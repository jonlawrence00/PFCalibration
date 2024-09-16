import os
import os.path as osp
import math

import numpy as np
import torch
import gc
import torch.nn as nn
import torch.nn.init as init
from torch.nn.functional import softplus
import torch_geometric.transforms as T

from torch.utils.checkpoint import checkpoint
from torch_cluster import knn_graph, graclus_cluster
from torch_scatter import scatter
from torch_sparse.storage import SparseStorage

from torch import Tensor
from torch_geometric.typing import OptTensor, Optional, Tuple


from torch_geometric.nn import EdgeConv, NNConv
from torch_geometric.nn.pool.pool import pool_batch
from torch_geometric.nn.pool.consecutive import consecutive_cluster
from torch_geometric.utils.num_nodes import maybe_num_nodes
from torch_geometric.utils import normalized_cut
from torch_geometric.utils import remove_self_loops
from torch_geometric.nn import (max_pool, max_pool_x, global_max_pool,
                                avg_pool, avg_pool_x, global_mean_pool, 
                                global_add_pool)

transform = T.Cartesian(cat=False)

def normalized_cut_2d(edge_index, pos):
    row, col = edge_index[0], edge_index[1]
    edge_attr = torch.norm(pos[row] - pos[col], p=2, dim=1)
    return normalized_cut(edge_index, edge_attr, num_nodes=pos.size(0))

# jit compatible version of coalesce
def coalesce(index, value: OptTensor, m: int, n: int, op: str = "add"):
    storage = SparseStorage(row=index[0], col=index[1], value=value,
                            sparse_sizes=(m, n), is_sorted=False)
    storage = storage.coalesce(reduce=op)
    return torch.stack([storage.row(), storage.col()], dim=0), storage.value()

# jit compatible version of to_undirected
def to_undirected(edge_index, num_nodes: Optional[int] = None) -> Tensor:
    num_nodes = maybe_num_nodes(edge_index, num_nodes)

    row, col = edge_index[0], edge_index[1]
    temp = torch.cat([row, col], dim=0), torch.cat([col, row], dim=0)
    row, col = temp[0], temp[1]
    edge_index = torch.stack([row, col], dim=0)
    edge_index, _ = coalesce(edge_index, None, num_nodes, num_nodes)
    return edge_index

# jit compatible version of pool_edge, depends on coalesce
def pool_edge(cluster, edge_index, edge_attr: Optional[torch.Tensor] = None):
    num_nodes = cluster.size(0)
    edge_index = cluster[edge_index.view(-1)].view(2, -1)
    edge_index, edge_attr = remove_self_loops(edge_index, edge_attr)
    if edge_index.numel() > 0:
        edge_index, edge_attr = coalesce(edge_index, edge_attr, num_nodes,
                                         num_nodes)
    return edge_index, edge_attr

def _aggr_pool_x(cluster, x, aggr: str, size: Optional[int] = None):
    """Call into scatter with configurable reduction op"""
    #print("aggr: ", aggr)
    return scatter(x, cluster, dim=0, dim_size=size, reduce=aggr)

def global_pool_aggr(x, batch: OptTensor, aggr: str, size: Optional[int] = None):
    """Global pool via passed aggregator: 'mean', 'add', 'max'"""
    if batch is None and size is None:
        raise Exception('Must provide at least one of "batch" or "size"')
    if batch is not None:
        size = int(batch.max().item() + 1)
    assert batch is not None
    return scatter(x, batch, dim=0, dim_size=size, reduce=aggr)

# this function is specialized compared to the more general non-jittable version
# in particular edge_attr can be removed since it is always None
def aggr_pool(cluster, x, batch: OptTensor, aggr: str) -> Tuple[Tensor, OptTensor]:
    """jit-friendly version of max/mean/add pool"""
    #print("cluster before; ", len(cluster))
    cluster, perm = consecutive_cluster(cluster)
    #print("cluster after; ", len(cluster))
    #print("perm: ,", len(perm))
    x = _aggr_pool_x(cluster, x, aggr)
    if batch is not None:
        batch = pool_batch(perm, batch)
    return x, batch

def aggr_pool_x(cluster, x, batch: OptTensor, aggr: str, size: Optional[int] = None) -> Tuple[Tensor, OptTensor]:
    """*_pool_x with configurable aggr method"""
    if batch is None and size is None:
        raise Exception('Must provide at least one of "batch" or "size"')
    if size is not None and batch is not None:
        batch_size = int(batch.max().item()) + 1
        return _aggr_pool_x(cluster, x, aggr, batch_size * size), None

    cluster, perm = consecutive_cluster(cluster)
    x = _aggr_pool_x(cluster, x, aggr)
    if batch is not None:
        batch = pool_batch(perm, batch)

    return x, batch
    
class DynamicReductionNetworkJit(nn.Module):
    '''
    This model iteratively contracts nearest neighbour graphs 
    until there is one output node.
    The latent space trained to group useful features at each level
    of aggregration.
    This allows single quantities to be regressed from complex point counts
    in a location and orientation invariant way.
    One encoding layer is used to abstract away the input features.

    @param input_dim: dimension of input features
    @param hidden_dim: dimension of hidden layers
    @param output_dim: dimensio of output
    
    @param k: size of k-nearest neighbor graphs
    @param aggr: message passing aggregation scheme. 
    @param norm: feature normaliztion. None is equivalent to all 1s (ie no scaling)
    @param loop: boolean for presence/absence of self loops in k-nearest neighbor graphs
    @param pool: type of pooling in aggregation layers. Choices are 'add', 'max', 'mean'
    
    @param agg_layers: number of aggregation layers. Must be >=0
    @param mp_layers: number of layers in message passing networks. Must be >=1
    @param in_layers: number of layers in inputnet. Must be >=1
    @param out_layers: number of layers in outputnet. Must be >=1
    '''
    latent_probe: Optional[int]
    def __init__(self, input_dim=4, hidden_dim=64, output_dim=1, k=16, aggr='add', norm=None, 
            loop=True, pool='max',
            agg_layers=2, mp_layers=2, in_layers=1, out_layers=3,
            graph_features = 0,
            latent_probe=None):
        super(DynamicReductionNetworkJit, self).__init__()

        self.graph_features = graph_features

        if latent_probe is not None and (latent_probe>agg_layers+1 or latent_probe<-1*agg_layers-1):
            print("Error: asked for invalid latent_probe layer")
            return
        
        if latent_probe is not None and latent_probe < 0:
            latent_probe = agg_layers+1 - latent_probe

        if latent_probe is not None:
            print("Probing latent features after %dth layer"%latent_probe)

        self.latent_probe = latent_probe

        self.loop = loop

        print("Pooling with",pool)
        print("Using self-loops" if self.loop else "Not using self-loops")
        print("There are",agg_layers,'aggregation layers')

        if norm is None:
            norm = torch.ones(input_dim)

        #normalization vector
        self.datanorm = nn.Parameter(norm)
        
        self.k = k

        #construct inputnet
        in_layers_l = []
        temp_layer = nn.Linear(input_dim, hidden_dim)
        init.kaiming_uniform_(temp_layer.weight)
        in_layers_l += [temp_layer,nn.ReLU()]

        for i in range(in_layers-1):
            temp_layer = nn.Linear(hidden_dim, hidden_dim)
            init.kaiming_uniform_(temp_layer.weight)
            in_layers_l += [temp_layer, 
                    nn.ReLU()]

        self.inputnet = nn.Sequential(*in_layers_l)
        #print("hi there big dog")
        #print("input net: ", self.inputnet)
        
        #construct aggregation layers
        self.agg_layers = nn.ModuleList()
        print("agg layers: ", self.agg_layers)
        for i in range(agg_layers):
            e = 2**(i+1)
            #print("hello there mate 1")
            #construct message passing network
            mp_layers_l = []

            for j in range(mp_layers-1):
                #print("hello there mate 2, He init is done")
                temp_layer = nn.Linear(2*hidden_dim, 2*hidden_dim)
                init.kaiming_uniform_(temp_layer.weight)
                mp_layers_l += [temp_layer, nn.BatchNorm1d(2*hidden_dim), nn.ReLU()]

            temp_layer = nn.Linear(2*hidden_dim, hidden_dim)
            init.kaiming_uniform_(temp_layer.weight)
            mp_layers_l += [temp_layer, nn.BatchNorm1d(hidden_dim), nn.ReLU()]
           
            convnn = nn.Sequential(*mp_layers_l)
            
            self.agg_layers.append(EdgeConv(nn=convnn, aggr=aggr).jittable())
            #print("aggr: ", aggr)
            #print("convnn: " , convnn)

        print("agg layers: ", self.agg_layers)
        #print("mp layers: ",mp_layers)
        
        #construct outputnet
        out_layers_l = []

        #temp_layer = nn.Linear(4*hidden_dim+self.graph_features, 2*hidden_dim+self.graph_features)
        #init.kaiming_uniform_(temp_layer.weight)
        #out_layers_l += [temp_layer, nn.BatchNorm1d(4*hidden_dim), nn.ReLU(), nn.Dropout(.1)]

        #temp_layer = nn.Linear(2*hidden_dim+self.graph_features, hidden_dim+self.graph_features)
        #init.kaiming_uniform_(temp_layer.weight)
        #out_layers_l += [temp_layer, nn.BatchNorm1d(1), nn.ReLU()]

        for i in range(out_layers-1):
            temp_layer = nn.Linear(hidden_dim+self.graph_features, hidden_dim+self.graph_features)
            init.kaiming_uniform_(temp_layer.weight)
            out_layers_l += [temp_layer, nn.ReLU(), nn.Dropout(.1)]


        temp_layer = nn.Linear(hidden_dim+self.graph_features, output_dim)
        init.kaiming_uniform_(temp_layer.weight)
        out_layers_l += [temp_layer]

        self.output = nn.Sequential(*out_layers_l)

        if pool not in {'max', 'mean', 'add'}:
            raise Exception("ERROR: INVALID POOLING")
        
        self.aggr_type = pool

        print("output: ", self.output)

        '''#construct inputnet
        in_layers_l = []
        in_layers_l += [nn.Linear(input_dim, hidden_dim),
                nn.ELU()]

        for i in range(in_layers-1):
            in_layers_l += [nn.Linear(hidden_dim, hidden_dim), 
                    nn.ELU()]

        self.inputnet = nn.Sequential(*in_layers_l)

        print("input net: ", self.inputnet)
        
        #construct aggregation layers
        self.agg_layers = nn.ModuleList()
        print("agg layers: ", self.agg_layers)
        for i in range(agg_layers):
            print("hello there mate 1")
            #construct message passing network
            mp_layers_l = []

            for j in range(mp_layers-1):
                print("hello there mate 2")
                mp_layers_l += [nn.Linear(2*hidden_dim, 2*hidden_dim), nn.ELU()]

            mp_layers_l += [nn.Linear(2*hidden_dim, hidden_dim),
                    nn.ELU()]
           
            convnn = nn.Sequential(*mp_layers_l)
            
            self.agg_layers.append(EdgeConv(nn=convnn, aggr=aggr).jittable())
            print("aggr: ", aggr)
            print("convnn: " , convnn)

        print("agg layers: ", self.agg_layers)
        #print("mp layers: ",mp_layers)
        
        #construct outputnet
        out_layers_l = []

        for i in range(out_layers-1):
            out_layers_l += [nn.Linear(hidden_dim+self.graph_features, hidden_dim+self.graph_features), nn.ELU()]
            #out_layers_l += [nn.Linear(hidden_dim+self.graph_features, hidden_dim+self.graph_features), nn.ELU(), nn.Dropout(.1)]

        out_layers_l += [nn.Linear(hidden_dim+self.graph_features, output_dim)]

        self.output = nn.Sequential(*out_layers_l)

        if pool not in {'max', 'mean', 'add'}:
            raise Exception("ERROR: INVALID POOLING")
        
        self.aggr_type = pool

        print("output: ", self.output)'''

        '''#construct inputnet
        in_layers_l = []
        in_layers_l += [nn.Linear(input_dim, hidden_dim),
                nn.ELU()]

        for i in range(in_layers-1):
            in_layers_l += [nn.Linear(hidden_dim, hidden_dim), 
                    nn.ELU()]

        self.inputnet = nn.Sequential(*in_layers_l)

        print("input net: ", self.inputnet)
        
        #construct aggregation layers
        self.agg_layers = nn.ModuleList()
        print("agg layers: ", self.agg_layers)
        for i in range(agg_layers):
            e = 2**(i+1)
            print("hello there mate 1")
            #construct message passing network
            mp_layers_l = []

            for j in range(mp_layers-1):
                print("hello there mate 2")
                mp_layers_l += [nn.Linear(e*hidden_dim, e*hidden_dim), nn.ELU()]

            #mp_layers_l += [nn.Linear(2*hidden_dim, hidden_dim),nn.ELU()]
           
            convnn = nn.Sequential(*mp_layers_l)
            
            self.agg_layers.append(EdgeConv(nn=convnn, aggr=aggr).jittable())
            print("aggr: ", aggr)
            print("convnn: " , convnn)

        print("agg layers: ", self.agg_layers)
        #print("mp layers: ",mp_layers)
        
        #construct outputnet
        out_layers_l = []

        out_layers_l += [nn.Linear(4*hidden_dim+self.graph_features, 2*hidden_dim+self.graph_features), nn.ELU()]

        out_layers_l += [nn.Linear(2*hidden_dim+self.graph_features, hidden_dim+self.graph_features), nn.ELU()]
        
        for i in range(out_layers-1):
            out_layers_l += [nn.Linear(hidden_dim+self.graph_features, hidden_dim+self.graph_features), nn.ELU()]
            #out_layers_l += [nn.Linear(hidden_dim+self.graph_features, hidden_dim+self.graph_features), nn.ELU(), nn.Dropout(.1)]

        out_layers_l += [nn.Linear(hidden_dim+self.graph_features, output_dim)]

        self.output = nn.Sequential(*out_layers_l)

        if pool not in {'max', 'mean', 'add'}:
            raise Exception("ERROR: INVALID POOLING")
        
        self.aggr_type = pool

        print("output: ", self.output)'''

    def forward(self, x: Tensor, batch: OptTensor, graph_x: OptTensor) -> Tensor:
        '''
        Push the batch 'data' through the network
        '''
        #print("################################")
        #print("initial size of x: ", x.size())
        #print("initial x: ",len(x))
        x = self.datanorm * x
        #print("self.datanorm: ", self.datanorm)
        x = self.inputnet(x)
        #print("input x", len(x))
        #print("initial size of x: ", x.size())
        latent_probe = self.latent_probe
        
        if graph_x is not None:
            #print("graph is none")
            graph_x = graph_x.view((-1, self.graph_features))

        # if there are no aggregation layers just leave x, batch alone
        nAgg = len(self.agg_layers)
        #counter = 0
        for i, edgeconv in enumerate(self.agg_layers):
            initial = x
            #counter += 1
            if latent_probe is not None and i == latent_probe:
                return x
            #print("batch len", len(batch))
            #print("size of batch: ", batch.size())
            #print(batch)
            knn = knn_graph(x, self.k, batch, loop=self.loop, flow=edgeconv.flow)
            #print("knn : " ,knn)
            #print("self.kL " , len(self.k))
            #print(self.k)
            #print("len of x:", len(x))
            #print("knn: ", knn[0].shape())
            edge_index = to_undirected(knn)
            #print("#i: ", i, " #edgeconv: ", edgeconv)
            #print("x[0] before edgeconv: ", x[0])
            #print("x.size before edgeconv: ", x.size())
            x = edgeconv(x, edge_index)
            #print("x[0] after edgeconv: ", x[0])
            x+=initial
            #print("x.size after edgeconv: ", x.size())
            #print("edge index: ", edge_index)


            weight = normalized_cut_2d(edge_index, x)
            #print("weight: ", weight.size())
            cluster = graclus_cluster(edge_index[0], edge_index[1], weight, x.size(0))
            #print("cluster: ", cluster.size())
            #print(cluster[0:100])

            #print("batch: ", batch)
            if i == nAgg - 1:
                x, batch = aggr_pool_x(cluster, x, batch, self.aggr_type)
            else:
                #print("hello?")
                x, batch = aggr_pool(cluster, x, batch, self.aggr_type)

            #print("batch: ", batch[0:10])
            #print("x.size after aggr_pool: ", x.size())
            #print("x[0] after aggregate: ", x[0])
            #print("batch.size after aggr_pool: ", batch)

        if latent_probe is not None and latent_probe == nAgg:
            #print("latent is none")
            return x

        # this xforms to batch-per-row so no need to return batch
        #print("after agg x.size", x.size())
        x = global_pool_aggr(x, batch, self.aggr_type)
        #print("after global x.size", x.size())

        if latent_probe is not None and latent_probe == nAgg + 1:
            #print("latent is none")
            return x

        #print("size of graph_x; ", graph_x)
        if graph_x is not None:
            #print("graph is not none 2")
            x = torch.cat((x, graph_x), 1)

        x = self.output(x).squeeze(-1)
        #print("size of x at the end: ", x.size())
        #print("####################################################")
        return x
