
import uproot
from numba import jit
import numpy as np
import awkward as ak
from time import time
import pickle
import tqdm
from torch_geometric.data import Data
import torch
objects=[]
with open('/home/rusack/shared/pickles/HGCAL_TestBeam/Test_alps/Hit_X.pickle', 'rb') as f:
    while True:
        try:
            objects.append(pickle.load(f))
            #            print(objects)
        except EOFError:
            break
            

print(len(objects))
