#!/usr/bin/env python3

import os
import subprocess
from datetime import date
import numpy as np

#  --in_layers IN_LAYERS
#                        Depth of input network
#  --agg_layers AGG_LAYERS
#                        Number of aggregation layers
#  --mp_layers MP_LAYERS
#                        Depth of message-passing networks
#  --out_layers OUT_LAYERS
#                        Depth of output network
#  --train_batches TRAIN_BATCHES
#                        Number of batches per training epoch. Incompatible with --train_batch_size
#  --n_epochs N_EPOCHS   
#                         Number of training epochs
#
#
#
#
#  --max_lr MAX_LR       Max learning rate
#  --min_lr MIN_LR       Min learning rate

import argparse
############################################################################
#ARGS
parser = argparse.ArgumentParser(description='So, you\'re staring at your hot new GNN thinking... I need to train this sucker. But like, there\'s like, a billion params to search over. I got you fam. This program has a few main modes for that: \t\n1. train, \t\n2. oops, just inference plz, \t\n3. make me plots... and a sandwich... \t\n4. No no. No mode 4... you need to do work too sometimes... go drink some water, and maybe call your parents every once and a while.')

parser.add_argument('--mode', type=str, choices=['train','inference','plots'], help='pick one: train,inference,plots')
parser.add_argument('--mlr_sweep', type=str, help='A list of max learning rates like: 1e-4,1e-5,1e-6')
parser.add_argument('--tb_sweep', type=str, help='A list of max of the number of training batches like: 500,1000,5000')
parser.add_argument('--ne_sweep', type=str, help='A list of the number of epochs like: 35,50,100')
parser.add_argument('--base_slurm', type=str, help='The relative path to whatever base slurm command you want to use i.e. submit_train_base.slurm')
parser.add_argument('--no_submit', action='store_true', help='When true, it does not call sbatch, it just makes the .slurm files')

args = parser.parse_args()
###########################################################################################################################
#INPUTS
mode = args.mode
max_lrs = args.mlr_sweep.split(',')
train_batches = args.tb_sweep.split(',')
n_epochs = args.ne_sweep.split(',')
no_submit = args.no_submit 


slurm_base = args.base_slurm

folder_prefix="/home/rusack/rothmans/Pho_Refined/grid_30M_EE/training_folder/"
submit_folder = "slurms_{}".format(date.today().strftime("%d_%m_%Y"))

#make permutations
sweeps = [max_lrs, train_batches, n_epochs]
#black magic
paramsets = np.array(np.meshgrid(*sweeps)).T.reshape(-1,len(sweeps))
print(paramsets)
print(no_submit)

if(mode == 'train'):
    if not os.path.exists(submit_folder):
        os.mkdir(submit_folder)
    
    base_comm = [line for line in open(slurm_base, "r")]
    os.chdir(submit_folder)    
    for max_lr, train_batch, n_epoch in paramsets:
        folder = folder_prefix+"maxlr_{}_trainbatches_{}_nepoch_{}".format(max_lr, train_batch, n_epoch)
        if not os.path.exists(folder):
            os.mkdir(folder)
        print(folder)
    
        train_comm = "train {} 2018_UL_Photon_Refined --idx_name EE --train_batches {} --max_lr {} --min_lr 1e-7 --n_epochs {} --lr_sched Cyclic --ES no --semi_dscb_sigmoid --thresh 1.5 --best_arch &>> {}/training.log".format(folder, train_batch, max_lr, n_epoch, folder)
        slurm_name = "submit_train_maxlr_{}_trainbatches_{}_nepoch_{}.slurm".format(max_lr, train_batch, n_epoch)
        submit_file = open(slurm_name, "w")
    
        submit_file.writelines(base_comm)
        submit_file.write("\n\n")
        submit_file.write(train_comm)
        submit_file.close()
        #quit before you submit
        if(no_submit):
           continue
        #or, to the grid!
        subprocess.run(["sbatch", slurm_name])
if(mode == 'inference'):
    print("Starting Inference")
    if not os.path.exists(submit_folder):
        os.mkdir(submit_folder)
    
    base_comm = [line for line in open(slurm_base, "r")]
    os.chdir(submit_folder)    
    for max_lr, train_batch, n_epoch in paramsets:
        folder = folder_prefix+"maxlr_{}_trainbatches_{}_nepoch_{}".format(max_lr, train_batch, n_epoch)
    
        train_comm = "train {} 2018_UL_Photon_Refined --idx_name EE --train_batches {} --max_lr {} --min_lr 1e-7 --n_epochs {} --lr_sched Cyclic --ES no --semi_dscb_sigmoid --thresh 1.5 --best_arch --predict_only &>> {}/predict.log".format(folder, train_batch, max_lr, n_epoch, folder)
        slurm_name = "submit_predict_maxlr_{}_trainbatches_{}_nepoch_{}.slurm".format(max_lr, train_batch, n_epoch)
        submit_file = open(slurm_name, "w")
    
        submit_file.writelines(base_comm)
        submit_file.write("\n\n")
        submit_file.write(train_comm)
        submit_file.close()
        #quit before you submit
        if(no_submit):
           print("NO SUBMIT")
           continue
        #or, to the grid!
        print("Submitting to the grid...")
        subprocess.run(["sbatch", slurm_name])


if(mode == 'plots'):
    folders = ""
    labels = ""
    for max_lr, train_batch, n_epoch in paramsets:
        folder = folder_prefix+"maxlr_{}_trainbatches_{}_nepoch_{}".format(max_lr, train_batch, n_epoch)
   #./validate convergence --folders ../training_folder ../training_folder --labels first second
        #check for a summaries file in the folder
        if not os.path.isfile(folder+"/summaries.npz"):
            print("WARNING: NO SUMMARY FOUND IN---\n\t"+folder)
            continue
        folders += (folder+" ")
        labels += "MLR_{}_TB_{}_NE_{} ".format(max_lr, train_batch, n_epoch)

    plot_comm = "validate convergence --folders {} --labels {} --semiparam".format(folders, labels)
    print(plot_comm)
    subprocess.run((plot_comm.split()))
