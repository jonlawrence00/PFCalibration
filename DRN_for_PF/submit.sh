#!/bin/bash -l
#SBATCH --time=48:00:00
#SBATCH --ntasks=4
#SBATCH --mem=100g
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mail-user=chirayu.gupta@gmail.com
#SBATCH -p gpu
#SBATCH --gres=gpu:1

#v2->trueE, no semi
#v1->semi_dscb, loss fun= dscb

#export PYTHONUNBUFFERED=1
#module load ohpc
#module load iiser/apps/cuda/11.4
#module load cmake/3.14.3
#module swap gnu8 cdac/compiler/gcc/10.2.0
#module load python/3.9.8

./train trained pickles_test --idx_name all --nosemi --target trueE --valid_batch_size 400 --train_batch_size 400 --acc_rate 2 --max_lr 0.0001 --n_epochs 100 --lr_sched Const  --in_layers 3 --mp_layers 2 --out_layers 2 --agg_layers 2 &>> trained/training.log
