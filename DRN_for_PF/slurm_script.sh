#!/bin/bash -l                                                                                                                                

#SBATCH --time=48:00:00     ## time limit                                                                                                            
#SBATCH --ntasks-per-node=10               ## nodes - can change it                                                                                             
#SBATCH --mail-type=ALL                                                                                                                         
#SBATCH --mail-user=bhumika.kansal@students.iiserpune.ac.in ## change it your mail ID
#SBATCH --job-name="ratio_full"                                                                                                               
#SBATCH --partition=gpu                                                                                                                      
#SBATCH --gres=gpu:4                                                                                                                         
#SBATCH --mem=100g                                                                                                                           
#SBATCH --output=pytorch_gpu_%j.out                                                                                                          
#SBATCH --error=pytorch_gpu_%j.err                                                                                                           

module load cdac/spack/0.17
source /home/apps/spack/share/spack/setup-env.sh
spack load python@3.8.2
source /home/apps/iiser/pytorch-venv/bin/activate
export CUDA_VISIBLE_DEVICES=0,1,2,3  ## this line is to use 4 GPU nodes

echo '-========================================================='
cd /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/

#./prepareHGCAL &>> /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/pickle_full.log

#echo '${PWD}'
#cd ${PWD}

#ec_out
#./train /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/Trained_eta2pt5to2pt75_pos /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/pickle_eta2pt5to2pt75_pos --idx_name all --nosemi --target trueE  --valid_batch_size 400 --train_batch_size 400 --acc_rate 2 --max_lr 0.0001 --n_epochs 20 --lr_sched Const  --in_layers 3 --mp_layers 2 --out_layers 2 --agg_layers 2 &>> /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/Trained_eta2pt5to2pt75_pos/training.log

#ec_in
#./train /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/Trained_eta1pt55to2pt5_pos /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/pickle_eta1pt55to2pt5_pos --idx_name all --nosemi --target trueE  --valid_batch_size 400 --train_batch_size 400 --acc_rate 2 --max_lr 0.0001 --n_epochs 20 --lr_sched Const  --in_layers 3 --mp_layers 2 --out_layers 2 --agg_layers 2 &>> /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/Trained_eta1pt55to2pt5_pos/training.log


#./train /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/Trained_eta2pt75_acc5_bs1000 /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/pickle_test --idx_name all --nosemi --target trueE  --valid_batch_size 1000 --train_batch_size 1000 --acc_rate 5 --n_epochs 10 --max_lr 0.0001 --lr_sched Const  --in_layers 3 --mp_layers 2 --out_layers 2 --agg_layers 2 &>> /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/Trained_eta2pt75_acc5_bs1000/training.log    

#./train /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/Trained_eta2pt75_epoch100 /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/pickle_eta2pt75 --idx_name all --nosemi --target trueE  --valid_batch_size 400 --train_batch_size 400 --acc_rate 2 --max_lr 0.0001 --n_epochs 100 --lr_sched Const  --in_layers 3 --mp_layers 2 --out_layers 2 --agg_layers 2 &>> /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/Trained_eta2pt75_epoch100/training.log

./train /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/Trained_eta2pt75_target_ratioflip_epoch100 /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/pickle_eta2pt75 --idx_name all --nosemi --target ratioflip  --valid_batch_size 400 --train_batch_size 400 --acc_rate 2 --max_lr 0.0001 --n_epochs 100 --lr_sched Const  --in_layers 3 --mp_layers 2 --out_layers 2 --agg_layers 2 &>> /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/Trained_eta2pt75_target_ratioflip_epoch100/training.log

#./train /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/Trained /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/pickle --idx_name all --nosemi --target ratioflip  --valid_batch_size 400 --train_batch_size 400 --acc_rate 2 --max_lr 0.0001 --lr_sched Const  --in_layers 3 --mp_layers 2 --out_layers 2 --agg_layers 2 --predict_only &>> /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/Trained/training.log





