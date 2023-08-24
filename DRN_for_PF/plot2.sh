#!/bin/bash                                                                                                                                    

#module load python/3.7      
module load iiser/apps/cuda/11.4
set umask 0002 
module load cdac/spack/0.17
source /home/apps/spack/share/spack/setup-env.sh
spack load python@3.8.2
source /home/apps/iiser/pytorch-venv/bin/activate
echo '-========================================================='
cd /home/bkansal/work/Bhumika/The_DRN_for_HGCAL
#cd <path to the directory where preparehGCAL & other codes are available>
#python tmp.py
#./prepareHGCAL &>> ./pickle_barrel.log
#./prepareHGCAL_extract
#./prepareHGCAL_extract &>> ./pickle_ec_in_pos.log

#<<EOF
module load python/3.7
module load iiser/apps/ROOT/6.24.02


#Trained_eta2pt75_target_ratio_epoch100
#python plotting_script/2d_binned_EH.py /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/Trained_eta2pt75_target_ratio_epoch100/ Output_eta2pt75_EH.root /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/pickle_eta2pt75/
#python plotting_script/2d_binned_H.py /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/Trained_eta2pt75_target_ratio_epoch100/ Output_eta2pt75_H.root /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/pickle_eta2pt75/

#Trained_2_350_eta2pt75_epoch100

#python plotting_script/2d_binned_EH.py /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/Trained_2_350_eta3_epoch100/  Output_eta2pt75_EH.root /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/pickle_2_350_eta3/
#python plotting_script/2d_binned_H.py /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/Trained_2_350_eta3_epoch100/ Output_eta2pt75_H.root /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/pickle_2_350_eta3/

#python plotting_script/2d_binned_EH.py /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/Trained_2_350_eta3_epoch100/  Output_eta2pt75_EH_trans.root /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/pickle_2_350_eta3/
#python plotting_script/2d_binned_H.py /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/Trained_2_350_eta3_epoch100/  Output_eta2pt75_H_trans.root /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/pickle_2_350_eta3/
python plotting_script/eta_binned_EH.py /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/Trained_2_350_eta3_epoch100/  Output_eta2pt75_EH_eta.root /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/pickle_2_350_eta3/
python plotting_script/eta_binned_H.py /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/Trained_2_350_eta3_epoch100/  Output_eta2pt75_H_eta.root /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/pickle_2_350_eta3/
