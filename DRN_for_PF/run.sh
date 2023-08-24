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
./prepareHGCAL
#./prepareHGCAL_extract &>> ./pickle_ec_in_pos.log

<<EOF
module load python/3.7
module load iiser/apps/ROOT/6.24.02
#python plotting_script/combined_pickle.py /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/pickle_eta1pt55to2pt5_pos/trueE_target.pickle /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/pickle_eta1pt55to2pt5_neg/trueE_target.pickle
#python plotting_script/full_binned.py /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/Trained tmp.root /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/pickle/
#python plotting_script/full_binned.py /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/Trained_eta1 Output_eta1.root /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/pickle_eta1/
#python plotting_script/full_binned.py /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/Trained_eta1pt55_trueE Output_eta1pt55_trueE.root /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/pickle_eta1pt55/
#python plotting_script/full_binned.py /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/Trained_eta2pt5 Output_eta2pt5.root /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/pickle_eta2pt5/
#python plotting_script/2d_binned.py /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/Trained_eta2pt75_acc5_bs1000 Output_eta2pt75_eta.root /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/pickle_test/
python plotting_script/2d_binned.py /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/Trained_eta1pt55 Output_eta1pt55_eta.root /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/pickle_eta1pt55/
#python plotting_script/full_binned.py /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/Trained_eta1pt55to2pt5_pos Output_eta1pt55to2pt5_pos.root /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/pickle_eta1pt55to2pt5_pos/
#python plotting_script/full_binned_v2.py /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/Trained_eta1pt55to2pt5_pos /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/Trained_eta1pt55to2pt5_neg Output_eta1pt55to2pt5.root /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/pickle_eta1pt55to2pt5_pos/ /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/pickle_eta1pt55to2pt5_neg/
#python plotting_script/eta_binned.py /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/Trained_eta1pt55_trueE Output_eta1pt55_trueE_eta.root /home/bkansal/work/Bhumika/The_DRN_for_HGCAL/pickle_eta1pt55/
