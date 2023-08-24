# PFcalibration for charged hadrons using DRN
There are three folders in this repository. Here are some following details regarding those:
1. ```PFChargedHadronAnalyzer```: This will be used for preparing ntuples using EDAnalyzer code which is ```PFChargedHadronAnalyzer/plugins/PFChargedHadronAnalyzer.cc```. In the mentioned code, some selections are applied on the reco file (after step 3). The basic summary of those selections are mentioned in [https://gitlab.cern.ch/tdr/notes/AN-22-015/uploads/003facf903bb534b57fc33d8a5ea2810/AN-22-015_temp-18Feb.pdf](url) except the details of ECAL and HCAL PF rechits.
   Detailed instructions for making ntuples on **lxplus** cluster can be found at the following link: [https://docs.google.com/document/d/14EFU3Mi--4UHU0bUn07zHIF0wTyqvyzkPYU4KNpZlKY/edit?usp=sharing](url). Some precise instructions for creating ntuples after login to lxplus are as follows:
   ```
   mkdir PFcalibration
   cd PFcalibration
   export SCRAM_ARCH=slc7_amd64_gcc900
   cmsrel CMSSW_12_6_4
   cd CMSSW_12_6_4/src
   cmsenv
   git clone https://github.com/bkansal/PFcalibration_DRN.git
   scram b
   cd PFCalibration/PFChargedHadronAnalyzer/test
   cmsRun step3_RAW2DIGI_L1Reco_RECO_RECOSIM_EI_PAT_VALIDATION_DQM.py
   ```
   The output files from above commands can be found at the following path on **lxplus** : ```/eos/user/b/bkansal/step3_ana/```
   
3. ```DRN_for_PF```: This will be used for training and validating the DRN model on **NSM** cluster as described in ```train``` and other important files. For this,run following commands ("\\\" indicates comments):
    ```
    mkdir < pickle_folder > 
    ./run.sh \\ preparing pickle files as mentioned in Extract.py
    mkdir < Trained_folder >
    sbatch slurm_script.sh \\ to submit slurm job which contains all training parameter as utilized by train file
    ./plot.sh \\ To make all output root files that will be stored in < Trained_folder >
    ```
    The input skimmed ntuple file can be found at the following path on **nsm** : ```/home/bkansal/work/Bhumika/The_DRN_for_HGCAL/root_file/step3_2_350_eta3.root```.
   
5. ```Plotting```: It contains all the scripts to compare the response and resolution results using chi square methods and DRN methods. The commands to run those file are mentioned in ```Plotting/run.sh```.

