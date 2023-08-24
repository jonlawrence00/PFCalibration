from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName ='PGun_step1_GEN_SIM_1200_2021_200_500_Sep17'
config.General.workArea = 'crab_projects'

#optional
#config.General.transferOutputs
#config.General.transferLogs
#config.General.failureLimit = 

#Expert use
#config.General.instance
#config.General.activity

config.section_("JobType")
config.JobType.pluginName = 'PrivateMC'
#config.JobType.psetName = 'SinglePiE50HCAL_pythia8_cfi_GEN_SIM.py'
config.JobType.psetName = 'SinglePi0E10_pythia8_cfi_GEN_SIM_200_500.py'
config.JobType.outputFiles = ['step1.root']
config.JobType.eventsPerLumi = 2000

config.section_("Data")
#config.Data.inputDataset = '/Single_Pion_gun_13TeV_pythia8/Fall14DR73-NoPU_MCRUN2_73_V9-v1/GEN-SIM-RAW-RECO'
#config.Data.primaryDataset = ''
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 2000
NJOBS = 5000
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.publication = False
#config.Data.publishDBS = '' default for the moment
#config.Data.outLFN = '/home/spandey/t3store/PF_PGun'
config.Data.outLFNDirBase = '/store/user/bkansal/step1/PGun_step1_GEN_SIM_1200_2021_200_500_Sep17/'

config.section_("Site")
config.Site.storageSite = 'T2_IN_TIFR'
#config.Site.blacklist = ['T3_US_UCR', 'T3_US_UMiss']
config.Site.whitelist = ['T2_US_Wisconsin','T2_US_MIT']                                              
#config.Site.whitelist = ['T2_BE_*','T2_HU_Budapest']
#config.Site.whitelist = ['T2_RU_ITEP','T2_RU_INR']

#config.section_("User")
#config.section_("Debug")
