from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")

config.General.requestName = 'PGun_step3_RECO_1264_200_500_rechit'
config.General.workArea = 'crab_projects'

#optional
#config.General.transferOutputs
#config.General.transferLogs
#config.General.failureLimit = 

#Expert use
#config.General.instance
#config.General.activity

config.section_("JobType")
#config.JobType.psetName = 'myEDAna.py'
#config.JobType.psetName = 'step3_RAW2DIGI_L1Reco_RECO_EI_PAT_VALIDATION_DQM.py'
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'step3_RAW2DIGI_L1Reco_RECO_RECOSIM_EI_PAT_VALIDATION_DQM.py'
config.JobType.outputFiles = ['step3.root']
#config.JobType.eventsPerLumi = 2000
config.JobType.maxMemoryMB=3000
config.section_("Data")

#config.Data.inputDataset = '/Single_Pion_gun_E_2to200_14TeV_pythia8/Run3Summer21DR-NoPURAWRECO_120X_mcRun3_2021_realistic_v6-v2/GEN-SIM-RAW'
#config.Data.inputDataset = '/SinglePionGun_E0p2to200/Run3Winter23Digi-NoPUGTv4_126X_mcRun3_2023_forPU65_v4-v2/GEN-SIM-RAW'
config.Data.inputDataset = '/SinglePionGun_E200to500/Run3Winter23Digi-NoPUGTv4_126X_mcRun3_2023_forPU65_v4-v2/GEN-SIM-RAW'
#config.Data.inputDataset = '/Single_Pion_gun_E_2to200_14TeV_pythia8/Run3Summer21DR-NoPURAWRECO_120X_mcRun3_2021_realistic_v6-v2/GEN-SIM-RAW'
#config.Data.inputDataset = '/Single_Pion_gun_E_2to200_14TeV_pythia8/Run3Summer21DR-NoPURAWRECO_120X_mcRun3_2021_realistic_v6-v2/GEN-SIM-RECO'
#config.Data.inputDataset = '/Single_Pion_gun_E_200to500_14TeV_pythia8/Run3Summer21DR-NoPURAWRECO_120X_mcRun3_2021_realistic_v6-v2/GEN-SIM-RECO'
#config.Data.inputDataset = '/Single_Pion_gun_E_2to200_13TeV_pythia8/RunIIWinter19PFCalibDR-2016ConditionsNoPU_105X_mcRun2_asymptotic_v2-v1/GEN-SIM-RECO'
#config.Data.primaryDataset = ''
#config.Data.splitting = 'EventBased'
# config.Data.userInputFiles = open('/afs/cern.ch/user/b/bkansal/work/PFcalibration_2019/CMSSW_10_6_0_pre2/src/PFCalibration/PFChargedHadronAnalyzer/test/step2_2_200_10million_1060_2018.txt').readlines()
# config.Data.ignoreLocality = True
#config.Data.partialDataset = False

config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#NJOBS = 1681
#config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.publication = False
#config.Data.publishDBS = '' default for the moment
#config.Data.outLFN = '/home/spandey/t3store/PF_PGun'
config.Data.outLFNDirBase = '/store/user/bkansal/step3_rechit/PGun_step3_RECO_1264_200_500_rechit/'


config.section_("Site")
config.Site.storageSite = 'T2_IN_TIFR'
# config.Site.blacklist = ['T2_CH_*', 'T2_US_*']
# # config.Site.whitelist = ['T2_UK_London_IC','T2_US_*']                                            
# config.Site.whitelist = ['T2_UA_KIPT','T2_IN_TIFR']
                                                                                                  
# config.Site.whitelist = ['T2_AT_Vienna','T2_UK_London_IC']


#config.section_("User")
#config.section_("Debug")
