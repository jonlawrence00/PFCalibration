from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'PGun_step3_RECO_1200_2021_2_200_Sep17_tmp1'
config.General.workArea = 'crab_projects'

#optional
#config.General.transferOutputs
#config.General.transferLogs
#config.General.failureLimit = 

#Expert use
#config.General.instance
#config.General.activity

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'step3_RAW2DIGI_L1Reco_RECO_RECOSIM_EI_PAT_VALIDATION_DQM.py'
#config.JobType.psetName = 'myEDAna.py'
config.JobType.outputFiles = ['step3.root']
#config.JobType.eventsPerLumi = 2000

config.section_("Data")
#config.Data.inputDataset = '/Single_Pion_gun_E_200to500_13TeV_pythia8/RunIIWinter19PFCalibDR-2018ConditionsNoPU_105X_upgrade2018_realistic_v4-v1/GEN-SIM-RECO'
#config.Data.primaryDataset = ''
#config.Data.splitting = 'EventBased'
config.Data.userInputFiles = open('/afs/cern.ch/user/b/bkansal/work/PFcalibration_run3/CMSSW_12_0_0/src/PFCalibration/PFChargedHadronAnalyzer/test/step2_2_200_10million_1200.txt').readlines()
# config.Data.ignoreLocality = True
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
NJOBS = 5000
# config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.publication = False
#config.Data.publishDBS = '' default for the moment
#config.Data.outLFN = '/home/spandey/t3store/PF_PGun'
config.Data.outLFNDirBase = '/store/user/bkansal/step3/PGun_step3_RECO_1200_2021_2_200_Sep17_tmp1/'

config.section_("Site")
config.Site.storageSite = 'T2_IN_TIFR'
# config.Site.blacklist = ['T2_CH_*', 'T2_US_*']
config.Site.whitelist = ['T2_US_Wisconsin','T2_US_MIT']
#config.Site.whitelist = ['T2_US_*','T2_UK_*']                                            
# config.Site.whitelist = ['T2_UA_KIPT','T2_IN_TIFR']
                                                                                                  
# config.Site.whitelist = ['T2_AT_Vienna','T2_UK_London_IC']


#config.section_("User")
#config.section_("Debug")
