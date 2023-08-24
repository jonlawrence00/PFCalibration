# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 --datatier GEN-SIM-RECO,MINIAODSIM,DQMIO --runUnscheduled --conditions auto:run1_mc -s RAW2DIGI,L1Reco,RECO,RECOSIM,EI,PAT,VALIDATION:@standardValidationNoHLT+@miniAODValidation,DQM:@standardDQMFakeHLT+@miniAODDQM --eventcontent RECOSIM,MINIAODSIM,DQM -n 100 --filein file:step2.root --fileout file:step3.root
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('ana',eras.Run3)

# import of standard configurations
# process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
# process.load('Configuration.EventContent.EventContent_cff')
# process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
# process.load('Configuration.StandardSequences.MagneticField_cff')
# process.load('Configuration.StandardSequences.RawToDigi_cff')
# process.load('Configuration.StandardSequences.L1Reco_cff')
# process.load('Configuration.StandardSequences.Reconstruction_cff')
# process.load('Configuration.StandardSequences.RecoSim_cff')
# process.load('CommonTools.ParticleFlow.EITopPAG_cff')
# process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
# process.load('Configuration.StandardSequences.PATMC_cff')
# process.load('Configuration.StandardSequences.Validation_cff')
# process.load('DQMOffline.Configuration.DQMOfflineMC_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('root://se01.indiacms.res.in//store/user/spandey/step2/PGun_step2_DIGI_1002_2_200_Feb_12/CRAB_UserFiles/crab_PGun_step2_DIGI_1002_2_200_Feb_12/180212_110432/0000/step2_2.root'),
#    fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/RunIIWinter19PFCalibDR/Single_Pion_gun_E_200to500_13TeV_pythia8/GEN-SIM-RECO/2016ConditionsNoPU_105X_mcRun2_asymptotic_v2-v1/270000/FF853C26-CDC1-4D44-95D9-924C3C3A482F.root'),
#    fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/Run3Summer21DR/Single_Pion_gun_E_200to500_14TeV_pythia8/GEN-SIM-RECO/NoPURAWRECO_120X_mcRun3_2021_realistic_v6-v2/260000/0788ebab-26eb-410e-97b0-868877d34e33.root'),
#    fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/Run3Summer21DR/Single_Pion_gun_E_2to200_14TeV_pythia8/GEN-SIM-RECO/NoPURAWRECO_120X_mcRun3_2021_realistic_v6-v2/270000/00617b79-69fd-4239-a0a6-cde116b7abe1.root'),
#    fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/Run3Winter23Reco/SinglePionGun_E0p2to200/GEN-SIM-RECO/EpsilonPU_126X_mcRun3_2023_forPU65_v1-v2/2550000/007001ee-f0d5-4161-99ea-f286ad7136e6.root'),
    fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/Run3Winter23Reco/SinglePionGun_E0p2to200/GEN-SIM-RECO/NoPUGTv4_GTv4_126X_mcRun3_2023_forPU65_v4-v2/2540000/02bc0753-d6e6-4626-9ec3-33ccd7501292.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
# process.configurationMetadata = cms.untracked.PSet(
#     annotation = cms.untracked.string('step3 nevts:100'),
#     name = cms.untracked.string('Applications'),
#     version = cms.untracked.string('$Revision: 1.19 $')
# )

# Output definition

# process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
#     dataset = cms.untracked.PSet(
#         dataTier = cms.untracked.string('GEN-SIM-RECO'),
#         filterName = cms.untracked.string('')
#     ),
#     fileName = cms.untracked.string('file:step3.root'),
#     outputCommands = process.RECOSIMEventContent.outputCommands,
#     splitLevel = cms.untracked.int32(0)
# )

# process.MINIAODSIMoutput = cms.OutputModule("PoolOutputModule",
#     compressionAlgorithm = cms.untracked.string('LZMA'),
#     compressionLevel = cms.untracked.int32(4),
#     dataset = cms.untracked.PSet(
#         dataTier = cms.untracked.string('MINIAODSIM'),
#         filterName = cms.untracked.string('')
#     ),
#     dropMetaData = cms.untracked.string('ALL'),
#     eventAutoFlushCompressedSize = cms.untracked.int32(-900),
#     fastCloning = cms.untracked.bool(False),
#     fileName = cms.untracked.string('file:step3_inMINIAODSIM.root'),
#     outputCommands = process.MINIAODSIMEventContent.outputCommands,
#     overrideBranchesSplitLevel = cms.untracked.VPSet(
#         cms.untracked.PSet(
#             branch = cms.untracked.string('patPackedCandidates_packedPFCandidates__*'),
#             splitLevel = cms.untracked.int32(99)
#         ), 
#         cms.untracked.PSet(
#             branch = cms.untracked.string('recoGenParticles_prunedGenParticles__*'),
#             splitLevel = cms.untracked.int32(99)
#         ), 
#         cms.untracked.PSet(
#             branch = cms.untracked.string('patTriggerObjectStandAlones_slimmedPatTrigger__*'),
#             splitLevel = cms.untracked.int32(99)
#         ), 
#         cms.untracked.PSet(
#             branch = cms.untracked.string('patPackedGenParticles_packedGenParticles__*'),
#             splitLevel = cms.untracked.int32(99)
#         ), 
#         cms.untracked.PSet(
#             branch = cms.untracked.string('patJets_slimmedJets__*'),
#             splitLevel = cms.untracked.int32(99)
#         ), 
#         cms.untracked.PSet(
#             branch = cms.untracked.string('recoVertexs_offlineSlimmedPrimaryVertices__*'),
#             splitLevel = cms.untracked.int32(99)
#         ), 
#         cms.untracked.PSet(
#             branch = cms.untracked.string('recoCaloClusters_reducedEgamma_reducedESClusters_*'),
#             splitLevel = cms.untracked.int32(99)
#         ), 
#         cms.untracked.PSet(
#             branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedEBRecHits_*'),
#             splitLevel = cms.untracked.int32(99)
#         ), 
#         cms.untracked.PSet(
#             branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedEERecHits_*'),
#             splitLevel = cms.untracked.int32(99)
#         ), 
#         cms.untracked.PSet(
#             branch = cms.untracked.string('recoGenJets_slimmedGenJets__*'),
#             splitLevel = cms.untracked.int32(99)
#         ), 
#         cms.untracked.PSet(
#             branch = cms.untracked.string('patJets_slimmedJetsPuppi__*'),
#             splitLevel = cms.untracked.int32(99)
#         ), 
#         cms.untracked.PSet(
#             branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedESRecHits_*'),
#             splitLevel = cms.untracked.int32(99)
#         )
#     ),
#     overrideInputFileSplitLevels = cms.untracked.bool(True),
#     splitLevel = cms.untracked.int32(0)
# )

# process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
#     dataset = cms.untracked.PSet(
#         dataTier = cms.untracked.string('DQMIO'),
#         filterName = cms.untracked.string('')
#     ),
#     fileName = cms.untracked.string('file:step3_inDQM.root'),
#     outputCommands = process.DQMEventContent.outputCommands,
#     splitLevel = cms.untracked.int32(0)
# )

# Additional output definition

# # Other statements
# process.mix.playback = True
# process.mix.digitizers = cms.PSet()
# for a in process.aliases: delattr(process, a)
# process.RandomNumberGeneratorService.restoreStateLabel=cms.untracked.string("randomEngineStateProducer")

from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run1_mc', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '100X_upgrade2018_realistic_v10', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '105X_mcRun2_asymptotic_v2', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '120X_mcRun3_2021_realistic_v6', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '126X_mcRun3_2023_forPU65_v1', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '126X_mcRun3_2023_forPU65_v4', '')

process.pfChargedHadronAnalyzer = cms.EDAnalyzer(
    "PFChargedHadronAnalyzer",
    PFCandidates = cms.InputTag("particleFlow"),
    PFSimParticles = cms.InputTag("particleFlowSimParticle"),
    EcalPFClusters = cms.InputTag("particleFlowClusterECAL"),
    HcalPFClusters = cms.InputTag("particleFlowClusterHCAL"),
    ptMin = cms.double(1.),                     # Minimum pt                                                                         
    pMin = cms.double(1.),                      # Minimum p                                                                          
    nPixMin = cms.int32(2),                     # Nb of pixel hits                                                                   
    nHitMin = cms.vint32(14,17,20,17,10),       # Nb of track hits                                                                   
    nEtaMin = cms.vdouble(1.4,1.6,2.0,2.4,2.6), # in these eta ranges                                                                
    hcalMin = cms.double(0.5),                   # Minimum hcal energy                                                               
    ecalMax = cms.double(1E9),                  # Maximum ecal energy                                                                
    verbose = cms.untracked.bool(True),         # not used.                                                                          
    #rootOutputFile = cms.string("PGun__2_200GeV__81X_upgrade2017_realistic_v22.root"),# the root tree                               
    rootOutputFile = cms.string("step3.root"),# the root tree                                                       
#    IsMinBias = cms.untracked.bool(False)                                                                                           
)




process.load("RecoParticleFlow.PFProducer.particleFlowSimParticle_cff")
#process.load("RecoParticleFlow.Configuration.HepMCCopy_cfi")                                                                        

process.particleFlowSimParticle.ParticleFilter = cms.PSet(
        # Allow *ALL* protons with energy > protonEMin                                                                               
        protonEMin = cms.double(5000.0),
        # Particles must have abs(eta) < etaMax (if close enough to 0,0,0)                                                           
        etaMax = cms.double(5.3),
        # Charged particles with pT < pTMin (GeV/c) are not simulated                                                                
        chargedPtMin = cms.double(0.0),
        # Particles must have energy greater than EMin [GeV]                                                                         
        EMin = cms.double(0.0),
        # half-length of the ECAL endcap inner surface                                                                                                                                                      
        rMax = cms.double(129.),
        zMax = cms.double(317.),
        # List of invisible particles (abs of pdgid)                                                                                                                                                        
        invisibleParticles = cms.vint32()
)
process.genReReco = cms.Sequence(#process.generator+                                                                                 
                                 #process.genParticles+                                                                              
                                 #process.genJetParticles+                                                                           
                                 #process.recoGenJets+                                                                               
                                 #process.genMETParticles+                                                                           
                                 #process.recoGenMET+                                                                                
                                 process.particleFlowSimParticle)



# Path and EndPath definitions


process.EDA = cms.EndPath(process.pfChargedHadronAnalyzer)
process.gRR = cms.EndPath(process.genReReco)



# process.raw2digi_step = cms.Path(process.RawToDigi)
# process.L1Reco_step = cms.Path(process.L1Reco)
# process.reconstruction_step = cms.Path(process.reconstruction)
# process.recosim_step = cms.Path(process.recosim)
# process.eventinterpretaion_step = cms.Path(process.EIsequence)
# process.EDA = cms.EndPath(process.pfChargedHadronAnalyzer)
# process.gRR = cms.EndPath(process.genReReco)




# process.Flag_trackingFailureFilter = cms.Path(process.goodVertices+process.trackingFailureFilter)
# process.Flag_goodVertices = cms.Path(process.primaryVertexFilter)
# process.Flag_CSCTightHaloFilter = cms.Path(process.CSCTightHaloFilter)
# process.Flag_trkPOGFilters = cms.Path(process.trkPOGFilters)
# process.Flag_HcalStripHaloFilter = cms.Path(process.HcalStripHaloFilter)
# process.Flag_trkPOG_logErrorTooManyClusters = cms.Path(~process.logErrorTooManyClusters)
# process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
# process.Flag_ecalLaserCorrFilter = cms.Path(process.ecalLaserCorrFilter)
# process.Flag_globalSuperTightHalo2016Filter = cms.Path(process.globalSuperTightHalo2016Filter)
# process.Flag_eeBadScFilter = cms.Path(process.eeBadScFilter)
# process.Flag_METFilters = cms.Path(process.metFilters)
# process.Flag_chargedHadronTrackResolutionFilter = cms.Path(process.chargedHadronTrackResolutionFilter)
# process.Flag_globalTightHalo2016Filter = cms.Path(process.globalTightHalo2016Filter)
# process.Flag_CSCTightHaloTrkMuUnvetoFilter = cms.Path(process.CSCTightHaloTrkMuUnvetoFilter)
# process.Flag_HBHENoiseIsoFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseIsoFilter)
# process.Flag_BadChargedCandidateSummer16Filter = cms.Path(process.BadChargedCandidateSummer16Filter)
# process.Flag_hcalLaserEventFilter = cms.Path(process.hcalLaserEventFilter)
# process.Flag_BadPFMuonFilter = cms.Path(process.BadPFMuonFilter)
# process.Flag_ecalBadCalibFilter = cms.Path()
# process.Flag_HBHENoiseFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseFilter)
# process.Flag_trkPOG_toomanystripclus53X = cms.Path(~process.toomanystripclus53X)
# process.Flag_EcalDeadCellBoundaryEnergyFilter = cms.Path(process.EcalDeadCellBoundaryEnergyFilter)
# process.Flag_BadChargedCandidateFilter = cms.Path(process.BadChargedCandidateFilter)
# process.Flag_trkPOG_manystripclus53X = cms.Path(~process.manystripclus53X)
# process.Flag_BadPFMuonSummer16Filter = cms.Path(process.BadPFMuonSummer16Filter)
# process.Flag_muonBadTrackFilter = cms.Path(process.muonBadTrackFilter)
# process.Flag_CSCTightHalo2015Filter = cms.Path(process.CSCTightHalo2015Filter)
# process.prevalidation_step = cms.Path(process.prevalidationNoHLT)
# process.prevalidation_step1 = cms.Path(process.prevalidationMiniAOD)
# process.validation_step = cms.EndPath(process.validationNoHLT)
# process.validation_step1 = cms.EndPath(process.validationMiniAOD)
# process.dqmoffline_step = cms.EndPath(process.DQMOfflineFakeHLT)
# process.dqmoffline_1_step = cms.EndPath(process.DQMOfflineMiniAOD)
# process.dqmofflineOnPAT_step = cms.EndPath(process.PostDQMOffline)
# process.dqmofflineOnPAT_1_step = cms.EndPath(process.PostDQMOfflineMiniAOD)
# process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)
# process.MINIAODSIMoutput_step = cms.EndPath(process.MINIAODSIMoutput)
# process.DQMoutput_step = cms.EndPath(process.DQMoutput)

# # Schedule definition
# #process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.recosim_step,process.eventinterpretaion_step,process.Flag_HBHENoiseFilter,process.Flag_HBHENoiseIsoFilter,process.Flag_CSCTightHaloFilter,process.Flag_CSCTightHaloTrkMuUnvetoFilter,process.Flag_CSCTightHalo2015Filter,process.Flag_globalTightHalo2016Filter,process.Flag_globalSuperTightHalo2016Filter,process.Flag_HcalStripHaloFilter,process.Flag_hcalLaserEventFilter,process.Flag_EcalDeadCellTriggerPrimitiveFilter,process.Flag_EcalDeadCellBoundaryEnergyFilter,process.Flag_ecalBadCalibFilter,process.Flag_goodVertices,process.Flag_eeBadScFilter,process.Flag_ecalLaserCorrFilter,process.Flag_trkPOGFilters,process.Flag_chargedHadronTrackResolutionFilter,process.Flag_muonBadTrackFilter,process.Flag_BadChargedCandidateFilter,process.Flag_BadPFMuonFilter,process.Flag_BadChargedCandidateSummer16Filter,process.Flag_BadPFMuonSummer16Filter,process.Flag_trkPOG_manystripclus53X,process.Flag_trkPOG_toomanystripclus53X,process.Flag_trkPOG_logErrorTooManyClusters,process.Flag_METFilters,process.prevalidation_step,process.prevalidation_step1,process.validation_step,process.validation_step1,process.dqmoffline_step,process.dqmoffline_1_step,process.dqmofflineOnPAT_step,process.dqmofflineOnPAT_1_step,process.RECOSIMoutput_step,process.MINIAODSIMoutput_step,process.DQMoutput_step)

# process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.eventinterpretaion_step,process.gRR,process.EDA)


# process.schedule = cms.Schedule(process.gRR,process.EDA)



# process.schedule.associate(process.patTask)



# from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
# associatePatAlgosToolsTask(process)

# # customisation of the process.

# # Automatic addition of the customisation function from SimGeneral.MixingModule.fullMixCustomize_cff
# from SimGeneral.MixingModule.fullMixCustomize_cff import setCrossingFrameOn 

# #call to customisation function setCrossingFrameOn imported from SimGeneral.MixingModule.fullMixCustomize_cff
# process = setCrossingFrameOn(process)

# # End of customisation functions
# #do not add changes to your config after this point (unless you know what you are doing)
# from FWCore.ParameterSet.Utilities import convertToUnscheduled
# process=convertToUnscheduled(process)

# # customisation of the process.

# # Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
# from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllMC 

# #call to customisation function miniAOD_customizeAllMC imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
# process = miniAOD_customizeAllMC(process)

# # End of customisation functions

# # Customisation from command line

# #Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
# from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
# process = customiseLogErrorHarvesterUsingOutputCommands(process)

# # Add early deletion of temporary data products to reduce peak memory need
# from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
# process = customiseEarlyDelete(process)
# # End adding early deletion
