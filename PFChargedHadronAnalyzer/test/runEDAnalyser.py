# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step2 --mc --eventcontent RECOSIM --datatier GEN-SIM-RECO --conditions 120X_mcRun3_2021_realistic_v6 --step RAW2DIGI,L1Reco,RECO,RECOSIM --nThreads 4 --geometry DB:Extended --era Run3 --filein file:JME-Run3Summer21DRPremix-00001_1_1.root --fileout file:JME-Run3Summer21DRPremix-00001_2_1.root -n -1
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_cff import Run3

process = cms.Process('RECO',Run3)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.RecoSim_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000), #NEvents
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
      #fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/Run3Winter23Digi/SinglePionGun_E200to500/GEN-SIM-RAW/NoPUGTv4_126X_mcRun3_2023_forPU65_v4-v2/50000/00161083-8bd0-4a00-b462-c5e5141ab85c.root'),
    fileNames = cms.untracked.vstring('file:step3_reco.root'), # input root file name
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    makeTriggerResults = cms.obsolete.untracked.bool,
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step2 nevts:-1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
#    fileName = cms.untracked.string('file:JME-Run3Summer21DRPremix-00001_2_1.root'),
    fileName = cms.untracked.string('file:step3_reco.root'),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '120X_mcRun3_2021_realistic_v6', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '126X_mcRun3_2023_forPU65_v4', '')

process.pfChargedHadronAnalyzer = cms.EDAnalyzer(
    "PFChargedHadronAnalyzer",
    PFCandidates = cms.InputTag("particleFlow"),
    PFSimParticles = cms.InputTag("particleFlowSimParticle"),
    EcalPFClusters = cms.InputTag("particleFlowClusterECAL"),
    HcalPFClusters = cms.InputTag("particleFlowClusterHCAL"),
    EcalPFrechit = cms.InputTag("particleFlowRecHitECAL"),
    HcalPFrechit = cms.InputTag("particleFlowRecHitHBHE"),
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
        rMax = cms.double(129.),
        # half-length of the ECAL endcap inner surface 
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
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.recosim_step = cms.Path(process.recosim)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)
process.EDA = cms.EndPath(process.pfChargedHadronAnalyzer)
process.gRR = cms.EndPath(process.genReReco)

# Schedule definition
process.schedule = cms.Schedule(process.EDA)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
process.options.numberOfThreads = 1
process.options.numberOfStreams = 0
process.options.numberOfConcurrentLuminosityBlocks = 2
process.options.eventSetup.numberOfConcurrentIOVs = 1
if hasattr(process, 'DQMStore'): process.DQMStore.assertLegacySafe=cms.untracked.bool(False)



# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
