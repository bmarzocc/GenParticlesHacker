import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("MODIFY")

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic25ns13TeVEarly2018Collision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag,'106X_upgrade2018_realistic_v4','')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( -1 ) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1 )
                                                                       
process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),                       
    fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/user/bmarzocc/ManyElectron_Pt-1To500-gun/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-MINIAODSIM_from_stdGEN/251108_225359/0001/step6_stdMINIAOD_1775.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.GENoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(1),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(20971520),
    fileName = cms.untracked.string('test/step1_newGEN.root'),
    outputCommands = cms.untracked.vstring(
        "keep *",
        "drop ints_genParticles__GEN",
        "drop recoGenParticles_genParticles__GEN",
        "drop *_TriggerResults_*_GEN",
        "drop *_ak4GenJets_*_GEN",
        "drop *_ak4GenJetsNoNu_*_GEN",
        "drop *_ak8GenJets_*_GEN",
        "drop *_ak8GenJetsNoNu_*_GEN",
        "drop *_genMetCalo_*_GEN",
        "drop *_genMetTrue_*_GEN"
    ),
    splitLevel = cms.untracked.int32(0)
)

# Redefine sequences to avoid collection conflicts
process.load("Configuration.StandardSequences.Generator_cff")
process.GeneInfoTask.remove(process.genParticles)
process.genParticles = cms.EDProducer("NewGenParticlesProducer",

    genInfo         = cms.InputTag("generator"),
    lheEvents       = cms.InputTag("externalLHEProducer"), 
    genInts         = cms.InputTag("genParticles","","GEN"),
    genCollection   = cms.InputTag("genParticles","","GEN"),
    
    motherPdgId     = cms.vint32(-999,-999), #PdgId of the mother
    motherStatus    = cms.vint32(-999,-999), #Status of the mother
    excludeMomPdgId = cms.vint32(-999),      #PdgId of the mother
    pdgIdIn         = cms.vint32(11,-11),    #Input pdgIds to change
    pdgIdOut        = cms.vint32(22, 22),    #Output new pdgIds
    statusIn        = cms.vint32(1,1),       #Input status to change
    statusOut       = cms.vint32(1,1),       #Output new status   
    chargeOut       = cms.vint32(0,0),       #Output new charge  
    massOut         = cms.vdouble(0.,0.),    #Output new mass   
    doLHEMatching   = cms.bool(False),        #Match with LHEparticles
    dRMaxLHEMatch   = cms.double(0.8),       #dRMax for LHEparticles-matching
    ignoreTauDecays = cms.bool(True),        #Ignore particles from tau-decays
    suppressFSR     = cms.bool(True),        #Suppress FSR
    preserveAngles  = cms.bool(True),        #If Suppress FSR, preserve the particle angles, therefore the mother invariante mass will change
    debug           = cms.bool(True)         #Use debug mode
    
)
process.genParticles = cms.EDProducer("NewGenParticlesProducer",

    genInts         = cms.InputTag("genParticles","","GEN"),
    genCollection   = cms.InputTag("genParticles","","GEN"),
    lheEvents       = cms.InputTag("externalLHEProducer"), 
    
    motherPdgId     = cms.vint32(-999,-999), #PdgId of the mother
    motherStatus    = cms.vint32(-999,-999), #Status of the mother
    excludeMomPdgId = cms.vint32(-999),      #PdgId of the mother
    pdgIdIn         = cms.vint32(11,-11),    #Input pdgIds to change
    pdgIdOut        = cms.vint32(22, 22),    #Output new pdgIds
    statusIn        = cms.vint32(1,1),       #Input status to change
    statusOut       = cms.vint32(1,1),       #Output new status   
    chargeOut       = cms.vint32(0,0),       #Output new charge  
    massOut         = cms.vdouble(0.,0.),    #Output new mass   
    doLHEMatching   = cms.bool(True),        #Match with LHEparticles
    dRMaxLHEMatch   = cms.double(0.8),       #dRMax for LHEparticles-matching
    ignoreTauDecays = cms.bool(True),        #Ignore particles from tau-decays
    suppressFSR     = cms.bool(True),        #Suppress FSR
    preserveAngles  = cms.bool(True),        #If Suppress FSR, preserve the particle angles, therefore the mother invariante mass will change
    debug           = cms.bool(False)         #Use debug mode
    
)
process.GeneInfoTask.add(process.genParticles)
process.producer_step = cms.Path(process.genParticles)
process.endjob_step = cms.EndPath(process.GENoutput)
process.schedule = cms.Schedule(process.producer_step,process.endjob_step)

