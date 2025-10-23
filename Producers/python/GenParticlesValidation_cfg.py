import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.parseArguments()

process = cms.Process("genAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(100))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    skipEvents = cms.untracked.uint32(98)
)

process.genParticleAnalyzer = cms.EDAnalyzer('GenParticleAnalyzer',
    #genParticles    = cms.InputTag("prunedGenParticles",""), 
    #genParticles    = cms.InputTag("genParticles","","GEN"), 
    genParticles    = cms.InputTag("genParticles","","NEWGEN"), 
    #genParticles    = cms.InputTag("genParticles","","DIGI2RAW"), 
    lheEvents       = cms.InputTag("externalLHEProducer"), 
    
    debug           = cms.bool(False),       #use DEBUG mode
    motherPdgId     = cms.vint32(-999),     #PdgId of the mother
    motherStatus    = cms.vint32(-999),     #Status of the mother
    excludeMomPdgId = cms.vint32(111,321,331,221),     #PdgId of the mother
    pdgIdBeforeFSR  = cms.vint32(11,-11),   #Input pdgIds
    #pdgIdIn         = cms.vint32(11,-11),   #Input pdgIds
    pdgIdIn         = cms.vint32(22),   #Input pdgIds
    #statusIn        = cms.vint32(1,1),      #Input status
    statusIn        = cms.vint32(1),      #Input status
    ignoreTauDecays = cms.bool(True),       #Ignore particles from tau-decays
    doLHEMatching   = cms.bool(True),       #Match with LHEparticles
    dRMaxLHEMatch   = cms.double(0.8)       #dRMax for LHEparticles-matching
    #doLHEMatching   = cms.bool(False),       #Match with LHEparticles
    #dRMaxLHEMatch   = cms.double(999.)       #dRMax for LHEparticles-matching
)
process.p = cms.Path(process.genParticleAnalyzer)
