import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("NEWGEN")

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
    fileNames = cms.untracked.vstring('file:test/output_stdgen.root'),
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
    fileName = cms.untracked.string('test/output_newgen.root'),
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
        "drop *_genMetTrue_*_GEN",
        "drop *_genParticlesForJets_*_NEWGEN",
        "drop *_genParticlesForJetsNoNu_*_NEWGEN",
        "drop *_ak4GenJets*_rho_NEWGEN",
        "drop *_ak4GenJets*_rhos_NEWGEN",
        "drop *_ak4GenJets*_sigma_NEWGEN",
        "drop *_ak4GenJets*_sigmas_NEWGEN",
        "drop *_ak8GenJets*_rho_NEWGEN",
        "drop *_ak8GenJets*_rhos_NEWGEN",
        "drop *_ak8GenJets*_sigma_NEWGEN",
        "drop *_ak8GenJets*_sigmas_NEWGEN",
        "drop *_genCandidatesForMET_*_NEWGEN",
        "drop *_genParticlesForMETAllVisible_*_NEWGEN",
    ),
    splitLevel = cms.untracked.int32(0)
)

# Redefine sequences to avoid collection conflicts
process.load("Configuration.StandardSequences.Generator_cff")
process.GeneInfoTask.remove(process.genParticles)

process.load('GenParticlesHacker.Producers.NewGenParticlesProducer_cfi')
process.GeneInfoTask.add(process.genParticles)

process.load("RecoJets.Configuration.RecoGenJets_cff")
process.genParticlesForJets.src = cms.InputTag("genParticles","modified","NEWGEN")
process.ak4GenJets.src = cms.InputTag("genParticlesForJets","","NEWGEN")
process.ak8GenJets.src = cms.InputTag("genParticlesForJets","","NEWGEN")
process.genParticlesForJetsNoNu.src = cms.InputTag("genParticles","modified","NEWGEN")
process.ak4GenJetsNoNu.src = cms.InputTag("genParticlesForJetsNoNu","","NEWGEN")
process.ak8GenJetsNoNu.src = cms.InputTag("genParticlesForJetsNoNu","","NEWGEN")
process.genJetSequence = cms.Sequence(
    process.genParticlesForJets +
    process.genParticlesForJetsNoNu +
    process.ak4GenJets +
    process.ak4GenJetsNoNu +
    process.ak8GenJets +
    process.ak8GenJetsNoNu
)

process.load("RecoMET.Configuration.RecoGenMET_cff")
process.genCandidatesForMET.src = cms.InputTag("genParticles","modified","NEWGEN")
process.genParticlesForMETAllVisible.src = cms.InputTag("genParticles","modified","NEWGEN")
process.genMetCalo.src = cms.InputTag("genCandidatesForMET","","NEWGEN")
process.genMetTrue.src = cms.InputTag("genParticlesForMETAllVisible","","NEWGEN")
process.genMETSequence = cms.Sequence(
    process.genCandidatesForMET +
    process.genParticlesForMETAllVisible +
    process.genMetCalo +
    process.genMetTrue
)

process.producer_step = cms.Path(process.genParticles)
process.genjet_step = cms.Path(process.genJetSequence)
process.genmet_step = cms.Path(process.genMETSequence)
process.endjob_step = cms.EndPath(process.GENoutput)
process.schedule = cms.Schedule(process.producer_step,process.genjet_step,process.genmet_step,process.endjob_step)
#process.schedule = cms.Schedule(process.producer_step,process.endjob_step)


