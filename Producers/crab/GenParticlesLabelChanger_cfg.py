import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('standard')
options.register('inputFile', 
                  'step1_newGEN_tmp.root', 
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,
                  "Input file")
options.register('outputFile', 
                  'step1_newGEN.root', 
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,
                  "Output file")
options.parseArguments()
print(options)

process = cms.Process("NEWGEN")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("file:"+options.inputFile)
)

process.GENoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(1),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(20971520),
    fileName = cms.untracked.string(options.outputFile),
    outputCommands = cms.untracked.vstring(
        "keep *",
        "drop *_genParticles_*_MODIFY",
        "drop *_TriggerResults_*_MODIFY",
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
        "drop *_genParticlesForMETAllVisible_*_NEWGEN"
    ),
    splitLevel = cms.untracked.int32(0)
)


process.genParticles = cms.EDProducer('GenParticlesLabelChanger',
    genInts         = cms.InputTag("genParticles","modified","MODIFY"),
    genCollection   = cms.InputTag("genParticles","modified","MODIFY")
)

process.load("RecoJets.Configuration.GenJetParticles_cff")
process.load("RecoJets.Configuration.RecoGenJets_cff")
process.genParticlesForJets.src = cms.InputTag("genParticles","","NEWGEN")
process.ak4GenJets.src = cms.InputTag("genParticlesForJets","","NEWGEN")
process.ak8GenJets.src = cms.InputTag("genParticlesForJets","","NEWGEN")
process.genParticlesForJetsNoNu.src = cms.InputTag("genParticles","","NEWGEN")
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

process.load("RecoMET.Configuration.GenMETParticles_cff")
process.load("RecoMET.Configuration.RecoGenMET_cff")
process.genCandidatesForMET.src = cms.InputTag("genParticles","","NEWGEN")
process.genParticlesForMETAllVisible.src = cms.InputTag("genParticles","","NEWGEN")
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

