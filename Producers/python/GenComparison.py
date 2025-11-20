import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.ioFilesOpt = cms.PSet(

    # #maxEvents
    skipEvents = cms.int32(-1),
    maxEvents = cms.int32(-1),
    
    ## genParticle-collection
    hepMC = cms.string("generatorSmeared"),
    #genCollection = cms.string("genParticles"),
    genCollection = cms.string("prunedGenParticles"),
    
    ## debug
    debug = cms.bool(True),
    
    ## isParticleGun
    isParticleGun = cms.bool(True),
    
    ## checkHepMC
    checkHepMC = cms.bool(False),
    
    ## doLHEMatching
    doLHEMatching = cms.bool(False),
    dRMaxLHEMatch = cms.double(0.8),
    dEFracMaxLHEMatch = cms.double(999.),

    
    ##input files
    phoInputFiles = cms.vstring('file:crab/step6_stdMINIAOD_21.root'),
    eleInputFiles = cms.vstring('root://cms-xrd-global.cern.ch//store/user/bmarzocc/ManyElectron_Pt-1To500-gun/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-MINIAODSIM_from_stdGEN/251108_225359/0000/step6_stdMINIAOD_21.root')
    
)
