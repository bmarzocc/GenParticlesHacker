import FWCore.ParameterSet.Config as cms

genParticles = cms.EDProducer("NewGenParticlesProducer",

    genInts         = cms.InputTag("genParticles","","GEN"),
    genCollection   = cms.InputTag("genParticles","","GEN"),
    
    motherPdgId     = cms.int32(23),      #PdgId of the mother
    motherStatus    = cms.int32(-999),    #Status of the mother
    pdgIdIn         = cms.vint32(11,-11), #Input pdgIds to change
    pdgIdOut        = cms.vint32(22, 22), #Output new pdgIds
    statusIn        = cms.vint32(1,1),    #Input status to change
    statusOut       = cms.vint32(1,1),    #Output new status   
    chargeOut       = cms.vint32(0,0),    #Output new charge  
    massOut         = cms.vdouble(0.,0.), #Output new mass   
    ignoreTauDecays = cms.bool(True),     #Ignore particles from tau-decays
    suppressFSR     = cms.bool(True),     #Suppress FSR
    preserveAngles  = cms.bool(True)      #If Suppress FSR, preserve the particle angles, therefore the mother invariante mass will change
    
)
