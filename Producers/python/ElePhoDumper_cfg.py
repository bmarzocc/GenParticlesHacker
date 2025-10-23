import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("ElePhoDumper")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff") # gives deprecated message in 80X but still runs
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag,'106X_upgrade2018_realistic_v16_L1v1','')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( -1 ) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1 )
                                                                       
process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),                       
    fileNames = cms.untracked.vstring("file:crab/step6_stdMINIAOD_19.root"),
    secondaryFileNames = cms.untracked.vstring()
    ) 

process.elephodumper = cms.EDAnalyzer("ElePhoDumper",

    rhoCollection                   = cms.InputTag("fixedGridRhoAll",""),
    pileupSummary                   = cms.InputTag("slimmedAddPileupInfo",""),
    vertexCollection                = cms.InputTag("offlineSlimmedPrimaryVertices",""),
    #genParticleCollection           = cms.InputTag("genParticles","","NEWGEN"),
    genParticleCollection           = cms.InputTag("prunedGenParticles"),
    ebRechitCollection              = cms.InputTag("reducedEgamma","reducedEBRecHits"),
    eeRechitCollection              = cms.InputTag("reducedEgamma","reducedEERecHits"),
    patElectronCollection           = cms.InputTag("slimmedElectrons",""), 
    patPhotonCollection             = cms.InputTag("slimmedPhotons",""), 
    
    isMC                            = cms.bool(True),     #isMC
    motherPdgId                     = cms.vint32(-999),      #PdgId of the mother
    motherStatus                    = cms.vint32(-999),    #Status of the mother
    pdgIdBeforeFSR                  = cms.vint32(11,-11), #genMatch pdgIds before FSR
    pdgId                           = cms.vint32(11,-11), #genMatch pdgIds
    #pdgId                           = cms.vint32(22),     #genMatch pdgIds
    status                          = cms.vint32(1,1),    #genMatch status
    #status                          = cms.vint32(1),    #genMatch status
    ignoreTauDecays                 = cms.bool(True),     #ignore gen-particles from tau-decays
    dRMin                           = cms.double(0.1),    #minDR for gen matching
   
    savePhotons                     = cms.bool(False),    #save patPhotons and patMET information
    saveElectrons                   = cms.bool(True),     #save patElectrons and patMET information
    
    egmCutBasedElectronIDVeto       = cms.string('cutBasedElectronID-Fall17-94X-V2-veto'),   #cutBasedEleID veto
    egmCutBasedElectronIDloose      = cms.string('cutBasedElectronID-Fall17-94X-V2-loose'),  #cutBasedEleID loose  
    egmCutBasedElectronIDmedium     = cms.string('cutBasedElectronID-Fall17-94X-V2-medium'), #cutBasedEleID medium 
    egmCutBasedElectronIDtight      = cms.string('cutBasedElectronID-Fall17-94X-V2-tight'),  #cutBasedEleID tight
    egmMVAElectronIDloose           = cms.string('mvaEleID-Fall17-iso-V2-wpLoose'),          #mvaEleID loose 
    egmMVAElectronIDmedium          = cms.string('mvaEleID-Fall17-iso-V2-wp90'),             #mvaEleID medium 
    egmMVAElectronIDtight           = cms.string('mvaEleID-Fall17-iso-V2-wp80'),             #mvaEleID tight  
    egmMVAElectronIDlooseNoIso      = cms.string('mvaEleID-Fall17-noIso-V2-wpLoose'),        #mvaEleIDNoIso loose 
    egmMVAElectronIDmediumNoIso     = cms.string('mvaEleID-Fall17-noIso-V2-wp90'),           #mvaEleIDNoIso medium 
    egmMVAElectronIDtightNoIso      = cms.string('mvaEleID-Fall17-noIso-V2-wp80'),           #mvaEleIDNoIso tight  
    heepElectronID                  = cms.string('heepElectronID-HEEPV70'),                  #mvaEleIDNoIso tight
    egmCutBasedPhotonIDloose        = cms.string('cutBasedPhotonID-Fall17-94X-V2-loose'),    #cutBasedPhoID loose  
    egmCutBasedPhotonIDmedium       = cms.string('cutBasedPhotonID-Fall17-94X-V2-medium'),   #cutBasedPhoID medium 
    egmCutBasedPhotonIDtight        = cms.string('cutBasedPhotonID-Fall17-94X-V2-tight'),    #cutBasedPhoID tight 
    egmMVAPhotonIDmedium            = cms.string('mvaPhoID-RunIIFall17-v2-wp90'),            #mvaPhoID medium 
    egmMVAPhotonIDtight             = cms.string('mvaPhoID-RunIIFall17-v2-wp80'),            #mvaPhoID tight     
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("test/output_elephoDumper_updated.root")
)

process.p = cms.Path(
    process.elephodumper
)
