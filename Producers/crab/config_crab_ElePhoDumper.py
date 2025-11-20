from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config
import FWCore.ParameterSet.Config as cms

config = Configuration()

config.section_('General')
#config.General.requestName       = 'DoubleElectron_Pt-1To300-gun_RunIISummer20UL18MiniAODv2-FlatPU0to70EdalIdealGT_EdalIdealGT'
#config.General.requestName       = 'DoublePhoton_Pt-5To300-gun_RunIISummer20UL18MiniAODv2-FlatPU0to70EdalIdealGT_EdalIdealGT'
#config.General.requestName       = 'DYJetsToEE_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8-stdGEN_Dumper'
#config.General.requestName       = 'DYJetsToGammaGamma_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8-newGEN_Dumper'
#config.General.requestName       = 'ManyElectron_Pt-1To500-gun-stdGEN_Dumper'
config.General.requestName       = 'ManyPhoton_Pt-1To500-gun-newGEN_Dumper'
#config.General.requestName       = 'GluGluHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8_Dumper'
#config.General.requestName       = 'DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_CENTRAL_Dumper'
config.General.transferLogs      = True
config.General.transferOutputs   = True

config.section_('JobType')
config.JobType.pluginName        = 'Analysis'
# Name of the CMSSW configuration file
#config.JobType.psetName          = 'ElePhoDumper_ele_cfg.py'
config.JobType.psetName          = 'ElePhoDumper_pho_cfg.py'
config.JobType.allowUndistributedCMSSW = True

config.section_('Data')
#config.Data.inputDataset         = '/DoubleElectron_Pt-1To300-gun/RunIISummer20UL18MiniAODv2-FlatPU0to70EdalIdealGT_EdalIdealGT_106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM'
#config.Data.inputDataset         = '/DoublePhoton_Pt-5To300-gun/RunIISummer20UL18MiniAODv2-FlatPU0to70EdalIdealGT_EdalIdealGT_106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM'
#config.Data.inputDataset         = '/DYJetsToEE_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/bmarzocc-RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v4-MINIAODSIM_from_stdGEN-v2-5994cd161b8b5688239399fbfeba9eac/USER'
#config.Data.inputDataset         = '/DYJetsToGammaGamma_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/fiemmi-RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v4-MINIAODSIM_from_newGEN_preserveAngles-5994cd161b8b5688239399fbfeba9eac/USER'
#config.Data.inputDataset         = '/ManyElectron_Pt-1To500-gun/bmarzocc-RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-MINIAODSIM_from_stdGEN-5994cd161b8b5688239399fbfeba9eac/USER'
config.Data.inputDataset         = '/ManyPhoton_Pt-1To500-gun/bmarzocc-RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-MINIAODSIM_from_newGEN-v2-5994cd161b8b5688239399fbfeba9eac/USER'
#config.Data.inputDataset         = '/GluGluHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM'
#config.Data.inputDataset         = '/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM'
config.Data.outputDatasetTag     = 'Run2018'
#config.Data.inputDBS             = 'global'
config.Data.inputDBS             = 'phys03'
config.Data.splitting            = 'FileBased' #'LumiBased'
config.Data.unitsPerJob          = 1 #30000
config.Data.outLFNDirBase        = '/store/group/phys_egamma/bmarzocc/ParticleGuns/'
config.Data.publication          = False
#config.Data.lumiMask             = 'Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt'

# GRID
config.section_('Site')
config.Site.storageSite   =  'T2_CH_CERN'
config.Site.ignoreGlobalBlacklist = True
