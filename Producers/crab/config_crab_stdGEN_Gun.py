from CRABClient.UserUtilities import config

config = config()

## General settings
config.General.requestName ='ManyElectron_Pt-15To250-gun_stdGEN'
config.General.transferOutputs = True
config.General.transferLogs = True

## PrivateMC type with a fake miniAOD step to circunvent crab requests (official data-tier for PrivateMC)
config.JobType.pluginName  = 'PrivateMC'
config.JobType.psetName    = 'ManyElectron_Pt-15To250-gun-RunIISummer20UL18MiniAODv2_stdGEN_cfg.py'
## To be executed on node with Arguments
config.JobType.scriptExe   = 'scriptExe_stdGEN_GUN.sh'
config.JobType.inputFiles  = ['scriptExe_stdGEN_GUN.sh','ManyElectron_Pt-15To250-gun-RunIISummer20UL18MiniAODv2_stdGEN_cfg.py']
## Output file to be collected
config.JobType.outputFiles = ["step1_stdGEN.root"]
config.JobType.disableAutomaticOutputCollection = True
## Memory, cores, cmssw
#config.JobType.allowUndistributedCMSSW = True
#config.JobType.maxMemoryMB = 2500
#config.JobType.numCores    = 1

## Data
config.Data.splitting   = 'EventBased'
config.Data.unitsPerJob = 10000
config.Data.totalUnits  = 100000000

config.Data.outputPrimaryDataset = 'ManyElectron_Pt-15To250-gun'
config.Data.publication          = True
config.Data.publishDBS           = 'phys03'
#config.Data.ignoreLocality       = True
config.Data.outputDatasetTag     = 'RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v4-stdGEN'
config.Data.outLFNDirBase        =  '/store/user/bmarzocc/'
#config.Data.publishWithGroupName = True

## Site
config.Site.storageSite         = 'T2_US_Wisconsin'
#config.Site.whitelist           = ['T2_CH_CERN']
