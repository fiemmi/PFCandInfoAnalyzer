from CRABClient.UserUtilities import config

config = config()
config.General.requestName = 'flatTree_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8'
config.General.workArea = 'crab_Sep03_2020_h1210'
config.General.transferOutputs = True
config.General.transferLogs = False
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ConfFile_cfg.py' #put here the name of analyzer config file
config.JobType.allowUndistributedCMSSW = True
config.Data.inputDataset = '/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer19UL17MiniAOD-106X_mc2017_realistic_v6-v2/MINIAODSIM' #dataset as it apperas in DAS query
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 20
config.Data.lumiMask = 'lumiMask_TTToSemileptonic.txt'
#config.Data.runRange = '297050-297100'
config.Data.outLFNDirBase = '/store/user/fiemmi'
config.Data.publication = False
config.Site.storageSite = 'T2_IT_Legnaro'

