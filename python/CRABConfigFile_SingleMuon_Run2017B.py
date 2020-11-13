from CRABClient.UserUtilities import config

config = config()
config.General.requestName = 'SingleMuon_Run2017B-09Aug2019_UL2017-v1'
config.General.workArea = 'crab_Sep05_2020_h1830_SingleMuon'
config.General.transferOutputs = True
config.General.transferLogs = False
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ConfFile_cfg.py' #put here the name of analyzer config file
config.JobType.allowUndistributedCMSSW = True
config.Data.inputDataset = '/SingleMuon/Run2017B-09Aug2019_UL2017-v1/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 20
config.Data.lumiMask = 'goldenJSON_Run2017B.txt'
config.Data.runRange = '297050-297100'
config.Data.outLFNDirBase = '/store/user/fiemmi'
config.Data.publication = False
config.Site.storageSite = 'T2_IT_Legnaro'

