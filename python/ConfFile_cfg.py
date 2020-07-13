import FWCore.ParameterSet.Config as cms

process = cms.Process("PFCandInfo")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(21000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    #fileNames = cms.untracked.vstring('file:root://xrootd-cms.infn.it///store/mc/RunIISummer19UL17MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/FlatPU0to70_106X_mc2017_realistic_v6-v3/40000/FFEE2656-43CF-4247-8A06-8B04E8FF00F5.root'),
    #fileNames = cms.untracked.vstring('file:root://xrootd-cms.infn.it///store/mc/RunIISummer19UL17MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/NoPU_106X_mc2017_realistic_v6-v3/40000/FDBFC4B3-F359-5B49-8F79-2075C6024381.root'),
    fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_1.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_2.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_3.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_4.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_5.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_6.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_7.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_8.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_9.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_10.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_11.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_12.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_13.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_14.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_15.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_16.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_17.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_18.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_19.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_20.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_21.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_22.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_23.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_24.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_25.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_26.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_27.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_28.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_29.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_30.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_31.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_32.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_33.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_34.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_35.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_36.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_37.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_38.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_39.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_40.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_41.root', 'file:/afs/cern.ch/user/f/fiemmi/JetMET/CMSSW_10_6_4/src/pickevents_PU_files_EXT_bis/pickevents_42.root',),
    

    
)

process.GetPFInfo = cms.EDAnalyzer('PFCandInfoAnalyzer',
                                   vertices  = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                   PUinfo  = cms.InputTag("slimmedAddPileupInfo"),
                                   AK4PUPPIJets  = cms.InputTag("slimmedJetsPuppi"),
                                   AK4CHSJets  = cms.InputTag("slimmedJets"),
                                   AK4GenJets = cms.InputTag("slimmedGenJets"),
                                   PFCands  = cms.InputTag("packedPFCandidates"),
                                   
)

process.TFileService = cms.Service("TFileService",
                                       #fileName = cms.string('flatTree_QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8_ext.root')
                                       #fileName = cms.string('flatTree_QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8_PU_EXT.root')
                                       #fileName = cms.string('flatTree_QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8_training_noPU_EXT.root')
                                       fileName = cms.string('flatTree_QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8_training_PU_EXT.root')
                                       #fileName = cms.string('try.root')
                                   )

process.p = cms.Path(process.GetPFInfo)
