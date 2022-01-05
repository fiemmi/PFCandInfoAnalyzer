import FWCore.ParameterSet.Config as cms

########### this is to get input files from command line ########## 
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing ('analysis')
options.parseArguments()

#################################################################

process = cms.Process("PFCandInfo")
#

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(80000) )

isMC = True
isminiAODv1 = True

#load globaltag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
if (isMC) : 
    process.GlobalTag.globaltag = "106X_mc2017_realistic_v8"
    print("Using globaltag {0}".format(process.GlobalTag.globaltag))
 
else :
    process.GlobalTag.globaltag = "106X_dataRun2_v32"
    print("Using globaltag {0}".format(process.GlobalTag.globaltag))

process.source = cms.Source("PoolSource",
    #FLAT PU SAMPLES (mind: they are v1 miniAODs and Summer19 MC campaign!)
    fileNames = cms.untracked.vstring('file:root://xrootd-cms.infn.it///store/mc/RunIISummer19UL17MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/FlatPU0to70_106X_mc2017_realistic_v6-v3/40000/FFEE2656-43CF-4247-8A06-8B04E8FF00F5.root'),
    #fileNames = cms.untracked.vstring('file:root://xrootd-cms.infn.it///store/mc/RunIISummer19UL17MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/EpsilonPU_106X_mc2017_realistic_v6-v2/100000/107F8FFF-6F20-A14D-BC89-E0C761BC300D.root'),
    #NEW miniAODv2 FILES WITH PUPPIv15 WEIGHTS
    #EpsiloPU
    #fileNames = cms.untracked.vstring('file:root://cmsxrootd.fnal.gov///store/mc/RunIISummer19UL17MiniAODv2/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/EpsilonPU_106X_mc2017_realistic_v9-v1/00000/1A7AF728-10B7-2345-A49A-FBAB937CEADF.root', 'file:root://cmsxrootd.fnal.gov///store/mc/RunIISummer19UL17MiniAODv2/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/EpsilonPU_106X_mc2017_realistic_v9-v1/00000/275CC5DF-0EA5-FD40-B3FB-BF2A6309A24E.root'),
    #with PU
    #fileNames = cms.untracked.vstring('file:root://cmsxrootd.fnal.gov///store/mc/RunIISummer19UL17MiniAODv2/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/106X_mc2017_realistic_v9-v1/00000/00379E19-FB20-7044-83F6-7A9E41F1BD1D.root'),
    #FILES FROM emdPickEvents
    #fileNames = cms.untracked.vstring('file:/eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/pickevents_PU_files_EXT40k/pickevents_1.root', 'file:/eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/pickevents_PU_files_EXT40k/pickevents_2.root', 'file:/eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/pickevents_PU_files_EXT40k/pickevents_3.root', 'file:/eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/pickevents_PU_files_EXT40k/pickevents_4.root', 'file:/eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/pickevents_PU_files_EXT40k/pickevents_5.root', 'file:/eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/pickevents_PU_files_EXT40k/pickevents_6.root', 'file:/eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/pickevents_PU_files_EXT40k/pickevents_7.root', 'file:/eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/pickevents_PU_files_EXT40k/pickevents_8.root',),
    #fileNames = cms.untracked.vstring('file:/eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/pickevents_PU_files_EXT80k/pickevents_1.root', 'file:/eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/pickevents_PU_files_EXT80k/pickevents_2.root', 'file:/eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/pickevents_PU_files_EXT80k/pickevents_3.root', 'file:/eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/pickevents_PU_files_EXT80k/pickevents_4.root', 'file:/eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/pickevents_PU_files_EXT80k/pickevents_5.root', 'file:/eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/pickevents_PU_files_EXT80k/pickevents_6.root', 'file:/eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/pickevents_PU_files_EXT80k/pickevents_7.root', 'file:/eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/pickevents_PU_files_EXT80k/pickevents_8.root', 'file:/eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/pickevents_PU_files_EXT80k/pickevents_9.root', 'file:/eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/pickevents_PU_files_EXT80k/pickevents_10.root', 'file:/eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/pickevents_PU_files_EXT80k/pickevents_11.root',  'file:/eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/pickevents_PU_files_EXT80k/pickevents_12.root',  'file:/eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/pickevents_PU_files_EXT80k/pickevents_13.root', 'file:/eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/pickevents_PU_files_EXT80k/pickevents_14.root', 'file:/eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/pickevents_PU_files_EXT80k/pickevents_15.root', 'file:/eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/pickevents_PU_files_EXT80k/pickevents_16.root', 'file:/eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/pickevents_PU_files_EXT80k/pickevents_17.root',),
    #TO BE USED WITH CONDOR
    #fileNames = cms.untracked.vstring(options.inputFiles),
    #----- SEMILEPTONIC TTBAR -------
    #SUMMER19
    #fileNames = cms.untracked.vstring('file:root://xrootd-cms.infn.it///store/mc/RunIISummer19UL17MiniAODv2/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_mc2017_realistic_v9-v1/020000/053108A0-62AD-E742-AC68-A88BABB476D6.root'),
    #SUMMER20
    #fileNames = cms.untracked.vstring('file:root://xrootd-cms.infn.it///store/mc/RunIISummer20UL17MiniAODv2/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_mc2017_realistic_v9-v1/00000/005708B7-331C-904E-88B9-189011E6C9DD.root'),
    #-------- DATA --------
    #JetHT 
    #fileNames = cms.untracked.vstring('file:root://xrootd-cms.infn.it///store/data/Run2017B/JetHT/MINIAOD/09Aug2019_UL2017-v1/270000/FC1A877A-9874-D143-B7D8-E16F1F1E2BB1.root'),
    #SingleMuon
    #fileNames = cms.untracked.vstring('file:root://xrootd-cms.infn.it///store/data/Run2017B/SingleMuon/MINIAOD/09Aug2019_UL2017-v1/270000/528F724C-F018-9B4E-B1DF-531A1DEDAAFA.root'),
    #SingleElectron
    #fileNames = cms.untracked.vstring('file:root://xrootd-cms.infn.it///store/data/Run2017B/SingleElectron/MINIAOD/09Aug2019_UL2017-v1/130000/005B3780-5B9B-3B4F-ABA9-DCC474BE15BE.root'),

)

print("options.inputFiles", options.inputFiles)

process.triggerSelectionQCD = cms.EDFilter("TriggerResultsFilter",
                                       triggerConditions = cms.vstring('HLT_PFHT1050_v*'),
                                       hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                       l1tResults = cms.InputTag( "" ),
                                       throw = cms.bool(False)
                                       )

process.triggerSelectionTTToSemileptonic = cms.EDFilter("TriggerResultsFilter",
                                       triggerConditions = cms.vstring(
                                       "HLT_IsoMu27_v*",
                                       "HLT_Ele35_WPTight_Gsf_v*",
                                       ),
                                       hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                       l1tResults = cms.InputTag( "" ),
                                       throw = cms.bool(False)
                                       )

if isminiAODv1:
    """
    Update PUPPI to v15
    """

    from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
    process.ak4PuppiJets  = ak4PFJets.clone (src = 'puppi', doAreaFastjet = True, jetPtMin = 2.)

    from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
    addJetCollection(process, labelName = 'Puppi', jetSource = cms.InputTag('ak4PuppiJets'), algo = 'AK', rParam=0.4, genJetCollection=cms.InputTag('slimmedGenJets'), jetCorrections = ('AK4PFPuppi', ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual'], 'None'), pfCandidates = cms.InputTag('packedPFCandidates'), pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'), svSource = cms.InputTag('slimmedSecondaryVertices'), muSource =cms.InputTag( 'slimmedMuons'), elSource = cms.InputTag('slimmedElectrons'), genParticles= cms.InputTag('prunedGenParticles'), getJetMCFlavour=isMC)                                                                                                        

    process.patJetsPuppi.addGenPartonMatch = cms.bool(isMC)                                                                                                          
    process.patJetsPuppi.addGenJetMatch = cms.bool(isMC)

    from CommonTools.PileupAlgos.customizePuppiTune_cff import UpdatePuppiTuneV15
    UpdatePuppiTuneV15(process, isMC)

process.GetPFInfo = cms.EDAnalyzer('PFCandInfoAnalyzer',
                                   vertices  = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                   PUinfo  = cms.InputTag("slimmedAddPileupInfo"),
                                   AK4PUPPIJets  = cms.InputTag("slimmedJetsPuppi"),
                                   AK4CHSJets  = cms.InputTag("slimmedJets"),
                                   AK4GenJets = cms.InputTag("slimmedGenJets"),
                                   AK8PUPPIJets = cms.InputTag('slimmedJetsAK8'),
                                   AK8PUPPIJetSDMass = cms.string('ak8PFJetsPuppiSoftDropMass'),
                                   AK8GenJets = cms.InputTag('slimmedGenJetsAK8'),
                                   genParticles = cms.InputTag('packedGenParticles'),
                                   btaggerCSVv2 = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
                                   CSVv2WP = cms.double(1.0), #not used for now
                                   btaggerDeepCSV = cms.string("pfDeepCSVJetTags:probb"),
                                   DeepCSVWP = cms.double(0.4506), #WP for UL2017
                                   PFCands  = cms.InputTag("packedPFCandidates"),
                                   muons = cms.InputTag("slimmedMuons"),
                                   electrons = cms.InputTag("slimmedElectrons"),
                                   missingEt = cms.InputTag("slimmedMETs"),
                                   PUPPImissingEt = cms.InputTag("slimmedMETsPuppi"),
                                   triggerNames = cms.vstring (
                                       "HLT_PFHT1050_v",
                                       "HLT_PFJet500_v",
                                       "HLT_PFJet550_v",
                                       "HLT_IsoMu27_v",
                                       "HLT_Ele35_WPTight_Gsf_v",
                                       
                                   ),
                                   triggerResults = cms.InputTag('TriggerResults','','HLT'),
                                   runOnMC = cms.untracked.bool(isMC)
)

process.TFileService = cms.Service("TFileService",
                                       #fileName = cms.string('flatTree_QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8_ext.root')
                                       #fileName = cms.string('flatTree_QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8_PU_EXT.root')
                                       #fileName = cms.string('flatTree_QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8_training_EpsilonPU_EXT350k_withPUPPIalpha.root')
                                       fileName = cms.string('flatTree_QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8_training_EpsilonPU_EXT80k_withPUPPIalpha_miniAODv1_withAK8.root')
                                       #fileName = cms.string('flatTree_QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8_training_PU_EXT80k_.root')
                                       #fileName = cms.string('flatTree_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_fixCHSbug.root')
                                       #fileName = cms.string('flatTree_JetHT_Run2017B.root')
                                       #fileName = cms.string('flatTree_SingleMuon_Run2017B.root')
                                       #fileName = cms.string('flatTree_SingleElectron_Run2017B.root')
                                       #fileName = cms.string('file.root')
                                       #TO BE USED WITH CONDOR
                                       #fileName = cms.string(options.outputFile),
                                   )
if (not isMC) :
    if isminiAODv1:
        process.p = cms.Path(process.triggerSelectionQCD*process.puppiSequence*process.GetPFInfo)
    else:
        process.p = cms.Path(process.triggerSelectionQCD*process.GetPFInfo)
else:
    if isminiAODv1:
        process.p = cms.Path(process.puppiSequence*process.GetPFInfo)
    else:
        process.p = cms.Path(process.GetPFInfo)

#print process.dumpPython() 
