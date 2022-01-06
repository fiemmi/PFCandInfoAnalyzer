#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TBenchmark.h>
#include <iostream>

using namespace std;

void createTreePU_noPU (int nfile) {

gBenchmark->Start("running time");    
    
 TFile *inputfile_noPU  = new TFile(Form("/afs/cern.ch/work/f/fiemmi/private/CMSSW_10_6_20/src/PFCandInfo/PFCandInfoAnalyzer/sliced_noPUfiles/files/EXT80k_miniAODv1_fixCHSbug/flatTree_QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8_training_EpsilonPU_EXT80k_withPUPPIalpha_miniAODv1_fixCHSbug_part%i.root", nfile), "READ" );

 TFile *inputfile_PU  = new TFile(Form("/afs/cern.ch/work/f/fiemmi/private/CMSSW_10_6_20/src/PFCandInfo/PFCandInfoAnalyzer/sorting_framework/sorted_files/EXT80k_miniAODv1_fixCHSbug/flatTree_QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8_training_PU_EXT80k_withPUPPIalpha_miniAODv1_fixCHSbug_%i.root", nfile), "READ" );

 int runNo, evtNo, lumiSec, nPUint_PU, nPFCands_PU, nPFCands_noPU, nAK4PUPPIJets_PU, nAK4PUPPIJets_noPU, nAK8PUPPIJets_PU, nAK8PUPPIJets_noPU, nAK4CHSJets_PU, nAK4CHSJets_noPU, nAK4GenJets_PU /*(same as noPU)*/, nAK8GenJets_PU /*(same as noPU)*/, nGenParticles_noPU /*(same as PU)*/, nLeptons_PU;
 float CHSMET_PU, CHSUnclMET_PU, RawCHSMET_PU, RawCHSUnclMET_PU, PUPPIMET_PU, PUPPIUnclMET_PU, RawPUPPIMET_PU, RawPUPPIUnclMET_PU, CHSMET_noPU, CHSUnclMET_noPU, RawCHSMET_noPU, RawCHSUnclMET_noPU, PUPPIMET_noPU, PUPPIUnclMET_noPU, RawPUPPIMET_noPU, RawPUPPIUnclMET_noPU, genMET, genUnclMET, VBFDijetGenMass, VBFDijetGenMass_PU, VBFDijetGenMass_noPU, VBFDijetCHSMass_PU, VBFDijetCHSMass_noPU, VBFDijetPUPPIMass_PU, VBFDijetPUPPIMass_noPU;
 //PF candidates
 vector <float> PFCandPt_PU, PFCandPt_noPU, PFCandEta_PU, PFCandEta_noPU, PFCandAbsEta_PU, PFCandPhi_PU, PFCandPhi_noPU, PFCandE_PU, PFCandE_noPU, PFCandpdgId_PU, PFCandpdgId_noPU, PFCandCharge_PU, PFCandCharge_noPU, PFCandPUPPIw_PU, PFCandPUPPIw_noPU, PFCandPUPPIalpha_PU, PFCandPUPPIalpha_noPU, PFCandHCalFrac_PU, PFCandHCalFrac_noPU, PFCandHCalFracCalib_PU, PFCandHCalFracCalib_noPU, PFCandVtxAssQual_PU, PFCandVtxAssQual_noPU, PFCandFromPV_PU, PFCandFromPV_noPU, PFCandLostInnerHits_PU, PFCandLostInnerHits_noPU, PFCandTrackHighPurity_PU, PFCandTrackHighPurity_noPU, PFCandDZ_PU, PFCandDZ_noPU, PFCandDXY_PU, PFCandDXY_noPU, PFCandDZsig_PU, PFCandDZsig_noPU, PFCandDXYsig_PU, PFCandDXYsig_noPU, PFCandNormChi2_PU, PFCandNormChi2_noPU, PFCandQuality_PU, PFCandQuality_noPU, PFCandNumHits_PU, PFCandNumHits_noPU, PFCandNumLayersHit_PU, PFCandNumLayersHit_noPU;
 //AK4 jets
 vector <float> AK4PUPPIJetPt_PU, AK4PUPPIJetRawPt_PU, AK4PUPPIJetEta_PU, AK4PUPPIJetPhi_PU, AK4PUPPIJetE_PU, AK4PUPPIJetRawE_PU, AK4PUPPIJetPt_noPU, AK4PUPPIJetRawPt_noPU, AK4PUPPIJetEta_noPU, AK4PUPPIJetPhi_noPU, AK4PUPPIJetE_noPU, AK4PUPPIJetRawE_noPU, AK4CHSJetPt_PU, AK4CHSJetRawPt_PU, AK4CHSJetEta_PU, AK4CHSJetPhi_PU, AK4CHSJetE_PU, AK4CHSJetRawE_PU, AK4CHSJetPt_noPU, AK4CHSJetRawPt_noPU, AK4CHSJetEta_noPU, AK4CHSJetPhi_noPU, AK4CHSJetE_noPU, AK4CHSJetRawE_noPU, AK4GenJetPt_PU, AK4GenJetEta_PU, AK4GenJetPhi_PU, AK4GenJetE_PU;
 //AK8 jets
 vector <float> AK8PUPPIJetPt_PU, AK8PUPPIJetRawPt_PU, AK8PUPPIJetEta_PU, AK8PUPPIJetPhi_PU, AK8PUPPIJetE_PU, AK8PUPPIJetRawE_PU, AK8PUPPIJetPt_noPU, AK8PUPPIJetRawPt_noPU, AK8PUPPIJetEta_noPU, AK8PUPPIJetPhi_noPU, AK8PUPPIJetE_noPU, AK8PUPPIJetRawE_noPU, AK8GenJetPt_PU, AK8GenJetEta_PU, AK8GenJetPhi_PU, AK8GenJetE_PU;
 //Gen particles
 vector <float> genParticlePt_noPU, genParticleEta_noPU, genParticlePhi_noPU, genParticleE_noPU, genParticleCharge_noPU, genParticlepdgId_noPU;
 vector <bool> triggerBit_PU, AK4CHSJetIsBtag_PU;
 
 TFile * outputfile  = new TFile( Form("flatTree_QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8_training_PU+EpsilonPU_EXT80k_withPUPPIalpha_miniAODv1_fixCHSbug_%i.root", nfile), "RECREATE" );
 TTree * flatTree = new TTree( "events", "events" );
 flatTree->Branch("runNo", &runNo, "runNo/I");
 flatTree->Branch("evtNo", &evtNo, "evtNo/I");
 flatTree->Branch("lumiSec", &lumiSec, "lumiSec/I");
 flatTree->Branch("nPUint", &nPUint_PU, "nPUint_PU/I");
 flatTree->Branch("triggerBit_PU", &triggerBit_PU);
 flatTree->Branch("nLeptons_PU", &nLeptons_PU, "nLeptons_PU/I");
 flatTree->Branch("nPFCands_noPU", &nPFCands_noPU, "nPFCands_noPU/I");
 flatTree->Branch("nPFCands_PU", &nPFCands_PU, "nPFCands_PU/I");
 flatTree->Branch("PFCandPt_noPU", &PFCandPt_noPU);
 flatTree->Branch("PFCandPt_PU", &PFCandPt_PU);
 flatTree->Branch("PFCandEta_noPU", &PFCandEta_noPU);
 flatTree->Branch("PFCandEta_PU", &PFCandEta_PU);
 flatTree->Branch("PFCandEta_noPU", &PFCandEta_noPU);
 flatTree->Branch("PFCandEta_PU", &PFCandEta_PU);
 flatTree->Branch("PFCandAbsEta_PU", &PFCandAbsEta_PU);
 flatTree->Branch("PFCandPhi_noPU", &PFCandPhi_noPU);
 flatTree->Branch("PFCandPhi_PU", &PFCandPhi_PU);
 flatTree->Branch("PFCandE_noPU", &PFCandE_noPU);
 flatTree->Branch("PFCandE_PU", &PFCandE_PU);
 flatTree->Branch("PFCandpdgId_PU", &PFCandpdgId_PU);
 flatTree->Branch("PFCandpdgId_noPU", &PFCandpdgId_noPU);
 flatTree->Branch("PFCandCharge_PU", &PFCandCharge_PU);
 flatTree->Branch("PFCandCharge_noPU", &PFCandCharge_noPU);
 flatTree->Branch("PFCandPUPPIw_PU", &PFCandPUPPIw_PU);
 flatTree->Branch("PFCandPUPPIw_noPU", &PFCandPUPPIw_noPU);
 flatTree->Branch("PFCandPUPPIalpha_PU", &PFCandPUPPIalpha_PU);
 flatTree->Branch("PFCandPUPPIalpha_noPU", &PFCandPUPPIalpha_noPU);
 flatTree->Branch("PFCandHCalFrac_PU", &PFCandHCalFrac_PU);
 flatTree->Branch("PFCandHCalFrac_noPU", &PFCandHCalFrac_noPU);
 flatTree->Branch("PFCandHCalFracCalib_PU", &PFCandHCalFracCalib_PU);
 flatTree->Branch("PFCandHCalFracCalib_noPU", &PFCandHCalFracCalib_noPU);
 flatTree->Branch("PFCandVtxAssQual_PU", &PFCandVtxAssQual_PU);
 flatTree->Branch("PFCandVtxAssQual_noPU", &PFCandVtxAssQual_noPU);
 flatTree->Branch("PFCandFromPV_PU", &PFCandFromPV_PU);
 flatTree->Branch("PFCandFromPV_noPU", &PFCandFromPV_noPU);
 flatTree->Branch("PFCandLostInnerHits_PU", &PFCandLostInnerHits_PU);
 flatTree->Branch("PFCandLostInnerHits_noPU", &PFCandLostInnerHits_noPU);
 flatTree->Branch("PFCandTrackHighPurity_PU", &PFCandTrackHighPurity_PU);
 flatTree->Branch("PFCandTrackHighPurity_noPU", &PFCandTrackHighPurity_noPU);
 flatTree->Branch("PFCandDZ_PU", &PFCandDZ_PU);
 flatTree->Branch("PFCandDZ_noPU", &PFCandDZ_noPU);
 flatTree->Branch("PFCandDXY_PU", &PFCandDXY_PU);
 flatTree->Branch("PFCandDXY_noPU", &PFCandDXY_noPU);
 flatTree->Branch("PFCandDZsig_PU", &PFCandDZsig_PU);
 flatTree->Branch("PFCandDZsig_noPU", &PFCandDZsig_noPU);
 flatTree->Branch("PFCandDXYsig_PU", &PFCandDXYsig_PU);
 flatTree->Branch("PFCandDXYsig_noPU", &PFCandDXYsig_noPU);
 flatTree->Branch("PFCandNormChi2_PU", &PFCandNormChi2_PU);
 flatTree->Branch("PFCandNormChi2_noPU", &PFCandNormChi2_noPU);
 flatTree->Branch("PFCandQuality_PU", &PFCandQuality_PU);
 flatTree->Branch("PFCandQuality_noPU", &PFCandQuality_noPU);
 flatTree->Branch("PFCandNumHits_PU", &PFCandNumHits_PU);
 flatTree->Branch("PFCandNumHits_noPU", &PFCandNumHits_noPU);
 flatTree->Branch("PFCandNumLayersHit_PU", &PFCandNumLayersHit_PU);
 flatTree->Branch("PFCandNumLayersHit_noPU", &PFCandNumLayersHit_noPU);
 flatTree->Branch("nAK4PUPPIJets_PU", &nAK4PUPPIJets_PU);
 flatTree->Branch("AK4PUPPIJetPt_PU", &AK4PUPPIJetPt_PU);
 flatTree->Branch("AK4PUPPIJetRawPt_PU", &AK4PUPPIJetRawPt_PU);
 flatTree->Branch("AK4PUPPIJetEta_PU", &AK4PUPPIJetEta_PU);
 flatTree->Branch("AK4PUPPIJetPhi_PU", &AK4PUPPIJetPhi_PU);
 flatTree->Branch("AK4PUPPIJetE_PU", &AK4PUPPIJetE_PU);
 flatTree->Branch("AK4PUPPIJetRawE_PU", &AK4PUPPIJetRawE_PU);
 flatTree->Branch("nAK4PUPPIJets_noPU", &nAK4PUPPIJets_noPU);
 flatTree->Branch("AK4PUPPIJetPt_noPU", &AK4PUPPIJetPt_noPU);
 flatTree->Branch("AK4PUPPIJetRawPt_noPU", &AK4PUPPIJetRawPt_noPU);
 flatTree->Branch("AK4PUPPIJetEta_noPU", &AK4PUPPIJetEta_noPU);
 flatTree->Branch("AK4PUPPIJetPhi_noPU", &AK4PUPPIJetPhi_noPU);
 flatTree->Branch("AK4PUPPIJetE_noPU", &AK4PUPPIJetE_noPU);
 flatTree->Branch("AK4PUPPIJetRawE_noPU", &AK4PUPPIJetRawE_noPU);
 flatTree->Branch("nAK8PUPPIJets_PU", &nAK8PUPPIJets_PU);
 flatTree->Branch("AK8PUPPIJetPt_PU", &AK8PUPPIJetPt_PU);
 flatTree->Branch("AK8PUPPIJetRawPt_PU", &AK8PUPPIJetRawPt_PU);
 flatTree->Branch("AK8PUPPIJetEta_PU", &AK8PUPPIJetEta_PU);
 flatTree->Branch("AK8PUPPIJetPhi_PU", &AK8PUPPIJetPhi_PU);
 flatTree->Branch("AK8PUPPIJetE_PU", &AK8PUPPIJetE_PU);
 flatTree->Branch("AK8PUPPIJetRawE_PU", &AK8PUPPIJetRawE_PU);
 flatTree->Branch("nAK8PUPPIJets_noPU", &nAK8PUPPIJets_noPU);
 flatTree->Branch("AK8PUPPIJetPt_noPU", &AK8PUPPIJetPt_noPU);
 flatTree->Branch("AK8PUPPIJetRawPt_noPU", &AK8PUPPIJetRawPt_noPU);
 flatTree->Branch("AK8PUPPIJetEta_noPU", &AK8PUPPIJetEta_noPU);
 flatTree->Branch("AK8PUPPIJetPhi_noPU", &AK8PUPPIJetPhi_noPU);
 flatTree->Branch("AK8PUPPIJetE_noPU", &AK8PUPPIJetE_noPU);
 flatTree->Branch("AK8PUPPIJetRawE_noPU", &AK8PUPPIJetRawE_noPU);
 flatTree->Branch("nAK4CHSJets_PU", &nAK4CHSJets_PU);
 flatTree->Branch("AK4CHSJetPt_PU", &AK4CHSJetPt_PU);
 flatTree->Branch("AK4CHSJetRawPt_PU", &AK4CHSJetRawPt_PU);
 flatTree->Branch("AK4CHSJetEta_PU", &AK4CHSJetEta_PU);
 flatTree->Branch("AK4CHSJetPhi_PU", &AK4CHSJetPhi_PU);
 flatTree->Branch("AK4CHSJetE_PU", &AK4CHSJetE_PU);
 flatTree->Branch("AK4CHSJetRawE_PU", &AK4CHSJetRawE_PU);
 flatTree->Branch("AK4CHSJetIsBtag_PU", &AK4CHSJetIsBtag_PU);
 flatTree->Branch("nAK4CHSJets_noPU", &nAK4CHSJets_noPU);
 flatTree->Branch("AK4CHSJetPt_noPU", &AK4CHSJetPt_noPU);
 flatTree->Branch("AK4CHSJetRawPt_noPU", &AK4CHSJetRawPt_noPU);
 flatTree->Branch("AK4CHSJetEta_noPU", &AK4CHSJetEta_noPU);
 flatTree->Branch("AK4CHSJetPhi_noPU", &AK4CHSJetPhi_noPU);
 flatTree->Branch("AK4CHSJetE_noPU", &AK4CHSJetE_noPU);
 flatTree->Branch("AK4CHSJetRawE_noPU", &AK4CHSJetRawE_noPU); 
 flatTree->Branch("nAK4GenJets_PU", &nAK4GenJets_PU);
 flatTree->Branch("AK4GenJetPt_PU", &AK4GenJetPt_PU);
 flatTree->Branch("AK4GenJetEta_PU", &AK4GenJetEta_PU);
 flatTree->Branch("AK4GenJetPhi_PU", &AK4GenJetPhi_PU);
 flatTree->Branch("AK4GenJetE_PU", &AK4GenJetE_PU);
 flatTree->Branch("nAK8GenJets_PU", &nAK8GenJets_PU);
 flatTree->Branch("AK8GenJetPt_PU", &AK8GenJetPt_PU);
 flatTree->Branch("AK8GenJetEta_PU", &AK8GenJetEta_PU);
 flatTree->Branch("AK8GenJetPhi_PU", &AK8GenJetPhi_PU);
 flatTree->Branch("AK8GenJetE_PU", &AK8GenJetE_PU);
 flatTree->Branch("nGenParticles_noPU", &nGenParticles_noPU, "nGenParticles_noPU/I");
 flatTree->Branch("genParticlePt_noPU", &genParticlePt_noPU);
 flatTree->Branch("genParticleEta_noPU", &genParticleEta_noPU);
 flatTree->Branch("genParticlePhi_noPU", &genParticlePhi_noPU);
 flatTree->Branch("genParticleE_noPU", &genParticleE_noPU);
 flatTree->Branch("genParticleCharge_noPU", &genParticleCharge_noPU);
 flatTree->Branch("genParticlepdgId_noPU", &genParticlepdgId_noPU);
 flatTree->Branch("CHSMET_PU", &CHSMET_PU);
 flatTree->Branch("CHSUnclMET_PU", &CHSUnclMET_PU);
 flatTree->Branch("RawCHSMET_PU", &RawCHSMET_PU);
 flatTree->Branch("RawCHSUnclMET_PU", &RawCHSUnclMET_PU);
 flatTree->Branch("PUPPIMET_PU", &PUPPIMET_PU);
 flatTree->Branch("PUPPIUnclMET_PU", &PUPPIUnclMET_PU);
 flatTree->Branch("RawPUPPIMET_PU", &RawPUPPIMET_PU);
 flatTree->Branch("RawPUPPIUnclMET_PU", &RawPUPPIUnclMET_PU);
 flatTree->Branch("CHSMET_noPU", &CHSMET_noPU);
 flatTree->Branch("CHSUnclMET_noPU", &CHSUnclMET_noPU);
 flatTree->Branch("RawCHSMET_noPU", &RawCHSMET_noPU);
 flatTree->Branch("RawCHSUnclMET_noPU", &RawCHSUnclMET_noPU);
 flatTree->Branch("PUPPIMET_noPU", &PUPPIMET_noPU);
 flatTree->Branch("PUPPIUnclMET_noPU", &PUPPIUnclMET_noPU);
 flatTree->Branch("RawPUPPIMET_noPU", &RawPUPPIMET_noPU);
 flatTree->Branch("RawPUPPIUnclMET_noPU", &RawPUPPIUnclMET_noPU);
 flatTree->Branch("genMET", &genMET);
 flatTree->Branch("genUnclMET", &genUnclMET);
 flatTree->Branch("VBFDijetGenMass", &VBFDijetGenMass);
 flatTree->Branch("VBFDijetGenMass_noPU", &VBFDijetGenMass_noPU);
 flatTree->Branch("VBFDijetGenMass_PU", &VBFDijetGenMass_PU);
 flatTree->Branch("VBFDijetPUPPIMass_noPU", &VBFDijetPUPPIMass_noPU);
 flatTree->Branch("VBFDijetCHSMass_noPU", &VBFDijetCHSMass_noPU);
 flatTree->Branch("VBFDijetPUPPIMass_PU", &VBFDijetPUPPIMass_PU);
 flatTree->Branch("VBFDijetCHSMass_PU", &VBFDijetCHSMass_PU);


 TTree *evt_noPU = (TTree*)inputfile_noPU->Get( "events" );

 UInt_t myrunNo_noPU = 0;
 evt_noPU->SetBranchAddress( "runNo", &myrunNo_noPU );

 UInt_t myevtNo_noPU = 0;
 evt_noPU->SetBranchAddress( "evtNo", &myevtNo_noPU );

 UInt_t mylumiSec_noPU = 0;
 evt_noPU->SetBranchAddress( "lumiSec", &mylumiSec_noPU );

 UInt_t mynPFCands_noPU = 0;
 evt_noPU->SetBranchAddress( "nPFCands", &mynPFCands_noPU );

 vector<float> *myPFCandPt_noPU = 0;
 evt_noPU->SetBranchAddress( "PFCandPt", &myPFCandPt_noPU );

 vector<float> *myPFCandEta_noPU = 0;
 evt_noPU->SetBranchAddress( "PFCandEta", &myPFCandEta_noPU );

 vector<float> *myPFCandPhi_noPU = 0;
 evt_noPU->SetBranchAddress( "PFCandPhi", &myPFCandPhi_noPU );

 vector<float> *myPFCandE_noPU = 0;
 evt_noPU->SetBranchAddress( "PFCandE", &myPFCandE_noPU );

 vector<float> *myPFCandpdgId_noPU = 0;
 evt_noPU->SetBranchAddress( "PFCandpdgId", &myPFCandpdgId_noPU );

 vector<float> *myPFCandCharge_noPU = 0;
 evt_noPU->SetBranchAddress( "PFCandCharge", &myPFCandCharge_noPU );

 vector<float> *myPFCandPUPPIw_noPU = 0;
 evt_noPU->SetBranchAddress( "PFCandPUPPIw", &myPFCandPUPPIw_noPU );

 vector<float> *myPFCandPUPPIalpha_noPU = 0;
 evt_noPU->SetBranchAddress( "PFCandPUPPIalpha", &myPFCandPUPPIalpha_noPU );

 vector<float> *myPFCandHCalFrac_noPU = 0;
 evt_noPU->SetBranchAddress( "PFCandHCalFrac", &myPFCandHCalFrac_noPU );

 vector<float> *myPFCandHCalFracCalib_noPU = 0;
 evt_noPU->SetBranchAddress( "PFCandHCalFracCalib", &myPFCandHCalFracCalib_noPU );

 vector<float> *myPFCandVtxAssQual_noPU = 0;
 evt_noPU->SetBranchAddress( "PFCandVtxAssQual", &myPFCandVtxAssQual_noPU );

 vector<float> *myPFCandFromPV_noPU = 0;
 evt_noPU->SetBranchAddress( "PFCandFromPV", &myPFCandFromPV_noPU );

 vector<float> *myPFCandLostInnerHits_noPU = 0;
 evt_noPU->SetBranchAddress( "PFCandLostInnerHits", &myPFCandLostInnerHits_noPU );

 vector<float> *myPFCandTrackHighPurity_noPU = 0;
 evt_noPU->SetBranchAddress( "PFCandTrackHighPurity", &myPFCandTrackHighPurity_noPU );

 vector<float> *myPFCandDZ_noPU = 0;
 evt_noPU->SetBranchAddress( "PFCandDZ", &myPFCandDZ_noPU );

 vector<float> *myPFCandDXY_noPU = 0;
 evt_noPU->SetBranchAddress( "PFCandDXY", &myPFCandDXY_noPU );

 vector<float> *myPFCandDZsig_noPU = 0;
 evt_noPU->SetBranchAddress( "PFCandDZsig", &myPFCandDZsig_noPU );

 vector<float> *myPFCandDXYsig_noPU = 0;
 evt_noPU->SetBranchAddress( "PFCandDXYsig", &myPFCandDXYsig_noPU );

 vector<float> *myPFCandNormChi2_noPU = 0;
 evt_noPU->SetBranchAddress( "PFCandNormChi2", &myPFCandNormChi2_noPU );

 vector<float> *myPFCandQuality_noPU = 0;
 evt_noPU->SetBranchAddress( "PFCandQuality", &myPFCandQuality_noPU );

 vector<float> *myPFCandNumHits_noPU = 0;
 evt_noPU->SetBranchAddress( "PFCandNumHits", &myPFCandNumHits_noPU );

 vector<float> *myPFCandNumLayersHit_noPU = 0;
 evt_noPU->SetBranchAddress( "PFCandNumLayersHit", &myPFCandNumLayersHit_noPU );

 UInt_t mynAK4PUPPIJets_noPU = 0;
 evt_noPU->SetBranchAddress( "nAK4PUPPIJets", &mynAK4PUPPIJets_noPU );

 vector<float> *myAK4PUPPIJetPt_noPU = 0;
 evt_noPU->SetBranchAddress( "AK4PUPPIJetPt", &myAK4PUPPIJetPt_noPU );

 vector<float> *myAK4PUPPIJetRawPt_noPU = 0;
 evt_noPU->SetBranchAddress( "AK4PUPPIJetRawPt", &myAK4PUPPIJetRawPt_noPU );

 vector<float> *myAK4PUPPIJetEta_noPU = 0;
 evt_noPU->SetBranchAddress( "AK4PUPPIJetEta", &myAK4PUPPIJetEta_noPU );

 vector<float> *myAK4PUPPIJetPhi_noPU = 0;
 evt_noPU->SetBranchAddress( "AK4PUPPIJetPhi", &myAK4PUPPIJetPhi_noPU );

 vector<float> *myAK4PUPPIJetE_noPU = 0;
 evt_noPU->SetBranchAddress( "AK4PUPPIJetE", &myAK4PUPPIJetE_noPU );

 vector<float> *myAK4PUPPIJetRawE_noPU = 0;
 evt_noPU->SetBranchAddress( "AK4PUPPIJetRawE", &myAK4PUPPIJetRawE_noPU );

 UInt_t mynAK8PUPPIJets_noPU = 0;
 evt_noPU->SetBranchAddress( "nAK8PUPPIJets", &mynAK8PUPPIJets_noPU );

 vector<float> *myAK8PUPPIJetPt_noPU = 0;
 evt_noPU->SetBranchAddress( "AK8PUPPIJetPt", &myAK8PUPPIJetPt_noPU );

 vector<float> *myAK8PUPPIJetRawPt_noPU = 0;
 evt_noPU->SetBranchAddress( "AK8PUPPIJetRawPt", &myAK8PUPPIJetRawPt_noPU );

 vector<float> *myAK8PUPPIJetEta_noPU = 0;
 evt_noPU->SetBranchAddress( "AK8PUPPIJetEta", &myAK8PUPPIJetEta_noPU );

 vector<float> *myAK8PUPPIJetPhi_noPU = 0;
 evt_noPU->SetBranchAddress( "AK8PUPPIJetPhi", &myAK8PUPPIJetPhi_noPU );

 vector<float> *myAK8PUPPIJetE_noPU = 0;
 evt_noPU->SetBranchAddress( "AK8PUPPIJetE", &myAK8PUPPIJetE_noPU );

 vector<float> *myAK8PUPPIJetRawE_noPU = 0;
 evt_noPU->SetBranchAddress( "AK8PUPPIJetRawE", &myAK8PUPPIJetRawE_noPU );

 UInt_t mynAK4CHSJets_noPU = 0;
 evt_noPU->SetBranchAddress( "nAK4CHSJets", &mynAK4CHSJets_noPU );

 vector<float> *myAK4CHSJetPt_noPU = 0;
 evt_noPU->SetBranchAddress( "AK4CHSJetPt", &myAK4CHSJetPt_noPU );

 vector<float> *myAK4CHSJetRawPt_noPU = 0;
 evt_noPU->SetBranchAddress( "AK4CHSJetRawPt", &myAK4CHSJetRawPt_noPU );

 vector<float> *myAK4CHSJetEta_noPU = 0;
 evt_noPU->SetBranchAddress( "AK4CHSJetEta", &myAK4CHSJetEta_noPU );

 vector<float> *myAK4CHSJetPhi_noPU = 0;
 evt_noPU->SetBranchAddress( "AK4CHSJetPhi", &myAK4CHSJetPhi_noPU );

 vector<float> *myAK4CHSJetE_noPU = 0;
 evt_noPU->SetBranchAddress( "AK4CHSJetE", &myAK4CHSJetE_noPU );

 vector<float> *myAK4CHSJetRawE_noPU = 0;
 evt_noPU->SetBranchAddress( "AK4CHSJetRawE", &myAK4CHSJetRawE_noPU );

 UInt_t mynGenParticles_noPU = 0;
 evt_noPU->SetBranchAddress( "nGenParticles", &mynGenParticles_noPU );

 vector<float> *mygenParticlePt_noPU = 0;
 evt_noPU->SetBranchAddress( "genParticlePt", &mygenParticlePt_noPU );

 vector<float> *mygenParticleEta_noPU = 0;
 evt_noPU->SetBranchAddress( "genParticleEta", &mygenParticleEta_noPU );

 vector<float> *mygenParticlePhi_noPU = 0;
 evt_noPU->SetBranchAddress( "genParticlePhi", &mygenParticlePhi_noPU );

 vector<float> *mygenParticleE_noPU = 0;
 evt_noPU->SetBranchAddress( "genParticleE", &mygenParticleE_noPU );

 vector<float> *mygenParticleCharge_noPU = 0;
 evt_noPU->SetBranchAddress( "genParticleCharge", &mygenParticleCharge_noPU );

 vector<float> *mygenParticlepdgId_noPU = 0;
 evt_noPU->SetBranchAddress( "genParticlepdgId", &mygenParticlepdgId_noPU );

 float myCHSMET_noPU = 0;
 evt_noPU->SetBranchAddress( "CHSMET", &myCHSMET_noPU );

 float myCHSUnclMET_noPU = 0;
 evt_noPU->SetBranchAddress( "CHSUnclMET", &myCHSUnclMET_noPU );

 float myRawCHSMET_noPU = 0;
 evt_noPU->SetBranchAddress( "RawCHSMET", &myRawCHSMET_noPU );

 float myRawCHSUnclMET_noPU = 0;
 evt_noPU->SetBranchAddress( "RawCHSUnclMET", &myRawCHSUnclMET_noPU );

 float myPUPPIMET_noPU = 0;
 evt_noPU->SetBranchAddress( "PUPPIMET", &myPUPPIMET_noPU );

 float myPUPPIUnclMET_noPU = 0;
 evt_noPU->SetBranchAddress( "PUPPIUnclMET", &myPUPPIUnclMET_noPU );

 float myRawPUPPIMET_noPU = 0;
 evt_noPU->SetBranchAddress( "RawPUPPIMET", &myRawPUPPIMET_noPU );

 float myRawPUPPIUnclMET_noPU = 0;
 evt_noPU->SetBranchAddress( "RawPUPPIUnclMET", &myRawPUPPIUnclMET_noPU );

 float mygenMET = 0; //take genMET from noPU sample, it is the same as in the PU sample
 evt_noPU->SetBranchAddress( "genMET", &mygenMET );

 float mygenUnclMET = 0; //take genUnclMET from noPU sample, it is the same as in the PU sample
 evt_noPU->SetBranchAddress( "genUnclMET", &mygenUnclMET );

 float myVBFDijetGenMass = 0; //take VBFDijetGenMass from noPU sample
 evt_noPU->SetBranchAddress( "VBFDijetGenMass", &myVBFDijetGenMass );

 float myVBFDijetGenMass_noPU = 0; //take VBFDijetGenMass from noPU sample
 evt_noPU->SetBranchAddress( "VBFDijetGenMass_noPU", &myVBFDijetGenMass_noPU );

 float myVBFDijetPUPPIMass_noPU = 0;
 evt_noPU->SetBranchAddress( "VBFDijetPUPPIMass", &myVBFDijetPUPPIMass_noPU );

 float myVBFDijetCHSMass_noPU = 0;
 evt_noPU->SetBranchAddress( "VBFDijetCHSMass", &myVBFDijetCHSMass_noPU );

 Long64_t nevents_noPU = evt_noPU->GetEntries();





 TTree *evt_PU = (TTree*)inputfile_PU->Get( "events" );

 UInt_t myrunNo_PU = 0;
 evt_PU->SetBranchAddress( "runNo", &myrunNo_PU );

 UInt_t myevtNo_PU = 0;
 evt_PU->SetBranchAddress( "evtNo", &myevtNo_PU);

 UInt_t mylumiSec_PU = 0;
 evt_PU->SetBranchAddress( "lumiSec", &mylumiSec_PU );

 UInt_t mynPUint_PU = 0;
 evt_PU->SetBranchAddress( "nPUint", &mynPUint_PU );

 UInt_t mynPFCands_PU = 0;
 evt_PU->SetBranchAddress( "nPFCands", &mynPFCands_PU );

 UInt_t mynLeptons_PU = 0;
 evt_PU->SetBranchAddress( "nLeptons", &mynLeptons_PU );

 vector<float> *myPFCandPt_PU = 0;
 evt_PU->SetBranchAddress( "PFCandPt", &myPFCandPt_PU );

 vector<float> *myPFCandEta_PU = 0;
 evt_PU->SetBranchAddress( "PFCandEta", &myPFCandEta_PU );

 vector<float> *myPFCandAbsEta_PU = 0;
 evt_PU->SetBranchAddress( "PFCandAbsEta", &myPFCandAbsEta_PU );

 vector<float> *myPFCandPhi_PU = 0;
 evt_PU->SetBranchAddress( "PFCandPhi", &myPFCandPhi_PU );

 vector<float> *myPFCandE_PU = 0;
 evt_PU->SetBranchAddress( "PFCandE", &myPFCandE_PU );

 vector<float> *myPFCandpdgId_PU = 0;
 evt_PU->SetBranchAddress( "PFCandpdgId", &myPFCandpdgId_PU );

 vector<float> *myPFCandCharge_PU = 0;
 evt_PU->SetBranchAddress( "PFCandCharge", &myPFCandCharge_PU );

 vector<float> *myPFCandPUPPIw_PU = 0;
 evt_PU->SetBranchAddress( "PFCandPUPPIw", &myPFCandPUPPIw_PU );

 vector<float> *myPFCandPUPPIalpha_PU = 0;
 evt_PU->SetBranchAddress( "PFCandPUPPIalpha", &myPFCandPUPPIalpha_PU );

 vector<float> *myPFCandHCalFrac_PU = 0;
 evt_PU->SetBranchAddress( "PFCandHCalFrac", &myPFCandHCalFrac_PU );

 vector<float> *myPFCandHCalFracCalib_PU = 0;
 evt_PU->SetBranchAddress( "PFCandHCalFracCalib", &myPFCandHCalFracCalib_PU );

 vector<float> *myPFCandVtxAssQual_PU = 0;
 evt_PU->SetBranchAddress( "PFCandVtxAssQual", &myPFCandVtxAssQual_PU );

 vector<float> *myPFCandFromPV_PU = 0;
 evt_PU->SetBranchAddress( "PFCandFromPV", &myPFCandFromPV_PU );

 vector<float> *myPFCandLostInnerHits_PU = 0;
 evt_PU->SetBranchAddress( "PFCandLostInnerHits", &myPFCandLostInnerHits_PU );

 vector<float> *myPFCandTrackHighPurity_PU = 0;
 evt_PU->SetBranchAddress( "PFCandTrackHighPurity", &myPFCandTrackHighPurity_PU );

 vector<float> *myPFCandDZ_PU = 0;
 evt_PU->SetBranchAddress( "PFCandDZ", &myPFCandDZ_PU );

 vector<float> *myPFCandDXY_PU = 0;
 evt_PU->SetBranchAddress( "PFCandDXY", &myPFCandDXY_PU );

 vector<float> *myPFCandDZsig_PU = 0;
 evt_PU->SetBranchAddress( "PFCandDZsig", &myPFCandDZsig_PU );

 vector<float> *myPFCandDXYsig_PU = 0;
 evt_PU->SetBranchAddress( "PFCandDXYsig", &myPFCandDXYsig_PU );

 vector<float> *myPFCandNormChi2_PU = 0;
 evt_PU->SetBranchAddress( "PFCandNormChi2", &myPFCandNormChi2_PU );

 vector<float> *myPFCandQuality_PU = 0;
 evt_PU->SetBranchAddress( "PFCandQuality", &myPFCandQuality_PU );

 vector<float> *myPFCandNumHits_PU = 0;
 evt_PU->SetBranchAddress( "PFCandNumHits", &myPFCandNumHits_PU );

 vector<float> *myPFCandNumLayersHit_PU = 0;
 evt_PU->SetBranchAddress( "PFCandNumLayersHit", &myPFCandNumLayersHit_PU );

 UInt_t mynAK4PUPPIJets_PU = 0;
 evt_PU->SetBranchAddress( "nAK4PUPPIJets", &mynAK4PUPPIJets_PU );

 vector<float> *myAK4PUPPIJetPt_PU = 0;
 evt_PU->SetBranchAddress( "AK4PUPPIJetPt", &myAK4PUPPIJetPt_PU );

 vector<float> *myAK4PUPPIJetRawPt_PU = 0;
 evt_PU->SetBranchAddress( "AK4PUPPIJetRawPt", &myAK4PUPPIJetRawPt_PU );

 vector<float> *myAK4PUPPIJetEta_PU = 0;
 evt_PU->SetBranchAddress( "AK4PUPPIJetEta", &myAK4PUPPIJetEta_PU );

 vector<float> *myAK4PUPPIJetPhi_PU = 0;
 evt_PU->SetBranchAddress( "AK4PUPPIJetPhi", &myAK4PUPPIJetPhi_PU );

 vector<float> *myAK4PUPPIJetE_PU = 0;
 evt_PU->SetBranchAddress( "AK4PUPPIJetE", &myAK4PUPPIJetE_PU );

 vector<float> *myAK4PUPPIJetRawE_PU = 0;
 evt_PU->SetBranchAddress( "AK4PUPPIJetRawE", &myAK4PUPPIJetRawE_PU );

  UInt_t mynAK8PUPPIJets_PU = 0;
 evt_PU->SetBranchAddress( "nAK8PUPPIJets", &mynAK8PUPPIJets_PU );

 vector<float> *myAK8PUPPIJetPt_PU = 0;
 evt_PU->SetBranchAddress( "AK8PUPPIJetPt", &myAK8PUPPIJetPt_PU );

 vector<float> *myAK8PUPPIJetRawPt_PU = 0;
 evt_PU->SetBranchAddress( "AK8PUPPIJetRawPt", &myAK8PUPPIJetRawPt_PU );

 vector<float> *myAK8PUPPIJetEta_PU = 0;
 evt_PU->SetBranchAddress( "AK8PUPPIJetEta", &myAK8PUPPIJetEta_PU );

 vector<float> *myAK8PUPPIJetPhi_PU = 0;
 evt_PU->SetBranchAddress( "AK8PUPPIJetPhi", &myAK8PUPPIJetPhi_PU );

 vector<float> *myAK8PUPPIJetE_PU = 0;
 evt_PU->SetBranchAddress( "AK8PUPPIJetE", &myAK8PUPPIJetE_PU );

 vector<float> *myAK8PUPPIJetRawE_PU = 0;
 evt_PU->SetBranchAddress( "AK8PUPPIJetRawE", &myAK8PUPPIJetRawE_PU );

 UInt_t mynAK4CHSJets_PU = 0;
 evt_PU->SetBranchAddress( "nAK4CHSJets", &mynAK4CHSJets_PU );

 vector<float> *myAK4CHSJetPt_PU = 0;
 evt_PU->SetBranchAddress( "AK4CHSJetPt", &myAK4CHSJetPt_PU );

 vector<float> *myAK4CHSJetRawPt_PU = 0;
 evt_PU->SetBranchAddress( "AK4CHSJetRawPt", &myAK4CHSJetRawPt_PU );

 vector<float> *myAK4CHSJetEta_PU = 0;
 evt_PU->SetBranchAddress( "AK4CHSJetEta", &myAK4CHSJetEta_PU );

 vector<float> *myAK4CHSJetPhi_PU = 0;
 evt_PU->SetBranchAddress( "AK4CHSJetPhi", &myAK4CHSJetPhi_PU );

 vector<float> *myAK4CHSJetE_PU = 0;
 evt_PU->SetBranchAddress( "AK4CHSJetE", &myAK4CHSJetE_PU );

 vector<float> *myAK4CHSJetRawE_PU = 0;
 evt_PU->SetBranchAddress( "AK4CHSJetRawE", &myAK4CHSJetRawE_PU );

 UInt_t mynAK4GenJets_PU = 0;
 evt_PU->SetBranchAddress( "nAK4GenJets", &mynAK4GenJets_PU );

 vector<float> *myAK4GenJetPt_PU = 0;
 evt_PU->SetBranchAddress( "AK4GenJetPt", &myAK4GenJetPt_PU );

 vector<float> *myAK4GenJetEta_PU = 0;
 evt_PU->SetBranchAddress( "AK4GenJetEta", &myAK4GenJetEta_PU );

 vector<float> *myAK4GenJetPhi_PU = 0;
 evt_PU->SetBranchAddress( "AK4GenJetPhi", &myAK4GenJetPhi_PU );

 vector<float> *myAK4GenJetE_PU = 0;
 evt_PU->SetBranchAddress( "AK4GenJetE", &myAK4GenJetE_PU );

 UInt_t mynAK8GenJets_PU = 0;
 evt_PU->SetBranchAddress( "nAK8GenJets", &mynAK8GenJets_PU );

 vector<float> *myAK8GenJetPt_PU = 0;
 evt_PU->SetBranchAddress( "AK8GenJetPt", &myAK8GenJetPt_PU );

 vector<float> *myAK8GenJetEta_PU = 0;
 evt_PU->SetBranchAddress( "AK8GenJetEta", &myAK8GenJetEta_PU );

 vector<float> *myAK8GenJetPhi_PU = 0;
 evt_PU->SetBranchAddress( "AK8GenJetPhi", &myAK8GenJetPhi_PU );

 vector<float> *myAK8GenJetE_PU = 0;
 evt_PU->SetBranchAddress( "AK8GenJetE", &myAK8GenJetE_PU );

 float myCHSMET_PU = 0;
 evt_PU->SetBranchAddress( "CHSMET", &myCHSMET_PU );

 float myCHSUnclMET_PU = 0;
 evt_PU->SetBranchAddress( "CHSUnclMET", &myCHSUnclMET_PU );

 float myRawCHSMET_PU = 0;
 evt_PU->SetBranchAddress( "RawCHSMET", &myRawCHSMET_PU );

 float myRawCHSUnclMET_PU = 0;
 evt_PU->SetBranchAddress( "RawCHSUnclMET", &myRawCHSUnclMET_PU );

 float myPUPPIMET_PU = 0;
 evt_PU->SetBranchAddress( "PUPPIMET", &myPUPPIMET_PU );

 float myPUPPIUnclMET_PU = 0;
 evt_PU->SetBranchAddress( "PUPPIUnclMET", &myPUPPIUnclMET_PU );

 float myRawPUPPIMET_PU = 0;
 evt_PU->SetBranchAddress( "RawPUPPIMET", &myRawPUPPIMET_PU );

 float myRawPUPPIUnclMET_PU = 0;
 evt_PU->SetBranchAddress( "RawPUPPIUnclMET", &myRawPUPPIUnclMET_PU );

 vector<bool> *mytriggerBit_PU = 0;
 evt_PU->SetBranchAddress( "triggerBit", &mytriggerBit_PU );

 vector<bool> *myAK4CHSJetIsBtag_PU = 0;
 evt_PU->SetBranchAddress( "AK4CHSJetIsBtag",  &myAK4CHSJetIsBtag_PU );

 float myVBFDijetPUPPIMass_PU = 0;
 evt_PU->SetBranchAddress( "VBFDijetPUPPIMass", &myVBFDijetPUPPIMass_PU );

 float myVBFDijetCHSMass_PU = 0;
 evt_PU->SetBranchAddress( "VBFDijetCHSMass", &myVBFDijetCHSMass_PU );

 float myVBFDijetGenMass_PU = 0; //take VBFDijetGenMass from noPU sample
 evt_PU->SetBranchAddress( "VBFDijetGenMass", &myVBFDijetGenMass_PU );

 Long64_t nevents_PU = evt_PU->GetEntries();

  
    for ( Long64_t ievent_noPU = 0; ievent_noPU < nevents_noPU; ++ievent_noPU ){
        
        evt_noPU->GetEntry( ievent_noPU );
        evt_PU->GetEntry( ievent_noPU ); //input files are sorted, so take the same entry for both files, it should correspond to the same event in both files
        
//         if (ievent_noPU > 10 ) break;
            
            if(myrunNo_noPU == myrunNo_PU && myevtNo_noPU == myevtNo_PU && mylumiSec_noPU == mylumiSec_PU) { //make a check on runNo, evtNo, lumiSec, just to be sure.
                
                if (ievent_noPU % 1000 == 0) {
                cout << "ievent = " << ievent_noPU << endl;
                std::cout << "*** FILE NO PU *** runNo: " << myrunNo_noPU << "; evtNo: " << myevtNo_noPU << "; lumiSec: " << mylumiSec_noPU <<endl;
                std::cout << "*** FILE PU *** runNo: " << myrunNo_PU << "; evtNo: " << myevtNo_PU << "; lumiSec: " << mylumiSec_PU <<endl;
                
                }
                
                runNo = myrunNo_noPU;
                evtNo = myevtNo_noPU;
                lumiSec = mylumiSec_noPU;
		nPUint_PU = mynPUint_PU;
		nLeptons_PU = mynLeptons_PU;
                nPFCands_noPU = mynPFCands_noPU;
                nPFCands_PU = mynPFCands_PU;
                nAK4PUPPIJets_PU = mynAK4PUPPIJets_PU;
		nAK4PUPPIJets_noPU = mynAK4PUPPIJets_noPU;
		nAK8PUPPIJets_PU = mynAK8PUPPIJets_PU;
		nAK8PUPPIJets_noPU = mynAK8PUPPIJets_noPU;
		nAK4CHSJets_PU = mynAK4CHSJets_PU;
		nAK4CHSJets_noPU = mynAK4CHSJets_noPU;
		nAK4GenJets_PU = mynAK4GenJets_PU;
		nAK8GenJets_PU = mynAK8GenJets_PU;
		nGenParticles_noPU = mynGenParticles_noPU;
		CHSMET_PU = myCHSMET_PU;
		CHSUnclMET_PU = myCHSUnclMET_PU;
                RawCHSMET_PU = myRawCHSMET_PU;
                RawCHSUnclMET_PU = myRawCHSUnclMET_PU;
		PUPPIMET_PU = myPUPPIMET_PU;
                PUPPIUnclMET_PU = myPUPPIUnclMET_PU;
		RawPUPPIMET_PU = myRawPUPPIMET_PU;
                RawPUPPIUnclMET_PU = myRawPUPPIUnclMET_PU;
		CHSMET_noPU = myCHSMET_noPU;
                CHSUnclMET_noPU = myCHSUnclMET_noPU;
		RawCHSMET_noPU = myRawCHSMET_noPU;
                RawCHSUnclMET_noPU = myRawCHSUnclMET_noPU;
		PUPPIMET_noPU = myPUPPIMET_noPU;
                PUPPIUnclMET_noPU = myPUPPIUnclMET_noPU;
		RawPUPPIMET_noPU = myRawPUPPIMET_noPU;
		RawPUPPIUnclMET_noPU = myRawPUPPIUnclMET_noPU;
		genMET = mygenMET;
		genUnclMET = mygenUnclMET;
		triggerBit_PU = *mytriggerBit_PU;
		AK4CHSJetIsBtag_PU = *myAK4CHSJetIsBtag_PU; 
		VBFDijetGenMass = myVBFDijetGenMass;
		VBFDijetGenMass_PU = myVBFDijetGenMass_PU;
		VBFDijetGenMass_noPU = myVBFDijetGenMass_noPU;
		VBFDijetCHSMass_PU = myVBFDijetCHSMass_PU; 
		VBFDijetCHSMass_noPU = myVBFDijetCHSMass_noPU; 
		VBFDijetPUPPIMass_PU = myVBFDijetPUPPIMass_PU; 
		VBFDijetPUPPIMass_noPU = myVBFDijetPUPPIMass_noPU;

                for (int i = 0; i < mynPFCands_noPU; i++) {
                    
                    PFCandPt_noPU.push_back(myPFCandPt_noPU->at(i));
                    PFCandEta_noPU.push_back(myPFCandEta_noPU->at(i));
                    PFCandPhi_noPU.push_back(myPFCandPhi_noPU->at(i));
                    PFCandE_noPU.push_back(myPFCandE_noPU->at(i));
                    PFCandpdgId_noPU.push_back(myPFCandpdgId_noPU->at(i));
                    PFCandCharge_noPU.push_back(myPFCandCharge_noPU->at(i));
                    PFCandPUPPIw_noPU.push_back(myPFCandPUPPIw_noPU->at(i));
		    PFCandPUPPIalpha_noPU.push_back(myPFCandPUPPIalpha_noPU->at(i));
                    PFCandHCalFrac_noPU.push_back(myPFCandHCalFrac_noPU->at(i));
                    PFCandHCalFracCalib_noPU.push_back(myPFCandHCalFracCalib_noPU->at(i));
                    PFCandVtxAssQual_noPU.push_back(myPFCandVtxAssQual_noPU->at(i));
                    PFCandFromPV_noPU.push_back(myPFCandFromPV_noPU->at(i));
                    PFCandLostInnerHits_noPU.push_back(myPFCandLostInnerHits_noPU->at(i));
                    PFCandTrackHighPurity_noPU.push_back(myPFCandTrackHighPurity_noPU->at(i));
                    PFCandDZ_noPU.push_back(myPFCandDZ_noPU->at(i));
                    PFCandDXY_noPU.push_back(myPFCandDXY_noPU->at(i));
                    PFCandDZsig_noPU.push_back(myPFCandDZsig_noPU->at(i));
                    PFCandDXYsig_noPU.push_back(myPFCandDXYsig_noPU->at(i));
                    PFCandNormChi2_noPU.push_back(myPFCandNormChi2_noPU->at(i));
                    PFCandQuality_noPU.push_back(myPFCandQuality_noPU->at(i));
		    PFCandNumHits_noPU.push_back(myPFCandNumHits_noPU->at(i));
		    PFCandNumLayersHit_noPU.push_back(myPFCandNumLayersHit_noPU->at(i));

		    
                }
                
                for (int i = 0; i < mynPFCands_PU; i++) {
                    
                    PFCandPt_PU.push_back(myPFCandPt_PU->at(i));
                    PFCandEta_PU.push_back(myPFCandEta_PU->at(i));
                    PFCandAbsEta_PU.push_back(myPFCandAbsEta_PU->at(i));
                    PFCandPhi_PU.push_back(myPFCandPhi_PU->at(i));
                    PFCandE_PU.push_back(myPFCandE_PU->at(i));
                    PFCandpdgId_PU.push_back(myPFCandpdgId_PU->at(i));
                    PFCandCharge_PU.push_back(myPFCandCharge_PU->at(i));
                    PFCandPUPPIw_PU.push_back(myPFCandPUPPIw_PU->at(i));
		    PFCandPUPPIalpha_PU.push_back(myPFCandPUPPIalpha_PU->at(i));
                    PFCandHCalFrac_PU.push_back(myPFCandHCalFrac_PU->at(i));
                    PFCandHCalFracCalib_PU.push_back(myPFCandHCalFracCalib_PU->at(i));
                    PFCandVtxAssQual_PU.push_back(myPFCandVtxAssQual_PU->at(i));
                    PFCandFromPV_PU.push_back(myPFCandFromPV_PU->at(i));
                    PFCandLostInnerHits_PU.push_back(myPFCandLostInnerHits_PU->at(i));
                    PFCandTrackHighPurity_PU.push_back(myPFCandTrackHighPurity_PU->at(i));
                    PFCandDZ_PU.push_back(myPFCandDZ_PU->at(i));
                    PFCandDXY_PU.push_back(myPFCandDXY_PU->at(i));
                    PFCandDZsig_PU.push_back(myPFCandDZsig_PU->at(i));
                    PFCandDXYsig_PU.push_back(myPFCandDXYsig_PU->at(i));
                    PFCandNormChi2_PU.push_back(myPFCandNormChi2_PU->at(i));
                    PFCandQuality_PU.push_back(myPFCandQuality_PU->at(i));
		    PFCandNumHits_PU.push_back(myPFCandNumHits_PU->at(i));
		    PFCandNumLayersHit_PU.push_back(myPFCandNumLayersHit_PU->at(i));
                    
                }

		for (int i = 0; i < mynGenParticles_noPU; i++) {
                    
                    genParticlePt_noPU.push_back(mygenParticlePt_noPU->at(i));
                    genParticleEta_noPU.push_back(mygenParticleEta_noPU->at(i));
                    genParticlePhi_noPU.push_back(mygenParticlePhi_noPU->at(i));
                    genParticleE_noPU.push_back(mygenParticleE_noPU->at(i));
                    genParticlepdgId_noPU.push_back(mygenParticlepdgId_noPU->at(i));
                    genParticleCharge_noPU.push_back(mygenParticleCharge_noPU->at(i));
                    
                }
                
                for (int i = 0; i < mynAK4PUPPIJets_PU; i++) {
                    
                    AK4PUPPIJetPt_PU.push_back(myAK4PUPPIJetPt_PU->at(i));
		    AK4PUPPIJetRawPt_PU.push_back(myAK4PUPPIJetRawPt_PU->at(i));
                    AK4PUPPIJetEta_PU.push_back(myAK4PUPPIJetEta_PU->at(i));
                    AK4PUPPIJetPhi_PU.push_back(myAK4PUPPIJetPhi_PU->at(i));
                    AK4PUPPIJetE_PU.push_back(myAK4PUPPIJetE_PU->at(i));
		    AK4PUPPIJetRawE_PU.push_back(myAK4PUPPIJetRawE_PU->at(i));
                    
                }

		for (int i = 0; i < mynAK4PUPPIJets_noPU; i++) {
                    
                    AK4PUPPIJetPt_noPU.push_back(myAK4PUPPIJetPt_noPU->at(i));
		    AK4PUPPIJetRawPt_noPU.push_back(myAK4PUPPIJetRawPt_noPU->at(i));
                    AK4PUPPIJetEta_noPU.push_back(myAK4PUPPIJetEta_noPU->at(i));
                    AK4PUPPIJetPhi_noPU.push_back(myAK4PUPPIJetPhi_noPU->at(i));
                    AK4PUPPIJetE_noPU.push_back(myAK4PUPPIJetE_noPU->at(i));
		    AK4PUPPIJetRawE_noPU.push_back(myAK4PUPPIJetRawE_noPU->at(i));
                    
                }

		for (int i = 0; i < mynAK8PUPPIJets_PU; i++) {
                    
                    AK8PUPPIJetPt_PU.push_back(myAK8PUPPIJetPt_PU->at(i));
		    AK8PUPPIJetRawPt_PU.push_back(myAK8PUPPIJetRawPt_PU->at(i));
                    AK8PUPPIJetEta_PU.push_back(myAK8PUPPIJetEta_PU->at(i));
                    AK8PUPPIJetPhi_PU.push_back(myAK8PUPPIJetPhi_PU->at(i));
                    AK8PUPPIJetE_PU.push_back(myAK8PUPPIJetE_PU->at(i));
		    AK8PUPPIJetRawE_PU.push_back(myAK8PUPPIJetRawE_PU->at(i));
                    
                }

		for (int i = 0; i < mynAK8PUPPIJets_noPU; i++) {
                    
                    AK8PUPPIJetPt_noPU.push_back(myAK8PUPPIJetPt_noPU->at(i));
		    AK8PUPPIJetRawPt_noPU.push_back(myAK8PUPPIJetRawPt_noPU->at(i));
                    AK8PUPPIJetEta_noPU.push_back(myAK8PUPPIJetEta_noPU->at(i));
                    AK8PUPPIJetPhi_noPU.push_back(myAK8PUPPIJetPhi_noPU->at(i));
                    AK8PUPPIJetE_noPU.push_back(myAK8PUPPIJetE_noPU->at(i));
		    AK8PUPPIJetRawE_noPU.push_back(myAK8PUPPIJetRawE_noPU->at(i));
                    
                }

		for (int i = 0; i < mynAK4CHSJets_PU; i++) {
                    
                    AK4CHSJetPt_PU.push_back(myAK4CHSJetPt_PU->at(i));
		    AK4CHSJetRawPt_PU.push_back(myAK4CHSJetRawPt_PU->at(i));
                    AK4CHSJetEta_PU.push_back(myAK4CHSJetEta_PU->at(i));
                    AK4CHSJetPhi_PU.push_back(myAK4CHSJetPhi_PU->at(i));
                    AK4CHSJetE_PU.push_back(myAK4CHSJetE_PU->at(i));
		    AK4CHSJetRawE_PU.push_back(myAK4CHSJetRawE_PU->at(i));
                    
                }

		for (int i = 0; i < mynAK4CHSJets_noPU; i++) {
                    
                    AK4CHSJetPt_noPU.push_back(myAK4CHSJetPt_noPU->at(i));
		    AK4CHSJetRawPt_noPU.push_back(myAK4CHSJetRawPt_noPU->at(i));
                    AK4CHSJetEta_noPU.push_back(myAK4CHSJetEta_noPU->at(i));
                    AK4CHSJetPhi_noPU.push_back(myAK4CHSJetPhi_noPU->at(i));
                    AK4CHSJetE_noPU.push_back(myAK4CHSJetE_noPU->at(i));
		    AK4CHSJetRawE_noPU.push_back(myAK4CHSJetRawE_noPU->at(i));
                    
                }

		for (int i = 0; i < mynAK4GenJets_PU; i++) {
                    
                    AK4GenJetPt_PU.push_back(myAK4GenJetPt_PU->at(i));
                    AK4GenJetEta_PU.push_back(myAK4GenJetEta_PU->at(i));
                    AK4GenJetPhi_PU.push_back(myAK4GenJetPhi_PU->at(i));
		    AK4GenJetE_PU.push_back(myAK4GenJetE_PU->at(i));
                    
                }

		for (int i = 0; i < mynAK8GenJets_PU; i++) {
                    
                    AK8GenJetPt_PU.push_back(myAK8GenJetPt_PU->at(i));
                    AK8GenJetEta_PU.push_back(myAK8GenJetEta_PU->at(i));
                    AK8GenJetPhi_PU.push_back(myAK8GenJetPhi_PU->at(i));
		    AK8GenJetE_PU.push_back(myAK8GenJetE_PU->at(i));
                    
                }
                
            }
            
            flatTree->Fill();
            PFCandPt_noPU.clear();
            PFCandEta_noPU.clear();
            PFCandPhi_noPU.clear();
            PFCandE_noPU.clear();
            PFCandPt_PU.clear();
            PFCandEta_PU.clear();
            PFCandAbsEta_PU.clear();
            PFCandPhi_PU.clear();
            PFCandE_PU.clear();
            PFCandpdgId_PU.clear();
            PFCandCharge_PU.clear();
            PFCandPUPPIw_PU.clear();
	    PFCandPUPPIalpha_PU.clear();
            PFCandHCalFrac_PU.clear();
            PFCandHCalFracCalib_PU.clear();
            PFCandVtxAssQual_PU.clear();
            PFCandFromPV_PU.clear();
            PFCandLostInnerHits_PU.clear();
            PFCandTrackHighPurity_PU.clear();
            PFCandDZ_PU.clear();
            PFCandDXY_PU.clear();
            PFCandDZsig_PU.clear();
            PFCandDXYsig_PU.clear();
            PFCandNormChi2_PU.clear();
            PFCandQuality_PU.clear();
	    PFCandNumHits_PU.clear();
	    PFCandNumLayersHit_PU.clear();
	    PFCandpdgId_noPU.clear();
            PFCandCharge_noPU.clear();
            PFCandPUPPIw_noPU.clear();
	    PFCandPUPPIalpha_noPU.clear();
            PFCandHCalFrac_noPU.clear();
            PFCandHCalFracCalib_noPU.clear();
            PFCandVtxAssQual_noPU.clear();
            PFCandFromPV_noPU.clear();
            PFCandLostInnerHits_noPU.clear();
            PFCandTrackHighPurity_noPU.clear();
            PFCandDZ_noPU.clear();
            PFCandDXY_noPU.clear();
            PFCandDZsig_noPU.clear();
            PFCandDXYsig_noPU.clear();
            PFCandNormChi2_noPU.clear();
            PFCandQuality_noPU.clear();
	    PFCandNumHits_noPU.clear();
	    PFCandNumLayersHit_noPU.clear();
            AK4PUPPIJetPt_PU.clear();
	    AK4PUPPIJetRawPt_PU.clear();
            AK4PUPPIJetEta_PU.clear();
            AK4PUPPIJetPhi_PU.clear();
            AK4PUPPIJetE_PU.clear();
	    AK4PUPPIJetRawE_PU.clear();
	    AK4PUPPIJetPt_noPU.clear();
	    AK4PUPPIJetRawPt_noPU.clear();
            AK4PUPPIJetEta_noPU.clear();
            AK4PUPPIJetPhi_noPU.clear();
            AK4PUPPIJetE_noPU.clear();
	    AK4PUPPIJetRawE_noPU.clear();
	    AK8PUPPIJetPt_PU.clear();
	    AK8PUPPIJetRawPt_PU.clear();
            AK8PUPPIJetEta_PU.clear();
            AK8PUPPIJetPhi_PU.clear();
            AK8PUPPIJetE_PU.clear();
	    AK8PUPPIJetRawE_PU.clear();
	    AK8PUPPIJetPt_noPU.clear();
	    AK8PUPPIJetRawPt_noPU.clear();
            AK8PUPPIJetEta_noPU.clear();
            AK8PUPPIJetPhi_noPU.clear();
            AK8PUPPIJetE_noPU.clear();
	    AK8PUPPIJetRawE_noPU.clear();
	    AK4CHSJetPt_PU.clear();
	    AK4CHSJetRawPt_PU.clear();
            AK4CHSJetEta_PU.clear();
            AK4CHSJetPhi_PU.clear();
            AK4CHSJetE_PU.clear();
	    AK4CHSJetRawE_PU.clear();
	    AK4CHSJetPt_noPU.clear();
	    AK4CHSJetRawPt_noPU.clear();
            AK4CHSJetEta_noPU.clear();
            AK4CHSJetPhi_noPU.clear();
            AK4CHSJetE_noPU.clear();
	    AK4CHSJetRawE_noPU.clear();
	    AK4GenJetPt_PU.clear();
            AK4GenJetEta_PU.clear();
            AK4GenJetPhi_PU.clear();
	    AK4GenJetE_PU.clear();
	    AK8GenJetPt_PU.clear();
            AK8GenJetEta_PU.clear();
            AK8GenJetPhi_PU.clear();
	    AK8GenJetE_PU.clear();
	    genParticlePt_noPU.clear();
            genParticleEta_noPU.clear();
            genParticlePhi_noPU.clear();
            genParticleE_noPU.clear();
	    genParticleCharge_noPU.clear();
	    genParticlepdgId_noPU.clear();
	    triggerBit_PU.clear();
	    AK4CHSJetIsBtag_PU.clear();
        
    }
    
flatTree->Write();
outputfile->Write();
gBenchmark->Show("running time");
outputfile->Close();
delete outputfile;

}
