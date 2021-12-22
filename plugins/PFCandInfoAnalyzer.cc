// Package:    PFCandInfo/PFCandInfoAnalyzer
// Class:      PFCandInfoAnalyzer
// 
/**\class PFCandInfoAnalyzer PFCandInfoAnalyzer.cc PFCandInfo/PFCandInfoAnalyzer/plugins/PFCandInfoAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Fabio Iemmi - University and INFN of Bologna
//         Created:  Fri, 20 Mar 2020 17:44:41 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/PatCandidates/interface/PFParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "DataFormats/PatCandidates/interface/MET.h"

//ROOT includes
#include "TH1.h"
#include "TTree.h"
#include <vector>
#include <TLorentzVector.h>
#include "TBranch.h"

//
// class declaration
//

class PFCandInfoAnalyzer : public edm::EDAnalyzer {
   public:
      explicit PFCandInfoAnalyzer(const edm::ParameterSet&);
      ~PFCandInfoAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
   
  // ----------member data ---------------------------
  TTree *outTree_;
  int nvtx;
  int nPUint = 0;
  int nAK4PUPPIJets = 0;
  int nAK4CHSJets = 0;
  int nAK4GenJets = 0;
  int nLeptons;
  int nPFCand;
  int nGenParticles;
  int run, evt, lumi;
  float CHSMET, CHSUnclusteredMET, RawCHSMET, RawCHSUnclusteredMET, PUPPIMET, PUPPIUnclusteredMET, RawPUPPIMET, RawPUPPIUnclusteredMET, genMET, genUnclusteredMET, VBFDijetCHSMass, VBFDijetPUPPIMass, VBFDijetGenMass;
  //next line is for debugging purposes
  std::vector <float> AK4PUPPIJetPt_fromConstituents, AK4PUPPIJetEta_fromConstituents, AK4PUPPIJetPhi_fromConstituents, AK4PUPPIJetE_fromConstituents; 
  std::vector <float> PFCandPt,PFCandPx, PFCandPy, PFCandPz, PFCandEta, PFCandAbsEta, PFCandPhi, PFCandE, PFCandpdgId, PFCandCharge, PFCandPUPPIw, PFCandPUPPIalpha, PFCandHCalFrac,
  PFCandHCalFracCalib, PFCandVtxAssQual, PFCandFromPV, PFCandLostInnerHits, PFCandTrackHighPurity, PFCandDZ, PFCandDXY, PFCandDZSig, PFCandDXYSig, PFCandNormChi2,
  PFCandQuality, PFCandNumHits, PFCandNumPixelHits, PFCandPixelLayersWithMeasurement, PFCandStripLayersWithMeasurement, PFCandTrackerLayersWithMeasurement, AK4PUPPIJetPt,
    AK4PUPPIJetEta, AK4PUPPIJetPhi, AK4PUPPIJetE, AK4PUPPIJetRawPt, AK4PUPPIJetRawE, AK4CHSJetPt, AK4CHSJetEta, AK4CHSJetPhi, AK4CHSJetE, AK4CHSJetRawPt, AK4CHSJetRawE, AK4GenJetPt, AK4GenJetEta, AK4GenJetPhi, AK4GenJetE, genParticlePt, genParticleEta, genParticlePhi, genParticleE, genParticlepdgId, genParticleCharge;
  
  std::vector<std::string> triggerNames_;
  std::string btaggerCSVv2_;
  std::string btaggerDeepCSV_;
  double DeepCSVWP_;
  TH1F* triggerNamesHisto_;
  std::vector<bool> triggerBit_;
  std::vector<bool> AK4CHSJetIsBtag;

  bool isMC;
  
  edm::Service<TFileService> fs_;
  edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > PUToken_;
  edm::EDGetTokenT<pat::JetCollection> PUPPIjetToken_;
  edm::EDGetTokenT<pat::JetCollection> CHSjetToken_;
  edm::EDGetTokenT<reco::GenJetCollection> GenjetToken_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> GenParticleToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfcandToken_;
  edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
  edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  edm::EDGetTokenT<pat::METCollection> METToken_;
  edm::EDGetTokenT<pat::METCollection> PUPPIMETToken_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
PFCandInfoAnalyzer::PFCandInfoAnalyzer(const edm::ParameterSet& iConfig) :
  triggerResultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  PUToken_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("PUinfo"))),
  PUPPIjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("AK4PUPPIJets"))),
  CHSjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("AK4CHSJets"))),
  GenjetToken_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("AK4GenJets"))),
  GenParticleToken_(consumes<pat::PackedGenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))),
  pfcandToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("PFCands"))),
  electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  METToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("missingEt"))),
  PUPPIMETToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("PUPPImissingEt")))
  
{
  
  isMC = iConfig.getUntrackedParameter<bool>("runOnMC");
  
  //now do what ever initialization is needed
  triggerNames_ = iConfig.getParameter<std::vector<std::string> > ("triggerNames");
  btaggerCSVv2_ = iConfig.getParameter<std::string> ("btaggerCSVv2");
  btaggerDeepCSV_ = iConfig.getParameter<std::string> ("btaggerDeepCSV");
  DeepCSVWP_ = iConfig.getParameter<double> ("DeepCSVWP");

  //--- booking the triggerNames histogram ---------                                                                                                            
  triggerNamesHisto_ = new TH1F("TriggerNames","TriggerNames",1,0,1);
  triggerNamesHisto_->SetCanExtend(TH1::kAllAxes);
  
  for(unsigned i = 0; i < triggerNames_.size(); i++) {

    triggerNamesHisto_->Fill(triggerNames_[i].c_str(), 1);
  
  }
  

  outTree_ = new TTree ("events","events");
  outTree_->Branch("runNo",&run,"run/i");
  outTree_->Branch("evtNo",&evt,"evt/i");
  outTree_->Branch("lumiSec",&lumi,"lumi/i");
  outTree_->Branch("triggerBit", &triggerBit_);
  outTree_->Branch("nVtx",&nvtx,"nvtx/i");
  outTree_->Branch("nPUint",&nPUint,"nPUint/i");
  outTree_->Branch("nLeptons", &nLeptons,"nLeptons/i");
  outTree_->Branch("nAK4PUPPIJets", &nAK4PUPPIJets,"nAK4PUPPIJets/i");
  outTree_->Branch("nAK4CHSJets", &nAK4CHSJets,"nAK4CHSJets/i");
  outTree_->Branch("nAK4GenJets", &nAK4GenJets,"nAK4GenJets/i");
  outTree_->Branch("nPFCands", &nPFCand,"nPFCand/i");
  outTree_->Branch("PFCandPt", &PFCandPt);
  //outTree_->Branch("PFCandPx", &PFCandPx);
  //outTree_->Branch("PFCandPy", &PFCandPy);
  //outTree_->Branch("PFCandPz", &PFCandPz);
  outTree_->Branch("PFCandEta", &PFCandEta);
  outTree_->Branch("PFCandAbsEta", &PFCandAbsEta);
  outTree_->Branch("PFCandPhi", &PFCandPhi);
  outTree_->Branch("PFCandE", &PFCandE);
  outTree_->Branch("PFCandpdgId", &PFCandpdgId);
  outTree_->Branch("PFCandCharge", &PFCandCharge);
  outTree_->Branch("PFCandPUPPIw", &PFCandPUPPIw);
  outTree_->Branch("PFCandPUPPIalpha", &PFCandPUPPIalpha);
  outTree_->Branch("PFCandHCalFrac", &PFCandHCalFrac);
  outTree_->Branch("PFCandHCalFracCalib", &PFCandHCalFracCalib);
  outTree_->Branch("PFCandVtxAssQual", &PFCandVtxAssQual);
  outTree_->Branch("PFCandFromPV", &PFCandFromPV);
  outTree_->Branch("PFCandLostInnerHits", &PFCandLostInnerHits);
  outTree_->Branch("PFCandTrackHighPurity", &PFCandTrackHighPurity);
  outTree_->Branch("PFCandDZ", &PFCandDZ);
  outTree_->Branch("PFCandDXY", &PFCandDXY);
  outTree_->Branch("PFCandDZsig", &PFCandDZSig);
  outTree_->Branch("PFCandDXYsig", &PFCandDXYSig);
  outTree_->Branch("PFCandNormChi2", &PFCandNormChi2);
  outTree_->Branch("PFCandQuality", &PFCandQuality);
  //outTree_->Branch("PFCandNumPixelHits", &PFCandNumPixelHits);
  outTree_->Branch("PFCandNumHits", &PFCandNumHits);
  //outTree_->Branch("PFCandPixelLayersWithMeasurement", &PFCandPixelLayersWithMeasurement);
  //outTree_->Branch("PFCandStripLayersWithMeasurement", &PFCandStripLayersWithMeasurement);
  outTree_->Branch("PFCandNumLayersHit", &PFCandTrackerLayersWithMeasurement);
  outTree_->Branch("nGenParticles", &nGenParticles,"nGenParticles/i");
  outTree_->Branch("genParticlePt", &genParticlePt);
  outTree_->Branch("genParticleEta", &genParticleEta);
  outTree_->Branch("genParticlePhi", &genParticlePhi);
  outTree_->Branch("genParticleE", &genParticleE);
  outTree_->Branch("genParticlepdgId", &genParticlepdgId);
  outTree_->Branch("genParticleCharge", &genParticleCharge);
  outTree_->Branch("AK4PUPPIJetPt", &AK4PUPPIJetPt);
  outTree_->Branch("AK4PUPPIJetEta", &AK4PUPPIJetEta);
  outTree_->Branch("AK4PUPPIJetPhi", &AK4PUPPIJetPhi);
  outTree_->Branch("AK4PUPPIJetE", &AK4PUPPIJetE);
  outTree_->Branch("AK4PUPPIJetRawPt", &AK4PUPPIJetRawPt);
  outTree_->Branch("AK4PUPPIJetRawE", &AK4PUPPIJetRawE);
  outTree_->Branch("AK4PUPPIJetPt_fromConstituents", &AK4PUPPIJetPt_fromConstituents);
  outTree_->Branch("AK4PUPPIJetEta_fromConstituents", &AK4PUPPIJetEta_fromConstituents);
  outTree_->Branch("AK4PUPPIJetPhi_fromConstituents", &AK4PUPPIJetPhi_fromConstituents);
  outTree_->Branch("AK4PUPPIJetE_fromConstituents", &AK4PUPPIJetE_fromConstituents);
  outTree_->Branch("AK4CHSJetPt", &AK4CHSJetPt);
  outTree_->Branch("AK4CHSJetEta", &AK4CHSJetEta);
  outTree_->Branch("AK4CHSJetPhi", &AK4CHSJetPhi);
  outTree_->Branch("AK4CHSJetE", &AK4CHSJetE);
  outTree_->Branch("AK4CHSJetRawPt", &AK4CHSJetRawPt);
  outTree_->Branch("AK4CHSJetRawE", &AK4CHSJetRawE);
  outTree_->Branch("AK4CHSJetIsBtag", &AK4CHSJetIsBtag);
  outTree_->Branch("AK4GenJetPt", &AK4GenJetPt);
  outTree_->Branch("AK4GenJetEta", &AK4GenJetEta);
  outTree_->Branch("AK4GenJetPhi", &AK4GenJetPhi);
  outTree_->Branch("AK4GenJetE", &AK4GenJetE);
  outTree_->Branch("CHSMET", &CHSMET);
  outTree_->Branch("CHSUnclMET", &CHSUnclusteredMET);
  outTree_->Branch("RawCHSMET", &RawCHSMET);
  outTree_->Branch("RawCHSUnclMET", &RawCHSUnclusteredMET);
  outTree_->Branch("PUPPIMET", &PUPPIMET);
  outTree_->Branch("PUPPIUnclMET", &PUPPIUnclusteredMET);
  outTree_->Branch("RawPUPPIMET", &RawPUPPIMET);
  outTree_->Branch("RawPUPPIUnclMET", &RawPUPPIUnclusteredMET);
  outTree_->Branch("genMET", &genMET);
  outTree_->Branch("genUnclMET", &genUnclusteredMET);
  outTree_->Branch("VBFDijetCHSMass", &VBFDijetCHSMass);
  outTree_->Branch("VBFDijetPUPPIMass", &VBFDijetPUPPIMass);
  outTree_->Branch("VBFDijetGenMass", &VBFDijetGenMass);

  triggerNamesHisto_->Write();

}


PFCandInfoAnalyzer::~PFCandInfoAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
PFCandInfoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;

  run    = iEvent.id().run();
  evt    = iEvent.id().event();
  lumi   = iEvent.id().luminosityBlock();
  //cout << run << ":" << lumi << ":" << evt << endl;
  //-------------- Handle vertices info ----------------------------------------------------------------------------------------------------------------------

  Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  
  if (vertices->empty()) {

    std::cout << "vertices is empty!!! skipping event..." << endl;
    return;
  
  }
  
  nvtx = vertices->size();

   //-------------- If Monte Carlo, handle pileup info ----------------------------------------------------------------------------------------------------------------------

  if (isMC) {
        
    Handle <std::vector<PileupSummaryInfo> > PUinfo;
    iEvent.getByToken(PUToken_, PUinfo);
  
    if (PUinfo->empty()) {
  
      std::cout << "PUinfo is empty!!! skipping event..." << endl;
      return;
  
    }

    std::vector<PileupSummaryInfo>::const_iterator PVI;

    for(PVI = PUinfo->begin(); PVI != PUinfo->end(); ++PVI) {

      int BX = PVI->getBunchCrossing();

      if(BX == 0) { 
     
	nPUint = PVI->getPU_NumInteractions();
    
      }

    }
 
  }
  
  //-------------- Handle trigger info ----------------------------------------------------------------------------------------------------------------------
 
  Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(triggerResultsToken_, triggerResults);
  
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerResults);
  
  for(unsigned int k = 0; k < triggerNames_.size(); k++) {//loop over the list of trigger names...
   
    bool bit(false);
    
    for(unsigned int itrig = 0; itrig < triggerResults->size(); itrig++) {//loop over triggerResults to see what triggers the event passed
      
      string trigger_name = string(names.triggerName(itrig)); //get the name of the itrig-th trigger that has been passed
      
      //--- erase the last characters, i.e. the version number---
      for (int a = 0; a < 2; a++) {

	char last_char = trigger_name.back();
	if ( isdigit( last_char ) ) trigger_name.pop_back();
      
      }

      if (trigger_name == triggerNames_[k]) {//if the name is in the list of trigger names...
	
	bit = triggerResults->accept(itrig);//...bit becomes true
     
      }
   
    }

    triggerBit_.push_back(bit);

  }

  //----------- Handle lepton info ----------------------------------------------------------------------------------------------------------------------------

  Handle<pat::ElectronCollection> electrons;
  iEvent.getByToken(electronToken_, electrons);

  Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);

  vector<const reco::Candidate *> myLeptons;
  
  for ( long unsigned int i = 0; i < muons->size(); i++ ) { //loop on muons imposing cuts to reject non-prompt muons: hadronic punch-throughs, in-flight decays from b/c-flavored hadrons
    
    if ( muons->at(i).passed(reco::Muon::PFIsoMedium) ) {//requre  an isolated muon
      /*
	+++ NOTE: SKIP THE QUALITY REQUIREMENTS FOR MUONS: WANT THE NON-PROMPT LEPTONS TO BE IN THE COUNT TO VETO THEM +++

      if ( muons->at(i).isGlobalMuon() && //track identified as a global muon
	   muons->at(i).globalTrack()->normalizedChi2() < 10.0 && //chi2/ndof of the global fit less than 10
	   muons->at(i).globalTrack()->hitPattern().numberOfValidMuonHits() > 0 && //at least one hit in muon detector
	   fabs(muons->at(i).globalTrack()->dxy(vertices->at(0).position())) < 0.002 && // transverse impact parameter less than 2mm
 	   muons->at(i).globalTrack()->numberOfValidHits() > 10 ) //more than ten valid hits on the track
      */
      myLeptons.push_back(&muons->at(i));
    
    }
  
  }
  
  for ( long unsigned int i = 0; i < electrons->size(); i++ ) {
    
    float eleRelIso = ( electrons->at(i).dr03TkSumPt() + max( 0., electrons->at(i).dr03EcalRecHitSumEt() -1.) + electrons->at(i).dr03HcalTowerSumEt() )/electrons->at(i).pt();
    
    if ( eleRelIso < 0.20) myLeptons.push_back(&electrons->at(i));
  
  }
  
  std::sort(myLeptons.begin(),myLeptons.end(),[](const reco::Candidate *a,const reco::Candidate *b){return a->pt() > b->pt();});
  
  nLeptons = (int)myLeptons.size();

   //-------------- Handle AK4 PUPPI jets info ----------------------------------------------------------------------------------------------------------------------
  
  Handle<pat::JetCollection> AK4PUPPIJets;
  iEvent.getByToken(PUPPIjetToken_, AK4PUPPIJets);
  
  
  /*if(AK4PUPPIJets->empty()) {
    
    std::cout << "AK4PUPPIJets is empty!!! skipping event..." << endl;
    return;
    
    }*/

  //nAK4PUPPIJets = AK4PUPPIJets->size(); count only jets with pt > 20 and |eta| < 2.4

   //-------------- Handle AK4 CHS jets info ----------------------------------------------------------------------------------------------------------------------

  Handle<pat::JetCollection> AK4CHSJets;
  iEvent.getByToken(CHSjetToken_, AK4CHSJets);
  
  /*if(AK4CHSJets->empty()) {
    
    std::cout << "AK4CHSJets is empty!!! skipping event..." << endl;
    return;
  
    }
  */
  //nAK4CHSJets = AK4CHSJets->size(); count only jets with pt > 20 and |eta| < 2.4

 //-------------- Handle particle-flow candidates info ----------------------------------------------------------------------------------------------------------------------

  Handle<pat::PackedCandidateCollection> PFCands;
  iEvent.getByToken(pfcandToken_, PFCands);
  
  if(PFCands->empty()) {
    
    std::cout << "PFCands is empty!!! skipping event..." << endl;
    return;
  
  }

  nPFCand = PFCands->size();
  
  //loop on PF candidates to store their Pt, Eta, Phi, E, pdgId, charge, PUPPI weight
  for (int i = 0; i < nPFCand; i++) {
    
    float alpha = 0.0;

    //compute PUPPIalpha
    alpha = 0.0;
    
    if(std::abs(PFCands->at(i).eta()) < 2.5) {//PF candidate is inside the tracking volume

      for (int j = 0; j < nPFCand; j++) {
	
	if (!(std::abs(PFCands->at(j).charge()) == 1 && PFCands->at(j).fromPV() > 1)) continue; //in the tracking volume, take only charged PFs from LV in alpha computation

	float deltaR = 0.0;

	if (j != i) {
	
	  float deltaEta = fabs( PFCands->at(i).eta() - PFCands->at(j).eta() );
	  float deltaPhi = fabs( PFCands->at(i).phi() - PFCands->at(j).phi() );
	  if(deltaPhi > 3.14159265358979323846) deltaPhi  = 2*3.14159265358979323846 - deltaPhi;         
	  deltaR = sqrt( pow(deltaEta, 2) + pow(deltaPhi, 2) );
      
	  if (deltaR < 0.4) {

	    alpha += pow(PFCands->at(j).pt()/deltaR ,2);
	
	  }
	  
	}//end j != i
    
      }//end inner loop on PF cands
    
    }//end abs(PFCandEta) < 2.5

    else {//PF candidate is outside the tracking volume

      for (int j = 0; j < nPFCand; j++) {
	
	float deltaR = 0.0;

	if (j != i) {
	
	  float deltaEta = fabs( PFCands->at(i).eta() - PFCands->at(j).eta() );
	  float deltaPhi = fabs( PFCands->at(i).phi() - PFCands->at(j).phi() );
	  if(deltaPhi > 3.14159265358979323846) deltaPhi  = 2*3.14159265358979323846 - deltaPhi;         
	  deltaR = sqrt( pow(deltaEta, 2) + pow(deltaPhi, 2) );
      
	  if (deltaR < 0.4) {

	    alpha += pow(PFCands->at(j).pt()/deltaR ,2);
	
	  }

	}
    
      }//end inner loop on PF cands
    
    }//end abs(PFCandEta) > 2.5   
    

    PFCandPt.push_back(PFCands->at(i).pt());
    PFCandPx.push_back(PFCands->at(i).px());
    PFCandPy.push_back(PFCands->at(i).py());
    PFCandPz.push_back(PFCands->at(i).pz());
    PFCandEta.push_back(PFCands->at(i).eta());
    PFCandAbsEta.push_back(std::abs(PFCands->at(i).eta()));
    PFCandPhi.push_back(PFCands->at(i).phi());
    PFCandE.push_back(PFCands->at(i).energy());
    PFCandpdgId.push_back(PFCands->at(i).pdgId());
    PFCandCharge.push_back(PFCands->at(i).charge());
    PFCandPUPPIw.push_back(PFCands->at(i).puppiWeight());
    PFCandPUPPIalpha.push_back(log(alpha));
    PFCandNumHits.push_back(PFCands->at(i).numberOfHits());
    PFCandNumPixelHits.push_back(PFCands->at(i).numberOfPixelHits());
    PFCandPixelLayersWithMeasurement.push_back(PFCands->at(i).pixelLayersWithMeasurement());
    PFCandStripLayersWithMeasurement.push_back(PFCands->at(i).stripLayersWithMeasurement());
    PFCandTrackerLayersWithMeasurement.push_back(PFCands->at(i).trackerLayersWithMeasurement());
   
    //for neutral hadrons or HF hadron, store fraction of energy recorded in the HCAL
    if (PFCands->at(i).pdgId() == 130 || PFCands->at(i).pdgId() == 1) {
      
      PFCandHCalFrac.push_back(PFCands->at(i).hcalFraction());
      
    }

    else if (PFCands->at(i).isIsolatedChargedHadron()) {

      PFCandHCalFrac.push_back(PFCands->at(i).rawHcalFraction());
    
    }

    else PFCandHCalFrac.push_back(0.0);

    PFCandHCalFracCalib.push_back(PFCands->at(i).rawHcalFraction());

    //for charged hadrons, store quality of PV-candidate association
    PFCandVtxAssQual.push_back(PFCands->at(i).pvAssociationQuality());

    PFCandFromPV.push_back(PFCands->at(i).fromPV());
    PFCandLostInnerHits.push_back(PFCands->at(i).lostInnerHits());
    PFCandTrackHighPurity.push_back(PFCands->at(i).trackHighPurity());
    PFCandDZ.push_back(PFCands->at(i).dz());
    PFCandDXY.push_back(PFCands->at(i).dxy());

    if (PFCands->at(i).bestTrack()) {
    
      PFCandDZSig.push_back(PFCands->at(i).dz()/PFCands->at(i).dzError());
      PFCandDXYSig.push_back(PFCands->at(i).dxy()/PFCands->at(i).dxyError());
      
      const auto *trk = PFCands->at(i).bestTrack();
      PFCandNormChi2.push_back(trk->normalizedChi2());
      PFCandQuality.push_back(trk->qualityMask());
      
    }

    else {

      PFCandDZSig.push_back(0.0);
      PFCandDXYSig.push_back(0.0);
      PFCandNormChi2.push_back(999.0);
      PFCandQuality.push_back(0.0);
      
    }

  }//end loop on PF candidates

  
  //loop on PUPPI jets to store their Pt, Eta, Phi

  TLorentzVector PUPPIJetsTLVector;
  PUPPIJetsTLVector.SetPtEtaPhiE(0,0,0,0);
  TLorentzVector RawPUPPIJetsTLVector;
  RawPUPPIJetsTLVector.SetPtEtaPhiE(0,0,0,0);

  
  for (long unsigned int i = 0; i < AK4PUPPIJets->size(); i++) {

    if ( AK4PUPPIJets->at(i).pt() > 20 && fabs(AK4PUPPIJets->at(i).eta()) < 4.7 ) {

      nAK4PUPPIJets++;
      AK4PUPPIJetPt.push_back(AK4PUPPIJets->at(i).pt());
      AK4PUPPIJetEta.push_back(AK4PUPPIJets->at(i).eta());
      AK4PUPPIJetPhi.push_back(AK4PUPPIJets->at(i).phi());
      AK4PUPPIJetE.push_back(AK4PUPPIJets->at(i).energy());
      AK4PUPPIJetRawPt.push_back(AK4PUPPIJets->at(i).correctedJet("Uncorrected").pt());
      AK4PUPPIJetRawE.push_back(AK4PUPPIJets->at(i).correctedJet("Uncorrected").energy());
      TLorentzVector myPUPPIJetTLVector;
      myPUPPIJetTLVector.SetPtEtaPhiE(AK4PUPPIJets->at(i).pt(), AK4PUPPIJets->at(i).eta(), AK4PUPPIJets->at(i).phi(), AK4PUPPIJets->at(i).energy());
      PUPPIJetsTLVector += myPUPPIJetTLVector;
      TLorentzVector myRawPUPPIJetTLVector;
      myRawPUPPIJetTLVector.SetPtEtaPhiE(AK4PUPPIJets->at(i).correctedJet("Uncorrected").pt(), AK4PUPPIJets->at(i).eta(), AK4PUPPIJets->at(i).phi(), AK4PUPPIJets->at(i).correctedJet("Uncorrected").energy());
      RawPUPPIJetsTLVector += myRawPUPPIJetTLVector;
      
      TLorentzVector myJetConstituents;
      myJetConstituents.SetPtEtaPhiE(0,0,0,0);
      int nconst = AK4PUPPIJets->at(i).numberOfDaughters();
      for (int j = 0; j < nconst; j++) {

	TLorentzVector myJetConstituent;
	pat::PackedCandidate * myConstituent  =  (pat::PackedCandidate*) AK4PUPPIJets->at(i).daughter(j);
	myJetConstituent.SetPtEtaPhiE(myConstituent->pt()*myConstituent->puppiWeight(), myConstituent->eta(), myConstituent->phi(), myConstituent->energy()*myConstituent->puppiWeight());
	myJetConstituents += myJetConstituent;

      }

      AK4PUPPIJetPt_fromConstituents.push_back(myJetConstituents.Pt());
      AK4PUPPIJetEta_fromConstituents.push_back(myJetConstituents.Eta());
      AK4PUPPIJetPhi_fromConstituents.push_back(myJetConstituents.Phi());
      AK4PUPPIJetE_fromConstituents.push_back(myJetConstituents.E());
    }

  }

  // Find pair among puppi jets with large deltaEta (>3.0) and largest dijet mass
  VBFDijetPUPPIMass = -1.;
  for (unsigned int i=0; i<AK4PUPPIJetPt.size(); i++){
    TLorentzVector j1;
    j1.SetPtEtaPhiE(AK4PUPPIJetRawPt[i],AK4PUPPIJetEta[i],AK4PUPPIJetPhi[i],AK4PUPPIJetRawE[i]);
    for (unsigned int j=i+1; j<AK4PUPPIJetPt.size(); j++){
      TLorentzVector j2;
      j2.SetPtEtaPhiE(AK4PUPPIJetRawPt[j],AK4PUPPIJetEta[j],AK4PUPPIJetPhi[j],AK4PUPPIJetRawE[j]);
      if (abs(j1.Eta()-j2.Eta()) < 3.0)
	continue;
      TLorentzVector combined = j1+j2;
      if (combined.M() > VBFDijetPUPPIMass)
	VBFDijetPUPPIMass = combined.M();
    }
  }
  
  //loop on CHS jets to store their Pt, Eta, Phi

  TLorentzVector CHSJetsTLVector;
  CHSJetsTLVector.SetPtEtaPhiE(0,0,0,0);
  TLorentzVector RawCHSJetsTLVector;
  RawCHSJetsTLVector.SetPtEtaPhiE(0,0,0,0);
  for (long unsigned int i = 0; i < AK4CHSJets->size(); i++) {

    if ( AK4CHSJets->at(i).pt() > 20 && fabs(AK4CHSJets->at(i).eta()) < 4.7 ) {

      nAK4CHSJets++;
      
      bool btagged = false;
      if ( AK4CHSJets->at(i).bDiscriminator(btaggerDeepCSV_) > DeepCSVWP_ ) {
      
	btagged = true;

      }

      AK4CHSJetIsBtag.push_back(btagged);

      AK4CHSJetPt.push_back(AK4CHSJets->at(i).pt());
      AK4CHSJetEta.push_back(AK4CHSJets->at(i).eta());
      AK4CHSJetPhi.push_back(AK4CHSJets->at(i).phi());
      AK4CHSJetE.push_back(AK4CHSJets->at(i).energy());
      AK4CHSJetRawPt.push_back(AK4CHSJets->at(i).correctedJet("Uncorrected").pt());
      AK4CHSJetRawE.push_back(AK4CHSJets->at(i).correctedJet("Uncorrected").energy());
      TLorentzVector myCHSJetTLVector;
      myCHSJetTLVector.SetPtEtaPhiE(AK4CHSJets->at(i).pt(), AK4CHSJets->at(i).eta(), AK4CHSJets->at(i).phi(), AK4CHSJets->at(i).energy());
      CHSJetsTLVector += myCHSJetTLVector;
      TLorentzVector myRawCHSJetTLVector;
      myRawCHSJetTLVector.SetPtEtaPhiE(AK4CHSJets->at(i).correctedJet("Uncorrected").pt(), AK4CHSJets->at(i).eta(), AK4CHSJets->at(i).phi(), AK4CHSJets->at(i).correctedJet("Uncorrected").energy());
      RawCHSJetsTLVector += myRawCHSJetTLVector;

    }
   
  }


  // Find pair among chs jets with large deltaEta (>3.0) and largest dijet mass
  VBFDijetCHSMass = -1.;
  for (unsigned int i=0; i<AK4CHSJetPt.size(); i++){
    TLorentzVector j1;
    j1.SetPtEtaPhiE(AK4CHSJetRawPt[i],AK4CHSJetEta[i],AK4CHSJetPhi[i],AK4CHSJetRawE[i]);
    for (unsigned int j=i+1; j<AK4CHSJetPt.size(); j++){
      TLorentzVector j2;
      j2.SetPtEtaPhiE(AK4CHSJetRawPt[j],AK4CHSJetEta[j],AK4CHSJetPhi[j],AK4CHSJetRawE[j]);
      if (abs(j1.Eta()-j2.Eta()) < 3.0)
	continue;
      TLorentzVector combined = j1+j2;
      if (combined.M() > VBFDijetCHSMass)
	VBFDijetCHSMass = combined.M();
    }
  }



  //-------------- If Monte Carlo, handle AK4 gen jets and gen particles info ---------------------------------------------------------------------------------------------------------
  TLorentzVector genJetsTLVector;
  if (isMC) {
    
    Handle<reco::GenJetCollection> AK4GenJets;
    iEvent.getByToken(GenjetToken_, AK4GenJets);
    /*
    if(AK4GenJets->empty()) {
    
      std::cout << "AK4GenJets is empty!!! skipping event..." << endl;
      return;
  
    }
  
    nAK4GenJets = AK4GenJets->size(); count only jets with pt > 20 and |eta| < 2.4
    */
      
    //loop on Gen jets to store their Pt, Eta, Phi
    for (long unsigned int i = 0; i < AK4GenJets->size(); i++) {

      if ( AK4GenJets->at(i).pt() > 20 && fabs(AK4GenJets->at(i).eta()) < 4.7 ) {

	nAK4GenJets++;
	AK4GenJetPt.push_back(AK4GenJets->at(i).pt());
	AK4GenJetEta.push_back(AK4GenJets->at(i).eta());
	AK4GenJetPhi.push_back(AK4GenJets->at(i).phi());
	AK4GenJetE.push_back(AK4GenJets->at(i).energy());
	TLorentzVector mygenJetTLVector;
	mygenJetTLVector.SetPtEtaPhiE(AK4GenJets->at(i).pt(), AK4GenJets->at(i).eta(), AK4GenJets->at(i).phi(), AK4GenJets->at(i).energy());
	genJetsTLVector += mygenJetTLVector;
      }
   
    }

    // Find pair among gen jets with large deltaEta (>3.0) and largest dijet mass
    VBFDijetGenMass = -1.;
    for (unsigned int i=0; i<AK4GenJetPt.size(); i++){
      TLorentzVector j1;
      j1.SetPtEtaPhiE(AK4GenJetPt[i],AK4GenJetEta[i],AK4GenJetPhi[i],AK4GenJetE[i]);
      for (unsigned int j=i+1; j<AK4GenJetPt.size(); j++){
	TLorentzVector j2;
	j2.SetPtEtaPhiE(AK4GenJetPt[j],AK4GenJetEta[j],AK4GenJetPhi[j],AK4GenJetE[j]);
	if (abs(j1.Eta()-j2.Eta()) < 3.0)
	  continue;
	TLorentzVector combined = j1+j2;
	if (combined.M() > VBFDijetGenMass)
	  VBFDijetGenMass = combined.M();
      }
    }


    Handle<pat::PackedGenParticleCollection> GenParticles;
    iEvent.getByToken(GenParticleToken_, GenParticles);
    
    if(GenParticles->empty()) {
    
    std::cout << "GenParticles is empty!!! skipping event..." << endl;
    return;
  
    }

    nGenParticles = GenParticles->size();
  
    //loop on GenParticles to store their Pt, Eta, Phi, E, pdgId, charge
    for (int i = 0; i < nGenParticles; i++) {

      genParticlePt.push_back(GenParticles->at(i).pt());
      genParticleEta.push_back(GenParticles->at(i).eta());
      genParticlePhi.push_back(GenParticles->at(i).phi());
      genParticleE.push_back(GenParticles->at(i).energy());
      genParticlepdgId.push_back(GenParticles->at(i).pdgId());
      genParticleCharge.push_back(GenParticles->at(i).charge());

    }
  
  }
  
  Handle<pat::METCollection> MET;
  iEvent.getByToken(METToken_, MET);
  
  CHSMET = MET->at(0).pt();
  CHSUnclusteredMET = sqrt( pow(MET->at(0).px() + CHSJetsTLVector.Px(), 2) + pow(MET->at(0).py() + CHSJetsTLVector.Py(), 2) );
  RawCHSMET = MET->at(0).uncorPt();
  RawCHSUnclusteredMET = sqrt( pow(MET->at(0).uncorPx() + RawCHSJetsTLVector.Px(), 2) + pow(MET->at(0).uncorPy() + RawCHSJetsTLVector.Py(), 2) );

  Handle<pat::METCollection> PuppiMET;
  iEvent.getByToken(PUPPIMETToken_, PuppiMET);
  
  PUPPIMET = PuppiMET->at(0).pt();
  PUPPIUnclusteredMET = sqrt( pow(PuppiMET->at(0).px() + PUPPIJetsTLVector.Px(), 2) + pow(PuppiMET->at(0).py() + PUPPIJetsTLVector.Py(), 2) );
  RawPUPPIMET = PuppiMET->at(0).uncorPt();
  RawPUPPIUnclusteredMET = sqrt( pow(PuppiMET->at(0).uncorPx() + RawPUPPIJetsTLVector.Px(), 2) + pow(PuppiMET->at(0).uncorPy() + RawPUPPIJetsTLVector.Py(), 2) );

  //if isMC, store genMET
  if (isMC) {

    genMET=MET->at(0).genMET()->pt();
    genUnclusteredMET = sqrt( pow(MET->at(0).genMET()->px() + genJetsTLVector.Px(), 2) + pow(MET->at(0).genMET()->py() + genJetsTLVector.Py(), 2) );

  }

  outTree_->Fill();
 
  nAK4PUPPIJets = 0;
  nAK4CHSJets = 0;
  nAK4GenJets = 0;
  triggerBit_.clear();
  myLeptons.clear();
  PFCandPt.clear();
  PFCandPx.clear();
  PFCandPy.clear();
  PFCandPz.clear();
  PFCandEta.clear();
  PFCandAbsEta.clear();
  PFCandPhi.clear();
  PFCandE.clear();
  PFCandpdgId.clear();
  PFCandCharge.clear();
  PFCandPUPPIw.clear();
  PFCandPUPPIalpha.clear();
  PFCandHCalFrac.clear();
  PFCandHCalFracCalib.clear();
  PFCandVtxAssQual.clear();
  PFCandFromPV.clear();
  PFCandLostInnerHits.clear();
  PFCandTrackHighPurity.clear();
  PFCandDZ.clear();
  PFCandDXY.clear();
  PFCandDZSig.clear();
  PFCandDXYSig.clear();
  PFCandNormChi2.clear();
  PFCandQuality.clear();
  PFCandNumHits.clear();
  PFCandNumPixelHits.clear();
  PFCandPixelLayersWithMeasurement.clear();
  PFCandStripLayersWithMeasurement.clear();
  PFCandTrackerLayersWithMeasurement.clear();
  genParticlePt.clear();
  genParticleEta.clear();
  genParticlePhi.clear();
  genParticleE.clear();
  genParticlepdgId.clear();
  genParticleCharge.clear();
  AK4PUPPIJetPt.clear();
  AK4PUPPIJetEta.clear();
  AK4PUPPIJetPhi.clear();
  AK4PUPPIJetE.clear();
  AK4PUPPIJetRawPt.clear();
  AK4PUPPIJetRawE.clear();
  AK4PUPPIJetPt_fromConstituents.clear();
  AK4PUPPIJetEta_fromConstituents.clear();
  AK4PUPPIJetPhi_fromConstituents.clear();
  AK4PUPPIJetE_fromConstituents.clear();
  AK4CHSJetPt.clear();
  AK4CHSJetEta.clear();
  AK4CHSJetPhi.clear();
  AK4CHSJetE.clear();
  AK4CHSJetRawPt.clear();
  AK4CHSJetRawE.clear();
  AK4CHSJetIsBtag.clear();
  AK4GenJetPt.clear();
  AK4GenJetEta.clear();
  AK4GenJetPhi.clear();
  AK4GenJetE.clear();
  


#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
PFCandInfoAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PFCandInfoAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
PFCandInfoAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
PFCandInfoAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
PFCandInfoAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
PFCandInfoAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PFCandInfoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PFCandInfoAnalyzer);
