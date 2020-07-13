
// Package:    PFCandInfo/PFCandInfoAnalyzer
// Class:      PFCandInfoAnalyzer
// 
/**\class PFCandInfoAnalyzer PFCandInfoAnalyzer.cc PFCandInfo/PFCandInfoAnalyzer/plugins/PFCandInfoAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Fabio Iemmi
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

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/PatCandidates/interface/PFParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

//ROOT includes
#include "TH1.h"
#include "TTree.h"
#include <vector>
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
  int nPUint;
  int nAK4PUPPIJets = 0;
  int nAK4CHSJets = 0;
  int nAK4GenJets = 0;
  int nPFCand;
  int run, evt, lumi;
  std::vector <float> PFCandPt,PFCandPx, PFCandPy, PFCandPz, PFCandEta, PFCandAbsEta, PFCandPhi, PFCandE, PFCandpdgId, PFCandCharge, PFCandPUPPIw, PFCandHCalFrac, PFCandHCalFracCalib, PFCandVtxAssQual, PFCandFromPV, PFCandLostInnerHits, PFCandTrackHighPurity, PFCandDZ, PFCandDXY, PFCandDZSig, PFCandDXYSig, PFCandNormChi2, PFCandQuality, AK4PUPPIJetPt,  AK4PUPPIJetEta, AK4PUPPIJetPhi, AK4PUPPIJetE, AK4PUPPIJetRawPt, AK4PUPPIJetRawE, AK4CHSJetPt, AK4CHSJetEta, AK4CHSJetPhi, AK4CHSJetE, AK4CHSJetRawPt, AK4CHSJetRawE, AK4GenJetPt, AK4GenJetEta, AK4GenJetPhi;
  edm::Service<TFileService> fs_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > PUToken_;
  edm::EDGetTokenT<pat::JetCollection> PUPPIjetToken_;
  edm::EDGetTokenT<pat::JetCollection> CHSjetToken_;
  edm::EDGetTokenT<reco::GenJetCollection> GenjetToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfcandToken_;
  
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
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  PUToken_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("PUinfo"))),
  PUPPIjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("AK4PUPPIJets"))),
  CHSjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("AK4CHSJets"))),
  GenjetToken_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("AK4GenJets"))),
  pfcandToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("PFCands")))
  
{
  //now do what ever initialization is needed
  outTree_ = new TTree ("events","events");
  outTree_->Branch("runNo",&run,"run/i");
  outTree_->Branch("evtNo",&evt,"evt/i");
  outTree_->Branch("lumiSec",&lumi,"lumi/i");
  outTree_->Branch("nVtx",&nvtx,"nvtx/i");
  outTree_->Branch("nPUint",&nPUint,"nPUint/i");
  outTree_->Branch("nAK4PUPPIJets", &nAK4PUPPIJets,"nAK4PUPPIJets/i");
  outTree_->Branch("nAK4CHSJets", &nAK4CHSJets,"nAK4CHSJets/i");
  //outTree_->Branch("nAK4GenJets", &nAK4GenJets,"nAK4GenJets/i");
  outTree_->Branch("nPFCands", &nPFCand,"nPFCand/i");
  outTree_->Branch("PFCandPt", &PFCandPt);
  // outTree_->Branch("PFCandPx", &PFCandPx);
  //outTree_->Branch("PFCandPy", &PFCandPy);
  //outTree_->Branch("PFCandPz", &PFCandPz);
  outTree_->Branch("PFCandEta", &PFCandEta);
  outTree_->Branch("PFCandAbsEta", &PFCandAbsEta);
  outTree_->Branch("PFCandPhi", &PFCandPhi);
  outTree_->Branch("PFCandE", &PFCandE);
  outTree_->Branch("PFCandpdgId", &PFCandpdgId);
  outTree_->Branch("PFCandCharge", &PFCandCharge);
  outTree_->Branch("PFCandPUPPIw", &PFCandPUPPIw);
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
  outTree_->Branch("AK4PUPPIJetPt", &AK4PUPPIJetPt);
  outTree_->Branch("AK4PUPPIJetEta", &AK4PUPPIJetEta);
  outTree_->Branch("AK4PUPPIJetPhi", &AK4PUPPIJetPhi);
  outTree_->Branch("AK4PUPPIJetE", &AK4PUPPIJetE);
  outTree_->Branch("AK4PUPPIJetRawPt", &AK4PUPPIJetRawPt);
  outTree_->Branch("AK4PUPPIJetRawE", &AK4PUPPIJetRawE);
  outTree_->Branch("AK4CHSJetPt", &AK4CHSJetPt);
  outTree_->Branch("AK4CHSJetEta", &AK4CHSJetEta);
  outTree_->Branch("AK4CHSJetPhi", &AK4CHSJetPhi);
  outTree_->Branch("AK4CHSJetE", &AK4CHSJetE);
  outTree_->Branch("AK4CHSJetRawPt", &AK4CHSJetRawPt);
  outTree_->Branch("AK4CHSJetRawE", &AK4CHSJetRawE);
  /*outTree_->Branch("AK4GenJetPt", &AK4GenJetPt);
  outTree_->Branch("AK4GenJetEta", &AK4GenJetEta);
  outTree_->Branch("AK4GenJetPhi", &AK4GenJetPhi);*/
  
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
  //run    = iEvent.getRun();
  evt    = iEvent.id().event();
  lumi   = iEvent.id().luminosityBlock();
  //lumi   = iEvent.getLuminosityBlock();

  Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()) {
    std::cout << "vertices is empty!!! skipping event..." << endl;
    return;
  }
  nvtx = vertices->size();

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
 
  
  Handle<pat::JetCollection> AK4PUPPIJets;
  iEvent.getByToken(PUPPIjetToken_, AK4PUPPIJets);
  /*if(AK4PUPPIJets->empty()) {
    std::cout << "AK4PUPPIJets is empty!!! skipping event..." << endl;
    return;
    }*/

  //nAK4PUPPIJets = AK4PUPPIJets->size(); count only jets with pt > 20 and |eta| < 2.4
  
  Handle<pat::JetCollection> AK4CHSJets;
  iEvent.getByToken(CHSjetToken_, AK4CHSJets);
  /*if(AK4CHSJets->empty()) {
    std::cout << "AK4CHSJets is empty!!! skipping event..." << endl;
    return;
  }
  //nAK4CHSJets = AK4CHSJets->size(); count only jets with pt > 20 and |eta| < 2.4

  Handle<reco::GenJetCollection> AK4GenJets;
  iEvent.getByToken(GenjetToken_, AK4GenJets);
  if(AK4GenJets->empty()) {
    std::cout << "AK4GenJets is empty!!! skipping event..." << endl;
    return;
  }
  //nAK4GenJets = AK4GenJets->size(); count only jets with pt > 20 and |eta| < 2.4
  
  */

  Handle<pat::PackedCandidateCollection> PFCands;
  iEvent.getByToken(pfcandToken_, PFCands);
  if(PFCands->empty()) {
    std::cout << "PFCands is empty!!! skipping event..." << endl;
    return;
  }
  nPFCand = PFCands->size();

  //loop on PF candidates to store their Pt, Eta, Phi, E, pdgId, charge, PUPPI weight
  for (int i = 0; i < nPFCand; i++) {

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

  }

  

  //loop on PUPPI jets to store their Pt, Eta, Phi
  for (long unsigned int i = 0; i < AK4PUPPIJets->size(); i++) {

    if ( AK4PUPPIJets->at(i).pt() > 20 && fabs(AK4PUPPIJets->at(i).eta()) < 2.4 ) {
      
      nAK4PUPPIJets++;
      AK4PUPPIJetPt.push_back(AK4PUPPIJets->at(i).pt());
      AK4PUPPIJetEta.push_back(AK4PUPPIJets->at(i).eta());
      AK4PUPPIJetPhi.push_back(AK4PUPPIJets->at(i).phi());
      AK4PUPPIJetE.push_back(AK4PUPPIJets->at(i).energy());
      AK4PUPPIJetRawPt.push_back(AK4PUPPIJets->at(i).correctedJet("Uncorrected").pt());
      AK4PUPPIJetRawE.push_back(AK4PUPPIJets->at(i).correctedJet("Uncorrected").energy());
      
    }

  }
   
  //loop on CHS jets to store their Pt, Eta, Phi
  for (long unsigned int i = 0; i < AK4CHSJets->size(); i++) {

    if ( AK4CHSJets->at(i).pt() > 20 && fabs(AK4CHSJets->at(i).eta()) < 2.4 ) {

      nAK4CHSJets++;
      AK4CHSJetPt.push_back(AK4CHSJets->at(i).pt());
      AK4CHSJetEta.push_back(AK4CHSJets->at(i).eta());
      AK4CHSJetPhi.push_back(AK4CHSJets->at(i).phi());
      AK4CHSJetE.push_back(AK4CHSJets->at(i).energy());
      AK4CHSJetRawPt.push_back(AK4CHSJets->at(i).correctedJet("Uncorrected").pt());
      AK4CHSJetRawE.push_back(AK4CHSJets->at(i).correctedJet("Uncorrected").energy());

    }
   
  }

  /*
  //loop on Gen jets to store their Pt, Eta, Phi
  for (long unsigned int i = 0; i < AK4GenJets->size(); i++) {

    if ( AK4GenJets->at(i).pt() > 20 && fabs(AK4GenJets->at(i).eta()) < 2.4 ) {

      nAK4GenJets++;
      AK4GenJetPt.push_back(AK4GenJets->at(i).pt());
      AK4GenJetEta.push_back(AK4GenJets->at(i).eta());
      AK4GenJetPhi.push_back(AK4GenJets->at(i).phi());

    }
   
  }
  
  */

  outTree_->Fill();
 
  nAK4PUPPIJets = 0;
  nAK4CHSJets = 0;
  nAK4GenJets = 0;
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
  AK4PUPPIJetPt.clear();
  AK4PUPPIJetEta.clear();
  AK4PUPPIJetPhi.clear();
  AK4PUPPIJetE.clear();
  AK4PUPPIJetRawPt.clear();
  AK4PUPPIJetRawE.clear();
  AK4CHSJetPt.clear();
  AK4CHSJetEta.clear();
  AK4CHSJetPhi.clear();
  AK4CHSJetE.clear();
  AK4CHSJetRawPt.clear();
  AK4CHSJetRawE.clear();
  AK4GenJetPt.clear();
  AK4GenJetEta.clear();
  AK4GenJetPhi.clear();

  












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
