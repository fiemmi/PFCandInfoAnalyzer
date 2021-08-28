#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TBenchmark.h>
#include <TTreeIndex.h>
#include <fstream>
#include <string>
#include <iostream>

using namespace std;

void sortTrees (int nfile) {
  //get the block of indices of no-PU events that you want to write in the correct order. They lay in a .txt file
  ifstream inputfile (Form("/eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/sorting_framework/lists_of_indices/EXT80k_v9-v1/list_%i.txt", nfile));
  cout << Form("/eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/sorting_framework/lists_of_indices/EXT80k_v9-v1/list_%i.txt", nfile) << endl;
  string line;
  std::vector<int> indices;
  if (inputfile.is_open()) {
    while ( getline (inputfile,line) ) {
      indices.push_back(stoi(line)); //store the indices in a std::vector
    }
    inputfile.close();
    //for(int i = 0; i < indices.size(); i++) cout << indices.at(i) << endl;

    TFile* fIn = new TFile ("/eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/flatTree_QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8_training_PU_EXT80k_withPUPPIalpha_v9-v1.root");
    TTree* tIn = (TTree*) fIn->Get("events");
    
    TFile* fOut = new TFile(Form("flatTree_QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8_training_PU_EXT80k_withPUPPIalpha_v9-v1_%i.root", nfile), "recreate");
    TTree* tOut = tIn->CloneTree(0);
    
    for (int i = 0; i < indices.size(); i++) {
      //if (i == 10) break; //make it shorter; debug
      if (i%1000 == 0) cout << "ievent = " << i << endl;
      tIn->GetEntry(indices.at(i));
      tOut->Fill();
    }
    fOut->cd();
    fOut->Write();
    
  }

  else cout << "Unable to open file";

}
