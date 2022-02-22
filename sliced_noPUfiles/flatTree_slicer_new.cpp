#include <iostream>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

using namespace std;

void flatTree_slicer_new(int n_evts_max, int n_evts_per_file) {

  gROOT->Reset();
  int n_output_files = n_evts_max/n_evts_per_file;
  TFile *outputfile[n_output_files];
  //Get old file, old tree and set top branch address
  TFile *oldfile = new TFile("/eos/user/f/fiemmi/JetMET/storage/PU_PFCandInfoAnalizer_files/FlatPU0to75_RunIISummer20UL17_MiniAODv2_EXT350k/flatTree_QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8_FlatPU0to75_RunIISummer20UL17_MiniAODv2_EXT350k.root");
  TTree *oldtree = (TTree*)oldfile->Get("events");
  TTree *newtree[n_output_files];
  
  for (int i = 0; i < n_output_files; i++){
    newtree[i] = oldtree->CloneTree(0);
    std::string str_file_counter = to_string(i);
    TString Str_file_counter = str_file_counter;
    outputfile[i] = new TFile("/eos/user/f/fiemmi/JetMET/storage/PU_sliced_files/FlatPU0to75_RunIISummer20UL17_MiniAODv2_EXT350k/flatTree_QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8_FlatPU0to75_RunIISummer20UL17_MiniAODv2_EXT350k_part"+Str_file_counter+".root","recreate");
    
    if (i > 0) {

      newtree[i-1]->AutoSave();
      delete newtree[i-1];
      outputfile[i-1]->Close();
      delete outputfile[i-1];
	
    }

    newtree[i]=oldtree->CopyTree("","",n_evts_per_file,n_evts_per_file*i);
    outputfile[i]->Write();
    cout << "Processed slice " << i << endl;
   
  }

}


