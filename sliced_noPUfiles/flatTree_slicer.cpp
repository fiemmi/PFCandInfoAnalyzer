#include <iostream>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

using namespace std;

void flatTree_slicer(int n_evts_max, int n_evts_per_file) {

  gROOT->Reset();
  TFile *outputfile[n_evts_max/n_evts_per_file];
  //Get old file, old tree and set top branch address
  TFile *oldfile = new TFile("/eos/user/f/fiemmi/JetMET/storage/noPU_crab_files/flatTree_QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8_training_EpsilonPU_EXT1kk_withPUPPIalpha.root");
  //TFile *oldfile = new TFile("test_file_part1.root");
  TTree *oldtree = (TTree*)oldfile->Get("events");
  TTree *newtree[n_evts_max/n_evts_per_file];
  Int_t nentries = (Int_t)oldtree->GetEntries();
 
  //int divide = n; //number of output files
  int file_counter = 0;

  for (int i = 0; i < n_evts_max+1; i++){

    if (i%1000 == 0) cout << "ievent = " << i << endl;
    
    //Create a new file + a clone of old tree in new file
    if (i%n_evts_per_file == 0) {

      if (i > 0) {

	file_counter++;
	newtree[file_counter-1]->AutoSave();
	delete newtree[file_counter-1];
	outputfile[file_counter-1]->Close();
	delete outputfile[file_counter-1];
	cout << "Processed slice " << file_counter << endl;
	if (i == n_evts_max) break;

      }

      std::string str_file_counter = to_string(file_counter);
      TString Str_file_counter = str_file_counter;
      outputfile[file_counter] = new TFile("/eos/user/f/fiemmi/JetMET/storage/noPU_sliced_files/flatTree_QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8_training_EpsilonPU_EXT1kk_withPUPPIalpha_"+Str_file_counter+".root","recreate");
      newtree[file_counter] = oldtree->CloneTree(0);
  
    }

    oldtree->GetEntry(i);
    newtree[file_counter]->Fill();

  }

}


