#include <iostream>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

using namespace std;

void flatTree_slicer(int n) {

gROOT->Reset();

//Get old file, old tree and set top branch address
TFile *oldfile = new TFile("../flatTree_QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8_training_EpsilonPU_EXT80k_withPUPPIalpha_v9-v1.root");
//TFile *oldfile = new TFile("test_file_part1.root");
TTree *oldtree = (TTree*)oldfile->Get("events");
Int_t nentries = (Int_t)oldtree->GetEntries();

int divide = n; //number of output files

for (int nfile = 0; nfile < divide; nfile++){

//Create a new file + a clone of old tree in new file
std::string str_file_counter = to_string(nfile);
TString Str_file_counter = str_file_counter;
TFile *newfile = new TFile("files/EXT80k_v9-v1/flatTree_QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8_training_EpsilonPU_EXT80k_withPUPPIalpha_v9-v1_"+Str_file_counter+".root","recreate");
TTree *newtree = oldtree->CloneTree(0);

 for (Int_t i=(nentries/divide)*nfile; i<(nentries/divide)*(1+nfile); i++) {
   oldtree->GetEntry(i);
   newtree->Fill();
 }

newtree->AutoSave();
delete newfile;
cout << "Processed slice " << nfile +1 << endl;
}

delete oldfile;

}


