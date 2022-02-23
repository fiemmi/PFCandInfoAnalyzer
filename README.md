# PFCandInfoAnalyzer
This analyzer is used to store information about PF candidates and jets in a flatTree format and preprocess it for ML tasks.

## Getting started
The analyzer will work with CMSSW_10_6_>=17. It's always good practice to check on the PdmV TWikis which CMSSW release is recommended for a given input miniAOD file. To get started do:

```shell
cmsrel CMSSW_10_6_17
cd CMSSW_10_6_17/src/
cmsenv
mkdir PFCandInfo
cd PFCandInfo
git clone https://github.com/fiemmi/PFCandInfoAnalyzer
scramv1 b -j 5
```
## General idea of this repository
Scripts contained in this repository were originally written to accomplish pileup (PU) mitigation studies inside the CMS Collaboration. The main idea is to start from two Monte Carlo QCD files containing the same simulated events, the first file being a no-PU file, while the second one having PU superimposed. The final goal is to have a `.root` output file storing, event-by-event, information about ParticleFlow (PF) candidates coming from PU and n-PU events in a coherent way.

## Step-by-step procedure
The procedure to get the wanted output `.root` file is made by several steps. They are listed below together with an explanation for each of them.
**CAVEAT**: the following steps won't work out-of-the-box. You'll need to make changes here and there: update the paths to macros, the number of jobs to be run and so on to your needs.

1. **Ntuplize a given number of events from no-PU files**. In the configuration file of PFCandInfoAnalyzer (`python/ConfFile_cfg.py`), ntuplize a given number of events from a no-PU file through the `fileNames` parameter. The number can be chosen by editing the `process.maxEvents` parameter.
```
cmsenv
voms-proxy-init --voms cms
cmsRun /your/path/PFCandInfo/PFCandInfoAnalyzer/python/ConfFile_cfg.py
```
2. **Save runNo, lumiSec, evtNo of no-PU events**. Run the script `printEvents.C` to save the run number, lumi section and event number of each no-PU event that you ntuplized in step 1 to a `.txt` file. 
3. **Extract the very same events from the PU file**. Let's suppose you ntuplized the first *n* events from the no-PU file. Running `cmsRun` on the first *n* events of the PU file won't ntuplize the same events (due to how events are generated). Thus, the very same events that have been ntuplized in step 1 must first be extracted from the PU sample. This is done through the `edmPickEvents.py ` script. It needs the DAS string for the file you want to extract events from and the `.txt` file that has been produced in step 2 as inputs. For example:
```
edmPickEvents.py "/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/RunIISummer19UL17MiniAOD-FlatPU0to70_106X_mc2017_realistic_v6-v3/MINIAODSIM" events_EpsilonPU_EXT80k.txt --crab
```
The `--crab` option will produce a crab configuration file that will do the job. You will have to slightly modify it to fit your needs (e.g., change the site where to write the output files). Launch CRAB jobs:
```
source /cvmfs/cms.cern.ch/crab3/crab.sh
crab submit -c pickevents_crab.py
```
4. **Launch condor jobs to ntuplize the PU events**. Edit the configuration file of PFCandInfoAnalyzer to work for condor jobs: `fileNames` should be equal to `cms.untracked.vstring(options.inputFiles)` and `fileName` should be equal to `cms.string(options.outputFile)`. Then, go to `condor/PFCandInfoAnalyzer` and edit `fileList.txt` with the path to the files that CRAB has produced. Finally, run the analyzer on them by doing:
```
condor_submit submit.sub
```
5. **Merge the output files of the condor jobs**.
```
hadd some_name_for_outputfile.root condorJob_*
```
6. **Save runNo, lumiSec, evtNo of PU events**. Repeat instructions in step 2 to save the run number, lumi section and event number of each PU event to a .txt file.
7. **Set up the sorting of TTrees**. You have now two lists of run number, lumi section and event number in the form of two `.txt` files. Feed them to the `sorting_framework/get_indices.py` script to find the correspondance between events in PU and no-PU files. Mind that there are many reasons why CRAB could fail in extracting all the no-PU events that you want from the PU file. Thus, you should feed the script **first** with the list of PU events and **then** with the list of no-PU events. If a PU event is not find in the no-PU set, it is simply skipped. This way you will get the PU--->noPU correspondance with no unmatched events. The script will produce a third list containing, for each event in the PU file, the index pointing at that event no-in the PU file. The script will also split this final list in *n* sublists, where *n* can be decided by the user, and save them.
8. **Sort no-PU events**. This is done by running condor jobs for `sorting_framework/sortTrees.cpp`. This macro takes a sublist produced in step 7 as input and uses the indices to sort the no-PU events in the same order of the PU events:
```
cd condor/sortTrees
condor_submit submit.sub
```
9. **Split the PU file in subfiles**. Take the file you got in step 5 and split it into the same number of files *n* that you got in step 7. For this you can use `flatTree_slicer.cpp` (or its faster version `flatTree_slicer_new.cpp`):
```
cd sliced_noPUfiles
root -l -q flatTree_slicer.cpp+(n_evts_max, n_evts_per_file)
// or root -l -q flatTree_slicer_new.cpp+(n_evts_max, n_evts_per_file)
```
10. **Merge information from PU and no-PU events**. Now that you have the same number of PU and no-PU files and events are in the same order, merge the files into a single file by running condor jobs executing `createTreePU_noPU_framework/createTreePU_noPU.cpp`:
```
cd condor/createTreePU_noPU
condor_submit submit.sub
```
11. **Perform matching between PF candidates**. Finally, match PF candidates from no-PU and PU events by running condor jobs for `PFmatching/PFmatching.cpp`:
```
cd condor/PFmatching
condor_submit submit.sub
```
