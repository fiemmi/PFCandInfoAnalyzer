#!/bin/bash
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node
source /cvmfs/cms.cern.ch/cmsset_default.sh  ## if a bash script, use .sh instead of .csh
cd /eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/PFCandInfo/PFCandInfoAnalyzer/createTreePU_noPU_framework/merged_files/EXT80k_v9-v1
eval `scramv1 runtime -sh` # cmsenv is an alias not on the workers
declare -i nfile=${1}
root -l -q -b '/eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/createTreePU_noPU_framework/createTreePU_noPU.cpp('$nfile')'

