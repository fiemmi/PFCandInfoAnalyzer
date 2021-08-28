#!/bin/bash
#set -x
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node
source /cvmfs/cms.cern.ch/cmsset_default.sh  ## if a bash script, use .sh instead of .csh
cd /eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/sorting_framework/sorted_files
eval `scramv1 runtime -sh` # cmsenv is an alias not on the workers
declare -i nfile=${1}
root -l -q -b '/eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/PFCandInfo/PFCandInfoAnalyzer/sorting_framework/sortTrees.cpp('$nfile')'
