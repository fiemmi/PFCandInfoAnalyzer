#!/bin/bash
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node
source /cvmfs/cms.cern.ch/cmsset_default.sh  ## if a bash script, use .sh instead of .csh
cd /afs/cern.ch/work/f/fiemmi/private/CMSSW_10_6_20/src/PFCandInfo/PFCandInfoAnalyzer/PFmatching/files/EXT80k_v9-v1
eval `scramv1 runtime -sh` # cmsenv is an alias not on the workers
declare -i nfile=${1}
root -l -q -b '/afs/cern.ch/work/f/fiemmi/private/CMSSW_10_6_20/src/PFCandInfo/PFCandInfoAnalyzer/PFmatching/PFmatching.cpp+('$nfile', 0.04, 0.1)'
