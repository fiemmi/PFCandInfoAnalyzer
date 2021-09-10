#!/bin/bash
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node
source /cvmfs/cms.cern.ch/cmsset_default.sh  ## if a bash script, use .sh instead of .csh
cd /afs/cern.ch/work/f/fiemmi/private/CMSSW_10_6_20/src/PFCandInfo/PFCandInfoAnalyzer/condor/PFCandInfoAnalyzer
eval `scramv1 runtime -sh` # cmsenv is an alias not on the workers
line=$((${1}+1))
inFile=$(awk "NR == ${line}" ${2})
echo $inFile
cmsRun /afs/cern.ch/work/f/fiemmi/private/CMSSW_10_6_20/src/PFCandInfo/PFCandInfoAnalyzer/python/ConfFile_cfg.py inputFiles=${inFile} outputFile=condorJob_${1}.root
cp condorJob_${1}.root /afs/cern.ch/work/f/fiemmi/private/CMSSW_10_6_20/src/PFCandInfo/PFCandInfoAnalyzer/condor/PFCandInfoAnalyzer/files/EXT80k_v9-v1
rm condorJob_${1}.root
