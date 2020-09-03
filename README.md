# PFCandInfoAnalyzer
This analyzer is used to store information about PF candidates and jets in flatTree format.

## Usage
The analyzer will work with CMSSW_10_6_>=17. To get started do:

``bash
cmsrel CMSSW_10_6_17
cd CMSSW_10_6_17/src/
cmsenv
mkdir PFCandInfo
cd PFCandInfo
git clone https://github.com/fiemmi/PFCandInfoAnalyzer
scramv1 b -j 5

