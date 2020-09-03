# PFCandInfoAnalyzer
This analyzer is used to store information about PF candidates and jets in flatTree format.

## Getting started
The analyzer will work with CMSSW_10_6_>=17. To get started do:

```shell
cmsrel CMSSW_10_6_17
cd CMSSW_10_6_17/src/
cmsenv
mkdir PFCandInfo
cd PFCandInfo
git clone https://github.com/fiemmi/PFCandInfoAnalyzer
scramv1 b -j 5
```
## Working on MC or data
As the MC true information is not available in data, to avoid ProductNotFound errors, you should tell the Config file if you run on MC or data. Change the variable isMC accordingly in the Config file.

## Trigger workflow
To know if an event passed a given trigger, the branch ```std::vector<bool>``` triggerBit is used. First, you should look at the TriggerNames histogram to know which triggers were saved in the ntuples. For example, you can see here:

![alt text](http://fiemmi.web.cern.ch/fiemmi/JetMET/TriggerNames.png)

that 5 triggers were saved in the ntuples. Then, if you want to see if an event passed the trigger which appears in the i-th bin of TriggerNames, look at the content of ```triggerBit->at(i-1)```, e.g., if you want to use HLT_IsoMu27_v, you can simply do 

```
if (triggerBit->at(3)) {

...do things

}
```
### MC
The Config file is written in such a way that no trigger is applied at ntuplization level. Thus, to apply triggers, you shall simply do what is described in the previous lines.
### Data
Due to the big size of data events, the Config file is written in such a way that triggers are applied at ntuplization level. This is done _via_ ```EDFilters```. For example, to get a data sample to be used to select semileptonic ttbar events, you may want to do:

```python
process.triggerSelectionTTToSemileptonic = cms.EDFilter("TriggerResultsFilter",
                                       triggerConditions = cms.vstring(
                                       "HLT_IsoMu27_v*",
                                       "HLT_Ele35_WPTight_Gsf_v*",
                                       ),
                                       hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                       l1tResults = cms.InputTag( "" ),
                                       throw = cms.bool(False)
                                       )
```
which results in the logical OR of HLT_IsoMu27_v and HLT_Ele35_WPTight_Gsf_v, and ad it to your path:

```python
process.p = cms.Path(process.triggerSelectionTTToSemileptonic*process.puppiSequence*process.GetPFInfo)
```

# PUPPI tune and b tagging
PUPPI tune is updated to v15 in the Config file. The b tagging information is not saved for PUPPIv15, so that b tagging information is stored for CHS jets. In particular, jets are tagged with DeepCSV using the medium WP.
