# Installing the Package:

Follow these instructions:   
```
mkdir myFolder   
cmsrel CMSSW_8_0_26   
cd CMSSW_8_0_26/src   
git clone https://github.com/lpernie/HhhAnalysis.git   
cd HhhAnalysis   
scram b -j 12   
```

# Running on data or simulation:   
Data:   
```
cd DataFormats/scripts/patifyData_Run2016/
python multicrab3.py
```
MC:
```
cd HhhAnalysis/CutFlowAnalyzer/test/DY_generation
vim README
```
# Notes on DiHiggsWWBBAnalyzer.cc
This analyzer is targeted to run on DiHiggs, TTbar, DY samples.     
Sampletype is used to mark which sample it is running on.   
If you are using data, set Sampletype=0 and analyzer will skip gen-level selection   
If you are using signal sample, e.g. B6, you can set Sampletype=0, or Sampletype=6 which will check gen-level information   
If you are using TTbar, you can set Sampletype=0, or Sampletype=13 which will check gen-level information   
If you are using DY sample, you can set Sampletype=0, or Sampletype=14 which will check gen-level information( not implemented yet)   
The event will be saved only if it pass preselection at reco levle but we also sometimes keep track of gen-level information   

# Physics objects from Muon POG:   
https://twiki.cern.ch/twiki/bin/viewauth/CMS/POGRecipesICHEP2016   
