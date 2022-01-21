# Installing the Package:

Follow these instructions to checkout CMSSW_10_2_0 branch   
```
mkdir myFolder   
cmsrel CMSSW_10_2_0   
cd CMSSW_10_2_0/src   
git clone https://github.com/rrabadan/nanoAOD-tools.git PhysicsTools/NanoAODTools   
git clone https://github.com/tahuang1991/HhhAnalysis.git   
cd HhhAnalysis
git checkout origin/CMSSW_10_2_0
scram b -j 12   
```
tahuang1991/nanoAOD-tools is built based on cms-nanoAOD/nanoAOD-tools with minor change
