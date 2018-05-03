import sys
sys.path.append("/home/taohuang/DiHiggsAnalysis/CMSSW_9_4_0_pre1/src/HhhAnalysis/python/NanoAOD")
import Samplelist as Slist
import os

full_local_samplelist = {}
localdir = "HHNtuple_20180412/"
#localdir = "/fdata/hepx/store/user/taohuang/HHNtuple_20180404/"
localdir = "/fdata/hepx/store/user/taohuang/HHNtuple_20180502_dataTT/"
for i,isample in enumerate(Slist.NumSample):
    
    if int(isample) <0:
        continue
    sampleName = Slist.sampleN_short[i]
    dataname =  Slist.Nanodatasets[i].split('/')[1]
    localfilepath = os.path.join(localdir, dataname+"_Friend.root")
    xsec =  Slist.MCxsections[i]
    if sampleName not in full_local_samplelist.keys():
        full_local_samplelist[sampleName] = {}

    full_local_samplelist[sampleName][dataname] = {}
    full_local_samplelist[sampleName][dataname]["path"] = localfilepath
    full_local_samplelist[sampleName][dataname]["cross_section"] = xsec
    #if sampleName == 'TT':
    #   full_local_samplelist[sampleName][dataname]["path"] = '/fdata/hepx/store/user/taohuang/HHNtuple_20180501_dataTT_v2/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8_Friend.root'
full_local_samplelist["Data"] = {}
datanames = ["DoubleMuon", "DoubleEG","MuonEG"]
for dataname in datanames:
    localdatadir = "/fdata/hepx/store/user/taohuang/HHNtuple_20180502_dataonly_HLT_v2/"
    #localdatadir = "/fdata/hepx/store/user/taohuang/HHNtuple_20180412/"
    full_local_samplelist["Data"][dataname] = {}
    full_local_samplelist["Data"][dataname]["path"] =  os.path.join(localdatadir,dataname+"Run2016.root")
#print full_local_samplelist
