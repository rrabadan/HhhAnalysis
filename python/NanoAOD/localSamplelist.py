import sys
import ROOT
#sys.path.append("/home/taohuang/DiHiggsAnalysis/CMSSW_9_4_0_pre1/src/HhhAnalysis/python/NanoAOD")
sys.path.append("/home/taohuang/HhhAnalysis/python/NanoAOD")
import Samplelist as Slist
import os

def get_event_weight_sum_file(filepath):
    #print "filepath ",filepath
    tfile = ROOT.TFile(filepath, "READ")
    hist = tfile.Get("h_cutflow")
    event_weight_sum = hist.GetBinContent(1)
    tfile.Close()
    return event_weight_sum

full_local_samplelist = {}
#localdir = "/fdata/hepx/store/user/taohuang/HHNtuple_20180502_dataTT/"
localdir = "/data/taohuang/HHNtuple_20180502_dataTT/"
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
    localdatadir = "/data/taohuang/HHNtuple_20180502_dataonly_HLT_v2/"
    #localdatadir = "/fdata/hepx/store/user/taohuang/HHNtuple_20180502_dataonly_HLT_v2/"
    full_local_samplelist["Data"][dataname] = {}
    full_local_samplelist["Data"][dataname]["path"] =  os.path.join(localdatadir,dataname+"Run2016.root")
#print full_local_samplelist

untagged_samplelist = {}
untagged_MCname = ["TT"]
untagged_localdir = '/data/taohuang/HHNtuple_20180418_DYestimation_withMbtagweight/'
for mcname in untagged_MCname:
    untagged_samplelist[mcname] = {}
    #for key in untagged_samplelist[mcname].keys():
    for key in full_local_samplelist[mcname].keys():
        untagged_samplelist[mcname][key] = {}
        #untagged_samplelist[mcname][key]['path'] = os.path.join(untagged_localdir, key)
        untagged_samplelist[mcname][key]['path'] = os.path.join(untagged_localdir, key+"_Friend_untagged.root")
        untagged_samplelist[mcname][key]['cross_section'] = full_local_samplelist[ mcname][key]['cross_section']
        untagged_samplelist[mcname][key]['event_weight_sum'] =  get_event_weight_sum_file( full_local_samplelist[ mcname][key]['path'])

untagged_samplelist["Data"] = {}
datas = ["DoubleMuon", "DoubleEG"]
for dataname in datas:
    untagged_localdatadir = untagged_localdir
    untagged_samplelist["Data"][dataname] = {}
    untagged_samplelist["Data"][dataname]["path"] =  os.path.join(untagged_localdatadir,dataname+"Run2016_untagged.root")
