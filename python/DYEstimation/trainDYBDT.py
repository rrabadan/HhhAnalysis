#! /usr/bin/env python
import ROOT

import sys, os
import copy

#sys.path.append("../../CommonTools/mvaTraining/")
from tmvaTools import *

#sys.path.append('/nfs/soft/python/python-2.7.5-sl6_amd64_gcc44/lib/python2.7/site-packages/storm-0.20-py2.7-linux-x86_64.egg')
#sys.path.append('/nfs/soft/python/python-2.7.5-sl6_amd64_gcc44/lib/python2.7/site-packages/MySQL_python-1.2.3-py2.7-linux-x86_64.egg')
#sys.path.append('/home/taohuang/DiHiggsAnalysis/CMSSW_9_4_0_pre1/src/HhhAnalysis/python/NanoAOD')
#
#CMSSW_BASE = os.environ['CMSSW_BASE']
#SCRAM_ARCH = os.environ['SCRAM_ARCH']
#sys.path.append(os.path.join(CMSSW_BASE,'bin', SCRAM_ARCH))

#import Samplelist
def get_sample(inFileDir, samplename):
    sampleinfo = {}
    sampleinfo["files"] = [ os.path.join(inFileDir, samplename+"_Friend.root") ]
    #allsamplenames = [x.split('/')[1] for x in Samplelist.datasets]
    #index = allsamplenames.index(samplename)
    #sampleinfo["cross_section"] = Samplelist.MCxsections[index]
    tfile = ROOT.TFile(sampleinfo["files"][0] ,"READ")
    h_cutflow = tfile.Get("h_cutflow")
    sampleinfo["cross_section"] = tfile.Get("cross_section").GetVal()
    sampleinfo["event_weight_sum"] = h_cutflow.GetBinContent(1)
    sampleinfo["relativeWeight"] = sampleinfo["cross_section"]/h_cutflow.GetBinContent(1)
    return sampleinfo

date = "2017_02_17"
suffix = "bb_cc_vs_rest_10var"
label_template = "DATE_BDTDY_SUFFIX"

inFileDir = "/fdata/hepx/store/user/taohuang/HHNtuple_20180405_DYestimation/"

# SAMPLES FOR THE TRAINING

#DYJetsToLL_M10to50 = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_extended_ext0_plus_ext1_v4.4.2+80X_HHAnalysis_2017-02-13.v2"
#DYJetsToLL_M10to50_db = get_sample(unicode(DYJetsToLL_M10to50))
#
#DYJetsToLL_M50_0J = "DYToLL_0J_13TeV-amcatnloFXFX-pythia8_extended_ext0_plus_ext1_v4.4.2+80X_HHAnalysis_2017-02-13.v2"
#DYJetsToLL_M50_0J_db = get_sample(unicode(DYJetsToLL_M50_0J))
#DYJetsToLL_M50_1J = "DYToLL_1J_13TeV-amcatnloFXFX-pythia8_extended_ext0_plus_ext1_v4.4.2+80X_HHAnalysis_2017-02-13.v2"
#DYJetsToLL_M50_1J_db = get_sample(unicode(DYJetsToLL_M50_1J))
#DYJetsToLL_M50_2J = "DYToLL_2J_13TeV-amcatnloFXFX-pythia8_extended_ext0_plus_ext1_v4.4.2+80X_HHAnalysis_2017-02-13.v2"
#DYJetsToLL_M50_2J_db = get_sample(unicode(DYJetsToLL_M50_2J))
#
##bkgFiles = { 
#        "DYJetsToLL_M-10to50": {
#                    "files": [ inFileDir + DYJetsToLL_M10to50 + "_histos.root" ],
#                    "relativeWeight": DYJetsToLL_M10to50_db.source_dataset.xsection / DYJetsToLL_M10to50_db.event_weight_sum
#                },
#        "DYJetsToLL_M-50_0J": {
#                    "files": [ inFileDir + DYJetsToLL_M50_0J + "_histos.root" ],
#                    "relativeWeight": DYJetsToLL_M50_0J_db.source_dataset.xsection / DYJetsToLL_M50_0J_db.event_weight_sum
#                },
#        "DYJetsToLL_M-50_1J": {
#                    "files": [ inFileDir + DYJetsToLL_M50_1J + "_histos.root" ],
#                    "relativeWeight": DYJetsToLL_M50_1J_db.source_dataset.xsection / DYJetsToLL_M50_1J_db.event_weight_sum
#                },
#        "DYJetsToLL_M-50_2J": {
#                    "files": [ inFileDir + DYJetsToLL_M50_2J + "_histos.root" ],
#                    "relativeWeight": DYJetsToLL_M50_2J_db.source_dataset.xsection / DYJetsToLL_M50_2J_db.event_weight_sum
#                },
#        }
#
 
DYJetsToLL_M10to50_db = get_sample(inFileDir, "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")
DYJetsToLL_M50_0J_db = get_sample(inFileDir, "DYToLL_0J_13TeV-amcatnloFXFX-pythia8")
DYJetsToLL_M50_1J_db = get_sample(inFileDir, "DYToLL_1J_13TeV-amcatnloFXFX-pythia8")
DYJetsToLL_M50_2J_db = get_sample(inFileDir, "DYToLL_2J_13TeV-amcatnloFXFX-pythia8")

bkgFiles = { 
        "DYJetsToLL_M-10to50": { 
	           "files": DYJetsToLL_M10to50_db["files"],
		   "relativeWeight": DYJetsToLL_M10to50_db["relativeWeight"]
                },
        "DYJetsToLL_M-50_0J": {
	           "files": DYJetsToLL_M50_0J_db["files"],
		   "relativeWeight": DYJetsToLL_M50_0J_db["relativeWeight"]
                },
        "DYJetsToLL_M-50_1J": {
	           "files": DYJetsToLL_M50_1J_db["files"],
		   "relativeWeight": DYJetsToLL_M50_1J_db["relativeWeight"]
                },
#        "DYJetsToLL_M-50_2J": {
#	           "files": DYJetsToLL_M50_2J_db["files"],
#		   "relativeWeight": DYJetsToLL_M50_2J_db["relativeWeight"]
#                },
        }

print "bgfiles  ",bkgFiles

discriList = [
        "jet1_pt",
        "jet1_eta",
        "jet2_pt",
        "jet2_eta",
        "jj_pt",
        "ll_pt",
        "ll_eta",
        "llmetjj_DPhi_ll_met",
        "ht",
        "nJetsL"
        ]
#discriList = [
#        "jet1_pt",
#        "jet1_eta",
#        "jet2_pt",
#        "jet2_eta",
#        "jj_DPhi_j_j",
#        "ll_pt",
#        "ll_eta",
#        "llmetjj_DPhi_ll_jj",
#        "llmetjj_DPhi_ll_met",
#        "lljj_pt"
#    ]

spectatorList = []
cut = "(isMuMu || isElEl)"
MVAmethods = ["kBDT"]
#weightExpr = "event_weight * trigeff * llidiso * pu * sample_weight"
weightExpr = "event_reco_weight * sample_weight"

sigFiles = copy.deepcopy(bkgFiles)

#sigSelection = "(gen_bb || gen_cc)"
#bkgSelection = "(!(gen_bb || gen_cc))"
genbb_selection = "genjet1_partonFlavour == 5 && genjet2_partonFlavour == 5"
gencc_selection = "genjet1_partonFlavour == 4 && genjet2_partonFlavour == 4"
sigSelection = "((%s) || (%s))"%(genbb_selection, gencc_selection)
bkgSelection = "(! %s)"%(sigSelection)

label = label_template.replace("DATE", date).replace("SUFFIX", suffix)

if __name__ == "__main__":
    trainMVA(bkgFiles, sigFiles, discriList, cut, weightExpr, MVAmethods, spectatorList, label, sigWeightExpr=sigSelection, bkgWeightExpr=bkgSelection, nSignal=1000000, nBkg=1000000)
