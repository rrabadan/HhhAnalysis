#! /bin/env python

from __future__ import division

import argparse
import os
import sys
from array import array

import numpy as np

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.Reset()


parser = argparse.ArgumentParser(description='Compute event category fraction on a given sample')
#group = parser.add_mutually_exclusive_group(required=True)
#
#group.add_argument('-d', '--json', type=str, metavar='FILE', help='JSON file describing the input')
#
#parser.add_argument('dummy', type=str, help="Dummy argument for compatibility with condorTools")
parser.add_argument('-o', '--outputdir',dest="outputdir", type=str, default= "./", help='output dir')
parser.add_argument('-i', '--inputdir',dest= "inputdir", type=str, metavar='STR', help='input dir')
parser.add_argument('-n', '--name',dest = "name", type=str, metavar='STR', help='Sample name')
#
options = parser.parse_args()
#


def get_sample(inFileDir, samplename):
    sampleinfo = {}
    sampleinfo["files"] = [ os.path.join(inFileDir, samplename+"_Friend.root") ]
    tfile = ROOT.TFile(sampleinfo["files"][0] ,"READ")
    h_cutflow = tfile.Get("h_cutflow")
    sampleinfo["cross_section"] = tfile.Get("cross_section").GetVal()
    sampleinfo["event_weight_sum"] = h_cutflow.GetBinContent(1)
    sampleinfo["relativeWeight"] = sampleinfo["cross_section"]/h_cutflow.GetBinContent(1)
    return sampleinfo

#files = ["/fdata/hepx/store/user/taohuang/HH_NanoAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIIFall17_NanoAOD_Friend.root"]
#samplelist = ['DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8','DYToLL_0J_13TeV-amcatnloFXFX-pythia8','DYToLL_1J_13TeV-amcatnloFXFX-pythia8','DYToLL_2J_13TeV-amcatnloFXFX-pythia8']
#inFileDir = '/fdata/hepx/store/user/taohuang/HHNtuple_20180410_DYestimation/'
inFileDir = options.inputdir

output = "dy_flavour_{jobname}.root".format(jobname = options.name)
output = os.path.join(options.outputdir, output)

print "=============================================================="
print "inputfileDir ", inFileDir," MC ",options.name
print "=============================================================="
print "outputFile ", output
print "=============================================================="

treename = "Friends"
sampleinfo = get_sample(inFileDir, options.name)
chain = ROOT.TChain( treename )
cross_section = sampleinfo['cross_section']
event_wgt_sum = sampleinfo['event_weight_sum']
for f in sampleinfo["files"]:
    chain.Add(f)



### 161218, 7 var + (real) ht, bb && cc vs. rest
#binning = np.asarray([-0.37849545503700767, -0.252662010953846, -0.21699735110544804, -0.19021412069054136, -0.171396453995042, -0.15676306613552415, -0.14321969667988593, -0.13153891023683148, -0.12107244096732181, -0.11105252860752612, -0.10194490420798973, -0.09332335631432828, -0.08536327877224573, -0.07801048636166992, -0.07120394344115692, -0.0647847651486761, -0.05993005923356425, -0.05408124411388808, -0.04856654729070964, -0.0428596004319479, -0.037361529766591905, -0.03188497839059719, -0.026731425996404594, -0.0216053159528186, -0.016547999887674154, -0.011556509911378651, -0.006637505166127072, -0.001911041855971231, 0.003017405742272187, 0.0076999820179795425, 0.012357098641527476, 0.017143242487650923, 0.021887924437848254, 0.026634181970421666, 0.03165959325886574, 0.03665486332869056, 0.04151823848231689, 0.04651083716266224, 0.051682726965222574, 0.056882360974048966, 0.06180807445815208, 0.06784876065619697, 0.07411968400707587, 0.08150640870709099, 0.09007916053735501, 0.10047784040389708, 0.11391604713844385, 0.12992068221277847, 0.1513985908285126, 0.18121531860439136, 0.3786900843132766],
#        dtype='float')
#binning = np.linspace(-0.38, 0.38, 25)

### 161220, 10 var
# 30 bins
binning = np.asarray([-0.4325139551124535, -0.2146539640268055, -0.17684879232551598, -0.1522156780133781, -0.13344360493544538, -0.1177783085968212, -0.10431773748076387, -0.09240803627202236, -0.08144732988778663, -0.07139562851774808, -0.06195872754019471, -0.053149265226606804, -0.044689436819594426, -0.036486494035769285, -0.028370020384749492, -0.02052289170780913, -0.01265119174726717, -0.004810595256756055, 0.003258152851774066, 0.01125285685430063, 0.019322492143167114, 0.02785483333896287, 0.03659553016370119, 0.04591206104108278, 0.05601279709011762, 0.06690819726322504, 0.07861467402378061, 0.09302953795299788, 0.11151410228370977, 0.13829367256021688, 0.333748766143408],
        dtype='float')

n_bins = len(binning) - 1

flavors = ["b", "c", "l"]
partonidmap  = {"b":5, "c":4, "l":0}

efficiencies = {}

for flav1 in flavors:
    for flav2 in flavors:
        name = "%s%s_frac" % (flav1, flav2)

        frac = ROOT.TEfficiency(name, name, n_bins, binning)
        frac.SetStatisticOption(ROOT.TEfficiency.kBUniform)
        frac.SetUseWeightedEvents()
        # Set global weight for sample
        frac.SetWeight(cross_section / event_wgt_sum)

        key = (flav1, flav2)
        efficiencies[key] = frac

ROOT.TH1.SetDefaultSumw2(True)

chain.SetBranchStatus("*", 0)

chain.SetBranchStatus("is*", 1)
chain.SetBranchStatus("lep*", 1)
chain.SetBranchStatus("jet*", 1)
chain.SetBranchStatus("ll*", 1)
chain.SetBranchStatus("jj*", 1)
chain.SetBranchStatus("ht*", 1)
chain.SetBranchStatus("nJetsL*", 1)
chain.SetBranchStatus("DPhi*", 1)
chain.SetBranchStatus("*partonFlavour*", 1)
chain.SetBranchStatus("event_reco_weight", 1)
chain.SetBranchStatus("sample_weight", 1)

bdt_tmva_variables = [
        "jet1_pt",
        "jet1_eta",
        "jet2_pt",
        "jet2_eta",
        "jj_pt",
        "ll_pt",
        #"ll_eta",
        #"llmetjj_DPhi_ll_met",
        "ht",
        "nJetsL"

]
#bdt_label = "2016_12_18_BDTDY_bb_cc_vs_rest_7var_ht_nJets"
bdt_label = "2017_02_17_BDTDY_bb_cc_vs_rest_10var"
#FIXME
bdt_xml_file = "/home/taohuang/DiHiggsAnalysis/CMSSW_9_4_0_pre1/src/HhhAnalysis/python/DYEstimation/DYBDTTraining/weights/{}_kBDT.weights.xml".format(bdt_label)
print "bdt_xml_file ",bdt_xml_file

dict_tmva_variables = { var: array('f', [0]) for var in bdt_tmva_variables }
m_reader = ROOT.TMVA.Reader("Silent=1")
for var in bdt_tmva_variables:
    m_reader.AddVariable(var, dict_tmva_variables[var])
m_reader.BookMVA(bdt_label, bdt_xml_file)
entries = None
print("Loading chain...")
if not entries:
    entries = chain.GetEntries()
print("Done.")

print("Computing jet flavor fraction using %d events." % entries)

for i in range(0, entries):
    chain.GetEntry(i)

    if (i % 100 == 0):
        print("Event %d over %d" % (i + 1, entries))


    if not (chain.isElEl or chain.isMuMu):
        continue

    if (chain.ll_M <= 12):
        continue

    # Weight: take into account? Also lepton ID SF?
    weight = chain.sample_weight * chain.event_reco_weight

    def pass_flavor_cut(flav1, flav2):
	return  chain.genjet1_partonFlavour == partonidmap[flav1] and chain.genjet2_partonFlavour == partonidmap[flav2]


    def get_value(object, val):
	return getattr(object, val)

    for var in bdt_tmva_variables:
        # Special treatment for variables not retrieved from the base object
	dict_tmva_variables[var][0] = get_value(chain, var)

    bdt_value = m_reader.EvaluateMVA(bdt_label)

    for flav1 in flavors:
        for flav2 in flavors:
            key = (flav1, flav2)
            efficiencies[key].FillWeighted(pass_flavor_cut(flav1, flav2), weight, bdt_value)

print("Done")
output = ROOT.TFile.Open(output, "recreate")

for key, value in efficiencies.items():
    value.Write()

output.Close()
