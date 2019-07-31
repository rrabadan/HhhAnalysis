#!/usr/bin/python
import os
import sys
sys.path.append('/home/taohuang/DiHiggsAnalysis/CMSSW_9_4_0_pre1/src/HhhAnalysis/python/NanoAOD/')
import ROOT
from localSamplelist import *

def getTotalEvents(filename, treename= "Friends"):
    h1 = ROOT.TH1F("h1","h1",10,0,10)
    ch = ROOT.TChain(treename)
    ch.Add(filename)
    ch.Draw("1>>h1","ll_M<76")
    return h1.GetEntries()

def getEventWeightSum(filename):
    tfile = ROOT.TFile(filename, "READ")
    h_cutflow = tfile.Get("h_cutflow")
    return h_cutflow.GetBinContent(1)

def getCrossSection(filename):
    tfile = ROOT.TFile(filename, "READ")
    xsec = tfile.Get("cross_section")
    return xsec.GetVal()


benchmarks =  ["ttV","Wjet","sT","VV","DY","TT","Data"]
#benchmarks =  ["ttV","Wjet","sT","VV","DY","Data"]
#benchmarks = ["ttV","Wjet","sT","VV","DY","TT","Data","radion","graviton"]
#benchmarks = ["radion_M400","graviton_M400"]
signals = []
masses = [260, 270, 300, 350, 400, 450, 500, 550, 600, 650,750, 800, 900]
#masses = [750]
for mass in masses:
	signals.append("RadionM%d"%mass)
	#benchmarks.append("graviton_M%d"%mass)
#alliters = [100, 1000, 10000, 100000, 1000000]
benchmarks = signals
#benchmarks = ['TT']
#benchmarks = signals
#benchmarks.append("sT_top")
#benchmarks.append("sT_antitop")
#benchmarks.append("DYM10to50")#69k
#benchmarks.append("DYToLL0J")#838k
#benchmarks.append("DYToLL1J")#5095k
#benchmarks.append("DYToLL2J")#15265k
#

#jobdir = "HMEJobs_Louvain_Test"
#jobdir = "TTJobs_Louvain"
#jobdir = "STJobs_Louvain"
#jobdir = "DYJobs_Louvain"
#jobdir = "20170814_DY_Louvain"
write_xsec_sumweight = True
jobdir = "20180412_HHbbWW_addHME_10k_final"
outputdir = "/fdata/hepx/store/user/taohuang/%s/"%(jobdir)
for job in benchmarks:
  for dataname in full_local_samplelist[job].keys():
	filename = full_local_samplelist[job][dataname]['path']
	totalevents = getTotalEvents(filename, "Friends")

	finalfile = outputdir + dataname + '_final.root'
	#if 'Radion' in job or job == 'TT':
	#    os.system("hadd -f %s"%finalfile +' '+ outputdir + job + '_ijob*.root')
	#else:
	#os.system("hadd -f %s"%finalfile +' '+ outputdir + dataname + '_ijob*.root')
	totalevents2 = getTotalEvents(finalfile, "Friends")
	if totalevents2 != totalevents :
	    print "the event number changed after adding HME: ",totalevents," ---> ",totalevents2
	if write_xsec_sumweight and job != "Data":
	    print "file  ",filename
	    event_weight_sum = getEventWeightSum(filename)
	    tfile = ROOT.TFile(finalfile,"UPDATE")
	    if 'Radion' not in job :
		cross_section = getCrossSection(filename)
		p1 = ROOT.TParameter(float)("cross_section", xsec)
		p1.Write()
	    p2 = ROOT.TParameter(float)("event_weight_sum", event_weight_sum)
	    p2.Write()
            tfile.Close()

	  
