import os, sys, random
import ROOT
from ROOT import TFile,TChain,TH1F,TH2F,TLegend
import Samplelist 


outroot = ROOT.TFile("allweight.root","RECREATE")
for key in Samplelist.Ntuplelist_Louvain:
    print "key ",key, " sample ",Samplelist.Ntuplelist_Louvain[key]
    
