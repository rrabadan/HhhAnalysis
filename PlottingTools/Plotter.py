import random
import ROOT
from ROOT import TFile,TChain,TH1F,TH2F,TLegend
import os
execfile("start.py")
execfile("functions.py")

#Creating folders and parameters
tree_name="DiHiggsWWBBAna/evtree"
os.system("mkdir -p Plots")
os.system("mkdir -p HADD")
Lumi=36.8#fb-1
#TT TChains
TT_ch = ROOT.TChain(tree_name)
os.system("find /fdata/hepx/store/user/lpernie/TTTo2L2Nu_13TeV-powheg/TTTo2L2Nu_powheg_miniAODv2_v0_ext1_01/170322_221256 | grep root | grep -v failed > HADD/TT0.txt")
with open('HADD/TT0.txt','r') as f:
  for line in f:
    if not line.isspace():
      TT_ch.Add(str(line[:-1]))
print "TT has", TT_ch.GetEntries(), "entries."
#DT TChains
DY_ch = ROOT.TChain(tree_name)
os.system("find /fdata/hepx/store/user/lpernie/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/170323_021221 | grep root | grep -v failed > HADD/DY0.txt")
os.system("find /fdata/hepx/store/user/lpernie/DYToLL_0J_13TeV-amcatnloFXFX-pythia8/crab_DYToLL_0J_13TeV-amcatnloFXFX-pythia8/170323_021233 | grep root | grep -v failed >> HADD/DY0.txt")
os.system("find /fdata/hepx/store/user/lpernie/DYToLL_1J_13TeV-amcatnloFXFX-pythia8/crab_DYToLL_1J_13TeV-amcatnloFXFX-pythia8/170323_021246 | grep root | grep -v failed >> HADD/DY0.txt")
os.system("find /fdata/hepx/store/user/lpernie/DYToLL_2J_13TeV-amcatnloFXFX-pythia8/crab_DYToLL_2J_13TeV-amcatnloFXFX-pythia8/170323_021258 | grep root | grep -v failed >> HADD/DY0.txt")
with open('HADD/DY0.txt','r') as f:
  for line in f:
    if not line.isspace():
      DY_ch.Add(str(line[:-1]))
print "DY has", DY_ch.GetEntries(), "entries."
#VV TChains
#Single top TChains
#W+jets TChains
#tt+V TChains

# Dataset to plot
filelist = [TT_ch,DY_ch]
#Cuts and Ordering
#cut = "mu1_pt>10 && mu2_pt>10 && fabs(mu1_eta)<2.4 && fabs(mu2_eta)<2.4"
cut = ""
benchmarks = ["TTbar","DY"]

#---Starting to plot histos here-------------------------------------------------------------------------------------------------------------------------------------
#NORM="unity"
NORM="lumi"
#GEN
#RECO Di-Leptons
draw1D(filelist, "muon1_pt", "(50,10,100)", "Pt [GeV]", cut, benchmarks, "h_MU1_pt.pdf", Lumi, NORM)
draw1D(filelist, "muon2_pt", "(50,10,100)", "Pt [GeV]", cut, benchmarks, "h_MU2_pt.pdf", Lumi, NORM)
draw1D(filelist, "muon1_eta", "(50,-3,3)", "#Eta", cut, benchmarks, "h_MU1_eta.pdf", Lumi, NORM)
draw1D(filelist, "muon2_eta", "(50,-3,3)", "#Eta", cut, benchmarks, "h_MU2_eta.pdf", Lumi, NORM)
draw1D(filelist, "mass_l1l2", "(36,10,150)", "m(l,l)", cut, benchmarks, "h_M_l1l2.pdf", Lumi, NORM)
draw1D(filelist, "dR_l1l2", "(50,0,5)", "#Delta R(l,l)", cut, benchmarks, "h_DR_l1l2.pdf", Lumi, NORM)
#RECO Di-Jets
draw1D(filelist, "b1jet_pt", "(50,0,100)", "Pt [GeV]", cut, benchmarks, "h_J1_pt.pdf", Lumi, NORM)
draw1D(filelist, "b2jet_pt", "(50,0,100)", "Pt [GeV]", cut, benchmarks, "h_J2_pt.pdf", Lumi, NORM)
draw1D(filelist, "b1jet_eta", "(50,-3.,3.)", "#Eta", cut, benchmarks, "h_J1_eta.pdf", Lumi, NORM)
draw1D(filelist, "b2jet_eta", "(50,-3.,3.)", "#Eta", cut, benchmarks, "h_J2_eta.pdf", Lumi, NORM)
draw1D(filelist, "mass_b1b2", "(36,10,150)", "m(j,j)", cut, benchmarks, "h_M_b1b2.pdf", Lumi, NORM)
draw1D(filelist, "dR_b1b2", "(50,0,5)", "#Delta R(j,j)", cut, benchmarks, "h_DR_b1b2.pdf", Lumi, NORM)
#RECO MIX
#draw1D(filelist, "mass_trans", "(50,5,100)", "M_trans [GeV]", cut, benchmarks, "h_MT.pdf", Lumi, NORM)
draw1D(filelist, "dR_l1l2b1b2", "(50,0,5)", "#Delta R(lljj)", cut, benchmarks, "h_dR_l1l2b1b2.pdf", Lumi, NORM)
draw1D(filelist, "dphi_llmet", "(50,-4,4)", "#Delta #Phi (ll,MET)", cut, benchmarks, "h_Dphi_llmet.pdf", Lumi, NORM)
