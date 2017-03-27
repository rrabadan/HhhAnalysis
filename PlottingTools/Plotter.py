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
os.system("rm -rf HADD/*txt")
Lumi=36.8#fb-1

TT_ch = ROOT.TChain(tree_name)
os.system("find /fdata/hepx/store/user/lpernie/TTTo2L2Nu_13TeV-powheg/crab_TTTo2L2Nu_13TeV-powheg/170327_023645 | grep root | grep -v failed > HADD/TT_TTTo2L2Nu13TeV-powheg.txt")
os.system("cat HADD/TT_* > HADD/TT.txt")
with open("HADD/TT.txt","r") as f:
  for line in f:
    if not line.isspace():
      TT_ch.Add(str(line[:-1]))
print "TT has", TT_ch.GetEntries(), "entries."
DY_ch = ROOT.TChain(tree_name)
os.system("find /fdata/hepx/store/user/lpernie/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/170327_023657 | grep root | grep -v failed > HADD/DY_DYJetsToLLM-10to50.txt")
os.system("find /fdata/hepx/store/user/lpernie/DYToLL_0J_13TeV-amcatnloFXFX-pythia8/crab_DYToLL_0J_13TeV-amcatnloFXFX-pythia8/170327_023709 | grep root | grep -v failed > HADD/DY_DYToLL0J.txt")
os.system("find /fdata/hepx/store/user/lpernie/DYToLL_1J_13TeV-amcatnloFXFX-pythia8/crab_DYToLL_1J_13TeV-amcatnloFXFX-pythia8/170327_023721 | grep root | grep -v failed > HADD/DY_DYToLL1J.txt")
os.system("find /fdata/hepx/store/user/lpernie/DYToLL_2J_13TeV-amcatnloFXFX-pythia8/crab_DYToLL_2J_13TeV-amcatnloFXFX-pythia8/170327_023733 | grep root | grep -v failed > HADD/DY_DYToLL2J.txt")
os.system("cat HADD/DY_* > HADD/DY.txt")
with open("HADD/DY.txt","r") as f:
  for line in f:
    if not line.isspace():
      DY_ch.Add(str(line[:-1]))
print "DY has", DY_ch.GetEntries(), "entries."
VV_ch = ROOT.TChain(tree_name)
os.system("find /fdata/hepx/store/user/lpernie/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/crab_ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/170327_023746 | grep root | grep -v failed > HADD/VV_ZZTo2L2Q13TeV.txt")
os.system("find /fdata/hepx/store/user/lpernie/ZZTo2L2Nu_13TeV_powheg_pythia8/crab_ZZTo2L2Nu_13TeV_powheg_pythia8/NEW | grep root | grep -v failed > HADD/VV_ZZTo2L2Nu13TeV.txt")
os.system("find /fdata/hepx/store/user/lpernie/ZZTo4L_13TeV_powheg_pythia8/crab_ZZTo4L_13TeV_powheg_pythia8/170327_023811 | grep root | grep -v failed > HADD/VV_ZZTo4L13TeV.txt")
os.system("find /fdata/hepx/store/user/lpernie/WWToLNuQQ_aTGC_13TeV-madgraph-pythia8/crab_WWToLNuQQ_aTGC_13TeV-madgraph-pythia8/170327_023824 | grep root | grep -v failed > HADD/VV_WWToLNuQQaTGC.txt")
os.system("find /fdata/hepx/store/user/lpernie/WWTo2L2Nu_MWW-600To800_aTGC_13TeV-amcatnloFXFX-madspin-pythia8/crab_WWTo2L2Nu_MWW-600To800_aTGC_13TeV-amcatnloFXFX-madspin-pythia8/170327_023836 | grep root | grep -v failed > HADD/VV_WWTo2L2NuMWW-600To800.txt")
os.system("find /fdata/hepx/store/user/lpernie/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/crab_WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/170327_023848 | grep root | grep -v failed > HADD/VV_WZTo2L2Q13TeV.txt")
os.system("find /fdata/hepx/store/user/lpernie/WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8/crab_WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8/170327_023901 | grep root | grep -v failed > HADD/VV_WZTo1L3Nu13TeV.txt")
os.system("find /fdata/hepx/store/user/lpernie/WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/crab_WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/170327_023913 | grep root | grep -v failed > HADD/VV_WZTo1L1Nu2Q13TeV.txt")
os.system("cat HADD/VV_* > HADD/VV.txt")
with open("HADD/VV.txt","r") as f:
  for line in f:
    if not line.isspace():
      VV_ch.Add(str(line[:-1]))
print "VV has", VV_ch.GetEntries(), "entries."
sT_ch = ROOT.TChain(tree_name)
os.system("find /fdata/hepx/store/user/lpernie/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/crab_ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/170327_025236 | grep root | grep -v failed > HADD/sT_STt-channel.txt")
os.system("find /fdata/hepx/store/user/lpernie/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/crab_ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/170327_025249 | grep root | grep -v failed > HADD/sT_STt-channel.txt")
os.system("find /fdata/hepx/store/user/lpernie/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/crab_ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/170327_025301 | grep root | grep -v failed > HADD/sT_STs-channel.txt")
os.system("find /fdata/hepx/store/user/lpernie/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/crab_ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/170327_025313 | grep root | grep -v failed > HADD/sT_STtW.txt")
os.system("find /fdata/hepx/store/user/lpernie/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/crab_ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/170327_025326 | grep root | grep -v failed > HADD/sT_STtW.txt")
os.system("cat HADD/sT_* > HADD/sT.txt")
with open("HADD/sT.txt","r") as f:
  for line in f:
    if not line.isspace():
      sT_ch.Add(str(line[:-1]))
print "sT has", sT_ch.GetEntries(), "entries."
Wjet_ch = ROOT.TChain(tree_name)
os.system("find /fdata/hepx/store/user/lpernie/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/170327_025339 | grep root | grep -v failed > HADD/Wjet_WJetsToLNuTuneCUETP8M1.txt")
os.system("find /fdata/hepx/store/user/lpernie/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/170327_025352 | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-100To200.txt")
os.system("find /fdata/hepx/store/user/lpernie/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/170327_025404 | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-200To400.txt")
os.system("find /fdata/hepx/store/user/lpernie/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/170327_025418 | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-400To600.txt")
os.system("find /fdata/hepx/store/user/lpernie/WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/170327_025430 | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-600To800.txt")
os.system("find /fdata/hepx/store/user/lpernie/WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/170327_025443 | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-800To1200.txt")
os.system("find /fdata/hepx/store/user/lpernie/WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/170327_025457 | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-1200To2500.txt")
os.system("find /fdata/hepx/store/user/lpernie/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/170327_025509 | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-2500ToInf.txt")
os.system("cat HADD/Wjet_* > HADD/Wjet.txt")
with open("HADD/Wjet.txt","r") as f:
  for line in f:
    if not line.isspace():
      Wjet_ch.Add(str(line[:-1]))
print "Wjet has", Wjet_ch.GetEntries(), "entries."
ttV_ch = ROOT.TChain(tree_name)
os.system("find /fdata/hepx/store/user/lpernie/TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/170327_025521 | grep root | grep -v failed > HADD/ttV_TTWJetsToQQTuneCUETP8M1.txt")
os.system("find /fdata/hepx/store/user/lpernie/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/170327_025533 | grep root | grep -v failed > HADD/ttV_TTWJetsToLNuTuneCUETP8M1.txt")
os.system("find /fdata/hepx/store/user/lpernie/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/170327_025546 | grep root | grep -v failed > HADD/ttV_TTZToQQTuneCUETP8M1.txt")
os.system("cat HADD/ttV_* > HADD/ttV.txt")
with open("HADD/ttV.txt","r") as f:
  for line in f:
    if not line.isspace():
      ttV_ch.Add(str(line[:-1]))
print "ttV has", ttV_ch.GetEntries(), "entries."
os.system("find /fdata/hepx/store/user/lpernie/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/170327_025558 | grep root | grep -v failed > HADD/ttV_TTZToLLNuNuM-10.txt")

# Dataset to plot
filelist = [TT_ch, DY_ch, VV_ch, sT_ch, Wjet_ch, ttV_ch]
#Cuts and Ordering
#cut = "mu1_pt>10 && mu2_pt>10 && fabs(mu1_eta)<2.4 && fabs(mu2_eta)<2.4"
cut = ""
benchmarks = ["TTbar","DY","VV","singTop","Wjet","ttV"]

#---Starting to plot histos here-------------------------------------------------------------------------------------------------------------------------------------
NORM="unity"
#NORM="lumi"
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
