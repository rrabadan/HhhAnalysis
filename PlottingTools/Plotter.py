import random
import ROOT
from ROOT import TFile,TChain,TH1F,TH2F,TLegend
import os
execfile("start.py")
execfile("functions.py")

#Creating folders and parameters
tree_name="DiHiggsWWBBAna/evtree"
os.system("mkdir -p Plots/C")
os.system("mkdir -p HADD")
os.system("rm -rf HADD/*txt")
#Lumi=36.42#fb-1
#Lumi=13.828#fb-1
Lumi=10.216#fb-1


# MC
MCprefix = "/fdata/hepx/store/user/tahuang/MC_Hhh_analysis_20171010/"
MCprefix_v2 = "/fdata/hepx/store/user/tahuang/MC_Hhh_analysis_20171010_v2/"
HLTntuples = MCprefix+"HLTNtuples/"
makeHadd = True
TT_ch = ROOT.TChain(tree_name)
#os.system("find /fdata/hepx/store/user/lpernie/TTTo2L2Nu_13TeV-powheg/crab_TTTo2L2Nu_13TeV-powheg/170401_035903 | grep root | grep -v failed > HADD/TT_TTTo2L2Nu13TeV-powheg.txt")
os.system("find "+MCprefix+"TTTo2L2Nu_13TeV-powheg/ | grep root | grep -v failed > HADD/TT_TTTo2L2Nu13TeV-powheg.txt")
print "find command ", "find "+MCprefix+"crab_TTTo2L2Nu_13TeV-powheg/ | grep root | grep -v failed > HADD/TT_TTTo2L2Nu13TeV-powheg.txt"
os.system("cat HADD/TT_* > HADD/TT.txt")
#if makeHadd: os.system("hadd -T -f -k /fdata/hepx/store/user/lpernie/TTTo2L2Nu_13TeV-powheg/crab_TTTo2L2Nu_13TeV-powheg.root @HADD/TT.txt")
if makeHadd: os.system("hadd -T -f -k "+ HLTntuples +"crab_TTTo2L2Nu_13TeV-powheg.root @HADD/TT.txt")
#N_tot_path_TT = "/fdata/hepx/store/user/lpernie/TTTo2L2Nu_13TeV-powheg/crab_TTTo2L2Nu_13TeV-powheg.root"
N_tot_path_TT = HLTntuples+"crab_TTTo2L2Nu_13TeV-powheg.root"
TT_file =  ROOT.TFile.Open(N_tot_path_TT,"read"); h_prehlt_TT =  ROOT.TH1F(TT_file.Get("TriggerResults/hevent_filter")); nTOT_prehlt_TT = h_prehlt_TT.GetBinContent(2)
h_posthlt_TT =  ROOT.TH1F(TT_file.Get("DiHiggsWWBBAna/hevent")); nTOT_posthlt_TT = h_posthlt_TT.GetBinContent(2);
with open("HADD/TT.txt","r") as f:
  for line in f:
    if not line.isspace():
      TT_ch.Add(str(line[:-1]))
print "TT has", TT_ch.GetEntries(), "entries."

DY_ch = ROOT.TChain(tree_name)
os.system("find "+MCprefix+"DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/ | grep root | grep -v failed > HADD/DY_DYJetsToLLM-10to50.txt")
os.system("find "+MCprefix+"DYToLL_0J_13TeV-amcatnloFXFX-pythia8/ | grep root | grep -v failed > HADD/DY_DYToLL0J.txt")
os.system("find "+MCprefix+"DYToLL_1J_13TeV-amcatnloFXFX-pythia8/ | grep root | grep -v failed > HADD/DY_DYToLL1J.txt")
os.system("find "+MCprefix+"DYToLL_2J_13TeV-amcatnloFXFX-pythia8/ | grep root | grep -v failed > HADD/DY_DYToLL2J.txt")
#os.system("find /fdata/hepx/store/user/lpernie/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/170401_035916 | grep root | grep -v failed > HADD/DY_DYJetsToLLM-10to50.txt")
#os.system("find /fdata/hepx/store/user/lpernie/DYToLL_0J_13TeV-amcatnloFXFX-pythia8/crab_DYToLL_0J_13TeV-amcatnloFXFX-pythia8/170401_035928 | grep root | grep -v failed > HADD/DY_DYToLL0J.txt")
#os.system("find /fdata/hepx/store/user/lpernie/DYToLL_1J_13TeV-amcatnloFXFX-pythia8/crab_DYToLL_1J_13TeV-amcatnloFXFX-pythia8/170401_035940 | grep root | grep -v failed > HADD/DY_DYToLL1J.txt")
#os.system("find /fdata/hepx/store/user/lpernie/DYToLL_2J_13TeV-amcatnloFXFX-pythia8/crab_DYToLL_2J_13TeV-amcatnloFXFX-pythia8/170401_035952 | grep root | grep -v failed > HADD/DY_DYToLL2J.txt")
os.system("cat HADD/DY_* > HADD/DY.txt")
#if makeHadd: os.system("hadd -T -f -k /fdata/hepx/store/user/lpernie/DYToLL_2J_13TeV-amcatnloFXFX-pythia8/crab_DYToLL_2J_13TeV-amcatnloFXFX-pythia8.root @HADD/DY.txt")
if makeHadd: os.system("hadd -T -f -k "+HLTntuples+"crab_DYToLL_2J_13TeV-amcatnloFXFX-pythia8.root @HADD/DY.txt")

#N_tot_path_DY = "/fdata/hepx/store/user/lpernie/DYToLL_2J_13TeV-amcatnloFXFX-pythia8/crab_DYToLL_2J_13TeV-amcatnloFXFX-pythia8.root"
N_tot_path_DY = HLTntuples+"crab_DYToLL_2J_13TeV-amcatnloFXFX-pythia8.root"
DY_file =  ROOT.TFile.Open(N_tot_path_DY,"read"); h_prehlt_DY =  ROOT.TH1F(DY_file.Get("TriggerResults/hevent_filter")); nTOT_prehlt_DY = h_prehlt_DY.GetBinContent(2)
h_posthlt_DY =  ROOT.TH1F(DY_file.Get("DiHiggsWWBBAna/hevent")); nTOT_posthlt_DY = h_posthlt_DY.GetBinContent(2);
with open("HADD/DY.txt","r") as f:
  for line in f:
    if not line.isspace():
      DY_ch.Add(str(line[:-1]))
print "DY has", DY_ch.GetEntries(), "entries."

sT_ch = ROOT.TChain(tree_name)
os.system("find "+ MCprefix_v2 + "ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/ | grep root | grep -v failed > HADD/sT_STt-channel.txt")
os.system("find "+ MCprefix_v2 + "ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/ | grep root | grep -v failed > sT_SantiTt-channel.txt")
os.system("find "+ MCprefix + "ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/ | grep root | grep -v failed > HADD/sT_STs-channel.txt")
os.system("find "+ MCprefix + "ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/ | grep root | grep -v failed > HADD/sT_SantiTtW.txt")
os.system("find "+ MCprefix + "ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/  | grep root | grep -v failed > HADD/sT_STtW.txt")
#os.system("find /fdata/hepx/store/user/lpernie/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/crab_ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/170401_040200 | grep root | grep -v failed > HADD/sT_STt-channel.txt")
#os.system("find /fdata/hepx/store/user/lpernie/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/crab_ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/170401_040213 | grep root | grep -v failed > HADD/sT_STt-channel.txt")
#os.system("find /fdata/hepx/store/user/lpernie/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/crab_ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/170401_040225 | grep root | grep -v failed > HADD/sT_STs-channel.txt")
#os.system("find /fdata/hepx/store/user/lpernie/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/crab_ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/170401_040237 | grep root | grep -v failed > HADD/sT_STtW.txt")
#os.system("find /fdata/hepx/store/user/lpernie/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/crab_ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/170401_040249 | grep root | grep -v failed > HADD/sT_STtW.txt")
os.system("cat HADD/sT_* > HADD/sT.txt")
#if makeHadd: os.system("hadd -T -f -k /fdata/hepx/store/user/lpernie/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/crab_ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1.root @HADD/sT.txt")
if makeHadd: os.system("hadd -T -f -k "+ HLTntuples +"crab_ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1.root @HADD/sT.txt")
N_tot_path_sT = HLTntuples+"crab_ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1.root"
sT_file =  ROOT.TFile.Open(N_tot_path_sT,"read"); h_prehlt_sT =  ROOT.TH1F(sT_file.Get("TriggerResults/hevent_filter")); nTOT_prehlt_sT = h_prehlt_sT.GetBinContent(2)
h_posthlt_sT =  ROOT.TH1F(sT_file.Get("DiHiggsWWBBAna/hevent")); nTOT_posthlt_sT = h_posthlt_sT.GetBinContent(2);
with open("HADD/sT.txt","r") as f:
  for line in f:
    if not line.isspace():
      sT_ch.Add(str(line[:-1]))
print "sT has", sT_ch.GetEntries(), "entries."

VV_ch = ROOT.TChain(tree_name)
os.system("find "+ MCprefix +"ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/ | grep root | grep -v failed > HADD/VV_ZZTo2L2Q13TeV.txt")
os.system("find "+ MCprefix +"ZZTo2L2Nu_13TeV_powheg_pythia8/ | grep root | grep -v failed >HADD/VV_ZZTo2L2Q13TeV.txt")
os.system("find "+ MCprefix +"ZZTo4L_13TeV_powheg_pythia8/ | grep root | grep -v failed >  HADD/VV_ZZTo4L13TeV.txt")
os.system("find "+ MCprefix_v2 +"WWToLNuQQ_aTGC_13TeV-madgraph-pythia8/ | grep root | grep -v failed > HADD/VV_WWToLNuQQaTGC.txt")
os.system("find "+ MCprefix +"WWTo2L2Nu_MWW-600To800_aTGC_13TeV-amcatnloFXFX-madspin-pythia8 | grep root | grep -v failed > HADD/VV_WWTo2L2NuMWW-600To800.txt")
os.system("find "+ MCprefix +"WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/ | grep root | grep -v failed >HADD/VV_WZTo2L2Q13TeV.txt")
os.system("find "+ MCprefix +"WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8/ | grep root | grep -v failed > HADD/VV_WZTo1L3Nu13TeV.txt")
os.system("find "+ MCprefix +"WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/ | grep root | grep -v failed > HADD/VV_WZTo1L3Nu13TeV.txt")
os.system("find "+ MCprefix +"WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/ | grep root | grep -v failed > VV_WZTo3LNuTuneCUETP8M1.txt")
#os.system("find /fdata/hepx/store/user/lpernie/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/crab_ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/170401_040005 | grep root | grep -v failed > HADD/VV_ZZTo2L2Q13TeV.txt")
#os.system("find /fdata/hepx/store/user/lpernie/ZZTo2L2Nu_13TeV_powheg_pythia8/crab_ZZTo2L2Nu_13TeV_powheg_pythia8/170401_040017 | grep root | grep -v failed > HADD/VV_ZZTo2L2Nu13TeV.txt")
#os.system("find /fdata/hepx/store/user/lpernie/ZZTo4L_13TeV_powheg_pythia8/crab_ZZTo4L_13TeV_powheg_pythia8/170401_040030 | grep root | grep -v failed > HADD/VV_ZZTo4L13TeV.txt")
#os.system("find /fdata/hepx/store/user/lpernie/WWToLNuQQ_aTGC_13TeV-madgraph-pythia8/crab_WWToLNuQQ_aTGC_13TeV-madgraph-pythia8/170401_040043 | grep root | grep -v failed > HADD/VV_WWToLNuQQaTGC.txt")
#os.system("find /fdata/hepx/store/user/lpernie/WWTo2L2Nu_MWW-600To800_aTGC_13TeV-amcatnloFXFX-madspin-pythia8/crab_WWTo2L2Nu_MWW-600To800_aTGC_13TeV-amcatnloFXFX-madspin-pythia8/170401_040055 | grep root | grep -v failed > HADD/VV_WWTo2L2NuMWW-600To800.txt")
#os.system("find /fdata/hepx/store/user/lpernie/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/crab_WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/170401_040108 | grep root | grep -v failed > HADD/VV_WZTo2L2Q13TeV.txt")
#os.system("find /fdata/hepx/store/user/lpernie/WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8/crab_WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8/170401_040121 | grep root | grep -v failed > HADD/VV_WZTo1L3Nu13TeV.txt")
#os.system("find /fdata/hepx/store/user/lpernie/WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/crab_WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/170401_040134 | grep root | grep -v failed > HADD/VV_WZTo1L1Nu2Q13TeV.txt")
#os.system("find /fdata/hepx/store/user/lpernie/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/crab_WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/170401_040146 | grep root | grep -v failed > HADD/VV_WZTo3LNuTuneCUETP8M1.txt")
os.system("cat HADD/VV_* > HADD/VV.txt")
#if makeHadd: os.system("hadd -T -f -k /fdata/hepx/store/user/lpernie/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/crab_WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8.root @HADD/VV.txt")
if makeHadd: os.system("hadd -T -f -k "+HLTntuples+"crab_WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8.root @HADD/VV.txt")
#N_tot_path_VV = "/fdata/hepx/store/user/lpernie/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/crab_WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8.root"
N_tot_path_VV = HLTntuples+"crab_WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8.root"
VV_file =  ROOT.TFile.Open(N_tot_path_VV,"read"); h_prehlt_VV =  ROOT.TH1F(VV_file.Get("TriggerResults/hevent_filter")); nTOT_prehlt_VV = h_prehlt_VV.GetBinContent(2)
h_posthlt_VV =  ROOT.TH1F(VV_file.Get("DiHiggsWWBBAna/hevent")); nTOT_posthlt_VV = h_posthlt_VV.GetBinContent(2);
with open("HADD/VV.txt","r") as f:
  for line in f:
    if not line.isspace():
      VV_ch.Add(str(line[:-1]))
print "VV has", VV_ch.GetEntries(), "entries."

Wjet_ch = ROOT.TChain(tree_name)
os.system("find "+ MCprefix +"WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/ | grep root | grep -v failed > HADD/Wjet_WJetsToLNuTuneCUETP8M1.txt")
os.system("find "+ MCprefix_v2 +"WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-100To200.txt")
os.system("find "+ MCprefix_v2 +"WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-200To400.txt")
os.system("find "+ MCprefix +"WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-400To600.txt")
os.system("find "+ MCprefix_v2 +"WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-600To800.txt")
os.system("find "+ MCprefix +"WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-800To1200.txt")
os.system("find "+ MCprefix +"WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-1200To2500.txt")
os.system("find "+ MCprefix +"WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-2500ToInf.txt")
#os.system("find /fdata/hepx/store/user/lpernie/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/170401_040303 | grep root | grep -v failed > HADD/Wjet_WJetsToLNuTuneCUETP8M1.txt")
#os.system("find /fdata/hepx/store/user/lpernie/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/170401_040315 | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-100To200.txt")
#os.system("find /fdata/hepx/store/user/lpernie/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/170401_040328 | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-200To400.txt")
#os.system("find /fdata/hepx/store/user/lpernie/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/170401_040341 | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-400To600.txt")
#os.system("find /fdata/hepx/store/user/lpernie/WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/170401_040353 | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-600To800.txt")
#os.system("find /fdata/hepx/store/user/lpernie/WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/170401_040408 | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-800To1200.txt")
#os.system("find /fdata/hepx/store/user/lpernie/WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/170401_040421 | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-1200To2500.txt")
#os.system("find /fdata/hepx/store/user/lpernie/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/170401_040433 | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-2500ToInf.txt")
os.system("cat HADD/Wjet_* > HADD/Wjet.txt")
#if makeHadd: os.system("hadd -T -f -k /fdata/hepx/store/user/lpernie/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root @HADD/Wjet.txt")
#N_tot_path_Wjet = "/fdata/hepx/store/user/lpernie/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"
if makeHadd: os.system("hadd -T -f -k "+ HLTntuples +"crab_WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root @HADD/Wjet.txt")
N_tot_path_Wjet = HLTntuples + "crab_WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"
Wjet_file =  ROOT.TFile.Open(N_tot_path_Wjet,"read"); h_prehlt_Wjet =  ROOT.TH1F(Wjet_file.Get("TriggerResults/hevent_filter")); nTOT_prehlt_Wjet = h_prehlt_Wjet.GetBinContent(2)
h_posthlt_Wjet =  ROOT.TH1F(Wjet_file.Get("DiHiggsWWBBAna/hevent")); nTOT_posthlt_Wjet = h_posthlt_Wjet.GetBinContent(2);
with open("HADD/Wjet.txt","r") as f:
  for line in f:
    if not line.isspace():
      Wjet_ch.Add(str(line[:-1]))
print "Wjet has", Wjet_ch.GetEntries(), "entries."

ttV_ch = ROOT.TChain(tree_name)
os.system("find "+ MCprefix +"TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/ | grep root | grep -v failed > HADD/ttV_TTWJetsToQQTuneCUETP8M1.txt")
os.system("find "+ MCprefix +"TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/ | grep root | grep -v failed > HADD/ttV_TTWJetsToLNuTuneCUETP8M1.txt")
os.system("find "+ MCprefix +"TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/ | grep root | grep -v failed > HADD/ttV_TTZToQQTuneCUETP8M1.txt")
os.system("find "+ MCprefix +"TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/ | grep root | grep -v failed > HADD/ttV_TTZToLLNuNuM-10.txt")
#os.system("find /fdata/hepx/store/user/lpernie/TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/170401_040446 | grep root | grep -v failed > HADD/ttV_TTWJetsToQQTuneCUETP8M1.txt")
#os.system("find /fdata/hepx/store/user/lpernie/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/170401_040458 | grep root | grep -v failed > HADD/ttV_TTWJetsToLNuTuneCUETP8M1.txt")
#os.system("find /fdata/hepx/store/user/lpernie/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/170401_040512 | grep root | grep -v failed > HADD/ttV_TTZToQQTuneCUETP8M1.txt")
#os.system("find /fdata/hepx/store/user/lpernie/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/170401_040525 | grep root | grep -v failed > HADD/ttV_TTZToLLNuNuM-10.txt")
os.system("cat HADD/ttV_* > HADD/ttV.txt")
if makeHadd: os.system("hadd -T -f -k "+ HLTntuples +"crab_TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root @HADD/ttV.txt")
N_tot_path_ttV = HLTntuples+"crab_TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root"
#if makeHadd: os.system("hadd -T -f -k /fdata/hepx/store/user/lpernie/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root @HADD/ttV.txt")
#N_tot_path_ttV = "/fdata/hepx/store/user/lpernie/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root"
ttV_file =  ROOT.TFile.Open(N_tot_path_ttV,"read"); h_prehlt_ttV =  ROOT.TH1F(ttV_file.Get("TriggerResults/hevent_filter")); nTOT_prehlt_ttV = h_prehlt_ttV.GetBinContent(2)
h_posthlt_ttV =  ROOT.TH1F(ttV_file.Get("DiHiggsWWBBAna/hevent")); nTOT_posthlt_ttV = h_posthlt_ttV.GetBinContent(2);
with open("HADD/ttV.txt","r") as f:
  for line in f:
    if not line.isspace():
      ttV_ch.Add(str(line[:-1]))
print "ttV has", ttV_ch.GetEntries(), "entries."


# Data
DATA_ch = ROOT.TChain(tree_name)
#Dataprefix = "/fdata/hepx/store/user/tahuang/DoubleMuon/"
#os.system("find "+Dataprefix+"crab_Hhh_Run2017B-PromptReco-v1/ | grep root | grep -v failed > HADD/DATA_Hhh_Run2017B-v1.txt")
#os.system("find "+Dataprefix+"crab_Hhh_Run2017C-PromptReco-v1/ | grep root | grep -v failed > HADD/DATA_Hhh_Run2017C-v1.txt")
#os.system("find "+Dataprefix+"crab_Hhh_Run2017C-PromptReco-v2/ | grep root | grep -v failed > HADD/DATA_Hhh_Run2017C-v2.txt")
#os.system("find "+Dataprefix+"crab_Hhh_Run2017C-PromptReco-v3/ | grep root | grep -v failed > HADD/DATA_Hhh_Run2017C-v3.txt")
#os.system("find /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Hhh_Run2016B-23Sep2016-v3/170329_223804 | grep root | grep -v failed > HADD/DATA_Hhh_Run2016B-23Sep2016-v3.txt")
Dataprefix = "/fdata/hepx/store/user/tahuang/Collisioins17/DoubleMuon/"
os.system("find "+Dataprefix+"crab_crabtest_2017Bv1 | grep root | grep -v failed > HADD/DATA_Hhh_Run2017B.txt")
os.system("find "+Dataprefix+"crab_crabtest_v2 | grep root | grep -v failed > HADD/DATA_Hhh_Run2017C.txt")
#os.system("find /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Hhh_Run2016C-23Sep2016-v1/170329_223816 | grep root | grep -v failed > HADD/DATA_Hhh_Run2016C-23Sep2016-v1.txt")
#os.system("find /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Hhh_Run2016D-23Sep2016-v1/170329_223828 | grep root | grep -v failed > HADD/DATA_Hhh_Run2016D-23Sep2016-v1.txt")
#os.system("find /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Hhh_Run2016E-23Sep2016-v1/170329_223841 | grep root | grep -v failed > HADD/DATA_Hhh_Run2016E-23Sep2016-v1.txt")
#os.system("find /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Hhh_Run2016F-23Sep2016-v1/170329_223852 | grep root | grep -v failed > HADD/DATA_Hhh_Run2016F-23Sep2016-v1.txt")
#os.system("find /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Hhh_Run2016G-23Sep2016-v1/170329_223905 | grep root | grep -v failed > HADD/DATA_Hhh_Run2016G-23Sep2016-v1.txt")
#os.system("find /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Hhh_Run2016H-03Feb2017_ver2-v1/170330_194716 | grep root | grep -v failed > HADD/DATA_Hhh_Run2016H-03Feb2017_ver2-v1.txt")
#os.system("find /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Hhh_Run2016H-03Feb2017_ver3-v1/170330_194729 | grep root | grep -v failed > HADD/DATA_Hhh_Run2016H-03Feb2017_ver3-v1.txt")
os.system("cat HADD/DATA_* > HADD/DATA.txt")
with open("HADD/DATA.txt","r") as f:
  for line in f:
    if not line.isspace():
      DATA_ch.Add(str(line[:-1]))
print "DATA has", DATA_ch.GetEntries(), "entries."

# Dataset to plot
#filelist   = [ttV_ch, Wjet_ch, sT_ch, VV_ch, DY_ch, TT_ch, DATA_ch] #If you draw a StackPlot place the smaller samples at the beginning and data as last one
#benchmarks = ["ttV", "Wjet", "singTop", "VV", "DY", "TTbar", "Data"]
#nTOT       = [nTOT_prehlt_ttV, nTOT_prehlt_Wjet, nTOT_prehlt_sT, nTOT_prehlt_VV, nTOT_prehlt_DY, nTOT_prehlt_TT, 1]
filelist   = [Wjet_ch, VV_ch, sT_ch, DY_ch, TT_ch, DATA_ch] #If you draw a StackPlot place the smaller samples at the beginning and data as last one
benchmarks = ["Wjet","VV","singleTop","DY", "TTbar", "Data"]
nTOT       = [nTOT_prehlt_Wjet, nTOT_prehlt_VV, nTOT_prehlt_sT, nTOT_prehlt_DY, nTOT_prehlt_TT, 1]
#Cuts and Ordering
cleaningcut = "mass_l1l2<75  && mass_l1l2>12 && dR_l1l2<3.8 && mass_b1b2>30 && dR_b1b2<3.8 && dR_l1l2b1b2>0.2 && dR_l1l2b1b2<4.0"
cut = "met_pt>20 && b1jet_pt>30 && TMath::Abs(b1jet_eta)<2.5 && b2jet_pt>30 && TMath::Abs(b2jet_eta)<2.5 && muon1_pt>20 && TMath::Abs(muon1_eta)<2.4 && muon2_pt>20 && TMath::Abs(muon2_eta)<2.4"+" && "+cleaningcut

#---Starting to plot histos here-------------------------------------------------------------------------------------------------------------------------------------
NORM=["lumi"] #NORM=["unity","lumi"]
DataOrMC="DataMC" # MC or DATA = all samples overimposed, DataMC = MC is stack, data is overimposed
LOG=["no"]#,"no"] # You can choose if you want the plots to be log, linear or both
Format=[".pdf",".C"]
Format=[".pdf"]
for this_LOG in LOG:
  for this_Norm in NORM:
    #GEN
    #RECO Di-Leptons
    draw1D(filelist, "muon1_pt", "(50,10.,210.)", "Lead. #mu P_{T} [GeV]", cut, benchmarks, "h_MU1_pt", Lumi, nTOT, this_Norm, DataOrMC, this_LOG, Format)
    draw1D(filelist, "muon2_pt", "(50,10.,210.)", "Sublead. #mu P_{T} [GeV]", cut, benchmarks, "h_MU2_pt", Lumi, nTOT, this_Norm, DataOrMC, this_LOG, Format)
    draw1D(filelist, "muon1_eta", "(50,-3.,3.)", "Lead. #mu #eta", cut, benchmarks, "h_MU1_eta", Lumi, nTOT, this_Norm, DataOrMC, this_LOG, Format)
    draw1D(filelist, "muon2_eta", "(50,-3.,3.)", "Sublead. #mu #eta", cut, benchmarks, "h_MU2_eta", Lumi, nTOT, this_Norm, DataOrMC, this_LOG, Format)
    draw1D(filelist, "mass_l1l2", "(40,10.,80.)", "m(l,l)", cut, benchmarks, "h_M_l1l2", Lumi, nTOT, this_Norm, DataOrMC, this_LOG, Format)
    draw1D(filelist, "pt_l1l2", "(40,0.,300.)", "p_{T}(l,l)", cut, benchmarks, "h_pt_l1l2", Lumi, nTOT, this_Norm, DataOrMC, this_LOG, Format)
    draw1D(filelist, "dR_l1l2", "(50,0.,5.)", "#Delta R(l,l)", cut, benchmarks, "h_DR_l1l2", Lumi, nTOT, this_Norm, DataOrMC, this_LOG, Format)
    #RECO Di-Jets
    draw1D(filelist, "b1jet_pt", "(50,20.,350.)", "Lead. jet P_{T} [GeV]", cut, benchmarks, "h_J1_pt", Lumi, nTOT, this_Norm, DataOrMC, this_LOG, Format)
    draw1D(filelist, "b2jet_pt", "(50,20.,350.)", "Sublead. jet P_{T} [GeV]", cut, benchmarks, "h_J2_pt", Lumi, nTOT, this_Norm, DataOrMC, this_LOG, Format)
    draw1D(filelist, "b1jet_eta", "(50,-3.,3.)", "Lead. jet #eta", cut, benchmarks, "h_J1_eta", Lumi, nTOT, this_Norm, DataOrMC, this_LOG, Format)
    draw1D(filelist, "b2jet_eta", "(50,-3.,3.)", "Sublead. jet #eta", cut, benchmarks, "h_J2_eta", Lumi, nTOT, this_Norm, DataOrMC, this_LOG, Format)
    draw1D(filelist, "mass_b1b2", "(40,40.,400.)", "m(j,j)", cut, benchmarks, "h_M_b1b2", Lumi, nTOT, this_Norm, DataOrMC, this_LOG, Format)
    draw1D(filelist, "pt_b1b2", "(50,0.,400.)", "p_{T}(j,j)", cut, benchmarks, "h_pt_b1b2", Lumi, nTOT, this_Norm, DataOrMC, this_LOG, Format)
    draw1D(filelist, "dR_b1b2", "(50,0.,5.)", "#Delta R(j,j)", cut, benchmarks, "h_DR_b1b2", Lumi, nTOT, this_Norm, DataOrMC, this_LOG, Format)
    #RECO MIX
    draw1D(filelist, "met_pt", "(50,10.,400.)", "MET [GeV]", cut, benchmarks, "h_MET", Lumi, nTOT, this_Norm, DataOrMC, this_LOG, Format)
    draw1D(filelist, "mass_trans", "(50,0.,250.)", "M_{trans} [GeV]", cut, benchmarks, "h_MT", Lumi, nTOT, this_Norm, DataOrMC, this_LOG, Format)
    draw1D(filelist, "dR_l1l2b1b2", "(50,0.,5.)", "#Delta R(ll,jj)", cut, benchmarks, "h_dR_l1l2b1b2", Lumi, nTOT, this_Norm, DataOrMC, this_LOG, Format)
    draw1D(filelist, "dR_minbl", "(50,0.,5.)", "min#Delta R(l,j)", cut, benchmarks, "h_dR_minbl", Lumi, nTOT, this_Norm, DataOrMC, this_LOG, Format)
    draw1D(filelist, "TMath::Abs(dphi_llmet)", "(40,0,3.2)", "#Delta #phi (ll,MET)", cut, benchmarks, "h_Dphi_llmet", Lumi, nTOT, this_Norm, DataOrMC, this_LOG, Format)
