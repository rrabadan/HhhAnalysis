import os, sys, random
import ROOT
from ROOT import TFile,TChain,TH1F,TH2F,TLegend
execfile("start.py")
execfile("functions.py")
#Creating folders and parameters
tree_name="DiHiggsWWBBAna/evtree"
os.system("rm -rf HADD/*txt")
print "Executing: python", sys.argv[0] , "-b", sys.argv[2], sys.argv[3], "(Arg1=makeHadd: HaddYes or HaddNo /|\ Arg2=whichSample: TT, DY, VV, singTop, Wjet, ttV, Data)"
makeHadd = sys.argv[2]
whichSample = sys.argv[3]
if( makeHadd!="HaddYes" and makeHadd!="HaddNo" ): print "WARNING! 1st argument has to be HaddYes or HaddNo"; sys.exit()
if( whichSample!="TT" and whichSample!="DY" and whichSample!="VV" and whichSample!="singTop" and whichSample!="Wjet" and whichSample!="ttV" and whichSample!="Data" ):  print "WARNING! 2nd argument have to be TT, DY, VV, singTop, Wjet, ttV, or Data"; sys.exit()

Find_str = []; this_cat = ""; this_hadd = ""; this_NtotPath = ""
# MC
if( whichSample == "TT" ):
  Find_str.append("find /fdata/hepx/store/user/lpernie/TTTo2L2Nu_13TeV-powheg/crab_TTTo2L2Nu_13TeV-powheg/170405_172215 | grep root | grep -v failed > HADD/TT_TTTo2L2Nu13TeV-powheg.txt")
  this_cat      = "cat HADD/TT_* > HADD/TT.txt"
  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/TTTo2L2Nu_13TeV-powheg/crab_TTTo2L2Nu_13TeV-powheg.root @HADD/TT.txt"
  this_NtotPath = "/fdata/hepx/store/user/lpernie/TTTo2L2Nu_13TeV-powheg/crab_TTTo2L2Nu_13TeV-powheg.root"
if( whichSample == "DY" ):
  Find_str.append("find /fdata/hepx/store/user/lpernie/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/170405_172231 | grep root | grep -v failed > HADD/DY_DYJetsToLLM-10to50.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/DYToLL_0J_13TeV-amcatnloFXFX-pythia8/crab_DYToLL_0J_13TeV-amcatnloFXFX-pythia8/170405_172248 | grep root | grep -v failed > HADD/DY_DYToLL0J.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/DYToLL_1J_13TeV-amcatnloFXFX-pythia8/crab_DYToLL_1J_13TeV-amcatnloFXFX-pythia8/170405_172303 | grep root | grep -v failed > HADD/DY_DYToLL1J.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/DYToLL_2J_13TeV-amcatnloFXFX-pythia8/crab_DYToLL_2J_13TeV-amcatnloFXFX-pythia8/170405_172320 | grep root | grep -v failed > HADD/DY_DYToLL2J.txt")
  this_cat      = "cat HADD/DY_* > HADD/DY.txt"
  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/DYToLL_2J_13TeV-amcatnloFXFX-pythia8/crab_DYToLL_2J_13TeV-amcatnloFXFX-pythia8.root @HADD/DY.txt"
  this_NtotPath = "/fdata/hepx/store/user/lpernie/DYToLL_2J_13TeV-amcatnloFXFX-pythia8/crab_DYToLL_2J_13TeV-amcatnloFXFX-pythia8.root"
if( whichSample == "VV" ):
  Find_str.append("find /fdata/hepx/store/user/lpernie/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/crab_ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/170405_172338 | grep root | grep -v failed > HADD/VV_ZZTo2L2Q13TeV.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/ZZTo2L2Nu_13TeV_powheg_pythia8/crab_ZZTo2L2Nu_13TeV_powheg_pythia8/170405_172352 | grep root | grep -v failed > HADD/VV_ZZTo2L2Nu13TeV.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/ZZTo4L_13TeV_powheg_pythia8/crab_ZZTo4L_13TeV_powheg_pythia8/170405_172407 | grep root | grep -v failed > HADD/VV_ZZTo4L13TeV.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/WWToLNuQQ_aTGC_13TeV-madgraph-pythia8/crab_WWToLNuQQ_aTGC_13TeV-madgraph-pythia8/170405_172423 | grep root | grep -v failed > HADD/VV_WWToLNuQQaTGC.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/WWTo2L2Nu_MWW-600To800_aTGC_13TeV-amcatnloFXFX-madspin-pythia8/crab_WWTo2L2Nu_MWW-600To800_aTGC_13TeV-amcatnloFXFX-madspin-pythia8/170405_172438 | grep root | grep -v failed > HADD/VV_WWTo2L2NuMWW-600To800.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/crab_WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/170405_172452 | grep root | grep -v failed > HADD/VV_WZTo2L2Q13TeV.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8/crab_WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8/170405_172508 | grep root | grep -v failed > HADD/VV_WZTo1L3Nu13TeV.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/crab_WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/170405_172522 | grep root | grep -v failed > HADD/VV_WZTo1L1Nu2Q13TeV.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/crab_WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/170405_172539 | grep root | grep -v failed > HADD/VV_WZTo3LNuTuneCUETP8M1.txt")
  this_cat      = "cat HADD/VV_* > HADD/VV.txt"
  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/crab_WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8.root @HADD/VV.txt"
  this_NtotPath = "/fdata/hepx/store/user/lpernie/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/crab_WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8.root"
if( whichSample == "sT" ):
  Find_str.append("find /fdata/hepx/store/user/lpernie/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/crab_ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/170405_172553 | grep root | grep -v failed > HADD/sT_STt-channel.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/crab_ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/170405_172612 | grep root | grep -v failed > HADD/sT_STt-channel.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/crab_ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/170405_172627 | grep root | grep -v failed > HADD/sT_STs-channel.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/crab_ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/170405_172642 | grep root | grep -v failed > HADD/sT_STtW.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/crab_ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/170405_172656 | grep root | grep -v failed > HADD/sT_STtW.txt")
  this_cat      = "cat HADD/sT_* > HADD/sT.txt"
  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/crab_ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1.root @HADD/sT.txt"
  this_NtotPath = "/fdata/hepx/store/user/lpernie/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/crab_ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1.root"
if( whichSample == "Wjet" ):
  Find_str.append("find /fdata/hepx/store/user/lpernie/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/170405_172710 | grep root | grep -v failed > HADD/Wjet_WJetsToLNuTuneCUETP8M1.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/170405_172724 | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-100To200.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/170405_172740 | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-200To400.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/170405_172757 | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-400To600.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/170405_172813 | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-600To800.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/170405_172829 | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-800To1200.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/170405_172846 | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-1200To2500.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/170405_172904 | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-2500ToInf.txt")
  this_cat      = "cat HADD/Wjet_* > HADD/Wjet.txt"
  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root @HADD/Wjet.txt"
  this_NtotPath = "/fdata/hepx/store/user/lpernie/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"
if( whichSample == "ttV" ):
  Find_str.append("find /fdata/hepx/store/user/lpernie/TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/170405_172925 | grep root | grep -v failed > HADD/ttV_TTWJetsToQQTuneCUETP8M1.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/170405_172942 | grep root | grep -v failed > HADD/ttV_TTWJetsToLNuTuneCUETP8M1.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/170405_173009 | grep root | grep -v failed > HADD/ttV_TTZToQQTuneCUETP8M1.txt")
  this_cat      = "cat HADD/ttV_* > HADD/ttV.txt"
  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root @HADD/ttV.txt"
  this_NtotPath = "/fdata/hepx/store/user/lpernie/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root"
  Find_str.append("find /fdata/hepx/store/user/lpernie/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/170405_173034 | grep root | grep -v failed > HADD/ttV_TTZToLLNuNuM-10.txt")
# Data
if( whichSample == "Data" ):
  Find_str.append("find /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Hhh_Run2016B-23Sep2016-v3/170405_170709 | grep root | grep -v failed > HADD/DATA_Hhh_Run2016B-23Sep2016-v3.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Hhh_Run2016C-23Sep2016-v1/170405_170722 | grep root | grep -v failed > HADD/DATA_Hhh_Run2016C-23Sep2016-v1.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Hhh_Run2016D-23Sep2016-v1/170405_170736 | grep root | grep -v failed > HADD/DATA_Hhh_Run2016D-23Sep2016-v1.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Hhh_Run2016E-23Sep2016-v1/170405_031502 | grep root | grep -v failed > HADD/DATA_Hhh_Run2016E-23Sep2016-v1.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Hhh_Run2016F-23Sep2016-v1/170405_170803 | grep root | grep -v failed > HADD/DATA_Hhh_Run2016F-23Sep2016-v1.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Hhh_Run2016G-23Sep2016-v1/170405_170817 | grep root | grep -v failed > HADD/DATA_Hhh_Run2016G-23Sep2016-v1.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Hhh_Run2016H-03Feb2017_ver2-v1/170405_170831 | grep root | grep -v failed > HADD/DATA_Hhh_Run2016H-03Feb2017_ver2-v1.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Hhh_Run2016H-03Feb2017_ver3-v1/170405_170845 | grep root | grep -v failed > HADD/DATA_Hhh_Run2016H-03Feb2017_ver3-v1.txt")
  this_cat      = "cat HADD/DATA_* > HADD/DATA.txt"

# Running
TCha = ROOT.TChain(tree_name)
for this_find in Find_str:
  os.system(this_find)
os.system(this_cat)
if (makeHadd=="HaddYes" and whichSample!="Data"): os.system(this_hadd)
Ntot_path = this_NtotPath
MyFile =  ROOT.TFile.Open(Ntot_path,"read");
h_prehlt  =  ROOT.TH1F(MyFile.Get("TriggerResults/hevent_filter")); nTOT_prehlt = h_prehlt.GetBinContent(2)
h_posthlt =  ROOT.TH1F(MyFile.Get("DiHiggsWWBBAna/hevent")); nTOT_posthlt = h_posthlt.GetBinContent(2);
print "HADD/TT.txt should be", this_cat.split(">")[1] 
with open(this_cat.split(">")[1],"r") as f:
  for line in f:
    if not line.isspace():
      TCha.Add(str(line[:-1]))
print whichSample, "TChain has", TCha.GetEntries(), "entries."
f = ROOT.TFile.Open('/fdata/hepx/store/user/lpernie/Hhh_For_Plotting/' + whichSample + '.root','w'); f.cd()

# Weights
h_Nev_preHLT  = ROOT.TH1F("h_Nev_preHLT","1,-0.5,0.5"); h_Nev_preHLT.GetXaxis().SetTitle("#Events (pre HLT)"); h_Nev_preHLT.SetBinContent(1,nTOT_prehlt)
h_Nev_posHLT  = ROOT.TH1F("h_Nev_posHLT","1,-0.5,0.5"); h_Nev_posHLT.GetXaxis().SetTitle("#Events (post HLT)"); h_Nev_posHLT.SetBinContent(1,nTOT_posthlt)
# Muons Reco
h_MU1_pt      = ROOT.TH1F("h_MU1_pt","",50,10.,350.); h_MU1_pt.GetXaxis().SetTitle("Lead. #mu P_{T} [GeV]"); h_MU1_pt.GetYaxis().SetTitle("Entries");
h_MU2_pt      = ROOT.TH1F("h_MU2_pt","",50,10.,350.); h_MU2_pt.GetXaxis().SetTitle("Sublead. #mu P_{T} [GeV]"); h_MU2_pt.GetYaxis().SetTitle("Entries");
h_MU1_eta     = ROOT.TH1F("h_MU1_eta","",50,-3.,3.); h_MU1_eta.GetXaxis().SetTitle("Lead. #mu #eta"); h_MU1_eta.GetYaxis().SetTitle("Entries");
h_MU2_eta     = ROOT.TH1F("h_MU2_eta","",50,-3.,3.); h_MU2_eta.GetXaxis().SetTitle("Sublead. #mu #eta"); h_MU2_eta.GetYaxis().SetTitle("Entries");
h_mass_l1l2   = ROOT.TH1F("h_mass_l1l2","",40,20.,150.); mass_l1l2.GetXaxis().SetTitle("m(l,l)"); h_mass_l1l2.GetYaxis().SetTitle("Entries");
h_dR_l1l2     = ROOT.TH1F("h_dR_l1l2","",50,0.,5.); dR_l1l2.GetXaxis().SetTitle("#Delta R(l,l)"); h_dR_l1l2.GetYaxis().SetTitle("Entries");
# B-jets Reco
h_J1_pt       = ROOT.TH1F("h_J1_pt","",50,20.,350.); b1jet_pt.GetXaxis().SetTitle("Lead. jet P_{T} [GeV]"); h_J1_pt.GetYaxis().SetTitle("Entries");
h_J2_pt       = ROOT.TH1F("h_J2_pt","",50,20.,350.); b2jet_pt.GetXaxis().SetTitle("Sublead. jet P_{T} [GeV]"); h_J2_pt.GetYaxis().SetTitle("Entries");
h_J1_eta      = ROOT.TH1F("h_J1_eta","",50,-3.,3.); b1jet_eta.GetXaxis().SetTitle("Lead. jet #eta"); h_J1_eta.GetYaxis().SetTitle("Entries");
h_J2_eta      = ROOT.TH1F("h_J2_eta","",50,-3.,3.); b2jet_eta.GetXaxis().SetTitle("Sublead. jet #eta"); h_J2_eta.GetYaxis().SetTitle("Entries");
h_mass_b1b2   = ROOT.TH1F("h_mass_b1b2","",40,40.,400.); mass_b1b2.GetXaxis().SetTitle("m(j,j)"); h_mass_b1b2.GetYaxis().SetTitle("Entries");
h_dR_b1b2     = ROOT.TH1F("h_dR_b1b2","",50,0.,6.); dR_b1b2.GetXaxis().SetTitle("#Delta R(j,j)"); h_dR_b1b2.GetYaxis().SetTitle("Entries");
# Mix Reco
h_met_pt      = ROOT.TH1F("h_met_pt","",50,10.,400.); met_pt.GetXaxis().SetTitle("MET [GeV]"); h_met_pt.GetYaxis().SetTitle("Entries");
h_mass_trans  = ROOT.TH1F("h_mass_trans","",50,0.,250.); mass_trans.GetXaxis().SetTitle("M_{trans} [GeV]"); h_mass_trans.GetYaxis().SetTitle("Entries");
h_dR_l1l2b1b2 = ROOT.TH1F("h_dR_l1l2b1b2","",50,0.,6.); dR_l1l2b1b2.GetXaxis().SetTitle("#Delta R(lljj)"); h_dR_l1l2b1b2.GetYaxis().SetTitle("Entries");
h_dphi_llmet  = ROOT.TH1F("h_dphi_llmet","",50,-4.,4.); dphi_llmet.GetXaxis().SetTitle("#Delta #phi (ll,MET)"); h_dphi_llmet.GetYaxis().SetTitle("Entries");

for ev in  TCha:
  MET_cut   = (ev.met_pt>20)
  MuMu_cut  = (ev.muon1_pt>20 and fabs(ev.muon1_eta)<2.4 and ev.muon2_pt>10 and fabs(ev.muon2_eta)<2.4 and ev.mass_l1l2>12)
  B1B2_cut  = (ev.b1jet_pt>20 and fabs(ev.b1jet_eta)<2.4 and ev.b2jet_pt>20 and fabs(ev.b2jet_eta)<2.4)
  preselection = (MET_cut and MuMu_cut and B1B2_cut)
  NoZ_cut   = (fabs(ev.mass_l1l2-91)>15)
  # Minimal Selection
  if( preselection ):
    h_MU1_pt.Fill( ev.muon1_pt )
    h_MU2_pt.Fill( ev.muon2_pt )
    h_MU2_eta.Fill( ev.muon2_eta )
    h_MU1_eta.Fill( ev.muon1_eta )
    h_mass_l1l2.Fill( ev.mass_l1l2 )
    h_dR_l1l2.Fill( ev.dR_l1l2 )
    h_J1_pt.Fill( ev.b1jet_pt )
    h_J2_pt.Fill( ev.b2jet_pt )
    h_J1_eta.Fill( ev.b1jet_eta )
    h_J2_eta.Fill( ev.b2jet_eta )
    h_mass_b1b2.Fill( ev.mass_b1b2 )
    h_dR_b1b2.Fill( ev.dR_b1b2 )
#    if( NoZ_cut ):

# Writing histograms
f.cd()
h_Nev_preHLT.Write()
h_Nev_posHLT.Write()
h_MU1_pt.Write()
h_MU2_pt.Write()
h_MU1_eta.Write()
h_MU2_eta.Write()
h_mass_l1l2.Write()
h_dR_l1l2.Write()
h_J1_pt.Write()
h_J2_pt.Write()
h_J1_eta.Write()
h_J2_eta.Write()
h_mass_b1b2.Write()
h_dR_b1b2.Write()
h_met_pt.Write()
h_mass_trans.Write()
h_dR_l1l2b1b2.Write()
h_dphi_llmet.Write()
f.Close()
