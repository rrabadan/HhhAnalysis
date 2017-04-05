import random
import ROOT
from ROOT import TFile,TChain,TH1F,TH2F,TLegend
import os
execfile("start.py")
execfile("functions.py")
#Creating folders and parameters
tree_name="DiHiggsWWBBAna/evtree"
os.system("rm -rf HADD/*txt")
print "Executing: python", sys.arg[0] , "-b", sys.arg[1], sys.arg[2], "(Arg1=makeHadd: True or False /|\ Arg2=whichSample: TT, DY, VV, singTop, Wjet, ttV, Data)"
makeHadd = sys.arg[1]
whichSample = sys.arg[2]

Find_str  = []; this_cat  = ""; this_hadd = ""; this_NtotPath = ""
# MC
if(whichSample=="TT"):
  Find_str.append("find /fdata/hepx/store/user/lpernie/TTTo2L2Nu_13TeV-powheg/crab_TTTo2L2Nu_13TeV-powheg/170401_035903 | grep root | grep -v failed > HADD/TT_TTTo2L2Nu13TeV-powheg.txt")
  this_cat  = "cat HADD/TT_* > HADD/TT.txt"
  this_hadd = "hadd -T -f -k /fdata/hepx/store/user/lpernie/TTTo2L2Nu_13TeV-powheg/crab_TTTo2L2Nu_13TeV-powheg.root @HADD/TT.txt"
  this_NtotPath = "/fdata/hepx/store/user/lpernie/TTTo2L2Nu_13TeV-powheg/crab_TTTo2L2Nu_13TeV-powheg.root"
if(whichSample=="DY"):
  Find_str.append("find /fdata/hepx/store/user/lpernie/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/170401_035916 | grep root | grep -v failed > HADD/DY_DYJetsToLLM-10to50.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/DYToLL_0J_13TeV-amcatnloFXFX-pythia8/crab_DYToLL_0J_13TeV-amcatnloFXFX-pythia8/170401_035928 | grep root | grep -v failed > HADD/DY_DYToLL0J.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/DYToLL_1J_13TeV-amcatnloFXFX-pythia8/crab_DYToLL_1J_13TeV-amcatnloFXFX-pythia8/170401_035940 | grep root | grep -v failed > HADD/DY_DYToLL1J.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/DYToLL_2J_13TeV-amcatnloFXFX-pythia8/crab_DYToLL_2J_13TeV-amcatnloFXFX-pythia8/170401_035952 | grep root | grep -v failed > HADD/DY_DYToLL2J.txt")
  this_cat = "cat HADD/DY_* > HADD/DY.txt")
  this_hadd = "hadd -T -f -k /fdata/hepx/store/user/lpernie/DYToLL_2J_13TeV-amcatnloFXFX-pythia8/crab_DYToLL_2J_13TeV-amcatnloFXFX-pythia8.root @HADD/DY.txt"
  this_NtotPath = "/fdata/hepx/store/user/lpernie/DYToLL_2J_13TeV-amcatnloFXFX-pythia8/crab_DYToLL_2J_13TeV-amcatnloFXFX-pythia8.root"

# Running
TCha = ROOT.TChain(tree_name)
for this_find in Find_str:
  os.system(this_find)
os.system(this_cat)
if makeHadd: os.system(this_hadd)
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
h_MU1_pt      = ROOT.TH1F("h_MU1_pt","",50,10.,350.); h_MU1_pt.GetXaxis().SetTitle("Lead. #mu P_{T} [GeV]"); .GetYaxis().SetTitle("Entries");
h_MU2_pt      = ROOT.TH1F("h_MU2_pt","",50,10.,350.); h_MU2_pt.GetXaxis().SetTitle("Sublead. #mu P_{T} [GeV]"); .GetYaxis().SetTitle("Entries");
h_MU1_eta     = ROOT.TH1F("h_MU1_eta","",50,-3.,3.); h_MU1_eta.GetXaxis().SetTitle("Lead. #mu #eta"); .GetYaxis().SetTitle("Entries");
h_MU2_eta     = ROOT.TH1F("h_MU2_eta","",50,-3.,3.); h_MU2_eta.GetXaxis().SetTitle("Sublead. #mu #eta"); .GetYaxis().SetTitle("Entries");
h_mass_l1l2   = ROOT.TH1F("h_mass_l1l2","",40,20.,150.); mass_l1l2.GetXaxis().SetTitle("m(l,l)"); .GetYaxis().SetTitle("Entries");
h_dR_l1l2     = ROOT.TH1F("h_dR_l1l2","",50,0.,5.); dR_l1l2.GetXaxis().SetTitle("#Delta R(l,l)"); .GetYaxis().SetTitle("Entries");
# B-jets Reco
h_J1_pt       = ROOT.TH1F("h_J1_pt","",50,20.,350.); b1jet_pt.GetXaxis().SetTitle("Lead. jet P_{T} [GeV]"); .GetYaxis().SetTitle("Entries");
h_J2_pt       = ROOT.TH1F("h_J2_pt","",50,20.,350.); b2jet_pt.GetXaxis().SetTitle("Sublead. jet P_{T} [GeV]"); .GetYaxis().SetTitle("Entries");
h_J1_eta      = ROOT.TH1F("h_J1_eta","",50,-3.,3.); b1jet_eta.GetXaxis().SetTitle("Lead. jet #eta"); .GetYaxis().SetTitle("Entries");
h_J2_eta      = ROOT.TH1F("h_J2_eta","",50,-3.,3.); b2jet_eta.GetXaxis().SetTitle("Sublead. jet #eta"); .GetYaxis().SetTitle("Entries");
h_mass_b1b2   = ROOT.TH1F("h_mass_b1b2","",40,40.,400.); mass_b1b2.GetXaxis().SetTitle("m(j,j)"); .GetYaxis().SetTitle("Entries");
h_dR_b1b2     = ROOT.TH1F("h_dR_b1b2","",50,0.,6.); dR_b1b2.GetXaxis().SetTitle("#Delta R(j,j)"); .GetYaxis().SetTitle("Entries");
# Mix Reco
h_met_pt      = ROOT.TH1F("h_met_pt","",50,10.,400.); met_pt.GetXaxis().SetTitle("MET [GeV]"); .GetYaxis().SetTitle("Entries");
h_mass_trans  = ROOT.TH1F("h_mass_trans","",50,0.,250.); mass_trans.GetXaxis().SetTitle("M_{trans} [GeV]"); .GetYaxis().SetTitle("Entries");
h_dR_l1l2b1b2 = ROOT.TH1F("h_dR_l1l2b1b2","",50,0.,6.); dR_l1l2b1b2.GetXaxis().SetTitle("#Delta R(lljj)"); .GetYaxis().SetTitle("Entries");
h_dphi_llmet  = ROOT.TH1F("h_dphi_llmet","",50,-4.,4.); dphi_llmet.GetXaxis().SetTitle("#Delta #phi (ll,MET)"); .GetYaxis().SetTitle("Entries");

for ev in TT_ch:
  MET_cut   = (ev.met_pt>20)
  MuMu_cut  = (ev.muon1_pt>20 && TMath::Abs(ev.muon1_eta)<2.4 && ev.muon2_pt>10 && TMath::Abs(ev.muon2_eta)<2.4 && ev.mass_l1l2>12)
  B1B2_cut  = (ev.b1jet_pt>20 && TMath::Abs(ev.b1jet_eta)<2.4 && ev.b2jet_pt>20 && TMath::Abs(ev.b2jet_eta)<2.4)
  NoZ_cut   = (fabs(ev.mass_l1l2-91)>15)
  # Minimal Selection
  if( MET_cut && MuMu_cut && B1B2_cut ):
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
    if( NoZ_cut ):


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
