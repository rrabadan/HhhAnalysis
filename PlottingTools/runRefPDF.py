import os, sys, random
import ROOT
from ROOT import TFile,TChain,TH1F,TH2F,TLegend
from math import *
from HeavyMassEstimator import *
import argparse
import numpy as np

import Samplelist 
execfile("start.py")
execfile("functions.py")
#Creating folders and parameters

doTest = False
doOnlyRefPDF = True
refPDF = ROOT.TFile("REFPDFPU40.root","READ")
onshellWmasspdf = refPDF.Get("onshellWmasspdf")
onshellnuptpdf = refPDF.Get("onshellnuptpdf")
recobjetrescalec1pdfPU40 = refPDF.Get("recobjetrescalec1pdfPU40")

tree_name="DiHiggsWWBBAna/evtree"
benchmarks = ["ttV", "Wjet", "singTop", "VV", "DY", "TTbar", "Data"]
for mass in [260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 900, 1000]:
    benchmarks.append("graviton_M%d"%mass)
    benchmarks.append("radion_M%d"%mass)

#whichSample = "Graviton_0419"

parser = argparse.ArgumentParser(description='runHME')
parser.add_argument("-n", "--njobs", dest="njobs", type=int, default=100, help="total splitted jobs. [Default: 100]")
parser.add_argument("-i", "--ijob", dest="ijob", type=int, default=0, help="ith job in splitting. [Default: 0]")
parser.add_argument("-jt", "--jobtype", dest="jobtype", default="radion_M260", help="sample type to run. [Default: radion]")
args = parser.parse_args()
whichSample = args.jobtype
makeHadd = "HaddYes"

"""
print "Executing: python", sys.argv[0] , "-b", sys.argv[2], sys.argv[3], "(Arg1=makeHadd: HaddYes or HaddNo /|\ Arg2=whichSample: TT, DY, VV, sT, Wjet, ttV, Data)"
makeHadd = sys.argv[2]
whichSample = sys.argv[3]
if( makeHadd!="HaddYes" and makeHadd!="HaddNo" ): print "WARNING! 1st argument has to be HaddYes or HaddNo"; sys.exit()
if( whichSample!="TT" and whichSample!="DY" and whichSample!="VV" and whichSample!="sT" and whichSample!="Wjet" and whichSample!="ttV" and whichSample!="Data" ):  print "WARNING! 2nd argument have to be TT, DY, VV, singTop, Wjet, ttV, or Data"; sys.exit()

"""
######################################
## define shell commands
######################################
Find_str = []; this_cat = ""; this_hadd = ""; this_NtotPath = ""
for i in range(len(Samplelist.outAnalist[whichSample])):
  Find_str.append("find %s | grep root | grep -v failed > HADD/%s_%d.txt"%(Samplelist.outAnalist[whichSample][i], Samplelist.sampleFullName[whichSample], i))
    #Find_str.append("find /fdata/hepx/store/user/lpernie/TTTo2L2Nu_13TeV-powheg/crab_TTTo2L2Nu_13TeV-powheg/170405_172215 | grep root | grep -v failed > HADD/TT_TTTo2L2Nu13TeV-powheg.txt")
this_cat = "cat HADD/%s_* > HADD/%s.txt"%(Samplelist.sampleFullName[whichSample], whichSample)
this_NtotDir  = "/fdata/hepx/store/user/%s/%s/"%(user,Samplelist.sampleFullName[whichSample])
if not(os.path.exists(this_NtotDir)):
    os.system("mkdir %s"%this_NtotDir)
if whichSample != "DATA":
    this_NtotPath = this_NtotDir+"crab_%s.root"%(Samplelist.sampleFullName[whichSample])
    this_hadd     = "hadd -T -f -k %s @HADD/%s.txt"%(this_NtotPath, whichSample)

######################################
# Running commands
######################################
for this_find in Find_str:
  os.system(this_find)
os.system(this_cat)
if (makeHadd=="HaddYes" and whichSample!="Data"): os.system(this_hadd)
######################################
#this_NtotPath = "/home/taohuang/DiHiggsAnalysis/CMSSW_8_0_26_patch2/src/HhhAnalysis/CutFlowAnalyzer/test/out_ana_10k_graviton.root"
#this_NtotPath = "/fdata/hepx/store/user/taohuang/DiHiggsAnalysisSample/out_ann_%s_20160411.root"%whichSample
if( whichSample != "Data" ):
  Ntot_path = this_NtotPath
  MyFile =  ROOT.TFile.Open(Ntot_path,"read");
  h_prehlt  =  ROOT.TH1F(MyFile.Get("TriggerResults/hevent_filter")); nTOT_prehlt = h_prehlt.GetBinContent(2)
  h_posthlt =  ROOT.TH1F(MyFile.Get("DiHiggsWWBBAna/hevent")); nTOT_posthlt = h_posthlt.GetBinContent(2);
else:
  nTOT_prehlt = 0.; nTOT_posthlt = 0.

TCha = ROOT.TChain(tree_name)
with open(this_cat.split(">")[1].split(" ")[1],"r") as f:
  for line in f:
    if not line.isspace():
      TCha.Add(str(line[:-1]))
#TCha.Add("/home/taohuang/DiHiggsAnalysis/CMSSW_8_0_26_patch2/src/HhhAnalysis/CutFlowAnalyzer/test/out_ana_10k_graviton.root")
#TCha.Add(this_NtotPath)
print whichSample, "TChain has", TCha.GetEntries(), "entries."
TotalEv = TCha.GetEntries()
nStart = TotalEv/args.njobs*args.ijob
nEnd = TotalEv/args.njobs*(args.ijob+1)
if args.ijob == args.njobs-1:
    nEnd = -1

histdir = '/fdata/hepx/store/user/%s/Hhh_For_Plotting/%s_REFPDF/'%(user, whichSample)
if not(os.path.exists(histdir)):
    os.system("mkdir %s"%histdir)
f = ROOT.TFile(histdir + whichSample + '_ijob%d.root'%args.ijob,'recreate'); f.cd()

ptbinSF= [20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 120.0, 200.0, 300.]
etabinSF = [0, 0.9, 1.2, 2.1, 2.4]
myptbinSF = np.asarray(ptbinSF)
myetabinSF = np.asarray(etabinSF)

# PRESELECTION
# Weights
h_pre_Nev_preHLT         = ROOT.TH1F("h_pre_Nev_preHLT","",1,-0.5,0.5);          h_pre_Nev_preHLT.GetXaxis().SetTitle("#Events (pre HLT)"); h_pre_Nev_preHLT.SetBinContent(1,nTOT_prehlt)
h_pre_Nev_posHLT         = ROOT.TH1F("h_pre_Nev_posHLT","",1,-0.5,0.5);          h_pre_Nev_posHLT.GetXaxis().SetTitle("#Events (post HLT)"); h_pre_Nev_posHLT.SetBinContent(1,nTOT_posthlt)
h_pre_XsecBr             = ROOT.TH1F("h_pre_XsecBr","",1,-0.5,0.5);              h_pre_XsecBr.GetXaxis().SetTitle("XsecBr");
h_pre_muon1_triggerSF    = ROOT.TH2F("h_pre_muon1_triggerSF","",len(etabinSF)-1,myetabinSF,len(ptbinSF)-1,myptbinSF);  h_pre_muon1_triggerSF.GetXaxis().SetTitle("#mu 1 Trigger SF");
h_pre_muon1_isoSF        = ROOT.TH2F("h_pre_muon1_isoSF","",len(etabinSF)-1,myetabinSF,len(ptbinSF)-1,myptbinSF); h_pre_muon1_isoSF.GetXaxis().SetTitle("#mu 1 Iso SF");
h_pre_muon1_idSF         = ROOT.TH2F("h_pre_muon1_idSF","",len(etabinSF)-1,myetabinSF,len(ptbinSF)-1,myptbinSF);  h_pre_muon1_idSF.GetXaxis().SetTitle("#mu 1 ID SF");
h_pre_muon1_trackingSF   = ROOT.TH1F("h_pre_muon1_trackingSF","",50,-2.4,2.4); h_pre_muon1_trackingSF.GetXaxis().SetTitle("#mu 1 tracking SF");
h_pre_muon1_SF_bg2       = ROOT.TH2F("h_pre_muon1_SF_bg2","",len(etabinSF)-1,myetabinSF,len(ptbinSF)-1,myptbinSF);  h_pre_muon1_SF_bg2.GetXaxis().SetTitle("#mu SF bg2");
h_pre_muon1_SF_bg1       = ROOT.TH1F("h_pre_moun1_SF_bg1","",50,-2.4,2.4); h_pre_muon1_SF_bg1.GetXaxis().SetTitle("#mu SF bg1");
h_pre_muon2_triggerSF    = ROOT.TH2F("h_pre_muon2_triggerSF","",len(etabinSF)-1,myetabinSF,len(ptbinSF)-1,myptbinSF);  h_pre_muon2_triggerSF.GetXaxis().SetTitle("#mu 1 Trigger SF");
h_pre_muon2_isoSF        = ROOT.TH2F("h_pre_muon2_isoSF","",len(etabinSF)-1,myetabinSF,len(ptbinSF)-1,myptbinSF);      h_pre_muon2_isoSF.GetXaxis().SetTitle("#mu 1 Iso SF");
h_pre_muon2_idSF         = ROOT.TH2F("h_pre_muon2_idSF","",len(etabinSF)-1,myetabinSF, len(ptbinSF)-1, myptbinSF);       h_pre_muon2_idSF.GetXaxis().SetTitle("#mu 1 ID SF");
h_pre_muon2_trackingSF   = ROOT.TH1F("h_pre_muon2_trackingSF","",50,-2.4,2.4); h_pre_muon2_trackingSF.GetXaxis().SetTitle("#mu 1 tracking SF");
h_pre_muon2_SF_bg2       = ROOT.TH2F("h_pre_muon2_SF_bg2","",len(etabinSF)-1,myetabinSF,len(ptbinSF)-1,myptbinSF);  h_pre_muon2_SF_bg2.GetXaxis().SetTitle("#mu SF bg2");
h_pre_muon2_SF_bg1       = ROOT.TH1F("h_pre_moun2_SF_bg1","",50,-2.4,2.4); h_pre_muon1_SF_bg1.GetXaxis().SetTitle("#mu SF bg1");
# Regression variables
h_pre_numOfVertices      = ROOT.TH1F("h_pre_numOfVertices","",50,0.,50.);      h_pre_numOfVertices.GetXaxis().SetTitle("");
h_pre_b1jet_mt           = ROOT.TH1F("h_pre_b1jet_mt","",50,0.,100.);          h_pre_b1jet_mt.GetXaxis().SetTitle("");
h_pre_b1jet_leadTrackPt  = ROOT.TH1F("h_pre_b1jet_leadTrackPt","",50,0.,100.); h_pre_b1jet_leadTrackPt.GetXaxis().SetTitle("");
h_pre_b1jet_leptonPtRel  = ROOT.TH1F("h_pre_b1jet_leptonPtRel","",50,0.,10.);  h_pre_b1jet_leptonPtRel.GetXaxis().SetTitle("");
h_pre_b1jet_leptonPt     = ROOT.TH1F("h_pre_b1jet_leptonPt","",50,0.,100.);    h_pre_b1jet_leptonPt.GetXaxis().SetTitle("");
h_pre_b1jet_leptonDeltaR = ROOT.TH1F("h_pre_b1jet_leptonDeltaR","",50,0.,5.);  h_pre_b1jet_leptonDeltaR.GetXaxis().SetTitle("");
h_pre_b1jet_neHEF        = ROOT.TH1F("h_pre_b1jet_neHEF","",50,0.,2.);         h_pre_b1jet_neHEF.GetXaxis().SetTitle("");
h_pre_b1jet_neEmEF       = ROOT.TH1F("h_pre_b1jet_neEmEF","",50,0.,2.);        h_pre_b1jet_neEmEF.GetXaxis().SetTitle("");
h_pre_b1jet_vtxNtracks   = ROOT.TH1F("h_pre_b1jet_vtxNtracks","",50,0.,50.);   h_pre_b1jet_vtxNtracks.GetXaxis().SetTitle("");
h_pre_b1jet_vtxPt        = ROOT.TH1F("h_pre_b1jet_vtxPt","",50,0.,100.);       h_pre_b1jet_vtxPt.GetXaxis().SetTitle("");
h_pre_b1jet_vtxMass      = ROOT.TH1F("h_pre_b1jet_vtxMass","",50,0.,50.);      h_pre_b1jet_vtxMass.GetXaxis().SetTitle("");
h_pre_b1jet_vtx3DSig     = ROOT.TH1F("h_pre_b1jet_vtx3DSig","",50,0.,5.);      h_pre_b1jet_vtx3DSig.GetXaxis().SetTitle("");
h_pre_b1jet_vtx3DVal     = ROOT.TH1F("h_pre_b1jet_vtx3DVal","",50,0.,5.);      h_pre_b1jet_vtx3DVal.GetXaxis().SetTitle("");

# Muons gen
h_gen_pre_MU1_pt             = ROOT.TH1F("h_gen_pre_MU1_pt","",50,10.,350.);    h_gen_pre_MU1_pt.GetXaxis().SetTitle("Lead. #mu P_{T} [GeV]");    
h_gen_pre_MU2_pt             = ROOT.TH1F("h_gen_pre_MU2_pt","",50,10.,350.);    h_gen_pre_MU2_pt.GetXaxis().SetTitle("Sublead. #mu P_{T} [GeV]"); 
h_gen_pre_MU1_eta            = ROOT.TH1F("h_gen_pre_MU1_eta","",50,-3.,3.);     h_gen_pre_MU1_eta.GetXaxis().SetTitle("Lead. #mu #eta");          
h_gen_pre_MU2_eta            = ROOT.TH1F("h_gen_pre_MU2_eta","",50,-3.,3.);     h_gen_pre_MU2_eta.GetXaxis().SetTitle("Sublead. #mu #eta");       
h_gen_pre_mass_l1l2          = ROOT.TH1F("h_gen_pre_mass_l1l2","",40,20.,150.); h_gen_pre_mass_l1l2.GetXaxis().SetTitle("m(l,l)");               
h_gen_pre_dR_l1l2            = ROOT.TH1F("h_gen_pre_dR_l1l2","",50,0.,5.);      h_gen_pre_dR_l1l2.GetXaxis().SetTitle("#Delta R(l,l)");           
# B-jets Reco
h_gen_pre_J1_pt              = ROOT.TH1F("h_gen_pre_J1_pt","",50,20.,350.);     h_gen_pre_J1_pt.GetXaxis().SetTitle("Lead. jet P_{T} [GeV]");     
h_gen_pre_J2_pt              = ROOT.TH1F("h_gen_pre_J2_pt","",50,20.,350.);     h_gen_pre_J2_pt.GetXaxis().SetTitle("Sublead. jet P_{T} [GeV]");  
h_gen_pre_J1_eta             = ROOT.TH1F("h_gen_pre_J1_eta","",50,-3.,3.);      h_gen_pre_J1_eta.GetXaxis().SetTitle("Lead. jet #eta");           
h_gen_pre_J2_eta             = ROOT.TH1F("h_gen_pre_J2_eta","",50,-3.,3.);      h_gen_pre_J2_eta.GetXaxis().SetTitle("Sublead. jet #eta");        
h_gen_pre_mass_b1b2          = ROOT.TH1F("h_gen_pre_mass_b1b2","",40,40.,400.); h_gen_pre_mass_b1b2.GetXaxis().SetTitle("m(j,j)");      
h_gen_pre_dR_b1b2            = ROOT.TH1F("h_gen_pre_dR_b1b2","",50,0.,6.);      h_gen_pre_dR_b1b2.GetXaxis().SetTitle("#Delta R(j,j)");           
# Mix Reco
h_gen_pre_met_pt             = ROOT.TH1F("h_gen_pre_met_pt","",50,10.,400.);    h_gen_pre_met_pt.GetXaxis().SetTitle("MET [GeV]");                
h_gen_pre_mass_trans         = ROOT.TH1F("h_gen_pre_mass_trans","",50,0.,250.); h_gen_pre_mass_trans.GetXaxis().SetTitle("M_{trans} [GeV]");     
h_gen_pre_dR_l1l2b1b2        = ROOT.TH1F("h_gen_pre_dR_l1l2b1b2","",50,0.,6.);  h_gen_pre_dR_l1l2b1b2.GetXaxis().SetTitle("#Delta R(lljj)");     
h_gen_pre_dphi_llmet         = ROOT.TH1F("h_gen_pre_dphi_llmet","",50,-4.,4.);  h_gen_pre_dphi_llmet.GetXaxis().SetTitle("#Delta #phi (ll,MET)");
h_gen_pre_onshellWmass       = ROOT.TH1F("h_gen_pre_onshellWmass","",70,30.,100.);  h_gen_pre_onshellWmass.GetXaxis().SetTitle("W^{onshell} mass [GeV]");
h_gen_pre_offshellWmass       = ROOT.TH1F("h_gen_pre_offshellWmass","",70,30.,100.);  h_gen_pre_offshellWmass.GetXaxis().SetTitle("W^{onshell} mass [GeV]");
h_gen_pre_onshellnupt        = ROOT.TH1F("h_gen_pre_onshellnupt","",500,0.,500.);  h_gen_pre_onshellnupt.GetXaxis().SetTitle("p_{T} of #nu^{onshellW} [GeV]");
h_gen_pre_onoffshellWmass    = ROOT.TH2F("h_gen_pre_onoffshellWmass","",70,30.,100.,70,0.0,70.0);  h_gen_pre_onoffshellWmass.GetXaxis().SetTitle(" W^{onshell} mass [GeV]"); h_gen_pre_onoffshellWmass.GetYaxis().SetTitle(" W^{offshell} mass [GeV]");
h_gen_pre_bjetrescalec1dR4   = ROOT.TH1F("h_gen_pre_bjetrescalec1dR4","",600,0.,6);  h_gen_pre_bjetrescalec1dR4.GetXaxis().SetTitle("#frac{p_{T}(b)}{p_{T}(bjet)}, leading bjet");
h_gen_pre_bjetrescalec2dR4   = ROOT.TH1F("h_gen_pre_bjetrescalec2dR4","",600,0.,6);  h_gen_pre_bjetrescalec2dR4.GetXaxis().SetTitle("#frac{p_{T}(b)}{p_{T}(bjet)}, subleading bjet");

# Muons Reco
h_pre_MU1_pt             = ROOT.TH1F("h_pre_MU1_pt","",50,10.,350.);    h_pre_MU1_pt.GetXaxis().SetTitle("Lead. #mu P_{T} [GeV]");    
h_pre_MU2_pt             = ROOT.TH1F("h_pre_MU2_pt","",50,10.,350.);    h_pre_MU2_pt.GetXaxis().SetTitle("Sublead. #mu P_{T} [GeV]"); 
h_pre_MU1_eta            = ROOT.TH1F("h_pre_MU1_eta","",50,-3.,3.);     h_pre_MU1_eta.GetXaxis().SetTitle("Lead. #mu #eta");          
h_pre_MU2_eta            = ROOT.TH1F("h_pre_MU2_eta","",50,-3.,3.);     h_pre_MU2_eta.GetXaxis().SetTitle("Sublead. #mu #eta");       
h_pre_mass_l1l2          = ROOT.TH1F("h_pre_mass_l1l2","",40,20.,150.); h_pre_mass_l1l2.GetXaxis().SetTitle("m(l,l)");               
h_pre_dR_l1l2            = ROOT.TH1F("h_pre_dR_l1l2","",50,0.,5.);      h_pre_dR_l1l2.GetXaxis().SetTitle("#Delta R(l,l)");           
# B-jets Reco
h_pre_J1_pt              = ROOT.TH1F("h_pre_J1_pt","",50,20.,350.);     h_pre_J1_pt.GetXaxis().SetTitle("Lead. jet P_{T} [GeV]");     
h_pre_J2_pt              = ROOT.TH1F("h_pre_J2_pt","",50,20.,350.);     h_pre_J2_pt.GetXaxis().SetTitle("Sublead. jet P_{T} [GeV]");  
h_pre_J1_eta             = ROOT.TH1F("h_pre_J1_eta","",50,-3.,3.);      h_pre_J1_eta.GetXaxis().SetTitle("Lead. jet #eta");           
h_pre_J2_eta             = ROOT.TH1F("h_pre_J2_eta","",50,-3.,3.);      h_pre_J2_eta.GetXaxis().SetTitle("Sublead. jet #eta");        
h_pre_mass_b1b2          = ROOT.TH1F("h_pre_mass_b1b2","",40,40.,400.); h_pre_mass_b1b2.GetXaxis().SetTitle("m(j,j)");      
h_pre_dR_b1b2            = ROOT.TH1F("h_pre_dR_b1b2","",50,0.,6.);      h_pre_dR_b1b2.GetXaxis().SetTitle("#Delta R(j,j)");           
# Mix Reco
h_pre_met_pt             = ROOT.TH1F("h_pre_met_pt","",50,10.,400.);    h_pre_met_pt.GetXaxis().SetTitle("MET [GeV]");                
h_pre_mass_trans         = ROOT.TH1F("h_pre_mass_trans","",50,0.,250.); h_pre_mass_trans.GetXaxis().SetTitle("M_{trans} [GeV]");     
h_pre_dR_l1l2b1b2        = ROOT.TH1F("h_pre_dR_l1l2b1b2","",50,0.,6.);  h_pre_dR_l1l2b1b2.GetXaxis().SetTitle("#Delta R(lljj)");     
h_pre_dphi_llmet         = ROOT.TH1F("h_pre_dphi_llmet","",50,-4.,4.);  h_pre_dphi_llmet.GetXaxis().SetTitle("#Delta #phi (ll,MET)");
# CLEANING CUT
# Muons gen
h_gen_cc_MU1_pt              = ROOT.TH1F("h_gen_cc_MU1_pt","",50,10.,350.);    h_gen_cc_MU1_pt.GetXaxis().SetTitle("Lead. #mu P_{T} [GeV]");    
h_gen_cc_MU2_pt              = ROOT.TH1F("h_gen_cc_MU2_pt","",50,10.,350.);    h_gen_cc_MU2_pt.GetXaxis().SetTitle("Sublead. #mu P_{T} [GeV]"); 
h_gen_cc_MU1_eta             = ROOT.TH1F("h_gen_cc_MU1_eta","",50,-3.,3.);     h_gen_cc_MU1_eta.GetXaxis().SetTitle("Lead. #mu #eta");          
h_gen_cc_MU2_eta             = ROOT.TH1F("h_gen_cc_MU2_eta","",50,-3.,3.);     h_gen_cc_MU2_eta.GetXaxis().SetTitle("Sublead. #mu #eta");       
h_gen_cc_mass_l1l2           = ROOT.TH1F("h_gen_cc_mass_l1l2","",40,20.,150.); h_gen_cc_mass_l1l2.GetXaxis().SetTitle("m(l,l)");             
h_gen_cc_dR_l1l2             = ROOT.TH1F("h_gen_cc_dR_l1l2","",50,0.,5.);      h_gen_cc_dR_l1l2.GetXaxis().SetTitle("#Delta R(l,l)");           
# B-jets gen
h_gen_cc_J1_pt               = ROOT.TH1F("h_gen_cc_J1_pt","",50,20.,350.);     h_gen_cc_J1_pt.GetXaxis().SetTitle("Lead. jet P_{T} [GeV]");     
h_gen_cc_J2_pt               = ROOT.TH1F("h_gen_cc_J2_pt","",50,20.,350.);     h_gen_cc_J2_pt.GetXaxis().SetTitle("Sublead. jet P_{T} [GeV]");  
h_gen_cc_J1_eta              = ROOT.TH1F("h_gen_cc_J1_eta","",50,-3.,3.);      h_gen_cc_J1_eta.GetXaxis().SetTitle("Lead. jet #eta");           
h_gen_cc_J2_eta              = ROOT.TH1F("h_gen_cc_J2_eta","",50,-3.,3.);      h_gen_cc_J2_eta.GetXaxis().SetTitle("Sublead. jet #eta");        
h_gen_cc_mass_b1b2           = ROOT.TH1F("h_gen_cc_mass_b1b2","",40,40.,400.); h_gen_cc_mass_b1b2.GetXaxis().SetTitle("m(j,j)");              
h_gen_cc_dR_b1b2             = ROOT.TH1F("h_gen_cc_dR_b1b2","",50,0.,6.);      h_gen_cc_dR_b1b2.GetXaxis().SetTitle("#Delta R(j,j)");           
# Mix gen
h_gen_cc_met_pt              = ROOT.TH1F("h_gen_cc_met_pt","",50,10.,400.);    h_gen_cc_met_pt.GetXaxis().SetTitle("MET [GeV]");                
h_gen_cc_mass_trans          = ROOT.TH1F("h_gen_cc_mass_trans","",50,0.,250.); h_gen_cc_mass_trans.GetXaxis().SetTitle("M_{trans} [GeV]");   
h_gen_cc_dR_l1l2b1b2         = ROOT.TH1F("h_gen_cc_dR_l1l2b1b2","",50,0.,6.);  h_gen_cc_dR_l1l2b1b2.GetXaxis().SetTitle("#Delta R(lljj)");
h_gen_cc_dphi_llmet          = ROOT.TH1F("h_gen_cc_dphi_llmet","",50,-4.,4.);  h_gen_cc_dphi_llmet.GetXaxis().SetTitle("#Delta #phi (ll,MET)");
h_gen_cc_onshellWmass       = ROOT.TH1F("h_gen_cc_onshellWmass","",70,30.,100.);  h_gen_cc_onshellWmass.GetXaxis().SetTitle("W^{onshell} mass [GeV]");
h_gen_cc_offshellWmass       = ROOT.TH1F("h_gen_cc_offshellWmass","",70,30.,100.);  h_gen_cc_offshellWmass.GetXaxis().SetTitle("W^{onshell} mass [GeV]");
h_gen_cc_onshellnupt        = ROOT.TH1F("h_gen_cc_onshellnupt","",500,0.,500.);  h_gen_cc_onshellnupt.GetXaxis().SetTitle("p_{T} of #nu^{onshellW} [GeV]");
h_gen_cc_onoffshellWmass    = ROOT.TH2F("h_gen_cc_onoffshellWmass","",70,30.,100.,70,0.0,70.0);  h_gen_cc_onoffshellWmass.GetXaxis().SetTitle(" W^{onshell} mass [GeV]"); h_gen_cc_onoffshellWmass.GetYaxis().SetTitle(" W^{offshell} mass [GeV]");
h_gen_cc_bjetrescalec1dR4   = ROOT.TH1F("h_gen_cc_bjetrescalec1dR4","",600,0.,6);  h_gen_cc_bjetrescalec1dR4.GetXaxis().SetTitle("#frac{p_{T}(b)}{p_{T}(bjet)}, leading bjet");
h_gen_cc_bjetrescalec2dR4   = ROOT.TH1F("h_gen_cc_bjetrescalec2dR4","",600,0.,6);  h_gen_cc_bjetrescalec2dR4.GetXaxis().SetTitle("#frac{p_{T}(b)}{p_{T}(bjet)}, subleading bjet");

# Muons Reco
h_cc_MU1_pt              = ROOT.TH1F("h_cc_MU1_pt","",50,10.,350.);    h_cc_MU1_pt.GetXaxis().SetTitle("Lead. #mu P_{T} [GeV]");    
h_cc_MU2_pt              = ROOT.TH1F("h_cc_MU2_pt","",50,10.,350.);    h_cc_MU2_pt.GetXaxis().SetTitle("Sublead. #mu P_{T} [GeV]"); 
h_cc_MU1_eta             = ROOT.TH1F("h_cc_MU1_eta","",50,-3.,3.);     h_cc_MU1_eta.GetXaxis().SetTitle("Lead. #mu #eta");          
h_cc_MU2_eta             = ROOT.TH1F("h_cc_MU2_eta","",50,-3.,3.);     h_cc_MU2_eta.GetXaxis().SetTitle("Sublead. #mu #eta");       
h_cc_mass_l1l2           = ROOT.TH1F("h_cc_mass_l1l2","",40,20.,150.); h_cc_mass_l1l2.GetXaxis().SetTitle("m(l,l)");             
h_cc_dR_l1l2             = ROOT.TH1F("h_cc_dR_l1l2","",50,0.,5.);      h_cc_dR_l1l2.GetXaxis().SetTitle("#Delta R(l,l)");           
# B-jets Reco
h_cc_J1_pt               = ROOT.TH1F("h_cc_J1_pt","",50,20.,350.);     h_cc_J1_pt.GetXaxis().SetTitle("Lead. jet P_{T} [GeV]");     
h_cc_J2_pt               = ROOT.TH1F("h_cc_J2_pt","",50,20.,350.);     h_cc_J2_pt.GetXaxis().SetTitle("Sublead. jet P_{T} [GeV]");  
h_cc_J1_eta              = ROOT.TH1F("h_cc_J1_eta","",50,-3.,3.);      h_cc_J1_eta.GetXaxis().SetTitle("Lead. jet #eta");           
h_cc_J2_eta              = ROOT.TH1F("h_cc_J2_eta","",50,-3.,3.);      h_cc_J2_eta.GetXaxis().SetTitle("Sublead. jet #eta");        
h_cc_mass_b1b2           = ROOT.TH1F("h_cc_mass_b1b2","",40,40.,400.); h_cc_mass_b1b2.GetXaxis().SetTitle("m(j,j)");              
h_cc_dR_b1b2             = ROOT.TH1F("h_cc_dR_b1b2","",50,0.,6.);      h_cc_dR_b1b2.GetXaxis().SetTitle("#Delta R(j,j)");           
# Mix Reco
h_cc_met_pt              = ROOT.TH1F("h_cc_met_pt","",50,10.,400.);    h_cc_met_pt.GetXaxis().SetTitle("MET [GeV]");                
h_cc_mass_trans          = ROOT.TH1F("h_cc_mass_trans","",50,0.,250.); h_cc_mass_trans.GetXaxis().SetTitle("M_{trans} [GeV]");   
h_cc_dR_l1l2b1b2         = ROOT.TH1F("h_cc_dR_l1l2b1b2","",50,0.,6.);  h_cc_dR_l1l2b1b2.GetXaxis().SetTitle("#Delta R(lljj)");
h_cc_dphi_llmet          = ROOT.TH1F("h_cc_dphi_llmet","",50,-4.,4.);  h_cc_dphi_llmet.GetXaxis().SetTitle("#Delta #phi (ll,MET)");


nEv = 0
for ev in TCha:
  if nEv < nStart or (nEnd > 0 and nEv >= nEnd):
      #print "nEv ",nEv," nstart ",nStart," nEnd ",nEnd
      nEv = nEv +1
      continue
  if (doTest and nEv%10 == 0 ):
      print "nEv ",nEv
  elif (nEv%10000 == 0):
      print "nEv ",nEv
  
  if (doTest and nEv>=10000):
      break
  # CUTS
  MET_cut             = (ev.met_pt>20)
  MuMu_cut            = (ev.muon1_pt>20 and fabs(ev.muon1_eta)<2.4 and ev.muon2_pt>10 and fabs(ev.muon2_eta)<2.4 and ev.mass_l1l2>12)
  B1B2_cut            = (ev.b1jet_pt>20 and fabs(ev.b1jet_eta)<2.4 and ev.b2jet_pt>20 and fabs(ev.b2jet_eta)<2.4)
  preselection        = (MET_cut and MuMu_cut and B1B2_cut)
  mt_cleancut         = (ev.mass_trans>10)
  NoZ_cleancut        = (ev.mass_l1l2<(91-15))
  mJJ_cleancut        = (ev.mass_b1b2>30)
  dR_lljj_cleancut    = (ev.dR_l1l2b1b2>0.2 and ev.dR_l1l2b1b2<4)
  dR_ll_cleancut      = (ev.dR_l1l2<3.8)
  dR_jj_cleancut      = (ev.dR_b1b2<3.8)
  cleaning_cuts       = (mt_cleancut and NoZ_cleancut and mJJ_cleancut and dR_lljj_cleancut and dR_ll_cleancut and dR_jj_cleancut)
  # WEIGHTS
  weight = 1.
  if( whichSample!="Data" ):
    # Cross section (We only weights to 1fb-1, the real normalization for xsec will be placed into FinalPlotter.py)
    xSec = ev.XsecBr
    if(nEv==0): h_pre_XsecBr.SetBinContent(1, xSec)
    Lumi    = 1. * 1000 # Convert from 1fb-1 to 1pb-1
    weight = weight * (Lumi*xSec)/nTOT_prehlt
    # Scale factors
    weight = weight * ev.muon1_pogSF * ev.muon2_pogSF
  # Minimal Selection
  if( preselection ):
    # SF
    h_pre_muon1_triggerSF.Fill(abs(ev.muon1_eta), ev.muon1_pt, ev.muon1_triggerSF)
    h_pre_muon1_isoSF.Fill(abs(ev.muon1_eta), ev.muon1_pt, ev.muon1_isoSF)
    h_pre_muon1_idSF.Fill(abs(ev.muon1_eta), ev.muon1_pt, ev.muon1_idSF)
    h_pre_muon1_trackingSF.Fill(ev.muon1_eta, ev.muon1_trackingSF)
    h_pre_muon1_SF_bg1.Fill(ev.muon1_eta)
    h_pre_muon1_SF_bg2.Fill(abs(ev.muon1_eta), ev.muon1_pt)
    h_pre_muon2_triggerSF.Fill(abs(ev.muon2_eta), ev.muon2_pt, ev.muon2_triggerSF)
    h_pre_muon2_isoSF.Fill(abs(ev.muon2_eta), ev.muon2_pt, ev.muon2_isoSF)
    h_pre_muon2_idSF.Fill(abs(ev.muon2_eta), ev.muon2_pt, ev.muon2_idSF)
    h_pre_muon2_trackingSF.Fill(ev.muon2_eta, ev.muon2_trackingSF)
    h_pre_muon2_SF_bg1.Fill(ev.muon2_eta)
    h_pre_muon2_SF_bg2.Fill(abs(ev.muon2_eta), ev.muon2_pt)
    # Regression Variables
    h_pre_numOfVertices.Fill( ev.numOfVertices, weight )   
    h_pre_b1jet_mt.Fill( ev.b1jet_mt, weight )
    h_pre_b1jet_leadTrackPt.Fill( ev.b1jet_leadTrackPt, weight )
    h_pre_b1jet_leptonPtRel.Fill( ev.b1jet_leptonPtRel, weight )
    h_pre_b1jet_leptonPt.Fill( ev.b1jet_leptonPt, weight )
    h_pre_b1jet_leptonDeltaR.Fill( ev.b1jet_leptonDeltaR, weight )
    h_pre_b1jet_neHEF.Fill( ev.b1jet_neHEF, weight )
    h_pre_b1jet_neEmEF.Fill( ev.b1jet_neEmEF, weight )
    h_pre_b1jet_vtxNtracks.Fill( ev.b1jet_vtxNtracks, weight )
    h_pre_b1jet_vtxPt.Fill( ev.b1jet_vtxPt, weight )
    h_pre_b1jet_vtxMass.Fill( ev.b1jet_vtxMass, weight )
    h_pre_b1jet_vtx3DSig.Fill( ev.b1jet_vtx3DSig, weight )
    h_pre_b1jet_vtx3DVal.Fill( ev.b1jet_vtx3DVal, weight )
    # Kinematic Variables
    if ev.findAllGenParticles: 
	# Gen
        if ev.w1_mass > ev.w2_mass:
		h_gen_pre_onshellWmass.Fill( ev.w1_mass )
		h_gen_pre_offshellWmass.Fill( ev.w2_mass )
		h_gen_pre_onoffshellWmass.Fill( ev.w1_mass, ev.w2_mass )
		h_gen_pre_onshellnupt.Fill( ev.nu1_pt )
        else:
		h_gen_pre_onshellWmass.Fill( ev.w2_mass )
		h_gen_pre_offshellWmass.Fill( ev.w1_mass )
		h_gen_pre_onoffshellWmass.Fill( ev.w2_mass, ev.w1_mass )
		h_gen_pre_onshellnupt.Fill( ev.nu2_pt )
	h_gen_pre_MU1_pt.Fill( ev.mu1_pt, weight )
	h_gen_pre_MU2_pt.Fill( ev.mu2_pt, weight )
	h_gen_pre_MU2_eta.Fill( ev.mu2_eta, weight )
	h_gen_pre_MU1_eta.Fill( ev.mu1_eta, weight )
	h_gen_pre_mass_l1l2.Fill( ev.mass_genl1l2, weight )
	h_gen_pre_dR_l1l2.Fill( ev.dR_genl1l2, weight )
	h_gen_pre_J1_pt.Fill( ev.b1_pt, weight )
	h_gen_pre_J2_pt.Fill( ev.b2_pt, weight )
	h_gen_pre_J1_eta.Fill( ev.b1_eta, weight )
	h_gen_pre_J2_eta.Fill( ev.b2_eta, weight )
	h_gen_pre_mass_b1b2.Fill( ev.mass_genb1b2, weight )
	h_gen_pre_dR_b1b2.Fill( ev.dR_genb1b2, weight )
        dR_jet11 = deltaR(ev.b1jet_eta, ev.b1jet_phi, ev.b1_eta, ev.b1_phi)
        dR_jet12 = deltaR(ev.b1jet_eta, ev.b1jet_phi, ev.b2_eta, ev.b2_phi)
        dR_jet21 = deltaR(ev.b2jet_eta, ev.b2jet_phi, ev.b1_eta, ev.b1_phi)
        dR_jet22 = deltaR(ev.b2jet_eta, ev.b2jet_phi, ev.b2_eta, ev.b2_phi)
        if dR_jet11 < dR_jet12 and dR_jet11 < 0.4:
	   h_gen_pre_bjetrescalec1dR4.Fill( ev.b1_pt/ev.b1jet_pt, weight ) 
        elif dR_jet11 > dR_jet12 and dR_jet12 < 0.4:
	   h_gen_pre_bjetrescalec1dR4.Fill( ev.b2_pt/ev.b1jet_pt, weight ) 
        if dR_jet22 < dR_jet21 and dR_jet22 < 0.4:
	   h_gen_pre_bjetrescalec2dR4.Fill( ev.b2_pt/ev.b2jet_pt, weight ) 
        elif dR_jet22 > dR_jet21 and dR_jet21 < 0.4:
	  h_gen_pre_bjetrescalec2dR4.Fill( ev.b1_pt/ev.b2jet_pt, weight ) 
        #print "dR_jet11 ",dR_jet11," dR_jet12 ",dR_jet12, " ev.dR_b1jet ",ev.dR_b1jet
        #print "dR_jet21 ",dR_jet21," dR_jet22 ",dR_jet22, " ev.dR_b2jet ",ev.dR_b2jet
	h_gen_pre_met_pt.Fill( ev.genmet_pt, weight )
	h_gen_pre_mass_trans.Fill( ev.mass_gentrans, weight )
	h_gen_pre_dR_l1l2b1b2.Fill( ev.dR_genl1l2b1b2, weight )
	h_gen_pre_dphi_llmet.Fill( ev.dphi_genllmet, weight )
    #Reco
    h_pre_MU1_pt.Fill( ev.muon1_pt, weight )
    h_pre_MU2_pt.Fill( ev.muon2_pt, weight )
    h_pre_MU2_eta.Fill( ev.muon2_eta, weight )
    h_pre_MU1_eta.Fill( ev.muon1_eta, weight )
    h_pre_mass_l1l2.Fill( ev.mass_l1l2, weight )
    h_pre_dR_l1l2.Fill( ev.dR_l1l2, weight )
    h_pre_J1_pt.Fill( ev.b1jet_pt, weight )
    h_pre_J2_pt.Fill( ev.b2jet_pt, weight )
    h_pre_J1_eta.Fill( ev.b1jet_eta, weight )
    h_pre_J2_eta.Fill( ev.b2jet_eta, weight )
    h_pre_mass_b1b2.Fill( ev.mass_b1b2, weight )
    h_pre_dR_b1b2.Fill( ev.dR_b1b2, weight )
    h_pre_met_pt.Fill( ev.met_pt, weight )
    h_pre_mass_trans.Fill( ev.mass_trans, weight )
    h_pre_dR_l1l2b1b2.Fill( ev.dR_l1l2b1b2, weight )
    h_pre_dphi_llmet.Fill( ev.dphi_llmet, weight )

    if( cleaning_cuts ):
      # Kinematic Variables
      if ev.findAllGenParticles: 
	  #Gen
          if ev.w1_mass > ev.w2_mass:
		h_gen_cc_onshellWmass.Fill( ev.w1_mass)
		h_gen_cc_offshellWmass.Fill( ev.w2_mass)
		h_gen_cc_onoffshellWmass.Fill( ev.w1_mass, ev.w2_mass )
		h_gen_cc_onshellnupt.Fill( ev.nu1_pt, weight)
          else:
		h_gen_cc_onshellWmass.Fill( ev.w2_mass)
		h_gen_cc_offshellWmass.Fill( ev.w1_mass)
		h_gen_cc_onoffshellWmass.Fill( ev.w2_mass, ev.w1_mass )
		h_gen_cc_onshellnupt.Fill( ev.nu2_pt )
	  h_gen_cc_MU1_pt.Fill( ev.mu1_pt, weight )
	  h_gen_cc_MU2_pt.Fill( ev.mu2_pt, weight )
	  h_gen_cc_MU2_eta.Fill( ev.mu2_eta, weight )
	  h_gen_cc_MU1_eta.Fill( ev.mu1_eta, weight )
	  h_gen_cc_mass_l1l2.Fill( ev.mass_genl1l2, weight )
	  h_gen_cc_dR_l1l2.Fill( ev.dR_genl1l2, weight )
	  h_gen_cc_J1_pt.Fill( ev.b1_pt, weight )
	  h_gen_cc_J2_pt.Fill( ev.b2_pt, weight )
	  h_gen_cc_J1_eta.Fill( ev.b1_eta, weight )
	  h_gen_cc_J2_eta.Fill( ev.b2_eta, weight )
	  h_gen_cc_mass_b1b2.Fill( ev.mass_genb1b2, weight )
	  h_gen_cc_dR_b1b2.Fill( ev.dR_genb1b2, weight )
	  h_gen_cc_met_pt.Fill( ev.genmet_pt, weight )
	  h_gen_cc_mass_trans.Fill( ev.mass_gentrans, weight )
	  h_gen_cc_dR_l1l2b1b2.Fill( ev.dR_genl1l2b1b2, weight )
	  h_gen_cc_dphi_llmet.Fill( ev.dphi_genllmet, weight )
          dR_jet11 = deltaR(ev.b1jet_eta, ev.b1jet_phi, ev.b1_eta, ev.b1_phi)
          dR_jet12 = deltaR(ev.b1jet_eta, ev.b1jet_phi, ev.b2_eta, ev.b2_phi)
          dR_jet21 = deltaR(ev.b2jet_eta, ev.b2jet_phi, ev.b1_eta, ev.b1_phi)
          dR_jet22 = deltaR(ev.b2jet_eta, ev.b2jet_phi, ev.b2_eta, ev.b2_phi)
          if dR_jet11 < dR_jet12 and dR_jet11 < 0.4:
	   h_gen_cc_bjetrescalec1dR4.Fill( ev.b1_pt/ev.b1jet_pt, weight ) 
          elif dR_jet11 > dR_jet12 and dR_jet12 < 0.4:
	   h_gen_cc_bjetrescalec1dR4.Fill( ev.b2_pt/ev.b1jet_pt, weight ) 
          if dR_jet22 < dR_jet21 and dR_jet22 < 0.4:
	   h_gen_cc_bjetrescalec2dR4.Fill( ev.b2_pt/ev.b2jet_pt, weight ) 
          elif dR_jet22 > dR_jet21 and dR_jet21 < 0.4:
	   h_gen_cc_bjetrescalec2dR4.Fill( ev.b1_pt/ev.b2jet_pt, weight ) 
      #REco
      h_cc_MU1_pt.Fill( ev.muon1_pt, weight )
      h_cc_MU2_pt.Fill( ev.muon2_pt, weight )
      h_cc_MU2_eta.Fill( ev.muon2_eta, weight )
      h_cc_MU1_eta.Fill( ev.muon1_eta, weight )
      h_cc_mass_l1l2.Fill( ev.mass_l1l2, weight )
      h_cc_dR_l1l2.Fill( ev.dR_l1l2, weight )
      h_cc_J1_pt.Fill( ev.b1jet_pt, weight )
      h_cc_J2_pt.Fill( ev.b2jet_pt, weight )
      h_cc_J1_eta.Fill( ev.b1jet_eta, weight )
      h_cc_J2_eta.Fill( ev.b2jet_eta, weight )
      h_cc_mass_b1b2.Fill( ev.mass_b1b2, weight )
      h_cc_dR_b1b2.Fill( ev.dR_b1b2, weight )
      h_cc_met_pt.Fill( ev.met_pt, weight )
      h_cc_mass_trans.Fill( ev.mass_trans, weight )
      h_cc_dR_l1l2b1b2.Fill( ev.dR_l1l2b1b2, weight )
      h_cc_dphi_llmet.Fill( ev.dphi_llmet, weight )
  nEv = nEv +1

h_pre_muon1_triggerSF.Divide(h_pre_muon1_SF_bg2)
h_pre_muon1_isoSF.Divide(h_pre_muon1_SF_bg2)
h_pre_muon1_idSF.Divide(h_pre_muon1_SF_bg2)
h_pre_muon1_trackingSF.Divide(h_pre_muon1_SF_bg1)
h_pre_muon2_triggerSF.Divide(h_pre_muon2_SF_bg2)
h_pre_muon2_isoSF.Divide(h_pre_muon2_SF_bg2)
h_pre_muon2_idSF.Divide(h_pre_muon2_SF_bg2)
h_pre_muon2_trackingSF.Divide(h_pre_muon2_SF_bg1)
#normalize1D(h_pre_muon1_trackingSF)
# Writing histograms
f.cd()
 
h_pre_XsecBr.Write()
h_pre_muon1_triggerSF.Write()
h_pre_muon1_isoSF.Write()
h_pre_muon1_idSF.Write()
h_pre_muon1_trackingSF.Write()
h_pre_muon2_triggerSF.Write()
h_pre_muon2_isoSF.Write()
h_pre_muon2_idSF.Write()
h_pre_muon2_trackingSF.Write()
h_pre_Nev_preHLT.Write()
h_pre_Nev_posHLT.Write()
h_pre_numOfVertices.Write()  
h_pre_b1jet_mt.Write()
h_pre_b1jet_leadTrackPt.Write()
h_pre_b1jet_leptonPtRel.Write()
h_pre_b1jet_leptonPt.Write()
h_pre_b1jet_leptonDeltaR.Write()
h_pre_b1jet_neHEF.Write()
h_pre_b1jet_neEmEF.Write()
h_pre_b1jet_vtxNtracks.Write()
h_pre_b1jet_vtxPt.Write()
h_pre_b1jet_vtxMass.Write()
h_pre_b1jet_vtx3DSig.Write()
h_pre_b1jet_vtx3DVal.Write()

if whichSample != "DATA":
    h_gen_pre_MU1_pt.Write()
    h_gen_pre_MU1_pt.Write()
    h_gen_pre_MU2_pt.Write()
    h_gen_pre_MU1_eta.Write()
    h_gen_pre_MU2_eta.Write()
    h_gen_pre_mass_l1l2.Write()
    h_gen_pre_dR_l1l2.Write()
    h_gen_pre_J1_pt.Write()
    h_gen_pre_J2_pt.Write()
    h_gen_pre_J1_eta.Write()
    h_gen_pre_J2_eta.Write()
    h_gen_pre_mass_b1b2.Write()
    h_gen_pre_dR_b1b2.Write()
    h_gen_pre_met_pt.Write()
    h_gen_pre_mass_trans.Write()
    h_gen_pre_dR_l1l2b1b2.Write()
    h_gen_pre_onshellWmass.Write()
    h_gen_pre_offshellWmass.Write()
    h_gen_pre_onshellnupt.Write()
    h_gen_pre_onoffshellWmass.Write()
    h_gen_pre_bjetrescalec1dR4.Write()
    h_gen_pre_bjetrescalec2dR4.Write()
h_pre_MU1_pt.Write()
h_pre_MU1_pt.Write()
h_pre_MU2_pt.Write()
h_pre_MU1_eta.Write()
h_pre_MU2_eta.Write()
h_pre_mass_l1l2.Write()
h_pre_dR_l1l2.Write()
h_pre_J1_pt.Write()
h_pre_J2_pt.Write()
h_pre_J1_eta.Write()
h_pre_J2_eta.Write()
h_pre_mass_b1b2.Write()
h_pre_dR_b1b2.Write()
h_pre_met_pt.Write()
h_pre_mass_trans.Write()
h_pre_dR_l1l2b1b2.Write()
h_pre_dphi_llmet.Write()
# Cleaning Cuts
if whichSample != "DATA":
    h_gen_cc_MU1_pt.Write()
    h_gen_cc_MU2_pt.Write()
    h_gen_cc_MU2_eta.Write()
    h_gen_cc_MU1_eta.Write()
    h_gen_cc_mass_l1l2.Write()
    h_gen_cc_dR_l1l2.Write()
    h_gen_cc_J1_pt.Write()
    h_gen_cc_J2_pt.Write()
    h_gen_cc_J1_eta.Write()
    h_gen_cc_J2_eta.Write()
    h_gen_cc_mass_b1b2.Write()
    h_gen_cc_dR_b1b2.Write()
    h_gen_cc_met_pt.Write()
    h_gen_cc_mass_trans.Write()
    h_gen_cc_dR_l1l2b1b2.Write()
    h_gen_cc_dphi_llmet.Write()
    h_gen_cc_onshellWmass.Write()
    h_gen_cc_offshellWmass.Write()
    h_gen_cc_onshellnupt.Write()
    h_gen_cc_onoffshellWmass.Write()
    h_gen_cc_bjetrescalec1dR4.Write()
    h_gen_cc_bjetrescalec2dR4.Write()
h_cc_MU1_pt.Write()
h_cc_MU2_pt.Write()
h_cc_MU2_eta.Write()
h_cc_MU1_eta.Write()
h_cc_mass_l1l2.Write()
h_cc_dR_l1l2.Write()
h_cc_J1_pt.Write()
h_cc_J2_pt.Write()
h_cc_J1_eta.Write()
h_cc_J2_eta.Write()
h_cc_mass_b1b2.Write()
h_cc_dR_b1b2.Write()
h_cc_met_pt.Write()
h_cc_mass_trans.Write()
h_cc_dR_l1l2b1b2.Write()
h_cc_dphi_llmet.Write()
f.Close()
print "Done."
