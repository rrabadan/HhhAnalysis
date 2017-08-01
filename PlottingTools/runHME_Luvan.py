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
def AddFiles(chain, inputDir):
    theInputFiles = []
    for d in range(len(inputDir)):
        my_dir = inputDir[d]
        if not os.path.isdir(my_dir):
            print "ERROR: This is not a valid directory: ", my_dir
            if d==len(inputDir)-1:
                print "ERROR: No input files were selected"
                exit()
            continue
        print "Proceed to next directory"
        ls = os.listdir(my_dir)
        chain.Add(my_dir[:] + x for x in ls if x.endswith('root'))

doTest = False
doHME = True
refPDF = ROOT.TFile("REFPDFPU40.root","READ")
onshellWmasspdf = refPDF.Get("onshellWmasspdf")
onshellnuptpdf = refPDF.Get("onshellnuptpdf")
recobjetrescalec1pdfPU40 = refPDF.Get("recobjetrescalec1pdfPU40")

tree_name="DiHiggsWWBBAna/evtree"
tree_name = "t"
benchmarks = ["ttV", "Wjet", "singTop", "VV", "DY", "TTbar", "Data"]
#whichSample = "Graviton_0419"

parser = argparse.ArgumentParser(description='runHME')
parser.add_argument("-n", "--njobs", dest="njobs", type=int, default=100, help="total splitted jobs. [Default: 100]")
parser.add_argument("-i", "--ijob", dest="ijob", type=int, default=0, help="ith job in splitting. [Default: 0]")
parser.add_argument("-jt", "--jobtype", dest="jobtype", default="radion_M300", help="sample type to run. [Default: radion]")
parser.add_argument("-it", "--iterations", dest="iterations", default=10000, help="iteration number used HME [Default: 100000]")
args = parser.parse_args()
whichSample = args.jobtype
makeHadd = "HaddNo"

#print "Executing: python", sys.argv[0] , "-b", sys.argv[2], sys.argv[3], "(Arg1=makeHadd: HaddYes or HaddNo /|\ Arg2=whichSample: TT, DY, VV, sT, Wjet, ttV, Data)"
#makeHadd = sys.argv[2]
#whichSample = sys.argv[3]
if( makeHadd!="HaddYes" and makeHadd!="HaddNo" ): print "WARNING! 1st argument has to be HaddYes or HaddNo"; sys.exit()
#if( whichSample!="TT" and whichSample!="DY" and whichSample!="VV" and whichSample!="sT" and whichSample!="Wjet" and whichSample!="ttV" and whichSample!="Data" ):  print "WARNING! 2nd argument have to be TT, DY, VV, singTop, Wjet, ttV, or Data"; sys.exit()
######################################
## define shell commands
######################################
print "sampletype ",whichSample, " ",Samplelist.outAnalist[whichSample]
Find_str = []; this_cat = ""; this_hadd = ""; this_NtotPath = ""
for i in range(len(Samplelist.outAnalist[whichSample])):
  Find_str.append("find %s | grep root | grep -v failed > HADD/%s_%d.txt"%(Samplelist.outAnalist[whichSample][i], Samplelist.sampleFullName[whichSample], i))
    #Find_str.append("find /fdata/hepx/store/user/lpernie/TTTo2L2Nu_13TeV-powheg/crab_TTTo2L2Nu_13TeV-powheg/170405_172215 | grep root | grep -v failed > HADD/TT_TTTo2L2Nu13TeV-powheg.txt")
this_cat = "cat HADD/%s_* > HADD/%s.txt"%(Samplelist.sampleFullName[whichSample], whichSample)
#if whichSample != "DATA":
    #this_NtotPath = "/fdata/hepx/store/user/%s/%s/crab_%s.root"%(user,Samplelist.sampleFullName[whichSample], Samplelist.sampleFullName[whichSample])
    #this_NtotPath = Samplelist.sampleFullName[whichSample]
    #this_hadd     = "hadd -T -f -k %s @HADD/%s.txt"%(this_NtotPath, whichSample)
    #outpath = "/fdata/hepx/store/user/%s/%s"%(user,Samplelist.sampleFullName[whichSample])
    #os.system("mkdir -p %s" % outpath)

######################################
# Running commands
######################################
for this_find in Find_str:
  print "this_find ",this_find
  os.system(this_find)

os.system(this_cat)
print "this_hadd ", this_hadd
#if (makeHadd=="HaddYes" and whichSample!="Data"): os.system(this_hadd)
######################################
#this_NtotPath = "/home/taohuang/DiHiggsAnalysis/CMSSW_8_0_26_patch2/src/HhhAnalysis/CutFlowAnalyzer/test/out_ana_10k_graviton.root"
#this_NtotPath = "/fdata/hepx/store/user/taohuang/DiHiggsAnalysisSample/out_ann_%s_20160411.root"%whichSample
if( whichSample != "Data" ):
  #Ntot_path = this_NtotPath
  #MyFile =  ROOT.TFile.Open(Ntot_path,"read");
  #h_prehlt  =  ROOT.TH1F(MyFile.Get("TriggerResults/hevent_filter")); nTOT_prehlt = h_prehlt.GetBinContent(2)
  #h_posthlt =  ROOT.TH1F(MyFile.Get("DiHiggsWWBBAna/hevent")); nTOT_posthlt = h_posthlt.GetBinContent(2);
  nTOT_prehlt = 1.; nTOT_posthlt = 1
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
f = ROOT.TFile('/fdata/hepx/store/user/%s/Hhh_For_Test/'%user + whichSample + '_ijob%d.root'%args.ijob,'recreate'); f.cd()
#TCha2 = TCha.Clone()
sample_weight       = np.zeros(1, dtype=float) 
jj_M                = np.zeros(1, dtype=float) 
lep1_pt             = np.zeros(1, dtype=float) 
lep2_pt             = np.zeros(1, dtype=float) 
jet1_pt             = np.zeros(1, dtype=float) 
jet2_pt             = np.zeros(1, dtype=float) 
met_pt              = np.zeros(1, dtype=float) 
ht                  = np.zeros(1, dtype=float) 
llmetjj_M           = np.zeros(1, dtype=float) 
cosThetaStar        = np.zeros(1, dtype=float) 
ll_M                = np.zeros(1, dtype=float) 
ll_DR_l_l           = np.zeros(1, dtype=float) 
jj_DR_j_j           = np.zeros(1, dtype=float) 
llmetjj_DPhi_ll_jj  = np.zeros(1, dtype=float) 
llmetjj_MTformula   = np.zeros(1, dtype=float) 
ll_pt               = np.zeros(1, dtype=float) 
jj_pt               = np.zeros(1, dtype=float) 
llmetjj_minDR_l_j   = np.zeros(1, dtype=float) 
lep1_eta            = np.zeros(1, dtype=float) 
lep1_phi            = np.zeros(1, dtype=float) 
lep1_Iso            = np.zeros(1, dtype=float) 
lep2_eta            = np.zeros(1, dtype=float) 
lep2_phi            = np.zeros(1, dtype=float) 
lep2_Iso            = np.zeros(1, dtype=float) 
ll_DPhi_l_l         = np.zeros(1, dtype=float) 
ll_DEta_l_l         = np.zeros(1, dtype=float) 
event_weight        = np.zeros(1, dtype=float) 
event_pu_weight     = np.zeros(1, dtype=float) 
isElEl              = np.zeros(1, dtype=float) 
isMuMu              = np.zeros(1, dtype=float) 
isElMu              = np.zeros(1, dtype=float) 
isMuEl              = np.zeros(1, dtype=float) 
event_number        = np.zeros(1, dtype=float) 
event_run           = np.zeros(1, dtype=float)
isSF                = np.zeros(1, dtype=float) 
total_weight        = np.zeros(1, dtype=float) 
jet1_cMVAv2         = np.zeros(1, dtype=float) 
jet2_cMVAv2         = np.zeros(1, dtype=float) 
nJetsL              = np.zeros(1, dtype=float) 
llidiso             = np.zeros(1, dtype=float)
mumuidiso           = np.zeros(1, dtype=float) 
elelidiso           = np.zeros(1, dtype=float) 
jjbtag_heavy        = np.zeros(1, dtype=float) 
jjbtag_light        = np.zeros(1, dtype=float) 
trigeff             = np.zeros(1, dtype=float)
pu                  = np.zeros(1, dtype=float) 



TCha2 = ROOT.Tree("evtreeHME","tree with HME")
TCha2.Branch("sample_weight",         sample_weight,       "sample_weight/F")
TCha2.Branch("jj_M",                  jj_M,                "jj_M/F")
TCha2.Branch("lep1_pt",               lep1_pt,             "lep1_pt/F")
TCha2.Branch("lep2_pt",               lep2_pt,             "lep2_pt/F")
TCha2.Branch("jet1_pt",               jet1_pt,             "jet1_pt/F")
TCha2.Branch("jet2_pt",               jet2_pt,             "jet2_pt/F")
TCha2.Branch("met_pt",                met_pt,              "met_pt/F")
TCha2.Branch("ht",                    ht,                  "ht/F")
TCha2.Branch("llmetjj_M",             llmetjj_M,           "llmetjj_M/F")
TCha2.Branch("cosThetaStar",          cosThetaStar,        "cosThetaStar/F")
TCha2.Branch("ll_M",                  ll_M,                "ll_M/F")
TCha2.Branch("ll_DR_l_l",             ll_DR_l_l,           "ll_DR_l_l/F")
TCha2.Branch("jj_DR_j_j",             jj_DR_j_j,           "jj_DR_j_j/F")
TCha2.Branch("llmetjj_DPhi_ll_jj",    llmetjj_DPhi_ll_jj,  "llmetjj_DPhi_ll_jj/F")
TCha2.Branch("llmetjj_MTformula",     llmetjj_MTformula,   "llmetjj_MTformula/F")
TCha2.Branch("ll_pt",                 ll_pt,                "ll_pt/F")
TCha2.Branch("jj_pt",                 jj_pt,                "ll_pt/F")
TCha2.Branch("llmetjj_minDR_l_j",     llmetjj_minDR_l_j,    "llmetjj_minDR_l_j/F")
TCha2.Branch("lep1_eta",              lep1_eta,             "lep1_eta/F")
TCha2.Branch("lep1_phi",              lep1_phi,             "lep1_phi/F")
TCha2.Branch("lep1_Iso",              lep1_Iso,             "lep1_Iso/F")
TCha2.Branch("lep2_eta",              lep2_eta,             "lep2_eta/F")
TCha2.Branch("lep2_phi",              lep2_phi,             "lep2_phi/F")
TCha2.Branch("lep2_Iso",              lep2_Iso,             "lep2_Iso/F")
TCha2.Branch("ll_DPhi_l_l",           ll_DPhi_l_l,          "ll_DPhi_l_l/F")
TCha2.Branch("ll_DEta_l_l",           ll_DEta_l_l,          "ll_DEta_l_l/F")
TCha2.Branch("event_weight",          event_weight,         "event_weight/F")
TCha2.Branch("event_pu_weight",       event_pu_weight,      "event_pu_weight/F")
TCha2.Branch("isElEl",                isElEl,               "isElEl/F")
TCha2.Branch("isMuMu",                isMuMu,               "isMuMu/F")
TCha2.Branch("isElMu",                isElMu,               "isElMu/F")
TCha2.Branch("isMuEl",                isMuEl,               "isMuEl/F")
TCha2.Branch("event_number",          event_number,         "event_number/F")
TCha2.Branch("event_run",             event_run,            "event_run/F")
TCha2.Branch("isSF",                  isSF,                 "isSF/F")
TCha2.Branch("total_weight",          total_weight,         "total_weight/F")
TCha2.Branch("jet1_cMVAv2",           jet1_cMVAv2,          "jet1_cMVAv2/F")
TCha2.Branch("jet2_cMVAv2",           jet2_cMVAv2,          "jet2_cMVAv2/F")
TCha2.Branch("nJetsL",                nJetsL,               "nJetsL/F")
TCha2.Branch("llidiso",               llidiso,              "llidiso/F")
TCha2.Branch("mumuidiso",             mumuidiso,            "mumuidiso/F")
TCha2.Branch("elelidiso",             elelidiso,            "elelidiso/F")
TCha2.Branch("jjbtag_heavy",          jjbtag_heavy,         "jjbtag_heavy/F")
TCha2.Branch("jjbtag_light",          jjbtag_light,         "jjbtag_light/F")
TCha2.Branch("trigeff",               trigeff,              "trigeff/F")
TCha2.Branch("pu",                    pu,                   "pu/F")



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

h_h2mass_weight_gen      = ROOT.TH1F("h_h2mass_weight_gen","",100,200.,1200.);  h_h2mass_weight_gen.GetXaxis().SetTitle("HME reco mass [GeV]");
h_h2mass_weight_gen_sum  = ROOT.TH1F("h_h2mass_weight_gen_sum","",1000,200.,1200.);  h_h2mass_weight_gen_sum.GetXaxis().SetTitle("HME reco mass [GeV]");
h_h2mass_weight1_gen     = ROOT.TH1F("h_h2mass_weight1_gen","",100,200.,1200.);  h_h2mass_weight1_gen.GetXaxis().SetTitle("HME reco mass [GeV]");
h_h2mass_weight_reco     = ROOT.TH1F("h_h2mass_weight_reco","",100,200.,1200.);  h_h2mass_weight_reco.GetXaxis().SetTitle("HME reco mass [GeV]");
h_h2mass_weight_reco_sum = ROOT.TH1F("h_h2mass_weight_reco_sum","",1000,200.,1200.);  h_h2mass_weight_reco_sum.GetXaxis().SetTitle("HME reco mass [GeV]");
h_h2mass_weight1_reco    = ROOT.TH1F("h_h2mass_weight1_reco","",100,200.,1200.);  h_h2mass_weight1_reco.GetXaxis().SetTitle("HME reco mass [GeV]");
offshellWmass_gen	 = ROOT.TH1F("hme_offshellWmass_gen","offshell W mass from HME", 60, 10.0, 70.0); offshellWmass_gen.GetXaxis().SetTitle("HME offshell W mass [GeV]");
offshellWmass_reco	 = ROOT.TH1F("hme_offshellWmass_reco","offshell W mass from HME", 60, 10.0, 70.0); offshellWmass_reco.GetXaxis().SetTitle("HME offshell W mass [GeV]");
# adding more branches to TTree
passCC         = np.zeros(1, dtype=bool)
passpre        = np.zeros(1, dtype=bool)
h_h2mass_gen   = np.zeros(1, dtype=float)
h_h2mass_reco  = np.zeros(1, dtype=float)
passCCbr       = TCha2.Branch("passCC", passCC,"passCC/I")
passprebr      = TCha2.Branch("passpre", passpre,"passpre/I")
h_h2mass_genbr = TCha2.Branch("h_h2mass_gen", h_h2mass_gen,"h_h2mass_gen/I")
h_h2mass_recobr= TCha2.Branch("h_h2mass_reco", h_h2mass_reco,"h_h2mass_reco/I")
def initbr():
    passCC[0] = False
    passpre[0] = False
    h_h2mass_gen[0] = -1
    h_h2mass_reco[0] = -1


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

  if ev.htoWW>0 and ev.w1_mass>0 and ev.w2_mass>0 :
      print "WWBB event ",nEv," w1_mass ",ev.w1_mass, " w2_mass ",ev.w2_mass
  elif ev.htoZZ>0 and ev.z1_mass>0 and ev.z2_mass>0:
      print "ZZBB event ",nEv, " z1_mass ",ev.z1_mass, " z2_mass ",ev.z2_mass
  else:
      print "NOTZZBBWWZ event"
  initbr()

  sample_weight[0]      =  ev.sample_weight      
  jj_M[0]               =  ev.jj_M               
  lep1_pt[0]            =  ev.lep1_pt            
  lep2_pt[0]            =  ev.lep2_pt            
  jet1_pt[0]            =  ev.jet1_pt            
  jet2_pt[0]            =  ev.jet2_pt            
  met_pt[0]             =  ev.met_pt             
  ht[0]                 =  ev.ht                 
  llmetjj_M[0]          =  ev.llmetjj_M          
  cosThetaStar[0]       =  ev.cosThetaStar       
  ll_M[0]               =  ev.ll_M               
  ll_DR_l_l[0]          =  ev.ll_DR_l_l          
  jj_DR_j_j[0]          =  ev.jj_DR_j_j          
  llmetjj_DPhi_ll_jj[0] =  ev.llmetjj_DPhi_ll_jj 
  llmetjj_MTformula[0]  =  ev.llmetjj_MTformula  
  ll_pt[0]              =  ev.ll_pt              
  jj_pt[0]              =  ev.jj_pt              
  llmetjj_minDR_l_j[0]  =  ev.llmetjj_minDR_l_j  
  lep1_eta[0]           =  ev.lep1_eta           
  lep1_phi[0]           =  ev.lep1_phi           
  lep1_Iso[0]           =  ev.lep1_Iso           
  lep2_eta[0]           =  ev.lep2_eta           
  lep2_phi[0]           =  ev.lep2_phi           
  lep2_Iso[0]           =  ev.lep2_Iso           
  ll_DPhi_l_l[0]        =  ev.ll_DPhi_l_l        
  ll_DEta_l_l[0]        =  ev.ll_DEta_l_l        
  event_weight[0]       =  ev.event_weight       
  event_pu_weight[0]    =  ev.event_pu_weight    
  isElEl[0]             =  ev.isElEl             
  isMuMu[0]             =  ev.isMuMu             
  isElMu[0]             =  ev.isElMu             
  isMuEl[0]             =  ev.isMuEl             
  event_number[0]       =  ev.event_number       
  event_run[0]          =  ev.event_run          
  isSF[0]               =  ev.isSF               
  total_weight[0]       =  ev.total_weight       
  jet1_cMVAv2[0]        =  ev.jet1_cMVAv2        
  jet2_cMVAv2[0]        =  ev.jet2_cMVAv2        
  nJetsL[0]             =  ev.nJetsL             
  llidiso[0]            =  ev.llidiso            
  mumuidiso[0]          =  ev.mumuidiso          
  elelidiso[0]          =  ev.elelidiso          
  jjbtag_heavy[0]       =  ev.jjbtag_heavy       
  jjbtag_light[0]       =  ev.jjbtag_light       
  trigeff[0]            =  ev.trigeff              
  pu[0]                 =  ev.pu                 
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
  if (preselection and cleaning_cuts ):
      print "event ",nEv," pass CC cut"
  if (preselection):
      passpre[0] = True
  if (cleaning_cuts):
      passCC[0] = True
  
  if (cleaning_cuts and ev.findAllGenParticles and doHME):
      
      mu1_p4 	      = ROOT.TLorentzVector(ev.mu1_px, ev.mu1_py, ev.mu1_pz, ev.mu1_energy)
      mu2_p4 	      = ROOT.TLorentzVector(ev.mu2_px, ev.mu2_py, ev.mu2_pz, ev.mu2_energy)
      b1jet_p4 	      = ROOT.TLorentzVector(ev.b1_px, ev.b1_py, ev.b1_pz, ev.b1_energy)
      b2jet_p4 	      = ROOT.TLorentzVector(ev.b2_px, ev.b2_py, ev.b2_pz, ev.b2_energy)
      met_vec2        = ROOT.TVector2(ev.nu1_px + ev.nu2_px, ev.nu1_py + ev.nu2_py)
      hme_gen = HeavyMassEstimator()
      onshellW_mu  = 1 
      if ev.w1_mass < ev.w2_mass:
      	  onshellW_mu  = 2
      hme_gen.setKinematic(mu1_p4, mu2_p4, b1jet_p4, b2jet_p4, met_vec2, onshellW_mu)
      hme_gen.setonshellWmasspdf(onshellWmasspdf)
      hme_gen.setonshellnuptpdf(onshellnuptpdf)
      #hme.showKinematic()
      hme_gen.setIterations(args.iterations)
      hme_gen.runHME()
      if hme_gen.hme_h2Mass.GetEntries()>0:
	  print "Gen Level most probable reco mass ",hme_gen.hme_h2Mass.GetXaxis().GetBinCenter(hme_gen.hme_h2Mass.GetMaximumBin())," entries ",hme_gen.hme_h2Mass.GetEntries()
	  hme_gen.hme_h2Mass.SetName("hme_h2Mass_ev%d_genlevel"%nEv)
          hme_gen.hme_h2MassWeight1.SetName("hme_h2MassWeight1_ev%d_genlevel"%nEv)
	  hme_gen.hmetree.SetName("hmetree_ev%d_genlevel"%nEv)
	  hme_gen.hme_h2Mass_correctmunupair.SetName("hme_h2Mass_correctmunupair_ev%d_genlevel"%nEv)
	  hme_gen.hme_h2Mass_incorrectmunupair.SetName("hme_h2Mass_incorrectmunupair_ev%d_genlevel"%nEv)
	  hme_gen.hme_h2Mass_correctmunupair.SetLineColor(ROOT.kRed)
	  hme_gen.hme_h2Mass_incorrectmunupair.SetLineColor(ROOT.kBlue)
    	  hme_gen.hme_offshellWmass.SetName("hme_offshellWmass_ev%d_genlevle"%nEv)
	  h_h2mass_weight_gen_sum.Add(hme_gen.hme_h2Mass)
          offshellWmass_gen.Add(hme_gen.hme_offshellWmass)
	  if nEv%100 == 0:
	      hme_gen.hme_h2Mass.Write()
	      hme_gen.hme_h2MassWeight1.Write()
    	      hme_gen.hme_h2Mass_correctmunupair.Write()
    	      hme_gen.hme_h2Mass_incorrectmunupair.Write()
    	      hme_gen.hme_offshellWmass.Write()
              hme_gen.hmetree.Write()
          h_h2mass_weight_gen.Fill(hme_gen.hme_h2Mass.GetXaxis().GetBinCenter(hme_gen.hme_h2Mass.GetMaximumBin()))
          h_h2mass_weight1_gen.Fill(hme_gen.hme_h2MassWeight1.GetXaxis().GetBinCenter(hme_gen.hme_h2MassWeight1.GetMaximumBin()))
          h_h2mass_gen[0] = hme_gen.hme_h2Mass.GetXaxis().GetBinCenter(hme_gen.hme_h2Mass.GetMaximumBin())
  
  if (cleaning_cuts and doHME):
      muon1_p4 	  = ROOT.TLorentzVector(ev.muon1_px, ev.muon1_py, ev.muon1_pz, ev.muon1_energy)
      muon2_p4 	  = ROOT.TLorentzVector(ev.muon2_px, ev.muon2_py, ev.muon2_pz, ev.muon2_energy)
      b1jet_p4 	  = ROOT.TLorentzVector(ev.b1jet_px, ev.b1jet_py, ev.b1jet_pz, ev.b1jet_energy)
      b2jet_p4 	  = ROOT.TLorentzVector(ev.b2jet_px, ev.b2jet_py, ev.b2jet_pz, ev.b2jet_energy)
      met_vec2    = ROOT.TVector2(ev.met_px, ev.met_py)
      hme = HeavyMassEstimator()
      hme.setKinematic(muon1_p4, muon2_p4, b1jet_p4, b2jet_p4, met_vec2, 0)
      hme.setonshellWmasspdf(onshellWmasspdf)
      hme.setonshellnuptpdf(onshellnuptpdf)
      hme.setrecobjetrescalec1pdf(recobjetrescalec1pdfPU40)
      #hme.showKinematic()
      hme.setIterations(args.iterations)
      hme.runHME()
      hme.hme_offshellWmass.SetName("hme_offshellWmass_ev%d_genlevle"%nEv)
      if hme.hme_h2Mass.GetEntries()>0:
	  print "Reco Level most probable reco mass ",hme.hme_h2Mass.GetXaxis().GetBinCenter(hme.hme_h2Mass.GetMaximumBin())," entries ",hme.hme_h2Mass.GetEntries()
	  hme.hme_h2Mass.SetName("hme_h2Mass_ev%d_recolevel"%nEv)
          hme.hme_h2MassWeight1.SetName("hme_h2MassWeight1_ev%d_recolevel"%nEv)
          hme.hme_offshellWmass.SetName("hme_offshellWmass_ev%d_recolevle"%nEv)
	  h_h2mass_weight_reco_sum.Add(hme.hme_h2Mass)
	  #if nEv%100 == 0:
	  #    hme.hme_h2Mass.Write()
	  #    hme.hme_h2MassWeight1.Write()
    	  #    hme.hme_offshellWmass.Write()
          h_h2mass_weight_reco.Fill(hme.hme_h2Mass.GetXaxis().GetBinCenter(hme.hme_h2Mass.GetMaximumBin()))
          h_h2mass_weight1_reco.Fill(hme.hme_h2MassWeight1.GetXaxis().GetBinCenter(hme.hme_h2MassWeight1.GetMaximumBin()))
    	  h_h2mass_reco[0] = hme.hme_h2Mass.GetXaxis().GetBinCenter(hme.hme_h2Mass.GetMaximumBin())
  # WEIGHTS
  passCCbr.Fill();   passprebr.Fill();   h_h2mass_genbr.Fill();  h_h2mass_recobr.Fill();
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
h_h2mass_weight_gen.Write()
h_h2mass_weight1_gen.Write()
h_h2mass_weight_reco.Write()
h_h2mass_weight1_reco.Write()
h_h2mass_weight_gen_sum.Write()
h_h2mass_weight_reco_sum.Write()
TCha2.Write()
f.Close()
print "Done."
