import os, sys, random
import ROOT
from ROOT import TFile,TChain,TH1F,TH2F,TLegend
from math import *
from HeavyMassEstimator import *
import argparse
import numpy as np
#sys.argv.append( '-b' ) 
from array import array
import mt2_bisect

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
offshellWmasspdf = refPDF.Get("offshellWmass")
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
parser.add_argument("-it", "--iterations", dest="iterations", type=int, default=10000, help="iteration number used HME [Default: 100000]")
parser.add_argument("-o", "--outputdir",dest="output", default="./",help="output root file directory, [Defualt:PWD] ")
args = parser.parse_args()
if args.output[-1] != '/':
   args.output = args.output+'/'
print " output ",args.output
whichSample = args.jobtype
makeHadd = "HaddNo"

maxn = 10 
sample_weight               = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
jj_M                        = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
lep1_pt                     = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
lep2_pt                     = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
jet1_pt                     = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
jet2_pt                     = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
met_pt                      = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
ht                          = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
llmetjj_M                   = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
cosThetaStar                = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
ll_M                        = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
ll_DR_l_l                   = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
jj_DR_j_j                   = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
llmetjj_DPhi_ll_jj          = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
llmetjj_MTformula           = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
ll_pt                       = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
jj_pt                       = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
llmetjj_minDR_l_j           = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
lep1_eta                    = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
lep1_phi                    = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
lep1_Iso                    = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
lep2_eta                    = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
lep2_phi                    = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
lep2_Iso                    = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
ll_DPhi_l_l                 = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
ll_DEta_l_l                 = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
jet1_eta                    = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
jet1_phi                    = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
jet2_eta                    = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
jet2_phi                    = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
event_weight                = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
event_pu_weight             = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
isElEl                      = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
isMuMu                      = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
isElMu                      = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
isMuEl                      = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
event_number                = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
event_run                   = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
isSF                        = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
total_weight                = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
jet1_cMVAv2                 = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
jet2_cMVAv2                 = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
nJetsL                      = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
llidiso                     = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
mumuidiso                   = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
elelidiso                   = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
jjbtag_heavy                = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
jjbtag_light                = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
trigeff                     = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
pu                          = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
DY_BDT_flat                 = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
dy_nobtag_to_btagM_weight   = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float) 
passCC                      = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
passpre                     = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
mt2                         = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_h2mass_gen              = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_h2mass_reco             = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_mean_reco               = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_stddev_reco             = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_entries_reco            = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_entry_peak_reco         = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_h2mass_weight2_gen      = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_h2mass_weight2_reco     = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_mean_weight2_reco       = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_stddev_weight2_reco     = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_entries_weight2_reco    = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_entry_peak_weight2_reco = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
#hme_offshellWmass_reco      = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_mostprob_offshellWmass_reco = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_mean_offshellWmass_reco = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_stddev_offshellWmass_reco = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_mostprob_offshellWmass_weight2_reco = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_mean_offshellWmass_weight2_reco = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_stddev_offshellWmass_weight2_reco = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hmecputime                  = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
cross_section               = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
event_weight_sum            = array( 'f', maxn*[ 0. ] ) #np.zeros(1, dtype=float)
hme_bins = []
for i in range(100):
    hme_bins.append(array( 'f', maxn*[ 0. ] ))
#print "Executing: python", sys.argv[0] , "-b", sys.argv[2], sys.argv[3], "(Arg1=makeHadd: HaddYes or HaddNo /|\ Arg2=whichSample: TT, DY, VV, sT, Wjet, ttV, Data)"
#makeHadd = sys.argv[2]
#whichSample = sys.argv[3]
if( makeHadd!="HaddYes" and makeHadd!="HaddNo" ): print "WARNING! 1st argument has to be HaddYes or HaddNo"; sys.exit()
#if( whichSample!="TT" and whichSample!="DY" and whichSample!="VV" and whichSample!="sT" and whichSample!="Wjet" and whichSample!="ttV" and whichSample!="Data" ):  print "WARNING! 2nd argument have to be TT, DY, VV, singTop, Wjet, ttV, or Data"; sys.exit()
######################################
## define shell commands
######################################
print "sampletype ",whichSample, " ",Samplelist.Ntuplelist_Louvain[whichSample]
Find_str = []; this_cat = ""; this_hadd = ""; this_NtotPath = ""
for i in range(len(Samplelist.Ntuplelist_Louvain[whichSample])):
  Find_str.append("find %s | grep root | grep -v failed > HADD/%s_%d.txt"%(Samplelist.Ntuplelist_Louvain[whichSample][i], Samplelist.sampleFullName[whichSample], i))
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
  #print "this_find ",this_find
  os.system(this_find)

os.system(this_cat)
#print "this_hadd ", this_hadd
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
      tfiletmp = ROOT.TFile(str(line[:-1]),"READ")
      cross_section[0] = tfiletmp.Get("cross_section").GetVal()
      event_weight_sum[0] = tfiletmp.Get("event_weight_sum").GetVal()
      tfiletmp.Close()
      print "cross_section[0] ",cross_section[0]," event_weight_sum[0]  ",event_weight_sum[0]
      TCha.Add(str(line[:-1]))
#TCha.Add("/home/taohuang/DiHiggsAnalysis/CMSSW_8_0_26_patch2/src/HhhAnalysis/CutFlowAnalyzer/test/out_ana_10k_graviton.root")
#TCha.Add(this_NtotPath)
#TCha.Add("/fdata/hepx/store/user/tahuang/HHNTuples/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-400_narrow_Summer16MiniAODv2_v5.0.1+80X_HHAnalysis_2017-03-01.v0_histos.root")
TotalEv = TCha.GetEntries()
nStart = TotalEv/args.njobs*args.ijob
nEnd = TotalEv/args.njobs*(args.ijob+1)
if args.ijob == args.njobs-1:
    nEnd = -1
print whichSample, "TChain has", TCha.GetEntries(), "entries."," job star from ",nStart," end at ",nEnd
#f = ROOT.TFile('/fdata/hepx/store/user/%s/Hhh_For_Test/'%user + whichSample + '_ijob%d.root'%args.ijob,'recreate'); f.cd()
f = ROOT.TFile(args.output + whichSample + '_ijob%d.root'%args.ijob, 'recreate'); f.cd()
#f = ROOT.TFile('test.root','recreate'); f.cd()
#TCha2 = TCha.Clone()


TCha2 = ROOT.TTree("evtreeHME","tree with HME")
TCha2.Branch("sample_weight",                 sample_weight,                "sample_weight/F")
TCha2.Branch("jj_M",                          jj_M,                         "jj_M/F")
TCha2.Branch("lep1_pt",                       lep1_pt,                      "lep1_pt/F")
TCha2.Branch("lep2_pt",                       lep2_pt,                      "lep2_pt/F")
TCha2.Branch("jet1_pt",                       jet1_pt,                      "jet1_pt/F")
TCha2.Branch("jet2_pt",                       jet2_pt,                      "jet2_pt/F")
TCha2.Branch("met_pt",                        met_pt,                       "met_pt/F")
TCha2.Branch("ht",                            ht,                           "ht/F")
TCha2.Branch("llmetjj_M",                     llmetjj_M,                    "llmetjj_M/F")
TCha2.Branch("cosThetaStar",                  cosThetaStar,                 "cosThetaStar/F")
TCha2.Branch("ll_M",                          ll_M,                         "ll_M/F")
TCha2.Branch("ll_DR_l_l",                     ll_DR_l_l,                    "ll_DR_l_l/F")
TCha2.Branch("jj_DR_j_j",                     jj_DR_j_j,                    "jj_DR_j_j/F")
TCha2.Branch("llmetjj_DPhi_ll_jj",            llmetjj_DPhi_ll_jj,           "llmetjj_DPhi_ll_jj/F")
TCha2.Branch("llmetjj_MTformula",             llmetjj_MTformula,            "llmetjj_MTformula/F")
TCha2.Branch("ll_pt",                         ll_pt,                        "ll_pt/F")
TCha2.Branch("jj_pt",                         jj_pt,                        "ll_pt/F")
TCha2.Branch("llmetjj_minDR_l_j",             llmetjj_minDR_l_j,            "llmetjj_minDR_l_j/F")
TCha2.Branch("lep1_eta",                      lep1_eta,                     "lep1_eta/F")
TCha2.Branch("lep1_phi",                      lep1_phi,                     "lep1_phi/F")
TCha2.Branch("lep1_Iso",                      lep1_Iso,                     "lep1_Iso/F")
TCha2.Branch("lep2_eta",                      lep2_eta,                     "lep2_eta/F")
TCha2.Branch("lep2_phi",                      lep2_phi,                     "lep2_phi/F")
TCha2.Branch("lep2_Iso",                      lep2_Iso,                     "lep2_Iso/F")
TCha2.Branch("ll_DPhi_l_l",                   ll_DPhi_l_l,                  "ll_DPhi_l_l/F")
TCha2.Branch("ll_DEta_l_l",                   ll_DEta_l_l,                  "ll_DEta_l_l/F")
TCha2.Branch("jet1_eta",                      jet1_eta,                     "jet1_eta/F")
TCha2.Branch("jet1_phi",                      jet1_phi,                     "jet1_phi/F")
TCha2.Branch("jet2_eta",                      jet2_eta,                     "jet2_eta/F")
TCha2.Branch("jet2_phi",                      jet2_phi,                     "jet2_phi/F")
TCha2.Branch("event_weight",                  event_weight,                 "event_weight/F")
TCha2.Branch("event_pu_weight",               event_pu_weight,              "event_pu_weight/F")
TCha2.Branch("isElEl",                        isElEl,                       "isElEl/F")
TCha2.Branch("isMuMu",                        isMuMu,                       "isMuMu/F")
TCha2.Branch("isElMu",                        isElMu,                       "isElMu/F")
TCha2.Branch("isMuEl",                        isMuEl,                       "isMuEl/F")
TCha2.Branch("event_number",                  event_number,                 "event_number/F")
TCha2.Branch("event_run",                     event_run,                    "event_run/F")
TCha2.Branch("isSF",                          isSF,                         "isSF/F")
TCha2.Branch("total_weight",                  total_weight,                 "total_weight/F")
TCha2.Branch("jet1_cMVAv2",                   jet1_cMVAv2,                  "jet1_cMVAv2/F")
TCha2.Branch("jet2_cMVAv2",                   jet2_cMVAv2,                  "jet2_cMVAv2/F")
TCha2.Branch("nJetsL",                        nJetsL,                       "nJetsL/F")
TCha2.Branch("llidiso",                       llidiso,                      "llidiso/F")
TCha2.Branch("mumuidiso",                     mumuidiso,                    "mumuidiso/F")
TCha2.Branch("elelidiso",                     elelidiso,                    "elelidiso/F")
TCha2.Branch("jjbtag_heavy",                  jjbtag_heavy,                 "jjbtag_heavy/F")
TCha2.Branch("jjbtag_light",                  jjbtag_light,                 "jjbtag_light/F")
TCha2.Branch("trigeff",                       trigeff,                      "trigeff/F")
TCha2.Branch("pu",                            pu,                           "pu/F")
TCha2.Branch("DY_BDT_flat",                   DY_BDT_flat,                  "DY_BDT_flat/F")
TCha2.Branch("dy_nobtag_to_btagM_weight",     dy_nobtag_to_btagM_weight,    "dy_nobtag_to_btagM_weight/F")
TCha2.Branch("passCC",                        passCC,                       "passCC/F")
TCha2.Branch("passpre",                       passpre,                      "passpre/F")
TCha2.Branch("mt2",                           mt2,                          "mt2/F")
TCha2.Branch("hme_h2mass_gen",                hme_h2mass_gen,               "hme_h2mass_gen/F")
TCha2.Branch("hme_h2mass_reco",               hme_h2mass_reco,              "hme_h2mass_reco/F")
TCha2.Branch("hme_mean_reco",                 hme_mean_reco,                "hme_mean_reco/F")
TCha2.Branch("hme_stddev_reco",               hme_stddev_reco,              "hme_stddev_reco/F")
TCha2.Branch("hme_entries_reco",              hme_entries_reco,             "hme_entries_reco/F")
TCha2.Branch("hme_entry_peak_reco",           hme_entry_peak_reco,          "hme_entry_peak_reco/F")
TCha2.Branch("hme_h2mass_weight2_gen",        hme_h2mass_weight2_gen,       "hme_h2mass_weight2_gen/F")
TCha2.Branch("hme_h2mass_weight2_reco",       hme_h2mass_weight2_reco,      "hme_h2mass_weight2_reco/F")
TCha2.Branch("hme_mean_weight2_reco",         hme_mean_weight2_reco,        "hme_mean_weight2_reco/F")
TCha2.Branch("hme_stddev_weight2_reco",       hme_stddev_weight2_reco,      "hme_stddev_weight2_reco/F")
TCha2.Branch("hme_entries_weight2_reco",      hme_entries_weight2_reco,     "hme_entries_weight2_reco/F")
TCha2.Branch("hme_entry_peak_weight2_reco",   hme_entry_peak_weight2_reco,  "hme_entry_peak_weight2_reco/F")
TCha2.Branch("hme_mostprob_offshellWmass_reco",   hme_mostprob_offshellWmass_reco,  "hme_mostprob_offshellWmass_reco/F")
TCha2.Branch("hme_mean_offshellWmass_reco",   hme_mean_offshellWmass_reco,  "hme_mean_offshellWmass_reco/F")
TCha2.Branch("hme_stddev_offshellWmass_reco", hme_stddev_offshellWmass_reco,"hme_stddev_offshellWmass_reco/F")
TCha2.Branch("hme_mean_offshellWmass_weight2_reco",   hme_mean_offshellWmass_weight2_reco,  "hme_mean_offshellWmass_weight2_reco/F")
TCha2.Branch("hme_mostprob_offshellWmass_weight2_reco",   hme_mostprob_offshellWmass_weight2_reco,  "hme_mostprob_offshellWmass_weight2_reco/F")
TCha2.Branch("hme_stddev_offshellWmass_weight2_reco", hme_stddev_offshellWmass_weight2_reco,"hme_stddev_offshellWmass_weight2_reco/F")
TCha2.Branch("hmecputime",                    hmecputime,                   "hmecputime/F")
TCha2.Branch("cross_section",                 cross_section,                "cross_section/F")
TCha2.Branch("event_weight_sum",              event_weight_sum,             "event_weight_sum/F")

for i in range(100):
    TCha2.Branch("hme_bin%d"%i, hme_bins[i], "hme_bin%d/F"%(i))



ptbinSF= [20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 120.0, 200.0, 300.]
etabinSF = [0, 0.9, 1.2, 2.1, 2.4]
myptbinSF = np.asarray(ptbinSF)
myetabinSF = np.asarray(etabinSF)

# PRESELECTION
# Weights
h_pre_Nev_preHLT         = ROOT.TH1F("h_pre_Nev_preHLT","",1,-0.5,0.5);          h_pre_Nev_preHLT.GetXaxis().SetTitle("#Events (pre HLT)"); h_pre_Nev_preHLT.SetBinContent(1,nTOT_prehlt)
h_pre_Nev_posHLT         = ROOT.TH1F("h_pre_Nev_posHLT","",1,-0.5,0.5);          h_pre_Nev_posHLT.GetXaxis().SetTitle("#Events (post HLT)"); h_pre_Nev_posHLT.SetBinContent(1,nTOT_posthlt)
h_pre_XsecBr             = ROOT.TH1F("h_pre_XsecBr","",1,-0.5,0.5);              h_pre_XsecBr.GetXaxis().SetTitle("XsecBr");
h_pre_lep1_triggerSF    = ROOT.TH2F("h_pre_lep1_triggerSF","",len(etabinSF)-1,myetabinSF,len(ptbinSF)-1,myptbinSF);  h_pre_lep1_triggerSF.GetXaxis().SetTitle("#mu 1 Trigger SF");
h_pre_lep1_isoSF        = ROOT.TH2F("h_pre_lep1_isoSF","",len(etabinSF)-1,myetabinSF,len(ptbinSF)-1,myptbinSF); h_pre_lep1_isoSF.GetXaxis().SetTitle("#mu 1 Iso SF");
h_pre_lep1_idSF         = ROOT.TH2F("h_pre_lep1_idSF","",len(etabinSF)-1,myetabinSF,len(ptbinSF)-1,myptbinSF);  h_pre_lep1_idSF.GetXaxis().SetTitle("#mu 1 ID SF");
h_pre_lep1_trackingSF   = ROOT.TH1F("h_pre_lep1_trackingSF","",50,-2.4,2.4); h_pre_lep1_trackingSF.GetXaxis().SetTitle("#mu 1 tracking SF");
h_pre_lep1_SF_bg2       = ROOT.TH2F("h_pre_lep1_SF_bg2","",len(etabinSF)-1,myetabinSF,len(ptbinSF)-1,myptbinSF);  h_pre_lep1_SF_bg2.GetXaxis().SetTitle("#mu SF bg2");
h_pre_lep1_SF_bg1       = ROOT.TH1F("h_pre_moun1_SF_bg1","",50,-2.4,2.4); h_pre_lep1_SF_bg1.GetXaxis().SetTitle("#mu SF bg1");
h_pre_lep2_triggerSF    = ROOT.TH2F("h_pre_lep2_triggerSF","",len(etabinSF)-1,myetabinSF,len(ptbinSF)-1,myptbinSF);  h_pre_lep2_triggerSF.GetXaxis().SetTitle("#mu 1 Trigger SF");
h_pre_lep2_isoSF        = ROOT.TH2F("h_pre_lep2_isoSF","",len(etabinSF)-1,myetabinSF,len(ptbinSF)-1,myptbinSF);      h_pre_lep2_isoSF.GetXaxis().SetTitle("#mu 1 Iso SF");
h_pre_lep2_idSF         = ROOT.TH2F("h_pre_lep2_idSF","",len(etabinSF)-1,myetabinSF, len(ptbinSF)-1, myptbinSF);       h_pre_lep2_idSF.GetXaxis().SetTitle("#mu 1 ID SF");
h_pre_lep2_trackingSF   = ROOT.TH1F("h_pre_lep2_trackingSF","",50,-2.4,2.4); h_pre_lep2_trackingSF.GetXaxis().SetTitle("#mu 1 tracking SF");
h_pre_lep2_SF_bg2       = ROOT.TH2F("h_pre_lep2_SF_bg2","",len(etabinSF)-1,myetabinSF,len(ptbinSF)-1,myptbinSF);  h_pre_lep2_SF_bg2.GetXaxis().SetTitle("#mu SF bg2");
h_pre_lep2_SF_bg1       = ROOT.TH1F("h_pre_moun2_SF_bg1","",50,-2.4,2.4); h_pre_lep1_SF_bg1.GetXaxis().SetTitle("#mu SF bg1");
# Regression variables
h_pre_numOfVertices      = ROOT.TH1F("h_pre_numOfVertices","",50,0.,50.);      h_pre_numOfVertices.GetXaxis().SetTitle("");
h_pre_jet1_mt           = ROOT.TH1F("h_pre_jet1_mt","",50,0.,100.);          h_pre_jet1_mt.GetXaxis().SetTitle("");
h_pre_jet1_leadTrackPt  = ROOT.TH1F("h_pre_jet1_leadTrackPt","",50,0.,100.); h_pre_jet1_leadTrackPt.GetXaxis().SetTitle("");
h_pre_jet1_leptonPtRel  = ROOT.TH1F("h_pre_jet1_leptonPtRel","",50,0.,10.);  h_pre_jet1_leptonPtRel.GetXaxis().SetTitle("");
h_pre_jet1_leptonPt     = ROOT.TH1F("h_pre_jet1_leptonPt","",50,0.,100.);    h_pre_jet1_leptonPt.GetXaxis().SetTitle("");
h_pre_jet1_leptonDeltaR = ROOT.TH1F("h_pre_jet1_leptonDeltaR","",50,0.,5.);  h_pre_jet1_leptonDeltaR.GetXaxis().SetTitle("");
h_pre_jet1_neHEF        = ROOT.TH1F("h_pre_jet1_neHEF","",50,0.,2.);         h_pre_jet1_neHEF.GetXaxis().SetTitle("");
h_pre_jet1_neEmEF       = ROOT.TH1F("h_pre_jet1_neEmEF","",50,0.,2.);        h_pre_jet1_neEmEF.GetXaxis().SetTitle("");
h_pre_jet1_vtxNtracks   = ROOT.TH1F("h_pre_jet1_vtxNtracks","",50,0.,50.);   h_pre_jet1_vtxNtracks.GetXaxis().SetTitle("");
h_pre_jet1_vtxPt        = ROOT.TH1F("h_pre_jet1_vtxPt","",50,0.,100.);       h_pre_jet1_vtxPt.GetXaxis().SetTitle("");
h_pre_jet1_vtxMass      = ROOT.TH1F("h_pre_jet1_vtxMass","",50,0.,50.);      h_pre_jet1_vtxMass.GetXaxis().SetTitle("");
h_pre_jet1_vtx3DSig     = ROOT.TH1F("h_pre_jet1_vtx3DSig","",50,0.,5.);      h_pre_jet1_vtx3DSig.GetXaxis().SetTitle("");
h_pre_jet1_vtx3DVal     = ROOT.TH1F("h_pre_jet1_vtx3DVal","",50,0.,5.);      h_pre_jet1_vtx3DVal.GetXaxis().SetTitle("");
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
h_h2mass_weight2_gen     = ROOT.TH1F("h_h2mass_weight2_gen","",100,200.,1200.);  h_h2mass_weight2_gen.GetXaxis().SetTitle("HME reco mass [GeV]");
h_h2mass_weight_reco     = ROOT.TH1F("h_h2mass_weight_reco","",100,200.,1200.);  h_h2mass_weight_reco.GetXaxis().SetTitle("HME reco mass [GeV]");
h_h2mass_weight_reco_sum = ROOT.TH1F("h_h2mass_weight_reco_sum","",1000,200.,1200.);  h_h2mass_weight_reco_sum.GetXaxis().SetTitle("HME reco mass [GeV]");
h_h2mass_weight1_reco    = ROOT.TH1F("h_h2mass_weight1_reco","",100,200.,1200.);  h_h2mass_weight1_reco.GetXaxis().SetTitle("HME reco mass [GeV]");
h_h2mass_weight2_reco     = ROOT.TH1F("h_h2mass_weight2_reco","",100,200.,1200.);  h_h2mass_weight2_reco.GetXaxis().SetTitle("HME reco mass [GeV]");
offshellWmass_gen	 = ROOT.TH1F("hme_offshellWmass_gen","offshell W mass from HME", 60, 10.0, 70.0); offshellWmass_gen.GetXaxis().SetTitle("HME offshell W mass [GeV]");
offshellWmass_reco	 = ROOT.TH1F("hme_offshellWmass_reco","offshell W mass from HME", 60, 10.0, 70.0); offshellWmass_reco.GetXaxis().SetTitle("HME offshell W mass [GeV]");
# adding more branches to TTree
def initbr():
    passCC[0] = 0
    passpre[0] = 0
    hme_h2mass_gen[0] = -1
    hme_h2mass_reco[0] = -1
    #print "hme_bins ",hme_bins
    hme_entries_reco[0] = 0.0
    for i in range(len(hme_bins)):
	hme_bins[i][0] = 0.0

stop_watch = ROOT.TStopwatch()
stop_watch2 = ROOT.TStopwatch()

stop_watch2.Start()
nEv = 0
for ev in TCha:
  if (doTest and nEv%10 == 0 ):
      print "nEv ",nEv
  elif (nEv%10000 == 0):
      print "nEv ",nEv
  if nEv < nStart:
      #print "nEv ",nEv," nstart ",nStart," nEnd ",nEnd
      nEv = nEv +1
      continue
  if (nEnd > 0 and nEv >= nEnd):
      break
  
  if (doTest and nEv>=100):
      break

  if nEv == nStart:
      print "First event to run: ",nEv, " last event ", nEnd," HME iterations ",args.iterations
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
  if whichSample.startswith("DY"):
      DY_BDT_flat[0]    =  ev.DY_BDT_flat
      dy_nobtag_to_btagM_weight[0] = ev.dy_nobtag_to_btagM_weight
  MET_cut               = (ev.met_pt>20)
  MuMu_cut              = (ev.lep1_pt>20 and fabs(ev.lep1_eta)<2.4 and ev.lep2_pt>10 and fabs(ev.lep2_eta)<2.4 and ev.ll_M>12)
  #B1B2_cut              = (ev.jet1_pt>20 and fabs(ev.jet1_eta)<2.4 and ev.jet2_pt>20 and fabs(ev.jet2_eta)<2.4)
  B1B2_cut              = (ev.jet1_pt>20 and ev.jet2_pt>20)
  preselection          = (MET_cut and MuMu_cut and B1B2_cut)
  mt_cleancut           = (ev.llmetjj_MTformula>10)
  NoZ_cleancut          = (ev.ll_M<(91-15))
  mJJ_cleancut          = (ev.jj_M>30)
  #dR_lljj_cleancut      = (ev.dR_l1l2b1b2>0.2 and ev.dR_l1l2b1b2<4)
  dR_lljj_cleancut      = 1
  dR_ll_cleancut        = (ev.ll_DR_l_l<3.8)
  dR_jj_cleancut        = (ev.jj_DR_j_j<3.8)
  cleaning_cuts         = (mt_cleancut and NoZ_cleancut and mJJ_cleancut and dR_lljj_cleancut and dR_ll_cleancut and dR_jj_cleancut)
  #if (preselection and cleaning_cuts ):
  #    print "event ",nEv," pass CC cut"
  if (preselection):
      passpre[0] = 1
  if (cleaning_cuts):
      passCC[0] = 1
  

  #if (cleaning_cuts and ev.findAllGenParticles and doHME and False):
  if (False):
      
      mu1_p4 	      = ROOT.TLorentzVector(ev.mu1_px, ev.mu1_py, ev.mu1_pz, ev.mu1_energy)
      mu2_p4 	      = ROOT.TLorentzVector(ev.mu2_px, ev.mu2_py, ev.mu2_pz, ev.mu2_energy)
      jet1_p4 	      = ROOT.TLorentzVector(ev.b1_px, ev.b1_py, ev.b1_pz, ev.b1_energy)
      jet2_p4 	      = ROOT.TLorentzVector(ev.b2_px, ev.b2_py, ev.b2_pz, ev.b2_energy)
      met_vec2        = ROOT.TVector2(ev.nu1_px + ev.nu2_px, ev.nu1_py + ev.nu2_py)
      hme_gen = HeavyMassEstimator()
      onshellW_mu  = 1 
      if ev.w1_mass < ev.w2_mass:
      	  onshellW_mu  = 2
      hme_gen.setKinematic(mu1_p4, mu2_p4, jet1_p4, jet2_p4, met_vec2, onshellW_mu)
      hme_gen.setonshellWmasspdf(onshellWmasspdf)
      hme_gen.setoffshellWmasspdf(offshellWmasspdf)
      hme_gen.setonshellnuptpdf(onshellnuptpdf)
      #hme.showKinematic()
      hme_gen.setIterations(args.iterations)
      hme_gen.runHME()
      if hme_gen.hme_h2Mass.GetEntries()>0 and hme_gen.hme_h2Mass.GetXaxis().GetBinCenter(hme_gen.hme_h2Mass.GetMaximumBin())>200:
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
          hme_h2mass_gen[0] = hme_gen.hme_h2Mass.GetXaxis().GetBinCenter(hme_gen.hme_h2Mass.GetMaximumBin())
  
  #if (cleaning_cuts and doHME):
  stop_watch.Start()
  if (doHME):
      #lep1_p4 	  = ROOT.TLorentzVector(ev.lep1_px, ev.lep1_py, ev.lep1_pz, ev.lep1_energy)
      #lep2_p4 	  = ROOT.TLorentzVector(ev.lep2_px, ev.lep2_py, ev.lep2_pz, ev.lep2_energy)
      #jet1_p4 	  = ROOT.TLorentzVector(ev.jet1_px, ev.jet1_py, ev.jet1_pz, ev.jet1_energy)
      #jet2_p4 	  = ROOT.TLorentzVector(ev.jet2_px, ev.jet2_py, ev.jet2_pz, ev.jet2_energy)
      #met_vec2    = ROOT.TVector2(ev.met_px, ev.met_py)
      lep1_p4 	  = ROOT.TLorentzVector(); lep1_p4.SetPtEtaPhiM(ev.lep1_pt, ev.lep1_eta, ev.lep1_phi, 0)
      lep2_p4 	  = ROOT.TLorentzVector(); lep2_p4.SetPtEtaPhiM(ev.lep2_pt, ev.lep2_eta, ev.lep2_phi, 0) 
      jet1_p4 	  = ROOT.TLorentzVector(); jet1_p4.SetPtEtaPhiM(ev.jet1_pt, ev.jet1_eta, ev.jet1_phi, 4.18)
      jet2_p4 	  = ROOT.TLorentzVector(); jet2_p4.SetPtEtaPhiM(ev.jet2_pt, ev.jet2_eta, ev.jet2_phi, 4.18)
      met_vec2    = ROOT.TVector2();       met_vec2.SetMagPhi(ev.met_pt, ev.met_phi)
      
      mt2_event1 =  mt2_bisect.mt2()
      sumesBl11_p4 = lep1_p4+jet1_p4
      sumesBl22_p4 = lep2_p4+jet2_p4
      mt2_event1.set_momenta(sumesBl11_p4.M(), sumesBl11_p4.Px(), sumesBl11_p4.Py(), sumesBl22_p4.M(), sumesBl22_p4.Px(), sumesBl22_p4.Py(), 0.0, met_vec2.Px(), met_vec2.Py())
      mt2_event1.set_mn(0.0)
      MT2_1 = mt2_event1.get_mt2()

      mt2_event2 =  mt2_bisect.mt2()
      sumesBl12_p4 = lep2_p4+jet1_p4
      sumesBl21_p4 = lep1_p4+jet2_p4
      mt2_event2.set_momenta(sumesBl12_p4.M(), sumesBl12_p4.Px(), sumesBl12_p4.Py(), sumesBl21_p4.M(), sumesBl21_p4.Px(), sumesBl21_p4.Py(), 0.0, met_vec2.Px(), met_vec2.Py())
      mt2_event2.set_mn(0.0)
      MT2_2 = mt2_event2.get_mt2()
      
      mt2 = min(MT2_1, MT2_2)
      hme = HeavyMassEstimator()
      hme.setKinematic(lep1_p4, lep2_p4, jet1_p4, jet2_p4, met_vec2, 0)
      hme.setonshellWmasspdf(onshellWmasspdf)
      hme.setoffshellWmasspdf(offshellWmasspdf)
      hme.setonshellnuptpdf(onshellnuptpdf)
      hme.setrecobjetrescalec1pdf(recobjetrescalec1pdfPU40)
      #hme.showKinematic()
      hme.setIterations(args.iterations)
      hme.runHME()
      #hme.hme_offshellWmass.SetName("hme_offshellWmass_ev%d_genlevle"%nEv)
      if hme.hme_h2Mass.GetEntries()>0 and hme.hme_h2Mass.GetXaxis().GetBinCenter(hme.hme_h2Mass.GetMaximumBin())>=250.0:
	  #print "Reco Level most probable reco mass ",hme.hme_h2Mass.GetXaxis().GetBinCenter(hme.hme_h2Mass.GetMaximumBin())," entries ",hme.hme_h2Mass.GetEntries()
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
          h_h2mass_weight2_reco.Fill(hme.hme_h2MassWeight2.GetXaxis().GetBinCenter(hme.hme_h2MassWeight2.GetMaximumBin()))

    	  hme_h2mass_reco[0] = hme.hme_h2Mass.GetXaxis().GetBinCenter(hme.hme_h2Mass.GetMaximumBin())
          hme_mean_reco[0] = hme.hme_h2Mass.GetMean()
          hme_stddev_reco[0] = hme.hme_h2Mass.GetStdDev(1)
          hme_entries_reco[0] = float(hme.hme_h2Mass.GetEntries())/args.iterations
	  hme_entry_peak_reco[0] = hme.hme_h2Mass.Integral(hme.hme_h2Mass.GetMaximumBin()-5, hme.hme_h2Mass.GetMaximumBin()+5)
    	  hme_h2mass_weight2_reco[0] = hme.hme_h2MassWeight2.GetXaxis().GetBinCenter(hme.hme_h2MassWeight2.GetMaximumBin())
          hme_mean_weight2_reco[0] = hme.hme_h2MassWeight2.GetMean()
          hme_stddev_weight2_reco[0] = hme.hme_h2MassWeight2.GetStdDev(1)
          hme_entries_weight2_reco[0] = float(hme.hme_h2MassWeight2.GetEntries())/args.iterations
	  hme_entry_peak_weight2_reco[0] = hme.hme_h2MassWeight2.Integral(hme.hme_h2MassWeight2.GetMaximumBin()-5, hme.hme_h2MassWeight2.GetMaximumBin()+5)
	  print "hme_h2mass_reco[0] ",hme_h2mass_reco[0]," hme_h2mass_weight2_reco[0] ",hme_h2mass_weight2_reco[0]
	  #offshell Wmass
	  h_offshellWmass_recoh2mass = hme.hme_h2MassAndoffshellWmass.ProjectionY("h_offshellWmass_recoh2mass",hme.hme_h2Mass.GetMaximumBin()-5, hme.hme_h2Mass.GetMaximumBin()+5)
	  h_offshellWmass_recoh2mass_weight2 = hme.hme_h2MassAndoffshellWmass_weight2.ProjectionY("h_offshellWmass_recoh2mass_weight2",hme.hme_h2MassWeight2.GetMaximumBin()-5, hme.hme_h2MassWeight2.GetMaximumBin()+5)
	  hme_mostprob_offshellWmass_reco[0] = h_offshellWmass_recoh2mass.GetXaxis().GetBinCenter(h_offshellWmass_recoh2mass.GetMaximumBin())
	  hme_mean_offshellWmass_reco[0] = h_offshellWmass_recoh2mass.GetMean()
	  hme_stddev_offshellWmass_reco[0] = h_offshellWmass_recoh2mass.GetStdDev(1)
	  hme_mostprob_offshellWmass_weight2_reco[0] = h_offshellWmass_recoh2mass_weight2.GetXaxis().GetBinCenter(h_offshellWmass_recoh2mass_weight2.GetMaximumBin())
	  hme_mean_offshellWmass_weight2_reco[0] = h_offshellWmass_recoh2mass_weight2.GetMean()
	  hme_stddev_offshellWmass_weight2_reco[0] = h_offshellWmass_recoh2mass_weight2.GetStdDev(1)

          hme.hme_h2Mass.Rebin(10)
          hme.hme_h2MassWeight2.Rebin(10)
          nbins = hme.hme_h2Mass.GetXaxis().GetNbins()
          #print "nbins ",nbins
	  if nbins != len(hme_bins):
	      print "error!!! nbins of hme_h2mass is not the same as the one written to TTree "
	  for i in range(1, len(hme_bins)+1):
	      hme_bins[i-1][0] =  hme.hme_h2Mass.GetBinContent(i)
        
  # WEIGHTS
  #passCCbr.Fill();   passprebr.Fill();   hme_h2mass_genbr.Fill();  hme_h2mass_recobr.Fill();
  hmecputime[0] =  stop_watch.CpuTime()
  #print "Cputime ",stop_watch.CpuTime()," realtime ",stop_watch.RealTime()
  stop_watch.Stop()
  TCha2.Fill()
  #weight = 1.
  #if( whichSample!="Data" ):
  #  # Cross section (We only weights to 1fb-1, the real normalization for xsec will be placed into FinalPlotter.py)
  #  xSec = ev.XsecBr
  #  if(nEv==0): h_pre_XsecBr.SetBinContent(1, xSec)
  #  Lumi    = 1. * 1000 # Convert from 1fb-1 to 1pb-1
  #  weight = weight * (Lumi*xSec)/nTOT_prehlt
  #  # Scale factors
  #  weight = weight * ev.lep1_pogSF * ev.lep2_pogSF
  ## Minimal Selection
  #if( preselection ):
  #  # SF
  #  h_pre_lep1_triggerSF.Fill(abs(ev.lep1_eta), ev.lep1_pt, ev.lep1_triggerSF)
  #  h_pre_lep1_isoSF.Fill(abs(ev.lep1_eta), ev.lep1_pt, ev.lep1_isoSF)
  #  h_pre_lep1_idSF.Fill(abs(ev.lep1_eta), ev.lep1_pt, ev.lep1_idSF)
  #  h_pre_lep1_trackingSF.Fill(ev.lep1_eta, ev.lep1_trackingSF)
  #  h_pre_lep1_SF_bg1.Fill(ev.lep1_eta)
  #  h_pre_lep1_SF_bg2.Fill(abs(ev.lep1_eta), ev.lep1_pt)
  #  h_pre_lep2_triggerSF.Fill(abs(ev.lep2_eta), ev.lep2_pt, ev.lep2_triggerSF)
  #  h_pre_lep2_isoSF.Fill(abs(ev.lep2_eta), ev.lep2_pt, ev.lep2_isoSF)
  #  h_pre_lep2_idSF.Fill(abs(ev.lep2_eta), ev.lep2_pt, ev.lep2_idSF)
  #  h_pre_lep2_trackingSF.Fill(ev.lep2_eta, ev.lep2_trackingSF)
  #  h_pre_lep2_SF_bg1.Fill(ev.lep2_eta)
  #  h_pre_lep2_SF_bg2.Fill(abs(ev.lep2_eta), ev.lep2_pt)
  #  # Regression Variables
  #  h_pre_numOfVertices.Fill( ev.numOfVertices, weight )   
  #  h_pre_jet1_mt.Fill( ev.jet1_mt, weight )
  #  h_pre_jet1_leadTrackPt.Fill( ev.jet1_leadTrackPt, weight )
  #  h_pre_jet1_leptonPtRel.Fill( ev.jet1_leptonPtRel, weight )
  #  h_pre_jet1_leptonPt.Fill( ev.jet1_leptonPt, weight )
  #  h_pre_jet1_leptonDeltaR.Fill( ev.jet1_leptonDeltaR, weight )
  #  h_pre_jet1_neHEF.Fill( ev.jet1_neHEF, weight )
  #  h_pre_jet1_neEmEF.Fill( ev.jet1_neEmEF, weight )
  #  h_pre_jet1_vtxNtracks.Fill( ev.jet1_vtxNtracks, weight )
  #  h_pre_jet1_vtxPt.Fill( ev.jet1_vtxPt, weight )
  #  h_pre_jet1_vtxMass.Fill( ev.jet1_vtxMass, weight )
  #  h_pre_jet1_vtx3DSig.Fill( ev.jet1_vtx3DSig, weight )
  #  h_pre_jet1_vtx3DVal.Fill( ev.jet1_vtx3DVal, weight )
  #  # Kinematic Variables
  #  h_pre_MU1_pt.Fill( ev.lep1_pt, weight )
  #  h_pre_MU2_pt.Fill( ev.lep2_pt, weight )
  #  h_pre_MU2_eta.Fill( ev.lep2_eta, weight )
  #  h_pre_MU1_eta.Fill( ev.lep1_eta, weight )
  #  h_pre_mass_l1l2.Fill( ev.mass_l1l2, weight )
  #  h_pre_dR_l1l2.Fill( ev.dR_l1l2, weight )
  #  h_pre_J1_pt.Fill( ev.jet1_pt, weight )
  #  h_pre_J2_pt.Fill( ev.jet2_pt, weight )
  #  h_pre_J1_eta.Fill( ev.jet1_eta, weight )
  #  h_pre_J2_eta.Fill( ev.jet2_eta, weight )
  #  h_pre_mass_b1b2.Fill( ev.mass_b1b2, weight )
  #  h_pre_dR_b1b2.Fill( ev.dR_b1b2, weight )
  #  h_pre_met_pt.Fill( ev.met_pt, weight )
  #  h_pre_mass_trans.Fill( ev.mass_trans, weight )
  #  h_pre_dR_l1l2b1b2.Fill( ev.dR_l1l2b1b2, weight )
  #  h_pre_dphi_llmet.Fill( ev.dphi_llmet, weight )
  #  if( cleaning_cuts ):
  #    # Kinematic Variables
  #    h_cc_MU1_pt.Fill( ev.lep1_pt, weight )
  #    h_cc_MU2_pt.Fill( ev.lep2_pt, weight )
  #    h_cc_MU2_eta.Fill( ev.lep2_eta, weight )
  #    h_cc_MU1_eta.Fill( ev.lep1_eta, weight )
  #    h_cc_mass_l1l2.Fill( ev.mass_l1l2, weight )
  #    h_cc_dR_l1l2.Fill( ev.dR_l1l2, weight )
  #    h_cc_J1_pt.Fill( ev.jet1_pt, weight )
  #    h_cc_J2_pt.Fill( ev.jet2_pt, weight )
  #    h_cc_J1_eta.Fill( ev.jet1_eta, weight )
  #    h_cc_J2_eta.Fill( ev.jet2_eta, weight )
  #    h_cc_mass_b1b2.Fill( ev.mass_b1b2, weight )
  #    h_cc_dR_b1b2.Fill( ev.dR_b1b2, weight )
  #    h_cc_met_pt.Fill( ev.met_pt, weight )
  #    h_cc_mass_trans.Fill( ev.mass_trans, weight )
  #    h_cc_dR_l1l2b1b2.Fill( ev.dR_l1l2b1b2, weight )
  #    h_cc_dphi_llmet.Fill( ev.dphi_llmet, weight )
  nEv = nEv +1

#h_pre_lep1_triggerSF.Divide(h_pre_lep1_SF_bg2)
#h_pre_lep1_isoSF.Divide(h_pre_lep1_SF_bg2)
#h_pre_lep1_idSF.Divide(h_pre_lep1_SF_bg2)
#h_pre_lep1_trackingSF.Divide(h_pre_lep1_SF_bg1)
#h_pre_lep2_triggerSF.Divide(h_pre_lep2_SF_bg2)
#h_pre_lep2_isoSF.Divide(h_pre_lep2_SF_bg2)
#h_pre_lep2_idSF.Divide(h_pre_lep2_SF_bg2)
#h_pre_lep2_trackingSF.Divide(h_pre_lep2_SF_bg1)
#normalize1D(h_pre_lep1_trackingSF)
# Writing histograms
f.cd()
#h_pre_XsecBr.Write()
#h_pre_lep1_triggerSF.Write()
#h_pre_lep1_isoSF.Write()
#h_pre_lep1_idSF.Write()
#h_pre_lep1_trackingSF.Write()
#h_pre_lep2_triggerSF.Write()
#h_pre_lep2_isoSF.Write()
#h_pre_lep2_idSF.Write()
#h_pre_lep2_trackingSF.Write()
#h_pre_Nev_preHLT.Write()
#h_pre_Nev_posHLT.Write()
#h_pre_numOfVertices.Write()  
#h_pre_jet1_mt.Write()
#h_pre_jet1_leadTrackPt.Write()
#h_pre_jet1_leptonPtRel.Write()
#h_pre_jet1_leptonPt.Write()
#h_pre_jet1_leptonDeltaR.Write()
#h_pre_jet1_neHEF.Write()
#h_pre_jet1_neEmEF.Write()
#h_pre_jet1_vtxNtracks.Write()
#h_pre_jet1_vtxPt.Write()
#h_pre_jet1_vtxMass.Write()
#h_pre_jet1_vtx3DSig.Write()
#h_pre_jet1_vtx3DVal.Write()
#h_pre_MU1_pt.Write()
#h_pre_MU2_pt.Write()
#h_pre_MU1_eta.Write()
#h_pre_MU2_eta.Write()
#h_pre_mass_l1l2.Write()
#h_pre_dR_l1l2.Write()
#h_pre_J1_pt.Write()
#h_pre_J2_pt.Write()
#h_pre_J1_eta.Write()
#h_pre_J2_eta.Write()
#h_pre_mass_b1b2.Write()
#h_pre_dR_b1b2.Write()
#h_pre_met_pt.Write()
#h_pre_mass_trans.Write()
#h_pre_dR_l1l2b1b2.Write()
#h_pre_dphi_llmet.Write()
## Cleaning Cuts
#h_cc_MU1_pt.Write()
#h_cc_MU2_pt.Write()
#h_cc_MU2_eta.Write()
#h_cc_MU1_eta.Write()
#h_cc_mass_l1l2.Write()
#h_cc_dR_l1l2.Write()
#h_cc_J1_pt.Write()
#h_cc_J2_pt.Write()
#h_cc_J1_eta.Write()
#h_cc_J2_eta.Write()
#h_cc_mass_b1b2.Write()
#h_cc_dR_b1b2.Write()
#h_cc_met_pt.Write()
#h_cc_mass_trans.Write()
#h_cc_dR_l1l2b1b2.Write()
#h_cc_dphi_llmet.Write()
#h_h2mass_weight_gen.Write()
#h_h2mass_weight1_gen.Write()
#h_h2mass_weight_reco.Write()
#h_h2mass_weight1_reco.Write()
#h_h2mass_weight_gen_sum.Write()
#h_h2mass_weight_reco_sum.Write()
TCha2.Write()
f.Close()
print "Done."
stop_watch2.Stop()
print "stop_watch2 Cputime ",stop_watch2.CpuTime()," realtime ",stop_watch2.RealTime()
