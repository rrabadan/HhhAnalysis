import os, sys, random
import ROOT
from ROOT import TFile,TChain,TH1F,TH2F,TLegend
from math import *
from HeavyMassEstimator import *
import numpy as np
execfile("start.py")
execfile("functions.py")
#Creating folders and parameters

doTest = False
tree_name="DiHiggsWWBBAna/evtree"
#benchmarks = ["ttV", "Wjet", "singTop", "VV", "DY", "TTbar", "Data"]
print "Executing: python", sys.argv[0] , "-b", sys.argv[2], sys.argv[3], "(Arg1=makeHadd: HaddYes or HaddNo /|\ Arg2=whichSample: TT, DY, VV, sT, Wjet, ttV, Data)"
makeHadd = sys.argv[2]
whichSample = sys.argv[3]
if( makeHadd!="HaddYes" and makeHadd!="HaddNo" ): print "WARNING! 1st argument has to be HaddYes or HaddNo"; sys.exit()
if( whichSample!="TT" and whichSample!="DY" and whichSample!="VV" and whichSample!="sT" and whichSample!="Wjet" and whichSample!="ttV" and whichSample!="Data" and not "Rad_" in whichSample ):  print "WARNING! 2nd argument have to be TT, DY, VV, singTop, Wjet, ttV, or Data"; sys.exit()

Find_str = []; this_cat = ""; this_hadd = ""; this_NtotPath = ""
# MC
if( whichSample == "Rad_260" ):
  Find_str.append("find /fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-260_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-260_narrow_13TeV-madgraph-v2/170420_232954 | grep root | grep -v failed > HADD/Rad_260_GluGluToRadionToHHTo2B2VTo2L2NuM-260.txt")
  this_cat      = "cat HADD/Rad_260_* > HADD/Rad_260.txt"
  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-260_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-260_narrow_13TeV-madgraph-v2.root @HADD/Rad_260.txt"
  this_NtotPath = "/fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-260_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-260_narrow_13TeV-madgraph-v2.root"
#if( whichSample == "Rad_270" ):
#  Find_str.append("find /fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-270_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-270_narrow_13TeV-madgraph-v2/170420_233008 | grep root | grep -v failed > HADD/Rad_270_GluGluToRadionToHHTo2B2VTo2L2NuM-270.txt")
#  this_cat      = "cat HADD/Rad_270_* > HADD/Rad_270.txt"
#  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-270_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-270_narrow_13TeV-madgraph-v2.root @HADD/Rad_270.txt"
#  this_NtotPath = "/fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-270_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-270_narrow_13TeV-madgraph-v2.root"
#if( whichSample == "Rad_300" ):
#  Find_str.append("find /fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-300_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-300_narrow_13TeV-madgraph-v2/170420_233022 | grep root | grep -v failed > HADD/Rad_300_GluGluToRadionToHHTo2B2VTo2L2NuM-300.txt")
#  this_cat      = "cat HADD/Rad_300_* > HADD/Rad_300.txt"
#  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-300_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-300_narrow_13TeV-madgraph-v2.root @HADD/Rad_300.txt"
#  this_NtotPath = "/fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-300_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-300_narrow_13TeV-madgraph-v2.root"
#if( whichSample == "Rad_350" ):
#  Find_str.append("find /fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-350_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-350_narrow_13TeV-madgraph-v2/170420_233035 | grep root | grep -v failed > HADD/Rad_350_GluGluToRadionToHHTo2B2VTo2L2NuM-350.txt")
#  this_cat      = "cat HADD/Rad_350_* > HADD/Rad_350.txt"
#  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-350_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-350_narrow_13TeV-madgraph-v2.root @HADD/Rad_350.txt"
#  this_NtotPath = "/fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-350_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-350_narrow_13TeV-madgraph-v2.root"
#if( whichSample == "Rad_400" ):
#  Find_str.append("find /fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-400_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-400_narrow_13TeV-madgraph-v2/170420_233048 | grep root | grep -v failed > HADD/Rad_400_GluGluToRadionToHHTo2B2VTo2L2NuM-400.txt")
#  this_cat      = "cat HADD/Rad_400_* > HADD/Rad_400.txt"
#  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-400_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-400_narrow_13TeV-madgraph-v2.root @HADD/Rad_400.txt"
#  this_NtotPath = "/fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-400_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-400_narrow_13TeV-madgraph-v2.root"
#if( whichSample == "Rad_450" ):
#  Find_str.append("find /fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-450_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-450_narrow_13TeV-madgraph-v2/170420_233102 | grep root | grep -v failed > HADD/Rad_450_GluGluToRadionToHHTo2B2VTo2L2NuM-450.txt")
#  this_cat      = "cat HADD/Rad_450_* > HADD/Rad_450.txt"
#  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-450_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-450_narrow_13TeV-madgraph-v2.root @HADD/Rad_450.txt"
#  this_NtotPath = "/fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-450_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-450_narrow_13TeV-madgraph-v2.root"
if( whichSample == "Rad_500" ):
  Find_str.append("find /fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-500_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-500_narrow_13TeV-madgraph-v2/170420_233115 | grep root | grep -v failed > HADD/Rad_500_GluGluToRadionToHHTo2B2VTo2L2NuM-500.txt")
  this_cat      = "cat HADD/Rad_500_* > HADD/Rad_500.txt"
  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-500_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-500_narrow_13TeV-madgraph-v2.root @HADD/Rad_500.txt"
  this_NtotPath = "/fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-500_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-500_narrow_13TeV-madgraph-v2.root"
#if( whichSample == "Rad_550" ):
#  Find_str.append("find /fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-550_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-550_narrow_13TeV-madgraph-v2/170420_233129 | grep root | grep -v failed > HADD/Rad_550_GluGluToRadionToHHTo2B2VTo2L2NuM-550.txt")
#  this_cat      = "cat HADD/Rad_550_* > HADD/Rad_550.txt"
#  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-550_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-550_narrow_13TeV-madgraph-v2.root @HADD/Rad_550.txt"
#  this_NtotPath = "/fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-550_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-550_narrow_13TeV-madgraph-v2.root"
#if( whichSample == "Rad_600" ):
#  Find_str.append("find /fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-600_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-600_narrow_13TeV-madgraph-v2/170420_233142 | grep root | grep -v failed > HADD/Rad_600_GluGluToRadionToHHTo2B2VTo2L2NuM-600.txt")
#  this_cat      = "cat HADD/Rad_600_* > HADD/Rad_600.txt"
#  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-600_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-600_narrow_13TeV-madgraph-v2.root @HADD/Rad_600.txt"
#  this_NtotPath = "/fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-600_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-600_narrow_13TeV-madgraph-v2.root"
#if( whichSample == "Rad_650" ):
#  Find_str.append("find /fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-650_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-650_narrow_13TeV-madgraph-v2/170420_233156 | grep root | grep -v failed > HADD/Rad_650_GluGluToRadionToHHTo2B2VTo2L2NuM-650.txt")
#  this_cat      = "cat HADD/Rad_650_* > HADD/Rad_650.txt"
#  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-650_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-650_narrow_13TeV-madgraph-v2.root @HADD/Rad_650.txt"
#  this_NtotPath = "/fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-650_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-650_narrow_13TeV-madgraph-v2.root"
#if( whichSample == "Rad_750" ):
#  Find_str.append("find /fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-750_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-750_narrow_13TeV-madgraph-v2/170420_233209 | grep root | grep -v failed > HADD/Rad_750_GluGluToRadionToHHTo2B2VTo2L2NuM-750.txt")
#  this_cat      = "cat HADD/Rad_750_* > HADD/Rad_750.txt"
#  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-750_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-750_narrow_13TeV-madgraph-v2.root @HADD/Rad_750.txt"
#  this_NtotPath = "/fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-750_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-750_narrow_13TeV-madgraph-v2.root"
#if( whichSample == "Rad_800" ):
#  Find_str.append("find /fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-800_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-800_narrow_13TeV-madgraph-v2/170420_233222 | grep root | grep -v failed > HADD/Rad_800_GluGluToRadionToHHTo2B2VTo2L2NuM-800.txt")
#  this_cat      = "cat HADD/Rad_800_* > HADD/Rad_800.txt"
#  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-800_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-800_narrow_13TeV-madgraph-v2.root @HADD/Rad_800.txt"
#  this_NtotPath = "/fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-800_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-800_narrow_13TeV-madgraph-v2.root"
if( whichSample == "Rad_900" ):
  Find_str.append("find /fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-900_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-900_narrow_13TeV-madgraph-v2/170420_233235 | grep root | grep -v failed > HADD/Rad_900_GluGluToRadionToHHTo2B2VTo2L2NuM-900.txt")
  this_cat      = "cat HADD/Rad_900_* > HADD/Rad_900.txt"
  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-900_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-900_narrow_13TeV-madgraph-v2.root @HADD/Rad_900.txt"
  this_NtotPath = "/fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-900_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-900_narrow_13TeV-madgraph-v2.root"
#if( whichSample == "Rad_1000" ):
#  Find_str.append("find /fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-1000_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-1000_narrow_13TeV-madgraph-v2/170420_233248 | grep root | grep -v failed > HADD/Rad_1000_GluGluToRadionToHHTo2B2VTo2L2NuM-1000.txt")
#  this_cat      = "cat HADD/Rad_1000_* > HADD/Rad_1000.txt"
#  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-1000_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-1000_narrow_13TeV-madgraph-v2.root @HADD/Rad_1000.txt"
#  this_NtotPath = "/fdata/hepx/store/user/lpernie/GluGluToRadionToHHTo2B2VTo2L2Nu_M-1000_narrow_13TeV-madgraph-v2/crab_GluGluToRadionToHHTo2B2VTo2L2Nu_M-1000_narrow_13TeV-madgraph-v2.root"
#if( whichSample == "Grav_260" ):
#  Find_str.append("find /fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-260_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-260_narrow_13TeV-madgraph-v2/170420_233302 | grep root | grep -v failed > HADD/Grav_260_GluGluToBulkGravitonToHHTo2B2VTo2L2NuM-260.txt")
#  this_cat      = "cat HADD/Grav_260_* > HADD/Grav_260.txt"
#  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-260_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-260_narrow_13TeV-madgraph-v2.root @HADD/Grav_260.txt"
#  this_NtotPath = "/fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-260_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-260_narrow_13TeV-madgraph-v2.root"
#if( whichSample == "Grav_270" ):
#  Find_str.append("find /fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-270_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-270_narrow_13TeV-madgraph-v2/170420_233315 | grep root | grep -v failed > HADD/Grav_270_GluGluToBulkGravitonToHHTo2B2VTo2L2NuM-270.txt")
#  this_cat      = "cat HADD/Grav_270_* > HADD/Grav_270.txt"
#  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-270_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-270_narrow_13TeV-madgraph-v2.root @HADD/Grav_270.txt"
#  this_NtotPath = "/fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-270_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-270_narrow_13TeV-madgraph-v2.root"
#if( whichSample == "Grav_300" ):
#  Find_str.append("find /fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-300_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-300_narrow_13TeV-madgraph-v2/170420_233329 | grep root | grep -v failed > HADD/Grav_300_GluGluToBulkGravitonToHHTo2B2VTo2L2NuM-300.txt")
#  this_cat      = "cat HADD/Grav_300_* > HADD/Grav_300.txt"
#  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-300_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-300_narrow_13TeV-madgraph-v2.root @HADD/Grav_300.txt"
#  this_NtotPath = "/fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-300_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-300_narrow_13TeV-madgraph-v2.root"
#if( whichSample == "Grav_350" ):
#  Find_str.append("find /fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-350_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-350_narrow_13TeV-madgraph-v2/170420_233342 | grep root | grep -v failed > HADD/Grav_350_GluGluToBulkGravitonToHHTo2B2VTo2L2NuM-350.txt")
#  this_cat      = "cat HADD/Grav_350_* > HADD/Grav_350.txt"
#  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-350_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-350_narrow_13TeV-madgraph-v2.root @HADD/Grav_350.txt"
#  this_NtotPath = "/fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-350_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-350_narrow_13TeV-madgraph-v2.root"
#if( whichSample == "Grav_400" ):
#  Find_str.append("find /fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-400_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-400_narrow_13TeV-madgraph-v2/170420_233355 | grep root | grep -v failed > HADD/Grav_400_GluGluToBulkGravitonToHHTo2B2VTo2L2NuM-400.txt")
#  this_cat      = "cat HADD/Grav_400_* > HADD/Grav_400.txt"
#  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-400_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-400_narrow_13TeV-madgraph-v2.root @HADD/Grav_400.txt"
#  this_NtotPath = "/fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-400_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-400_narrow_13TeV-madgraph-v2.root"
#if( whichSample == "Grav_450" ):
#  Find_str.append("find /fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-450_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-450_narrow_13TeV-madgraph-v2/170420_233408 | grep root | grep -v failed > HADD/Grav_450_GluGluToBulkGravitonToHHTo2B2VTo2L2NuM-450.txt")
#  this_cat      = "cat HADD/Grav_450_* > HADD/Grav_450.txt"
#  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-450_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-450_narrow_13TeV-madgraph-v2.root @HADD/Grav_450.txt"
#  this_NtotPath = "/fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-450_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-450_narrow_13TeV-madgraph-v2.root"
#if( whichSample == "Grav_500" ):
#  Find_str.append("find /fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-500_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-500_narrow_13TeV-madgraph-v2/170420_233422 | grep root | grep -v failed > HADD/Grav_500_GluGluToBulkGravitonToHHTo2B2VTo2L2NuM-500.txt")
#  this_cat      = "cat HADD/Grav_500_* > HADD/Grav_500.txt"
#  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-500_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-500_narrow_13TeV-madgraph-v2.root @HADD/Grav_500.txt"
#  this_NtotPath = "/fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-500_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-500_narrow_13TeV-madgraph-v2.root"
#if( whichSample == "Grav_550" ):
#  Find_str.append("find /fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-550_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-550_narrow_13TeV-madgraph-v2/170420_233436 | grep root | grep -v failed > HADD/Grav_550_GluGluToBulkGravitonToHHTo2B2VTo2L2NuM-550.txt")
#  this_cat      = "cat HADD/Grav_550_* > HADD/Grav_550.txt"
#  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-550_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-550_narrow_13TeV-madgraph-v2.root @HADD/Grav_550.txt"
#  this_NtotPath = "/fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-550_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-550_narrow_13TeV-madgraph-v2.root"
#if( whichSample == "Grav_600" ):
#  Find_str.append("find /fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-600_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-600_narrow_13TeV-madgraph-v2/170420_233449 | grep root | grep -v failed > HADD/Grav_600_GluGluToBulkGravitonToHHTo2B2VTo2L2NuM-600.txt")
#  this_cat      = "cat HADD/Grav_600_* > HADD/Grav_600.txt"
#  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-600_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-600_narrow_13TeV-madgraph-v2.root @HADD/Grav_600.txt"
#  this_NtotPath = "/fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-600_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-600_narrow_13TeV-madgraph-v2.root"
#if( whichSample == "Grav_650" ):
#  Find_str.append("find /fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-650_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-650_narrow_13TeV-madgraph-v2/170420_233504 | grep root | grep -v failed > HADD/Grav_650_GluGluToBulkGravitonToHHTo2B2VTo2L2NuM-650.txt")
#  this_cat      = "cat HADD/Grav_650_* > HADD/Grav_650.txt"
#  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-650_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-650_narrow_13TeV-madgraph-v2.root @HADD/Grav_650.txt"
#  this_NtotPath = "/fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-650_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-650_narrow_13TeV-madgraph-v2.root"
#if( whichSample == "Grav_700" ):
#  Find_str.append("find /fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-700_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-700_narrow_13TeV-madgraph-v2/170420_233519 | grep root | grep -v failed > HADD/Grav_700_GluGluToBulkGravitonToHHTo2B2VTo2L2NuM-700.txt")
#  this_cat      = "cat HADD/Grav_700_* > HADD/Grav_700.txt"
#  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-700_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-700_narrow_13TeV-madgraph-v2.root @HADD/Grav_700.txt"
#  this_NtotPath = "/fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-700_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-700_narrow_13TeV-madgraph-v2.root"
#if( whichSample == "Grav_800" ):
#  Find_str.append("find /fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-800_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-800_narrow_13TeV-madgraph-v2/170420_233532 | grep root | grep -v failed > HADD/Grav_800_GluGluToBulkGravitonToHHTo2B2VTo2L2NuM-800.txt")
#  this_cat      = "cat HADD/Grav_800_* > HADD/Grav_800.txt"
#  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-800_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-800_narrow_13TeV-madgraph-v2.root @HADD/Grav_800.txt"
#  this_NtotPath = "/fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-800_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-800_narrow_13TeV-madgraph-v2.root"
#if( whichSample == "Grav_900" ):
#  Find_str.append("find /fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-900_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-900_narrow_13TeV-madgraph-v2/170420_233546 | grep root | grep -v failed > HADD/Grav_900_GluGluToBulkGravitonToHHTo2B2VTo2L2NuM-900.txt")
#  this_cat      = "cat HADD/Grav_900_* > HADD/Grav_900.txt"
#  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-900_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-900_narrow_13TeV-madgraph-v2.root @HADD/Grav_900.txt"
#  this_NtotPath = "/fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-900_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-900_narrow_13TeV-madgraph-v2.root"
#if( whichSample == "Grav_1000" ):
#  Find_str.append("find /fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-1000_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-1000_narrow_13TeV-madgraph-v2/170420_233601 | grep root | grep -v failed > HADD/Grav_1000_GluGluToBulkGravitonToHHTo2B2VTo2L2NuM-1000.txt")
#  this_cat      = "cat HADD/Grav_1000_* > HADD/Grav_1000.txt"
#  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-1000_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-1000_narrow_13TeV-madgraph-v2.root @HADD/Grav_1000.txt"
#  this_NtotPath = "/fdata/hepx/store/user/lpernie/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-1000_narrow_13TeV-madgraph-v2/crab_GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-1000_narrow_13TeV-madgraph-v2.root"
if( whichSample == "TT" ):
  Find_str.append("find /fdata/hepx/store/user/lpernie/TTTo2L2Nu_13TeV-powheg/crab_TTTo2L2Nu_13TeV-powheg/170421_015437 | grep root | grep -v failed > HADD/TT_TTTo2L2Nu13TeV-powheg.txt")
  this_cat      = "cat HADD/TT_* > HADD/TT.txt"
  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/TTTo2L2Nu_13TeV-powheg/crab_TTTo2L2Nu_13TeV-powheg.root @HADD/TT.txt"
  this_NtotPath = "/fdata/hepx/store/user/lpernie/TTTo2L2Nu_13TeV-powheg/crab_TTTo2L2Nu_13TeV-powheg.root"
if( whichSample == "DY" ):
  Find_str.append("find /fdata/hepx/store/user/lpernie/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/170421_015450 | grep root | grep -v failed > HADD/DY_DYJetsToLLM-10to50.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/DYToLL_0J_13TeV-amcatnloFXFX-pythia8/crab_DYToLL_0J_13TeV-amcatnloFXFX-pythia8/170421_015503 | grep root | grep -v failed > HADD/DY_DYToLL0J.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/DYToLL_1J_13TeV-amcatnloFXFX-pythia8/crab_DYToLL_1J_13TeV-amcatnloFXFX-pythia8/170421_015517 | grep root | grep -v failed > HADD/DY_DYToLL1J.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/DYToLL_2J_13TeV-amcatnloFXFX-pythia8/crab_DYToLL_2J_13TeV-amcatnloFXFX-pythia8/170421_015531 | grep root | grep -v failed > HADD/DY_DYToLL2J.txt")
  this_cat      = "cat HADD/DY_* > HADD/DY.txt"
  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/DYToLL_2J_13TeV-amcatnloFXFX-pythia8/crab_DYToLL_2J_13TeV-amcatnloFXFX-pythia8.root @HADD/DY.txt"
  this_NtotPath = "/fdata/hepx/store/user/lpernie/DYToLL_2J_13TeV-amcatnloFXFX-pythia8/crab_DYToLL_2J_13TeV-amcatnloFXFX-pythia8.root"
if( whichSample == "VV" ):
  Find_str.append("find /fdata/hepx/store/user/lpernie/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/crab_ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/170421_015545 | grep root | grep -v failed > HADD/VV_ZZTo2L2Q13TeV.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/ZZTo2L2Nu_13TeV_powheg_pythia8/crab_ZZTo2L2Nu_13TeV_powheg_pythia8/170421_015558 | grep root | grep -v failed > HADD/VV_ZZTo2L2Nu13TeV.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/ZZTo4L_13TeV_powheg_pythia8/crab_ZZTo4L_13TeV_powheg_pythia8/170421_015611 | grep root | grep -v failed > HADD/VV_ZZTo4L13TeV.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/WWToLNuQQ_aTGC_13TeV-madgraph-pythia8/crab_WWToLNuQQ_aTGC_13TeV-madgraph-pythia8/170421_015624 | grep root | grep -v failed > HADD/VV_WWToLNuQQaTGC.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/WWTo2L2Nu_MWW-600To800_aTGC_13TeV-amcatnloFXFX-madspin-pythia8/crab_WWTo2L2Nu_MWW-600To800_aTGC_13TeV-amcatnloFXFX-madspin-pythia8/170421_015637 | grep root | grep -v failed > HADD/VV_WWTo2L2NuMWW-600To800.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/crab_WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/170421_015651 | grep root | grep -v failed > HADD/VV_WZTo2L2Q13TeV.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8/crab_WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8/170421_015705 | grep root | grep -v failed > HADD/VV_WZTo1L3Nu13TeV.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/crab_WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/170421_015719 | grep root | grep -v failed > HADD/VV_WZTo1L1Nu2Q13TeV.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/crab_WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/170421_015733 | grep root | grep -v failed > HADD/VV_WZTo3LNuTuneCUETP8M1.txt")
  this_cat      = "cat HADD/VV_* > HADD/VV.txt"
  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/crab_WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8.root @HADD/VV.txt"
  this_NtotPath = "/fdata/hepx/store/user/lpernie/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/crab_WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8.root"
if( whichSample == "sT" ):
  Find_str.append("find /fdata/hepx/store/user/lpernie/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/crab_ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/170421_015747 | grep root | grep -v failed > HADD/sT_STt-channel.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/crab_ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/170421_015801 | grep root | grep -v failed > HADD/sT_STt-channel.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/crab_ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/170421_015814 | grep root | grep -v failed > HADD/sT_STs-channel.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/crab_ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/170421_015828 | grep root | grep -v failed > HADD/sT_STtW.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/crab_ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/170421_015842 | grep root | grep -v failed > HADD/sT_STtW.txt")
  this_cat      = "cat HADD/sT_* > HADD/sT.txt"
  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/crab_ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1.root @HADD/sT.txt"
  this_NtotPath = "/fdata/hepx/store/user/lpernie/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/crab_ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1.root"
if( whichSample == "Wjet" ):
  Find_str.append("find /fdata/hepx/store/user/lpernie/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/170421_015855 | grep root | grep -v failed > HADD/Wjet_WJetsToLNuTuneCUETP8M1.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/170421_015910 | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-100To200.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/170421_015924 | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-200To400.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/170421_015937 | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-400To600.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/170421_015951 | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-600To800.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/170421_020004 | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-800To1200.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/170421_020018 | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-1200To2500.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/170421_020035 | grep root | grep -v failed > HADD/Wjet_WJetsToLNuHT-2500ToInf.txt")
  this_cat      = "cat HADD/Wjet_* > HADD/Wjet.txt"
  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root @HADD/Wjet.txt"
  this_NtotPath = "/fdata/hepx/store/user/lpernie/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"
if( whichSample == "ttV" ):
  Find_str.append("find /fdata/hepx/store/user/lpernie/TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/170421_020049 | grep root | grep -v failed > HADD/ttV_TTWJetsToQQTuneCUETP8M1.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/170421_020103 | grep root | grep -v failed > HADD/ttV_TTWJetsToLNuTuneCUETP8M1.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/170421_020117 | grep root | grep -v failed > HADD/ttV_TTZToQQTuneCUETP8M1.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/170421_020133 | grep root | grep -v failed > HADD/ttV_TTZToLLNuNuM-10.txt")
  this_cat      = "cat HADD/ttV_* > HADD/ttV.txt"
  this_hadd     = "hadd -T -f -k /fdata/hepx/store/user/lpernie/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root @HADD/ttV.txt"
  this_NtotPath = "/fdata/hepx/store/user/lpernie/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root"
# Data
if( whichSample == "Data" ):
  Find_str.append("find /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Hhh_Run2016B-23Sep2016-v3/170421_020242 | grep root | grep -v failed > HADD/DATA_Hhh_Run2016B-23Sep2016-v3.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Hhh_Run2016C-23Sep2016-v1/170421_020255 | grep root | grep -v failed > HADD/DATA_Hhh_Run2016C-23Sep2016-v1.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Hhh_Run2016D-23Sep2016-v1/170421_020311 | grep root | grep -v failed > HADD/DATA_Hhh_Run2016D-23Sep2016-v1.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Hhh_Run2016E-23Sep2016-v1/170421_020325 | grep root | grep -v failed > HADD/DATA_Hhh_Run2016E-23Sep2016-v1.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Hhh_Run2016F-23Sep2016-v1/170421_020338 | grep root | grep -v failed > HADD/DATA_Hhh_Run2016F-23Sep2016-v1.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Hhh_Run2016G-23Sep2016-v1/170421_020352 | grep root | grep -v failed > HADD/DATA_Hhh_Run2016G-23Sep2016-v1.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Hhh_Run2016H-03Feb2017_ver2-v1/170421_020405 | grep root | grep -v failed > HADD/DATA_Hhh_Run2016H-03Feb2017_ver2-v1.txt")
  Find_str.append("find /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Hhh_Run2016H-03Feb2017_ver3-v1/170421_020420 | grep root | grep -v failed > HADD/DATA_Hhh_Run2016H-03Feb2017_ver3-v1.txt")
  this_cat      = "cat HADD/DATA_* > HADD/DATA.txt"

# Running
TCha = ROOT.TChain(tree_name)
for this_find in Find_str:
  os.system(this_find)
os.system(this_cat)
if (makeHadd=="HaddYes" and whichSample!="Data"): os.system(this_hadd)
if( whichSample != "Data" ):
  Ntot_path = this_NtotPath
  MyFile =  ROOT.TFile.Open(Ntot_path,"read");
  h_prehlt  =  ROOT.TH1F(MyFile.Get("TriggerResults/hevent_filter")); nTOT_prehlt = h_prehlt.GetBinContent(2)
  h_posthlt =  ROOT.TH1F(MyFile.Get("DiHiggsWWBBAna/hevent")); nTOT_posthlt = h_posthlt.GetBinContent(2);
else:
  nTOT_prehlt = 0.; nTOT_posthlt = 0.
with open(this_cat.split(">")[1].split(" ")[1],"r") as f:
  for line in f:
    if not line.isspace():
      TCha.Add(str(line[:-1]))
print whichSample, "TChain has", TCha.GetEntries(), "entries."
f = ROOT.TFile('/fdata/hepx/store/user/%s/Hhh_For_Plotting/'%user + whichSample + '.root','recreate'); f.cd()

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

nEv = 0
for ev in TCha:
  if (doTest and nEv%1000 == 0 ):
      print "ev ",nEv
  elif (nEv%10000 == 0):
      print "ev ",nEv
  
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
f.Close()
print "Done."
