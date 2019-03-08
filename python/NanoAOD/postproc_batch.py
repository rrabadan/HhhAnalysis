#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
#sys.argv.append( '-b' )
#sys.argv.append( '-q' )

import warnings

from RunConfiguration import *

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor

sys.path.append('/home/taohuang/DiHiggsAnalysis/CMSSW_9_4_0_pre1/src/HhhAnalysis/python/NanoAOD')
from countHistogramProducer import *
from genParticleProducer import *
from HHbbWWProducer import *


from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import *
btagSF2016_cMVA = lambda : btagSFProducer("2016",  algo = 'cmva', sfFileName='cMVAv2_Moriond17_B_H.csv')## file is under NanoAODTools
btagSF2017_deepCSV = lambda : btagSFProducer("2017", algo = 'deepcsv', sfFileName='DeepCSV_94XSF_V3_B_F.csv')
def btagSFyear(year):
    return{
	2016: btagSF2016_cMVA(),
	2017: btagSF2017_deepCSV()
    }[year]

from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jecUncertainties import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetUncertainties import *
def jetmetUncertaintiesyear(year):
    """jet uncertainties and recalibration for MC"""
    return{
	2016: jetmetUncertainties2016(),
	2017: jetmetUncertainties2017()
    }[year]

from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetRecalib import *

from PhysicsTools.NanoAODTools.postprocessing.modules.common.lepSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.muonScaleResProducerTao import *
muonScaleRes =  lambda year : muonScaleResProducer('roccor.Run2.v3', 'RoccoR%d.txt'%year) 

from PhysicsTools.NanoAODTools.postprocessing.modules.jme.mht import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import *
#from  PhysicsTools.NanoAODTools.postprocessing.examples.puWeightProducer import *
def puWeightyear(year):
    return {
      2016: puWeight(),
      2017: puAutoWeight()
    }[year]

import argparse

parser = argparse.ArgumentParser(description='postprocess_batch')
parser.add_argument("-i", "--inputdir", dest="inputdir", default="none", help="inputdir")
parser.add_argument("-json", "--json", dest="json", default=None, help="json file ")
parser.add_argument("-o", "--outputdir",dest="output", default="/fdata/hepx/store/user/taohuang/HH_NanoAOD/",help="output root file directory, [Defualt:fdata] ")
parser.add_argument("-j", "--jobtype", dest="jobtype", default="", help="jobtype ")
args = parser.parse_args()
jsonfile = args.json
#jsonfile = "Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt"
jobtype = args.jobtype
#outputdir = "/fdata/hepx/store/user/taohuang/HHNtuples_20180322/"
outputdir  = args.output
inputdir = args.inputdir
inputfiles = []
if os.path.isdir(inputdir):
    #print "This is a valid directory: ",inputdir
    ls = os.listdir(inputdir)
    inputfiles = [os.path.join(inputdir,x) for x in ls if x.endswith("root")]
elif os.path.isfile(inputdir):
    inputfiles = [inputdir]
elif inputdir.startswith("root://cms-xrd-global.cern.ch") or inputdir.startswith("root://cmsxrootd.fnal.gov/"):
    inputfiles = [inputdir]
else:
    print "ERROR: This is NOT a valid directory: ",inputdir
    exit()
if not os.path.isdir(outputdir):
    os.system("mkdir "+outputdir)

modules = []
mht_hh = lambda : mhtProducer( lambda j : j.pt > 20 and abs(j.eta) < 2.4,
                            lambda mu : mu.pt > 10 and abs(mu.eta) < 2.4,
                            lambda el : el.pt > 10 and abs(el.eta) < 2.5 )

print "=============================================================="
if "DoubleMuon" in jobtype:
    Runperiod = jobtype.replace("DoubleMuonRun","")[:5]
    print "Dataset DoubleMuon ",Runperiod
    modules = [mht_hh(),  muonScaleRes(Runyear), jetRecalibData(Runperiod), HHbbWWProducer(False, triggertype = "DoubleMuon", verbose=1)]
elif "DoubleEG" in jobtype:
    Runperiod = jobtype.replace("DoubleEGRun","")[:5]
    print "Dataset DoubleEG ",Runperiod
    modules = [mht_hh(),  muonScaleRes(Runyear), jetRecalibData(Runperiod), HHbbWWProducer(False, triggertype = "DoubleEG", verbose=1)]
elif "MuonEG" in jobtype:
    Runperiod = jobtype.replace("MuonEGRun","")[:5]
    print "Dataset MuonEG ",Runperiod
    modules = [mht_hh(),  muonScaleRes(Runyear), jetRecalibData(Runperiod), HHbbWWProducer(False, triggertype = "MuonEG", verbose=1)]
elif jobtype != "":
    print "MC samples "
    jsonfile = None
    #modules = [genHHAndTTbar(), puWeight(), countHistogramAll_2016(), jetmetUncertainties2016(), btagSF2016_cMVA(),  muonScaleRes2016(), mht_hh(), HHbbWWProducer(True, verbose = 1) ]
    modules = [ puWeight(), countHistogramAll_2016(), jetmetUncertaintiesyear(Runyear), btagSFyear(Runyear),  muonScaleRes(Runyear), mht_hh(), HHbbWWProducer(True, verbose = 1) ]

else:
    print "jobtype to run is not found, exit "
    exit()
print "=============================================================="
print "inputfiles ", inputfiles
print "=============================================================="
print "outputdir ", outputdir
print "=============================================================="
p=PostProcessor(outputdir, inputfiles,"1","keep_and_drop.txt", modules, friend = True, jsonInput = jsonfile, provenance=True)

p.run()


