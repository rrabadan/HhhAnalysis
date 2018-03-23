#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor

from HHbbWWProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jecUncertainties import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetUncertainties import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.lepSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.mht import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import *
#from  PhysicsTools.NanoAODTools.postprocessing.examples.puWeightProducer import *

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
else:
    print "ERROR: This is NOT a valid directory: ",inputdir
    exit()
if not os.path.isdir(outputdir):
    os.system("mkdir "+outputdir)

modules = []
#selection='''(Sum$(Electron_pt > 20 && Electron_mvaSpring16GP_WP90) >= 2  ||
# Sum$(Muon_pt > 20) >= 2 ||
# Sum$(Electron_pt > 20 && Electron_mvaSpring16GP_WP80) >= 1   ||
# Sum$(Muon_pt > 20 && Muon_tightId) >= 1 ||
# (Sum$(Muon_pt > 20) == 0 && Sum$(Electron_pt > 20 && Electron_mvaSpring16GP_WP90) == 0 && MET_pt > 150 ) ) 
# &&  Sum$((abs(Jet_eta)<2.5 && Jet_pt > 20 && Jet_jetId)) >= 2 && Entry$ < 1000 
###mhtProducer(jetSelection, muonSelection, electronSelection), adding the selected obj and filling MHT_pt and MHT_phi
mht_hh = lambda : mhtProducer( lambda j : j.pt > 20 and abs(j.eta) < 2.4,
                            lambda mu : mu.pt > 10 and abs(mu.eta) < 2.4,
                            lambda el : el.pt > 10 and abs(el.eta) < 2.5 )

print "=============================================================="
if "DoubleMuon" in jobtype:
    print "Dataset DoubleMuon "
    modules = [mht_hh(), HHbbWWProducer(False, triggertype = "DoubleMuon", verbose=1)]
elif "DoubleEG" in jobtype:
    print "Dataset DoubleEG "
    modules = [mht_hh(), HHbbWWProducer(False, triggertype = "DoubleEG", verbose=1)]
elif "MuonEG" in jobtype:
    print "Dataset MuonEG "
    modules = [mht_hh(), HHbbWWProducer(False, triggertype = "MuonEG", verbose=1)]
elif jobtype != "":
    print "MC samples "
    jsonfile = None
    modules = [puAutoWeight(), btagSF2016(), mht_hh(), HHbbWWProducer(True, verbose = 1) ]
else:
    print "jobtype to run is not found, exit "
    exit()
print "=============================================================="
print "inputfiles ", inputfiles
print "=============================================================="
print "outputdir ", outputdir
print "=============================================================="
#p=PostProcessor(".",files,selection.replace('\n',' '),"keep_and_drop.txt",[puAutoWeight(),jetmetUncertainties2016All(), btagSF2016(), hhbbWW()],provenance=True)
p=PostProcessor(outputdir, inputfiles,"1","keep_and_drop.txt", modules, jsonInput = jsonfile, provenance=True)
#p=PostProcessor(".",filesdata_MuEl,"1","keep_and_drop.txt",[mht_hh(), hhbbWW_data("MuonEG")],provenance=True)
#p=PostProcessor(".",filesdata_ElEl,"1","keep_and_drop.txt",[mht_hh(), hhbbWW_data("DoubleEG")],provenance=True)

p.run()


