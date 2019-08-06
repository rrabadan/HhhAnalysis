#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

import warnings

#from RunConfiguration import *
Runyear = 2017

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor

sys.path.append('/home/taohuang/DiHiggsAnalysis/CMSSW_9_4_0_pre1/src/HhhAnalysis/python/NanoAOD')
#from countHistogramProducer import *
from genParticleProducer import *
from HHbbWWProducer_sync import *


from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import *
btagSF2016_cMVA = lambda : btagSFProducer("2016",  algo = 'cmva', sfFileName='cMVAv2_Moriond17_B_H.csv')## file is under NanoAODTools
btagSF2017_deepCSV = lambda : btagSFProducer("2017", algo = 'deepcsv', sfFileName='DeepCSV_94XSF_V3_B_F.csv')
btagSF2018_deepCSV = lambda : btagSFProducer("2018", algo = 'deepcsv', sfFileName='DeepCSV_102XSF_V1.csv')
def btagSFyear(year):
    return{
	2016: btagSF2016_cMVA(),
	2017: btagSF2017_deepCSV(),
	2018: btagSF2018_deepCSV()
    }[year]

from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jecUncertainties import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetUncertainties import *
def jetmetUncertaintiesyear(year):
    """jet uncertainties and recalibration for MC"""
    return{
	2016: jetmetUncertainties2016(),
	2017: jetmetUncertainties2017(),
	2018: jetmetUncertainties2018()
    }[year]

from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetRecalib import *

from PhysicsTools.NanoAODTools.postprocessing.modules.common.lepSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.muonScaleResProducerTao import *
#from PhysicsTools.NanoAODTools.postprocessing.modules.common.muonScaleResProducer import *
muonScaleRes =  lambda year : muonScaleResProducer('roccor.Run2.v3', 'RoccoR%d.txt'%year) 
#def muonScaleRes(year):
#    return {
#	2016:muonScaleRes2016(),
#	2017:muonScaleRes2017(),
#	2018:muonScaleRes2018()
#    }[year]

from PhysicsTools.NanoAODTools.postprocessing.modules.jme.mht import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import *
#from  PhysicsTools.NanoAODTools.postprocessing.examples.puWeightProducer import *
def puWeightyear(year):
    return {
      2016: puAutoWeight_2016(),
      2017: puAutoWeight_2017(),
      2018: puAutoWeight_2018()
    }[year]

###mhtProducer(jetSelection, muonSelection, electronSelection), adding the selected obj and filling MHT_pt and MHT_phi
mht_hh = lambda : mhtProducer( lambda j : j.pt > 20 and abs(j.eta) < 2.4,
                            lambda mu : mu.pt > 10 and abs(mu.eta) < 2.4,
                            lambda el : el.pt > 10 and abs(el.eta) < 2.5 )

#outputdir = "/fdata/hepx/store/user/taohuang/HH_NanoAOD/"
outputdir = "/eos/uscms/store/user/tahuang/HHbbWW_sync/"
filesTTbar = ["/fdata/hepx/store/user/taohuang/HH_NanoAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIIFall17_NanoAOD.root"]
filesSignal = ["/fdata/hepx/store/user/taohuang/HH_NanoAOD/GravitonToHHTo2B2VTo2L2Nu_M-400_narrow_13TeV_NanoAOD_2.root"]
filesTTbar2017 = ["/fdata/hepx/store/mc/RunIIFall17NanoAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/100000/2E9913D6-BAA9-E811-8ABE-0CC47A4DEF3E.root"]
fileSignal2017 = ["/fdata/hepx/store/mc/RunIIFall17NanoAOD/GluGluToRadionToHHTo2B2VTo2L2Nu_M-350_narrow_13TeV-madgraph_correctedcfg/NANOAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/70B64607-EEB2-E811-B1D7-A0369FD0B268.root"]
filesDY1J2017 = ["/fdata/hepx/store/mc/RunIIFall17NanoAODv4/DYJetsToLL_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/90000/62B5DB6B-F42C-D648-A6EF-5F26495985BF.root"]
#filesRadion750 = ["/fdata/hepx/store/user/taohuang/Radion750_NanoAOD_RunIIFall17MiniAODv2-PU2017/myNanoProdMc_NANOAOD_11_sync.root"]
#filesRadion750 = ["/fdata/hepx/store/user/taohuang/Radion750_NanoAOD_RunIIFall17MiniAODv2-PU2017/myNanoProdMc_NANOAOD_1.root"]
#fileRadion750 = ["/eos/uscms/store/user/tahuang/GluGluToRadionToHHTo2B2Tau_M-750_narrow_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018/190728_042406/0000/myNanoProdMcv2_NANO_10.root"]
#filesRadion750 = ["/eos/uscms/store/user/tahuang/GluGluToRadionToHHTo2B2VTo2L2Nu_M-750_narrow_13TeV-madgraph-v2/RunIISummer16MiniAODv3-PU2017_NanoAOD/190728_045153/0000/myNanoProdMc2016_NANO_2.root"]
#filesRadion750 = ["/eos/uscms/store/user/tahuang/GluGluToRadionToHHTo2B2VTo2L2Nu_M-750_narrow_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_NanoAOD_10010X_M-750_CMSSW10215_v2/190627_205303/0000/myNanoProdMcv2_NANO_11.root"]
filesRadion750 = ["/eos/uscms/store/user/tahuang/GluGluToRadionToHHTo2B2VTo2L2Nu_M-750_narrow_13TeV-madgraph_correctedcfg/RunIIFall17MiniAODv2-PU2017_NanoAOD_20190730/190731_060502/0000/myNanoProdMc2017_NANO_11.root"]


#modules = [ puWeightyear(Runyear), countHistogramAll_2016(), jetmetUncertaintiesyear(Runyear), btagSFyear(Runyear),  mht_hh(), HHbbWWProducer(True, verbose = 1) ]
modules = [ puWeightyear(Runyear), muonScaleRes(Runyear),  mht_hh(), HHbbWWProducer(True, verbose = 1) ]
#p=PostProcessor(outputdir, filesTTbar2017,"1","keep_and_drop.txt", modules, friend = True, jsonInput = None, provenance=True)
#p=PostProcessor(outputdir, filesDY1J2017,"1","keep_and_drop.txt", modules, friend = True, jsonInput = None, provenance=True)
p=PostProcessor(outputdir, filesRadion750,"1","keep_and_drop.txt", modules, friend = True, jsonInput = None, provenance=True, outtreeName="syncTree")
print "run Postprocessor here "
p.run()

