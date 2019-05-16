#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True



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

sys.path.append("/home/taohuang/DiHiggsAnalysis/CMSSW_9_4_0_pre1/src/HhhAnalysis/python/NanoAOD/")
from countHistogramProducer import *
from genParticleProducer import *

files = ["root://cms-xrd-global.cern.ch//store/mc/RunIISummer16NanoAOD/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-400_narrow_13TeV-madgraph-v2/NANOAODSIM/PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/00000/82738046-5711-E811-ADAF-0CC47A4D75F0.root"]
filesTTbar= [
'root://cms-xrd-global.cern.ch//store/group/cmst3/group/nanoAOD/NanoTestProd006/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer17MiniAOD-92X-NanoCrabProd006/171006_155430/0000/nanolzma_1.root',
]
filesTTbar = ["/fdata/hepx/store/user/taohuang/HH_NanoAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIIFall17_NanoAOD.root"]
filesSignal = ["/fdata/hepx/store/user/taohuang/HH_NanoAOD/GravitonToHHTo2B2VTo2L2Nu_M-400_narrow_13TeV_NanoAOD_2.root"]
filesdata_MuMu = ["/fdata/hepx/store/user/taohuang/HH_NanoAOD/Run2016B_DoubleMuon_EC6E35C0-630C-E811-AB01-00215A450982.root"]
filesdata_ElEl = ["/fdata/hepx/store/user/taohuang/HH_NanoAOD/Run2016B_DoubleEG_E8EE52A1-4D0C-E811-BFE1-90E2BACC5EEC.root"]
filesdata_MuEl = ["/fdata/hepx/store/user/taohuang/HH_NanoAOD/Run2016B_MuonEG_E480DF59-890C-E811-995A-0242AC130002.root"]
filesDY0j = ["/fdata/hepx/store/user/taohuang/NANOAOD/DYToLL_0J_13TeV-amcatnloFXFX-pythia8/84C5AEED-8E19-E811-9948-F04DA2757237.root"]
filesDY2017 = ["/fdata/hepx/store/mc/RunIIFall17NanoAOD/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/46ADA924-A2AC-E811-B855-003048C45C8A.root"]


#selection='''(Sum$(Electron_pt > 20 && Electron_mvaSpring16GP_WP90) >= 2  ||
mht_hh = lambda : mhtProducer( lambda j : j.pt > 20 and abs(j.eta) < 2.4,
                            lambda mu : mu.pt > 10 and abs(mu.eta) < 2.4,
                            lambda el : el.pt > 10 and abs(el.eta) < 2.5 )

##cp leptonSF/cMVAv2_Moriond17_B_H.csv   ../../../PhysicsTools/NanoAODTools/data/btagSF/
outputdir = "/fdata/hepx/store/user/taohuang/HH_NanoAOD"
#p=PostProcessor(outputdir, filesTTbar, cut = "1", branchsel = "keep_and_drop_pre.txt", modules = [ ], friend = False, provenance=True, outputbranchsel="keep_and_drop_out.txt")
modules = [puWeightyear(Runyear), countHistogramAll_2016(), jetmetUncertaintiesyear(Runyear), btagSFyear(Runyear),  muonScaleRes(Runyear), mht_hh(), HHbbWWProducer(True, DYestimation = True, verbose = 1)]
p=PostProcessor(outputdir, filesDY2017,"1","keep_and_drop.txt", modules, friend = True, provenance=True)


p.run()


