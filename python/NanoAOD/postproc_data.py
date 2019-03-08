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



files = ["root://cms-xrd-global.cern.ch//store/mc/RunIISummer16NanoAOD/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-400_narrow_13TeV-madgraph-v2/NANOAODSIM/PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/00000/82738046-5711-E811-ADAF-0CC47A4D75F0.root"]
filesTTbar= [
'root://cms-xrd-global.cern.ch//store/group/cmst3/group/nanoAOD/NanoTestProd006/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer17MiniAOD-92X-NanoCrabProd006/171006_155430/0000/nanolzma_1.root',
]
filesTTbar = ["/fdata/hepx/store/user/taohuang/HH_NanoAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIIFall17_NanoAOD.root"]
filesSignal = ["/fdata/hepx/store/user/taohuang/HH_NanoAOD/GravitonToHHTo2B2VTo2L2Nu_M-400_narrow_13TeV_NanoAOD_2.root"]
#filesdata_MuMu = ["/fdata/hepx/store/user/taohuang/HH_NanoAOD/Run2016B_DoubleMuon_EC6E35C0-630C-E811-AB01-00215A450982.root"]
filesdata_ElEl = ["/fdata/hepx/store/user/taohuang/HH_NanoAOD/Run2016B_DoubleEG_E8EE52A1-4D0C-E811-BFE1-90E2BACC5EEC.root"]
filesdata_MuEl = ["/fdata/hepx/store/user/taohuang/HH_NanoAOD/Run2016B_MuonEG_E480DF59-890C-E811-995A-0242AC130002.root"]
filesdata_MuMu = ["/fdata/hepx/store/user/taohuang/NANOAOD/DoubleMuonRun2016Bver2/DEF277D4-550C-E811-8262-001E675A6653.root"]


filesdata_ElEl2017 = ["/fdata/hepx/store/data/Run2017B/DoubleEG/NANOAOD/31Mar2018-v1/70000/8254B705-BB44-E811-A808-FA163E95DCD9.root"]
filesdata_MuMu2017 = ["/fdata/hepx/store/data/Run2017B/DoubleMuon/NANOAOD/31Mar2018-v1/70000/A41543AB-C444-E811-9599-FA163E9592B6.root"]
filesdata_MuEl2017 = ["/fdata/hepx/store/data/Run2017B/MuonEG/NANOAOD/31Mar2018-v1/70000/B2D34F84-0D45-E811-AAA3-FA163E8A6597.root"]


#selection='''(Sum$(Electron_pt > 20 && Electron_mvaSpring16GP_WP90) >= 2  ||
# Sum$(Muon_pt > 20) >= 2 ||
# Sum$(Electron_pt > 20 && Electron_mvaSpring16GP_WP80) >= 1   ||
# Sum$(Muon_pt > 20 && Muon_tightId) >= 1 ||
# (Sum$(Muon_pt > 20) == 0 && Sum$(Electron_pt > 20 && Electron_mvaSpring16GP_WP90) == 0 && MET_pt > 150 ) ) 
# &&  Sum$((abs(Jet_eta)<2.5 && Jet_pt > 20 && Jet_jetId)) >= 2 && Entry$ < 1000 
#'''
#
#selectionALL='''Sum$(Electron_pt > 20 && Electron_mvaSpring16GP_WP90) >= 2  ||
# Sum$(Electron_pt > 20 && Electron_mvaSpring16GP_WP80) >= 1   ||
# Sum$(Jet_pt > 40 && Jet_jetId) >= 4   || 
#Sum$(Jet_pt *(abs(Jet_eta)<2.5 && Jet_pt > 20 && Jet_jetId)) > 160  || 
#MET_pt > 100  || Sum$(Muon_pt > 20 && Muon_tightId) >= 1
#'''
###mhtProducer(jetSelection, muonSelection, electronSelection), adding the selected obj and filling MHT_pt and MHT_phi
mht_hh = lambda : mhtProducer( lambda j : j.pt > 20 and abs(j.eta) < 2.4,
                            lambda mu : mu.pt > 10 and abs(mu.eta) < 2.4,
                            lambda el : el.pt > 10 and abs(el.eta) < 2.5 )

jsonfile = "Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt"
#jobtype = ""
outputdir = "/fdata/hepx/store/user/taohuang/HH_NanoAOD/"
#inputdir = 

#p=PostProcessor(".",files,selection.replace('\n',' '),"keep_and_drop.txt",[puAutoWeight(),jetmetUncertainties2016All(), btagSF2016, hhbbWW()],provenance=True)
#p=PostProcessor(".",files,selection.replace('\n',' '),"keep_and_drop.txt",[puAutoWeight(),jetmetUncertainties2016All(), btagSF2016(), hhbbWW()],provenance=True)
#p=PostProcessor(".",files,selection.replace('\n',' '),"keep_and_drop.txt",[puAutoWeight(), hhbbWW()],provenance=True)
#p=PostProcessor(".",filesSignal,"1","keep_and_drop.txt",[puAutoWeight(), lepSF(), btagSF2016(), mht_hh(), hhbbWW()],provenance=True)
#p=PostProcessor(outputdir, filesdata_MuMu,"1","keep_and_drop.txt",[muonScaleRes2016(),mht_hh(), HHbbWWProducer(False, triggertype = "DoubleMuon", verbose=4)], friend =True, jsonInput = jsonfile, provenance=True)
#p=PostProcessor(outputdir, filesdata_MuEl,"1","keep_and_drop.txt",[muonScaleRes2016(),mht_hh(), HHbbWWProducer(False, triggertype = "MuonEG", verbose=4)], friend = True, jsonInput = jsonfile, provenance=True)
#p=PostProcessor(outputdir,filesdata_ElEl2017,"1","keep_and_drop.txt",[ muonScaleRes(Runyear),  HHbbWWProducer(False, triggertype = "DoubleEG", verbose=4)], friend = True, provenance=True)
p=PostProcessor(outputdir,filesdata_MuMu2017,"1","keep_and_drop.txt",[ muonScaleRes(Runyear),  HHbbWWProducer(False, triggertype = "DoubleMuon", verbose=0)], friend = True, provenance=True)
#p=PostProcessor(outputdir,filesdata_MuEl2017,"1","keep_and_drop.txt",[ muonScaleRes(Runyear),  HHbbWWProducer(False, triggertype = "MuonEG", verbose=4)], friend = True, provenance=True)

p.run()


