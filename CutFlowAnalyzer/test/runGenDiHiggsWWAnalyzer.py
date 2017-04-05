import FWCore.ParameterSet.Config as cms
import os
import sys

process = cms.Process("DiHiggsAnalyzer")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring(
        #'file:/eos/uscms/store/user/tahuang/DiHiggs/out_sim.root'
	#'/store/mc/RunIISpring16MiniAODv1/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext3-v1/00000/02F64C80-990E-E611-A2FE-842B2B185476.root'
	'/store/mc/RunIISpring16MiniAODv2/DYBBJetsToLL_M-10To70_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/premix_withHLT_80X_mcRun2_asymptotic_v14-v3/100000/0072F5C1-5089-E611-8340-D8D385AF8AE4.root'
	#'/store/data/Run2016H/DoubleMuon/MINIAOD/03Feb2017_ver2-v1/110000/08ADA6AA-D3EC-E611-AF17-B083FED42488.root'
	#'file:/eos/uscms/store/user/tahuang/DiHiggs/out_miniaod.root'
    )
)

from HhhAnalysis.MCProduction.InputFileHelpers import *
#inputdir = ['/fdata/hepx/store/user/tahuang/xSM_HeavyHiggs2DiHiggs2bbWW_B3_leptonW_CMSSW80X_13TeV_1M/xSM_HeavyHiggs2DiHiggs2bbWW_B3_leptonW_CMSSW80X_13TeV_1M/170329_023747/0000/']
#inputdir = ['/fdata/hepx/store/user/tahuang/xSM_HeavyHiggs2DiHiggs2bbWW_B3_leptonW_CMSSW80X_13TeV_10k/xSM_HeavyHiggs2DiHiggs2bbWW_B3_leptonW_CMSSW80X_13TeV_10k/170330_023219/0000/']
#inputdir = ['/fdata/hepx/store/user/tahuang/xSM_HeavyHiggs2DiHiggs2bbWW2bbll_B3_CMSSW80X_13TeV_10k_nobbwwfilter/xSM_HeavyHiggs2DiHiggs2bbWW2bbll_B3_CMSSW80X_13TeV_10k_nobbwwfilter/170331_052059/0000/']
inputdir = ['/fdata/hepx/store/user/tahuang/xSM_HeavyHiggs2DiHiggs2bbWW2bbll_B3_CMSSW80X_13TeV_10k_nobbwwfilter_updatehwidth/xSM_HeavyHiggs2DiHiggs2bbWW2bbll_B3_CMSSW80X_13TeV_10k_nobbwwfilter_updatehwidth/170331_053038/0000/']
process = useInputDir(process, inputdir)

process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(10000) 
)

process.MessageLogger = cms.Service("MessageLogger", 
    destinations = cms.untracked.vstring("cout"), 
    cout = cms.untracked.PSet(threshold = cms.untracked.string("ERROR"))
)

process.eventCounterFilter = cms.EDFilter("EventCounterFilter")
import HLTrigger.HLTfilters.triggerResultsFilter_cfi as hlt
#print hlt," hlt.triggerResultsFilter ",hlt.triggerResultsFilter
# accept if any path succeeds (explicit OR)
"""
process.hltfilter = hlt.triggerResultsFilter.clone(
	HLTPaths = ('HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*','HLT_Mu17*'),
	l1tResults = '',#not use L1t results
	throw = cms.bool(False) 
"""
process.hltfilter = cms.EDFilter( "TriggerResultsFilter",
	triggerConditions = cms.vstring( 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*',
	    				 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*',
					 'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*',
					 'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*',
	    ),
	hltResults = cms.InputTag( "TriggerResults","","HLT"),
	#l1tResults = cms.InputTag( "hltGtDigis" ),
	l1tResults = cms.InputTag( "" ),
	l1tIgnoreMask = cms.bool( False ),
	l1techIgnorePrescales = cms.bool( False ),
	daqPartitions = cms.uint32( 1 ),
	throw = cms.bool(True)    
	)
muonPOGSFdir = os.getenv( "CMSSW_BASE" ) +"/src/HhhAnalysis/CutFlowAnalyzer/test/MuonEffAndSF_2016Data/"
process.DiHiggsWWBBAna = cms.EDAnalyzer('DiHiggsWWBBAnalyzer',
    verbose = cms.untracked.int32(0),
    SampleType = cms.untracked.int32(3), #enum {Data = 0, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, TTbar};//add other background
    sampleName = cms.untracked.int32(3), #B1 = 1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, tt, DYJets, DY0Jets, DY1Jets, DY2Jets
    onlyGenLevel = cms.bool(True),
    genParticles = cms.InputTag("genParticles"),
    #genParticles = cms.InputTag("prunedGenParticles"),#minAOD
    #genjets = cms.InputTag("slimmedGenJets"),
    #genjets = cms.InputTag("ak4GenJetsNoNu"),
    genjets = cms.InputTag("ak4GenJets"),
    #muons = cms.InputTag("cleanPatPFMuonsTriggerMatch"),
    muons = cms.InputTag("slimmedMuons"),
    triggerSFFile = cms.string(muonPOGSFdir+"EfficienciesAndSF_BCDEF_trigger.root"),
    isoSFFile = cms.string(muonPOGSFdir+"EfficienciesAndSF_BCDEF_ISO.root"),
    idSFFile = cms.string(muonPOGSFdir+"EfficienciesAndSF_BCDEF_ID.root"),
    trackingSFFile = cms.string(muonPOGSFdir+"EfficienciesAndSF_BCDEFGH_Tracking.root"),
    triggerSFhist = cms.string("IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio"),
    isoSFhist = cms.string("TightISO_MediumID_pt_eta/abseta_pt_ratio"),
    idSFhist = cms.string("MC_NUM_MediumID2016_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio"),
    trackingSFhist = cms.string("ratio_eff_eta3_dr030e030_corr"),
    electrons = cms.InputTag("slimmedElectrons"),
    jets = cms.InputTag("slimmedJets"),
    mets = cms.InputTag("slimmedMETs"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    triggerEvent = cms.InputTag("patTriggerEvent"),
    tracks = cms.InputTag("generalTracks"),
    TriggerResults = cms.InputTag("TriggerResults","","RECO"),
    TrackRefitter = cms.InputTag("TrackRefitter"),
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    Traj = cms.InputTag("TrackRefitter"),
    debug = cms.untracked.bool(False),
    #simulation = cms.bool(True),
    runMMC = cms.bool(False)
)
#print "process.DiHiggsWWBBAna ",process.DiHiggsWWBBAna
process.dump=cms.EDAnalyzer('EventContentAnalyzer')

process.TFileService = cms.Service("TFileService",
    #fileName = cms.string("out_ana.root")
    fileName = cms.string("file:/fdata/hepx/store/user/taohuang/DiHiggsAnalysisSample/out_ana_genlevel_B3_10k_nobbwwfilter_updatehiggswidth.root")
)

#process.phlt = cms.Path(process.hltfilter)
process.pDiHiggsWWBBAna = cms.Path(
    #process.hltfilter*
    process.eventCounterFilter *
    process.DiHiggsWWBBAna
)

process.pdump = cms.Path(process.dump)	

#process.schedule = cms.Schedule(process.phlt*process.pDiHiggsWWBBAna)
process.schedule = cms.Schedule( process.pDiHiggsWWBBAna)
print "----------------------------------------------\n"
print "Inputfiles: \n"
print process.source.fileNames
print "Outputfile: \n"
print process.TFileService.fileName
