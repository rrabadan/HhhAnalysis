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
#  'file:/fdata/hepx/store/user/lpernie/TEST_LOCALLY/DYJETS_7A385961-C6D9-E611-85B2-0025905B85BC.root'
  'file:/fdata/hepx/store/user/lpernie/TEST_LOCALLY/Run2016D-23Sep2016_MINIAOD_12B2DEA9-B68C-E611-99A4-0CC47A1DF810.root'
  )
)

from HhhAnalysis.MCProduction.InputFileHelpers import *
#inputdir = ['/fdata/hepx/store/user/tahuang/xSM_HeavyHiggs2DiHiggs2bbWW_B3_leptonW_CMSSW80X_13TeV_1M/xSM_HeavyHiggs2DiHiggs2bbWW_B3_leptonW_CMSSW80X_13TeV_1M/170329_023747/0000/']
#inputdir = ['/fdata/hepx/store/user/tahuang/xSM_HeavyHiggs2DiHiggs2bbWW_B3_leptonW_CMSSW80X_13TeV_10k/xSM_HeavyHiggs2DiHiggs2bbWW_B3_leptonW_CMSSW80X_13TeV_10k/170330_023219/0000/']
#process = useInputDir(process, inputdir)

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
triggerPaths = cms.vstring( 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*',
			     'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*',
			     'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*',
			     'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*',
	)
process.hltfilter = cms.EDFilter( "TriggerResultsFilter",
        triggerConditions = triggerPaths,
	hltResults = cms.InputTag( "TriggerResults","","HLT"),
	#l1tResults = cms.InputTag( "hltGtDigis" ),
	l1tResults = cms.InputTag( "" ),
	l1tIgnoreMask = cms.bool( False ),
	l1techIgnorePrescales = cms.bool( False ),
	daqPartitions = cms.uint32( 1 ),
	throw = cms.bool(True)    
)
print sys.argv
sample = 0#int(sys.argv[1])
print "sample",sample

process.DiHiggsWWBBAna = cms.EDAnalyzer('DiHiggsWWBBAnalyzer',
  verbose = cms.untracked.int32(0),
  #enum {Data = 0, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, TTbar, DYJets, DY0Jets, DY1Jets, DY2Jets, ZZTo2L2Q, ZZTo2L2Nu, ZZTo4L, WWToLNuQQ, WWTo2L2Nu, WZTo2L2Q, WZTo1L3Nu, WZTo1L1Nu2Q, WZTo3LNu, ST_tchannel_top, ST_tchannel_antitop, ST_schannel, ST_tW_antitop, ST_tW_top, WJetsToLNu, WJetsToLNu_HT100To200, WJetsToLNu_HT200To400, WJetsToLNu_HT400To600, WJetsToLNu_HT600To800, WJetsToLNu_HT800To1200, WJetsToLNu_HT1200To2500, WJetsToLNu_HT2500ToInf, TTWJetsToQQ, TTWJetsToLNu, TTZToQQ, TTZToLLNuNu };//add other background
  SampleType = cms.untracked.int32(sample),
  #genParticles = cms.InputTag("genParticles"),
  genParticles = cms.InputTag("prunedGenParticles"),#minAOD
  genjets = cms.InputTag("slimmedGenJets"),
  #genjets = cms.InputTag("ak4GenJetsNoNu"),
  #genjets = cms.InputTag("ak4GenJets"),

  #trigger matching
  doTriggerMatching = cms.bool(True),
  hltPaths = triggerPaths,
  deltaPtRel_trigger = cms.untracked.double(.5),
  deltaR_trigger  = cms.untracked.double(.1),

  #reco 
  #muons = cms.InputTag("cleanPatPFMuonsTriggerMatch"),
  muons = cms.InputTag("slimmedMuons"),
  #2016data: Run BCDEF use 2016Medium, GH use Medium
  mu_id = cms.untracked.string("2016Medium"),
  mu_PFIso = cms.untracked.double(0.15),#tight iso
  triggerSFFile = cms.string("EfficienciesAndSF_BCDEF_trigger.root"),
  isoSFFile = cms.string("EfficienciesAndSF_BCDEF_ISO.root"),
  idSFFile = cms.string("EfficienciesAndSF_BCDEF_ID.root"),
  trackingSFFile = cms.string("EfficienciesAndSF_BCDEFGH_Tracking.root"),
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
  TriggerResults = cms.InputTag("TriggerResults","","HLT"),
  TriggerObjects = cms.InputTag("selectedPatTrigger"),
  TrackRefitter = cms.InputTag("TrackRefitter"),
  primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
  Traj = cms.InputTag("TrackRefitter"),

  debug = cms.untracked.bool(False),
  onlyGenLevel = cms.bool(False),
  runMMC = cms.bool(False)
)
process.dump=cms.EDAnalyzer('EventContentAnalyzer')

process.TFileService = cms.Service("TFileService",
  fileName = cms.string("out_ana.root")
)

process.phlt = cms.Path(process.hltfilter)
process.pDiHiggsWWBBAna = cms.Path(
  process.eventCounterFilter*
  process.hltfilter*
  process.DiHiggsWWBBAna
)

process.pdump = cms.Path(process.dump)

process.schedule = cms.Schedule(process.pDiHiggsWWBBAna)
