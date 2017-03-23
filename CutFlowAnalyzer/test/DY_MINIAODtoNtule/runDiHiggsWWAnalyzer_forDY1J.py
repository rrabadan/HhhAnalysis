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
	'file:/fdata/hepx/store/user/lpernie/TEST_LOCALLY/DYJETS_7A385961-C6D9-E611-85B2-0025905B85BC.root'
    )
)

process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(100000) 
)

process.MessageLogger = cms.Service("MessageLogger", 
    destinations = cms.untracked.vstring("cout"), 
    cout = cms.untracked.PSet(threshold = cms.untracked.string("ERROR"))
)

process.DiHiggsWWBBAna = cms.EDAnalyzer('DiHiggsWWBBAnalyzer',
    verbose = cms.untracked.int32(0),
    #enum {Data = 0, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, TTbar, DYJets, DY0Jets, DY1Jets, DY2Jets};//add other background
    SampleType = cms.untracked.int32(16),
    #genParticles = cms.InputTag("genParticles"),
    genParticles = cms.InputTag("prunedGenParticles"),#minAOD
    #muons = cms.InputTag("cleanPatPFMuonsTriggerMatch"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    genjets = cms.InputTag("slimmedGenJets"), # For miniAOD
    #genjets = cms.InputTag("ak4GenJetsNoNu"),
    #genjets = cms.InputTag("ak4GenJets"),
    jets = cms.InputTag("slimmedJets"),
    mets = cms.InputTag("slimmedMETs"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    triggerEvent = cms.InputTag("patTriggerEvent"),
    tracks = cms.InputTag("generalTracks"),
    TriggerResults = cms.InputTag("TriggerResults","","RECO"),
    TrackRefitter = cms.InputTag("TrackRefitter"),
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    Traj = cms.InputTag("TrackRefitter"),
    finalStates = cms.bool(False),
    simulation = cms.bool(True),
    runMMC = cms.bool(False)
)
process.dump=cms.EDAnalyzer('EventContentAnalyzer')

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("out_ana.root")
)

process.pDiHiggsWWBBAna = cms.Path(
    process.DiHiggsWWBBAna
)

process.pdump = cms.Path(process.dump)

process.schedule = cms.Schedule(process.pDiHiggsWWBBAna)
