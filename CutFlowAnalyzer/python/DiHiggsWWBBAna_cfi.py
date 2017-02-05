import FWCore.ParameterSet.Config as cms

print "Deprecated"
DiHiggsWWBBAna = cms.EDAnalyzer('DiHiggsWWBBAnalyzer',
	verbose = cms.untracked.int32(0),
        #enum {Data = 0, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, TTbar};//add other background
	SampleType = cms.untracked.int32(0),
    	genParticles = cms.InputTag("genParticles"),
	muons = cms.InputTag("cleanPatPFMuonsTriggerMatch"),
	beamSpot = cms.InputTag("offlineBeamSpot"),
	triggerEvent = cms.InputTag("patTriggerEvent"),
	tracks = cms.InputTag("generalTracks"),
	TriggerResults = cms.InputTag("TriggerResults","","RECO"),
	TrackRefitter = cms.InputTag("TrackRefitter"),
	primaryVertices = cms.InputTag("offlinePrimaryVertices"),
	PATJet = cms.InputTag("patJets"),
	Traj = cms.InputTag("TrackRefitter"),
	finalStates = cms.bool(False),
	simulation = cms.bool(True),
        runMMC = cms.bool(False)
        )
