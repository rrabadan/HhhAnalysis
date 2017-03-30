import FWCore.ParameterSet.Config as cms

#from Configuration.Generator.PythiaUEZ2starSettings_cfi import *
#from Configuration.GenProduction.PythiaUESettings_cfi import *

process = cms.Process("TEST")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
##process.MessageLogger.destinations = ['cerr']
##process.MessageLogger.statistics = []
##process.MessageLogger.fwkJobReports = []
##process.MessageLogger.cerr.threshold = "Warning"


process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
    fileNames = cms.untracked.vstring('file:out_sim_generatoronly_100.root')
    #fileNames = cms.untracked.vstring('file:/store/user/tahuang/DiHiggs/out_sim.root'),
    )

#w1  Candidate id: 24 mass: 62.4283 (P,E)= (-6.40087, 52.4484, -49.2489, 95.4702)(Pt,E) = (52.8375, -0.832524, 1.69224, 95.4702) status: 11
#w2  Candidate id: -24 mass: 59.3055 (P,E)= (15.3628, 53.444, -71.1402, 108.029)(Pt,E) = (55.6082, -1.06577, 1.29089, 108.029) status: 11
#b1  Candidate id: 5 mass: 29.5901 (P,E)= (-27.5392, -16.2606, -114.209, 122.238)(Pt,E) = (31.9815, -1.98507, -2.60822, 122.238) status: 11
#b2  Candidate id: -5 mass: 5 (P,E)= (32.998, -108.667, -72.0992, 134.613)(Pt,E) = (113.567, -0.598487, -1.27598, 134.613) status: 11
process.bbWWFilter = cms.EDFilter("MCHhhMultiParticleFilter",
  src         = cms.untracked.InputTag("generator", "unsmeared"),
  #src         = cms.InputTag("generatorSmeared"),
  NumRequired = cms.int32(4), #
  AcceptMore  = cms.bool(True),
  Status      = cms.vint32(11, 11, 11, 11),
  ParticleID  = cms.vint32(5, -5, 24, -24),
  EtaMax      = cms.vdouble(10, 10, 10, 10), # no cut on eta
  PtMin       = cms.vdouble(0., 0.,0, 0)
)

process.filter_step = cms.Path(process.bbWWFilter)
process.output = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
       SelectEvents = cms.vstring('filter_step')
       ),
    fileName = cms.untracked.string('out_filter.root'),
    outputCommands = cms.untracked.vstring('keep *'),
    splitLevel = cms.untracked.int32(0)
)

process.outpath = cms.EndPath(process.output)
process.schedule = cms.Schedule( process.filter_step, process.outpath)
