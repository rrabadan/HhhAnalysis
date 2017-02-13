import FWCore.ParameterSet.Config as cms

process = cms.Process("PATANA")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")

# verbose flags for the PF2PAT modules
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(False)

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("HhhAnalysis.DataFormats.miniAODtoPAT_cff")
process.load("HhhAnalysis.CutFlowAnalyzer.CutFlowAnalyzer_cff")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'file:022B1BC6-A789-E611-808B-B499BAAB427C.root'
	'/store/mc/RunIISpring16MiniAODv1/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext3-v1/00000/02F64C80-990E-E611-A2FE-842B2B185476.root'
    )
)

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('patTuple.root')
)

### Add MuJet Dataformats
from HhhAnalysis.DataFormats.EventContent_cff import *
process = customizePatOutput(process)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.p = cms.Path(
#    process.patifyMC *
    process.patifyData 
#   * process.cutFlowAnalyzers
)

process.outpath = cms.EndPath(process.out)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("out_ana.root")
)
