#this fake PSET is needed for local test and for crab to figure the output filename
#you do not need to edit it unless you want to do a local test using a different input file than
#the one marked below
import FWCore.ParameterSet.Config as cms
process = cms.Process('NANO')
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(),
#	lumisToProcess=cms.untracked.VLuminosityBlockRange("254231:1-254231:24")
)
process.source.fileNames = [
    	#'file:/fdata/hepx/store/user/tahuang/JPsiToMuMu_Pt20to100-pythia8-gun/RAW2DIGI_RECO_Muons_JPsiToMuMu_v2/170712_091000/0000/out_reco_1.root'
    	'file:/fdata/hepx/store/user/taohuang/HH_NanoAOD/nanolzma_1.root'
	#'../lzma_1.root' ##you can change only this line
	#'../../NanoAOD/test/lzma.root' ##you can change only this line
]
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))
process.output = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string('tree.root'))
process.out = cms.EndPath(process.output)

