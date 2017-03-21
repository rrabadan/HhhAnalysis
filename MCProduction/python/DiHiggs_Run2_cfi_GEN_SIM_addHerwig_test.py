# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: DiHiggs_Run2_cfi -s GEN,SIM --mc --datatier GEN-SIM --beamspot Realistic25ns13TeV2016Collision --conditions 80X_mcRun2_asymptotic_2016_miniAODv2_v1 --eventcontent RAWSIM --era Run2_2016 --filetype LHE --filein file:output_from_herwig.hepmc --fileout file:out_sim.root -n 10 --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('SIM',eras.Run2_2016)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.Geometry.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic25ns13TeV2016Collision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

# Input source
process.source = cms.Source("EmptySource")
#process.source = cms.Source("MCFileSource",
#    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
#    #fileNames = cms.untracked.vstring('file:output_from_herwig.hepmc')
#   fileNames = cms.untracked.vstring('file:/fdata/hepx/store/user/taohuang/Pheno/HH-bbWW-B6-20160518-leptonW-1000000.hepmc')
    #fileNames = cms.untracked.vstring('file:/home/taohuang/DiHiggsAnalysis/CMSSW_8_0_19/src/Crabjobs/HH-bbWW-B3-13TeV-leptonW-OnlyME-10k.hepmc')
    #fileNames = cms.untracked.vstring('file:HH-bbWW-B3-NoHadronization-1k.hepmc')
   
#)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('DiHiggs_Run2_cfi nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('file:/fdata/hepx/store/user/taohuang/DiHiggsAnalysisSample/out_sim_hadronization_10k.root'),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_miniAODv2_v1', '')

from Configuration.Generator.HerwigppDefaults_cfi import *
from Configuration.Generator.HerwigppUE_EE_5C_cfi import *
from Configuration.Generator.HerwigppPDF_CTEQ6_LO_cfi import *
from Configuration.Generator.HerwigppEnergy_13TeV_cfi import *

process.generator = cms.EDFilter("ThePEGGeneratorFilter",
	herwigDefaultsBlock,
	herwigppUESettingsBlock,
	crossSection = cms.untracked.double(2.0e+10),
	filterEfficiency = cms.untracked.double(1.0),

	configFiles = cms.vstring(),
	parameterSets = cms.vstring(
		'productionParameters',
		'hwpp_cmsDefaults',
		#'hwpp_ue_EE5C',
		#'hwpp_pdf_CTEQ6L1',
		#'hwpp_cm_13TeV',
		'hwpp_ue_EE5C',
	),
	productionParameters = cms.vstring(
		
		#############
		### Copy here the content of the Herwig++ input card, each line must be in quotation marks and separated by a comma.
		
		###############
		
    'cd /Herwig/Generators',
    'set LHCGenerator:NumberOfEvents 1000',
    'set LHCGenerator:RandomNumberGenerator:Seed 31122001',
    'set LHCGenerator:DebugLevel 1',
    'set LHCGenerator:PrintEvent 10',
    'set LHCGenerator:MaxErrors 10000',                

    'cd /Herwig/Particles/',
    'create ThePEG::ParticleData boxon',
    'setup boxon 99926 boxon 0.0 0.0 0.0 0.0 0 0 0 1',
    'create ThePEG::ParticleData triangon',
    'setup triangon 99927 triangon 0.0 0.0 0.0 0.0 0 0 0 1',
    'create ThePEG::ParticleData h0',
    'setup h0 25 h0 125 0.003196 0.03196 0 0 0 1 0',
    'create ThePEG::ParticleData H',
    'setup H 35 H 352.75 2.17 21.7 0 0 0 1 0',
    'create ThePEG::ParticleData t',
    'setup t 6 t 174.3 1.4 14 0 2 3 2 0', 
    'create ThePEG::ParticleData tbar',
    'setup tbar -6 tbar 174.3 1.4 14 0 -2 -3 2 0', 
    'makeanti tbar t',
    'create ThePEG::ConstituentParticleData b',
    'setup b 5 b 4.7 0.  0 0 -1 3 2 1  5',
    'create ThePEG::ConstituentParticleData bbar',
    'setup bbar -5 bbar 4.7 0.  0 0 1 -3 2 1  5',
    'makeanti bbar b',
    # Disable/Enable Higgs decays
    'set /Herwig/Particles/h0/h0->b,bbar;:OnOff On',
    'set /Herwig/Particles/h0/h0->W+,W-;:OnOff On',
    'set /Herwig/Particles/h0/h0->tau-,tau+;:OnOff Off',
    'set /Herwig/Particles/h0/h0->g,g;:OnOff Off',
    'set /Herwig/Particles/h0/h0->c,cbar;:OnOff Off',
    'set /Herwig/Particles/h0/h0->Z0,Z0;:OnOff Off',
    'set /Herwig/Particles/h0/h0->gamma,gamma;:OnOff Off',
    'set /Herwig/Particles/h0/h0->mu-,mu+;:OnOff Off',
    'set /Herwig/Particles/h0/h0->t,tbar;:OnOff Off',
    #W decay
    #W+->bbar,c;
    #W+->c,dbar;
    #W+->c,sbar;
    #W+->nu_e,e+;
    #W+->nu_mu,mu+;
    #W+->nu_tau,tau+;
    #W+->sbar,u;
    #W+->u,dbar
    'set /Herwig/Particles/W+/W+->bbar,c;:OnOff Off',
    'set /Herwig/Particles/W+/W+->c,dbar;:OnOff Off',
    'set /Herwig/Particles/W+/W+->c,sbar;:OnOff Off',
    'set /Herwig/Particles/W+/W+->sbar,u;:OnOff Off',
    'set /Herwig/Particles/W+/W+->u,dbar;:OnOff Off',
    'set /Herwig/Particles/W+/W+->nu_e,e+;:OnOff On',
    'set /Herwig/Particles/W+/W+->nu_mu,mu+;:OnOff On',
    'set /Herwig/Particles/W+/W+->nu_tau,tau+;:OnOff On',
    'set /Herwig/Particles/W-/W-->b,cbar;:OnOff Off',
    'set /Herwig/Particles/W-/W-->cbar,d;:OnOff Off',
    'set /Herwig/Particles/W-/W-->cbar,s;:OnOff Off',
    'set /Herwig/Particles/W-/W-->s,ubar;:OnOff Off',
    'set /Herwig/Particles/W-/W-->ubar,d;:OnOff Off',
    'set /Herwig/Particles/W-/W-->nu_ebar,e-;:OnOff On',
    'set /Herwig/Particles/W-/W-->nu_mubar,mu-;:OnOff On',
    'set /Herwig/Particles/W-/W-->nu_taubar,tau-;:OnOff On',
    #
    # Set B-mesons stable (?)
    #set /Herwig/Particles/B+:Stable Stable
    #set /Herwig/Particles/B-:Stable Stable
    #set /Herwig/Particles/B0:Stable Stable
    #set /Herwig/Particles/Bbar0:Stable Stable
    #set /Herwig/Particles/B_s0:Stable Stable
    #set /Herwig/Particles/B_sbar0:Stable Stable
    #set /Herwig/Particles/B_c+:Stable Stable
    #set /Herwig/Particles/B_c-:Stable Stable
    #set /Herwig/Particles/Upsilon:Stable Stable

    # Set b-baryons stable
    #set /Herwig/Particles/Sigma_b+:Stable Stable
    #set /Herwig/Particles/Lambda_b0:Stable Stable
    #set /Herwig/Particles/Sigma_b-:Stable Stable
    #set /Herwig/Particles/Xi_b0:Stable Stable
    #set /Herwig/Particles/Xi_b-:Stable Stable
    #set /Herwig/Particles/Omega_b-:Stable Stable
    #set /Herwig/Particles/Sigma_bbar-:Stable Stable
    #set /Herwig/Particles/Lambda_bbar0:Stable Stable
    #set /Herwig/Particles/Sigma_bbar+:Stable Stable
    #set /Herwig/Particles/Xi_bbar0:Stable Stable
    #set /Herwig/Particles/Xi_bbar+:Stable Stable
    #set /Herwig/Particles/Omega_bbar+:Stable Stable

    #set /Herwig/Particles/h0:Stable Stable

    'create ThePEG::LHAPDF /Herwig/Partons/cmsPDFSet ThePEGLHAPDF.so',	# cmsPDFSet Default name for shower PDF
    'set /Herwig/Partons/cmsPDFSet:PDFName cteq6ll.LHpdf',	
    'set /Herwig/Partons/cmsPDFSet:RemnantHandler /Herwig/Partons/HadronRemnants',
    'set /Herwig/Particles/p+:PDF /Herwig/Partons/cmsPDFSet',		# Use PDF in shower
    'set /Herwig/Particles/pbar-:PDF /Herwig/Partons/cmsPDFSet',

    # Colour reconnection settings
    'set /Herwig/Hadronization/ColourReconnector:ColourReconnection Yes',
    'set /Herwig/Hadronization/ColourReconnector:ReconnectionProbability 0.49',
    # Colour Disrupt settings
    'set /Herwig/Partons/RemnantDecayer:colourDisrupt 0.80',
    # inverse hadron radius
    'set /Herwig/UnderlyingEvent/MPIHandler:InvRadius 2.30',
    # MPI model settings
    'set /Herwig/UnderlyingEvent/MPIHandler:softInt Yes',
    'set /Herwig/UnderlyingEvent/MPIHandler:twoComp Yes',
    'set /Herwig/UnderlyingEvent/MPIHandler:DLmode 2',
    'set /Herwig/UnderlyingEvent/MPIHandler:EnergyExtrapolation Power',
    'set /Herwig/UnderlyingEvent/MPIHandler:ReferenceScale 7000.*GeV',
    'set /Herwig/UnderlyingEvent/MPIHandler:Power 0.33',
    'set /Herwig/UnderlyingEvent/MPIHandler:pTmin0 3.91*GeV',

    'set /Herwig/Generators/LHCGenerator:EventHandler:LuminosityFunction:Energy 13000.0',
    'set /Herwig/Shower/Evolver:IntrinsicPtGaussian 2.2*GeV',

    'insert LHCGenerator:AnalysisHandlers 0 /Herwig/Analysis/HepMCFile',
    'set /Herwig/Analysis/HepMCFile:PrintEvent 100'
    'set /Herwig/Analysis/HepMCFile:Format GenEvent'
    'set /Herwig/Analysis/HepMCFile:Units GeV_mm'

##################################################
# Matrix Elements for hadron-hadron collisions    
# PLUGIN: gg -> Jpsi Jpsi                         
##################################################
    'cd /',
    'cd /Herwig/MatrixElements/',
    #'library ./MEgg2JpsiJpsi.so',
    #'create PLUGIN::MEgg2JpsiJpsi MEgg2JpsiJpsi',
    'library ./MEHiggsPair.so',
    #'create PLUGIN::MEHiggsPair MEHiggsPair',
    'create Herwig::MEHiggsPair MEHiggsPair MEHiggsPair.so',
    'insert SimpleQCD:MatrixElements[0] MEHiggsPair',
    'set MEHiggsPair:Process ggToHTohh',
    #set MEHiggsPair:Process All
    'set MEHiggsPair:SelfCoupling 7.1514',
    'set MEHiggsPair:hhHCoupling -0.5681',

##################################################
# Change default parton level kinematic cuts
# (default in share/Herwig++/defaults/Cuts.in)
##################################################
     'set /Herwig/Cuts/MassCut:MinM 2.0*GeV',
     'set /Herwig/Cuts/MassCut:MaxM 6.0*GeV',
     'set /Herwig/Cuts/QCDCuts:MHatMin 0.0*GeV',
     'set /Herwig/Cuts/QCDCuts:MHatMax 1000.0*GeV',
     'set /Herwig/Cuts/QCDCuts:X1Min 1.0e-8',
     'set /Herwig/Cuts/QCDCuts:X2Min 1.0e-8',
     'set /Herwig/UnderlyingEvent/MPIHandler:IdenticalToUE 0',

# Save run for later usage with 'Herwig++ run'
##################################################
     'cd /Herwig/Generators',
     'saverun LHC-HH-bbWW-B3-LeptonW LHCGenerator'
	)
)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.RAWSIMoutput_step)
for path in process.paths:
        getattr(process,path)._seq = process.generator *  getattr(process,path)._seq 
        #getattr(process,path)._seq = process.generator * process.fourMuonFilter  *  getattr(process,path)._seq 



"""
def MassReplaceInputTag(aProcess,oldT="rawDataCollector",newT="rawDataRepacker"):
	from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag
	for s in aProcess.paths_().keys():
		massSearchReplaceAnyInputTag(getattr(aProcess,s),oldT,newT)
MassReplaceInputTag(process, "generator", "source")
process.VtxSmeared.src = "source"
"""
