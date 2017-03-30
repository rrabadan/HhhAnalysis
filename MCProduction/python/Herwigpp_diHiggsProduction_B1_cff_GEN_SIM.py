# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: HhhAnalysis/MCProduction/Herwigpp_diHiggsProduction_B1_cff -s GEN,SIM --mc --datatier GEN-SIM --beamspot Realistic25ns13TeV2016Collision --conditions 80X_mcRun2_asymptotic_2016_miniAODv2_v1 --eventcontent RAWSIM --era Run2_2016 --fileout file:out_sim.root -n 1 --no_exec
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
    input = cms.untracked.int32(1)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('HhhAnalysis/MCProduction/Herwigpp_diHiggsProduction_B1_cff nevts:1'),
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
    fileName = cms.untracked.string('file:out_sim.root'),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_miniAODv2_v1', '')

process.bbWWFilter = cms.EDFilter("MCHhhMultiParticleFilter",
    AcceptMore = cms.bool(True),
    EtaMax = cms.vdouble(10, 10, 10, 10),
    NumRequired = cms.int32(4),
    ParticleID = cms.vint32(5, -5, 24, -24),
    PtMin = cms.vdouble(0.0, 0.0, 0, 0),
    Status = cms.vint32(11, 11, 11, 11),
    src = cms.untracked.InputTag("generator","unsmeared")
)


process.generator = cms.EDFilter("ThePEGGeneratorFilter",
    configFiles = cms.vstring(),
    crossSection = cms.untracked.double(-1),
    dataLocation = cms.string('${HERWIGPATH}'),
    diHiggsProduction = cms.vstring('cd /Herwig/Particles/', 
        'create ThePEG::ParticleData boxon', 
        'setup boxon 99926 boxon 0.0 0.0 0.0 0.0 0 0 0 1', 
        'create ThePEG::ParticleData triangon', 
        'setup triangon 99927 triangon 0.0 0.0 0.0 0.0 0 0 0 1', 
        'create ThePEG::ParticleData h0', 
        'setup h0 25 h0 125 0.003196 0.03196 0 0 0 1 0', 
        'create ThePEG::ParticleData H', 
        'setup H 35 H 258.0126 0.68 6.8 0 0 0 1 0', 
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
        'set /Herwig/Particles/h0/h0->b,bbar;:OnOff On', 
        'set /Herwig/Particles/h0/h0->W+,W-;:OnOff On', 
        'set /Herwig/Particles/h0/h0->tau-,tau+;:OnOff Off', 
        'set /Herwig/Particles/h0/h0->g,g;:OnOff Off', 
        'set /Herwig/Particles/h0/h0->c,cbar;:OnOff Off', 
        'set /Herwig/Particles/h0/h0->Z0,Z0;:OnOff Off', 
        'set /Herwig/Particles/h0/h0->gamma,gamma;:OnOff Off', 
        'set /Herwig/Particles/h0/h0->mu-,mu+;:OnOff Off', 
        'set /Herwig/Particles/h0/h0->t,tbar;:OnOff Off', 
        'set /Herwig/Particles/W+/W+->bbar,c;:OnOff Off', 
        'set /Herwig/Particles/W+/W+->c,dbar;:OnOff Off', 
        'set /Herwig/Particles/W+/W+->c,sbar;:OnOff Off', 
        'set /Herwig/Particles/W+/W+->sbar,u;:OnOff Off', 
        'set /Herwig/Particles/W+/W+->u,dbar;:OnOff Off', 
        'set /Herwig/Particles/W+/W+->nu_e,e+;:OnOff On', 
        'set /Herwig/Particles/W+/W+->nu_mu,mu+;:OnOff On', 
        'set /Herwig/Particles/W+/W+->nu_tau,tau+;:OnOff Off', 
        'set /Herwig/Particles/W-/W-->b,cbar;:OnOff Off', 
        'set /Herwig/Particles/W-/W-->cbar,d;:OnOff Off', 
        'set /Herwig/Particles/W-/W-->cbar,s;:OnOff Off', 
        'set /Herwig/Particles/W-/W-->s,ubar;:OnOff Off', 
        'set /Herwig/Particles/W-/W-->ubar,d;:OnOff Off', 
        'set /Herwig/Particles/W-/W-->nu_ebar,e-;:OnOff On', 
        'set /Herwig/Particles/W-/W-->nu_mubar,mu-;:OnOff On', 
        'set /Herwig/Particles/W-/W-->nu_taubar,tau-;:OnOff Off', 
        'cd /', 
        'cd /Herwig/MatrixElements/', 
        'library ./MEHiggsPair.so', 
        'create Herwig::MEHiggsPair MEHiggsPair MEHiggsPair.so', 
        'insert SimpleQCD:MatrixElements[0] MEHiggsPair', 
        'set MEHiggsPair:Process ggToHTohh', 
        'set MEHiggsPair:SelfCoupling 3.4387', 
        'set MEHiggsPair:hhHCoupling -0.7025'),
    eventHandlers = cms.string('/Herwig/EventHandlers'),
    filterEfficiency = cms.untracked.double(1.0),
    generatorModule = cms.string('/Herwig/Generators/LHCGenerator'),
    hwpp_basicSetup = cms.vstring('create ThePEG::RandomEngineGlue /Herwig/RandomGlue', 
        'set /Herwig/Generators/LHCGenerator:RandomNumberGenerator /Herwig/RandomGlue', 
        'set /Herwig/Generators/LHCGenerator:NumberOfEvents 10000000', 
        'set /Herwig/Generators/LHCGenerator:DebugLevel 1', 
        'set /Herwig/Generators/LHCGenerator:UseStdout 0', 
        'set /Herwig/Generators/LHCGenerator:PrintEvent 0', 
        'set /Herwig/Generators/LHCGenerator:MaxErrors 10000'),
    hwpp_cm_13TeV = cms.vstring('set /Herwig/Generators/LHCGenerator:EventHandler:LuminosityFunction:Energy 13000.0', 
        'set /Herwig/Shower/Evolver:IntrinsicPtGaussian 2.2*GeV'),
    hwpp_cmsDefaults = cms.vstring('+hwpp_basicSetup', 
        '+hwpp_setParticlesStableForDetector'),
    hwpp_pdf_CTEQ6L1 = cms.vstring('+hwpp_pdf_CTEQ6L1_Common', 
        '+hwpp_ue_EE5C'),
    hwpp_pdf_CTEQ6L1_CUETHS1 = cms.vstring('+hwpp_pdf_CTEQ6L1_Common', 
        '+hwpp_ue_CUETHS1'),
    hwpp_pdf_CTEQ6L1_Common = cms.vstring('create ThePEG::LHAPDF /Herwig/Partons/cmsPDFSet ThePEGLHAPDF.so', 
        'set /Herwig/Partons/cmsPDFSet:PDFName cteq6ll.LHpdf', 
        'set /Herwig/Partons/cmsPDFSet:RemnantHandler /Herwig/Partons/HadronRemnants', 
        'set /Herwig/Particles/p+:PDF /Herwig/Partons/cmsPDFSet', 
        'set /Herwig/Particles/pbar-:PDF /Herwig/Partons/cmsPDFSet'),
    hwpp_pdf_CTEQ6L1_Hard = cms.vstring('+hwpp_pdf_CTEQ6L1_Hard_Common', 
        '+hwpp_ue_EE5C'),
    hwpp_pdf_CTEQ6L1_Hard_CUETHS1 = cms.vstring('+hwpp_pdf_CTEQ6L1_Hard_Common', 
        '+hwpp_ue_CUETHS1'),
    hwpp_pdf_CTEQ6L1_Hard_Common = cms.vstring('create ThePEG::LHAPDF /Herwig/Partons/cmsHardPDFSet ThePEGLHAPDF.so', 
        'set /Herwig/Partons/cmsHardPDFSet:PDFName cteq6ll.LHpdf', 
        'set /Herwig/Partons/cmsHardPDFSet:RemnantHandler /Herwig/Partons/HadronRemnants'),
    hwpp_pdf_CTEQ6LL = cms.vstring('+hwpp_pdf_CTEQ6L1'),
    hwpp_pdf_CTEQ6LL_CUETHS1 = cms.vstring('+hwpp_pdf_CTEQ6L1_CUETHS1'),
    hwpp_pdf_CTEQ6LL_Hard = cms.vstring('+hwpp_pdf_CTEQ6L1_Hard'),
    hwpp_pdf_CTEQ6LL_Hard_CUETHS1 = cms.vstring('+hwpp_pdf_CTEQ6L1_Hard_CUETHS1'),
    hwpp_setParticlesStableForDetector = cms.vstring('set /Herwig/Particles/mu-:Stable Stable', 
        'set /Herwig/Particles/mu+:Stable Stable', 
        'set /Herwig/Particles/Sigma-:Stable Stable', 
        'set /Herwig/Particles/Sigmabar+:Stable Stable', 
        'set /Herwig/Particles/Lambda0:Stable Stable', 
        'set /Herwig/Particles/Lambdabar0:Stable Stable', 
        'set /Herwig/Particles/Sigma+:Stable Stable', 
        'set /Herwig/Particles/Sigmabar-:Stable Stable', 
        'set /Herwig/Particles/Xi-:Stable Stable', 
        'set /Herwig/Particles/Xibar+:Stable Stable', 
        'set /Herwig/Particles/Xi0:Stable Stable', 
        'set /Herwig/Particles/Xibar0:Stable Stable', 
        'set /Herwig/Particles/Omega-:Stable Stable', 
        'set /Herwig/Particles/Omegabar+:Stable Stable', 
        'set /Herwig/Particles/pi+:Stable Stable', 
        'set /Herwig/Particles/pi-:Stable Stable', 
        'set /Herwig/Particles/K+:Stable Stable', 
        'set /Herwig/Particles/K-:Stable Stable', 
        'set /Herwig/Particles/K_S0:Stable Stable', 
        'set /Herwig/Particles/K_L0:Stable Stable'),
    hwpp_ue_EE5C = cms.vstring('+hwpp_ue_EE5CEnergyExtrapol', 
        'set /Herwig/Hadronization/ColourReconnector:ColourReconnection Yes', 
        'set /Herwig/Hadronization/ColourReconnector:ReconnectionProbability 0.49', 
        'set /Herwig/Partons/RemnantDecayer:colourDisrupt 0.80', 
        'set /Herwig/UnderlyingEvent/MPIHandler:InvRadius 2.30', 
        'set /Herwig/UnderlyingEvent/MPIHandler:softInt Yes', 
        'set /Herwig/UnderlyingEvent/MPIHandler:twoComp Yes', 
        'set /Herwig/UnderlyingEvent/MPIHandler:DLmode 2'),
    hwpp_ue_EE5CEnergyExtrapol = cms.vstring('set /Herwig/UnderlyingEvent/MPIHandler:EnergyExtrapolation Power', 
        'set /Herwig/UnderlyingEvent/MPIHandler:ReferenceScale 7000.*GeV', 
        'set /Herwig/UnderlyingEvent/MPIHandler:Power 0.33', 
        'set /Herwig/UnderlyingEvent/MPIHandler:pTmin0 3.91*GeV'),
    parameterSets = cms.vstring('hwpp_cmsDefaults', 
        'hwpp_ue_EE5C', 
        'hwpp_pdf_CTEQ6L1', 
        'hwpp_cm_13TeV', 
        'diHiggsProduction'),
    repository = cms.string('HerwigDefaults.rpo'),
    run = cms.string('LHC')
)


process.ProductionFilterSequence = cms.Sequence(process.generator+process.bbWWFilter)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.RAWSIMoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 


