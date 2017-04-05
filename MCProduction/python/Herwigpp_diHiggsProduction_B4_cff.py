import FWCore.ParameterSet.Config as cms

from Configuration.Generator.HerwigppDefaults_cfi import *
from Configuration.Generator.HerwigppUE_EE_5C_cfi import *
from Configuration.Generator.HerwigppPDF_CTEQ6_LO_cfi import *
from Configuration.Generator.HerwigppEnergy_13TeV_cfi import *

generator = cms.EDFilter("ThePEGGeneratorFilter",
	herwigDefaultsBlock,
	herwigppUESettingsBlock,
	herwigppPDFSettingsBlock,
	herwigppEnergySettingsBlock,

	configFiles = cms.vstring(),
	parameterSets = cms.vstring(
		'hwpp_cmsDefaults',
		'hwpp_ue_EE5C',
		'hwpp_pdf_CTEQ6L1',
		'hwpp_cm_13TeV',
		'diHiggsProduction',

	),

	diHiggsProduction = cms.vstring(
		#############
		### Copy here the content of the Herwig++ input card, each line must be in quotation marks and separated by a comma.
		
		###############
    'cd /Herwig/Particles/',
    'create ThePEG::ParticleData boxon',
    'setup boxon 99926 boxon 0.0 0.0 0.0 0.0 0 0 0 1',
    'create ThePEG::ParticleData triangon',
    'setup triangon 99927 triangon 0.0 0.0 0.0 0.0 0 0 0 1',
    'create ThePEG::ParticleData h0',
    'setup h0 25 h0 125 0.004088 0.04088 0 0 0 1 0',
    'create ThePEG::ParticleData H',
    'setup H 35 H 415.488 1.593 15.93 0 0 0 1 0',
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
    'set /Herwig/Particles/W+/W+->nu_tau,tau+;:OnOff Off',#should keep tau or not ??
    'set /Herwig/Particles/W-/W-->b,cbar;:OnOff Off',
    'set /Herwig/Particles/W-/W-->cbar,d;:OnOff Off',
    'set /Herwig/Particles/W-/W-->cbar,s;:OnOff Off',
    'set /Herwig/Particles/W-/W-->s,ubar;:OnOff Off',
    'set /Herwig/Particles/W-/W-->ubar,d;:OnOff Off',
    'set /Herwig/Particles/W-/W-->nu_ebar,e-;:OnOff On',
    'set /Herwig/Particles/W-/W-->nu_mubar,mu-;:OnOff On',
    'set /Herwig/Particles/W-/W-->nu_taubar,tau-;:OnOff Off',

##################################################
# Matrix Elements for hadron-hadron collisions    
# gg->H->hh               
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
    'set MEHiggsPair:SelfCoupling 1.4152',
    'set MEHiggsPair:hhHCoupling 0.2354',


##################################################
# Save run for later usage with 'Herwig++ run'
##################################################
#     'cd /Herwig/Generators',
#     'saverun LHC-HH-bbWW-B3-LeptonW LHCGenerator'
	    ),

        crossSection = cms.untracked.double(-1),
        filterEfficiency = cms.untracked.double(1.0),
)

#w1  Candidate id: 24 mass: 62.4283 (P,E)= (-6.40087, 52.4484, -49.2489, 95.4702)(Pt,E) = (52.8375, -0.832524, 1.69224, 95.4702) status: 11
#w2  Candidate id: -24 mass: 59.3055 (P,E)= (15.3628, 53.444, -71.1402, 108.029)(Pt,E) = (55.6082, -1.06577, 1.29089, 108.029) status: 11
#b1  Candidate id: 5 mass: 29.5901 (P,E)= (-27.5392, -16.2606, -114.209, 122.238)(Pt,E) = (31.9815, -1.98507, -2.60822, 122.238) status: 11
#b2  Candidate id: -5 mass: 5 (P,E)= (32.998, -108.667, -72.0992, 134.613)(Pt,E) = (113.567, -0.598487, -1.27598, 134.613) status: 11
#process.bbWWFilter = cms.EDFilter("MCMultiParticleFilter",
bbWWFilter = cms.EDFilter("MCHhhMultiParticleFilter",
  src         = cms.untracked.InputTag("generator", "unsmeared"),
  debug       = cms.untracked.bool(False)
)


ProductionFilterSequence = cms.Sequence(generator * bbWWFilter)
