// -*- C++ -*-
//
// Package:    DiHiggsWWBBAnalyzer
// Class:      DiHiggsWWBBAnalyzer
// 
/**\class DiHiggsWWBBAnalyzer DiHiggsWWBBAnalyzer.cc DiHiggsWW/DiHiggsWWBBAnalyzer/plugins/DiHiggsWWBBAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  tao huang
//         Created:  Wed, 26 Nov 2014 17:58:07 GMT
// $Id$
//
//


// system include files
#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "math.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//headers from root lib
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TH1F.h"



//#include "HhhAnalysis/CutFlowAnalyzer/src/MMC.h"

float dxy(const reco::Candidate *cand, const reco::Vertex *point){
    return (-(cand->vx()-point->x())*cand->py()+(cand->vy()+point->y())*cand->px())/cand->pt();
};

float dz(const reco::Candidate *cand, const reco::Vertex *point){
    return ((cand->vz() - point->z()) - ((cand->vx() - point->x()) * cand->px() + (cand->vy() - point->y()) * cand->py()) /cand->pt()) * cand->pz() /cand->pt();
};

class MMC;
//#define WMass 80.385   // W mass
//#define SMHMass 125.03 // SM module higgs mass


//
// class declaration
//

using namespace reco;

typedef std::pair<float, float> EtaPhi;

class DiHiggsWWBBAnalyzer : public edm::EDAnalyzer {
   public:
      explicit DiHiggsWWBBAnalyzer(const edm::ParameterSet&);
      ~DiHiggsWWBBAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
    
   private:
      edm::ParameterSet cfg_;
      edm::ParameterSet mmcset_;
      // Labels to access
      edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;        // reconstructed muons
      edm::EDGetTokenT<pat::ElectronCollection> electronToken_;        // reconstructed electron
      edm::EDGetTokenT<pat::JetCollection> jetToken_;
      edm::EDGetTokenT<std::vector<reco::GenJet>> genjetToken_;
      edm::EDGetTokenT<pat::METCollection> metToken_;
      edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
      edm::EDGetTokenT<pat::TriggerEvent> triggerEventToken_;
      edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
      edm::EDGetTokenT<edm::TriggerResults> trigResToken_;
      edm::EDGetTokenT<reco::TrackCollection> trackRefToken_;
      edm::EDGetTokenT< std::vector<Trajectory> > trajToken_;
      edm::EDGetTokenT<reco::VertexCollection> primaryVerticesToken_;
      //configuration
      //data: 0 
      //or MC:signal(B1-B12)+background
      enum {Data = 0, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, TTbar};//add other background
      int sampleType_;
      //bool runMCMatch;//select physics objects by matching them to gen information
      int jetSelectionAlgo_;
      int muonSelectionAlgo_;
      int jetCorrectionAlgo_;
      int metCorrectionAlgo_;
      float mu_eta_;
      float el_eta_;
      float leadingpt_mumu_;//mu-mu events, leading
      float trailingpt_mumu_;//
      float leadingpt_muel_;//
      float trailingpt_muel_;//
      float leadingpt_elmu_;//
      float trailingpt_elmu_;//
      float leadingpt_elel_;
      float trailingpt_elel_;
      float jet_eta_;
      float jet_leadingpt_;
      float jet_trailingpt_;
      std::string bjetDiscrName_;
      float bjetDiscrCut_loose_;
      float bjetDiscrCut_medium_;
      float bjetDiscrCut_tight_;
      float met_;
      float jetleptonDeltaR_;
      float muIso_;
      float elIso_;
      float iterations_;
      //gen matching 
      float leptonsDeltaR_;//dR(gen, reco)
      float jetsDeltaR_;//gen matching 
      bool finalStates_;
      
      // debuglevel constrol 
      int verbose_; 
      void print();
      void printHtoWWChain();
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    private:
      TTree *evtree;
      TFile *output;
      void initBranches();
      void fillbranches(); 
      edm::Service< TFileService > fs;
      
      //----------branches of tree ---------------------------
      int ievent;
   //******************************
   //GEN Level
   //******************************
   private:
      const reco::Candidate* findmudaughter(const reco::Candidate* );
      const reco::Candidate* findnudaughter(const reco::Candidate* );
      const reco::Candidate* findmudescendants(const reco::Candidate*, int& );
      const reco::Candidate* findnudescendants(const reco::Candidate*, int& );
     //decendants and ancestor
      const reco::Candidate* stabledecendant(const reco::Candidate* cand, int id);
     //const reco::Candidate* stablehtoWWdecendant(Particle p, PdgId id);
      const reco::Candidate* finddecendant(const reco::Candidate* cand, int id, bool first=false);
      const reco::Candidate* findancestor(const reco::Candidate* cand, int id, bool first=false);
      bool hasMother(const reco::Candidate* cand, int id);
      bool hasDaughter(const reco::Candidate* cand, int id);
      
      void printCandidate(const reco::Candidate* );
      void printallDecendants(const reco::Candidate* );
      void printallAncestors(const reco::Candidate* );
      void checkGenParticlesSignal(edm::Handle<reco::GenParticleCollection> genParticleColl);
      void checkGenParticlesTTbar(edm::Handle<reco::GenParticleCollection> genParticleColl);
      void matchGenJet2Parton(edm::Handle<std::vector<reco::GenJet>> genjetColl);
      void matchmuon2Gen();//match pat:::Muon to gen muon 
      void matchBjets2Gen();//match genjet to gen b and then match pat::Jet to genjet


    private:
      // ---------- Candidates in signal channel ---------------------------
      const reco::Candidate* mu1_W1_cand;
      const reco::Candidate* nu1_W1_cand;
      const reco::Candidate* mu2_W2_cand;
      const reco::Candidate* nu2_W2_cand;
      const reco::Candidate* mu1_htoWW_cand;
      const reco::Candidate* mu2_htoWW_cand;
      const reco::Candidate* b1_htobb_cand;
      const reco::Candidate* b2_htobb_cand;
      const reco::Candidate* h2tohh_cand;

      const reco::Candidate* mu1cand;
      const reco::Candidate* nu1cand;
      const reco::Candidate* mu2cand;
      const reco::Candidate* nu2cand;
      const reco::Candidate* b1cand;
      const reco::Candidate* b2cand;
      const reco::Candidate* w1cand;
      const reco::Candidate* w2cand;
      const reco::Candidate* htoWWcand;
      const reco::Candidate* htoBBcand;
      const reco::Candidate* h2tohhcand;
      const reco::Candidate* t1cand;
      const reco::Candidate* t2cand;
      const reco::GenJet* b1genjet;
      const reco::GenJet* b2genjet;



    private:
      
      TLorentzVector mu1_lorentz;
      TLorentzVector mu2_lorentz;
      TLorentzVector bbar_lorentz;
      TLorentzVector nu1_lorentz;
      TLorentzVector nu2_lorentz;
      TLorentzVector met_lorentz;
      TLorentzVector bjet_lorentz;
      TLorentzVector bbarjet_lorentz;

    private:
      TLorentzVector stableDecendantsLorentz(const reco::Candidate* cand); 
      TLorentzVector calculateMET(); 

    private:
      //gen particles
      bool findAllGenParticles;
      float mu1_energy;
      float mu1_pt;
      float mu1_px;
      float mu1_py;
      float mu1_pz;
      float mu1_eta;
      float mu1_phi;
      float w1_energy;
      float w1_pt;
      float w1_eta;
      float w1_phi;
      float w1_px;
      float w1_py;
      float w1_pz;
      float w1_mass;
      float nu1_energy;
      float nu1_pt;
      float nu1_px;
      float nu1_py;
      float nu1_pz;
      float nu1_eta;
      float nu1_phi;
      bool Wtomu1nu1;

      float mu2_energy;
      float mu2_pt;
      float mu2_px;
      float mu2_py;
      float mu2_pz;
      float mu2_eta;
      float mu2_phi;
      float w2_energy;
      float w2_pt;
      float w2_eta;
      float w2_phi;
      float w2_px;
      float w2_py;
      float w2_pz;
      float w2_mass;
      float nu2_energy;
      float nu2_pt;
      float nu2_px;
      float nu2_py;
      float nu2_pz;
      float nu2_eta;
      float nu2_phi;
      bool Wtomu2nu2;
     
      float htoWW_energy;
      float htoWW_pt;
      float htoWW_eta;
      float htoWW_phi;
      float htoWW_px;
      float htoWW_py;
      float htoWW_pz; 
      float htoWW_mass;

      float b1_energy;
      float b1_pt;
      float b1_eta;
      float b1_phi;
      float b1_px;
      float b1_py;
      float b1_pz;
      float b2_energy;
      float b2_pt;
      float b2_eta;
      float b2_phi;
      float b2_px;
      float b2_py;
      float b2_pz;

      float htobb_energy;
      float htobb_pt;
      float htobb_eta;
      float htobb_phi;
      float htobb_px;
      float htobb_py;
      float htobb_pz;
      float htobb_mass;
      //genjet
      float b1genjet_energy;
      float b1genjet_pt;
      float b1genjet_eta;
      float b1genjet_phi;
      float b1genjet_px;
      float b1genjet_py;
      float b1genjet_pz;
      float b1genjet_mass;
      float dR_b1genjet;
      float b2genjet_energy;
      float b2genjet_pt;
      float b2genjet_eta;
      float b2genjet_phi;
      float b2genjet_px;
      float b2genjet_py;
      float b2genjet_pz;
      float b2genjet_mass;
      float dR_b2genjet;
      bool hastwogenjets;


      float genmet_pt;
      float genmet_phi;
      float genmet_px;
      float genmet_py;
      
      float h2tohh_energy;
      float h2tohh_px;
      float h2tohh_py;
      float h2tohh_pz;
      float h2tohh_mass;

      //for ttbar
      float t1_energy;
      float t1_px;
      float t1_py;
      float t1_pz;
      float t1_mass;
      float t2_energy;
      float t2_px;
      float t2_py;
      float t2_pz;
      float t2_mass;


      float dR_genbl;
      float dR_genb1l1;
      float dR_genb1l2;
      float dR_genb2l1;
      float dR_genb2l2;
      float dR_genl1l2;
      float dR_genl1l2b1b2;
      float dphi_genl1l2b1b2;
      float dR_genb1b2;
      float dR_genminbl;
      float mass_genl1l2;
      float energy_genl1l2;
      float pt_genl1l2;
      float phi_genl1l2;
      float eta_genl1l2;
      float mass_genb1b2;
      float energy_genb1b2;
      float pt_genb1b2;
      float phi_genb1b2;
      float eta_genb1b2;
      float dphi_genllbb;
      float dphi_genllmet;
      float mass_gentrans;

      //reco leve
      float numberOfmuon1;
      float numberOfmuon2;
      float muon1_energy;
      float muon1_pt;
      float muon1_eta;
      float muon1_phi;
      float muon1_px;
      float muon1_py;
      float muon1_pz;
      float muon1_isoVar;
      float muon1_dxy;
      float muon1_dz;
      //muonid loose, medium, tight
      int muon1_id;
      float dR_mu1;
      float muon2_energy;
      float muon2_pt;
      float muon2_eta;
      float muon2_phi;
      float muon2_px;
      float muon2_py;
      float muon2_pz;
      float muon2_isoVar;
      float muon2_dxy;
      float muon2_dz;
      int muon2_id;
      float dR_mu2;
      bool hastwomuons;
    
      float dR_b1jet;
      float dR_b2jet;
      float b1jet_px;
      float b1jet_py;
      float b1jet_pz;
      float b1jet_eta;
      float b1jet_phi;
      float b1jet_pt;
      float b1jet_energy;
      float b1jet_mass;
      unsigned int b1jet_btag;
      float b1jet_bDiscVar;
      float b2jet_px;
      float b2jet_py;
      float b2jet_pz;
      float b2jet_eta;
      float b2jet_phi;
      float b2jet_pt;
      float b2jet_energy;
      float b2jet_mass;
      unsigned int b2jet_btag;
      float b2jet_bDiscVar;
      bool hastwojets;
      bool hasb1jet;
      bool hasb2jet;

      float met_pt;
      float met_phi;
      float met_px;
      float met_py;
    
      float minMass;
      float MINdR_bl;
      float dR_bl;
      float dR_b1l1;
      float dR_b1l2;
      float dR_b2l1;
      float dR_b2l2;
      float dR_l1l2;
      float dR_l1l2b1b2;
      float dphi_l1l2b1b2;
      float dR_b1b2;
      float dR_minbl;
      float mass_l1l2;
      float energy_l1l2;
      float pt_l1l2;
      float phi_l1l2;
      float eta_l1l2;
      float mass_b1b2;
      float energy_b1b2;
      float pt_b1b2;
      float phi_b1b2;
      float eta_b1b2;
      float dphi_llbb;
      float dphi_llmet;
      float mass_trans;

     
    private:
      bool runMMC_;
      bool simulation_;
      MMC* thismmc;
     // MMC tree branches
};

    //
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DiHiggsWWBBAnalyzer::DiHiggsWWBBAnalyzer(const edm::ParameterSet& iConfig)
{
  //****************************************************************************
  //                 SET GEN LEVEL VARIABLES AND MATCHING                      
  //****************************************************************************
    genParticlesToken_    = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));
    jetsDeltaR_ = iConfig.getUntrackedParameter<double>("jetsDeltaR",0.4);
    leptonsDeltaR_ = iConfig.getUntrackedParameter<double>("leptonsDeltaR",0.1);

  //****************************************************************************
  //                 SET RECO LEVEL VARIABLES AND COUNTERS                       
  //****************************************************************************
	
    muonToken_           = consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"));
    electronToken_       = consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"));
    jetToken_          = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"));
    genjetToken_          = consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genjets"));
    metToken_          = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"));

    beamSpotToken_        = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"));
    triggerEventToken_    = consumes<pat::TriggerEvent>(iConfig.getParameter<edm::InputTag>("triggerEvent"));
    tracksToken_          = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"));
    trigResToken_         = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"));
    trackRefToken_        = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("TrackRefitter"));
    trajToken_            = consumes< std::vector<Trajectory> >(iConfig.getParameter<edm::InputTag>("Traj"));
    primaryVerticesToken_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"));

  //****************************************************************************
  //                 CLEAR_UP CUTS                       
  //****************************************************************************
    verbose_ = iConfig.getUntrackedParameter<int>("verbose",0);
    mu_eta_ = iConfig.getUntrackedParameter<double>("mu_eta",2.4);
    el_eta_ = iConfig.getUntrackedParameter<double>("el_eta",2.5);
    leadingpt_mumu_ = iConfig.getUntrackedParameter<double>("leadingpt_mumu",10);
    trailingpt_mumu_ = iConfig.getUntrackedParameter<double>("trailingpt_mumu",10);
    leadingpt_muel_ = iConfig.getUntrackedParameter<double>("leadingpt_muel",10);
    trailingpt_muel_ = iConfig.getUntrackedParameter<double>("trailingpt_muel",10);
    leadingpt_elmu_ = iConfig.getUntrackedParameter<double>("leadingpt_elmu",10);
    trailingpt_elmu_ = iConfig.getUntrackedParameter<double>("trailingpt_elmu",10);
    leadingpt_elel_ = iConfig.getUntrackedParameter<double>("leadingpt_elel",10);
    trailingpt_elel_ = iConfig.getUntrackedParameter<double>("trailingpt_elel",10);
    jet_eta_ = iConfig.getUntrackedParameter<double>("jet_eta",2.5);
    jet_leadingpt_ = iConfig.getUntrackedParameter<double>("jet_leadingpt",20);
    jet_trailingpt_ = iConfig.getUntrackedParameter<double>("jet_trailingpt",20);
    bjetDiscrName_ = iConfig.getUntrackedParameter<std::string>("bjetDiscrName","pfCombinedMVAV2BJetTags");
    bjetDiscrCut_loose_ = iConfig.getUntrackedParameter<double>("bjetDiscrCut_loose",0.5);
    bjetDiscrCut_medium_ = iConfig.getUntrackedParameter<double>("bjetDiscrCut_medium",0.7);
    bjetDiscrCut_tight_ = iConfig.getUntrackedParameter<double>("bjetDiscrCut_tight",0.9);
    jetleptonDeltaR_ = iConfig.getUntrackedParameter<double>("jetleptonDeltaR",0.3);
    muIso_ = iConfig.getUntrackedParameter<double>("muIso",0.15);
    elIso_ = iConfig.getUntrackedParameter<double>("elIso",0.04);
    
	//mmcset_ = iConfig.getParameter<edm::ParameterSet>("mmcset"); 
    sampleType_ = iConfig.getUntrackedParameter<int>("SampleType",0);
    finalStates_ = iConfig.getParameter<bool>("finalStates");
    runMMC_ = iConfig.getParameter<bool>("runMMC");
    simulation_ = iConfig.getParameter<bool>("simulation");
/*
     iterations_ = iConfig.getUntrackedParameter<int>("iterations",100000);
     seed_ = iConfig.getParameter<int>("seed");
     RefPDFfile_ = iConfig.getParameter<std::string>("RefPDFfile");
     std::cout <<" RefPDFfile_ " << RefPDFfile_ << std::endl;
     */
    // initilize candidates pointer
    //now do what ever initialization is needed
    ievent = 0;
    mu1_W1_cand = NULL;
    nu1_W1_cand = NULL;
    mu2_W2_cand = NULL;
    nu2_W2_cand = NULL;
    mu1_htoWW_cand = NULL;
    mu2_htoWW_cand = NULL;
    b1_htobb_cand = NULL;
    b2_htobb_cand = NULL;
    h2tohh_cand = NULL;
    //  
}


void 
DiHiggsWWBBAnalyzer::initBranches(){


    findAllGenParticles = false;
    mu1cand = NULL;
    nu1cand = NULL;
    mu2cand = NULL;
    nu2cand = NULL;
    w1cand =NULL;
    w2cand = NULL;
    b1cand = NULL;
    b2cand = NULL;
    htoWWcand = NULL;
    htoBBcand = NULL;
    h2tohhcand = NULL;
    t1cand = NULL;
    t2cand = NULL;
    
    mu1_energy = -1;
    mu1_pt = -1;
    mu1_eta = -9;
    mu1_phi = -9;
    mu1_px = 0.0;
    mu1_py = 0.0;
    mu1_pz = 0.0;
    nu1_energy = -1;
    nu1_pt = -1;
    nu1_eta = -9;
    nu1_phi = -9;
    nu1_px = 0.0;
    nu1_py = 0.0;
    nu1_pz = 0.0;
    w1_energy = -1;
    w1_pt = -1;
    w1_eta = -9;
    w1_phi = -9;
    w1_px = 0.0;
    w1_py = 0.0;
    w1_pz = 0.0;
    w1_mass = 0.0;

    mu2_energy = -1;
    mu2_eta =-9;
    mu2_phi =-9;
    mu2_pt =-1;
    mu2_px = 0.0;
    mu2_py = 0.0;
    mu2_pz = 0.0;
    nu2_energy = -1;
    nu2_eta =-9;
    nu2_phi =-9;
    nu2_pt =-1;
    nu2_px = 0.0;
    nu2_py = 0.0;
    nu2_pz = 0.0;
    w2_energy = -1;
    w2_eta =-9;
    w2_phi =-9;
    w2_pt =-1;
    w2_px = 0.0;
    w2_py = 0.0;
    w2_pz = 0.0;
    w2_mass = 0.0;

    htoWW_energy = -1;
    htoWW_pt = -1;
    htoWW_px = 0.0;
    htoWW_py = 0.0;
    htoWW_pz = 0.0;
    htoWW_mass = 0.0;

    b1_energy = -1;
    b1_pt = -1;
    b1_eta = -9;
    b1_phi = -9;
    b1_px = 0.0;
    b1_py = 0.0;
    b1_pz = 0.0;
    b2_energy = -1;
    b2_pt = -1;
    b2_eta = -9;
    b2_phi = -9;
    b2_px = 0.0;
    b2_py = 0.0;
    b2_pz = 0.0;

    htobb_energy = -1;
    htobb_px = 0.0;
    htobb_py = 0.0;
    htobb_pz = 0.0;
    htobb_mass = 0.0;

    h2tohh_energy = -1;
    h2tohh_px = 0.0;
    h2tohh_py = 0.0;
    h2tohh_pz = 0.0;
    h2tohh_mass = 0.0;

    b1genjet_px=0;
    b1genjet_py=0;
    b1genjet_pz=0;
    b1genjet_eta=-9;
    b1genjet_phi=-9;
    b1genjet_pt=-1;
    b1genjet_energy=-1;
    b2genjet_px=0;
    b2genjet_py=0;
    b2genjet_pz=0;
    b2genjet_eta=-9;
    b2genjet_phi=-9;
    b2genjet_pt=-1;
    b2genjet_energy=-1;
    dR_b1genjet=jetsDeltaR_;
    dR_b2genjet=jetsDeltaR_;
    hastwogenjets = false;

    genmet_pt = -999.;
    genmet_phi = -999.;
    genmet_px = -999.;
    genmet_py = -999.;

    //ttar
    t1_px = 0.0;
    t1_py = 0.0;
    t1_pz = 0.0;
    t1_energy = -1;
    t1_mass = -1;
    t2_px = 0.0;
    t2_py = 0.0;
    t2_pz = 0.0;
    t2_energy = -1;
    t2_mass = -1;

    dR_genbl=-1.0;
    dR_genb1l1=-1.0;
    dR_genb1l2=-1.0;
    dR_genb2l1=-1.0;
    dR_genb2l2=-1.0;
    dR_genl1l2=-1.0;
    dR_genl1l2b1b2=-1.0;
    dphi_genl1l2b1b2=-1.0;
    dR_genb1b2=-1.0;
    dR_genminbl = -1.0;
    mass_genl1l2 = -1.0;
    energy_genl1l2 = -1;
    pt_genl1l2 = -1;
    phi_genl1l2 = -1;
    eta_genl1l2 = -1;
    mass_genb1b2 = -1.0;
    energy_genb1b2 = -1;
    pt_genb1b2 = -1;
    phi_genb1b2 = -9;
    eta_genb1b2 = -9;
    dphi_genllbb = -10;
    dphi_genllmet = -10;
    mass_gentrans = -1;

    //reco level
    numberOfmuon1 = -1;
    numberOfmuon2 = -1;
    muon1_px = 0;
    muon1_py = 0;
    muon1_pz = 0;
    muon1_eta = -9;
    muon1_phi = -9;
    muon1_pt = -1;
    muon1_energy = -1;
    muon1_isoVar = 10.0;
    muon1_dxy = -999;
    muon1_dz = -999;
    muon2_px = 0;
    muon2_py = 0;
    muon2_pz = 0;
    muon2_eta = -9;
    muon2_phi = -9;
    muon2_pt = -1;
    muon2_energy = -1;
    muon2_isoVar = 10.0;
    muon2_dxy = -999;
    muon2_dz = -999;
    dR_mu1 = 2.0;
    dR_mu2 = 2.0;
    hastwomuons = false;

    dR_b1jet = jetsDeltaR_;
    dR_b2jet = jetsDeltaR_;
    b1jet_px=0;
    b1jet_py=0;
    b1jet_pz=0;
    b1jet_eta=-1;
    b1jet_phi=-9;
    b1jet_pt= -1;
    b1jet_energy=-1;
    b1jet_btag = 0;
    b1jet_bDiscVar = 0;
    b2jet_px=0;
    b2jet_py=0;
    b2jet_pz=0;
    b2jet_eta=-9;
    b2jet_phi=-9;
    b2jet_pt=-1;
    b2jet_energy=-1;
    b2jet_btag = 0;
    b2jet_bDiscVar = 0;
    hastwojets = false;
    hasb1jet=false;
    hasb2jet=false;

    met_pt = -999.;
    met_phi = -999.;
    met_px = -999.;
    met_py = -999.;

    minMass=999.;
    dR_bl=-1.0;
    dR_b1l1=-1.0;
    dR_b1l2=-1.0;
    dR_b2l1=-1.0;
    dR_b2l2=-1.0;
    MINdR_bl=-1.0;
    dR_l1l2=-1.0;
    dR_l1l2b1b2=-1.0;
    dphi_l1l2b1b2=-1.0;
    dR_b1b2=-1.0;
    dR_minbl = -1.0;
    dR_genl1l2 = -1;
    mass_l1l2 = -1.0;
    energy_l1l2 = -1;
    pt_l1l2 = -1;
    phi_l1l2 = -1;
    eta_l1l2 = -1;
    mass_b1b2 = -1.0;
    energy_b1b2 = -1;
    pt_b1b2 = -1;
    phi_b1b2 = -1;
    eta_b1b2 = -1;
    dphi_llbb = -10;
    dphi_llmet = -10;
    mass_trans = -1;
}


DiHiggsWWBBAnalyzer::~DiHiggsWWBBAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//
//step 1 find gen particles if it is MC sample
//step 2 select pat objects: matched by gen or not, and apply clear-up cuts at the same time 
//step 3 run MMC on selected objects: gen level or reco level
// ------------ method called for each event  ------------
void
DiHiggsWWBBAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    using namespace edm;
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
//   if(iEvent.isRealData()) std::cout << " Not a real Data " << std::endl;
    initBranches(); 
    ievent++;
    std::cout << "event  " << iEvent.id().event()<<" ievent "<< ievent << std::endl;
    /*
    edm::Handle<pat::METCollection> mets;
    iEvent.getByToken(metToken_, mets);
    const pat::MET &met = mets->front();
    const reco::GenMET *genmet = met.genMET();
    genmet_px = genmet->px(); genmet_py = genmet->py(); genmet_phi = genmet->phi(); genmet_pt = genmet->pt();
    printf("MET: pt %5.1f, phi %+4.2f, sumEt (%.1f). genMET %.1f. MET with JES up/down: %.1f/%.1f\n",
   	 met.pt(), met.phi(), met.sumEt(), met.genMET()->pt(),met.shiftedPt(pat::MET::JetEnUp), met.shiftedPt(pat::MET::JetEnDown));
    met_px = met.px(); met_py = met.py(); met_phi = met.phi(); met_pt = met.pt();
	*/

  //****************************************************************************
  //                GENERATOR LEVEL                       
  //****************************************************************************
    edm::Handle<reco::GenParticleCollection> genParticleColl;
    iEvent.getByToken(genParticlesToken_, genParticleColl);
    edm::Handle<std::vector<reco::GenJet>> genjetColl;
    iEvent.getByToken(genjetToken_, genjetColl);
    if (sampleType_<=12 and sampleType_>0)
	checkGenParticlesSignal(genParticleColl);
    else if (sampleType_==13)
	checkGenParticlesTTbar(genParticleColl);
    if (sampleType_>0 and findAllGenParticles) matchGenJet2Parton( genjetColl );
    if (findAllGenParticles) 
	fillbranches();//fill Gen info into tree

  //****************************************************************************
  //                RECO LEVEL
  //****************************************************************************
  // Cut on primary vertex in event
    edm::Handle<reco::VertexCollection> primaryVertices;
    iEvent.getByToken(primaryVerticesToken_, primaryVertices);
    if (primaryVertices->empty()) return; // skip the event if no PV found
    const reco::Vertex &PV = primaryVertices->front();
   
    edm::Handle<pat::MuonCollection> muons;
    iEvent.getByToken(muonToken_, muons);
    edm::Handle<pat::ElectronCollection> electrons;
    iEvent.getByToken(electronToken_, electrons);
    std::vector<const reco::Candidate *> pleptons;
    std::vector<const reco::Candidate *> nleptons;

    
  //****************************************************************************
  //                Triggering matching 
  //****************************************************************************

  //****************************************************************************
  //                di-Leptons selection
  //****************************************************************************
    //for (pat::MuonCollection::const_iterator iMuon = muons->begin();  iMuon != muons->end();  ++iMuon) {
    //or for (const pat::Muon &mu : *muons) {
    for (const pat::Muon &mu : *muons) {
	const MuonPFIsolation& muonIso = mu.pfIsolationR03();
	float isoVar = (muonIso.sumChargedHadronPt + muonIso.sumNeutralHadronEt + muonIso.sumPhotonEt)/mu.pt();
	if (fabs(mu.eta())<mu_eta_ and mu.pt()>10 and fabs(mu.muonBestTrack()->dz(PV.position()))<0.1 and 
		((mu.pt()>20 and fabs(mu.muonBestTrack()->dxy(PV.position()))<0.02) or 
		(mu.pt()<20 and fabs(mu.muonBestTrack()->dxy(PV.position()))<0.01)) and isoVar<0.15){
	   if (mu.charge()>0) pleptons.push_back(&mu);
	   else if (mu.charge()<0) nleptons.push_back(&mu);
	   std::cout <<"get one muon passed selection eta "<< mu.eta() <<" pt "<< mu.pt()<<" isovar "<< isoVar <<" charge "<< mu.charge()<< std::endl;
	   const reco::GenParticle * genp = mu.genParticle();
	   if (genp)
	       std::cout <<"matched genParticle: id "<< genp->pdgId()<<" px "<< genp->px() <<" py "<< genp->py()<<" pz "<< genp->pz() << std::endl;
	}
	printf("muon with pt %4.1f, charge %d, IsoVar %5.3f, dz(PV) %+5.3f, dxy(PV)%+5.3f, POG loose id %d, tight id %d\n",
		mu.pt(), mu.charge(), isoVar, mu.muonBestTrack()->dz(PV.position()), fabs(mu.muonBestTrack()->dxy(PV.position())),mu.isLooseMuon(), mu.isTightMuon(PV));
    }
    //ignore electron now
    for (const pat::Electron &el : *electrons) {
	if (el.pt()>1000000) continue;
    	
    }
    const reco::Candidate * selectedPlep = NULL;
    const reco::Candidate * selectedNlep = NULL;
    TLorentzVector dilep_p4;
    float sumPt=0.0;
    for (const reco::Candidate *plep : pleptons) {
	for (const reco::Candidate *nlep : nleptons){
	    //select lepton pairs with larger sumPt
	    dilep_p4.SetPxPyPzE(plep->px()+nlep->px(), plep->py()+nlep->py(), plep->pz()+nlep->pz(), plep->energy()+nlep->energy());
	    if(((plep->pt()>10 and nlep->pt()>20) or (nlep->pt()>10 and plep->pt()>20))and dilep_p4.M()>12 and (plep->pt()+nlep->pt())>sumPt){
		selectedPlep = plep;
		selectedNlep = nlep;
		sumPt = plep->pt()+nlep->pt();
	    }
	}
    }
    //bool hastwomuons = false;
    if (sumPt>=30){
	//muon1, mu1: positive charge
	muon1_px = selectedPlep->px(); muon1_py = selectedPlep->py(); muon1_pz = selectedPlep->pz(); muon1_energy = selectedPlep->energy();
	muon1_pt = selectedPlep->pt(); muon1_eta = selectedPlep->eta(); muon1_phi = selectedPlep->phi();
	muon1_dxy = dxy(selectedPlep, &PV); muon1_dz = dz(selectedPlep, &PV);
	muon2_px = selectedNlep->px(); muon2_py = selectedNlep->py(); muon2_pz = selectedNlep->pz(); muon2_energy = selectedNlep->energy();
	muon2_pt = selectedNlep->pt(); muon2_eta = selectedNlep->eta(); muon2_phi = selectedNlep->phi();
	muon2_dxy = dxy(selectedNlep, &PV); muon2_dz = dz(selectedNlep, &PV);
	//if (selectedPlep->genParticle()){
	//    const reco::GenParticle * genp = selectedPlep->genParticle();
	//    std::cout <<"selectedPlep, matched genParticle: id "<< genp->pdgId()<<" px "<< genp->px() <<" py "<< genp->py()<<" pz "<< genp->pz() << std::endl;
	//}
	//if (selectedNlep->genParticle()){
	 //   const reco::GenParticle * genp = selectedPlep->genParticle();
	 //   std::cout <<"selectedPNlep, matched genParticle: id "<< genp->pdgId()<<" px "<< genp->px() <<" py "<< genp->py()<<" pz "<< genp->pz() << std::endl;
	//}
	hastwomuons = true;
    	std::cout <<"Found two leptons " << std::endl;
	printf("plepton: pt %5.1f, eta %+4.2f \n", selectedPlep->pt(), selectedPlep->eta());
	printf("nlepton: pt %5.1f, eta %+4.2f \n", selectedNlep->pt(), selectedNlep->eta());
    }
      
  //****************************************************************************
  //                di-Jets selection
  //****************************************************************************
    edm::Handle<pat::JetCollection> jets;
    iEvent.getByToken(jetToken_, jets);
    std::vector<pat::Jet> allbjets;
    for (const pat::Jet &j : *jets) {
	if (j.pt() < jet_trailingpt_ or fabs(j.eta()) > jet_eta_) continue;
  	if (hastwomuons){
	    TLorentzVector jet_p4(j.px(), j.py(), j.pz(), j.energy());    
	    float dR1 = jet_p4.DeltaR(TLorentzVector(selectedPlep->px(), selectedPlep->py(), selectedPlep->pz(), selectedPlep->energy()));
	    float dR2 = jet_p4.DeltaR(TLorentzVector(selectedNlep->px(), selectedNlep->py(), selectedNlep->pz(), selectedNlep->energy()));
	    if (dR1 < jetleptonDeltaR_ or dR2 < jetleptonDeltaR_) continue;
	}	
	float bDiscVar = j.bDiscriminator(bjetDiscrName_);
	if (bDiscVar < bjetDiscrCut_medium_)  continue;
	allbjets.push_back(j);
	printf("Jet with pt %6.1f, eta %+4.2f, pileup mva disc %+.2f, btag CSV %.3f, CISV %.3f\n",
		j.pt(),j.eta(), j.userFloat("pileupJetId:fullDiscriminant"), std::max(0.f,j.bDiscriminator("combinedSecondaryVertexBJetTags")), std::max(0.f,j.bDiscriminator("combinedInclusiveSecondaryVertexBJetTags")));
	const reco::GenParticle * genp = j.genParticle();
	if (genp)
	    std::cout <<"matched genParticle: id "<< genp->pdgId()<<" px "<< genp->px() <<" py "<< genp->py()<<" pz "<< genp->pz() << std::endl;
    }

    // sort jets by pt
    std::sort(allbjets.begin(), allbjets.end(), [](pat::Jet& jet1, pat::Jet& jet2) { return jet1.pt() > jet2.pt(); });
    std::cout <<"allbjets size "<< allbjets.size() << std::endl;
    unsigned int jet1=0, jet2=0;
    float diff_higgsmass = 9999;
    for (unsigned int i=0; i<allbjets.size(); i++){
    	for (unsigned int j=i+1; j<allbjets.size(); j++){
	    TLorentzVector dijet_p4(allbjets[i].px()+allbjets[j].px(), allbjets[i].py()+allbjets[j].py(), 
		    allbjets[i].pz()+allbjets[j].pz(),allbjets[i].energy()+allbjets[j].energy());
	    if (fabs(dijet_p4.M()-125)<diff_higgsmass){
	       jet1 = i;
	       jet2 = j;
	       diff_higgsmass = fabs(dijet_p4.M()-125); 
	    }
	}
    }
    if (allbjets.size()>2){
	std::cout <<"found two bjets "<< std::endl;
	b1jet_px = allbjets[jet1].px(); b1jet_py = allbjets[jet1].py(); b1jet_pz = allbjets[jet1].pz(); b1jet_energy = allbjets[jet1].energy();
	b1jet_pt = allbjets[jet1].pt(); b1jet_eta = allbjets[jet1].eta(); b1jet_phi = allbjets[jet1].phi();
	b1jet_bDiscVar = allbjets[jet1].bDiscriminator(bjetDiscrName_);
	b2jet_px = allbjets[jet2].px(); b2jet_py = allbjets[jet1].py(); b2jet_pz = allbjets[jet2].pz(); b2jet_energy = allbjets[jet2].energy();
	b2jet_pt = allbjets[jet2].pt(); b2jet_eta = allbjets[jet1].eta(); b2jet_phi = allbjets[jet2].phi();
	b2jet_bDiscVar = allbjets[jet2].bDiscriminator(bjetDiscrName_);
	hastwojets = true;
    }

   if (hastwomuons and hastwojets){
          TLorentzVector Muon1_p4(muon1_px, muon1_py, muon1_pz, muon1_energy); 
          TLorentzVector Muon2_p4(muon2_px, muon2_py, muon2_pz, muon2_energy); 
          TLorentzVector b1jet_p4(b1jet_px, b1jet_py, b1jet_pz, b1jet_energy); 
          TLorentzVector b2jet_p4(b2jet_px, b2jet_py, b2jet_pz, b2jet_energy); 
	  dR_bl   = (b1jet_p4.Pt()>b2jet_p4.Pt()) ? (b1jet_p4.DeltaR( (Muon1_p4.Pt()>Muon2_p4.Pt()) ? Muon1_p4 : Muon2_p4 )) : (b2jet_p4.DeltaR( (Muon1_p4.Pt()>Muon2_p4.Pt()) ? Muon1_p4 : Muon2_p4 ));
	  dR_b1l1 = b1jet_p4.DeltaR(Muon1_p4);
	  dR_b1l2 = b1jet_p4.DeltaR(Muon2_p4);
	  dR_b2l1 = b2jet_p4.DeltaR(Muon1_p4);
	  dR_b2l2 = b2jet_p4.DeltaR(Muon2_p4);
	  dR_b1b2 = b1jet_p4.DeltaR(b2jet_p4);
	  dR_l1l2 = Muon1_p4.DeltaR(Muon2_p4);
	  dR_l1l2b1b2 = (Muon1_p4+Muon2_p4).DeltaR(b1jet_p4+b2jet_p4);
	  dphi_l1l2b1b2 = fabs(TVector2::Phi_mpi_pi( (Muon1_p4+Muon2_p4).Phi()-(b1jet_p4+b2jet_p4).Phi() ));
	  TLorentzVector ll_p4 = Muon1_p4+Muon2_p4;
	  TLorentzVector bjets_p4 = b1jet_p4+b2jet_p4;
	  dR_minbl = std::min(std::min(dR_b1l1,dR_b1l2),std::min(dR_b2l1,dR_b2l2));
	  dphi_llbb = TVector2::Phi_mpi_pi(ll_p4.Phi()-bjets_p4.Phi());
	  dphi_llmet = TVector2::Phi_mpi_pi(ll_p4.Phi()-met_phi);
	  mass_l1l2 = ll_p4.M(); energy_l1l2 = ll_p4.Energy(); pt_l1l2 = ll_p4.Pt(); eta_l1l2 = ll_p4.Eta(); phi_l1l2 = ll_p4.Phi();
	  mass_b1b2 = bjets_p4.M(); energy_b1b2 = bjets_p4.Energy(); pt_b1b2 = bjets_p4.Pt(); eta_b1b2 = bjets_p4.Eta(); phi_b1b2 = bjets_p4.Phi();
	  mass_trans = sqrt(2*ll_p4.Pt()*met_pt*(1-cos(dphi_llmet)));
	  //if (dR_b1l1 > jetleptonDeltaR_ and dR_b1l2 > jetleptonDeltaR_ and dR_b2l1 > jetleptonDeltaR_ and dR_b2l2 > jetleptonDeltaR_) hasdRljet =true;
	  //MINdR_bl = dR_b1l1*(dR_b1l1<dR_b1l2 && dR_b1l1<dR_b2l1 && dR_b1l1<dR_b2l2) + dR_b2l1*(dR_b2l1<dR_b2l2 && dR_b2l1<dR_b1l1 && dR_b2l1<dR_b1l2) + dR_b1l2*(dR_b1l2<dR_b1l1 && dR_b1l2<dR_b2l1 && dR_b1l2<dR_b2l2) + dR_b2l2*(dR_b2l2<dR_b1l1 && dR_b2l2<dR_b1l2 && dR_b2l2<dR_b2l1);
     
	  //MT2: In order to construct MT2 for either the t tbar -> bW bW system or our signal H -> h h -> bb WW,
	  //we group each pair b_jet-lepton into an object: we then get "Particle A" and "Particle B" (each one given by a b_jet-lepton object),
	  //whose Kinematics is to be fed into the MT2 variable.
	  //There are two possible Lepton-Bquark pairings. We compute MT2 for both and pick the smallest value.
   
   }*/

   if (findAllGenParticles) evtree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
DiHiggsWWBBAnalyzer::beginJob()
{
    evtree = fs->make<TTree>("evtree", "evtree");
 //  output = new TFile("output.root","recreate");
      // output->cd();
    evtree->Branch("ievent",&ievent);
 //  evtree = new TTree("evtree","event tree");
 
   //evtree->Branch("h2tohh",&findAllGenParticles);
    evtree->Branch("mu1_px",&mu1_px, "mu1_px/F");
    evtree->Branch("mu1_py",&mu1_py, "mu1_py/F");
    evtree->Branch("mu1_pz",&mu1_pz, "mu1_pz/F");
    evtree->Branch("mu1_eta",&mu1_eta, "mu1_eta/F");
    evtree->Branch("mu1_phi",&mu1_phi, "mu1_phi/F");
    evtree->Branch("mu1_pt",&mu1_pt, "mu1_pt/F");
    evtree->Branch("mu1_energy",&mu1_energy, "mu1_energy/F");
    evtree->Branch("nu1_px",&nu1_px, "nu1_px/F");
    evtree->Branch("nu1_py",&nu1_py, "nu1_py/F");
    evtree->Branch("nu1_pz",&nu1_pz, "nu1_pz/F");
    evtree->Branch("nu1_eta",&nu1_eta, "nu1_eta/F");
    evtree->Branch("nu1_phi",&nu1_phi, "nu1_phi/F");
    evtree->Branch("nu1_pt",&nu1_pt, "nu1_pt/F");
    evtree->Branch("nu1_energy",&nu1_energy, "nu1_energy/F");
    evtree->Branch("mu2_px",&mu2_px, "mu2_px/F");
    evtree->Branch("mu2_py",&mu2_py, "mu2_py/F");
    evtree->Branch("mu2_pz",&mu2_pz, "mu2_pz/F");
    evtree->Branch("mu2_eta",&mu2_eta, "mu2_eta/F");
    evtree->Branch("mu2_phi",&mu2_phi, "mu2_phi/F");
    evtree->Branch("mu2_pt",&mu2_pt, "mu2_pt/F");
    evtree->Branch("mu2_energy",&mu2_energy, "mu2_energy/F");
    evtree->Branch("nu2_px",&nu2_px, "nu2_px/F");
    evtree->Branch("nu2_py",&nu2_py, "nu2_py/F");
    evtree->Branch("nu2_pz",&nu2_pz, "nu2_pz/F");
    evtree->Branch("nu2_eta",&nu2_eta, "nu2_eta/F");
    evtree->Branch("nu2_phi",&nu2_phi, "nu2_phi/F");
    evtree->Branch("nu2_pt",&nu2_pt, "nu2_pt/F");
    evtree->Branch("nu2_energy",&nu2_energy, "nu2_energy/F");
    //evtree->Branch("nu1and2_pt",&nu1and2_pt, "nu1and2_pt/F");
    //evtree->Branch("nu1and2_px",&nu1and2_px, "nu1and2_px/F");
    //evtree->Branch("nu1and2_py",&nu1and2_py, "nu1and2_py/F");
    //evtree->Branch("nu1and2_diBaxis_p",&nu1and2_diBaxis_p, "nu1and2_diBaxis_p/F");
    //evtree->Branch("nu1and2_diBaxis_t",&nu1and2_diBaxis_t, "nu1and2_diBaxis_t/F");

    evtree->Branch("w1_mass",&w1_mass, "w1_mass/F");
    evtree->Branch("w1_px",&w1_px, "w1_px/F");
    evtree->Branch("w1_py",&w1_py, "w1_py/F");
    evtree->Branch("w1_pz",&w1_pz, "w1_pz/F");
    evtree->Branch("w1_energy",&w1_energy, "w1_energy/F");
    evtree->Branch("w1_pt",&w1_pt, "w1_pt/F");
    evtree->Branch("w1_eta",&w1_eta, "w1_eta/F");
    evtree->Branch("w1_phi",&w1_phi, "w1_phi/F");
    evtree->Branch("w2_mass",&w2_mass, "w2_mass/F");
    evtree->Branch("w2_px",&w2_px, "w2_px/F");
    evtree->Branch("w2_py",&w2_py, "w2_py/F");
    evtree->Branch("w2_pz",&w2_pz, "w2_pz/F");
    evtree->Branch("w2_energy",&w2_energy, "w2_energy/F");
    evtree->Branch("w2_pt",&w2_pt, "w2_pt/F");
    evtree->Branch("w2_eta",&w2_eta, "w2_eta/F");
    evtree->Branch("w2_phi",&w2_phi, "w2_phi/F");
   
    evtree->Branch("htoWW_energy",&htoWW_energy);
    evtree->Branch("htoWW_px",&htoWW_px,"htoWW_px/F");
    evtree->Branch("htoWW_py",&htoWW_py,"htoWW_px/F");
    evtree->Branch("htoWW_pz",&htoWW_pz,"htoWW_pz/F");
    evtree->Branch("htoWW_mass",&htoWW_mass,"htoWW_mass/F");
    //evtree->Branch("Wtomu1nu1",&Wtomu1nu1,"Wtomu1nu1/B");
    //evtree->Branch("Wtomu2nu2",&Wtomu2nu2,"Wtomu2nu2/B");
    //evtree->Branch("htoWW",&htoWW,"htoWW/B");

    evtree->Branch("b1_px",&b1_px, "b1_px/F");
    evtree->Branch("b1_py",&b1_py, "b1_py/F");
    evtree->Branch("b1_pz",&b1_pz, "b1_pz/F");
    evtree->Branch("b1_eta",&b1_eta, "b1_eta/F");
    evtree->Branch("b1_phi",&b1_phi, "b1_phi/F");
    evtree->Branch("b1_pt",&b1_pt, "b1_pt/F");
    evtree->Branch("b1_energy",&b1_energy, "b1_energy/F");
    evtree->Branch("b2_px",&b2_px, "b2_px/F");
    evtree->Branch("b2_py",&b2_py, "b2_py/F");
    evtree->Branch("b2_pz",&b2_pz, "b2_pz/F");
    evtree->Branch("b2_eta",&b2_eta, "b2_eta/F");
    evtree->Branch("b2_phi",&b2_phi, "b2_phi/F");
    evtree->Branch("b2_pt",&b2_pt, "b2_pt/F");
    evtree->Branch("b2_energy",&b2_energy, "b2_energy/F");
    evtree->Branch("htobb_px",&htobb_px, "htobb_px/F");
    evtree->Branch("htobb_py",&htobb_py, "htobb_py/F");
    evtree->Branch("htobb_pz",&htobb_pz, "htobb_pz/F");
    evtree->Branch("htobb_energy",&htobb_energy, "htobb_energy/F");
    evtree->Branch("htobb_mass",&htobb_mass, "htobb_mass/F");

    evtree->Branch("b1genjet_px",&b1genjet_px, "b1genjet_px/F");
    evtree->Branch("b1genjet_py",&b1genjet_py, "b1genjet_py/F");
    evtree->Branch("b1genjet_pz",&b1genjet_pz, "b1genjet_pz/F");
    evtree->Branch("b1genjet_eta",&b1genjet_eta, "b1genjet_eta/F");
    evtree->Branch("b1genjet_phi",&b1genjet_phi, "b1genjet_phi/F");
    evtree->Branch("b1genjet_pt",&b1genjet_pt, "b1genjet_pt/F");
    evtree->Branch("b1genjet_energy",&b1genjet_energy, "b1genjet_energy/F");
    evtree->Branch("b1genjet_mass",&b1genjet_mass, "b1genjet_mass/F");
    evtree->Branch("b2genjet_px",&b2genjet_px, "b2genjet_px/F");
    evtree->Branch("b2genjet_py",&b2genjet_py, "b2genjet_py/F");
    evtree->Branch("b2genjet_pz",&b2genjet_pz, "b2genjet_pz/F");
    evtree->Branch("b2genjet_eta",&b2genjet_eta, "b2genjet_eta/F");
    evtree->Branch("b2genjet_phi",&b2genjet_phi, "b2genjet_phi/F");
    evtree->Branch("b2genjet_pt",&b2genjet_pt, "b2genjet_pt/F");
    evtree->Branch("b2genjet_energy",&b2genjet_energy, "b2genjet_energy/F");
    evtree->Branch("b2genjet_mass",&b2genjet_mass, "b2genjet_mass/F");
    evtree->Branch("dR_b1genjet", &dR_b1genjet,"dR_b1genjet/F");  
    evtree->Branch("dR_b2genjet", &dR_b2genjet,"dR_b2genjet/F");  
    evtree->Branch("hastwogenjets", &hastwogenjets,"hastwogenjets/B");  

    evtree->Branch("genmet_pt",&genmet_pt,"genmet_pt/F");
    evtree->Branch("genmet_phi",&genmet_phi,"genmet_phi/F");
    evtree->Branch("genmet_px",&genmet_px,"genmet_px/F");
    evtree->Branch("genmet_py",&genmet_py,"genmet_py/F");

    evtree->Branch("t1_px",&t1_px,"t1_px/F");
    evtree->Branch("t1_py",&t1_py,"t1_py/F");
    evtree->Branch("t1_pz",&t1_pz,"t1_pz/F");
    evtree->Branch("t1_energy",&t1_energy,"t1_energy/F");
    evtree->Branch("t1_mass",&t1_mass,"t1_mass/F");
    evtree->Branch("t2_px",&t2_px,"t2_px/F");
    evtree->Branch("t2_py",&t2_py,"t2_py/F");
    evtree->Branch("t2_pz",&t2_pz,"t2_pz/F");
    evtree->Branch("t2_energy",&t2_energy,"t2_energy/F");
    evtree->Branch("t2_mass",&t2_mass,"t2_mass/F");

    evtree->Branch("dR_genbl",&dR_genbl, "dR_genbl/F");
    evtree->Branch("dR_genb1l1",&dR_genb1l1, "dR_genb1l1/F");
    evtree->Branch("dR_genb1l2",&dR_genb1l2, "dR_genb1l2/F");
    evtree->Branch("dR_genb2l1",&dR_genb2l1, "dR_genb2l1/F");
    evtree->Branch("dR_genb2l2",&dR_genb2l2, "dR_genb2l2/F");
    evtree->Branch("dR_genb1b2",&dR_genb1b2, "dR_genb1b2/F");
    evtree->Branch("dR_genl1l2",&dR_genl1l2, "dR_genl1l2/F");
    evtree->Branch("dR_genl1l2b1b2",&dR_genl1l2b1b2, "dR_genl1l2b1b2/F");
    evtree->Branch("dphi_genl1l2b1b2",&dphi_genl1l2b1b2, "dphi_genl1l2b1b2/F");
    evtree->Branch("dR_genminbl",&dR_genminbl, "dR_genminbl/F");
    evtree->Branch("mass_genb1b2",&mass_genb1b2, "mass_genb1b2/F");
    evtree->Branch("energy_genb1b2",&energy_genb1b2, "energy_genb1b2/F");
    evtree->Branch("pt_genb1b2",&pt_genb1b2, "pt_genb1b2/F");
    evtree->Branch("phi_genb1b2",&phi_genb1b2, "phi_genb1b2/F");
    evtree->Branch("eta_genb1b2",&eta_genb1b2, "eta_genb1b2/F");
    evtree->Branch("mass_genl1l2",&mass_genl1l2, "mass_genl1l2/F");
    evtree->Branch("energy_genl1l2",&energy_genl1l2, "energy_genl1l2/F");
    evtree->Branch("pt_genl1l2",&pt_genl1l2, "pt_genl1l2/F");
    evtree->Branch("phi_genl1l2",&phi_genl1l2, "phi_genl1l2/F");
    evtree->Branch("eta_genl1l2",&eta_genl1l2, "eta_genl1l2/F");
    evtree->Branch("dphi_genllbb",&dphi_genllbb, "dphi_genllbb/F");
    evtree->Branch("dphi_genllmet",&dphi_genllmet, "dphi_genllmet/F");
    evtree->Branch("mass_gentrans",&mass_gentrans, "mass_gentrans/F");


    //reco level
    evtree->Branch("muon1_px",&muon1_px, "muon1_px/F");
    evtree->Branch("muon1_py",&muon1_py, "muon1_py/F");
    evtree->Branch("muon1_pz",&muon1_pz, "muon1_pz/F");
    evtree->Branch("muon1_eta",&muon1_eta, "muon1_eta/F");
    evtree->Branch("muon1_phi",&muon1_phi, "muon1_phi/F");
    evtree->Branch("muon1_pt",&muon1_pt, "muon1_pt/F");
    evtree->Branch("muon1_energy",&muon1_energy, "muon1_energy/F");
    evtree->Branch("muon1_isoVar",&muon1_isoVar, "muon1_isoVar/F");
    evtree->Branch("muon1_dxy",&muon1_dxy, "muon1_dxy/F");
    evtree->Branch("muon1_dz",&muon1_dz, "muon1_dz/F");
    evtree->Branch("muon2_px",&muon2_px, "muon2_px/F");
    evtree->Branch("muon2_py",&muon2_py, "muon2_py/F");
    evtree->Branch("muon2_pz",&muon2_pz, "muon2_pz/F");
    evtree->Branch("muon2_eta",&muon2_eta, "muon2_eta/F");
    evtree->Branch("muon2_phi",&muon2_phi, "muon2_phi/F");
    evtree->Branch("muon2_pt",&muon2_pt, "muon2_pt/F");
    evtree->Branch("muon2_energy",&muon2_energy, "muon2_energy/F");
    evtree->Branch("muon2_dxy",&muon1_dxy, "muon2_dxy/F");
    evtree->Branch("muon2_dz",&muon1_dz, "muon2_dz/F");
    evtree->Branch("muon2_isoVar",&muon2_isoVar, "muon2_isoVar/F");
    evtree->Branch("dR_mu1",&dR_mu1, "dR_mu1/F");
    evtree->Branch("dR_mu2",&dR_mu2, "dR_mu2/F");
    evtree->Branch("hastwomuons",&hastwomuons, "hastwomuons/B");
    //evtree->Branch("hasmuon1",&hasmuon1, "hasmuon1/B");
    //evtree->Branch("hasmuon2",&hasmuon2, "hasmuon2/B");
    //evtree->Branch("hasRecomuon1",&hasRecomuon1, "hasRecomuon1/B");
    //evtree->Branch("hasRecomuon2",&hasRecomuon2, "hasRecomuon2/B");
    //evtree->Branch("muon1_hasgenMu",&muon1_hasgenMu, "muon1_hasgenMu/B");
    //evtree->Branch("muon2_hasgenMu",&muon2_hasgenMu, "muon2_hasgenMu/B");

    evtree->Branch("b1jet_px",&b1jet_px, "b1jet_px/F");
    evtree->Branch("b1jet_py",&b1jet_py, "b1jet_py/F");
    evtree->Branch("b1jet_pz",&b1jet_pz, "b1jet_pz/F");
    evtree->Branch("b1jet_eta",&b1jet_eta, "b1jet_eta/F");
    evtree->Branch("b1jet_phi",&b1jet_phi, "b1jet_phi/F");
    evtree->Branch("b1jet_pt",&b1jet_pt, "b1jet_pt/F");
    evtree->Branch("b1jet_energy",&b1jet_energy, "b1jet_energy/F");
    evtree->Branch("b1jet_mass",&b1jet_mass, "b1jet_mass/F");
    evtree->Branch("b1jet_btag",&b1jet_btag, "b1jet_btag/i");//unsigned int
    evtree->Branch("b1jet_bDiscVar",&b1jet_bDiscVar, "b1jet_bDiscVar/F");
    evtree->Branch("b2jet_px",&b2jet_px, "b2jet_px/F");
    evtree->Branch("b2jet_py",&b2jet_py, "b2jet_py/F");
    evtree->Branch("b2jet_pz",&b2jet_pz, "b2jet_pz/F");
    evtree->Branch("b2jet_eta",&b2jet_eta, "b2jet_eta/F");
    evtree->Branch("b2jet_phi",&b2jet_phi, "b2jet_phi/F");
    evtree->Branch("b2jet_pt",&b2jet_pt, "b2jet_pt/F");
    evtree->Branch("b2jet_energy",&b2jet_energy, "b2jet_energy/F");
    evtree->Branch("b2jet_mass",&b2jet_mass, "b2jet_mass/F");
    evtree->Branch("b2jet_btag",&b2jet_btag, "b2jet_btag/i");
    evtree->Branch("b2jet_bDiscVar",&b2jet_bDiscVar, "b2jet_bDiscVar/F");
    evtree->Branch("dR_b1jet", &dR_b1jet,"dR_b1jet/F");  
    evtree->Branch("dR_b2jet", &dR_b2jet,"dR_b2jet/F");  
    evtree->Branch("hastwojets",&hastwojets, "hastwojets/B");
    evtree->Branch("hasb1jet",&hasb1jet, "hasb1jet/B");
    evtree->Branch("hasb2jet",&hasb2jet, "hasb2jet/B");

    evtree->Branch("met_pt",&met_pt,"met_pt/F");
    evtree->Branch("met_phi",&met_phi,"met_phi/F");
    evtree->Branch("met_px",&met_px,"met_px/F");
    evtree->Branch("met_py",&met_py,"met_py/F");

    evtree->Branch("dR_bl",&dR_bl, "dR_bl/F");
    evtree->Branch("dR_b1l1",&dR_b1l1, "dR_b1l1/F");
    evtree->Branch("dR_b1l2",&dR_b1l2, "dR_b1l2/F");
    evtree->Branch("dR_b2l1",&dR_b2l1, "dR_b2l1/F");
    evtree->Branch("dR_b2l2",&dR_b2l2, "dR_b2l2/F");
    evtree->Branch("dR_b1b2",&dR_b1b2, "dR_b1b2/F");
    evtree->Branch("MINdR_bl",&MINdR_bl, "MINdR_bl/F");
    evtree->Branch("dR_l1l2",&dR_l1l2, "dR_l1l2/F");
    evtree->Branch("dR_l1l2b1b2",&dR_l1l2b1b2, "dR_l1l2b1b2/F");
    evtree->Branch("dphi_l1l2b1b2",&dphi_l1l2b1b2, "dphi_l1l2b1b2/F");
    evtree->Branch("dR_minbl",&dR_minbl, "dR_minbl/F");
    evtree->Branch("dR_genl1l2",&dR_genl1l2, "dR_genl1l2/F");
    evtree->Branch("mass_b1b2",&mass_b1b2, "mass_b1b2/F");
    evtree->Branch("energy_b1b2",&energy_b1b2, "energy_b1b2/F");
    evtree->Branch("pt_b1b2",&pt_b1b2, "pt_b1b2/F");
    evtree->Branch("phi_b1b2",&phi_b1b2, "phi_b1b2/F");
    evtree->Branch("eta_b1b2",&eta_b1b2, "eta_b1b2/F");
    evtree->Branch("mass_l1l2",&mass_l1l2, "mass_l1l2/F");
    evtree->Branch("energy_l1l2",&energy_l1l2, "energy_l1l2/F");
    evtree->Branch("pt_l1l2",&pt_l1l2, "pt_l1l2/F");
    evtree->Branch("phi_l1l2",&phi_l1l2, "phi_l1l2/F");
    evtree->Branch("eta_l1l2",&eta_l1l2, "eta_l1l2/F");
    evtree->Branch("dphi_llbb",&dphi_llbb, "dphi_llbb/F");
    evtree->Branch("dphi_llmet",&dphi_llmet, "dphi_llmet/F");
    evtree->Branch("mass_trans",&mass_trans, "mass_trans/F");

    
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DiHiggsWWBBAnalyzer::endJob() 
{
//   std::cout << "endJob, ievent  " << ievent << std::endl;
    //output->Write();
   // output->Close();
/*
   // release space 
   delete mu_onshellW_lorentz;
   delete mu_offshellW_lorentz;
   delete nu_onshellW_lorentz;
   delete nu_offshellW_lorentz;
   delete onshellW_lorentz;
   delete offshellW_lorentz;
   delete htoWW_lorentz;
   delete htoBB_lorentz;
   delete h2tohh_lorentz;*/
   //delete jets_lorentz;
}


// ------------ method called to check singal genParticles   ------------
void 
DiHiggsWWBBAnalyzer::checkGenParticlesSignal(edm::Handle<reco::GenParticleCollection> genParticleColl){

    std::vector<reco::GenParticle*> b1Coll; 
    std::vector<reco::GenParticle*> b2Coll;
    std::vector<reco::GenParticle*> W1Coll;
    std::vector<reco::GenParticle*> W2Coll;
    std::vector<const reco::Candidate*> htoWWColl;
    std::vector<const reco::Candidate*> htoBBColl;

    std::cout <<"*********** start to check GenParticles for Signal sample ***********"<< std::endl;
    bool h2tohh = false;
    for (reco::GenParticleCollection::const_iterator it = genParticleColl->begin(); it != genParticleColl->end(); ++it) {
    //particle id, (electron12)(muon13),(b5),(W+24),(SM higgs25)
    // particle id  it->pdgId()
    //std::cout << "Gen paticles: id " << it->pdgId() << std::endl; 
	//if (it->pdgId() == 24 && hasDaughter(it->clone(), -13) ){
      if (it->pdgId() == 24){
        const reco::Candidate* tmpw1 = it->mother();
        while (tmpw1->pdgId() == 24 && tmpw1->numberOfMothers() == 1) tmpw1 = tmpw1->mother();
        if (tmpw1->pdgId() == 25)  W1Coll.push_back(it->clone());
	  }
      else if (it->pdgId() == -24){
        const reco::Candidate* tmpw2 = it->mother();
        while (tmpw2->pdgId() == -24 && tmpw2->numberOfMothers() == 1) tmpw2 = tmpw2->mother();
        if (tmpw2->pdgId() == 25)  W2Coll.push_back(it->clone());
      }
      else if (it->pdgId() == 5 && it->mother()->pdgId() == 25 )
      {
        if (it->numberOfMothers() != 1) std::cout << "bquark has more than one mother particle" << std::endl;
        //   std::cout << "find bquark candidate" << std::endl;
        b1Coll.push_back(it->clone());
      }
      else if (it->pdgId() == -5 && it->mother()->pdgId() == 25 )
      {
        if (it->numberOfMothers() != 1) std::cout << "bquark has more than one mother particle" << std::endl;
        // std::cout << "find bbarquark candidate" << std::endl;
        b2Coll.push_back(it->clone());
      }

    }// all Gen particles

    std::cout <<"size W1Coll "<< W1Coll.size() <<" W2Coll "<< W2Coll.size()<<" b1Coll "<< b1Coll.size()<<" b2Coll "<< b2Coll.size()<< std::endl;
    //htoWW
    if (W1Coll.size() && W2Coll.size()){
	for (auto W1_cand : W1Coll)
            for (auto W2_cand : W2Coll){
                const reco::Candidate* W1_mother = W1_cand->mother();
                const reco::Candidate* W2_mother = W2_cand->mother();
                while (W1_mother->pdgId() == -24) W1_mother = W1_mother->mother();
                while (W2_mother->pdgId() ==  24) W2_mother = W2_mother->mother();
                if (W1_mother == W2_mother && W1_mother->pdgId() == 25) {
		    htoWWColl.push_back(W1_mother);
		    break;
		}
	    }
    }
    //htoBB
    if (b1Coll.size() && b2Coll.size()){
	for(auto b1_cand : b1Coll)
	    for (auto b2_cand : b2Coll) {
                const reco::Candidate* b1_mother = b1_cand->mother();
                const reco::Candidate* b2_mother = b2_cand->mother();
                if (b1_mother == b2_mother && b1_mother->pdgId() == 25) {
		    htoBBColl.push_back(b1_mother);
		    break;
		}
	    }
    }
    //h2tohh
    if (htoWWColl.size() && htoBBColl.size()){
	for (auto htoWW_cand : htoWWColl){
            for (auto htoBB_cand : htoBBColl){
                const reco::Candidate* htoWW_mother = htoWW_cand->mother();
                const reco::Candidate* htoBB_mother = htoBB_cand->mother();
                while (htoWW_mother->pdgId() == 25)  htoWW_mother = htoWW_mother->mother();
                while (htoBB_mother->pdgId() == 25)  htoBB_mother = htoBB_mother->mother();
                if (htoWW_mother == htoBB_mother && (htoBB_mother->pdgId()==99927 || htoBB_mother->pdgId()==99926)){ 
		    h2tohhcand = htoWW_mother;
		    htoWWcand = htoWW_cand;
		    htoBBcand = htoBB_cand;
		    h2tohh = true;
		    break;
		}
            }
	    if (h2tohh) break;
		
	}
	if (h2tohh)
	    std::cout << "find h2 candidate " << std::endl;
	else std::cout <<"failed to find h2tohh "<< std::endl;
    }

    if (h2tohh){
        std::cout << "h2 candidate id " << h2tohhcand->pdgId() << " mass " << h2tohhcand->mass() << std::endl;
	b1cand = finddecendant(htoBBcand, 5, false);
	b2cand = finddecendant(htoBBcand, -5, false);
	w1cand = finddecendant(htoWWcand, 24, false);
	w2cand = finddecendant(htoWWcand, -24, false);   
	if (hasDaughter(w1cand, -13) and hasDaughter(w2cand, 13)){
	    mu1cand = findmudaughter(w1cand);
	    nu1cand = findnudaughter(w1cand);
	    mu2cand = findmudaughter(w2cand);
	    nu2cand = findnudaughter(w2cand);
	    //make sure all candiates are in same frame
	    while (mu1cand->numberOfDaughters()==1 and  mu1cand->daughter(0)->pdgId()==mu1cand->pdgId())
		    mu1cand = mu1cand->daughter(0);
	    while (nu1cand->numberOfDaughters()==1 and  nu1cand->daughter(0)->pdgId()==nu1cand->pdgId())
		    nu1cand = nu1cand->daughter(0);
	    while (mu2cand->numberOfDaughters()==1 and  mu2cand->daughter(0)->pdgId()==mu2cand->pdgId())
		    mu2cand = mu2cand->daughter(0);
	    while (nu2cand->numberOfDaughters()==1 and  nu2cand->daughter(0)->pdgId()==nu2cand->pdgId())
		    nu2cand = nu2cand->daughter(0);
	    findAllGenParticles = true;
	    std::cout <<" mu1 " ; printCandidate(mu1cand);
	    std::cout <<" nu1 " ; printCandidate(nu1cand);
	    std::cout <<" mu2 " ; printCandidate(mu2cand);
	    std::cout <<" nu2 " ; printCandidate(nu2cand);
	}else{
	    findAllGenParticles = false;
	    std::cout <<"failed to two leptons from  W decays "<< std::endl;
	}
        std::cout <<" w1 " ; printCandidate(w1cand);
        std::cout <<" w2 " ; printCandidate(w2cand);
        std::cout <<" b1 " ; printCandidate(b1cand);
        std::cout <<" b2 " ; printCandidate(b2cand);
    }

}



// ------------ method called to check ttbar genParticles   ------------
void 
DiHiggsWWBBAnalyzer::checkGenParticlesTTbar(edm::Handle<reco::GenParticleCollection> genParticleColl){

    std::vector<reco::GenParticle*> b1Coll; 
    std::vector<reco::GenParticle*> b2Coll;
    std::vector<reco::GenParticle*> W1Coll;
    std::vector<reco::GenParticle*> W2Coll;
    std::vector<const reco::Candidate*> tColl;
    std::vector<const reco::Candidate*> tbarColl;

    std::cout <<"*********** start to check GenParticles for TTbar sample ***********"<< std::endl;
    for (reco::GenParticleCollection::const_iterator it = genParticleColl->begin(); it != genParticleColl->end(); ++it) {
    //particle id, (electron12)(muon13),(b5),(W+24),(SM higgs25)
    // particle id  it->pdgId()
    //std::cout << "Gen paticles: id " << it->pdgId() << std::endl; 
	//if (it->pdgId() == 24 && hasDaughter(it->clone(), -13) ){
      if (it->pdgId() == 24 ){
	const reco::Candidate* tmpw1 = it->mother();
	//while (tmpw1->pdgId() == 24 && tmpw1->numberOfMothers() == 1) tmpw1 = tmpw1->mother();
	if (tmpw1->pdgId() == 6)  W1Coll.push_back(it->clone());
      }
      else if (it->pdgId() == -24 ){
        const reco::Candidate* tmpw2 = it->mother();
        //while (tmpw2->pdgId() == -24 && tmpw2->numberOfMothers() == 1) tmpw2 = tmpw2->mother();
        if (tmpw2->pdgId() == -6)  W2Coll.push_back(it->clone());
      }
      else if (it->pdgId() == 5 && it->mother()->pdgId() == 6 ){
          if (it->numberOfMothers() != 1) std::cout << "bquark has more than one mother particle" << std::endl;
          b1Coll.push_back(it->clone());
      }
      else if (it->pdgId() == -5 && it->mother()->pdgId() == -6 ){
	  if (it->numberOfMothers() != 1) std::cout << "bquark has more than one mother particle" << std::endl;
	  b2Coll.push_back(it->clone());
      }
    }// all Gen particles


    std::cout <<"size W1Coll "<< W1Coll.size() <<" W2Coll "<< W2Coll.size()<<" b1Coll "<< b1Coll.size()<<" b2Coll "<< b2Coll.size()<< std::endl;
    //t->bW+ 6->5,24
    if (W1Coll.size() && b1Coll.size()){
         for (auto W1_cand : W1Coll)
              for (auto b1_cand : b1Coll){
		  const reco::Candidate* b1_mother = b1_cand->mother();
		  const reco::Candidate* W1_mother = W1_cand->mother();
	  //        while (W1_mother->pdgId() == 24) W1_mother = W1_mother->mother();
		  if (W1_mother == b1_mother && W1_mother->pdgId() == 6) {
			tColl.push_back(W1_mother);
			//std::cout <<" find t->bW+" << std::endl;
			break;
		    }
              }
    }

     //tbar->bbarW- -6->-5,-24
    if (W2Coll.size() && b2Coll.size()){
         for(auto W2_cand : W2Coll)
             for (auto b2_cand : b2Coll) {
		  const reco::Candidate* W2_mother = W2_cand->mother();
		  const reco::Candidate* b2_mother = b2_cand->mother();
	       //   while (W2_mother->pdgId() == -24) W2_mother = W2_mother->mother();
		   if (W2_mother == b2_mother && W2_mother->pdgId() == -6) {
			tbarColl.push_back(W2_mother);
			//std::cout <<" find tbar->bbarW- " << std::endl;
			break;
		    }
              }
    }
    std::cout <<"tColl size " << tColl.size() <<" tbarColl " << tbarColl.size() << std::endl;
    //bool ttbar =false;
    if (tColl.size()==1 and tbarColl.size() ==1){
	//ttbar = true;
	t1cand = tColl.at(0);
	t2cand = tbarColl.at(0);
	//while (t1cand->numberOfMothers()==1 and (t1cand->mother())->pdgId()==6) t1cand = t1cand->mother();
	//while (t2cand->numberOfMothers()==1 and (t2cand->mother())->pdgId()==-6) t2cand = t2cand->mother();
        std::cout << "find t and tbar candidate " << std::endl;
	b1cand = finddecendant(t1cand, 5, false);
	w1cand = finddecendant(t1cand, 24, false);
	b2cand = finddecendant(t2cand, -5, false);
	w2cand = finddecendant(t2cand, -24, false);   
   	//for (unsigned int i=0; i < b1cand->numberOfDaughters(); i++)
   	   //std::cout <<"candidate id "<< b1cand->pdgId()<<" daughter i "<< i <<" id "<<(b1cand->daughter(i))->pdgId()<< std::endl;
	if (hasDaughter(w1cand, -13) and hasDaughter(w2cand, 13)){
	    mu1cand = findmudaughter(w1cand);
	    nu1cand = findnudaughter(w1cand);
	    mu2cand = findmudaughter(w2cand);
	    nu2cand = findnudaughter(w2cand);
	    //make sure all candiates are in same frame
	    while (mu1cand->numberOfDaughters()==1 and  mu1cand->daughter(0)->pdgId()==mu1cand->pdgId())
		    mu1cand = mu1cand->daughter(0);
	    while (nu1cand->numberOfDaughters()==1 and  nu1cand->daughter(0)->pdgId()==nu1cand->pdgId())
		    nu1cand = nu1cand->daughter(0);
	    while (mu2cand->numberOfDaughters()==1 and  mu2cand->daughter(0)->pdgId()==mu2cand->pdgId())
		    mu2cand = mu2cand->daughter(0);
	    while (nu2cand->numberOfDaughters()==1 and  nu2cand->daughter(0)->pdgId()==nu2cand->pdgId())
		    nu2cand = nu2cand->daughter(0);
	    findAllGenParticles = true;
	    std::cout <<" mu1 " ; printCandidate(mu1cand);
	    std::cout <<" nu1 " ; printCandidate(nu1cand);
	    std::cout <<" mu2 " ; printCandidate(mu2cand);
	    std::cout <<" nu2 " ; printCandidate(nu2cand);
	}else{
	    findAllGenParticles = false;
	    std::cout <<"failed to two leptons from  W decays "<< std::endl;
	}
        std::cout <<" w1 " ; printCandidate(w1cand);
        std::cout <<" w2 " ; printCandidate(w2cand);
        std::cout <<" b1 " ; printCandidate(b1cand);
        std::cout <<" b2 " ; printCandidate(b2cand);
    }

}


void 
DiHiggsWWBBAnalyzer::matchGenJet2Parton(edm::Handle<std::vector<reco::GenJet>> genjetColl){
    TLorentzVector b1_p4(b1cand->px(), b1cand->py(), b1cand->pz(), b1cand->energy());
    TLorentzVector b2_p4(b2cand->px(), b2cand->py(), b2cand->pz(), b2cand->energy());
   // float dR_b1genjet = 9;
   // float dR_b2genjet = 9;

    bool hasb1genjet = false, hasb2genjet = false;
    for (const reco::GenJet& genjet: *genjetColl){
	TLorentzVector genjet_p4(genjet.px(), genjet.py(), genjet.pz(), genjet.energy());
	float dR_b1 = genjet_p4.DeltaR(b1_p4);
	float dR_b2 = genjet_p4.DeltaR(b2_p4);
	//std::cout <<"dR_b1 "<<dR_b1 <<" dR_b2 "<< dR_b2 <<" genjet px "<< genjet.px() << " py "<< genjet.py()<<" pz "<< genjet.pz()<<" E "<<genjet.energy() <<std::endl;
	if (dR_b1 < dR_b2 and dR_b1 < dR_b1genjet){
	    b1genjet = &genjet;
	    dR_b1genjet = dR_b1;
	    hasb1genjet = true;
	}else if (dR_b1 > dR_b2 and dR_b2 < dR_b2genjet){
	    b2genjet = &genjet;
	    dR_b2genjet = dR_b2;
	    hasb2genjet = true;
	}
    }
    if (hasb1genjet and hasb2genjet){
    	hastwogenjets = true;
	std::cout <<" b1genjet: genparticle components " << std::endl;
	for (auto genp : b1genjet->getGenConstituents())
		std::cout <<"genp id "<< genp->pdgId()<<" px "<< genp->px()<<" py "<< genp->py()<<" pt "<< genp->pt() << std::endl;
    }

}
// ------------ method called when starting to processes a run  ------------
/*
void 
DiHiggsWWBBAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
DiHiggsWWBBAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
DiHiggsWWBBAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
DiHiggsWWBBAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/


//---------- method called to fill branches --------------------------------------------
void 
DiHiggsWWBBAnalyzer::fillbranches(){
    using namespace std;//include min, max
    if (findAllGenParticles){
      mu1_energy = mu1cand->energy();
      mu1_pt = mu1cand->pt();
      mu1_eta = mu1cand->eta();
      mu1_phi = mu1cand->phi();
      mu1_px = mu1cand->px();
      mu1_py = mu1cand->py();
      mu1_pz = mu1cand->pz();
      nu1_energy = nu1cand->energy();
      nu1_pt = nu1cand->pt();
      nu1_eta = nu1cand->eta();
      nu1_phi = nu1cand->phi();
      nu1_px = nu1cand->px();
      nu1_py = nu1cand->py();
      nu1_pz = nu1cand->pz();
      w1_energy = w1cand->energy();
      w1_pt = w1cand->pt();
      w1_eta = w1cand->eta();
      w1_phi = w1cand->phi();
      w1_px = w1cand->px();
      w1_py = w1cand->py();
      w1_pz = w1cand->pz();
      w1_mass = w1cand->mass();

      mu2_energy = mu2cand->energy();
      mu2_pt = mu2cand->pt();
      mu2_eta = mu2cand->eta();
      mu2_phi = mu2cand->phi();
      mu2_px = mu2cand->px();
      mu2_py = mu2cand->py();
      mu2_pz = mu2cand->pz();
      nu2_energy = nu2cand->energy();
      nu2_pt = nu2cand->pt();
      nu2_px = nu2cand->px();
      nu2_py = nu2cand->py();
      nu2_pz = nu2cand->pz();
      nu2_eta = nu2cand->eta();
      nu2_phi = nu2cand->phi();
      w2_energy = w2cand->energy();
      w2_pt = w2cand->pt();
      w2_eta = w2cand->eta();
      w2_phi = w2cand->phi();
      w2_px = w2cand->px();
      w2_py = w2cand->py();
      w2_pz = w2cand->pz();
      w2_mass = w2cand->mass();

      b1_energy = b1cand->energy();
      b1_pt = b1cand->pt();
      b1_eta = b1cand->eta();
      b1_phi = b1cand->phi();
      b1_px = b1cand->px();
      b1_py = b1cand->py();
      b1_pz = b1cand->pz();
      b2_energy = b2cand->energy();
      b2_pt = b2cand->pt();
      b2_eta = b2cand->eta();
      b2_phi = b2cand->phi();
      b2_px = b2cand->px();
      b2_py = b2cand->py();
      b2_pz = b2cand->pz();
      if (sampleType_<=12 and sampleType_>0){
	  htoWW_energy = htoWWcand->energy();
	  htoWW_px = htoWWcand->px();
	  htoWW_py = htoWWcand->py();
	  htoWW_pz = htoWWcand->pz();
	  htoWW_mass = htoWWcand->mass();


	  htobb_energy = htoBBcand->energy();
	  htobb_px = htoBBcand->px();
	  htobb_py = htoBBcand->py();
	  htobb_pz = htoBBcand->pz();
	  htobb_mass = htoBBcand->mass();
	 
	  
	  h2tohh_energy = h2tohhcand->energy();
	  h2tohh_px = h2tohhcand->px();
	  h2tohh_py = h2tohhcand->py();
	  h2tohh_pz = h2tohhcand->pz();
	  h2tohh_mass = h2tohhcand->mass();
   

      }else if (sampleType_==13){
	  t1_energy = t1cand->energy();
	  t1_px = t1cand->px();
	  t1_py = t1cand->py();
	  t1_pz = t1cand->pz();
	  t2_energy = t2cand->energy();
	  t2_px = t2cand->px();
	  t2_py = t2cand->py();
	  t2_pz = t2cand->pz();
      }
       
      if (hastwogenjets){
	  b1genjet_energy = b1genjet->energy();
	  b1genjet_pt = b1genjet->pt();
	  b1genjet_eta = b1genjet->eta();
	  b1genjet_phi = b1genjet->phi();
	  b1genjet_px = b1genjet->px();
	  b1genjet_py = b1genjet->py();
	  b1genjet_pz = b1genjet->pz();
	  b2genjet_energy = b2genjet->energy();
	  b2genjet_pt = b2genjet->pt();
	  b2genjet_eta = b2genjet->eta();
	  b2genjet_phi = b2genjet->phi();
	  b2genjet_px = b2genjet->px();
	  b2genjet_py = b2genjet->py();
	  b2genjet_pz = b2genjet->pz();
      }
      TLorentzVector mu1_p4(mu1cand->px(), mu1cand->py(), mu1cand->pz(), mu1cand->energy());
      TLorentzVector mu2_p4(mu2cand->px(), mu2cand->py(), mu2cand->pz(), mu2cand->energy());
      TLorentzVector b1_p4(b1cand->px(), b1cand->py(), b1cand->pz(), b1cand->energy());
      TLorentzVector b2_p4(b2cand->px(), b2cand->py(), b2cand->pz(), b2cand->energy());
      dR_genbl   = (b1_p4.Pt()>b2_p4.Pt()) ? (b1_p4.DeltaR( (mu1_p4.Pt()>mu2_p4.Pt()) ? mu1_p4 : mu2_p4 )) : (b2_p4.DeltaR( (mu1_p4.Pt()>mu2_p4.Pt()) ? mu1_p4 : mu2_p4 ));
      dR_genb1l1 = b1_p4.DeltaR(mu1_p4);
      dR_genb1l2 = b1_p4.DeltaR(mu2_p4);
      dR_genb2l1 = b2_p4.DeltaR(mu1_p4);
      dR_genb2l2 = b2_p4.DeltaR(mu2_p4);
      dR_genb1b2 = b1_p4.DeltaR(b2_p4);
      dR_genl1l2 = mu1_p4.DeltaR(mu2_p4);
      dR_genl1l2b1b2 = (mu1_p4+mu2_p4).DeltaR(b1_p4+b2_p4);
      dphi_genl1l2b1b2 = TVector2::Phi_mpi_pi( (mu1_p4+mu2_p4).Phi()-(b1_p4+b2_p4).Phi() );
      TLorentzVector genll_p4 = mu1_p4+mu2_p4;
      TLorentzVector genbb_p4 = b1_p4+b2_p4;
      dR_genminbl = min(min(dR_genb1l1,dR_genb1l2),min(dR_genb2l1,dR_genb2l2));
      dphi_genllbb = TVector2::Phi_mpi_pi(genll_p4.Phi()-genbb_p4.Phi());
      dphi_genllmet = TVector2::Phi_mpi_pi(genll_p4.Phi()-genmet_phi);
      mass_genl1l2 = genll_p4.M(); energy_genl1l2 = genll_p4.Energy(); pt_genl1l2 = genll_p4.Pt(); eta_genl1l2 = genll_p4.Eta(); phi_genl1l2 = genll_p4.Phi();
      mass_genb1b2 = genbb_p4.M(); energy_genb1b2 = genbb_p4.Energy(); pt_genb1b2 = genbb_p4.Pt(); eta_genb1b2 = genbb_p4.Eta(); phi_genb1b2 = genbb_p4.Phi();
      mass_gentrans = sqrt(2*genll_p4.Pt()*genmet_pt*(1-cos(dphi_genllmet)));
   }
}


//-------------- method called to find a muon daughter for given candidate -------------------------
const reco::Candidate* 
DiHiggsWWBBAnalyzer::findmudaughter(const reco::Candidate *Wcand){

    const reco::Candidate *tmpcand = NULL;
    int count = 0;
    for (unsigned int i = 0; i<Wcand->numberOfDaughters(); i++)
		if (abs(Wcand->daughter(i)->pdgId()) == 13) {
			count++;
            tmpcand = Wcand->daughter(i);
    	}
    if (count != 1) std::cout << " this candidate has more one mu daughter " << std::endl;
    return tmpcand;
}


//------------ method called to find muon descendants for given candidate -------------------------
// return one muon descendant and the number of muon descendants by "count"
const reco::Candidate*
DiHiggsWWBBAnalyzer::findmudescendants(const reco::Candidate *cand, int& count){

   const reco::Candidate * tmpcand = NULL;
   for (unsigned int i = 0; i<cand->numberOfDaughters(); i++){
         tmpcand = findmudescendants(cand->daughter(i), count);
         if (abs(cand->daughter(i)->pdgId()) == 13) count++;
         if (tmpcand == NULL && abs(cand->daughter(i)->pdgId()) == 13)  tmpcand = cand->daughter(i);                            
   }

   return tmpcand;
}


//-------------- method called to find a neutrino daughter for given candidate -------------------------
const reco::Candidate* 
DiHiggsWWBBAnalyzer::findnudaughter(const reco::Candidate *Wcand){

    const reco::Candidate *tmpcand = NULL;
    int count = 0;
    for (unsigned int i = 0; i<Wcand->numberOfDaughters(); i++)
        if (abs(Wcand->daughter(i)->pdgId()) == 14) {
			count++;
            tmpcand = Wcand->daughter(i);
        }
    if (count != 1) std::cout << " this candidate has more one nu daughter " << std::endl;
    return tmpcand;
}

//------------ method called to find neutrino descendants for given candidate -------------------------
// return one neutrino descendant and the number of neutrno descendants by "count"
const reco::Candidate*
DiHiggsWWBBAnalyzer::findnudescendants(const reco::Candidate *cand, int& count){

   const reco::Candidate * tmpcand = NULL;
   for (unsigned int i = 0; i<cand->numberOfDaughters(); i++){
         tmpcand = findmudescendants(cand->daughter(i), count);
         if (abs(cand->daughter(i)->pdgId()) == 14) count++;
         if (tmpcand == NULL && abs(cand->daughter(i)->pdgId()) == 14)  tmpcand = cand->daughter(i);                            
   }
   return tmpcand;
}


//------------- method called to find stable decendant with pdgid = id
const reco::Candidate* 
DiHiggsWWBBAnalyzer::stabledecendant(const reco::Candidate* cand, int id){
   const reco::Candidate* tmp = NULL;
   for (unsigned int i=0; i < cand->numberOfDaughters(); i++){
	if ((cand->daughter(i))->pdgId() == id && (cand->daughter(i))->status() == 1)
		return 	tmp=cand->daughter(i);
        else if (stabledecendant(cand->daughter(i),id)) 
		return tmp=stabledecendant(cand->daughter(i),id);
   }   
   return tmp;
}


//------------ method called to find decendant with pdgid = id, 
//if first it true, then return the candidate closest to seed
//if first it false, then return the candidate farthest to seed
const reco::Candidate* 
DiHiggsWWBBAnalyzer::finddecendant(const reco::Candidate* cand, int id, bool first){
   const reco::Candidate* tmp = NULL;
   for (unsigned int i=0; i < cand->numberOfDaughters(); i++){
        
	if ((cand->daughter(i))->pdgId() == id && first && cand->pdgId() != id)
		return 	tmp=cand->daughter(i);
	else if ((cand->daughter(i))->pdgId() == id && !first && !hasDaughter(cand->daughter(i), id)) 
		return  tmp=cand->daughter(i); // tmp does not has daughter with pdgid = id
	else if ((cand->daughter(i))->pdgId() == id && !first && (cand->daughter(i))->numberOfDaughters()>1) 
		return  tmp=cand->daughter(i);// tmp has more one daughters therefore it is final-states
        else if (finddecendant(cand->daughter(i),id, first)) 
		return tmp=finddecendant(cand->daughter(i),id);
   }
    return tmp;
}

//---------- method called to find a ancestor with pdgid = id, 
//if first is true, then return the candidate closest to seed
//if first is false, then return the candidate furthest to seed
const reco::Candidate*
DiHiggsWWBBAnalyzer::findancestor(const reco::Candidate* cand, int id, bool first){

   const reco::Candidate* tmp = NULL;
   for (unsigned int i=0; i < cand->numberOfMothers(); i++){
        
	if ((cand->mother(i))->pdgId() == id && first && cand->pdgId() != id)
		return 	tmp=cand->mother(i);
	else if ((cand->mother(i))->pdgId() == id && !first && !hasMother(cand->mother(i), id)) 
		return  tmp=cand->mother(i);
        else if (findancestor(cand->mother(i),id, first)) 
		return tmp=findancestor(cand->mother(i),id, first);
   }
   return tmp;

}

//---------- method called to check whether cand has mother with pdgid = id -----------------------------------
bool
DiHiggsWWBBAnalyzer::hasMother(const reco::Candidate* cand, int id){

   for (unsigned int i=0; i < cand->numberOfMothers(); i++)
        if ((cand->mother(i))->pdgId() == id) return true;
   return false;

}

//-------- method called to check whether cand has daughter with pdgid = id ------------------------------------
bool
DiHiggsWWBBAnalyzer::hasDaughter(const reco::Candidate* cand, int id){
 
   //for (unsigned int i=0; i < cand->numberOfDaughters(); i++)
   //	   std::cout <<"candidate id "<< cand->pdgId()<<" daughter i "<< i <<" id "<<(cand->daughter(i))->pdgId()<< std::endl;
   for (unsigned int i=0; i < cand->numberOfDaughters(); i++)
        if ((cand->daughter(i))->pdgId() == id) return true;
   return false;

}


//----------- method called to calculate total lorentz vector from b bar jets ---------------
//
TLorentzVector
DiHiggsWWBBAnalyzer::stableDecendantsLorentz(const reco::Candidate *cand){
  
   TLorentzVector bjets = TLorentzVector();
   for (unsigned i = 0; i < cand->numberOfDaughters(); i++){
	if ((cand->daughter(i))->status() == 1){
		TLorentzVector tmp((cand->daughter(i))->px(), (cand->daughter(i))->py(), 
				(cand->daughter(i))->pz(),(cand->daughter(i))->energy());
                //printCandidate(cand->daughter(i));
        	bjets = bjets+tmp; 
	}else bjets = bjets+stableDecendantsLorentz(cand->daughter(i));
   }

   return bjets;
}

//-------------  method called to calculate MET in simuation ------------------ 
TLorentzVector 
DiHiggsWWBBAnalyzer::calculateMET(){

   TLorentzVector METlorentz = TLorentzVector();
   TVector2 met_pxpy(nu1cand->px()+nu2cand->px(), nu1cand->py()+nu2cand->py());
   METlorentz.SetPxPyPzE(nu1cand->px()+nu2cand->px(), nu1cand->py()+nu2cand->py(),0,met_pxpy.Mod());

   return METlorentz;
}




// ------------ method called for printing additional information, useful for debugging  ------------
void
DiHiggsWWBBAnalyzer::print() {
   //todo: combine all necessary print out for debug 
    std::cout << "print() " << std::endl;
}

//------------- method called for printing h->WW->mumununu chain -------------------------
/*
void 
DiHiggsWWBBAnalyzer::printHtoWWChain(){
          
    float h_energy = mu1_energy+mu2_energy+nu1_energy+nu2_energy;
    float h_px = mu1_px+mu2_px+nu1_px+nu2_px;
    float h_py = mu1_py+mu2_py+nu1_py+nu2_py;
    float h_pz = mu1_pz+mu2_pz+nu1_pz+nu2_pz;
          
    TLorentzVector final_p4(h_px, h_py, h_pz, h_energy);
    
    float h_energy_bb = b1_energy+b2_energy; 
    float h_px_bb = b1_px+b2_px; 
    float h_py_bb = b1_py+b2_py;
    float h_pz_bb = b1_pz+b2_pz;
    TLorentzVector h_bb_p4(h_px_bb,h_py_bb,h_pz_bb,h_energy_bb);

    TLorentzVector h2_final_tot(final_p4+h_bb_p4); 
    if (htoWW){
             
               TLorentzVector h_p4(htoWW_px, htoWW_py, htoWW_pz, htoWW_energy);
               TLorentzVector WW_p4(mu1_W1_cand->px()+mu2_W2_cand->px(),mu1_W1_cand->py()+mu2_W2_cand->py(), 
                                    mu1_W1_cand->pz()+mu2_W2_cand->pz(),mu1_W1_cand->energy()+mu2_W2_cand->energy());
               std::cout << "invariant mass from h_p4: " << h_p4.M() 
                         << "	, from WW_p4 " << WW_p4.M() 
                         << "	, from final_p4 " << final_p4.M() << std::endl;
               if (abs(WW_p4.M()-h_p4.M())>1)  std::cout << "h->WW, invariant mass reconstruction discrepancy ? " << std::endl; 
               std::cout <<  " H -> WW " << std::endl;
               const reco::Candidate* tmp_cand1 = NULL;
               const reco::Candidate* tmp_cand2 = NULL;
               
               for (unsigned int n = 0; n<mu1_htoWW_cand->numberOfDaughters(); n++){                                  
                      if ((mu1_htoWW_cand->daughter(n))->pdgId()==-24) tmp_cand1 =  mu1_htoWW_cand->daughter(n);
                      if ((mu1_htoWW_cand->daughter(n))->pdgId()== 24) tmp_cand2 =  mu1_htoWW_cand->daughter(n);
                      if (n >= 2) std::cout << "h has more 2 daughters, id " << (mu1_htoWW_cand->daughter(n))->pdgId() << std::endl;
                  }
               while (tmp_cand1 != mu1_W1_cand || tmp_cand2 != mu2_W2_cand){
                TLorentzVector W1_p4(tmp_cand1->px(),tmp_cand1->py(),tmp_cand1->pz(), tmp_cand1->energy());
                TLorentzVector W2_p4(tmp_cand2->px(),tmp_cand2->py(),tmp_cand2->pz(), tmp_cand2->energy());
                TLorentzVector tmp_WW_p4(W1_p4+W2_p4);
                std::cout <<"W- num of daughters "<< tmp_cand1->numberOfDaughters() << " W-mass " << W1_p4.M() << " W1 four momentum "; W1_p4.Print();
                std::cout <<"W+ num of daughters "<< tmp_cand2->numberOfDaughters() << " W+mass " << W2_p4.M() << " W2 four momentum "; W2_p4.Print();
                std::cout << "Total invariant mass " << tmp_WW_p4.M() <<" tmp_WW four momentum "; tmp_WW_p4.Print(); 
                //if (tmp_cand1 != mu1_W1_cand) {
                     for (unsigned int i = 0; i<tmp_cand1->numberOfDaughters(); i++){
                          std::cout << " daughter of W- , id " << tmp_cand1->daughter(i)->pdgId() << "  status " << tmp_cand1->daughter(i)->status() <<std::endl; 
                          TLorentzVector dau_W1_p4(tmp_cand1->daughter(i)->px(),tmp_cand1->daughter(i)->py(),tmp_cand1->daughter(i)->pz(),tmp_cand1->daughter(i)->energy());
                          std::cout << " four momentum "; dau_W1_p4.Print();
                          if (tmp_cand1->daughter(i)->pdgId() == -24) tmp_cand1 = tmp_cand1->daughter(i);  
                          }
                // }
                //if (tmp_cand2 != mu2_W2_cand) {
                     for (unsigned int j = 0; j<tmp_cand2->numberOfDaughters(); j++){
                          std::cout << " daughter of W+ , id " << tmp_cand2->daughter(j)->pdgId() << "  status " << tmp_cand2->daughter(j)->status() <<std::endl; 
                          TLorentzVector dau_W2_p4(tmp_cand2->daughter(j)->px(),tmp_cand2->daughter(j)->py(),tmp_cand2->daughter(j)->pz(),tmp_cand2->daughter(j)->energy());
                          std::cout << " four momentum "; dau_W2_p4.Print();
                          if (tmp_cand2->daughter(j)->pdgId() == 24) tmp_cand2 = tmp_cand2->daughter(j);  
                          }
                // }
               }
               std::cout << "WW -> mumununu" << std::endl;
               while (tmp_cand1->status() != 1 || tmp_cand2->status() != 1){
                     
                     std::cout << "-------------  begin of this loop ----------------------" << std::endl;
     
                     int size1 = tmp_cand1->numberOfDaughters();
                     int size2 = tmp_cand2->numberOfDaughters();
                     int muon1 = -1;
                     int muon2 = -1;
                     float px=0;
                     float py=0;
                     float pz=0;
                     float energy=0;
                      if (tmp_cand1->pdgId() == 13)  {
                             //std::cout << "cand1 reaches final states, status of particle " << tmp_cand1->status() << std::endl;
                             px += nu1_px;
                             py += nu1_py;
                             pz += nu1_pz;
                             energy += nu1_energy;
                             }
                     std::cout << "cand1, id"<< tmp_cand1->pdgId() << " status " << tmp_cand1->status() <<" size of daughters " << size1 << std::endl; 
                     std::cout << " daughters of " << ((tmp_cand1->pdgId()==-24)?"W-  ":"muon- ") << std::endl;
                       for (int i = 0; i < size1; i++){
                             const Candidate * d1 = tmp_cand1->daughter(i); 
                             std::cout << "daughter id " << d1->pdgId() << "  status " << d1->status() << std::endl;
                             if (d1->pdgId() == 13 ) muon1 = i;
                             printCandidate(d1);
                             px += d1->px();
                             py += d1->py();
                             pz += d1->pz();
                             energy += d1->energy(); 
                   }
                      TLorentzVector cand1_lorentz(px, py, pz, energy); 
                      std::cout << " W- mass from W- Candidate " << w1_mass << " from mu-,nu " << cand1_lorentz.M() << std::endl;
                      if (muon1 != -1 && tmp_cand1->status() != 1) tmp_cand1 = tmp_cand1->daughter(muon1);
                      float px2 = 0.0;
                      float py2 = 0.0;
                      float pz2 = 0.0;
                      float energy2 = 0.0;
                      if (tmp_cand2->pdgId() == -13)  {
                           //  std::cout << "cand2 reaches final states, status of particle "<< tmp_cand2->status() << std::endl;
                             px2 += nu2_px;
                             py2 += nu2_py;
                             pz2 += nu2_pz;
                             energy2 += nu2_energy;
                             }
                     std::cout << "cand2, id" << tmp_cand2->pdgId() <<" status " << tmp_cand2->status() <<" size of daughters " << size2 << std::endl; 
                     std::cout << " daughters of " << ((tmp_cand2->pdgId()==24)?"W+  ":"muon+ ") << std::endl;
                       for (int j = 0; j < size2; j++){
                             const Candidate * d2 = tmp_cand2->daughter(j); 
                             std::cout << "daughter id " << d2->pdgId() << "  status " << d2->status() << std::endl;
                             if (d2->pdgId() == -13 ) muon2 = j;
                             printCandidate(d2);
                             px2 += d2->px();
                             py2 += d2->py();
                             pz2 += d2->pz();
                             energy2 += d2->energy(); 
                   }
                      TLorentzVector cand2_lorentz(px2, py2, pz2, energy2); 
                      std::cout << " W+ mass from W+ Candidate " << w2_mass << " from mu+,nu " << cand2_lorentz.M() << std::endl;
                      if (muon2 != -1 && tmp_cand2->status() != 1) tmp_cand2 = tmp_cand2->daughter(muon2);
                     TLorentzVector tmp = cand1_lorentz+cand2_lorentz;
                     std::cout <<"Total px " << tmp.Px() << " py " << tmp.Py() << " pz " << tmp.Pz() << " E " << tmp.Energy() << std::endl;
                     std::cout << " invariant mass from daughters of WW " << tmp.M() << std::endl;  
                     std::cout << "For Next loop status of cand1 " << tmp_cand1->status()  << "	cand2 " << tmp_cand2->status() << std::endl; 
                     std::cout << "-------------  end of this loop ----------------------" << std::endl;
              } 
            if (findAllGenParticles){	
			std::cout <<"htoWW invariant mass " << final_p4.M() <<" four momentum " << std::endl; 
                        final_p4.Print();
                        std::cout <<"htobb invariant mass " << h_bb_p4.M() << " four momentum " << std::endl;
                        h_bb_p4.Print();
                        std::cout <<"h2tohh invariant mass " << h2_final_tot.M() <<" total momentum " << std::endl;
                        h2_final_tot.Print();
			} 

        }//end if (htoWW)

}*/



//---------- method called to print candidates for debug ---------------------
void
DiHiggsWWBBAnalyzer::printCandidate(const reco::Candidate* cand){

   std::cout <<" Candidate id: "<< cand->pdgId() << " mass: " << cand->mass() <<" (P,E)= ("<< cand->px() <<", "<< cand->py()<<", "<< cand->pz()<<", "<< cand->energy() <<")" <<"(Pt,E) = ("<< cand->pt() <<", "<< cand->eta() <<", "<< cand->phi()<<", "<<cand->energy()<< ")" <<" status: " << cand->status() << std::endl;

}

//--------- method called to print all decendants for cand -------------------
void 
DiHiggsWWBBAnalyzer::printallDecendants(const reco::Candidate* cand){
   
   if (cand->status() != 0 && cand->numberOfDaughters() > 0){
        std::cout << "******************  children of id "<< cand->pdgId() <<"      *********************" << std::endl;
   	for (unsigned int i=0; i < cand->numberOfDaughters(); i++)
        	printCandidate(cand->daughter(i));
        std::cout << "***********************************************************" << std::endl;
   	for (unsigned int i=0; i < cand->numberOfDaughters(); i++)
		printallDecendants(cand->daughter(i));

    }
}


//--------- method called to print all Ancestors for cand -------------------
void 
DiHiggsWWBBAnalyzer::printallAncestors(const reco::Candidate* cand){
   
   if (cand->status() != 0 && cand->numberOfMothers() > 0){
        std::cout << "******************  mothers of id "<< cand->pdgId() <<"      *********************" << std::endl;
   	for (unsigned int i=0; i < cand->numberOfMothers(); i++)
        	printCandidate(cand->mother(i));
        std::cout << "***********************************************************" << std::endl;
   	for (unsigned int i=0; i < cand->numberOfMothers(); i++)
		printallAncestors(cand->mother(i));

    }
}




// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DiHiggsWWBBAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiHiggsWWBBAnalyzer);


