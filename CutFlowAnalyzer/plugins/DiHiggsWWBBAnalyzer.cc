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
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

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
    
   //runMMC
   /*private: 
      void runMMC();
      void initTree(TTree* mmctree);
        
     
      float genEtaGuass(float mean, float rms);
      float genPhiFlat();
      EtaPhi generatenu1_etaphi();
      float nu1pt_onshellW(EtaPhi nu1_etaphi, TLorentzVector* mu1lorentz, float wMass);
      bool  nulorentz_offshellW(TLorentzVector* jetlorentz, TLorentzVector* mu1lorentz, 
			       TLorentzVector* mu2lorentz, TLorentzVector* nu1lorentz, 
 			       TLorentzVector* nu2lorentz, int control, float hMass);
      bool checkSolution(TLorentzVector* jetslorentz,
                          TLorentzVector* mu1lorentz,
                          TLorentzVector* mu2lorentz,
                          TLorentzVector* nu1lorentz, int control, float hMass); 
      bool cutsCheck();
      void assignMuLorentzVec(int control);  
        */  
   private:
      edm::ParameterSet cfg_;
      edm::ParameterSet mmcset_;
      // Labels to access
      edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;        // reconstructed muons
      edm::EDGetTokenT<pat::ElectronCollection> electronToken_;        // reconstructed electron
      edm::EDGetTokenT<pat::JetCollection> jetToken_;
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
      float mu1_eta_;
      float mu2_eta_;
      float mu1_pt_;
      float mu2_pt_;
      float jet1_eta_;
      float jet2_eta_;
      float jet1_pt_;
      float jet2_pt_;
      float met_;
      float jetleptonDeltaR_;
      float leptonIso_;
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
      void checkGenParticlesSingal(edm::Handle<reco::GenParticleCollection> genParticleColl);
      void checkGenParticlesTTbar(edm::Handle<reco::GenParticleCollection> genParticleColl);
      void matchMuon2Gen();//match pat:::Muon to gen muon 
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
      float mu1_energy;
      float mu1_px;
      float mu1_py;
      float mu1_pz;
      float mu1_eta;
      float mu1_phi;
      float w1_energy;
      float w1_px;
      float w1_py;
      float w1_pz;
      float w1_mass;
      float nu1_energy;
      float nu1_px;
      float nu1_py;
      float nu1_pz;
      float nu1_eta;
      float nu1_phi;
      bool Wtomu1nu1;

      float mu2_energy;
      float mu2_px;
      float mu2_py;
      float mu2_pz;
      float mu2_eta;
      float mu2_phi;
      float w2_energy;
      float w2_px;
      float w2_py;
      float w2_pz;
      float w2_mass;
      float nu2_energy;
      float nu2_px;
      float nu2_py;
      float nu2_pz;
      float nu2_eta;
      float nu2_phi;
      bool Wtomu2nu2;
     
      float htoWW_energy;
      float htoWW_px;
      float htoWW_py;
      float htoWW_pz; 
      float htoWW_mass;
    //  float w1_mass;
    //  float w2_mass;
      float b1_energy;
      float b1_px;
      float b1_py;
      float b1_pz;
      int b1_motherid;
      float b2_energy;
      float b2_px;
      float b2_py;
      float b2_pz;
      int b2_motherid;
      float bjet_energy;
      float bjet_px;
      float bjet_py;
      float bjet_pz;
      float bjet_mass;
      float bbarjet_energy;
      float bbarjet_px;
      float bbarjet_py;
      float bbarjet_pz;
      float bbarjet_mass;

      float met;
      float met_phi;
      float met_px;
      float met_py;

      float htobb_energy;
      float htobb_px;
      float htobb_py;
      float htobb_pz;
      float htobb_mass;
      
      float h2tohh_energy;
      float h2tohh_px;
      float h2tohh_py;
      float h2tohh_pz;
      float h2tohh_mass;
      //cuts for higgstoWWbb
      bool mu_positive;
      bool mu_negative;
      bool nu_positive;
      bool nu_negative;
      bool bquark;
      bool bbarquark;
      bool htobb;
      bool htoWW;
      bool findAllGenParticles;
      float virtualW_lowM;
      float virtualW_highM;

     
    private:
      bool runMMC_;
      bool simulation_;
      MMC* thismmc;
     // MMC tree branches
     /*
    private:
      float onshellWMassRandomWalk(float x0, float step, float random);
      float onshellWMassRandomWalk(float x0, float step, float random, TH1F* hist);
      float onshellWMassPDF(float wmass);
  
    private:
      TH1F* readoutonshellWMassPDF();
      TH1F* readoutoffshellWMassPDF();
      TH1F* readoutonshellnuptPDF();
 
    private:
      float weightfromhist(TH1F* pdf, float x); 
      float weightfromonshellnupt(float nupt); 
   
    private:
      bool weightfromonshellnupt_func_;
      bool weightfromonshellnupt_hist_;
      bool weightfromoffshellWmass_hist_;

    private:
      int iterations_;
      int seed_;
      std::string RefPDFfile_;
      float eta_mean;
      float eta_rms;
      float eta_gen; 
      float phi_gen;
      float wmass_gen;
      float hmass_gen;
      TLorentzVector* mu_onshellW_lorentz;
      TLorentzVector* mu_offshellW_lorentz;
      TVector2* MMCmet_vec2;
      TLorentzVector* nu_onshellW_lorentz;
      TLorentzVector* nu_offshellW_lorentz;
      TLorentzVector* offshellW_lorentz;
      TLorentzVector* onshellW_lorentz;
      TLorentzVector* htoWW_lorentz;
      TLorentzVector* htoBB_lorentz;
      TLorentzVector* h2tohh_lorentz;
      
      int control;
      float weight;
      float weight1;//extra weight
      float weight2;//extra weight
 
      float mu_onshellW_Eta;
      float mu_onshellW_Phi;
      float mu_onshellW_Pt;
      float mu_onshellW_E;
      float mu_offshellW_Eta;
      float mu_offshellW_Phi;
      float mu_offshellW_Pt;
      float mu_offshellW_E;
      float nu_onshellW_Eta;
      float nu_onshellW_Phi;
      float nu_onshellW_Pt;
      float nu_onshellW_E;
      float nu_offshellW_Eta;
      float nu_offshellW_Phi;
      float nu_offshellW_Pt;
      float nu_offshellW_E;
      
      float onshellW_Eta;
      float onshellW_Phi;
      float onshellW_Pt;
      float onshellW_E;
      float onshellW_Mass;
      float offshellW_Eta;
      float offshellW_Phi;
      float offshellW_Pt;
      float offshellW_E;
      float offshellW_Mass;
     
      float htoBB_Eta;
      float htoBB_Phi;
      float htoBB_Pt;
      float htoBB_E;
      float htoBB_Mass;
      float htoWW_Eta;
      float htoWW_Phi;
      float htoWW_Pt;
      float htoWW_E;
      float htoWW_Mass;

      float MMCmet_E;
      float MMCmet_Phi;
      float MMCmet_Px;
      float MMCmet_Py;

      float h2tohh_Eta;
      float h2tohh_Phi;
      float h2tohh_Pt;
      float h2tohh_E;
      float h2tohh_Mass;


      float eta_nuoffshellW_true;
      float phi_nuoffshellW_true;
      float pt_nuoffshellW_true;
      float eta_nuonshellW_true;
      float phi_nuonshellW_true;
      float pt_nuonshellW_true;
      float mass_offshellW_true;
      float mass_onshellW_true;
      float pt_h2tohh_true;
*/
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
    mu1_eta_ = iConfig.getUntrackedParameter<double>("mu1_eta",2.4);
    mu2_eta_ = iConfig.getUntrackedParameter<double>("mu2_eta",2.4);
    mu1_pt_ = iConfig.getUntrackedParameter<double>("mu1_pt",10);
    mu2_pt_ = iConfig.getUntrackedParameter<double>("mu2_pt",10);
    jet1_eta_ = iConfig.getUntrackedParameter<double>("jet1_eta",2.5);
    jet2_eta_ = iConfig.getUntrackedParameter<double>("jet2_eta",2.5);
    jet1_pt_ = iConfig.getUntrackedParameter<double>("jet1_pt",20);
    jet2_pt_ = iConfig.getUntrackedParameter<double>("jet2_pt",20);
	//mmcset_ = iConfig.getParameter<edm::ParameterSet>("mmcset"); 
    sampleType_ = iConfig.getUntrackedParameter<int>("SampleType",0);
    finalStates_ = iConfig.getParameter<bool>("finalStates");
    runMMC_ = iConfig.getParameter<bool>("runMMC");
    simulation_ = iConfig.getParameter<bool>("simulation");
/*
weightfromonshellnupt_func_ = iConfig.getParameter<bool>("weightfromonshellnupt_func");
weightfromonshellnupt_hist_ = iConfig.getParameter<bool>("weightfromonshellnupt_hist");
     weightfromoffshellWmass_hist_ = iConfig.getParameter<bool>("weightfromoffshellWmass_hist");
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
    mu1_energy = 0.0;
    mu1_px = 0.0;
    mu1_py = 0.0;
    mu1_pz = 0.0;
    w1_energy = 0.0;
    w1_px = 0.0;
    w1_py = 0.0;
    w1_pz = 0.0;
    w1_mass = 0.0;
    nu1_energy = 0.0;
    nu1_px = 0.0;
    nu1_py = 0.0;
    nu1_pz = 0.0;
    Wtomu1nu1 = false;

    mu2_energy = 0.0;
    mu2_px = 0.0;
    mu2_py = 0.0;
    mu2_pz = 0.0;
    w2_energy = 0.0;
    w2_px = 0.0;
    w2_py = 0.0;
    w2_pz = 0.0;
    w2_mass = 0.0;
    nu2_energy = 0.0;
    nu2_px = 0.0;
    nu2_py = 0.0;
    nu2_pz = 0.0;
    Wtomu2nu2 = false;

    htoWW_energy = 0.0;
    htoWW_px = 0.0;
    htoWW_py = 0.0;
    htoWW_pz = 0.0;
    htoWW_mass = 0.0;

    b1_energy = 0.0;
    b1_px = 0.0;
    b1_py = 0.0;
    b1_pz = 0.0;
    b1_motherid = 0;
    b2_energy = 0.0;
    b2_px = 0.0;
    b2_py = 0.0;
    b2_pz = 0.0;
    b2_motherid = 0;
    bjet_energy = 0.0;
    bjet_px = 0.0;
    bjet_py = 0.0;
    bjet_pz = 0.0;
    bjet_mass = 0.0;
    bbarjet_energy = 0.0;
    bbarjet_px = 0.0;
    bbarjet_py = 0.0;
    bbarjet_pz = 0.0;
    bbarjet_mass = 0.0;

    htobb_energy = 0.0;
    htobb_px = 0.0;
    htobb_py = 0.0;
    htobb_pz = 0.0;
    htobb_mass = 0.0;

    h2tohh_energy = 0.0;
    h2tohh_px = 0.0;
    h2tohh_py = 0.0;
    h2tohh_pz = 0.0;
    h2tohh_mass = 0.0;

    mu_positive = false;
    mu_negative = false;
    bquark = false;
    bbarquark = false;
    htobb = false;
    htoWW = false;
    virtualW_lowM = 25;
    virtualW_highM = 45;
    findAllGenParticles = false;

     /* 
      mu1_lorentz = new TLorentzVector();
      nu1_lorentz = new TLorentzVector();
      mu2_lorentz = new TLorentzVector();
      nu2_lorentz = new TLorentzVector();
      bbar_lorentz = new TLorentzVector();
      met_lorentz = new TLorentzVector();

     // runMMC
      mu_onshellW_lorentz = new TLorentzVector();
      mu_offshellW_lorentz = new TLorentzVector();
      MMCmet_vec2 = new TVector2();
      nu_onshellW_lorentz = new TLorentzVector();
      nu_offshellW_lorentz = new TLorentzVector();
      offshellW_lorentz = new TLorentzVector();
      onshellW_lorentz = new TLorentzVector();
      htoWW_lorentz = new TLorentzVector();
      htoBB_lorentz = new TLorentzVector();
      h2tohh_lorentz = new TLorentzVector();

      mu_onshellW_Eta =0;
      mu_onshellW_Phi =0;
      mu_onshellW_Pt =0;
      mu_onshellW_E =0;
      mu_offshellW_Eta =0;
      mu_offshellW_Phi =0;
      mu_offshellW_Pt =0;
      mu_offshellW_E =0;
      nu_onshellW_Eta =0;
      nu_onshellW_Phi =0;
      nu_onshellW_Pt =0;
      nu_onshellW_E =0 ;
      nu_offshellW_Eta =0;
      nu_offshellW_Phi =0; 
      nu_offshellW_Pt =0;
      nu_offshellW_E =0;
      
      onshellW_Eta =0;
      onshellW_Phi =0;
      onshellW_Pt =0;
      onshellW_E =0;
      onshellW_Mass =0;
      offshellW_Eta =0;
      offshellW_Phi =0;
      offshellW_Pt =0;
      offshellW_E =0;
      offshellW_Mass =0;
     
      htoBB_Eta =0;
      htoBB_Phi =0;
      htoBB_Pt =0;
      htoBB_E =0;
      htoBB_Mass =0;
      htoWW_Eta =0;
      htoWW_Phi =0;
      htoWW_Pt =0;
      htoWW_E =0;
      htoWW_Mass =0;

      h2tohh_Eta =0;
      h2tohh_Phi =0;
      h2tohh_Pt =0;
      h2tohh_E =0;
      h2tohh_Mass =0;
      eta_nuoffshellW_true = 0.0;
      pt_nuoffshellW_true = 0.0;
      eta_nuonshellW_true = 0.0;
      pt_nuonshellW_true = 0.0;
      mass_offshellW_true = 0.0;
      mass_onshellW_true =0.0;
      pt_h2tohh_true = 0;
      
  */  

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
    
    std::cout << "event  " << iEvent.id().event() << std::endl;
    ievent = iEvent.id().event();
    edm::Handle<reco::GenParticleCollection> genParticleColl;
    iEvent.getByToken(genParticlesToken_, genParticleColl);
    if (sampleType_<=12 and sampleType_>0)
	checkGenParticlesSingal(genParticleColl);
    else if (sampleType_==13)
	checkGenParticlesTTbar(genParticleColl);

  // Cut on primary vertex in event
    edm::Handle<reco::VertexCollection> primaryVertices;
    iEvent.getByToken(primaryVerticesToken_, primaryVertices);
    if (primaryVertices->empty()) return; // skip the event if no PV found
    const reco::Vertex &PV = primaryVertices->front();
   
    edm::Handle<pat::MuonCollection> muons;
    iEvent.getByToken(muonToken_, muons);
    edm::Handle<pat::ElectronCollection> electrons;
    iEvent.getByToken(electronToken_, electrons);
    std::vector<const reco::Candidate *> leptons;
    for (const pat::Muon &mu : *muons) {
	leptons.push_back(&mu);
	printf("muon with pt %4.1f, dz(PV) %+5.3f, POG loose id %d, tight id %d\n",
		mu.pt(), mu.muonBestTrack()->dz(PV.position()), mu.isLooseMuon(), mu.isTightMuon(PV));
    }
    for (const pat::Electron &el : *electrons) leptons.push_back(&el);
    for (const reco::Candidate *lep : leptons) {
	if (lep->pt() < 5) continue;
	printf("lepton: pt %5.1f, eta %+4.2f \n", lep->pt(), lep->eta());
    }

    //for (pat::MuonCollection::const_iterator iMuon = muons->begin();  iMuon != muons->end();  ++iMuon) {
    //or for (const pat::Muon &mu : *muons) {
      
    edm::Handle<pat::JetCollection> jets;
    iEvent.getByToken(jetToken_, jets);
    for (const pat::Jet &j : *jets) {
	if (j.pt() < 30 || fabs(j.eta()) > 2.4) continue;
	printf("Jet with pt %6.1f, eta %+4.2f, pileup mva disc %+.2f, btag CSV %.3f, CISV %.3f\n",
		j.pt(),j.eta(), j.userFloat("pileupJetId:fullDiscriminant"), std::max(0.f,j.bDiscriminator("combinedSecondaryVertexBJetTags")), std::max(0.f,j.bDiscriminator("combinedInclusiveSecondaryVertexBJetTags")));
    }


    edm::Handle<pat::METCollection> mets;
    iEvent.getByToken(metToken_, mets);
    const pat::MET &met = mets->front();
    printf("MET: pt %5.1f, phi %+4.2f, sumEt (%.1f). genMET %.1f. MET with JES up/down: %.1f/%.1f\n",
   	 met.pt(), met.phi(), met.sumEt(), met.genMET()->pt(),met.shiftedPt(pat::MET::JetEnUp), met.shiftedPt(pat::MET::JetEnDown));
   /* 
   //if (h2tohh && runMMC_) runMMC();
   if (h2tohh && runMMC_){
   	mu1_lorentz.SetPtEtaPhiM(mu1cand->pt(), mu1cand->eta(), mu1cand->phi(), 0);
   	nu1_lorentz.SetPtEtaPhiM(nu1cand->pt(), nu1cand->eta(), nu1cand->phi(), 0);
        mu2_lorentz.SetPtEtaPhiM(mu2cand->pt(), mu2cand->eta(), mu2cand->phi(), 0); 
        nu2_lorentz.SetPtEtaPhiM(nu2cand->pt(), nu2cand->eta(), nu2cand->phi(), 0); 
        bbar_lorentz.SetXYZT(b1cand->px()+b2cand->px(), b1cand->py()+b2cand->py(), b1cand->pz()+b2cand->pz(), b1cand->energy()+b2cand->energy());
        int onshellMarker = -1;
        if (w1cand->mass() > w2cand->mass()) onshellMarker = 1;
        else onshellMarker = 2;
         std::cout <<" mu1 lorenz "; mu1_lorentz.Print(); 
        //thismmc = new MMC();
        //std::cout << "onshellMarkder  " << onshellMarker << std::endl;
	thismmc = new MMC(&mu1_lorentz, &mu2_lorentz, &bbar_lorentz, &met_lorentz, &nu1_lorentz, &nu2_lorentz, onshellMarker, 
	simulation_, ievent, mmcset_, fs, verbose_);
        //thismmc->printTrueLorentz();
        thismmc->runMMC();	
        delete thismmc;

    }   */ 
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
   evtree->Branch("mu1_energy",&mu1_energy);
   evtree->Branch("mu1_px",&mu1_px);
   evtree->Branch("mu1_py",&mu1_py);
   evtree->Branch("mu1_pz",&mu1_pz);
   evtree->Branch("mu1_eta",&mu1_eta);
   evtree->Branch("mu1_phi",&mu1_phi);
   evtree->Branch("w1_energy",&w1_energy);
   evtree->Branch("w1_px",&w1_px);
   evtree->Branch("w1_py",&w1_py);
   evtree->Branch("w1_pz",&w1_pz);
   evtree->Branch("w1_mass",&w1_mass);
   evtree->Branch("nu1_energy",&nu1_energy);
   evtree->Branch("nu1_px",&nu1_px);
   evtree->Branch("nu1_py",&nu1_py);
   evtree->Branch("nu1_pz",&nu1_pz);
   evtree->Branch("nu1_eta",&nu1_eta);
   evtree->Branch("nu1_phi",&nu1_phi);
   evtree->Branch("Wtomu1nu1",&Wtomu1nu1);

   evtree->Branch("mu2_energy",&mu2_energy);
   evtree->Branch("mu2_px",&mu2_px);
   evtree->Branch("mu2_py",&mu2_py);
   evtree->Branch("mu2_pz",&mu2_pz);
   evtree->Branch("mu2_eta",&mu2_eta);
   evtree->Branch("mu2_phi",&mu2_phi);
   evtree->Branch("w2_energy",&w2_energy);
   evtree->Branch("w2_px",&w2_px);
   evtree->Branch("w2_py",&w2_py);
   evtree->Branch("w2_pz",&w2_pz);
   evtree->Branch("w2_mass",&w2_mass);
   evtree->Branch("nu2_energy",&nu2_energy);
   evtree->Branch("nu2_px",&nu2_px);
   evtree->Branch("nu2_py",&nu2_py);
   evtree->Branch("nu2_pz",&nu2_pz);
   evtree->Branch("nu2_eta",&nu2_eta);
   evtree->Branch("nu2_phi",&nu2_phi);
   evtree->Branch("Wtomu2nu2",&Wtomu2nu2);

   evtree->Branch("htoWW_energy",&htoWW_energy);
   evtree->Branch("htoWW_px",&htoWW_px);
   evtree->Branch("htoWW_py",&htoWW_py);
   evtree->Branch("htoWW_pz",&htoWW_pz);
   evtree->Branch("htoWW_mass",&htoWW_mass);
   
   evtree->Branch("b1_energy",&b1_energy);
   evtree->Branch("b1_px",&b1_px);
   evtree->Branch("b1_py",&b1_py);
   evtree->Branch("b1_pz",&b1_pz);
   evtree->Branch("b1_motherid",&b1_motherid);
   evtree->Branch("b2_energy",&b2_energy);
   evtree->Branch("b2_px",&b2_px);
   evtree->Branch("b2_py",&b2_py);
   evtree->Branch("b2_pz",&b2_pz);
   evtree->Branch("b2_motherid",&b2_motherid);
   evtree->Branch("bjet_energy",&bjet_energy);
   evtree->Branch("bjet_px",&bjet_px);
   evtree->Branch("bjet_py",&bjet_py);
   evtree->Branch("bjet_pz",&bjet_pz);
   evtree->Branch("bjet_mass",&bjet_mass);
   evtree->Branch("bbarjet_energy",&bbarjet_energy);
   evtree->Branch("bbarjet_px",&bbarjet_px);
   evtree->Branch("bbarjet_py",&bbarjet_py);
   evtree->Branch("bbarjet_pz",&bbarjet_pz);
   evtree->Branch("bbarjet_mass",&bbarjet_mass);
   
   evtree->Branch("htobb_energy",&htobb_energy);
   evtree->Branch("htobb_px",&htobb_px);
   evtree->Branch("htobb_py",&htobb_py);
   evtree->Branch("htobb_pz",&htobb_pz);
   evtree->Branch("htobb_mass",&htobb_mass);
   
   evtree->Branch("h2tohh_energy",&h2tohh_energy);
   evtree->Branch("h2tohh_px",&h2tohh_px);
   evtree->Branch("h2tohh_py",&h2tohh_py);
   evtree->Branch("h2tohh_pz",&h2tohh_pz);
   evtree->Branch("h2tohh_mass",&h2tohh_mass);

   evtree->Branch("met",&met);
   evtree->Branch("met_phi",&met_phi);
   evtree->Branch("met_px",&met_px);
   evtree->Branch("met_py",&met_py);
   
   evtree->Branch("htobb",&htobb);
   evtree->Branch("htoWW",&htoWW);
   evtree->Branch("h2tohh",&findAllGenParticles);
    
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
DiHiggsWWBBAnalyzer::checkGenParticlesSingal(edm::Handle<reco::GenParticleCollection> genParticleColl){

    std::vector<reco::GenParticle*> b1Coll; 
    std::vector<reco::GenParticle*> b2Coll;
    std::vector<reco::GenParticle*> W1Coll;
    std::vector<reco::GenParticle*> W2Coll;
    std::vector<const reco::Candidate*> htoWWColl;
    std::vector<const reco::Candidate*> htoBBColl;

    std::cout <<"*********** start to check GenParticles for Singal sample ***********"<< std::endl;
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
		//fillbranches();
        //evtree->Fill();
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
    const reco::Candidate* t1cand;
    const reco::Candidate* t2cand;

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
	//fillbranches();
        //evtree->Fill();
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
      mu1_energy = mu1cand->energy();
      mu1_px = mu1cand->px();
      mu1_py = mu1cand->py();
      mu1_pz = mu1cand->pz();
      w1_energy = w1cand->energy();
      w1_px = w1cand->px();
      w1_py = w1cand->py();
      w1_pz = w1cand->pz();
      w1_mass = w1cand->mass();
      nu1_energy = nu1cand->energy();
      nu1_px = nu1cand->px();
      nu1_py = nu1cand->py();
      nu1_pz = nu1cand->pz();

      mu1_eta = mu1cand->eta();
      mu1_phi = mu1cand->phi();
      nu1_eta = nu1cand->eta();
      nu1_phi = nu1cand->phi();

      mu2_energy = mu2cand->energy();
      mu2_px = mu2cand->px();
      mu2_py = mu2cand->py();
      mu2_pz = mu2cand->pz();
      w2_energy = w2cand->energy();
      w2_px = w2cand->px();
      w2_py = w2cand->py();
      w2_pz = w2cand->pz();
      w2_mass = w2cand->mass();
      nu2_energy = nu2cand->energy();
      nu2_px = nu2cand->px();
      nu2_py = nu2cand->py();
      nu2_pz = nu2cand->pz();

      mu2_eta = mu2cand->eta();
      mu2_phi = mu2cand->phi();
      nu2_eta = nu2cand->eta();
      nu2_phi = nu2cand->phi();

      htoWW_energy = htoWWcand->energy();
      htoWW_px = htoWWcand->px();
      htoWW_py = htoWWcand->py();
      htoWW_pz = htoWWcand->pz();
      htoWW_mass = htoWWcand->mass();

      b1_energy = b1cand->energy();
      b1_px = b1cand->px();
      b1_py = b1cand->py();
      b1_pz = b1cand->pz();
      b2_energy = b2cand->energy();
      b2_px = b2cand->px();
      b2_py = b2cand->py();
      b2_pz = b2cand->pz();

      bjet_energy = bjet_lorentz.Energy();
      bjet_px = bjet_lorentz.Px();
      bjet_py = bjet_lorentz.Py();
      bjet_pz = bjet_lorentz.Pz();
      bjet_mass = bjet_lorentz.M();
      bbarjet_energy = bbarjet_lorentz.Energy();
      bbarjet_px = bbarjet_lorentz.Px();
      bbarjet_py = bbarjet_lorentz.Py();
      bbarjet_pz = bbarjet_lorentz.Pz();
      bbarjet_mass = bbarjet_lorentz.M();
      
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
   
      met = met_lorentz.Energy();
      met_phi = met_lorentz.Phi();
      met_px = met_lorentz.Px();
      met_py = met_lorentz.Py();
   
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

}



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


