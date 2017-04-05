#ifndef HhhAnalysis_CutFlowAnalyzer_POGRecipesRun2
#define HhhAnalysis_CutFlowAnalyzer_POGRecipesRun2


//#include "DataFormats/JetReco/interface/PFJet.h"
//#include "DataFormats/PatCandidates/interface/Jet.h"
//#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
//standard selector
//#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "TH2F.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"


namespace POGRecipesRun2{

    bool isMediumMuon(const pat::Muon & recoMu);//Run2016 GH data
    bool is2016MediumMuon(const pat::Muon & recoMu);//Run2016 BCDEF data
    float MuonIsoPFbased(const pat::Muon & recoMu);
    float MuonIsoTrackerbased(const pat::Muon & recoMu);
    float getMuonTriggerSF(float mueta, float mupt, std::string filename, std::string histname);
    float getMuonISOSF(float mueta, float mupt, std::string filename, std::string histname);
    float getMuonIDSF(float mueta, float mupt, std::string filename, std::string histname);
    float getMuonTrackingSF(float mueta, std::string filename, std::string histname);

}


#endif
