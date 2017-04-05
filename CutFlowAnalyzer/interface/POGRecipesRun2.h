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
#include <sstream>
#include <string>
#include <vector>


namespace POGRecipesRun2{

    bool isMediumMuon(const pat::Muon & recoMu);//Run2016 GH data
    bool is2016MediumMuon(const pat::Muon & recoMu);//Run2016 BCDEF data
    float MuonIsoPFbased(const pat::Muon & recoMu);
    float MuonIsoTrackerbased(const pat::Muon & recoMu);
    float getMuonTriggerSF(float mueta, float mupt, std::string filename, std::string histname);
    float getMuonISOSF(float mueta, float mupt, std::string filename, std::string histname);
    float getMuonIDSF(float mueta, float mupt, std::string filename, std::string histname);
    float getMuonTrackingSF(float mueta, std::string filename, std::string histname);
  class MuonPOGSFManager {
    public:
	MuonPOGSFManager(std::vector<std::string> descriptions, std::vector<std::string> files, std::vector<std::string> histnames);
	~MuonPOGSFManager();
	float getMuonTriggerMCSF(float mueta, float mupt);
	float getMuonISOMCSF(float mueta, float mupt);
	float getMuonIDMCSF(float mueta, float mupt);
	float getMuonTrackingMCSF(float mueta);
    private:
	std::string filename_trigger;
	std::string filename_iso;
	std::string filename_id;
	std::string filename_tracking;
	std::string histname_trigger;
	std::string histname_iso;
	std::string histname_id;
	std::string histname_tracking;
	TH2F* hist_trigger;
	TH2F* hist_iso;
	TH2F* hist_id;
	TGraphAsymmErrors* graph_tracking;
	TFile* file_trigger;
	TFile* file_iso;
	TFile* file_id;
	TFile* file_tracking;
	bool trigger_valid;
	bool iso_valid;
	bool id_valid;
	bool tracking_valid;
  };

}


#endif
