#include "HhhAnalysis/CutFlowAnalyzer/interface/POGRecipesRun2.h"
//===================================
//Muon POG twiki: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Medium_Muon
//Electron POG: 
//Jet POG:
//MET
//===================================


bool POGRecipesRun2::is2016MediumMuon(const pat::Muon & recoMu) 
{
    bool goodGlob = recoMu.isGlobalMuon() && 
		     recoMu.globalTrack()->normalizedChi2() < 3 && 
		     recoMu.combinedQuality().chi2LocalPosition < 12 && 
		     recoMu.combinedQuality().trkKink < 20; 
    bool isMedium = muon::isLooseMuon(recoMu) && 
		   recoMu.innerTrack()->validFraction() > 0.49 && 
		   muon::segmentCompatibility(recoMu) > (goodGlob ? 0.303 : 0.451); 
    return isMedium; 
}


bool POGRecipesRun2::isMediumMuon(const pat::Muon & recoMu) 
{
    bool goodGlob = recoMu.isGlobalMuon() && 
		     recoMu.globalTrack()->normalizedChi2() < 3 && 
		     recoMu.combinedQuality().chi2LocalPosition < 12 && 
		     recoMu.combinedQuality().trkKink < 20; 
    bool isMedium = muon::isLooseMuon(recoMu) && 
		   recoMu.innerTrack()->validFraction() > 0.8 && 
		   muon::segmentCompatibility(recoMu) > (goodGlob ? 0.303 : 0.451); 
    return isMedium; 
}


float POGRecipesRun2::MuonIsoPFbased(const pat::Muon & recoMu)
{
    return (recoMu.pfIsolationR04().sumChargedHadronPt + std::max(0., recoMu.pfIsolationR04().sumNeutralHadronEt + recoMu.pfIsolationR04().sumPhotonEt - 0.5*recoMu.pfIsolationR04().sumPUPt))/recoMu.pt();
}


float POGRecipesRun2::MuonIsoTrackerbased(const pat::Muon & recoMu)
{
    return recoMu.isolationR03().sumPt/recoMu.pt();
}


float POGRecipesRun2::getMuonTriggerSF(float mueta, float mupt, std::string filename, std::string histname)
{
     TFile* file = new TFile(filename.c_str());
     TH2F* abseta_pt_ratio = (TH2F*)file->Get(histname.c_str());
     int bin1 = abseta_pt_ratio->GetXaxis()->FindBin(mueta);
     int bin2 = abseta_pt_ratio->GetYaxis()->FindBin(mupt);
     if (bin1==0 or bin1==abseta_pt_ratio->GetNbinsX()+1 or bin2==0 or bin2==abseta_pt_ratio->GetNbinsY()+1)
	 return 1.0;//not find corresponding bin
     float sf = abseta_pt_ratio->GetBinContent(bin1, bin2);
     delete file;
     return sf;
}


float POGRecipesRun2::getMuonISOSF(float mueta, float mupt, std::string filename, std::string histname)
{
     TFile* file = new TFile(filename.c_str());
     TH2F* abseta_pt_ratio = (TH2F*)file->Get(histname.c_str());
     int bin1 = abseta_pt_ratio->GetXaxis()->FindBin(mueta);
     int bin2 = abseta_pt_ratio->GetYaxis()->FindBin(mupt);
     if (bin1==0 or bin1==abseta_pt_ratio->GetNbinsX()+1 or bin2==0 or bin2==abseta_pt_ratio->GetNbinsY()+1)
	 return 1.0;//not find corresponding bin
     float sf = abseta_pt_ratio->GetBinContent(bin1, bin2);
     delete file;
     return sf;
}

float POGRecipesRun2::getMuonIDSF(float mueta, float mupt, std::string filename, std::string histname)
{
     TFile* file = new TFile(filename.c_str());
     TH2F* abseta_pt_ratio = (TH2F*)file->Get(histname.c_str());
     int bin1 = abseta_pt_ratio->GetXaxis()->FindBin(mueta);
     int bin2 = abseta_pt_ratio->GetYaxis()->FindBin(mupt);
     if (bin1==0 or bin1==abseta_pt_ratio->GetNbinsX()+1 or bin2==0 or bin2==abseta_pt_ratio->GetNbinsY()+1)
	 return 1.0;//not find corresponding bin
     float sf = abseta_pt_ratio->GetBinContent(bin1, bin2);
     delete file;
     return sf;
}

float POGRecipesRun2::getMuonTrackingSF(float mueta, std::string filename, std::string histname)
{
     TFile* file = new TFile(filename.c_str());
     TGraphAsymmErrors* ratio_eff_eta3_dr030e030_corr = (TGraphAsymmErrors*)file->Get(histname.c_str());
     int n = ratio_eff_eta3_dr030e030_corr->GetN();
     double eta_up = 0.0;
     double sf_up = 0.0;
     double eta_low = 0.0;
     double sf_low = 0.0;
     for (int i =0; i < n-1; i++ ){
	 ratio_eff_eta3_dr030e030_corr->GetPoint(i, eta_low, sf_low);
	 ratio_eff_eta3_dr030e030_corr->GetPoint(i+1, eta_up, sf_up);
	 if (float(eta_low) <= mueta and mueta < float(eta_up))
	     break;
     }
     if (mueta > float(eta_up)) return 1.0;//not found eta bin
     float sf = sf_low + (sf_up-sf_low)/(eta_up-eta_low)*(mueta-eta_low);
     delete file;
     return sf;
}

POGRecipesRun2::MuonPOGSFManager::MuonPOGSFManager(std::vector<std::string> descriptions, std::vector<std::string> files, std::vector<std::string> histnames)
{
    int i=0;
    for (auto type : descriptions){
	if (type == "Trigger"){
	    filename_trigger = files[i];
	    histname_trigger = histnames[i];
	    trigger_valid  = true;
	}else if (type == "ISO"){
	    filename_iso = files[i];
	    histname_iso = histnames[i];
	    iso_valid  = true;
	}else if (type == "ID"){
	    filename_id = files[i];
	    histname_id = histnames[i];
	    id_valid  = true;
	}else if (type == "Tracking"){
	    filename_tracking = files[i];
	    histname_tracking = histnames[i];
	    tracking_valid  = true;
	}
	++i;
    }
    if (trigger_valid){
	file_trigger = new TFile(filename_trigger.c_str());
	hist_trigger = (TH2F*)file_trigger->Get(histname_trigger.c_str());
    }
    if (iso_valid){
	file_iso = new TFile(filename_iso.c_str());
	hist_iso = (TH2F*)file_iso->Get(histname_iso.c_str());
    }
    if (id_valid){
	file_id = new TFile(filename_id.c_str());
	hist_id = (TH2F*)file_id->Get(histname_id.c_str());
    }
    if (tracking_valid){
	file_tracking = new TFile(filename_tracking.c_str());
	graph_tracking = (TGraphAsymmErrors*)file_tracking->Get(histname_tracking.c_str());
    }
}

POGRecipesRun2::MuonPOGSFManager::~MuonPOGSFManager()
{
    if (trigger_valid) delete file_trigger;
    if (iso_valid) delete file_iso;
    if (id_valid) delete file_id;
    if (tracking_valid) delete file_tracking;
}


	
float POGRecipesRun2::MuonPOGSFManager::getMuonTriggerMCSF(float mueta, float mupt)
{
    if (not trigger_valid) return 1.0;
     int bin1 = hist_trigger->GetXaxis()->FindBin(mueta);
     int bin2 = hist_trigger->GetYaxis()->FindBin(mupt);
     if (bin1==0 or bin1==hist_trigger->GetNbinsX()+1 or bin2==0 or bin2==hist_trigger->GetNbinsY()+1)
	 return 1.0;//not find corresponding bin
     float sf = hist_trigger->GetBinContent(bin1, bin2);
     return sf;

}

float POGRecipesRun2::MuonPOGSFManager::getMuonISOMCSF(float mueta, float mupt)
{
    if (not iso_valid) return 1.0;
     int bin1 = hist_iso->GetXaxis()->FindBin(mueta);
     int bin2 = hist_iso->GetYaxis()->FindBin(mupt);
     if (bin1==0 or bin1==hist_iso->GetNbinsX()+1 or bin2==0 or bin2==hist_iso->GetNbinsY()+1)
	 return 1.0;//not find corresponding bin
     float sf = hist_iso->GetBinContent(bin1, bin2);
     return sf;

}

float POGRecipesRun2::MuonPOGSFManager::getMuonIDMCSF(float mueta, float mupt)
{
    if (not id_valid) return 1.0;
     int bin1 = hist_id->GetXaxis()->FindBin(mueta);
     int bin2 = hist_id->GetYaxis()->FindBin(mupt);
     if (bin1==0 or bin1==hist_id->GetNbinsX()+1 or bin2==0 or bin2==hist_id->GetNbinsY()+1)
	 return 1.0;//not find corresponding bin
     float sf = hist_id->GetBinContent(bin1, bin2);
     return sf;

}

float POGRecipesRun2::MuonPOGSFManager::getMuonTrackingMCSF(float mueta)
{
    if (not tracking_valid) return 1.0;
     int n = graph_tracking->GetN();
     double eta_up = 0.0;
     double sf_up = 0.0;
     double eta_low = 0.0;
     double sf_low = 0.0;
     for (int i =0; i < n-1; i++ ){
	 graph_tracking->GetPoint(i, eta_low, sf_low);
	 graph_tracking->GetPoint(i+1, eta_up, sf_up);
	 if (float(eta_low) <= mueta and mueta < float(eta_up))
	     break;
     }
     if (mueta > float(eta_up)) return 1.0;//not found eta bin
     float sf = sf_low + (sf_up-sf_low)/(eta_up-eta_low)*(mueta-eta_low);
     return sf;
}
