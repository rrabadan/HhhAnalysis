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
     if (bin1==0 or bin1==abseta_pt_ratio->GetNbinsX() or bin2==0 or bin2==abseta_pt_ratio->GetNbinsY())
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
     if (bin1==0 or bin1==abseta_pt_ratio->GetNbinsX() or bin2==0 or bin2==abseta_pt_ratio->GetNbinsY())
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
     if (bin1==0 or bin1==abseta_pt_ratio->GetNbinsX() or bin2==0 or bin2==abseta_pt_ratio->GetNbinsY())
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

