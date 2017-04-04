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



