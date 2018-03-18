
def electronImpactParameterCut( electron):
    """ check electron impact parameter, https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2 """
    dxy_endcap = 0.1 #cm
    dxy_barrel = 0.05
    dz_endcap  = 0.2
    dz_barrel = 0.1
    if abs(electron.eta) < 1.479: ##barrel
	return (abs(electron.dz) < dz_barrel and abs(electron.dxy) < dxy_barrel)
    elif abs(electron.eta) >= 1.479 and abs(electron.eta) < 2.5:
	return (abs(electron.dz) < dz_endcap and abs(electron.dxy) < dxy_endcap)
    else:
	return False

def electronID( electron):
    """ check electron ID, https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2 """
    ### we use medium Electron
    return electron.cutBased >= 3
    
def electronIso( electron):
    """ check electron isolation, https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2 """
    ### electron iso cut
    return True

def electronHLTSafeID( electron):
    """ check electron isolation, https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2 """
    """https://twiki.cern.ch/twiki/bin/view/CMS/ChangesEGMHLTAlgo2014"""
    ### extrac HLT safe ID cut is applied for electron selection for 2016 data only
    return True
    #if abs(electron.eta) < 1.479: #barrel
    #return (abs(electron.sieie) < 0.11 and 
    
def muonID( muon):
    """https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#Muon_Isolation """
    #goodGlob = (muon.isGlobalMuon() and muon.globalTrack()->normalizedChi2() < 3 and muon.combinedQuality().chi2LocalPosition < 12 and muon.combinedQuality().trkKink < 20)
    ###check Muon ID
    return muon.mediumId

def muonIso(muon):
    ###check Muon iso 
    return muon.pfRelIso03_all < 0.15

def jetMediumBtagging(jet):
    """https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80X """ 
    ### use medium btagging 
    return jet.btagCMVA >= 0.185
