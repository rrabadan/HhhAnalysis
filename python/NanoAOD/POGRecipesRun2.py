import ROOT


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
    return electron.pfRelIso03_all < .04

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
    ###check Muon iso , tight isolation: 0.15; loose iso: 0.25 
    return muon.pfRelIso04_all < 0.15

def leptonImpactParameter(lepton):
    return ((abs(lepton.pdgId) == 11 and electronImpactParameterCut(lepton)) or abs(lepton.pdgId) == 13 )
def leptonpairImpactParameter(leptonpair):
    return (leptonImpactParameter(leptonpair[0]) and leptonImpactParameter(leptonpair[1]))

def leptonIso(lepton):
    #print "lepton id ",lepton.pdgId
    return ((abs(lepton.pdgId) == 11 and electronIso(lepton)) or (abs(lepton.pdgId) == 13 and muonIso(lepton)))
def leptonpairIso(leptonpair):
    return (leptonIso(leptonpair[0]) and leptonIso(leptonpair[1]))

def leptonID(lepton):
    return ((abs(lepton.pdgId) == 11 and electronID(lepton)) or (abs(lepton.pdgId) == 13 and muonID(lepton)))
def leptonpairID(leptonpair):
    return (leptonID(leptonpair[0]) and leptonID(leptonpair[1]))

def leptonHLTSafeID(lepton):
    return ((abs(lepton.pdgId) == 11 and electronHLTSafeID(lepton)) or (abs(lepton.pdgId) == 13))
def leptonpairHLTSafeID(leptonpair):
    return (leptonHLTSafeID(leptonpair[0]) and leptonHLTSafeID(leptonpair[1]))

def jetMediumBtagging(jet):
    """https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80X """ 
    ### use medium btagging 
    return jet.btagCMVA >= 0.185



class LeptonSFManager():
    def __init__(self):
	self.EGSF_filename  = "leptonSF/EGM2D_eleGSF.root"
	self.MuonIDSF_filename = "leptonSF/Mu_ID.root"
	self.MuonIsoSF_filename = "leptonSF/Mu_Iso.root"
	self.MuonTrgSF_filename = "leptonSF/Mu_Trg.root"
	
	### x-axis: eta,  y-axis: pt
	self.EGSF_histname = "EGamma_SF2D"
	self.MuonIDSF_histname = "MC_NUM_MediumID2016_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio"
	self.MuonIsoSF_histname = "TightISO_MediumID_pt_eta/abseta_pt_ratio"
	self.MuonTrgSF_histname = "IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio"

	self.EGSF_tfile = ROOT.TFile(self.EGSF_filename,"READ")
	self.MuonIDSF_tfile = ROOT.TFile(self.MuonIDSF_filename,"READ")
	self.MuonIsoSF_tfile = ROOT.TFile(self.MuonIsoSF_filename,"READ")
	self.MuonTrgSF_tfile = ROOT.TFile(self.MuonTrgSF_filename,"READ")

	self.EGSF_th2 = self.EGSF_tfile.Get(self.EGSF_histname)
        #print "self.EGSF_th2 ",self.EGSF_th2.Print("ALL")
	self.MuonIDSF_th2 = self.MuonIDSF_tfile.Get(self.MuonIDSF_histname)
	self.MuonIsoSF_th2 = self.MuonIsoSF_tfile.Get(self.MuonIsoSF_histname)
	self.MuonTrgSF_th2 = self.MuonTrgSF_tfile.Get(self.MuonTrgSF_histname)


    #### FIXME, add uncertainty in next version
    def getSF(self,  th2, eta, pt):
	 bin1 = th2.GetXaxis().FindBin(eta);
	 bin2 = th2.GetYaxis().FindBin(pt);
	 #if th2.GetName() == "EGamma_SF2D":
	 #    print "EGamma_SF2D bin1 ".bin1," bin2 ",bin2
	 if (bin1==0 or bin1== th2.GetNbinsX()+1 or bin2==0 or bin2== th2.GetNbinsY()+1):
	     return 1.0,0.0,0.0
	 return th2.GetBinContent(bin1, bin2),th2.GetBinErrorUp(bin1, bin2),th2.GetBinErrorLow(bin1, bin2)

    def getEGSF(self, eta, pt):##final one ?
	SF = self.getSF(self.EGSF_th2, eta, pt)
        #print "Electron SFs ",SF
	return SF

    def getMuonIDSF(self, eta, pt):
	return self.getSF(self.MuonIDSF_th2, abs(eta), pt)

    def getMuonIsoSF(self, eta, pt):
	return self.getSF(self.MuonIsoSF_th2, abs(eta), pt)

    def getMuonTrgSF(self, eta, pt):
	return self.getSF(self.MuonTrgSF_th2, abs(eta), pt)

    def getleptonTrgSF(self, lep):
	if abs(lep.pdgId) == 11:
	    return 1.0,0.0,0.0
	elif abs(lep.pdgId) == 13:
	    return self.getMuonTrgSF(lep.eta, lep.pt)
	else:
	    raise ValueError('Getting lepton SFs, not electron nor muon ', lep.pdgId)

    def getleptonpairTrgSF(self, leptonpair):
	SF1 = self.getleptonTrgSF(leptonpair[0])
       	SF2 = self.getleptonTrgSF(leptonpair[1])
        return SF1[0]*SF2[0],  SF1[1]*SF2[1], SF1[2]*SF2[2]

    def getleptonIsoSF(self, lep):
	if abs(lep.pdgId) == 11:
	    return self.getEGSF(lep.eta, lep.pt)
	elif abs(lep.pdgId) == 13:
	    return self.getMuonIsoSF(lep.eta, lep.pt)
	else:
	    raise ValueError('Getting lepton SFs, not electron nor muon ', lep.pdgId)

    def getleptonpairIsoSF(self, leptonpair):
	SF1 = self.getleptonIsoSF(leptonpair[0])
       	SF2 = self.getleptonIsoSF(leptonpair[1])
        return SF1[0]*SF2[0],  SF1[1]*SF2[1], SF1[2]*SF2[2]

    def getleptonIDSF(self, lep):
	if abs(lep.pdgId) == 11:
	    return 1.0,0.0,0.0
	elif abs(lep.pdgId) == 13:
	    return self.getMuonIDSF(lep.eta, lep.pt)
	else:
	    raise ValueError('Getting lepton SFs, not electron nor muon ', lep.pdgId)

    def getleptonpairIDSF(self, leptonpair):
	SF1 = self.getleptonIDSF(leptonpair[0])
       	SF2 = self.getleptonIDSF(leptonpair[1])
        return SF1[0]*SF2[0],  SF1[1]*SF2[1], SF1[2]*SF2[2]




	  
