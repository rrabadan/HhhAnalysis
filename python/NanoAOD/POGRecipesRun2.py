import ROOT
import string
from math import sqrt
import os


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
    return electron.pfRelIso03_all < 0.04

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
    """https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco"""
    ### use medium btagging 
    ### SF should be from  leptonSF/cMVAv2_Moriond17_B_H.csv 
    return (jet.btagCMVA > 0.4432)

def combinedError(err1, err2, weight1):
    return sqrt(err1*err1*weight1+err2*err2*(1.0-weight1));

class LeptonSFManager():
    ### to take the lumi into consideration?
    ### Electron Triggering, IP, Isolation, tracking?
    ### Muon tracking ?
    def __init__(self):
	##brilcalc lumi -u /pb  --normtag normtag_PHYSICS.json -i json.txt
	## used delivered lumi for normalization

	##self.Lumi_BCDEF = 5.750+2.573+4.242+4.025+3.105 ## recorded,23Sep2016ReReco
	self.Lumi_BCDEF = 5.991+2.685+4.411+4.222+3.303 ##delievered
	##self.Lumi_GH = 7.576+8.651 ## recorded, 23Sep2016ReReco
	self.Lumi_GH = 7.865+8.985
	self.totalLumi = self.Lumi_BCDEF  + self.Lumi_GH#fb-1
	HhhPath = "/home/taohuang/DiHiggsAnalysis/CMSSW_9_4_0_pre1/src/HhhAnalysis/python/NanoAOD/"
	#self.EGSF_filename  = "leptonSF/EGM2D_eleGSF.root" ## for electron ID
	self.EGSF_filename  =             os.path.join(HhhPath , "leptonSF/egammaEffi.txt_EGM2D.root" )## for electron ID
	self.MuonIDSF_filename =          os.path.join(HhhPath , "leptonSF/EfficienciesAndSF_BCDEF_ID.root")
	self.MuonIsoSF_filename =         os.path.join(HhhPath , "leptonSF/EfficienciesAndSF_BCDEF_ISO.root")
	self.MuonTrgSF_filename =         os.path.join(HhhPath , "leptonSF/EfficienciesAndSF_BCDEF_trigger.root")
	self.MuonTrackingSF_filename =    os.path.join(HhhPath , "leptonSF/EfficienciesAndSF_BCDEF_Tracking.root")
	self.MuonIDSF_GH_filename =       os.path.join(HhhPath , "leptonSF/EfficienciesAndSF_GH_ID.root")
	self.MuonIsoSF_GH_filename =      os.path.join(HhhPath , "leptonSF/EfficienciesAndSF_GH_ISO.root")
	self.MuonTrgSF_GH_filename =      os.path.join(HhhPath , "leptonSF/EfficienciesAndSF_GH_trigger.root")
	self.MuonTrackingSF_GH_filename = os.path.join(HhhPath , "leptonSF/EfficienciesAndSF_GH_Tracking.root")

	
	### x-axis: eta,  y-axis: pt
	self.EGSF_histname = "EGamma_SF2D"
	self.MuonIDSF_histname = "MC_NUM_MediumID2016_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio"
	self.MuonIsoSF_histname = "TightISO_MediumID_pt_eta/abseta_pt_ratio"
	self.MuonTrgSF_histname = "IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio"
	self.MuonTrackingSF_tgraphname = "ratio_eff_eta3_dr030e030_corr"

	self.EGSF_tfile = ROOT.TFile(self.EGSF_filename,"READ")
	self.MuonIDSF_tfile = ROOT.TFile(self.MuonIDSF_filename,"READ")
	self.MuonIsoSF_tfile = ROOT.TFile(self.MuonIsoSF_filename,"READ")
	#self.MuonTrgSF_tfile = ROOT.TFile(self.MuonTrgSF_filename,"READ")
	self.MuonTrackingSF_tfile = ROOT.TFile(self.MuonTrackingSF_filename,"READ")
	self.MuonIDSF_GH_tfile = ROOT.TFile(self.MuonIDSF_GH_filename,"READ")
	self.MuonIsoSF_GH_tfile = ROOT.TFile(self.MuonIsoSF_GH_filename,"READ")
	#self.MuonTrgSF_GH_tfile = ROOT.TFile(self.MuonTrgSF_GH_filename,"READ")
	self.MuonTrackingSF_GH_tfile = ROOT.TFile(self.MuonTrackingSF_GH_filename,"READ")

	self.EGSF_th2 = self.EGSF_tfile.Get(self.EGSF_histname)
        #print "self.EGSF_th2 ",self.EGSF_th2.Print("ALL")
	self.MuonIDSF_th2 = self.MuonIDSF_tfile.Get(self.MuonIDSF_histname)
	self.MuonIsoSF_th2 = self.MuonIsoSF_tfile.Get(self.MuonIsoSF_histname)
	#self.MuonTrgSF_th2 = self.MuonTrgSF_tfile.Get(self.MuonTrgSF_histname)
	self.MuonTrackingSF_tgraph= self.MuonTrackingSF_tfile.Get(self.MuonTrackingSF_tgraphname)
        ##### GH
	self.MuonIDSF_GH_th2 = self.MuonIDSF_GH_tfile.Get(self.MuonIDSF_histname)
	self.MuonIsoSF_GH_th2 = self.MuonIsoSF_GH_tfile.Get(self.MuonIsoSF_histname)
	#self.MuonTrgSF_GH_th2 = self.MuonTrgSF_GH_tfile.Get(self.MuonTrgSF_histname)
	self.MuonTrackingSF_GH_tgraph= self.MuonTrackingSF_GH_tfile.Get(self.MuonTrackingSF_tgraphname)

    	self.MuonTrackingSF_nbins = self.MuonTrackingSF_tgraph.GetN()
    	eta = ROOT.Double(0.0); trackingSF =  ROOT.Double(0.0);  trackingSF_GH =  ROOT.Double(0.0);
	self.MuonTrackingSF_allbins = []
	self.MuonTrackingSF_GH_lumiratio = self.Lumi_GH/self.totalLumi
	## how to weight SF based on Lumi ?
	for i in range(0, self.MuonTrackingSF_nbins):
	    self.MuonTrackingSF_tgraph.GetPoint(i, eta, trackingSF)
	    self.MuonTrackingSF_GH_tgraph.GetPoint(i, eta, trackingSF_GH)
	    thisbin = {}
	    xlow = self.MuonTrackingSF_tgraph.GetErrorXlow(i)
	    xhigh = self.MuonTrackingSF_tgraph.GetErrorXhigh(i)
	    ylow = self.MuonTrackingSF_tgraph.GetErrorYlow(i)
	    yhigh = self.MuonTrackingSF_tgraph.GetErrorYhigh(i)
	    ylow_GH = self.MuonTrackingSF_GH_tgraph.GetErrorYlow(i)
	    yhigh_GH = self.MuonTrackingSF_GH_tgraph.GetErrorYhigh(i)
    	    thisbin["etalow"]  = eta - xlow
    	    thisbin["etahigh"]  = eta + xhigh
	    thisbin["SF"] = trackingSF*(1-self.MuonTrackingSF_GH_lumiratio) +  trackingSF_GH*self.MuonTrackingSF_GH_lumiratio
	    
    	    thisbin["SFerrlow"]  = (-1.0)*combinedError(ylow_GH, ylow, self.MuonTrackingSF_GH_lumiratio) + thisbin["SF"]
    	    thisbin["SFerrhigh"]  = combinedError(yhigh_GH, yhigh, self.MuonTrackingSF_GH_lumiratio) + thisbin["SF"]
	    self.MuonTrackingSF_allbins.append(thisbin)


    	self.legs = ["DoubleEleLegHigPt","DoubleEleLegLowPt","DoubleMuLegHigPt","DoubleMuLegLowPt","MuEleLegHigPt", "MuEleLegLowPt","EleMuLegHigPt","EleMuLegLowPt"]
	self.TriggerSFs = {}
	for leg in self.legs:
	    filepath = os.path.join(HhhPath, "leptonSF/HLT_"+leg+".text")
	    trigfile = open(filepath,"r")
	    ##in string
	    trigdata = [line.rstrip().split() for line in trigfile  if '#' not in line]
	    self.TriggerSFs[leg] = [[float(x) for x in bins ] for bins in trigdata]
	    #print "leg  ",leg,"\n",self.TriggerSFs[leg]
	    
	##from HHbbWW analysis note: AN_HIG-17-006
	self.DZEffs = {"DoubleEle": 0.983, "DoubleMu":0.993, "MuEle":0.988, "EleMu":0.982}	

    #### FIXME, add uncertainty in next version
    def getSF(self,  th2, eta, pt):
	 bin1 = th2.GetXaxis().FindBin(eta);
	 bin2 = th2.GetYaxis().FindBin(pt);
	 #if th2.GetName() == "EGamma_SF2D":
	 #    print "EGamma_SF2D bin1 ".bin1," bin2 ",bin2
	 if (bin1==0 or bin1== th2.GetNbinsX()+1 or bin2==0 or bin2== th2.GetNbinsY()+1):
	     return 1.0,1.0,1.0
	 errup = th2.GetBinErrorUp(bin1, bin2)
	 errlow = th2.GetBinErrorLow(bin1, bin2)
	 SF = th2.GetBinContent(bin1, bin2)
	 return SF, SF+errup, SF-errlow

    def getEGSF(self, eta, pt):##final one ?
	SF = self.getSF(self.EGSF_th2, eta, pt)
        #print "Electron SFs ",SF
	return SF

    def getMuonIDSF(self, eta, pt):
	SF_BCDEF, errlow_BCDEF, errhigh_BCDEF = self.getSF(self.MuonIDSF_th2, abs(eta), pt)
	SF_GH, errlow_GH, errhigh_GH = self.getSF(self.MuonIDSF_GH_th2, abs(eta), pt)
	#print "MuonID SF_BCDEF ",SF_BCDEF, errlow_BCDEF, errhigh_BCDEF, " SF_GH ",SF_GH, errlow_GH, errhigh_GH
	SF = SF_BCDEF*(1-self.MuonTrackingSF_GH_lumiratio) + SF_GH*self.MuonTrackingSF_GH_lumiratio
	SF_low = SF - combinedError(errlow_GH-SF_GH, errlow_BCDEF-SF_BCDEF, self.MuonTrackingSF_GH_lumiratio)
	SF_high = SF + combinedError(errhigh_GH-SF_GH, errhigh_BCDEF-SF_BCDEF, self.MuonTrackingSF_GH_lumiratio)
	return SF, SF_high, SF_low

    def getMuonIsoSF(self, eta, pt):
	SF_BCDEF, errlow_BCDEF, errhigh_BCDEF = self.getSF(self.MuonIsoSF_th2, abs(eta), pt)
	SF_GH, errlow_GH, errhigh_GH = self.getSF(self.MuonIsoSF_GH_th2, abs(eta), pt)
	#print "MuonIso SF_BCDEF ",SF_BCDEF, errlow_BCDEF, errhigh_BCDEF, " SF_GH ",SF_GH, errlow_GH, errhigh_GH
	SF = SF_BCDEF*(1-self.MuonTrackingSF_GH_lumiratio) + SF_GH*self.MuonTrackingSF_GH_lumiratio
	SF_low = SF - combinedError(errlow_GH-SF_GH, errlow_BCDEF-SF_BCDEF, self.MuonTrackingSF_GH_lumiratio)
	SF_high = SF + combinedError(errhigh_GH-SF_GH, errhigh_BCDEF-SF_BCDEF, self.MuonTrackingSF_GH_lumiratio)
	return SF, SF_high, SF_low


    #def getMuonTrgSF(self, eta, pt):
	#return self.getSF(self.MuonTrgSF_th2, abs(eta), pt)

    def getleptonTrgSF(self, pt, eta, leg):
	for bins in self.TriggerSFs[leg]:
	    if eta > bins[0] and eta <= bins[1] and pt > bins[2] and pt <= bins[3]:
		errhigh = 0.0; errlow = 0.0
		if len(bins) == 6:
		    errhigh = bins[5]; errlow = bins[5]
		else:
		    errhigh = bins[5]; errlow = bins[6]
		#print "leg ",leg," pt ",pt," eta ",eta, " weight ",bins[4]
		return bins[4],bins[4]+errhigh, bins[4]-errlow
	#print "leg ",leg," pt ",pt," eta ",eta ," noweight found, use 1.0"
	return 1.0, 1.0, 1.0


    def getleptonpairTrgSF(self, leptonpair):

	legs = []
	Dzeff = 1.0
	if abs(leptonpair[0].pdgId) == 11 and abs(leptonpair[1].pdgId) == 11:
	    legs.append("DoubleEleLegHigPt")
	    legs.append("DoubleEleLegLowPt")
	    Dzeff = self.DZEffs["DoubleEle"]
	elif  abs(leptonpair[0].pdgId) == 13 and abs(leptonpair[1].pdgId) == 13:
	    legs.append("DoubleMuLegHigPt")
	    legs.append("DoubleMuLegLowPt")
	    Dzeff = self.DZEffs["DoubleMu"]
	elif  abs(leptonpair[0].pdgId) == 11 and abs(leptonpair[1].pdgId) == 13:
	    legs.append("EleMuLegHigPt")
	    legs.append("EleMuLegLowPt")
	    Dzeff = self.DZEffs["EleMu"]
	elif  abs(leptonpair[0].pdgId) == 13 and abs(leptonpair[1].pdgId) == 11:
	    legs.append("MuEleLegHigPt")
	    legs.append("MuEleLegLowPt")
	    Dzeff = self.DZEffs["MuEle"]
	else:
	    raise ValueError('Getting lepton SFs, lepton pair is true! ')

	SF1 = self.getleptonTrgSF(leptonpair[0].pt, leptonpair[0].eta, legs[0])
       	SF2 = self.getleptonTrgSF(leptonpair[1].pt, leptonpair[1].eta, legs[1])
        return SF1[0]*SF2[0]*Dzeff,  SF1[1]*SF2[1], SF1[2]*SF2[2]

    def getleptonIsoSF(self, lep):
	if abs(lep.pdgId) == 11:
	    return 1.0,1.0,1.0
	elif abs(lep.pdgId) == 13:
	    return self.getMuonIsoSF(abs(lep.eta), lep.pt)
	else:
	    raise ValueError('Getting lepton Iso SFs, not electron nor muon ', lep.pdgId)

    def getleptonpairIsoSF(self, leptonpair):
	SF1 = self.getleptonIsoSF(leptonpair[0])
       	SF2 = self.getleptonIsoSF(leptonpair[1])
        return SF1[0]*SF2[0],  SF1[1]*SF2[1], SF1[2]*SF2[2]

    def getleptonIDSF(self, lep):
	if abs(lep.pdgId) == 11:
	    return self.getEGSF(lep.eta, lep.pt)
	elif abs(lep.pdgId) == 13:
	    return self.getMuonIDSF(abs(lep.eta), lep.pt)
	else:
	    raise ValueError('Getting lepton ID SFs, not electron nor muon ', lep.pdgId)

    def getleptonpairIDSF(self, leptonpair):
	SF1 = self.getleptonIDSF(leptonpair[0])
       	SF2 = self.getleptonIDSF(leptonpair[1])
        return SF1[0]*SF2[0],  SF1[1]*SF2[1], SF1[2]*SF2[2]

    def getMuonTrackingSF(self, eta):
	 for thisbin in self.MuonTrackingSF_allbins:
	     if eta > thisbin["etalow"] and eta <= thisbin["etahigh"]:
		 return thisbin["SF"],thisbin["SFerrhigh"], thisbin["SFerrlow"]
	 return  1.0, 1.0, 1.0


    def getleptonTrackingSF(self, lep):
	if abs(lep.pdgId) == 13:
	    return self.getMuonTrackingSF(abs(lep.eta))
	else:
	    return 1.0, 1.0, 1.0
    def getleptonpairTrackingSF(self, leptonpair):
	SF1 = self.getleptonTrackingSF(leptonpair[0])
       	SF2 = self.getleptonTrackingSF(leptonpair[1])
        return SF1[0]*SF2[0],  SF1[1]*SF2[1], SF1[2]*SF2[2]


	  
