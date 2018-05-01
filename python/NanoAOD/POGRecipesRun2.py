import ROOT
import string
from math import sqrt, pi, degrees
import os
import json
import numpy as np

def printObject(obj):
    print " Id ",obj.pdgId, " pt ",obj.pt, " eta ", obj.eta," phi ",obj.phi


def jetLooseID(jet):
    """https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_data"""
    ##only for abs(jet.eta)<2.7 , does not care about jet.eta > 2.7     
    #CHM = jet.nElectrons + jet.nMuons
    CHM = 1 ## this cut not available in NANoAOD?
    if abs(jet.eta)<2.7:
       ##looseJetID = (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4) && abs(eta)<=2.7
       return  (jet.neHEF < 0.99 and jet.neEmEF<0.99 and jet.nConstituents>1 and ((abs(jet.eta)<= 2.4 and jet.chHEF>0 and CHM>0 and jet.chEmEF<0.99) or (abs(jet.eta)>2.4)))
    else:
      return False
def electronImpactParameterCut( electron):
    """ check electron impact parameter, https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2 """
    dxy_endcap = 0.1 #cm
    dxy_barrel = 0.05
    dz_endcap  = 0.2
    dz_barrel = 0.1
    superCluster_eta = electron.eta + electron.deltaEtaSC
    if abs(superCluster_eta) <= 1.479: ##barrel
	return (abs(electron.dz) < dz_barrel and abs(electron.dxy) < dxy_barrel)
    elif abs(superCluster_eta) > 1.479 and abs(electron.eta) < 2.5:
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
    #return electron.pfRelIso03_all < 0.04
    return True

def electronHLTSafeID( electron, jetRhoCalo):
    """ check electron isolation, https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2 """
    """https://twiki.cern.ch/twiki/bin/view/CMS/ChangesEGMHLTAlgo2014"""
    """https://github.com/cms-sw/cmssw/blob/CMSSW_9_4_X/RecoEgamma/ElectronIdentification/python/Identification/cutBasedElectronHLTPreselecition_Summer16_V1_cff.py"""
    ### extrac HLT safe ID cut is applied for electron selection for 2016 data only
    ##
    #GsfEleDEtaInSeedCut 
    #ptCutOff  = 20.0
    barrelCutOff  = 1.479
    #isRelativeIso_Ecal = True
    #isRelativeIso_Hcal = True
    superCluster_eta = electron.eta + electron.deltaEtaSC
    sieie = electron.sieie
    hoe = electron.hoe
    eInvMinusPInv = abs(electron.eInvMinusPInv)
    ##Iso value here is not exact same as used in officialy analysis, Eff to be validated 
    TrkPtIsoRel =  electron.dr03TkSumPt/electron.pt
    ##https://indico.cern.ch/event/491507/contributions/2192817/attachments/1285452/1911768/EGM_HLTsafeCuts_31May16.pdf
    ##quote from HWW https://github.com/latinos/LatinoAnalysis/blob/master/Gardener/python/variables/l2Sel.py#L188

    eA_ecal = 0.165
    eA_hcal = 0.06
    if abs(superCluster_eta) > barrelCutOff:
        eA_ecal = 0.132
	eA_hcal = 0.131
        
    EcalPFClusterIsoRel = (electron.dr03EcalRecHitSumEt - eA_ecal * jetRhoCalo )/electron.pt
    HcalPFClusterIsoRel = (electron.dr03HcalDepth1TowerSumEt - eA_hcal * jetRhoCalo )/electron.pt
    normalizedGsfChi2cut = electron.isPFcand
    missingHits = electron.lostHits
    #print "electron pt ",electron.pt," supercluter_eta ",superCluster_eta, " sieie ",sieie, " hoe ",hoe," eInvMinusPInv ",eInvMinusPInv, " EcalPFClusterIsoRel ",EcalPFClusterIsoRel, " HcalPFClusterIsoRel ",HcalPFClusterIsoRel, " TrkPtIsoRel ",TrkPtIsoRel," jetRhoCalo ",jetRhoCalo

    if abs(superCluster_eta) <= barrelCutOff: #barrel
        return (abs(sieie) < 0.11 and hoe < 0.06 and eInvMinusPInv < 0.013 and EcalPFClusterIsoRel < 0.16 and HcalPFClusterIsoRel < 0.12 and TrkPtIsoRel < 0.08)
    else:
        return (abs(sieie) < 0.31 and hoe < 0.065 and eInvMinusPInv < 0.013 and EcalPFClusterIsoRel < 0.12 and HcalPFClusterIsoRel < 0.12 and TrkPtIsoRel < 0.08 and normalizedGsfChi2cut) 

def muonID( muon):
    """https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#Muon_Isolation """
    #goodGlob = (muon.isGlobalMuon() and muon.globalTrack()->normalizedChi2() < 3 and muon.combinedQuality().chi2LocalPosition < 12 and muon.combinedQuality().trkKink < 20)
    ###check Muon ID
    #return muon.mediumId
    return muon.tightId

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

def leptonHLTSafeID(lepton, jetRhoCalo):
    return ((abs(lepton.pdgId) == 11 and electronHLTSafeID(lepton, jetRhoCalo)) or (abs(lepton.pdgId) == 13))
def leptonpairHLTSafeID(leptonpair, jetRhoCalo):
    return (leptonHLTSafeID(leptonpair[0], jetRhoCalo) and leptonHLTSafeID(leptonpair[1], jetRhoCalo))

def jetMediumBtagging(jet):
    """https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco"""
    ### use medium btagging 
    ### SF should be from  leptonSF/cMVAv2_Moriond17_B_H.csv 
    return (jet.btagCMVA > 0.4432)

def combinedError(err1, err2, weight1):
    return sqrt(err1*err1*weight1+err2*err2*(1.0-weight1));

def phiconversion_radintodegree(phi):
    if phi < 0:
       phi = phi + pi*2
    if phi > 2*pi:
       phi = phi - pi*2 
    return degrees(phi)

def propagateToME2X_phi(lepton):
    ##From Louvain
    charge = lepton.charge
    pt = lepton.pt
    lep_p4 = ROOT.TLorentzVector()
    lep_p4.SetPtEtaPhiM(lepton.pt, lepton.eta, lepton.phi, lepton.mass)
    theta = degrees(lep_p4.Theta()) ##should be between 0, 180
    if theta > 90:
        theta = 180 - theta
    return  lep_p4.Phi() + pi/180.0 * charge * (1. / pt) * (10.48 - 5.1412 * theta + 0.02308 * theta * theta)

def checkMuonSector(lepton):
    eta = lepton.eta
    side = 1
    if eta < 0.0:	side = -1
    if(abs(eta)<1.25):		return 0
    elif(abs(eta)<1.8):
    #Code is: 1=> first 50 degrees of the sector, 11=> last 10 deg of previous sector
	phi_ME2x = propagateToME2X_phi(lepton)
	phiindegree = phiconversion_radintodegree( phi_ME2x )
        print "phi ",lepton.phi," phi in ME2x ", phi_ME2x, " phi in degree ",phiindegree
	if(phiindegree>=5   and phiindegree<15):     	return 11*side
	if(phiindegree>=15  and phiindegree<65):    	return 1*side
	if(phiindegree>=65  and phiindegree<75):    	return 12*side
	if(phiindegree>=75  and phiindegree<125):   	return 2*side
	if(phiindegree>=125 and phiindegree<135):  	return 13*side
	if(phiindegree>=135 and phiindegree<185):  	return 3*side
	if(phiindegree>=185 and phiindegree<195):  	return 14*side
	if(phiindegree>=195 and phiindegree<245):  	return 4*side
	if(phiindegree>=245 and phiindegree<255):  	return 15*side
	if(phiindegree>=255 and phiindegree<305):  	return 5*side
	if(phiindegree>=305 and phiindegree<315):  	return 16*side
	if(phiindegree>=315 or phiindegree < 5): 	return 6*side
    else:
    #Code is: 1=> first 40 degrees of the sector, 11=> last 20 deg of previous sector
	phi_ME2x = propagateToME2X_phi(lepton)
	phiindegree = phiconversion_radintodegree( phi_ME2x )
        print "phi ",lepton.phi," phi in ME2x ", phi_ME2x, " phi in degree ",phiindegree
	if(phiindegree>=355 or  phiindegree<15):   	return 11*side
	if(phiindegree>=15  and phiindegree<55):     	return 1*side
	if(phiindegree>=55  and phiindegree<75):     	return 12*side
	if(phiindegree>=75  and phiindegree<115):    	return 2*side
	if(phiindegree>=115 and phiindegree<135):   	return 13*side
	if(phiindegree>=135 and phiindegree<175):   	return 3*side
	if(phiindegree>=175 and phiindegree<195):   	return 14*side
	if(phiindegree>=195 and phiindegree<235):   	return 4*side
	if(phiindegree>=235 and phiindegree<255):   	return 15*side
	if(phiindegree>=255 and phiindegree<295):   	return 5*side
	if(phiindegree>=295 and phiindegree<315):   	return 16*side;
	if(phiindegree>=315 and phiindegree < 355): 	return 6*side;
    return 0


def checkMuonPairSectors(muon1, muon2):
    ### 0, good event, no SF
    ### 1, both in overlap region or non-overlap region
    ### 2, one in overlap region and one in non-overlap region
    if abs(muon1.pdgId) != 13 or abs(muon2.pdgId) != 13:
        print "error in checkMuonPairSectors!! not both inputs are muon "
        return 0
    if muon1.eta*muon2.eta < 0 or abs(muon1.eta)<1.25 or abs(muon2.eta)<1.25:
        return 0
    #print "muon1 ", printObject(muon1)
    sector1 = checkMuonSector(muon1)
    #print "muon2 ", printObject(muon2)
    sector2 = checkMuonSector(muon2)
    #print "sector1 ",sector1," sector2 ",sector2
    if sector1*sector2 < 0: ## not same endcap
        return 0
    if sector1 == sector2: 
        return 1
    if abs(sector1 - sector2) == 10:
        return 2
    return 0

def isMuonPairSameCSCRegion(muon1, muon2):
    return checkMuonPairSectors(muon1, muon2) == 1


def loadJsonFile(filename):
    return json.loads(open(filename).read())

class LeptonSFManager():
    ### to take the lumi into consideration?
    ### Electron Triggering, IP, Isolation, tracking?
    ### Muon tracking ?
    def __init__(self):
	##brilcalc lumi -u /pb  --normtag normtag_PHYSICS.json -i json.txt
	## used delivered lumi for normalization

	self.useJsonSFs = False
	##self.Lumi_BCDEF = 5.750+2.573+4.242+4.025+3.105 ## recorded,23Sep2016ReReco
	self.Lumi_BCDEF = 5.991+2.685+4.411+4.222+3.303 ##delievered
	##self.Lumi_GH = 7.576+8.651 ## recorded, 23Sep2016ReReco
	self.Lumi_GH = 7.865+8.985
	self.totalLumi = self.Lumi_BCDEF  + self.Lumi_GH#fb-1
	HhhPath = "/home/taohuang/DiHiggsAnalysis/CMSSW_9_4_0_pre1/src/HhhAnalysis/python/NanoAOD/"
	#self.EGSF_filename  = "leptonSF/EGM2D_eleGSF.root" ## for electron ID
	self.EGIDSF_filename  =             os.path.join(HhhPath , "leptonSF/egammaEffi.txt_EGM2D_CutBasedId.root" )## for electron ID
	self.EGRecoSF_filename  =             os.path.join(HhhPath , "leptonSF/egammaEffi.txt_EGM2D_reco.root" )## for electron ID
	self.MuonIDSF_filename =          os.path.join(HhhPath , "leptonSF/EfficienciesAndSF_BCDEF_ID.root")
	self.MuonIsoSF_filename =         os.path.join(HhhPath , "leptonSF/EfficienciesAndSF_BCDEF_ISO.root")
	self.MuonTrgSF_filename =         os.path.join(HhhPath , "leptonSF/EfficienciesAndSF_BCDEF_trigger.root")
	self.MuonTrackingSF_filename =    os.path.join(HhhPath , "leptonSF/EfficienciesAndSF_BCDEF_Tracking.root")
	self.MuonIDSF_GH_filename =       os.path.join(HhhPath , "leptonSF/EfficienciesAndSF_GH_ID.root")
	self.MuonIsoSF_GH_filename =      os.path.join(HhhPath , "leptonSF/EfficienciesAndSF_GH_ISO.root")
	self.MuonTrgSF_GH_filename =      os.path.join(HhhPath , "leptonSF/EfficienciesAndSF_GH_trigger.root")
	self.MuonTrackingSF_GH_filename = os.path.join(HhhPath , "leptonSF/EfficienciesAndSF_GH_Tracking.root")

	
	### x-axis: eta,  y-axis: pt
	self.EGIDSF_histname = "EGamma_SF2D"
	self.EGRecoSF_histname = "EGamma_SF2D"
	#self.MuonIDSF_histname = "MC_NUM_MediumID2016_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio"
	#self.MuonIsoSF_histname = "TightISO_MediumID_pt_eta/abseta_pt_ratio"
	self.MuonIDSF_histname = "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio"
	self.MuonIsoSF_histname = "TightISO_TightID_pt_eta/abseta_pt_ratio"
	self.MuonTrgSF_histname = "IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio"
	self.MuonTrackingSF_tgraphname = "ratio_eff_eta3_dr030e030_corr"

	self.EGIDSF_tfile = ROOT.TFile(self.EGIDSF_filename,"READ")
	self.EGRecoSF_tfile = ROOT.TFile(self.EGRecoSF_filename,"READ")
	self.MuonIDSF_tfile = ROOT.TFile(self.MuonIDSF_filename,"READ")
	self.MuonIsoSF_tfile = ROOT.TFile(self.MuonIsoSF_filename,"READ")
	#self.MuonTrgSF_tfile = ROOT.TFile(self.MuonTrgSF_filename,"READ")
	self.MuonTrackingSF_tfile = ROOT.TFile(self.MuonTrackingSF_filename,"READ")
	self.MuonIDSF_GH_tfile = ROOT.TFile(self.MuonIDSF_GH_filename,"READ")
	self.MuonIsoSF_GH_tfile = ROOT.TFile(self.MuonIsoSF_GH_filename,"READ")
	#self.MuonTrgSF_GH_tfile = ROOT.TFile(self.MuonTrgSF_GH_filename,"READ")
	self.MuonTrackingSF_GH_tfile = ROOT.TFile(self.MuonTrackingSF_GH_filename,"READ")

	self.EGIDSF_th2 = self.EGIDSF_tfile.Get(self.EGIDSF_histname)
	self.EGRecoSF_th2 = self.EGRecoSF_tfile.Get(self.EGRecoSF_histname)
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
	##apply for muon pairs in same endcap(abs(eta)>1.25) due to bug  in EMTF, 2016 Run
	## case1 both two muons in CSC overlap region or non overlap region => if EMTFBug, eff = 0
	## case2 one muon in overlap region and another in non overlap region => if EMTFBug, eff = 0.5
	## should use phi at CSC station2 !!!
	##https://twiki.cern.ch/twiki/bin/view/CMS/EndcapHighPtMuonEfficiencyProblem
	self.EMTFBug_run2016_sameOverlap_or_SameNonOverlap = 0.564474
	self.EMTFBug_run2016_oneOverlap_oneNonOverlap = 0.782237

	###SFs in Json files
	self.Ele_HLTSafeID_file = 'leptonSF/Electron_MediumPlusHLTSafeID_moriond17.json'
	self.Electron_MediumPlusHLTSafeID_moriond17 = loadJsonFile(self.Ele_HLTSafeID_file)

        self.Muon_id_jsonfile = os.path.join(HhhPath, "leptonSF/Muon_TightID_genTracks_id_BCDEFGH_weighted.json")
        self.Muon_iso_jsonfile = os.path.join(HhhPath, "leptonSF/Muon_TightISO_TightID_iso_BCDEFGH_weighted.json")
        self.Muon_reco_jsonfile = os.path.join(HhhPath, "leptonSF/Muon_tracking_BCDEFGH.json")

	self.Muon_id_dict = loadJsonFile(self.Muon_id_jsonfile)
	self.Muon_iso_dict = loadJsonFile(self.Muon_iso_jsonfile)
	self.Muon_reco_dict = loadJsonFile(self.Muon_reco_jsonfile)

        self.Ele_id_jsonfile = os.path.join(HhhPath, 'leptonSF/Electron_EGamma_SF2D_medium_moriond17.json')
        self.Ele_reco_jsonfile = os.path.join(HhhPath, 'leptonSF/Electron_EGamma_SF2D_reco_moriond17.json')
	 
	self.Ele_id_dict = loadJsonFile(self.Ele_id_jsonfile)
	self.Ele_reco_dict = loadJsonFile(self.Ele_reco_jsonfile)

    def useJsonFiles(self, x):
	self.useJsonSFs = x


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
	
    def getSF_json(self, SFs_dict, eta, pt):
	for etabin in SFs_dict:
	    etas = re.findall(r"[-+]?\d*\.\d+|\d+", etabin)
	    #print "etabin ",etabin," range ",etas
	    if not ((eta >= etas[0] and eta < etas[1]) or (eta >= etas[1] and eta < etas[0])):
		continue
	    for ptbin in SFs_dict[etabin]:
		pts = re.findall(r"[-+]?\d*\.\d+|\d+", ptbin)
		#print "ptbin ", ptbin, " range ",pts
		if (pt >= pts[0] and pt<pts[1]) or (pt <pts[0] and pt >= pts[1]):
		    SF = SFs_dict[etabin][ptbin]['value']
		    err = SFs_dict[etabin][ptbin]['error']
		    return SF, SF+error, SF-error
	return 1.0,1.0,1.0


    def getSF_json_v2(self, SFs_dict, eta, pt):
	#print "eta ",eta, " pt ",pt
	for etabin in SFs_dict:
	    etas = etabin['bin']
	    #print "etabin  ",etabin
	    if not ((eta >= etas[0] and eta < etas[1]) or (eta >= etas[1] and eta < etas[0])):
		continue
	    for ptbin in etabin['values']:
	        pts = ptbin['bin']
	        #print "ptbin ",ptbin
	        if (pt >= pts[0] and pt<pts[1]) or (pt <pts[0] and pt >= pts[1]):
		    SF = ptbin['value']
		    SF_errup = ptbin['error_high']
		    SF_errlow = ptbin['error_low']
		    #print " getSF_Ele_HLTSafeID ", SF, " ele pt ",pt, " eta ",eta
		    return SF, SF + SF_errup, SF - SF_errlow
	return 1.0,1.0,1.0
  
    def getleptonHLTSafeIDSF(self, lep):
        if abs(lep.pdgId) == 13:
            return  1.0,1.0,1.0
        elif  abs(lep.pdgId) == 11:	
	    superCluster_eta = lep.eta + lep.deltaEtaSC
	    return self.getSF_json_v2(self.Electron_MediumPlusHLTSafeID_moriond17['data'], superCluster_eta, lep.pt)

    def getleptonpairHTLSafeIDSF(self, leptonpair):
	SF1 = self.getleptonHLTSafeIDSF(leptonpair[0])
       	SF2 = self.getleptonHLTSafeIDSF(leptonpair[1])
        return SF1[0]*SF2[0],  SF1[1]*SF2[1], SF1[2]*SF2[2]


    def getEGIDSF(self, eta, pt):##final one ?
	SF = [1.0, 1.0, 1.0]
        if self.useJsonSFs:
	    SF = self.getSF_json_v2(self.Ele_id_dict['data'], eta, pt)
	else:
	    SF = self.getSF(self.EGIDSF_th2, eta, pt)
        #print "Electron SFs ",SF
	return SF

    def getEGRecoSF(self, eta, pt):##final one ?
	SF = [1.0, 1.0, 1.0]
        if self.useJsonSFs:
	    SF = self.getSF_json_v2(self.Ele_reco_dict['data'], eta, pt)
	else:
	    SF = self.getSF(self.EGRecoSF_th2, eta, pt)
        #print "Electron SFs ",SF
	return SF

    def getMuonIDSF(self, eta, pt):
	SF = [1.0, 1.0, 1.0]
        if self.useJsonSFs:
	    SF = self.getSF_json_v2(self.Muon_id_dict['data'], eta, pt)
	else:
	    SF_BCDEF, errlow_BCDEF, errhigh_BCDEF = self.getSF(self.MuonIDSF_th2, abs(eta), pt)
	    SF_GH, errlow_GH, errhigh_GH = self.getSF(self.MuonIDSF_GH_th2, abs(eta), pt)
	    #print "MuonID SF_BCDEF ",SF_BCDEF, errlow_BCDEF, errhigh_BCDEF, " SF_GH ",SF_GH, errlow_GH, errhigh_GH
	    SF = []
	    SF.append( SF_BCDEF*(1-self.MuonTrackingSF_GH_lumiratio) + SF_GH*self.MuonTrackingSF_GH_lumiratio )
	    SF.append( SF - combinedError(errlow_GH-SF_GH, errlow_BCDEF-SF_BCDEF, self.MuonTrackingSF_GH_lumiratio))
	    SF.append( SF + combinedError(errhigh_GH-SF_GH, errhigh_BCDEF-SF_BCDEF, self.MuonTrackingSF_GH_lumiratio))
	return SF

    def getMuonIsoSF(self, eta, pt):
	SF = [1.0, 1.0, 1.0]
        if self.useJsonSFs:
	    SF = self.getSF_json_v2(self.Muon_iso_dict['data'], eta, pt)
	else:
	    SF_BCDEF, errlow_BCDEF, errhigh_BCDEF = self.getSF(self.MuonIsoSF_th2, abs(eta), pt)
	    SF_GH, errlow_GH, errhigh_GH = self.getSF(self.MuonIsoSF_GH_th2, abs(eta), pt)
	    #print "MuonIso SF_BCDEF ",SF_BCDEF, errlow_BCDEF, errhigh_BCDEF, " SF_GH ",SF_GH, errlow_GH, errhigh_GH
	    SF.append(SF_BCDEF*(1-self.MuonTrackingSF_GH_lumiratio) + SF_GH*self.MuonTrackingSF_GH_lumiratio)
	    SF.append(SF - combinedError(errlow_GH-SF_GH, errlow_BCDEF-SF_BCDEF, self.MuonTrackingSF_GH_lumiratio))
	    SF.append(SF + combinedError(errhigh_GH-SF_GH, errhigh_BCDEF-SF_BCDEF, self.MuonTrackingSF_GH_lumiratio))
	return SF


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
	    cscsector_case  = checkMuonPairSectors(leptonpair[0], leptonpair[1])
    	    EMTFeff = 1.0
            if cscsector_case == 1:
	       EMTFeff = self.EMTFBug_run2016_sameOverlap_or_SameNonOverlap 
	    elif cscsector_case == 2:
	       EMTFeff =  self.EMTFBug_run2016_oneOverlap_oneNonOverlap
	    Dzeff = self.DZEffs["DoubleMu"]*EMTFeff
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
        #print "IsoSF lep1 ",printObject(leptonpair[0]), " SF1 ", SF1, " lep2 ",printObject(leptonpair[1])," SF2 ",SF2 
        return SF1[0]*SF2[0],  SF1[1]*SF2[1], SF1[2]*SF2[2]

    def getleptonIDSF(self, lep):
	if abs(lep.pdgId) == 11:
	    return self.getEGIDSF(lep.eta + lep.deltaEtaSC, lep.pt)#use super cluster eta
	elif abs(lep.pdgId) == 13:
	    return self.getMuonIDSF(abs(lep.eta), lep.pt)
	else:
	    raise ValueError('Getting lepton ID SFs, not electron nor muon ', lep.pdgId)

    def getleptonpairIDSF(self, leptonpair):
	SF1 = self.getleptonIDSF(leptonpair[0])
       	SF2 = self.getleptonIDSF(leptonpair[1])
        #print "IDSF lep1 ",printObject(leptonpair[0]), " SF1 ", SF1, " lep2 ",printObject(leptonpair[1])," SF2 ",SF2 
        return SF1[0]*SF2[0],  SF1[1]*SF2[1], SF1[2]*SF2[2]

    def getMuonTrackingSF(self, eta, pt):
        if self.useJsonSFs:
	    SF = self.getSF_json_v2(self.Muon_reco_dict['data'], eta, pt)
	    return SF
	else:
	    eta = abs(eta)
	    for thisbin in self.MuonTrackingSF_allbins:
		if eta > thisbin["etalow"] and eta <= thisbin["etahigh"]:
		    return thisbin["SF"],thisbin["SFerrhigh"], thisbin["SFerrlow"]
	return  1.0, 1.0, 1.0


    def getleptonTrackingSF(self, lep):
	if abs(lep.pdgId) == 13:
	    return self.getMuonTrackingSF(lep.eta, lep.pt)
	elif abs(lep.pdgId) == 11:
	    return self.getEGRecoSF(lep.eta + lep.deltaEtaSC, lep.pt)
	else:
	    raise ValueError('Getting lepton ID SFs, not electron nor muon ', lep.pdgId)

    def getleptonpairTrackingSF(self, leptonpair):
	SF1 = self.getleptonTrackingSF(leptonpair[0])
       	SF2 = self.getleptonTrackingSF(leptonpair[1])
        #print "RecoSF lep1 ",printObject(leptonpair[0]), " SF1 ", SF1, " lep2 ",printObject(leptonpair[1])," SF2 ",SF2 
        return SF1[0]*SF2[0],  SF1[1]*SF2[1], SF1[2]*SF2[2]


	  
