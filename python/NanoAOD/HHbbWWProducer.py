import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
print "ROOT version ",ROOT.gROOT.GetVersion()
from math import sqrt,cos
import random
#random.randint(0,100)

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection,Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import * #deltaR, matching etc..
import sys
sys.path.append('/home/taohuang/DiHiggsAnalysis/CMSSW_9_4_0_pre1/src/HhhAnalysis/python/NanoAOD')

from RunConfiguration import *
if Runyear == 2016:
    import POGRecipesRun2016 as POGRecipesRun2
elif Runyear == 2017:
    import POGRecipesRun2017 as POGRecipesRun2
else:
    print "wrong year run: %d !! "%Runyear
    sys.exit()

print "HHbbWWProducer, import finished here"

class HHbbWWProducer(Module):
    ## data or MC, which L1 trigger, HLT?
    ###kwargs: triggertype, verbose, run_lumi
    def __init__(self,isMC, **kwargs):
	print "init HHbbWWProduer"
	self.writeHistFile=True
        self.isMC = isMC ## two mode: data or MC
	print "kwargs ",kwargs
	self.triggertype  = ''##"DoubleMuon, DoubleEG, MuonEG"
	self.DYestimation = False
	self.CheckBtaggingEff = False
	self.deltaR_trigger_reco = 0.1; self.deltaPtRel_trigger_reco = 0.5
	self.verbose = 3
	self.muonEta = 2.4; self.EGEta = 2.5
	self.leadingMuonPt = {"DoubleMuon": 20.0, "MuonEG":25.0}
	self.subleadingMuonPt = {"DoubleMuon": 10.0, "MuonEG":10.0}
	self.leadingEGPt = {"DoubleEG": 25.0, "MuonEG":25.0}
	self.subleadingEGPt = {"DoubleEG": 15.0, "MuonEG":15.0}
	self.jetPt = 20; self.jetEta = 2.4
	self.deltaR_j_l = 0.3 #jet,lepton seperation
	self.maxnjets = 20

	self.run_lumi = None
	self.__dict__.update((key, kwargs[key]) for key in kwargs.keys())
	self.CheckBtaggingEff = (self.CheckBtaggingEff and self.DYestimation and self.isMC)
        print "self.run_lumi  ",self.run_lumi," trigger type ",self.triggertype," verbose ",self.verbose, " DYestimation ",self.DYestimation

	self.addGenToTree = False
	self.addSystematic = True

	self.ll_sys_branchname = ["Isosf","IDsf","trgsf","trackingsf","HLTsafeIDsf"]
	self.jjbtag_sys_branchname = ["jjbtag_light","jjbtag_heavy"]
	self.lep1_branchname_sysvalue = {}
	self.lep2_branchname_sysvalue = {}
	self.jjbtag_branchname_sysvalue = {}


	### for debug
	self.nEv_DoubleEG = 0
	self.nEv_DoubleMuon = 0
	self.nEv_MuonEG = 0
	self.lepSFmanager = POGRecipesRun2.LeptonSFManager()
        self.lepSFmanager.useJsonFiles(True)

    def beginJob(self, histFile=None,histDirName=None):
	print "BeginJob "
	#Module.beginJob(self,histFile,histDirName)
	#self.addObject(self.h_cutflow)
	#self.addObject(self.h_cutflow_weight)
	#self.addObject(self.h_cutflowlist["DoubleMuon"])
	#self.addObject(self.h_cutflowlist["DoubleEG"])
	#self.addObject(self.h_cutflowlist["MuonEG"])

    def endJob(self):
        pass


    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
	print "BeginFiles "
	self.h_eventcounter = ROOT.TH1F("h_raweventcounter","h_raweventcounter", 20, 0, 20.0)
	self.h_cutflow = ROOT.TH1F("h_cutflow","h_cutflow", 20, 0, 20.0)
	self.h_cutflow_weight = ROOT.TH1F("h_cutflow_weight","h_cutflow_weight", 200, 0, 200.0)
	self.h_cutflowlist = {"DoubleMuon":ROOT.TH1F("h_cutflow_DoubleMuon","h_cutflow_DoubleMuon", 20, 0, 20.0), 
			      "DoubleEG":ROOT.TH1F("h_cutflow_DoubleEG","h_cutflow_DoubleEG", 20, 0, 20.0), 
			      "MuonEG":ROOT.TH1F("h_cutflow_MuonEG","h_cutflow_MuonEG", 20, 0, 20.0)}

	self.h_cutflowlist["DoubleMuon"].SetLineColor(ROOT.kRed)
	self.h_cutflowlist["DoubleEG"].SetLineColor(ROOT.kBlue)
	self.h_cutflowlist["MuonEG"].SetLineColor(ROOT.kBlack)
	
        self.out = wrappedOutputTree
        self.out.branch("jet1_idx",  "I");
        self.out.branch("jet1_pt",  "F");
        self.out.branch("jet1_E",  "F");
        self.out.branch("jet1_eta",  "F");
        self.out.branch("jet1_phi",  "F");
        self.out.branch("jet1_cMVAv2",  "F");
        self.out.branch("jet1_partonFlavour",  "I");
        self.out.branch("jet1_hadronFlavour",  "I");
        self.out.branch("jet2_idx",  "I");
        self.out.branch("jet2_pt",  "F");
        self.out.branch("jet2_E",  "F");
        self.out.branch("jet2_eta",  "F");
        self.out.branch("jet2_phi",  "F");
        self.out.branch("jet2_cMVAv2",  "F");
        self.out.branch("jet2_partonFlavour",  "I");
        self.out.branch("jet2_hadronFlavour",  "I");
	self.out.branch("isElEl",  "I") # 0 or 1
	self.out.branch("isElMu",  "I") # 0 or 1, mu_pt<el_pt
	self.out.branch("isMuEl",  "I") # 0 or 1
	self.out.branch("isMuMu",  "I") # 0 or 1
	self.out.branch("isSF",  "F")
	#self.out.branch("lepstype",  "F")
        self.out.branch("lep1_pt",  "F");
        self.out.branch("lep1_E",  "F");
        self.out.branch("lep1_eta",  "F");
        self.out.branch("lep1_phi",  "F");
        self.out.branch("lep1_pdgid",  "I");
        self.out.branch("lep1_id",  "F");
        self.out.branch("lep1_iso",  "F");
        self.out.branch("lep2_pt",  "F");
        self.out.branch("lep2_E",  "F");
        self.out.branch("lep2_eta",  "F");
        self.out.branch("lep2_phi",  "F");
        self.out.branch("lep2_pdgid",  "I");
        self.out.branch("lep2_id",  "F");
        self.out.branch("lep2_iso",  "F");
	self.out.branch("met_pt",  "F")
	self.out.branch("met_phi",  "F")
	self.out.branch("nJetsL",  "F")
	self.out.branch("jj_M",  "F")
	self.out.branch("el_hltsafeid", "F")
	self.out.branch("llreco_weight",  "F")
	self.out.branch("llid_weight",  "F")
	self.out.branch("lliso_weight",  "F")
	self.out.branch("lltrigger_weight",  "F")
	self.out.branch("ht_jets",  "F")
	self.out.branch("ht",  "F")
	self.out.branch("ll_M",  "F")
	self.out.branch("llmet_M",  "F")
	self.out.branch("llmetjj_MT2",  "F")
	self.out.branch("llmetjj_M",  "F")
	self.out.branch("lljj_M",  "F")
	self.out.branch("cosThetaStar",  "F")
	self.out.branch("ll_DR_l_l",  "F")
	self.out.branch("jj_DR_j_j",  "F")
	self.out.branch("llmetjj_DPhi_ll_jj",  "F")
	self.out.branch("llmetjj_DPhi_ll_met",  "F")
	self.out.branch("llmetjj_DPhi_llmet_jj",  "F")
	self.out.branch("ll_pt",  "F")
	self.out.branch("ll_eta",  "F")
	self.out.branch("jj_pt",  "F")
	self.out.branch("jj_eta",  "F")
	self.out.branch("llmetjj_minDR_l_j",  "F")
	self.out.branch("llmetjj_MTformula",  "F")
	self.out.branch("ll_DPhi_l_l",  "F")
	self.out.branch("ll_DEta_l_l",  "F")
	self.out.branch("event_pu_weight",  "F")
	self.out.branch("event_lep_weight",  "F")
	self.out.branch("event_btag_weight",  "F")
	self.out.branch("jjbtag_heavy",  "F")
	self.out.branch("jjbtag_light",  "F")
        self.out.branch("genweight",  "F");
        self.out.branch("sample_weight",  "F");
	self.out.branch("event_reco_weight",  "F")
	self.out.branch("total_weight",  "F")
	self.out.branch("pu",  "F")
	self.out.branch("event_run",  "I")
	self.out.branch("event_lumiblock",  "I")
	#self.out.branch("event_weight",  "F")
	#self.out.branch("DY_BDT_flat",  "F")
	#self.out.branch("dy_nobtag_to_btagM_weight",  "F")
	#self.out.branch("mt2",  "F")
	#self.out.branch("mt2_bb",  "F")
	#self.out.branch("mt2_ll",  "F")
	#self.out.branch("event_number",  "I")

	if self.isMC and self.addSystematic:
	    for name in self.ll_sys_branchname:
		self.out.branch("lep1%s_up"%name,"F")
		self.out.branch("lep1%s_down"%name,"F")
		self.out.branch("lep2%s_up"%name,"F")
		self.out.branch("lep2%s_down"%name,"F")
	    for name in self.jjbtag_sys_branchname:
		self.out.branch("%s_up"%name,"F")
		self.out.branch("%s_down"%name,"F")

	    self.out.branch("event_pu_weight_up",  "F")
	    self.out.branch("event_pu_weight_down",  "F")
	    self.out.branch("event_pdf_weight_up",  "F")
	    self.out.branch("event_pdf_weight_down",  "F")
	    #self.out.branch("LHEScaleWeight", "F", lenVar = "nLHEScaleWeight")
	    for i in range(0, 9):
		self.out.branch("LHEScaleWeight_%d"%i, "F")
	if self.DYestimation or self.addGenToTree:	
	    self.out.branch("genjet1_partonFlavour",  "I");
	    self.out.branch("genjet2_partonFlavour",  "I");
	    self.out.branch("genmet_pt",  "F");
	    self.out.branch("genmet_phi",  "F");
	##how to add gen information???
	if self.isMC and self.addGenToTree:
	    self.out.branch("genb1_pt",  "F");##quark level
	    self.out.branch("genjet1_pt",  "F");
	    self.out.branch("genjet1_E",  "F");
	    self.out.branch("genjet1_eta",  "F");
	    self.out.branch("genjet1_phi",  "F");
	    self.out.branch("genb2_pt",  "F");
	    self.out.branch("genjet2_pt",  "F");
	    self.out.branch("genjet2_E",  "F");
	    self.out.branch("genjet2_eta",  "F");
	    self.out.branch("genjet2_phi",  "F");
	    #self.out.branch("genW1_pt",  "F");
	    #self.out.branch("genW1_mass",  "F");
	    #self.out.branch("genW2_pt",  "F");
	    #self.out.branch("genW2_mass",  "F");
	    #self.out.branch("gennu1_pt",  "F");##nu1 from W1
	    #self.out.branch("gennu1_eta",  "F");
	    #self.out.branch("gennu1_phi",  "F");
	    #self.out.branch("gennu2_pt",  "F");
	    #self.out.branch("gennu2_eta",  "F");
	    #self.out.branch("gennu2_phi",  "F");
	    self.out.branch("genl1_pt",  "F");##genl1 from lep1, leading 
	    self.out.branch("genl1_E",  "F");
	    self.out.branch("genl1_eta",  "F");
	    self.out.branch("genl1_phi",  "F");
	    self.out.branch("genl1_id",  "I");
	    self.out.branch("genl2_pt",  "F");
	    self.out.branch("genl2_E",  "F");
	    self.out.branch("genl2_eta",  "F");
	    self.out.branch("genl2_phi",  "F");
	    self.out.branch("genl2_id",  "I");
	    self.out.branch("dR_genb1_genjet1",  "F");
	    self.out.branch("dR_genb1_genjet2",  "F");
	    self.out.branch("dR_genb2_genjet1",  "F");
	    self.out.branch("dR_genb2_genjet2",  "F");
	    self.out.branch("dR_genl1_lepFromW1",  "F");
	    self.out.branch("dR_genl1_lepFromW2",  "F");
	    self.out.branch("dR_genl2_lepFromW1",  "F");
	    self.out.branch("dR_genl2_lepFromW2",  "F");
	    self.out.branch("met_diNuetrino_pt",  "F");
	    self.out.branch("met_diNuetrino_phi",  "F");
	    if self.CheckBtaggingEff:
		self.out.branch("alljets_pt", "F", n=self.maxnjets)
		self.out.branch("alljets_eta", "F", n=self.maxnjets)
		self.out.branch("alljets_cMVAv2", "F", n=self.maxnjets)
		#self.out.branch("alljets_partonFlavour", "I", n=self.maxnjets)
		#self.out.branch("alljets_hadronFlavour", "I", n=self.maxnjets)
		self.out.branch("alljets_genpartonFlavour", "I", n=self.maxnjets)
	
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
	self.h_eventcounter.Write()
	self.h_cutflow.Write()
	self.h_cutflowlist["DoubleMuon"].Write()
	self.h_cutflowlist["DoubleEG"].Write()
	self.h_cutflowlist["MuonEG"].Write()
	self.h_cutflow_weight.Write()
	#inputTree.SetBranchStatus("Jet_btagSF_*", 0)
	#c1 = ROOT.TCanvas("h_cutflow","h_cutflow")
	#self.h_cutflowlist["MuonEG"].Draw("hist")
	#self.h_cutflowlist["DoubleEG"].Draw("histsame")
	#self.h_cutflowlist["DoubleMuon"].Draw("histsame")
	#c1.Write()
	#c2 = ROOT.TCanvas("h_cutflow_weight","h_cutflow_weight")
	#c2.cd()
	#self.h_cutflow_weight.Draw("hist")
	#c2.Write()

    #def goldenJason(self, run, lumi):
    #    run_str = '%d'%run
    #    if run_str in self.run_lumi.keys():
    #        alllumis = self.run_lumi[run_str]
    #        for lumirange in alllumis:
    #            if lumi >= lumirange[0] and lumi <= lumirange[1]:
    #                return True
    #    else:
    #        return False


    def findingLeptonPairs(self, leptons_mu, leptons_el):

	leptonpairs = []
	alltriggertypes = []
	if self.isMC:
	    alltriggertypes = ["DoubleMuon", "DoubleEG", "MuonEG"]
	else:
	    alltriggertypes = [self.triggertype]


	for triggertype in alltriggertypes:
	    if triggertype == "DoubleMuon" and len(leptons_mu) >= 2:
		nmuons = len(leptons_mu)
		for imu1 in range(0, nmuons):
		    for imu2 in range(imu1+1, nmuons):
			if leptons_mu[imu1].charge * leptons_mu[imu2].charge < 0 and leptons_mu[imu1].pt > self.leadingMuonPt["DoubleMuon"] and leptons_mu[imu2].pt > self.subleadingMuonPt["DoubleMuon"]:
			    leptonpairs.append(( leptons_mu[imu1],  leptons_mu[imu2] ))

	    elif triggertype == "DoubleEG" and len(leptons_el) >= 2:
		nelectrons = len(leptons_el)
		for iel1 in range(0, nelectrons):
		    for iel2 in range(iel1+1, nelectrons):
			if leptons_el[iel1].charge * leptons_el[iel2].charge < 0 and leptons_el[iel1].pt > self.leadingEGPt["DoubleEG"] and leptons_el[iel2].pt > self.subleadingEGPt["DoubleEG"]:
			    leptonpairs.append(( leptons_el[iel1],  leptons_el[iel2] ))

	    elif triggertype == "MuonEG" and len(leptons_mu) >= 1 and len(leptons_el) >= 1:
		nmuons = len(leptons_mu)
		nelectrons = len(leptons_el)
		for imu1 in range(0, nmuons):
		    for iel1 in range(0, nelectrons):
			if leptons_mu[imu1].charge * leptons_el[iel1].charge > 0:
			    continue
			if (leptons_mu[imu1].pt >= self.leadingMuonPt["MuonEG"] and leptons_el[iel1].pt >= self.subleadingEGPt["MuonEG"]) or (leptons_mu[imu1].pt >= self.subleadingMuonPt["MuonEG"] and leptons_el[iel1].pt >= self.leadingEGPt["MuonEG"]):
			    if leptons_mu[imu1].pt > leptons_el[iel1].pt:
				leptonpairs.append(( leptons_mu[imu1],  leptons_el[iel1] ))
			    else:
				leptonpairs.append(( leptons_el[iel1],  leptons_mu[imu1] ))

	if self.verbose > 2:
	    for leps in leptonpairs:
		print "Lepton pair leading ",leps[0].pdgId, " pt ",leps[0].pt," subleading ",leps[1].pdgId," pt ",leps[1].pt
        return leptonpairs


    def pt(self, jet, isMC):
        ## the MC has JER smearing applied which has output branch Jet_pt_smeared which should be compared 
        ## with data branch Jet_pt. This essentially aliases the two branches to one common jet pt variable.
        if isMC:
            return jet.pt_smeared
        else:
            return jet.pt		    
    
    def met(self, met, isMC):
        ## the MC has JER smearing applied which has output branch met_[pt/phi]_smeared which should be compared 
        ## the MC has JER smearing applied which has output branch met_[pt/phi]_smeared which should be compared 
        ## the MC has JER smearing applied which has output branch met_[pt/phi]_smeared which should be compared 
        ## with data branch MET_[pt/phi]. This essentially aliases the two branches to one common variable.
        if isMC:
            return (met.pt_smeared,met.phi_smeared)
        else:
            return (met.pt,met.phi)
	   
    def matchTriggerObjwithHLT(self, path, trigobjs):
	alllegs = path.split('_')[1:-1]
	leg_trigobj_map  = {};
        leg_trigobj_map['leg1'] = []
        leg_trigobj_map['leg2'] = []
	foundleg1 = False; foundleg2 = False
	for x in alllegs:
	    legobjmap = None
	    if "Mu" in x:
	    	pt = 8.0
	    	if "TkMu" not in x:
		    pt = int(x[2:])
		#matchedobjs = list(filter(lambda obj: obj.id == 13 and obj.pt >= pt and obj.l1pt > pt/2.0 and (obj.filterBits & 3) > 0, trigobjs))
		matchedobjs = list(filter(lambda obj: obj.id == 13 and obj.pt >= pt and obj.filterBits > 0 , trigobjs))

	        if not(foundleg1):
		    leg_trigobj_map['leg1'] = matchedobjs
		    foundleg1 = True
		elif not (foundleg2):
		    leg_trigobj_map['leg2'] = matchedobjs
		    foundleg2 = True
		else:
		    print "error in matching HLT to trigger objects ", path
	    elif "Ele" in x:
		pt = int(x[3:])
		matchedobjs = list(filter(lambda obj: obj.id == 11 and obj.pt >= pt and (obj.filterBits & 1) > 0, trigobjs))
	        if not(foundleg1):
		    leg_trigobj_map['leg1'] = matchedobjs
		    foundleg1 = True
		elif not (foundleg2):
		    leg_trigobj_map['leg2'] = matchedobjs
		    foundleg2 = True
		else:
		    print "error in matching HLT to trigger objects ", path

	if self.verbose > 3: 
	    print "firedpath ", path, " matched l1objs, leg1 ",len(leg_trigobj_map['leg1']), " leg2 ",len(leg_trigobj_map['leg2'])

	return leg_trigobj_map 
         
	#selectedObjs = []
	#for x in alllegs:
	#    if "Mu" in x:
	#    	pt = 8.0
	#    	if "TkMu" not in x:
	#	    pt = int(x[2:])
	#	for obj in trigobjs:
	#	    #if obj.id == 13 and (obj.filterBits & 1)>0 and obj.pt >= pt and (obj.l1pt >= 4.0 or obj.l1pt_2 >= 4.0):
	#	    #if obj.id == 13 and obj.pt >= pt:
	#	    #if obj.id == 13 and obj.pt >= pt and obj.l1pt >= 4.0 and (obj.filterBits & 3) > 0:
	#	    if obj.id == 13 and obj.pt >= pt and (obj.filterBits & 3) > 0:
	#	    	selectedObjs.append(obj)
	#	    elif  obj.id == 13 and  self.verbose >3 :
	#	        print "path ",path, " trig pt ",pt, " trig obj pt ",obj.pt, " l1pt ",obj.l1pt," filterbits ",obj.filterBits
	#    elif "Ele" in x:
	#	pt = int(x[3:])
	#	for obj in trigobjs:
	#	    #if obj.id == 11 and (obj.filterBits & 1)>0 and obj.pt >= pt and (obj.l1pt >= 10.0 or obj.l1pt_2 >= 10.0):##mini pt for L1 seed? 10GeV?
	#	    #if obj.id == 11 and obj.pt >= pt:##mini pt for L1 seed? 10GeV?
	#	    if obj.id == 11 and obj.pt >= pt and  (obj.filterBits & 1) > 0:##mini pt for L1 seed? 10GeV?
	#	    	selectedObjs.append(obj)
	#	    elif  obj.id == 11 and  self.verbose >3 :
	#	        print "path ",path, " trig pt ",pt, " trig obj pt ",obj.pt, " l1pt ",obj.l1pt," filterbits ",obj.filterBits

	#return set(selectedObjs)

    def findTriggerType(self, leptons):
	if abs(leptons[0].pdgId) == 13 and abs(leptons[1].pdgId) == 13:
	    return "DoubleMuon"
	elif abs(leptons[0].pdgId) == 11 and  abs(leptons[1].pdgId) == 11:
	    return "DoubleEG"
	elif (abs(leptons[0].pdgId) == 13 and  abs(leptons[1].pdgId) == 11) or (abs(leptons[0].pdgId) == 11 and  abs(leptons[1].pdgId) == 13):
	    return "MuonEG"
	else:
	    print "in findTriggerType, lepton id is 11 or 13, error!!!"
	    return ""

    def matchHLTPath(self, hlt, trigobjs, leptons):
	##triggerobj filterBits: 
	##1 = CaloIdL_TrackIdL_IsoVL, 2 = WPLoose, 4 = WPTight, 8 = OverlapFilter PFTau for Electron (PixelMatched e/gamma); 
	##1 = TrkIsoVVL, 2 = Iso, 4 = OverlapFilter PFTau for Muon

	##whether cut on l1pt for triggerobjects?
	## DoubleMu, L1 seed: DoubleMu_11_4 or 12_5
	## MuonEG: L1 seed: Mu20_EG10, L1_Mu5_IsoEG18..
	##DoubleEG: L1_SingleIsoEG22, L1_DoubleEG_15_10...


	""" check which HLT is fired """
	L1t_hlt = {"DoubleMuon": ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
				  "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ"],
	    	   "MuonEG": ["HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL","HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 
		   		"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL","HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ"],
		   "DoubleEG": ["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"],
			}
        def findfiredPaths(l1trigger):
	    return filter(lambda path : hasattr(hlt, path) and getattr(hlt,  path, False), L1t_hlt[l1trigger])
	def printLepton(lep, string):
	    print string," pdgid ",lep.pdgId," pt ",lep.pt, " eta ",lep.eta, " phi ",lep.phi
	def printTrigObjs(objs, string):
	    print string
	    for obj in objs:
	        print "trigger obj id ",obj.id, " pt ",obj.pt, " eta ",obj.eta, " phi ", obj.phi 
	def l1ptcutForPath(path):
	    l1ptcut = {}
	    if "Mu17" in path and "Mu8" in path:
	        l1ptcut["leg1"] = 11; l1ptcut["leg2"] = 4
	    elif "Mu23" in path and "Ele12" in path:
	        l1ptcut["leg1"] = 20; l1ptcut["leg2"] = 10
	    elif "Mu8" in path and "Ele23" in path:
	        l1ptcut["leg1"] = 5; l1ptcut["leg2"] = 18
	    elif "Ele23" in path and "Ele12" in path:
	        l1ptcut["leg1"] = 15; l1ptcut["leg2"] = 10
	    return l1ptcut
	def l1ptcut(path, leg1Objs, leg2Objs):
	    if len(leg1Objs)  == 0 or len(leg2Objs) == 0:
	        return False
	    l1ptcut = l1ptcutForPath(path)
	    leg1Objs_l1ptcut = list(filter(lambda obj : obj.l1pt >= l1ptcut["leg1"], leg1Objs))
	    leg2Objs_l1ptcut = list(filter(lambda obj : obj.l1pt >= l1ptcut["leg2"], leg2Objs))
	    passl1pt = (len(leg1Objs_l1ptcut)  >= 1 and len(leg2Objs_l1ptcut) >= 1 and not(leg1Objs_l1ptcut == leg2Objs_l1ptcut and len(leg1Objs_l1ptcut) == 1))
	    if path == "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ":
		leg1Objs_l1ptcut = list(filter(lambda obj : obj.l1pt >= 22, leg1Objs))
		leg2Objs_l1ptcut = list(filter(lambda obj : obj.l1pt >= 22, leg2Objs))
		passl1pt = passl1pt or ((len(leg1Objs_l1ptcut)>=1  or len(leg2Objs_l1ptcut) >= 1) and not(leg1Objs == leg2Objs and len(leg1Objs) == 1))
	    return passl1pt


        
	trigobjs.sort(key=lambda x: x.pt,reverse=True)##sort trigobjs by pt increasing order
	l1t = self.findTriggerType(leptons)	
        if not  self.isMC and l1t != self.triggertype:
	     print  "error!!! for data, the lepton's type is not the same as trigger type "
	     return False
        allfiredHLTs = findfiredPaths(l1t)
        if self.verbose > 3:
	    print "all fired paths ",allfiredHLTs
	    printLepton(leptons[0], "leading lep"); printLepton(leptons[1], "subleading lep")

        matchTolep1 = lambda tobj : tobj.id == abs(leptons[0].pdgId) and deltaR(tobj.eta, tobj.phi, leptons[0].eta, leptons[0].phi) < self.deltaR_trigger_reco and abs(tobj.pt-leptons[0].pt)/leptons[0].pt < self.deltaPtRel_trigger_reco 
        matchTolep2 = lambda tobj : tobj.id == abs(leptons[1].pdgId) and deltaR(tobj.eta, tobj.phi, leptons[1].eta, leptons[1].phi) < self.deltaR_trigger_reco and abs(tobj.pt-leptons[1].pt)/leptons[1].pt < self.deltaPtRel_trigger_reco 
        for  path in allfiredHLTs:
	    leg_trigobj_map = self.matchTriggerObjwithHLT(path, trigobjs)
            goodLeg1objs_lep1 = list(filter(matchTolep1, leg_trigobj_map['leg1']))
            goodLeg2objs_lep1 = list(filter(matchTolep1, leg_trigobj_map['leg2']))
            goodLeg1objs_lep2 = list(filter(matchTolep2, leg_trigobj_map['leg1']))
            goodLeg2objs_lep2 = list(filter(matchTolep2, leg_trigobj_map['leg2']))
	    if self.verbose > 3: 
		printTrigObjs( leg_trigobj_map['leg1'], "triggerobj matched to leg1")
		printTrigObjs( leg_trigobj_map['leg2'], "triggerobj matched to leg2")
		print "goodLeg1objs_lep1 ",len(goodLeg1objs_lep1), " goodLeg2objs_lep1 ",len(goodLeg2objs_lep1)," goodLeg1objs_lep2 ",len(goodLeg1objs_lep2)," goodLeg2objs_lep2 ",len(goodLeg2objs_lep2)
	         
	    #fired_option1 = (len(goodLeg1objs_lep1) > 0 and len(goodLeg2objs_lep2)>0 and not(goodLeg1objs_lep1 == goodLeg2objs_lep2 and len(goodLeg1objs_lep1) == 1))
	    #fired_option2 = (len(goodLeg2objs_lep1) > 0 and len(goodLeg1objs_lep2)>0 and not(goodLeg2objs_lep1 == goodLeg1objs_lep2 and len(goodLeg2objs_lep1) == 1))
	    ##add l1pt cut:
	    fired_option1 = l1ptcut(path, goodLeg1objs_lep1, goodLeg2objs_lep2)
	    fired_option2 = l1ptcut(path, goodLeg1objs_lep2, goodLeg2objs_lep1)
	    
            if fired_option1  or fired_option2:
    		return True

	    #ismatched_leps = [False, False]
	    #selectedObjs = self.matchTriggerObjwithHLT(path, trigobjs)
	    #if self.verbose > 3: 
	    #    print "triggertype ",l1t, " firedpath ", path, " matched l1objs ", len(selectedObjs)
            #if len(selectedObjs) <= 1:
	    #    continue
            #for tobj in selectedObjs:	
	    #   for i in range(0, len(leptons)):
	    #        if tobj.id == abs(leptons[i].pdgId):
	    #    	dR = deltaR(tobj.eta, tobj.phi, leptons[i].eta, leptons[i].phi)
	    #    	dPtRel = abs(tobj.pt-leptons[i].pt)/leptons[i].pt 
	    #    	if self.verbose > 3:
	    #    	    print "lepton pdgId ",leptons[i].pdgId," obj id ",tobj.id, " dR ",dR, " dPtRel ",dPtRel
	    #    	ismatched_leps[i] = ismatched_leps[i] or (dR < self.deltaR_trigger_reco and dPtRel < self.deltaPtRel_trigger_reco)
	    #if ismatched_leps[0] and ismatched_leps[1]:
	    #    return True
	        
	return False


    def trigger_reco_matching(self, muons, trigobjs, leptons_mu):
	""" match reco object to triggger obj """
	for imu, mu in enumerate(muons):
	    triggermatch = False
	    for tobj in trigobjs:
		if abs(mu.pdgId) != tobj.id:
		    continue
		dR = deltaR(tobj.eta, tobj.phi, mu.eta, mu.phi)
		dPtRel = abs(tobj.pt-mu.pt)/mu.pt 
		##check id, dR, dpt
		if dR < self.deltaR_trigger_reco and dPtRel < self.deltaPtRel_trigger_reco:
		    if self.verbose > 3:
			print "trigger obj l1 pt ", tobj.l1pt," l2pt ", tobj.l2pt," pt ",tobj.pt, " eta ",tobj.eta, " bits ", tobj.filterBits, " offlep pt ",mu.pt," eta ",mu.eta," id ",mu.pdgId
		    triggermatch = True
		    break;
	    if triggermatch:
		leptons_mu.append(mu)


    def fillCutFlow(self, cutflow_bin, weight):
	#cutflow_bin += 1.0
	if self.verbose > 4:
	    print "final filling cutflow ",cutflow_bin," triggertype ",self.triggertype," weight ",weight
	self.h_cutflow_weight.Fill(weight)##used to mornitor the weight
	while cutflow_bin > 0:
	    self.h_cutflow.Fill( cutflow_bin, weight)
	    self.h_cutflowlist[self.triggertype].Fill( cutflow_bin, weight)
	    cutflow_bin = cutflow_bin - 1

    def fillCutFlow_leptons(self, cutflow_bin, weight, leptonpairs):
	allL1triggertype = []
	if self.isMC:
	    allL1triggertype = [self.findTriggerType(x) for x in leptonpairs]
	else:
	    allL1triggertype = [self.triggertype]
	
	#only fill for first case
	#print "allL1triggertype ",allL1triggertype
	#l1t = allL1triggertype[0]
	while cutflow_bin > 0:
	    self.h_cutflow.Fill( cutflow_bin, weight)
            for l1t in allL1triggertype:
		self.h_cutflowlist[l1t].Fill( cutflow_bin, weight)
	    cutflow_bin = cutflow_bin - 1
	

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
	#hlts = Collection(event, "HLT")

        self.h_eventcounter.Fill(1)	
        ### PV	
	PV = Object(event, "PV")
	event_pu_weight = 1.0; event_pu_weight_up = 1.0; event_pu_weight_down = 1.0
	LHEPdfweight_stddev = 0.0; LHEScaleWeight = None
	pu = PV.npvsGood
        run = getattr(event,"run", False)	
	luminosityBlock = getattr(event, "luminosityBlock", False)
	ievent =  getattr(event,"event", False)
	RhoFastjetCentralCalo = getattr(event, "fixedGridRhoFastjetCentralCalo", 0.0)
	event_reco_weight = 1.0 ## for pu_weight*btag_SF*lepSF
	sample_weight = 1.0
	genweight = 1.0
	if self.isMC:
	    event_pu_weight = event.puWeight
	    event_pu_weight_up = event.puWeightUp/event_pu_weight
	    event_pu_weight_down = event.puWeightDown/event_pu_weight
	    genweight = getattr(event, "genWeight", 1)
	    sample_weight = genweight * event.puWeight
	    
	    if self.addSystematic and hasattr(event, "LHEPdfWeight"):
		LHEPdfWeight   = event.LHEPdfWeight
	        LHEPdfWeight_sum = sum(LHEPdfWeight)
		nLHEPdfWeight = len(LHEPdfWeight)
		LHEPdfWeight_ave = LHEPdfWeight_sum/nLHEPdfWeight
		LHEPdfweight_tmplist = [(x - LHEPdfWeight_ave)*(x-LHEPdfWeight_ave) for x in LHEPdfWeight]
		LHEPdfweight_stddev = sqrt(sum(LHEPdfweight_tmplist)/(len(LHEPdfWeight)-1.0))
    		#print "sum ",sum(LHEPdfWeight)," ave ",LHEPdfWeight_ave," LHEPdfweight_stddev ",LHEPdfweight_stddev
	    if self.addSystematic and hasattr(event, "LHEScaleWeight"):
		LHEScaleWeight = list(event.LHEScaleWeight)
    		#print "LHEScaleWeight ",LHEScaleWeight, " type ", type(LHEScaleWeight)

	if self.verbose > 1:
	    print "run ",run," luminosityBlock ",luminosityBlock," ievent ",ievent," sample_weight ",sample_weight
	cutflow_bin = 0
	##to get even weight sum
	##first time to fill cutflow histogram at the begining 
	self.h_cutflowlist["DoubleMuon"].Fill( cutflow_bin, sample_weight )
	self.h_cutflowlist["DoubleEG"].Fill( cutflow_bin, sample_weight )
	self.h_cutflowlist["MuonEG"].Fill( cutflow_bin, sample_weight )
	self.h_cutflow.Fill( cutflow_bin, sample_weight)

	### MET
        met = Object(event, "MET")
	metPt = met.pt; metPhi = met.phi
	if hasattr(event, "MET_pt_nom"):
	    metPt = getattr(event, "MET_pt_nom")
	    metPhi = getattr(event, "MET_phi_nom")
	    #print "use calibrated MET pt ",metPt," phi ",metPhi," old pt ",met.pt," phi ",met.phi

        electrons = list(Collection(event, "Electron"))
        muons = list(Collection(event, "Muon"))
	muon_pt_corrected = None 
	if hasattr(event, "Muon_pt_corrected"):
	    muon_pt_corrected = getattr(event, "Muon_pt_corrected")
	    for muon in muons:
	        imu = muons.index(muon)
    		#print "imu ",imu, " pt ",muon.pt," corrected pt ",muon_pt_corrected[imu]
		muon.pt = muon_pt_corrected[imu]
	
        jets = list(Collection(event, "Jet"))
	mht = Object(event, "MHT")
	jet_btagSF = None; jet_btagSF_up = None; jet_btagSF_down = None
	if self.isMC:
	    jet_btagSF      = event.Jet_btagSF
	    jet_btagSF_up   = event.Jet_btagSF_up
	    jet_btagSF_down = event.Jet_btagSF_down
	

        #####################################
        ## di-leps selection
        #####################################

        #muons.sort(key=lambda x:x.pt,reverse=True)	
        #electrons.sort(key=lambda x:x.pt,reverse=True)	


	###cutstep: dilepton pt, and eta
	for mu in muons:
	    if  self.verbose > 3:
	        print "Muon id ",mu.pdgId, " pt ",mu.pt," eta ",mu.eta
	lepton_muons = list(filter(lambda x : abs(x.eta)<self.muonEta, muons))
   
	for el in electrons:
	    if  self.verbose > 3:
	        print "Electron id ", el.pdgId, " pt ", el.pt," eta ", el.eta
	## check eta and whether it passes the conversion veto(not from photon conversion), convVeto = True: good electron
	lepton_electrons = list(filter(lambda x :  abs(x.eta)<self.EGEta, electrons))

        leptonpairs = self.findingLeptonPairs(lepton_muons, lepton_electrons)
	leptonpairs.sort(key=lambda x:x[0].pt+x[1].pt,reverse=True)	
        #ptsumlist = [x[0].pt+x[1].pt for x in leptonpairs]
	#print "ptsumlist after sorting ", ptsumlist

        passdileptonPtEta = (len(leptonpairs) >= 1)
	if not passdileptonPtEta:
	    if self.verbose > 3:
	       print "cutflow_bin ",cutflow_bin," failed , weight ",event_reco_weight * sample_weight," dileptonPtEta"
	    ##cutflow filled already at the beginning 
	    return False
        if self.isMC:
	    self.triggertype = self.findTriggerType(leptonpairs[0]) ## temparory one  for MC
	cutflow_bin += 1


       
	###cutstep: dilepton dz, dxy
	leptonpairs = [x for x in leptonpairs if POGRecipesRun2.leptonpairImpactParameter(x)]
	passdileptonIP = (len(leptonpairs) >= 1)
	if not passdileptonIP:
	    if self.verbose > 3:
	       print "cutflow_bin ",cutflow_bin," failed , weight ",event_reco_weight * sample_weight," dileptonIP "
	    self.fillCutFlow(cutflow_bin, event_reco_weight * sample_weight)
	    #self.fillCutFlow_leptons(cutflow_bin, event_reco_weight * sample_weight, leptonpairs)
	    return False
        if self.isMC: ##update event type, 
	    self.triggertype = self.findTriggerType(leptonpairs[0]) ## temparory one  for MC
	cutflow_bin += 1

		
	###cutstep: dilepton, ID
	leptonpairs = [x for x in leptonpairs if POGRecipesRun2.leptonpairID(x)]
        passdileptonID = (len(leptonpairs) >= 1)
	if not passdileptonID:
	    if self.verbose > 3:
	       print "cutflow_bin ",cutflow_bin," failed , weight ",event_reco_weight * sample_weight," dileptonID"
	    self.fillCutFlow(cutflow_bin, event_reco_weight * sample_weight)
	    #self.fillCutFlow_leptons(cutflow_bin, event_reco_weight * sample_weight, leptonpairs)
	    return False
	cutflow_bin += 1
	llIDsf = 1.0; llIsosf = 1.0; lltrackingsf, lltrgsf, llHLTsafeIDsf  = 1.0, 1.0, 1.0
	########
	##for lep uncertainty
	#######
	l1Isosf       = [1.0, 1.0, 1.0]; l2Isosf       = [1.0, 1.0, 1.0]
	l1IDsf        = [1.0, 1.0, 1.0]; l2IDsf        = [1.0, 1.0, 1.0] 
	l1trgsf       = [1.0, 1.0, 1.0]; l2trgsf       = [1.0, 1.0, 1.0] 
        l1trackingsf  = [1.0, 1.0, 1.0]; l2trackingsf  = [1.0, 1.0, 1.0] 
	l1HLTsafeIDsf = [1.0, 1.0, 1.0]; l2HLTsafeIDsf = [1.0, 1.0, 1.0] 
	if self.isMC: ##may need to update the SFs for trigging since leading lepton pair  may change
	    self.triggertype = self.findTriggerType(leptonpairs[0]) ## temparory one  for MC
	    l1IDsf, l2IDsf                         = self.lepSFmanager.getleptonpairIDSF(leptonpairs[0])
            #print "l1IDsf  ",l1IDsf, " l2IDsf ",l2IDsf
	    llIDsf  = l1IDsf[0] * l2IDsf[0]
	    #llIDsf, llIDsf_up, llIDsf_low = self.lepSFmanager.getleptonpairIDSF(leptonpairs[0])


	###cutstep: dilepton, ISO
	leptonpairs = [x for x in leptonpairs if POGRecipesRun2.leptonpairIso(x)]
        passdileptonIso = (len(leptonpairs) >= 1)
	if not passdileptonIso:
	    if self.verbose > 3:
	       print "cutflow_bin ",cutflow_bin," failed , weight ",event_reco_weight * sample_weight," dileptonIso"
	    if self.isMC:
		event_reco_weight = event_reco_weight * llIDsf## already pass the ID cut, apply ID SF
	    self.fillCutFlow(cutflow_bin, event_reco_weight * sample_weight)
	    #self.fillCutFlow_leptons(cutflow_bin, event_reco_weight * sample_weight, leptonpairs)
	    return False
	cutflow_bin += 1
	if self.isMC:
	    self.triggertype = self.findTriggerType(leptonpairs[0]) ## temparory one  for MC
	    l1Isosf, l2Isosf                       = self.lepSFmanager.getleptonpairIsoSF(leptonpairs[0])
	    l1IDsf, l2IDsf                         = self.lepSFmanager.getleptonpairIDSF(leptonpairs[0])
            llIsosf = l1Isosf[0] * l2Isosf[0]
	    llIDsf  = l1IDsf[0] * l2IDsf[0]
	    #llIsosf, llIsosf_up, llIsosf_low = self.lepSFmanager.getleptonpairIsoSF(leptonpairs[0])
	    #llIDsf, llIDsf_up, llIDsf_low = self.lepSFmanager.getleptonpairIDSF(leptonpairs[0])

	###cutstep: dilepton, HLTSafeID
	leptonpairs = [x for x in leptonpairs if POGRecipesRun2.leptonpairHLTSafeID(x, RhoFastjetCentralCalo)]
        passdileptonHLTSafeID = (len(leptonpairs) >= 1)
	if not passdileptonHLTSafeID:
	    if self.verbose > 3:
	       print "cutflow_bin ",cutflow_bin," failed , weight ",event_reco_weight * sample_weight," dileptonHLTSafeID "
	    if self.isMC:
		event_reco_weight = event_reco_weight * llIsosf * llIDsf  ## already pass the ID+Iso cut, apply ll ID+Iso SF
	    self.fillCutFlow(cutflow_bin, event_reco_weight * sample_weight)
	    #self.fillCutFlow_leptons(cutflow_bin, event_reco_weight * sample_weight, leptonpairs)
	    return False
	cutflow_bin += 1

	
        ### Trigger matching here: data, 
	### we apply trigger SFs to MC and do not cut on trigger matching 
	if not self.isMC:
	    trigobjs = list(Collection(event, "TrigObj"))
	    leptonpairs = [ x for x in leptonpairs if  self.matchHLTPath(event, trigobjs, x)]
	    if len(leptonpairs) == 0:
		if self.verbose > 3:
		   print "cutflow_bin ",cutflow_bin," failed, HTL matching for real data only "
		self.fillCutFlow(cutflow_bin, event_reco_weight * sample_weight)
		return False
	cutflow_bin += 1

	### we should NOT see muon pairs in same CSC region with ETMF bug in real data, if we saw it, for safety, kill it 
	### kill the muon pairs in same CSC region, either both in overlap or both in "non-overlap", for real datea, run < 278167
	### efficiency correction applied for MC is in self.lepSFmanager.getleptonpairTrg
	if not self.isMC and self.triggertype == "DoubleMuon" and run < 277166:
	    leptonpairs = [ x for x in leptonpairs if not(POGRecipesRun2.isMuonPairSameCSCRegion(x[0], x[1])) ]
	    if len(leptonpairs) == 0:
		if self.verbose > 3:
		   print "cutflow_bin ",cutflow_bin," failed , dimuon in same CSC region and due to EMTF bug , run ", run
		self.fillCutFlow(cutflow_bin, event_reco_weight * sample_weight)
		return False
	cutflow_bin += 1



	def dileptonMass(leptons):
	    lep1_p4 = ROOT.TLorentzVector()
	    lep1_p4.SetPtEtaPhiM(leptons[0].pt, leptons[0].eta, leptons[0].phi, leptons[0].mass)
	    lep2_p4 = ROOT.TLorentzVector()
	    lep2_p4.SetPtEtaPhiM(leptons[1].pt, leptons[1].eta, leptons[1].phi, leptons[1].mass)
	    ll_p4 = lep1_p4 + lep2_p4
	    return ll_p4.M()

	###cut: ll_M > 12
	leptonpairs = [x for x in leptonpairs if dileptonMass(x) >12 ]
	passdileptonLowMass = (len(leptonpairs) >= 1)
	if not(passdileptonLowMass):
	    if self.verbose > 3:
	       print "cutflow_bin ",cutflow_bin," failed , weight ",event_reco_weight * sample_weight," llM >12 "
	    self.fillCutFlow(cutflow_bin, event_reco_weight * sample_weight)
	    #self.fillCutFlow_leptons(cutflow_bin, event_reco_weight * sample_weight, leptonpairs)
	    return False
	cutflow_bin += 1


	### not apply right now
	### add lepton pair veto 
	#if len(leptonpairs) > 1:
	#    if self.verbose > 3:
	#       print "cutflow_bin ",cutflow_bin," failed , weight ",event_reco_weight * sample_weight," dileptonpair>1 vetoed "
	#    self.fillCutFlow(cutflow_bin, event_reco_weight * sample_weight)
	#    return False
	#cutflow_bin += 1
	    
        ## select final two leptons, sort the pair by pt sum again 
	leptonpairs.sort(key=lambda x:x[0].pt+x[1].pt, reverse=True)	
	leptons = leptonpairs[0]##leading lepton first
	event_lep_weight = 1.0
	if self.isMC:
	    #llIsosf, llIsosf_up, llIsosf_low = self.lepSFmanager.getleptonpairIsoSF(leptonpairs[0])
	    #llIDsf, llIDsf_up, llIDsf_low = self.lepSFmanager.getleptonpairIDSF(leptonpairs[0])
	    #lltrgsf, lltrgsf_up, lltrgsf_low = self.lepSFmanager.getleptonpairTrgSF(leptonpairs[0])
	    #lltrackingsf, lltrackingsf_up, lltrackingsf_low = self.lepSFmanager.getleptonpairTrackingSF(leptonpairs[0])
	    #llHLTsafeIDsf, llHLTsafeIDsf_up, llHLTsafeIDsf_low = self.lepSFmanager.getleptonpairHTLSafeIDSF(leptonpairs[0])
	    l1Isosf, l2Isosf                       = self.lepSFmanager.getleptonpairIsoSF(leptonpairs[0])
	    l1IDsf, l2IDsf                         = self.lepSFmanager.getleptonpairIDSF(leptonpairs[0])
	    l1trgsf, l2trgsf                       = self.lepSFmanager.getleptonpairTrgSF(leptonpairs[0])
	    l1trackingsf, l2trackingsf             = self.lepSFmanager.getleptonpairTrackingSF(leptonpairs[0])
	    l1HLTsafeIDsf, l2HLTsafeIDsf           = self.lepSFmanager.getleptonpairHTLSafeIDSF(leptonpairs[0])
            llIsosf = l1Isosf[0] * l2Isosf[0]
	    llIDsf  = l1IDsf[0] * l2IDsf[0]
	    lltrgsf  = l1trgsf[0] * l2trgsf[0]
	    lltrackingsf  = l1trackingsf[0] * l2trackingsf[0]
	    llHLTsafeIDsf  = l1HLTsafeIDsf[0] * l2HLTsafeIDsf[0]
    	    event_lep_weight = llIsosf * llIDsf * lltrgsf * llHLTsafeIDsf * lltrackingsf
	    ##finished the lepton selection, add lep weight to event_reco_weight
	    event_reco_weight = event_reco_weight * event_lep_weight
	    #################
	    self.lep1_branchname_sysvalue["Isosf"]        = l1Isosf
	    self.lep1_branchname_sysvalue["IDsf"]         = l1IDsf
	    self.lep1_branchname_sysvalue["trgsf"]        = l1trgsf
	    self.lep1_branchname_sysvalue["HLTsafeIDsf"]  = l1HLTsafeIDsf
	    self.lep1_branchname_sysvalue["trackingsf"]   = l1trackingsf
	    self.lep2_branchname_sysvalue["Isosf"]        = l2Isosf
	    self.lep2_branchname_sysvalue["IDsf"]         = l2IDsf
	    self.lep2_branchname_sysvalue["trgsf"]        = l2trgsf
	    self.lep2_branchname_sysvalue["HLTsafeIDsf"]  = l2HLTsafeIDsf
	    self.lep2_branchname_sysvalue["trackingsf"]   = l2trackingsf


        isMuMu = 0; isMuEl   = 0; isElMu = 0; isElEl = 0; isSF = 0;
	if abs(leptons[0].pdgId) == 13 and abs(leptons[1].pdgId) == 13:
	    isMuMu = 1
	    isSF = 1
	elif abs(leptons[0].pdgId) == 11 and  abs(leptons[1].pdgId) == 11:
	    isElEl = 1
	    isSF = 1
	elif  abs(leptons[0].pdgId) == 13 and  abs(leptons[1].pdgId) == 11:
	    isMuEl = 1
	    isSF = 0
	elif  abs(leptons[0].pdgId) == 11 and  abs(leptons[1].pdgId) == 13:
	    isElMu = 1
	    isSF = 0 


	lep1_iso = 0.0; lep2_iso = 0.0
	lep1_id = 0.0; lep2_id = 0.0
	if abs(leptons[0].pdgId) == 13:
	    lep1_iso = leptons[0].pfRelIso04_all
	    lep1_id = 1 + leptons[0].softId + leptons[0].mediumId + leptons[0].tightId
	else:
	    lep1_iso = leptons[0].pfRelIso03_all
	    lep1_id = leptons[0].cutBased

	if abs(leptons[1].pdgId) == 13:
	    lep2_iso = leptons[1].pfRelIso04_all
	    lep2_id = 1 + leptons[1].softId + leptons[1].mediumId + leptons[1].tightId
	else:
	    lep2_iso = leptons[1].pfRelIso03_all
	    lep2_id = leptons[1].cutBased

        # SetPtEtaPhiE(pt,eta,phi,e); and SetPtEtaPhiM(pt,eta,phi,m);
        lep1_p4 = ROOT.TLorentzVector()
        lep1_p4.SetPtEtaPhiM(leptons[0].pt, leptons[0].eta, leptons[0].phi, leptons[0].mass)
        lep2_p4 = ROOT.TLorentzVector()
        lep2_p4.SetPtEtaPhiM(leptons[1].pt, leptons[1].eta, leptons[1].phi, leptons[1].mass)
        ll_p4 = ROOT.TLorentzVector()
        ll_p4 = lep1_p4 + lep2_p4
	ll_pt = ll_p4.Pt(); ll_M = ll_p4.M(); ll_eta = ll_p4.Eta(); ll_phi = ll_p4.Phi()
        ll_DR_l_l = lep1_p4.DeltaR(lep2_p4)
	ll_DPhi_l_l = deltaPhi(lep1_p4.Phi(), lep2_p4.Phi())
	ll_DEta_l_l = lep1_p4.Eta() - lep2_p4.Eta()
        
        #####################################
        ## di-jets selection
        #####################################
	if not (len(jets) >= 2):
	    if self.verbose > 3:
	       print "cutflow_bin ",cutflow_bin," failed , weight ",event_reco_weight * sample_weight," njets ",len(jets)
	    self.fillCutFlow(cutflow_bin, event_reco_weight * sample_weight)
	    #self.fillCutFlow_leptons(cutflow_bin, event_reco_weight * sample_weight, [leptons])
	    return False
	cutflow_bin += 1

	### jet recalibration
	if hasattr(event, "Jet_pt_nom"):
	    jet_pt_corrected = getattr(event, "Jet_pt_nom")
	    for jet in jets:
                ijet = jets.index(jet)
                #print "Jet pt from reCalibration : ",jet_pt_corrected[ijet]," old pt ",jet.pt
                jets[ijet].pt = jet_pt_corrected[ijet]
	
	###cut: JetPt, JetEta
        bjets = [x for x in jets if x.puId>0 and x.jetId>0 and x.pt>self.jetPt and abs(x.eta)<self.jetEta]
        #bjets = [x for x in jets if x.puId>0 and x.pt>self.jetPt and abs(x.eta)<self.jetEta] ## remove jetid > 0
	if not(len(bjets) >= 2) :
	    if self.verbose > 3:
	       print "cutflow_bin ",cutflow_bin," failed , weight ",event_reco_weight * sample_weight," jetPtEta "
	    self.fillCutFlow(cutflow_bin, event_reco_weight * sample_weight)
	    #self.fillCutFlow_leptons(cutflow_bin, event_reco_weight * sample_weight, [leptons])
	    return False
	cutflow_bin += 1

	###cut: the seperation between jet and lepton
        #bjets.sort(key=lambda x:x.pt,reverse=True)	
    	deltaRl1j = lambda jet :  deltaR(jet.eta, jet.phi, leptons[0].eta, leptons[0].phi) > self.deltaR_j_l
    	deltaRl2j = lambda jet :  deltaR(jet.eta, jet.phi, leptons[1].eta, leptons[1].phi) > self.deltaR_j_l
        bjets = [x for x in bjets if deltaRl1j(x) and deltaRl2j(x)]
	nJetsL = len(bjets)
        if not(len(bjets) >= 2): 
	    if self.verbose > 3:
	       print "cutflow_bin ",cutflow_bin," failed , weight ",event_reco_weight * sample_weight," bl seperation"
	    self.fillCutFlow(cutflow_bin, event_reco_weight * sample_weight)
	    #self.fillCutFlow_leptons(cutflow_bin, event_reco_weight * sample_weight, [leptons])
	    return False
	cutflow_bin += 1

	#### select final two jets: jets with maximum pT , for DY estimation 
        ## sort alljets with a cMVAv2 descreasing order 
        bjets.sort(key=lambda x:x.btagCMVA, reverse=True)	
	ht_jets = sum([x.pt for x in bjets])
	hJets = bjets[0:2]
	hJets.sort(key=lambda x:x.pt, reverse=True)
	hJidx = [jets.index(x) for x in hJets]
	jet1 = hJets[0]; jet2 = hJets[1]

	##include all jets before btagging for eff measurement 
	##only select jets up to self.maxnjets
	def resize_alljets(thislist):
	    if len(thislist)> self.maxnjets:	 
    		return thislist[0: self.maxnjets] 
	    else: 
	    	return thislist+[0]*(self.maxnjets-len(alljets))
	alljets = bjets
	alljets_pt = [jet.pt for jet in alljets]
	alljets_eta = [jet.eta for jet in alljets]
	alljets_cMVAv2 = [jet.btagCMVA for jet in alljets]
	flavour = lambda x : abs(x)*(abs(x) == 4 or abs(x) == 5)
        alljets_pt = resize_alljets(alljets_pt)
        alljets_eta = resize_alljets(alljets_eta)
        alljets_cMVAv2 = resize_alljets(alljets_cMVAv2)
	
            		
	###cut: jet btagging , normal selection, skip it for DY estimation
	if not self.DYestimation:
	    bjets = [x for x in bjets if POGRecipesRun2.jetMediumBtagging(x)]
	    if not(len(bjets) >= 2): 
		if self.verbose > 3:
		   print "cutflow_bin ",cutflow_bin," failed , weight ",event_reco_weight * sample_weight, " btagging "
		self.fillCutFlow(cutflow_bin, event_reco_weight * sample_weight)
		#self.fillCutFlow_leptons(cutflow_bin, event_reco_weight * sample_weight, [leptons])
		return False
	    cutflow_bin += 1

	    #### select final two bjets: jets with maximum btagging  and then fill cutflow hist
	    hJets = sorted(bjets, key = lambda jet : jet.btagCMVA, reverse=True)[0:2]
	    hJets.sort(key=lambda x:x.pt, reverse=True)
	    hJidx = [jets.index(x) for x in hJets]
	    jet1 = hJets[0]; jet2 = hJets[1]

	### add btagging weight 
	hJets_BtagSF      = [1.0, 1.0]
	hJets_BtagSF_up   = [1.0, 1.0]
	hJets_BtagSF_down = [1.0, 1.0]

	if self.isMC:
	    hJets_BtagSF       = [jet_btagSF[x] for x in hJidx]
	    hJets_BtagSF_up    = [jet_btagSF_up[x] for x in hJidx]
	    hJets_BtagSF_down  = [jet_btagSF_down[x] for x in hJidx]
	    ## FIXME, how to apply SFs:  https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods
	    event_reco_weight = event_reco_weight * hJets_BtagSF[0] * hJets_BtagSF[1]

	if ll_M < 76:
	    cutflow_bin += 1
	else:
	    self.fillCutFlow(cutflow_bin, event_reco_weight * sample_weight)
	    #self.fillCutFlow_leptons(cutflow_bin, event_reco_weight * sample_weight, [leptons])
	    if not self.DYestimation:##only keep ll_M >76 for DY estimation
		return False
	self.fillCutFlow(cutflow_bin, event_reco_weight * sample_weight)
	if self.verbose > 3:
	    print "cutflow_bin ",cutflow_bin," Fine , weight ",event_reco_weight * sample_weight


        ## Save a few basic reco. H kinematics
        hj1_p4 = ROOT.TLorentzVector()
        hj2_p4 = ROOT.TLorentzVector()
        hj1_p4.SetPtEtaPhiM(jets[hJidx[0]].pt,jets[hJidx[0]].eta,jets[hJidx[0]].phi,jets[hJidx[0]].mass)
        hj2_p4.SetPtEtaPhiM(jets[hJidx[1]].pt,jets[hJidx[1]].eta,jets[hJidx[1]].phi,jets[hJidx[1]].mass)
        hbb_p4 = hj1_p4 + hj2_p4
	jj_pt = hbb_p4.Pt(); jj_M = hbb_p4.M(); jj_eta = hbb_p4.Eta(); jj_phi = hbb_p4.Phi()
        jj_DR_j_j = hj1_p4.DeltaR(hj2_p4)
	jj_DPhi_j_j = deltaPhi(hj1_p4.Phi(), hj2_p4.Phi())
        jj_DEta_j_j = hj1_p4.Eta() - hj2_p4.Eta()

        met_p4 = ROOT.TLorentzVector()## mass=0, eta=0
        met_p4.SetPtEtaPhiM(metPt,0.,metPhi,0.) # only use met vector to derive transverse quantities
	dR_b1l1 = hj1_p4.DeltaR(lep1_p4)
	dR_b2l1 = hj2_p4.DeltaR(lep1_p4)
	dR_b1l2 = hj1_p4.DeltaR(lep2_p4)
	dR_b2l2 = hj2_p4.DeltaR(lep2_p4)
        llmetjj_minDR_l_j = min([dR_b1l1, dR_b2l1, dR_b1l2, dR_b2l2])
        llmetjj_DR_ll_jj = ll_p4.DeltaR(hbb_p4)
        llmetjj_DPhi_ll_met =  deltaPhi(ll_p4.Phi(), met_p4.Phi())
        llmetjj_DPhi_ll_jj =  deltaPhi(ll_p4.Phi(), hbb_p4.Phi())
	llmetjj_MTformula = sqrt(2*ll_p4.Pt()*met_p4.Pt()*(1-cos(llmetjj_DPhi_ll_met)))
        lljj_p4 = hbb_p4 + ll_p4
	lljj_M = lljj_p4.M()
        llmet_p4 = met_p4 + ll_p4
	llmetjj_p4 = met_p4 + ll_p4 + hbb_p4
	llmetjj_M = llmetjj_p4.M()
        llmet_M = llmet_p4.M()
        llmetjj_DPhi_llmet_jj =  deltaPhi(llmet_p4.Phi(), hbb_p4.Phi())
    
        cosThetaStar = cos(llmet_p4.Theta() - hbb_p4.Theta())
    	ht = ht_jets + lep1_p4.Pt() + lep2_p4.Pt() 
        ###

        #jet selection, old
        #jetsForHiggs = [x for x in jets if x.lepFilter and x.puId>0 and x.jetId>0 and x.Pt>20 and abs(x.eta)<2.5]
        #if (len(jetsForHiggs) < 2): return False
        self.out.fillBranch("jet1_idx",  hJidx[0]);
        self.out.fillBranch("jet1_pt",  hj1_p4.Pt());
        self.out.fillBranch("jet1_E",  hj1_p4.E());
        self.out.fillBranch("jet1_eta",  hj1_p4.Eta());
        self.out.fillBranch("jet1_phi",  hj1_p4.Phi());
        self.out.fillBranch("jet1_cMVAv2",  jet1.btagCMVA);
        self.out.fillBranch("jet2_idx",  hJidx[1]);
        self.out.fillBranch("jet2_pt", hj2_p4.Pt());
        self.out.fillBranch("jet2_E",  hj2_p4.E());
        self.out.fillBranch("jet2_eta", hj2_p4.Eta());
        self.out.fillBranch("jet2_phi",  hj2_p4.Phi());
        self.out.fillBranch("jet2_cMVAv2",  jet2.btagCMVA);
	self.out.fillBranch("isElEl", isElEl) # 0 or 1
	self.out.fillBranch("isElMu",  isElMu) # 0 or 1, mu_pt<el_pt
	self.out.fillBranch("isMuEl",  isMuEl) # 0 or 1
	self.out.fillBranch("isMuMu",  isMuMu) # 0 or 1
	self.out.fillBranch("isSF",  isSF)
	#self.out.fillBranch("lepstype", Lepstype)
        self.out.fillBranch("lep1_pt",  lep1_p4.Pt());
        self.out.fillBranch("lep1_E",  lep1_p4.E());
        self.out.fillBranch("lep1_eta",  lep1_p4.Eta());
        self.out.fillBranch("lep1_phi",  lep1_p4.Phi());
        self.out.fillBranch("lep1_pdgid", leptons[0].pdgId );
        self.out.fillBranch("lep1_id",  lep1_id);
        self.out.fillBranch("lep1_iso",  lep1_iso);
        self.out.fillBranch("lep2_pt",  lep2_p4.Pt());
        self.out.fillBranch("lep2_E",   lep2_p4.E());
        self.out.fillBranch("lep2_eta",  lep2_p4.Eta());
        self.out.fillBranch("lep2_phi",  lep2_p4.Phi());
        self.out.fillBranch("lep2_pdgid",  leptons[1].pdgId);
        self.out.fillBranch("lep2_id",  lep2_id);
        self.out.fillBranch("lep2_iso",  lep2_iso);
	self.out.fillBranch("met_pt",  met_p4.Pt())
	self.out.fillBranch("met_phi", met_p4.Phi())
	self.out.fillBranch("ht",  ht)
	self.out.fillBranch("ht_jets",  ht_jets)
	self.out.fillBranch("nJetsL",  nJetsL)
	self.out.fillBranch("jj_M",  jj_M)
	self.out.fillBranch("el_hltsafeid", llHLTsafeIDsf)
	self.out.fillBranch("llreco_weight",  lltrackingsf)
	self.out.fillBranch("llid_weight",  llIDsf)
	self.out.fillBranch("lliso_weight",  llIsosf)
	self.out.fillBranch("lltrigger_weight",  lltrgsf)
	self.out.fillBranch("ll_M",  ll_M)
	self.out.fillBranch("llmet_M",  llmet_M)
	#self.out.fillBranch("llmetjj_MT2",  "F")
	self.out.fillBranch("llmetjj_M", llmetjj_M)
	self.out.fillBranch("lljj_M",  lljj_M)
	self.out.fillBranch("cosThetaStar",  cosThetaStar)
	self.out.fillBranch("ll_DR_l_l",  ll_DR_l_l)
	self.out.fillBranch("jj_DR_j_j",  jj_DR_j_j)
	self.out.fillBranch("llmetjj_DPhi_ll_jj",  llmetjj_DPhi_ll_jj)
	self.out.fillBranch("llmetjj_DPhi_ll_met",  llmetjj_DPhi_ll_met)
	self.out.fillBranch("llmetjj_DPhi_llmet_jj",  llmetjj_DPhi_llmet_jj)
	self.out.fillBranch("ll_pt",  ll_pt)
	self.out.fillBranch("ll_eta",  ll_eta)
	self.out.fillBranch("jj_pt",  jj_pt)
	self.out.fillBranch("jj_eta",  jj_eta)
	self.out.fillBranch("llmetjj_minDR_l_j",  llmetjj_minDR_l_j)
	self.out.fillBranch("llmetjj_MTformula",  llmetjj_MTformula)
	self.out.fillBranch("ll_DPhi_l_l",  ll_DPhi_l_l)
	self.out.fillBranch("pu",  pu)
	self.out.fillBranch("event_pu_weight",  event_pu_weight)
	self.out.fillBranch("event_lep_weight",  event_lep_weight)
	self.out.fillBranch("event_btag_weight",  hJets_BtagSF[0] * hJets_BtagSF[1])
        self.out.fillBranch("sample_weight",  sample_weight)
        self.out.fillBranch("genweight",  genweight)
	self.out.fillBranch("event_reco_weight",  event_reco_weight)
	self.out.fillBranch("total_weight",  event_reco_weight * sample_weight)
	self.out.fillBranch("event_run",  run)
	self.out.fillBranch("event_lumiblock",  luminosityBlock)
	#self.out.fillBranch("event_weight",  "F")
	#self.out.fillBranch("DY_BDT_flat",  "F")
	#self.out.fillBranch("dy_nobtag_to_btagM_weight",  "F")
	#self.out.fillBranch("mt2",  "F")
	#self.out.fillBranch("mt2_bb",  "F")
	#self.out.fillBranch("mt2_ll",  "F")
	#self.out.fillBranch("event_number",  ievent)
	if self.isMC and self.addSystematic:
	    for name in self.ll_sys_branchname:
		self.out.fillBranch("lep1%s_up"%name,     self.lep1_branchname_sysvalue[name][1])
		self.out.fillBranch("lep1%s_down"%name,   self.lep1_branchname_sysvalue[name][2])
		self.out.fillBranch("lep2%s_up"%name,     self.lep2_branchname_sysvalue[name][1])
		self.out.fillBranch("lep2%s_down"%name,   self.lep2_branchname_sysvalue[name][2])
	    self.out.fillBranch("event_pu_weight_up",    event_pu_weight_up)
	    self.out.fillBranch("event_pu_weight_down",  event_pu_weight_down)
	    self.out.fillBranch("event_pdf_weight_up",   1.0 + LHEPdfweight_stddev)
	    self.out.fillBranch("event_pdf_weight_down", 1.0 - LHEPdfweight_stddev)
	    #self.out.fillBranch("LHEScaleWeight", LHEScaleWeight)
	    #print "len(LHEScaleWeight) ",len(LHEScaleWeight)
	    for i in range(0, len(LHEScaleWeight)):
		self.out.fillBranch("LHEScaleWeight_%d"%i, LHEScaleWeight[i])

	if self.isMC:
	    genmet = Object(event, "GenMET")
	    genParticles = Collection(event, "GenPart")
	    genjets = Collection(event, "GenJet")
	    nGenPart = getattr(event, "nGenPart", False)
	    nGenJet = getattr(event, "nGenJet", False)

	    nGenW1FromHiggs=0; nGenW2FromHiggs=0; nGenLepFromW1FromHiggs=0; nGenLepFromW2FromHiggs=0; nGenNuFromW1FromHiggs =0; nGenNuFromW2FromHiggs=0
	    nGenBQuarkFromHiggs=0
	    if hasattr(event, "nGenX"):
	        #nGenX   =  getattr(event, "nGenX", 0) 
		nGenW1FromHiggs = getattr(event, "nGenW1FromHiggs", 0) 
		nGenW2FromHiggs = getattr(event, "nGenW2FromHiggs", 0) 
		nGenLepFromW1FromHiggs = getattr(event, "nGenLepFromW1FromHiggs", 0)
		nGenNuFromW1FromHiggs = getattr(event, "nGenNuFromW1FromHiggs", 0)
		nGenLepFromW2FromHiggs = getattr(event, "nGenLepFromW2FromHiggs", 0)
		nGenNuFromW2FromHiggs = getattr(event, "nGenNuFromW2FromHiggs", 0)
		nGenBQuarkFromHiggs = getattr(event, "nGenBQuarkFromHiggs", 0)
	    if nGenW1FromHiggs  == 1 and nGenW2FromHiggs  == 1:
	        GenW1FromHiggs_mass = getattr(event, "GenW1FromHiggs_mass")
	        GenW1FromHiggs_pt = getattr(event, "GenW1FromHiggs_pt")
	        GenW2FromHiggs_mass = getattr(event, "GenW2FromHiggs_mass")
	        GenW2FromHiggs_pt = getattr(event, "GenW2FromHiggs_pt")
	        #self.out.fillBranch("genW1_pt",    GenW1FromHiggs_pt[0])
	        #self.out.fillBranch("genW1_mass",  GenW1FromHiggs_mass[0])
	        #self.out.fillBranch("genW2_pt",    GenW2FromHiggs_pt[0])
	        #self.out.fillBranch("genW2_mass",  GenW2FromHiggs_mass[0])
		if self.verbose > 3:
		  print "Gen W1 from Higgs, mass ", GenW1FromHiggs_mass[0], " pt ",GenW1FromHiggs_pt[0],\
		    "Gen W2 from Higgs, mass ",GenW2FromHiggs_mass[0], " pt ",GenW2FromHiggs_pt[0]
	    genNuFromW1_p4 = ROOT.TLorentzVector()
	    genNuFromW2_p4 = ROOT.TLorentzVector()
	    genLepFromW1_p4 = ROOT.TLorentzVector()
	    genLepFromW2_p4 = ROOT.TLorentzVector()
            dR_genl1_lepFromW1 = 99.0            
            dR_genl1_lepFromW2 = 99.0            
            dR_genl2_lepFromW1 = 99.0            
            dR_genl2_lepFromW2 = 99.0            
            foundLepsfromWs = False
	    if nGenLepFromW1FromHiggs == 1 and nGenNuFromW1FromHiggs == 1 and nGenLepFromW2FromHiggs == 1 and nGenNuFromW2FromHiggs == 1:
	        foundLepsfromWs = True
	        GenLepFromW1FromHiggs_pt   = getattr(event, "GenLepFromW1FromHiggs_pt")
	        GenLepFromW1FromHiggs_eta  = getattr(event, "GenLepFromW1FromHiggs_eta")
	        GenLepFromW1FromHiggs_phi  = getattr(event, "GenLepFromW1FromHiggs_phi")
	        GenNuFromW1FromHiggs_pt   = getattr(event, "GenNuFromW1FromHiggs_pt")
	        GenNuFromW1FromHiggs_eta  = getattr(event, "GenNuFromW1FromHiggs_eta")
	        GenNuFromW1FromHiggs_phi  = getattr(event, "GenNuFromW1FromHiggs_phi")
	        GenLepFromW2FromHiggs_pt   = getattr(event, "GenLepFromW2FromHiggs_pt")
	        GenLepFromW2FromHiggs_eta  = getattr(event, "GenLepFromW2FromHiggs_eta")
	        GenLepFromW2FromHiggs_phi  = getattr(event, "GenLepFromW2FromHiggs_phi")
	        GenNuFromW2FromHiggs_pt   = getattr(event, "GenNuFromW2FromHiggs_pt")
	        GenNuFromW2FromHiggs_eta  = getattr(event, "GenNuFromW2FromHiggs_eta")
	        GenNuFromW2FromHiggs_phi  = getattr(event, "GenNuFromW2FromHiggs_phi")
		genLepFromW1_p4.SetPtEtaPhiM(GenLepFromW1FromHiggs_pt[0], GenLepFromW1FromHiggs_eta[0], GenLepFromW1FromHiggs_phi[0], 0.0)
		genNuFromW1_p4.SetPtEtaPhiM(GenNuFromW1FromHiggs_pt[0], GenNuFromW1FromHiggs_eta[0], GenNuFromW1FromHiggs_phi[0], 0.0)
		genLepFromW2_p4.SetPtEtaPhiM(GenLepFromW2FromHiggs_pt[0], GenLepFromW2FromHiggs_eta[0], GenLepFromW2FromHiggs_phi[0], 0.0)
		genNuFromW2_p4.SetPtEtaPhiM(GenNuFromW2FromHiggs_pt[0], GenNuFromW2FromHiggs_eta[0], GenNuFromW2FromHiggs_phi[0], 0.0)
                diNu_p4 = genNuFromW1_p4+genNuFromW2_p4
		met_diNu_pt = diNu_p4.Pt()
		met_diNu_phi = diNu_p4.Phi()
		if  self.addGenToTree:
		    self.out.fillBranch("met_diNuetrino_pt",  met_diNu_pt);
		    self.out.fillBranch("met_diNuetrino_phi",  met_diNu_phi);
		    #self.out.fillBranch("gennu1_pt",  GenNuFromW1FromHiggs_pt[0]) #nu1 from W1
		    #self.out.fillBranch("gennu1_eta", GenNuFromW1FromHiggs_eta[0]) 
		    #self.out.fillBranch("gennu1_phi", GenNuFromW1FromHiggs_phi[0]) 
		    #self.out.fillBranch("gennu2_pt",  GenNuFromW2FromHiggs_pt[0]) 
		    #self.out.fillBranch("gennu2_eta", GenNuFromW2FromHiggs_eta[0]) 
		    #self.out.fillBranch("gennu2_phi", GenNuFromW2FromHiggs_phi[0]) 
		if self.verbose > 3:
                  print "GenLepFromW1 pt ",GenLepFromW1FromHiggs_pt[0], " eta ",GenLepFromW1FromHiggs_eta[0], " phi ",GenLepFromW1FromHiggs_phi[0],\
		      " GenLepFromW2 pt ",GenLepFromW2FromHiggs_pt[0], " eta ",GenLepFromW2FromHiggs_eta[0], " phi ",GenLepFromW2FromHiggs_phi[0]
		  print "MET from two nuetrinos pt ",met_diNu_pt, " phi ",met_diNu_phi
		  print "genmet_pt ",genmet.pt," phi ",genmet.phi

	    genb1FromHiggs_p4 = ROOT.TLorentzVector()
	    genb2FromHiggs_p4 = ROOT.TLorentzVector()
            dR_genb1_genjet1 = 99.0
            dR_genb1_genjet2 = 99.0
            dR_genb2_genjet1 = 99.0
            dR_genb2_genjet2 = 99.0
            foundBQuarkFromHiggs = False
	    if nGenBQuarkFromHiggs == 2 :
	        foundBQuarkFromHiggs = True
	        GenBQuarkFromHiggs_pt = getattr(event, "GenBQuarkFromHiggs_pt")
	        GenBQuarkFromHiggs_mass = getattr(event, "GenBQuarkFromHiggs_mass")
	        GenBQuarkFromHiggs_eta = getattr(event, "GenBQuarkFromHiggs_eta")
	        GenBQuarkFromHiggs_phi = getattr(event, "GenBQuarkFromHiggs_phi")
		genb1FromHiggs_p4.SetPtEtaPhiM(GenBQuarkFromHiggs_pt[0], GenBQuarkFromHiggs_eta[0], GenBQuarkFromHiggs_phi[0], GenBQuarkFromHiggs_mass[0])
		genb2FromHiggs_p4.SetPtEtaPhiM(GenBQuarkFromHiggs_pt[1], GenBQuarkFromHiggs_eta[1], GenBQuarkFromHiggs_phi[1], GenBQuarkFromHiggs_mass[1])
		if  self.addGenToTree:
		    self.out.fillBranch("genb1_pt",  GenBQuarkFromHiggs_pt[0]) #b1
		    self.out.fillBranch("genb2_pt",  GenBQuarkFromHiggs_pt[1]) #b2
		    #self.out.fillBranch("genb1_eta",  GenBQuarkFromHiggs_eta[0]) #b1
		    #self.out.fillBranch("genb1_phi",  GenBQuarkFromHiggs_phi[0]) #b1
		    #self.out.fillBranch("genb2_eta",  GenBQuarkFromHiggs_eta[1]) #b2
		    #self.out.fillBranch("genb2_phi",  GenBQuarkFromHiggs_phi[1]) #b2
		if self.verbose > 3:
		    print "GEn B quark from Higgs pt ",GenBQuarkFromHiggs_pt, " eta ",GenBQuarkFromHiggs_eta, " phi ",GenBQuarkFromHiggs_phi
	    
	    lep1_genindex = leptons[0].genPartIdx
	    lep2_genindex = leptons[1].genPartIdx
	    if lep1_genindex >= 0 and lep1_genindex < nGenPart and self.addGenToTree:
		genl1 = genParticles[lep1_genindex]
		if self.verbose > 3:
		    print "lep1 pt ",genl1.pt," eta ",genl1.eta, " phi ",genl1.phi
		genl1_p4 = ROOT.TLorentzVector(); genl1_p4.SetPtEtaPhiM(genl1.pt, genl1.eta, genl1.phi, genl1.mass)
		self.out.fillBranch("genl1_pt",  genl1_p4.Pt());
		self.out.fillBranch("genl1_E",  genl1_p4.E());
		self.out.fillBranch("genl1_eta",  genl1_p4.Eta());
		self.out.fillBranch("genl1_phi",  genl1_p4.Phi());
		self.out.fillBranch("genl1_id",  abs(genl1.pdgId));
		if foundLepsfromWs:
		    dR_genl1_lepFromW1 = genl1_p4.DeltaR(genLepFromW1_p4)
		    dR_genl1_lepFromW2 = genl1_p4.DeltaR(genLepFromW2_p4)
	    if lep2_genindex >= 0 and lep2_genindex < nGenPart and self.addGenToTree:
		genl2 = genParticles[lep2_genindex]
		if self.verbose > 3:
		    print "lep2 pt ",genl2.pt," eta ",genl2.eta, " phi ",genl2.phi
		genl2_p4 = ROOT.TLorentzVector(); genl2_p4.SetPtEtaPhiM(genl2.pt, genl2.eta, genl2.phi, genl2.mass)
		self.out.fillBranch("genl2_pt",  genl2_p4.Pt());
		self.out.fillBranch("genl2_E",  genl2_p4.E());
		self.out.fillBranch("genl2_eta",  genl2_p4.Eta());
		self.out.fillBranch("genl2_phi",  genl2_p4.Phi());
		self.out.fillBranch("genl2_id",  abs(genl2.pdgId));
		if foundLepsfromWs:
		    dR_genl2_lepFromW1 = genl2_p4.DeltaR(genLepFromW1_p4)
		    dR_genl2_lepFromW2 = genl2_p4.DeltaR(genLepFromW2_p4)
	    ## jet matching  
            genjet1_partonflavour = 0; genjet2_partonflavour = 0
	    jjbtag_heavy = 1.0; jjbtag_light = 1.0
    	    if jet1.genJetIdx >= 0 and jet1.genJetIdx < nGenJet:
		genjet1 = genjets[jet1.genJetIdx]
		genjet1_p4 = ROOT.TLorentzVector(); 
		if self.verbose > 3:
		    print "genjet1 pt ",genjet1.pt," eta ",genjet1.eta, " phi ",genjet1.phi
		genjet1_p4.SetPtEtaPhiM(genjet1.pt, genjet1.eta, genjet1.phi, genjet1.mass)
		genjet1_partonflavour = flavour(genjet1.partonFlavour)
    		if self.addGenToTree:
		    self.out.fillBranch("genjet1_pt", genjet1_p4.Pt());
		    self.out.fillBranch("genjet1_E",  genjet1_p4.E());
		    self.out.fillBranch("genjet1_eta",  genjet1_p4.Eta());
		    self.out.fillBranch("genjet1_phi",  genjet1_p4.Phi());
		if self.DYestimation or self.addGenToTree:	
		    self.out.fillBranch("genjet1_partonFlavour",  genjet1_partonflavour);
		if foundBQuarkFromHiggs:
		    dR_genb1_genjet1 = genjet1_p4.DeltaR(genb1FromHiggs_p4)
		    dR_genb2_genjet1 = genjet1_p4.DeltaR(genb2FromHiggs_p4)
		    
    	    if jet2.genJetIdx >= 0 and jet2.genJetIdx < nGenJet :
		genjet2 = genjets[jet2.genJetIdx]
		genjet2_p4 = ROOT.TLorentzVector(); 
		if self.verbose > 3:
		    print "genjet2 pt ",genjet2.pt," eta ",genjet2.eta, " phi ",genjet2.phi
		genjet2_partonflavour = flavour(genjet2.partonFlavour)
		genjet2_p4.SetPtEtaPhiM(genjet2.pt, genjet2.eta, genjet2.phi, genjet2.mass)
    		if self.addGenToTree:
		    self.out.fillBranch("genjet2_pt", genjet2_p4.Pt());
		    self.out.fillBranch("genjet2_E",  genjet2_p4.E());
		    self.out.fillBranch("genjet2_eta",  genjet2_p4.Eta());
		    self.out.fillBranch("genjet2_phi",  genjet2_p4.Phi());
		if self.DYestimation or self.addGenToTree:	
		    self.out.fillBranch("genjet2_partonFlavour",  genjet2_partonflavour);

		if foundBQuarkFromHiggs:
		    dR_genb1_genjet2 = genjet2_p4.DeltaR(genb1FromHiggs_p4)
		    dR_genb2_genjet2 = genjet2_p4.DeltaR(genb2FromHiggs_p4) 
	    #self.jjbtag_sys_branchname = ["jjbtag_light","jjbtag_heavy"]
            self.jjbtag_branchname_sysvalue["jjbtag_heavy"] = [1.0, 1.0, 1.0]
            self.jjbtag_branchname_sysvalue["jjbtag_light"] = [1.0, 1.0, 1.0]
	    if genjet1_partonflavour > 0:
	        jjbtag_heavy = jjbtag_heavy * hJets_BtagSF[0] 
		self.jjbtag_branchname_sysvalue["jjbtag_heavy"][0] = jjbtag_heavy
		self.jjbtag_branchname_sysvalue["jjbtag_heavy"][1] = self.jjbtag_branchname_sysvalue["jjbtag_heavy"][1]* hJets_BtagSF_up[0]/hJets_BtagSF[0]
		self.jjbtag_branchname_sysvalue["jjbtag_heavy"][2] = self.jjbtag_branchname_sysvalue["jjbtag_heavy"][2]* hJets_BtagSF_down[0]/hJets_BtagSF[0]
	    else:
	        jjbtag_light = jjbtag_light * hJets_BtagSF[0] 
		self.jjbtag_branchname_sysvalue["jjbtag_light"][0] = jjbtag_light
		self.jjbtag_branchname_sysvalue["jjbtag_light"][1] = self.jjbtag_branchname_sysvalue["jjbtag_light"][1]* hJets_BtagSF_up[0]/hJets_BtagSF[0]
		self.jjbtag_branchname_sysvalue["jjbtag_light"][2] = self.jjbtag_branchname_sysvalue["jjbtag_light"][2]* hJets_BtagSF_down[0]/hJets_BtagSF[0]
	    if genjet2_partonflavour > 0:
	        jjbtag_heavy = jjbtag_heavy * hJets_BtagSF[1] 
		self.jjbtag_branchname_sysvalue["jjbtag_heavy"][0] = jjbtag_heavy
		self.jjbtag_branchname_sysvalue["jjbtag_heavy"][1] = self.jjbtag_branchname_sysvalue["jjbtag_heavy"][1]* hJets_BtagSF_up[1]/hJets_BtagSF[1]
		self.jjbtag_branchname_sysvalue["jjbtag_heavy"][2] = self.jjbtag_branchname_sysvalue["jjbtag_heavy"][2]* hJets_BtagSF_down[1]/hJets_BtagSF[1]
	    else:
	        jjbtag_light = jjbtag_light * hJets_BtagSF[1] 
		self.jjbtag_branchname_sysvalue["jjbtag_light"][0] = jjbtag_light
		self.jjbtag_branchname_sysvalue["jjbtag_light"][1] = self.jjbtag_branchname_sysvalue["jjbtag_light"][1]* hJets_BtagSF_up[1]/hJets_BtagSF[1]
		self.jjbtag_branchname_sysvalue["jjbtag_light"][2] = self.jjbtag_branchname_sysvalue["jjbtag_light"][2]* hJets_BtagSF_down[1]/hJets_BtagSF[1]
	    self.out.fillBranch("jjbtag_heavy", jjbtag_heavy);
	    self.out.fillBranch("jjbtag_light", jjbtag_light);
	    if self.addSystematic:
		for name in self.jjbtag_branchname_sysvalue.keys():
		    self.out.fillBranch(name+"_up",   self.jjbtag_branchname_sysvalue[name][1])
		    self.out.fillBranch(name+"_down", self.jjbtag_branchname_sysvalue[name][2])


	    #genflavour = lambda jet : (jet.genJetIdx < nGenJet and jet.genJetIdx >=0)*flavour(genjets[abs(jet.genJetIdx)].partonFlavour)
            if self.CheckBtaggingEff:
		alljets_genpartonFlavour = []
		for jet in alljets:
		    if jet.genJetIdx < nGenJet and jet.genJetIdx >=0:
			alljets_genpartonFlavour.append(flavour(genjets[abs(jet.genJetIdx)].partonFlavour))
		    else:
			alljets_genpartonFlavour.append(0)
			
		alljets_genpartonFlavour = resize_alljets(alljets_genpartonFlavour)
		alljets_partonFlavour = [jet.partonFlavour for jet in alljets]
		alljets_hadronFlavour = [jet.hadronFlavour for jet in alljets]
		alljets_partonFlavour = resize_alljets(alljets_partonFlavour)
		alljets_hadronFlavour = resize_alljets(alljets_hadronFlavour)
		self.out.fillBranch("alljets_pt", alljets_pt)
		self.out.fillBranch("alljets_eta", alljets_eta)
		self.out.fillBranch("alljets_cMVAv2", alljets_cMVAv2)
		self.out.fillBranch("alljets_genpartonFlavour", alljets_genpartonFlavour)
		#self.out.fillBranch("alljets_partonFlavour", alljets_partonFlavour)
		#self.out.fillBranch("alljets_hadronFlavour", alljets_hadronFlavour)

	    if self.DYestimation or self.addGenToTree:	
		self.out.fillBranch("jet1_partonFlavour",  flavour(jet1.partonFlavour));
		self.out.fillBranch("jet1_hadronFlavour",  flavour(jet1.hadronFlavour));
		self.out.fillBranch("jet2_partonFlavour",  flavour(jet2.partonFlavour));
		self.out.fillBranch("jet2_hadronFlavour",  flavour(jet2.hadronFlavour));
		self.out.fillBranch("genmet_pt",  genmet.pt);
		self.out.fillBranch("genmet_phi", genmet.phi);

    	    if  self.addGenToTree:
		self.out.fillBranch("dR_genb1_genjet1",  dR_genb1_genjet1);
		self.out.fillBranch("dR_genb1_genjet2",  dR_genb1_genjet2);
		self.out.fillBranch("dR_genb2_genjet1",  dR_genb2_genjet1);
		self.out.fillBranch("dR_genb2_genjet2",  dR_genb2_genjet2);
		self.out.fillBranch("dR_genl1_lepFromW1",  dR_genl1_lepFromW1);
		self.out.fillBranch("dR_genl1_lepFromW2",  dR_genl1_lepFromW2);
		self.out.fillBranch("dR_genl2_lepFromW1",  dR_genl2_lepFromW1);
		self.out.fillBranch("dR_genl2_lepFromW2",  dR_genl2_lepFromW2);



        return True
                

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

#hhbbWW = lambda : HHbbWWProducer(True, "") 
#hhbbWW_data = lambda x : HHbbWWProducer(False, x) 
