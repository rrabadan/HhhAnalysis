import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from math import sqrt,cos
import random
#random.randint(0,100)

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection,Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import * #deltaR, matching etc..
import POGRecipesRun2
lepSFmanager = POGRecipesRun2.LeptonSFManager()

class HHbbWWProducer(Module):
    ## data or MC, which L1 trigger, HLT?
    ###kwargs: triggertype, verbose, run_lumi
    def __init__(self,isMC, **kwargs):
        self.isMC = isMC ## two mode: data or MC
	self.triggertype  = ''##"DoubleMuon, DoubleEG, MuonEG"
	self.deltaR_trigger_reco = 0.1; self.deltaPtRel_trigger_reco = 0.5
	self.verbose = 3
	self.muonEta = 2.4; self.EGEta = 2.5
	self.leadingMuonPt = {"DoubleMuon": 20.0, "MuonEG":25.0}
	self.subleadingMuonPt = {"DoubleMuon": 10.0, "MuonEG":10.0}
	self.leadingEGPt = {"DoubleEG": 25.0, "MuonEG":25.0}
	self.subleadingEGPt = {"DoubleEG": 15.0, "MuonEG":15.0}
	self.jetPt = 20; self.jetEta = 2.4
	self.deltaR_j_l = 0.3 #jet,lepton seperation
	self.h_cutflow = ROOT.TH1F("h_cutflow","h_cutflow", 20, 0, 20.0)
	self.h_cutflowlist = {"DoubleMuon":ROOT.TH1F("h_cutflow_DoubleMuon","h_cutflow_DoubleMuon", 20, 0, 20.0), 
			      "DoubleEG":ROOT.TH1F("h_cutflow_DoubleEG","h_cutflow_DoubleEG", 20, 0, 20.0), 
			      "MuonEG":ROOT.TH1F("h_cutflow_MuonEG","h_cutflow_MuonEG", 20, 0, 20.0)}

	self.h_cutflowlist["DoubleMuon"].SetLineColor(ROOT.kRed)
	self.h_cutflowlist["DoubleEG"].SetLineColor(ROOT.kBlue)
	self.h_cutflowlist["MuonEG"].SetLineColor(ROOT.kBlack)

	self.run_lumi = None
	self.__dict__.update((key, kwargs[key]) for key in kwargs.keys())
        print "self.run_lumi  ",self.run_lumi," trigger type ",self.triggertype," verbose ",self.verbose

	### for debug
	self.nEv_DoubleEG = 0
	self.nEv_DoubleMuon = 0
	self.nEv_MuonEG = 0
        pass
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("jet1_pt",  "F");
        self.out.branch("jet1_E",  "F");
        self.out.branch("jet1_eta",  "F");
        self.out.branch("jet1_phi",  "F");
        self.out.branch("jet1_cMVAv2",  "F");
        self.out.branch("jet2_pt",  "F");
        self.out.branch("jet2_E",  "F");
        self.out.branch("jet2_eta",  "F");
        self.out.branch("jet2_phi",  "F");
        self.out.branch("jet2_cMVAv2",  "F");
	self.out.branch("isElEl",  "F") # 0 or 1
	self.out.branch("isElMu",  "F") # 0 or 1, mu_pt<el_pt
	self.out.branch("isMuEl",  "F") # 0 or 1
	self.out.branch("isMuMu",  "F") # 0 or 1
	self.out.branch("isSF",  "F")
	self.out.branch("lepstype",  "F")
        self.out.branch("lep1_pt",  "F");
        self.out.branch("lep1_E",  "F");
        self.out.branch("lep1_eta",  "F");
        self.out.branch("lep1_phi",  "F");
        self.out.branch("lep1_id",  "F");
        self.out.branch("lep1_iso",  "F");
        self.out.branch("lep2_pt",  "F");
        self.out.branch("lep2_E",  "F");
        self.out.branch("lep2_eta",  "F");
        self.out.branch("lep2_phi",  "F");
        self.out.branch("lep2_id",  "F");
        self.out.branch("lep2_iso",  "F");
	self.out.branch("met_pt",  "F")
	self.out.branch("met_phi",  "F")
	self.out.branch("nJetsL",  "F")
	self.out.branch("jjbtag_heavy",  "F")
	self.out.branch("jjbtag_light",  "F")
	self.out.branch("jj_M",  "F")
	self.out.branch("llidiso",  "F")
	self.out.branch("trigeff",  "F")
	self.out.branch("ll_M",  "F")
	self.out.branch("ht",  "F")
	self.out.branch("llmetjj_MT2",  "F")
	self.out.branch("llmetjj_M",  "F")
	self.out.branch("lljj_M",  "F")
	self.out.branch("cosThetaStar",  "F")
	self.out.branch("ll_DR_l_l",  "F")
	self.out.branch("jj_DR_j_j",  "F")
	self.out.branch("llmetjj_DPhi_ll_jj",  "F")
	self.out.branch("ll_pt",  "F")
	self.out.branch("jj_pt",  "F")
	self.out.branch("llmetjj_minDR_l_j",  "F")
	self.out.branch("llmetjj_MTformula",  "F")
	self.out.branch("ll_DPhi_l_l",  "F")
	self.out.branch("ll_DEta_l_l",  "F")
	self.out.branch("event_pu_weight",  "F")
	self.out.branch("event_lep_weight",  "F")
	self.out.branch("event_btag_weight",  "F")
        self.out.branch("sample_weight",  "F");
	self.out.branch("event_reco_weight",  "F")
	self.out.branch("event_weight",  "F")
	self.out.branch("pu",  "F")
	self.out.branch("DY_BDT_flat",  "F")
	self.out.branch("dy_nobtag_to_btagM_weight",  "F")
	self.out.branch("mt2",  "F")
	self.out.branch("mt2_bb",  "F")
	self.out.branch("mt2_ll",  "F")
	self.out.branch("event_number",  "I")
	self.out.branch("event_run",  "I")
	##how to add gen information???
	
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
	self.h_cutflow.Write()
	self.h_cutflowlist["DoubleMuon"].Write()
	self.h_cutflowlist["DoubleEG"].Write()
	self.h_cutflowlist["MuonEG"].Write()
        pass

    def goldenJason(self, run, lumi):
	run_str = '%d'%run
	if run_str in self.run_lumi.keys():
	    alllumis = self.run_lumi[run_str]
	    for lumirange in alllumis:
	        if lumi >= lumirange[0] and lumi <= lumirange[1]:
	            return True
	else:
	    return False


    def findingLeptonPairs(self, leptons_mu, leptons_el):

	leptonpairs = []
	if self.triggertype == "DoubleMuon" and len(leptons_mu) >= 2:
	    nmuons = len(leptons_mu)
	    for imu1 in range(0, nmuons):
		for imu2 in range(imu1+1, nmuons):
		    if leptons_mu[imu1].charge * leptons_mu[imu2].charge < 0 and leptons_mu[imu1].pt > self.leadingMuonPt["DoubleMuon"] and leptons_mu[imu2].pt > self.subleadingMuonPt["DoubleMuon"]:
			leptonpairs.append(( leptons_mu[imu1],  leptons_mu[imu2] ))

	elif self.triggertype == "DoubleEG" and len(leptons_el) >= 2:
	    nelectrons = len(leptons_el)
	    for iel1 in range(0, nelectrons):
		for iel2 in range(iel1+1, nelectrons):
		    if leptons_el[iel1].charge * leptons_el[iel2].charge < 0 and leptons_el[iel1].pt > self.leadingEGPt["DoubleEG"] and leptons_el[iel2].pt > self.subleadingEGPt["DoubleEG"]:
			leptonpairs.append(( leptons_el[iel1],  leptons_el[iel2] ))

	elif self.triggertype == "MuonEG" and len(leptons_mu) >= 1 and len(leptons_el) >= 1:
	    nmuons = len(leptons_mu)
	    nelectrons = len(leptons_el)
	    for imu1 in range(0, nmuons):
		for iel1 in range(0, nelectrons):
		    if leptons_mu[imu1].charge * leptons_el[iel1].charge > 0:
		    	continue
		    if leptons_mu[imu1].pt >= self.leadingMuonPt["MuonEG"] and leptons_el[iel1].pt >= self.subleadingEGPt["MuonEG"]:
			leptonpairs.append(( leptons_mu[imu1],  leptons_el[iel1] ))
		    elif leptons_mu[imu1].pt >= self.subleadingMuonPt["MuonEG"] and leptons_el[iel1].pt >= self.leadingEGPt["MuonEG"]:
			leptonpairs.append(( leptons_el[iel1],  leptons_mu[imu1] ))

	if self.verbose > 3:
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
        ## with data branch MET_[pt/phi]. This essentially aliases the two branches to one common variable.
        if isMC:
            return (met.pt_smeared,met.phi_smeared)
        else:
            return (met.pt,met.phi)
	   
    def matchTriggerObjwithHLT(self, path, trigobjs):
	allhltlist = path.split('_')[1:-1]
	selectedObjs = []
	for x in allhltlist:
	    if "Mu" in x:##TkMu is not inclued in HLT path, FIXME??
		pt = int(x[2:])
		for obj in trigobjs:
		    if obj.id == 13 and obj.pt >= pt and (obj.l1pt >= 4.0 or obj.l1pt_2 >= 4.0):
		    	selectedObjs.append(obj)
	    elif "Ele" in x:
		pt = int(x[3:])
		for obj in trigobjs:
		    if obj.id == 11 and obj.pt >= pt and (obj.l1pt >= 10.0 or obj.l1pt_2 >= 10.0):##mini pt for L1 seed? 10GeV?
		    	selectedObjs.append(obj)

	return set(selectedObjs)

    def HLTPath(self, hlt, trigobjs, isMC):
	""" check which HLT is fired """
	L1t_hlt = {"DoubleMuon": ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ"],
				  #"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ"], ##FIXME
	    	   "MuonEG": ["HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL","HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 
		   		"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL","HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ"],
		   "DoubleEG": ["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"]
			}
        paths_objs = {}
	l1t_objs = {}
        if isMC:
	    for l1 in L1t_hlt:
	    	matchedobjs = []
	    	for path in L1t_hlt[l1]:
		    firepath = getattr(hlt,  path, False)
    		    if firepath:
		        if self.verbose > 4:
			    if l1 == "DoubleEG":
				self.nEv_DoubleEG += 1.0
				print "firedpath ",path," nEv_DoubleEG ",self.nEv_DoubleEG
			    elif l1 == "MuonEG":
				self.nEv_MuonEG += 1.0
				print "firedpath ",path," nEv_MuonEG ",self.nEv_MuonEG
			    elif l1 == "DoubleMuon":
				self.nEv_DoubleMuon += 1.0
				print "firedpath ",path," nEv_DoubleMuon ",self.nEv_DoubleMuon

		    	selectedObjs = self.matchTriggerObjwithHLT(path, trigobjs)
    			if len(selectedObjs) >= 2:
			    #self.triggertype = l1
			    paths_objs[path] = selectedObjs
			    matchedobjs.extend(selectedObjs)
    			    break

	   	if len(matchedobjs)>0 :
		    matchedobjs = set(matchedobjs)
		    l1t_objs[l1] = matchedobjs

	    l1ts = l1t_objs.keys()
	    if len(l1ts) > 1:
	    	randint = random.randint(0, len(l1ts)-1)
    		if self.verbose > 1:
		    print "all fired paths_objs ",paths_objs.keys()," selected trigger type ",l1ts[randint]
    		return True, l1ts[randint], l1t_objs[l1ts[randint]]
	    elif len(l1ts) == 1:
    		return True, l1ts[0], l1t_objs[l1ts[0]]

	else:
	    matchedobjs = []
	    l1 = self.triggertype
	    for path in L1t_hlt[self.triggertype]:
		firepath = getattr(hlt,  path, False)
		if firepath:
		    selectedObjs = self.matchTriggerObjwithHLT(path, trigobjs)
		    if len(selectedObjs) >= 2:
			#paths_objs[path] = selectedObjs
			matchedobjs.extend(selectedObjs)
	    if len(matchedobjs)>0 :
		matchedobjs = set(matchedobjs)
		l1t_objs[l1] = matchedobjs
		return True, l1, matchedobjs
		    #return True, path

	return False," ",[]
    def trigger_reco_matching(self, muons, trigobjs, leptons_mu):
	""" match reco object to triggger obj """
	for imu, mu in enumerate(muons):
	    triggermatch = False
	    for tobj in trigobjs:
		if abs(mu.pdgId) != tobj.id:
		    continue
		dR = deltaR(tobj.eta, tobj.phi, mu.eta, mu.phi)
		dPtRel = abs(tobj.pt-mu.pt)/mu.pt ###which trigger object pt should be used here, FIXME
		##check id, dR, dpt
		if dR < self.deltaR_trigger_reco and dPtRel < self.deltaPtRel_trigger_reco:
		    if self.verbose > 3:
			print "trigger obj l1 pt ", tobj.l1pt," l2pt ", tobj.l2pt," pt ",tobj.pt, " eta ",tobj.eta, " bits ", tobj.filterBits, " offlep pt ",mu.pt," eta ",mu.eta," id ",mu.pdgId
		    triggermatch = True
		    break;
	    if triggermatch:
		leptons_mu.append(mu)


    def fillMCinfo(self, event):
	""" fill all gen information here """
	if not self.isMC:
		return 

	genmet = Object(event)

    def fillCutFlow(self, cutflow_bin, weight):
	#cutflow_bin += 1.0
	if self.verbose > 4:
	    print "filling cutflow ",cutflow_bin," triggertype ",self.triggertype," weight ",weight
	self.h_cutflow.Fill( cutflow_bin, weight)
	self.h_cutflowlist[self.triggertype].Fill( cutflow_bin, weight)

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
	#hlts = Collection(event, "HLT")

	
        ### PV	
	PV = Object(event, "PV")
	event_pu_weight = 1.0
	pu = PV.npvs
        run = getattr(event,"run", False)	
	luminosityBlock = getattr(event, "luminosityBlock", False)
	ievent =  getattr(event,"event", False)
	event_reco_weight = 1.0 ## for pu_weight*btag_SF*lepSF
	sample_weight = 1.0
	if not(self.isMC) and self.run_lumi != None:
            if not(self.goldenJason(run, luminosityBlock)):
		### skip events failed the goldenjason check for real data 
		return False
	if self.isMC:
	    #sample_weight = 1.0
	    event_pu_weight = event.puWeight
	    sample_weight = getattr(event, "genWeight", 1) * event.puWeight
	if self.verbose > 1:
	    print "run ",run," luminosityBlock ",luminosityBlock," ievent ",ievent," sample_weight ",sample_weight
	cutflow_bin = 0.0
	##first time to fill cutflow histogram 
	if self.isMC:
	    self.h_cutflowlist["DoubleMuon"].Fill( cutflow_bin, event_reco_weight * sample_weight *.25)
	    self.h_cutflowlist["DoubleEG"].Fill( cutflow_bin, event_reco_weight * sample_weight * .25)
	    self.h_cutflowlist["MuonEG"].Fill( cutflow_bin, event_reco_weight * sample_weight * .5)
	else:
	    self.h_cutflowlist["DoubleMuon"].Fill( cutflow_bin, event_reco_weight * sample_weight )
	    self.h_cutflowlist["DoubleEG"].Fill( cutflow_bin, event_reco_weight * sample_weight )
	    self.h_cutflowlist["MuonEG"].Fill( cutflow_bin, event_reco_weight * sample_weight )
	self.h_cutflow.Fill( cutflow_bin, event_reco_weight * sample_weight)

	###cutstep1: initial events, no selection, 
	trigobjs = Collection(event, "TrigObj")
	fired,triggertype, matchedobjs =  self.HLTPath(event, trigobjs, self.isMC)
        
	##cutstep2: HLT is fired or not
	if not fired:
	    if self.verbose > 3:
	        print "cutflow_bin ",cutflow_bin," failed , weight ",event_reco_weight * sample_weight
	    return False
	if self.isMC:
	    self.triggertype = triggertype
	    if self.verbose > 2:
		print "Final fired trigger type ",triggertype
	cutflow_bin += 1.0
	self.fillCutFlow(cutflow_bin, event_reco_weight * sample_weight)
	
        Lepstype = -1##MuMu:1, MuEl:2, ElMu:3, ElEl:4
	if self.triggertype == "DoubleMuon":
        	Lepstype = 1	
	elif self.triggertype == "DoubleEG":
		Lepstype = 4
	elif self.triggertype == "MuonEG":
		Lepstype = 2
	#elif self.triggertype == "MuonEG" and "Ele23" in path:
		#Lepstype = 3
	else:
		raise ValueError('Triggertype is not acceptable , check self.triggertype ')
	        

	### MET
        met = Object(event, "MET")
        #metPt,metPhi = self.met(met,self.isMC)
	metPt = met.pt; metPhi = met.phi

        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        jets = list(Collection(event, "Jet"))
	jet_btagSF = 0
	#mu_effSF = 0
	#el_effSF = 0
	if self.isMC:
	    #mu_effSF = event.Muon_effSF
	    #el_effSF = event.Electron_effSF
	    jet_btagSF = event.Jet_btagSF
	

        #####################################
        ## di-leps selection
        #####################################
	###match trigger obj to HLT?
	leptons_mu = []
	leptons_el = []
	
	if self.triggertype == "DoubleMuon" or self.triggertype == "MuonEG":
            self.trigger_reco_matching(muons, matchedobjs, leptons_mu)
	if self.triggertype == "DoubleEG" or self.triggertype == "MuonEG":
            self.trigger_reco_matching(electrons, matchedobjs, leptons_el)
	

        leptons_mu.sort(key=lambda x:x.pt,reverse=True)	
        leptons_el.sort(key=lambda x:x.pt,reverse=True)	
	###cutstep3: L1_HLT_RECO matching
        passdileptonTrig = False
	if (self.triggertype == "DoubleMuon" and len(leptons_mu) >= 2):
	    if self.isMC:
	        lltrgsf, lltrgsf_up, lltrgsf_low = lepSFmanager.getleptonpairTrgSF(leptons_mu[:2])
	        event_reco_weight = lltrgsf * event_reco_weight
	    passdileptonTrig = True
	elif (self.triggertype == "DoubleEG" and len(leptons_el) >= 2):
	    if self.isMC:
	        lltrgsf, lltrgsf_up, lltrgsf_low = lepSFmanager.getleptonpairTrgSF(leptons_el[:2])
	        event_reco_weight = lltrgsf * event_reco_weight
	    passdileptonTrig = True
	elif (self.triggertype == "MuonEG" and len(leptons_mu)>=1 and len(leptons_el) >= 1):
	    if self.isMC:
	        lltrgsf, lltrgsf_up, lltrgsf_low = lepSFmanager.getleptonpairTrgSF([leptons_mu[0], leptons_el[0]])
	        event_reco_weight = lltrgsf * event_reco_weight
	    passdileptonTrig = True
	else:
	    if self.verbose > 3:
	        print " cutflow_bin ",cutflow_bin," failed , weight ",event_reco_weight * sample_weight," Triggermatching "
	    return False

	cutflow_bin += 1.0
	self.fillCutFlow(cutflow_bin, event_reco_weight * sample_weight)
	


	###cutstep4: dilepton pt, and eta, dz, dxy??
	for mu in leptons_mu:
	    if  self.verbose > 3:
	        print "Muon id ",mu.pdgId, " pt ",mu.pt," eta ",mu.eta
	    if abs(mu.eta) >= self.muonEta:
	        leptons_mu.remove(mu)
	for el in leptons_el:
	    if  self.verbose > 3:
	        print "Electron id ", el.pdgId, " pt ", el.pt," eta ", el.eta
	    if abs(el.eta) >= self.EGEta or not(el.convVeto): 
		## check eta and whether it passes the conversion veto(not from photon conversion), convVeto = True: good electron
	        leptons_el.remove(el)

        leptonpairs = self.findingLeptonPairs(leptons_mu, leptons_el)

        for pair in leptonpairs:
	    if not(POGRecipesRun2.leptonpairImpactParameter(pair)):
		leptonpairs.remove(pair)

        passdileptonPtEta = (len(leptonpairs) >= 1)
	if not passdileptonPtEta:
	    if self.verbose > 3:
	       print "cutflow_bin ",cutflow_bin," failed , weight ",event_reco_weight * sample_weight," dileptoPtEta"
	    return False
	cutflow_bin += 1.0
	self.fillCutFlow(cutflow_bin, event_reco_weight * sample_weight)
		
	###cutstep6: dilepton, ID
        for pair in leptonpairs:
	    if not(POGRecipesRun2.leptonpairID(pair)):
		leptonpairs.remove(pair)

        passdileptonID = (len(leptonpairs) >= 1)
	if not passdileptonID:
	    if self.verbose > 3:
	       print "cutflow_bin ",cutflow_bin," failed , weight ",event_reco_weight * sample_weight," dileptonID"
	    return False
	cutflow_bin += 1.0
	if self.isMC:
	    llIDsf, llIDsf_up, llIDsf_low = lepSFmanager.getleptonpairIDSF(leptonpairs[0])
	    event_reco_weight = event_reco_weight*llIDsf
	self.fillCutFlow(cutflow_bin, event_reco_weight * sample_weight)


	###cutstep5: dilepton, ISO
        for pair in leptonpairs:
	    if not(POGRecipesRun2.leptonpairIso(pair)):
		leptonpairs.remove(pair)

        passdileptonIso = (len(leptonpairs) >= 1)
	if not passdileptonIso:
	    if self.verbose > 3:
	       print "cutflow_bin ",cutflow_bin," failed , weight ",event_reco_weight * sample_weight," dileptonIso"
	    return False
	cutflow_bin += 1.0
	if self.isMC:
	    llIsosf, llIsosf_up, llIsosf_low = lepSFmanager.getleptonpairIsoSF(leptonpairs[0])
	    event_reco_weight = event_reco_weight*llIsosf
	self.fillCutFlow(cutflow_bin, event_reco_weight * sample_weight)

	###cutstep6: dilepton, HLTSafeID
        for pair in leptonpairs:
	    if not(POGRecipesRun2.leptonpairHLTSafeID(pair)):
		leptonpairs.remove(pair)

        passdileptonHLTSafeID = (len(leptonpairs) >= 1)
	if not passdileptonHLTSafeID:
	    if self.verbose > 3:
	       print "cutflow_bin ",cutflow_bin," failed , weight ",event_reco_weight * sample_weight," dileptonHLTSafeID "
	    return False
	cutflow_bin += 1.0
	self.fillCutFlow(cutflow_bin, event_reco_weight * sample_weight)

        ## select final two leptons	
	leptons = leptonpairs[0]##leading lepton first
	if self.triggertype == "MuonEG":
	    if abs(leptons[0].pdgId) == 13:
    		Lepstype = 2
    	    else:
    		Lepstype = 3

        # SetPtEtaPhiE(pt,eta,phi,e); and SetPtEtaPhiM(pt,eta,phi,m);
        lep1_p4 = ROOT.TLorentzVector()
        lep1_p4.SetPtEtaPhiM(leptons[0].pt, leptons[0].eta, leptons[0].phi, leptons[0].mass)
        lep2_p4 = ROOT.TLorentzVector()
        lep2_p4.SetPtEtaPhiM(leptons[1].pt, leptons[1].eta, leptons[1].phi, leptons[1].mass)
        ll_p4 = ROOT.TLorentzVector()
        ll_p4 = lep1_p4 + lep2_p4
	ll_pt = ll_p4.Pt(); ll_M = ll_p4.M(); ll_Eta = ll_p4.Eta(); ll_Phi = ll_p4.Phi()
        ll_DR_l_l = lep1_p4.DeltaR(lep2_p4)
	ll_DPhi_l_l = deltaPhi(lep1_p4.Phi(), lep2_p4.Phi())
	ll_DEta_l_l = lep1_p4.Eta() - lep2_p4.Eta()
        
	###cut: ll_M > 12
	if not(ll_p4.M() > 12):
	    if self.verbose > 3:
	       print "cutflow_bin ",cutflow_bin," failed , weight ",event_reco_weight * sample_weight," llM ",ll_p4.M()
	    return False
	cutflow_bin += 1.0
	self.fillCutFlow(cutflow_bin, event_reco_weight * sample_weight)


        #####################################
        ## di-jets selection
        #####################################
	if len(jets) >= 2:
	    cutflow_bin += 1.0
	    self.fillCutFlow(cutflow_bin, event_reco_weight * sample_weight)
	else:
	    if self.verbose > 3:
	       print "cutflow_bin ",cutflow_bin," failed , weight ",event_reco_weight * sample_weight," njets ",len(jets)
	    return False

	###cut: JetPt, JetEta
	ht = sum([x.pt for x in jets])
        bjets = [x for x in jets if x.puId>0 and x.jetId>0 and x.pt>self.jetPt and abs(x.eta)<self.jetEta]
	if not(len(bjets) >= 2) :
	    if self.verbose > 3:
	       print "cutflow_bin ",cutflow_bin," failed , weight ",event_reco_weight * sample_weight," jetPtEta "
	    return False
            		
	cutflow_bin += 1.0
	self.fillCutFlow(cutflow_bin, event_reco_weight * sample_weight)
	###cut: the seperation between jet and lepton
        bjets.sort(key=lambda x:x.pt,reverse=True)	
        for jet in bjets:
	    dR1 = deltaR(jet.eta, jet.phi, leptons[0].eta, leptons[0].phi)
	    dR2 = deltaR(jet.eta, jet.phi, leptons[1].eta, leptons[1].phi)
            if dR1 < self.deltaR_j_l or dR2 < self.deltaR_j_l:
	        bjets.remove(jet)
	nJetsL = len(bjets)
        if not(len(bjets) >= 2): 
	    if self.verbose > 3:
	       print "cutflow_bin ",cutflow_bin," failed , weight ",event_reco_weight * sample_weight," bl seperation"
	    return False
	cutflow_bin += 1.0
	self.fillCutFlow(cutflow_bin, event_reco_weight * sample_weight)
            		
	###cut: jet btagging 
        for jet in bjets:
	    if not (POGRecipesRun2.jetMediumBtagging(jet)):
		bjets.remove(jet)
        if not(len(bjets) >= 2): 
	    if self.verbose > 3:
	       print "cutflow_bin ",cutflow_bin," failed , weight ",event_reco_weight * sample_weight, " btagging "
	    return False
		
        #### select final two bjets: jets with maximum btagging  and then fill cutflow hist
        hJets = sorted(bjets, key = lambda jet : jet.btagCMVA, reverse=True)[0:2]
        hJidx = [jets.index(x) for x in hJets]
	jet1 = hJets[0]; jet2 = hJets[1]
	hJets_BtagSF = [1.0, 1.0]

	if self.isMC:
	    hJets_BtagSF = [jet_btagSF[x] for x in hJidx]
	    event_reco_weight = event_reco_weight * hJets_BtagSF[0] * hJets_BtagSF[1]
	cutflow_bin += 1.0
	self.fillCutFlow(cutflow_bin, event_reco_weight * sample_weight)
	if self.verbose > 3:
	    print "cutflow_bin ",cutflow_bin," Fine , weight ",event_reco_weight * sample_weight


        ## Save a few basic reco. H kinematics
        hj1_p4 = ROOT.TLorentzVector()
        hj2_p4 = ROOT.TLorentzVector()
        hj1_p4.SetPtEtaPhiM(jets[hJidx[0]].pt,jets[hJidx[0]].eta,jets[hJidx[0]].phi,jets[hJidx[0]].mass)
        hj2_p4.SetPtEtaPhiM(jets[hJidx[1]].pt,jets[hJidx[1]].eta,jets[hJidx[1]].phi,jets[hJidx[1]].mass)
        hbb_p4 = hj1_p4 + hj2_p4
	jj_pt = hbb_p4.Pt(); jj_M = hbb_p4.M(); jj_Eta = hbb_p4.Eta(); jj_Phi = hbb_p4.Phi()
        jj_DR_j_j = hj1_p4.DeltaR(hj2_p4)
	jj_DPhi_j_j = deltaPhi(hj1_p4.Phi(), hj2_p4.Phi())
        jj_DEta_j_j = hj1_p4.Eta() - hj2_p4.Eta()

        met_p4 = ROOT.TLorentzVector()## mass=0, eta=0
        met_p4.SetPtEtaPhiM(met.pt,0.,met.phi,0.) # only use met vector to derive transverse quantities
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
        ###

        #jetsForHiggs = [x for x in jets if x.lepFilter and x.puId>0 and x.jetId>0 and x.Pt>20 and abs(x.eta)<2.5]
        #if (len(jetsForHiggs) < 2): return False
        isMuMu = 0; isMuEl   = 0; isElMu = 0; isElEl = 0; isSF = 0;
	if Lepstype == 1:
	    isMuMu = 1
	    isSF = 1
	elif Lepstype == 2:
	    isMuEl = 1
	    isSF = 0
	elif Lepstype == 3:
	    isElMu = 1
	    isSF = 0
	else:
	    isElEl = 1
	    isSF = 1 

        self.out.fillBranch("jet1_pt",  hj1_p4.Pt());
        self.out.fillBranch("jet1_E",  hj1_p4.E());
        self.out.fillBranch("jet1_eta",  hj1_p4.Eta());
        self.out.fillBranch("jet1_phi",  hj1_p4.Phi());
        self.out.fillBranch("jet1_cMVAv2",  jet1.btagCMVA);
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
	self.out.fillBranch("lepstype", Lepstype)
        self.out.fillBranch("lep1_pt",  lep1_p4.Pt());
        self.out.fillBranch("lep1_E",  lep1_p4.E());
        self.out.fillBranch("lep1_eta",  lep1_p4.Eta());
        self.out.fillBranch("lep1_phi",  lep1_p4.Phi());
        #self.out.fillBranch("lep1_id", "F" );
        #self.out.fillBranch("lep1_iso",  "F");
        self.out.fillBranch("lep2_pt",  lep2_p4.Pt());
        self.out.fillBranch("lep2_E",   lep2_p4.E());
        self.out.fillBranch("lep2_eta",  lep2_p4.Eta());
        self.out.fillBranch("lep2_phi",  lep2_p4.Phi());
        #self.out.fillBranch("lep2_id",  "F");
        #self.out.fillBranch("lep2_iso",  "F");
	self.out.fillBranch("met_pt",  met_p4.Pt())
	self.out.fillBranch("met_phi", met_p4.Phi())
	self.out.fillBranch("ht",  ht)
	self.out.fillBranch("nJetsL",  nJetsL)
	#self.out.fillBranch("jjbtag_heavy",  "F")
	#self.out.fillBranch("jjbtag_light",  "F")
	self.out.fillBranch("jj_M",  jj_M)
	#self.out.fillBranch("llidiso",  "F")
	#self.out.fillBranch("trigeff",  "F")
	self.out.fillBranch("ll_M",  ll_M)
	#self.out.fillBranch("llmetjj_MT2",  "F")
	#self.out.fillBranch("llmetjj_M",  "F")
	self.out.fillBranch("lljj_M",  lljj_M)
	#self.out.fillBranch("cosThetaStar",  "F")
	self.out.fillBranch("ll_DR_l_l",  ll_DR_l_l)
	self.out.fillBranch("jj_DR_j_j",  jj_DR_j_j)
	self.out.fillBranch("llmetjj_DPhi_ll_jj",  llmetjj_DPhi_ll_jj)
	self.out.fillBranch("ll_pt",  ll_pt)
	self.out.fillBranch("jj_pt",  jj_pt)
	self.out.fillBranch("llmetjj_minDR_l_j",  llmetjj_minDR_l_j)
	self.out.fillBranch("llmetjj_MTformula",  llmetjj_MTformula)
	self.out.fillBranch("ll_DPhi_l_l",  ll_DPhi_l_l)
	self.out.fillBranch("pu",  pu)
	self.out.fillBranch("event_pu_weight",  event_pu_weight)
	#self.out.fillBranch("event_lep_weight",  "F")
	self.out.fillBranch("event_btag_weight",  hJets_BtagSF[0] * hJets_BtagSF[1])
        self.out.fillBranch("sample_weight",  sample_weight)
	self.out.fillBranch("event_reco_weight",  event_reco_weight)
	#self.out.fillBranch("event_weight",  "F")
	#self.out.fillBranch("DY_BDT_flat",  "F")
	#self.out.fillBranch("dy_nobtag_to_btagM_weight",  "F")
	#self.out.fillBranch("mt2",  "F")
	#self.out.fillBranch("mt2_bb",  "F")
	#self.out.fillBranch("mt2_ll",  "F")
	#self.out.fillBranch("event_number",  "I")
	#self.out.fillBranch("event_run",  "I")
        #self.out.fillBranch("H_pt",hbb.Pt())
        #self.out.fillBranch("H_phi",hbb.Phi())
        #self.out.fillBranch("H_eta",hbb.Eta())
        #self.out.fillBranch("H_mass",hbb.M())

	### Compute soft activity vetoing Higgs jets
	##find signal footprint
	#matchedSAJets=self.matchSoftActivity(hJets,sa)
	## update SA variables 


	#softActivityJetHT=event.SoftActivityJetHT2-sum([x.pt for x in matchedSAJets])
        #self.out.fillBranch("SA_Ht",softActivityJetHT)

	#matchedSAJetsPt5=[x for x in matchedSAJets if x.pt>5]
        #softActivityJetNjets5=event.SoftActivityJetNjets5-len(matchedSAJetsPt5)
        #self.out.fillBranch("SA5",softActivityJetNjets5)

        return True
                

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

#hhbbWW = lambda : HHbbWWProducer(True, "") 
#hhbbWW_data = lambda x : HHbbWWProducer(False, x) 
