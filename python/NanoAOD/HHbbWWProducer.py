import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection,Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import * #deltaR, matching etc..

class HHbbWWProducer(Module):
    ## data or MC, which L1 trigger, HLT?
    def __init__(self, isMC, L1trigger, verbose = 4):
        self.isMC = isMC ## two mode: data or MC
	self.triggertype  = L1trigger##"DoubleMuon, DoubleEG, MuonEG"
	self.deltaR_trigger_reco = 0.1; self.deltaPtRel_trigger_reco = 0.5
	self.ievent = 0
	self.verbose = verbose
	self.leadingMuonPt = {"DoubleMuon": 20; "MuonEG":25}
	self.subleadingMuonPt = {"DoubleMuon": 10; "MuonEG":10}
	self.leadingEGPt = {"DoubleEG": 25; "MuonEG":25}
	self.subleadingEGPt = {"DoubleEG": 15; "MuonEG":15}
	self.h_cutflow = ROOT.TH1F("h_cutflow","h_cutflow", 20, 0, 20.0)

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
        pass

    def electronImpactParameterCut(self, electron):
	""" check electron impact parameter, https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2 """
	dxy_endcap = 0.1 #cm
	dxy_barrel = 0.05
	dz_endcap  = 0.2
	dz_barrel = 0.1
	if abs(electron.eta) < 1.479: ##barrel
	    return (abs(electron.dz) < dz_barrel and abs(electron.dxy) < dxy_barrel)
        else abs(electron.eta) >= 1.479 and abs(electron.eta) < 2.5:
	    return (abs(electron.dz) < dz_endcap and abs(electron.dxy) < dxy_endcap)
	else:
	    return False

    def electronID(self, electron):
	""" check electron ID, https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2 """
	### we use medium Electron
	return electron.cutBased >= 3
	
    def electronIso(self, electron):
	""" check electron isolation, https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2 """
	### electron iso cut
	return True
    
    def electronHLTSafeID(self, electron):
	""" check electron isolation, https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2 """
	"""https://twiki.cern.ch/twiki/bin/view/CMS/ChangesEGMHLTAlgo2014"""
	### extrac HLT safe ID cut is applied for electron selection for 2016 data only
	return True
	#if abs(electron.eta) < 1.479: #barrel
	#return (abs(electron.sieie) < 0.11 and 
	
    def muonID(self, muon):
        """https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#Muon_Isolation """
	#goodGlob = (muon.isGlobalMuon() and muon.globalTrack()->normalizedChi2() < 3 and muon.combinedQuality().chi2LocalPosition < 12 and muon.combinedQuality().trkKink < 20)
        ###check Muon ID
	return muon.mediumId

    def muonIso(self. muon):
        ###check Muon iso 
        return muon.pfRelIso03_all < 0.15
    
    def checkTwoleptons(self, leptons_mu_keys, leptons_el_keys):
        passdilepton = False
	if self.triggertype == "DoubleMuon":
	    if len(leptons_mu_keys) >= 2:
	        passdilepton = True
	elif self.triggertype == "DoubleEG":
	    if len(leptons_el_keys) >= 2:
	        passdilepton = True
	elif self.triggertype == "MuonEG":
	    if len(leptons_mu_keys) >= 1 and len(leptons_el_keys) >= 1:
	        passdilepton = True
        return passdilepton
	

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
	   
    def HLTPath(self, hlt, isMC):
	""" check which HLT is fired """
	L1t_hlt = {"DoubleMuon": ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ"], 
	    	   "MuonEG": ["HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL","HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 
		   		"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL","HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ"],
		   "DoubleEG": ["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL","HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"]
			}
        if isMC:
	    for l1 in L1t_hlt:
	    	for path in L1t_hlt[l1]:
		    firepath = getattr(hlt,  path, False)
    		    if firepath:
		    	self.triggertype = l1
		    	return True, path
	    
	else:
	    for path in L1t_hlt[self.triggertype]:
		firepath = getattr(hlt,  path, False)
		if firepath:
		    return True, path
	return False,"" 

    def trigger_reco_matching(self, muons, mu_effSF, trigobjs, leptons):
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
			    print "trigger obj l1 pt ", tobj.l1pt," l2pt ", tobj.l2pt," pt ",tobj.pt, " eta ",tobj.eta, " bits ", tobj.filterBits, " offlep pt ",mu.pt," eta ",mu.eta
			triggermatch = True
			break;
		if triggermatch:
		    leptons[mu] = mu_effSF[imu]

    def fillMCinfo(self, event):
	""" fill all gen information here """
	if not self.isMC:
		return 

	genmet = Object(event)

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
	#hlts = Collection(event, "HLT")

	
	event_reco_weight = 1.0 ## for pu_weight*btag_SF*lepSF
	sample_weight = 1.0
	###cutstep1: initial events, no selection, 
	fired, path =  self.HLTPath(event, self.isMC)
	##cutstep2: HLT is fired or not
	cutflow_bin = 1.0
	if not fired:
	    self.h_cutflow.Fill( cutflow_bin, event_reco_weight * sample_weight)
	    return False
	else:
	    cutflow_bin += 1.0
	
	self.ievent += 1
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
	        
	if self.verbose > 2:
		print "iEvent ",self.ievent," HLT path ",path," leptype ",Lepstype

        ### PV	
	PV = Object(event, "PV")
	pu_weight = event.puWeight
        #self.out.fillBranch("pu", PV.npvs)
        self.out.fillBranch("event_pu_weight", pu_weight)

	### MET
        met = Object(event, "MET")
        #metPt,metPhi = self.met(met,self.isMC)
	metPt = met.pt; metPhi = met.phi
        #self.out.fillBranch("met_pt",metPt)
        #self.out.fillBranch("met_phi",metPhi) 

        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
	mu_effSF = event.Muon_effSF
	el_effSF = event.Electron_effSF
        jets = Collection(event, "Jet")
	jet_btagSF = event.Jet_btagSF
	njets = len(jets)
	trigobjs = Collection(event, "TrigObj")
	

	###match trigger obj to HLT?
	trigobjs_lep = []
  	for tobj in trigobjs:
	    if (tobj.id == 11 and tobj.l1pt >= 12) or (tobj.id == 13 and tobj.l1pt >= 8):
		trigobjs_lep.append(tobj)
		if self.verbose > 3:
	            print "triggerobj id  ",tobj.id," l1pt ",tobj.l1pt," l2pt ",tobj.l2pt," pt ",tobj.pt," eta ",tobj.eta," phi ",tobj.phi," bits ",tobj.filterBits	
	leptons_mu = {}
	leptons_el = {}
	
	if self.triggertype == "DoubleMuon" or self.triggertype == "MuonEG":
            self.trigger_reco_matching(muons, mu_effSF, trigobjs_lep, leptons_mu)
	if self.triggertype == "DoubleEG" or self.triggertype == "MuonEG":
            self.trigger_reco_matching(electrons, el_effSF, trigobjs_lep, leptons_el)
	

	###cutstep3: L1_HLT_RECO matching
	if (self.triggertype == "DoubleMuon" and len(leptons_mu) >= 2) or (self.triggertype == "DoubleEG" and len(leptons_el) >= 2) or (self.triggertype == "MuonEG" and len(leptons_mu)>=1 and len(leptons_el) >= 1):
	    cutflow_bin += 1.0
	else:
	    self.h_cutflow.Fill( cutflow_bin, event_reco_weight * sample_weight)
	    return False
		
	


	leptons_el_keys = leptons_el.keys()
        leptons_el_keys.sort(key=lambda x:x.pt,reverse=True)	
	for el in leptons_el_keys:
	    if not(self.electronImpactParameterCut(el)):
		leptons_el_keys.remove(el)
        leptons_mu_keys = leptons_mu.keys()
        leptons_mu_keys.sort(key=lambda x:x.pt,reverse=True)	
	passdileptonPtEta = False
	###cutstep4: dilepton pt, and eta, 
	if self.triggertype == "DoubleMuon":
	    for mu in leptons_mu_keys:
	        if mu.pt < self.subleadingMuonPt["DoubleMuon"] :
		    leptons_mu_keys.remove( mu )
            if len(leptons_mu_keys) >= 2 and leptons_mu_keys[0].pt >= self.leadingMuonPt["DoubleMuon"]:
	        cutflow_bin += 1.0
	    else:
	        self.h_cutflow.Fill( cutflow_bin, event_reco_weight * sample_weight)
	        return False
	elif self.triggertype == "DoubleEG":
	    for el in leptons_el_keys:
	        if el.pt < self.subleadingEGPt["DoubleEG"] :
		    leptons_el_keys.remove( el )
            if len(leptons_el_keys) >= 2 and leptons_el_keys[0].pt >= self.leadingEGPt["DoubleEG"]:
	        cutflow_bin += 1.0
	    else:
	        self.h_cutflow.Fill( cutflow_bin, event_reco_weight * sample_weight)
	        return False
	elif self.triggertype == "MuonEG":
            if len(leptons_mu_keys) >= 1 and len(leptons_el_keys) >= 1:
	        if (leptons_mu_keys[0].pt >= self.leadingMuonPt["MuonEG"] and leptons_el_keys[0].pt >= self.subleadingEGPt["MuonEG"]):
		     cutflow_bin += 1.0
		elif (leptons_mu_keys[0].pt >= self.subleadingMuonPt["MuonEG"] and leptons_el_keys[0].pt >= self.leadingEGPt["MuonEG"]):
		     cutflow_bin += 1.0
		else:
	            self.h_cutflow.Fill( cutflow_bin, event_reco_weight * sample_weight)
	            return False
	    else:
	        self.h_cutflow.Fill( cutflow_bin, event_reco_weight * sample_weight)
	        return False
	    
		
	###cutstep5: dilepton, ISO
	passdileptonIso = False
	for el in leptons_el_keys:
	    if not(self.electronIso(el)):
		leptons_el_keys.remove(el)
	for mu in leptons_mu_keys:
	    if not(self.muonIso(mu)):
		leptons_mu_keys.remove(mu)
	
	passdileptonIso =  self.checkTwoleptons(leptons_mu_keys, leptons_el_keys)

	if passdileptonIso:
	     cutflow_bin += 1.0
	else:
	     self.h_cutflow.Fill( cutflow_bin, event_reco_weight * sample_weight)
	     return False
	
	###cutstep6: dilepton, ID
	for el in leptons_el_keys:
	    if not(self.electronID(el)):
		leptons_el_keys.remove(el)
	for mu in leptons_mu_keys:
	    if not(self.muonID(mu)):
		leptons_mu_keys.remove(mu)
	
	passdileptonID =  self.checkTwoleptons(leptons_mu_keys, leptons_el_keys)
	if passdileptonID:
	     cutflow_bin += 1.0
	else:
	     self.h_cutflow.Fill( cutflow_bin, event_reco_weight * sample_weight)
	     return False


	###cutstep6: dilepton, ID
	for el in leptons_el_keys:
	    if not(self.electronHLTSafeID(el)):
		leptons_el_keys.remove(el)

	passdileptonHLTSafeID =  self.checkTwoleptons(leptons_mu_keys, leptons_el_keys)

	if passdileptonHLTSafeID:
	     cutflow_bin += 1.0
	else:
	     self.h_cutflow.Fill( cutflow_bin, event_reco_weight * sample_weight)
	     return False

        ## selection final two leptons	

	##Lepton Selection
        #wMuons = [x for x in muons if x.pt > 25 and x.tightId >= 1 and x.pfRelIso04_all < 0.15 and x.dxy < 0.05 and x.dz < 0.2]
        #wElectrons = [x for x in electrons if x.mvaSpring16GP_WP80 and x.pt > 25 and x.pfRelIso03_all < 0.12]      

        #wMuons.sort(key=lambda x:x.pt,reverse=True)
        #wElectrons.sort(key=lambda x:x.pt,reverse=True)
      
	self.out.fillBranch("lepstype", Lepstype)
        ## add branches for some basic kinematics
        #ll = ROOT.TLorentzVector()



        #vLeptons = [] # decay products of H->WW

        #for vLepton in vLeptons:
        #    vLepton_4vec = ROOT.TLorentzVector()
        #    vLepton_4vec.SetPtEtaPhiM(vLepton.pt,vLepton.eta,vLepton.phi,vLepton.mass)
        #    V = V + vLepton_4vec
        #if Vtype >=2 and Vtype<=4:
        #    met_4vec = ROOT.TLorentzVector()
        #    met_4vec.SetPtEtaPhiM(met.pt,0.,met.phi,0.) # only use met vector to derive transverse quantities
        #    V = V + met_4vec
        #self.out.fillBranch("V_pt",V.Pt())
        #self.out.fillBranch("V_eta",V.Eta())
        #self.out.fillBranch("V_phi",V.Phi())
        #self.out.fillBranch("V_mass",V.M())
        #self.out.fillBranch("V_mt",V.Mt())

        ### filter jets that overlap with any of the selected leptons
        #allLeptons = zElectrons[:]
        #allLeptons.extend(zMuons)
        #allLeptons.extend(wElectrons)
        #allLeptons.extend(wMuons)
        #jetFilterFlags = [True]*len(jets)
        ##for jet in jets:
        ##    jet.jetFilter = True
        #for lepton in allLeptons:
        #    jetInd = lepton.jetIdx
        #    if jetInd >= 0:
        #        jetFilterFlags[jetInd] = False
        #        #jets[jetInd].jetFilter = False
        #self.out.fillBranch("Jet_lepFilter",jetFilterFlags)

        ### alias JER-smeared MC jet pT and data jet pT to the same
        ### branch name
        #jetPts = [-99.]*len(jets)
        #for i in xrange(len(jets)):
        #    jetPts[i] = self.pt(jets[i],self.isMC)

        #self.out.fillBranch("Jet_Pt",jetPts)

        ### Add explicit indices for selected H(bb) candidate jets
        #jetsForHiggs = [x for x in jets if x.lepFilter and x.puId>0 and x.jetId>0 and x.Pt>20 and abs(x.eta)<2.5]
        #if (len(jetsForHiggs) < 2): return False
        #hJets = sorted(jetsForHiggs, key = lambda jet : jet.btagCMVA, reverse=True)[0:2]
        #hJidx = [jets.index(x) for x in hJets]
        #self.out.fillBranch("hJidx",hJidx)

        ### Save a few basic reco. H kinematics
        #hj1 = ROOT.TLorentzVector()
        #hj2 = ROOT.TLorentzVector()
        #hj1.SetPtEtaPhiM(jets[hJidx[0]].Pt,jets[hJidx[0]].eta,jets[hJidx[0]].phi,jets[hJidx[0]].mass)
        #hj2.SetPtEtaPhiM(jets[hJidx[1]].Pt,jets[hJidx[1]].eta,jets[hJidx[1]].phi,jets[hJidx[1]].mass)
        #hbb = hj1 + hj2
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

hhbbWW = lambda : HHbbWWProducer(True, "") 
hhbbWW_data = lambda x : HHbbWWProducer(False, x) 
