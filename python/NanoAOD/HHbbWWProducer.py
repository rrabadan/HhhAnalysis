import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection,Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import * #deltaR, matching etc..

class HHbbWWProducer(Module):
    ## data or MC, which L1 trigger, HLT?
    def __init__(self, isMC, L1trigger):
        self.isMC = isMC ## two mode: data or MC
	self.triggertype  = L1trigger##"DoubleMuon, DoubleEG, MuonEG"
        pass
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("sample_weight",  "F");
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
	self.out.branch("event_weight",  "F")
	self.out.branch("event_pu_weight",  "F")
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

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
	#hlts = Collection(event, "HLT")
	fired, path =  self.HLTPath(event, self.isMC)
	if not fired:
	    return False
	
        Lepstype = -1##MuMu:1, MuEl:2, ElMu:3, ElEl:4
	if self.triggertype == "DoubleMuon":
        	Lepstype = 1	
	elif self.triggertype == "DoubleEG":
		Lepstype = 4
	elif self.triggertype == "MuonEG" and "Mu23" in path:
		Lepstype = 2
	elif self.triggertype == "MuonEG" and "Ele23" in path:
		Lepstype = 3
	else:
		return False

	
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        jets = list(Collection(event, "Jet"))
        met = Object(event, "MET")

        #metPt,metPhi = self.met(met,self.isMC)
	metPt = met.pt; metPhi = met.phi
	print "Lepstype ",Lepstype," metpt ",metPt, " metphi ",metPhi
        self.out.fillBranch("met_pt",metPt)
        self.out.fillBranch("met_phi",metPhi) 
      
	self.out.fillBranch("lepstype", Lepstype)
        ## add branches for some basic kinematics
        #ll = ROOT.TLorentzVector()

        #wElectrons = [x for x in electrons if x.mvaSpring16GP_WP80 and x.pt > 25 and x.pfRelIso03_all < 0.12]      
        #wMuons = [x for x in muons if x.pt > 25 and x.tightId >= 1 and x.pfRelIso04_all < 0.15 and x.dxy < 0.05 and x.dz < 0.2]

        #wMuons.sort(key=lambda x:x.pt,reverse=True)
        #wElectrons.sort(key=lambda x:x.pt,reverse=True)

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
