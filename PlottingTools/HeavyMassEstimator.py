
import ROOT
import numpy as np
from math import *


class HeavyMassEstimator(object):
    """ Heavy  Mass Estimator 
    Attributes: 
    	two muons lorentz vectors
        two bjet lorentz vectors 
        missing ET
        root fils contains all PDFs
    hmetree should contain full information from  HME
    """
    iterations = 100
    onshellnuptpdf = ROOT.TH1F()
    onshellWmasspdf = ROOT.TH1F()
    onshellnuptpdf_flag = False
    onshellWmasspdf_flag = False
    hmetree = ROOT.TTree("hmetree","HME Tree") 
    hme_h2Mass = ROOT.TH1F("hme_h2Mass","h2 mass from HME",1000, 200.0,1200.0)
    hme_h2MassWeight1 = ROOT.TH1F("hme_h2MassWeight1","h2 mass from HME",1000, 200.0,1200.0)
    hme_h2MassWeight4 = ROOT.TH1F("hme_h2MassWeight4","h2 mass from HME",1000, 200.0,1200.0)
    eta_gen  = np.zeros(1, dtype=float);   phi_gen  = np.zeros(1, dtype=float)
    wmass_gen =  np.zeros(1, dtype=float); hmass_gen = np.zeros(1, dtype=float)
    metpx_gen = np.zeros(1, dtype=float);  metpy_gen = np.zeros(1, dtype=float);
    nsolutions = np.zeros(1, dtype=int)
    weight = np.zeros(1, dtype=float)
    weight1 = np.zeros(1, dtype=float)
    weight2 = np.zeros(1, dtype=float)
    weight3 = np.zeros(1, dtype=float)
    weight4 = np.zeros(1, dtype=float)
    l_onshellW_eta = np.zeros(1, dtype=float); l_onshellW_phi = np.zeros(1, dtype=float); l_onshellW_pt = np.zeros(1, dtype=float); l_onshellW_energy = np.zeros(1, dtype=float)
    l_offshellW_eta = np.zeros(1, dtype=float); l_offshellW_phi = np.zeros(1, dtype=float); l_offshellW_pt = np.zeros(1, dtype=float); l_offshellW_energy = np.zeros(1, dtype=float)
    nu_onshellW_eta = np.zeros(1, dtype=float); nu_onshellW_phi = np.zeros(1, dtype=float); nu_onshellW_pt = np.zeros(1, dtype=float); nu_onshellW_energy = np.zeros(1, dtype=float)
    nu_offshellW_eta = np.zeros(1, dtype=float); nu_offshellW_phi = np.zeros(1, dtype=float); nu_offshellW_pt = np.zeros(1, dtype=float); nu_offshellW_energy = np.zeros(1, dtype=float)
    onshellW_eta = np.zeros(1, dtype=float); onshellW_phi = np.zeros(1, dtype=float); onshellW_pt = np.zeros(1, dtype=float); onshellW_energy = np.zeros(1, dtype=float); onshellW_mass = np.zeros(1, dtype=float)
    offshellW_eta = np.zeros(1, dtype=float); offshellW_phi = np.zeros(1, dtype=float); offshellW_pt = np.zeros(1, dtype=float); offshellW_energy = np.zeros(1, dtype=float); offshellW_mass = np.zeros(1, dtype=float)
    met_pt = np.zeros(1, dtype=float); met_phi = np.zeros(1, dtype=float); met_px = np.zeros(1, dtype=float); met_py = np.zeros(1, dtype=float) 
    b1jet_eta = np.zeros(1, dtype=float); b1jet_phi = np.zeros(1, dtype=float); b1jet_pt = np.zeros(1, dtype=float); b1jet_energy = np.zeros(1, dtype=float)
    b2jet_eta = np.zeros(1, dtype=float); b2jet_phi = np.zeros(1, dtype=float); b2jet_pt = np.zeros(1, dtype=float); b2jet_energy = np.zeros(1, dtype=float)
    htoBB_eta = np.zeros(1, dtype=float); htoBB_phi = np.zeros(1, dtype=float); htoBB_pt = np.zeros(1, dtype=float); htoBB_energy = np.zeros(1, dtype=float)
    htoWW_eta = np.zeros(1, dtype=float); htoWW_phi = np.zeros(1, dtype=float); htoWW_pt = np.zeros(1, dtype=float); htoWW_energy = np.zeros(1, dtype=float)
    htoBB_mass = np.zeros(1, dtype=float); htoWW_mass = np.zeros(1, dtype=float);
    h2tohh_eta = np.zeros(1, dtype=float); h2tohh_phi = np.zeros(1, dtype=float); h2tohh_pt = np.zeros(1, dtype=float); h2tohh_energy = np.zeros(1, dtype=float); h2tohh_mass = np.zeros(1, dtype=float)

    hmetree.Branch('nsolutions', nsolutions, 'nsolutions/I')
    hmetree.Branch('l_onshellW_eta', l_onshellW_eta, 'l_onshellW_eta/D')  
    hmetree.Branch('l_onshellW_phi', l_onshellW_phi, 'l_onshellW_phi/D')  
    hmetree.Branch('l_onshellW_pt', l_onshellW_pt, 'l_onshellW_pt/D')  
    hmetree.Branch('l_onshellW_energy', l_onshellW_energy, 'l_onshellW_energy/D')  
    hmetree.Branch('l_offshellW_eta', l_offshellW_eta, 'l_offshellW_eta/D')  
    hmetree.Branch('l_offshellW_phi', l_offshellW_phi, 'l_offshellW_phi/D')  
    hmetree.Branch('l_offshellW_pt', l_offshellW_pt, 'l_offshellW_pt/D')  
    hmetree.Branch('l_offshellW_energy', l_offshellW_energy, 'l_offshellW_energy/D')  
    hmetree.Branch('nu_onshellW_eta', nu_onshellW_eta, 'nu_onshellW_eta/D')  
    hmetree.Branch('nu_onshellW_phi', nu_onshellW_phi, 'nu_onshellW_phi/D')  
    hmetree.Branch('nu_onshellW_pt', nu_onshellW_pt, 'nu_onshellW_pt/D')  
    hmetree.Branch('nu_onshellW_energy', nu_onshellW_energy, 'nu_onshellW_energy/D')  
    hmetree.Branch('nu_offshellW_eta', nu_offshellW_eta, 'nu_offshellW_eta/D')  
    hmetree.Branch('nu_offshellW_phi', nu_offshellW_phi, 'nu_offshellW_phi/D')  
    hmetree.Branch('nu_offshellW_pt', nu_offshellW_pt, 'nu_offshellW_pt/D')  
    hmetree.Branch('nu_offshellW_energy', nu_offshellW_energy, 'nu_offshellW_energy/D')  
    hmetree.Branch('onshellW_eta', onshellW_eta, 'onshellW_eta/D')  
    hmetree.Branch('onshellW_phi', onshellW_phi, 'onshellW_phi/D')  
    hmetree.Branch('onshellW_pt', onshellW_pt, 'onshellW_pt/D')  
    hmetree.Branch('onshellW_energy', onshellW_energy, 'onshellW_energy/D')  
    hmetree.Branch('onshellW_mass', onshellW_mass, 'onshellW_mass/D')  
    hmetree.Branch('offshellW_eta', offshellW_eta, 'offshellW_eta/D')  
    hmetree.Branch('offshellW_phi', offshellW_phi, 'offshellW_phi/D')  
    hmetree.Branch('offshellW_pt', offshellW_pt, 'offshellW_pt/D')  
    hmetree.Branch('offshellW_energy', offshellW_energy, 'offshellW_energy/D')  
    hmetree.Branch('offshellW_mass', offshellW_mass, 'offshellW_mass/D')  
    hmetree.Branch('met_pt', met_pt, 'met_pt')
    hmetree.Branch('met_phi', met_phi, 'met_phi')
    hmetree.Branch('met_px', met_px, 'met_px')
    hmetree.Branch('met_py', met_py, 'met_py')
    hmetree.Branch('b1jet_eta', b1jet_eta, 'b1jet_eta/D')
    hmetree.Branch('b1jet_phi', b1jet_phi, 'b1jet_phi/D')
    hmetree.Branch('b1jet_pt', b1jet_pt, 'b1jet_pt/D')
    hmetree.Branch('b1jet_energy', b1jet_energy, 'b1jet_energy/D')
    hmetree.Branch('b2jet_eta', b2jet_eta, 'b2jet_eta/D')
    hmetree.Branch('b2jet_phi', b2jet_phi, 'b2jet_phi/D')
    hmetree.Branch('b2jet_pt', b2jet_pt, 'b2jet_pt/D')
    hmetree.Branch('b2jet_energy', b2jet_energy, 'b2jet_energy/D')
    hmetree.Branch('htoBB_eta', htoBB_eta, 'htoBB_eta/D')
    hmetree.Branch('htoBB_phi', htoBB_phi, 'htoBB_phi/D')
    hmetree.Branch('htoBB_pt', htoBB_pt, 'htoBB_pt/D')
    hmetree.Branch('htoBB_energy', htoBB_energy, 'htoBB_energy/D')
    hmetree.Branch('htoBB_mass', htoBB_mass, 'htoBB_mass/D')
    hmetree.Branch('htoWW_eta', htoWW_eta, 'htoWW_eta/D')
    hmetree.Branch('htoWW_phi', htoWW_phi, 'htoWW_phi/D')
    hmetree.Branch('htoWW_pt', htoWW_pt, 'htoWW_pt/D')
    hmetree.Branch('htoWW_energy', htoWW_energy, 'htoWW_energy/D')
    hmetree.Branch('htoWW_mass', htoWW_mass, 'htoWW_mass/D')
    hmetree.Branch('h2tohh_eta', h2tohh_eta, 'h2tohh_eta/D')
    hmetree.Branch('h2tohh_phi', h2tohh_phi, 'h2tohh_phi/D')
    hmetree.Branch('h2tohh_pt', h2tohh_pt, 'h2tohh_pt/D')
    hmetree.Branch('h2tohh_energy', h2tohh_energy, 'h2tohh_energy/D')
    hmetree.Branch('h2tohh_mass', h2tohh_mass, 'h2tohh_mass/D')


    lepton1_p4  = ROOT.TLorentzVector()
    lepton2_p4  = ROOT.TLorentzVector()
    b1jet_p4  = ROOT.TLorentzVector()
    b2jet_p4  = ROOT.TLorentzVector()
    met = ROOT.TVector2()
    wmasshist = ROOT.TH1F()
    onshellnupthist = ROOT.TH1F()

    lepton1_onshellW_p4 = ROOT.TLorentzVector()
    lepton1_offshellW_p4 = ROOT.TLorentzVector()
    nu_onshellW_p4 = ROOT.TLorentzVector()
    nu_offshellW_p4 = ROOT.TLorentzVector()
    onshellW_p4 = ROOT.TLorentzVector()
    onshellW_p4 = ROOT.TLorentzVector()
    htoWW_p4 =  ROOT.TLorentzVector()
    htoBB_p4 =  ROOT.TLorentzVector()
    h2tohh_p4 = ROOT.TLorentzVector()



    def __init__(self):
	print "  create a HeavyMassEstimator object "
	#self.lepton1_p4  = ROOT.TLorentzVector()
    def setKinematic(self, lepton1_p4, lepton2_p4, jet1_p4, jet2_p4, met):
	self.lepton1_p4 = lepton1_p4
	self.lepton2_p4 = lepton2_p4
	self.b1jet_p4 = jet1_p4
	self.b2jet_p4 = jet2_p4
	self.met = met
	#self.RefPDFFileName = RefPDFFileName
	#self.RefPDFFile = ROOT.TFile(RefPDFFileName,"READ")

    def setonshellWmasspdf(self, hist):
	self.onshellWmasspdf = hist
	self.onshellWmasspdf_flag = True
    
    def setonshellnuptpdf(self, hist):
	self.onshellnuptpdf = hist
	self.onshellnuptpdf_flag = True

    def showKinematic(self):
	print "lepton1 ",self.lepton1_p4.Print()
	print "lepton2 ",self.lepton2_p4.Print()
	print "b1jet ",self.b1jet_p4.Print()
	print "b2jet ",self.b2jet_p4.Print()
	print "Met ",self.met.Print()

    def initHMETree(self):
	""" intialize the HME tree """
	self.l_onshellW_eta[0] = -9.0; self.l_onshellW_phi[0] = -9.0; self.l_onshellW_pt[0] = -1.0; self.l_onshellW_energy[0] = -1.0
	self.l_offshellW_eta[0] = -9.0; self.l_offshellW_phi[0] = -9.0; self.l_offshellW_pt[0] = -1.0; self.l_offshellW_energy[0] = -1.0
	self.nu_onshellW_eta[0] = -9.0; self.nu_onshellW_phi[0] = -9.0; self.nu_onshellW_pt[0] = -1.0; self.nu_onshellW_energy[0] = -1.0
	self.nu_offshellW_eta[0] = -9.0; self.nu_offshellW_phi[0] = -9.0; self.nu_offshellW_pt[0] = -1.0; self.nu_offshellW_energy[0] = -1.0
	self.onshellW_eta[0] = -9.0; self.onshellW_phi[0] = -9.0; self.onshellW_pt[0] = -1.0; self.onshellW_energy[0] = -1.0; self.onshellW_mass[0] = -1.0
	self.offshellW_eta[0] = -9.0; self.offshellW_phi[0] = -9.0; self.offshellW_pt[0] = -1.0; self.offshellW_energy[0] = -1.0;self.offshellW_mass[0] = -1.0
	self.htoWW_eta[0] = -9.0; self.htoWW_phi[0] = -9.0; self.htoWW_pt[0] = -1.0; self.htoWW_energy[0] = -1.0; self.htoWW_mass[0] = -1.0
	self.b1jet_eta[0] = -9.0; self.b1jet_phi[0] = -9.0; self.b1jet_pt[0] = -1.0; self.b1jet_energy[0] = -1.0
	self.b2jet_eta[0] = -9.0; self.b2jet_phi[0] = -9.0; self.b2jet_pt[0] = -1.0; self.b2jet_energy[0] = -1.0
	self.htoBB_eta[0] = -9.0; self.htoBB_phi[0] = -9.0; self.htoBB_pt[0] = -1.0; self.htoBB_energy[0] = -1.0; self.htoBB_mass[0] = -1.0
	self.h2tohh_eta[0] = -9.0; self.h2tohh_phi[0] = -9.0; self.h2tohh_pt[0] = -1.0; self.h2tohh_energy[0] = -1.0; self.h2tohh_mass[0] = -1.0
	self.met_pt[0] = -1.0; self.met_px[0] = -99999.0; self.met_py[0] = -99999.0; self.met_phi[0] = -99999.0
	self.weight[0] = 1.0; self.weight1[0] = 1.0;  self.weight2[0] = 1.0; self.weight3[0] = 1.0; self.weight4[0] = 1.0
	self.nsolutions[0] = 0
   
    def getOnshellWMass(self, x0, step, random):
	
	xmin = 50.0; xmax = 90.0
        while (x0 > xmax or x0 < xmin):
	    if x0 > xmax:
	    	x0 = x0 - xmax + xmin
	    if x0 < xmin:
	        x0 = xmax - (xmin - x0)
	x1 = x0 + step
	while (x1 > xmax or x1 < xmin):
	    if x1 > xmax:
	        x1 = x1 - xmax + xmin
	    if x1 < xmin:
	        x1 = xmax - (xmin - x1)
	binx00 = self.onshellWmasspdf.FindBin(x0)
	binx01 = binx00 + 1
	binx10 = self.onshellWmasspdf.FindBin(x1)
	binx11 = binx10 + 1
        if self.onshellWmasspdf.GetBinCenter(binx00) > x0:
     	    binx00 = binx00 - 1
     	    binx01 = binx01 - 1
        if self.onshellWmasspdf.GetBinCenter(binx10) > x1:
     	    binx10 = binx10 - 1
     	    binx11 = binx11 - 1
        y00 = self.onshellWmasspdf.GetBinContent( binx00 )
        x00 = self.onshellWmasspdf.GetBinCenter( binx00 )
        y01 = self.onshellWmasspdf.GetBinContent( binx01 )
        x01 = self.onshellWmasspdf.GetBinCenter( binx01 )
        y10 = self.onshellWmasspdf.GetBinContent( binx10 )
        x10 = self.onshellWmasspdf.GetBinCenter( binx10 )
        y11 = self.onshellWmasspdf.GetBinContent( binx11 )
        x11 = self.onshellWmasspdf.GetBinCenter( binx11 )
        w0  = (x0 - x00)*(y00 - y01)/(x00 - x01) + y00;
        w1  = (x1 - x10)*(y10 - y11)/(x10 - x11) + y10;
	#w1/w0: transition probability 
	if (w1/w0 >= random):
	    return x1
	elif (w1/w0 < random):
	    return x0
        else:	
            print "error in getOnshellWMass "
            return 80.3
            
    
    def assignMuP4(self, case):
	"""lepton+nu permutation
	    in simulation:
	    case 0: assgin lepton from onshell W to l_onshellW_p4
	    case 1: assgin lepton from offshell W to l_onshellW_p4
	   not in simuation: what ever comes first is assign to l_onshellW_p4
	"""
	if case == 0:
	    self.lepton_onshellW_p4 = self.lepton1_p4
	    self.lepton_offshellW_p4 = self.lepton2_p4
	elif case == 1:
	    self.lepton_offshellW_p4 = self.lepton1_p4
	    self.lepton_onshellW_p4 = self.lepton2_p4
	
    def nuPtFromOnshellW(self, nu_eta, nu_phi, lepton_p4, wMass):
	deta = nu_eta-lepton_p4.Eta()
	dphi = nu_phi-lepton_p4.Phi()
	nuPt = wMass*wMass/(2*lepton_p4.Pt()*(cosh(deta)-cosh(dphi)))
        return nuPt

    def nuP4FromOffshellW(self, met, lepton1_p4, lepton2_p4, nu1_p4, nu2_p4, case, hMass):
	tmp_p4 = lepton1_p4 + lepton2_p4 + nu1_p4
        tmp_nu_px = met.Px() - nu1_p4.Px()
        tmp_nu_py = met.Py() - nu1_p4.Py()
        nu_pxpy =  ROOT.TVector2(tmp_nu_px, tmp_nu_py)
        tmp_nu_pt = nu_pxpy.Mod()
        tmp_p4_v2 = ROOT.TLorentzVector(sqrt(pow(tmp_p4.Pt(), 2) + pow(tmp_p4.M(), 2)), 0, tmp_p4.Pz(), tmp_p4.Energy())
        chdeta = (pow(hMass, 2) + 2*(nu_pxpy.Px()*tmp_p4.Px() + nu_pxpy.Py()*tmp_p4.Py()) - pow(tmp_p4.M(), 2))/(2.0*tmp_p4_v2.Pt()*tmp_nu_pt)
        if chdeta < 1.0:
        #no solution if chdeta<1.0
    	    #print "no solution since chdeta<1.0, chdeta ",chdeta
	    nu2_p4.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0)
            return False 
	tmp_nu_phi = nu_pxpy.Phi_mpi_pi(nu_pxpy.Phi())
        deta = acosh( chdeta )
        tmp_nu_eta = 0.0
        if case == 1:
	    tmp_nu_eta = tmp_p4_v2.Eta() - deta
	else :
	    tmp_nu_eta = tmp_p4_v2.Eta() + deta
        if (abs(tmp_nu_eta) > 6.0):
        #very unlikely solution
	    print "tmp_nu_eta ",tmp_nu_eta, " very unlikely solution, pass"
	    nu2_p4.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0)
            return False 
	nu2_p4.SetPtEtaPhiM(tmp_nu_pt, tmp_nu_eta, tmp_nu_phi, 0.0)	
        htoWW_tmp = tmp_p4 + nu2_p4
	if abs(htoWW_tmp.M() - hMass) > 1.0:
	    print "Warning!!! gen hmass ", hMass, " HME htoWW mass ", htoWW_tmp.M()
	return True

    def runHME(self):
	if not(self.onshellWmasspdf_flag):
	    print "no onshellWmasspdf, error!!! "
	    return  False
	self.htoBB_p4 = self.b1jet_p4 + self.b2jet_p4
	self.eta_gen[0] = 0.0 ; self.phi_gen[0] = 0.0
	#initial wmass_gen
	self.wmass_gen[0] = 70.0
	self.met_px[0] = self.met.Px(); self.met_py[0] = self.met.Py()
	it = 0
	genRandom = ROOT.TRandom3(0)
        #PUSample: 25.2, PU0: 14.8
        met_sigma = 25.2
        #genRandom.SetSeed()
	while (it < self.iterations ):
	    it += 1
	    self.initHMETree()
	    self.eta_gen[0] = genRandom.Uniform(-6, 6)
	    self.phi_gen[0] = genRandom.Uniform(-3.1415926, 3.1415926)
	    self.hmass_gen[0] = genRandom.Gaus(125.03, 0.004)
	    rand01 = genRandom.Uniform(0., 1.0)
	    step = genRandom.Uniform(-10.0, 10.0)
	    self.wmass_gen[0] = self.getOnshellWMass(self.wmass_gen[0], step, rand01)
	    #print "it ",it," self.eta_gen[0] ",self.eta_gen[0]," wmass_gen ",self.wmass_gen[0]
	    #update met 
	    met_dpx = genRandom.Gaus(0.0, met_sigma)
	    met_dpy = genRandom.Gaus(0.0, met_sigma)
            met_corr = self.met + ROOT.TVector2(met_dpx, met_dpy)
	    self.metpx_gen[0]= met_corr.Px()
	    self.metpy_gen[0] = met_corr.Py()
            #print "met_px ",self.met.Px()," met_py ",self.met.Py()," after correction px ",met_corr.Px()," py ",met_corr.Py()
	    self.nsolutions[0] = 0
	    isolution = 0
	    solutions = [False, False, False, False]
	    #1. permutation 
	    #2. check nu_onshell_W pt
	    #3. solve the kinematics
	    #4. mark solution is ture if it is solved
	    #5. dump information into tree
	    while isolution < len(solutions):
		case = isolution/2
		#print "case ",case," isolution ",isolution
		self.assignMuP4(case)
	        nu_onshellW_pt_tmp = self.nuPtFromOnshellW(self.eta_gen[0], self.phi_gen[0], self.lepton_onshellW_p4, self.wmass_gen[0])
	        nu_onshellW_p4_tmp = ROOT.TLorentzVector()
	        nu_onshellW_p4_tmp.SetPtEtaPhiM(nu_onshellW_pt_tmp, self.eta_gen[0], self.phi_gen[0], 0)
	        nu_offshellW_p4_tmp = ROOT.TLorentzVector()
	        solutions[isolution] = self.nuP4FromOffshellW(met_corr, self.lepton_onshellW_p4, self.lepton_offshellW_p4, nu_onshellW_p4_tmp, nu_offshellW_p4_tmp, isolution%2, self.hmass_gen[0])
	 	if solutions[isolution]:
		    self.nsolutions[0] += 1
	 	isolution += 1 
	    isolution = 0
	    if self.nsolutions[0] == 0:#no soution in this iterations
	 	continue
	    self.weight[0] = 1.0/float(self.nsolutions[0])
	    while isolution < len(solutions):
		if not(solutions[isolution]):	
		    isolution += 1
		    continue
		case = isolution/2
		self.assignMuP4(case)
		self.nu_onshellW_pt[0] = self.nuPtFromOnshellW(self.eta_gen[0], self.phi_gen[0], self.lepton_onshellW_p4, self.wmass_gen[0])
		self.nu_onshellW_p4.SetPtEtaPhiM(self.nu_onshellW_pt[0], self.eta_gen[0], self.phi_gen[0], 0)
		self.nuP4FromOffshellW(met_corr, self.lepton_onshellW_p4, self.lepton_offshellW_p4, self.nu_onshellW_p4, self.nu_offshellW_p4, isolution%2, self.hmass_gen[0])
	        self.onshellW_p4 = self.lepton_onshellW_p4+self.nu_onshellW_p4
	        self.offshellW_p4 = self.lepton_offshellW_p4+self.nu_offshellW_p4
		self.htoWW_p4 = self.onshellW_p4 + self.offshellW_p4
		self.h2tohh_p4 = self.htoWW_p4 + self.htoBB_p4

		if (fabs(self.htoWW_p4.M() - self.hmass_gen[0])>1.0):
		    print "Error!! hmass_gen ", self.hmass_gen[0], " higgs mass from HME htoWW_p4 ", self.htoWW_p4.M()
		#print "get this h2tohh_mass ",self.h2tohh_p4.M()," iter ",it
		self.l_onshellW_eta[0] = self.lepton_onshellW_p4.Eta()
		self.l_onshellW_phi[0] = self.lepton_onshellW_p4.Phi()
		self.l_onshellW_pt[0] = self.lepton_onshellW_p4.Pt()
		self.l_onshellW_energy[0] = self.lepton_onshellW_p4.Energy()
		self.l_offshellW_eta[0] = self.lepton_offshellW_p4.Eta()
		self.l_offshellW_phi[0] = self.lepton_offshellW_p4.Phi()
		self.l_offshellW_pt[0] = self.lepton_offshellW_p4.Pt()
		self.l_offshellW_energy[0] = self.lepton_offshellW_p4.Energy()

		self.nu_onshellW_eta[0] = self.nu_onshellW_p4.Eta()
		self.nu_onshellW_phi[0] = self.nu_onshellW_p4.Phi()
		self.nu_onshellW_pt[0] = self.nu_onshellW_p4.Pt()
		self.nu_onshellW_energy[0] = self.nu_onshellW_p4.Energy()
		self.nu_offshellW_eta[0] = self.nu_offshellW_p4.Eta()
		self.nu_offshellW_phi[0] = self.nu_offshellW_p4.Phi()
		self.nu_offshellW_pt[0] = self.nu_offshellW_p4.Pt()
		self.nu_offshellW_energy[0] = self.nu_offshellW_p4.Energy()

		self.onshellW_eta[0] = self.onshellW_p4.Eta()
		self.onshellW_phi[0] = self.onshellW_p4.Phi()
		self.onshellW_pt[0] = self.onshellW_p4.Pt()
		self.onshellW_energy[0] = self.onshellW_p4.Energy()
		self.offshellW_eta[0] = self.offshellW_p4.Eta()
		self.offshellW_phi[0] = self.offshellW_p4.Phi()
		self.offshellW_pt[0] = self.offshellW_p4.Pt()
		self.offshellW_energy[0] = self.offshellW_p4.Energy()

    		self.htoWW_eta[0] = self.htoWW_p4.Eta()
    		self.htoWW_phi[0] = self.htoWW_p4.Phi()
    		self.htoWW_pt[0] = self.htoWW_p4.Pt()
    		self.htoWW_energy[0] = self.htoWW_p4.Energy()
    		self.htoWW_mass[0] = self.htoWW_p4.M()
    		self.htoBB_eta[0] = self.htoBB_p4.Eta()
    		self.htoBB_phi[0] = self.htoBB_p4.Phi()
    		self.htoBB_pt[0] = self.htoBB_p4.Pt()
    		self.htoBB_energy[0] = self.htoBB_p4.Energy()
    		self.htoBB_mass[0] = self.htoBB_p4.M()
    		self.h2tohh_pt[0] = self.h2tohh_p4.Pt()
    		self.h2tohh_energy[0] = self.h2tohh_p4.Energy()
    		self.h2tohh_mass[0] = self.h2tohh_p4.M()

		if (self.h2tohh_p4.Pt()/self.h2tohh_p4.E() <.00000001):
		    print "Strange case: h2tohh pt ", self.h2tohh_p4.Pt(), " energy ",self.h2tohh_p4.E()
		    self.h2tohh_eta[0] = 1000000.0; self.h2tohh_phi[0] = 0.0
		else:
		    self.h2tohh_eta[0] = self.h2tohh_p4.Eta(); self.h2tohh_phi[0] = self.h2tohh_p4.Phi()
		self.hme_h2Mass.Fill(self.h2tohh_mass[0], self.weight[0])
	 	isolution += 1 
	##### end of iteration



