import ROOT
import string
from math import sqrt, pi, degrees
import os
import json
import re
import numpy as np
from PhysicsTools.NanoAODTools.postprocessing.tools import * #deltaR

#Information about these cuts can be found at
#https://gitlab.cern.ch/cms-hh-bbww/cms-hh-to-bbww/-/blob/master/Legacy/objects.md#hh-h-bbww-analysis-documentation-with-full-run2-data

def btagDeepFalvB_MediumWP(runyear):
    if runyear == 2016: return 0.3093
    if runyear == 2017: return 0.3033
    if runyear == 2018: return 0.2770

def ak4jetIdCut(runyear):
    #Require loose for 2016, tight for 2017 and 2018
    if runyear == 2016: return 0
    if runyear == 2017: return 2
    if runyear == 2018: return 2

def ak8jetIdCut(runyear):
    #Require loose for 2016, tight for 2017 and 2018
    if runyear == 2016: return 0
    if runyear == 2017: return 2
    if runyear == 2018: return 2

def ak8lsjetIdCut(runyear):
    #Require loose for 2016, tight for 2017 and 2018
    if runyear == 2016: return 0
    if runyear == 2017: return 2
    if runyear == 2018: return 2

def ak8subjetMediumBtag(runyear):
    if runyear == 2016: return 0.6321
    if runyear == 2017: return 0.4941
    if runyear == 2018: return 0.4184

def conept(lep):
    """
    #https://github.com/CERN-PH-CMG/cmgtools-lite/blob/f8a34c64a4489d94ff9ac4c0d8b0b06dad46e521/TTHAnalysis/python/tools/conept.py#L74
    if (abs(lep.pdgId) != 11 and abs(lep.pdgId) != 13): return lep.pt
    if (abs(lep.pdgId) != 13 or lep.mediumMuonId > 0) and lep.mvaTTH > 0.90: return lep.pt
    return lep.pt #Currently do not have jetPtRatiov2, looking for a fix
    #else: return 0.90* lep.pt / lep.jetPtRatiov2
    """
    #Fix taken from Florian
    #https://github.com/FlorianBury/HHbbWWAnalysis/blob/84204a8d8c31eb67d7f1a8e4bd77ce00d7232bd6/BaseHHtobbWW.py#L951-L962Le
    if (abs(lep.pdgId) != 11 and abs(lep.pdgId) != 13): return lep.pt
    elif (abs(lep.pdgId) == 11 and lep.mvaTTH > 0.30): return lep.pt
    elif (abs(lep.pdgId) == 13 and lep.mediumId and lep.mvaTTH > 0.50): return lep.pt
    else: return 0.9 * lep.pt * (1.0 + lep.jetRelIso)



def get_jet_from_lepton(lep, jets, runyear):
    jetId = lep.jetIdx
    if jetId < 0 or jetId > len(jets):
        return 0
    return jets[jetId]

def leptonCleaning(ele, muons, runyear):
    for mu in muons:
      if deltaR(ele.eta, ele.phi, mu.eta, mu.phi) < 0.3:
        return False
    return True

def muonPreselection(muon, runyear):
    return abs(muon.eta)<2.4 and muon.pt>5 and abs(muon.dxy) <= 0.05 and abs(muon.dz) <= 0.1 and muon.looseId and muon.miniPFRelIso_all < 0.4 and muon.sip3d < 8

def muonFakeable(muon, runyear):
    jetDeepJet_cut = False
    lepMVA_cut = True
    if muon.mvaTTH <= 0.85:
        lepMVA_cut = False
        if muon.jetRelIso <= 0.5 and 
    return conept(muon) >= 10 and jetDeepJet_cut and lepMVA_cut

def muonTight(muon, runyear):
    return True

def electronPreselection(ele, runyear):
    return abs(ele.eta)<2.5 and ele.pt>7 and abs(ele.dxy) <= 0.05 and abs(ele.dz) <= 0.1 and ele.miniPFRelIso_all < 0.4 and ele.lostHits <=1 and ele.mvaFall17V2noIso_WPL and ele.sip3d < 8

def electronFakeable(ele, jets, runyear):
    sieie_cut = False
    if ele.eta > 1.479 and ele.eta < 2.5:
        if ele.sieie <= 0.030:
            sieie_cut = True
    if ele.eta < 1.479:
        if ele.sieie <= 0.011:
            sieie_cut = True
    jetDeepJet_cut = False
    if get_jet_from_lepton(ele, jets, runyear) != 0:
        if get_jet_from_lepton(ele, jets, runyear).btagDeepFlavB <= btagDeepFalvB_MediumWP(runyear):
            jetDeepJet_cut = True
    lepMVA_cut = True
    if ele.mvaTTH <= 0.80:
        lepMVA_cut = False
        if ele.jetRelIso <= 0.7 and ele.mvaFall17V2noIso_WP80:
            lepMVA_cut = True
    return conept(ele) >= 10 and sieie_cut and ele.hoe <= 0.10 and ele.eInvMinusPInv >= -0.04 and jetDeepJet_cut and lepMVA_cut and ele.lostHits == 0 and ele.convVeto

def electronTight(ele, runyear):
    return True

def ak4jetPreselection(jet, runyear):
    return abs(jet.eta)<2.4 and jet.pt>25 and jet.jetId >= ak4jetIdCut(runyear)

def ak4jetCleaning(jet, muons, electrons, runyear):
    for mu in muons:
        if deltaR(jet.eta, jet.phi, mu.eta, mu.phi) < 0.4:
            return False
    for ele in electrons:
        if deltaR(jet.eta, jet.phi, ele.eta, ele.phi) < 0.4:
            return False
    return True

def ak4jetBtagging(jet, runyear):
    return True

def ak8jetPreselection(jet, subjets, runyear):
    twoSubJets = False
    subjet1_idx = jet.subJetIdx1
    subjet2_idx = jet.subJetIdx2
    if subjet1_idx < 0 or subjet2_idx < 0 or subjet1_idx >= len(subjets) or subjet2_idx >= len(subjets):
        return False
    subjet1 = subjets[subjet1_idx]
    subjet2 = subjets[subjet2_idx]
    if subjet1.pt > 20 and abs(subjet1.eta) < 2.4 and subjet2.pt > 20 and abs(subjet2.eta) < 2.4:
      twoSubJets = True
    return jet.jetId >= ak8jetIdCut(runyear) and jet.pt >= 200 and abs(jet.eta) <= 2.4 and twoSubJets and jet.msoftdrop > 30 and jet.msoftdrop < 210 and jet.tau2/jet.tau1 <= 0.75

def ak8jetCleaning(jet, muons, electrons, runyear):
    for mu in muons:
        if deltaR(jet.eta, jet.phi, mu.eta, mu.phi) < 0.8:
            return False
    for ele in electrons:
        if deltaR(jet.eta, jet.phi, ele.eta, ele.phi) < 0.8:
            return False
    return True

def ak8jetBtagging(jet, subjets, runyear):
    subjet1_idx = jet.subJetIdx1
    subjet2_idx = jet.subJetIdx2
    if subjet1_idx < 0 or subjet2_idx < 0 or subjet1_idx >= len(subjets) or subjet2_idx >= len(subjets):
        return False
    subjet1 = subjets[subjet1_idx]
    subjet2 = subjets[subjet2_idx]
    subjet_pt_btagging_pass = False
    if (subjet1.btagDeepB > ak8subjetMediumBtag(runyear) and subjet1.pt > 30) or (subjet2.btagDeepB > ak8subjetMediumBtag(runyear) and subjet2.pt > 30):
        subjet_pt_btagging_pass = True
    return subjet_pt_btagging_pass

def ak8lsjetPreselection(jet, subjets, muons, electrons, runyear):
    twoSubJets = False
    subjet1_idx = jet.subJetIdx1
    subjet2_idx = jet.subJetIdx2
    if subjet1_idx < 0 or subjet2_idx < 0 or subjet1_idx >= len(subjets) or subjet2_idx >= len(subjets):
        return False
    subjet1 = subjets[subjet1_idx]
    subjet2 = subjets[subjet2_idx]
    if subjet1.pt > 20 and abs(subjet1.eta) < 2.4 and subjet2.pt > 20 and abs(subjet2.eta) < 2.4:
      twoSubJets = True
    mindR_lep_jet = 999.0
    mindR_lep_subjet = 999.0
    for mu in muons:
        thisdR_lep_jet = deltaR(jet.eta, jet.phi, mu.eta, mu.phi)
        if thisdR_lep_jet < mindR_lep_jet:
          mindR_lep_jet = thisdR_lep_jet
          if thisdR_lep_jet < 1.2:
              mindR_lep_subjet = min(deltaR(subjet1.eta, subjet1.phi, mu.eta, mu.phi), deltaR(subjet2.eta, subjet2.phi, mu.eta, mu.phi))
    for ele in electrons:
        thisdR_lep_jet = deltaR(jet.eta, jet.phi, ele.eta, ele.phi)
        if thisdR_lep_jet < mindR_lep_jet:
          mindR_lep_jet = thisdR_lep_jet
          if thisdR_lep_jet < 1.2:
              mindR_lep_subjet = min(deltaR(subjet1.eta, subjet1.phi, ele.eta, ele.phi), deltaR(subjet2.eta, subjet2.phi, ele.eta, ele.phi))

    return jet.jetId >= ak8lsjetIdCut(runyear) and jet.pt >= 100 and abs(jet.eta) <= 2.4 and twoSubJets and jet.tau2/jet.tau1 <= 0.75 and mindR_lep_jet < 1.2 and  mindR_lep_subjet > 0.1

def ak8lsjetCleaning(jet, ak4jets, ak8jets, runyear):
    return True

#Need to add all the VBF stuff
