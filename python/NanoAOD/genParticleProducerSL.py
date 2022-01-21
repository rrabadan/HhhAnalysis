import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import math
import logging


sign = lambda x: int(math.copysign(1, x) if x != 0 else 0)


statusFlagsMap = {
  # comments taken from:
  # DataFormats/HepMCCandidate/interface/GenParticle.h
  # PhysicsTools/HepMCCandAlgos/interface/MCTruthHelper.h
  #
  # nomenclature taken from:
  # PhysicsTools/NanoAOD/python/genparticles_cff.py
  #
  #TODO: use this map in other gen-lvl particle selectors as well
  # GenLepFromTauFromTop -> isDirectPromptTauDecayProduct &&
  #                         isDirectHardProcessTauDecayProduct &&
  #                         isLastCopy &&
  #                         ! isDirectHadronDecayProduct
  # GenLepFromTau -> isDirectTauDecayProduct (or isDirectPromptTauDecayProduct?) &&
  #                  isLastCopy &&
  #                  ! isDirectHadronDecayProduct
  #                  (&& maybe isHardProcessTauDecayProduct?)
  #
  # GenLepFromTop -> isPrompt &&
  #                  isHardProcess &&
  #                  (isLastCopy || isLastCopyBeforeFSR) &&
  #                  ! isDirectHadronDecayProduct
  #
  # Not sure if to choose (isLastCopy or isLastCopyBeforeFSR) or just isFirstCopy:
  # GenWZQuark, GenHiggsDaughters, GenVbosons
  #
  # Have no clue what exactly to require from GenTau
  #
  #
  'isPrompt'                           : 0,  # any decay product NOT coming from hadron, muon or tau decay
  'isDecayedLeptonHadron'              : 1,  # a particle coming from hadron, muon, or tau decay
                                             # (does not include resonance decays like W,Z,Higgs,top,etc)
                                             # equivalent to status 2 in the current HepMC standard
  'isTauDecayProduct'                  : 2,  # a direct or indirect tau decay product
  'isPromptTauDecayProduct'            : 3,  # a direct or indirect decay product of a prompt tau
  'isDirectTauDecayProduct'            : 4,  # a direct tau decay product
  'isDirectPromptTauDecayProduct'      : 5,  # a direct decay product from a prompt tau
  'isDirectHadronDecayProduct'         : 6,  # a direct decay product from a hadron
  'isHardProcess'                      : 7,  # part of the hard process
  'fromHardProcess'                    : 8,  # the direct descendant of a hard process particle of the same pdg id
  'isHardProcessTauDecayProduct'       : 9,  # a direct or indirect decay product of a tau from the hard process
  'isDirectHardProcessTauDecayProduct' : 10, # a direct decay product of a tau from the hard process
  'fromHardProcessBeforeFSR'           : 11, # the direct descendant of a hard process particle of the same pdg id
                                             # for outgoing particles the kinematics are those before QCD or QED FSR
  'isFirstCopy'                        : 12, # the first copy of the particle in the chain with the same pdg id
  'isLastCopy'                         : 13, # the last copy of the particle in the chain with the same pdg id
                                             # (and therefore is more likely, but not guaranteed,
                                             # to carry the final physical momentum)
  'isLastCopyBeforeFSR'                : 14, # the last copy of the particle in the chain with the same pdg id
                                             # before QED or QCD FSR (and therefore is more likely,
                                             # but not guaranteed, to carry the momentum after ISR;
                                             # only really makes sense for outgoing particles
}

def printGenParticles(genparts, genParticles, string=''):
    print "debugging ",string
    genParticles = list(genParticles)
    for genp in genparts:
        idx = None
        if hasattr(genp, 'idx'):
	   idx = genp.idx
        print "---- id ",genp.pdgId," idx ",idx," pt ",genp.pt," eta ",genp.eta," mass ",genp.mass, " status ",genp.status, " its motheridx ",genp.genPartIdxMother
        if genp.genPartIdxMother >0:
	    genp_mother = genParticles[genp.genPartIdxMother]
	    print "\t its  mother id ",genp_mother.pdgId, " mass ",genp_mother.mass, " status ",genp_mother.status

class MassTable:
  def __init__(self):
    self.pdgTable = ROOT.TDatabasePDG()

  def getMass(self, mass, pdgId):
    if mass > 10. or (pdgId == 22 and mass > 1.) or abs(pdgId) == 24 or pdgId == 23:
      return mass
    else:
      genParticleInstance = self.pdgTable.GetParticle(pdgId)
      if not genParticleInstance:
        # Since most of the common low-mass particles are defined in ROOT's PDG table,
        # and that it's more than likely we don't need such generator-level information,
        # we can safely set the masses of such particles to 0 GeV
        logging.debug("Setting the mass to 0 GeV for a particle with PDG id of %d" % pdgId)
        return 0.
      return genParticleInstance.Mass()

  def getCharge(self, pdgId):
    genParticleInstance = self.pdgTable.GetParticle(pdgId)
    if not genParticleInstance:
      # It's more than likely that we don't need to know the charges of generator-level particles
      # that are not defined in ROOT's PDG id table. Therefore, we assign neutral charges to
      # these particles.
      logging.debug("Setting the charge to neutral for a particle with PDG id of %d" % pdgId)
      return 0
    return sign(genParticleInstance.Charge())


class GenPartAux:
  def __init__(self, genPart, idx, massTable):
    self.pt               = genPart.pt
    self.eta              = genPart.eta
    self.phi              = genPart.phi
    self.mass             = massTable.getMass(genPart.mass, genPart.pdgId)
    self.pdgId            = genPart.pdgId
    self.charge           = massTable.getCharge(genPart.pdgId)
    self.status           = genPart.status
    self.statusFlags      = genPart.statusFlags
    self.genPartIdxMother = genPart.genPartIdxMother
    self.idx              = idx

  def __str__(self):
    return "pt = %.3f eta = %.3f phi = %.3f mass = %.3f pdgId = %i charge = %i status = %i " \
           "statusFlags = %i mom = %i idx = %i" % \
      (self.pt, self.eta, self.phi, self.mass, self.pdgId, self.charge, self.status, \
       self.statusFlags, self.genPartIdxMother, self.idx)

  def __repr__(self):
    return self.__str__()

  def checkIf(self, condition):
    assert(condition in statusFlagsMap)
    return (self.statusFlags & (1 << statusFlagsMap[condition]) != 0)


class SelectionOptions:
  SAVE_TAU                      = 0
  SAVE_LEPTONIC_TAU             = 1
  SAVE_HADRONIC_TAU             = 2
  SAVE_LEPTON_FROM_TAU          = 3
  SAVE_LEPTONIC_NU_FROM_TAU     = 4
  SAVE_TAU_NU_FROM_LEPTONIC_TAU = 5
  SAVE_TAU_NU_FROM_HADRONIC_TAU = 6

  SAVE_TOP                               = 9
  SAVE_W_FROM_TOP                        = 10
  SAVE_BQUARK_FROM_TOP                   = 11
  SAVE_LEPTON_FROM_W_FROM_TOP            = 12
  SAVE_LEPTONIC_NU_FROM_W_FROM_TOP       = 13
  SAVE_TAU_FROM_TOP                      = 14
  SAVE_TAU_NU_FROM_TOP                   = 15
  SAVE_LEPTON_FROM_TAU_FROM_TOP          = 16
  SAVE_LEPTON_NU_FROM_TAU_FROM_TOP       = 17
  SAVE_TAU_NU_FROM_LEPTONIC_TAU_FROM_TOP = 18
  SAVE_TAU_NU_FROM_HADRONIC_TAU_FROM_TOP = 19
  SAVE_NU_FROM_TAU_FROM_TOP              = 20
  SAVE_QUARK_FROM_W_FROM_TOP             = 21

class SelectionXToHHTobbWWTobblvqqOptions:
  SAVE_HIG_FROM_X                  = 0
  SAVE_BQUARK_FROM_HIG             = 1
  SAVE_W1_FROM_HIG                 = 2
  SAVE_W2_FROM_HIG                 = 3
  SAVE_LEPTON_FROM_W_FROM_HIG      = 4
  SAVE_LEPTONIC_NU_FROM_W_FROM_HIG = 5
  SAVE_QUARK_FROM_W_FROM_HIG       = 6
  SAVE_QUARKBAR_FROM_W_FROM_HIG    = 7

def genLeptonSelection(genParticles):
  return filter(lambda genPart: abs(genPart.pdgId) in [11, 13] and genPart.status == 1, genParticles)

def genPromptLeptonSelection(genParticles):
  return filter(
    lambda genLepton:
      genLepton.checkIf('isLastCopy') and
      not genLepton.checkIf('isDirectHadronDecayProduct') and
      (
        genLepton.checkIf('isPrompt') or
        genLepton.checkIf('isDirectPromptTauDecayProduct')
      ),
    genLeptonSelection(genParticles)
  )

def genHiggsSelection(genParticles):
  return filter(
    lambda genPart:
      genPart.pdgId == 25 and \
      (genParticles[genPart.genPartIdxMother].pdgId != 25 if genPart.genPartIdxMother >= 0 else True),
    genParticles
  )

def genRadionGravitonSelection(genParticles): ## X->HH
  #print "genRadionGravitonSelection type(genParticles) ",type(genParticles)
  XCandidates = []
  pdgIds =[35, 39]## add potential id for X
     
  XCandidates = filter (lambda genPart: abs(genPart.pdgId) in pdgIds, genParticles)
  #printGenParticles(XCandidates, genParticles, "genRadionGravitonSelection before final filter")

  return filter(
    lambda genPart:
      (genParticles[genPart.genPartIdxMother].pdgId != genPart.pdgId if genPart.genPartIdxMother >= 0 else True),
    XCandidates
  )

def genHiggsDaughtersSelection(genParticles):
  return filter(
    lambda genPart:
      genPart.pdgId != 25 and \
      (genParticles[genPart.genPartIdxMother].pdgId == 25 if genPart.genPartIdxMother >= 0 else False),
    genParticles
  )

def genWZquarkSelection(genParticles):
  return filter(
    lambda genPart:
      abs(genPart.pdgId) in [1, 2, 3, 4, 5, 6] and genPart.genPartIdxMother >= 0 and \
      abs(genParticles[genPart.genPartIdxMother].pdgId) in [23, 24],
    genParticles
  )

def genVbosonSelection(genParticles):
  return filter(
    lambda genPart:
      abs(genPart.pdgId) in [23, 24] and genPart.genPartIdxMother >= 0 and \
      genParticles[genPart.genPartIdxMother].pdgId != genPart.pdgId,
    genParticles
  )

def genNuSelection(genParticles):
  return filter(lambda genPart: abs(genPart.pdgId) in [12, 14, 16], genParticles)

def findLastDecendantWithSameId(cand, genParticles):
  nextgeneration = filter(lambda genPart: genPart.pdgId == cand.pdgId and genPart.genPartIdxMother == cand.idx, genParticles)
  if len(nextgeneration) == 1:
      return findLastDecendantWithSameId(nextgeneration[0])
  elif len(nextgeneration) == 0:
      return cand
  else:
      raise ValueError("Invalid number of nextgeneration(%s) of cand (%s) in findLastDecendantWithSameId %i:" % \
	(', '.join(map(lambda genPart: str(genPart), nextgeneration)), str(cand), len(nextgeneration))
      )

def genXToHHTobbWWTobblvqqSelection(genParticles, choice, enable_consistency_checks = True):
  genBQuarkFromHIG = []
  genW1FromHIG = None
  genW2FromHIG = None
  genLepFromWFromHIG = None
  genNuFromWFromHIG = None
  genQuarkFromWFromHIG = None
  debugXToHHTobbWW = True
  #XCandidates  = genRadionGravitonSelection(genParticles)
  XCandidates = []
  pdgIds =[35, 39]## add potential id for X
     
  XCandidates = filter (lambda genPart: abs(genPart.pdgId) in pdgIds, genParticles)
  XCandidatesidx = [X.idx  for X in XCandidates]
  finalXCandidates = filter(lambda X : X.genPartIdxMother not in XCandidatesidx, XCandidates)
     
  if len(finalXCandidates) != 1: 
    #print "Xcandidates is not one ",len(finalXCandidates)
    #printGenParticles(finalXCandidates, genParticles, "genXToHHTobbWWSelection")
    ##avoid printout for non-signal MC
    return []
  
  ##step1 find X->HH
  genHiggs = filter(lambda genPart: genPart.pdgId == 25, genParticles)
  genHiggsidx = [Higgs.idx for Higgs in genHiggs]
  #printGenParticles(genHiggs, genParticles, "genHiggs")
  genHiggsFromX = filter(lambda genPart: genPart.genPartIdxMother in XCandidatesidx, genHiggs)
  if len(genHiggsFromX) != 2:
    raise ValueError("Invalid number of X->HH (%s) decay products (%s): %i" % \
          (XCandidates[0], ', '.join(map(lambda genPart: str(genPart), genHiggsFromX)), len(genHiggsFromX))
    )

  if choice == SelectionXToHHTobbWWTobblvqqOptions.SAVE_HIG_FROM_X:
    return genHiggsFromX

  ##step2 find HH->bbWW
  ##possible combination: HH->bbbb; HH->WWWW;  HH->bbbb?
  ##possible HH->bbZZ
  ## no intermediate state in NanoAOD sample ?
  genBQuarkFromHiggs = filter(lambda genPart: abs(genPart.pdgId) == 5 and genPart.genPartIdxMother in genHiggsidx, genParticles)
  genWFromHiggs = filter(lambda genPart: abs(genPart.pdgId) == 24 and genPart.genPartIdxMother in genHiggsidx, genParticles)
  genZFromHiggs = filter(lambda genPart: abs(genPart.pdgId) == 23 and genPart.genPartIdxMother in genHiggsidx, genParticles)
  if not(len(genBQuarkFromHiggs) == 2 and len(genWFromHiggs) == 2) :
    if len(genZFromHiggs) != 2 and debugXToHHTobbWW:
      print "Invalid number of H->bb (%s):%d or H->WW (%s):%d in HH (%s) decay" % \
        (', '.join(map(lambda genPart: str(genPart), genBQuarkFromHiggs)), len(genBQuarkFromHiggs), \
         ', '.join(map(lambda genPart: str(genPart), genWFromHiggs)), len(genWFromHiggs), \
         ', '.join(map(lambda genPart: str(genPart), genHiggsFromX)))
      #decendantsFromHiggs = filter(lambda genPart: genPart.genPartIdxMother in genHiggsidx, genParticles)
      #printGenParticles(decendantsFromHiggs, genParticles, "All decendants from HH")
    return []
  
  if choice == SelectionXToHHTobbWWTobblvqqOptions.SAVE_BQUARK_FROM_HIG:
    if len(genBQuarkFromHiggs)  != 2 or genBQuarkFromHiggs[0].genPartIdxMother != genBQuarkFromHiggs[1].genPartIdxMother:
      if debugXToHHTobbWW:
        print ("Invalid number of H->bb decay products (%s): %i" % \
        ( ', '.join(map(lambda genPart: str(genPart), genBQuarkFromHiggs)), len(genBQuarkFromHiggs)))
      return []
    return genBQuarkFromHiggs

  if len(genWFromHiggs)  != 2 or genWFromHiggs[0].genPartIdxMother != genWFromHiggs[1].genPartIdxMother:
    if debugXToHHTobbWW:
      print ("Invalid number of H->WW decay products (%s): %i" % \
      ( ', '.join(map(lambda genPart: str(genPart), genWFromHiggs)), len(genWFromHiggs)))
    return []
       
  genWFromHiggs.sort(key=lambda x:x.mass,reverse=True)	
  genW1FromHIG = genWFromHiggs[0]
  genW2FromHIG = genWFromHiggs[1]
  if choice == SelectionXToHHTobbWWTobblvqqOptions.SAVE_W1_FROM_HIG:
      return [genW1FromHIG]
  if choice == SelectionXToHHTobbWWTobblvqqOptions.SAVE_W2_FROM_HIG:
      return [genW2FromHIG]
      
  ##step3 find W->lv, W->qq'
  def findDaughtersFromW(genW):
    genWDaughters = filter(lambda genPart: genPart.genPartIdxMother == genW.idx, genParticles)
    daughters = {}
    if len(genWDaughters) != 2:
      raise ValueError("Invalid number (%i) of W (%s) daughters decay: %s" % \
        (len(genWDaughters), genW, ', '.join(map(str, genWDaughters)))
      )

    if any(map(lambda genPart: abs(genPart.pdgId) in [11, 13, 15], genWDaughters)):
      lepsFromW = filter(lambda genPart: -sign(genW.pdgId) * genPart.pdgId in [11, 13, 15], genWDaughters)
      if len(lepsFromW) != 1:
        raise ValueError("Inconsistent W (%s) decay products decay: %s" % \
          (genW, ', '.join(map(str, genWDaughters)))
        )
      genLepFromW = lepsFromW[0]
      nusLepFromW = filter(lambda genPart: genPart.pdgId == sign(genW.pdgId) * (abs(genLepFromW.pdgId) + 1), genWDaughters)
      if len(nusLepFromW) != 1:
        raise ValueError("Inconsistent W (%s) decay products decay: %s" % \
          (genW, ', '.join(map(str, genWDaughters)))
        )
      daughters['Lep'] =  [genLepFromW]
      daughters['Nu']  =  [nusLepFromW[0]]
      return daughters
    elif all(map(lambda genPart: abs(genPart.pdgId) in [1, 2, 3, 4, 5], genWDaughters)):
      quarksFromWFromH_IdSorted = list(sorted(genWDaughters, key = lambda genPart: abs(genPart.pdgId), reverse = True))
      daughters['Qbar'] = [quarksFromWFromH_IdSorted[0]]
      daughters['Q'] = [quarksFromWFromH_IdSorted[1]]
      return daughters
    else:
      raise ValueError("Not implemented W decay.")

  #genW1Daughters = filter(lambda genPart: genPart.genPartIdxMother == genW1FromHIG.idx, genParticles)
  #genW2Daughters = filter(lambda genPart: genPart.genPartIdxMother == genW2FromHIG.idx, genParticles)

  W1Daughters = findDaughtersFromW(genW1FromHIG)
  W2Daughters = findDaughtersFromW(genW2FromHIG)

  #print(W1Daughters)
  #print(W2Daughters)
  
  #if choice == SelectionXToHHTobbWWTobblvqqOptions.SAVE_W1_FROM_HIG:
  #    return [genW1FromHIG]
  #if choice == SelectionXToHHTobbWWTobblvqqOptions.SAVE_W2_FROM_HIG:
  #    return [genW2FromHIG]

  WDaughters = W1Daughters.copy()
  WDaughters.update(W2Daughters)

  if choice in [ SelectionXToHHTobbWWTobblvqqOptions.SAVE_LEPTON_FROM_W_FROM_HIG ,\
    SelectionXToHHTobbWWTobblvqqOptions.SAVE_LEPTONIC_NU_FROM_W_FROM_HIG ]:
    #genW1Daughters = filter(lambda genPart: genPart.genPartIdxMother == genW1FromHIG.idx, genParticles)
    #W1Daughters = findLepandNufromW(genW1FromHIG)
    #W1Daughters = findLepandNufromW(findLastDecendantWithSameId(genW1FromHIG, genParticles))
    if choice == SelectionXToHHTobbWWTobblvqqOptions.SAVE_LEPTON_FROM_W_FROM_HIG:
        return WDaughters['Lep']
    else:
        return WDaughters['Nu']
    
  elif choice in [ SelectionXToHHTobbWWTobblvqqOptions.SAVE_QUARK_FROM_W_FROM_HIG ,\
    SelectionXToHHTobbWWTobblvqqOptions.SAVE_QUARKBAR_FROM_W_FROM_HIG ]:
    #W2Daughters = findLepandNufromW(genW2FromHIG)
    #W2Daughters = findLepandNufromW(findLastDecendantWithSameId(genW2FromHIG, genParticles))
    if choice == SelectionXToHHTobbWWTobblvqqOptions.SAVE_QUARK_FROM_W_FROM_HIG:
        return WDaughters['Q']
    else:
        return WDaughters['Qbar']
  else:
    raise ValueError("Choice %i not implemented in  XToHHTobbWWTobblvqq" % choice)

class genParticleProducer(Module):

  def __init__(self, genEntry, verbose = False):
    self.massTable = MassTable()
    self.branchLenNames  = {}
    self.selections      = {}
    self.branchBaseNames = []

    self.genBranches = {
        "pt"          : "F",
        "eta"         : "F",
        "phi"         : "F",
        "mass"        : "F",
        "pdgId"       : "I",
        "charge"      : "I",
        "status"      : "I",
        "statusFlags" : "I",
        "idx"         : "I",
      }

    for branchBaseName, selection in genEntry.items():
      self.branchBaseNames.append(branchBaseName)
      self.selections[branchBaseName]     = selection
      self.branchLenNames[branchBaseName] = "n%s" % branchBaseName

    if verbose:
      logging.getLogger().setLevel(logging.DEBUG)

  def beginJob(self):
    pass

  def endJob(self):
    pass

  def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    self.out = wrappedOutputTree
    for branchBaseName in self.branchBaseNames:
      for branchName, branchType in self.genBranches.items():
        self.out.branch(
          "%s_%s" % (branchBaseName, branchName),
          branchType,
          lenVar = self.branchLenNames[branchBaseName]
        )

  def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    pass

  def analyze(self, event):
    genParticles  = map(lambda genPartIdx: GenPartAux(genPartIdx[1], genPartIdx[0], self.massTable), enumerate(Collection(event, "GenPart")))
    allgenParticles = Collection(event, "GenPart")
    #print " type(genParticles) ",type(genParticles)  ," type(allgenParticles) ",type(allgenParticles)

    for branchBaseName in self.branchBaseNames:
      gen_arr = self.selections[branchBaseName](genParticles)
      gen_arr = list(sorted(gen_arr, key = lambda genPart: genPart.pt, reverse = True)) # sort by pT
      #printGenParticles(gen_arr, allgenParticles,branchBaseName)

      for branchName, branchType in self.genBranches.items():
        self.out.fillBranch(
          "%s_%s" % (branchBaseName, branchName),
          map(lambda genPart: getattr(genPart, branchName), gen_arr)
        )
    return True


genLeptonEntry = ("GenLep", genPromptLeptonSelection)
genLeptonAllEntry = ("GenLepAll", genLeptonSelection)
genXEntry = ("GenX", genRadionGravitonSelection)
genHiggsEntry = ("GenHiggs", genHiggsSelection)
genHiggsDaughtersEntry = ("GenHiggsDaughters", genHiggsDaughtersSelection)
genNuEntry = ("GenNu", genNuSelection)
genWZquarkEntry = ("GenWZQuark", genWZquarkSelection)
genVbosonEntry = ("GenVbosons", genVbosonSelection)
genBQuarkFromHIGEntry = ("GenBQuarkFromHiggs", (lambda genParticles : genXToHHTobbWWTobblvqqSelection(genParticles, SelectionXToHHTobbWWTobblvqqOptions.SAVE_BQUARK_FROM_HIG)))
genW1FromHIGEntry = ("GenW1FromHiggs", (lambda genParticles : genXToHHTobbWWTobblvqqSelection(genParticles, SelectionXToHHTobbWWTobblvqqOptions.SAVE_W1_FROM_HIG)))
genW2FromHIGEntry = ("GenW2FromHiggs", (lambda genParticles : genXToHHTobbWWTobblvqqSelection(genParticles, SelectionXToHHTobbWWTobblvqqOptions.SAVE_W2_FROM_HIG)))
genLepFromWFromHIGEntry = ("GenLepFromWFromHiggs", (lambda genParticles : genXToHHTobbWWTobblvqqSelection(genParticles, SelectionXToHHTobbWWTobblvqqOptions.SAVE_LEPTON_FROM_W_FROM_HIG)))
genNuFromWFromHIGEntry = ("GenNuFromWFromHiggs",  (lambda genParticles : genXToHHTobbWWTobblvqqSelection(genParticles, SelectionXToHHTobbWWTobblvqqOptions.SAVE_LEPTONIC_NU_FROM_W_FROM_HIG)))
genQuarkFromWFromHIGEntry = ("GenQuarkFromWFromHiggs", (lambda genParticles : genXToHHTobbWWTobblvqqSelection(genParticles, SelectionXToHHTobbWWTobblvqqOptions.SAVE_QUARK_FROM_W_FROM_HIG)))
genQuarkBarFromWFromHIGEntry = ("GenQuarkBarFromWFromHiggs",  (lambda genParticles : genXToHHTobbWWTobblvqqSelection(genParticles, SelectionXToHHTobbWWTobblvqqOptions.SAVE_QUARKBAR_FROM_W_FROM_HIG)))

# provide these variables as the 2nd arguments to the import option for the nano_postproc.py script
genLepton                      = lambda : genParticleProducer(dict([genLeptonEntry]))                      # all prompt stable leptons
genLeptonAll                   = lambda : genParticleProducer(dict([genLeptonAllEntry]))                   # all stable leptons
genHiggs                       = lambda : genParticleProducer(dict([genHiggsEntry]))                       # all Higgs (first in the decay chain)
genHiggsDaughters              = lambda : genParticleProducer(dict([genHiggsDaughtersEntry]))              # all Higgs daughters
genBQuarkFromHiggs             = lambda : genParticleProducer(dict([genBQuarkFromHIGEntry]))
genW1FromHiggs                 = lambda : genParticleProducer(dict([genW1FromHIGEntry]))
genW2FromHiggs                 = lambda : genParticleProducer(dict([genW2FromHIGEntry]))
genLepFromWFromHiggs           = lambda : genParticleProducer(dict([genLepFromWFromHIGEntry]))
genNuFromWFromHiggs            = lambda : genParticleProducer(dict([genNuFromWFromHIGEntry]))
genQuarkFromWFromHiggs         = lambda : genParticleProducer(dict([genQuarkFromWFromHIGEntry]))
genQuarkBarFromWFromHiggs      = lambda : genParticleProducer(dict([genQuarkBarFromWFromHIGEntry]))
#genTau                         = lambda : genParticleProducer(dict([genTauEntry]))                         # all taus
#genNu                          = lambda : genParticleProducer(dict([genNuEntry]))                          # all neutrinos
#genWZquark                     = lambda : genParticleProducer(dict([genWZquarkEntry]))                     # all quarks coming from W or Z decay
#genVboson                      = lambda : genParticleProducer(dict([genVbosonEntry]))                      # all W and Z bosons (first in the decay chain)

genAll = lambda : genParticleProducer(dict([
    genLeptonEntry,
    genLeptonAllEntry,
    genXEntry,
    genHiggsEntry,
    #genHiggsDaughtersEntry,
    genBQuarkFromHIGEntry,
    genW1FromHIGEntry,
    genW2FromHIGEntry,
    #genLepFromW1FromHIGEntry,
    #genNuFromW1FromHIGEntry,
    #genLepFromW2FromHIGEntry,
    #genNuFromW2FromHIGEntry,
    #genNuEntry,
    #genWZquarkEntry,
    #genVbosonEntry,
  ]))

genHH = lambda : genParticleProducer(dict([
    genXEntry,
    genHiggsEntry,
    genBQuarkFromHIGEntry,
    genW1FromHIGEntry,
    genW2FromHIGEntry,
    genLepFromWFromHIGEntry,
    genNuFromWFromHIGEntry,
    genQuarkFromWFromHIGEntry,
    genQuarkBarFromWFromHIGEntry,
  ]))
