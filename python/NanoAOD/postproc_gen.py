#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor

#from HhhAnalysis.NanoAOD.HHbbWWProducer import *
from genParticleProducerSL import genHH

filesSignal = [
  "/eos/home-r/rrabadan/CMS/Hhh/mc/GluGluToRadionToHHTo2B2WToLNu2J_M-280_narrow_13TeV/tree_1.root",
  "/eos/home-r/rrabadan/CMS/Hhh/mc/GluGluToRadionToHHTo2B2WToLNu2J_M-280_narrow_13TeV/tree_2.root",
  #"/eos/home-r/rrabadan/CMS/Hhh/mc/GluGluToRadionToHHTo2B2WToLNu2J_M-280_narrow_13TeV/tree_3.root",
  #"/eos/home-r/rrabadan/CMS/Hhh/mc/GluGluToRadionToHHTo2B2WToLNu2J_M-280_narrow_13TeV/tree_4.root",
  #"/eos/home-r/rrabadan/CMS/Hhh/mc/GluGluToRadionToHHTo2B2WToLNu2J_M-280_narrow_13TeV/tree_5.root"
]

modules = [genHH()]

p=PostProcessor(".", filesSignal, "1","keep_and_drop_gen.txt", modules, provenance=True, jsonInput=None, maxEntries=20000)
p.run()
