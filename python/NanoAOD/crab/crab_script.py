#!/usr/bin/env python
import os
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import * 

from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jecUncertainties import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetUncertainties import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import *

#this takes care of converting the input files from CRAB
from PhysicsTools.NanoAODTools.postprocessing.framework.crabhelper import inputFiles,runsAndLumis

from  PhysicsTools.NanoAODTools.postprocessing.examples.mhtProducer import *
mht_hh = lambda : mhtProducer( lambda j : j.pt > 20 and abs(j.eta) < 2.4,
                            lambda mu : mu.pt > 10 and abs(mu.eta) < 2.4,
                            lambda el : el.pt > 10 and abs(el.eta) < 2.5 )

jsonfile = None
##cp leptonSF/cMVAv2_Moriond17_B_H.csv   ../../../PhysicsTools/NanoAODTools/data/btagSF/
btagSF2016_cMVA = lambda : btagSFProducer("2016",  algo = 'cmva', sfFileName='cMVAv2_Moriond17_B_H.csv')
#btagSF2017_cMVA = lambda : btagSFProducer("2017",  algo = 'cmva')## for 2017Data?
modules = [puAutoWeight(), btagSF2016_cMVA(), mht_hh(), HHbbWWProducer(True, verbose = 1) ]
p=PostProcessor(".",inputFiles(),"1", modules=modules, friend = True, provenance=True,fwkJobReport=True,jsonInput=jsonfile)
p.run()

print "DONE"
os.system("ls -lR")

