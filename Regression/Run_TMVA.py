#!/usr/bin/python
import os, ROOT, sys
from ROOT import gSystem, gROOT, gApplication, TFile, TTree, TCut, TH1F
from optparse import OptionParser
print "===> Optimizing cuts for selecting Heavy Higgs vs tt."

folder_basic="./files/"
endtag=".root"
signalFile           = [folder_basic+'signal_B1'+endtag]
backgroundFile       = folder_basic+'background_tt'+endtag

presel = 'MMC_h2massweight1_prob>200 && hasRECOjet1 && hasRECOjet1 && hasMET && hastwomuons && (((b1jet_btag&2)>0 && (b2jet_btag&3)>0) || ((b1jet_btag&3)>0 && (b 2jet_btag&2)>0)) && dR_l1l2<3.3 && dR_l1l2>0.07 && dR_b1b2<5. && mass_l1l2<100 && mass_l1l2>5. && mass_b1b2>150 && dR_bl<5 && dR_l1l2b1b2<6 && MINdR_bl<3.2 && MINdR_bl>0.4 && mass_b1b2<700 && mass_trans<250 && MT2<400 && pt_b1b2<300'
cuts  = [' && MMC_h2massweight1_prob>200 && MMC_h2massweight1_prob<380']
Tag       = '_prova_'
outputs   = ['TMVA_B1_RECO'+Tag]
weightDir = ['weights_B1_RECO'+Tag]
MVAS      = ['ALL'] # = ['Likelihood','LikelihoodMIX','KNN','MLP','MLPBFGS','BDT','BDTD']

# Loop on the category to be optimized
for i in range(len(signalFile)):
    # Check consistency:
    print "DOING ITERATION NUMBER: " + str(i)
    if len(signalFile) != len(outputs) or len(signalFile) != len(cuts) or len(signalFile) != len(weightDir) :
        raise RuntimeError('ERROR::Problem in the lenght of your array!')
    # In case you want to test different MVAs in different TMVA jobs
    for iM in range(len(MVAS)):
        # Create Selection and weight folder
        fullcut = presel + cuts[i]
        print "Preselection is: " + fullcut
        OutPut = outputs[i] + "_" + MVAS[iM] + '.root'
        Weight = weightDir[i] + "_" + MVAS[iM]
        print "Output will be in: " + OutPut + ", and xml file in " + str(Weight)
        os.system('rm -rf ' + Weight)
        # Loading files and entries before MVA
        inputSig = TFile.Open( signalFile[i] )
        inputBkg = TFile.Open( backgroundFile )
        treeS    = inputSig.Get("evtree")
        treeB    = inputBkg.Get("evtree")
        hs = TH1F("hs","",100,0.,100)
        hb = TH1F("hb","",100,0.,100)
        treeS.Draw("dR_l1l2>>hs","reweighting * (" + fullcut + ")","goff");
        treeB.Draw("dR_l1l2>>hb","reweighting * (" + fullcut + ")","goff");
        NEntries_S = hs.Integral()
        NEntries_B = hb.Integral()
        # Get the corresponding weights
        hW = TH1F("hW","", 1000, 0., 50.) 
        treeS.Draw("weight>>hW",fullcut,"goff")
        signalWeight =  hW.GetMean()
        treeB.Draw("weight>>hW",fullcut,"goff")
        backgroundWeight = hW.GetMean()
        print "---> Signal: " + str(NEntries_S*signalWeight) + " the weight used is: " + str(signalWeight)
        print "---> Backgr: " + str(NEntries_B*backgroundWeight) + " the weight used is: " + str(backgroundWeight)
        #Now the MVA
        print "Now starting serious things..."
        if( MVAS[iM]=="ALL" and len(MVAS)==1 ):
            print 'Executing: python CutsOptimization.py -a "' + str(fullcut) + '" -o ' + str(OutPut) + ' -i ' + str(signalFile[i]) + ' -j ' + str(backgroundFile) + ' -w ' + str(Weight)
            os.system('python CutsOptimization.py -a "' + fullcut + '" -o ' + OutPut + ' -i ' + signalFile[i] + ' -j ' + backgroundFile + ' -w ' + str(Weight) )
        elif( MVAS[iM]=="ALL" and len(MVAS)!=1 ):
            print "WARNING:: you want to run an all method and you have len(MVAS)!=1"; exit;
        else:
            print 'Executing: python CutsOptimization.py -a "' + str(fullcut) + '" -o ' + str(OutPut) + ' -i ' + str(signalFile[i]) + ' -j ' + str(backgroundFile) + ' -m ' + str(MVAS[iM]) + ' -w ' + str(Weight)
            os.system('python CutsOptimization.py -a "' + fullcut + '" -o ' + OutPut + ' -i ' + signalFile[i] + ' -j ' + backgroundFile + ' -m ' + MVAS[iM] + ' -w ' + str(Weight) )
print "DONE!!\n \n"
