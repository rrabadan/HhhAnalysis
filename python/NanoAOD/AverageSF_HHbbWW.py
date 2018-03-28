#!/usr/bin/env python

import json
import sys
from sys import exit
import ROOT
import optparse
import os.path
import string
import logging
import traceback
from array import array
from collections import OrderedDict
import math

import os

from Trigger_cfg import *

def AverageLeg(leg, triggerMaps):
    #parser = optparse.OptionParser(usage)

    #parser.add_option('--outputFile'     , dest='outputFile'     , help='output file with averaged values'           , default='test.txt')
          
    # read default parsing options as well
    #hwwtools.addOptions(parser)
    #hwwtools.loadOptDefaults(parser)
    #(opt, args) = parser.parse_args()


    sys.argv.append( '-b' )
    ROOT.gROOT.SetBatch()


    outputFile = "HLT_"+leg+".text"
    print " outputFile       = ", outputFile

     
      
    print "list of trigger files"
    
    files_names_trigger = {}

    #files_names_trigger["ICHEP2016fullLumi/PartialRuns/HLT_DoubleMuLegHigPt_Run_274094_275000.txt"] = 3.426 
    #files_names_trigger["ICHEP2016fullLumi/PartialRuns/HLT_DoubleMuLegHigPt_Run_275001_275783.txt"] = 3.191 
    #files_names_trigger["ICHEP2016fullLumi/PartialRuns/HLT_DoubleMuLegHigPt_Run_275784_276500.txt"] = 3.191 
    #files_names_trigger["ICHEP2016fullLumi/PartialRuns/HLT_DoubleMuLegHigPt_Run_276501_276811.txt"] = 3.191 
    #opt.outputFile = "ICHEP2016fullLumi/HLT_DoubleMuLegHigPt.txt"
    print "triggerMaps.keys ",triggerMaps.keys()
    for i in triggerMaps.keys():
	filename = triggerMaps[i]["LegEff"][leg]
	print "i ",i," file ",triggerMaps[i]["LegEff"][leg], " lumi ",triggerMaps[i]["lumi"]
	if filename in files_names_trigger.keys():
	    print "add lumi"
	    files_names_trigger[filename] += triggerMaps[i]["lumi"]
	else:
	    print "just find this file"
	    files_names_trigger[filename] = triggerMaps[i]["lumi"]



    
    print " files_names_trigger = ", files_names_trigger
    
    files_trigger = {}    
    for names in files_names_trigger :
       #print names
       ##print files_names_trigger[names]       
       
       files_trigger [names] = {
         'file' : open (names),
         'lumi' : files_names_trigger[names]
         }

      
   
    list_triggers = {}

    for names in files_trigger :
      list_triggers[names]  = {  
             'values' : [line.rstrip().split() for line in files_trigger[names]['file']        if '#' not in line] ,
             'lumi'   : files_trigger[names]['lumi']
          }
      




    outFile = open(outputFile, 'w')
    outFile.write('# etamin  etamax  ptmin   ptmax   eff     deff_high       deff_low\n')

    first_key = list_triggers.keys()[0] 
    
    for point in list_triggers[first_key]['values']:
      outFile.write( '  {0:.2f}  {1:.2f}  {2:.2f}  {3:.2f} '.format(float(point[0]), float(point[1]), float(point[2]), float(point[3])) )
      
      eff = float(point[4])
      error_eff = float(point[5])
      
      error_eff_up = eff + error_eff
      error_eff_lo = eff - error_eff

      # muons and electrons have provided different formats!
      if len(point) > 6 :
         error_eff_lo = float(point[6])
         error_eff_lo = eff - error_eff_lo                  
      
      luminosity = list_triggers[first_key]['lumi']
     
      all_lumi       =   luminosity
      eff            =   eff          * luminosity
      error_eff_up   =   error_eff_up * luminosity
      error_eff_lo   =   error_eff_lo * luminosity
      
      
      for names in list_triggers :
        if names != first_key:
          for point_other in list_triggers[names]['values']:
            if    float(point_other[0]) == float(point[0]) \
              and float(point_other[1]) == float(point[1]) \
              and float(point_other[2]) == float(point[2]) \
              and float(point_other[3]) == float(point[3]) :

              new_eff = float(point_other[4])
              new_error_eff = float(point_other[5])
              
              new_error_eff_up = new_eff + new_error_eff
              new_error_eff_lo = new_eff - new_error_eff

              # muons and electrons have provided different formats!
              if len(point_other) > 6 :
                 new_error_eff_lo = float(point_other[6])
                 new_error_eff_lo = new_eff - new_error_eff_lo                  
        
        
              luminosity = list_triggers[names]['lumi']

              all_lumi       =   all_lumi     + luminosity
              eff            =   eff          +  new_eff          * luminosity
              error_eff_up   =   error_eff_up + new_error_eff_up * luminosity
              error_eff_lo   =   error_eff_lo + new_error_eff_lo * luminosity
                      
      eff            =   eff /  all_lumi
      error_eff_up   =   error_eff_up / all_lumi
      error_eff_lo   =   error_eff_lo / all_lumi
              
      outFile.write( '        {0:.4f}  {1:.4f}  {2:.4f}  \n'.format(float(eff), float(error_eff_up - eff), float(eff-error_eff_lo)) )
    
          
  
    os.system("cat "+outputFile)


def AverageDZEff(triggerMaps):
    Channels =["DoubleEle","DoubleMu","MuEle","EleMu"]
    DZefffile = open("DZEff2016.text","w")
    for ch in Channels:
    	totallumi = 0.0
    	totaleff = 0.0
	for i in triggerMaps.keys():
	    lumi = triggerMaps[i]["lumi"]
	    eff = triggerMaps[i]["DZEff"][ch]
	    #print " ch ", ch," eff ",eff," lumi ",lumi
	    totaleff += eff*lumi
	    totallumi +=  lumi
	totaleff = totaleff/totallumi
	#print "total eff ",totaleff
	DZefffile.write(ch+" %.3f\n"%eff)
    
	


	    



if __name__ == '__main__':
    print '''
--------------------------------------------------------------------------------------------------

             \                                            \ \     /       |                    
            _ \ \ \   /  _ \   __|  _` |   _` |   _ \      \ \   /  _` |  |  |   |   _ \   __| 
           ___ \ \ \ /   __/  |    (   |  (   |   __/       \ \ /  (   |  |  |   |   __/ \__ \ 
         _/    _\ \_/  \___| _|   \__,_| \__, | \___|        \_/  \__,_| _| \__,_| \___| ____/ 
                                         |___/                                                 
 
--------------------------------------------------------------------------------------------------
'''    

    AverageDZEff( Trigger['Full2016'])
    legs = ["DoubleEleLegHigPt","DoubleEleLegLowPt","DoubleMuLegHigPt","DoubleMuLegLowPt","MuEleLegHigPt", "MuEleLegLowPt","EleMuLegHigPt","EleMuLegLowPt"]
    #for leg in legs:
	#AverageLeg(leg, Trigger['Full2016'])
