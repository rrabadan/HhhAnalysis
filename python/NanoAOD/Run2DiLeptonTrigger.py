import re
DiLeptonTrigger = {}

errorcolor1 = '\x1b[1;31m'
errorcolor2 = '\x1b[0m'
#------------------------------------  Run 2016 -------------------------------------------
"""
period  runrange
Run2016A  271036-271658
Run2016B  272007-275376
Run2016C  275657-276283
Run2016D  276315-276811
Run2016E  276831-277420
Run2016F  277772-278808
Run2016G  278820-280385
Run2016H  280919-284044
"""


DiLeptonTrigger["Run2016"] = {}
Run2016DoubleMuon          = {}
Run2016DoubleEG            = {}
Run2016MuonEG              = {}
DiLeptonTrigger["Run2016"]["DoubleMuon"] = Run2016DoubleMuon
DiLeptonTrigger["Run2016"]["DoubleEG"]   = Run2016DoubleEG
DiLeptonTrigger["Run2016"]["MuonEG"]     = Run2016MuonEG


Run2016DoubleMuon = {
    "Runrange:[273158,281612]" : {
				  "HLTPaths":["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL"],
				  "L1Pts": { ##each element contains l1pt for leg1 and lep2
					     "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL":[[11, 4]], 
					     "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL":[[11, 4]]
				            },
				  "IntLumi": 27.655
                                  },
    "Runrange:[281613,284044]" : {
				  "HLTPaths":["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ"],
				  "L1Pts": { ##each element contains l1pt for leg1 and lep2
					     "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL":[[11, 4]], 
					     "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL":[[11, 4]]
				            },
				  "IntLumi": 8.263
    				  }
}

Run2016DoubleEG = {
    "Runrange:[273158,284044]" : {
				  "HLTPaths": ["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"],
				  "L1Pts": { ##each element contains l1pt for leg1 and lep2
					     "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ":[[24, 0],[22,12],[18, 17]] ## or logic
				            },
				  "IntLumi": 35.918
                                  }
}

Run2016MuonEG   = {## double check
    "Runrange:[273158,278272]" : { ##B-F
				  "HLTPaths": ["HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL","HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL"],
				  "L1Pts": { ##each element contains l1pt for leg1 and lep2
					     "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL":[[12, 10]],
					     "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL":[[5, 15]]
				            },
				  "IntLumi": 17.641
                                  },
    "Runrange:[278273,284044]" : {##G,H
				  "HLTPaths": ["HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ","HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ"],
				  "L1Pts": { ##each element contains l1pt for leg1 and lep2
					     "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ":[[20, 10]],
					     "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ":[[5, 18]]
				            },
				  "IntLumi": 18.277
                                  },
}

#------------------------------------  Run 2017 -------------------------------------------
"""
#period, absolute runrange, collision runrange, sqrt(s), Lumi(no normtag), Lumi(normtag) 
Run2017B  297020-299329   297046-29932913 4.792 4.823
Run2017C  299337-302029   299368-30202913 9.755 9.664
Run2017D  302030-303434   302030-30343413 4.319 4.252
Run2017E  303435-304826   303824-30479713 9.424 9.278
Run2017F  304911-306462   305040-30646213 13.50 13.540
"""

DiLeptonTrigger["Run2017"] = {}

Run2017DoubleMuon = {
    "Runrange:[297046,299329]" : {# B
				  "HLTPaths":["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8"],
				  "L1Pts": { ##each element contains l1pt for leg1 and lep2
					     "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8":[[12, 5]] ## or logic
				            },
				  "IntLumi": 4.794
                                  },
    "Runrange:[299368,306462]" : {#C-F
				  "HLTPaths":["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"],
				  "L1Pts": { ##each element contains l1pt for leg1 and lep2
					     "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8":[[12, 5]], ## or logic
					     "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8":[[12, 5]] ## or logic
				            },
				  "IntLumi": 36.733
    				  }
}

Run2017DoubleEG = {
    "Runrange:[297046,306462]" : {#B-F
				  "HLTPaths": ["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"],
				  "L1Pts": { ##each element contains l1pt for leg1 and lep2
					     "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL":[[22,12],[23, 10],[24,0]] ## or logic
				            },
				  "IntLumi":41.557
                                  }
}

##question: why HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ is used while HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL is not prescaled 
Run2017MuonEG   = {
    "Runrange:[297046,299329]" : {# B
				  "HLTPaths": ["HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ","HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ","HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ"],
				  "L1Pts": { ##each element contains l1pt for leg1 and lep2
					     "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ":[[20, 17],[23, 10]], ## or logic
					     "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ":[[5, 18],[7,20]],
					     "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ":[[5,18],[7,20]]
				            },
				  "IntLumi": 4.794
                                  },
    "Runrange:[299368,306462]" : { #C-F
				  "HLTPaths": ["HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL","HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ","HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ"],
				  "L1Pts": { ##each element contains l1pt for leg1 and lep2
					     "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL":[[20, 17],[23, 10]], ## or logic
					     "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ":[[5, 18],[7,20]],
					     "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ":[[20, 17],[23, 10]],
					     "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ":[[5,18],[7,20]]
				            },
				  "IntLumi": 36.733
                                  },
}
DiLeptonTrigger["Run2017"]["DoubleMuon"] = Run2017DoubleMuon
DiLeptonTrigger["Run2017"]["DoubleEG"]   = Run2017DoubleEG
DiLeptonTrigger["Run2017"]["MuonEG"]     = Run2017MuonEG

#------------------------------------  Run 2018 -------------------------------------------
"""
#period, absolute runrange, collision runrange, sqrt(s), Lumi(no normtag), Lumi(normtag) 
Run2018A  315252-316995  315252-316995 13 13.48  14.00 
Run2018B  316998-319312  317080-319310 13 6.785  7.10 
Run2018C  319313-320393  319337-320065 13 6.612  6.94 
Run2018D  320394-325273  320673-325175 13 31.95  31.93 
"""

DiLeptonTrigger["Run2018"] = {}

Run2018DoubleMuon = {
    "Runrange:[315252,325175]" : {# A-D
				  "HLTPaths": ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8"],
				  "L1Pts": {
				            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8": [[12, 5]]
				             },
				  "IntLumi": 59.735
                                  }
}

Run2018DoubleEG   = {
    "Runrange:[315252,325175]" : {# A-D
				  "HLTPaths": ["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"],
				  "L1Pts": {
				      	    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL": [[22, 10],[30, 0]] ##or logic
				            },
				  "IntLumi": 59.735
                                  }
}

Run2018MuonEG     = {
    "Runrange:[315252,325175]" : {# A-D
				  "HLTPaths": ["HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ","HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ"],
				  "L1Pts": {
				            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL": [[23, 10],[20,17],[22, 0]],
					    "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ":[[5, 20]],
					    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ":[[5, 20]]
				            },
				  "IntLumi": 59.735
                                  }
}
DiLeptonTrigger["Run2018"]["DoubleMuon"] = Run2018DoubleMuon
DiLeptonTrigger["Run2018"]["DoubleEG"]   = Run2018DoubleEG
DiLeptonTrigger["Run2018"]["MuonEG"]     = Run2018MuonEG

IntLumiyear = {
    2016 : 35.92,
    2017 : 41.53,
    2018 : 59.74
}

DiLeptonTrigger_allHLTPaths_prescale = {}
def addHLTPathprescale(Runyear):
    allHLTPaths = []
    allHLTPaths_prescale = {}
    for triggertype in DiLeptonTrigger["Run%d"%Runyear]:
        for key in DiLeptonTrigger["Run%d"%Runyear][triggertype].keys():
	    #allHLTPaths.extend(x for x in  DiLeptonTrigger["Run%d"%Runyear][triggertype][key]["HLTPaths"] if x not in allHLTPaths)
            for path in  DiLeptonTrigger["Run%d"%Runyear][triggertype][key]["HLTPaths"]:
	        if path not in allHLTPaths_prescale.keys():
		    allHLTPaths_prescale[path] =  DiLeptonTrigger["Run%d"%Runyear][triggertype][key]["IntLumi"]
		else:
		    allHLTPaths_prescale[path] =  DiLeptonTrigger["Run%d"%Runyear][triggertype][key]["IntLumi"] + allHLTPaths_prescale[path]

    for path in allHLTPaths_prescale.keys():
	allHLTPaths_prescale[path] = round(allHLTPaths_prescale[path]/IntLumiyear[Runyear], 3)
    DiLeptonTrigger_allHLTPaths_prescale["Run%d"%Runyear] = allHLTPaths_prescale
    
##ignore 2016 results
#addHLTPathprescale(2016)
addHLTPathprescale(2017)
addHLTPathprescale(2018)
	        
#print "DiLeptonTrigger_allHLTPaths_prescale \n ",DiLeptonTrigger_allHLTPaths_prescale        


def findHLTPathsAndL1Pts(Runyear, datatype, run, isMC = False):
     
    thisdict = DiLeptonTrigger["Run%d"%Runyear][datatype]
    #print "Run%d"%Runyear, " ",datatype," table ",thisdict
    hltpaths = []
    hltpath_l1pts = {}
    for key in thisdict.keys():
	if not isMC:
	    runrange = re.findall(r"[-+]?\d*\.\d+|\d+", key)
	    #print "key in DiLeptonTrigger dict ",key," runrange ",runrange
	    if run >= int(runrange[0]) and  run <= int(runrange[1]):
		#print "HLT paths and l1pts ",thisdict[key]["L1Pts"]
		return thisdict[key]["HLTPaths"], thisdict[key]["L1Pts"]
	else:
	    for path in thisdict[key]["HLTPaths"]:
	        if path not in hltpaths:
		    hltpaths.append(path)
                    hltpath_l1pts[path] = thisdict[key]["L1Pts"][path]
    if isMC and len(hltpaths)>0:
        return hltpaths, hltpath_l1pts
    else:
	print errorcolor1+"warning!!! no HLT path is found for %d%s, run number %d"%(Runyear, datatype, run)+errorcolor2
	return None

