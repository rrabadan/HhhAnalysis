#!/usr/bin/python
import os
import sys

fdatadir = "/fdata/hepx/store/user/taohuang/NANOAOD/"
Nanodatasets  = []; NumSample = []; sampleN_short = []
localdirs = {}
MCxsections = []
#doTT=True; doDY=True; doVV=True; doSingleT=True; doWjets=True; dottV=True

##DoubleEG
Nanodatasets.append('/DoubleEG/Run2017B-31Mar2018-v1/NANOAOD')
NumSample.append('-1'); sampleN_short.append('DoubleEGRun2017B')
MCxsections.append(-1.0)
Nanodatasets.append('/DoubleEG/Run2017C-31Mar2018-v1/NANOAOD')
NumSample.append('-3'); sampleN_short.append('DoubleEGRun2017C')
MCxsections.append(-1.0)
Nanodatasets.append('/DoubleEG/Run2017D-31Mar2018-v1/NANOAOD')
NumSample.append('-4'); sampleN_short.append('DoubleEGRun2017D')
MCxsections.append(-1.0)
Nanodatasets.append('/DoubleEG/Run2017E-31Mar2018-v1/NANOAOD')
NumSample.append('-5'); sampleN_short.append('DoubleEGRun2017E')
MCxsections.append(-1.0)
Nanodatasets.append('/DoubleEG/Run2017F-31Mar2018-v1/NANOAOD')
NumSample.append('-6'); sampleN_short.append('DoubleEGRun2017F')
MCxsections.append(-1.0)
##DoubleMuon
Nanodatasets.append('/DoubleMuon/Run2017B-31Mar2018-v1/NANOAOD')
NumSample.append('-10'); sampleN_short.append('DoubleMuonRun2017B')
MCxsections.append(-1.0)
Nanodatasets.append('/DoubleMuon/Run2017C-31Mar2018-v1/NANOAOD')
NumSample.append('-12'); sampleN_short.append('DoubleMuonRun2017C')
MCxsections.append(-1.0)
Nanodatasets.append('/DoubleMuon/Run2017D-31Mar2018-v1/NANOAOD')
NumSample.append('-13'); sampleN_short.append('DoubleMuonRun2017D')
MCxsections.append(-1.0)
Nanodatasets.append('/DoubleMuon/Run2017E-31Mar2018-v1/NANOAOD')
NumSample.append('-14'); sampleN_short.append('DoubleMuonRun2017E')
MCxsections.append(-1.0)
Nanodatasets.append('/DoubleMuon/Run2017F-31Mar2018-v1/NANOAOD')
NumSample.append('-15'); sampleN_short.append('DoubleMuonRun2017F')
MCxsections.append(-1.0)
#MuonEG
Nanodatasets.append('/MuonEG/Run2017B-31Mar2018-v1/NANOAOD')
NumSample.append('-19'); sampleN_short.append('MuonEGRun2017B')
MCxsections.append(-1.0)
Nanodatasets.append('/MuonEG/Run2017C-31Mar2018-v1/NANOAOD')
NumSample.append('-21'); sampleN_short.append('MuonEGRun2017C')
MCxsections.append(-1.0)
Nanodatasets.append('/MuonEG/Run2017D-31Mar2018-v1/NANOAOD')
NumSample.append('-22'); sampleN_short.append('MuonEGRun2017D')
MCxsections.append(-1.0)
Nanodatasets.append('/MuonEG/Run2017E-31Mar2018-v1/NANOAOD')
NumSample.append('-23'); sampleN_short.append('MuonEGRun2017E')
MCxsections.append(-1.0)
Nanodatasets.append('/MuonEG/Run2017F-31Mar2018-v1/NANOAOD')
NumSample.append('-24'); sampleN_short.append('MuonEGRun2017F')
MCxsections.append(-1.0)


masspoints = [270, 280, 300, 320, 350, 400, 450, 500,  600, 650, 750, 800, 900]
for mass in masspoints:
    Nanodatasets.append("/GluGluToRadionToHHTo2B2VTo2L2Nu_M-%d_narrow_13TeV-madgraph_correctedcfg/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM"%mass)
    NumSample.append(masspoints.index(mass)); sampleN_short.append('RadionM%d'%mass)
    MCxsections.append(5.0)#by default, assume the cross section for signal is 5pb
#Nanodatasets.append("/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-*_narrow_13TeV-madgraph-v2/RunIISummer16NanoAOD-PUMoriond17_31Mar2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM")
#NumSample.append('2'); sampleN_short.append('Graviton')
##/GluGluToRadionToHHTo2B2VTo2L2Nu_M-1000_narrow_13TeV-madgraph-v2/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2017_TrancheIV_v6-v1/MINIAODSIM

##/GluGluToRadionToHHTo2B2VTo2L2Nu_M-1000_narrow_13TeV-madgraph-v2/RunIISummer16NanoAOD-PUMoriond17_31Mar2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM

# TT#
Nanodatasets.append('/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/NANOAODSIM')
NumSample.append('13'); sampleN_short.append('TT')
MCxsections.append(87.31)
# DY
Nanodatasets.append('/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM')
NumSample.append('14'); sampleN_short.append('DY')
MCxsections.append(18610.0)
Nanodatasets.append('/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM')
NumSample.append('15'); sampleN_short.append('DY')
MCxsections.append(4758.9)
Nanodatasets.append('/DYJetsToLL_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM')
NumSample.append('16'); sampleN_short.append('DY')
MCxsections.append(929.1)
Nanodatasets.append('/DYJetsToLL_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM')
NumSample.append('17'); sampleN_short.append('DY')
MCxsections.append(337.1)

# VV

#Nanodatasets.append('/WZTo3LNu_3Jets_MLL-4to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM')
#Nanodatasets.append('/WZTo3LNu_2Jets_MLL-4to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM')
#Nanodatasets.append('/WZTo3LNu_1Jets_MLL-4to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM')
#Nanodatasets.append('/WZTo3LNu_0Jets_MLL-4to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM')


Nanodatasets.append('/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM')
NumSample.append('18'); sampleN_short.append('VV')
MCxsections.append(3.22)
Nanodatasets.append('/ZZTo2L2Nu_13TeV_powheg_pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM')
NumSample.append('19'); sampleN_short.append('VV')
MCxsections.append(0.564)
Nanodatasets.append('/ZZTo4L_13TeV_powheg_pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM')
NumSample.append('20'); sampleN_short.append('VV')
MCxsections.append(1.256)
#/WWToLNuQQ_NNPDF31_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/NANOAODSIM 
#/WWTo2L2Nu_NNPDF31_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/NANOAODSIM 
Nanodatasets.append('/WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/NANOAODSIM')
NumSample.append('21'); sampleN_short.append('VV')
MCxsections.append(49.997)#
Nanodatasets.append('/WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM')
NumSample.append('22'); sampleN_short.append('VV') ### not available now
MCxsections.append(12.178)
Nanodatasets.append('/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM')
NumSample.append('23'); sampleN_short.append('VV')
MCxsections.append(5.595)
Nanodatasets.append('/WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8_v2/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM')
NumSample.append('24'); sampleN_short.append('VV') ### not available now 
MCxsections.append(3.033)
Nanodatasets.append('/WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/NANOAODSIM')
NumSample.append('25'); sampleN_short.append('VV')
MCxsections.append(10.71)
Nanodatasets.append('/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM')
NumSample.append('26'); sampleN_short.append('VV')
MCxsections.append(4.42965)
##sT
## production rate t-channel(~218pb@13TeV) > tW > s-channel 
#Nanodatasets.append('/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/NANOAODSIMi')                   
#Nanodatasets.append('/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/NANOAODSIM')
#subleading in sT
#/ST_s-channel_antitop_leptonDecays_13TeV-PSweights_powheg-pythia/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM 
#/ST_s-channel_top_leptonDecays_13TeV-PSweights_powheg-pythia/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM 
Nanodatasets.append('/ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM')
NumSample.append('27'); sampleN_short.append('sT')
MCxsections.append(136.02)
Nanodatasets.append('/ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/NANOAODSIM') 
NumSample.append('28'); sampleN_short.append('sT')
MCxsections.append(80.95)
Nanodatasets.append('/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/NANOAODSIM')
NumSample.append('29'); sampleN_short.append('sT')
MCxsections.append(3.36)
Nanodatasets.append('/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/NANOAODSIM')
NumSample.append('30'); sampleN_short.append('sT')
MCxsections.append(19.5545)
Nanodatasets.append('/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM')
NumSample.append('31'); sampleN_short.append('sT')
MCxsections.append(19.5545)


# W + Jets

#Nanodatasets.append('/WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM')


Nanodatasets.append('/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/NANOAODSIM')
NumSample.append('32'); sampleN_short.append('Wjet')
MCxsections.append(61526.7)
Nanodatasets.append('/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM')
NumSample.append('33'); sampleN_short.append('Wjet')
MCxsections.append(1627.45)

Nanodatasets.append('/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM')
NumSample.append('34'); sampleN_short.append('Wjet')
MCxsections.append(435.237)
Nanodatasets.append('/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM')
NumSample.append('35'); sampleN_short.append('Wjet')
MCxsections.append(59.181)
Nanodatasets.append('/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM')
NumSample.append('36'); sampleN_short.append('Wjet')### not available now
MCxsections.append(14.58)
Nanodatasets.append('/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM')
NumSample.append('37'); sampleN_short.append('Wjet')
MCxsections.append(6.656)
Nanodatasets.append('/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM')
NumSample.append('38'); sampleN_short.append('Wjet')
MCxsections.append(1.608)
Nanodatasets.append('/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM')
NumSample.append('39'); sampleN_short.append('Wjet')
MCxsections.append(0.0389)
# tt + V

Nanodatasets.append('/TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM')
NumSample.append('40'); sampleN_short.append('ttV')
MCxsections.append(0.4062)
Nanodatasets.append('/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM')
NumSample.append('41'); sampleN_short.append('ttV')
MCxsections.append(0.2043)
Nanodatasets.append('/TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM')
NumSample.append('42'); sampleN_short.append('ttV')
MCxsections.append(0.5297)
Nanodatasets.append('/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM')
NumSample.append('43'); sampleN_short.append('ttV')
MCxsections.append(0.2529)
#Nanodatasets.append('/TTZToLL_M-1to10_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM')

"""
alljobtypes = set(sampleN_short)
for job in alljobtypes:
    localdirs[job] = []

for ijob, job in enumerate(Nanodatasets):
    nsample = int(NumSample[ijob])
    jobtype = sampleN_short[ijob]
    dataname = ""
    datadir = " "
    #print "nsample ",nsample, " jobtype ",jobtype
    if nsample < 0:
        datadir = sampleN_short[ijob]
	dataname = job
        #print "real data nsample ",nsample, " datadir ",datadir
    elif nsample > 0:
        datadir = job.split('/')[1]
        #print "MC nsample ",nsample, " datadir ",datadir, "MiniAOD dataset ",job.split('/')
	#query = "dataset dataset=/%s/*/NANOAODSIM"%(datadir)
        #pdata = os.popen("dasgoclient -limit=0 -query='{query}'".format(query = query))	
        #founddataset = False
	#for line in pdata:
	#    #print "dataset ",line," datatype ",datadir
	#    if datadir in line:
	#        founddataset = True
	#        dataname = line[:-1]	
	#if not(founddataset): 
	#    print "WARNING!!!!! no dataset found for ",datadir
    localdirs[jobtype].append(os.path.join(fdatadir, datadir))


outAnalist = {}
outAnadir = "/fdata/hepx/store/user/taohuang/HHNtuple_20180328_fixedleptonDZeff/"
for i,datasetname in enumerate( NanoNanodatasets ):
    sampleName = sampleN_short[i]
    if NumSample[i] < 0:
    	sampleName = "Data"
    outAnafile = os.path.join(outAnadir, NanoNanodatasets[i].split('/')[1])
    if hasattr(outAnalist, sampleName):
	outAnalist[sampleName].append(outAnafile)
    else:
	outAnalist[sampleName] = []
	outAnalist[sampleName].append(outAnafile)

"""
