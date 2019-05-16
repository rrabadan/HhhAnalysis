#!/usr/bin/python                                                                                                                                        
import os
import sys

datasets  = []; NumSample = []; sampleN_short = []
Nanodatasets = []; localdirs = {}
MCxsections = []        

DateGlobalTag= 'Nano14Dec2018'
MCNanoAODGlobalTag= 'RunIIFall17NanoAODv4'

### DoubleEGRun2017B
Nanodatasets.append('/DoubleEG/Run2017B-Nano14Dec2018-v1/NANOAOD')
NumSample.append('-1'); sampleN_short.append('DoubleEGRun2017B'); MCxsections.append(-1.000000)

### DoubleEGRun2017C
Nanodatasets.append('/DoubleEG/Run2017C-Nano14Dec2018-v1/NANOAOD')
NumSample.append('-3'); sampleN_short.append('DoubleEGRun2017C'); MCxsections.append(-1.000000)

### DoubleEGRun2017D
Nanodatasets.append('/DoubleEG/Run2017D-Nano14Dec2018-v1/NANOAOD')
NumSample.append('-4'); sampleN_short.append('DoubleEGRun2017D'); MCxsections.append(-1.000000)

### DoubleEGRun2017E
Nanodatasets.append('/DoubleEG/Run2017E-Nano14Dec2018-v1/NANOAOD')
NumSample.append('-5'); sampleN_short.append('DoubleEGRun2017E'); MCxsections.append(-1.000000)

### DoubleEGRun2017F
Nanodatasets.append('/DoubleEG/Run2017F-Nano14Dec2018-v1/NANOAOD')
NumSample.append('-6'); sampleN_short.append('DoubleEGRun2017F'); MCxsections.append(-1.000000)

### DoubleMuonRun2017B
Nanodatasets.append('/DoubleMuon/Run2017B-Nano14Dec2018-v1/NANOAOD')
NumSample.append('-10'); sampleN_short.append('DoubleMuonRun2017B'); MCxsections.append(-1.000000)

### DoubleMuonRun2017C
Nanodatasets.append('/DoubleMuon/Run2017C-Nano14Dec2018-v1/NANOAOD')
NumSample.append('-12'); sampleN_short.append('DoubleMuonRun2017C'); MCxsections.append(-1.000000)

### DoubleMuonRun2017D
Nanodatasets.append('/DoubleMuon/Run2017D-Nano14Dec2018-v1/NANOAOD')
NumSample.append('-13'); sampleN_short.append('DoubleMuonRun2017D'); MCxsections.append(-1.000000)

### DoubleMuonRun2017E
Nanodatasets.append('/DoubleMuon/Run2017E-Nano14Dec2018-v1/NANOAOD')
NumSample.append('-14'); sampleN_short.append('DoubleMuonRun2017E'); MCxsections.append(-1.000000)

### DoubleMuonRun2017F
Nanodatasets.append('/DoubleMuon/Run2017F-Nano14Dec2018-v1/NANOAOD')
NumSample.append('-15'); sampleN_short.append('DoubleMuonRun2017F'); MCxsections.append(-1.000000)

### MuonEGRun2017B
Nanodatasets.append('/MuonEG/Run2017B-Nano14Dec2018-v1/NANOAOD')
NumSample.append('-19'); sampleN_short.append('MuonEGRun2017B'); MCxsections.append(-1.000000)

### MuonEGRun2017C
Nanodatasets.append('/MuonEG/Run2017C-Nano14Dec2018-v1/NANOAOD')
NumSample.append('-21'); sampleN_short.append('MuonEGRun2017C'); MCxsections.append(-1.000000)

### MuonEGRun2017D
Nanodatasets.append('/MuonEG/Run2017D-Nano14Dec2018-v1/NANOAOD')
NumSample.append('-22'); sampleN_short.append('MuonEGRun2017D'); MCxsections.append(-1.000000)

### MuonEGRun2017E
Nanodatasets.append('/MuonEG/Run2017E-Nano14Dec2018-v1/NANOAOD')
NumSample.append('-23'); sampleN_short.append('MuonEGRun2017E'); MCxsections.append(-1.000000)

### MuonEGRun2017F
Nanodatasets.append('/MuonEG/Run2017F-Nano14Dec2018-v1/NANOAOD')
NumSample.append('-24'); sampleN_short.append('MuonEGRun2017F'); MCxsections.append(-1.000000)

### RadionM270
Nanodatasets.append('/GluGluToRadionToHHTo2B2VTo2L2Nu_M-270_narrow_13TeV-madgraph_correctedcfg/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('0'); sampleN_short.append('RadionM270'); MCxsections.append(5.000000)

### RadionM280
#Nanodatasets.append('/GluGluToRadionToHHTo2B2VTo2L2Nu_M-280_narrow_13TeV-madgraph_correctedcfg/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
#NumSample.append('1'); sampleN_short.append('RadionM280'); MCxsections.append(5.000000)

### RadionM300

### RadionM320
Nanodatasets.append('/GluGluToRadionToHHTo2B2VTo2L2Nu_M-320_narrow_13TeV-madgraph_correctedcfg/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('3'); sampleN_short.append('RadionM320'); MCxsections.append(5.000000)

### RadionM350
Nanodatasets.append('/GluGluToRadionToHHTo2B2VTo2L2Nu_M-350_narrow_13TeV-madgraph_correctedcfg/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('4'); sampleN_short.append('RadionM350'); MCxsections.append(5.000000)

### RadionM400
Nanodatasets.append('/GluGluToRadionToHHTo2B2VTo2L2Nu_M-400_narrow_13TeV-madgraph_correctedcfg/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('5'); sampleN_short.append('RadionM400'); MCxsections.append(5.000000)

### RadionM450
Nanodatasets.append('/GluGluToRadionToHHTo2B2VTo2L2Nu_M-450_narrow_13TeV-madgraph_correctedcfg/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('6'); sampleN_short.append('RadionM450'); MCxsections.append(5.000000)

### RadionM500
Nanodatasets.append('/GluGluToRadionToHHTo2B2VTo2L2Nu_M-500_narrow_13TeV-madgraph_correctedcfg/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('7'); sampleN_short.append('RadionM500'); MCxsections.append(5.000000)

### RadionM600
Nanodatasets.append('/GluGluToRadionToHHTo2B2VTo2L2Nu_M-600_narrow_13TeV-madgraph_correctedcfg/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('8'); sampleN_short.append('RadionM600'); MCxsections.append(5.000000)

### RadionM650
Nanodatasets.append('/GluGluToRadionToHHTo2B2VTo2L2Nu_M-650_narrow_13TeV-madgraph_correctedcfg/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('9'); sampleN_short.append('RadionM650'); MCxsections.append(5.000000)

### RadionM750
#Nanodatasets.append('/GluGluToRadionToHHTo2B2VTo2L2Nu_M-750_narrow_13TeV-madgraph_correctedcfg/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
#NumSample.append('10'); sampleN_short.append('RadionM750'); MCxsections.append(5.000000)

### RadionM800
#Nanodatasets.append('/GluGluToRadionToHHTo2B2VTo2L2Nu_M-800_narrow_13TeV-madgraph_correctedcfg/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
#NumSample.append('11'); sampleN_short.append('RadionM800'); MCxsections.append(5.000000)

### RadionM900
Nanodatasets.append('/GluGluToRadionToHHTo2B2VTo2L2Nu_M-900_narrow_13TeV-madgraph_correctedcfg/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('12'); sampleN_short.append('RadionM900'); MCxsections.append(5.000000)

### TT
Nanodatasets.append('/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_new_pmx_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('13'); sampleN_short.append('TT'); MCxsections.append(87.310000)

### DY
Nanodatasets.append('/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('14'); sampleN_short.append('DY'); MCxsections.append(18610.000000)
Nanodatasets.append('/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('15'); sampleN_short.append('DY'); MCxsections.append(4758.900000)
Nanodatasets.append('/DYJetsToLL_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('16'); sampleN_short.append('DY'); MCxsections.append(929.100000)
Nanodatasets.append('/DYJetsToLL_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('17'); sampleN_short.append('DY'); MCxsections.append(337.100000)

### VV
Nanodatasets.append('/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('18'); sampleN_short.append('VV'); MCxsections.append(3.220000)
Nanodatasets.append('/ZZTo2L2Nu_13TeV_powheg_pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('19'); sampleN_short.append('VV'); MCxsections.append(0.564000)
#Nanodatasets.append('/ZZTo4L_13TeV_powheg_pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
#NumSample.append('20'); sampleN_short.append('VV'); MCxsections.append(1.256000)
#Nanodatasets.append('/ZZTo4L_13TeV_powheg_pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6_ext1-v1/NANOAODSIM')
#NumSample.append('20'); sampleN_short.append('VV'); MCxsections.append(1.256000)
Nanodatasets.append('/ZZTo4L_13TeV_powheg_pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_new_pmx_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('20'); sampleN_short.append('VV'); MCxsections.append(1.256000)
#Nanodatasets.append('/WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6_ext1-v1/NANOAODSIM')
#NumSample.append('21'); sampleN_short.append('VV'); MCxsections.append(49.997000)
Nanodatasets.append('/WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_new_pmx_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('21'); sampleN_short.append('VV'); MCxsections.append(49.997000)
Nanodatasets.append('/WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('22'); sampleN_short.append('VV'); MCxsections.append(12.178000)
Nanodatasets.append('/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('23'); sampleN_short.append('VV'); MCxsections.append(5.595000)
Nanodatasets.append('/WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8_v2/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('24'); sampleN_short.append('VV'); MCxsections.append(3.033000)
Nanodatasets.append('/WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('25'); sampleN_short.append('VV'); MCxsections.append(10.710000)
Nanodatasets.append('/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_new_pmx_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('26'); sampleN_short.append('VV'); MCxsections.append(4.429650)

### sT
Nanodatasets.append('/ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_new_pmx_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('27'); sampleN_short.append('sT'); MCxsections.append(136.020000)
Nanodatasets.append('/ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('28'); sampleN_short.append('sT'); MCxsections.append(80.950000)
Nanodatasets.append('/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('29'); sampleN_short.append('sT'); MCxsections.append(3.360000)
Nanodatasets.append('/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6_ext1-v1/NANOAODSIM')
NumSample.append('30'); sampleN_short.append('sT'); MCxsections.append(19.554500)
#Nanodatasets.append('/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_new_pmx_102X_mc2017_realistic_v6-v1/NANOAODSIM')
#NumSample.append('30'); sampleN_short.append('sT'); MCxsections.append(19.554500)
Nanodatasets.append('/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_new_pmx_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('31'); sampleN_short.append('sT'); MCxsections.append(19.554500)

### Wjet
Nanodatasets.append('/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('32'); sampleN_short.append('Wjet'); MCxsections.append(61526.700000)
#Nanodatasets.append('/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6_ext1-v1/NANOAODSIM')
#NumSample.append('32'); sampleN_short.append('Wjet'); MCxsections.append(61526.700000)
Nanodatasets.append('/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('33'); sampleN_short.append('Wjet'); MCxsections.append(1627.450000)
Nanodatasets.append('/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('34'); sampleN_short.append('Wjet'); MCxsections.append(435.237000)
Nanodatasets.append('/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('35'); sampleN_short.append('Wjet'); MCxsections.append(59.181000)
Nanodatasets.append('/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('36'); sampleN_short.append('Wjet'); MCxsections.append(14.580000)
Nanodatasets.append('/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('37'); sampleN_short.append('Wjet'); MCxsections.append(6.656000)
Nanodatasets.append('/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('38'); sampleN_short.append('Wjet'); MCxsections.append(1.608000)
Nanodatasets.append('/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('39'); sampleN_short.append('Wjet'); MCxsections.append(0.038900)

### ttV
Nanodatasets.append('/TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('40'); sampleN_short.append('ttV'); MCxsections.append(0.406200)
Nanodatasets.append('/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('41'); sampleN_short.append('ttV'); MCxsections.append(0.204300)
Nanodatasets.append('/TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('42'); sampleN_short.append('ttV'); MCxsections.append(0.529700)
Nanodatasets.append('/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM')
NumSample.append('43'); sampleN_short.append('ttV'); MCxsections.append(0.252900)
