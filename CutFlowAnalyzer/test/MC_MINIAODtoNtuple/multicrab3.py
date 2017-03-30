from CRABClient.UserUtilities import config
config = config()

config.section_("General")
config.General.workArea     = 'crab_projects'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName   = 'Analysis'
config.JobType.outputFiles  = ['out_ana.root']
config.JobType.psetName     = "runDiHiggsWWAnalyzer.py"

config.section_("Data")
config.Data.inputDBS        = 'global'
config.Data.splitting       = 'FileBased'
config.Data.unitsPerJob     = 1
config.Data.outLFNDirBase   =  '/store/user/lpernie/'
config.Data.publication     = False

config.section_("Site")
config.Site.storageSite     = 'T3_US_TAMU'

import os
import glob
import operator
def findNewestDir(directory):
  os.chdir(directory)
  dirs = {}
  for dir in glob.glob('*'):
    if os.path.isdir(dir):
      dirs[dir] = os.path.getctime(dir)
  lister = sorted(dirs.iteritems(), key=operator.itemgetter(1))
  return lister[-1][0]

OnlySubmitCRAB=True
datasets  = []; NumSample = []; sampleN_short = []
doTT=True; doDY=True; doVV=True; doSingleT=True; doWjets=True; dottV=True
#doTT=False; doDY=False; doVV=False; doSingleT=False; doWjets=False; dottV=False
# TT
if doTT:
  datasets.append('/TTTo2L2Nu_13TeV-powheg/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/MINIAODSIM')
  NumSample.append('13'); sampleN_short.append('TT')
# DY
if doDY:
  datasets.append('/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('14'); sampleN_short.append('DY')
  datasets.append('/DYToLL_0J_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
  NumSample.append('15'); sampleN_short.append('DY')
  datasets.append('/DYToLL_1J_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
  NumSample.append('16'); sampleN_short.append('DY')
  datasets.append('/DYToLL_2J_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
  NumSample.append('17'); sampleN_short.append('DY')
# VV
if doVV:
  datasets.append('/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('18'); sampleN_short.append('VV')
  datasets.append('/ZZTo2L2Nu_13TeV_powheg_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('19'); sampleN_short.append('VV')
  datasets.append('/ZZTo4L_13TeV_powheg_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('20'); sampleN_short.append('VV')
  datasets.append('/WWToLNuQQ_aTGC_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('21'); sampleN_short.append('VV')
  datasets.append('/WWTo2L2Nu_MWW-600To800_aTGC_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('22'); sampleN_short.append('VV')
  datasets.append('/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('23'); sampleN_short.append('VV')
  datasets.append('/WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('24'); sampleN_short.append('VV')
  datasets.append('/WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v3/MINIAODSIM')
  NumSample.append('25'); sampleN_short.append('VV')
  datasets.append('/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('26'); sampleN_short.append('VV')
# Single-t
if doSingleT:
  datasets.append('/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('27'); sampleN_short.append('sT')
  datasets.append('/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('28'); sampleN_short.append('sT')
  datasets.append('/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('29'); sampleN_short.append('sT')
  datasets.append('/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
  NumSample.append('30'); sampleN_short.append('sT')
  datasets.append('/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
  NumSample.append('31'); sampleN_short.append('sT')
# W + Jets
if doWjets:
  datasets.append('/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('32'); sampleN_short.append('Wjet')
  datasets.append('/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM')
  NumSample.append('33'); sampleN_short.append('Wjet')
  datasets.append('/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM')
  NumSample.append('34'); sampleN_short.append('Wjet')
  datasets.append('/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
  NumSample.append('35'); sampleN_short.append('Wjet')
  datasets.append('/WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
  NumSample.append('36'); sampleN_short.append('Wjet')
  datasets.append('/WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
  NumSample.append('37'); sampleN_short.append('Wjet')
  datasets.append('/WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
  NumSample.append('38'); sampleN_short.append('Wjet')
  datasets.append('/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
  NumSample.append('39'); sampleN_short.append('Wjet')
# tt + V
if dottV:
  datasets.append('/TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('40'); sampleN_short.append('ttV')
  datasets.append('/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM')
  NumSample.append('41'); sampleN_short.append('ttV')
  datasets.append('/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
  NumSample.append('42'); sampleN_short.append('ttV')
  datasets.append('/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext3-v1/MINIAODSIM')
  NumSample.append('43'); sampleN_short.append('ttV')

plotter_f = open("for_plotter.py",'w')
check_f   = open("check_crab.sh",'w')
if not OnlySubmitCRAB:
  check_f.write("#!/bin/bash\n")
args = []; args.append("-99")
lastSampleShort="NOSAMPLESHORT"
lastSample="NOSAMPLE"
i=0
if __name__ == '__main__':
  from CRABAPI.RawCommand import crabCommand
  for dataset in datasets:
    # To check the status
    check_f.write("crab status -d crab_projects/crab_"+dataset.split('/')[1]+"\n")
    # To plot easily the datasets
    if not OnlySubmitCRAB:
      sampleN   = dataset.split('/')[1].split('_')[0]+dataset.split('/')[1].split('_')[1]
      path      = "/fdata/hepx/store/user/lpernie/" + dataset.split('/')[1] + "/crab_" + dataset.split('/')[1] + "/"
      NewestDir = findNewestDir(path)
      path      = path + NewestDir
      newSample = True
      if(sampleN_short[i]==lastSampleShort): newSample = False
      if((newSample and lastSampleShort!="NOSAMPLESHORT") or i==int(len(NumSample)-1) ):
        plotter_f.write('os.system("cat HADD/' + sampleN_short[i-1] + '_*' ' > HADD/' + sampleN_short[i-1] + '.txt")\n')
        plotter_f.write('os.system("hadd -T -f -k /fdata/hepx/store/user/lpernie/' + oldDataset + '/crab_' + oldDataset + '.root @HADD/' + sampleN_short[i-1] + '.txt")\n')
        plotter_f.write('N_tot_path_' + sampleN_short[i-1] + ' = "/fdata/hepx/store/user/lpernie/' + oldDataset + '/crab_' + oldDataset + '.root"\n')
        plotter_f.write(sampleN_short[i-1] + '_file =  ROOT.TFile.Open(N_tot_path_' + sampleN_short[i-1] + ',"read"); h_' + sampleN_short[i-1] + ' =  ROOT.TH1F(' + sampleN_short[i-1] + '_file.Get("DiHiggsWWBBAna/hevent")); nTOT_' + sampleN_short[i-1] + ' = h_' + sampleN_short[i-1] + '.GetBinContent(2);\n')
        plotter_f.write('with open("HADD/' + sampleN_short[i-1] +'.txt","r") as f:\n')
        plotter_f.write('  for line in f:\n')
        plotter_f.write('    if not line.isspace():\n')
        plotter_f.write('      '+sampleN_short[i-1]+'_ch.Add(str(line[:-1]))\n')
        plotter_f.write('print "'+sampleN_short[i-1]+' has", '+sampleN_short[i-1]+'_ch.GetEntries(), "entries."\n')
      if(newSample):  plotter_f.write(sampleN_short[i]+'_ch = ROOT.TChain(tree_name)\n')
      plotter_f.write('os.system("find ' + path + ' | grep root | grep -v failed > HADD/' + sampleN_short[i] + "_" + sampleN + '.txt")\n')
      lastSampleShort = sampleN_short[i]
      lastSample      = sampleN
    if OnlySubmitCRAB:
      config.Data.inputDataset = dataset
      config.General.requestName = dataset.split('/')[1]
      args[0] = NumSample[i]
      config.JobType.pyCfgParams = args
      crabCommand('submit', config = config)
    i=i+1
    oldDataset = dataset.split('/')[1]
print "bash check_crab.sh > RESULTS.txt"
