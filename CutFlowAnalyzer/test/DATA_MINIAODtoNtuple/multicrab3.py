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
config.Data.splitting       = 'LumiBased'
config.Data.unitsPerJob     = 10
config.Data.lumiMask        = 'Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt'
config.Data.outLFNDirBase   = '/store/user/lpernie/'
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
datasets  = []; 
datasets.append("/DoubleMuon/Run2016B-23Sep2016-v3/MINIAOD")
datasets.append("/DoubleMuon/Run2016C-23Sep2016-v1/MINIAOD")
datasets.append("/DoubleMuon/Run2016D-23Sep2016-v1/MINIAOD")
datasets.append("/DoubleMuon/Run2016E-23Sep2016-v1/MINIAOD")
datasets.append("/DoubleMuon/Run2016F-23Sep2016-v1/MINIAOD")
datasets.append("/DoubleMuon/Run2016G-23Sep2016-v1/MINIAOD")

check_f   = open("check_crab.sh",'w')
check_f.write("#!/bin/bash\n")

plotter_f = open("for_plotter.sh",'w')
if not OnlySubmitCRAB: plotter_f.write('DATA_ch = ROOT.TChain(tree_name)\n')

if __name__ == '__main__':
  from CRABAPI.RawCommand import crabCommand
  i=0
  args = []; args.append("-99")
  for dataset in datasets:
    args[0] = 0
    # To check the status
    check_f.write("crab status -d crab_projects/crab_Hhh_"+dataset.split('/')[2]+"\n")
    # To plot easily the datasets
    if not OnlySubmitCRAB:
      sampleN   = "Hhh_"+dataset.split('/')[2]
      path      = "/fdata/hepx/store/user/lpernie/" + dataset.split('/')[1] + "/crab_Hhh_" + dataset.split('/')[2] + "/"
      NewestDir = findNewestDir(path)
      path      = path + NewestDir
      plotter_f.write('os.system("find ' + path + ' | grep root | grep -v failed > HADD/DATA_"' + sampleN + '.txt")\n')
    if OnlySubmitCRAB:
      config.Data.inputDataset = dataset
      config.General.requestName = "Hhh_"+dataset.split('/')[2]
      #config.JobType.pyCfgParams = args
      crabCommand('submit', config = config)
    i=i+1
print "bash check_crab.sh > RESULTS.txt"
if not OnlySubmitCRAB:
  plotter_f.write('os.system("cat HADD/DATA_* > HADD/DATA.txt")\n')
  plotter_f.write('with open("HADD/HADD/DATA.txt","r") as f:\n')
  plotter_f.write('  for line in f:\n')
  plotter_f.write('    if not line.isspace():\n')
  plotter_f.write('      DATA_ch.Add(str(line[:-1]))\n')
  plotter_f.write('print "DATA has", DATA_ch.GetEntries(), "entries."\n')