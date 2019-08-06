from CRABClient.UserUtilities import config
config = config()

config.section_("General")
config.General.workArea     = 'crab_projects'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName   = 'Analysis'
config.JobType.outputFiles  = ['out_ana.root']
config.JobType.psetName     = "runDiHiggsWWAnalyzer.py"
config.JobType.inputFiles   = ["Files/EfficienciesAndSF_BCDEFGH_Tracking.root","Files/EfficienciesAndSF_BCDEF_Tracking.root","Files/EfficienciesAndSF_GH_ISO.root","Files/EfficienciesAndSF_BCDEF_ID.root","Files/EfficienciesAndSF_BCDEF_trigger.root","Files/EfficienciesAndSF_GH_Tracking.root","Files/EfficienciesAndSF_BCDEF_ISO.root","Files/EfficienciesAndSF_GH_ID.root","Files/EfficienciesAndSF_GH_trigger.root"]

config.section_("Data")
config.Data.inputDBS        = 'global'
#config.Data.splitting       = 'FileBased'
#config.Data.unitsPerJob     = 5
config.Data.splitting       = 'LumiBased'
config.Data.unitsPerJob     = 20
config.Data.lumiMask        = 'Cert_294927-304120_13TeV_PromptReco_Collisions17_JSON.txt'
config.Data.outLFNDirBase   = '/store/user/tahuang/Collisioins17'
config.Data.publication     = False

config.section_("Site")
config.Site.storageSite     = 'T3_US_TAMU'
OnlySubmitCRAB = True

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

datasets  = []; 
#datasets.append("/DoubleMuon/Run2016B-23Sep2016-v3/MINIAOD")
#datasets.append("/DoubleMuon/Run2016C-23Sep2016-v1/MINIAOD")
#datasets.append("/DoubleMuon/Run2016D-23Sep2016-v1/MINIAOD")
#datasets.append("/DoubleMuon/Run2016E-23Sep2016-v1/MINIAOD")
#datasets.append("/DoubleMuon/Run2016F-23Sep2016-v1/MINIAOD")
#datasets.append("/DoubleMuon/Run2016G-23Sep2016-v1/MINIAOD")
#datasets.append("/DoubleMuon/Run2016H-03Feb2017_ver2-v1/MINIAOD")
#datasets.append("/DoubleMuon/Run2016H-03Feb2017_ver3-v1/MINIAOD")
#datasets.append("/DoubleMuon/Run2017A-PromptReco-v1/MINIAOD")
#datasets.append("/DoubleMuon/Run2017A-PromptReco-v2/MINIAOD")
#datasets.append("/DoubleMuon/Run2017A-PromptReco-v3/MINIAOD")
#datasets.append("/DoubleMuon/Run2017B-PromptReco-v1/MINIAOD")
#datasets.append("/DoubleMuon/Run2017B-PromptReco-v2/MINIAOD")
#datasets.append("/DoubleMuon/Run2017C-PromptReco-v1/MINIAOD")
#datasets.append("/DoubleMuon/Run2017C-PromptReco-v2/MINIAOD")
#datasets.append("/DoubleMuon/Run2017C-PromptReco-v3/MINIAOD")
#datasets.append("/DoubleMuon/Run2017D-PromptReco-v1/MINIAOD")
#datasets.append("/DoubleMuon/Run2017E-PromptReco-v1/MINIAOD")
datasets.append("/DoubleMuon/Run2017B-12Sep2017-v1/MINIAOD")
datasets.append("/DoubleMuon/Run2017C-12Sep2017-v1/MINIAOD")


check_f = open("check_crab.sh",'w'); check_f.write("#!/bin/bash\n")
resub_f = open("resub_crab.sh",'w'); resub_f.write("#!/bin/bash\n")

plotter_f = open("for_plotter.py",'w')
plotter_f.write('if( whichSample == "Data" ):\n')
if __name__ == '__main__':
  from CRABAPI.RawCommand import crabCommand
  i=0
  args = []; args.append("-99")
  for dataset in datasets:
    args[0] = 0
    # To check the status
    check_f.write("crab status -d crab_projects/crab_Hhh_"+dataset.split('/')[2]+"\n")
    resub_f.write("crab resubmit -d crab_projects/crab_Hhh_"+dataset.split('/')[2]+"\n")
    # To plot easily the datasets
    if not OnlySubmitCRAB:
      sampleN   = "Hhh_"+dataset.split('/')[2]
      path      = "/fdata/hepx/store/user/tahuang/" + dataset.split('/')[1] + "/crab_Hhh_" + dataset.split('/')[2] + "/"
      NewestDir = findNewestDir(path)
      path      = path + NewestDir
      plotter_f.write('  Find_str.append("find ' + path + ' | grep root | grep -v failed > HADD/DATA_' + sampleN + '.txt")\n')
    if OnlySubmitCRAB:
      config.Data.inputDataset = dataset
      config.General.requestName = "Hhh_"+dataset.split('/')[2]
      #config.JobType.pyCfgParams = args
      crabCommand('submit', config = config)
    i=i+1
print "bash check_crab.sh > RESULTS.txt"
if not OnlySubmitCRAB:
        plotter_f.write('  this_cat      = "cat HADD/DATA_* > HADD/DATA.txt"\n')
