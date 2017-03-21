from CRABClient.UserUtilities import config
config = config()

config.section_("General")
config.General.requestName = 'DYJetsToLL_M-10to50_TuneCUETP8M1_amcatnloFXFX_pythia8_01'
config.General.workArea = 'crab_projects'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runDiHiggsWWAnalyzer_forDY.py'
config.JobType.outputFiles = ['out_pat.root','out_ana.root']

config.section_("Data")
config.Data.inputDataset = '/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase =  '/store/user/lpernie/'
config.Data.publication = True
config.Data.outputDatasetTag = 'DYJetsToLL_M-10to50_TuneCUETP8M1_amcatnloFXFX_pythia8_01'

config.section_("Site")
config.Site.storageSite = 'T3_US_TAMU'

#NJOBS = 1
#config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
