from CRABClient.UserUtilities import config
config = config()

config.section_("General")
config.General.requestName = 'TTTo2L2Nu_powheg_miniAODv2_v0_ext1_01'
config.General.workArea = 'crab_projects'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runDiHiggsWWAnalyzer_forTT.py'
config.JobType.outputFiles = ['out_ana.root']

config.section_("Data")
config.Data.inputDataset = '/TTTo2L2Nu_13TeV-powheg/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase =  '/store/user/lpernie/'
config.Data.publication = True
config.Data.outputDatasetTag = 'TTTo2L2Nu_powheg_miniAODv2_v0_ext1_01'

config.section_("Site")
config.Site.storageSite = 'T3_US_TAMU'

#NJOBS = 1
#config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
