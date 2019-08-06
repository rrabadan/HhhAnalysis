from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'NanoPost_March18_TTbar_v2'
config.General.workArea = 'crab_projects'
config.General.transferLogs=True
config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'crab_script.sh'
config.JobType.inputFiles = ['../keep_and_drop.txt','../postproc.py','../POGRecipesRun2.py','../leptonSF.zip','../HHbbWWProducer.py','../../../../PhysicsTools/NanoAODTools/scripts/haddnano.py'] #hadd nano will not be needed once nano tools are in cmssw
config.JobType.sendPythonFolder	 = True
config.section_("Data")
config.Data.inputDataset = '/GluGluToRadionToHHTo2B2VTo2L2Nu_M-260_narrow_13TeV-madgraph-v2/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 1
#config.Data.unitsPerJob = 100
#config.Data.totalUnits = 2000
#config.Data.inputDBS='phys03'
username = "tahuang"
config.Data.outLFNDirBase = '/store/user/%s/%s/'%(username, config.General.requestName)
config.Data.publication = False
config.Data.outputDatasetTag = 'NanoTestPost'
config.section_("Site")
config.Site.storageSite = "T3_US_TAMU"

