from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'NanoPost_March18_TTbar'
config.General.workArea = 'crab_projects'
config.General.transferLogs=True
config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'crab_script.sh'
config.JobType.inputFiles = ['../keep_and_drop.txt','../postproc.py','../POGRecipesRun2.py','../HHbbWWProducer.py','../../../../PhysicsTools/NanoAODTools/scripts/haddnano.py'] #hadd nano will not be needed once nano tools are in cmssw
config.JobType.sendPythonFolder	 = True
config.section_("Data")
#config.Data.inputDataset = '/ZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8/arizzi-NanoCrabXmasRunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asym___v6_ext1-v1__315-f64d1fc6d0aff52acf7debc448857e96/USER'
#config.Data.inputDataset = '/ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8/arizzi-NanoCrabXmasRunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asym___heIV_v6-v1__307-f64d1fc6d0aff52acf7debc448857e96/USER'
#config.Data.inputDataset = '/WH_HToBB_WToLNu_M125_13TeV_amcatnloFXFX_madspin_pythia8/arizzi-NanoCrabXmasRunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asym___heIV_v6-v1__215-f64d1fc6d0aff52acf7debc448857e96/USER'
config.Data.inputDataset = '/TT_TuneCUETP8M2T4up_13TeV-powheg-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM'
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

