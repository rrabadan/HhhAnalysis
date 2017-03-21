from WMCore.Configuration import Configuration

config = Configuration()

config.section_("General")

config.General.workArea = 'crab_DiHiggs_GEN_SIM'
config.General.requestName='crab_DiHigs_B3_addherwig'
config.General.transferOutputs = True

config.section_("JobType")

config.JobType.pluginName = 'PrivateMC'
#config.JobType.psetName = 'DiHiggs_Run2_cfi_GEN_SIM.py'
config.JobType.psetName = 'DiHiggs_Run2_cfi_GEN_SIM_addHerwig_test.py'
#config.JobType.inputFiles = ['/fdata/hepx/store/user/taohuang/Pheno/HH-bbWW-B6-20160518-leptonW-1000000.hepmc']
#config.JobType.inputFiles = ['HH-bbWW-B3-13TeV-leptonW-OnlyME-10k.hepmc']
config.JobType.inputFiles = ['MEHiggsPair.so']


config.section_("Data")

#config.Data.inputDBS = 'global'
config.Data.outputPrimaryDataset = 'xSM_HeavyHiggs2DiHiggs2bbWW_B3_leptonW_CMSSW80X_13TeV_1M'
#config.Data.splitting = 'FileBased'
config.Data.splitting = 'EventBased'

config.Data.unitsPerJob = 100
config.Data.totalUnits = 10

config.Data.outLFNDirBase = '/store/user/tahuang/'

config.Data.publication = True
#config.Data.inputDataset='/MinBias_TuneCUETP8M1_13TeV-pythia8/RunIIWinter15GS-MCRUN2_71_V1-v1/GEN-SIM'
#config.Data.inputDataset = '/MinBias_TuneCUETP8M1_14TeV-pythia8/PhaseIIFall16GS82-90X_upgrade2023_realistic_v1-v1/GEN-SIM'
#config.Data.userInputFiles = ['/fdata/hepx/store/user/taohuang/Pheno/HH-bbWW-B6-20160518-leptonW-1000000.hepmc']

config.Data.outputDatasetTag= config.Data.outputPrimaryDataset

config.section_("Site")

config.Site.storageSite = 'T3_US_TAMU'

config.Site.whitelist = ['T3_US_TAMU']


