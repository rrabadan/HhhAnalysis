from CRABClient.UserUtilities import config
config = config()

config.section_("General")
config.General.workArea = 'crab_projects'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
#config.JobType.psetName = 'runDiHiggsWWAnalyzer_forDY.py'
config.JobType.outputFiles = ['out_ana.root']

config.section_("Data")
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase =  '/store/user/lpernie/'
config.Data.publication = False

config.section_("Site")
config.Site.storageSite = 'T3_US_TAMU'


psetName=['runDiHiggsWWAnalyzer_forDYJetsToLL.py','runDiHiggsWWAnalyzer_forDY0J.py','runDiHiggsWWAnalyzer_forDY1J.py','runDiHiggsWWAnalyzer_forDY2J.py']
i=0
if __name__ == '__main__':
    from CRABAPI.RawCommand import crabCommand
    for dataset in ['/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM','/DYToLL_0J_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM','/DYToLL_1J_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM','/DYToLL_2J_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM']:
        config.JobType.psetName = psetName[i]
        config.Data.inputDataset = dataset
        config.General.requestName = dataset.split('/')[1]
        crabCommand('submit', config = config)
        i=i+1
