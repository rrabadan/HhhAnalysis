from CRABClient.UserUtilities import config
config = config()

config.section_("General")
config.General.workArea = 'crab_projects'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'patTuple_cutana_mujets_80X_cfg.py'
config.JobType.outputFiles = ['out_ana.root']

config.section_("Data")
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 10
config.Data.lumiMask = 'Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt'
#config.Data.outLFNDirBase = '/store/group/alca_muonalign/lpernie/DoubleMuonRun2015D_PromptReco_AOD246908-258750_7412p4_patch1_patTuple_cutana_mujets'
#config.Data.outLFNDirBase = '/store/group/lpcdarksusy/'
config.Data.outLFNDirBase = '/store/user/lpernie/'
config.Data.publication = False
#config.Data.outputDatasetTag = 'DoubleMuon_2016BC_AOD_808_patana_01'

config.section_("Site")
config.Site.storageSite = 'T3_US_TAMU' #'T2_CH_CERN', 'T3_US_FNALLPC'

if __name__ == '__main__':
    from CRABAPI.RawCommand import crabCommand
    for dataset in ['/DoubleMuon/Run2016B-23Sep2016-v1/MINIAOD','/DoubleMuon/Run2016B-23Sep2016-v3/MINIAOD','/DoubleMuon/Run2016C-23Sep2016-v1/MINIAOD']:
        config.Data.inputDataset = dataset
        config.General.requestName = dataset.split('/')[2]
        crabCommand('submit', config = config)
