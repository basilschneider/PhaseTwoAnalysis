from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'TChiWZ_300_250_run05'
config.General.workArea = 'crab_tasks/'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'scripts/produceNtuples_cfg.py'
config.JobType.pyCfgParams= ['skim=False','inputFormat=PAT','outFilename=MiniEvents.root']
# Uncomment the following line when running on PAT events
#config.JobType.inputFiles = ['TMVAClassification_BDT.weights.xml']
config.JobType.outputFiles = ['MiniEvents.root']

config.section_("Data")
config.Data.inputDataset = '/SMS-TChiWZ_ZToLL_mChargino-300_mLsp-250_TuneCUETP8M1_14TeV-madgraphMLM-pythia8/PhaseIITDRSpring17MiniAOD-PU140_91X_upgrade2023_realistic_v3-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
# Uncomment to run on a fraction of the dataset
config.Data.totalUnits = 2
config.Data.outLFNDirBase = '/store/user/bschneid/upgradeFS/'
config.Data.publication = False

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
config.Site.ignoreGlobalBlacklist = True
