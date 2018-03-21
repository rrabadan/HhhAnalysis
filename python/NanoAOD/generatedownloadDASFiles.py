#!/usr/bin/python
import os
import sys

datasets  = []; NumSample = []; sampleN_short = []
#doTT=True; doDY=True; doVV=True; doSingleT=True; doWjets=True; dottV=True


##DoubleEG
datasets.append('/DoubleEG/Run2016B-05Feb2018_ver1-v1/NANOAOD')
NumSample.append('-1'); sampleN_short.append('DoubleEGRun2016Bver1')
datasets.append('/DoubleEG/Run2016B-05Feb2018_ver2-v1/NANOAOD')
NumSample.append('-2'); sampleN_short.append('DoubleEGRun2016Bver2')
datasets.append('/DoubleEG/Run2016C-05Feb2018-v1/NANOAOD')
NumSample.append('-3'); sampleN_short.append('DoubleEGRun2016C')
datasets.append('/DoubleEG/Run2016D-05Feb2018-v1/NANOAOD')
NumSample.append('-4'); sampleN_short.append('DoubleEGRun2016')
datasets.append('/DoubleEG/Run2016E-05Feb2018-v1/NANOAOD')
NumSample.append('-5'); sampleN_short.append('DoubleEGRun2016E')
datasets.append('/DoubleEG/Run2016F-05Feb2018-v1/NANOAOD')
NumSample.append('-6'); sampleN_short.append('DoubleEGRun2016F')
datasets.append('/DoubleEG/Run2016G-05Feb2018-v1/NANOAOD')
NumSample.append('-7'); sampleN_short.append('DoubleEGRun2016G')
datasets.append('/DoubleEG/Run2016H-05Feb2018_ver2-v1/NANOAOD')
NumSample.append('-8'); sampleN_short.append('DoubleEGRun2016Hver2')
datasets.append('/DoubleEG/Run2016H-05Feb2018_ver3-v1/NANOAOD')
NumSample.append('-9'); sampleN_short.append('DoubleEGRun2016Hver3')
##DoubleMuon
datasets.append('/DoubleMuon/Run2016B-05Feb2018_ver1-v1/NANOAOD')
NumSample.append('-10'); sampleN_short.append('DoubleMuonRun2016Bver1')
datasets.append('/DoubleMuon/Run2016B-05Feb2018_ver2-v1/NANOAOD')
NumSample.append('-11'); sampleN_short.append('DoubleMuonRun2016Bver2')
datasets.append('/DoubleMuon/Run2016C-05Feb2018-v1/NANOAOD')
NumSample.append('-12'); sampleN_short.append('DoubleMuonRun2016C')
datasets.append('/DoubleMuon/Run2016D-05Feb2018-v1/NANOAOD')
NumSample.append('-13'); sampleN_short.append('DoubleMuonRun2016D')
datasets.append('/DoubleMuon/Run2016E-05Feb2018-v1/NANOAOD')
NumSample.append('-14'); sampleN_short.append('DoubleMuonRun2016E')
datasets.append('/DoubleMuon/Run2016F-05Feb2018-v1/NANOAOD')
NumSample.append('-15'); sampleN_short.append('DoubleMuonRun2016F')
datasets.append('/DoubleMuon/Run2016G-05Feb2018-v1/NANOAOD')
NumSample.append('-16'); sampleN_short.append('DoubleMuonRun2016G')
datasets.append('/DoubleMuon/Run2016H-05Feb2018_ver2-v1/NANOAOD')
NumSample.append('-17'); sampleN_short.append('DoubleMuonRun2016Hver2')
datasets.append('/DoubleMuon/Run2016H-05Feb2018_ver3-v1/NANOAOD')
NumSample.append('-18'); sampleN_short.append('DoubleMuonRun2016Hver3')
#MuonEG
datasets.append('/MuonEG/Run2016B-05Feb2018_ver1-v1/NANOAOD')
NumSample.append('-19'); sampleN_short.append('MuonEGRun2016Bver2')
datasets.append('/MuonEG/Run2016B-05Feb2018_ver2-v1/NANOAOD')
NumSample.append('-20'); sampleN_short.append('MuonEGRun2016Bver2')
datasets.append('/MuonEG/Run2016C-05Feb2018-v1/NANOAOD')
NumSample.append('-21'); sampleN_short.append('MuonEGRun2016C')
datasets.append('/MuonEG/Run2016D-05Feb2018-v1/NANOAOD')
NumSample.append('-22'); sampleN_short.append('MuonEGRun2016D')
datasets.append('/MuonEG/Run2016E-05Feb2018-v1/NANOAOD')
NumSample.append('-23'); sampleN_short.append('MuonEGRun2016E')
datasets.append('/MuonEG/Run2016F-05Feb2018-v1/NANOAOD')
NumSample.append('-24'); sampleN_short.append('MuonEGRun2016F')
datasets.append('/MuonEG/Run2016G-05Feb2018-v1/NANOAOD')
NumSample.append('-25'); sampleN_short.append('MuonEGRun2016G')
datasets.append('/MuonEG/Run2016H-05Feb2018_ver2-v1/NANOAOD')
NumSample.append('-26'); sampleN_short.append('MuonEGRun2016Hver2')
datasets.append('/MuonEG/Run2016H-05Feb2018_ver3-v1/NANOAOD')
NumSample.append('-27'); sampleN_short.append('MuonEGRun2016Hver3')


#doTT=False; doDY=False; doVV=False; doSingleT=False; doWjets=False; dottV=False
# TT
#datasets.append('/TTTo2L2Nu_13TeV-powheg/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/MINIAODSIM')
#NumSample.append('13'); sampleN_short.append('TT')
# DY
#datasets.append('/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
datasets.append('/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
NumSample.append('14'); sampleN_short.append('DY')
datasets.append('/DYToLL_0J_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
NumSample.append('15'); sampleN_short.append('DY')
datasets.append('/DYToLL_1J_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
NumSample.append('16'); sampleN_short.append('DY')
datasets.append('/DYToLL_2J_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
NumSample.append('17'); sampleN_short.append('DY')
# VV
datasets.append('/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
NumSample.append('18'); sampleN_short.append('VV')
datasets.append('/ZZTo2L2Nu_13TeV_powheg_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
NumSample.append('19'); sampleN_short.append('VV')
datasets.append('/ZZTo4L_13TeV_powheg_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
NumSample.append('20'); sampleN_short.append('VV')
datasets.append('/WWToLNuQQ_aTGC_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
NumSample.append('21'); sampleN_short.append('VV')
#datasets.append('/WWTo2L2Nu_MWW-600To800_aTGC_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
#NumSample.append('22'); sampleN_short.append('VV') ### not available now
datasets.append('/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
NumSample.append('23'); sampleN_short.append('VV')
#datasets.append('/WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
#NumSample.append('24'); sampleN_short.append('VV') ### not available now 
datasets.append('/WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v3/MINIAODSIM')
NumSample.append('25'); sampleN_short.append('VV')
datasets.append('/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
NumSample.append('26'); sampleN_short.append('VV')
datasets.append('/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
##sT
NumSample.append('27'); sampleN_short.append('sT')
datasets.append('/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
NumSample.append('28'); sampleN_short.append('sT')
datasets.append('/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
NumSample.append('29'); sampleN_short.append('sT')
datasets.append('/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
NumSample.append('30'); sampleN_short.append('sT')
datasets.append('/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
NumSample.append('31'); sampleN_short.append('sT')
# W + Jets
datasets.append('/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
NumSample.append('32'); sampleN_short.append('Wjet')
datasets.append('/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM')
NumSample.append('33'); sampleN_short.append('Wjet')
datasets.append('/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM')
NumSample.append('34'); sampleN_short.append('Wjet')
datasets.append('/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
NumSample.append('35'); sampleN_short.append('Wjet')
#datasets.append('/WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
#NumSample.append('36'); sampleN_short.append('Wjet')### not available now
datasets.append('/WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
NumSample.append('37'); sampleN_short.append('Wjet')
datasets.append('/WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
NumSample.append('38'); sampleN_short.append('Wjet')
datasets.append('/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
NumSample.append('39'); sampleN_short.append('Wjet')
# tt + V
datasets.append('/TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
NumSample.append('40'); sampleN_short.append('ttV')
datasets.append('/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM')
NumSample.append('41'); sampleN_short.append('ttV')
datasets.append('/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
NumSample.append('42'); sampleN_short.append('ttV')
datasets.append('/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext3-v1/MINIAODSIM')
NumSample.append('43'); sampleN_short.append('ttV')



jobdir = "checkoutAllSamples"
outdir = "/fdata/hepx/store/user/taohuang/NANOAOD/"
os.system("mkdir -p %s" % jobdir)
#os.system("mkdir -p %s" % outputdir)
submitscript = open("submitall_%s.sh"%(jobdir),"w")
submitscript.write("""#!/bin/bash
cd $CMSSW_BASE/src/HhhAnalysis/python/NanoAOD/
	""")

#stakeholder , background
squeues = ["stakeholder-4g","background","background-4g"]
queue = "background"
# kepperror > /dev/null
#drop error &> /dev/null
#for job in benchmarks:
outfile = "Samplelists.txt"
for ijob, job in enumerate(datasets):
    nsample = int(NumSample[ijob])
    jobtype = sampleN_short[ijob]
    dataname = ""
    datadir = " "
    #print "nsample ",nsample, " jobtype ",jobtype
    if nsample < 0:
        datadir = sampleN_short[ijob]
	dataname = job
        #print "real data nsample ",nsample, " datadir ",datadir
    else:
        datadir = job.split('/')[1]
        #print "MC nsample ",nsample, " datadir ",datadir, "MiniAOD dataset ",job.split('/')
	query = "dataset dataset=/%s/*/NANOAODSIM"%(datadir)
        pdata = os.popen("dasgoclient -limit=0 -query='{query}'".format(query = query))	
        founddataset = False
	for line in pdata:
	    #print "dataset ",line," datatype ",datadir
	    if datadir in line:
	        founddataset = True
	        dataname = line[:-1]	
	if not(founddataset): 
	    print "WARNING!!!!! no dataset found for ",datadir
    	
    outputdir = os.path.join(outdir, datadir) 
    if not os.path.exists(outputdir):
	    os.makedirs(outputdir)

    print "finaldataset ",dataname, " outputdir ",outputdir
    #os.system("echo {datadir} >> {outfilename}".format(datadir = datadir,  outfilename = outfile))
    query = "file dataset="+dataname
    #flist = os.peopn("dasgoclient -limit=0 -query='{query}' >> {outfilename}".format(query = query, outfilename = outfile))
    flist = os.popen("dasgoclient -limit=0 -query='{query}'".format(query = query), 'r')
    for line in flist:
	print "line ",line[:-1]
	if ".root" in line and "NANOAOD" in line:
	    print "download... "
	    os.system("xrdcp root://cms-xrd-global.cern.ch/"+line[:-1]+" "+outputdir)
    jobscript = open("{0}/Send_downloadDAS_{1}.slrm".format(jobdir, datadir), "w")
    jobscript.write("""#!/bin/bash
#SBATCH -J {jobtype}
#SBATCH -p {partition}
#SBATCH -n1
#SBATCH --mem-per-cpu=2000
#SBATCH --time=72:00:00
#SBATCH -o {jobdir}/batchjobs_{jobtype}-%A-%a.out
#SBATCH -e {jobdir}/batchjobs_{jobtype}-%A-%a.err
#SBATCH --ntasks-per-core=1

echo "starting at `date` on `hostname`"
echo "SLURM_JOBID=$SLURM_JOBID"
jobid=$SLURM_JOBID
source ~/.bashrc
. /etc/profile.d/modules.sh
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd $CMSSW_BASE/src/HhhAnalysis/python/NanoAOD/
eval `scramv1 runtime -sh`
export X509_USER_PROXY=$HOME/.x509up_u1468
export XRD_NETWORKSTACK=IPv4
voms-proxy-info -all
echo $X509_USER_PROXY
python downloadDASFiles.py -d {dataset} -o {output}
echo "job$jobid starts, `date`"
echo "job$jobid is done, `date`"
exit 0""".format(dataset = dataname, output=outputdir, partition=queue, jobtype = jobtype, jobdir =jobdir))
    jobscript.close()

    #os.system("cat %s/Send_PlotterProducer_%s_%d_%d.slrm"%(jobdir, job, ijob, iters))
    submitscript.write("""
sbatch {0}/Send_downloadDAS_{1}.slrm""".format(jobdir, datadir))
submitscript.close()
os.system("chmod +x submitall_%s.sh"%(jobdir))

#python PlotterProducer.py -b HaddNo {jobtype}
#os.system("./submitallHME.sh")

