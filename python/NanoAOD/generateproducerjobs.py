#!/usr/bin/python
import os
import sys
from Samplelist import *


jobdir = "producerbbWW"
outdir = "/fdata/hepx/store/user/taohuang/HHNtuple_20180322/"
Nanoaoddir = "/fdata/hepx/store/user/taohuang/NANOAOD/"
os.system("mkdir -p %s" % jobdir)


#stakeholder , background
squeues = ["stakeholder-4g","background","background-4g"]
queue = "background"
# kepperror > /dev/null
#drop error &> /dev/null
#for job in benchmarks:
##DoubleMuon
torun_datasets = []
torun_datasets.append('/DoubleMuon/Run2016B-05Feb2018_ver1-v1/NANOAOD')
torun_datasets.append('/DoubleMuon/Run2016B-05Feb2018_ver2-v1/NANOAOD')
torun_datasets.append('/DoubleMuon/Run2016C-05Feb2018-v1/NANOAOD')
torun_datasets.append('/DoubleMuon/Run2016D-05Feb2018-v1/NANOAOD')
torun_datasets.append('/DoubleMuon/Run2016E-05Feb2018-v1/NANOAOD')
torun_datasets.append('/DoubleMuon/Run2016F-05Feb2018-v1/NANOAOD')
torun_datasets.append('/DoubleMuon/Run2016G-05Feb2018-v1/NANOAOD')
torun_datasets.append('/DoubleMuon/Run2016H-05Feb2018_ver2-v1/NANOAOD')
torun_datasets.append('/DoubleMuon/Run2016H-05Feb2018_ver3-v1/NANOAOD')
masspoints = [260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 750, 800, 900]
for mass in masspoints:
    torun_datasets.append("/GluGluToRadionToHHTo2B2VTo2L2Nu_M-%d_narrow_13TeV-madgraph-v2/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM"%mass)
torun_datasets.append('/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/arizzi-RunIIFall17MiniAOD-94X-Nano01Fall17-e273b12d9f89d622a34e4bc98b05ee29/USER')

def generateslrm(torun_datasets):
    submitscript = open("submitall_%s.sh"%(jobdir),"w")
    submitscript.write("""#!/bin/bash
cd $CMSSW_BASE/src/HhhAnalysis/python/NanoAOD/
	    """)

    for ijob, job in enumerate(torun_datasets):
	index = datasets.index(job)
	nsample = int(NumSample[index])
	jobtype = sampleN_short[index]
	print "nsample ",nsample, " jobtype ",jobtype, "dataset ", job
	if nsample < 0:
	    datadir = sampleN_short[index]
	    dataname = job
	    #print "real data nsample ",nsample, " datadir ",datadir
	elif nsample>= 0:
	    datadir = job.split('/')[1]
	    #print "MC nsample ",nsample, " datadir ",datadir, "MiniAOD dataset ",job.split('/')

	inputdir = os.path.join(Nanoaoddir, datadir)
	inputfiles = []
	if os.path.isdir(inputdir):
	    #print "This is a valid directory: ",inputdir
	    ls = os.listdir(inputdir)
	    inputfiles = [os.path.join(inputdir,x) for x in ls if x.endswith("root")]
	elif os.path.isfile(inputdir):
	    inputfiles = [inputdir]

	outputdir = os.path.join(outdir, datadir)
	if not os.path.isdir(outputdir):
	    os.system("mkdir "+outputdir)

	for ifile, infile in enumerate(inputfiles):

	    jobscript = open("{0}/Send_producerbbWW_{1}_{2}.slrm".format(jobdir, jobtype, ifile), "w")
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
python postproc_batch.py -i {inputdir} -o {outputdir} -j {jobtype}
echo "job$jobid starts, `date`"
echo "job$jobid is done, `date`"
exit 0""".format( inputdir=infile, outputdir=outputdir, partition=queue, jobtype = jobtype, jobdir =jobdir))
	    jobscript.close()

	    submitscript.write("""
sbatch {0}/Send_producerbbWW_{1}_{2}.slrm""".format(jobdir, datadir, ifile))
    submitscript.close()
    os.system("chmod +x submitall_%s.sh"%(jobdir))

    ##python PlotterProducer.py -b HaddNo {jobtype}
    ##os.system("./submitallHME.sh")



def mergeoutputNtuples(torun_datasets):
    for ijob, job in enumerate(torun_datasets):
	index = datasets.index(job)
	nsample = int(NumSample[index])
	jobtype = sampleN_short[index]
	print "nsample ",nsample, " jobtype ",jobtype
	if nsample < 0:
	    datadir = sampleN_short[index]
	    dataname = job
	    #print "real data nsample ",nsample, " datadir ",datadir
	elif nsample>= 0:
	    datadir = job.split('/')[1]
	    #print "MC nsample ",nsample, " datadir ",datadir, "MiniAOD dataset ",job.split('/')
	outputdir = os.path.join(outdir, datadir)
    	finalfile = os.path.join(outdir, datadir+"_Skim.root")
    	os.system("hadd -f "+finalfile+ " " +outputdir+"/*root ")



generateslrm(torun_datasets)
