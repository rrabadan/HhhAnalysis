#!/usr/bin/python
import os
import sys
#Runyear, Slist from Samplelist
from RunConfiguration  import *
    
import ROOT


jobdir = "producerbbWW%d"%Runyear
#jobdir = "producerbbWWcounter"
#jobdir = "producerbbWWDYestimation"
outdir = "/fdata/hepx/store/user/taohuang/HHNtuple_Run%d_20190304_addSys/"%(Runyear)
#outdir = "/fdata/hepx/store/user/taohuang/HHNtuple_20180618_DYEstimation/"
#outdir = "/fdata/hepx/store/user/taohuang/HHNtuple_20180620_MC_eventcounter/"
Nanoaoddir_tao = "/fdata/hepx/store/user/taohuang/NANOAOD/"
os.system("mkdir -p %s" % jobdir)
os.system("mkdir -p %s" % outdir)


#stakeholder , background
squeues = ["stakeholder-4g","background","background-4g"]
queue = "background-4g"
#queue = "stakeholder"
dodata = True
# kepperror > /dev/null
#drop error &> /dev/null
#for job in benchmarks:
##DoubleMuon
torun_datasets = []
#torun_datasets.append("/DoubleEG/Run2016B-05Feb2018_ver1-v1/NANOAOD")
#torun_datasets.append("/DoubleEG/Run2016B-05Feb2018_ver2-v1/NANOAOD")
#torun_datasets.append("/DoubleEG/Run2016C-05Feb2018-v1/NANOAOD")
#torun_datasets.append("/DoubleEG/Run2016D-05Feb2018-v1/NANOAOD")
#torun_datasets.append("/DoubleEG/Run2016E-05Feb2018-v1/NANOAOD")
#torun_datasets.append("/DoubleEG/Run2016F-05Feb2018-v1/NANOAOD")
#torun_datasets.append("/DoubleEG/Run2016G-05Feb2018-v1/NANOAOD")
#torun_datasets.append("/DoubleEG/Run2016H-05Feb2018_ver2-v1/NANOAOD")
#torun_datasets.append("/DoubleEG/Run2016H-05Feb2018_ver3-v1/NANOAOD")
#torun_datasets.append('/DoubleMuon/Run2016B-05Feb2018_ver1-v1/NANOAOD')
#torun_datasets.append('/DoubleMuon/Run2016B-05Feb2018_ver2-v1/NANOAOD')
#torun_datasets.append('/DoubleMuon/Run2016C-05Feb2018-v1/NANOAOD')
#torun_datasets.append('/DoubleMuon/Run2016D-05Feb2018-v1/NANOAOD')
#torun_datasets.append('/DoubleMuon/Run2016E-05Feb2018-v1/NANOAOD')
#torun_datasets.append('/DoubleMuon/Run2016F-05Feb2018-v1/NANOAOD')
#torun_datasets.append('/DoubleMuon/Run2016G-05Feb2018-v1/NANOAOD')
#torun_datasets.append('/DoubleMuon/Run2016H-05Feb2018_ver2-v1/NANOAOD')
#torun_datasets.append('/DoubleMuon/Run2016H-05Feb2018_ver3-v1/NANOAOD')
#torun_datasets.append("/MuonEG/Run2016B-05Feb2018_ver1-v1/NANOAOD")
#torun_datasets.append("/MuonEG/Run2016B-05Feb2018_ver2-v1/NANOAOD")
#torun_datasets.append("/MuonEG/Run2016C-05Feb2018-v1/NANOAOD")
#torun_datasets.append("/MuonEG/Run2016D-05Feb2018-v1/NANOAOD")
#torun_datasets.append("/MuonEG/Run2016E-05Feb2018-v1/NANOAOD")
#torun_datasets.append("/MuonEG/Run2016F-05Feb2018-v1/NANOAOD")
#torun_datasets.append("/MuonEG/Run2016G-05Feb2018-v1/NANOAOD")
#torun_datasets.append("/MuonEG/Run2016H-05Feb2018_ver2-v1/NANOAOD")
#torun_datasets.append("/MuonEG/Run2016H-05Feb2018_ver3-v1/NANOAOD")
#for mass in masspoints:
#    torun_datasets.append("/GluGluToRadionToHHTo2B2VTo2L2Nu_M-%d_narrow_13TeV-madgraph-v2/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM"%mass)
####TTbar
###torun_datasets.append("/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/arizzi-RunIIFall17MiniAOD-94X-Nano01Fall17-e273b12d9f89d622a34e4bc98b05ee29/USER")
#torun_datasets.append('/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
#torun_datasets.append('/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
#torun_datasets.append('/DYToLL_0J_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
#torun_datasets.append('/DYToLL_1J_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
#torun_datasets.append('/DYToLL_2J_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
#torun_datasets.append('/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
#torun_datasets.append('/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
#torun_datasets.append('/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
#torun_datasets.append('/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
#torun_datasets.append('/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
##no LHEPDFweight
#torun_datasets.append('/WWToLNuQQ_aTGC_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
#torun_datasets.append('/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
#torun_datasets.append('/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM')
torun_datasets = Slist.Nanodatasets
#for mass in Slist.masspoints:
#    torun_datasets.append("/GluGluToRadionToHHTo2B2VTo2L2Nu_M-%d_narrow_13TeV-madgraph_correctedcfg/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM"%mass)
#print "=============================================================="
#print "outputdir ",outdir
#print "=============================================================="

hepxdir = "/fdata/hepx/"

todonanoaod = open("tod0nanoaod2016.txt","w")
def generateslrm(torun_datasets):
    submitscript = open("submitall_%s.sh"%(jobdir),"w")
    submitscript.write("""#!/bin/bash
cd $CMSSW_BASE/src/HhhAnalysis/python/NanoAOD/
	    """)

    for ijob, job in enumerate(torun_datasets):
	index = Slist.Nanodatasets.index(job)
	nsample = int(Slist.NumSample[index])
	jobtype = Slist.sampleN_short[index]
	Nanodataset =  Slist.Nanodatasets[index]
	#print "nsample ",nsample, " jobtype ",jobtype, "dataset ", job," NanoAOD ",Nanodataset
	if job.split('/')[1] != Nanodataset.split('/')[1]:
	    print "dataset matching is wrong!! job is ",job," NanoAOD is ",Nanodataset
	if nsample < 0:
	    datadir = Slist.sampleN_short[index]
	    dataname = job
	    if not dodata:
	       continue
	    print "real data nsample ",nsample, " datadir ",datadir
	elif nsample>= 0:
	    datadir = job.split('/')[1]
	    print "MC nsample ",nsample, " datadir ",datadir, "dataset ",job.split('/')

	localdir = os.path.join(Nanoaoddir_tao, datadir)
        if os.path.isdir(localdir):
	    print "localdir ",localdir
	    todonanoaod.write(Nanodataset+"\n")
	else:
	    print "warning!! no localdir ",Nanodataset

	outputdir = os.path.join(outdir, datadir)
	if not os.path.isdir(outputdir):
	    os.system("mkdir "+outputdir)


	query = "file dataset="+Nanodataset
	#print "query ",query
	#os.system("dasgoclient -limit=0 -query='{query}' >> {outfilename}".format(query = query, outfilename = outfile))
	#flist = os.peopn("dasgoclient -limit=0 -query='{query}' >> {outfilename}".format(query = query, outfilename = outfile))
	flist = os.popen("dasgoclient -limit=0 -query='{query}'".format(query = query), 'r')
	#if jobtype == "TT":
	#    print "inputdir ",inputdir
	#    flist = os.popen("ls "+inputdir)
	ifile = 0
	for line in flist:
	    ifile += 1
	    print "line ",line
	    if ".root" in line :#and "NANOAOD" in line:
	    	infile = ""
	    	file_das = line[:-1]
		rfilename = file_das.split('/')[-1]
		infile = os.path.join(hepxdir, file_das[1:])
	    	#print "torunjob ",datadir, " ",ifile, " filename  ",line[:-1]," infile ",infile
              
		if not os.path.isfile(infile):
		    infile = os.path.join(localdir,  rfilename)
		    if not os.path.isfile(infile):
			print "Error, file is not found on /fdata/hepx! not transfered ? ",infile
			print "datasetname ", Nanodataset

    		print "final infile ",infile
		#print "file on brazos ",infile
		    #else:
		 	#print "find file on brazos by transfer ",infile

	#for ifile, infile in enumerate(inputfiles):
	#	print "infile ",infile

		jobscript = open("{0}/Send_producerbbWW_{1}_{2}.slrm".format(jobdir, datadir, ifile), "w")
		jobscript.write("""#!/bin/bash
#SBATCH -J {jobtype}
#SBATCH -p {partition}
#SBATCH -n1
#SBATCH --mem-per-cpu=2000
#SBATCH --time=12:00:00
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
#export X509_USER_PROXY=$HOME/x509up_u1468
#voms-proxy-info -all
#echo $X509_USER_PROXY
python postproc_batch.py -i {inputdir} -o {outputdir} -j {jobtype}
#python postproc_eventcounter.py -i {inputdir} -o {outputdir} -j {jobtype}
#python postproc_DY_batch.py -i {inputdir} -o {outputdir} -j {jobtype}

echo "job$jobid starts, `date`"
echo "job$jobid is done, `date`"
exit 0""".format( inputdir=infile, outputdir=outputdir, partition=queue, jobtype = jobtype, jobdir =jobdir))
		jobscript.close()

		submitscript.write("""
sbatch {0}/Send_producerbbWW_{1}_{2}.slrm""".format(jobdir, datadir, ifile))
    submitscript.close()
    #os.system("chmod +x submitall_%s.sh"%(jobdir))


def mergeoutputNtuples(torun_datasets):
    overwrite = True
    for ijob, job in enumerate(torun_datasets):
	index = Slist.Nanodatasets.index(job)
	nsample = int(NumSample[index])
	jobtype = sampleN_short[index]
	print "nsample ",nsample, " jobtype ",jobtype
	if nsample < 0:
	    if not dodata:
	       continue
	    datadir = sampleN_short[index]
	    dataname = job
	    #print "real data nsample ",nsample, " datadir ",datadir
	elif nsample>= 0:
	    datadir = job.split('/')[1]
	    #print "MC nsample ",nsample, " datadir ",datadir, "MiniAOD dataset ",job.split('/')
	#faileddatasets = ["WWToLNuQQ_aTGC_13TeV-madgraph-pythia8", "ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1"]
	#faileddatasets.append("ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1")
	#if datadir in faileddatasets:
	#    continue

	outputdir = os.path.join(outdir, datadir)
    	finalfile = os.path.join(outdir, datadir+"_Run2016.root")
	if os.path.isfile(finalfile) and not(overwrite):
	    continue
	xsec = MCxsections[index]
	print "sample ",datadir, " cross section ",xsec
	haddfiles = True; write_xsec = True
	if haddfiles:
	    os.system("python haddnano.py "+finalfile+ " " +outputdir+"/*Friend.root ")
	if write_xsec and nsample > 0:
	    tfile = ROOT.TFile(finalfile,"UPDATE")
	    #h_cutflow = tfile.Get("h_cutflow")
	    h_cutflow = tfile.Get("CountFullWeighted")
	    event_weight_sum = h_cutflow.GetBinContent(1)
	    p = ROOT.TParameter(float)("cross_section", xsec)
	    p2 = ROOT.TParameter(float)("event_weight_sum", event_weight_sum)
	    p2.Write()
	    p.Write()
            tfile.Close()
    if dodata:
        for ch in ["DoubleMuon","DoubleEG","MuonEG"]:
	    finalfile = os.path.join(outdir, ch+"_Run2016.root")
	    alldatafiles = ch+"*_Friend.root"
	    os.system("python haddnano.py "+finalfile +" "+outdir+alldatafiles)



generateslrm(torun_datasets)
#mergeoutputNtuples(torun_datasets)
