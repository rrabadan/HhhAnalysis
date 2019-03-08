#!/usr/bin/python
import os
import sys
from Samplelist import *


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
masspoints = [260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900]
#outfile = open("Samplelists_test.txt","w")
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
    elif nsample > 0 and nsample < 14:#signal, skip now and TTbar, FIXME
        datadir = sampleN_short[ijob]
        dataname = job
    elif nsample>= 14:
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
    
    #outfile.write('Nanodatasets.append("'+dataname+'")\n')
    site = "root://cmsxrootd.fnal.gov/"##or root://cms-xrd-global.cern.ch/
    
    #if ijob == 0 or nsample != int(NumSample[ijob-1]):
    #	os.system("echo {datadir} >> {outfilename}".format(datadir = datadir,  outfilename = outfile))
    query = "file dataset="+dataname
    #os.system("dasgoclient -limit=0 -query='{query}' >> {outfilename}".format(query = query, outfilename = outfile))
    #flist = os.peopn("dasgoclient -limit=0 -query='{query}' >> {outfilename}".format(query = query, outfilename = outfile))
    flist = os.popen("dasgoclient -limit=0 -query='{query}'".format(query = query), 'r')
    ifile = 0
    for line in flist:
	print "line ",line[:-1]
	ifile += 1
	if ".root" in line and "NANOAOD" in line:
	    #print "download... "
	    #os.system("xrdcp root://cms-xrd-global.cern.ch/"+line[:-1]+" "+outputdir)
	    fileonfdata = os.path.join(outputdir, line[:-1].split("/")[-1])	
	    if os.path.isfile(fileonfdata):
		print "file on fdata ",fileonfdata
		continue
	    jobscript = open("{0}/Send_downloadDAS_{1}_{2}.slrm".format(jobdir, datadir, ifile), "w")
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
export X509_USER_PROXY=$HOME/x509up_u1468
#export XRD_NETWORKSTACK=IPv4
voms-proxy-info -all
echo $X509_USER_PROXY
echo dataset, {dataset}
echo "file to download "
echo {site}{filename}
xrdcp {site}{filename} {output}
#python downloadDASFiles.py -d {dataset} -o {output}
echo "job$jobid starts, `date`"
echo "job$jobid is done, `date`"
exit 0""".format(dataset = dataname, site = site, filename = line[:-1], output=outputdir, partition=queue, jobtype = jobtype, jobdir =jobdir))
	    jobscript.close()

    #os.system("cat %s/Send_PlotterProducer_%s_%d_%d.slrm"%(jobdir, job, ijob, iters))
	    submitscript.write("""
sbatch {0}/Send_downloadDAS_{1}_{2}.slrm""".format(jobdir, datadir, ifile))
submitscript.close()
os.system("chmod +x submitall_%s.sh"%(jobdir))

##python PlotterProducer.py -b HaddNo {jobtype}
##os.system("./submitallHME.sh")
#
