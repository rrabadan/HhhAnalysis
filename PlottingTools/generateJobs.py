#!/usr/bin/python
import os
import sys


benchmarks = ["ttV", "Wjet", "singTop", "VV", "DY", "TTbar", "Data"]
jobdir = "slr"
os.system("mkdir -p %s" % jobdir)
submitscript = open("submitallPlotProducer.sh","w")
submitscript.write("""#!/bin/bash
cd $CMSSW_BASE/src/HhhAnalysis/PlottingTools/
rm -rf HADD/*txt	
	""")


for job in benchmarks:
    jobscript = open("{0}/Send_PlotterProducer_{1}.slrm".format(jobdir, job), "w")
    jobscript.write("""#!/bin/bash
#SBATCH -J runsplit
#SBATCH -p hepx
#SBATCH -n1
#SBATCH --mem-per-cpu=12000
#SBATCH -o batchjobs_{jobtype}-%A-%a.out
#SBATCH -e batchjobs_{jobtype}-%A-%a.err
#SBATCH --ntasks-per-core=1

echo "starting at `date` on `hostname`"
echo "SLURM_JOBID=$SLURM_JOBID"
jobid=$SLURM_JOBID

source ~/.bashrc
. /etc/profile.d/modules.sh
cd $CMSSW_BASE/src/HhhAnalysis/PlottingTools/
echo "job$jobid starts, `date`"
python PlotterProducer.py -b HaddYes {jobtype}
echo "job$jobid is done, `date`"
exit 0""".format(jobtype=job))
    jobscript.close()

    submitscript.write("""
sbatch {0}/submit_{1}.slrm""".format(jobdir, job))
submitscript.close()
os.system("chmod +x submitallPlotProducer.sh")

#os.system("./submitallPlotProducer.sh")

