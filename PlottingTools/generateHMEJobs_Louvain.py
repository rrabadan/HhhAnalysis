#!/usr/bin/python
import os
import sys


#benchmarks =  ["ttV","Wjet","sT","VV","DY","TT","Data"]
#benchmarks = ["ttV","Wjet","sT","VV","DY","TT","Data","radion","graviton"]
#benchmarks = ["radion_M400","graviton_M400"]
benchmarks = []
masses = [260, 270, 300, 350, 400, 450, 500, 550, 600, 650,700, 800, 900]
masses = [400]
for mass in masses:
	benchmarks.append("radion_M%d"%mass)
	benchmarks.append("graviton_M%d"%mass)
#alliters = [100, 1000, 10000, 100000, 1000000]
#benchmarks = []
#benchmarks.append("TT")
outputdir = "/fdata/hepx/store/user/taohuang/HHNtuple_Louvain_Test/"
alliters = [10000]
jobdir = "HMEJobs_Louvain_Test"
#jobdir = "TTJobs_Louvain"
totaljobs = 100
os.system("mkdir -p %s" % jobdir)
os.system("mkdir -p %s" % outputdir)
submitscript = open("submitall%s.sh"%jobdir,"w")
submitscript.write("""#!/bin/bash
cd $CMSSW_BASE/src/HhhAnalysis/PlottingTools/
rm -rf HADD/*txt	
	""")

#stakeholder , background
for job in benchmarks:
  for iters in alliters:
    totaljobs = iters/100
    if totaljobs>200 or job == "TT":
    	totaljobs = 2000
    print "totaljobs ",totaljobs," job ",job
    for ijob in range(0, totaljobs):
	jobscript = open("{0}/Send_PlotterProducer_{1}_{2}_{3}.slrm".format(jobdir, job, ijob, iters), "w")
	jobscript.write("""#!/bin/bash
#SBATCH -J {ijob}{jobtype}
#SBATCH -p stakeholder
#SBATCH -n1
#SBATCH --mem-per-cpu=2000
#SBATCH --time=72:00:00
#SBATCH -o batchjobs_{jobtype}_{ijob}-%A-%a.out
#SBATCH -e batchjobs_{jobtype}_{ijob}-%A-%a.err
#SBATCH --ntasks-per-core=1

echo "starting at `date` on `hostname`"
echo "SLURM_JOBID=$SLURM_JOBID"
jobid=$SLURM_JOBID
source ~/.bashrc
. /etc/profile.d/modules.sh
cd $CMSSW_BASE/src/HhhAnalysis/PlottingTools/
python runHME_Louvain.py -n {njobs} -i {ijob} -jt {jobtype} -it {iters} -o {output} 
echo "job$jobid starts, `date`"
echo "job$jobid is done, `date`"
exit 0""".format(jobtype=job, njobs=totaljobs, ijob=ijob, iters=iters, output=outputdir))
	jobscript.close()

	os.system("cat %s/Send_PlotterProducer_%s_%d_%d.slrm"%(jobdir, job, ijob, iters))
	submitscript.write("""
sbatch {0}/Send_PlotterProducer_{1}_{2}_{3}.slrm""".format(jobdir, job, ijob, iters))
submitscript.close()
os.system("chmod +x submitall%s.sh"%jobdir)

#python PlotterProducer.py -b HaddNo {jobtype}
#os.system("./submitallHME.sh")

