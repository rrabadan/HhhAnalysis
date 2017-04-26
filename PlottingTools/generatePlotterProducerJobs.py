#!/usr/bin/python
import os
import sys

# Jobs to send
benchmarks = ["ttV","Wjet","sT","VV","DY","TT", "Rad_260","Rad_270","Rad_300","Rad_350","Rad_400","Rad_450","Rad_500","Rad_550","Rad_600","Rad_650","Rad_750","Rad_800","Rad_900","Rad_1000","Grav_260","Grav_270","Grav_300","Grav_350","Grav_400","Grav_450","Grav_500","Grav_550","Grav_600","Grav_650","Grav_700","Grav_800","Grav_900","Grav_1000","Data"]
# Wanna do the HADD? The first time you run it is mandatory
HaddYes = "HaddYes"#HaddNo

# Folder containing the slrm jobs, one per sample
jobdir = "slrm_PlotterProducer"
os.system("mkdir -p %s" % jobdir)

# File thatyou have to source
submitscript = open("Send_PlotterProducer.sh","w")
submitscript.write("""#!/bin/bash
cd $CMSSW_BASE/src/HhhAnalysis/PlottingTools/
rm -rf HADD/*txt
""")


for sample in benchmarks:
  jobscript = open("{0}/Send_PlotterProducer_{1}.slrm".format(jobdir, sample), "w")
  jobscript.write("""#!/bin/bash
#SBATCH -J run{0}
#SBATCH -p stakeholder
#SBATCH -n1
#SBATCH --mem-per-cpu=2000
#SBATCH -o batchjobs_{0}-%A-%a.out
#SBATCH -e batchjobs_{0}-%A-%a.err
#SBATCH --ntasks-per-core=1

echo "starting at `date` on `hostname`"
echo "SLURM_JOBID=$SLURM_JOBID"
jobid=$SLURM_JOBID
source ~/.bashrc
. /etc/profile.d/modules.sh
cd $CMSSW_BASE/src/HhhAnalysis/PlottingTools/
python PlotterProducer.py -b {1} {0}
echo "job$jobid starts, `date`"
echo "job$jobid is done, `date`"
exit 0""".format(sample, HaddYes))
  jobscript.close()
  submitscript.write("""sbatch {0}/Send_PlotterProducer_{1}.slrm\n""".format(jobdir, sample))

submitscript.close()
os.system("chmod +x Send_PlotterProducer.sh")

#os.system("./Send_PlotterProducer.sh")
