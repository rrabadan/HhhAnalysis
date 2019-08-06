#!/usr/bin/python
import os
import sys
import ROOT


jobdir = "Btaggingeff"
outdir = "/fdata/hepx/store/user/taohuang/HHNtuple_20180410_DYestimation_btaggingEff/"
#outdir = "/fdata/hepx/store/user/taohuang/HHNtuple_20180411_TTbaronly/"
inputdir =  "/fdata/hepx/store/user/taohuang/HHNtuple_20180410_DYestimation/"
os.system("mkdir -p %s" % outdir)


#stakeholder , background
squeues = ["stakeholder-4g","background","background-4g"]
queue = "background-4g"
# kepperror > /dev/null
#drop error &> /dev/null
#for job in benchmarks:
##DoubleMuon
torun_datasets = []
torun_datasets.append("TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8")
torun_datasets.append('DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8')
torun_datasets.append('DYToLL_0J_13TeV-amcatnloFXFX-pythia8')
torun_datasets.append('DYToLL_1J_13TeV-amcatnloFXFX-pythia8')
torun_datasets.append('DYToLL_2J_13TeV-amcatnloFXFX-pythia8')
torun_datasets.append('ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1')
torun_datasets.append('ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1')
torun_datasets.append('ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1')
torun_datasets.append('ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1')
torun_datasets.append('ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1')
#torun_datasets.append('/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/arizzi-RunIIFall17MiniAOD-94X-Nano01Fall17-e273b12d9f89d622a34e4bc98b05ee29/USER')
#torun_datasets.append('/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM')
print "=============================================================="
print "outputdir ",outdir
print "=============================================================="

#todonanoaod = open("todaonanoaod.txt","w")
def generateslrm(torun_datasets, scriptname, jobdir):
    os.system("mkdir -p %s" % jobdir)
    submitscript = open("submitall_%s.sh"%(jobdir),"w")
    submitscript.write("""#!/bin/bash
cd $CMSSW_BASE/src/HhhAnalysis/python/DYEstimation/
	    """)

    for ijob, job in enumerate(torun_datasets):

	jobscript = open("{0}/Send_btaggingeff_{1}.slrm".format(jobdir, job), "w")
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
cd $CMSSW_BASE/src/HhhAnalysis/python/DYEstimation/
eval `scramv1 runtime -sh`
#export X509_USER_PROXY=$HOME/x509up_u1468
#voms-proxy-info -all
#echo $X509_USER_PROXY
python {scriptname} -n {jobtype}  -o {outputdir}  -i {inputdir}

echo "job$jobid starts, `date`"
echo "job$jobid is done, `date`"
exit 0""".format( inputdir=inputdir, outputdir=outdir, partition=queue, jobtype = job, jobdir =jobdir, scriptname = scriptname))
	jobscript.close()

	submitscript.write("""
sbatch {0}/Send_btaggingeff_{1}.slrm """.format(jobdir, job))
    submitscript.close()
    os.system("chmod +x submitall_%s.sh"%(jobdir))
    #c1.SaveAs("%s.eps"%(picname))

    

def plotallNtuples(torun_datasets, plotdir):
    for ijob, job in enumerate(torun_datasets):
    	finalfile = os.path.join(outdir, "btagging_efficiency_"+job+".root")
	os.system(" python drawBTaggingEfficiency.py -i {infile} -o {outdir} -n {name}".format(infile = finalfile, outdir = plotdir, name=job))



#generateslrm(torun_datasets, "computeBTaggingEfficiency_SkimNano.py", "")
generateslrm(torun_datasets, "computeFlavorFractionsOnBDT_SkimNano.py", "DY_flavour_fraction")
#plotdir = "plots_btaggingEff/"
#plotallNtuples(torun_datasets, plotdir)
