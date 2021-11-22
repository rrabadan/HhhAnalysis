import os
import ROOT
import re
import numpy as np
import datetime
import sys
import glob

mainSys = ["scale_mu","PS","pile"]
eventcat = ["boosted_DY_VVV", "boosted_GGF", "boosted_other", "resolved1b_GGF", "resolved2b_GGF","resolved_DY_VVV","resolved_other"]
runyears = ["2016", "2017","2018"]

folder = "/afs/cern.ch/work/t/tahuang/CombinedLimit/ForGithub/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/Rebin_Florian/"
	    

def generate_scripts_onebinning(spin, masspoints, inputfolder, scriptsuffix,  runscript):
    #scriptsuffix = "s1fb_tm1_autoMC10" 
    condordir  = os.path.join(folder, "combinecondor")
    if not os.path.exists(condordir):
        print("creating condor dir", condordir)
        os.system("mkdir "+condordir)
    condor_all_fname = os.path.join(condordir, "condor_submit_reprocess_florian_allmass_allcat_FR_spin%d.sh"%spin)
    condor_all = open(condor_all_fname, "write")
    for mass in masspoints:
        filefolder = "/eos/user/t/tahuang/FR2_Florian_20211105/spin%d"%spin
        #allfiles = glob.glob(os.path.join(inputfolder,"*.txt"))
        #print("mass ",mass, " all txt ",allfiles)
        script_name = "spin{spin}_{mass}_FR_{scriptsuffix}".format(mass=mass,scriptsuffix=scriptsuffix, spin=spin)
        file_script = os.path.join(inputfolder, script_name+".sh")
        script = open(file_script, "write")	    
        script.write("#!/bin/bash\n")
        script.write("echo 'start to process "+file_script +"'\n")
        script.write("pushd "+inputfolder+"\n")
        txt_fname = os.path.join(filefolder, "datacard_hh_bbww_dl_all_{mass}_spin{spin}_FR2.txt".format(mass=mass,
        spin=spin))
        ws_fname = "%s_combine_workspace.root"%(script_name)
        script.write("eval `scramv1 runtime -sh`\n")
        script.write("# If workspace does not exist, create it once\n")
        #script.write("if [ ! -f  %s ]; then\n"%(ws_fname))
        script.write("text2workspace.py {txtfile} -m {mass} -o {wsname}_combine_workspace.root &> {wsname}_text2ws.log\n".format(mass = mass,
        txtfile=txt_fname, wsname = script_name))
        #script.write("fi\n\n")
        script.write("ValidateDatacards.py {txtfile} --mass {mass} --printLevel 3 --jsonFile validateDatacard_{wsname}.json \n".format(mass = mass, txtfile =txt_fname, wsname = script_name))
        script.write("# Run limit\n")
        #script.write("combine -M AsymptoticLimits --rMax 500 -m {mass} -n {prefix}_M{mass}_{ch} {prefix}_M{mass}_{ch}_combine_workspace.root &> {prefix}_M{mass}_{ch}_{suffix}.log\n".format(mass = mass, ch = channel,  suffix =scriptsuffix,  prefix = prefix))
        #script.write("combine -M AsymptoticLimits --rMax 500 -m {mass} -n {prefix}_M{mass}_{ch} {prefix}_M{mass}_{ch}_combine_workspace.root -s 1 -t 100   &> {prefix}_M{mass}_{ch}_{suffix}.log\n".format(mass = mass, ch = channel,  suffix =scriptsuffix,  prefix = prefix))
        #script.write("combine  --rMax 500 -m {mass} -n {prefix}_M{mass}_{ch} {prefix}_M{mass}_{ch}_combine_workspace.root \n".format(mass = mass, ch = channel,  suffix =scriptsuffix,  prefix = prefix))
        #script.write("combine -M Asymptotic --rMax 500 -m {mass} -n {prefix}_M{mass}_{ch} {prefix}_M{mass}_{ch}_combine_workspace.root -S 1 &> {prefix}_M{mass}_{ch}_{suffix}.log\n".format(mass = mass, ch = channel,  suffix =scriptsuffix,  prefix = prefix))

        ### add --cminDefaultMinimizerTolerance 0.5 and --cminDefaultMinimizerStrategy 0
        #script.write("combine -M AsymptoticLimits -m {mass} -n {prefix}_M{mass}_{ch} {prefix}_M{mass}_{ch}_combine_workspace.root -s $i --cminDefaultMinimizerTolerance 0.5 --cminDefaultMinimizerStrategy 0 --verbose 10  &> {prefix}_M{mass}_{ch}_{suffix}.log\n".format(mass = mass, ch = channel,  suffix =scriptsuffix,  prefix = prefix))

        #### add --picky / --rMax/--rMin --strickBounds option in fit to avoid minimization error
        #script.write("combine -M AsymptoticLimits -m {mass} -n {prefix}_M{mass}_{ch} {prefix}_M{mass}_{ch}_combine_workspace.root -s $i -t -1 --picky --verbose 1  &> {prefix}_M{mass}_{ch}_{suffix}.log\n".format(mass = mass, ch = channel,  suffix =scriptsuffix,  prefix = prefix))
        
        ### asimov fit
        script.write("combine -M AsymptoticLimits -m {mass} -n {wsfile} {wsfile}_combine_workspace.root -s 1 -t -1 --verbose 1  &> {wsfile}_{suffix}_limits.log\n".format(mass = mass, wsfile = script_name, suffix =scriptsuffix))

        script.write("popd\n")
        script.write("echo 'finish channel %s'\n"% script_name)
        script.close()
        os.system("chmod 775 "+file_script)
        if runscript:
            batchfname = os.path.join(condordir, "Batch_%s_florian.cmd"%(script_name))
            batchscript = open(batchfname, "write")
            rqstime = "testmatch" 
            batchscript.write("Notification          = Complete\n")
            batchscript.write("""universe              = vanilla 
executable            = {script}
arguments             = no
output                = limitsout/{suffix}.$(ClusterId).$(ProcId).out
error                 = limitserr/{suffix}.$(ClusterId).$(ProcId).err
log                   = limitslog/{suffix}.$(ClusterId).$(ProcId).log
request_memory        = 4000M                                                                                                                        
+JobFlavour           = "{rqstime}"
queue
    """.format(script = file_script, suffix = script_name, rqstime = rqstime))
            batchscript.close()
            condor_all.write("condor_submit "+batchfname+"\n")
            #Notification          = Complete
            #notify_user           = taohuang@email.tamu.edu
            #workday=8h, longlunch= 2h, tomorrow=24h, testmatch=3d, nextweek=7d 
            ### lxplus batch job, 8nh = 8hour queue,
            ###  -R "pool>30000"
            ### bjobs  chekc status 
        ### end of mass loop
    condor_all.close()
    os.system("chmod u+x "+condor_all_fname)
    #os.system("source "+fname_all)


masspoints = [260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 900]
scriptsuffix = "s1fb_tm1_autoMC10" 
runscript = True
inputfolder = os.path.join(folder, "Reproduce_Florian_Run2_DL_20211108/")
if not os.path.exists(inputfolder):
    print("creating inputfolder dir", inputfolder)
    os.system("mkdir "+inputfolder)
for spin in [0]:
    generate_scripts_onebinning(spin, masspoints, inputfolder, scriptsuffix,  runscript)
