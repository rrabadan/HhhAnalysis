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
def produce_datacard(infile, outfile, path, bkg_qtlbinning=True, mode=""):
    binname = infile.split("/")[-1][:-4]
    bkgnode = "GGF" not in binname
    qtlbinning = True
    #print("binname ",binname)
    fw = open(outfile, 'w')
    with open(infile, 'r') as f:
        startSys=False
        for line in f:
            systematic_name = line.split(' ')[0]
            mainSys_true = any(sysname in systematic_name for sysname in mainSys)
            if '.root' in line:
                if bkg_qtlbinning and bkgnode:
                    line = line.replace(binname+'.root', path+binname+'_rebin1d_qtlbin.root')
                elif mode == "lnN":
                    line = line.replace('.root', '_rebin1d_lnN.root')
                elif qtlbinning:
                    line = line.replace('.root', '_rebin1d_qtlbin.root')
                else:
                    line = line.replace('.root', '_rebin1d.root')
                #print("root file line ", line)
            if line.startswith('nuisance'):
                line = line[:-1] + " ifexists \n"
                #print("nuisance line ",line)
            if (startSys and ('lnN' in line  or (mode == "mainSys" and mainSys_true) or mode == "allsys")) or not startSys or line.startswith('---------'):
                fw.write(line) 
            if line.startswith('rate'):
                startSys = True
                #print("starting to process sys", line)
            elif line.startswith('nuisance'):
                startSys = False
        
            #print(line)
    fw.close()


def produce_datacard_onebinning(binsuffix, masspoints, bkg_qtlbinning=True, mode="lnN"):
    for mass in masspoints:
        inputfolder = None
        originpath = "/eos/user/t/tahuang/FR2_Florian/datacard_fit_Resonant_LowMass_Graviton_2D_syst_M_%d_raw_FR2/"%(mass)
        if mass < 550:
            inputfolder = folder+"datacard_fit_Resonant_LowMass_Graviton_2D_syst_M_%d_raw_FR2%s"%(mass, binsuffix)
        else:
            inputfolder = folder+"datacard_fit_Resonant_HighMass_Graviton_2D_syst_M_%d_raw_FR2%s"%(mass, binsuffix)
            originpath =  "/eos/user/t/tahuang/FR2_Florian/datacard_fit_Resonant_HighMass_Graviton_2D_syst_M_%d_raw_FR2/"%(mass)
        #allfiles = glob.glob(os.path.join(inputfolder,"*.txt"))
        #print("mass ",mass, " all txt ",allfiles)
        combine_FR = "combineCards.py " 
        os.chdir(inputfolder)
        for year in runyears:
            fname_list = []
            combine_str = "combineCards.py "
            for cat in eventcat:
            #if f.endswith("2016.txt") or f.endswith("2017.txt") or f.endswith("2018.txt"):
                #filename = f.split('/')[-1]
                filename = "HH_{mass}_{cat}_{year}".format(mass = mass, cat=cat, year=year)
                if filename == "HH_800_boosted_GGF_2016":
                   ## broken file
                   continue 
                oldf = os.path.join(inputfolder, filename+".txt")
                new_fname = filename+"_rebin1d_%s"%mode
                combine_str += " {fname}={newfname}.txt ".format(fname = filename, newfname=new_fname)
                newf = os.path.join(inputfolder, new_fname+".txt")
                produce_datacard(oldf, newf, originpath, bkg_qtlbinning, mode)
                combine_str += " > HH_{mass}_allcat_{year}_rebin1d_{mode}.txt".format(mass = mass, year=year, mode = mode)
            os.system(combine_str)
            combine_FR += " bbWW=HH_{mass}_allcat_{year}_rebin1d_{mode}.txt ".format(mass = mass, year=year, mode = mode)
        combine_FR += " > HH_{mass}_allcat_FR_rebin1d_{mode}.txt".format(mass = mass, mode = mode)
        os.system(combine_FR)
	    

def generate_scripts_onebinning(binsuffix, masspoints, scriptsuffix,  runscript,  mode="lnN"):
    eospath = "/eos/user/t/tahuang/FR2_Florian_v1_rebin/"
    #scriptsuffix = "s1fb_tm1_autoMC10" 
    runyears_all = runyears + ["FR"]
    #eventcat_all = eventcat+["allcat"]
    eventcat_all = ["boosted_GGF", "resolved1b_GGF", "resolved2b_GGF", "allcat"]
    condordir  = os.path.join(folder, "combinecondor")
    if not os.path.exists(condordir):
        print("creating condor dir", condordir)
        os.system("mkdir "+condordir)
    condor_all_fname = os.path.join(condordir, "condor_submit_allmass_allcat_FR{mode}{binsuffix}.sh".format(mode=mode, binsuffix=binsuffix))
    condor_all = open(condor_all_fname, "write")
    inputfolder = "Graviton_syst_allbenchmarks_FR2{binsuffix}".format(binsuffix=binsuffix)
    if not os.path.exists(inputfolder):
        print("creating folder for rebin analysis", inputfolder)
        os.mkdir(inputfolder)
    eosfolder = os.path.join(eospath, inputfolder) 
    for mass in masspoints:
        #inputfolder = None
        #if mass < 550:
        #    inputfolder = folder+"datacard_fit_Resonant_LowMass_Graviton_2D_syst_M_%d_raw_FR2%s"%(mass, binsuffix)
        #else:
        #    inputfolder = folder+"datacard_fit_Resonant_HighMass_Graviton_2D_syst_M_%d_raw_FR2%s"%(mass, binsuffix)
        #allfiles = glob.glob(os.path.join(inputfolder,"*.txt"))
        #print("mass ",mass, " all txt ",allfiles)
        for year in runyears_all:
            script_year_name = "Graviton_{mass}_{year}_{scriptsuffix}{mode}".format(mass=mass, year=year, scriptsuffix=scriptsuffix, mode=mode)
            file_script_year = os.path.join(inputfolder, script_year_name+".sh")
            script_year = open(file_script_year, "write")	    
            script_year.write("#!/bin/bash\n")
            script_year.write("echo 'start to process "+file_script_year +"'\n")
            script_year.write("pushd "+inputfolder+"\n")
            for cat in eventcat_all:
                if year == "FR" and cat != "allcat":
                    continue
                filename = "HH_{mass}_{cat}_{year}".format(mass = mass, cat=cat, year=year)
                if filename == "HH_800_boosted_GGF_2016":
                   ## broken file
                   continue 
                new_fname = filename+"_rebin1d_%s"%mode
                dcfile = os.path.join(eosfolder, new_fname+".txt")
                #########################################3
                ## script to run data card
                #ws_fname = os.path.join(inputfolder, "%s_combine_workspace.root"%(new_fname))
                ws_fname = "%s_combine_workspace.root"%(new_fname)
                fname = os.path.join(inputfolder, "Graviton_%s_%s.sh"%(new_fname, scriptsuffix))
                script_year.write("source "+fname+"\n")
                script = open(fname, "write")
                script.write("#!/bin/bash\n")
                script.write("echo 'start  %s'\n"% new_fname)
                script.write("pushd "+inputfolder+"\n")
                script.write("eval `scramv1 runtime -sh`\n")
                script.write("# If workspace does not exist, create it once\n")
                #script.write("if [ ! -f  %s ]; then\n"%(ws_fname))
                script.write("text2workspace.py {dcfile} -m {mass} -o {wsfile} &> {dcfile}_text2ws.log\n".format(mass =
                mass, dcfile=dcfile, wsfile=ws_fname))
                #script.write("fi\n\n")
                script.write("ValidateDatacards.py {dcfile} --mass {mass} --printLevel 3 --jsonFile  validateDatacard_{fname}.json \n".format(mass = mass, dcfile = dcfile, fname=new_fname))
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
                script.write("combine -M AsymptoticLimits -m {mass} -n {fname} {fname}_combine_workspace.root -s 1 -t -1 --verbose 1  &> {fname}_{suffix}_limits.log\n".format(mass = mass, fname = new_fname, suffix =scriptsuffix))

                script.write("popd\n")
                script.write("echo 'finish channel %s'\n"% new_fname)
                script.close()
                os.system("chmod 775 "+fname)
            ############## end of loop over all event categories
            script_year.write("popd\n")
            script_year.close()
            os.system("chmod 775 "+file_script_year)
            if runscript:
                batchfname = os.path.join(condordir, "Batch_%s%s.cmd"%(script_year_name, binsuffix))
                batchscript = open(batchfname, "write")
                rqstime = "workday"
                if  year == "FR" or mode == "allsys":
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
		""".format(script = file_script_year, suffix = script_year_name+binsuffix, rqstime = rqstime))
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
nnbins_v = [5, 8, 10, 20]
hmbins_v = [0,1,2,3]
#masspoints = [550, 600, 650, 700, 800, 900]

#masspoints = [260]
nnbins_v = [10, 20]
#hmbins_v = [1]

scriptsuffix = "s1fb_tm1_autoMC10" 
runscript = True
#mode = "mainSys"
#mode = "lnN"
mode = "allsys"
bkg_qtlbinning = True
qtlbinning = True

for mode in ["lnN","allsys"]:
    for nnbin in nnbins_v:
        for hmebin in hmbins_v:
            binsuffix = "_rebin_nnbin%sHMEv%d"%(nnbin, hmebin)
            if qtlbinning:
                binsuffix = binsuffix+"qltbin"
            #produce_datacard_onebinning(binsuffix, masspoints, bkg_qtlbinning, mode)
            generate_scripts_onebinning(binsuffix, masspoints, scriptsuffix, runscript, mode)
