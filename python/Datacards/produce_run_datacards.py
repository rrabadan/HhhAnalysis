import os
import ROOT
import re
import numpy as np
import datetime
import sys
import glob
import limithelper

mainSys = ["scale_mu","PS","pile"]
#eventcat = ["boosted_DY_VVV", "boosted_GGF", "boosted_other", "resolved1b_GGF", "resolved2b_GGF","resolved_DY_VVV","resolved_other"]
eventcat = ["inclusive_DY_VVV", "boosted_GGF", "boosted_other", "resolved1b_GGF", "resolved2b_GGF","resolved_other"]

runyears = ["2016", "2017","2018"]

localfolder = "/afs/cern.ch/work/t/tahuang/CombinedLimit/ForGithub/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/Rebin_Florian/"
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
                    line = line.replace(binname+'.root', path+binname+'_rebin1d.root')
                #elif mode == "lnN":
                #    line = line.replace('.root', '_rebin1d_lnN.root')
                #elif qtlbinning:
                #    line = line.replace('.root', '_rebin1d_qtlbin.root')
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
        #originpath = "/eos/user/t/tahuang/FR2_Florian/datacard_fit_Resonant_LowMass_Graviton_2D_syst_M_%d_raw_FR2/"%(mass)
        #if mass < 550:
        #    inputfolder = folder+"datacard_fit_Resonant_LowMass_Graviton_2D_syst_M_%d_raw_FR2%s"%(mass, binsuffix)
        #else:
        #    inputfolder = folder+"datacard_fit_Resonant_HighMass_Graviton_2D_syst_M_%d_raw_FR2%s"%(mass, binsuffix)
        #    originpath =  "/eos/user/t/tahuang/FR2_Florian/datacard_fit_Resonant_HighMass_Graviton_2D_syst_M_%d_raw_FR2/"%(mass)
        #originpath = "/eos/user/t/tahuang/HHbbww_resonant_DL_Rebin_v2_202111/Graviton_syst_allbenchmarks_FR2_rebin_backgroundnode_8qtlbins/"
        originpath = "../Graviton_syst_allbenchmarks_FR2_rebin_backgroundnode_8qtlbins/"
        inputfolder = "/eos/user/t/tahuang/HHbbww_resonant_DL_Rebin_v2_202111/Graviton_syst_allbenchmarks_FR2_rebin_%s"%binsuffix
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
                #filename = "HH_{mass}_{cat}_{year}".format(mass = mass, cat=cat, year=year)
                filename = "HH_DL_{mass}_{cat}_{year}".format(mass = mass, cat=cat, year=year)
                #if filename == "HH_800_boosted_GGF_2016":
                #   ## broken file
                #   continue 
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
        os.chdir(localfolder)
	    
def condorjob(scriptname, batchfname, jobsuffix, rqstime):
    batchscript = open(batchfname, "write")
    if rqstime == "testmatch":
        batchscript.write("Notification          = Complete\n")
    batchscript.write("""universe              = vanilla 
executable            = {script}
arguments             = no
output                = limitsout/{suffix}.$(ClusterId).$(ProcId).out
error                 = limitserr/{suffix}.$(ClusterId).$(ProcId).err
log                   = limitslog/{suffix}.$(ClusterId).$(ProcId).log
request_cpus          = 4
request_memory        = 8000M                                                                                                                        
+JobFlavour           = "{rqstime}"
queue
""".format(script = scriptname, suffix = jobsuffix, rqstime = rqstime))
    batchscript.close()


def checkcondorjob_limits_onebinning(masspoints, binsuffix, scriptsuffix, mode):
    """
    check whether limits are produced from condorjob, if not, resubmit the jobs
    """
    xtitle = "Graviton mass [GeV]"
    runyears_all = runyears+["FR"]
    #eventcat_all = eventcat+["allcat"]
    eventcat_all = ["boosted_GGF","resolved1b_GGF", "resolved2b_GGF", "allcat"]
    eventcat_all = ["allcat"]
    condordir  = os.path.join(localfolder, "combinecondor")
    if not os.path.exists(condordir):
        print("creating condor dir", condordir)
        os.system("mkdir "+condordir)
    condor_all_fname = os.path.join(condordir, "condor_resubmit_allmass_allcat_FR{mode}{binsuffix}.sh".format(mode=mode, binsuffix=binsuffix))
    condor_all = open(condor_all_fname, "write")
    inputfolder = os.path.join(localfolder, "Graviton_syst_allbenchmarks_FR2{binsuffix}".format(binsuffix=binsuffix))
    if not os.path.exists(inputfolder):
        print("Error!! script to run datacard are not found ", inputfolder)
        exit()
    for mass in masspoints:
        for year in runyears_all:
            script_year_name = "Graviton_{mass}_{year}_{scriptsuffix}{mode}".format(mass=mass, year=year, scriptsuffix=scriptsuffix, mode=mode)
            file_script_year = os.path.join(inputfolder, script_year_name+".sh")
            script_year = open(file_script_year, "write")	    
            script_year.write("#!/bin/bash\n")
            script_year.write("echo 'start to process "+file_script_year +"'\n")
            script_year.write("pushd "+inputfolder+"\n")
            allcat_found = True
            for cat in eventcat_all:
                if year == "FR" and cat != "allcat":
                    continue
                #if "HH_{mass}_{cat}_{year}".format(mass = mass, cat=cat, year=year) == "HH_800_boosted_GGF_2016":
                #   ## broken file
                #   continue 
                new_fname = "HH_DL_{mass}_{cat}_{year}".format(mass = mass, cat=cat, year=year) + "_rebin1d_%s"%mode
                fname = os.path.join(inputfolder, "Graviton_%s_%s.sh"%(new_fname, scriptsuffix))
                limitfilename = new_fname+"_rebin1d_{mode}_{scriptsuffix}_limits.log".format(mass = mass, cat=cat, year=year, mode=mode, scriptsuffix=scriptsuffix)
                #inputfolder = "datacard_fit_Resonant_LowMass_Graviton_2D_syst_M_{mass}_raw_FR2_rebin_{binsuffix}".format(binsuffix=binsuffix, mass=mass)
                #if mass >= 550:
                #    inputfolder = "datacard_fit_Resonant_HighMass_Graviton_2D_syst_M_{mass}_raw_FR2_rebin_{binsuffix}".format(binsuffix=binsuffix, mass=mass)
                logfile = os.path.join(inputfolder, limitfilename)
                if not os.path.exists(logfile):
                    print("Error!! no logfile is found ", logfile)
                    script_year.write("source "+fname+"\n")
                    allcat_found = False
                    continue
                limits = limithelper.extractlimitfromtxtfile(logfile)
                if len(limits.keys()) != 6:
                    print("Warning! not all limits for this mass found ", logfile, limits)
                    script_year.write("source "+fname+"\n")
                    allcat_found = False
                    continue
            script_year.write("popd\n")
            script_year.close()
            os.system("chmod 775 "+file_script_year)
            if not allcat_found:
                print("regenerate condor job for  M ", mass, " runyear ", year)
                rqstime = "workday"
                if  year == "FR" or mode == "allsys":
                    rqstime = "testmatch" 
                batchfname = os.path.join(condordir, "Batch_resubmit_%s%s.cmd"%(script_year_name, binsuffix))
                condorjob(file_script_year, batchfname, script_year_name+binsuffix, rqstime)
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

    return 0


def generate_scripts_onebinning(binsuffix, masspoints, scriptsuffix,  runscript,  mode="lnN"):
    #eospath = "/eos/user/t/tahuang/FR2_Florian_v1_rebin/"
    eospath = "/eos/user/t/tahuang/HHbbww_resonant_DL_Rebin_v2_202111/"
    #scriptsuffix = "s1fb_tm1_autoMC10" 
    runyears_all = runyears + ["FR"]
    #eventcat_all = eventcat+["allcat"]
    #eventcat_all = ["boosted_GGF", "resolved1b_GGF", "resolved2b_GGF", "allcat"]
    eventcat_all = ["allcat"]
    condordir  = os.path.join(localfolder, "combinecondorv2")
    if not os.path.exists(condordir):
        print("creating condor dir", condordir)
        os.system("mkdir "+condordir)
    condor_all_fname = os.path.join(condordir, "condor_submit_allmass_allcat_FR_{binsuffix}_{mode}.sh".format(mode=mode, binsuffix=binsuffix))
    condor_all = open(condor_all_fname, "write")
    dirname = "Graviton_syst_allbenchmarks_FR2_rebin_{binsuffix}".format(binsuffix=binsuffix)
    inputfolder = os.path.join(localfolder, dirname)
    if not os.path.exists(inputfolder):
        print("creating folder for rebin analysis", inputfolder)
        os.mkdir(inputfolder)
    eosfolder = os.path.join(eospath, dirname) 
    if not os.path.exists(eosfolder):
        print("Error!! rebinned shape folder is not found", eosfolder)
        exit()
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
            #script_year.write("pushd "+inputfolder+"\n")
            for cat in eventcat_all:
                if year == "FR" and cat != "allcat":
                    continue
                filename = "HH_{mass}_{cat}_{year}".format(mass = mass, cat=cat, year=year)
                #if filename == "HH_800_boosted_GGF_2016":
                #   ## broken file
                #   continue 
                new_fname = filename+"_rebin1d_%s"%mode
                dcfile = os.path.join(eosfolder, new_fname+".txt")
                #########################################3
                ## script to run data card
                fname = os.path.join(inputfolder, "Graviton_%s_%s.sh"%(new_fname, scriptsuffix))
                script_year.write("source "+fname+"\n")
                script = open(fname, "write")
                script.write("#!/bin/bash\n")
                script.write("echo 'start  %s'\n"% new_fname)
                script.write("pushd "+inputfolder+"\n")
                script.write("eval `scramv1 runtime -sh`\n")
                script.write("ValidateDatacards.py {dcfile} --mass {mass} --printLevel 3 --jsonFile  validateDatacard_{fname}.json \n".format(mass = mass, dcfile = dcfile, fname=new_fname))
                script.write("# If workspace does not exist, create it once\n")
                #script.write("if [ ! -f  %s ]; then\n"%(ws_fname+"_combine_workspace.root"))
                script.write("text2workspace.py {dcfile} -m {mass} -o {ws_fname}_combine_workspace.root &> {ws_fname}_text2ws.log\n".format(mass =
                mass, dcfile=dcfile, ws_fname=new_fname))
                #script.write("fi\n\n")
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
            #script_year.write("popd\n")
            script_year.close()
            os.system("chmod 775 "+file_script_year)
            if runscript:
                rqstime = "workday"
                if  year == "FR" or mode == "allsys":
                    rqstime = "testmatch" 
                batchfname = os.path.join(condordir, "Batch_%s_rebin_%s.cmd"%(script_year_name, binsuffix))
                condorjob(file_script_year, batchfname, script_year_name+binsuffix, rqstime)
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

def generate_scripts(masspoints, scriptsuffix, runscript, mode, nn_bins, HME_bins, threholds, HMEqtl, quadraticthresh):
    bkg_qtlbinning = True
    condordir = "combinecondorv2"
    ##condor_submit_allmass_allcat_FR_nnbin10HME10qtlthresh10_lnN.sh 
    condorname = "condjor_submit"
    condorname += "i_HMqtlbin" if HMEqtl else "_HMEnorm"
    condorname += "_quadraticthresh" if quadraticthresh else "_fixedthresh"
    condorall = os.path.join(condordir, condorname+".sh")
    condorall_w = open(condorall, "w")
    condorall_w.write("echo 'starting to submit job "+condorname+"...'\n")
    for hmev in HME_bins:
      for nnv in nn_bins:
          for thresh in thresholds:
              binsuffix = "nnbin%dHME%.2fthresh%d"%(nnv, hmev, thresh)
              if HMEqtl and not quadraticthresh:
                  binsuffix = "nnbin%dHME%dqtlthresh%d"%(nnv, hmev, thresh)
              elif HMEqtl and quadraticthresh:
                  binsuffix = "nnbin%dHME%dqtlthreshquadratic"%(nnv, hmev)
              elif not HMEqtl and quadraticthresh:
                  binsuffix = "nnbin%dHME%.2fthreshquadratic"%(nnv, hmev)
              binsuffix = binsuffix.replace('.','p')
              lastjob = hmev == HME_bins[-1] and nnv == nn_bins[-1] and thresh == thresholds[-1]
              if lastjob:
                  print("creating condor job ",binsuffix)
              condorall_w.write("source condor_submit_allmass_allcat_FR_"+binsuffix+"_"+mode+".sh\n")
              #produce_datacard_onebinning(binsuffix, masspoints, bkg_qtlbinning, mode)
              #generate_scripts_onebinning(binsuffix, masspoints, scriptsuffix, runscript, mode)
    condorall_w.close()
    os.system("chmod u+x "+condorall)

def generate_scripts_v1(masspoints, scriptsuffix, runscript, mode, nn_bins, HME_bins, qtlbinning):
    bkg_qtlbinning = True
    for nnbin in nn_bins:
        for hmebin in HME_bins:
            binsuffix = "nnbin%sHMEv%d"%(nnbin, hmebin)
            if qtlbinning:
                binsuffix = binsuffix+"qltbin"
            #produce_datacard_onebinning(binsuffix, masspoints, bkg_qtlbinning, mode)
            generate_scripts_onebinning(binsuffix, masspoints, scriptsuffix, runscript, mode)
            #checkcondorjob_limits_onebinning(masspoints, binsuffix, scriptsuffix, mode)

masspoints = [260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 900]
nn_bins = [5, 8, 10, 20]
HME_bins = [0,1,2,3]
#masspoints = [550, 600, 650, 700, 800, 900]
bkg_qtlbinning = True
qtlbinning = False
if qtlbinning:
    nn_bins = [10, 20]

masspoints = [260]

scriptsuffix = "s1fb_tm1_autoMC10" 
runscript = True
#mode = "mainSys"
mode = "lnN"
#mode = "allsys"

for i in range(4):
    HMEquantile_binning = i/2
    quadraticthreshalgo = i%2
    HME_bins = [1.1, 1.15, 1.2]
    nn_bins = [10,20, 40]
    thresholds = [2, 5, 10]
    #HME_bins = [1.1]
    #nn_bins = [10]
    #thresholds = [2]
    if HMEquantile_binning:
       HME_bins = [5, 10, 15]
    if quadraticthreshalgo:
       thresholds = [10]

    if HMEquantile_binning and quadraticthreshalgo:
        print("Running rebin mode with HME quantile binning and quadratic threshold algo, and systematics: ", mode)
    elif HMEquantile_binning and not quadraticthreshalgo:
        print("Running rebin mode with HME quantile binning and fixed threshold algo, and systematics: ", mode)
    elif not HMEquantile_binning and quadraticthreshalgo:
        print("Running rebin mode with customized HME binning and quadratic threshold algo, and systematics: ", mode)
    else:
        print("Running rebin mode with customized HME binning and fixed threshold algo, and systematics: ", mode)

    generate_scripts(masspoints, scriptsuffix, runscript, mode, nn_bins, HME_bins, thresholds, HMEquantile_binning, quadraticthreshalgo)
    
