#!/usr/bin/python
import os
import sys
import ROOT


def generate_condorjob(mass, binsuffix, hmeversion, nnversion, thresh, filename, shellname, pyfile, addnotification=False):
    #condorscript = open("condor_rebin_mass%d_nnv%dHMEbinv%d.cmd"%(mass, nnv, hmev),"w")
    condorscript = open(filename, "w")
    if addnotification:
        condorscript.write("Notification          = Complete\n")
        condorscript.write("notify_user           = taohuang@tamu.edu\n")
    condorscript.write("""universe              = vanilla 
executable            = {shellname}
arguments             = {pyfile} {mass} {hmeversion} {nnversion} {thresh}
should_transfer_files = YES
transfer_input_files  = {pyfile}
output                = condorout/rebin_M{mass}_{binsuffix}.$(ClusterId).$(ProcId).out
error                 = condorerr/rebin_M{mass}_{binsuffix}.$(ClusterId).$(ProcId).err
log                   = condorlog/rebin_M{mass}_{binsuffix}.$(ClusterId).$(ProcId).log
request_cpus          = 4
request_memory        = 8000M
+JobFlavour           = "workday"
queue""".format(mass = mass, binsuffix=binsuffix, hmeversion=hmeversion, nnversion=nnversion, thresh=thresh, shellname=shellname,pyfile=pyfile))
    condorscript.close()

def genereate_condorjobs_multiqueue(mass, inputdir, inputlist, binsuffix, hmeversion, nnversion, thresh, filename, shellname, pyfile, HMEqtl, quadraticthresh, addnotification):
    #runrebin_condorsig.sh
    ##echo 'start running rebin python  ', $1, ' with HMT version ',$2, 'NN version ', $3," thresh ", $4, "quadraticthreshalgo", $5, " HMEquantile_binning ",$6, " infile ",$7, " outfile ",$8
# If workspace does not exist, create it once
    condorscript = open(filename, "w")
    if addnotification:
        condorscript.write("Notification          = Complete\n")
        condorscript.write("notify_user           = taohuang@tamu.edu\n")
    condorscript.write("""universe              = vanilla 
executable            = {shellname}
arguments             = {pyfile} {hmeversion} {nnversion} {thresh} {quadraticthresh} {HMEqtl} $(infile)
input                 = {pyfile}
output                = condorout/rebin_M{mass}_{binsuffix}.$(ClusterId).$(ProcId).out
error                 = condorerr/rebin_M{mass}_{binsuffix}.$(ClusterId).$(ProcId).err
log                   = condorlog/rebin_M{mass}_{binsuffix}.$(ClusterId).$(ProcId).log
request_cpus          = 2
request_memory        = 4000M
+JobFlavour           = "workday"
should_transfer_files = YES
transfer_input_files  = $(infile) 
""".format(mass=mass, binsuffix=binsuffix, hmeversion=hmeversion, nnversion=nnversion, thresh=thresh, shellname=shellname,pyfile=pyfile, HMEqtl=HMEqtl, quadraticthresh=quadraticthresh))
    for i,infile in enumerate(inputlist):
        infile_path = os.path.join(inputdir, infile)
        condorscript.write("infile="+infile_path+"\n")
        condorscript.write("queue 1"+"\n")

    condorscript.close()



#notify_user           = taohuang@email.tamu.edu
#Notification          = Complete

def generatejob_bkgnode(masspoints):
    pyfile = "rebin_florian_bkgnode.py"
    shellname = "runrebin_general.sh"
    submitscript = open("submit_rebincondorjobs_bkgnode.sh","w")
    #submitscript.write("rm condorout/* \n")
    #submitscript.write("rm condorerr/* \n")
    #submitscript.write("rm condorlog/* \n")
    binsuffix = "8qltbin"
    for mass in masspoints:
      filename = "condor_rebin_mass%d_bkgnode.cmd"%(mass)
      
      generate_condorjob(mass, binsuffix, 1, 1, 1, filename, shellname)
      submitscript.write("condor_submit "+filename+" \n")
    submitscript.close()
    os.system("chmod +x submit_rebincondorjobs_bkgnode.sh")



def generatejob_signode(masspoints, HME_bins, nn_bins, thresholds, HMEqtl, quadraticthresh):
    #shellname = "runrebin_general.sh"
    #pyfile = "process_rebin_florian_condorsig.py"
    shellname = "runrebin_condorsig.sh"
    pyfile = "rebin_florian.py"
    submitscriptname = "submit_rebincondorjobs_signode.sh"
    if HMEqtl and not quadraticthresh: 
        #pyfile = "rebin_florian_hmeqtl_fixedthresh.py"
        submitscriptname = "submit_rebincondorjobs_signode_hmeqtl_fixedthresh.sh"
    elif HMEqtl and quadraticthresh:
        #pyfile = "rebin_florian_hmeqtl_quadraticthresh.py"
        submitscriptname = "submit_rebincondorjobs_signode_hmeqtl_quadraticthresh.sh"
    elif not HMEqtl and quadraticthresh:
        #pyfile = "rebin_florian_hmenorm_quadraticthresh.py"
        submitscriptname = "submit_rebincondorjobs_signode_hmenorm_quadraticthresh.sh"
    elif thresholds == [0]:
        #pyfile = "rebin_florian.py"
        submitscriptname = "submit_rebincondorjobs_signode_equalNNbins.sh"
    submitscript = open(submitscriptname,"w")
    outputfolder = "/eos/user/t/tahuang/HHbbww_resonant_DL_Rebin_v2_202111/" 
    submitscript.write("echo 'starting to submit job "+submitscriptname +"'\n")
    #submitscript.write("rm condorout/* \n")
    #submitscript.write("rm condorerr/* \n")
    #submitscript.write("rm condorlog/* \n")
    eventcat = ["boosted_GGF", "resolved1b_GGF", "resolved2b_GGF"]
    runyears = ["2016", "2017","2018"]
    for mass in masspoints:
        folder = "/eos/user/t/tahuang/FR2_Florian_v2/datacard_fit_Resonant_HighMass_Graviton_2D_syst_M_%d_raw_FR2"%mass
        if mass <= 500:
            folder = "/eos/user/t/tahuang/FR2_Florian_v2/datacard_fit_Resonant_LowMass_Graviton_2D_syst_M_%d_raw_FR2"%mass
        filelist = []
        for year in runyears:
            for cat in eventcat:
                filelist.append("HH_DL_{mass}_{cat}_{year}.root".format(mass = mass, cat=cat, year=year))
        for hmev in HME_bins:
          for nnv in nn_bins:
              for thresh in thresholds:
                  binsuffix = "nnbin%dHME%.2fthresh%d"%(nnv, hmev, thresh)
                  if HMEqtl and not quadraticthresh:
                      binsuffix = "nnbin%dHME%dqtlthresh%d"%(nnv, hmev, thresh)
                  elif HMEqtl and quadraticthresh:
                      #binsuffix = "nnbin%dHME%dqtlthreshquadratic"%(nnv, hmev)
                      binsuffix = "nnbin%dHME%dqtlthreshquadraticpoisson"%(nnv, hmev)
                  elif not HMEqtl and quadraticthresh:
                      #binsuffix = "nnbin%dHME%.2fthreshquadratic"%(nnv, hmev)
                      binsuffix = "nnbin%dHME%.2fthreshquadraticpoisson"%(nnv, hmev)
                  if thresh == 0:##equalNN binnings
                       binsuffix = "nnbin%dHME%.2fNothresh"%(nnv, hmev) 
                  binsuffix = binsuffix.replace('.','p')
                  #newname = "Graviton_syst_allbenchmarks_FR2_rebin_%s"%binsuffix
                  #newfolder = os.path.join(outputfolder, newname)
                  #if not os.path.exists(newfolder):
                  #    os.mkdir(newfolder)
                  #os.system("cp "+folder+"/*txt "+newfolder+"/")
                  filename = "condor_rebin_mass%d_%s.cmd"%(mass, binsuffix)
                  lastjob = hmev == HME_bins[-1] and nnv == nn_bins[-1] and thresh == thresholds[-1]
                  if lastjob:
                      print("creating condor job ",mass, binsuffix)
                  #generate_condorjob(mass, binsuffix, hmev, nnv, thresh, filename, shellname, pyfile, lastjob)
                  genereate_condorjobs_multiqueue(mass, folder, filelist, binsuffix, hmev, nnv, thresh, filename, shellname, pyfile, HMEqtl, quadraticthresh,lastjob)
                  submitscript.write("condor_submit "+filename+" \n")
    submitscript.write("echo 'finishing submitting job "+submitscriptname +" ... '  \n")
    submitscript.close()
    os.system("chmod +x "+submitscriptname)

def checkfiledone(mass, binsuffix):
    eventcat = ["boosted_GGF", "resolved1b_GGF", "resolved2b_GGF"]
    runyears = ["2016", "2017","2018"]
    folder = "/eos/user/t/tahuang/FR2_Florian_v2/datacard_fit_Resonant_HighMass_Graviton_2D_syst_M_%d_raw_FR2"%mass
    if mass <= 500:
        folder = "/eos/user/t/tahuang/FR2_Florian_v2/datacard_fit_Resonant_LowMass_Graviton_2D_syst_M_%d_raw_FR2"%mass
    outfolder = "/eos/user/t/tahuang/HHbbww_resonant_DL_Rebin_v2_202111/Graviton_syst_allbenchmarks_FR2_rebin_{binsuffix}".format(binsuffix=binsuffix)
    filedone = True
    filelist = []
    for year in runyears:
        for cat in eventcat:
            fname = "HH_DL_{mass}_{cat}_{year}".format(mass=mass, cat=cat, year=year)
            infile  = os.path.join(folder, fname +".root")
            outfile = os.path.join(outfolder, fname +"_rebin1d.root")
            #infile = os.path.join(folder, "HH_DL_{mass}_{cat}_{year}.root".format(mass = mass, cat=cat, year=year))
            #outfile = os.path.join(folder, "HH_DL_{mass}_{cat}_{year}_rebin1d.root".format(mass = mass, cat=cat, year=year)) 
            if not os.path.exists(outfile):
                #print("file not found ", outfile)
                filedone = False
                filelist.append(infile)
                continue
            rinfile = ROOT.TFile(infile, "READ")
            keylist_in = rinfile.GetListOfKeys()
            if len(keylist_in) == 0:
                print("Error! no histogram in infile ", infile)
            routfile = ROOT.TFile(outfile, "READ")
            keylist_out = routfile.GetListOfKeys()
            if len(keylist_in)  != len(keylist_out)-1:
                print("numbers of hists old", len(keylist_in)," new ", len(keylist_out), outfile)
                filedone = False
                filelist.append(infile)
            rinfile.Close()
            routfile.Close()
    return filedone,filelist

def checkjob_signode(masspoints, HME_bins, nn_bins, thresholds, HMEqtl, quadraticthresh):
    #shellname = "runrebin_general.sh"
    shellname = "runrebin_condorsig.sh"
    pyfile = "rebin_florian.py"
    submitscriptname = "submit_checkcondorjobs_signode.sh"
    if HMEqtl and not quadraticthresh: 
        #pyfile = "rebin_florian_hmeqtl_fixedthresh.py"
        submitscriptname = "submit_checkcondorjobs_signode_hmeqtl_fixedthresh.sh"
    elif HMEqtl and quadraticthresh:
        #pyfile = "rebin_florian_hmeqtl_quadraticthresh.py"
        submitscriptname = "submit_checkcondorjobs_signode_hmeqtl_quadraticthresh.sh"
    elif not HMEqtl and quadraticthresh:
        #pyfile = "rebin_florian_hmenorm_quadraticthresh.py"
        submitscriptname = "submit_checkcondorjobs_signode_hmenorm_quadraticthresh.sh"
    elif thresholds == [0]:
        #pyfile = "rebin_florian.py"
        submitscriptname = "submit_checkcondorjobs_signode_equalNNbins.sh"
    else:
        #pyfile = "rebin_florian.py"
        submitscriptname = "submit_checkcondorjobs_signode.sh"
    submitscript = open(submitscriptname,"w")
    submitscript.write("echo 'starting to submit job "+submitscriptname +"'\n")
    #submitscript.write("rm condorout/* \n")
    #submitscript.write("rm condorerr/* \n")
    #submitscript.write("rm condorlog/* \n")
    for mass in masspoints:
        folder = "/eos/user/t/tahuang/FR2_Florian_v2/datacard_fit_Resonant_HighMass_Graviton_2D_syst_M_%d_raw_FR2"%mass
        if mass <= 500:
            folder = "/eos/user/t/tahuang/FR2_Florian_v2/datacard_fit_Resonant_LowMass_Graviton_2D_syst_M_%d_raw_FR2"%mass
        nfiles = 0
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
                  if thresh == 0:##equalNN binnings
                       binsuffix = "nnbin%dHME%.2fNothresh"%(nnv, hmev) 
                  binsuffix = binsuffix.replace('.','p')
                  lastjob = hmev == HME_bins[-1] and nnv == nn_bins[-1] and thresh == thresholds[-1]
                  outfilename = "Graviton_syst_allbenchmarks_FR2_rebin_{binsuffix}".format(binsuffix = binsuffix)
                  filename = "condor_rebin_M%d_%s.cmd"%(mass, binsuffix)
                  if not os.path.exists(filename):
                      print("Error!! this condor job is not created before! ", filename)
                  filedone, filelist = checkfiledone(mass, binsuffix)
                  nfiles += len(filelist)
                  if lastjob:
                      print("creating condor job for resubmission ",mass, binsuffix, " total num files ", nfiles)
                  if not filedone:
                      #print("previous job failed ,resubmit it again ",filename)
                      #generate_condorjob(mass, binsuffix, hmev, nnv, thresh, filename, shellname, pyfile, lastjob)
                      genereate_condorjobs_multiqueue(mass, folder, filelist, binsuffix, hmev, nnv, thresh, filename, shellname, pyfile, HMEqtl, quadraticthresh,lastjob)
                      submitscript.write("condor_submit "+filename+" \n")
    submitscript.write("echo 'finishing submitting job "+submitscriptname +" ... '  \n")
    submitscript.close()
    os.system("chmod +x "+submitscriptname)




EqualNNbin = False
#masspoints = [260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 900]
masspoints = [260, 270, 300, 350, 400, 450, 500, 600, 650, 700, 800, 900]
mode = "allsys"
#masspoints = [700]
#HME_bins = [1, 3]
#for i in range(4):
for i in [1, 3]:
    HMEquantile_binning = i/2
    quadraticthreshalgo = i%2
    HME_bins = [1.1, 1.15, 1.2]
    nn_bins = [10,20, 40]
    thresholds = [2, 5, 10]
    if HMEquantile_binning:
        #HME_bins = [5, 10, 15]
        HME_bins = [5,8,10]
    if quadraticthreshalgo:
        thresholds = [10]

    if EqualNNbin:
        thresholds = [0]
        nn_bins = [5, 10, 20]


    if HMEquantile_binning and quadraticthreshalgo:
        print("Running rebin mode with HME quantile binning and quadratic threshold algo, and systematics: ", mode)
    elif HMEquantile_binning and not quadraticthreshalgo:
        print("Running rebin mode with HME quantile binning and fixed threshold algo, and systematics: ", mode)
    elif not HMEquantile_binning and quadraticthreshalgo:
        print("Running rebin mode with customized HME binning and quadratic threshold algo, and systematics: ", mode)
    else:
        print("Running rebin mode with customized HME binning and fixed threshold algo, and systematics: ", mode)

    print("NN binning options ", nn_bins, " HME binning options ", HME_bins," thresholds ", thresholds)

    #generatejob_bkgnode(masspoints)
    generatejob_signode(masspoints, HME_bins, nn_bins, thresholds, HMEquantile_binning, quadraticthreshalgo)
    #checkjob_signode(masspoints, HME_bins, nn_bins, thresholds, HMEquantile_binning, quadraticthreshalgo)
