import os
import sys



finalpath= "/eos/user/t/tahuang/FR2_Florian_v1_rebin/"
intpath = "/afs/cern.ch/work/t/tahuang/CombinedLimit/ForGithub/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/Rebin_Florian"
masspoints = [260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 900]
nnbins_v = [5, 8, 10, 20]
hmbins_v = [0,1,2,3]

#masspoints = [ 300, 400, 450, 500, 550, 600, 650, 700, 900]
#masspoints = [900]

#nnbins_v = [10]
#hmbins_v = [3]

scriptsuffix = "s1fb_tm1_autoMC10" 
#mode = "mainSys"
mode = "lnN"
#mode = "allsys"
compareplotdir = "compare_limitsplots_allsys_FR2_rebin_2021106"
bkg_qltbinning = True
qtlbin = False

for nnbin in nnbins_v:
    for hmebin in hmbins_v:
        binsuffix = "nnbin%sHMEv%d"%(nnbin, hmebin)
        if qtlbin:
            binsuffix = binsuffix+"qltbin"
        eosfolder_name = "Graviton_syst_allbenchmarks_FR2_rebin_{binsuffix}".format(binsuffix=binsuffix)
        finalfolder = os.path.join(finalpath, eosfolder_name) 
        if not os.path.exists(finalfolder):
            os.mkdir(finalfolder)
        folder_results = os.path.join(intpath, eosfolder_name)
        if not os.path.exists(folder_results):
            os.mkdir(folder_results)
        for mass in masspoints:
            thisfolder = "datacard_fit_Resonant_LowMass_Graviton_2D_syst_M_{mass}_raw_FR2_rebin_{binsuffix}".format(binsuffix=binsuffix, mass=mass)
            if mass >= 550:
                thisfolder = "datacard_fit_Resonant_HighMass_Graviton_2D_syst_M_{mass}_raw_FR2_rebin_{binsuffix}".format(binsuffix=binsuffix, mass=mass)
            inputfolder = os.path.join(intpath, thisfolder)
            os.system("cp -f " + inputfolder + "/*.root " + finalfolder )
            os.system("cp -f " + inputfolder + "/*.txt " + finalfolder)
            os.system("cp -f " + inputfolder + "/*.log " + folder_results)
            os.system("cp -f " + inputfolder + "/*.sh " + folder_results)




