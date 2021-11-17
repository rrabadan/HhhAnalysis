import os
import ROOT
import re
import numpy as np
from datetime import datetime
import sys 
import csv
sys.argv.append( '-b' )
#or ROOT.gROOT.SetBatch(1)

import limithelper

mainSys = ["scale_mu","PS","pile"]
eventcat = ["boosted_DY_VVV", "boosted_GGF", "boosted_other", "resolved1b_GGF", "resolved2b_GGF","resolved_DY_VVV","resolved_other"]
eventcat = ["boosted_GGF", "resolved1b_GGF", "resolved2b_GGF"]
runyears = ["2016", "2017","2018"]

folder = "/afs/cern.ch/work/t/tahuang/CombinedLimit/ForGithub/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/Rebin_Florian/"

def get_florian_limits(spin, masspoints, scriptsuffix):
    folder_florian = "Reproduce_Florian_Run2_DL_20211108/"
    limits_dict = {}
    allmassfound = True
    for mass in masspoints:
        ##spin0_900_FR_s1fb_tm1_autoMC10_s1fb_tm1_autoMC10_limits.log
        logfile = os.path.join(folder_florian,
        "spin{spin}_{mass}_FR_{scriptsuffix}_{scriptsuffix}_limits.log".format(mass=mass, spin=spin,
        scriptsuffix=scriptsuffix))
        if not os.path.exists(logfile):
            print("Error!! no logfile is found ", logfile)
            allmassfound = False
            continue
        limits_dict[mass] = limithelper.extractlimitfromtxtfile(logfile)
        if len(limits_dict[mass].keys()) != 6:
            print("Warning! not all limits for this mass found ", logfile, limits_dict[mass])
            allmassfound = False

    #print("florian limits ", limits_dict)
    return limits_dict

def get_limits_onebinning(masspoints, binsuffix, scriptsuffix, mode, plotdir, makeplot):
    """ Return limits for one binning case, dic[Runyear][cat] = limits
    HH_260_resolved_other_2018_rebin1d_lnN_s1fb_tm1_autoMC10_limits.log
    """
    xtitle = "Graviton mass [GeV]"
    runyears_all = runyears+["FR"]
    #eventcat_all = eventcat+["allcat"]
    eventcat_all = ["allcat"]
    limits_dict = {}
    #thisfolder = "/eos/user/t/tahuang/HHbbww_resonant_DL_Rebin_v2_202111/Graviton_syst_allbenchmarks_FR2_rebin_%s"%binsuffix
    thisfolder = "Graviton_syst_allbenchmarks_FR2_rebin_%s"%binsuffix
    for year in runyears_all:
        limits_dict[year] = {}
        for cat in eventcat_all:
            if year == "FR" and cat != "allcat":
                continue
            limits_dict[year][cat] = {}
            allmassfound = True
            for mass in masspoints:
                #if "HH_{mass}_{cat}_{year}".format(mass = mass, cat=cat, year=year) == "HH_800_boosted_GGF_2016":
                #   ## broken file
                #   continue 
                filename = "HH_{mass}_{cat}_{year}_rebin1d_{mode}_{scriptsuffix}_limits.log".format(mass = mass, cat=cat, year=year, mode=mode, scriptsuffix=scriptsuffix)
                #thisfolder = "datacard_fit_Resonant_LowMass_Graviton_2D_syst_M_{mass}_raw_FR2_rebin_{binsuffix}".format(binsuffix=binsuffix, mass=mass)
                #if mass >= 550:
                #    thisfolder = "datacard_fit_Resonant_HighMass_Graviton_2D_syst_M_{mass}_raw_FR2_rebin_{binsuffix}".format(binsuffix=binsuffix, mass=mass)
                logfile = os.path.join(thisfolder, filename)
                if not os.path.exists(logfile):
                    print("Error!! no logfile is found ", logfile)
                    allmassfound = False
                    continue
                limits_dict[year][cat][mass] = limithelper.extractlimitfromtxtfile(logfile)
                if len(limits_dict[year][cat][mass].keys()) != 6:
                    print("Warning! not all limits for this mass found ", logfile, limits_dict[year][cat][mass])
                    allmassfound = False
            text = "Run"+year+"_"+cat+"_"+binsuffix
            plotname = os.path.join(plotdir, text+"_"+binsuffix+"_"+scriptsuffix+"_"+mode)
            #if makeplot or (year == "FR" and cat == "allcat"):
            if makeplot:
                limithelper.makeBrazilPlot(masspoints, limits_dict[year][cat], xtitle, text, plotname, False, year)

    return limits_dict 

def get_binsuffix(nnv, hmev, thresh, HMEqtl, quadraticthresh):
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
    return binsuffix

def get_limits_allbinnings(masspoints, nnbins_v, hmebins_v, thresholds, scriptsuffix, mode, compareplotdir, HMEqtl, quadraticthresh):
    
    EqualNNbins = thresholds == [0]
    addFlorianCompare = True
    florian_text = "Florian_FR_allcat_allsys"
    florian_limits = {}
    if addFlorianCompare:
        florian_limits = get_florian_limits(2, masspoints, scriptsuffix)

    xtitle = "Graviton mass [GeV]"
    makesingleplot = False
    runyears_all = runyears + ["FR"]
    eventcat_all = eventcat+["allcat"]
    ##all_limits_dict [binsuffix][Runyear][cat]
    all_limits_dict = {}
    for nnbin in nnbins_v:
        for hmebin in hmebins_v:
            for thresh in thresholds:
                binsuffix = get_binsuffix(nnbin, hmebin, thresh, HMEqtl, quadraticthresh)
                #print("hmebin  ", hmebin, " binsuffix ", binsuffix)
                lastjob = hmebin == HME_bins[-1] and nnbin == nn_bins[-1] and thresh == thresholds[-1]
                if lastjob:
                    print("creating condor job ",binsuffix)
                #plotdir = "limitsplots_graviton_{binsuffix}_FR2".format(binsuffix=binsuffix)
                #if not os.path.exists(plotdir):
                #    os.mkdir(plotdir)
                all_limits_dict[binsuffix] = get_limits_onebinning(masspoints, binsuffix, scriptsuffix, mode, compareplotdir, makesingleplot)

    #nominalbin = "nnbin10HMEv1"
    nominalbin = "nnbin20HME1p15thresh5"
    if HMEqtl and not quadraticthresh:
        nominalbin = "nnbin20HME10qtlthresh5"
    elif HMEqtl and quadraticthresh:
        nominalbin = "nnbin20HME10qtlthreshquadratic"
    elif not HMEqtl and quadraticthresh:
        nominalbin = "nnbin20HME1p15threshquadratic"
    if EqualNNbin:
        nominalbin = "nnbin20HME1p15Nothresh"

    ##compare within one year, different categories
    plot_catergories = False
    reference = "allcat"
    cats_to_plot = ["allcat","boosted_DY_VVV", "boosted_other", "resolved1b_GGF",
    "resolved2b_GGF","resolved_DY_VVV","resolved_other"]
    if plot_catergories:
        for year in runyears:
            cat_to_limits = {}
            for cat in cats_to_plot:
                cat_to_limits[cat] = all_limits_dict[nominalbin][year][cat]
            if addFlorianCompare and year == "FR":
                reference = florian_text
                cat_to_limits["Florian_FR_allcat_allsys"] = florian_limits
            text="Evtcategories_"+year+"_"+nominalbin
            plotname = os.path.join(compareplotdir, text)
            reference = florian_text
            limithelper.makeComparePlot(masspoints, cat_to_limits, reference, xtitle, text, plotname, False, year)

    ##compare different years
    plot_years = True
    #cats_plot_years = ["allcat", "resolved1b_GGF", "resolved2b_GGF"]
    cats_plot_years = ["allcat"]
    years_to_plot = ["FR"] + runyears
    reference = "FR"
    if plot_years:
        for cat in cats_plot_years:
            years_to_limits = {}
            for year in years_to_plot:
                years_to_limits[year] = all_limits_dict[nominalbin][year][cat]
            if addFlorianCompare and cat == "allcat":
                reference = florian_text
                years_to_limits["Florian_FR_allcat_allsys"] = florian_limits
            text="Runyears_"+cat+"_"+nominalbin
            plotname = os.path.join(compareplotdir, text)
            limithelper.makeComparePlot(masspoints, years_to_limits, reference, xtitle, text, plotname, False, "Full Run2")


    ##compare different HME binnings
    years_bincompare = ["2016", "FR"]
    cats_bincompare = ["allcat"]
    nnbins_bincompare = [10, 20]
    hmebins_bincompare = [1.1, 1.15, 1.2]
    threshs_bincompare = [2, 5, 10]
    if HMEqtl:
        hmebins_bincompare = [5, 10]
    
    plot_HMTbinnings = True
    reference = "HMEquantilebin10" if HMEqtl else "HME1p15" 
    if plot_HMTbinnings:
        for year in years_bincompare:
            for cat in cats_bincompare:
                for nnbin in nnbins_bincompare:
                    for thresh in thresholds:
                        hmebins_to_limits= {}
                        for hmebin in hmebins_v:
                            binsuffix = get_binsuffix(nnbin, hmebin, thresh, HMEqtl, quadraticthresh)
                            hmebinsuffix = "HMEquantilebin%d"%hmebin if HMEqtl else ("HME%.2f"%hmebin).replace(".","p")
                            hmebins_to_limits[hmebinsuffix] = all_limits_dict[binsuffix][year][cat]
                        if addFlorianCompare and year == "FR" and cat == "allcat":
                            reference = florian_text
                            hmebins_to_limits["Florian_FR_allcat_allsys"] = florian_limits
                        text="HMEquantilebins_" if HMEqtl else "HMEbins_"
                        text += year+"_"+cat+"_nnbin%d"%nnbin
                        if thresh > 0:
                            text += "_quadraticthresh%d"%thresh if quadraticthresh else "_thresh%d"%thresh
                        else:
                            text += "_equalNNbins"
                        plotname = os.path.join(compareplotdir, text)
                        limithelper.makeComparePlot(masspoints, hmebins_to_limits, reference, xtitle, text, plotname, False, year)

    plot_nnbinnings = True
    reference = "nnbin10"
    if plot_nnbinnings:
        for year in years_bincompare:
            for cat in cats_bincompare:
                for hmebin in hmebins_bincompare:
                    for thresh in thresholds:
                        nnbins_to_limits= {}
                        for nnbin in nnbins_v:
                            binsuffix = get_binsuffix(nnbin, hmebin, thresh, HMEqtl, quadraticthresh)
                            nnbins_to_limits["nnbin%d"%nnbin] = all_limits_dict[binsuffix][year][cat]
                        if addFlorianCompare and year == "FR" and cat == "allcat":
                            reference = florian_text
                            nnbins_to_limits["Florian_FR_allcat_allsys"] = florian_limits
                        text="NNbins_"+year+"_"+cat
                        text += "_hmeqtlbin%d"%hmebin if HMEqtl else ("HME%.2f"%hmebin).replace(".","p") 
                        if thresh > 0:
                            text += "_quadraticthresh%d"%thresh if quadraticthresh else "_thresh%d"%thresh
                        else:
                            text += "_equalNNbins"
                        plotname = os.path.join(compareplotdir, text)
                        limithelper.makeComparePlot(masspoints, nnbins_to_limits, reference, xtitle, text, plotname, False, year)

    plot_thresholds = not EqualNNbins
    reference = ""
    if plot_nnbinnings:
        for year in years_bincompare:
            for cat in cats_bincompare:
                for hmebin in hmebins_bincompare:
                    for nnbin in nnbins_bincompare:
                        threshs_to_limits= {}
                        for thresh in thresholds:
                            binsuffix = get_binsuffix(nnbin, hmebin, thresh, HMEqtl, quadraticthresh)
                            threshbinsuffix = "NNbin%dthresh%d"%(nnbin, thresh)
                            threshs_to_limits[threshbinsuffix] = all_limits_dict[binsuffix][year][cat]
                            reference = threshbinsuffix
                        if addFlorianCompare and year == "FR" and cat == "allcat":
                            reference = florian_text
                            threshs_to_limits["Florian_FR_allcat_allsys"] = florian_limits
                        text="ThreshAlgo_"+year+"_"+cat+"_nnbin%d"%nnbin
                        text += "_hmeqtlbin%d"%hmebin if HMEqtl else ("_HME%.2f"%hmebin).replace(".","p") 
                        plotname = os.path.join(compareplotdir, text)
                        limithelper.makeComparePlot(masspoints, threshs_to_limits, reference, xtitle, text, plotname, False, year)
		        
		        

#masspoints = [260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 900]
masspoints = [260, 270, 300, 350, 400, 450, 500, 600, 650, 700, 800, 900]
nn_bins = [10, 20, 40]
HME_bins = [1.1, 1.15, 1.2]
thresholds = [2, 5, 10]
dateTimeObj = datetime.now()


scriptsuffix = "s1fb_tm1_autoMC10" 
#mode = "mainSys"
mode = "lnN"
#mode = "allsys"
compareplotdir = "compare_limitsplots_FR2_rebin_%s_"%mode+dateTimeObj.strftime("%d-%b-%Y")
EqualNNbin = False
irange = 4
if EqualNNbin:
    irange =1 
    compareplotdir += "_equalNNbins"
if not os.path.exists(compareplotdir):
    os.makedirs(compareplotdir)

#get_florian_limits(2, masspoints, scriptsuffix)
#get_limits_allbinnings(masspoints, nnbins_v, hmebins_v, thresholds, scriptsuffix, mode, compareplotdir, threshalgo)

for i in range(irange):
#for i in [2, 3]:
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
    if EqualNNbin:
        thresholds = [0]
        nn_bins = [5, 10,20]

    if HMEquantile_binning and quadraticthreshalgo:
        print("Running rebin mode with HME quantile binning and quadratic threshold algo, and systematics: ", mode)
    elif HMEquantile_binning and not quadraticthreshalgo:
        print("Running rebin mode with HME quantile binning and fixed threshold algo, and systematics: ", mode)
    elif not HMEquantile_binning and quadraticthreshalgo:
        print("Running rebin mode with customized HME binning and quadratic threshold algo, and systematics: ", mode)
    else:
        print("Running rebin mode with customized HME binning and fixed threshold algo, and systematics: ", mode)
    print("NN binning options ", nn_bins, " HME binning options ", HME_bins," thresholds ", thresholds)

    get_limits_allbinnings(masspoints, nn_bins, HME_bins, thresholds, scriptsuffix, mode, compareplotdir, HMEquantile_binning, quadraticthreshalgo)

