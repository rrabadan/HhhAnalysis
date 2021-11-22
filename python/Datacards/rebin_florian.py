import os
from array import array
import datetime
from math import sqrt
import numpy.lib.recfunctions as recfunctions
import numpy as np
import ROOT

import argparse
import sys
#sys.argv.append( '-b' )
#sys.argv.append( '-q' )

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

mainSys = ["scale_mu","PS","pile"]
eventcat = ["boosted_DY_VVV", "boosted_GGF", "boosted_other", "resolved1b_GGF", "resolved2b_GGF","resolved_DY_VVV","resolved_other"]
runyears = ["2016", "2017","2018"]

##========================================================
## origin hist: x-axis: 0-1.0, 400bins, 0-2000.0, 400bins
##========================================================
def print_xybins(hist):
    """
    print out the 2d histogram binning
    """
    hist_nbin_x = hist.GetNbinsX()
    hist_nbin_y = hist.GetNbinsY()
    xaxis = hist.GetXaxis() 
    yaxis = hist.GetYaxis() 
    for i in range(1,hist_nbin_x):
        print("xbin ", i, " lowEdge ", xaxis.GetBinLowEdge(i))
    print("xbin ", hist_nbin_x, " UpEdge ", xaxis.GetBinUpEdge(hist_nbin_x))
    for i in range(1,hist_nbin_y):
        print("ybin ", i, " lowEdge ", yaxis.GetBinLowEdge(i))
    print("ybin ", i, " UpEdge ", yaxis.GetBinUpEdge(hist_nbin_y))

    for i in range(1,hist_nbin_x+1):
        for j in range(1,hist_nbin_y+1):
            print("binx ",i," biny ",j, " bincontent ", hist.GetBinContent(i,j), " binerror ", hist.GetBinError(i,j))

def generateHMEbins(mass, version=0):
    """
    generate HME binning with different version number
    """
    lowM = 250.0; highM = 1250.0
    xbins = [lowM]
    x = lowM
    step = 25.0 
    if mass >= 400 and mass < 500:
        step = 30.0
    elif mass>=500 and mass < 700:
        step = 40.0
    elif mass>=700:
        step = 50.0
    step = step + 10*version
    gap1 = 50 + mass*.1 #120.0
    gap2 = 50 + mass*.2 #200.0
    gap3 = 60 + mass*.3 #
    while x <= highM-150.0:
        if abs(x-mass) > mass:
            x = x+100.0
            xbins.append(x)
        elif abs(x-mass)<= gap2 and abs(x-mass)>gap1:
            x = x+step*1.6
            xbins.append(x)
        elif abs(x-mass)> gap3:
            x = x+step*3
            xbins.append(x)
        #elif abs(x-mass)> 200.:
        #    x = x+step*2.5
        #    xbins.append(x)
        else:
            x = x+step
            xbins.append(x)
    if mass<800:
        xbins = xbins[:-1]
    xbins.append(highM)
    xbins[0]  = 0.0
    xbins.append(2000.0)
    print("Benchmark ", mass , " v ", version, " HME mass bins ",xbins)
    return np.asarray(xbins)
    #return xbins

def generateHMEbinsv2(mass):
    lowM = 250.0; highM = 1200.0
    #if mass > 700: ## from
    ## within 100 GeV, use 50 Gev binsize, otherwise use 100 GeV binsize
    xbins = [lowM]
    step = 50
    x = lowM
    while x < highM-100.0:
        if abs(x-mass)<100.0:
            x = xbins[-1]+40.0
        elif abs(x-mass)<250.0:
            x = xbins[-1]+70.0
        elif abs(x-mass)<400.0:
            x = xbins[-1]+100.0
        else:
            x = xbins[-1]+150.0
        xbins.append(x)

    xbins[0] = 0.0
    xbins.append(2000.0)
    print("Benchmark ", mass ," HME mass bins ",xbins)
    return np.asarray(xbins)
	
def generateHMEbinsv3(mass):
    lowM = 250.0; highM = 1200.0
    #if mass > 700: ## from
    ## within 100 GeV, use 50 Gev binsize, otherwise use 100 GeV binsize
    xbins = [lowM]
    step = 50
    x = lowM
    while x < 1150.0:
        x = xbins[-1]+90.0
        xbins.append(x)

    xbins[0]  = 0.0
    xbins.append(2000.0)
    print("Benchmark ", mass ," HME mass bins ",xbins)
    return np.asarray(xbins)
	
def generateHMEbinsv4(a):
    """
    a must be > 1.0. : 1.1, 1.15, 1.2, 1.25
    """
    lowM = 250.0; highM = 1200.0
    x = lowM
    xbins = [lowM]
    while x <= highM-100.0:
      x = int(xbins[-1]*a)
      xbins.append(x)
    xbins[0]  = 0.0
    xbins.append(2000.0)
    print("generateHMEbinsv4: ",a,xbins)
    return np.asarray(xbins)


def getBinsArray(hist, ibin=0):
    """Return the numpy array of histogram binning
    ibin = 0 => return xaxis binning
    ibin = 1 => return yaxis binning
    """
    bins = []
    hist_nbin = None
    xaxis = None
    if ibin == 0:
        hist_nbin = hist.GetNbinsX()
        xaxis = hist.GetXaxis()
    else:
        hist_nbin = hist.GetNbinsY()
        xaxis = hist.GetYaxis()

    for i in range(1, hist_nbin+1):
        bins.append(xaxis.GetBinLowEdge(i))
    bins.append(xaxis.GetBinUpEdge(hist_nbin))
    return np.asarray(bins)

def find_quantile_binning(hist, n):
    """Return the quantile binning for 1D histogram
    """
    nbin = hist.GetNbinsX()
    fraction = 1.0/n
    totalint = hist.Integral()
    #print("hist integral ", totalint, " number of final bins ", n," quantile fraction ", fraction)
    xbin_list = [hist.GetXaxis().GetBinUpEdge(nbin)]
    i  = nbin
    tot_bincontent = 0.0
    while i >= 1:
        tot_bincontent += hist.GetBinContent(i)
        if tot_bincontent >= totalint*fraction:
            binedge = hist.GetXaxis().GetBinLowEdge(i)
            xbin_list.insert(0, binedge)
            #print("found new binnning ", binedge, " bincontent ", tot_bincontent)
            tot_bincontent = 0.0
        i -= 1

    binedge = hist.GetXaxis().GetBinLowEdge(1)
    if tot_bincontent < totalint*fraction/2.0:
        xbin_list[0] = binedge
    else:
        xbin_list.insert(0, binedge)
    #print("new binning ", xbin_list)
    return np.asarray(xbin_list)

def find_quantile_HMEbinning(hist, n):
    """Return the quantile binning for 1D histogram
    """
    projHistY = hist.ProjectionY()
    return find_quantile_binning(projHistY, n)


def find_threshold_binning_linearized2d(totalhist, mainbkg_list,  nbin_y, thresh, use_quadratic_thresh=False):
    """Return the quantile binning for linearized 2D histogram
    """
    nbin = int(totalhist.GetNbinsX()/nbin_y)
    #print("totalhist.GetNbinsX() ",totalhist.GetNbinsX(), " nbiny ", nbin_y, " totalhist.GetNbinsX()/nbin_y ",totalhist.GetNbinsX()/nbin_y)
    thresh_array = None
    if use_quadratic_thresh:
        thresh_array = thresh + np.arange(0, nbin) ** 2
        
    axis = totalhist.GetXaxis()
    mainbkg_sum = [0.0]*len(mainbkg_list)
    xbin_list = []
    newbin_index = 0
    for i in range(nbin_y):
        #upbound = axis.GetBinUpEdge()
        upbin  = (nbin_y-i)*nbin
        lowbin = (nbin_y-i-1)*nbin+1
        xbin_list.insert(0, axis.GetBinUpEdge(upbin))
        j = upbin
        mainbkg_sum = [0.0]*len(mainbkg_list)
        bincontent_sum = 0.0
        binerr2_sum = 0.0
        newbin_index = 0
        #print("i ",i, " lowbin ",lowbin, " upbin ",upbin, " current xbinlist ", xbin_list)
        while j >= lowbin:
            for ih,h in enumerate(mainbkg_list):
                mainbkg_sum[ih] += h.GetBinContent(j)
            bincontent_sum += totalhist.GetBinContent(j)
            binerr2_sum  += totalhist.GetBinError(j)**2
            nonzero = all(x>0.0 for x in mainbkg_sum)
            thisthresh = thresh_array[newbin_index] if use_quadratic_thresh else  thresh
            if bincontent_sum-sqrt(binerr2_sum) >= thisthresh and nonzero:
                if axis.GetBinUpEdge(j) not in xbin_list:
                    xbin_list.insert(0, axis.GetBinUpEdge(j))
                    #print("inserting new binedge ", axis.GetBinUpEdge(j))
                #print("biny i ",i," newbin_index ", newbin_index, " thresh ", thresh, xbin_list)
                mainbkg_sum = [0.0]*len(mainbkg_list)
                bincontent_sum = 0.0
                binerr2_sum = 0.0
                newbin_index += 1
            j -= 1

    xbin_list.insert(0, axis.GetBinLowEdge(1))
    #print("old threshold algo for quantile binning for linearized 2D ", xbin_list)
    return np.asarray(xbin_list)

def find_threshold_binning_linearized2d_new(totalhist, mainbkg_list,  nbin_y, thresh, use_quadratic_thresh=False):
    """Return the quantile binning for linearized 2D histogram
    reserved for new threshold algorithm
    """
    nbin = int(totalhist.GetNbinsX()/nbin_y)
    #print("totalhist.GetNbinsX() ",totalhist.GetNbinsX(), " nbiny ", nbin_y, " totalhist.GetNbinsX()/nbin_y ",totalhist.GetNbinsX()/nbin_y)
    thresh_array = np.arange(1, nbin) ** 2
    print("thresh array ", thresh_array)
        
    axis = totalhist.GetXaxis()
    mainbkg_sum = [0.0]*len(mainbkg_list)
    xbin_list = []
    newbin_index = 0
    poissonerror = 1.86
    for i in range(nbin_y):
        #upbound = axis.GetBinUpEdge()
        upbin  = (nbin_y-i)*nbin
        lowbin = (nbin_y-i-1)*nbin+1
        xbin_list.insert(0, axis.GetBinUpEdge(upbin))
        j = upbin
        mainbkg_sum = [0.0]*len(mainbkg_list)
        bincontent_sum = 0.0
        binerr2_sum = 0.0
        newbin_index = 0
        #print("i ",i, " lowbin ",lowbin, " upbin ",upbin, " current xbinlist ", xbin_list)
        while j >= lowbin:
            for ih,h in enumerate(mainbkg_list):
                mainbkg_sum[ih] += h.GetBinContent(j)
            if totalhist.GetBinContent(j) > 0.0:
                bincontent_sum += totalhist.GetBinContent(j)
                binerr2_sum  += totalhist.GetBinError(j)**2
            else:
                #bincontent_sum += totalhist.GetBinContent(j)
                binerr2_sum  += poissonerror**2
            nonzero = all(x>0.0 for x in mainbkg_sum)
            thisthresh = thresh_array[newbin_index] if use_quadratic_thresh else  thresh
            if bincontent_sum-sqrt(binerr2_sum) >= thisthresh and nonzero:
                if axis.GetBinUpEdge(j) not in xbin_list:
                    xbin_list.insert(0, axis.GetBinUpEdge(j))
                    #print("inserting new binedge ", axis.GetBinUpEdge(j))
                #print("biny i ",i," newbin_index ", newbin_index, " thresh ", thresh, xbin_list)
                mainbkg_sum = [0.0]*len(mainbkg_list)
                bincontent_sum = 0.0
                binerr2_sum = 0.0
                newbin_index += 1
            j -= 1

    xbin_list.insert(0, axis.GetBinLowEdge(1))
    #print("New threshold algo for quantile binning for linearized 2D ", xbin_list)
    return np.asarray(xbin_list)


def rebinx_2d(hist, xbins):
    """Return the x-axis rebinned histogram for 2D histogram
    """
    ybins = getBinsArray(hist, 1)
    #print("rebinx, keep ybins ", ybins)
    hist_nbin_x = hist.GetNbinsX()
    if np.array_equal(xbins, getBinsArray(hist, 0)):
        return hist
    #print("rebinx ", len(xbins), " getNbinsX ", hist_nbin_x)
    xaxis = hist.GetXaxis() 
    xmin = xaxis.GetBinLowEdge(1)
    xmax = xaxis.GetBinUpEdge(hist_nbin_x)
    if xmin >= xbins[-1] or xmax <= xbins[0]:
        print("ignore this rebinx since xmin is larger than max of new bins or xmax is smaller than min of new bins")
        return hist
    if xmin  < xbins[0]:
        xbins[0] =  xmin
    if xmax  > xbins[-1]:
        xbins[-1] =  xmax

    newhist = ROOT.TH2F(hist.GetName()+"_rebinx", hist.GetTitle(),  len(xbins)-1, xbins, len(ybins)-1, ybins)
    hist_nbin_y = hist.GetNbinsY()
    klow  = 1
    for i, xlow in enumerate(xbins[:-1]):
        xup = xbins[i+1]
        #print("klow ", klow, " klowEdge ", xaxis.GetBinLowEdge(klow), "new binx ", i, " xlow ", xlow, " xup ", xup)
        for j in range(1, hist_nbin_y+1):
            bincontent = 0.0
            binerr2  = 0.0
            k = klow
            while xaxis.GetBinUpEdge(k) <= xup:
                bincontent +=  hist.GetBinContent(k, j)
                binerr2  += hist.GetBinError(k, j)*hist.GetBinError(k, j)
                k += 1
                #print("mering ", k," for xup ",xup)
            newhist.SetBinContent(i+1, j, bincontent)
            newhist.SetBinError(i+1, j, sqrt(binerr2))
        klow = k

    newhist.SetName(hist.GetName())
    return newhist


def rebiny_2d(hist, ybins):
    """Return the y-axis rebinned histogram for 2D histogram
    """
    xbins = getBinsArray(hist, 0)
    #print("rebiny, keep xbins ", xbins)
    if np.array_equal(ybins, getBinsArray(hist, 1)):
        return hist
    hist_nbin_y = hist.GetNbinsY()
    xaxis = hist.GetYaxis() 
    ymin = xaxis.GetBinLowEdge(1)
    ymax = xaxis.GetBinUpEdge(hist_nbin_y)
    #print("old bins ymin ",ymin, " ymax ",ymax, " new bins ymin ", ybins[0]," ymax ", ybins[-1])
    if ymin >= ybins[-1] or ymax <= ybins[0]:
        print("ignore this rebiny since xmin is larger than max of new bins or xmax is smaller than min of new bins")
        return hist
    if ymin  < ybins[0]:
        ybins[0] =  ymin
    if ymax  > ybins[-1]:
        ybins[-1] =  ymax

    #print("new ybins used ",ybins)
    newhist = ROOT.TH2F(hist.GetName()+"_rebiny", hist.GetTitle(),  len(xbins)-1, xbins, len(ybins)-1, ybins)
    hist_nbin_x = hist.GetNbinsX()
    klow  = 1
    for i, ylow in enumerate(ybins[:-1]):
        yup = ybins[i+1]
        #print("klow ", klow, " klowEdge ", xaxis.GetBinLowEdge(klow), "new biny ", i, " ylow ", ylow, " yup ", yup)
        for j in range(1, hist_nbin_x+1):
            #print("new biny ", i, " ylow ", xlow, " yup ", xup," binx ", j)
            bincontent = 0.0
            binerr2  = 0.0
            k = klow
            while xaxis.GetBinUpEdge(k) <= yup:
                bincontent +=  hist.GetBinContent(j, k)
                binerr2  += hist.GetBinError(j,k)*hist.GetBinError(j,k)
                k += 1
                #print("mering ", k," for yup ",xup)
            newhist.SetBinContent(j, i+1, bincontent)
            newhist.SetBinError(j, i+1, sqrt(binerr2))
        klow = k

    newhist.SetName(hist.GetName())
    return newhist

def linearized2D(hist2d):
    """Return the linearized 2D histogram
    """
    binX = hist2d.GetNbinsX()
    binY = hist2d.GetNbinsY()
    xlow = hist2d.GetXaxis().GetBinLowEdge(1)
    xhigh = hist2d.GetXaxis().GetBinUpEdge(binX)
    ##need to update this to unequal binning
    xbins = getBinsArray(hist2d, ibin=0) 
    ybins = getBinsArray(hist2d, ibin=1) 
    #print("2d histogram, xbins ",xbins, " yins ",ybins)
    newbins = [xbins[0]]
    for k in range(binY):
        newbins += [x+k for x in xbins[1:]]
    newbins = np.asarray(newbins) 
    #print("new 1d bins ",len(newbins)-1," expected ", binX*binY," ",newbins)
    #hist1d = ROOT.TH1F(histname, hist2d.GetTitle(), binX*binY, xlow, xlow+(xhigh - xlow)*binY)
    hist1d = ROOT.TH1F(hist2d.GetName(), hist2d.GetTitle(), binX*binY, newbins)
    for i in range(binX):
        for j in range(binY):
            bincontent = hist2d.GetBinContent(i+1, j+1)
            binerr  = hist2d.GetBinError(i+1, j+1)
            hist1d.SetBinContent(i+1+j*binX, bincontent)
            hist1d.SetBinError(i+1+j*binX, binerr)
    hist1d.SetDirectory(0)
    return hist1d


def rebin_hist2d(hist, xbins, ybins):
    """Return the both x and y rebinned and linearized 2D histogram
    """
    old_integral = hist.Integral()
    #print_xybins(hist)
    #oldxbins = getBinsArray(hist, ibin=0) 
    #oldybins = getBinsArray(hist, ibin=1) 
    #print("Old xbins", oldxbins, " yins", oldybins)
    #print("New xbins", xbins, " yins", ybins)
    if not isinstance(ybins, type(None)):
        hist = rebinx_2d(hist, xbins)
    if not isinstance(xbins, type(None)):
        hist = rebiny_2d(hist, ybins)

    hist = linearized2D(hist) 
    new_integral = hist.Integral()

    if old_integral > 0 and abs(old_integral - new_integral) > 0.01 * max(old_integral, 1.0):
        print("Warning!! integral of rebinned histogram %s has changed, old integral "%histname, old_integral, " new ", new_integral)
    #print_xybins(hist)
    return hist


def rebin_hist1d(hist, xbins):
    """Return the rebinned 1D histogram
    """
    histname = hist.GetName()
    old_integral = hist.Integral()
    xaxis = hist.GetXaxis()
    oldxmin = xaxis.GetBinLowEdge(1)
    oldxmax = xaxis.GetBinUpEdge(hist.GetNbinsX())
    if len(xbins) == 3:
        nbins = xbins[0]; xmin = xbins[1]; xmax =  xbins[2]
        xbins = []
        binwidth = (xmax-xmin)*1.0/nbins
        for i in range(0, nbins+1):
            xbins.append(xmin + i*binwidth)
        xbins = np.asarray(xbins)
    if oldxmin  < xbins[0]:
        xbins[0] =  oldxmin
    if oldxmax  > xbins[-1]:
        xbins[-1] = oldxmax

    newhist = ROOT.TH1F(hist.GetName()+"_rebinx", hist.GetTitle(),  len(xbins)-1, xbins)
    hist_nbin_y = hist.GetNbinsY()
    k = 1
    for i, xlow in enumerate(xbins[:-1]):
        xup = xbins[i+1]
        bincontent = 0.0; binerr2 = 0.0
        while xaxis.GetBinUpEdge(k) <= xup:
            bincontent +=  hist.GetBinContent(k)
            binerr2  += hist.GetBinError(k)*hist.GetBinError(k)
            k += 1
        newhist.SetBinContent(i+1, bincontent)
        newhist.SetBinError(i+1, sqrt(binerr2))
    newhist.SetName(histname)
    new_integral = newhist.Integral()
    if old_integral > 0 and abs(old_integral - new_integral) > 0.01 * max(old_integral, 1.0):
        print("Warning!! integral of rebinned histogram %s has changed, old integral "%histname, old_integral, " new ", new_integral)
    return newhist



def rebin_onefile(rootfile, newfile, xbins, ybins, mode=""):
    rfile = ROOT.TFile(rootfile,"READ") 
    keylist = rfile.GetListOfKeys()
    #print("keylist ",keylist)
    if len(keylist) == 0:
        print("no histogram in file ", rootfile)
        return 0

    ihist = 0; maxhist = -1
    #histname_list = ["data_obs"]
    histlist = []
    hist_totexpect = None
    for key in keylist:
        if ihist >= maxhist and maxhist > 0:
            print("only rebins %d histograms now"%maxhist)
            break
        #isbkg_nominal = "__" not in key.ReadObj().GetName() and "data_obs" not in key.ReadObj().GetName() and "ggHH_" not in key.ReadObj().GetName()  
        isbkg_nominal = "__" not in key.ReadObj().GetName() and "data_obs" not in key.ReadObj().GetName() and "signal_ggf" not in key.ReadObj().GetName()  
        if not hist_totexpect and isbkg_nominal:
            hist_totexpect = key.ReadObj().Clone()
        elif isbkg_nominal:
            hist_totexpect.Add(key.ReadObj())

        if mode == "lnN" and "__" in key.ReadObj().GetName():
            continue
        mainSys_true = any(sysname in key.ReadObj().GetName() for sysname in mainSys)  
        if mode == "mainSys" and "__" in key.ReadObj().GetName() and not mainSys_true: 
            continue
        #histlist.append( key.ReadObj().GetName() )
        histlist.append( key.ReadObj())
        ihist += 1

    if not isinstance(hist_totexpect, type(None)):
        hist_totexpect.SetName("total_expect")
        hist_totexpect.SetTitle("Total background MC")
        histlist.append(hist_totexpect)

    rfile.Close()
    rfile2 = ROOT.TFile(newfile,"RECREATE") 
    rfile2.cd()
    #directory = rfile2.GetDirectory("")
    for hist in histlist:
        #hist = rfile.Get(histname)
        if type(hist)  == ROOT.TH2D:
            hist = rebin_hist2d(hist, xbins, ybins)
        elif type(hist)  == ROOT.TH1D:
            hist = rebin_hist1d(hist, xbins)
        else:
            print("hist is not 2D nor 1D, rebin is ignore")
        hist.SetDirectory(rfile2) 
        #c1 = ROOT.TCanvas("c", "canvas", 800, 800)
        #c1.Clear()
        #hist.Draw("hist")
        #c1.SaveAs(newfile[:-5]+hist.GetName()+"_1d.png")
        hist.Write()
    rfile2.Close()

def rebin_onefile_bkgnode(rootfile, newfile, xbins, ybins, mode="allsys", quantile_binning=True):
    if not (type(xbins).__module__ == np.__name__) and not quantile_binning:
        print("input new bins are not array type")
        return 0
    rfile = ROOT.TFile(rootfile,"READ") 
    keylist = rfile.GetListOfKeys()
    if len(keylist) == 0:
        print("no histogram in file ", rootfile)
        return 0

    signal = []
    ihist = 0; maxhist = -1
    histlist = []
    hist_totexpect = None
    for key in keylist:
        if ihist >= maxhist and maxhist > 0:
            print("only rebins %d histograms now"%maxhist)
            break
        #isbkg_nominal = "__" not in key.ReadObj().GetName() and "data_obs" not in key.ReadObj().GetName() and "ggHH_" not in key.ReadObj().GetName()  
        isbkg_nominal = "__" not in key.ReadObj().GetName() and "data_obs" not in key.ReadObj().GetName() and "signal_ggf" not in key.ReadObj().GetName()  
        if not hist_totexpect and isbkg_nominal:
            hist_totexpect = key.ReadObj().Clone()
        elif isbkg_nominal:
            hist_totexpect.Add(key.ReadObj())

        if mode == "lnN" and "__" in key.ReadObj().GetName():
            continue
        mainSys_true = any(sysname in key.ReadObj().GetName() for sysname in mainSys)  
        if mode == "mainSys" and "__" in key.ReadObj().GetName() and not mainSys_true: 
            continue
        histlist.append( key.ReadObj())
        ihist += 1

    if not isinstance(hist_totexpect, type(None)):
        hist_totexpect.SetName("total_expect")
        hist_totexpect.SetTitle("Total background MC")
        histlist.append(hist_totexpect)
    if quantile_binning:
        xbins = find_quantile_binning(hist_totexpect, 8) ##total 8 bins
        
    rfile.Close()
    rfile2 = ROOT.TFile(newfile,"RECREATE") 
    rfile2.cd()
    #directory = rfile2.GetDirectory("")
    for hist in histlist:
        #hist = rfile.Get(histname)
        #print("hist type ", type(hist))
        if type(hist)  == ROOT.TH2D:## signal node
            hist = rebin_hist2d(hist, xbins, ybins)
        elif  type(hist)  == ROOT.TH1D: ## background node
            hist = rebin_hist1d(hist, xbins)
        else:
            print("hist type is not 1D or 2D, no rebin")
        hist.SetDirectory(rfile2) 
        #c1 = ROOT.TCanvas("c", "canvas", 800, 800)
        #c1.Clear()
        #hist.SetMinimum(0.0)
        #hist.Draw("hist")
        #c1.SaveAs(newfile[:-5]+hist.GetName()+"_1d.png")
        hist.Write()
    rfile2.Close()

def rebin_onefile_signode(rootfile, newfile, xbins, ybins, mode,  thresh, quadraticthreshalgo=False,
HMEqtlbinning=False):
    if not (type(xbins).__module__ == np.__name__ ) and not (type(ybins).__module__ == np.__name__ or HMEqtlbinning):
        print("input new bins are not array type")
        return 0
    rfile = ROOT.TFile(rootfile,"READ") 
    keylist = rfile.GetListOfKeys()
    #print("keylist ",keylist)
    if len(keylist) == 0:
        print("no histogram in file ", rootfile)
        return 0

    signal = []
    ihist = 0; maxhist = -1
    histlist = []
    hist_totexpect = None
    hist_signal = None
    for key in keylist:
        #print("processing hist ",key.ReadObj().GetName())
        if ihist >= maxhist and maxhist > 0:
            print("only rebins %d histograms now"%maxhist)
            break
        if "signal_ggf_" in key.ReadObj().GetName() and "hbbhww" in key.ReadObj().GetName() and "__" not in key.ReadObj().GetName():
            hist_signal = key.ReadObj()
        #isbkg_nominal = "__" not in key.ReadObj().GetName() and "data_obs" not in key.ReadObj().GetName() and "ggHH_" not in key.ReadObj().GetName()  
        isbkg_nominal = "__" not in key.ReadObj().GetName() and "data_obs" not in key.ReadObj().GetName() and "signal_ggf" not in key.ReadObj().GetName()  
        if isinstance(hist_totexpect, type(None)) and isbkg_nominal:
            hist_totexpect = key.ReadObj().Clone()
        elif isbkg_nominal:
            hist_totexpect.Add(key.ReadObj())

        if mode == "lnN" and "__" in key.ReadObj().GetName():
            continue
        #mainSys_true = any(sysname in key.ReadObj().GetName() for sysname in mainSys)  
        #if mode == "mainSys" and "__" in key.ReadObj().GetName() and not mainSys_true: 
        #    continue
        histlist.append( key.ReadObj())
        ihist += 1

    if not isinstance(hist_totexpect, type(None)):
        hist_totexpect.SetName("total_expect")
        hist_totexpect.SetTitle("Total background MC")
        histlist.append(hist_totexpect)
        print("adding total MC")
    ### use quantile HME binning
    if HMEqtlbinning and not isinstance(hist_signal, type(None)):
        nybin = len(ybins)-1 if len(ybins)>3 else ybins[0]
        ybins = find_quantile_HMEbinning(hist_signal, nybin)
        print("use quantile HME binning for ",newfile, ybins)
    elif HMEqtlbinning: 
        print("ERROR!!! signal histogram is not found !!, HMEqtlbinning is not working")


    histlist_rebin = []
    mainbkgs = ["TT","ST","DY"]
    mainbkg_list = []
    hist_totexpect_rebin = None
    for hist in histlist:
        if type(hist)  == ROOT.TH2D:## signal node
            hist = rebin_hist2d(hist, xbins, ybins)
            if hist.GetName() in mainbkgs:
                mainbkg_list.append(hist)
            if hist.GetName() == "total_expect":
                hist_totexpect_rebin = hist
            histlist_rebin.append(hist)
        elif  type(hist)  == ROOT.TH1D: ## background node
            histlist_rebin.append(hist)
        else:
            print("hist type is not 1D or 2D, no rebin")
            histlist_rebin.append(hist)

    if thresh > 0:
        xbins = find_threshold_binning_linearized2d(hist_totexpect_rebin, mainbkg_list, len(ybins)-1, thresh,
        quadraticthreshalgo) ##total 8 bins
        print("use threshold algorithm binning for ",newfile," xbins ",xbins)
        xbins = find_threshold_binning_linearized2d_new(hist_totexpect_rebin, mainbkg_list, len(ybins)-1, thresh,
        quadraticthreshalgo) ##total 8 bins
        print("use new threshold algorithm binning for ",newfile," xbins ",xbins)

    rfile.Close()
    rfile2 = ROOT.TFile(newfile,"RECREATE") 
    rfile2.cd()
    for hist in histlist_rebin:
        if thresh>0:
            hist = rebin_hist1d(hist, xbins)  
        hist.SetDirectory(rfile2) 
        #c1 = ROOT.TCanvas("c", "canvas", 800, 800)
        #c1.Clear()
        #hist.SetMinimum(0.0)
        #hist.Draw("hist")
        #c1.SaveAs(newfile[:-5]+hist.GetName()+"_1d.png")
        hist.Write()
    rfile2.Close()

def rebin_onefolder(oldfolder, newfolder, xbins, ybins, runyears, mode,  thresh, quadraticthreshalgo, HMEqtlbinning):
    if len(xbins) == 3:
        nbins = xbins[0]; xmin = xbins[1]; xmax =  xbins[2]
        xbins = []
        binwidth = (xmax-xmin)*1.0/nbins
        for i in range(0, nbins+1):
            xbins.append(round(xmin + i*binwidth, 4))
        xbins = np.asarray(xbins)

    if len(ybins) == 3:
        nbins = ybins[0]; ymin = ybins[1]; ymax =  ybins[2]
        ybins = []
        binwidth = (ymax-ymin)*1.0/nbins
        for i in range(0, nbins+1):
            ybins.append(round(ymin + i*binwidth, 4))
        ybins = np.asarray(ybins)

    print("old folder ", oldfolder, " new folder ", newfolder)
    print("new binning before quantile , x ",xbins, " y ",ybins)
    if not os.path.exists(newfolder):
        os.makedirs(newfolder)
    allfiles = os.listdir(oldfolder)
    ##datacard_fit_Resonant_HighMass_Graviton_2D_syst_M_550_raw_FR2
    node = "GGF"
    #node = "boosted_GGF"
    #node = "resolved1b_GGF"
    print("process file in ",node)
    for filename in allfiles:
        for year in runyears:
            if filename.endswith(".root") and year in filename and node in filename:
            #if filename.endswith(".root") and year in filename:
                signalnode = ("GGF" in filename)
                rootfile = os.path.join(oldfolder, filename)
                newfilename = None
                if mode == "lnN":
                    newfilename = filename[:-5]+"_rebin1d_%s.root"%mode
                else:
                    newfilename = filename[:-5]+"_rebin1d.root"
                newfile  = os.path.join(newfolder, newfilename)
                print("oldfile ", rootfile, " new file ", newfile)
                #if os.path.exsits(newfile):
                #    rinfile = ROOT.TFile(rootfile, "READ")
                #    keylist_in = rinfile.GetListOfKeys()
                #    routfile = ROOT.TFile(newfile, "READ")
                #    keylist_out = routfile.GetListOfKeys()
                #    rinfile.Close()
                #    routfile.Close()
                #    if len(keylist_in)  == len(keylist_out)-1:
                #        print("Rebinning is already processed, remove this condition if rebinning is forced again")
                #        continue
            
                if signalnode:
                    rebin_onefile_signode(rootfile, newfile, xbins, ybins, mode, thresh, quadraticthreshalgo,
                    HMEqtlbinning)
                else:
                    rebin_onefile_bkgnode(rootfile, newfile, xbins, ybins, mode)## use quantile binning for bkg node
        
    
def checkbrokenfile(oldfolder):
    ### note: found HH_800_boosted_GGF_2016.root is broken for v1 
    print("checking folder ", oldfolder)
    allfiles = os.listdir(oldfolder)
    for filename in allfiles:
        if filename.endswith(".root") and "rebin1d" not in filename:
            rootfile = os.path.join(oldfolder, filename)
            rfile = ROOT.TFile(rootfile, "READ")
            rfile.Close()
    
def main():
    xbins= [10,  0.0, 1.0]
    masspoints = [260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 900]
    ybins = [250.0, 287, 330, 379, 435, 500, 575, 661, 760, 873, 1003, 1153, 1200.0, 2000.0]
    #HME_bins = [0, 1, 2, 3]
    HME_bins= [1.1, 1.15, 1.2]
    nn_bins = [10, 20, 40]
    threshlist = [2, 5, 10] 
    mode = "allsys"
    quadraticthreshalgo = False
    HMEquantilebinning = False

    if HMEquantilebinning:
        HME_bins = [5, 10, 15]
    if quadraticthreshalgo:
        threshlist = [10]


    if HMEquantilebinning and quadraticthreshalgo:
        print("Running rebin mode with HME quantile binning and quadratic threshold algo, and systematics: ", mode)
    elif HMEquantilebinning and not quadraticthreshalgo:
        print("Running rebin mode with HME quantile binning and fixed threshold algo, and systematics: ", mode)
    elif not HMEquantilebinning and quadraticthreshalgo:
        print("Running rebin mode with customized HME binning and quadratic threshold algo, and systematics: ", mode)
    else:
        print("Running rebin mode with customized HME binning and fixed threshold algo, and systematics: ", mode)

    #masspoints = [600]
    ##masspoints = [800]
    #HME_bins=[3]
    #nn_bins=[10]
    print("arguments number ",len(sys.argv), sys.argv)
    ## arguments: 3 + '-b' +'-q'
    condormode = len(sys.argv) >= 3+2
    if condormode:
        print("runnning rebin with condor mode", sys.argv)
        if len(sys.argv) == 7:
            nn_bins = [int(sys.argv[3])]
            threshlist = [int(sys.argv[4])]
        elif len(sys.argv) == 6:
            nn_bins = [int(sys.argv[3])]
            threshlist = [2, 5, 10] 
        elif len(sys.argv) == 5:
            nn_bins =  [10, 20, 40, 400]
            threshlist = [2, 5, 10] 

        masspoints = [int(sys.argv[1])]
        HME_bins = [float(sys.argv[2])]

    folder = None
    #if mass <= 500:
    #    folder = "/eos/user/t/tahuang/FR2_Florian/datacard_fit_Resonant_LowMass_Graviton_2D_syst_M_%d_raw_FR2"%mass
    #else:
    #    folder = "/eos/user/t/tahuang/FR2_Florian/datacard_fit_Resonant_HighMass_Graviton_2D_syst_M_%d_raw_FR2"%mass
    #newfolder = folder+"_rebin%HMEv%d"%hmebin_version
    #rebin_onefolder(mass, hmebin_version, folder, newfolder, xbins, ["2016", "2017","2018"])
        
    #outputfolder = "/afs/cern.ch/work/t/tahuang/CombinedLimit/ForGithub/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/Rebin_Florian"
    ## put rebinned root file to eos
    eospath = "/eos/user/t/tahuang/"
    outputfolder = "HHbbww_resonant_DL_Rebin_v2_202111"
    if condormode:
        outputfolder = os.path.join(eospath, outputfolder)

    for mass in masspoints:
      if mass <= 500:
          folder = "datacard_fit_Resonant_LowMass_Graviton_2D_syst_M_%d_raw_FR2"%mass
      else:
          folder = "datacard_fit_Resonant_HighMass_Graviton_2D_syst_M_%d_raw_FR2"%mass
      if condormode:
          folder = os.path.join(eospath+"FR2_Florian_v2", folder)
      #checkbrokenfile(folder)
      #for hmebin_version in HME_bins:
      for hmebin_version in HME_bins:
          for nnbin_version in nn_bins:
              for thresh in threshlist:
                  #newname = folder.split('/')[-1] + "_rebin_nnbin%dHMEv%d"%(nnbin_version, hmebin_version)
                  binsuffix = "nnbin%dHME%.2fthresh%d"%(nnbin_version, hmebin_version, thresh)
                  if HMEquantilebinning and not quadraticthreshalgo:
                      binsuffix = "nnbin%dHME%dqtlthresh%d"%(nnbin_version, hmebin_version, thresh)
                  elif HMEquantilebinning and quadraticthreshalgo:
                      binsuffix = "nnbin%dHME%dqtlthreshquadratic"%(nnbin_version, hmebin_version)
                  elif not HMEquantilebinning and quadraticthreshalgo:
                      binsuffix = "nnbin%dHME%.2fthreshquadratic"%(nnbin_version, hmebin_version)
                  binsuffix = binsuffix.replace('.','p')
                      
                  newname = "Graviton_syst_allbenchmarks_FR2_rebin_%s"%binsuffix
                  #newfolder = folder+"_rebin_nnbin%dHMEv%d"%(nnbin_version, hmebin_version)
                  newfolder = os.path.join(outputfolder, newname)
                  #newfolder = folder
                  xbins = [nnbin_version, 0.0, 1.0]
                  #ybins = generateHMEbins(int(mass), hmebin_version) for v1 only
                  if HMEquantilebinning:
                      ybins = [int(hmebin_version), 0.0, 2000.0]
                  else:
                      ybins = generateHMEbinsv4(hmebin_version)## mass-independent  for v2
                  rebin_onefolder(folder, newfolder, xbins, ybins, ["2016", "2017","2018"], mode, thresh, quadraticthreshalgo,
                  HMEquantilebinning)
                  #os.system("cp "+folder+"/*txt "+newfolder+"/")
                  #rebin_onefolder(folder, newfolder, xbins, ybins, ["2016"], mode, thresh, quadraticthreshalgo,
                  #HMEquantilebinning)
                  #continue
if __name__ == '__main__':
    print("starting to main ")
    parser = argparse.ArgumentParser(description='Process rebin_florian arguments')
    #parser.add_argument("outfile", help="output file to rebin ")
    parser.add_argument("--infile", help="input file to rebin ")
    parser.add_argument("--nnbin_version", help="number of nn output bins", type=int)
    parser.add_argument("--hme_version", help="number of HME binning or HME binning version", type=float)
    parser.add_argument("--thresh", help="threshold for threshold algorithm", type=int)
    parser.add_argument("--quadraticthreshalgo", help="use quadratic threshold or not for threshold algo", type=int)
    parser.add_argument("--HMEquantilebinning", help="use quantile binning for HME or not", type=int)
    args = parser.parse_args()
    print("infile ", args.infile)
    print("nnbin_version", args.nnbin_version)
    print("hme_version", args.hme_version)
    print("thresh", args.thresh)
    print("quadraticthreshalgo ", args.quadraticthreshalgo)
    print("HMEquantilebinning ", args.HMEquantilebinning)

    mode = "allsys"
    infile = args.infile.split('/')[-1]
    #outfile = args.outfile
    outfile = infile[:-5]+"_rebin1d.root"
    nnbin_version = args.nnbin_version
    hmebin_version = args.hme_version
    thresh = args.thresh
    quadraticthreshalgo = args.quadraticthreshalgo != 0
    HMEquantilebinning = args.HMEquantilebinning != 0

    if HMEquantilebinning and quadraticthreshalgo:
        print("Running rebin mode with HME quantile binning and quadratic threshold algo, and systematics: ", mode)
    elif HMEquantilebinning and not quadraticthreshalgo:
        print("Running rebin mode with HME quantile binning and fixed threshold algo, and systematics: ", mode)
    elif not HMEquantilebinning and quadraticthreshalgo:
        print("Running rebin mode with customized HME binning and quadratic threshold algo, and systematics: ", mode)
    elif thresh > 0:
        print("Running rebin mode with customized HME binning and fixed threshold algo, and systematics: ", mode)
    else:
        print("Running rebin mode with customized HME binning and Equal NN binning, and systematics: ", mode)


    #newname = folder.split('/')[-1] + "_rebin_nnbin%dHMEv%d"%(nnbin_version, hmebin_version)
    binsuffix = "nnbin%dHME%.2fthresh%d"%(nnbin_version, hmebin_version, thresh)
    if HMEquantilebinning and not quadraticthreshalgo:
      binsuffix = "nnbin%dHME%dqtlthresh%d"%(nnbin_version, hmebin_version, thresh)
    elif HMEquantilebinning and quadraticthreshalgo:
      #binsuffix = "nnbin%dHME%dqtlthreshquadratic"%(nnbin_version, hmebin_version)
      binsuffix = "nnbin%dHME%dqtlthreshquadraticpoisson"%(nnbin_version, hmebin_version)
    elif not HMEquantilebinning and quadraticthreshalgo:
      #binsuffix = "nnbin%dHME%.2fthreshquadratic"%(nnbin_version, hmebin_version)
      binsuffix = "nnbin%dHME%.2fthreshquadraticpoisson"%(nnbin_version, hmebin_version)
    if thresh == 0:##equalNN binnings
        binsuffix = "nnbin%dHME%.2fNothresh"%(nnbin_version, hmebin_version) 
    binsuffix = binsuffix.replace('.','p')
                      
    xbins = [nnbin_version, 0.0, 1.0]
    #ybins = generateHMEbins(int(mass), hmebin_version) for v1 only
    ybins = None
    if HMEquantilebinning:
      ybins = [int(hmebin_version), 0.0, 2000.0]
    else:
      ybins = generateHMEbinsv4(hmebin_version)## mass-independent  for v2

    if len(xbins) == 3:
        nbins = xbins[0]; xmin = xbins[1]; xmax =  xbins[2]
        xbins = []
        binwidth = (xmax-xmin)*1.0/nbins
        for i in range(0, nbins+1):
            xbins.append(round(xmin + i*binwidth, 4))
        xbins = np.asarray(xbins)

    if len(ybins) == 3:
        nbins = ybins[0]; ymin = ybins[1]; ymax =  ybins[2]
        ybins = []
        binwidth = (ymax-ymin)*1.0/nbins
        for i in range(0, nbins+1):
            ybins.append(round(ymin + i*binwidth, 4))
        ybins = np.asarray(ybins)
    print("Before quantile binning, xbins ", xbins, " ybins ",ybins," thresh ", thresh, " binsuffix ", binsuffix)
    rebin_onefile_signode(infile, outfile, xbins, ybins, mode,  thresh, quadraticthreshalgo, HMEquantilebinning)

    dirname = "Graviton_syst_allbenchmarks_FR2_rebin_%s"%binsuffix
    eospath = os.path.join("/eos/user/t/tahuang/HHbbww_resonant_DL_Rebin_v2_202111", dirname)
    if not os.path.exists(eospath):
        os.mkdir(eospath)
    print("infile ", infile, " now copying the output ", outfile, " to ", eospath)
    os.system("xrdcp "+outfile +" "+eospath)
