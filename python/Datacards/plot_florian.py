import os
from array import array
from datetime import datetime
#import keras
from math import sqrt
import numpy.lib.recfunctions as recfunctions
import numpy as np
import ROOT

import sys
sys.argv.append( '-b' )
sys.argv.append( '-q' )

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
## ggHH_M_550_hbbhwwdl,ggHH_M_550_hbbhtt
bkg_dict = {
        "TT": ["TT"],
        "ST": ["ST"],
        "DY": ["DY"],
        "Fakes": ["Fakes"],
        "SM Higgs":
        ["ggH_hbb","ggH_hgg","ggH_hmm","ggH_htt","ggH_hww","ggH_hzz","qqH_hbb","qqH_hgg","qqH_hmm","qqH_htt","qqH_hww","WH_hbb","ZH_hbb","ZH_htt","ZH_hww","tHq_hww","tHW_hww","VH_hww"],
        "VV(V)" : ["VV","VVV"],
        "W+jets":["WJets"],
        "Rares" : ["ttZ", "ttW","Other_bbWW"]
        }

eventcat = ["boosted_DY_VVV", "boosted_GGF", "boosted_other", "resolved1b_GGF", "resolved2b_GGF","resolved_DY_VVV","resolved_other"]
eventcat = ["boosted_GGF", "resolved1b_GGF", "resolved2b_GGF"]
runyears = ["2016", "2017","2018"]
intLumi_years = {
 "2016" : "35.9 fb^{-1}",
 "2017" : "41.53 fb^{-1}",
 "2018" : "59.74 fb^{-1}",
 "FR2" : "147.17 fb^{-1}"
}

def plot1d(mass, year, infiles, category, bkgnames,  plotdir, plotsuffix, adduncertainty, adddata):

    colors = [900+2, 800-4, 860-9, 632-9, 616-10, 432-6, 600-6, 400-6]
    bkg_hists = {}
    totalbkg = None
    hist_bbww = None
    hist_bbtt = None
    hist_bbww_name = "signal_ggf_spin2_%d_hbbhww"%mass
    hist_bbtt_name = "signal_ggf_spin2_%d_hbbhtt"%mass
    print("histname ", hist_bbww_name, hist_bbtt_name)
    for f in infiles:
        rf = ROOT.TFile(f, "READ")

        if isinstance(hist_bbww, type(None)):
            hist_bbww = rf.Get(hist_bbww_name).Clone()
            hist_bbww.SetDirectory(0)
        else:
            hist_bbww.Add(rf.Get(hist_bbww_name).Clone())
        if isinstance(hist_bbtt, type(None)):
            hist_bbtt = rf.Get(hist_bbtt_name).Clone()
            hist_bbtt.SetDirectory(0)
        else:
            hist_bbtt.Add(rf.Get(hist_bbtt_name).Clone())

        for bkg in bkgnames:
            histname_list = bkg_dict[bkg]
            for name in histname_list:
                hist = rf.Get(name)
                if not(type(hist) == ROOT.TH1F or type(hist) == ROOT.TH1D):
                    print("Failed to find histogram ", hist ," in file ", f)
                    continue
                if bkg not in bkg_hists.keys():
                    bkg_hists[bkg] = hist.Clone()
                    bkg_hists[bkg].SetDirectory(0)
                else:
                    bkg_hists[bkg].Add(hist.Clone())

    hs = ROOT.THStack("bbWW_"+category, "  ")
    legend = ROOT.TLegend(0.6,0.65,0.88,0.68+len(bkgnames)*.024);
    legend.SetNColumns(2);
    legend.SetTextSize(0.03); legend.SetTextFont(42)
    legend.SetBorderSize(0)
    num_bkg = len(bkgnames)
    for i, bkg in enumerate(bkgnames):
        bkg_hists[bkg].SetFillColor(colors[i])

        hs.Add(bkg_hists[bkg])
        if isinstance(totalbkg, type(None)):
            totalbkg = bkg_hists[bkg].Clone()
        else:
            totalbkg.Add(bkg_hists[bkg].Clone())

        legend.AddEntry(bkg_hists[bkgnames[num_bkg-1-i]], bkgnames[num_bkg-1-i], "f")

    maxbgbin = totalbkg.GetBinContent(totalbkg.GetMaximumBin())
    max_bbww = hist_bbww.GetBinContent(hist_bbww.GetMaximumBin())
    max_bbtt = hist_bbtt.GetBinContent(hist_bbtt.GetMaximumBin())
    maxbin = max(max_bbww, max_bbtt, maxbgbin)

    hs.SetMaximum(maxbin*2.0)

    hist_data = ROOT.TH1D()
    hist_data.SetMarkerStyle(20)
    legend2 = ROOT.TLegend(0.15,0.73,0.5,0.73+len(bkgnames)*.012);
    legend2.AddEntry(hist_data, "Data", "p")
    legend2.AddEntry(hist_bbww, "Graviton M=%d GeV, HH->bbWW [1pb]"%mass,"l") 
    legend2.AddEntry(hist_bbtt, "Graviton M=%d GeV, HH->bb#tau#tau [1pb]"%mass,"l") 
    legend2.SetTextSize(0.025)
    legend2.SetBorderSize(0)
    hist_bbww.SetLineColor(416+1)
    hist_bbww.SetLineWidth(3)
    hist_bbww.SetFillStyle(3315)
    hist_bbtt.SetLineColor(632+2)
    hist_bbtt.SetLineWidth(3)
    hist_bbtt.SetFillStyle(3351)

    c1 = ROOT.TCanvas("c", "canvas", 800, 800)
    c1.SetLogy()
    #c1.Clear()
    #pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    #pad1.SetBottomMargin(.0)
    #pad1.SetLogy()
    #pad1.SetLogx()
    #pad1.Draw()
    #pad1.cd()
    hs.Draw("hist")
    hist_bbww.Draw("samehist")
    hist_bbtt.Draw("samehist")

    legend.Draw("same")
    legend2.Draw("same")
    
    tex0 = ROOT.TLatex(0.1,0.92, " #scale[1.4]{#font[61]{CMS}} Preliminary"+"  "*10 + intLumi_years[year]+"(13TeV)")
    #tex0.SetNDC(); tex0.SetTextSize(.05); tex0.SetTextFont(42)
    tex0.SetNDC(); tex0.SetTextSize(.04); tex0.SetTextFont(42)
    tex0.Draw("same")

    hs.GetHistogram().GetYaxis().SetTitle("Events")
    #hs.GetHistogram().GetYaxis().SetTitleSize(.05)
    hs.GetHistogram().GetYaxis().SetTitleSize(.04)
    hs.GetHistogram().GetYaxis().SetLabelFont(42)
    hs.GetHistogram().GetYaxis().SetLabelSize(.045)
    hs.GetHistogram().GetXaxis().SetTitle("Linearized 1D bin")
    hs.GetHistogram().GetYaxis().SetTitleOffset(1.1)
    #hs.GetHistogram().GetXaxis().SetTitleSize(.06)
    hs.GetHistogram().GetXaxis().SetTitleSize(.04)
    hs.GetHistogram().GetXaxis().SetLabelFont(42)
    #c1.cd()
    #c1.Update()

    #pad2 = ROOT.TPad("pad2", "pad2", 0, 0.0, 1, .29)
    #pad2.SetTopMargin(0.)
    #pad2.SetBottomMargin(.35)
    #pad2.SetTicks(1,1)#draw x, y axis on both side (left right for y, and bottom up for x)
    #pad2.Draw()
    #pad2.cd()

    #tex_pad2 = ROOT.TLatex(0.2,0.35, "Maximum #frac{S}{#sqrt{B}} = %.1f @ %.3f"%(bestS, bestWP))
    #tex_pad2.SetNDC()
    #tex_pad2.SetTextSize(.035)
    #tex_pad2.Draw("same")
 
    c1.SaveAs(plotdir+"Graviton_"+category+"_"+plotsuffix+".C")
    c1.SaveAs(plotdir+"Graviton_"+category+"_"+plotsuffix+".png")
    c1.SaveAs(plotdir+"Graviton_"+category+"_"+plotsuffix+".pdf")

def plot_all(masspoints, folder, binsuffix, mode, bkgnames, plotdir, adduncertainty, adddata):
     # 
    for mass in masspoints:
        for cat in eventcat:
            FR2_file = []
            for year in runyears:
                category = "HH_DL_{mass}_{cat}_{year}".format(mass = mass, cat = cat, year = year)
                #filename = category+"_rebin1d_%s.root"%mode
                filename = category+"_rebin1d.root"
                infiles = []
                infiles.append(os.path.join(folder, filename))
                FR2_file.append(os.path.join(folder, filename))
                #print("infiles ", infiles, " cat ", cat, " year ", year)
                plot1d(mass, year, infiles, category, bkgnames,  plotdir, binsuffix, adduncertainty, adddata)
            category = "HH_{mass}_{cat}_{year}".format(mass = mass, cat = cat, year = "FR2")
            print("FR2 infiles ", FR2_file, " cat ", cat)
            #plot1d(mass, "FR2", FR2_file, category, bkgnames,  plotdir, binsuffix, adduncertainty, adddata)



runyears_all = ["2016", "2017","2018","FR2"]
bkgnames = ["Fakes","Rares","W+jets","VV(V)","SM Higgs","DY","ST","TT"]

adduncertainty=False
adddata = False
mode = "lnN"

dateTimeObj = datetime.now()
eospath = "/eos/user/t/tahuang/"
rebinfolder="HHbbww_resonant_DL_Rebin_v2_202111/"
rebinfolder = os.path.join(eospath, rebinfolder)
plotdir = "HHbbWW_final_shapes_rebin_%s_"%mode+dateTimeObj.strftime("%d-%b-%Y")

#masspoints = [260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 900]
masspoints = [260, 270, 300, 350, 400, 450, 500, 600, 650, 700, 800, 900]
#masspoints = [600]

EqualNNbin = False
irange = 4
if EqualNNbin:
    irange =1
    plotdir += "_equalNNbins"
if not os.path.exists(plotdir):
    os.makedirs(plotdir)
fullpath = "/afs/cern.ch/work/t/tahuang/CombinedLimit/ForGithub/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/Rebin_Florian/"
plotdir = os.path.join(fullpath, plotdir)

for i in range(irange):
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
    for hmebin_version in HME_bins:
        for nnbin_version in nn_bins:
            for thresh in thresholds:
                binsuffix = "nnbin%dHME%.2fthresh%d"%(nnbin_version, hmebin_version, thresh)
                if HMEquantile_binning and not quadraticthreshalgo:
                    binsuffix = "nnbin%dHME%dqtlthresh%d"%(nnbin_version, hmebin_version, thresh)
                elif HMEquantile_binning and quadraticthreshalgo:
                    binsuffix = "nnbin%dHME%dqtlthreshquadratic"%(nnbin_version, hmebin_version)
                elif not HMEquantile_binning and quadraticthreshalgo:
                    binsuffix = "nnbin%dHME%.2fthreshquadratic"%(nnbin_version, hmebin_version)
                if thresh == 0:##equalNN binnings
                    binsuffix = "nnbin%dHME%.2fNothresh"%(nnbin_version, hmebin_version) 
                binsuffix = binsuffix.replace('.','p')
                    
                newname = "Graviton_syst_allbenchmarks_FR2_rebin_%s"%binsuffix
                thisfolder = os.path.join(rebinfolder, newname)
                if plotdir[-1] != "/":
                    plotdir += "/"
                plot_all(masspoints, thisfolder, binsuffix, mode, bkgnames, plotdir, adduncertainty, adddata)





