#! /bin/env python

import argparse
import os


import ROOT as R

#R.PyConfig.IgnoreCommandLineOptions = True
#R.gROOT.Reset()
#R.gROOT.SetBatch(True)

R.gStyle.SetOptStat(0)
#R.gStyle.SetOptTitle(0)
R.gStyle.SetPaintTextFormat(".2f")
R.gStyle.SetTextFont(42)

title = " #scale[1.4]{#font[61]{CMS}} Simulation Preliminary"+"  "*6+"35.87 fb^{-1} (13 TeV),2016"
def drawEfficiency(obj, plotname):
    # Ignore systematics for the moment
    # if "__" in obj.GetName():
        # return

    c = R.TCanvas("c", "c", 600, 800)

    h = obj.GetCopyPassedHisto()
    h.Reset()
    h.SetDirectory(0)

    for x in range(1, h.GetNbinsX() + 1):
        #for y in range(1, h.GetNbinsY() + 1):
            bin = h.GetBin(x)
            eff = obj.GetEfficiency(bin)
            err = max(obj.GetEfficiencyErrorLow(bin), obj.GetEfficiencyErrorUp(bin))

            h.SetBinContent(x, eff)
            h.SetBinError(x, err)

    h.GetXaxis().SetTitle("BDT output")
    h.GetYaxis().SetTitle("Fraction")
    h.GetYaxis().SetTitleOffset(1.2)
    h.SetTitle(obj.GetName())

    #c.SetLogx()
    #h.SetTitle(title)
    #h.SetContour(100)
    #h.Draw("colz text err")
    h.Draw("e")
    c.SaveAs(plotname + ".pdf")



parser = argparse.ArgumentParser(description='Draw b-tagging efficiency')
parser.add_argument('-i','--input', dest='input', type=str, metavar='FILE', help='Input file containing efficiencies')
parser.add_argument('-o','--output', dest='outputdir', type=str,default='./', help='out put dir')
parser.add_argument('-n','--name', dest='name', type=str, default='', help='plotname for this file')


options = parser.parse_args()

flavours = ['b','c','l']
f = R.TFile(options.input,"READ")
print "f ",f, " inputs ",options.input
def drawEach():
    for key in f.GetListOfKeys():
        clazz = R.gROOT.GetClass(key.GetClassName())
        if not clazz.InheritsFrom("TEfficiency"):
            continue
        obj = key.ReadObj()
        if "vs_pt" in obj.GetName() or "vs_eta" in obj.GetName():
            continue
        plotname  = os.path.join(options.outputdir, obj.GetName()+"_"+options.name)
        drawEfficiency(obj, plotname)

c = R.TCanvas("c", "c", 900, 900)
c.Clear()
c.Divide(3,3)
ipad = 1
for flav1 in flavours:
    for flav2 in flavours:
        c.cd(ipad)
        name = "%s%s_frac" % (flav1, flav2)
	eff = f.Get(name)
	#b1 = eff.GetCopyPassedHisto()
	#b1.Reset()
	#b1.SetDirectory(0)
        #b1.SetTitle(name)
    	#b1.Set
        #b1.Draw()
        eff.Draw()
        c.Update()
        #print eff.GetPaintedGraph().GetTitle()
        title = R.gPad.GetPrimitive("title")
        title.SetTextFont(62)
        title.SetTextSize(.07)
	ipad += 1
c.cd()
c.Update()
plotname  = os.path.join(options.outputdir, "DY_fraction_all_"+options.name)
c.SaveAs(plotname + ".pdf")

