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

    c = R.TCanvas("c", "c", 1600, 800)

    h = obj.GetCopyPassedHisto()
    h.Reset()
    h.SetDirectory(0)

    for x in range(1, h.GetNbinsX() + 1):
        for y in range(1, h.GetNbinsY() + 1):
            bin = h.GetBin(x, y)
            eff = obj.GetEfficiency(bin)
            err = max(obj.GetEfficiencyErrorLow(bin), obj.GetEfficiencyErrorUp(bin))

            h.SetBinContent(x, y, eff)
            h.SetBinError(x, y, err)

    h.GetXaxis().SetTitle("p_{T} of jet")
    h.GetYaxis().SetTitle("|#eta| of jet")
    h.GetYaxis().SetTitleOffset(1.2)

    c.SetLogx()
    if "eff_on_q" in obj.GetName() or "eff_on_l" in obj.GetName() or "eff_on_g" in obj.GetName() or "eff_on_n" in obj.GetName():
        c.SetLogz()
    h.SetTitle(title)
    h.SetContour(100)
    h.Draw("colz text err")
    c.SaveAs(plotname + ".pdf")



def plotTeff1D(Teffs,  xtitle, title, legs, text, setlogx,  picname):
    color = [R.kRed, R.kBlue, R.kGreen+2, R.kBlack]
    marker = [25, 21, 26, 20, 24]

    #b1 = ROOT.TH1F("b1","b1",len(xbinning)-1, xbinning)
    b1 = Teffs[0].GetCopyPassedHisto()
    b1.Reset()
    b1.SetDirectory(0)
    b1.GetYaxis().SetRangeUser(0.001,1.2)
    b1.GetYaxis().SetTitleOffset(1.0)
    b1.GetYaxis().SetNdivisions(520)
    b1.GetYaxis().SetTitle("Medium B-tagging efficiency")
    b1.GetXaxis().SetTitle(xtitle )
    b1.GetXaxis().SetTitleSize(0.05)
    b1.GetXaxis().SetLabelSize(0.05)
    b1.GetYaxis().SetTitleSize(0.05)
    b1.GetYaxis().SetLabelSize(0.05)
    #b1.SetTitle(" #scale[1.4]{#font[61]{CMS}} #font[52]{Simulation preliminary}"+"  "*8+" 14 TeV, 200 PU")
    #b1.SetTitle(" #scale[1.4]{#font[61]{CMS}} Phase-2 Simulation"+"  "*4+" #sqrt{s}=14 TeV, <PU>=200")
    #b1.SetTitle(title)
    b1.SetTitle("")

    c1 = R.TCanvas()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetTickx()
    c1.SetTicky()
    c1.SetLogy()
    if setlogx: 
	c1.SetLogx()

    dy_leg = len(Teffs)*0.05
    legend = R.TLegend(0.65,0.4,0.85,0.4+dy_leg)
    legend.SetFillColor(R.kWhite)
    legend.SetTextFont(42)
    #legend.SetTextSize(.05)
    #legend.SetHeader("%s"%legheader)
    b1.SetStats(0)
    b1.Draw()
    for n in range(len(Teffs)):
	Teffs[n].SetLineColor(color[n])
	Teffs[n].SetMarkerColor(color[n])
	Teffs[n].SetLineWidth(2)
	Teffs[n].SetMarkerStyle(marker[n])
        Teffs[n].SetTitle("")
	Teffs[n].Draw("samepZ")

	legend.AddEntry(Teffs[n],"%s"%legs[n],"pl")
    legend.Draw("same")

    tex = R.TLatex(0.4,0.65,"%s"%text)
    tex.SetTextSize(0.04)
    tex.SetTextFont(62)
    tex.SetNDC()
    tex.Draw("same")
    tex0 = R.TLatex(0.12,0.92, title)
    tex0.SetTextSize(0.04)
    #tex0.SetTextFont(62)
    tex0.SetNDC()
    tex0.Draw("same")
    #c1.Update()
    c1.SaveAs("%s.pdf"%(picname))

parser = argparse.ArgumentParser(description='Draw b-tagging efficiency')
parser.add_argument('-i','--input', dest='input', type=str, metavar='FILE', help='Input file containing efficiencies')
parser.add_argument('-o','--output', dest='outputdir', type=str,default='./', help='out put dir')
parser.add_argument('-n','--name', dest='name', type=str, default='', help='plotname for this file')


options = parser.parse_args()

flavours = ['b','c','light']
effvseta_list = []
effvspt_list = []
legs = []
f = R.TFile(options.input,"READ")
print "f ",f, " inputs ",options.input
for flavour in flavours:
    effvseta_list.append(f.Get('btagging_eff_on_%s_vs_eta'%flavour))
    effvspt_list.append(f.Get('btagging_eff_on_%s_vs_pt'%flavour))
    legs.append(flavour+"-jet")

effetaplotname = os.path.join(options.outputdir, "btagging_eff_vs_eta_1D_%s"%options.name)
effptplotname = os.path.join(options.outputdir, "btagging_eff_vs_pt_1D_%s"%options.name)
plotTeff1D(effvseta_list, "|#eta| of jet", title, legs, "", False,  effetaplotname)
plotTeff1D(effvspt_list, "|p_{T}| of jet", title, legs, "", True,  effetaplotname)
##draw2D
for key in f.GetListOfKeys():
    clazz = R.gROOT.GetClass(key.GetClassName())
    if not clazz.InheritsFrom("TEfficiency"):
        continue
    obj = key.ReadObj()
    if "vs_pt" in obj.GetName() or "vs_eta" in obj.GetName():
	continue
    plotname  = os.path.join(options.outputdir, obj.GetName()+"_"+options.name)
    drawEfficiency(obj, plotname)
