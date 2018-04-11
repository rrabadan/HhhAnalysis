import sys
import ROOT 
ROOT.gROOT.SetBatch(1)
#from Helpers import *
ROOT.gErrorIgnoreLevel=1001
from ROOT import * 
import os
import numpy as np
from math import *

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetStatW(0.07)
ROOT.gStyle.SetStatH(0.06)

#ROOT.gStyle.SetOptStat(111110)
ROOT.gStyle.SetOptStat(0)

#ROOT.gStyle.SetErrorX(0)
#ROOT.gStyle.SetErrorY(0)

ROOT.gStyle.SetTitleStyle(0)
ROOT.gStyle.SetTitleAlign(13) ## coord in top left
ROOT.gStyle.SetTitleX(0.)
ROOT.gStyle.SetTitleY(1.)
ROOT.gStyle.SetTitleW(1)
ROOT.gStyle.SetTitleH(0.058)
ROOT.gStyle.SetTitleBorderSize(0)

ROOT.gStyle.SetPadLeftMargin(0.126)
ROOT.gStyle.SetPadRightMargin(0.10)
ROOT.gStyle.SetPadTopMargin(0.06)
ROOT.gStyle.SetPadBottomMargin(0.13)

ROOT.gStyle.SetMarkerStyle(1)
etabin = [1.0, 1.2, 1.3,  1.4, 1.5, 1.6, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95,
	      2.0, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5]
myetabin = np.asarray(etabin)
binLow = [0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,12.0,14.0,16.0,18.0,20.0,24.0,28.0,32.0,36.0,42.0,50.0]
ptbins = np.asarray(binLow)
rateptbin = [5.0, 6.0, 7.0,8.0, 10.0, 12.0, 15.0, 20.0, 25.0, 30.0, 35.0]
myrateptbin = np.asarray(rateptbin)

ratesample = "RateTree_GE11GE21_displacedMu.root"
ratesample = "RateTree_GE11GE21_displacedmu_noveto_ME0eta21.root"
ratesample = "RateTree_GE11GE21_displacedmu_ME0eta215_20171009.root"

#signalsample = "../20170712/SingleMu_EMTF_ME0_PU200_GEMCSCANA_20170724_GE11GE21_ME0eta205.root"
#signalsample = "GEMCSCAna_PU200_192_TenMuE0To50_20171114.root"
#signalsample = "GEMCSCAna_PU200_192_DisplacedMuonGun_20171122.root"
signalsample = "GEMCSCAna_PU200_192_DaskSUSY_ctau1000_20171122.root"
#signalsampleGE11 = "../20170712/SingleMu_EMTF_ME0_PU200_GEMCSCANA_20170724_GE11only_ME0eta205.root"
#signalsampleGE11 = "../20170712/SingleMu_EMTF_ME0_PU200_GEMCSCANA_20170724_GE11GE21_ME0eta205.root"
#signalsampleGE11 = "GEMCSCAna_PU200_192_DisplacedMuonGun_20171122.root"
signalsampleGE11 = "GEMCSCAna_PU200_192_DaskSUSY_ctau1000_20171122.root"
rfile = TFile(ratesample, "READ")
#total = rfile.Get("h_eventcount").GetEntries()
total = rfile.Get("h_eventcountME0Segment192").GetEntries()
#total = 284900.000000
SF = 30.0*1000.0/total#khz
SF140 = 30.0*1000.0/total*0.7#khz
pu = 200
print "Rate sample, total events ",total," SF ",SF
#color = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen+2, ROOT.kMagenta+2,ROOT.kCyan+2]
#maker = [20,21,22,23,33]
color = [ROOT.kBlack, ROOT.kBlue, ROOT.kGreen+2, ROOT.kRed, ROOT.kOrange+8] #ROOT.kMagenta+2,ROOT.kOrange+8]
color = [ROOT.kBlue, ROOT.kGreen+2, ROOT.kRed, ROOT.kOrange+8, ROOT.kBlack]
maker =[24, 21,25,20,22]
evenodds = ["odd,even","odd,odd","even,even","even,odd","all pairs"]

def DecompresZipFiles(path):

    if os.path.isdir(path):
      ls = os.listdir(path)
      for x in ls:
	if not(x.endswith(".zip")):
	    #print "x.endswith(.root) ", x.endswith(".root")
	    continue
	os.system("unzip %s"%(path+x))
    elif os.path.isfile(path):
      os.system("unzip "+path)
    else:
      print " it is not file or dir ", filedir[0]
def addTwoHistogram(hist1, hist2):
    nbins = hist1.GetXaxis().GetNbins()
    for i in range(1, nbins+1):
	binwidth = hist2.GetXaxis().GetBinWidth(1)
	content1 = hist1.GetBinContent(i)
	upedge1 = hist1.GetXaxis().GetBinUpEdge(i)
	lowedge1 = hist1.GetXaxis().GetBinLowEdge(i)
	bincenter2 = upedge1-binwidth
	xbin2 = hist2.FindBin(upedge1-binwidth)
	content2 = hist2.GetBinContent(xbin2)
	content = content1 + content2
	while (bincenter2 > lowedge1):
	    bincenter2 = bincenter2-binwidth
	    xbin2 = hist2.FindBin(bincenter2)
	    content2 = hist2.GetBinContent(xbin2)
	    content += content2
	hist1.SetBinContent(i, content)
	print "i ",i," bincenter ",bincenter2, " xbin2 ",xbin2," content1 ",content1," content2 ",content2

def addFilesToChain(ch, filedir):
    if os.path.isdir(filedir):
    	  ls = os.listdir(filedir)
    	  for x in ls:
	      	if not(x.endswith(".root")):
			#print "x.endswith(.root) ", x.endswith(".root")
			continue
		x = filedir[:]+x
    		if os.path.isdir(x):
			continue
		ch.Add(x)
    elif os.path.isfile(filedir):
	  ch.Add(filedir)
    else:
	  print " it is not file or dir ", filedir
	  
def gethist1D(chain, den, todraw="pt"):
    h1 = ROOT.TH1F("h1","h1",len(ptbins)-1, ptbins)
    chain.Draw("fabs(%s) >> h1"%todraw,den)
    return h1

def gethist1DEta(chain, den, todraw="eta"):
    h1 = ROOT.TH1F("h1","h1",50, 0.0, 2.5)
    chain.Draw("fabs(%s) >> h1"%todraw,den)
    return h1

def gethist1DDxy(chain, den, todraw="dxy"):
    h1 = ROOT.TH1F("h1","h1",50, 0.0, 50.0)
    chain.Draw("fabs(%s) >> h1"%todraw,den)
    return h1

def getOverallEff(Teff, xinit):
    hden = Teff.GetCopyTotalHisto()
    hnum = Teff.GetCopyPassedHisto()
    nbins = hden.GetXaxis().GetNbins()
    xbin = hden.FindBin(xinit)
    den = hden.Integral(xbin, nbins)
    num = hnum.Integral(xbin, nbins)
    eff = 0.0; efferr  = 0.0
    if den == 0:
	    print "xinit ",xinit, " den ",den," num ",num," eff ",eff
    else:
	    eff = float(num)/float(den) 
	    efferr = sqrt(eff*(1-eff)/den)
    return eff, efferr

def getPartOfHist(hist, xinit, xend):
    nbins = hist.GetXaxis().GetNbins()
    xbin1 = hist.FindBin(xinit)
    xbin2 = hist.FindBin(xend)
    hist2 = hist.Clone()
    for x in range(0, xbin1):
	hist2.SetBinContent(x, 0.0)
	hist2.SetBinError(x, 0.0)

    for x in range(xbin2, nbins):
	hist2.SetBinContent(x, 0.0)
	hist2.SetBinError(x, 0.0)
    return hist2

def getPartOfEff(Teff, xinit, xend):
    hden = Teff.GetCopyTotalHisto()
    hnum = Teff.GetCopyPassedHisto()
    hden2 =  getPartOfHist(hden, xinit, xend)
    hnum2 =  getPartOfHist(hnum, xinit, xend)
    return ROOT.TEfficiency(hnum2, hden2)


def maximizeSensitivity(efflist, ratelist, nsignals, fraction):
    fraction = float(fraction)/100.0
    signal_tot = 0.0
    for n in nsignals:
	signal_tot += n
    weightlist = []
    selectedPoints = [0,0,0,0]
    Eff = 0.0
    for npar in range(0, 4):
	weightlist.append(float(nsignals[npar])/float(signal_tot))
    def getEffAndRate():
	efftmp = 0.0
	ratetmp = 0
	for npar in range(0, 4):
	    x = selectedPoints[npar]
	    efftmp += efflist[npar][x]*weightlist[npar]
	    ratetmp += ratelist[npar][x]
	return efftmp, ratetmp
    def getNextPoint():
	nextpar = -1
	drate = 999999
	for  npar in range(0, 4):
	    x = selectedPoints[npar]
	    if x>= len(ratelist[npar])-1:
		continue#impossible value
	    if ratelist[npar][x+1] < drate:
		nextpar = npar
		drate = ratelist[npar][x+1]
	return nextpar

    Eff,Nrate = getEffAndRate()
    #print "ratelist ",ratelist, " efflist ",efflist, " weightlist ",weightlist," nsignals ",nsignals
    #print "initial Eff ",Eff, " Nrate ",Nrate
    while (Eff<fraction):
	nextpar = getNextPoint()
	if nextpar < 0 :
	    break
	selectedPoints[nextpar] += 1
	Eff,Nrate = getEffAndRate()		
    ##varying final selected points
    #print "final Eff ",Eff, " Rate ",Nrate, " selectedPoints ",selectedPoints
    return Eff,Nrate,selectedPoints
		
def getEvents(ch, cut):
    hist = TH2F("hist","hist", 1,-1, 1, 1,-1, 1)
    ch.Draw("1:1>>hist",cut)
    return hist.GetEntries()


def getCut(ch, todraw, basecut, dphimin, dphimax, step, fraction):
    frac = 1.0
    total = getEvents(ch, basecut)
    naccept = getEvents(ch, basecut+"&& abs(%s)<%f"%(todraw, dphimax))
    naccept2 = getEvents(ch, basecut+"&& abs(%s)<%f"%(todraw, dphimin))
    #print "ch ",ch, " todraw ",todraw, " basecut ",basecut," tot ",total," naccept ",naccept," naccept2 ",naccept2
    frac = float(naccept)/float(total)
    frac2 = float(naccept2)/float(total)
    #print "basecut ",basecut," low side ",frac2," up side ",frac
    if frac<fraction:
	print "error! , with init range max ",dphimax," efficiency is ",frac," required eff: ",fraction
	return dphimax
    elif frac2>fraction:
	print "error! , with init range min ",dphimin," efficiency is ",frac2," required eff: ",fraction
	return dphimin
    xcut = (dphimin+dphimax)/2.0
    naccept = getEvents(ch, basecut+"&& abs(%s)<%f"%(todraw, xcut))
    frac = float(naccept)/float(total)
    #print "naccept ",naccept," frac ",frac, " xcut ",xcut," dphimin ",dphimin, " dphimax ",dphimax," required eff: ",fraction
    while(dphimax-dphimin > step):
	if frac>fraction:
	    dphimax = xcut
	elif frac<fraction:
	    dphimin = xcut
	else:
	    break
	xcut = (dphimin+dphimax)/2.0
	naccept = getEvents(ch, basecut+"&& abs(%s)<%f"%(todraw, xcut))
	frac = float(naccept)/float(total)
        #print "naccept ",naccept, " frac ",frac, " xcut ",xcut," dphimin ",dphimin, " dphimax ",dphimax," required eff: ",fraction
    return xcut

    
def getFinalLUTValue(pt, fraction, etamin, etamax, npar, suffix, filename):
    algos_lut = {"Position": -1, "Direction"+suffix: -1,"Hybrid"+suffix:[]}
    readinLUT = open(filename,'r')
    startread = [False, False, False]
    endread = [False, False, False]
    
    #special case: Directionbased: eta:2.1
    if etamax >= 2.05 and etamax < 2.2:
	    etamax = 2.1
    if etamin<2.1 and etamin>=2.0:
	    etamin = 2.0
    if etamin>=2.1 and etamin<2.2:
	    etamin = 2.1
    if etamin>= 1.6 and etamin<1.8:
	    etamin = 1.65
    inpar = 0
    for line in readinLUT:
	#Positioneta16to18pt5fraction76  
	getall = (endread[0] and endread[1] and endread[2])
	if getall:
	    break
	ialgo = -1
	for algo in algos_lut:
	    ialgo += 1
	    #print "ialgo ",ialgo, " algo ",algo
	    if line.startswith("%seta%dto%dpt%dfraction%d"%(algo, int(etamin*10), int(etamax*10),pt, fraction)):
		startread[ialgo] = True
		continue
	    if startread[ialgo] and line=="}\n":
		endread[ialgo] = True
		continue
	    if line=="{\n":
		continue
	    if startread[ialgo] and not(endread[ialgo]):
		if algo=="Position" or algo[:9]=="Direction":
		    allnums = line.split(",")
		    algos_lut[algo] = float(allnums[npar])
		    endread[ialgo] = True
		elif algo[:6]=="Hybrid" and inpar == npar:
		    allnums = line[1:-3].split(",")
		    for ipara in range(0,5):
			algos_lut[algo].append(float(allnums[ipara]))
		    endread[ialgo] = True
		elif algo[:6]=="Hybrid":
		    inpar += 1
	  
    #print "all lut for pt ",pt," fraction ",fraction, " etamin ", etamin, " etamax ",etamax, " npar ",npar
    #for key in algos_lut:
    #print "algo ",key, " value ",algos_lut[key]
    return algos_lut

def getDirectiondPhiME0LUTValue(pt, fraction, filename):
    readinLUT = open(filename,'r')
    startread = False; endread = False

    for line in readinLUT:
	if endread:
		break
	#ME0pt5fraction96
	if  line.startswith("ME0pt%dfraction%d"%(pt, fraction)):
		startread = True
		continue
	if startread and line=="}\n":
		endread = True
		continue
	if line=="{\n":
		continue
	if startread and not(endread):
		#print "line ",line ," float ",float(line)
		return float(line)

		    
#_____________________________________________________________
def makeEffplot_eta(chainlist,todraw, den, num, etamin, etamax, xtitle,ytitle,leg, legheader, txt,picname):
	
    b1 = ROOT.TH1F("b1","b1", 50, 0.0, 2.5)
    b1.GetYaxis().SetRangeUser(0.0,1.05)
    b1.GetYaxis().SetTitleOffset(1.1)
    b1.GetYaxis().SetNdivisions(520)
    b1.GetYaxis().SetTitle(xtitle)
    b1.GetXaxis().SetTitle(ytitle)
    b1.GetXaxis().SetTitleSize(0.05)
    b1.GetXaxis().SetLabelSize(0.05)
    b1.GetYaxis().SetTitleSize(0.05)
    b1.GetYaxis().SetLabelSize(0.05)
    #b1.SetTitle(" #scale[1.4]{#font[61]{CMS}} #font[52]{Simulation preliminary}"+"  "*8+" 14 TeV, 200 PU")
    b1.SetTitle(" #scale[1.4]{#font[61]{CMS}} Phase-2 Simulation"+"  "*4+" #sqrt{s}=14 TeV, <PU>=200")
    
    
    Teffs = []
    hdens = []
    hnums = []
    npars = 4
    Upperlimits = [0 for x in range(npars)]
    #print "chainlist ",chainlist
    
    if (etamin<1.6):
	    ring = 2
    if (etamin>=1.6):
	    ring = 1
    #legend = ROOT.TLegend(0.45,0.15,0.9,0.5)
    legend = ROOT.TLegend(0.45,0.20,0.8,0.50)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetTextFont(42)
    legend.SetTextSize(.05)
    legend.SetHeader("%s"%legheader)
    muons = ["Prompt Muons","Displaced Muons, 5<|d_{xy}|<50cm"]
    for n in range(len(chainlist)):
	#	print "n ",n," den ",den[n]," num ",num[n]
	hden = ROOT.TH1F("hden%d"%n,"hden%d"%n,50, 0.0, 2.5)
	hnum = ROOT.TH1F("hnum%d"%n,"hnum%d"%n,50, 0.0, 2.5)
	for i in range(0, 50):
		hden.SetBinContent(i, 0)
		hnum.SetBinContent(i,.0)
	#Teffs.append(ROOT.TEfficiency(hnum, hden))
	hdens.append(hden)
	hnums.append(hnum)
	chain = chainlist[n]
	#print "den ",den[n]," todraw ",todraw
	hdens[n].Add(gethist1DEta(chain, den[n],todraw))
	hnums[n].Add(gethist1DEta(chain, den[n]+" && %s"%(num[n]), todraw))

    c1 = ROOT.TCanvas()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetTickx()
    c1.SetTicky()

    b1.SetStats(0)
    b1.Draw()
    for n in range(len(chainlist)):
	Teffs.append(ROOT.TEfficiency(hnums[n],hdens[n]))	
	SetOwnership(Teffs[n], False)
	m = len(Teffs)-1
	Teffs[m].SetLineColor(color[n])
	Teffs[m].SetMarkerColor(color[n])
	Teffs[m].SetMarkerStyle(maker[n])
	Teffs[m].Draw("same")
	legend.AddEntry(Teffs[n],"%s"%leg[n],"pl")
    #print "Teffs ",Teffs
    #Teffs[0].Print("ALL")
    legend.Draw("same")

    tex = ROOT.TLatex(0.35,0.57,"%s"%txt)
    #tex = ROOT.TLatex(0.45,0.57,"#splitline{%s}{%d%% eff at %d [GeV]}"%(txt,fractionToKeep,pt))
    #tex = ROOT.TLatex(0.45,0.57,"#splitline{%s}{check the sign of #Delta Y_{12} and #Delta Y_{23}}"%(txt))
    tex.SetTextSize(0.05)
    tex.SetTextFont(62)
    tex.SetNDC()
    tex.Draw("same")
    #c1.Update()
    c1.SaveAs("%s.png"%(picname))
    c1.SaveAs("%s.pdf"%(picname))
    c1.SaveAs("%s.C"%(picname))
    return Teffs
	


#_____________________________________________________________
def makeEffplot_pt(chainlist,todraw, den, num, etamin, etamax, xtitle,ytitle,leg, legheader, txt,picname):
	
    b1 = ROOT.TH1F("b1","b1",len(ptbins)-1, ptbins)
    b1.GetYaxis().SetRangeUser(0.0,1.05)
    b1.GetYaxis().SetTitleOffset(1.1)
    b1.GetYaxis().SetNdivisions(520)
    b1.GetYaxis().SetTitle("Trigger efficiency")
    b1.GetXaxis().SetTitle("True muon p_{T} [GeV]")
    b1.GetXaxis().SetTitleSize(0.05)
    b1.GetXaxis().SetLabelSize(0.05)
    b1.GetYaxis().SetTitleSize(0.05)
    b1.GetYaxis().SetLabelSize(0.05)
    #b1.SetTitle(" #scale[1.4]{#font[61]{CMS}} #font[52]{Simulation preliminary}"+"  "*8+" 14 TeV, 200 PU")
    b1.SetTitle(" #scale[1.4]{#font[61]{CMS}} Phase-2 Simulation"+"  "*4+" #sqrt{s}=14 TeV, <PU>=200")
    
    
    Teffs = []
    hdens = []
    hnums = []
    npars = 4
    Upperlimits = [0 for x in range(npars)]
    #print "chainlist ",chainlist
    
    if (etamin<1.6):
	    ring = 2
    if (etamin>=1.6):
	    ring = 1
    #legend = ROOT.TLegend(0.45,0.15,0.9,0.5)
    legend = ROOT.TLegend(0.45,0.20,0.8,0.50)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetTextFont(42)
    legend.SetTextSize(.05)
    legend.SetHeader("%s"%legheader)
    muons = ["Prompt Muons","Displaced Muons, 5<|d_{xy}|<50cm"]
    for n in range(len(chainlist)):
	#	print "n ",n," den ",den[n]," num ",num[n]
	hden = ROOT.TH1F("hden%d"%n,"hden%d"%n,len(ptbins)-1, ptbins)
	hnum = ROOT.TH1F("hnum%d"%n,"hnum%d"%n,len(ptbins)-1, ptbins)
	for i in range(len(ptbins)):
		hden.SetBinContent(i, 0)
		hnum.SetBinContent(i,.0)
	#Teffs.append(ROOT.TEfficiency(hnum, hden))
	hdens.append(hden)
	hnums.append(hnum)
	chain = chainlist[n]
	#print "den ",den[n]
	hdens[n].Add(gethist1D(chain, den[n],todraw))
	hnums[n].Add(gethist1D(chain, den[n]+" && %s"%(num[n]), todraw))

    c1 = ROOT.TCanvas()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetTickx()
    c1.SetTicky()

    b1.SetStats(0)
    b1.Draw()
    for n in range(len(chainlist)):
	Teffs.append(ROOT.TEfficiency(hnums[n],hdens[n]))	
	SetOwnership(Teffs[n], False)
	m = len(Teffs)-1
	Teffs[m].SetLineColor(color[n])
	Teffs[m].SetLineWidth(2)
	
	Teffs[m].SetMarkerColor(color[n])
	Teffs[m].SetMarkerStyle(maker[n])
	Teffs[m].Draw("sameZ")
	legend.AddEntry(Teffs[n],"%s"%leg[n],"pl")
    #print "Teffs ",Teffs
    #Teffs[0].Print("ALL")
    legend.Draw("same")

    tex = ROOT.TLatex(0.35,0.57,"%s"%txt)
    #tex = ROOT.TLatex(0.45,0.57,"#splitline{%s}{%d%% eff at %d [GeV]}"%(txt,fractionToKeep,pt))
    #tex = ROOT.TLatex(0.45,0.57,"#splitline{%s}{check the sign of #Delta Y_{12} and #Delta Y_{23}}"%(txt))
    tex.SetTextSize(0.05)
    tex.SetTextFont(62)
    tex.SetNDC()
    tex.Draw("same")
    #c1.Update()
    c1.SaveAs("%s.png"%(picname))
    c1.SaveAs("%s.pdf"%(picname))
    c1.SaveAs("%s.C"%(picname))
    return Teffs



#____________________________________________________________________________________________
def makeEffplot_dxy(chainlist,todraw, den, num, etamin, etamax, xtitle,ytitle,leg, legheader, txt,picname):
	
    b1 = ROOT.TH1F("b1","b1",50, 0.0, 50.0)
    b1.GetYaxis().SetRangeUser(0.0,1.05)
    b1.GetYaxis().SetTitleOffset(1.1)
    b1.GetYaxis().SetNdivisions(520)
    b1.GetYaxis().SetTitle("Trigger efficiency")
    b1.GetXaxis().SetTitle("|dxy| [cm]")
    b1.GetXaxis().SetTitleSize(0.05)
    b1.GetXaxis().SetLabelSize(0.05)
    b1.GetYaxis().SetTitleSize(0.05)
    b1.GetYaxis().SetLabelSize(0.05)
    #b1.SetTitle(" #scale[1.4]{#font[61]{CMS}} #font[52]{Simulation preliminary}"+"  "*8+" 14 TeV, 200 PU")
    b1.SetTitle(" #scale[1.4]{#font[61]{CMS}} Phase-2 Simulation"+"  "*4+" #sqrt{s}=14 TeV, <PU>=200")
    
    
    Teffs = []
    hdens = []
    hnums = []
    npars = 4
    Upperlimits = [0 for x in range(npars)]
    #print "chainlist ",chainlist
    
    if (etamin<1.6):
	    ring = 2
    if (etamin>=1.6):
	    ring = 1
    #legend = ROOT.TLegend(0.45,0.15,0.9,0.5)
    legend = ROOT.TLegend(0.45,0.20,0.8,0.50)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetTextFont(42)
    legend.SetTextSize(.05)
    legend.SetHeader("%s"%legheader)
    muons = ["Prompt Muons","Displaced Muons, 5<|d_{xy}|<50cm"]
    for n in range(len(chainlist)):
	#	print "n ",n," den ",den[n]," num ",num[n]
	hden = ROOT.TH1F("hden%d"%n,"hden%d"%n, 50, 0.0, 50.0)
	hnum = ROOT.TH1F("hnum%d"%n,"hnum%d"%n, 50, 0.0, 50.0)
	for i in range(50):
		hden.SetBinContent(i, 0)
		hnum.SetBinContent(i,.0)
	#Teffs.append(ROOT.TEfficiency(hnum, hden))
	hdens.append(hden)
	hnums.append(hnum)
	chain = chainlist[n]
	#print "den ",den[n]
	hdens[n].Add(gethist1DDxy(chain, den[n],todraw))
	hnums[n].Add(gethist1DDxy(chain, den[n]+" && %s"%(num[n]), todraw))

    c1 = ROOT.TCanvas()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetTickx()
    c1.SetTicky()

    b1.SetStats(0)
    b1.Draw()
    for n in range(len(chainlist)):
	Teffs.append(ROOT.TEfficiency(hnums[n],hdens[n]))	
	SetOwnership(Teffs[n], False)
	m = len(Teffs)-1
	Teffs[m].SetLineColor(color[n])
	Teffs[m].SetLineWidth(2)
	
	Teffs[m].SetMarkerColor(color[n])
	Teffs[m].SetMarkerStyle(maker[n])
	Teffs[m].Draw("sameZ")
	legend.AddEntry(Teffs[n],"%s"%leg[n],"pl")
    #print "Teffs ",Teffs
    #Teffs[0].Print("ALL")
    legend.Draw("same")

    tex = ROOT.TLatex(0.35,0.57,"%s"%txt)
    #tex = ROOT.TLatex(0.45,0.57,"#splitline{%s}{%d%% eff at %d [GeV]}"%(txt,fractionToKeep,pt))
    #tex = ROOT.TLatex(0.45,0.57,"#splitline{%s}{check the sign of #Delta Y_{12} and #Delta Y_{23}}"%(txt))
    tex.SetTextSize(0.05)
    tex.SetTextFont(62)
    tex.SetNDC()
    tex.Draw("same")
    #c1.Update()
    c1.SaveAs("%s.png"%(picname))
    c1.SaveAs("%s.pdf"%(picname))
    c1.SaveAs("%s.C"%(picname))
    return Teffs


#_____________________________________________________________
def makeplotsEta(Teffs, legs, legheader, text, picname):
	
    color = [ROOT.kBlue,  ROOT.kBlack, ROOT.kRed, ROOT.kGreen+2, kBlack]
    marker = [25, 21, 26, 20, 24]
    b1 = ROOT.TH1F("b1","b1",18, 1.6, 2.5)
    b1.GetYaxis().SetRangeUser(0.0,1.05)
    b1.GetYaxis().SetTitleOffset(1.1)
    b1.GetYaxis().SetNdivisions(520)
    b1.GetYaxis().SetTitle("Trigger efficiency")
    b1.GetXaxis().SetTitle("|#eta| at 2nd CSC station")
    b1.GetXaxis().SetTitleSize(0.05)
    b1.GetXaxis().SetLabelSize(0.05)
    b1.GetYaxis().SetTitleSize(0.05)
    b1.GetYaxis().SetLabelSize(0.05)
    #b1.SetTitle(" #scale[1.4]{#font[61]{CMS}} #font[52]{Simulation preliminary}"+"  "*8+" 14 TeV, 200 PU")
    b1.SetTitle(" #scale[1.4]{#font[61]{CMS}} Phase-2 Simulation"+"  "*4+" #sqrt{s}=14 TeV, <PU>=200")

    c1 = ROOT.TCanvas()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetTickx()
    c1.SetTicky()

    dy_leg = 0.06*(len(Teffs)+1)
    
    legend = ROOT.TLegend(0.4,0.2,0.88, 0.2+dy_leg)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetTextFont(42)
    #legend.SetTextSize(.05)
    legend.SetHeader("%s"%legheader)
    b1.SetStats(0)
    b1.Draw()
    for n in range(len(Teffs)):
	Teffs[n].SetLineColor(color[n])
	Teffs[n].SetMarkerColor(color[n])
	Teffs[n].SetLineWidth(2)
	Teffs[n].SetMarkerStyle(marker[n])
	Teffs[n].Draw("samepZ")
	legend.AddEntry(Teffs[n],"%s"%legs[n],"pl")
    legend.Draw("same")

    tex = ROOT.TLatex(0.4,0.65,"%s"%text)
    tex.SetTextSize(0.04)
    tex.SetTextFont(62)
    tex.SetNDC()
    tex.Draw("same")
    #c1.Update()
    c1.SaveAs("%s.png"%(picname))
    #c1.SaveAs("%s.eps"%(picname))
    c1.SaveAs("%s.pdf"%(picname))
    c1.SaveAs("%s.C"%(picname))


def makeplotsPt(Teffs, legs, legheader, text, picname):

    color = [ROOT.kBlue, ROOT.kBlack, ROOT.kRed, ROOT.kGreen+2]
    marker = [25, 21, 26, 20, 24]

    b1 = ROOT.TH1F("b1","b1",len(ptbins)-1, ptbins)
    b1.GetYaxis().SetRangeUser(0.0,1.05)
    b1.GetYaxis().SetTitleOffset(1.1)
    b1.GetYaxis().SetNdivisions(520)
    b1.GetYaxis().SetTitle("Trigger efficiency")
    b1.GetXaxis().SetTitle("True muon p_{T} [GeV]")
    b1.GetXaxis().SetTitleSize(0.05)
    b1.GetXaxis().SetLabelSize(0.05)
    b1.GetYaxis().SetTitleSize(0.05)
    b1.GetYaxis().SetLabelSize(0.05)
    #b1.SetTitle(" #scale[1.4]{#font[61]{CMS}} #font[52]{Simulation preliminary}"+"  "*8+" 14 TeV, 200 PU")
    b1.SetTitle(" #scale[1.4]{#font[61]{CMS}} Phase-2 Simulation"+"  "*4+" #sqrt{s}=14 TeV, <PU>=200")

    c1 = ROOT.TCanvas()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetTickx()
    c1.SetTicky()

    dy_leg = 0.06*(len(Teffs)+1)
    
    legend = ROOT.TLegend(0.45,0.2,0.9,0.2+dy_leg)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetTextFont(42)
    #legend.SetTextSize(.05)
    legend.SetHeader("%s"%legheader)
    b1.SetStats(0)
    b1.Draw()
    for n in range(len(Teffs)):
	Teffs[n].SetLineColor(color[n])
	Teffs[n].SetLineWidth(2)
	Teffs[n].SetMarkerColor(color[n])
	Teffs[n].SetMarkerStyle(marker[n])
	Teffs[n].Draw("samepZ")
	legend.AddEntry(Teffs[n],"%s"%legs[n],"pl")
    legend.Draw("same")

    tex = ROOT.TLatex(0.4,0.6,"%s"%text)
    tex.SetTextSize(0.04)
    tex.SetTextFont(62)
    tex.SetNDC()
    tex.Draw("same")
    #c1.Update()
    c1.SaveAs("%s.png"%(picname))
    c1.SaveAs("%s.pdf"%(picname))
    c1.SaveAs("%s.C"%(picname))

def makeplotsDxy(Teffs, legs, legheader, text, picname):

    color = [ROOT.kBlue, ROOT.kBlack, ROOT.kRed, ROOT.kGreen+2]
    marker = [25, 21, 26, 20, 24]

    b1 = ROOT.TH1F("b1","b1", 50, 0.0, 50.0)
    b1.GetYaxis().SetRangeUser(0.0,1.05)
    b1.GetYaxis().SetTitleOffset(1.1)
    b1.GetYaxis().SetNdivisions(520)
    b1.GetYaxis().SetTitle("Trigger efficiency")
    b1.GetXaxis().SetTitle("|dxy| [cm]")
    b1.GetXaxis().SetTitleSize(0.05)
    b1.GetXaxis().SetLabelSize(0.05)
    b1.GetYaxis().SetTitleSize(0.05)
    b1.GetYaxis().SetLabelSize(0.05)
    #b1.SetTitle(" #scale[1.4]{#font[61]{CMS}} #font[52]{Simulation preliminary}"+"  "*8+" 14 TeV, 200 PU")
    b1.SetTitle(" #scale[1.4]{#font[61]{CMS}} Phase-2 Simulation"+"  "*4+" #sqrt{s}=14 TeV, <PU>=200")

    c1 = ROOT.TCanvas()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetTickx()
    c1.SetTicky()

    dy_leg = 0.06*(len(Teffs)+1)
    
    legend = ROOT.TLegend(0.45,0.2,0.9,0.2+dy_leg)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetTextFont(42)
    #legend.SetTextSize(.05)
    legend.SetHeader("%s"%legheader)
    b1.SetStats(0)
    b1.Draw()
    for n in range(len(Teffs)):
	Teffs[n].SetLineColor(color[n])
	Teffs[n].SetLineWidth(2)
	Teffs[n].SetMarkerColor(color[n])
	Teffs[n].SetMarkerStyle(marker[n])
	Teffs[n].Draw("samepZ")
	legend.AddEntry(Teffs[n],"%s"%legs[n],"pl")
    legend.Draw("same")

    tex = ROOT.TLatex(0.4,0.6,"%s"%text)
    tex.SetTextSize(0.04)
    tex.SetTextFont(62)
    tex.SetNDC()
    tex.Draw("same")
    #c1.Update()
    c1.SaveAs("%s.png"%(picname))
    c1.SaveAs("%s.pdf"%(picname))
    c1.SaveAs("%s.C"%(picname))


def plotshistograms(hists, legs,legheader, text, picname, RateVsPt=False):

    color = [ROOT.kBlue, ROOT.kBlack,  ROOT.kRed, ROOT.kGreen+2, kBlack]
    marker = [25, 21, 26, 20, 24]
    #hs = ROOT.THStack("hs"," #scale[1.4]{#font[61]{CMS}} #font[52]{Simulation preliminary}"+"  "*8+" 14 TeV,%d PU"%pu)
    hs = ROOT.THStack("hs"," #scale[1.4]{#font[61]{CMS}} Phase-2 Simulation"+"  "*4+" #sqrt{s}=14 TeV, <PU>=200")
    c1 = ROOT.TCanvas()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetTickx()
    c1.SetTicky()
    c1.SetLogy()
    if RateVsPt:
	hs.SetMinimum(1)
	hs.SetMaximum(2000.)
    dy_leg = 0.04*(len(hists)+1)

    legend = ROOT.TLegend(0.46,0.45,0.88,0.450+dy_leg)
    legend.SetHeader(legheader)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetTextFont(42)
    #legend.SetTextSize(.05)
    nhist = -1
    for hist in hists:
	nhist +=1
	hist.SetMarkerColor(color[nhist])
	hist.SetLineColor(color[nhist])
	hist.SetLineWidth(2)
	hist.SetMarkerStyle(marker[nhist])
	#for rate
	#hist.Scale(SF)
	hs.Add(hist)
	#legend.AddEntry(hist, legs[nhist]+",event: %d"%(hist.GetEntries()),"pl")
	legend.AddEntry(hist, legs[nhist],"pl")
    hs.Draw("nostacke")
    legend.Draw("same")
    if RateVsPt:
	c1.SetLogy()
	#c2.SetLogx()
	c1.SetGridx()
	c1.SetGridy()
	
    c1.Update()
    c1.Modified()
    #print "hists[0] ",hists[0]
    #hs.GetHistogram().GetXaxis().SetTitle(hists[0].GetXaxis().GetTitle())
    #hs.GetHistogram().GetYaxis().SetTitle(hists[0].GetYaxis().GetTitle())
    #hs.GetHistogram().GetXaxis().SetTitle("L1Mu |#eta|")
    hs.GetHistogram().GetXaxis().SetTitle("|#eta| at 2nd CSC station")
    texh=0.2
    if RateVsPt:
	hs.GetHistogram().GetXaxis().SetTitle("Trigger p_{T} threshold [GeV]")
	texh = 0.2
    hs.GetHistogram().GetYaxis().SetTitle("Trigger rate [kHz]")
	
    tex = ROOT.TLatex(0.2,texh,"%s"%text)
    tex.SetTextSize(0.045)
    tex.SetTextFont(62)
    tex.SetNDC()
    tex.Draw("same")
    #c1.Update()
    c1.SaveAs("%s.png"%(picname))
    c1.SaveAs("%s.pdf"%(picname))
    c1.SaveAs("%s.eps"%(picname))
    c1.SaveAs("%s.C"%(picname))


def plotsgraphs(graphs, legs, legheader, text, picname):

    c1 = TCanvas()
    c1.SetGridx()
    c1.SetGridy()
    mgerr =TMultiGraph()
    #print "size ", len(efflist_p[npar]), " ",efflist_p[npar]
    dy_leg = 0.06*(len(graphs)+1)
    legend = ROOT.TLegend(0.15,0.65,0.5,0.65+dy_leg)
    legend.SetHeader(legheader)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetTextFont(42)
    #legend.SetTextSize(.05)
    n= -1
    for graph in graphs:
	n += 1
	graph.SetMarkerStyle(maker[n])
	graph.SetMarkerColor(color[n])
	graph.SetLineColor(color[n])
	graph.SetFillColor(color[n])
    	legend.AddEntry(graph, legs[n], "p")
	mgerr.Add(graph)

    mgerr.Draw("AP")
    legend.Draw("same")
    mgerr.SetMinimum(0)
    #mgerr.SetTitle(" #scale[1.4]{#font[61]{CMS}} #font[52]{Simulation preliminary}"+"  "*10+" 14 TeV, %d PU"%pu)
    mgerr.SetTitle(" #scale[1.4]{#font[61]{CMS}} Phase-2 Simulation"+"  "*4+" #sqrt{s}=14 TeV, <PU>=200")
    mgerr.GetXaxis().SetTitle("Trigger plateau efficiency")
    mgerr.GetYaxis().SetTitle("Trigger rate [kHz]")

    tex = ROOT.TLatex(0.2,0.6,"%s"%text)
    tex.SetTextSize(0.045)
    tex.SetTextFont(62)
    tex.SetNDC()
    tex.Draw("same")
    #c1.Update()
    c1.Update()
    c1.SaveAs("%s.png"%(picname))
    c1.SaveAs("%s.pdf"%(picname))
    c1.SaveAs("%s.C"%(picname))
    return mgerr


def getRateVsEta(L1MuTriggerRate, histname, myetabin, cut, vetotype):
    if vetotype == 1:
	    cut = cut + "&& !L1Mu_looseveto"
    elif vetotype == 2:
	    cut = cut + "&& !L1Mu_mediumveto"
    elif vetotype == 3:
	    cut = cut + "&& !L1Mu_tightveto"
    hist = TH1F(histname, histname, len(myetabin)-1, myetabin)
    L1MuTriggerRate.Draw("abs(L1Mu_eta_st2)>>"+histname, cut)
    hist.GetXaxis().SetTitle("L1Mu |#eta|")
    hist.GetYaxis().SetTitle("Trigger Rate [kHz]")
    hist.Sumw2()
    print "SF ",SF," cut ",cut
    hist.Scale(SF)
    return hist


def getEffAndRate(pt, fractions, etamin, etamax, wfilename, filedir, condoroutroot, SBPlots, useGE21, doPlots=True, doLUTs=False, doRates=True):#uncertainty?

    efflist_p = []
    ratelist_p = []
    effErrlist_p = []
    rateErrlist_p = []
    efflist_d = []
    ratelist_d = []
    effErrlist_d = []
    rateErrlist_d = []
    efflist_h = []
    ratelist_h = []
    effErrlist_h = []
    rateErrlist_h = []
    nsignals = []

    extrasuffix = ""
    extraleg = "with GE21"
    if not(useGE21):
	extrasuffix = "_noGE21"
	extraleg = "w/o GE21"
    if etamin<1.6:
	extraleg = ""

    #nstarts = {1.8 : 1, 2.0: 2}
    #nstarts = {1.2 : 1, 1.4: 2}
    nstarts = {1.6: 1, 1.8 : 2, 2.0: 3, 1.2 : 1, 1.4: 2}
    for npar in range(0, 5):
	efflist_p.append([])
	ratelist_p.append([])
	effErrlist_p.append([])
	rateErrlist_p.append([])
	efflist_d.append([])
	ratelist_d.append([])
	effErrlist_d.append([])
	rateErrlist_d.append([])
	efflist_h.append([])
	ratelist_h.append([])
	effErrlist_h.append([])
	rateErrlist_h.append([])

    for frac in fractions:
	filename  = filedir+"EventsNum_Pt_%d_fraction_%d.log"%(pt, frac)
	if not(os.path.isfile(filename)):
	    print "failed to open file ",filename
	    continue
	fileopen = open(filename,'r')
	npar = 0
	startreading = False
	endreading = False
	istart = 0
	iend = 0
	print "fileopen ",fileopen, " nstarts ",nstarts[etamin]," etamin ",etamin
	for line in fileopen:
	    #v1: {284900.000000 1501.000000 64.000000 286.000000 982.000000 877.000000 876.000000 42.000000 840.000000
	    #v2: 284900.000000 1641.000000 29.000000 167.000000 18.000000 19122.000000 14524.000000 14528.000000 14601.000000  
	    #print "line ",line," istart ",istart," iend ",iend
	    if line.startswith("{"):
		istart += 1
		continue
	    if line.startswith("},"):
		iend += 1
		continue
	    if istart >= nstarts[etamin]:
		startreading = True
	    if iend >= (nstarts[etamin]):
		endreading = True
	    #print "line ",line," istart ",istart," iend ",iend," starreading ",startreading," endreading ",endreading
	    if (not (startreading and not(endreading))):
		continue
	    allnums = line.split(" ")
	    #print "filename ", filename, "all nums ", allnums
	    n_signal_tot = float(allnums[5])
	    n_signal_p = float(allnums[6])
	    n_signal_d = float(allnums[7])
	    n_signal_h = float(allnums[8])
	    n_rate_tot = float(allnums[1])
	    n_rate_p = float(allnums[2])*SF
	    n_rate_d = float(allnums[3])*SF
	    n_rate_h = float(allnums[4])*SF
	    eff_p = n_signal_p/n_signal_tot
	    eff_d = n_signal_d/n_signal_tot
	    eff_h = n_signal_h/n_signal_tot
	    effErr_p = eff_p * sqrt((1-eff_p)/n_signal_p)
	    effErr_d = eff_d * sqrt((1-eff_d)/n_signal_p)
	    effErr_h = eff_h * sqrt((1.0-eff_h)/n_signal_p)
	    RateErr_p = sqrt(float(allnums[2]))*SF
	    RateErr_d = sqrt(float(allnums[3]))*SF
	    RateErr_h = sqrt(float(allnums[4]))*SF
	    
	    if (len(nsignals)<4):
		nsignals.append(n_signal_tot)

	    efflist_p[npar].append(eff_p)
	    ratelist_p[npar].append(n_rate_p)
	    effErrlist_p[npar].append(effErr_p)
	    rateErrlist_p[npar].append(RateErr_p)
	    efflist_d[npar].append(eff_d)
	    ratelist_d[npar].append(n_rate_d)
	    effErrlist_d[npar].append(effErr_d)
	    rateErrlist_d[npar].append(RateErr_d)
	    efflist_h[npar].append(eff_h)
	    ratelist_h[npar].append(n_rate_h)
	    effErrlist_h[npar].append(effErr_h)
	    rateErrlist_h[npar].append(RateErr_h)
	    npar += 1 

    totsignals = 0.0
    for n in nsignals:
	totsignals += n
    #print "nsignals ",nsignals
    rfilename = filedir+condoroutroot
    rfile = TFile(rfilename, "READ")
    print "filename ",filename
    def getCombinedEffandRatePlot(algo, selectedPoint, fractions, thisfrac):
	hratelist = []
	Tefflistprompt = []
	Tefflistdisplaced = []
	for npar in range(0, 4):
	    ncount = -1
	    fraction = 0
	    for frac in fractions: 
		filename  = filedir+"EventsNum_Pt_%d_fraction_%d.log"%(pt, frac)
		if os.path.isfile(filename):
		    ncount += 1
		if ncount == selectedPoint[npar]:
		    fraction = frac
		    break
	    if ncount == -1:
		print "failed to find point ",selectedPoint[npar]," in par ",npar," pt ",pt
	    #ratehist_hybrid_eta20to21_npar1_frac96_pt7;1
	    #hybrideta18to20npar2promptmuonpt7ptbg0fraction92;1
	    #hybrideta18to20npar3displacedmuonpt7ptbg0fraction92;1
	    algo2 = algo
	    if algo=="position":
		algo2 = "poistion"
	    histname = "ratehist_%s_eta%dto%d_npar%d_frac%d_pt%d"%(algo, int(etamin*10), int(etamax*10), npar, fraction, pt)
	    teffname1 = "%seta%dto%dnpar%dpromptmuonpt%dptbg0fraction%d"%(algo2, int(etamin*10), int(etamax*10), npar,pt, fraction)
	    teffname2 = "%seta%dto%dnpar%ddisplacedmuonpt%dptbg0fraction%d"%(algo2, int(etamin*10), int(etamax*10), npar, pt, fraction)
	    hratelist.append(rfile.Get(histname).Clone())
	    Tefflistprompt.append(rfile.Get(teffname1).Clone())
	    Tefflistdisplaced.append(rfile.Get(teffname2).Clone())
	    #getOverallEff(Tefflistprompt[npar], pt)
	    #getOverallEff(Tefflistdisplaced[npar], pt)
	    #print "selected rate ",histname
	for npar in range(1, 4):
	    #if npar == 0:
		#hist = hratelist[0].Clone()
	    #else:
	    #hist.Add(hratelist[npar])
	    hratelist[0].Add(hratelist[npar])
	    Tefflistprompt[0].Add(Tefflistprompt[npar])
	    Tefflistdisplaced[0].Add(Tefflistdisplaced[npar])
	#rfile.Close()
	hratelist[0].SetName("ratehist_%s_eta%dto%d_npar%d_frac%d_pt%d"%(algo, int(etamin*10), int(etamax*10), 4, thisfrac, pt)+extrasuffix)
	Tefflistprompt[0].SetName("Teff_prompt_%s_eta%dto%d_npar%d_frac%d_pt%d"%(algo, int(etamin*10), int(etamax*10), 4, thisfrac, pt)+extrasuffix)
	Tefflistdisplaced[0].SetName("Teff_displaced_%s_eta%dto%d_npar%d_frac%d_pt%d"%(algo, int(etamin*10), int(etamax*10), 4, thisfrac, pt)+extrasuffix)
	#SetOwnership(hratelist[0], False)
	#hratelist[0].Print("ALL")
	return Tefflistprompt[0], Tefflistdisplaced[0], hratelist[0]
	#return 
	    
    def getLUTsValue(algo, selectedPoint, fractions, lut, tag):
	cutlist = []
	for npar in range(0, 4):
	    #npoint = selectedPoint[npar]
	    ncount = -1
	    fraction = 0
	    for frac in fractions: 
		filename  = filedir+"EventsNum_Pt_%d_fraction_%d.log"%(pt, frac)
		if os.path.isfile(filename):
		    ncount += 1
		if ncount == selectedPoint[npar]:
		    fraction = frac
		    break
	    if ncount == -1:
		print "failed to find point ",selectedPoint[npar]," in par ",npar," pt ",pt
	    readinLUT = open(filedir+"%sLUT_Pt%d_fraction%d.log"%(algo, pt, fraction),'r')
	    startreading = False
	    endreading = False
	    istart = 0
	    iend = 0
	    iline = -1
	    print "fileopen ",readinLUT
	    for line in readinLUT:
		 #{284900.000000 1501.000000 64.000000 286.000000 982.000000 877.000000 876.000000 42.000000 840.000000
		#print "line ",line," istart ",istart," iend ",iend
		if line.startswith("{") and not(line.endswith(",\n")):
		    istart += 1
		    continue
		if line.startswith("},"):
		    iend += 1
		    continue
		if istart >= nstarts[etamin]:
		    startreading = True
		if iend >= nstarts[etamin]:
		    endreading = True
		print "line ",line," istart ",istart," iend ",iend," starreading ",startreading," endreading ",endreading," nstart ", nstarts[etamin]
		if (not (startreading and not(endreading))):
		    continue
		iline += 1
		if (algo == "Position" or algo == "Direction"):
		    allnums = line.split(",")
		    cutlist.append(allnums[npar])
		    break
		if algo == "Hybrid" and npar == iline:
		    cutlist.append(line)
		    break
	print "algo ",algo, " cuts ",cutlist
	lut.write(algo+extrasuffix+tag+"\n{\n")
	if algo == "Position" or algo == "Direction":
	    for x in cutlist[:-1]:
		lut.write(x+",")
	    lut.write(cutlist[-1]+"}\n")
	else:
	    lut.writelines(x for x in cutlist)
	    lut.write("}\n")
		#print "filename ", filename, "all nums ", allnums

    finalLUT = open("finalLUT_20170303.log",'a')
    hratelist_p = []
    hratelist_d = []
    hratelist_h = []
    nfrac = -1
    for frac in fractions:
	nfrac += 1
	print "pt ",pt," etamin ",etamin, "etamax ",etamax, " fraction ",frac
	peff,prate,pselected = maximizeSensitivity(efflist_p, ratelist_p, nsignals, frac)
	deff,drate,dselected = maximizeSensitivity(efflist_d, ratelist_d, nsignals, frac)
	heff,hrate,hselected = maximizeSensitivity(efflist_h, ratelist_h, nsignals, frac)
	print "from log file position eff ",peff, " prate ",prate ," hybrid eff ",heff," hrate ",hrate
	efflist_p[4].append(peff)
	ratelist_p[4].append(prate)
	effErrlist_p[4].append(sqrt((1-peff)*peff/totsignals))
	rateErrlist_p[4].append(sqrt(prate*SF))
	efflist_d[4].append(deff)
	ratelist_d[4].append(drate)
	effErrlist_d[4].append(sqrt((1-deff)*deff/totsignals))
	rateErrlist_d[4].append(sqrt(drate*SF))
	efflist_h[4].append(heff)
	ratelist_h[4].append(hrate)
	effErrlist_h[4].append(sqrt((1-heff)*heff/totsignals))
	rateErrlist_h[4].append(sqrt(hrate*SF))
	if doLUTs:
	    #finalLUT.write("eta%dto%dpt%dfraction%d"%(int(etamin*10), int(etamax*10), pt, frac)+"\n")	
	    tag = "eta%dto%dpt%dfraction%d"%(int(etamin*10), int(etamax*10), pt, frac)
	    #getLUTsValue("Position", pselected, fractions, finalLUT, tag)
	    getLUTsValue("Direction", dselected, fractions, finalLUT, tag)
	    getLUTsValue("Hybrid", hselected, fractions, finalLUT, tag)
	if doRates:
	    wfile = TFile(wfilename, "UPDATE") 
	    teff1_p, teff2_p, hrate_p = getCombinedEffandRatePlot("position", pselected, fractions, frac)
	    teff1_d, teff2_d, hrate_d = getCombinedEffandRatePlot("direction", dselected, fractions, frac)
	    teff1_h, teff2_h, hrate_h = getCombinedEffandRatePlot("hybrid", hselected, fractions, frac)
	    eff1_p, err1_p = getOverallEff(teff1_p, pt)
	    eff2_p, err2_p = getOverallEff(teff2_p, pt)
	    eff1_h, err1_h = getOverallEff(teff1_h, pt)
	    eff2_h, err2_h = getOverallEff(teff2_h, pt)
	    print "eff1 p ",eff1_p, " eff2_p ",eff2_p, " eff1_h ",eff1_h, " eff2_h ",eff2_h, " rate entries position ",hrate_p.GetEntries()," hyrid ",hrate_h.GetEntries()
	    teff1_p.Write();teff2_p.Write();hrate_p.Write()
	    teff1_d.Write();teff2_d.Write();hrate_h.Write()
	    teff1_h.Write();teff2_h.Write();hrate_d.Write()
	    wfile.Close()
	    
    wfile = TFile(wfilename, "UPDATE") 
    #print "total rate ",ratelist_h 
    legs = ["Position-based algorithm","Stub alignment algorithm %s"%extraleg,"Hybrid algorithm %s"%extraleg]

    if not(doPlots):	
	wfile.Close()
	return
    graphlist_p = []
    graphlist_d = []
    graphlist_h = []
    for npar in range(0, 5):
	legend = ROOT.TLegend(0.15,0.75,0.63,0.9)
	legend.SetFillColor(ROOT.kWhite)
	legend.SetTextFont(42)
	legend.SetTextSize(.05)
	c1 = TCanvas()
	mgerr =TMultiGraph()
	#print "size ", len(efflist_p[npar]), " ",efflist_p[npar]

	graphlist_p.append(TGraphErrors(len(efflist_p[npar]), np.array(efflist_p[npar]), np.array(ratelist_p[npar]), np.array(effErrlist_p[npar]), np.array(rateErrlist_p[npar])))
	graphlist_d.append(TGraphErrors(len(efflist_d[npar]), np.array(efflist_d[npar]), np.array(ratelist_d[npar]), np.array(effErrlist_d[npar]), np.array(rateErrlist_d[npar])))
	graphlist_h.append(TGraphErrors(len(efflist_h[npar]), np.array(efflist_h[npar]), np.array(ratelist_h[npar]), np.array(effErrlist_h[npar]), np.array(rateErrlist_h[npar])))
	graphlist_p[npar].SetMarkerStyle(maker[0])
	graphlist_p[npar].SetMarkerColor(color[0])
	#graphlist_p[npar].SetMarkerSize(.8)
	graphlist_d[npar].SetMarkerStyle(maker[1])
	graphlist_d[npar].SetMarkerColor(color[1])
	graphlist_h[npar].SetMarkerStyle(maker[2])
	graphlist_h[npar].SetMarkerColor(color[2])
	
	graphlist_p[npar].SetName("TriggerEffVsRate_position_pt%d_eta%dto%d_npar%d%s"%(pt, int(etamin*10), int(etamax*10), npar, extrasuffix))
	graphlist_d[npar].SetName("TriggerEffVsRate_direction_pt%d_eta%dto%d_npar%d%s"%(pt, int(etamin*10), int(etamax*10), npar, extrasuffix))
	graphlist_h[npar].SetName("TriggerEffVsRate_hybrid_pt%d_eta%dto%d_npar%d%s"%(pt, int(etamin*10), int(etamax*10), npar, extrasuffix))
	
	mgerr.Add(graphlist_p[npar])
	#mgerr.Add(graphlist_d[npar])
	mgerr.Add(graphlist_h[npar])
	legend.AddEntry(graphlist_p[npar],legs[0],"pl")
	#legend.AddEntry(graphlist_d[npar],legs[1],"pl")
	legend.AddEntry(graphlist_h[npar],legs[2],"pl")
	mgerr.Draw("AP")
	legend.Draw("same")
	mgerr.SetTitle("Trigger Eff Vs Rate, p_{T}>%d GeV, %.1f<|#eta|<%.2f, parity=%d"%(pt, etamin, etamax, npar))
	mgerr.SetMinimum(0)
	if npar==4:
	    mgerr.SetTitle("Trigger Eff Vs Rate, p_{T}>%d GeV, %.1f<|#eta|<%.2f, all parities"%(pt, etamin, etamax))

	mgerr.GetXaxis().SetTitle("Trigger Efficiency")
	mgerr.GetYaxis().SetTitle("Trigger rate [kHz]")
	c1.SaveAs(SBPlots+"TriggerEff_Rate_pt%d_eta%dto%d_npar%d%s.png"%(pt, int(etamin*10), int(etamax*10), npar, extrasuffix))
	c1.SaveAs(SBPlots+"TriggerEff_Rate_pt%d_eta%dto%d_npar%d%s.pdf"%(pt, int(etamin*10), int(etamax*10), npar, extrasuffix))
	c1.SaveAs(SBPlots+"TriggerEff_Rate_pt%d_eta%dto%d_npar%d%s.C"%(pt, int(etamin*10), int(etamax*10), npar, extrasuffix))
	mgerr.SetName("TriggerEffVsRate_all_pt%d_eta%dto%d_npar%d%s"%(pt, int(etamin*10), int(etamax*10), npar, extrasuffix))
	graphlist_p[npar].Write()
	graphlist_d[npar].Write()
	graphlist_h[npar].Write()
	mgerr.Write()
    wfile.Close()
    return totsignals, efflist_p[4], ratelist_p[4], efflist_d[4], ratelist_d[4], efflist_h[4], ratelist_h[4]

def CombinedGE21_SB(filename, pts, netas, SBPlots):
    tfile = TFile(filename, "UPDATE")
    color = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen+2, ROOT.kMagenta+2,ROOT.kCyan+2]
    maker = [20,21,22,23,33]
    legs = ["Position-based algorithm","Hybrid algorithm without GE21","Hybrid algorithm with GE21"]
    for pt in pts:
      for neta in range(len(netas)):
	if neta < len(netas)-1:
	    etamin = netas[neta]
	    etamax = netas[neta+1]
	    allnpar = 5
	else:
	    etamin = netas[0]
	    etamax = netas[-1]
	    allnpar = 1
	
	for npar in range(allnpar):
	    #print "pt ",pt, " etamin ",etamin, " etamax ",etamax," npar ",npar,"TriggerEffVsRate_hybrid_pt%d_eta%dto%d_npar%d"%(pt, int(etamin*10), int(etamax*10), npar)
	    legend = ROOT.TLegend(0.15,0.7,0.65,0.92)
	    legend.SetFillColor(ROOT.kWhite)
	    legend.SetTextFont(42)
	    legend.SetTextSize(.05)
	    c1 = TCanvas()
	    mgerr =TMultiGraph()
	    if allnpar>1:
		graph_GE21 = tfile.Get("TriggerEffVsRate_hybrid_pt%d_eta%dto%d_npar%d"%(int(pt), int(etamin*10), int(etamax*10), npar))
		graph_noGE21 = tfile.Get("TriggerEffVsRate_hybrid_pt%d_eta%dto%d_npar%d_noGE21"%(int(pt), int(etamin*10), int(etamax*10), npar))
		graph_p = tfile.Get("TriggerEffVsRate_position_pt%d_eta%dto%d_npar%d"%(int(pt), int(etamin*10), int(etamax*10), npar))
	    else:
		graph_GE21 = tfile.Get("TriggerEffVsRate_hybrid_pt%d_eta%dto%d"%(int(pt), int(etamin*10), int(etamax*10)))
		graph_noGE21 = tfile.Get("TriggerEffVsRate_hybrid_pt%d_eta%dto%d_noGE21"%(int(pt), int(etamin*10), int(etamax*10)))
		graph_p = tfile.Get("TriggerEffVsRate_position_pt%d_eta%dto%d"%(int(pt), int(etamin*10), int(etamax*10)))

	    graph_p.SetMarkerStyle(maker[0])
	    graph_p.SetMarkerColor(color[0])
	    graph_p.SetLineColor(color[0])
	    graph_noGE21.SetMarkerStyle(maker[3])
	    graph_noGE21.SetMarkerColor(color[3])
	    graph_noGE21.SetLineColor(color[3])
	    graph_GE21.SetMarkerStyle(maker[2])
	    graph_GE21.SetMarkerColor(color[2])
	    graph_GE21.SetLineColor(color[2])

	    mgerr.Add(graph_p)
	    mgerr.Add(graph_noGE21)
	    mgerr.Add(graph_GE21)

	    legend.AddEntry(graph_p,legs[0],"p")
	    legend.AddEntry(graph_noGE21,legs[1],"p")
	    legend.AddEntry(graph_GE21,legs[2],"p")
	    mgerr.Draw("AP")
	    legend.Draw("same")
	    #mgerr.SetTitle("Trigger Eff Vs Rate, p_{T}>%d GeV, %.2f<|#eta|<%.2f, parity=%d"%(pt, etamin, etamax, npar))
	    mgerr.SetTitle(" #scale[1.4]{#font[61]{CMS}} #font[52]{Simulation preliminary}"+"  "*10+" 14 TeV, 200 PU")
	    mgerr.SetMinimum(0)
	    tex = TLatex()
	    if npar==4 or allnpar==1:
		tex.Draw(0.2, 0.88,"p_{T}>%d GeV, %.2f<|#eta|<%.2f, all parities"%(pt, etamin, etamax))
	    else:
		tex.Draw(0.2, 0.88,"p_{T}>%d GeV, %.2f<|#eta|<%.2f, parity=%d"%(pt, etamin, etamax, npar))
	    mgerr.GetXaxis().SetTitle("Trigger Efficiency")
	    mgerr.GetYaxis().SetTitle("Trigger rate [kHz]")
	    if allnpar==5:
		c1.SaveAs(SBPlots+"TriggerEff_Rate_pt%d_eta%dto%d_npar%d_all.png"%(pt, int(etamin*10), int(etamax*10), npar))
		c1.SaveAs(SBPlots+"TriggerEff_Rate_pt%d_eta%dto%d_npar%d_all.pdf"%(pt, int(etamin*10), int(etamax*10), npar))
		c1.SaveAs(SBPlots+"TriggerEff_Rate_pt%d_eta%dto%d_npar%d_all.C"%(pt, int(etamin*10), int(etamax*10), npar))
		mgerr.SetName("TriggerEffVsRate_all_pt%d_eta%dto%d_npar%d_combined"%(pt, int(etamin*10), int(etamax*10), npar))
	    else:
		c1.SaveAs(SBPlots+"TriggerEff_Rate_pt%d_eta%dto%d_all.png"%(pt, int(etamin*10), int(etamax*10)))
		c1.SaveAs(SBPlots+"TriggerEff_Rate_pt%d_eta%dto%d_all.pdf"%(pt, int(etamin*10), int(etamax*10)))
		c1.SaveAs(SBPlots+"TriggerEff_Rate_pt%d_eta%dto%d_all.C"%(pt, int(etamin*10), int(etamax*10)))
		mgerr.SetName("TriggerEffVsRate_all_pt%d_eta%dto%d_combined"%(pt, int(etamin*10), int(etamax*10)))
	    mgerr.Write()
    tfile.Close()
	

def CombinedGE21_effAndrate(filename, pts, netas, SBPlots):
    tfile = TFile(filename, "READ")
    color = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen+2, ROOT.kMagenta+2,ROOT.kCyan+2]
    maker = [20,21,22,23,33]
    legs = ["Position-based algorithm","Hybrid algorithm without GE21","Hybrid algorithm with GE21"]
    pu = 200
    for pt in pts:
	efflist_p = []
	ratelist_p = []
	effErrlist_p = []
	rateErrlist_p = []
	efflist_h_noGE21 = []
	ratelist_h_noGE21 = []
	effErrlist_h_noGE21 = []
	rateErrlist_h_noGE21 = []
	efflist_h = []
	ratelist_h = []
	effErrlist_h = []
	rateErrlist_h = []
	for frac in fractions:
	    teff_h_prompt = []
	    teff_h_noGE21_prompt = []
	    teff_p_prompt = []
	    teff_h_displaced = []
	    teff_h_noGE21_displaced = []
	    teff_p_displaced = []
	    histrate_h = []
	    histrate_p = []
	    histrate_h_noGE21 = []
	    for neta in range(len(netas)-1):
		etamin = netas[neta]
		etamax = netas[neta+1]
		suffix = ""
		if etamin>=1.6:
		    suffix = "_noGE21"
		print "pt ",pt, " etamin ",etamin, " etamax ",etamax," frac ",frac
		teff_p_prompt.append(tfile.Get("Teff_displaced_position_eta%dto%d_npar4_frac%d_pt%d"%(int(etamin*10), int(etamax*10), frac, pt)))
		teff_h_prompt.append(tfile.Get("Teff_displaced_hybrid_eta%dto%d_npar4_frac%d_pt%d"%(int(etamin*10), int(etamax*10), frac, pt)))
		teff_h_noGE21_prompt.append(tfile.Get("Teff_displaced_hybrid_eta%dto%d_npar4_frac%d_pt%d"%(int(etamin*10), int(etamax*10), frac, pt)+suffix))
		getOverallEff(teff_p_prompt[neta], pt)
		getOverallEff(teff_h_prompt[neta], pt)
		getOverallEff(teff_h_noGE21_prompt[neta], pt)
		histrate_p.append(tfile.Get("ratehist_position_eta%dto%d_npar%d_frac%d_pt%d"%(int(etamin*10), int(etamax*10), 4, frac, pt)))
		histrate_h.append(tfile.Get("ratehist_hybrid_eta%dto%d_npar%d_frac%d_pt%d"%(int(etamin*10), int(etamax*10), 4, frac, pt)))
		histrate_h_noGE21.append(tfile.Get("ratehist_hybrid_eta%dto%d_npar%d_frac%d_pt%d"%(int(etamin*10), int(etamax*10), 4, frac, pt)+suffix))
	    for neta in range(1, len(netas)-1):
		teff_p_prompt[0].Add(teff_p_prompt[neta])
		teff_h_prompt[0].Add(teff_h_prompt[neta])
		teff_h_noGE21_prompt[0].Add(teff_h_noGE21_prompt[neta])
		if netas[neta] != 1.6:		
		    histrate_p[0].Add(histrate_p[neta])
		    histrate_h[0].Add(histrate_h[neta])
		    histrate_h_noGE21[0].Add(histrate_h_noGE21[neta])
		else:
		    addTwoHistogram(histrate_p[neta], histrate_p[0]); histrate_p[0]= histrate_p[neta].Clone()
		    addTwoHistogram(histrate_h[neta], histrate_h[0]); histrate_h[0] = histrate_h[neta].Clone()
		    addTwoHistogram(histrate_h_noGE21[neta], histrate_h_noGE21[0]); histrate_h_noGE21[0] = histrate_h_noGE21[neta].Clone()

	    etamin = netas[0]
	    etamax = netas[-1]
	    peff,peffErr = getOverallEff(teff_p_prompt[0], pt)
	    prate = histrate_p[0].GetEntries()*SF200
	    efflist_p.append(peff)
	    ratelist_p.append(prate)
	    effErrlist_p.append(peffErr)
	    rateErrlist_p.append(sqrt(prate*SF200))
	    hnoGE21eff,hnoGE21effErr = getOverallEff(teff_h_noGE21_prompt[0], pt)
	    hnoGE21rate = histrate_h_noGE21[0].GetEntries()*SF200
	    efflist_h_noGE21.append(hnoGE21eff)
	    ratelist_h_noGE21.append(hnoGE21rate)
	    effErrlist_h_noGE21.append(hnoGE21effErr)
	    rateErrlist_h_noGE21.append(sqrt(hnoGE21rate*SF200))
	    heff,heffErr = getOverallEff(teff_h_prompt[0], pt)
	    hrate = histrate_h[0].GetEntries()*SF200
	    efflist_h.append(heff)
	    ratelist_h.append(hrate)
	    effErrlist_h.append(heffErr)
	    rateErrlist_h.append(sqrt(hrate*SF200))
	    print "pt ",pt," peff ",peff," prate ",prate, " hnoGE21eff ",hnoGE21eff," hnoGE21rate ",hnoGE21rate," heff ",heff," hrate ",hrate," frac ",frac 
	     
	    teff_p_prompt[0].SetName("Teff_prompt_position_eta%dto%d_npar4_frac%d_pt%d_PU%d"%(int(etamin*10), int(etamax*10), frac, pt, pu))
	    teff_h_prompt[0].SetName("Teff_prompt_hybrid_eta%dto%d_npar4_frac%d_pt%d_PU%d"%(int(etamin*10), int(etamax*10), frac, pt, pu))
	    teff_h_noGE21_prompt[0].SetName("Teff_prompt_hybrid_eta%dto%d_npar4_frac%d_pt%d_noGE21_PU%d"%(int(etamin*10), int(etamax*10), frac, pt, pu))
	    histrate_p[0].SetName("ratehist_position_eta%dto%d_npar%d_frac%d_pt%d_PU%d"%(int(etamin*10), int(etamax*10), 4, frac, pt, pu))
	    histrate_h[0].SetName("ratehist_hybrid_eta%dto%d_npar%d_frac%d_pt%d_PU%d"%(int(etamin*10), int(etamax*10), 4, frac, pt, pu))
	    histrate_h_noGE21[0].SetName("ratehist_hybrid_eta%dto%d_npar%d_frac%d_pt%d_noGE21_PU%d"%(int(etamin*10), int(etamax*10), 4, frac, pt, pu))
	    teff_p_prompt[0].SetMarkerStyle(maker[0])
	    teff_p_prompt[0].SetMarkerColor(color[0])
	    teff_p_prompt[0].SetLineColor(color[0])
	    histrate_p[0].SetMarkerStyle(maker[0])
	    histrate_p[0].SetMarkerColor(color[0])
	    histrate_p[0].SetLineColor(color[0])
	    teff_h_noGE21_prompt[0].SetMarkerStyle(maker[3])
	    teff_h_noGE21_prompt[0].SetMarkerColor(color[3])
	    teff_h_noGE21_prompt[0].SetLineColor(color[3])
	    histrate_h_noGE21[0].SetMarkerStyle(maker[3])
	    histrate_h_noGE21[0].SetMarkerColor(color[3])
	    histrate_h_noGE21[0].SetLineColor(color[3])
	    teff_h_prompt[0].SetMarkerStyle(maker[2])
	    teff_h_prompt[0].SetMarkerColor(color[2])
	    teff_h_prompt[0].SetLineColor(color[2])
	    histrate_h[0].SetMarkerStyle(maker[2])
	    histrate_h[0].SetMarkerColor(color[2])
	    histrate_h[0].SetLineColor(color[2])

	    tex1 = TLatex(.4, .6,"p_{T}^{L1}>%d GeV, %.2f<|#eta|<%.2f"%(pt, etamin, etamax))
	    tex1.SetNDC()
	    c1 = TCanvas()
	    c1.SetGridx()
	    c1.SetGridy()
	    b1 = ROOT.TH1F("b1","b1",len(ptbins)-1, ptbins)
	    b1.GetYaxis().SetRangeUser(0.0,1.05)
	    b1.GetYaxis().SetTitleOffset(1.1)
	    b1.GetYaxis().SetNdivisions(520)
	    b1.GetYaxis().SetTitle("Trigger efficiency")
	    b1.GetXaxis().SetTitle("True muon p_{T} [GeV]")
	    b1.GetXaxis().SetTitleSize(0.05)
	    b1.GetXaxis().SetLabelSize(0.05)
	    b1.GetYaxis().SetTitleSize(0.05)
	    b1.GetYaxis().SetLabelSize(0.05)
	    b1.SetTitle(" #scale[1.4]{#font[61]{CMS}} #font[52]{Simulation preliminary}"+"  "*10+" 14 TeV, 200 PU")
	
	    legend1 = ROOT.TLegend(0.4,0.2,0.88,0.5)
	    legend1.SetFillColor(ROOT.kWhite)
	    legend1.SetTextFont(42)
	    #legend1.SetTextSize(.05)
	    legend1.AddEntry(teff_p_prompt[0], legs[0], "p")
	    legend1.AddEntry(teff_h_noGE21_prompt[0], legs[1], "p")
	    legend1.AddEntry(teff_h_prompt[0], legs[2], "p")
	    b1.Draw()
	    teff_p_prompt[0].Draw("samep")
	    teff_h_noGE21_prompt[0].Draw("samep")
	    teff_h_prompt[0].Draw("samep")
	    legend1.Draw("same")
	    tex1.Draw("same")
	    #c1.SaveAs(SBPlots+"TriggerEff_turnon_pt%d_eta%dto%d_frac%d_PU%d.png"%(pt, int(etamin*10), int(etamax*10), frac, pu))
	    c1.SaveAs(SBPlots+"TriggerEff_turnon_pt%d_eta%dto%d_frac%d_PU%d.pdf"%(pt, int(etamin*10), int(etamax*10), frac, pu))
	    c1.SaveAs(SBPlots+"TriggerEff_turnon_pt%d_eta%dto%d_frac%d_PU%d.C"%(pt, int(etamin*10), int(etamax*10), frac, pu))

	    c2 = TCanvas()
	    c2.Clear()
	    b2 = ROOT.TH1F("b2","b2",len(myetabin)-1, myetabin)
	    b2.GetYaxis().SetRangeUser(0.0, 50)
	    b2.GetYaxis().SetTitleOffset(1.1)
	    b2.GetYaxis().SetNdivisions(510)
	    b2.GetYaxis().SetTitle("Trigger Rate [kHz]")
	    b2.GetXaxis().SetTitle("L1Mu |#eta|")
	    b2.GetXaxis().SetTitleSize(0.05)
	    b2.GetXaxis().SetLabelSize(0.05)
	    b2.GetYaxis().SetTitleSize(0.05)
	    b2.GetYaxis().SetLabelSize(0.05)
	    b2.SetTitle(" #scale[1.4]{#font[61]{CMS}} #font[52]{Simulation preliminary}"+"  "*10+" 14 TeV, %d PU"%pu)
	    c2.SetGridx()
	    c2.SetGridy()

	    tex2 = TLatex(.15, .85,"p_{T}^{L1}>%d GeV, %.2f<|#eta|<%.2f, overall efficiency %d%%"%(pt, etamin, etamax, frac))
	    tex2.SetNDC()
	    legend2 = ROOT.TLegend(0.15,0.6,0.6,0.8)
	    legend2.SetFillColor(ROOT.kWhite)
	    legend2.SetTextFont(42)
	    #legend2.SetTextSize(.04)
	    legend2.AddEntry(histrate_p[0], legs[0], "p")
	    legend2.AddEntry(histrate_h_noGE21[0], legs[1], "p")
	    legend2.AddEntry(histrate_h[0], legs[2], "p")
	    b2.Draw()
	    histrate_p[0].Scale(SF200)
	    histrate_h_noGE21[0].Scale(SF200)
	    histrate_h[0].Scale(SF200)

	    histrate_p[0].Draw("samee")
	    histrate_h_noGE21[0].Draw("samee")
	    histrate_h[0].Draw("samee")
	    legend2.Draw("same")
	    tex2.Draw("same")
	    #c2.SaveAs(SBPlots+"TriggerRate_eta_pt%d_eta%dto%d_frac%d_PU%d.png"%(pt, int(etamin*10), int(etamax*10), frac, pu))
	    c2.SaveAs(SBPlots+"TriggerRate_eta_pt%d_eta%dto%d_frac%d_PU%d.pdf"%(pt, int(etamin*10), int(etamax*10), frac, pu))
	    c2.SaveAs(SBPlots+"TriggerRate_eta_pt%d_eta%dto%d_frac%d_PU%d.C"%(pt, int(etamin*10), int(etamax*10), frac, pu))
	etamin = netas[0]
	etamax = netas[-1]
	tex3 = TLatex(.2, .88,"p_{T}^{L1}>%d GeV, %.2f<|#eta|<%.2f"%(pt, etamin, etamax))
	tex3.SetNDC()
	legend = ROOT.TLegend(0.15,0.7,0.65,0.85)
	legend.SetFillColor(ROOT.kWhite)
	legend.SetTextFont(42)
	#legend.SetTextSize(.04)
	c3 = TCanvas()
	mgerr =TMultiGraph()
	graph_p = TGraphErrors(len(efflist_p), np.array(efflist_p), np.array(ratelist_p), np.array(effErrlist_p), np.array(rateErrlist_p))
	graph_noGE21 = TGraphErrors(len(efflist_h_noGE21), np.array(efflist_h_noGE21), np.array(ratelist_h_noGE21), np.array(effErrlist_h_noGE21), np.array(rateErrlist_h_noGE21))
	graph_GE21 = TGraphErrors(len(efflist_h), np.array(efflist_h), np.array(ratelist_h), np.array(effErrlist_h), np.array(rateErrlist_h))
	graph_p.SetMarkerStyle(maker[0])
	graph_p.SetMarkerColor(color[0])
	graph_p.SetLineColor(color[0])
	graph_noGE21.SetMarkerStyle(maker[3])
	graph_noGE21.SetMarkerColor(color[3])
	graph_noGE21.SetLineColor(color[3])
	graph_GE21.SetMarkerStyle(maker[2])
	graph_GE21.SetMarkerColor(color[2])
	graph_GE21.SetLineColor(color[2])

	mgerr.Add(graph_p)
	mgerr.Add(graph_noGE21)
	mgerr.Add(graph_GE21)

	legend.AddEntry(graph_p,legs[0],"p")
	legend.AddEntry(graph_noGE21,legs[1],"p")
	legend.AddEntry(graph_GE21,legs[2],"p")
	mgerr.Draw("AP")
	legend.Draw("same")
	tex3.Draw("same")
	mgerr.SetMinimum(0)
	mgerr.SetTitle(" #scale[1.4]{#font[61]{CMS}} #font[52]{Simulation preliminary}"+"  "*10+" 14 TeV, %d PU"%pu)
	mgerr.GetXaxis().SetTitle("Trigger Efficiency")
	mgerr.GetYaxis().SetTitle("Trigger Rate [kHz]")
	#c3.SaveAs(SBPlots+"TriggerEff_Rate_pt%d_eta%dto%d_all_v2.png"%(pt, int(etamin*10), int(etamax*10)))
	c3.SaveAs(SBPlots+"TriggerEff_Rate_pt%d_eta%dto%d_all_PU%d.pdf"%(pt, int(etamin*10), int(etamax*10), pu))
	c3.SaveAs(SBPlots+"TriggerEff_Rate_pt%d_eta%dto%d_all_PU%d.C"%(pt, int(etamin*10), int(etamax*10), pu))
	mgerr.SetName("TriggerEffVsRate_all_pt%d_eta%dto%d_combined_PU%d"%(pt, int(etamin*10), int(etamax*10), pu))
	



def getallPtEffAndRate(pts, fractions, etas, rootfile, filedir, condoroutroot, SBPlots, useGE21, doPlots, doLUTs, doRates):
    for pt in pts:
	peffalleta = []
	pratealleta = []
	deffalleta = []
	dratealleta = []
	heffalleta = []
	hratealleta = []
	nalleta = []
	for neta in range(len(etas)-1):
	    etamin = etas[neta]
	    etamax = etas[neta+1]
	    #all 4 parities efflist_p[nfrac]
	    ntot, efflist_p, ratelist_p, efflist_d, ratelist_d,efflist_h, ratelist_h = getEffAndRate(pt, fractions, etamin, etamax, rootfile, filedir, condoroutroot, SBPlots, useGE21, doPlots, doLUTs, doRates)
	    nalleta.append(ntot)
	    peffalleta.append(efflist_p), pratealleta.append(ratelist_p)
	    deffalleta.append(efflist_d), dratealleta.append(ratelist_d)
	    heffalleta.append(efflist_h), hratealleta.append(ratelist_h)

	if len(etas)<3:	continue
	def getfinalEffAndRate(nalleta, effalleta,  ratealleta, nfrac):
	    totden = 0.0
	    totnum = 0.0
	    totrate = 0.0
	    for neta in range(len(etas)-1):
		totden += nalleta[neta]
		totnum += effalleta[neta][nfrac]*nalleta[neta]
		totrate += ratealleta[neta][nfrac]
	    return totnum, totden, totrate
		

	efflist_p_alleta = []
	ratelist_p_alleta = []
	effErrlist_p_alleta = []
	rateErrlist_p_alleta = []
	efflist_d_alleta = []
	ratelist_d_alleta = []
	effErrlist_d_alleta = []
	rateErrlist_d_alleta = []
	efflist_h_alleta = []
	ratelist_h_alleta = []
	effErrlist_h_alleta = []
	rateErrlist_h_alleta = []
	nfrac = -1
	etamin = netas[0]; etamax = netas[-1]
	for frac in fractions:
	    nfrac += 1
	    ptotnum, ptotden, ptotrate  = getfinalEffAndRate(nalleta, peffalleta,  pratealleta, nfrac)
	    dtotnum, dtotden, dtotrate  = getfinalEffAndRate(nalleta, deffalleta,  dratealleta, nfrac)
	    htotnum, htotden, htotrate  = getfinalEffAndRate(nalleta, heffalleta,  hratealleta, nfrac)
	    peff = float(ptotnum)/float(ptotden)
	    efflist_p_alleta.append(peff)
	    ratelist_p_alleta.append(ptotrate)
	    effErrlist_p_alleta.append(sqrt((1-peff)*peff/ptotden))
	    rateErrlist_p_alleta.append(sqrt(ptotrate*SF))
	    deff = float(dtotnum)/float(dtotden)
	    efflist_d_alleta.append(deff)
	    ratelist_d_alleta.append(dtotrate)
	    effErrlist_d_alleta.append(sqrt((1-deff)*deff/dtotden))
	    rateErrlist_d_alleta.append(sqrt(dtotrate*SF))
	    heff = float(htotnum)/float(htotden)
	    efflist_h_alleta.append(heff)
	    ratelist_h_alleta.append(htotrate)
	    effErrlist_h_alleta.append(sqrt((1-heff)*heff/htotden))
	    rateErrlist_h_alleta.append(sqrt(htotrate*SF))
	    print "pt ",pt, " etamin ",etamin," etamax ",etamax, " frac ",frac," ptotden ",ptotden," ptotnum ",ptotnum," htotnum ",htotnum," nalleta ",nalleta," peff ",peff, " heff ",heff
	extrasuffix = ""
	extraleg = "with GE21"
	if not(useGE21):
	    extrasuffix = "_noGE21"
	    extraleg = "w/o GE21"
	if etamin<1.6:
	    extraleg = ""
	wfile = TFile(rootfile, "UPDATE") 
	#print "total rate ",ratelist_h 
	legs = ["Position-based algorithm","Stub alignment algorithm %s"%extraleg,"Hybrid algorithm %s"%extraleg]
	tex = TLatex(0.2, 0.88,"p_{T}>%d GeV, %.2f<|#eta|<%.2f, all parities"%(pt, etamin, etamax))
	tex.SetNDC()
	legend = ROOT.TLegend(0.15,0.7,0.63,0.85)
	legend.SetFillColor(ROOT.kWhite)
	legend.SetTextFont(42)
	#legend.SetTextSize(.05)
	c1 = TCanvas()
	mgerr =TMultiGraph()
	graph_p = TGraphErrors(len(efflist_p_alleta), np.array(efflist_p_alleta), np.array(ratelist_p_alleta), np.array(effErrlist_p_alleta), np.array(rateErrlist_p_alleta))
	graph_d = TGraphErrors(len(efflist_d_alleta), np.array(efflist_d_alleta), np.array(ratelist_d_alleta), np.array(effErrlist_d_alleta), np.array(rateErrlist_d_alleta))
	graph_h = TGraphErrors(len(efflist_h_alleta), np.array(efflist_h_alleta), np.array(ratelist_h_alleta), np.array(effErrlist_h_alleta), np.array(rateErrlist_h_alleta))
	graph_p.SetMarkerStyle(maker[0])
	graph_p.SetMarkerColor(color[0])
	graph_d.SetMarkerStyle(maker[1])
	graph_d.SetMarkerColor(color[1])
	graph_h.SetMarkerStyle(maker[2])
	graph_h.SetMarkerColor(color[2])
	etamin = netas[0]; etamax = netas[-1]
	graph_p.SetName("TriggerEffVsRate_position_pt%d_eta%dto%d%s"%(pt, int(etamin*10), int(etamax*10), extrasuffix))
	graph_d.SetName("TriggerEffVsRate_direction_pt%d_eta%dto%d%s"%(pt, int(etamin*10), int(etamax*10), extrasuffix))
	graph_h.SetName("TriggerEffVsRate_hybrid_pt%d_eta%dto%d%s"%(pt, int(etamin*10), int(etamax*10), extrasuffix))
	mgerr.SetName("TriggerEffVsRate_all_pt%d_eta%dto%d%s"%(pt, int(etamin*10), int(etamax*10), extrasuffix))
	
	mgerr.Add(graph_p)
	#mgerr.Add(graph_d)
	mgerr.Add(graph_h)
	legend.AddEntry(graph_p,legs[0],"pl")
	#legend.AddEntry(graph_d,legs[1],"pl")
	legend.AddEntry(graph_h,legs[2],"pl")
	mgerr.Draw("AP")
	legend.Draw("same")
	tex.Draw("same")
	mgerr.SetMinimum(0)
	mgerr.SetTitle(" #scale[1.4]{#font[61]{CMS}} #font[52]{Simulation preliminary}"+"  "*10+" 14 TeV, 200 PU")
	#mgerr.SetTitle("Trigger Eff Vs Rate, p_{T}>%d GeV, %.2f<|#eta|<%.2f, all parities"%(pt, etamin, etamax))
	mgerr.GetXaxis().SetTitle("Trigger Efficiency")
	mgerr.GetYaxis().SetTitle("Trigger rate [kHz]")
	c1.SaveAs(SBPlots+"TriggerEff_Rate_pt%d_eta%dto%d%s.png"%(pt, int(etamin*10), int(etamax*10),  extrasuffix))
	c1.SaveAs(SBPlots+"TriggerEff_Rate_pt%d_eta%dto%d%s.pdf"%(pt, int(etamin*10), int(etamax*10),  extrasuffix))
	c1.SaveAs(SBPlots+"TriggerEff_Rate_pt%d_eta%dto%d%s.C"%(pt, int(etamin*10), int(etamax*10),  extrasuffix))
	graph_p.Write()
	graph_d.Write()
	graph_h.Write()
	mgerr.Write()

	wfile.Close()

	     

def plotalleta(pt, netas, fraction, outputrootfile, ratesample, signalsample, allLUTfile, SBPlots):
    #fraction_d = 99
    rate_y = "dphi_dir12_GE21"
    rate_x = "ddY123"
    xaxis = "(deltay123_test)"
    yaxis = "(dphi_dir_st1_st2_L1)"
    yaxis_noGE21 = "(dphi_dir_st1_st2_L1_csc)"
    rate_y_noGE21 = "dphi_dir12_noGE21"
    legs3_v2 = ["Position-based algorithm","Stub alignment algorithm","Hybrid algorithm"]
    legs3 = ["Position-based algorithm","Hybrid algorithm w/o GE21","Hybrid algorithm w GE21"]

    ratefile = TFile(ratesample)
    h_total = ratefile.Get("h_eventcountME0Segment192")

    rateEvents = h_total.GetEntries()
    L1MuTriggerRate = TChain("L1MuTriggerRateME0Segment192")
    addFilesToChain(L1MuTriggerRate, ratesample)
    GEMCSCTrackCh0GE11 = TChain("GEMCSCAnalyzer/trk_eff_CSC_ALL")
    #GEMCSCTrackCh1GE11 = TChain("GEMCSCAnalyzer/trk_eff_CSC_ME11")
    GEMCSCTrackCh2GE11 = TChain("GEMCSCAnalyzer/trk_eff_CSC_ME32")
    GEMCSCTrackCh3GE11 = TChain("GEMCSCAnalyzer/trk_eff_CSC_ME31")
    GEMCSCTrackCh0 = TChain("GEMCSCAnalyzer/trk_eff_CSC_ALL")
    #GEMCSCTrackCh1 = TChain("GEMCSCAnalyzer/trk_eff_CSC_ME11")
    GEMCSCTrackCh2 = TChain("GEMCSCAnalyzer/trk_eff_CSC_ME32")
    GEMCSCTrackCh3 = TChain("GEMCSCAnalyzer/trk_eff_CSC_ME31")
    #GEMCSCTrackCh1 = TChain("GEMCSCAnalyzer/trk_eff_CSC_ALL")
    #addFilesToChain(GEMCSCTrackCh1, signalsample) 
    addFilesToChain(GEMCSCTrackCh0GE11, signalsampleGE11) 
    #addFilesToChain(GEMCSCTrackCh1GE11, signalsampleGE11) 
    addFilesToChain(GEMCSCTrackCh2GE11, signalsampleGE11) 
    addFilesToChain(GEMCSCTrackCh3GE11, signalsampleGE11) 
    addFilesToChain(GEMCSCTrackCh0, signalsample) 
    #addFilesToChain(GEMCSCTrackCh1, signalsample) 
    addFilesToChain(GEMCSCTrackCh2, signalsample) 
    addFilesToChain(GEMCSCTrackCh3, signalsample) 
    #GEMCSCTrackCh0GE11.AddFriend(GEMCSCTrackCh1GE11)
    GEMCSCTrackCh0GE11.AddFriend(GEMCSCTrackCh2GE11)
    GEMCSCTrackCh0GE11.AddFriend(GEMCSCTrackCh3GE11)
    #GEMCSCTrackCh0.AddFriend(GEMCSCTrackCh1)
    GEMCSCTrackCh0.AddFriend(GEMCSCTrackCh2)
    GEMCSCTrackCh0.AddFriend(GEMCSCTrackCh3)

    Teff_out = ROOT.TFile(outputrootfile,"UPDATE")
    Teffall = []
    Teffall2 = []
    Teffalleta = []
    Teffalleta2 = []
    Teffalldxy = []
    Teffalldxy2 = []
    hrateall = []
    hrateall2 = []
    quality = 12
    for neta in range(len(netas)-1):
       if (netas[neta]<0 or netas[neta+1]<0):
	    continue
       ring = 1
       if (netas[neta]<1.6):
	    ring = 2
    #if (netas[0]==1.6):
       etamin = netas[neta]
       etamax = netas[neta+1]
       Teffallnparseta = []
       Teffallnparseta_veto2 = []
       Teffallnparsdxy = []
       Teffallnparsdxy_veto2 = []
       Teffallnpars = []
       Teffallnpars_veto2 = []
       hrateallnpars = []
       hrateallnpars_veto2 = []
       #Teffallnpars_p = []
       #Teffallnpars_d = []
       legs = ["#splitline{Displaced muon}{10<|d_{xy}|<50 cm}"]
       for npar in range(0, 4):
	    me11 = evenodds[npar].split(',')[0]
	    me21 = evenodds[npar].split(',')[1]
	    chambers = "ME1%d %s,ME2%d %s,ME3%d %s"%(ring, me11, ring, me21, ring, me21)
	    me31_n = 1
	    if me21 == "even":
		    me31_n = 2
	    ME31cut = " && ((trk_eff_CSC_ME31.has_lct&%d)>0 || (trk_eff_CSC_ME32.has_lct&%d)>0)"%(me31_n, me31_n)
	    #st_title = ["Prompt muon, 2<p_{T}<%d"%pt1, "Displaced Muon, 10<|d_{xy}|<50, p_{T}>%d"%pt]
	    algos_lut = getFinalLUTValue(pt, fraction, etamin, etamax, npar, "", allLUTfile)
	    algos_lut_noGE21 = getFinalLUTValue(pt, fraction, etamin, etamax, npar, "_noGE21", allLUTfile)
	    ddY_cut = algos_lut["Position"]
	    dphi_dir_cut = algos_lut["Direction"]
	    maxa = algos_lut["Hybrid"][0] 
	    maxb = algos_lut["Hybrid"][1] 
	    alpha= algos_lut["Hybrid"][2]
	    x0 =  algos_lut["Hybrid"][3]
	    y0 =  algos_lut["Hybrid"][4]
            print "Position ddY_cut ",ddY_cut," Direction ",dphi_dir_cut," hybrid ",algos_lut["Hybrid"]
	    maxa1 = maxa; maxb1 = maxb; alpha1 = alpha; x01 = x0; y01= y0
	    if ring==1:
		maxa1 = algos_lut_noGE21["Hybrid_noGE21"][0] 
		maxb1 = algos_lut_noGE21["Hybrid_noGE21"][1] 
		alpha1= algos_lut_noGE21["Hybrid_noGE21"][2]
		x01 =  algos_lut_noGE21["Hybrid_noGE21"][3]
		y01 =  algos_lut_noGE21["Hybrid_noGE21"][4]
                print " Hybrid no GE21 ",algos_lut_noGE21
	    
            #L1Mucut = "&& L1Mu_quality>=4 && hasME2"
            L1Mucut = "&& 1"
	    ratecut = "maxPositionPt&& abs(L1Mu_eta_st2)>=%f && abs(L1Mu_eta_st2)<%f && npar==%d"%(etamin, etamax, npar)+L1Mucut
	    #etacuts ="&& meRing==%d"%ring+"&& fabs(trk_eff_CSC_ME2%d.eta_lct_%s)>%f && fabs(trk_eff_CSC_ME2%d.eta_lct_%s)<%f"%(ring,me21,etamin,ring,me21,etamax)
	    checkvalue = "&& fabs(deltay23_test)<50 && fabs(deltay12_test)<50 && fabs(%s)>0 && fabs(%s)>0  && abs(dphi_dir_st1_st2_L1)<2.0 && meRing==%d"%(xaxis, yaxis,ring)
	    checkx = "&& fabs(%s)<60 && fabs(%s)>0 "%(xaxis, xaxis)
	    checky = "&& fabs(%s)<2.0 && fabs(%s)>0 "%(yaxis, yaxis)
	    checky_noGE21 = "&& fabs(%s)<2.0 && fabs(%s)>0 "%(yaxis_noGE21, yaxis_noGE21)
	    if ring==2:
		checky_noGE21 = checky
	    checkx = ""; checky = ""; checky_noGE21 = ""
            dxycut = "&& fabs(genGdMu_dxy)>10 && fabs(genGdMu_dxy)<50 && fabs(genGdMu_dR)<0.1"
            #dxycut = " && fabs(dxy)>10 && fabs(dxy)<50"
            
	    cutsbase = "has_gmtCand>0 && meRing==%d && "%ring+"fabs(eta_st2_sim)>%f && fabs(eta_st2_sim)<%f && npar_sim==%d"%(etamin,etamax, npar)+ME31cut+dxycut+checkvalue+" && abs(dxy)>1.0"
            cutsbase_nodxycut = "has_gmtCand>0 && meRing==%d && "%ring+"fabs(eta_st2_sim)>%f && fabs(eta_st2_sim)<%f && npar_sim==%d"%(etamin,etamax, npar)+ME31cut+checkvalue + "&& abs(dxy)>1.0"
	    cuts_p =  [cutsbase+checkx]
	    cuts_d =  [cutsbase+checky]

	    text = "#splitline{%s}{%.1f<|#eta|<%.2f, p_{T}^{L1}>%d GeV}"%(chambers, etamin,etamax,pt)
	    text_position = "#splitline{%s}{%.1f<|#eta|<%.2f, p_{T}^{L1}>%d GeV}"%(chambers,etamin,etamax,pt)
	    #Teffs_ddY = makeEffplot_pt([GEMCSCTrackCh0], "pt", cuts_p, ["fabs(%s)<%f"%(xaxis, ddY_cut)], etamin, etamax,"true muon p_{T} GeV","Trigger efficiency",legs,"Position-based algo" ,text_position,SBPlots+"GEMCSC_dxy10_50_position_eff_20170131_pt%d_fraction%d_st2eta%dto%d_npar%d"%(pt, fraction, int(etamin*10),int(etamax*10), npar))
	    #Teffs_dphi_dir = makeEffplot_pt([GEMCSCTrackCh0], "pt", cuts_d, ["fabs(%s)<%f"%(yaxis, dphi_dir_cut)], etamin, etamax,"true muon p_{T} GeV","Trigger efficiency",legs,"Stub alignment algo",text,SBPlots+"GEMCSC_dxy10_50_direction_eff_20170131_pt%d_fraction%d_st2eta%dto%d_npar%d"%(pt, fraction, int(etamin*10),int(etamax*10), npar))
	    #(maxa, maxb, alpha, x0, y0) = loopEllipse(L1MuTriggerRate, GEMCSCTrackCh1, nEvent_list, ddYDPhiCuts, fraction, npar, astart, bstart, xaxislist, yaxislist,x_bins, y_bins,xtitle, ytitle,st_title, etamin, etamax, cut_list, text,outputdir1+"GEMCSC_ctau0andctau1000_hyrid_profile_20170131_pt%d_ptbg%d_fraction%d_st2eta%dto%d_npar%d%s"%(pt, pt1, fraction, int(etamin*10),int(etamax*10), npar, extrapic))
	    xaxis1_rate = "(%s*TMath::Cos(%f)+%s*TMath::Sin(%f)-%f)"%(rate_x, alpha, rate_y, alpha, x0)
	    yaxis1_rate = "(%s*TMath::Sin(%f)-%s*TMath::Cos(%f)-%f)"%(rate_x, alpha, rate_y, alpha, y0)
	    ellipsecut_rate = "(%s*%s/(%f*%f)+%s*%s/(%f*%f))<=1.0"%(xaxis1_rate, xaxis1_rate, maxa, maxa, yaxis1_rate, yaxis1_rate, maxb, maxb)
	    xaxis1_rate_noGE21 = "(%s*TMath::Cos(%f)+%s*TMath::Sin(%f)-%f)"%(rate_x, alpha1, rate_y_noGE21, alpha1, x01)
	    yaxis1_rate_noGE21 = "(%s*TMath::Sin(%f)-%s*TMath::Cos(%f)-%f)"%(rate_x, alpha1, rate_y_noGE21, alpha1, y01)
	    ellipsecut_rate_noGE21 = "(%s*%s/(%f*%f)+%s*%s/(%f*%f))<=1.0"%(xaxis1_rate_noGE21, xaxis1_rate_noGE21, maxa1, maxa1, yaxis1_rate_noGE21, yaxis1_rate_noGE21, maxb1, maxb1)
	    if ring ==2:
		ellipsecut_rate_noGE21 = ellipsecut_rate
	    etacut_h = "&& abs(L1Mu_eta_st2)>1.65 && abs(L1Mu_eta_st2)<2.15 "
	    ratehist_p_veto0 = getRateVsEta(L1MuTriggerRate, "ratehist_position_eta%dto%d_npar%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), npar, fraction, pt), myetabin, ratecut + "&& abs(%s)<=%f"%(rate_x, ddY_cut), 0)
	    ratehist_h_veto0 = getRateVsEta(L1MuTriggerRate, "ratehist_hybrid_eta%dto%d_npar%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), npar, fraction, pt), myetabin, ratecut + "&& "+ellipsecut_rate+etacut_h, 0)
	    print "hrate cut ",ratecut + "&& "+ellipsecut_rate+etacut_h
	    ratehist_h_noGE21_veto0 = getRateVsEta(L1MuTriggerRate, "ratehist_hybrid_noGE21_eta%dto%d_npar%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), npar, fraction, pt), myetabin, ratecut + "&& "+ellipsecut_rate_noGE21, 0)
	    ratehist_p_veto2 = getRateVsEta(L1MuTriggerRate, "ratehist_position_mediumveto_eta%dto%d_npar%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), npar, fraction, pt), myetabin, ratecut + "&& abs(%s)<=%f"%(rate_x, ddY_cut), 2)
	    ratehist_h_veto2 = getRateVsEta(L1MuTriggerRate, "ratehist_hybrid_mediumveto_eta%dto%d_npar%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), npar, fraction, pt), myetabin, ratecut + "&& "+ellipsecut_rate+etacut_h, 2)
	    ratehist_h_noGE21_veto2 = getRateVsEta(L1MuTriggerRate, "ratehist_hybrid_noGE21_mediumveto_eta%dto%d_npar%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), npar, fraction, pt), myetabin, ratecut + "&& "+ellipsecut_rate_noGE21, 2)
	    hrateallnpars.append([ratehist_p_veto0, ratehist_h_noGE21_veto0, ratehist_h_veto0])
	    hrateallnpars_veto2.append([ratehist_p_veto2, ratehist_h_noGE21_veto2, ratehist_h_veto2])
	    ratehist_p_veto0.Write()
	    ratehist_h_noGE21_veto0.Write()
	    ratehist_h_veto0.Write()
	    ratehist_p_veto2.Write()
	    ratehist_h_noGE21_veto2.Write()
	    ratehist_h_veto2.Write()
	    #print "max_a ",maxa," maxb ",maxb," alpha ",alpha," x0 ",x0, " y0 ",y0
	    #xaxis1 = "(%s*TMath::Cos(%f)*charge+%s*TMath::Sin(%f)*charge-%f)"%(xaxis, alpha, yaxis, alpha, x0)
	    #yaxis1 = "(%s*TMath::Sin(%f)*charge-%s*TMath::Cos(%f)*charge-%f)"%(xaxis, alpha, yaxis, alpha, y0)
	    xaxis1 = "(%s*TMath::Cos(%f)+%s*TMath::Sin(%f)-%f)"%(xaxis, alpha, yaxis, alpha, x0)
	    yaxis1 = "(%s*TMath::Sin(%f)-%s*TMath::Cos(%f)-%f)"%(xaxis, alpha, yaxis, alpha, y0)
	    ellipse = "(%s*%s/(%f*%f)+%s*%s/(%f*%f))<=1.0"%(xaxis1, xaxis1, maxa, maxa, yaxis1, yaxis1, maxb, maxb)
	    xaxis1_noGE21 = "(%s*TMath::Cos(%f)+%s*TMath::Sin(%f)-%f)"%(xaxis, alpha1, yaxis_noGE21, alpha1, x01)
	    yaxis1_noGE21 = "(%s*TMath::Sin(%f)-%s*TMath::Cos(%f)-%f)"%(xaxis, alpha1, yaxis_noGE21, alpha1, y01)
	    ellipse_noGE21 = "(%s*%s/(%f*%f)+%s*%s/(%f*%f))<=1.0"%(xaxis1_noGE21, xaxis1_noGE21, maxa1, maxa1, yaxis1_noGE21, yaxis1_noGE21, maxb1, maxb1)
	    if ring ==2 :
		ellipse_noGE21 = ellipse
            print "position den cut ",cutsbase+checkx," num cut ","fabs(%s)<%f"%(xaxis, ddY_cut)+L1Mucut
            print "hybrid nogem den cut ",cutsbase+checkx+checky_noGE21," num cut ", ellipse_noGE21+L1Mucut
            print "hybrid with gem den cut ",cutsbase+checkx+checky," num cut ",ellipse+L1Mucut
	    Teffs = makeEffplot_pt([GEMCSCTrackCh0GE11, GEMCSCTrackCh0GE11, GEMCSCTrackCh0], "pt", [cutsbase+checkx, cutsbase+checkx+checky_noGE21, cutsbase+checkx+checky], ["fabs(%s)<%f"%(xaxis, ddY_cut)+L1Mucut, ellipse_noGE21+L1Mucut, ellipse+L1Mucut], etamin, etamax,"true muon p_{T} GeV","Trigger efficiency",legs3,"no track veto",text,SBPlots+"GEMCSC_dxy10_50_combined_effvspt_20170131_pt%d_fraction%d_st2eta%dto%d_npar%d"%(pt, fraction, int(etamin*10),int(etamax*10), npar))
	    Teffs_veto2 = makeEffplot_pt([GEMCSCTrackCh0GE11, GEMCSCTrackCh0GE11, GEMCSCTrackCh0], "pt", [cutsbase+checkx, cutsbase+checkx+checky_noGE21, cutsbase+checkx+checky], ["fabs(%s)<%f"%(xaxis, ddY_cut)+"&& !isL1MediumVeto"+L1Mucut, ellipse_noGE21+"&& !isL1MediumVeto"+L1Mucut, ellipse+"&& !isL1MediumVeto"+L1Mucut], etamin, etamax,"true muon p_{T} GeV","Trigger efficiency",legs3,"medium track veto",text,SBPlots+"GEMCSC_dxy10_50_mediumveto_combined_effvspt_20170131_pt%d_fraction%d_st2eta%dto%d_npar%d"%(pt, fraction, int(etamin*10),int(etamax*10), npar))
	    Teffseta = makeEffplot_eta([GEMCSCTrackCh0, GEMCSCTrackCh0GE11, GEMCSCTrackCh0], "abs(eta_st2_sim)", [cutsbase+checkx+"&& pt>=%d"%(pt+5), cutsbase+checkx+checky_noGE21+"&& pt>=%d"%(pt+5), cutsbase+checkx+checky+"&& pt>=%d"%(pt+5)], ["fabs(%s)<%f"%(xaxis, ddY_cut)+L1Mucut, ellipse_noGE21+L1Mucut, ellipse+L1Mucut], etamin, etamax,"|#eta| at 2nd CSC station","Trigger efficiency",legs3,"no track veto",text,SBPlots+"GEMCSC_dxy10_50_combined_effvseta_20170131_pt%d_fraction%d_st2eta%dto%d_npar%d"%(pt, fraction, int(etamin*10),int(etamax*10), npar))
	    Teffseta_veto2 = makeEffplot_eta([GEMCSCTrackCh0GE11, GEMCSCTrackCh0GE11, GEMCSCTrackCh0], "abs(eta_st2_sim)", [cutsbase+checkx+"&& pt>=%d"%(pt+5), cutsbase+checkx+checky_noGE21+"&& pt>=%d"%(pt+5), cutsbase+checkx+checky+"&& pt>=%d"%(pt+5)], ["fabs(%s)<%f"%(xaxis, ddY_cut)+"&& !isL1MediumVeto"+L1Mucut, ellipse_noGE21+"&& !isL1MediumVeto"+L1Mucut, ellipse+"&& !isL1MediumVeto"+L1Mucut], etamin, etamax,"|#eta| at 2nd CSC station","Trigger efficiency",legs3,"medium track veto",text,SBPlots+"GEMCSC_dxy10_50_mediumveto_combined_effvseta_20170131_pt%d_fraction%d_st2eta%dto%d_npar%d"%(pt, fraction, int(etamin*10),int(etamax*10), npar))
	    Teffsdxy = makeEffplot_dxy([GEMCSCTrackCh0, GEMCSCTrackCh0GE11, GEMCSCTrackCh0], "abs(genGdMu_dxy)", [cutsbase_nodxycut+checkx+"&& pt>=%d"%(pt+5), cutsbase_nodxycut+checkx+checky_noGE21+"&& pt>=%d"%(pt+5), cutsbase_nodxycut+checkx+checky+"&& pt>=%d"%(pt+5)], ["fabs(%s)<%f"%(xaxis, ddY_cut)+L1Mucut, ellipse_noGE21+L1Mucut, ellipse+L1Mucut], etamin, etamax,"|#eta| at 2nd CSC station","Trigger efficiency",legs3,"no track veto",text,SBPlots+"GEMCSC_dxy10_50_combined_effvsdxy_20170131_pt%d_fraction%d_st2eta%dto%d_npar%d"%(pt, fraction, int(etamin*10),int(etamax*10), npar))
	    Teffsdxy_veto2 = makeEffplot_dxy([GEMCSCTrackCh0GE11, GEMCSCTrackCh0GE11, GEMCSCTrackCh0], "abs(genGdMu_dxy)", [cutsbase_nodxycut+checkx+"&& pt>=%d"%(pt+5), cutsbase_nodxycut+checkx+checky_noGE21+"&& pt>=%d"%(pt+5), cutsbase_nodxycut+checkx+checky+"&& pt>=%d"%(pt+5)], ["fabs(%s)<%f"%(xaxis, ddY_cut)+"&& !isL1MediumVeto"+L1Mucut, ellipse_noGE21+"&& !isL1MediumVeto"+L1Mucut, ellipse+"&& !isL1MediumVeto"+L1Mucut], etamin, etamax,"|#eta| at 2nd CSC station","Trigger efficiency",legs3,"medium track veto",text,SBPlots+"GEMCSC_dxy10_50_mediumveto_combined_effvsdxy_20170131_pt%d_fraction%d_st2eta%dto%d_npar%d"%(pt, fraction, int(etamin*10),int(etamax*10), npar))
	    peff, perr = getOverallEff(Teffs[0], pt)
	    heff, herr = getOverallEff(Teffs[2], pt)
	    heff_noGE21, herr_noGE21 = getOverallEff(Teffs[1], pt)
	    peff_v2, perr_v2 = getOverallEff(Teffs_veto2[0], pt)
	    heff_v2, herr_v2 = getOverallEff(Teffs_veto2[2], pt)
	    heff_noGE21_v2, herr_noGE21_v2 = getOverallEff(Teffs_veto2[1], pt)
	    print "pt ",pt," fraction ",fraction," etamin ",etamin, " etamax ",etamax," npar ",npar, "peff ",peff, " heff ",heff," after veto peff2 ",peff_v2," heff2 ",heff_v2," prate ",ratehist_p_veto0.GetEntries(), " hrate ",ratehist_h_veto0.GetEntries()," after veto prate2 ",ratehist_p_veto2.GetEntries()," hrate2 ",ratehist_h_veto2.GetEntries()
	    #makeplotsPt([Teffs_ddY[0], Teffs_dphi_dir[0], Teffs[0]], legs3, text,SBPlots+"GEMCSC_ctau0andctau1000_combined3_eff_20170131_pt%d_ptbg%d_fraction%d_fractiond%d_St2eta%dto%d_npar%d_dirinterstation%s"%(pt, pt1, fraction, fraction_d, int(etamin*10),int(etamax*10), npar, extrapic)) 
	    Teffallnpars.append(Teffs)
	    Teffallnpars_veto2.append(Teffs_veto2)
	    Teffallnparseta.append(Teffseta)
	    Teffallnparseta_veto2.append(Teffseta_veto2)
	    Teffallnparsdxy.append(Teffsdxy)
	    Teffallnparsdxy_veto2.append(Teffsdxy_veto2)
	    #Teffallnpars_d.append(Teffs)
	    #Teffallnpars_d.append(Teffs_dphi_dir)
	    Teffs[0].SetName("positioneta%dto%dnpar%d"%(int(etamin*10),int(etamax*10), npar)+"displacedmuonpt%d"%(pt)+"fraction%d"%fraction)
	    Teffs[2].SetName("hybrideta%dto%dnpar%d"%(int(etamin*10),int(etamax*10), npar)+"displacedmuonpt%d"%(pt)+"fraction%d"%fraction)
	    Teffs[1].SetName("hybrid_noGE21eta%dto%dnpar%d"%(int(etamin*10),int(etamax*10), npar)+"displacedmuonpt%d"%(pt)+"fraction%d"%fraction)
	    Teffs_veto2[0].SetName("positionmediumvetoeta%dto%dnpar%d"%(int(etamin*10),int(etamax*10), npar)+"displacedmuonpt%d"%(pt)+"fraction%d"%fraction)
	    Teffs_veto2[2].SetName("hybridmediumvetoeta%dto%dnpar%d"%(int(etamin*10),int(etamax*10), npar)+"displacedmuonpt%d"%(pt)+"fraction%d"%fraction)
	    Teffs_veto2[1].SetName("hybrid_noGE21mediumvetoeta%dto%dnpar%d"%(int(etamin*10),int(etamax*10), npar)+"displacedmuonpt%d"%(pt)+"fraction%d"%fraction)
	    Teffseta[0].SetName("positioneta%dto%dnpar%d"%(int(etamin*10),int(etamax*10), npar)+"displacedmuonpt%d"%(pt)+"fraction%deta"%fraction)
	    Teffseta[2].SetName("hybrideta%dto%dnpar%d"%(int(etamin*10),int(etamax*10), npar)+"displacedmuonpt%d"%(pt)+"fraction%deta"%fraction)
	    Teffseta[1].SetName("hybrid_noGE21eta%dto%dnpar%d"%(int(etamin*10),int(etamax*10), npar)+"displacedmuonpt%d"%(pt)+"fraction%deta"%fraction)
	    Teffseta_veto2[0].SetName("positionmediumvetoeta%dto%dnpar%d"%(int(etamin*10),int(etamax*10), npar)+"displacedmuonpt%d"%(pt)+"fraction%deta"%fraction)
	    Teffseta_veto2[2].SetName("hybridmediumvetoeta%dto%dnpar%d"%(int(etamin*10),int(etamax*10), npar)+"displacedmuonpt%d"%(pt)+"fraction%deta"%fraction)
	    Teffseta_veto2[1].SetName("hybrid_noGE21mediumvetoeta%dto%dnpar%d"%(int(etamin*10),int(etamax*10), npar)+"displacedmuonpt%d"%(pt)+"fraction%deta"%fraction)
	    #Teffs_ddY[0].SetName("poistioneta%dto%dnpar%d"%(int(etamin*10),int(etamax*10), npar)+"promptmuonpt%dptbg%d"%(pt,pt1)+"fraction%d"%fraction+extrapic)
	    #Teffs_ddY[0].SetName("poistioneta%dto%dnpar%d"%(int(etamin*10),int(etamax*10), npar)+"promptmuonpt%dptbg%d"%(pt,pt1)+"fraction%d"%fraction+extrapic)
	    #Teffs_dphi_dir[0].SetName("directioneta%dto%dnpar%d"%(int(etamin*10),int(etamax*10), npar)+"promptmuonpt%dptbg%d"%(pt,pt1)+"fraction%d"%fraction+extrapic)
	    #Teffs[1].SetName("hybrideta%dto%dnpar%d"%(int(etamin*10),int(etamax*10), npar)+"displacedmuonpt%dptbg%d"%(pt,pt1)+"fraction%d"%fraction+extrapic)
	    #Teffs_ddY[1].SetName("poistioneta%dto%dnpar%d"%(int(etamin*10),int(etamax*10), npar)+"displacedmuonpt%dptbg%d"%(pt,pt1)+"fraction%d"%fraction+extrapic)
	    #Teffs_dphi_dir[1].SetName("directioneta%dto%dnpar%d"%(int(etamin*10),int(etamax*10), npar)+"displacedmuonpt%dptbg%d"%(pt,pt1)+"fraction%d"%fraction+extrapic)
	    Teffs[0].Write();Teffs[1].Write();Teffs[2].Write()
	    Teffs_veto2[0].Write();Teffs_veto2[1].Write();Teffs_veto2[2].Write()
	    Teffseta[0].Write();Teffseta[1].Write();Teffseta[2].Write()
	    Teffseta_veto2[0].Write();Teffseta_veto2[1].Write();Teffseta_veto2[2].Write()
	    #Teffs_ddY[0].Write()
	    #Teffs_dphi_dir[0].Write()
	    #Teffs_ddY[1].Write()
	    #Teffs_dphi_dir[1].Write()
       #print "Teffallnpars len ",len(Teffallnpars),Teffallnpars
       Teffeta0_p = Teffallnparseta[0][0]; Teffeta0_h = Teffallnparseta[0][2];Teffeta0_h_noGE21 = Teffallnparseta[0][1]
       Teffeta2_p = Teffallnparseta_veto2[0][0]; Teffeta2_h = Teffallnparseta_veto2[0][2];Teffeta2_h_noGE21 = Teffallnparseta_veto2[0][1]
       Teffdxy0_p = Teffallnparsdxy[0][0]; Teffdxy0_h = Teffallnparsdxy[0][2];Teffdxy0_h_noGE21 = Teffallnparsdxy[0][1]
       Teffdxy2_p = Teffallnparsdxy_veto2[0][0]; Teffdxy2_h = Teffallnparsdxy_veto2[0][2];Teffdxy2_h_noGE21 = Teffallnparsdxy_veto2[0][1]
       Teff0_p = Teffallnpars[0][0]; Teff0_h = Teffallnpars[0][2];Teff0_h_noGE21 = Teffallnpars[0][1]
       Teff2_p = Teffallnpars_veto2[0][0]; Teff2_h = Teffallnpars_veto2[0][2];Teff2_h_noGE21 = Teffallnpars_veto2[0][1]
       hrate0_p = hrateallnpars[0][0]; hrate0_h = hrateallnpars[0][2]; hrate0_h_noGE21 = hrateallnpars[0][1]
       hrate2_p = hrateallnpars_veto2[0][0]; hrate2_h = hrateallnpars_veto2[0][2]; hrate2_h_noGE21 = hrateallnpars_veto2[0][1]
       for xpar in range(len(Teffallnpars)-1):
	       
	    Teffeta0_p.Add(Teffallnparseta[xpar+1][0])		   
	    Teffeta0_h.Add(Teffallnparseta[xpar+1][2])		   
	    Teffeta0_h_noGE21.Add(Teffallnparseta[xpar+1][1])		   
	    Teffeta2_p.Add(Teffallnparseta_veto2[xpar+1][0])		   
	    Teffeta2_h.Add(Teffallnparseta_veto2[xpar+1][2])		   
	    Teffeta2_h_noGE21.Add(Teffallnparseta_veto2[xpar+1][1])		   
	       
	    Teffdxy0_p.Add(Teffallnparsdxy[xpar+1][0])		   
	    Teffdxy0_h.Add(Teffallnparsdxy[xpar+1][2])		   
	    Teffdxy0_h_noGE21.Add(Teffallnparsdxy[xpar+1][1])		   
	    Teffdxy2_p.Add(Teffallnparsdxy_veto2[xpar+1][0])		   
	    Teffdxy2_h.Add(Teffallnparsdxy_veto2[xpar+1][2])		   
	    Teffdxy2_h_noGE21.Add(Teffallnparsdxy_veto2[xpar+1][1])		   
	       
	    Teff0_p.Add(Teffallnpars[xpar+1][0])		   
	    Teff0_h.Add(Teffallnpars[xpar+1][2])		   
	    Teff0_h_noGE21.Add(Teffallnpars[xpar+1][1])		   
	    Teff2_p.Add(Teffallnpars_veto2[xpar+1][0])		   
	    Teff2_h.Add(Teffallnpars_veto2[xpar+1][2])		   
	    Teff2_h_noGE21.Add(Teffallnpars_veto2[xpar+1][1])		   

	    hrate0_p.Add(hrateallnpars[xpar+1][0])
	    hrate0_h.Add(hrateallnpars[xpar+1][2])
	    hrate0_h_noGE21.Add(hrateallnpars[xpar+1][1])
	    hrate2_p.Add(hrateallnpars_veto2[xpar+1][0])
	    hrate2_h.Add(hrateallnpars_veto2[xpar+1][2])
	    hrate2_h_noGE21.Add(hrateallnpars_veto2[xpar+1][1])

       Teffalleta.append([Teffeta0_p, Teffeta0_h_noGE21, Teffeta0_h])
       Teffalleta2.append([Teffeta2_p, Teffeta2_h_noGE21, Teffeta2_h])
       Teffalldxy.append([Teffdxy0_p, Teffdxy0_h_noGE21, Teffdxy0_h])
       Teffalldxy2.append([Teffdxy2_p, Teffdxy2_h_noGE21, Teffdxy2_h])
       Teffall.append([Teff0_p, Teff0_h_noGE21, Teff0_h])
       Teffall2.append([Teff2_p, Teff2_h_noGE21, Teff2_h])
       hrateall.append([hrate0_p, hrate0_h_noGE21, hrate0_h])
       hrateall2.append([hrate2_p, hrate2_h_noGE21, hrate2_h])
       peff0, perr0 = getOverallEff(Teff0_p, pt)
       peff2, perr2 = getOverallEff(Teff2_p, pt)
       heff0_noGE21, herr0_noGE21 = getOverallEff(Teff0_h_noGE21, pt)
       heff2_noGE21, herr2_noGE21 = getOverallEff(Teff2_h_noGE21, pt)
       heff0, herr0 = getOverallEff(Teff0_h, pt)
       heff2, herr2 = getOverallEff(Teff2_h, pt)
       #print "pt ",pt," fraction ",fraction," etamin ",etamin, " etamax ",etamax," allnpar peff0 ",peff0, " heff0 ",heff0," heff0_noGE21 ",heff0_noGE21
       #print "Rate entreis p ",hrate0_p.GetEntries()," h_noGE21 ",hrate0_h_noGE21.GetEntries()," h ",hrate0_h.GetEntries()
       #print "all eff, medium veto p ",peff2," heff2_noGE21 ",heff2_noGE21," heff2 ",heff2
       #print "Rate, medium veto,  entreis p ",hrate2_p.GetEntries()," h_noGE21 ",hrate2_h_noGE21.GetEntries()," h ",hrate2_h.GetEntries()
       Teffeta0_p.SetName("positioneta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%dvseta"%(fraction, pt))
       Teffeta0_h.SetName("hybrideta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%dvseta"%(fraction, pt))
       Teffeta0_h_noGE21.SetName("hybrid_noGE21eta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%dvseta"%(fraction, pt))
       Teffeta2_p.SetName("positionmediumvetoeta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%dvseta"%(fraction, pt))
       Teffeta2_h.SetName("hybridmediumvetoeta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%dvseta"%(fraction, pt))
       Teffeta2_h_noGE21.SetName("hybrid_noGE21mediumvetoeta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%dvseta"%(fraction, pt))
       Teffeta0_p.Write(); Teffeta0_h.Write(); Teffeta0_h_noGE21.Write()
       Teffeta2_p.Write(); Teffeta2_h.Write(); Teffeta2_h_noGE21.Write()

       Teffdxy0_p.SetName("positioneta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%dvsdxy"%(fraction, pt))
       Teffdxy0_h.SetName("hybrideta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%dvsdxy"%(fraction, pt))
       Teffdxy0_h_noGE21.SetName("hybrid_noGE21eta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%dvsdxy"%(fraction, pt))
       Teffdxy2_p.SetName("positionmediumvetoeta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%dvsdxy"%(fraction, pt))
       Teffdxy2_h.SetName("hybridmediumvetoeta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%dvsdxy"%(fraction, pt))
       Teffdxy2_h_noGE21.SetName("hybrid_noGE21mediumvetoeta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%dvsdxy"%(fraction, pt))
       Teffdxy0_p.Write(); Teffdxy0_h.Write(); Teffdxy0_h_noGE21.Write()
       Teffdxy2_p.Write(); Teffdxy2_h.Write(); Teffdxy2_h_noGE21.Write()

       Teff0_p.SetName("positioneta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%d"%(fraction, pt))
       Teff0_h.SetName("hybrideta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%d"%(fraction, pt))
       Teff0_h_noGE21.SetName("hybrid_noGE21eta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%d"%(fraction, pt))
       Teff2_p.SetName("positionmediumvetoeta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%d"%(fraction, pt))
       Teff2_h.SetName("hybridmediumvetoeta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%d"%(fraction, pt))
       Teff2_h_noGE21.SetName("hybrid_noGE21mediumvetoeta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%d"%(fraction, pt))
       Teff0_p.Write(); Teff0_h.Write(); Teff0_h_noGE21.Write()
       Teff2_p.Write(); Teff2_h.Write(); Teff2_h_noGE21.Write()

       hrate0_p.SetName("ratehist_position_eta%dto%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10),  fraction, pt))
       hrate2_p.SetName("ratehist_position_mediumveto_eta%dto%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), fraction, pt))
       hrate0_h_noGE21.SetName("ratehist_hybrid_noGE21_eta%dto%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), fraction, pt))
       hrate2_h_noGE21.SetName("ratehist_hybrid_noGE21_mediumveto_eta%dto%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), fraction, pt))
       hrate0_h.SetName("ratehist_hybrid_eta%dto%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), fraction, pt))
       hrate2_h.SetName("ratehist_hybrid_mediumveto_eta%dto%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), fraction, pt))
       hrate0_p.Write(); hrate0_h_noGE21.Write(); hrate0_h.Write()
       hrate2_p.Write(); hrate2_h_noGE21.Write(); hrate2_h.Write()
       text_h = "#splitline{Hybrid algorithm w/ GE21}{%.2f<|#eta|<%.2f, p_{T}^{L1}>%d GeV}"%(etamin,etamax, pt)
       text_p = "#splitline{Position-based algorithm}{%.2f<|#eta|<%.2f, p_{T}^{L1}>%d GeV}"%(etamin,etamax, pt)
       text_h_noGE21 = "#splitline{Hybrid algorithm w/o GE21}{%.2f<|#eta|<%.2f, p_{T}^{L1}>%d GeV}"%(etamin,etamax, pt)
       legsveto = ["no track veto","medium track veto"]
       #makeplotsPt([Teff0_p, Teff2_p], legsveto, text_p,SBPlots+"GEMCSC_dxy10_50_position_vetos_eff_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10))) 
       #makeplotsPt([Teff0_h, Teff2_h], legsveto, text_h,SBPlots+"GEMCSC_dxy10_50_hybrid_vetos_eff_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10))) 
       #makeplotsPt([Teff0_h_noGE21, Teff2_h_noGE21], legsveto, text_h_noGE21,SBPlots+"GEMCSC_dxy10_50_hybrid_noGE21_vetos_eff_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10))) 
       #plotshistograms([hrate0_p, hrate2_p], legsveto, text_p, SBPlots+"TriggerRate_position_vetos_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10)))
       #plotshistograms([hrate0_h, hrate2_h], legsveto, text_h, SBPlots+"TriggerRate_hybrid_vetos_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10)))
       #plotshistograms([hrate0_h_noGE21, hrate2_h_noGE21], legsveto, text_h_noGE21, SBPlots+"TriggerRate_hybrid_noGE21_vetos_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10)))
       text_noveto =  "%.2f<|#eta|<%.2f, p_{T}^{L1}>%d GeV, no track veto"%(etamin,etamax, pt)
       text_mveto =  "%.2f<|#eta|<%.2f, p_{T}^{L1}>%d GeV, medium track veto"%(etamin,etamax, pt)
       #makeplotsPt([Teff0_p, Teff0_h_noGE21, Teff0_h], legs3, text_noveto,SBPlots+"GEMCSC_dxy10_50_combined3_noveto_eff_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10))) 
       #makeplotsPt([Teff2_p, Teff2_h_noGE21, Teff2_h], legs3, text_mveto,SBPlots+"GEMCSC_dxy10_50_combined3_mediumveto_eff_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10))) 
       #plotshistograms([hrate0_p, hrate0_h_noGE21, hrate0_h], legs3, text_noveto, SBPlots+"TriggerRate_combiend_noveto_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10)))
       #plotshistograms([hrate2_p, hrate2_h_noGE21, hrate2_h], legs3, text_mveto, SBPlots+"TriggerRate_combiend_mediumveto_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10)))
    etamin = netas[0]; etamax = netas[-1]
    Teffalleta0_p = Teffalleta[0][0]; Teffalleta0_h = Teffalleta[0][2];Teffalleta0_h_noGE21 = Teffalleta[0][1]
    Teffalleta2_p = Teffalleta2[0][0]; Teffalleta2_h = Teffalleta2[0][2];Teffalleta2_h_noGE21 = Teffalleta2[0][1]
    Teffalldxy0_p = Teffalldxy[0][0]; Teffalldxy0_h = Teffalldxy[0][2];Teffalldxy0_h_noGE21 = Teffalldxy[0][1]
    Teffalldxy2_p = Teffalldxy2[0][0]; Teffalldxy2_h = Teffalldxy2[0][2];Teffalldxy2_h_noGE21 = Teffalldxy2[0][1]
    Teffall0_p = Teffall[0][0]; Teffall0_h = Teffall[0][2];Teffall0_h_noGE21 = Teffall[0][1]
    Teffall2_p = Teffall2[0][0]; Teffall2_h = Teffall2[0][2];Teffall2_h_noGE21 = Teffall2[0][1]
    hrateall0_p = hrateall[0][0]; hrateall0_h = hrateall[0][2]; hrateall0_h_noGE21 = hrateall[0][1]
    hrateall2_p = hrateall2[0][0]; hrateall2_h = hrateall2[0][2]; hrateall2_h_noGE21 = hrateall2[0][1]

    for xeta in range(len(Teffall)-1):

	Teffalleta0_p.Add(Teffalleta[xeta+1][0])		   
	Teffalleta0_h_noGE21.Add(Teffalleta[xeta+1][1])		   
	Teffalleta0_h.Add(Teffalleta[xeta+1][2])		   
	Teffalleta2_p.Add(Teffalleta2[xeta+1][0])		   
	Teffalleta2_h_noGE21.Add(Teffalleta2[xeta+1][1])		   
	Teffalleta2_h.Add(Teffalleta2[xeta+1][2])		   

	Teffalldxy0_p.Add(Teffalldxy[xeta+1][0])		   
	Teffalldxy0_h_noGE21.Add(Teffalldxy[xeta+1][1])		   
	Teffalldxy0_h.Add(Teffalldxy[xeta+1][2])		   
	Teffalldxy2_p.Add(Teffalldxy2[xeta+1][0])		   
	Teffalldxy2_h_noGE21.Add(Teffalldxy2[xeta+1][1])		   
	Teffalldxy2_h.Add(Teffalldxy2[xeta+1][2])		   

	Teffall0_p.Add(Teffall[xeta+1][0])		   
	Teffall0_h_noGE21.Add(Teffall[xeta+1][1])		   
	Teffall0_h.Add(Teffall[xeta+1][2])		   
	Teffall2_p.Add(Teffall2[xeta+1][0])		   
	Teffall2_h_noGE21.Add(Teffall2[xeta+1][1])		   
	Teffall2_h.Add(Teffall2[xeta+1][2])		   

        hrateall0_p.Add(hrateall[xeta+1][0])
        hrateall0_h_noGE21.Add(hrateall[xeta+1][1])
        hrateall0_h.Add(hrateall[xeta+1][2])
        hrateall2_p.Add(hrateall2[xeta+1][0])
        hrateall2_h_noGE21.Add(hrateall2[xeta+1][1])
        hrateall2_h.Add(hrateall2[xeta+1][2])

    peff0, perr0 = getOverallEff(Teffall0_p, pt+5)
    peff2, perr2 = getOverallEff(Teffall2_p, pt+5)
    heff0_noGE21, herr0_noGE21 = getOverallEff(Teffall0_h_noGE21, pt+5)
    heff2_noGE21, herr2_noGE21 = getOverallEff(Teffall2_h_noGE21, pt+5)
    heff0, herr0 = getOverallEff(Teffall0_h, pt+5)
    heff2, herr2 = getOverallEff(Teffall2_h, pt+5)
     
    print "pt ",pt," fraction ",fraction," etamin ",etamin, " etamax ",etamax," allnpar peff0 ",peff0, " heff0 ",heff0," heff0_noGE21 ",heff0_noGE21
    print "Rate entreis p ",hrateall0_p.GetEntries()," h_noGE21 ",hrateall0_h_noGE21.GetEntries()," h ",hrateall0_h.GetEntries()
    print "all eff, medium veto p ",peff2," heff2_noGE21 ",heff2_noGE21," heff2 ",heff2
    print "Rate, medium veto,  entreis p ",hrateall2_p.GetEntries()," h_noGE21 ",hrateall2_h_noGE21.GetEntries()," h ",hrateall2_h.GetEntries()

    Teffalleta0_p.SetName("positioneta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%dvseta"%(fraction, pt))
    Teffalleta0_h.SetName("hybrideta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%dvseta"%(fraction, pt))
    Teffalleta0_h_noGE21.SetName("hybrid_noGE21eta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%dvseta"%(fraction, pt))
    Teffalleta2_p.SetName("positionmediumvetoeta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%dvseta"%(fraction, pt))
    Teffalleta2_h.SetName("hybridmediumvetoeta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%dvseta"%(fraction, pt))
    Teffalleta2_h_noGE21.SetName("hybrid_noGE21mediumvetoeta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%dvseta"%(fraction, pt))
    Teffalleta0_p.Write(); Teffalleta0_h.Write(); Teffalleta0_h_noGE21.Write()
    Teffalleta2_p.Write(); Teffalleta2_h.Write(); Teffalleta2_h_noGE21.Write()

    Teffalldxy0_p.SetName("positioneta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%dvsdxy"%(fraction, pt))
    Teffalldxy0_h.SetName("hybrideta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%dvsdxy"%(fraction, pt))
    Teffalldxy0_h_noGE21.SetName("hybrid_noGE21eta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%dvsdxy"%(fraction, pt))
    Teffalldxy2_p.SetName("positionmediumvetoeta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%dvsdxy"%(fraction, pt))
    Teffalldxy2_h.SetName("hybridmediumvetoeta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%dvsdxy"%(fraction, pt))
    Teffalldxy2_h_noGE21.SetName("hybrid_noGE21mediumvetoeta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%dvsdxy"%(fraction, pt))
    Teffalldxy0_p.Write(); Teffalldxy0_h.Write(); Teffalldxy0_h_noGE21.Write()
    Teffalldxy2_p.Write(); Teffalldxy2_h.Write(); Teffalldxy2_h_noGE21.Write()

    Teffall0_p.SetName("positioneta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%d"%(fraction, pt))
    Teffall0_h.SetName("hybrideta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%d"%(fraction, pt))
    Teffall0_h_noGE21.SetName("hybrid_noGE21eta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%d"%(fraction, pt))
    Teffall2_p.SetName("positionmediumvetoeta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%d"%(fraction, pt))
    Teffall2_h.SetName("hybridmediumvetoeta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%d"%(fraction, pt))
    Teffall2_h_noGE21.SetName("hybrid_noGE21mediumvetoeta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%d"%(fraction, pt))
    Teffall0_p.Write(); Teffall0_h.Write(); Teffall0_h_noGE21.Write()
    Teffall2_p.Write(); Teffall2_h.Write(); Teffall2_h_noGE21.Write()

    hrateall0_p.SetName("ratehist_position_eta%dto%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), fraction, pt))
    hrateall2_p.SetName("ratehist_position_mediumveto_eta%dto%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), fraction, pt))
    hrateall0_h_noGE21.SetName("ratehist_hybrid_noGE21_eta%dto%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), fraction, pt))
    hrateall2_h_noGE21.SetName("ratehist_hybrid_noGE21_mediumveto_eta%dto%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), fraction, pt))
    hrateall0_h.SetName("ratehist_hybrid_eta%dto%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), fraction, pt))
    hrateall2_h.SetName("ratehist_hybrid_mediumveto_eta%dto%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), fraction, pt))
    hrateall0_p.Write(); hrateall0_h_noGE21.Write(); hrateall0_h.Write()
    hrateall2_p.Write(); hrateall2_h_noGE21.Write(); hrateall2_h.Write()

    text_h = "#splitline{Hybrid algorithm w/ GE21}{%.2f<|#eta|<%.2f, p_{T}^{L1}>%d GeV}"%(etamin,etamax, pt)
    text_p = "#splitline{Position-based algorithm}{%.2f<|#eta|<%.2f, p_{T}^{L1}>%d GeV}"%(etamin,etamax, pt)
    text_h_noGE21 = "#splitline{Hybrid algorithm w/o GE21}{%.2f<|#eta|<%.2f, p_{T}^{L1}>%d GeV}"%(etamin,etamax, pt)
    legsveto = ["no track veto","medium track veto"]
    legheader = "L1Mu with Q>=4"
    makeplotsPt([Teffall0_p, Teffall2_p], legsveto, legheader, text_p,SBPlots+"GEMCSC_dxy10_50_position_vetos_effvspt_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10))) 
    makeplotsPt([Teffall0_h, Teffall2_h], legsveto, legheader,text_h,SBPlots+"GEMCSC_dxy10_50_hybrid_vetos_effvspt_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10))) 
    makeplotsPt([Teffall0_h_noGE21, Teffall2_h_noGE21], legsveto, legheader,text_h_noGE21,SBPlots+"GEMCSC_dxy10_50_hybrid_noGE21_vetos_effvspt_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10))) 
    makeplotsEta([Teffalleta0_p, Teffalleta2_p], legsveto, legheader, text_p,SBPlots+"GEMCSC_dxy10_50_position_vetos_effvseta_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10))) 
    makeplotsEta([Teffalleta0_h, Teffalleta2_h], legsveto, legheader,text_h,SBPlots+"GEMCSC_dxy10_50_hybrid_vetos_effvseta_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10))) 
    makeplotsDxy([Teffalldxy0_p, Teffalldxy2_p], legsveto, legheader, text_p,SBPlots+"GEMCSC_dxy10_50_position_vetos_effvsdxy_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10))) 
    makeplotsDxy([Teffalldxy0_h, Teffalldxy2_h], legsveto, legheader,text_h,SBPlots+"GEMCSC_dxy10_50_hybrid_vetos_effvsdxy_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10))) 

    plotshistograms([hrateall0_p, hrateall2_p], legsveto,"L1Mu with Q>=4", text_p, SBPlots+"TriggerRate_position_vetos_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10)))
    plotshistograms([hrateall0_h, hrateall2_h], legsveto, "L1Mu with Q>=4", text_h, SBPlots+"TriggerRate_hybrid_vetos_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10)))
    plotshistograms([hrateall0_h_noGE21, hrateall2_h_noGE21], legsveto, "L1Mu with Q>=4", text_h_noGE21, SBPlots+"TriggerRate_hybrid_noGE21_vetos_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10)))

    text_noveto =  "%.2f<|#eta|<%.2f, p_{T}^{L1}>%d GeV, no track veto"%(etamin,etamax, pt)
    text_mveto =  "%.2f<|#eta|<%.2f, p_{T}^{L1}>%d GeV, medium track veto"%(etamin,etamax, pt)
    legheader = "L1Mu with Q>=4, no track veto"
    legheaderm = "L1Mu with Q>=4, medium track veto"
    legs2 = ["CSC+GE11: Position based","CSC+GE11+GE21: Hybrid algorithm"]
    makeplotsEta([Teffalleta0_p,  Teffalleta0_h], legs2, legheader,text_noveto,SBPlots+"GEMCSC_dxy10_50_combined2_noveto_effvseta_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10))) 
    makeplotsEta([Teffalleta2_p,  Teffalleta2_h], legs2, legheaderm,text_mveto,SBPlots+"GEMCSC_dxy10_50_combined2_mediumveto_effvseta_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10))) 
    makeplotsDxy([Teffalldxy0_p,  Teffalldxy0_h], legs2, legheader,text_noveto,SBPlots+"GEMCSC_dxy10_50_combined2_noveto_effvsdxy_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10))) 
    makeplotsDxy([Teffalldxy2_p,  Teffalldxy2_h], legs2, legheaderm,text_mveto,SBPlots+"GEMCSC_dxy10_50_combined2_mediumveto_effvsdxy_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10))) 
    makeplotsPt([Teffall0_p, Teffall0_h], legs2, legheader,text_noveto,SBPlots+"GEMCSC_dxy10_50_combined2_noveto_effvspt_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10))) 
    makeplotsPt([Teffall2_p, Teffall2_h], legs2, legheaderm,text_mveto,SBPlots+"GEMCSC_dxy10_50_combined2_mediumveto_effvspt_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10))) 

    plotshistograms([hrateall0_p, hrateall0_h], legs2, legheader, text_noveto, SBPlots+"TriggerRate_combined2_noveto_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10)))
    plotshistograms([hrateall2_p, hrateall2_h], legs2, legheaderm, text_mveto, SBPlots+"TriggerRate_combined2_mediumveto_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10)))
    return peff0, perr0, heff0_noGE21, herr0_noGE21, heff0, herr0, peff2, perr2, heff2_noGE21, herr2_noGE21,heff2,herr2, hrateall0_p.Integral(),hrateall0_h_noGE21.Integral(), hrateall0_h.Integral(), hrateall2_p.Integral(), hrateall2_h_noGE21.Integral(),hrateall2_h.Integral()

    Teff_out.Close()


def plotalletaME0GE21(pt, netas, fraction, outputrootfile, ratesample, signalsample, allLUTfile, SBPlots):
    rate_y1 = "dphi_dir12_noME0_noGE21"
    rate_y2 = "dphi_dir12_noME0_GE21"
    rate_y3 = "dphi_dir12_GE21"
    rate_x = "ddY123"
    xaxis = "(deltay123_test)"
    yaxis1 = "(dphi_dir_ME11only_ME21only_L1)"
    yaxis2 = "(dphi_dir_ME11only_st2_L1)"
    yaxis3 = "(dphi_dir_st1_st2_L1)"
    legs3_v2 = ["Position-based algorithm","Stub alignment algorithm","Hybrid algorithm"]
    legs3 = ["Position-based algorithm","Hybrid algorithm w/o GE21","Hybrid algorithm w GE21"]
    legsall = ["Position-based algorithm","Hybrid algirthm+ w/o ME11+ w/o GE21","Hybrid algorithm+ w/o ME0 + w GE21"," Hybrid algorithm+ w ME0 + w GE21"]
    legsall = ["Position-based algorithm","Stub alignment algorithm+ w/o ME0 + w/o GE21","Stub alignment algorithm+ w/o ME0 + w GE21","Stub alignment algorithm+ w ME0 + w GE21"]


    ratefile = TFile(ratesample)
    #h_total = ratefile.Get("h_eventcount")

    rateEvents = total
    L1MuTriggerRate = TChain("L1MuTriggerRateME0Segment192")
    #L1MuTriggerRate = TChain("L1MuTriggerRate")
    addFilesToChain(L1MuTriggerRate, ratesample)
    GEMCSCTrackCh0 = TChain("GEMCSCAnalyzer/trk_eff_CSC_ALL")
    GEMCSCME31Ch0 = TChain("GEMCSCAnalyzer/trk_eff_CSC_ME31")
    GEMCSCTrackCh0GE11 = TChain("GEMCSCAnalyzer/trk_eff_CSC_ALL")
    GEMCSCME31Ch0GE11 = TChain("GEMCSCAnalyzer/trk_eff_CSC_ME31")
    #GEMCSCTrackCh1 = TChain("GEMCSCAnalyzer/trk_eff_CSC_ALL")
    #addFilesToChain(GEMCSCTrackCh1, signalsample) 
    addFilesToChain(GEMCSCTrackCh0, signalsample) 
    addFilesToChain(GEMCSCME31Ch0, signalsample) 
    addFilesToChain(GEMCSCTrackCh0GE11, signalsampleGE11) 
    addFilesToChain(GEMCSCME31Ch0GE11, signalsampleGE11) 
    GEMCSCTrackCh0.AddFriend(GEMCSCME31Ch0)
    GEMCSCTrackCh0GE11.AddFriend(GEMCSCME31Ch0GE11)

    Teff_out = ROOT.TFile(outputrootfile,"UPDATE")
    Teffall = []
    Teffall2 = []
    Teffalleta = []
    Teffalleta2 = []
    hrateall = []
    hrateall2 = []
    quality = 4
    L1Mucut = "&& L1Mu_quality>=%d"%(quality)
    etamin = netas[0]; etamax = netas[-1]
    ratecut0 = "maxPromptPt && abs(L1Mu_eta_st2)>=%f && abs(L1Mu_eta_st2)<%f"%(etamin, etamax)
    ratecut1 = "maxPositionPt&& npar>=0 && abs(L1Mu_eta_st2)>=%f && abs(L1Mu_eta_st2)<%f"%(etamin, etamax)
    ratecut2 = "maxPositionPt && L1Mu_me0_dR<0.2 && npar>=0 && abs(L1Mu_eta_st2)>=%f && abs(L1Mu_eta_st2)<%f"%(etamin, etamax)
    ratehist_step0 = getRateVsEta(L1MuTriggerRate, "ratehist_step0_eta%dto%d"%(int(etamin*10),int(etamax*10)), myetabin, ratecut0 + L1Mucut , 0)
    ratehist_step1 = getRateVsEta(L1MuTriggerRate, "ratehist_step1_eta%dto%d"%(int(etamin*10),int(etamax*10)), myetabin, ratecut1 + L1Mucut , 0)
    ratehist_step2 = getRateVsEta(L1MuTriggerRate, "ratehist_step2_eta%dto%d"%(int(etamin*10),int(etamax*10)), myetabin, ratecut2 + L1Mucut , 0)
    print "hrate events, step0 ",ratehist_step0.GetEntries(), " step1 ",ratehist_step1.GetEntries(), " step2 ",ratehist_step2.GetEntries()
    for neta in range(len(netas)-1):
       if (netas[neta]<0 or netas[neta+1]<0):
	    continue
       ring = 1
       if (netas[neta]<1.6):
	    ring = 2
    #if (netas[0]==1.6):
       etamin = netas[neta]
       etamax = netas[neta+1]
       Teffallnpars = []
       Teffallnpars_veto2 = []
       Teffallnparseta = []
       Teffallnparseta_veto2 = []
       hrateallnpars = []
       hrateallnpars_veto2 = []
       #Teffallnpars_p = []
       #Teffallnpars_d = []
       legs = ["#splitline{Displaced muon}{10<|d_{xy}|<25 cm}"]
       for npar in range(0, 4):
	    me11 = evenodds[npar].split(',')[0]
	    me21 = evenodds[npar].split(',')[1]
	    chambers = "ME1%d %s,ME2%d %s,ME3%d %s"%(ring, me11, ring, me21, ring, me21)
	    me31_n = 1
	    if me21 == "even":
		    me31_n = 2
	    ME31cut = " && ((trk_eff_CSC_ME31.has_lct&%d)>0)"%(me31_n)
	    #st_title = ["Prompt muon, 2<p_{T}<%d"%pt1, "Displaced Muon, 10<|d_{xy}|<50, p_{T}>%d"%pt]
	    algos_lut = getFinalLUTValue(pt, fraction, etamin, etamax, npar, "", allLUTfile)
	    ddY_cut = algos_lut["Position"]
	    dphi_cut3 = algos_lut["Direction"]
	    """
	    algos_lut_noGE21 = getFinalLUTValue(pt, fraction, etamin, etamax, npar, "_noGE21", allLUTfile)
	    dphi_dir_cut = algos_lut["Direction"]
	    maxa = algos_lut["Hybrid"][0] 
	    maxb = algos_lut["Hybrid"][1] 
	    alpha= algos_lut["Hybrid"][2]
	    x0 =  algos_lut["Hybrid"][3]
	    y0 =  algos_lut["Hybrid"][4]
	    maxa1 = maxa; maxb1 = maxb; alpha1 = alpha; x01 = x0; y01= y0
	    if ring==1:
		maxa1 = algos_lut_noGE21["Hybrid_noGE21"][0] 
		maxb1 = algos_lut_noGE21["Hybrid_noGE21"][1] 
		alpha1= algos_lut_noGE21["Hybrid_noGE21"][2]
		x01 =  algos_lut_noGE21["Hybrid_noGE21"][3]
		y01 =  algos_lut_noGE21["Hybrid_noGE21"][4]
	    """

	    #ratecut = "maxPositionPt && L1Mu_quality>= %d && abs(L1Mu_eta)>=%f && abs(L1Mu_eta)<%f && npar==%d"%(quality, etamin, etamax, npar)
	    L1Mucut = "&& L1Mu_quality>=4  && hasME2"
	    ratecut = "maxPositionPt&& abs(L1Mu_eta_st2)>=%f && abs(L1Mu_eta_st2)<%f && npar==%d"%(etamin, etamax, npar)+L1Mucut
	    #etacuts ="&& meRing==%d"%ring+"&& fabs(trk_eff_CSC_ME2%d.eta_lct_%s)>%f && fabs(trk_eff_CSC_ME2%d.eta_lct_%s)<%f"%(ring,me21,etamin,ring,me21,etamax)
	    checkvalue = "&& fabs(deltay23_test)<50 && fabs(deltay12_test)<50 && fabs(%s)>0 && fabs(%s)>0 && meRing==%d"%(xaxis, yaxis1,ring)
	    checkx = "&& fabs(%s)<60 && fabs(%s)>0 "%(xaxis, xaxis)
	    checky = "&& fabs(%s)<2.0 && fabs(%s)>0 "%(yaxis1, yaxis1)
	    dxycut = "&& fabs(genGdMu_dxy)>5 && fabs(genGdMu_dxy)<25 && genGdMu_lxy<50 && fabs(genGdMu_dR)<0.1"
	    #checky_noGE21 = "&& fabs(%s)<2.0 && fabs(%s)>0 "%(yaxis1, yaxis1)
	    #cutsbase =  ["meRing==%d && "%ring+"hasSt1St2St3 && fabs(eta_st2_L1)>%f && fabs(eta_st2_L1)<%f && npar_L1==%d"%(etamin,etamax, npar), "meRing==%d && "%ring+"hasSt1St2St3 && fabs(eta_st2_L1)>%f && fabs(eta_st2_L1)<%f && fabs(genGdMu_dxy)>10 && fabs(genGdMu_dxy)<50 && fabs(genGdMu_dR)<0.1 && npar_L1==%d"%(etamin,etamax, npar)]
	    #cutsbase = "isMe0Region && meRing==%d && "%ring+"hasSt1St2St3 && fabs(eta_st2_L1)>%f && fabs(eta_st2_L1)<%f && fabs(genGdMu_dxy)>5 && fabs(genGdMu_dxy)<25 && genGdMu_lxy<50 && fabs(genGdMu_dR)<0.1 && npar_L1==%d"%(etamin,etamax, npar) + L1Mucut

	    #cutsbase = "meRing==%d && "%ring+"hasSt1St2St3 && fabs(eta_st2_L1)>%f && fabs(eta_st2_L1)<%f && fabs(genGdMu_dxy)>10 && fabs(genGdMu_dxy)<50 && fabs(genGdMu_dR)<0.1 && npar_L1==%d"%(etamin,etamax, npar)
	    checkx = ""; checky = ""
	    cutsbase = "has_gmtCand>0 && isMe0Region && fabs(eta_st2_sim)>%f && fabs(eta_st2_sim)<%f && npar_sim==%d"%(etamin,etamax, npar)+ME31cut
	    cuts_p =  [cutsbase+checkx]
	    cuts_d =  [cutsbase+checky]
	    #cuts_d =  ["meRing==%d && "%ring+"hasSt1St2St3 && fabs(deltay12_fit)<50 &&  fabs(deltay23_fit)<50 && fabs(eta_st2_sh)>%f && fabs(eta_st2_sh)<%f && npar==%d"%(etamin,etamax, npar)+checky, "meRing==%d && "%ring+" hasSt1St2St3 && fabs(deltay12_fit)<50 &&  fabs(deltay23_fit)<50 && fabs(eta_st2_sh)>%f && fabs(eta_st2_sh)<%f && fabs(genGdMu_dxy)>10 && fabs(genGdMu_dxy)<50 && fabs(genGdMu_dR)<0.1 && npar==%d"%(etamin,etamax, npar)+checky]

	    text = "#splitline{%s}{%.1f<|#eta|<%.2f, p_{T}^{L1}>%d GeV}"%(chambers, etamin,etamax,pt)
	    text_position = "#splitline{%s}{%.1f<|#eta|<%.2f, p_{T}^{L1}>%d GeV}"%(chambers,etamin,etamax,pt)
	    #ddYDPhiCuts = [ddY_cut, dphi_dir_cut]

	    ytitle = "#Delta#Phi_{dir}"; xtitle = "#Delta#Delta Y"
	    
	    dphi_cut1 =10; dphi_cut2 = 10
	    #dphi_cut1 = getCut(GEMCSCTrackCh0, "abs("+ yaxis1 +")", cuts_d[0]+"&& pt>%f-1 && pt<%f+1"%(pt, pt), 0.0, 1.5, 0.001, fraction*1.0/100.0)
	    #dphi_cut2 = getCut(GEMCSCTrackCh0, "abs("+yaxis2+")", cuts_d[0]+"&& pt>%f-1 && pt<%f+1"%(pt, pt), 0.0, 1.5, 0.001, fraction*1.0/100.0)
	    #dphi_cut2 = getCut(GEMCSCTrackCh0, "abs("+yaxis2+")", cuts_d[0]+"&& pt>%f-1 && pt<%f+1"%(pt, pt), 0.0, 1.5, 0.001, fraction*1.0/100.0)
	    #dphi_cut3 = getCut(GEMCSCTrackCh0, "abs("+yaxis3+")", cuts_d[0]+"&& pt>%f-1 && pt<%f+1"%(pt, pt), 0.0, 1.5, 0.001, fraction*1.0/100.0)
	    #dphi_cut1 = getCut(GEMCSCTrackCh0, "abs("+ yaxis1 +")", cuts_d[0]+"&& pt>%f"%(pt)+"&& pt<%f"%(pt+10), 0.0, 1.5, 0.005, fraction*1.0/100.0)
	    #dphi_cut2 = getCut(GEMCSCTrackCh0, "abs("+yaxis2+")", cuts_d[0]+"&& pt>%f"%(pt)+"&& pt<%f"%(pt+10), 0.0, 1.5, 0.005, fraction*1.0/100.0)
	    #dphi_cut3 = getCut(GEMCSCTrackCh0, "abs("+yaxis3+")", cuts_d[0]+"&& pt>%f"%(pt)+"&& pt<%f"%(pt+10), 0.0, 1.5, 0.005, fraction*1.0/100.0)
	    #ddY_cut2 = getCut(GEMCSCTrackCh0, "abs("+xaxis+")", cuts_p[0]+"&& pt>%f"%(pt), 0.0, 10, 0.01, 0.99)
	   
	    
	    if ( ddY_cut<0):
		    ddY_cut = getCut(GEMCSCTrackCh0, "abs("+xaxis+")", cuts_p[0]+"&& pt>%f"%(pt)+"&& pt<%f"%(pt+10), 0.0, 40, 0.01, fraction*1.0/100.0)
	    ddY_cut2 = 2*ddY_cut
	    #dphicut1 = getcuts(GEMCSCTrackCh0, cuts_d[0], yaxis1, pt, npar,etamin, etamax, fraction, "(3000,0,1.50)", ytitle,"",SBPlots+"GEMCSC_ctau0_directionOnly_getcuts_recopt_20170611_ptcut%d_st2eta%dto%d_npar%d%s"%(pt, int(etamin*10), int(etamax*10), npar,"noME0_noGE21"))
	    #dphicut2 = getcuts(GEMCSCTrackCh0, cuts_d[0], yaxis2, pt, npar,etamin, etamax, fraction, "(3000,0,1.50)", ytitle,"",SBPlots+"GEMCSC_ctau0_directionOnly_getcuts_recopt_20170611_ptcut%d_st2eta%dto%d_npar%d%s"%(pt, int(etamin*10), int(etamax*10), npar, "noME0_GE21"))
	    #dphicut3 = getcuts(GEMCSCTrackCh0, cuts_d[0], yaxis3, pt, npar,etamin, etamax, fraction, "(3000,0,1.50)", ytitle,"",SBPlots+"GEMCSC_ctau0_directionOnly_getcuts_recopt_20170611_ptcut%d_st2eta%dto%d_npar%d%s"%(pt, int(etamin*10), int(etamax*10), npar, "ME0_GE21"))
	    #ddY_cut = getcuts(GEMCSCTrackCh1, cuts_p[0], xaxis, pt, npar,etamin, etamax, fraction, "(4000,0,40.0)", xtitle,"",SBPlots+"GEMCSC_ctau0_positionOnly_getcuts_recopt_20170611_ptcut%d_st2eta%dto%d_npar%d%s"%(pt, int(etamin*10), int(etamax*10), npar,extrapic))
	    #ddY_cut2 = getcuts(GEMCSCTrackCh1, cuts_p[0], xaxis, pt, npar,etamin, etamax, fraction+2, "(4000,0,40.0)", xtitle,"",SBPlots+"GEMCSC_ctau0_positionOnly_getcuts_recopt_20170131_ptcut%d_st2eta%dto%d_npar%d%s"%(pt, int(etamin*10), int(etamax*10), npar,extrapic))
	    print "dphicut1 ",dphi_cut1," dphicut2 ",dphi_cut2, " dphicut3 ",dphi_cut3, " ddY_cut ",ddY_cut," ddY_cut2 ",ddY_cut2
	    
	    ratehist_p_veto0 = getRateVsEta(L1MuTriggerRate, "ratehist_position_eta%dto%d_npar%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), npar, fraction, pt), myetabin, ratecut + L1Mucut + "&& abs(%s)<=%f"%(rate_x, ddY_cut), 0)
	    ratehist_h1_veto0 = getRateVsEta(L1MuTriggerRate, "ratehist_hybrid1_eta%dto%d_npar%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), npar, fraction, pt), myetabin, ratecut + L1Mucut + "&& abs(%s)<=%f"%(rate_y1, dphi_cut1)+ "&& abs(%s)<=%f"%(rate_x, ddY_cut2), 0)
	    ratehist_h2_veto0 = getRateVsEta(L1MuTriggerRate, "ratehist_hybrid2_eta%dto%d_npar%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), npar, fraction, pt), myetabin, ratecut + L1Mucut + "&& abs(%s)<=%f"%(rate_y2, dphi_cut2)+ "&& abs(%s)<=%f"%(rate_x, ddY_cut2), 0)
	    ratehist_h3_veto0 = getRateVsEta(L1MuTriggerRate, "ratehist_hybrid3_eta%dto%d_npar%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), npar, fraction, pt), myetabin, ratecut + L1Mucut + "&& abs(%s)<=%f"%(rate_y3, dphi_cut3)+ "&& abs(%s)<=%f"%(rate_x, ddY_cut2), 0)
	    #ratehist_p_veto2 = getRateVsEta(L1MuTriggerRate, "ratehist_position_mediumveto_eta%dto%d_npar%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), npar, fraction, pt), myetabin, ratecut + "&& abs(%s)<=%f"%(rate_x, ddY_cut), 2)
	    #ratehist_h_veto2 = getRateVsEta(L1MuTriggerRate, "ratehist_hybrid_mediumveto_eta%dto%d_npar%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), npar, fraction, pt), myetabin, ratecut + "&& "+ellipsecut_rate, 2)
	    #ratehist_h_noGE21_veto2 = getRateVsEta(L1MuTriggerRate, "ratehist_hybrid_noGE21_mediumveto_eta%dto%d_npar%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), npar, fraction, pt), myetabin, ratecut + "&& "+ellipsecut_rate_noGE21, 2)
	    hrateallnpars.append([ratehist_p_veto0, ratehist_h1_veto0, ratehist_h2_veto0, ratehist_h3_veto0])
	    #hrateallnpars_veto2.append([ratehist_p_veto2, ratehist_h_noGE21_veto2, ratehist_h_veto2])
	    ratehist_p_veto0.Write()
	    ratehist_h1_veto0.Write()
	    ratehist_h2_veto0.Write()
	    ratehist_h3_veto0.Write()
	    #ratehist_p_veto2.Write()
	    #ratehist_h_noGE21_veto2.Write()
	    #ratehist_h_veto2.Write()

	    numcuts = ["fabs(%s)<%f"%(xaxis, ddY_cut)+ L1Mucut, "fabs(%s)<%f"%(yaxis1, dphi_cut1)+"&& fabs(%s)<%f"%(xaxis, ddY_cut2)+ L1Mucut, "fabs(%s)<%f"%(yaxis2, dphi_cut2)+" && fabs(%s)<%f"%(xaxis, ddY_cut2)+ L1Mucut, "fabs(%s)<%f"%(yaxis3, dphi_cut3)+"&& fabs(%s)<%f"%(xaxis, ddY_cut2)+ L1Mucut]
	    Teffs = makeEffplot_pt([GEMCSCTrackCh0GE11, GEMCSCTrackCh0GE11, GEMCSCTrackCh0GE11, GEMCSCTrackCh0], "pt", [cutsbase+checkx, cutsbase+checkx+checky, cutsbase+checkx+checky, cutsbase+checkx+checky], numcuts, etamin, etamax,"true muon p_{T} GeV","Trigger efficiency",legsall,"no track veto",text,SBPlots+"GEMCSC_dxy5_25_combined_eff_20170611_pt%d_fraction%d_st2eta%dto%d_npar%d"%(pt, fraction, int(etamin*10),int(etamax*10), npar))
	    #print "Teffs dencut ",cutsbase+checkx," numcut ","fabs(%s)<%f"%(xaxis, ddY_cut)+ L1Mucut
	    Teffseta = makeEffplot_eta([GEMCSCTrackCh0GE11, GEMCSCTrackCh0GE11, GEMCSCTrackCh0GE11, GEMCSCTrackCh0], "abs(eta_st2_sim)", [cutsbase+checkx+"&& pt>%d"%(pt+5), cutsbase+checkx+checky+"&& pt>%d"%(pt+5), cutsbase+checkx+checky+"&& pt>%d"%(pt+5), cutsbase+checkx+checky+"&& pt>%d"%(pt+5)], numcuts, etamin, etamax,"|#eta| at 2nd CSC station","Trigger efficiency",legsall,"no track veto",text,SBPlots+"GEMCSC_dxy5_25_combined_effvseta_20170611_pt%d_fraction%d_st2eta%dto%d_npar%d"%(pt, fraction, int(etamin*10),int(etamax*10), npar))
	    #Teffs_veto2 = makeEffplot_pt([GEMCSCTrackCh0, GEMCSCTrackCh0, GEMCSCTrackCh0], "pt", [cutsbase+checkx, cutsbase+checkx+checky_noGE21, cutsbase+checkx+checky_noGE21], ["fabs(%s)<%f"%(xaxis, ddY_cut)+"&& !isL1MediumVeto", ellipse_noGE21+"&& !isL1MediumVeto", ellipse+"&& !isL1MediumVeto"], etamin, etamax,"true muon p_{T} GeV","Trigger efficiency",legs3,"medium track veto",text,SBPlots+"GEMCSC_dxy10_50_mediumveto_combined_eff_20170131_pt%d_fraction%d_st2eta%dto%d_npar%d"%(pt, fraction, int(etamin*10),int(etamax*10), npar))


	    #peff2, perr2 = getOverallEff(Teffs_ddY[1], pt)
	    peff, perr = getOverallEff(Teffs[0], pt+5)
	    h1eff, h1err = getOverallEff(Teffs[1], pt+5)
	    h2eff, h2err = getOverallEff(Teffs[2], pt+5)
	    h3eff, h3err = getOverallEff(Teffs[3], pt+5)
	    #heff_noGE21, herr_noGE21 = getOverallEff(Teffs[1], pt)
	    #peff_v2, perr_v2 = getOverallEff(Teffs_veto2[0], pt)
	    #heff_v2, herr_v2 = getOverallEff(Teffs_veto2[2], pt)
	    #heff_noGE21_v2, herr_noGE21_v2 = getOverallEff(Teffs_veto2[1], pt)
	    print "pt ",pt," fraction ",fraction," etamin ",etamin, " etamax ",etamax," npar ",npar, "peff ",peff," h1eff ",h1eff ," h2eff ",h2eff," h3eff ",h3eff," p_rate ", ratehist_p_veto0.GetEntries() ," h1rate ",ratehist_h1_veto0.GetEntries()," h1rate ",ratehist_h2_veto0.GetEntries()," h3rate ",ratehist_h3_veto0.GetEntries()
	    #makeplotsPt([Teffs_ddY[0], Teffs_dphi_dir[0], Teffs[0]], legs3, text,SBPlots+"GEMCSC_ctau0andctau1000_combined3_eff_20170131_pt%d_ptbg%d_fraction%d_fractiond%d_St2eta%dto%d_npar%d_dirinterstation%s"%(pt, pt1, fraction, fraction_d, int(etamin*10),int(etamax*10), npar, extrapic)) 
	    Teffallnpars.append(Teffs)
	    Teffallnparseta.append(Teffseta)
	    #Teffallnpars_veto2.append(Teffs_veto2)
	    #Teffallnpars_d.append(Teffs)
	    #Teffallnpars_d.append(Teffs_dphi_dir)
	    Teffs[0].SetName("positioneta%dto%dnpar%d"%(int(etamin*10),int(etamax*10), npar)+"displacedmuonpt%d"%(pt)+"fraction%d"%fraction)
	    Teffs[1].SetName("hybrid1eta%dto%dnpar%d"%(int(etamin*10),int(etamax*10), npar)+"displacedmuonpt%d"%(pt)+"fraction%d"%fraction)
	    Teffs[2].SetName("hybrid2eta%dto%dnpar%d"%(int(etamin*10),int(etamax*10), npar)+"displacedmuonpt%d"%(pt)+"fraction%d"%fraction)
	    Teffs[3].SetName("hybrid3eta%dto%dnpar%d"%(int(etamin*10),int(etamax*10), npar)+"displacedmuonpt%d"%(pt)+"fraction%d"%fraction)
	    Teffseta[0].SetName("positioneta%dto%dnpar%d"%(int(etamin*10),int(etamax*10), npar)+"displacedmuonpt%d"%(pt)+"fraction%dvseta"%fraction)
	    Teffseta[1].SetName("hybrid1eta%dto%dnpar%d"%(int(etamin*10),int(etamax*10), npar)+"displacedmuonpt%d"%(pt)+"fraction%dvseta"%fraction)
	    Teffseta[2].SetName("hybrid2eta%dto%dnpar%d"%(int(etamin*10),int(etamax*10), npar)+"displacedmuonpt%d"%(pt)+"fraction%dvseta"%fraction)
	    Teffseta[3].SetName("hybrid3eta%dto%dnpar%d"%(int(etamin*10),int(etamax*10), npar)+"displacedmuonpt%d"%(pt)+"fraction%dvseta"%fraction)
	    #Teffs_veto2[0].SetName("positionmediumvetoeta%dto%dnpar%d"%(int(etamin*10),int(etamax*10), npar)+"displacedmuonpt%d"%(pt)+"fraction%d"%fraction)
	    Teffs[0].Write();Teffs[1].Write();Teffs[2].Write(); Teffs[3].Write()
	    Teffseta[0].Write();Teffseta[1].Write();Teffseta[2].Write(); Teffseta[3].Write()
	    #Teffs_veto2[0].Write();Teffs_veto2[1].Write();Teffs_veto2[2].Write()
	    #Teffs_ddY[0].Write()
	    #Teffs_dphi_dir[0].Write()
	    #Teffs_ddY[1].Write()
	    #Teffs_dphi_dir[1].Write()
       #print "Teffallnpars len ",len(Teffallnpars),Teffallnpars
       Teff0_p = Teffallnpars[0][0]; Teff0_h1 = Teffallnpars[0][1];Teff0_h2 = Teffallnpars[0][2]; Teff0_h3 = Teffallnpars[0][3]
       Teffeta0_p = Teffallnparseta[0][0]; Teffeta0_h1 = Teffallnparseta[0][1];Teffeta0_h2 = Teffallnparseta[0][2]; Teffeta0_h3 = Teffallnparseta[0][3]
       #Teff2_p = Teffallnpars_veto2[0][0]; Teff2_h = Teffallnpars_veto2[0][2];Teff2_h_noGE21 = Teffallnpars_veto2[0][1]
       hrate0_p = hrateallnpars[0][0]; hrate0_h1 = hrateallnpars[0][1]; hrate0_h2 = hrateallnpars[0][2]; hrate0_h3 = hrateallnpars[0][3]
       #hrate2_p = hrateallnpars_veto2[0][0]; hrate2_h = hrateallnpars_veto2[0][2]; hrate2_h_noGE21 = hrateallnpars_veto2[0][1]
       for xpar in range(len(Teffallnpars)-1):
	    Teff0_p.Add(Teffallnpars[xpar+1][0])		   
	    Teff0_h1.Add(Teffallnpars[xpar+1][1])		   
	    Teff0_h2.Add(Teffallnpars[xpar+1][2])		   
	    Teff0_h3.Add(Teffallnpars[xpar+1][3])		   
	    Teffeta0_p.Add(Teffallnparseta[xpar+1][0])		   
	    Teffeta0_h1.Add(Teffallnparseta[xpar+1][1])		   
	    Teffeta0_h2.Add(Teffallnparseta[xpar+1][2])		   
	    Teffeta0_h3.Add(Teffallnparseta[xpar+1][3])		   
	    #Teff2_p.Add(Teffallnpars_veto2[xpar+1][0])		   
	    #Teff2_h.Add(Teffallnpars_veto2[xpar+1][2])		   
	    #Teff2_h_noGE21.Add(Teffallnpars_veto2[xpar+1][1])		   
	    hrate0_p.Add(hrateallnpars[xpar+1][0])
	    hrate0_h1.Add(hrateallnpars[xpar+1][1])
	    hrate0_h2.Add(hrateallnpars[xpar+1][2])
	    hrate0_h3.Add(hrateallnpars[xpar+1][3])
	    #hrate2_p.Add(hrateallnpars_veto2[xpar+1][0])
	    #hrate2_h.Add(hrateallnpars_veto2[xpar+1][2])
	    #hrate2_h_noGE21.Add(hrateallnpars_veto2[xpar+1][1])
       Teffall.append([Teff0_p, Teff0_h1, Teff0_h2, Teff0_h3])
       Teffalleta.append([Teffeta0_p, Teffeta0_h1, Teffeta0_h2, Teffeta0_h3])
       #Teffall2.append([Teff2_p, Teff2_h_noGE21, Teff2_h])
       hrateall.append([hrate0_p, hrate0_h1, hrate0_h2, hrate0_h3])
       #hrateall2.append([hrate2_p, hrate2_h_noGE21, hrate2_h])
       peff0, perr0 = getOverallEff(Teff0_p, pt)
       #peff2, perr2 = getOverallEff(Teff2_p, pt)
       h1eff0, h1err0 = getOverallEff(Teff0_h1, pt)
       #heff2_noGE21, herr2_noGE21 = getOverallEff(Teff2_h_noGE21, pt)
       h2eff0, h2err0 = getOverallEff(Teff0_h2, pt)
       h3eff0, h3err0 = getOverallEff(Teff0_h3, pt)
       #heff2, herr2 = getOverallEff(Teff2_h, pt)
       #print "pt ",pt," fraction ",fraction," etamin ",etamin, " etamax ",etamax," allnpar peff0 ",peff0, " heff0 ",heff0," heff0_noGE21 ",heff0_noGE21
       #print "Rate entreis p ",hrate0_p.GetEntries()," h_noGE21 ",hrate0_h_noGE21.GetEntries()," h ",hrate0_h.GetEntries()
       #print "all eff, medium veto p ",peff2," heff2_noGE21 ",heff2_noGE21," heff2 ",heff2
       #print "Rate, medium veto,  entreis p ",hrate2_p.GetEntries()," h_noGE21 ",hrate2_h_noGE21.GetEntries()," h ",hrate2_h.GetEntries()
       Teff0_p.SetName("positioneta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%d"%(fraction, pt))
       Teff0_h1.SetName("hybrid1eta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%d"%(fraction, pt))
       Teff0_h2.SetName("hybrid2eta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%d"%(fraction, pt))
       Teff0_h3.SetName("hybrid3eta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%d"%(fraction, pt))
       Teffeta0_p.SetName("positioneta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%dvseta"%(fraction, pt))
       Teffeta0_h1.SetName("hybrid1eta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%dvseta"%(fraction, pt))
       Teffeta0_h2.SetName("hybrid2eta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%dvseta"%(fraction, pt))
       Teffeta0_h3.SetName("hybrid3eta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%dvseta"%(fraction, pt))
       #Teff2_p.SetName("positionmediumvetoeta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonpt%d"%(pt))
       #Teff2_h.SetName("hybridmediumvetoeta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonpt%d"%(pt))
       #Teff2_h_noGE21.SetName("hybrid_noGE21mediumvetoeta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonpt%d"%(pt))
       Teff0_p.Write(); Teff0_h1.Write(); Teff0_h2.Write(); Teff0_h3.Write()
       Teffeta0_p.Write(); Teffeta0_h1.Write(); Teffeta0_h2.Write(); Teffeta0_h3.Write()
       #Teff2_p.Write(); Teff2_h.Write(); Teff2_h_noGE21.Write()
       #hrate2_p.SetName("ratehist_position_mediumveto_eta%dto%d_npar%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), npar, fraction, pt))
       hrate0_p.SetName("ratehist_position_eta%dto%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), fraction, pt))
       hrate0_h1.SetName("ratehist_hybrid1_eta%dto%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), fraction, pt))
       hrate0_h2.SetName("ratehist_hybrid2_eta%dto%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), fraction, pt))
       hrate0_h3.SetName("ratehist_hybrid3_eta%dto%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), fraction, pt))
       #hrate2_h_noGE21.SetName("ratehist_hybrid_noGE21_mediumveto_eta%dto%d_npar%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), npar, fraction, pt))
       hrate0_p.Write(); hrate0_h1.Write(); hrate0_h2.Write(); hrate0_h3.Write();
       #hrate2_p.Write(); hrate2_h_noGE21.Write(); hrate2_h.Write()
       text_h = "#splitline{Hybrid algorithm w/ GE21}{%.2f<|#eta|<%.2f, p_{T}^{L1}>%d GeV}"%(etamin,etamax, pt)
       text_p = "#splitline{Position-based algorithm}{%.2f<|#eta|<%.2f, p_{T}^{L1}>%d GeV}"%(etamin,etamax, pt)
       text_h_noGE21 = "#splitline{Hybrid algorithm w/o GE21}{%.2f<|#eta|<%.2f, p_{T}^{L1}>%d GeV}"%(etamin,etamax, pt)
       legsveto = ["no track veto","medium track veto"]
       #makeplotsPt([Teff0_p, Teff2_p], legsveto, text_p,SBPlots+"GEMCSC_dxy10_50_position_vetos_eff_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10))) 
       #makeplotsPt([Teff0_h, Teff2_h], legsveto, text_h,SBPlots+"GEMCSC_dxy10_50_hybrid_vetos_eff_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10))) 
       #makeplotsPt([Teff0_h_noGE21, Teff2_h_noGE21], legsveto, text_h_noGE21,SBPlots+"GEMCSC_dxy10_50_hybrid_noGE21_vetos_eff_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10))) 
       #plotshistograms([hrate0_p, hrate2_p], legsveto, text_p, SBPlots+"TriggerRate_position_vetos_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10)))
       #plotshistograms([hrate0_h, hrate2_h], legsveto, text_h, SBPlots+"TriggerRate_hybrid_vetos_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10)))
       #plotshistograms([hrate0_h_noGE21, hrate2_h_noGE21], legsveto, text_h_noGE21, SBPlots+"TriggerRate_hybrid_noGE21_vetos_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10)))
       text_noveto =  "%.2f<|#eta|<%.2f, p_{T}^{L1}>%d GeV, no track veto"%(etamin,etamax, pt)
       text_mveto =  "%.2f<|#eta|<%.2f, p_{T}^{L1}>%d GeV, medium track veto"%(etamin,etamax, pt)
       text_all =  "#splitline{%.2f<|#eta|<%.2f, p_{T}^{L1}>%d GeV}{5<|dxy|<25cm, Lxy<50cm}"%(etamin,etamax, pt)
       legheader = "L1Mu with Q>=%d, no track veto"%(quality)
       makeplotsPt([Teff0_p, Teff0_h1, Teff0_h2, Teff0_h3], legsall, legheader, text_all, SBPlots+"GEMCSC_dxy5_25_combined3_noveto_eff_20170611_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10))) 
       makeplotsEta([Teffeta0_p, Teffeta0_h1, Teffeta0_h2, Teffeta0_h3], legsall, legheader, text_all, SBPlots+"GEMCSC_dxy5_25_combined3_noveto_effvseta_20170611_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10))) 
       #makeplotsPt([Teff2_p, Teff2_h_noGE21, Teff2_h], legs3, text_mveto,SBPlots+"GEMCSC_dxy10_50_combined3_mediumveto_eff_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10))) 
       #plotshistograms([hrate0_p, hrate0_h_noGE21, hrate0_h], legs3, text_noveto, SBPlots+"TriggerRate_combiend_noveto_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10)))
       #plotshistograms([hrate2_p, hrate2_h_noGE21, hrate2_h], legs3, text_mveto, SBPlots+"TriggerRate_combiend_mediumveto_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10)))
    etamin = netas[0]; etamax = netas[-1]
    Teffall0_p = Teffall[0][0]; Teffall0_h2 = Teffall[0][2];Teffall0_h1 = Teffall[0][1]; Teffall0_h3 = Teffall[0][3]
    Teffalleta0_p = Teffalleta[0][0]; Teffalleta0_h2 = Teffalleta[0][2];Teffalleta0_h1 = Teffalleta[0][1]; Teffalleta0_h3 = Teffalleta[0][3]
    #Teffall2_p = Teffall2[0][0]; Teffall2_h = Teffall2[0][2];Teffall2_h_noGE21 = Teffall2[0][1]
    hrateall0_p = hrateall[0][0]; hrateall0_h2 = hrateall[0][2]; hrateall0_h1 = hrateall[0][1]; hrateall0_h3 = hrateall[0][3]
    #hrateall2_p = hrateall2[0][0]; hrateall2_h = hrateall2[0][2]; hrateall2_h_noGE21 = hrateall2[0][1]
    for xeta in range(len(Teffall)-1):
	Teffall0_p.Add(Teffall[xeta+1][0])		   
	Teffall0_h1.Add(Teffall[xeta+1][1])		   
	Teffall0_h2.Add(Teffall[xeta+1][2])		   
	Teffall0_h2.Add(Teffall[xeta+1][3])		   
	Teffalleta0_p.Add(Teffalleta[xeta+1][0])		   
	Teffalleta0_h1.Add(Teffalleta[xeta+1][1])		   
	Teffalleta0_h2.Add(Teffalleta[xeta+1][2])		   
	Teffalleta0_h3.Add(Teffalleta[xeta+1][3])		   
	#Teffall2_p.Add(Teffall2[xeta+1][0])		   
	#Teffall2_h_noGE21.Add(Teffall2[xeta+1][1])		   
	#Teffall2_h.Add(Teffall2[xeta+1][2])		   
        hrateall0_p.Add(hrateall[xeta+1][0])
        hrateall0_h1.Add(hrateall[xeta+1][1])
        hrateall0_h2.Add(hrateall[xeta+1][2])
        hrateall0_h3.Add(hrateall[xeta+1][3])
	#hrateall2_p.Add(hrateall2[xeta+1][0])
	#hrateall2_h_noGE21.Add(hrateall2[xeta+1][1])
	#hrateall2_h.Add(hrateall2[xeta+1][2])
    peff0, perr0 = getOverallEff(Teffall0_p, pt+5)
    #peff2, perr2 = getOverallEff(Teffall2_p, pt+5)
    h1eff0, h1err0 = getOverallEff(Teffall0_h1, pt+5)
    #heff2_noGE21, herr2_noGE21 = getOverallEff(Teffall2_h_noGE21, pt+5)
    h2eff0, h2err0 = getOverallEff(Teffall0_h2, pt+5)
    h3eff0, h3err0 = getOverallEff(Teffall0_h3, pt+5)
    #heff2, herr2 = getOverallEff(Teffall2_h, pt+5)
     
    print "pt ",pt," fraction ",fraction," etamin ",etamin, " etamax ",etamax," allnpar peff0 ",peff0, " h1eff0 ",h1eff0," h2eff0 ",h2eff0, " h3eff0 ",h3eff0
    print "Rate entreis p ",hrateall0_p.GetEntries()," h1 ",hrateall0_h1.GetEntries()," h2 ",hrateall0_h2.GetEntries(), " h3 ",hrateall0_h3.GetEntries()
    #print "all eff, medium veto p ",peff2," heff2_noGE21 ",heff2_noGE21," heff2 ",heff2
    #print "Rate, medium veto,  entreis p ",hrateall2_p.GetEntries()," h_noGE21 ",hrateall2_h_noGE21.GetEntries()," h ",hrateall2_h.GetEntries()
    Teffall0_p.SetName("positioneta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%d"%(fraction,pt))
    Teffall0_h1.SetName("hybrid1eta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%d"%(fraction,pt))
    Teffall0_h2.SetName("hybrid2eta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%d"%(fraction,pt))
    Teffall0_h3.SetName("hybrid3eta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%d"%(fraction,pt))
    Teffall0_p.Write(); Teffall0_h1.Write(); Teffall0_h2.Write(); Teffall0_h3.Write()
    Teffalleta0_p.SetName("positioneta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%dvseta"%(fraction, pt))
    Teffalleta0_h1.SetName("hybrid1eta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%dvseta"%(fraction, pt))
    Teffalleta0_h2.SetName("hybrid2eta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%dvseta"%(fraction,pt))
    Teffalleta0_h3.SetName("hybrid3eta%dto%d"%(int(etamin*10),int(etamax*10))+"displacedmuonfrac%dpt%dvseta"%(fraction, pt))
    Teffalleta0_p.Write(); Teffalleta0_h1.Write(); Teffalleta0_h2.Write(); Teffalleta0_h3.Write()
    #Teffall2_p.Write(); Teffall2_h.Write(); Teffall2_h_noGE21.Write()
    #hrateall2_p.SetName("ratehist_position_mediumveto_eta%dto%d_npar%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), npar, fraction, pt))
    #hrateall2_h_noGE21.SetName("ratehist_hybrid_noGE21_mediumveto_eta%dto%d_npar%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), npar, fraction, pt))
    hrateall0_p.SetName("ratehist_position_eta%dto%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), fraction, pt))
    hrateall0_h1.SetName("ratehist_hybrid1_eta%dto%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), fraction, pt))
    hrateall0_h2.SetName("ratehist_hybrid2_eta%dto%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), fraction, pt))
    hrateall0_h3.SetName("ratehist_hybrid3_eta%dto%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), fraction, pt))
    #hrateall2_h.SetName("ratehist_hybrid_mediumveto_eta%dto%d_npar%d_frac%d_pt%d"%(int(etamin*10),int(etamax*10), npar, fraction, pt))
    hrateall0_p.Write(); hrateall0_h1.Write(); hrateall0_h2.Write();  hrateall0_h3.Write()
    #hrateall2_p.Write(); hrateall2_h_noGE21.Write(); hrateall2_h.Write()
    text_h = "#splitline{Hybrid algorithm w/ GE21}{%.2f<|#eta|<%.2f, p_{T}^{L1}>%d GeV}"%(etamin,etamax, pt)
    text_p = "#splitline{Position-based algorithm}{%.2f<|#eta|<%.2f, p_{T}^{L1}>%d GeV}"%(etamin,etamax, pt)
    text_h_noGE21 = "#splitline{Hybrid algorithm w/o GE21}{%.2f<|#eta|<%.2f, p_{T}^{L1}>%d GeV}"%(etamin,etamax, pt)
    legsveto = ["no track veto","medium track veto"]
    legheader = "L1Mu with Q>=12"
    #makeplotsPt([Teffall0_p, Teffall2_p], legsveto, legheader, text_p,SBPlots+"GEMCSC_dxy10_50_position_vetos_eff_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10))) 
    #makeplotsPt([Teffall0_h, Teffall2_h], legsveto, legheader,text_h,SBPlots+"GEMCSC_dxy10_50_hybrid_vetos_eff_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10))) 
    #makeplotsPt([Teffall0_h_noGE21, Teffall2_h_noGE21], legsveto, legheader,text_h_noGE21,SBPlots+"GEMCSC_dxy10_50_hybrid_noGE21_vetos_eff_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10))) 
    #plotshistograms([hrateall0_p, hrateall2_p], legsveto,"L1Mu with Q>=4", text_p, SBPlots+"TriggerRate_position_vetos_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10)))
    #plotshistograms([hrateall0_h, hrateall2_h], legsveto, "L1Mu with Q>=4", text_h, SBPlots+"TriggerRate_hybrid_vetos_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10)))
    #plotshistograms([hrateall0_h_noGE21, hrateall2_h_noGE21], legsveto, "L1Mu with Q>=4", text_h_noGE21, SBPlots+"TriggerRate_hybrid_noGE21_vetos_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10)))
    text_noveto =  "%.2f<|#eta|<%.2f, p_{T}^{L1}>%d GeV, no track veto"%(etamin,etamax, pt)
    text_mveto =  "%.2f<|#eta|<%.2f, p_{T}^{L1}>%d GeV, medium track veto"%(etamin,etamax, pt)
    #text_all =  "#splitline{%.2f<|#eta|<%.2f, p_{T}^{L1}>%d GeV}{ 5<|dxy|<25cm, Lxy<50cm}"%(etamin,etamax, pt)
    text_all =  "%.2f<|#eta|<%.2f, p_{T}^{L1}>%d GeV"%(etamin,etamax, pt)
    #text_rate =  "%.2f<|#eta|<%.2f, p_{T}^{L1}>%d GeV"%(etamin,etamax, pt)
    text_rate =  "Displaced Muon Trigger, %.2f<|#eta|<%.2f"%(etamin,etamax)
    legheader = "L1Mu, no track veto"
    legheaderm = "L1Mu Q>=%d, medium track veto"%(quality)
    legs2all = ["CSC+GE11: position based","CSC+GE11+GE21+ME0: hybrid algorithm"]
    makeplotsPt([Teffall0_h3, Teffall0_p], [legsall[1], legs2all[0]], legheader,text_all,SBPlots+"GEMCSC_dxy5_25_combined3_noveto_eff_v2_20170611_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10))) 
    makeplotsPt([Teffall0_p, Teffall0_h1, Teffall0_h2 ,Teffall0_h3], legsall, legheader,text_all,SBPlots+"GEMCSC_dxy5_25_combined3_noveto_eff_20170611_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10))) 
    makeplotsEta([Teffalleta0_h3, Teffalleta0_p], [legsall[1], legs2all[0]], legheader,text_all,SBPlots+"GEMCSC_dxy5_25_combined3_noveto_effvseta_v2_20170611_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10))) 
    makeplotsEta([Teffalleta0_p, Teffalleta0_h1, Teffalleta0_h2 ,Teffalleta0_h3], legsall, legheader,text_all,SBPlots+"GEMCSC_dxy5_25_combined3_noveto_effvseta_20170611_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10))) 
    #makeplotsPt([Teffall2_p, Teffall2_h_noGE21, Teffall2_h], legs3, legheaderm,text_mveto,SBPlots+"GEMCSC_dxy10_50_combined3_mediumveto_eff_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10))) 
    plotshistograms([hrateall0_p, hrateall0_h3], legs2all, legheader, text_rate, SBPlots+"TriggerRatevseta_v2_combiend_noveto_20170611_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10)))
    plotshistograms([hrateall0_p, hrateall0_h1, hrateall0_h2, hrateall0_h3], legsall, legheader, text_rate, SBPlots+"TriggerRate_combiend_noveto_20170611_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10)))
    legssteps = ["L1Mu, Q>=4, no pt cut","L1Mu,Q>=4 + hits in ME11,ME21,ME31","Position based, Q>=4","hybrid with ME0, Q>=4"]
    legsstepsv2 = ["L1Mu, Q>=4, no pt cut","L1Mu,Q>=4 + hits in ME11,ME21,ME31","L1Mu, Q>=4, hits in ME11, ME21, ME21, ME0"]
    legsstepsv3 = ["No pt cut","No pt cut, hits in ME11,ME21,ME31","No pt cut, hits in ME11, ME21, ME21, ME0","Position based, using CSC only, p_{T}^{Displaced}>%d GeV"%pt,"Hybrid using CSC+GE11+GE21+ME0,p_{T}^{Displaced}>%d GeV"%pt]
    plotshistograms([ratehist_step0, ratehist_step1, hrateall0_p, hrateall0_h3], legssteps, legheader, text_rate, SBPlots+"TriggerRate_steps_noveto_20170611_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10)))
    plotshistograms([ratehist_step0, ratehist_step1, ratehist_step2], legsstepsv2, legheader, text_rate, SBPlots+"TriggerRate_stepsv2_noveto_20170611_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10)))
    plotshistograms([ratehist_step0, ratehist_step1, ratehist_step2 ,hrateall0_p, hrateall0_h3], legsstepsv3, legheader, text_rate, SBPlots+"TriggerRate_stepsv3_noveto_20170611_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10)))
    #plotshistograms([hrateall2_p, hrateall2_h_noGE21, hrateall2_h], legs3, legheaderm, text_mveto, SBPlots+"TriggerRate_combiend_mediumveto_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10)))
    return peff0, perr0, h1eff0, h1err0, h2eff0, h2err0, h3eff0, h3err0, hrateall0_p.Integral(), hrateall0_h1.Integral(), hrateall0_h2.Integral(), hrateall0_h3.Integral()
    Teff_out.Close()


def plotalletaallptallfractionsME0GE21(pts, netas, fractions, fraction_rate, outputrootfile, ratesample, signalsample, allLUTfile, SBPlots):
    hrate_template = ROOT.TH1F("hrate_ME0_fractionrate%d"%fraction_rate,"hrate_ME0",len(rateptbin)-1, myrateptbin)
    hrate_p = ROOT.TH1F("hrate_p_fractionrate%d"%fraction_rate,"hrate_p_fractionrate%d"%fraction_rate,len(rateptbin)-1, myrateptbin)
    hrate_h1 = ROOT.TH1F("hrate_h1_fractionrate%d"%fraction_rate,"hrate_h1_fractionrate%d"%fraction_rate,len(rateptbin)-1, myrateptbin)
    hrate_h2 = ROOT.TH1F("hrate_h2_fractionrate%d"%fraction_rate,"hrate_h2_fractionrate%d"%fraction_rate,len(rateptbin)-1, myrateptbin)
    hrate_h3 = ROOT.TH1F("hrate_h3_fractionrate%d"%fraction_rate,"hrate_h3_fractionrate%d"%fraction_rate,len(rateptbin)-1, myrateptbin)
    #hrate2_p = ROOT.TH1F("hrate2_p_fractionrate%d"%fraction_rate,"hrate2_p_fractionrate%d"%fraction_rate,len(rateptbin)-1, myrateptbin)
    #hrate2_h_noGE21 = ROOT.TH1F("hrate2_h_noGE21_fractionrate%d"%fraction_rate,"hrate2_h_noGE21_fractionrate%d"%fraction_rate,len(rateptbin)-1, myrateptbin)
    #hrate2_h = ROOT.TH1F("hrate2_h_fractionrate%d"%fraction_rate,"hrate2_h_fractionrate%d"%fraction_rate,len(rateptbin)-1, myrateptbin)
    etamin = netas[0]; etamax = netas[-1]
    #legs3 = ["Position-based algorithm","Hybrid algorithm w/o GE21","Hybrid algorithm w GE21"]
    legsall = ["Position-based algorithm","Hybrid algirthm+ w/o ME11+ w/o GE21","Hybrid algorithm+ w/o ME0 + w GE21"," Hybrid algorithm+ w ME0 + w GE21"]
    legsall = ["Position-based algorithm","Stub alignment algorithm+ w/o ME0 + w/o GE21","Stub alignment algorithm+ w/o ME0 + w GE21","Stub alignment algorithm+ w ME0 + w GE21"]
    for pt in pts:
	
	ratelist_p = []
	rateErrlist_p = []
	ratelist_h1 = []
	rateErrlist_h1 = []
	ratelist_h2 = []
	rateErrlist_h2 = []
	ratelist_h3 = []
	rateErrlist_h3 = []
	efflist_p = [];effErrlist_p = [];efflist_h1 = [];effErrlist_h1 = [];efflist_h2 = [];effErrlist_h2 = [];efflist_h3 = [];effErrlist_h3 = []
	#efflist2_p = [];effErrlist2_p = [];efflist2_h_noGE21 = [];effErrlist2_h_noGE21 = [];efflist2_h = [];effErrlist2_h = []
	for fraction in fractions:
	    #peff0, perr0, heff0_noGE21, herr0_noGE21, heff0, herr0, peff2, perr2, heff2_noGE21, herr2_noGE21, heff2, herr2, nrate0_p, nrate0_h_noGE21, nrate0_h ,nrate2_p, nrate2_h_noGE21, nrate2_h = plotalleta(pt, netas, fraction, outputrootfile, ratesample, signalsample, allLUTfile, SBPlots)
	    peff0, perr0, h1eff0, h1err0, h2eff0, h2err0, h3eff0, h3err0, nrate0_p, nrate0_h1, nrate0_h2 ,nrate0_h3 = plotalletaME0GE21(pt, netas, fraction, outputrootfile, ratesample, signalsample, allLUTfile, SBPlots)
	    if fraction==fraction_rate:
		ipt = hrate_p.FindBin(pt)
		hrate_p.SetBinContent(ipt, nrate0_p)
		hrate_h1.SetBinContent(ipt, nrate0_h1)
		hrate_h2.SetBinContent(ipt, nrate0_h2)
		hrate_h3.SetBinContent(ipt, nrate0_h3)
		hrate_p.SetBinError(ipt, sqrt(nrate0_p*SF))
		hrate_h1.SetBinError(ipt, sqrt(nrate0_h1*SF))
		hrate_h2.SetBinError(ipt, sqrt(nrate0_h2*SF))
		hrate_h3.SetBinError(ipt, sqrt(nrate0_h3*SF))
		#print "rate value prate ",nrate0_p," hrate ",nrate0_h3
		#hrate2_p.Fill(pt, nrate2_p)
		#hrate2_h_noGE21.Fill(pt, nrate2_h_noGE21)
		#hrate2_h.Fill(pt, nrate2_h)
	    efflist_p.append(peff0); effErrlist_p.append(perr0)
	    efflist_h1.append(h1eff0); effErrlist_h1.append(h1err0)
	    efflist_h2.append(h2eff0); effErrlist_h2.append(h2err0)
	    efflist_h3.append(h3eff0); effErrlist_h3.append(h3err0)
	    #efflist2_p.append(peff2); effErrlist2_p.append(perr2)
	    #efflist2_h_noGE21.append(heff2_noGE21); effErrlist2_h_noGE21.append(herr2_noGE21)
	    #efflist2_h.append(heff2); effErrlist2_h.append(herr2)
	    ratelist_p.append(nrate0_p); rateErrlist_p.append(sqrt(nrate0_p*SF))
	    ratelist_h1.append(nrate0_h1); rateErrlist_h1.append(sqrt(nrate0_h1*SF))
	    ratelist_h2.append(nrate0_h2); rateErrlist_h2.append(sqrt(nrate0_h2*SF))
	    ratelist_h3.append(nrate0_h3); rateErrlist_h3.append(sqrt(nrate0_h3*SF))
	    #ratelist2_p.append(nrate2_p); rateErrlist2_p.append(sqrt(nrate2_p*SF))
	    #ratelist2_h.append(nrate2_h); rateErrlist2_h.append(sqrt(nrate2_h*SF))
	    #ratelist2_h_noGE21.append(nrate2_h_noGE21); rateErrlist2_h_noGE21.append(sqrt(nrate2_h_noGE21*SF))
	Teff_out = ROOT.TFile(outputrootfile,"UPDATE")
	graph_p = TGraphErrors(len(efflist_p), np.array(efflist_p), np.array(ratelist_p), np.array(effErrlist_p), np.array(rateErrlist_p))
	graph_h1 = TGraphErrors(len(efflist_h1), np.array(efflist_h1), np.array(ratelist_h1), np.array(effErrlist_h1), np.array(rateErrlist_h1))
	graph_h2 = TGraphErrors(len(efflist_h2), np.array(efflist_h2), np.array(ratelist_h2), np.array(effErrlist_h2), np.array(rateErrlist_h2))
	graph_h3 = TGraphErrors(len(efflist_h3), np.array(efflist_h3), np.array(ratelist_h3), np.array(effErrlist_h3), np.array(rateErrlist_h3))
	#graph2_p = TGraphErrors(len(efflist2_p), np.array(efflist2_p), np.array(ratelist2_p), np.array(effErrlist2_p), np.array(rateErrlist2_p))
	#graph2_h = TGraphErrors(len(efflist2_h), np.array(efflist2_h), np.array(ratelist2_h), np.array(effErrlist2_h), np.array(rateErrlist2_h))
	#graph2_h_noGE21 = TGraphErrors(len(efflist2_h_noGE21), np.array(efflist2_h_noGE21), np.array(ratelist2_h_noGE21), np.array(effErrlist2_h_noGE21), np.array(rateErrlist2_h_noGE21))
	legend = TLegend(0.45,0.75,0.90,0.90)
	legend.SetFillColor(ROOT.kWhite)
	legend.SetTextFont(42)
	tex = TLatex(0.2,.88, "%.2f<|#eta|<%.2f, L1Mu with p_{T}^{L1}>%d GeV, Q>=12"%(etamin, etamax, pt))
        graph_p.SetName("TriggerEffVsRate_position_pt%d_eta%dto%d"%(pt,int(etamin*10), int(etamax*10)))
        graph_h1.SetName("TriggerEffVsRate_hybrid1_pt%d_eta%dto%d"%(pt,int(etamin*10), int(etamax*10)))
        graph_h2.SetName("TriggerEffVsRate_hybrid2_pt%d_eta%dto%d"%(pt,int(etamin*10), int(etamax*10)))
        graph_h3.SetName("TriggerEffVsRate_hybrid3_pt%d_eta%dto%d"%(pt,int(etamin*10), int(etamax*10)))
	#graph2_p.SetName("TriggerEffVsRate_position_mediumveto_pt%d_eta%dto%d"%(pt,int(etamin*10), int(etamax*10)))
	#graph2_h.SetName("TriggerEffVsRate_hybrid_mediumveto_pt%d_eta%dto%d"%(pt,int(etamin*10), int(etamax*10)))
	#graph2_h_noGE21.SetName("TriggerEffVsRate_hybrid_noGE21_mediumveto_pt%d_eta%dto%d"%(pt,int(etamin*10), int(etamax*10)))
	graph_p.Write();graph_h1.Write(); graph_h2.Write(); graph_h3.Write()
	#graph2_p.Write();graph2_h_noGE21.Write();graph2_h.Write()
	text_noveto =  "%.2f<|#eta|<%.2f, p_{T}^{L1}>%d GeV"%(etamin,etamax,pt)
	text_mveto =  "%.2f<|#eta|<%.2f, p_{T}^{L1}>%d GeV"%(etamin,etamax,pt)
	#ratelist2_p = [];rateErrlist2_p = [];ratelist2_h_noGE21 = [];rateErrlist2_h_noGE21 = [];ratelist2_h = [];rateErrlist2_h = []
	legheader = "L1Mu with Q>=12, no track veto"
	legheaderm = "L1Mu with medium track veto, Q>=12"
	mgerr = plotsgraphs([graph_p, graph_h1, graph_h2, graph_h3], legsall, legheader, text_noveto, SBPlots+"TriggerRateVsEff_combiend_noveto_20170611_pt%d_St2eta%dto%d_allnpar"%(pt, int(etamin*10),int(etamax*10)))
	#mgerr2 = plotsgraphs([graph2_p, graph2_h_noGE21, graph2_h], legs3, legheaderm, text_mveto, SBPlots+"TriggerRateVsEff_combiend_mediumveto_20170131_pt%d_St2eta%dto%d_allnpar"%(pt, int(etamin*10),int(etamax*10)))
        mgerr.SetName("TriggerEffVsRate_combined_noveto_pt%d_eta%dto%d"%(pt,int(etamin*10), int(etamax*10)))
	#mgerr2.SetName("TriggerEffVsRate_combined_mediumveto_pt%d_eta%dto%d"%(pt,int(etamin*10), int(etamax*10)))
	mgerr.Write()
	#mgerr2.Write()
	Teff_out.Close()
    #print "Rate events, p ",hrate_p.GetEntries()," h1 ", hrate_h1.GetEntries(), " h2 ",hrate_h2.GetEntries()," h3 ",hrate_h3.GetEntries()
    Teff_out = ROOT.TFile(outputrootfile,"UPDATE")
    #c2 = TCanvas()
    #c2.SetLogy()
    #c2.SetLogx()
    #c2.SetGridx()
    #c2.SetGridy()
    #legend2 = TLegend(0.4,0.7,0.90,0.85)
    #legend2.SetFillColor(ROOT.kWhite)
    #legend2.SetTextFont(42)
    #hrate_p.Sumw2(False)
    #hrate_h1.Sumw2(False)
    #hrate_h2.Sumw2(False)
    #hrate_h3.Sumw2(False)
    #hrate2_p.Sumw2(False)
    #hrate2_h.Sumw2(False)
    #hrate2_h_noGE21.Sumw2(False)
    hrate_template.SetMinimum(.2)
    hrate_template.SetMaximum(2000)
    #hrate_template.SetTitle(" #scale[1.4]{#font[61]{CMS}} #font[52]{Simulation preliminary}"+"  "*10+" 14 TeV, 200 PU")
    hrate_template.SetTitle(" #scale[1.4]{#font[61]{CMS}} Phase-2 Simulation"+"  "*4+" #sqrt{s}=14 TeV, <PU>=200")
    hrate_template.GetXaxis().SetTitle("L1Mu p_{T} [GeV]")
    hrate_template.GetYaxis().SetTitle("Trigger rate [kHz]")
    hrate_template.Write()
    hrate_p.Write();hrate_h1.Write();hrate_h2.Write(); hrate_h3.Write()
    #hrate2_p.Write();hrate2_h.Write();hrate2_h_noGE21.Write()
    text_noveto =  "%.2f<|#eta|<%.2f"%(etamin,etamax)
    text_mveto =  "%.2f<|#eta|<%.2f"%(etamin,etamax)
    legheader = "L1Mu, no track veto"
    legheaderm = "L1Mu, medium track veto"
    legs2all = ["CSC+GE11: position based","CSC+GE11+GE21+ME0: hybrid algorithm"]
    plotshistograms([hrate_p, hrate_h3], legs2all,legheader, text_noveto, SBPlots+"TriggerRate_pt_combined_noveto_v2_20170611_ratefraction%d_St2eta%dto%d_allnpar"%(fraction_rate, int(etamin*10),int(etamax*10)), True)
    #plotshistograms([hrate_p, hrate_h1, hrate_h2, hrate_h3], legsall,legheader, text_noveto, SBPlots+"TriggerRate_pt_combined_noveto_20170611_ratefraction%d_St2eta%dto%d_allnpar"%(fraction_rate, int(etamin*10),int(etamax*10)), True)

def plotalletaallptallfractions(pts, netas, fractions, fraction_rate, outputrootfile, ratesample, signalsample, allLUTfile, SBPlots):
    hrate_template = ROOT.TH1F("hrate_ge11_fractionrate%d"%fraction_rate,"hrate_ge11",len(rateptbin)-1, myrateptbin)
    hrate_p = ROOT.TH1F("hrate_p_fractionrate%d"%fraction_rate,"hrate_p_fractionrate%d"%fraction_rate,len(rateptbin)-1, myrateptbin)
    hrate_h_noGE21 = ROOT.TH1F("hrate_h_noGE21_fractionrate%d"%fraction_rate,"hrate_h_noGE21_fractionrate%d"%fraction_rate,len(rateptbin)-1, myrateptbin)
    hrate_h = ROOT.TH1F("hrate_h_fractionrate%d"%fraction_rate,"hrate_h_fractionrate%d"%fraction_rate,len(rateptbin)-1, myrateptbin)
    hrate2_p = ROOT.TH1F("hrate2_p_fractionrate%d"%fraction_rate,"hrate2_p_fractionrate%d"%fraction_rate,len(rateptbin)-1, myrateptbin)
    hrate2_h_noGE21 = ROOT.TH1F("hrate2_h_noGE21_fractionrate%d"%fraction_rate,"hrate2_h_noGE21_fractionrate%d"%fraction_rate,len(rateptbin)-1, myrateptbin)
    hrate2_h = ROOT.TH1F("hrate2_h_fractionrate%d"%fraction_rate,"hrate2_h_fractionrate%d"%fraction_rate,len(rateptbin)-1, myrateptbin)
    etamin = netas[0]; etamax = netas[-1]
    legs3 = ["Position-based algorithm","Hybrid algorithm w/o GE21","Hybrid algorithm w GE21"]
    for pt in pts:
	
	ratelist_p = []
	rateErrlist_p = []
	ratelist_h_noGE21 = []
	rateErrlist_h_noGE21 = []
	ratelist_h = []
	rateErrlist_h = []
	ratelist2_p = [];rateErrlist2_p = [];ratelist2_h_noGE21 = [];rateErrlist2_h_noGE21 = [];ratelist2_h = [];rateErrlist2_h = []
	efflist_p = [];effErrlist_p = [];efflist_h_noGE21 = [];effErrlist_h_noGE21 = [];efflist_h = [];effErrlist_h = []
	efflist2_p = [];effErrlist2_p = [];efflist2_h_noGE21 = [];effErrlist2_h_noGE21 = [];efflist2_h = [];effErrlist2_h = []
	for fraction in fractions:
	    peff0, perr0, heff0_noGE21, herr0_noGE21, heff0, herr0, peff2, perr2, heff2_noGE21, herr2_noGE21, heff2, herr2, nrate0_p, nrate0_h_noGE21, nrate0_h ,nrate2_p, nrate2_h_noGE21, nrate2_h = plotalleta(pt, netas, fraction, outputrootfile, ratesample, signalsample, allLUTfile, SBPlots)
	    if fraction==fraction_rate:
		ipt = hrate_p.FindBin(pt)
		hrate_p.SetBinContent(ipt, nrate0_p)
		hrate_h_noGE21.SetBinContent(ipt, nrate0_h_noGE21)
		hrate_h.SetBinContent(ipt, nrate0_h)
		hrate_p.SetBinError(ipt, sqrt(nrate0_p*SF))
		hrate_h_noGE21.SetBinError(ipt, sqrt(nrate0_p*SF))
		hrate_h.SetBinError(ipt, sqrt(nrate0_p*SF))

		hrate2_p.SetBinContent(ipt, nrate2_p)
		hrate2_h_noGE21.SetBinContent(ipt, nrate2_h_noGE21)
		hrate2_h.SetBinContent(ipt, nrate2_h)
		hrate2_p.SetBinError(ipt, sqrt(nrate2_p*SF))
		hrate2_h_noGE21.SetBinError(ipt, sqrt(nrate2_p*SF))
		hrate2_h.SetBinError(ipt, sqrt(nrate2_p*SF))
	    efflist_p.append(peff0); effErrlist_p.append(perr0)
	    efflist_h_noGE21.append(heff0_noGE21); effErrlist_h_noGE21.append(herr0_noGE21)
	    efflist_h.append(heff0); effErrlist_h.append(herr0)
	    efflist2_p.append(peff2); effErrlist2_p.append(perr2)
	    efflist2_h_noGE21.append(heff2_noGE21); effErrlist2_h_noGE21.append(herr2_noGE21)
	    efflist2_h.append(heff2); effErrlist2_h.append(herr2)
	    ratelist_p.append(nrate0_p); rateErrlist_p.append(sqrt(nrate0_p*SF))
	    ratelist_h.append(nrate0_h); rateErrlist_h.append(sqrt(nrate0_h*SF))
	    ratelist_h_noGE21.append(nrate0_h_noGE21); rateErrlist_h_noGE21.append(sqrt(nrate0_h_noGE21*SF))
	    ratelist2_p.append(nrate2_p); rateErrlist2_p.append(sqrt(nrate2_p*SF))
	    ratelist2_h.append(nrate2_h); rateErrlist2_h.append(sqrt(nrate2_h*SF))
	    ratelist2_h_noGE21.append(nrate2_h_noGE21); rateErrlist2_h_noGE21.append(sqrt(nrate2_h_noGE21*SF))
	Teff_out = ROOT.TFile(outputrootfile,"UPDATE")
	graph_p = TGraphErrors(len(efflist_p), np.array(efflist_p), np.array(ratelist_p), np.array(effErrlist_p), np.array(rateErrlist_p))
	graph_h = TGraphErrors(len(efflist_h), np.array(efflist_h), np.array(ratelist_h), np.array(effErrlist_h), np.array(rateErrlist_h))
	graph_h_noGE21 = TGraphErrors(len(efflist_h_noGE21), np.array(efflist_h_noGE21), np.array(ratelist_h_noGE21), np.array(effErrlist_h_noGE21), np.array(rateErrlist_h_noGE21))
	graph2_p = TGraphErrors(len(efflist2_p), np.array(efflist2_p), np.array(ratelist2_p), np.array(effErrlist2_p), np.array(rateErrlist2_p))
	graph2_h = TGraphErrors(len(efflist2_h), np.array(efflist2_h), np.array(ratelist2_h), np.array(effErrlist2_h), np.array(rateErrlist2_h))
	graph2_h_noGE21 = TGraphErrors(len(efflist2_h_noGE21), np.array(efflist2_h_noGE21), np.array(ratelist2_h_noGE21), np.array(effErrlist2_h_noGE21), np.array(rateErrlist2_h_noGE21))
	legend = TLegend(0.45,0.75,0.90,0.90)
	legend.SetFillColor(ROOT.kWhite)
	legend.SetTextFont(42)
	tex = TLatex(0.2,.88, "%.2f<|#eta|<%.2f, L1Mu with p_{T}^{L1}>%d GeV, Q>=4"%(etamin, etamax, pt))
        graph_p.SetName("TriggerEffVsRate_position_pt%d_eta%dto%d"%(pt,int(etamin*10), int(etamax*10)))
        graph_h.SetName("TriggerEffVsRate_hybrid_pt%d_eta%dto%d"%(pt,int(etamin*10), int(etamax*10)))
        graph_h_noGE21.SetName("TriggerEffVsRate_hybrid_noGE21_pt%d_eta%dto%d"%(pt,int(etamin*10), int(etamax*10)))
        graph2_p.SetName("TriggerEffVsRate_position_mediumveto_pt%d_eta%dto%d"%(pt,int(etamin*10), int(etamax*10)))
        graph2_h.SetName("TriggerEffVsRate_hybrid_mediumveto_pt%d_eta%dto%d"%(pt,int(etamin*10), int(etamax*10)))
        graph2_h_noGE21.SetName("TriggerEffVsRate_hybrid_noGE21_mediumveto_pt%d_eta%dto%d"%(pt,int(etamin*10), int(etamax*10)))
	graph_p.Write();graph_h_noGE21.Write();graph_h.Write()
	graph2_p.Write();graph2_h_noGE21.Write();graph2_h.Write()
	text_noveto =  "%.2f<|#eta|<%.2f, p_{T}^{L1}>%d GeV"%(etamin,etamax,pt)
	text_mveto =  "%.2f<|#eta|<%.2f, p_{T}^{L1}>%d GeV"%(etamin,etamax,pt)
	legheader = "L1Mu with no track veto, Q>=4"
	legheaderm = "L1Mu with medium track veto, Q>=4"
	mgerr = plotsgraphs([graph_p, graph_h_noGE21, graph_h], legs3, legheader, text_noveto, SBPlots+"TriggerRateVsEff_combiend_noveto_20170131_pt%d_St2eta%dto%d_allnpar"%(pt, int(etamin*10),int(etamax*10)))
	mgerr2 = plotsgraphs([graph2_p, graph2_h_noGE21, graph2_h], legs3, legheaderm, text_mveto, SBPlots+"TriggerRateVsEff_combiend_mediumveto_20170131_pt%d_St2eta%dto%d_allnpar"%(pt, int(etamin*10),int(etamax*10)))
        mgerr.SetName("TriggerEffVsRate_combined_noveto_pt%d_eta%dto%d"%(pt,int(etamin*10), int(etamax*10)))
        mgerr2.SetName("TriggerEffVsRate_combined_mediumveto_pt%d_eta%dto%d"%(pt,int(etamin*10), int(etamax*10)))
	mgerr.Write()
	mgerr2.Write()
	Teff_out.Close()
    Teff_out = ROOT.TFile(outputrootfile,"UPDATE")
    #c2 = TCanvas()
    #c2.SetLogy()
    #c2.SetLogx()
    #c2.SetGridx()
    #c2.SetGridy()
    #legend2 = TLegend(0.4,0.7,0.90,0.85)
    #legend2.SetFillColor(ROOT.kWhite)
    #legend2.SetTextFont(42)
    #hrate_p.Sumw2(False)
    #hrate_h.Sumw2(False)
    #hrate_h_noGE21.Sumw2(False)
    #hrate2_p.Sumw2(False)
    #hrate2_h.Sumw2(False)
    #hrate2_h_noGE21.Sumw2(False)
    hrate_template.SetMinimum(.2)
    hrate_template.SetMaximum(2000)
    #hrate_template.SetTitle(" #scale[1.4]{#font[61]{CMS}} #font[52]{Simulation preliminary}"+"  "*10+" 14 TeV, %d PU"%pu)
    hrate_template.SetTitle(" #scale[1.4]{#font[61]{CMS}} Phase-2 Simulation"+"  "*4+" #sqrt{s}=14 TeV, <PU>=%d"%pu)
    hrate_template.GetXaxis().SetTitle("L1Mu p_{T} [GeV]")
    hrate_template.GetYaxis().SetTitle("Trigger rate [kHz]")
    hrate_template.Write()
    hrate_p.Write();hrate_h.Write();hrate_h_noGE21.Write()
    hrate2_p.Write();hrate2_h.Write();hrate2_h_noGE21.Write()
    text_noveto =  "%.2f<|#eta|<%.2f"%(etamin,etamax)
    text_mveto =  "%.2f<|#eta|<%.2f"%(etamin,etamax)
    legheader = "L1Mu, no track veto"
    legheaderm = "L1Mu, medium track veto"
    plotshistograms([hrate_p, hrate_h], legs3,legheader, text_noveto, SBPlots+"TriggerRate_pt_combined_noveto_20170131_ratefraction%d_St2eta%dto%d_allnpar"%(fraction_rate, int(etamin*10),int(etamax*10)), True)
    plotshistograms([hrate2_p, hrate2_h], legs3,legheaderm, text_mveto, SBPlots+"TriggerRate_pt_combined_mediumveto_20170131_ratefraction%d_St2eta%dto%d_allnpar"%(fraction_rate, int(etamin*10),int(etamax*10)), True)

def getallTriggerRateAndEffVsEta(outputrootfile, fractions, pt, netas, SBPlots):
	
    quality = 12
    #plotalleta(pt, 0, netas, fractions[0], outputrootfile, badGE11)
    #plotalleta(pt, 0, netas, fractions[1], outputrootfile, badGE11)
    rfile1 = TFile(outputrootfile1,"READ")
    rfile2 = TFile(outputrootfile2,"READ")

    etamin = netas[0]; etamax = netas[-1];
    
    legsall = ["CSC + GE11","CSC + GE11 + GE21", "CSC + GE11 + GE21 + ME0"]
    ratehist_CSC = rfile.Get("ratehist_position_eta12to21_frac%d_pt%d"%(fractions[0], pt))
    ratehist_CSC2 = rfile.Get("ratehist_position_eta20to24_frac%d_pt%d"%(fractions[0], pt))
    ratehist_CSC2_v2= getPartOfHist(ratehist_CSC2, 2.16, 2.4)
    ratehist_CSC.Add(ratehist_CSC2_v2)
    ratehist_h_all = rfile.Get("ratehist_hybrid_eta12to21_frac%d_pt%d"%(fractions[0], pt))
    ratehist_st1 = rfile.Get("ratehist_hybrid_noGE21_eta12to21_frac%d_pt%d"%(fractions[0], pt))
    ratehist_h_all =  getPartOfHist(ratehist_h_all, 1.6, 2.15)
    ratehist_st1 =  getPartOfHist(ratehist_st1, 1.6, 2.15)
    #### ME0 region
    ratehist_ME0 = rfile.Get("ratehist_hybrid3_eta20to24_frac%d_pt%d"%(fractions[0], pt))
    ratehist_ME0.Add(getPartOfHist(ratehist_h_all, 1.6, 2.06))
    ratehist_ME0.Sumw2()

    eff_CSC = rfile.Get("positioneta12to21displacedmuonfrac%dpt%dvseta"%(fractions[0], pt))
    eff_h_all = rfile.Get("hybrideta12to21displacedmuonfrac%dpt%dvseta"%(fractions[0], pt))
    eff_st1 = rfile.Get("hybrid_noGE21eta12to21displacedmuonfrac%dpt%dvseta"%(fractions[0], pt))
    eff_ME0 = rfile.Get("hybrid3eta20to24displacedmuonfrac%dpt%dvseta"%(fractions[0], pt))
    
    eff_CSC2 = rfile.Get("positioneta20to24displacedmuonfrac%dpt%dvseta"%(fractions[0], pt))
    eff_CSC2_v2 = getPartOfEff(eff_CSC2, 2.16, 2.4)
    eff_CSC.Add(eff_CSC2_v2)
    eff_h_all = getPartOfEff(eff_h_all, 1.6, 2.15)
    eff_st1 = getPartOfEff(eff_st1, 1.6, 2.15)
    eff_ME0.Add(getPartOfEff(eff_h_all, 1.6, 2.06))
    

    
    #ratehist_CSC.Scale(SF); ratehist_st1.Scale(SF); ratehist_h_all.Scale(SF); ratehist_ME0.Scale(SF); ratehist_CSC_q12.Scale(SF)
    #SF140 = 0.7
    #ratehist_CSC.Scale(SF140); ratehist_st1.Scale(SF140); ratehist_h_all.Scale(SF140); ratehist_ME0.Scale(SF140); ratehist_CSC_q12.Scale(SF140)
    print "Rate number CSC ",ratehist_CSC.Integral()," GE11 ",ratehist_st1.Integral()," GE11GE21 ",ratehist_h_all.Integral()," GE11GE21ME0 ",ratehist_ME0.Integral()
    text4 = "p_{T}^{L1} >%d GeV, before applying track veto"%(pt)
    plotshistograms([ratehist_CSC, ratehist_st1, ratehist_h_all, ratehist_ME0], legsall, "L1 Muon with Q>=%d"%quality, text4, SBPlots+"TriggerRate_combiend_20170625_pt%d_St2eta%dto%d_withME0_q%d_allnpar_PU200"%(pt, int(etamin*10),int(etamax*10), quality), False)



    simpt = pt+5
    #texteff = "p_{T}^{L1} >%d GeV, p_{T}^{Sim}>%d GeV"%(pt, simpt)
    texteff = "#splitline{p_{T}^{L1} >%d GeV, p_{T}^{Sim}>%d GeV}{Before applying track veto}"%(pt, simpt)
    #makeplotsEta([eff_CSC, eff_st1, eff_h_all, eff_ME0], legsall, "L1 Muon with Q>=12", texteff, SBPlots+"GEMCSC_promptMu_combined3_effvsetaf_20170621_pt%d_simpt%d_eta%dto%d_allnpar_v3"%(pt, pt+5, int(1.6*10),int(2.4*10))) 

    #makeplotsPt([Teff_prompt1[0], Teff_prompt2[0], Teff_displaced], legs, legheader, text,SBPlots+"GEMCSC_dxy10to20_30to50_CurrentL1andDisplaced_eff_20170131_pt%d_fraction%d_St2eta%dto%d_allnpar"%(pt, fraction, int(etamin*10),int(etamax*10))) 


def produceME0GE21LUTvalue(pts, netas, fractions, signalsample, SBPlots):
    xaxis = "(deltay123_test)"
    yaxis1 = "(dphi_dir_ME11only_ME21only_L1)"
    yaxis2 = "(dphi_dir_ME11only_st2_L1)"
    yaxis3 = "(dphi_dir_st1_st2_L1)"
    legs3_v2 = ["Position-based algorithm","Stub alignment algorithm","Hybrid algorithm"]
    legs3 = ["Position-based algorithm","Hybrid algorithm w/o GE21","Hybrid algorithm w GE21"]
    legsall = ["Position-based algorithm","Hybrid algirthm+ w/o ME11+ w/o GE21","Hybrid algorithm+ w/o ME0 + w GE21"," Hybrid algorithm+ w ME0 + w GE21"]
    legsall = ["Position-based algorithm","Stub alignment algorithm+ w/o ME0 + w/o GE21","Stub alignment algorithm+ w/o ME0 + w GE21","Stub alignment algorithm+ w ME0 + w GE21"]


    GEMCSCTrackCh0 = TChain("GEMCSCAnalyzer/trk_eff_CSC_ALL")
    GEMCSCME31Ch0 = TChain("GEMCSCAnalyzer/trk_eff_CSC_ME31")
    #GEMCSCTrackCh1 = TChain("GEMCSCAnalyzer/trk_eff_CSC_ALL")
    #addFilesToChain(GEMCSCTrackCh1, signalsample) 
    addFilesToChain(GEMCSCTrackCh0, signalsample) 
    addFilesToChain(GEMCSCME31Ch0, signalsample) 
    GEMCSCTrackCh0.AddFriend(GEMCSCME31Ch0)
    alldPhiLUT = {}
    allddYLUT = {}
    ipt = -1
    for pt in pts:
     ipt += 1
     for fraction in fractions:
      for neta in range(len(netas)-1):
       if (netas[neta]<0 or netas[neta+1]<0):
	    continue
       ring = 1
       if (netas[neta]<1.6):
	    ring = 2
    #if (netas[0]==1.6):
       etamin = netas[neta]
       etamax = netas[neta+1]
       ddYcutlist = []; dphicutlist = []
       for npar in range(0, 4):
	    me11 = evenodds[npar].split(',')[0]
	    me21 = evenodds[npar].split(',')[1]
	    chambers = "ME1%d %s,ME2%d %s,ME3%d %s"%(ring, me11, ring, me21, ring, me21)
	    me31_n = 1
	    if me21 == "even":
		    me31_n = 2
	    ME31cut = " && ((trk_eff_CSC_ME31.has_lct&%d)>0)"%(me31_n)
	    #st_title = ["Prompt muon, 2<p_{T}<%d"%pt1, "Displaced Muon, 10<|d_{xy}|<50, p_{T}>%d"%pt]
	    ptcuts = "&& pt>%f && pt<%f"%(pt-1, pt+1)
	    cutsbase = "has_gmtCand>0 && isMe0Region && fabs(eta_st2_sim)>%f && fabs(eta_st2_sim)<%f && npar_sim==%d"%(etamin,etamax, npar)+ME31cut

	    #dphi_cut1 = getCut(GEMCSCTrackCh0, "abs("+ yaxis1 +")", cutsbase+ptcuts, 0.0, 1.5, 0.005, fraction*1.0/100.0)
	    #dphi_cut2 = getCut(GEMCSCTrackCh0, "abs("+yaxis2+")", cutsbase+ptcuts, 0.0, 1.5, 0.005, fraction*1.0/100.0)
	    dphi_cut3 = getCut(GEMCSCTrackCh0, "abs("+yaxis3+")", cutsbase+ptcuts, 0.0, 1.5, 0.005, fraction*1.0/100.0)
	    ddY_cut = getCut(GEMCSCTrackCh0, "abs("+xaxis+")", cutsbase+ptcuts, 0.0, 10, 0.01, fraction*1.0/100.0)
	    if (pt == pts[0]):
		    ddYcutlist.append(ddY_cut); dphicutlist.append(dphi_cut3)
	    else:
		    if dphi_cut3>alldPhiLUT["pt%dfrac%dnpar%d"%(pts[ipt-1], fraction, npar)]:
			    dphicutlist.append(alldPhiLUT["pt%dfrac%dnpar%d"%(pts[ipt-1], fraction, npar)])
	            else:
			    dphicutlist.append(dphi_cut3)
		    if ddY_cut>allddYLUT["pt%dfrac%dnpar%d"%(pts[ipt-1], fraction, npar)]:
			    ddYcutlist.append(allddYLUT["pt%dfrac%dnpar%d"%(pts[ipt-1], fraction, npar)])
	            else:
			    ddYcutlist.append(ddY_cut)
	    alldPhiLUT["pt%dfrac%dnpar%d"%(pt, fraction, npar)] = dphicutlist[-1]
	    allddYLUT["pt%dfrac%dnpar%d"%(pt, fraction, npar)] = ddYcutlist[-1]

			    
	    #print "pt ",pt," fraction ", fraction, " dphicut1 ",dphi_cut1," dphicut2 ",dphi_cut2, " dphicut3 ",dphi_cut3, " ddY_cut ",ddY_cut
#Positioneta12to14pt5fraction76
#{
#2.489333,1.020008,2.058920,1.162629
#}
#Directioneta12to14pt5fraction76
       #print "pt ",pt," fraction ", fraction," etamin ",etamin, " etamax ",etamax
       #print "ddY cut ",ddYcutlist," dphicut ",dphicutlist
       print "Directioneta%dto%dpt%dfraction%d"%(int(etamin*10), int(etamax*10),pt,fraction)
       print "{\n", str(dphicutlist)[1:-1],"\n}"
       if (etamin>=2.2):
	  print "Positioneta%dto%dpt%dfraction%d"%(int(etamin*10), int(etamax*10),pt,fraction)
	  print "{\n", str(ddYcutlist)[1:-1],"\n}"

       	
	    
	   

#DecompresZipFiles("./")
#SBPlots = "eta165to21SBPlots_20171120_oldpromptMu_PU200_allwithGE21_checkddY_checkdphi/"
#SBPlots = "eta165to21SBPlots_20171120_displaced_dxy10_50_PU200_allwithGE21_checkddY_checkdphi/"
#SBPlots = "eta165to21SBPlots_20171122_displacedMuonGun_allwithGE21_checkxy/"
SBPlots = "eta165to21SBPlots_20171122_darkSUSY_allwithGE21_checkxyanddxy/"
if len(sys.argv)>=2:
        SBPlots = sys.argv[1]
	if SBPlots[-1] != "/":
	    SBPlots = SBPlots+"/"
	print SBPlots
if not os.path.exists(SBPlots):
        os.makedirs(SBPlots)

filedir = "condor_20170303_rate_fraction90_95_99_eta16to215_noGE21/"
fractions = [x for x in range(76, 97)]
pts = [5, 7, 10,12, 15, 20]
netas = [1.6,1.8,2.0,2.15]
netas = [1.2,1.4,1.6]
netas = [1.2,1.4,1.6,1.8,2.0,2.1]
netas = [1.65, 1.8,2.0, 2.1]
#netas = [2.1,2.2,2.4]
ME0minEta = 2.1
#netas = [1.2,2.4]
useGE21 = False
condoroutroot = "Teffout_eta16to215_noGE21.root"
outputroot = "SB_eta16to215_test.root"
#getallPtEffAndRate(pts, fractions, netas, outputroot, filedir, condoroutroot, SBPlots, useGE21, True, True, True)
filename = "SB_eta12to215_20170306.root"
#CombinedGE21_SB(filename, pts, netas, SBPlots)
#netas = [1.2,1.4,1.6,1.8,2.0,2.15]
#CombinedGE21_effAndrate(filename, pts, netas, SBPlots)
#allLUTfile = "finalLUT_20170312_all.log"
allLUTfile = "finalLUT_20170619_ME0.log"
#allLUTfile = "finalLUT_20170727.log"
outputrootfile = "SB_eta%dto%d_20170625_test.root"%(int(10*netas[0]), int(10*netas[-1]))
outputrootfile = SBPlots+outputrootfile
#signalsample = "/Users/taohuang/Documents/CSCTrigger/GEMCSCbending_20170221/GEMCSC_Ana_ctau0_Pt2_50_PU140_20170712_GE11GE21/out_ana_1.root"
#for frac in fractions:
#    plotalletaallptallfractions(pts, netas, fractions, frac, outputrootfile, ratesample, signalsample, allLUTfile, SBPlots)
frac = 96
fractions = [96]
#getallTriggerRateAndEffVsEta(outputrootfile, [frac], 10, netas, SBPlots)
#plotalletaallptallfractions([10], netas, [96], frac, outputrootfile, ratesample, signalsample, allLUTfile, SBPlots)
plotalleta(10, netas, frac, outputrootfile, ratesample, signalsample, allLUTfile, SBPlots)


#netas = [2.1, 2.2, 2.4]
#fractions = [x for x in range(76, 97)]

#produceME0GE21LUTvalue(pts, netas, [frac], signalsample, SBPlots)
#plotalletaME0GE21(10, netas, frac, outputrootfile, ratesample, signalsample, allLUTfile, SBPlots)
#pts = [12, 15, 20]
#plotalletaallptallfractionsME0GE21(pts, netas, [frac], frac, outputrootfile, ratesample, signalsample, allLUTfile, SBPlots)
#plotalletaallptallfractions(pts, netas, [frac], frac, outputrootfile, ratesample, signalsample, allLUTfile, SBPlots)
#plotalletaallptallfractions([10], netas, fractions, frac, outputrootfile, ratesample, signalsample, allLUTfile, SBPlots)
#CompareEffPromptAndDisplaced([10, 20], netas, frac, outputrootfile, signalsample, SBPlots)

