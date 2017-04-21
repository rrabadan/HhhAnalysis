import ROOT
from ROOT import TFile,TChain,TH1F,TH2F,TLegend,gDirectory
import sys, os, random
execfile("start.py")

#Creating folders and parameters
Filefolder   = "/fdata/hepx/store/user/%s/Hhh_For_Plotting/"%user
Lumi     = 36.42#fb-1
Samples  = ["ttV","Wjet","sT","VV","DY","TT","Rad_260","Rad_500","Rad_900","Data"]
DataSignal_startAt = 6 # The index where you want to start not doing a StackPlot (for Signal and Data).

DataMC   = True
Norm     = "lumi"#"uni"
log      = "log"#"lin"
Format   = [".pdf",".png",".C"]
VetoList = ["h_XsecBr","h_muon1_triggerSF","h_muon1_isoSF","h_muon1_idSF","h_muon1_trackingSF","h_muon2_triggerSF","h_muon2_isoSF","h_muon2_idSF","h_muon2_trackingSF","h_Nev_preHLT","h_Nev_posHLT"]
SFtodraw = ["h_pre_muon1_triggerSF","h_pre_muon1_isoSF","h_pre_muon1_idSF","h_pre_muon1_trackingSF","h_pre_muon2_triggerSF","h_pre_muon2_isoSF","h_pre_muon2_idSF","h_pre_muon2_trackingSF"]
if( (Norm!="uni" and Norm!="lumi") or (log!="log" and log!="lin") ): print "WARNING!!! Wrong paramters."; sys.exit()
DataOrMC = "DataMC"
if (not DataMC): DataOrMC = "OverImposed"
if DataMC: print "You are producing Starck plots in", log, "scale. Normalized to:", Norm, ". Format:", Format
else:      print "You are producing OverImposed plots in", log, "scale. Normalized to:", Norm, ". Format:", Format

# List of files
print "Loading files in:", Filefolder
Folder   = "Plots_" + DataOrMC + "_" + Norm + "_" + log + "/"
os.system("mkdir -p " + Folder + "/C")
File_List = []
for sample in Samples:
  Myfile = ROOT.TFile.Open(Filefolder+"/"+sample+".root","read");
  File_List.append(Myfile)
  os.system("mkdir -p " + Folder + "/" + sample)

# For each root file you loop over the histograms it contains
print "Looping over Files and Histograms."
ROOT.gStyle.SetOptStat(0)

def drawVariabes():
    c1 = ROOT.TCanvas()
    c1.SetGridx(); c1.SetGridy(); c1.SetTickx(); c1.SetTicky(); c1.cd()
    if(log=="log"): c1.SetLogy(1)
    else:           c1.SetLogy(0)
    hdata = []
    hstack = []; hstack_max = []; Xaxis = []; Yaxis = []
    nFile = 0
    for ThisFile in File_List:
      ThisFile.cd()
      nHist = 0
      for h in ThisFile.GetListOfKeys():
	h = h.ReadObj()
	if( h.ClassName()=="TH1F"):
	  if(log=="log"): h.SetMinimum(0.00001)
	  if(Norm=="uni" and not h.GetName() in VetoList):  h.Scale(1./h.GetEntries())
	  h.Draw(); c1.SaveAs(Folder + "/" + Samples[nFile] + "/" + h.GetName() + ".pdf")
	  if( not h.GetName() in VetoList ):
	    # If is Data (ALWAYS the last one)
	    if( nFile==int(len(File_List)-1) and DataMC ):
	      hdata.append(h)
	      hdata[nHist].SetLineColor(1); hdata[nHist].SetMarkerStyle(20); hdata[nHist].SetMarkerColor(1);
	      if( hdata[nHist].GetMaximum() > hstack_max[nHist] ): hstack_max[nHist] = hdata[nHist].GetMaximum()
	    # If is MC
	    else:
	      if(Norm=="lumi"): h.Scale(Lumi)
	      if( nHist==0 ): legend.AddEntry(h, Samples[nFile], "l") # One entry per file
	      if( nFile==0 ): # One entry per Histo
		hstack.append(ROOT.THStack(h.GetName(),""))
		hstack_max.append(-1)
		Xaxis.append(h.GetXaxis().GetTitle())
		Yaxis.append(h.GetYaxis().GetTitle())
	      h.SetLineColor(color[Samples[nFile]]); h.SetFillColor(color[Samples[nFile]]); h.SetMarkerColor(color[Samples[nFile]]); h.SetMarkerStyle(marker[Samples[nFile]]); h.SetLineWidth(2);
	      hstack[nHist].Add(h)
	      if( h.GetMaximum() > hstack_max[nHist] ): hstack_max[nHist] = h.GetMaximum()
	    nHist = nHist + 1
      nFile = nFile +1

    # Now Draw The histos
    ROOT.gStyle.SetOptStat(0)
    print "Drawing the histograms."
    nHist = 0
    for this_stack in hstack:
      this_stack.SetMaximum(hstack_max[nHist])
      pic_name = this_stack.GetName()
      if not DataMC: this_stack.Draw("nostack") # Histos overimposed
      else:          this_stack.Draw("hist"); hdata[nHist].Draw("PEsame")
      legend.Draw("same")
      this_stack.GetXaxis().SetTitle(Xaxis[nHist])
      this_stack.GetYaxis().SetTitle(Yaxis[nHist])
      for this_format in Format:
	nameH = Folder + pic_name + this_format
	if this_format == ".C":
	  nameH = Folder + "/C/" + pic_name + this_format
	c1.SaveAs( nameH )
      nHist = nHist +1

def drawSF():
    c1 = ROOT.TCanvas()
    c1.SetGridx(); c1.SetGridy(); c1.SetTickx(); c1.SetTicky(); c1.cd()
    #draw SF
    nFile  = 0
    for ThisFile in File_List:
      ThisFile.cd()
      for sfname in SFtodraw:
	    hist = ThisFile.Get(sfname)
	    hist.SetMaximum(1.01)
	    hist.SetMinimum(.9)

	    if hist.ClassName() == "TH2F":
	       hist.Draw("colztext")
	    elif hist.ClassName() == "TH1F":
	       hist.SetMarkerSize(2)
	       hist.Draw("histp")
	    for this_format in Format:
	       nameH = Folder + hist.GetName() + "_" + Samples[nFile] + this_format
	       c1.SaveAs(nameH)
      nFile +=1
#drawVariabes()
drawSF()
print "The end! Thanks for choosing: ", sys.argv[0]
