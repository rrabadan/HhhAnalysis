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
Format   = [".pdf",".C"]
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

def drawSF():
  c1 = ROOT.TCanvas()
  c1.SetGridx(); c1.SetGridy(); c1.SetTickx(); c1.SetTicky(); c1.cd()
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
        if this_format == ".C":
          nameH = Folder + "/C/" + hist.GetName() + "_" + Samples[nFile] + this_format
        c1.SaveAs(nameH)
    nFile +=1

c0 = ROOT.TCanvas("c0","",800, 800)
c0.SetGridx(); c0.SetGridy(); c0.SetTickx(); c0.SetTicky(); c0.cd()
if(log=="log"): c0.SetLogy(1)
else:           c0.SetLogy(0)
hdata = []; hSignAll = []
hstack = []; hstack_max = []; Xaxis = []; Yaxis = []
nFile = 0
LastHisto = 0
# Find the total number of histos to plot
for h in File_List[0].GetListOfKeys():
  if( h.ClassName()=="TH1F" and not h.GetName() in VetoList ): LastHisto += 1
# For each file
for ThisFile in File_List:
  ThisFile.cd()
  nHist = 0
  hSign = []
  # Loop over its object
  for h in ThisFile.GetListOfKeys():
    h = h.ReadObj()
    # And only take care of TH1F
    if( h.ClassName()=="TH1F"):
      if(Norm=="uni" and not h.GetName() in VetoList):  h.Scale(1./h.GetEntries())
      # Each histo is drawn alone in the Sample Folder
      c0.cd(); h.Draw(); c0.SaveAs(Folder + "/" + Samples[nFile] + "/" + h.GetName() + ".pdf")
      # The not vetoed histos are drawn all together
      if( not h.GetName() in VetoList ):
        # IF is MC (or you do not want a stack plot)
        if( nFile<DataSignal_startAt or not DataMC ):
          if(Norm=="lumi"): h.Scale(Lumi)
          if( nHist==0 ): legend.AddEntry(h, Samples[nFile], "l") # One entry per File
          if( nFile==0 ):                                         # One per Histo
            hstack.append(ROOT.THStack(h.GetName(),""))
            hstack_max.append(-1)
          Xaxis.append(h.GetXaxis().GetTitle())
          Yaxis.append(h.GetYaxis().GetTitle())
          h.SetLineColor(color[Samples[nFile]]); h.SetFillColor(color[Samples[nFile]]); h.SetMarkerColor(color[Samples[nFile]]); h.SetMarkerStyle(marker[Samples[nFile]]); h.SetLineWidth(2);
          hstack[nHist].Add(h)
          if( h.GetMaximum() > hstack_max[nHist] ): hstack_max[nHist] = h.GetMaximum()
        # IF is SIGNAL (and you want a stack plot)
        if( nFile>=DataSignal_startAt and nFile!=int(len(File_List)-1) and DataMC ):
          if(Norm=="lumi"): h.Scale(Lumi*50)
          if( nHist==0 ): legend.AddEntry(h, Samples[nFile], "l") # One entry per file
          hSign.append(h)
          hSign[nHist].SetLineColor(color[Samples[nFile]]); hSign[nHist].SetLineStyle(2); hSign[nHist].SetLineWidth(2);
          if( hSign[nHist].GetMaximum() > hstack_max[nHist] ): hstack_max[nHist] = hSign[nHist].GetMaximum()
          if (nHist==LastHisto-1) : hSignAll.append(hSign); 
        # IF is Data (ALWAYS the last one, and you want a stack plot)
        if( nFile==int(len(File_List)-1) and DataMC ):
          hdata.append(h)
          hdata[0].SetLineColor(1); hdata[0].SetMarkerStyle(20); hdata[0].SetMarkerColor(1);
        nHist += 1
  nFile += 1

# Now Draw The histos
ROOT.gStyle.SetOptStat(0)
c1 = ROOT.TCanvas("c1","",800, 800)
c1.SetGridx(); c1.SetGridy(); c1.SetTickx(); c1.SetTicky(); c1.cd()
if(log=="log"): c1.SetLogy(1)
else:           c1.SetLogy(0)
pad1 = ROOT.TPad("pad1", "", 0, 0.2, 1, 1.0)
pad1.SetBottomMargin(0) # Upper and lower plot are joined
pad1.SetGridx();        # Vertical grid
pad1.Draw();            # Draw the upper pad: pad1
pad1.cd();              # pad1 becomes the current pad
nHist = 0
for this_stack in hstack:
  if(log=="log"): this_stack.SetMinimum(0.00001)
  this_stack.SetMaximum(hstack_max[nHist])
  pic_name = this_stack.GetName()
  if not DataMC: this_stack.Draw("nostack") # Histos overimposed
  else:
    this_stack.Draw("hist");
    hdata[nHist].SetMarkerStyle(0); hdata[nHist].SetMarkerColor(ROOT.kBlack);
    hdata[nHist].Draw("PEsame");
    for iS in range(len(hSignAll)): hSignAll[iS][nHist].Draw("same")
  legend.Draw("same")
  this_stack.GetXaxis().SetTitle(Xaxis[nHist]); this_stack.GetXaxis().SetTitleSize(1.);
  this_stack.GetYaxis().SetTitle(Yaxis[nHist]); this_stack.GetYaxis().SetTitleSize(1.);
  c1.cd(); # Go back to the main canvas before defining pad2
  pad2 = ROOT.TPad("pad2", "", 0, 0.05, 1, 0.2)
  pad2.SetTopMargin(0)
  pad2.SetBottomMargin(0.2)
  pad2.SetGridx(); pad2.SetGridy();
  pad2.Draw()
  pad2.cd() # pad2 becomes the current pad
  # Define the ratio plot
  h3 = hdata[nHist].Clone("h3")
  h3.SetLineColor(ROOT.kBlack)
  h3.SetMinimum(0.5); h3.SetMaximum(1.5);
  h3.Sumw2()
  h3.SetStats(0);
  h3.GetXaxis().SetTitle(Xaxis[nHist]); h3.GetXaxis().SetTitleOffset(0.5); h3.GetXaxis().SetTitleSize(0.2);
  h3.GetYaxis().SetTitle("Data/MC ratio"); h3.GetXaxis().SetTitleOffset(0.5); h3.GetYaxis().SetTitleSize(0.2);
  h3.GetYaxis().SetLabelSize(0.08) 
  h3.GetXaxis().SetLabelSize(0.1)
  h_tot = ROOT.TH1F()
  index=0
  for h_temp in this_stack.GetHists():
    if(index==0): h_tot = h_temp
    else:         h_tot.Add(h_temp)
    index += 1
  h3.Divide(h_tot)
  h3.SetMarkerStyle(2)
  h3.Draw("ep")       # Draw the ratio plot
  for this_format in Format:
    nameH = Folder + pic_name + this_format
    if this_format == ".C":
      nameH = Folder + "/C/" + pic_name + this_format
    c1.SaveAs( nameH )
  nHist = nHist +1
  del h3

#drawSF()
print "The end! Thanks for choosing: ", sys.argv[0]
