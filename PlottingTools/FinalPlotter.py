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

def drawVariabes():
  c1 = ROOT.TCanvas()
  c1.SetGridx(); c1.SetGridy(); c1.SetTickx(); c1.SetTicky(); c1.cd()
  if(log=="log"): c1.SetLogy(1)
  else:           c1.SetLogy(0)
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
        if(log=="log"): h.SetMinimum(0.00001)
        if(Norm=="uni" and not h.GetName() in VetoList):  h.Scale(1./h.GetEntries())
        # Each histo is drawn alone in the Sample Folder
        h.Draw(); c1.SaveAs(Folder + "/" + Samples[nFile] + "/" + h.GetName() + ".pdf")
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
            if(Norm=="lumi"): h.Scale(Lumi)
            if( nHist==0 ): legend.AddEntry(h, Samples[nFile], "l") # One entry per file
            hSign.append(h)
            hSign[nHist].SetLineColor(color[Samples[nFile]]); hSign[nHist].SetLineStyle(2);
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
  print "Drawing the histograms."
  nHist = 0
  for this_stack in hstack:
    this_stack.SetMaximum(hstack_max[nHist])
    pic_name = this_stack.GetName()
    if not DataMC: this_stack.Draw("nostack") # Histos overimposed
    else:
      this_stack.Draw("hist");
      hdata[nHist].Draw("PEsame");
      for iS in range(len(hSignAll)): hSignAll[iS][nHist].Draw("Lsame")
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

drawVariabes()
drawSF()
print "The end! Thanks for choosing: ", sys.argv[0]
