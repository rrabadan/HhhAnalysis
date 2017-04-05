import random
import ROOT
from ROOT import TFile,TChain,TH1F,TH2F,TLegend,gDirectory
import os
execfile("start.py")

#Creating folders and parameters
Filefolder   = "/fdata/hepx/store/user/lpernie/Hhh_For_Plotting/"
Lumi     = 36.42#fb-1
#Samples  = ["TT"]#,"DY","VV","singTop","Wjet","ttV","Data"]
Samples  = ["TT","Wjet"]
DataOrMC = "DataMC"
Norm     = "lumi"
Format   = [".pdf",".C"]

# List of files
File_List = []
for sample in Samples:
  Myfile = ROOT.TFile.Open(Filefolder+"/"+sample+".root","read");
  File_List.append(Myfile)

# For each root file you loop over the histograms it contains
hdata = ROOT.TH1F()
hdata.SetLineColor(1); hdata.SetMarkerStyle(20); hdata.SetMarkerColor(1);
hstack = []
hstack_max = []
nFile = 0
for ThisFile in File_List:
  ThisFile.cd()
  nHist = 0
  for h in Myfile.GetListOfKeys():
    h = h.ReadObj()
    if( h.ClassName()=="TH1F"):
      # If is Data (ALWAYS the last one)
      if( nFile==int(len(File_List)-1) ):
        hdata = h
      # If is MC
      else:
        if( nHist==0 ): legend.AddEntry(h, Samples[nFile], "l") # One entry per file
        if( nFile==0 ): # One entry per Histo
          hstack.append(ROOT.THStack("hs_"+h.GetName(),""));
          hstack_max.append(-1)
        hstack[nHist].Add(h)
        if( h.GetMaximum() > hstack_max[nHist] ): hstack_max[nHist] = h.GetMaximum()
      nHist = nHist + 1
  nFile = nFile +1

# Now Draw The histos
nHist = 0
for this_stack in hstack:
  this_stack.SetMaximum(hstack_max[nHist])
  pic_name = this_stack.GetName()
  if DataOrMC!="DataMC": this_stack.Draw("nostack") # Histos overimposed
  else:                  this_stack.Draw("hist"); hdata.Draw("same")
  legend.Draw("same")
  for this_format in Format:
    Folder   = "Plots_" + DataOrMC + "_" + Norm + "/"
    os.system("mkdir -p " + Folder + "/C")
    if this_format == ".C": Folder = Folder + "C/"
    c1.SaveAs( Folder + pic_name + this_format )

  nHist = nHist +1
