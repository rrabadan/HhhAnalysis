import random
import ROOT
from ROOT import TFile,TChain,TH1F,TH2F,TLegend,gDirectory
import os
execfile("start.py")

#Creating folders and parameters
Filefolder   = "/fdata/hepx/store/user/lpernie/Hhh_For_Plotting/"
Lumi     = 36.42#fb-1
Samples  = ["ttV","Wjet","sT","VV","DY","TTbar","Data"]
#Samples  = ["ttV","Wjet","VV","DY","Data"]

DataOrMC = "DataMC"
Norm     = "uni"#"lumi"
log      = "log"#"lin"
Format   = [".pdf",".C"]
VetoList = ["h_XsecBr","h_muon1_triggerSF","h_muon1_isoSF","h_muon1_idSF","h_muon1_trackingSF","h_muon2_triggerSF","h_muon2_isoSF","h_muon2_idSF","h_muon2_trackingSF","h_Nev_preHLT","h_Nev_posHLT"]
print "You are producing ",Norm,log,Format,DataOrMC

# List of files
File_List = []
for sample in Samples:
  Myfile = ROOT.TFile.Open(Filefolder+"/"+sample+".root","read");
  File_List.append(Myfile)

# For each root file you loop over the histograms it contains
c1 = ROOT.TCanvas()
c1.SetGridx(); c1.SetGridy(); c1.SetTickx(); c1.SetTicky(); c1.cd()
if(log=="log"): c1.SetLogy(1)
else:           c1.SetLogy(0)
hdata = ROOT.TH1F()
hstack = []
hstack_max = []
nFile = 0
for ThisFile in File_List:
  ThisFile.cd()
  nHist = 0
  for h in Myfile.GetListOfKeys():
    h = h.ReadObj()
    if( h.ClassName()=="TH1F" and not h.GetName() in VetoList ):
      # If is Data (ALWAYS the last one)
      if( nFile==int(len(File_List)-1) ):
        hdata = h
        hdata.SetLineColor(1); hdata.SetMarkerStyle(20); hdata.SetMarkerColor(1);
        if(Norm=="uni"): h.Scale(1./h.GetEntries())
      # If is MC
      else:
        if(Norm=="lumi"): h.Scale(Lumi)
        if(Norm=="uni"):  h.Scale(1./h.GetEntries())
        if( nHist==0 ): legend.AddEntry(h, Samples[nFile], "l") # One entry per file
        if( nFile==0 ): # One entry per Histo
          hstack.append(ROOT.THStack("hs_"+h.GetName(),""))
          hstack_max.append(-1)
        h.SetLineColor(color[Samples[nFile]]); h.SetFillColor(color[Samples[nFile]]); h.SetMarkerColor(color[Samples[nFile]]); h.SetMarkerStyle(marker[Samples[nFile]])
        h.SetLineWidth(2); h.SetMinimum(0.00001)
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
  else:                  this_stack.Draw("hist"); hdata.Draw("Psame")
  legend.Draw("same")
  this_stack.GetXaxis().SetTitle(ROOT.TH1F(this_stack.GetHistogram()).GetXaxis().GetTitle())
  this_stack.GetYaxis().SetTitle("Y")
  Folder   = "Plots_" + DataOrMC + "_" + Norm + "_" + log + "/"
  os.system("mkdir -p " + Folder + "/C")
  for this_format in Format:
    if this_format == ".C": Folder = Folder + "C/"
    nameH = Folder + pic_name + this_format
    print "Saving", nameH
    c1.SaveAs( nameH )
    print "Saved"
  nHist = nHist +1
