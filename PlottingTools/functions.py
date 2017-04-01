def hist1D(tree, todraw, x_bins, cut, B, Lumi, nTOT, isMC):
  if cut=="" or cut==" ": cut="1"
  Lumi    = Lumi * 1000 # Convert from fb-1 to pb-1
  xBins   = int(x_bins[1:-1].split(',')[0])
  xminBin = float(x_bins[1:-1].split(',')[1])
  xmaxBin = float(x_bins[1:-1].split(',')[2])
  b1      = ROOT.TH1F("%s_%s"%(B,todraw), "%s"%B, xBins, xminBin, xmaxBin)
  if isMC:
    cut_and_weight = str(Lumi) + "*(XsecBr/" + str(nTOT) + ")*(" + str(cut) + ")"
  else:
    if nTOT!=1: print 'WARNING!!! DATA have "nTOT" that is different from 1. nTOT is the totoal number of MC event geenrated, that should be set to 1 for DATA.'
    #cut_and_weight = str(Lumi) + "*((1.05*XsecBr)/" + str(nTOT) + ")*(" + str(cut) + ")" #Use this if you want to make a test and use a MC sample as if it was DATA
    cut_and_weight = str(cut)
  tree.Draw("%s>>%s_%s"%(todraw,B,todraw), cut_and_weight)
  ROOT.SetOwnership(b1, False)
  return b1

def draw1D_v2(filelist,x_bins,x_title,cut,benchmarks, pic_name):
  xBins = int(x_bins[1:-1].split(',')[0])
  xminBin = float(x_bins[1:-1].split(',')[1])
  xmaxBin = float(x_bins[1:-1].split(',')[2])
  b1 = ROOT.TH1F("b1","b1",xBins,xminBin,xmaxBin)
  b1.SetTitle("h2#rightarrow hh#rightarrow WWWW"+" "*24 + "14TeV")
  b1.GetYaxis().SetTitle("Events")
  b1.GetXaxis().SetTitle("%s"%x_title)
  b1.SetStats(0)
  color = [ROOT.kRed, ROOT.kBlue, ROOT.kMagenta+2, ROOT.kGreen+2, ROOT.kCyan]
  marker = [20,21,22,23,34]
  legend = ROOT.TLegend(0.75,0.6,0.86,0.94)
  legend.SetFillColor(ROOT.kWhite)
  legend.SetTextSize(0.05)
  legend.SetTextFont(62)
  hs1 = ROOT.THStack("hs1","%s distribution"%x_title)
  hs2 = ROOT.THStack("hs2","%s distribution"%x_title)
  hs3 = ROOT.THStack("hs3","%s distribution"%x_title)
  hs4 = ROOT.THStack("hs4","%s distribution"%x_title)
  hists_1 = []
  hists_2 = []
  hists_3 = []
  hists_4 = []
  for nfile in range(len(filelist)):
    hists_1.append(ROOT.TH1F("hist1_%d"%nfile,"hist1_%d"%nfile, 190, 200, 4000))
    hists_2.append(ROOT.TH1F("hist2_%d"%nfile,"hist2_%d"%nfile, 190, 200, 4000))
    hists_3.append(ROOT.TH1F("hist3_%d"%nfile,"hist3_%d"%nfile,30,0,300))
    hists_4.append(ROOT.TH1F("hist4_%d"%nfile,"hist4_%d"%nfile,30,0,900))
    for nfile in range(len(filelist)):
      rootfile = filelist[nfile]
      B = benchmarks[nfile]
      f = ROOT.TFile(rootfile)
      lists = f.GetListOfKeys()
      for x in range(len(lists)):
        subkey = lists.At(x)
        obj = subkey.ReadObj()
        if obj.GetName()=="evtree":
          continue
      print " title ",obj.GetTitle()," Name ",obj.GetName()
      maxbin = obj.GetMaximumBin()
      hists_1[nfile].Fill(obj.GetXaxis().GetBinCenter(maxbin))
      hists_3[nfile].Fill(obj.GetBinContent(maxbin))
      hists_4[nfile].Fill(obj.Integral())
      obj.Scale(1.0/obj.Integral())
      hists_2[nfile].Add(obj.Rebin(20))
      hists_1[nfile].SetLineColor(color[nfile])
      hists_1[nfile].SetMarkerColor(color[nfile])
      hists_1[nfile].SetMarkerStyle(marker[nfile])
      
      hists_2[nfile].SetLineColor(color[nfile])
      hists_2[nfile].SetMarkerColor(color[nfile])
      hists_2[nfile].SetMarkerStyle(marker[nfile])
      
      hists_3[nfile].SetLineColor(color[nfile])
      hists_3[nfile].SetMarkerColor(color[nfile])
      hists_3[nfile].SetMarkerStyle(marker[nfile])
      
      hists_4[nfile].SetLineColor(color[nfile])
      hists_4[nfile].SetMarkerColor(color[nfile])
      hists_4[nfile].SetMarkerStyle(marker[nfile])
      if (nfile==len(filelist)-1):
        hists_1[nfile].Scale(1.0/5.0)
        hists_2[nfile].Scale(1.0/5.0)
      hists_3[nfile].Scale(1.0/hists_3[nfile].Integral())
      hists_4[nfile].Scale(1.0/hists_4[nfile].Integral())
      hs1.Add(hists_1[nfile])
      hs2.Add(hists_2[nfile])
      hs3.Add(hists_3[nfile])
      hs4.Add(hists_4[nfile])
      legend.AddEntry(hists_1[nfile], "%s"%B, "pl")
      c1 = ROOT.TCanvas()
      c1.SetGridx()
      c1.SetGridy()
      c1.SetTickx()
      c1.SetTicky()
      c1.Divide(2,2)
      c1.cd(1)
      hs1.Draw("nostack+p")
      hs1.GetHistogram().GetXaxis().SetTitle("%s"%x_title)
      hs1.GetHistogram().GetXaxis().SetRangeUser(xminBin, xmaxBin)
      legend.Draw("same")
      tex1 = ROOT.TLatex(0.25,.50,"Most probable mass")
      tex1.SetTextSize(0.05)
      tex1.SetTextFont(62)
      tex1.SetNDC()
      tex1.Draw("same")
      c1.cd(2)
      hs2.Draw("nostack+p")
      hs2.GetHistogram().GetXaxis().SetTitle("%s"%x_title)
      hs2.GetHistogram().GetXaxis().SetRangeUser(xminBin, xmaxBin)
      legend.Draw("same")
      tex2 = ROOT.TLatex(0.25,.50,"add all survival solution(normalized in each event)")
      tex2.SetTextSize(0.05)
      tex2.SetTextFont(62)
      tex2.SetNDC()
      tex2.Draw("same")
     
      c1.cd(3)
      hs3.Draw("nostack+p")
      hs3.SetTitle("bincontent of maximum bin, normalized  distribution")
      hs3.GetHistogram().GetXaxis().SetTitle("maximum bincontent")
      legend.Draw("same")
      tex3 = ROOT.TLatex(0.25,.50,"maximum  bin content from MMC")
      tex3.SetTextSize(0.05)
      tex3.SetTextFont(62)
      tex3.SetNDC()
      tex3.Draw("same")
     
      c1.cd(4)
      hs4.SetTitle("total number of survival solutions, normalized distribution")
      hs4.Draw("nostack+p")
      hs4.GetHistogram().GetXaxis().SetTitle("total number of survival solutions")
      legend.Draw("same")
      tex4 = ROOT.TLatex(0.25,.50,"total number of survival solutions")
      tex4.SetTextSize(0.05)
      tex4.SetTextFont(62)
      tex4.SetNDC()
      tex4.Draw("same")
  c1.cd()
  c1.SaveAs("Hhh_PDFvalidation_%s_combined.png"%pic_name)
    
def draw1D(filelist, todraw, x_bins, x_title,cut, benchmarks, pic_name, Lumi, nTOT, Norm, DataOrMC, LOG, Format):
  #Check parameters
  if (not(".pdf" in Format) and not(".C" in Format) ): print "WARNING: 'Format' has a wrong inzialization!"
  if( (Norm=="unity") and DataOrMC=="DataMC"):
    print "You are asking to normalize to unity and you are overimposing MC and data. This make no sense."; return;
  #Loop over the samples
  c1 = ROOT.TCanvas()
  c1.SetGridx(); c1.SetGridy(); c1.SetTickx(); c1.SetTicky()
  if(LOG=="yes"): c1.SetLogy(1)
  if(LOG=="no"):  c1.SetLogy(0)
  color = [ROOT.kRed, ROOT.kBlue, ROOT.kMagenta+2, ROOT.kGreen+2, ROOT.kCyan, ROOT.kOrange+2, ROOT.kViolet-1, ROOT.kPink+2]
  marker = [20,21,22,23,34,33,29,28]
  legend = ROOT.TLegend(0.65,0.65,0.8,0.94); legend.SetFillColor(ROOT.kWhite); legend.SetTextSize(0.05); legend.SetTextFont(62)
  hs = ROOT.THStack("hs"," ")
  hdata = ROOT.TH1F()
  hists = []
  i = 0
  for nfile in range(len(filelist)):
    print "TEST: printing", nfile, benchmarks[nfile]
    isMC = True
    if( DataOrMC=="DataMC" and (nfile==int(len(filelist)-1)) ): isMC=False
    B = benchmarks[nfile]
    hist = hist1D(filelist[nfile], todraw, x_bins, cut, B, Lumi, nTOT[i], isMC)
    hist.SetLineColor(color[nfile])
    hist.SetLineWidth(2)
    hist.SetMarkerColor(color[nfile])
    hist.SetMarkerStyle(marker[nfile])
    hist.SetMinimum(0.00001)
    print ' ->', hist.Integral()
    if(hist.Integral()<=0):
      print "-> NO events for",B,"using the selection:"
      print "   ",cut
    else:
      if isMC:
        print "IS MC"
        if(DataOrMC=="DataMC"): hist.SetFillColor(color[nfile])
        if(Norm=="unity"): hist.Scale(1./hist.Integral())
        hs.Add(hist)
        legend.AddEntry(hist, "%s"%B, "l")
        hists.append(hist)
      else: # Data is expected to be the last item of filelist
        print "IS DATA"
        hdata = hist
        hdata.SetLineColor(1); hdata.SetMarkerStyle(20); hdata.SetMarkerColor(1);
        legend.AddEntry(hist, "%s"%B, "l")
    i = i+1
  print 'Done with LOOP'
  if DataOrMC!="DataMC": hs.Draw("nostack") # One on top of the others
  else:                  hs.Draw("hist"); hdata.Draw("same")
  hs.GetHistogram().GetXaxis().SetTitle("%s"%x_title)
  if(Norm=="unity"):  hs.GetHistogram().GetYaxis().SetTitle("Normalized to unity")
  if(Norm=="lumi"):   hs.GetHistogram().GetYaxis().SetTitle("Normalized to "+str(Lumi)+" fb^{-1}")
  legend.Draw("same")
  logTXT = ""
  if(LOG=="yes"): logTXT = "log_"
  for this_format in Format:
    Folder = "Plots/"
    if this_format == ".C": Folder = "Plots/C/"
    c1.SaveAs( Folder + logTXT + Norm + "_" + DataOrMC + "_" + pic_name + this_format )

def findNewestDir(directory):
  os.chdir(directory)
  dirs = {}
  for dir in glob.glob('*'):
    if os.path.isdir(dir):
      dirs[dir] = os.path.getctime(dir)
  lister = sorted(dirs.iteritems(), key=operator.itemgetter(1))
  return lister[-1][0]
