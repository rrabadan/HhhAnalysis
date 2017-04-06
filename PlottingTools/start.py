ROOT.gROOT.SetBatch(1)

ROOT.gStyle.SetStatW(0.17)
ROOT.gStyle.SetStatH(0.15)

ROOT.gStyle.SetTitleStyle(0)
ROOT.gStyle.SetTitleAlign(13) ## coord in top left
ROOT.gStyle.SetTitleX(0.)
ROOT.gStyle.SetTitleY(1.)
ROOT.gStyle.SetTitleW(1)

ROOT.gStyle.SetTitleXSize(0.05)
ROOT.gStyle.SetTitleYSize(0.05)
ROOT.gStyle.SetTitleH(0.058)
ROOT.gStyle.SetTitleBorderSize(0)

ROOT.gStyle.SetPadLeftMargin(0.126)
ROOT.gStyle.SetPadRightMargin(0.14)
ROOT.gStyle.SetPadTopMargin(0.06)
ROOT.gStyle.SetPadBottomMargin(0.13)

legend = ROOT.TLegend(0.75,0.6,0.86,0.94)
legend.SetFillColor(ROOT.kWhite)
legend.SetTextSize(0.05)
legend.SetTextFont(62)

color  = {"TT":ROOT.kRed, "DY":ROOT.kBlue, "VV":ROOT.kGreen+2, "sT":ROOT.kMagenta+2, "Wjet":ROOT.kOrange+2, "ttV":ROOT.kViolet-1, "Data": 0 }
marker = {"TT":21, "DY":22, "VV":23, "sT":34, "Wjet":33, "ttV":29, "Data": 0 }
