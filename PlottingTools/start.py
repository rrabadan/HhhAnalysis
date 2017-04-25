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
ROOT.gStyle.SetPadRightMargin(0.04)
ROOT.gStyle.SetPadTopMargin(0.06)
ROOT.gStyle.SetPadBottomMargin(0.13)

ROOT.gStyle.SetLegendBorderSize(0)
ROOT.gStyle.SetLegendFillColor(0)
legend = ROOT.TLegend(0.85,0.6,0.95,0.94)
legend.SetFillColor(ROOT.kWhite)
legend.SetTextSize(0.02)
legend.SetTextFont(62)

color  = {"TT":ROOT.kRed, "DY":ROOT.kBlue, "VV":ROOT.kGreen+2, "sT":ROOT.kMagenta+2, "Wjet":ROOT.kOrange+2, "ttV":ROOT.kViolet-1, "Rad_260":9, "Rad_500":28, "Rad_900":12, "Data": 0 }
marker = {"TT":21, "DY":22, "VV":23, "sT":34, "Wjet":33, "ttV":29 }

import getpass
user = getpass.getuser()
