import random
import ROOT
from ROOT import TFile,TChain,TH1F,TH2F,TLegend
import os
execfile("start.py")
execfile("functions.py")

#Creating folders and parameters
tree_name="DiHiggsWWBBAna/evtree"
os.system("mkdir -p Plots")
os.system("mkdir -p HADD")
Lumi=40#fb-1
#TChains
TT_ch = ROOT.TChain(tree_name)
os.system("find /fdata/hepx/store/user/lpernie/TTTo2L2Nu_13TeV-powheg/TTTo2L2Nu_powheg_miniAODv2_v0_ext1_01/170322_221256 | grep root | grep -v failed > HADD/TT0.txt")
with open('HADD/TT0.txt','r') as f:
  for line in f:
    if not line.isspace():
      TT_ch.Add(str(line[:-1]))
print "TT has", TT_ch.GetEntries(), "entries."
#
DY_ch = ROOT.TChain(tree_name)
os.system("find /fdata/hepx/store/user/lpernie/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/170323_021221 | grep root | grep -v failed > HADD/DY0.txt")
os.system("find /fdata/hepx/store/user/lpernie/DYToLL_0J_13TeV-amcatnloFXFX-pythia8/crab_DYToLL_0J_13TeV-amcatnloFXFX-pythia8/170323_021233 | grep root | grep -v failed >> HADD/DY0.txt")
os.system("find /fdata/hepx/store/user/lpernie/DYToLL_1J_13TeV-amcatnloFXFX-pythia8/crab_DYToLL_1J_13TeV-amcatnloFXFX-pythia8/170323_021246 | grep root | grep -v failed >> HADD/DY0.txt")
os.system("find /fdata/hepx/store/user/lpernie/DYToLL_2J_13TeV-amcatnloFXFX-pythia8/crab_DYToLL_2J_13TeV-amcatnloFXFX-pythia8/170323_021258 | grep root | grep -v failed >> HADD/DY0.txt")
with open('HADD/DY0.txt','r') as f:
  for line in f:
    if not line.isspace():
      DY_ch.Add(str(line[:-1]))
print "DY has", DY_ch.GetEntries(), "entries."

filelist = [TT_ch,DY_ch]
#Cuts and Ordering
#cut = "mu1_pt>10 && mu2_pt>10 && fabs(mu1_eta)<2.4 && fabs(mu2_eta)<2.4"
cut = ""
benchmarks = ["TTbar","DY"]

#---Starting to plot histos here-------------------------------------------------------------------------------------------------------------------------------------
pic_name = "dR_l1l2"
todraw   = "dR_l1l2"
x_bins   = "(50,0.0,5)"
x_title  = "#Delta R(l,l)"
draw1D(filelist, todraw, x_bins, x_title, cut, benchmarks, pic_name, Lumi, "unity")

pic_name = "mass_l1l2"
todraw   = "mass_l1l2"
x_bins   = "(36,10,150)"
x_title  = "m(l,l)"
draw1D(filelist, todraw, x_bins, x_title, cut, benchmarks, pic_name, Lumi, "unity")

pic_name = "mass_b1b2"
todraw   = "mass_b1b2"
x_bins   = "(36,10,150)"
x_title  = "m(j,j)"
draw1D(filelist, todraw, x_bins, x_title, cut, benchmarks, pic_name, Lumi, "unity")

pic_name = "MMC_h2mass"
x_title = "Most Probable mass"
x_bins = "(60,200,1400)"
#draw1D_v2(filelist,x_bins,x_title,cut,benchmarks, pic_name)
pic_name = "EtaonshellW"
todraw = "w1_eta*(w1_mass>w2_mass)+w2_eta*(w1_mass<w2_mass)"
x_title = "#eta of onshell W"
x_bins = "(100,-4,4)"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)

pic_name = "PhionshellW"
todraw = "w1_phi*(w1_mass>w2_mass)+w2_phi*(w1_mass<w2_mass)"
x_title = "#Phi of onshell W"
x_bins = "(100,-3.1415,3.1415)"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)

pic_name = "PtonshellW"
todraw = "w1_pt*(w1_mass>w2_mass)+w2_pt*(w1_mass<w2_mass)"
x_title = "p_{T} of onshell W"
x_bins = "(100,0,400)"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)

pic_name = "EtanuonshellW"
todraw = "nu1_eta*(w1_mass>w2_mass)+nu2_eta*(w1_mass<w2_mass)"
x_title = "#eta of neutrinos from onshell W"
x_bins = "(100,-4,4)"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)

pic_name = "PhinuonshellW"
todraw = "nu1_phi*(w1_mass>w2_mass)+nu2_phi*(w1_mass<w2_mass)"
x_title = "#phi of neutrinos from onshell W"
x_bins = "(100,-3.1415,3.1415)"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)

pic_name = "PtnuonshellW"
todraw = "nu1_pt*(w1_mass>w2_mass)+nu2_pt*(w1_mass<w2_mass)"
x_title = "p_{T} of neutrinos from onshell W"
x_bins = "(100,0,400)"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)

pic_name = "EtaoffshellW"
todraw = "w1_eta*(w1_mass<w2_mass)+w2_eta*(w1_mass>w2_mass)"
x_title = "#eta of offshell W"
x_bins = "(100,-4,4)"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)

pic_name = "PhioffshellW"
todraw = "w1_phi*(w1_mass<w2_mass)+w2_phi*(w1_mass>w2_mass)"
x_title = "#Phi of offshell W"
x_bins = "(100,-3.1415,3.1415)"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)

pic_name = "PtoffshellW"
todraw = "w1_pt*(w1_mass<w2_mass)+w2_pt*(w1_mass>w2_mass)"
x_title = "p_{T} of offshell W"
x_bins = "(100,0,400)"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)

pic_name = "EtanuoffshellW"
todraw = "nu1_eta*(w1_mass<w2_mass)+nu2_eta*(w1_mass>w2_mass)"
x_title = "#eta of neutrinos from offshell W"
x_bins = "(100,-4,4)"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)

pic_name = "PhinuoffshellW"
todraw = "nu1_phi*(w1_mass<w2_mass)+nu2_phi*(w1_mass>w2_mass)"
x_title = "#phi of neutrinos from offshell W"
x_bins = "(100,-3.1415,3.1415)"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)

pic_name = "PtnuoffshellW"
todraw = "nu1_pt*(w1_mass<w2_mass)+nu2_pt*(w1_mass>w2_mass)"
x_title = "p_{T} of neutrinos from offshell W"
x_bins = "(100,0,400)"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)

pic_name = "leadingbjetratio"
todraw = "(b1_pt/b1jet_pt)*(b1jet_pt>b2jet_pt)+(b2_pt/b2jet_pt)*(b1jet_pt<b2jet_pt)"
x_title = "#frac{p_{T}(b parton)}{p_{T}(bjet)}, leading bjet"
x_bins = "(600,0.0,3.0)"
cut = "h2tohh && hasRECOjet1 && hasRECOjet2 && dR_b1jet<0.4 && dR_b2jet<0.4 && hastwomuons"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)

pic_name = "MT2_reco"
todraw = "MT2_reco"
x_bins = "(80,0.0,400)"
x_title = "M_{T2} [GeV]"
cut = "(h2tohh || ttbar) && hasRECOjet1 && hasRECOjet2 && hastwomuons && dR_b1jet<0.4 && dR_b2jet<0.4 && hasdRljet"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)

pic_name = "dR_l1l2"
todraw = "dR_l1l2"
x_bins = "(50,0.0,5)"
x_title = "#Delta R(l,l)"
cut = "(h2tohh || ttbar) && hasRECOjet1 && hasRECOjet2 && hastwomuons && dR_b1jet<0.4 && dR_b2jet<0.4 && hasdRljet"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)
pic_name = "dR_l1l2gen"
todraw = "dR_genl1l2"
x_title = "#Delta R(l,l), gen level"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)

pic_name = "dR_b1b2"
todraw = "dR_b1b2"
x_bins = "(50,0.0,5)"
x_title = "#Delta R(j,j)"
cut = "(h2tohh || ttbar) && hasRECOjet1 && hasRECOjet2 && hastwomuons && dR_b1jet<0.4 && dR_b2jet<0.4 && hasdRljet"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)
pic_name = "dR_b1b2gen"
todraw = "dR_genb1b2"
x_title = "#Delta R(j,j), gen level"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)

pic_name = "dR_l1l2b1b2"
todraw = "dR_l1l2b1b2"
x_bins = "(50,0.0,5)"
x_title = "#Delta R(ll,jj)"
cut = "(h2tohh || ttbar) && hasRECOjet1 && hasRECOjet2 && hastwomuons && dR_b1jet<0.4 && dR_b2jet<0.4 && hasdRljet"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)
pic_name = "dR_l1l2b1b2gen"
todraw = "dR_genl1l2b1b2"
x_title = "#Delta R(ll,jj), gen level"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)

pic_name = "dphi_l1l2b1b2"
todraw = "dphi_l1l2b1b2"
x_bins = "(50,0.0,3.5)"
x_title = "#Delta #phi(ll,jj)"
cut = "(h2tohh || ttbar) && hasRECOjet1 && hasRECOjet2 && hastwomuons && dR_b1jet<0.4 && dR_b2jet<0.4 && hasdRljet"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)
pic_name = "dphi_l1l2b1b2gen"
todraw = "dphi_genl1l2b1b2"
x_title = "#Delta #phi(ll,jj), gen level"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)

pic_name = "dphi_llmet"
todraw = "dphi_llmet"
x_bins = "(50,0.0,3.5)"
x_title = "#Delta #phi(ll, #slash{E}_{T})"
cut = "(h2tohh || ttbar) && hasRECOjet1 && hasRECOjet2 && hastwomuons && dR_b1jet<0.4 && dR_b2jet<0.4 && hasdRljet"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)
pic_name = "dphi_llmetgen"
todraw = "dphi_genllmet"
x_title = "#Delta #phi(ll, #slash{E}_{T}), gen level"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)

pic_name = "dR_minbl"
todraw = "dR_minbl"
x_bins = "(50,0.0,5)"
x_title = "min #Delta R(l,j)"
cut = "(h2tohh || ttbar) && hasRECOjet1 && hasRECOjet2 && hastwomuons && dR_b1jet<0.4 && dR_b2jet<0.4 && hasdRljet"#hasdRljet
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)
pic_name = "dR_minblgen"
todraw = "dR_genminbl"
x_title = "min #Delta R(l,j), gen level"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)

pic_name = "mass_l1l2"
todraw = "mass_l1l2"
x_bins = "(50,.0,400)"
x_title = "M(ll) [GeV]"
cut = "(h2tohh || ttbar) && hasRECOjet1 && hasRECOjet2 && hastwomuons && dR_b1jet<0.4 && dR_b2jet<0.4 && hasdRljet"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)
pic_name = "mass_l1l2gen"
todraw = "mass_genl1l2"
x_title = "M(ll), gen level"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)

pic_name = "mass_b1b2"
todraw = "mass_b1b2"
x_bins = "(80,0.0,400)"
x_title = "M(jj) [GeV]"
cut = "(h2tohh || ttbar) && hasRECOjet1 && hasRECOjet2 && hastwomuons && dR_b1jet<0.4 && dR_b2jet<0.4 && hasdRljet"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)
pic_name = "mass_b1b2gen"
todraw = "mass_genb1b2"
x_title = "M(jj), gen level"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)

pic_name = "mass_trans"
todraw = "mass_trans"
x_bins = "(60,0.0,300)"
x_title = "M_{trans} [GeV]"
cut = "(h2tohh || ttbar) && hasRECOjet1 && hasRECOjet2 && hastwomuons && dR_b1jet<0.4 && dR_b2jet<0.4 && hasdRljet"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)
pic_name = "mass_transgen"
todraw = "mass_gentrans"
x_title = "M_{trans}, gen level"
#draw1D(filelist,"evtree",todraw,x_bins,x_title,cut,benchmarks, pic_name)
