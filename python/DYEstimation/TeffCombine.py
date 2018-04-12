import ROOT
import os

def addHistos(lst):
    if len(lst) == 1:
        return
    for h in lst[1:]:
        lst[0].Add(h)

def combineTeffs(files, Teffnames, outname, outfile):
    #print "files ",files, " teffname ",Teffnames
    tfiles = []
    for f in files:
        tfiles.append(ROOT.TFile(f, "READ"))

    totalHists = []
    passHists = []
    tefflist = []
    #den = 0.0; num = 0.0
    #den2 = 0.0; num2 = 0.0
    for i, tf in enumerate(tfiles):
        teffname = Teffnames[i]
        tf.cd()
        teff = tf.Get(teffname)
        totHist = teff.GetCopyTotalHisto()
        passHist = teff.GetCopyPassedHisto()
        totHist.SetDirectory(0)
        passHist.SetDirectory(0)
        teff.SetDirectory(0)
        #den =totHist.GetBinContent(11, 5)
        #num =passHist.GetBinContent(11, 5)
        #den2= totHist.GetBinContent(11, 3)
        #num2= passHist.GetBinContent(11, 3)
        totalHists.append(totHist)
        passHists.append(passHist)
        tefflist.append(teff)
        #print "totHist ",totHist.Integral()," passHist ",passHist.Integral()
        #print " den ",den ," num ",num," (11, 3) den ",den2," num ",num2



    outtf = ROOT.TFile.Open(outfile, "update")
    addHistos(totalHists)
    addHistos(passHists)
    addHistos(tefflist)
    totalHists[0].SetTitle(outname)
    passHists[0].SetTitle(outname)
    #for x in range(1, totalHists[0].GetNbinsX() + 1):
    #    for y in range(1, totalHists[0].GetNbinsY() + 1):
    #        den = totalHists[0].GetBinContent(x, y)
    #        num = passHists[0].GetBinContent(x, y)
    #        print " x ",x," y ",y," den ",den, " num ",num
    #        if den < num :
    #            print "error !!!! "
    #print "totalHists[0] ",totalHists[0].Integral(), " passHists[0] ",passHists[0].Integral()
    thisEff = tefflist[0]
    #thisEff = ROOT.TEfficiency(passHists[0], totalHists[0])
    #thisEff = passHists[0].Clone()
    #thisEff.Divide(totalHists[0])
    thisEff.SetName(outname)
    thisEff.SetStatisticOption(ROOT.TEfficiency.kBUniform)

    thisEff.Write()
    outtf.Close()

def combineAllFractions(inputdir, todo_datasets, outfile):


    files = []
    for job in todo_datasets:
        files.append(os.path.join(inputdir, "dy_flavour_"+job+".root"))

    allTeffnames = []
    flavours = ["b","c","l"]
    for f1 in flavours:
        for f2 in flavours:
            allTeffnames.append(f1+f2+"_frac")
    temptf = ROOT.TFile(outfile, "recreate")
    temptf.Close()
    for teff in allTeffnames:
        tefflist = [teff]*len(files)
        outname = teff
        combineTeffs(files, tefflist, outname, outfile)


def combineBtaggingEff(inputdir, todo_datasets, outfile):
    files = []
    for job in todo_datasets:
        files.append(os.path.join(inputdir, "btagging_efficiency_"+job+".root"))

    
    allTeffnames = []
    flavours = ["b","c","light"]
    for f in flavours:
        allTeffnames.append("btagging_eff_on_%s_vs_pt"%f)
        allTeffnames.append("btagging_eff_on_%s_vs_eta"%f)
    for f in ['b','c']:
        allTeffnames.append("btagging_eff_on_"+f)
    allTeffnames.append("mistagging_eff_on_light")
    temptf = ROOT.TFile(outfile, "recreate")
    temptf.Close()
    for teff in allTeffnames:
        tefflist = [teff]*len(files)
        outname = teff
        combineTeffs(files, tefflist, outname, outfile)


todo_datasets = []
todo_datasets.append('DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8')
todo_datasets.append('DYToLL_0J_13TeV-amcatnloFXFX-pythia8')
todo_datasets.append('DYToLL_1J_13TeV-amcatnloFXFX-pythia8')
todo_datasets.append('DYToLL_2J_13TeV-amcatnloFXFX-pythia8')
inputdir = "BtaggingEff_fraction_cmstca/"
outfile_frac = "dy_flavour_frac_combined.root"
#combineAllFractions(inputdir, todo_datasets, outfile_frac)
outfile_btag = "dy_flavour_btagging_combined.root"
combineBtaggingEff(inputdir, todo_datasets, outfile_btag)
