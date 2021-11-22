import os
import ROOT
import re
import numpy as np
#import sys 
#sys.argv.append( '-b' )
#or ROOT.gROOT.SetBatch(1)

def rescalesignal(rootfile, rdir, prefix, factor):
    rfile = ROOT.TFile(rootfile, "UPDATE")
    #rdir = "mjj_vs_NN_M400"
    a = rfile.Get(rdir)
    #print "a ",a, " ",a.GetListOfKeys()
    keys = a.GetListOfKeys()
    #print "keys ",keys
    for i, key in enumerate(keys):
        #obj = ROOT.TH1F(key.ReadObj())
        obj = key.ReadObj()
        #print "obj ",obj
        if prefix in obj.GetName():
            obj.Scale(factor)
        #print "hist name ", obj.GetName()
    rfile.Write("",ROOT.TObject.kOverwrite)
    rfile.Close()

def rescalesignalall(masslist, workdir, factor):
    channels =   ["ElEl","MuEl","MuMu"]
    for mass in masslist:
       thisdir = workdir + "M%d.r7526/"%mass
       for ch in channels:
    	   rootfile = thisdir + "GGToX0ToHHTo2B2L2Nu_M%d_%s_shapes.root"%(mass, ch)
    	   newfile = thisdir + "GGToX0ToHHTo2B2L2Nu_M%d_%s_shapes_signalscale%d.root"%(mass, ch, int(factor))
    	   os.system("cp %s %s "%(rootfile, newfile))
    	   rdir = "mjj_vs_NN_M%d"%mass
    	   prefix = "ggX0HH%d"%mass
    	   rescalesignal(newfile, rdir, prefix, factor)
    

def extranumber(s):
    num = []
    for t in s.split():
        try:
            num.append(float(t))
        except ValueError:
            pass

    return num


def extractlimitfromtxtfile(logfile):
    logopen = open(logfile, "read")
    signal_xsec = 1.0
    #signal_xsec = 5000.0
    percents = [2.5, 16.0, 50.0, 84.0, 97.5, -1]
    limits_lines = {}
    limits = {}
    #for key in ["50.0%", " 2.5%", "16.0%", "84.0%", "97.5%", "Observed"]:
    for key in percents:
	#print(key, "Expected %4.1f%%:"%key)
        limits_lines[key] = []
    for line in logopen:
        if line.startswith("Expected "):
            for key in [2.5, 16.0, 50.0, 84.0, 97.5]:
                keystr = "Expected %4.1f"%key
                if keystr in line:
                    limits_lines[key].append(line)
        elif line.startswith("Observed Limit:"):
            limits_lines[-1].append(line)
    #print("limits_lines ", limits_lines)
    for key in percents:
        if len(limits_lines[key]) == 0:
            continue
        line = limits_lines[key][-1] 
        if key != -1:
            #print(key, "Expected %4.1f%%:"%key)
            line = line.replace("Expected %4.1f%%:"%key, "")
        nums = extranumber(line)
        if len(nums)>0:
            limits[key]  = nums[0] * signal_xsec
    #for line in logopen:
    #    #if line.startswith ("median expected limit: "):
    #	if line.startswith("Expected 50.0%:"):
    #        line = line.replace("Expected 50.0%:", "")
    #        nums = extranumber(line)
    #        limits[50.0]  = nums[0] * signal_xsec
    #    elif line.startswith("Expected  2.5%:"):
    #        line = line.replace("Expected  2.5%:", "")
    #        nums = extranumber(line)
    #        limits[2.50] = nums[0] * signal_xsec
    #    elif line.startswith("Expected 16.0%:"):
    #        line = line.replace("Expected 16.0%:", "")
    #        nums = extranumber(line)
    #        limits[16.0] = nums[0] * signal_xsec
    #    elif line.startswith("Expected 84.0%:"):
    #        line = line.replace("Expected 84.0%:", "")
    #        nums = extranumber(line)
    #        limits[84.0] = nums[0] * signal_xsec
    #    elif line.startswith("Expected 97.5%:"):
    #        line = line.replace("Expected 97.5%:", "")
    #        nums = extranumber(line)
    #        limits[97.5] = nums[0] * signal_xsec
    #    elif line.startswith("Observed Limit:"):
    #        nums = extranumber(line)
    #        limits[-1] = nums[0] * signal_xsec
    #    else :
    #        pass
    #          
    #print "limits ",limits
    logopen.close()
    return limits

def extractlimitfromtxtfile_t100(logfile):
    logopen = open(logfile, "read")
    signal_xsec = 1.0
    #signal_xsec = 5000.0
    percents = [2.5, 16.0, 50.0, 84.0, 97.5, -1]
    limits_lines = {}
    limits = {}
    observed_list = []
    ntoys = 0; observed_tot = 0.0
    for line in logopen:
        #if line.startswith ("median expected limit: "):
    	if line.startswith("median expected limit: "):
            #line = line.replace("Expected 50.0%:", "")
            nums = extranumber(line)
            limits[50.0]  = nums[0] * signal_xsec
        elif line.startswith("   68% expected band :"):
            line = line.replace("   68% expected band :", "")
            nums = extranumber(line)
            limits[16.0] = nums[0] * signal_xsec
            limits[84.0] = nums[1] * signal_xsec
        elif line.startswith("   95% expected band :"):
            line = line.replace("   95% expected band :", "")
            nums = extranumber(line)
            limits[2.5] = nums[0] * signal_xsec
            limits[97.5] = nums[1] * signal_xsec
        elif line.startswith("Observed Limit:"):
            nums = extranumber(line)
            observed_list.append(nums[0])
	    observed_tot = observed_tot+nums[0]
	    ntoys += 1
    from numpy import median
    limits[-1] = median(observed_list)
    print limits
    logopen.close()
    return limits



def CombineLimitplots(filelist, histlist, masspoints, xtitle, legends, text, plotname, drawData=True, drawUncertainty=True, runyear="2016"):
    #drawData = True
    #drawUncertainty = True
    colors = [ROOT.kRed, ROOT.kMagenta+2,ROOT.kBlue+1, ROOT.kOrange+2, ROOT.kAzure, ROOT.kCyan, ROOT.kViolet, ROOT.kBlack]
    markers = [21,22,23, 20, 24, 33, 34, 32, 29]
    tfilelist = []
    for f in filelist:
        tfilelist.append( ROOT.TFile(f, "READ"))


    c1 = ROOT.TCanvas("c1","c1",600, 800)
    c1.SetLogy()
    c1.SetGridx()  
    c1.SetGridy()  
    c1.SetTickx()  
    c1.SetTicky()  
    minx = min(masspoints)*0.8
    maxx = max(masspoints)*1.3
    b1 = ROOT.TH1F("b2","b2", len(masspoints)*2, minx, maxx)
    #b1.SetTitle(" #scale[1.4]{#font[61]{CMS}} Simulation Preliminary"+"  "*14+"35.87 fb^{-1} (13 TeV),2016")
    b1.SetTitle(" ")
    b1.GetYaxis().SetTitle("95% C.L. limits on production rate (fb)")
    b1.GetXaxis().SetTitle(xtitle)
    #yhigh = max(twosigma_up)*1.2
    #ylow = min(twosigma_low)*.9
    yhigh = 2000.0 #2000.0
    ylow = 1.0#1.0
    b1.GetYaxis().SetRangeUser(ylow, yhigh)
    b1.SetStats(0)
    b1.Draw()


    leg0 = ROOT.TLegend(0.65,0.8,0.9,0.90)
    leg0.SetFillColor(ROOT.kWhite)
    leg0.SetTextFont(42)
    leg0.SetHeader( legends[-1] )
    #le0g.AddEntry(grdata,"observed","pl")

    ## draw error band
    if drawUncertainty:
        i =  len(tfilelist)- 1
        tf = tfilelist[i]
        tf.cd()
        g_central = tf.Get(histlist[i]+"_central")
        g_central.SetLineColor(colors[i])
        g_central.SetMarkerColor(colors[i])
        g_central.SetMarkerStyle(markers[i])
        g_onesigma = tf.Get(histlist[i]+"_onesigma")
        g_twosigma = tf.Get(histlist[i]+"_twosigma")
        g_twosigma.SetFillColorAlpha(ROOT.kOrange, 0.5)
        g_twosigma.SetLineStyle(2)
        g_onesigma.SetFillColorAlpha(ROOT.kGreen+2, 0.5)
        g_onesigma.SetLineStyle(2)
        #leg0.AddEntry(g_data,"Observed","l")
        leg0.AddEntry(g_central,"Expected 95% upper limit","l")
        leg0.AddEntry(g_onesigma,"1 std. deviation, Expected","f")
        leg0.AddEntry(g_twosigma,"2 std. deviation, Expected","f")
        #if not drawData:
        #    g_central.Draw("lsame")
        g_twosigma.Draw("fe3same")
        g_onesigma.Draw("fe3same")


    leg = ROOT.TLegend(0.15,0.15,0.4,0.15+0.045*len(tfilelist))
    leg.SetFillColor(ROOT.kWhite)
    leg.SetTextFont(42)
    #if drawData:
    #    leg.SetHeader("Observed")
    #else:
    ## use asimov dataset, always expected 
    leg.SetHeader("Expected Limits")

    for i, tf in enumerate(tfilelist):
        tf.cd()
        g_data = tf.Get(histlist[i]+"_data")
        g_data.SetLineColor(colors[i])
        g_data.SetMarkerColor(colors[i])
        g_data.SetMarkerStyle(markers[i])
        g_central = tf.Get(histlist[i]+"_central")
        g_central.SetLineColor(colors[i])
        g_central.SetMarkerColor(colors[i])
        g_central.SetMarkerStyle(markers[i])
        if drawData:
            g_data.Draw("lpsame")
            g_central.SetLineStyle(2)
            g_central.Draw("lsame")
            thisleg = leg.AddEntry(g_data, legends[i],"p")
            thisleg.SetTextColor(colors[i])
        else:
            g_central.Draw("lpsame")
            thisleg = leg.AddEntry(g_central, legends[i],"p")
            thisleg.SetTextColor(colors[i])
    tex0 = ROOT.TLatex(0.08,0.91, "#scale[1.4]{#font[61]{CMS}} Internal"+" "*10+ runyear)
    tex0.SetNDC(); tex0.SetTextSize(.04); tex0.SetTextFont(42)
    tex0.Draw("same")
    tex1 = ROOT.TLatex(0.2,0.21, text)
    tex1.SetNDC(); tex1.SetTextSize(.04); tex1.SetTextFont(42)
    #tex1.Draw("same")
    
    leg0.Draw("same")
    leg.Draw("same")
    if drawData:
        c1.SaveAs(plotname+"_95Upperlmit_data_HHbbWW_combined.pdf")
        c1.SaveAs(plotname+"_95Upperlmit_data_HHbbWW_combined.C")
    else:
        c1.SaveAs(plotname+"_95Upperlmit_nodata_HHbbWW_combined.pdf")
        c1.SaveAs(plotname+"_95Upperlmit_nodata_HHbbWW_combined.C")

def makeBrazilPlot(masspoints_v0, alllimits, xtitle, text, plotname, drawData=False, runyear="2016"):
    onesigma_up = []
    twosigma_up = []
    onesigma_low = []
    twosigma_low = []
    central = []
    data = []
    masspoints = []
    for mass in masspoints_v0:
        if mass in alllimits.keys() and len(alllimits[mass].keys()) < 6:
            print("warning!!!!!!, not all limits found on mass %d"%mass, limits)
            continue
        limits = alllimits[mass]	
        central.append(limits[50.0])
        twosigma_low.append(limits[2.5])
        onesigma_low.append(limits[16.0])
        onesigma_up.append(limits[84.0])
        twosigma_up.append(limits[97.5])
    	data.append(limits[-1])
        masspoints.append(mass)
    fakeerrors = [0.0]*len(masspoints)
    #g_onesigma = TGraphAsymmErrors(len(masspoints),  np.array(masspoints),  np.array(central), np.array(fakeerrors), np.array(fakeerrors), np.array(onesigma_low), np.array(onesigma_up))
    #g_twosigma = TGraphAsymmErrors(len(masspoints),  np.array(masspoints),  np.array(central), np.array(fakeerrors), np.array(fakeerrors), np.array(twosigma_low), np.array(twosigma_up))
    c1 = ROOT.TCanvas("c1","c1",600, 800)
    outfilename = plotname+".root"
    #tfile = ROOT.TFile(outfilename,"UPDATE")
    tfile = ROOT.TFile(outfilename,"RECREATE")
    
    c1.SetLogy()
    c1.SetGridx()  
    c1.SetGridy()  
    c1.SetTickx()  
    c1.SetTicky()  
    onesigma_low.reverse()
    twosigma_low.reverse()
    onesigma_all = onesigma_up + onesigma_low
    twosigma_all = twosigma_up + twosigma_low
    masspoints_all = masspoints + list(reversed(masspoints))
    masspoints_f =  np.array(masspoints)+0.0
    #print "allXpoints ",masspoints_all," onesigma ",onesigma_all," float masspoints ",masspoints_f
    g_data = ROOT.TGraph(len(masspoints),  np.array(masspoints)+0.0,  np.array(data))
    g_central = ROOT.TGraph(len(masspoints),  np.array(masspoints)+0.0,  np.array(central))
    g_onesigma = ROOT.TGraph(len(masspoints)*2,  np.array(masspoints_all)+0.0,  np.array(onesigma_all))
    g_twosigma = ROOT.TGraph(len(masspoints)*2,  np.array(masspoints_all)+0.0,  np.array(twosigma_all))


    #g_twosigma.SetFillColor(ROOT.kYellow)
    g_twosigma.SetFillColorAlpha(ROOT.kOrange, 0.5)
    g_twosigma.SetLineStyle(2)
    g_onesigma.SetFillColorAlpha(ROOT.kGreen+2, 0.5)
    g_onesigma.SetLineStyle(2)
    g_central.SetLineWidth(2)
    g_central.SetLineStyle(7)
    g_central.SetLineColor(9)
    g_central.SetMarkerStyle(20)
    g_central.SetMarkerSize(1)
    g_central.SetMarkerColor(9)
    g_data.SetLineWidth(2)
    g_data.SetLineStyle(1)
    g_data.SetLineColor(ROOT.kBlack)
    g_data.SetMarkerStyle(20)
    g_data.SetMarkerSize(1)
    g_data.SetMarkerColor(ROOT.kBlack)

    #b1 = ROOT.TH1F("b2","b2",14, 250.0, 950.0)
    minx = min(masspoints)*0.9
    maxx = max(masspoints)*1.3
    b1 = ROOT.TH1F("b2","b2", len(masspoints)*2, minx, maxx)
    #b1.SetTitle(" #scale[1.4]{#font[61]{CMS}} Simulation Preliminary"+"  "*14+"35.87 fb^{-1} (13 TeV),2016")
    b1.SetTitle(" ")
    b1.GetYaxis().SetTitle("95% C.L. limits on production rate (fb)")
    b1.GetXaxis().SetTitle(xtitle)
    #yhigh = max(twosigma_up)*1.2
    #ylow = min(twosigma_low)*.9
    yhigh = 2000.0 #10000.0
    ylow = 0.01#10.0
    b1.GetYaxis().SetRangeUser(ylow, yhigh)
    b1.SetStats(0)

    b1.Draw()
    g_twosigma.Draw("fe3same")
    g_onesigma.Draw("fe3same")
    g_central.Draw("lpsame")
    if drawData :
        g_data.Draw("lpsame")
    suffixname =  plotname.split("/")[-1]
    g_data.SetName("%s_data"%suffixname)
    g_central.SetName("%s_central"%suffixname)
    g_onesigma.SetName("%s_onesigma"%suffixname)
    g_twosigma.SetName("%s_twosigma"%suffixname)
    g_central.Write()
    g_data.Write()
    g_onesigma.Write()
    g_twosigma.Write()

    
    leg = ROOT.TLegend(0.65,0.8,0.9,0.90)
    leg.SetFillColor(ROOT.kWhite)
    leg.SetTextFont(42)
    if drawData :
        leg.AddEntry(g_data,"Observed","pl")
    leg.AddEntry(g_central,"Expected 95% upper limit","l")
    leg.AddEntry(g_onesigma,"1 std. deviation, Expected","f")
    leg.AddEntry(g_twosigma,"2 std. deviation, Expected","f")
    tex0 = ROOT.TLatex(0.08,0.91, "#scale[1.4]{#font[61]{CMS}} Internal"+" "*10+runyear)
    tex0.SetNDC(); tex0.SetTextSize(.04); tex0.SetTextFont(42)
    tex0.Draw("same")
    tex1 = ROOT.TLatex(0.13,0.21, text)
    tex1.SetNDC(); tex1.SetTextSize(.04); tex1.SetTextFont(42)
    tex1.Draw("same")
    
    leg.Draw("same")
    if drawData:
        c1.SaveAs(plotname+"_95Upperlmit_data_HHbbWW.pdf")
        c1.SaveAs(plotname+"_95Upperlmit_data_HHbbWW.C")
    else:
        c1.SaveAs(plotname+"_95Upperlmit_expected_HHbbWW.pdf")
        c1.SaveAs(plotname+"_95Upperlmit_expected_HHbbWW.C")

    tfile.Close()

def makeComparePlot(masspoints_v0, key_to_limits, reference, xtitle, text, plotname, drawData, year):
    colors = [ROOT.kBlack, ROOT.kRed, ROOT.kMagenta+2,ROOT.kBlue+1, ROOT.kOrange+2, ROOT.kAzure, ROOT.kCyan, ROOT.kViolet]
    markers = [20, 21,22,23, 24, 33, 34, 32, 29]
    onesigma_up = []
    twosigma_up = []
    onesigma_low = []
    twosigma_low = []
    central = {}
    data = {}
    key_masspoints = {}
    for key in key_to_limits.keys(): 
        central[key] = []
        data[key] = []
        key_masspoints[key] =[]
    #key0 = key_to_limits.keys()[0]
    #for mass in masspoints_v0:
    #    limits = key_to_limits[key0][mass]	
    #    if len(limits.keys()) < 6:
    #        print("warning!!!!!!, not all limits found on mass %d"%mass, limits)
    #        continue
    #    twosigma_low.append(limits[2.5])
    #    onesigma_low.append(limits[16.0])
    #    onesigma_up.append(limits[84.0])
    #    twosigma_up.append(limits[97.5])
    #    for key in key_to_limits.keys(): 
    #        #print("mass ", mass, " key ", key, " median limit ", key_to_limits[key][mass][50.0])
    #        central[key].append(key_to_limits[key][mass][50.0])
    #        data[key].append(key_to_limits[key][mass][-1])
    #    masspoints.append(mass)
    for key in key_to_limits.keys():
        for mass in masspoints_v0:
            if mass in key_to_limits[key].keys() and len(key_to_limits[key][mass].keys()) == 6:
                limits = key_to_limits[key][mass]
                central[key].append(limits[50.0])
                data[key].append(limits[-1])
                key_masspoints[key].append(mass)
                if key == reference:
                    twosigma_low.append(limits[2.5])
                    onesigma_low.append(limits[16.0])
                    onesigma_up.append(limits[84.0])
                    twosigma_up.append(limits[97.5])
            elif mass in key_to_limits[key].keys():
                print("warning!!!!!!, not all limits found on mass %d for %s"%(mass, key), key_to_limits[key][mass])
                continue
            else:
                print("warning!!!!!!, benchmark mass %d is not found for %s "%(mass, key))
                continue
            
    ## first set of limits. also ploting the error band usually
    fakeerrors = [0.0]*len(key_masspoints[reference])
    masspoints = key_masspoints[reference]
    #g_onesigma = TGraphAsymmErrors(len(masspoints),  np.array(masspoints),  np.array(central), np.array(fakeerrors), np.array(fakeerrors), np.array(onesigma_low), np.array(onesigma_up))
    #g_twosigma = TGraphAsymmErrors(len(masspoints),  np.array(masspoints),  np.array(central), np.array(fakeerrors), np.array(fakeerrors), np.array(twosigma_low), np.array(twosigma_up))
    c1 = ROOT.TCanvas("c1","c1",600, 800)
    outfilename = plotname+".root"
    tfile = ROOT.TFile(outfilename,"RECREATE")
    
    c1.SetLogy()
    c1.SetGridx()  
    c1.SetGridy()  
    c1.SetTickx()  
    c1.SetTicky()  
    onesigma_low.reverse()
    twosigma_low.reverse()
    onesigma_all = onesigma_up + onesigma_low
    twosigma_all = twosigma_up + twosigma_low
    masspoints_all = masspoints + list(reversed(masspoints))
    masspoints_f =  np.array(masspoints)+0.0
    #print "allXpoints ",masspoints_all," onesigma ",onesigma_all," float masspoints ",masspoints_f
    g_data = ROOT.TGraph(len(masspoints),  np.array(masspoints)+0.0,  np.array(data[reference]))
    g_central = ROOT.TGraph(len(masspoints),  np.array(masspoints)+0.0,  np.array(central[reference]))
    g_onesigma = ROOT.TGraph(len(masspoints)*2,  np.array(masspoints_all)+0.0,  np.array(onesigma_all))
    g_twosigma = ROOT.TGraph(len(masspoints)*2,  np.array(masspoints_all)+0.0,  np.array(twosigma_all))


    #g_twosigma.SetFillColor(ROOT.kYellow)
    g_twosigma.SetFillColorAlpha(ROOT.kOrange, 0.5)
    g_twosigma.SetLineStyle(2)
    g_onesigma.SetFillColorAlpha(ROOT.kGreen+2, 0.5)
    g_onesigma.SetLineStyle(2)

    g_central.SetLineWidth(2)
    g_central.SetLineStyle(7)
    g_central.SetLineColor(ROOT.kBlack)
    g_central.SetMarkerStyle(20)
    g_central.SetMarkerSize(1)
    g_central.SetMarkerColor(ROOT.kBlack)
    g_data.SetLineWidth(2)
    g_data.SetLineStyle(1)
    g_data.SetLineColor(ROOT.kBlack)
    g_data.SetMarkerStyle(20)
    g_data.SetMarkerSize(1)
    g_data.SetMarkerColor(ROOT.kBlack)

    ### compare limits
    #leg2 = ROOT.TLegend(0.15,0.15,0.4,0.15+0.045*len(key_to_limits.keys()))
    leg2 = ROOT.TLegend(0.65,0.76-0.04*len(key_to_limits.keys()), 0.9, 0.76)
    leg2.SetFillColor(ROOT.kWhite)
    leg2.SetTextFont(42)
    leg2.AddEntry(g_central, reference, "lp")
    #if drawData:
    i = 1
    g_data_rest = {}
    g_central_rest = {}
    for key in key_to_limits.keys(): 
        if key == reference:
            continue
        if drawData :
            g_data_rest[key] = ROOT.TGraph(len(key_masspoints[key]),  np.array(key_masspoints[key])+0.0,  np.array(data[key]))
            g_data_rest[key].SetLineWidth(2)
            g_data_rest[key].SetLineStyle(1)
            g_data_rest[key].SetLineColor(colors[i])
            g_data_rest[key].SetMarkerStyle(markers[i])
            g_data_rest[key].SetMarkerSize(1)
            g_data_rest[key].SetMarkerColor(colors[i])

        g_central_rest[key]  = ROOT.TGraph(len(key_masspoints[key]),  np.array(key_masspoints[key])+0.0,  np.array(central[key]))
        g_central_rest[key].SetLineWidth(2)
        g_central_rest[key].SetLineStyle(7)
        g_central_rest[key].SetLineColor(colors[i])
        g_central_rest[key].SetMarkerStyle(markers[i])
        g_central_rest[key].SetMarkerSize(1)
        g_central_rest[key].SetMarkerColor(colors[i])
        leg2.AddEntry(g_central_rest[key], key, "lp")
        i += 1

    #b1 = ROOT.TH1F("b2","b2",14, 250.0, 950.0)
    minx = min(masspoints)*0.9
    maxx = max(masspoints)*1.3
    b1 = ROOT.TH1F("b2","b2", len(masspoints)*2, minx, maxx)
    #b1.SetTitle(" #scale[1.4]{#font[61]{CMS}} Simulation Preliminary"+"  "*14+"35.87 fb^{-1} (13 TeV),2016")
    b1.SetTitle(" ")
    b1.GetYaxis().SetTitle("95% C.L. limits on production rate (fb)")
    b1.GetXaxis().SetTitle(xtitle)
    #yhigh = max(twosigma_up)*1.2
    #ylow = min(twosigma_low)*.9
    yhigh = 2000.0 #10000.0
    ylow = 0.01#10.0
    b1.GetYaxis().SetRangeUser(ylow, yhigh)
    b1.SetStats(0)

    b1.Draw()
    g_twosigma.Draw("fe3same")
    g_onesigma.Draw("fe3same")
    g_central.Draw("lpsame")
    if drawData :
        g_data.Draw("lpsame")
    suffixname =  reference
    g_data.SetName("%s_data"%suffixname)
    g_central.SetName("%s_central"%suffixname)
    g_onesigma.SetName("%s_onesigma"%suffixname)
    g_twosigma.SetName("%s_twosigma"%suffixname)
    g_central.Write()
    g_data.Write()
    g_onesigma.Write()
    g_twosigma.Write()
    
    leg = ROOT.TLegend(0.65,0.8,0.9,0.90)
    leg.SetFillColor(ROOT.kWhite)
    leg.SetTextFont(42)
    if drawData :
        leg.AddEntry(g_data,"Observed","pl")
    leg.AddEntry(g_central,"Expected 95% upper limit","l")
    leg.AddEntry(g_onesigma,"1 std. deviation, Expected","f")
    leg.AddEntry(g_twosigma,"2 std. deviation, Expected","f")
    leg.Draw("same")

    for key in g_central_rest.keys(): 
        g_central_rest[key].SetName(key+"_central")
        g_central_rest[key].Draw("samelp")
        g_central_rest[key].Write()
        if drawData :
            g_data_rest[key].Draw("samelp")
            g_data_rest[key].SetName(key+"_data")
            g_data_rest[key].Write()
    leg2.Draw("same")

    tex0 = ROOT.TLatex(0.08,0.91, "#scale[1.4]{#font[61]{CMS}} Internal"+" "*20+year)
    tex0.SetNDC(); tex0.SetTextSize(.04); tex0.SetTextFont(42)
    tex0.Draw("same")
    tex1 = ROOT.TLatex(0.13,0.21, text)
    tex1.SetNDC(); tex1.SetTextSize(.04); tex1.SetTextFont(42)
    tex1.Draw("same")
    

    if drawData:
        c1.SaveAs(plotname+"_95Upperlmit_data_HHbbWW.pdf")
        c1.SaveAs(plotname+"_95Upperlmit_data_HHbbWW.C")
    else:
        c1.SaveAs(plotname+"_95Upperlmit_expected_HHbbWW.pdf")
        c1.SaveAs(plotname+"_95Upperlmit_expected_HHbbWW.C")
    tfile.Close()


if __name__ == '__main__':
   print("helper.py")
