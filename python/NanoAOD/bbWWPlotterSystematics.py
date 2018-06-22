import os
import sys
import ROOT
from math import fabs, sqrt

## ------------------------------ ##
## Setup Matplotlib Configuration ##



class bbWWPlotterSystematics(object):
    """Class for handling systematic uncertainties to put in plots"""
    def __init__(self,  inputfiledict, nominalWeight, systematiclist, DYdatadriven):
        """Initialize class for each sample you want systematics"""
        self.channels = ["MuMu", "MuEl","ElEl"]
        self.channelcuts = {"MuMu":"isMuMu", "MuEl":"(isMuEl || isElMu)","ElEl":"isElEl"}
        #self.outfile = outfile
        #self.addStatistics = True
        #self.writeSystematichists = False
        self.SysDict = {
                "CMS_eff_b_heavy":{  "weight":"jjbtag_heavy_nominal",
                                     "description":"b-tagging uncertainty from b,c-jets", 
                                     "event_weight_sum": "nominal",
                                     },
                "CMS_eff_b_light":{  "weight":"jjbtag_light_nominal",
                                     "description":"b-tagging uncertainty from light-jets",
                                     "event_weight_sum": "nominal",
                                     },
                "CMS_pu":{ "weight":"event_pu_weight_nominal",
                           "description":"pielup uncertainty",
                           "event_weight_sum":"nominal",
                           },
                "CMS_pdf":{"weight":"event_pdf_weight_nominal",
                            "description":"PDF uncertainty",
                            "event_weight_sum":"nominal",
                            },
                "CMS_eff_trigger":{"weight":"lep1trgsf_nominal * lep2trgsf_nominal", 
                                   "description": "trigger eff uncertainty",
                                   "event_weight_sum":"nominal",
                                   },
                #"CMS_eff_e":{"weight": "(isMuEl*lep2trackingsf_nominal*lep2HLTsafeIDsf_nominal*lep2IDsf_nominal + isElMu*lep1trackingsf_nominal*lep1HLTsafeIDsf_nominal*lep1IDsf_nominal+ isElEl*lep1trackingsf_nominal*lep1HLTsafeIDsf_nominal*lep1IDsf_nominal*lep2trackingsf_nominal*lep2HLTsafeIDsf_nominal*lep2IDsf_nominal+isMuMu)", 
                "CMS_eff_e":{"weight": "(isMuEl*lep2trackingsf_nominal*lep2IDsf_nominal + isElMu*lep1trackingsf_nominal*lep1IDsf_nominal+ isElEl*lep1trackingsf_nominal*lep1IDsf_nominal*lep2trackingsf_nominal*lep2IDsf_nominal+isMuMu)", 
                            "description":"Electron eff",
                            "event_weight_sum":"nominal",
                            },
                "CMS_eff_mu":{"weight": "(isElMu*lep2trackingsf_nominal*lep2IDsf_nominal + isMuEl*lep1trackingsf_nominal*lep1IDsf_nominal+ isMuMu*lep1trackingsf_nominal*lep1IDsf_nominal*lep2trackingsf_nominal*lep2IDsf_nominal+isElEl)",
                            "description":"Muon eff",
                            },
                "CMS_iso_mu":{"weight": "(lep1Isosf_nominal * lep2Isosf_nominal)", 
                              "description": "Muon Isolation",
                              "event_weight_sum":"nominal",
                              },
                "QCDscale":{"weight":None,
                            "description":"QCD scale uncertainty",
                            "event_weight_sum":None,
                            },
                }

        
        self.filedict = inputfiledict

        self.nominalWeight = nominalWeight

        self.systematiclist = systematiclist

        self.DYdatadriven = DYdatadriven
        self.shortnames_datadriven = ["TT_untagged", "sT_untagged"]

        self.treename = "evtreeHME_nn"
        #self.treename = "Friends"

    #def initialize1D(self, inputfiledict, nominalWeight, systematiclist, todraw, xtitle, xbins, cuts, DYdatadriven, untagged_filedict):
    def initialize1D(self, todraw, xtitle, xbins, cuts):
        """Setup some initial variables"""
        self.todraw = todraw
        self.xtitle = xtitle
        self.xbins = xbins
        self.cuts = cuts
        #self.shortname = shortname 

        self.sample_systematic_hist = {}
        self.systematic_hist = {}
        self.finalhist = {}
        self.channel_shortname_systematic_hist = {}
        for channel in self.channels:
            self.channel_shortname_systematic_hist[channel] = {}

    def get_xsection_eventweightsum_file(self, filename):
        tfile = ROOT.TFile(filename, "READ")
        xsec = tfile.Get("cross_section")
        event_weight_sum = tfile.Get("event_weight_sum")
        return xsec.GetVal(),event_weight_sum.GetVal()


    def GetQCDScaleHist(self, filename):
        tfile = ROOT.TFile(filename, "READ")
        hist = tfile.Get("CountWeightedLHEWeightScale")
        finalhist = hist.Clone()
        finalhist.SetDirectory(0)
        return finalhist


    #def QCDScaleSystematic(self, chain, weight, nominal_event_weight_sum, QCDscalehist, plotname):
    def QCDScaleSystematic(self, chain, weight, plotname):

        
        #####
        ####LHE scale variation weights (w_var / w_nominal); [0] is renscfact=0.5d0 facscfact=0.5d0 ; [1] is renscfact=0.5d0 facscfact=1d0 ; [2] is renscfact=0.5d0 facscfact=2d0 ; [3] is renscfact=1d0 facscfact=0.5d0 ; [4] is renscfact=1d0 facscfact=1d0 ; [5] is renscfact=1d0 facscfact=2d0 ; [6] is renscfact=2d0 facscfact=0.5d0 ; [7] is renscfact=2d0 facscfact=1d0 ; [8] is renscfact=2d0 facscfact=2d0
        ####
        indexlist = [0, 1, 3, 5, 7, 8]
        QCDscale_allshapes = []
        for index in indexlist:
            #this_event_weight_sum = QCDscalehist.GetBinContent(index + 1)
            #finalcut = "LHEScaleWeight[%d]"%index +"*"+weight 
            finalcut = "LHEScaleWeight_%d"%index +"*"+weight 
            hist = ROOT.TH1F("QCDScale_shape%d"%index, "QCDScale_shape%d"%index, len(self.xbins)-1, self.xbins)
            QCDscale_allshapes.append(hist)
            chain.Draw(self.todraw + ">> QCDScale_shape%d"%index, finalcut)
            #print "finalcut ",finalcut, " hist integral ", hist.Integral()

        hist_up = ROOT.TH1F("QCDscaleup","QCDscaleup", len(self.xbins)-1, self.xbins)
        hist_down = ROOT.TH1F("QCDscaledown", "QCDscaledown", len(self.xbins)-1, self.xbins)
        for bin in xrange( hist_up.GetNbinsX()):
            binvalues = []
            for shape in QCDscale_allshapes:
                binvalues.append(shape.GetBinContent(bin + 1))
            hist_up.SetBinContent(bin+1, max(binvalues))
            hist_down.SetBinContent(bin+1, min(binvalues))
            #print "bin ",bin," binvalues ",binvalues, " max ", max(binvalues), " min ",min(binvalues)


        #print "final integral up ", hist_up.Integral(), " down ", hist_down.Integral(), " nominal ",QCDscale_allshapes[3].Integral()
        plotQCDScale = False
        if plotQCDScale:
            samplename = plotname.split('/')[-1]
            hs = ROOT.THStack(samplename, "QCD scale uncertainty, intermedian plot")
            #hist_up.SetMarkerColor(2)
            #hist_up.SetMarkerStyle(22)
            #hist_down.SetMarkerColor(4)
            #hist_down.SetMarkerStyle(23)
            hist_up.SetLineColor(1)
            #hist_up.SetLineStyle(2)
            hist_up.SetLineWidth(2)
            hist_down.SetLineColor(1)
            #hist_down.SetLineStyle(2)
            hist_down.SetLineWidth(2)
            colors = [800-4, 820-3, 900-3, 860-3, 616-7, 432+2, 400+2]
            legend = ROOT.TLegend(0.74,0.5,0.84,0.5+7*.05); 
            legend.SetTextSize(0.04); legend.SetTextFont(42)
            for i, shape in enumerate(QCDscale_allshapes):
                shape.SetLineColor(colors[i])
                hs.Add(shape)
                legend.AddEntry(shape, "QCDscale_%d"%indexlist[i],"l")
            legend.SetBorderSize(0)
            legend.AddEntry(hist_up, "QCDscale Up","l")
            legend.AddEntry(hist_down, "QCDscale down","l")
     
            scale_c = ROOT.TCanvas("scale", "scale",800, 600)
            hs.Draw("hist nostack")
            hist_up.Draw("histsame")
            hist_down.Draw("histsame")
            legend.Draw("same")
            hs.GetHistogram().GetXaxis().SetTitle(self.xtitle)
            hs.GetHistogram().GetYaxis().SetTitle("Events")
            tex1 = ROOT.TLatex(0.13,0.87, samplename)
            tex1.SetNDC(); tex1.SetTextSize(.035)
            tex1.Draw("same")
            #plotdir = "DataDriven_DY_plots/"
            scale_c.SaveAs(plotname+"_QCDscale.pdf")
          


        hist_up.SetDirectory(0)
        hist_down.SetDirectory(0)
        return hist_up,hist_down



    #################################################################################################
    #### do systematics for one file
    #################################################################################################
    def plotSystematics_singleFile(self, key, channel, samplename, filename, nominalWeight): 

        
        self.sample_systematic_hist[samplename] = {}

        #xsec,event_weight_sum = self.get_xsection_eventweightsum_file(filename)
        #weight = nominalWeight+"*{cross_section}*1000.0/{event_weight_sum}".format(cross_section = xsec, event_weight_sum = event_weight_sum)
        weight = nominalWeight+"*cross_section*1000.0/event_weight_sum"
        if "Radion" in key or "Graviton" in key:
            xsec = 5.0
        chain = ROOT.TChain(self.treename)
        chain.Add( filename )
        self.sample_systematic_hist[samplename]["nominal"] = ROOT.TH1F(key+"_"+channel+"_%s"%samplename+"_nominal", key+"_"+channel+"_%s"%samplename+"_nominal", len(self.xbins)-1, self.xbins)
        finalcut = "("+ self.cuts +"&&"+ self.channelcuts[channel] +")*"+weight
        #print "sample ",samplename, " xsec ",xsec, "  event_weight_sum ",event_weight_sum, " weight ",weight , " finalcut ",finalcut
        chain.Draw(self.todraw + ">> " + self.sample_systematic_hist[samplename]["nominal"].GetName(), finalcut)
        #print "nominal integral ", self.sample_systematic_hist[samplename]["nominal"].Integral()
        self.sample_systematic_hist[samplename]["nominal"].SetDirectory(0)


        ###now plotting systematics
        for sys in self.systematiclist:
            self.sample_systematic_hist[samplename][sys] = {}

            if sys == "QCDscale":
                faileddatasets =  ["WWToLNuQQ_aTGC_13TeV-madgraph-pythia8", "ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1"]
                faileddatasets.append("ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1")
                if (samplename in faileddatasets) or (self.DYdatadriven and samplename.replace("_untagged","") in faileddatasets):
                    self.sample_systematic_hist[samplename][sys]["up"] = None
                    self.sample_systematic_hist[samplename][sys]["down"] = None
                    continue
                plotdir = "QCDScale_intermedianplots/"
                #QCDscalehist = self.GetQCDScaleHist(filename)
                plotname = os.path.join(plotdir, samplename+"_"+channel)
                weight_cuts = "("+ self.cuts +"&&"+ self.channelcuts[channel] +")*"+weight
                #hist_up, hist_down = self.QCDScaleSystematic(chain, weight_cuts, event_weight_sum, QCDscalehist, plotname)
                hist_up, hist_down = self.QCDScaleSystematic(chain, weight_cuts, plotname)
                hist_up.SetName(key+"_"+channel+"_%s"%samplename+"_QCDscaleup")
                hist_down.SetName(key+"_"+channel+"_%s"%samplename+"_QCDscaledown")
                #print "samplename ",samplename, " up ",hist_up.Print("ALL")," down ",hist_down.Print("ALL")," nominal ",self.sample_systematic_hist[samplename]["nominal"].Print("ALL")
                self.sample_systematic_hist[samplename][sys]["up"] = hist_up
                self.sample_systematic_hist[samplename][sys]["down"] = hist_down
                continue

            for plottype in ["up","down"]:

                suffix = sys+"_"+plottype
                hist = ROOT.TH1F(key+"_"+channel+"_%s"%samplename+"_%s"%(suffix), key+"_"+channel+"_%s"%samplename+"_%s"%(suffix), len(self.xbins)-1, self.xbins)
                self.sample_systematic_hist[samplename][sys][plottype] = hist
                
                sysweight = self.SysDict[sys]["weight"].replace("nominal", plottype)
                thisweight = sysweight+"*"+weight
                finalcut = "(" + self.cuts +"&&"+ self.channelcuts[channel] + ")*"+thisweight
                #chlist[samplename].Draw(todraw + ">> " + self.sample_systematic_hist[samplename][sys][plottype].GetName(), finalcut)
                #print "sample ",samplename," file ",self.filedict[key][samplename]['path']," weight ",thisweight," todraw ",self.todraw," finalcut ",finalcut
                chain.Draw(self.todraw + ">> " + self.sample_systematic_hist[samplename][sys][plottype].GetName(), finalcut)
                self.sample_systematic_hist[samplename][sys][plottype].SetDirectory(0)
                

    #################################################################################################
    ##plot systematcs of all sub-processes in one process like TT, sT, DY...
    #################################################################################################
    def plotSystematics(self, key, channel):
        
        self.shortname = key
        self.hist_data_untagged = None
        self.sample_systematic_hist = {}
       
        #print "plotSystematics ",key, " channel ", channel
        if self.shortname == "DY" and self.DYdatadriven and (channel == "MuMu" or channel == "ElEl"):
            ##plot data
            self.hist_data_untagged = ROOT.TH1F(self.shortname+"_"+channel+"_data_untagged", self.shortname+"_"+channel+"_untagged", len(self.xbins)-1, self.xbins)
            Mbtag_weight = "dy_Mbtag_weight"
            untagged_suffix = "_untagged"
            data_untagged_name = "DoubleMuon"
            if channel == "ElEl":
                data_untagged_name = "DoubleEG"
            datafile_untagged =  self.filedict["Data"+untagged_suffix][data_untagged_name]['path']
            ch_d = ROOT.TChain(self.treename)
            ch_d.AddFile(datafile_untagged)
            cut_data = "("+ self.cuts +"&&"+ self.channelcuts[channel] +")*"+Mbtag_weight
            ch_d.Draw(self.todraw + ">> " + self.hist_data_untagged.GetName(), cut_data)
            self.hist_data_untagged.SetDirectory(0)

            weight_datadriven = self.nominalWeight+"* "+ Mbtag_weight
            for key_datadriven in self.shortnames_datadriven:
                for samplename_datadriven in self.filedict[key_datadriven].keys():
                    self.plotSystematics_singleFile(key, channel, samplename_datadriven, self.filedict[key_datadriven][samplename_datadriven]['path'], weight_datadriven)
        else:
            for iname, samplename in enumerate(self.filedict[self.shortname].keys()):
                self.plotSystematics_singleFile(self.shortname, channel, samplename, self.filedict[self.shortname][samplename]['path'], self.nominalWeight)

    
    #################################################################################################
    ##combine subprocess systematcs for one process: like TT, sT, ttV, TT_untagged, sT_untagged
    #################################################################################################
    def combineSystematics(self, shortname, channel):
        ##combine to final plots
        self.shortname = shortname
        suffix = "nominal"
        #print "combineSystematics ", shortname, " channel ", channel
        self.finalhist["nominal"] = ROOT.TH1F(self.shortname+"_"+channel, self.shortname+"_"+channel, len(self.xbins)-1, self.xbins)
        allsamples_combine = None
        if "untagged" in self.shortname and self.DYdatadriven and (channel == "MuMu" or channel == "ElEl"):
            allsamples_combine = self.filedict[shortname]
        else:
            allsamples_combine = self.filedict[self.shortname].keys()

        for iname, samplename in enumerate(allsamples_combine):
            self.finalhist["nominal"].Add( self.sample_systematic_hist[samplename]["nominal"] )
        self.finalhist["nominal"].SetDirectory(0)

        dummyplots = {}
        totalSys = {}
        for sys in self.systematiclist:
            self.systematic_hist[sys] = {}
            totalSys[sys] = {}
        for plottype in ["up","down"]:

            dummyplots[plottype] = self.cloneDummyHistogram( self.finalhist["nominal"] )
            dummyplots[plottype].SetName("dummy"+plottype)
            for sys in self.systematiclist:
                totalSys[sys][plottype] = 0.0
                suffix = "_"+sys+"_"+plottype
                #hist = ROOT.TH1F(self.shortname+"_"+channel+suffix, self.shortname+"_"+channel+suffix, len(self.xbins)-1, self.xbins)
                self.systematic_hist[sys][plottype] = self.cloneDummyHistogram( self.finalhist["nominal"] )
                self.systematic_hist[sys][plottype].SetName(self.shortname+"_"+channel+suffix)
                #print "systematic hist name ",self.shortname+"_"+channel+suffix
                for iname, samplename in enumerate(allsamples_combine):
                    if self.sample_systematic_hist[samplename][sys][plottype]:
                        self.systematic_hist[sys][plottype].Add( self.sample_systematic_hist[samplename][sys][plottype] )
                    else:
                        ##if sys shift not available, still need to add the nominal 
                        self.systematic_hist[sys][plottype].Add( self.sample_systematic_hist[samplename]["nominal"] )
                #print "plottype ",plottype, self.systematic_hist[sys][plottype].Print("ALL")," nominal ",self.finalhist["nominal"].Print("ALL")
                self.systematic_hist[sys][plottype].SetDirectory(0)
                for bin in xrange(self.systematic_hist[sys][plottype].GetNbinsX()):
                    sysvalue = self.systematic_hist[sys][plottype].GetBinContent(bin+1) - self.finalhist["nominal"].GetBinContent(bin+1) ##error: hist_up - hist_nominal
                    #print "plottype ",plottype, " bin ",bin," systematic_hist ",self.systematic_hist[sys][plottype].GetBinContent(bin+1)," nominal ",self.finalhist["nominal"].GetBinContent(bin+1)  ," sysvalue ",sysvalue ," dummyplot bin content ", dummyplots[plottype].GetBinContent(bin+1)
                    totalSys[sys][plottype] +=  sysvalue
                    bincontent = dummyplots[plottype].GetBinContent(bin+1) + sysvalue * sysvalue
                    dummyplots[plottype].SetBinContent(bin+1, bincontent)
                totalSys[sys][plottype] = abs(totalSys[sys][plottype])/self.finalhist["nominal"].Integral()

            for bin in xrange( dummyplots[plottype].GetNbinsX()):
                bincontent = dummyplots[plottype].GetBinContent(bin+1)
                dummyplots[plottype].SetBinContent(bin+1, sqrt(bincontent))


            self.finalhist[plottype] = self.finalhist["nominal"].Clone()
            self.finalhist[plottype].SetName(self.shortname+"_"+channel+plottype)
            addweight = 1.0
            #print "plottype ", plottype, " dummpyplot ",dummyplots[plottype].Print("ALL")
            if plottype == "down":
                addweight = -1.0
            self.finalhist[plottype].Add(dummyplots[plottype], addweight)
            self.finalhist[plottype].SetDirectory(0)


        #print "nominal ", self.finalhist["nominal"].Print("ALL")
        ##add one side 
        self.finalhist["one_sided"] = self.cloneDummyHistogram( self.finalhist["nominal"] )
        self.finalhist["one_sided"].SetName(self.shortname+"_"+channel+"one_sided")
        for bin in xrange(self.finalhist["one_sided"].GetNbinsX()):
            err_up = self.finalhist["up"].GetBinContent(bin + 1)
            err_down = self.finalhist["down"].GetBinContent(bin + 1)
            nominal = self.finalhist["nominal"].GetBinContent(bin + 1)
            onesided = (abs(nominal - err_up) + abs(nominal - err_down))/2.0
            #print "err_up ",err_up, " err_down ",err_down, " nominal ",nominal," oneside ",onesided
            #sys_stat_total = sqrt(onesided*onesided + nominal)
            self.finalhist["nominal"].SetBinError(bin+1, onesided)
            self.finalhist["one_sided"].SetBinContent(bin+1, nominal + onesided)


        #print "Integral nominal ", self.finalhist["nominal"].Integral(), " up ",self.finalhist["up"].Integral()," down ",self.finalhist["down"].Integral(), " one sided ",self.finalhist["one_sided"].Integral()
        #print "totalSys ",totalSys," self.systematic_hist ",self.systematic_hist
        for sys in self.systematiclist:
            totalSys[sys]["oneside"] = (totalSys[sys]["up"] + totalSys[sys]["down"])/2.0
            #print "totalSys ",sys," samplename ",self.shortname, " up ",totalSys[sys]["up"], " down " ,totalSys[sys]["down"]," oneside ",totalSys[sys]["oneside"]," Up histname ",self.systematic_hist[sys]["up"].GetName()


    #################################################################################################
    ### 
    #################################################################################################
    def runSystematics(self, shortname, channel):
        self.plotSystematics(shortname, channel)
        if shortname == "DY" and self.DYdatadriven and (channel == "MuMu" or channel == "ElEl"):
            datadriven_finalhist ={}
            for name in self.shortnames_datadriven:
                self.combineSystematics(name, channel)
                datadriven_finalhist[name] = self.finalhist.copy()
                self.channel_shortname_systematic_hist[channel][name] = self.systematic_hist.copy()
                #print "Data-driven ", datadriven_finalhist[name]["nominal"].Print("ALL")

            ### calculate final hist for data-driven DY    
            self.finalhist = {}
            for plottype in ["nominal","up","down","one_sided"]:
                self.finalhist[plottype] = self.cloneDummyHistogram(self.hist_data_untagged)

            for  bin in xrange(self.hist_data_untagged.GetNbinsX()):
                err_up_square = 0.0
                err_down_square = 0.0
                nominal_mc = 0.0
                data_binvalue = self.hist_data_untagged.GetBinContent(bin+1)
                for name in self.shortnames_datadriven:
                    thisnominal = datadriven_finalhist[name]["nominal"].GetBinContent(bin+1)
                    thisup = datadriven_finalhist[name]["up"].GetBinContent(bin+1)
                    thisdown = datadriven_finalhist[name]["down"].GetBinContent(bin+1)
                    err_up_square = err_up_square + (thisnominal-thisup)*(thisnominal-thisup)
                    err_down_square = err_down_square + (thisnominal-thisdown)*(thisnominal-thisdown)
                    nominal_mc  = nominal_mc+thisnominal
                one_sided_err = (sqrt(err_up_square)+sqrt(err_down_square))/2.0
                self.finalhist["nominal"].SetBinContent(bin+1, data_binvalue - nominal_mc)
                self.finalhist["up"].SetBinContent(bin+1, data_binvalue - nominal_mc + sqrt(err_up_square))
                self.finalhist["down"].SetBinContent(bin+1, data_binvalue - nominal_mc - sqrt(err_down_square))
                self.finalhist["one_sided"].SetBinContent(bin+1, data_binvalue - nominal_mc + one_sided_err)
                self.finalhist["nominal"].SetBinError(bin+1, one_sided_err)
            #print "Data-Driven DY final ", self.finalhist["nominal"].Print("ALL")


        else:
            self.combineSystematics(shortname, channel)
            self.channel_shortname_systematic_hist[channel][shortname] = self.systematic_hist.copy()
        #print "shortname ",shortname, " all process systematics ", self.channel_shortname_systematic_hist

    def addTH1withError(self, hist1, hist2, c2=1.0):
        
        h_dummy = hist1.Clone()
        for bin in xrange(h_dummy.GetNbinsX()):
            value1 = hist1.GetBinContent(bin+1)
            err1 = hist1.GetBinError(bin+1)
            value2 = hist2.GetBinContent(bin+1)
            err2 = hist2.GetBinError(bin+1)
            h_dummy.SetBinContent(bin+1, value1+ value2*c2)
            h_dummy.SetBinError(bin+1, sqrt(err1*err1 + err2*err2*c2*c2))

        #h_dummy.Print("ALL")
        return h_dummy
            
    def cloneDummyHistogram(self,rHistogram):
        """Clone a histogram and set content/error to 0 -- for 1-sided systematics"""
        h_dummy = rHistogram.Clone()
        h_dummy.SetName("dummy")
        for bin in xrange(h_dummy.GetNbinsX()):
            h_dummy.SetBinContent(bin+1,0.0)
            h_dummy.SetBinError(bin+1,0.0)

        return h_dummy

    def writeSystematicsToFile(self, directory):

        ## how to write histogram into file 
        self.outfile = os.path.join(directory, self.shortname+"_systematics.root")
        print "self.outfile ", self.outfile
        tfile = ROOT.TFile(self.outfile, "recreate")
        for samplename in self.sample_systematic_hist.keys():
            #self.sample_systematic_hist[samplename]["nominal"].Print("ALL")
            self.sample_systematic_hist[samplename]["nominal"].SetDirectory(tfile)
            self.sample_systematic_hist[samplename]["nominal"].Write()
            for sys in self.sample_systematic_hist[samplename].keys():
                for plottype in ["up","down"]:
                    self.sample_systematic_hist[samplename][sys][plottype].SetDirectory(tfile)
                    self.sample_systematic_hist[samplename][sys][plottype].Write()

        self.systematic_hist["nominal"].SetDirectory(tfile)
        self.systematic_hist["nominal"].Write()
        for sys in self.systematic_hist.keys():
            for plottype in ["up","down"]:
                self.systematic_hist[sys][plottype].SetDirectory(tfile)
                self.systematic_hist[sys][plottype].Write()

        """
        for plotytype in ["nominal", "up","down"]
            self.finalhist[plottype].SetDirectory(tfile)
            self.finalhist[plottype].Write()
        """

        tfile.Close()




### THE END ###


   


