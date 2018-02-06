import ROOT
import os

masses = [260, 270, 300, 350, 400, 450, 500, 550, 600, 650,750, 800, 900]
filecounts = {"TT":0, "sT_top":0, "sT_antitop":0, "DYM10to50":0, "DYToLL0J":0, "DYToLL1J":0,"DYToLL2J":0}
for mass in masses:
   filecounts["radion_M%d"%mass] = 0

def checkSampletype(filename ):
    for name in filecounts:
        if filename.startswith(name):
	    filecounts[name] += 1
	    return name
def checkOneSampletype(filedir, sampletype, num):
    for i in range(0, num):
        f = sampletype + "_ijob%d.root"%i
        if os.path.isfile(filedir+f):
	   ch = ROOT.TChain("t")
	   ch.Add(filedir+f)
	   if not ch.GetEntries():
	       print "f ",f
	else:
	     print "f not exist ", f
##sbatch 20180205_10k_Signal_TT_Louvain/Send_PlotterProducer_radion_M260_1_10000.slrm
def checkAllFiles(filedir, outfile ):
    allfiles = os.listdir(filedir)
    jobscriptfile = open(outfile, "write")
    for f in allfiles:
    #for f in ['TT_ijob294.root']:
	if not f.endswith(".root"):
	    continue
	samplename = checkSampletype(f)
	ch = ROOT.TChain("t")
	ch.Add(filedir+f)
	#print "file  ",f," entries ",ch.GetEntries()
	if not ch.GetEntries():
	   print "f ",f
	   jobname = f[:-5].split('_')[-1]
	   jobnum = jobname[4:]
	   jobscriptfile.write("sbatch 20180205_10k_Signal_TT_Louvain/Send_PlotterProducer_%s_%s_10000.slrm\n"%(samplename, jobnum))
	
    print "files counting ",filecounts
filedir = "/fdata/hepx/store/user/taohuang/20180205_10k_Signal_TT_Louvain/"
outfile = "resubmitall20180205_10k_Louvain.sh"
#checkAllFiles(filedir, outfile)

checkOneSampletype(filedir, "TT",2000)
