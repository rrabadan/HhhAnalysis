#!/usr/bin/python
import os
import sys
#Runyear, Slist from Samplelist
from RunConfiguration  import *
import ROOT

MCNanoAODGTlib = {
    2016 : "RunIISummer16NanoAODv3",
    2017 : "RunIIFall17NanoAODv4",
    2018 : "RunIIAutumn18NanoAODv4"
}
#DateGlobalTag = "05Feb2018"
#MCGlobalTag = "RunIISummer16NanoAOD"
newRunyear = 2016
newDateGT = "Nano14Dec2018"
newMCNanoAODGT = MCNanoAODGTlib[newRunyear]
oldDateGT = Slist.DateGlobalTag
oldMCNanoAODGT = Slist.MCNanoAODGlobalTag 



torun_datasets = Slist.Nanodatasets
#for mass in Slist.masspoints:
#    torun_datasets.append("/GluGluToRadionToHHTo2B2VTo2L2Nu_M-%d_narrow_13TeV-madgraph_correctedcfg/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM"%mass)
#print "=============================================================="
#print "outputdir ",outdir
#print "=============================================================="


newnanoaod = open("NewNanoAODSample%d.txt"%newRunyear,"w")
newSamplelist = open("newSamplelist%d.py"%newRunyear, "w")
newSamplelist.write("""#!/usr/bin/python                                                                                                                 import os
import sys

datasets  = []; NumSample = []; sampleN_short = []
Nanodatasets = []; localdirs = {}
MCxsections = []        

""")

newSamplelist.write("DateGlobalTag= '%s'\n"%newDateGT)
newSamplelist.write("MCNanoAODGlobalTag= '%s'\n\n"%newMCNanoAODGT)

def updateSamplelist(torun_datasets):
    for ijob, job in enumerate(torun_datasets):
	index = Slist.Nanodatasets.index(job)
	nsample = int(Slist.NumSample[index])
	jobtype = Slist.sampleN_short[index]
	Nanodataset =  Slist.Nanodatasets[index]
	mc_cs = Slist.MCxsections[index]
	if index ==0 or jobtype != Slist.sampleN_short[index-1]:
	    newSamplelist.write("\n### "+jobtype+"\n")
	sampletype = ""
	sampletag = ""
	#print "nsample ",nsample, " jobtype ",jobtype, "dataset ", job," NanoAOD ",Nanodataset
	if job.split('/')[1] != Nanodataset.split('/')[1]:
	    print "dataset matching is wrong!! job is ",job," NanoAOD is ",Nanodataset
	if nsample < 0:
	    datadir = Slist.sampleN_short[index]
	    dataname = job
	    sampletype = "NANOAOD"
	    #print "real data nsample ",nsample, " datadir ",datadir
	elif nsample>= 0:
	    datadir = job.split('/')[1]
	    sampletype = "NANOAODSIM"
	    sampletag = "%s*%s*"%(newMCNanoAODGT, newDateGT)
	    #print "MC nsample ",nsample, " datadir ",datadir," sampletag ",sampletag


	
	query = "dataset=/%s/%s/%s"%(datadir, sampletag, sampletype)
	if nsample < 0 and oldDateGT in Nanodataset:
	    query = Nanodataset.replace(oldDateGT, newDateGT)
    	    if newRunyear != Runyear and 'Run%d'%Runyear in query:
	          querylist = query.split('/')
		  query = '/'+querylist[1]+'/'+'Run%d'%newRunyear+querylist[2][7]+"*"+newDateGT+"*"+'/NANOAOD'
        elif  nsample < 0:
            print "warning!!! no good query for data ", query
	

	#print "query ",query
	#os.system("dasgoclient -limit=0 -query='{query}' >> {outfilename}".format(query = query, outfilename = outfile))
	#flist = os.peopn("dasgoclient -limit=0 -query='{query}' >> {outfilename}".format(query = query, outfilename = outfile))
	flist = os.popen("dasgoclient -limit=0 -query='{query}'".format(query = query), 'r')
	ifile = 0
	for line in flist:
	    ifile += 1
	    #print "line ",line[:-1]
	    if "NANOAOD" in line :#and "NANOAOD" in line:
	    	newdas = line[:-1]
		newnanoaod.write(line)
    		newSamplelist.write("Nanodatasets.append('%s')\n"%newdas)
		newSamplelist.write("NumSample.append('%d'); sampleN_short.append('%s'); MCxsections.append(%f)\n"%(nsample, jobtype, mc_cs))
	if ifile > 1 or ifile == 0:
	   print "Error!! found #samples :%d"%ifile," check sample by hand ",Nanodataset," query ",query
	   newSamplelist.write("##Error!! found #samples :%d,check the sample by hand with query %s"%(ifile, query))
	   if ifile == 0:
	       newSamplelist.write("##NumSample.append('%d'); sampleN_short.append('%s'); MCxsections.append(%f)\n"%(nsample, jobtype, mc_cs))
                



updateSamplelist(torun_datasets)
