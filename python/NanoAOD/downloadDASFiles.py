import os,sys
import argparse

parser = argparse.ArgumentParser(description='downloadDASfiles')
parser.add_argument("-d", "--dataset", dest="dataset", default="none", help="dataset to download from DAS")
parser.add_argument("-o", "--outputdir",dest="output", default="/fdata/hepx/store/user/taohuang/HH_NanoAOD/",help="output root file directory, [Defualt:fdata] ")
args = parser.parse_args()
#2016
#dataset = "/GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-*_narrow_13TeV-madgraph-v2/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM"
#2017
#dataset = "/GluGluToRadionToHHTo2B2VTo2L2Nu_M-280_narrow_13TeV-madgraph_correctedcfg/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM"
dataset = args.dataset
query = "file dataset="+dataset
print "To download ",dataset, " to ", args.output
flist = os.popen("dasgoclient -limit=0 -query='{query}'".format(query = query))
print "find datafile "
#os.system("dasgoclient -limit=0 -query='{query}' > {filelist}".format(query = query, filelist = filelist))
#output_dir = "/fdata/hepx/store/user/taohuang/HH_NanoAOD/"

output_folder = args.output
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
#flist = open(filelist, "read")
masslist = [260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900]
for line in flist:
    print "line ",line[:-1]
    if ".root" in line and "NANOAOD" in line:
#        for mass in masslist:
#	    if "M-%d"%mass in line:
#		#output_folder = args.output
#	        output_folder = os.path.join(output_dir, "GluGluToBulkGravitonToHHTo2B2VTo2L2Nu_M-%d_narrow_13TeV-madgraph_correctedcfg"%mass)
#		if not os.path.exists(output_folder):
#		    os.makedirs(output_folder)
		print "download... ",output_folder
		os.system("xrdcp root://cms-xrd-global.cern.ch/"+line[:-1]+" "+output_folder)
