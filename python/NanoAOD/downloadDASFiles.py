import os,sys
import argparse

parser = argparse.ArgumentParser(description='downloadDASfiles')
parser.add_argument("-d", "--dataset", dest="dataset", default="none", help="dataset to download from DAS")
parser.add_argument("-o", "--outputdir",dest="output", default="/fdata/hepx/store/user/taohuang/HH_NanoAOD/",help="output root file directory, [Defualt:fdata] ")
args = parser.parse_args()
#dataset = "/GluGluToRadionToHHTo2B2VTo2L2Nu_M-*_narrow_13TeV-madgraph-v2/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM"
dataset = args.dataset
query = "file dataset="+dataset
print "To download ",dataset, " to ", args.output
flist = os.popen("dasgoclient -limit=0 -query='{query}'".format(query = query))
print "find datafile "
#os.system("dasgoclient -limit=0 -query='{query}' > {filelist}".format(query = query, filelist = filelist))
#output_folder = "/fdata/hepx/store/user/taohuang/HH_NanoAOD/"
output_folder = args.output
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
#flist = open(filelist, "read")
for line in flist:
    print "line ",line[:-1]
    if ".root" in line and "NANOAOD" in line:
    	print "download... "
	os.system("xrdcp root://cms-xrd-global.cern.ch/"+line[:-1]+" "+output_folder)
