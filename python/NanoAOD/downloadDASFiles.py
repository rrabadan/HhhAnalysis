import os,sys

dataset = "/GluGluToRadionToHHTo2B2VTo2L2Nu_M-*_narrow_13TeV-madgraph-v2/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM"
query = "file dataset="+dataset
filelist = "filelist.txt"
flist = os.popen("dasgoclient -limit=0 -query='{query}'".format(query = query))
#os.system("dasgoclient -limit=0 -query='{query}' > {filelist}".format(query = query, filelist = filelist))
output_folder = "/fdata/hepx/store/user/taohuang/HH_NanoAOD/"
#flist = open(filelist, "read")
masslist = [260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 750, 800, 900]
for line in flist:
    print "line ",line[:-1]
    for mass in masslist:
        massstr = "GluGluToRadionToHHTo2B2VTo2L2Nu_M-%d_narrow_13TeV-madgraph-v2"%mass
	thisfolder = os.path.join(output_folder, massstr)
	if not os.path.exists(thisfolder):
	    os.makedirs(thisfolder)
	if massstr in line:
	    print "mass ",mass," file  ",line[:-1]
	    os.system("xrdcp root://cms-xrd-global.cern.ch/"+line[:-1]+" "+thisfolder)
	    break
#for i in range(1, 20):
    #os.system("xrdcp root://cms-xrd-global.cern.ch//store/user/arizzi/Nano01Fall17/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAOD-94X-Nano01Fall17/180205_181145/0000/test94X_NANO_{nfile}.root /fdata/hepx/store/user/taohuang/HH_NanoAOD/".format(nfile = i))
