#!/usr/bin/python
import os
import sys


def generate_condorjob(mass, hmev, nnv):
    condorscript = open("condor_rebin_mass%d_nnv%dHMEbinv%d.cmd"%(mass, nnv, hmev),"w")

    condorscript.write("""universe              = vanilla 
executable            = runrebin.sh
arguments             = {mass} {hmeversion} {nnversion}
output                = output/rebin_M{mass}_hmev{hmeversion}_nnv{nnversion}.$(ClusterId).$(ProcId).out
error                 = error/rebin_M{mass}_hmev{hmeversion}_nnv{nnversion}.$(ClusterId).$(ProcId).err
log                   = log/rebin_M{mass}_hmev{hmeversion}_nnv{nnversion}.$(ClusterId).log
request_cpus          = 2
request_memory        = 4000M
+JobFlavour           = "workday"
Notification          = Complete
queue""".format(mass = mass, hmeversion = hmev, nnversion=nnv))
    condorscript.close()

#notify_user           = taohuang@email.tamu.edu

masspoints = [260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 900]
#masspoints = [700]
HME_bins = [0, 1, 2, 3]
#nn_bins = [5, 8, 10, 20]
nn_bins = [10,20]
submitscript = open("submit_rebincondorjobs.sh","w")
submitscript.write("rm output/* \n")
submitscript.write("rm error/* \n")
submitscript.write("rm log/* \n")
for mass in masspoints:
    for hmev in HME_bins:
      for nnv in nn_bins:
          #generate_condorjob(mass, hmev, nnv)
          filename = "condor_rebin_mass%d_nnv%dHMEbinv%d.cmd"%(mass, nnv, hmev)
          submitscript.write("condor_submit "+filename+" \n")

submitscript.close()
os.system("chmod +x submit_rebincondorjobs.sh")
