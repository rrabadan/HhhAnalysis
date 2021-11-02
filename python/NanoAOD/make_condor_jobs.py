import os
pwd = os.getcwd()+'/'

for year in [2016]:
  for mass in [260, 270, 300, 350, 400, 450, 500, 550, 600, 650, 750, 800, 900, 1000]:
    os.system("mkdir -p condor")
    os.system("mkdir -p condor/{year}".format(year = year))
    os.system("mkdir -p condor/{year}/m{mass}".format(mass = mass, year = year))
    os.system("mkdir -p condor/{year}/m{mass}/output".format(mass = mass, year = year))
    os.system("mkdir -p condor/{year}/m{mass}/error".format(mass = mass, year = year))
    os.system("mkdir -p condor/{year}/m{mass}/log".format(mass = mass, year = year))
    os.system("mkdir -p condor/{year}/m{mass}/scripts".format(mass = mass, year = year))


    c_name = "condor/{year}/m{mass}/scripts/condor_m{mass}.sh".format(mass = mass, year = year)
    c_script = open(c_name, "write")
    c_script.write("#!/bin/bash\n")
    c_script.write("date\n")
    c_script.write("cd "+pwd+"\n")
    c_script.write("eval `scramv1 runtime -sh`\n")
    c_script.write("python postproc_radion750_2016.py {mass} condor/{year}/m{mass}/\n".format(mass = mass, year = year))
    c_script.write("echo 'done'\n")
    c_script.write("date\n")

    os.system("chmod 755 "+c_name)

    condor_sub = open("condor/{year}/m{mass}/scripts/submit_m{mass}.sh".format(mass = mass, year = year), "write")
    condor_sub.write("""universe                = vanilla
executable              = {pwd}/{fname}
arguments               = no
output                  = {pwd}/condor/{year}/m{mass}/output/out_m{mass}.$(ClusterId).$(ProcId).out
error                   = {pwd}/condor/{year}/m{mass}/error/err_m{mass}.$(ClusterId).$(ProcId).err
log                     = {pwd}/condor/{year}/m{mass}/log/log_m{mass}.$(ClusterId).log
request_memory          = 4000M
+JobFlavour             = "workday"
queue""".format(fname = c_name, pwd = pwd, mass = mass, year = year))
