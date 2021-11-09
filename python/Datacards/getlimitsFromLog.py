
import re

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
    print("limits_lines ", limits_lines)
    for key in percents:
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
    print limits
    return limits


def extractlimitfromtxtfile_t100(logfile):
    logopen = open(logfile, "read")
    signal_xsec = 1.0
    #signal_xsec = 5000.0
    percents = [2.5, 16.0, 50.0, 84.0, 97.5, -1]
    limits_lines = {}
    limits = {}
    ntoys = 0; observed_tot = 0.0
    observed_list = []
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
    #limits[-1] = observed_tot/ntoys
    print "observed_list ",observed_list
    from numpy import median
    limits[-1] = median(observed_list)
    print limits
    return limits



logfile = "GGToX0ToHHTo2B2L2Nu_nnout_MTandMT2_MJJ_nnstep0p04_nncut0p12_limits/M270/GGToX0ToHHTo2B2L2Nu_M270_MuMu_ElEl_MuEl_signalxsec1fb.log"
#extractlimitfromtxtfile(logfile)
logfile = "GGToX0ToHHTo2B2L2Nu_M900_MuMu_ElEl_MuEl_signalxsec1fb.log"
extractlimitfromtxtfile_t100(logfile)
