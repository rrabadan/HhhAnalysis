#!/usr/bin/python
import os
from optparse import OptionParser

def getid(line):
    parts=line.split("\"")
    return eval(parts[1])

def getiso(line):
    parts=line.split("\"")
    return eval(parts[7])

def main():

    wplines = [251,281,293,311,323]
    wp = [70,80,85,90,95]
    # You can give a List if you have more
    files = ['weights/OptimizedCuts.xml']

    etabins = [0,1,2,0,1,2]
    ptbins =  [0,0,0,1,1,1]
    etalabels = ['/&eta;/<0.8','0.8</&eta;/<1.479','1.479</&eta;/<2.5']

    # 2 pt bins x 3 eta bins x 5 WPs
    idtable = [ [ [ 999 for i in range(5) ] for j in range(3) ] for k in range(2) ]
    isotable = [ [ [ 999 for i in range(5) ] for j in range(3) ] for k in range(2) ]

    bin=0
    for ifile in files:
        etabin = etabins[bin]
        ptbin = ptbins[bin]

        eleid = []
        eleiso = []
        f = open(ifile)
        counter=0
        for line in f:
            counter = counter+1
            for i in range(len(wp)):
                if(counter==wplines[i]):
                    print str(getid(line))
                    print str(getiso(line))
                    print "---"
                    eleid.append(getid(line))
                    eleiso.append(getiso(line))

        for i in range(len(eleid)):
            idtable[ptbin][etabin][i] = eleid[i]
            isotable[ptbin][etabin][i] = eleiso[i]
        bin = bin+1


    # write the table now
    for ipt in range (2):
        print "pt bin = "+str(ipt)
        print "| *WP* | ",
        for ieta in range(3):
            print etalabels[ieta]+" | ",
            if(ieta==2):
                print ""
        for iwp in range(5):
            print "| *WP "+str(wp[iwp])+"* | ",
            for ieta in range(3):
                print "%.3f" % idtable[ipt][ieta][iwp]+" / %.3f" % isotable[ipt][ieta][iwp]+" | ",
                if(ieta==2):
                    print ""


if __name__ == "__main__":
    main()
