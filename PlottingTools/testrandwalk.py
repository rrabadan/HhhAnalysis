
import numpy as np
import ROOT
refPDF = ROOT.TFile("REFPDFPU40_test.root","READ")
onshellWmasspdf = refPDF.Get("onshellWmasspdf")
def getOnshellWMass(onshellWmasspdf, x0, step, random):
    
    xmin = 50.0; xmax = 90.0
    while (x0 > xmax or x0 < xmin):
	if x0 > xmax:
	    x0 = x0 - xmax + xmin
	if x0 < xmin:
	    x0 = xmax - (xmin - x0)
    x1 = x0 + step
    while (x1 > xmax or x1 < xmin):
	if x1 > xmax:
	    x1 = x1 - xmax + xmin
	if x1 < xmin:
	    x1 = xmax - (xmin - x1)
    """
    binx00 = onshellWmasspdf.FindBin(x0)
    binx01 = binx00 + 1
    binx10 = onshellWmasspdf.FindBin(x1)
    binx11 = binx10 + 1
    if onshellWmasspdf.GetBinCenter(binx00) > x0:
	binx00 = binx00 - 1
	binx01 = binx01 - 1
    if onshellWmasspdf.GetBinCenter(binx10) > x1:
	binx10 = binx10 - 1
	binx11 = binx11 - 1
    y00 = onshellWmasspdf.GetBinContent( binx00 )
    x00 = onshellWmasspdf.GetBinCenter( binx00 )
    y01 = onshellWmasspdf.GetBinContent( binx01 )
    x01 = onshellWmasspdf.GetBinCenter( binx01 )
    y10 = onshellWmasspdf.GetBinContent( binx10 )
    x10 = onshellWmasspdf.GetBinCenter( binx10 )
    y11 = onshellWmasspdf.GetBinContent( binx11 )
    x11 = onshellWmasspdf.GetBinCenter( binx11 )
    w0  = (x0 - x00)*(y00 - y01)/(x00 - x01) + y00;
    w1  = (x1 - x10)*(y10 - y11)/(x10 - x11) + y10;
    """
    w0  = onshellWmasspdf.Interpolate(x0)
    w1  = onshellWmasspdf.Interpolate(x1)
    #w1/w0: transition probability 
    if (w1/w0 >= random):
	return x1
    elif (w1/w0 < random):
	return x0
    else:	
	print "error in getOnshellWMass "
	return 80.3

iterations = 100000
wfile = ROOT.TFile("randwalk.root","recreate")
tree = ROOT.TTree("wmasstree","Wmass tree")
wmass_gen1 = np.zeros(1, dtype=float)
wmass_gen2 = np.zeros(1, dtype=float)
wmass_gen3 = np.zeros(1, dtype=float)
wmass_gen4 = np.zeros(1, dtype=float)
wmass_gen5 = np.zeros(1, dtype=float)
tree.Branch('wmass_gen1', wmass_gen1, 'wmass_gen1/D')
tree.Branch('wmass_gen2', wmass_gen2, 'wmass_gen2/D')
tree.Branch('wmass_gen3', wmass_gen3, 'wmass_gen3/D')
tree.Branch('wmass_gen4', wmass_gen4, 'wmass_gen4/D')
tree.Branch('wmass_gen5', wmass_gen5, 'wmass_gen4/D')
wmass_gen1[0] = 70.0
wmass_gen2[0] = 70.0
wmass_gen3[0] = 70.3
wmass_gen4[0] = 70.3
wmass_gen5[0] = 70.3
genRandom = ROOT.TRandom3(0)
it = 0
while it<iterations:
    it += 1 
    rand01 = genRandom.Uniform(0., 1.0)
    step1 = genRandom.Uniform(-4.0, 4.0)
    step2 = genRandom.Uniform(-8.0, 8.0)
    step3 = genRandom.Uniform(-16.0, 16.0)
    step4 = genRandom.Uniform(-32.0, 32.0)
    step5 = genRandom.Uniform(-64.0, 64.0)
    #print "step1 ",step1," step2 ",step2," step3 ",step3," step4 ",step4
    wmass_gen1[0] = getOnshellWMass(onshellWmasspdf, wmass_gen1[0], step1, rand01)
    wmass_gen2[0] = getOnshellWMass(onshellWmasspdf, wmass_gen2[0], step2, rand01)
    wmass_gen3[0] = getOnshellWMass(onshellWmasspdf, wmass_gen3[0], step3, rand01)
    wmass_gen4[0] = getOnshellWMass(onshellWmasspdf, wmass_gen4[0], step4, rand01)
    wmass_gen5[0] = getOnshellWMass(onshellWmasspdf, wmass_gen5[0], step5, rand01)
    tree.Fill() 
    
tree.Write()
onshellWmasspdf.Write()
wfile.Close()
