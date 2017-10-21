
import os
#write a interface file:mt2_bisect.i
mtfile = "mt2_bisect"
pythonpath = "/home/taohuang/DiHiggsAnalysis/CMSSW_8_0_26_patch2/python"#use python from CMSSW
print "Compiling mt2_bisect with python, using Swig"
os.system("rm  %s.o & rm _%s.so & rm %s_wrap.*"%(mtfile, mtfile, mtfile))
print "generate _wrap.cxx" 
os.system("swig -c++ -python %s.i"%(mtfile))
os.system("g++ -O2 -fPIC -c %s.cpp"%(mtfile))
os.system("g++ -c -m64 -fPIC %s_wrap.cxx -I%s"%(mtfile, pythonpath))
os.system("g++ -shared %s.o %s_wrap.o -o _%s.so "%(mtfile, mtfile, mtfile))
print "Done"
#swig -c++ -python mt2_bisect.i 
#g++ -O2 -fPIC -c mt2_bisect.cpp
#g++ -c -m64 -fPIC mt2_bisect.cpp mt2_bisect_wrap.cxx -I/home/taohuang/DiHiggsAnalysis/CMSSW_8_0_26_patch2/python
#g++ -shared mt2_bisect.o mt2_bisect_wrap.o -o _mt2_bisect.so 
