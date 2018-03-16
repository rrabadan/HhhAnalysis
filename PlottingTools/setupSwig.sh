#!/bin/bash


swig -c++ -python mt2_bisect.i 
g++ -O2 -fPIC -c mt2_bisect.cpp
g++ -c -m64 -fPIC mt2_bisect.cpp mt2_bisect_wrap.cxx -I$CMSSW_BASE/python
g++ -shared mt2_bisect.o mt2_bisect_wrap.o -o _mt2_bisect.so 
