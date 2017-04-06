import FWCore.ParameterSet.Config as cms
from HhhAnalysis.MCProduction.fileNamesPU import fileNamesPU

def addPileup(process):
    process.mix.input.fileNames = fileNamesPU
    return process

