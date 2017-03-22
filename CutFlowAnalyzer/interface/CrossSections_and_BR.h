#ifndef CrossSections_and_BR_h
#define CrossSections_and_BR_h
#include <iostream>
#include <utility>
#include <vector>
#include "TString.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TChain.h"

class CrossSections_and_BR{
  public:
    CrossSections_and_BR();
    ~CrossSections_and_BR();
    float GetWeight( int proc ); //Return Xsec * BR
  private:
    enum {B1 = 1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, tt, DYJets, DY0Jets, DY1Jets, DY2Jets};
    Float_t XsecBr_B1;
    Float_t XsecBr_B2;
    Float_t XsecBr_B3;
    Float_t XsecBr_B4;
    Float_t XsecBr_B5;
    Float_t XsecBr_B6;
    Float_t XsecBr_B7;
    Float_t XsecBr_B8;
    Float_t XsecBr_B9;
    Float_t XsecBr_B10;
    Float_t XsecBr_B11;
    Float_t XsecBr_B12;
    Float_t XsexBr_tt_To_2l2n;
    Float_t XsecBr_DYJets_ToLL;
    Float_t XsecBr_DY0Jets_ToLL;
    Float_t XsecBr_DY1Jets_ToLL;
    Float_t XsecBr_DY2Jets_ToLL;
};
#endif
