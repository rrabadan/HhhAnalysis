#include <iostream>
#include <utility>
#include <vector>
#include <cmath>
#include "TTree.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TRefArray.h"
#include "TObject.h"
#include "TString.h"
#include "HhhAnalysis/CutFlowAnalyzer/interface/CrossSections_and_BR.h"

using namespace std;
// This is in pb
CrossSections_and_BR::CrossSections_and_BR(){
  //BR
  float BR_h_bb   = 0.577;
  float BR_h_WW   = 0.215;
  float BR_W_lnu  = 0.3257;
  //float BR_W_munu = 0.1057;
  //float BR_W_munu_taunu = 0.2182;
  //float BR_t_Wb   = 1.;
  // 0.577^2 + 0.215^2*0.3257^4 + 2*0.577*0.215*0.3257^2 -> Correct for our samples: bbbb + bbWW + WWWW
  // 2*0.577*0.215*0.3257^2 -> Just to know muv Br: bbWW -> bblvlv
  // 2*0.577*0.215*0.1057^2 -> Just to know muv Br: bbWW -> bbmuvmuv
  float BR_final = pow(BR_h_bb,2) + pow(BR_h_WW,2)*pow(BR_W_lnu,4) + 2*BR_h_bb*BR_h_WW*pow(BR_W_lnu,2);
  //Background (BR already included)
  XsexBr_tt_To_2l2n   = 87.31; //L=e, mu only?
  XsecBr_DYJets_ToLL  = 18610; //L=e, mu only?
  XsecBr_DY0Jets_ToLL = 4758.9;
  XsecBr_DY1Jets_ToLL = 929.1;
  XsecBr_DY2Jets_ToLL = 337.1;
  //Optimistic Xsecs
  XsecBr_B1  = 1.1948649279633416   * 0.5035043253778491  * BR_final;
  XsecBr_B2  = 0.5900264559172156   * 0.7401476862292035  * BR_final;
  XsecBr_B3  = 0.4431392872406074   * 0.7550948295752251  * BR_final;
  XsecBr_B4  = 0.3629024495635394   * 0.3271838865481416  * BR_final;
  XsecBr_B5  = 0.26198989157860825  * 0.30714428242103997 * BR_final;
  XsecBr_B6  = 0.1505457803080438   * 0.23554145150771982 * BR_final;
  XsecBr_B7  = 0.08778859336272271  * 0.23284593958553115 * BR_final;
  XsecBr_B8  = 0.044958592682492124 * 0.30234102603499013 * BR_final;
  XsecBr_B9  = 0.02330230741512451  * 0.3283405583212012  * BR_final;
  XsecBr_B10 = 0.016878791231549145 * 0.19895423048929178 * BR_final;
  XsecBr_B11 = 0.008205926642337763 * 0.22349301758113208 * BR_final;
  XsecBr_B12 = 0.006805317561386914 * 0.07869395182066141 * BR_final;
}

CrossSections_and_BR::~CrossSections_and_BR(){
}

// By definition is nitialize to 1 fb-1 (so you can multiply later to the correct lumi) 
float CrossSections_and_BR::GetWeight( int proc ){
  float xsecbr = -999.;
  if(proc==B1)           xsecbr = XsecBr_B1;  
  else if(proc==B2)      xsecbr = XsecBr_B2;  
  else if(proc==B3)      xsecbr = XsecBr_B3;  
  else if(proc==B4)      xsecbr = XsecBr_B4;  
  else if(proc==B5)      xsecbr = XsecBr_B5;  
  else if(proc==B6)      xsecbr = XsecBr_B6;  
  else if(proc==B7)      xsecbr = XsecBr_B7;  
  else if(proc==B8)      xsecbr = XsecBr_B8;  
  else if(proc==B9)      xsecbr = XsecBr_B9;  
  else if(proc==B10)     xsecbr = XsecBr_B10; 
  else if(proc==B11)     xsecbr = XsecBr_B11; 
  else if(proc==B12)     xsecbr = XsecBr_B12; 
  else if(proc==tt)      xsecbr = XsexBr_tt_To_2l2n;
  else if(proc==DYJets)  xsecbr = XsecBr_DYJets_ToLL;
  else if(proc==DY0Jets) xsecbr = XsecBr_DY0Jets_ToLL;
  else if(proc==DY1Jets) xsecbr = XsecBr_DY1Jets_ToLL;
  else if(proc==DY2Jets) xsecbr = XsecBr_DY2Jets_ToLL;
  else cout<<"WARNING::CrossSections_and_BR -> No process specified!"<<endl;
  //Lumi*=1000; //convert to pb-1
  //cout<<"CROSS SECTION: "<<xsecbr<<" * "<<Lumi<<" / "<<InitEv<<" = "<<xsecbr*Lumi/InitEv<<endl;
  return xsecbr;//*Lumi/InitEv;
}
