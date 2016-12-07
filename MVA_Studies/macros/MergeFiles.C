#include <iostream>
#include <vector>
#include "TPad.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TSystem.h"
#include "TImage.h"
#include "TKey.h"
#include "TH1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TObjArray.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TString.h"

using namespace std;

void CopyDir( TFile * thisfile, TFile * outFile, TString CommonDir );

//.x MergeFiles.C+
void MergeFiles(){
  bool isB3 = true;
  TString BaseName = "RECO_1btag_40PU";
  TString Sig = isB3 ? "B3" : "B6";
  bool is_RuleFit=false, is_MLP=true;
  TString OutName = "TMVA_" + Sig  + "_" + BaseName + "_ALL.root";
  std::vector<TString> FilesToMerge;
  FilesToMerge.clear();
  FilesToMerge.push_back("TMVA_" + Sig + "_" + BaseName + "_Likelihood.root");
  FilesToMerge.push_back("TMVA_" + Sig + "_" + BaseName + "_LikelihoodD.root");
  FilesToMerge.push_back("TMVA_" + Sig + "_" + BaseName + "_LikelihoodMIX.root");
  FilesToMerge.push_back("TMVA_" + Sig + "_" + BaseName + "_BDT.root");
  FilesToMerge.push_back("TMVA_" + Sig + "_" + BaseName + "_HMatrix.root");
  if(is_RuleFit) FilesToMerge.push_back("TMVA_" + Sig + "_" + BaseName + "_RuleFit.root");
  if(is_MLP)     FilesToMerge.push_back("TMVA_" + Sig + "_" + BaseName + "_MLP.root");
  TFile *outFile=TFile::Open(OutName.Data(), "recreate");
  TTree TestTree("TestTree","TestTree");
  TTree TrainTree("TrainTree","TrainTree");
  TTree t1("t1","a simple Tree with simple variables");
  std::vector< std::pair<TString,TString> > NamesDir; NamesDir.clear();
  std::pair<TString,TString> Pair;
  Pair.first = "Method_Likelihood"; Pair.second = "Likelihood";
  NamesDir.push_back(Pair); 
  Pair.first = "Method_Likelihood"; Pair.second = "LikelihoodD";
  NamesDir.push_back(Pair);
  Pair.first = "Method_Likelihood"; Pair.second = "LikelihoodMIX";
  NamesDir.push_back(Pair);
  Pair.first = "Method_BDT"; Pair.second = "BDT";
  NamesDir.push_back(Pair);
  Pair.first = "Method_HMatrix"; Pair.second = "HMatrix";
  NamesDir.push_back(Pair);
  Pair.first = "Method_RuleFit"; Pair.second = "RuleFit";
  if(is_RuleFit) NamesDir.push_back(Pair);
  Pair.first = "Method_MLP"; Pair.second = "MLP";
  if(is_MLP) NamesDir.push_back(Pair);

  //Copy common files
  TFile *f = TFile::Open( FilesToMerge[0].Data() );
  if( !f ){ cout<<"WARNING1! "<<FilesToMerge[0].Data()<<" cannot be opened!"<<endl; return; }
  TH2F *h1 = (TH2F*) f->Get("CorrelationMatrixS");
  TH2F *h2 = (TH2F*) f->Get("CorrelationMatrixB");
  outFile->cd();
  h1->Write();
  h2->Write();
  f->Close();
  vector< vector<float> > MyTest;
  vector< vector<float> > MyTrain;
  //Copy specific folder and take Braches
  for(unsigned int i=0; i<FilesToMerge.size(); i++){
    cout<<"OPENING: "<<FilesToMerge[i].Data()<<endl;
    TFile *thisfile = TFile::Open( FilesToMerge[i].Data() );
    if( !thisfile ){ cout<<"WARNING2! "<<FilesToMerge[i].Data()<<" cannot be opened!"<<endl; continue; }
    //Specific stuff
    CopyDir( thisfile, outFile, NamesDir[i].first );
    outFile->cd();
    float Bname = -999.;
    TTree *TteA = (TTree*) thisfile->Get("TestTree");
    TteA->SetBranchAddress(NamesDir[i].second, &Bname);
    Int_t nentriesA = (Int_t)TteA->GetEntries();
    vector<float> MyTest_tmp; MyTest_tmp.clear();
    for (Int_t iT=0; iT<nentriesA; iT++) {
	TteA->GetEntry(iT);
	MyTest_tmp.push_back(Bname);
    }
    MyTest.push_back(MyTest_tmp);
    float BnameT = -999.;
    TTree *TttA = (TTree*) thisfile->Get("TrainTree");
    TttA->SetBranchAddress(NamesDir[i].second, &BnameT);
    nentriesA = (Int_t)TttA->GetEntries();
    vector<float> MyTrain_tmp; MyTrain_tmp.clear();
    for (Int_t iT=0; iT<nentriesA; iT++) {
	TttA->GetEntry(iT);
	MyTrain_tmp.push_back(Bname);
    }
    MyTrain.push_back(MyTrain_tmp);
    if( i==FilesToMerge.size()-1 ){
	std::vector<TString> CommonDir;
	CommonDir.clear();
	CommonDir.push_back("InputVariables_Id");
	CommonDir.push_back("InputVariables_Deco");
	CommonDir.push_back("InputVariables_PCA");
	CommonDir.push_back("InputVariables_Gauss_Deco");
	//Loop on Each common directory
	for(unsigned int j=0; j<CommonDir.size(); j++){
	  CopyDir( thisfile, outFile, CommonDir[j] );
	}//Common dir
cout<<"AA0"<<endl;
//	outFile->cd();
//	TTree *Tte = (TTree*) thisfile->Get("TestTree");
//	TTree *Ttt = (TTree*) thisfile->Get("TrainTree");
//	int classID, classID2;
//	char className, className2;
//	float dR_l1l2, dR_b1b2, dR_bl, mass_l1l2, mass_b1b2, MT2, weight;
//	float dR_l1l22, dR_b1b22, dR_bl2, mass_l1l22, mass_b1b22, MT22, weight2;
//	float Likelihood2, LikelihoodD2, LikelihoodMIX2, BDT2, HMatrix2, RuleFit2, MLP2;
//	Tte->SetBranchAddress("classID", &classID);
//	Tte->SetBranchAddress("className", &className);
//	Tte->SetBranchAddress("dR_l1l2", &dR_l1l2);
//	Tte->SetBranchAddress("dR_b1b2", &dR_b1b2);
//	Tte->SetBranchAddress("dR_bl", &dR_bl);
//	Tte->SetBranchAddress("mass_l1l2", &mass_l1l2);
//	Tte->SetBranchAddress("mass_b1b2", &mass_b1b2);
//	Tte->SetBranchAddress("MT2", &MT2);
//	Tte->SetBranchAddress("weight", &weight);
//	TestTree.Branch("classID", &classID2,"classID/I");
//	TestTree.Branch("className", &className2,"className/C");
//	TestTree.Branch("dR_l1l2", &dR_l1l22,"dR_l1l2/F");
//	TestTree.Branch("dR_b1b2", &dR_b1b22,"dR_b1b2/F");
//	TestTree.Branch("dR_bl", &dR_bl2,"dR_bl/F");
//	TestTree.Branch("mass_l1l2", &mass_l1l22,"mass_l1l2/F");
//	TestTree.Branch("mass_b1b2", &mass_b1b22,"mass_b1b2/F");
//	TestTree.Branch("MT2", &MT22,"MT2/F");
//	TestTree.Branch("weight", &weight2,"weight/F");
//	TestTree.Branch("Likelihood", &Likelihood2,"Likelihood/F");
//	TestTree.Branch("LikelihoodD", &LikelihoodD2,"LikelihoodD/F");
//	TestTree.Branch("LikelihoodMIX", &LikelihoodMIX2,"LikelihoodMIX/F");
//	TestTree.Branch("BDT", &BDT2,"BDT/F");
//	TestTree.Branch("HMatrix", &HMatrix2,"HMatrix/F");
//	if(is_RuleFit) TestTree.Branch("RuleFit", &RuleFit2,"RuleFit/F");
//	if(is_MLP)     TestTree.Branch("MLP", &MLP2,"MLP/F");
//	int TclassID, TclassID2;
//	char TclassName, TclassName2;
//	float TdR_l1l2, TdR_b1b2, TdR_bl, Tmass_l1l2, Tmass_b1b2, TMT2, Tweight;
//	float TdR_l1l22, TdR_b1b22, TdR_bl2, Tmass_l1l22, Tmass_b1b22, TMT22, Tweight2;
//	float TLikelihood2, TLikelihoodD2, TLikelihoodMIX2, TBDT2, THMatrix2, TRuleFit2, TMLP2;
//	Ttt->SetBranchAddress("classID", &TclassID);
//	Ttt->SetBranchAddress("className", &TclassName);
//	Ttt->SetBranchAddress("dR_l1l2", &TdR_l1l2);
//	Ttt->SetBranchAddress("dR_b1b2", &TdR_b1b2);
//	Ttt->SetBranchAddress("dR_bl", &TdR_bl);
//	Ttt->SetBranchAddress("mass_l1l2", &Tmass_l1l2);
//	Ttt->SetBranchAddress("mass_b1b2", &Tmass_b1b2);
//	Ttt->SetBranchAddress("MT2", &TMT2);
//	Ttt->SetBranchAddress("weight", &Tweight);
//	TrainTree.Branch("classID", &TclassID2,"classID/I");
//	TrainTree.Branch("className", &TclassName2,"className/C");
//	TrainTree.Branch("dR_l1l2", &TdR_l1l22,"dR_l1l2/F");
//	TrainTree.Branch("dR_b1b2", &TdR_b1b22,"dR_b1b2/F");
//	TrainTree.Branch("dR_bl", &TdR_bl2,"dR_bl/F");
//	TrainTree.Branch("mass_l1l2", &Tmass_l1l22,"mass_l1l2/F");
//	TrainTree.Branch("mass_b1b2", &Tmass_b1b22,"mass_b1b2/F");
//	TrainTree.Branch("MT2", &TMT22,"MT2/F");
//	TrainTree.Branch("weight", &Tweight2,"weight/F");
//	TrainTree.Branch("Likelihood", &TLikelihood2,"Likelihood/F");
//	TrainTree.Branch("LikelihoodD", &TLikelihoodD2,"LikelihoodD/F");
//	TrainTree.Branch("LikelihoodMIX", &TLikelihoodMIX2,"LikelihoodMIX/F");
//	TrainTree.Branch("BDT", &TBDT2,"BDT/F");
//	TrainTree.Branch("HMatrix", &THMatrix2,"HMatrix/F");
//	if(is_RuleFit) TestTree.Branch("RuleFit", &TRuleFit2,"RuleFit/F");
//	if(is_MLP)     TestTree.Branch("MLP", &TMLP2,"MLP/F");
//	Int_t nentries = (Int_t)Tte->GetEntries();
//	for (Int_t iT=0; iT<nentries; iT++) {
//	  Tte->GetEntry(iT);
//	  classID2 = classID;
//	  className2 = className;
//	  dR_l1l22 = dR_l1l2;
//	  dR_b1b22 = dR_b1b2;
//	  dR_bl2 = dR_bl;
//	  mass_l1l22 = mass_l1l2;
//	  mass_b1b22 = mass_b1b2;
//	  MT22 = MT2;
//	  weight2 = weight;
//	  Likelihood2 = MyTest[0][iT];
//	  LikelihoodD2 = MyTest[1][iT];
//	  LikelihoodMIX2 = MyTest[2][iT];
//	  BDT2 = MyTest[3][iT];
//	  HMatrix2 = MyTest[4][iT];
//	  if(is_RuleFit) RuleFit2 = MyTest[5][iT];
//	  if(is_MLP)     MLP2 = MyTest[is_RuleFit ? 6:5][iT];
//	  //TestTree.Fill();
//	}
//	nentries = (Int_t)Ttt->GetEntries();
//	for (Int_t iT=0; iT<nentries; iT++) {
//	  Ttt->GetEntry(iT);
//	  TclassID2 = TclassID;
//	  TclassName2 = TclassName;
//	  TdR_l1l22 = TdR_l1l2;
//	  TdR_b1b22 = TdR_b1b2;
//	  TdR_bl2 = TdR_bl;
//	  Tmass_l1l22 = Tmass_l1l2;
//	  Tmass_b1b22 = Tmass_b1b2;
//	  TMT22 = TMT2;
//	  Tweight2 = Tweight;
//	  TLikelihood2 = MyTrain[0][iT];
//	  TLikelihoodD2 = MyTrain[1][iT];
//	  TLikelihoodMIX2 = MyTrain[2][iT];
//	  TBDT2 = MyTrain[3][iT];
//	  THMatrix2 = MyTrain[4][iT];
//	  if(is_RuleFit) TRuleFit2 = MyTrain[5][iT];
//	  if(is_MLP)     TMLP2 = MyTrain[is_RuleFit ? 6:5][iT];
//	  //TrainTree.Fill();
//	}
    }//i==Last
cout<<"AA1"<<endl;
  }//All Files
cout<<"AA2"<<endl;
  //TestTree.Write();
  //TrainTree.Write();
  outFile->cd();
  outFile->Close();
cout<<"AA3"<<endl;
}

void CopyDir( TFile * thisfile, TFile * outFile, TString CommonDir ){
  outFile->cd();
  TDirectory *plotDirectory = thisfile->GetDirectory(CommonDir);
  TDirectory *dir1;
  if (!outFile->GetDirectory(CommonDir)) dir1 = outFile->mkdir(CommonDir);  
  else dir1 = outFile->GetDirectory(CommonDir);
  dir1->cd();
  TIter listPlots(plotDirectory->GetListOfKeys());
  TKey *key;
  cout<<"Loop on "<<CommonDir<<endl;
  while((key = (TKey*)listPlots())) {
    TString typeDirectory = key->GetClassName();
    TString nameDirectory = key->GetTitle();
    if(typeDirectory == "TDirectoryFile") {
	TDirectory *plotFinalDirectory = thisfile->GetDirectory(CommonDir+"/"+nameDirectory);
	dir1->cd();
	TDirectory *dir2;
	if (!outFile->GetDirectory(CommonDir+"/"+nameDirectory)) dir2 = dir1->mkdir(nameDirectory);
	else dir2 = dir1;
	dir2->cd();
	TIter listPlots2(plotFinalDirectory->GetListOfKeys());
	TKey *key2;
	while((key2 = (TKey*)listPlots2())) {
	  TString typeDirectory2 = key2->GetClassName();
	  TString nameDirectory2 = key2->GetTitle();
	  if( typeDirectory2 != "TDirectoryFile" ) {
	    dir2->cd();
	    TObject *obj2 = key2->ReadObj();
	    obj2->Write();
	    outFile->cd();
	  }
	  else{
	    TDirectory *plotFinalDirectory2 = thisfile->GetDirectory(CommonDir+"/"+nameDirectory+"/"+nameDirectory2);
	    dir2->cd();
	    TDirectory *dir3;
	    if (!outFile->GetDirectory(CommonDir+"/"+nameDirectory+"/"+nameDirectory2)) dir3 = dir2->mkdir(nameDirectory2);
	    else dir3 = dir2;
	    dir3->cd();
	    TIter listPlots3(plotFinalDirectory2->GetListOfKeys());
	    TKey *key3;
	    while((key3 = (TKey*)listPlots3())) {
		TString typeDirectory3 = key3->GetClassName();
		TString nameDirectory3 = key3->GetTitle();
		if( typeDirectory3 != "TDirectoryFile" ) {
		  dir3->cd();
		  TObject *obj3 = key3->ReadObj();
		  obj3->Write();
		  outFile->cd();
		}
		else{
		  cout<<"WARNING!!!! ANOTHER LAYER!!!"<<endl;
		}
	    }
	  }
	}
    }
    else{
	outFile->cd(CommonDir);
	TObject *obj = key->ReadObj();
	obj->Write();
	outFile->cd("");
    }
  }
}
