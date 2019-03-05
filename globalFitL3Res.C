// File: globalFitL3Res.C
// Created by Mikko Voutilainen, on June 9th, 2014
// Updated on Aug 12, 2015 (Run 2 global fit, patch hard-coded 0.98)
// Purpose: Use JEC combination file to fit L3Res globally
//          including scale uncertainties, FSR eigenvectors etc.
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TF1.h"
#include "TMultiGraph.h"
#include "TLatex.h"
#include "TFitter.h"
#include "TMinuit.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMath.h"
#include "TLine.h"

//#include "tools.h"
#include "tdrstyle_mod14.C"

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <functional>
#include <iterator>

//jec-fit-protype adds
#include "FitBase.hpp"
//#include "MultijetBinnedSum.hpp"
//#include "MultijetCrawlingBins.hpp"
#include "JetCorrDefinitions.hpp"
#include <list>
#include <Minuit2/Minuit2Minimizer.h>
#include <Math/Functor.h>

using namespace std;
int iMultijet =0;

bool dofsr = true; // correct for FSR
double ptreco_gjet = 15.; // min jet pT when evaluating alphamax for gamma+jet
double ptreco_zjet = 5.; // same for Z+jet
bool dol1bias = false; // correct MPF for L1L2L3-L1 (instead of L1L2L3-RC)
bool _paper = false;//true;
bool _useZoom = true;
double _cleanUncert = 0.05; // for eta>2
//double _cleanUncert = 0.020; // Clean out large uncertainty points from PR plot
//bool _g_dcsonly = false;
string scalingForL2OrL3Fit = "None";// was "ApplyL3ResDontScaleDijets";
//"None" - for inpunt combination files without any residual applied
//"PutBackL2Res" - put L2res back in for gamma/Z+jet for vs eta studies
//"ApplyL3ResDontScaleDijets" - apply barrel JES (use case: check closure when only L2Res is applied to the inputs and L3Res didn't change)
//N.B.: Barrel JES from input text file is always applied to dijet results

bool useNewMultijet = false;//true; //use MultijetCrawlingBins now
bool verboseGF = false;

unsigned int _nsamples(0);
unsigned int _nmethods(0);

int cnt(0); int Nk(0);
TF1 *_jesFit(0);
vector<TGraphErrors*> *_vdt(0);
vector<TGraphErrors*> *_vpt(0);
vector<TGraphErrors*> *_vdt2(0);
vector<TGraphErrors*> *_vpt2(0);
vector<TGraphErrors*> *_vdt3(0);
vector<TH1D*> *_vsrc;
void jesFitter(Int_t &npar, Double_t *grad, Double_t &chi2, Double_t *par,
	       Int_t flag);

// Helper functions to draw fit uncertainty band for arbitrary TF1
TF1 *_fitError_func(0);
TMatrixD *_fitError_emat(0);
Double_t fitError(Double_t *xx, Double_t *p);

// Future improvements:
// 1) add more parameters, for SPRE, SPRH, offset, tracker, ECAL *but*
//   also add a priori Gaussian constraints on these parameters based on
//   these parameters (e.g. +/-sqrt(2)*3% for SPRE, SPRH, +/-1% for ECAL etc.)
// => this should allow for more flexible shape, with easy factorization
//    of sources without huge correlations between parameters
// 2) add JER uncertainty for multijets

// Alternative parameterizations
bool useOff = true;//false; // pT-dependent offset
bool useTDI = false; // tracker dynamic inefficiency for 3p fit
//double fixTDI = 1; // fix TDI for BCD+EF and G+H
double fixTDI = 0; // do NOT fix TDI for BCD+EF and G+H
bool useEG = false; // ECAL gain shift for 3p fit

//const int njesFit = 1; // scale only
const int njesFit = 2; // scale(ECAL)+HB
//const int njesFit = 3; //useOff=true; // scale(ECAL)+HB+offset
//const int njesFit = 3; //useTDI=true; // scale(ECAL)+HB+tracker (switchable)
//const int njesFit = 3; //useEG=true; // scale(ECAL)+HB+ECALgain
//const int njesFit = 4; // scale(ECAL)+HB+offset+ECALgain

std::vector<double> parTransForMultiJet(njesFit, 0.001);
//Nuisances DummyMultijetNuisances(NuisanceDefinitions);
//new multijet sources
int nsrcnmj= 0;// dummy value to be overwritten automatically

//std::unique_ptr<MeasurementBase> _lossFunc;
//MeasurementBaseCombLossFunction* _lossFunc=0;
//CombLossFunction* _lossFunc=0;
//CombLossFunction* _inputFunc=0;
//NuisanceDefinitions* nuisanceDefs_ =0;
//Nuisances* nuisances_ =0;
const double ptminJesFit = 30;
TF1 *fhb(0), *fl1(0), *fl1mc(0), *ftr(0), *feg(0); double _etamin(0);
Double_t jesFit(Double_t *x, Double_t *p) {
  
  double pt = *x;

  // Initialize SinglePionHCAL shape
  if (!fhb) fhb = new TF1("fhb","max(0.,[0]+[1]*pow(x,[2]))",10,3500);
  fhb->SetParameters(1.03091e+00, -5.11540e-02, -1.54227e-01); // SPRH
 
  // Initialize L1FastJet-L1RC difference
  // Values from fitting ratio/eta00-13/hl1bias (JEC set in reprocess.C)
  //if (!fl1) fl1 = new TF1("fl1","1-([0]+[1]*log(x))/x",10,3500);
  if (!fl1) fl1 = new TF1("fl1","1-([0]+[1]*log(x)+[2]*pow(log(x),2))/x",
			  10,3500);
  //fl1->SetParameters(2.60382e-01, 1.96664e-01); // Sum16V6G hl1bias
  //fl1->SetParameters(5.71298e-01, 1.59635e-01);
  //fl1->SetParameters(-1.96332e-01, 3.07378e-01); // Sum16_07AugBCDEFGHV6 hl1bias
  fl1->SetParameters(3.906, -1.652, 0.2257); // Sum16_07AugBCDEFGHV16 hl1bias
  
  if (!fl1mc) fl1mc = new TF1("fl1mc","1-([0]+[1]*log(x))/x",10,3500);
  fl1mc->SetParameters(1.414,-0.131); // Sum16_07AugBCDEFGHV16 hl1diff

  // Initialize tracker inefficiency shape
  // Values from drawAvsB.C for EvsG
  // The p3 is turned off for flatter low pT and steeper high pT, which
  // is more consistent with multijet data (not used in the fit, yet)
  //if (!ftr) ftr = new TF1("ftr","1-[0]-[1]*pow(x,[2]) + ([3]+[4]*log(x))/x",10,3500);
  //ftr->SetParameters(-0.04432, 1.304, -0.4624, 0, 1.724);

  // Sigmoid from drawAvsB.C for (BCD+EF)vs(G+H) 2016 Legacy (add -1)
  if (!ftr) ftr = new TF1("ftr","[0]+(1-[0])/(1. + exp(-(log(x)-[1])/[2]))",
			  30,3500);
  //ftr->SetParameters(0.9763, 5.04, 0.3695); // ENDCAP
  //ftr->SetParameters(0.979, 5.030, 0.395); // BCD/G+H
  //ftr->SetParameters(0.980, 5.055, 0.324); // BCD/GH
  //ftr->SetParameters(0.979, 5.036, 0.391); // BCD/GH g180; p1,p2 fix BCD+EF/GH
  ftr->SetParameters(0.982, 5.036, 0.391); // BCD/GH g100; p1,p2 fix BCD+EF/GH

  //if (!feg) feg = new TF1("feg","[0]*TMath::Gaus(x,[1],[2]*sqrt(x))",
  if (!feg) feg = new TF1("feg","[0]+[1]*log(x)+"
			  "[2]*TMath::Gaus(x,[3],[4]*sqrt(x))",
			  10,3500);
  //feg->SetParameters(-1.45, 0.17, -1.4, 366, 13.2);
  //feg->SetParameters(1.86, -0.17, -1.75, 410, 14.0); // Sep23V1 RunG SMP-J
  feg->SetParameters(0,0, -1.75, 410, 14.0); // Sep23V1 RunG SMP-J, bump only

  // As it is in L2L3Residuals:
  //[8]*((1+0.04432-1.304*pow(max(30,min(6500,x)),-0.4624)+(0+1.724*TMath::Log(max(30,min(6500,x))))/x)-(1+0.04432-1.304*pow(208.,-0.4624)+(0+1.724*TMath::Log(208.))/x))

  // Just fitting a constant scale factor for RelativePt uncertainty
  if (njesFit==1) {
    return p[0];
  }

  // Directly using fitted SPR HCAL shape (from JECUncertainty.cpp)
  if (njesFit==2 && !fixTDI) {

    // p[0]: overall scale shift, p[1]: HCAL shift in % (full band +3%)
    double ptref = 208; // pT that minimizes correlation in p[0] and p[1]
    //double ptref = 265; // For Run II, 208.*sqrt(13/8)
    // 208: 0.10, 338: 0.19, 265: 0.15

    return (p[0] + p[1]/3.*100*(fhb->Eval(pt)-fhb->Eval(ptref)));
  } // njesFit==2

  if (njesFit==3 && useOff) {

    // p[0]: overall scale shift, p[1]: HCAL shift in % (full band +3%)
    // p[2]: L1FastJet-L1RC difference
    double ptref = 208; // pT that minimizes correlation in p[0] and p[1]
    
    return (p[0] + p[1]/3.*100*(fhb->Eval(pt)-fhb->Eval(ptref))
	    + p[2]*(fl1->Eval(pt)-fl1->Eval(ptref)));
  } // njesFit==3 && useOff

  if (njesFit==3 && useTDI) {

    // p[0]: overall scale shift, p[1]: HCAL shift in % (full band +3%)
    double ptref = 208; // pT that minimizes correlation in p[0] and p[1]
    
    return (p[0] + p[1]/3.*100*(fhb->Eval(pt)-fhb->Eval(ptref))
	    //+ p[2]*(ftr->Eval(pt)-ftr->Eval(ptref)));
	    + p[2]*(ftr->Eval(pt)-1));
  } // njesFit==3 && TDI

  if (njesFit==2 && fixTDI) {

    double ptref = 208; // pT that minimizes correlation in p[0] and p[1]
    return (p[0] + p[1]/3.*100*(fhb->Eval(pt)-fhb->Eval(ptref))
	    + fixTDI*(ftr->Eval(pt)-1));
  } // njesFit==3 && TDI && fixTDI

  if (njesFit==3 && useEG) {

    // p[0]: overall scale shift, p[1]: HCAL shift in % (full band +3%)
    // p[2]: tracker dynamic inefficiency, p[2]: ECAL gain shift
    double ptref = 208; // pT that minimizes correlation in p[0] and p[1]
    
    return (p[0] + p[1]/3.*100*(fhb->Eval(pt)-fhb->Eval(ptref))
	    //+ 0.01*p[2]*(feg->Eval(pt)-feg->Eval(ptref)));
	    //+ 0.01*p[2]*feg->Eval(pt));
	    + 0.01*p[2]*feg->Eval(pt));
  } // njesFit==3 && EG

  if (njesFit==4 && useOff && useEG) {

    // p[0]: overall scale shift, p[1]: HCAL shift in % (full band +3%)
    // p[2]: tracker dynamic inefficiency, p[2]: ECAL gain shift
    double ptref = 208; // pT that minimizes correlation in p[0] and p[1]
    
    return (p[0] + p[1]/3.*100*(fhb->Eval(pt)-fhb->Eval(ptref))
	    //+ p[2]*(ftr->Eval(pt)-ftr->Eval(ptref))
	    + p[2]*(fl1->Eval(pt)-fl1->Eval(ptref))
	    + 0.01*p[3]*feg->Eval(pt));
  } // njesFit==4
  if (njesFit==4 && useOff && !useEG) {

    // p[0]: overall scale shift, p[1]: HCAL shift in % (full band +3%)
    // p[2]: L1FastJet-L1RC difference, p[3]: L1Data-L1MC difference
    double ptref = 208;
    
    return (p[0] + p[1]/3.*100*(fhb->Eval(pt)-fhb->Eval(ptref))
	    + p[2]*(fl1->Eval(pt)-fl1->Eval(ptref))
	    + p[3]*(fl1mc->Eval(pt)-fl1mc->Eval(ptref)));
  } // njesFit==4 && useOff
  assert(0);
  exit(0);
}

void globalFitL3Res(double etamin = 0, double etamax = 1.3,
		    string epoch="", string selectSample="Standard_MJDJ_gam_zee_zmm", string selectMethods="PtBalMPF") {
  if(verboseGF)cout << Form("Running globalFitL3Res(etamin=%02.2f,etamax=%02.2f,epoch=%s, selectSample=%s, selectMethods=%s",etamin,etamax,epoch.c_str(),selectSample.c_str(),selectMethods.c_str()) << endl << flush;
  _etamin = etamin;
  const char *cep = epoch.c_str();
  //njesFit = (njesFit==3 && useTDI && (epoch=="G"||epoch=="H") ? 2 : njesFit);
  iMultijet=0;
//  if(_lossFunc!=0){
//    delete _lossFunc;
//  }
//  if(_inputFunc!=0){
//    delete _inputFunc;
//  }
  //  auto jetCorr2 = make_unique<JetCorrStd2P>();
  //  auto jetCorr3 = make_unique<JetCorrStd3P>();

  //  auto nuisanceDefs = make_unique<NuisanceDefinitions>;
  //  if(nuisanceDefs_!=0){
  //    delete nuisanceDefs_;
  //  }
  //  NuisanceDefinitions nuisanceDefs;//N.B.: assign _lossFunc only once the measurements and corrsponding nuisances have been registered
  //  nuisanceDefs_ = &nuisanceDefs;

  //  cout << "going to use _lossFunc" << _lossFunc << " " << _lossFunc->GetNumParams()<< endl;
  
  // This is for drawing multijet in _raw.pdf and _orig.pdf around input JEC,
  // to validate the global fit steps (kFSR, shifting by nuisances)
  //  NuisanceDefinitions nuisanceDefsDummy;
  //  _inputFunc = new CombLossFunction(move(jetCorr3),nuisanceDefsDummy);

  // // Andrey Popov, March 19, 2018
  // // https://indico.cern.ch/event/713034/#4-residuals-with-multijet-2016
  // Andrey Popov, April 10, 2018 (v2 of above)
  // https://indico.cern.ch/event/720429/#7-unhealthy-high-pt-electrons
  map<string,const char*> fm_files;
  //dummy files, temp
  fm_files["A"] = "All";
  fm_files["B"] = "All";
  fm_files["C"] = "All";
  fm_files["ABC"] = "All";

//  list<unique_ptr<MeasurementBase>> measurements;//multijet_180507_2016All.root
//  
//  MultijetCrawlingBins* MPF_MJ = new MultijetCrawlingBins(Form("rootfiles/multijet_181217b_2016%s.root", fm_files[cep]), MultijetCrawlingBins::Method::MPF, nuisanceDefs, {"JER"});
//  // JER uncertainty is not considered as the corresponding L2Res variations are buggy (from jec-fit-prototype documentation)
//  MPF_MJ->SetPtLeadRange(0.,1600.);
//  //  measurements.emplace_back(MPF_MJ);

//  MultijetBinnedSum* MPF_MJ = new MultijetBinnedSum(Form("rootfiles/multijet_180507_2016%s.root", fm_files[cep]), MultijetBinnedSum::Method::MPF );
//  if(dropFirstXNewMultijetTriggerBins>0)MPF_MJ->SetTriggerBinRange(dropFirstXNewMultijetTriggerBins);
//  measurements.emplace_back(MPF_MJ);
//  MultijetBinnedSum* MJB_MJ = new MultijetBinnedSum(Form("rootfiles/multijet_180507_2016%s.root", fm_files[cep]), MultijetBinnedSum::Method::PtBal );
//  if(dropFirstXNewMultijetTriggerBins>0)MJB_MJ->SetTriggerBinRange(dropFirstXNewMultijetTriggerBins);
//  measurements.emplace_back(MJB_MJ);



//  if(njesFit==1){
//    _lossFunc = new CombLossFunction(move(jetCorr2),nuisanceDefs);
//    std::cout << "FORCING TO DEACTIVATE NEW MULTIJET: 'njesFit==1' not supported right now" << std::endl;
//    useNewMultijet=false;
//  }
//  if(njesFit==2){
//    _lossFunc = new CombLossFunction(move(jetCorr2),nuisanceDefs);
//    //    std::cout << "FORCING TO DEACTIVATE NEW MULTIJET: Before reactivation make sure that parametrization is compatible" << std::endl;
//    //    useNewMultijet=false;
//  }    
//  //else if(njesFit==3 && useOff){
//  else if(njesFit==3 && useOff){
//    Double_t temp_x;
//    std::vector<Double_t> temp_p = {1.0,0.0,0.0};
//    jesFit(&temp_x,&temp_p[0]);
//    jetCorr3->SetParamsL1({fl1->GetParameter(0),fl1->GetParameter(1)}); //fl1
//    _lossFunc = new CombLossFunction(move(jetCorr3),nuisanceDefs);
//    std::cout << "FORCING TO DEACTIVATE NEW MULTIJET: Before reactivation make sure that parametrization is compatible" << std::endl;
//    useNewMultijet=false;
//    assert(0);
//  }
//  else if (njesFit!=1) {
//    cout << "not defined configuration for multijet..." << endl;
//    useNewMultijet=false;
//    assert(0);
//  }



//  for (auto const &measurement: measurements)
//    _lossFunc->AddMeasurement(measurement.get());
//  //  nsrcnmj = _lossFunc->GetNuisances().GetNumParams();
//  nsrcnmj = nuisanceDefs.GetNumParams();
//
//  std::cout << nuisanceDefs.GetNumParams() << std::endl;
//
//  for ( unsigned i=0; i!=nuisanceDefs.GetNumParams(); i++){
//    std::cout << nuisanceDefs.GetName(i).c_str() << std::endl;
//  }
  
  TDirectory *curdir = gDirectory;
  setTDRStyle();
  writeExtraText = true; // for JEC paper CWR

  TFile *f = new TFile(Form("rootfiles/jecdata%s.root",cep),"READ");
  assert(f && !f->IsZombie());

  // Patch photon E/p (moved to reprocess.C in v7)
  //TFile *feop = new TFile("rootfiles/EoverP_dataMCratio/EoverP_vsRegrCorrEnergy_dataMCRatio.root","READ");
  //assert(feop && !feop->IsZombie());
  //TH1F *heop = (TH1F*)feop->Get("EoverP_vsRegrCorrEnergy_dataMCRatio");
  //assert(heop);
  //TH1F *heope = (TH1F*)heop->Clone("heope");
  //for (int i = 1; i != heop->GetNbinsX()+1; ++i) {
  //heope->SetBinContent(i, heop->GetBinError(i));
  //}
  //TF1 *frot = new TF1("frot","[0]+[1]*pow(x,[2])",10,6500);
  //frot->SetParameters(1.014,-0.105,-0.4);

  // Settings
  const double alpha = 0.30;
  const double ptmj = 30.;
  const char *ct = "ratio";
  //
  map<string, std::vector<string> > methodsmap;
  methodsmap["PtBalMPF"] = {"ptchs","mpfchs1"};
  methodsmap["MPF"] = {"mpfchs1"};
  methodsmap["PtBal"] = {"ptchs"};

  assert(methodsmap.find(selectMethods)!=methodsmap.end());
  const unsigned int nmethods = methodsmap[selectMethods].size();
  vector<const char*> methodsvec;
  for(unsigned int i = 0; i < nmethods; ++i)methodsvec.push_back(methodsmap[selectMethods].at(i).c_str());
  const char** methods = &methodsvec.front();
  
  //  const int nmethods = 2;
  //  const char* methods[nmethods] = {"ptchs","mpfchs1"};
  //const int nmethods = 1;//MPFOnlyTest for FineEtaBins
  //const char* methods[nmethods] = {"mpfchs1"};
  //const char* methods[nmethods] = {"ptchs"};
  _nmethods = nmethods; // for multijets in global fit

  // Global fit with multijets, gamma+jet, Z+jet

  bool isl3 = (etamin==0 && ((epoch!="L4" && fabs(etamax-1.3)<0.1) ||
			     (epoch=="L4" && fabs(etamax-2.4)<0.1)));

  
  map<string, std::vector<string> > samplesmap;
  map<string, int > nsample0map;
  map<string, int > igjmap;
  map<string, int > izeemap;
  map<string, int > izmmmap;
  map<string, int > izllmap;
  // Normal global fit with all four samples (multijet/dijet, gamma+jet, Z+jets)
  samplesmap["Standard_MJDJ_gam_zee_zmm"] =   {(isl3 ? "multijet" : "dijet"),
					     "gamjet", "zeejet", "zmmjet"};
  nsample0map["Standard_MJDJ_gam_zee_zmm"] = 1;
  igjmap["Standard_MJDJ_gam_zee_zmm"] = 0;
  izeemap["Standard_MJDJ_gam_zee_zmm"] = 1;
  izmmmap["Standard_MJDJ_gam_zee_zmm"] = 2;
  izllmap["Standard_MJDJ_gam_zee_zmm"] = -1;

  // Global fit with only dijet, Z+jets
  samplesmap["DJ_zee_zmm"] =   {"dijet", "zeejet", "zmmjet"};
  nsample0map["DJ_zee_zmm"] = 1;
  igjmap["DJ_zee_zmm"] = -1;
  izeemap["DJ_zee_zmm"] = 0;
  izmmmap["DJ_zee_zmm"] = 1;
  izllmap["DJ_zee_zmm"] = -1;

  // Global fit without multijets/dijets
  samplesmap["gam_zee_zmm"] =   {"gamjet", "zeejet", "zmmjet"};
  nsample0map["gam_zee_zmm"] = 0;
  igjmap["gam_zee_zmm"] = 0;
  izeemap["gam_zee_zmm"] = 1;
  izmmmap["gam_zee_zmm"] = 2;
  izllmap["gam_zee_zmm"] = -1;

  // Global fit with only dijets, merged Z+jet
  samplesmap["MJDJ_zll"] =   {(isl3 ? "multijet" : "dijet"), "zlljet"};
  nsample0map["MJDJ_zll"] = 1;
  igjmap["MJDJ_zll"] = -1;
  izeemap["MJDJ_zll"] = -1;
  izmmmap["MJDJ_zll"] = -1;
  izllmap["MJDJ_zll"] = 0;

  // Global fit without multijets/dijets and with merged Z+jet
  samplesmap["gam_zll"] =   {"gamjet", "zlljet"};
  nsample0map["gam_zll"] = 0;
  igjmap["gam_zll"] = 0;
  izeemap["gam_zll"] = -1;
  izmmmap["gam_zll"] = -1;
  izllmap["gam_zll"] = 1;

  // Global fit with  merged Z+jet only
  samplesmap["zll"] =   { "zlljet"};
  nsample0map["zll"] = 0;
  igjmap["zll"] = -1;
  izeemap["zll"] = -1;
  izmmmap["zll"] = -1;
  izllmap["zll"] = 0;
  
  // Global fit with gamma+jet only
  samplesmap["gam"] =   { "gamjet"};
  nsample0map["gam"] = 0;
  igjmap["gam"] = 0;
  izeemap["gam"] = -1;
  izmmmap["gam"] = -1;
  izllmap["gam"] = -1;
  
  // Global fit with only dijets
  samplesmap["DJ"] =   { "dijet"};
  nsample0map["DJ"] = 1;
  igjmap["DJ"] = -1;
  izeemap["DJ"] = -1;
  izmmmap["DJ"] = -1;
  izllmap["DJ"] = -1;
  
  // Global fit with all samples: multijets/dijets, gamma+jet, merged Z+jet
  samplesmap["MJDJ_gam_zll"] =   { (isl3 ? "multijet" : "dijet"),"gamjet", "zlljet"};
  nsample0map["MJDJ_gam_zll"] = 1;
  igjmap["MJDJ_gam_zll"] = 0;
  izeemap["MJDJ_gam_zll"] = -1;
  izmmmap["MJDJ_gam_zll"] = -1;
  izllmap["MJDJ_gam_zll"] = 1;


  // Global fit without photon+jet
  samplesmap["MJDJ_zee_zmm"] =   {(etamin==0 && (etamax==1.3||etamax==2.4)? "multijet" : "dijet"),
				   "zeejet", "zmmjet"};
  nsample0map["MJDJ_zee_zmm"] = 1;
  igjmap["MJDJ_zee_zmm"] = -1;
  izeemap["MJDJ_zee_zmm"] = 0;
  izmmmap["MJDJ_zee_zmm"] = 1;
  izllmap["MJDJ_zee_zmm"] = -1;
  
  // Global fit with Z+jet only
  samplesmap["zee_zmm"] =   {"zeejet", "zmmjet"};
  nsample0map["zee_zmm"] = 0;
  igjmap["zee_zmm"] = -1;
  izeemap["zee_zmm"] = 0;
  izmmmap["zee_zmm"] = 1;
  izllmap["zee_zmm"] = -1;
  
  // Global fit with Z+jet only
  samplesmap["zmm"] =   {"zmmjet"};
  nsample0map["zmm"] = 0;
  igjmap["zmm"] = -1;
  izeemap["zmm"] = -1;
  izmmmap["zmm"] = 0;
  izllmap["zmm"] = -1;
  
  //  string selectSample="Standard_MJDJ_gam_zee_zmm";
  cout<<"Available global fit sample configs, accessible by selectSample in globalFit-function call:" <<endl;
  for(auto elem : samplesmap){
    for (auto i: elem.second)
      std::cout << i << ", ";
    cout << "... choose key: \"" <<elem.first << "\"\n";
  }
  const int nsamples = samplesmap[selectSample].size();
  vector<const char*> samplevec;
  for(int i = 0; i < nsamples; ++i)samplevec.push_back(samplesmap[selectSample].at(i).c_str());
  const char** samples = &samplevec.front();
  //const int nsample0 = 1; // first Z/gamma+jet sample
  const int nsample0 = nsample0map[selectSample];
  const int igj = igjmap[selectSample];
  const int izee = izeemap[selectSample];
  const int izmm = izmmmap[selectSample];
  const int izll = izllmap[selectSample];


  // old default
  /*
  const int nsamples = 3;
  const int nsample0 = 1; // first Z/gamma+jet sample
  const char* samples[3] = {(etamin==0 ? "multijet" : "dijet",
			    "gamjet", "zmmjet"};
  const int igj = 0;
  const int izee = 1; // TEMP
  const int izmm = 1;
  */
  /*
  const int nsamples = 2;
  const int nsample0 = 0; // first Z/gamma+jet sample
  const char* samples[2] = {"gamjet", "zmmjet"};
  */
  /*
  const int nsamples = 1;
  const int nsample0 = 0; // first Z/gamma+jet sample
  const char* samples[1] = {"gamjet"};
  */
  /*
  const int nsamples = 2;
  const int nsample0 = 1; // first Z/gamma+jet sample
  const char* samples[2] = {(etamin==0 ? "multijet" : "dijet"), "gamjet"};
  */
  //
  _nsamples = nsamples; // for multijets in global fit
  // Global fit without multijets
  //const int nsamples = 3;
  //const int nsample0 = 0; // first Z/gamma+jet sample
  //const char* samples[4] = {"gamjet", "zeejet", "zmmjet"};

  
  string sbin = Form("eta%02.0f-%02.0f",10*etamin,10*etamax);
  const char* bin = sbin.c_str();

  assert(f->cd(ct));
  TDirectory *din1 = f->GetDirectory(ct); assert(din1);
  assert(din1->cd(bin));
  TDirectory *d = din1->GetDirectory(bin); assert(d);
  d->cd();


  ///////////////////////////////////////////////
  // Load and setup data and uncertainty sources
  ///////////////////////////////////////////////

  // Load data points
  vector<TGraphErrors*> gs;
  vector<TGraphErrors*> gs2;
  vector<TGraphErrors*> gs3;
  for (unsigned int imethod = 0; imethod != nmethods; ++imethod) {
    for (int isample = 0; isample != nsamples; ++isample) {
      const char *cm = methods[imethod];
      const char *cs = samples[isample];

      //string s = Form("%s_%s_a%02.0f",cm,cs,alpha*100);
      string s = Form("%s_%s_a%02.0f",cm,cs,
		      string(cs)=="multijet" ? ptmj : alpha*100);
      TGraphErrors *g = (TGraphErrors*)d->Get(s.c_str());
      if (!g) cout << "Graph "<<s<<" not found!" << endl << flush;
      assert(g);

      if (string(cs)=="multijet") {
	g->SetMarkerStyle(imethod==0 ? kOpenTriangleUp : kFullTriangleUp);
      }
      g->SetLineWidth(2);

      // Clean out empty points
      // as well as trigger-biased ones for dijets
      for (int i = g->GetN()-1; i != -1; --i) {
	if (g->GetY()[i]==0 || g->GetEY()[i]==0 ||
	    (string(cs)=="dijet" && g->GetX()[i]<70.))
	  g->RemovePoint(i);
      }
      
      gs.push_back(g);
      gs2.push_back((TGraphErrors*)g->Clone(Form("%s_shifted",g->GetName())));
      gs3.push_back((TGraphErrors*)g->Clone(Form("%s_raw",g->GetName())));
    } // for isample
  } // for imethod

  _vdt = &gs;
  _vdt2 = &gs2;
  _vdt3 = &gs3;
  
  // Load pT fractions for multijet balancing
  vector<TGraphErrors*> gfs, gfs2;
  if (string(samples[0])=="multijet") {

    // Fractions for MJB (pT>30 GeV)
    //string s = "ptf30_multijet_a30";
    //string s = "crecoil_multijet_a30";
    string s = Form("crecoil_multijet_a%1.0f",ptmj);
    TGraphErrors *g = (TGraphErrors*)d->Get(s.c_str());
    if (!g) cout << "Graph "<<s<<" not found!" << endl << flush;
    assert(g);
    
    g->SetLineWidth(2);
    gfs.push_back((TGraphErrors*)g->Clone());
    g->SetMarkerStyle(kOpenTriangleDown);
    gfs2.push_back((TGraphErrors*)g->Clone()); // for MJB
    
    // Fractions for MPF (pT>10 GeV)
    //s = "ptf10_multijet_a30";
    //s = "crecoil_multijet_a10"; // 80XV1rereco
    s = "crecoil_multijet_a15"; // 80XV1rereco+Sum16
    g = (TGraphErrors*)d->Get(s.c_str());
    if (!g) cout << "Graph "<<s<<" not found!" << endl << flush;
    assert(g);

    g->SetLineWidth(2);
    gfs.push_back((TGraphErrors*)g->Clone());
    g->SetMarkerStyle(kFullTriangleDown);
    gfs2.push_back((TGraphErrors*)g->Clone()); // for MPF
  }
  
  _vpt = &gfs;
  _vpt2 = &gfs2;

  // Load reference barrel JES for dijets
  TH1D *hjes0 = (TH1D*)f->Get(epoch=="L4" ? "ratio/eta00-24/hjes" :
			      "ratio/eta00-13/hjes"); assert(hjes0);
  hjes0->SetName("hjes0"); // L2Res only (usually)
  TH1D *herr0 = (TH1D*)f->Get(epoch=="L4" ? "ratio/eta00-24/herr" :
			      "ratio/eta00-13/herr"); assert(herr0);
  herr0->SetName("herr0"); // L2L3Res (usually)

  // Load JEC uncertainty band
  TH1D *herr = (TH1D*)d->Get("herr"); assert(herr);
  TH1D *hrun1 = (TH1D*)d->Get("hrun1"); assert(hrun1);
  TH1D *hjes = (TH1D*)d->Get("hjes"); assert(hjes);
  TH1D *herr_ref = (TH1D*)d->Get("herr_ref"); assert(herr_ref);
  TH1D *herr_noflv = (TH1D*)d->Get("herr_noflv"); assert(herr_noflv);
  TH1D *herr_spr = (TH1D*)d->Get("herr_spr"); assert(herr_spr);
  TH1D *herr_pu = (TH1D*)d->Get("herr_pu"); assert(herr_pu);
  
  // Load FSR (JESref) corrections, and correct data
  vector<TH1D*> hks(nsamples*nmethods,0);
  for (unsigned int imethod = 0; imethod != nmethods; ++imethod) {
    for (int isample = 0; isample != nsamples; ++isample) {

      const char *cm = methods[imethod];
      const char *cs = samples[isample];

      string s = Form("fsr/hkfsr_%s_%s",cm,cs);
      TH1D *h = (TH1D*)d->Get(s.c_str());
      if (!h) cout << "Histo "<<s<<" not found!" << endl << flush;
      assert(h);

      // For correcting L1L2L3-L1 to L1L2L3-RC
      TH1D *hl1 = (TH1D*)d->Get("hl1bias"); assert(hl1);

      int ibm = isample + nsamples*imethod;
      h->SetName(Form("bm%d_%s",(1<<ibm),h->GetName()));

      hks[ibm] = h;

      // Correct data for FSR
      TGraphErrors *g = (*_vdt)[ibm];
      TGraphErrors *g2 = (*_vdt2)[ibm];
      TGraphErrors *g3 = (*_vdt3)[ibm];
      TGraphErrors *gf = (_vpt->size()==nmethods ? (*_vpt)[imethod] : 0);
      TGraphErrors *gf2 = (_vpt2->size()==nmethods ? (*_vpt2)[imethod] : 0);
      for (int i = 0; i != g->GetN(); ++i) {
	
	double pt = g->GetX()[i];
	double r = g->GetY()[i];
	double ptreco = (string(cs)=="zeejet"||string(cs)=="zmmjet"
			 ||string(s)=="zlljet" ?
			 ptreco_zjet : ptreco_gjet);
	double aeff = max(alpha, ptreco/pt);
	double kfsr = (dofsr ? 1./(1+aeff*h->GetBinContent(h->FindBin(pt))):1);
	double l1 = (dol1bias && string(cm)=="mpfchs1" && string(cs)=="gamjet" ?
		     1./hl1->GetBinContent(hl1->FindBin(pt)) : 1);

	double scale = 1.00; // correct out previous L3Res
	double escale = 0; // E/p moved to reprocess.C

	// put L2res back in for gamma/Z+jet for vs eta studies
	if (!(etamin==0 && fabs(etamax-1.3)<0.1) &&
	    string(cs)!="dijet" && string(cs)!="multijet") {
	  double jes = hjes->GetBinContent(hjes->FindBin(pt)); // L2Res only
	  double jesref = hjes0->GetBinContent(hjes0->FindBin(pt)); // barrel
	  // divide by jesref(=1) in case hjes had L2L3Res instead of L2Res
          if(scalingForL2OrL3Fit=="None")scale= 1.0;//no residual input
          else if(scalingForL2OrL3Fit=="PutBackL2Res")scale = jes / jesref;
          else if(scalingForL2OrL3Fit=="ApplyL3ResDontScaleDijets")scale = 1 / jesref;
          else assert(false);
	}
	g->SetPoint(i, pt, scale*l1*r*kfsr);
	g2->SetPoint(i, pt, scale*l1*r*kfsr);
	g3->SetPoint(i, pt, scale*l1*r);
	//
	double ex = g->GetErrorX(i);
	double ey = g->GetErrorY(i);
	g->SetPointError(i, ex, sqrt(ey*ey + escale*escale));
	g2->SetPointError(i, ex, sqrt(ey*ey + escale*escale));
	g3->SetPointError(i, ex, sqrt(ey*ey + escale*escale));

	
	// For multijet, correct instead for reference JES
	if (string(cs)=="multijet" && gf) {
	  double ptref = gf->GetY()[i] * pt;
	  double jesref = herr->GetBinContent(herr->FindBin(ptref));
	  // Don't correct g, we need raw input for global fit
	  //g->SetPoint(i, pt, jesref * r);
	  g2->SetPoint(i, pt, jesref * r);
	  g3->SetPoint(i, pt, jesref * r);
	  double jes = herr->GetBinContent(herr->FindBin(pt));
	  gf2->SetPoint(i, ptref, jes / r);
	}
	// For dijet, multiply by barrel JES
	if (string(cs)=="dijet" && scalingForL2OrL3Fit!="ApplyL3ResDontScaleDijets") {
	  double jesref = herr0->GetBinContent(herr0->FindBin(pt));
	  g->SetPoint(i, pt, r*kfsr * jesref);
	  g2->SetPoint(i, pt, r*kfsr * jesref);
	  g3->SetPoint(i, pt, r * jesref);
	}

      } // for i
      
    } // for isample
  } // for imethod


  // Container for all uncertainty sources (one per nuisance parameter)
  vector<TH1D*> hs;

  // FSR uncertainty sources (fit uncertainty eigenvectors)
  // 3x2x(3-4) (we have FSR for multijets as well?)
  const int neig = 3;
  for (int ieig = 0; ieig != neig; ++ieig) {
    for (unsigned int imethod = 0; imethod != nmethods; ++imethod) {
      for (int isample = 0; isample != nsamples; ++isample) {
	
	const char *cm = methods[imethod];
	const char *cs = samples[isample];

	string s = Form("fsr/hkfsr_%s_%s_eig%d",cm,cs,ieig);
	TH1D *h = (TH1D*)d->Get(s.c_str());
	if (!h) cout << "Histo "<<s<<" not found!" << endl << flush;
	assert(h);

	int ibm = isample + nsamples*imethod;	
	for (int i = 1; i != h->GetNbinsX()+1; ++i) {
	  double pt = h->GetBinCenter(i);
	  double ptreco = (string(cs)=="zeejet"||string(cs)=="zmmjet"
			   ||string(cs)=="zlljet" ?
			 ptreco_zjet : ptreco_gjet);
	  double aeff = max(alpha, ptreco/pt);
	  h->SetBinContent(i, aeff*h->GetBinContent(i));
	}
	h->SetName(Form("bm%d_%s",(1<<ibm),h->GetName()));

	hs.push_back(h);
      } // for isample
    } // for imethod
  } // for ieig

  // Uncertainty sources for e, gamma and mu scale uncertainties
  // (and extra source for multijets?)
  int is(0);
  int is_gj(0), is_zee(0), is_zmm(0), is_zll(0);
  for (unsigned int i = 0; i != _vdt->size(); ++i) {
      
    string s = samples[i%nsamples]; const char *cs = s.c_str();
    unsigned int imethod = i/nsamples;
    string m = methods[imethod];    const char *cm = m.c_str();

    TH1D *h = hs[i]; assert(h);
    TH1D *h2 = (TH1D*)h->Clone(Form("bm%d_inactive_%s_%s_%d",1<<i,cm,cs,i));

    double escale(0);
    //bool gain6(true), gain1(true);
    const int n0 = nsample0;
    const int n1 = nsamples+nsample0;
    if (s=="multijet" && m=="ptchs") { // imethod==0) {
      escale = 1.0; // no constraint on absolute scale
      h2->SetName(Form("bm%d_scale_multijet_%d",(1<<(n0-1) | (1<<(n1-1))), i));
    }
    if (s=="dijet" && m=="ptchs") { //imethod==0) {
      escale = 0.005; // for JER bias and whatnot
      h2->SetName(Form("bm%d_scale_dijet_%02.0f_%d",
		       (1<<(n0-1) | (1<<(n1-1))), escale*1000., i));
    }
    if (s=="gamjet" && m=="ptchs") { //imethod==0) {
      //escale = 0.020; // Legacy2016
      // BCDEGHG 2%: 47.7/56 0.989,0.028, 1%: 48.2/56(0.5369) same,
      // BCDEFGH 0.5%: 49.8/56 0.988,0.027(0.5904)
      // BCDEFGH 0.2%: 54.5/56 0.986,0.023(0.6199) => some tension
      //if (epoch=="GH" || epoch=="BCDEFGH") escale = 0.005;
      escale = 0.005;
      //gain6 = true; gain1 = true;//false;
      // Use same source for both MPF and pT balance
      h2->SetName(Form("bm%d_scale_gamjet_%02.0f_%d",
		       (1<<(n0+igj) | (1<<(n1+igj))), escale*1000., i));
      is = hs.size();
      is_gj = hs.size();
    }
    // Separate photon scale uncertainty for gains 1 and 6 (pT>400 GeV)
    // (turned off for time being)
    //if (s=="gamjet" && m=="mpfchs1") { //imethod==1{
    //escale = 0.005;
    //gain6 = false; gain1 = false;//true;
      // Use same source for both MPF and pT balance 
      //h2->SetName(Form("bm%d_scale2_gamjet_%03.0f_%d",
    //	       (1<<(n0+igj) | (1<<(n1+igj))), escale*10000., i));
    //is = hs.size();
    //is_gj = hs.size();
    //}
    if (s=="zeejet" && m=="ptchs") { //imethod==0) {
      escale = 0.002; // Legacy2016G, Zee mass fit within 0.2% up to 300 GeV
      //escale = 0.0005; // Legacy2016BCDEFGH after minitools/drawZmass.C fit
      // Use same source for both MPF and pT balance
      h2->SetName(Form("bm%d_scale_zeejet_%03.0f_%d",
		       (1<<(n0+izee) | (1<<(n1+izee))), escale*10000, i));
      is_zee = hs.size();
    }
    if (s=="zmmjet" && m=="ptchs") { // imethod==0) {
      escale = 0.002; // Legacy2016G, Zmm mass fit within 0.2% up to 300 GeV
      //escale = 0.0005; // Legacy2016BCDEFGH minitools/drawZmass.C fit
      // Use same source for both MPF and pT balance
      h2->SetName(Form("bm%d_scale_zmmjet_%03.0f_%d",
		       (1<<(n0+izmm) | (1<<(n1+izmm))), escale*10000, i));
      is_zmm = hs.size();
    }
    if (s=="zlljet" && m=="ptchs") {
      escale = 0.0020;
      //escale = 0.0005;
      // Use same source for both MPF and pT balance
      h2->SetName(Form("bm%d_scale_zlljet_%03.0f_%d",
		       (1<<(n0+izll) | (1<<(n1+izll))), escale*10000, i));
      is_zll = hs.size();
    }

    // Same scale uncertainty applies to all pT bins
    // UPDATE: now separated by ECAL gain
    for (int j = 1; j != h2->GetNbinsX()+1; ++j) {
      //if ((gain6 && h2->GetBinCenter(j)<400.) ||
      //  (gain1 && h2->GetBinCenter(j)>=400.)) {
      h2->SetBinContent(j, escale);
      h2->SetBinError(j, escale);
      //}
      //else {
      //h2->SetBinContent(j, 0.);
      //h2->SetBinError(j, 0.);
      //}
    } // for j

    hs.push_back(h2);
  } // for i

  // Create additional sources for MPF uncertainties with e/gamma
  // (one for each sample x method, but all except two are empty)
  for (unsigned int i = 0; i != _vdt->size(); ++i) {

    string s = samples[i%nsamples]; const char *cs = s.c_str();
    string m = methods[i/nsamples]; const char *cm = m.c_str();

    double escale(0);
    if (s=="gamjet" && m=="mpfchs1") { escale = 0.005; }
    if (s=="zeejet" && m=="mpfchs1") { escale = 0.005; }
    if (s=="zlljet" && m=="mpfchs1") { escale = 0.002; }

    TH1D *h = hs[i]; assert(h);
    TH1D *h2 = (TH1D*)h->Clone(Form("bm%d_%s_%02.0f_%s_%s_%d",
				    1<<i,
				    escale!=0 ? "mpfscale" : "inactive",
				    escale*1000.,cm,cs,i));
    
    for (int j = 1; j != h2->GetNbinsX()+1; ++j) {
      h2->SetBinContent(j, escale);
      h2->SetBinError(j, escale);
    } // for j
    hs.push_back(h2);
  } // for i

  // *** NEW for 2018 V5M, V5M2 ***
  // Additional uncertainties for pT balance shape so that
  // it is not given more weight than MPF with footprint uncertainty
  // (one for each sample x method, but all except two are empty)
  for (unsigned int i = 0; i != _vdt->size(); ++i) {

    string s = samples[i%nsamples]; const char *cs = s.c_str();
    string m = methods[i/nsamples]; const char *cm = m.c_str();

    double escale(0);
    if (s=="gamjet" && m=="ptchs") { escale = 0.005; }
    if (s=="zeejet" && m=="ptchs") { escale = 0.005; }
    if (s=="zmmjet" && m=="ptchs") { escale = 0.005; }
    if (s=="zlljet" && m=="ptchs") { escale = 0.005; }

    TH1D *h = hs[i]; assert(h);
    TH1D *h2 = (TH1D*)h->Clone(Form("bm%d_%s_%02.0f_%s_%s_%d",
				    1<<i,
				    escale!=0 ? "ptbalscale" : "inactive",
				    escale*1000.,cm,cs,i));
    
    for (int j = 1; j != h2->GetNbinsX()+1; ++j) {
      h2->SetBinContent(j, escale);
      h2->SetBinError(j, escale);
    } // for j
    hs.push_back(h2);
  } // for i

  // MPF uncertainties from offset pT dependence and type-I
  /*
  for (unsigned int i = 0; i != _vdt->size(); ++i) {

    TH1D *hl1 = (TH1D*)d->Get("hl1bias"); assert(hl1);
    TH1D *h2 = (TH1D*)hl1->Clone(Form("bm%d_mpftype1_%d",1<<i,i));
    for (int j = 1; j != h2->GetNbinsX()+1; ++j) {
      int n1 = nsamples+nsample0;
      if (i==n1+0) {
	//if (i==n1+1) {
	h2->SetBinContent(j, 1-hl1->GetBinContent(j));
	h2->SetBinError(j, fabs(1-hl1->GetBinContent(j)));
	//h2->SetName(Form("bm%d_mpfscale_%d",(1<<(n1+0)|1<<(n1+1)|1<<(n1+2)),i));
	// New gamma+jet fixed so turn this off for that (Z+jet MPF only)
	h2->SetName(Form("bm%d_mpfscale_%d",(1<<(n1+1)|1<<(n1+2)),i));
      }
      else {
	h2->SetBinContent(j, 0);
	h2->SetBinError(j, 0);
      }
    } // for j
    hs.push_back(h2);
  } // for i
  */

  // With type-I fix, PileUpPtBB systematic now shared by MPF and pTbal
  // We could therefore fit it, but it could bias high pT end too much
  // Current compromise is to include the lower end shift of -9% in L3Res
  // and to keep this as a systematic so fit uses more high pT information
  if (false && njesFit<=2) { // if njesFit==3, this systematic is included in the fit

    for (unsigned int i = 0; i != _vdt->size(); ++i) {
      
      TH1D *hl1 = (TH1D*)d->Get("hl1bias"); assert(hl1);
      TH1D *h2 = (TH1D*)hl1->Clone(Form("bm%d_l1bias_%d",1<<i,i));
      
      for (int j = 1; j != h2->GetNbinsX()+1; ++j) {
	int n0 = nsample0;
	int n1 = nsamples+nsample0;
	if (int(i)==n0+0) {
	  h2->SetBinContent(j, 1-hl1->GetBinContent(j));
	  h2->SetBinError(j, fabs(1-hl1->GetBinContent(j)));
	  // PileUpPt source applies to all samples and both MPF and pT balance
	  h2->SetName(Form("bm%d_l1bias_%d",(1<<(n0+izee)|1<<(n0+izmm)|
					     1<<(n1+izee)|1<<(n1+izmm)),i));
	}
	else {
	  h2->SetBinContent(j, 0);
	  h2->SetBinError(j, 0);
	}
      } // for j
      hs.push_back(h2);
    } // for i
  }

  // New uncertainty sources for multijets
  if (false && string (samples[0])=="multijet") {

    const int isample = 0;
    for (unsigned int imethod = 0; imethod != nmethods; ++imethod) {
	
      const char *cm = methods[imethod];
      const char *cs = "multijet";
      
      const int npt = 24;
      double vpt[npt+1]= {200, 250, 300, 330, 370, 410, 450, 510, 530, 550, 575,
			  600, 650, 700, 800, 900, 1000, 1100, 1200, 1400, 1600,
			  1800, 2000, 2200};
      TH1D *h = new TH1D(Form("hmj_%s",cm),"",npt,vpt);
      
      // Pt10/Pt30 fits with drawMultijet.C
      // Using 10/30 rather than 20/30, because this shows more effect
      // at low pT where MPF and MJB differ for Pt30
      //MJB:
      double vmjb[4] = {1.8054e+06, -3.2847, -1.0056e+09, -4.3954};
      //MPF:
      double vmpf[4] = {1.221e+05, -2.7955, -2.5519e+08, -4.1888};
      TF1 *f1 = new TF1(Form("f1_%s",cm),
			"1+[0]*pow(x,[1])+[2]*pow(x,[3])",200,1600);
      f1->SetParameters(string(cm)=="ptchs" ? vmjb : vmpf);
      for (int i = 1; i != h->GetNbinsX()+1; ++i) {
	h->SetBinContent(i, f1->Eval(h->GetBinCenter(i))-1);
      }
      
      int ibm = isample + nsamples*imethod;
      h->SetName(Form("bm%d_%s",(1<<ibm),h->GetName()));
      
      hs.push_back(h);
    } // for imethod
  } // for ieig
  
  // Uncertainty sources for multijets
  if (true && string(samples[0])=="multijet") {

    const int isample = 0;
    const int nsrcmj = 3;
    for (unsigned int imethod = 0; imethod != nmethods; ++imethod) {
      for (int isrc = 0; isrc != nsrcmj; ++isrc) {
	
	const char *cm = methods[imethod];
	const char *cs = "multijet";
	
	// Don't remember why these sources were stored under fsr/, but ok
	string s = Form("fsr/%s_%s_src%d",cm,cs,isrc);
	TH1D *h = (TH1D*)d->Get(s.c_str());
	if (!h) cout << "Histo "<<s<<" not found!" << endl << flush;
	assert(h);
	
	int ibm = isample + nsamples*imethod;
	h->SetName(Form("bm%d_%s",(1<<ibm),h->GetName()));
	
	hs.push_back(h);
      } // for imethod
    } // for ieig
  }

  // use global pointer to give outside access
  _vsrc = &hs;

  
  /////////////////////////
  // Draw original data  //
  /////////////////////////
  cout << "Draw original data" << endl;

  const int maxpt = 3500;//1600;
  const int minpt = 30;
  TH1D *h = new TH1D("h",";p_{T} (GeV);Jet response (ratio)",
		     maxpt-minpt,minpt,maxpt);
  h->SetMinimum(etamin>=3 ? 0.50 : (etamin>=2.5 ? 0.70 : 0.91));
  h->SetMaximum(etamin>=3 ? 1.75 : (etamin>=2.5 ? 1.45 : 1.15));
  if (_useZoom) {
    //h->SetMinimum(0.97); // GH
    //h->SetMaximum(1.03); // GH
    h->SetMinimum(0.9501); // GH
    h->SetMaximum(1.045); // GH
  }
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  h->DrawClone("AXIS");

  //map<string, const char*> lumimap;
  //lumimap["BCD"] = "Run2016BCD Legacy, 12.9 fb^{-1}";
  //lumimap["BCD"] = "Run2016BCD 12.9 fb^{-1}"; // DPNote
  //lumimap["E"] = "Run2016E re-mAOD, 4.0 fb^{-1}";
  //lumimap["F"] = "Run2016F re-mAOD, 2.8 fb^{-1}";//3.1 fb^{-1}";
  //lumimap["EF"] = "Run2016EF Legacy, 6.8 fb^{-1}";
  //lumimap["EF"] = "Run2016EF 6.8 fb^{-1}"; // DPNote
  //lumimap["G"] = "Run2016fG Legacy, 8.0 fb^{-1}";
  //lumimap["H"] = "Run2016H Legacy, 8.8 fb^{-1}";
  //lumimap["GH"] = "Run2016FGH re-mAOD, 16.8 fb^{-1}";
  //lumimap["GH"] = "Run2016fGH Legacy, 16.8 fb^{-1}";
  //lumimap["GH"] = "Run2016GH 16.8 fb^{-1}"; // DPNote
  //lumimap["BCDEF"] = "Run2016BCDEF re-mAOD, 19.7 fb^{-1}";
  //lumimap["BCDEFGH"] = "Run2016BCDEFGH Legacy, 36.5 fb^{-1}";
  //lumimap["BCDEFGH"] = "Run2016BCDEFGH 36.5 fb^{-1}"; // DPNote
  //lumimap["L4"] = "Run2016BCDEFGH closure, 36.5 fb^{-1}";
  map<string, const char*> lumimap;
  lumimap["A"] = "Run2018A 14.0 fb^{-1}"; //PdmV Analysis TWiki
  lumimap["B"] = "Run2018B 7.1 fb^{-1}"; //PdmV Analysis TWiki
  lumimap["C"] = "Run2018C 6.9 fb^{-1}"; //PdmV Analysis TWiki
  lumimap["D"] = "Run2018D 31.9 fb^{-1}"; //PdmV Analysis TWiki
  lumimap["ABC"] = "Run2018ABC 28.0 fb^{-1}"; //PdmV Analysis TWiki
  lumimap["ABCD"] = "Run2018ABCD 59.9 fb^{-1}"; //PdmV Analysis TWiki

  lumi_13TeV = lumimap[epoch];

  TCanvas *c0 = tdrCanvas("c0",h,4,11,true);
  gPad->SetLogx();
  
  // multijets up/down
  TLegend *legpf = tdrLeg(0.58,0.79,0.88,0.83);
  legpf->AddEntry(gs[0]," ","PL");
  if (_vpt2->size()!=0) legpf->AddEntry((*_vpt2)[0]," ","PL");
  TLegend *legmf = tdrLeg(0.66,0.79,0.96,0.83);
  if(nmethods==2)legmf->AddEntry(gs[nsamples]," ","PL");
  if (_vpt2->size()>1) legmf->AddEntry((*_vpt2)[1]," ","PL");

  TLegend *legp = tdrLeg(0.58,_useZoom ? 0.70 : 0.55,0.88,0.90);
  if( (nmethods==1||nmethods==2) && strcmp(methods[0],"ptchs")  ==0 ){
    legp->SetHeader("p_{T}^{bal}");
    for (int i = 0; i != nsamples; ++i)
      legp->AddEntry(gs[i]," ",i==0 ? "" : "PL");
  }
  map<string, const char*> texlabel;
  texlabel["multijet"] = "Multijet";
  texlabel["dijet"] = "Dijet";
  texlabel["gamjet"] = "#gamma+jet";
  texlabel["zeejet"] = "Zee+jet";
  texlabel["zmmjet"] = "Z#mu#mu+jet";
  texlabel["zlljet"] = "Z+jet";//"Zl^{+}l{-}+jet";

  TLegend *legm = tdrLeg(0.66,_useZoom ? 0.70 : 0.55,0.96,0.90);
  if( (nmethods==1&&strcmp(methods[0],"mpfchs1")  ==0) || (nmethods==2 && strcmp(methods[1],"mpfchs1")  ==0 )){
    legm->SetHeader("MPF");
    for (int i = 0; i != nsamples; ++i)
      legm->AddEntry(gs[i+(strcmp(methods[0],"mpfchs1")==0 ? 0 : nsamples)],texlabel[samples[i]],i==0 ? "" : "PL");
  }
  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  if (etamin==0 && (fabs(etamax-1.3)<0.1 || fabs(etamax-2.4)<0.1)) {
    if (epoch=="L4") {
      if (dofsr) tex->DrawLatex(0.20,0.73,"|#eta|<2.4, #alpha<0.3#rightarrow0");
      else       tex->DrawLatex(0.20,0.73,"|#eta|<2.4, #alpha<0.3");
    } else {
      if (dofsr) tex->DrawLatex(0.20,0.73,"|#eta|<1.3, #alpha<0.3#rightarrow0");
      else       tex->DrawLatex(0.20,0.73,"|#eta|<1.3, #alpha<0.3");
    }
  }
  else {
    assert(dofsr);
    tex->DrawLatex(0.20,0.73,Form("%1.3f#leq|#eta|<%1.3f",etamin,etamax));
  }
  tex->DrawLatex(0.20, _useZoom ?  0.17 : 0.22, "Before global fit");

  hrun1->SetLineWidth(2);
  hrun1->SetLineColor(kCyan+3);
  hrun1->SetLineStyle(kDashed);

  herr->SetLineWidth(2);
  herr->SetLineColor(kCyan+3);//kGray+2);//kRed+1);
  herr->SetLineStyle(kDashed);
  herr->SetFillColor(kCyan+1);//kGray);//kRed-9);

  herr_ref->SetLineWidth(2);
  herr_ref->SetLineColor(kYellow+3);
  herr_ref->SetLineStyle(kDashed);

  if (!_useZoom) {
    hrun1->DrawClone("SAME E5");
    (new TGraph(hrun1))->DrawClone("SAMEL");
  }
  else {
    herr->DrawClone("SAME E5");
    (new TGraph(herr))->DrawClone("SAMEL");
  }
  herr_ref->DrawClone("SAME E5");
  (new TGraph(herr_ref))->DrawClone("SAMEL");
  if (!_useZoom) {  
    hrun1->SetFillStyle(kNone);
    hrun1->DrawClone("SAME E5");
    (new TGraph(hrun1))->DrawClone("SAMEL");
    hrun1->SetFillStyle(1001);
  }
  else {
    herr->SetFillStyle(kNone);
    herr->DrawClone("SAME E5");
    (new TGraph(herr))->DrawClone("SAMEL");
    herr->SetFillStyle(1001);
  }
  herr_ref->SetFillStyle(kNone);
  herr_ref->DrawClone("SAME E5");
  (new TGraph(herr_ref))->DrawClone("SAMEL");
  herr_ref->SetFillStyle(1001);


  for (unsigned int i = 0; i != _vdt->size(); ++i) {

    TGraphErrors *g2 = (*_vdt2)[i];
    
    // Add multijet downward points
    if (string(samples[i%nsamples])=="multijet") {
      int imethod = i/nsamples;
      TGraphErrors *gf2 = (_vpt2->size()==nmethods ? (*_vpt2)[imethod] : 0);
      if (gf2) {
	gf2->SetMarkerColor(kGray+2);
	gf2->SetLineColor(kGray+2);
	gf2->DrawClone("SAMEPz");
      }
    }
    
    g2->DrawClone("SAMEPz");
  }


//  JetCorrStd2P  jetCorr2Dummy;
//  Nuisances DummyMultijetNuisances(nuisanceDefsDummy); 
//  MultijetCrawlingBins* MPFMultijet = new MultijetCrawlingBins(Form("rootfiles/multijet_181217b_2016%s.root", fm_files[cep]), MultijetCrawlingBins::Method::MPF, nuisanceDefsDummy, {"JER","L1Res","L2Res"});
//  MPFMultijet->SetPtLeadRange(0.,1600.);
//  // MPF
//  TH1D MPF_Raw =  MPFMultijet->RecomputeBalanceData(jetCorr2Dummy, DummyMultijetNuisances);
//  TH1D MPF_SimRaw =  MPFMultijet->RecomputeBalanceSim(jetCorr2Dummy, DummyMultijetNuisances);
////  TH1D* MPF_RatioRaw = (TH1D*) MPF_SimRaw.Clone("MPF_RatioRaw");
////  MPF_RatioRaw->Divide(&MPF_Raw);
//  TGraphErrors MPF_RatioRaw = MPFMultijet->ComputeResiduals(jetCorr2Dummy, DummyMultijetNuisances);
//  MPF_RatioRaw.SetMarkerSize(0.5);
//
//
//  MultijetCrawlingBins* MJBMultijet = new MultijetCrawlingBins(Form("rootfiles/multijet_181217b_2016%s.root", fm_files[cep]), MultijetCrawlingBins::Method::PtBal, nuisanceDefsDummy, {"JER","L1Res","L2Res"});
//  MJBMultijet->SetPtLeadRange(0.,1600.);
//  // MJB
//  TH1D MJB_Raw =  MJBMultijet->RecomputeBalanceData(jetCorr2Dummy, DummyMultijetNuisances);
//  TH1D MJB_SimRaw =  MJBMultijet->RecomputeBalanceSim(jetCorr2Dummy, DummyMultijetNuisances);
//  //  TH1D* MJB_RatioRaw = (TH1D*) MJB_SimRaw.Clone("MJB_RatioRaw");
//  //  MJB_RatioRaw->Divide(&MJB_Raw);
//  TGraphErrors MJB_RatioRaw = MJBMultijet->ComputeResiduals(jetCorr2Dummy, DummyMultijetNuisances);
//  MJB_RatioRaw.SetMarkerSize(0.5);
//  MJB_RatioRaw.SetMarkerStyle(kOpenCircle);
//  
//  // Center points around input JES
//
//  //SetPoint
//  for (int i = 0; i != MPF_RatioRaw.GetN(); ++i) {
//    double ptlead = 0;
//    double residual = 0;
//    MPF_RatioRaw.GetPoint(i,ptlead,residual);
//    double jesfit = herr_ref->GetBinContent(herr_ref->FindBin(ptlead));
//    MPF_RatioRaw.SetPoint(i, ptlead, 1/(residual+1) * jesfit);
//    //    cout << residual << " " << ptlead << " " << jesfit << " " << (residual+1) * jesfit << endl;
//  }  
////  for (int i = 1; i != MPF_RatioRaw->GetNbinsX()+1; ++i) {
////    double ptlead = MPF_RatioRaw->GetBinCenter(i);
////    double mpf = MPF_RatioRaw->GetBinContent(i);
////    double jesfit = herr_ref->GetBinContent(herr_ref->FindBin(ptlead));
////    MPF_RatioRaw->SetBinContent(i, mpf * jesfit);
////    cout << mpf << " " << ptlead << " " << jesfit << " " << mpf * jesfit << endl;
////  }
//  for (int i = 0; i != MJB_RatioRaw.GetN(); ++i) {
//    double ptlead = 0;
//    double residual = 0;
//    MJB_RatioRaw.GetPoint(i,ptlead,residual);
//    double jesfit = herr_ref->GetBinContent(herr_ref->FindBin(ptlead));
//    MJB_RatioRaw.SetPoint(i, ptlead, 1/(residual+1) * jesfit);
//    //    cout << residual << " " << ptlead << " " << jesfit << " " << (residual+1) * jesfit << endl;
//  }  
////  for (int i = 1; i != MJB_RatioRaw->GetNbinsX()+1; ++i) {
////    double ptlead = MJB_RatioRaw->GetBinCenter(i);
////    double mjb = MJB_RatioRaw->GetBinContent(i);
////    double jesfit = herr_ref->GetBinContent(herr_ref->FindBin(ptlead));
////    MJB_RatioRaw->SetBinContent(i, mjb * jesfit);
////  }
//  
//  if(useNewMultijet){
//    MPF_RatioRaw.Draw("PSAME");
//    MJB_RatioRaw.Draw("PSAME");
//  
//    legp->AddEntry(&MJB_RatioRaw," ","PL");
//    legm->AddEntry(&MPF_RatioRaw,"Multijet","PL");
//  }

  if (!_useZoom) {
    legp->AddEntry(hrun1," ","");
    legm->AddEntry(hrun1,"Run I","FL");
    legp->AddEntry(herr_ref," ","");
  }
  else legp->AddEntry(herr," ","FL");
  //legm->AddEntry(herr_ref,"07AugV7","FL");
  //legm->AddEntry(herr_ref,"Run II","FL");
  //legm->AddEntry(herr_ref,"V10","FL");
  //legm->AddEntry(herr_ref,"V16M","FL");
  legm->AddEntry(herr_ref,"V5M2","FL");


  ///////////////////////
  // Draw raw response
  //////////////////////
  cout << "Draw raw response" << endl;

  h = (TH1D*)h->Clone("h0b");
  h->SetYTitle("Raw jet response (ratio)");  
  TCanvas *c0b = tdrCanvas("c0b",h,4,11,true);
  gPad->SetLogx();

  legpf->Draw();
  legmf->Draw();
  legp->Draw();
  legm->Draw();

  if (etamin==0 && (fabs(etamax-1.3)<0.1 || fabs(etamax-2.4)<0.1)) {
    if (epoch!="L4") tex->DrawLatex(0.20,0.73,"|#eta|<1.3, #alpha<0.3");
    if (epoch=="L4") tex->DrawLatex(0.20,0.73,"|#eta|<2.4, #alpha<0.3");
  }
  else tex->DrawLatex(0.20,0.73,Form("%1.3f#leq|#eta|<%1.3f",etamin,etamax));

  tex->DrawLatex(0.20, _useZoom ? 0.17 : 0.22, "Before FSR correction");

  if (!_useZoom) {
    hrun1->DrawClone("SAME E5");
    (new TGraph(hrun1))->DrawClone("SAMEL");
  }
  else {
    herr->DrawClone("SAME E5");
    (new TGraph(herr))->DrawClone("SAMEL");
  }

  herr_ref->DrawClone("SAME E5");
  (new TGraph(herr_ref))->DrawClone("SAMEL");

  if (!_useZoom) {
    hrun1->SetFillStyle(kNone);
    hrun1->DrawClone("SAME E5");
    (new TGraph(hrun1))->DrawClone("SAMEL");
    hrun1->SetFillStyle(1001);
  }
  else {
    herr->SetFillStyle(kNone);
    herr->DrawClone("SAME E5");
    (new TGraph(herr))->DrawClone("SAMEL");
    herr->SetFillStyle(1001);
  }

  herr_ref->SetFillStyle(kNone);
  herr_ref->DrawClone("SAME E5");
  (new TGraph(herr_ref))->DrawClone("SAMEL");
  herr_ref->SetFillStyle(1001);

  for (unsigned int i = 0; i != _vdt3->size(); ++i) {
    
    TGraphErrors *g3 = (*_vdt3)[i];
    // Skip multijet downward points
    g3->DrawClone("SAMEPz");
  }

//  if(useNewMultijet){
//    MPF_RatioRaw.Draw("PSAME");
//    MJB_RatioRaw.Draw("PSAME");
//  }
  
  ///////////////////////
  // Perform global fit
  //////////////////////
  cout << "Perform global fit" << endl;
  
  // Fit function
  // 2015-01-17: use only range visible on plots (one g+jet point at >2.5 TeV)
  TF1 *jesfit = new TF1("jesfit",jesFit,minpt,maxpt,njesFit);
  jesfit->SetLineColor(kBlack);
  if (njesFit==1) {
    jesfit->SetParameter(0, 0.98);
  }
  if (njesFit==2) {
    jesfit->SetParameters(0.990, 0.043);
  }
  if (njesFit==3 && useOff) {
    jesfit->SetParameters(0.989, 0.053, -0.370);
    // per era settings to make converge faster with new multijets (EGM1 values)
    //if (epoch=="BCDEFGH") jesfit->SetParameters(0.988, 0.096, -0.338);
    //if (epoch=="BCD")     jesfit->SetParameters(0.983, 0.138, -0.714);
    //if (epoch=="EF")      jesfit->SetParameters(0.972, 0.135, -0.578);
    //if (epoch=="GH")      jesfit->SetParameters(0.989, 0.055, -0.377);
    // per era settings to make converge faster with new multijets (EGM2 values)
    //if (epoch=="BCDEFGH") jesfit->SetParameters(0.985, 0.088, -0.389);
    //if (epoch=="BCD")     jesfit->SetParameters(0.983, 0.138, -0.785);
    //if (epoch=="EF")      jesfit->SetParameters(0.973, 0.133, -0.530);
    //if (epoch=="GH")      jesfit->SetParameters(0.989, 0.053, -0.370);
    // per era settings to make converge faster with new multijets (EGM3 values)
    //if (epoch=="BCDEFGH") jesfit->SetParameters(0.989, 0.095, -0.378);
    //if (epoch=="BCD")     jesfit->SetParameters(0.984, 0.137, -0.798);
    //if (epoch=="EF")      jesfit->SetParameters(0.973, 0.132, -0.608);
    //if (epoch=="GH")      jesfit->SetParameters(0.990, 0.054, -0.413);
    // per era for EGM3 and Sum16_07Aug2017BCDEFGH_V6 L1 bias
    if (epoch=="BCD")     jesfit->SetParameters(0.9845, 0.1803, -1.682);
    if (epoch=="EF")      jesfit->SetParameters(0.9726, 0.1643, -1.172);
    if (epoch=="GH")      jesfit->SetParameters(0.9894, 0.0748, -0.783);
    if (epoch=="BCDEFGH") jesfit->SetParameters(0.9902, 0.1152, -0.868);
  }
  if (njesFit==3 && useTDI) {
    jesfit->SetParameters(0.985, 0.001, 0.5);
  }
  if (njesFit==2 && fixTDI) {
    double ftdi = 0.5;
    //if (epoch=="BCD" || epoch=="EF") fixTDI = 1;
    if (epoch=="BCD") fixTDI = 1;
    //if (epoch=="EF") fixTDI = 1.524;//=(1-0.968)/(1-0.979); // g180
    if (epoch=="EF") fixTDI = 1.500;//=(1-0.968)/(1-0.979); // g100
    if (epoch=="G" || epoch=="H" || epoch=="GH") fixTDI = 0;
    //if (epoch=="BCDEFGH") fixTDI = 19.7/(19.7+16.8);
    //if (epoch=="BCDEFGH") fixTDI = (12.9*1+6.8*1.524)/(19.7+16.8); // g180
    if (epoch=="BCDEFGH") fixTDI = (12.9*1+6.8*1.500)/(19.7+16.8); // g100
    jesfit->SetParameters(0.985, 0.001);
  }

  if (njesFit==3 && useEG) {
    jesfit->SetParameters(0.995, -0.025, 1.0);
  }
  if (njesFit==4 && useOff && useEG) {
    jesfit->SetParameters(0.995, -0.025, 1.0, 1.0);
  }
  if (njesFit==4 && useOff && !useEG) {
    jesfit->SetParameters(0.989, 0.053, -0.370, 0);
    if (epoch=="BCD")     jesfit->SetParameters(0.9845, 0.1803, -1.682, 0);
    if (epoch=="EF")      jesfit->SetParameters(0.9726, 0.1643, -1.172, 0);
    if (epoch=="GH")      jesfit->SetParameters(0.9894, 0.0748, -0.783, 0);
    if (epoch=="BCDEFGH") jesfit->SetParameters(0.9902, 0.1152, -0.868, 0);
  }

  _jesFit = jesfit;
  

  // Get linear equations and solve them to get initial values for fitter
  const int np = _jesFit->GetNpar();
  //const int np = _jesFit->GetNumberFreeParameters();
  const int nsrc = _vsrc->size();
  Int_t Npar = np+nsrc+nsrcnmj;

  cout << "Global fit has " << np << " fit parameters, "
       << nsrc << " nuisance parameters, and " <<nsrcnmj << " 'new multijet' nuisance parameters" << endl;

  vector<double> a(Npar, 0);
  for (int i = 0; i != np; ++i) a[i] = jesfit->GetParameter(i);

  // Setup global chi2 fit (jesFitter is our function)
  TFitter *fitter = new TFitter(Npar);
  fitter->SetFCN(jesFitter);

  // Set parameters
  vector<string> parnames(Npar);
//  for(int i = Npar-nsrcnmj; i<Npar; ++i){
//    parnames.at(i)=nuisanceDefs.GetName(i-(Npar-nsrcnmj)).c_str();
//  }
  for (int i = 0; i != Npar; ++i)
    fitter->SetParameter(i, parnames[i].c_str(), a[i], (i<np ? 0.01 : 1),
			 -100, 100);

  // Run fitter (multiple times if needed)
  const int nFit = 1;
  cnt = 0;
  for (int i = 0; i != nFit; ++i)
    fitter->ExecuteCommand("MINI", 0, 0);
  TMatrixD emat(Npar, Npar);
  gMinuit->mnemat(emat.GetMatrixArray(), Npar);

  // Retrieve the chi2 the hard way
  Double_t tmp_par[Npar], tmp_err[Npar];
  Double_t chi2_gbl(0), chi2_src(0), chi2_data(0);
  vector<double> vchi2_data(gs2.size(),0);
  vector<int> vndata(gs2.size(),0);
  int nsrc_true(0);
  Double_t grad[Npar];
  Int_t flag = 1;
  TH1D *hsrc = new TH1D("hsrc",";Nuisance parameter;",12,-3,3);//30,-3,3);

  for (int i = 0; i != Npar; ++i) {
    tmp_par[i] = fitter->GetParameter(i);
    tmp_err[i] = fitter->GetParError(i);
  }
  jesFitter(Npar, grad, chi2_gbl, tmp_par, flag);
  cout << "List of fitted sources that count (nonzero):" << endl;
  for (int i = np; i != Npar-nsrcnmj; ++i) {
    if (fabs(tmp_par[i])!=0 || fabs(tmp_err[i]-1)>1e-2) {
      ++nsrc_true;
      hsrc->Fill(tmp_par[i]);
    }
    chi2_src += pow(tmp_par[i],2);
  }
  for (unsigned int i = 0; i != gs2.size(); ++i) {
    for (int j = 0; j != gs2[i]->GetN(); ++j) {
      double x = gs2[i]->GetX()[j];
      double y = gs2[i]->GetY()[j];
      double ey = gs2[i]->GetEY()[j];
      chi2_data += pow((y-jesfit->Eval(x))/ey,2);
      vchi2_data[i] += pow((y-jesfit->Eval(x))/ey,2);
      ++vndata[i];
    }
  }

  cout << endl;
  cout << "*** Processing eta bin " << etamin<<" - "<<etamax << " ***" << endl;
  cout << "Global chi2/ndf = " << chi2_gbl
       << " / " << Nk-np << " for " << Nk <<" data points, "
       << np << " fit parameters, "
       << nsrc << " ("<<nsrc_true<<") uncertainty sources, and " << nsrcnmj << "new multijet nuisances" << endl;
  cout << endl;
  //  cout <<"New multijet:" <<endl;
  //  cout << "Dimensions: "<< _lossFunc->GetNDF() + _lossFunc->GetNumParams() << endl;
//  if(useNewMultijet){
//    cout << "chi2(multijet) term: "<< _lossFunc->Eval(parTransForMultiJet,nuisanceDefs) << endl;
//    cout << "mutltijet chi2 retrieved " << iMultijet << " times."<< endl;
//  }
  cout << endl;
  cout << "For data chi2/ndf = " << chi2_data << " / " << Nk << endl;
  cout << "For sources chi2/ndf = " << chi2_src << " / " << nsrc_true << endl;
  cout << "Per data set:" << endl;
  for (unsigned int i = 0; i != vchi2_data.size(); ++i) {
    cout << "  " << Form("%1.1f",vchi2_data[i]) << " / " << vndata[i]
	 << "  ("<<gs2[i]->GetName()<<")"<< endl;
  }

  ofstream txtL2L3("txt2/GlobalFitOutput_L2L3Residuals.txt",ios_base::app);
  if(njesFit==2){
    if(etamin==0.&&etamax==0.261&&np==2)txtL2L3 << "{ 1 JetEta 1 JetPt 1./([0]+[1]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[2]*((-2.36997+0.413917*TMath::Log(x))/x-(-2.36997+0.413917*TMath::Log(208))/208)) Correction L2Relative}";
    if(np==2&& !(etamin==0.&&etamax==1.3)){
      txtL2L3 << Form("\n %7.4f  %7.4f  5 10 6500 %7.4f %7.4f 0.0 ", etamin, etamax, tmp_par[0], tmp_par[1]);
    }
  }
  if(njesFit==1){
    if(etamin==0.&&etamax==0.261&&np==1)txtL2L3 << "{ 1 JetEta 1 JetPt 1./[0] Correction L2Relative}";
    if(np==1&& !(etamin==0.&&etamax==1.3)){
      txtL2L3 << Form("\n %7.4f  %7.4f  5 10 6500 %7.4f 0.0 0.0 ", etamin, etamax, tmp_par[0]);
    }
  }

  ofstream txtChi2_NDF("txt2/GlobalFitOutput_L2L3Residuals_Chi2OverNDF.txt",ios_base::app);
  if(etamin==0.&&etamax==0.261)txtChi2_NDF << "{ 1 JetEta 1 JetPt [0] Correction L2Relative}";
  if(np==2&& !(etamin==0.&&etamax==1.3)){
    txtChi2_NDF << Form("\n %7.4f  %7.4f  5 10 6500 %7.4f 0.0 0.0 ", etamin, etamax, chi2_gbl/(Nk-np));
  }

  for (int isample = 0; isample != nsamples; ++isample) {
    for (unsigned int imethod = 0; imethod != nmethods; ++imethod) {
      const char *cm = methods[imethod];
      const char *cs = samples[isample];
      string s = Form("fsr/hkfsr_%s_%s",cm,cs);
      TH1D *h = (TH1D*)d->Get(s.c_str());
      if (!h) cout << "Histo "<<s<<" not found!" << endl << flush;
      assert(h);

      ofstream txtFSRDiJet(Form("txt2/GlobalFitOutput_FSR_%s_%s.txt",cs,cm),ios_base::app);
      if(etamin==0.&&etamax==0.261)txtFSRDiJet << "{ 2 JetEta  JetPt 1 JetPt [0] Correction L2Relative}";
      for (int i = 1; i != h->GetNbinsX()+1; ++i) {
        //        double pt = h->GetBinCenter(i);
        double ptlow = h->GetBinLowEdge(i);
        double ptup = h->GetBinLowEdge(i+1);
        double FSR = h->GetBinContent(i);
        //        ofstream txtFSRDiJetPtBin(Form("GlobalFitOutput_FSR_%s_%s_%04.0f.txt",cs,cm,ptlow),ios_base::app);
        //        if(etamin==0.)txtFSRDiJetPtBin << "{ 2 JetEta  JetPt 1 JetPt [0] Correction L2Relative}";
        if(np==2&& !(etamin==0.&&etamax==1.3)){
          //          txtFSRDiJetPtBin << Form("\n %7.4f   %7.4f  %7.0f   %7.0f  3 10 6500 %7.4f ", etamin, etamax, ptlow, ptup, FSR*0.3); // dr/dalpha x dalpha with dalpha = 0.3 to be on same page as dijet standalone
          txtFSRDiJet      << Form("\n %7.4f   %7.4f  %7.0f   %7.0f  3 10 6500 %7.4f ", etamin, etamax, ptlow, ptup, FSR*0.3);
        }
      }
    }
  }

  
  cout << "Fit parameters:" << endl;
  for (int i = 0; i != np; ++i) {
    cout << Form("%2d: %9.4f +/- %6.4f,   ",
		 i+1,tmp_par[i],sqrt(emat[i][i]));// << endl;
  }
  cout << endl;
  // For l2l3res.C
  cout << Form("    // eta bin %1.1f - %1.1f", etamin, etamax) << endl;
  cout << "    {";
  for (int i = 0; i != np; ++i) {
    cout << Form("%7.4f%s",tmp_par[i],i==np-1 ? "}," : ",");
  }
  cout << endl;

  cout << "Constant scale @ HCAL=0: ";
  cout << jesfit->Eval(fhb->GetX(1,1,30)) << endl;
  cout << endl;

  cout << "Error matrix:" << endl;
  for (int i = 0; i != np; ++i) {
    for (int j = 0; j != np; ++j) {
      cout << Form("%10.4g ", emat[i][j]);
    } // for j
    cout << endl;
  } // for i
  cout << endl;

  cout << "Error matrix square root:" << endl;
  for (int i = 0; i != np; ++i) {
    for (int j = 0; j != np; ++j) {
      cout << Form("%7.4f ", sqrt(fabs(emat[i][j]))*TMath::Sign(1.,emat[i][j]));
    } // for j
    cout << endl;
  } // for i
  cout << endl;

  cout << "Correlation matrix:" << endl;
  for (int i = 0; i != np; ++i) {
    for (int j = 0; j != np; ++j) {
      cout << Form("%5.2f ", emat[i][j]/sqrt(emat[i][i]*emat[j][j]));
    } // for j
    cout << endl;
  } // for i
  cout << endl;

  cout << "Uncertainty sources:" << endl;
  for (int i = np; i != Npar-nsrcnmj; ++i) {
    //if (tmp_par[i]==0&&fabs(tmp_err[i]-1)<1e-2)
    if (tmp_par[i]==0&&fabs(tmp_err[i]-1)<2e-2)
      cout << Form("%2d - %35s:  ----------\n",
		   (i-np+1), (*_vsrc)[i-np]->GetName());
    else
      cout << Form("%2d - %35s:  %5.2f+/-%4.2f\n",
		   i-np+1, (*_vsrc)[i-np]->GetName(),
		   tmp_par[i],tmp_err[i]);
  }
  cout << endl;
//  cout << "Nuisances new multijet:" << endl;
//  for (int i = Npar-nsrcnmj; i != Npar; ++i) {
//    //    cout << i-(Npar-nsrcnmj) << endl;
//    //    cout << i << endl;
//    if (tmp_par[i]==0&&fabs(tmp_err[i]-1)<2e-2)
//      cout << Form("%2d - %35s:  ----------\n",
//                   (i-np+1), nuisanceDefs.GetName(i-(Npar-nsrcnmj)).c_str()); 
//    else
//      cout << Form("%2d - %35s:  %5.2f+/-%4.2f\n",
//		   i-np+1, nuisanceDefs.GetName(i-(Npar-nsrcnmj)).c_str(), 
//		   tmp_par[i],tmp_err[i]);
//  }
//  cout << endl;



  ///////////////////////
  // Draw shifted data
  ///////////////////////
  cout << "Draw shifted data" <<endl;
  
  h = (TH1D*)h->Clone("h1");
  h->SetYTitle("Post-fit jet response (ratio)");  

  TCanvas *c1 = tdrCanvas("c1",h,4,11,true);
  gPad->SetLogx();

  legpf->Draw();
  legmf->Draw();
  legp->Draw();
  legm->Draw();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(h->GetMinimum(),
	      6500./cosh(etamin),
	      h->GetMinimum()+0.5*(h->GetMaximum()-h->GetMinimum()),
	      6500./cosh(etamin));

  if (etamin==0 && (fabs(etamax-1.3)<0.1 || fabs(etamax-2.4)<0.1)) {
    if (epoch!="L4") tex->DrawLatex(0.20,0.73,"|#eta|<1.3, #alpha<0.3#rightarrow0");
    if (epoch=="L4") tex->DrawLatex(0.20,0.73,"|#eta|<2.4, #alpha<0.3#rightarrow0");
  }
  else tex->DrawLatex(0.20,0.73,Form("%1.3f#leq|#eta|<%1.3f",etamin,etamax));

  if (!_useZoom) tex->DrawLatex(0.20,0.22,"After global fit");
  tex->DrawLatex(0.20,0.17,Form("#chi^{2} / NDF = %1.1f / %d",
				chi2_gbl, Nk-np));

  is += np;

  // Fitted lepton/photon scales too much detailed for the paper
  tex->SetTextSize(0.030); tex->SetTextColor(kBlue-9);

  if (!_useZoom) {
    hrun1->DrawClone("SAME E5");
    (new TGraph(hrun1))->DrawClone("SAMEL");
  }
  else {
    herr->DrawClone("SAME E5");
    (new TGraph(herr))->DrawClone("SAMEL");
  }

  herr_ref->DrawClone("SAME E5");
  (new TGraph(herr_ref))->DrawClone("SAMEL");

  // Return SPR uncertainty to unnormalized
  for (int i = 1; i != herr_spr->GetNbinsX()+1; ++i) {
    double pt = herr_spr->GetBinCenter(i);
    double err = (pt<200 ? -1 : +1)*herr_spr->GetBinError(i);
    double newerr = (err + 0.008*sqrt(2)) / sqrt(2);
    herr_spr->SetBinContent(i, herr_spr->GetBinContent(i)+0.006);
    herr_spr->SetBinError(i, newerr);
  }
  herr_spr->SetFillColor(kCyan+1);
  herr_spr->SetFillStyle(1001);

  if (!_useZoom) {
    hrun1->SetFillStyle(kNone);
    hrun1->DrawClone("SAME E5");
    (new TGraph(hrun1))->DrawClone("SAMEL");
    hrun1->SetFillStyle(1001);
  }
  else {
    herr->SetFillStyle(kNone);
    herr->DrawClone("SAME E5");
    (new TGraph(herr))->DrawClone("SAMEL");
    herr->SetFillStyle(1001);
  }

  herr_ref->SetFillStyle(kNone);
  herr_ref->DrawClone("SAME E5");
  (new TGraph(herr_ref))->DrawClone("SAMEL");
  herr_ref->SetFillStyle(1001);

  _jesFit->SetNpx(1000);
  _jesFit->DrawClone("SAME");
  for (unsigned int i = 0; i != _vdt2->size(); ++i) {
    
    TGraphErrors *g2 = (*_vdt2)[i];

    // Add multijet downward points
    if (string(samples[i%nsamples])=="multijet") {
      unsigned int imethod = i/nsamples;
      TGraphErrors *gf2 = (_vpt2->size()==nmethods ? (*_vpt2)[imethod] : 0);
      if (gf2) {
	gf2->SetMarkerColor(kGray+2);
	gf2->SetLineColor(kGray+2);
	gf2->DrawClone("SAMEPz");
      }
    }
    
    if (string(samples[i%nsamples])=="gamjet")
      g2->SetLineWidth(3);

    // Clean out large uncertainties
    if (_cleanUncert) {
      for (int i = g2->GetN()-1; i != 0; --i) {
	if (g2->GetEY()[i]>_cleanUncert) g2->RemovePoint(i);
      }
    }

    g2->DrawClone("SAMEPz");
  }

  TF1 *fke = new TF1("fitError", fitError, minpt, maxpt, 1);
  _fitError_func = jesfit;
  TMatrixD emat2(np, np);
  for (int i = 0; i != np; ++i) {
    for (int j = 0; j != np; ++j) {
      emat2[i][j] = emat[i][j];
    } // for i
  } // for j
  _fitError_emat = &emat2;

  fke->SetNpx(1000);
  fke->SetLineStyle(kDotted);
  fke->SetLineColor(jesfit->GetLineColor());
  fke->SetParameter(0,-1);
  fke->DrawClone("SAME");
  fke->SetParameter(0,+1);
  fke->DrawClone("SAME");

//  TH1D MPF_PostFit =  MPFMultijet->RecomputeBalanceData(*(_lossFunc->GetCorrector()), *nuisances_);
//  TH1D MPF_SimPostFit =  MPFMultijet->RecomputeBalanceSim(*(_lossFunc->GetCorrector()), *nuisances_);
//  //  TH1D* MPF_RatioPostFit = (TH1D*) MPF_SimPostFit.Clone("MPF_RatioPostFit");
//  //  MPF_RatioPostFit->Divide(&MPF_PostFit);
//  TGraphErrors MPF_RatioPostFit = MPFMultijet->ComputeResiduals(*(_lossFunc->GetCorrector()), *nuisances_);
//  MPF_RatioPostFit.SetMarkerSize(0.5);
//  TGraphErrors *MPF_RatioPostFitScaled = (TGraphErrors*)MPF_RatioPostFit.Clone("NewMultiJet_RatioPostFitScaled");
//  for (int i = 0; i != MPF_RatioPostFitScaled->GetN(); ++i) {
//    double ptlead = 0;
//    double residual = 0;
//    MPF_RatioPostFitScaled->GetPoint(i,ptlead,residual);
//    MPF_RatioPostFitScaled->SetPoint(i, ptlead, 1/(residual+1) * _jesFit->Eval(ptlead));
//  }  
//  if(useNewMultijet)MPF_RatioPostFitScaled->Draw("PSAME");
//
//  TH1D MJB_PostFit =  MJBMultijet->RecomputeBalanceData(*(_lossFunc->GetCorrector()), *nuisances_);
//  TH1D MJB_SimPostFit =  MJBMultijet->RecomputeBalanceSim(*(_lossFunc->GetCorrector()), *nuisances_);
//  //  TH1D* MJB_RatioPostFit = (TH1D*) MJB_SimPostFit.Clone("MJB_RatioPostFit");
//  //  MJB_RatioPostFit->Divide(&MJB_PostFit);
//  TGraphErrors MJB_RatioPostFit = MJBMultijet->ComputeResiduals(*(_lossFunc->GetCorrector()), *nuisances_);
//  MJB_RatioPostFit.SetMarkerSize(0.5);
//  MJB_RatioPostFit.SetMarkerStyle(kOpenCircle);
//  TGraphErrors *MJB_RatioPostFitScaled = (TGraphErrors*)MJB_RatioPostFit.Clone("NewMultiJet_RatioPostFitScaled");
//  // Center points around fitted JES
//  for (int i = 0; i != MJB_RatioPostFitScaled->GetN(); ++i) {
//    double ptlead = 0;
//    double residual = 0;
//    MJB_RatioPostFitScaled->GetPoint(i,ptlead,residual);
//    MJB_RatioPostFitScaled->SetPoint(i, ptlead, 1/(residual+1) * _jesFit->Eval(ptlead));
//  }  
//  if(useNewMultijet)MJB_RatioPostFitScaled->Draw("PSAME");
//
//  TFile *fout = new TFile(Form("rootfiles/jecdata%s_Multijet.root",cep),"RECREATE");
//  MPF_Raw      	  .SetName("MPF_Raw"       );	
//  MPF_RatioRaw    .SetName("MPF_RatioRaw"       );	
//  MPF_SimRaw   	  .SetName("MPF_SimRaw"    );	
//  MPF_PostFit  	  .SetName("MPF_PostFit"	  );
//  MPF_RatioPostFit.SetName("MPF_RatioPostFit"	  );
//  MPF_SimPostFit  .SetName("MPF_SimPostFit");
//
//  MPF_Raw      	.Write();
//  MPF_RatioRaw      	.Write();
//  MPF_SimRaw   	.Write();
//  MPF_PostFit  	.Write();
//  MPF_RatioPostFit  	.Write();
//  MPF_RatioPostFitScaled->Write();
//  MPF_SimPostFit.Write();
//
//  MJB_Raw      	  .SetName("MJB_Raw"       );	
//  MJB_RatioRaw    .SetName("MJB_RatioRaw"       );	
//  MJB_SimRaw   	  .SetName("MJB_SimRaw"    );	
//  MJB_PostFit  	  .SetName("MJB_PostFit"	  );
//  MJB_RatioPostFit.SetName("MJB_RatioPostFit"	  );
//  MJB_SimPostFit  .SetName("MJB_SimPostFit");
//
//  MJB_Raw      	.Write();
//  MJB_RatioRaw      	.Write();
//  MJB_SimRaw   	.Write();
//  MJB_PostFit  	.Write();
//  MJB_RatioPostFit  	.Write();
//  MJB_RatioPostFitScaled->Write();
//  MJB_SimPostFit.Write();
//  fout->Close();

  if (!_paper) {
 
    //tex->SetTextColor(kWhite); // hide from view
    tex->DrawLatex(0.32,0.64,Form("Npar=%d",_jesFit->GetNpar()));
    //tex->DrawLatex(0.20,0.61,Form("N_{par}=%d (%d)",_jesFit->GetNpar(),
    //				_jesFit->GetNumberFreeParameters()));
    tex->DrawLatex(0.32,0.61,Form("p0=%1.4f #pm %1.4f",
				  _jesFit->GetParameter(0),
				  sqrt(emat[0][0])));
    if (njesFit>=2)
      tex->DrawLatex(0.32,0.58,Form("p1=%1.4f #pm %1.4f",
				    _jesFit->GetParameter(1),
				    sqrt(emat[1][1])));
    if (njesFit==2 && fixTDI)
      tex->DrawLatex(0.32,0.55,Form("p2=%1.3f (TDI fix)",fixTDI));
    if (njesFit>=3)
      tex->DrawLatex(0.32,0.55,Form("p2=%1.3f #pm %1.3f",
				  _jesFit->GetParameter(2),
				    sqrt(emat[2][2])));
    if (njesFit>=4)
      tex->DrawLatex(0.32,0.52,Form("p3=%1.3f #pm %1.3f",
				    _jesFit->GetParameter(3),
				    sqrt(emat[3][3])));
    tex->SetTextSize(0.045); tex->SetTextColor(kBlack);
  }
  
  // Factorize error matrix into eigenvectors
  TVectorD eigvec(np);
  TMatrixD eigmat = emat2.EigenVectors(eigvec);

  for (int ieig = 0; ieig != np; ++ieig) {

    TF1 *fkeig = new TF1(Form("fitEig_feig%d",ieig), jesFit, minpt, maxpt, np);
    fkeig->SetNpx(1000);
    fkeig->SetLineStyle(kDashed);
    for (int i = 0; i != np; ++i) {
      fkeig->SetParameter(i, _jesFit->GetParameter(i)
			  + eigmat[i][ieig] * sqrt(eigvec[ieig]));
    }
    //fkeig->Draw("SAME"); // For visualizing uncertainty correlations vs pT
  } // for ieig

  c0->RedrawAxis();
  c0b->RedrawAxis();
  c1->RedrawAxis();

  if (etamin==0 && (fabs(etamax-1.3)<0.1 || fabs(etamax-2.4)<0.1)) {
    if (dofsr) {
      c0b->SaveAs(Form("pdf/%s/globalFitL3res_raw.pdf",cep));
      //c0b->SaveAs("pdfC/globalFitL3res_raw.C");
      c0->SaveAs(Form("pdf/%s/globalFitL3res_orig.pdf",cep));
      //c0->SaveAs("pdfC/globalFitL3res_orig.C");
      c1->SaveAs(Form("pdf/%s/globalFitL3res_shifted.pdf",cep));
      //c1->SaveAs("pdfC/globalFitL3res_shifted.C");
    }
  }
  else {
    c0b->SaveAs(Form("pdf/%s/globalFitL3res_raw_eta%02.0f-%02.0f.pdf",
		     cep,10*etamin,10*etamax));
    c0->SaveAs(Form("pdf/%s/globalFitL3res_orig_eta%02.0f-%02.0f.pdf",
		    cep,10*etamin,10*etamax));
    c1->SaveAs(Form("pdf/%s/globalFitL3res_shifted_eta%02.0f-%02.0f.pdf",
		    cep,10*etamin,10*etamax));
  }

  // Draw nuisance parameter distribution
  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  hsrc->SetFillColor(kBlue-10);
  hsrc->SetFillStyle(1001);
  hsrc->Draw();
  gStyle->SetOptStat(111111);

  if (etamin==0 && (fabs(etamax-1.3)<0.1 || fabs(etamax-2.4)<0.1))
    c2->SaveAs(Form("pdf/%s/globalFitL3res_hsrc.pdf",cep));
  else {
    c2->SaveAs(Form("pdf/%s/globalFitL3res_hsrc_eta%02.0f-%02.0f.pdf",
		    cep,10*etamin, 10*etamax));
  }
  gStyle->SetOptStat(0);

  // Draw FSR corrections redone

  int colorsdark[4]  = {kBlack,kBlue, kGreen+2,kRed};
  int colorslight[4] = {kGray, kBlue-10, kGreen+2-10,kRed-10};
  for (unsigned int imethod = 0; imethod != nmethods; ++imethod) {

    h->SetMinimum(-0.06);
    h->SetMaximum(+0.09);
    h->SetYTitle("dR / d#alpha (ratio)");
    h = (TH1D*)h->Clone(Form("h3_%d",imethod));

    TCanvas *c3 = tdrCanvas(Form("c3_%d",imethod),h,4,11,true);
    c3->SetLogx();

    TLegend *leg1 = tdrLeg(0.58, _useZoom ? 0.75 : 0.65, 0.88, 0.90);
    leg1->SetHeader("In");
    TLegend *leg2 = tdrLeg(0.65, _useZoom ? 0.75 : 0.65, 0.95, 0.90);
    leg2->SetHeader("Out");

    tex->DrawLatex(0.55,0.20,
		   imethod==0 ? "p_{T} balance method" : "MPF method");

    //    for (int isample = nsample0; isample != nsamples; ++isample) {
    int offsetsample=0;//nsample0; ==>patch to also show dijet fit
    for (int isample = offsetsample; isample != nsamples; ++isample) {

      int ibm = isample + nsamples*imethod;
      TH1D *hk = hks[ibm]; assert(hk);

      //double minx = (isample==nsample0 ? 40  : 30);
      double minx = 30;//(isample==igj ? 40  : 30);
      double maxx = (isample==nsample0 ? 1500 : 1000); // 80X
      hk->GetXaxis()->SetRangeUser(minx,maxx);

      hk->SetLineColor(colorslight[isample-offsetsample]);
      hk->SetFillColor(colorslight[isample-offsetsample]);
      hk->SetFillStyle(kNone);
      hk->SetMarkerSize(0.);
      hk->DrawClone("SAME E5");
      hk->SetFillStyle(3002);
      hk->DrawClone("SAME E5");

      //(new TGraph(hk))->DrawClone("SAMEL");
      TGraph *gk = new TGraph(hk);
      for (int i = gk->GetN()-1; i != -1; --i) {
        if(verboseGF)cout << Form("%s %s pt: %4.0f fsr: %09.6f", samples[isample], methods[imethod], gk->GetX()[i], gk->GetY()[i]) <<endl <<flush;
	//if (gk->GetX()[i] < hk->GetXaxis()->GetXmin())
	if (gk->GetX()[i] < minx || gk->GetX()[i] > maxx)
	  gk->RemovePoint(i);
        if (gk->GetY()[i]==0.)gk->RemovePoint(i);
	//if (gk->GetY()[i]==0 || gk->GetEY()[i]==0) gk->RemovePoint(i);
      }
      gk->Draw("SAMEL");

      leg1->AddEntry(hk, " ", "FL");
    
      // Sum over eigenvectors
      TH1D *hke = (TH1D*)hk->Clone(); hke->Reset();
      for (int ipt = 1; ipt != hke->GetNbinsX()+1; ++ipt) {
	double yeig(0);
	vector<double> df(Npar,0);
	//const int neig = 3;
	for (int ieig = 0; ieig != neig; ++ieig) {

	  int ie = ibm+ieig*nsamples*nmethods;
	  TH1D *hi = hs[ie]; assert(hi);
	
	  yeig += hi->GetBinContent(ipt) * tmp_par[np+ie];
	  df[np+ie] = hi->GetBinContent(ipt);
	} // for ieig

	// FSR is sum of three eigenfunctions so use error propagation for that
	double err2(0);
	for (int i = 0; i != Npar; ++i) {
	  for (int j = 0; j != Npar; ++j) {
	    err2 += emat[i][j]*df[i]*df[j];
	  }
	} // for i

	// Double-check the minus sign for yeig
	// This is probably because kfsr = 1/(1+alpha*hk) so ~ 1 - alpha*hk
	double pt = hke->GetBinCenter(ipt);
	const char *cs = samples[isample];
	double ptreco = (string(cs)=="zeejet"||string(cs)=="zmmjet"
			 ||string(cs)=="zlljet" ?
			 ptreco_zjet : ptreco_gjet);
	double aeff = max(alpha, ptreco/pt);
	hke->SetBinContent(ipt, (aeff*hk->GetBinContent(ipt) - yeig)/aeff);
	hke->SetBinError(ipt, sqrt(err2));
      } // for ipt
      
      //if (isample==nsample0) hke->GetXaxis()->SetRangeUser(40,800);
      //if (isample!=nsample0) hke->GetXaxis()->SetRangeUser(30,300);
      hke->GetXaxis()->SetRangeUser(minx,maxx);
      hke->SetFillStyle(1001);
      hke->SetLineColor(colorsdark[isample-offsetsample]);// hk->GetLineColor()+10);
      hke->DrawClone("SAME E5");
      //(new TGraph(hke))->DrawClone("SAMEL");
      TGraph *gke = new TGraph(hke);
      for (int i = gke->GetN()-1; i != -1; --i) {
        if (gke->GetX()[i] < minx || gke->GetX()[i] > maxx)
          gke->RemovePoint(i);
        if (gke->GetY()[i]==0.)gke->RemovePoint(i);
      }
      gke->Draw("SAMEL");

      leg2->AddEntry(hke, texlabel[samples[isample]], "FL");
    } // for isample

    if (etamin==0 && (fabs(etamax-1.3)<0.1 || fabs(etamax-2.4)<0.1)) {
      c3->SaveAs(Form("pdf/%s/globalFitL3res_%s_kfsr.pdf",
		      cep,methods[imethod]));
      //c3->SaveAs(Form("pdfC/globalFitL3res_%s_kfsr.C", methods[imethod]));
    }
    else {
      c3->SaveAs(Form("pdf/%s/globalFitL3res_%s_kfsr_eta%02.0f-%02.0f.pdf",
		      cep,methods[imethod],10*etamin,10*etamax));
    }
  } // for imethod
  
  //if (etamin==0) c3->SaveAs("pdf/globalFitL3res_kfsr.pdf");

  curdir->cd();
  f->Close();
} // globalFitL3Res


void jesFitter(Int_t& npar, Double_t* grad, Double_t& chi2, Double_t* par,
	       Int_t flag) {

  // Basic checks
  assert(_vdt);
  assert(_vsrc);
  assert(int(_vsrc->size())==npar-_jesFit->GetNpar()-nsrcnmj); // nuisance per source
  //assert(int(_vsrc->size())==npar-_jesFit->GetNumberFreeParameters());
  assert(_vdt->size()<=sizeof(int)*8); // bitmap large enough

//  cout << _lossFunc << endl;
//  assert(_lossFunc!=0);
//  //cout << _lossFunc->GetNDF() << " and " << _lossFunc->GetNumParams() << " and " << (unsigned int)_jesFit->GetNpar() << endl;
//  if(useNewMultijet)assert(_lossFunc->GetNumParams() - nsrcnmj==(unsigned int)_jesFit->GetNpar());
  assert(parTransForMultiJet.size()==(unsigned int)_jesFit->GetNpar());
  parTransForMultiJet[0] = par[0]-1;
  parTransForMultiJet[1] = par[1];
  if(_jesFit->GetNpar()==3)parTransForMultiJet[2] = par[2];

  Double_t *ps = &par[_jesFit->GetNpar()];
  int ns = npar - _jesFit->GetNpar();// fine to include new multijet parameters here... // - nsrcnmj;
  //int ns = npar - _jesFit->GetNumberFreeParameters();

  if (flag) {
    
    // do the following calculation:
    // chi2 = sum_i( (x_i+sum_s(a_s y_si) -fit)^2/sigma_i^2) + sum_s(a_s^2)
    chi2 = 0;
    Nk = 0;

    // Loop over input data (graphs x bins)
    // - for each point, add up source eigenvectors x nuisance parameters
    // - then calculate chi2 adding up residuals + nuisance parameters
    for (unsigned int ig = 0; ig != _vdt->size(); ++ig) {
      TGraphErrors *g = (*_vdt)[ig]; assert(g);
      for (int i = 0; i != g->GetN(); ++i) {

	// Retrieve central value and uncertainty for this point
	double pt = g->GetX()[i];
	double data = g->GetY()[i];
	double sigma = g->GetEY()[i];

	if (pt < ptminJesFit) continue;

	// Calculate fit value at this point
	for (int ipar = 0; ipar != _jesFit->GetNpar(); ++ipar) {
      //for (int ipar = 0; ipar != _jesFit->GetNumberFreeParameters(); ++ipar) {
	  _jesFit->SetParameter(ipar, par[ipar]);
	}
	double fit = _jesFit->EvalPar(&pt,par);

	// For multijet balancing, multiply data by reference JES
	if (TString(g->GetName()).Contains("multijet")) {
	  assert(_vpt);
	  assert(_vpt->size()==_nmethods);
	  unsigned int imethod = ig/_nsamples; assert(imethod<_nmethods);
	  TGraphErrors *gf = (*_vpt)[imethod]; assert(gf);
	  double ptref = gf->GetY()[i] * pt;
	  double fitRef = _jesFit->EvalPar(&ptref,par);
	  data *= fitRef;
	  sigma *= fitRef;
	}

	// Calculate total shift caused by all nuisance parameters
	double shifts = 0;
	for (unsigned int is = 0; is != _vsrc->size(); ++is) {

	  TH1D *hsrc = (*_vsrc)[is]; assert(hsrc);
	  // Bitmap controls which datasets this eigenvector applies to
	  int bitmap; char foo[256];
	  sscanf(hsrc->GetName(),"bm%d_%s", &bitmap,foo);
	  if (!(bitmap>0)) {
	    cout << "Source "<<hsrc->GetName()<<" bitmap "<<bitmap<<endl<<flush;
	  }
	  assert(bitmap>0);
	  bool ison = (bitmap & (1<<ig));
	  //bool ison = (is==ig);//(bitmap & (1<<ig));
	  //bool ison = (is%_vdt->size()==ig);
	  int ipt = hsrc->FindBin(pt);

	  //if (ison) shifts += hsrc->GetBinContent(ipt)
	  //      + ps[is] * hsrc->GetBinError(ipt);
	  if (ison) shifts += ps[is] * hsrc->GetBinContent(ipt);
	}

	// Add chi2 from residual
	double chi = (data + shifts - fit) / sigma;
	chi2 += chi * chi;
	++Nk;

	// if _vdt2 is provided, store shifted data there
	{//if (_vdt2 && _vdt2->size()==_vdt->size()) {
	  assert(_vdt2);
	  assert(_vdt2->size()==_vdt->size());

	  TGraphErrors *g2 = (*_vdt2)[ig];
	  if (g2 && g2->GetN()==g->GetN() && g2->GetX()[i]==pt) {
	    g2->SetPoint(i, pt, data + shifts);

	    // For multijets, store also downward extrapolation
	    if (TString(g->GetName()).Contains("multijet")) {
		assert(_vpt->size()==_nmethods);
		assert(_vpt2->size()==_nmethods);

	      unsigned int imethod = ig/_nsamples; assert(imethod<_nmethods);
	      TGraphErrors *gf = (*_vpt)[imethod]; assert(gf);
	      TGraphErrors *gf2 = (*_vpt2)[imethod]; assert(gf2);
	      double jes = _jesFit->EvalPar(&pt, par);
	      double ptref = pt * gf->GetY()[i];
	      double jesref = _jesFit->EvalPar(&ptref, par);
	      // MJB = jes / jesref
	      // data = MJB*jesref => "jesref" = jes / MJB = jes * jesref/data
	      gf2->SetPoint(i, ptref, jes * jesref / (data + shifts));
	    }
	  }
	}
      } // for ipt
    } // for ig

    // Add chi2 from nuisance parameters //including new multijet
    for (int is = 0; is != ns; ++is) {
      
      chi2 += ps[is]*ps[is];
    } // for ipar

    //cout << Form("standard chi2: %2.4f  multijet chi2: %2.4f \n", chi2, _lossFunc->Eval(parTransForMultiJet));
    //simpy add multijet MPF chi2-term for now and add bining dimensionality from multijet
//    if(useNewMultijet){
//      //      nuisanceDefs_->SetValues(&ps[ns-nsrcnmj]); // set correct index, no copy involved
//      if(nuisances_==0)nuisances_ = new Nuisances(*nuisanceDefs_);
//      nuisances_->SetValues(&ps[ns-nsrcnmj]); // set correct index, no copy involved
//      //      _lossFunc->SetParams(&ps[ns-nsrcnmj]); // set correct index, no copy involved
//      chi2 += _lossFunc->Eval(parTransForMultiJet,*nuisances_);
//      iMultijet++;
//      Nk += _lossFunc->GetNDF() + _lossFunc->GetNumParams();
//      //does this need to be written back?
//      
//    }
    // Give some feedback on progress in case loop gets stuck
    if ((++cnt)%100==0) cout << "." << flush;
  } // if flag
  else {
    if (grad) {}; // suppress warning;
    return;
  }

} // jesFitter


// Calculate fit uncertainty for an arbitrary function
// using standard error propagation with differentials
Double_t fitError(Double_t *xx, Double_t *pp) {

  // Sanity checks on inputs
  assert(_fitError_func);
  assert(_fitError_emat);
  double x = *xx;
  double k = pp[0];
  TF1 *f = _fitError_func;
  int n = f->GetNpar();
  //int n = f->GetNumberFreeParameters();
  TMatrixD &emat = (*_fitError_emat);
  assert(emat.GetNrows()==n);
  assert(emat.GetNcols()==n);
  
  // Partial derivatives as differentials with 10% step size
  vector<double> df(n);
  for (int i = 0; i != n; ++i) {

    double p = f->GetParameter(i);
    double dp = 0.1*sqrt(emat[i][i]); // step size 10% of uncertainty
    f->SetParameter(i, p+dp);
    double fup = f->Eval(x);
    f->SetParameter(i, p-dp);
    double fdw = f->Eval(x);
    f->SetParameter(i, p);

    // Calculate partial derivative as a simple differential
    df[i] = (dp ? (fup - fdw) / (2.*dp) : 0);
  }

  // Perform standard error propagation
  double sumerr2(0);
  for (int i = 0; i != n; ++i) {
    for (int j = 0; j != n; ++j) {
      sumerr2 += emat[i][j]*df[i]*df[j];
    }
  }

  double err = sqrt(sumerr2);

  return (f->Eval(x) + k*err);
}

