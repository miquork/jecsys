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
#include "tdrstyle_mod15.C"

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <functional>
#include <iterator>
#include <utility>
#include <set>

#include <list>
#include <Minuit2/Minuit2Minimizer.h>
#include <Math/Functor.h>

using namespace std;

// plot in DP2021 style instead of Sep 18, 2020 style:
// https://indico.cern.ch/event/955359/contributions/4014344/attachments/2104773/3539752/UL18V3M1_L3Res_2020_09_18.pdf
bool DP2021b = true;

// Global fit settings
double alpha = 0.30; // reference alpha value for zlljet, gamjet (def:0.3)
double ptmj = 30.; // reference pTmin value for multijet (def:30)
bool dofsr = true; // correct for FSR central value (def:true)
bool dofsrcombo1 = true; // use input refIOV for FSR in IOVs (def: true)
bool dofsrcombo2 = false; // use output refIOV for FSR in IOVs (def:false)
bool fixfsrcombo = false; // do not refit FSR (except refIOV) (def:false)
bool useIncJetRefIOVUncert = false; // Include refIOV fit uncertainty? (def:false)
// Settings for 7-parameter fit with PF fractions
bool usePF = true;//true; // Include PF fractions in the fit (def:true for now)
bool usePFZ = true;//false; // use Z+jet PF instead of dijet PF (def:false for now)
bool fitPF = false;
bool fitPFZ = true;
double minPFpt = 74;//114;//49; // Minimum pT where PF fractions used (def:114)
double maxPFpt = 400;//400;//1000;//1500; // Maximum pT where PF fractions used (def:1000)
double minPFZpt = 25;
double maxPFZpt = 600;
int penalizeFitPars = 9; // Add penalty for 9-par fit parameters (def:9)
bool useP0Trk = true; // ok2
bool useP1Gam = true; // ok1
bool useP2CalX = false; // ok1 (292.9/173) (def:false)
bool useP2HadX = true; // ok1 (295.7/173)
bool useP3HadH = true;//false; // over3
bool useP4HadE = true;//false; // nn3
bool useP5Frag = true;//false; // nn3
bool useP6L1RC = true;//false;
bool useP7TrkD = true;
bool useP8MB80pf = false; // (def:false)
bool useP8MB80fit = false; // (def:false)
//
bool useZJet50 = false; // Use Z+jet with alpha<0.50 (def:false)
bool useZJet100 = false; // Use Z+jet with alpha<1.00 (def:false)
double alphaeffmax = 0.45; // max alphaeff for FSR corrections (def:0.45)
double ptreco_gjet = 15.; // min jet pT when evaluating alphamax (def:15)
double ptreco_zlljet = 5.; // same for Zll+jet (def:5, tbu)
double ptreco_zjet = 15.; // same for Z+jet (def:15)
bool dol1bias = false; // correct MPF to L1L2L3-L1 from L1L2L3-RC (def:false)
// see also ScalingForL2OrL3Fit below (left there since long explanation)

// Plotting settings
bool _paper = false; // switch of plotting fit parameters
bool _useZoom = true; // also affects the kind of uncertainty band plotted: useZoom=true comes by default with AbsoluteScale+TotalNoFlavorNoTime; false--> Run1 and reference AbsoluteScale
bool plotIncJet = true; // plot inclusive jet data
bool plotHadW = true; // plot hadronic W data
bool plotMultijetDown = true; // plot gray downward points for multijets
double ptmaxMultijetDown = 300; // max pT for downward multijet points
double shiftPtBal = 0.975; // move x-axis for pTbal, 0 or 1 for none
bool storeErrorMatrix = true; // error matrix for globalFitL3Pulls.C
bool storeJesFit = true; // save jesFit and eigenvectors
bool writeTextFile = true; // textfile for minitools/createL2L3ResTextFile.C
double _cleanUncert = 0.05; // Clean out large uncertainty points from PR plot, for eta>2 mostly

// Scaling options listed below are "None", "putBackL2Res", "DontScaleDijets", "ApplyL3ResDontScaleDijets"
string scalingForL2OrL3Fit = "None";
//"None" - for input combination files without any residual applied (dijets are still scaled, see below)
//"PutBackL2Res" - put L2res back in for gamma/Z+jet for vs eta studies
//"DontScaleDijets" - don't modify inputs at all (good for full closure tests with L2L3Res applied)  -- not needed in case all the reference JEC bands are moved to 1.0 in reprocess.C 
//"ApplyL3ResDontScaleDijets" - apply barrel JES (use case: check closure when only L2Res is applied to the inputs and L3Res didn't change)
//N.B.: Barrel JES from input text file is always applied to dijet results (unless "ApplyL3ResDontScaleDijets" or "DontScaleDijets" is chosen)

// Printing options
bool verboseGF = false;

unsigned int _nsamples(0);
unsigned int _nmethods(0);

int cnt(0); int Nk(0);
TF1 *_jesFit(0);
vector<TGraphErrors*> *_vdt1(0); // original data
vector<TGraphErrors*> *_vmjf(0); // multijet creoil
vector<TGraphErrors*> *_vpt1(0); // original crecoil => change 
vector<TGraphErrors*> *_vdt2(0); // shifted data
vector<TGraphErrors*> *_vpt2(0); // shifted multijet down
vector<TGraphErrors*> *_vdt3(0); // raw data
vector<TGraphErrors*> *_vpt3(0); // shifted multijet up
vector<TH1D*> *_vsrc; // uncertainty sources
vector< pair<string,TGraphErrors*> > _vpf; // active dijet PF fractions
vector< pair<string,TGraphErrors*> > _vpf2; // fitted dijet PF fractions
vector< pair<string,TGraphErrors*> > _vpfz; // active Z+jet PF fractions
vector< pair<string,TGraphErrors*> > _vpfz2; // fitted Z+jet PF fractions
map<string, map<int, TF1*> > _mpf; // map of PF fraction variations

void jesFitter(Int_t &npar, Double_t *grad, Double_t &chi2, Double_t *par,
	       Int_t flag);

// Helper functions to draw fit uncertainty band for arbitrary TF1
TF1 *_fitError_func(0);
TMatrixD *_fitError_emat(0);
Double_t fitError(Double_t *xx, Double_t *p);

// Alternative parameterizations
// 1: EM scale only
// 2: EM+HB (w/ fixOff option)
// 3: EM+HB+offset
// 20200515: V4 is the best of the day, minimum for all IOVs
const int njesFit = 9;//4;//8;//7;//6;//3;//6;//7;//6;//4;//3;//4;//5;
//double fixOff = -0.1239; // E-p3fx 64.1/60
double fixOff = 0;//-0.2220; //CE-p3fx 64.1/60
//double fixOff = -0.1813; // p3 77.4/62
//double fixOff = 0.0190; // p4 76.0/62
  // -0.1655; // p3 75.7/62
  //-0.1464;//-0.6035; // p3 185.9/62
const double ptminJesFit = 15;//30;

//TF1 *fhb(0), *fl1(0), *fl1b(0), *fl1mc(0), *ftr(0), *feg(0); double _etamin(0);
//double _etamin(0);
TF1 *fhb(0), *fl1(0); // Run I functions
//TF1 *fx(0), *fc(0), *fcx(0), *fhh(0), *feh(0), *fch(0), *fp(0), *ft(0), *fd(0);
TF1 *ft(0), *fc(0), *fcx(0); // Deprecating shapes
TF1 *ftd(0), *ftm(0), *fhx(0), *fhh(0), *feh(0), *fp(0), *fhw(0); // Fit shapes
TF1 *fm80(0); // new and experimental
string _epoch ="";
Double_t jesFit(Double_t *x, Double_t *p) {
  
  double pt = *x;

  // Initialize SinglePionHCAL shape (Run I shape)
  if (!fhb) fhb = new TF1("fhb","max(0.,[0]+[1]*pow(x,[2]))",10,3500);
  fhb->SetParameters(1.03091e+00, -5.11540e-02, -1.54227e-01); // SPRH
 
  /*
  // Fits from minitools/varPlots.C
  // -3% to +3% cross variation for SPR calorimeter scale (ECAL+HCAL)
  if (!fx) fx = new TF1("fx","[p0]+[p1]*pow(x/[p2],[p3])/(1+pow(x/[p2],[p3]))*(1-pow(x/[p2],-[p3]))",15,4500);
  fx->SetParameters(1.184, 1.428, 1402, 1.225); // toyPF
  // => different parameters, but essentially same fit as fcx below
  */

  // Regular +3% SPR calorimeter scale variation (ECAL+HCAL)
  if (!fc) fc = new TF1("fc","[p0]+[p1]*pow(x/[p2],[p3])/(1+pow(x/[p2],[p3]))*(1-pow(x/[p2],-[p3]))",15,4500);
  fc->SetParameters(0.335, 0.8171, 406.8, 1.34); // toyPF

  // Fits from minitools/varPlots.C
  // SPR -3% to +3% cross variation (ECAL+HCAL)
  if (!fcx) fcx = new TF1("fcx","[p0]+[p1]*pow(x/[p2],[p3])/(1+pow(x/[p2],[p3]))*(1-pow(x/[p2],-[p3]))",15,4500);
  fcx->SetParameters(1.177, 1.42, 1388, 1.231); // toyPF

  /*
  // SPR +3% variation
  if (!fch) fch = new TF1("fch","[p0]+[p1]*pow(x/[p2],[p3])/(1+pow(x/[p2],[p3]))*(1-pow(x/[p2],-[p3]))",15,4500);
  fch->SetParameters(0.321, 0.7894, 394.5, 1.401); // toyPF
  // => essentially same as fc above
  */

  // SPRH -3% to +3% cross variation (HCAL only)
  if (!fhx) fhx = new TF1("fhx","[p0]+[p1]*pow(x/[p2],[p3])/(1+pow(x/[p2],[p3]))*(1-pow(x/[p2],-[p3]))",15,4500);
  fhx->SetParameters(0.8904, 1.082, 1408, 1.204); // toyPF

  /*
  // SPRH +3% variation
  if (!fhh) fhh = new TF1("fhh","[p0]+[p1]*pow(x/[p2],[p3])/(1+pow(x/[p2],[p3]))*(1-pow(x/[p2],-[p3]))",15,4500);
  fhh->SetParameters(0.7898, 0.5758, 392.2, 1.428); // toyPF
  */
  // SPRH -3% variation
  if (!fhh) fhh = new TF1("fhh","[p0]+[p1]*pow(x/[p2],[p3])/(1+pow(x/[p2],[p3]))*(1-pow(x/[p2],-[p3]))",15,4500);
  fhh->SetParameters(-0.7938, -0.5798, 396.1, 1.412); // toyPF

  // SPRE -3% variation
  if (!feh) feh = new TF1("feh","[p0]+[p1]*pow(x/[p2],[p3])/(1+pow(x/[p2],[p3]))*(1-pow(x/[p2],-[p3]))",15,4500);
  feh->SetParameters(-0.2603, -0.2196, 409.4, 1.276); // toyPF

  // Tracking -3% variation
  if (!ft) ft = new TF1("ft","[p0]+[p1]*pow(x/208.,[p2])",15,4500);
  ft->SetParameters(0.057, -0.3845, -0.3051); // toyPF

  /*
  // Tracking -1% variation in data
  if (!fd) fd = new TF1("fd","[p0]+[p1]*pow(x/208.,[p2])+[p3]/x",15,4500);
  fd->SetParameters(28.16, -28.73, -0.005094, 14.24); // Data E or F
  */

  /*
  // Tracking -1% variation in data
  if (!fd) fd = new TF1("fd","[p0]+[p1]*pow(x/208.,[p2])+[p3]/x",15,4500);
  fd->SetParameters(5.37, -5.811, -0.01197, -0.8343); // Data BCDEF
  */

  // Tracking '-1%' variation in UL2017 data
  if (!ftd) ftd = new TF1("ftd","[p0]+[p1]*pow(x/208.,[p2])+[p3]/x",15,4500);
  ftd->SetParameters(-0.116, -0.6417, -0.3051, 23.63); // Data

  // Tracking '-1%' variation in UL2017 MC
  if (!ftm) ftm = new TF1("ftm","[p0]+[p1]*pow(x/208.,[p2])+[p3]/x",15,4500);
  ftm->SetParameters(0.2683, -0.6994, -0.3051, 18.49); // MC

  // Photon -3% variation
  if (!fp) fp = new TF1("fp","[p0]",15,4500);
  //fp->SetParameter(0,-0.8292); // toyPF (20200511)
  fp->SetParameter(0,-0.8295); // toyPF (20200514)

  // sigmaMB 69.2 mb to 80 mb variation
  if (!fm80) fm80 = new TF1("fm80","[p0]+[p1]*pow(x/208.,[p2])+[p3]*exp(-[p4]*x)",15,4500);
  fm80->SetParameters(0.03904,-0.01612,-0.942, -7.145,0.09907); // toyPF

  // H++ vs P8CP5
  if (!fhw) fhw = new TF1("fhw","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]/x+[p5]*log(x)/x",15,4500);
  fhw->SetParameters(0.9526,-0.3883,1285,2.46,18.1,-2.062); // FullMC

  // Initialize L1FastJet_Simple-L1RC difference
  // Values from fitting ratio/eta00-13/hl1bias (JEC set in reprocess.C)
  if (!fl1) fl1 = new TF1("fl1","1-([0]+[1]*log(x)+[2]*pow(log(x),2))/x",10,3500);
  //fl1->SetParameters(1.72396, 0, 0); // UL17 V1S hl1bias (p1=p2=0) (UL17_V1)
  fl1->SetParameters(0.350077, 0.553560, -0.0527681); // BCDEF hl1rcos (RC vs Simple)

  /*
  if (!fl1b) fl1b = new TF1("fl1b","1-([0]+[1]*log(x)+[2]*pow(log(x),2))/x",10,3500);
  fl1b->SetParameters(-2.62921, 2.34488, -0.429717); // BCDEF hl1cos
  */

  double ptref = 208; // pT that minimizes correlation in p[0] and p[1]
  // (or used to minimize in Run I: for Run II, 208.*sqrt(13/8)=265)

  // Constant scale factor (EM scale) only
  if (njesFit==1) {

    return p[0];
  } // njesFit==2

  // Directly using fitted SPR HCAL shape (from JECUncertainty.cpp)
  if (njesFit==2) {

    // p[0]: overall scale shift, p[1]: HCAL shift in % (full band +3%)
    //return (p[0] + p[1]/3.*100*(fhb->Eval(pt)-fhb->Eval(ptref))); // Run I
    return (p[0]
	    + p[1]*0.01*fcx->Eval(pt)
	    + fixOff*(fl1->Eval(pt)-fl1->Eval(ptref)));
  } // njesFit==2

  if (njesFit==3) {

    // p[0]: overall scale shift, p[1]: HCAL shift in % (full band +3%)
    // p[2]: L1FastJet-L1RC difference
    return (p[0]
	    + p[1]*0.01*fcx->Eval(pt)
	    + p[2]*(fl1->Eval(pt)-fl1->Eval(ptref)));
  } // njesFit==3

  if (njesFit==4) {

    // p[0]: overall scale shift, p[1]: HCAL shift in % (full band +3%)
    // p[2]: L1FastJet-L1RC difference, p[3]: Hadron scale (in HCAL+ECAL)
    return (p[0]
	    //+ p[1]*0.01*fcx->Eval(pt)
	    + p[1]*0.01*fcx->Eval(pt)
	    + p[2]*(fl1->Eval(pt)-fl1->Eval(ptref))
	    + p[3]*0.01*(fc->Eval(pt)-fc->Eval(ptref)));
  } // njesFit==4

  // Full-blown detector and fragmentation model:
  // - Photon, hadrons from toyPF
  // - Tracker inefficiency from data (phi>2.7 in UL2017), compared to toyPF
  // - FullMC fragmentation from H++ vs P8
  // - Offset from L1RC vs L1Simple (update to L1SemiSimple from EpsilonMC)
  if (njesFit==9) {
    assert(!(useP2HadX && useP2CalX));
    return (1
	    + (useP0Trk  ? p[0]*0.01*ftd->Eval(pt) : 0) // Tracker eff. (data)
	    + (useP1Gam  ? p[1]*0.01*fp->Eval(pt) : 0) // photon scale
	    + (useP2CalX ? p[2]*0.01*fcx->Eval(pt) : 0) // Hadron cross-over
	    + (useP2HadX ? p[2]*0.01*fhx->Eval(pt) : 0) // Hadron cross-over
	    + (useP3HadH ? p[3]*0.01*fhh->Eval(pt) : 0) // Hadron in HCAL (SPRH)
	    + (useP4HadE ? p[4]*0.01*feh->Eval(pt) : 0) // Hadron in ECAL (SPRE)
	    //+ (useP5Frag ? p[5]*0.01*fhw->Eval(pt) : 0) // Herwig fragmentation
	    + (useP5Frag ? 0.3*p[5]*0.01*fhw->Eval(pt) : 0) // Herwig fragmentation (rough estimate of impact on Z+jet)
	    + (useP6L1RC ? p[6]*(fl1->Eval(pt)-1) : 0) // L1RC-L1Simple diff.
	    + (useP7TrkD ? p[7]*0.01*3*(ftd->Eval(pt)-ftm->Eval(pt)) : 0) // Tracker Data-MC ('3%' vs '1%' for useP0Trk)
	    + (useP8MB80fit ? p[8]*0.01*fm80->Eval(pt) : 0)
	    );
  } // nJesFit==9

  assert(0);
  exit(0);
}

void globalFitL3Res(double etamin = 0, double etamax = 1.3,
		    string epoch="", string selectSample="Standard_MJDJ_gam_zee_zmm", string selectMethods="PtBalMPF") {
  if(verboseGF)cout << Form("Running globalFitL3Res(etamin=%02.2f,etamax=%02.2f,epoch=%s, selectSample=%s, selectMethods=%s",etamin,etamax,epoch.c_str(),selectSample.c_str(),selectMethods.c_str()) << endl << flush;
  //_etamin = etamin;
  const char *cep = epoch.c_str();
  _epoch = epoch;

  TDirectory *curdir = gDirectory;
  setTDRStyle();

  if (DP2021b) _paper = true;

  bool isUL18 = (epoch=="2018ABCD" || epoch=="2018A" || 
		 epoch=="2018B" || epoch=="2018C" || epoch=="2018D");
  string srefIOV = (isUL18 ? "2018ABCD" : "BCDEF");
  const char *refIOV = srefIOV.c_str();
  bool isRefIOV = (epoch=="2018ABCD" || epoch=="BCDEF");

  TFile *f = new TFile(Form("rootfiles/jecdata%s.root",cep),"READ");
  assert(f && !f->IsZombie());

  TFile *fsys = new TFile(Form("rootfiles/jecdata%s.root",refIOV),"READ");
  assert(fsys && !fsys->IsZombie());

  //TFile *fjes = new TFile(Form("rootfiles/jecdata%s.root",cep),"UPDATE");
  //assert(fjes && !fjes->IsZombie());

  TFile *ffsr = new TFile(Form("rootfiles/jecdata%s.root",refIOV),"READ");
  //isRefIOV ? "UPDATE" : "READ");
  TFile *fout = new TFile(Form("rootfiles/jecdata%s.root",cep),//"UPDATE");
			  storeJesFit || isRefIOV ? "UPDATE" : "READ");
  assert(fout && !fout->IsZombie());

  bool fsrcombo = (dofsrcombo1||dofsrcombo2);
  assert((ffsr && !ffsr->IsZombie()) || (!fsrcombo && !isRefIOV));
  if (dofsrcombo1) cout << Form("Use input FSR from %s for all IOVs\n",refIOV);
  if (dofsrcombo2 && !isRefIOV) cout << Form("Use output FSR from %s\n",refIOV);
  if (fixfsrcombo && !isRefIOV) cout << "Do not refit FSR" << endl;
  
  TFile *fij = (isUL18 ?
		new TFile("rootfiles/drawDeltaJEC_18UL_JECv3.root","READ") :
		new TFile("rootfiles/drawDeltaJEC_17UL_V2M4res_cp2_all_v4.root","READ"));
  assert((fij && !fij->IsZombie()) || ! plotIncJet);
  //TFile *fw = new TFile("rootfiles/hadW.root","READ");
  //assert((fw && !fw->IsZombie()) || ! plotHadW);

  curdir->cd();

  // Settings
  const char *ct = "ratio";

  map<string, std::vector<string> > methodsmap;
  methodsmap["PtBalMPF"] = {"ptchs","mpfchs1"};
  methodsmap["MPF"] = {"mpfchs1"};
  methodsmap["PtBal"] = {"ptchs"};

  assert(methodsmap.find(selectMethods)!=methodsmap.end());
  const unsigned int nmethods = methodsmap[selectMethods].size();
  vector<const char*> methodsvec;
  for(unsigned int i = 0; i < nmethods; ++i)methodsvec.push_back(methodsmap[selectMethods].at(i).c_str());
  const char** methods = &methodsvec.front();
  
  _nmethods = nmethods; // for multijets in global fit

  // Is this wide-bin L3Res fit (as opposed to narrow-bin L2Res)?
  bool isl3 = (etamin==0 && ((epoch!="L4" && fabs(etamax-1.3)<0.1) ||
			     (epoch=="L4" && fabs(etamax-2.4)<0.1)));
  
  map<string, std::vector<string> > samplesmap;
  map<string, int > nsample0map;
  map<string, int > igjmap;
  map<string, int > izeemap;
  map<string, int > izmmmap;
  map<string, int > izllmap;
  map<string, int > iwmhmap;
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

  // Global fit with all samples: multijets/dijets, gamma+jet, merged Z+jet
  samplesmap["MJDJ_gam_z"] =   { (isl3 ? "multijet" : "dijet"),"gamjet", "zjet"};
  nsample0map["MJDJ_gam_z"] = 1;
  igjmap["MJDJ_gam_z"] = 0;
  izeemap["MJDJ_gam_z"] = -1;
  izmmmap["MJDJ_gam_z"] = -1;
  izllmap["MJDJ_gam_z"] = 1;

  // Global fit with all samples and inclusive jets as extra
  samplesmap["MJDJ_inc_gam_zll"] =   { (isl3 ? "multijet" : "dijet"),"incjet","gamjet", "zlljet"};
  nsample0map["MJDJ_inc_gam_zll"] = 2;
  igjmap["MJDJ_inc_gam_zll"] = 0;
  izeemap["MJDJ_inc_gam_zll"] = -1;
  izmmmap["MJDJ_inc_gam_zll"] = -1;
  izllmap["MJDJ_inc_gam_zll"] = 1;

  // Global fit with all samples and inclusive jets + hadw as extra
  samplesmap["MJDJ_inc_gam_zll_hadw"] =   { (isl3 ? "multijet" : "dijet"),"incjet","gamjet", "zlljet", "hadw"};
  nsample0map["MJDJ_inc_gam_zll_hadw"] = 2;
  igjmap["MJDJ_inc_gam_zll_hadw"] = 0;
  izeemap["MJDJ_inc_gam_zll_hadw"] = -1;
  izmmmap["MJDJ_inc_gam_zll_hadw"] = -1;
  izllmap["MJDJ_inc_gam_zll_hadw"] = 1;
  iwmhmap["MJDJ_inc_gam_zll_hadw"] = 2;

  // Global fit with: multijets/dijets, gamma+jet, merged Z+jet, hadronic w
  samplesmap["MJDJ_gam_zll_hadw"] =   { (isl3 ? "multijet" : "dijet"),"gamjet", "zlljet","hadw"};
  nsample0map["MJDJ_gam_zll_hadw"] = 1;
  igjmap["MJDJ_gam_zll_hadw"] = 0;
  izeemap["MJDJ_gam_zll_hadw"] = -1;
  izmmmap["MJDJ_gam_zll_hadw"] = -1;
  izllmap["MJDJ_gam_zll_hadw"] = 1;
  iwmhmap["MJDJ_gam_zll_hadw"] = 2;

  // Global fit with: multijets/dijets, gamma+jet, merged Z+jet, hadronic w
  samplesmap["MJDJ_gam_z_hadw"] =   { (isl3 ? "multijet" : "dijet"),"gamjet", "zjet","hadw"};
  nsample0map["MJDJ_gam_z_hadw"] = 1;
  igjmap["MJDJ_gam_z_hadw"] = 0;
  izeemap["MJDJ_gam_z_hadw"] = -1;
  izmmmap["MJDJ_gam_z_hadw"] = -1;
  izllmap["MJDJ_gam_z_hadw"] = 1;
  iwmhmap["MJDJ_gam_z_hadw"] = 2;

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
  cout << "Selected sample config: " << selectSample << endl;
  const int nsamples = samplesmap[selectSample].size();
  vector<const char*> samplevec;
  set<string> sampleset;
  for(int i = 0; i < nsamples; ++i) {
    samplevec.push_back(samplesmap[selectSample].at(i).c_str());
    sampleset.insert(samplesmap[selectSample].at(i));
  }
  const char** samples = &samplevec.front();
  //const int nsample0 = 1; // first Z/gamma+jet sample
  const int nsample0 = nsample0map[selectSample];
  const int igj = igjmap[selectSample];
  const int izee = izeemap[selectSample];
  const int izmm = izmmmap[selectSample];
  const int izll = izllmap[selectSample];
  const int iwmh = iwmhmap[selectSample];

  _nsamples = nsamples; // for multijets in global fit
  
  string sbin = Form("eta%02.0f-%02.0f",10*etamin,10*etamax);
  const char* bin = sbin.c_str();

  assert(fsys->cd(ct));
  TDirectory *dsys1 = fsys->GetDirectory(ct); assert(dsys1);
  assert(dsys1->cd(bin));
  TDirectory *dsys = dsys1->GetDirectory(bin); assert(dsys);
  curdir->cd();

  //assert(fjes->cd(ct));
  //TDirectory *djes1 = fjes->GetDirectory(ct); assert(djes1);
  //assert(djes1->cd(bin));
  //TDirectory *djes = djes1->GetDirectory(bin); assert(djes);
  //curdir->cd();

  TDirectory *dfsr(0);
  if (fsrcombo || isRefIOV) {
    assert(ffsr->cd(ct));
    TDirectory *dfsrin1 = ffsr->GetDirectory(ct); assert(dfsrin1);
    assert(dfsrin1->cd(bin));
    dfsr = dfsrin1->GetDirectory(bin); assert(dfsr);
    curdir->cd();
  } // fsrcombo

  assert(fout->cd(ct));
  TDirectory *dout1 = fout->GetDirectory(ct); assert(dout1);
  assert(dout1->cd(bin));
  TDirectory *dout = dout1->GetDirectory(bin); assert(dout);
  curdir->cd();

  assert(f->cd(ct));
  TDirectory *din1 = f->GetDirectory(ct); assert(din1);
  assert(din1->cd(bin));
  TDirectory *d = din1->GetDirectory(bin); assert(d);
  //d->cd();
  curdir->cd();

  //////////////////////////////////////////////////
  // Load and setup fractions and their variations
  //////////////////////////////////////////////////

  //vector<pair<string, TGraphErrors> > _vpf;
  //map<string, map<int, TF1*> > _mpf;
  if (usePF || usePFZ) {
    //if (true) {
    //assert(njesFit==6);
    //assert(njesFit==7);
    assert(njesFit==9);

    TGraphErrors *gchf(0), *gnhf(0), *gnef(0);
    TGraphErrors *gchfz(0), *gnhfz(0), *gnefz(0);
    if (usePFZ) {
      gchfz = (TGraphErrors*)d->Get("chf_zlljet_a30");
      gnhfz = (TGraphErrors*)d->Get("nhf_zlljet_a30");
      gnefz = (TGraphErrors*)d->Get("nef_zlljet_a30");
      /*
      gchf = (TGraphErrors*)d->Get("chf_zlljet_a100");
      gnhf = (TGraphErrors*)d->Get("nhf_zlljet_a100");
      gnef = (TGraphErrors*)d->Get("nef_zlljet_a100");
      */
      assert(gchfz);
      assert(gnhfz);
      assert(gnefz);
      // Use open markers for Z+jet if both present and also fitting dijet
      if (usePF && fitPF) {
	gchfz->SetMarkerStyle(kOpenCircle);
	gnhfz->SetMarkerStyle(kOpenDiamond);
	gnefz->SetMarkerStyle(kOpenSquare);
      }
      _vpfz.push_back(make_pair<string,TGraphErrors*>("chf", move(gchfz)));
      _vpfz.push_back(make_pair<string,TGraphErrors*>("nhf", move(gnhfz)));
      _vpfz.push_back(make_pair<string,TGraphErrors*>("nef", move(gnefz)));
      
      // Graphs for output fractions
      curdir->cd();
      TGraphErrors *gchfz2 = (TGraphErrors*)gchfz->Clone("gchfz2");
      TGraphErrors *gnhfz2 = (TGraphErrors*)gnhfz->Clone("gnhfz2");
      TGraphErrors *gnefz2 = (TGraphErrors*)gnefz->Clone("gnefz2");
      _vpfz2.push_back(make_pair<string,TGraphErrors*>("chf", move(gchfz2)));
      _vpfz2.push_back(make_pair<string,TGraphErrors*>("nhf", move(gnhfz2)));
      _vpfz2.push_back(make_pair<string,TGraphErrors*>("nef", move(gnefz2)));
    }
    if (usePF) {
      gchf = (TGraphErrors*)d->Get("chf_pfjet_a30");
      gnhf = (TGraphErrors*)d->Get("nhf_pfjet_a30");
      gnef = (TGraphErrors*)d->Get("nef_pfjet_a30");
      assert(gchf);
      assert(gnhf);
      assert(gnef);

      // Use open markers for dijet if both present, and not fitting Z
      if (fitPFZ && !fitPF) {
	gchf->SetMarkerStyle(kOpenCircle);
	gnhf->SetMarkerStyle(kOpenDiamond);
	gnef->SetMarkerStyle(kOpenSquare);
      }

      //vector<pair<string, TGraphErrors*> > _vpf;
      // std::move needed to turn lvalue to rvalue, as std::make_pair
      // only accepts T2&&, which is rvalue in C++11 standard
      _vpf.push_back(make_pair<string,TGraphErrors*>("chf", move(gchf)));
      _vpf.push_back(make_pair<string,TGraphErrors*>("nhf", move(gnhf)));
      _vpf.push_back(make_pair<string,TGraphErrors*>("nef", move(gnef)));
      
      // Graphs for output fractions
      curdir->cd();
      TGraphErrors *gchf2 = (TGraphErrors*)gchf->Clone("gchf2");
      TGraphErrors *gnhf2 = (TGraphErrors*)gnhf->Clone("gnhf2");
      TGraphErrors *gnef2 = (TGraphErrors*)gnef->Clone("gnef2");
      _vpf2.push_back(make_pair<string,TGraphErrors*>("chf", move(gchf2)));
      _vpf2.push_back(make_pair<string,TGraphErrors*>("nhf", move(gnhf2)));
      _vpf2.push_back(make_pair<string,TGraphErrors*>("nef", move(gnef2)));
    }
    
    // Fits from minitools/varPlots.C
    TF1 *ft3_nhf = new TF1("ft3_nhf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,-0.3)",15,4500);
    ft3_nhf->SetParameters(1.312, -0.2902, 951.2, 1.289, -0.6573); // 50.9/18
    TF1 *fp_nhf = new TF1("fp_nhf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500);
    //fp_nhf->SetParameters(0.08365, 0, 1000, 252.7, 0.0001695, 0.9006); // 9.2/18 (was ECAL, not photon; fix 20200610)
    fp_nhf->SetParameters(0.07395, 0, 1000, 258.2, 1.223e-05, 1.158); // 7.6/18
    TF1 *fcx_nhf = new TF1("fcx_nhf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500);
    fcx_nhf->SetParameters(-1.54, -0.4506, 156.3, 0.5089, 1.234, 0.09256); // 3.6/17
    TF1 *fhx_nhf = new TF1("fhx_nhf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500);
    fhx_nhf->SetParameters(-0.295, 0.09444, 2713, 2.292, 0.06437, 0.2845); // 3.6/17
    //TF1 *fh_nhf = new TF1("fh_nhf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500);
    //fh_nhf->SetParameters(-0.4011, -3.571, 7980, 0.3416, 0.7455, 0.2044); // 7.1/17
    TF1 *fhm_nhf = new TF1("fhm_nhf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500);
    fhm_nhf->SetParameters(-0.2746, -0.6358, 9664, 0.6547, 0.05559, 0.1816); // 10.1/17
    //TF1 *fe_nhf = new TF1("fe_nhf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500);
    //fe_nhf->SetParameters(-2.43, 0, 1705, 275.2, 2.375, 0.006368); // 6.7/18
    TF1 *fem_nhf = new TF1("fem_nhf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500);
    fem_nhf->SetParameters(-0.03458, 0, 1713, 274.8, 0.01665, 0.2426); // 3.4/18
    //TF1 *fd_nhf = new TF1("fd_nhf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,+0.3)",15,4500);
    //fd_nhf->SetParameters(0.4208, 6.059, 7174, 0.4467, -0.4718); // 84.4/50
    //TF1 *fd_nhf = new TF1("fd_nhf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,+0.3)",15,4500);
    //fd_nhf->SetParameters(0.1631, 1.035, 1427, 0.5606, -0.1133); // 77.9/52
    //TF1 *ftmg_nhf = new TF1("ftmg_nhf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,+0.3)+[p5]/x",15,4500);
    //ftmg_nhf->SetParameters(0.9795, 20.45, 5.299e+04, 0.408, -0.9336, -3.579); // 256.2/115
    TF1 *ftmg_nhf = new TF1("ftmg_nhf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,+0.3)+[p5]/x",15,4500);
    ftmg_nhf->SetParameters(-0.01022, -0.1962, 4000, 3.071, 0.04211, 0.01005); // 253.1/116
    TF1 *fm80_chf = new TF1("fm80_chf","[p0]+[p1]*pow(x,[p2])",15,4500);
    fm80_chf->SetParameters(1.848, -1.685, -0.006643); // 70.3/53
    TF1 *fhw_nhf = new TF1("fhw_nhf","[p0]+[p1]*pow(x,[p2])+[p3]/x",15,4500);
    fhw_nhf->SetParameters(-5.151, 4.495, 0.03335, -12.3); // 65.4/42

    //<TCanvas::Print>: pdf file pdf/varPlotsComp_nhf.pdf has been created

    // Fits from minitools/varPlots.C
    TF1 *ft3_nef = new TF1("ft3_nef","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,-0.3)",15,4500);
    ft3_nef->SetParameters(-0.8857, 1.481, 2223, -0.07145, -1.411); // 35.2/18
    TF1 *fp_nef = new TF1("fp_nef","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500);
    //fp_nef->SetParameters(0.8609, 0, 1000, 1.302, -1.325, 0.01247); // 19.2/18 (was ECAL, not photon; fix 20200610)
    fp_nef->SetParameters(2.283, 0, 1000, 1.302, -2.738, 0.002452); // 21.3/18
    TF1 *fcx_nef = new TF1("fcx_nef","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500);
    fcx_nef->SetParameters(0.01659, -0.005997, 754.8, 121.2, -0.0001236, 0.9665); // 7.9/17
    TF1 *fhx_nef = new TF1("fhx_nef","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500);
    fhx_nef->SetParameters(0.05474, -0.003141, 798.6, 78.84, -0.000957, 0.7676); // 3.7/17
    //TF1 *fh_nef = new TF1("fh_nef","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500);
    //fh_nef->SetParameters(-0.1723, 1.334, 9774, 0.1691, -0.3062, 0.1916); // 7.7/17
    TF1 *fhm_nef = new TF1("fhm_nef","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500);
    fhm_nef->SetParameters(0.4158, -2.14, 9426, 0.1723, 0.4111, 0.1937); // 9.5/17
    //TF1 *fe_nef = new TF1("fe_nef","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500);
    //fe_nef->SetParameters(-0.02182, 0, 1573, 274.7, -0.01082, 0.2467); // 5.3/18
    TF1 *fem_nef = new TF1("fem_nef","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500);
    fem_nef->SetParameters(-0.02364, 0, 1481, 246.2, -0.009737, 0.2576); // 5.3/18
    //TF1 *fd_nef = new TF1("fd_nef","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,+0.3)",15,4500);
    //fd_nef->SetParameters(-12.22, 5.872, 1.487e+04, -0.5201, 0.3729); // 80.7/50
    //TF1 *fd_nef = new TF1("fd_nef","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,+0.3)",15,4500);
    //fd_nef->SetParameters(-8.455, 4.092, 2.45e+04, -0.6258, 0.1758); // 61.6/52
    //TF1 *ftmg_nef = new TF1("ftmg_nef","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,+0.3)+[p5]/x",15,4500);
    //ftmg_nef->SetParameters(-6.473, 3.525, 3.595e+05, -6.328, -0.04752, -8.047); // 300.0/115
    //TF1 *ftmg_nef = new TF1("ftmg_nef","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,+0.3)+[p5]/x",15,4500);
    //ftmg_nef->SetParameters(-7.798, 8.892, 472.9, 0.04415, -0.1297, -1.122); // 236.2/115 (V2)
    //TF1 *ftmg_nef = new TF1("ftmg_nef","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,+0.3)+[p5]/x",15,4500);
    //ftmg_nef->SetParameters(0.2257, -0.168, 1180, 3.338, 0.01264, -3.072); // 151.4/115 (V3)
    TF1 *ftmg_nef = new TF1("ftmg_nef","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])+[p6]/x",15,4500);
    ftmg_nef->SetParameters(0.07453, 0.1457, 1131, -3.68, -0.4155, -0.3, -1.878); // 149.7/115
    TF1 *fm80_nef = new TF1("fm80_nef","[p0]+[p1]*pow(x,[p2])",15,4500);
    fm80_nef->SetParameters(3.611, -3.909, -0.007482); // 70.2/53
    TF1 *fhw_nef = new TF1("fhw_nef","[p0]+[p1]*pow(x,[p2])+[p3]/x",15,4500);
    fhw_nef->SetParameters(0.8417, -0.2605, 0.2289, 2.426); // 35.9/42

    //<TCanvas::Print>: pdf file pdf/varPlotsComp_nef.pdf has been created

    // Fits from minitools/varPlots.C
    TF1 *ft3_chf = new TF1("ft3_chf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,-0.3)",15,4500);
    ft3_chf->SetParameters(-1.792, 0.2976, 1107, 1.559, 1.039); // 64.1/18
    TF1 *fp_chf = new TF1("fp_chf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500);
    //fp_chf->SetParameters(-0.8522, 1.708, 904.6, 0.09591, -0.05902, 0.2962); // 13.2/17 (was ECAL, not photon; fix 20200610)
    fp_chf->SetParameters(0.3333, 0.7433, 1023, 0.3926, -0.09446, 0.2883); // 10.3/17
    TF1 *fcx_chf = new TF1("fcx_chf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500);
    fcx_chf->SetParameters(-0.1188, -0.3705, 408.1, -0.2583, 1.39, -0.1831); // 2.7/17
    TF1 *fhx_chf = new TF1("fhx_chf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500);
    fhx_chf->SetParameters(-0.0637, -0.2811, 4531, -0.3172, 1.071, -0.153); // 1.7/17
    //TF1 *fh_chf = new TF1("fh_chf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500);
    //fh_chf->SetParameters(-0.131, 0.03603, 332.2, 3.025, 0.005807, -3.233); // 4.4/17
    TF1 *fhm_chf = new TF1("fhm_chf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500);
    fhm_chf->SetParameters(0.1552, -0.04221, 315.4, 2.787, -0.06628, -0.2572); // 6.2/17
    //TF1 *fe_chf = new TF1("fe_chf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500);
    //fe_chf->SetParameters(0.06035, 0, 1000, 1.3, -0.007813, 0.2174); // 2.7/18
    TF1 *fem_chf = new TF1("fem_chf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500);
    fem_chf->SetParameters(0.06085, 0, 1000, 1.3, -0.008137, 0.2135); // 2.7/18
    //TF1 *fd_chf = new TF1("fd_chf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,+0.3)",15,4500);
    //fd_chf->SetParameters(0.0799, -1.954, 1142, 0.393, 0.1761); // 98.9/50
    //TF1 *fd_chf = new TF1("fd_chf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,+0.3)",15,4500);
    //fd_chf->SetParameters(-0.2958, -5.025, 2400, 0.4452, 0.5037); // 68.3/52
    TF1 *ftmg_chf = new TF1("ftmg_chf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,+0.3)+[p5]/x",15,4500);
    ftmg_chf->SetParameters(1.982, -2.678, 47.02, 0.262, 0.1494, -3.097); // 214.9/115
    TF1 *fm80_nhf = new TF1("fm80_nhf","[p0]+[p1]*pow(x,[p2])",15,4500);
    fm80_nhf->SetParameters(-0.05047, -0.0008452, 0.6402); // 56.7/53
    TF1 *fhw_chf = new TF1("fhw_chf","[p0]+[p1]*pow(x,[p2])+[p3]/x",15,4500);
    fhw_chf->SetParameters(-0.2176, 1.064e-05, 1.373, 0); // 41.8/43


    //<TCanvas::Print>: pdf file pdf/varPlotsComp_chf.pdf has been created
    //<TCanvas::Print>: pdf file pdf/varPlotsComp_all.pdf has been created

    //map<string, map<int, TF1*> > _mpf;
    //_mpf["chf"][0] = ft3_chf; // tracking in MC
    _mpf["chf"][0] = (useP0Trk  ? ftmg_chf : 0); // tracking in data+MC
    _mpf["chf"][1] = (useP1Gam  ? fp_chf : 0); // photon
    _mpf["chf"][2] = (useP2CalX ? fcx_chf : 0); // SPR cross (HCAL+ECAL)
    _mpf["chf"][2] = (useP2HadX ? fhx_chf : 0); // SPR cross (HCAL)
    _mpf["chf"][3] = (useP3HadH ? fhm_chf : 0); // SPR-hcal
    _mpf["chf"][4] = (useP4HadE ? fem_chf : 0); // SPR-ecal
    _mpf["chf"][5] = (useP5Frag ? fhw_chf : 0); // fragmentation
    _mpf["chf"][6] = 0;//(useP6L1RC ? 0 : 0);
    _mpf["chf"][7] = 0;//(useP7TrkD ? 0 : 0);
    _mpf["chf"][8] = (useP8MB80pf ? fm80_chf : 0); // fragmentation
    //
    //_mpf["nhf"][0] = ft3_nhf; // tracking in MC
    _mpf["nhf"][0] = (useP0Trk  ? ftmg_nhf : 0); // tracking in data+MC
    _mpf["nhf"][1] = (useP1Gam  ? fp_nhf : 0); // photon
    _mpf["nhf"][2] = (useP2CalX ? fcx_nhf : 0); // SPR cross (HCAL+ECAL)
    _mpf["nhf"][2] = (useP2HadX ? fhx_nhf : 0); // SPR cross (HCAL)
    _mpf["nhf"][3] = (useP3HadH ? fhm_nhf : 0); // SPR-hcal
    _mpf["nhf"][4] = (useP4HadE ? fem_nhf : 0); // SPR-ecal
    _mpf["nhf"][5] = (useP5Frag ? fhw_nhf : 0);
    _mpf["nhf"][6] = 0;//(useP6L1RC ? 0 : 0);
    _mpf["nhf"][7] = 0;//(useP7TrkD ? 0 : 0);
    _mpf["nhf"][8] = (useP8MB80pf ? fm80_nhf : 0);
    //
    //_mpf["nef"][0] = ft3_nef; // tracking in MC
    _mpf["nef"][0] = (useP0Trk ? ftmg_nef : 0); // tracking in data+MC
    _mpf["nef"][1] = (useP1Gam ? fp_nef : 0); // photon
    _mpf["nef"][2] = (useP2CalX ? fcx_nef : 0); // SPR cross (HCAL+ECAL)
    _mpf["nef"][2] = (useP2HadX ? fhx_nef : 0); // SPR cross (HCAL)
    _mpf["nef"][3] = (useP3HadH ? fhm_nef : 0); // SPR-hcal
    _mpf["nef"][4] = (useP4HadE ? fem_nef : 0); // SPR-ecal
    _mpf["nef"][5] = (useP5Frag ? fhw_nef : 0);
    _mpf["nef"][6] = 0;//(useP6L1RC ? 0 : 0);
    _mpf["nef"][7] = 0;//(useP7TrkD ? 0 : 0);
    _mpf["nef"][8] = (useP8MB80pf ? fm80_nef : 0);
  } // usePF


  ///////////////////////////////////////////////
  // Load and setup data and uncertainty sources
  ///////////////////////////////////////////////

  // Load data points
  vector<TGraphErrors*> gs1; // original
  vector<TGraphErrors*> gs2; // shifted
  vector<TGraphErrors*> gs3; // raw
  for (unsigned int imethod = 0; imethod != nmethods; ++imethod) {
    for (int isample = 0; isample != nsamples; ++isample) {

      string sm = methods[imethod];
      const char *cm = sm.c_str();
      string ss = samples[isample];
      const char *cs = ss.c_str();

      string s = Form("%s_%s_a%02.0f",cm,cs,
		      ss=="multijet" ? ptmj : alpha*100);
      if (useZJet50  && (ss=="zjet"||ss=="zeejet"||ss=="zmmjet"||ss=="zlljet"))
	s = Form("%s_%s_a%02.0f",cm,cs,50.);
      if (useZJet100 && (ss=="zjet"||ss=="zeejet"||ss=="zmmjet"||ss=="zlljet"))
	s = Form("%s_%s_a%02.0f",cm,cs,100.);
      TGraphErrors *g = (TGraphErrors*)d->Get(s.c_str());
      if (!g) cout << "Graph "<<s<<" not found!" << endl << flush;
      assert(g);

      // Clone to keep original data in file intact
      g = (TGraphErrors*)g->Clone(Form("g_%s",g->GetName()));
      
      //if (ss=="multijet") {
      //g->SetMarkerStyle(imethod==0 ? kOpenTriangleUp : kFullTriangleUp);
      //}
      //g->SetLineWidth(2);

      // Clean out empty points
      // as well as trigger-biased ones for dijets
      for (int i = g->GetN()-1; i != -1; --i) {
	if (g->GetY()[i]==0 || g->GetEY()[i]==0 ||
	    (ss=="dijet" && g->GetX()[i]<70.))
	  g->RemovePoint(i);
      }
      
      gs1.push_back((TGraphErrors*)g->Clone(Form("%s_original",g->GetName())));
      gs2.push_back((TGraphErrors*)g->Clone(Form("%s_shifted",g->GetName())));
      gs3.push_back((TGraphErrors*)g->Clone(Form("%s_raw",g->GetName())));
    } // for isample
  } // for imethod

  _vdt1 = &gs1;
  _vdt2 = &gs2;
  _vdt3 = &gs3;
  
  // Load pT fractions for multijet balancing
  vector<TGraphErrors*> gfs, gfs1, gfs2, gfs3;
  //if (string(samples[0])=="multijet") {
  if (sampleset.find("multijet")!=sampleset.end()) {

    // Fractions for MJB (pT>30 GeV)
    string s = Form("crecoil_multijet_a%1.0f",ptmj);
    const char *c = s.c_str();
    TGraphErrors *g = (TGraphErrors*)d->Get(s.c_str());
    if (!g) cout << "Graph "<<s<<" not found!" << endl << flush;
    assert(g);
    
    g->SetLineWidth(2);
    g->SetMarkerStyle(kOpenTriangleDown);
    gfs .push_back((TGraphErrors*)g->Clone(Form("%s_mjb_crecoil",c)));
    gfs1.push_back((TGraphErrors*)g->Clone(Form("%s_mjb_original",c)));
    gfs2.push_back((TGraphErrors*)g->Clone(Form("%s_mjb_shifted",c)));
    gfs3.push_back((TGraphErrors*)g->Clone(Form("%s_mjb_raw",c)));
    
    // Fractions for MPF (pT>15 GeV)
    s = "crecoil_multijet_a15";
    //if (isUL18) s = "crecoil_multijet_a30"; // PATCH
    g = (TGraphErrors*)d->Get(s.c_str());
    if (!g) cout << "Graph "<<s<<" not found!" << endl << flush;
    assert(g);

    g->SetLineWidth(2);
    g->SetMarkerStyle(kFullTriangleDown);
    gfs .push_back((TGraphErrors*)g->Clone(Form("%s_mpf_crecoil",c)));
    gfs1.push_back((TGraphErrors*)g->Clone(Form("%s_mpf_original",c)));
    gfs2.push_back((TGraphErrors*)g->Clone(Form("%s_mpf_shifted",c)));
    gfs3.push_back((TGraphErrors*)g->Clone(Form("%s_mpf_raw",c)));
  }

  _vmjf = &gfs;
  _vpt1 = &gfs1;
  _vpt2 = &gfs2;
  _vpt3 = &gfs3;

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

  // Load inclusive jet reference JES
  TH1D *hjesfitref = (TH1D*)dsys->Get(isRefIOV?"herr_ref":"sys/hjesfit");
  assert(hjesfitref);
  
  // Load inclusive jet data
  TH1D *hij(0);
  TGraphErrors *gij(0);
  if (plotIncJet) {
    map<string,const char*> fij_files;
    fij_files["2018A"] = "A";
    fij_files["2018B"] = "B";
    fij_files["2018C"] = "C";
    fij_files["2018D"] = "D";
    fij_files["2018ABCD"] = "";
    hij = (isUL18 ?
	   (TH1D*)fij->Get(Form("jet_Run18UL%s_det",fij_files[epoch])) :
	   (TH1D*)fij->Get(Form("jet_Run17UL%s_fwd3",
				isRefIOV ? "" : epoch.c_str())));
    assert(hij);
    for (int i = 1; i != hij->GetNbinsX()+1; ++i) {
      double pt = hij->GetBinCenter(i);
      double jes = hjesfitref->GetBinContent(hjesfitref->FindBin(pt));
      hij->SetBinContent(i, hij->GetBinContent(i) * jes);
    } // for i

    gij = (TGraphErrors*)d->Get("mpfchs1_incjet_a30");
    assert(gij);
    for (int i = 0; i != gij->GetN(); ++i) {
      double pt = gij->GetX()[i];
      double jes = hjesfitref->GetBinContent(hjesfitref->FindBin(pt));
      gij->SetPoint(i, gij->GetX()[i], gij->GetY()[i] * jes);
    } // for i
  }

  // Load hadronic W data (if being used)
  TH1D *hw(0);
  if (sampleset.find("hadw")!=sampleset.end()) {// && plotHadW) {
    //assert(fw);
    //hw = (TH1D*)fw->Get("data_nevents_hadw_fitprob02_L1L2L3");
    hw = (TH1D*)d->Get("counts_hadw_a30");
    assert(hw);
  }
  
  // Load FSR (JESref) corrections, and correct data
  vector<TH1D*> hks(nsamples*nmethods,0);
  for (unsigned int imethod = 0; imethod != nmethods; ++imethod) {
    for (int isample = 0; isample != nsamples; ++isample) {

      string sm = methods[imethod];
      const char *cm = sm.c_str();
      string ss = samples[isample];
      const char *cs = ss.c_str();

      int ibm = isample + nsamples*imethod;
      TH1D *hfsr(0);

      // No FSR corrections for inclusive jets
      if (ss=="incjet") {
	hfsr = (TH1D*)hij->Clone(Form("bm%d_%s",(1<<ibm),hij->GetName()));
	hfsr->Reset();
	hks[ibm] = hfsr;
      }
      // fitProb-based FSR uncertaity for hadronic W, but first zero
      else if (ss=="hadw") {
	assert(hw);
	hfsr = (TH1D*)hw->Clone(Form("bm%d_%s",(1<<ibm),hw->GetName()));
	hfsr->Reset();
	hks[ibm] = hfsr;
      }
      else {
	// For correcting FSR bias
	string s = Form("fsr/hkfsr_%s_%s",cm,cs);
	//if (ss=="zjet") s = Form("fsr/hkfsr_%s_zlljet",cm);
	hfsr = (TH1D*)d->Get(s.c_str());
	if (dofsrcombo2 && !isRefIOV) s = Form("fsr/hkfsr2_%s_%s",cm,cs);
	if (fsrcombo && !isRefIOV) hfsr = (TH1D*)dfsr->Get(s.c_str());
	if (!hfsr) cout << "Histo "<<s<<" not found!" << endl << flush;
      }
      assert(hfsr);

      // For correcting L1L2L3-L1 to L1L2L3-RC
      TH1D *hl1 = (TH1D*)d->Get("hl1bias"); assert(hl1);
      
      hfsr->SetName(Form("bm%d_%s",(1<<ibm),hfsr->GetName()));

      hks[ibm] = hfsr;

      // Correct data for FSR
      TGraphErrors *g1 = (*_vdt1)[ibm];
      TGraphErrors *g2 = (*_vdt2)[ibm];
      TGraphErrors *g3 = (*_vdt3)[ibm];
      TGraphErrors *gf1 = (_vpt1->size()==nmethods ? (*_vpt1)[imethod] : 0);
      TGraphErrors *gf2 = (_vpt2->size()==nmethods ? (*_vpt2)[imethod] : 0);
      TGraphErrors *gf3 = (_vpt3->size()==nmethods ? (*_vpt3)[imethod] : 0);
      for (int i = 0; i != g1->GetN(); ++i) {
	
	double pt = g1->GetX()[i];
	double r = g1->GetY()[i];
	double ptreco(0);
	if (ss=="zeejet"||ss=="zmmjet"||ss=="zlljet") ptreco = ptreco_zlljet;
	if (ss=="zjet") ptreco = ptreco_zjet;
	if (ss=="gamjet") ptreco = ptreco_gjet;
	double aeff = min(max(alpha, ptreco/pt),alphaeffmax);
	if (ss=="multijet") aeff = alpha;
	double kfsr(1);
	if (dofsr) kfsr = 1./(1+aeff*hfsr->GetBinContent(hfsr->FindBin(pt)));
	double l1(1);
	if (dol1bias && sm=="mpfchs1" && ss=="gamjet")
	  l1 = 1./hl1->GetBinContent(hl1->FindBin(pt));

	double scale = 1.00; // correct out previous L3Res

	if (ss=="incjet") {
	  //scale = herr_ref->GetBinContent(herr_ref->FindBin(pt));
	  scale = hjesfitref->GetBinContent(hjesfitref->FindBin(pt));
	}

	// put L2res back in for gamma/Z+jet for vs eta studies
	if (!(isl3 && ss!="dijet" && ss!="multijet")) {

	  double jes = hjes->GetBinContent(hjes->FindBin(pt)); // L2Res only
	  double jesref = hjes0->GetBinContent(hjes0->FindBin(pt)); // barrel
	  // divide by jesref(=1) in case hjes had L2L3Res instead of L2Res
          if (scalingForL2OrL3Fit=="None" ||
	      scalingForL2OrL3Fit=="DontScaleDijets") scale = 1.0;
          else if (scalingForL2OrL3Fit=="PutBackL2Res") scale = jes / jesref;
          else if (scalingForL2OrL3Fit=="ApplyL3ResDontScaleDijets") scale = 1. / jesref;
          else assert(false);
          cout << "scale..." << scale << endl;
	}
	g1->SetPoint(i, pt, scale*l1*r*kfsr);
	g2->SetPoint(i, pt, scale*l1*r*kfsr);
	g3->SetPoint(i, pt, scale*l1*r);
	
	// For multijet, correct instead for reference JES (and FSR)
	if (ss=="multijet") {
	  assert(gf1);
	  double ptref = gf1->GetY()[i] * pt;
	  double jesref = herr->GetBinContent(herr->FindBin(ptref));
	  // Don't correct g1, we need raw input for global fit
	  //g1->SetPoint(i, pt, jesref * r);
	  g1->SetPoint(i, pt, r * kfsr);
	  g2->SetPoint(i, pt, jesref * r * kfsr);
	  g3->SetPoint(i, pt, jesref * r);
	  double jes = herr->GetBinContent(herr->FindBin(pt));
	  //gf2->SetPoint(i, ptref, jes / r);
	  gf2->SetPoint(i, ptref, jes / (r * kfsr));
	  gf2->SetPointError(i, g1->GetEX()[i], g1->GetEY()[i]);
	  gf3->SetPoint(i, ptref, jes / r);
	  gf3->SetPointError(i, g1->GetEX()[i], g1->GetEY()[i]);
	}

	// For dijet, multiply by barrel JES
	if (ss=="dijet" &&
	    !(scalingForL2OrL3Fit=="ApplyL3ResDontScaleDijets" ||
	      scalingForL2OrL3Fit=="DontScaleDijets")) {
	  double jesref = herr0->GetBinContent(herr0->FindBin(pt));
	  g1->SetPoint(i, pt, r*kfsr * jesref);
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
  const int neigmax = 2;//4;//3;
  for (int ieig = 0; ieig != neigmax; ++ieig) {
    for (unsigned int imethod = 0; imethod != nmethods; ++imethod) {
      for (int isample = 0; isample != nsamples; ++isample) {
	
	string sm = methods[imethod];
	const char *cm = sm.c_str();
	string ss = samples[isample];
	const char *cs = ss.c_str();

	// bitmap for which methods the sources applies to
	int ibm = isample + nsamples*imethod;	
	
	string s = Form("fsr/hkfsr_%s_%s_eig%d",cm,cs,ieig);
	//if (ss=="zjet") s = Form("fsr/hkfsr_%s_zlljet_eig%d",cm,ieig);
	TH1D *h = (TH1D*)d->Get(s.c_str());
	if (dofsrcombo2 && !isRefIOV)
	  s = Form("fsr/hkfsr2_%s_%s_eig%d",cm,cs,ieig);
	if (fsrcombo && !isRefIOV) h = (TH1D*)dfsr->Get(s.c_str());
	if (ieig>0 && !h) {
	  h = hs[ibm]; assert(h); // use src0
	  h = (TH1D*)h->Clone(Form("bm%d_inactive_hkfsr_%s_%s_eig%d",
	  			   (1<<ibm),cm,cs,ieig));
	  h->Reset();
	  // Store also empty ones for now to have FSR plots work
	  hs.push_back(h);
	  continue;
	}
	if (ss=="incjet" && !h && hij) {
	  h = (TH1D*)hij->Clone(Form("bm%d_inactive_hkfsr_%s_%s_eig%d",
	  			   (1<<ibm),cm,cs,ieig));
	  h->Reset();
	  // Store also empty ones for now to have FSR plots work
	  hs.push_back(h);
	  continue;
	}
	if (ss=="hadw" && !h && hw) {
	  h = (TH1D*)hw->Clone(Form("bm%d_inactive_hkfsr_%s_%s_eig%d",
				    (1<<ibm),cm,cs,ieig));
	  h->Reset();
	  // Store also empty ones for now to have FSR plots work
	  hs.push_back(h);
	  continue;
	}
	if (!h) cout << "Histo "<<s<<" not found!" << endl << flush;
	assert(h);

	// Clone to not change the original in file
	h = (TH1D*)h->Clone(Form("bm%d_%s",(1<<ibm),h->GetName()));

	// Forbid refitting if fixfsrcombo and not BCDEF
	if (fixfsrcombo && !isRefIOV)
	  h->Reset();
	
	// Scale dR/dalpha to dR for given alpha setting
	for (int i = 1; i != h->GetNbinsX()+1; ++i) {
	  double pt = h->GetBinCenter(i);
	  double ptreco(0);
	  if (ss=="zeejet"||ss=="zmmjet"||ss=="zlljet") ptreco = ptreco_zlljet;
	  if (ss=="zjet") ptreco = ptreco_zjet;
	  if (ss=="gamjet") ptreco = ptreco_gjet;
	  double aeff = min(max(alpha, ptreco/pt),alphaeffmax);
	  if (ss=="multijet") aeff = alpha;
	  h->SetBinContent(i, aeff*h->GetBinContent(i));
	}
	//h->SetName(Form("bm%d_%s",(1<<ibm),h->GetName()));

	hs.push_back(h);
      } // for isample
    } // for imethod
  } // for ieig

  // Uncertainty sources for e, gamma and mu scale uncertainties
  // (and extra source for multijets?)
  int is(0);
  int is_gj(0), is_zee(0), is_zmm(0), is_zll(0), is_z(0);
  for (unsigned int i = 0; i != _vdt1->size(); ++i) {
    
    unsigned int isample = i%nsamples;
    string ss = samples[isample];
    const char *cs = ss.c_str();
    unsigned int imethod = i/nsamples;
    string sm = methods[imethod];
    const char *cm = sm.c_str();

    TH1D *h = hs[i]; assert(h);
    TH1D *h2 = (TH1D*)h->Clone(Form("bm%d_inactive_%s_%s_%d",1<<i,cm,cs,i));

    double escale(0);
    const int n0 = nsample0;
    const int n1 = nsamples+nsample0;
    if (ss=="dijet" && sm=="ptchs") {
      escale = 0.005; // for JER bias and whatnot
      // Use same source for both MPF and pT balance
      h2->SetName(Form("bm%d_scale_%03.0f_dijet_%d",
		       (1<<(n0-1) | (1<<(n1-1))), escale*10000., i));
    }
    if (ss=="gamjet" && sm=="ptchs") {
      escale = 0.005; // UL17-v2 (after 1.0% rescale)
      // Use same source for both MPF and pT balance
      //h2->SetName(Form("bm%d_scale_%03.0f_gamjet_%d",
      //(1<<(n0+igj) | (1<<(n1+igj))), escale*10000., i));
      h2->SetName(Form("bm%d_scale_%03.0f_gamjet",
		       (1<<(n0+igj) | (1<<(n1+igj))), escale*10000.));
      is = hs.size();
      is_gj = hs.size();
    }
    if (ss=="zeejet" && sm=="ptchs") {
      escale = 0.002; // Consistent with Zmm to 0.2% after mass fit and fix
      // Use same source for both MPF and pT balance
      h2->SetName(Form("bm%d_scale_%03.0f_zeejet_%d",
		       (1<<(n0+izee) | (1<<(n1+izee))), escale*10000, i));
      is_zee = hs.size();
    }
    if (ss=="zmmjet" && sm=="ptchs") {
      escale = 0.002; // Z mass fit within 0.2% of unity, flat vs pT
      // Use same source for both MPF and pT balance
      h2->SetName(Form("bm%d_scale_%03.0f_zmmjet_%d",
		       (1<<(n0+izmm) | (1<<(n1+izmm))), escale*10000, i));
      is_zmm = hs.size();
    }
    if (ss=="zlljet" && sm=="ptchs") {
      escale = 0.0020; // Zmm mass within 0.2% of unity, Zee fixed to it
      // Use same source for both MPF and pT balance
      //h2->SetName(Form("bm%d_scale_%03.0f_zlljet_%d",
      //	       (1<<(n0+izll) | (1<<(n1+izll))), escale*10000, i));
      h2->SetName(Form("bm%d_scale_%03.0f_zlljet",
		       (1<<(n0+izll) | (1<<(n1+izll))), escale*10000));
      is_zll = hs.size();
    }
    if (ss=="zjet" && sm=="ptchs") {
      escale = 0.0020; // Zmm mass within 0.2% of unity, Zee fixed to it
      // Use same source for both MPF and pT balance
      //h2->SetName(Form("bm%d_scale_%03.0f_zjet_%d",
      //	       (1<<(n0+izll) | (1<<(n1+izll))), escale*10000, i));
      h2->SetName(Form("bm%d_scale_%03.0f_zjet",
		       (1<<(n0+izll) | (1<<(n1+izll))), escale*10000));
      is_zll = hs.size();
    }

    // Same scale uncertainty applies to all pT bins
    for (int j = 1; j != h2->GetNbinsX()+1; ++j) {
      h2->SetBinContent(j, escale);
      h2->SetBinError(j, escale);
    } // for j

    if (escale!=0) hs.push_back(h2);
  } // for i

  // Create additional sources for MPF uncertainties with e/gamma
  // (one for each sample x method, but all except two are empty)
  for (unsigned int i = 0; i != _vdt1->size(); ++i) {

    unsigned int isample = i%nsamples;
    string s = samples[isample];
    const char *cs = s.c_str();
    unsigned int imethod = i/nsamples;
    string m = methods[imethod];
    const char *cm = m.c_str();

    double escale(0);
    if (s=="gamjet" && m=="mpfchs1") { escale = 0.002;} 
    if (s=="zeejet" && m=="mpfchs1") { escale = 0.002; }
    if (s=="zlljet" && m=="mpfchs1") { escale = 0.002; }
    if (s=="zjet"   && m=="mpfchs1") { escale = 0.002; }

    TH1D *h = hs[i]; assert(h);
    //TH1D *h2 = (TH1D*)h->Clone(Form("bm%d_%s_%02.0f_%s_%s_%d",
    TH1D *h2 = (TH1D*)h->Clone(Form("bm%d_%s_%02.0f_%s_%s",
				    1<<i,
				    escale!=0 ? "mpfscale" : "inactive",
				    //escale*1000.,cm,cs,i));
				    escale*1000.,cm,cs));
    
    for (int j = 1; j != h2->GetNbinsX()+1; ++j) {
      h2->SetBinContent(j, escale);
      h2->SetBinError(j, escale);
    } // for j
    if (escale!=0) hs.push_back(h2);
  } // for i

  // Additional uncertainties for pT balance shape so that
  // it is not given more weight than MPF with footprint uncertainty
  // (one for each sample x method, but all except two are empty)
  // With dofsrcombo: switch off again
  if (!fsrcombo) {
    for (unsigned int i = 0; i != _vdt1->size(); ++i) {

      unsigned int isample = i%nsamples;
      string ss = samples[isample];
      const char *cs = ss.c_str();
      unsigned int imethod = i/nsamples;
      string sm = methods[imethod];
      const char *cm = sm.c_str();
      
      double escale(0);
      if (ss=="gamjet" && sm=="ptchs") { escale = 0.005; }
      if (ss=="zeejet" && sm=="ptchs") { escale = 0.005; }
      if (ss=="zmmjet" && sm=="ptchs") { escale = 0.005; }
      if (ss=="zlljet" && sm=="ptchs") { escale = 0.005; }
      if (ss=="zjet"   && sm=="ptchs") { escale = 0.005; }
      
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
  } // !fsrcombo

  // UL17-V2: add l2cos x <rho>/<rho>_ref as low pT systematic
  if (false) {

    // Derived with help of minitools/rhoBias.C
    TF1 *fl2 = new TF1("fl1","1-([0]+[1]*log(x)+[2]*pow(log(x),2))/x"
		       "*log(x)*pow(1+x/[3],[4])",10,3500);
    fl2->SetParameters(1.57597,-1.89194,0.383362, 158.4,-4.812);
    
    for (unsigned int i = 0; i != _vdt1->size(); ++i) {
      
      //TH1D *hl1 = (TH1D*)d->Get("hl1bias"); assert(hl1);
      TH1D *hl1 = (TH1D*)d->Get("hl2cos"); assert(hl1);
      TH1D *h2 = (TH1D*)hl1->Clone(Form("bm%d_l1bias_%d",1<<i,i));
      
      bool usesrc(false);
      for (int j = 1; j != h2->GetNbinsX()+1; ++j) {
	int n0 = nsample0;
	int n1 = nsamples+nsample0;
	if (int(i)==n0+izll) {
	  double pt = h2->GetBinCenter(j);
	  h2->SetBinContent(j, 1-fl2->Eval(pt));
	  h2->SetBinError(j, fabs(1-fl2->Eval(pt)));
	  
	  // PileUpPt source applies to all samples and both MPF and pT balance
	  h2->SetName(Form("bm%d_l1bias_zlljet_%d",(1<<(n0+izll)|1<<(n1+izll)),i));
	  usesrc = true;
	}
	else {
	  h2->SetBinContent(j, 0);
	  h2->SetBinError(j, 0);
	}
      } // for j
      if (usesrc) hs.push_back(h2);
    } // for i
  }

  // Uncertainty sources for multijets
  if (sampleset.find("multijet")!=sampleset.end()) {

    const int nsrcmj = 0;//1;//3;
    for (unsigned int imethod = 0; imethod != nmethods; ++imethod) {
      for (int isrc = 0; isrc != nsrcmj; ++isrc) {

	string ss = "multijet";
	const char *cs = ss.c_str();
	string sm = methods[imethod];
	const char *cm = sm.c_str();
	
	string s;
	if (isrc==0) s = "sys/jer_multijet";
	TH1D *h = (TH1D*)d->Get(s.c_str());
	if (fsrcombo && !isRefIOV) h = (TH1D*)dfsr->Get(s.c_str());
	if (!h) cout << "Histo "<<s<<" not found!" << endl << flush;
	assert(h);
	
	// JER uncertainty applies for both MPF and MJB
	if (isrc==0 && sm=="ptchs") { // jer
	  h->SetName(Form("bm%d_multijet_jer_src%d_%d",
			  ((1<<0) | (1<<nsamples)),isrc,int(hs.size()+1)));
	  hs.push_back(h);
	}
      } // for imethod
    } // for ieig
  } // multijet syst.

  // Uncertainty sources for inclusive jets
  if (sampleset.find("incjet")!=sampleset.end() && useIncJetRefIOVUncert) {
    
    const int iij(1);
    assert(string(samples[iij])=="incjet");

    TH1D *h(0); int neig(0);
    while ( (h = (TH1D*)dsys->Get(Form("sys/hjesfit_eig%d",neig))) ) {
      h->SetName(Form("bm%d_incjet_jer_src%d_%d",
		      ((1<<(0+iij)) | (1<<(nsamples+iij))),neig,
		      int(hs.size()+1)));
      hs.push_back(h);
      ++neig;
      /*
      // Add separate ones for fwd2 and dag (different JER)
      h->SetName(Form("bm%d_incjet_jer_src%d_%d",
		      (1<<(0+iij)),neig,hs.size()+1));
      hs.push_back(h);
      ++neig;
      h->SetName(Form("bm%d_incjet_jer_src%d_%d",
		      ((1<<(nsamples+iij))),neig,hs.size()+1));
      hs.push_back(h);
      ++neig;
      */
    }
  } // incjet syst.

  // Uncertainty sources for hadronic W
  if (sampleset.find("hadw")!=sampleset.end()) {

    const int n0 = nsample0;
    const int n1 = nsamples+nsample0;

    assert(string(samples[n0+iwmh])=="hadw");

    string sa = "sys/hadw_ptave_fitprob";
    TH1D *ha = (TH1D*)dsys->Get(sa.c_str());
    if (!ha) cout << "Histo "<<sa<<" not found!" << endl << flush;
    assert(ha);

    string sb = "sys/hadw_ptboth_fitprob";
    TH1D *hb = (TH1D*)dsys->Get(sb.c_str());
    if (!hb) cout << "Histo "<<sb<<" not found!" << endl << flush;
    assert(hb);

    string sa2 = "sys/hadw_ptave_fitprob2";
    TH1D *ha2 = (TH1D*)dsys->Get(sa2.c_str());
    if (!ha2) cout << "Histo "<<sa2<<" not found!" << endl << flush;
    assert(ha2);

    string sb2 = "sys/hadw_ptboth_fitprob2";
    TH1D *hb2 = (TH1D*)dsys->Get(sb2.c_str());
    if (!hb2) cout << "Histo "<<sb2<<" not found!" << endl << flush;
    assert(hb2);

    // JER uncertainty applies for both MPF and MJB
    //h->SetName(Form("bm%d_hadw_fitprob",
    //		    ((1<<0) | (1<<nsamples)))); // bug!
    //   	    ((1<<0+iwmh) | (1<<nsamples+iwmh)))); // fixed

    // Separate systematics for pTave and pTboth
    ha->SetName(Form("bm%d_hadw_ptave_fitprob",1<<(n0+iwmh)));
    hs.push_back(ha);
    hb->SetName(Form("bm%d_hadw_ptboth_fitprob",1<<(n1+iwmh)));
    hs.push_back(hb);

    // Separate systematics for pTave and pTboth
    ha2->SetName(Form("bm%d_hadw_ptave_fitprob2",1<<(n0+iwmh)));
    hs.push_back(ha2);
    hb2->SetName(Form("bm%d_hadw_ptboth_fitprob2",1<<(n1+iwmh)));
    hs.push_back(hb2);
  } // hadw syst.

  // Uncertainty sources for gamma+jet EM energy scale from Zee fit
  if (true) {

    const int n0 = nsample0;
    const int n1 = nsamples+nsample0;

    const int nees = 3;
    for (int ieig = 0; ieig != nees; ++ieig) {
	
	string s = Form("sys/zee_gamjet_eig%d",ieig);
	TH1D *h = (TH1D*)d->Get(s.c_str());
	if (!h) cout << "Histo "<<s<<" not found!" << endl << flush;
	assert(h);
	
	// EES uncertainty applies for both MPF and MJB
	h->SetName(Form("bm%d_eesfromzee_gamjet_eig%d",
			((1<<(n0+igj)) | (1<<(n1+igj))),ieig));
	hs.push_back(h);
    } // for ieig
  } // Photon EES

  // use global pointer to give outside access
  _vsrc = &hs;

  
  /////////////////////////
  // Draw original data  //
  /////////////////////////
  cout << "Draw original data" << endl;

  const int maxpt = 4500;
  const int minpt = 15;
  TH1D *h = new TH1D("h",";p_{T} (GeV);Jet response (ratio)",
		     maxpt-minpt,minpt,maxpt);
  h->SetMinimum(etamin>=3 ? 0.50 : (etamin>=2.5 ? 0.70 : 0.9301));
  h->SetMaximum(etamin>=3 ? 1.75 : (etamin>=2.5 ? 1.45 : 1.07));
  if (_useZoom) {
    //less aggressive zoom, depending on detector region
    h->SetMinimum(etamin>=3 ? 0.80 : (etamin>=2.5 ? 0.80 : 0.9401));
    h->SetMaximum(etamin>=3 ? 1.30 : (etamin>=2.5 ? 1.30 : 1.0599));
  }
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  h->DrawClone("AXIS");

  map<string, const char*> lumimap;
  lumimap["BCDEF"] = "2017, 41.5 fb^{-1}"; // for DP note
  lumimap["B"] = "Run2017B, 4.8 fb^{-1}";
  lumimap["C"] = "Run2017C, 9.6 fb^{-1}";
  lumimap["D"] = "Run2017D, 4.2 fb^{-1}";
  lumimap["E"] = "Run2017E, 9.3 fb^{-1}";
  lumimap["F"] = "Run2017F, 13.4 fb^{-1}";
  lumimap["2018A"] = "2018A"; // placeholder
  lumimap["2018B"] = "2018B"; // placeholder
  lumimap["2018C"] = "2018C"; // placeholder
  lumimap["2018D"] = "2018D"; // placeholder
  lumimap["2018ABCD"] = "2018, 59.8 fb^{-1}"; // placeholder
  lumi_13TeV = lumimap[epoch];

  TCanvas *c0 = tdrCanvas("c0",h,4,11,true);
  gPad->SetLogx();
  
  // multijets up/down
  TLegend *legpf = tdrLeg(0.54,0.805,0.81,0.845);
  legpf->AddEntry(gs1[0]," ","PL");
  if (_vpt2->size()!=0) legpf->AddEntry((*_vpt2)[0]," ","PL");
  TLegend *legmf = tdrLeg(0.62,0.805,0.90,0.845);
  if(nmethods==2)legmf->AddEntry(gs1[nsamples]," ","PL");
  if (_vpt2->size()>1) legmf->AddEntry((*_vpt2)[1]," ","PL");

  TLegend *legp = tdrLeg(0.54,_useZoom ? 0.65 : 0.60,0.81,0.90);
  if( (nmethods==1||nmethods==2) && strcmp(methods[0],"ptchs")  ==0 ){
    legp->SetHeader("p_{T}^{bal}");
    for (int i = 0; i != nsamples; ++i)
      legp->AddEntry(gs1[i]," ",i==0 ? "" : "PL");
  }
  map<string, const char*> texlabel;
  texlabel["multijet"] = "Multijet";
  texlabel["incjet"] = "Incl. jet";
  texlabel["hadw"] = "W>qq'";
  texlabel["dijet"] = "Dijet";
  texlabel["gamjet"] = "#gamma+jet";
  texlabel["zeejet"] = "Zee+jet";
  texlabel["zmmjet"] = "Z#mu#mu+jet";
  texlabel["zlljet"] = "Z+jet";
  texlabel["zjet"] = "Z+jet(SL)";

  TLegend *legm = tdrLeg(0.62,_useZoom ? 0.65 : 0.60,0.90,0.90);
  if( (nmethods==1&&strcmp(methods[0],"mpfchs1")  ==0) || (nmethods==2 && strcmp(methods[1],"mpfchs1")  ==0 )){
    legm->SetHeader("MPF");
    for (int i = 0; i != nsamples; ++i)
      legm->AddEntry(gs1[i+(strcmp(methods[0],"mpfchs1")==0 ? 0 : nsamples)],texlabel[samples[i]],i==0 ? "" : "PL");
  }

  TLegend *legi = tdrLeg(0.62,0.17,0.90,0.21);

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  if (isl3) {
    if (useZJet100) 
      tex->DrawLatex(0.20,0.73,Form("|#eta|<%1.1f, #alpha<1.0%s",etamax,
				    dofsr ? "#rightarrow0" : ""));
    else if (useZJet50) 
      tex->DrawLatex(0.20,0.73,Form("|#eta|<%1.1f, #alpha<0.5%s",etamax,
				    dofsr ? "#rightarrow0" : ""));
    else
      tex->DrawLatex(0.20,0.73,Form("|#eta|<%1.1f, #alpha<0.3%s",etamax,
				    dofsr ? "#rightarrow0" : ""));
  }
  else {
    assert(dofsr);
    tex->DrawLatex(0.20,0.73,Form("%1.3f#leq|#eta|<%1.3f",etamin,etamax));
  }
  tex->DrawLatex(0.20, _useZoom ?  0.17 : 0.22, "Before global fit");

  hrun1->SetLineWidth(2);
  hrun1->SetLineColor(kCyan+3);
  hrun1->SetLineStyle(kDashed);
  hrun1->SetFillColorAlpha(hrun1->GetFillColor(),0.8);

  herr->SetLineWidth(2);
  herr->SetLineColor(kCyan+3);
  herr->SetLineStyle(kDashed);
  herr->SetFillColor(kCyan+1);
  herr->SetFillColorAlpha(herr->GetFillColor(),0.8);

  herr_ref->SetLineWidth(2);
  herr_ref->SetLineColor(kYellow+3);
  herr_ref->SetLineStyle(kDashed);
  herr_ref->SetFillColorAlpha(herr_ref->GetFillColor(),0.8);

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

  if (plotIncJet) {
    //tdrDraw(hij,"Pz",kFullDiamond,kOrange+2);//kBlack);
    //legi->AddEntry(hij,"Incl. jet","PLE");
    //hij->SetMarkerSize(0.5);
    tdrDraw(gij,"Pz",kFullDiamond,kOrange+2);//kBlack);
    gij->SetMarkerSize(0.8);
    legi->AddEntry(gij,"Incl. jet","PLE");
  }

  for (unsigned int i = 0; i != _vdt1->size(); ++i) {

    TGraphErrors *g2 = (*_vdt2)[i];

    string ss = samples[i%nsamples];
    string sm = methods[i/nsamples];
    
    // Add multijet downward points
    if (ss=="multijet") {
      int imethod = i/nsamples;
      TGraphErrors *gf2 = (_vpt2->size()==nmethods ? (*_vpt2)[imethod] : 0);
      if (gf2) {
	gf2->SetMarkerColor(kGray+2);
	gf2->SetLineColor(kGray+2);
	for (int j = gf2->GetN()-1; j != -1; --j)
	  if (gf2->GetX()[j]>ptmaxMultijetDown) gf2->RemovePoint(j);
	if (plotMultijetDown) {
	  TGraphErrors *gf2tmp = (TGraphErrors*)gf2->DrawClone("SAMEPz");
	  gf2tmp->SetMarkerSize(0.5);
	}
      }
    }

    // Shift pT balance points for better visibility of error bars
    if (sm=="ptchs" && shiftPtBal) {
      for (int j = 0; j != g2->GetN(); ++j) {
	g2->SetPoint(j, g2->GetX()[j]*shiftPtBal, g2->GetY()[j]);
      }
    }
    // And also make point smaller for better visibility of error bars
    if (ss=="multijet") { g2->SetMarkerSize(0.5); }
    if (ss=="gamjet")   { g2->SetMarkerSize(0.7); }
    if (ss=="zlljet")   { g2->SetMarkerColor(kRed); g2->SetLineColor(kRed); }
    if (ss=="zjet")     { g2->SetMarkerColor(kRed); g2->SetLineColor(kRed); }
    if (ss=="hadw")     { g2->SetMarkerSize(0.5); }
    
    g2->DrawClone("SAMEPz");
  }


  if (!_useZoom) {
    legp->AddEntry(hrun1," ","");
    legm->AddEntry(hrun1,"Run I","FL");
    legp->AddEntry(herr_ref," ","");
    //if (!_paper) legm->AddEntry(herr_ref,"UL17V2M5 IOV","FL");
    if (!_paper) legm->AddEntry(herr_ref,"UL17V4 IOV","FL");
    if ( _paper) legm->AddEntry(herr_ref,"Run II","FL");
  }
  else {
    legp->AddEntry(herr," ","FL");
    //if (!_paper) legm->AddEntry(herr_ref,"UL17V2M5","FL");
    if (!_paper) {
      if (isUL18)
	legm->AddEntry(herr_ref,"UL17V4CDE","FL");
      else
	legm->AddEntry(herr_ref,"UL17V4","FL");
    }
    if ( _paper) legm->AddEntry(herr_ref,"Syst. (tot,abs)","FL");
   }


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
  legi->Draw();

  if (etamin==0 && (fabs(etamax-1.3)<0.1 || fabs(etamax-2.4)<0.1)) {
    if (epoch!="L4") {
      if (useZJet100) 
	tex->DrawLatex(0.20,0.73,"|#eta|<1.3, #alpha<1.0");
      else if (useZJet50) 
	tex->DrawLatex(0.20,0.73,"|#eta|<1.3, #alpha<0.5");
      else
	tex->DrawLatex(0.20,0.73,"|#eta|<1.3, #alpha<0.3");
    }
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

  if (plotIncJet) {
    //tdrDraw(hij,"Pz",kFullDiamond,kOrange+2);//kBlack);
    //hij->SetMarkerSize(0.5);
    tdrDraw(gij,"Pz",kFullDiamond,kOrange+2);//kBlack);
  }

  for (unsigned int i = 0; i != _vdt3->size(); ++i) {
    
    TGraphErrors *g3 = (*_vdt3)[i];

    string ss = samples[i%nsamples];
    string sm = methods[i/nsamples];

    // Add multijet downward points
    if (ss=="multijet") {
      int imethod = i/nsamples;
      TGraphErrors *gf3 = (_vpt3->size()==nmethods ? (*_vpt3)[imethod] : 0);
      if (gf3) {
	gf3->SetMarkerColor(kGray+2);
	gf3->SetLineColor(kGray+2);
	gf3->SetMarkerSize(0.5); // UL17
	for (int j = gf3->GetN()-1; j != -1; --j)
	  if (gf3->GetX()[j]>ptmaxMultijetDown) gf3->RemovePoint(j);
	gf3->DrawClone("SAMEPz");
      }
    }

    // Shift pT balance points for better visibility of error bars
    if (sm=="ptchs" && shiftPtBal) {
      for (int j = 0; j != g3->GetN(); ++j) {
	g3->SetPoint(j, g3->GetX()[j]*shiftPtBal, g3->GetY()[j]);
      }
    }
    // And also make them smaller for better visibility of error bars
    if (ss=="multijet") { g3->SetMarkerSize(0.5); }
    if (ss=="gamjet")   { g3->SetMarkerSize(0.7); }
    if (ss=="zlljet")   { g3->SetMarkerColor(kRed); g3->SetLineColor(kRed); }
    if (ss=="zjet")     { g3->SetMarkerColor(kRed); g3->SetLineColor(kRed); }
    if (ss=="hadw")     { g3->SetMarkerSize(0.5); }

    // Skip multijet downward points
    g3->DrawClone("SAMEPz");
  }

  
  ///////////////////////
  // Perform global fit
  //////////////////////
  cout << "Perform global fit" << endl;
  
  // Fit function
  TF1 *jesfit = new TF1("jesfit",jesFit,minpt,maxpt,njesFit);
  jesfit->SetLineColor(kBlack);
  if (njesFit==1)  jesfit->SetParameter (0, 0.98);
  if (njesFit==2)  jesfit->SetParameters(0.99, 0.05);
  //if (njesFit==3)  jesfit->SetParameters(0.981, 0.044, -0.519);
  //if (njesFit==4)  jesfit->SetParameters(0.981, 0.044, -0.519,0 );
  //if (njesFit==3)  jesfit->SetParameters(0.9811, 1.0868, -0.1813);
  if (njesFit==3)  jesfit->SetParameters(0.98, 1.0, 0.);
  if (njesFit==4)  jesfit->SetParameters(0.9811, 1.0868, -0.1813, 0);
  if (njesFit==5)  jesfit->SetParameters(2., 0.,2., 3., 0.);
  //if (njesFit==6)  jesfit->SetParameters(2., 1.,epoch=="B"?1.:0.,1., 1., 0.);
  //if (njesFit==6)  jesfit->SetParameters(1, 2., 1.,epoch=="B"?1.:0.,1., 0.);
  //if (njesFit==6)  jesfit->SetParameters(0.572,1.889,1.213,-0.029,1.243,-0.128);
  //if (njesFit==6)  jesfit->SetParameters(0.1,0.1,0.1,0.1,0.1,0.1);
  //if (njesFit==6)  jesfit->SetParameters(0.01,0.01,0.01,0.01,0.01,0.01);
  if (njesFit==6)  jesfit->SetParameters(1,1,1,0,0,0);
  //if (njesFit==7)  jesfit->SetParameters(1., 1.,0.,0.,1., 1., 0.);
  //if (njesFit==7)  jesfit->SetParameters(1.5, 1.1,0.1,-0.1,1.4, 1.4, -0.1);
  //if (njesFit==7)  jesfit->SetParameters(1.57,1.10,0.10,-0.08,1.40,1.33,-0.18);
  //if (njesFit==7)  jesfit->SetParameters(1,1,1,0,0,-1,0);
  //if (njesFit==7)  jesfit->SetParameters(1,1,1,0,1,-1,0);

  //if (njesFit==8)  jesfit->SetParameters(1,1,1,0,1,-1,0,1);
  //if (njesFit==9)  jesfit->SetParameters(0,0,0,0,0,0,0,0,0);
  if (njesFit==9)  jesfit->SetParameters(1,1,1,1,1,0,0,0,0);

  _jesFit = jesfit;
  

  // Get linear equations and solve them to get initial values for fitter
  const int np = _jesFit->GetNpar();
  const int nsrc = _vsrc->size();
  Int_t Npar = np+nsrc;

  cout << "Global fit has " << np << " fit parameters and "
       << nsrc << " nuisance parameters." << endl;

  vector<double> a(Npar, 0);
  for (int i = 0; i != np; ++i) a[i] = jesfit->GetParameter(i);

  // Setup global chi2 fit (jesFitter is our function)
  TFitter *fitter = new TFitter(Npar);
  fitter->SetFCN(jesFitter);

  // Set parameters
  vector<string> parnames(Npar);
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
  TVectorD vpar(Npar);
  TVectorD verr(Npar);
  Double_t chi2_gbl(0), chi2_src(0), chi2_data(0);
  vector<double> vchi2_data(gs2.size(),0);
  vector<int> vndata(gs2.size(),0);
  int nsrc_true(0);
  Double_t grad[Npar];
  Int_t flag = 1;
  TH1D *hsrc = new TH1D("hsrc",";Nuisance parameter;",12,-3,3);

  for (int i = 0; i != Npar; ++i) {
    tmp_par[i] = fitter->GetParameter(i);
    tmp_err[i] = fitter->GetParError(i);
    vpar[i] = fitter->GetParameter(i);
    verr[i] = fitter->GetParError(i);
  }
  jesFitter(Npar, grad, chi2_gbl, tmp_par, flag);
  cout << "List of fitted sources that count (nonzero):" << endl;
  for (int i = np; i != Npar; ++i) {
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

  if (storeErrorMatrix) {
    TFile *femat = new TFile("rootfiles/globalFitL3Res_emat.root","RECREATE");
    TH2D *h2emat = new TH2D("h2emat",";i_{par};j_{par}",
			    Npar,-0.5,Npar-0.5, Npar,-0.5,Npar-0.5);
    TH2D *h2cov = new TH2D("h2cov",";i_{par};j_{par}",
			   Npar,-0.5,Npar-0.5, Npar,-0.5,Npar-0.5);
    TH1D *h1par = new TH1D("h1par",";i_{par};",Npar,-0.5,Npar-0.5);
    for (int i = 0; i != Npar; ++i) {

      h1par->SetBinContent(i+1, vpar[i]);
      h1par->SetBinError(i+1, verr[i]);
      for (int j = 0; j != Npar; ++j) {
	h2emat->SetBinContent(i+1, j+1, emat[i][j]);
	h2cov->SetBinContent(i+1, j+1, emat[i][i]*emat[j][j]==0 ? 0 :
			     emat[i][j] / sqrt(emat[i][i]*emat[j][j]));
      } // for j
      string sp = (i<np ? Form("p%d",i) : (*_vsrc)[i-np]->GetName());
      const char *cp = sp.c_str();
      //h2emat->GetYaxis()->SetBinLabel(i+1,Form("%s (%d)",cp,i));
      //h2cov->GetYaxis()->SetBinLabel(i+1,Form("%s (%d)",cp,i));
      //h1par->GetXaxis()->SetBinLabel(i+1,Form("%s (%d)",cp,i));
      h2emat->GetYaxis()->SetBinLabel(i+1,Form("%s",cp));
      h2cov->GetYaxis()->SetBinLabel(i+1,Form("%s",cp));
      h1par->GetXaxis()->SetBinLabel(i+1,Form("%s",cp));
    } // for i
    // Store to /sys folder
    if (storeJesFit || isRefIOV) {
      dout->cd("sys");
      emat.Write("emat",TObject::kOverwrite);
      vpar.Write("vpar",TObject::kOverwrite);
      verr.Write("verr",TObject::kOverwrite);
      h2emat->Write("h2emat",TObject::kOverwrite);
      h2cov->Write("h2cov",TObject::kOverwrite);
      h1par->Write("h1par",TObject::kOverwrite);
    }
    // And to separate file
    femat->cd();
    emat.Write("emat");
    vpar.Write("vpar");
    verr.Write("verr");
    h2emat->Write();
    h2cov->Write();
    h1par->Write();
    femat->Close();
    curdir->cd();
  }

  cout << endl;
  cout << "*** Processing eta bin " << etamin<<" - "<<etamax << " ***" << endl;
  cout << "Global chi2/ndf = " << chi2_gbl
       << " / " << Nk-np << " for " << Nk <<" data points, "
       << np << " fit parameters and "
       << nsrc << " ("<<nsrc_true<<") uncertainty sources." << endl;
  cout << endl;

  cout << endl;
  cout << "For data chi2/ndf = " << chi2_data << " / " << Nk << endl;
  cout << "For sources chi2/ndf = " << chi2_src << " / " << nsrc_true << endl;
  cout << "Per data set:" << endl;
  for (unsigned int i = 0; i != vchi2_data.size(); ++i) {
    cout << "  " << Form("%1.1f",vchi2_data[i]) << " / " << vndata[i]
	 << "  ("<<gs2[i]->GetName()<<")"<< endl;
  }

  /*
  // Printing code not recently updated
  // Can use minitools/createL2L3ResTextFile.C also
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
      if (string(cs)=="incjet") continue;
      string s = Form("fsr/hkfsr_%s_%s",cm,cs);
      TH1D *h = (TH1D*)d->Get(s.c_str());
      if (dofsrcombo2 && epoch!="BCDEF") s = Form("fsr/hkfsr2_%s_%s",cm,cs);
      if (fsrcombo) h = (TH1D*)dfsr->Get(s.c_str());
      if (!h) cout << "Histo "<<s<<" not found!" << endl << flush;
      assert(h);

      ofstream txtFSRDiJet(Form("txt2/GlobalFitOutput_FSR_%s_%s.txt",cs,cm),ios_base::app);
      if(etamin==0.&&etamax==0.261)txtFSRDiJet << "{ 2 JetEta  JetPt 1 JetPt [0] Correction L2Relative}";
      for (int i = 1; i != h->GetNbinsX()+1; ++i) {
        double ptlow = h->GetBinLowEdge(i);
        double ptup = h->GetBinLowEdge(i+1);
        double FSR = h->GetBinContent(i);
        if(np==2&& !(etamin==0.&&etamax==1.3)){
          txtFSRDiJet      << Form("\n %7.4f   %7.4f  %7.0f   %7.0f  3 10 6500 %7.4f ", etamin, etamax, ptlow, ptup, FSR*0.3);
        }
      }
    }
  }
  */
  
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
  for (int i = np; i != Npar; ++i) {
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
  legi->Draw();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(h->GetMinimum(),
	      6500./cosh(etamin),
	      h->GetMinimum()+0.5*(h->GetMaximum()-h->GetMinimum()),
	      6500./cosh(etamin));

  if (etamin==0 && (fabs(etamax-1.3)<0.1 || fabs(etamax-2.4)<0.1)) {
    if (epoch!="L4") {
      if (useZJet100)
	tex->DrawLatex(0.20,0.73,"|#eta|<1.3, #alpha<1.0#rightarrow0");
      else if (useZJet50)
	tex->DrawLatex(0.20,0.73,"|#eta|<1.3, #alpha<0.5#rightarrow0");
      else
	tex->DrawLatex(0.20,0.73,"|#eta|<1.3, #alpha<0.3#rightarrow0");
    }
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

  if (plotIncJet) {
    //tdrDraw(hij,"Pz",kFullDiamond,kOrange+2);//kBlack);
    //hij->SetMarkerSize(0.5);
    tdrDraw(gij,"Pz",kFullDiamond,kOrange+2);//kBlack);
  }

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
    string ss = samples[i%nsamples];
    string sm = methods[i/nsamples];

    // Add multijet downward points
    if (ss=="multijet") {
      unsigned int imethod = i/nsamples;
      TGraphErrors *gf2 = (_vpt2->size()==nmethods ? (*_vpt2)[imethod] : 0);
      if (gf2) {
	gf2->SetMarkerColor(kGray+2);
	gf2->SetLineColor(kGray+2);
	for (int j = gf2->GetN()-1; j != -1; --j)
	  if (gf2->GetX()[j]>ptmaxMultijetDown) gf2->RemovePoint(j);
	if (sm=="ptchs" && shiftPtBal) {
	  for (int j = 0; j != gf2->GetN(); ++j) {
	    gf2->SetPoint(j, gf2->GetX()[j]*shiftPtBal, gf2->GetY()[j]);
	  }
	}
	if (plotMultijetDown) {
	  TGraphErrors *gf2tmp = (TGraphErrors*)gf2->DrawClone("SAMEPz");
	  gf2tmp->SetMarkerSize(0.5); // UL17
	}
      }
    }
    
    if (ss=="gamjet")
      g2->SetLineWidth(3);

    // Clean out large uncertainties
    if (_cleanUncert) {
      for (int i = g2->GetN()-1; i != 0; --i) {
	if (g2->GetEY()[i]>_cleanUncert) g2->RemovePoint(i);
      }
    }
    
    // Shift pT balance points for better visibility of error bars
    if (sm=="ptchs" && shiftPtBal) {
      for (int j = 0; j != g2->GetN(); ++j) {
    	g2->SetPoint(j, g2->GetX()[j]*shiftPtBal, g2->GetY()[j]);
      }
    }
    // And also make point smaller for better visibility of error bars
    if (ss=="multijet") { g2->SetMarkerSize(0.5); }
    if (ss=="gamjet")   { g2->SetMarkerSize(0.7); }
    if (ss=="zlljet")   { g2->SetMarkerColor(kRed); g2->SetLineColor(kRed); }
    if (ss=="zjet")     { g2->SetMarkerColor(kRed); g2->SetLineColor(kRed); }
    if (ss=="hadw")     { g2->SetMarkerSize(0.5); }

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

  if (!_paper && njesFit==9) {

    const char *name[] = {"trk","ph","hx","hh","he",
			  "hw","rc","td","mb"};

    //tex->DrawLatex(0.25,0.70,Form("Full 9-par model"));
    tex->DrawLatex(0.54,0.70-0.09,Form("Full 9-par model"));

    double ts_orig = tex->GetTextSize();
    double ts_new = 0.025;
    tex->SetTextSize(ts_new);
    for (int i = 0; i != njesFit; ++i) {
      tex->DrawLatex(0.30,0.70-i*ts_new,
		     Form("p%d(%s)=%1.3f #pm %1.3f",
			  i,name[i],
			  _jesFit->GetParameter(i),
			  sqrt(emat[i][i])));
    }
    tex->SetTextSize(ts_orig);

  }
  else if (!_paper) {
 
    //tex->SetTextColor(kWhite); // hide from view
    tex->DrawLatex(0.25,0.70,Form("Npar=%d",_jesFit->GetNpar()));
    tex->DrawLatex(0.25,0.67,Form("p0=%1.4f #pm %1.4f",
				  _jesFit->GetParameter(0),
				  sqrt(emat[0][0])));
    if (njesFit>=2)
      tex->DrawLatex(0.25,0.64,Form("p1=%1.4f #pm %1.4f",
				    _jesFit->GetParameter(1),
				    sqrt(emat[1][1])));
    if (njesFit==2 && fixOff)
      tex->DrawLatex(0.25,0.61,Form("p2=%1.4f (fixOff)",fixOff));
    if (njesFit>=3)
      tex->DrawLatex(0.25,0.61,Form("p2=%1.4f #pm %1.4f",
				  _jesFit->GetParameter(2),
				    sqrt(emat[2][2])));
    if (njesFit==3 && fixOff)
      tex->DrawLatex(0.25,0.58,Form("p3=%1.4f (fixOff)",fixOff));
    if (njesFit>=4)
      tex->DrawLatex(0.25,0.58,Form("p3=%1.4f #pm %1.4f",
				    _jesFit->GetParameter(3),
				    sqrt(emat[3][3])));
    if (njesFit>=5)
      tex->DrawLatex(0.25,0.55,Form("p4=%1.4f #pm %1.4f",
				    _jesFit->GetParameter(4),
				    sqrt(emat[4][4])));
    if (njesFit>=6)
      tex->DrawLatex(0.25,0.52,Form("p5=%1.4f #pm %1.4f",
				    _jesFit->GetParameter(5),
				    sqrt(emat[5][5])));
    if (njesFit>=7)
      tex->DrawLatex(0.25,0.49,Form("p6=%1.4f #pm %1.4f",
				    _jesFit->GetParameter(6),
				    sqrt(emat[6][6])));
    tex->SetTextSize(0.045); tex->SetTextColor(kBlack);
  }
  
  // Factorize error matrix into eigenvectors
  TVectorD eigvec(np);
  TMatrixD eigmat = emat2.EigenVectors(eigvec);

  // Eigenvectors
  vector<TH1D*> vhkeig(np);
  vector<TF1*> vfkeig(np);
  // Original components
  vector<TH1D*> vhcomp(np);
  vector<TF1*> vfcomp(np);
  for (int ieig = 0; ieig != np; ++ieig) {

    TF1 *fkeig = new TF1(Form("fitEig_feig%d",ieig), jesFit, minpt, maxpt, np);
    fkeig->SetNpx(1000);
    fkeig->SetLineStyle(kDashed);
    for (int i = 0; i != np; ++i) {
      fkeig->SetParameter(i, _jesFit->GetParameter(i)
			  + eigmat[i][ieig] * sqrt(eigvec[ieig]));
    }
    //fkeig->Draw("SAME"); // For visualizing uncertainty correlations vs pT

    TF1 *fcomp = new TF1(Form("fitComp_fpar%d",ieig), jesFit, minpt, maxpt, np);
    TF1 *fdiff = new TF1(Form("fitDiff_fpar%d",ieig), jesFit, minpt, maxpt, np);
    fcomp->SetNpx(1000);
    fcomp->SetLineStyle(kDotted);
    for (int i = 0; i != np; ++i) {
      fcomp->SetParameter(i, (i==ieig ? _jesFit->GetParameter(i) : 0));
      fdiff->SetParameter(i, (i==ieig ?
			      _jesFit->GetParameter(i)+sqrt(emat[i][i]) : 0));
    }

    TH1D *h = (TH1D*)hij->Clone(Form("hkeig_%d",ieig));
    h->Reset();
    for (int i = 1; i != h->GetNbinsX()+1; ++i) {
      double pt = h->GetBinCenter(i);
      h->SetBinContent(i, fkeig->Eval(pt) - _jesFit->Eval(pt));
    }

    TH1D *hc = (TH1D*)hij->Clone(Form("hcomp_%d",ieig));
    hc->Reset();
    for (int i = 1; i != hc->GetNbinsX()+1; ++i) {
      double pt = hc->GetBinCenter(i);
      hc->SetBinContent(i, fcomp->Eval(pt));
      hc->SetBinError(i, fabs(fdiff->Eval(pt)-fcomp->Eval(pt)));
    }

    vhkeig[ieig] = h;
    vhcomp[ieig] = hc;
    vfkeig[ieig] = fkeig;
    vfcomp[ieig] = fcomp;
  } // for ieig

  // Store fit with total uncertainty from eigenvectors
  TH1D *hjesfit = (TH1D*)hij->Clone("hjesfit");
  hjesfit->Reset();
  for (int i = 1; i != hjesfit->GetNbinsX()+1; ++i) {
    double sume2(0);
    for (int ieig = 0; ieig != np; ++ieig) {
      sume2 += pow(vhkeig[ieig]->GetBinContent(i),2);
    }
    double pt = hjesfit->GetBinCenter(i);
    hjesfit->SetBinContent(i, _jesFit->Eval(pt));
    hjesfit->SetBinError(i, sqrt(sume2));
  } // for i

  c0->RedrawAxis();
  c0b->RedrawAxis();
  c1->RedrawAxis();

  if (etamin==0 && (fabs(etamax-1.3)<0.1 || fabs(etamax-2.4)<0.1)) {
    if (dofsr) {
      c0b->SaveAs(Form("pdf/%s/globalFitL3res_raw.pdf",cep));
      c0->SaveAs(Form("pdf/%s/globalFitL3res_orig.pdf",cep));
      c1->SaveAs(Form("pdf/%s/globalFitL3res_shifted.pdf",cep));
      c1->SaveAs(Form("pdf/%s/globalFitL3res_shifted.root",cep));
    }
    if (writeTextFile) {
      ofstream fout("textFiles/globalFitL3Res.txt",ios::app);
      fout << Form("  if (set==\"%s SimpleL1\") { "
		   "p[0] = %1.4g, p[1] = %1.4g, p[2] = %1.4g; }\n",
		   cep,
		   _jesFit->GetParameter(0), _jesFit->GetParameter(1),
		   fixOff&&njesFit==2 ? fixOff : _jesFit->GetParameter(2));
      fout.close();
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
  // Also store output fit and eigenvectors

  int colorsdark[4]  = {kBlack,kBlue, kGreen+2,kRed};
  int colorslight[4] = {kGray, kBlue-10, kGreen+2-10,kRed-10};
  for (unsigned int imethod = 0; imethod != nmethods; ++imethod) {

    string sm = methods[imethod];
    const char *cm = sm.c_str();

    h->SetMinimum(imethod==0 ? -0.10 : -0.06);
    h->SetMaximum(imethod==0 ? +0.05 : +0.09);
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

    for (int isample = 0; isample != nsamples; ++isample) {

      string ss = samples[isample];
      const char *cs = ss.c_str();

      // No FSR correction for inclusive jets
      if (ss=="incjet") continue;
      // No FSR correction for hadronic W
      if (ss=="hadw") continue;

      int ibm = isample + nsamples*imethod;
      TH1D *hk = hks[ibm]; assert(hk);

      double minx(30);
      double maxx(1500);
      if (ss=="zeejet"||ss=="zmmjet"||ss=="zlljet") {  maxx = 700; }
      if (ss=="zjet") {  maxx = 700; }
      if (ss=="gamjet") { minx=230; maxx = 1500; }
      if (ss=="multijet") { minx=114; maxx = 2116; }//2640; }
      if (ss=="incjet") { minx=15; maxx = 2116; }//2640; }
      if (ss=="hadw") { minx=30; maxx = 200; }
      hk->GetXaxis()->SetRangeUser(minx,maxx);

      /*
      hk->SetLineColor(colorslight[isample]);
      hk->SetFillColor(colorslight[isample]);
      hk->SetFillStyle(kNone);
      hk->SetMarkerSize(0.);
      hk->DrawClone("SAME E5");
      hk->SetFillStyle(3002);
      hk->DrawClone("SAME E5");
      */
      hk->SetLineColor(colorslight[isample]);
      hk->SetFillColorAlpha(colorslight[isample],0.3);
      hk->SetFillStyle(kNone);
      hk->SetMarkerSize(0.);
      hk->DrawClone("SAME E5");
      hk->SetFillStyle(1001);//3002);
      hk->DrawClone("SAME E5");

      TGraph *gk = new TGraph(hk);
      for (int i = gk->GetN()-1; i != -1; --i) {
        if (verboseGF) {
	  cout << Form("%s %s pt: %4.0f fsr: %09.6f", samples[isample],
		       methods[imethod], gk->GetX()[i], gk->GetY()[i])
	       << endl << flush;
	}
	if (gk->GetX()[i] < minx || gk->GetX()[i] > maxx)
	  gk->RemovePoint(i);
        if (gk->GetY()[i]==0) gk->RemovePoint(i);
      }
      gk->Draw("SAMEL");

      leg1->AddEntry(hk, " ", "FL");
    
      // Sum over eigenvectors
      TH1D *hke = (TH1D*)hk->Clone(); hke->Reset();
      for (int ipt = 1; ipt != hke->GetNbinsX()+1; ++ipt) {
	double yeig(0);
	vector<double> df(Npar,0);
	for (int ieig = 0; ieig != neigmax; ++ieig) {

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
	string ss = samples[isample];
	const char *cs = ss.c_str();
	double ptreco(0);
	if (ss=="zeejet"||ss=="zmmjet"||ss=="zlljet") ptreco = ptreco_zlljet;
	if (ss=="zjet") ptreco = ptreco_zjet;
	if (ss=="gamjet") ptreco = ptreco_gjet;
	double aeff = min(max(alpha, ptreco/pt),alphaeffmax);
	if (ss=="multijet") aeff = alpha;
	hke->SetBinContent(ipt, (aeff*hk->GetBinContent(ipt) - yeig)/aeff);
	// BUG FIX 20200506: err2 is for dR at alpha, scale to dR/dalpha
	//hke->SetBinError(ipt, sqrt(err2));
	hke->SetBinError(ipt, sqrt(err2)/aeff);
      } // for ipt
      
      hke->GetXaxis()->SetRangeUser(minx,maxx);
      hke->SetFillStyle(1001);
      hke->SetLineColor(colorsdark[isample]);
      hke->SetFillColorAlpha(colorslight[isample],0.7);
      hke->DrawClone("SAME E5");

      TGraph *gke = new TGraph(hke);
      for (int i = gke->GetN()-1; i != -1; --i) {
        if (gke->GetX()[i] < minx || gke->GetX()[i] > maxx)
          gke->RemovePoint(i);
        if (gke->GetY()[i]==0.)gke->RemovePoint(i);
      }
      gke->Draw("SAMEL");

      leg2->AddEntry(hke, texlabel[samples[isample]], "FL");

      // Store scaled FSR eigenvectors and central value for BCDEF
      if (isRefIOV) {
	dout->cd("fsr");
	hke->Write(Form("hkfsr2_%s_%s",cm,cs),TObject::kOverwrite);
	for (int ieig = 0; ieig != neigmax; ++ieig) {
	  int ie = ibm+ieig*nsamples*nmethods;
	  TH1D *hi = hs[ie]; assert(hi);
	  hi = (TH1D*)hi->Clone(Form("hkfsr2_%s_%s_eig%d",cm,cs,ieig));
	  hi->Scale(tmp_par[np+ie]);
	  hi->Write(Form("hkfsr2_%s_%s_eig%d",cm,cs,ieig),TObject::kOverwrite);
	} // for ieig
	//
	curdir->cd();
      } // BCDEF
    } // for isample

    // Store also jesFit and it's eigenvectors for BCDEF
    // (using BCDEF also as input; update code later to also store output IOVs)
    if (storeJesFit || isRefIOV) {
      dout->cd("sys");
      hjesfit->Write("hjesfit",TObject::kOverwrite);
      for (int ieig = 0; ieig != np; ++ieig) {
	vhkeig[ieig]->Write(Form("hjesfit_eig%d",ieig),TObject::kOverwrite);
      }
      for (int ipar = 0; ipar != np; ++ipar) {
	vhcomp[ipar]->Write(Form("hjesfit_par%d",ipar),TObject::kOverwrite);
      }
      _jesFit->Write("fjesfit",TObject::kOverwrite);
      for (int ieig = 0; ieig != np; ++ieig) {
	vfkeig[ieig]->Write(Form("fjesfit_eig%d",ieig),TObject::kOverwrite);
      }
      for (int ipar = 0; ipar != np; ++ipar) {
	vfcomp[ipar]->Write(Form("fjesfit_par%d",ipar),TObject::kOverwrite);
      }
      eigvec.Write("eigvec",TObject::kOverwrite);
      eigmat.Write("eigmat",TObject::kOverwrite);
      curdir->cd();
    } // store jesFit and eigenvectors

    //if (etamin==0 && (fabs(etamax-1.3)<0.1 || fabs(etamax-2.4)<0.1)) {
    if (isl3) {
      c3->SaveAs(Form("pdf/%s/globalFitL3res_%s_kfsr.pdf",
		      cep,methods[imethod]));
    }
    else {
      c3->SaveAs(Form("pdf/%s/globalFitL3res_%s_kfsr_eta%02.0f-%02.0f.pdf",
		      cep,methods[imethod],10*etamin,10*etamax));
    }
  } // for imethod
  

  ///////////////////////////////////////////
  // Draw PF compositions and their fit
  //////////////////////////////////////////

  double ptmin(15), ptmax(4500);
  TH1D *h4 = tdrHist("hpfjet","PF composition changes (10^{-2})",
		     -1.2,1.3,"p_{T} (GeV)",ptmin,ptmax);
  if (isUL18) {
    h4->SetMinimum(-2.4);
    h4->SetMaximum(+2.6);
  }
  TCanvas *c4 = tdrCanvas("c4",h4,4,11,kSquare);
  gPad->SetLogx();
  TLegend *leg4 = tdrLeg(0.60,0.72,0.90,0.90);
  leg4->SetHeader("Dijet");
  TLegend *leg4z = tdrLeg(0.40,0.72,0.70,0.90);
  leg4z->SetHeader("Z+jet");
  TLine *l4 = new TLine(); l4->SetLineStyle(kDashed);
  l4->DrawLine(ptmin,0,ptmax,0);

  // Draw output fits with bands firsts
  //bool isrefz = (usePFZ&&fitPFZ || usePFZ&&!usePF)
  bool isrefz = (usePFZ && !usePF);
  vector< pair<string,TGraphErrors*> > &_vpfx = (isrefz ? _vpfz : _vpf);
  vector< pair<string,TGraphErrors*> > &_vpfx2 = (isrefz ? _vpfz2 : _vpf2);
  for (unsigned int ic = 0; ic != _vpfx.size(); ++ic) {
    string sf = _vpfx[ic].first;
    const char *cf = sf.c_str();
    TGraphErrors *g2 = _vpfx2[ic].second;
    TGraphErrors *gc2 = (TGraphErrors*)g2->Clone("gc2");
    // Turn units to percentages (10^-2)
    for (int i = 0; i != gc2->GetN(); ++i) {
      gc2->SetPoint(i, gc2->GetX()[i], 100.*gc2->GetY()[i]);
      gc2->SetPointError(i, gc2->GetEX()[i], 100.*gc2->GetEY()[i]);
    } // for i
    TGraphErrors *gc3 = (TGraphErrors*)gc2->Clone("gc3");
    for (int i = 0; i != gc3->GetN(); ++i) {
      gc3->SetPointError(i, gc3->GetEX()[i], 0);
      // Evaluate post-fit uncertainty band for composition
      if (njesFit==9) {
	double pt = gc3->GetX()[i];
	// Calculate composition bias at this pT	
	double df[njesFit];
	for (int ipar = 0; ipar != njesFit; ++ipar) {
	  if (_mpf[sf][ipar]!=0) { // parameter maps to composition change
	    TF1 *f1 = _mpf[sf][ipar]; // effect shape vs pt
	    df[ipar] = f1->Eval(pt);//*sqrt(emat[ipar][ipar]);
	  }
	  else
	    df[ipar] = 0;
	} // for ipar
	double sumerr2(0);
	for (int ipar = 0; ipar != njesFit; ++ipar) {
	  for (int jpar = 0; jpar != njesFit; ++jpar) {
	    sumerr2 += df[ipar]*df[jpar]*emat[ipar][jpar];
	  } // for jpar
	} // for ipar
	gc3->SetPointError(i, gc3->GetEX()[i], sqrt(sumerr2));
      } // njesFit==8
    } // for i
    // Remove points outside fit range
    TGraphErrors *gc4 = (TGraphErrors*)gc3->Clone("gc4");
    for (int i = gc4->GetN()-1; i != -1; --i) {
      double pt = gc4->GetX()[i];
      //if ((usePF && (pt<minPFpt || pt>maxPFpt)) ||
      //  (usePFZ && (pt<minPFZpt || pt>maxPFZpt)))
      if ((fitPF && (pt<minPFpt || pt>maxPFpt)) ||
	  (fitPFZ && (pt<minPFZpt || pt>maxPFZpt)))
	gc4->RemovePoint(i);
    } // for i

    tdrDraw(gc2,"E3",g2->GetMarkerStyle(),g2->GetMarkerColor()+1,kSolid,-1,
	    1001,g2->GetMarkerColor()-9);
    gc2->SetFillColorAlpha(g2->GetMarkerColor()-9,0.2);//0.7);
    tdrDraw(gc3,"E3",g2->GetMarkerStyle(),g2->GetMarkerColor()+1,kSolid,-1,
    	    1001,g2->GetMarkerColor()-9);
    gc3->SetFillColorAlpha(g2->GetMarkerColor()-9,0.4);
    tdrDraw(gc4,"E3",g2->GetMarkerStyle(),g2->GetMarkerColor()+1,kSolid,-1,
	    1001,g2->GetMarkerColor());
    gc4->SetFillColorAlpha(g2->GetMarkerColor(),0.7);
    tdrDraw((TGraphErrors*)gc2->Clone(),"LX",
	    kNone,g2->GetMarkerColor(),kSolid,-1,kNone);
  } // for ic

  // Draw input graphs with points on top
  for (unsigned int ic = 0; ic != _vpfz.size(); ++ic) {
    string sf = _vpfz[ic].first;
    const char *cf = sf.c_str();
    TGraphErrors *g = _vpfz[ic].second;
    TGraphErrors *gc = (TGraphErrors*)g->Clone();
    for (int i = 0; i != gc->GetN(); ++i) {
      gc->SetPoint(i, g->GetX()[i], 100.*g->GetY()[i]);
      gc->SetPointError(i, g->GetEX()[i], 100.*g->GetEY()[i]);
    }
    tdrDraw(gc,"Pz",g->GetMarkerStyle(),g->GetMarkerColor());
    gc->SetMarkerSize(0.8);
    leg4z->AddEntry(gc,cf,"PLE");
  }
  for (unsigned int ic = 0; ic != _vpf.size(); ++ic) {
    string sf = _vpf[ic].first;
    const char *cf = sf.c_str();
    TGraphErrors *g = _vpf[ic].second;
    TGraphErrors *gc = (TGraphErrors*)g->Clone();
    for (int i = 0; i != gc->GetN(); ++i) {
      gc->SetPoint(i, g->GetX()[i], 100.*g->GetY()[i]);
      gc->SetPointError(i, g->GetEX()[i], 100.*g->GetEY()[i]);
    }
    tdrDraw(gc,"Pz",g->GetMarkerStyle(),g->GetMarkerColor());
    gc->SetMarkerSize(0.8);
    leg4->AddEntry(gc,cf,"PLE");
  }

  if (isl3) {
    c4->SaveAs(Form("pdf/%s/globalFitL3res_pfjet.pdf",cep));
  }
  else
    c4->SaveAs(Form("pdf/%s/globalFitL3res_pfjet_eta%02.0f-%02.0f.pdf",
		    cep,10*etamin,10*etamax));

  curdir->cd();
  if (ffsr) {
    ffsr->Close();
  }
  //if (fjes) {
  //fjes->Close();
  //}
  if (fsys) {
    fsys->Close();
  }
  if (fout) {
    fout->Close();
  }
  f->Close();

} // globalFitL3Res


void jesFitter(Int_t& npar, Double_t* grad, Double_t& chi2, Double_t* par,
	       Int_t flag) {

  // Basic checks
  assert(_vdt1);
  assert(_vsrc);
  assert(int(_vsrc->size())==npar-_jesFit->GetNpar()); // nuisance per source
  assert(_vdt1->size()<=sizeof(int)*8); // bitmap large enough

  Double_t *ps = &par[_jesFit->GetNpar()];
  int ns = npar - _jesFit->GetNpar();

  if (flag) {
    
    // do the following calculation:
    // chi2 = sum_i( (x_i+sum_s(a_s y_si) -fit)^2/sigma_i^2) + sum_s(a_s^2)
    chi2 = 0;
    Nk = 0;

    // Loop over input data (graphs x bins)
    // - for each point, add up source eigenvectors x nuisance parameters
    // - then calculate chi2 adding up residuals + nuisance parameters
    for (unsigned int ig = 0; ig != _vdt1->size(); ++ig) {
      TGraphErrors *g = (*_vdt1)[ig]; assert(g);
      for (int i = 0; i != g->GetN(); ++i) {

	// Retrieve central value and uncertainty for this point
	double pt = g->GetX()[i];
	double data = g->GetY()[i];
	double sigma = g->GetEY()[i];

	if (pt < ptminJesFit) continue;

	// Calculate fit value at this point
	for (int ipar = 0; ipar != _jesFit->GetNpar(); ++ipar) {
	  _jesFit->SetParameter(ipar, par[ipar]);
	}
	double fit = _jesFit->EvalPar(&pt,par);

	// For multijet balancing, multiply data by reference JES
	if (TString(g->GetName()).Contains("multijet")) {
	  assert(_vpt1);
	  assert(_vpt1->size()==_nmethods);
	  unsigned int imethod = ig/_nsamples; assert(imethod<_nmethods);
	  TGraphErrors *gf = (*_vpt1)[imethod]; assert(gf);
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
	  int ipt = hsrc->FindBin(pt);

	  if (ison) shifts += ps[is] * hsrc->GetBinContent(ipt);
	}
	
	// Add chi2 from residual
	double chi = (data + shifts - fit) / sigma;
	chi2 += chi * chi;
	++Nk;

	// if _vdt2 is provided, store shifted data there
	{
	  assert(_vdt2);
	  assert(_vdt2->size()==_vdt1->size());

	  TGraphErrors *g2 = (*_vdt2)[ig];
	  if (g2 && g2->GetN()==g->GetN() && (g2->GetX()[i]==pt ||
					      g2->GetX()[i]==pt*shiftPtBal)) {
	    g2->SetPoint(i, pt, data + shifts);

	    // For multijets, store also downward extrapolation
	    if (TString(g->GetName()).Contains("multijet")) {
		assert(_vpt1->size()==_nmethods);
		assert(_vpt2->size()==_nmethods);

	      unsigned int imethod = ig/_nsamples; assert(imethod<_nmethods);
	      TGraphErrors *gf = (*_vpt1)[imethod]; assert(gf);
	      TGraphErrors *gf2 = (*_vpt2)[imethod]; assert(gf2);
	      double jes = _jesFit->EvalPar(&pt, par);
	      double ptref = pt * gf->GetY()[i];
	      double jesref = _jesFit->EvalPar(&ptref, par);
	      // MJB = jes / jesref
	      // data = MJB*jesref => "jesref" = jes / MJB = jes * jesref/data
	      gf2->SetPoint(i, ptref, jes * jesref / (data + shifts));
	      gf2->SetPointError(i, g2->GetEX()[i], g2->GetEY()[i]);
	    }
	  }
	}
      } // for ipt
    } // for ig
  
    // Add "pfjet" compositions (chf, nhf, nef) also into the fit
    if (usePF) {
      for (unsigned int ic = 0; ic != _vpf.size(); ++ic) {
	string sf = _vpf[ic].first;
	TGraphErrors *g = _vpf[ic].second; // input fraction
	TGraphErrors *g2 = _vpf2[ic].second; // output fraction
	for (int i = 0; i != g->GetN(); ++i) {
	  
	  // Retrieve central value and uncertainty for this point
	  double pt = g->GetX()[i];
	  double data = g->GetY()[i];
	  double sigma = g->GetEY()[i];
	  
	  if (pt < ptminJesFit) continue;
	  
	  // Calculate composition bias at this pT	
	  double df(0);
	  for (int ipar = 0; ipar != _jesFit->GetNpar(); ++ipar) {
	    //assert(_mpf[sf][ipar]!=0 || ipar==_jesFit->GetNpar()-1);
	    if (_mpf[sf][ipar]!=0) { // parameter maps to composition change
	      TF1 *f1 = _mpf[sf][ipar]; // effect shape vs pt
	      df += 0.01*f1->Eval(pt)*par[ipar];
	    }
	  } // for ipar
	  // Later: enforce that sum df = 0? (then need to add muf+nef)
	  
	  // Store fitted composition for plotting
	  g2->SetPoint(i, pt, df);

	  // Add chi2 from composition (nuisances not yet considered)
	  //if (usePF) {
	  //if (usePF && pt>28) {
	  if (fitPF && pt > minPFpt && pt < maxPFpt) {
	    //if (usePF && pt>28 && pt<686) {
	    //if (usePF && pt>28 && pt<686 && sf!="nef") {
	  //if (usePF && pt>49 && pt<686) {
	    double chi = (data - df) / sigma;
	    //double chi = 0.1 * (data - df) / sigma; // downscale pfjet weight
	    chi2 += chi * chi;
	    ++Nk;
	  } // usePF
	} // for i
      } // for ic
    } // usePF

    if (usePFZ) {
      for (unsigned int ic = 0; ic != _vpfz.size(); ++ic) {
	string sf = _vpfz[ic].first;
	TGraphErrors *g = _vpfz[ic].second; // input fraction
	TGraphErrors *g2 = _vpfz2[ic].second; // output fraction
	for (int i = 0; i != g->GetN(); ++i) {
	  
	  // Retrieve central value and uncertainty for this point
	  double pt = g->GetX()[i];
	  double data = g->GetY()[i];
	  double sigma = g->GetEY()[i];
	  
	  if (pt < ptminJesFit) continue;
	  
	  // Calculate composition bias at this pT	
	  double df(0);
	  for (int ipar = 0; ipar != _jesFit->GetNpar(); ++ipar) {
	    if (_mpf[sf][ipar]!=0) { // parameter maps to composition change
	      TF1 *f1 = _mpf[sf][ipar]; // effect shape vs pt
	      df += 0.01*f1->Eval(pt)*par[ipar];
	    }
	  } // for ipar
	  
	  // Store fitted composition for plotting
	  g2->SetPoint(i, pt, df);

	  // Add chi2 from composition (nuisances not yet considered)
	  if (fitPFZ && pt > minPFZpt && pt < maxPFZpt) {
	    double chi = (data - df) / sigma;
	    chi2 += chi * chi;
	    ++Nk;
	  } // pt range
	} // for i
      } // for ic
    } // usePFZ

    // Add chi2 from nuisance parameters //including new multijet
    for (int is = 0; is != ns; ++is) {
      
      chi2 += ps[is]*ps[is];
    } // for ipar

    // Penalize fit parameters (outside [0,1] range)
    if (penalizeFitPars>0 && njesFit==9) {
      for (int ipar = 0; ipar != _jesFit->GetNpar(); ++ipar) {
	if (ipar < penalizeFitPars) {
	  chi2 += par[ipar] * par[ipar];
	  //double chi = 0;
	  //if (par[ipar]<0) chi = fabs(par[ipar]);
	  //if (par[ipar]>1) chi = fabs(par[ipar]-1);
	  //chi2 += chi * chi;
	  ++Nk;
	}
      }
    } // penalizeFitPars

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

