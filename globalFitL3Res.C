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

// Global fit settings
double alpha = 0.30; // reference alpha value for zlljet, gamjet (def:0.3)
double ptmj = 30.; // reference pTmin value for multijet (def:30)

// FSR settings
bool dofsr = true;         // correct for FSR central value (def:true)
bool dofsrcombo1 = false;  // use input refIOV for FSR in IOVs (def: true)
bool dofsrcombo2 = false;  // use output refIOV for FSR in IOVs (def:false)
bool dofsrhadw = false; // DP_2021;     // secondary FSR switch for hadW
bool dofsr3 = true;        // use softrad3.C corrections instead of softrad.C
bool dofsr3zj = true;      // secondary softrad3.C switch for Z+jet
bool dofsr3gj = true;      // secondary softrad3.C switch for gamma+jet
bool dofsr3mj = true;      // secondary softrad3.C switch for multijet
bool fixfsrcombo = false;  // do not refit FSR (except refIOV) (def:false)
bool useZJet50 = false;    // Use Z+jet MPF with alpha<0.50 (def:false)
bool useGJet50 = false;    // Use gamma+jet MPF with alpha<0.50 (def:false)
bool alsoPtBal50 = false;  // Use Z/gamma+jet DB with alpha<0.50 (def:false)
bool useZJet100 = true;    // Use Z+jet MPF with alpha<1.00 (def:true)
bool useGJet100 = true;    // Use gamma+jet MPF with alpha<1.00 (def:true)
bool alsoPtBal100 = dofsr3;// Use Z/gamma+jet DB with alpha<1.00 (def:dofsr3)
//bool usePFJet100 = true; // Use PFJet compsition with alpha<1.00
double alphaeffmax = 0.45; // max alphaeff for FSR corrections (def:0.45)

// Further FSR-related settings
double ptreco_gjet = 15.;  // min jet pT when evaluating alphamax (def:15)
double ptreco_zlljet = 5.; // same for Zll+jet (def:5, tbu)
double ptreco_zjet = 15.;  // same for Z+jet (def:15)
bool dol1bias = false;     // correct MPF to L1L2L3-L1 from L1L2L3-RC (def:false)

// Inclusive jet settings (for IOV fits vs reference)
bool useIncJetRefIOVUncert = false;//true; // Include refIOV fit uncertainty? (def:no)

// PF composition settings
bool useSpecialLowPU = false; // DP_2021 // Special 2017H settings if useLowPU
bool usePF = true;     // Load dijet PF fractions in the fit (def:true for now)
bool usePFZ = true;    // Load Z+jet PF fraction in the fit (def:true for now)
bool usePFG = false; // DP_2021    // Load G+jet PF fraction in the fit (def:true for now)
bool usePFC = false;//true;    // Load Combo-PF fraction in the fit (def:true for now)
bool fitPF = true;     // Use dijet PF fractions in the fit (def:true for now)
bool fitPFZ = true;    // Use Z+jet PF fractions in the fit (def:true for now)
bool fitPFG = false; // DP_2021    // Use G+jet PF fractions in the fit (def:true for now)
bool fitPFC = false;//true;    // Use Combo-PF fractions in the fit (def:true for now)
double fitPF_delta = 0.20; // DP_2021 //0.25;  // "free" shift in % for dijet PF (def:0.25)
double fitPFZ_delta = 0.10; // DP_2021 //0.25; // "free" shift in % for Z+jet PF (def:0.25)
double fitPFG_delta = 0.25; // "free" shift in % for Z+jet PF (def:0.25)
double fitPFC_delta = 0.10; // "free" shift in % for Combo-PF (def:0.25)
double minPFpt = 43;   // Minimum pT where dijet PF fractions used (def:114)
double maxPFpt = 600;  // Maximum pT where PF fractions used (def:1000)
double minPFZpt = 25;  // Minimum pT where Z+jet PF fractions used (def:25)
double maxPFZpt = 300; // Maximum pT where Z+jet PF fractios used (def:300)
double minPFGpt = 40;  // Minimum pT where G+jet PF fractions used (def:40)
double maxPFGpt = 600; // Maximum pT where G+jet PF fractios used (def:600)
double minPFCpt = 40;  // Minimum pT where Combo-PF fractions used (def:25)
double maxPFCpt = 600; // Maximum pT where Combo-PF fractios used (def:600)

// Correction level
extern string CorLevel; // defined in reprocess.C ("L1L2Res" or "L1L2L3Res")
bool useL2L3ResForClosure = true;//(fitPF||fitPFZ); // scale out L2L3Res for closure

// Fit parameterization settings
int penalizeFitPars = 9; // Add penalty for 9-par fit parameters (def:9)
bool useP0Trk = true;    // Tracking with CHF change (def:true)
bool useP1Gam = true;    // Photon scale with NEF change (def:true)
bool useP2CalX = false;     // Hadron ECAL+HCAL scale cross (def:false)
bool useP2HadX = true;   // Hadron HCAL scale cross with NHF change (def:true)
bool useP3HadH = true;   // Hadron HCAL scale -3% with NHF change (def:true)
bool useP4HadE = true;   // Hadron ECAL scale -3% with NEF change (def:true)
bool useP5Frag = true;   // Herwig vs Pythia (def:true)
bool useP6L1RC = true;   // Offset bias (def:true)
bool useP7TrkD = true;   // Tracking without CHF change (def:true)
bool useP8MB80pf = false;   // MinBias xsec=80 mb PF fractions (def:false)
bool useP8MB80fit = false;  // MinBias xsec=80 mb (def:false)
bool useP8Flv = true;    // Z+jet / Z+q response shape (def:true)

// Other fit settings
const int njesFit = 9;   // Number of fit parameters (def:9)
double fixOff = 0;
const double ptminJesFit = 15;

// Plotting settings: axis ranges
bool _paper = false;     // Switch plotting to paper style (e.g. no pars)
bool _useZoom = true;    // Zoom y-axis range and change error bands and ref.:
// true--> comes by default with AbsoluteScale+TotalNoFlavorNoTime
// false--> Run1 and reference AbsoluteScale
double shiftPtBal = 0.975;  // Move x-axis for pTbal plots, 0 or 1 for none
double shiftYrange = +0.010;// Move y-axis range for closure plots

// Plottig settings: data points to be shown
bool useIncJet = true;  // Load inclusive jet data, for hjesfit (def:true)
bool plotIncJet = true; // Plot inclusive jet data, if not loaded (def:true)
bool plotHadW = true;   // Plot hadronic W data
bool plotMultijetDown = true; // Plot also (gray) downward points for multijets
double ptmaxMultijetDown = 300; // Max pT for downward multijet points
double _cleanUncert = 0.05; // Max uncertainty for points on the plot (def:0.05)

// Output settings
bool storeErrorMatrix = true; // Error matrix for globalFitL3Pulls.C
bool storeJesFit = true;      // Save jesFit and eigenvectors
bool writeTextFile = true;    // Textfile for minitools/createL2L3ResTextFile.C

// Scaling for inputs: obsolete. Previous L2(L3)Res now stored in reprocess.C
// => See useL2L3ResForClosure above
string scalingForL2OrL3Fit = "None";
// Scaling options listed below are:
// "None", "putBackL2Res", "DontScaleDijets", "ApplyL3ResDontScaleDijets"
// "None" - for input combination files without any residual applied (dijets are still scaled, see below)
// "PutBackL2Res" - put L2res back in for gamma/Z+jet for vs eta studies
// "DontScaleDijets" - don't modify inputs at all (good for full closure tests with L2L3Res applied)  -- not needed in case all the reference JEC bands are moved to 1.0 in reprocess.C 
// "ApplyL3ResDontScaleDijets" - apply barrel JES (use case: check closure when only L2Res is applied to the inputs and L3Res didn't change)
// N.B.: Barrel JES from input text file is always applied to dijet results (unless "ApplyL3ResDontScaleDijets" or "DontScaleDijets" is chosen)

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
vector< pair<string,TGraphErrors*> > _vpfg; // active G+jet PF fractions
vector< pair<string,TGraphErrors*> > _vpfg2; // fitted G+jet PF fractions
vector< pair<string,TGraphErrors*> > _vpfc; // active Combo-PF fractions
vector< pair<string,TGraphErrors*> > _vpfc2; // fitted Combo-PF fractions
map<string, map<int, TF1*> > _mpf; // map of PF fraction variations
// Set PF composition shapes in _mpf, also used by minitools/fullSimShapes.C
void setToyShapeFuncs();
void setFullShapeFuncs();

// Reference for L1L2L3Res closure test
TH1D *_hjesref(0);
// Reference for new JEC (DP2021)
TH1D *_hjesrefnew(0);

void jesFitter(Int_t &npar, Double_t *grad, Double_t &chi2, Double_t *par,
	       Int_t flag);

// Helper functions to draw fit uncertainty band for arbitrary TF1
TF1 *_fitError_func(0);
TMatrixD *_fitError_emat(0);
Double_t fitError(Double_t *xx, Double_t *p);
Double_t fitError1d(Double_t *xx, Double_t *p);
Double_t origJES(Double_t *x, Double_t *);
Double_t newJES(Double_t *x, Double_t *);
void multiplyGraph(TGraphErrors *g, TF1 *f);
						
TF1 *fhb(0), *fl1(0); // Run I functions
TF1 *ft(0), *fc(0), *fcx(0); // Deprecating shapes
TF1 *ftd(0), *ftm(0), *fhx(0), *fhh(0), *feh(0), *fp(0), *fhw(0); // Fit shapes
TF1 *fm80(0), *f1q3(0); // new and experimental
string _epoch ="";
Double_t jesFit(Double_t *x, Double_t *p) {
  
  double pt = *x;

  // Initialize SinglePionHCAL shape (Run I shape)
  if (!fhb) fhb = new TF1("fhb","max(0.,[0]+[1]*pow(x,[2]))",10,3500);
  fhb->SetParameters(1.03091e+00, -5.11540e-02, -1.54227e-01); // SPRH
 
  // Regular +3% SPR calorimeter scale variation (ECAL+HCAL)
  //if (!fc) fc = new TF1("fc","[p0]+[p1]*pow(x/[p2],[p3])/(1+pow(x/[p2],[p3]))*(1-pow(x/[p2],-[p3]))",15,4500);
  //fc->SetParameters(0.335, 0.8171, 406.8, 1.34); // toyPF
  assert(fc);

  // Fits from minitools/varPlots.C
  // SPR -3% to +3% cross variation (ECAL+HCAL)
  //if (!fcx) fcx = new TF1("fcx","[p0]+[p1]*pow(x/[p2],[p3])/(1+pow(x/[p2],[p3]))*(1-pow(x/[p2],-[p3]))",15,4500);
  //fcx->SetParameters(1.177, 1.42, 1388, 1.231); // toyPF
  assert(fcx || !useP2CalX);

  // SPRH -3% to +3% cross variation (HCAL only)
  //if (!fhx) fhx = new TF1("fhx","[p0]+[p1]*pow(x/[p2],[p3])/(1+pow(x/[p2],[p3]))*(1-pow(x/[p2],-[p3]))",15,4500);
  //fhx->SetParameters(0.8904, 1.082, 1408, 1.204); // toyPF
  assert(fhx || !useP2HadX); 

  // SPRH -3% variation
  //if (!fhh) fhh = new TF1("fhh","[p0]+[p1]*pow(x/[p2],[p3])/(1+pow(x/[p2],[p3]))*(1-pow(x/[p2],-[p3]))",15,4500);
  //fhh->SetParameters(-0.7938, -0.5798, 396.1, 1.412); // toyPF
  assert(fhh);

  // SPRE -3% variation
  //if (!feh) feh = new TF1("feh","[p0]+[p1]*pow(x/[p2],[p3])/(1+pow(x/[p2],[p3]))*(1-pow(x/[p2],-[p3]))",15,4500);
  //feh->SetParameters(-0.2603, -0.2196, 409.4, 1.276); // toyPF
  assert(feh);

  // Tracking -3% variation
  //if (!ft) ft = new TF1("ft","[p0]+[p1]*pow(x/208.,[p2])",15,4500);
  //ft->SetParameters(0.057, -0.3845, -0.3051); // toyPF
  assert(ft);

  // Tracking '-1%' variation in UL2017 data
  //if (!ftd) ftd = new TF1("ftd","[p0]+[p1]*pow(x/208.,[p2])+[p3]/x",15,4500);
  //ftd->SetParameters(-0.116, -0.6417, -0.3051, 23.63); // Data
  assert(ftd || !useP7TrkD);

  // Tracking '-1%' variation in UL2017 MC
  //if (!ftm) ftm = new TF1("ftm","[p0]+[p1]*pow(x/208.,[p2])+[p3]/x",15,4500);
  //ftm->SetParameters(0.2683, -0.6994, -0.3051, 18.49); // MC
  assert(ftm || !useP7TrkD);

  // Photon -3% variation
  //if (!fp) fp = new TF1("fp","[p0]",15,4500);
  //fp->SetParameter(0,-0.8295); // toyPF
  assert(fp);

  // sigmaMB 69.2 mb to 80 mb variation
  //if (!fm80) fm80 = new TF1("fm80","[p0]+[p1]*pow(x/208.,[p2])+[p3]*exp(-[p4]*x)",15,4500);
  //fm80->SetParameters(0.03904,-0.01612,-0.942, -7.145,0.09907); // toyPF
  assert(fm80);

  // H++ vs P8CP5
  //if (!fhw) fhw = new TF1("fhw","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]/x+[p5]*log(x)/x",15,4500);
  //fhw->SetParameters(0.9526,-0.3883,1285,2.46,18.1,-2.062); // FullMC
  assert(fhw);

  // Initialize L1FastJet_Simple-L1RC difference
  // Values from fitting ratio/eta00-13/hl1bias (JEC set in reprocess.C)
  if (!fl1) fl1 = new TF1("fl1","1-([0]+[1]*log(x)+[2]*pow(log(x),2))/x",10,3500);
  //fl1->SetParameters(1.72396, 0, 0); // UL17 V1S hl1bias (p1=p2=0) (UL17_V1)
  fl1->SetParameters(0.350077, 0.553560, -0.0527681); // BCDEF hl1rcos (RC vs Simple)

  // Z+q vs Z+jet for |eta|<1.3 from minitools/Zflavor.C through Flavor.C
  // Extrapolated from [45,300] GeV range (1.02@45 to 1.06@15, so quite far)
  if (!f1q3) f1q3 = new TF1("f1q3","1+0.01*([0]+[1]*(pow(0.01*x,[2])-1))",15,3500);
  f1q3->SetParameters(0.7966, 0.9311, -1);

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

    // Don't use both ECAL+HCAL cross and HCAL cross, too similar
    assert(!(useP2HadX && useP2CalX));

    // Scale out L2L3Res so can fit original parameters and use PF composition
    assert(_hjesref || !(useL2L3ResForClosure && CorLevel=="L1L2L3Res"));
    double jesref = (useL2L3ResForClosure && CorLevel=="L1L2L3Res" ?
		     _hjesref->Interpolate(pt) : 1);

    // Calculate parameterization with each component switchable on/off
    return (1
	    + (useP0Trk  ? p[0]*0.01*ftd->Eval(pt) : 0) // Tracker eff. (data) // DP_2021
	    //+ (useP0Trk  ? p[0]*0.01*ft->Eval(pt) : 0) // Tracker eff. (data) // latest
	    + (useP1Gam  ? p[1]*0.01*fp->Eval(pt) : 0)  // photon scale
	    + (useP2CalX ? p[2]*0.01*fcx->Eval(pt) : 0) // Hadron cross-over
	    + (useP2HadX ? p[2]*0.01*fhx->Eval(pt) : 0) // Hadron cross-over
	    + (useP3HadH ? p[3]*0.01*fhh->Eval(pt) : 0) // Hadron in HCAL (SPRH)
	    + (useP4HadE ? p[4]*0.01*feh->Eval(pt) : 0) // Hadron in ECAL (SPRE)
	    + (useP5Frag ? 0.3*p[5]*0.01*fhw->Eval(pt) : 0) // Herwig fragmentation (rough estimate of impact on Z+jet)
	    + (useP6L1RC ? p[6]*(fl1->Eval(pt)-1) : 0) // L1RC-L1Simple diff.
	    + (useP7TrkD ? p[7]*0.01*3*(ftd->Eval(pt)-ftm->Eval(pt)) : 0) // Tracker Data-MC ('3%' vs '1%' for useP0Trk)
	    + (useP8MB80fit ? p[8]*0.01*fm80->Eval(pt) : 0)
	    + (useP8Flv  ? p[8]*(1-f1q3->Eval(pt)) : 0)
	    ) / jesref; // Divide by reference JEC for closure
  } // nJesFit==9

  // Not supported number of parameters, exit code
  assert(0);
  exit(0);
}

void globalFitL3Res(double etamin = 0, double etamax = 1.3,
		    string epoch="",
		    string selectSample="Standard_MJDJ_gam_zee_zmm",
		    string selectMethods="PtBalMPF") {

  if(verboseGF) cout << Form("Running globalFitL3Res"
			     "(etamin=%02.2f,etamax=%02.2f,epoch=%s,"
			     " selectSample=%s, selectMethods=%s",
			     etamin,etamax,epoch.c_str(),
			     selectSample.c_str(),selectMethods.c_str())
		     << endl << flush;

  const char *cep = epoch.c_str();
  _epoch = epoch;

  TDirectory *curdir = gDirectory;
  setTDRStyle();

  bool isUL18 = (epoch=="2018ABCD" || epoch=="2018A" || 
		 epoch=="2018B" || epoch=="2018C" || epoch=="2018D");
  bool isUL17 = (epoch=="2017BCDEF" || epoch=="2017B" || epoch=="2017C" ||
		 epoch=="2017D" || epoch=="2017E" || epoch=="2017F");
  bool isUL16 = (epoch=="2016BCDEF" || epoch=="2016BCD" ||
		 epoch=="2016EF" || epoch=="2016GH");
  bool isAPV = (epoch=="2016BCDEF" || epoch=="2016BCD" ||
		epoch=="2016EF");
  bool isRun2 = (epoch=="Run2Test");
  bool isLowPU = (epoch=="2017H");
  string srefIOV = (isUL18 ? "2018ABCD" :
		    (isUL17 ? "2017BCDEF" :
		     (isUL16 ? "2016BCDEF" :
		      (isLowPU ? "Run2Test" :
		       (isRun2 ? "Run2Test" :
			"2018ABCD")))));
  const char *refIOV = srefIOV.c_str();
  bool isRefIOV = (epoch=="2018ABCD" ||
		   epoch=="2017BCDEF" ||
		   epoch=="2016BCDEF" ||
		   epoch=="Run2Test");
  
  //if (isLowPU) { plotIncJet = usePF = usePFG = usePFC = false; }
  if (isLowPU && useSpecialLowPU) {
    usePFC = false; usePF = true; usePFZ = true; usePFG = false;
    fitPFC = false; fitPF = true; fitPFZ = true; fitPFG = false;
    // Only PF 393.9 / 390
    // Also PFZ 415.3 / 432
    maxPFpt = 2000;
    minPFpt = minPFZpt = 15.;
    fitPF_delta = fitPFZ_delta = 0.25;
    useIncJetRefIOVUncert = true;
  }

  // Input data
  TFile *f = new TFile(Form("rootfiles/jecdata%s.root",cep),"READ");
  assert(f && !f->IsZombie());

  // Input systematics from reference IOV
  TFile *fsys = new TFile(Form("rootfiles/jecdata%s.root",refIOV),"READ");
  assert(fsys && !fsys->IsZombie());

  // FSR corrections from reference IOV, in case fixed
  TFile *ffsr = new TFile(Form("rootfiles/jecdata%s.root",refIOV),"READ");

  // Output JES fit results, if requested or if reference IOV
  TFile *fout = new TFile(Form("rootfiles/jecdata%s.root",cep),
			  storeJesFit || isRefIOV ? "UPDATE" : "READ");
  assert(fout && !fout->IsZombie());
  
  // Figure out which variant of FSR corrections used, if any
  bool fsrcombo = (dofsrcombo1||dofsrcombo2);
  assert((ffsr && !ffsr->IsZombie()) || (!fsrcombo && !isRefIOV));
  if (dofsrcombo1) cout << Form("Use input FSR from %s for all IOVs\n",refIOV);
  if (dofsrcombo2 && !isRefIOV) cout << Form("Use output FSR from %s\n",refIOV);
  if (fixfsrcombo && !isRefIOV) cout << "Do not refit FSR" << endl;
  
  // Inclusive jets input file, for plotting as a reference (also from 'f')
  TFile *fij = (isUL16 ? 
		new TFile("../JERCProtoLab/Summer19UL16/L3Residual_incljet/drawDeltaJEC_16ULJECV5_vsBCDEF.root","READ") :
		isUL18 ?
		new TFile("rootfiles/drawDeltaJEC_18UL_JECv3.root","READ") :
		isRun2 ?
		new TFile("rootfiles/drawDeltaJEC_Run2_vsUL4X.root","READ") :
		isLowPU ? // until new 5X with custom 16APV and 16GH
		new TFile("rootfiles/drawDeltaJEC_Run2_vsUL4X.root","READ") :
		new TFile("rootfiles/drawDeltaJEC_17UL_V2M4res_cp2_all_v4.root","READ"));
  assert(fij && !fij->IsZombie());

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
  for(unsigned int i = 0; i != nmethods; ++i) 
    methodsvec.push_back(methodsmap[selectMethods].at(i).c_str());
  const char** methods = &methodsvec.front();
  
  _nmethods = nmethods; // for multijets in global fit

  // Is this wide-bin L3Res fit (as opposed to narrow-bin L2Res)?
  bool isl3 = (etamin==0 && fabs(etamax-1.3)<0.1);
  
  // Samples mapping. NB: this is too complex code, should generalize
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

  // Global fit without multijets/dijets and with merged Z+jet
  samplesmap["gam_zll_hadw"] =   {"gamjet", "zlljet","hadw"};
  nsample0map["gam_zll_hadw"] = 0;
  igjmap["gam_zll_hadw"] = 0;
  izeemap["gam_zll_hadw"] = -1;
  izmmmap["gam_zll_hadw"] = -1;
  izllmap["gam_zll_hadw"] = 1;
  iwmhmap["gam_zll_hadw"] = 2;

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
  samplesmap["MJDJ_gam_zll"] = {(isl3 ? "multijet":"dijet"),"gamjet","zlljet"};
  nsample0map["MJDJ_gam_zll"] = 1;
  igjmap["MJDJ_gam_zll"] = 0;
  izeemap["MJDJ_gam_zll"] = -1;
  izmmmap["MJDJ_gam_zll"] = -1;
  izllmap["MJDJ_gam_zll"] = 1;

  // Global fit with all samples: multijets/dijets, gamma+jet, merged Z+jet
  samplesmap["MJDJ_gam_zmm"] = {(isl3 ? "multijet":"dijet"),"gamjet","zmmjet"};
  nsample0map["MJDJ_gam_zmm"] = 1;
  igjmap["MJDJ_gam_zmm"] = 0;
  izeemap["MJDJ_gam_zmm"] = -1;
  izmmmap["MJDJ_gam_zmm"] = 1;
  izllmap["MJDJ_gam_zll"] = -1;

  // Global fit with all samples: multijets/dijets, gamma+jet, merged Z+jet
  samplesmap["MJDJ_gam_z"] = { (isl3 ? "multijet" : "dijet"),"gamjet", "zjet"};
  nsample0map["MJDJ_gam_z"] = 1;
  igjmap["MJDJ_gam_z"] = 0;
  izeemap["MJDJ_gam_z"] = -1;
  izmmmap["MJDJ_gam_z"] = -1;
  izllmap["MJDJ_gam_z"] = 1;

  // Global fit with all samples and inclusive jets as extra
  samplesmap["MJDJ_inc_gam_zll"] = {(isl3 ? "multijet" : "dijet"),
				    "incjet", "gamjet", "zlljet"};
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

  // Global fit with all samples and inclusive jets + hadw as extra
  samplesmap["MJDJ_inc_gam_z_hadw"] =   { (isl3 ? "multijet" : "dijet"),"incjet","gamjet", "zjet", "hadw"};
  nsample0map["MJDJ_inc_gam_z_hadw"] = 2;
  igjmap["MJDJ_inc_gam_z_hadw"] = 0;
  izeemap["MJDJ_inc_gam_z_hadw"] = -1;
  izmmmap["MJDJ_inc_gam_z_hadw"] = -1;
  izllmap["MJDJ_inc_gam_z_hadw"] = 1;
  iwmhmap["MJDJ_inc_gam_z_hadw"] = 2;

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
  
  samplesmap["zmm"] =   {"zmmjet"};
  nsample0map["zmm"] = 0;
  igjmap["zmm"] = -1;
  izeemap["zmm"] = -1;
  izmmmap["zmm"] = 0;
  izllmap["zmm"] = -1;

  samplesmap["z"] =   {"zjet"};
  nsample0map["z"] = 0;
  igjmap["z"] = -1;
  izeemap["z"] = -1;
  izmmmap["z"] = -1;
  izllmap["z"] = 0;

  samplesmap["inc_z"] =   {"incjet","zjet"};
  nsample0map["inc_z"] = 1;
  igjmap["inc_z"] = -1;
  izeemap["inc_z"] = -1;
  izmmmap["inc_z"] = -1;
  izllmap["inc_z"] = 1;
  
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
  const int nsample0 = nsample0map[selectSample]; // first Z/gamma+jet sample
  const int imj = 0; // imjmap later
  const int igj = igjmap[selectSample];
  const int izee = izeemap[selectSample];
  const int izmm = izmmmap[selectSample];
  const int izll = izllmap[selectSample];
  const int iwmh = iwmhmap[selectSample];

  _nsamples = nsamples; // for multijets in global fit
  
  // Switch off duplicate plotting of inclusive jets
  if (sampleset.find("incjet")!=sampleset.end()) plotIncJet = false;

  string sbin = Form("eta%02.0f-%02.0f",10*etamin,10*etamax);
  const char* bin = sbin.c_str();

  assert(fsys->cd(ct));
  TDirectory *dsys1 = fsys->GetDirectory(ct); assert(dsys1);
  assert(dsys1->cd(bin));
  TDirectory *dsys = dsys1->GetDirectory(bin); assert(dsys);
  curdir->cd();

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

  if (usePF || usePFZ || usePFG || usePFC) {
    assert(njesFit==9);

    TGraphErrors *gchf(0), *gnhf(0), *gnef(0);
    TGraphErrors *gchfz(0), *gnhfz(0), *gnefz(0);
    TGraphErrors *gchfg(0), *gnhfg(0), *gnefg(0);
    TGraphErrors *gchfc(0), *gnhfc(0), *gnefc(0);
    if (usePFZ) {
      if (useZJet100) {
	gchfz = (TGraphErrors*)d->Get("chf_zlljet_a100");
	gnhfz = (TGraphErrors*)d->Get("nhf_zlljet_a100");
	gnefz = (TGraphErrors*)d->Get("nef_zlljet_a100");
	if (!gchfz && !gnhfz && !gnefz) {
	  gchfz = (TGraphErrors*)d->Get("chf_zmmjet_a100");
	  gnhfz = (TGraphErrors*)d->Get("nhf_zmmjet_a100");
	  gnefz = (TGraphErrors*)d->Get("nef_zmmjet_a100");
	}
	if (sampleset.find("zjet")!=sampleset.end()) {
	  gchfz = (TGraphErrors*)d->Get("chf_zjet_a100");
	  gnhfz = (TGraphErrors*)d->Get("nhf_zjet_a100");
	  gnefz = (TGraphErrors*)d->Get("nef_zjet_a100");
	}
      }
      else {
	gchfz = (TGraphErrors*)d->Get("chf_zlljet_a30");
	gnhfz = (TGraphErrors*)d->Get("nhf_zlljet_a30");
	gnefz = (TGraphErrors*)d->Get("nef_zlljet_a30");
	if (sampleset.find("zjet")!=sampleset.end()) {
	  gchfz = (TGraphErrors*)d->Get("chf_zjet_a30");
	  gnhfz = (TGraphErrors*)d->Get("nhf_zjet_a30");
	  gnefz = (TGraphErrors*)d->Get("nef_zjet_a30");
	}
      }
      assert(gchfz);
      assert(gnhfz);
      assert(gnefz);
      // Use open markers for Z+jet if both present and also fitting dijet
      if (usePF && fitPF) {
	gchfz->SetMarkerStyle(kOpenCircle);
	gnhfz->SetMarkerStyle(kOpenDiamond);
	gnefz->SetMarkerStyle(kOpenSquare);
      }
      //vector<pair<string, TGraphErrors> > _vpfz;
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
      //if (usePFJet100) {
      //gchf = (TGraphErrors*)d->Get("chf_pfjet_a100");
      //gnhf = (TGraphErrors*)d->Get("nhf_pfjet_a100");
      //gnef = (TGraphErrors*)d->Get("nef_pfjet_a100");
      //}
      //else {
      // Only a30 available for now; should add a100?
      gchf = (TGraphErrors*)d->Get("chf_pfjet_a30");
      gnhf = (TGraphErrors*)d->Get("nhf_pfjet_a30");
      gnef = (TGraphErrors*)d->Get("nef_pfjet_a30");
      //}
      assert(gchf);
      assert(gnhf);
      assert(gnef);

      // Use open markers for dijet if both present, but only fitting Z
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
    if (usePFG) {
      gchfg = (TGraphErrors*)d->Get("chf_gamjet_a100");
      gnhfg = (TGraphErrors*)d->Get("nhf_gamjet_a100");
      gnefg = (TGraphErrors*)d->Get("nef_gamjet_a100");

      assert(gchfg);
      assert(gnhfg);
      assert(gnefg);

      //vector<pair<string, TGraphErrors*> > _vpf;
      // std::move needed to turn lvalue to rvalue, as std::make_pair
      // only accepts T2&&, which is rvalue in C++11 standard
      _vpfg.push_back(make_pair<string,TGraphErrors*>("chf", move(gchfg)));
      _vpfg.push_back(make_pair<string,TGraphErrors*>("nhf", move(gnhfg)));
      _vpfg.push_back(make_pair<string,TGraphErrors*>("nef", move(gnefg)));
      
      // Graphs for output fractions
      curdir->cd();
      TGraphErrors *gchfg2 = (TGraphErrors*)gchfg->Clone("gchfg2");
      TGraphErrors *gnhfg2 = (TGraphErrors*)gnhfg->Clone("gnhfg2");
      TGraphErrors *gnefg2 = (TGraphErrors*)gnefg->Clone("gnefg2");
      _vpfg2.push_back(make_pair<string,TGraphErrors*>("chf", move(gchfg2)));
      _vpfg2.push_back(make_pair<string,TGraphErrors*>("nhf", move(gnhfg2)));
      _vpfg2.push_back(make_pair<string,TGraphErrors*>("nef", move(gnefg2)));
    }
    if (usePFC) {
      gchfc = (TGraphErrors*)d->Get("chf_cmb_ren");
      gnhfc = (TGraphErrors*)d->Get("nhf_cmb_ren");
      gnefc = (TGraphErrors*)d->Get("nef_cmb_ren");

      assert(gchfc);
      assert(gnhfc);
      assert(gnefc);

      //vector<pair<string, TGraphErrors*> > _vpf;
      // std::move needed to turn lvalue to rvalue, as std::make_pair
      // only accepts T2&&, which is rvalue in C++11 standard
      _vpfc.push_back(make_pair<string,TGraphErrors*>("chf", move(gchfc)));
      _vpfc.push_back(make_pair<string,TGraphErrors*>("nhf", move(gnhfc)));
      _vpfc.push_back(make_pair<string,TGraphErrors*>("nef", move(gnefc)));
      
      // Graphs for output fractions
      curdir->cd();
      TGraphErrors *gchfc2 = (TGraphErrors*)gchfc->Clone("gchfc2");
      TGraphErrors *gnhfc2 = (TGraphErrors*)gnhfc->Clone("gnhfc2");
      TGraphErrors *gnefc2 = (TGraphErrors*)gnefc->Clone("gnefc2");
      _vpfc2.push_back(make_pair<string,TGraphErrors*>("chf", move(gchfc2)));
      _vpfc2.push_back(make_pair<string,TGraphErrors*>("nhf", move(gnhfc2)));
      _vpfc2.push_back(make_pair<string,TGraphErrors*>("nef", move(gnefc2)));
    }
    
  } // usePF
  setToyShapeFuncs(); // DP_2021
  //setFullShapeFuncs(); // DP_2021 => off


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
		      ss=="multijet" ? ptmj :
		      (ss=="incjet" ? 100 : alpha*100));
      if (useZJet50  && (ss=="zjet" || ss=="zeejet" || ss=="zmmjet"||
			 ss=="zlljet") &&
	  (sm=="mpfchs1" || alsoPtBal50))
	s = Form("%s_%s_a%02.0f",cm,cs,50.);
      if (useGJet50  && ss=="gamjet" && (sm=="mpfchs1" || alsoPtBal50))
	s = Form("%s_%s_a%02.0f",cm,cs,50.);
      if (useZJet100 && (ss=="zjet" || ss=="zeejet" || ss=="zmmjet" ||
			 ss=="zlljet") &&
	  (sm=="mpfchs1" || alsoPtBal100))
	s = Form("%s_%s_a%02.0f",cm,cs,100.);
      if (useGJet100  && ss=="gamjet" && (sm=="mpfchs1" || alsoPtBal100))
	s = Form("%s_%s_a%02.0f",cm,cs,100.);
      TGraphErrors *g = (TGraphErrors*)d->Get(s.c_str());
      if (!g) cout << "Graph "<<s<<" not found!" << endl << flush;
      assert(g);

      // Clone to keep original data in file intact
      g = (TGraphErrors*)g->Clone(Form("g_%s_%s_%s",g->GetName(),cs,cm));
      
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
    
    // Fractions for MPF (pT>15 GeV), only without HDM
    if (!dofsr3 && !dofsr3mj) s = "crecoil_multijet_a15";
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
  TH1D *hjes0 = (TH1D*)f->Get("ratio/eta00-13/hjes"); assert(hjes0);
  hjes0->SetName("hjes0"); // L2Res only (usually)
  TH1D *herr0 = (TH1D*)f->Get("ratio/eta00-13/herr"); assert(herr0);
  herr0->SetName("herr0"); // L2L3Res (usually)

  // Load JEC uncertainty band
  TH1D *herr = (TH1D*)d->Get("herr"); assert(herr);
  TH1D *hrun1 = (TH1D*)d->Get("hrun1"); assert(hrun1);
  TH1D *hjes = (TH1D*)d->Get("hjes"); assert(hjes);
  TH1D *herr_l2l3res = (TH1D*)d->Get("herr_l2l3res"); assert(herr_l2l3res);
  TH1D *herr_ref = (TH1D*)d->Get("herr_ref"); assert(herr_ref);
  TH1D *herr_noflv = (TH1D*)d->Get("herr_noflv"); assert(herr_noflv);
  TH1D *herr_spr = (TH1D*)d->Get("herr_spr"); assert(herr_spr);
  TH1D *herr_pu = (TH1D*)d->Get("herr_pu"); assert(herr_pu);
  
  // Set reference for L1L2L3Res closure
  _hjesref = herr_l2l3res;
  
  // Load inclusive jet reference JES
  TH1D *hjesfitref = (TH1D*)dsys->Get(isRefIOV ? "herr_ref" : "sys/hjesfit");
  assert(hjesfitref);
  
  // Load inclusive jet data
  // => should update this to use jecdata.root instead? 
  TH1D *hij(0);
  TGraphErrors *gij(0);
  if (plotIncJet || useIncJet) {
    map<string,const char*> fij_files;
    fij_files["2018A"] = "A";
    fij_files["2018B"] = "B";
    fij_files["2018C"] = "C";
    fij_files["2018D"] = "D";
    fij_files["2018ABCD"] = "";
    fij_files["2017BCDEF"] = "";
    fij_files["2016BCD"] = "BCD";
    fij_files["2016EF"] = "EF";
    fij_files["2016GH"] = "GH";
    fij_files["2016BCDEF"] = "BCDEF";

    fij_files["Run2Test"] = "UL4X";

    if (isUL18)
      hij = (TH1D*)fij->Get(Form("jet_Run18UL%s_det",fij_files[epoch]));
    else if (isUL17)
      hij = (TH1D*)fij->Get(Form("jet_Run17UL%s_fwd3",fij_files[epoch]));
    else if (isUL16)
      hij = (TH1D*)fij->Get(Form("jet_Run16UL%s_fwd",fij_files[epoch]));
    else if (isRun2)
      hij = (TH1D*)fij->Get(Form("jet_Run%s_fwd",fij_files[epoch]));
    else if (isLowPU)
      hij = (TH1D*)fij->Get(Form("jet_Run%s_fwd","UL17H"));
    else
      assert(false);

    assert(hij);
    for (int i = 1; i != hij->GetNbinsX()+1; ++i) {
      double pt = hij->GetBinCenter(i);
      double jes = hjesfitref->GetBinContent(hjesfitref->FindBin(pt));
      hij->SetBinContent(i, hij->GetBinContent(i) * jes);
    } // for i

    if (plotIncJet) {
      gij = (TGraphErrors*)d->Get("mpfchs1_incjet_a30");
      if (!gij) gij = (TGraphErrors*)d->Get("mpfchs1_incjet_a100");
      assert(gij);
      for (int i = 0; i != gij->GetN(); ++i) {
	double pt = gij->GetX()[i];
	double jes = hjesfitref->GetBinContent(hjesfitref->FindBin(pt));
	gij->SetPoint(i, gij->GetX()[i], gij->GetY()[i] * jes);
      } // for i
    }
  }

  // Load hadronic W data (if being used)
  TH1D *hw(0);
  if (sampleset.find("hadw")!=sampleset.end()) {
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
      // fitProb-based FSR uncertainty for hadronic W, but first zero
      else if (ss=="hadw") {
	//assert(hw);
	//hfsr = (TH1D*)hw->Clone(Form("bm%d_%s",(1<<ibm),hw->GetName()));
	//hfsr->Reset();
	//hks[ibm] = hfsr;

	string s = "fsr/hadw_fsr";
	hfsr = (TH1D*)d->Get(s.c_str()); assert(hfsr);
      }
      else {
	// For correcting FSR bias
	string s = Form("fsr/hkfsr_%s_%s",cm,cs);
	hfsr = (TH1D*)d->Get(s.c_str());
	if (dofsrcombo2 && !isRefIOV) s = Form("fsr/hkfsr2_%s_%s",cm,cs);
	if (fsrcombo && !isRefIOV) hfsr = (TH1D*)dfsr->Get(s.c_str());
	if (dofsr3 && ((dofsr3zj && (ss=="zjet" || ss=="zlljet" ||
				     ss=="zmmjet")) ||
		       (dofsr3mj && ss=="multijet") ||
		       (dofsr3gj && ss=="gamjet"))) {
	  s = Form("fsr/hkfsr3_%s_%s",cm,cs);
	  hfsr = (TH1D*)d->Get(s.c_str());
	  if (fsrcombo && !isRefIOV) hfsr = (TH1D*)dfsr->Get(s.c_str());
	}
	if (!hfsr) cout << "Histo "<<s<<" not found!" << endl << flush;
      }
      assert(hfsr);

      // For correcting L1L2L3-L1 to L1L2L3-RC
      TH1D *hl1 = (TH1D*)d->Get("hl1bias"); assert(hl1);
      
      hfsr->SetName(Form("bm%d_%s_%s",(1<<ibm),cs,hfsr->GetName()));

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
	if (ss=="zeejet" || ss=="zmmjet" || ss=="zlljet")
	  ptreco = ptreco_zlljet;
	if (ss=="zjet") ptreco = ptreco_zjet;
	if (ss=="gamjet") ptreco = ptreco_gjet;
	double aeff = min(max(alpha, ptreco/pt),alphaeffmax);
	if (ss=="multijet") aeff = alpha;
	double kfsr(1);
	if (dofsr) kfsr = 1./(1+aeff*hfsr->GetBinContent(hfsr->FindBin(pt)));
	if (dofsr3 && ((dofsr3zj && (ss=="zjet" || ss=="zlljet" ||
				     ss=="zmmjet")) ||
		       (dofsr3mj && ss=="multijet") ||
		       (dofsr3gj && ss=="gamjet"))) // override regular dofsr
	  kfsr = (r + hfsr->GetBinContent(hfsr->FindBin(pt))) / r;
	if (dofsr && dofsrhadw && ss=="hadw")
	  kfsr = 1./hfsr->GetBinContent(hfsr->FindBin(pt));
	if (ss=="hadw" && !dofsrhadw) kfsr = 1; // DP_2021
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
  map<int, const char*> meig;
  meig[0] = "mpfu1";
  meig[1] = "mpfn1";
  assert(meig.size()>=neigmax);
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
	if (dofsr3 && ((dofsr3zj && (ss=="zjet" || ss=="zlljet" ||
				     ss=="zmmjet")) ||
		       (dofsr3mj && ss=="multijet") ||
		       (dofsr3gj && ss=="gamjet"))) {
	  s = Form("fsr/hkfsr3_%s_%s_%s",cm,cs,meig[ieig]);
	}

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
	h = (TH1D*)h->Clone(Form("bm%d_%s_%s",(1<<ibm),cs,h->GetName()));

	// Forbid refitting if fixfsrcombo and not BCDEF
	if (fixfsrcombo && !isRefIOV)
	  h->Reset();
	
	// Scale dR/dalpha to dR for given alpha setting
	for (int i = 1; i != h->GetNbinsX()+1; ++i) {
	  double pt = h->GetBinCenter(i);
	  double ptreco(0);
	  if (ss=="zeejet" || ss=="zmmjet" || ss=="zlljet")
	    ptreco = ptreco_zlljet;
	  if (ss=="zjet") ptreco = ptreco_zjet;
	  if (ss=="gamjet") ptreco = ptreco_gjet;
	  if (dofsr3 && ((dofsr3zj && (ss=="zjet" || ss=="zlljet" ||
				       ss=="zmmjet" || ss=="zeejet")) ||
			 (dofsr3mj && ss=="multijet") ||
			 (dofsr3gj && ss=="gamjet"))) {
	    h->SetBinContent(i, h->GetBinContent(i));
	  }
	  else {
	    double aeff = min(max(alpha, ptreco/pt),alphaeffmax);
	    if (ss=="multijet") aeff = alpha;
	    h->SetBinContent(i, aeff*h->GetBinContent(i));
	  }
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
    // Use same scale uncertainty source for MPF and pT balance, hence to ibm's
    if (ss=="dijet" && sm=="ptchs") {
      escale = 0.005; // for JER bias and whatnot
      h2->SetName(Form("bm%d_scale_%03.0f_dijet",
		       (1<<(n0-1) | (1<<(n1-1))), escale*10000.));
    }
    if (ss=="gamjet" && sm=="ptchs") {
      escale = 0.005; // UL17-v2 (after 1.0% rescale)
      if (isUL16) escale = 0.0020; // to better match scale from Zee
      h2->SetName(Form("bm%d_scale_%03.0f_gamjet",
		       (1<<(n0+igj) | (1<<(n1+igj))), escale*10000.));
      is = hs.size();
      is_gj = hs.size();
    }
    if (ss=="zeejet" && sm=="ptchs") {
      escale = 0.002; // Consistent with Zmm to 0.2% after mass fit and fix
      if (isUL16) escale = 0.0010;
      h2->SetName(Form("bm%d_scale_%03.0f_zeejet",
		       (1<<(n0+izee) | (1<<(n1+izee))), escale*10000));
      is_zee = hs.size();
    }
    if (ss=="zmmjet" && sm=="ptchs") {
      escale = 0.002; // Z mass fit within 0.2% of unity, flat vs pT
      if (isUL16) escale = 0.0010;
      h2->SetName(Form("bm%d_scale_%03.0f_zmmjet",
		       (1<<(n0+izmm) | (1<<(n1+izmm))), escale*10000));
      is_zmm = hs.size();
    }
    if (ss=="zlljet" && sm=="ptchs") {
      escale = 0.0020; // Zmm mass within 0.2% of unity, Zee fixed to it
      if (isUL16) escale = 0.0010;
      h2->SetName(Form("bm%d_scale_%03.0f_zlljet",
		       (1<<(n0+izll) | (1<<(n1+izll))), escale*10000));
      is_zll = hs.size();
    }
    if (ss=="zjet" && sm=="ptchs") {
      escale = 0.0020; // Zmm mass within 0.2% of unity, Zee fixed to it
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
  // In case of softrad3.C, apply same uncertainty to both DB and MPF
  for (unsigned int i = 0; i != _vdt1->size(); ++i) {

    unsigned int isample = i%nsamples;
    string s = samples[isample];
    const char *cs = s.c_str();
    unsigned int imethod = i/nsamples;
    string m = methods[imethod];
    const char *cm = m.c_str();

    double escale(0); int iref(0);
    const int n0 = nsample0;
    const int n1 = nsamples+nsample0;
    if (s=="gamjet" && m=="mpfchs1") { escale = 0.002; iref = igj;  } 
    if (s=="zeejet" && m=="mpfchs1") { escale = 0.002; iref = izee; }
    if (s=="zmmjet" && m=="mpfchs1") { escale = 0.002; iref = izmm; }
    if (s=="zlljet" && m=="mpfchs1") { escale = 0.002; iref = izll; }
    if (s=="zjet"   && m=="mpfchs1") { escale = 0.002; iref = izll; }

    TH1D *h = hs[i]; assert(h);
    TH1D *h2 = (TH1D*)h->Clone(Form("bm%d_%s_%03.0f_%s_%s",
				    dofsr3 ? 
				    ((1<<(n0+iref)) | (1<<(n1+iref))) :
				    1<<i,
				    escale!=0 ?
				    (dofsr3 ? "hdmscale" : "mpfscale") :
				    "inactive",
				    escale*10000.,cm,cs));
    
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
  if (!fsrcombo && !dofsr3) {
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
  } // !fsrcombo && !dofsr3

  // UL17-V2: add l2cos x <rho>/<rho>_ref as low pT systematic
  if (false) {

    // Derived with help of minitools/rhoBias.C
    TF1 *fl2 = new TF1("fl2","1-([0]+[1]*log(x)+[2]*pow(log(x),2))/x"
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
    const int n0 = 0; // multijet scheme, imj = 0
    const int n1 = nsamples; // multijet scheme, imj = 0
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
			  (1<<(n0+imj) | (1<<(n1+imj))),isrc,int(hs.size()+1)));
			  //((1<<0) | (1<<nsamples)),isrc,int(hs.size()+1)));
	  hs.push_back(h);
	}
      } // for imethod
    } // for ieig
  } // multijet syst.

  // Uncertainty sources for inclusive jets
  if (sampleset.find("incjet")!=sampleset.end() && useIncJetRefIOVUncert) {
    
    const int iij(isLowPU ? 0 : 1);
    const int n0 = 0; // multijet scheme, imj = 0
    const int n1 = nsamples; // multijet scheme, imj = 0
    assert(string(samples[iij])=="incjet");

    // Automatically load unknown number of systematic sources
    TH1D *h(0); int neig(0);
    while ( (h = (TH1D*)dsys->Get(Form("sys/hjesfit_eig%d",neig))) ) {
      h->SetName(Form("bm%d_incjet_jer_src%d_%d",
		      ((1<<(n0+iij)) | (1<<(n1+iij))),neig,
		      int(hs.size()+1)));
      hs.push_back(h);
      ++neig;
    }
  } // incjet syst.

  // Additional pileup uncertainty for inclusive jets
  // NB: low-PU may fit a difference similar to PileUpMuZero in
  // JECUncert_PileUp_AK4PFchs_Eta00.pdf
  if (useSpecialLowPU && sampleset.find("incjet")!=sampleset.end()) {
    
    const int iij(isLowPU ? 0 : 1);
    const int n0 = 0;
    const int n1 = nsamples;
    assert(string(samples[iij])=="incjet");
    TH1D *hpu = (TH1D*)dsys->Get("herr_pu"); assert(hpu);
    TH1D *h = (TH1D*)hpu->Clone(Form("bm%d_incjet_pu_%d",
				     ((1<<(n0+iij)) | (1<<(n1+iij))),
				     int(hs.size()+1)));
    for (int i = 1; i != hpu->GetNbinsX()+1; ++i) {
      double pt = h->GetBinCenter(i);
      //if (pt < 100) h->SetBinContent(i, hpu->GetBinError(i));
      //if (pt >= 100 && pt<600) h->SetBinContent(i, -hpu->GetBinError(i));
      //if (pt >= 600) h->SetBinContent(i, hpu->GetBinError(i));
      //h->SetBinContent(i, hpu->GetBinError(i));
      //if (pt < 100) h->SetBinContent(i, hpu->GetBinError(i));
      if (pt < 80) h->SetBinContent(i, hpu->GetBinError(i));
      else          h->SetBinContent(i, h->GetBinContent(i-1)*0.5);
      //else          h->SetBinContent(i, 0.);
      //assert(false); // check where PU should change sign
    }
    hs.push_back(h);
    // without this source 415.3 / 432
    // with signed 415.3 / 432, fits 0.03+/-0.19
    // with unsigned 399.5 / 432, fits 0.86+/-0.22
    // with pt<100+ 413.5 / 432, fits 0.37+/-0.27
    // with pt<100x0 410.4 / 432, fits 0.54+/-0.24
    // with pt<80x0 402.4 / 432, fits 0.81+/-0.23
    // with pt<80+ 410.0 / 432, fits 0.60+/-0.26
  } // incjet PU

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
  if (sampleset.find("gamjet")!=sampleset.end()) {

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
  h->SetMinimum(h->GetMinimum()+shiftYrange);
  h->SetMaximum(h->GetMaximum()+shiftYrange);
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  h->DrawClone("AXIS");

  map<string, const char*> lumimap;
  lumimap["BCDEF"] = "2017, 41.5 fb^{-1}"; // for DP note
  lumimap["2017BCDEF"] = "2017, 41.5 fb^{-1}"; // for DP note
  lumimap["B"] = "Run2017B, 4.8 fb^{-1}";
  lumimap["C"] = "Run2017C, 9.6 fb^{-1}";
  lumimap["D"] = "Run2017D, 4.2 fb^{-1}";
  lumimap["E"] = "Run2017E, 9.3 fb^{-1}";
  lumimap["F"] = "Run2017F, 13.4 fb^{-1}";
  lumimap["2018A"] = "2018A"; // placeholder
  lumimap["2018B"] = "2018B"; // placeholder
  lumimap["2018C"] = "2018C"; // placeholder
  lumimap["2018D"] = "2018D"; // placeholder
  lumimap["2018ABCD"] = "2018"; // placeholder
  lumimap["2016BCD"] = "Run2016BCD, 12.9 fb^{-1}";
  lumimap["2016EF"] = "Run2016EF, 6.8 fb^{-1}";
  lumimap["2016GH"] = "Run2016GH, 16.8 fb^{-1}";
  lumimap["2016BCDEF"] = "Run2016BCDEF, 19.8 fb^{-1}";
  lumimap["2016BCDEFGH"] = "Run2016BCDEFGH, 36.5 fb^{-1}";
  lumimap["Run2Test"] = "Run 2, 137.6 fb^{-1}";
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
    //legp->SetHeader("p_{T}^{bal}");
    legp->SetHeader("DB"); // ATLAS-style name, better typesetting
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
  //texlabel["zlljet"] = "Z+jet(KIT)";
  texlabel["zjet"] = "Z+jet(UH)";

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
  herr_ref->SetFillColor(kYellow+1);
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
    tdrDraw(gij,"Pz",kFullDiamond,kOrange+2);
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
    if (ss=="zmmjet")   { g2->SetMarkerColor(kRed); g2->SetLineColor(kRed); }
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
    if (!_paper) {
      if (isUL18)
	legm->AddEntry(herr_ref,"UL18 V5","FL");
      else if (isUL17) 
	legm->AddEntry(herr_ref,"UL17 V5","FL");
      else if (isUL16)
	//legm->AddEntry(herr_ref,"UL16 V5M1","FL");
	legm->AddEntry(herr_ref,"UL16 V5","FL");
      else if (isRun2)
	legm->AddEntry(herr_ref,"Run2","FL");
      else if (isLowPU)
	legm->AddEntry(herr_ref,"2017H","FL");
      else
	assert(false);
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
    tdrDraw(gij,"Pz",kFullDiamond,kOrange+2);
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
    if (ss=="zmmjet")   { g3->SetMarkerColor(kRed); g3->SetLineColor(kRed); }
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
  if (njesFit==3)  jesfit->SetParameters(0.98, 1.0, 0.);
  if (njesFit==4)  jesfit->SetParameters(0.9811, 1.0868, -0.1813, 0);
  if (njesFit==5)  jesfit->SetParameters(2., 0.,2., 3., 0.);
  if (njesFit==6)  jesfit->SetParameters(1,1,1,0,0,0);

  //if (njesFit==9)  jesfit->SetParameters(1,1,1,1,1,0,0,0,0);
  if (njesFit==9)  jesfit->SetParameters(1,1,1,1,1,0,0,1,0);
  // Setting originally for L1L2L3Res, obsolete with useL2L3ResForClosure
  if (njesFit==9 && !useL2L3ResForClosure &&
      fabs(herr_ref->GetBinContent(herr_ref->FindBin(100.))-1)<1e-3)
    jesfit->SetParameters(0,0,0,0,0,0,0,0,0);

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
    //if (fabs(tmp_par[i])!=0 || fabs(tmp_err[i]-1)>1e-2) {
    if (!(fabs(tmp_par[i])<1.5e-2 && fabs(tmp_err[i]-1)<1.5e-2)) {
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
    TFile *femat = new TFile(Form("rootfiles/globalFitL3Res_emat_%s.root",
				  epoch.c_str()),"RECREATE");
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

  cout << "Fit parameters:" << endl;
  for (int i = 0; i != np; ++i) {
    cout << Form("%2d: %9.4f +/- %6.4f,   ",
		 i+1,tmp_par[i],sqrt(emat[i][i]));
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
    //if (tmp_par[i]==0&&fabs(tmp_err[i]-1)<2e-2)
    if (fabs(tmp_par[i])<1e-5&&fabs(tmp_err[i]-1)<1.5e-2)
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
    tdrDraw(gij,"Pz",kFullDiamond,kOrange+2);
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
  vector<TGraphErrors*> vgc1ds;
  for (unsigned int i = 0; i != _vdt2->size(); ++i) {
    
    TGraphErrors *g2 = (*_vdt2)[i];
    string ss = samples[i%nsamples];
    string sm = methods[i/nsamples];
    const char *cs = ss.c_str();

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
	  vgc1ds.push_back((TGraphErrors*)gf2tmp->Clone(Form("%sd2",cs)));
	}
      }
    }
    
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
    if (ss=="zmmjet")   { g2->SetMarkerColor(kRed); g2->SetLineColor(kRed); }
    if (ss=="zjet")     { g2->SetMarkerColor(kRed); g2->SetLineColor(kRed); }
    if (ss=="hadw")     { g2->SetMarkerSize(0.5); }

    g2->DrawClone("SAMEPz");
    vgc1ds.push_back((TGraphErrors*)g2->Clone(Form("%s2",cs)));
  } // for i in _vdt2

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
  ((TF1*)fke->DrawClone("SAME"))->SetName("fitError_down");
  fke->SetParameter(0,+1);
  ((TF1*)fke->DrawClone("SAME"))->SetName("fitError_down");

  if (!_paper && njesFit==9) {

    const char *name[] = {"trk","ph","hx","hh","he",
			  "hw","rc","td","flv"};//"mb"};

    tex->DrawLatex(0.54,0.70-0.09,Form("Full 9-par model"));

    double ts_orig = tex->GetTextSize();
    double ts_new = 0.018;
    tex->SetTextSize(ts_new);
    for (int i = 0; i != njesFit; ++i) {
      tex->DrawLatex(0.36,0.70-i*ts_new,
		     Form("p%d(%s)=%+1.3f#pm%1.3f",
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
  TH1D *hjesfit2 = (TH1D*)hij->Clone("hjesfit2");
  hjesfit2->Reset(); // After L2Res (if using L1L2L3Res)
  for (int i = 1; i != hjesfit->GetNbinsX()+1; ++i) {
    double sume2(0);
    for (int ieig = 0; ieig != np; ++ieig) {
      sume2 += pow(vhkeig[ieig]->GetBinContent(i),2);
    }
    double pt = hjesfit->GetBinCenter(i);
    hjesfit->SetBinContent(i, _jesFit->Eval(pt));
    hjesfit->SetBinError(i, sqrt(sume2));
    
    double jesref = (useL2L3ResForClosure && CorLevel=="L1L2L3Res" ?
		     _hjesref->Interpolate(pt) : 1);
    hjesfit2->SetBinContent(i, _jesFit->Eval(pt) * jesref);
    hjesfit2->SetBinError(i, sqrt(sume2));
  } // for i


  ////////////////////////////
  // Draw doubly-shifted data
  ///////////////////////////

  cout << "Draw doubly-shifted data" <<endl;
  
  h = (TH1D*)h->Clone("h2");
  h->SetYTitle("Post-fit jet response (ratio)");  
  h->GetYaxis()->SetRangeUser(0.94+1e-5,1.06-1e-5);

  TCanvas *c1d = tdrCanvas("c1d",h,4,11,true);
  gPad->SetLogx();

  gPad->SetLogx();

  legpf->Draw();
  legmf->Draw();
  legp->Draw();
  legm->Draw();
  legi->Draw();

  l->SetLineStyle(kDashed);
  l->DrawLine(h->GetMinimum(),
	      6500./cosh(etamin),
	      h->GetMinimum()+0.5*(h->GetMaximum()-h->GetMinimum()),
	      6500./cosh(etamin));

  tex->SetTextColor(kBlack);
  tex->SetTextSize(0.045);
  if (etamin==0 && fabs(etamax-1.3)<0.1 && useZJet100)
    tex->DrawLatex(0.20,0.73,"|#eta|<1.3, #alpha<1.0#rightarrow0");
  tex->DrawLatex(0.20,0.17,Form("#chi^{2} / NDF = %1.1f / %d",chi2_gbl,Nk-np));

  // DP2021: error as in herr_ref, expand for herr
  TH1D *herr_dp = (TH1D*)herr_l2l3res->Clone("herr_dp");
  herr_dp->SetLineWidth(2);
  herr_dp->SetLineColor(kYellow+3);
  herr_dp->SetLineStyle(kDashed);
  herr_dp->SetFillColorAlpha(kYellow+1,0.8);
  ((TH1D*)herr_dp->DrawClone("SAME E5"))->SetName("herr_dp_band");
  ((TGraph*)(new TGraph(herr_dp))->DrawClone("SAMEL"))
    ->SetName("herr_dp_line");
  TF1 *forigJES = new TF1("origJES",origJES,minpt,maxpt,0);
  TF1 *fnewJES = new TF1("newJES",newJES,minpt,maxpt,0);
  _hjesrefnew = herr_dp; // DP2021

  TH1D *herr1d = (TH1D*)herr->Clone("herr1d");
  herr1d->Multiply(fnewJES);
  ((TH1D*)herr1d->DrawClone("SAME E5"))->SetName("herr1d_band");
  ((TGraph*)(new TGraph(herr1d))->DrawClone("SAMEL"))
    ->SetName("herr1d_line");

  TH1D *herr_ref1d = (TH1D*)herr_ref->Clone("herr_ref1d");
  herr_ref1d->Multiply(fnewJES);
  herr_ref1d->SetFillColorAlpha(kYellow+1,0.8);
  ((TH1D*)herr_ref1d->DrawClone("SAME E5"))->SetName("herr_ref1d_band");
  ((TGraph*)(new TGraph(herr_ref1d))->DrawClone("SAMEL"))
    ->SetName("herr_ref1d_line");

  //if (epoch=="2016BCDEF") {
  if (plotIncJet) {
    TGraphErrors *gij1d = (TGraphErrors*)gij->Clone("gij1d");
    multiplyGraph(gij1d,fnewJES);
    tdrDraw(gij1d,"Pz",kFullDiamond,kOrange+2);
  }
  //}

  herr1d->SetFillStyle(kNone);
  ((TH1D*)herr1d->DrawClone("SAME E5"))->SetName("herr1d_band2");
  ((TGraph*)(new TGraph(herr1d))->DrawClone("SAMEL"))
    ->SetName("herr1d_line2");
  herr1d->SetFillStyle(1001);

  herr_ref1d->SetFillStyle(kNone);
  ((TH1D*)herr_ref1d->DrawClone("SAME E5"))->SetName("herr_ref1d_band2");
  ((TGraph*)(new TGraph(herr_ref1d))->DrawClone("SAMEL"))
    ->SetName("herr_ref1d_line2");
  herr_ref1d->SetFillStyle(1001);

  //_jesFit->SetNpx(1000);
  //_jesFit->DrawClone("SAME");

  for (int i = 0; i != vgc1ds.size(); ++i) {
    multiplyGraph(vgc1ds[i],forigJES);
    vgc1ds[i]->Draw("SAMEPz");
  }

  TF1 *fke1d = new TF1("fitError1d", fitError1d, minpt, maxpt, 1);
  _fitError_func = jesfit;
  _fitError_emat = &emat2;
  //_origJES = new TF1("_origJES","[0]",minpt,maxpt);
  //_origJES->SetParameter(0,1);

  fke1d->SetNpx(1000);
  fke1d->SetLineStyle(kDotted);
  fke1d->SetLineColor(jesfit->GetLineColor());
  fke1d->SetParameter(0,-1);
  ((TF1*)fke1d->DrawClone("SAME"))->SetName("fitError1d_down");
  fke1d->SetParameter(0,+1);
  ((TF1*)fke1d->DrawClone("SAME"))->SetName("fitError1d_up");
  fke1d->SetLineStyle(kSolid);
  fke1d->SetParameter(0,0);
  ((TF1*)fke1d->DrawClone("SAME"))->SetName("fitError1d_mid");

  //////////////////////////////////
  // End drawing doubly-shifted data
  //////////////////////////////////


  //TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  c0->cd();  l->DrawLine(minpt,1,maxpt,1);
  c0b->cd(); l->DrawLine(minpt,1,maxpt,1);
  c1->cd();  l->DrawLine(minpt,1,maxpt,1);
  c1d->cd(); l->DrawLine(minpt,1,maxpt,1);

  c0->RedrawAxis();
  c0b->RedrawAxis();
  c1->RedrawAxis();
  c1d->RedrawAxis();

  if (etamin==0 && (fabs(etamax-1.3)<0.1 || fabs(etamax-2.4)<0.1)) {
    if (dofsr) {
      c0b->SaveAs(Form("pdf/%s/globalFitL3res_raw.pdf",cep));
      c0->SaveAs(Form("pdf/%s/globalFitL3res_orig.pdf",cep));
      c1->SaveAs(Form("pdf/%s/globalFitL3res_shifted.pdf",cep));
      c1->SaveAs(Form("pdf/%s/globalFitL3res_shifted.root",cep));
      c1d->SaveAs(Form("pdf/%s/globalFitL3res_doublyshifted.pdf",cep));
      c1d->SaveAs(Form("pdf/%s/globalFitL3res_doublyshifted.root",cep));
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
      if (ss=="zeejet" || ss=="zmmjet" || ss=="zlljet")
	{  minx=15; maxx = 700; }
      if (ss=="zjet") {  minx=15; maxx = 700; }
      if (ss=="gamjet") { minx=230; maxx = 1500; }
      if (ss=="multijet") { minx=114; maxx = (isUL18 ? 2640 : 2116); }//2640; }
      if (ss=="incjet") { minx=15; maxx = 2116; }//2640; }
      if (ss=="hadw") { minx=30; maxx = 200; }
      hk->GetXaxis()->SetRangeUser(minx,maxx);

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
	if (ss=="zeejet" || ss=="zmmjet" || ss=="zlljet")
	  ptreco = ptreco_zlljet;
	if (ss=="zjet") ptreco = ptreco_zjet;
	if (ss=="gamjet") ptreco = ptreco_gjet;
	if (dofsr3 && ((dofsr3zj && (ss=="zjet" || ss=="zlljet" ||
				     ss=="zmmjet")) ||
		       (dofsr3mj && ss=="multijet") ||
		       (dofsr3gj && ss=="gamjet"))) {
	  hke->SetBinContent(ipt, hk->GetBinContent(ipt) - yeig);
	  hke->SetBinError(ipt, sqrt(err2));
	}
	else {
	  double aeff = min(max(alpha, ptreco/pt),alphaeffmax);
	  if (ss=="multijet") aeff = alpha;
	  hke->SetBinContent(ipt, (aeff*hk->GetBinContent(ipt) - yeig)/aeff);
	  hke->SetBinError(ipt, sqrt(err2)/aeff);
	}
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
      } // isRefIOV
    } // for isample

    // Store also jesFit and it's eigenvectors for refIOV
    // (using refIOV also as input; update code later to also store output IOVs)
    if (storeJesFit || isRefIOV) {
      dout->cd("sys");
      hjesfit->Write("hjesfit",TObject::kOverwrite);
      hjesfit2->Write("hjesfit2",TObject::kOverwrite);
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
		     -2.0+1e-5,2.5-1e-5,"p_{T} (GeV)",ptmin,ptmax); // Run2Test
		     //-1.2,1.3,"p_{T} (GeV)",ptmin,ptmax);
  if (isUL18) {
    h4->SetMinimum(-2.4);
    h4->SetMaximum(+2.6);
  }
  if (isUL16) {
    h4->SetMinimum(-4.0);
    h4->SetMaximum(+4.4);
  }
  TCanvas *c4 = tdrCanvas("c4",h4,4,11,kSquare);
  gPad->SetLogx();

  TLegend *leg4 = tdrLeg(0.56,0.72,0.86,0.90);
  //leg4->SetHeader(fitPF_delta==0? "Dijet" : 
  //		  Form("Dijet (#Delta=%1.2f%%)",fitPF_delta));
  if (usePF) leg4->SetHeader("Dijet");

  TLegend *leg4z = tdrLeg(0.38,0.72,0.68,0.90);
  if (usePFZ) leg4z->SetHeader("Z+jet");
  //leg4z->SetHeader(fitPFZ_delta==0? "Z+jet" : 
  //		   Form("Z+jet (#Delta=%1.2f%%)",fitPFZ_delta));

  TLegend *leg4g = tdrLeg(0.74,0.72,1.04,0.90);
  if (usePFG) leg4g->SetHeader("#gamma+jet");
  //leg4z->SetHeader(fitPFG_delta==0? "#gamma+jet" : 
  //		   Form("#gamma+jet (#Delta=%1.2f%%)",fitPFG_delta));

  TLegend *leg4c = tdrLeg(0.38,0.72,0.68,0.90);
  if (usePFC) leg4c->SetHeader("Combined PF");

  if (usePF||usePFZ||usePFG)
    tex->DrawLatex(0.18,0.17,Form("#Delta=%1.2f%%(z), %1.2f%%(d),"
				  " %1.2f%%(#gamma)",
				  fitPFZ_delta, fitPF_delta, fitPFG_delta));
  if (usePFC)
    tex->DrawLatex(0.18,0.17,Form("#Delta=%1.2f%%(cmb)",fitPFC_delta));
				

  TLine *l4 = new TLine(); l4->SetLineStyle(kDashed);
  l4->DrawLine(ptmin,0,ptmax,0);

  // Draw output fits with bands firsts
  bool isrefz = (usePFZ && !usePF);
  vector< pair<string,TGraphErrors*> > &_vpfx = (isrefz ? _vpfz : 
						 (usePFC ? _vpfc : _vpf));
  vector< pair<string,TGraphErrors*> > &_vpfx2 = (isrefz ? _vpfz2 :
						  (usePFC ? _vpfc2 : _vpf2));
  for (unsigned int ic = 0; ic != _vpfx.size(); ++ic) {
    string sf = _vpfx[ic].first;
    const char *cf = sf.c_str();
    TGraphErrors *g2 = _vpfx2[ic].second;
    TGraphErrors *gc2 = (TGraphErrors*)g2->Clone(Form("gc2_%d",ic));
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
      } // njesFit==9
    } // for i
    // Remove points outside fit range
    TGraphErrors *gc4 = (TGraphErrors*)gc3->Clone(Form("gc4_%d",ic));
    for (int i = gc4->GetN()-1; i != -1; --i) {
      double pt = gc4->GetX()[i];
      if ((!fitPF && !fitPFZ && !fitPFG && !fitPFC) ||
	  (fitPF && !fitPFZ && (pt<minPFpt || pt>maxPFpt)) ||
	  (fitPFZ && !fitPF && (pt<minPFZpt || pt>maxPFZpt)) ||
	  (fitPFZ && fitPF && (pt<min(minPFZpt,minPFpt) ||
			       pt>max(maxPFZpt,maxPFpt))) ||
	  (fitPFC && !fitPF && !fitPFZ && !fitPFG &&
	   (pt<minPFCpt || pt>maxPFCpt)))
	gc4->RemovePoint(i);
    } // for i

    tdrDraw(gc2,"E3",g2->GetMarkerStyle(),g2->GetMarkerColor()+1,kSolid,-1,
	    1001,g2->GetMarkerColor()-9);
    gc2->SetFillColorAlpha(g2->GetMarkerColor()-9,0.2);
    tdrDraw(gc3,"E3",g2->GetMarkerStyle(),g2->GetMarkerColor()+1,kSolid,-1,
    	    1001,g2->GetMarkerColor()-9);
    gc3->SetFillColorAlpha(g2->GetMarkerColor()-9,0.4);
    tdrDraw(gc4,"E3",g2->GetMarkerStyle(),g2->GetMarkerColor()+1,kSolid,-1,
	    1001,g2->GetMarkerColor());
    gc4->SetFillColorAlpha(g2->GetMarkerColor(),0.7);
    tdrDraw((TGraphErrors*)gc2->Clone(Form("gc5_%d",ic)),"LX",
	    kNone,g2->GetMarkerColor(),kSolid,-1,kNone);
  } // for ic

  // Draw input graphs with points on top
  for (unsigned int ic = 0; ic != _vpfz.size(); ++ic) {
    string sf = _vpfz[ic].first;
    const char *cf = sf.c_str();
    TGraphErrors *g = _vpfz[ic].second;
    TGraphErrors *gc = (TGraphErrors*)g->Clone(Form("gc6_%d",ic));
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
    TGraphErrors *gc = (TGraphErrors*)g->Clone(Form("gc7_%d",ic));
    for (int i = 0; i != gc->GetN(); ++i) {
      gc->SetPoint(i, g->GetX()[i], 100.*g->GetY()[i]);
      gc->SetPointError(i, g->GetEX()[i], 100.*g->GetEY()[i]);
    }
    tdrDraw(gc,"Pz",g->GetMarkerStyle(),g->GetMarkerColor());
    gc->SetMarkerSize(0.8);
    leg4->AddEntry(gc,cf,"PLE");
  }
  for (unsigned int ic = 0; ic != _vpfg.size(); ++ic) {
    string sf = _vpfg[ic].first;
    const char *cf = sf.c_str();
    TGraphErrors *g = _vpfg[ic].second;
    TGraphErrors *gc = (TGraphErrors*)g->Clone(Form("gc8_%d",ic));
    for (int i = 0; i != gc->GetN(); ++i) {
      gc->SetPoint(i, g->GetX()[i], 100.*g->GetY()[i]);
      gc->SetPointError(i, g->GetEX()[i], 100.*g->GetEY()[i]);
    }
    tdrDraw(gc,"Pz",g->GetMarkerStyle(),g->GetMarkerColor()+2);
    gc->SetMarkerSize(0.6);
    leg4g->AddEntry(gc,cf,"PLE");
  }
  for (unsigned int ic = 0; ic != _vpfc.size(); ++ic) {
    string sf = _vpfc[ic].first;
    const char *cf = sf.c_str();
    TGraphErrors *g = _vpfc[ic].second;
    TGraphErrors *gc = (TGraphErrors*)g->Clone(Form("gc9_%d",ic));
    for (int i = 0; i != gc->GetN(); ++i) {
      gc->SetPoint(i, g->GetX()[i], 100.*g->GetY()[i]);
      gc->SetPointError(i, g->GetEX()[i], 100.*g->GetEY()[i]);
    }
    tdrDraw(gc,"Pz",g->GetMarkerStyle(),g->GetMarkerColor());
    gc->SetMarkerSize(0.8);
    leg4c->AddEntry(gc,cf,"PLE");
  }

  if (isl3) {
    c4->SaveAs(Form("pdf/%s/globalFitL3res_pfjet.pdf",cep));
    c4->SaveAs(Form("pdf/%s/globalFitL3res_pfjet.root",cep));
  }
  else
    c4->SaveAs(Form("pdf/%s/globalFitL3res_pfjet_eta%02.0f-%02.0f.pdf",
		    cep,10*etamin,10*etamax));

  curdir->cd();
  if (ffsr) {
    ffsr->Close();
  }
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
  
    // Add dijet PF compositions (chf, nhf, nef) into the chi2 fit
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
	    if (_mpf[sf][ipar]!=0) { // parameter maps to composition change
	      TF1 *f1 = _mpf[sf][ipar]; // effect shape vs pt
	      df += 0.01*f1->Eval(pt)*par[ipar];
	    }
	  } // for ipar
	  // Later: enforce that sum df = 0? (then need to add muf+nef)
	  
	  // Store fitted composition for plotting
	  g2->SetPoint(i, pt, df);

	  // Add chi2 from composition (nuisances not yet considered)
	  if (fitPF && pt > minPFpt && pt < maxPFpt) {
	    double chi = max(fabs(data - df) - 0.01*fitPF_delta,0.) / sigma;
	    chi2 += chi * chi;
	    ++Nk;
	  } // usePF
	} // for i
      } // for ic
    } // usePF

    // Add Z+jet PF compositions (chf, nhf, nef) into the chi2 fit
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
	  // Later: enforce that sum df = 0? (then need to add muf+nef)
	  
	  // Store fitted composition for plotting
	  g2->SetPoint(i, pt, df);

	  // Add chi2 from composition (nuisances not yet considered)
	  if (fitPFZ && pt > minPFZpt && pt < maxPFZpt) {
	    //double chi = (data - df) / sigma;
	    // Later: add -0.01*fitPF_delta?
	    double chi = max(fabs(data - df) - 0.01*fitPFZ_delta,0.) / sigma;
	    chi2 += chi * chi;
	    ++Nk;
	  } // pt range
	} // for i
      } // for ic
    } // usePFZ

    // Add G+jet PF compositions (chf, nhf, nef) into the chi2 fit
    if (usePFG) {
      for (unsigned int ic = 0; ic != _vpfg.size(); ++ic) {
	string sf = _vpfg[ic].first;
	TGraphErrors *g = _vpfg[ic].second; // input fraction
	TGraphErrors *g2 = _vpfg2[ic].second; // output fraction
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
	  // Later: enforce that sum df = 0? (then need to add muf+nef)
	  
	  // Store fitted composition for plotting
	  g2->SetPoint(i, pt, df);

	  // Add chi2 from composition (nuisances not yet considered)
	  if (fitPFG && pt > minPFGpt && pt < maxPFGpt) {
	    //double chi = (data - df) / sigma;
	    // Later: add -0.01*fitPF_delta?
	    double chi = max(fabs(data - df) - 0.01*fitPFG_delta,0.) / sigma;
	    chi2 += chi * chi;
	    ++Nk;
	  } // pt range
	} // for i
      } // for ic
    } // usePFG

    // Add G+jet PF compositions (chf, nhf, nef) into the chi2 fit
    if (usePFC) {
      for (unsigned int ic = 0; ic != _vpfc.size(); ++ic) {
	string sf = _vpfc[ic].first;
	TGraphErrors *g = _vpfc[ic].second; // input fraction
	TGraphErrors *g2 = _vpfc2[ic].second; // output fraction
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
	  // Later: enforce that sum df = 0? (then need to add muf+nef)
	  
	  // Store fitted composition for plotting
	  g2->SetPoint(i, pt, df);

	  // Add chi2 from composition (nuisances not yet considered)
	  if (fitPFC && pt > minPFCpt && pt < maxPFCpt) {
	    //double chi = (data - df) / sigma;
	    // Later: add -0.01*fitPF_delta?
	    double chi = max(fabs(data - df) - 0.01*fitPFC_delta,0.) / sigma;
	    chi2 += chi * chi;
	    ++Nk;
	  } // pt range
	} // for i
      } // for ic
    } // usePFC

    // Add chi2 from nuisance parameters //including new multijet
    for (int is = 0; is != ns; ++is) {
      
      chi2 += ps[is]*ps[is];
    } // for ipar

    // Add penalty for first N fit parameters (Bayesian prior, essentially)
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

// DP2021
Double_t fitError1d(Double_t *xx, Double_t *pp) {
  return (fitError(xx,pp) * origJES(xx,pp));
}
Double_t origJES(Double_t *x, Double_t *p) {
  return (_hjesref->Interpolate(*x));
}
Double_t newJES(Double_t *x, Double_t *p) {
  return (_hjesrefnew->Interpolate(*x));
}
void multiplyGraph(TGraphErrors *g, TF1 *f) {

  for (int i = 0; i != g->GetN(); ++i) {
    g->SetPoint(i, g->GetX()[i], g->GetY()[i]*f->Eval(g->GetX()[i]));
    g->SetPointError(i, g->GetEX()[i], g->GetEY()[i]*f->Eval(g->GetX()[i]));
  }
}


void setToyShapeFuncs() {

  // Regular +3% SPR calorimeter scale variation (ECAL+HCAL)
  if (!fc) fc = new TF1("fc","[p0]+[p1]*pow(x/[p2],[p3])/(1+pow(x/[p2],[p3]))*(1-pow(x/[p2],-[p3]))",15,4500);
  fc->SetParameters(0.335, 0.8171, 406.8, 1.34); // toyPF

  // Fits from minitools/varPlots.C
  // SPR -3% to +3% cross variation (ECAL+HCAL)
  if (!fcx) fcx = new TF1("fcx","[p0]+[p1]*pow(x/[p2],[p3])/(1+pow(x/[p2],[p3]))*(1-pow(x/[p2],-[p3]))",15,4500);
  fcx->SetParameters(1.177, 1.42, 1388, 1.231); // toyPF

  // SPRH -3% to +3% cross variation (HCAL only)
  if (!fhx) fhx = new TF1("fhx","[p0]+[p1]*pow(x/[p2],[p3])/(1+pow(x/[p2],[p3]))*(1-pow(x/[p2],-[p3]))",15,4500); // DP_2021
  fhx->SetParameters(0.8904, 1.082, 1408, 1.204); // toyPF // DP_2021
  //if (!fhx) fhx = new TF1("fhc3_Rjet","[p0]+log(x)*([p1]+log(x)*([p2]+log(x)*([p3]+log(x)*([p4]+log(x)*([p5]+log(x)*[p6])))))",15,4500); // latest
  //fhx->SetParameters(-6.575,3.888,-0.4755,-0.06304,0.0063,0.001629,-0.0001509); // 0.0/16 // latest

  // SPRH -3% variation
  if (!fhh) fhh = new TF1("fhh","[p0]+[p1]*pow(x/[p2],[p3])/(1+pow(x/[p2],[p3]))*(1-pow(x/[p2],-[p3]))",15,4500);
  fhh->SetParameters(-0.7938, -0.5798, 396.1, 1.412); // toyPF

  // SPRE -3% variation
  if (!feh) feh = new TF1("feh","[p0]+[p1]*pow(x/[p2],[p3])/(1+pow(x/[p2],[p3]))*(1-pow(x/[p2],-[p3]))",15,4500);
  feh->SetParameters(-0.2603, -0.2196, 409.4, 1.276); // toyPF

  // Tracking -3% variation
  //if (!ft) ft = new TF1("ft","[p0]+[p1]*pow(x/208.,[p2])",15,4500);
  //ft->SetParameters(0.057, -0.3845, -0.3051); // toyPF

  // Tracking '-1%' variation in UL2017 data
  if (!ftd) ftd = new TF1("ftd","[p0]+[p1]*pow(x/208.,[p2])+[p3]/x",15,4500);
  ftd->SetParameters(-0.116, -0.6417, -0.3051, 23.63); // Data
  ft = ftd;

  // Tracking '-1%' variation in UL2017 MC
  if (!ftm) ftm = new TF1("ftm","[p0]+[p1]*pow(x/208.,[p2])+[p3]/x",15,4500);
  ftm->SetParameters(0.2683, -0.6994, -0.3051, 18.49); // MC

  // Photon -3% variation
  if (!fp) fp = new TF1("fp","[p0]",15,4500);
  fp->SetParameter(0,-0.8295); // toyPF

  // sigmaMB 69.2 mb to 80 mb variation
  if (!fm80) fm80 = new TF1("fm80","[p0]+[p1]*pow(x/208.,[p2])+[p3]*exp(-[p4]*x)",15,4500);
  fm80->SetParameters(0.03904,-0.01612,-0.942, -7.145,0.09907); // toyPF

  // H++ vs P8CP5
  if (!fhw) fhw = new TF1("fhw","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]/x+[p5]*log(x)/x",15,4500);
  fhw->SetParameters(0.9526,-0.3883,1285,2.46,18.1,-2.062); // FullMC

    // Fits from minitools/varPlots.C
    //////////////////////////////////

    // Neutral hadron fraction (NHF)
    // Fits from minitools/varPlots.C
    //TF1 *ft3_nhf = new TF1("ft3_nhf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,-0.3)",15,4500);
    //ft3_nhf->SetParameters(1.312, -0.2902, 951.2, 1.289, -0.6573); // 50.9/18
    TF1 *fp_nhf = new TF1("fp_nhf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500);
    fp_nhf->SetParameters(0.07395, 0, 1000, 258.2, 1.223e-05, 1.158); // 7.6/18
    TF1 *fcx_nhf = new TF1("fcx_nhf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500);
    fcx_nhf->SetParameters(-1.54, -0.4506, 156.3, 0.5089, 1.234, 0.09256); // 3.6/17
    TF1 *fhx_nhf = new TF1("fhx_nhf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500); // DP_2021
    fhx_nhf->SetParameters(-0.295, 0.09444, 2713, 2.292, 0.06437, 0.2845); // 3.6/17 // DP_2021
    //TF1 *fhx_nhf = new TF1("fhc3_nhf","[p0]+log(x)*([p1]+log(x)*([p2]+log(x)*([p3]+log(x)*([p4]+log(x)*([p5]+log(x)*[p6])))))",15,4500);
    //fhx_nhf->SetParameters(-2.475,1.441,-0.177,-0.02317,0.002375,0.0006077,-5.567e-05); // 0.0/16

    TF1 *fhm_nhf = new TF1("fhm_nhf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500);
    fhm_nhf->SetParameters(-0.2746, -0.6358, 9664, 0.6547, 0.05559, 0.1816); // 10.1/17
    TF1 *fem_nhf = new TF1("fem_nhf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500);
    fem_nhf->SetParameters(-0.03458, 0, 1713, 274.8, 0.01665, 0.2426); // 3.4/18
    TF1 *ftmg_nhf = new TF1("ftmg_nhf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,+0.3)+[p5]/x",15,4500);
    ftmg_nhf->SetParameters(-0.01022, -0.1962, 4000, 3.071, 0.04211, 0.01005); // 253.1/116
    TF1 *fm80_chf = new TF1("fm80_chf","[p0]+[p1]*pow(x,[p2])",15,4500);
    fm80_chf->SetParameters(1.848, -1.685, -0.006643); // 70.3/53
    TF1 *fhw_nhf = new TF1("fhw_nhf","[p0]+[p1]*pow(x,[p2])+[p3]/x",15,4500);
    fhw_nhf->SetParameters(-5.151, 4.495, 0.03335, -12.3); // 65.4/42

    //<TCanvas::Print>: pdf file pdf/varPlotsComp_nhf.pdf has been created

    // Photon fraction (NEF)
    // Fits from minitools/varPlots.C
    //TF1 *ft3_nef = new TF1("ft3_nef","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,-0.3)",15,4500);
    //ft3_nef->SetParameters(-0.8857, 1.481, 2223, -0.07145, -1.411); // 35.2/18
    TF1 *fp_nef = new TF1("fp_nef","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500);
    fp_nef->SetParameters(2.283, 0, 1000, 1.302, -2.738, 0.002452); // 21.3/18
    TF1 *fcx_nef = new TF1("fcx_nef","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500);
    fcx_nef->SetParameters(0.01659, -0.005997, 754.8, 121.2, -0.0001236, 0.9665); // 7.9/17

    TF1 *fhx_nef = new TF1("fhx_nef","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500); // DP_2021
    fhx_nef->SetParameters(0.05474, -0.003141, 798.6, 78.84, -0.000957, 0.7676); // 3.7/17 // DP_2021
    //TF1 *fhx_nef = new TF1("fhc3_nef","[p0]+log(x)*([p1]+log(x)*([p2]+log(x)*([p3]+log(x)*([p4]+log(x)*([p5]+log(x)*[p6])))))",15,4500);
    //fhx_nef->SetParameters(-1.027,0.6699,-0.09337,-0.01124,0.001406,0.0003225,-3.572e-05); // 0.0/16

    TF1 *fhm_nef = new TF1("fhm_nef","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500);
    fhm_nef->SetParameters(0.4158, -2.14, 9426, 0.1723, 0.4111, 0.1937); // 9.5/17
    TF1 *fem_nef = new TF1("fem_nef","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500);
    fem_nef->SetParameters(-0.02364, 0, 1481, 246.2, -0.009737, 0.2576); // 5.3/18
    TF1 *ftmg_nef = new TF1("ftmg_nef","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])+[p6]/x",15,4500);
    ftmg_nef->SetParameters(0.07453, 0.1457, 1131, -3.68, -0.4155, -0.3, -1.878); // 149.7/115
    TF1 *fm80_nef = new TF1("fm80_nef","[p0]+[p1]*pow(x,[p2])",15,4500);
    fm80_nef->SetParameters(3.611, -3.909, -0.007482); // 70.2/53
    TF1 *fhw_nef = new TF1("fhw_nef","[p0]+[p1]*pow(x,[p2])+[p3]/x",15,4500);
    fhw_nef->SetParameters(0.8417, -0.2605, 0.2289, 2.426); // 35.9/42

    //<TCanvas::Print>: pdf file pdf/varPlotsComp_nef.pdf has been created

    // Charged hadron fraction (CHF)
    // Fits from minitools/varPlots.C
    //TF1 *ft3_chf = new TF1("ft3_chf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,-0.3)",15,4500);
    //ft3_chf->SetParameters(-1.792, 0.2976, 1107, 1.559, 1.039); // 64.1/18
    TF1 *fp_chf = new TF1("fp_chf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500);
    fp_chf->SetParameters(0.3333, 0.7433, 1023, 0.3926, -0.09446, 0.2883); // 10.3/17
    TF1 *fcx_chf = new TF1("fcx_chf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500);
    fcx_chf->SetParameters(-0.1188, -0.3705, 408.1, -0.2583, 1.39, -0.1831); // 2.7/17

    TF1 *fhx_chf = new TF1("fhx_chf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500); // DP_2021
    fhx_chf->SetParameters(-0.0637, -0.2811, 4531, -0.3172, 1.071, -0.153); // 1.7/17 // DP_2021
    //TF1 *fhx_chf = new TF1("fhc3_chf","[p0]+log(x)*([p1]+log(x)*([p2]+log(x)*([p3]+log(x)*([p4]+log(x)*([p5]+log(x)*[p6])))))",15,4500);
    //fhx_chf->SetParameters(3.42,-2.059,0.2634,0.03367,-0.003686,-0.0009108,8.937e-05); // 0.0/16

    TF1 *fhm_chf = new TF1("fhm_chf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500);
    fhm_chf->SetParameters(0.1552, -0.04221, 315.4, 2.787, -0.06628, -0.2572); // 6.2/17
    TF1 *fem_chf = new TF1("fem_chf","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])",15,4500);
    fem_chf->SetParameters(0.06085, 0, 1000, 1.3, -0.008137, 0.2135); // 2.7/18
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
    _mpf["chf"][0] = (useP0Trk  ? ftmg_chf : 0);// tracking in data+MC
    _mpf["chf"][1] = (useP1Gam  ? fp_chf : 0);  // photon
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
    _mpf["nhf"][0] = (useP0Trk  ? ftmg_nhf : 0);// tracking in data+MC
    _mpf["nhf"][1] = (useP1Gam  ? fp_nhf : 0);  // photon
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
    _mpf["nef"][1] = (useP1Gam ? fp_nef : 0);   // photon
    _mpf["nef"][2] = (useP2CalX ? fcx_nef : 0); // SPR cross (HCAL+ECAL)
    _mpf["nef"][2] = (useP2HadX ? fhx_nef : 0); // SPR cross (HCAL)
    _mpf["nef"][3] = (useP3HadH ? fhm_nef : 0); // SPR-hcal
    _mpf["nef"][4] = (useP4HadE ? fem_nef : 0); // SPR-ecal
    _mpf["nef"][5] = (useP5Frag ? fhw_nef : 0);
    _mpf["nef"][6] = 0;//(useP6L1RC ? 0 : 0);
    _mpf["nef"][7] = 0;//(useP7TrkD ? 0 : 0);
    _mpf["nef"][8] = (useP8MB80pf ? fm80_nef : 0);

} // setToyShapeFuncs

void setFullShapeFuncs() {
  
      // Fits from minitools/fullSimShapes.C
    //////////////////////////////////////

    // Jet response (Rjet)
    // Fits from minitools/fullSimShapes.C
    TF1 *fhp3_Rjet = new TF1("fhp3_Rjet","[p0]+log(x)*([p1]+log(x)*([p2]+log(x)*([p3]+log(x)*[p4])))",15,4500);
    fhp3_Rjet->SetParameters(1.802,-1.629,0.496,-0.05204,0.001873); // 20.5/18
    TF1 *fhc3_Rjet = new TF1("fhc3_Rjet","[p0]+log(x)*([p1]+log(x)*([p2]+log(x)*([p3]+log(x)*([p4]+log(x)*([p5]+log(x)*[p6])))))",15,4500);
    fhc3_Rjet->SetParameters(-6.575,3.888,-0.4755,-0.06304,0.0063,0.001629,-0.0001509); // 0.0/16
    TF1 *fhm3_Rjet = new TF1("fhm3_Rjet","[p0]+log(x)*([p1]+log(x)*([p2]+log(x)*([p3]+log(x)*[p4])))",15,4500);
    fhm3_Rjet->SetParameters(-6.241,4.701,-1.277,0.1386,-0.005413); // 20.3/18
    //TF1 *fem3_Rjet = new TF1("fem3_Rjet","min(-(0.6+[p0])+(0.4+[p1])*pow(x,-(0.3+[p2])),-abs(0.06+[p3]))",15,4500);
    //fem3_Rjet->SetParameters(-0.00729,15.08,0.4843,-0.06); // 19.0/19
    TF1 *fem3_Rjet = new TF1("fem3_Rjet","[p0]+log(x)*([p1]+log(x)*([p2]+log(x)*([p3]+log(x)*([p4]+log(x)*([p5]+log(x)*[p6])))))",15,4500);
    fem3_Rjet->SetParameters(4.991,-1.863,0.07648,0.02939,-0.0003534,-0.0006707,4.937e-05); // 7.0/16
    TF1 *fpm3_Rjet = new TF1("fpm3_Rjet","[p0]+log(x)*([p1]+log(x)*[p2])",15,4500);
    fpm3_Rjet->SetParameters(-0.5754,-0.08557,0.007061); // 0.4/20
    TF1 *ftm3_Rjet = new TF1("ftm3_Rjet","[p0]-(1.0+[p1])*pow(x,-(0.3+[p2]))",15,4500);
    ftm3_Rjet->SetParameters(0.4822,3.867,-0.01323); // 34.4/20

    // Fits from minitools/fullSimShapes.C
    //////////////////////////////////////

    // Charged hadron fraction (CHF)
    // Fits from minitools/fullSimShapes.C
    TF1 *fhp3_chf = new TF1("fhp3_chf","[p0]+log(x)*([p1]+log(x)*([p2]+log(x)*([p3]+log(x)*[p4])))",15,4500);
    fhp3_chf->SetParameters(-4.448,3.116,-0.7725,0.07301,-0.002181); // 130.1/18
    TF1 *fhc3_chf = new TF1("fhc3_chf","[p0]+log(x)*([p1]+log(x)*([p2]+log(x)*([p3]+log(x)*([p4]+log(x)*([p5]+log(x)*[p6])))))",15,4500);
    fhc3_chf->SetParameters(3.42,-2.059,0.2634,0.03367,-0.003686,-0.0009108,8.937e-05); // 0.0/16
    TF1 *fhm3_chf = new TF1("fhm3_chf","[p0]+log(x)*([p1]+log(x)*([p2]+log(x)*([p3]+log(x)*[p4])))",15,4500);
    fhm3_chf->SetParameters(2.493,-1.765,0.4322,-0.03585,0.0006989); // 187.8/18
    TF1 *fem3_chf = new TF1("fem3_chf","[p0]+log(x)*([p1]+log(x)*([p2]+log(x)*([p3]+log(x)*[p4])))",15,4500);
    fem3_chf->SetParameters(0.04343,-0.1138,0.02241,0.005785,-0.0008414); // 32.9/18
    TF1 *fpm3_chf = new TF1("fpm3_chf","[p0]+log(x)*([p1]+log(x)*[p2])",15,4500);
    fpm3_chf->SetParameters(0.05363,0.1341,-0.01254); // 12.6/20
    TF1 *ftm3_chf = new TF1("ftm3_chf","[p0]+log(x)*([p1]+log(x)*([p2]+log(x)*([p3]+log(x)*[p4])))",15,4500);
    ftm3_chf->SetParameters(-9.824,6.347,-1.589,0.1567,-0.004988); // 5287.9/16

    // Fits from minitools/fullSimShapes.C
    //////////////////////////////////////

    // Neutral hadron fraction (NHF)
    // Fits from minitools/fullSimShapes.C
    TF1 *fhp3_nhf = new TF1("fhp3_nhf","[p0]+log(x)*([p1]+log(x)*([p2]+log(x)*([p3]+log(x)*([p4]+log(x)*([p5]+log(x)*[p6])))))",15,4500);
    fhp3_nhf->SetParameters(11.9,-10.71,3.249,-0.2538,-0.05074,0.01027,-0.0005032); // 361.3/16
    TF1 *fhc3_nhf = new TF1("fhc3_nhf","[p0]+log(x)*([p1]+log(x)*([p2]+log(x)*([p3]+log(x)*([p4]+log(x)*([p5]+log(x)*[p6])))))",15,4500);
    fhc3_nhf->SetParameters(-2.475,1.441,-0.177,-0.02317,0.002375,0.0006077,-5.567e-05); // 0.0/16
    TF1 *fhm3_nhf = new TF1("fhm3_nhf","[p0]+log(x)*([p1]+log(x)*([p2]+log(x)*([p3]+log(x)*([p4]+log(x)*([p5]+log(x)*[p6])))))",15,4500);
    fhm3_nhf->SetParameters(-9.651,9.03,-2.85,0.2346,0.04518,-0.009451,0.0004705); // 486.9/16
    TF1 *fem3_nhf = new TF1("fem3_nhf","[p0]+log(x)*([p1]+log(x)*([p2]+log(x)*([p3]+log(x)*[p4])))",15,4500);
    fem3_nhf->SetParameters(-0.3027,0.1342,-0.005781,-0.003766,0.0003855); // 73.3/18
    TF1 *fpm3_nhf = new TF1("fpm3_nhf","[p0]+log(x)*([p1]+log(x)*[p2])",15,4500);
    fpm3_nhf->SetParameters(0.2565,-0.07617,0.008004); // 9.9/20
    TF1 *ftm3_nhf = new TF1("ftm3_nhf","[p0]+log(x)*([p1]+log(x)*([p2]+log(x)*([p3]+log(x)*[p4])))",15,4500);
    ftm3_nhf->SetParameters(-3.537,4.217,-1.418,0.1935,-0.0092); // 9175.7/18

    // Fits from minitools/fullSimShapes.C
    //////////////////////////////////////

    // Photon fraction (NEF)
    // Fits from minitools/fullSimShapes.C
    TF1 *fhp3_nef = new TF1("fhp3_nef","[p0]+log(x)*([p1]+log(x)*([p2]+log(x)*([p3]+log(x)*([p4]+log(x)*([p5]+log(x)*[p6])))))",15,4500);
    fhp3_nef->SetParameters(-13.07,12.29,-3.883,0.3423,0.05416,-0.01187,0.0005907); // 182.2/16
    TF1 *fhc3_nef = new TF1("fhc3_nef","[p0]+log(x)*([p1]+log(x)*([p2]+log(x)*([p3]+log(x)*([p4]+log(x)*([p5]+log(x)*[p6])))))",15,4500);
    fhc3_nef->SetParameters(-1.027,0.6699,-0.09337,-0.01124,0.001406,0.0003225,-3.572e-05); // 0.0/16
    TF1 *fhm3_nef = new TF1("fhm3_nef","[p0]+log(x)*([p1]+log(x)*([p2]+log(x)*([p3]+log(x)*([p4]+log(x)*([p5]+log(x)*[p6])))))",15,4500);
    fhm3_nef->SetParameters(15.37,-14.03,4.332,-0.3683,-0.06156,0.01315,-0.0006485); // 418.9/16
    TF1 *fem3_nef = new TF1("fem3_nef","[p0]+log(x)*([p1]+log(x)*([p2]+log(x)*([p3]+log(x)*[p4])))",15,4500);
    fem3_nef->SetParameters(0.5027,-0.1245,-0.01303,0.0006432,0.0002053); // 117.3/18
    TF1 *fpm3_nef = new TF1("fpm3_nef","[p0]+log(x)*([p1]+log(x)*[p2])",15,4500);
    fpm3_nef->SetParameters(-0.2875,-0.06739,0.005384); // 16.9/20
    TF1 *ftm3_nef = new TF1("ftm3_nef","[p0]+log(x)*([p1]+log(x)*([p2]+log(x)*([p3]+log(x)*[p4])))",15,4500);
    ftm3_nef->SetParameters(21.13,-15.9,4.355,-0.4987,0.02021); // 5355.5/18

    ft  = ftm3_Rjet;
    fhx = fhc3_Rjet;
    fhh = fhm3_Rjet;
    feh = fem3_Rjet;
    fp  = fpm3_Rjet;

  // Copied from setToyShapeFuncs()
  //////////////////////////////////
  TF1 *fm80_chf = new TF1("fm80_chf","[p0]+[p1]*pow(x,[p2])",15,4500);
  fm80_chf->SetParameters(1.848, -1.685, -0.006643); // 70.3/53
  TF1 *fhw_nhf = new TF1("fhw_nhf","[p0]+[p1]*pow(x,[p2])+[p3]/x",15,4500);
  fhw_nhf->SetParameters(-5.151, 4.495, 0.03335, -12.3); // 65.4/42

  TF1 *fm80_nef = new TF1("fm80_nef","[p0]+[p1]*pow(x,[p2])",15,4500);
  fm80_nef->SetParameters(3.611, -3.909, -0.007482); // 70.2/53
  TF1 *fhw_nef = new TF1("fhw_nef","[p0]+[p1]*pow(x,[p2])+[p3]/x",15,4500);
  fhw_nef->SetParameters(0.8417, -0.2605, 0.2289, 2.426); // 35.9/42
  
  TF1 *fm80_nhf = new TF1("fm80_nhf","[p0]+[p1]*pow(x,[p2])",15,4500);
  fm80_nhf->SetParameters(-0.05047, -0.0008452, 0.6402); // 56.7/53
  TF1 *fhw_chf = new TF1("fhw_chf","[p0]+[p1]*pow(x,[p2])+[p3]/x",15,4500);
  fhw_chf->SetParameters(-0.2176, 1.064e-05, 1.373, 0); // 41.8/43

  // Regular +3% SPR calorimeter scale variation (ECAL+HCAL)
  if (!fc) fc = new TF1("fc","[p0]+[p1]*pow(x/[p2],[p3])/(1+pow(x/[p2],[p3]))*(1-pow(x/[p2],-[p3]))",15,4500);
  fc->SetParameters(0.335, 0.8171, 406.8, 1.34); // toyPF

  // Tracking '-1%' variation in UL2017 data
  if (!ftd) ftd = new TF1("ftd","[p0]+[p1]*pow(x/208.,[p2])+[p3]/x",15,4500);
  ftd->SetParameters(-0.116, -0.6417, -0.3051, 23.63); // Data

  // Tracking '-1%' variation in UL2017 MC
  if (!ftm) ftm = new TF1("ftm","[p0]+[p1]*pow(x/208.,[p2])+[p3]/x",15,4500);
  ftm->SetParameters(0.2683, -0.6994, -0.3051, 18.49); // MC

  // sigmaMB 69.2 mb to 80 mb variation
  if (!fm80) fm80 = new TF1("fm80","[p0]+[p1]*pow(x/208.,[p2])+[p3]*exp(-[p4]*x)",15,4500);
  fm80->SetParameters(0.03904,-0.01612,-0.942, -7.145,0.09907); // toyPF

  // H++ vs P8CP5
  if (!fhw) fhw = new TF1("fhw","[p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]/x+[p5]*log(x)/x",15,4500);
  fhw->SetParameters(0.9526,-0.3883,1285,2.46,18.1,-2.062); // FullMC

  assert(!useP2CalX);
  //useP6L1RC = false;
  //useP7TrkD = false;
  //useP8Flv = false;
  _mpf["chf"][0] = (useP0Trk  ? ftm3_chf : 0);// tracking in data+MC
  _mpf["chf"][1] = (useP1Gam  ? fpm3_chf : 0);  // photon
  _mpf["chf"][2] = (useP2HadX ? fhc3_chf : 0); // SPR cross (custom HCAL #3)
  _mpf["chf"][3] = (useP3HadH ? fhm3_chf : 0); // SPR-hcal
  _mpf["chf"][4] = (useP4HadE ? fem3_chf : 0); // SPR-ecal
  _mpf["chf"][5] = (useP5Frag ? fhw_chf : 0); // fragmentation
  _mpf["chf"][6] = 0;//(useP6L1RC ? 0 : 0);
  _mpf["chf"][7] = 0;//(useP7TrkD ? 0 : 0);
  _mpf["chf"][8] = (useP8MB80pf ? fm80_chf : 0); // fragmentation
  //
  _mpf["nhf"][0] = (useP0Trk  ? ftm3_nhf : 0);// tracking in data+MC
  _mpf["nhf"][1] = (useP1Gam  ? fpm3_nhf : 0);  // photon
  _mpf["nhf"][2] = (useP2HadX ? fhc3_nhf : 0); // SPR cross (HCAL)
  _mpf["nhf"][3] = (useP3HadH ? fhm3_nhf : 0); // SPR-hcal
  _mpf["nhf"][4] = (useP4HadE ? fem3_nhf : 0); // SPR-ecal
  _mpf["nhf"][5] = (useP5Frag ? fhw_nhf : 0);
  _mpf["nhf"][6] = 0;//(useP6L1RC ? 0 : 0);
  _mpf["nhf"][7] = 0;//(useP7TrkD ? 0 : 0);
  _mpf["nhf"][8] = (useP8MB80pf ? fm80_nhf : 0);
  //
  _mpf["nef"][0] = (useP0Trk ? ftm3_nef : 0); // tracking in data+MC
  _mpf["nef"][1] = (useP1Gam ? fpm3_nef : 0);   // photon
  _mpf["nef"][2] = (useP2HadX ? fhc3_nef : 0); // SPR cross (HCAL)
  _mpf["nef"][3] = (useP3HadH ? fhm3_nef : 0); // SPR-hcal
  _mpf["nef"][4] = (useP4HadE ? fem3_nef : 0); // SPR-ecal
  _mpf["nef"][5] = (useP5Frag ? fhw_nef : 0);
  _mpf["nef"][6] = 0;//(useP6L1RC ? 0 : 0);
  _mpf["nef"][7] = 0;//(useP7TrkD ? 0 : 0);
  _mpf["nef"][8] = (useP8MB80pf ? fm80_nef : 0);
} // setFullShapeFuncs
