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

using namespace std;

bool dofsr = true; // correct for FSR
double ptreco_gjet = 15.; // min jet pT when evaluating alphamax for gamma+jet
double ptreco_zjet = 5.; // same for Z+jet
bool dol1bias = false; // correct MPF for L1L2L3-L1 (instead of L1L2L3-RC)
bool _paper = true;
double _cleanUncert = 0.05; // for eta>2
//double _cleanUncert = 0.020; // Clean out large uncertainty points from PR plot
//bool _g_dcsonly = false;
string scalingForL2OrL3Fit = "None"; //"None" - for inpunt combination files without any residual applied
//"PutBackL2Res" - put L2res back in for gamma/Z+jet for vs eta studies
//N.B.: Barrel JES from input text file is always applied to dijet results

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
bool useOff = true; // pT-dependent offset
bool useTDI = false; // tracker dynamic inefficiency for 3p fit
bool useEG = false; // ECAL gain shift for 3p fit

//const int njesFit = 1; // scale only
const int njesFit = 2; // scale(ECAL)+HB
//const int njesFit = 3; //useOff=true; // scale(ECAL)+HB+offset
//const int njesFit = 3; //useTDI=true; // scale(ECAL)+HB+tracker (switchable)
//const int njesFit = 3; //useEG=true; // scale(ECAL)+HB+ECALgain
//const int njesFit = 4; // scale(ECAL)+HB+offset+ECALgain

const double ptminJesFit = 30;
TF1 *fhb(0), *fl1(0), *ftr(0), *feg(0); double _etamin(0);
Double_t jesFit(Double_t *x, Double_t *p) {

  double pt = *x;

  // Initialize SinglePionHCAL shape
  if (!fhb) fhb = new TF1("fhb","max(0.,[0]+[1]*pow(x,[2]))",10,3500);
  fhb->SetParameters(1.03091e+00, -5.11540e-02, -1.54227e-01); // SPRH
 
  // Initialize L1FastJet-L1RC difference
  if (!fl1) fl1 = new TF1("fl1","1-([0]+[1]*log(x))/x",10,3500);
  //fl1->SetParameters(0.553999, -8.10939e-03);
  //fl1->SetParameters(0., 1.89377e-01); // 80XV1re+Sum16 hl1bias (log(pT) only)
  //fl1->SetParameters(0.5, 0.); // Z+jet UE
  //fl1->SetParameters(1.79089e-01, 1.89377e-01); // 80XV1re+Sum16 hl1bias
  fl1->SetParameters(2.60382e-01, 1.96664e-01); // Sum16V6G hl1bias

  // Initialize tracker inefficiency shape
  //if (!ftr) ftr = new TF1("ftr","0.85-[0] + (0.15+[0])*([3]-[1]*pow(x,[2]))"
  //		  " + [4]/x",10,3500);
  //ftr->SetParameters(0.131,1.348,-0.241,1.261,2.649);
  // Values from drawAvsB.C for FvsG
  //ftr->SetParameters(0.3389, 1.313, -0.03876, 2.012, 1.747);
  // Values from drawAvsB.C for EvsG
  // The p3 is turned off for flatter low pT and steeper high pT, which
  // is more consistent with multijet data (not used in the fit, yet)
  if (!ftr) ftr = new TF1("ftr","1-[0]-[1]*pow(x,[2]) + ([3]+[4]*log(x))/x",10,3500);
  ftr->SetParameters(-0.04432, 1.304, -0.4624, 0, 1.724);

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
  if (njesFit==2) {

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
  } // njesFit==3 && TDI

  if (njesFit==3 && useTDI) {

    // p[0]: overall scale shift, p[1]: HCAL shift in % (full band +3%)
    double ptref = 208; // pT that minimizes correlation in p[0] and p[1]
    
    return (p[0] + p[1]/3.*100*(fhb->Eval(pt)-fhb->Eval(ptref))
	    + p[2]*(ftr->Eval(pt)-ftr->Eval(ptref)));
  } // njesFit==3 && TDI

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

}

void globalFitL3Res(double etamin = 0, double etamax = 1.3,
		    string epoch="") {
  if(verboseGF)cout << Form("Running globalFitL3Res(etamin=%02.2f,etamax=%02.2f,epoch=%s",etamin,etamax,epoch.c_str()) << endl << flush;
  _etamin = etamin;
  const char *cep = epoch.c_str();
  //njesFit = (njesFit==3 && useTDI && (epoch=="G"||epoch=="H") ? 2 : njesFit);

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
  const int nmethods = 2;
  const char* methods[nmethods] = {"ptchs","mpfchs1"};
  //  const int nmethods = 1;//MPFOnlyTest for FineEtaBins
  //  const char* methods[nmethods] = {"mpfchs1"};
  _nmethods = nmethods; // for multijets in global fit

  // Global fit with multijets, gamma+jet, Z+jet

  bool isl3 = (etamin==0 && ((epoch!="L4" && fabs(etamax-1.3)<0.1) ||
			     (epoch=="L4" && fabs(etamax-2.4)<0.1)));

  
  // Normal global fit with all four samples (multijet/dijet, gamma+jet, Z+jets)
  const int nsamples = 4;
  const int nsample0 = 1; // first Z/gamma+jet sample
  const char* samples[4] = {(isl3 ? "multijet" : "dijet"),
			    "gamjet", "zeejet", "zmmjet"};
  const int igj = 0;
  const int izee = 1;
  const int izmm = 2;
  

  /*
  // Global fit with only dijet, Z+jets
  const int nsamples = 3;
  const int nsample0 = 1; // first Z/gamma+jet sample
  const char* samples[4] = {"dijet",
			    "zeejet", "zmmjet"};
  const int igj = -1;
  const int izee = 1;
  const int izmm = 2;
  */

  /*
  // Global fit with only dijet
  const int nsamples = 1;
  const int nsample0 = 0; // first Z/gamma+jet sample
  const char* samples[1] = {"dijet"};
  const int igj = -1;
  const int izee = -1;
  const int izmm = -1;
  */

  /*
  // Global fit without multijets/dijets
  const int nsamples = 3;
  const int nsample0 = 0; // first Z/gamma+jet sample
  const char* samples[3] = {"gamjet", "zeejet", "zmmjet"};
  const int igj = 0;
  const int izee = 1;
  const int izmm = 2;
  */

  /*
  // Global fit without photon+jet
  const int nsamples = 3;
  const int nsample0 = 1; // first Z/gamma+jet sample
  const char* samples[3] = {(etamin==0 && (etamax==1.3||etamax==2.4)? "multijet" : "dijet"),
			    "zeejet", "zmmjet"};
  const int igj = -1;
  const int izee = 1;
  const int izmm = 2;
  */
  
  /*
  // Global fit with Z+jet only
  const int nsamples = 2;
  const int nsample0 = 0; // first Z/gamma+jet sample
  const char* samples[2] = {"zeejet", "zmmjet"};
  const int igj = -1;
  const int izee = 0;
  const int izmm = 1;
  */
  

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
  for (int imethod = 0; imethod != nmethods; ++imethod) {
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
	    (string(cs)=="dijet" && g->GetX()[i]<70.)) // ||
	    //(string(cs)=="gamjet" && g->GetX()[i]>600. && etamin!=0))
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
  for (int imethod = 0; imethod != nmethods; ++imethod) {
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
	double ptreco = (string(cs)=="zeejet"||string(cs)=="zmmjet" ?
			 ptreco_zjet : ptreco_gjet);
	double aeff = max(alpha, ptreco/pt);
	double kfsr = (dofsr ? 1./(1+aeff*h->GetBinContent(h->FindBin(pt))):1);
	double l1 = (dol1bias && string(cm)=="mpfchs1" && string(cs)=="gamjet" ?
		     1./hl1->GetBinContent(hl1->FindBin(pt)) : 1);

	// 0.98 was on until Aug 12, 10:12am :(
	double scale = 1.00;//0.98; // correct out previous L3Res
	//double scale = (string(cs)=="gamjet" ? 1.000 : 1);
	// Scales determined by hand:
	// BCD (+0.00-274.6)
	// EF (+0.50-214.9)
	// G (+1.00-224.0)
	// H (weird chi2, keep 0.00-291.4)
	// BCD: +0.00-274.6, -0.50-275.8, +0.50-275.0, -0.25-275.0, +0.25-274.6
	// EF: +1.50-218.0, +1.00-215.7, +0.50-214.9, +0.00-215.7, +0.25-215.1
	// G: +1.00-224.0, +0.50-224.5, +1.50-225.1, +0.75-224.1, +1.25-224.4
	// H: +0.00-291.4, -0.50-294.8, +0.50-289.6, +1.0-289.4

	double escale = 0; // E/p moved to reprocess.C
	// Add photon E/p correction (and it's uncertianty) on the fly
	//double scale = (string(cs)=="gamjet" ?
	//		heop->Interpolate(pt) : 1);
	//double escale = (string(cs)=="gamjet" ?
	//		heope->Interpolate(pt) : 0.);

	// put L2res back in for gamma/Z+jet for vs eta studies
	if (!(etamin==0 && fabs(etamax-1.3)<0.1) &&
	    string(cs)!="dijet" && string(cs)!="multijet") {
	  //int ijes = min(max(1,hjes->FindBin(pt)),hjes->GetNbinsX());
	  //scale = hjes->GetBinContent(ijes);
	  double jes = hjes->GetBinContent(hjes->FindBin(pt)); // L2Res only
	  double jesref = hjes0->GetBinContent(hjes0->FindBin(pt)); // barrel
	  // divide by jesref(=1) in case hjes had L2L3Res instead of L2Res
          if(scalingForL2OrL3Fit=="None")scale= 1.0;//no residual input
          else if(scalingForL2OrL3Fit=="PutBackL2Res")scale = jes / jesref;
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
	if (string(cs)=="dijet") {
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
    for (int imethod = 0; imethod != nmethods; ++imethod) {
      for (int isample = 0; isample != nsamples; ++isample) {
	
	const char *cm = methods[imethod];
	const char *cs = samples[isample];

	string s = Form("fsr/hkfsr_%s_%s_eig%d",cm,cs,ieig);
	TH1D *h = (TH1D*)d->Get(s.c_str());
	if (!h) cout << "Histo "<<s<<" not found!" << endl << flush;
	assert(h);

	int ibm = isample + nsamples*imethod;	
	//h->Scale(alpha); // set which nominal points this applies to
	// Sum16V6: take into account pTmin threshold effect on alphamax
	for (int i = 1; i != h->GetNbinsX()+1; ++i) {
	  double pt = h->GetBinCenter(i);
	  double ptreco = (string(cs)=="zeejet"||string(cs)=="zmmjet" ?
			 ptreco_zjet : ptreco_gjet);
	  double aeff = max(alpha, ptreco/pt);
	  h->SetBinContent(i, aeff*h->GetBinContent(i));
	}
	h->SetName(Form("bm%d_%s",(1<<ibm),h->GetName()));

	// TEMP TEMP TEMP
	// Scale FSR uncertainty by x10 to see the impact on Run 2 fit
	//h->Scale(10.0);
	//if (!(string(cs)=="gamjet" && string(cm)=="ptchs"))
	//if (!(string(cs)=="gamjet"))
	//h->Scale(4.0);
	// TEMP TEMP TEMP

	hs.push_back(h);
      } // for isample
    } // for imethod
  } // for ieig

  // Uncertainty sources for e, gamma and mu scale uncertainties
  // (and extra source for multijets?)
  int is(0);
  int is_gj(0), is_zee(0), is_zmm(0);
  for (unsigned int i = 0; i != _vdt->size(); ++i) {
      
    string s = samples[i%nsamples]; const char *cs = s.c_str();
    unsigned int imethod = i/nsamples;
    string m = methods[imethod];    const char *cm = m.c_str();

    TH1D *h = hs[i]; assert(h);
    //TH1D *h2 = (TH1D*)h->Clone(Form("bm%d_%s_%s_scale_%d",1<<i,cm,cs,i));
    TH1D *h2 = (TH1D*)h->Clone(Form("bm%d_inactive_%s_%s_%d",1<<i,cm,cs,i));

    double escale(0);
    const int n0 = nsample0;
    const int n1 = nsamples+nsample0;
    if (s=="multijet" && imethod==0) {
      escale = 1.0; // no constraint on absolute scale
      h2->SetName(Form("bm%d_scale_multijet_%d",(1<<(n0-1) | (1<<(n1-1))), i));
    }
    if (s=="dijet" && imethod==0) {
      escale = 0.005; // for JER bias and whatnot
      h2->SetName(Form("bm%d_scale_dijet_%02.0f_%d",
		       (1<<(n0-1) | (1<<(n1-1))), escale*1000., i));
    }
    if (s=="gamjet" && imethod==0) {
      //escale = 0.005; // photon (ECAL) scale without regression
      //escale = 0.002; // photon (ECAL) scale with regression (Run I)
      //escale = 0.010; // Run II guess
      //escale = 0.020; // 80X_V7 until photon scale explicitly fixed
      //escale = 0.03;//0.005; // 80XV1re+Sum16
      //escale = 0.005; // 80XV1re+Sum16
      escale = 0.020; // 03Feb2016V2 TESTB
      //escale = 0.010; // Sep23V1
      // Use same source for both MPF and pT balance
      h2->SetName(Form("bm%d_scale_gamjet_%02.0f_%d",
		       (1<<(n0+igj) | (1<<(n1+igj))), escale*1000., i));
      is = hs.size();
      is_gj = hs.size();
    }
    /*
    // reduce photon scale uncertainty to 0.2% with regression corrections for pTbal
    // but keep 0.5% uncorrelated uncertainty for MPF with footprint
    // (until fully fixed)
    if (s=="gamjet" && imethod==0) { // ptchs
      escale = 0.002; // photon (ECAL) scale in pTbal
      h2->SetName(Form("bm%d_scale_%d",(1<<(n0+0)), i));
      is = hs.size();
    }
    if (s=="gamjet" && imethod==1) { // mpfchs1
      escale = 0.005; // photon (ECAL) scale in MPF
      h2->SetName(Form("bm%d_scale_%d",(1<<(n1+0)), i));
      is = hs.size();
    }
    */
    if (s=="zeejet" && imethod==0) {
      //escale = 0.005; // electron (ECAL) scale (Run I)
      //escale = 0.002; // with Zee mass scale fix
      //escale = 0.010; // Run II guess
      //escale = 0.005; // 76X guess
      //escale = 0.005; // 80X-BCD,E,F,G guess
      escale = 0.002; // 80X BCDEF+GH: flat -0.15% shift on pol0/pol1 pT<400
      // Use same source for both MPF and pT balance
      h2->SetName(Form("bm%d_scale_zeejet_%02.0f_%d",
		       (1<<(n0+izee) | (1<<(n1+izee))), escale*1000, i));
      is_zee = hs.size();
    }
    if (s=="zmmjet" && imethod==0) {
      // Use same source for both MPF and pT balance
      //escale = 0.002; // muon (tracker) scale (Run I)
      //escale = 0.010; // Run II guess
      //escale = 0.005; // 76X guess
      //escale = 0.002; // 80X-BCD,E,F,G guess
      //escale = 0.005; // 80XV1re+Sum16
      escale = 0.002; // 80X BCDEF+GH: flat -0.15% shift on pol0/pol1 pT<400
      //escale = 0.005; // Sep23V1-BCD,E,F,G guess
      h2->SetName(Form("bm%d_scale_zmmjet_%02.0f_%d",
		       (1<<(n0+izmm) | (1<<(n1+izmm))), escale*1000, i));
      is_zmm = hs.size();
    }

    // Same scale uncertainty applies to all pT bins
    for (int j = 1; j != h2->GetNbinsX()+1; ++j) {
      h2->SetBinContent(j, escale);
      h2->SetBinError(j, escale);
    } // for j

    hs.push_back(h2);
  } // for i

  // Create additional sources for MPF uncertainties with e/gamma
  // (one for each sample x method, but all except two are empty)
  for (unsigned int i = 0; i != _vdt->size(); ++i) {

    string s = samples[i%nsamples]; const char *cs = s.c_str();
    string m = methods[i/nsamples]; const char *cm = m.c_str();

    double escale(0);
    //int n1 = nsamples+nsample0;
    //if (i==n1+0) { escale = 0.005; } // gamma+jet
    //if (i==n1+1) { escale = 0.005; } // Zee+jet
    //if (s=="gamjet" && m=="mpfchs1") { escale = 0.005; } // GT
    // EM footprint uncertainties for photons and electrons
    //if (s=="gamjet" && m=="mpfchs1") { escale = 0.002; } // V8
    //if (s=="gamjet" && m=="mpfchs1") { escale = 0.010; } // 2015-10-21
    if (s=="gamjet" && m=="mpfchs1") { escale = 0.005; } // 76X
    //if (s=="zeejet" && m=="mpfchs1") { escale = 0.005; }
    //if (s=="zeejet" && m=="mpfchs1") { escale = 0.010; } // 2015-12-10
    if (s=="zeejet" && m=="mpfchs1") { escale = 0.005; } // 76X

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
	  //h2->SetName(Form("bm%d_l1bias_%d",(1<<(n0+0)|1<<(n0+1)|1<<(n0+2)|
	  //			     1<<(n1+0)|1<<(n1+1)|1<<(n1+2)),i));
	  h2->SetName(Form("bm%d_l1bias_%d",(1<<(n0+izee)|1<<(n0+izmm)|
					     1<<(n1+izee)|1<<(n1+izmm)),i));
	  // TEMP TEMP TEMP
	  // Only switch PileUpPt source for MPF (due to L1L2L3-L1 method)
	  //h2->SetName(Form("bm%d_l1bias_%d",(1<<(n1+0)|1<<(n1+1)|1<<(n1+2)),i));
	  // Hack to fix pT. dep of missing gamma+jet MPF footprint correction
	  //h2->SetName(Form("bm%d_l1bias_%d",1<<(n1+igj),i));
	  // TEMP TEMP TEMP
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
    for (int imethod = 0; imethod != nmethods; ++imethod) {
	
      const char *cm = methods[imethod];
      const char *cs = "multijet";
      
      //const int npt = 14; // 20161019
      //double vpt[npt+1] = {200, 250, 300, 370, 450, 510, 550, 600, 700, 800,
      //		   900, 1000, 1200, 1400, 1600};
      //const int npt = 23; // 20161021, 20161202
      //double vpt[npt+1] = {200, 250, 300, 330, 370, 410, 450, 510, 530, 550,
      //		   575, 600, 650, 700, 800, 900, 1000, 1100, 1200,
      //		   1400, 1600, 1800, 2000, 2200}; // V8
      const int npt = 24; // 20161222, 80XV1re+Sum16
      double vpt[npt+1]= {200, 250, 300, 330, 370, 410, 450, 510, 530, 550, 575,
			  600, 650, 700, 800, 900, 1000, 1100, 1200, 1400, 1600,
			  1800, 2000, 2200};// 2016Early (has 3000 also?)
      //const int npt = 67;//68; // 20161123, 80X rereco
      //double vpt[npt+1] = {200, 225, 250, 275, 300, 310, 320, 330, 340, 350,
      //		   360, 370, 380, 390, 400, 410, 420, 430, 440, 450,
      //		   460, 470, 480, 490, 500, 510, 520, 530, 540, 550,
      //		   560, 570, 580, 590, 600, 610, 620, 630, 640, 650,
      //		   660, 670, 680, 690, 700, 720, 740, 760, 780, 800,
      //		   820, 840, 860, 880, 900, 920, 940, 960, 980, 1000,
      //		   1050, 1100, 1200, 1400, 1600, 1800, 2000, 2200};
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
    for (int imethod = 0; imethod != nmethods; ++imethod) {
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

  const int maxpt = 3500;//1600;
  const int minpt = 30;
  TH1D *h = new TH1D("h",";p_{T} (GeV);Jet response (ratio)",
		     maxpt-minpt,minpt,maxpt);
  //h->SetMinimum(etamin==0 ? 0.93 : 0.80);
  //h->SetMaximum(etamin==0 ? 1.08 : 1.20);
  //h->SetMinimum(etamin>=3 ? 0.00 : (etamin>=1.9 ? 0.50 : 0.79)); // 0.82
  //h->SetMaximum(etamin>=3 ? 2.00 : (etamin>=1.9 ? 1.50 : 1.24)); // 1.20
  //h->SetMinimum(etamin>=3 ? 0.50 : (etamin>=2.5 ? 0.70 : 0.85)); // 74X V7
  //h->SetMaximum(etamin>=3 ? 1.75 : (etamin>=2.5 ? 1.45 : 1.20)); // 74X V7
  //h->SetMinimum(etamin>=3 ? 0.50 : (etamin>=2.5 ? 0.70 : 0.93)); // 74X V7
  //h->SetMaximum(etamin>=3 ? 1.75 : (etamin>=2.5 ? 1.45 : 1.08)); // 74X V7
  //h->SetMinimum(etamin>=3 ? 0.50 : (etamin>=2.5 ? 0.70 : 0.92)); // 76X V1
  //h->SetMaximum(etamin>=3 ? 1.75 : (etamin>=2.5 ? 1.45 : 1.12)); // 76X V1
  h->SetMinimum(etamin>=3 ? 0.50 : (etamin>=2.5 ? 0.70 : 0.91)); // 76X V1
  h->SetMaximum(etamin>=3 ? 1.75 : (etamin>=2.5 ? 1.45 : 1.15)); // 76X V1
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  h->DrawClone("AXIS");

  //////////////////////////////////////////////
  // for pas-v6, recenter uncertainty around fit
  //TF1 *jesshift = new TF1("jesshift",jesFit,minpt,maxpt,njesFit);
  //jesshift->SetParameters(0.9784, -0.0362, 0.054); // AN2015_023_v2 Table 1
  //for (int i = 1; i != herr_ref->GetNbinsX()+1; ++i) {
  //double x = herr_ref->GetBinCenter(i);
  //if (etamin==0 && _paper) herr_ref->SetBinContent(i, jesshift->Eval(x));
  //}
  //delete jesshift;
  ///////////////////////////////////////////////

  //TCanvas *c0 = new TCanvas("c0","c0",600,600);
  //lumi_13TeV = "42 pb^{-1}";
  //lumi_13TeV = "25 ns - 16 pb^{-1}";
  //lumi_13TeV = "Run2015D - 25 ns - 106.5 pb^{-1}";
  //lumi_13TeV = "Run2015D - 25 ns - 121 pb^{-1}";
  //lumi_13TeV = (_g_dcsonly ? 
  //		"Run2015D-DCSOnly - 25 ns - 980 pb^{-1}" :
  //		"Run2015D-Cert - 25 ns - 552 pb^{-1}");
  //lumi_13TeV = "Run2015D - Oct 19 - 1.28 fb^{-1}";
  //lumi_13TeV = "Run2015D - Golden - 2.1 fb^{-1}";
  //lumi_13TeV = "Run2016, 590 pb^{-1}";
  //lumi_13TeV = "Run2016, 0.6 fb^{-1}";
  //lumi_13TeV = "Run2016, 2.1 fb^{-1}";
  //lumi_13TeV = "Run2016, 2.6 fb^{-1}";
  //lumi_13TeV = "Run2016G, 4.4 fb^{-1}";
  //lumi_13TeV = "Run2016BCD, 13 fb^{-1}";
  //lumi_13TeV = "Run2016, 2.6+1.5 fb^{-1}";
  map<string, const char*> lumimap;
  lumimap["BCD"] = "Run2016BCD re-mAOD, 12.9 fb^{-1}";
  lumimap["E"] = "Run2016E re-mAOD, 4.0 fb^{-1}";
  lumimap["F"] = "Run2016F re-mAOD, 2.8 fb^{-1}";//3.1 fb^{-1}";
  //lumimap["G"] = "Run2016G re-mAOD, 7.1 fb^{-1}";
  lumimap["G"] = "Run2016FG re-mAOD, 8.0 fb^{-1}";
  lumimap["H"] = "Run2016H re-mAOD, 8.8 fb^{-1}";
  lumimap["GH"] = "Run2016FGH re-mAOD, 16.8 fb^{-1}";
  lumimap["BCDEF"] = "Run2016BCDEF re-mAOD, 19.7 fb^{-1}";
  lumimap["BCDEFGH"] = "Run2016BCDEFGH re-mAOD, 36.5 fb^{-1}";
  lumimap["EF"] = "Run2016EF re-mAOD, 6.8 fb^{-1}";
  lumimap["L4"] = "Run2016BCDEFGH closure, 36.5 fb^{-1}";
  lumi_13TeV = lumimap[epoch];

  TCanvas *c0 = tdrCanvas("c0",h,4,11,true);
  gPad->SetLogx();
  
  //cmsPrel(19800.);
  //CMS_lumi(c0, 2, 0);

  // multijets up/down
  TLegend *legpf = tdrLeg(0.58,0.79,0.88,0.83);
  legpf->AddEntry(gs[0]," ","PL");
  if (_vpt2->size()!=0) legpf->AddEntry((*_vpt2)[0]," ","PL");
  TLegend *legmf = tdrLeg(0.66,0.79,0.96,0.83);
  legmf->AddEntry(gs[nsamples]," ","PL");
  if (_vpt2->size()>1) legmf->AddEntry((*_vpt2)[1]," ","PL");

  TLegend *legp = tdrLeg(0.58,0.55,0.88,0.90);
  legp->SetHeader("p_{T}^{bal}");
  for (int i = 0; i != nsamples; ++i)
    legp->AddEntry(gs[i]," ",i==0 ? "" : "PL");

  map<string, const char*> texlabel;
  texlabel["multijet"] = "Multijet";
  texlabel["dijet"] = "Dijet";
  texlabel["gamjet"] = "#gamma+jet";
  texlabel["zeejet"] = "Zee+jet";// TB";
  //texlabel["zeejet"] = "Zee+jet RS";
  texlabel["zmmjet"] = "Z#mu#mu+jet";// TB";
  //texlabel["zmmjet"] = "Z#mu#mu+jet RS";
  TLegend *legm = tdrLeg(0.66,0.55,0.96,0.90);
  legm->SetHeader("MPF");
  for (int i = 0; i != nsamples; ++i)
    legm->AddEntry(gs[i+nsamples],texlabel[samples[i]],i==0 ? "" : "PL");

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
    tex->DrawLatex(0.20,0.73,Form("%1.1f#leq|#eta|<%1.1f",etamin,etamax));
  }
  tex->DrawLatex(0.20, 0.22, "Before global fit");

  hrun1->SetLineWidth(2);
  hrun1->SetLineColor(kCyan+3);
  hrun1->SetLineStyle(kDashed);

  herr_ref->SetLineWidth(2);
  herr_ref->SetLineColor(kYellow+3);
  herr_ref->SetLineStyle(kDashed);

  hrun1->DrawClone("SAME E5");
  (new TGraph(hrun1))->DrawClone("SAMEL");
  herr_ref->DrawClone("SAME E5");
  (new TGraph(herr_ref))->DrawClone("SAMEL");

  legp->AddEntry(hrun1," ","");
  //legm->AddEntry(hrun1,"Run 1","FL");
  //legm->AddEntry(hrun1,"80X V6","FL");
  legm->AddEntry(hrun1,"Run I","FL");

  legp->AddEntry(herr_ref," ","");
  //legm->AddEntry(herr_ref,"JES unc.","FL");
  //legm->AddEntry(herr_ref,"80X V7","FL");
  //legm->AddEntry(herr_ref,"80X V8","FL");
  //legm->AddEntry(herr_ref,"80XreV1","FL");
  //legm->AddEntry(herr_ref,"Sum16V3","FL");
  //  legm->AddEntry(herr_ref,"03FebV3","FL");
  legm->AddEntry(herr_ref,"03FebVXX","FL");

  hrun1->SetFillStyle(kNone);
  hrun1->DrawClone("SAME E5");
  (new TGraph(hrun1))->DrawClone("SAMEL");
  hrun1->SetFillStyle(1001);

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

  ///////////////////////
  // Draw raw response
  //////////////////////

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
  else tex->DrawLatex(0.20,0.73,Form("%1.1f#leq|#eta|<%1.1f",etamin,etamax));

  tex->DrawLatex(0.20, 0.22, "Before FSR correction");

  hrun1->DrawClone("SAME E5");
  (new TGraph(hrun1))->DrawClone("SAMEL");

  herr_ref->DrawClone("SAME E5");
  (new TGraph(herr_ref))->DrawClone("SAMEL");

  hrun1->SetFillStyle(kNone);
  hrun1->DrawClone("SAME E5");
  (new TGraph(hrun1))->DrawClone("SAMEL");
  hrun1->SetFillStyle(1001);

  herr_ref->SetFillStyle(kNone);
  herr_ref->DrawClone("SAME E5");
  (new TGraph(herr_ref))->DrawClone("SAMEL");
  herr_ref->SetFillStyle(1001);

  for (unsigned int i = 0; i != _vdt3->size(); ++i) {
    
    TGraphErrors *g3 = (*_vdt3)[i];
    // Skip multijet downward points
    g3->DrawClone("SAMEPz");
  }


  ///////////////////////
  // Perform global fit
  //////////////////////

  // Fit function
  // 2015-01-17: use only range visible on plots (one g+jet point at >2.5 TeV)
  TF1 *jesfit = new TF1("jesfit",jesFit,minpt,maxpt,njesFit);
  jesfit->SetLineColor(kBlack);
  if (njesFit==1) {
    jesfit->SetParameter(0, 0.98);
  }
  if (njesFit==2) {
    jesfit->SetParameters(0.98, 0.0);
  }
  if (njesFit==3 && useOff) {
    jesfit->SetParameters(0.98, 0.0, 0.0);
  }
  if (njesFit==3 && useTDI) {
    //jesfit->SetParameters(0.98, 0.0, 0.0);
    //jesfit->SetParameters(0.99, -0.03, 0.0);
    jesfit->SetParameters(0.985, 0.001, 0.5);
  }
  if (njesFit==3 && useEG) {
    //jesfit->SetParameters(0.98, 0.0, 1);
    jesfit->SetParameters(0.995, -0.025, 1.0);
  }
  if (njesFit==4 && useOff && useEG) {
    //jesfit->SetParameters(1.00, -0.01, 0.0, 1.0);
    jesfit->SetParameters(0.995, -0.025, 1.0, 1.0);
  }
  _jesFit = jesfit;
  

  // Get linear equations and solve them to get initial values for fitter
  const int np = _jesFit->GetNpar();
  const int nsrc = _vsrc->size();
  Int_t Npar = np+nsrc;

  cout << "Global fit has " << np << " fit parameters and "
       << nsrc << " nuisance parameters" << endl;

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
  Double_t chi2_gbl(0), chi2_src(0), chi2_data(0);
  vector<double> vchi2_data(gs2.size(),0);
  vector<int> vndata(gs2.size(),0);
  int nsrc_true(0);
  Double_t grad[Npar];
  Int_t flag = 1;
  TH1D *hsrc = new TH1D("hsrc",";Nuisance parameter;",30,-3,3);

  for (int i = 0; i != Npar; ++i) {
    tmp_par[i] = fitter->GetParameter(i);
    tmp_err[i] = fitter->GetParError(i);
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

  cout << endl;
  cout << "*** Processing eta bin " << etamin<<" - "<<etamax << " ***" << endl;
  cout << "Global chi2/ndf = " << chi2_gbl
       << " / " << Nk-np << " for " << Nk <<" data points, "
       << np << " fit parameters and "
       << nsrc << " ("<<nsrc_true<<") uncertainty sources" << endl;
  cout << endl;
  cout << "For data chi2/ndf = " << chi2_data << " / " << Nk << endl;
  cout << "For sources chi2/ndf = " << chi2_src << " / " << nsrc_true << endl;
  cout << "Per data set:" << endl;
  for (unsigned int i = 0; i != vchi2_data.size(); ++i) {
    cout << "  " << Form("%1.1f",vchi2_data[i]) << " / " << vndata[i]
	 << "  ("<<gs2[i]->GetName()<<")"<< endl;
  }

  ofstream txtL2L3("GlobalFitOutput_L2L3Residuals.txt",ios_base::app);
  if(etamin==0.&&etamax==0.261)txtL2L3 << "{ 1 JetEta 1 JetPt 1./([0]+[1]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[2]*((-2.36997+0.413917*TMath::Log(x))/x-(-2.36997+0.413917*TMath::Log(208))/208)) Correction L2Relative}";
  if(np==2&& !(etamin==0.&&etamax==1.3)){
    txtL2L3 << Form("\n %7.4f  %7.4f  5 10 6500 %7.4f %7.4f 1.0 ", etamin, etamax, tmp_par[0], tmp_par[1]);
    //    txtL2L3 << Form("%7.4f  %7.4f  5 10 6500 %7.4f %7.4f 1.0 \n ", -etamax, -etamin, tmp_par[0], tmp_par[1]);
  }


  for (int isample = 0; isample != nsamples; ++isample) {
    for (int imethod = 0; imethod != nmethods; ++imethod) {
      const char *cm = methods[imethod];
      const char *cs = samples[isample];
      string s = Form("fsr/hkfsr_%s_%s",cm,cs);
      TH1D *h = (TH1D*)d->Get(s.c_str());
      if (!h) cout << "Histo "<<s<<" not found!" << endl << flush;
      assert(h);

      ofstream txtFSRDiJet(Form("GlobalFitOutput_FSR_%s_%s.txt",cs,cm),ios_base::app);
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

  /*
  cout << "Uncertainty sources:" << endl;
  for (int i = np; i != Npar; ++i) {
    if ((i-np)%(2*nsamples)==0) cout << (*_vsrc)[i-np]->GetName() << endl;
    if (tmp_par[i]==0&&fabs(tmp_err[i]-1)<1e-2)
      cout << Form("%2d:   ----------,    ",i-np+1);
    else
      cout << Form("%2d: %5.2f+/-%4.2f,   ",i-np+1,tmp_par[i],tmp_err[i]);
    if ((i-np)%nsamples==nsamples-1) cout << endl;
  }
  cout << endl;
  */
  cout << "Uncertainty sources:" << endl;
  for (int i = np; i != Npar; ++i) {
    //if ((i-np)%(2*nsamples)==0) cout << (*_vsrc)[i-np]->GetName() << endl;
    if (tmp_par[i]==0&&fabs(tmp_err[i]-1)<1e-2)
      cout << Form("%2d - %35s:  ----------\n",
		   (i-np+1), (*_vsrc)[i-np]->GetName());
    else
      cout << Form("%2d - %35s:  %5.2f+/-%4.2f\n",
		   i-np+1, (*_vsrc)[i-np]->GetName(),
		   tmp_par[i],tmp_err[i]);
  //if ((i-np)%nsamples==nsamples-1) cout << endl;
  }
  cout << endl;


  ///////////////////////
  // Draw shifted data
  ///////////////////////

  h = (TH1D*)h->Clone("h1");
  //h->SetYTitle("Shifted jet response (ratio)");  
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
  else tex->DrawLatex(0.20,0.73,Form("%1.1f#leq|#eta|<%1.1f",etamin,etamax));

  //tex->DrawLatex(0.20,0.18,Form("(data %1.1f / %d, sources %1.1f / %d)",
  //			chi2_data,Nk,chi2_src,nsrc_true));

  tex->DrawLatex(0.20,0.22,"After global fit");
  tex->DrawLatex(0.20,0.17,Form("#chi^{2} / NDF = %1.1f / %d",
				chi2_gbl, Nk-np));

  is += np;

  // Fitted lepton/photon scales too much detailed for the paper
  tex->SetTextSize(0.030); tex->SetTextColor(kBlue-9); // 0.68,0.64,0.60
  /*
  tex->DrawLatex(0.20,0.70,Form("#gamma #times (%1.3f#pm%1.3f)",
				1+0.01*tmp_par[is_gj],
				0.01*sqrt(emat[is_gj][is_gj])));
				//1+0.01*tmp_par[is],
				//0.01*sqrt(emat[is][is])));
  tex->DrawLatex(0.20,0.67,Form("e #times (%1.3f#pm%1.3f)",
				1+0.01*tmp_par[is_zee],
				0.01*sqrt(emat[is_zee][is_zee])));
  //1+0.005*tmp_par[is+1],
  //0.005*sqrt(emat[is+1][is+1])));
  tex->DrawLatex(0.20,0.64,Form("#mu #times (%1.3f#pm%1.3f)",
				1+0.01*tmp_par[is_zmm],
				0.01*sqrt(emat[is_zmm][is_zmm])));
  //1+0.005*tmp_par[is+2],
  //0.005*sqrt(emat[is+2][is+2])));
  */			

  tex->DrawLatex(0.20,0.61,Form("N_{par}=%d",_jesFit->GetNpar()));
  tex->DrawLatex(0.32,0.61,Form("p_{0}=%1.3f #pm %1.3f",
				_jesFit->GetParameter(0),
				sqrt(emat[0][0])));
  if (njesFit>=2)
    tex->DrawLatex(0.32,0.58,Form("p_{1}=%1.3f #pm %1.3f",
				  _jesFit->GetParameter(1),
				  sqrt(emat[1][1])));
  if (njesFit>=3)
    tex->DrawLatex(0.32,0.55,Form("p_{2}=%1.3f #pm %1.3f",
				  _jesFit->GetParameter(2),
				  sqrt(emat[2][2])));
  if (njesFit>=4)
    tex->DrawLatex(0.32,0.52,Form("p_{3}=%1.3f #pm %1.3f",
				  _jesFit->GetParameter(3),
				  sqrt(emat[3][3])));
  //tex->DrawLatex(0.32,0.55,Form("p_{2}=%1.3f #pm %1.3f",
  //			_jesFit->GetParameter(2),
  //			sqrt(emat[2][2])));
  tex->SetTextSize(0.045); tex->SetTextColor(kBlack);

  hrun1->DrawClone("SAME E5");
  (new TGraph(hrun1))->DrawClone("SAMEL");

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

  hrun1->SetFillStyle(kNone);
  hrun1->DrawClone("SAME E5");
  (new TGraph(hrun1))->DrawClone("SAMEL");
  hrun1->SetFillStyle(1001);

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

  // Draw constant fit as reference as well (Run 2)
  if (false) {

    double vp[6][2] = // 
      {{ 0.9409 , 0.0073},  // 0-1.3
       { 0.9352 , 0.0046},  // 1.3-1.9
       { 0.9160 , 0.0049},  // 1.9.2.5
       { 0.8755 , 0.0056},  // 2.5-3.0
       { 0.7568 , 0.0085},  // 3.0-3.2
       { 0.7556 , 0.0067}}; // 3.2-5.2

    double ix = int(etamin*10+0.5);
    double val(0), err(0);
    if (ix==0)  { val = vp[0][0]; err = vp[0][1]; }
    if (ix==13) { val = vp[1][0]; err = vp[1][1]; }
    if (ix==19) { val = vp[2][0]; err = vp[2][1]; }
    if (ix==25) { val = vp[3][0]; err = vp[3][1]; }
    if (ix==30) { val = vp[4][0]; err = vp[4][1]; }
    if (ix==32) { val = vp[5][0]; err = vp[5][1]; }

    TLine *l = new TLine();
    l->SetLineStyle(kDashed);
    l->SetLineColor(kGray+1);
    l->DrawLine(minpt, val, maxpt, val);
    l->SetLineStyle(kDashDotted);
    l->DrawLine(minpt, val+err, maxpt, val+err);
    l->DrawLine(minpt, val-err, maxpt, val-err);
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

  if (false) { // Overlay effect of time slew change

    TDirectory *dir = gDirectory;
    TFile *fts = new TFile("rootfiles/JEClosure_757p1_with753JEC.root","READ");
    assert(fts && !fts->IsZombie());

    TDirectory *d = fts->GetDirectory("ak4PFJetAnalyzer");
    assert(d);
    TGraphErrors *g0 = (TGraphErrors*)d->Get("grMean_ak4PF");
    assert(g0);
    
    TGraphErrors *g = (TGraphErrors*)g0->Clone();
    //for (int i = 0; i != g->GetN(); ++i) {
    //g->SetPoint(i, g->GetX()[i], 0.952+(1-g->GetY()[i]));
    //}

    c0->cd();
    tdrDraw(g,"PL",kFullStar,kBlack);
    
    c1->cd();
    tdrDraw(g,"PL",kFullStar,kBlack);

    fts->Close();
    dir->cd();
  }

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
  //int colors[nsamples] = {kBlue, kGreen+2,kRed};
  //int colors[nsamples] = {kBlue, kRed};
  //int colors[nsamples] = {kBlue};
  for (int imethod = 0; imethod != nmethods; ++imethod) {

    h->SetMinimum(-0.06);//imethod==0 ? -0.20 : -0.20);//-0.04);
    h->SetMaximum(+0.09);//imethod==0 ? +0.30 : +0.30);//+0.06);
    //h->SetYTitle("k_{FSR} = dR/d#alpha (ratio)");
    h->SetYTitle("dR / d#alpha (ratio)");
    h = (TH1D*)h->Clone(Form("h3_%d",imethod));

    TCanvas *c3 = tdrCanvas(Form("c3_%d",imethod),h,4,11,true);
    c3->SetLogx();

    TLegend *leg1 = tdrLeg(0.58,0.65,0.88,0.90);
    leg1->SetHeader("In");
    TLegend *leg2 = tdrLeg(0.65,0.65,0.95,0.90);
    leg2->SetHeader("Out");

    tex->DrawLatex(0.55,0.20,
		   imethod==0 ? "p_{T} balance method" : "MPF method");

    //    for (int isample = nsample0; isample != nsamples; ++isample) {
    int offsetsample=0;//nsample0; ==>patch to also show dijet fit
    for (int isample = offsetsample; isample != nsamples; ++isample) {

      int ibm = isample + nsamples*imethod;
      TH1D *hk = hks[ibm]; assert(hk);

      double minx = (isample==nsample0 ? 40  : 30);
      //double maxx = (isample==nsample0 ? 900 : 600);//800 : 300);
      double maxx = (isample==nsample0 ? 1500 : 1000); // 80X
      hk->GetXaxis()->SetRangeUser(minx,maxx);

      //if (isample==nsample0) {
      //hk->GetXaxis()->SetRangeUser(40,800);
      //}
      //if (isample!=nsample0) {
      //hk->GetXaxis()->SetRangeUser(30,300);
      //}

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
	double ptreco = (string(cs)=="zeejet"||string(cs)=="zmmjet" ?
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
} // globalFitL3Res


void jesFitter(Int_t& npar, Double_t* grad, Double_t& chi2, Double_t* par,
	       Int_t flag) {

  // Basic checks
  assert(_vdt);
  assert(_vsrc);
  assert(int(_vsrc->size())==npar-_jesFit->GetNpar()); // nuisance per source
  assert(_vdt->size()<=sizeof(int)*8); // bitmap large enough

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

    // Add chi2 from nuisance parameters
    for (int is = 0; is != ns; ++is) {
      
      chi2 += ps[is]*ps[is];
    } // for ipar
    
    // Give some feedback on progress in case loop gets stuck
    if ((++cnt)%1000==0) cout << "." << flush;
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
