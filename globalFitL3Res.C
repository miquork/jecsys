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
#include <iostream>
#include <vector>
#include <map>

using namespace std;

bool dofsr = true; // correct for FSR
bool dol1bias = false; // correct MPF for L1L2L3-L1 (instead of L1L2L3-RC)
bool _paper = true;
//bool _g_dcsonly = false;
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

//const int njesFit = 1; // scale only
//const int njesFit = 2; // scale+HB+fixedPileUpPt
//const int njesFit = 3; // scale+HB+PileUpPt
bool hbpt = false;//true; // pT dependent HB scale instead of PileUpPt
const int njesFit = 4; // scale+HBxPt+PileUpPt
bool erffit = false;//true; // error function fit
bool stepfit = true;
bool basicfit = false;
const double ptminJesFit = 30;//300;
TF1 *fhb(0), *fl1(0); double _etamin(0);
Double_t jesFit(Double_t *x, Double_t *p) {

  double pt = *x;

  // Initialize SinglePionHCAL and PileUpPt shapes
  if (!fhb) fhb = new TF1("fhb","max(0.,[0]+[1]*pow(x,[2]))",10,3500);
  fhb->SetParameters(1.03091e+00, -5.11540e-02, -1.54227e-01); // SPRH
  if (!fl1) fl1 = new TF1("fl1","1+([0]+[1]*log(x))/x",30,2000);
  //fl1->SetParameters(-2.36997, 0.413917); // Run I
  fl1->SetParameters(-0.7866069, 0.1330742); // Run II
 
  // Just fitting a constant scale factor for RelativePt uncertainty
  if (njesFit==1) {
    return p[0];
  }

  // Directly using fitted SPR HCAL shape (from JECUncertainty.cpp)
  if (njesFit==2) {

    // p[0]: overall scale shift, p[1]: HCAL shift in % (full band +3%)
    double ptref = 208; // pT that minimizes correlation in p[0] and p[1]
    double ptx = (stepfit ? max(150.,min(340.,pt)) : pt); // 147.4/87 V2M2

    return (p[0] + p[1]/3.*100*(fhb->Eval(ptx)-fhb->Eval(ptref)));
  } // njesFit==2


  // HB scale + const + residual log-pT offset
  if (njesFit==3 && !hbpt) {

    // p[0]: overall scale shift, p[1]: HCAL shift in % (full band +3%)
    // p[2]: fraction of PileUpPtBB uncertainty
    //double ptref = 225.; // GT
    double ptref = 208.; // newL1
    //return (p[0] + p[1]/3.*100*(fhb->Eval(pt)-fhb->Eval(ptref))
    //    + p[2]*(fl1->Eval(pt)-fl1->Eval(ptref)));
    double ptx = max(220.,min(300.,pt));
    return (p[0] + p[1]/3.*100*(fhb->Eval(ptx)-fhb->Eval(ptref))
	    + p[2]*(fl1->Eval(ptx)-fl1->Eval(ptref)));
  } // njesFit==3

  // pT-dependent HB scale + const
  if (njesFit==3 && hbpt) {

    // p[0]: overall scale shift, p[1]: HCAL shift in % (full band +3%)
    // p[3]: HCAL scale log-linear pT dependence
    double ptref = 208.;
    double hbscale = p[1] + p[2]*log(pt/ptref);
    return (p[0] + hbscale/3.*100*(fhb->Eval(pt)-fhb->Eval(ptref)));
  } // njesFit==3

  // Less physical error function fit
  if (njesFit==4 && erffit) {
    
    //return (p[0] + p[1]*TMath::Erf((log(pt/p[2]) / log(p[3]/p[2]))));
    return (p[0] + p[1]*TMath::Erf((log(pt) - log(p[2])) / p[3]));
    //return (p[0] + p[1]*TMath::Erf((pt - p[2]) / p[3]));
  } // njesFit==4

  // Directly using fitted SPR HCAL shape (from JECUncertainty.cpp)
  if (njesFit==4 && stepfit) {

    // p[0]: overall scale shift, p[1]: HCAL shift in % (full band +3%)
    double ptref = 208; // pT that minimizes correlation in p[0] and p[1]
    //double x0(150), x1(340);
    double x0(150), x1(400);
    double r(1), k(100./3.);
    double f = fhb->Eval(pt);
    double f0 = fhb->Eval(x0);
    double f1 = fhb->Eval(x1);
    double fref = fhb->Eval(ptref);
//     if (pt<x0) {
//       r = (p[0] + (p[1]-p[2])*k*(f0-fref)) + p[2]*k*(f - fref);
//     }
//     if (pt>=x0 && pt<x1) {
//       r = p[0] + p[1]*k*(f - fref);
//     }
//     if (pt>=x1) {
//       r = (p[0] + (p[1]-p[3])*k*(f1-fref)) + p[3]*k*(f - fref);
//     }
    double r3p = p[0] + p[1]*k*(f-fref);
    double r1p = p[2] + p[3]*k*(f-fref);
    double w = (pt-x0)/(x1-x0); // linear version
    //double w = (log(pt)-log(x0))/(log(x1)-log(x0)); // log-linear version
    if (pt<x0) r = r3p;
    if (pt>=x0 && pt<x1) r = r3p*(1-w) + r1p*w;
    if (pt>=x1) r = r1p;

    return r;
  } // njesFit==4

  // pT-dependent HB scale + const + residual log-pT offset
  if (njesFit==4 && basicfit) {

    // p[0]: overall scale shift, p[1]: HCAL shift in % (full band +3%)
    // p[2]: fraction of PileUpPtBB uncertainty
    // p[3]: HCAL scale log-linear pT dependence
    double ptref = 208.;
    double hbscale = p[1] + p[3]*log(pt/ptref);
    return (p[0] + hbscale/3.*100*(fhb->Eval(pt)-fhb->Eval(ptref))
	    + p[2]*(fl1->Eval(pt)-fl1->Eval(ptref)));
  } // njesFit==4

}

void globalFitL3Res(double etamin = 0, double etamax = 1.3) {
  
  _etamin = etamin;
  
  TDirectory *curdir = gDirectory;
  setTDRStyle();
  writeExtraText = false; // for JEC paper CWR

  TFile *f = new TFile("rootfiles/jecdata.root","READ");
  assert(f && !f->IsZombie());

  // Settings
  const double alpha = 0.30;
  const char *ct = "ratio";
  //
  const int nmethods = 2;
  const char* methods[nmethods] = {"ptchs","mpfchs1"};
  _nmethods = nmethods; // for multijets in global fit

  // Global fit with multijets
  const int nsamples = 4;
  const int nsample0 = 1; // first Z/gamma+jet sample
  const char* samples[4] = {(etamin==0 ? "multijet" : "dijet"),
			    "gamjet", "zeejet", "zmmjet"};
  const int igj = 0;
  const int izee = 1;
  const int izmm = 2;


  /*
  // Global fit without multijets
  const int nsamples = 3;
  const int nsample0 = 0; // first Z/gamma+jet sample
  const char* samples[3] = {"gamjet", "zeejet", "zmmjet"};
  const int igj = 0;
  const int izee = 1;
  const int izmm = 2;
  */
  // old default
  /*
  const int nsamples = 3;
  const int nsample0 = 1; // first Z/gamma+jet sample
  const char* samples[3] = {(etamin==0 ? "multijet" : "dijet"),
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
  //
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

      string s = Form("%s_%s_a%02.0f",cm,cs,alpha*100);
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
	    (string(cs)=="dijet" && g->GetX()[i]<70.) ||
	    (string(cs)=="gamjet" && g->GetX()[i]>600. && etamin!=0))
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
    string s = "ptf30_multijet_a30";
    TGraphErrors *g = (TGraphErrors*)d->Get(s.c_str());
    if (!g) cout << "Graph "<<s<<" not found!" << endl << flush;
    assert(g);
    
    g->SetLineWidth(2);
    gfs.push_back((TGraphErrors*)g->Clone());
    g->SetMarkerStyle(kOpenTriangleDown);
    gfs2.push_back((TGraphErrors*)g->Clone()); // for MJB
    
    // Fractions for MPF (pT>30 GeV)
    s = "ptf10_multijet_a30";
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
  TH1D *hjes0 = (TH1D*)f->Get("ratio/eta00-13/hjes"); assert(hjes0);
  hjes0->SetName("hjes0");
  TH1D *herr0 = (TH1D*)f->Get("ratio/eta00-13/herr"); assert(herr0);
  herr0->SetName("herr0");

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
	double kfsr = (dofsr ? 1./(1+alpha*h->GetBinContent(h->FindBin(pt))):1);
	double l1 = (dol1bias && string(cm)=="mpfchs1" && string(cs)=="gamjet" ?
		     1./hl1->GetBinContent(hl1->FindBin(pt)) : 1);
	// 0.98 was on until Aug 12, 10:12am :(
	double scale = 1.00;//0.98; // correct out previous L3Res
	if (etamin!=0) {
	  int ijes = min(max(1,hjes->FindBin(pt)),hjes->GetNbinsX());
	  scale = hjes->GetBinContent(ijes);
	}
	g->SetPoint(i, pt, scale*l1*r*kfsr);
	g2->SetPoint(i, pt, scale*l1*r*kfsr);
	g3->SetPoint(i, pt, scale*l1*r);
	
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
	h->Scale(alpha); // set which nominal points this applies to
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
  for (unsigned int i = 0; i != _vdt->size(); ++i) {
      
    TH1D *h = hs[i]; assert(h);
    TH1D *h2 = (TH1D*)h->Clone(Form("bm%d_scale_%d",1<<i,i));

    string s = samples[i%nsamples];
    unsigned int imethod = i/nsamples;

    double escale(0);
    const int n0 = nsample0;
    const int n1 = nsamples+nsample0;
    if (s=="multijet" && imethod==0) {
      escale = 1.0; // no constraint on absolute scale
      h2->SetName(Form("bm%d_scale_%d",(1<<(n0-1) | (1<<(n1-1))), i));
    }
    if (s=="dijet" && imethod==0) {
      escale = 0.005; // for JER bias and whatnot
      h2->SetName(Form("bm%d_scale_%d",(1<<(n0-1) | (1<<(n1-1))), i));
    }
    if (s=="gamjet" && imethod==0) {
      //escale = 0.005; // photon (ECAL) scale without regression
      //escale = 0.002; // photon (ECAL) scale with regression (Run I)
      escale = 0.010; // Run II guess
      // Use same source for both MPF and pT balance
      h2->SetName(Form("bm%d_scale_%d",(1<<(n0+igj) | (1<<(n1+igj))), i));
      is = hs.size();
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
      escale = 0.005; // 76X guess
      // Use same source for both MPF and pT balance
      h2->SetName(Form("bm%d_scale_%d",(1<<(n0+izee) | (1<<(n1+izee))), i));
    }
    if (s=="zmmjet" && imethod==0) {
      // Use same source for both MPF and pT balance
      //escale = 0.002; // muon (tracker) scale (Run I)
      //escale = 0.010; // Run II guess
      escale = 0.005; // 76X guess
      h2->SetName(Form("bm%d_scale_%d",(1<<(n0+izmm) | (1<<(n1+izmm))), i));
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

    TH1D *h = hs[i]; assert(h);
    TH1D *h2 = (TH1D*)h->Clone(Form("bm%d_mpfscale_%d",1<<i,i));

    string s = samples[i%nsamples];
    string m = methods[i/nsamples];

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
	  //				     1<<(n1+0)|1<<(n1+1)|1<<(n1+2)),i));
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

  // Uncertainty sources for multijets
  if (string(samples[0])=="multijet") {

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
  lumi_13TeV = "Run2015D - Golden - 2.1 fb^{-1}";
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
  texlabel["zeejet"] = "Zee+jet";
  texlabel["zmmjet"] = "Z#mu#mu+jet";
  TLegend *legm = tdrLeg(0.66,0.55,0.96,0.90);
  legm->SetHeader("MPF");
  for (int i = 0; i != nsamples; ++i)
    legm->AddEntry(gs[i+nsamples],texlabel[samples[i]],i==0 ? "" : "PL");

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  if (etamin==0) {
    if (dofsr) tex->DrawLatex(0.20,0.79,"|#eta|<1.3, #alpha<0.3#rightarrow0");
    else       tex->DrawLatex(0.20,0.79,"|#eta|<1.3, #alpha<0.3");
  }
  else {
    assert(dofsr);
    tex->DrawLatex(0.20,0.79,Form("%1.1f#leq|#eta|<%1.1f",etamin,etamax));
  }
  tex->DrawLatex(0.20, 0.22, "Before global fit");

  herr_ref->SetLineWidth(2);
  herr_ref->SetLineColor(kYellow+3);
  herr_ref->SetLineStyle(kDashed);
  herr_ref->DrawClone("SAME E3");
  (new TGraph(herr_ref))->DrawClone("SAMEL");
  legp->AddEntry(herr_ref," ","");
  legm->AddEntry(herr_ref,"JES unc.","FL");

  hrun1->SetLineWidth(2);
  hrun1->SetLineColor(kCyan+3);
  hrun1->SetLineStyle(kDashed);
  hrun1->DrawClone("SAME E3");
  (new TGraph(hrun1))->DrawClone("SAMEL");
  legp->AddEntry(hrun1," ","");
  legm->AddEntry(hrun1,"Run 1","FL");


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

  if (etamin==0) tex->DrawLatex(0.20,0.79,"|#eta|<1.3, #alpha<0.3");
  else tex->DrawLatex(0.20,0.79,Form("%1.1f#leq|#eta|<%1.1f",etamin,etamax));

  tex->DrawLatex(0.20, 0.22, "Before FSR correction");

  herr_ref->DrawClone("SAME E3");
  (new TGraph(herr_ref))->DrawClone("SAMEL");

  hrun1->DrawClone("SAME E3");
  (new TGraph(hrun1))->DrawClone("SAMEL");

  herr_ref->SetFillStyle(kNone);
  herr_ref->DrawClone("SAME E3");
  (new TGraph(herr_ref))->DrawClone("SAMEL");
  herr_ref->SetFillStyle(1001);

  hrun1->SetFillStyle(kNone);
  hrun1->DrawClone("SAME E3");
  (new TGraph(hrun1))->DrawClone("SAMEL");
  hrun1->SetFillStyle(1001);

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
    jesfit->SetParameter(0, 0.95);
  }
  if (njesFit==2) {
    //jesfit->SetParameters(0.9828, -0.00608);
    jesfit->SetParameters(0.9350, -0.0922);
  }
  if (njesFit==3 && !hbpt) {
    jesfit->SetParameters(0.99, -0.035, -0.09);
  }
  if (njesFit==3 && hbpt) {
    jesfit->SetParameters(0.98, -0.035, 0.01);
  }
  if (njesFit==4 && basicfit) {
    jesfit->SetParameters(0.98, -0.035, -0.09, 0.01);
  }
  if (njesFit==4 && stepfit) {
    //jesfit->SetParameters( 0.9708, 0.4081, 250., 400.);
    //jesfit->SetParLimits(2,100,300);
    //jesfit->SetParLimits(3,350,600);
    //jesfit->FixParameter(2, 250.);
    //jesfit->FixParameter(3, 400.);
    //jesfit->SetParameters( 0.9708, 0.4081, 0.001, 0.001);
    jesfit->SetParameters( 0.985, 0.03, 1.03, -0.09 );
  }
  if (njesFit==4 && erffit) {
    jesfit->SetParameters(0.985, 0.02, 200, 0.5); // log(pt)
    //jesfit->SetParameters(0.975, 0.02, 250, 50);
    //jesfit->SetParameters(0.9785, +0.0199, 100, 100);
    //jesfit->FixParameter(3,100);
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
    cout << "  " << vchi2_data[i] << " / " << vndata[i]
	 << "  ("<<gs2[i]->GetName()<<")"<< endl;
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
  for (int i = np; i != Npar; ++i) {
    if ((i-np)%(2*nsamples)==0) cout << (*_vsrc)[i-np]->GetName() << endl;
    if (tmp_par[i]==0&&fabs(tmp_err[i]-1)<1e-2)
      cout << Form("%2d:   ----------,    ",i-np+1);
    else
      cout << Form("%2d: %5.2f+/-%4.2f,   ",i-np+1,tmp_par[i],tmp_err[i]);
    if ((i-np)%nsamples==nsamples-1) cout << endl;
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

  if (etamin==0) tex->DrawLatex(0.20,0.79,"|#eta|<1.3, #alpha<0.3#rightarrow0");
  else tex->DrawLatex(0.20,0.79,Form("%1.1f#leq|#eta|<%1.1f",etamin,etamax));

  //tex->DrawLatex(0.20,0.18,Form("(data %1.1f / %d, sources %1.1f / %d)",
  //			chi2_data,Nk,chi2_src,nsrc_true));

  tex->DrawLatex(0.20,0.22,"After global fit");
  tex->DrawLatex(0.20,0.17,Form("#chi^{2} / NDF = %1.1f / %d",
				chi2_gbl, Nk-np));

  is += np;

  // Fitted lepton/photon scales too much detailed for the paper
  //tex->DrawLatex(0.20,0.73,Form("#gamma #times (%1.3f#pm%1.3f)",
  //			1+0.01*tmp_par[is],
  //			0.01*sqrt(emat[is][is])));
  //tex->DrawLatex(0.20,0.73,Form("e #times (%1.3f#pm%1.3f)",
  //			1+0.005*tmp_par[is+1],
  //			0.005*sqrt(emat[is+1][is+1])));
  //tex->DrawLatex(0.20,0.67,Form("#mu #times (%1.3f#pm%1.3f)",
  //			1+0.01*tmp_par[is+2],
  //			0.01*sqrt(emat[is+2][is+2])));

  herr_ref->DrawClone("SAME E3");
  (new TGraph(herr_ref))->DrawClone("SAMEL");

  hrun1->DrawClone("SAME E3");
  (new TGraph(hrun1))->DrawClone("SAMEL");

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

  herr_ref->SetFillStyle(kNone);
  herr_ref->DrawClone("SAME E3");
  (new TGraph(herr_ref))->DrawClone("SAMEL");
  herr_ref->SetFillStyle(1001);

  hrun1->SetFillStyle(kNone);
  hrun1->DrawClone("SAME E3");
  (new TGraph(hrun1))->DrawClone("SAMEL");
  hrun1->SetFillStyle(1001);

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

  if (etamin==0) {
    if (dofsr) {
      c0b->SaveAs("pdf/globalFitL3res_raw.pdf");
      c0b->SaveAs("pdfC/globalFitL3res_raw.C");
      c0->SaveAs("pdf/globalFitL3res_orig.pdf");
      c0->SaveAs("pdfC/globalFitL3res_orig.C");
      c1->SaveAs("pdf/globalFitL3res_shifted.pdf");
      c1->SaveAs("pdfC/globalFitL3res_shifted.C");
    }
  }
  else {
    c0b->SaveAs(Form("pdf/globalFitL3res_raw_eta%1.0f-%1.0f.pdf",
		     10*etamin,10*etamax));
    c0->SaveAs(Form("pdf/globalFitL3res_orig_eta%1.0f-%1.0f.pdf",
		    10*etamin,10*etamax));
    c1->SaveAs(Form("pdf/globalFitL3res_shifted_eta%1.0f-%1.0f.pdf",
		    10*etamin,10*etamax));
  }

  // Draw nuisance parameter distribution
  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  hsrc->SetFillColor(kBlue-10);
  hsrc->SetFillStyle(1001);
  hsrc->Draw();
  gStyle->SetOptStat(111111);

  if (etamin==0) c2->SaveAs("pdf/globalFitL3res_hsrc.pdf");
  else {
    c2->SaveAs(Form("pdf/globalFitL3res_hsrc_eta%1.0f-%1.0f.pdf",
		    10*etamin, 10*etamax));
  }
  gStyle->SetOptStat(0);

  // Draw FSR corrections redone

  int colors[3] = {kBlue, kGreen+2,kRed};
  //int colors[nsamples] = {kBlue, kGreen+2,kRed};
  //int colors[nsamples] = {kBlue, kRed};
  //int colors[nsamples] = {kBlue};
  for (int imethod = 0; imethod != nmethods; ++imethod) {

    h->SetMinimum(-0.06);//imethod==0 ? -0.20 : -0.20);//-0.04);
    h->SetMaximum(+0.09);//imethod==0 ? +0.30 : +0.30);//+0.06);
    h->SetYTitle("k_{FSR} = dR/d#alpha (ratio)");
    h = (TH1D*)h->Clone(Form("h3_%d",imethod));

    TCanvas *c3 = tdrCanvas(Form("c3_%d",imethod),h,4,11,true);
    c3->SetLogx();

    TLegend *leg1 = tdrLeg(0.58,0.65,0.88,0.90);
    leg1->SetHeader("In");
    TLegend *leg2 = tdrLeg(0.65,0.65,0.95,0.90);
    leg2->SetHeader("Out");

    tex->DrawLatex(0.55,0.20,
		   imethod==0 ? "p_{T} balance method" : "MPF method");

    for (int isample = nsample0; isample != nsamples; ++isample) {

      int ibm = isample + nsamples*imethod;
      TH1D *hk = hks[ibm]; assert(hk);

      double minx = (isample==nsample0 ? 40  : 30);
      double maxx = (isample==nsample0 ? 900 : 600);//800 : 300);
      hk->GetXaxis()->SetRangeUser(minx,maxx);

      //if (isample==nsample0) {
      //hk->GetXaxis()->SetRangeUser(40,800);
      //}
      //if (isample!=nsample0) {
      //hk->GetXaxis()->SetRangeUser(30,300);
      //}

      hk->SetLineColor(colors[isample-nsample0]-10);
      hk->SetFillColor(colors[isample-nsample0]-10);
      hk->SetFillStyle(kNone);
      hk->SetMarkerSize(0.);
      hk->DrawClone("SAME E3");
      hk->SetFillStyle(3002);
      hk->DrawClone("SAME E3");

      //(new TGraph(hk))->DrawClone("SAMEL");
      TGraph *gk = new TGraph(hk);
      for (int i = gk->GetN()-1; i != -1; --i) {
	//if (gk->GetX()[i] < hk->GetXaxis()->GetXmin())
	if (gk->GetX()[i] < minx || gk->GetX()[i] > maxx)
	  gk->RemovePoint(i);
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
	hke->SetBinContent(ipt, (alpha*hk->GetBinContent(ipt) - yeig)/alpha);
	hke->SetBinError(ipt, sqrt(err2));
      } // for ipt
      
      //if (isample==nsample0) hke->GetXaxis()->SetRangeUser(40,800);
      //if (isample!=nsample0) hke->GetXaxis()->SetRangeUser(30,300);
      hke->GetXaxis()->SetRangeUser(minx,maxx);
      hke->SetFillStyle(1001);
      hke->SetLineColor(hk->GetLineColor()+10);
      hke->DrawClone("SAME E3");
      //(new TGraph(hke))->DrawClone("SAMEL");
      TGraph *gke = new TGraph(hke);
      for (int i = gke->GetN()-1; i != -1; --i) {
        if (gke->GetX()[i] < minx || gke->GetX()[i] > maxx)
          gke->RemovePoint(i);
      }
      gke->Draw("SAMEL");

      leg2->AddEntry(hke, texlabel[samples[isample]], "FL");
    } // for isample

    if (etamin==0) {
      c3->SaveAs(Form("pdf/globalFitL3res_%s_kfsr.pdf",methods[imethod]));
      c3->SaveAs(Form("pdfC/globalFitL3res_%s_kfsr.C", methods[imethod]));
    }
    else {
      c3->SaveAs(Form("pdf/globalFitL3res_%s_kfsr_eta%1.0f-%1.0f.pdf",
		      methods[imethod],10*etamin,10*etamax));
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
