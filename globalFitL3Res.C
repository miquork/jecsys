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

#include <list>
#include <Minuit2/Minuit2Minimizer.h>
#include <Math/Functor.h>

using namespace std;

// Global fit settings
double alpha = 0.30; // reference alpha value for zjet, gamjet
double ptmj = 30.; // reference pTmin value for multijet
bool dofsr = true; // correct for FSR central value
bool dofsrcombo = true; // use BCDEF for FSR, not individual IOVs
double alphaeffmax = 0.45; // max alphaeff for FSR corrections
double ptreco_gjet = 15.; // min jet pT when evaluating alphamax for gamma+jet
double ptreco_zjet = 5.; // same for Z+jet
bool dol1bias = false; // correct MPF for L1L2L3-L1 (instead of L1L2L3-RC)
// see also ScalingForL2OrL3Fit below (left there since long explanation)

// Plotting settings
bool _paper = true; // switch of plotting fit parameters
bool _useZoom = true; // also affects the kind of uncertainty band plotted: useZoom=true comes by default with AbsoluteScale+TotalNoFlavorNoTime; false--> Run1 and reference AbsoluteScale
bool plotMultijetDown = true; // plot gray downward points for multijets
double ptmaxMultijetDown = 300; // max pT for downward multijet points
double shiftPtBal = 0.975; // move x-axis for pTbal, 0 or 1 for none
bool storeErrorMatrix = true; // error matrix for globalFitL3Pulls.C
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
const int njesFit = 2;
double fixOff = -0.6035; // p3 185.9/62
const double ptminJesFit = 30;

TF1 *fhb(0), *fl1(0), *fl1b(0), *fl1mc(0), *ftr(0), *feg(0); double _etamin(0);
Double_t jesFit(Double_t *x, Double_t *p) {
  
  double pt = *x;

  // Initialize SinglePionHCAL shape
  if (!fhb) fhb = new TF1("fhb","max(0.,[0]+[1]*pow(x,[2]))",10,3500);
  fhb->SetParameters(1.03091e+00, -5.11540e-02, -1.54227e-01); // SPRH
 
  // Initialize L1FastJe_Simple-L1RC difference
  // Values from fitting ratio/eta00-13/hl1bias (JEC set in reprocess.C)
  if (!fl1) fl1 = new TF1("fl1","1-([0]+[1]*log(x)+[2]*pow(log(x),2))/x",10,3500);
  //fl1->SetParameters(1.72396, 0, 0); // UL17 V1S hl1bias (p1=p2=0) (UL17_V1)
  fl1->SetParameters(0.350077, 0.553560, -0.0527681); // BCDEF hl1rcos (RC vs Simple)

  if (!fl1b) fl1b = new TF1("fl1b","1-([0]+[1]*log(x)+[2]*pow(log(x),2))/x",10,3500);
  fl1b->SetParameters(-2.62921, 2.34488, -0.429717); // BCDEF hl1cos

  double ptref = 208; // pT that minimizes correlation in p[0] and p[1]
  // (or used to minimize in Run I: for Run II, 208.*sqrt(13/8)=265)

  // Constant scale factor (EM scale) only
  if (njesFit==1) {

    return p[0];
  } // njesFit==2

  // Directly using fitted SPR HCAL shape (from JECUncertainty.cpp)
  if (njesFit==2) {

    // p[0]: overall scale shift, p[1]: HCAL shift in % (full band +3%)
    return (p[0] + p[1]/3.*100*(fhb->Eval(pt)-fhb->Eval(ptref))
	    + fixOff*(fl1->Eval(pt)-fl1->Eval(ptref)));
  } // njesFit==2

  if (njesFit==3) {

    // p[0]: overall scale shift, p[1]: HCAL shift in % (full band +3%)
    // p[2]: L1FastJet-L1RC difference
    return (p[0] + p[1]/3.*100*(fhb->Eval(pt)-fhb->Eval(ptref))
	    + p[2]*(fl1->Eval(pt)-fl1->Eval(ptref)));
  } // njesFit==3

  if (njesFit==4) {

    // p[0]: overall scale shift, p[1]: HCAL shift in % (full band +3%)
    // p[2]: L1FastJet-L1RC difference
    return (p[0] + p[1]/3.*100*(fhb->Eval(pt)-fhb->Eval(ptref))
	    + p[2]*(fl1->Eval(pt)-fl1->Eval(ptref))
	    + p[3]*(fl1b->Eval(pt)-fl1b->Eval(ptref)));
  } // njesFit==4

  assert(0);
  exit(0);
}

void globalFitL3Res(double etamin = 0, double etamax = 1.3,
		    string epoch="", string selectSample="Standard_MJDJ_gam_zee_zmm", string selectMethods="PtBalMPF") {
  if(verboseGF)cout << Form("Running globalFitL3Res(etamin=%02.2f,etamax=%02.2f,epoch=%s, selectSample=%s, selectMethods=%s",etamin,etamax,epoch.c_str(),selectSample.c_str(),selectMethods.c_str()) << endl << flush;
  _etamin = etamin;
  const char *cep = epoch.c_str();

  TDirectory *curdir = gDirectory;
  setTDRStyle();

  TFile *f = new TFile(Form("rootfiles/jecdata%s.root",cep),"READ");
  assert(f && !f->IsZombie());

  TFile *ffsr = new TFile(Form("rootfiles/jecdata%s.root","BCDEF"),"READ");
  assert((ffsr && !ffsr->IsZombie()) || !dofsrcombo);
  if (dofsrcombo) cout << "Use FSR from jecdataBCDEF.root for all IOVs\n";

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

  _nsamples = nsamples; // for multijets in global fit
  
  string sbin = Form("eta%02.0f-%02.0f",10*etamin,10*etamax);
  const char* bin = sbin.c_str();

  TDirectory *dfsr(0);
  if (dofsrcombo) {
    assert(ffsr->cd(ct));
    TDirectory *dfsrin1 = ffsr->GetDirectory(ct); assert(dfsrin1);
    assert(dfsrin1->cd(bin));
    dfsr = dfsrin1->GetDirectory(bin); assert(dfsr);
  } // usefsrcombo

  assert(f->cd(ct));
  TDirectory *din1 = f->GetDirectory(ct); assert(din1);
  assert(din1->cd(bin));
  TDirectory *d = din1->GetDirectory(bin); assert(d);
  d->cd();


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
      TGraphErrors *g = (TGraphErrors*)d->Get(s.c_str());
      if (!g) cout << "Graph "<<s<<" not found!" << endl << flush;
      assert(g);
      
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
  
  // Load FSR (JESref) corrections, and correct data
  vector<TH1D*> hks(nsamples*nmethods,0);
  for (unsigned int imethod = 0; imethod != nmethods; ++imethod) {
    for (int isample = 0; isample != nsamples; ++isample) {

      string sm = methods[imethod];
      const char *cm = sm.c_str();
      string ss = samples[isample];
      const char *cs = ss.c_str();

      // For correcting FSR bias
      string s = Form("fsr/hkfsr_%s_%s",cm,cs);
      TH1D *hfsr = (TH1D*)d->Get(s.c_str());
      if (dofsrcombo) hfsr = (TH1D*)dfsr->Get(s.c_str());
      if (!hfsr) cout << "Histo "<<s<<" not found!" << endl << flush;
      assert(hfsr);

      // For correcting L1L2L3-L1 to L1L2L3-RC
      TH1D *hl1 = (TH1D*)d->Get("hl1bias"); assert(hl1);
      
      int ibm = isample + nsamples*imethod;
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
	if (ss=="zeejet"||ss=="zmmjet"||ss=="zlljet") ptreco = ptreco_zjet;
	if (ss=="gamjet") ptreco = ptreco_gjet;
	double aeff = min(max(alpha, ptreco/pt),alphaeffmax);
	if (ss=="multijet") aeff = alpha;
	double kfsr(1);
	if (dofsr) kfsr = 1./(1+aeff*hfsr->GetBinContent(hfsr->FindBin(pt)));
	double l1(1);
	if (dol1bias && sm=="mpfchs1" && ss=="gamjet")
	  l1 = 1./hl1->GetBinContent(hl1->FindBin(pt));

	double scale = 1.00; // correct out previous L3Res

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
  const int neigmax = 4;//3;
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
	TH1D *h = (TH1D*)d->Get(s.c_str());
	if (dofsrcombo) h = (TH1D*)dfsr->Get(s.c_str());
	if (ieig>0 && !h) {
	  h = hs[ibm]; assert(h); // use src0
	  h = (TH1D*)h->Clone(Form("bm%d_inactive_hkfsr_%s_%s_eig%d",
	  			   (1<<ibm),cm,cs,ieig));
	  h->Reset();
	  // Store also empty ones for now to have FSR plots work
	  hs.push_back(h);
	  continue;
	}
	if (!h) cout << "Histo "<<s<<" not found!" << endl << flush;
	assert(h);

	for (int i = 1; i != h->GetNbinsX()+1; ++i) {
	  double pt = h->GetBinCenter(i);
	  double ptreco(0);
	  if (ss=="zeejet"||ss=="zmmjet"||ss=="zlljet") ptreco = ptreco_zjet;
	  if (ss=="gamjet") ptreco = ptreco_gjet;
	  double aeff = min(max(alpha, ptreco/pt),alphaeffmax);
	  if (ss=="multijet") aeff = alpha;
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
      h2->SetName(Form("bm%d_scale_%03.0f_gamjet_%d",
		       (1<<(n0+igj) | (1<<(n1+igj))), escale*10000., i));
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
      h2->SetName(Form("bm%d_scale_%03.0f_zlljet_%d",
		       (1<<(n0+izll) | (1<<(n1+izll))), escale*10000, i));
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

    TH1D *h = hs[i]; assert(h);
    TH1D *h2 = (TH1D*)h->Clone(Form("bm%d_%s_%02.0f_%s_%s_%d",
				    1<<i,
				    escale!=0 ? "mpfscale" : "inactive",
				    escale*1000.,cm,cs,i));
    
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
  if (!dofsrcombo) {
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
  } // !dofsrcombo

  // UL17-V2: add l2cos x <rho>/<rho>_ref as low pT systematic
  if (true) {

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

    const int nsrcmj = 1;//3;
    for (unsigned int imethod = 0; imethod != nmethods; ++imethod) {
      for (int isrc = 0; isrc != nsrcmj; ++isrc) {

	string ss = "multijet";
	const char *cs = ss.c_str();
	string sm = methods[imethod];
	const char *cm = sm.c_str();
	
	// Don't remember why these sources were stored under fsr/, but ok
	string s;
	if (isrc==0) s = "sys/jer_multijet";
	TH1D *h = (TH1D*)d->Get(s.c_str());
	if (dofsrcombo) h = (TH1D*)dfsr->Get(s.c_str());
	if (!h) cout << "Histo "<<s<<" not found!" << endl << flush;
	assert(h);
	
	// JER uncertainty applies for both MPF and MJB
	if (isrc==0 && sm=="ptchs") { // jer
	  h->SetName(Form("bm%d_multijet_jer_src%d_%d",
			  (1<<0 | 1<<nsamples),isrc,hs.size()+1));
	  hs.push_back(h);
	}
      } // for imethod
    } // for ieig
  } // multijet syst.

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
			(1<<(n0+igj) | 1<<(n1+igj)),ieig));
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
  texlabel["dijet"] = "Dijet";
  texlabel["gamjet"] = "#gamma+jet";
  texlabel["zeejet"] = "Zee+jet";
  texlabel["zmmjet"] = "Z#mu#mu+jet";
  texlabel["zlljet"] = "Z+jet";

  TLegend *legm = tdrLeg(0.62,_useZoom ? 0.65 : 0.60,0.90,0.90);
  if( (nmethods==1&&strcmp(methods[0],"mpfchs1")  ==0) || (nmethods==2 && strcmp(methods[1],"mpfchs1")  ==0 )){
    legm->SetHeader("MPF");
    for (int i = 0; i != nsamples; ++i)
      legm->AddEntry(gs1[i+(strcmp(methods[0],"mpfchs1")==0 ? 0 : nsamples)],texlabel[samples[i]],i==0 ? "" : "PL");
  }
  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  if (isl3) {
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
    
    g2->DrawClone("SAMEPz");
  }


  if (!_useZoom) {
    legp->AddEntry(hrun1," ","");
    legm->AddEntry(hrun1,"Run I","FL");
    legp->AddEntry(herr_ref," ","");
    if (!_paper) legm->AddEntry(herr_ref,"UL17V2M3 IOV","FL");
    if ( _paper) legm->AddEntry(herr_ref,"Run II","FL");
  }
  else {
    legp->AddEntry(herr," ","FL");
    if (!_paper) legm->AddEntry(herr_ref,"UL17V2M3","FL");
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
  if (njesFit==3)  jesfit->SetParameters(0.981, 0.044, -0.519);
  if (njesFit==4)  jesfit->SetParameters(0.981, 0.044, -0.519,0 );

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
      h2emat->GetYaxis()->SetBinLabel(i+1,Form("%s (%d)",cp,i));
      h2cov->GetYaxis()->SetBinLabel(i+1,Form("%s (%d)",cp,i));
      h1par->GetXaxis()->SetBinLabel(i+1,Form("%s (%d)",cp,i));
    } // for i
    emat.Write("emat");
    vpar.Write("vpar");
    verr.Write("verr");
    h2emat->Write();
    h2cov->Write();
    h1par->Write();
    femat->Close();
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
      if (dofsrcombo) h = (TH1D*)dfsr->Get(s.c_str());
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

  if (!_paper) {
 
    //tex->SetTextColor(kWhite); // hide from view
    tex->DrawLatex(0.25,0.67,Form("Npar=%d",_jesFit->GetNpar()));
    tex->DrawLatex(0.25,0.64,Form("p0=%1.4f #pm %1.4f",
				  _jesFit->GetParameter(0),
				  sqrt(emat[0][0])));
    if (njesFit>=2)
      tex->DrawLatex(0.25,0.61,Form("p1=%1.4f #pm %1.4f",
				    _jesFit->GetParameter(1),
				    sqrt(emat[1][1])));
    if (njesFit==2 && fixOff)
      tex->DrawLatex(0.25,0.58,Form("p2=%1.4f (fixOff)",fixOff));
    if (njesFit>=3)
      tex->DrawLatex(0.25,0.58,Form("p2=%1.4f #pm %1.4f",
				  _jesFit->GetParameter(2),
				    sqrt(emat[2][2])));
    if (njesFit>=4)
      tex->DrawLatex(0.25,0.54,Form("p3=%1.4f #pm %1.4f",
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

  int colorsdark[4]  = {kBlack,kBlue, kGreen+2,kRed};
  int colorslight[4] = {kGray, kBlue-10, kGreen+2-10,kRed-10};
  for (unsigned int imethod = 0; imethod != nmethods; ++imethod) {

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

      int ibm = isample + nsamples*imethod;
      TH1D *hk = hks[ibm]; assert(hk);

      string ss = samples[isample];
      const char *cs = ss.c_str();

      double minx(30);
      double maxx(1500);
      if (ss=="zeejet"||ss=="zmmjet"||ss=="zlljet") {  maxx = 700; }
      if (ss=="gamjet") { minx=230; maxx = 1500; }
      if (ss=="multijet") { minx=114; maxx = 2116; }//2640; }
      hk->GetXaxis()->SetRangeUser(minx,maxx);

      hk->SetLineColor(colorslight[isample]);
      hk->SetFillColor(colorslight[isample]);
      hk->SetFillStyle(kNone);
      hk->SetMarkerSize(0.);
      hk->DrawClone("SAME E5");
      hk->SetFillStyle(3002);
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
	if (ss=="zeejet"||ss=="zmmjet"||ss=="zlljet") ptreco = ptreco_zjet;
	if (ss=="gamjet") ptreco = ptreco_gjet;
	double aeff = min(max(alpha, ptreco/pt),alphaeffmax);
	if (ss=="multijet") aeff = alpha;
	hke->SetBinContent(ipt, (aeff*hk->GetBinContent(ipt) - yeig)/aeff);
	hke->SetBinError(ipt, sqrt(err2));
      } // for ipt
      
      hke->GetXaxis()->SetRangeUser(minx,maxx);
      hke->SetFillStyle(1001);
      hke->SetLineColor(colorsdark[isample]);
      hke->DrawClone("SAME E5");

      TGraph *gke = new TGraph(hke);
      for (int i = gke->GetN()-1; i != -1; --i) {
        if (gke->GetX()[i] < minx || gke->GetX()[i] > maxx)
          gke->RemovePoint(i);
        if (gke->GetY()[i]==0.)gke->RemovePoint(i);
      }
      gke->Draw("SAMEL");

      leg2->AddEntry(hke, texlabel[samples[isample]], "FL");
    } // for isample

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
  
  curdir->cd();
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

    // Add chi2 from nuisance parameters //including new multijet
    for (int is = 0; is != ns; ++is) {
      
      chi2 += ps[is]*ps[is];
    } // for ipar

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

