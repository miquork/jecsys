// File: softrad.C
// Created by Mikko Voutilainen, on June 5th, 2014
// Updated on Oct 25, 2014 (cleanup for 8 TeV JEC paper)
// Purpose: Use JEC combination file to derive FSR+ISR corrections,
//          including uncertainty eigenvectors for global refit
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TF1.h"
#include "TMultiGraph.h"
#include "TLatex.h"
#include "TMinuit.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMath.h"

//#include "../tools.h"
#include "tdrstyle_mod14.C"

#include <string>
#include <iostream>
#include <vector>
#include <map>

const double _lumi = 19800.;

bool dodijet = true;//false;
bool _s_dcsonly = false;

using namespace std;

TF1 *_sr_fitError_func(0);
TMatrixD *_sr_fitError_emat(0);
Double_t sr_fitError(Double_t *xx, Double_t *p);

// Soft radiation corrections for L3Res
void softrad(double etamin=0.0, double etamax=1.3, bool dodijet=false) {

  cout << "Calling softrad("<<etamin<<","<<etamax<<","<<dodijet<<");/n"<<flush;

  setTDRStyle();
  writeExtraText = false; // for JEC paper CWR

  TDirectory *curdir = gDirectory;

  // Open jecdata.root produced by reprocess.C
  //TFile *fin = new TFile("rootfiles/jecdata_runI_20fb.root","READ");
  //TFile *fin = new TFile("rootfiles/jecdata.root","READ");
  //assert(fin && !fin->IsZombie());
  //
  //TFile *fout = new TFile("rootfiles/jecdata.root","UPDATE");
  //assert(fout && !fout->IsZombie());
  TFile *finout = new TFile("rootfiles/jecdata.root","UPDATE");
  assert(finout && !finout->IsZombie());
  
  const int ndirs = 3;
  const char* dirs[ndirs] = {"data", "mc", "ratio"};
  const int nmethods = 2;
  const char* methods[nmethods] = {"mpfchs1", "ptchs"};
  const int nsamples = (dodijet ? 4 : 3);
  const char* samples[4] = {"gamjet", "zeejet", "zmmjet", "dijet"};
  //const int nsamples = (dodijet ? 3 : 2);
  //const char* samples[3] = {"gamjet", "zmmjet", "dijet"};
  const int idj = (dodijet ? nsamples-1 : -1);
  string sbin = Form("eta%02.0f-%02.0f",10*etamin,10*etamax);
  const char* bin = sbin.c_str();
  const int nalphas = 4;
  const int alphas[nalphas] = {30, 20, 15, 10};
  // TEMP TEMP TEMP
  // For Run 2 gamma+jet sample compatibility
  //const int nalphas = 3;
  //const int alphas[nalphas] = {30, 20, 10};
  //

  // new^2: 30 40 *50* 60 85 *105* 130 175 *230* 300 400 500 700 1000 1500

  // Z+jet bins
  //const double ptbins1[] = {30, 40, 50, 60, 75, 95, 125, 180, 300, 1000};
  //const double ptbins1[] = {30, 40, 50, 60, 75, 95, 125, 180, 300, 400, 1000};
  const double ptbins1[] = {30, 40, 50, 60, 85, 105, 130, 175, 230,
			    300, 400, 500, 700, 1000, 1500};
  const int npt1 = sizeof(ptbins1)/sizeof(ptbins1[0])-1;
  TH1D *hpt1 = new TH1D("hpt1","",npt1,&ptbins1[0]);
  TProfile *ppt1 = new TProfile("ppt1","",npt1,&ptbins1[0]);

  // gamma+jet bins
  //const double ptbins2[] = {30, 40, 50, 60, 75, 100, 125, 155, 180,
  //		    210, 250, 300, 350, 400, 500, 600, 800};
  //const double ptbins2[] = {40, 60, 85, 100, 130, 175, 250, 300, 2000};
  //const double ptbins2[] = {40, 60, 85, 100, 130, 175, 250, 300, 400,
  //		    500, 700, 1000, 1500};
  const double ptbins2[] = {40, 50, 60, 85, 105, 130, 175, 230, 300, 400,
			    500, 700, 1000, 1500};
  const int npt2 = sizeof(ptbins2)/sizeof(ptbins2[0])-1;
  TH1D *hpt2 = new TH1D("hpt2","",npt2,&ptbins2[0]);
  TProfile *ppt2 = new TProfile("ppt2","",npt2,&ptbins2[0]);

  // dijet bins
  const double ptbins4[] = {20, 62, 107, 175, 242, 310, 379, 467,
			    628, 839, 1121, 1497, 2000};
  const int npt4 = sizeof(ptbins4)/sizeof(ptbins4[0])-1;
  TH1D *hpt4 = new TH1D("hpt4","",npt4,&ptbins4[0]);
  TProfile *ppt4 = new TProfile("ppt4","",npt4,&ptbins4[0]);

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);

  map<string,const char*> texlabel;
  texlabel["gamjet"] = "#gamma+jet";
  texlabel["zeejet"] = "Z#rightarrowee+jet";
  texlabel["zmmjet"] = "Z#rightarrow#mu#mu+jet";
  texlabel["dijet"] = "Dijet";
  texlabel["ptchs"] = "p_{T} balance (CHS)";
  texlabel["mpfchs"] = "MPF raw (CHS)";
  texlabel["mpfchs1"] = "MPF type-I (CHS)";

  // overlay of various alpha values
  TCanvas *c1 = new TCanvas("c1","c1",ndirs*400,nmethods*400);
  c1->Divide(ndirs,nmethods);

  TH1D *h1 = new TH1D("h1",";p_{T} (GeV);Response",1270,30,1300);

  // extrapolation vs alpha for each pT bin
  vector<TCanvas*> c2s(ndirs*nmethods);
  for (unsigned int icanvas = 0; icanvas != c2s.size(); ++icanvas) {
    TCanvas *c2 = new TCanvas(Form("c2_%d",icanvas),Form("c2_%d",icanvas),
			      1200,1200);
    c2->Divide(3,3);
    c2s[icanvas] = c2;
  }

  TH1D *h2 = new TH1D("h2",";#alpha;Response",10,0.,0.4);
  h2->SetMaximum(1.10);//1.08);
  h2->SetMinimum(0.85);//0.88);

  // krad corrections
  TCanvas *c3 = new TCanvas("c3","c3",ndirs*400,nmethods*400);
  c3->Divide(ndirs,nmethods);

  TH1D *h3 = new TH1D("h3",";p_{T,ref} (GeV);FSR sensitivity: -dR/d#alpha [%]",
		      1270,30,1300);

  cout << "Reading in data" << endl << flush;
  // Read in plots vs pT (and alpha)
  map<string, map<string, map<string, map<int, TGraphErrors*> > > > gemap;
  map<string, map<string, map<string, map<int, TGraphErrors*> > > > gamap;
  for (int idir = 0; idir != ndirs; ++idir) {
    for (int  imethod = 0; imethod != nmethods; ++imethod) {
      for (int  isample = 0; isample != nsamples; ++isample) {

	for (int  ialpha = 0; ialpha != nalphas; ++ialpha) {

	  finout->cd();
	  //assert(gDirectory->cd(dirs[idir]));
	  //assert(gDirectory->cd(bin));
	  //TDirectory *d = gDirectory;
	  assert(finout->cd(dirs[idir]));
	  TDirectory *din1 = finout->GetDirectory(dirs[idir]); assert(din1);
	  assert(din1->cd(bin));
	  TDirectory *d = din1->GetDirectory(bin); assert(d);
	  d->cd();

	  const char *cd = dirs[idir];
	  const char *cm = methods[imethod];
	  const char *cs = samples[isample];
	  const int a = alphas[ialpha];
	  // Get graph made vs pT
	  string s = Form("%s/%s/%s_%s_a%d",dirs[idir],bin,cm,cs,a);
	  TGraphErrors *g = (TGraphErrors*)finout->Get(s.c_str());
	  // TEMP TEMP TEMP
	  // Use Run 2 gamma+jet sample for pT balance	  
	  //TGraphErrors *g = 0;
	  //if (string(cs)=="gamjet" && string(cm)=="ptchs") 
	  //if (string(cs)=="gamjet") 
	    //if (string(cs)=="gamjet" || string(cs)=="zmmjet") 
	    //g = (TGraphErrors*)fout->Get(s.c_str());
	  //else
	  //g = (TGraphErrors*)finout->Get(s.c_str());
	  // TEMP TMP TMP
	  if (!g) cout << "Missing " << s << endl << flush;
	  assert(g);

	  // Clean out empty points
	  // as well as trigger-biased ones for dijets
	  // as well as weird gamma+jet high pT point
	  for (int i = g->GetN()-1; i != -1; --i) {
	    if (g->GetY()[i]==0 || g->GetEY()[i]==0 ||
		(string(cs)=="dijet" && g->GetX()[i]<70.)) // ||
		//(string(cs)=="gamjet" && g->GetX()[i]>600. && etamin!=0))
	      g->RemovePoint(i);
	  }

	  gemap[cd][cm][cs][a] = g;
	  
	  // Sort points into new graphs vs alpha
	  TH1D *hpt = (isample==0 ? hpt2 : hpt1);
	  TProfile *ppt = (isample==0 ? ppt2 : ppt1);
	  //if (isample==3) { hpt = hpt4; ppt = ppt4; } // pas-v6
	  if (isample==idj) { hpt = hpt4; ppt = ppt4; } // pas-v6
	  for (int i = 0; i != g->GetN(); ++i) {
	    
	    double pt = g->GetX()[i];
	    ppt->Fill(pt, pt);
	    int ipt = int(hpt->GetBinLowEdge(hpt->FindBin(pt))+0.5);
	    //int ipt = int(pt+0.5);
	    TGraphErrors *ga = gamap[cd][cm][cs][ipt];
	    if (!ga) {
	      ga = new TGraphErrors(0);
	      ga->SetMarkerStyle(g->GetMarkerStyle());
	      ga->SetMarkerColor(g->GetMarkerColor());
	      ga->SetLineColor(g->GetLineColor());
	      gamap[cd][cm][cs][ipt] = ga;
	    }
	    int n = ga->GetN();
	    ga->SetPoint(n, 0.01*a, g->GetY()[i]);
	    ga->SetPointError(n, 0, g->GetEY()[i]);
	  } // for i 

	} // for ialpha

      } // for isample
    } // for imethod
  } // for idir

  cout << "Drawing plots vs pT for each alpha" << endl << flush;

  // 2x6 plots
  for (int idir = 0; idir != ndirs; ++idir) {
    for (int  imethod = 0; imethod != nmethods; ++imethod) {

      const char *cd = dirs[idir];
      const char *cm = methods[imethod];

      int ipad = ndirs*imethod + idir + 1; assert(ipad<=6);
      c1->cd(ipad);
      gPad->SetLogx();
      h1->SetMaximum(1.15);//idir<2 ? 1.15 : 1.08);
      h1->SetMinimum(0.80);//idir<2 ? 0.85 : 0.93);
      h1->SetYTitle(Form("Response (%s)",cd));
      h1->DrawClone("AXIS");
      tex->DrawLatex(0.20,0.85,texlabel[cm]);
      tex->DrawLatex(0.20,0.80,"|#eta| < 1.3, #alpha=0.1--0.3");
      TLegend *leg = tdrLeg(0.60,0.75,0.90,0.90);

      for (int  isample = 0; isample != nsamples; ++isample) {
	for (int  ialpha = 0; ialpha != nalphas; ++ialpha) {

	  const char *cs = samples[isample];
	  if (etamin==0 && string(cs)=="dijet") continue;

	  const int a = alphas[ialpha];
	  TGraphErrors *g = gemap[cd][cm][cs][a]; assert(g);

	  // Clean out points with very large uncertainty for plot readability
	  for (int i = g->GetN()-1; i != -1; --i) {
	    if (g->GetEY()[i]>0.02) g->RemovePoint(i);
	  }

	  g->Draw("SAME Pz");

	  if (ialpha==0) leg->AddEntry(g,texlabel[cs],"P");
	}
      } // for isample

      // Individual plots for JEC paper
      if ( true ) { // paper

	TH1D *h = new TH1D(Form("h_5%s_%s",cd,cm),
			   Form(";p_{T} (GeV);Response (%s)",cd),
			   1270,30,1300);
	h->GetXaxis()->SetMoreLogLabels();
	h->GetXaxis()->SetNoExponent();
	h->SetMinimum(0.83);//0.88);
	h->SetMaximum(1.13);//1.13);
	
	//lumi_13TeV = (_s_dcsonly ?
	//	      "Run2015D-DCSOnly - 25 ns - 980 pb^{-1}" :
	//	      "Run2015D-Cert - 25 ns - 552 pb^{-1}");
	//lumi_13TeV = "Run2015D - Oct 19 - 1.28 fb^{-1}";
	lumi_13TeV = "Run2015D - Golden - 2.1 fb^{-1}";
	TCanvas *c0 = tdrCanvas(Form("c0_%s_%s",cm,cd), h, 4, 11, true);
	c0->SetLogx();
	

	TLegend *leg = tdrLeg(0.55,0.68,0.85,0.83);
	tex->DrawLatex(0.55,0.85,texlabel[cm]);
	tex->DrawLatex(0.55,0.18,"|#eta| < 1.3, #alpha<0.3");
	//tex->DrawLatex(0.55,0.18,"Anti-k_{T} R=0.5");

	// Loop over Z+jet and gamma+jet (only, no dijet/multijet)
	//for (int  isample = 0; isample != min(3,nsamples); ++isample) {
	for (int  isample = 0; isample != min(idj,nsamples); ++isample) {
	  
	  const char *cs = samples[isample];
	  TGraphErrors *g = gemap[cd][cm][cs][30]; assert(g);
	  g->Draw("SAME Pz");
	  
	  leg->AddEntry(g,texlabel[cs],"P");
	} // for isample

	if (etamin==0) {
	  c0->SaveAs(Form("pdf/paper_softrad_%s_%s_vspt.pdf",cd,cm));
	  c0->SaveAs(Form("pdfC/paper_softrad_%s_%s_vspt.C",cd,cm));
	}
	else {
	  c0->SaveAs(Form("pdf/an_softrad_%s_%s_eta%1.0f-%1.0f_vspt.pdf",
			  cd,cm,10*etamin,10*etamax));
	}
      } // paper

    } // for imethod
  } // for idir
  
  c1->cd(0);
  //cmsPrel(_lumi, true);
  //CMS_lumi(c1, 2, 33);
  c1->SaveAs("pdf/softrad_2x6_vspt.pdf");


  cout << "Drawing plots vs alpha for each pT" << endl << flush;
  cout << "...and fitting slope vs alpha" << endl << flush;

  map<string, map<string, map<string, TGraphErrors* > > > gkmap;

  // 2x6 plots
  for (int idir = 0; idir != ndirs; ++idir) {
    for (int  imethod = 0; imethod != nmethods; ++imethod) {
      
      int icanvas = nmethods*imethod + idir; assert(icanvas<=6);
      TCanvas *c2 = c2s[icanvas]; assert(c2);

      const char *cd = dirs[idir];
      string dd = cd;
      const char *cm = methods[imethod];
      string mm = cm;

      const int npads = 9;
      for (int ipad = 0; ipad != npads; ++ipad) {
	c2->cd(ipad+1);
	h2->SetYTitle(Form("Response (%s)",cd));
	h2->DrawClone("AXIS");
	tex->DrawLatex(0.20,0.85,texlabel[cm]);
	tex->DrawLatex(0.20,0.80,"|#eta| < 1.3");
	tex->DrawLatex(0.20,0.75,Form("%1.0f < p_{T} < %1.0f GeV",
				      hpt1->GetBinLowEdge(ipad+1),
				      hpt1->GetBinLowEdge(ipad+2)));
	TLegend *leg = tdrLeg(0.65,0.75,0.90,0.90);
	leg->AddEntry(gemap[cd][cm]["gamjet"][30], texlabel["gamjet"], "P");
	//leg->AddEntry(gemap[cd][cm]["zeejet"][30], texlabel["zeejet"], "P");
	leg->AddEntry(gemap[cd][cm]["zmmjet"][30], texlabel["zmmjet"], "P");
	leg->AddEntry(gemap[cd][cm]["dijet"][30], texlabel["dijet"], "P");
      }

      for (int  isample = 0; isample != nsamples; ++isample) {

	const char *cs = samples[isample];
	string ss = cs;

	map<int, TGraphErrors*> &gam = gamap[cd][cm][cs];
	map<int, TGraphErrors*>::iterator itpt;
	for (itpt = gam.begin(); itpt != gam.end(); ++itpt) {

	  int ipt = itpt->first;
	  int jpt = hpt1->FindBin(ipt);
	  if (jpt>npads) continue;
	  assert(jpt<=npads);
	  c2->cd(jpt);
	  
	  TGraphErrors *ga = itpt->second; assert(ga);
	  
	  ga->Draw("SAME Pz");

	  // Fit slope
	  TF1 *f1 = new TF1(Form("f1_%s_%s_%s_%d",cd,cm,cs,ipt),
			    "(x<1)*([0]+[1]*x) + (x>1 && x<2)*[0] +"
			    "(x>2)*[1]",-1,1);
	  f1->SetLineColor(ga->GetLineColor());
	  f1->SetParameters(1,0);
	  //const double minalpha = (isample==0 ? 10./ipt : 5./ipt);
	  const double minalpha = (isample==0 ? 10./ipt : 5./ipt); // 2015-10-20
	  // Constrain slope to within reasonable values
	  // in the absence of sufficient data using priors
	  // TEMP TEMP TEMP
	  //if (!(ss=="gamjet" && mm=="ptchs")) {//true) { // use priors
	  //if (!(ss=="gamjet")) {//true) { // use priors
	  if (true) { // use priors
	  // turn off priors
	  //if (false) { // use priors
	  //if (!(string(cs)=="gamjet" && string(cm)=="ptchs")) {
	  //if (!(string(cs)=="gamjet")) {
	  //if (false) {
	  // TEMP TEMP TEMP
	    int n = ga->GetN();
	    // For response, limit to 1+/-0.02 (we've corrected for L3Res
	    //ga->SetPoint(n, 1.5, 1);
	    //ga->SetPointError(n, 0, 0.02);
	    //n = ga->GetN();
	    //if (imethod==1) { // pT balance
	    if (mm=="ptchs") {
	      // For pT balance, estimate slope of <vecpT2>/alpha from data
	      // => 7.5%/0.30 = 25%
	      // Approximate uncertainty on this to be
	      // 0.5%/0.30 ~ 1.5% for data, 0.5%/0.30 ~ 1.5% for Z+jet MC, and
	      // 2%/0.30 ~ 6% for gamma+jet MC (same as slope)
	      /*
	      if (idir==0 && isample!=idj) ga->SetPoint(n, 2.5, -0.300); //DT
	      if (idir==0 && isample==idj) ga->SetPoint(n, 2.5, -0.075);
	      if (idir==1 && isample!=idj) ga->SetPoint(n, 2.5, -0.300); //MC
	      if (idir==1 && isample==idj) ga->SetPoint(n, 2.5, -0.075);
	      if (idir==2 && isample!=idj) ga->SetPoint(n, 2.5, -0.000); //rt
	      if (idir==2 && isample==idj) ga->SetPoint(n, 2.5, -0.000); 
	      */
	      // x=2.5 constraints the slope dR/dalpha; set expexted value
	      if (dd=="data"  && ss!="dijet") ga->SetPoint(n, 2.5, -0.300);
	      if (dd=="mc"    && ss!="dijet") ga->SetPoint(n, 2.5, -0.300);
	      if (dd=="ratio" && ss!="dijet") ga->SetPoint(n, 2.5, -0.000);
	      //
	      if (dd=="data"  && ss=="dijet") ga->SetPoint(n, 2.5, -0.075);
	      if (dd=="mc"    && ss=="dijet") ga->SetPoint(n, 2.5, -0.075);
	      if (dd=="ratio" && ss=="dijet") ga->SetPoint(n, 2.5, -0.000); 
	      // For V7, consider different slope for gamma+jet
	      if (dd=="data"  && ss=="gamjet") ga->SetPoint(n, 2.5, -0.275);
	      if (dd=="mc"    && ss=="gamjet") ga->SetPoint(n, 2.5, -0.300);
	      if (dd=="ratio" && ss=="gamjet") ga->SetPoint(n, 2.5, +0.025); 
	      //
	      // BUG: found 2015-01-08 (no effect on ratio)
	      // BUG: found 2015-10-20 (negative errors made them inapplicable)
	      //if (idir==1)               ga->SetPointError(n, 0, -0.015);
	      /*
	      if (idir==0 && isample!=idj) ga->SetPointError(n, 0, 0.050); //DT
	      if (idir==0 && isample==idj) ga->SetPointError(n, 0, 0.025);
	      if (idir==1 && isample!=idj) ga->SetPointError(n, 0, 0.050); //MC
	      if (idir==1 && isample==idj) ga->SetPointError(n, 0, 0.025);
	      if (idir==2 && isample!=idj) ga->SetPointError(n, 0, 0.025); //rt
	      if (idir==2 && isample==idj) ga->SetPointError(n, 0, 0.0125); 
	      */
	      // x=2.5 for slope of dR/dalpha; set the a priori uncertainty
	      if (dd=="data"  && ss!="dijet") ga->SetPointError(n, 0, 0.050);
	      if (dd=="mc"    && ss!="dijet") ga->SetPointError(n, 0, 0.050);
	      if (dd=="ratio" && ss!="dijet") ga->SetPointError(n, 0, 0.025);
	      //
	      if (dd=="data"  && ss=="dijet") ga->SetPointError(n, 0, 0.025);
	      if (dd=="mc"    && ss=="dijet") ga->SetPointError(n, 0, 0.025);
	      if (dd=="ratio" && ss=="dijet") ga->SetPointError(n, 0, 0.0125);
	      // For V7, consider different slope for gamma+jet
	      if (dd=="data"  && ss=="gamjet") ga->SetPointError(n, 0, 0.050);
	      if (dd=="mc"    && ss=="gamjet") ga->SetPointError(n, 0, 0.050);
	      if (dd=="ratio" && ss=="gamjet") ga->SetPointError(n, 0, 0.050);
	    }
	    //if (imethod==0) { // MPF
	    if (mm=="mpfchs1") { // MPF
	      // For MPF, expectation is no slope
	      // Maximal slope would be approximately
	      // (<vecpT2>/alpha ~ 25% from pT balance) times
	      // (response difference between pT1 and vecpT2~10%)
	      // => 0.25*0.10 = 2.5%
	      // For data/MC, estimate uncertainty as half of this
	      // => 1.25%
	      /*
	      ga->SetPoint(n, 2.5, 0.);
	      if (idir!=2) ga->SetPointError(n, 0, 0.025);
	      if (idir==2) ga->SetPointError(n, 0, 0.0125);
	      */
	      // V6 constraint: no slope
	      /*
	      ga->SetPoint(n, 2.5, 0.);
	      if (dd!="ratio") ga->SetPointError(n, 0, 0.025);
	      if (dd=="ratio") ga->SetPointError(n, 0, 0.0125);
	      */
	      // V7 constraint: HF down by -20%, negative slope from gamjet
	      ga->SetPoint(n, 2.5, dd=="mc" ? 0. : -0.010);
	      if (dd!="mc") ga->SetPointError(n, 0, 0.0125);
	      if (dd=="mc") ga->SetPointError(n, 0, 0.0125);
	    } // MPF
	  } // use priors

	  if (ga->GetN()>2) {

	    f1->SetRange(minalpha, 3.);
	    ga->Fit(f1,"QRN");

	    if (f1->GetNDF()>=0) {
	      f1->DrawClone("SAME");
	      f1->SetRange(0,0.4);
	      f1->SetLineStyle(kDashed);
	      f1->DrawClone("SAME");

	      // Store results
	      TGraphErrors *gk = gkmap[cd][cm][cs];
	      if (!gk) {
		gk = new TGraphErrors(0);
		gk->SetMarkerStyle(ga->GetMarkerStyle());
		gk->SetMarkerColor(ga->GetMarkerColor());
		gk->SetLineColor(ga->GetLineColor());
		gkmap[cd][cm][cs] = gk;
	      }
	      int n = gk->GetN();
	      TProfile *ppt = (isample==0 ? ppt2 : ppt1);
	      //if (isample==3) { ppt = ppt4; } // pas-v6
	      if (isample==idj) { ppt = ppt4; } // pas-v6
	      double pt = ppt->GetBinContent(ppt->FindBin(ipt));
	      gk->SetPoint(n, pt, f1->GetParameter(1));
	      gk->SetPointError(n, 0, f1->GetParError(1));
	    } // f1->GetNDF()>=0
	  } // ga->GetN()>2
	} // for itpt
	
      } // for isample
      
      c2->SaveAs(Form("pdf/softrad_3x3_%s_%s_vsalpha.pdf",cd,cm));
      
    }
  }


  cout << "Drawing plots of kFSR vs pT" << endl;

  // 2x6 plots
  for (int idir = 0; idir != ndirs; ++idir) {
    for (int  imethod = 0; imethod != nmethods; ++imethod) {

      const char *cd = dirs[idir];
      const char *cm = methods[imethod];

      TMultiGraph *mgk = new TMultiGraph();

      int ipad = ndirs*imethod + idir + 1; assert(ipad<=6);
      c3->cd(ipad);
      gPad->SetLogx();
      h3->SetMaximum(imethod==0 ? 0.05 : (idir!=2 ? 0.1 : 0.25));
      h3->SetMinimum(imethod==0 ? -0.05 : (idir!=2 ? -0.4 : -0.25));
      h3->SetYTitle(Form("k_{FSR} = dR/d#alpha (%s)",cd));
      h3->DrawClone("AXIS");
      tex->DrawLatex(0.20,0.85,texlabel[cm]);
      tex->DrawLatex(0.20,0.80,"|#eta| < 1.3");
      TLegend *leg = tdrLeg(0.60,0.75,0.90,0.90);

      for (int  isample = 0; isample != nsamples; ++isample) {

	const char *cs = samples[isample];
	TGraphErrors *gk = gkmap[cd][cm][cs]; assert(gk);
	
	leg->AddEntry(gk,texlabel[cs],"P");

	// Fit each sample separately for pT balance
	if (true) {

	  TF1 *fk = new TF1(Form("fk_%s_%s_%s",cd,cm,cs),
			    "[0]+[1]*log(0.01*x)+[2]*pow(log(0.01*x),2)",
			    30,1300);
	  // TEMP TEMP TEMP
	  fk->SetParameters(-0.05,+0.01,0);//-0.25,-0.5,0);
	  //fk->FixParameter(2,0);
	  //if (etamin>=2.5) fk->FixParameter(1,0);
	  // Change parameters for Run 2
	  //fk->SetParameters(0,0,0);
	  // TEMP TEMP TEMP
	  fk->SetLineColor(gk->GetLineColor());
	  gk->Fit(fk, "QRN");

	  tex->SetTextColor(fk->GetLineColor());
	  tex->DrawLatex(0.55,0.27-0.045*isample,
			 Form("#chi^{2}/NDF = %1.1f / %d",
			      fk->GetChisquare(), fk->GetNDF()));
	  tex->SetTextColor(kBlack);

	  // Error band
	  const int n = fk->GetNpar();
	  TMatrixD emat(n,n);
	  gMinuit->mnemat(emat.GetMatrixArray(), n);
	  TF1 *fke = new TF1(Form("fk_%s_%s_%s",cd,cm,cs),
			     sr_fitError, 30, 1300, 1);
	  _sr_fitError_func = fk;
	  _sr_fitError_emat = &emat;

	  fke->SetLineStyle(kSolid);
	  fke->SetLineColor(fk->GetLineColor()-10);
	  fke->SetParameter(0,-1);
	  fke->DrawClone("SAME");
	  fke->SetParameter(0,+1);
	  fke->DrawClone("SAME");

	  fk->DrawClone("SAME");
	  gk->DrawClone("SAME Pz");

	  // Store soft radiation corrections in fsr subdirectory
	  assert(finout->cd(cd));
	  TDirectory *dout1 = finout->GetDirectory(cd); assert(dout1);
	  //assert(gDirectory->cd(bin)); // broke in ROOT 6.04/08
	  //if (!gDirectory->FindObject("fsr")) gDirectory->mkdir("fsr");
	  //assert(gDirectory->cd("fsr"));
	  assert(dout1->cd(bin));
	  TDirectory *dout2 = dout1->GetDirectory(bin); assert(dout2);
	  if (!dout2->FindObject("fsr")) dout2->mkdir("fsr");
	  assert(dout2->cd("fsr"));
	  TDirectory *dout3 = dout2->GetDirectory("fsr"); assert(dout3);
	  dout3->cd();

	  TH1D *hk = (TH1D*)(isample==0 ? hpt2->Clone() : hpt1->Clone());
	  hk->SetName(Form("hkfsr_%s_%s",cm,cs));
	  TProfile *ppt = (isample==0 ? ppt2 : ppt1);
	  //if (isample==3) { ppt = ppt4; } // pas-v6
	  if (isample==2) { ppt = ppt4; } // pas-v6
	  for (int i = 1; i != hk->GetNbinsX()+1; ++i) {
	    double pt = ppt->GetBinContent(i);
	    if (pt>30 && pt<1300) {
	      hk->SetBinContent(i, fk->Eval(pt));
	      hk->SetBinError(i, fabs(fke->Eval(pt)-fk->Eval(pt)));
	    }
	    else {
	      hk->SetBinContent(i, 0);
	      hk->SetBinError(i, 0);
	    }
	  }
	  
	  hk->Write(hk->GetName(), TObject::kOverwrite);

	  // Factorize error matrix into eigenvectors
	  // Remember: A = Q*Lambda*Q^-1, where
	  // A is emat, Q is eigmat, and Lambda is a diagonal matrix with
	  // eigenvalues from eigvec on the diagonal. For eigenmatrix
	  // Q^-1 = Q^T, i.e. inverse matrix is the original transposed
	  TVectorD eigvec(n);
	  TMatrixD eigmat = emat.EigenVectors(eigvec);

	  // Eigenvectors are the columns and sum of eigenvectors squared
	  // equals original uncertainty. Calculate histograms from the
	  // eigenvectors and store them
	  TF1 *fkeig = (TF1*)fk->Clone(Form("%s_eig",fk->GetName()));
	  fkeig->SetLineStyle(kDotted);
	  for (int ieig = 0; ieig != n; ++ieig) {

	    // Eigenvector functions
	    for (int i = 0; i != n; ++i) {
	      fkeig->SetParameter(i, fk->GetParameter(i)
				  + eigmat[i][ieig] * sqrt(eigvec[ieig]));
	    }
	    fkeig->DrawClone("SAMEL");

	    // Eigenvector histograms evaluated at bin mean pT
	    TH1D *hke = (TH1D*)hk->Clone(Form("%s_eig%d",hk->GetName(),ieig));
	    hke->Reset();

	    for (int i = 0; i != gk->GetN(); ++i) {

	      double pt = gk->GetX()[i];
	      int ipt = hke->FindBin(pt);
	      // Need to store central value as well, because
	      // uncertainty sources are signed
	      hke->SetBinContent(ipt, fkeig->Eval(pt)-fk->Eval(pt));
	      hke->SetBinError(ipt, fabs(fkeig->Eval(pt)-fk->Eval(pt)));
	    }
	    hke->Write(hke->GetName(), TObject::kOverwrite);
	  }

	  cout << "." << flush;
	} // if tree
      } // for isample
    } // for imethod
  } // for idir
  
  c3->cd(0);
  //cmsPrel(_lumi, true);
  CMS_lumi(c3, 2, 33);
  c3->SaveAs("pdf/softrad_2x6_kfsr.pdf");

  finout->Close();
  //fout->Close();
  curdir->cd();
} // softrad

Double_t sr_fitError(Double_t *xx, Double_t *p) {

  assert(_sr_fitError_func);
  assert(_sr_fitError_emat);
  double x = *xx;
  double k = p[0];
  TF1 *f = _sr_fitError_func;
  int n = f->GetNpar();
  TMatrixD &emat = (*_sr_fitError_emat);
  assert(emat.GetNrows()==n);
  assert(emat.GetNcols()==n);
  
  vector<double> df(n);
  for (int i = 0; i != n; ++i) {

    double p = f->GetParameter(i);
    double dp = 0.1*f->GetParError(i);
    f->SetParameter(i, p+dp);
    double fup = f->Eval(x);
    f->SetParameter(i, p-dp);
    double fdw = f->Eval(x);
    f->SetParameter(i, p);
    df[i] = (dp ? (fup - fdw) / (2.*dp) : 0);
  }

  double sumerr2(0);
  for (int i = 0; i != n; ++i) {
    for (int j = 0; j != n; ++j) {
      sumerr2 += emat[i][j]*df[i]*df[j];
    }
  }

  double err = sqrt(sumerr2);

  return (f->Eval(x) + k*err);
}
