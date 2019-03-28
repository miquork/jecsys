// Purpose: Set of tools to derive MC truth JER and JER SF from SMP-J tuples
//          to be used for inclusive jet analysis (d|y|=0.5 binning)
// Related: ptresolution.h stores results and interfaces JME official JER
//          mk_resolution.C runs this code and interfaces with JER packages
//          ../pdf/resolutionslides.tex to plot results
// author: mikko.voutilainen at cern.ch

// Specific tools:
// - resolution(): derive MC truth JER and compare own JER to JME official
// - redoJER(): average JME JER SF over eta for inclusive jet rapidity bins
// - resolution_datamc(): derive JER SF (to be maintained)
//
// run with 'root -l -b -q minitools/mk_resolution.C'

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile2D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"

#include "ptresolution.h"
#include "../tdrstyle_mod15.C"
#//include "settings12.h"
#include "tools.h"

#include <iostream>

using namespace std;

bool _closejer = true;
double _recopt = 15.;

Double_t ptreso(Double_t *x, Double_t *p) {
  return ptresolution(x[0], p[0])*p[1];
}

// http://en.wikipedia.org/wiki/Crystal_Ball_function
Double_t fCrystalBall(Double_t *xx, Double_t *p) {
  
  double x = xx[0];
  double N0 = p[0]; // Overall normalization
  double alpha = p[1]; // start of non-Gaussian tail in units of sigma
  double a = fabs(alpha);
  double n = p[2]; // ??
  double xbar = p[3]; // core Gaussian mean
  double sigma = p[4]; // core Gaussian width

  double A = pow(n / a, n) * exp(-a*a/2);
  double B = n / a - a;
  double C = n / a * 1 / (n-1) * exp(-a*a/2);
  double D = sqrt(TMath::Pi()/2) * (1 + TMath::Erf(a/sqrt(2)));
  double N = 1. / (sigma * (C + D));
		   
  if ( (x-xbar)/sigma > -alpha) {
    return ( N0 * N * exp(-pow(x-xbar,2)/(2*sigma*sigma)) );
  }
  if ( (x-xbar)/sigma <= -alpha) {
    return ( N0 * N * A * pow(B - (x-xbar)/sigma, -n) );
  }
  
  //assert(false);
  return 0;
}
// Double-sided version of fCrystalBall (see above)
Double_t fCrystalBall2(Double_t *xx, Double_t *p) {
  
  double x = xx[0];
  double N0 = p[0]; // Overall normalization
  double alpha1 = p[1]; // start of non-Gaussian tail in units of sigma
  double a1 = fabs(alpha1);
  double n1 = p[2]; // ??
  double alpha2 = p[3];
  double a2 = fabs(alpha2);
  double n2 = p[4];
  double xbar = p[5]; // core Gaussian mean
  double sigma = p[6]; // core Gaussian width

  double A1 = pow(n1 / a1, n1) * exp(-a1*a1/2);
  double B1 = n1 / a1 - a1;
  double C1 = n1 / a1 * 1 / (n1-1) * exp(-a1*a1/2);
  double D1 = sqrt(TMath::Pi()/2) * (1 + TMath::Erf(a1/sqrt(2)));
  double N1 = 1. / (sigma * (C1 + D1));

  double A2 = pow(n2 / a2, n2) * exp(-a2*a2/2);
  double B2 = n2 / a2 - a2;
  double C2 = n2 / a2 * 1 / (n2-1) * exp(-a2*a2/2);
  double D2 = sqrt(TMath::Pi()/2) * (1 + TMath::Erf(a2/sqrt(2)));
  double N2 = 1. / (sigma * (C2 + D2));
		   
  if ( (x-xbar)/sigma > -alpha1 && x-xbar<=0) {
    return ( N0 * N1 * exp(-pow(x-xbar,2)/(2*sigma*sigma)) );
  }
  if ( (x-xbar)/sigma <= -alpha1) {
    return ( N0 * N1 * A1 * pow(B1 - (x-xbar)/sigma, -n1) );
  }
  if ( (x-xbar)/sigma < alpha2 && x-xbar>0) {
    return ( N0 * N2 * exp(-pow(x-xbar,2)/(2*sigma*sigma)) );
  }
  if ( (x-xbar)/sigma >= alpha2) {
    return ( N0 * N2 * A2 * pow(B2 - (xbar-x)/sigma, -n2) );
  }
  
  //assert(false);
  return 0;
}

// D0 Jet function
// http://lib.tkk.fi/Diss/2008/isbn9789521037238/isbn9789521037238.pdf
// page 146 (171), Eq.(7.18)
Double_t fD0Jet(Double_t *xx, Double_t *p) {
  
  double x = xx[0];
  double N = p[0];
  double mu = p[1];
  double sigma = p[2];
  double P = p[3];
  double lambda = p[4];

  // we have <x> = mu - P/lambda, RMS(x) = sqrt(sigma^2 + P(2-P)/lambda^2)

  double g = (1-P)*TMath::Gaus(x,mu,sigma,kTRUE)
    + P*lambda/2 * exp( lambda * (x - mu + lambda*sigma*sigma/2) )
    * TMath::Erfc( (x - mu + lambda*sigma*sigma) / (sqrt(2)*sigma));

  return ( N * g );
}
// Double-sided version
Double_t fD0Jet2(Double_t *xx, Double_t *p) {
  
  double x = xx[0];
  double N = p[0];
  double mu = p[1];
  double sigma = p[2];
  double P1 = p[3];
  double lambda1 = p[4];
  double P2 = p[5];
  double lambda2 = p[6];

  // we have <x> = mu - P/lambda, RMS(x) = sqrt(sigma^2 + P(2-P)/lambda^2)
  // => how did this change with two sides?

  double g = (1-P1-P2)*TMath::Gaus(x,mu,sigma,kTRUE)
    + P1*lambda1/2 * exp( lambda1 * (x - mu + lambda1*sigma*sigma/2) )
    * TMath::Erfc( (x - mu + lambda1*sigma*sigma) / (sqrt(2)*sigma))
    + P2*lambda2/2 * exp( -lambda2 * (x - mu - lambda2*sigma*sigma/2) )
    * TMath::Erfc( -(x - mu - lambda2*sigma*sigma) / (sqrt(2)*sigma));

  return ( N * g );
}

Double_t ptresPlusN(Double_t *x, Double_t *p) {
  double res = ptresolution(x[0],p[1]);
  return sqrt(res*res + p[0]*p[0]/(x[0]*x[0]))/res;
}

Double_t ptresPlusS(Double_t *x, Double_t *p) {
  double res = ptresolution(x[0],p[1]);
  return sqrt(res*res + p[0]*p[0]/x[0])/res;
}

Double_t ptresPlusC(Double_t *x, Double_t *p) {
  double res = ptresolution(x[0],p[1]);
  return sqrt(res*res + p[0]*p[0])/res;
}

Double_t ptresPlusNSC(Double_t *x, Double_t *p) {
  double res = ptresolution(x[0],p[3]);
  return sqrt(res*res + p[0]*p[0]/(x[0]*x[0])
	      + p[1]*p[1]/x[0] + p[2]*p[2])/res;
}

void resolution(string type="MC",string file="") {

  //_ismcjer = true;
  _jer_iov = run1;
  if (TString(file.c_str()).Contains("Fall18")) _jer_iov = run2018;
  if (TString(file.c_str()).Contains("17nov17")) _jer_iov = run2017;
  if (TString(file.c_str()).Contains("Legacy16")) _jer_iov = run2016;
  jer_iov jer_ref = _jer_iov;

  TDirectory *curdir = gDirectory;
  setTDRStyle();

  string filename = type + "-1-" + file;
  const char *cff = filename.c_str();
  string idname = type + "-" + file;
  const char *cf = idname.c_str();

  const char *c5 = "R=0.4, 13 TeV";
  //TFile *fin5 = new TFile("rootfiles/output-MC-1-Fall18V8-D.root", "READ");
  //TFile *fin5 = new TFile(Form("rootfiles/output-MC-1-%s.root",cf),"READ");
  TFile *fin5 = new TFile(Form("rootfiles/output-%s.root",cff),"READ");
  assert(fin5 && !fin5->IsZombie());
  assert(fin5->cd("Standard"));
  TDirectory *din5 = gDirectory;

  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextSize(0.045);

  const int ny0 = 0;
  const int ny = 8;//6;
  const double ybins[ny+1] = {0,0.5,1,1.5,2,2.5,3,3.2,4.7};

  double vpar5[ny][3];
  double vchi5[ny];
  int vndf5[ny];

  // Summary plots of AK5 and AK7 JEC
  TCanvas *c0b = new TCanvas("c0b","c0b",1200,600);
  c0b->Divide(2);

  TCanvas *c3 = new TCanvas("c3","c3",600,600);
  gPad->SetLogx();

  TH1D *h3 = new TH1D("h3",";p_{T}^{ptcl} (GeV);Resolution",2480,20,2500);
  h3->SetMinimum(0);
  h3->SetMaximum(0.25);
  h3->GetXaxis()->SetMoreLogLabels();
  h3->GetXaxis()->SetNoExponent();
  h3->Draw("AXIS");

  TCanvas *c3e = new TCanvas("c3e","c3e",600,600);
  gPad->SetLogx();

  TLegend *leg3 = new TLegend(0.6,0.6,0.8,0.9,"","brNDC");
  leg3->SetFillStyle(kNone);
  leg3->SetBorderSize(0);
  leg3->SetTextSize(0.045);
  leg3->Draw();

  TH1D *h3e = new TH1D("h3e",";E (GeV);Resolution",3480,20,3500);
  h3e->SetMinimum(0);
  h3e->SetMaximum(0.25);
  h3e->GetXaxis()->SetMoreLogLabels();
  h3e->GetXaxis()->SetNoExponent();
  h3e->Draw("AXIS");

  TLegend *leg3e = new TLegend(0.6,0.6,0.8,0.9,"","brNDC");
  leg3e->SetFillStyle(kNone);
  leg3e->SetBorderSize(0);
  leg3e->SetTextSize(0.045);
  leg3e->Draw();

  int colors[8] = {kBlack, kGray+1, kRed, kRed+1, kBlue, kBlue+1,
		   kOrange+1, kOrange+3};
  int styles[8] = {kSolid, kSolid, kSolid, kDashed, kDashed, kDotted,
		   kDashDotted, kDashDotted};

  for (int iy = ny0; iy != ny; ++iy) {


    //const double etamin = 0.5*iy;
    //const double etamax = 0.5*(iy+1);
    const double etamin = ybins[iy];
    const double etamax = ybins[iy+1];
    const double etamid = 0.5*(etamin+etamax);
    const double ptmin = 15.;//20;
    const double fitxmin = 20;//37.;
    const double emax = 3000.;
    const double ptmax = emax/cosh(etamin);

    assert(din5->cd(Form("Eta_%1.1f-%1.1f", etamin, etamax)));
    assert(gDirectory->cd("mc"));
    TDirectory *d5 = gDirectory;

    curdir->cd();
    
    TH1D *hpt5 = (TH1D*)d5->Get("hpt_g0"); assert(hpt5);
    TH2D *h2r5 = (TH2D*)d5->Get("h2r_g"); assert(h2r5);

    // Plot and fit resolution
    TCanvas *c1 = new TCanvas(Form("c1_%d",iy),Form("c1_%d",iy),8*150,5*150);
    c1->Divide(8,5,-1,-1);
    const int nmax = 8*5;
    
    TH1D *h = new TH1D(Form("h_%d",iy),";p_{T}^{reco} / p_{T}^{gen};"
		       "1/N/dR_{jet}",200,0,2);
    h->SetMinimum(1e-5);
    h->SetMaximum(50.);

    TF1 *fit5 = new TF1(Form("fit5_%d",iy),"TMath::Gaus(x,[1],[0],1)",0.,2.);
    fit5->SetLineColor(kBlue);
    
    TGraphErrors *gr5 = new TGraphErrors();
    TGraphErrors *gr5a = new TGraphErrors();
    TGraphErrors *gr5b = new TGraphErrors();
    //
    TGraphErrors *gs5 = new TGraphErrors();
    TGraphErrors *gs5a = new TGraphErrors();

    //const int nskip = h2r5->GetXaxis()->FindBin(21.)-1;
    const int nskip = h2r5->GetXaxis()->FindBin(ptmin)-1;
    assert(nskip<=h2r5->GetNbinsX());
    for (int ibin = 1+nskip; ibin != h2r5->GetNbinsX()+1 && ibin-nskip<nmax+1;
	 ++ibin) {

      assert(ibin-nskip<nmax+1);
      c1->cd(ibin-nskip);
      gPad->SetLogy();
      h->Draw("AXIS");

      TH1D *hpr5 = h2r5->ProjectionY(Form("hpr5_%d_%d",iy,ibin),ibin,ibin);
      hpr5->Rebin(4);
      hpr5->Scale(1./hpr5->Integral()/hpr5->GetBinWidth(1));
      hpr5->SetLineColor(kRed+2);
      hpr5->SetFillColor(kRed-9);
      hpr5->SetFillStyle(1001);
      hpr5->DrawClone("SAMEHIST");

      TLatex *tex = new TLatex();
      tex->SetNDC();
      tex->DrawLatex(0.20,0.90,Form("%1.0f<p_{T}<%1.0f GeV",
				    h2r5->GetXaxis()->GetBinLowEdge(ibin),
				    h2r5->GetXaxis()->GetBinLowEdge(ibin+1)));
      double pt = h2r5->GetXaxis()->GetBinLowEdge(ibin);
      double ptave = 0.5*(pt + h2r5->GetXaxis()->GetBinLowEdge(ibin+1));
      double errpt = 0.5*0.5*(h2r5->GetXaxis()->GetBinLowEdge(ibin+1) - pt);
      //double minx = ptmin*1.02/pt;
      double minx = max(0.66, ptmin*1.02/pt);
      double maxx = 1.5;

      // Two-sided Crystal Ball fit
      TF1 *fcb2 = new TF1("fcb2",fCrystalBall2,
			  //0.35, 1.65, 7);
			  minx, maxx, 7);
      //fcb2->SetParameters(1., 1.7, 30, 1.5, 30, 1.008, 0.108);
      fcb2->SetParameters(1., 2.0, 30, 3.0, 30, 1., ptresolution(pt,etamin));
      fcb2->SetParLimits(1, 0.5, 5);//3);
      fcb2->SetParLimits(3, 0.5, 5);//3);
      fcb2->SetLineStyle(kSolid);
      fcb2->SetLineWidth(2);
      fcb2->SetLineColor(kBlack);
      hpr5->Fit(fcb2,"QRN");

      fcb2->SetLineStyle(kSolid);
      fcb2->SetLineWidth(1);
      fcb2->DrawClone("SAME");
      fcb2->SetRange(0.2,1.8);
      fcb2->SetLineStyle(kDotted);
      fcb2->DrawClone("SAME");

      const double srange_up = 2.0;//3.0;
      const double srange_dw = 2.0;
      TF1 *fg5 = new TF1(Form("fg5_%d_%d",iy,ibin),"gaus",
			 max(minx,1-srange_dw*ptresolution(pt,etamin)),
			 min(maxx,1+srange_up*ptresolution(pt,etamin)));
      //fg5->SetParLimits(1,0.98,1.02);
      hpr5->Fit(fg5,"QRN");
      fg5->SetLineWidth(1);
      fg5->SetLineColor(kRed+1);
      fg5->SetLineStyle(kSolid);
      fg5->DrawClone("SAME");
      fg5->SetRange(minx,maxx);
      fg5->SetLineStyle(kDotted);
      fg5->DrawClone("SAME");

      int n = gr5->GetN();
      double k5 = max(1.,sqrt(fg5->GetChisquare()/fg5->GetNDF()));
      gr5->SetPoint(n, ptave, fg5->GetParameter(1));
      gr5->SetPointError(n, errpt, fg5->GetParError(1)*k5);
      gr5a->SetPoint(n, ptave, hpr5->GetMean());
      gr5a->SetPointError(n, errpt, hpr5->GetMeanError());
      gr5b->SetPoint(n, ptave, 0.5*(hpr5->GetMean()+fg5->GetParameter(1)));
      gr5b->SetPointError(n, errpt, tools::oplus(0.5*hpr5->GetMeanError(),
						 0.5*fg5->GetParError(1)*k5));

      gs5->SetPoint(n, ptave, fg5->GetParameter(2));
      gs5->SetPointError(n, errpt, fg5->GetParError(2)*k5);
      gs5a->SetPoint(n, ptave, hpr5->GetRMS());
      gs5a->SetPointError(n, errpt, hpr5->GetRMSError());

      TLegend *leg1 = new TLegend(0.05,1.0,0.6,5,"","br");
      leg1->SetFillStyle(kNone);
      leg1->SetBorderSize(0);
      leg1->SetTextSize(0.045);
      leg1->AddEntry(hpr5,Form("%s (%1.1f)",
			       c5, fg5->GetChisquare()/fg5->GetNDF()),"F");
      leg1->Draw();

      // JER paper plots
      if ((iy==0 && ibin==nskip+12) || // eta 0-1.3, pT 97-114
	  (iy==0 && ibin==nskip+37)) { // eta 0-1.3, pT 967-1032

	bool b1 = (ibin==nskip+12);

	TH1D *h0 = new TH1D(Form("h0_%d",ibin),
			    ";R_{ptcl} = p_{T} / p_{T,ptcl};"
			    "1 / (N dR_{ptcl})",200,0,2); // paper v5
	h0->SetMaximum(39.999);
	h0->SetMinimum(1e-4);

	lumi_8TeV = "";
	extraText = "Simulation"; // paper v5
	//extraText = "Simulation Preliminary";
	//extraText2 = "Preliminary";
	TCanvas *c0 = tdrCanvas("c0",h0,2,0,kSquare);
	h0->GetYaxis()->SetTitleOffset(1.1);
	h0->GetXaxis()->SetTitleOffset(0.95);
	gPad->SetLogy();

	// Crystal Ball fit
	double mu = fg5->GetParameter(1);
	double sigma = fg5->GetParameter(2);
	TF1 *fcb = new TF1("fcb",fCrystalBall,
			   //0.35, mu+2*sigma, 5);
			   minx, min(maxx,mu+2*sigma), 5);
        // f(x; alpha, n, xbar, sigma)
	fcb->SetParameters(1, 1.8, 30, 1.008, 0.107);
	fcb->SetParLimits(1, 0.5, 3);
	fcb->SetLineColor(kBlack);
	hpr5->Fit(fcb,"QRN");

	// Two-sided Crystal Ball fit moved earlier

	// D0 function fit
	TF1 *fd02 =  new TF1("fd0",fD0Jet2,0.35,1.65, 7);
	// f(x; N, mu, sigma, P, lambda, P2, lambda2)
	if ( b1) fd02->SetParameters(1, 1.014, 0.107, 0.044, 9.5, 0.06, 9.5);
	if (!b1) fd02->SetParameters(0.99, 1.008, 0.043, 0.0055, 9.1);
	fd02->SetParLimits(0,0.9,1.1);
	fd02->SetParLimits(3,0.0,0.2);
	fd02->SetParLimits(5,0.0,0.2);
	fd02->SetParLimits(4,3,100);
	fd02->SetParLimits(6,3,100);
	fd02->SetLineColor(kBlack);
	hpr5->Fit(fd02,"QRN");

	hpr5->SetMarkerColor(kRed+2);
	hpr5->SetLineColor(kRed+2);
	hpr5->SetFillColor(kRed-9);
	hpr5->SetFillStyle(1001);

	hpr5->DrawClone("SAMEHIST");
	hpr5->DrawClone("SAME");

	fcb->SetRange(0.2,1.6);
	fcb->SetLineStyle(kSolid);

	fcb2->SetLineStyle(kSolid);
	fcb2->SetLineWidth(2);
	fcb2->SetRange(minx,maxx);
	fcb2->DrawClone("SAME");
	fcb2->SetRange(0.2,1.8);
	fcb2->SetLineStyle(kDotted);
	fcb2->SetLineWidth(1);
	fcb2->DrawClone("SAME");

	//fd02->DrawClone("SAME");
	fd02->SetRange(0.2,1.8);
	fd02->SetLineStyle(kDotted);
	//fd02->DrawClone("SAME");
	fd02->SetLineStyle(kSolid);
	
	fg5->SetLineStyle(kSolid);
	fg5->SetLineColor(kRed+1);
	fg5->SetLineWidth(2);
	fg5->DrawClone("SAME");
	fg5->SetLineStyle(kDotted);
	fg5->SetRange(0.5,1.5);
	fg5->DrawClone("SAME");
	fg5->SetLineStyle(kSolid);

	tex->SetTextSize(0.045);
	if (b1)
	  //tex->DrawLatex(0.20,0.86,"97 #leq p_{T}^{ptcl} < 114 GeV");
	  tex->DrawLatex(0.20,0.86,"97 #leq p_{T,ptcl} < 114 GeV"); // paper v5
	if (!b1)
	  //tex->DrawLatex(0.20,0.86,"967 #leq p_{T}^{ptcl} < 1032 GeV");
	  tex->DrawLatex(0.20,0.86,"967#leq p_{T,ptcl}< 1032 GeV"); // v5
	tex->DrawLatex(0.20,0.78,"|#eta| < 0.5, PF");
	//tex->DrawLatex(0.20,0.74,"R_{cone} = 0.7");
	//tex->DrawLatex(0.20,0.74,"R_{cone} = 0.5");
	tex->DrawLatex(0.20,0.73,"R = 0.5"); // paper v5
	//tex->DrawLatex(0.20,0.68,"#DeltaR < R_{cone}/2");
	tex->DrawLatex(0.20,0.68,"#DeltaR < R / 2"); // paper v5

	//TLegend *leg = tdrLeg(0.65,0.66,0.90,0.90);
	TLegend *leg = tdrLeg(0.63,0.72,0.90,0.90);
	//leg->AddEntry(hpr7,"R_{cone} = 0.7","FL");
	leg->AddEntry(hpr5,"MC","FLP");
	leg->AddEntry(fg5,"Gaussian","L");
	//leg->AddEntry(fcb2,"Crystal Ball^{2}","L");
	//leg->AddEntry(fd02,"D0 Jet^{2}","L");
	leg->AddEntry(fcb2,"Crystal Ball","L");
	//leg->AddEntry(fd02,"D0 Jet","L");

	//cmsPrel(0);
	gPad->RedrawAxis();
	//if ( b1) c0->SaveAs("pdf/resolution_jertails_97-114GeV.pdf");
	//if (!b1) c0->SaveAs("pdf/resolution_jertails_967-1032GeV.pdf");
	if ( b1) c0->SaveAs(Form("pdf/resolution_jertails_97-114GeV_%s.pdf",
				 cf));
	if (!b1) c0->SaveAs(Form("pdf/resolution_jertails_967-1032GeV_%s.pdf",
				 cf));
      }

    } // for ibin
    
    //c1->SaveAs(Form("pdf/resolution_Rap%d_gausfits.pdf",iy));
    c1->SaveAs(Form("pdf/resolution_Rap%d_gausfits_%s.pdf",iy,cf));

    TCanvas *c2 = new TCanvas(Form("c2_%d",iy),Form("c2_%d",iy),500,750);
    c2->Divide(1,3);

    TH1D *h2 = new TH1D(Form("h2_%d",iy),";p_{T,ptcl} (GeV);Response",
			//int(ptmax-ptmin),ptmin,ptmax);
			int(emax-ptmin),ptmin,emax);

    c2->cd(1);
    gPad->SetLogx();
    h2->SetMinimum(0.98);//0.98);
    h2->SetMaximum(1.05);//1.02);
    h2->GetXaxis()->SetMoreLogLabels();
    h2->GetXaxis()->SetNoExponent();
    h2->DrawClone("AXIS");

    TLine *l = new TLine();
    l->SetLineStyle(kDashed);
    //l->DrawLine(ptmin,1,ptmax,1);
    l->DrawLine(ptmin,1,emax,1);
    l->SetLineStyle(kDotted);
    //l->DrawLine(ptmin,1.005,ptmax,1.005);
    //l->DrawLine(ptmin,0.995,ptmax,0.995);
    l->DrawLine(ptmin,1.005,emax,1.005);
    l->DrawLine(ptmin,0.995,emax,0.995);

    gr5->SetMarkerColor(kBlue);
    gr5->Draw("SAMEP");
    gr5a->SetMarkerColor(kBlue);
    gr5a->SetMarkerStyle(kOpenCircle);
    gr5a->Draw("SAMEP");
    gr5b->SetLineColor(kBlue);
    gr5b->Draw("SAMELx");

    TLegend *leg1 = new TLegend(0.25,0.70,0.45,0.90,"","brNDC");
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(kNone);
    leg1->SetTextSize(0.045);
    leg1->AddEntry(gr5,Form("#mu(%s)",c5),"P");
    leg1->AddEntry(gr5a,Form("MEAN(%s)",c5),"P");
    leg1->Draw();

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->DrawLatex(0.50,0.85,Form("%1.1f < |#eta| < %1.1f",etamin,etamax));
    tex->SetTextColor(kBlue);
    tex->DrawLatex(0.75,0.85,cf);

    c2->cd(2);
    gPad->SetLogx();
    h2->SetMinimum(0.0);
    h2->SetMaximum(0.25);
    h2->GetXaxis()->SetMoreLogLabels();
    h2->GetXaxis()->SetNoExponent();
    h2->GetYaxis()->SetTitle("Resolution");
    h2->DrawClone("AXIS");

    gs5->SetMarkerColor(kBlue);
    gs5->Draw("SAMEP");
    gs5a->SetMarkerColor(kBlue);
    gs5a->SetMarkerStyle(kOpenCircle);
    gs5a->Draw("SAMEP");

    //TF1 *fs = new TF1("fs",ptreso,ptmin,ptmax,2);
    _usejme = false;
    //_run2012 = true; _run2018 = false;
    _jer_iov = run1;
    TF1 *fs = new TF1("fs",ptreso,ptmin,emax,2);
    //fs->SetNpx(int(emax-ptmin));
    fs->SetNpx(int((emax-ptmin)/5.));
    fs->SetParameters(etamin,1.);
    fs->SetLineColor(kGreen+2);
    fs->SetLineStyle(kDashed);
    fs->SetLineWidth(2);
    fs->DrawClone("SAME");

    // Save into TGraph for plotting ratio
    TGraph *gsrun1 = new TGraph(0);
    gsrun1->SetLineColor(fs->GetLineColor());
    gsrun1->SetLineStyle(fs->GetLineStyle());
    gsrun1->SetLineWidth(fs->GetLineWidth());
    for (int i = 0; i != gs5->GetN(); ++i) {
      gsrun1->SetPoint(i, gs5->GetX()[i], fs->Eval(gs5->GetX()[i]));
    }

    _usejme = true;
    //_run2012 = false; _run2018 = true;
    _jer_iov = jer_ref;//run2016;
    TF1 *fsjme = new TF1("fsjme",ptreso,ptmin,emax,2);
    //fsjme->SetNpx(int(emax-ptmin));
    fsjme->SetNpx(int((emax-ptmin)/5.));
    fsjme->SetParameters(etamid,1.);
    fsjme->SetLineColor(kRed+2);
    fsjme->SetLineStyle(kDotted);
    fsjme->SetLineWidth(2);
    fsjme->DrawClone("SAME");
    //fsjme->Draw("SAME");

    // Save into TGraph for plotting ratio
    TGraph *gsjme = new TGraph(0);
    gsjme->SetLineColor(fsjme->GetLineColor());
    gsjme->SetLineStyle(fsjme->GetLineStyle());
    gsjme->SetLineWidth(fsjme->GetLineWidth());
    for (int i = 0; i != gs5->GetN(); ++i) {
      gsjme->SetPoint(i, gs5->GetX()[i], fsjme->Eval(gs5->GetX()[i]));
    }

    _usejme = false;
    //_run2012 = false; _run2018 = true;
    _jer_iov = jer_ref;//run2016;
    TF1 *fsjer = new TF1("fsjer",ptreso,ptmin,emax,2);
    fsjer->SetNpx(int((emax-ptmin)/5.));
    fsjer->SetParameters(etamid,1.);
    fsjer->SetLineColor(kBlue+2);
    fsjer->SetLineStyle(kDashDotted);
    fsjer->SetLineWidth(2);
    fsjer->DrawClone("SAME");

    // Save into TGraph for plotting ratio
    TGraph *gsjer = new TGraph(0);
    gsjer->SetLineColor(fsjer->GetLineColor());
    gsjer->SetLineStyle(fsjer->GetLineStyle());
    gsjer->SetLineWidth(fsjer->GetLineWidth());
    for (int i = 0; i != gs5->GetN(); ++i) {
      gsjer->SetPoint(i, gs5->GetX()[i], fsjer->Eval(gs5->GetX()[i]));
    }

    //TF1 *fs5 = new TF1(Form("fs5_%d",iy),"sqrt([0]*[0]/(x*x)"
    TF1 *fs5 = new TF1(Form("fs5_%d",iy),"sqrt(abs([0])*[0]/(x*x)"
		       "+ [1]*[1]/x*pow(x,[2])"
		       "+ [3]*[3])",
		       //fitxmin, ptmax);
		       fitxmin, emax);
    fs5->SetParameters(1,0.8,0.2,0.04);
    fs5->FixParameter(2,0.);
    gs5->Fit(fs5,"QRN");
    fs5->SetLineColor(kBlue);
    fs5->DrawClone("SAME");
    //fs5->SetRange(ptmin,ptmax);
    fs5->SetRange(ptmin,emax);
    fs5->SetLineStyle(kDashed);
    fs5->DrawClone("SAME");

    TLegend *leg2 = new TLegend(0.20,0.20,0.40,0.45,"","brNDC");
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(kNone);
    leg2->SetTextSize(0.045);
    //leg2->AddEntry(fs,"QCD-11-004","L");
    leg2->AddEntry(gs5,Form("#sigma(%s)",c5),"P");
    leg2->AddEntry(gs5a,Form("RMS(%s)",c5),"P");
    leg2->Draw();

    tex->SetTextColor(kRed+2);
    tex->DrawLatex(0.50,0.81,Form("JME official JER for #LT#rho#GT=%1.2f,"
				  " |#eta|=%1.2f",_rho,etamid));
    tex->SetTextColor(kGreen+2);
    tex->DrawLatex(0.50,0.73,"Run 1 AK5 JER scaled to AK4");
    tex->SetTextColor(kBlue);
    tex->DrawLatex(0.50,0.65,Form("#chi^{2}/ndf = %1.1f/%d (%s) %s",
				  fs5->GetChisquare(),
				  fs5->GetNDF(), "R=0.4", cf));
    tex->DrawLatex(0.50,0.57,Form("N=%1.2f, S=%1.3f, C=%1.4f",
				  fs5->GetParameter(0),
				  fs5->GetParameter(1),
				  fs5->GetParameter(3)));

    c0b->cd(1);
    gPad->SetLogx();
    if (iy==0) h2->DrawClone("AXIS");
    TGraphErrors *gs5b = (TGraphErrors*)gs5->DrawClone("SAMEP");
    gs5b->SetLineColor(colors[iy]);
    gs5b->SetMarkerColor(colors[iy]);
    TF1 *fs5b = (TF1*)fs5->DrawClone("SAME");
    fs5b->SetLineColor(colors[iy]);

    c0b->cd(2);
    gPad->SetLogx();
    if (iy==0) h2->DrawClone("AXIS");

    c2->cd(3);
    gPad->SetLogx();
    h2->SetMinimum(0.80);
    h2->SetMaximum(1.30);//1.20);
    h2->GetYaxis()->SetTitle("Resolution / Fit");
    h2->DrawClone("AXIS");
    //l->DrawLine(ptmin,1,ptmax,1);
    l->DrawLine(ptmin,1,emax,1);
    l->SetLineStyle(kDotted);
    //l->DrawLine(ptmin,1.05,ptmax,1.05);
    //l->DrawLine(ptmin,0.95,ptmax,0.95);
    l->DrawLine(ptmin,1.05,emax,1.05);
    l->DrawLine(ptmin,0.95,emax,0.95);

    TGraphErrors *gs5ar = (TGraphErrors*)gs5a->Clone();
    for (int i = 0; i != gs5a->GetN(); ++i) {
      double sfit = fs5->Eval(gs5a->GetX()[i]);
      gs5ar->SetPoint(i, gs5a->GetX()[i], gs5a->GetY()[i]/sfit);
      gs5ar->SetPointError(i, gs5a->GetEX()[i], gs5a->GetEY()[i]/sfit);
    }
    gs5ar->Draw("SAMEP");

    TGraph *gsrun1r = (TGraph*)gsrun1->Clone();
    for (int i = 0; i != gsrun1->GetN(); ++i) {
      double sfit = fs5->Eval(gsrun1->GetX()[i]);
      gsrun1r->SetPoint(i, gsrun1->GetX()[i], gsrun1->GetY()[i]/sfit);
    }
    gsrun1r->Draw("SAMEL");

    TGraph *gsjmer = (TGraph*)gsjme->Clone();
    for (int i = 0; i != gsjme->GetN(); ++i) {
      double sfit = fs5->Eval(gsjme->GetX()[i]);
      gsjmer->SetPoint(i, gsjme->GetX()[i], gsjme->GetY()[i]/sfit);
    }
    gsjmer->Draw("SAMEL");

    TGraph *gsjerr = (TGraph*)gsjer->Clone();
    for (int i = 0; i != gsjer->GetN(); ++i) {
      double sfit = fs5->Eval(gsjer->GetX()[i]);
      gsjerr->SetPoint(i, gsjer->GetX()[i], gsjer->GetY()[i]/sfit);
    }
    gsjerr->Draw("SAMEL");


    TGraphErrors *gs5r = (TGraphErrors*)gs5->Clone();
    for (int i = 0; i != gs5->GetN(); ++i) {
      double sfit = fs5->Eval(gs5->GetX()[i]);
      gs5r->SetPoint(i, gs5->GetX()[i], gs5->GetY()[i]/sfit);
      gs5r->SetPointError(i, gs5->GetEX()[i], gs5->GetEY()[i]/sfit);
    }
    gs5r->Draw("SAMEP");

    tex->DrawLatex(0.75,0.85,cf);

    //c2->SaveAs(Form("pdf/resolution_Rap%d.pdf",iy));
    c2->SaveAs(Form("pdf/resolution_Rap%d_%s.pdf",iy,cf));


    c3->cd();
    //fs5->SetLineStyle(kSolid);
    //fs5->SetLineColor(iy+1);
    fs5->SetLineStyle(styles[iy]);
    fs5->SetLineColor(colors[iy]);
    fs5->SetLineWidth(3);
    fs5->DrawClone("SAME");
    leg3->AddEntry(fs5, Form("%1.1f<|#eta|<%1.1f",0.5*iy,0.5*(iy+1)),"L");

    vpar5[iy][0] = fs5->GetParameter(0);
    vpar5[iy][1] = fs5->GetParameter(1);
    vpar5[iy][2] = fs5->GetParameter(3);
    vchi5[iy] = fs5->GetChisquare(); 
    vndf5[iy] = fs5->GetNDF();

    c3e->cd();

    TF1 *fs5e = new TF1(Form("fs5e_%d",iy),"sqrt([0]*[0]/(x/[4]*x/[4])"
			"+ [1]*[1]/(x/[4])*pow(x/[4],[2])"
			"+ [3]*[3])",
			//fs5->GetXmin()*cosh(0.5*iy),
			//fs5->GetXmax()*cosh(0.5*iy));
			fitxmin*cosh(0.5*iy), ptmax*cosh(0.5*iy));
    fs5e->SetParameters(fs5->GetParameter(0), fs5->GetParameter(1),
			fs5->GetParameter(2), fs5->GetParameter(3),
			cosh(0.5*iy));
    fs5e->SetLineStyle(styles[iy]);
    fs5e->SetLineColor(colors[iy]);
    fs5e->SetLineWidth(3);
    fs5e->DrawClone("SAME");
    leg3e->AddEntry(fs5e, Form("%1.1f<|#eta|<%1.1f",0.5*iy,0.5*(iy+1)),"L");
  }
  
  //c0b->SaveAs("pdf/resolution_summary.pdf");
  c0b->SaveAs(Form("pdf/resolution_summary_%s.pdf",cf));

  // Print out fit parameters
  cout << "// Fit of JER for " << c5 << " file " << cf << endl;
  for (int iy = 0; iy != ny; ++iy) {
    cout << Form("  %s{%1.2f, %1.3f, %1.4f}%s // y %1.1f-%1.1f, chi2 %1.1f/%d",
		 iy==0 ? "{" : " ",
		 vpar5[iy][0], vpar5[iy][1], vpar5[iy][2],
		 iy==ny-1 ? "};" : ", ", ybins[iy],ybins[iy+1],
		 //0.5*iy, 0.5*(iy+1),
		 vchi5[iy], vndf5[iy]) << endl;
  }
}

// redoJER() is copied from qcdjet/ak7ak5resolution.C (24 April 2015 version)
// Assumes JER SF are independent of pT, and symmetric over +/-eta
void redoJER(string run="") {

  TDirectory *curdir = gDirectory;
  setTDRStyle();

  const char *cf = run.c_str();

  // Run I
  //const int nbins1 = 6;
  //const double bins1[nbins1+1] = {0,0.5,1.1,1.7,2.3,3.0,5.0};
  //const double vals1[nbins1] = {1.052, 1.057, 1.096, 1.134, 1.288, 1.288};
  //const double errs1[nbins1] = {0.063, 0.057, 0.065, 0.094, 0.200, 0.200};
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#JER_Scaling_factors_and_Unce_AN1
  const int nbins1 = 7;
  const double bins1[nbins1+1] = {0.0,0.5,1.1,1.7,2.3,2.8,3.2,5.0};
  const double vals1[nbins1] = {1.079,1.099,1.121,1.208,1.254,1.395,1.056};
  const double errs1[nbins1] =
    {1.105-1.053,1.127-1.071,1.150-1.092,1.254-1.162,1.316-1.192,1.458-1.332,
     1.247-0.865};
  //const int nbins = 6;
  //const double bins[nbins+1] = {0,0.5,1,1.5,2,2.5,3};

  // Reference binning and scale factors manually copied using
  // [cat file.txt | awk '{print $2","}'] and $4, and $6"-"$5 from
  // Autumn18_V1_MC_SF_AK4PFchs.txt
  const int nbins18 = 13;
  const double bins18[nbins18+1] =
    {0, 0.522, 0.783, 1.131, 1.305, 1.74, 1.93, 2.043, 2.322,
     2.5, 2.853, 2.964, 3.139, 5.191};
  const double vals18[nbins18+1] = 
    {1.15, 1.134, 1.102, 1.134, 1.104, 1.149, 1.148, 1.114, 1.347,
     2.137, 1.65, 1.225, 1.082};
  const double errs18[nbins18+1] =
    {1.193-1.107, 1.214-1.054, 1.154-1.05, 1.246-1.022, 1.315-0.893,
     1.308-0.99, 1.357-0.939, 1.305-0.923, 1.621-1.073,
     2.661-1.613, 2.591-0.709, 1.419-1.031, 1.28-0.884};

  // Fall17_V3_MC_SF_AK4PFchs.txt
  const double vals17[nbins18+1] =
    {1.1432, 1.1815, 1.0989, 1.1137, 1.1307, 1.1600, 1.2393, 1.2604,
     1.4085, 1.9909, 2.2923, 1.2696, 1.1542};
  const double errs17[nbins18+1] =
    {1.1654-1.1210, 1.2299-1.1332, 1.1444-1.0533, 1.2533-0.9740, 1.2778-0.9837,
     1.2576-1.0623, 1.4301-1.0484, 1.4105-1.1103, 1.6105-1.2066,
     2.5593-1.4225, 2.6665-1.9180, 1.3785-1.1607, 1.3066-1.0019};

  // Summer16_25nsV1_MC_SF_AK4PFchs.txt
  const double vals16[nbins18+1] =
    {1.1595, 1.1948, 1.1464, 1.1609, 1.1278, 1.1000, 1.1426, 1.1512,
     1.2963, 1.3418, 1.7788, 1.1869, 1.1922};
  const double errs16[nbins18+1] =
    {1.224-1.095, 1.26-1.1296, 1.2096-1.0832, 1.2634-1.0584, 1.2264-1.0292,
     1.2079-0.9921, 1.264-1.0212, 1.2652-1.0372, 1.5334-1.0592,
     1.5509-1.1327, 1.9796-1.578, 1.3112-1.0626, 1.341-1.0434};

  // Inclusive jet analysis bins
  const int nbins = 8;
  const double bins[nbins+1] = {0,0.5,1,1.5,2,2.5,3,3.2,4.7};

  int nbins0(0);
  const double *bins0(0), *vals0(0), *errs0(0);
  //jer_iov jer_ref(none);
  _ismcjer = true;
  _usejme = false;
  if (run=="Run2018") {
    _jer_iov = run2018;
    nbins0 = nbins18;
    bins0 = &bins18[0];
    vals0 = &vals18[0];
    errs0 = &errs18[0];
  }
  if (run=="Run2017") {
    _jer_iov = run2017;
    nbins0 = nbins18;
    bins0 = &bins18[0];
    vals0 = &vals17[0];
    errs0 = &errs17[0];
  }
  if (run=="Run2016") {
    _jer_iov = run2016;
    nbins0 = nbins18;
    bins0 = &bins18[0];
    vals0 = &vals16[0];
    errs0 = &errs16[0];
  }
  if (run=="Run1") {
    _jer_iov = run1;
    nbins0 = nbins1;
    bins0 = &bins1[0];
    vals0 = &vals1[0];
    errs0 = &errs1[0];
  }
  assert(bins0);

  //const double refpt = 50;
  const double refpt = 95.5; // Eta_0.0-1.3 jt60 hselpt->GetMean()

  TGraphErrors *gk = new TGraphErrors(0);
  for (int i = 0; i != nbins0; ++i) {
    gk->SetPoint(i, 0.5*(bins0[i]+bins0[i+1]), vals0[i]);
    //gk->SetPointError(i, 0.5*(bins0[i+1]-bins0[i]), errs0[i]);
    gk->SetPointError(i, 0.5*(bins0[i+1]-bins0[i]), 0.5*errs0[i]);
  }

  // Take y histogram from data and average JER SF over y bin using that
  TFile *f = new TFile("rootfiles/output-DATA-1-Fall18V8-D.root","READ");
  assert(f && !f->IsZombie());
  TH1D *heta(0);
  TGraphErrors *gjer2 = new TGraphErrors(0);
  TGraphErrors *gjer20 = new TGraphErrors(0);
  TGraphErrors *gjer3 = new TGraphErrors(0);
  TGraphErrors *gk2 = new TGraphErrors(0);
  TGraphErrors *gk20 = new TGraphErrors(0);
  for (int i = 0; i != nbins; ++i) {
    double y1 = bins[i]; double y2 = bins[i+1];
    TH1D *hy = (TH1D*)f->Get(Form("Standard/Eta_%1.1f-%1.1f/jt40/hy",y1,y2));
    assert(hy);
    curdir->cd();
    if (!heta) heta = (TH1D*)hy->Clone("heta");
    else heta->Add(hy);

    double sumv2(0), sumw2(0), sumu2(0), sumw(0);
    for (int j = 1; j != hy->GetNbinsX()+1; ++j) {
      double y = hy->GetBinCenter(j);
      double nj = hy->GetBinContent(j); 
      //if (y>0 && nj>0) {
      if (nj>0) {
	//_ak7 = false;
	//_ismcjer = true;
	//_usejme = false;
	//_jer_iov = jer_ref;
	double jermc = ptresolution(refpt, y);
	sumv2 += nj*pow(jermc,2);
	//int jj = TMath::BinarySearch(nbins0, bins0, y);
	int jj = TMath::BinarySearch(nbins0, bins0, fabs(y));
	double k = vals0[jj];
	double jerdt = k * jermc;
	sumw2 += nj*pow(jerdt,2);
	sumu2 += nj*pow(jermc*(k+0.5*errs0[jj]),2);
	sumw += nj;
      }
    }
    double jermc = sqrt(sumv2)/sqrt(sumw);
    double jerdt = sqrt(sumw2)/sqrt(sumw);
    double jerdtu = sqrt(sumu2)/sqrt(sumw);
    double errdt = jerdtu - jerdt;
    gjer2->SetPoint(i, 0.5*(y1+y2), jerdt);
    gjer2->SetPointError(i, 0.5*(y2-y1), errdt);
    gjer20->SetPoint(i, 0.5*(y1+y2), jerdt);
    gjer20->SetPointError(i, 0.5*(y2-y1), 0);
    gjer3->SetPoint(i, 0.5*(y1+y2), jermc);
    gjer3->SetPointError(i, 0.5*(y2-y1), 0.);
    gk2->SetPoint(i, 0.5*(y1+y2), jerdt / jermc);
    gk2->SetPointError(i, 0.5*(y2-y1), errdt / jermc);
    gk20->SetPoint(i, 0.5*(y1+y2), jerdt / jermc);
    gk20->SetPointError(i, 0.5*(y2-y1), 0);
  } // for i

  heta->Scale(10./heta->Integral());
  heta->GetXaxis()->SetRangeUser(0.,3.0);
  heta->SetMaximum(0.30);//0.50);
  heta->SetMinimum(0.00);//-0.10);

  TGraphErrors *gjer = new TGraphErrors(0);
  TGraphErrors *gjer1 = new TGraphErrors(0);
  for (int i = 1; i != heta->GetNbinsX()+1; ++i) {
    double x = heta->GetBinCenter(i);
    if (i>0.) {
      int n = gjer->GetN();
      //_ak7 = false;
      //_ismcjer = true;
      //_usejme = false;
      //_jer_iov = jer_ref;
      double jermc = ptresolution(refpt, x);
      gjer->SetPoint(n, x, jermc);
      gjer->SetPointError(n, 0.5*heta->GetBinWidth(i), 0.);
      int j = TMath::BinarySearch(nbins0, bins0, x);
      double k = vals0[j];
      double jerdt = k*jermc;
      gjer1->SetPoint(n, x, jerdt);
      gjer1->SetPointError(n, 0.5*heta->GetBinWidth(i), jermc*0.5*errs0[j]);
    }
  }
  gjer1->SetFillStyle(1001);
  gjer1->SetFillColor(kYellow+1);

  TH1D *h = new TH1D("h",";Jet |#eta|;Resolution (#sigma_{p_{T}} / p_{T})",
		     47,0,4.7);
  h->SetMaximum(0.3);
  h->SetMinimum(0.0);

  lumi_13TeV = run;
  lumi_8TeV = "Run1";
  TCanvas *c1 = tdrCanvas("c1", h, run=="Run1" ? 2: 4, 11, kSquare);

  tdrDraw(gjer2,"E2");
  tdrDraw(gjer20,"PZ",kNone,kYellow+3,kSolid,kYellow+3);
  gjer20->SetLineWidth(3);
  tdrDraw(gjer,"PZ",kNone);
  gjer->SetLineWidth(3);

  TLegend *leg = tdrLeg(0.18,0.65,0.48,0.77);
  leg->AddEntry(gjer20,"Data","FL");
  leg->AddEntry(gjer,"MC","L");
  leg->Draw();

  gPad->RedrawAxis();

  cout << "// JER SF produced with minitools/resolution.C:redoJER"
       << "(\""<<run<<"\")" << endl;
  for (int i = 0; i != gk2->GetN(); ++i) {

    double y1 = gk2->GetX()[i] - gk2->GetEX()[i];
    double y2 = gk2->GetX()[i] + gk2->GetEX()[i];
    double val =  gk2->GetY()[i];
    double eval = gk2->GetEY()[i];
    cout << Form(" {%1.3f, %1.3f}, // %1.1f-%1.1f",
		 val, eval, y1, y2) << endl;
  }

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);
  tex->DrawLatex(0.20,0.59,Form("p_{T} = %1.1f GeV",refpt));

  c1->SaveAs(Form("pdf/resolution_interpol_%s.pdf",cf));


  TH1D *h2 = new TH1D("h2",";Jet |#eta|;JER data/MC scale factor",47,0,4.7);
  h2->SetMaximum(2.5);
  h2->SetMinimum(0.7);

  TCanvas *c2 = tdrCanvas("c2", h2, run=="Run1" ? 2 : 4, 11, kSquare);

  tdrDraw(gk2,"E2");
  tdrDraw(gk20,"PZ",kNone,kYellow+3,kSolid,kYellow+3);
  gk20->SetLineWidth(3);
  tdrDraw(gk,"Pz");

  TLine *l = new TLine();
  l->SetLineStyle(kSolid);
  l->DrawLine(0,1,4.7,1);
  l->SetLineStyle(kDotted);
  l->DrawLine(0,1.15,4.7,1.15);
  //l->DrawLine(0,0.85,4.7,0.85);

  TLegend *leg2 = tdrLeg(0.18,0.65,0.48,0.77);
  leg2->AddEntry(gk20,"Interpolated SF","FL");
  leg2->AddEntry(gk,"Original SF","PL");
  leg2->Draw();


  //cmsPrel(800);

  gPad->RedrawAxis();

  c2->SaveAs(Form("pdf/resolution_interpol_datamc_%s.pdf",cf));
} // redoJER

// (resolution_datamc has not yet been maintained, is a copy from qcdjet)
// Purpose is to derive JER data/MC SF from SMP-J tuples,
// to compare to JME official results and double-check them 
void resolution_datamc() {

  //_ismcjer = true;

  setTDRStyle();

  TFile *fd = new TFile("GR13_AK7_20fb/output-DATA-2b.root","READ");
  assert(fd && !fd->IsZombie());
  assert(gDirectory->cd("Standard"));
  TDirectory *dind = gDirectory;

  //TFile *fm = new TFile("output-MC-2b.root","READ");
  //TFile *fm = new TFile("GRV23_AK5_42x_v2/output-MC-2b.root","READ");
  TFile *fm = new TFile("GR13_AK7_20fb/output-MC-2b.root","READ");
  assert(fm && !fm->IsZombie());
  assert(gDirectory->cd("Standard"));
  TDirectory *dinm = gDirectory;

  TH1D *h1 = new TH1D("h1",";p_{T,3}/p_{T,ave};Dijet asymmetry RMS",100,0,1.0);
  h1->SetMinimum(0.+0.0001);
  h1->SetMaximum(0.5-0.0001);
  TH1D *h1jes = new TH1D("h1jes",";p_{T,3}/p_{T,ave};Dijet asymmetry mean",
			 100,0,1.0);
  h1jes->SetMinimum(-0.3+0.0001);
  h1jes->SetMaximum(0.3-0.0001);

  TH1D *h2 = new TH1D("h2",";p_{T,3}/p_{T,ave};Resolution",100,0,0.5);
  h2->SetMinimum(-0.3+0.0001);
  h2->SetMaximum(0.3-0.0001);
  TH1D *h2jes = new TH1D("h2jes",";p_{T,3}/p_{T,ave};Response",100,0,0.5);
  //h2jes->SetMinimum(-0.3+0.0001);
  //h2jes->SetMaximum(0.3-0.0001);
  h2jes->SetMinimum(0.9+0.0001);
  h2jes->SetMaximum(1.1-0.0001);

  TH1D *h3 = new TH1D("h3",";p_{T,3}/p_{T,ave};Data / MC (JER)",100,0,0.5);
  h3->SetMinimum(0.0+0.0001);
  h3->SetMaximum(2.5-0.0001);
  TH1D *h3jes = new TH1D("h3jes",";p_{T,3}/p_{T,ave};Data / MC (JES)",
			 100,0,0.5);
  h3jes->SetMinimum(0.8+0.0001);
  h3jes->SetMaximum(1.2-0.0001);

  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextSize(0.045);

  TCanvas *c4_0 = new TCanvas("c4_0","c4_0",1200,1200);
  c4_0->SetTopMargin(0.17);
  c4_0->Divide(3,3,-1,-1);
  //cmsPrel(_lumi,true);

  TCanvas *c4jes_0 = new TCanvas("c4jes_0","c4jes_0",1200,1200);
  c4jes_0->SetTopMargin(0.17);
  c4jes_0->Divide(3,3,-1,-1);
  //cmsPrel(_lumi,true);

  TCanvas *c4fsr_0 = new TCanvas("c4fsr_0","c4fsr_0",1200,1200);
  c4fsr_0->SetTopMargin(0.17);
  c4fsr_0->Divide(3,3,-1,-1);
  //cmsPrel(_lumi,true);

  double cterm013 = 0;
  // average constant term for |y|<1.5
  const double cterm = 0;//0.0241;//0.0341;//36pb-1,Z2

  const int nybins = 8;
  const double ybins[nybins+1] = {0,0.5,1,1.5,2,2.5,3,3.2,4.7};

  for (int iy = 0; iy != nybins; ++iy) {

    const double ymin = ybins[iy];
    const double ymax = ybins[iy+1];

    const int nx=6, ny=4;
    TCanvas *c1_0 = new TCanvas(Form("c1_%d",iy),Form("c1_%d",iy),
				nx*150,ny*150);
    c1_0->Divide(nx,ny,-1,-1);
    TCanvas *c1jes_0 = new TCanvas(Form("c1jes_%d",iy),Form("c1jes_%d",iy),
				nx*150,ny*150);
    c1jes_0->Divide(nx,ny,-1,-1);
    TCanvas *c2_0 = new TCanvas(Form("c2_%d",iy),Form("c2_%d",iy),
				nx*150,ny*150);
    c2_0->Divide(nx,ny,-1,-1);
    TCanvas *c2jes_0 = new TCanvas(Form("c2jes_%d",iy),Form("c2jes_%d",iy),
				nx*150,ny*150);
    c2jes_0->Divide(nx,ny,-1,-1);
    TCanvas *c3_0 = new TCanvas(Form("c3_%d",iy),Form("c3_%d",iy),
				nx*150,ny*150);
    c3_0->Divide(nx,ny,-1,-1);
    TCanvas *c3jes_0 = new TCanvas(Form("c3jes_%d",iy),Form("c3jes_%d",iy),
				   nx*150,ny*150);
    c3jes_0->Divide(nx,ny,-1,-1);

    assert(dind->cd(Form("Eta_%1.1f-%1.1f",ymin,ymax)));
    TDirectory *dd = gDirectory;
    TH3D *hd3 = (TH3D*)dd->Get("hdjasymm"); assert(hd3);

    assert(dinm->cd(Form("Eta_%1.1f-%1.1f",ymin,ymax)));
    TDirectory *dm = gDirectory;
    TH3D *hm3 = (TH3D*)dm->Get("hdjasymm"); assert(hm3);

    TGraphErrors *gr = new TGraphErrors(0);
    TGraphErrors *grjes = new TGraphErrors(0);
    TGraphErrors *grfsr = new TGraphErrors(0);

    const int nskip = 12;//10;
    const int ndrop = 1;
    for (int ipt = 1+nskip; ipt != hd3->GetNbinsX()+1-ndrop; ++ipt) {

      // Project RMS vs pT3,rd
      TH1D *hd = hd3->ProjectionY(Form("hd_y%d_pt%d",iy,ipt));
      TH1D *hdjes = hd3->ProjectionY(Form("hdjes_y%d_pt%d",iy,ipt));
      for (int i = 0; i != hd->GetNbinsX()+1; ++i) {
	TH1D *hd1 = hd3->ProjectionZ("tmp",ipt,ipt,i,i);
	hd->SetBinContent(i, hd1->GetRMS()/(1+hd1->GetMean()));
	hd->SetBinError(i, hd1->GetRMSError()/(1+hd1->GetMean()));
	//
	hdjes->SetBinContent(i, hd1->GetMean());
	hdjes->SetBinError(i, hd1->GetMeanError());
	delete hd1;
      } // for i
      //
      TH1D *hm = hm3->ProjectionY(Form("hm_y%d_pt%d",iy,ipt));
      TH1D *hmjes = hm3->ProjectionY(Form("hmjes_y%d_pt%d",iy,ipt));
      for (int i = 0; i != hm->GetNbinsX()+1; ++i) {
	TH1D *hm1 = hm3->ProjectionZ("tmp",ipt,ipt,i,i);
	hm->SetBinContent(i, hm1->GetRMS()/(1+hm1->GetMean()));
	hm->SetBinError(i, hm1->GetRMSError()/(1+hm1->GetMean()));
	//
	hmjes->SetBinContent(i, hm1->GetMean());
	hmjes->SetBinError(i, hm1->GetMeanError());
	delete hm1;
      } // for i
      
      // Move the 0-bin closer to true value (was to the first empty)
      double pt = hd3->GetXaxis()->GetBinCenter(ipt);
      double pt1 = hd3->GetXaxis()->GetBinLowEdge(ipt);
      double pt2 = hd3->GetXaxis()->GetBinLowEdge(ipt+1);
      double ptmax = hd3->GetXaxis()->GetBinLowEdge(ipt+1);
      int k = hd->FindBin(_recopt/ptmax)-1;
      int k2 = k;
      if (hd->GetBinContent(k2)==0) {
	hd->SetBinContent(k2, hd->GetBinContent(1));
	hd->SetBinError(k2, hd->GetBinError(1));
	hd->SetBinContent(1, 0.);
	hd->SetBinError(1, 0.);
      }
      //
      if (hdjes->GetBinContent(k2)==0) {
	hdjes->SetBinContent(k2, hdjes->GetBinContent(1));
	hdjes->SetBinError(k2, hdjes->GetBinError(1));
	hdjes->SetBinContent(1, 0.);
	hdjes->SetBinError(1, 0.);
      }
      //
      k = hm->FindBin(_recopt/ptmax)-1;
      k2 = k;
      if (hm->GetBinContent(k2)==0) {
	hm->SetBinContent(k2, hm->GetBinContent(1));
	hm->SetBinError(k2, hm->GetBinError(1));
	hm->SetBinContent(1, 0.);
	hm->SetBinError(1, 0.);
      }
      //
      if (hmjes->GetBinContent(k2)==0) {
	hmjes->SetBinContent(k2, hmjes->GetBinContent(1));
	hmjes->SetBinError(k2, hmjes->GetBinError(1));
	hmjes->SetBinContent(1, 0.);
	hmjes->SetBinError(1, 0.);
      }
	
      // Figure out sensitivity to intrinsic resolution
      TF1 *f1 = new TF1(Form("f1_y%d_pt%d",iy,ipt),
			"sqrt([0]*[0]+[1]*[1]*x*x)",0.0,0.5);
      double res0 = sqrt(1/3.*pow(ptresolution(pt, 0.0+0.001),2)
			 +1/3.*pow(ptresolution(pt, 0.5+0.001),2)
			 +1/3.*pow(ptresolution(pt, 1.0+0.001),2));
      double res = sqrt(pow(ptresolution(pt, ymin+0.001),2) + res0*res0)/2.;
      f1->SetParameters(res,1.);
      // Adapt free parameters and fit range to gain sensitivity
      f1->SetRange(0., _recopt/ptmax>0.3 ? 0.5 : 0.4);
      if (_recopt/ptmax>0.15) f1->FixParameter(0, res);
      // Fit with JER, FSR to extract resolution
      hm->Fit(f1,"QRN");

      // Figure out sensitivity to intrinsic bias
      TF1 *f1jes = new TF1(Form("f1jes_y%d_pt%d",iy,ipt),
			   //"[0]+[1]*x",0.0,0.5);
			   "[0]+[1]*x+[1]*x*x",0.0,0.5);
      f1jes->SetParameters(0.,0.,0.);
      // Adapt free parameters and fit range to gain sensitivity
      f1jes->SetRange(0., _recopt/ptmax>0.3 ? 0.5 : 0.4);
      if (_recopt/ptmax>0.15) {
	f1jes->FixParameter(0, 0.);
	f1jes->FixParameter(1, 0.);
      }
      // Fit
      hmjes->Fit(f1jes,"QRN");

      TH1D *hdi = (TH1D*)hd->Clone(Form("hdi_y%d_pt%d",iy,ipt));
      for (int i = 1; i != hd->GetNbinsX()+1; ++i) {
	if (hd->GetBinContent(i)!=0) {
	  double bal = f1->GetParameter(1)*hd->GetBinCenter(i);
	  double y = hd->GetBinContent(i);
	  double in = sqrt(fabs(y*y - bal*bal));
	  double in2 = (4*in*in - res0*res0 - cterm*cterm);
	  double sigma = TMath::Sign(1.,in2)*sqrt(fabs(in2));
	  if (sigma>0) {
	    hdi->SetBinContent(i, sigma);
	    hdi->SetBinError(i, hd->GetBinError(i)*y/sigma*4);
	  }
	  else {
	    hdi->SetBinContent(i, 0);
	    hdi->SetBinError(i, 0);
	  }
	}
      } // for i
      //
      TH1D *hmi = (TH1D*)hm->Clone(Form("hmi_y%d_pt%d",iy,ipt));
      for (int i = 1; i != hm->GetNbinsX()+1; ++i) {
	if (hm->GetBinContent(i)!=0) {
	  double bal = f1->GetParameter(1)*hm->GetBinCenter(i);
	  double y = hm->GetBinContent(i);
	  double in = sqrt(fabs(y*y - bal*bal));
	  double in2 = fabs(4*in*in - res0*res0);
	  double sigma = TMath::Sign(1.,in2)*sqrt(fabs(in2));
	  if (sigma>0) {
	    hmi->SetBinContent(i, sigma);
	    hmi->SetBinError(i, hm->GetBinError(i)*y/sigma*4);
	  }
	  else {
	    hmi->SetBinContent(i, 0);
	    hmi->SetBinError(i, 0);
	  }
	}
      } // for i

      TH1D *hdijes = (TH1D*)hdjes->Clone(Form("hdijes_y%d_pt%d",iy,ipt));
      for (int i = 1; i != hdjes->GetNbinsX()+1; ++i) {
	if (hdjes->GetBinContent(i)!=0) {
	  double jecbias = f1jes->GetParameter(1)*hdjes->GetBinCenter(i);
	  double jec = 1 + hdjes->GetBinContent(i) - jecbias;
	  hdijes->SetBinContent(i, jec);
	  hdijes->SetBinError(i, sqrt(pow(f1jes->GetParError(0),2)+
				      pow(f1jes->GetParError(1)
					  *hdjes->GetBinCenter(i),2)));
	}
      } // for i
      //
      TH1D *hmijes = (TH1D*)hmjes->Clone(Form("hmijes_y%d_pt%d",iy,ipt));
      for (int i = 1; i != hmjes->GetNbinsX()+1; ++i) {
	if (hmjes->GetBinContent(i)!=0) {
	  double jecbias = f1jes->GetParameter(1)*hmjes->GetBinCenter(i);
	  double jec = 1 + hmjes->GetBinContent(i) - jecbias;
	  hmijes->SetBinContent(i, jec);
	  hmijes->SetBinError(i, sqrt(pow(f1jes->GetParError(0),2)+
				      pow(f1jes->GetParError(1)
					  *hmjes->GetBinCenter(i),2)));
	}
      } // for i


      // RMS(asymmetry) vs pT3/pTave
      c1_0->cd(ipt-nskip);
      h1->Draw("AXIS");
      hm->SetMarkerSize(0.5);
      hm->SetMarkerStyle(kFullSquare);
      hm->SetMarkerColor(kBlue);
      hm->SetLineColor(kBlue);
      hm->Draw("SAMEP");
      hd->SetMarkerSize(0.5);
      hd->SetMarkerStyle(kFullCircle);
      hd->SetMarkerColor(kRed);
      hd->SetLineColor(kRed);
      hd->Draw("SAMEP");
      f1->SetLineColor(kGreen+2);
      f1->SetLineWidth(1);
      f1->DrawClone("SAME");
      
      t->DrawLatex(0.20,0.85,Form("%1.1f#leq|#eta|<%1.1f,"
				  " %1.0f#leqp_{T,ave}<%1.0f GeV",
				  ymin, ymax, pt1, pt2));

      // Mean(asymmetry) vs pT3/pTave
      c1jes_0->cd(ipt-nskip);
      h1jes->Draw("AXIS");
      hmjes->SetMarkerSize(0.5);
      hmjes->SetMarkerStyle(kFullSquare);
      hmjes->SetMarkerColor(kBlue);
      hmjes->SetLineColor(kBlue);
      hmjes->Draw("SAMEP");
      hdjes->SetMarkerSize(0.5);
      hdjes->SetMarkerStyle(kFullCircle);
      hdjes->SetMarkerColor(kRed);
      hdjes->SetLineColor(kRed);
      hdjes->Draw("SAMEP");
      f1jes->SetLineColor(kGreen+2);
      f1jes->SetLineWidth(1);
      f1jes->DrawClone("SAME");
      
      t->DrawLatex(0.20,0.85,Form("%1.1f#leq|#eta|<%1.1f,"
				  " %1.0f#leqp_{T,ave}<%1.0f GeV",
				  ymin, ymax, pt1, pt2));

      // sigma=sqrt(RMS^2 - (FSR*pT3/pTave)^2) vs pT3/pTave
      c2_0->cd(ipt-nskip);
      h2->Draw("AXIS");
      hmi->GetXaxis()->SetRangeUser(0.,0.5);
      hmi->SetMarkerSize(0.5);
      hmi->SetMarkerStyle(kFullSquare);
      hmi->SetMarkerColor(kBlue);
      hmi->SetLineColor(kBlue);
      hmi->Draw("SAMEP");
      hdi->SetMarkerSize(0.5);
      hdi->SetMarkerStyle(kFullCircle);
      hdi->SetMarkerColor(kRed);
      hdi->SetLineColor(kRed);
      hdi->Draw("SAMEP");
      
      f1->SetLineWidth(1);      
      f1->ReleaseParameter(0);
      f1->FixParameter(1,0.);
      hmi->Fit(f1,"QRN");
      f1->SetLineColor(kBlue+2);
      f1->DrawClone("SAME");
      hdi->Fit(f1,"QRN");
      f1->SetLineColor(kRed+2);
      f1->SetLineWidth(1);
      f1->DrawClone("SAME");
      
      TH1D *hr = (TH1D*)hd->Clone(Form("hr_y%d_pt%d",iy,ipt));
      hr->Divide(hd, hm, 1, 1);
      TH1D *hri = (TH1D*)hdi->Clone(Form("hri_y%d_pt%d",iy,ipt));
      hri->Divide(hdi, hmi, 1, 1);

      t->DrawLatex(0.20,0.85,Form("%1.1f#leq|#eta|<%1.1f,"
				  " %1.0f#leqp_{T,ave}<%1.0f GeV",
				  ymin, ymax, pt1, pt2));

      // JEC=sqrt(Mean^2 - FSR*pT3/pTave) vs pT3/pTave
      c2jes_0->cd(ipt-nskip);
      h2jes->Draw("AXIS");
      hmijes->GetXaxis()->SetRangeUser(0.,0.5);
      hmijes->SetMarkerSize(0.5);
      hmijes->SetMarkerStyle(kFullSquare);
      hmijes->SetMarkerColor(kBlue);
      hmijes->SetLineColor(kBlue);
      hmijes->Draw("SAMEP");
      hdijes->SetMarkerSize(0.5);
      hdijes->SetMarkerStyle(kFullCircle);
      hdijes->SetMarkerColor(kRed);
      hdijes->SetLineColor(kRed);
      hdijes->Draw("SAMEP");
      
      f1jes->SetLineWidth(1);      
      f1jes->ReleaseParameter(0);
      f1jes->FixParameter(1,0.);
      hmijes->Fit(f1jes,"QRN");
      f1jes->SetLineColor(kBlue+2);
      f1jes->DrawClone("SAME");
      hdijes->Fit(f1jes,"QRN");
      f1jes->SetLineColor(kRed+2);
      f1jes->SetLineWidth(1);
      f1jes->DrawClone("SAME");
      
      TH1D *hrjes = (TH1D*)hd->Clone(Form("hrjes_y%d_pt%d",iy,ipt));
      hrjes->Divide(hdjes, hmjes, 1, 1);
      TH1D *hrijes = (TH1D*)hdijes->Clone(Form("hrijes_y%d_pt%d",iy,ipt));
      hrijes->Divide(hdijes, hmijes, 1, 1);

      t->DrawLatex(0.20,0.85,Form("%1.1f#leq|#eta|<%1.1f,"
				  " %1.0f#leqp_{T,ave}<%1.0f GeV",
				  ymin, ymax, pt1, pt2));


      // sigma_data / sigma_MC vs pT3/pTave
      c3_0->cd(ipt-nskip);
      h3->Draw("AXIS");
      hri->GetXaxis()->SetRangeUser(0.,0.5);
      hri->SetMarkerSize(0.5);
      hri->SetMarkerStyle(kFullCircle);
      hri->Draw("SAMEP");

      TF1 *fp0 = new TF1(Form("fp0_y%d_pt%d",iy,ipt),"[0]+[1]*x",0.,0.4);
      if (_recopt/ptmax>0.3) fp0->SetRange(0.,0.5);
      hri->Fit(fp0,"QRN");
      fp0->SetLineColor(kBlue);
      fp0->SetLineWidth(1);
      fp0->DrawClone("SAME");
      fp0->FixParameter(1,0.);
      hri->Fit(fp0,"QRN");
      fp0->SetLineColor(kBlack);
      fp0->SetLineWidth(2);
      fp0->DrawClone("SAME");

      if (fp0->GetParError(0)<0.30) {
	int n = gr->GetN();
	gr->SetPoint(n, pt, fp0->GetParameter(0));
	gr->SetPointError(n, 0., fp0->GetParError(0));
      }

      t->DrawLatex(0.20,0.85,Form("%1.1f#leq|#eta|<%1.1f,"
				  " %1.0f#leqp_{T,ave}<%1.0f GeV",
				  ymin, ymax, pt1, pt2));


      // JES_data / JES_MC vs pT3/pTave
      c3jes_0->cd(ipt-nskip);
      h3jes->Draw("AXIS");
      hrijes->GetXaxis()->SetRangeUser(0.,0.5);
      hrijes->SetMarkerSize(0.5);
      hrijes->SetMarkerStyle(kFullCircle);
      hrijes->Draw("SAMEP");

      TF1 *fp0jes = new TF1(Form("fp0jes_y%d_pt%d",iy,ipt),"[0]+[1]*x",0.,0.4);
      if (_recopt/ptmax>0.3) fp0jes->SetRange(0.,0.5);
      fp0jes->SetParameters(1,0.);

      fp0jes->FixParameter(1,0.);
      hrijes->Fit(fp0jes,"QRN");
      fp0jes->SetLineColor(kBlack);
      fp0jes->SetLineWidth(2);
      fp0jes->DrawClone("SAME");

      fp0jes->ReleaseParameter(1);
      hrijes->Fit(fp0jes,"QRN");
      fp0jes->SetLineColor(kBlue);
      fp0jes->SetLineWidth(1);
      fp0jes->DrawClone("SAME");


      {//if (fp0jes->GetParError(0)<0.30) {
	int n = grjes->GetN();
	grjes->SetPoint(n, pt, fp0jes->GetParameter(0));
	grjes->SetPointError(n, 0., fp0jes->GetParError(0));
      }
      {//if (fp0jes->GetParError(1)<0.30) {
	int n = grfsr->GetN();
	grfsr->SetPoint(n, pt, fp0jes->GetParameter(1));
	grfsr->SetPointError(n, 0., fp0jes->GetParError(1));
      }

      t->DrawLatex(0.20,0.85,Form("%1.1f#leq|#eta|<%1.1f,"
				  " %1.0f#leqp_{T,ave}<%1.0f GeV",
				  ymin, ymax, pt1, pt2));
      
    } // for ipt
  
    c4_0->cd(iy+1);
    gPad->SetLogx();

    int xmin = 37.;
    int xmax = 2500;
    TH1D *h4 = new TH1D(Form("h4_y%d",iy),";p_{T} (GeV);Data / MC (JER)",
			int(xmax-xmin),xmin,xmax);
    h4->GetXaxis()->SetMoreLogLabels();
    h4->GetXaxis()->SetNoExponent();
    h4->GetYaxis()->SetRangeUser(0.3,1.7999);
    h4->Draw("AXIS");

    TLine *l = new TLine();
    l->SetLineStyle(kDashed);
    l->DrawLine(xmin, 1, xmax, 1);

    // Clean out points where mean+stat<1
    for (int i = gr->GetN()-1; i != -1; --i) {
      if (gr->GetY()[i]+gr->GetEY()[i]<1.00)
	;//gr->RemovePoint(i);
    }

    gr->SetMarkerSize(0.7);
    gr->SetMarkerStyle(kFullCircle);
    gr->SetMarkerColor(kBlue);
    gr->SetLineColor(kBlue);
    gr->Draw("SAMEPz");

    TF1 *fk = new TF1(Form("fk_y%d",iy),"[0]",xmin,xmax);
    fk->SetLineColor(kBlue);
    fk->SetParameter(0,1.);
    gr->Fit(fk,"QRN");
    fk->Draw("SAME");
    

    TF1 *fnsc = new TF1(Form("fnsc_y%d",iy),ptresPlusNSC,xmin,xmax,4);
    fnsc->SetLineColor(kBlack);
    fnsc->SetParameters(0.,0.,0.,ymin+0.001);
    fnsc->FixParameter(3,ymin+0.001);
    gr->Fit(fnsc,"QRN");
    fnsc->Draw("SAME");

    TF1 *fn = new TF1(Form("fn_y%d",iy),ptresPlusN,xmin,xmax,2);
    fn->SetLineColor(kMagenta+1);
    fn->SetParameters(0.,ymin+0.001);
    fn->FixParameter(1,ymin+0.001);
    gr->Fit(fn,"QRN");
    fn->Draw("SAME");

    TF1 *fs = new TF1(Form("fs_y%d",iy),ptresPlusS,xmin,xmax,2);
    fs->SetLineColor(kGreen+2);
    fs->SetParameters(0.,ymin+0.001);
    fs->FixParameter(1,ymin+0.001);
    gr->Fit(fs,"QRN");
    fs->Draw("SAME");

    TF1 *fc = new TF1(Form("fc_y%d",iy),ptresPlusC,xmin,xmax,2);
    fc->SetLineColor(kRed);
    fc->SetParameters(0.,ymin+0.001);
    fc->FixParameter(1,ymin+0.001);
    gr->Fit(fc,"QRN");
    fc->Draw("SAME");


    cout << Form("k=%1.3g+/-%1.3g (chi2/NDF=%1.1f/%d) / "
		 "N=%1.2g, S=%1.2g, C=%1.2g (chi2/NDF=%1.1f/%d) / "
		 "N=%1.2g+/-%1.2g (chi2/NDF=%1.1f/%d) / "
		 "S=%1.2g+/-%1.2g (chi2/NDF=%1.1f/%d) / "
		 "C=%1.2g+/-%1.2g (chi2/NDF=%1.1f/%d)",
		 fk->GetParameter(0), fk->GetParError(0),
		 fk->GetChisquare(), fk->GetNDF(),
		 fnsc->GetParameter(0),fnsc->GetParameter(1),
		 fnsc->GetParameter(2),fn->GetChisquare(), fn->GetNDF(),
		 fn->GetParameter(0), fn->GetParError(0),
		 fn->GetChisquare(), fn->GetNDF(),
		 fs->GetParameter(0), fs->GetParError(0),
		 fs->GetChisquare(), fs->GetNDF(),
		 fc->GetParameter(0), fc->GetParError(0),
		 fc->GetChisquare(), fc->GetNDF()) << endl;
    if (iy<2) cterm013 += pow(fc->GetParameter(0),2);
    if (iy==2) cterm013 += 3./5.*pow(fc->GetParameter(0),2);
  
    double x0 = (iy%3==0 ? 0.20 : 0.03);
    double y0 = (iy/3==2 ? 0.20 : 0.03);
    t->SetTextColor(kBlue);
    t->DrawLatex(x0,y0+0.20,Form("k=%1.3f#pm%1.3f (#chi^{2}/NDF = %1.1f/%d)",
				fk->GetParameter(0), fk->GetParError(0),
				fk->GetChisquare(), fk->GetNDF()));
    t->SetTextColor(kMagenta+1);
    t->DrawLatex(x0,y0+0.15,Form("N=%1.1f#pm%1.1f (#chi^{2}/NDF = %1.1f/%d)",
				fn->GetParameter(0), fn->GetParError(0),
				fn->GetChisquare(), fn->GetNDF()));
    t->SetTextColor(kGreen+1);
    t->DrawLatex(x0,y0+0.10,Form("S=%1.2f#pm%1.2f (#chi^{2}/NDF = %1.1f/%d)",
				 fs->GetParameter(0), fs->GetParError(0),
				 fs->GetChisquare(), fs->GetNDF()));
    t->SetTextColor(kRed);
    t->DrawLatex(x0,y0+0.05,Form("C=%1.3f#pm%1.3f (#chi^{2}/NDF = %1.1f/%d)",
				fc->GetParameter(0), fc->GetParError(0),
				fc->GetChisquare(), fc->GetNDF()));
    t->SetTextColor(kBlack);
    t->DrawLatex(x0,y0+0.00,Form("N=%1.1f, S=%1.2f, C=%1.3f"
				 " (#chi^{2}/NDF = %1.1f/%d)",
				 fnsc->GetParameter(0), fnsc->GetParameter(1),
				 fnsc->GetParameter(2),
				 fnsc->GetChisquare(), fnsc->GetNDF()));

    t->SetTextColor(kBlack);
    double y1 = (iy/3==2 ? 0.95 : 0.92);
    t->DrawLatex(x0,y1,Form("%1.2g #leq |y| < %1.2g",ymin,ymax));
    t->DrawLatex(x0,y1-0.05,"Data vs Pythia Z2");

    //if (iy==0) {
    c1_0->cd(0);
    //c1_0->SaveAs(Form("eps/resolution_datamc_c1_%d.eps",iy));
    c1_0->SaveAs(Form("pdf/resolution_datamc_c1_%d.pdf",iy));
    c2_0->cd(0);
    //c2_0->SaveAs(Form("eps/resolution_datamc_c2_%d.eps",iy));
    c2_0->SaveAs(Form("pdf/resolution_datamc_c2_%d.pdf",iy));
    c3_0->cd(0);
    //c3_0->SaveAs(Form("eps/resolution_datamc_c3_%d.eps",iy));
    c3_0->SaveAs(Form("pdf/resolution_datamc_c3_%d.pdf",iy));
    //}
    if (_closejer) {
      delete c1_0;
      delete c2_0;
      delete c3_0;
    }


    c4jes_0->cd(iy+1);
    gPad->SetLogx();

    TH1D *h4jes = new TH1D(Form("h4jes_y%d",iy),";p_{T} (GeV);Data / MC (JES)",
			int(xmax-xmin),xmin,xmax);
    h4jes->GetXaxis()->SetMoreLogLabels();
    h4jes->GetXaxis()->SetNoExponent();
    h4jes->GetYaxis()->SetRangeUser(0.8+0.0001,1.2-0.0001);
    h4jes->Draw("AXIS");

    l->SetLineStyle(kDashed);
    l->DrawLine(xmin, 1, xmax, 1);

    // Clean out points where mean+stat<1
    for (int i = grjes->GetN()-1; i != -1; --i) {
      if (grjes->GetY()[i]+grjes->GetEY()[i]<1.00)
	;//grjes->RemovePoint(i);
    }

    grjes->SetMarkerSize(0.7);
    grjes->SetMarkerStyle(kFullCircle);
    grjes->SetMarkerColor(kBlue);
    grjes->SetLineColor(kBlue);
    //grjes->Draw("SAMEPz");
    grjes->DrawClone("SAMEPz");

    t->SetTextColor(kBlack);
    //double y1 = (iy/3==2 ? 0.95 : 0.92);
    t->DrawLatex(x0,y1,Form("%1.2g #leq |y| < %1.2g",ymin,ymax));
    t->DrawLatex(x0,y1-0.05,"Data vs Pythia Z2");


    c4fsr_0->cd(iy+1);
    gPad->SetLogx();

    TH1D *h4fsr = new TH1D(Form("h4fsr_y%d",iy),";p_{T} (GeV);Data / MC (FSR)",
			int(xmax-xmin),xmin,xmax);
    h4fsr->GetXaxis()->SetMoreLogLabels();
    h4fsr->GetXaxis()->SetNoExponent();
    h4fsr->GetYaxis()->SetRangeUser(-0.2+0.0001,0.2-0.0001);
    h4fsr->Draw("AXIS");

    l->SetLineStyle(kDashed);
    l->DrawLine(xmin, 0, xmax, 0);

    // Clean out points where mean+stat<1
    for (int i = grfsr->GetN()-1; i != -1; --i) {
      if (grfsr->GetY()[i]+grfsr->GetEY()[i]<1.00)
	;//grfsr->RemovePoint(i);
    }

    grfsr->SetMarkerSize(0.7);
    grfsr->SetMarkerStyle(kFullCircle);
    grfsr->SetMarkerColor(kBlue);
    grfsr->SetLineColor(kBlue);
    grfsr->Draw("SAMEPz");

    t->SetTextColor(kBlack);
    //double y1 = (iy/3==2 ? 0.95 : 0.92);
    t->DrawLatex(x0,y1,Form("%1.2g #leq |y| < %1.2g",ymin,ymax));
    t->DrawLatex(x0,y1-0.05,"Data vs Pythia Z2");

    c1jes_0->cd(0);
    c1jes_0->SaveAs(Form("pdf/response_datamc_c1_%d.pdf",iy));
    c2jes_0->cd(0);
    c2jes_0->SaveAs(Form("pdf/response_datamc_c2_%d.pdf",iy));
    c3jes_0->cd(0);
    c3jes_0->SaveAs(Form("pdf/response_datamc_c3_%d.pdf",iy));
  } // for iy

  cout << Form("Input cterm = C=%1.4f",cterm) << endl;
  cterm013 = sqrt(cterm013/(13./5.));
  cout << Form("Averaged out C=%1.4f (|y|<1.3)",cterm013) << endl;
  double ctermnew = sqrt((cterm013*cterm013+cterm*cterm)/2.);
  cout << Form("New input C=%1.4f (|y|<1.3)",ctermnew) << endl;

  c4_0->cd(0);
  c4_0->SaveAs("eps/resolution_datamc.eps");
  c4_0->SaveAs("pdf/resolution_datamc.pdf");
  if (_closejer) {
    delete c4_0;
  }

  c4jes_0->cd(0);
  c4jes_0->SaveAs("pdf/resolution_response_datamc.pdf");
  c4fsr_0->cd(0);
  c4fsr_0->SaveAs("pdf/resolution_responsebias_datamc.pdf");
}

