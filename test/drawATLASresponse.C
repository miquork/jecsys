#include "TCanvas.h"
#include "TH1D.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLine.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "tdrstyle_mod14.C"

#include <vector>

using namespace std;

// Move to mu-based mapping, which is better for comparing
// different PU scenarios as it considers both IT and OOT PU,
// plus we have a number directly comparable to ATLAS
double rhoFromMu(double mu) {
  // Eta_0.0-1.3, jt320
  return (1.01272 + 0.551183*mu + 0.000362936*mu*mu);
}

// Use this for solving ptmeas from ptgen = JEC(ptmeas) * ptmeas
FactorizedJetCorrector *_jec(0);
Double_t fJECPt(Double_t *x, Double_t *p) {

  double ptmeas = x[0];
  double eta = p[0];
  double jeta = p[1];
  double rho = p[2];

  _jec->setJetPt(ptmeas);
  _jec->setJetEta(eta);
  _jec->setJetA(jeta);
  _jec->setRho(rho);

  double jec = _jec->getCorrection();
  
  return (ptmeas * jec);
}

// Response as a function of pTprime (instead of JEC vs pTmeas)
TF1 *_jecpt(0);
double getResp(double ptgen, double eta, double jeta, double mu) {

  _jecpt->SetParameters(eta, jeta, rhoFromMu(mu));
  double ptmeas = _jecpt->GetX(ptgen, 1, 4000);
  double resp = ptmeas / _jecpt->Eval(ptmeas); // 1/jec

  return resp;
}

void drawATLASresponse() {

  // New style settings
  // #include "tdrstyle_mod14.C"
  setTDRStyle();
  // call cmsPrel(iPeriod, iPos);
  // int iPeriod = 2;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV 
  // second parameter in example_plot is iPos, which drives the position of the CMS logo in the plot
  // iPos=11 : top-left, left-aligned
  // iPos=33 : top-right, right-aligned
  // iPos=22 : center, centered
  // mode generally : 
  //   iPos = 10*(alignement 1/2/3) + position (1/2/3 = left/center/right)

  if (!_jec) {

    const char *sd = "CondFormats/JetMETObjects/data";
    const char *st = "Winter14_V5_MC";
    const char *s;

    //s = Form("%s/%s_L1FastJet_AK5PFchs.txt",sd,st); cout << s << endl;
    //JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
    s = Form("%s/%s_L2Relative_AK5PFchs.txt",sd,st); cout << s << endl;
    JetCorrectorParameters *l2 = new JetCorrectorParameters(s);
    s = Form("%s/%s_L3Absolute_AK5PFchs.txt",sd,st); cout << s << endl;
    JetCorrectorParameters *l3 = new JetCorrectorParameters(s);

    vector<JetCorrectorParameters> v;
    //v.push_back(l1);
    v.push_back(*l2);
    v.push_back(*l3);
    _jec = new FactorizedJetCorrector(v);
  }
  if (!_jecpt) {
    _jecpt = new TF1("jecpt",fJECPt,0,4000,3);
  }

  double ergs[] = {30, 60, 110, 400, 2000};
  const int ne = sizeof(ergs)/sizeof(ergs[0]);
  const int neta = 48;//52;
  const int jeta = TMath::Pi()*0.5*0.5;
  const int mu = 0;

  TGraph *gs[ne];
  for (int ie = 0; ie != ne; ++ie) {

    double energy = ergs[ie];

    TGraph *g = new TGraph(0); gs[ie] = g;
    for (int ieta = 0; ieta != neta; ++ieta) {
      
      double eta = (ieta+0.5)*0.1;
      double pt = energy / cosh(eta);
      if (pt > 10.) {

	double jes = getResp(pt, eta, jeta, mu);
	int n = g->GetN();
	g->SetPoint(n, eta, jes);
      }
    } // for ie
  } // for ieta


  // Draw results
  //TCanvas *c1 = new TCanvas("c1","c1",600,600);
  //TCanvas *c1 = new TCanvas("c1","c1",800,600); // ATLAS shape
  TH1D *h = new TH1D("h",";Jet |#eta|;Jet response at PF scale",40,0,4.8);
  h->SetMaximum(1.25);
  h->SetMinimum(0.5);
  //h->Draw("AXIS");
  TCanvas *c1 = tdrCanvas("c1",h,2,0);

  TLegend *leg1 = tdrLeg(0.25,0.25,0.55,0.30);
  TLegend *leg2 = tdrLeg(0.25,0.20,0.55,0.25);
  TLegend *leg3 = tdrLeg(0.25,0.15,0.55,0.20);
  TLegend *leg4 = tdrLeg(0.55,0.25,0.85,0.30);
  TLegend *leg5 = tdrLeg(0.55,0.20,0.85,0.25);
  TLegend *legs[ne] = {leg1, leg2, leg3, leg4, leg5};

  int colors[] = {kGreen+2, kBlack, kOrange+1, kBlue, kRed+1};
  int markers[] = {kFullCircle, kOpenCircle, kFullSquare, kOpenSquare,
		   kFullTriangleUp};

  for (int ie = 0; ie != ne; ++ie) {
    
    TGraph *g = gs[ie];
    g->SetMarkerColor(colors[ie]);
    g->SetMarkerStyle(markers[ie]);
    g->Draw("SAMEP");

    //TLegend *leg = (ie<3 ? leg1 : leg2);
    TLegend *leg = legs[ie];
    leg->SetTextColor(colors[ie]);
    leg->AddEntry(g, Form("E = %1.0f GeV",ergs[ie]), "P");
  }


  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);
  
  TLine *l = new TLine();
  l->DrawLine(1.3,0.7,1.3,1.1);
  l->DrawLine(2.5,0.7,2.5,1.1);
  l->DrawLine(3.0,0.7,3.0,1.1);
  l->DrawLine(4.5,0.7,4.5,1.1);
  l->SetLineStyle(kDashed);
  l->DrawLine(3.2,0.7,3.2,1.1);

  tex->DrawLatex(0.35,0.86,"2012 JES: Anti-k_{t} R = 0.5, PF");

  tex->DrawLatex(0.19,0.78,"Barrel");
  tex->DrawLatex(0.47,0.78,"Endcap"); //0.42
  tex->DrawLatex(0.73,0.78,"Forward");

  tex->DrawLatex(0.21,0.73,"BB");
  tex->DrawLatex(0.43,0.73,"EC1");
  tex->DrawLatex(0.57,0.73,"EC2");
  tex->DrawLatex(0.77,0.73,"HF");

  c1->SaveAs("pdf/drawATLASresponse.pdf");
} // drawATLASresponse
