#include "TCanvas.h"
#include "TH1D.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLine.h"

#include "../CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "../CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

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
const double emax = 6500;//3800; //6500
const double ptmin = 15.;
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
  double ptmeas = _jecpt->GetX(ptgen, 1, emax);
  double resp = ptmeas / _jecpt->Eval(ptmeas); // 1/jec

  return resp;
}

void drawCMSresponse() {

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

    const char *sd = "../CondFormats/JetMETObjects/data";
    //const char *st = "Winter14_V5_MC";
    //const char *st = "Winter14_V8_MC"; // 2012
    //const char *st = "Winter14_V8_DATA";
    //const char *st = "Summer15_25nsV3_DATA";
    //const char *st = "Summer15_25nsV6_DATA";
    //const char *st = "Summer15_25nsV7_DATA";
    //const char *st = "Fall15_25nsV1_DATA";
    //const char *st = "Spring16_25nsV3_MC";
    //const char *st = "Summer16_23Sep2016GV3_DATA"; // 2017
    //const char *st = "Summer16_03Feb2017G_V3_DATA"; // 2017 03FebV3
    //const char *st = "Summer16_03Feb2017BCD_V7_DATA"; // 2017 03FebV7
    const char *st = "Summer16_07Aug2017GH_V10_DATA"; // 2017 03FebV
    //const char *st = "Summer16_07Aug2017BCD_V10_DATA"; // 2017 03FebV7
    const char *s;

    //s = Form("%s/%s_L1FastJet_AK5PFchs.txt",sd,st); cout << s << endl;
    //JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
    s = Form("%s/%s_L1FastJet_AK4PFchs.txt",sd,st); cout << s << endl;
    JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
    //s = Form("%s/%s_L2Relative_AK5PFchs.txt",sd,st); cout << s << endl;
    s = Form("%s/%s_L2Relative_AK4PFchs.txt",sd,st); cout << s << endl;
    JetCorrectorParameters *l2 = new JetCorrectorParameters(s);
    //s = Form("%s/%s_L3Absolute_AK5PFchs.txt",sd,st); cout << s << endl;
    s = Form("%s/%s_L3Absolute_AK4PFchs.txt",sd,st); cout << s << endl;
    JetCorrectorParameters *l3 = new JetCorrectorParameters(s);
    s = Form("%s/%s_L2L3Residual_AK4PFchs.txt",sd,st); cout << s << endl;
    JetCorrectorParameters *l2l3 = new JetCorrectorParameters(s);

    vector<JetCorrectorParameters> v;
    v.push_back(*l1);
    //v.push_back(*l2);
    //v.push_back(*l3);
    //v.push_back(*l2l3);
    _jec = new FactorizedJetCorrector(v);
  }
  if (!_jecpt) {
    _jecpt = new TF1("jecpt",fJECPt,ptmin,emax,3);
  }

  const double etabins[] = 
    //{-5.191, -3.839, -3.489, -3.139, -2.964, -2.853, -2.65, -2.5, -2.322,
    //-2.172, -1.93, -1.653, -1.479, -1.305, -1.044, -0.783, -0.522, -0.261,
    //0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322,
    //2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};
    {0, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879,
     0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.74, 1.83,
     1.93, 2.043, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489,
     3.664, 3.839, 4.013, 4.191, 4.363, 4.538, 4.716, 4.889, 5.191};
  const int neta = sizeof(etabins)/sizeof(etabins[0])-1; 

  //double ergs[] = {30, 60, 110, 400, 2000};
  //const int ne = sizeof(ergs)/sizeof(ergs[0]);
  //double pts[] = {10, 30, 100, 400, 2000};
  double pts[] = {15, 30, 60, 110, 400, 2000};
  const int npt = sizeof(pts)/sizeof(pts[0]);
  //const int neta = 48;//52;
  const int jeta = TMath::Pi()*0.4*0.4;
  const int mu = 25;

  TGraph *gs[npt];
  //for (int ie = 0; ie != ne; ++ie) {
  for (int ipt = 0; ipt != npt; ++ipt) {

    //double energy = ergs[ie];
    double pt = pts[ipt];

    TGraph *g = new TGraph(0); gs[ipt] = g;
    for (int ieta = 0; ieta != neta; ++ieta) {
      
      //double eta = (ieta+0.5)*0.1;
      double eta = 0.5*(etabins[ieta]+etabins[ieta+1]);
      //double pt = energy / cosh(eta);
      double energy = pt * cosh(eta);
      if (pt >= ptmin && energy < emax) { // 13 TeV

	//double jes = getResp(pt, eta, jeta, mu);
	double jes_plus = getResp(pt, eta, jeta, mu);
	double jes_minus = getResp(pt, -eta, jeta, mu);
	double jes = 0.5 * (jes_plus + jes_minus);
	//double jes = jes_plus;
	//double jes = jes_minus;
	int n = g->GetN();
	g->SetPoint(n, eta, jes);
      }
    } // for ie
  } // for ieta


  // Draw results
  TH1D *h = new TH1D("h",";Jet |#eta|;Simulated jet response",40,0,4.8);
  //TH1D *h = new TH1D("h",";Jet |#eta|;Data jet response",40,0,4.8);
  //TH1D *h = new TH1D("h",";Jet |#eta|;Data response+offset",40,0,4.8);
  h->SetMaximum(1.25);
  h->SetMinimum(0.5);
  //extraText = "Simulation";
  extraText = "Simulation Preliminary";
  //extraText = "Preliminary";
  lumi_8TeV = "";
  lumi_13TeV = "";
  //lumi_13TeV = "2.1 fb^{-1}";
  //TCanvas *c1 = tdrCanvas("c1",h,2,0,kSquare);
  TCanvas *c1 = tdrCanvas("c1",h,4,0,kSquare);

  TLegend *leg1 = tdrLeg(0.25,0.25,0.55,0.30);
  TLegend *leg2 = tdrLeg(0.25,0.20,0.55,0.25);
  TLegend *leg3 = tdrLeg(0.25,0.15,0.55,0.20);
  TLegend *leg4 = tdrLeg(0.55,0.25,0.85,0.30);
  TLegend *leg5 = tdrLeg(0.55,0.20,0.85,0.25);
  TLegend *leg6 = tdrLeg(0.55,0.15,0.85,0.20);
  TLegend *legs[npt] = {leg1, leg2, leg3, leg4, leg5, leg6};

  int colors[] = {kGray+2, kGreen+2, kBlack, kOrange+1, kBlue, kRed+1};
  int markers[] = {kOpenDiamond, kFullCircle, kOpenCircle, kFullSquare,
		   kOpenSquare, kFullTriangleUp};

  for (int ipt = 0; ipt != npt; ++ipt) {
    
    TGraph *g = gs[ipt];
    g->SetMarkerColor(colors[ipt]);
    g->SetMarkerStyle(markers[ipt]);
    g->Draw("SAMEP");

    //TLegend *leg = (ie<3 ? leg1 : leg2);
    TLegend *leg = legs[ipt];
    leg->SetTextColor(colors[ipt]);
    //leg->AddEntry(g, Form("E = %1.0f GeV",ergs[ie]), "P");
    leg->AddEntry(g, Form("p_{T} = %1.0f GeV",pts[ipt]), "P");
  }


  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);
  
  TLine *l = new TLine();
  // Run I
  /*
  l->DrawLine(1.3,0.7,1.3,1.1);
  l->DrawLine(2.5,0.7,2.5,1.1);
  l->DrawLine(3.0,0.7,3.0,1.1);
  l->DrawLine(4.5,0.7,4.5,1.1);
  l->SetLineStyle(kDashed);
  l->DrawLine(3.2,0.7,3.2,1.1);
  */
  // Run II
  l->DrawLine(1.305,0.66,1.305,1.06); // BB
  l->DrawLine(2.5,0.66,2.5,1.06); // EC1
  l->DrawLine(2.853,0.66,2.853,1.06); // EC2
  l->DrawLine(3.139,0.66,3.139,1.06); // HF
  l->DrawLine(4.538,0.66,4.538,1.06); // HF
  l->SetLineStyle(kDashed);
  //l->DrawLine(3.139,0.66,3.139,1.06); // HF
  l->DrawLine(2.964,0.66,2.964,1.00); // HF

  //tex->DrawLatex(0.23,0.86,"2012 JES: Anti-k_{t} R = 0.5, PF+CHS");
  //tex->DrawLatex(0.30,0.86,"53X JES: Anti-k_{t} R = 0.5, PF+CHS");
  //tex->DrawLatex(0.30,0.86,"74X JES: Anti-k_{t} R = 0.4, PF+CHS");
  //tex->DrawLatex(0.30,0.86,"76X JES: Anti-k_{t} R = 0.4, PF+CHS");
  tex->DrawLatex(0.19,0.86,"2016 JES: Anti-k_{T} R = 0.4, PF + CHS");
  //tex->DrawLatex(0.23,0.86,"2017 JES: Anti-k_{t} R = 0.4, PF+CHS");
  //tex->DrawLatex(0.23,0.86,"2017 03FebV3: Anti-k_{t} R = 0.4, PF+CHS");
 
  tex->DrawLatex(0.19,0.78,"Barrel");
  tex->DrawLatex(0.47,0.78,"Endcap"); //0.42
  tex->DrawLatex(0.73,0.78,"Forward");

  tex->DrawLatex(0.21,0.73,"BB");
  tex->DrawLatex(0.43,0.73,"EC1");
  //tex->DrawLatex(0.57,0.73,"EC2");
  tex->DrawLatex(0.56,0.73,"EC2");
  tex->DrawLatex(0.635,0.68,"X");
  tex->DrawLatex(0.625,0.63,"1");
  tex->DrawLatex(0.65,0.63,"2");
  tex->DrawLatex(0.77,0.73,"HF");
  
  c1->RedrawAxis();
  c1->SaveAs("../pdf/drawCMSresponse.pdf");
} // drawCMSresponse
