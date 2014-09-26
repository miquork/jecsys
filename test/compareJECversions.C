#include "TCanvas.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TLine.h"

#include "TLatex.h"
#include "tdrstyle_mod12.C"
//#include "settings12.h"

#include <fstream>

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

// symmetric if both below false
// positive true overrides negative
bool _usenegative = false;
bool _usepositive = false;//true;

bool _useptgen = false; // iterate to produce JEC vs pTgen
bool _compare = false; // not used yet?

const double _mu = 19.83; // 20/fb at 8 TeV (htrpu)
const double _lumi = 19800.;
bool _pdf = true;
bool _mc(false);
string _alg("");

// Mapping from mu to NPV and rho
// Derived from inclusive jet data, |eta|<1.3, jt320: prhovstrpu, pnpvvstrpu
// TF1 *f1 = new TF1("f1","[0]+[1]*x+[2]*x*x",7,33);
// for (int i=0; i!=f1->GetNpar(); ++i) cout<<Form("%1.4g, ",f1->GetParameter(i)); cout<<endl;
//
// prhovstrpu->Fit(f1,"R")
double getRho(double mu) {

  double p[3] = {1.009, 0.5515, 0.0003597}; // DATA

  return (p[0]+p[1]*mu+p[2]*mu*mu);

}
// pnpvvstrpu->Fit(f1,"R")
double getNPV(double mu) {

  double p[3] = {1.032, 0.7039, -0.001319}; // DATA

  return (p[0]+p[1]*mu+p[2]*mu*mu);
}

void setEtaPtE(FactorizedJetCorrector *jec, double eta, double pt, double e,
	       int mu) {

  int npv = getNPV(mu);
  double rho = getRho(mu);

  jec->setJetEta(eta);
  jec->setJetPt(pt);
  // L1Offset
  jec->setJetE(e);
  jec->setNPV(npv);
  // L1FastJEt
  bool is5 = (_alg=="AK5PF"||_alg=="AK5PFchs"||_alg=="AK5CALO");
  double jeta = TMath::Pi()*(is5 ? 0.5*0.5 : 0.7*0.7);
  jec->setJetA(jeta);
  jec->setRho(rho);

  return;
}

double getEtaPtE(FactorizedJetCorrector *jec, double eta, double pt, double e,
		 int mu = _mu) {

  setEtaPtE(jec, eta, pt, e, mu);

  // if using pTgen, need to iterate to solve ptreco
  if (_useptgen) {
    double ptgen = pt;
    assert(false); // need to develop this further
  }

  return (jec->getCorrection());
}

double getEtaPtUncert(JetCorrectionUncertainty *unc, double eta, double pt) {
  
  unc->setJetEta(eta);
  unc->setJetPt(pt);
  return (unc->getUncertainty(true));
  
}


void compareJECversions(string algo="AK5PFchs",
			bool l1=true, bool l2l3=true, bool res=true,
			string type="DATA") {

  gROOT->ProcessLine(".L tdrstyle_mod12.C");
  setTDRStyle();
  
  assert(type=="DATA" || type=="MC");
  const bool mc = (type=="MC");
  _mc = mc;
  _alg = algo;
  assert(!mc || !res);
  const char *a = algo.c_str();
  const char *str;

  // Legends
  const char *cm0 = "DATA";
  const char *cm = type.c_str();
  const char *s2 = "Winter14_V5";
  const char *s2s = "V5";
  const char *s1 = "Summer13";
  const char *s1s = "GT";
  string sgen = (_useptgen ? "ptcl" : "raw");
  const char *cgen = sgen.c_str();

  // Summer13
  string sid1 = (_mc ? "START53_V26_MC" : "FT_53_V21_AN6_DATA");
  const char *cid1 = sid1.c_str();
  const char *a1 = a;
  str=Form("CondFormats/JetMETObjects/data/%s_L1FastJet_%s.txt",cid1,a1);
  cout << str << endl << flush;
  JetCorrectorParameters *JetCorPar1L1 = new JetCorrectorParameters(str);
  str=Form("CondFormats/JetMETObjects/data/%s_L2Relative_%s.txt",cid1,a1);
  cout << str << endl << flush;
  JetCorrectorParameters *JetCorPar1L2 = new JetCorrectorParameters(str);
  str=Form("CondFormats/JetMETObjects/data/%s_L3Absolute_%s.txt",cid1,a1);
  cout << str << endl << flush;
  JetCorrectorParameters *JetCorPar1L3 = new JetCorrectorParameters(str);
  str=Form("CondFormats/JetMETObjects/data/%s_L2L3Residual_%s.txt",cid1,a1);
  if (!mc) cout << str << endl << flush;
  JetCorrectorParameters *JetCorPar1 = (mc ? 0 : new JetCorrectorParameters(str));
  str=Form("CondFormats/JetMETObjects/data/%s_Uncertainty_%s.txt",cid1,a1);
  cout << str << endl << flush;
  JetCorrectionUncertainty *jecUnc1 = new JetCorrectionUncertainty(str);

  // JEC in 2010 JINST paper
  /*
    str=Form("CondFormats/JetMETObjects/data/Fall10_L1Offset_%s_SCALED.txt",a); cout << str << endl << flush;
    JetCorrectorParameters *JetCorPar1L1 = new JetCorrectorParameters(str);
    str=Form("CondFormats/JetMETObjects/data/Fall10_L2Relative_%s.txt",a); cout << str << endl << flush;
    JetCorrectorParameters *JetCorPar1L2 = new JetCorrectorParameters(str);
    str=Form("CondFormats/JetMETObjects/data/Fall10_L3Absolute_%s.txt",a); cout << str << endl << flush;
    JetCorrectorParameters *JetCorPar1L3 = new JetCorrectorParameters(str);
    str=Form("CondFormats/JetMETObjects/data/Fall10_L2L3Residual_%s.txt",a); if (mc) cout << str << endl << flush;
    JetCorrectorParameters *JetCorPar1 = (mc ? 0 : new JetCorrectorParameters(str));
    str=Form("CondFormats/JetMETObjects/data/Jec10V1_Uncertainty_%s.txt",a); cout << str << endl << flush;
    JetCorrectionUncertainty *jecUnc1 = new JetCorrectionUncertainty(str);
  */

  // Latest 53X V5 (Winter14_V5)
  string sid2 = (_mc ? "Winter14_V4_MC" : "Winter14_V4_DATA");
  const char *cid2 = sid2.c_str();
  const char *a2 = a;
  str=Form("CondFormats/JetMETObjects/data/%s_L1FastJet_%s.txt",cid2,a2);
  cout << str << endl << flush;
  JetCorrectorParameters *JetCorPar2L1 = new JetCorrectorParameters(str);
  str=Form("CondFormats/JetMETObjects/data/%s_L2Relative_%s.txt",cid2,a2);
  cout << str << endl << flush;
  JetCorrectorParameters *JetCorPar2L2 = new JetCorrectorParameters(str);
  str=Form("CondFormats/JetMETObjects/data/%s_L3Absolute_%s.txt",cid2,a2);
  cout << str << endl << flush;
  JetCorrectorParameters *JetCorPar2L3 = new JetCorrectorParameters(str);
  str=Form("CondFormats/JetMETObjects/data/%s_L2L3Residual_%s.txt",cid2,a2);
  if (!mc) cout << str << endl << flush;
  JetCorrectorParameters *JetCorPar2 = (mc ? 0 : new JetCorrectorParameters(str));
  //str=Form("CondFormats/JetMETObjects/data/%s_Uncertainty_%s.txt",cid2,a2);
  str=Form("../txt/Winter14_V8M_DATA_Uncertainty_%s.txt",a2);
  cout << str << endl << flush;
  JetCorrectionUncertainty *jecUnc2 = new JetCorrectionUncertainty(str);

  vector<JetCorrectorParameters> vParam1;
  if (l1)   vParam1.push_back(*JetCorPar1L1);
  if (l2l3) vParam1.push_back(*JetCorPar1L2);
  if (l2l3) vParam1.push_back(*JetCorPar1L3);
  if (res && !mc && JetCorPar1)  vParam1.push_back(*JetCorPar1);
  vector<JetCorrectorParameters> vParam2;
  if (l1)   vParam2.push_back(*JetCorPar2L1);
  if (l2l3) vParam2.push_back(*JetCorPar2L2);
  if (l2l3) vParam2.push_back(*JetCorPar2L3);
  if (res && !mc && JetCorPar2)  vParam2.push_back(*JetCorPar2);

  FactorizedJetCorrector *JEC1 = new FactorizedJetCorrector(vParam1);
  FactorizedJetCorrector *JEC2 = new FactorizedJetCorrector(vParam2);
  //_JEC1 = JEC1;

  TCanvas *c0 = new TCanvas(Form("c0_%s",a),Form("c0_%s",a),600,600);
  TCanvas *c1a = new TCanvas(Form("c1a_%s",a),Form("c1a_%s",a),600,600);
  TCanvas *c1b = new TCanvas(Form("c1b_%s",a),Form("c1b_%s",a),600,600);
  TCanvas *c1c = new TCanvas(Form("c1c_%s",a),Form("c1c_%s",a),600,600);

  TCanvas *c2 = new TCanvas(Form("c2_%s",a),Form("c2_%s",a),600,600);

  TH1D *h = new TH1D(Form("h_%s",a),Form(";|#eta|;%s L2L3 residual",a),
		     50,0,5);
  if (_usenegative) h->GetXaxis()->SetTitle("-|#eta|");
  if (_usepositive) h->GetXaxis()->SetTitle("+|#eta|");
  const char *cl1 = (l1 ? "L1" : "");
  const char *cl2l3 = (l2l3 ? "L2L3" : "");
  const char *cpl = (res&&(l1||l2l3) ? "+" : "");
  const char *cplus = (res&&(l1||l2l3) ? "Plus" : "");
  const char *cres = (res ? "L2L3res" : "");

  //h->GetYaxis()->SetTitle(Form("%s %s%s%s%s",a,cl1,cl2l3,cpl,cres));
  /*
  h->SetMinimum(l1 ? 0.85 : 1.00);
  h->SetMaximum(l1 ? 1.30 : 1.45);
  if ((res && !l1 && !l2l3) || (l2l3 && !res && !l1)) {
    h->SetMinimum(0.95);
    h->SetMaximum(1.25);
  }
  if (l2l3 && algo=="CALO") h->SetMaximum(2.2);
  */
  h->GetYaxis()->SetTitle(Form("%s%s%s%s",cl1,cl2l3,cpl,cres));
  h->SetMinimum(0.3);
  h->SetMaximum(2.0);//1.6);

  TGraph *g1a = new TGraph(0);
  TGraph *g1b = new TGraph(0);
  TGraph *g1c = new TGraph(0);
  //
  TGraph *g2a = new TGraph(0);
  TGraph *g2b = new TGraph(0);
  TGraph *g2c = new TGraph(0);
  //
  TGraph *g21a = new TGraph(0);
  TGraph *g21b = new TGraph(0);
  TGraph *g21c = new TGraph(0);

  TGraphErrors *g1a_e = new TGraphErrors(0);
  TGraphErrors *g1b_e = new TGraphErrors(0);
  TGraphErrors *g1c_e = new TGraphErrors(0);
  TGraph *g1a_pl = new TGraph(0);
  TGraph *g1a_mn = new TGraph(0);
  TGraph *g1b_pl = new TGraph(0);
  TGraph *g1b_mn = new TGraph(0);
  TGraph *g1c_pl = new TGraph(0);
  TGraph *g1c_mn = new TGraph(0);

  TGraphErrors *g2a_e = new TGraphErrors(0);
  TGraphErrors *g2b_e = new TGraphErrors(0);
  TGraphErrors *g2c_e = new TGraphErrors(0);
  TGraph *g2a_pl = new TGraph(0);
  TGraph *g2a_mn = new TGraph(0);
  TGraph *g2b_pl = new TGraph(0);
  TGraph *g2b_mn = new TGraph(0);
  TGraph *g2c_pl = new TGraph(0);
  TGraph *g2c_mn = new TGraph(0);

  
  // Different projections in one loop
  /*
  for (int icase = 0; icase != 3; ++icase) {

    TGraph *g1 = new TGraph(0);
    TGraph *g2 = new TGraph(0);
    TGraph *g21 = new TGraph(0);

    TGraphErrors *g1_e = new TGraphErrors(0);
    TGraph *g1_pl = new TGraph(0);
    TGraph *g1_mn = new TGraph(0);
    TGraphErrors *g2_e = new TGraphErrors(0);
    TGraph *g2_pl = new TGraph(0);
    TGraph *g2_mn = new TGraph(0);

    for (int i = 1; i != h->GetNbinsX()+1; ++i) {
      double eta = h->GetBinCenter(i);
    } // for i

    TGraphErrors *g1_e = new TGraphErrors(0);
    TGraph *g1_pl = new TGraph(0);
    TGraph *g1_mn = new TGraph(0);
    TGraphErrors *g2_e = new TGraphErrors(0);
    TGraph *g2_pl = new TGraph(0);
    TGraph *g2_mn = new TGraph(0);
  } // icase
  */

  const int npt = 6;
  //double ptbins[npt] = {30, 40, 50, 80, 140,500};
  double ptbins[npt] = {30, 60, 120, 240, 480, 960};
  TGraphErrors *g21s[npt];
  for (int i = 0; i != npt; ++i) {
    g21s[i] = new TGraphErrors(0);
  }
  
  for (int i = 1; i != h->GetNbinsX()+1; ++i) {

    double eta = h->GetBinCenter(i);

    // ***** Pt = 30, 50, 80, 120, 200, 500 *****
    {
      for (int j = 0; j != npt; ++j) {

	TGraphErrors *g21 = g21s[j];
	double pt = ptbins[j];
	double energy = pt*cosh(eta);

	if (energy < 4000.) {
	  // Asymmetric corrections now
	  double y1 = 0.5*(getEtaPtE(JEC1, +eta, pt, energy)
			    + getEtaPtE(JEC1, -eta, pt, energy));
	  double y2 = 0.5*(getEtaPtE(JEC2, +eta, pt, energy)
			    + getEtaPtE(JEC2, -eta, pt, energy));

	  g21->SetPoint(g21->GetN(), eta, y2/y1);
	} // energy < 4000
      } // for j
    } // pt bins
 
    // ***** Pt = 30 
    {
      double pt = 30.;
      double energy = pt*cosh(eta);
      
      if (energy < 3500.) {
	// Asymmetric corrections now
	double y1 = 0.5*(getEtaPtE(JEC1, +eta, pt, energy)
			  + getEtaPtE(JEC1, -eta, pt, energy));
	double y2 = 0.5*(getEtaPtE(JEC2, +eta, pt, energy)
			  + getEtaPtE(JEC2, -eta, pt, energy));
	// negative side
	if (_usenegative) {
	  y1 = getEtaPtE(JEC1, -eta, pt, energy);
	  y2 = getEtaPtE(JEC2, -eta, pt, energy);
	}
	// positive side
	if (_usepositive) {
	  y1 = getEtaPtE(JEC1, +eta, pt, energy);
	  y2 = getEtaPtE(JEC2, +eta, pt, energy);
	}

	double e1 = getEtaPtUncert(jecUnc1, eta, pt);
	double e2 = getEtaPtUncert(jecUnc2, eta, pt);
	
	g1a->SetPoint(g1a->GetN(), eta, y1);
	g2a->SetPoint(g2a->GetN(), eta, y2);
	g21a->SetPoint(g21a->GetN(),eta, y2/y1);
	//
	g1a_pl->SetPoint(g1a_pl->GetN(), eta, y1*(1+e1));
	g1a_mn->SetPoint(g1a_mn->GetN(), eta, y1*(1-e1));
	g1a_e->SetPoint(i-1, eta, y1);
	g1a_e->SetPointError(i-1, 0., y1*e1);
	//
	g2a_pl->SetPoint(g2a_pl->GetN(), eta, y2*(1+e2));
	g2a_mn->SetPoint(g2a_mn->GetN(), eta, y2*(1-e2));
	g2a_e->SetPoint(i-1, eta, y2);
	g2a_e->SetPointError(i-1, 0., y2*e2);
      }
    }

    // ***** Pt = 100 
    {
      double pt = 100.;
      double energy = pt*cosh(eta);
      
      if (energy < 3500.) {
	// Asymmetric corrections now
	double y1 = 0.5*(getEtaPtE(JEC1, +eta, pt, energy)
			  + getEtaPtE(JEC1, -eta, pt, energy));
	double y2 = 0.5*(getEtaPtE(JEC2, +eta, pt, energy)
			  + getEtaPtE(JEC2, -eta, pt, energy));
	// negative side
	if (_usenegative) {
	  y1 = getEtaPtE(JEC1, -eta, pt, energy);
	  y2 = getEtaPtE(JEC2, -eta, pt, energy);
	}
	// positive side
	if (_usepositive) {
	  y1 = getEtaPtE(JEC1, +eta, pt, energy);
	  y2 = getEtaPtE(JEC2, +eta, pt, energy);
	}

	double e1 = getEtaPtUncert(jecUnc1, eta, pt);
	double e2 = getEtaPtUncert(jecUnc2, eta, pt);
	
	g1b->SetPoint(g1b->GetN(), eta, y1);
	g2b->SetPoint(g2b->GetN(), eta, y2);
	g21b->SetPoint(g21b->GetN(),eta, y2/y1);
	//
	g1b_pl->SetPoint(g1b_pl->GetN(), eta, y1*(1+e1));
	g1b_mn->SetPoint(g1b_mn->GetN(), eta, y1*(1-e1));
	g1b_e->SetPoint(i-1, eta, y1);
	g1b_e->SetPointError(i-1, 0., y1*e1);
	//
	g2b_pl->SetPoint(g2b_pl->GetN(), eta, y2*(1+e2));
	g2b_mn->SetPoint(g2b_mn->GetN(), eta, y2*(1-e2));
	g2b_e->SetPoint(i-1, eta, y2);
	g2b_e->SetPointError(i-1, 0., y2*e2);
      }
    }

    // ***** E = 1000 
    {
      double energy = 1000.;
      double pt = energy/cosh(eta);
     
      if (pt > 10.) {
	// Asymmetric corrections now
	double y1 = 0.5*(getEtaPtE(JEC1, +eta, pt, energy)
			  + getEtaPtE(JEC1, -eta, pt, energy));
	double y2 = 0.5*(getEtaPtE(JEC2, +eta, pt, energy)
			  + getEtaPtE(JEC2, -eta, pt, energy));
	// negative side
	if (_usenegative) {
	  y1 = getEtaPtE(JEC1, -eta, pt, energy);
	  y2 = getEtaPtE(JEC2, -eta, pt, energy);
	}
	// positive side
	if (_usepositive) {
	  y1 = getEtaPtE(JEC1, +eta, pt, energy);
	  y2 = getEtaPtE(JEC2, +eta, pt, energy);
	}
	double e1 = getEtaPtUncert(jecUnc1, eta, pt);
	double e2 = getEtaPtUncert(jecUnc2, eta, pt);
	
	g1c->SetPoint(g1c->GetN(), eta, y1);
	g2c->SetPoint(g2c->GetN(), eta, y2);
	g21c->SetPoint(g21c->GetN(),eta, y2/y1);
	//
	g1c_pl->SetPoint(g1c_pl->GetN(), eta, y1*(1+e1));
	g1c_mn->SetPoint(g1c_mn->GetN(), eta, y1*(1-e1));
	g1c_e->SetPoint(i-1, eta, y1);
	g1c_e->SetPointError(i-1, 0., y1*e1);
	//
	g2c_pl->SetPoint(g2c_pl->GetN(), eta, y2*(1+e2));
	g2c_mn->SetPoint(g2c_mn->GetN(), eta, y2*(1-e2));
	g2c_e->SetPoint(i-1, eta, y2);
	g2c_e->SetPointError(i-1, 0., y2*e2);
      }
    }
  }

  // Generic legend
  TLegend *leg = new TLegend(0.20,0.75,0.40,0.85,"","brNDC");
  leg->SetTextSize(0.045);
  leg->SetBorderSize(0);
  leg->SetFillStyle(kNone);
  leg->AddEntry(g2a,s2,"LPF");
  leg->AddEntry(g1a,s1,"LPF");
  leg->Draw();

  g2a->SetFillStyle(3003);
  g2a->SetFillColor(kRed);

  g1a->SetFillStyle(3003);
  g1a->SetFillColor(kBlue);
  
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);

  // ***** Pt = 30
  {
    c1a->cd();
    h->DrawClone("AXIS");

    g1a_e->SetFillStyle(3003);
    g1a_e->SetFillColor(kBlue);
    g1a_e->Draw("SAME E3");
    g1a_pl->SetLineColor(kBlue);
    g1a_pl->SetLineStyle(kDotted);
    g1a_pl->Draw("SAMEL");
    g1a_mn->SetLineColor(kBlue);
    g1a_mn->SetLineStyle(kDotted);
    g1a_mn->Draw("SAMEL");

    g2a_e->SetFillStyle(3003);
    g2a_e->SetFillColor(kRed);
    g2a_e->Draw("SAME E3");
    g2a_pl->SetLineColor(kRed);
    g2a_pl->SetLineStyle(kDotted);
    g2a_pl->Draw("SAMEL");
    g2a_mn->SetLineColor(kRed);
    g2a_mn->SetLineStyle(kDotted);
    g2a_mn->Draw("SAMEL");
        
    g1a->SetMarkerStyle(kFullSquare);
    g1a->SetMarkerColor(kBlue);
    g1a->SetLineColor(kBlue);
    g1a->Draw("SAMEPL");

    g2a->SetMarkerStyle(kFullCircle);
    g2a->SetMarkerColor(kRed);
    g2a->SetLineColor(kRed);
    g2a->Draw("SAMEPL");

    tex->DrawLatex(0.20,0.88,Form("p_{T,%s} = 30 GeV%s, %s",cgen,
				  l1||l2l3 ? Form(", #LT#mu#GT = %1.1f",_mu) : "",
				  cm));
    tex->DrawLatex(0.65,0.80,a);
    leg->Draw();
    if (!mc) cmsPrel(_lumi);
    if (mc)  cmsPrel(0);
    gPad->RedrawAxis();
  }

  // ***** Pt = 100
  {
    c1b->cd();
    h->DrawClone("AXIS");

    g1b_e->SetFillStyle(3003);
    g1b_e->SetFillColor(kBlue);
    g1b_e->Draw("SAME E3");
    g1b_pl->SetLineColor(kBlue);
    g1b_pl->SetLineStyle(kDotted);
    g1b_pl->Draw("SAMEL");
    g1b_mn->SetLineColor(kBlue);
    g1b_mn->SetLineStyle(kDotted);
    g1b_mn->Draw("SAMEL");    

    g2b_e->SetFillStyle(3003);
    g2b_e->SetFillColor(kRed);
    g2b_e->Draw("SAME E3");
    g2b_pl->SetLineColor(kRed);
    g2b_pl->SetLineStyle(kDotted);
    g2b_pl->Draw("SAMEL");
    g2b_mn->SetLineColor(kRed);
    g2b_mn->SetLineStyle(kDotted);
    g2b_mn->Draw("SAMEL");    

    g1b->SetMarkerStyle(kFullSquare);
    g1b->SetMarkerColor(kBlue);
    g1b->SetLineColor(kBlue);
    g1b->Draw("SAMEPL");
    
    g2b->SetMarkerStyle(kFullCircle);
    g2b->SetMarkerColor(kRed);
    g2b->SetLineColor(kRed);
    g2b->Draw("SAMEPL");
    
    tex->DrawLatex(0.20,0.88,Form("p_{T,%s} = 100 GeV%s, %s",cgen,
				  l1||l2l3 ? Form(", #LT#mu#GT = %1.1f",_mu) : "",
				  cm));
    tex->DrawLatex(0.65,0.80,a);
    leg->Draw();
    if (!mc) cmsPrel(_lumi);
    if (mc)  cmsPrel(0);
    gPad->RedrawAxis();
  }

  // ***** E = 1000
  {
    c1c->cd();
    h->DrawClone("AXIS");

    g1c_e->SetFillStyle(3003);
    g1c_e->SetFillColor(kBlue);
    g1c_e->Draw("SAME E3");
    g1c_pl->SetLineColor(kBlue);
    g1c_pl->SetLineStyle(kDotted);
    g1c_pl->Draw("SAMEL");
    g1c_mn->SetLineColor(kBlue);
    g1c_mn->SetLineStyle(kDotted);
    g1c_mn->Draw("SAMEL");

    g2c_e->SetFillStyle(3003);
    g2c_e->SetFillColor(kRed);
    g2c_e->Draw("SAME E3");
    g2c_pl->SetLineColor(kRed);
    g2c_pl->SetLineStyle(kDotted);
    g2c_pl->Draw("SAMEL");
    g2c_mn->SetLineColor(kRed);
    g2c_mn->SetLineStyle(kDotted);
    g2c_mn->Draw("SAMEL");
    
    g1c->SetMarkerStyle(kFullSquare);
    g1c->SetMarkerColor(kBlue);
    g1c->SetLineColor(kBlue);
    g1c->Draw("SAMEPL");
    
    g2c->SetMarkerStyle(kFullCircle);
    g2c->SetMarkerColor(kRed);
    g2c->SetLineColor(kRed);
    g2c->Draw("SAMEPL");
    
    tex->DrawLatex(0.20,0.88,Form("E_{%s} = 1000 GeV%s, %s",cgen,
				  l1||l2l3 ? Form(", #LT#mu#GT = %1.1f",_mu) : "",
				  cm));
    tex->DrawLatex(0.65,0.80,a);
    leg->Draw();
    if (!mc) cmsPrel(_lumi);
    if (mc)  cmsPrel(0);
    gPad->RedrawAxis();
  }

  string ctype = string(cl1)+string(cl2l3)+string(cplus)+string(cres);
  const char *cs = ctype.c_str();
  if(_pdf) c1a->SaveAs(Form("pdf/compareJECversions_%s_%s_%s_Pt030.pdf",a,cm,cs));
  if(_pdf) c1b->SaveAs(Form("pdf/compareJECversions_%s_%s_%s_Pt100.pdf",a,cm,cs));
  if(_pdf) c1c->SaveAs(Form("pdf/compareJECversions_%s_%s_%s_Q1000.pdf",a,cm,cs));


  // ***** Multiple pT bins for ratio only
  {
    int colors[] = {kBlack, kBlue, kCyan+2, kGreen+2, kOrange+2, kRed};
    int styles[] = {kSolid, kDashed, kDotted, kDashDotted, kDashed, kSolid};
    c0->cd();
    h->SetMinimum(0.85);//0.90);
    h->SetMaximum(1.25);//1.20);
    TH1D *h0 = (TH1D*)h->DrawClone("AXIS");
    h0->GetYaxis()->SetTitle(Form("%s%s%s%s (%s / %s)",cl1,cl2l3,cpl,cres,
				  s2s, s1s));

    TLine *l = new TLine();
    l->SetLineStyle(kDashed);
    l->DrawLine(0,1,5,1);
    l->DrawLine(0,1.05,5,1.05);
    l->SetLineStyle(kDotted);
    l->DrawLine(0,1.02,5,1.02);
    l->DrawLine(0,0.98,5,0.98);

    TLegend *leg = new TLegend(0.20,0.62,0.60,0.92,"","brNDC");
    leg->SetTextSize(0.045);
    leg->SetBorderSize(0);
    leg->SetFillStyle(kNone);
    leg->Draw();

    for (int i = 0; i != npt; ++i) {
      TGraphErrors *gr = g21s[i];
      gr->SetLineColor(colors[i]);//i+1);
      gr->SetLineStyle(styles[i]);//i+1);
      gr->SetLineWidth(3);
      gr->Draw("SAME L");

      //leg->AddEntry(gr,Form("%s / %s (p_{T,%s}=%1.0f GeV)",
      //s2s, s1s, cgen,ptbins[i]),"L");
      leg->AddEntry(gr,Form("p_{T,%s} = %1.0f GeV",
			    cgen,ptbins[i]),"L");
    } // for i
    tex->DrawLatex(0.70,0.85,a);
    cmsPrel(_lumi);
  }
  c0->SaveAs(Form("pdf/compareJECversions_%s_%s_%s_%sover%s.pdf",a,cm,cs,s2s,s1s));

  
  {// Ratio plots
    c2->cd();

    h->SetMinimum(l1||l2l3 ? 0.80 : 0.95);
    h->SetMaximum(l1||l2l3 ? 1.20 : 1.12);
    if (res && (algo=="CALO" || algo=="JPT")) h->SetMinimum(0.80);
    h->DrawClone("AXIS");

    TLine *l = new TLine();
    l->SetLineStyle(kDashed);
    l->DrawLine(0,1,5,1);
    l->SetLineColor(kGreen+2);
    l->SetLineStyle(kDotted);
    l->DrawLine(0,1.02,5,1.02);
    l->DrawLine(0,0.98,5,0.98);
    l->SetLineColor(kBlue);
    l->SetLineStyle(kDotted);
    l->DrawLine(0,1.05,5,1.05);
    l->DrawLine(0,0.95,5,0.95);

    g21a->SetMarkerStyle(kOpenCircle);
    g21a->SetMarkerColor(kBlue);
    g21a->SetLineColor(kBlue);
    g21a->SetLineStyle(kDashed);
    g21a->SetLineWidth(3);
    g21a->Draw("SAMEL");

    g21c->SetMarkerStyle(kOpenSquare);
    g21c->SetMarkerColor(kRed);
    g21c->SetLineColor(kRed);
    g21c->SetLineStyle(kDotted);
    g21c->SetLineWidth(3);
    g21c->Draw("SAMEL");

    g21b->SetMarkerStyle(kFullCircle);
    g21b->SetMarkerColor(kGreen+2);
    g21b->SetLineColor(kGreen+2);
    g21b->SetLineStyle(kSolid);
    g21b->SetLineWidth(3);
    g21b->Draw("SAMEL");
    
    TLegend *leg = new TLegend(0.20,0.77,0.60,0.92,"","brNDC");
    leg->SetTextSize(0.045);
    leg->SetBorderSize(0);
    leg->SetFillStyle(kNone);
    leg->AddEntry(g21a,Form("%s / %s (p_{T,%s}=30 GeV)",s2s,s1s,cgen),"LP");
    leg->AddEntry(g21b,Form("%s / %s (p_{T,%s}=100 GeV)",s2s,s1s,cgen),"LP");
    leg->AddEntry(g21c,Form("%s / %s (E_{%s}=1000 GeV)",s2s,s1s,cgen),"LP");
    leg->Draw();

    if (!mc) cmsPrel(_lumi);
    if (mc)  cmsPrel(0);
    gPad->RedrawAxis();
    
    if(_pdf) c2->SaveAs(Form("pdf/compareJECversions_%s_%s_%s_Ratios.pdf",a,cm,cs));
  } // Ratio plots
} // compareJECversions


