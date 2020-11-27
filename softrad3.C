// File: softrad3.C
// Created by Mikko Voutilainen, on November 11, 2020
// Purpose: Use JEC combination file to derive jet1+jetn+jetu corrections
//          including uncertainty eigenvectors for global refit
//          These replace the FSR+ISR corrections based on alpha extrapolation
//          in the previous generation softrad.C
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
#include "TLine.h"

#include "tdrstyle_mod15.C"

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>

// Find entry in graph correspondign to this histogram bin
double getY(double pt, double ptmin, double ptmax, TGraphErrors *g,
	    double &ey) {
  
  int j(-1);
  double dx(0);
  for (int i = 0; i != g->GetN(); ++i) {
    double x = g->GetX()[i];
    if ((j<0 || fabs(x-pt)<dx) && x>ptmin && x<ptmax) {
      j = i;
      x = g->GetX()[i];
    }
  } // for i

  ey = (j<0 ? 0 : g->GetEY()[j]);
  return (j<0 ? 0 : g->GetY()[j]);
}

void getEntry(int i, double &pt, double &p0, double &p1, double &pn, double &pu,
	      TH1D *h, TGraphErrors *g0, TGraphErrors *g1,
	      TGraphErrors *gn, TGraphErrors *gu,
	      double &ep0, double &ep1, double &epn, double &epu) {

  pt = h->GetBinCenter(i);
  double ptmin = h->GetBinLowEdge(i);
  double ptmax = h->GetBinLowEdge(i+1);

  p0 = getY(pt, ptmin, ptmax, g0, ep0);
  p1 = getY(pt, ptmin, ptmax, g1, ep1);
  pn = getY(pt, ptmin, ptmax, gn, epn);
  pu = getY(pt, ptmin, ptmax, gu, epu);
} // getEntry

// root -l -b -q softrad3.C+g\(0\,1\.3\,false\,\"2018ABCD\"\)

// 3-point soft radiation corrections for L3Res
void softrad3(double etamin=0.0, double etamax=1.3, bool dodijet=false,
	      string epoch="") {
  
  //if (dodijet && fabs(etamin-0.0)<0.1 && fabs(etamax-1.3)<0.1)
  //domultijet = true;
  
  cout << "Calling softrad3("<<etamin<<","<<etamax<<","<<dodijet<<");\n"<<flush;
  const char *cep = epoch.c_str();
  cout << "For epoch " << epoch << endl;
  
  bool isUL18 = (epoch=="2018ABCD" || epoch=="2018A" || 
		 epoch=="2018B" || epoch=="2018C" || epoch=="2018D");
  bool isUL17 = (epoch=="2017BCDEF");
  
  setTDRStyle();
  
  TDirectory *curdir = gDirectory;

  // Open jecdata.root produced by reprocess.C
  TFile *finout = new TFile(Form("rootfiles/jecdata%s.root",epoch.c_str()),
			    "UPDATE");
  assert(finout && !finout->IsZombie());

  curdir->cd();

  vector<string> dirs;
  dirs.push_back("data");
  dirs.push_back("mc");
  dirs.push_back("ratio");

  vector<string> dirs2;
  dirs2.push_back("ratio");

  vector<string> methods;
  methods.push_back("counts");
  methods.push_back("mpf1");
  methods.push_back("mpfn");
  methods.push_back("mpfu");
  methods.push_back("mpfchs1");
  methods.push_back("ptchs");

  vector<string> methods2;
  methods2.push_back("mpfchs1");
  methods2.push_back("ptchs");

  vector<string> samples;
  //samples.push_back("gamjet");
  //samples.push_back("zeejet");
  //samples.push_back("zmmjet");
  if (isUL18) samples.push_back("multijet");
  if (isUL18) samples.push_back("zlljet");
  samples.push_back("zjet");

  //vector<int> alphas;
  //alphas.push_back(100);
  const int aref = 100; // Z/gamma+jet
  const int mptref = 30; // multijet

  map<string, map<string, map<string, TGraphErrors*> > > gs;
  map<string, map<string, TH1D*> > hs;
 
  // Load all data sets first
  for (unsigned int idir = 0; idir != dirs.size(); ++idir) {
    string d = dirs[idir];
    const char *cd = d.c_str();

    for (unsigned int imet = 0; imet != methods.size(); ++imet) {
      string m = methods[imet];
      const char *cm = m.c_str();

      for (unsigned int isam = 0; isam != samples.size(); ++isam) {
	string s = samples[isam];
	const char *cs = s.c_str();
	
	int xref = (s=="multijet" ? mptref : aref);
	string sg = Form("%s/eta%02d-%02d/%s_%s_a%d",
			 cd,int(10*etamin),int(10*etamax),cm,cs,xref);
	const char *cg = sg.c_str();

	if (m=="counts") {
	  if (s=="zlljet") { // PATCH
	    sg = Form("%s/eta%02d-%02d/%s_%s_a%d",
		      cd,int(10*etamin),int(10*etamax),cm,"zmmjet",xref);
	    cg = sg.c_str();
	  }
	  TH1D *h = (TH1D*)finout->Get(cg); //assert(h);
	  if (!h) cout << cg << flush << endl; assert(h);
	  hs[d][s] = h;
	}
	else {
	  TGraphErrors *g = (TGraphErrors*)finout->Get(cg); //assert(g);
	  if (!g) cout << cg << flush << endl; assert(g);
	  gs[d][m][s] = g;
	}
      } // isam
    } // imet
  } // idir

  // Helper function for solving non-linear equations
  TF1 *fm = new TF1("fm","x - ([0] + (x/[3]-1)*[1] + (x/[4]-1)*[2])",0,13000);
  TF1 *fp = new TF1("fp","x - ([0] + (x/[3]  )*[1] + (x/[4]  )*[2])",0,13000);

  // Then start working on corrections
  for (unsigned int idir = 0; idir != dirs2.size(); ++idir) {
    string d = dirs2[idir];
    const char *cd = d.c_str();

    for (unsigned int imet = 0; imet != methods2.size(); ++imet) {
      string m = methods2[imet];
      const char *cm = m.c_str();

      for (unsigned int isam = 0; isam != samples.size(); ++isam) {
	string s = samples[isam];
	const char *cs = s.c_str();

	// mpf1, mpfn and mpfu must be loaded before correcting mpfchs1
	// statistics needed for checking pT binning
	if ((m=="mpfchs1" || m=="ptchs") && d=="ratio") {
	  assert(hs["data"][s]!=0);// || s=="multijet");
	  assert(gs["mc"]["mpf1"][s]!=0);
	  assert(gs["mc"]["mpfn"][s]!=0);
	  assert(gs["mc"]["mpfu"][s]!=0);
	  assert(gs["mc"]["mpfchs1"][s]!=0);
	  assert(gs["data"]["mpf1"][s]!=0);
	  assert(gs["data"]["mpfn"][s]!=0);
	  assert(gs["data"]["mpfu"][s]!=0);
	  assert(gs["data"]["mpfchs1"][s]!=0);
	  assert(gs["ratio"]["mpfchs1"][s]!=0);
	  assert(gs["ratio"]["ptchs"][s]!=0);

	  TH1D *h = hs["data"][s];
	  TGraphErrors *g0 = gs["data"]["mpfchs1"][s];
	  TGraphErrors *g1 = gs["data"]["mpf1"][s];
	  TGraphErrors *gn = gs["data"]["mpfn"][s];
	  TGraphErrors *gu = gs["data"]["mpfu"][s];
	  TGraphErrors *h0 = gs["mc"]["mpfchs1"][s];
	  TGraphErrors *h1 = gs["mc"]["mpf1"][s];
	  TGraphErrors *hn = gs["mc"]["mpfn"][s];
	  TGraphErrors *hu = gs["mc"]["mpfu"][s];
	  //
	  TGraphErrors *gm = gs["ratio"]["mpfchs1"][s];
	  TGraphErrors *gp = gs["ratio"]["ptchs"][s];


	  TH1D *hfsr   = (TH1D*)h->Clone("fsr");   hfsr->Reset();
	  TH1D *hfsr0s = (TH1D*)h->Clone("fsr0s"); hfsr0s->Reset();
	  TH1D *hfsr1s = (TH1D*)h->Clone("fsr1s"); hfsr1s->Reset();
	  TH1D *hfsrn1 = (TH1D*)h->Clone("fsrn1"); hfsrn1->Reset();
	  TH1D *hfsrn2 = (TH1D*)h->Clone("fsrn2"); hfsrn2->Reset();
	  TH1D *hfsrns = (TH1D*)h->Clone("fsrns"); hfsrns->Reset();
	  TH1D *hfsru1 = (TH1D*)h->Clone("fsru1"); hfsru1->Reset();
	  TH1D *hfsru2 = (TH1D*)h->Clone("fsru2"); hfsru2->Reset();
	  TH1D *hfsrus = (TH1D*)h->Clone("fsrus"); hfsrus->Reset();
	  
	  double pt(0), p0(0), p1(0), pn(0), pu(0);
	  double pm(0), q0(0), q1(0), qn(0), qu(0);
	  double eq0(0), eq1(0), eqn(0), equ(0);
	  double pd(0), r0(0), r1(0), rn(0), ru(0);
	  double er0(0), er1(0), ern(0), eru(0);
	  double mpf(0), ptbal(0);
	  for (int i = 1; i != h->GetNbinsX()+1; ++i) {

	    getEntry(i, pt, mpf, ptbal, qn, qu, h, gm, gp, hn, hu,
		     eq0, eq1, eqn, equ);
	    getEntry(i, pm, q0, q1,     qn, qu, h, h0, h1, hn, hu,
		     eq0, eq1, eqn, equ);
	    getEntry(i, pd, r0, r1,     rn, ru, h, g0, g1, gn, gu,
		     er0, er1, ern, eru);
	    if (q0==0 && r0==0 && m=="mpfchs1") continue;
	    if (q1==0 && r1==0 && m=="ptchs") continue;
	    double pt = 0.5*(pm+pd);
	    cout << s << " pt="<<pt<<endl;
	    if (fabs(r1+rn+ru-r0)>1e-3) {
	      cout << "i="<<i<<" pt="<<pt<<" r1+rn+ru-r0="<<r1+rn+ru-r0<<endl
 		   << " r0="<<r0 << " r1="<<r1<<" rn="<<rn<<" ru="<<ru<<endl
		   <<flush;
	      //assert(fabs(r1+rn+ru-r0)<1e-3);
	    }
	    if (fabs(q1+qn+qu-q0)>1e-3) {
	      cout << "i="<<i<<" pt="<<pt<<" q1+qn+qu-q0="<<q1+qn+qu-q0<<endl
		   << " q0="<<q0 << " q1="<<q1<<" qn="<<qn<<" qu="<<qu<<endl
		   <<flush;
	      //assert(fabs(q1+qn+qu-q0)<1e-3);
	    }


	    // Master equations:
	    // MPF = r0 = r1 + rn + ru (1)
	    //        1 = p1 + pn + pu (2)
	    // r1 = p1 * R1 (3a)
	    // rn = pn * Rn (3b)
	    // ru = pu * Ru (3c)
	    // 
	    // Substitute (3x) to (1), then subtract (1) from (2)
	    // 1-r0 =    (1-R1)*p1       + (1-Rn)*pn + (1-Ru)*pu | +R1-1+r0
	    // R1 = r0 + (1-R1)*(p1-1)   + (1-Rn)*pn + (1-Ru)*pu | (2) for p1-1
	    // R1 = r0 + (1-R1)*(-pn-pu) + (1-Rn)*pn + (1-Ru)*pu | reorder
	    // R1 = r0 + (R1-Rn)*pn + (R1-Ru)*pu                 | (3b,3c)   
	    // R1 = r0 + [(R1-Rn)/Rn] * rn + [(R1-Ru)/Rn] * ru
	    //
	    // R1 = r0 + Dn*pn + Du*pu, where
	    // Dn = R1-Rn and Du = R1-Ru
	    // These are typically positive, as are pn and pu, so R1>r0
	    // R1 and Rn are known to 1%, so p1 and pn have uncertainty of 1%
	    // Ru on the other hand is much less precise (5%?), but probably
	    // more precise than pu (10%?) => pu = ru/Ru

	    // Z+jet and Z+g responses vs Z+q from drawZflavor.C
	    //double Rq = 1;
	    //double Ri_d = 1 + 0.01*(-1.48 - 0.40*pow(pt/100., -2.432));
	    //double Rg_d = 1 + 0.01*(-2.50 - 2.16*pow(pt/100., -1.586));
	    //double Rn_d = min(Ri_d, 0.5*(1+Rg_d)); // guess 50% q + 50% g
	    double Rn_d = 1;
	    //double R1 = 1.01; // Guesstimate for Z+jet MC truth response
	    //double Ru = 0.72; // R_UE(FullSim) for QCD UL17
	    //double Ru = 0.75; // R_UE(Detector) for QCD UL17
	    //double Ru = 0.80;//0.5*(0.85+0.75); // pT15/eta=0 and UE average
	    //double Ru = 0.85; // pT15/eta=0 jet response
	    //double Ru = 1.2; // test
	    //double Ru_d(0.72), Ru_m(0.72); // RUE(DY,HS1) or RUE(QCD,CP5)
	    //double Ru_d(0.65), Ru_m(0.65); // RUE(DY,CP5)
	    double Ru_d(0.685), Ru_m(0.685); // Ru(DY,CP5,UL17)
	    // Flavor response differences data-MC
	    //double dRq = 0;
	    //double dRi_d = 0.01*(-0.15);
	    //double dRg_d = 0.01*(-0.64);
	    
	    //double Ri_m = Ri_d - dRi_d;
	    //double Rg_m = Ri_d - dRi_d;
	    //double Rn_m = min(Ri_m, 0.5*(1+Rg_m)); // guess 50% q + 50% g
	    double Rn_m = 1;
	    /*
	    double Dn_d = Ri_d - Rn_d;
	    double Du_d = R1 - Ru_d;
	    double Pn_d = rn / (R1 - Dn_d);
	    double Pu_d = ru / (R1 - Du_d);
	    double P1_d = 1 - Pn_d - Pu_d;
	    //
	    double Dn_m = Ri_m - Rn_m;
	    double Du_m = R1 - Ru_m;
	    double Pn_m = qn / (R1 - Dn_m);
	    double Pu_m = qu / (R1 - Du_m);
	    double P1_m = 1 - Pn_m - Pu_m;
	    */

	    // Systematic variations for Rn_m, Ru_m numerically
	    double dRnd = -0.01; // guesstimate for data-MC
	    double dRnm = -0.02; // quark/gluon response variation
	    //double dRud = -0.05; // guesstimate for data-MC
	    //double dRum = -0.10; // RUE(DY,HS1)=0.72 / UE(DY,CP5)=0.65
	    //double dRud = +0.05; // guesstimate for data-MC
	    //double dRum = +0.10; // RUE(DY,HS1)=0.72 / UE(DY,CP5)=0.65
	    double dRud = +0.08; // QCD data=2.7 / FullSim(UL17)=2.5
	    double dRum = +0.05; // Ru(DY,HS1)=0.72 / Ru(DY,CP5)=0.65 vs 0.685

	    if (s=="multijet") {

	      Rn_m = 0.92; // low pT gluon jets
	      Rn_d = 0.90; // low pT gluon jets
	      Ru_d = Ru_m - 0.08;
	    }

	    if (m=="mpfchs1") {

	      // Ratio
	      //double dR = ((r0 + Dn_d*Pn_d + Du_d*Pu_d) /
	      //		 (q0 + Dn_m*Pn_m + Du_m*Pu_m)) - r0 / q0;
	      //double dRn = (Dn_d-Dn_m)*Pn_d + Dn_m*(Pn_d-Pn_m);
	      //double dRu = (Du_d-Du_m)*Pn_d + Du_m*(Pu_d-Pu_m);
	      
	      // Difference is better than ratio because factorization is exact
	      /*
	      double dR = ((r0 + Dn_d*Pn_d + Du_d*Pu_d) -
			   (q0 + Dn_m*Pn_m + Du_m*Pu_m)) - (r0 - q0);
	      */

	      // Solve master equation numerically
	      // R1 = r0 + [(R1-Rn)/Rn] * rn + [(R1-Ru)/Rn] * ru
	      // => fm = "x - ([0] + (x/[3]-1)*[1] + (x/[4]-1)*[2])"
	      fm->SetParameters(r0, rn, ru, Rn_d, Ru_d);
	      double R1_d = fm->GetX(0,0.50,1.50);
	      fm->SetParameters(q0, qn, qu, Rn_m, Ru_m);
	      double R1_m = fm->GetX(0,0.50,1.50);
	      double dR = R1_d / R1_m - r0 / q0;

	      fm->SetParameters(r0, rn, ru, Rn_d+dRnd, Ru_d);
	      double R1_dn1 = fm->GetX(0,0.50,1.50);
	      double dRn1 = (R1_dn1 - R1_d) / R1_m;
	      //
	      fm->SetParameters(r0, rn, ru, Rn_d+dRnm, Ru_d);
	      double R1_dn2 = fm->GetX(0,0.50,1.50);
	      fm->SetParameters(q0, qn, qu, Rn_m+dRnm, Ru_m);
	      double R1_mn2 = fm->GetX(0,0.50,1.50);
	      double dRn2 = R1_dn2/R1_mn2 - R1_d/R1_m; 
	      //
	      fm->SetParameters(r0, rn, ru, Rn_d, Ru_d+dRud);
	      double R1_du1 = fm->GetX(0,0.50,1.50);
	      double dRu1 = (R1_du1 - R1_d) / R1_m;
	      //
	      fm->SetParameters(r0, rn, ru, Rn_d, Ru_d+dRum);
	      double R1_du2 = fm->GetX(0,0.50,1.50);
	      fm->SetParameters(q0, qn, qu, Rn_m, Ru_m+dRum);
	      double R1_mu2 = fm->GetX(0,0.50,1.50);
	      double dRu2 = R1_du2/R1_mu2 - R1_d/R1_m; 

	      // Statistical variations with error propagation from r1, rn, ru
	      // (r0,r1,rn,ru correlations tricky, though)
	      double dR0s = sqrt(pow(er0,2) + pow(eq0,2));
	      double dR1s = sqrt(pow(er1,2) + pow(eq1,2));
	      double dRns = (R1_m/Rn_m - 1)*sqrt(pow(ern,2) + pow(eqn,2));
	      double dRus = (R1_m/Ru_m - 1)*sqrt(pow(eru,2) + pow(equ,2));

	      if (fabs(r0/q0-mpf)>1e-3) {
		cout <<s<<": r0/q0="<<r0/q0<<" mpf="<<mpf<<endl<<flush;
		assert(fabs(r0/q0-mpf)<1e-3);
	      }

	      /*
	      double Dn_d = R1_d - Rn_d;
	      double Du_d = R1_d - Ru_d;
	      double P1_d = r1 / R1_d;
	      double Pn_d = rn / Rn_d;
	      double Pu_d = ru / Ru_d;
	      //
	      double Dn_m = R1_m - Rn_m;
	      double Du_m = R1_m - Ru_m;
	      double P1_m = q1 / R1_m;
	      double Pn_m = qn / Rn_m;
	      double Pu_m = qu / Ru_m;	      

	      if (fabs(P1_d+Pn_d+Pu_d-1)>1e-3) {
		cout << "i="<<i<<" pt="<<pt
		     <<" P1+Pn+Pu-1="<<P1_d+Pn_d+Pu_d-1<<endl
		     << " (data) P1="<<P1_d<<" Pn="<<Pn_d<<" Pu="<<Pu_d<<endl
		     <<flush;
		//assert(fabs(P1_d+Pn_d+Pu_d-1)<1e-3);
	      }
	      if (fabs(P1_m+Pn_m+Pu_m-1)>1e-3) {
		cout << "i="<<i<<" pt="<<pt
		     <<" P1+Pn+Pu-1="<<P1_m+Pn_m+Pu_m-1<<endl
		     << " (MC) P1="<<P1_m<<" Pn="<<Pn_m<<" Pu="<<Pu_m<<endl
		     <<flush;
		//assert(fabs(P1_m+Pn_m+Pu_m-1)<1e-3);
	      }
	      
	      double dRn = (Dn_d-Dn_m)*Pn_d + Dn_m*(Pn_d-Pn_m);
	      double dRu = (Du_d-Du_m)*Pu_d + Du_m*(Pu_d-Pu_m);

	      hfsr->SetBinContent(i, dR);
	      hfsrn->SetBinContent(i, dRn);
	      hfsru->SetBinContent(i, dRu);
	      */

	      double ddR = sqrt(pow(dRn1,2) + pow(dRn2,2) +
				pow(dRu1,2)+pow(dRu2,2));
	      hfsr->SetBinContent(i, dR);
	      hfsr->SetBinError(i, ddR);
	      hfsr0s->SetBinContent(i, dR0s);
	      hfsr1s->SetBinContent(i, dR1s);
	      hfsrn1->SetBinContent(i, dRn1);
	      hfsrn2->SetBinContent(i, dRn2);
	      hfsrns->SetBinContent(i, dRns);
	      hfsru1->SetBinContent(i, dRu1);
	      hfsru2->SetBinContent(i, dRu2);
	      hfsrus->SetBinContent(i, dRus);
	   

	      cout << " MPF pt="<<pt<<" R1_d="<<R1_d<<" R1_m="<<R1_m
		   << " R1="<<R1_d/R1_m<<" raw="<<mpf<<endl;
	      cout << "  - r0="<<r0<<" r1="<<r1<<" rn="<<rn<<" ru="<<ru<<endl;
	      cout << "  - q0="<<q0<<" q1="<<q1<<" qn="<<qn<<" qu="<<qu<<endl;
	    } // mpfchs1

	    // For pT balance, correct estimate p1 based on pn and pu
	    if (m=="ptchs") {

	      //double dR = (r1/P1_d) / (q1/P1_m) - (r1 / q1);
	      // r1/q1*( (P1_m - P1_d)/P1_d )
	      // r1/q1*( (1 - Pn_m - Pu_m - (1 - Pn_d - Pu_d)) / P1_d)
	      // r1/q1*( ((Pn_d - Pn_m) + (Pu_d - Pu_m)) / P1_d)
	      //double dR = dRn + dRu;

	      // r1 = R1*p1 = R1*(1-pn-pu) = R1*(1-rn/Rn-ru/Ru)
	      // R1 = r1 + [R1/Rn] * rn + [R1/Ru] * ru
	      // => fp = "x - ([0] + x/[3]*[1] + x/[4]*[2])"
	      double xmin(0.1), xmax(1.9);
	      fp->SetParameters(r1, rn, ru, Rn_d, Ru_d);
	      double R1_d = fp->GetX(0,xmin,xmax);
	      fp->SetParameters(q1, qn, qu, Rn_m, Ru_m);
	      double R1_m = fp->GetX(0,xmin,xmax);
	      double dR = R1_d / R1_m - r1 / q1;

	      fp->SetParameters(r1, rn, ru, Rn_d+dRnd, Ru_d);
	      double R1_dn1 = fp->GetX(0,xmin,xmax);
	      double dRn1 = (R1_dn1 - R1_d) / R1_m;
	      //
	      fp->SetParameters(r1, rn, ru, Rn_d+dRnm, Ru_d);
	      double R1_dn2 = fp->GetX(0,xmin,xmax);
	      fp->SetParameters(q1, qn, qu, Rn_m+dRnm, Ru_m);
	      double R1_mn2 = fp->GetX(0,xmin,xmax);
	      double dRn2 = R1_dn2/R1_mn2 - R1_d/R1_m; 
	      //
	      fp->SetParameters(r1, rn, ru, Rn_d, Ru_d+dRud);
	      double R1_du1 = fp->GetX(0,xmin,xmax);
	      double dRu1 = (R1_du1 - R1_d) / R1_m;
	      //
	      fp->SetParameters(r1, rn, ru, Rn_d, Ru_d+dRum);
	      double R1_du2 = fp->GetX(0,xmin,xmax);
	      fp->SetParameters(q1, qn, qu, Rn_m, Ru_m+dRum);
	      double R1_mu2 = fp->GetX(0,xmin,xmax);
	      double dRu2 = R1_du2/R1_mu2 - R1_d/R1_m; 

	      // Statistical variations with error propagation from r1, rn, ru
	      // (r0,r1,rn,ru correlations tricky, though)
	      double dR0s = sqrt(pow(er0,2) + pow(eq0,2));
	      double dR1s = sqrt(pow(er1,2) + pow(eq1,2));
	      double dRns = R1_m/Rn_m*sqrt(pow(ern,2) + pow(eqn,2));
	      double dRus = R1_m/Ru_m*sqrt(pow(eru,2) + pow(equ,2));

	      if (fabs(r1/q1-ptbal)>1e-3) {
		cout << "*i="<<i<<" pt="<<pt
		     <<" r1/q1="<<r1/q1<<" mpf="<<mpf<<endl<<flush;
		//assert(fabs(r1/q1-ptbal)<1e-3);
	      }

	      // temp check
	      //fm->SetParameters(r0, rn, ru, Rn_d, Ru_d);
	      //double R1m_d = fm->GetX(0,0.50,1.50);
	      //fm->SetParameters(q0, qn, qu, Rn_m, Ru_m);
	      //double R1m_m = fm->GetX(0,0.50,1.50);
	      //dR = R1_d / R1_m - r1 / q1;
	      // temp check
	      /*
	      double P1_d = r1 / R1_d;
	      double Pn_d = rn / Rn_d;
	      double Pu_d = ru / Ru_d;
	      //
	      double Pn_m = qn / Rn_m;
	      double Pu_m = qu / Ru_m;
	      
	      double dRn = r1/q1 * (Pn_d - Pn_m) / P1_d;
	      double dRu = r1/q1 * (Pu_d - Pu_m) / P1_d;
	      */
	      //hfsr->SetBinContent(i, dR);
	      //hfsrn->SetBinContent(i, dRn);
	      //hfsru->SetBinContent(i, dRu);

	      double ddR = sqrt(pow(dRn1,2) + pow(dRn2,2) +
				pow(dRu1,2)+pow(dRu2,2));
	      hfsr->SetBinContent(i, dR);
	      hfsr->SetBinError(i, ddR);
	      hfsr0s->SetBinContent(i, dR0s);
	      hfsr1s->SetBinContent(i, dR1s);
	      hfsrn1->SetBinContent(i, dRn1);
	      hfsrn2->SetBinContent(i, dRn2);
	      hfsrns->SetBinContent(i, dRns);
	      hfsru1->SetBinContent(i, dRu1);
	      hfsru2->SetBinContent(i, dRu2);
	      hfsrus->SetBinContent(i, dRus);

	      cout << " pTbal pt="<<pt<<" R1_d="<<R1_d<<" R1_m="<<R1_m
		   << " R1="<<R1_d/R1_m<<" raw="<<ptbal<<endl;
	      cout << "  - r0="<<r0<<" r1="<<r1<<" rn="<<rn<<" ru="<<ru<<endl;
	      cout << "  - q0="<<q0<<" q1="<<q1<<" qn="<<qn<<" qu="<<qu<<endl;
	    } // ptchs
	  } // for i

	  finout->cd(Form("%s/eta%02d-%02d/fsr",
			  cd,int(10*etamin),int(10*etamax)));
	  hfsr->Write(Form("hkfsr3_%s_%s",cm,cs),TObject::kOverwrite);
	  hfsr0s->Write(Form("hkfsr3_%s_%s_mpf0s",cm,cs),TObject::kOverwrite);
	  hfsr1s->Write(Form("hkfsr3_%s_%s_mpf1s",cm,cs),TObject::kOverwrite);
	  hfsrn1->Write(Form("hkfsr3_%s_%s_mpfn1",cm,cs),TObject::kOverwrite);
	  hfsrn2->Write(Form("hkfsr3_%s_%s_mpfn2",cm,cs),TObject::kOverwrite);
	  hfsrns->Write(Form("hkfsr3_%s_%s_mpfns",cm,cs),TObject::kOverwrite);
	  hfsru1->Write(Form("hkfsr3_%s_%s_mpfu1",cm,cs),TObject::kOverwrite);
	  hfsru2->Write(Form("hkfsr3_%s_%s_mpfu2",cm,cs),TObject::kOverwrite);
	  hfsrus->Write(Form("hkfsr3_%s_%s_mpfus",cm,cs),TObject::kOverwrite);
	  curdir->cd();

	} // mpfchs

	// mpfchs1 needs to be loaded (and corrected) before ptchs
	//if (m=="ptchs" && d=="ratio) {
	//assert(gs[d]["mpfchs1"][s]!=0);

	//} // ptchs
      } // isam
    } // imet
  } // idir

  finout->Close();
  curdir->cd();
} // softrad2
