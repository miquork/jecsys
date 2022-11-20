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

#include "Flavor.h"
#include "tdrstyle_mod15.C"

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>

// Scale soft jets and unclustered energy with gJES from Z+flavor (Flavor.C)
// gJESpt=15 => 0.947, gJESpt=45 => 0.981
//bool useGluonJES = false;// DP_2021 // (def:true)
bool useGluonJES = true; // (def:true)
double gJESpt = 15.; // reference pT for gluonJES (def:45)

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
Flavor *_flv(0);
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
  bool isUL16 = (epoch=="2016BCDEF" || epoch=="2016BCD" ||
		 epoch=="2016EF" || epoch=="2016GH");
  bool isAPV = (epoch=="2016BCDEF" || epoch=="2016BCD" || epoch=="2016EF");
  bool isRun2 = (epoch=="Run2Test");

  if (!_flv) {
    cout << "Initializing Flavor" << endl;
    _flv = new Flavor();
    _flv->getResp(gJESpt,0.,"MultijetRecoil25",-1,"25"); // Initialization
    cout << Form("gluonJES(pT=%1.0f,|eta|<2.5) = ",gJESpt)
	 << _flv->getResp(gJESpt,0.,"MultijetRecoil25",-1,"25") << endl;
  }

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
  dirs2.push_back("data");
  dirs2.push_back("mc");

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
  if (etamax==1.3) {
    //if (!isUL16) 
    if (epoch!="2017H") {
      samples.push_back("multijet");
      samples.push_back("zlljet");
      samples.push_back("zeejet");
      samples.push_back("zmmjet");
    }
    if (epoch=="2016EF" || epoch=="2016GH" || epoch=="2016BCD" ||
	epoch=="2018ABCD" || epoch=="2017BCDEF" || epoch=="2016BCDEF" ||
	epoch=="Run2Test")
      samples.push_back("gamjet");
    if (isUL18||isUL17||isUL16||isRun2) {
      samples.push_back("gi");
      samples.push_back("gb");
      samples.push_back("gc");
      samples.push_back("gq");
      samples.push_back("gg");
      samples.push_back("gn");
      //
      samples.push_back("gii");
      samples.push_back("gib");
      samples.push_back("gic");
      samples.push_back("giq");
      samples.push_back("gig");
      samples.push_back("gin");
      //
      samples.push_back("gbi");
      samples.push_back("gbb");
      samples.push_back("gbc");
      samples.push_back("gbq");
      samples.push_back("gbg");
      samples.push_back("gbn");
      //
      samples.push_back("gci");
      samples.push_back("gcb");
      samples.push_back("gcc");
      samples.push_back("gcq");
      samples.push_back("gcg");
      samples.push_back("gcn");
      //
      samples.push_back("gqi");
      samples.push_back("gqb");
      samples.push_back("gqc");
      samples.push_back("gqq");
      samples.push_back("gqg");
      samples.push_back("gqn");
      //
      samples.push_back("ggi");
      samples.push_back("ggb");
      samples.push_back("ggc");
      samples.push_back("ggq");
      samples.push_back("ggg");
      samples.push_back("ggn");
      //
      samples.push_back("gni");
      samples.push_back("gnb");
      samples.push_back("gnc");
      samples.push_back("gnq");
      samples.push_back("gng");
      samples.push_back("gnn");
    }
  }
  samples.push_back("zjet");
  if (isUL18||isUL17||isUL16||isRun2) {
  samples.push_back("zi");
  samples.push_back("zb");
  samples.push_back("zc");
  samples.push_back("zq");
  samples.push_back("zg");
  samples.push_back("zn");
  //
  samples.push_back("zii");
  samples.push_back("zib");
  samples.push_back("zic");
  samples.push_back("ziq");
  samples.push_back("zig");
  samples.push_back("zin");
  //
  samples.push_back("zbi");
  samples.push_back("zbb");
  samples.push_back("zbc");
  samples.push_back("zbq");
  samples.push_back("zbg");
  samples.push_back("zbn");
  //
  samples.push_back("zci");
  samples.push_back("zcb");
  samples.push_back("zcc");
  samples.push_back("zcq");
  samples.push_back("zcg");
  samples.push_back("zcn");
  //
  samples.push_back("zqi");
  samples.push_back("zqb");
  samples.push_back("zqc");
  samples.push_back("zqq");
  samples.push_back("zqg");
  samples.push_back("zqn");
  //
  samples.push_back("zgi");
  samples.push_back("zgb");
  samples.push_back("zgc");
  samples.push_back("zgq");
  samples.push_back("zgg");
  samples.push_back("zgn");
  //
  samples.push_back("zni");
  samples.push_back("znb");
  samples.push_back("znc");
  samples.push_back("znq");
  samples.push_back("zng");
  samples.push_back("znn");
  }

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

	bool isflavormc = 
	  (s=="zii" || s=="zbi" || s=="zci" || s=="zqi" || s=="zgi" ||s=="zni"||
	   s=="zib" || s=="zbb" || s=="zcb" || s=="zqb" || s=="zgb" ||s=="znb"||
	   s=="zic" || s=="zbc" || s=="zcc" || s=="zqc" || s=="zgc" ||s=="znc"||
	   s=="ziq" || s=="zbq" || s=="zcq" || s=="zqq" || s=="zgq" ||s=="znq"||
	   s=="zig" || s=="zbg" || s=="zcg" || s=="zqg" || s=="zgg" ||s=="zng"||
	   s=="zin" || s=="zbn" || s=="zcn" || s=="zqn" || s=="zgn" ||s=="znn"||
	   s=="gii" || s=="gbi" || s=="gci" || s=="gqi" || s=="ggi" ||s=="gni"||
	   s=="gib" || s=="gbb" || s=="gcb" || s=="gqb" || s=="ggb" ||s=="gnb"||
	   s=="gic" || s=="gbc" || s=="gcc" || s=="gqc" || s=="ggc" ||s=="gnc"||
	   s=="giq" || s=="gbq" || s=="gcq" || s=="gqq" || s=="ggq" ||s=="gnq"||
	   s=="gig" || s=="gbg" || s=="gcg" || s=="gqg" || s=="ggg" ||s=="gng"||
	   s=="gin" || s=="gbn" || s=="gcn" || s=="gqn" || s=="ggn" ||s=="gnn");
	if (isflavormc && d!="mc") continue;

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

  // Helper function for solving jet response x from MPF decomposition
  TF1 *fm = new TF1("fm","x - ([0] + (x/[3]-1)*[1] + (x/[4]-1)*[2])",0,13000);
  TF1 *fp = new TF1("fp","x - ([0] + (x/[3]  )*[1] + (x/[4]  )*[2])",0,13000);

  // Helper functions for multijet with recoil binning, x=R1
  // [0]=MPF or MJB with lead in denomator, [1]=n,[2]=u, [3]=Rn,[4]=Ru,[5]=R2
  /*
  TF1 *fml = new TF1("fml","[0] - ([5]/x + (1-[5]/[3])*[1] + (1-[5]/[4])*[2])",
		     0,13000);
  TF1 *fpl = new TF1("fpl","[0] - ([5]/x + ( -[5]/[3])*[1] + ( -[5]/[4])*[2])",
		     0,13000);
  */
  // Helper function for solving jet response ratio x from multijet MPF/MJB
  // in 'ptave' binning
  // [0]=MPF or MJB, [1]=n, [2]=u, [3]=Rn, [4]=Ru, [5]=R2
  /*
  TF1 *fma = new TF1("fma","(x-1)/(x+1)-[0]"
		     " + ((1./[5])/(x+1) - 1./(2.*[3])) * [1]"
		     " + ((1./[5])/(x+1) - 1./(2.*[4])) * [2]",0,13000);
  TF1 *fpa = new TF1("fma","(x-1)/(x+1)-[0]"
		     " + (0              - 1./(2.*[3])) * [1]"
		     " + (0              - 1./(2.*[4])) * [2]",0,13000);
  */
  TF1 *fma = new TF1("fma","(x-1)/(x+1)-[0]"
		     " + (1 - (x+1)/(2.*[3]/[5])) * [1]"
		     " + (1 - (x+1)/(2.*[4]/[5])) * [2]",0,13000);
  TF1 *fpa = new TF1("fma","(x-1)/(x+1)-[0]"
		     " + (0 - (x+1)/(2.*[3]/[5])) * [1]"
		     " + (0 - (x+1)/(2.*[4]/[5])) * [2]",0,13000);

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

	bool isflavormc =
	  (s=="zii" || s=="zbi" || s=="zci" || s=="zqi" || s=="zgi" ||s=="zni"||
	   s=="zib" || s=="zbb" || s=="zcb" || s=="zqb" || s=="zgb" ||s=="znb"||
	   s=="zic" || s=="zbc" || s=="zcc" || s=="zqc" || s=="zgc" ||s=="znc"||
	   s=="ziq" || s=="zbq" || s=="zcq" || s=="zqq" || s=="zgq" ||s=="znq"||
	   s=="zig" || s=="zbg" || s=="zcg" || s=="zqg" || s=="zgg" ||s=="zng"||
	   s=="zin" || s=="zbn" || s=="zcn" || s=="zqn" || s=="zgn" ||s=="znn"||
	   s=="gii" || s=="gbi" || s=="gci" || s=="gqi" || s=="ggi" ||s=="gni"||
	   s=="gib" || s=="gbb" || s=="gcb" || s=="gqb" || s=="ggb" ||s=="gnb"||
	   s=="gic" || s=="gbc" || s=="gcc" || s=="gqc" || s=="ggc" ||s=="gnc"||
	   s=="giq" || s=="gbq" || s=="gcq" || s=="gqq" || s=="ggq" ||s=="gnq"||
	   s=="gig" || s=="gbg" || s=="gcg" || s=="gqg" || s=="ggg" ||s=="gng"||
	   s=="gin" || s=="gbn" || s=="gcn" || s=="gqn" || s=="ggn" ||s=="gnn");

	// mpf1, mpfn and mpfu must be loaded before correcting mpfchs1
	// statistics needed for checking pT binning
	if ((m=="mpfchs1" || m=="ptchs")) { // && d=="ratio") {
	  assert(hs["data"][s]!=0 || isflavormc);// || s=="multijet");
	  assert(gs["mc"]["mpf1"][s]!=0);
	  assert(gs["mc"]["mpfn"][s]!=0);
	  assert(gs["mc"]["mpfu"][s]!=0);
	  assert(gs["mc"]["mpfchs1"][s]!=0);
	  assert(gs["data"]["mpf1"][s]!=0 || isflavormc);
	  assert(gs["data"]["mpfn"][s]!=0 || isflavormc);
	  assert(gs["data"]["mpfu"][s]!=0 || isflavormc);
	  assert(gs["data"]["mpfchs1"][s]!=0 || isflavormc);
	  assert(gs["ratio"]["mpfchs1"][s]!=0 || isflavormc);
	  assert(gs["ratio"]["ptchs"][s]!=0 || isflavormc);

	  TH1D *h = hs["data"][s];
	  TGraphErrors *g0 = gs["data"]["mpfchs1"][s];
	  TGraphErrors *g1 = gs["data"]["mpf1"][s];
	  TGraphErrors *gn = gs["data"]["mpfn"][s];
	  TGraphErrors *gu = gs["data"]["mpfu"][s];
	  if (isflavormc) { // placeholder for missing data
	    h = hs["mc"][s];
	    g0 = (TGraphErrors*)gs["mc"]["mpfchs1"][s]->Clone();
	    g1 = (TGraphErrors*)gs["mc"]["mpf1"][s]->Clone();
	    gn = (TGraphErrors*)gs["mc"]["mpfn"][s]->Clone();
	    gu = (TGraphErrors*)gs["mc"]["mpfu"][s]->Clone();
	  }
	  TGraphErrors *h0 = gs["mc"]["mpfchs1"][s];
	  TGraphErrors *h1 = gs["mc"]["mpf1"][s];
	  TGraphErrors *hn = gs["mc"]["mpfn"][s];
	  TGraphErrors *hu = gs["mc"]["mpfu"][s];
	  //
	  TGraphErrors *gm = gs["ratio"]["mpfchs1"][s];
	  TGraphErrors *gp = gs["ratio"]["ptchs"][s];
	  if (isflavormc) { // placeholder for missing data ratio
	    gm = (TGraphErrors*)gs["mc"]["mpfchs1"][s]->Clone();
	    gp = (TGraphErrors*)gs["mc"]["ptchs"][s]->Clone();
	    for (int i = 0; i != gm->GetN(); ++i)
	      gm->SetPoint(i, gm->GetX()[i], 1);
	    for (int i = 0; i != gp->GetN(); ++i)
	      gp->SetPoint(i, gp->GetX()[i], 1);
	  }


	  TH1D *hdm    = (TH1D*)h->Clone("hdm");   hdm->Reset();
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
	    //if ((q0==0 || r0==0) && m=="mpfchs1") continue;
	    //if ((q1==0 || r1==0) && m=="ptchs") continue;
	    double pt = (isflavormc ? pm : 0.5*(pm+pd));
	    cout << s << " " << m << " pt="<<pt<<endl;
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

	    // start latest ==>
	    double Rn_d(1.000), Rn_m(1.000);
	    // v1: Runcl = 0.685
	    //double Ru_d(0.685), Ru_m(0.685); // Ru(DY,CP5,UL17)
	    // v2: Runcl = 0.92
	    double Ru_d(0.92), Ru_m(0.92); // soft gluon jets in MC
	    // v3: Runcl_data * gJES, Rn_data * gJES (optional)
	    double R2_m(1.000), R2_d(1.000); // multijet recoil?
	    if (useGluonJES) {
	      double gJES = _flv->getResp(gJESpt,0.,"MultijetRecoil25",-1,"25");
	      Ru_d *= gJES;
	    }
	    
	    // Systematic variations for Rn_m, Ru_m numerically
	    double dRnd = -0.01; // guesstimate for data-MC
	    double dRnm = -0.02; // quark/gluon response variation
	    // rewind two:
	    double dRud = +0.08; // QCD data=2.7 / FullSim(UL17)=2.5 //
	    //double dRum = +0.05; // Ru(DY,HS1)=0.72 / Ru(DY,CP5)=0.65 vs 0.685
	    // latest two:
	    //double dRud = -0.02; // -gJES or gJESpt=45->15 extrapolation
	    double dRum = -0.10; // ca. (R_ue-Ru)/2 or Ru_m-1 towards R_ue
	    // ==> end latest

	    /*
	    // start DP_2021 ==>
	    double R2_m(1.000), R2_d(1.000); // multijet recoil?
	    double Rn_d(1.000), Rn_m(1.000);
            double Ru_d(0.685), Ru_m(0.685); // Ru(DY,CP5,UL17)

            // Systematic variations for Rn_m, Ru_m numerically
            double dRnd = -0.01; // guesstimate for data-MC
            double dRnm = -0.02; // quark/gluon response variation
            double dRud = +0.08; // QCD data=2.7 / FullSim(UL17)=2.5
            double dRum = +0.05; // Ru(DY,HS1)=0.72 / Ru(DY,CP5)=0.65 vs 0.685
	    // ==> end DP2021
	    */

	    if (s=="multijet") {

	      //Rn_m = 0.92; // low pT gluon jets
	      //Rn_d = 0.90; // low pT gluon jets
	      //Ru_d = Ru_m - 0.08;
	      //Rn_m = 1.000; // DP_2021
	      //Rn_d = 1.000; // DP_2021
	      // v1: Runcl = 0.65
	      //Ru_m = 0.65; // DP_2021 // from minitools/drawMultijetMPB_ISRonly.pdf
	      //Ru_d = 0.65; // DP_2021
	      // v2: Runcl = 0.92
	      //Ru_m = Ru_d = 0.92; // low pT gluon jets in MC
	      //double R2_m = 1.000; // DP_2021
	      //double R2_d = 1.000; // DP_2021

	      // turn <Rlead>/<Rrecoil> ratio back to <A> and <B>
	      double Ad = (r0-1)/(r0+1);
	      double Bd = (r1-1)/(r1+1);
	      double Am = (q0-1)/(q0+1);
	      double Bm = (q1-1)/(q1+1);

	      // same for N-jet and unclustered, which are
	      // rn = 2*Pn / (1 - Pn) <=> rn-rn*Pn=+2*Pn <=> Pn=rn/(2+rn)
	      // ru = -2*Pu / (1 - Pu) <=> ru-ru*Pu=-2*Pu <=> Pu=-ru/(2-ru)
	      double k = 1; // extra fudge factor => fixed in fma, fpa
	      double Pnd = k*rn/(2+rn);
	      //double Pud = +k*ru/(2-ru); // Orig. PU had extra minus sign
	      double Pud = k*ru/(2+ru); // Now fixed?
	      double Pnm = k*qn/(2+qn);
	      //double Pum = +k*qu/(2-qu); // Orig. PU had extra minus sign
	      double Pum = k*qu/(2+qu); // Now fixed?

	      double Rd(0), Rm(0), dR(0), dR_d(0), dR_m(0);
	      double xmin(0.1), xmax(1.9);
	      if (m=="mpfchs1") { // multijet

		fma->SetParameters(Ad,Pnd,Pud,Rn_d,Ru_d,R2_d);
		Rd = fma->GetX(0,xmin,xmax);
		fma->SetParameters(Am,Pnm,Pum,Rn_m,Ru_m,R2_m);
		Rm = fma->GetX(0,xmin,xmax);
		dR = Rd / Rm - r0 / q0;
		
		// Approximate parts for data and MC
		dR_d = Rd - r0;
		dR_m = Rm - q0;

		// vary N-jet Rn in data (data/MC difference)
		fma->SetParameters(Ad,Pnd,Pud,Rn_d+dRnd,Ru_d,R2_d);
		double Rd_n1 = fma->GetX(0,xmin,xmax);
		double dRn1 = (Rd_n1 - Rd) / Rm;
		// vary N-jet Rn in both data and MC (MC reference scale)
		fma->SetParameters(Ad,Pnd,Pud,Rn_d+dRnm,Ru_d,R2_d);
		double Rd_n2 = fma->GetX(0,xmin,xmax);
		fma->SetParameters(Am,Pnm,Pum,Rn_m+dRnm,Ru_m,R2_m);
		double Rm_n2 = fma->GetX(0,xmin,xmax);
		double dRn2 = Rd_n2/Rm_n2 - Rd/Rm; 
		// vary unclustered Ru in data (data/MC difference)
		fma->SetParameters(Ad,Pnd,Pud,Rn_d,Ru_d+dRud,R2_d);
		double Rd_u1 = fma->GetX(0,xmin,xmax);
		double dRu1 = (Rd_u1 - Rd) / Rm;
		// vary unclustered Ru in both data and MC (MC reference scale)
		fma->SetParameters(Ad,Pnd,Pud,Rn_d,Ru_d+dRum,R2_d);
		double Rd_u2 = fma->GetX(0,xmin,xmax);
		fma->SetParameters(Am,Pnm,Pum,Rn_m,Ru_m+dRum,R2_m);
		double Rm_u2 = fma->GetX(0,xmin,xmax);
		double dRu2 = Rd_u2/Rm_u2 - Rd/Rm; 

		// Statistical variations with error propagation from r1, rn, ru
		// (r0,r1,rn,ru correlations tricky, though)
		// ...Not yet adapted from Z+jet to multijet...
		double dR0s = sqrt(pow(er0,2) + pow(eq0,2));
		double dR1s = sqrt(pow(er1,2) + pow(eq1,2));
		double dRns = (Rm/Rn_m - 1)*sqrt(pow(ern,2) + pow(eqn,2));
		double dRus = (Rm/Ru_m - 1)*sqrt(pow(eru,2) + pow(equ,2));

		if (fabs(r0/q0-mpf)>1e-3) {
		  cout <<s<<": r0/q0="<<r0/q0<<" mpf="<<mpf<<endl<<flush;
		  assert(fabs(r0/q0-mpf)<1e-3);
		}
		
		double ddR = sqrt(pow(dRn1,2) + pow(dRn2,2) +
				  pow(dRu1,2)+pow(dRu2,2));

		if (d=="data")    hdm->SetBinContent(i, Rd);
		if (d=="data")    hdm->SetBinError(i, er0);
		if (Rm>0 || Rm<0) { // Safety againt Inf/NaN
		  if (d=="ratio") hdm->SetBinContent(i, Rd / Rm);
		  if (d=="mc")    hdm->SetBinContent(i, Rm);
		  // Statistical uncertainty will need more careful handling
		  if (d=="ratio") hdm->SetBinError(i, sqrt(er0*er0 + eq0*eq0));
		  if (d=="mc")    hdm->SetBinError(i, eq0);
		}

		if (d=="data") hfsr->SetBinContent(i, dR_d);
		if (dR_m>0 || dR_m<0) {// Safety against Inf/NaN
		  if (d=="ratio") hfsr->SetBinContent(i, dR);
		  if (d=="mc")  hfsr->SetBinContent(i, dR_m);
		  hfsr->SetBinError(i, ddR);
		  hfsr0s->SetBinContent(i, dR0s);
		  hfsr1s->SetBinContent(i, dR1s);
		  hfsrn1->SetBinContent(i, dRn1);
		  hfsrn2->SetBinContent(i, dRn2);
		  hfsrns->SetBinContent(i, dRns);
		  hfsru1->SetBinContent(i, dRu1);
		  hfsru2->SetBinContent(i, dRu2);
		  hfsrus->SetBinContent(i, dRus);
		}		

	      } // mpfchs1 multijet

	      if (m=="ptchs") { // multijet

		fpa->SetParameters(Bd,Pnd,Pud,Rn_d,Ru_d,R2_d);
		Rd = fpa->GetX(0,xmin,xmax);
		fpa->SetParameters(Bm,Pnm,Pum,Rn_m,Ru_m,R2_m);
		Rm = fpa->GetX(0,xmin,xmax);
		dR = Rd / Rm - r1 / q1;

		dR_d = Rd - r1;
		dR_m = Rm - q1;

		// vary N-jet Rn in data (data/MC difference)
		fpa->SetParameters(Bd,Pnd,Pud,Rn_d+dRnd,Ru_d,R2_d);
		double Rd_n1 = fpa->GetX(0,xmin,xmax);
		double dRn1 = (Rd_n1 - Rd) / Rm;
		// vary N-jet Rn in both data and MC (MC reference scale)
		fpa->SetParameters(Bd,Pnd,Pud,Rn_d+dRnm,Ru_d,R2_d);
		double Rd_n2 = fpa->GetX(0,xmin,xmax);
		fpa->SetParameters(Bm,Pnm,Pum,Rn_m+dRnm,Ru_m,R2_m);
		double Rm_n2 = fpa->GetX(0,xmin,xmax);
		double dRn2 = Rd_n2/Rm_n2 - Rd/Rm; 
		// vary unclustered Ru in data (data/MC difference)
		fpa->SetParameters(Bd,Pnd,Pud,Rn_d,Ru_d+dRud,R2_d);
		double Rd_u1 = fpa->GetX(0,xmin,xmax);
		double dRu1 = (Rd_u1 - Rd) / Rm;
		// vary unclustered Ru in both data and MC (MC reference scale)
		fpa->SetParameters(Bd,Pnd,Pud,Rn_d,Ru_d+dRum,R2_d);
		double Rd_u2 = fpa->GetX(0,xmin,xmax);
		fpa->SetParameters(Bm,Pnm,Pum,Rn_m,Ru_m+dRum,R2_m);
		double Rm_u2 = fpa->GetX(0,xmin,xmax);
		double dRu2 = Rd_u2/Rm_u2 - Rd/Rm; 

		// Statistical variations with error propagation from r1, rn, ru
		// (r0,r1,rn,ru correlations tricky, though)
		// ...Not yet adapted from Z+jet to multijet...
		double dR0s = sqrt(pow(er0,2) + pow(eq0,2));
		double dR1s = sqrt(pow(er1,2) + pow(eq1,2));
		double dRns = (Rm/Rn_m - 1)*sqrt(pow(ern,2) + pow(eqn,2));
		double dRus = (Rm/Ru_m - 1)*sqrt(pow(eru,2) + pow(equ,2));

		if (fabs(r0/q0-mpf)>1e-3) {
		  cout <<s<<": r0/q0="<<r0/q0<<" mpf="<<mpf<<endl<<flush;
		  assert(fabs(r0/q0-mpf)<1e-3);
		}
		
		double ddR = sqrt(pow(dRn1,2) + pow(dRn2,2) +
				  pow(dRu1,2)+pow(dRu2,2));

		if (d=="data")  hdm->SetBinContent(i, Rd);
		if (d=="data")  hdm->SetBinError(i, er0);
		if (Rm>0 || Rm<0) { // Safety against Inf/NaN
		  if (d=="ratio") hdm->SetBinContent(i, Rd / Rm);
		  if (d=="mc")    hdm->SetBinContent(i, Rm);
		  // Statistical uncertainty will need more careful handling
		  if (d=="ratio") hdm->SetBinError(i, sqrt(er0*er0 + eq0*eq0));
		  if (d=="mc")    hdm->SetBinError(i, eq0);
		}

		if (d=="data")  hfsr->SetBinContent(i, dR_d);
		if (dR>0 || dR<0) { // Safety against Inf/NaN
		  if (d=="ratio") hfsr->SetBinContent(i, dR);
		  if (d=="mc")    hfsr->SetBinContent(i, dR_m);
		  hfsr->SetBinError(i, ddR);
		  hfsr0s->SetBinContent(i, dR0s);
		  hfsr1s->SetBinContent(i, dR1s);
		  hfsrn1->SetBinContent(i, dRn1);
		  hfsrn2->SetBinContent(i, dRn2);
		  hfsrns->SetBinContent(i, dRns);
		  hfsru1->SetBinContent(i, dRu1);
		  hfsru2->SetBinContent(i, dRu2);
		  hfsrus->SetBinContent(i, dRus);
		}
	      } // ptchs multijet

	      cout << "r0="<<r0<<" r1="<<r1<<" q0="<<q0<<" q1="<<q1<<endl;
	      cout << "Ad="<<Ad<<" Bd="<<Bd<<" Am="<<Am<<" Bm="<<Bm<<endl;
	      cout << "Pnd="<<Pnd<<" Pud="<<Pud
		   <<" Pnm="<<Pnm<<" Pum="<<Pum<<endl;
	      cout << "Rd="<<Rd<<" Rm="<<Rm<<endl;
	    } // multijet
	    

	    // xmin(0.5), xmax(1.5)
	    double xmin(0.1), xmax(1.9);
	    if (m=="mpfchs1" && s!="multijet") {

	      // Solve master equation numerically
	      // R1 = r0 + [(R1-Rn)/Rn] * rn + [(R1-Ru)/Rn] * ru
	      // => fm = "x - ([0] + (x/[3]-1)*[1] + (x/[4]-1)*[2])"
	      fm->SetParameters(r0, rn, ru, Rn_d, Ru_d);
	      double R1_d = fm->GetX(0,xmin,xmax);
	      fm->SetParameters(q0, qn, qu, Rn_m, Ru_m);
	      double R1_m = fm->GetX(0,xmin,xmax);
	      double dR = R1_d / R1_m - r0 / q0;
	      double dR_d = R1_d - r0;
	      double dR_m = R1_m - q0;

	      fm->SetParameters(r0, rn, ru, Rn_d+dRnd, Ru_d);
	      double R1_dn1 = fm->GetX(0,xmin,xmax);
	      double dRn1 = (R1_dn1 - R1_d) / R1_m;
	      //
	      fm->SetParameters(r0, rn, ru, Rn_d+dRnm, Ru_d);
	      double R1_dn2 = fm->GetX(0,xmin,xmax);
	      fm->SetParameters(q0, qn, qu, Rn_m+dRnm, Ru_m);
	      double R1_mn2 = fm->GetX(0,xmin,xmax);
	      double dRn2 = R1_dn2/R1_mn2 - R1_d/R1_m; 
	      //
	      fm->SetParameters(r0, rn, ru, Rn_d, Ru_d+dRud);
	      double R1_du1 = fm->GetX(0,xmin,xmax);
	      double dRu1 = (R1_du1 - R1_d) / R1_m;
	      //
	      fm->SetParameters(r0, rn, ru, Rn_d, Ru_d+dRum);
	      double R1_du2 = fm->GetX(0,xmin,xmax);
	      fm->SetParameters(q0, qn, qu, Rn_m, Ru_m+dRum);
	      double R1_mu2 = fm->GetX(0,xmin,xmax);
	      double dRu2 = R1_du2/R1_mu2 - R1_d/R1_m; 

	      // Statistical variations with error propagation from r1, rn, ru
	      // (r0,r1,rn,ru correlations tricky, though)
	      double dR0s = sqrt(pow(er0,2) + pow(eq0,2));
	      double dR1s = sqrt(pow(er1,2) + pow(eq1,2));
	      double dRns = (R1_m/Rn_m - 1)*sqrt(pow(ern,2) + pow(eqn,2));
	      double dRus = (R1_m/Ru_m - 1)*sqrt(pow(eru,2) + pow(equ,2));

	      if (fabs(r0/q0-mpf)>1e-3) {
		cout <<s<<": r0/q0="<<r0/q0<<" mpf="<<mpf
		     << " diff="<<fabs(r0/q0-mpf)<<endl<<flush;
		//assert(fabs(r0/q0-mpf)<1e-3);
		assert(fabs(r0/q0-mpf)<3.9e-3 || mpf==0); // Run2Test
		//assert(fabs(r0/q0-mpf)<3.1e-3 || mpf==0); // Eta13ee_v2
		//assert(fabs(r0/q0-mpf)<2.3e-3); // Eta13ee_v1
		//assert(fabs(r0/q0-mpf)<1.2e-3); // Eta25ee
		//assert(fabs(r0/q0-mpf)<1e-3 ||
		//     (fabs(r0/q0-mpf)<8e-3  &&s=="zlljet"));
	      }

	      double ddR = sqrt(pow(dRn1,2) + pow(dRn2,2) +
				pow(dRu1,2)+pow(dRu2,2));

	      if (d=="data")  hdm->SetBinContent(i, R1_d);
	      if (d=="data")  hdm->SetBinError(i, er0);
	      if (R1_m>0 || R1_m<0) { // Safety against Inf/NaN
		if (d=="ratio") hdm->SetBinContent(i, R1_d / R1_m);
		if (d=="mc")    hdm->SetBinContent(i, R1_m);
		// Statistical uncertainty will need more careful handling
		if (d=="ratio") hdm->SetBinError(i, sqrt(er0*er0 + eq0*eq0));
		if (d=="mc")    hdm->SetBinError(i, eq0);
	      }

	      if (d=="data")  hfsr->SetBinContent(i, dR_d);
	      if (dR>0 || dR<0) { // Safety against Inf/Nan
		if (d=="ratio") hfsr->SetBinContent(i, dR);
		if (d=="mc")    hfsr->SetBinContent(i, dR_m);
		hfsr->SetBinError(i, ddR);
		hfsr0s->SetBinContent(i, dR0s);
		hfsr1s->SetBinContent(i, dR1s);
		hfsrn1->SetBinContent(i, dRn1);
		hfsrn2->SetBinContent(i, dRn2);
		hfsrns->SetBinContent(i, dRns);
		hfsru1->SetBinContent(i, dRu1);
		hfsru2->SetBinContent(i, dRu2);
		hfsrus->SetBinContent(i, dRus);
	      }

	      cout << " MPF pt="<<pt<<" R1_d="<<R1_d<<" R1_m="<<R1_m
		   << " R1="<<R1_d/R1_m<<" raw="<<mpf<<endl;
	      cout << "  - r0="<<r0<<" r1="<<r1<<" rn="<<rn<<" ru="<<ru<<endl;
	      cout << "  - q0="<<q0<<" q1="<<q1<<" qn="<<qn<<" qu="<<qu<<endl;
	    } // mpfchs1 z

	    // For pT balance, correct estimate p1 based on pn and pu
	    if (m=="ptchs" && s!="multijet") {

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
	      double dR_d = R1_d - r1;
	      double dR_m = R1_m - q1;

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
	      //double R1m_d = fm->GetX(0,xmin,xmax);
	      //fm->SetParameters(q0, qn, qu, Rn_m, Ru_m);
	      //double R1m_m = fm->GetX(0,xmin,xmax);
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

	      if (d=="data")  hdm->SetBinContent(i, R1_d);
	      if (d=="data")  hdm->SetBinError(i, er0);
	      if (R1_m>0 || R1_m<0) { // Safety against Inf/NaN
		if (d=="ratio") hdm->SetBinContent(i, R1_d / R1_m);
		if (d=="mc")    hdm->SetBinContent(i, R1_m);
		// Statistical uncertainty will need more careful handling
		if (d=="ratio") hdm->SetBinError(i, sqrt(er0*er0 + eq0*eq0));
		if (d=="mc")    hdm->SetBinError(i, eq0);
	      }

	      if (d=="data")  hfsr->SetBinContent(i, dR_d);
	      if (dR_m>0 | dR_m<0) { // Safety against Inf/NaN
		if (d=="ratio") hfsr->SetBinContent(i, dR);
		if (d=="mc")    hfsr->SetBinContent(i, dR_m);
		hfsr->SetBinError(i, ddR);
		hfsr0s->SetBinContent(i, dR0s);
		hfsr1s->SetBinContent(i, dR1s);
		hfsrn1->SetBinContent(i, dRn1);
		hfsrn2->SetBinContent(i, dRn2);
		hfsrns->SetBinContent(i, dRns);
		hfsru1->SetBinContent(i, dRu1);
		hfsru2->SetBinContent(i, dRu2);
		hfsrus->SetBinContent(i, dRus);
	      }

	      cout << " pTbal pt="<<pt<<" R1_d="<<R1_d<<" R1_m="<<R1_m
		   << " R1="<<R1_d/R1_m<<" raw="<<ptbal<<endl;
	      cout << "  - r0="<<r0<<" r1="<<r1<<" rn="<<rn<<" ru="<<ru<<endl;
	      cout << "  - q0="<<q0<<" q1="<<q1<<" qn="<<qn<<" qu="<<qu<<endl;
	    } // ptchs z
	  } // for i

	  finout->mkdir(Form("%s/eta%02d-%02d/fsr",
			     cd,int(10*etamin),int(10*etamax)));
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
	  finout->cd(Form("%s/eta%02d-%02d",
			  cd,int(10*etamin),int(10*etamax)));
	  hdm->Write(Form("hdm_%s_%s",cm,cs),TObject::kOverwrite);
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
