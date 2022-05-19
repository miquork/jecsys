// Purpose: Draw <rho> vs <npv> (or <rho> vs mu) to extract UE
//          Analyze Z+jet, gamma+jet, dijet, ttbarj
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMultiGraph.h"
#include "TLine.h"
#include "TProfile2D.h"

#include "tdrstyle_mod15.C"

//string epoch = "UL16GH";
///string epoch = "UL16BCDEF";
string sdj = "jt450";//"jt400";//"jt40";//"jt400";
bool plotZB = true;
bool plotSJ = false;//true;
bool plotNG = true;//true;
bool plotDJ = true;
bool plotZmm = false;
bool plotZJ = true;//false;
bool plotGam = true;
bool plotTT = true;
bool plotMCup = false;//true;//false;
bool plotZBup = true;//false;//true;
//bool useMu = false;

// ZB MC seems to be weirdly off in UL16GH, while other MC agree well with data
bool normalizeMCbyDataZB = false;//true;

// Dijet profiles generated with 80mb setting? Fix it to 69.2mb
bool fixDijetMu = true;

// Reference mu for 80mb from
// https://twiki.cern.ch/twiki/bin/view/CMSPublic/LumiPublicResults#Multi_year_plots
// These settings can be changed per epoch
double muRef80 = 27;
double sigmaMB = 69.2;
double uefitmin = 8;
double uefitmax = 22;
double zbfitmin = 2;//8;
double zbfitmax = 25;//23;
double kmuzb = 1.00; // <mu>_ZB / <mu>_ref

void cleanGraph(TGraphErrors *g) {
  for (int i = g->GetN(); i != -1; --i) {
    if (g->GetEY()[i]==0 || (g->GetEY()[i]>1 || g->GetEX()[i]>1))
      g->RemovePoint(i);
  }
}

// map TProfiles
TGraphErrors* makeGraph(TProfile *p1, TProfile *p2,
			double dnpv=0, double due=0) {

  int nbinsx = max(p1->GetNbinsX(),p2->GetNbinsX());
  TGraphErrors *g = new TGraphErrors(0);
  for (int i = 1; i != nbinsx+1; ++i) {
    //if (p1->GetBinCenter(i)>4 && p1->GetBinCenter(i)<33) {
    if (true) {
      int n = g->GetN();
      g->SetPoint(n, p1->GetBinContent(i)+dnpv, p2->GetBinContent(i)+due);
      g->SetPointError(n, p1->GetBinError(i), p2->GetBinError(i));
    }
  }  

  return g;
}

// map TGraphErrors
TGraphErrors* makeGraph(TGraphErrors *g1, TGraphErrors *g2,
			double dnpv=0, double due=0) {

  TGraphErrors *g = new TGraphErrors(0);
  for (int i = 1; i != g1->GetN(); ++i) {
    for (int j = 1; j != g2->GetN(); ++j) {
    
      if (fabs(g1->GetX()[i]-g2->GetX()[j])<0.5) {
	double mu = 0.5*(g1->GetX()[i]+g2->GetX()[j]);
	double npv = g1->GetY()[i];
	double rho = g2->GetY()[j];
	//if (mu>4 && mu<33) {
	if (true) {
	  int n = g->GetN();
	  g->SetPoint(n, npv+dnpv, rho+due);
	  g->SetPointError(n, g1->GetEY()[i], g2->GetEY()[j]);
	}
      }
    }
  }

  return g;
}

void scaleGraph(TGraphErrors *g, TF1 *f) {
  for (int i = 0; i != g->GetN(); ++i) {
    g->SetPoint(i, g->GetX()[i], g->GetY()[i]/f->Eval(g->GetX()[i]));
    g->SetPointError(i, g->GetEX()[i], g->GetEY()[i]/f->Eval(g->GetX()[i]));
  }
}
void shiftGraph(TGraphErrors *g, TF1 *f) {
  for (int i = 0; i != g->GetN(); ++i) {
    g->SetPoint(i, g->GetX()[i], g->GetY()[i] - f->Eval(g->GetX()[i]));
    g->SetPointError(i, g->GetEX()[i], g->GetEY()[i]);
  }
}
void shiftGraph(TGraphErrors *g, double dy) {
  for (int i = 0; i != g->GetN(); ++i) {
    g->SetPoint(i, g->GetX()[i], g->GetY()[i] - dy);
    g->SetPointError(i, g->GetEX()[i], g->GetEY()[i]);
  }
}

// Fix profile mu from 80mb to 69.2mb
TGraphErrors *fixProf(TProfile *p) {
  TGraphErrors *g = new TGraphErrors(p);
  cleanGraph(g);
  for (int i = 0; i != g->GetN(); ++i) {
    g->SetPoint(i, g->GetX()[i]*69.2/80., g->GetY()[i]);
  }
  return g;
}

// Plot <Rho> vs <NPV>, both from bins of TruePU
void RhoVsNPV(string epoch = "UL16GH") {

  TDirectory *curdir = gDirectory;
  setTDRStyle();

  double rhoDJ = 2.20;//1.528;
  double rhoDJMC = rhoDJ+0.30;
  double rhoSJ = 1.70;//1.60;//1.;
  double rhoSJMC = rhoSJ+0.00;
  //double rhoSJNR = rhoSJ;
  double rhoZmm = 1.68;//1.64;//1.2;
  double rhoZmmMC = rhoZmm-0.12;//+0.19;
  double rhoZJ = 1.68;//1.64;//1.2;
  double rhoZJMC = rhoZJ-0.12;//+0.19;
  //double rhoZmmHW = rhoZmm;
  double rhoGJ = 1.8;//1.2;
  double rhoGJMC = rhoGJ+1.00;
  double rhoTT = 2.20;//1.2;
  double rhoTTMC = rhoTT+0.30;
  double rhoZB = 0.;
  double rhoZBMC = rhoZB;
  //double rhoZBNR = rhoZB;
  double kRho = 1.2; // Control slope in rhoUE
  double eDJ(0),eDJMC(0),eSJ(0),eSJMC(0),eZmm(0),eZmmMC(0),eGJ(0),eGJMC(0);
  double eZJ(0), eZJMC(0),eTT(0),eTTMC(0),eZB(0),eZBMC(0);

  //double vtxeff = 0.7050; // Run I
  //double vtxeff2 = -0.001334; // Run I
  double vtxeff = 0.764614; // Run II fG
  double vtxeff2 = -0.0030028; // Run II fG
  double rholin = 0.573921; // Run II fG
  double rhoquad = -0.000112284; // Run II fG
  double nfakeDJ = 0.3;

  map<string, int> mcol;
  map<string, int> mmcm;
  map<string, int> mdtm;
  mcol["sj"]=kBlue-9; mmcm["sj"]=kOpenDiamond; mdtm["sj"]=kFullDiamond;
  mcol["zmm"]=kRed; mmcm["zmm"]=kOpenCircle; mdtm["zmm"]=kFullCircle;
  mcol["zj"]=kRed+1; mmcm["zj"]=kOpenCircle; mdtm["zj"]=kFullCircle;
  mcol["gj"]=kBlue; mmcm["gj"]=kOpenDiamond; mdtm["gj"]=kFullDiamond;
  mcol["tt"]=kGreen+2; mmcm["tt"]=kOpenCross; mdtm["tt"]=kFullCross;
  mcol["zb"]=kGray+2; mmcm["zb"]=kFullDotSmall; mdtm["zb"]=kFullDotMedium;
  mcol["ng"] = kBlue-9;
  mmcm["ng"] = kOpenDiamond; mdtm["ng"] = kFullDiamond;
  mcol["dj"] = kBlack;
  mmcm["dj"] = kOpenSquare; mdtm["dj"] = kFullSquare;
  mcol["jt450"] = kBlack;
  mmcm["jt450"] = kOpenSquare; mdtm["jt450"] = kFullSquare;
  mcol["jt500"] = kBlack;
  mmcm["jt500"] = kOpenSquare; mdtm["jt500"] = kFullSquare;
  mcol["jt400"] = kBlack;
  mmcm["jt400"] = kOpenSquare; mdtm["jt400"] = kFullSquare;
  mcol["jt40"] = kGray+3;
  mmcm["jt40"] = kOpenSquare; mdtm["jt40"] = kFullSquare;
  mcol["jt0"] = kBlue-9;
  mmcm["jt0"] = kOpenDiamond; mdtm["jt0"] = kFullDiamond;

  // muRef from rootfiles/output-DATA-2a.root:Standard/Eta_0.0-1.3/jtX/htrpu
  if (epoch=="UL16BCDEF") { // V7V3_BCDEF, jt450
    lumi_13TeV="UL16APV, 19.7 fb^{-1}";
    muRef80 = 27.42;//24.15;
    uefitmin = 8; uefitmax = 19; // <NPV>~14
    zbfitmin = 2; zbfitmax = 22; // <NPV>~14
    sdj = "jt450";
  }
  else if (epoch=="UL16GH") { // V7V3_GH, jt450
    lumi_13TeV="UL16GH, 16.8 fb^{-1}";
    muRef80 = 28.14;
    uefitmin = 9; uefitmax = 21; // <NPV>~16
    zbfitmin = 2; zbfitmax = 25; // <NPV>~14
    sdj = "jt450";
  }
  else if (epoch=="UL17") { // V5V3_BCDEF, jt500
    lumi_13TeV = "UL17, 41.5 fb^{-1}";
    muRef80 = 37.38;
    uefitmin = 9; uefitmax = 21; // <NPV>~16?
    zbfitmin = 2; zbfitmax = 25; // <NPV>~14?
    sdj = "jt500";
  }
  else if (epoch=="UL18") { // V5V2_ABCD, jt500
    //assert(false);
    lumi_13TeV = "UL18, 59.9 fb^{-1}";
    muRef80 = 36.58; // <NPV> ~23
    uefitmin = 12; uefitmax = 29;
    zbfitmin = 2; zbfitmax = 40;
    sdj = "jt500";
  }
  else assert(false);
  extraText = "Preliminary";

  TLine *l = new TLine();
  l->SetLineStyle(kDotted);

  const int mumin = 10;
  //TH1D *h3 = new TH1D("h3",";#mu;#LTN_{PV,good}-N_{HS}#GT / #epsilon_{vtx}(#mu)",50-mumin,mumin,50);
  TH1D *h3 = new TH1D("h3",";#mu;#LTN_{PV,good}-1#GT / #epsilon_{vtx}(#mu)",50-mumin,mumin,50);
  h3->SetMaximum(1.30);//2.0);//3.00);
  h3->SetMinimum(0.75);//0.50);
  TCanvas *c3 = tdrCanvas("c3",h3,4,11,kSquare);
  l->DrawLine(10,1,50,1);
  double muRef = muRef80*sigmaMB/80.;
  l->DrawLine(muRef,0.9,muRef,1.1);

  //TF1 *fs = new TF1("fz","([0]+[1]*x)*x",0,50);
  //fs->SetParameters(vtxeff,vtxeff2);
  //TF1 *fsmc = fs;
  TF1 *fs = new TF1("fs","[0]+[1]*x+[2]*x*x",mumin,50);
  fs->SetParameters(0, 0.7396, -0.002563); // chi^2/NDF=59.1/38
  TF1 *fsmc = new TF1("fsmc","[0]+[1]*x+[2]*x*x",mumin,50);
  TF1 *fsmczj = new TF1("fsmczj","[0]+[1]*x+[2]*x*x",mumin,50);
  fsmc->SetParameters(0, 0.6973, -0.001342); // chi^2/NDF=1987.2/38
  if (epoch=="UL16BCDEF") {
    fs->SetParameters(-0.1523, 0.7469, -0.00362); // chi^2/NDF=57.5/35
    fsmc->SetParameters(-0.09664, 0.7308, -0.003001); // chi^2/NDF=2590.2/37
  }
  else if (epoch=="UL16GH") {
    fs->SetParameters(0, 0.7396, -0.002563); // chi^2/NDF=59.1/38
    fsmc->SetParameters(0, 0.6973, -0.001342); // chi^2/NDF=1987.2/38
  }
  else if (epoch=="UL17") {
    fs->SetParameters(-0.1122, 0.8017, -0.001248); // chi^2/NDF=1006.4/97
    fsmc->SetParameters(-0.3749, 0.7868, -0.001729); // chi^2/NDF=275.8/97
  }
  else if (epoch=="UL18") {
    fs->SetParameters(0.05813, 0.722, 0.0003573); // chi^2/NDF=644.0/92
    fsmc->SetParameters(-0.1223, 0.7585, -0.001076); // chi^2/NDF=193.3/92
    fsmczj->SetParameters(-0.1223, 1.85*0.7585, -0.001076); // chi^2/NDF=193.3/92
  }
  else assert(false);

  TH1D *h4 = new TH1D("h4",";#mu;#LT#rho-#rho_{UE}#GT / #rho_{PU}(#mu)",
		      50-mumin,mumin,50);
  h4->SetMaximum(1.30);//2.00);//3.00);
  h4->SetMinimum(0.75);//0.50);
  TCanvas *c4 = tdrCanvas("c4",h4,4,11,kSquare);
  l->DrawLine(10,1,50,1);
  l->DrawLine(muRef,0.9,muRef,1.1);

  //TF1 *fpu = new TF1("fpu","([0]+[1]*x)*x",0,50);
  //fpu->SetParameters(rholin,rhoquad);
  //TF1 *fpumc = fpu;
  TF1 *fpu = new TF1("fpu","[0]+[1]*x+[2]*x*x",mumin,50);
  fpu->SetParameters(-0.4931, 0.6138, -0.0005717); // chi^2/NDF=42.8/35
  TF1 *fpumc = new TF1("fpumc","[0]+[1]*x+[2]*x*x",mumin,50);
  TF1 *fpumczj = new TF1("fpumczj","[0]+[1]*x+[2]*x*x",mumin,50);
  fpumc->SetParameters(-0.8043, 0.6066, -0.0004275); // chi^2/NDF=611.2/35
  if (epoch=="UL16BCDEF") {
    fpu->SetParameters(-0.8563, 0.6232, -0.001411); // chi^2/NDF=52.1/34
    fpumc->SetParameters(-0.8514, 0.6195, -0.001314); // chi^2/NDF=999.8/35
  }
  else if (epoch=="UL16GH") {
    fpu->SetParameters(-0.4931, 0.6138, -0.0005717); // chi^2/NDF=42.8/35
    fpumc->SetParameters(-0.8043, 0.6066, -0.0004275); // chi^2/NDF=611.2/35
  }
  else if (epoch=="UL17") {
    fpu->SetParameters(-0.6534, 0.6344, -0.0005665); // chi^2/NDF=607.7/93
    fpumc->SetParameters(-1.123, 0.6474, -0.001016); // chi^2/NDF=234.3/93
  }
  else if (epoch=="UL18") {
    fpu->SetParameters(-0.7638, 0.6152, -0.0006248); // chi^2/NDF=1078.8/92
    fpumc->SetParameters(-0.6688, 0.6039, -0.0003828); // chi^2/NDF=213.0/92
    fpumczj->SetParameters(-0.6688, 1.85*0.6039, -0.0003828); // chi^2/NDF=213.0/92
  }
  else assert(false);

  map<string, TGraphErrors* > gds;
  map<string, TGraphErrors* > gms; // RD MC
  map<string, TGraphErrors* > gns; // non-RD MC
  if (true) { // Zmm+jet

    TGraphErrors *gmn(0), *gmr(0), *gdn(0), *gdr(0);
    string dir="";
    TFile *fz(0);
    if (epoch=="UL16GH") {
      //fz = new TFile("../JERCProtoLab/Summer19UL16/L3Residual_Z/nonAPV/JEC_Combination_Zmm/ZJetCombination_Zmm_DYJets_amc_21Feb2020_Summer19UL16_V5old_L1L2L3Res.root","READ");
      fz = new TFile("../JERCProtoLab/Summer19UL16/L3Residual_Z/nonAPV/JEC_Combination_Zmm/ZJetCombination_Zmm_DYJets_amc_21Feb2020_Summer19UL16_V7_L1L2L3Res.root","READ");
      //assert(fz && !fz->IsZombie());
      dir = "Run2016postVFPFlateGH";
      /*
      gmn = (TGraphErrors*)fz->Get("Run2016postVFPFlateGH/MC_NPV_vs_npumean_CHS_a100_eta_00_13_L1L2L3Res"); assert(gmn);
      gmr = (TGraphErrors*)fz->Get("Run2016postVFPFlateGH/MC_Rho_vs_npumean_CHS_a100_eta_00_13_L1L2L3Res"); assert(gmr);
      gdn = (TGraphErrors*)fz->Get("Run2016postVFPFlateGH/Data_NPV_vs_npumean_CHS_a100_eta_00_13_L1L2L3Res"); assert(gdn);
      gdr = (TGraphErrors*)fz->Get("Run2016postVFPFlateGH/Data_Rho_vs_npumean_CHS_a100_eta_00_13_L1L2L3Res"); assert(gdr);
      */
    }
    else if (epoch=="UL16BCDEF") {
      //fz = new TFile("../JERCProtoLab/Summer19UL16/L3Residual_Z/APV/JEC_Combination_Zmm/ZJetCombination_Zmm_DYJets_amc_21Feb2020_Summer19UL16_V5old_L1L2L3Res.root","READ");
      fz = new TFile("../JERCProtoLab/Summer19UL16/L3Residual_Z/APV/JEC_Combination_Zmm/ZJetCombination_Zmm_DYJets_amc_21Feb2020_Summer19UL16_V7_L1L2L3Res.root","READ");
      //assert(fz && !fz->IsZombie());
      dir = "Run2016preVFPBCDEFearly";
      /*
      gmn = (TGraphErrors*)fz->Get("Run2016preVFPBCDEFearly/MC_NPV_vs_npumean_CHS_a100_eta_00_13_L1L2L3Res"); assert(gmn);
      gmr = (TGraphErrors*)fz->Get("Run2016preVFPBCDEFearly/MC_Rho_vs_npumean_CHS_a100_eta_00_13_L1L2L3Res"); assert(gmr);
      gdn = (TGraphErrors*)fz->Get("Run2016preVFPBCDEFearly/Data_NPV_vs_npumean_CHS_a100_eta_00_13_L1L2L3Res"); assert(gdn);
      gdr = (TGraphErrors*)fz->Get("Run2016preVFPBCDEFearly/Data_Rho_vs_npumean_CHS_a100_eta_00_13_L1L2L3Res"); assert(gdr);
      */
    }
    else if (epoch=="UL17") {
      fz = new TFile("../JERCProtoLab/Summer19UL17/L3Residual_Z/JEC_Combination_Zmm/ZJetCombination_Zmm_DYJets_Madgraph_09Aug2019_Summer19UL17_JECV6_L1L2L3Res.root","READ");
      dir = "Run2017BCDEF";
    }
    else if (epoch=="UL18") {
      fz = new TFile("../JERCProtoLab/Summer19UL18/L3Residual_Z/JEC_Combination_Zmm/splitZPtBin70/ZJetCombination_Zmm_DYJets_Madgraph_12Nov2019_Summer19UL18_JECV4_L1L2L3Res.root","READ");
      dir = "Run2018ABCD";
    }
    else
      assert(false);

    assert(fz && !fz->IsZombie());

    const char *cd = dir.c_str(); 
    gmn = (TGraphErrors*)fz->Get(Form("%s/MC_NPV_vs_npumean_CHS_a100_eta_00_13_L1L2L3Res",cd)); assert(gmn);
    gmr = (TGraphErrors*)fz->Get(Form("%s/MC_Rho_vs_npumean_CHS_a100_eta_00_13_L1L2L3Res",cd)); assert(gmr);
    gdn = (TGraphErrors*)fz->Get(Form("%s/Data_NPV_vs_npumean_CHS_a100_eta_00_13_L1L2L3Res",cd)); assert(gdn);
    gdr = (TGraphErrors*)fz->Get(Form("%s/Data_Rho_vs_npumean_CHS_a100_eta_00_13_L1L2L3Res",cd)); assert(gdr);

    TGraphErrors *gm = makeGraph(gmn,gmr, -1, 0);
    TGraphErrors *gd = makeGraph(gdn,gdr, -1, 0);

    cleanGraph(gd);
    cleanGraph(gm);

    gds["zmm"] = gd;
    gms["zmm"] = gm;

    c3->cd();
    TGraphErrors *gd1 = (TGraphErrors*)gdn->Clone("gd1z"); 
    cleanGraph(gd1); shiftGraph(gd1,1.); scaleGraph(gd1,fs);
    if (plotZmm) tdrDraw(gd1,"Pz",kFullCircle,kRed);
    TGraphErrors *gm1 = (TGraphErrors*)gmn->Clone("gm1z"); 
    cleanGraph(gm1); shiftGraph(gm1,1.); scaleGraph(gm1,fsmc);
    if (plotZmm) tdrDraw(gm1,"Pz",kOpenCircle,kRed);

    c4->cd();
    TGraphErrors *gd2 = (TGraphErrors*)gdr->Clone("gd2z"); 
    cleanGraph(gd2); shiftGraph(gd2,rhoZmm); scaleGraph(gd2,fpu);
    if (plotZmm) tdrDraw(gd2,"Pz",kFullCircle,kRed);
    TGraphErrors *gm2 = (TGraphErrors*)gmr->Clone("gm2z"); 
    cleanGraph(gm2); shiftGraph(gm2,rhoZmmMC); scaleGraph(gm2,fpumc);
    if (plotZmm) tdrDraw(gm2,"Pz",kOpenCircle,kRed);
  } // Zmm+jet

  if (true) { // Z+jet

    TGraphErrors *gmn(0), *gmr(0), *gdn(0), *gdr(0);
    TFile *fz(0);
    if (epoch=="UL18") {
      fz = new TFile("rootfiles/jme_bplusZ_merged_v44_muon2018.root","READ");
    }
    else
      assert(false);

    assert(fz && !fz->IsZombie());
    curdir->cd();

    TProfile2D *p2mn = (TProfile2D*)fz->Get("mc/eta_00_13/h_Zpt_Mu_NPVgood_alpha100"); assert(p2mn); p2mn->SetName("p2mn");
    TProfile2D *p2mr = (TProfile2D*)fz->Get("mc/eta_00_13/h_Zpt_Mu_Rho_alpha100"); assert(p2mr); p2mr->SetName("p2mr");
    TProfile2D *p2dn = (TProfile2D*)fz->Get("data/eta_00_13/h_Zpt_Mu_NPVgood_alpha100"); assert(p2dn); p2dn->SetName("p2dn");
    TProfile2D *p2dr = (TProfile2D*)fz->Get("data/eta_00_13/h_Zpt_Mu_Rho_alpha100"); assert(p2dr); p2dr->SetName("p2dr");

    
    p2mn->Rebin2D(7,1);
    p2mr->Rebin2D(7,1);
    p2dn->Rebin2D(7,1);
    p2dr->Rebin2D(7,1);
    int i1 = p2mn->GetXaxis()->FindBin(90.);
    int i2 = p2mn->GetXaxis()->FindBin(270.);
    cout << "Z+jet : " << p2mn->GetXaxis()->GetBinLowEdge(i1) << " - "
	 << p2mn->GetXaxis()->GetBinLowEdge(i1+1) << " GeV" << endl;
    //gmn = new TGraphErrors(p2mn->ProfileY("p1mn",i1,i2));//->ProjectionX("h1mn"));
    //gmr = new TGraphErrors(p2mr->ProfileY("p1mr",i1,i2));//->ProjectionX("h1mr"));
    //gdn = new TGraphErrors(p2dn->ProfileY("p1dn",i1,i2));//->ProjectionX("h1dn"));
    //gdr = new TGraphErrors(p2dr->ProfileY("p1dr",i1,i2));//->ProjectionX("h1dr"));
    gmn = new TGraphErrors(p2mn->ProjectionY("p1mn",i1,i1));
    gmr = new TGraphErrors(p2mr->ProjectionY("p1mr",i1,i1));
    gdn = new TGraphErrors(p2dn->ProjectionY("p1dn",i1,i1));
    gdr = new TGraphErrors(p2dr->ProjectionY("p1dr",i1,i1));

    /*
    TProfile *pmn = (TProfile*)fz->Get("mc/NpvVsMu");
    TProfile *pmr = (TProfile*)fz->Get("mc/RhoVsMu");
    TProfile *pdn = (TProfile*)fz->Get("data/NpvVsMu");
    TProfile *pdr = (TProfile*)fz->Get("data/RhoVsMu");
    gmn = new TGraphErrors(pmn->ProjectionX("p1mn"));
    gmr = new TGraphErrors(pmr->ProjectionX("p1mr"));
    gdn = new TGraphErrors(pdn->ProjectionX("p1dn"));
    gdr = new TGraphErrors(pdr->ProjectionX("p1dr"));
    */

    // Patch Z+jet MC
    double k = 1.87;//1.85;
    for (int i = 0; i != gmn->GetN(); ++i) {
      gmn->SetPoint(i, gmn->GetX()[i], (gmn->GetY()[i]-1)/k+1);
    }
    for (int i = 0; i != gmr->GetN(); ++i) {
      gmr->SetPoint(i, gmn->GetX()[i], gmr->GetY()[i]/k);
    }

    TGraphErrors *gm = makeGraph(gmn,gmr, -1, 0);
    TGraphErrors *gd = makeGraph(gdn,gdr, -1, 0);

    cleanGraph(gd);
    cleanGraph(gm);

    gds["zj"] = gd;
    gms["zj"] = gm;

    c3->cd();
    TGraphErrors *gd1 = (TGraphErrors*)gdn->Clone("gd1zj"); 
    cleanGraph(gd1); shiftGraph(gd1,1.); scaleGraph(gd1,fs);
    if (plotZJ) tdrDraw(gd1,"Pz",kFullCircle,kRed+1);
    TGraphErrors *gm1 = (TGraphErrors*)gmn->Clone("gm1zj"); 
    cleanGraph(gm1); shiftGraph(gm1,1.); scaleGraph(gm1,fsmc);//zj);
    if (plotZJ) tdrDraw(gm1,"Pz",kOpenCircle,kRed+1);

    c4->cd();
    TGraphErrors *gd2 = (TGraphErrors*)gdr->Clone("gd2zj"); 
    cleanGraph(gd2); shiftGraph(gd2,rhoZJ); scaleGraph(gd2,fpu);
    if (plotZJ) tdrDraw(gd2,"Pz",kFullCircle,kRed+1);
    TGraphErrors *gm2 = (TGraphErrors*)gmr->Clone("gm2zj"); 
    cleanGraph(gm2); shiftGraph(gm2,rhoZJMC); scaleGraph(gm2,fpumc);//zj);
    if (plotZJ) tdrDraw(gm2,"Pz",kOpenCircle,kRed+1);
  } // Z+jet

  if (true) { // gamma+jet, missing for 2016 Legacy re-reco

    TFile *fd(0), *fm(0);
    if (epoch=="UL16GH") {
      //fd = new TFile("rootfiles/PhotonJet_DATA_L2L3Res_1_0_2021-05-18_FGH_L2L3Res_PFlowAK4chs.root","READ");
      //fm = new TFile("rootfiles/PhotonJet_MC_L2L3Res_1_0_2021-05-18_FGH_L2L3Res_PFlowAK4chs.root","READ");
      //fd = new TFile("rootfiles/PhotonJet_DATA_L2L3Res_0_3_2021-05-18_FGH_L2L3Res_PFlowAK4chs.root","READ");
      //fm = new TFile("rootfiles/PhotonJet_MC_L2L3Res_0_3_2021-05-18_FGH_L2L3Res_PFlowAK4chs.root","READ");
      //fd = new TFile("../gamjet/files/GamHistosFill_data_2016FGH_v10.root");
      //fm = new TFile("../gamjet/files/GamHistosFill_mc_2016P8_v10.root");
      //fd = new TFile("../gamjet/files/GamHistosFill_data_2016FGH_v17.root");
      //fm = new TFile("../gamjet/files/GamHistosFill_mc_2016P8_v17.root");
      fd = new TFile("../gamjet/files/GamHistosFill_data_2016FGH_v20.root");
      fm = new TFile("../gamjet/files/GamHistosFill_mc_2016P8_v20.root");
    }
    else if (epoch=="UL16BCDEF") {
      //fd = new TFile("rootfiles/PhotonJet_DATA_L2L3Res_1_0_2021-05-18_BCDEI_L2L3Res_PFlowAK4chs.root","READ");
      //fm = new TFile("rootfiles/PhotonJet_MC_L2L3Res_1_0_2021-05-18_BCDEI_L2L3Res_PFlowAK4chs.root","READ");
      //fd = new TFile("rootfiles/PhotonJet_DATA_L2L3Res_0_3_2021-05-18_BCDEI_L2L3Res_PFlowAK4chs.root","READ");
      //fm = new TFile("rootfiles/PhotonJet_MC_L2L3Res_0_3_2021-05-18_BCDEI_L2L3Res_PFlowAK4chs.root","READ");
      //fd = new TFile("../gamjet/files/GamHistosFill_data_2016BCDEF_v10.root");
      //fm = new TFile("../gamjet/files/GamHistosFill_mc_2016P8APV_v10.root");
      fd = new TFile("../gamjet/files/GamHistosFill_data_2016BCDEF_v20.root");
      fm = new TFile("../gamjet/files/GamHistosFill_mc_2016APVP8_v20.root");
    }
    else if (epoch=="UL17") {
      fd = new TFile("../gamjet/files/GamHistosFill_data_2017BCDEF_v20.root");
      fm = new TFile("../gamjet/files/GamHistosFill_mc_2017P8_v20.root");
      //fm = new TFile("../gamjet/files/GamHistosMix_mc_2017P8QCD_v20.root");
      // consider also adding 2018QCD, which has bit higher rho vs mu
    }
    else if (epoch=="UL18") {
      fd = new TFile("../gamjet/files/GamHistosFill_data_2018ABCD_v20.root");
      fm = new TFile("../gamjet/files/GamHistosFill_mc_2018P8_v20.root");
      //fm = new TFile("../gamjet/files/GamHistosMix_mc_2018P8QCD_v20.root");
      // consider also adding 2018QCD, which has bit higher rho vs mu
    }
    else
      assert(false);

    //TFile *f = new TFile("files/gjet_rho_vs_mu_info-nobadruns.root","READ");
    //TFile *f = new TFile("rootfiles/gjet_rho_vs_mu_2016BCDEFGH.root","READ");
    //TFile *fd = new TFile("rootfiles/gjet_rho_vs_mu_2016fG.root","READ");
    //TFile *fd = new TFile("rootfiles/PhotonJet_DATA_L2L3Res_1_0_2021-05-18_BCDEI_L2L3Res_PFlowAK4chs.root","READ");
    assert(fd && !fd->IsZombie());
    //TFile *fm = new TFile("rootfiles/gjet_rho_vs_mu_2016MC.root","READ");
    //TFile *fm = new TFile("rootfiles/PhotonJet_MC_L2L3Res_1_0_2021-05-18_BCDEI_L2L3Res_PFlowAK4chs.root","READ");
    assert(fm && !fm->IsZombie());
    //assert(f->cd("histos"));
    //TDirectory *fp = gDirectory;
    // ptgt180, ptgt70, allpt
    //TProfile *pp1 = (TProfile*)fp->Get("ptgt180_npv_vs_mu_MCgjetANDqcd");
    //TProfile *ppm1 = (TProfile*)fm->Get("analysis/MUDir/NpvGood_vs_mu_eta0013_ptPhot_230_300");
    TProfile *ppm1 = (TProfile*)fm->Get("control/pnpvgoodvsmu");
    assert(ppm1);
    //TProfile *pp2 = (TProfile*)fp->Get("ptgt180_rho_vs_mu_MCgjetANDqcd");
    //TProfile *ppm2 = (TProfile*)fm->Get("analysis/MUDir/rho_vs_mu_eta0013_ptPhot_230_300");
    TProfile *ppm2 = (TProfile*)fm->Get("control/prhovsmu");
    assert(ppm2);
    //TProfile *ppd1 = (TProfile*)fp->Get("ptgt180_npv_vs_mu_DATA");
    //TProfile *ppd1 = (TProfile*)fd->Get("analysis/MUDir/NpvGood_vs_mu_eta0013_ptPhot_230_300");
    TProfile *ppd1 = (TProfile*)fd->Get("control/pnpvgoodvsmu");
    //175_230");//"40_50");
    assert(ppd1);
    //TProfile *ppd2 = (TProfile*)fp->Get("ptgt180_rho_vs_mu_DATA");
    //TProfile *ppd2 = (TProfile*)fd->Get("analysis/MUDir/rho_vs_mu_eta0013_ptPhot_230_300");
    TProfile *ppd2 = (TProfile*)fd->Get("control/prhovsmu");
    //175_230");//40_50");
    assert(ppd2);

    ppm1->Scale(0.96); // for nAOD?

    TGraphErrors *gm = makeGraph(ppm1, ppm2, -1, 0);
    TGraphErrors *gd = makeGraph(ppd1, ppd2, -1, 0);

    cleanGraph(gd);
    cleanGraph(gm);

    gds["gj"] = gd;
    gms["gj"] = gm;
    
    if (plotGam) {
      c3->cd();
      TGraphErrors *gd1 = new TGraphErrors(ppd1);
      cleanGraph(gd1); shiftGraph(gd1,1.); scaleGraph(gd1,fs);
      tdrDraw(gd1,"Pz",kFullDiamond,kBlue);
      TGraphErrors *gm1 = new TGraphErrors(ppm1);
      cleanGraph(gm1); shiftGraph(gm1,1.); scaleGraph(gm1,fsmc);
      tdrDraw(gm1,"Pz",kOpenDiamond,kBlue);

      c4->cd();
      TGraphErrors *gd2 = new TGraphErrors(ppd2);
      cleanGraph(gd2); shiftGraph(gd2,rhoGJ); scaleGraph(gd2,fpu);
      tdrDraw(gd2,"Pz",kFullDiamond,kBlue);
      TGraphErrors *gm2 = new TGraphErrors(ppm2);
      cleanGraph(gm2); shiftGraph(gm2,rhoGJMC); scaleGraph(gm2,fpumc);
      tdrDraw(gm2,"Pz",kOpenDiamond,kBlue);
    }
  } // gamma+jet

  if (true) { // ttbar

    /*
    //TFile *f = new TFile("files/ttbarue_obspu.root","READ");
    TFile *f = new TFile("files/ttbarue_truepu.root","READ");
    assert(f && !f->IsZombie());

    TProfile *pd1 = (TProfile*)f->Get("Profile_nvx_mu_dilep_data");
    assert(pd1);
    TProfile *pd2 = (TProfile*)f->Get("Profile_rho_mu_dilep_data");
    assert(pd2);    

    TProfile *pm1 = (TProfile*)f->Get("Profile_nvx_mu_dilep_mc");
    assert(pm1);
    TProfile *pm2 = (TProfile*)f->Get("Profile_rho_mu_dilep_mc");
    assert(pm2);

    TProfile *pn1 = (TProfile*)f->Get("Profile_nvx_mu_dilep_mc");
    assert(pn1);
    TProfile *pn2 = (TProfile*)f->Get("Profile_rho_mu_dilep_mc");
    assert(pn2);
    */

    TFile *fd(0), *fm(0);
    if (epoch=="UL16GH") {
      fd = new TFile("rootfiles/hadWUL16GHV5_Glu.root","READ");
      fm = new TFile("rootfiles/hadWMC16GHV5_Glu.root","READ");
    }
    else if (epoch=="UL16BCDEF") {
      fd = new TFile("rootfiles/hadWUL16APVV5_Glu.root","READ");
      fm = new TFile("rootfiles/hadWMC16APVV5_Glu.root","READ");
    }
    else if (epoch=="UL17") {
      fd = new TFile("rootfiles/hadWUL17V5_JECv2.root","READ");
      fm = new TFile("rootfiles/hadWMC17V5_JECv2.root","READ");
    }
    else if (epoch=="UL18") {
      fd = new TFile("rootfiles/hadWUL18V5_JECv2.root","READ");
      fm = new TFile("rootfiles/hadWMC18V5_JECv2.root","READ");
    }
    else assert(false);
    assert(fd && !fd->IsZombie());
    assert(fm && !fm->IsZombie());

    TProfile *pd1 = (TProfile*)fd->Get("pnpvvsmu");
    assert(pd1);
    TProfile *pd2 = (TProfile*)fd->Get("prhovsmu");
    assert(pd2);    

    TProfile *pm1 = (TProfile*)fm->Get("pnpvvsmu");
    assert(pm1);
    TProfile *pm2 = (TProfile*)fm->Get("prhovsmu");
    assert(pm2);

    // Already vs <Npv-NHS> ?
    TGraphErrors *gd = makeGraph(pd1, pd2, -1, 0);
    TGraphErrors *gm = makeGraph(pm1, pm2, -1, 0);
    //TGraphErrors *gn = makeGraph(pn1, pn2, 0, 0);

    cleanGraph(gd);
    cleanGraph(gm);
    //cleanGraph(gn);

    gds["tt"] = gd;
    gms["tt"] = gm;
    //gns["tt"] = gn;

    if (plotTT) {
      
      c3->cd();
      TGraphErrors *gd1 = new TGraphErrors(pd1);
      cleanGraph(gd1); shiftGraph(gd1,1); scaleGraph(gd1,fs);
      tdrDraw(gd1,"Pz",kFullCross,kGreen+2);
      TGraphErrors *gm1 = new TGraphErrors(pm1);
      cleanGraph(gm1); shiftGraph(gm1,1); scaleGraph(gm1,fsmc);
      tdrDraw(gm1,"Pz",kOpenCross,kGreen+2);

      c4->cd();
      TGraphErrors *gd2 = new TGraphErrors(pd2);
      cleanGraph(gd2); shiftGraph(gd2,rhoTT); scaleGraph(gd2,fpu);
      tdrDraw(gd2,"Pz",kFullCross,kGreen+2);
      TGraphErrors *gm2 = new TGraphErrors(pm2);
      cleanGraph(gm2); shiftGraph(gm2,rhoTTMC); scaleGraph(gm2,fpumc);
      tdrDraw(gm2,"Pz",kOpenCross,kGreen+2);
    }
  } // ttbar


  const int ntrg = 11;//10;
  const char *trg[ntrg] = {"jt0","jt40","jt60","jt80","jt140","jt200","jt260",
			   "jt320","jt400","jt450","jt500"};
  if (true) { // dijet --- from Helsinki SMP-HAD files

    for (int itrg = 0; itrg != ntrg; ++itrg) {

      string s = trg[itrg];
      const char *ct = trg[itrg];
      TFile *fd(0), *fm(0);
      if (epoch=="UL16GH") {
	fd = new TFile("rootfiles/output-DATA-2a-UL16V5V2_GH.root","READ");
	// To be updated to MCNU-2a-UL16V5V2
	fm = new TFile("rootfiles/output-MCNU-2a-UL16V2V1_GH.root","READ");
      }
      else if (epoch=="UL16BCDEF") {
	fd = new TFile("rootfiles/output-DATA-2a-UL16V5V2_BCDEF.root","READ");
	// To be updated to MCNU-2a-UL16V5V2x
	fm = new TFile("rootfiles/output-MCNU-2a-UL16V3V1_BCDEF.root","READ");
      }
      else if (epoch=="UL17") {
	fd = new TFile("rootfiles/output-DATA-2a-UL17V5V3_BCDEF.root","READ");
	fm = new TFile("rootfiles/output-MC-2a-UL17V5V3_BCDEF.root","READ");
      }
      else if (epoch=="UL18") {
	fd = new TFile("rootfiles/output-DATA-2a-UL18V5V2_ABCD.root","READ");
	//fm = new TFile("rootfiles/output-MC-2a-UL18V5V2_ABCD.root","READ");
	//fd = new TFile("rootfiles/output-DATA-2a-UL18V5V2_ABCD-pThat.root","READ");
	fm = new TFile("rootfiles/output-MC-2a-UL18V5V2_ABCD-pThat.root","READ");
      }
      else assert(false);
      assert(fd && !fd->IsZombie());
      assert(fm && !fm->IsZombie());

      if ((epoch=="UL16BCDEF"||epoch=="UL16GH") && s=="jt500") continue;


      TProfile *pd1 = (TProfile*)fd->Get(Form("Standard/Eta_0.0-1.3/%s"
					      "/pnpvvstrpu",ct));
      assert(pd1);
      TProfile *pd2 = (TProfile*)fd->Get(Form("Standard/Eta_0.0-1.3/%s"
					      "/prhovstrpu",ct));
      assert(pd2);    
      
      TProfile *pm1 = (TProfile*)fm->Get(Form("Standard/Eta_0.0-1.3/%s"
					      "/pnpvvstrpu",ct));
      assert(pm1);
      TProfile *pm2 = (TProfile*)fm->Get(Form("Standard/Eta_0.0-1.3/%s"
					      "/prhovstrpu",ct));
      assert(pm2);

      //TFile *fh = new TFile("files/output-HW-1.root","READ");
      //assert(fh && !fh->IsZombie());
      //TProfile *ph1 = (TProfile*)fh->Get(Form("Standard/Eta_0.0-1.3/%s"
      //				      "/pnpvvstrpu",ct));
      //assert(ph1);
      //TProfile *ph2 = (TProfile*)fh->Get(Form("Standard/Eta_0.0-1.3/%s"
      //				      "/prhovstrpu",ct));
      //assert(ph2);
      
      TGraphErrors *gd = makeGraph(pd1, pd2, -1, 0);
      TGraphErrors *gm = makeGraph(pm1, pm2, -1, 0);
      //TGraphErrors *gh = makeGraph(ph1, ph2, -1, 0);
      
      if (fixDijetMu) {
	gd = makeGraph(fixProf(pd1), fixProf(pd2), -1, 0);
	//gm = makeGraph(fixProf(pm1), fixProf(pm2), -1, 0);
      }


      cleanGraph(gd);
      cleanGraph(gm);
      //cleanGraph(gh);

      gds[Form("dj_%s",ct)] = gd;
      gms[Form("dj_%s",ct)] = gm;
      //gns[Form("dj_%s",ct)] = gh;
    
      if ((plotNG && s=="jt0") || (plotDJ && s==sdj)) {

	c3->cd();
	TGraphErrors *gd1 = new TGraphErrors(pd1);
	if (fixDijetMu) gd1 = fixProf(pd1);
	cleanGraph(gd1); shiftGraph(gd1,1); scaleGraph(gd1,fs);
	//if (plotDJ) tdrDraw(gd1,"P",kFullStar,kBlue-9);
	tdrDraw(gd1,"Pz",mdtm[s],mcol[s]);
	TGraphErrors *gm1 = new TGraphErrors(pm1);
	//if (fixDijetMu) gm1 = fixProf(pm1);
	cleanGraph(gm1); shiftGraph(gm1,1); scaleGraph(gm1,fsmc);
	tdrDraw(gm1,"Pz",mmcm[s],mcol[s]);

	c4->cd();
	TGraphErrors *gd2 = new TGraphErrors(pd2);
	if (fixDijetMu) gd2 = fixProf(pd2);
	cleanGraph(gd2); shiftGraph(gd2,rhoDJ); scaleGraph(gd2,fpu);
	tdrDraw(gd2,"Pz",mdtm[s],mcol[s]);
	TGraphErrors *gm2 = new TGraphErrors(pm2);
	//if (fixDijetMu) gm2 = fixProf(pm2);
	cleanGraph(gm2); shiftGraph(gm2,rhoDJMC); scaleGraph(gm2,fpumc);
	tdrDraw(gm2,"Pz",mmcm[s],mcol[s]);
      }
    } // for itrg
  
    gds["dj"] = gds[Form("dj_%s",sdj.c_str())]; assert(gds["dj"]);
    gms["dj"] = gms[Form("dj_%s",sdj.c_str())]; assert(gms["dj"]);
    gds["ng"] = gds["dj_jt0"]; assert(gds["ng"]);
    gms["ng"] = gms["dj_jt0"]; assert(gms["ng"]);
  } // dijet

  if (true) { // ZB soft jet (sj) --- (later) from Buffalo files

    //TFile *fd = new TFile("rootfiles/zerobiasL1DataGH.root","READ");
    //TFile *fd = new TFile("rootfiles/rhoVsMu2017DE.root","READ");
    //TFile *fd = new TFile("rootfiles/output-DATA-1-UL17V4_E.root","READ");
    //TFile *fd = new TFile("rootfiles/output-DATA-2a-UL16V2V1_GH.root","READ");
    TFile *fd(0), *fm(0);
    if (epoch=="UL16GH") {
      fd = new TFile("rootfiles/output-DATA-2a-UL16V7V3_GH.root","READ");
      //fm = new TFile("rootfiles/output-DATA-2a-UL16V7V3_GH.root","READ");
      fm = new TFile("rootfiles/output-MC-2a-UL16V2V1_GH_oldv1.root","READ");
    }
    else if (epoch=="UL16BCDEF") {
      fd = new TFile("rootfiles/output-DATA-2a-UL16V7V3_BCDEF.root","READ");
      fm = new TFile("rootfiles/output-MC-2a-UL16V3V1_BCDEF_oldv1.root","READ");
    }
    else if (epoch=="UL17") {
      fd = new TFile("rootfiles/output-DATA-2a-UL17V5V3_BCDEF.root","READ");
      fm = new TFile("rootfiles/output-MC-2a-UL17V5V3_BCDEF.root","READ");
    }
    else if (epoch=="UL18") {
      fd = new TFile("rootfiles/output-DATA-2a-UL18V5V2_ABCD.root","READ");
      fm = new TFile("rootfiles/output-MC-2a-UL18V5V2_ABCD.root","READ");
    }
    else assert(false);
    assert(fd && !fd->IsZombie());
        //TProfile *pd1 = (TProfile*)fd->Get("pnpv30");
    TProfile *pd1 = (TProfile*)fd->Get("Standard/Eta_0.0-1.3/jt0/pnpvvstrpu");
    assert(pd1);
    //TProfile *pd2 = (TProfile*)fd->Get("prho30");
    TProfile *pd2 = (TProfile*)fd->Get("Standard/Eta_0.0-1.3/jt0/prhovstrpu");
    assert(pd2);    
      
    //TFile *fm = new TFile("rootfiles/zerobiasL1MCMC.root","READ");
    //TFile *fm = new TFile("rootfiles/output-MC-1-UL17V4_E.root","READ");
    //TFile *fm = new TFile("rootfiles/output-MC-2a-UL16V2V1_GH.root","READ");
    //TFile *fm = new TFile("rootfiles/output-MCNU-2a-UL16V2V1_GH.root","READ");
    assert(fm && !fm->IsZombie());
    //TProfile *pm1 = (TProfile*)fm->Get("pnpv30");
    TProfile *pm1 = (TProfile*)fm->Get("Standard/Eta_0.0-1.3/jt40/pnpvvstrpu");
    assert(pm1);
    //TProfile *pm2 = (TProfile*)fm->Get("prho30");
    TProfile *pm2 = (TProfile*)fm->Get("Standard/Eta_0.0-1.3/jt40/prhovstrpu");
    assert(pm2);

    TGraphErrors *gd = makeGraph(pd1, pd2, -1, 0);
    TGraphErrors *gm = makeGraph(pm1, pm2, -1, 0);
      
    cleanGraph(gd);
    cleanGraph(gm);

    gds["sj"] = gd;
    gms["sj"] = gm;

    if (plotSJ) {
      c3->cd();
      TGraphErrors *gd1 = new TGraphErrors(pd1);
      cleanGraph(gd1); shiftGraph(gd1,1.); scaleGraph(gd1,fs);
      tdrDraw(gd1,"Pz",kFullDiamond,kBlue-9); gd1->SetMarkerSize(0.7);
      TGraphErrors *gm1 = new TGraphErrors(pm1);
      cleanGraph(gm1); shiftGraph(gm1,1.); scaleGraph(gm1,fsmc);
      tdrDraw(gm1,"Pz",kOpenDiamond,kBlue-9); gm1->SetMarkerSize(0.7);
      
      c4->cd();
      TGraphErrors *gd2 = new TGraphErrors(pd2);
      cleanGraph(gd2); shiftGraph(gd2,rhoSJ); scaleGraph(gd2,fpu);
      tdrDraw(gd2,"Pz",kFullDiamond,kBlue-9); gd2->SetMarkerSize(0.7);
      TGraphErrors *gm2 = new TGraphErrors(pm2);
      cleanGraph(gm2); shiftGraph(gm2,rhoSJMC); scaleGraph(gm2,fpumc);
      tdrDraw(gm2,"Pz",kOpenDiamond,kBlue-9); gm2->SetMarkerSize(0.7);
    }
  } // ZB soft jets

  if (true) { // ZeroBias (zb)

    //const char *p = "../jecsys/rootfiles/Legacy_zerobias2016_Aug072017_"
    //"SingleNeutrinoRunIISumer16/";
    //const char *p = "rootfiles/L1OffsetUL16_VFP/";
    //TFile *fd = new TFile(Form("%sLegacy_FG_R4.root",p),"READ");
    //TFile *fd = new TFile(Form("%sLegacy_H_R4.root",p),"READ");
    //TFile *fd = new TFile(Form("%sLegacy_FGH_R4.root",p),"READ");
    const char *p = "rootfiles/";
    // File produced with minitools/drawZBrhoVsNpv.C
    //TFile *fd = new TFile(Form("%sL1RC_RhoVsNpv_S19.root",p),"READ");

    TFile *fd(0), *fm(0);
    if (epoch=="UL16GH" || epoch=="UL16BCDEF") {
      fd = new TFile(Form("%sL1RC_RhoVsNpv_S20.root",p),"READ");
      fm = fd;
    }
    else if (epoch=="UL17") {
      fd = new TFile("rootfiles/GarvitaZB_Data_UL2017BCDEF_R4.root","READ");
      fm = new TFile("rootfiles/GarvitaZBv2_MC_UL2017BCDEF_R4.root","READ");
    }
    else if (epoch=="UL18") {
      fd = new TFile("rootfiles/GarvitaZB_Data_UL2018ABCD_R4.root","READ");
      fm = new TFile("rootfiles/GarvitaZB_MC_UL2018ABCD_R4.root","READ");
    }
    assert(fd && !fd->IsZombie());
    assert(fm && !fm->IsZombie());

    TProfile *pd1(0), *pd2(0), *pd3(0), *pm1(0), *pm2(0), *pm3(0);
    if (epoch=="UL16GH") {
      pd1 = (TProfile*)fd->Get("p_nPV_nPU_UL16GH");
      pd2 = (TProfile*)fd->Get("p_rho_nPU_UL16GH");
      pd3 = (TProfile*)fd->Get("p_nPU_nPU_UL16GH");
      pm1 = (TProfile*)fm->Get("p_nPV_nPU_MC16GH");
      pm2 = (TProfile*)fm->Get("p_rho_nPU_MC16GH");
      pm3 = (TProfile*)fm->Get("p_nPU_nPU_MC16GH");
    }
    else if (epoch=="UL16BCDEF") {
      pd1 = (TProfile*)fd->Get("p_nPV_nPU_UL16BCDEF");
      pd2 = (TProfile*)fd->Get("p_rho_nPU_UL16BCDEF");
      pd3 = (TProfile*)fd->Get("p_nPU_nPU_UL16BCDEF");
      pm1 = (TProfile*)fm->Get("p_nPV_nPU_MC16BCDEF");
      pm2 = (TProfile*)fm->Get("p_rho_nPU_MC16BCDEF");
      pm3 = (TProfile*)fm->Get("p_nPU_nPU_MC16BCDEF");
    }
    else if (epoch=="UL17" || epoch=="UL18") {
      pd1 = (TProfile*)fd->Get("p_nPV_nPU"); pd1->SetName("p_nPV_nPU_DATA");
      pd2 = (TProfile*)fd->Get("p_rho_nPU"); pd2->SetName("p_rho_nPU_DATA");
      //pd3 = (TProfile*)fd->Get("p_nPU_nPU");
      pm1 = (TProfile*)fm->Get("p_nPV_nPU"); pm1->SetName("p_nPV_nPU_MC");
      pm2 = (TProfile*)fm->Get("p_rho_nPU"); pm2->SetName("p_rho_nPU_MC");
      //pm3 = (TProfile*)fm->Get("p_nPU_nPU");

      // Patch missing profiles
      TH1D *hd3 = pd1->ProjectionX("p_nPU_nPU_DATA");
      TH1D *hm3 = pm1->ProjectionX("p_nPU_nPU_MC");
      for (int i = 1; i != hd3->GetNbinsX()+1; ++i) {
	hd3->SetBinContent(i, hd3->GetBinCenter(i) * kmuzb);
      }
      for (int i = 1; i != hm3->GetNbinsX()+1; ++i) {
	hm3->SetBinContent(i, hm3->GetBinCenter(i));
      }
      pd3 = (TProfile*)hd3;
      pm3 = (TProfile*)hm3;
    }
    else assert(false);

    assert(pd1);
    assert(pd2);    
    assert(pd3);    

    //TFile *fm = new TFile(Form("%sSingleNeutrino_MC_R4.root",p),"READ");
    //TFile *fm = new TFile(Form("%sL1RC_RhoVsNpv_S19.root",p),"READ");
    assert(pm1);
    assert(pm2);
    assert(pm3);

    TGraphErrors *gd = makeGraph(pd1, pd2, +0, 0);
    TGraphErrors *gm = makeGraph(pm1, pm2, +0, 0);

    cleanGraph(gd);
    cleanGraph(gm);

    gds["zb"] = gd;
    gms["zb"] = gm;

    /*
    TF1 *fsx = new TF1("fsx","[0]+([1]+[2]*x)*x",5,30);
    fsx->FixParameter(0,0);
    pd1->Fit(fsx,"QRN");
    cout << "fsx: " << fsx->GetParameter(0) << " + "
	 << fsx->GetParameter(1) << "*x + "
	 << fsx->GetParameter(2) << "*x*x" << endl;
    cout << "fsx_err: " << fsx->GetParError(0) << " + "
	 << fsx->GetParError(1) << "*x + "
	 << fsx->GetParError(2) << "*x*x" << endl;
    */

    if (plotZB) {
      c3->cd();
      TGraphErrors *gd1 = new TGraphErrors(pd1);
      cleanGraph(gd1); scaleGraph(gd1,fs); 
      tdrDraw(gd1,"Pz",kFullDotMedium,kGray+2);
      TGraphErrors *gm1 = new TGraphErrors(pm1);
      cleanGraph(gm1); scaleGraph(gm1,fsmc);
      tdrDraw(gm1,"Pz",kFullDotSmall,kGray+2);

      c4->cd();
      TGraphErrors *gd2 = new TGraphErrors(pd2);
      cleanGraph(gd2); scaleGraph(gd2,fpu);
      tdrDraw(gd2,"Pz",kFullDotMedium,kGray+2);
      TGraphErrors *gm2 = new TGraphErrors(pm2);
      cleanGraph(gm2); scaleGraph(gm2,fpumc);
      tdrDraw(gm2,"Pz",kFullDotSmall,kGray+2);
    }
    if (true) { // do <Vtx> vs mu and <rho> vs mu for Zero Bias
      //TH1D *h0 = tdrHist("h0","#LT#rho#GT or #LTN_{PV}#GT",0,35,"#mu",0,50);
      TH1D *h0 = tdrHist("h0","#LT#rho#GT or #LTN_{PV}#GT",0,45,
			 "#LT#mu#GT or #LTN_{PV}#GT",0,56);
      TCanvas *c0 = tdrCanvas("c0",h0,4,11,kSquare);
      l->SetLineColor(kBlue+2);
      //l->DrawLine(0,0,50,35.0); // 0.70 vertex efficiency per interaction
      l->DrawLine(0,0,56,56*0.70); // 0.70 vertex efficiency per interaction
      l->SetLineColor(kRed+2);
      //l->DrawLine(0,0,50,50.*0.56); // 0.59 GeV per interaction
      l->DrawLine(0,0,56,56*0.57); // 0.57 GeV per interaction
      l->SetLineColor(kBlack);

      TGraphErrors *gd1 = makeGraph(pd3, pd1, +0, 0);
      TGraphErrors *gm1 = makeGraph(pm3, pm1, +0, 0);
      TGraphErrors *gd2 = makeGraph(pd3, pd2, +0, 0);
      TGraphErrors *gm2 = makeGraph(pm3, pm2, +0, 0);
      //
      TGraphErrors *gd3 = makeGraph(pd1, pd2, +0, 0);
      TGraphErrors *gm3 = makeGraph(pm1, pm2, +0, 0);

      tdrDraw(gd3,"Pz",kFullDotMedium,kGray+2);  //gd3->SetMarkerSize(0.7);
      tdrDraw(gm3,"Pz",kFullDotSmall,kGray+1);  //gm3->SetMarkerSize(0.7);
      tdrDraw(gd1,"Pz",kFullSquare,kBlue);   gd1->SetMarkerSize(0.5);
      tdrDraw(gm1,"Pz",kOpenSquare,kBlue-9); gm1->SetMarkerSize(0.5);
      tdrDraw(gd2,"Pz",kFullCircle,kRed);    gd2->SetMarkerSize(0.5);
      tdrDraw(gm2,"Pz",kOpenCircle,kRed-9);  gm2->SetMarkerSize(0.5);

      //cleanGraph(gd1);
      //cleanGraph(gm2);
      //cleanGraph(gd2);
      //cleanGraph(gm2);

      //TLegend *leg0 = tdrLeg(0.20,0.50,0.50,0.74);
      TLegend *leg0 = tdrLeg(0.20,0.44,0.50,0.80);
      leg0->AddEntry(gd3,"#LT#rho#GT vs #LTN_{PV}#GT data","PLE");
      leg0->AddEntry(gm3,"#LT#rho#GT  vs #LTN_{PV}#GT MC","PLE");
      leg0->AddEntry(gd1,"#LTN_{PV}#GT data","PLE");
      leg0->AddEntry(gm1,"#LTN_{PV}#GT MC","PLE");
      leg0->AddEntry(gd2,"#LT#rho#GT data","PLE");
      leg0->AddEntry(gm2,"#LT#rho#GT MC","PLE");

      if (epoch=="UL16BCDEF" || epoch=="UL16GH")
	l->DrawLine(muRef,5,muRef,20);
      else if (epoch=="UL17" || epoch=="UL18")
	l->DrawLine(muRef,10,muRef,30);

      TF1 *fs = new TF1("fs","[0]+[1]*x+[2]*x*x",0,40);
      if (epoch=="UL16BCDEF" || epoch=="UL16GH") fs->SetRange(0,40);
      else if (epoch=="UL17" || epoch=="UL18")  fs->SetRange(0,50);
      else assert(false);

      fs->SetLineColor(kBlue);
      gd1->Fit(fs,"QRN");
      bool fixP0 = (fabs(fs->GetParameter(0))-2*fs->GetParError(0)<0);
      if (fixP0) {
	fs->FixParameter(0,0);
	gd1->Fit(fs,"QRN");
      }
      fs->Draw("SAME");
      cout << Form("fs: %1.4g + %1.4g*x + %1.4g*x*x // chi^2/NDF=%1.1f/%d",
		   fs->GetParameter(0),fs->GetParameter(1),
		   fs->GetParameter(2),
		   fs->GetChisquare(),fs->GetNDF()) << endl;
      cout << Form("fs->SetParameters(%1.4g, %1.4g, %1.4g);"
		   " // chi^2/NDF=%1.1f/%d",
		   fs->GetParameter(0),fs->GetParameter(1),
		   fs->GetParameter(2),
		   fs->GetChisquare(),fs->GetNDF()) << endl;

      TF1 *fsmc = new TF1("fsmc","[0]+[1]*x+[2]*x*x",0,40);
      if (epoch=="UL16BCDEF" || epoch=="UL16GH") fsmc->SetRange(0,40);
      else if (epoch=="UL17" || epoch=="UL18")  fsmc->SetRange(0,50);
      else assert(false);
      fsmc->SetLineColor(kBlue-9); //fsmc->SetLineStyle(kDotted);
      gm1->Fit(fsmc,"QRN");
      if (fixP0) {
	fsmc->FixParameter(0,0);
	gm1->Fit(fsmc,"QRN");
      }
      fsmc->Draw("SAME");
      cout << Form("fsmc: %1.4g + %1.4g*x + %1.4g*x*x // chi^2/NDF=%1.1f/%d",
		   fsmc->GetParameter(0),fsmc->GetParameter(1),
		   fsmc->GetParameter(2),
		   fsmc->GetChisquare(),fsmc->GetNDF()) << endl;
      cout << Form("fsmc->SetParameters(%1.4g, %1.4g, %1.4g);"
		   " // chi^2/NDF=%1.1f/%d",
		   fsmc->GetParameter(0),fsmc->GetParameter(1),
		   fsmc->GetParameter(2),
		   fsmc->GetChisquare(),fsmc->GetNDF()) << endl;

      //TF1 *fpu = new TF1("fpu","[0]+[1]*x+[2]*x*x+[3]*x*x*x",2,40);
      TF1 *fpu = new TF1("fpu","[0]+[1]*x+[2]*x*x",2,40);
      if (epoch=="UL16BCDEF" || epoch=="UL16GH") fpu->SetRange(2,40);
      else if (epoch=="UL17" || epoch=="UL18")  fpu->SetRange(2,50);
      //else if (epoch=="UL18")  fpu->SetRange(2,50);
      else assert(false);
      fpu->SetLineColor(kRed);
      gd2->Fit(fpu,"QRN");
      fpu->Draw("SAME");
      cout << Form("fpu: %1.4g + %1.4g*x + %1.4g*x*x // chi^2/NDF=%1.1f/%d",
		   fpu->GetParameter(0),fpu->GetParameter(1),
		   fpu->GetParameter(2),
		   fpu->GetChisquare(),fpu->GetNDF()) << endl;
      cout << Form("fpu->SetParameters(%1.4g, %1.4g, %1.4g);"
		   " // chi^2/NDF=%1.1f/%d",
		   fpu->GetParameter(0),fpu->GetParameter(1),
		   fpu->GetParameter(2),
		   fpu->GetChisquare(),fpu->GetNDF()) << endl;

      //TF1 *fpumc = new TF1("fpumc","[0]+[1]*x+[2]*x*x+[3]*x*x*x",2,40);
      TF1 *fpumc = new TF1("fpumc","[0]+[1]*x+[2]*x*x",2,40);
      if (epoch=="UL16BCDEF" || epoch=="UL16GH") fpumc->SetRange(2,40);
      else if (epoch=="UL17" || epoch=="UL18")  fpumc->SetRange(2,50);
      //else if (epoch=="UL18")  fpumc->SetRange(2,50);
      fpumc->SetLineColor(kRed-9); //fpumc->SetLineStyle(kDotted);
      gm2->Fit(fpumc,"QRN");
      fpumc->Draw("SAME");
      cout << Form("fpumc: %1.4g + %1.4g*x + %1.4g*x*x // chi^2/NDF=%1.1f/%d",
		   fpumc->GetParameter(0),fpumc->GetParameter(1),
		   fpumc->GetParameter(2),
		   fpumc->GetChisquare(),fpumc->GetNDF()) << endl;
      cout << Form("fpumc->SetParameters(%1.4g, %1.4g, %1.4g);"
		   " // chi^2/NDF=%1.1f/%d",
		   fpumc->GetParameter(0),fpumc->GetParameter(1),
		   fpumc->GetParameter(2),
		   fpumc->GetChisquare(),fpumc->GetNDF()) << endl;


      c0->SaveAs(Form("pdf/RhoVsNPV/RhoVsNPV_ZBvsMu_%s.pdf",epoch.c_str()));
    }
  } // ZB

  curdir->cd();

  //TH1D *h2 = new TH1D("h2",";#LTN_{PV,good} - N_{HS}#GT;"
  TH1D *h2 = new TH1D("h2",";#LTN_{PV,good} - 1#GT;"
		      "#LT#rho#GT - ZB (GeV)",50,0,50);
  h2->SetMaximum(6.0);//4.5);
  h2->SetMinimum(0.0);//0.5);
  
  TH1D *h = new TH1D("h",";#LTN_{PV,good}#GT;#LT#rho#GT (GeV)",50,0,50);
  h->SetMaximum(47.); // 30->50
  h->SetMinimum(0.);

  TCanvas *c1 = tdrDiCanvas("c1",h,h2,4,0);
  h->GetYaxis()->SetTitleOffset(1.10);
  h2->GetXaxis()->SetTitleOffset(0.85);

  c1->cd(1);
  gPad->Update();
  
  if (normalizeMCbyDataZB) {
    TLatex *tex = new TLatex();
    tex->SetNDC(); tex->SetTextSize(0.045); tex->SetTextColor(kRed);
    tex->DrawLatex(0.35,0.87,"NB! MC normalized by ZB *DATA*");
  }

  // Original range 0,30. UL16GH range 2,25, UE 8,22
  //TF1 *f1mzb = new TF1("f1mzb","[0] + [1]*x + [2]*x*x",zbfitmin,zbfitmax);
  //TF1 *f1dzb = new TF1("f1dzb","[0] + [1]*x + [2]*x*x",zbfitmin,zbfitmax);
  TF1 *f1mzb = new TF1("f1mzb","[0]+[1]*x+[2]*x*x+[3]*x*x*x",zbfitmin,zbfitmax);
  TF1 *f1dzb = new TF1("f1dzb","[0]+[1]*x+[2]*x*x+[3]*x*x*x",zbfitmin,zbfitmax);
  //TF1 *f3 = new TF1("f2","[0] + [1]*x + [2]*x*x",0,30);
  TF1 *fue = new TF1("fue","[0]+[1]*x",uefitmin,uefitmax);
  TF1 *fuemc = new TF1("fuemc","[0]+[1]*x",uefitmin,uefitmax);
  fuemc->SetLineStyle(kSolid);//kDashed);

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.030);

  vector<string> v;
  if (plotDJ) v.push_back("dj");
  if (plotZmm) v.push_back("zmm");
  if (plotZJ) v.push_back("zj");
  if (plotGam) v.push_back("gj");
  if (plotTT) v.push_back("tt");
  if (plotSJ) v.push_back("sj");
  if (plotNG) v.push_back("ng");
  if (plotZB) v.push_back("zb");

  vector<int> colors(v.size());
  vector<int> mcmarker(v.size());
  vector<int> dtmarker(v.size());
  for (int i = 0; i != v.size(); ++i) {
    colors[i] = mcol[v[i]];
    mcmarker[i] = mmcm[v[i]];
    dtmarker[i] = mdtm[v[i]];
  }

  gds["zb"]->Fit(f1dzb,"QRN");
  f1dzb->SetLineWidth(4);
  f1dzb->SetLineColor(kGray);
  f1dzb->Draw("SAME");

  gms["zb"]->Fit(f1mzb,"QRN");
  f1mzb->SetLineWidth(2);//4);
  f1mzb->SetLineColor(kGray+1);
  if (plotZBup) f1mzb->Draw("SAME");

  for (unsigned int i = 0; i != v.size(); ++i) {

    string s = v[i];
    const char *c = v[i].c_str();
    bool plotX = ((s=="zb" && plotZB) || (s=="sj" && plotSJ) ||
		  (s=="ng" && plotNG) || (s=="dj" && plotDJ) ||
		  (s=="zmm" && plotZmm) || (s=="gj" && plotGam) ||
		  (s=="zj" && plotZJ) || (s=="tt" && plotTT));

    TGraphErrors *gd = gds[c];
    if (plotX) tdrDraw(gd,"Pz",dtmarker[i],colors[i],kSolid,colors[i]);
    gd->SetMarkerSize(0.7);

    // Extra check
    //if (s=="gj") tdrDraw(gms[c],"Pz",mcmarker[i],colors[i],kSolid,colors[i]);

    tex->SetTextColor(gd->GetMarkerColor());

    TGraphErrors *gm = gms[c];
    //if (gm && v[i]!="tt") {
    if (gm) {
      gm->SetMarkerStyle(mcmarker[i]);
      gm->SetMarkerColor(colors[i]);
      gm->SetLineColor(colors[i]);
      gm->SetMarkerSize(0.7);
    }
    if (plotX && ((plotMCup && s!="zb") || (plotZBup && s=="zb")))
      tdrDraw(gm,"Pz",mcmarker[i],colors[i],kSolid,colors[i]);

    c1->cd(2);
    TGraphErrors *gdr = (TGraphErrors*)gd->Clone();
    shiftGraph(gdr, f1dzb);

    // Clean bigger error bars
    for (int i = gdr->GetN()-1; i != -1; --i) {
      if (gdr->GetEY()[i]>0.2) gdr->RemovePoint(i);
    }
    gdr->Draw("SAMEPz");
    gdr->SetMarkerSize(0.7);

    TGraphErrors *gmr = (TGraphErrors*)gm->Clone();
    if (normalizeMCbyDataZB) shiftGraph(gmr, f1dzb);
    else shiftGraph(gmr, f1mzb);
    //if (v[i]=="sj" || v[i]=="zmm" || v[i]=="gj") {
    if (plotX) {
      // Clean bigger error bars
      for (int i = gmr->GetN()-1; i != -1; --i) {
	if (gmr->GetEY()[i]>0.2) gmr->RemovePoint(i);
      }
      if (plotX) tdrDraw(gmr,"Pz",mcmarker[i],colors[i],kSolid,colors[i]);
      gmr->SetMarkerSize(0.7);
    }

    gdr->Fit(fue,"QRN");
    TF1 *fuec = (TF1*)fue->Clone(Form("fue_%s",c));
    fuec->SetLineColor(colors[i]);
    fuec->Draw("SAME");

    gmr->Fit(fuemc,"QRN");
    TF1 *fuemcc = (TF1*)fuemc->Clone(Form("fuemc_%s",c));
    fuemcc->SetLineColor(colors[i]);
    //if (v[i]!="zb") 
    fuemcc->Draw("SAME");

    double muRef = muRef80 * sigmaMB / 80.; 
    double NpvDT = fs->Eval(muRef);
    double NpvMC = fsmc->Eval(muRef);
    if (i==0) {
      l->DrawLine(NpvDT,h2->GetMinimum(),NpvDT,h2->GetMaximum()-1);
      l->DrawLine(NpvMC,h2->GetMinimum(),NpvMC,h2->GetMaximum()-1);
    }
    if (v[i]=="sj" || v[i]=="ng") {
      rhoSJ = fue->Eval(NpvDT);
      rhoSJMC = fuemc->Eval(NpvMC);
      //rhoSJNR = fuemcc->Eval(NpvMC);
      eSJ = (fue->Eval(NpvDT*kRho)-rhoSJ);
      eSJMC = (fuemc->Eval(NpvMC*kRho)-rhoSJMC);
    }
    if (v[i]=="dj") {
      rhoDJ = fue->Eval(NpvDT);
      rhoDJMC = fuemc->Eval(NpvMC);
      eDJ = (fue->Eval(NpvDT*kRho)-rhoDJ);
      eDJMC = (fuemc->Eval(NpvMC*kRho)-rhoDJMC);
    }
    if (v[i]=="zmm") {
      rhoZmm = fue->Eval(NpvDT);
      rhoZmmMC = fuemc->Eval(NpvMC);
      eZmm = (fue->Eval(NpvDT*kRho)-rhoZmm);
      eZmmMC = (fuemc->Eval(NpvMC*kRho)-rhoZmmMC);
    }
    if (v[i]=="zj") {
      rhoZJ = fue->Eval(NpvDT);
      rhoZJMC = fuemc->Eval(NpvMC);
      eZJ = (fue->Eval(NpvDT*kRho)-rhoZJ);
      eZJMC = (fuemc->Eval(NpvMC*kRho)-rhoZJMC);
    }
    if (v[i]=="gj") {
      rhoGJ = fue->Eval(NpvDT);
      rhoGJMC = fuemc->Eval(NpvMC);
      eGJ = (fue->Eval(NpvDT*kRho)-rhoGJ);
      eGJMC = (fuemc->Eval(NpvMC*kRho)-rhoGJMC);
    }
    if (v[i]=="tt") {
      rhoTT = fue->Eval(NpvDT);
      rhoTTMC = fuemc->Eval(NpvMC);
      eTT = (fue->Eval(NpvDT*kRho)-rhoTT);
      eTTMC = (fuemc->Eval(NpvMC*kRho)-rhoTTMC);
    }
    if (v[i]=="zb") {
      rhoZB = f1dzb->Eval(NpvDT)/muRef;//NpvDT)/20.;
      rhoZBMC = f1mzb->Eval(NpvMC)/muRef;//NpvMC)/20.;
      //rhoZBNR = f3->Eval(NpvMC)/20.;
      eZB = (f1dzb->Eval(NpvDT*kRho)/(muRef*kRho)-rhoZB);
      eZBMC = (f1mzb->Eval(NpvMC*kRho)/(muRef*kRho)-rhoZBMC);
    }

    c1->cd(1);
  }

  TLegend *leg2 = tdrLeg(0.19,0.57,0.45,0.87);
  if (plotDJ) leg2->AddEntry(gds["dj"],"Incl. jet (p_{T}>450 GeV)","PL");
  if (plotTT) leg2->AddEntry(gds["tt"],"t#bar{t} l+jet","PL");
  if (plotGam) leg2->AddEntry(gds["gj"],"#gamma+jet (p_{T}>230 GeV)","PL");
  if (plotZmm) leg2->AddEntry(gds["zmm"],"Z+jet (p_{T,Z}>15 GeV)","PL");
  if (plotZJ) leg2->AddEntry(gds["zj"],"Z+jet (p_{T,Z}>45 GeV)","PL");
  if (plotSJ) leg2->AddEntry(gds["sj"],"MB jet (p_{T}>15 GeV)","PL");
  if (plotNG) leg2->AddEntry(gds["ng"],"MB jet (p_{T}>15 GeV)","PL");
  if (plotZB) leg2->AddEntry(gds["zb"],"Zero Bias","PL");

  tex->SetTextColor(kBlack);
  tex->SetTextSize(0.045);
  tex->DrawLatex(0.20,0.50,"|#eta| < 1.3");
  tex->SetTextSize(0.035);
  tex->SetTextColor(kBlack);
  if (plotDJ) {
    tex->SetTextColor(kBlack);
    //tex->DrawLatex(0.60,0.32,Form("#rho_{UE}(Incl. Jet) = %1.2f GeV",
    //				  rhoDJ));
    tex->DrawLatex(0.55,0.32,Form("#rho_{UE}(Incl. Jet) = %1.2f (%1.2f) GeV",
				  rhoDJ,rhoDJMC));
  }
  if (plotNG || plotSJ) {
    tex->SetTextColor(kBlue-9);
    //tex->DrawLatex(0.60,0.27,Form("#rho_{UE}(ZB Jet) = %1.2f GeV",
    //				  rhoSJ));
    tex->DrawLatex(0.55,0.12,Form("#rho_{UE}(MB jet) = %1.2f (%1.2f) GeV",
				  rhoSJ,rhoSJMC));
  }

  if (plotGam) {
    tex->SetTextColor(kBlue);
    //tex->DrawLatex(0.60,0.22,Form("#rho_{UE}(#gamma+jet) = %1.2f GeV",
    //				  rhoGJ));
    tex->DrawLatex(0.55,0.22,Form("#rho_{UE}(#gamma+jet) = %1.2f (%1.2f) GeV",
				  rhoGJ,rhoGJMC));
  }
  if (plotZmm) {
    tex->SetTextColor(kRed);
    //tex->DrawLatex(0.60,0.17,Form("#rho_{UE}(Z+jet) = %1.2f GeV",
    //				  rhoZmm));
    tex->DrawLatex(0.55,0.17,Form("#rho_{UE}(Z+jet) = %1.2f (%1.2f) GeV",
				  rhoZmm,rhoZmmMC));
  }
  if (plotZJ) {
    tex->SetTextColor(kRed+1);
    tex->DrawLatex(0.55,0.17,Form("#rho_{UE}(Z+jet) = %1.2f (%1.2f) GeV",
				   rhoZJ,rhoZJMC));
  }
  if (plotTT) {
    tex->SetTextColor(kGreen+2);
    //tex->DrawLatex(0.60,0.12,Form("#rho_{UE}(t#bar{t}) = %1.2f GeV",
    //				  rhoTT));
    tex->DrawLatex(0.55,0.27,Form("#rho_{UE}(t#bar{t}) = %1.2f (%1.2f) GeV",
				  rhoTT,rhoTTMC));
  }
  if (plotZB) {
    tex->SetTextColor(kGray+2);
    //tex->DrawLatex(0.60,0.07,Form("#rho_{PU}(ZB) = %1.2f GeV",
    //				  rhoZB));
    tex->DrawLatex(0.55,0.07,Form("#rho_{PU}(ZB) = %1.3f (%1.3f) GeV",
				  rhoZB,rhoZBMC));
  }

  c1->cd(2);

  TLegend *leg2d = tdrLeg(0.19,0.90-2*0.12,0.93,0.90);
  leg2d->SetNColumns(3);
  leg2d->SetTextSize(0.045*2.);
  if (plotDJ) leg2d->AddEntry(gms["dj"],"Dijet MC","PL");
  if (plotGam) leg2d->AddEntry(gms["gj"],"#gamma+jet MC","PL");
  if (plotZmm) leg2d->AddEntry(gms["zmm"],"Z+jet MC","PL");
  if (plotZJ) leg2d->AddEntry(gms["zj"],"Z+jet MC","PL");
  if (plotTT) leg2d->AddEntry(gms["tt"],"TT MC","PL");
  if (plotSJ) leg2d->AddEntry(gms["sj"],"NuGun MC","PL");
  if (plotNG) leg2d->AddEntry(gms["ng"],"NuGun MC","PL");

  curdir->cd();

  c1->Update();
  c1->SaveAs(Form("pdf/RhoVsNPV/RhovsNPV_%s.pdf",epoch.c_str()));


  // Dijet study
  if (false) {

    //TH1D *h22 = new TH1D("h22",";#LTN_{PV,good} - N_{HS}#GT;"
    TH1D *h22 = new TH1D("h22",";#LTN_{PV,good} - 1#GT;"
			"#LT#rho#GT- ZB (GeV)",30,0,30);
    h22->SetMaximum(3.5);
    h22->SetMinimum(-0.5);
    TCanvas *c22 = tdrDiCanvas("c22",h,h22,2,0); 
    h2->GetXaxis()->SetTitleOffset(0.85);
   
    int jcol[ntrg] = {kBlack, kRed, kOrange+2, kGreen+2, kCyan+2,
		      kCyan+3, kBlue, kBlue+1, kMagenta+1, kMagenta+2};
    
    c22->cd(1);

    TF1 *f12 = new TF1("f12","[0] + [1]*x + [2]*x*x",0,30);
    TF1 *fue = new TF1("fue2","[0]+[1]*x",0,30);
    gds["zb"]->Fit(f12,"QRN");
    f12->SetLineWidth(4);
    f12->SetLineColor(kGray);
    f12->Draw("SAME");
    if (plotZB) tdrDraw((TGraphErrors*)gds["zb"]->Clone(), "P", kOpenCircle, kBlack);

    // Treat one ZB interaction as "hard scatter"
    TGraphErrors *gzb = (TGraphErrors*)gds["zb"]->Clone("zbue");
    for (int i = 0; i != gzb->GetN(); ++i) {
      gzb->SetPoint(i, gzb->GetX()[i]-0.70, gzb->GetY()[i]);
    }

    if (plotZB) tdrDraw(gzb,"P",kFullCircle,kGray+1);
    TF1 *fuezb = (TF1*)fue->Clone("fue_zb");
    fuezb->SetLineColor(kGray+1);
    fuezb->Draw("SAME");

    c22->cd(2);

    TGraphErrors *gdrzb = (TGraphErrors*)gzb->Clone();
    shiftGraph(gdrzb, f12);
    gdrzb->Draw("SAMEP");

    gdrzb->Fit(fue,"QRN");
    TF1 *fuzb = (TF1*)fue->Clone("fue_zb");
    fuzb->SetLineColor(kGray+1);
    fuzb->Draw("SAME");

    for (int itrg = 0; itrg != ntrg; ++itrg) {
      
      c22->cd(1);
      
      const char *ct = Form("dj_%s",trg[itrg]);
      TGraphErrors *gd = gds[ct];
      if (plotDJ) tdrDraw(gd,"P",kFullCircle,jcol[itrg]);
      
      c22->cd(2);
      
      TGraphErrors *gdr = (TGraphErrors*)gd->Clone();
      shiftGraph(gdr, f12);
      gdr->Draw("SAMEP");

      gdr->Fit(fue,"QRN");
      TF1 *fuec = (TF1*)fue->Clone(Form("fue_%s",ct));
      fuec->SetLineColor(jcol[itrg]);
      fuec->Draw("SAME");
    } // for itrg

    c22->SaveAs(Form("pdf/RhoVsNPV/RhoVsNPV_dijet_%s.pdf",epoch.c_str()));
    //delete c2;
  } // jet trigger comparisons


  // Add legends
  c3->cd();
  TLegend *leg3a = tdrLeg(0.49,0.67,0.75,0.92);
  leg3a->SetHeader("MC");
  if (plotSJ) leg3a->AddEntry(gms["sj"]," ","PL");
  if (plotNG) leg3a->AddEntry(gms["ng"]," ","PL");
  if (plotDJ) leg3a->AddEntry(gms["dj"]," ","PL");
  if (plotGam) leg3a->AddEntry(gms["gj"]," ","PL");
  if (plotZmm) leg3a->AddEntry(gms["zmm"]," ","PL");
  if (plotZJ) leg3a->AddEntry(gms["zj"]," ","PL");
  if (plotTT) leg3a->AddEntry(gms["tt"]," ","PL");
  if (plotZB) leg3a->AddEntry(gms["zb"]," ","PL");
  
  TLegend *leg3b = tdrLeg(0.56,0.67,0.82,0.92);
  leg3b->SetHeader("Data");
  if (plotSJ) leg3b->AddEntry(gds["sj"],"Jet (p_{T} > 15 GeV)","PL");
  if (plotNG) leg3b->AddEntry(gds["ng"],"Jet (p_{T} > 15 GeV)","PL");
  if (plotDJ) leg3b->AddEntry(gds["dj"],"Jet (p_{T} > 450 GeV)","PL");
  if (plotGam) leg3b->AddEntry(gds["gj"],"#gamma+jet (p_{T} > 230 GeV)","PL");
  if (plotZmm) leg3b->AddEntry(gds["zmm"],"Z+jet","PL");
  if (plotZJ) leg3b->AddEntry(gds["zj"],"Z+jet","PL");
  if (plotTT) leg3b->AddEntry(gds["tt"],"t#bar{t}","PL");
  if (plotZB) leg3b->AddEntry(gds["zb"],"Zero Bias","PL");

  c4->cd();
  TLegend *leg4a = tdrLeg(0.49,0.67,0.75,0.92);
  leg4a->SetHeader("MC");
  if (plotSJ) leg4a->AddEntry(gms["sj"]," ","PL");
  if (plotNG) leg4a->AddEntry(gms["ng"]," ","PL");
  if (plotDJ)  leg4a->AddEntry(gms["dj"]," ","PL");
  if (plotGam) leg4a->AddEntry(gms["gj"]," ","PL");
  if (plotZmm) leg4a->AddEntry(gms["zmm"]," ","PL");
  if (plotZJ) leg4a->AddEntry(gms["zj"]," ","PL");
  if (plotTT) leg4a->AddEntry(gms["tt"]," ","PL");
  if (plotZB) leg4a->AddEntry(gms["zb"]," ","PL");
  
  TLegend *leg4b = tdrLeg(0.56,0.67,0.82,0.92);
  leg4b->SetHeader("Data");
  if (plotSJ) leg4b->AddEntry(gds["sj"],"Jet (p_{T} > 15 GeV)","PL");
  if (plotNG) leg4b->AddEntry(gds["ng"],"Jet (p_{T} > 15 GeV)","PL");
  if (plotDJ) leg4b->AddEntry(gds["dj"],"Jet (p_{T} > 450 GeV)","PL");
  if (plotGam) leg4b->AddEntry(gds["gj"],"#gamma+jet (p_{T} > 230 GeV)","PL");
  if (plotZmm) leg4b->AddEntry(gds["zmm"],"Z+jet","PL");
  if (plotZJ) leg4b->AddEntry(gds["zj"],"Z+jet","PL");
  if (plotTT) leg4b->AddEntry(gds["tt"],"t#bar{t}","PL");
  if (plotZB) leg4b->AddEntry(gds["zb"],"Zero Bias","PL");
  
  c3->SaveAs(Form("pdf/RhoVsNPV/RhoVsNPV_NpvVsMu_%s.pdf",epoch.c_str()));
  c4->SaveAs(Form("pdf/RhoVsNPV/RhoVsNPV_RhoVsMu_%s.pdf",epoch.c_str()));

  // Plot summary of measured offsets
  TH1D *h5u = tdrHist("h5u","#LT#rho#GT_{UE,meas} (GeV)",0,3.5,"Sample",0,7);
  TH1D *h5ud = tdrHist("h5ud","#LT#rho#GT_{UE,meas} (GeV)",0,3.5,"Sample",0,7);
  TH1D *h5um = tdrHist("h5um","#LT#rho#GT_{UE,meas,MC} (GeV)",0,3.5,"Sample",0,7);
  TH1D *h5d = tdrHist("h5d","Data-MC (GeV)",-0.7,0.2,"Sample",0,7);
  TCanvas *c5 = tdrDiCanvas("c5",h5u,h5d,4,0);

  const int ns = 7;
  double rho[ns] = {rhoZB, rhoSJ, rhoZmm, rhoZJ, rhoGJ, rhoTT, rhoDJ};
  double rhoMC[ns] = {rhoZBMC, rhoSJMC, rhoZmmMC, rhoZJMC, rhoGJMC, rhoTTMC, rhoDJMC};
  double err[ns] = {eZB, eSJ, eZmm, eZJ, eGJ, eTT, eDJ};
  double errMC[ns] = {eZBMC, eSJMC, eZmmMC, eZJMC, eGJMC, eTTMC, eDJMC};
  const char *name[ns] = {"ZeroBias","MinBias","Z+jet","Z+jet","#gamma+jet","t#bar{t}","QCD"};
  for (int i = 0; i != ns; ++i) {
    h5ud->SetBinContent(i+1,rho[i]);
    h5ud->SetBinError(i+1,err[i]);
    h5um->SetBinContent(i+1,rhoMC[i]);
    h5um->SetBinError(i+1,errMC[i]);
    h5d->SetBinContent(i+1,rho[i]-rhoMC[i]);
    h5d->SetBinError(i+1,fabs(err[i]-errMC[i]));
    h5d->GetXaxis()->SetBinLabel(i+1,name[i]);
  }

  c5->cd(1);
  tdrDraw(h5um,"H",kNone,kBlack);
  tdrDraw(h5ud,"EP",kFullCircle,kBlack,kSolid,-1,kNone);
  gPad->RedrawAxis();

  TLegend *leg5 = tdrLeg(0.70,0.88-0.06*2,1.00,0.88);
  leg5->AddEntry(h5ud,"Data","PL");
  leg5->AddEntry(h5um,"MC","F");

  c5->cd(2);
  tdrDraw(h5d,"HP",kFullCircle,kBlack);
  gPad->RedrawAxis();
  
  c5->SaveAs(Form("pdf/RhoVsNPV/RhoVsNPV_RhoUE_%s.pdf",epoch.c_str()));
}
