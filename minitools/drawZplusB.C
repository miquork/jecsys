// Purpose: draw inputs from Z+b analysis, compare to default Z+jet
// - check event counts vs pT in bins of alpha
// - check Z+jet MPF and pTbal vs pT in bins of alpha
// - check data, check MC, check ratio for Z+jet
// + purity of Z+b?
// run with 'root -l minitools/drawZplusB.C+g'
//
// To actually extract bJSF, use minitools/extractBSJF.C
#include "TFile.h"
#include "TGraphErrors.h"
#include "TF1.h"

#include "../tdrstyle_mod15.C"

#include <map>

using namespace std;

// Switches
//string sm = "mpfchs"; // method: mpfchs, ptchs (before v25)
string sm = "rmpf"; // method: mpfchs, ptchs (after v25)
//string sm = "ptchs";
string sb = "none"; // b-tag: none, loose, medium, tight
//string sb = "tight";
//string sb = "medium";

//void drawZplusB() {
void drawZplusB(string year="2018mm") {

  string sy = year;

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  // Sami's new results (nAOD with full JEC)
  //TFile *fzb = new TFile("rootfiles/jme_bplusZ_v1_2016GH.root","READ"); //X
  //TFile *fzb = new TFile("rootfiles/jme_bplusZ_v2_2016GH.root","READ");
  //TFile *fzb = new TFile("rootfiles/jme_bplusZ_v3_2016GH.root","READ");
  //TFile *fzb = new TFile("rootfiles/jme_bplusZ_v6_2016GH.root","READ");
  //TFile *fzb = new TFile("rootfiles/jme_bplusZ_v8_2016GH.root","READ");
  //TFile *fzb = new TFile("rootfiles/jme_bplusZ_v9_2016GH.root","READ");
  //TFile *fzb = new TFile("rootfiles/jme_bplusZ_v10_2016GH.root","READ");
  //TFile *fzb = new TFile("rootfiles/jme_bplusZ_2016BCDEFGH_v11.root","READ");
  //TFile *fzb = new TFile("rootfiles/jme_bplusZ_2016GH_v12.root","READ");
  //TFile *fzb = new TFile("rootfiles/jme_bplusZ_Run2Muons_v13.root","READ");
  //
  //TFile *fzb = new TFile("rootfiles/jme_bplusZ_2016Muons_v13.root","READ");
  //TFile *fzb = new TFile("rootfiles/jme_bplusZ_2017Muons_v14.root","READ");
  //TFile *fzb = new TFile("rootfiles/jme_bplusZ_2018Muons_v14.root","READ");
  //TFile *fzb = new TFile("rootfiles/jme_bplusZ_Run2Muons_v14.root","READ");
  //TFile *fzb = new TFile("rootfiles/jme_bplusZ_Run2_v18.root","READ");
  //TFile *fzb = new TFile("rootfiles/jme_bplusZ_Run2_v20.root","READ");
  //TFile *fzb = new TFile("rootfiles/jme_bplusZ_Run2_v25.root","READ");
  //TFile *fzb = new TFile("rootfiles/jme_bplusZ_Run2_Ver2.root","READ");
  TFile *fzb = new TFile("rootfiles/jme_bplusZ_Run2_Ver3.root","READ");
  assert(fzb && !fzb->IsZombie());
  curdir->cd();

  // Daniel's results as reference
  //TFile *fz = new TFile("../jecsys2016/rootfiles/zjet_combination_07Aug2017_Summer16_JECV15_Zmm_GH_2018-09-14.root","READ"); // L1L2Res
  //TFile *fz = new TFile("../jecsys2016/rootfiles/zjet_combination_07Aug2017_Summer16_JECV15_Zmm_BCDEFGH_2018-09-14.root","READ"); // L1L2Res
  //TFile *fz = new TFile("../jecsys2016/rootfiles/zjet_combination_07Aug2017_Summer16_JECV6_Zmm_GH_2018-03-06.root","READ"); // L1L2L3
  //TFile *fz = new TFile("rootfiles/zjet_combination_Autumn18_JECV16_Zmm_2019-08-19.root","READ"); { assert(fz && !fz->IsZombie()); fz->cd("Run2018ABCD"); fz = (TFile*)gDirectory; } // L1L2Res, L1L2L3Res
  //TFile *fz = new TFile("rootfiles/FullCombination_Zmm_17Sep2018_Autumn18_JECv17.root","READ"); { assert(fz && !fz->IsZombie()); fz->cd("Run2018ABCD"); fz = (TFile*)gDirectory; } // L1L2Res, L1L2L3Res
  //
  // [1] https://indico.cern.ch/event/758051/contributions/3143597/ Summer16 JEC V15 (nearest official release V11)
  // [2] https://indico.cern.ch/event/767001/contributions/3192422/ Fall17 JEC V31 (nearest official release V32)
  // [3] https://indico.cern.ch/event/837707/contributions/3512864/ Autumn18 JEC V16 (nearest official release V19)
  TFile *fz(0); string sjec("");
  if (sy=="2016mm") { // Missing L1L2L3res, V15 vs V11
    fz = new TFile("rootfiles/zjet_combination_07Aug2017_Summer16_JECV15_Zmm_BCDEFGH_2018-09-14.root","READ"); sjec="L1L2res"; assert(fz && !fz->IsZombie());
  }
  if (sy=="2016ee") { // Missing L1L2L3Res, V15 vs V11 
    fz = new TFile("rootfiles/zjet_combination_07Aug2017_Summer16_JECV15_Zee_BCDEFGH_2018-09-14.root","READ"); sjec="L1L2Res"; assert(fz && !fz->IsZombie());
  }
  if (sy=="2017mm") { // V32 vs V31
    fz = new TFile("rootfiles/zjet_combination_Fall17_JECV31_Zmm_BCDEF_2018-10-26.root","READ"); sjec="L1L2L3res"; assert(fz && !fz->IsZombie());
  }
  if (sy=="2017ee") { // V32 vs V31
    fz = new TFile("rootfiles/zjet_combination_Fall17_JECV31_Zee_BCDEF_2018-10-26.root","READ"); sjec="L1L2L3res"; assert(fz && !fz->IsZombie());
  }
  if (sy=="2018mm") { // V16 vs V19
    fz = new TFile("rootfiles/zjet_combination_Autumn18_JECV16_Zmm_2019-08-19.root","READ"); sjec="L1L2L3Res"; assert(fz && !fz->IsZombie()); fz->cd("Run2018ABCD"); fz = (TFile*)gDirectory;
  } // L1L2Res, L1L2L3Res
  if (sy=="2018ee") { // V16 vs V19
    fz = new TFile("rootfiles/zjet_combination_Autumn18_JECV16_Zmm_2019-08-19.root","READ"); sjec="L1L2L3Res"; assert(fz && !fz->IsZombie()); fz->cd("Run2018ABCD"); fz = (TFile*)gDirectory;
  } // L1L2Res, L1L2L3Res

  //{ assert(fz && !fz->IsZombie()); fz->cd("Run2018ABCD"); fz = (TFile*)gDirectory; } // L1L2Res, L1L2L3Res


  assert(fz && !fz->IsZombie());

  // L3Res response to add it to Daniel's L2res results for full JEC
  TF1 *fl3 = new TF1("fl3","[0]+[1]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[2]*((1-(-0.196332+0.307378*TMath::Log(x))/x)-(1-(-0.196332+0.307378*TMath::Log(208.))/208.))",30,2000);
  // From  Summer16_07Aug2017GH_V11_DATA_L2L3Residual_AK4PFchs.txt
  fl3->SetParameters(0.9894, 0.0747, -0.783);

  //int va[] = {10, 15, 20, 30, 40, 50, 60, 80, 100};
  //int va[] = {100, 80, 60, 50, 40, 30, 20, 15, 10};
  //int va[] = {30, 20, 15, 10};
  //int va[] = {100, 80, 60, 50, 40, 30};
  int va[] = {100, 30};
  const int na = sizeof(va)/sizeof(va[0]);

  string ssb = "X";
  if (sb=="none")  ssb = "";
  if (sb=="loose") ssb = "_btagCSVV2loose";
  if (sb=="medium") ssb = "_btagCSVV2medium";
  if (sb=="tight") ssb = "_btagCSVV2tight";
  const char *cb = ssb.c_str();

  // https://root.cern.ch/doc/master/classTColor.html
  //gStyle->SetPalette(kDeepSea);
  map<int,int> mcolor;
  mcolor[10] = kAzure-1;
  mcolor[15] = kAzure-2;
  mcolor[20] = kAzure-3;
  mcolor[30] = kAzure-4;
  mcolor[40] = kAzure-5;
  mcolor[50] = kAzure-6;
  mcolor[60] = kAzure-7;
  mcolor[80] = kAzure-8;
  mcolor[100] = kAzure-9;
  map<int,int> mcolor2;
  mcolor2[10] = kOrange-1;
  mcolor2[15] = kOrange-2;
  mcolor2[20] = kOrange-3;
  mcolor2[30] = kOrange-4;
  mcolor2[40] = kOrange-5;
  mcolor2[50] = kOrange-6;
  mcolor2[60] = kOrange-7;
  mcolor2[80] = kOrange-8;
  mcolor2[100] = kOrange-9;

  // ### Event counts ###

  TH1D *h1 = new TH1D("h1",";p_{T,Z} (GeV);Events",1670,30,1700);
  h1->GetXaxis()->SetMoreLogLabels();
  h1->GetXaxis()->SetNoExponent();
  h1->SetMinimum(0.5);//1e2);
  h1->SetMaximum(3e6);//1e9);

  TH1D *h2 = new TH1D("h2",";p_{T,Z} (GeV);Ratio to #alpha<1",1670,30,1700);
  h2->GetXaxis()->SetMoreLogLabels();
  h2->GetXaxis()->SetNoExponent();
  h2->SetMinimum(0.0);
  h2->SetMaximum(1.2);

  //lumi_13TeV = "2016GH, 16.2 fb^{-1}";
  //lumi_13TeV = "2016BCDEFGH, 36.5 fb^{-1}";
  //lumi_13TeV = "Run2, 137.5 fb^{-1}"; // 36.5+41.4+59.9
  lumi_13TeV = "2018ABCD, 59.9 fb^{-1}";
  TCanvas *c1 = tdrCanvas("c1",h1,4,11,kSquare);
  c1->SetLogx();
  c1->SetLogy();
  TLegend *leg1 = tdrLeg(0.70,0.90-0.06*na,0.95,0.90);

  TCanvas *c2 = tdrCanvas("c2",h2,4,0,kSquare);
  c2->SetLogx();
  //TLegend *leg2 = tdrLeg(0.70,0.45,0.95,0.90);
  TLegend *leg2 = tdrLeg(0.20,0.81,0.90,0.90);
  leg2->SetNColumns(5);

  TH1D *h100(0), *h0100(0);
  for (int i = 0; i != na; ++i) {

    int a = va[i];
    // mc = data (v1)
    //TH1D *h = (TH1D*)fzb->Get(Form("mc/eta13/statistics_zmmjet%s_a%d",cb,a));
    //TH1D *h = (TH1D*)fzb->Get(Form("data/eta13/statistics_zmmjet%s_a%d",cb,a));
    TH1D *h = (TH1D*)fzb->Get(Form("data/eta_00_13/statistics_zmmjet%s_a%d",cb,a));
    assert(h);
    //TH1D *h0 = (a<40 ? (TH1D*)fz->Get(Form("Data_RawNEvents_CHS_a%d_eta_00_13_L1L2Res",a)) : 0);
    //TH1D *h0 = (a<40 ? (TH1D*)fz->Get(Form("Data_RawNEvents_CHS_a%d_eta_00_13_L1L2L3",a)) : 0);
    TH1D *h0 = (a<40 ? (TH1D*)fz->Get(Form("Data_RawNEvents_CHS_a%d_eta_00_13_L1L2L3Res",a)) : 0); // L1L2L3Res


    if (!h100) h100 = (TH1D*)h->Clone("h100");
    h = (TH1D*)h->Clone(Form("hstat_%d",a));
    if (!h0100 && h0) h0100 = (TH1D*)h0->Clone("h0100");
    if (h0) h0 = (TH1D*)h0->Clone(Form("h0stat_%d",a));

    // Zere out bins where alpha*pT,Z<15 GeV
    for (int j = 1; j != h->GetNbinsX()+1; ++j) {
      double ptzmin = h->GetBinLowEdge(j);
      if (0.01*a*ptzmin<15.) {
	//h->SetBinContent(j, 0);
	//h->SetBinError(j, 0);
      }
    }

    // Ratio to alpha<1
    TH1D *hr = (TH1D*)h->Clone(Form("hr_%d",a));
    hr->Sumw2(); // created for some, but not all?
    hr->Divide(h,h100,1,1,"B");

    TH1D *hr0 = (h0 ? (TH1D*)h0->Clone(Form("hr_%d",a)) : 0);
    if (hr0 && h0) {
      hr0->Sumw2(); // created for some, but not all?
      hr0->Divide(h0,h0100,1,1,"B");
    }

    c1->cd();
    tdrDraw(h,"HE",kNone,mcolor[a],kSolid,-1,1001,mcolor[a]);
    if (h0) tdrDraw(h0,"H",kOpenCircle,mcolor[a]+1,kDashed,-1,kNone);
    leg1->AddEntry(h,Form("#alpha < %1.2f",0.01*a),"FL");
    if (h0 && i==na-1) {
      leg1->AddEntry(h0,"KIT 0.3","L");
      leg1->SetY1NDC(leg1->GetY1NDC()-0.06);
    }

    c2->cd();
    tdrDraw(hr,"HE",kNone,mcolor[a],kSolid,-1,1001,mcolor[a]);
    if (hr0) tdrDraw(hr0,"H",kOpenCircle,mcolor[a]+1,kDashed,-1,kNone);
    //if (a==100) leg2->AddEntry(h,"#alpha<1");
    //else        leg2->AddEntry(h,Form("<%1.2f",0.01*a),"FL");
    leg2->AddEntry(h,Form("<%1.2f",0.01*a),"FL");
  } // for i

  c1->cd();
  c1->RedrawAxis();
  leg1->Draw();

  c2->RedrawAxis();
  leg2->Draw();

  //c2->cd();
  //TF1 *f1 = new TF1("f1","30./x",30,1700);
  //f1->Draw("SAME");


  // ### MPF response or pT balance###

  TH1D *h3 = new TH1D("h3",";p_{T,Z} (GeV);MPF response",1670,30,1700);
  h3->GetXaxis()->SetMoreLogLabels();
  h3->GetXaxis()->SetNoExponent();
  h3->SetMinimum(0.950-1e-5); // mpfchs
  h3->SetMaximum(1.050+1e-5); // mpfchs
  if (sm == "ptchs") {
    h3->SetYTitle("p_{T} balance");
    h3->SetMinimum(0.82-1e-5);
    h3->SetMaximum(1.05+1e-5);
  }
  if (sb != "none") {
    h3->SetYTitle(Form("%s (%s b-tag)",h3->GetYaxis()->GetTitle(),sb.c_str()));
    h3->SetMinimum(0.920+1e-5);
    h3->SetMaximum(1.020+2e-5);
  }

  TH1D *h3dw = new TH1D("h3dw",";p_{T,Z} (GeV);Data/MC-1 (%)",1670,30,1700);
  h3dw->GetXaxis()->SetMoreLogLabels();
  h3dw->GetXaxis()->SetNoExponent();
  h3dw->SetMinimum(-3.0+1e-5); // mpfhcs
  h3dw->SetMaximum(+6.0-1e-5); // mpfchs
  if (sm == "ptchs") {
    h3dw->SetMinimum(-3.0+1e-5);
    h3dw->SetMaximum(+6.0-1e-5);
  }
  if (sb != "none") {
    h3dw->SetMinimum(-1.0+1e-5);
    h3dw->SetMaximum(+3.0-1e-5);
  }

  //TCanvas *c3 = tdrCanvas("c3",h3,4,0,kSquare);
  TCanvas *c3 = tdrDiCanvas("c3",h3,h3dw,4,0);
  c3->cd(2);
  gPad->SetLogx();
  c3->cd(1);
  gPad->SetLogx();

  TLegend *leg3 = tdrLeg(0.2,0.70,0.9,0.9);
  leg3->SetHeader("MC Data      #alpha_{max}              #alpha_{max}");
  leg3->SetNColumns(6);

  for (int i = 0; i != na; ++i) {

    int a = va[i];

    // From Sami
    TGraphErrors *gd = (TGraphErrors*)fzb->Get(Form("data/eta_00_13/%s_zmmjet%s_a%d",sm.c_str(),cb,a));
    assert(gd);
    // mc errors miscalculated?
    TGraphErrors *gm = (TGraphErrors*)fzb->Get(Form("mc/eta_00_13/%s_zmmjet%s_a%d",sm.c_str(),cb,a));
    assert(gm);
    // ratio
    //TGraphErrors *gr = (TGraphErrors*)fzb->Get(Form("ratio/eta_00_13/%s_zmmjet%s_a%d",sm.c_str(),cb,a)); // before v25
    TGraphErrors *gr = (TGraphErrors*)gd->Clone(); // after v25
    assert(gd->GetN() == gm->GetN());
    for (int i = 0; i != gr->GetN(); ++i) {
      assert(fabs(gd->GetX()[i] / gm->GetX()[i]-1)<0.05);
      gr->SetPoint(i, gd->GetX()[i], gd->GetY()[i] / gm->GetY()[i]);
    } // for i
    assert(gr);

    // From Daniel
    //string skit = (sm=="mpfchs" ? "MPF_CHS" : "PtBal_CHS"); // before v25
    string skit = (sm=="rmpf" ? "MPF_CHS" : "PtBal_CHS"); // after v25
    //string skit = (sm=="mpfchs" ? "MPF-notypeI_CHS" : "PtBal_CHS");
    const char *ck = skit.c_str();
    //TGraphErrors *gd0 = (a<40 ? (TGraphErrors*)fz->Get(Form("Data_%s_a%d_eta_00_13_L1L2L3",ck,a)) : 0);
    //TGraphErrors *gm0 = (a<40 ? (TGraphErrors*)fz->Get(Form("MC_%s_a%d_eta_00_13_L1L2L3",ck,a)) : 0);
    //TH1D *hgr0 = (a<50 ? (TH1D*)fz->Get(Form("Ratio_%s_a%d_eta_00_13_L1L2L3",ck,a)) : 0);
    //TGraphErrors *gd0 = (a<40 ? (TGraphErrors*)fz->Get(Form("Data_%s_a%d_eta_00_13_L1L2Res",ck,a)) : 0);
    //TGraphErrors *gm0 = (a<40 ? (TGraphErrors*)fz->Get(Form("MC_%s_a%d_eta_00_13_L1L2Res",ck,a)) : 0);
    // //TGraphErrors *gr0 = (a<40 ? (TGraphErrors*)fz->Get(Form("Ratio_MPF_CHS_a%d_eta_00_13_L1L2Res",a)) : 0);
    //TH1D *hgr0 = (a<40 ? (TH1D*)fz->Get(Form("Ratio_%s_a%d_eta_00_13_L1L2Res",ck,a)) : 0);
    //
    TGraphErrors *gd0 = (a<40 ? (TGraphErrors*)fz->Get(Form("Data_%s_a%d_eta_00_13_L1L2L3Res",ck,a)) : 0);
    TGraphErrors *gm0 = (a<40 ? (TGraphErrors*)fz->Get(Form("MC_%s_a%d_eta_00_13_L1L2L3Res",ck,a)) : 0);
    TH1D *hgr0 = (a<40 ? (TH1D*)fz->Get(Form("Ratio_%s_a%d_eta_00_13_L1L2L3Res",ck,a)) : 0);
    assert(a>30 || gd0);
    assert(a>30 || gm0);
    //assert(a>30 || gr0);
    assert(a>30 || hgr0);


    // Clean points with alpha*pT,z<15
    for (int j = gd->GetN()-1; j != -1; --j) {
      double ptzmin = h100->GetBinLowEdge(h100->FindBin(gd->GetX()[j]));
      if (0.01*a*ptzmin < 15.) gd->RemovePoint(j);
      else {
	double x = h100->GetBinCenter(h100->FindBin(gd->GetX()[j])-1);
	gd->SetPoint(j, x, gd->GetY()[j]);
      }
    }
    for (int j = gm->GetN()-1; j != -1; --j) {
      double ptzmin = h100->GetBinLowEdge(h100->FindBin(gm->GetX()[j]));
      if (0.01*a*ptzmin < 15.) gm->RemovePoint(j);
      else {
	double x = h100->GetBinCenter(h100->FindBin(gm->GetX()[j])-1);
	gm->SetPoint(j, x, gm->GetY()[j]);
      }
    }
    for (int j = gr->GetN()-1; j != -1; --j) {
      double ptzmin = h100->GetBinLowEdge(h100->FindBin(gr->GetX()[j]));
      if (0.01*a*ptzmin < 15.) gr->RemovePoint(j);
      else {
	//gr->SetPoint(j, gr->GetX()[j], 100.*(gr->GetY()[j]-1));
	double x = h100->GetBinCenter(h100->FindBin(gr->GetX()[j])-1);
	gr->SetPoint(j, x, 100.*(gr->GetY()[j]-1));
	gr->SetPointError(j, gr->GetEX()[j], 100.*gr->GetEY()[j]);
      }
    }
    // Clone ratio (is ok)
//     TGraphErrors *gr2 = (TGraphErrors*)gd->Clone("gr2");
//     for (int j = 0; j != gr2->GetN(); ++j) {
//     gr2->SetPoint(j, 0.5*(gd->GetX()[j]+gm->GetX()[j]),
//     	    100.*(gd->GetY()[j] / gm->GetY()[j] - 1));
//     gr2->SetPointError(j, sqrt(pow(gd->GetEX()[j],2)+
//     			 pow(0.5*(gd->GetX()[j]-gm->GetX()[j]),2)),
//     		 100.*sqrt(pow(gd->GetEY()[j],2)+pow(gm->GetEY()[j],2)));
//     }
    
    // Change Daniel's result to percentages
    //if (gr0) {
      //for (int j = gr0->GetN()-1; j != -1; --j) {
      //for (int j = 0; j != gr0->GetN(); ++j) {
	//double ptzmin = h0100->GetBinLowEdge(h0100->FindBin(gr0->GetX()[j]));
	//if (0.01*a*ptzmin < 15.) true;//gr0->RemovePoint(j);
	//else {
	//gr0->SetPoint(j, gr0->GetX()[j], 100.*(gr0->GetY()[j]-1));
	//gr0->SetPointError(j, gr0->GetEX()[j], 100.*gr0->GetEY()[j]);
	//}
    //}
    //}

    // Fix Daniel's results for missing L3Res
    if (gd0) {
      for (int j = gd0->GetN()-1; j != -1; --j) {
	double x = gd0->GetX()[j];
	double r3 = fl3->Eval(x);
	gd0->SetPoint(j, x, gd0->GetY()[j]/r3);
      }
    }
    if (hgr0) {
      for (int j = 1; j != hgr0->GetNbinsX()+1; ++j) {
	double x = hgr0->GetBinCenter(j);
	double r3 = fl3->Eval(x);
	hgr0->SetBinContent(j, 100.*(hgr0->GetBinContent(j)/r3-1));
	hgr0->SetBinError(j, 100.*hgr0->GetBinError(j));
      }
    }

    c3->cd(1);
    if (gd0) tdrDraw(gd0,"Pz",kFullSquare,mcolor2[a]);
    if (gm0) tdrDraw(gm0,"Pz",kOpenSquare,mcolor2[a]);
    tdrDraw(gd,"P",kFullCircle,mcolor[a]);//,kSolid,-1,1001,mcolor[a]);
    tdrDraw(gm,"P",kOpenCircle,mcolor[a]);//,kSolid,-1,1001,mcolor[a]);
    
    leg3->AddEntry(gm," ","PL");
    leg3->AddEntry(gd,Form("<%1.2f",0.01*a),"PL");

    c3->cd(2);
    //if (gr0) tdrDraw(gr0,"PL",kOpenSquare,mcolor2[a]);
    if (hgr0) tdrDraw(hgr0,"Pz",kOpenSquare,mcolor2[a]);
    tdrDraw(gr,"P",kFullCircle,mcolor[a]);//,kSolid,-1,1001,mcolor[a]);
    //tdrDraw(gr2,"PL",kOpenCircle,mcolor[a]);//,kSolid,-1,1001,mcolor[a]);
  } // for i

  c1->SaveAs(Form("pdf/drawZplusB_events_%s.pdf",sb.c_str()));
  c2->SaveAs(Form("pdf/drawZplusB_alphafraction_%s.pdf",sb.c_str()));
  c3->SaveAs(Form("pdf/drawZplusB_%s_%s.pdf",sm.c_str(),sb.c_str()));

  // §14 bins in pT, from 30 to 1.5 TeV
  // b-jet purity 5.7% (no tag), 24.7% (loose), 58.5% (medium), 86.5% (tight)
}
