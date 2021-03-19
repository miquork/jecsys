// Purpose: compare multijets from 2016 and 2017 to 2018, data and MC
#include "TFile.h"
#include "TLine.h"

#include "../tdrstyle_mod15.C"

bool plotMC = false;//true;

TH1D *diffY(TH1D *h1, TH1D *h2, const char* name) {
  
  TH1D *h = (TH1D*)h1->Clone(name);
  h->Add(h1,h2,1,-1);
  h->Divide(h2);
  h->Scale(100.);
  
  return h;
}

void invertY(TH1D *h) {
  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    double y = h->GetBinContent(i);
    double ey = h->GetBinError(i);
    h->SetBinContent(i, (y!=0 ? 1./y : 0.));
    h->SetBinError(i, (y!=0 ? 1./y * ey/y : 0.));
  } // for i
} // invertY

void drawMultijet() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  //TFile *f16 = new TFile("../JERCProtoLab/Summer19UL16/L3Residual_multijet/multijet_UL2016FGH_jecV2_jerV1_v2.root","READ");
  TFile *f16 = new TFile("../JERCProtoLab/Summer19UL16/L3Residual_multijet/multijet_UL2016FGH_jecV2_jerV1_v3.root","READ"); // Laura
  assert(f16 && !f16->IsZombie());
  TFile *f16d = new TFile("../JERCProtoLab/Summer19UL16/L3Residual_multijet/multijet_Rebin2_20210312_UL2016FGH_jecV2_jerV1.root","READ"); // Minusk
  assert(f16d && !f16d->IsZombie());

  TFile *f17 = new TFile("rootfiles/multijet_Rebin2_20200526_UL2017BCDEF_SimpleL1V4JRV2_jecV4_jerV2.root","READ");
  assert(f17 && !f17->IsZombie());
	 
  //TFile *f18 = new TFile("rootfiles/multijet_UL2018ABCD_jecDTV4MCV2_jerUL17V2.root","READ");
  TFile *f18 = new TFile("rootfiles/multijet_Rebin2_20201207_UL2018ABCD_jecV5_jerV2.root","READ");
  assert(f18 && !f18->IsZombie());
  
  curdir->cd();

  //TH1D *hld18 = (TH1D*)f18->Get("Data/Pt30/MPF_leading_L1L2res");
  //TH1D *had18 = (TH1D*)f18->Get("Data/Pt30/MPF_ptave_L1L2res");
  //TH1D *hrd18 = (TH1D*)f18->Get("Data/Pt30/MPF_recoil_L1L2res");
  TH1D *hld18 = (TH1D*)f18->Get("Data/Pt30/MPF_leading_L1L2Res");
  TH1D *had18 = (TH1D*)f18->Get("Data/Pt30/MPF_ptave_L1L2Res");
  TH1D *hrd18 = (TH1D*)f18->Get("Data/Pt30/MPF_recoil_L1L2Res");
  assert(hld18);
  assert(had18);
  assert(hrd18);  

  TH1D *hlm18 = (TH1D*)f18->Get("MC/Pt30/MPF_leading_MG");
  TH1D *ham18 = (TH1D*)f18->Get("MC/Pt30/MPF_ptave_MG");
  TH1D *hrm18 = (TH1D*)f18->Get("MC/Pt30/MPF_recoil_MG");
  assert(hlm18);
  assert(ham18);
  assert(hrm18);

  //invertY(hrm18);

  TH1D *hld17 = (TH1D*)f17->Get("Data/Pt30/MPF_leading_L1L2Res");
  TH1D *had17 = (TH1D*)f17->Get("Data/Pt30/MPF_ptave_L1L2Res");
  TH1D *hrd17 = (TH1D*)f17->Get("Data/Pt30/MPF_recoil_L1L2Res");
  assert(hld17);
  assert(had17);
  assert(hrd17);

  TH1D *hlm17 = (TH1D*)f17->Get("MC/Pt30/MPF_leading_MG");
  TH1D *ham17 = (TH1D*)f17->Get("MC/Pt30/MPF_ptave_MG");
  TH1D *hrm17 = (TH1D*)f17->Get("MC/Pt30/MPF_recoil_MG");
  assert(hlm17);
  assert(ham17);
  assert(hrm17);

  /*
  TH1D *hld16 = (TH1D*)f16->Get("Data/Pt30/MPF_leading_L1L2res");
  TH1D *had16 = (TH1D*)f16->Get("Data/Pt30/MPF_ptave_L1L2res");
  TH1D *hrd16 = (TH1D*)f16->Get("Data/Pt30/MPF_recoil_L1L2res");
  assert(hld16);
  assert(had16);
  assert(hrd16);
  */

  TH1D *hld16 = (TH1D*)f16d->Get("Data/Pt30/MPF_leading_L1L2Res");
  TH1D *had16 = (TH1D*)f16d->Get("Data/Pt30/MPF_ptave_L1L2Res");
  TH1D *hrd16 = (TH1D*)f16d->Get("Data/Pt30/MPF_recoil_L1L2Res");
  assert(hld16);
  assert(had16);
  assert(hrd16);

  TH1D *hlm16 = (TH1D*)f16->Get("MC/Pt30/MPF_leading_MG");
  TH1D *ham16 = (TH1D*)f16->Get("MC/Pt30/MPF_ptave_MG");
  TH1D *hrm16 = (TH1D*)f16->Get("MC/Pt30/MPF_recoil_MG");
  assert(hlm16);
  assert(ham16);
  assert(hrm16);

  TH1D *hl18 = diffY(hld18,hlm18,"hl18");
  TH1D *ha18 = diffY(had18,ham18,"ha18");
  TH1D *hr18 = diffY(hrd18,hrm18,"hr18");

  TH1D *hl17 = diffY(hld17,hlm17,"hl17");
  TH1D *ha17 = diffY(had17,ham17,"ha17");
  TH1D *hr17 = diffY(hrd17,hrm17,"hr17");

  TH1D *hl16 = diffY(hld16,hlm16,"hl16");
  TH1D *ha16 = diffY(had16,ham16,"ha16");
  TH1D *hr16 = diffY(hrd16,hrm16,"hr16");
 

  double ptmin = 35;
  double ptmax = 4000;
  TH1D *hu = tdrHist("hu","#LTp_{T,lead}#GT / #LTp_{T,recoil}#GT",
		     0.90+1e-5,1.05,
		     //0.80,1.15, //0.89+1e-5,1.12-1e-5,
		     "p_{T,recoil}, p_{T,ave} or p_{T,lead} (GeV)",ptmin,ptmax);
  TH1D *hd = tdrHist("hd","Data/MC-1 (%)",-1.5,+2,//-3,+5,//-3.3,+3.3,
		     "p_{T,recoil}, p_{T,ave} or p_{T,lead} (GeV)",ptmin,ptmax);
  
  lumi_13TeV = "UL 16GH+17+18, 16.8+41.5+59.9 fb^{-1}";
  TCanvas *c1 = tdrDiCanvas("c1",hu,hd,4,11);

  TLine *l = new TLine();

  c1->cd(1);
  gPad->SetLogx();

  l->SetLineStyle(kDashed);
  l->DrawLine(ptmin,1,ptmax,1);
  l->SetLineStyle(kDotted);
  l->DrawLine(ptmin,1.01,ptmax,1.01);
  l->DrawLine(ptmin,0.99,ptmax,0.99);
  l->SetLineStyle(kDashDotted);
  l->DrawLine(114,0.93,114,1.02);

  if (plotMC) {
    tdrDraw(hlm18,"Pz",kFullCircle,kRed);
    tdrDraw(ham18,"PZ",kFullCircle,kGreen+2);
    tdrDraw(hrm18,"Pz",kFullCircle,kBlue);
    
    tdrDraw(hlm17,"Pz",kOpenCircle,kRed);
    tdrDraw(ham17,"PZ",kOpenCircle,kGreen+2);
    tdrDraw(hrm17,"Pz",kOpenCircle,kBlue);
    
    tdrDraw(hlm16,"Pz",kOpenDiamond,kRed);
    tdrDraw(ham16,"PZ",kOpenDiamond,kGreen+2);
    tdrDraw(hrm16,"Pz",kOpenDiamond,kBlue);
  }
  else {
    tdrDraw(hld18,"Pz",kFullCircle,kRed);
    tdrDraw(had18,"PZ",kFullCircle,kGreen+2);
    tdrDraw(hrd18,"Pz",kFullCircle,kBlue);
    
    tdrDraw(hld17,"Pz",kOpenCircle,kRed);
    tdrDraw(had17,"PZ",kOpenCircle,kGreen+2);
    tdrDraw(hrd17,"Pz",kOpenCircle,kBlue);
    
    tdrDraw(hld16,"Pz",kOpenDiamond,kRed);
    tdrDraw(had16,"PZ",kOpenDiamond,kGreen+2);
    tdrDraw(hrd16,"Pz",kOpenDiamond,kBlue);
  }

  TLegend *leg18 = tdrLeg(0.60,0.10,0.80,0.34);
  TLegend *leg17 = tdrLeg(0.55,0.10,0.75,0.34);
  TLegend *leg16 = tdrLeg(0.50,0.10,0.70,0.34);

  if (plotMC) { // plot MC
    leg18->SetHeader("18 (MC)");
    leg18->AddEntry(hlm18,"Leading","PLE");
    leg18->AddEntry(ham18,"Average","PLE");
    leg18->AddEntry(hrm18,"Recoil","PLE");
    leg17->SetHeader("17,");
    leg17->AddEntry(hlm17," ","PLE");
    leg17->AddEntry(ham17," ","PLE");
    leg17->AddEntry(hrm17," ","PLE");
    leg16->SetHeader("16,");
    leg16->AddEntry(hlm16," ","PLE");
    leg16->AddEntry(ham16," ","PLE");
    leg16->AddEntry(hrm16," ","PLE");
  }
  else {
    leg18->SetHeader("18 (DATA)");
    leg18->AddEntry(hld18,"Leading","PLE");
    leg18->AddEntry(had18,"Average","PLE");
    leg18->AddEntry(hrd18,"Recoil","PLE");
    leg17->SetHeader("17,");
    leg17->AddEntry(hld17," ","PLE");
    leg17->AddEntry(had17," ","PLE");
    leg17->AddEntry(hrd17," ","PLE");
    leg16->SetHeader("16,");
    leg16->AddEntry(hld16," ","PLE");
    leg16->AddEntry(had16," ","PLE");
    leg16->AddEntry(hrd16," ","PLE");
  }



  c1->cd(2);
  gPad->SetLogx();

  l->SetLineStyle(kDashed);
  l->DrawLine(ptmin,0,ptmax,0);
  l->SetLineStyle(kDotted);
  l->DrawLine(ptmin,+1,ptmax,+1);
  l->DrawLine(ptmin,-1,ptmax,-1);
  l->SetLineStyle(kDashDotted);
  l->DrawLine(114,-0.5,114,+1.5);

  /*
  tdrDraw(hl18,"Pz",kFullCircle,kRed);
  tdrDraw(ha18,"PZ",kFullCircle,kGreen+2);
  tdrDraw(hr18,"Pz",kFullCircle,kBlue);

  tdrDraw(hl17,"Pz",kOpenCircle,kRed);
  tdrDraw(ha17,"PZ",kOpenCircle,kGreen+2);
  tdrDraw(hr17,"Pz",kOpenCircle,kBlue);

  tdrDraw(hl16,"Pz",kOpenDiamond,kRed);
  tdrDraw(ha16,"PZ",kOpenDiamond,kGreen+2);
  tdrDraw(hr16,"Pz",kOpenDiamond,kBlue);
  */


  tdrDraw(ha18,"PZ",kFullCircle,kGreen+2);
  tdrDraw(ha17,"PZ",kOpenCircle,kGreen+2);
  tdrDraw(ha16,"PZ",kOpenDiamond,kGreen+2);

  /*
  tdrDraw(hl18,"Pz",kFullCircle,kRed);
  tdrDraw(hl17,"Pz",kOpenCircle,kRed);
  tdrDraw(hl16,"Pz",kOpenDiamond,kRed);
  */
  /*
  tdrDraw(hr18,"Pz",kFullCircle,kBlue);
  tdrDraw(hr17,"Pz",kOpenCircle,kBlue);
  tdrDraw(hr16,"Pz",kOpenDiamond,kBlue);
  */
  c1->SaveAs("pdf/drawMultijet_UL16vs17vsUL18.pdf");
}
