// Created by Mikko Voutilainen, on Jan 15th, 2020
// (initially copied from drawL1.C)
// Purpose is to compare simple and complex parameterizations of offset
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TF1.h"
#include "TMultiGraph.h"
#include "TLatex.h"
#include "TMinuit.h"
#include "TMatrixD.h"
#include "TMath.h"
#include "TLine.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "tools.h"
#include "tdrstyle_mod14.C"

#include <string>

using namespace std;

// global constants used also in getCorrPt
//const double gRho = 20.47;
const double gRho = 22.49;
const double gRho0 = 1.519;
const double gPt = 100;//30;//15;//500;//100;//30;//15.;//8;//5
const double gJetArea = TMath::Pi()*0.4*0.4;

FactorizedJetCorrector *getJEC(string version, string version2="") {

  string s;
  const char *cd = "CondFormats/JetMETObjects/data";
  //const char *cd = "rootfiles/JECULMC2017v2";
  //const char *gt = "Summer16_07Aug2017";
  //const char *gt = "Autumn18";
  const char *gt = "";

  FactorizedJetCorrector *jecl1;
  s = Form("%s/%s%s.txt",cd,gt,version.c_str());
  cout << s << endl << flush;
  JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
  vector<JetCorrectorParameters> v;
  v.push_back(*l1);

  if (version2!="") {
    s = Form("%s/%s%s.txt",cd,gt,version2.c_str());
    cout << s << endl << flush;
    JetCorrectorParameters *l2 = new JetCorrectorParameters(s);
    v.push_back(*l2);
  }

  jecl1 = new FactorizedJetCorrector(v);

  return jecl1;
}

double getCorrPt(double ptraw, double eta, FactorizedJetCorrector *jec) {

  jec->setJetPt(ptraw);
  jec->setJetEta(eta);
  jec->setRho(gRho);
  jec->setJetA(gJetArea);
  double c = jec->getCorrection();  

  return (c*ptraw);
}
FactorizedJetCorrector *_jec(0);
Double_t fGetCorrPt(Double_t *x, Double_t *p) {
  return getCorrPt(x[0], p[0],_jec);
}
TF1 *_f1(0);
double getRawPt(double ptcorr, double eta, FactorizedJetCorrector *jec) {
  if (!_f1) _f1 = new TF1("f1",fGetCorrPt,0,6500,1);
  _f1->SetParameter(0,eta);
  _jec = jec;
  double ptraw = _f1->GetX(ptcorr,1,6500);
  
  return ptraw;
} 

void drawSimpleL1() {

  setTDRStyle();

  // L1 corrections only
  FactorizedJetCorrector *jecrcmc = getJEC("Fall17_17Nov2017_V32_MC_L1RC_AK4PFchs"); // subtract 1.519 GeV
  FactorizedJetCorrector *jecl1mc = getJEC("Fall17_17Nov2017_V32_MC_L1FastJet_AK4PFchs");
  //
  FactorizedJetCorrector *jecrcnmc = getJEC("Summer19UL17_V1_SimpleL1_MC_L1RC_AK4PFchs"); // subtracts 2.00 GeV => 0.5 GeV difference to EOY17
  FactorizedJetCorrector *jecl1smc = getJEC("Summer19UL17_V1_SimpleL1_MC_L1FastJet_AK4PFchs"); 
  FactorizedJetCorrector *jecl1cmc = getJEC("Summer19UL17_V1_ComplexL1_MC_L1FastJet_AK4PFchs"); 

  // L1L2 corrections for mapping pTgen -> pTraw for L1 alone
  FactorizedJetCorrector *jec2rcmc = getJEC("Fall17_17Nov2017_V32_MC_L1RC_AK4PFchs","Fall17_17Nov2017_V32_MC_L2Relative_AK4PFchs");
  FactorizedJetCorrector *jec2l1mc = getJEC("Fall17_17Nov2017_V32_MC_L1FastJet_AK4PFchs","Fall17_17Nov2017_V32_MC_L2Relative_AK4PFchs");
  //
  FactorizedJetCorrector *jec2rcnmc = getJEC("Summer19UL17_V1_SimpleL1_MC_L1RC_AK4PFchs","Summer19UL17_V1_SimpleL1_MC_L2Relative_AK4PFchs");
  FactorizedJetCorrector *jec2l1smc = getJEC("Summer19UL17_V1_SimpleL1_MC_L1FastJet_AK4PFchs","Summer19UL17_V1_SimpleL1_MC_L2Relative_AK4PFchs");
  FactorizedJetCorrector *jec2l1cmc = getJEC("Summer19UL17_V1_ComplexL1_MC_L1FastJet_AK4PFchs","Summer19UL17_V1_ComplexL1_MC_L2Relative_AK4PFchs");

  TH1D *h = new TH1D("h",";#eta;off(#LT#rho#GT) / A#times(#LT#rho#GT-#rho_{0})",
		     94,-4.7,4.7);
  TH1D *hrcmc = (TH1D*)h->Clone("hrcmc");
  TH1D *hl1mc = (TH1D*)h->Clone("hl1mc");
  //
  TH1D *hrcnmc = (TH1D*)h->Clone("hrcnmc");
  TH1D *hl1smc = (TH1D*)h->Clone("hl1smc");
  TH1D *hl1cmc = (TH1D*)h->Clone("hl1cmc");

  const int njec = 5;
  FactorizedJetCorrector *vjec[njec] =
    {jecrcmc, jecl1mc, jecrcnmc, jecl1smc, jecl1cmc};
  FactorizedJetCorrector *vjec2[njec] =
    {jec2rcmc, jec2l1mc, jec2rcnmc, jec2l1smc, jec2l1cmc};
  TH1D *vh[njec] = {hrcmc, hl1mc, hrcnmc, hl1smc, hl1cmc};
  //int colors[njec] = {kGray, kBlack, kRed, kOrange+2, kBlue};
  //int style[njec] = {kSolid, kSolid, kSolid, kSolid, kSolid};

  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    
    double eta = h->GetBinCenter(i);
    for (int j = 0; j != njec; ++j) {

      FactorizedJetCorrector *jec2 = vjec2[j];
      double ptraw = getRawPt(gPt, eta, jec2);

      FactorizedJetCorrector *jec = vjec[j];
      jec->setJetEta(eta);
      jec->setJetPt(ptraw);//gPt);
      jec->setRho(gRho);
      
      jec->setJetA(gJetArea);
      double c = jec->getCorrection();
      // c = 1-off/pt => off = (1-c)*pt
      double off = (1-c)*ptraw;//gPt;
    
      TH1D *h = vh[j];
      h->SetBinContent(i, off / (gJetArea*(gRho-gRho0)));
    }
  }


  double ptmax = 1200.*2.;
  TH1D *hup = new TH1D("hup",";#eta;offset(#LT#rho#GT) / A#times(#LT#rho#GT-#rho_{0})", 94,-4.7,4.7);
  hup->SetMinimum(-0.2+1e-4);
  hup->SetMaximum(1.7-1e-4);

  TH1D *hdw = new TH1D("hdw",";#eta;Data / MC",94,-4.7,4.7);
  hdw->SetMinimum(-0.3);
  hdw->SetMaximum(+0.5);

  lumi_13TeV = "UL2017";

  TCanvas *c1 = tdrDiCanvas("c1",hup,hdw,4,11);

  c1->cd(1);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kBlack);
  l->DrawLine(-4.7,0,4.7,0);

  double ymc = hl1mc->Integral()/hl1mc->GetNbinsX();
  double yrc = hrcmc->Integral()/hrcmc->GetNbinsX();

  l->SetLineStyle(kDotted);

  // MC red, data black, RC blue
  tdrDraw(hrcmc,"HIST",0,kRed,kDotted,-1,0,0); // EOY17
  tdrDraw(hrcnmc,"HIST",0,kRed,kSolid,-1,0,0);
  tdrDraw(hl1mc,"HIST",0,kOrange+2,kDotted,-1,0,0); // EOY17
  tdrDraw(hl1cmc,"HIST",0,kBlue,kSolid,-1,0,0);
  tdrDraw(hl1smc,"HIST",0,kGreen+2,kSolid,-1,0,0);

  hrcmc->SetLineWidth(3);
  hl1mc->SetLineWidth(3);
  hrcnmc->SetLineWidth(3);
  hl1cmc->SetLineWidth(2);
  hl1smc->SetLineWidth(2);

  TLegend *leg1 = tdrLeg(0.43,0.50,0.63,0.75);
  leg1->AddEntry(hrcmc,"RC EOY17","L");
  leg1->AddEntry(hrcnmc,"RC ","L");
  leg1->AddEntry(hl1cmc,"L1 complex","L");
  leg1->AddEntry(hl1smc,"L1 simple","L");
  leg1->AddEntry(hl1mc,"L1 EOY17","L");

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);
  tex->SetTextColor(kGray+2);
  tex->DrawLatex(0.40,0.85,Form("p_{T,corr}=%1.0f GeV, #LT#rho#GT = %1.1f GeV",gPt,gRho));
  tex->DrawLatex(0.40,0.79,Form("Anti-k_{T} R=0.4, #rho_{0}= %1.3f GeV",gRho0));

  gPad->RedrawAxis();


  c1->cd(2);

  l->SetLineStyle(kDashed);
  l->SetLineColor(kBlack);
  l->DrawLine(-4.7,0,4.7,0);

  TH1D *hrc = (TH1D*)hrcmc->Clone("hrc");
  TH1D *hl1 = (TH1D*)hl1mc->Clone("hl1");
  TH1D *hrcn = (TH1D*)hrcnmc->Clone("hrcn");
  TH1D *hl1s = (TH1D*)hl1smc->Clone("hl1s");
  TH1D *hl1c = (TH1D*)hl1cmc->Clone("hl1c");

  hrc->Add(hl1s,-1);
  hrcn->Add(hl1s,-1);
  hl1->Add(hl1s,-1);
  hl1c->Add(hl1s,-1);

  hdw->SetYTitle("Minus L1S");

  tdrDraw(hrc,"HIST",0,kRed,kDotted,-1,0,0);
  tdrDraw(hrcn,"HIST",0,kRed,kSolid,-1,0,0);
  tdrDraw(hl1,"HIST",0,kOrange+2,kDotted,-1,0,0);
  tdrDraw(hl1c,"HIST",0,kBlue,kSolid,-1,0,0);
    
  hrc->SetLineWidth(3);
  hrcn->SetLineWidth(3);
  hl1s->SetLineWidth(2);
  hl1c->SetLineWidth(2);
    
  c1->SaveAs(Form("pdf/drawSimpleL1_UL2017_pt%d.pdf",int(gPt+0.5)));


  // Solve offset difference as a fraction of jet pT
  c1->cd(2);

  
}
