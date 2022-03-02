// Purpose: test JER SF text file produced by minitools/jerSF
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "TMath.h"

#include "tdrstyle_mod15.C"

#include <iostream>
#include <vector>

using namespace std;

const bool debug = true;

void testJERSF(string filename) {

  setTDRStyle();

  //filename = "../JECDatabase/textFiles/Summer19UL18_V5_MC/Summer19UL18_V5_MC_L2Relative_AK4PFchs.txt"; // L2Res example, runs
  //filename = "../JRDatabase/textFiles/Summer19UL18_JRV2_MC/Summer19UL18_JRV2_MC_PtResolution_AK4PFchs.txt"; // JER example, doesn't run
  cout << "Open " << filename << endl << flush;
  JetCorrectorParameters *p = new JetCorrectorParameters(filename.c_str());
  vector<JetCorrectorParameters> v;
  v.push_back(*p);
  FactorizedJetCorrector *jer = new FactorizedJetCorrector(v);

  if (debug) {
    jer->setJetEta(0.);
    jer->setJetPt(80.);
    jer->setRho(20.85); // UL18 Z+jet vs pT at [60,130] GeV, |eta|<1.3
    cout << "JER(eta=0,pt=80,rho=20.85)="<<jer->getCorrection()<<endl;
    jer->setJetEta(0.);
    jer->setJetPt(80.);
    jer->setRho(5.); // UL18 Z+jet vs pT at [60,130] GeV, |eta|<1.3
    cout << "JER(eta=0,pt=80,rho=5)="<<jer->getCorrection()<<endl;
  }

  const int ndiv = 1;//5;
  const int nbins = 60*ndiv;
  const double deta = TMath::TwoPi()/72./ndiv;
  const double maxeta = nbins*deta; // 5.236
  const double rho = 20.85;
  TH1D *h8 = new TH1D("h8",";#eta_{jet};JER SF",nbins,0.,maxeta);
  TH1D *h80 = new TH1D("h80",";#eta_{jet};JER SF",nbins,0.,maxeta);
  TH1D *hxmax = new TH1D("hxmax",";#eta_{jet};JER SF",nbins,0.,maxeta);
  for (int i = 0; i != h80->GetNbinsX()+1; ++i) {

    double eta = h80->GetBinCenter(i);
    double etamin = h80->GetBinLowEdge(i);

    jer->setJetEta(eta);
    jer->setJetPt(8.);
    jer->setRho(rho);
    h8->SetBinContent(i, jer->getCorrection());

    jer->setJetEta(eta);
    jer->setJetPt(80.);
    jer->setRho(rho);
    h80->SetBinContent(i, jer->getCorrection());

    double ptmax = min(3500.,6500./cosh(etamin));
    jer->setJetEta(eta);
    jer->setJetPt(ptmax);
    jer->setRho(rho);
    hxmax->SetBinContent(i, jer->getCorrection());
  }

  TH1D *h = tdrHist("h","JER SF",0.85,1.75,"#eta_{jet}",0,5.2);
  lumi_13TeV = "UL2018, 59.9 fb^{-1}";
  TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);
  tdrDraw(h8,"HPz",kOpenDiamond,kRed,kSolid,-1,kNone);
  tdrDraw(h80,"HPz",kFullCircle,kBlack,kSolid,-1,kNone);
  tdrDraw(hxmax,"HPz",kOpenDiamond,kBlue,kSolid,-1,kNone);
  
  c1->SaveAs("pdf/jersf/testJERSF.pdf");
} // testJERSF
