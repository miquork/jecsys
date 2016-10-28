#include "TF1.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "tdrstyle_mod15.C"

#include <string>
#include <vector>
#include <iostream>

using namespace std;

FactorizedJetCorrector *initJEC(string s) {
  cout << s << endl << flush;
  JetCorrectorParameters *p = new JetCorrectorParameters(s.c_str());
  vector<JetCorrectorParameters> v;
  v.push_back(*p);
  return (new FactorizedJetCorrector(v));
}

// Helper function to solve pTgen from pTgen = pTreco - offset
FactorizedJetCorrector *_jec(0);
Double_t fPtCorr(Double_t *x, Double_t *p) {
  assert(_jec);
  
  _jec->setJetPt(*x);
  _jec->setJetEta(p[0]);
  _jec->setRho(p[1]);
  _jec->setJetA(p[2]);
  double c = _jec->getCorrection();
    
  // c = 1 - off/pt => pt - c*pt = off =>
  //offset = (1-c)*pt;

  return (c * (*x));
}

TF1 *_fptc(0);
double getOffset(double ptgen, double eta, double rho, double r,
		 FactorizedJetCorrector *jec) {
  assert(jec);
  _jec = jec;
  if (!_fptc) {
    _fptc = new TF1("fptc",fPtCorr,5,6500,3);
  }
  _fptc->SetParameter(0, eta);
  _fptc->SetParameter(1, rho);
  _fptc->SetParameter(2, TMath::Pi()*r*r);
  double ptreco = _fptc->GetX(ptgen, ptgen, ptgen+30);
  double offset = ptreco - ptgen;

  /*
  // Iterate to solve offset for ptgen
  double offset_old(0), offset(0);
  double doffset(10);
  while (fabs(doffset)>0.01) {

    double pt = ptgen + offset;
    jec->setJetPt(pt);
    jec->setJetEta(eta);
    jec->setRho(rho);
    jec->setJetA(TMath::Pi()*r*r);
    double c = jec->getCorrection();
    
  // c = 1 - off/pt => pt - c*pt = off =>
    
    offset = (1-c)*pt;
    doffset = offset - offset_old;
    offset_old = offset;
  }
  */

  return offset;
}


// draw offset pT dependendence and differences between PF and PFchs
// call mk_drawOffsetPt() to run
void drawOffsetPt() {
  
  setTDRStyle();

  double rho = 10.;
  //double r = 0.5;
  //const char *cd = "CondFormats/JetMETObjects/data/Winter14_V8_";
  double r = 0.4;
  const char *cd = "CondFormats/JetMETObjects/data/Spring16_25nsV7G_";
    
  vector<string> s;
  /*
  s.push_back("MC_RC_AK5PF.txt");
  s.push_back("DATA_RC_AK5PF.txt");
  s.push_back("MC_L1FastJet_AK5PF.txt");
  s.push_back("DATA_L1FastJet_AK5PF.txt");
  //
  s.push_back("MC_RC_AK5PFchs.txt");
  s.push_back("DATA_RC_AK5PFchs.txt");
  s.push_back("MC_L1FastJet_AK5PFchs.txt");
  s.push_back("DATA_L1FastJet_AK5PFchs.txt");
  */
  s.push_back("MC_L1RC_AK4PF.txt");
  s.push_back("DATA_L1RC_AK4PF.txt");
  s.push_back("MC_L1FastJet_AK4PF.txt");
  s.push_back("DATA_L1FastJet_AK4PF.txt");
  //
  s.push_back("MC_L1RC_AK4PFchs.txt");
  s.push_back("DATA_L1RC_AK4PFchs.txt");
  s.push_back("MC_L1FastJet_AK4PFchs.txt");
  s.push_back("DATA_L1FastJet_AK4PFchs.txt");

  const int ns = 8;
  assert(ns==s.size());
  double markers[ns] = {kOpenSquare, kFullSquare, kOpenCircle, kFullCircle,
			kOpenSquare, kFullSquare, kOpenCircle, kFullCircle};
  double colors[ns] = {kGreen+2, kOrange+2, kBlue, kRed,
		       kGreen+2, kOrange+2, kBlue, kRed};
  const char *labels[ns] =
    {"MC RC", "Data RC", "MC true", "Data true",
     "MC RC CHS", "Data RC CHS", "MC true CHS", "Data true CHS"};

  
  vector<Double_t> ptbins;
  for (int i = 0; 5.*pow(1.2,i)<250.; ++i) {
    ptbins.push_back(5.*pow(1.2,i));
  }

  //TH1D *h = new TH1D("h",";p_{T,gen};offset / (#rho#timesA)",245,5,250);
  TH1D *h = new TH1D("h",";p_{T,corr} (GeV);Offset / (#rho#timesA)",
		     ptbins.size()-1, &ptbins[0]);
  h->SetMinimum(0);//2);
  h->SetMaximum(2);//12);
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  //h->Draw();

  lumi_13TeV = "19.7 fb^{-1}";
  TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);
  //c1->SetLogx();

  TLegend *leg = tdrLeg(0.5,0.70,0.8,0.9);
  //TLegend *leg1 = tdrLeg(0.3,0.70,0.6,0.9);
  //TLegend *leg2 = tdrLeg(0.6,0.70,0.9,0.9);

  for (unsigned int is = 0; is != s.size(); ++is) {

    const char *cs = s[is].c_str();
    FactorizedJetCorrector *jec = initJEC(Form("%s%s",cd,cs));

    TH1D *ho = (TH1D*)h->Clone(Form("ho_%s",cs));
    for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    
      //ho->SetBinContent(i, getOffset(h->GetBinCenter(i),0,10.,0.5,jec));
      // average over |eta|<1.3
      double pt = h->GetBinCenter(i);

      double sumo(0);
      const int neta = 26;//1;//26;
      for (int ieta = 0; ieta != neta; ++ieta) {
	double eta = -1.3 + 2.6/neta*(ieta+0.5);
	sumo += 1./neta * getOffset(pt,eta,rho,r,jec);
      }
      ho->SetBinContent(i, sumo / (rho * TMath::Pi()*r*r));
    }
    ho->SetMarkerStyle(markers[is]);
    ho->SetMarkerColor(colors[is]);
    ho->SetLineColor(colors[is]);
    ho->Draw("SAMEPL");
    //if (is==0||is==1||is==4||is==5) leg->AddEntry(ho,labels[is],"PL");
    if (is<4) leg->AddEntry(ho,labels[is],"PL");
    //if (is<4)  leg1->AddEntry(ho,labels[is],"PL");
    //if (is>=4) leg2->AddEntry(ho,labels[is],"PL");
  } // for i

  c1->RedrawAxis();
  // For some reason calling SetLogx earlier drops result precision?
  c1->SetLogx();
  
  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  tex->DrawLatex(0.20,0.75,"Anti-k_{T} R=0.5");
  tex->DrawLatex(0.20,0.70,Form("#rho = %1.1f GeV/A",rho)); 
  tex->DrawLatex(0.20,0.65,"|#eta| < 1.3"); 

  c1->SaveAs("pdf/drawOffsetPt.pdf");
}
