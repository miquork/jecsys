// Purpose: Estimate minimum safe jet pT as twice the current offset level
//          (or use offset*(1+2*RMS), where RMS=sqrt(offset/GeV)?)
//          Compare L1RC to L1MC to check stability
//          Plot vs eta to see regions approaching safe limits
// run with 'root -l -b -q minitools/mk_estimateMinJetPt.C'
#include "TMath.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
//#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "../tdrstyle_mod15.C"

FactorizedJetCorrector *getJEC(string version) {

  string s;
  const char *cd = "CondFormats/JetMETObjects/data";
  const char *gt = "Autumn18";

  FactorizedJetCorrector *jecl1;
  s = Form("%s/%s%s.txt",cd,gt,version.c_str());
  cout << s << endl << flush;
  JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
  vector<JetCorrectorParameters> v;
  v.push_back(*l1);
  jecl1 = new FactorizedJetCorrector(v);

  return jecl1;
}


void estimateMinJetPt() {

  setTDRStyle();

  double veta[] =
    {-5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664,
     -3.489, -3.314, -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172,
     -2.043, -1.93, -1.83, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305,
     -1.218, -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522,
     -0.435, -0.348, -0.261, -0.174, -0.087, 0, 0.087, 0.174, 0.261, 0.348,
     0.435, 0.522, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305,
     1.392, 1.479, 1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5,
     2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191,
     4.363, 4.538, 4.716, 4.889, 5.191};
  const int neta = sizeof(veta)/sizeof(veta[0])-1;

  const double mu = 35.39; // 2018D jt500 htrpu average 35.39
  const double rho = 0.55*mu; // 2018D jt500 hrho average 19.31

  TH1D *hrc = new TH1D("hrc",";Jet #eta;RC offset (GeV)",neta,veta);
  TH1D *hl1 = new TH1D("hl1",";Jet #eta;L1MC offset (GeV)",neta,veta);

  FactorizedJetCorrector *jecrc = getJEC("_V8_MC_L1RC_AK4PFchs");
  FactorizedJetCorrector *jecl1 = getJEC("_V8_MC_L1FastJet_AK4PFchs");

  for (int i = 1; i != hrc->GetNbinsX()+1; ++i) {
    
    double pt = 30.; // value not really relevant for L1RC
    double eta = hrc->GetBinCenter(i);

    jecrc->setJetPt(pt);
    jecrc->setJetEta(eta);
    jecrc->setJetA(TMath::Pi()*0.4*0.4);
    jecrc->setRho(rho);
    double corrrc = jecrc->getCorrection(); // = (pt-off)/pt
    double offrc = pt - pt*corrrc;
    hrc->SetBinContent(i, 2*offrc);

    double ptl1 = 2*offrc;
    jecl1->setJetPt(ptl1);
    jecl1->setJetEta(eta);
    jecl1->setJetA(TMath::Pi()*0.4*0.4);
    jecl1->setRho(rho);
    double corrl1 = jecl1->getCorrection(); // = (pt-off)/pt
    double offl1 = ptl1 - ptl1*corrl1;
    hl1->SetBinContent(i, 2*offl1);

  } // for i

  TH1D *h = new TH1D("h",";Jet #eta;2 #times offset (GeV)",neta,veta);
  h->SetMaximum(30.);

  TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);
  tdrDraw(hrc,"HIST");
  tdrDraw(hl1,"P");

  gPad->RedrawAxis();

  c1->SaveAs("pdf/estimateMinJetPt.pdf");
}
