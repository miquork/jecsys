// Draw ratio of gamma+jet with two different EGamma corrections
#include "../tdrstyle_mod14.C"
#include "../tools.h"
#define STANDALONE
#include "../CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "../CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "../CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "../ErrorTypes.hpp"
#include "../JetDefs.hpp"

#include "TFile.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TLine.h"
#include "TF1.h"
#include "TH1D.h"
#include "TMultiGraph.h"


TH1D* getHistoFromTxtFile(string label, const char* txtFile){
  TH1D *histo = new TH1D(label.c_str(),";#eta;chi2/ndf",52,0,5.2);

  cout << txtFile << endl << flush;
  JetCorrectorParameters *l2l3res = new JetCorrectorParameters(txtFile);
  vector<JetCorrectorParameters> v;
  v.push_back(*l2l3res);
  FactorizedJetCorrector* jec = new FactorizedJetCorrector(v);

  jec->setJetA(0.1);
  jec->setRho(20);

  for (unsigned binx=1;binx<=histo->GetNbinsX()+1;binx++) {
    jec->setJetPt(20);
    jec->setJetEta(histo->GetBinCenter(binx));
    double corr = jec->getCorrection();
    histo->SetBinContent(binx,corr);
  }
  return histo;
}

void drawL2L3ResChi2Variants() {
  
  setTDRStyle();

  // Save background histogram for easy drawing from root file
  TH1D *h = new TH1D("h",";#eta;chi2/ndf",
		     52,0,5.2);
  h->SetMinimum(0.00);
  h->SetMaximum(15);

  //  string directory = "CondFormats/JetMETObjects/data/";
  string directory = "EGM2_2par/";
  const char *d = directory.c_str();
  const char *a = "AK4PFchs";

  
  map<string, const char*> txtFileMap;
  //  txtFileMap["1par_MJDJ_gam_zll_PtBalMPF"]=Form("%sCollectL2Output_1par_MJDJ_gam_zll_PtBalMPF/Summer16LegacyBCDEFGH_VXXX_DATA_L2L3Residual_Chi2OverNDF_%s.txt",d,a);
  txtFileMap["2par_MJDJ_gam_zll_PtBalMPF"]=Form("%sCollectL2Output_MJDJ_gam_zll_PtBalMPF/Summer16LegacyBCDEFGH_VXXX_DATA_L2L3Residual_Chi2OverNDF_%s.txt",d,a);
  txtFileMap["2par_MJDJ_gam_zll_MPF"]=Form("%sCollectL2Output_MJDJ_gam_zll_MPF/Summer16LegacyBCDEFGH_VXXX_DATA_L2L3Residual_Chi2OverNDF_%s.txt",d,a);
  txtFileMap["2par_MJDJ_gam_zll_PtBal"]=Form("%sCollectL2Output_MJDJ_gam_zll_PtBal/Summer16LegacyBCDEFGH_VXXX_DATA_L2L3Residual_Chi2OverNDF_%s.txt",d,a);
  txtFileMap["2par_gam_zll_PtBalMPF"]=Form("%sCollectL2Output_gam_zll_PtBalMPF/Summer16LegacyBCDEFGH_VXXX_DATA_L2L3Residual_Chi2OverNDF_%s.txt",d,a);
  txtFileMap["2par_DJ_PtBalMPF"]=Form("%sCollectL2Output_DJ_PtBalMPF/Summer16LegacyBCDEFGH_VXXX_DATA_L2L3Residual_Chi2OverNDF_%s.txt",d,a);
  txtFileMap["2par_zll_PtBalMPF"]=Form("%sCollectL2Output_zll_PtBalMPF/Summer16LegacyBCDEFGH_VXXX_DATA_L2L3Residual_Chi2OverNDF_%s.txt",d,a);
  txtFileMap["2par_gam_PtBalMPF"]=Form("%sCollectL2Output_gam_PtBalMPF/Summer16LegacyBCDEFGH_VXXX_DATA_L2L3Residual_Chi2OverNDF_%s.txt",d,a);
  
  int colors[8]  = {kBlack,kBlue, kGreen+2, kRed, kGray, kBlue-10, kGreen+2-10,kRed-10};
  int markerstyle[8]  = {kOpenCircle, kFullCircle, kOpenTriangleDown, kFullTriangleDown, kOpenTriangleUp, kFullTriangleUp, kOpenSquare, kFullSquare};

  TCanvas *c1 = tdrCanvas("c1",h);
  TLegend *legdw = tdrLeg(0.15,0.70,0.50,0.90);
  h->Draw();
  for (auto it=txtFileMap.begin(); it!=txtFileMap.end(); ++it){
    int idx = std::distance(txtFileMap.begin(),it);
    //    std::cout << std::distance(aMap.begin(), it_map) << '\n';
    std::cout << idx << ": " << it->first << " => " << it->second << '\n';
    TH1D* histo = getHistoFromTxtFile(it->first,it->second);
    tdrDraw(histo,"PL",markerstyle[idx],colors[idx],kSolid,-1,0);
    legdw->AddEntry(histo,it->first.c_str(),"PL");

  }
  legdw->Draw();
  c1->SaveAs("pdf/GlobalFitChi2_Summary.pdf");

}
