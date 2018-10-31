// File: drawUncEtaVsPt.C
// Created by Henning Kirschenmann, on 5 October 2018
// Purpose: Grab uncertainty txt-files and produce 
//          2D uncertanty map as .root and .pdf 
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
#include <iostream>
#include <vector>
#include <map>

using namespace std;


// put all the different methods in a single file for easy access for everybody
void drawUncEtaVsPt(string JEC="") {

  // Set TDR style to have correct graphical setttings when storing graphs
  setTDRStyle();
 
  TDirectory *curdir = gDirectory;

  const char *s, *s2;
  // Total uncertainty
  s = Form("%s",JEC.c_str());
  s2 = "Total";
  cout << s << ":" << s2 << endl << flush;
  JetCorrectorParameters *p_unc = new JetCorrectorParameters(s,s2);
  JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p_unc);



  const double x_pt[] =
  //    {8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
    {15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
     97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 362, 430,
     507, 592, 686, 790, 905, 1032, 1172, 1327, 1497, 1684, 1890, //1999};
     //2000, 2238, 2500, 2787, 3103, 3450,
     2116, 2366, 2640, 2941, 3273, 3637, 
     4037, 4477, 4961, 5492, 6076, 7000};
  const int ndiv_pt = sizeof(x_pt)/sizeof(x_pt[0])-1;
  
  const double x_eta[] =
    {-5.4,-5.0,-4.4,-4,-3.5,-3,-2.8,-2.6,-2.4,-2.2,-2.0,
     -1.8,-1.6,-1.4,-1.2,-1.0, -0.8,-0.6,-0.4,-0.2,0.,
     0.2,0.4,0.6,0.8,1.0,1.2,1.4,
     1.6,1.8,2.0,2.2,2.4,2.6,2.8,3,3.5,4,4.4,5.0,5.4};
  const int ndiv_eta = sizeof(x_eta)/sizeof(x_eta[0])-1;     

  // Uncertainty bands
  TH2D *herr = new TH2D("herr",";p_{T} (GeV);#eta;JEC uncertainty",
                        ndiv_pt, &x_pt[0],ndiv_eta,&x_eta[0]);


  for (int i = 1; i != herr->GetNbinsX()+1; ++i) {
    for (int j = 1; j != herr->GetNbinsY()+1; ++j) {
      double pt = herr->GetXaxis()->GetBinCenter(i);
      double eta = herr->GetYaxis()->GetBinCenter(j);

      //      cout << "eta" << eta << " pt: " << pt  << endl;

      // JEC uncertainties
      unc->setJetEta(eta);
      unc->setJetPt(pt);
      double err = unc->getUncertainty(true);
      //      cout << "eta" << eta << " pt: " << pt  << " err: " << err<< endl;
      herr->SetBinContent(i,j,err);
      herr->SetBinError(i,j,err);
    }
  }



  TCanvas *c1 = tdrCanvas("c1",(TH1D*)herr,4,11,kSquare);
  //c1->SetLogx();


  herr->Draw("colz");
  c1->SetRightMargin( 0.22 );
  herr->GetZaxis()->SetTitleOffset(1.4);
  CMS_lumi( c1, 4, 11 );  
  c1->RedrawAxis();
  // For some reason calling SetLogx earlier drops result precision?
  c1->SetLogx();
  
  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  tex->DrawLatex(0.20,0.75,"Anti-k_{T} R=0.4");


  
  c1->SaveAs("IdealizedUncertainty.pdf");
  TFile *fout = new TFile(Form("outHistos.root"), "RECREATE");
  herr->Write();
  fout->Close();

} // reprocess



