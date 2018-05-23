// File: drawJetFlavor.C
// Created by Mikko Voutilainen, on May 22nd, 2018
// Purpose: Draw the flavor fractions of Z+jet, gamma+jet and dijet samples,
//          which are added into the global fit combination file by jetflavor.C

#include "TFile.h"
//#include "TGraphErrors.h"
#include "TH1D.h"
//#include "TProfile.h"
//#include "TCanvas.h"
//#include "TLegend.h"
//#include "TF1.h"
//#include "TMultiGraph.h"
#include "TLatex.h"
//#include "TMinuit.h"
//#include "TMatrixD.h"
//#include "TVectorD.h"
//#include "TMath.h"

//#include "../tools.h"
#include "tdrstyle_mod14.C"

#include <string>
//#include <iostream>
//#include <fstream>
//#include <vector>
#include <map>

using namespace std;

bool samex = true;//false;

void drawJetFlavor(double etamin=0.0, double etamax=1.3,
		   //string sample="zlljet",  
		   string sample="gamjet",  
		   //string sample="dijet",  
		   string epoch="BCDEFGH") {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  const int nflavors = 4;
  string flavors[nflavors] = {"uds","g","c","b"};
  map<string, string> labels;
  labels["uds"] = "Quark";
  labels["c"] = "Charm";
  labels["b"] = "Bottom";
  labels["g"] = "Gluon";

  double kx = 13000./8000.;
  double xmin(30.), xmax(400.);
  if (sample=="zlljet") { xmin=30.; xmax=400.; }
  if (sample=="gamjet") { xmin=30.; xmax=800.; }
  if (sample=="dijet")  { xmin=50.; xmax=2000.; }
  if (samex) { xmin *= kx; xmax *= kx; }

  TH1D *h = new TH1D("h",";p_{T,Z} (GeV);Flavor fraction", 100, xmin, xmax);
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  h->SetMaximum(1.2-1e-4);

  lumi_13TeV = "";
  extraText = "Simulation";
  TCanvas *c1 = tdrCanvas("c1",h,4,0,kSquare);
  c1->SetLogx();

  TLegend *leg = tdrLeg(0.20,0.65,0.50,0.90);
  leg->SetHeader("ME parton flavor");

  TLatex *tex = new TLatex();
  tex->SetTextSize(0.045);
  tex->SetNDC();
  if (sample=="zlljet") {
    tex->DrawLatex(0.47,0.81,"Z+jet, |#eta| < 1.3, #alpha < 0.3");
    tex->DrawLatex(0.47,0.76,"MG+Pythia8 CUETP8M1");
    h->SetXTitle("p_{T,Z} (GeV)");
  }
  if (sample=="gamjet") {
    tex->DrawLatex(0.47,0.81,"#gamma+jet, |#eta| < 1.3, #alpha < 0.3");
    tex->DrawLatex(0.47,0.76,"Pythia8 CUETP8M1");
    h->SetXTitle("p_{T,#gamma} (GeV)");
  }
  if (sample=="dijet") {
    tex->DrawLatex(0.47,0.81,"Dijet, |#eta| < 1.3, #alpha < 0.3");
    tex->DrawLatex(0.47,0.76,"Pythia8 CUETP8M1");
    h->SetXTitle("p_{T,jet} (GeV)");
  }

  TFile *f = new TFile("rootfiles/jecdataBCDEFGH.root","READ");
  assert(f && !f->IsZombie());

  f->cd("ratio");
  gDirectory->cd("eta00-13");
  gDirectory->cd("flv");
  TDirectory *fd = gDirectory;

  for (int i = 0 ; i != nflavors; ++i) {

    string sf = Form("%s_%s",sample.c_str(),flavors[i].c_str());
    cout << sf << endl;
    TH1D *hf = (TH1D*)fd->Get(sf.c_str());
    assert(hf);

    //tdrDraw(hf,"P",hf->GetLineColor());
    hf->GetXaxis()->SetRangeUser(xmin, xmax);
    hf->Draw("SAMEP");
    leg->AddEntry(hf,labels[flavors[i]].c_str(),"PL");
  }

  c1->SaveAs(Form("pdf/drawJetFlavor_%s%s.pdf",
		  sample.c_str(),samex ? "_samexrange" : ""));
}
