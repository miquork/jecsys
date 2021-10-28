// Purpose: derive JER C-term scale factor from modified hotzone maps
#include "TFile.h"
#include "TProfile2D.h"

#include "../tdrstyle_mod15.C"

void jerCterm() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  // Open input file
  TFile *f = new TFile("rootfiles/output-DATA-1-Ctest.root","READ");
  assert(f && !f->IsZombie());
  curdir->cd();

  // Retrieve input 2D profile
  TProfile2D *p2d = (TProfile2D*)f->Get("FullEta_Reco/jt0/p2djmpf");
  assert(p2d);
  TH2D *h2d = p2d->ProjectionXY("h2dfmpf");

  // Setup plotting for 2D profile
  TH1D *h0 = new TH1D("h0",";#eta;#phi",120,-5.236,5.236);
  h0->GetYaxis()->SetRangeUser(-TMath::Pi(),TMath::Pi());
  h2d->GetZaxis()->SetRangeUser(-1,+1);
  lumi_13TeV = "ZeroBias 10k - 2D profile of MPF";
  TCanvas *c0 = tdrCanvas("c0",h0,4,11,kSquare);
  h2d->Draw("SAME COLZ");
  gPad->SetRightMargin(0.15);
  gPad->RedrawAxis();

  // Setup plotting for eta strip scan
  TH1D *h = new TH1D("h",";#phi;MPF mean for probe vs tag",
		     72,-TMath::Pi(),TMath::Pi());
  h->GetYaxis()->SetRangeUser(-1,1.5);
  lumi_13TeV = "ZeroBias 10k - #eta scan";
  TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);
  
  // Loop over barrel eta bins
  int ieta1 = h2d->GetXaxis()->FindBin(-1.3+1e-4);
  int ieta2 = h2d->GetXaxis()->FindBin(+1.3-1e-4);
  TH1D *hrms = new TH1D("hrms",";Mean of probe vs tag;Entries",250,-1,1.5);
  TH1D *hrmsw = new TH1D("hrmsw",";Mean of probe vs tag;Weights",250,-1,1.5);
  for (int ieta = ieta1; ieta != ieta2+1; ++ieta) {

    // Project eta bin to phi
    TH1D *h1 = h2d->ProjectionY(Form("h1_%d",ieta),ieta,ieta);
    tdrDraw(h1,"HISTE",kNone,ieta-ieta1+1,kSolid,-1,kNone);

    // Fill scatter of means (non-empty bins only for now)
    for (int j = 1; j != h1->GetNbinsX()+1; ++j) {
      if (h1->GetBinError(j)!=0) {
	hrms->Fill(h1->GetBinCenter(j));
	hrmsw->Fill(h1->GetBinCenter(j), 1./pow(h1->GetBinError(j),2));
      }
    }
  }

  // Setup plotting for RMS
  TH1D *h2 = new TH1D("h",";MPF mean for probe vs tag;Normalized entries",
		     72,-TMath::Pi(),TMath::Pi());
  h2->GetYaxis()->SetRangeUser(0,0.1);
  lumi_13TeV = "ZeroBias 10k - RMS distribution";
  TCanvas *c2 = tdrCanvas("c2",h2,4,11,kSquare);

  hrms->Scale(1./hrms->Integral());
  hrms->Draw("SAME");

  hrms->SetLineColor(kRed);
  hrms->SetMarkerColor(kRed);
  hrmsw->Scale(1./hrmsw->Integral());
  hrmsw->Draw("SAME");

} // jerCterm
