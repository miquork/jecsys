// Purpose: plot top pT variants to derive and/or check reweighing
#include "TFile.h"
#include "TLine.h"
#include "TF1.h"

#include "../tdrstyle_mod15.C"

TF1 *_f1(0);
Double_t der1(Double_t *x, Double_t *p) {
  assert(_f1);
  return (_f1->Derivative(*x));
}

void hadW_toppt() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TFile *fd = new TFile("rootfiles/hadWUL1718V5_MPDGcorrNoW.root","READ");
  assert(fd && !fd->IsZombie());
  TFile *fm = new TFile("rootfiles/hadWMC1718V5_MPDGcorrNoW.root","READ");
  assert(fm && !fm->IsZombie());

  curdir->cd();

  TH1D *htd = (TH1D*)fd->Get("htoppt"); assert(htd);
  TH1D *htm = (TH1D*)fm->Get("htoppt"); assert(htm);
  TH1D *htm_fsr = (TH1D*)fm->Get("htoppt_fsr"); assert(htm_fsr);
  TH1D *htm_q2qg = (TH1D*)fm->Get("htoppt_q2qg"); assert(htm_q2qg);
  TH1D *htm_x2xg = (TH1D*)fm->Get("htoppt_x2xg"); assert(htm_x2xg);
  TH1D *htm_isr = (TH1D*)fm->Get("htoppt_isr"); assert(htm_isr);
  TH1D *htg = (TH1D*)fm->Get("htopptgen"); assert(htg);

  htd->Scale(1./htd->Integral());
  double kmc = 1./htm->Integral();
  htm->Scale(kmc);
  htm_fsr->Scale(kmc); //1./htm_fsr->Integral());
  htm_q2qg->Scale(kmc); //1./htm_q2qg->Integral());
  htm_x2xg->Scale(kmc); //1./htm_x2xg->Integral());
  htm_isr->Scale(kmc); //1./htm_isr->Integral());
  htg->Scale(kmc); //1./htg->Integral());

  TH1D *hu = tdrHist("hu","Fraction of events per bin",1.5e-3,1.5e-1,
		     "Top p_{T} (GeV)",0,400);
  TH1D *hd = tdrHist("hd","Data/MC",0.7,1.2,
		     "Top p_{T} (GeV)",0,400);

  TH1D *htr = (TH1D*)htd->Clone("htr");
  htr->Divide(htm);
  TH1D *htr_fsr = (TH1D*)htm_fsr->Clone("htr_fsr");
  htr_fsr->Divide(htm);
  TH1D *htr_q2qg = (TH1D*)htm_q2qg->Clone("htr_q2qg");
  htr_q2qg->Divide(htm);
  TH1D *htr_x2xg = (TH1D*)htm_x2xg->Clone("htr_x2xg");
  htr_x2xg->Divide(htm);
  TH1D *htr_isr = (TH1D*)htm_isr->Clone("htr_isr");
  htr_isr->Divide(htm);
  TH1D *htr_g = (TH1D*)htg->Clone("htr_g");
  htr_g->Divide(htm);

  lumi_13TeV = "UL17+18, 101.4 fb^{-1}";
  TCanvas *c1 = tdrDiCanvas("c1",hu,hd,4,11);

  c1->cd(1);
  gPad->SetLogx();
  gPad->SetLogy();

  tdrDraw(htg,"HPz",kOpenDiamond,kBlack,kSolid,-1,kNone);
  tdrDraw(htm,"HPz",kOpenCircle,kBlack,kSolid,-1,kNone);
  tdrDraw(htm_q2qg,"HPz",kOpenTriangleUp,kBlue,kSolid,-1,kNone);
  tdrDraw(htm_x2xg,"HPz",kOpenTriangleUp,kRed,kSolid,-1,kNone);
  tdrDraw(htm_fsr,"HPz",kOpenTriangleUp,kBlack,kSolid,-1,kNone);
  tdrDraw(htm_isr,"HPz",kOpenTriangleDown,kBlack,kSolid,-1,kNone);
  tdrDraw(htd,"HPz",kFullCircle,kBlack,kSolid,-1,kNone);

  gPad->RedrawAxis();

  c1->cd(2);
  gPad->SetLogx();

  tdrDraw(htr_g,"HPz",kOpenDiamond,kBlack,kSolid,-1,kNone);
  tdrDraw(htr,"HPz",kFullCircle,kBlack,kSolid,-1,kNone);
  tdrDraw(htr_q2qg,"HPz",kOpenTriangleUp,kBlue,kSolid,-1,kNone);
  tdrDraw(htr_x2xg,"HPz",kOpenTriangleUp,kRed,kSolid,-1,kNone);
  tdrDraw(htr_fsr,"HPz",kOpenTriangleUp,kBlack,kSolid,-1,kNone);
  tdrDraw(htr_isr,"HPz",kOpenTriangleDown,kBlack,kSolid,-1,kNone);

  htr->SetMarkerSize(0.7);

  TLine *l = new TLine();
  l->SetLineStyle(kDotted);
  l->DrawLine(1,1,400,1);

  TF1 *f1pt = new TF1("f1pt","[0]+[1]*pow(x,[2])",1,400);
  f1pt->SetParameters(1,-0.1,0.1);

  htr->Fit(f1pt,"RN");
  f1pt->Draw("SAME");

  cout << "// Data/MC ratio from minitools/hadW_toppt.C"<< endl;
  cout << Form("// Chisquare / NDF = %1.1f / %d (%1.3f)",
	       f1pt->GetChisquare(), f1pt->GetNDF(),
	       f1pt->GetChisquare() / f1pt->GetNDF()) << endl;
  cout << Form("TF1 *f1pt = new TF1(\"f1pt\",\"%s\",1,400);",
	       f1pt->GetExpFormula().Data()) << endl;
  cout << Form("f1pt->SetParameters(%1.4g,%1.4g,%1.4g);",
	       f1pt->GetParameter(0),f1pt->GetParameter(1),
	       f1pt->GetParameter(2)) << endl;

  c1->SaveAs("pdf/hadW_toppt.pdf");


  // Why is top pT distribution in data such a pure log-normal
  // with peak value exactly at mW?
  TCanvas *c2 = tdrCanvas("c2",hu,4,11,kSquare);
  gPad->SetLogx();
  gPad->SetLogy();

  tdrDraw(htd,"HPz",kFullCircle,kBlack,kSolid,-1,kNone);

  TF1 *f1 = new TF1("f1","[0]*TMath::LogNormal(x,[1],[2],[3])",
		    //1,400); // chi2/NDF poor
		    //10,180); // fairly good fit
		    40,140); // good fit, mode 80.46

  //f1->SetParameters(0.06,80.4,0,1);
  f1->SetParameters(10,0.46,-4.2,151);
  //f1->FixParameter(2,0); // not so good fit
  htd->Fit(f1,"RN");
  f1->Draw("SAME");

  // Show extended fit range
  TF1 *f1ext = new TF1("f1ext","[0]*TMath::LogNormal(x,[1],[2],[3])",1,400);
  f1ext->SetParameters(f1->GetParameter(0),f1->GetParameter(1),
		       f1->GetParameter(2),f1->GetParameter(3));
  f1ext->SetLineStyle(kDotted);
  f1ext->DrawClone("SAME");
  htd->Fit(f1ext,"RN");
  f1ext->SetLineStyle(kDashed);
  f1ext->Draw("SAME");

  cout << "// Top pT fit in data from minitools/hadW_toppt.C"<< endl;
  cout << Form("// Chisquare / NDF = %1.1f / %d (%1.3f)",
	       f1->GetChisquare(), f1->GetNDF(),
	       f1->GetChisquare() / f1->GetNDF()) << endl;
  cout << Form("TF1 *f1 = new TF1(\"f1\",\"%s\",1,400);",
	       f1->GetExpFormula().Data()) << endl;
  cout << Form("f1->SetParameters(%1.4g,%1.4g,%1.4g,%1.4g);",
	       f1->GetParameter(0),f1->GetParameter(1),
	       f1->GetParameter(2),f1->GetParameter(3)) << endl;

  cout << Form("// Chisquare / NDF = %1.1f / %d (%1.3f)",
	       f1ext->GetChisquare(), f1ext->GetNDF(),
	       f1ext->GetChisquare() / f1ext->GetNDF()) << endl;
  cout << Form("TF1 *f1ext = new TF1(\"f1ext\",\"%s\",1,400);",
	       f1ext->GetExpFormula().Data()) << endl;
  cout << Form("f1ext->SetParameters(%1.4g,%1.4g,%1.4g,%1.4g);",
	       f1ext->GetParameter(0),f1ext->GetParameter(1),
	       f1ext->GetParameter(2),f1ext->GetParameter(3)) << endl;

  // Double_t LogNormal(Double_t x, Double_t sigma, Double_t theta = 0, Double_t m = 1)
  // lognormal_pdf uses same definition of http://en.wikipedia.org/wiki/Log-normal_distribution
  // where mu = log(m)
  // return ::ROOT::Math::lognormal_pdf(x, TMath::Log(m), sigma, theta);
  // But this has
  // 1/(x*sqrt(2pi s^2)) * exp(-(log(x)-m)^2/(2 s^2))
  //double mu = TMath::Log(f1->GetParameter(3));
  //double sigma = f1->GetParameter(1);
  //double mode = exp(mu - sigma*sigma);
  //double mode = 1./exp(sigma*sigma);

  TF1 *df1 = new TF1("df1",der1,0,1,400);
  _f1 = f1;
  double mode = df1->GetX(0.,60,100);

  cout << Form("mode = %1.4g ",mode) << endl;

  c2->SaveAs("pdf/hadW_toppt_LogNormal.pdf");
} // hadW_toppt
