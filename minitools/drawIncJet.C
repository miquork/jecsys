// Purpose: draw inclusive jets at |y|<1.3 in 2016 paper style
// with option to flip it vs NLO+NLL-to-NNLO interpolation
// (assuming NLO+NLL for gluon jets, NNLO for quark jets)
#include "TFile.h"
#include "TF1.h"

#include "../tdrstyle_mod15.C"

void drawIncJet() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  // Load data
  TFile *f = new TFile("rootfiles/commonUL2017_V4_V2M5res_hotzone.root","READ");
  assert(f && !f->IsZombie());
  curdir->cd();

  TH1D *hd = (TH1D*)f->Get("ak4/Eta_0.0-1.3/hpt_data_2017_all_ptcl_fwd");
  assert(hd);
  TH1D *hd1 = (TH1D*)f->Get("ak4/Eta_0.0-0.5/hpt_data_2017_all_ptcl_fwd");
  assert(hd1);
  TH1D *hd2 = (TH1D*)f->Get("ak4/Eta_0.5-1.0/hpt_data_2017_all_ptcl_fwd");
  assert(hd2);
  TH1D *hd3 = (TH1D*)f->Get("ak4/Eta_1.0-1.5/hpt_data_2017_all_ptcl_fwd");
  assert(hd3);
  TH1D *hstat = (TH1D*)hd->Clone("hstat");
  hstat->Divide(hd);

  // Load JEC uncertainty
  TFile *fe = new TFile("rootfiles/jecdataBCDEF.root","READ");
  assert(f && !f->IsZombie());
  curdir->cd();  

  TH1D *he = (TH1D*)fe->Get("ratio/eta00-13/herr_ref"); assert(he);
  he = (TH1D*)he->Clone("he");
  for (int i = 1; i != he->GetNbinsX()+1; ++i) {
    he->SetBinContent(i,1);
    he->SetBinError(i, he->GetBinError(i)*5);
  }

  // Load theory predictions
  TFile *f1 = new TFile("rootfiles/cmsJetsNLO_AK4.root","READ");
  assert(f1 && !f1->IsZombie());
  TFile *f2 = new TFile("rootfiles/Ratio_Theory/Ratio_NLO_Fig5b.root","READ");
  assert(f2 && !f2->IsZombie());
  TFile *f3 = new TFile("rootfiles/np_ew.root","READ");
  assert(f3 && !f3->IsZombie());
  curdir->cd();

  // NLO prediction
  TH1D *hnlo1 = (TH1D*)f1->Get("histCT14nlo_PDF_Cnt_y0");
  assert(hnlo1);
  TH1D *hnlo2 = (TH1D*)f1->Get("histCT14nlo_PDF_Cnt_y1");
  assert(hnlo2);
  TH1D *hnlo3 = (TH1D*)f1->Get("histCT14nlo_PDF_Cnt_y2");
  assert(hnlo3);
  TH1D *hnlo = (TH1D*)hnlo3->Clone("hnlo");
  TH1D *href = (TH1D*)hnlo3->Clone("href");

  // Corrections for R-dependence, (+NP, EW). Only for |y|<0.5 but flat vs y
  TH1D *h8lo = (TH1D*)f2->Get("Cross_section_ratio_ak8_by_ak4_ybin1_LO");
  TH1D *h8nlo = (TH1D*)f2->Get("Cross_section_ratio_ak8_by_ak4_ybin1_NLO");
  TH1D *hr = (TH1D*)h8lo->Clone("hr");
  hr->Divide(h8nlo);

  // Smooth out R-dependence, NLO has low stats and up to 2% fluctuations
  TF1 *fr = new TF1("fr","[p0]+log(0.01*x)*([p1]+[p2]*log(0.01*x))",5,7000);
  fr->SetParameters(8.45577e-01, 5.51566e-02, -7.10071e-03);
  hr->Fit(fr,"QRNW");
      
  // EW corrections, NP corrections
  TH1D *hew1 = (TH1D*)f3->Get("ew16_ak4_y0");
  assert(hew1);
  TH1D *hew2 = (TH1D*)f3->Get("ew16_ak4_y1");
  assert(hew2);
  TH1D *hew3 = (TH1D*)f3->Get("ew16_ak4_y2");
  assert(hew3);

  // NP corrections
  TH1D *hnp1 = (TH1D*)f3->Get("np16_ak4_y0");
  assert(hnp1);
  TH1D *hnp2 = (TH1D*)f3->Get("np16_ak4_y1");
  assert(hnp2);
  TH1D *hnp3 = (TH1D*)f3->Get("np16_ak4_y2");
  assert(hnp3);

  // NLL kFactor
  TH1D *hnll1 = (TH1D*)f3->Get("kFactorNLL_ak4_y0");
  assert(hnll1);
  TH1D *hnll2 = (TH1D*)f3->Get("kFactorNLL_ak4_y1");
  assert(hnll2);
  TH1D *hnll3 = (TH1D*)f3->Get("kFactorNLL_ak4_y2");
  assert(hnll3);
  TH1D *hnll = (TH1D*)hnll3->Clone("hnll");

  // NNLO kFactor
  TH1D *hnnlo1 = (TH1D*)f3->Get("kFactorNNLO_ak4_y0");
  assert(hnnlo1);
  TH1D *hnnlo2 = (TH1D*)f3->Get("kFactorNNLO_ak4_y1");
  assert(hnnlo2);
  TH1D *hnnlo3 = (TH1D*)f3->Get("kFactorNNLO_ak4_y2");
  assert(hnnlo3);
  TH1D *hnnlo = (TH1D*)hnnlo3->Clone("hnnlo");

      
  for (int i = 1; i != hnlo3->GetNbinsX()+1; ++i) {

    double x = hnlo3->GetBinCenter(i);
    double dy = 0.6;

    double yd(0), eyd2(0);
    double y(0), ey2(0);
    int j1 = hd->FindBin(hnlo3->GetBinLowEdge(i)+1e-1);
    int j2 = hd->FindBin(hnlo3->GetBinLowEdge(i+1)-1e-1);
    for (int j = j1; j != j2+1; ++j) {
      yd += hd->GetBinContent(j)*hd->GetBinWidth(j);
      eyd2 += pow(hd->GetBinError(j)*hd->GetBinWidth(j),2);
      y += (1.0*hd1->GetBinContent(j)*hd1->GetBinWidth(j) +
	    1.0*hd2->GetBinContent(j)*hd2->GetBinWidth(j) +
	    dy*hd3->GetBinContent(j)*hd3->GetBinWidth(j)) / (2+dy);
      ey2 += (1.0*pow(hd1->GetBinError(j)*hd->GetBinWidth(j),2) +
	      1.0*pow(hd1->GetBinError(j)*hd->GetBinWidth(j),2) +
	      dy*pow(hd1->GetBinError(j)*hd->GetBinWidth(j),2)) / (2+dy);
    }
    yd /= (hd->GetBinLowEdge(j2+1)-hd->GetBinLowEdge(j1));
    double eyd = sqrt(eyd2) / (hd->GetBinLowEdge(j2+1)-hd->GetBinLowEdge(j1));
    y /= (hd1->GetBinLowEdge(j2+1)-hd1->GetBinLowEdge(j1));
    double ey = sqrt(ey2) / (hd1->GetBinLowEdge(j2+1)-hd1->GetBinLowEdge(j1));

    double xrh = min(max(x,hr->GetBinLowEdge(1)+1e-3),
		    hr->GetBinLowEdge(hr->GetNbinsX()+1));

    double xew = min(max(x,hew3->GetBinLowEdge(1)),
		     hew3->GetBinLowEdge(hew3->GetNbinsX()+1));

    // First bin set to 1.0 for AK4 NPF
    double xnp = min(max(x,hnp3->GetBinLowEdge(2)+1e-2),
		     hnp3->GetBinLowEdge(hnp3->GetNbinsX()+1));

    // Only [97,2238] GeV available for kfactor y2
    double xnll = min(max(x,hnll3->GetBinLowEdge(3)+1e-2),
		      hnll3->GetBinLowEdge(22));

    // Only [97,2238] GeV available for kfactor y2
    double xnnlo = min(max(x,hnnlo3->GetBinLowEdge(3)+1e-2),
		       hnnlo3->GetBinLowEdge(22));
    
    double krh = hr->GetBinContent(hr->FindBin(xrh)); // (hist)
    double krf = fr->Eval(x); // radius correction wrt AK8 (fit)
    double kew = (1.0*hew1->GetBinContent(hew1->FindBin(xew)) +
		  1.0*hew2->GetBinContent(hew2->FindBin(xew)) +
		  dy*hew3->GetBinContent(hew3->FindBin(xew))) / (2+dy);
    double knp = (1.0*hnp1->GetBinContent(hnp1->FindBin(xnp)) +
		  1.0*hnp2->GetBinContent(hnp2->FindBin(xnp)) +
		  dy*hnp3->GetBinContent(hnp3->FindBin(xnp))) / (2+dy);
    double knll = (1.0*hnll1->GetBinContent(hnll1->FindBin(xnll)) +
		   1.0*hnll2->GetBinContent(hnll2->FindBin(xnll)) +
		   dy*hnll3->GetBinContent(hnll3->FindBin(xnll))) / (2+dy);

    double knnlo = (1.0*hnnlo1->GetBinContent(hnnlo1->FindBin(xnnlo)) +
		    1.0*hnnlo2->GetBinContent(hnnlo2->FindBin(xnnlo)) +
		    dy*hnnlo3->GetBinContent(hnnlo3->FindBin(xnnlo))) / (2+dy);

    double nlo = (1.0*hnlo1->GetBinContent(hnlo1->FindBin(x)) +
		  1.0*hnlo2->GetBinContent(hnlo2->FindBin(x)) +
		  dy*hnlo3->GetBinContent(hnlo3->FindBin(x))) /  (2+dy);

    /*
    // With NP
    href->SetBinContent(i,  nlo * krf * kew * knp / yd);
    hnlo->SetBinContent(i,  nlo * kew * knp / yd);
    hnll->SetBinContent(i,  nlo * knll * kew * knp / yd);
    hnnlo->SetBinContent(i,  nlo * knnlo * kew * knp / yd);
    */
    // Without NP
    href->SetBinContent(i,  nlo * krf * kew / y);
    hnlo->SetBinContent(i,  nlo * kew / y);
    hnll->SetBinContent(i,  nlo * knll * kew  / y);
    hnnlo->SetBinContent(i,  nlo * knnlo * kew / y);
  } // for i

  double ptmin = 15;
  double ptmax = 3500;
  TH1D *h = tdrHist("h","Theory / data",0.92,1.22,//0.7,1.5,
		    "Jet p_{T} (GeV)",ptmin,ptmax);
  lumi_13TeV = "UL2017, 41.5 fb^{-1}";
  TCanvas *c1 = tdrCanvas("c1",h,4,11);//,kSquare);
  gPad->SetLogx();
  //gPad->SetLogy();

  tdrDraw(he,"LE3",kNone,kYellow+2,kSolid,-1,1001,kYellow+1);
  tdrDraw(hstat,"Pz",kFullCircle,kBlack);
  tdrDraw(href,"H",kNone,kBlack,kSolid,-1,kNone);
  tdrDraw(hnlo,"H",kNone,kRed,kSolid,-1,kNone);
  tdrDraw(hnll,"H",kNone,kBlue,kSolid,-1,kNone);
  tdrDraw(hnnlo,"H",kNone,kMagenta+1,kSolid,-1,kNone);

  hstat->SetMarkerSize(0.5);
  hnll->SetLineWidth(3);
  hnnlo->SetLineWidth(3);

  TLegend *leg = tdrLeg(0.65,0.63,0.85,0.88);
  leg->AddEntry(hstat,"Data (|y|<1.3'ish)","PLE");
  leg->AddEntry(hnlo,"NLO","L");
  leg->AddEntry(hnnlo,"NNLO","L");
  leg->AddEntry(href,"NLO #times RF8","L");
  leg->AddEntry(hnll,"NLO #times NLL","L");

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  tex->DrawLatex(0.5,0.17,"All theories #times NP #times EW");

  gPad->RedrawAxis();

  c1->SaveAs("pdf/drawIncJet.pdf");
} // drawIncJet
