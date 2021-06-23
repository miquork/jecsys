// Purpose: Redraw B hadron semileptonic branching ratios, add PDG 2020 b buark
// Original code: https://gitlab.cern.ch/CMS-TOPPAG/BFragmentationAnalyzer/-/blob/master/test/drawBrs.py
// Original plot: https://gitlab.cern.ch/CMS-TOPPAG/BFragmentationAnalyzer/-/blob/master/data/semilepbr_unc.png
// LHCb results: https://lhcbproject.web.cern.ch/lhcbproject/Publications/LHCbProjectPublic/Directory_LHCb-PAPER-2018-050/Table_1.pdf
// at https://lhcbproject.web.cern.ch/lhcbproject/Publications/LHCbProjectPublic/LHCb-PAPER-2018-050.html
#include "TGraphErrors.h"
#include "TBox.h"
#include "TLatex.h"
#include "TFile.h"
#include "TProfile.h"

#include "../tdrstyle_mod15.C"

void drawBRs() {

  setTDRStyle();

  // Input data for B0, B+, B0s and LambdaB BR(X->lnuX)
  const int nb = 4;
  double ahfdt[nb] = {0.400,  0.400,  0.12,  0.08}; 
  double ahfp8[nb] = {0.425,  0.425,  0.11,  0.04}; 
  double abrp8[nb] = {0.104, 0.113, 0.093, 0.077};
  double abr16[nb] = {0.104, 0.110, 0.096, 0.1025};
  double ebr16[nb] = {0.003, 0.003, 0.008, 0.002};

  // New values from TOP 2020 graph
  double abrc5[nb] = {0.1040, 0.1127, 0.093, 0.077};
  double abr20[nb] = {0.1035, 0.1100, 0.096, 0.109}; // LambdaB mean change
  double ebr20[nb] = {0.0027, 0.0027, 0.008, 0.022}; // LambdaB error x10

  // New values from LHCb paper
  double alhcb[nb+1] = {0.1030, 0.1108, 0.1024, 0.1026, 0.1070};
  double elhcb[nb+1] = {0.0019, 0.0020, 0.0021, 0.0025, 0.0019};

  // Data from TOP mt l+jet analysis
  //TFile *f = new TFile("rootfiles/hadWMC17V5New_Glu.root","READ");
  //TFile *f = new TFile("rootfiles/hadWMC17V5New_JEC.root","READ");
  //TFile *f = new TFile("rootfiles/hadWMC17V5New_JECV2.root","READ");
  //TFile *f = new TFile("rootfiles/hadWMC17V5New_SLBRWeight.root","READ");
  TFile *f = new TFile("rootfiles/hadWMC17V5New_SLBRWeightV3.root","READ");
  assert(f && !f->IsZombie());
  TProfile *pbr = (TProfile*)f->Get("pbhslbr"); assert(pbr);
  TProfile *pfr = (TProfile*)f->Get("pbhfrac"); assert(pfr);
  int aidx[nb] = {2,1,3,4};

  // NB: could high muon content in Z+b and mt analyses be missing
  // b-tag SF of ~0.90 and mistag rates above unity? Bkgs have few mu's?

  // Create graphs
  TGraphErrors *gp8 = new TGraphErrors(nb);
  TGraphErrors *g16 = new TGraphErrors(nb);
  TGraphErrors *geb = new TGraphErrors(nb);
  TGraphErrors *gc5 = new TGraphErrors(nb);
  TGraphErrors *g20 = new TGraphErrors(nb);
  TGraphErrors *g20a = new TGraphErrors(nb);
  TGraphErrors *glhcb = new TGraphErrors(nb);
  TGraphErrors *gmt = new TGraphErrors(nb);
  double brdt(0), brd8(0), brp8(0), brd5(0), brc5(0);
  double brdmt(0), brmt(0), br20(0), br20e(0), brlhcb(0), brlhcbe(0);
  for (int i = 0; i != nb; ++i) {
    
    brdt += ahfdt[i]*abr20[i];
    brd8 += ahfdt[i]*abrp8[i];
    brp8 += ahfp8[i]*abrp8[i];
    brd5 += ahfdt[i]*abrc5[i];
    brc5 += ahfp8[i]*abrc5[i];
    brdmt += ahfdt[i]*0.5*pbr->GetBinContent(aidx[i]);
    brmt  += pfr->GetBinContent(aidx[i])/pfr->Integral(1,4)
      *0.5*pbr->GetBinContent(aidx[i]);
    br20  += ahfdt[i]*abr20[i];
    br20e += ahfdt[i]*ebr20[i];
    brlhcb  += ahfdt[i]*alhcb[i];
    brlhcbe += ahfdt[i]*elhcb[i];

    cout << " BR = " << 0.5*pbr->GetBinContent(aidx[i]) << endl;

    // Set points
    gp8->SetPoint(i, abrp8[i], i+1+0.25);
    gp8->SetPointError(i, 0., 0.);
    g16->SetPoint(i, abr16[i], i+1-0.25);
    g16->SetPointError(i, ebr16[i], 0.);
    gc5->SetPoint(i, abrc5[i], i+1+0.25);//+0.125);
    gc5->SetPointError(i, 0., 0.);
    g20->SetPoint(i, abr20[i], i+1-0.125);
    g20->SetPointError(i, ebr20[i], 0.);
    glhcb->SetPoint(i, alhcb[i], i+1-0.25);
    glhcb->SetPointError(i, elhcb[i], 0.);
    gmt->SetPoint(i, 0.5*pbr->GetBinContent(aidx[i]), i+1+0.125);
    gmt->SetPointError(i, 0.5*pbr->GetBinError(aidx[i]), 0.);

    // Calculate envelope
    //double eup = max(abrp8[i], abr16[i]+ebr16[i]);
    //double edw = min(abrp8[i], abr16[i]-ebr16[i]);
    double eup = max(abrc5[i], abr20[i]+ebr20[i]);
    double edw = min(abrc5[i], abr20[i]-ebr20[i]);
    double mid = 0.5*(eup+edw);
    double err = 0.5*(eup-edw);

    geb->SetPoint(i, mid, i+1);
    geb->SetPointError(i, err, 0.5);

  } // for i
  gp8->SetPoint(nb, brp8, 0.25);
  gp8->SetPointError(nb, 2*0.5*fabs(brd8-brp8), 0.);
  //gc5->SetPoint(nb, brc5, 0.25+0.125);
  gc5->SetPoint(nb, brc5, 0.40);
  gc5->SetPointError(nb, 2*0.5*fabs(brd5-brc5), 0.);
  //glhcb->SetPoint(nb, alhcb[nb], 0.10);
  //glhcb->SetPointError(nb, elhcb[nb], 0.);
  glhcb->SetPoint(nb, brlhcb, 0.10);
  glhcb->SetPointError(nb, brlhcbe, 0.);
  gmt->SetPoint(nb, brmt, 0.30);
  gmt->SetPointError(nb, 2*0.5*(brdmt-brmt), 0.);
  g20->SetPoint(nb, br20, 0.20);
  g20->SetPointError(nb, br20e, 0.);

  cout << "SL BR = " << brmt << endl;

  TH1D *h = new TH1D("h",";BR(b#rightarrowl#nuX);",120,0.05,0.16);
  lumi_13TeV = "";
  TCanvas *c1 = tdrCanvas("c1",h,4,0,kSquare);
  h->SetMinimum(0.);//0.5);
  h->SetMaximum(nb+0.5);

  //TBox *Abox = new TBox(0.1069-0.0035,0.5,0.1069+0.0035,nb+0.5);
  TBox *Abox = new TBox(0.1069-0.0035,0.0,0.1069+0.0035,nb+0.5);
  Abox->SetFillColor(kGreen-9);
  Abox->Draw("SAME");

  TLine *l = new TLine();
  l->SetLineColor(kGreen+2);
  //l->DrawLine(0.1069,0.5,0.1069,nb+0.5);
  l->DrawLine(0.1069,0.0,0.1069,nb+0.5);

  for (int i = 0; i != nb; ++i) {
    TBox *Bbox = new TBox(geb->GetX()[i]-geb->GetEX()[i],
			  geb->GetY()[i]-geb->GetEY()[i],
			  geb->GetX()[i]+geb->GetEX()[i],
			  geb->GetY()[i]+geb->GetEY()[i]);
    Bbox->SetFillColor(kGray);
    Bbox->SetFillColorAlpha(kGray,0.7);
    Bbox->Draw("SAME");
  }

  TLatex *tex = new TLatex();
  tex->SetTextSize(0.045);
  /*
  tex->DrawLatex(0.135,0.25,"All");
  tex->DrawLatex(0.135,1.,"B^{0}");
  tex->DrawLatex(0.135,2.,"B^{+}");
  tex->DrawLatex(0.135,3.,"B^{0}_{s}");
  tex->DrawLatex(0.135,4.,"#Lambda_{b}");
  */
  tex->DrawLatex(0.132,0.25,Form("All (%1.3f)",brlhcb/gmt->GetX()[nb]));
  tex->DrawLatex(0.132,1.,Form("B^{0} (%1.3f)",alhcb[0]/gmt->GetX()[0]));
  tex->DrawLatex(0.132,2.,Form("B^{+} (%1.3f)",alhcb[1]/gmt->GetX()[1]));
  tex->DrawLatex(0.132,3.,Form("B^{0}_{s} (%1.3f)",alhcb[2]/gmt->GetX()[2]));
  tex->DrawLatex(0.132,4.,Form("#Lambda_{b} (%1.3f)",alhcb[3]/gmt->GetX()[3]));

  g20a->SetFillColor(kGreen-9);
  g20a->SetLineColor(kGreen+2);
  geb->SetFillColor(kGray);
  geb->SetLineColor(kGray);

  //tdrDraw(geb,"E3",kNone,kGray,kSolid,-1,1001,kGray);//->Draw("SAME E3");
  //tdrDraw(g16,"P",kFullCircle,kBlack);//->Draw("SAME P");
  tdrDraw(glhcb,"P",kFullCircle,kBlack);//->Draw("SAME P");
  //tdrDraw(gp8,"P",kOpenCircle,kBlack);//->Draw("SAME P");
  tdrDraw(g20,"P",kFullCircle,kBlue);//->Draw("SAME P");
  //tdrDraw(gc5,"P",kOpenCircle,kBlue);//->Draw("SAME P");
  tdrDraw(gc5,"P",kOpenCircle,kBlack);//->Draw("SAME P");
  tdrDraw(gmt,"P",kOpenCircle,kBlue);//->Draw("SAME P");

  TLegend *leg = tdrLeg(0.18,0.52-0.06*6,0.38,0.52);
  leg->AddEntry(g20a,"PDG 2020 avg.","LF");
  //leg->AddEntry(geb,"TOP syst. (2016)","F");
  leg->AddEntry(geb,"TOP syst.","F");
  leg->AddEntry(g20,"PDG 2020","PE");
  //leg->AddEntry(g16,"PDG 2016","PE");
  leg->AddEntry(glhcb,"LHCb 2019","PE");
  //leg->AddEntry(gc5,"Pythia8 (2020)","P");
  //leg->AddEntry(gp8,"Pythia8 (2018)","P");
  leg->AddEntry(gc5,"Pythia8","P");
  leg->AddEntry(gmt,"Pythia8 (m_{t})","P");

  gPad->RedrawAxis();

  c1->SaveAs("pdf/drawBRs.pdf");
} // drawBRs
