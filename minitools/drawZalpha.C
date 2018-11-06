#include "../tdrstyle_mod15.C"
#include "../tools.C"

#include "TFile.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TLine.h"

void drawZalpha() {

  setTDRStyle();

  //TFile *f = new TFile("../rootfiles/jecdataBCDEFGH.root","READ");
  TFile *f = new TFile("../../jecsys/rootfiles/jecdataBCDEFGH.root","READ");
  assert(f && !f->IsZombie());

  // eta00-13
  // zmmjet (could also use zlljet)
  TGraphErrors *gd30 = (TGraphErrors*)f->Get("data/eta00-13/ptchs_zmmjet_a30");
  TGraphErrors *gd20 = (TGraphErrors*)f->Get("data/eta00-13/ptchs_zmmjet_a20");
  TGraphErrors *gd15 = (TGraphErrors*)f->Get("data/eta00-13/ptchs_zmmjet_a15");
  TGraphErrors *gd10 = (TGraphErrors*)f->Get("data/eta00-13/ptchs_zmmjet_a10");
  assert(gd30);
  assert(gd20);
  assert(gd15);
  assert(gd10);

  TGraphErrors *gd3 = (TGraphErrors*)f->Get("data/eta00-13/mpfchs1_zmmjet_a30");
  TGraphErrors *gd2 = (TGraphErrors*)f->Get("data/eta00-13/mpfchs1_zmmjet_a20");
  TGraphErrors *gd1 = (TGraphErrors*)f->Get("data/eta00-13/mpfchs1_zmmjet_a15");
  TGraphErrors *gd0 = (TGraphErrors*)f->Get("data/eta00-13/mpfchs1_zmmjet_a10");
  assert(gd3);
  assert(gd2);
  assert(gd1);
  assert(gd0);

  TGraphErrors *gm30 = (TGraphErrors*)f->Get("mc/eta00-13/ptchs_zmmjet_a30");
  TGraphErrors *gm20 = (TGraphErrors*)f->Get("mc/eta00-13/ptchs_zmmjet_a20");
  TGraphErrors *gm15 = (TGraphErrors*)f->Get("mc/eta00-13/ptchs_zmmjet_a15");
  TGraphErrors *gm10 = (TGraphErrors*)f->Get("mc/eta00-13/ptchs_zmmjet_a10");
  assert(gm30);
  assert(gm20);
  assert(gm15);
  assert(gm10);

  TGraphErrors *gm3 = (TGraphErrors*)f->Get("mc/eta00-13/mpfchs1_zmmjet_a30");
  TGraphErrors *gm2 = (TGraphErrors*)f->Get("mc/eta00-13/mpfchs1_zmmjet_a20");
  TGraphErrors *gm1 = (TGraphErrors*)f->Get("mc/eta00-13/mpfchs1_zmmjet_a15");
  TGraphErrors *gm0 = (TGraphErrors*)f->Get("mc/eta00-13/mpfchs1_zmmjet_a10");
  assert(gm3);
  assert(gm2);
  assert(gm1);
  assert(gm0);

  //int ipt = 2;//3;//0;//3;//3; // jecsys2016/
  int ipt = 3; // jecsys/
  // pT~200 GeV is first point where alpha_max=0.10 is (almost) fully efficient
  double pt = gd30->GetX()[ipt];
  TGraphErrors *gd = new TGraphErrors(0);
  tools::SetPoint(gd, 0, 0.10, gd10->GetY()[ipt], 0, gd10->GetEY()[ipt]);
  tools::SetPoint(gd, 1, 0.15, gd15->GetY()[ipt], 0, gd15->GetEY()[ipt]);
  tools::SetPoint(gd, 2, 0.20, gd20->GetY()[ipt], 0, gd20->GetEY()[ipt]);
  tools::SetPoint(gd, 3, 0.30, gd30->GetY()[ipt], 0, gd30->GetEY()[ipt]);
  TGraphErrors *gm = new TGraphErrors(0);
  tools::SetPoint(gm, 0, 0.10, gm10->GetY()[ipt], 0, gm10->GetEY()[ipt]);
  tools::SetPoint(gm, 1, 0.15, gm15->GetY()[ipt], 0, gm15->GetEY()[ipt]);
  tools::SetPoint(gm, 2, 0.20, gm20->GetY()[ipt], 0, gm20->GetEY()[ipt]);
  tools::SetPoint(gm, 3, 0.30, gm30->GetY()[ipt], 0, gm30->GetEY()[ipt]);

  int impf = ipt+4;
  double mpf = gd3->GetX()[impf];
  TGraphErrors *gdm = new TGraphErrors(0);
  tools::SetPoint(gdm, 0, 0.10, gd0->GetY()[impf], 0, gd0->GetEY()[impf]);
  tools::SetPoint(gdm, 1, 0.15, gd1->GetY()[impf], 0, gd1->GetEY()[impf]);
  tools::SetPoint(gdm, 2, 0.20, gd2->GetY()[impf], 0, gd2->GetEY()[impf]);
  tools::SetPoint(gdm, 3, 0.30, gd3->GetY()[impf], 0, gd3->GetEY()[impf]);
  TGraphErrors *gmm = new TGraphErrors(0);
  tools::SetPoint(gmm, 0, 0.10, gm0->GetY()[impf], 0, gm0->GetEY()[impf]);
  tools::SetPoint(gmm, 1, 0.15, gm1->GetY()[impf], 0, gm1->GetEY()[impf]);
  tools::SetPoint(gmm, 2, 0.20, gm2->GetY()[impf], 0, gm2->GetEY()[impf]);
  tools::SetPoint(gmm, 3, 0.30, gm3->GetY()[impf], 0, gm3->GetEY()[impf]);

  // Scale to %
  for (int i = 0; i != gd->GetN(); ++i) {
    tools::SetPoint(gd, i, gd->GetX()[i], 100.*(gd->GetY()[i]-1),
		    0, gd->GetEY()[i]*100);
    tools::SetPoint(gm, i, gm->GetX()[i], 100.*(gm->GetY()[i]-1),
		    0, gm->GetEY()[i]*100);
    tools::SetPoint(gdm, i, gdm->GetX()[i], 100.*(gdm->GetY()[i]-1),
		    0, gdm->GetEY()[i]*100);
    tools::SetPoint(gmm, i, gmm->GetX()[i], 100.*(gmm->GetY()[i]-1),
		    0, gmm->GetEY()[i]*100);
  }

  double ptmax = 1200.*2.;
  TH1D *hup = new TH1D("hup",";#alpha_{max};#LTp_{T,jet}#GT / p_{T,Z}^{ } - 1 (%)",
		       10,0,0.34);
  hup->SetMinimum(-10.5);//-26);//-10.5);
  hup->SetMaximum(+3);//+7);//+3);

  TH1D *hdw = new TH1D("hdw",";#alpha_{max};Data/MC-1 (%)",10,0,0.34);
  hdw->SetMinimum(-0.4);//-10);//-0.4);
  hdw->SetMaximum(+0.6);//0);//+0.6);
  hdw->GetXaxis()->SetMoreLogLabels();
  hdw->GetXaxis()->SetNoExponent();

  //lumi_13TeV = "Run2016BCDEFGH 36.5 fb^{-1}";
  lumi_13TeV = "Run2016 36.5 fb^{-1}";
  extraText = "Private Work";
  TCanvas *c1 = tdrDiCanvas("c1",hup,hdw,4,11);
  hup->GetYaxis()->SetTitleOffset(1.00);

  c1->cd(1);
  tdrDraw(gdm,"Pz",kFullSquare,kBlue);//-9);
  tdrDraw(gmm,"Pz",kOpenSquare,kBlue);//-9);
  tdrDraw(gd,"Pz",kFullCircle,kRed);//-9);
  tdrDraw(gm,"Pz",kOpenCircle,kRed);//-9);

  TF1 *fdm = new TF1("fdm","[0]+[1]*x",0,0.35);
  fdm->FixParameter(1,0);
  fdm->SetLineColor(kBlue-0);
  fdm->SetLineStyle(kSolid);
  gdm->Fit(fdm,"QRN");
  fdm->Draw("SAME");

  TF1 *fmm = new TF1("fmm","[0]+[1]*x",0,0.35);
  fmm->FixParameter(1,0);
  gmm->Fit(fmm,"QRN");
  fmm->SetLineColor(kBlue);//-9);
  fmm->SetLineStyle(kDashed);
  fmm->Draw("SAME");

  // Scale data by MPF data/MC ratio
  double  k = (fmm->GetParameter(0)/100.+1) / (fdm->GetParameter(0)/100.+1);
  for (int i = 0; i != gd->GetN(); ++i) {
    gd->SetPoint(i, gd->GetX()[i], ((gd->GetY()[i]/100.+1)*k-1)*100.);
    gdm->SetPoint(i, gdm->GetX()[i], ((gdm->GetY()[i]/100.+1)*k-1)*100.);
  }
  gdm->Fit(fdm,"QRN");


  TF1 *fm = new TF1("fm","[0]+[1]*x",0,0.35);
  double ue = TMath::Pi()*0.4*0.4*2.00;//1.52*log(13)/log(8); // Z+jet UE
  //double ue = TMath::Pi()*0.4*0.4*2.00*log(13)/log(8); // dijet UE
  //double ue = TMath::Pi()*0.4*0.4*2.00; // dijet UE??
  fm->FixParameter(0, ue/pt * 100.);
  //fm->SetParameter(0, ue/pt * 100.);
  gm->Fit(fm,"QRN");
  ue = fm->GetParameter(0)*pt/100.;
  double due = fm->GetParError(0)*pt/100.;
  gm->SetPoint(gm->GetN(), 0, ue/pt * 100.); // UE/pT,Z-1
  //gmm->SetPoint(gmm->GetN(), 0, 1.50*ue/pt * 100.); // R_UE*UE/pT,Z-1
  fm->SetLineColor(kRed);//-9);
  fm->SetLineStyle(kDashed);
  fm->Draw("SAME");

  TF1 *fd = new TF1("fd","[0]+[1]*x",0,0.35);
  fd->FixParameter(0, ue/pt * 100.);
  fd->SetLineColor(kRed);//-9);
  //gd->SetPoint(gd->GetN(),0,fd->GetParameter(0)); // UE/pT,Z-1 (before fit)
  fd->SetLineStyle(kSolid);
  gd->Fit(fd,"QRN");
  gd->SetPoint(gd->GetN(),0,fd->GetParameter(0)); // UE/pT,Z-1 (after fit)
  fd->Draw("SAME");


  TLatex *tex = new TLatex();
  tex->SetTextSize(0.045);
  tex->SetNDC();

  tex->DrawLatex(0.20,0.11,"|#eta| < 1.3, R = 0.4 PF+CHS");
  //tex->DrawLatex(0.20,0.11,"2.853 < |#eta|< 2.964, R = 0.4 PF+CHS");
  //tex->DrawLatex(0.20,0.05,Form("#LTp_{T}#GT = %1.0f (%1.0f) GeV",pt,mpf));
  tex->DrawLatex(0.20,0.05,Form("#LTp_{T}#GT = %1.0f GeV",pt));//,mpf));

  tex->SetTextColor(kBlue);//-9);
  tex->DrawLatex(0.50,0.70,Form("#chi^{2} / NDF = %1.1f / %d (%1.1f / %d)",
				fdm->GetChisquare(), fdm->GetNDF(),
				fmm->GetChisquare(), fmm->GetNDF()));

  tex->SetTextColor(kRed);//-9);
  tex->DrawLatex(0.20,0.17,Form("#chi^{2} / NDF = %1.1f / %d (%1.1f / %d)",
				fd->GetChisquare(), fd->GetNDF(),
				fm->GetChisquare(), fm->GetNDF()));
  tex->DrawLatex(0.20,0.24,Form("k_{FSR} = %1.3f #pm %1.3f",
				//tex->DrawLatex(0.20,0.17,Form("k_{FSR} = %1.3f #pm %1.3f",
				fd->GetParameter(1)*0.01,
				fd->GetParError(1)*0.01));
  //tex->DrawLatex(0.20,0.30,Form("#rho_{UE} = %1.2f #pm %1.2f GeV/A",
  tex->DrawLatex(0.20,0.30,Form("#rho_{UE} = %1.2f GeV/A",
				ue/(TMath::Pi()*0.4*0.4)));
  //			ue/(TMath::Pi()*0.4*0.4),
  //			due/(TMath::Pi()*0.4*0.4)));
  
  gPad->RedrawAxis();
  
  TLegend *leg = tdrLeg(0.65,0.35,0.85,0.65);
  leg->SetHeader("Z(#rightarrow#mu#mu)+jet");
  //leg->SetHeader("Z(#rightarrowl^{+}l^{-})+jet");
  leg->AddEntry(gdm,"Data MPF","PL");
  leg->AddEntry(gmm,"MC MPF","PL");
  leg->AddEntry(gd,"Data p_{T}^{bal}","PL");
  leg->AddEntry(gm,"MC p_{T}^{bal}","PL");

  c1->cd(2);

  TGraphErrors *grm = (TGraphErrors*)gdm->Clone();
  for (int i = 0; i != grm->GetN(); ++i) {
    tools::SetPoint(grm, i, gdm->GetX()[i],
		    100*((gdm->GetY()[i]/100.+1)/(gmm->GetY()[i]/100.+1)-1),
		    0, sqrt(pow(gdm->GetEY()[i],2)+pow(gmm->GetEY()[i],2)));
  }

  TGraphErrors *gr = (TGraphErrors*)gd->Clone();
  for (int i = 0; i != gr->GetN(); ++i) {
    tools::SetPoint(gr, i, gd->GetX()[i],
		    100*((gd->GetY()[i]/100.+1)/(gm->GetY()[i]/100.+1)-1),
		    0, sqrt(pow(gd->GetEY()[i],2)+pow(gm->GetEY()[i],2)));
  }

  tdrDraw(grm,"Pz",kFullSquare,kBlue);//-9);
  tdrDraw(gr,"Pz",kFullCircle,kRed);//-9);

  TF1 *frm = new TF1("frm","((([0]+[1]*x)/100.+1)"
		     " / (([2]+[3]*x)/100.+1) - 1)*100",0,0.35);
  frm->SetParameters(fdm->GetParameter(0), fdm->GetParameter(1),
		     fmm->GetParameter(0), fmm->GetParameter(1));
  frm->SetLineColor(kBlue);//-9);
  frm->Draw("SAME");

  TF1 *fr = new TF1("fr","((([0]+[1]*x)/100.+1)"
		    " / (([2]+[3]*x)/100.+1) - 1)*100",0,0.35);
  fr->SetParameters(fd->GetParameter(0), fd->GetParameter(1),
		    fm->GetParameter(0), fm->GetParameter(1));
  fr->SetLineColor(kRed);//-9);
  fr->Draw("SAME");

  gPad->RedrawAxis();

  c1->Update();
  c1->SaveAs("../pdf/drawZalpha.pdf");
}
