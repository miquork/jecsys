// Draw ratio of Zee+jet and Zmm+jet
#include "../tdrstyle_mod15.C"
#include "../tools.C"

#include "TFile.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TLine.h"
#include "TF1.h"
#include "TMultiGraph.h"
#include "TMatrixD.h"

void drawZeeVsZmm(string run = "BCDEF") {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TFile *f = new TFile(Form("../rootfiles/jecdata%s.root",run.c_str()),"READ");
  assert(f && !f->IsZombie());
  gDirectory->cd("ratio");
  gDirectory->cd("eta00-13");
  TDirectory *d = gDirectory;
  TGraphErrors *gem = (TGraphErrors*)d->Get("mpfchs1_zeejet_a30"); assert(gem);
  TGraphErrors *gep = (TGraphErrors*)d->Get("ptchs_zeejet_a30"); assert(gep);
  TGraphErrors *gmm = (TGraphErrors*)d->Get("mpfchs1_zmmjet_a30"); assert(gmm);
  TGraphErrors *gmp = (TGraphErrors*)d->Get("ptchs_zmmjet_a30"); assert(gmp);

  TGraphErrors *gpm = (TGraphErrors*)d->Get("mpfchs1_gamjet_a30"); assert(gpm);
  TGraphErrors *gpp = (TGraphErrors*)d->Get("ptchs_gamjet_a30"); assert(gpp);

  TH1D *hpm = (TH1D*)d->Get("fsr/hkfsr_mpfchs1_gamjet"); assert(hpm);
  TH1D *hpp = (TH1D*)d->Get("fsr/hkfsr_ptchs_gamjet"); assert(hpp);
  for (int i = 0; i != gpm->GetN(); ++i) {
    double kfsr  = hpm->GetBinContent(hpm->FindBin(gpm->GetX()[i]));
    gpm->SetPoint(i, gpm->GetX()[i], gpm->GetY()[i]*(1. - 0.3*kfsr));
  }
  for (int i = 0; i != gpp->GetN(); ++i) {
    double kfsr  = hpp->GetBinContent(hpp->FindBin(gpp->GetX()[i]));
    gpp->SetPoint(i, gpp->GetX()[i], gpp->GetY()[i]*(1. - 0.3*kfsr));
  }


  curdir->cd();

  double ptmax = 2.*1200.;
  TH1D *hup = new TH1D("hup",";p_{T,Z} (GeV);Data / MC",670,30,ptmax);
  hup->SetMinimum(0.96);
  hup->SetMaximum(1.01);
  hup->GetXaxis()->SetMoreLogLabels();
  hup->GetXaxis()->SetNoExponent();

  TH1D *hdw = new TH1D("hdw",";p_{T} (GeV);Zee / Z#mu#mu - 1 (%)",670,30,ptmax);
  hdw->SetMinimum((0.970+1e-5-1)*100);
  hdw->SetMaximum((1.030+1e-5-1)*100);
  hdw->GetXaxis()->SetMoreLogLabels();
  hdw->GetXaxis()->SetNoExponent();

  map<string, const char*> lumimap;
  lumimap["B"] = "Run2017B, 4.8 fb^{-1}";
  lumimap["C"] = "Run2017C, 9.6 fb^{-1}";
  lumimap["D"] = "Run2017D, 4.2 fb^{-1}";
  lumimap["E"] = "Run2017E, 9.3 fb^{-1}";
  lumimap["F"] = "Run2017F, 13.4 fb^{-1}";
  lumimap["BCDEF"] = "Run20167BCDEF, 41.4 fb^{-1}";
  lumi_13TeV = lumimap[run];
  TCanvas *c1 = tdrDiCanvas("c1",hdw,hup,4,11);

  c1->cd(2);
  gPad->SetLogx();

  tdrDraw(gmm,"Pz",kFullSquare,kRed);
  tdrDraw(gmp,"Pz",kOpenSquare,kRed);

  tdrDraw(gem,"Pz",kFullCircle,kGreen+2);
  tdrDraw(gep,"Pz",kOpenCircle,kGreen+2);

  tdrDraw(gpm,"Pz",kFullDiamond,kBlue);
  // Don't draw gam+jet pTbal, because different MC => kFSR does not cancel
  tdrDraw(gpp,"Pz",kOpenDiamond,kBlue);

  c1->cd(1);
  gPad->SetLogx();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(30,0,ptmax,0);

  // Ratio new/old

  TGraphErrors *grm = tools::ratioGraphs(gem,gmm);
  TGraphErrors *grp = tools::ratioGraphs(gep,gmp);

  TGraphErrors *grpm = tools::ratioGraphs(gpm,gmm);
  TGraphErrors *grpp = tools::ratioGraphs(gpp,gmp);

  // Change result to (Zee/Zmm-1)*100
  for (int i = 0; i != grm->GetN(); ++i) {
    grm->SetPoint(i, grm->GetX()[i], (grm->GetY()[i]-1)*100);
    grm->SetPointError(i, grm->GetEX()[i], (grm->GetEY()[i])*100);
  }
  for (int i = 0; i != grp->GetN(); ++i) {
    grp->SetPoint(i, grp->GetX()[i], (grp->GetY()[i]-1)*100);
    grp->SetPointError(i, grp->GetEX()[i], (grp->GetEY()[i])*100);
  }
  // Change result to (gam/Zmm-1)*100
  for (int i = 0; i != grpm->GetN(); ++i) {
    grpm->SetPoint(i, grpm->GetX()[i], (grpm->GetY()[i]-1)*100);
    grpm->SetPointError(i, grpm->GetEX()[i], (grpm->GetEY()[i])*100);
  }
  for (int i = 0; i != grpp->GetN(); ++i) {
    grpp->SetPoint(i, grpp->GetX()[i], (grpp->GetY()[i]-1)*100);
    grpp->SetPointError(i, grpp->GetEX()[i], (grpp->GetEY()[i])*100);
  }



  TLegend *legdw = tdrLeg(0.5,0.60,0.70,0.90);
  legdw->AddEntry(grm,"MPF","PL");
  legdw->AddEntry(grp,"p_{T} balance","PL");
  legdw->AddEntry(grpm,"#gamma / Z#mu#mu MPF","PL");
  legdw->AddEntry(grpp,"#gamma / Z#mu#mu p_{T}^{bal}","PL");

  tdrDraw(grm,"Pz",kFullCircle,kGreen+2);
  tdrDraw(grp,"Pz",kOpenSquare,kGreen+2);

  tdrDraw(grpm,"Pz",kFullDiamond,kBlue);
  // Don't draw gam+jet pTbal, because different MC => kFSR does not cancel
  tdrDraw(grpp,"Pz",kOpenDiamond,kBlue);

  TMultiGraph *gm = new TMultiGraph();
  gm->Add(grm);
  gm->Add(grp);

  TF1 *fs = new TF1("fs","([0]+[1]*log(0.01*x)+[2]*pow(log(0.01*x),2)-1)*100",
		    30,ptmax);
  fs->SetParameters(1,0,0);//1.7,-1,1);
  //fs->FixParameter(2,0);

  gm->Fit(fs,"QRN");
  fs->SetLineStyle(kSolid);
  fs->SetLineColor(kGreen+2);
  fs->Draw("SAME");

  TMatrixD emat(fs->GetNpar(),fs->GetNpar());
  gMinuit->mnemat(&emat[0][0],fs->GetNpar());

  // d/dp0=1, d/dp1=log(x), d/dp2=log(x)^2
  // eps =   (d/dp0)^2 m_00 +   (d/dp1)^2 m_11   +   (d/dp2)^2 m_22
  //     + 2*(d/dp0/p1)m_01 + 2*(d/dp0/p2)m_02 + 2*(d/dp1/p2)m_12
  TF1 *fse = new TF1("fse","([0]+[1]*log(0.01*x)+[2]*pow(log(0.01*x),2)"
		     "+ [3]*sqrt([4] + pow(log(0.01*x),2)*[5]"
		     "+ pow(log(0.01*x),4)*[6]"
		     "+ 2*log(0.01*x)*[7]+2*pow(log(0.01*x),2)*[8]"
		     "+ 2*pow(log(0.01*x),3)*[9])"
		     " - 1)*100",30,ptmax);
  fse->SetParameters(fs->GetParameter(0),fs->GetParameter(1),
		     fs->GetParameter(2), +1,
		     emat[0][0], emat[1][1], emat[2][2],
		     emat[0][1], emat[0][2], emat[1][2]);

  fse->SetLineColor(kGreen-8);
  fse->DrawClone("SAME");
  fse->SetParameter(3, -1);
  fse->DrawClone("SAME");
  
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);
  tex->DrawLatex(0.20,0.72,"#alpha < 0.30");
  tex->DrawLatex(0.20,0.67,"Data / MC ratio");
  //
  tex->DrawLatex(0.20,0.20,Form("#chi^{2}/NDF = %1.1f/%d",
				fs->GetChisquare(), fs->GetNDF()));
  tex->DrawLatex(0.20,0.15,Form("p_{0} = %1.2f #pm %1.2f",
				100*(fs->GetParameter(0)-1),
				100*fs->GetParError(0)));
  tex->DrawLatex(0.20,0.10,Form("p_{1} = %1.2f #pm %1.2f",
				100.*fs->GetParameter(1),
				100.*fs->GetParError(1)));
  tex->DrawLatex(0.20,0.05,Form("p_{2} = %1.2f #pm %1.2f",
				100.*fs->GetParameter(2),
				100.*fs->GetParError(2)));


  // Compare to Zee mass corrections
  TF1 *f1mzee = new TF1("f1mzee","([3]/([0]+[1]*log(0.01*x)"
			"+[2]*pow(log(0.01*x),2))"
			"-1)*100", 30, ptmax);
  // BCDEFGH fit with minitools/drawZmass.C
  //f1mzee->SetParameters(0.99885, 0.00176, 0.00135, 0.99868);//EGM1
  //f1mzee->SetParameters(1.00017, 0.00166, 0.00114, 0.99868);//EGM2
  //f1mzee->SetParameters(1.00279, 0.00166, 0.00112, 0.99868);//EGM3
  f1mzee->SetParameters(1.00246, 0.00214, 0.00116, 0.99854); // Nov17V10BCDEF

  // Zee mass applied to gamma+jet at pT,Z=2*pT,gamma
  TF1 *f1mgam = new TF1("f1mgam","([3]/([0]+[1]*log(0.01*x*2)"
			"+[2]*pow(log(0.01*x*2),2))"
			"-1)*100", 30, ptmax);
  f1mgam->SetParameters(1.00246, 0.00214, 0.00116, 0.99854); // Nov17V10BCDEF
      
  c1->cd(1);
  f1mzee->DrawClone("SAME");
  f1mzee->SetLineStyle(kSolid);

  f1mgam->SetLineStyle(kDashed);
  f1mgam->DrawClone("SAME");

  legdw->AddEntry(f1mzee,"#mu#mu/ee mass","L");
  legdw->AddEntry(f1mgam,"#mu#mu/#gamma mass","L");

  gPad->Update();

  c1->SaveAs(Form("../pdf/drawZeeVsZmm_Run%s.pdf",run.c_str()));
}
