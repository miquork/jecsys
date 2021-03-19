// Draw ratio of Zee+jet and Zmm+jet
// run with 'root -l -b -q minitools/drawZeeVsZmm.C+g'
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

void drawZeeVsZmms(string run = "BCDEFGH") {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TFile *f = new TFile(Form("rootfiles/jecdata%s.root",run.c_str()),"READ");
  assert(f && !f->IsZombie());
  gDirectory->cd("ratio");
  gDirectory->cd("eta00-13");
  TDirectory *d = gDirectory;

  TFile *ffsr = new TFile(Form("rootfiles/jecdata%s.root",
			       run=="2018ABCD" ? "2018ABCD" :
			       run=="2016GH" ? "2016GH" : "BCDEF"),"READ");
  assert(ffsr && !ffsr->IsZombie());
  gDirectory->cd("ratio");
  gDirectory->cd("eta00-13");
  TDirectory *dfsr = gDirectory;

  /*
  TGraphErrors *gem = (TGraphErrors*)d->Get("mpfchs1_zeejet_a30"); assert(gem);
  TGraphErrors *gep = (TGraphErrors*)d->Get("ptchs_zeejet_a30"); assert(gep);
  TGraphErrors *gmm = (TGraphErrors*)d->Get("mpfchs1_zmmjet_a30"); assert(gmm);
  TGraphErrors *gmp = (TGraphErrors*)d->Get("ptchs_zmmjet_a30"); assert(gmp);

  TGraphErrors *gpm = (TGraphErrors*)d->Get("mpfchs1_gamjet_a30"); assert(gpm);
  TGraphErrors *gpp = (TGraphErrors*)d->Get("ptchs_gamjet_a30"); assert(gpp);
  */
  TGraphErrors *gem = (TGraphErrors*)d->Get("mpfchs1_zeejet_a100"); assert(gem);
  TGraphErrors *gep = (TGraphErrors*)d->Get("ptchs_zeejet_a100"); assert(gep);
  TH1D *heh = (TH1D*)d->Get("hdm_mpfchs1_zeejet"); assert(heh);
  TGraphErrors *gmm = (TGraphErrors*)d->Get("mpfchs1_zmmjet_a100"); assert(gmm);
  TGraphErrors *gmp = (TGraphErrors*)d->Get("ptchs_zmmjet_a100"); assert(gmp);
  TH1D *hmh = (TH1D*)d->Get("hdm_mpfchs1_zmmjet"); assert(hmh);

  TGraphErrors *gpm = (TGraphErrors*)d->Get("mpfchs1_gamjet_a100"); assert(gpm);
  TGraphErrors *gpp = (TGraphErrors*)d->Get("ptchs_gamjet_a100"); assert(gpp);

  //TH1D *hpm = (TH1D*)d->Get("fsr/hkfsr_mpfchs1_gamjet"); assert(hpm);
  //TH1D *hpp = (TH1D*)d->Get("fsr/hkfsr_ptchs_gamjet"); assert(hpp);
  TH1D *hpm = (TH1D*)dfsr->Get("fsr/hkfsr_mpfchs1_gamjet"); assert(hpm);
  TH1D *hpp = (TH1D*)dfsr->Get("fsr/hkfsr_ptchs_gamjet"); assert(hpp);
  for (int i = 0; i != gpm->GetN(); ++i) {
    double kfsr  = hpm->GetBinContent(hpm->FindBin(gpm->GetX()[i]));
    gpm->SetPoint(i, gpm->GetX()[i], gpm->GetY()[i]*(1. - 0.3*kfsr));
  }
  for (int i = 0; i != gpp->GetN(); ++i) {
    double kfsr  = hpp->GetBinContent(hpp->FindBin(gpp->GetX()[i]));
    gpp->SetPoint(i, gpp->GetX()[i], gpp->GetY()[i]*(1. - 0.3*kfsr));
  }

  // Use Zll+jet FSR for Zee+jet and Zmm+jet to keep them on same footing
  TH1D *hzm = (TH1D*)d->Get("fsr/hkfsr_mpfchs1_zlljet"); assert(hzm);
  TH1D *hzp = (TH1D*)d->Get("fsr/hkfsr_ptchs_zlljet"); assert(hzp);
  //TH1D *hzm = (TH1D*)d->Get("fsr/hkfsr3_mpfchs1_zlljet"); assert(hzm);
  //TH1D *hzp = (TH1D*)d->Get("fsr/hkfsr3_ptchs_zlljet"); assert(hzp);
  // Correct Zee
  for (int i = 0; i != gem->GetN(); ++i) {
    double kfsr  = hzm->GetBinContent(hzm->FindBin(gem->GetX()[i]));
    gem->SetPoint(i, gem->GetX()[i], gem->GetY()[i]*(1. - 0.3*kfsr));
    //gem->SetPoint(i, gem->GetX()[i], gem->GetY()[i]*(1. - 1.0*kfsr));
  }
  for (int i = 0; i != gep->GetN(); ++i) {
    double kfsr  = hzp->GetBinContent(hzp->FindBin(gep->GetX()[i]));
    gep->SetPoint(i, gep->GetX()[i], gep->GetY()[i]*(1. - 0.3*kfsr));
    //gep->SetPoint(i, gep->GetX()[i], gep->GetY()[i]*(1. - 1.0*kfsr));
  }
  // Correct Zmm
  for (int i = 0; i != gmm->GetN(); ++i) {
    double kfsr  = hzm->GetBinContent(hzm->FindBin(gmm->GetX()[i]));
    gmm->SetPoint(i, gmm->GetX()[i], gmm->GetY()[i]*(1. - 0.3*kfsr));
    //gmm->SetPoint(i, gmm->GetX()[i], gmm->GetY()[i]*(1. - 1.0*kfsr));
  }
  for (int i = 0; i != gmp->GetN(); ++i) {
    double kfsr  = hzp->GetBinContent(hzp->FindBin(gmp->GetX()[i]));
    gmp->SetPoint(i, gmp->GetX()[i], gmp->GetY()[i]*(1. - 0.3*kfsr));
    //gmp->SetPoint(i, gmp->GetX()[i], gmp->GetY()[i]*(1. - 1.0*kfsr));
  }


  curdir->cd();

  double ptmax = 2.*1200.;
  double ptmin = 15.;
  TH1D *hup = new TH1D("hup",";p_{T,Z} (GeV);Data / MC",ptmax-ptmin,ptmin,ptmax);
  /*
  hup->SetMinimum(0.94);//0.96);
  hup->SetMaximum(1.01);
  */
  hup->SetMinimum(0.96);
  hup->SetMaximum(1.03);
  hup->GetXaxis()->SetMoreLogLabels();
  hup->GetXaxis()->SetNoExponent();

  TH1D *hdw = new TH1D("hdw",";p_{T} (GeV);Zee / Z#mu#mu - 1 (%)",ptmax-ptmin,ptmin,ptmax);
  hdw->SetMinimum((0.970+1e-5-1)*100);
  hdw->SetMaximum((1.030+1e-5-1)*100);
  hdw->GetXaxis()->SetMoreLogLabels();
  hdw->GetXaxis()->SetNoExponent();

  map<string, const char*> lumimap;
  //lumimap["A"] = "Run2018A 14.0 fb^{-1}"; //PdmV Analysis TWiki
  //lumimap["B"] = "Run2018B 7.1 fb^{-1}"; //PdmV Analysis TWiki
  //lumimap["C"] = "Run2018C 6.9 fb^{-1}"; //PdmV Analysis TWiki
  //lumimap["D"] = "Run2018D 31.9 fb^{-1}"; //PdmV Analysis TWiki
  //lumimap["ABC"] = "Run2018ABC 28.0 fb^{-1}"; //PdmV Analysis TWiki
  //lumimap["ABCD"] = "Run2018ABCD 59.9 fb^{-1}"; //PdmV Analysis TWiki
  lumimap["BCDEF"] = "2017, 41.5 fb^{-1}"; // for DP note
  lumimap["B"] = "Run2017B, 4.8 fb^{-1}";
  lumimap["C"] = "Run2017C, 9.6 fb^{-1}";
  lumimap["D"] = "Run2017D, 4.2 fb^{-1}";
  lumimap["E"] = "Run2017E, 9.3 fb^{-1}";
  lumimap["F"] = "Run2017F, 13.4 fb^{-1}";
  lumimap["2018ABCD"] = "2018";
  lumimap["2016BCD"] = "Run2016BCD, 12.9 fb^{-1}";
  lumimap["2016EF"] = "Run2016EF, 6.8 fb^{-1}";
  lumimap["2016GH"] = "Run2016GH, 16.8 fb^{-1}";
  lumimap["2016BCDEFGH"] = "Run2016BCDEFGH, 36.5 fb^{-1}";
  lumi_13TeV = lumimap[run];
  TCanvas *c1 = tdrDiCanvas("c1",hdw,hup,4,11);

  c1->cd(2);
  gPad->SetLogx();

  tdrDraw(hmh,"Pz",kFullSquare,kRed+1);
  tdrDraw(heh,"Pz",kFullSquare,kGreen+3);

  //tdrDraw(gmm,"Pz",kFullSquare,kRed);
  //tdrDraw(gmp,"Pz",kOpenSquare,kRed);

  //tdrDraw(gem,"Pz",kFullCircle,kGreen+2);
  //tdrDraw(gep,"Pz",kOpenCircle,kGreen+2);

  tdrDraw(gpm,"Pz",kFullDiamond,kBlue);
  // Don't draw gam+jet pTbal, because different MC => kFSR does not cancel
  tdrDraw(gpp,"Pz",kOpenDiamond,kBlue);

  c1->cd(1);
  gPad->SetLogx();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(ptmin,0,ptmax,0);

  // Ratio new/old

  TGraphErrors *grm = tools::ratioGraphs(gem,gmm);
  TGraphErrors *grp = tools::ratioGraphs(gep,gmp);
  TH1D *hrh = (TH1D*)heh->Clone("hrh"); hrh->Divide(hmh);

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
  TGraphErrors *grh = new TGraphErrors(0);
  for (int i = 1; i != hrh->GetNbinsX()+1; ++i) {
    hrh->SetBinContent(i, (hrh->GetBinContent(i)-1)*100);
    hrh->SetBinError(i, hrh->GetBinError(i)*100);
    if (hrh->GetBinError(i)!=0) {
      int n = grh->GetN();
      grh->SetPoint(n, hrh->GetBinCenter(i), hrh->GetBinContent(i));
      grh->SetPointError(n, hrh->GetBinWidth(i)*0.5, hrh->GetBinError(i));
    }
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



  //TLegend *legdw = tdrLeg(0.5,0.55,0.70,0.90);
  TLegend *legdw = tdrLeg(0.5,0.65,0.70,0.90);
  legdw->AddEntry(hrh,"HDM","PL");
  //legdw->AddEntry(grm,"MPF","PL");
  //legdw->AddEntry(grp,"p_{T} balance","PL");
  legdw->AddEntry(grpm,"#gamma / Z#mu#mu MPF","PL");
  legdw->AddEntry(grpp,"#gamma / Z#mu#mu p_{T}^{bal}","PL");

  tdrDraw(hrh,"Pz",kFullSquare,kGreen+3);
  //tdrDraw(grm,"Pz",kFullCircle,kGreen+2);
  //tdrDraw(grp,"Pz",kOpenSquare,kGreen+2);

  tdrDraw(grpm,"Pz",kFullDiamond,kBlue);
  // Don't draw gam+jet pTbal, because different MC => kFSR does not cancel
  tdrDraw(grpp,"Pz",kOpenDiamond,kBlue);

  TMultiGraph *gm = new TMultiGraph();
  gm->Add(grh);
  //gm->Add(grm);
  //gm->Add(grp);

  TF1 *fs = new TF1("fs","([0]+[1]*log(0.01*x)+[2]*pow(log(0.01*x),2)-1)*100",
		    ptmin,ptmax);
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
		     " - 1)*100",ptmin,ptmax);
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
  //tex->DrawLatex(0.20,0.72,"#alpha < 0.30");
  tex->DrawLatex(0.20,0.72,"#alpha < 1.00");
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
			"-1)*100", ptmin, ptmax);
  // ABCD fit with minitools/drawZmass.C (f1mzee x 3 + f1mzmm x 1)
  //f1mzee->SetParameters(1.00276, 0.00243, 0.00021, 0.99844); // V5M2
  if (run=="2018ABCD")
    // UL18 Run2018ABCD fit with minitools/drawZmass.C
    f1mzee->SetParameters(1.00153, 0.00214, -0.00012, 0.99839); //p3=Zmm
  if (run=="2016GH")
    // UL16 Run2016GH fit with minitools/drawZmass.C                        
    //f1mzee->SetParameters(1.00128, 0.00287, 0.00062, 1.00017); //p3=Zmm,eta25
    f1mzee->SetParameters(1.00175, 0.00349, 0.00159, 1.00017); //p3=Zmm,eta13
  else
    // UL17 RunBCDEF fit with minitools/drawZmass.C
    f1mzee->SetParameters(0.99780, 0.00225, 0.00031, 0.99821); //p3=Zmm
			  //0.99781, 0.00244, 0.00014);  

  // Zee mass applied to gamma+jet at pT,Z=2*pT,gamma
  TF1 *f1mgam = new TF1("f1mgam","([3]/([0]+[1]*log(0.01*x*2)"
                        "+[2]*pow(log(0.01*x*2),2))"
                        "-1)*100", ptmin, ptmax);
  //f1mgam->SetParameters(1.00276, 0.00243, 0.00021, 0.99844); // V5M2
  if (run=="2018ABCD")
    // UL18 Run2018ABCD fit with minitools/drawZmass.C
    f1mgam->SetParameters(1.00153, 0.00214, -0.00012, 0.99839); //p3=Zmm
  if (run=="2016GH")
    // UL16 Run2016GH fit with minitools/drawZmass.C                        
    f1mgam->SetParameters(1.00128, 0.00287, 0.00062, 1.00017); //p3=Zmm
  else
    // UL17 RunBCDEF fit with minitools/drawZmass.C
    //f1mzee->SetParameters(0.99781, 0.00244, 0.00014);
    f1mgam->SetParameters(0.99780, 0.00225, 0.00031, 0.99821); //p3=Zmm

  c1->cd(1);
  f1mzee->SetLineStyle(kSolid);
  f1mzee->DrawClone("SAME");

  f1mgam->SetLineStyle(kDashed);
  f1mgam->DrawClone("SAME");

  legdw->AddEntry(f1mzee,"#mu#mu/ee mass","L");
  legdw->AddEntry(f1mgam,"#mu#mu/#gamma mass","L");

  gPad->Update();

  //c1->SaveAs(Form("pdf/drawZeeVsZmm_Run%s.pdf",run.c_str()));
  //c1->SaveAs(Form("pdf/drawZeeVsZmm_Run%s_nomasscorr.pdf",run.c_str()));
  c1->SaveAs(Form("pdf/drawZeeVsZmm_Run%s_masscorr.pdf",run.c_str()));
} // drawZeeVsZmms

void drawZeeVsZmm() {

  /*
  drawZeeVsZmms("B");
  drawZeeVsZmms("C");
  drawZeeVsZmms("D");
  drawZeeVsZmms("E");
  drawZeeVsZmms("F");
  drawZeeVsZmms("BCDEF");
  */
  //drawZeeVsZmms("2018ABCD");
  drawZeeVsZmms("2016GH");
}
