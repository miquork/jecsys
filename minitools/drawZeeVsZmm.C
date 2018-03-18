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

void drawZeeVsZmm(string run = "BCDEFGH") {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  /*
  map<string,const char*> fz_files;
  fz_files["BCD"] = "BCD_2017-12-08";
  fz_files["EF"] = "EF_2017-12-08";
  fz_files["G"] = "G_2017-12-08";
  fz_files["H"] = "H_2017-12-08";
  fz_files["GH"] = "GH_2017-12-19";
  fz_files["BCDEFGH"] = "BCDEFGH_2017-12-08"; 

  TFile *fm = new TFile(Form("../rootfiles/combination_ZJet_Zmm_%s_wide.root", fz_files[run]),"READ");
  TFile *fe = new TFile(Form("../rootfiles/combination_ZJet_Zee_%s_wide.root",fz_files[run]),"READ");

  assert(fm && !fm->IsZombie());
  assert(fe && !fe->IsZombie());

  // NB: Ratio is TH1D inoriginal files...
  TGraphErrors *gem = (TGraphErrors*)fe->Get("Data_MPF_CHS_a30_eta_00_13_L1L2L3");
  assert(gem);
  TGraphErrors *gep = (TGraphErrors*)fe->Get("Data_PtBal_CHS_a30_eta_00_13_L1L2L3");
  assert(gep);
  TGraphErrors *gmm = (TGraphErrors*)fm->Get("Data_MPF_CHS_a30_eta_00_13_L1L2L3");
  assert(gmm);
  TGraphErrors *gmp = (TGraphErrors*)fm->Get("Data_PtBal_CHS_a30_eta_00_13_L1L2L3");
  assert(gmp);
  */
  
  TFile *f = new TFile(Form("../rootfiles/jecdata%s.root",run.c_str()),"READ");
  assert(f && !f->IsZombie());
  gDirectory->cd("ratio");
  gDirectory->cd("eta00-13");
  TDirectory *d = gDirectory;
  TGraphErrors *gem = (TGraphErrors*)d->Get("mpfchs1_zeejet_a30"); assert(gem);
  TGraphErrors *gep = (TGraphErrors*)d->Get("ptchs_zeejet_a30"); assert(gep);
  TGraphErrors *gmm = (TGraphErrors*)d->Get("mpfchs1_zmmjet_a30"); assert(gmm);
  TGraphErrors *gmp = (TGraphErrors*)d->Get("ptchs_zmmjet_a30"); assert(gmp);
  curdir->cd();

  double ptmax = 2.*1200.;
  TH1D *hup = new TH1D("hup",";p_{T,Z} (GeV);Data / MC",670,30,ptmax);
  hup->SetMinimum(0.97);
  hup->SetMaximum(1.02);
  hup->GetXaxis()->SetMoreLogLabels();
  hup->GetXaxis()->SetNoExponent();

  TH1D *hdw = new TH1D("hdw",";p_{T} (GeV);Zee / Z#mu#mu - 1 (%)",670,30,ptmax);
  hdw->SetMinimum((0.970+1e-5-1)*100);
  hdw->SetMaximum((1.030+1e-5-1)*100);
  hdw->GetXaxis()->SetMoreLogLabels();
  hdw->GetXaxis()->SetNoExponent();

  map<string, const char*> lumimap;
  lumimap["BCD"] = "Run2016BCD Legacy, 12.9 fb^{-1}";
  lumimap["EF"] = "Run2016EF Legacy, 6.8 fb^{-1}";
  lumimap["FG"] = "Run2016fG Legacy, 8.0 fb^{-1}";
  lumimap["H"] = "Run2016H Legacy, 8.8 fb^{-1}";
  lumimap["GH"] = "Run2016fGH Legacy, 16.8 fb^{-1}";
  lumimap["BCDEFGH"] = "Run2016BCDEFGH Legacy, 36.5 fb^{-1}";
  lumi_13TeV = lumimap[run];
  TCanvas *c1 = tdrDiCanvas("c1",hdw,hup,4,11);

  c1->cd(2);
  gPad->SetLogx();

  tdrDraw(gmm,"Pz",kFullSquare,kRed);
  tdrDraw(gmp,"Pz",kOpenSquare,kRed);

  tdrDraw(gem,"Pz",kFullCircle,kGreen+2);
  tdrDraw(gep,"Pz",kOpenCircle,kGreen+2);

  c1->cd(1);
  gPad->SetLogx();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(30,0,ptmax,0);

  // Ratio new/old

  TGraphErrors *grm = tools::ratioGraphs(gem,gmm);
  TGraphErrors *grp = tools::ratioGraphs(gep,gmp);

  // Change result to (Zee/Zmm-1)*100
  for (int i = 0; i != grm->GetN(); ++i) {
    grm->SetPoint(i, grm->GetX()[i], (grm->GetY()[i]-1)*100);
    grm->SetPointError(i, grm->GetEX()[i], (grm->GetEY()[i])*100);
  }
  for (int i = 0; i != grp->GetN(); ++i) {
    grp->SetPoint(i, grp->GetX()[i], (grp->GetY()[i]-1)*100);
    grp->SetPointError(i, grp->GetEX()[i], (grp->GetEY()[i])*100);
  }

  TLegend *legdw = tdrLeg(0.5,0.70,0.70,0.90);
  legdw->AddEntry(grm,"MPF","PL");
  legdw->AddEntry(grp,"p_{T} balance","PL");

  tdrDraw(grm,"Pz",kFullCircle,kGreen+2);
  tdrDraw(grp,"Pz",kOpenSquare,kGreen+2);

  TMultiGraph *gm = new TMultiGraph();
  gm->Add(grm);
  gm->Add(grp);

  TF1 *fs = new TF1("fs","[0]+[1]*log(0.01*x)+[2]*pow(log(0.01*x),2)",30,ptmax);
  fs->SetParameters(0,0,0);//1.7,-1,1);
  fs->FixParameter(2,0);

  gm->Fit(fs,"QRN");
  fs->SetLineStyle(kSolid);
  fs->SetLineColor(kGreen+2);
  fs->Draw("SAME");
  
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);
  tex->DrawLatex(0.20,0.72,"#alpha < 0.30");
  tex->DrawLatex(0.20,0.67,"Data / MC ratio");
  //
  tex->DrawLatex(0.20,0.20,Form("#chi^{2}/NDF = %1.1f/%d",
				fs->GetChisquare(), fs->GetNDF()));
  tex->DrawLatex(0.20,0.15,Form("p_{0} = %1.2f #pm %1.2f",
				fs->GetParameter(0),fs->GetParError(0)));
  tex->DrawLatex(0.20,0.10,Form("p_{1} = %1.2f #pm %1.2f",
				fs->GetParameter(1),fs->GetParError(1)));
  tex->DrawLatex(0.20,0.05,Form("p_{2} = %1.2f #pm %1.2f",
				fs->GetParameter(2),fs->GetParError(2)));


  // Compare to Zee mass corrections
  /*
  TF1 *f1mzee = new TF1("f1mzee","([3]/([0]+[1]*log(2*x)+[2]*log(2*x)*log(2*x))"
			"-1)*100", 30, ptmax);
  TF1 *f1ezee = new TF1("f1ezee","sqrt([0]+pow(log(2*x),2)*[1]"
			"+pow(log(2*x),4)*[2]+2*log(2*x)*[3]"
			"+2*pow(log(2*x),2)*[4]+2*pow(log(2*x),3)*[5])"
			,30,ptmax);
  */
  TF1 *f1mzee = new TF1("f1mzee","([3]/([0]+[1]*log(0.01*x)"
			"+[2]*pow(log(0.01*x),2))"
			"-1)*100", 30, ptmax);
  //TF1 *f1ezee = new TF1("f1ezee","sqrt([0]+pow(log(2*x),2)*[1]"
  //			"+pow(log(2*x),4)*[2]+2*log(2*x)*[3]"
  //			"+2*pow(log(2*x),2)*[4]+2*pow(log(2*x),3)*[5])"
  //			,30,ptmax);
  // BCDEFGH fit with minitools/drawZmass.C
  //f1mzee->SetParameters(1.01732, -0.00909, 0.00116, 0.99855);
  f1mzee->SetParameters(1.00017, 0.00166, 0.00114, 0.99868);
  //f1ezee->SetParameters(7.54e-05, 1.41e-05, 1.63e-07,
  //			-3.26e-05, 3.47e-06, -1.51e-06);

  c1->cd(1);
  f1mzee->DrawClone("SAME");
  f1mzee->SetLineStyle(kSolid);

  legdw->AddEntry(f1mzee,"#mu#mu/ee mass","L");

  gPad->Update();

  c1->SaveAs(Form("../pdf/drawZeeVsZmm_Run%s.pdf",run.c_str()));
}
