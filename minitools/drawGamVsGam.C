// Draw ratio of gamma+jet with two different EGamma corrections
#include "../tdrstyle_mod15.C"
#include "../tools.C"

#include "TFile.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TLine.h"
#include "TF1.h"
#include "TMultiGraph.h"

void drawGamVsGam(string run = "GH") {

  setTDRStyle();

  
  //string sn = Form("V6_%s",run.c_str());
  //string so = Form("V2_%s",run.c_str());
  //string so = Form("V1-Jan15_%s",run.c_str());

  //TFile *fn = new TFile(Form("../rootfiles/Gjet_combinationfile_07Aug17_nores_%s.root", sn.c_str()),"READ");
  //TFile *fn = new TFile(Form("../../jecsys2016/rootfiles/Gjet_combinationfile_07Aug17_L2res_%s_2016.root", sn.c_str()),"READ");
  //
  TFile *fn = new TFile(Form("../../jecsys2016/rootfiles/Gjet_combinationfile_07Aug17_V15_L2res_run%s.root", (run=="GH" ? "FGH" : run.c_str())),"READ");
  assert(fn && !fn->IsZombie());

  //TFile *fo = new TFile(Form("../rootfiles/Gjet_combinationfile_07Aug17_nores_%s.root", so.c_str()),"READ");
  //TFile *fo = new TFile(Form("../../jecsys2016/rootfiles/Gjet_combinationfile_07Aug17_nores_%s.root", sn.c_str()),"READ");
  //TFile *fo = new TFile(Form("../../jecsys/rootfiles/Gjet_combinationfile_07Aug17_nores_%s.root", so.c_str()),"READ");
  //
  //TFile *fo = new TFile(Form("../../jecsys/rootfiles/Gjet_combinationfile_07Aug17_nores_%s_%s.root", (run=="GH" ? "V1" : "V2"), (run=="GH" ? "FGH" : run.c_str())),"READ");
  TFile *fo = new TFile(Form("../../jecsys2016/rootfiles/Gjet_combinationfile_07Aug17_L2Res_V6_%s_EGM2.root", (run=="GH" ? "FGH" : run.c_str())),"READ");
assert(fo && !fo->IsZombie());  

  TGraphErrors *gn = (TGraphErrors*)fn->Get("resp_MPFchs_DATA_a30_eta00_13");
  assert(gn);
  TGraphErrors *gn2 = (TGraphErrors*)fn->Get("resp_PtBalchs_DATA_a30_eta00_13");
  assert(gn2);
  TGraphErrors *go = (TGraphErrors*)fo->Get("resp_MPFchs_DATA_a30_eta00_13");
  assert(go);
  TGraphErrors *go2 = (TGraphErrors*)fo->Get("resp_PtBalchs_DATA_a30_eta00_13");
  assert(go2);


  // Same for Z
  TFile *fzn = new TFile(Form("../../jecsys2016/rootfiles/zjet_combination_07Aug2017_Summer16_JECV15_Zmm_%s_2018-09-14.root",run.c_str()),"READ");

  TFile *fzo = new TFile(Form("../../jecsys2016/rootfiles/zjet_combination_07Aug2017_Summer16_JECV6_Zmm_%s_2018-03-06.root",run.c_str()),"READ");
  //zjet_pileupSummary_07Aug2017_Summer16_JECV6_Zmm_%s_2018-03-27.root

  TGraphErrors *gzn = (TGraphErrors*)fzn->Get("Data_MPF_CHS_a30_eta_00_13_L1L2Res");
  assert(gzn);
  TGraphErrors *gzn2 = (TGraphErrors*)fzn->Get("Data_PtBal_CHS_a30_eta_00_13_L1L2Res");
  assert(gzn2);
  TGraphErrors *gzo = (TGraphErrors*)fzo->Get("Data_MPF_CHS_a30_eta_00_13_L1L2L3");//Res");
  assert(gzo);
  TGraphErrors *gzo2 = (TGraphErrors*)fzo->Get("Data_PtBal_CHS_a30_eta_00_13_L1L2L3");//Res");
  assert(gzo2);

//gn->SetMarkerStyle(kFullCircle);
//gn->Draw("AP");

//go->SetMarkerStyle(kOpenCircle);
//go->Draw("SAMEP");
									
  double ptmax = 2.*1200.;
  TH1D *hup = new TH1D("hup",";p_{T,#gamma} (GeV);MPF or p_{T}^{bal}",670,30,ptmax);
  hup->SetMinimum(0.83);//0.84);
  hup->SetMaximum(1.03);//1.00);
  hup->GetXaxis()->SetMoreLogLabels();
  hup->GetXaxis()->SetNoExponent();

  TH1D *hdw = new TH1D("hdw",";p_{T} (GeV);New / Old - 1 (%)",670,30,ptmax);
  hdw->SetMinimum(-2.5+1e-4);//(0.985+1e-5-1)*100);
  hdw->SetMaximum(+4.5);//(1.025+1e-5-1)*100);
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

  // Z+jet 
  tdrDraw(gzn2,"Pz",kFullSquare,kRed-9);
  tdrDraw(gzo2,"Pz",kOpenSquare,kRed-9);

  tdrDraw(gzn,"Pz",kFullCircle,kRed);
  tdrDraw(gzo,"Pz",kOpenCircle,kRed);

  // gamma+jet
  tdrDraw(gn2,"Pz",kFullSquare,kGray+1);
  tdrDraw(go2,"Pz",kOpenSquare,kGray+1);

  tdrDraw(gn,"Pz",kFullCircle,kBlue);
  tdrDraw(go,"Pz",kOpenCircle,kBlue);

  c1->cd(1);
  gPad->SetLogx();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(30,0,ptmax,0);

  // Ratio new/old
  TGraphErrors *gr = tools::ratioGraphs(gn,go);
  TGraphErrors *gr2 = tools::ratioGraphs(gn2,go2);
  TGraphErrors *gzr = tools::ratioGraphs(gzn,gzo);
  TGraphErrors *gzr2 = tools::ratioGraphs(gzn2,gzo2);

  // Change result to (new/old-1)*100
  for (int i = 0; i != gr->GetN(); ++i) {
    gr->SetPoint(i, gr->GetX()[i], (gr->GetY()[i]-1)*100);
    gr->SetPointError(i, gr->GetEX()[i], (gr->GetEY()[i])*100);
  }
  for (int i = 0; i != gr2->GetN(); ++i) {
    gr2->SetPoint(i, gr2->GetX()[i], (gr2->GetY()[i]-1)*100);
    gr2->SetPointError(i, gr2->GetEX()[i], (gr2->GetEY()[i])*100);
  }
  //
  for (int i = 0; i != gzr->GetN(); ++i) {
    gzr->SetPoint(i, gzr->GetX()[i], (gzr->GetY()[i]-1)*100);
    gzr->SetPointError(i, gzr->GetEX()[i], (gzn->GetEY()[i])*100);
  }
  for (int i = 0; i != gzr2->GetN(); ++i) {
    gzr2->SetPoint(i, gzr2->GetX()[i], (gzr2->GetY()[i]-1)*100);
    gzr2->SetPointError(i, gzr2->GetEX()[i], (gzn2->GetEY()[i])*100);
  }

  //TLegend *legdw = tdrLeg(0.5,0.70,0.70,0.90);
  TLegend *legdw = tdrLeg(0.5,0.78,0.70,0.90);
  legdw->AddEntry(gr,"MPF","PL");
  legdw->AddEntry(gr2,"p_{T} balance","PL");
  TLegend *legz = tdrLeg(0.45,0.78,0.65,0.90);
  legz->AddEntry(gzr," ","PL");
  legz->AddEntry(gzr2," ","PL");

  tdrDraw(gzr2,"Pz",kFullSquare,kRed-9);
  tdrDraw(gzr,"Pz",kFullCircle,kRed);
  tdrDraw(gr2,"Pz",kFullSquare,kGray+1);
  tdrDraw(gr,"Pz",kFullCircle,kBlue);

  TMultiGraph *gm = new TMultiGraph();
  gm->Add(gr);
  gm->Add(gr2);

  //TF1 *fs = new TF1("fs","[0]*(x<450)+[1]*(x>450)",30,ptmax);
  TF1 *fs = new TF1("fs","[0]*(x<400)+[1]*(x>500)"
		    "+0.5*([0]+[1])*(x>400 && x<500)",30,ptmax);
  fs->SetParameters(0.5,-0.5);
  
  gm->Fit(fs,"QRN");
  fs->SetLineStyle(kSolid);//Dotted);
  //fs->SetLineWidth(2);
  fs->SetLineColor(kBlue);//kBlack);
  fs->Draw("SAME");
  
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);
  tex->DrawLatex(0.20,0.72,"#alpha < 0.30");
  tex->DrawLatex(0.20,0.67,"Data only");


  tex->SetTextColor(kBlue);
  tex->DrawLatex(0.20,0.23,Form("#Delta = %1.2f #pm %1.2f%%",
				fs->GetParameter(1)-fs->GetParameter(0),
				sqrt(pow(fs->GetParError(0),2)
				     +pow(fs->GetParError(1),2))));
  tex->DrawLatex(0.20,0.17,Form("p_{T} < 400 GeV: %1.2f #pm %1.2f%%",
				fs->GetParameter(0), fs->GetParError(0)));
  tex->DrawLatex(0.20,0.11,Form("p_{T} > 500 GeV: %1.2f #pm %1.2f%%",
				fs->GetParameter(1), fs->GetParError(1)));
  tex->DrawLatex(0.20,0.05,Form("#chi^{2} / NDF = %1.1f / %d",
				fs->GetChisquare(), fs->GetNDF()));


  // Compare to Zee mass corrections
  TF1 *f1mzee = new TF1("f1mzee","(1./([0]+[1]*log(2*x)+[2]*log(2*x)*log(2*x))"
			"-1)*100", 30, ptmax);
  TF1 *f1ezee = new TF1("f1ezee","sqrt([0]+pow(log(2*x),2)*[1]"
			"+pow(log(2*x),4)*[2]+2*log(2*x)*[3]"
			"+2*pow(log(2*x),2)*[4]+2*pow(log(2*x),3)*[5])"
			,30,ptmax);
  // BCDEFGH fit with minitools/drawZmass.C
  f1mzee->SetParameters(1.01732, -0.00909, 0.00116);
  f1ezee->SetParameters(7.54e-05, 1.41e-05, 1.63e-07,
			-3.26e-05, 3.47e-06, -1.51e-06);

  //f1mzee->SetParameter(0,f1mzee->GetParameter(0)-fs->GetParameter(0)/100.);
  f1mzee->SetParameter(0,f1mzee->GetParameter(0)-0.005);
  //f1mzee->DrawClone("SAME");
  f1mzee->SetLineStyle(kDotted);
  f1mzee->SetParameter(0,f1mzee->GetParameter(0)-0.005);
  //f1mzee->DrawClone("SAME");
  f1mzee->SetParameter(0,f1mzee->GetParameter(0)+0.010);
  //f1mzee->DrawClone("SAME");

  f1mzee->SetLineStyle(kSolid);
  //legdw->AddEntry(f1mzee,"Z#rightarrow ee mass + 0.5%","L");

  c1->SaveAs(Form("../pdf/drawGamVsGam_Run%s.pdf",run.c_str()));
}
