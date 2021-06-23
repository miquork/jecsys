// Purpose: draw ATLAS style plots for 7 TeV l+jet 3D method
// (https://arxiv.org/abs/1503.05427)
// Main difference is requiring all jets in barrel at |eta|<1.3
// Kinematic fit is run in full acceptance to select W>qq' and two b's
// Uses output from minitools/mk_hadW.C
#include "TFile.h"
#include "TProfile.h"
#include "TLine.h"
#include "TLatex.h"

#include "../tdrstyle_mod15.C"

//void drawATLASmt(string mode="18V5") {
//void drawATLASmt(string mode="1718V5") {
//void drawATLASmt(string mode="MUF") {
//void drawATLASmt(string mode="EMUF") {
//void drawATLASmt(string mode="Glu_v6") {
//void drawATLASmt(string mode="16V2") {
bool isJEC = true;
//void drawATLASmt(string mode="161718V5") {
void drawATLASmt(string mode="17V5New") {

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  const char *cm = mode.c_str();

  //TFile *fd = new TFile(Form("rootfiles/hadWUL%s_MPDGcorrNoW.root",cm),"READ");
  //TFile *fd = new TFile(Form("rootfiles/hadWUL1718V5_%s.root",cm),"READ");
  TFile *fd(0);
  if (mode=="17V5New") {
    //fd = new TFile(Form("rootfiles/hadWUL%s_SLBRWeight.root",cm),"READ");
    fd = new TFile(Form("rootfiles/hadWUL%s_JECV2.root",cm),"READ");
  }
  else if (mode=="16V2" || mode=="161718V5") {
    fd = new TFile(Form("rootfiles/hadWUL%s_Glu.root",cm),"READ");
    if (isJEC) fd = new TFile(Form("rootfiles/hadWUL%s_JEC.root",cm),"READ");
  }
  else 
    fd = new TFile(Form("rootfiles/hadWUL1718V5_%s.root",cm),"READ");
  assert(fd && !fd->IsZombie());

  //TFile *fm = new TFile(Form("rootfiles/hadWMC%s_MPDGcorrNoW.root",cm),"READ");
  //TFile *fm = new TFile(Form("rootfiles/hadWMC1718V5_%s.root",cm),"READ");
  TFile *fm(0);
  if (mode=="17V5New") {
    //fm = new TFile(Form("rootfiles/hadWMC%s_SLBRWeight.root",cm),"READ");
    fm = new TFile(Form("rootfiles/hadWMC%s_JECV2.root",cm),"READ");
  }
  else if (mode=="16V2" || mode=="161718V5") {
    fm = new TFile(Form("rootfiles/hadWMC%s_Glu.root",cm),"READ");
    if (isJEC) fm = new TFile(Form("rootfiles/hadWMC%s_JEC.root",cm),"READ");
  }
  else fm = new TFile(Form("rootfiles/hadWMC1718V5_%s.root",cm),"READ");
  assert(fm && !fm->IsZombie());

  curdir->cd();

  TH1D *hatmwd = (TH1D*)fd->Get("atlas_mw"); assert(hatmwd);
  TH1D *hatbqd = (TH1D*)fd->Get("atlas_rbq"); assert(hatbqd);
  TH1D *hatmtd = (TH1D*)fd->Get("atlas_mt"); assert(hatmtd);
  TH1D *hatlbd = (TH1D*)fd->Get("atlas_mlb"); assert(hatlbd);
  TH1D *hatmld = (TH1D*)fd->Get("atlas_mlbmet"); assert(hatmld);
  //
  TH1D *hatmwm = (TH1D*)fm->Get("atlas_mw"); assert(hatmwm);
  TH1D *hatbqm = (TH1D*)fm->Get("atlas_rbq"); assert(hatbqm);
  TH1D *hatmtm = (TH1D*)fm->Get("atlas_mt"); assert(hatmtm);
  TH1D *hatlbm = (TH1D*)fm->Get("atlas_mlb"); assert(hatlbm);
  TH1D *hatmlm = (TH1D*)fm->Get("atlas_mlbmet"); assert(hatmlm);
  //
  TH1D *hatmwsw = (TH1D*)fm->Get("atlas_mw_signal"); assert(hatmwsw);
  TH1D *hatmws1 = (TH1D*)fm->Get("atlas_mw_signal_noweight"); assert(hatmws1);
  TH1D *hatbqsw = (TH1D*)fm->Get("atlas_rbq_signal"); assert(hatbqsw);
  TH1D *hatbqs1 = (TH1D*)fm->Get("atlas_rbq_signal_noweight");assert(hatbqs1);
  
  // Summary histograms (mW,Rbq,mt,ml,mlb) x (all, signal, signal+no weights)
  TH1D *hall = new TH1D("hall",";Result;Data/MC-1 (%)",5,0,5);
  //TH1D *hsig = new TH1D("hsig",";Result;Data/MC-1 (%)",5,0,5);
  //TH1D *hsnw = new TH1D("hsnw",";Result;Data/MC-1 (%)",5,0,5);

  //if (mode=="17V5")   lumi_13TeV = "2017, 43.9 fb^{-1}";
  if (mode=="16V2")   lumi_13TeV = "2016GH, 16.8 fb^{-1}";
  if (mode=="17V5")   lumi_13TeV = "2017, 41.5 fb^{-1}";
  if (mode=="18V5")   lumi_13TeV = "2018, 59.9 fb^{-1}";
  //if (mode=="1718V5") lumi_13TeV = "UL17+18, 103.8 fb^{-1}";
  if (mode=="1718V5") lumi_13TeV = "UL17+18, 101.4 fb^{-1}";
  if (mode=="161718V5") lumi_13TeV = "Run 2, 137.9 fb^{-1}";
  if (mode=="MUF") lumi_13TeV = "UL17+18, 101.4 fb^{-1}";
  if (mode=="EMUF") lumi_13TeV = "[TT l+jet] UL17+18, 101.4 fb^{-1}";
  if (mode=="Glu_v6") lumi_13TeV = "[TT l+jet 21.3.1] UL17+18, 101.4 fb^{-1}";
  if (mode=="17V5New")   lumi_13TeV = "2017New, 41.5 fb^{-1}";

  double k(1.7);
  if (mode=="17V5") k = 0.6;
  if (mode=="18V5") k = 0.6*1.355;//*0.333;//0.8*59.9/43.3;
  if (mode=="1718V5") k = 0.6*2.355;
  if (mode=="161718V5") {
    if (isJEC) k = 0.45*2.355;
    else       k = 0.85*2.355;
  }
  if (mode=="MUF") k = 0.6*2.355;
  if (mode=="16V2") k = 0.02;
  if (mode=="17V5New") k = 0.35;
  k *= 2.5;

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);

  TH1D *hmw = tdrHist("hmw","1000 events / GeV",0,0.001*5300*k,
		      "m_{W}^{reco} (GeV)",55,110);
  TCanvas *cmw = tdrCanvas("cmw",hmw,4,11,kSquare);
  //gPad->SetLeftMargin(0.17);//0.15);
  //hmw->SetTitleOffset(1.6,"Y");//1.3,"Y");
  //hmw->GetYaxis()->SetNoExponent(kFALSE);

  hatmwm->Scale(0.001*hatmwd->Integral()/hatmwm->Integral());
  hatmwd->Scale(0.001);
  tdrDraw(hatmwm,"HIST",kNone,kOrange+2,kSolid,-1,1001,kOrange);
  tdrDraw(hatmwd,"PEz",kFullCircle,kBlack);
  hatmwd->SetMarkerSize(0.7);

  if (isJEC) {
    tex->DrawLatex(0.40,0.90-0.045*1,"|#eta_{jets}| < 1.3");
    tex->DrawLatex(0.40,0.90-0.055*3,"fitProb>0.01");
  }
  else {
    tex->DrawLatex(0.40,0.90-0.045*1,"|#eta_{jets}| < 2.5");
    tex->DrawLatex(0.40,0.90-0.055*3,"fitProb>0.2");
  }
  tex->DrawLatex(0.40,0.90-0.045*2,"p_{T,jets}| > 30 GeV");
  //tex->DrawLatex(0.40,0.90-0.045*3,"Anti-k_{T} R=0.4");

  double mw = 100.*(hatmwd->GetMean()/hatmwm->GetMean()-1);
  double mwe = 100.*sqrt(pow(hatmwd->GetMeanError()/hatmwd->GetMean(),2) +
			 pow(hatmwm->GetMeanError()/hatmwm->GetMean(),2));
  tex->DrawLatex(0.70,0.55,"Data/MC-1 =");
  tex->DrawLatex(0.70,0.50,Form("%1.2f#pm%1.2f%%",mw,mwe));
  hall->SetBinContent(1, mw);
  hall->SetBinError(1, mwe);

  if (hatmwsw) {
    tex->DrawLatex(0.69,0.45,Form("(%1.2f#pm%1.2f%%)",
				  100.*(hatmwd->GetMean()/
					hatmwsw->GetMean()-1),
				  100.*sqrt(pow(hatmwd->GetMeanError()/
						hatmwd->GetMean(),2) +
					    pow(hatmwsw->GetMeanError()/
						hatmwsw->GetMean(),2))));
    tex->DrawLatex(0.69,0.40,Form("(%1.2f#pm%1.2f%%)",
				  100.*(hatmwd->GetMean()/
					hatmws1->GetMean()-1),
				  100.*sqrt(pow(hatmwd->GetMeanError()/
						hatmwd->GetMean(),2) +
					    pow(hatmws1->GetMeanError()/
						hatmws1->GetMean(),2))));
  }

  TLegend *legmw = tdrLeg(0.70,0.90-0.05*2,1.00,0.90);
  legmw->AddEntry(hatmwd,"Data","PE");
  legmw->AddEntry(hatmwm,"MC","F");  

  gPad->RedrawAxis();
  gPad->Update();

  TH1D *hbq = tdrHist("hbq","1000 events / 0.03",0,0.001*0.9*3500*k,
		      "R_{bq}^{reco} (GeV)",0.3,3);
  TCanvas *cbq = tdrCanvas("cbq",hbq,4,11,kSquare);
  //gPad->SetLeftMargin(0.15);
  //hbq->SetTitleOffset(1.3,"Y");

  hatbqm->Scale(0.001*hatbqd->Integral()/hatbqm->Integral());
  hatbqd->Scale(0.001);
  tdrDraw(hatbqm,"HISTE",kNone,kOrange+2,kSolid,-1,1001,kOrange);
  tdrDraw(hatbqd,"PEz",kFullCircle,kBlack);
  hatbqd->SetMarkerSize(0.7);

  if (isJEC) {
    tex->DrawLatex(0.40,0.90-0.045*1,"|#eta_{jets}| < 1.3");
    tex->DrawLatex(0.40,0.90-0.055*3,"fitProb>0.01");
  }
  else {
    tex->DrawLatex(0.40,0.90-0.045*1,"|#eta_{jets}| < 2.5");
    tex->DrawLatex(0.40,0.90-0.055*3,"fitProb>0.2");
  }
  tex->DrawLatex(0.40,0.90-0.045*2,"p_{T,jets}| > 30 GeV");
  //tex->DrawLatex(0.40,0.90-0.045*3,"Anti-k_{T} R=0.4");

  double bqd = 100.*(hatbqd->GetMean()/hatbqm->GetMean()-1);
  double bqde = 100.*sqrt(pow(hatbqd->GetMeanError()/hatbqd->GetMean(),2) +
			  pow(hatbqm->GetMeanError()/hatbqm->GetMean(),2));
  tex->DrawLatex(0.70,0.55,"Data/MC-1 = ");
  tex->DrawLatex(0.70,0.50,Form("%1.2f#pm%1.2f%%",bqd,bqde));
  hall->SetBinContent(2,bqd);
  hall->SetBinError(2,bqde);

  if (hatbqsw) {
    tex->DrawLatex(0.69,0.45,Form("(%1.2f#pm%1.2f%%)",
				  100.*(hatbqd->GetMean()/
					hatbqsw->GetMean()-1),
				  100.*sqrt(pow(hatbqd->GetMeanError()/
						hatbqd->GetMean(),2) +
					    pow(hatbqsw->GetMeanError()/
						hatbqsw->GetMean(),2))));
    tex->DrawLatex(0.69,0.40,Form("(%1.2f#pm%1.2f%%)",
				  100.*(hatbqd->GetMean()/
					hatbqs1->GetMean()-1),
				  100.*sqrt(pow(hatbqd->GetMeanError()/
						hatbqd->GetMean(),2) +
					    pow(hatbqs1->GetMeanError()/
						hatbqs1->GetMean(),2))));
  }
  // double chi2 = 0;
  //tex->DrawLatex(0.70,0.45,"#chi^{2}/NDF = ");
  //tex->DrawLatex(0.70,0.40,Form("%1.1f / %d",
				


  TLegend *legbq = tdrLeg(0.70,0.90-0.05*2,1.00,0.90);
  legbq->AddEntry(hatbqd,"Data","PE");
  legbq->AddEntry(hatbqm,"MC","F");  

  gPad->RedrawAxis();
  gPad->Update();

  TH1D *hmt = tdrHist("hmt","1000 events / GeV",0,0.001*3100*k,
		      "m_{top}^{reco} (GeV)",130,220);
  TCanvas *cmt = tdrCanvas("cmt",hmt,4,11,kSquare);
  gPad->SetLeftMargin(0.15);
  hmt->SetTitleOffset(1.3,"Y");

  if (isJEC) {
    tex->DrawLatex(0.40,0.90-0.045*1,"|#eta_{jets}| < 1.3");
    tex->DrawLatex(0.40,0.90-0.055*3,"fitProb>0.01");
  }
  else {
    tex->DrawLatex(0.40,0.90-0.045*1,"|#eta_{jets}| < 2.5");
    tex->DrawLatex(0.40,0.90-0.055*3,"fitProb>0.2");
  }
  tex->DrawLatex(0.40,0.90-0.045*2,"p_{T,jets}| > 30 GeV");
  //tex->DrawLatex(0.40,0.90-0.045*3,"Anti-k_{T} R=0.4");

  hatmtm->Scale(0.001*hatmtd->Integral()/hatmtm->Integral());
  hatmtd->Scale(0.001);
  tdrDraw(hatmtm,"HISTE",kNone,kOrange+2,kSolid,-1,1001,kOrange);
  tdrDraw(hatmtd,"PEz",kFullCircle,kBlack);
  hatmtd->SetMarkerSize(0.7);

  double mt = 100.*(hatmtd->GetMean()/hatmtm->GetMean()-1);
  double mte = 100.*sqrt(pow(hatmtd->GetMeanError()/hatmtd->GetMean(),2) +
			 pow(hatmtm->GetMeanError()/hatmtm->GetMean(),2));
  tex->DrawLatex(0.70,0.55,"Data/MC-1 = ");
  tex->DrawLatex(0.70,0.50,Form("%1.2f#pm%1.2f%%",mt,mte));
  hall->SetBinContent(3, mt);
  hall->SetBinError(3, mte);
				
  TLegend *legmt = tdrLeg(0.70,0.90-0.05*2,1.00,0.90);
  legmt->AddEntry(hatmtd,"Data","PE");
  legmt->AddEntry(hatmtm,"MC","F");  

  gPad->RedrawAxis();
  gPad->Update();


  TH1D *hml = tdrHist("hml","1000 events / GeV",0,0.001*1400*k,
		      "m_{lbmet}^{reco} (GeV)",35,320);
  TCanvas *cml = tdrCanvas("cml",hml,4,11,kSquare);
  gPad->SetLeftMargin(0.15);
  hml->SetTitleOffset(1.3,"Y");

  if (isJEC) {
    tex->DrawLatex(0.40,0.90-0.045*1,"|#eta_{jets}| < 1.3");
    tex->DrawLatex(0.40,0.90-0.055*3,"fitProb>0.01");
  }
  else {
    tex->DrawLatex(0.40,0.90-0.045*1,"|#eta_{jets}| < 2.5");
    tex->DrawLatex(0.40,0.90-0.055*3,"fitProb>0.2");
  }
  tex->DrawLatex(0.40,0.90-0.045*2,"p_{T,jets}| > 30 GeV");
  //tex->DrawLatex(0.40,0.90-0.045*3,"Anti-k_{T} R=0.4");

  hatmld->Rebin(5);  hatmld->Scale(1./5.);
  hatmlm->Rebin(5);  hatmlm->Scale(1./5.);

  hatmlm->Scale(0.001*hatmld->Integral()/hatmlm->Integral());
  hatmld->Scale(0.001);
  tdrDraw(hatmlm,"HISTE",kNone,kOrange+2,kSolid,-1,1001,kOrange);
  tdrDraw(hatmld,"PEz",kFullCircle,kBlack);
  hatmld->SetMarkerSize(0.7);

  double ml = 100.*(hatmld->GetMean()/hatmlm->GetMean()-1);
  double mle = 100.*sqrt(pow(hatmld->GetMeanError()/hatmld->GetMean(),2) +
			 pow(hatmlm->GetMeanError()/hatmlm->GetMean(),2));
  tex->DrawLatex(0.70,0.55,"Data/MC-1 = ");
  tex->DrawLatex(0.70,0.50,Form("%1.2f#pm%1.2f%%",ml,mle));
  hall->SetBinContent(4, ml);
  hall->SetBinError(4, mle);
				
  TLegend *legml = tdrLeg(0.70,0.90-0.05*2,1.00,0.90);
  legml->AddEntry(hatmld,"Data","PE");
  legml->AddEntry(hatmlm,"MC","F");  

  gPad->RedrawAxis();
  gPad->Update();



  TH1D *hlb = tdrHist("hlb","Events / 3 GeV",0,2*2500*k/0.8,
		      "m_{lb}^{reco} (GeV)",35,170);
  TCanvas *clb = tdrCanvas("clb",hlb,4,11,kSquare);
  gPad->SetLeftMargin(0.15);
  hlb->SetTitleOffset(1.3,"Y");

  if (isJEC) {
    tex->DrawLatex(0.40,0.90-0.045*1,"|#eta_{jets}| < 1.3");
    tex->DrawLatex(0.40,0.90-0.055*3,"fitProb>0.01");
  }
  else {
    tex->DrawLatex(0.40,0.90-0.045*1,"|#eta_{jets}| < 2.5");
    tex->DrawLatex(0.40,0.90-0.055*3,"fitProb>0.2");
  }
  tex->DrawLatex(0.40,0.90-0.045*2,"p_{T,jets}| > 30 GeV");
  //tex->DrawLatex(0.40,0.90-0.045*3,"Anti-k_{T} R=0.4");

  double mlb = 100.*(hatlbd->GetMean()/hatlbm->GetMean()-1);
  double mlbe = 100.*sqrt(pow(hatlbd->GetMeanError()/hatlbd->GetMean(),2) +
			  pow(hatlbm->GetMeanError()/hatlbm->GetMean(),2));
  tex->DrawLatex(0.70,0.65,"Data/MC-1 =");
  tex->DrawLatex(0.70,0.60,Form("%1.2f#pm%1.2f%%",mlb,mlbe));
  hall->SetBinContent(5, mlb);
  hall->SetBinError(5, mlbe);
				
  hatlbm->Scale(hatlbd->Integral()/hatlbm->Integral());
  tdrDraw(hatlbm,"HISTE",kNone,kOrange+2,kSolid,-1,1001,kOrange);
  tdrDraw(hatlbd,"PEz",kFullCircle,kBlack);
  hatlbd->SetMarkerSize(0.7);

  gPad->RedrawAxis();
  gPad->Update();



  TH1D *h = tdrHist("h","Data/MC-1 (%)",-0.5,1,"Result",0,5);
  h->GetXaxis()->SetBinLabel(1,"m_{W}");
  h->GetXaxis()->SetBinLabel(2,"R_{bq}");
  h->GetXaxis()->SetBinLabel(3,"m_{t}");
  h->GetXaxis()->SetBinLabel(4,"m_{lbmet}");
  h->GetXaxis()->SetBinLabel(5,"m_{lb}");

  TCanvas *call = tdrCanvas("call",h,4,11,kSquare);
  gPad->SetLeftMargin(0.15);
  h->SetTitleOffset(1.3,"Y");

  tdrDraw(hall,"Pz",kFullCircle);

  if (isJEC) {
    tex->DrawLatex(0.40,0.90-0.045*1,"|#eta_{jets}| < 1.3");
    tex->DrawLatex(0.40,0.90-0.055*3,"fitProb>0.01");
  }
  else {
    tex->DrawLatex(0.40,0.90-0.045*1,"|#eta_{jets}| < 2.5");
    tex->DrawLatex(0.40,0.90-0.055*3,"fitProb>0.2");
  }
  tex->DrawLatex(0.40,0.90-0.045*2,"p_{T,jets}| > 30 GeV");
  //tex->DrawLatex(0.40,0.90-0.045*3,"Anti-k_{T} R=0.4");

  if (isJEC) {
    cmw->SaveAs(Form("pdf/drawATLASmt_mw_%s_eta13_fp001.pdf",cm));
    cbq->SaveAs(Form("pdf/drawATLASmt_rbq_%s_eta13_fp001.pdf",cm));
    cmt->SaveAs(Form("pdf/drawATLASmt_mt_%s_eta13_fp001.pdf",cm));
    cml->SaveAs(Form("pdf/drawATLASmt_ml_%s_eta13_fp001.pdf",cm));
    clb->SaveAs(Form("pdf/drawATLASmt_mlb_%s_eta13_fp001.pdf",cm));

    call->SaveAs(Form("pdf/drawATLASmt_all_%s_eta13_fp001.pdf",cm));
  }
  else {
    cmw->SaveAs(Form("pdf/drawATLASmt_mw_%s.pdf",cm));
    cbq->SaveAs(Form("pdf/drawATLASmt_rbq_%s.pdf",cm));
    cmt->SaveAs(Form("pdf/drawATLASmt_mt_%s.pdf",cm));
    cml->SaveAs(Form("pdf/drawATLASmt_ml_%s.pdf",cm));
    clb->SaveAs(Form("pdf/drawATLASmt_mlb_%s.pdf",cm));

    call->SaveAs(Form("pdf/drawATLASmt_all_%s.pdf",cm));
  }
} // drawATLASmt
