// Purpose: draw stability of Wb (top candidate) mass
// Uses output from minitools/mk_hadW.C
#include "TFile.h"
#include "TProfile.h"
#include "TLine.h"
#include "TLatex.h"
#include "TGraphErrors.h"

#include "../tdrstyle_mod15.C"

void drawWb(string mode = "18V5") {

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  const char *cm = mode.c_str();

  TFile *fd = new TFile(Form("rootfiles/hadWUL%s_MPDGcorrNoW.root",cm),"READ");
  assert(fd && !fd->IsZombie());

  TFile *fm = new TFile(Form("rootfiles/hadWMC%s_MPDGcorrNoW.root",cm),"READ");
  assert(fm && !fm->IsZombie());

  curdir->cd();

  // bJES+JER, <pTgenb>/<pTrecob>
  TProfile *pgb = (TProfile*)fm->Get("pgb"); assert(pgb);

  TProfile *pmtd = (TProfile*)fd->Get("pmt"); assert(pmtd);
  TProfile *pmtm = (TProfile*)fm->Get("pmt"); assert(pmtm);
  TH1D *hr = pmtd->ProjectionX("pr");
  hr->Add(pmtm,-1);
  hr->Divide(hr,pmtm);
  hr->Scale(100.);

  TProfile *pmtupd = (TProfile*)fd->Get("pmtup"); assert(pmtupd);
  //TProfile *pmtupm = (TProfile*)fm->Get("pmtup"); assert(pmtupm);
  TH1D *hrup = pmtupd->ProjectionX("prup");
  //hrup->Add(pmtupm,-1);
  //hrup->Divide(hrup,pmtupm);
  //hrup->Scale(100.);

  TProfile *pmtdwd = (TProfile*)fd->Get("pmtdw"); assert(pmtdwd);
  //TProfile *pmtdwm = (TProfile*)fm->Get("pmtdw"); assert(pmtdwm);
  TH1D *hrdw = pmtdwd->ProjectionX("prdw");
  //hrdw->Add(pmtdwm,-1);
  //hrdw->Divide(hrdw,pmtdwm);
  //hrdw->Scale(100.);

  TH1D *herr = (TH1D*)hrup->Clone("herr"); herr->Reset();
  TH1D *hrerr = (TH1D*)hrup->Clone("hrerr"); hrerr->Reset();
  for (int i = 1; i != herr->GetNbinsX()+1; ++i) {
    herr->SetBinContent(i, 0.5*(hrup->GetBinContent(i)+hrdw->GetBinContent(i)));
    herr->SetBinError(i, 0.5*(hrup->GetBinContent(i)-hrdw->GetBinContent(i)));
    //hrerr->SetBinContent(i, 100.*(herr->GetBinContent(i)/pmtm->GetBinContent(i)-1));
    hrerr->SetBinContent(i, 0);
    hrerr->SetBinError(i, 100.*herr->GetBinError(i)/pmtm->GetBinContent(i));
  }

  TProfile *pmt1d = (TProfile*)fd->Get("pmt1"); assert(pmt1d);
  TProfile *pmt1m = (TProfile*)fm->Get("pmt1"); assert(pmt1m);
  TH1D *hr1 = pmt1d->ProjectionX("pr1");
  hr1->Add(pmt1m,-1);
  hr1->Divide(hr1,pmt1m);
  hr1->Scale(100.);
  /*
  // Divide manually to check uncertainty => ok, is the same
  for (int i = 1; i != hr1->GetNbinsX()+1; ++i) {
    double x = pmt1d->GetBinContent(i);
    double ex = pmt1d->GetBinError(i);
    double y = pmt1m->GetBinContent(i);
    double ey = pmt1m->GetBinError(i);
    hr1->SetBinContent(i, y!=0 ? 100.*(x/y-1) : 0);
    hr1->SetBinError(i, x*y!=0 ? 100.*x/y*sqrt(pow(ex/x,2)+pow(ey/y,2)) : 0);
  }
  */

  TProfile *pmt2d = (TProfile*)fd->Get("pmt2"); assert(pmt2d);
  TProfile *pmt2m = (TProfile*)fm->Get("pmt2"); assert(pmt2m);
  TH1D *hr2 = pmt2d->ProjectionX("pr2");
  hr2->Add(pmt2m,-1);
  hr2->Divide(hr2,pmt2m);
  hr2->Scale(100.);

  TProfile *pmtmaxd = (TProfile*)fd->Get("pmtmax"); assert(pmtmaxd);
  TProfile *pmtmaxm = (TProfile*)fm->Get("pmtmax"); assert(pmtmaxm);
  TH1D *hrmax = pmtmaxd->ProjectionX("prmax");
  hrmax->Add(pmtmaxm,-1);
  hrmax->Divide(hrmax,pmtmaxm);
  hrmax->Scale(100.);

  TProfile *pmtmind = (TProfile*)fd->Get("pmtmin"); assert(pmtmind);
  TProfile *pmtminm = (TProfile*)fm->Get("pmtmin"); assert(pmtminm);
  TH1D *hrmin = pmtmind->ProjectionX("prmin");
  hrmin->Add(pmtminm,-1);
  hrmin->Divide(hrmin,pmtminm);
  hrmin->Scale(100.);

  TH1D *hup = tdrHist("hup","m_{q#bar{q}'b} (GeV)",145.1,194.9,
		      "p_{T,b} (GeV)",30,240);
  double djes = (mode=="17V5" ? 0.5 : 0.);
  TH1D *hdw = tdrHist("hdw","Data/MC-1 (%)",-1.2+djes,+0.6+djes,
		      "p_{T,b} (GeV)",30,240);

  if (mode=="17V5") lumi_13TeV = "2017, 43.9 fb^{-1}";
  if (mode=="18V5") lumi_13TeV = "2018, 59.9 fb^{-1}";
  TCanvas *c1 = tdrDiCanvas("c1",hup,hdw,4,11);

  TLine *l = new TLine();
  l->SetLineStyle(kDotted);

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);

  c1->cd(1);
  gPad->SetLogx();
  l->DrawLine(30,172.5,240,172.5);

  tdrDraw(herr,"E3",kNone,kYellow+2,kSolid,-1,1001,kYellow+1);

  tdrDraw(pmtmaxm,"Pz",kOpenTriangleUp,kRed);
  tdrDraw(pmtmaxd,"Pz",kFullTriangleUp,kRed);

  tdrDraw(pmtminm,"Pz",kOpenTriangleDown,kBlue);
  tdrDraw(pmtmind,"Pz",kFullTriangleDown,kBlue);

  tdrDraw(pmt2m,"Pz",kOpenDiamond,kOrange-9);
  tdrDraw(pmt2d,"Pz",kFullDiamond,kOrange-9);

  tdrDraw(pmt1m,"Pz",kOpenDiamond,kGreen-9);
  tdrDraw(pmt1d,"Pz",kFullDiamond,kGreen-9);

  //tdrDraw(pmtupm,"PL",kNone,kCyan+1,kSolid,-1);
  //tdrDraw(pmtupd,"PL",kNone,kCyan+1,kSolid,-1);
  //tdrDraw(pmtdwm,"PL",kNone,kCyan+1,kDashed,-1);
  //tdrDraw(pmtdwd,"PL",kNone,kCyan+1,kDashed,-1);

  tdrDraw(pmtm,"Pz",kOpenCircle,kBlack);
  tdrDraw(pmtd,"Pz",kFullCircle,kBlack);

  tex->DrawLatex(0.35,0.06+2*0.06,"|#eta_{jets}|<1.3");
  tex->DrawLatex(0.275,0.06+1*0.06,"60<|m_{q#bar{q}'}|<100 GeV");
  tex->DrawLatex(0.25,0.06+0*0.06,"125<|m_{q#bar{q}'b}|<220 GeV");

  TLegend *legm = tdrLeg(0.61,0.05,0.91,0.10+0.06*6);
  legm->SetHeader("MC");
  legm->AddEntry(pmtmaxm," ","PLE");
  legm->AddEntry(pmt1m," ","PLE");
  legm->AddEntry(pmtm," ","PLE");
  legm->AddEntry(pmt2m," ","PLE");
  legm->AddEntry(pmtminm," ","PLE");

  TLegend *legd = tdrLeg(0.69,0.05,0.99,0.10+0.06*6);
  legd->SetHeader("Data");
  legd->AddEntry(pmtmaxd,"max m_{q#bar{q'}b}","PLE");
  legd->AddEntry(pmt1d,"b-jet 1","PLE");
  legd->AddEntry(pmtd,"b-jet 1+2","PLE");
  legd->AddEntry(pmt2d,"b-jet 2","PLE");
  legd->AddEntry(pmtmind,"min m_{q#bar{q'}b}","PLE");

  c1->cd(2);
  gPad->SetLogx();
  l->DrawLine(30,0,240,0);

  tdrDraw(hrerr,"E3",kNone,kYellow+2,kSolid,-1,1001,kYellow+1);

  tdrDraw(hrmax,"Pz",kFullTriangleUp,kRed);
  tdrDraw(hrmin,"Pz",kFullTriangleDown,kBlue);
  tdrDraw(hr2,"Pz",kFullDiamond,kOrange-9);
  tdrDraw(hr1,"Pz",kFullDiamond,kGreen-9);

  //tdrDraw(hrup,"HL",kNone,kCyan+1,kSolid,-1);
  //tdrDraw(hrdw,"HL",kNone,kCyan+1,kSolid,-1);
  tdrDraw(hr,"Pz",kFullCircle,kBlack);

  c1->SaveAs(Form("pdf/drawWb_%s.pdf",cm));

  // Check stability of Wb quantitatively for all options
  TH1D *h1mtd = (TH1D*)fd->Get("h1mt"); assert(h1mtd);
  TH1D *h1mt1d = (TH1D*)fd->Get("h1mt1"); assert(h1mt1d);
  TH1D *h1mt2d = (TH1D*)fd->Get("h1mt2"); assert(h1mt2d);
  TH1D *h1mtmind = (TH1D*)fd->Get("h1mtmin"); assert(h1mtmind);
  TH1D *h1mtmaxd = (TH1D*)fd->Get("h1mtmax"); assert(h1mtmaxd);
  //
  TH1D *h1mtm = (TH1D*)fm->Get("h1mt"); assert(h1mtm);
  TH1D *h1mt1m = (TH1D*)fm->Get("h1mt1"); assert(h1mt1m);
  TH1D *h1mt2m = (TH1D*)fm->Get("h1mt2"); assert(h1mt2m);
  TH1D *h1mtminm = (TH1D*)fm->Get("h1mtmin"); assert(h1mtminm);
  TH1D *h1mtmaxm = (TH1D*)fm->Get("h1mtmax"); assert(h1mtmaxm);

  TGraphErrors *gad = new TGraphErrors(5);
  gad->SetPoint(0,1,h1mtmind->GetMean());
  gad->SetPointError(0,1,h1mtmind->GetMeanError());

} // drawWb
