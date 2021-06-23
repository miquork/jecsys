// Purpose: draw stability of mW, Rbq, mt, mlb etc. vs IOV
// Can also draw rho, NPV, mu with or without normalizing by mu
// (add later ttbar xsec?)
#include "TFile.h"
#include "TProfile.h"
#include "TLine.h"

#include "../tdrstyle_mod15.C"

// Type of files:
// TOP = (|eta|<2.5 + fitProb>0.2);
// JEC =  (|eta|<1.3 + fitProb>0.01)
string mode = "TOP";
//string mode = "JEC";


void hadW_IOV() {
  
  setTDRStyle();
  TDirectory *curdir = gDirectory;

  //TFile *f = new TFile("rootfiles/hadWUL161718V5b_JEC.root","READ");
  TFile *f = new TFile(Form("rootfiles/hadWUL161718V7_%s.root",mode.c_str()),
		       "READ");
  assert(f && !f->IsZombie());
  
  //TFile *fm = new TFile("rootfiles/hadWMC161718V5b_JEC.root","READ");
  TFile *fm = new TFile(Form("rootfiles/hadWMC161718V7_%s.root",mode.c_str()),
			"READ");
  assert(fm && !fm->IsZombie());

  curdir->cd();

  TProfile *pw = (TProfile*)f->Get("pmwiov"); assert(pw);
  TProfile *pt = (TProfile*)f->Get("pmtiov"); assert(pt);
  TProfile *pl = (TProfile*)f->Get("pmlbiov"); assert(pl);
  TProfile *p2 = (TProfile*)f->Get("pmlbmetiov"); assert(p2);
  TProfile *pb = (TProfile*)f->Get("prbqiov"); assert(pb);
  TProfile *pbb = (TProfile*)f->Get("ppbbiov"); assert(pbb);
  TProfile *pqq = (TProfile*)f->Get("ppqqiov"); assert(pqq);
  TH1D *hw = pw->ProjectionX("hw");
  TH1D *ht = pt->ProjectionX("ht");
  TH1D *hl = pl->ProjectionX("hl");
  TH1D *h2 = p2->ProjectionX("h2");
  TH1D *hb = pb->ProjectionX("hb");
  TH1D *hbb = pbb->ProjectionX("hbb");
  //TH1D *hb = pbb->ProjectionX("hb");
  TH1D *hqq = pqq->ProjectionX("hqq");
  //TH1D *hb = (TH1D*)hbb->Clone("hb");
  hb->Divide(hbb,hqq);

  TProfile *pwm = (TProfile*)fm->Get("pmwiov"); assert(pwm);
  TProfile *ptm = (TProfile*)fm->Get("pmtiov"); assert(ptm);
  TProfile *plm = (TProfile*)fm->Get("pmlbiov"); assert(plm);
  TProfile *p2m = (TProfile*)fm->Get("pmlbmetiov"); assert(p2m);
  TProfile *pbm = (TProfile*)fm->Get("prbqiov"); assert(pbm);
  TProfile *pbbm = (TProfile*)fm->Get("ppbbiov"); assert(pbbm);
  TProfile *pqqm = (TProfile*)fm->Get("ppqqiov"); assert(pqqm);
  TH1D *hwm = pwm->ProjectionX("hwm");
  TH1D *htm = ptm->ProjectionX("htm");
  TH1D *hlm = plm->ProjectionX("hlm");
  TH1D *h2m = p2m->ProjectionX("h2m");
  TH1D *hbm = pbm->ProjectionX("hbm");
  TH1D *hbbm = pbbm->ProjectionX("hbbm");
  //TH1D *hbm = pbbm->ProjectionX("hbm");
  TH1D *hqqm = pqqm->ProjectionX("hqqm");
  //TH1D *hbm = (TH1D*)hbbm->Clone("hbm");
  hbm->Divide(hbbm,hqqm);
  
  hwm->Scale(1/pw->GetMean(2));//80.4 /0.985);
  htm->Scale(1/pt->GetMean(2));//1/172.5 /0.980);
  hlm->Scale(1/pl->GetMean(2));
  h2m->Scale(1/p2->GetMean(2));
  //hbm->Scale(1/pb->GetMean(2));//1/1.35);
  //hbm->Scale(1/pb->GetMean(2));//1/1.35);
  hbm->Scale(pqq->GetMean(2)/pbb->GetMean(2));//1/1.35);

  hw->Scale(1/pw->GetMean(2));//80.4 /0.985);
  ht->Scale(1/pt->GetMean(2));//1/172.5 /0.980);
  hl->Scale(1/pl->GetMean(2));
  h2->Scale(1/p2->GetMean(2));
  //hb->Scale(1/pb->GetMean(2));//1/1.35);
  //double meanhb = hb->GetMean(2);
  //hb->Scale(1/meanhb);//1/1.35);
  //hb->Scale(1/pb->GetMean(2));//1/1.35);
  hb->Scale(pqq->GetMean(2)/pbb->GetMean(2));//1/1.35);

  // Copy MC into the visible range
  hwm->SetBinContent(pwm->FindBin(290000), hwm->GetBinContent(1));
  hwm->SetBinError(pwm->FindBin(290000), //pwm->GetBinError(1));
		   sqrt(pow(hwm->GetBinError(1),2)+
			pow(pw->GetMeanError(2)/pw->GetMean(2),2)));
  htm->SetBinContent(ptm->FindBin(290000), htm->GetBinContent(1));
  htm->SetBinError(ptm->FindBin(290000), //ptm->GetBinError(1));
		   sqrt(pow(htm->GetBinError(1),2)+
			pow(pt->GetMeanError(2)/pt->GetMean(2),2)));
  hlm->SetBinContent(plm->FindBin(290000), hlm->GetBinContent(1));
  hlm->SetBinError(plm->FindBin(290000), //plm->GetBinError(1));
		   sqrt(pow(hlm->GetBinError(1),2)+
			pow(pl->GetMeanError(2)/pl->GetMean(2),2)));
  h2m->SetBinContent(p2m->FindBin(290000), h2m->GetBinContent(1));
  h2m->SetBinError(p2m->FindBin(290000), //p2m->GetBinError(1));
		   sqrt(pow(h2m->GetBinError(1),2)+
			pow(p2->GetMeanError(2)/p2->GetMean(2),2)));
  hbm->SetBinContent(pbm->FindBin(290000), hbm->GetBinContent(1));
  hbm->SetBinError(pbm->FindBin(290000), //pbm->GetBinError(1));
		   sqrt(pow(hbm->GetBinError(1),2)+
			pow(pb->GetMeanError(2)/pb->GetMean(2),2)));

  //TH1D *h = (TH1D*)hw->ProjectionX("h");
  //h->Reset();
  //h->SetMinimum(0.98);//0.5);
  //h->SetMaximum(1.03);//2.0);
  //h->GetXaxis()->SetRangeUser(273000, 325000);
  //TH1D *h = tdrHist("h","Ratio to mean",0.98,1.03,"Run",272760,326001);
  //TH1D *h = tdrHist("h","Ratio to mean",0.985,1.015,"Run",272760,326001);
  TH1D *h = tdrHist("h","Ratio to mean",0.980,1.020,"Run",272760,326001);
  //h->GetXaxis()->SetNoExponent();

  lumi_13TeV = "UL16V7 + UL17V5 + UL18V5";
  TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);
  gPad->SetLeftMargin(0.17);
  h->SetTitleOffset(1.5,"Y");


  tdrDraw(ht,"Pz",kFullSquare,kBlue);
  tdrDraw(hl,"Pz",kFullSquare,kCyan+2); hl->SetMarkerSize(0.8);
  tdrDraw(h2,"Pz",kFullSquare,kCyan+4); h2->SetMarkerSize(0.6);
  tdrDraw(hw,"Pz",kFullCircle,kGreen+2);
  tdrDraw(hb,"Pz",kFullDiamond,kRed);

  /*
  tdrDraw(htm,"Pz",kOpenSquare,kBlue);
  tdrDraw(hlm,"Pz",kOpenSquare,kCyan+2); hl->SetMarkerSize(0.8);
  tdrDraw(h2m,"Pz",kOpenSquare,kCyan+4); h2->SetMarkerSize(0.6);
  tdrDraw(hwm,"Pz",kOpenCircle,kGreen+2);
  tdrDraw(hbm,"Pz",kOpenDiamond,kRed);
  */

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(272760,1,326001,1);
  l->SetLineStyle(kDotted);
  //l->DrawLine(272760,1.005,326001,1.005);
  //l->DrawLine(272760,0.995,326001,0.995);
  //l->DrawLine(272760,1.002,326001,1.002);
  //l->DrawLine(272760,0.998,326001,0.998);
  l->DrawLine(272760,1.001,326001,1.001);
  l->DrawLine(272760,0.999,326001,0.999);

  l->DrawLine(272760,0.98,272760,1.010); // 16BCDE
  l->DrawLine(284046,0.98,284046,1.010); // 16GH
  l->DrawLine(297020,0.98,297020,1.005); // 17B
  l->DrawLine(306461,0.98,306461,1.005); // 17F
  l->DrawLine(315000,0.98,315000,1.005); // 18A
  l->DrawLine(326001,0.98,326001,1.005); // 18D

  //h->GetXaxis()->SetBinLabel(h->GetXaxis()->FindBin(278000),"2016");
  //h->GetXaxis()->SetBinLabel(h->GetXaxis()->FindBin(302000),"2017");
  //h->GetXaxis()->SetBinLabel(h->GetXaxis()->FindBin(321500),"2018");
  TLatex *tex = new TLatex();
  //tex->SetNDC(); 
  tex->SetTextSize(0.045);
  tex->DrawLatex(276000,0.982,"2016");
  tex->DrawLatex(285000,0.982,"MC/Data");
  tex->DrawLatex(299000,0.982,"2017");
  tex->DrawLatex(318000,0.982,"2018");
  if (mode=="TOP") {
    tex->SetNDC();
    tex->DrawLatex(0.41,0.86,"|#eta| < 2.5");
    tex->DrawLatex(0.41,0.81,"fitProb > 0.2");
  }

  TLegend *leg = tdrLeg(0.68,0.9-5*0.05,0.88,0.90);

  leg->AddEntry(ht,"m_{bqq}","PLE");
  leg->AddEntry(hl,"m_{lb}","PLE");
  leg->AddEntry(h2,"m_{lbmet}","PLE");
  leg->AddEntry(hw,"m_{qq}","PLE");
  leg->AddEntry(hb,"R_{bq}","PLE");
  /*
  leg->AddEntry(htm,"m_{bqq}","PLE");
  leg->AddEntry(hlm,"m_{lb}","PLE");
  leg->AddEntry(h2m,"m_{lbmet}","PLE");
  leg->AddEntry(hwm,"m_{qq}","PLE");
  leg->AddEntry(hbm,"R_{bq}","PLE");
  */

  gPad->Update();

  c1->SaveAs(Form("pdf/hadW_IOV_%s.pdf",mode.c_str()));
} // hadW_IOV
