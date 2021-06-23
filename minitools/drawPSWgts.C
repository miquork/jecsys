// Purpose: draw PSWeights variations in TT events
#include "TFile.h"
#include "TProfile.h"
#include "TLine.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TF1.h"

#include "../tdrstyle_mod15.C"

void drawPSWgts() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  //TFile *fps = new TFile("rootfiles/HadW/UL18V5/WMass_Muo18_PowHeg.root",
  TFile *fps = new TFile("rootfiles/HadW/UL18V5_WithEMu/Muo18_MC.root",
			 "READ");
  assert(fps && !fps->IsZombie());

  TH1D *hps = (TH1D*)fps->Get("PSWeights"); assert(hps);

  //TFile *f = new TFile("rootfiles/hadWMC18V5_MPDGcorrNoW.root","READ");
  //TFile *f = new TFile("rootfiles/hadWMC1718V5_MPDGcorrNoW.root","READ");
  TFile *f = new TFile("rootfiles/hadWMC161718V5_Glu.root","READ");
  assert(f && !f->IsZombie());

  //string obs[] = {"mt","mw","mlb","rbq","rbb","rqq","mwb","ptqq","ptw"};
  //string obs[] = {"mt","mw","mlb","mwb","mwbk","rq","rb"};//
		  //"ptqq","ptw","ptb","ptbk"};
  string obs[] = {"gmt","gmw","gmlb","gmwb","gmwbk"};
  //string obs[] = {"rbq","rbb","rqq","rb","rq"};
  const int nobs = sizeof(obs)/sizeof(obs[0]);

  map<string, string> label;
  label["mt"] = "m_{qqb} _{(0D)}";
  label["mw"] = "m_{qq}";
  label["mwb"] = "m_{Wb} _{(1D)}";
  label["mwbk"] = "m_{WB} _{(2D)}";
  label["mwB"] = "m_{WB} _{(2D)}";
  label["mlb"] = "m_{lb}";
  label["rbq"] = "R_{bq}";
  label["rbb"] = "R_{bb}";
  label["rqq"] = "R_{qq}";
  label["ptb"] = "p_{T,b}";
  label["ptbk"] = "p_{T,B} _{(2D)}";
  label["ptqq"] = "p_{T,qq}";
  label["ptw"] = "p_{T,W} _{(1D)}";
  label["rb"] = "R_{b}";
  label["rq"] = "R_{q}";

  label["gmt"] = "M_{qqb} _{(0D)}";
  label["gmw"] = "M_{qq}";
  label["gmlb"] = "M_{lb}";
  label["gmwb"] = "M_{Wb} _{(1D)}";
  label["gmwbk"] = "M_{WB} _{(2D)}";

  map<string,int> color;
  color["mt"] = kBlack;
  color["mwb"] = kGray+2;
  color["mwbk"] = kGray+1;
  color["mwB"] = kGray+1;
  color["mlb"] = kBlue;
  color["mw"] = kRed;
  color["rbq"] = kGreen+2;
  color["rbb"] = kCyan+2;
  color["rqq"] = kMagenta+2;
  color["ptb"] = kCyan+2;
  color["ptbk"] = kCyan+2;
  color["ptqq"] = kGreen+2;//kRed;//kCyan+2;
  color["ptw"] = kGreen+2;//kMagenta+2;
  color["rb"] = kCyan+2;
  color["rq"] = kGreen+2;

  color["gmt"] = kBlack;
  color["gmw"] = kRed;
  color["gmlb"] = kBlue;
  color["gmwb"] = kGray+2;
  color["gmwbk"] = kGray+1;

  map<string,int> marker;
  marker["mt"] = kFullCircle;
  marker["mwb"] = kFullCircle;
  marker["mwbk"] = kFullCircle;
  marker["mwB"] = kFullCircle;
  marker["mlb"] = kOpenCircle;
  marker["mw"] = kOpenSquare;
  marker["rbq"] = kOpenSquare;
  marker["rbb"] = kOpenDiamond;
  marker["rqq"] = kOpenDiamond;
  marker["ptb"] = kOpenDiamond;
  marker["ptbk"] = kFullDiamond;
  marker["ptqq"] = kOpenDiamond;//kFullDiamond;
  marker["ptw"] = kFullDiamond;
  marker["rb"] = kFullDiamond;//kOpenDiamond;
  marker["rq"] = kFullDiamond;

  marker["gmt"] = kFullCircle;
  marker["gmw"] = kOpenSquare;
  marker["gmlb"] = kOpenCircle;
  marker["gmwb"] = kFullCircle;
  marker["gmwbk"] = kFullCircle;

  TF1 *f1 = new TF1("f1","-100",0.5,45.5);
  vector<TProfile*> ps(nobs);
  vector<TH1D*> hs(nobs);
  for (int i = 0; i != nobs; ++i) {
    const char *c = obs[i].c_str();
    TProfile *p = (TProfile*)f->Get(Form("ps%s",c)); assert(p);
    TH1D *h = p->ProjectionX(Form("hps%s",c));
    h->Scale(100./h->GetBinContent(1));
    h->Add(f1);

    ps[i] = p;
    hs[i] = h;
  }

  /*
  TProfile *psmw = (TProfile*)f->Get("psmw"); assert(psmw);
  TProfile *psmt = (TProfile*)f->Get("psmt"); assert(psmt);
  TProfile *psmlb = (TProfile*)f->Get("psmlb"); assert(psmlb);
  TProfile *psrbq = (TProfile*)f->Get("psrbq"); assert(psrbq);
  TProfile *psrbb = (TProfile*)f->Get("psrbb"); assert(psrbb);
  TProfile *psrqq = (TProfile*)f->Get("psrqq"); assert(psrqq);

  TH1D *hpsmw = psmw->ProjectionX("hpsmw");
  TH1D *hpsmt = psmt->ProjectionX("hpsmt");
  TH1D *hpsmlb = psmlb->ProjectionX("hpsmlb");
  TH1D *hpsrbq = psrbq->ProjectionX("hpsrbq");
  TH1D *hpsrbb = psrbb->ProjectionX("hpsrbb");
  TH1D *hpsrqq = psrqq->ProjectionX("hpsrqq");

  hpsmw->Scale(100./hpsmw->GetBinContent(1));
  hpsmt->Scale(100./hpsmt->GetBinContent(1));
  hpsmlb->Scale(100./hpsmlb->GetBinContent(1));
  hpsrbq->Scale(100./hpsrbq->GetBinContent(1));
  hpsrbb->Scale(100./hpsrbb->GetBinContent(1));
  hpsrqq->Scale(100./hpsrqq->GetBinContent(1));

  TF1 *f1 = new TF1("f1","-100",0.5,45.5);
  hpsmw->Add(f1);
  hpsmt->Add(f1);
  hpsmlb->Add(f1);
  hpsrbq->Add(f1);
  hpsrbb->Add(f1);
  hpsrqq->Add(f1);
  */

  //TH1D *h = tdrHist("h","PS variation",0.992,1.008,"iPS",0.5,51.5);//45.5);
  TH1D *h = tdrHist("h","PS variation (%)",-0.8,+1.0,"iPS",0.5,51.5);//45.5);
  
  // Clone to get axis labels upright
  //TH1D *h = (TH1D*)hps->Clone("h"); h->Reset();
  //h->SetMaximum(1.008);
  //h->SetMinimum(0.992);

  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    h->GetXaxis()->SetBinLabel(i, hps->GetXaxis()->GetBinLabel(i+1));
  }
  h->GetXaxis()->SetLabelFont(hps->GetXaxis()->GetLabelFont());

  //lumi_13TeV = "Powheg+Pythia8, UL2018";
  //lumi_13TeV = "Powheg+Pythia8, UL17+18";
  lumi_13TeV = "Powheg+Pythia8, UL16+17+18";
  TCanvas *c1 = tdrCanvas("c1",h,4,11);
  gPad->SetBottomMargin(0.30);

  /*
  tdrDraw(hpsrbb,"HP",kOpenDiamond,kCyan+2,kSolid,-1,kNone);
  tdrDraw(hpsrqq,"HP",kOpenDiamond,kMagenta+2,kSolid,-1,kNone);
  tdrDraw(hpsrbq,"HP",kOpenSquare,kGreen+2,kSolid,-1,kNone);
  tdrDraw(hpsmw,"HP",kOpenSquare,kRed,kSolid,-1,kNone);
  tdrDraw(hpsmlb,"HP",kOpenCircle,kBlue,kSolid,-1,kNone);
  tdrDraw(hpsmt,"HP",kFullCircle,kBlack,kSolid,-1,kNone);
  */
  for (int i = 0; i != nobs; ++i) {
    const char *c = obs[i].c_str();
    tdrDraw(hs[i],"HP",marker[c],color[c],kSolid,-1,kNone);
  } // for i


  TLine *l = new TLine();
  l->SetLineStyle(kDotted);

  l->DrawLine(7.5,0.6,7.5,-0.8);
  l->DrawLine(11.5,0.6,11.5,-0.8);
  l->DrawLine(13.5,0.6,13.5,-0.8);
  l->DrawLine(15.5,0.6,15.5,-0.8);

  gPad->Update();


  TH1D *h2 = tdrHist("h2","PS variation (%)",-0.28,0.42,
		     "iPS",0.5,6.5);
  TCanvas *c2 = tdrCanvas("c2",h2,4,11);
  //gPad->SetBottomMargin(0.30);

  l->DrawLine(0.5,0,6.5,0);

  // Outline the rough effect from FSR: -0.2% on mW and mt, -0.1% from parts
  l->DrawLine(0.5,-0.2,4.5,-0.2); // -0.2%
  l->DrawLine(0.5,-0.1,5.5,-0.1); // -0.1%
  l->DrawLine(3.5,+0.1,6.5,+0.1); // +0.1%

  l->SetLineColor(kGray);
  l->DrawLine(0.5,0.05,5.5,0.05); // +0.05%
  l->DrawLine(0.5,-.15,2.5,-.15); // -0.15%

  
  /*
  TGraphErrors *gmw = new TGraphErrors(6);
  TGraphErrors *gmt = new TGraphErrors(6);
  TGraphErrors *gmlb = new TGraphErrors(6);
  TGraphErrors *grbq = new TGraphErrors(6);
  TGraphErrors *grbb = new TGraphErrors(6);
  TGraphErrors *grqq = new TGraphErrors(6);
  TGraphErrors *gs[] = {gmw, gmt, gmlb, grbq, grbb, grqq};
  TH1D *hs[] = {hpsmw, hpsmt, hpsmlb, hpsrbq, hpsrbb, hpsrqq};
  const int ng = sizeof(gs)/sizeof(gs[0]);
  */

  //TLegend *leg2 = tdrLeg(0.79,0.89-nobs*0.05,0.99,0.89);
TLegend *leg2 = tdrLeg(0.79,0.89-nobs*0.040,0.99,0.89);

  //for (int i = 0; i != ng; ++i) {
    //TGraphErrors *g = gs[i];
  for (int i = 0; i != nobs; ++i) {
    TGraphErrors *g = new TGraphErrors(6);
    TH1D *h = hs[i];
    double ke = (i<2 ? 1 : 0);
    double d(0);
    //if (nobs==9) d = (i<6 ? 0 : 0.25+0.125);
    if (nobs==9) d = (i<5 ? 0 : 0.25);
    // Baseline
    h2->GetXaxis()->SetBinLabel(1,"Baseline");
    g->SetPoint     (0, 1+d,   h->GetBinContent(1));
    g->SetPointError(0, 0.5-d, h->GetBinError  (1));
    // FSR overall ...
    h2->GetXaxis()->SetBinLabel(2,"FSR");
    g->SetPoint     (1, 2+d,   h->GetBinContent(4));
    g->SetPointError(1, 0.5-d, h->GetBinError  (4)*ke);
    //g->SetPoint     (2, 2,   h->GetBinContent(5));
    //g->SetPointError(2, 0.5, h->GetBinError  (5)*0);
    // ...and its three main components
    h2->GetXaxis()->SetBinLabel(3,"FSR:g2gg");
    g->SetPoint     (3, 3+d,   h->GetBinContent(8));
    g->SetPointError(3, 0.25,  h->GetBinError  (8)*ke);
    //g->SetPoint     (4, 3,   h->GetBinContent(9));
    //g->SetPointError(4, 0.25,h->GetBinError  (9)*0);
    h2->GetXaxis()->SetBinLabel(4,"FSR:q2qg");
    g->SetPoint     (5, 4+d,   h->GetBinContent(12));
    g->SetPointError(5, 0.25,  h->GetBinError  (12)*ke);
    //g->SetPoint     (6, 4,   h->GetBinContent(13));
    //g->SetPointError(6, 0.25,h->GetBinError  (13)*0);
    h2->GetXaxis()->SetBinLabel(5,"FSR:x2xg");
    g->SetPoint     (7, 5+d,   h->GetBinContent(14));
    g->SetPointError(7, 0.25,  h->GetBinError  (14)*ke);
    //g->SetPoint     (8, 5,   h->GetBinContent(15));
    //g->SetPointError(8, 0.25,h->GetBinError  (15)*0);
    // ISR overall ... +11-12
    h2->GetXaxis()->SetBinLabel(6,"ISR");
    g->SetPoint     (9, 6+d,   h->GetBinContent(26));
    g->SetPointError(9, 0.5,   h->GetBinError  (26)*ke);
    //g->SetPoint     (10,6,   h->GetBinContent(27));
    //g->SetPointError(10,0.5, h->GetBinError  (27)*ke);

    const char *c = obs[i].c_str();
    tdrDraw(g,"Pz",marker[c],color[c],kSolid,-1,kNone);
    leg2->AddEntry(g,label[c].c_str(),"PLE");

    if (marker[c]==kFullDiamond||marker[c]==kOpenDiamond)
      g->SetMarkerSize(1.5);
  } // for i

  /*
  tdrDraw(grbb,"Pz",kOpenDiamond,kCyan+2,kSolid,-1,kNone);
  tdrDraw(grqq,"Pz",kOpenDiamond,kMagenta+2,kSolid,-1,kNone);
  tdrDraw(grbq,"Pz",kOpenSquare,kGreen+2,kSolid,-1,kNone);
  tdrDraw(gmw, "Pz",kOpenSquare,kRed,kSolid,-1,kNone);
  tdrDraw(gmlb,"Pz",kOpenCircle,kBlue,kSolid,-1,kNone);
  tdrDraw(gmt, "Pz",kFullCircle,kBlack,kSolid,-1,kNone);
  */  

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  //tex->DrawLatex(0.34,0.85,"|#eta_{jets}|<1.3, 1l2b2^{+}q, fitProb>0.2");
  //tex->DrawLatex(0.34,0.85,"|#eta_{jets}|<2.5, 1l2b2^{+}q, fitProb>0.2");
  tex->DrawLatex(0.34,0.85,"|#eta_{jets}|<1.3, 1l2b2^{+}q, fitProb>0.01");
  tex->SetTextSize(0.040);
  tex->DrawLatex(0.34,0.79,"  60<m_{qq} <100 GeV");
  tex->DrawLatex(0.34,0.75,"130<m_{qqb}<220 GeV");
  tex->DrawLatex(0.34,0.71,"  30<m_{lb}  <170 GeV");
  tex->DrawLatex(0.34,0.67," 0.3<R_{bq}  <3.0,  #mu_{Rfac} = 0.5");
  // Is b included in x, or only t? See
  // http://home.thep.lu.se/~torbjorn/pythia82html/Variations.html
  //tex->DrawLatex(0.67,0.02,"x=t or b+t?");
  tex->DrawLatex(0.72,0.03,"x=b,t"); // According to slides Hannu linked to

  /*
  TLegend *leg2 = tdrLeg(0.79,0.89-6*0.05,0.99,0.89);
  //leg2->SetHeader("#mu=0.5 (2.0)");
  //leg2->SetHeader("#mu_{Rfac} = 0.5");
  leg2->AddEntry(gmt,"m_{qqb}","PLE");
  leg2->AddEntry(gmlb,"m_{lb}","PLE");
  leg2->AddEntry(gmw,"m_{W}","PLE");
  leg2->AddEntry(grbq,"R_{bq}","PLE");
  leg2->AddEntry(grbb,"p_{T,bb}","PLE");
  leg2->AddEntry(grqq,"p_{T,qq}","PLE");
  */

  c1->SaveAs("pdf/drawPSWgts_Run2_all.pdf");
  c2->SaveAs("pdf/drawPSWgts_Run2_top5.pdf");
}
