// Purpose: draw TT time stability
// 
#include "TFile.h"
#include "TProfile.h"
#include "TF1.h"
#include "TLatex.h"
#include "TLine.h"

#include <vector>

#include "../tdrstyle_mod15.C"

bool doBlind = true;

// helper function to retrieve profile and turn it into a histogram
TH1D *getIOV(string name, vector<TFile*> fs, bool scale = true,
	     double delta=0) {

  const char *cname = name.c_str();
  TH1D *h = new TH1D(Form("h%siov",cname),Form(";IOV;#LT%s#GT (%%)",cname),
		     11,0.5+delta,11.5+delta);
  
  h->GetXaxis()->SetBinLabel(1,"17B");
  h->GetXaxis()->SetBinLabel(2,"C");
  h->GetXaxis()->SetBinLabel(3,"D");
  h->GetXaxis()->SetBinLabel(4,"E");
  h->GetXaxis()->SetBinLabel(5,"F");
  h->GetXaxis()->SetBinLabel(6,"18A");
  h->GetXaxis()->SetBinLabel(7,"B");
  h->GetXaxis()->SetBinLabel(8,"C");
  h->GetXaxis()->SetBinLabel(9,"D");
  h->GetXaxis()->SetBinLabel(10,"DT");
  h->GetXaxis()->SetBinLabel(11,"MC");

  for (unsigned int ifile = 0; ifile != fs.size(); ++ifile) {
    TFile *f = fs[ifile];
    TProfile *p = (TProfile*)f->Get(Form("p%siov",cname)); assert(p);

    for (int i = 1; i != p->GetNbinsX()+1; ++i) {
      if (p->GetBinContent(i)!=0) {
	int k = (i>10 ? i-5 : i>1 ? i-4 : 11);
	assert(h->GetBinContent(k)==0);

	// Blind MC (for data/MC ratio) for now)
	if (k==11 && doBlind && (name=="mlb"||name=="mt"||name=="mt1"))
	  continue;

	h->SetBinContent(k, p->GetBinContent(i));
	h->SetBinError(k, p->GetBinError(i));
      }
    } // for i
    if (!scale) h->SetYTitle(p->GetYaxis()->GetTitle());
    if ( scale) h->SetYTitle("#LTobs#GT (%)");
  } // for ifile

  if (scale) {
    TF1 *f1 = new TF1("f1","[0]",0.5,9.5);
    h->Fit(f1,"QRN");
    double ref = f1->GetParameter(0);
    double eref = f1->GetParError(0);
    delete f1;

    h->SetBinContent(10,ref);
    h->SetBinError(10,eref);

    for (int i = 1; i != h->GetNbinsX()+1; ++i) {
      h->SetBinContent(i, 100. * (h->GetBinContent(i) / ref - 1));
      h->SetBinError(i, 100. * h->GetBinError(i) / ref);
    } // for i
  } // scale

  return h;
} // getIOV

void drawTTtime() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TFile *f17 = new TFile("rootfiles/hadWUL17V5_MPDGcorrNoW.root","READ");
  assert(f17 && !f17->IsZombie());
  TFile *f18 = new TFile("rootfiles/hadWUL18V5_MPDGcorrNoW.root","READ");
  assert(f18 && !f18->IsZombie());
  TFile *f1718 = new TFile("rootfiles/hadWUL1718V5_MPDGcorrNoW.root","READ");
  assert(f1718 && !f1718->IsZombie());
  TFile *fmc = new TFile("rootfiles/hadWMC18V5_MPDGcorrNoW.root","READ");
  assert(fmc && !fmc->IsZombie());

  vector<TFile*> fs;
  //fs.push_back(f17);
  //fs.push_back(f18);
  fs.push_back(f1718);
  fs.push_back(fmc);
  
  bool scale = true;
  TH1D *hmw = getIOV("mw",fs,scale,0);
  TH1D *hmt = getIOV("mt",fs,scale,+0.05);
  TH1D *hmt1 = getIOV("mt1",fs,scale,+0.10);
  TH1D *hrbq = getIOV("rbq",fs,scale,-0.05);
  TH1D *hmlb = getIOV("mlb",fs,scale,-0.10);

  //lumi_13TeV = "UL17+UL18, 101.4 fb^{-1}";
  lumi_13TeV = "UL17+UL18, 103.4 fb^{-1}";
  TCanvas *c1 = tdrCanvas("c1",hmw,4,11,kSquare);

  //hmw->Draw();
  hmw->SetMinimum(-0.6);//-0.3);
  hmw->SetMaximum(+1.2);//+0.6);

  TLine *l = new TLine();
  l->SetLineStyle(kDotted);
  l->DrawLine(0.5,0,11.5,0);

  tdrDraw(hmlb,"P",kOpenSquare,kRed);
  tdrDraw(hrbq,"P",kOpenDiamond,kGreen+2);
  tdrDraw(hmt,"P",kOpenCircle,kCyan+2);
  tdrDraw(hmt1,"P",kOpenCircle,kBlue);
  tdrDraw(hmw,"P",kFullCircle);

  hmt1->SetLineWidth(2);
  hmw->SetLineWidth(3);

  TLegend *leg = tdrLeg(0.70,0.90-5*0.05,1.00,0.90);
  leg->AddEntry(hmlb,"m_{lb}","PLE");
  leg->AddEntry(hrbq,"R_{bq}","PLE");
  leg->AddEntry(hmt,"m_{t,12}","PLE");
  leg->AddEntry(hmt1,"m_{t,1}","PLE");
  leg->AddEntry(hmw,"m_{W}","PLE");

  TLatex *tex = new TLatex();
  tex->SetTextSize(0.045); tex->SetNDC(true);

  tex->DrawLatex(0.45,0.85,"|#eta_{jets}| < 1.3");
  tex->DrawLatex(0.45,0.78,"fitProb > 0.2");

  tex->SetTextSize(0.035);
  tex->DrawLatex(0.45,0.65,Form("b/q: %+1.2f#pm%1.2f%%",
				hrbq->GetBinContent(10)-hrbq->GetBinContent(11),
				sqrt(pow(hrbq->GetBinError(10),2)+
				     pow(hrbq->GetBinError(11),2))));

  c1->SaveAs("pdf/drawTTtime.pdf");
} // drawTTtime
