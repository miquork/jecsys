// Purpose: plots mw, mt etc. vs mu to check PU stability
#include "TFile.h"
#include "TProfile.h"
#include "TLine.h"
#include "TF1.h"

#include "../tdrstyle_mod15.C"

void hadW_mu() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TFile *fd = new TFile("rootfiles/hadWUL161718V7_JEC.root","READ");
  assert(fd && !fd->IsZombie());
  
  TFile *fm = new TFile("rootfiles/hadWMC161718V7_JEC.root","READ");
  assert(fm && !fm->IsZombie());
  
  curdir->cd();

  // Select profiles to be plotted
  string ah[] = {"mw","mt","mlb","mlbmet","rbq"};
  const int nh = sizeof(ah)/sizeof(ah[0]);
  map<string, TH1D*> mhd;
  map<string, TH1D*> mhm;

  // Set marker, color, size and label
  map<string,int> color;
  color["mw"] =     kGreen+2;
  color["mt"] =     kBlue;
  color["mlb"] =    kCyan+2;
  color["mlbmet"] = kCyan+4;
  color["rbq"] =    kRed;
  map<string,int> markerd;
  markerd["mw"] =     kFullCircle;
  markerd["mt"] =     kFullSquare;
  markerd["mlb"] =    kFullSquare;
  markerd["mlbmet"] = kFullSquare;
  markerd["rbq"] =    kFullDiamond;
  map<string,int> markerm;
  markerm["mw"] =     kOpenCircle;
  markerm["mt"] =     kOpenSquare;
  markerm["mlb"] =    kOpenSquare;
  markerm["mlbmet"] = kOpenSquare;
  markerm["rbq"] =    kOpenDiamond;
  map<string,double> size;
  size["mw"] =     1.0;
  size["mt"] =     1.0;
  size["mlb"] =    0.8;
  size["mlbmet"] = 0.6;
  size["rbq"] =    1.0;
  map<string,const char*> label;
  label["mw"] =     "m_{qq}";
  label["mt"] =     "m_{bqq}";
  label["mlb"] =    "m_{lb}";
  label["mlbmet"] = "m_{lbmet}";
  label["rbq"] =    "R_{bq}";

  // Load profiles
  for (int i = 0; i != nh; ++i) {
  
    const char *c = ah[i].c_str();
    TProfile *pd = (TProfile*)fd->Get(Form("p%svsmu",c));
    assert(pd);
    TH1D *hd = pd->ProjectionX(Form("h%sd",c));

    TProfile *pm = (TProfile*)fm->Get(Form("p%svsmu",c));
    assert(pm);
    TH1D *hm = pm->ProjectionX(Form("h%sm",c));

    hm->Scale(1./pd->GetMean(2));
    hd->Scale(1./pd->GetMean(2));

    mhd[c] = hd;
    mhm[c] = hm;
  }
  
  //TH1D *h = tdrHist("h","Ratio to mean",0.980,1.020,"#mu",0,60);
  TH1D *h = tdrHist("h","Ratio to mean",0.950,1.050,"#mu",0,60);
  
  lumi_13TeV = "UL16V7 + UL17V5 + UL18V5";
  TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);
  gPad->SetLeftMargin(0.17);
  h->SetTitleOffset(1.5,"Y");

  TLegend *leg = tdrLeg(0.50,0.9-mhd.size()*0.05,0.70,0.90);
  
  // Draw histograms, determine slopes vs mu
  double muref = 30;
  TH1D *h2 =  tdrHist("h2", "Slope vs #mu (%)",-1,5,"Observable",-0.5,nh-0.5);
  TH1D *hsd = tdrHist("hsd","Slope vs #mu (%)",-1,5,"Observable",-0.5,nh-0.5);
  TH1D *hsm = tdrHist("hsm","Slope vs #mu (%)",-1,5,"Observable",-0.5,nh-0.5);

  for (int i = 0; i != nh; ++i) {
    
    const char *c = ah[i].c_str();
    TH1D *hd = mhd[c];
    TH1D *hm = mhm[c];
    tdrDraw(hd,"Pz",markerd[c],color[c]); hd->SetMarkerSize(size[c]);
    tdrDraw(hm,"Pz",markerm[c],color[c]); hm->SetMarkerSize(size[c]);
    leg->AddEntry(hd,label[c],"PLE");
    
    TF1 *f1d = new TF1(Form("f1d_%s",c),"[0]+[1]*x",0,60);
    hd->Fit(f1d,"QRN");
    f1d->SetLineColor(color[c]);
    f1d->Draw("SAME");

    h2->GetXaxis()->SetBinLabel(i+1, label[c]);
    hsd->SetBinContent(i+1, 100.*muref*f1d->GetParameter(1));
    hsd->SetBinError(i+1, 100.*muref*f1d->GetParError(1));

    TF1 *f1m = new TF1(Form("f1m_%s",c),"[0]+[1]*x",0,60);
    hm->Fit(f1m,"QRN");
    f1m->SetLineColor(color[c]);
    f1m->SetLineStyle(kDotted);
    f1m->Draw("SAME");

    hsm->SetBinContent(i+1, 100.*muref*f1m->GetParameter(1));
    hsm->SetBinError(i+1, 100.*muref*f1m->GetParError(1));
  } // for ih

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(0,1,60,1);
  l->SetLineStyle(kDotted);

  TLatex *tex = new TLatex();
  tex->SetNDC(); 
  tex->SetTextSize(0.045);

  gPad->Update();

  c1->SaveAs("pdf/hadW_mu.pdf");

  // Create plot to show slopes in data and MC in single plot with uncertainty
  TCanvas *c2 = tdrCanvas("c2",h2,4,11,kSquare);

  l->SetLineStyle(kDashed);
  l->DrawLine(-0.5,0,nh-0.5,0);

  tdrDraw(hsd,"Pz",kFullCircle);
  tdrDraw(hsm,"Pz",kOpenCircle);

  c2->SaveAs("pdf/hadW_mu_slopes.pdf");
} // hadW_mu
