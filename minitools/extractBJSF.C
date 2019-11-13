
// Purpose: extract b-JSF given Z+b and Z+jet data and MC
//            R_incl = R_b*p0 + R_qg*(1-p0),  p0~0.04
//            R_meas = R_b*p  + R_qg*(1-p),   p~0.85(t) - 0.55(m) - 0.25(l)
//            DeltaR = R_b - R_incl = (1-p0)/(p-p0)*(R_meas - R_incl)
//          and then data/MC for deltaR, with same p, p0 assumed (for now)
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TLatex.h"

#include <string>

#include "../tdrstyle_mod15.C"

using namespace std;

void extractBJSFs(string sm="mpfchs", string sb="tight", int a=30,
		  double ptmin=30., double ptmax=250.);//1000.);

// Run all variants at once
void extractBJSF() {

  extractBJSFs("mpfchs","tight",30);
  extractBJSFs("ptchs","tight",30);
  extractBJSFs("mpfchs","medium",30);
  extractBJSFs("ptchs","medium",30);
  extractBJSFs("mpfchs","loose",30);
  extractBJSFs("ptchs","loose",30);

  extractBJSFs("mpfchs","tight",100);
  extractBJSFs("ptchs","tight",100);
  extractBJSFs("mpfchs","medium",100);
  extractBJSFs("ptchs","medium",100);
  extractBJSFs("mpfchs","loose",100);
  extractBJSFs("ptchs","loose",100);
}


//void extractBJSFs() {
void extractBJSFs(string sm, string sb, int a, double ptmin, double ptmax) {

  setTDRStyle();

  TDirectory *curdir = gDirectory;

  // Input data from Sami
  //TFile *f = new TFile("rootfiles/jme_bplusZ_v8_2016GH.root","READ");
  //TFile *f = new TFile("rootfiles/jme_bplusZ_v10_2016GH.root","READ");
  //TFile *f = new TFile("rootfiles/jme_bplusZ_2016BCDEFGH_v11.root","READ");
  //TFile *f = new TFile("rootfiles/jme_bplusZ_2016GH_v12.root","READ");
  //TFile *f = new TFile("rootfiles/jme_bplusZ_Run2Muons_v14.root","READ");
  TFile *f = new TFile("rootfiles/jme_bplusZ_Run2_14.root","READ");
  assert(f && !f->IsZombie());

  curdir->cd();

  //double ptmin = 30;//50; //30;
  //double ptmax = 1000.;//500;//700.; //1000; //700;

  // Reference alpha value (loop later?)
  //int a = 100;
  //int a = 30;

  // Reference method (mpfchs)
  //string sm = "mpfchs"; // method: mpfchs, ptchs
  //string sm = "ptchs"; // method: mpfchs, ptchs

  // Reference b-tagger (loop later?)
  //string sb = "tight"; // short name: loose, medium, tight
  //string sb = "medium"; // short name: loose, medium, tight
  string ssb; // long name
  if (sb=="tight") ssb = "_btagCSVV2tight";
  if (sb=="medium") ssb = "_btagCSVV2medium";
  if (sb=="loose") ssb = "_btagCSVV2loose";
  const char *cb = ssb.c_str();
  string seta = "eta13";
  const char *ceta = seta.c_str();

  // Load R_incl, R_meas, purities
  TH1D *hds = (TH1D*)f->Get(Form("data/%s/statistics_zmmjet%s_a%d",ceta,cb,a));
  assert(hds);
  TGraphErrors *gdri = (TGraphErrors*)f->Get(Form("data/%s/%s_zmmjet_a%d",ceta,sm.c_str(),a));
  assert(gdri);
  TGraphErrors *gdrb = (TGraphErrors*)f->Get(Form("data/%s/%s_zmmjet%s_a%d",ceta,sm.c_str(),cb,a));
  assert(gdrb);
  TGraphErrors *gmri = (TGraphErrors*)f->Get(Form("mc/%s/%s_zmmjet_a%d",ceta,sm.c_str(),a));
  assert(gmri);
  TGraphErrors *gmrb = (TGraphErrors*)f->Get(Form("mc/%s/%s_zmmjet%s_a%d",ceta,sm.c_str(),cb,a));
  assert(gmrb);
  TGraphErrors *gdpi = (TGraphErrors*)f->Get(Form("mc/%s/bpurity_zmmjet_a%d",ceta,a));
  assert(gdpi);
  TGraphErrors *gdpb = (TGraphErrors*)f->Get(Form("mc/%s/bpurity_zmmjet%s_a%d",ceta,cb,a));
  assert(gdpb);

  curdir->cd();

  TGraphErrors *gdr = new TGraphErrors(gdrb->GetN());
  TGraphErrors *gmr = new TGraphErrors(gmrb->GetN());
  TGraphErrors *gr = new TGraphErrors(gmrb->GetN());
  for (int i = 0; i != gdrb->GetN(); ++i) {
    
    double pt = gdrb->GetX()[i];
    int ibin = hds->FindBin(pt);
    int nbin = hds->GetBinContent(ibin);
    double dpt = (nbin ? 0.5*hds->GetBinWidth(ibin)/sqrt(nbin) : 0);

    // Need to do some more checks that pT always the same
    double p0 = gdpi->GetY()[i];
    double p = gdpb->GetY()[i];
    //
    double Ri = gdri->GetY()[i];
    double Rb = gdrb->GetY()[i];
    double dRb = gdrb->GetEY()[i];
    double dRd = (p!=p0 ? (1-p0)/(p-p0)*(Rb - Ri) : 0);
    double ddRd = (p!=p0 ? (1-p0)/(p-p0)*dRb : 0);
    //
    double Rim = gmri->GetY()[i];
    double Rbm = gmrb->GetY()[i];
    double dRbm = gmrb->GetEY()[i];
    double dRm = (p!=p0 ? (1-p0)/(p-p0)*(Rbm - Rim) : 0);
    double ddRm = (p!=p0 ? (1-p0)/(p-p0)*dRbm : 0);
    //
    double dR = dRd-dRm;
    double ddR = sqrt(pow(ddRd,2)+pow(ddRm,2));

    gdr->SetPoint(i, pt, 100.*dRd);
    gdr->SetPointError(i, dpt, 100.*ddRd);
    gmr->SetPoint(i, pt, 100.*dRm);
    gmr->SetPointError(i, dpt, 100.*ddRm);
    gr->SetPoint(i, pt, 100.*dR);
    gr->SetPointError(i, dpt, 100.*ddR);
  } // for i

  /*
  gdr->SetMarkerStyle(kFullCircle);
  //gdr->Draw("AP");
  gmr->SetMarkerStyle(kOpenCircle);
  //gmr->Draw("SAMEP");

  gr->SetMarkerStyle(kFullCircle);
  gr->Draw("AP");*/


  TH1D *hup = new TH1D("hup",";p_{T,Z} (GeV);Data - MC (%)",1470,30,1500);
  hup->GetXaxis()->SetMoreLogLabels();
  hup->GetXaxis()->SetNoExponent();
  hup->SetMinimum(-7);
  hup->SetMaximum(14);

  TH1D *hdw = new TH1D("hdw",";p_{T,Z} (GeV);Z+b - Z+j (%)",1470,30,1500);
  hdw->GetXaxis()->SetMoreLogLabels();
  hdw->GetXaxis()->SetNoExponent();
  hdw->SetMinimum(-10);
  hdw->SetMaximum(+5);

  //lumi_13TeV = "2016GH, 16.2 fb^{-1}";
  //lumi_13TeV = "2016BCDEFGH, 36.5 fb^{-1}";
  lumi_13TeV = "Run2, 137.5 fb^{-1}";
  TCanvas *c1 = tdrDiCanvas("c1",hup,hdw,4,11);

  TLatex *tex = new TLatex();
  tex->SetTextSize(0.045); tex->SetNDC();
  
  c1->cd(1);
  gPad->SetLogx();
  tdrDraw(gr,"Pz",kFullCircle);
  
  TF1 *fr = new TF1("fr","[0]",ptmin,ptmax);
  fr->SetLineColor(kBlack);
  gr->Fit(fr,"QRN");
  fr->Draw("SAME");
  tex->DrawLatex(0.55,0.80,Form("#DeltaR_{b} = %1.2f #pm %1.2f (%%)",
				fr->GetParameter(0),fr->GetParError(0)));
  tex->DrawLatex(0.55,0.74,Form("#chi^{2} / NDF = %1.1f / %d",
				fr->GetChisquare(), fr->GetNDF()));

  tex->DrawLatex(0.20,0.72,Form("%s",sm.c_str()));
  if (seta=="eta13")
    tex->DrawLatex(0.20,0.67,Form("#alpha<%1.2f, |#eta|<1.3",a*0.01));
  else 
    tex->DrawLatex(0.20,0.67,Form("#alpha<%1.2f, other eta",a*0.01));
  tex->DrawLatex(0.20,0.62,Form("%s b-id",sb.c_str()));

  c1->cd(2);
  gPad->SetLogx();
  tdrDraw(gdr,"P",kFullCircle);
  tdrDraw(gmr,"P",kOpenCircle);
  
  TF1 *fd = new TF1("fd","[0]",ptmin,ptmax);
  fd->SetLineColor(kBlack);
  gdr->Fit(fd,"QRN");
  fd->Draw("SAME");
  TF1 *fm = new TF1("fm","[0]",30,ptmax);
  fm->SetLineColor(kBlack);
  fm->SetLineStyle(kDashed);
  gmr->Fit(fm,"QRN");
  fm->Draw("SAME");

  TLegend *legdw = tdrLeg(0.20,0.65,0.40,0.90);
  legdw->SetTextSize(0.045*2.2);
  legdw->AddEntry(gdr,"Data","PL");
  legdw->AddEntry(gmr,"MC","PL");

  tex->DrawLatex(0.78,0.66,Form("R_{b,data} = %1.2f#pm%1.2f (%%)",
				fd->GetParameter(0),fd->GetParError(0)));
  tex->DrawLatex(0.78,0.60,Form("#chi^{2} / NDF = %1.1f / %d",
				fd->GetChisquare(), fd->GetNDF()));
  tex->DrawLatex(0.78,0.50,Form("R_{b,MC} = %1.2f#pm%1.2f (%%)",
				fm->GetParameter(0),fm->GetParError(0)));
  tex->DrawLatex(0.78,0.44,Form("#chi^{2} / NDF = %1.1f / %d",
				fm->GetChisquare(), fm->GetNDF()));

  c1->SaveAs(Form("pdf/extractBJSF_%s_%s_a%d.pdf",sm.c_str(),sb.c_str(),a));
} // void extractbJSF
