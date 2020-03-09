// Purpose: extract b-JSF given Z+b and Z+jet data and MC
//            R_incl = R_b*p0 + R_qg*(1-p0),  p0~0.04
//            R_meas = R_b*p  + R_qg*(1-p),   p~0.85(t) - 0.55(m) - 0.25(l)
//            DeltaR = R_b - R_incl = (1-p0)/(p-p0)*(R_meas - R_incl)
//          and then data/MC for deltaR, with same p, p0 assumed (for now)
//
//
// Update: gluon fraction sizable for inclusive and response different from q,
//         charm fraction sizable for Z+b and balance different from b,
//         so account for these two effects to leading order
//            R_i = R_b*pb + R_c*pc + R_g*pg + R_q*pq  [pb+pc+pq+pg=1]
//            R_g = R_q * rg  [include FSR for ptchs; R_g later from data]
//            R_c = R_b * rc  [includes neutrinos in R_c and R_b]
//            R_m = R_b*p  + R_c*(1-p),   pc~(1-p)~0.07-0.15
//            R_b/R_q = R_m/R_i * (pq+rg*pg) / ((p+rc*(1-p))-R_m/R_i*(pb+rc*pc))
//
// To check inputs, use minitools/drawZplusB.C
//
// To-do:
// 0) Check compatibility with KIT Z+jet inputs
// 1) Update MC with ttbar background (for Z+b 75 GeV)
// 2) Update calculation with R_qg, p0 => R_uds, R_c, R_g, p0uds, p0c, p0g
//    2b) Need R_uds, R_c, R_g separately, or as ratio to inclusive R_incl
// 3) Add FSR weighted variants to constrain FSR
#include "TFile.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TLatex.h"

#include <string>

#include "../tdrstyle_mod15.C"

using namespace std;
int cnt(0);

void extractBJSFs(string sm="mpfchs", string sq="q", string sb="tight",
		  int a=30, double ptmin=30., double ptmax=400);
//250.);//1000.);

// Run all variants at once
void extractBJSF() {

  //extractBJSFs("mpfchs","deepflav",30);
  //extractBJSFs("ptchs","deepflav",30);
  extractBJSFs("rmpf","incl","deepb",30);
  extractBJSFs("ptchs","incl","deepb",30);
  extractBJSFs("rmpf","qtag","deepb",30);
  extractBJSFs("ptchs","qtag","deepb",30);
  //extractBJSFs("mpfchs","tight",30);
  //extractBJSFs("ptchs","tight",30);
  //extractBJSFs("mpfchs","medium",30);
  //extractBJSFs("ptchs","medium",30);
  //extractBJSFs("mpfchs","loose",30);
  //extractBJSFs("ptchs","loose",30);

  /*
  extractBJSFs("rmpf","qtag","deepb",20);
  extractBJSFs("ptchs","qtag","deepb",20);
  extractBJSFs("rmpf","qtag","deepb",40);
  extractBJSFs("ptchs","qtag","deepb",40);
  extractBJSFs("rmpf","qtag","deepb",50);
  extractBJSFs("ptchs","qtag","deepb",50);
  extractBJSFs("rmpf","qtag","deepb",60);
  extractBJSFs("ptchs","qtag","deepb",60);

  extractBJSFs("rmpf","qtag","deepb",100);
  extractBJSFs("ptchs","qtag","deepb",100);
  //extractBJSFs("mpfchs","deepb",100);
  //extractBJSFs("ptchs","deepb",100);
  //extractBJSFs("mpfchs","tight",100);
  //extractBJSFs("ptchs","tight",100);
  //extractBJSFs("mpfchs","medium",100);
  //extractBJSFs("ptchs","medium",100);
  //extractBJSFs("mpfchs","loose",100);
  //extractBJSFs("ptchs","loose",100);
  */
}


//void extractBJSFs() {
void extractBJSFs(string sm, string sq, string sb, int a, double ptmin, double ptmax) {

  setTDRStyle();

  TDirectory *curdir = gDirectory;

  // Input data from Sami
  //TFile *f = new TFile("rootfiles/jme_bplusZ_v8_2016GH.root","READ");
  //TFile *f = new TFile("rootfiles/jme_bplusZ_v10_2016GH.root","READ");
  //TFile *f = new TFile("rootfiles/jme_bplusZ_2016BCDEFGH_v11.root","READ");
  //TFile *f = new TFile("rootfiles/jme_bplusZ_2016GH_v12.root","READ");
  //TFile *f = new TFile("rootfiles/jme_bplusZ_Run2Muons_v14.root","READ");
  //TFile *f = new TFile("rootfiles/jme_bplusZ_Run2_14.root","READ");
  //TFile *f = new TFile("rootfiles/jme_bplusZ_Run2_v15.root","READ");
  //TFile *f = new TFile("rootfiles/jme_bplusZ_Run2_v16.root","READ");
  //TFile *f = new TFile("rootfiles/jme_bplusZ_Run2_Muons_v17.root","READ");
  //TFile *f = new TFile("rootfiles/jme_bplusZ_Run2_v18.root","READ");
  //TFile *f = new TFile("rootfiles/jme_bplusZ_Run2_v20.root","READ");
  //TFile *f = new TFile("rootfiles/jme_bplusZ_Run2_v20.root","READ");
  //TFile *f = new TFile("rootfiles/jme_bplusZ_Run2_v20_noTTBar.root","READ");
  //TFile *f = new TFile("rootfiles/jme_bplusZ_Run2_Ver2.root","READ"); // tt?
  //TFile *f = new TFile("rootfiles/jme_bplusZ_Run2_Ver3.root","READ"); // tt?
  //TFile *f = new TFile("rootfiles/jme_bplusZ_Run2_Ver4.root","READ"); // tt?
  //TFile *f = new TFile("rootfiles/jme_bplusZ_Run2_Ver5.root","READ"); // tt?
  TFile *f = new TFile("rootfiles/jme_bplusZ_Run2_Ver6.root","READ"); // tt?
  assert(f && !f->IsZombie());

  // Markus Seidel, 12 Feb 2014
  TFile *fj = new TFile("rootfiles/jetPtEta-Markus.root","READ");
  assert(fj && !fj->IsZombie());

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
  string ssq, ssb, ssb2; // long name
  if (sq=="incl") ssq = "";
  if (sq=="qtag") ssq = "_quarktag";
  //if (sb=="deepb") ssb = "_btagDeepBtight";
  if (sb=="deepb") ssb = "_btagdeepbtight";
  if (sb=="deepb" && sm=="ptchs") ssb = "_btagDeepBtight";
  if (sb=="deepb") ssb2 = "_btagDeepBtight";
  if (sb=="deepflav") ssb = "_btagDeepFlavBtight";
  if (sb=="tight") ssb = "_btagCSVV2tight";
  if (sb=="medium") ssb = "_btagCSVV2medium";
  if (sb=="loose") ssb = "_btagCSVV2loose";
  const char *cq = ssq.c_str();
  const char *cb = ssb.c_str();
  const char *cb2 = ssb2.c_str();
  //string seta = "eta13";
  //string seta = "eta_00_13";
  string seta = "eta_00_25";
  const char *ceta = seta.c_str();
  const char *cm = sm.c_str();

  // ATLAS B,D hadron fractions and muon branching fractions (PDG and P8)
  // Use these to estimate neutrino fraction
  // ATLAS-CONF-2019-046
  //double rb = 0.85; // <pTB/ptjet>
  //const int nb = 4; // B0, B+, B0s, bb
  //doubld vb0[nb] = {0.404, 0.404, 0.103+0.001, 0.088};
  //const int nc; // D+, D0, D0s, cb
  //double vc0[nc] = {0.226, 0.564, 0.080, 0.109}; // 0.021 left (uds?)
  //double nmu = 5; // bmu, btau, bcmu, bcbarmu, cmu
  //double vmu[nmu] = {0.1095, 0.0042, 0.0802, 0.016, 0.082};
  //double fracnu = 2*0.85*(vmu[0]
  // ATLAS bJES and neutrino fraction refrence at the bottom here
  // https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/PERF-2012-01/

  // Load R_incl, R_meas, purities
  TH1D *hds = (TH1D*)f->Get(Form("data/%s/statistics_zmmjet%s_a%d",ceta,cb,a));
  assert(hds);
  TGraphErrors *gdri = (TGraphErrors*)f->Get(Form("data/%s/%s_zmmjet%s_a%d",ceta,cm,cq,a));
  assert(gdri);
  TGraphErrors *gdrb = (TGraphErrors*)f->Get(Form("data/%s/%s_zmmjet%s_a%d",ceta,cm,cb,a));
  assert(gdrb);
  TGraphErrors *gmri = (TGraphErrors*)f->Get(Form("mc/%s/%s_zmmjet%s_a%d",ceta,cm,cq,a));
  assert(gmri);
  TGraphErrors *gmrb = (TGraphErrors*)f->Get(Form("mc/%s/%s_zmmjet%s_a%d",ceta,cm,cb,a));
  assert(gmrb);
  TGraphErrors *gdpi = (TGraphErrors*)f->Get(Form("mc/%s/bpurity_zmmjet%s_a%d",ceta,cq,a));
  assert(gdpi);
  TGraphErrors *gdpb = (TGraphErrors*)f->Get(Form("mc/%s/bpurity_zmmjet%s_a%d",ceta,cb2,a));
  assert(gdpb);

  // Get b-jet pT distribution at |eta|<1.3 in 10 GeV bins
  TH2F *h2j = (TH2F*)fj->Get("hBTagged"); assert(h2j);
  int i1 = h2j->GetXaxis()->FindBin(-1.3+0.05);
  int i2 = h2j->GetXaxis()->FindBin(+1.3-0.05);
  TH1D *hj1 = h2j->ProjectionY("_px",i1,i2);
  // Remap to Z+jet pT bins

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


  TH1D *hup = new TH1D(Form("hup%d",++cnt),
		       ";p_{T,Z} (GeV);Data - MC (%)",1470,30,1500);
  hup->GetXaxis()->SetMoreLogLabels();
  hup->GetXaxis()->SetNoExponent();
  hup->SetMinimum(-7);
  hup->SetMaximum(14);

  TH1D *hdw = new TH1D(Form("hdw%d",++cnt),
		       ";p_{T,Z} (GeV);Z+b - Z+j (%)",1470,30,1500);
  hdw->GetXaxis()->SetMoreLogLabels();
  hdw->GetXaxis()->SetNoExponent();
  hdw->SetMinimum(-15);//-10);
  hdw->SetMaximum(+5);
  if (sq=="incl") hdw->SetYTitle("Z+b - Z+j (%)");
  if (sq=="qtag") hdw->SetYTitle("Z+b - Z+q (%)");

  //lumi_13TeV = "2016GH, 16.2 fb^{-1}";
  //lumi_13TeV = "2016BCDEFGH, 36.5 fb^{-1}";
  lumi_13TeV = "Run2, 137.5 fb^{-1}";
  TCanvas *c1 = tdrDiCanvas(Form("c1_%d",++cnt),hup,hdw,4,11);

  TLatex *tex = new TLatex();
  tex->SetTextSize(0.045); tex->SetNDC();
  
  c1->cd(1);
  gPad->SetLogx();
  tdrDraw(gr,"Pz",kFullCircle);
  
  // Alternative shape
  //TF1 *frb = new TF1("fr","[0]+[1]*TMath::Gaus(x,[2],[3])",ptmin,ptmax);
  //frb->SetParameters(0,0.5,75,10); frb->FixParameter(2,75.);
  TF1 *frb = new TF1("fr","[0]+[1]*TMath::Gaus(log(x/[2]),0,log(1+[3]/[2]))",ptmin,ptmax);
  frb->SetParameters(0,0.5,75.,15.); frb->FixParameter(2,75.);
  frb->SetParLimits(3,7.5,22.5);
  frb->SetLineColor(kBlue-9);
  gr->Fit(frb,"QRN");
  frb->Draw("SAME");

  TF1 *fr = new TF1("fr","[0]",ptmin,ptmax);
  fr->SetLineColor(kBlack);
  gr->Fit(fr,"QRN");
  fr->Draw("SAME");
  tex->DrawLatex(0.55,0.80,Form("#DeltaR_{b} = %1.2f #pm %1.2f (%%)",
				fr->GetParameter(0),fr->GetParError(0)));
  tex->DrawLatex(0.55,0.74,Form("#chi^{2} / NDF = %1.1f / %d",
				fr->GetChisquare(), fr->GetNDF()));

  tex->DrawLatex(0.20,0.72,Form("%s",sm.c_str()));
  if (seta=="eta13" || seta=="eta_00_13")
    tex->DrawLatex(0.20,0.67,Form("#alpha<%1.2f, |#eta|<1.3",a*0.01));
  else if (seta=="eta_00_25")
    tex->DrawLatex(0.20,0.67,Form("#alpha<%1.2f, |#eta|<2.5",a*0.01));
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

  //TLegend *legdw = tdrLeg(0.20,0.65,0.40,0.90);
  TLegend *legdw = tdrLeg(0.40,0.75,0.70,0.90);
  legdw->SetNColumns(2);
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

  //c1->SaveAs(Form("pdf/extractBJSF_%s_%s_a%d.pdf",sm.c_str(),sb.c_str(),a));
  c1->SaveAs(Form("pdf/extractBJSF_%s_%s_%s_a%d.pdf",
		  cm,sq.c_str(),sb.c_str(),a));
} // void extractbJSF
