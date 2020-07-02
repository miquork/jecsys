// Purpose: study Z+jet for jet flavor (b, c, uds, g)
//          draw <pT,reco>/<pT,gen> vs pTgen(?)
//
#include "TFile.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TLine.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMultiGraph.h"

#include <map>
#include <string>
#include <iostream>

#include "../tdrstyle_mod15.C"
#include "../tools.C"

using namespace std;

void addOff(TGraphErrors *g, double doff=-0.734) {
  for (int i = 0; i != g->GetN(); ++i) {
    g->SetPoint(i, g->GetX()[i], g->GetY()[i]+doff/g->GetX()[i]);
  } // for i
} // addOff

void subtractGraph(TGraphErrors *g1, TGraphErrors *g2,
		   double k1=1, double k2=1) {

  assert(g1->GetN()==g2->GetN());
  for (int i = 0; i != g1->GetN(); ++i) {
    double x = g1->GetX()[i];
    assert(g2->GetX()[i]==x);
    g1->SetPoint(i, x, k1*g1->GetY()[i]-k2*g2->GetY()[i]);
    double ey1 = k1*g1->GetEY()[i];
    double ey2 = k2*g2->GetEY()[i];
    g1->SetPointError(i, g1->GetEX()[i], sqrt(ey1*ey1+ey2*ey2));
  }
} // subtractGraph

void subtractGraphFit(TGraphErrors *g1, TGraphErrors *g2,
		   double k1=1, double k2=1) {

  assert(g1->GetN()==g2->GetN());
  for (int i = 0; i != g1->GetN(); ++i) {
    double x = g1->GetX()[i];
    //assert(g2->GetX()[i]==x);
    double y2 = g2->Eval(x);
    g1->SetPoint(i, x, k1*g1->GetY()[i]-k2*y2);
    double ey1 = k1*g1->GetEY()[i];
    double ey2 = k2*g2->GetEY()[i];
    g1->SetPointError(i, g1->GetEX()[i], sqrt(ey1*ey1+ey2*ey2));
  }
} // subtractGraph

// Map <pTreco/pTgen> vs pTZ (ptchs) and <pTgen/pTZ> vs pTZ (bal) to
// <pTreco>/<pTgen> vs <pTgen>
TGraphErrors* mapZ(TGraphErrors *gr, TGraphErrors *gg) {
  assert(gr->GetN()==gg->GetN());

  TGraphErrors *g = (TGraphErrors*)gr->Clone(Form("g_%s",gr->GetName()));
  for (int i = 0; i != g->GetN(); ++i) {
    g->SetPoint(i, g->GetX()[i]*gg->GetY()[i], g->GetY()[i]/gg->GetY()[i]);
  } // for i

  return g;
} // mapZ

void drawZtruths(unsigned long int bitmap = 1023L) {
  
  // bitmap encodes which variants are to be drawn
  // shifts 0-4: pTgen (a,b,c,q,g)
  // shifts 5-9: pTZ (a,b,c,q,g)
  // a,q,g both = 825L
  // q,g both = 792L

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  //TFile *f = new TFile("rootfiles/jme_bplusZ_Run2_Ver6.root","READ");
  //TFile *f = new TFile("rootfiles/jme_bplusZ_merged_v11.root","READ");
  //TFile *f = new TFile("rootfiles/jme_bplusZ_merged_v12.root","READ");
  //TFile *f = new TFile("rootfiles/jme_bplusZ_merged_v13.root","READ");
  //TFile *f = new TFile("rootfiles/jme_bplusZ_merged_vX.root","READ");
  //TFile *f = new TFile("rootfiles/jme_bplusZ_merged_v24.root","READ");
  TFile *f = new TFile("rootfiles/jme_bplusZ_merged_v25.root","READ");
  assert(f && !f->IsZombie());

  // jet1=response?
  // ptchs vs genjet1? what is bal?
  // mpfchs vs gmpf? and vs 1-guncl

  // (mpfjet1, mpfjetn, mpfuncl) vs (genjet1, genjetn, genuncl)

  //why didn't this work:
  //f->cd("profile_Run2017");
  //TDirectory *d = curdir;
  //cout << d->GetName() << endl << flush;

  const char *cd = "mc/eta_00_13";
  // Response <pTreco/pTgen> vs pTgen 
  TGraphErrors *pag = (TGraphErrors*)f->Get(Form("%s/response_ptgen_zmmjet_a30",cd));
  TGraphErrors *pbg = (TGraphErrors*)f->Get(Form("%s_genb/response_ptgen_zmmjet_a30",cd));
  TGraphErrors *pcg = (TGraphErrors*)f->Get(Form("%s_genc/response_ptgen_zmmjet_a30",cd));
  TGraphErrors *pqg = (TGraphErrors*)f->Get(Form("%s_genuds/response_ptgen_zmmjet_a30",cd));
  TGraphErrors *pgg = (TGraphErrors*)f->Get(Form("%s_geng/response_ptgen_zmmjet_a30",cd));

  assert(pag);
  assert(pbg);
  assert(pcg);
  assert(pqg);
  assert(pgg);

  // Response <pTreco/pTgen> vs pTZ (biased by nominator and x-axis)
  // => also variant from ptchs / bal
  TGraphErrors *paz = (TGraphErrors*)f->Get(Form("%s/response_zmmjet_a30",cd));
  TGraphErrors *pbz = (TGraphErrors*)f->Get(Form("%s_genb/response_zmmjet_a30",cd));
  TGraphErrors *pcz = (TGraphErrors*)f->Get(Form("%s_genc/response_zmmjet_a30",cd));
  TGraphErrors *pqz = (TGraphErrors*)f->Get(Form("%s_genuds/response_zmmjet_a30",cd));
  TGraphErrors *pgz = (TGraphErrors*)f->Get(Form("%s_geng/response_zmmjet_a30",cd));

  assert(paz);
  assert(pbz);
  assert(pcz);
  assert(pqz);
  assert(pgz);

  // ptchs: <pTreco/pTZ> vs pTZ
  // bal: <pTgen/pTZ> vs pTZ
  // => ptchs / bal vs x*bal = <pTreco>/<pTgen> vs <pTgen>
  // TMP PATCH: ptchs => mpfjet1
  TGraphErrors *pazr = (TGraphErrors*)f->Get(Form("%s/mpfjet1_zmmjet_a30",cd));
  TGraphErrors *pbzr = (TGraphErrors*)f->Get(Form("%s_genb/mpfjet1_zmmjet_a30",cd));
  TGraphErrors *pczr = (TGraphErrors*)f->Get(Form("%s_genc/mpfjet1_zmmjet_a30",cd));
  TGraphErrors *pqzr = (TGraphErrors*)f->Get(Form("%s_genuds/mpfjet1_zmmjet_a30",cd));
  TGraphErrors *pgzr = (TGraphErrors*)f->Get(Form("%s_geng/mpfjet1_zmmjet_a30",cd));
  //
  // TMP PATCH: bal => genjet1
  TGraphErrors *pazg = (TGraphErrors*)f->Get(Form("%s/genjet1_zmmjet_a30",cd));
  TGraphErrors *pbzg = (TGraphErrors*)f->Get(Form("%s_genb/genjet1_zmmjet_a30",cd));
  TGraphErrors *pczg = (TGraphErrors*)f->Get(Form("%s_genc/genjet1_zmmjet_a30",cd));
  TGraphErrors *pqzg = (TGraphErrors*)f->Get(Form("%s_genuds/genjet1_zmmjet_a30",cd));
  TGraphErrors *pgzg = (TGraphErrors*)f->Get(Form("%s_geng/genjet1_zmmjet_a30",cd));

  assert(pazr);
  assert(pbzr);
  assert(pczr);
  assert(pqzr);
  assert(pgzr);
  //
  assert(pazg);
  assert(pbzg);
  assert(pczg);
  assert(pqzg);
  assert(pgzg);

  //const double off = -0.734;
  //const double off = -0.734 / 0.92;
  //const double off = -0.5 / 0.92; // PF-like
  //const double off = -0.5 / 0.92 * 0.5; // PFchs-like
  //const double off = -0.2 / 0.92 * 0.5; // PFchs-like
  const double off = 0; // dijet-like
  if (off!=0) {
    addOff(pag,off);
    addOff(pbg,off);
    addOff(pcg,off);
    addOff(pqg,off);
    addOff(pgg,off);
    //
    addOff(paz,off);
    addOff(pbz,off);
    addOff(pcz,off);
    addOff(pqz,off);
    addOff(pgz,off);
    //
    addOff(pazr,off);
    addOff(pbzr,off);
    addOff(pczr,off);
    addOff(pqzr,off);
    addOff(pgzr,off);
  }

  TGraphErrors *paz2 = mapZ(pazr,pazg);
  TGraphErrors *pbz2 = mapZ(pbzr,pbzg);
  TGraphErrors *pcz2 = mapZ(pczr,pczg);
  TGraphErrors *pqz2 = mapZ(pqzr,pqzg);
  TGraphErrors *pgz2 = mapZ(pgzr,pgzg);

  curdir->cd();

  const double ptmin = 25;
  const double ptmax = 1500;//400;
  TH1D *hup = new TH1D("hup",
		       ";#LTp_{T,gen}#GT (GeV);"
		       "#LTp_{T,reco}#GT / #LTp_{T,gen}#GT",
		       //";p_{T,Z} or p_{T,gen} (GeV);"
		       //"#LTp_{T,reco}/p_{T,gen}#GT",
		       int(ptmax-ptmin),ptmin,ptmax);
  hup->SetMinimum(0.975+1e-4);
  hup->SetMaximum(1.115-1e-4);
  hup->GetXaxis()->SetMoreLogLabels();
  hup->GetXaxis()->SetNoExponent();

  TH1D *hdw = new TH1D("hdw",
		       ";#LTp_{T,gen}#GT (GeV);"
		       //";p_{T,Z} or p_{T,gen} (GeV);"
		       "R_{f} - R_{q} (%)",
		       int(ptmax-ptmin),ptmin,ptmax);
  hdw->SetMinimum(-5.5+1e-4);
  hdw->SetMaximum(+0.5-1e-4);
  hdw->GetXaxis()->SetMoreLogLabels();
  hdw->GetXaxis()->SetNoExponent();

  lumi_13TeV = "Run2, 136.5 fb^{-1}";
  TCanvas *c1 = tdrDiCanvas("c1",hup,hdw,4,11);
  hup->GetYaxis()->SetTitleOffset(1.20);
  hdw->GetYaxis()->SetTitleOffset(0.50);

  c1->cd(1);
  gPad->SetLeftMargin(0.18);
  //gPad->SetRightMargin(0.02);
  gPad->SetLogx();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(ptmin,1,ptmax,1);

  // Add dijet to the back
  bool addDijet = true;
  if (addDijet) {
    TFile *fd = new TFile("rootfiles/output-MC-1-UL17V4_BCDEF.root","READ");
    assert(fd && !fd->IsZombie());
    fd->cd("Standard/Eta_0.0-1.3/mc");
    TProfile *pa = (TProfile*)gDirectory->Get("p2r_g"); assert(pa);
    TProfile *pq = (TProfile*)gDirectory->Get("p2r_q_g"); assert(pq);
    TProfile *pg = (TProfile*)gDirectory->Get("p2r_g_g"); assert(pg);
    TH1D *hrq = pq->ProjectionX("hrq");
    hrq->Add(pq,pq,100,-100);
    TH1D *hrg = pg->ProjectionX("hrg");
    hrg->Add(pg,pq,100,-100);
    
    c1->cd(2);
    hrq->GetXaxis()->SetRangeUser(25,1500);
    hrg->GetXaxis()->SetRangeUser(25,1500);
    tdrDraw(hrq,"Pz",kFullCircle,kMagenta-9); hrq->SetMarkerSize(0.8);
    tdrDraw(hrg,"Pz",kFullCircle,kBlue-9); hrg->SetMarkerSize(0.8);

    c1->cd(1);

    tdrDraw(pa,"Pz",kFullCircle,kGray); pa->SetMarkerSize(0.8);
    tdrDraw(pq,"Pz",kFullCircle,kMagenta-9); pq->SetMarkerSize(0.8);
    tdrDraw(pg,"Pz",kFullCircle,kBlue-9); pg->SetMarkerSize(0.8);

    TLegend *legd = tdrLeg(0.70,0.32,1.00,0.52);
    legd->SetHeader("Dijet");
    legd->AddEntry(pa,"All","PL");
    legd->AddEntry(pq,"q","PL");
    legd->AddEntry(pg,"g","PL");
  }

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  tex->DrawLatex(0.22,0.70,"|#eta| < 1.3");
  if (off!=0 )tex->DrawLatex(0.22,0.64,Form("#DeltaO = %1.0f MeV",1000*off));

  /*
  if (bitmap&(1L<<0)) tdrDraw(paz,"Pz",kOpenCircle,kBlack);
  if (bitmap&(1L<<1)) tdrDraw(pbz,"Pz",kOpenCircle,kRed);
  if (bitmap&(1L<<2)) tdrDraw(pcz,"Pz",kOpenCircle,kGreen+2);
  if (bitmap&(1L<<3)) tdrDraw(pqz,"Pz",kOpenDiamond,kMagenta+2);
  if (bitmap&(1L<<4)) tdrDraw(pgz,"Pz",kOpenSquare,kBlue);
  */
  if (bitmap&(1L<<0)) tdrDraw(paz2,"Pz",kOpenCircle,kBlack);
  if (bitmap&(1L<<1)) tdrDraw(pbz2,"Pz",kOpenCircle,kRed);
  if (bitmap&(1L<<2)) tdrDraw(pcz2,"Pz",kOpenCircle,kGreen+2);
  if (bitmap&(1L<<3)) tdrDraw(pqz2,"Pz",kOpenDiamond,kMagenta+2);
  if (bitmap&(1L<<4)) tdrDraw(pgz2,"Pz",kOpenSquare,kBlue);

  //tdrDraw(pag,"Pz",kOpenCircle,kBlack);
  if (bitmap&(1L<<5)) tdrDraw(pag,"Pz",kFullCircle,kBlack);
  if (bitmap&(1L<<6)) tdrDraw(pbg,"Pz",kFullCircle,kRed);
  //tdrDraw(pcg,"Pz",kOpenCircle,kGreen+2);
  if (bitmap&(1L<<7)) tdrDraw(pcg,"Pz",kFullCircle,kGreen+2);
  if (bitmap&(1L<<8)) tdrDraw(pqg,"Pz",kFullDiamond,kMagenta+2);
  if (bitmap&(1L<<9)) tdrDraw(pgg,"Pz",kFullSquare,kBlue);

  TLegend *leg1 = tdrLeg(0.70,0.55,1.00,0.85); int ng(0);
  if (bitmap&(1L<<0)) { leg1->AddEntry(pag,"All","PL"); ++ng; }
  if (bitmap&(1L<<1)) { leg1->AddEntry(pbg,"b","PL"); ++ng; }
  if (bitmap&(1L<<2)) { leg1->AddEntry(pcg,"c","PL"); ++ng; }
  if (bitmap&(1L<<3)) { leg1->AddEntry(pqg,"q","PL"); ++ng; }
  if (bitmap&(1L<<4)) { leg1->AddEntry(pgg,"g","PL"); ++ng; }
  if (ng>0) leg1->SetHeader("vs genjet");
  leg1->SetY1NDC(0.80-0.05*ng); // why not working?

  /*
  TLegend *leg2 = tdrLeg(0.50,0.55,0.80,0.85); int nz(0);
  if (bitmap&(1L<<0)) { leg2->AddEntry(paz,"All","PL"); ++nz; }
  if (bitmap&(1L<<1)) { leg2->AddEntry(pbz,"b","PL"); ++nz; }
  if (bitmap&(1L<<2)) { leg2->AddEntry(pcz,"c","PL"); ++nz; }
  if (bitmap&(1L<<3)) { leg2->AddEntry(pqz,"q","PL"); ++nz; }
  if (bitmap&(1L<<4)) { leg2->AddEntry(pgz,"g","PL"); ++nz; }
  leg2->SetY1NDC(0.80-0.05*nz); // why not working?
  if (nz>0) leg2->SetHeader("vs gen-Z");
  */
  TLegend *leg2 = tdrLeg(0.50,0.55,0.80,0.85); int nz(0);
  if (bitmap&(1L<<0)) { leg2->AddEntry(paz2,"All","PL"); ++nz; }
  if (bitmap&(1L<<1)) { leg2->AddEntry(pbz2,"b","PL"); ++nz; }
  if (bitmap&(1L<<2)) { leg2->AddEntry(pcz2,"c","PL"); ++nz; }
  if (bitmap&(1L<<3)) { leg2->AddEntry(pqz2,"q","PL"); ++nz; }
  if (bitmap&(1L<<4)) { leg2->AddEntry(pgz2,"g","PL"); ++nz; }
  leg2->SetY1NDC(0.80-0.05*nz); // why not working?
  if (nz>0) leg2->SetHeader("from vs-Z");

  // Fit responses with constant offset bias
  if (bitmap&(1L<<0) || bitmap&(1L<<5)) {
    TF1 *fa = new TF1("fa","[0]+[1]*pow(x,[2])+[3]/x",ptmin,ptmax);
    fa->SetParameters(1.00,-0.1,-0.5,0.);
    pag->Fit(fa,"QRN");
    fa->SetLineColor(kBlack);
    fa->Draw("SAME");
  }
  if (bitmap&(1L<<1) || bitmap&(1L<<6)) {
    TF1 *fb = new TF1("fb","[0]+[1]*pow(x,[2])+[3]/x",ptmin,ptmax);
    fb->SetParameters(0.990,-0.1,-0.5,0.);
    pbg->Fit(fb,"QRN");
    fb->SetLineColor(kRed);
    fb->Draw("SAME");
  }
  if (bitmap&(1L<<2) || bitmap&(1L<<7)) {
    TF1 *fc = new TF1("fc","[0]+[1]*pow(x,[2])+[3]/x",ptmin,ptmax);
    fc->SetParameters(1.00,-0.1,-0.5,0.);
    pcg->Fit(fc,"QRN");
    fc->SetLineColor(kGreen+2);
    fc->Draw("SAME");
  }
  if (bitmap&(1L<<3) || bitmap&(1L<<8)) {
    TF1 *fq = new TF1("fq","[0]+[1]*pow(x,[2])+[3]/x",ptmin,ptmax);
    fq->SetParameters(1.01,-0.1,-0.5,0.);
    pqg->Fit(fq,"QRN");
    fq->SetLineColor(kMagenta+2);
    fq->Draw("SAME");
  }
  if (bitmap&(1L<<4) || bitmap&(1L<<9)) {
    TF1 *fg = new TF1("fg","[0]+[1]*pow(x,[2])+[3]/x",ptmin,ptmax);
    fg->SetParameters(0.985,-0.1,-0.5,0.);
    pgg->Fit(fg,"QRN");
    fg->SetLineColor(kBlue);
    fg->Draw("SAME");
  }

  c1->cd(2);
  gPad->SetLeftMargin(0.18);
  //gPad->SetRightMargin(0.02);
  gPad->SetLogx();

  l->DrawLine(ptmin,0,ptmax,0);

  TGraphErrors *gqz = new TGraphErrors(*pqz);
  TGraphErrors *graz = new TGraphErrors(*paz); subtractGraph(graz, gqz, 100,100);
  TGraphErrors *grbz = new TGraphErrors(*pbz); subtractGraph(grbz, gqz, 100,100);
  TGraphErrors *grcz = new TGraphErrors(*pcz); subtractGraph(grcz, gqz, 100,100);
  TGraphErrors *grgz = new TGraphErrors(*pgz); subtractGraph(grgz, gqz, 100,100);
  //
  TGraphErrors *grqz = new TGraphErrors(*pqz); subtractGraph(grqz, gqz, 100,100);

  // For mapped Z, need to subtract interpolated quark
  TGraphErrors *gqz2 = new TGraphErrors(*pqz2);
  TGraphErrors *graz2 = new TGraphErrors(*paz2); subtractGraphFit(graz2, gqz2, 100,100);
  TGraphErrors *grbz2 = new TGraphErrors(*pbz2); subtractGraphFit(grbz2, gqz2, 100,100);
  TGraphErrors *grcz2 = new TGraphErrors(*pcz2); subtractGraphFit(grcz2, gqz2, 100,100);
  TGraphErrors *grgz2 = new TGraphErrors(*pgz2); subtractGraphFit(grgz2, gqz2, 100,100);
  //
  TGraphErrors *grqz2 = new TGraphErrors(*pqz2); subtractGraphFit(grqz2, gqz2, 100,100);

  /*
  if (bitmap&(1L<<0)) tdrDraw(graz,"Pz",kOpenCircle,kBlack);
  if (bitmap&(1L<<1)) tdrDraw(grbz,"Pz",kOpenCircle,kRed);
  if (bitmap&(1L<<2)) tdrDraw(grcz,"Pz",kOpenCircle,kGreen+2);
  if (bitmap&(1L<<3)) tdrDraw(grqz,"Pz",kOpenDiamond,kMagenta+2);
  if (bitmap&(1L<<4)) tdrDraw(grgz,"Pz",kOpenSquare,kBlue);
  */
  if (bitmap&(1L<<0)) tdrDraw(graz2,"Pz",kOpenCircle,kBlack);
  if (bitmap&(1L<<1)) tdrDraw(grbz2,"Pz",kOpenCircle,kRed);
  if (bitmap&(1L<<2)) tdrDraw(grcz2,"Pz",kOpenCircle,kGreen+2);
  if (bitmap&(1L<<3)) tdrDraw(grqz2,"Pz",kOpenDiamond,kMagenta+2);
  if (bitmap&(1L<<4)) tdrDraw(grgz2,"Pz",kOpenSquare,kBlue);

  TGraphErrors *gqg = new TGraphErrors(*pqg);
  TGraphErrors *grag = new TGraphErrors(*pag); subtractGraph(grag, gqg, 100,100);
  TGraphErrors *grbg = new TGraphErrors(*pbg); subtractGraph(grbg, gqg, 100,100);
  TGraphErrors *grcg = new TGraphErrors(*pcg); subtractGraph(grcg, gqg, 100,100);
  TGraphErrors *grgg = new TGraphErrors(*pgg); subtractGraph(grgg, gqg, 100,100);
  //
  TGraphErrors *grqg = new TGraphErrors(*pqg); subtractGraph(grqg, gqg, 100,100);

  if (bitmap&(1L<<5)) tdrDraw(grag,"Pz",kFullCircle,kBlack);
  if (bitmap&(1L<<6)) tdrDraw(grbg,"Pz",kFullCircle,kRed);
  if (bitmap&(1L<<7)) tdrDraw(grcg,"Pz",kFullCircle,kGreen+2);
  if (bitmap&(1L<<8)) tdrDraw(grqg,"Pz",kFullDiamond,kMagenta+2);
  if (bitmap&(1L<<9)) tdrDraw(grgg,"Pz",kFullSquare,kBlue);

  c1->cd();
  c1->Update();

  c1->SaveAs(Form("pdf/drawZtruth_%d.pdf",(int)bitmap));
} // drawZtruths

// draw true Z+jet response, MPF and pTbal per flavor tag region
//void drawZtruthPerTag() {

//} // drawZtruthsPerTag();

// Demonstration of MPF vs pTbal biases vs truth information
// balance = genbalance
// MPF = genmet + unclust
// response = flavor
// UE?
// Later: update to per flavor
void drawZtruthMethod(string flavor="") {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  //TFile *f = new TFile("rootfiles/jme_bplusZ_merged_v12.root","READ");
  TFile *f = new TFile("rootfiles/jme_bplusZ_merged_v13.root","READ");
  assert(f && !f->IsZombie());  

  curdir->cd();

  // graph nomenclature
  // - 'jes' : pTreco/pTgen type of variables (true response)
  // - 'jet' : pTreco/pTZ type of variables (balancing)
  // - 'met' : MPF observables (full MET projections)
  // [- 'jme' : jet1 projections, MPF component]
  
  string sd = "mc/eta_00_13";
  const char *cf = flavor.c_str();
  if (flavor!="") {
    sd = Form("mc/eta_00_13_gen%s",cf);
  }
  const char *cd = sd.c_str();

  //const char *cd = "mc/eta_00_13";
  //const char *cd = "mc/eta_00_13_geng";
  //const char *cd = "mc/eta_00_13_genuds";
  // <pTreco/pTgen> vs pTz (truth / flavor)
  // "(p_{T}^{Z},ljet_pt/lgenjet_pt)"
  TGraphErrors *pgjes = (TGraphErrors*)f->Get(Form("%s/response_zmmjet_a30",cd));  // same as jet1
  assert(pgjes);
  // <pTreco/pTgen> vs pTgen (truth / flavor)
  // "(p_{T}^{Z},ljet_pt/lgenjet_pt)" => should be (pTgen,...)
  TGraphErrors *pgjes4 = (TGraphErrors*)f->Get(Form("%s/response_ptgen_zmmjet_a30",cd));
  assert(pgjes4);
  // <pTgen*cos(dphi)/pTZ> vs pTZ (genbalance, except for cos(dphi))
  // "(p_{T}^{Z},-leadingGenJet.Pt()*(cos(leadingGenJet.Phi())*Zboson.Px() + sin(leadingGenJet.Phi())*Zboson.Py())/(Zboson.Pt()*Zboson.Pt()))"x
  TGraphErrors *pgjet = (TGraphErrors*)f->Get(Form("%s/genjet1_zmmjet_a30",cd));
  assert(pgjet);
  // <pTreco/pTZ> vs pTZ (regular balance)
  // "(p_{T}^{Z},ljet_pt/Zpt)" => should this say lgenjet_pt/Zpt?
  TGraphErrors *pgjet2 = (TGraphErrors*)f->Get(Form("%s/bal_zmmjet_a30",cd));
  if (flavor!="") // TMP PATCH
    pgjet2 = (TGraphErrors*)f->Get(Form("%s/genjet1_zmmjet_a30",cd)); // TMP PATCH
  assert(pgjet2);
  // <genmet/pTZ> vs pTZ (genmet)
  // "(p_{T}^{Z},1 + GenMET_pt*(cos(GenMET_phi)*Zboson.Px() + sin(GenMET_phi)*Zboson.Py())/(Zboson.Pt()*Zboson.Pt()))"
  TGraphErrors *pgmet = (TGraphErrors*)f->Get(Form("%s/gmpf_zmmjet_a30",cd));
  assert(pgmet);
  //
  // <pTreco/pTZ> vs pTZ (pT balance)
  // "(p_{T}^{Z},RpT)" => what's Rpt? Not pTreco/pTZ (bal is a bit lower) or pTreco/pTgen (response & jet1 are much higher)?
  // => likely pTreco/pTZ, just that bal was mislabeled
  TGraphErrors *pjet = (TGraphErrors*)f->Get(Form("%s/ptchs_zmmjet_a30",cd));
  if (flavor!="") // TMP PATCH
    pjet = (TGraphErrors*)f->Get(Form("%s/mpfjet1_zmmjet_a30",cd)); // TMP PATCH
  assert(pjet);
  // <met/pTZ> vs pTZ (MPF)
  // "(p_{T}^{Z},RMPF)"
  TGraphErrors *pmet = (TGraphErrors*)f->Get(Form("%s/mpfchs_zmmjet_a30",cd));
  assert(pmet);
  // <met1/pTZ> vs pTZ (MPF)
  // "(p_{T}^{Z},RMPFjet1)" => compare to genjet1
  TGraphErrors *pmet1 = (TGraphErrors*)f->Get(Form("%s/mpfjet1_zmmjet_a30",cd));
  assert(pmet1);

  TH1D *h = new TH1D("h",";p_{T,Z} (GeV);Observable",1470,30,1500);
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  h->SetMinimum(0.83);
  h->SetMaximum(1.13);
  if (flavor!="") {
    h->SetYTitle(Form("Observable (%s)",cf));
    h->SetMinimum(0.75);
    h->SetMaximum(1.19);
  }

  lumi_13TeV = "UL2017 Z+jet MC (v13)";
  TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);
  c1->SetLogx();

  tdrDraw(pgjes,"Pz",kOpenSquare,kGreen+1);
  tdrDraw(pgjes4,"Pz",kFullSquare,kGreen+4); pgjes4->SetMarkerSize(0.6);
  //tdrDraw(pgjet,"Pz",kOpenCircle,kRed+2); pgjet->SetMarkerSize(0.8);
  tdrDraw(pgjet2,"Pz",kOpenCircle,kRed); //pgjet2->SetMarkerSize(0.8);
  tdrDraw(pgmet,"Pz",kOpenDiamond,kBlue);
  //
  tdrDraw(pjet,"Pz",kFullCircle,kRed);
  tdrDraw(pmet,"Pz",kFullDiamond,kBlue);
  //
  //tdrDraw(pmet1,"Pz",kFullCircle,kRed+2); pmet1->SetMarkerSize(0.8);

  TGraphErrors *pgjes2 = (TGraphErrors*)pjet->Clone("pgjes2");
  assert(pgjes2->GetN()==pgjet->GetN());
  for (int i = 0; i != pgjes2->GetN(); ++i) {
    pgjes2->SetPoint(i, pjet->GetX()[i], pjet->GetY()[i]/pgjet->GetY()[i]);
  } // for i

  TGraphErrors *pgjes3 = (TGraphErrors*)pjet->Clone("pgjes3");
  assert(pgjes3->GetN()==pgjet2->GetN());
  for (int i = 0; i != pgjes3->GetN(); ++i) {
    pgjes3->SetPoint(i, pjet->GetX()[i], pjet->GetY()[i]/pgjet2->GetY()[i]);
  } // for i

  TGraphErrors *pgjes5 = (TGraphErrors*)pjet->Clone("pgjes5");
  assert(pgjes5->GetN()==pgjet2->GetN());
  for (int i = 0; i != pgjes5->GetN(); ++i) {
    pgjes5->SetPoint(i, pjet->GetX()[i] * pgjet2->GetY()[i],
		     pjet->GetY()[i]/pgjet2->GetY()[i]);
  } // for i

  //tdrDraw(pgjes2,"Pz",kFullSquare,kGreen+2); pgjes2->SetMarkerSize(0.8);
  tdrDraw(pgjes3,"Pz",kOpenSquare,kGreen+2); pgjes3->SetMarkerSize(0.8);
  tdrDraw(pgjes5,"Pz",kOpenSquare,kGreen+3); pgjes5->SetMarkerSize(0.7);
  
  TF1 *f4 = new TF1("f4","1+[0]/x",30,1500.);
  pgjes4->Fit(f4,"QRN");
  f4->SetLineColor(kGreen+3);
  f4->Draw("SAME");

  TF1 *f5 = new TF1("f5","1+[0]/x",30,1500.);
  pgjes5->Fit(f5,"QRN");
  f5->SetLineColor(kGreen+4);
  f5->Draw("SAME");

  TLegend *legjes = tdrLeg(0.36,0.67,0.66,0.83);
  legjes->SetTextSize(0.030); legjes->SetTextColor(kGreen+2);
  //legjes->AddEntry(pgjes2,"#LTp_{T,rec}#GT / #LTp_{T,gen}cos(#Delta#phi)#GT (ptchs / genjet1)","PL");
  legjes->AddEntry(pgjes,"#LTp_{T,rec} / p_{T,gen}#GT (response vs p_{T,Z})","PL");
  legjes->AddEntry(pgjes3,"#LTp_{T,rec}#GT / #LTp_{T,gen}#GT (ptchs / bal vs p_{T,Z})","PL");
  legjes->AddEntry(pgjes5,"#LTp_{T,rec}#GT / #LTp_{T,gen}#GT (ptchs / bal vs #LTp_{T,gen}#GT)","PL");
  legjes->AddEntry(pgjes4,"#LTp_{T,rec} / p_{T,gen}#GT (response vs p_{T,gen})","PL");

  TLegend *legmet = tdrLeg(0.42,0.415,0.72,0.495);
  legmet->SetTextSize(0.030); legmet->SetTextColor(kBlue);
  legmet->AddEntry(pgmet,"#LTp_{T,genMET}#GT / #LTp_{T,Z}#GT (gmpf)","PL");
  legmet->AddEntry(pmet,"#LTp_{T,recMET}#GT / #LTp_{T,Z}#GT (mpfchs)","PL");

  //TLegend *legjet = tdrLeg(0.48,0.16,0.78,0.32);
  TLegend *legjet = tdrLeg(0.48,0.24,0.78,0.32);
  legjet->SetTextSize(0.030); legjet->SetTextColor(kRed);
  legjet->AddEntry(pjet,"#LTp_{T,rec}#GT / #LTp_{T,Z}#GT (ptchs)","PL");
  legjet->AddEntry(pgjet2,"#LTp_{T,gen}#GT / #LTp_{T,Z}#GT (bal)","PL");
  //legjet->AddEntry(pmet1,"#LTp_{T,rec}cos(#Delta#phi)#GT / #LTp_{T,Z}#GT (mpfjet1)","PL");
  //legjet->AddEntry(pgjet,"#LTp_{T,gen}cos(#Delta#phi)#GT / #LTp_{T,Z}#GT (genjet1)","PL");


  if (flavor=="") c1->SaveAs("pdf/drawZtruth_method.pdf");
  else            c1->SaveAs(Form("pdf/drawZtruth_method_%s.pdf",cf));
} // drawZtruthMethod

void drawZtruth() {
  
  drawZtruths(1023);
  drawZtruths(792);
  drawZtruths(825);
  drawZtruthMethod();
  drawZtruthMethod("uds");
  drawZtruthMethod("g");
  drawZtruthMethod("b");
}
