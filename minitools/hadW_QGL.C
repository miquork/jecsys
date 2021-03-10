// Purpose: Derive quark and b-jet (also c-jet?) scale factors for QGL templates
//          These are needed for Zflavor.C

#include "TFile.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TLine.h"
#include "TF1.h"
#include "TMultiGraph.h"

#include "../tdrstyle_mod15.C"

#include <map>
#include <string>

// Plot MC on top pad?
bool plotMChist = false; // c1 QGL template
bool plotMC = false;//true; // c0 efficiency vs pT
bool debug = false;

void getDijetQGL(TFile *fjd, TFile *fjm, TFile *fjd2, TFile *fjm2,
		 TH2D **h2jd, TH2D **h2jm, TH2D **h2jmq, TH2D **h2jmg);

void hadW_QGL() {
  
  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TLine *l = new TLine();
  l->SetLineStyle(kDotted);

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);

  //TFile *fd = new TFile("rootfiles/hadWUL1718V5_EMUF.root","READ");
  TFile *fd = new TFile("rootfiles/hadWUL1718V5_Glu_v4.root","READ");
  //TFile *fd = new TFile("rootfiles/hadWUL1718V5_Glu.root","READ");
  assert(fd && !fd->IsZombie());
  //TFile *fm = new TFile("rootfiles/hadWMC1718V5_EMUF.root","READ");
  TFile *fm = new TFile("rootfiles/hadWMC1718V5_Glu_v4.root","READ");
  //TFile *fm = new TFile("rootfiles/hadWMC1718V5_Glu.root","READ");
  assert(fm && !fm->IsZombie());
  TFile *fz = new TFile("rootfiles/jecdata2018ABCD.root","READ");
  assert(fz && !fz->IsZombie());
  //TFile *fz2 = new TFile("rootfiles/jme_bplusZ_merged_v33.root","READ");
  TFile *fz2 = new TFile("rootfiles/jme_bplusZ_merged_v34.root","READ");
  assert(fz2 && !fz2->IsZombie());

  // Later: combine 0.5 bins up to |eta|<2.5, merge jet triggers
  // Need to also consider bias from JEC of QGL>0.5 vs QGL<0.5 with steep pT
  // (Z+jet suggets at least tagged quark responses differ)
  TFile *fjd = new TFile("rootfiles/output-DATA-2a-UL18V2V3_ABCD.root","READ");
  assert(fjd && !fjd->IsZombie());
  TFile *fjd2 = new TFile("rootfiles/output-DATA-2a-UL17V4_BCDEF.root","READ");
  assert(fjd2 && !fjd2->IsZombie());
  TFile *fjm = new TFile("rootfiles/output-MC-1-UL18V2V3_ABCD.root","READ");
  assert(fjm && !fjm->IsZombie());
  TFile *fjm2 = new TFile("rootfiles/output-MC-1-UL17V4_BCDEF.root","READ");
  assert(fjm2 && !fjm2->IsZombie());


  // More accurate pT binning per jet, but spiky
  TH2D *h2qd = (TH2D*)fd->Get("h2qglqb"); assert(h2qd);
  TH2D *h2qm = (TH2D*)fm->Get("h2qglqb"); assert(h2qm);
  TH2D *h2bd = (TH2D*)fd->Get("h2qglbb"); assert(h2bd);
  TH2D *h2bm = (TH2D*)fm->Get("h2qglbb"); assert(h2bm);
  /*
  // Smoother QGL distributions, but also more smeared pT range
  // Could be lower bias from JER, though?
  TH2D *h2qd = (TH2D*)fd->Get("h2qglqa"); assert(h2qd);
  TH2D *h2qm = (TH2D*)fm->Get("h2qglqa"); assert(h2qm);
  TH2D *h2bd = (TH2D*)fd->Get("h2qglba"); assert(h2bd);
  TH2D *h2bm = (TH2D*)fm->Get("h2qglba"); assert(h2bm);
  */
  TH2D *h2gd = (TH2D*)fd->Get("h2qglg"); assert(h2gd);
  TH2D *h2gm = (TH2D*)fm->Get("h2qglg"); assert(h2gm);

  // Split by true jet flavor
  TH2D *h2qum = (TH2D*)fm->Get("h2qglqb_u"); assert(h2qum);
  TH2D *h2qdm = (TH2D*)fm->Get("h2qglqb_d"); assert(h2qdm);
  TH2D *h2qsm = (TH2D*)fm->Get("h2qglqb_s"); assert(h2qsm);
  TH2D *h2qcm = (TH2D*)fm->Get("h2qglqb_c"); assert(h2qcm);
  TH2D *h2qbm = (TH2D*)fm->Get("h2qglqb_b"); assert(h2qbm);
  TH2D *h2qgm = (TH2D*)fm->Get("h2qglqb_g"); assert(h2qgm);

  TH2D *h2ggm = (TH2D*)fm->Get("h2qglg_g"); assert(h2ggm);
  TH2D *h2gum = (TH2D*)fm->Get("h2qglg_u"); assert(h2gum);
  TH2D *h2gdm = (TH2D*)fm->Get("h2qglg_d"); assert(h2gdm);
  TH2D *h2gsm = (TH2D*)fm->Get("h2qglg_s"); assert(h2gsm);
  TH2D *h2gcm = (TH2D*)fm->Get("h2qglg_c"); assert(h2gcm);
  TH2D *h2gbm = (TH2D*)fm->Get("h2qglg_b"); assert(h2gbm);

  // Inclusive QGL distributions, no pT binning involved; smoothest of all
  TH1D *h1q0d = (TH1D*)fd->Get("hqglq"); assert(h1q0d);
  TH1D *h1q0m = (TH1D*)fm->Get("hqglq"); assert(h1q0m);
  TH1D *h1b0d = (TH1D*)fd->Get("hqglb"); assert(h1b0d);
  TH1D *h1b0m = (TH1D*)fm->Get("hqglb"); assert(h1b0m);
  TH1D *h1g0d = (TH1D*)fd->Get("hqglg"); assert(h1g0d);
  TH1D *h1g0m = (TH1D*)fm->Get("hqglg"); assert(h1g0m);

  // Z+jet
  TH2D *h2zd = (TH2D*)fz2->Get("data/eta_00_25/h_JetPt_QGL_alpha100");
  TH2D *h2zm = (TH2D*)fz2->Get("mc/eta_00_25/h_JetPt_QGL_alpha100");
  assert(h2zd);
  assert(h2zm);
  // Dijet
  TH2D *h2jd(0), *h2jm(0), *h2jmq(0), *h2jmg(0);
  getDijetQGL(fjd, fjm, fjd2, fjm2, &h2jd, &h2jm, &h2jmq, &h2jmg);
  //TH2D *h2jd = (TH2D*)fjd->Get("Standard/Eta_0.0-1.3/jt0/hqgl2");
  //TH2D *h2jm = (TH2D*)fjm->Get("Standard/Eta_0.0-1.3/jt0/hqgl2");
  assert(h2jd);
  assert(h2jm);
  assert(h2jmq);
  assert(h2jmg);

  TH1D *hnj  = h2jm->ProjectionX("hnj");
  TH1D *hfjq = h2jmq->ProjectionX("hfjq");
  hfjq->Divide(hnj);
  TH1D *hfjg = h2jmg->ProjectionX("hfjg");
  hfjg->Divide(hnj);

  curdir->cd();

  TH1D *hqed = h2qd->ProjectionY("hqed",1,1); hqed->Reset();
  TH1D *hqem = h2qd->ProjectionY("hqem",1,1); hqem->Reset();
  TH1D *hqe  = h2qd->ProjectionY("hqe",1,1);  hqe->Reset();
  TH1D *hbed = h2bd->ProjectionY("hbed",1,1); hbed->Reset();
  TH1D *hbem = h2bd->ProjectionY("hbem",1,1); hbem->Reset();
  TH1D *hbe  = h2bd->ProjectionY("hbe",1,1);  hbe->Reset();
  TH1D *hged = h2gd->ProjectionY("hged",1,1); hged->Reset();
  TH1D *hgem = h2gd->ProjectionY("hgem",1,1); hgem->Reset();
  TH1D *hge  = h2gd->ProjectionY("hge",1,1);  hge->Reset();

  TH1D *hzed = h2zd->ProjectionX("hzed",1,1); hzed->Reset();
  TH1D *hzem = h2zd->ProjectionX("hzem",1,1); hzem->Reset();
  TH1D *hze  = h2zd->ProjectionX("hze",1,1);  hze->Reset();

  TH1D *hjed = h2jd->ProjectionX("hjed",1,1); hjed->Reset();
  TH1D *hjem = h2jd->ProjectionX("hjem",1,1); hjem->Reset();
  TH1D *hjemq = h2jd->ProjectionX("hjemq",1,1); hjemq->Reset();
  TH1D *hjemg = h2jd->ProjectionX("hjemg",1,1); hjemg->Reset();
  TH1D *hje  = h2jd->ProjectionX("hje",1,1);  hje->Reset();
  TH1D *hjeq  = h2jd->ProjectionX("hjeq",1,1);  hjeq->Reset();
  TH1D *hjeg  = h2jd->ProjectionX("hjeg",1,1);  hjeg->Reset();


  TH1D *hziim = (TH1D*)fz->Get("mc/eta00-25/counts_zii_a100");
  assert(hziim);

  TH1D *hzium = (TH1D*)fz->Get("mc/eta00-25/counts_ziq_a100");
  assert(hzium);
  TH1D *hzuum = (TH1D*)fz->Get("mc/eta00-25/counts_zqq_a100");
  assert(hzuum);
  TH1D *hzicm = (TH1D*)fz->Get("mc/eta00-25/counts_zic_a100");
  assert(hzicm);
  TH1D *hzucm = (TH1D*)fz->Get("mc/eta00-25/counts_zqc_a100");
  assert(hzucm);
  TH1D *hzccm = (TH1D*)fz->Get("mc/eta00-25/counts_zcc_a100");
  assert(hzccm);
  // Add charm to uds to stay compatible with W>qq' and dijet
  TH1D *hzqqm = (TH1D*)hzuum->Clone("hzqqm");
  hzqqm->Add(hzucm);
  hzqqm->Add(hzccm,0.5);//approximation
  TH1D *hziqm = (TH1D*)hzium->Clone("hziqm");
  hziqm->Add(hzicm);
  TH1D *hzqm = (TH1D*)hzqqm->Clone("hzqm");
  hzqm->Divide(hziqm);

  TH1D *hzigm = (TH1D*)fz->Get("mc/eta00-25/counts_zig_a100");
  assert(hzigm);
  TH1D *hzggm = (TH1D*)fz->Get("mc/eta00-25/counts_zgg_a100");
  assert(hzggm);
  TH1D *hzgm = (TH1D*)hzggm->Clone("hzgm");
  hzgm->Divide(hzigm);
  for (int i = 1; i != hzgm->GetNbinsX()+1; ++i) 
    hzgm->SetBinContent(i, 1-hzgm->GetBinContent(i));

  TH1D *hfzq = (TH1D*)hziqm->Clone("hfzq");
  hfzq->Add(hzicm);
  hfzq->Divide(hziim);
  TH1D *hfzg = (TH1D*)hzigm->Clone("hfzq");
  hfzg->Divide(hziim);

  //for (int j = 0; j != h2qd->GetNbinsY()+1; ++j) {
  for (int j = 1; j != h2qd->GetNbinsY()+1; ++j) {

    if (debug) cout << j<<","<<endl<<flush;

    TH1D *hu = new TH1D(Form("hu_%d",j),";QGL;Fraction",100,0,1);
    hu->SetMinimum(5e-3);//0);
    hu->SetMaximum(0.4);//0.26);
    TH1D *hd = new TH1D(Form("hd_%d",j),";QGL;Data/MC",100,0,1);
    hd->SetMinimum(0.0);//0.20);//0.25);//0.35);
    hd->SetMaximum(2.0);//1.80);//1.75);//1.65);

    TH1D *hu2 = new TH1D(Form("hu2_%d",j),";QGL;Fraction",100,0,1);
    hu2->SetMinimum(5e-3);
    hu2->SetMaximum(0.4);
    TH1D *hd2 = new TH1D(Form("hd2_%d",j),";QGL;Data/MC",100,0,1);
    hd2->SetMinimum(0.0);//0.4);//0.6);//0.20);
    hd2->SetMaximum(2.0);//1.6);//1.4);//1.80);

    //lumi_13TeV = "[TT l+jet, W>qq'] UL17+18, 101.4 fb^{-1}";
    //lumi_13TeV = "[TT l+jet] UL17+18, 101.4 fb^{-1}";
    lumi_13TeV = "UL17+18, 101.4 fb^{-1}";
    TCanvas *c1 = tdrDiCanvas(Form("c1_%d",j),hu,hd,4,11);
    TCanvas *c2 = tdrDiCanvas(Form("c2_%d",j),hu2,hd2,4,11);
    
    c1->cd(1);
    gPad->SetLogy();
    
    tex->DrawLatex(0.19,0.70,Form("|#eta| < 2.5, %1.0f--%1.0f GeV",
				  h2qd->GetYaxis()->GetBinLowEdge(j),
				  h2qd->GetYaxis()->GetBinLowEdge(j+1)));
    if (!plotMChist) tex->DrawLatex(0.19,0.65,"Data only");

    TLegend *leg1 = tdrLeg(0.65,0.90-5*0.05,0.85,0.90);

    c2->cd(1);
    gPad->SetLogy();

    c1->cd(2);
    l->DrawLine(0,1,1,1);
  
  //for (int j =0; j != h2qd->GetNbinsY()+1; ++j) {

    double pt = h2qd->GetYaxis()->GetBinCenter(j);
    double ptmin = h2qd->GetYaxis()->GetBinLowEdge(j);
    double ptmax = h2qd->GetYaxis()->GetBinLowEdge(j+1);
    //if (!(pt>80 && pt<110)) continue;
    //if (!(pt>45 && pt<50)) continue;
    //if (!(pt>110 && pt<150)) continue;
    //if (!(pt>150 && pt<200)) continue;
    if (j!=0 && pt<28) continue;

    TH1D *h1qd = h2qd->ProjectionX(Form("h1qd_%d",j),j,j);
    TH1D *h1qm = h2qm->ProjectionX(Form("h1qm_%d",j),j,j);
    TH1D *h1bd = h2bd->ProjectionX(Form("h1bd_%d",j),j,j);
    TH1D *h1bm = h2bm->ProjectionX(Form("h1bm_%d",j),j,j);
    TH1D *h1gd = h2gd->ProjectionX(Form("h1gd_%d",j),j,j);
    TH1D *h1gm = h2gm->ProjectionX(Form("h1gm_%d",j),j,j);

    TH1D *h1qum = h2qum->ProjectionX(Form("h1qum_%d",j),j,j);
    TH1D *h1qdm = h2qdm->ProjectionX(Form("h1qdm_%d",j),j,j);
    TH1D *h1qsm = h2qsm->ProjectionX(Form("h1qsm_%d",j),j,j);
    TH1D *h1qcm = h2qcm->ProjectionX(Form("h1qcm_%d",j),j,j);
    TH1D *h1qbm = h2qbm->ProjectionX(Form("h1qbm_%d",j),j,j);
    TH1D *h1qgm = h2qgm->ProjectionX(Form("h1qgm_%d",j),j,j);

    TH1D *h1ggm = h2ggm->ProjectionX(Form("h1ggm_%d",j),j,j);
    TH1D *h1gum = h2gum->ProjectionX(Form("h1gum_%d",j),j,j);
    TH1D *h1gdm = h2gdm->ProjectionX(Form("h1gdm_%d",j),j,j);
    TH1D *h1gsm = h2gsm->ProjectionX(Form("h1gsm_%d",j),j,j);
    TH1D *h1gcm = h2gcm->ProjectionX(Form("h1gcm_%d",j),j,j);
    TH1D *h1gbm = h2gbm->ProjectionX(Form("h1gbm_%d",j),j,j);

    int kmin = h2zd->GetXaxis()->FindBin(ptmin+1e-3);
    int kmax = h2zd->GetXaxis()->FindBin(ptmax-1e-3);
    TH1D *h1zd = h2zd->ProjectionY(Form("h1zd_%d",j),kmin,kmax);
    TH1D *h1zm = h2zm->ProjectionY(Form("h1zm_%d",j),kmin,kmax);
    int lmin = h2jd->GetXaxis()->FindBin(ptmin+1e-3);
    int lmax = h2jd->GetXaxis()->FindBin(ptmax-1e-3);
    TH1D *h1jd = h2jd->ProjectionY(Form("h1jd_%d",j),lmin,lmax);
    TH1D *h1jm = h2jm->ProjectionY(Form("h1jm_%d",j),lmin,lmax);
    TH1D *h1jmq = h2jmq->ProjectionY(Form("h1jmq_%d",j),lmin,lmax);
    TH1D *h1jmg = h2jmg->ProjectionY(Form("h1jmg_%d",j),lmin,lmax);

    if (j==0) { // inclusive case
      h1qd = h1q0d;
      h1qm = h1q0m;
      h1bd = h1b0d;
      h1bm = h1b0m;
      h1gd = h1g0d;
      h1gm = h1g0m;
    }

    int nbin = 2;//4;
    h1qd->Rebin(nbin);
    h1qm->Rebin(nbin);
    h1bd->Rebin(nbin);
    h1bm->Rebin(nbin);
    h1gd->Rebin(nbin);
    h1gm->Rebin(nbin);

    h1qum->Rebin(nbin);
    h1qdm->Rebin(nbin);
    h1qsm->Rebin(nbin);
    h1qcm->Rebin(nbin);
    h1qbm->Rebin(nbin);
    h1qgm->Rebin(nbin);

    h1ggm->Rebin(nbin);
    h1gum->Rebin(nbin);
    h1gdm->Rebin(nbin);
    h1gsm->Rebin(nbin);
    h1gcm->Rebin(nbin);
    h1gbm->Rebin(nbin);

    if (debug) cout << "Rebin Z+jet and dijet" << endl << flush;

    h1zd->Rebin(nbin);
    h1zm->Rebin(nbin);
    h1jd->Rebin(nbin);
    h1jm->Rebin(nbin);
    h1jmq->Rebin(nbin);
    h1jmg->Rebin(nbin);

    h1qd->Scale(1./h1qd->Integral());
    h1qm->Scale(1./h1qm->Integral());
    h1bd->Scale(1./h1bd->Integral());
    h1bm->Scale(1./h1bm->Integral());
    h1gd->Scale(1./h1gd->Integral());
    h1gm->Scale(1./h1gm->Integral());

    h1qum->Scale(1./h1qum->Integral());
    h1qdm->Scale(1./h1qdm->Integral());
    h1qsm->Scale(1./h1qsm->Integral());
    h1qcm->Scale(1./h1qcm->Integral());
    h1qbm->Scale(1./h1qbm->Integral());
    h1qgm->Scale(1./h1qgm->Integral());

    h1ggm->Scale(1./h1ggm->Integral());
    h1gum->Scale(1./h1gum->Integral());
    h1gdm->Scale(1./h1gdm->Integral());
    h1gsm->Scale(1./h1gsm->Integral());
    h1gcm->Scale(1./h1gcm->Integral());
    h1gbm->Scale(1./h1gbm->Integral());

    if (debug) cout << "Scale Z+jet and dijet" << endl << flush;

    h1zd->Scale(1./h1zd->Integral());
    h1zm->Scale(1./h1zm->Integral());
    h1jd->Scale(1./h1jd->Integral());
    h1jm->Scale(1./h1jm->Integral());
    h1jmq->Scale(1./h1jmq->Integral());
    h1jmg->Scale(1./h1jmg->Integral());
    
    // Calculate efficiencies
    int k0 = h1qd->FindBin(0.);
    int k1 = h1qd->FindBin(0.5+1e-3);
    int k2 = h1qd->GetNbinsX()+1;
    double qede(0), qeme(0);
    double qed = h1qd->IntegralAndError(k1,k2,qede) /  h1qd->Integral(k0,k2);
    double qem = h1qm->IntegralAndError(k1,k2,qeme) /  h1qm->Integral(k0,k2);
    hqed->SetBinContent(j, qed);
    hqem->SetBinContent(j, qem);
    hqe ->SetBinContent(j, qed / qem);
    hqed->SetBinError(j, qede);
    hqem->SetBinError(j, qeme);
    hqe ->SetBinError(j, qed / qem * sqrt(pow(qede/qed,2)+pow(qeme/qem,2)));
    /*
    double bed =
      h1bd->Integral(h1bd->FindBin(0.5+1e-3),h1bd->GetNbinsX()) /
      h1bd->Integral();
    double bem =
      h1bm->Integral(h1bm->FindBin(0.5+1e-3),h1bm->GetNbinsX()) /
      h1bm->Integral();
    hbed->SetBinContent(j, bed);
    hbem->SetBinContent(j, bem);
    hbe ->SetBinContent(j, bed / bem);
    */
    double bede(0), beme(0);
    double bed = h1bd->IntegralAndError(k1,k2,bede) /  h1bd->Integral(k0,k2);
    double bem = h1bm->IntegralAndError(k1,k2,beme) /  h1bm->Integral(k0,k2);
    hbed->SetBinContent(j, bed);
    hbem->SetBinContent(j, bem);
    hbe ->SetBinContent(j, bed / bem);
    hbed->SetBinError(j, bede);
    hbem->SetBinError(j, beme);
    hbe ->SetBinError(j, bed / bem * sqrt(pow(bede/bed,2)+pow(beme/bem,2)));

    double gede(0), geme(0);
    double ged = h1gd->IntegralAndError(k1,k2,gede) /  h1gd->Integral(k0,k2);
    double gem = h1gm->IntegralAndError(k1,k2,geme) /  h1gm->Integral(k0,k2);
    hged->SetBinContent(j, ged);
    hgem->SetBinContent(j, gem);
    hge ->SetBinContent(j, ged / gem);
    hged->SetBinError(j, gede);
    hgem->SetBinError(j, geme);
    hge ->SetBinError(j, ged / gem * sqrt(pow(gede/ged,2)+pow(geme/gem,2)));

    int m0 = h1zd->FindBin(0.);
    int m1 = h1zd->FindBin(0.5+1e-3);
    int m2 = h1zd->GetNbinsX()+1;
    double zede(0), zeme(0);
    double zed = h1zd->IntegralAndError(m1,m2,zede) /  h1zd->Integral(m0,m2);
    double zem = h1zm->IntegralAndError(m1,m2,zeme) /  h1zm->Integral(m0,m2);
    hzed->SetBinContent(kmin, zed);
    hzem->SetBinContent(kmin, zem);
    hze ->SetBinContent(kmin, zed / zem);
    hzed->SetBinError(kmin, zede);
    hzem->SetBinError(kmin, zeme);
    hze ->SetBinError(kmin, zed / zem * sqrt(pow(zede/zed,2)+pow(zeme/zem,2)));

    int l0 = h1jd->FindBin(0.);
    int l1 = h1jd->FindBin(0.5+1e-3);
    int l2 = h1jd->GetNbinsX()+1;
    double jede(0), jeme(0), jemqe(0), jemge(0);
    double jed = h1jd->IntegralAndError(l1,l2,jede) /  h1jd->Integral(l0,l2);
    double jem = h1jm->IntegralAndError(l1,l2,jeme) /  h1jm->Integral(l0,l2);
    double jemq = h1jmq->IntegralAndError(l1,l2,jemqe) / h1jmq->Integral(l0,l2);
    double jemg = h1jmg->IntegralAndError(l1,l2,jemge) / h1jmg->Integral(l0,l2);
    hjed->SetBinContent(lmin, jed);
    hjem->SetBinContent(lmin, jem);
    hjemq->SetBinContent(lmin, jemq);
    hjemg->SetBinContent(lmin, jemg);
    hje ->SetBinContent(lmin, jed / jem);
    hjeq ->SetBinContent(lmin, jemq / jem);
    hjeg ->SetBinContent(lmin, jemg / jem);
    hjed->SetBinError(lmin, jede);
    hjem->SetBinError(lmin, jeme);
    hjemq->SetBinError(lmin, jemqe);
    hjemg->SetBinError(lmin, jemge);
    hje ->SetBinError(lmin, jed / jem * sqrt(pow(jede/jed,2)+pow(jeme/jem,2)));

    c1->cd(1);
    
    if (debug) cout << "Draw c1 pad1" << endl << flush;

    if (plotMChist) {
      tdrDraw(h1qm,"HISTE",kNone,kBlue,kSolid,-1,1001,kBlue-9);
      tdrDraw(h1bm,"HISTE",kNone,kRed,kSolid,-1,1001,kRed-9);
      h1bm->SetFillColorAlpha(kRed-9,0.5);
      tdrDraw(h1gm,"HISTE",kNone,kOrange+2,kSolid,-1,1001,kOrange-9);
      h1gm->SetFillColorAlpha(kOrange-9,0.5);

      tdrDraw(h1zm,"HISTE",kNone,kMagenta+2,kSolid,-1,1001,kMagenta-9);
      h1zm->SetFillColorAlpha(kMagenta-9,0.5);
      tdrDraw(h1jm,"HISTE",kNone,kGreen+2,kSolid,-1,1001,kGreen-9);
      h1jm->SetFillColorAlpha(kGreen-9,0.5);

      tdrDraw(h1jmq,"HISTE",kNone,kGreen+2,kSolid,-1,3001,kGreen-9);
      h1jmq->SetFillColorAlpha(kGreen-9,0.5);
      tdrDraw(h1jmg,"HISTE",kNone,kGreen+2,kSolid,-1,1003,kGreen-9);
      h1jmg->SetFillColorAlpha(kGreen-9,0.5);
    }

    tdrDraw(h1bd,"Pz",kFullSquare,kRed);
    tdrDraw(h1qd,"Pz",kFullCircle,kBlue);
    tdrDraw(h1gd,"Pz",kFullDiamond,kOrange+2);

    tdrDraw(h1zd,"Pz",kOpenCircle,kMagenta+2);
    tdrDraw(h1jd,"Pz",kOpenDiamond,kGreen+2);

    leg1->AddEntry(h1qd,"W>qq' (q rich)","PLE");
    leg1->AddEntry(h1zd,"Z+jet (q rich)","PLE");
    leg1->AddEntry(h1bd,"t>Wb (b rich)","PLE");
    leg1->AddEntry(h1jd,"Dijet (g rich)","PLE");
    leg1->AddEntry(h1gd,"TT+jet (g rich)","PLE");
    
    gPad->RedrawAxis();

    c2->cd(1);

    if (debug) cout << "Draw c2 pad1" << endl << flush;

    tdrDraw(h1qm,"HISTE",kNone,kBlue,kSolid,-1,1001,kBlue-9);
    tdrDraw(h1gm,"HISTE",kNone,kOrange+2,kSolid,-1,1001,kOrange-9);
    h1gm->SetFillColorAlpha(kRed-9,0.5);

    if (false) {
      tdrDraw(h1qum,"Pz",kOpenCircle,kBlue);
      tdrDraw(h1qdm,"Pz",kOpenCircle,kCyan+2);
      tdrDraw(h1qsm,"Pz",kOpenCircle,kMagenta+2);
      tdrDraw(h1qcm,"Pz",kOpenCircle,kGreen+2);
      tdrDraw(h1qbm,"Pz",kOpenCircle,kRed);
      tdrDraw(h1qgm,"Pz",kOpenCircle,kOrange+2);
    }

    tdrDraw(h1ggm,"Pz",kOpenDiamond,kOrange+2);
    tdrDraw(h1gum,"Pz",kOpenDiamond,kBlue);
    tdrDraw(h1gdm,"Pz",kOpenDiamond,kCyan+2);
    tdrDraw(h1gsm,"Pz",kOpenDiamond,kMagenta+2);
    tdrDraw(h1gcm,"Pz",kOpenDiamond,kGreen+2);
    tdrDraw(h1gbm,"Pz",kOpenDiamond,kRed);

    gPad->RedrawAxis();

    c1->cd(2);

    if (debug) cout << "Draw c1 pad2" << endl << flush;

    TH1D *h1q = (TH1D*)h1qd->Clone(Form("h1q_%d",j));
    h1q->Divide(h1qm);
    TH1D *h1b = (TH1D*)h1bd->Clone(Form("h1b_%d",j));
    h1b->Divide(h1bm);
    TH1D *h1g = (TH1D*)h1gd->Clone(Form("h1g_%d",j));
    h1g->Divide(h1gm);

    TH1D *h1z = (TH1D*)h1zd->Clone(Form("h1z_%d",j));
    h1z->Divide(h1zm);
    TH1D *h1j = (TH1D*)h1jd->Clone(Form("h1j_%d",j));
    h1j->Divide(h1jm);

    tdrDraw(h1b,"Pz",kFullSquare,kRed);
    tdrDraw(h1q,"Pz",kFullCircle,kBlue);
    tdrDraw(h1g,"Pz",kFullDiamond,kOrange+2);

    tdrDraw(h1z,"Pz",kOpenCircle,kMagenta+2);
    tdrDraw(h1j,"Pz",kOpenDiamond,kGreen+2);

    TF1 *fqgl = new TF1(Form("fqll_%d",j),
			"2.5626*x^3 - 3.2240*x^2 + 1.8687*x + 0.6770",0,1);
    fqgl->Draw("SAME");

    tex->SetTextColor(kRed);
    tex->SetTextSize(0.08);
    tex->DrawLatex(0.17,0.80,"---- EOY16 QGL shape for gluons");
    tex->SetTextColor(kBlack);
    tex->SetTextSize(0.045);

    c2->cd(2);
    l->DrawLine(0,1,1,1);

    if (debug) cout << "Draw c2 pad2" << endl << flush;

    TH1D *h1qu = (TH1D*)h1qm->Clone(Form("h1qu_%d",j));
    h1qu->Divide(h1qum);
    TH1D *h1qw = (TH1D*)h1qdm->Clone(Form("h1qw_%d",j));
    h1qw->Divide(h1qum);
    TH1D *h1qs = (TH1D*)h1qsm->Clone(Form("h1qs_%d",j));
    h1qs->Divide(h1qum);
    TH1D *h1qc = (TH1D*)h1qcm->Clone(Form("h1qc_%d",j));
    h1qc->Divide(h1qum);
    TH1D *h1qb = (TH1D*)h1qbm->Clone(Form("h1qb_%d",j));
    h1qb->Divide(h1qum);
    TH1D *h1qg = (TH1D*)h1qgm->Clone(Form("h1qg_%d",j));
    //h1qg->Divide(h1qum);
    h1qg->Divide(h1ggm);

    TH1D *h1gg = (TH1D*)h1gm->Clone(Form("h1gg_%d",j));
    h1gg->Divide(h1ggm);
    TH1D *h1gu = (TH1D*)h1gum->Clone(Form("h1gu_%d",j));
    h1gu->Divide(h1qum);
    TH1D *h1gw = (TH1D*)h1gdm->Clone(Form("h1gd_%d",j));
    h1gw->Divide(h1qum);
    TH1D *h1gs = (TH1D*)h1gsm->Clone(Form("h1gs_%d",j));
    h1gs->Divide(h1qum);
    TH1D *h1gc = (TH1D*)h1gcm->Clone(Form("h1gc_%d",j));
    h1gc->Divide(h1qum);
    TH1D *h1gb = (TH1D*)h1gbm->Clone(Form("h1gb_%d",j));
    h1gb->Divide(h1qum);

    if (false) {
      //tdrDraw(h1qu,"Pz",kOpenCircle,kBlue);
      tdrDraw(h1qu,"HIST",kNone,kBlue,kSolid,-1,kNone,kNone);
      tdrDraw(h1qw,"Pz",kOpenCircle,kCyan+2);
      tdrDraw(h1qs,"Pz",kOpenCircle,kMagenta+2);
      tdrDraw(h1qc,"Pz",kOpenCircle,kGreen+2);    
      tdrDraw(h1qb,"Pz",kOpenCircle,kRed);
      tdrDraw(h1qg,"Pz",kOpenCircle,kOrange+2);    
    }

    //tdrDraw(h1gg,"Pz",kOpenDiamond,kOrange+2);
    tdrDraw(h1gg,"HIST",kNone,kOrange+2,kSolid,-1,kNone,kNone);
    tdrDraw(h1gu,"Pz",kOpenDiamond,kBlue);
    tdrDraw(h1gw,"Pz",kOpenDiamond,kCyan+2);
    tdrDraw(h1gs,"Pz",kOpenDiamond,kMagenta+2);
    tdrDraw(h1gc,"Pz",kOpenDiamond,kGreen+2);
    tdrDraw(h1gb,"Pz",kOpenDiamond,kRed);

    fqgl->Draw("SAME");

    //if (pt<80 || pt>110) delete c1;
    //else c1->SaveAs("pdf/hadW_QGL_80_110.pdf");
    //if (pt<97 || pt>114) delete c1;
    //else c1->SaveAs("pdf/hadW_QGL_97_114.pdf");
    if (pt<84 || pt>97) {
      delete c1;
      delete c2;
    }
    else {
      c1->SaveAs("pdf/hadW_QGL_84_97.pdf");
      c2->SaveAs("pdf/hadW_QGL_84_97_f.pdf");
    }

  } // for j


  TH1D *hu5 = tdrHist("hu5","Efficiency (QGL>0.5)",0.2+1e-5,1.1,"p_{T} (GeV)",
		      30,200);
  //TH1D *hd5 =tdrHist("hd5","Data/MC",0.95,1.75,"p_{T} (GeV)",30,200);
  TH1D *hd5 = tdrHist("hd5","Data/MC",0.8,2.0,"p_{T} (GeV)",30,200);

  TH1D *hu7 = tdrHist("hu7","Efficiency (QGL>0.5)",0.2+1e-5,1.1,"p_{T} (GeV)",
		      30,200);
  //TH1D *hd7 =tdrHist("hd7","Data/Sim",0.95,1.75,"p_{T} (GeV)",30,200);
  //TH1D *hd7 = tdrHist("hd7","Data/Sim",0.8,2.0,"p_{T} (GeV)",30,200);
  TH1D *hd7 = tdrHist("hd7","Data/Sim",0.95,1.05,"p_{T} (GeV)",30,200);

  TH1D *hu6 = tdrHist("hu6","Efficiency (QGL>0.5)",0.2+1e-5,1.1,"p_{T} (GeV)",
		      30,200);
  TH1D *hd6 = tdrHist("hd6","Data/MC",0.8,2.0,"p_{T} (GeV)",30,200);

  TH1D *hu3 = tdrHist("hu3","Efficiency (QGL>0.5)",0.2,1.1,"p_{T} (GeV)",
		      30,200);
  TH1D *hd3 = tdrHist("hd3","Quark/Gluon",1.0,4.0,"p_{T} (GeV)",30,200);

  //TH1D *hu8 = tdrHist("hu8","Efficiency (QGL>0.5)",0.35,0.95,"p_{T} (GeV)",
  TH1D *hu8 = tdrHist("hu8","Efficiency (QGL>0.5)",0.2,1.1,"p_{T} (GeV)",
		      30,200);
  //TH1D *hd8 = tdrHist("hd8","Quark/Gluon",1.0,2.0,"p_{T} (GeV)",30,200);
  TH1D *hd8 = tdrHist("hd8","Quark/Gluon",1.0,4.0,"p_{T} (GeV)",30,200);

  TH1D *hu4 = tdrHist("hu4","Flavor fraction",0.5,1.25,"p_{T} (GeV)",30,200);
  TH1D *hd4 = tdrHist("hd4","|Quark-Gluon|",0.0,1.0,"p_{T} (GeV)",30,200);

  if (!plotMC) {
    //hu->SetMinimum(0.35);
    //hu->SetMaximum(1.10);
  }

  //lumi_13TeV = "[TT l+jet, W>qq'] UL17+18, 101.4 fb^{-1}";
  lumi_13TeV = "UL17+18, 101.4 fb^{-1}";
  TCanvas *c3 = tdrDiCanvas("c3",hu3,hd3,4,11);
  TCanvas *c8 = tdrDiCanvas("c8",hu8,hd8,4,11);
  TCanvas *c4 = tdrDiCanvas("c4",hu4,hd4,4,11);
  TCanvas *c5 = tdrDiCanvas("c5",hu5,hd5,4,11);
  TCanvas *c7 = tdrDiCanvas("c7",hu7,hd7,4,11);
  TCanvas *c6 = tdrDiCanvas("c6",hu6,hd6,4,11);

  TProfile *peqq = (TProfile*)fm->Get("peqq"); assert(peqq);
  TProfile *pequ = (TProfile*)fm->Get("pequ"); assert(pequ);
  TProfile *peqg = (TProfile*)fm->Get("peqg"); assert(peqg);
  TProfile *pegg = (TProfile*)fm->Get("pegg"); assert(pegg);
  TProfile *pegu = (TProfile*)fm->Get("pegu"); assert(pegu);
  TProfile *pegq = (TProfile*)fm->Get("pegq"); assert(pegq);
  TH1D *hequ = pequ->ProjectionX("hequ");
  TH1D *heqq = peqq->ProjectionX("heqq");
  TH1D *heqg = peqg->ProjectionX("heqg");
  TH1D *hegg = pegg->ProjectionX("hegg");
  TH1D *hegu = pegu->ProjectionX("hegu");
  TH1D *hegq = pegq->ProjectionX("hegq");

  TProfile *pfqq = (TProfile*)fm->Get("pfqq"); assert(pfqq);
  TProfile *pfqg = (TProfile*)fm->Get("pfqg"); assert(pfqg);
  TProfile *pfgg = (TProfile*)fm->Get("pfgg"); assert(pfgg);
  TProfile *pfgq = (TProfile*)fm->Get("pfgq"); assert(pfgq);
  TH1D *hfqq = pfqq->ProjectionX("hfqq");
  TH1D *hfqg = pfqg->ProjectionX("hfqg");
  TH1D *hfgg = pfgg->ProjectionX("hfgg");
  TH1D *hfgq = pfgq->ProjectionX("hfgq");


  // Infer (solve) pure quark and gluon responses from mixed ones
  // Assume MC purities. Correct for lower quark efficiency in gluon sample
  // (cause of this in low efficiency not yet understood, quite big effect)
  // Also correct for potentially different gluon efficiency in quark sample
  // (relatively small effect)
  // And correct for sizable b-jet contamination in gluon sample at low pT
  TH1D *hqemx = (TH1D*)hqem->Clone("hqemx"); //hqemx->Reset();
  TH1D *hgemx = (TH1D*)hgem->Clone("hgemx"); //hgemx->Reset();
  TH1D *hqedx = (TH1D*)hqed->Clone("hqedx"); //hqemx->Reset();
  TH1D *hgedx = (TH1D*)hged->Clone("hgedx"); //hgemx->Reset();
  for (int i = 1; i != hqem->GetNbinsX()+1; ++i) {
    double pq = hfqq->GetBinContent(i);
    double pg = hfgg->GetBinContent(i);
    double eqm = hqem->GetBinContent(i);
    double egm = hgem->GetBinContent(i);
    double kq = hegq->GetBinContent(i)/heqq->GetBinContent(i);
    double kg = heqg->GetBinContent(i)/hegg->GetBinContent(i);
    // eqm =     pq *  eq + (1-pq)*kg*eg
    // egm = (1-pg)*kq*eq +     pg *  eg
    // => pg*eqm - (1-pq)*kg*egm = pg*pq*eq - (1-pq)*kg*(1-pg)*kq*eq
    double eq = (pg*eqm - (1-pq)*kg*egm) / (pg*pq - (1-pq)*kg*(1-pg)*kq);
    // => (1-pg)*kq*eqm - pq*egm = (1-pg)*kq*(1-pq)*kg*eg - pq*pg*eg
    double eg = ((1-pg)*kq*eqm - pq*egm) / ((1-pg)*kq*(1-pq)*kg - pq*pg);
    hqemx->SetBinContent(i, eq);
    hgemx->SetBinContent(i, eg);

    double eqd = hqed->GetBinContent(i);
    double egd = hged->GetBinContent(i);
    double eqx = (pg*eqd - (1-pq)*kg*egd) / (pg*pq - (1-pq)*kg*(1-pg)*kq);
    double egx = ((1-pg)*kq*eqd - pq*egd) / ((1-pg)*kq*(1-pq)*kg - pq*pg);
    hqedx->SetBinContent(i, eqx);
    hgedx->SetBinContent(i, egx);
  } // for i

  TH1D *hqerx = (TH1D*)hqedx->Clone("hqerx");
  hqerx->Divide(hqemx);
  TH1D *hgerx = (TH1D*)hgedx->Clone("hgerx");
  hgerx->Divide(hgemx);


  // Calculate expected dijet and Z+jet efficiencies using QGL SF from TT
  // xxx  
  //TH1D* hqes = (TH1D*)hqem->Clone("hqes");
  TH1D* hzes = (TH1D*)hzem->Clone("hzes");
  TH1D* hjes = (TH1D*)hjem->Clone("hjes");
  TH1D* hges = (TH1D*)hgem->Clone("hgses");
  const int ns = 4;
  TH1D* hs[ns][6] = {{hqem,hqed,hfqq,hfqg,heqq,heqg},
		     {hgem,hged,hfgq,hfgg,hegq,hegg},
		     {hzem,hzed,hfzq,hfzg,hzqm,hzgm},
		     {hjem,hjed,hfjq,hfjg,hjemq,hjemg}};
  const char *ss[ns] = {"q","g","z","j"};
  map<string,TH1D*> mhs;
  map<string,TH1D*> mhsq;
  map<string,TH1D*> mhsg;
  map<string,TH1D*> mhr;
  for (int n = 0; n != ns; ++n) {
    /*
    TH1D *hes = hqes;
    TH1D *hfq = hfqq;
    TH1D *hfg = hfqg;
    TH1D *heq = heqq;
    TH1D *heg = heqg;
    TH1D *hem = hqem;
    */
    //int n = 1;//0;
    TH1D *hem(hs[n][0]), *hed(hs[n][1]);
    TH1D *hfq(hs[n][2]), *hfg(hs[n][3]), *heq(hs[n][4]), *heg(hs[n][5]);

    assert(hem->GetNbinsX()==hed->GetNbinsX());
    //assert(hem->GetNbinsX()==hfq->GetNbinsX());
    assert(hfq->GetNbinsX()==hfg->GetNbinsX());
    assert(hfg->GetNbinsX()==heq->GetNbinsX());
    assert(heq->GetNbinsX()==heg->GetNbinsX());

    const char *s = ss[n];
    //TH1D *hes = (TH1D*)hem->Clone(Form("h%ses",s));
    TH1D *hes = (TH1D*)heq->Clone(Form("h%ses",s));
    mhs[s] = hes;
    TH1D *hesq = (TH1D*)heq->Clone(Form("h%sesq",s));
    mhsq[s] = hesq;
    TH1D *hesg = (TH1D*)heg->Clone(Form("h%sesg",s));
    mhsg[s] = hesg;

    //for (int i = 1; i != hem->GetNbinsX()+1; ++i) {
    //double pt = hem->GetBinCenter(i);
    //int j = hfq->FindBin(pt);
    for (int j = 1; j != heq->GetNbinsX()+1; ++j) {
      double pt = heq->GetBinCenter(j);
      int i = hem->FindBin(pt);

      double pq = hfq->GetBinContent(j);
      double pg = hfg->GetBinContent(j);
      double eq = heq->GetBinContent(j);
      double eg = heg->GetBinContent(j);
      double ex = hem->GetBinContent(i);
      double kq = 1.084;
      double kg = 1.723;
      
      // Correct efficiency assuming MC fractions
      double es = ex + eq*(kq-1)*pq + eg*(kg-1)*pg;
      double esq = eq*kq;
      double esg = eg*kg;
      // Except for dijet: gluon JES+FSR changes reconstructed gluon fraction
      // by 25% for 5% change in JES (1-2%) and FSR (3-4%)
      if (string(s)=="j" && true) {
	double pqd = pq + 0.25*pg;
	double pgd = pg - 0.25*pg;
	es = ex + (pqd-pq)*eq + (pgd-pg)*eg + eq*(kq-1)*pqd + eg*(kg-1)*pgd;
	//esq = eq*kq;
	//esg = eg*kg;
      }
      // Assume Z+jet MC flavor effiencies are 5% higher than W>qq' and dijet
      // due to some extra simulation artefact, and are the same in data
      // Maybe gluon-tagged quark jets radiate more FSR?
      if (string(s)=="z" && true) {
	double rz = 0.95;
	double kqd = kq*rz;
	double kgd = kg*rz;
	es = ex + eq*(kqd-1)*pq + eg*(kgd-1)*pg;
	esq = eq*kqd;
	esg = eg*kgd;
      }

      //hes->SetBinContent(i, es);
      //hesq->SetBinContent(i, esq);
      //hesg->SetBinContent(i, esg);
      hes->SetBinContent(j, es);
      hesq->SetBinContent(j, esq);
      hesg->SetBinContent(j, esg);
    } // for i

    TH1D *hrs = (TH1D*)hed->Clone(Form("h%srs",s));
    hrs->Divide(hes);
    mhr[s] = hrs;
  } 
  //TH1D *hqs = (TH1D*)hqed->Clone("hqs");
  //hqs->Divide(hqes);

  c5->cd(1);
  gPad->SetLogx();

  TLegend *leg5 = tdrLeg(0.65,0.90-5*0.05,0.85,0.90);

  tex->DrawLatex(0.19,0.70,"|#eta| < 2.5");
  if (!plotMC) tex->DrawLatex(0.19,0.65,"Data only");

  l->DrawLine(30,0.5,200,0.5);

  tdrDraw(hzed,"Pz",kFullSquare,kMagenta+2);
  if (plotMC) tdrDraw(hzem,"Pz",kOpenDiamond,kMagenta+2);
  tdrDraw(hjed,"Pz",kFullStar,kGreen+2);
  if (plotMC) tdrDraw(hjem,"Pz",kOpenDiamond,kGreen+2);

  tdrDraw(hbed,"Pz",kFullCross,kRed);
  if (plotMC) tdrDraw(hbem,"Pz",kOpenSquare,kRed);
  tdrDraw(hqed,"Pz",kFullCircle,kBlue);
  if (plotMC) tdrDraw(hqem,"Pz",kOpenCircle,kBlue);
  tdrDraw(hged,"Pz",kFullDiamond,kOrange+2);
  if (plotMC) tdrDraw(hgem,"Pz",kOpenDiamond,kOrange+2);

  leg5->AddEntry(hqed,"W>qq' (q rich)","PLE");
  leg5->AddEntry(hzed,"Z+jet   (q rich)","PLE");
  leg5->AddEntry(hbed,"t>Wb   (b rich)","PLE");
  leg5->AddEntry(hjed,"Dijet    (g rich)","PLE");
  leg5->AddEntry(hged,"TT+jet (g rich)","PLE");

  c5->cd(2);
  gPad->SetLogx();

  //l->DrawLine(30,1.10,200,1.10);
  l->DrawLine(30,1.0,200,1.0);
  //tdrDraw(hzq, "Pz",kFullSquare,kMagenta+2);
  tdrDraw(hbe, "Pz",kFullCross,kRed);
  tdrDraw(hqe, "Pz",kFullCircle,kBlue);
  tdrDraw(hge, "Pz",kFullDiamond,kOrange+2);

  tdrDraw(hze, "Pz",kFullSquare,kMagenta+2);
  tdrDraw(hje, "Pz",kFullStar,kGreen+2);

  c5->SaveAs("pdf/hadW_QGL_DataVsPt.pdf");

  c7->cd(1);
  gPad->SetLogx();

  TLegend *leg7 = tdrLeg(0.65,0.90-4*0.05,0.85,0.90);

  tex->DrawLatex(0.19,0.70,"|#eta| < 2.5");
  tex->DrawLatex(0.19,0.65,"Sim only");

  l->DrawLine(30,0.5,200,0.5);

  tex->DrawLatex(0.19,0.20,"Simulation settings:");
  tex->DrawLatex(0.19,0.15,"All samples, eff(g)#times1.723, eff(q)#times1.084");
  tex->DrawLatex(0.19,0.10,"+ for Z+jet, eff(g')#times0.95, eff(q')#times0.95");
  tex->DrawLatex(0.19,0.05,"+ for dijet, frac(g)#times0.75");

  //tdrDraw(hzes,"Pz",kFullSquare,kMagenta+2);
  //tdrDraw(hjes,"Pz",kFullStar,kGreen+2);
  //tdrDraw(hbed,"Pz",kFullCross,kRed);
  //tdrDraw(hqes,"Pz",kFullCircle,kBlue);
  //tdrDraw(hges,"Pz",kFullDiamond,kOrange+2);

  //leg7->AddEntry(hqes,"W>qq' (q rich)","PLE");
  //leg7->AddEntry(hzes,"Z+jet   (q rich)","PLE");
  //leg7->AddEntry(hbes,"t>Wb   (b rich)","PLE");
  //leg7->AddEntry(hjes,"Dijet    (g rich)","PLE");
  //leg7->AddEntry(hges,"TT+jet (g rich)","PLE");

  tdrDraw(mhs["z"],"Pz",kFullSquare,kMagenta+2);
  tdrDraw(mhs["j"],"Pz",kFullStar,kGreen+2);
  tdrDraw(mhs["q"],"Pz",kFullCircle,kBlue);
  tdrDraw(mhs["g"],"Pz",kFullDiamond,kOrange+2);

  leg7->AddEntry(mhs["q"],"W>qq' (q rich)","PLE");
  leg7->AddEntry(mhs["z"],"Z+jet   (q rich)","PLE");
  leg7->AddEntry(mhs["j"],"Dijet    (g rich)","PLE");
  leg7->AddEntry(mhs["g"],"TT+jet (g rich)","PLE");

  c7->cd(2);
  gPad->SetLogx();

  //l->DrawLine(30,1.10,200,1.10);
  l->DrawLine(30,1.0,200,1.0);

  //tdrDraw(hbe, "Pz",kFullCross,kRed);
  //tdrDraw(hqs, "Pz",kFullCircle,kBlue);
  /*
  tdrDraw(hgs, "Pz",kFullDiamond,kOrange+2);
  tdrDraw(hzs, "Pz",kFullSquare,kMagenta+2);
  tdrDraw(hjs, "Pz",kFullStar,kGreen+2);
  */

  tdrDraw(mhr["q"], "Pz",kFullCircle,kBlue);
  tdrDraw(mhr["g"], "Pz",kFullDiamond,kOrange+2);
  tdrDraw(mhr["z"], "Pz",kFullSquare,kMagenta+2);
  tdrDraw(mhr["j"], "Pz",kFullStar,kGreen+2);

  c7->SaveAs("pdf/hadW_QGL_SimVsPt.pdf");


  c6->cd(1);
  gPad->SetLogx();

  TLegend *leg6 = tdrLeg(0.65,0.90-6*0.05,0.85,0.90);

  tex->DrawLatex(0.19,0.70,"|#eta| < 2.5");
  tex->DrawLatex(0.19,0.65,"MC only");

  l->DrawLine(30,0.5,200,0.5);
  //tex->DrawLatex(0.50,0.33,"#uparrow Quark jets (udsc)");
  //tex->DrawLatex(0.50,0.26,"#downarrow Gluon jets (g)");
  tex->DrawLatex(0.35,0.33,"#uparrow Quark jets (udsc) --- W>qq'");
  tex->DrawLatex(0.35,0.26,"#downarrow Gluon jets (g) --- TT+jet");
  
  // MC reco, truth and "inferred" references
  tdrDraw(hqem,"Pz",kOpenCircle,kBlue);
  tdrDraw(hgem,"Pz",kOpenDiamond,kOrange+2);
  tdrDraw(hqemx,"Pz",kFullSquare,kBlue);
  tdrDraw(hgemx,"Pz",kFullSquare,kOrange+2);
  tdrDraw(heqq,"Pz",kOpenCross,kBlue+2);
  tdrDraw(hegg,"Pz",kOpenCross,kOrange+4);
  /*
  tdrDraw(hqed,"Pz",kFullCircle,kBlue);
  tdrDraw(hged,"Pz",kFullCircle,kOrange+2);
  tdrDraw(hqedx,"Pz",kFullSquare,kBlue);
  tdrDraw(hgedx,"Pz",kFullSquare,kOrange+2);
  */

  //tex->DrawLatex(0.42,0.77,"Quark jets #uparrow");
  //tex->DrawLatex(0.42,0.71,"Gluon jets #downarrow");
  leg6->AddEntry(hqemx,"Inferred q","PLE");
  leg6->AddEntry(heqq,"True q","PLE");
  leg6->AddEntry(hqem,"Reco q rich","PLE");
  leg6->AddEntry(hgem,"Reco g rich","PLE");
  leg6->AddEntry(hegg,"True g","PLE");
  leg6->AddEntry(hgemx,"Inferred g","PLE");

  c6->cd(2);
  gPad->SetLogx();
  l->DrawLine(30,1.0,200,1.0);

  tdrDraw(hqe,"Pz",kOpenCircle,kBlue);
  tdrDraw(hge,"Pz",kOpenDiamond,kOrange+2);
  tdrDraw(hqerx,"Pz",kFullSquare,kBlue);
  tdrDraw(hgerx,"Pz",kFullSquare,kOrange+2);

  TF1 *fg = new TF1("fg","[0]",32,200);
  fg->SetLineColor(kOrange+2);
  TF1 *fq = new TF1("fq","[0]",32,200);
  fq->SetLineColor(kBlue);

  hqerx->Fit(fq,"QRN");
  hgerx->Fit(fg,"QRN");

  fq->Draw("SAME");
  fg->Draw("SAME");

  tex->SetTextSize(1.6*0.045);
  tex->SetTextColor(kOrange+2);
  //tex->DrawLatex(0.30,0.77,Form("Eff(g) = %1.3f #pm %1.3f" // 2.2
  tex->DrawLatex(0.30,0.84,Form("Eff(g) = %1.3f #pm %1.3f"
				" (#chi^{2}/NDF = %1.1f / %d)",
				fg->GetParameter(0), fg->GetParError(0),
				fg->GetChisquare(), fg->GetNDF()));
  tex->SetTextColor(kBlue);
  //tex->DrawLatex(0.30,0.48,Form("Eff(q) = %1.3f #pm %1.3f" // 2.2
  tex->DrawLatex(0.30,0.52,Form("Eff(q) = %1.3f #pm %1.3f"
				" (#chi^{2}/NDF = %1.1f / %d)",
				fq->GetParameter(0), fq->GetParError(0),
				fq->GetChisquare(), fq->GetNDF()));
  tex->SetTextColor(kBlack);
  tex->SetTextSize(0.045);

  c6->SaveAs("pdf/hadW_QGL_MCvsPt.pdf");

  c3->cd(1);
  gPad->SetLogx();
  l->DrawLine(30,0.5,200,0.5);

  TLegend *leg3 = tdrLeg(0.65,0.90-4*0.05,0.85,0.90);

  tex->DrawLatex(0.19,0.70,"|#eta| < 2.5");
  tex->DrawLatex(0.19,0.65,"MC only");
  tex->DrawLatex(0.50,0.33,"#uparrow Quark jets (udsc)");
  tex->DrawLatex(0.50,0.26,"#downarrow Gluon jets (g)");
  
  hzqm->Scale(0.95);
  hzgm->Scale(0.95);
  tdrDraw(hzqm,"Pz",kOpenSquare,kMagenta+2);
  tdrDraw(hzgm,"Pz",kOpenSquare,kMagenta+2);

  tdrDraw(hjemq,"Pz",kOpenStar,kGreen+2);
  tdrDraw(hjemg,"Pz",kOpenStar,kGreen+2);

  //tdrDraw(hequ,"Pz",kOpenCircle,kBlue+2);//extra
  tdrDraw(heqq,"Pz",kOpenCircle,kBlue);
  tdrDraw(heqg,"Pz",kOpenCircle,kBlue);
  tdrDraw(hegg,"Pz",kOpenDiamond,kOrange+2);
  //tdrDraw(hegu,"Pz",kOpenDiamond,kOrange+4);//extra
  tdrDraw(hegq,"Pz",kOpenDiamond,kOrange+2);
  
  leg3->AddEntry(hzqm,"Z+jet (#times0.95)","PLE");
  leg3->AddEntry(heqq,"W>qq'","PLE");
  leg3->AddEntry(hjemq,"Dijet","PLE");
  leg3->AddEntry(hegq,"TT+jet","PLE");

  // Fit MC quark and gluon efficiencies
  TMultiGraph *mgmq = new TMultiGraph();
  //mgmq->Add(new TGraphErrors(hzqm));//mhsq["z"]));
  mgmq->Add(new TGraphErrors(heqq));//mhsq["q"]));
  mgmq->Add(new TGraphErrors(hjemq));//mhsq["j"]));

  TMultiGraph *mgmg = new TMultiGraph();
  //mgmg->Add(new TGraphErrors(hzgm));//mhsg["z"]));
  mgmg->Add(new TGraphErrors(hjemg));//mhsg["j"]));
  mgmg->Add(new TGraphErrors(hegg));//mhsg["g"]));
  
  TF1 *f1mq = new TF1("f1mq","[0]+[1]*pow(x,[2])",34,200);
  f1mq->SetParameters(0.8,-0.5,-0.5);
  mgmq->Fit(f1mq,"QRNW");
  f1mq->Draw("SAME");

  TF1 *f1mg = new TF1("f1mg","[0]+[1]*pow(x,[2])",34,200);
  f1mg->SetParameters(0.25,+0.5,-0.5);
  mgmg->Fit(f1mg,"QRNW");
  f1mg->Draw("SAME");

  TF1 *f1m = new TF1("f1m","([0]+[1]*pow(x,[2]))/([3]+[4]*pow(x,[5]))",34,200);
  f1m->SetParameters(f1mq->GetParameter(0),f1mq->GetParameter(1),
		     f1mq->GetParameter(2),f1mg->GetParameter(0),
		     f1mg->GetParameter(1),f1mg->GetParameter(2));

  c3->cd(2);
  gPad->SetLogx();

  TH1D *hzrm = (TH1D*)hzqm->Clone("hzrm");
  hzrm->Divide(hzgm);
  TH1D *hjrm = (TH1D*)hjemq->Clone("hjrm");
  hjrm->Divide(hjemg);
  TH1D *hqrm = (TH1D*)heqq->Clone("hqrm");
  hqrm->Divide(heqg);
  TH1D *hgrm = (TH1D*)hegq->Clone("hgrm");
  hgrm->Divide(hegg);

  tdrDraw(hzrm,"Pz",kOpenSquare,kMagenta+2);
  tdrDraw(hjrm,"Pz",kOpenStar,kGreen+2);
  tdrDraw(hqrm,"Pz",kOpenCircle,kBlue);
  tdrDraw(hgrm,"Pz",kOpenDiamond,kOrange+2);

  f1m->Draw("SAME");

  c3->SaveAs("pdf/hadW_QGL_MCefficiencies.pdf");

  c8->cd(1);
  gPad->SetLogx();
  l->DrawLine(30,0.5,200,0.5);

  TLegend *leg8 = tdrLeg(0.65,0.90-4*0.05,0.85,0.90);

  tex->DrawLatex(0.19,0.70,"|#eta| < 2.5");
  tex->DrawLatex(0.19,0.65,"Sim only");
  //tex->DrawLatex(0.50,0.33,"#uparrow Quark jets (udsc)");
  //tex->DrawLatex(0.50,0.26,"#downarrow Gluon jets (g)");
  tex->DrawLatex(0.50,0.48,"#uparrow Quark jets (udsc)");
  tex->DrawLatex(0.50,0.38,"#downarrow Gluon jets (g)");
  
  tdrDraw(mhsq["z"],"Pz",kOpenSquare,kMagenta+2);
  tdrDraw(mhsg["z"],"Pz",kOpenSquare,kMagenta+2);
  tdrDraw(mhsq["j"],"Pz",kOpenStar,kGreen+2);
  tdrDraw(mhsg["j"],"Pz",kOpenStar,kGreen+2);

  tdrDraw(mhsq["q"],"Pz",kOpenCircle,kBlue);
  tdrDraw(mhsg["q"],"Pz",kOpenCircle,kBlue);
  tdrDraw(mhsq["g"],"Pz",kOpenDiamond,kOrange+2);
  tdrDraw(mhsg["g"],"Pz",kOpenDiamond,kOrange+2);
 
  leg8->AddEntry(mhsq["z"],"Z+jet","PLE");
  leg8->AddEntry(mhsq["q"],"W>qq'","PLE");
  leg8->AddEntry(mhsq["j"],"Dijet","PLE");
  leg8->AddEntry(mhsq["g"],"TT+jet","PLE");

  // Fit Sim quark and gluon efficiencies
  TMultiGraph *mgsq = new TMultiGraph();
  //mgsq->Add(new TGraphErrors(mhsq["z"]));
  mgsq->Add(new TGraphErrors(mhsq["q"]));
  mgsq->Add(new TGraphErrors(mhsq["j"]));

  TMultiGraph *mgsg = new TMultiGraph();
  //mgsg->Add(new TGraphErrors(mhsg["z"]));
  mgsg->Add(new TGraphErrors(mhsg["j"]));
  mgsg->Add(new TGraphErrors(mhsg["g"]));
  
  TF1 *f1sq = new TF1("f1sq","[0]+[1]*pow(x,[2])",34,200);
  f1sq->SetParameters(0.8,-0.5,-0.5);
  mgsq->Fit(f1sq,"QRNW");
  f1sq->Draw("SAME");

  TF1 *f1sg = new TF1("f1sg","[0]+[1]*pow(x,[2])",34,200);
  f1sg->SetParameters(0.45,+0.5,-0.5);
  mgsg->Fit(f1sg,"QRNW");
  f1sg->Draw("SAME");

  c8->cd(2);
  gPad->SetLogx();


  if (false) {

    // Plot Quark/Gluon for Sim
    TH1D *hzrs = (TH1D*)mhsq["z"]->Clone("hzrs");
    hzrs->Divide(mhsg["z"]);
    TH1D *hjrs = (TH1D*)mhsq["j"]->Clone("hjrs");
    hjrs->Divide(mhsg["j"]);
    TH1D *hqrs = (TH1D*)mhsq["q"]->Clone("hqrs");
    hqrs->Divide(mhsg["q"]);
    TH1D *hgrs = (TH1D*)mhsq["g"]->Clone("hgrs");
    hgrs->Divide(mhsg["g"]);
    
    tdrDraw(hzrs,"Pz",kOpenSquare,kMagenta+2);
    tdrDraw(hjrs,"Pz",kOpenStar,kGreen+2);
    tdrDraw(hqrs,"Pz",kOpenCircle,kBlue);
    tdrDraw(hgrs,"Pz",kOpenDiamond,kOrange+2);
    
    TF1 *f1s = new TF1("f1s","([0]+[1]*pow(x,[2]))/([3]+[4]*pow(x,[5]))",34,200);
    f1s->SetParameters(f1sq->GetParameter(0),f1sq->GetParameter(1),
		       f1sq->GetParameter(2),f1sg->GetParameter(0),
		       f1sg->GetParameter(1),f1sg->GetParameter(2));
    
    f1s->Draw("SAME");
  }
  else {
  /*
  // These are plotted for MC efficiencies
  tdrDraw(hzqm,"Pz",kOpenSquare,kMagenta+2);
  tdrDraw(hzgm,"Pz",kOpenSquare,kMagenta+2);

  tdrDraw(hjemq,"Pz",kOpenStar,kGreen+2);
  tdrDraw(hjemg,"Pz",kOpenStar,kGreen+2);

  //tdrDraw(hequ,"Pz",kOpenCircle,kBlue+2);//extra
  tdrDraw(heqq,"Pz",kOpenCircle,kBlue);
  tdrDraw(heqg,"Pz",kOpenCircle,kBlue);
  tdrDraw(hegg,"Pz",kOpenDiamond,kOrange+2);
  //tdrDraw(hegu,"Pz",kOpenDiamond,kOrange+4);//extra
  tdrDraw(hegq,"Pz",kOpenDiamond,kOrange+2);
  */

    hd8->SetMinimum(0.9);
    hd8->SetMaximum(1.9);
    hd8->SetYTitle("Sim/MC");

    // Plot Sim/MC for SimEfficiencies for quarks and gluons
    TH1D *hzrs = (TH1D*)mhsq["z"]->Clone("hzrs");
    hzrs->Divide(hzqm);
    TH1D *hjrs = (TH1D*)mhsg["j"]->Clone("hjrs");
    hjrs->Divide(hjemg);
    TH1D *hqrs = (TH1D*)mhsq["q"]->Clone("hqrs");
    hqrs->Divide(heqq);
    TH1D *hgrs = (TH1D*)mhsg["g"]->Clone("hgrs");
    hgrs->Divide(hegg);

    tdrDraw(hzrs,"Pz",kOpenSquare,kMagenta+2);
    tdrDraw(hjrs,"Pz",kOpenStar,kGreen+2);
    tdrDraw(hqrs,"Pz",kOpenCircle,kBlue);
    tdrDraw(hgrs,"Pz",kOpenDiamond,kOrange+2);
  }

  c8->SaveAs("pdf/hadW_QGL_SimEfficiencies.pdf");
  
  
  c4->cd(1);
  gPad->SetLogx();

  TLegend *leg4 = tdrLeg(0.65,0.90-4*0.05,0.85,0.90);

  tex->DrawLatex(0.19,0.70,"|#eta| < 2.5");
  tex->DrawLatex(0.19,0.65,"MC only");

  l->DrawLine(30,1.0,200,1.0);
  l->DrawLine(30,0.5,200,0.5);

  tdrDraw(hfqq,"Pz",kOpenCircle,kBlue);
  tdrDraw(hfzq,"Pz",kOpenSquare,kMagenta+2);
  //tdrDraw(hfjq,"Pz",kOpenStar,kGreen+2);
  //tdrDraw(hfgq,"Pz",kOpenDiamond,kOrange+2);

  //tdrDraw(hfqg,"Pz",kOpenCircle,kBlue);
  //tdrDraw(hfzg,"Pz",kOpenSquare,kMagenta+2);
  tdrDraw(hfjg,"Pz",kOpenStar,kGreen+2);
  tdrDraw(hfgg,"Pz",kOpenDiamond,kOrange+2);

  leg4->AddEntry(hzqm,"Z+jet (q)","PLE");
  leg4->AddEntry(heqq,"W>qq' (q)","PLE");
  leg4->AddEntry(hjemq,"Dijet (g)","PLE");
  leg4->AddEntry(hegq,"TT+jet (g)","PLE");

  c4->cd(2);
  gPad->SetLogx();

  TH1D *hdfq = (TH1D*)hfqq->Clone("hdfq");
  hdfq->Add(hfqg,-1);
  TH1D *hdfz = (TH1D*)hfzq->Clone("hdfz");
  hdfz->Add(hfzg,-1);
  TH1D *hdfj = (TH1D*)hfjg->Clone("hdfj");
  hdfj->Add(hfjq,-1);
  TH1D *hdfg = (TH1D*)hfgg->Clone("hdfg");
  hdfg->Add(hfgq,-1);

  tdrDraw(hdfq,"Pz",kOpenCircle,kBlue);
  tdrDraw(hdfz,"Pz",kOpenSquare,kMagenta+2);
  tdrDraw(hdfj,"Pz",kOpenStar,kGreen+2);
  tdrDraw(hdfg,"Pz",kOpenDiamond,kOrange+2);

  c4->SaveAs("pdf/hadW_QGL_MCfractions.pdf");

  // Print out QGL>0.5 efficiency fits
  cout << endl;
  cout << "  // hadW_QGL.C QGL>0.5 efficiency fits\n";
  cout << "  // For Z+jet, multiply MC results by 1.05\n";
  cout << "  TF1 *f1mq = new TF1(\"f1mq\",\"[0]+[1]*pow(x,[2])\",34,200);\n";
  cout << "  TF1 *f1mg = new TF1(\"f1mg\",\"[0]+[1]*pow(x,[2])\",34,200);\n";
  cout << "  TF1 *f1dq = new TF1(\"f1dq\",\"[0]+[1]*pow(x,[2])\",34,200);\n";
  cout << "  TF1 *f1dg = new TF1(\"f1dg\",\"[0]+[1]*pow(x,[2])\",34,200);\n";
  cout << Form("  f1mq->SetParameters(%1.4g,%1.4g,%1.4g);\n",
	       f1mq->GetParameter(0),f1mq->GetParameter(1),
	       f1mq->GetParameter(2));
  cout << Form("  f1mg->SetParameters(%1.4g,%1.4g,%1.4g);\n",
	       f1mg->GetParameter(0),f1mg->GetParameter(1),
	       f1mg->GetParameter(2));
  cout << Form("  f1dq->SetParameters(%1.4g,%1.4g,%1.4g);\n",
	       f1sq->GetParameter(0),f1sq->GetParameter(1),
	       f1sq->GetParameter(2));
  cout << Form("  f1dg->SetParameters(%1.4g,%1.4g,%1.4g);\n",
	       f1sg->GetParameter(0),f1sg->GetParameter(1),
	       f1sg->GetParameter(2));
  cout << endl;
} // hadW_QGL


void  getDijetQGL(TFile *fjd, TFile *fjm, TFile *fjd2, TFile *fjm2,
		  TH2D **h2jd, TH2D **h2jm, TH2D **h2jmq, TH2D **h2jmg) {

  double eta[] = {0,0.5,1.0,1.5,2.0,2.5};
  const int neta = sizeof(eta)/sizeof(eta[0])-1;

  //for (int ieta = 0; ieta != neta; ++ieta) {
  for (int ieta = neta-1; ieta != -1; --ieta) {

    //TH2D *h2d = (TH2D*)fjd->Get("Standard/Eta_0.0-1.3/jt0/hqgl2");
    //TH2D *h2m = (TH2D*)fjm->Get("Standard/Eta_0.0-1.3/jt0/hqgl2");
    TH2D *h2d = (TH2D*)fjd->Get(Form("Standard/Eta_%1.1f-%1.1f/jt0/hqgl2",
    				     eta[ieta],eta[ieta+1]));
    assert(h2d);
    TH2D *h2m = (TH2D*)fjm->Get(Form("Standard/Eta_%1.1f-%1.1f/jt0/hqgl2",
    				     eta[ieta],eta[ieta+1]));
    assert(h2m);
    TH2D *h2mq = (TH2D*)fjm->Get(Form("Standard/Eta_%1.1f-%1.1f/jt0/hqgl2_q",
				      eta[ieta],eta[ieta+1]));
    assert(h2mq);
    TH2D *h2mg = (TH2D*)fjm->Get(Form("Standard/Eta_%1.1f-%1.1f/jt0/hqgl2_g",
				      eta[ieta],eta[ieta+1]));
    assert(h2mg);
    TH2D *h2mu = (TH2D*)fjm->Get(Form("Standard/Eta_%1.1f-%1.1f/jt0/hqgl2_u",
				      eta[ieta],eta[ieta+1]));
    assert(h2mu);

    TH2D *h2d2 = (TH2D*)fjd2->Get(Form("Standard/Eta_%1.1f-%1.1f/jt0/hqgl2",
				      eta[ieta],eta[ieta+1]));
    assert(h2d2);
    TH2D *h2m2 = (TH2D*)fjm2->Get(Form("Standard/Eta_%1.1f-%1.1f/jt0/hqgl2",
				       eta[ieta],eta[ieta+1]));
    assert(h2m2);
    TH2D *h2mq2 = (TH2D*)fjm2->Get(Form("Standard/Eta_%1.1f-%1.1f/jt0/hqgl2_q",
					eta[ieta],eta[ieta+1]));
    assert(h2mq2);
    TH2D *h2mg2 = (TH2D*)fjm2->Get(Form("Standard/Eta_%1.1f-%1.1f/jt0/hqgl2_g",
					eta[ieta],eta[ieta+1]));
    assert(h2mg2);
    TH2D *h2mu2 = (TH2D*)fjm2->Get(Form("Standard/Eta_%1.1f-%1.1f/jt0/hqgl2_u",
					eta[ieta],eta[ieta+1]));
    assert(h2mu2);

    /*
    if (!h2jd) h2jd = h2d;
    else        h2jd->Add(h2d);
    if (!h2jm) h2jm = h2m;
    else        h2jm->Add(h2m);
    */
    if (!(*h2jd) && !(*h2jm)) {
      (*h2jd) = (TH2D*)h2d->Clone("h2jd");   (*h2jd)->Add(h2d2);
      (*h2jm) = (TH2D*)h2m->Clone("h2jm");   (*h2jm)->Add(h2m2);
      (*h2jmq) = (TH2D*)h2mq->Clone("h2mq"); (*h2jmq)->Add(h2mq2);
      (*h2jmg) = (TH2D*)h2mg->Clone("h2mg"); (*h2jmg)->Add(h2mg2);
      (*h2jmg)->Add(h2mu); (*h2jmg)->Add(h2mu2); 
      cout << "QGL bins (dijet):\n  {" << endl;
      for (int i = 1; i != h2d->GetNbinsX()+1; ++i) {
	double ptmin = h2d->GetXaxis()->GetBinLowEdge(i);
	if (ptmin>=15 && ptmin<400) {
	  cout << ptmin << ", ";
	}
      }
      cout << endl;
    }
    else { // hadd
      assert((*h2jd)->GetNbinsX()==(*h2jm)->GetNbinsX());
      assert((*h2jd)->GetNbinsX()==(*h2jmq)->GetNbinsX());
      assert((*h2jd)->GetNbinsX()==(*h2jmg)->GetNbinsX());
      assert(h2d->GetNbinsX()>=(*h2jd)->GetNbinsX());

      for (int i = 1; i != (*h2jd)->GetNbinsX()+1; ++i) {
	for (int j = 1; j != (*h2jd)->GetNbinsY()+1; ++j) {

	    (*h2jd)->SetBinContent(i,j, (*h2jd)->GetBinContent(i,j)
				    + h2d->GetBinContent(i,j)
				    + h2d2->GetBinContent(i,j));
	    (*h2jm)->SetBinContent(i,j, (*h2jm)->GetBinContent(i,j)
				    + h2m->GetBinContent(i,j)
				    + h2m2->GetBinContent(i,j));
	    (*h2jmq)->SetBinContent(i,j, (*h2jmq)->GetBinContent(i,j)
				    + h2mq->GetBinContent(i,j)
				    + h2mq2->GetBinContent(i,j));
	    (*h2jmg)->SetBinContent(i,j, (*h2jmg)->GetBinContent(i,j)
				    + h2mg->GetBinContent(i,j)
				    + h2mg2->GetBinContent(i,j)
				    + h2mu->GetBinContent(i,j)
				    + h2mu2->GetBinContent(i,j));
	} // for j
      } // for i
    } // hadd
  } // for ieta
} // getDijetQGL
