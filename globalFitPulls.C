// Purpose: plot the pulls and post-fit parameter values
//          produced by globalFitL3Res.C and
//          stored in rootfiles/globalFitL3Res_emat.root
//          in a nice Higgs Combine tool style
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH2D.h"
#include "TLine.h"

#include "tdrstyle_mod15.C"


void globalFitPulls() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TFile *f = new TFile("rootfiles/globalFitL3Res_emat.root","READ");
  assert(f && !f->IsZombie());

  curdir->cd();
  
  TH1D *hp = (TH1D*)f->Get("h1par"); assert(hp);
  TH1D *h2cov = (TH1D*)f->Get("h2cov"); assert(hp);

  // Map names to better ones for graph
  map<string,const char*> mpar;
  mpar["bm1_fsru1"] = "Multijet HDM uncl.";
  mpar["bm4_fsru1"] = "Z(ll)+jet HDM uncl.";
  mpar["bm16_fsru1"] = "Multijet HDM uncl. v2";
  mpar["bm32_hkfsr_mpfchs1_gamjet_eig0"] = "#gamma+jet MPF s_{0}";
  mpar["bm64_fsru1"] = "Z(ll)+jet HDM uncl. v2";
  mpar["bm1_fsrn1"] = "Multijet HDM jets";
  mpar["bm4_fsrn1"] = "Z(ll)+jet HDM jets";
  mpar["bm16_fsrn1"] = "Multijet HDM jets v2";
  mpar["bm32_hkfsr_mpfchs1_gamjet_eig1"] = "s_{1}";//"#gamma+jet MPF s_{1}";
  mpar["bm64_fsrn1"] = "Z(ll)+jet HDM jets v2";
  //
  mpar["bm1_hkfsr_ptchs_multijet_eig0"] = "Multijet p_{T} s_{0}";
  mpar["bm2_hkfsr_ptchs_gamjet_eig0"] = "#gamma+jet p_{T} s_{0}";
  mpar["bm4_hkfsr_ptchs_zlljet_eig0"] = "Z(ll)+jet p_{T} s_{0}";
  mpar["bm8_hkfsr_mpfchs1_multijet_eig0"] = "Multijet MPF s_{0}";
  mpar["bm16_hkfsr_mpfchs1_gamjet_eig0"] = "#gamma+jet MPF s_{0}";
  mpar["bm32_hkfsr_mpfchs1_zlljet_eig0"] = "Z(ll)+jet MPF s_{0}";
  mpar["bm1_hkfsr_ptchs_multijet_eig1"] = "s_{1}";//"multijet p_{T} s_{1}";
  mpar["bm2_hkfsr_ptchs_gamjet_eig1"] = "s_{1}";//"#gamma+jet p_{T} s_{1}";
  mpar["bm4_hkfsr_ptchs_zlljet_eig1"] = "s_{1}";//"Z(ll)+jet p_{T} s_{1}";
  mpar["bm8_hkfsr_mpfchs1_multijet_eig1"] = "s_{1}";//"multijet MPF s_{1}";
  mpar["bm16_hkfsr_mpfchs1_gamjet_eig1"] = "s_{1}";//"#gamma+jet MPF s_{1}";
  mpar["bm32_hkfsr_mpfchs1_zlljet_eig1"] = "s_{1}";//"Z(ll)+jet MPF s_{1}";
  //
  mpar["bm34_scale_050_gamjet"] = "#gamma scale (0.5%)";
  mpar["bm68_scale_010_zlljet"] = "Z(ll) scale (0.1%)";
  mpar["bm32_mpfscale_02_mpfchs1_gamjet"] ="#gamma+jet MPF scale (0.2%)";
  mpar["bm64_mpfscale_02_mpfchs1_zlljet"] = "Z(ll)+jet MPF scale (0.2%)";
  //mpar["bm9_multijet_jer_src0"] = "multijet JER";
  mpar["bm34_eesfromzee_gamjet_eig0"] = "EM scale s_{0}";
  mpar["bm34_eesfromzee_gamjet_eig1"] = "s_{1}";//"EM scale s_{1}";
  mpar["bm34_eesfromzee_gamjet_eig2"] = "s_{2}";//"EM scale s_{2}";
  //
  mpar["bm8_hadw_ptave_fitprob"] = "W>qq' fitProb s_{0}";
  mpar["bm128_hadw_ptboth_fitprob"] = "W>qq' fitProb s_{2}";
  mpar["bm8_hadw_ptave_fitprob2"] = "W>qq' fitProb s_{1}";
  mpar["bm128_hadw_ptboth_fitprob2"] = "W>qq' fitProb s_{3}";
  /*
  mpar["bm1_hkfsr_ptchs_multijet_eig2"] = "s_{2}";//"multijet p_{T} s_{2}";
  mpar["bm2_hkfsr_ptchs_gamjet_eig2"] = "s_{2}";//"#gamma+jet p_{T} s_{2}";
  mpar["bm4_hkfsr_ptchs_zlljet_eig2"] = "s_{2}";//"Z(ll)+jet p_{T} s_{2}";
  mpar["bm8_hkfsr_mpfchs1_multijet_eig2"] = "s_{2}";//"multijet MPF s_{2}";
  mpar["bm16_hkfsr_mpfchs1_gamjet_eig2"] = "s_{2}";//"#gamma+jet MPF s_{2}";
  mpar["bm32_hkfsr_mpfchs1_zlljet_eig2"] = "s_{2}";//"Z(ll)+jet MPF s_{2}";
  mpar["bm1_hkfsr_ptchs_multijet_eig3"] = "s_{3}";//"multijet p_{T} s_{3}";
  mpar["bm8_hkfsr_mpfchs1_multijet_eig3"] = "s_{3}";//"multijet MPF s_{3}";
  mpar["bm16_hkfsr_mpfchs1_gamjet_eig3"] = "s_{3}";//"#gamma+jet MPF s_{3}";
  mpar["bm32_hkfsr_mpfchs1_zlljet_eig3"] = "s_{3}";//"Z(ll)+jet MPF s_{3}";
  */
  mpar["bm18_scale_050_gamjet"] = "#gamma scale (0.5%)";
  mpar["bm36_scale_020_zlljet"] = "Z(ll) scale (0.2%)";
  mpar["bm16_mpfscale_02_mpfchs1_gamjet"] ="#gamma+jet MPF scale (0.2%)";
  mpar["bm32_mpfscale_02_mpfchs1_zlljet"] = "Z(ll)+jet MPF scale (0.2%)";
  mpar["bm9_multijet_jer_src0"] = "multijet JER";
  mpar["bm18_eesfromzee_gamjet_eig0"] = "EM scale s_{0}";
  mpar["bm18_eesfromzee_gamjet_eig1"] = "s_{1}";//"EM scale s_{1}";
  mpar["bm18_eesfromzee_gamjet_eig2"] = "s_{2}";//"EM scale s_{2}";
  //
  mpar["p0"] = "Tracks (p0)";
  mpar["p1"] = "Photons (p1)";
  mpar["p2"] = "Hadrons (p2)";
  mpar["p3"] = "HCAL (p3)";
  mpar["p4"] = "ECAL (p4)";
  mpar["p5"] = "Herwig (p5)";
  mpar["p6"] = "L1RC (p6)";
  mpar["p7"] = "Trk data (p7)";



  //TH1D *h = new TH1D("h",";(#hat{#theta}-#theta_{0})/#Delta#theta;",
  //		     100,-2, 2);
  TH2D *h2 = new TH2D("h2",";(#hat{#theta}-#theta_{0})/#Delta#theta;",
		      100,-3,+3,hp->GetNbinsX(),-0.5,hp->GetNbinsX()-0.5);

  //TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);
  lumi_13TeV = "UL17";
  TCanvas *c2 = tdrCanvas("c2",(TH1D*)h2,4,0,kSquare);
  c2->SetWindowSize(500,1000);
  c2->SetLeftMargin(0.50);

  TGraphErrors *g = new TGraphErrors(0);
  vector<string> vpar;
  for (int i = 1; i != hp->GetNbinsX()+1; ++i) {
    const char *cpar = hp->GetXaxis()->GetBinLabel(i);
    TString par = cpar;
    //if (!(par.Contains("inactive") || par.Contains("p0") || 
    //	  par.Contains("p1") || par.Contains("p2"))) {
    //if ( !par.Contains("inactive") && par.Contains("bm") ) { // check these are not constrained
    if ( !par.Contains("inactive")
	 && !par.Contains("p8")) { // check these are not constrained
      int j = g->GetN();
      g->SetPoint(j, hp->GetBinContent(i), j);
      g->SetPointError(j, hp->GetBinError(i), 0);
      //h->GetYaxis()->SetBinLabel(j, par);
      h2->GetYaxis()->SetBinLabel(j+1, mpar[cpar]);
      vpar.push_back(cpar);
      cout << Form("  mpar[\"%s\"] = \"\";",cpar) << endl;
    }
  }
  const int np = g->GetN();
  h2->GetYaxis()->SetRangeUser(-0.5,np-0.5);
  //g->Draw("SAMEP");


  // Re-order sources for better readability (numbering from the top)
  map<string,int> midx;
  //const int n = 7; // non-FSR sources
  //const int n = np - 22;
  const int nfsr = 6+6; // FSR sources
  const int nf = 8 - 1; // number of fit parameters -1
  const int n = np - nfsr - 1;  // number of non-FSR sources -1
  midx["bm1_fsru1"] = n+1;
  midx["bm2_hkfsr_ptchs_gamjet_eig0"] = n+3;//2;
  midx["bm2_hkfsr_ptchs_gamjet_eig1"] = n+4;//3;
  midx["bm4_fsru1"] = n+5;//4;
  midx["bm16_fsru1"] = n+7;//5;
  midx["bm32_hkfsr_mpfchs1_gamjet_eig0"] = n+9;//6;
  midx["bm64_fsru1"] = n+11;//7;
  midx["bm1_fsrn1"] = n+2;//8;
  midx["bm4_fsrn1"] = n+6;//9;
  midx["bm16_fsrn1"] = n+8;//10;
  midx["bm32_hkfsr_mpfchs1_gamjet_eig1"] = n+10;//11;
  midx["bm64_fsrn1"] = n+12;//12;
  //
  midx["bm34_scale_050_gamjet"] = nf+2;
  midx["bm68_scale_010_zlljet"] = nf+1;
  midx["bm32_mpfscale_02_mpfchs1_gamjet"] = nf+7;
  midx["bm64_mpfscale_02_mpfchs1_zlljet"] = nf+6;
  midx["bm34_eesfromzee_gamjet_eig0"] = nf+3;
  midx["bm34_eesfromzee_gamjet_eig1"] = nf+4;
  midx["bm34_eesfromzee_gamjet_eig2"] = nf+5;
  //
  midx["bm8_hadw_ptave_fitprob"] = nf+8;
  midx["bm128_hadw_ptboth_fitprob"] = nf+9;
  midx["bm8_hadw_ptave_fitprob2"] = nf+10;
  midx["bm128_hadw_ptboth_fitprob2"] = nf+11;
  /*
  // 
  midx["bm32_hkfsr_mpfchs1_zlljet_eig0"] = n+1;
  midx["bm32_hkfsr_mpfchs1_zlljet_eig1"] = n+2;
  midx["bm16_hkfsr_mpfchs1_gamjet_eig0"] = n+3;//5;
  midx["bm16_hkfsr_mpfchs1_gamjet_eig1"] = n+4;//6;
  midx["bm8_hkfsr_mpfchs1_multijet_eig0"] = n+5;//9;
  midx["bm8_hkfsr_mpfchs1_multijet_eig1"] = n+6;//10;
  midx["bm4_hkfsr_ptchs_zlljet_eig0"] = n+7;//13;
  midx["bm4_hkfsr_ptchs_zlljet_eig1"] = n+8;//14;
  midx["bm2_hkfsr_ptchs_gamjet_eig0"] = n+9;//16;
  midx["bm2_hkfsr_ptchs_gamjet_eig1"] = n+10;//17;
  midx["bm1_hkfsr_ptchs_multijet_eig0"] = n+11;//19;
  midx["bm1_hkfsr_ptchs_multijet_eig1"] = n+12;//20;
  */
  //midx["bm1_hkfsr_ptchs_multijet_eig2"] = n+21;
  //midx["bm2_hkfsr_ptchs_gamjet_eig2"] = n+18;
  //midx["bm4_hkfsr_ptchs_zlljet_eig2"] = n+15;
  //midx["bm8_hkfsr_mpfchs1_multijet_eig2"] = n+11;
  //midx["bm16_hkfsr_mpfchs1_gamjet_eig2"] = n+7;
  //midx["bm32_hkfsr_mpfchs1_zlljet_eig2"] = n+3;
  //midx["bm1_hkfsr_ptchs_multijet_eig3"] = n+22;
  //midx["bm8_hkfsr_mpfchs1_multijet_eig3"] = n+12;
  //midx["bm16_hkfsr_mpfchs1_gamjet_eig3"] = n+8;
  //midx["bm32_hkfsr_mpfchs1_zlljet_eig3"] = n+4;
  midx["bm18_scale_050_gamjet"] = nf+2;
  midx["bm36_scale_020_zlljet"] = nf+1;
  midx["bm16_mpfscale_02_mpfchs1_gamjet"] = nf+7;
  midx["bm32_mpfscale_02_mpfchs1_zlljet"] = nf+6;
  //midx["bm9_multijet_jer_src0"] = 7;
  midx["bm18_eesfromzee_gamjet_eig0"] = nf+3;
  midx["bm18_eesfromzee_gamjet_eig1"] = nf+4;
  midx["bm18_eesfromzee_gamjet_eig2"] = nf+5;
  midx["p0"] = 0;
  midx["p1"] = 1;
  midx["p2"] = 2;
  midx["p3"] = 3;
  midx["p4"] = 4;
  midx["p5"] = 5;
  midx["p6"] = 6;
  midx["p7"] = 7;

  bool reorder = true;
  if (reorder) {

    TGraphErrors *g2 = new TGraphErrors(np);
    for (int i = 0; i != np; ++i) {
      string spar = vpar[i].c_str();
      int j = np-1 - midx[spar];
      if (midx[spar]==0) {
	cout << Form("par j=%d not found? (n=%d,np=%d,midx=%d,%s)",
		     j,n,np,midx[spar],spar.c_str())
	     << endl <<flush;
	assert(j>=0);
	assert(j<np);
	assert(midx[spar]!=0 || spar=="p0");
      }
      if (g2->GetX()[j]!=0 || j<0 || j>=np) {
	cout << Form("par j=%d taken (n=%d,np=%d,midx=%d,%s)",
		     j,n,np,midx[spar],spar.c_str())
	     << endl <<flush;
	assert(j>=0);
	assert(j<np);
	assert(g2->GetX()[j]==0);
      }
      g2->SetPoint(j, g->GetX()[i], j);
      g2->SetPointError(j, g->GetEX()[i], 0);
      h2->GetYaxis()->SetBinLabel(j+1, mpar[spar]);
    }
    g2->Draw("SAMEP");
  }

  TLatex *tex = new TLatex();
  tex->SetTextSize(0.045); tex->SetNDC();
  tex->DrawLatex(0.05,.005,Form("n_{s} = %d",g->GetN()));

  TLine *l = new TLine();
  l->SetNDC();
  l->SetLineStyle(kDashed);
  l->DrawLine(-1,-0.5,-1,g->GetN()-0.5);
  l->DrawLine(+1,-0.5,+1,g->GetN()-0.5);
  l->SetLineStyle(kDotted);
  l->DrawLine(0,-0.5,0,g->GetN()-0.5);

  c2->Update();
  c2->SaveAs("pdf/globalFitPulls_pulls.pdf");

}
