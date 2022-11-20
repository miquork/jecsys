// Purpose: plot the pulls and post-fit parameter values
//          produced by globalFitL3Res.C and
//          stored in rootfiles/globalFitL3Res_emat.root
//          in a nice Higgs Combine tool style
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH2D.h"
#include "TLine.h"
#include "TArrow.h"

#include "tdrstyle_mod15.C"

#include <string>

void globalFitPull(string set);

void globalFitPulls() {
  /*
  globalFitPull("2016BCD");
  globalFitPull("2016EF");
  globalFitPull("2016BCDEF");
  globalFitPull("2016GH");
  */
  globalFitPull("Run2Test");
}

//void globalFitPulls(string set="2016BCDEF") {
//void globalFitPulls(string set="2016GH") {
//void globalFitPulls(string set="2016BCD") {
  //void globalFitPulls(string set="2016EF") {
void globalFitPull(string set) { // singular

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  const char *cs = set.c_str();

  TFile *f = new TFile(Form("rootfiles/globalFitL3Res_emat_%s.root",cs),"READ");
  assert(f && !f->IsZombie());

  curdir->cd();
  
  TH1D *hp = (TH1D*)f->Get("h1par"); assert(hp);
  TH1D *h2cov = (TH1D*)f->Get("h2cov"); assert(hp);

  // Map names to better ones for graph
  map<string,const char*> mpar;
  mpar["p0"] = "Tracks (p0)";
  mpar["p1"] = "Photons (p1)";
  mpar["p2"] = "Hadrons (p2)";
  mpar["p3"] = "Had. HCAL (p3)";
  mpar["p4"] = "Had. ECAL (p4)";
  mpar["p5"] = "Shower (p5)";
  mpar["p6"] = "Offset (p6)";
  mpar["p7"] = "MC tracks (p7)";
  mpar["p8"] = "Flavor (p8)";
  //
  //mpar["bm68_scale_010_zlljet"] = "Z(ll) scale (0.10%)";
  mpar["bm68_scale_020_zjet"] = "Z scale (0.20%)";
  //mpar["bm34_scale_020_gamjet"] = "#gamma scale (0.20%)";
  mpar["bm34_scale_050_gamjet"] = "#gamma scale (0.50%)";
  mpar["bm34_eesfromzee_gamjet_eig0"] = "EM scale s_{0}";
  mpar["bm34_eesfromzee_gamjet_eig1"] = "s_{1}";//"EM scale s_{1}";
  mpar["bm34_eesfromzee_gamjet_eig2"] = "s_{2}";//"EM scale s_{2}";
  mpar["bm34_hdmscale_020_mpfchs1_gamjet"] ="#gamma+jet HDM scale (0.20%)";
  //mpar["bm68_hdmscale_020_mpfchs1_zlljet"] = "Z(ll)+jet HDM scale (0.20%)";
  mpar["bm68_hdmscale_020_mpfchs1_zjet"] = "Z+jet HDM scale (0.20%)";
  //mpar["bm9_multijet_jer_src0"] = "multijet JER";
  //
  mpar["bm8_hadw_ptave_fitprob"] = "W>qq' fitProb (pTave) s_{0}";
  mpar["bm8_hadw_ptave_fitprob2"] = "s_{1}";//W>qq' fitProb (pTave) s_{1}";
  mpar["bm128_hadw_ptboth_fitprob"] = "W>qq' fitProb (pTboth) s_{0}";
  mpar["bm128_hadw_ptboth_fitprob2"] = "s_{1}";//W>qq' fitProb (pTboth) s_{1}";
  //
  //mpar["bm4_zlljet_fsru1"] = "Z(ll)+jet FSR uncl. (DB)";
  //mpar["bm4_zlljet_fsrn1"] = "Z(ll)+jet FSR jets (DB)";
  mpar["bm4_zjet_fsru1"] = "Z+jet FSR uncl. (DB)";
  mpar["bm4_zjet_fsrn1"] = "Z+jet FSR jets (DB)";
  mpar["bm1_multijet_fsru1"] = "Multijet FSR uncl. (DB)";
  mpar["bm1_multijet_fsrn1"] = "Multijet FSR jets (DB)";
  mpar["bm2_gamjet_fsru1"] = "#gamma+jet FSR uncl. (DB)";
  mpar["bm2_gamjet_fsrn1"] = "#gamma+jet FSR jets (DB)";
  //mpar["bm2_gamjet_hkfsr_ptchs_gamjet_eig0"] = "#gamma+jet FSR (DB) s_{0}";
  //mpar["bm2_gamjet_hkfsr_ptchs_gamjet_eig1"] = "s_{1}";//#gamma+jet FSR (DB) s_{0}";
  //
  //mpar["bm32_gamjet_hkfsr_mpfchs1_gamjet_eig0"] = "#gamma+jet FSR (MPF) s_{0}";
  //mpar["bm32_gamjet_hkfsr_mpfchs1_gamjet_eig1"] = "s_{1}";
  mpar["bm32_gamjet_fsru1"] = "#gamma+jet FSR uncl. (MPF)";
  mpar["bm32_gamjet_fsrn1"] = "#gamma+jet HDM jets (MPF)";
  mpar["bm16_multijet_fsru1"] = "Multijet FSR uncl. (MPF)";
  mpar["bm16_multijet_fsrn1"] = "Multijet FSR jets (MPF)";
  //mpar["bm64_zlljet_fsru1"] = "Z(ll)+jet FSR uncl. (MPF)";
  //mpar["bm64_zlljet_fsrn1"] = "Z(ll)+jet HDM jets (MPF)";
  mpar["bm64_zjet_fsru1"] = "Z+jet FSR uncl. (MPF)";
  mpar["bm64_zjet_fsrn1"] = "Z+jet HDM jets (MPF)";


  TH2D *h2 = new TH2D("h2",";(#hat{#theta}-#theta_{0})/#Delta#theta;",
		      100,-3,+3,hp->GetNbinsX(),-0.5,hp->GetNbinsX()-0.5);

  lumi_13TeV = set.c_str();
  if (set=="2016BCDEF") lumi_13TeV = "2016BCDEF";
  if (set=="2016GH")    lumi_13TeV = "2016GH";
  if (set=="Run2Test")  lumi_13TeV = "Run 2, 136 fb^{-1}";

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
    //if ( !par.Contains("inactive")
    //	 && !par.Contains("p8")) { // check these are not constrained
    if ( !par.Contains("inactive") ) { // check these are not constrained
      int j = g->GetN();
      g->SetPoint(j, hp->GetBinContent(i), j);
      g->SetPointError(j, hp->GetBinError(i), 0);
      //h->GetYaxis()->SetBinLabel(j, par);
      h2->GetYaxis()->SetBinLabel(j+1, mpar[cpar]);
      vpar.push_back(cpar);
      cout << Form("  mpar[\"%s\"] = \"%s\";",cpar,mpar[cpar]) << endl;
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
  const int nf = 9 - 1; // number of fit parameters -1
  const int n = np - nfsr - 1;  // number of non-FSR sources -1
  midx["p0"] = 0;
  midx["p1"] = 1;
  midx["p2"] = 2;
  midx["p3"] = 3;
  midx["p4"] = 4;
  midx["p5"] = 5;
  midx["p6"] = 6;
  midx["p7"] = 7;
  midx["p8"] = 8;
  //
  //midx["bm68_scale_010_zlljet"] = nf+1;
  midx["bm68_scale_020_zjet"] = nf+1;
  //midx["bm34_scale_020_gamjet"] = nf+2;
  midx["bm34_scale_050_gamjet"] = nf+2;
  midx["bm34_eesfromzee_gamjet_eig0"] = nf+3;
  midx["bm34_eesfromzee_gamjet_eig1"] = nf+4;
  midx["bm34_eesfromzee_gamjet_eig2"] = nf+5;
  midx["bm34_hdmscale_020_mpfchs1_gamjet"] = nf+7;
  //midx["bm68_hdmscale_020_mpfchs1_zlljet"] = nf+6;
  midx["bm68_hdmscale_020_mpfchs1_zjet"] = nf+6;
  //
  midx["bm8_hadw_ptave_fitprob"] = nf+8;
  midx["bm8_hadw_ptave_fitprob2"] = nf+9;
  midx["bm128_hadw_ptboth_fitprob"] = nf+10;
  midx["bm128_hadw_ptboth_fitprob2"] = nf+11;
  //
  //midx["bm4_zlljet_fsru1"] = n+1;
  //midx["bm4_zlljet_fsrn1"] = n+2;
  midx["bm4_zjet_fsru1"] = n+1;
  midx["bm4_zjet_fsrn1"] = n+2;
  midx["bm1_multijet_fsru1"] = n+3;
  midx["bm1_multijet_fsrn1"] = n+4;
  midx["bm2_gamjet_fsru1"] = n+5;
  midx["bm2_gamjet_fsrn1"] = n+6;
  //midx["bm2_gamjet_hkfsr_ptchs_gamjet_eig0"] = n+5;
  //midx["bm2_gamjet_hkfsr_ptchs_gamjet_eig1"] = n+6;
  //
  //midx["bm64_zlljet_fsru1"] = n+7;
  //midx["bm64_zlljet_fsrn1"] = n+8;
  midx["bm64_zjet_fsru1"] = n+7;
  midx["bm64_zjet_fsrn1"] = n+8;
  midx["bm16_multijet_fsru1"] = n+9;
  midx["bm16_multijet_fsrn1"] = n+10;
  midx["bm32_gamjet_fsru1"] = n+11;
  midx["bm32_gamjet_fsrn1"] = n+12;
  //midx["bm32_gamjet_hkfsr_mpfchs1_gamjet_eig0"] = n+11;
  //midx["bm32_gamjet_hkfsr_mpfchs1_gamjet_eig1"] = n+12;

  bool reorder = true;
  if (reorder) {

    TGraphErrors *g2 = new TGraphErrors(np);
    for (int i = 0; i != np; ++i) {
      string spar = vpar[i].c_str();
      int j = np-1 - midx[spar];
      if (midx[spar]==0 && spar!="p0") {
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
      if (mpar[spar]==0) {
	cout << "Label missing for " << spar << endl << flush;
      }
      if (fabs(g->GetX()[i])>3.0) {

	double x = g->GetX()[i];	
	TArrow *arr = new TArrow(TMath::Sign(2.2,x),j,TMath::Sign(3,x),j,0.02);
	arr->SetLineWidth(2);
	arr->SetLineColor(kRed);
	arr->Draw();
	TLatex *tex = new TLatex();
	tex->SetTextSize(0.035); tex->SetTextColor(kRed);
	tex->DrawLatex(x<0 ? -1.8: 0.1, j-0.25,
		       Form("%+1.1f#pm%1.1f",x,g->GetEX()[i]));
      }
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
  //c2->SaveAs("pdf/globalFitPulls_pulls.pdf");
  c2->SaveAs(Form("pdf/globalFitPulls_pulls_%s.pdf",cs));

}
