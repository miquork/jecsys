// Purpose: Estimate bias on <rho> in Z+jet
//          Combine with bias from L2C ("NoPU") vs L2S ("WithPU")
//          to estimate bias or uncertainty on L3Res global fit
#include "TFile.h"
#include "TGraphErrors.h"
#include "TF1.h"

#include "../tdrstyle_mod15.C"

#include <vector>
using namespace std;

TGraphErrors *extrapolate(const int na, const int *va,
			  vector<TGraphErrors*> &vg) {

  TF1 *f1 = new TF1("f1","[0]+[1]*x",0,100);
  assert(vg[0]);
  TGraphErrors *g = (TGraphErrors*)vg[0]->Clone("g0");
  // Loop over pT bins
  for (int i = 0; i != g->GetN(); ++i) {
    if (g->GetX()[i]<30.) continue;
    if (g->GetX()[i]>1000.) continue;
    // Copy values vs alpha
    TGraphErrors *ga = new TGraphErrors(0);
    for (int ia = 0; ia != na; ++ia) {
      // Find matching point in pT
      for (int j = 0; j != vg[ia]->GetN(); ++j) {
	if (fabs(g->GetX()[i]/vg[ia]->GetX()[j]-1)<0.1) {
	  int n = ga->GetN();
	  ga->SetPoint(n, va[ia], vg[ia]->GetY()[j]);
	  ga->SetPointError(n, 0., vg[ia]->GetEY()[j]);
	}
      } // for j
    } // for ia
    ga->Fit(f1,"QRN");
    g->SetPoint(i, g->GetX()[i], f1->GetParameter(0));
    g->SetPointError(i, g->GetEX()[i], f1->GetParError(0));
    //delete ga;
  } // for i

  delete f1;
  return g;
} // extrapolate

void rhoBias() {

  setTDRStyle();

  TDirectory *curdir = gDirectory;
  
  TFile *fmm = new TFile("rootfiles/ZJetCombination_Zmm_DYJets_Madgraph_09Aug2019_Summer19UL17_JECV1_SimpleL1.root","READ");
  assert(fmm && !fmm->IsZombie());
  
  TFile *fjec = new TFile("rootfiles/jecdataBCDEF.root","READ");
  assert(fjec && !fjec->IsZombie());

  curdir->cd();
  
  TH1D *h = new TH1D("h",";p_{T,Z} (GeV);#LT#rho(Z+jet)#GT",3485,15,3500);
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  h->SetMinimum(-10.);
  h->SetMaximum(30.);
  
  lumi_13TeV = "UL17 BCDEF";
  TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);
  gPad->SetLogx();
  
  const int va[] = {30,20,15,10};
  const int na = sizeof(va)/sizeof(va[0]);
  const int color[na] = {kGreen+2, kYellow+2, kOrange+2, kRed+1};

  TLegend *legm = tdrLeg(0.50,0.20,0.75,0.50);
  legm->SetHeader("MC");
  TLegend *legd = tdrLeg(0.57,0.20,0.82,0.50);
  legd->SetHeader("Data");
  
  const char *cd = "Run2017BCDEF/";
  vector<TGraphErrors*> vgd(na);
  vector<TGraphErrors*> vgm(na);
  for (int i = 0; i != na; ++i) {
    const int a = va[i];
    TGraphErrors *gd = (TGraphErrors*)fmm->Get(Form("%sData_Rho_CHS_a%d_eta_00_13_L1L2Res",cd,a));
    assert(gd);
    TGraphErrors *gm = (TGraphErrors*)fmm->Get(Form("%sMC_Rho_CHS_a%d_eta_00_13_L1L2Res",cd,a));
    assert(gm);
    
    tdrDraw(gm,"Pz",kOpenCircle,color[i]); gm->SetMarkerSize(1.0-0.1*i);
    tdrDraw(gd,"Pz",kFullCircle,color[i]); gd->SetMarkerSize(1.0-0.1*i);
    
    legm->AddEntry(gm," ","PL");
    legd->AddEntry(gd,Form("#alpha < %1.2f",0.01*a),"PL");

    vgd[i] = gd;
    vgm[i] = gm;
  } // for i

  TGraphErrors *gd0 = extrapolate(na,va,vgd);
  tdrDraw(gd0,"Pz",kFullStar,kBlack);
  TGraphErrors *gm0 = extrapolate(na,va,vgm);
  tdrDraw(gm0,"Pz",kOpenStar,kBlack);

  legm->AddEntry(gm0," ","PL");
  legd->AddEntry(gd0,"#alpha_{max} #rightarrow 0","PL");

  //TF1 *f0 = new TF1("f0","[0]*min((1+[1]*log(x/[2]),1.)",30,1000);
  //f0->SetParameters(21.5,1,200);
  //TF1 *f0 = new TF1("f0","[0]*(log(x/[2])/[1])/(1+(log(x/[2])/[1]))",30,1000);
  //f0->SetParameters(21.5,log(100),1);
  //TF1 *f0 = new TF1("f0","[0]*(1+[1]*log(x)*pow(1+x/[2],[3]))",30,1000);
  //f0->SetParameters(21.06,-0.93,218.,-6.05);
  //f0->FixParameter(1,-1);
  TF1 *f0 = new TF1("f0","[0]*(1-log(x)*pow(1+x/[1],[2]))",30,1000);
  f0->SetParameters(21.13,158.4,-4.812);
  gd0->Fit(f0,"RNW");
  f0->Draw("SAME");

  TF1 *f0ref = new TF1("f0ref","[0]*(1-log(x)*pow(1+x/[1],[2]))",30,1000);
  f0ref->SetParameters(21.13,158.4,-4.812);
  f0ref->SetLineColor(kGreen+2);
  vgd[0]->Fit(f0ref,"RNW");
  f0ref->Draw("SAME");

  c1->SaveAs("pdf/rhoBias_rho.pdf");


  TH1D *h2 = new TH1D("h2",";p_{T,Z} (GeV);L1 bias",3485,15,3500);
  h2->GetXaxis()->SetMoreLogLabels();
  h2->GetXaxis()->SetNoExponent();
  h2->SetMinimum(0.98);
  h2->SetMaximum(1.08);


  TCanvas *c2 = tdrCanvas("c2",h2,4,11,kSquare);
  c2->SetLogx();

  TH1D *hl2cos = (TH1D*)fjec->Get("ratio/eta00-13/hl2cos");
  assert(hl2cos);
  
  double rhoref = 21.13;
  TH1D *hl2cos0 = (TH1D*)hl2cos->Clone("hl2cos0");
  for (int i = 1; i != hl2cos0->GetNbinsX()+1; ++i) {
    int pt = hl2cos0->GetBinCenter(i);
    double dpt(1000); int k(0);
    for (int j = 0; j != gd0->GetN(); ++j) {
      if (fabs(gd0->GetX()[j]-pt)<dpt) {
	dpt = fabs(gd0->GetX()[j]-pt);
      	k = j;
      }
    } // for j
    double w = 1 - gd0->GetY()[k] / rhoref;
    hl2cos0->SetBinContent(i, 1 + w*(hl2cos->GetBinContent(i)-1));
			   
  } // for i

  tdrDraw(hl2cos,"HIST",kNone,kBlack,kSolid,-1,kNone,0);
  tdrDraw(hl2cos0,"HP",kFullCircle,kBlack,kSolid,-1,kNone,0);

  TF1 *fl1 = new TF1("fl1","1-([0]+[1]*log(x)+[2]*pow(log(x),2))"
  			       "*max(1-[3]*log(x)/log(21.5),0.)/x",10,3500);
  fl1->SetParameters(1.57597,-1.89194,0.383362,0); // UL17 hl2cos
  fl1->SetLineColor(kGreen+2);
  fl1->Draw("SAME");

  TF1 *fl2 = new TF1("fl1","1-([0]+[1]*log(x)+[2]*pow(log(x),2))/x"
		     "*log(x)*pow(1+x/[3],[4])",10,3500);
  fl2->SetParameters(1.57597,-1.89194,0.383362, 158.4,-4.812);
  fl2->SetLineColor(kRed);
  fl2->DrawClone("SAME");
  fl2->SetLineColor(kBlue);
  hl2cos0->Fit(fl2,"WRN");
  fl2->Draw("SAME");

  c2->RedrawAxis();

  c2->SaveAs("pdf/rhoBias_cosxrho.pdf");
} // rhoBias
