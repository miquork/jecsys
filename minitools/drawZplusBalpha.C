// Purpose: Draw Z+b (and Z+j) vs alpha to study difference in FSR
//          between light flavor jets and b quark jets
//          Does the "dead cone" effect increase b quark jet FSR?
#include "TFile.h"
#include "TGraphErrors.h"
#include "TF1.h"

#include "../tdrstyle_mod15.C"

#include <map>
#include <vector>

using namespace std;

// Switches
//string sm = "mpfchs"; // method: mpfchs, ptchs, ratio
//string sm = "ptchs";
string sb = "none"; // b-tag: none, loose, medium, tight
//string sb = "tight";
//string sb = "medium";
//string sd = "data";
string sd = "mc";
//string sd = "ratio";
//double ptref = 170;
//double ptref = 100;
//double ptref = 70;
double ptref = 50;

void drawZplusBalpha() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TFile *f0 = new TFile("rootfiles/jme_bplusZ_Run2_14.root","READ");
  assert(f0 && !f0->IsZombie());
  f0->cd(sd.c_str());
  gDirectory->cd("eta13");
  TFile *f = (TFile*)gDirectory;
  curdir->cd();

  int aa[] = {10, 15, 20, 30, 40, 50, 60, 80, 100};
  const int na = sizeof(aa)/sizeof(aa[0]);

  // Get reference pT binning and find the reference bin
  TH1D *href = (TH1D*)f0->Get("data/eta13/statistics_zmmjet_a100");
  assert(href);
  const int iptref = href->FindBin(ptref);

  // Create canvas for drawing results
  TH1D *h = new TH1D("h",Form(";#alpha;p_{T} balance (%s)",sd.c_str()),
		     11,0,1.1);
  h->SetMinimum(sd=="ratio" ? 0.990 : 0.85);
  h->SetMaximum(sd=="ratio" ? 1.015 : 1.05);

  TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);


  for (int ipt = 1; ipt != href->GetNbinsX()+1; ++ipt) {

    double ptmin = href->GetBinLowEdge(ipt);
    if (ptmin<30 || ptmin>300) continue;

    // Load MPF and pT balance results, map vs alpha
    // (turn TGraphErrors into TH1D for easier handling)
    vector<TH1D*> vm;
    vector<TH1D*> vp;
    TGraphErrors *gam = new TGraphErrors(na);
    TGraphErrors *gap = new TGraphErrors(na);
    TGraphErrors *gar = new TGraphErrors(na);
    for (int ia = 0; ia != na; ++ia) {
      
      const int a = aa[ia];
      TGraphErrors *gm = (TGraphErrors*)f->Get(Form("mpfchs_zmmjet_a%d",a));
      assert(gm);
      TGraphErrors *gp = (TGraphErrors*)f->Get(Form("ptchs_zmmjet_a%d",a));
      assert(gp);
      
      // Copy over contents
      TH1D *hm = (TH1D*)href->Clone(Form("h_%s",gm->GetName()));
      hm->Reset();
      for (int i = 0; i != gm->GetN(); ++i) {
	int j = hm->FindBin(gm->GetX()[i]);
	hm->SetBinContent(j, gm->GetY()[i]);
	hm->SetBinError(j, gm->GetEY()[i]);
      } // for i
      //
      TH1D *hp = (TH1D*)href->Clone(Form("h_%s",gp->GetName()));
      hp->Reset();
      for (int i = 0; i != gp->GetN(); ++i) {
	int j = hp->FindBin(gp->GetX()[i]);
	hp->SetBinContent(j, gp->GetY()[i]);
	hp->SetBinError(j, gp->GetEY()[i]);
      } // for i
      
      if (hm->GetBinLowEdge(ipt)*0.01*a>=15.) {

	double ym = hm->GetBinContent(ipt);
	double eym = hm->GetBinError(ipt);
	//gam->SetPoint(ia, 0.01*a, ym);
	//gam->SetPointError(ia, 0, eym);
	
	double yp = hp->GetBinContent(ipt);
	double eyp = hp->GetBinError(ipt);
	gap->SetPoint(ia, 0.01*a, yp);
	gap->SetPointError(ia, 0, eyp);
	
	gar->SetPoint(ia, 0.01*a, yp/ym);
	gar->SetPointError(ia, 0, yp/ym*sqrt(pow(eyp/yp,2) + pow(eym/ym,2)));
      }
    } // for ia

    if (ipt!=iptref) {
      gam->SetMarkerColor(kBlue-9);//kGray);
      gap->SetMarkerColor(kBlue-9);//kGray);
      gar->SetMarkerColor(kBlue-9);//kGray);

      gam->SetMarkerSize(0.95-0.05*ipt);
      gap->SetMarkerSize(0.95-0.05*ipt);
      gar->SetMarkerSize(0.95-0.05*ipt);
    }

    gam->SetMarkerStyle(kOpenSquare);//kOpenCircle);
    gam->Draw("SAMEP");
    gap->SetMarkerStyle(kOpenCircle);//kFullCircle);
    gap->Draw("SAMEP");
    gar->SetMarkerStyle(kFullCircle);
    //gar->Draw("SAMEP");
  } // for ipt
    
}
