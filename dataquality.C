#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TMath.h"
#include "tdrstyle_mod.C"

// Draw 2D plot of jet rates in (eta,phi) to spot issues
void dataquality() {

  TDirectory *curdir = gDirectory;
  setTDRStyle();
  gStyle->SetOptTitle();
  gStyle->SetPalette(1);

  //TFile *f = new TFile("rootfiles/output_RunFlateG_Feb03V0/nol2l3res/output-DATA-1.root","READ");
  //TFile *fd = new TFile("rootfiles/output_RunFlateG_Feb03V0/output-DATA-1.root","READ");
  //TFile *fd = new TFile("rootfiles/exclusiontests/normalH/output-DATA-1.root","READ");
  //TFile *fd = new TFile("rootfiles/exclusiontests/normalfG/output-DATA-1.root","READ");
  //TFile *fd = new TFile("rootfiles/exclusion_cmsweek/normal/BCD/output-DATA-1.root","READ");
  //TFile *fd = new TFile("rootfiles/exclusion_cmsweek/normal/EF/output-DATA-1.root","READ");
  //TFile *fd = new TFile("rootfiles/exclusion_cmsweek/normal/BCDEF/output-DATA-1.root","READ");
  TFile *fd = new TFile("rootfiles/exclusion_cmsweek/mikko/GH/output-DATA-1.root","READ");
  assert(fd && !fd->IsZombie());

  TFile *fd2 = new TFile("rootfiles/exclusion_cmsweek/mikko/BCDEF/output-DATA-1.root","READ");
  assert(fd2 && !fd2->IsZombie());

  //TFile *fm = new TFile("rootfiles/exclusiontests/normalMC_fG/output-MC-1.root","READ"); 
  //TFile *fm = new TFile("rootfiles/exclusion_cmsweek/mc/BCD/output-MC-1.root","READ"); 
  //TFile *fm = new TFile("rootfiles/exclusion_cmsweek/mc/EF/output-MC-1.root","READ"); 
  //TFile *fm = new TFile("rootfiles/exclusion_cmsweek/mc/BCDEF/output-MC-1.root","READ");
  TFile *fm = new TFile("rootfiles/exclusion_cmsweek/mc/GH/output-MC-1.root","READ"); 
  //TFile *fm = new TFile("rootfiles/output_RunFlateG_Feb03V0/output-MC-1.root","READ"); 
  assert(fm && !fm->IsZombie());

  //TFile *fh = new TFile("rootfiles/hotjets-RunG.root","READ");
  //TFile *fh = new TFile("rootfiles/hotjets-RunBCD.root","READ");
  //TFile *fh = new TFile("rootfiles/hotjets-RunEF.root","READ");
  //TFile *fh = new TFile("rootfiles/hotjets-RunBCDEF.root","READ");
  //TFile *fh = new TFile("rootfiles/hotjets-RunGH.root","READ");
  //TFile *fh = new TFile("rootfiles/hotjets-RunBCDEFGH.root","READ");
  TFile *fh = new TFile("rootfiles/coldjets-RunBCDEFGH.root","READ");
  TH2D *h2jet(0);
  double minsumsig(0);
  if (fh && !fh->IsZombie()) {
    h2jet = (TH2D*)fh->Get("h2jet");
    assert(h2jet);
  }


  const int ntrg = 9;//10;
  string triggers[ntrg] = {"jt40", "jt60", "jt80", "jt140", "jt200",
			   "jt260", "jt320", "jt400", /*"data",*/ "mc"};
  TH2D *h2hots[ntrg], *h2hotr(0), *h2hotm(0);// *h2data(0);
  TH2D *h2colds[ntrg];

  //const char *ctrg = "jt400";
  for (int itrg = 0; itrg != ntrg; ++itrg) {

  //const char *ctrg = "jt400";
  string strg = triggers[itrg];
  const char *ctrg = triggers[itrg].c_str();
  TFile *f = (strg=="mc" ? fm : fd);
  TFile *f2 = (strg=="mc" ? fm : fd2);

  assert(f->cd("Standard"));
  TDirectory *din = gDirectory;

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  gPad->SetLeftMargin(0.10);
  gPad->SetRightMargin(0.15);
  gPad->SetTopMargin(0.10);
  gPad->SetBottomMargin(0.10);

  TH2D *h2 = 0;
  const int neta = 8;
  double etabins[neta+1] = {0,0.5,1.0,1.5,2.0,2.5,3.0,3.2,4.7};
  for (int ieta = 0; ieta != neta; ++ieta) {

    double etamin = etabins[ieta]; double etamax = etabins[ieta+1];
    assert(din->cd(Form("Eta_%1.1f-%1.1f",etamin,etamax)));
    //assert(gDirectory->cd(ctrg));
    //assert(gDirectory->cd(strg=="mc"||strg=="data" ? "jt450" : ctrg));
    assert(gDirectory->cd(strg=="mc"||strg=="data" ? "jt40" : ctrg));
    //assert(gDirectory->cd("jt400"));
    //assert(gDirectory->cd("jt40")); // large correlated regions
    //assert(gDirectory->cd("jt60")); // noisy
    //assert(gDirectory->cd("jt80")); // noisy
    //assert(gDirectory->cd("jt140")); // noisy, maxeta 3.5
    //assert(gDirectory->cd("jt200")); // maxeta 3.2
    //assert(gDirectory->cd("jt260")); // maxeta 3.0
    //assert(gDirectory->cd("jt320")); // maxeta 2.8
    TDirectory *d = gDirectory;
    
    TH2D *h = (TH2D*)d->Get("hetaphi"); assert(h);
    if (!h2) h2 = (TH2D*)h->Clone("h2");
    else h2->Add(h);

    if (strg!="mc" && f2) {
	assert(f2->cd("Standard"));
	assert(gDirectory->cd(Form("Eta_%1.1f-%1.1f",etamin,etamax)));
	assert(gDirectory->cd(strg=="mc"||strg=="data" ? "jt450" : ctrg));
	TH2D *hb = (TH2D*)gDirectory->Get("hetaphi"); assert(hb);
	h2->Add(hb);
    }
  } // for ieta

  // if (!h2data) {
  //   h2data = (TH2D*)h2->Clone("h2data");
  //   for (int i = 1; i != h2->GetNbinsX()+1; ++i) {
  //     for (int j = 1; j != h2->GetNbinsY()+1; ++j) {
  // 	h2data->SetBinContent(i, j, 0);
  //     }
  //   }
  // }
  // if (strg=="data") h2 = h2data;
  // if (strg!="data"&&strg!="mc"&&strg!="jt40"&&h2data) h2data->Add(h2);

  // Create map of known hot ECAL regions from Robert Schoefbeck:
  // https://github.com/schoef/JetMET/blob/master/JEC/python/L2res/jet_cleaning.py#L2-L8
  const int nhot = 7;
  double thr[nhot][nhot] = {{ -2.650, -2.500,  -1.35, -1.05 },
			    { -2.964, -2.650,  -1.10, -0.80 },
			    { -2.964, -2.650,  -0.25, 0.1 },
			    { -2.964, -2.650,  -3.14159, -2.8 },
			    { -2.964, -2.650,  2.9, 3.14159 },
			    { 2.650, 2.964,  -2., -1.6 },
			    { 2.650, 3.139,  0, 0.25 }};
  const int nhot2 = 10;
  double dx = TMath::TwoPi()/72.;
  double thr2[nhot2][4] =
    {// Mikko's extra regions first
      {2.6, 2.6+dx,   -2.7-dx, -2.7+dx},
      {1.4, 1.4+2*dx, -0.8-dx, -0.8+dx},
      {-3.2, -3.2+dx,  0.2, 0.2+2*dx},
      {2.6, 2.6+dx,    0.35, 0.35+3*dx},
      {2.4, 2.4+2*dx, -0.1, -0.1+2*dx },
      // Robert's sub-regions next
      {2.8, 2.8+2*dx, 0, 0+2*dx},
      {2.7, 2.7+dx, -1.8, -1.8+2*dx},
      {-2.8-dx, -2.8+dx, 3.0, 3.0+2*dx},
      {-2.7, -2.7+dx, -1.2, -1.2+2*dx},
      {-2.8, -2.8+dx, -3.15, -3.15+2*dx}
    };
    

  TH2D *h2hot = (TH2D*)h2->Clone("h2hot");
  TH2D *h2cold = (TH2D*)h2->Clone("h2cold");
  h2hotr = (TH2D*)h2->Clone("h2hotr");
  h2hotm = (TH2D*)h2->Clone("h2hotm");
  for (int i = 1; i != h2hot->GetNbinsX()+1; ++i) {
    for (int j = 1; j != h2hot->GetNbinsY()+1; ++j) {

      double eta = h2hot->GetXaxis()->GetBinCenter(i);
      double phi = h2hot->GetYaxis()->GetBinCenter(j);

      // Reset hot region map (fill later from h2a)
      h2hot->SetBinContent(i, j, 0);
      h2hot->SetBinError(i, j, 0);
      h2cold->SetBinContent(i, j, 0);
      h2cold->SetBinError(i, j, 0);

      // Produce map for Robert's hot spots
      h2hotr->SetBinContent(i, j, -10);
      h2hotr->SetBinError(i, j, 0);
      for (int k = 0; k != nhot; ++k) {
	if (eta >= thr[k][0] && eta <= thr[k][1] &&
	    phi >= thr[k][2] && phi <= thr[k][3]) {
	  h2hotr->SetBinContent(i, j, 10);
	} // hot region
      } // for k

      // Produce map for Mikko's manual hot spots
      h2hotm->SetBinContent(i, j, -10);
      h2hotm->SetBinError(i, j, 0);
      for (int k = 0; k != nhot2; ++k) {
	if (eta >= thr2[k][0] && eta <= thr2[k][1] &&
	    phi >= thr2[k][2] && phi <= thr2[k][3]) {
	  h2hotm->SetBinContent(i, j, 10);
	} // hot region 2
      } // for k

    } // for j
  } // for i
  h2hotr->GetZaxis()->SetRange(0,10);
  h2hotr->SetFillStyle(0);
  h2hotr->SetLineColor(kGray);
  h2hotm->GetZaxis()->SetRange(0,10);
  h2hotm->SetFillStyle(0);
  h2hotm->SetLineColor(kRed+1);

  h2->SetTitle("Number of jets;#eta_{jet};#phi_{jet}");
  h2->GetYaxis()->SetRangeUser(-TMath::Pi(),TMath::Pi());
  h2->DrawClone("COLZ");   
  h2hotr->DrawClone("BOX SAME");
  //h2hotm->DrawClone("BOX SAME");

  TH2D *h2a = (TH2D*)h2->Clone("h2a"); // stat significance
  TH2D *h2b = (TH2D*)h2->Clone("h2b"); // relative fluctuation
  TH1D *hrms = new TH1D("hrms","Relative RMS of #eta_{jet} strip;#eta_{jet};"
			"Relative RMS",
			72,-TMath::Pi(),TMath::Pi());
  TH1D *hrms2 = new TH1D("hrms2","Relative RMS of #eta_{jet} strip;#eta_{jet};"
			 "Expected relative RMS",
			 72,-TMath::Pi(),TMath::Pi());

  for (int i = 1; i != h2->GetNbinsX()+1; ++i) {

    // Calculate pedestal
    double sum(0), esum(0);
    int nbins(0), nbinsgood(0);
    for (int j = 1; j != h2->GetNbinsY()+1; ++j) {
      sum += h2->GetBinContent(i,j);
      esum += (h2->GetBinError(i, j) ?
	       pow(h2->GetBinContent(i, j)/h2->GetBinError(i, j),2) : 0.);
      if (h2->GetBinContent(i,j)!=0) ++nbins;
    }
    double ped = sum/nbins;

    // Subtract pedestal, calculate raw rms
    double rms2(0);
    for (int j = 1; j != h2->GetNbinsY()+1; ++j) {
      if (h2->GetBinContent(i,j)!=0) {
	h2a->SetBinContent(i, j, (h2->GetBinContent(i,j)-ped) / sqrt(ped));
	h2b->SetBinContent(i, j, (h2->GetBinContent(i,j)-ped) / ped);
	rms2 += pow(h2b->GetBinContent(i,j), 2);
      }
      int k = hrms->FindBin(h2->GetXaxis()->GetBinCenter(i));
    } // for j
    double rms = sqrt(rms2/nbins);
    hrms->SetBinContent(i, rms);
    hrms2->SetBinContent(i, 1./sqrt(sum/nbins));

    // Calculate cleaned RMS and pedestal with only [-maxsig,+maxsig] range
    double rms2good(0), sumgood(0);
    int nbinsbood(0);
    const int maxsig = 2.;
    for (int j = 1; j != h2->GetNbinsY()+1; ++j) {
      if (h2->GetBinContent(i,j)!=0 &&
	  fabs((h2->GetBinContent(i,j)-ped)/(ped*rms)) < maxsig) {
	sumgood += h2->GetBinContent(i,j);
	rms2good += pow(h2b->GetBinContent(i,j), 2);
	++nbinsgood;
      }
    }
    double pedgood = sumgood/nbinsgood;
    double rmsgood = sqrt(rms2good/nbinsgood);

    // For MC, replace sqrt(ped) with rms, because not using event counts
    //if (triggers[itrg]=="mc") {
      for (int j = 1; j != h2->GetNbinsY()+1; ++j) {
	if (h2->GetBinContent(i,j)!=0) {
	  h2a->SetBinContent(i, j, (h2->GetBinContent(i,j)-pedgood)/
			     (pedgood*rmsgood));
	}
      } // for j
      //} // mc
  } // for i


  TCanvas *c2a = new TCanvas("c2a","c2a",600,600);
  gPad->SetLeftMargin(0.10);
  gPad->SetRightMargin(0.15);
  gPad->SetTopMargin(0.10);
  gPad->SetBottomMargin(0.10);

  // For drawing, mark deficit <-8 as -8 (so blue box instead of white)
  //for (int i = 1; i != h2a->GetNbinsX()+1; ++i) {
  //for (int j = 1; j != h2a->GetNbinsY()+1; ++j) {
  //  if (h2a->GetBinContent(i,j)<-8)
  //h2a->SetBinContent(i,j,-8);
  //}
  //}
  h2a->SetTitle(Form("Significance of excess/deficit (%s);"
		     "#eta_{jet};#phi_{jet}",ctrg));

  if (strg=="mc" || f2) {
  //   h2a->Scale(8./10.);
  //   h2a->SetTitle(Form("Significance of excess/deficit #times 8/10 (%s);"
  // 		       "#eta_{jet};#phi_{jet}",ctrg));
    for (int i = 1; i != h2a->GetNbinsX()+1; ++i) {
      for (int j = 1; j != h2a->GetNbinsY()+1; ++j) {
	if (h2a->GetBinContent(i, j) < -8)
	  h2a->SetBinContent(i, j, -8);
      }
    }
  }
  // if (strg=="mc") {
  //   h2a->Scale(8./18.);
  //   h2a->SetTitle("Significance of excess/deficit #times 8/18 (MC);"
  // 		  "#eta_{jet};#phi_{jet}");
  // }
  // if (strg=="data") {
  //   //h2a->Scale(0.1);
  //   h2a->SetTitle("Significance of excess/deficit #times 1 (data);"
  // 		  "#eta_{jet};#phi_{jet}");
  // }

  //h2a->SetMinimum(strg=="mc" ? -10 : -8);//-4);
  //h2a->SetMaximum(strg=="mc" ? +10 : +8);//+4);
  h2a->SetMinimum(-8);//-4);
  h2a->SetMaximum(+8);//+4);

  h2a->DrawClone("COLZ");
  h2hotr->Draw("SAMEBOX");
  //h2hotm->Draw("SAMEBOX");

  if (h2jet) {
    h2jet->GetZaxis()->SetRange(0,40);
    h2jet->SetFillStyle(0);
    h2jet->SetLineColor(kBlack);
    h2jet->DrawClone("BOXSAME");
  }

  // Keep track of found hot regions
  h2hots[itrg] = h2hot;
  h2colds[itrg] = h2cold;
  const double minsig = 3.;
  const double maxsig = 8.5;
  minsumsig = 9.;
  for (int i = 1; i != h2hot->GetNbinsX()+1; ++i) {
    for (int j = 1; j != h2hot->GetNbinsY()+1; ++j) {  
      if (h2a->GetBinContent(i,j)>=minsig) {
	h2hot->SetBinContent(i, j, min(maxsig, h2a->GetBinContent(i,j)));
      }
      if (h2a->GetBinContent(i,j)<=-minsig) {
	h2cold->SetBinContent(i, j, min(maxsig, -h2a->GetBinContent(i,j)));
      }
    }
  }
  //if (strg!="data"&&strg!="mc") h2data->Add(h2a);


  TCanvas *c2b = new TCanvas("c2b","c2b",600,600);
  gPad->SetLeftMargin(0.10);
  gPad->SetRightMargin(0.15);
  gPad->SetTopMargin(0.10);
  gPad->SetBottomMargin(0.10);

  h2b->SetTitle("Relative fluctuation;#eta_{jet};#phi_{jet}");
  h2b->SetMinimum(-0.25);
  h2b->SetMaximum(+0.25);
  h2b->DrawClone("COLZ");

  TCanvas *c3 = new TCanvas("c3","c3",600,600);
  gPad->SetTopMargin(0.10);
  hrms->SetMinimum(0.);
  hrms->SetMaximum(0.5);
  hrms->SetLineColor(kRed);
  hrms->Draw();
  hrms2->SetLineColor(kBlue);
  hrms2->Draw("SAME");

  if (itrg==ntrg-1) {
    c1->SaveAs("pdf/dataquality_njet.pdf");
    c2a->SaveAs("pdf/dataquality_significance.pdf");
    c2b->SaveAs("pdf/dataquality_relfluctuation.pdf");
    c3->SaveAs("pdf/dataquality_relRMS.pdf");
  }

  c2a->SaveAs(Form("pdf/dataquality_significance_%s.pdf",ctrg));
  } // itrg

  // Sum up hot region maps
  TH2D *h2hot = (TH2D*)h2hots[0]->Clone("h2hot");
  TH2D *h2hot2 = (TH2D*)h2hots[0]->Clone("h2hot2");
  TH2D *h2cold = (TH2D*)h2colds[0]->Clone("h2cold");
  TH2D *h2cold2 = (TH2D*)h2colds[0]->Clone("h2cold2");
  for (int i = 1; i != h2hot->GetNbinsX()+1; ++i) {
    for (int j = 1; j != h2hot->GetNbinsY()+1; ++j) {

      // Hot regions
      h2hot->SetBinContent(i, j, 0);
      h2hot2->SetBinContent(i, j, 0);
      for (int itrg = 0; itrg != ntrg; ++itrg) {
	if (triggers[itrg]!="jt40" && triggers[itrg]!="mc") {
	  h2hot->SetBinContent(i, j, h2hot->GetBinContent(i,j)
			       + h2hots[itrg]->GetBinContent(i,j));
	}
      } // for itrg

      // Select regions with excesses in at least two triggers
      if (h2hot->GetBinContent(i,j)>=minsumsig)
	h2hot2->SetBinContent(i,j,10);
      else
	h2hot2->SetBinContent(i,j,-10);

      // Cold regions
      h2cold->SetBinContent(i, j, 0);
      h2cold2->SetBinContent(i, j, 0);
      for (int itrg = 0; itrg != ntrg; ++itrg) {
	if (triggers[itrg]!="jt40" && triggers[itrg]!="mc") {
	  h2cold->SetBinContent(i, j, h2cold->GetBinContent(i,j)
				+ h2colds[itrg]->GetBinContent(i,j));
	}
      } // for itrg

      // Select regions with deficits in at least two triggers
      if (h2cold->GetBinContent(i,j)>=minsumsig)
	h2cold2->SetBinContent(i,j,10);
      else
	h2cold2->SetBinContent(i,j,-10);
      
    } // for j
  } // for i

  //TFile *fout = new TFile("rootfiles/hotjets-runG.root","RECREATE");
  //TFile *fout = new TFile("rootfiles/hotjets-runBCD.root","RECREATE");
  //TFile *fout = new TFile("rootfiles/hotjets-runEF.root","RECREATE");
  //TFile *fout = new TFile("rootfiles/hotjets-runBCDEF.root","RECREATE");
  //TFile *fout = new TFile("rootfiles/hotjets-runGH.root","RECREATE");
  TFile *fout = new TFile("rootfiles/coldjets-runBCDEFGH.root","RECREATE");
  h2hot->Write("h2hot");
  h2hot2->Write("h2jet");
  //h2hotr->Write("h2hotr");
  //h2hotm->Write("h2hotm");
  h2cold->Write("h2cold");
  h2cold2->Write("h2hole");
  fout->Close();
}
