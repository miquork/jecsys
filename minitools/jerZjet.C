// Purpose: use Z+jet RMS to estimate JER SF vs pT
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TH2D.h"

#include "../tdrstyle_mod15.C"

void jerZjet() {
  
  setTDRStyle();
  TDirectory *curdir = gDirectory;
  


  //TFile *f = new TFile("rootfiles/jme_bplusZ_merged_v40_muon2018.root");
  TFile *f = new TFile("rootfiles/jme_bplusZ_merged_v41_muon2018.root");
  assert(f && !f->IsZombie());
  
  f->cd("data");
  gDirectory->cd("eta_00_13");
  TDirectory *dd = gDirectory;

  f->cd("mc");
  gDirectory->cd("eta_00_13");
  TDirectory *dm = gDirectory;

  curdir->cd();
  
  // Missing j2pt_pt15 variants for MPFx, RMPFjet1 and RMPFjet1x in v40
  /*
  TH2D *h2dr = (TH2D*)dd->Get("h_Zpt_RMPF_alpha100_eta_00"); assert(h2dr);
  TH2D *h2dx = (TH2D*)dd->Get("h_Zpt_RMPFx_alpha100_eta_00"); assert(h2dx);
  TH2D *h2dm = (TH2D*)dd->Get("h_Zpt_mZ_alpha100_eta_00"); assert(h2dm);
  //
  TH2D *h2mr = (TH2D*)dm->Get("h_Zpt_RMPF_alpha100_eta_00"); assert(h2mr);
  TH2D *h2mx = (TH2D*)dm->Get("h_Zpt_RMPFx_alpha100_eta_00"); assert(h2mx);
  TH2D *h2mm = (TH2D*)dm->Get("h_Zpt_mZ_alpha100_eta_00"); assert(h2mm);
  */
  TH2D *h2dr = (TH2D*)dd->Get("h_Zpt_RMPF_alpha30_eta_00"); assert(h2dr);
  TH2D *h2dx = (TH2D*)dd->Get("h_Zpt_RMPFx_alpha30_eta_00"); assert(h2dx);
  TH2D *h2dm = (TH2D*)dd->Get("h_Zpt_mZ_alpha30_eta_00"); assert(h2dm);
  //
  TH2D *h2mr = (TH2D*)dm->Get("h_Zpt_RMPF_alpha30_eta_00"); assert(h2mr);
  TH2D *h2mx = (TH2D*)dm->Get("h_Zpt_RMPFx_alpha30_eta_00"); assert(h2mx);
  TH2D *h2mm = (TH2D*)dm->Get("h_Zpt_mZ_alpha30_eta_00"); assert(h2mm);

  // j2pt_pt15 variants added in v41 (poor statistics)
  /*
  TH2D *h2dr = (TH2D*)dd->Get("h_Zpt_RMPFj2pt_pt15_eta_00"); assert(h2dr);
  TH2D *h2dx = (TH2D*)dd->Get("h_Zpt_RMPFj2ptx_pt15_eta_00"); assert(h2dr);
  TH2D *h2dm = (TH2D*)dd->Get("h_Zpt_mZ_alpha100_eta_00"); assert(h2dm);
  //
  TH2D *h2mr = (TH2D*)dm->Get("h_Zpt_RMPFj2pt_pt15_eta_00"); assert(h2mr);
  TH2D *h2mx = (TH2D*)dm->Get("h_Zpt_RMPFj2ptx_pt15_eta_00"); assert(h2mr);
  TH2D *h2mm = (TH2D*)dm->Get("h_Zpt_mZ_alpha100_eta_00"); assert(h2mm);
  */
  /*
  // j2pt_pt30 variants added in v41 (poor statistics)
  TH2D *h2dr = (TH2D*)dd->Get("h_Zpt_RMPFj2pt_pt30_eta_00"); assert(h2dr);
  TH2D *h2dx = (TH2D*)dd->Get("h_Zpt_RMPFj2ptx_pt30_eta_00"); assert(h2dr);
  TH2D *h2dm = (TH2D*)dd->Get("h_Zpt_mZ_alpha100_eta_00"); assert(h2dm);
  //
  TH2D *h2mr = (TH2D*)dm->Get("h_Zpt_RMPFj2pt_pt30_eta_00"); assert(h2mr);
  TH2D *h2mx = (TH2D*)dm->Get("h_Zpt_RMPFj2ptx_pt30_eta_00"); assert(h2mr);
  TH2D *h2mm = (TH2D*)dm->Get("h_Zpt_mZ_alpha100_eta_00"); assert(h2mm);
  */

  TH1D *h1dr = h2dr->ProjectionX("h1dr"); h1dr->Reset();
  TH1D *h1dx = h2dx->ProjectionX("h1dx"); h1dx->Reset();
  TH1D *h1dm = h2dm->ProjectionX("h1dm"); h1dm->Reset();
  TH1D *h1d = h2dr->ProjectionX("h1d"); h1d->Reset();
  assert(h1dx->GetNbinsX()==h1dr->GetNbinsX());
  assert(h1dm->GetNbinsX()==h1dr->GetNbinsX());
  TH1D *h1mr = h2mr->ProjectionX("h1mr"); h1mr->Reset();
  TH1D *h1mx = h2mx->ProjectionX("h1mx"); h1mx->Reset();
  TH1D *h1mm = h2mm->ProjectionX("h1mm"); h1mm->Reset();
  TH1D *h1m = h2mr->ProjectionX("h1m"); h1m->Reset();
  assert(h1mx->GetNbinsX()==h1mr->GetNbinsX());
  assert(h1mm->GetNbinsX()==h1mr->GetNbinsX());
  TH1D *htmp(0);
  //double wZ = 1.9*2.4952; // Z boson width / Z boson mass (0%)
  double wZ = 1.85*2.4952; // Z boson width / Z boson mass (~2%)
  //double wZ = 1.8*2.4952; // Z boson width / Z boson mass (~2%)
  //double wZ = 1.7*2.4952; // Z boson width / Z boson mass (~2.5%)
  double kpt = 1.5; // Z pT smearing vs Z mass smearing
  for (int i = 1; i != h1dr->GetNbinsX()+1; ++i) {
    htmp = h2dr->ProjectionY("htmp",i,i);
    h1dr->SetBinContent(i, htmp->GetRMS());
    h1dr->SetBinError(i, htmp->GetRMSError());
    delete htmp;
    htmp = h2dx->ProjectionY("htmp",i,i);
    h1dx->SetBinContent(i, htmp->GetRMS());
    h1dx->SetBinError(i, htmp->GetRMSError());
    delete htmp;
    htmp = h2dm->ProjectionY("htmp",i,i);
    if (htmp->GetMean()>0) {
      h1dm->SetBinContent(i, sqrt(max(pow(htmp->GetRMS(),2) - wZ*wZ, 0.)) /
			  htmp->GetMean() * kpt);
      h1dm->SetBinError(i, htmp->GetRMSError() / htmp->GetMean());
    }
    delete htmp;
    h1d->SetBinContent(i, sqrt(max(pow(h1dr->GetBinContent(i),2) -
				   pow(h1dx->GetBinContent(i),2) -
				   pow(h1dm->GetBinContent(i),2), 0.)));
    // y = sqrt(a^2 + b^2 + c^2)
    // => dy = 2*a*0.5/y (oplus) 2*b*0.5/y (oplus 2*c*0.5/y
    if (h1d->GetBinContent(i)!=0) {
      h1d->SetBinError(i, sqrt(pow(h1dr->GetBinError(i),2) +
			       pow(h1dx->GetBinError(i),2) +
			       pow(h1dm->GetBinError(i),2))
		       / h1d->GetBinContent(i));
    }

    htmp = h2mr->ProjectionY("htmp",i,i);
    h1mr->SetBinContent(i, htmp->GetRMS());
    h1mr->SetBinError(i, htmp->GetRMSError());
    delete htmp;
    htmp = h2mx->ProjectionY("htmp",i,i);
    h1mx->SetBinContent(i, htmp->GetRMS());
    h1mx->SetBinError(i, htmp->GetRMSError());
    delete htmp;
    htmp = h2mm->ProjectionY("htmp",i,i);
    if (htmp->GetMean()>0) {
      h1mm->SetBinContent(i, sqrt(max(pow(htmp->GetRMS(),2) - wZ*wZ, 0.)) /
			  htmp->GetMean() * kpt);
      h1mm->SetBinError(i, htmp->GetRMSError() / htmp->GetMean());
    }
    delete htmp;
    h1m->SetBinContent(i, sqrt(max(pow(h1mr->GetBinContent(i),2) -
				   pow(h1mx->GetBinContent(i),2) -
				   pow(h1mm->GetBinContent(i),2), 0.)));
    if (h1m->GetBinContent(i)!=0) {
      h1m->SetBinError(i, sqrt(pow(h1mr->GetBinError(i),2) +
			       pow(h1mx->GetBinError(i),2) +
			       pow(h1mm->GetBinError(i),2))
		       / h1d->GetBinContent(i));
    }
  } // for i in bins

  TH1D *h = tdrHist("h","RMS",0,1,"p_{T,Z} (GeV)",15,1500);
  lumi_13TeV = "UL2018, 59.9 fb^{-1}";
  TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);
  gPad->SetLogx();

  tdrDraw(h1dr,"Pz",kFullCircle,kRed);
  tdrDraw(h1dx,"Pz",kFullCircle,kBlue);  h1dx->SetMarkerSize(0.8);
  tdrDraw(h1dm,"Pz",kFullStar,kMagenta+2);
  tdrDraw(h1d,"Pz",kFullDiamond,kGreen+2);  h1d->SetMarkerSize(1.5);

  tdrDraw(h1mr,"Pz",kOpenCircle,kRed);
  tdrDraw(h1mx,"Pz",kOpenCircle,kBlue);  h1mx->SetMarkerSize(0.8);
  tdrDraw(h1mm,"Pz",kOpenStar,kMagenta+2);
  tdrDraw(h1m,"Pz",kFullDiamond,kGreen+1);  h1m->SetMarkerSize(1.5);

  TLegend *leg = tdrLeg(0.50,0.90-8*0.05,0.70,0.90);
  leg->AddEntry(h1dr,"Z+jet MPF Data","PLE");
  leg->AddEntry(h1mr,"Z+jet MPF MC","PLE");
  leg->AddEntry(h1dx,"Z+jet MPF-X Data","PLE");
  leg->AddEntry(h1mx,"Z+jet MPF-X MC","PLE");
  leg->AddEntry(h1d,"Z+jet JER Data","PLE");
  leg->AddEntry(h1m,"Z+jet JER MC","PLE");
  leg->AddEntry(h1dm,"Z+jet Mass Data","PLE");
  leg->AddEntry(h1mm,"Z+jet Mass MC","PLE");

  gPad->RedrawAxis();
  c1->SaveAs("pdf/jerZjet/jerZjet_c1.pdf");

  // Retrieve dijet data for comparison
  TFile *fd = new TFile("rootfiles/jerCombo/dijet.root","READ");
  assert(fd && !fd->IsZombie());
  TH1D *hdj = (TH1D*)fd->Get("data_JER_standard_SM1"); assert(hdj);
  TH1D *hmj = (TH1D*)fd->Get("MC_JER_standard_SM1"); assert(hmj);

  // Retrieve RC noise for adding to Z+jet S terms
  TFile *fn = new TFile("rootfiles/jerCombo/RC.root");
  assert(fn && !fn->IsZombie());
  TGraphAsymmErrors *gdn = (TGraphAsymmErrors*)fn->Get("Data/RMS"); assert(gdn);
  TGraphAsymmErrors *gmn = (TGraphAsymmErrors*)fn->Get("MC/RMS");   assert(gmn);
  double dtn(0), mcn(0);
  for (int i = 0; i != gdn->GetN(); ++i) {
    double eta = gdn->GetX()[i];
    double dtrms = gdn->GetY()[i];
    double mcrms = gmn->GetY()[i];
    if (eta < 1.3) {
      dtn = sqrt((dtn*dtn*i + dtrms*dtrms) / (i+1));
      mcn = sqrt((mcn*mcn*i + mcrms*mcrms) / (i+1));
    }
  } // for i in gdn

  h1m->Scale(1./1.14); // Undo UL18_V2 JER SF (roughly 1.14 for |eta|<1.3)
  TH1D *h1ds = (TH1D*)h1d->Clone("h1ds");
  TH1D *h1ms = (TH1D*)h1m->Clone("h1ms");
  for (int i = 1; i != h1ds->GetNbinsX()+1; ++i) {
    double pt = h1d->GetBinCenter(i);
    h1ds->SetBinContent(i, sqrt(pow(h1d->GetBinContent(i),2)+dtn*dtn/(pt*pt)));
    h1ms->SetBinContent(i, sqrt(pow(h1m->GetBinContent(i),2)+mcn*mcn/(pt*pt)));
  } // for i in h1ds

  TH1D *h2 = tdrHist("h2","RMS",0,0.35,"p_{T,ref} (GeV)",15,1500);
  TCanvas *c2 = tdrCanvas("c2",h2,4,11,kSquare);
  gPad->SetLogx();

  if (hdj->GetListOfFunctions()!=0) hdj->GetListOfFunctions()->Delete();
  if (hmj->GetListOfFunctions()!=0) hmj->GetListOfFunctions()->Delete();
  tdrDraw(hdj,"Pz",kFullCircle,kGreen+2);
  tdrDraw(hmj,"Pz",kOpenCircle,kGreen+2);  

  //h1d->GetXaxis()->SetRangeUser(15,130.); // v40 a100
  //h1m->GetXaxis()->SetRangeUser(15,130.); // v40 a100
  h1d->GetXaxis()->SetRangeUser(15,500.); // v41 pt15
  h1m->GetXaxis()->SetRangeUser(15,500.); // v41 pt15
  tdrDraw(h1d,"Pz",kFullDiamond,kRed); h1d->SetMarkerSize(1.5);
  tdrDraw(h1m,"Pz",kOpenDiamond,kRed); h1m->SetMarkerSize(1.5);

  //h1ds->GetXaxis()->SetRangeUser(15,130.); // v40 a100
  //h1ms->GetXaxis()->SetRangeUser(15,130.); // v40 a100
  h1ds->GetXaxis()->SetRangeUser(15,500.); // v41 pt15
  h1ms->GetXaxis()->SetRangeUser(15,500.); // v41 pt15
  tdrDraw(h1ds,"Pz",kFullDiamond,kMagenta+2); h1ds->SetMarkerSize(1.5);
  tdrDraw(h1ms,"Pz",kOpenDiamond,kMagenta+2); h1ms->SetMarkerSize(1.5);

  TLegend *leg2 = tdrLeg(0.42,0.90-6*0.05,0.62,0.90);
  leg2->AddEntry(h1ds,"Z+jet JER (S #oplus N) Data","PLE");
  leg2->AddEntry(h1ms,"Z+jet JER (S #oplus N) MC","PLE");
  leg2->AddEntry(h1d,"Z+jet JER (S only) Data","PLE");
  leg2->AddEntry(h1m,"Z+jet JER (S only) MC","PLE");
  leg2->AddEntry(hdj,"Dijet JER (S #oplus N) Data","PLE");
  leg2->AddEntry(hmj,"Dijet JER (S #oplus N) MC","PLE");

  gPad->RedrawAxis();
  c2->SaveAs("pdf/jerZjet/jerZjet_c2.pdf");

  TFile *fout = new TFile("rootfiles/jerZjet.root","RECREATE");
  assert(fout && !fout->IsZombie());
  h1d->Write("zjet_jer_s_data",TObject::kOverwrite);
  h1m->Write("zjet_jer_s_mc",TObject::kOverwrite);
  h1ds->Write("zjet_jer_sn_data",TObject::kOverwrite);
  h1ms->Write("zjet_jer_sn_mc",TObject::kOverwrite);
  fout->Write();
  fout->Close();

} // jerZjet

void jerZjet_old() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  //TFile *f = new TFile("rootfiles/jecdata2018ABCD.root","READ");
  TFile *f = new TFile("rootfiles/jecdata2017BCDEF.root","READ");
  //TFile *f = new TFile("rootfiles/jecdataRun2Test.root","READ");
  assert(f && !f->IsZombie());

  f->cd("data");
  gDirectory->cd("eta00-13");
  TDirectory *dd = gDirectory;
  f->cd("mc");
  gDirectory->cd("eta00-13");
  TDirectory *dm = gDirectory;
  
  TH1D *hds = (TH1D*)dd->Get("counts_zmmjet_a100"); assert(hds);
  TH1D *hms = (TH1D*)dm->Get("counts_zmmjet_a100"); assert(hms);

  TGraphErrors *gdm0 = (TGraphErrors*)dd->Get("mpfchs1_zmmjet_a100");assert(gdm0);
  TGraphErrors *gdm1 = (TGraphErrors*)dd->Get("mpf1_zmmjet_a100"); assert(gdm1);
  TGraphErrors *gdmn = (TGraphErrors*)dd->Get("mpfn_zmmjet_a100"); assert(gdmn);
  TGraphErrors *gdmu = (TGraphErrors*)dd->Get("mpfu_zmmjet_a100"); assert(gdmu);

  TGraphErrors *gds0 = new TGraphErrors(gdm1->GetN());
  TGraphErrors *gds1 = new TGraphErrors(gdm1->GetN());
  TGraphErrors *gdsn = new TGraphErrors(gdm1->GetN());
  TGraphErrors *gdsu = new TGraphErrors(gdm1->GetN());
  for (int i = 0; i != gdm1->GetN(); ++i) {
    int j = hds->FindBin(gdm1->GetX()[i]);
    // epsilon = rms/sqrt(N_eff = >rms = epsilon*sqrt(N_eff)
    double rms0 = gdm0->GetEY()[i] * sqrt(hds->GetBinContent(j));
    gds0->SetPoint(i, gdm0->GetX()[i], rms0);
    double rms1 = gdm1->GetEY()[i] * sqrt(hds->GetBinContent(j));
    gds1->SetPoint(i, gdm1->GetX()[i], rms1);
    double rmsn = gdmn->GetEY()[i] * sqrt(hds->GetBinContent(j));
    gdsn->SetPoint(i, gdm1->GetX()[i], rmsn);
    double rmsu = gdmu->GetEY()[i] * sqrt(hds->GetBinContent(j));
    gdsu->SetPoint(i, gdm1->GetX()[i], rmsu);
  } // for i

  TGraphErrors *gmm0 = (TGraphErrors*)dm->Get("mpfchs1_zmmjet_a100");assert(gmm0);
  TGraphErrors *gmm1 = (TGraphErrors*)dm->Get("mpf1_zmmjet_a100"); assert(gmm1);
  TGraphErrors *gmmn = (TGraphErrors*)dm->Get("mpfn_zmmjet_a100"); assert(gmmn);
  TGraphErrors *gmmu = (TGraphErrors*)dm->Get("mpfu_zmmjet_a100"); assert(gmmu);

  TGraphErrors *gms0 = new TGraphErrors(gdm1->GetN());
  TGraphErrors *gms1 = new TGraphErrors(gdm1->GetN());
  TGraphErrors *gmsn = new TGraphErrors(gdm1->GetN());
  TGraphErrors *gmsu = new TGraphErrors(gdm1->GetN());
  for (int i = 0; i != gdm1->GetN(); ++i) {
    int j = hds->FindBin(gdm1->GetX()[i]);
    double k = 1;//0.5; // N_eff / N estimate
    // epsilon = rms/sqrt(N_eff = >rms = epsilon*sqrt(N_eff)
    double rms0 = gmm0->GetEY()[i] * sqrt(k*hms->GetBinContent(j));
    gms0->SetPoint(i, gmm1->GetX()[i], rms0);
    double rms1 = gmm1->GetEY()[i] * sqrt(k*hms->GetBinContent(j));
    gms1->SetPoint(i, gmm1->GetX()[i], rms1);
    double rmsn = gmmn->GetEY()[i] * sqrt(k*hms->GetBinContent(j));
    gmsn->SetPoint(i, gmm1->GetX()[i], rmsn);
    double rmsu = gmmu->GetEY()[i] * sqrt(k*hms->GetBinContent(j));
    gmsu->SetPoint(i, gmm1->GetX()[i], rmsu);
  } // for i

  gds0->Draw("AP");    gds0->SetMarkerColor(kRed);
  gds1->Draw("SAMEP"); gds1->SetMarkerColor(kGreen+1);
  gdsn->Draw("SAMEP"); gdsn->SetMarkerColor(kBlue);
  gdsu->Draw("SAMEP"); gdsu->SetMarkerColor(kCyan+1);

  gms0->Draw("SAMEP");    gms0->SetMarkerColor(kRed);
  gms1->Draw("SAMEP"); gms1->SetMarkerColor(kGreen+1);
  gmsn->Draw("SAMEP"); gmsn->SetMarkerColor(kBlue);
  gmsu->Draw("SAMEP"); gmsu->SetMarkerColor(kCyan+1);
  gms0->SetMarkerStyle(kOpenCircle);
  gms1->SetMarkerStyle(kOpenCircle);
  gmsn->SetMarkerStyle(kOpenCircle);
  gmsu->SetMarkerStyle(kOpenCircle);
  
  gPad->SetLogx();
} // jerZjet_old
