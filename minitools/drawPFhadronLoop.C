// Purpose draw and analyze results produced by minitools/PFhadronLoop.C
// (called by minitools/mk_PFhadronLoop.C)
#include "TFile.h"
#include "TH2D.h"

#include "../tdrstyle_mod15.C"

void drawPFhadronLoop() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  //TFile *f = new TFile("rootfiles/Hadrons.root","READ");
  //TFile *f = new TFile("rootfiles/Hadrons_18MC.root","READ");
  //TFile *f = new TFile("rootfiles/Hadrons_18DT.root","READ");
  TFile *f = new TFile("rootfiles/Hadrons_18.root","READ");
  assert(f && !f->IsZombie());
  f->cd("MC18");
  TDirectory *dm = gDirectory;
  f->cd("DT18");
  TDirectory *dd = gDirectory;

  // First figure out background from high pT sideband extrapolated down
  TH2D *h2pd = (TH2D*)dd->Get("h2p");
  TH2D *h2pm = (TH2D*)dm->Get("h2p");

  // Functional shapes for some features seen in data and MC
  TF1 *fp5 = new TF1("fp5","[0]/x",3,50);
  fp5->SetParameter(0,1.2);

  TF1 *fmin = new TF1("fmin","[0]/x",3,50);
  fmin->SetParameter(0,1.2);
  TF1 *fmax = new TF1("fmax","[0]/sqrt(x)+[1]",3,50);
  fmax->SetParameters(2.0,1);

  TH1D *h0 = tdrHist("h0","(ECAL + HCAL) / Track",-0.05,2.495,
		     "Charged particle p_{T} (GeV)",3,50);
  TCanvas *c0 = tdrCanvas("c0",h0,4,11,kSquare);
  gPad->SetLogx();
  gPad->SetLogz();
  h2pd->Draw("SAME COLZ");
  //h2pm->Draw("SAME BOX");
  fp5->Draw("SAME");
  fmin->Draw("SAME");
  fmax->Draw("SAME");
  gPad->RedrawAxis();

  TF1 *f1 = new TF1("f1","gaus(0)+fabs(gaus(3))+fabs(gaus(6))",0,1.5);
  TF1 *f1_l = new TF1("f1_l","gaus(0)",0,1.5);
  TF1 *f1_m = new TF1("f1_m","gaus(0)",0,1.5);
  TF1 *f1_h = new TF1("f1_h","gaus(0)",0,1.5);
  TF1 *f1_mc = new TF1("f1","gaus(0)+fabs(gaus(3))+fabs(gaus(6))",0,1.5);
  double pt0 = 4.27;
  double ptmin = 2.5;  
  vector<TGraphErrors*> vgd(f1->GetNpar()+1);
  vector<TGraphErrors*> vgm(f1->GetNpar()+1);
  for (int i = 0; i != f1->GetNpar()+1; ++i) {
    vgd[i] = new TGraphErrors(h2pd->GetNbinsX());
    vgd[i]->SetNameTitle(Form("pd%d",i),Form("pd%d",i));
    vgm[i] = new TGraphErrors(h2pm->GetNbinsX());
    vgm[i]->SetNameTitle(Form("pm%d",i),Form("pm%d",i));
  }

  for (int i = h2pd->GetNbinsX(); i != 0; --i) {

    double pt = h2pd->GetXaxis()->GetBinCenter(i);

    TH1D *h1pd = h2pd->ProjectionY(Form("h1pd_%d",i),i,i);
    TH1D *h1pm = h2pm->ProjectionY(Form("h1pm_%d",i),i,i);

    // Patch x10 entries in 0.99 for now; figure out what's causing it later
    int j = h1pd->FindBin(0.99);
    h1pd->SetBinContent(j,0);
    h1pd->SetBinError(j,0);
    h1pd->Scale(h1pd->Integral()>0 ? 1./h1pd->Integral() : 1);

    h1pm->SetBinContent(j,0);
    h1pm->SetBinError(j,0);
    h1pm->Scale(h1pm->Integral()>0 ? 1./h1pm->Integral() : 1);

    double mean = h1pd->GetRMS();
    double rms = h1pd->GetRMS();

    f1->SetRange(ptmin/pt,1+2*rms);
    f1->SetParameters(0.014,0.94,0.24,
		      0.008,pt0/pt,//0.25*pt0/pt,
		      fp5->Eval(pt),
		      0.003,1.345,//1.5/sqrt(pt));//0.17);
		      fp5->Eval(pt));
    f1->FixParameter(4, 4./pt);
    f1->FixParameter(5, fp5->Eval(pt));
    f1->FixParameter(7, 1.345);
    f1->FixParameter(8, fp5->Eval(pt));
    //h1p->Fit(f1,"R");
    h1pd->Fit(f1,"QRN");

    for (int j = 0; j != f1->GetNpar(); ++j) {
      f1_mc->SetParameter(j, f1->GetParameter(j));
      if (j==4 || j==5 || j==7 || j==8)
	f1_mc->FixParameter(j, f1->GetParameter(j));
    }
    h1pm->Fit(f1_mc,"QRN");

    //delete h1p;
    for (int ip = 0; ip != f1->GetNpar(); ++ip) {
      double k = (ip%3==0 ? 0.5*100 : 1);
      vgd[ip]->SetPoint(i-1, pt, fabs(f1->GetParameter(ip)) * k);
      vgd[ip]->SetPointError(i-1, 0., f1->GetParError(ip) * k);
      vgm[ip]->SetPoint(i-1, pt, fabs(f1_mc->GetParameter(ip)) * k);
      vgm[ip]->SetPointError(i-1, 0., f1_mc->GetParError(ip) * k);
    }

    h1pd->Fit(f1,"QRN");
    for (int j = 0; j != 3; ++j) {
      f1_m->SetParameter(j, f1->GetParameter(0+j));
      f1_l->SetParameter(j, f1->GetParameter(3+j));
      f1_h->SetParameter(j, f1->GetParameter(6+j));
    }

    int n = f1->GetNpar();
    // mean = ((all-tail)*meancore + tail*meantail) / all
    double nd = h1pd->Integral();
    double md = h1pd->GetMean();
    double nl = f1_l->Integral(0,1.5);
    double ml = f1_l->GetParameter(1);
    vgd[n]->SetPoint(i-1, pt, (nd*md-nl*ml)/(nd-nl));

    //if (fabs(pt-32.5)>1) continue;
    //if (fabs(pt-21)>1) continue;
    //if (fabs(pt-15)>1) continue;
    //if (fabs(pt-12)>1) continue;
    //if (fabs(pt-10)>1) continue;
    //if (fabs(pt-8)>1) continue;
    if (fabs(pt-6)>1) continue;
    //if (fabs(pt-5)>1) continue;
    //if (fabs(pt-4.5)>1) continue;

    // Draw example distribution
    TH1D *h1 = tdrHist("h1","Fraction of particles",0,0.025,
		       "(ECAL + HCAL) / Track",-0.005,2.495);
    TCanvas *c1 = tdrCanvas("c1",h1,4,11,kSquare);

    //h1pm->Draw("SAMEHIST");
    tdrDraw(h1pm,"HISTE",kNone,kBlue-9,kSolid,-1,1001,kBlue-9);
    //h1pd->Draw("SAMEPz");
    tdrDraw(h1pd,"Pz",kFullCircle,kBlack,kSolid,-1);
    f1->SetLineColor(kRed);
    f1->DrawClone("SAME");

    f1_mc->SetLineStyle(kDashed);
    f1_mc->SetLineColor(kRed);
    f1_mc->DrawClone("SAME");
    f1_l->SetLineColor(kBlue);
    f1_l->DrawClone("SAME");
    f1_m->SetLineColor(kGreen+2);
    f1_m->DrawClone("SAME");
    f1_h->SetLineColor(kOrange+1);
    f1_h->DrawClone("SAME");

    gPad->RedrawAxis();
  } // for i

  TH1D *h = tdrHist("h","Parameter",0,2.,"p_{T} (GeV)",3,50);
  lumi_13TeV = "UL18";
  TCanvas *c2 = tdrCanvas("c2",h,4,11,kSquare);
  gPad->SetLogx();
  for (int i = 0; i != f1->GetNpar(); ++i) {
    //vg[i]->SetMarkerColor(i/3+1);
    if (i/3==0) vgd[i]->SetMarkerColor(kBlack); // Core
    if (i/3==1) vgd[i]->SetMarkerColor(kBlue);  // Low side
    if (i/3==2) vgd[i]->SetMarkerColor(kRed);   // High side

    if (i%3==0) vgd[i]->SetMarkerStyle(kFullSquare);  // Height
    if (i%3==1) vgd[i]->SetMarkerStyle(kFullCircle);  // Mean
    if (i%3==2) vgd[i]->SetMarkerStyle(kFullDiamond); // Width
    vgd[i]->Draw("SAMEPz");

    vgm[i]->SetMarkerColor(vgd[i]->GetMarkerColor());
    if (i%3==0) vgm[i]->SetMarkerStyle(kOpenSquare);  // Height
    if (i%3==1) vgm[i]->SetMarkerStyle(kOpenCircle);  // Mean
    if (i%3==2) vgm[i]->SetMarkerStyle(kOpenDiamond); // Width
    vgm[i]->Draw("SAMEPz");
  }

  int n = f1->GetNpar();
  vgd[n]->SetMarkerStyle(kFullCircle);
  vgd[n]->SetMarkerColor(kGray+1);
  vgd[n]->Draw("SAMEPz");

  fp5->Draw("SAME");

} // drawPFhadronLoop()
