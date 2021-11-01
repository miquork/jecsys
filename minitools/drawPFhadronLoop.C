// Purpose draw and analyze results produced by minitools/PFhadronLoop.C
// (called by minitools/mk_PFhadronLoop.C)
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TLine.h"

#include "../tdrstyle_mod15.C"
#include "../tools.C"

void drawPFhadronLoop() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  //TFile *f = new TFile("rootfiles/Hadrons.root","READ");
  //TFile *f = new TFile("rootfiles/Hadrons_18MC.root","READ");
  //TFile *f = new TFile("rootfiles/Hadrons_18DT.root","READ");
  //TFile *f = new TFile("rootfiles/Hadrons_18.root","READ");
  //TFile *f = new TFile("rootfiles/Hadrons_18_v2.root","READ");
  TFile *f = new TFile("rootfiles/Hadrons_1718.root","READ");
  //TFile *f = new TFile("rootfiles/Hadrons_161718.root","READ");
  //TFile *f = new TFile("rootfiles/Hadrons_16GH1718.root","READ");
  assert(f && !f->IsZombie());
  //f->cd("MC18");
  f->cd("MC1718");
  //f->cd("MC161718");
  //f->cd("MC16GH1718");
  TDirectory *dm = gDirectory;
  //f->cd("DT18");
  f->cd("DT1718");
  //f->cd("DT161718");
  //f->cd("DT16GH1718");
  TDirectory *dd = gDirectory;
  curdir->cd();

  // First figure out background from high pT sideband extrapolated down
  TH2D *h2pd = (TH2D*)dd->Get("h2p");
  TH2D *h2pm = (TH2D*)dm->Get("h2p");

  // Functional shapes for some features seen in data and MC
  double pt0 = 4.27;
  double ptmin = 2.5;  

  TF1 *fp5 = new TF1("fp5","[0]/x",3,130);
  fp5->SetParameter(0,1.2);

  TF1 *fmin = new TF1("fmin","[0]/x",3,130); // window min
  fmin->SetParameter(0,1.2);
  TF1 *flo = new TF1("flo","[0]/x",3,130); // soft particle peak
  flo->SetParameter(0,pt0);
  TF1 *fhi = new TF1("fhi","[0]/sqrt(x)+[1]",3,130); // hard particle peak
  fhi->SetParameters(1.0,1.0);
  TF1 *fmax = new TF1("fmax","[0]/sqrt(x)+[1]",3,130); // window max
  fmax->SetParameters(+2.0,1);
  TF1 *fmin2 = new TF1("fmin2","[0]/sqrt(x)+[1]",3,130); // second window min
  fmin2->SetParameters(-2.0,1);
  TF1 *fmin3 = new TF1("fmin3","[0]/sqrt(x)+[1]",3,130); // second window min
  //fmin3->SetParameters(-0.5,0.55); // third window min
  fmin3->SetParameters(+0.5,0.25); // third window min
  //TF1 *fmax2 = new TF1("fmax2","2-[0]/x",3,130); // second window max
  //fmax2->SetParameter(0,1.2);
  TF1 *fmax3 = new TF1("fmax3","[0]+[1]*x",3,130);
  fmax3->SetParameters(1.1,0.01);

  TLatex *tex = new TLatex();
  tex->SetTextSize(0.045*1.5);
  tex->SetNDC();  

  TH1D *h0 = tdrHist("h0","(ECAL + HCAL) / Track",-0.05,2.495,
		     "Charged particle p_{T} (GeV)",3,130);
  lumi_13TeV = "UL17+UL18 ZeroBias";
  TCanvas *c0 = tdrCanvas("c0",h0,4,11,kSquare);
  gPad->SetLogx();
  gPad->SetLogz();
  h2pd->Draw("SAME COLZ");
  //h2pm->Draw("SAME BOX");
  fp5->Draw("SAME"); fp5->SetLineStyle(kDashDotted);
  fmin->Draw("SAME"); fmin->SetLineStyle(kDashDotted);
  flo->Draw("SAME"); flo->SetLineStyle(kDashed);
  fhi->Draw("SAME"); fhi->SetLineStyle(kDashed);
  fmax->Draw("SAME"); fmax->SetLineStyle(kDotted);
  fmin2->Draw("SAME"); fmin2->SetLineStyle(kDotted);
  fmin3->Draw("SAME");
  //fmax2->Draw("SAME");
  fmax3->Draw("SAME");
  gPad->RedrawAxis();

  TF1 *f0 = new TF1("f0","gaus(0)+fabs(gaus(3))+fabs(gaus(6))",0,1.5);
  /*
  TF1 *f1 = new TF1("f1","gaus(0)+fabs(gaus(3))+fabs(gaus(6))",0,1.5);
  TF1 *f1_l = new TF1("f1_l","gaus(0)",0,1.5);
  TF1 *f1_m = new TF1("f1_m","gaus(0)",0,1.5);
  TF1 *f1_h = new TF1("f1_h","gaus(0)",0,1.5);
  TF1 *f1_mc = new TF1("f1","gaus(0)+fabs(gaus(3))+fabs(gaus(6))",0,1.5);
  */
  vector<TGraphErrors*> vgd(f0->GetNpar()+1);
  vector<TGraphErrors*> vgm(f0->GetNpar()+1);
  //for (int i = 0; i != f0->GetNpar()+1; ++i) {
  for (int i = 0; i != vgd.size(); ++i) {
    vgd[i] = new TGraphErrors(h2pd->GetNbinsX());
    vgd[i]->SetNameTitle(Form("pd%d",i),Form("pd%d",i));
    vgm[i] = new TGraphErrors(h2pm->GetNbinsX());
    vgm[i]->SetNameTitle(Form("pm%d",i),Form("pm%d",i));
  }
  TGraphErrors *gd0 = new TGraphErrors(h2pd->GetNbinsX());
  TGraphErrors *gm0 = new TGraphErrors(h2pm->GetNbinsX());

  // Array of 31 bins => 7x5 (35) or drop one for 6x5 (30)
  const int nx = 7, ny = 5, width = 200;
  TCanvas *c1s = new TCanvas("c1s","c1s",nx*width,ny*width);
  c1s->Divide(nx,ny,0.,0.);

  for (int i = h2pd->GetNbinsX(); i != 0; --i) {

    double pt = h2pd->GetXaxis()->GetBinCenter(i);
    double minpt = h2pd->GetXaxis()->GetBinLowEdge(i);
    double maxpt = h2pd->GetXaxis()->GetBinLowEdge(i+1);

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

    TF1 *f1 = new TF1(Form("f1_%d",i),
		      "gaus(0)+fabs(gaus(3))+fabs(gaus(6))",0,1.5);
    TF1 *f1_l = new TF1(Form("f1_l_%d",i),"gaus(0)",0,1.5);
    TF1 *f1_m = new TF1(Form("f1_m_%d",i),"gaus(0)",0,1.5);
    TF1 *f1_h = new TF1(Form("f1_h_%d",i),"gaus(0)",0,1.5);
    TF1 *f1_mc = new TF1(Form("f1_mc_%d",i),
			 "gaus(0)+fabs(gaus(3))+fabs(gaus(6))",0,1.5);

    //f1->SetRange(ptmin/pt,1+2*rms);
    f1->SetRange(fmin->Eval(pt), fmax->Eval(pt));
    f1->SetParameters(0.014,0.94,0.24,
		      0.008, //pt0/pt,//0.25*pt0/pt,
		      flo->Eval(pt),
		      fp5->Eval(pt),
		      0.003,1.345,//1.5/sqrt(pt));//0.17);
		      fp5->Eval(pt));
    f1->FixParameter(4, flo->Eval(pt));//4./pt);
    f1->FixParameter(5, fp5->Eval(pt));
    f1->FixParameter(7, fhi->Eval(pt));//1.345);
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

    /*
    int n = f1->GetNpar();
    // mean = ((all-tail)*meancore + tail*meantail) / all
    double nd = h1pd->Integral();
    double md = h1pd->GetMean();
    double nl = f1_l->Integral(0,1.5);
    double ml = f1_l->GetParameter(1);
    vgd[n]->SetPoint(i-1, pt, (nd*md-nl*ml)/(nd-nl));
    */

    //if (fabs(pt-32.5)>1) continue;
    //if (fabs(pt-21)>1) continue;
    //if (fabs(pt-15)>1) continue;
    //if (fabs(pt-12)>1) continue;
    //if (fabs(pt-10)>1) continue;
    //if (fabs(pt-8)>1) continue;
    //if (fabs(pt-6)>1) continue;
    //if (fabs(pt-5)>1) continue;
    //if (fabs(pt-4.5)>1) continue;

    // Select data and MC only within window
    TH1D *h1pm2 = (TH1D*)h1pm->Clone(Form("h1pm2_%d",i));
    TH1D *h1pd2 = (TH1D*)h1pd->Clone(Form("h1pd2_%d",i));
    for (int i = 1; i != h1pm->GetNbinsX()+1; ++i) {
      double xd = h1pd->GetBinCenter(i);
      double xm = h1pm->GetBinCenter(i);
      if (/*xd<fmin->Eval(pt) ||*/ xd<fmin3->Eval(pt) ||
	  /*xd>fmax->Eval(pt) ||*/ /*xd>fmax2->Eval(pt) || */
	  xd>fmax3->Eval(pt)) {
	h1pd2->SetBinContent(i, 0.);
	h1pd2->SetBinError(i, 0.);
      }
      if (/*xm<fmin->Eval(pt) ||*/ xm<fmin3->Eval(pt) ||
	  /*xm>fmax->Eval(pt) ||*/ /*xm>fmax2->Eval(pt) ||*/
	  xm>fmax3->Eval(pt)) {
	h1pm2->SetBinContent(i, 0.);
	h1pm2->SetBinError(i, 0.);
      }
    }
    
    // Draw example distribution
    TH1D *h1 = tdrHist(Form("h1_%d",i),
		       "Fraction of particles",1e-5,0.035,//0.025,
		       "(ECAL + HCAL) / Track",-0.005,2.495);
    //TCanvas *c1 = tdrCanvas("c1",h1,4,11,kSquare);
    c1s->cd(i);
    if (i<=7 || i>= 29) h1->SetMaximum(0.045-1e-5);
    if (i>=8 && i<= 28) h1->SetMaximum(0.025-1e-5);
    h1->Draw();

    //h1pm->Draw("SAMEHIST");
    tdrDraw(h1pm,"HISTE",kNone,kRed-9,kSolid,-1,1001,kRed-9);
    tdrDraw(h1pm2,"HISTE",kNone,kBlue-9,kSolid,-1,1001,kBlue-9);
    //h1pd->Draw("SAMEPz");
    tdrDraw(h1pd,"Pz",kOpenCircle,kGray+2,kSolid,-1);
    tdrDraw(h1pd2,"Pz",kFullCircle,kBlack,kSolid,-1);
    h1pd->SetMarkerSize(0.25);
    h1pd2->SetMarkerSize(0.25);
    f1->SetLineColor(kRed);
    //f1->DrawClone("SAME");
    f1->Draw("SAME");

    f1_mc->SetLineStyle(kDashed);
    f1_mc->SetLineColor(kRed);
    //f1_mc->DrawClone("SAME");
    f1_mc->Draw("SAME");
    f1_l->SetLineColor(kBlue);
    //f1_l->DrawClone("SAME");
    f1_l->Draw("SAME");
    f1_m->SetLineColor(kGreen+2);
    //f1_m->DrawClone("SAME");
    f1_m->Draw("SAME");
    f1_h->SetLineColor(kOrange+2);
    //f1_h->DrawClone("SAME");
    f1_h->Draw("SAME");

    tex->DrawLatex(0.60,0.85,Form("%1.3g#leqp_{T}<%1.3g",minpt,maxpt));

    gPad->RedrawAxis();

    int n = f1->GetNpar();
    // Restricted range
    vgd[n]->SetPoint(i-1, pt, h1pd2->GetMean());
    vgd[n]->SetPointError(i-1, 0., h1pd2->GetMeanError());
    vgm[n]->SetPoint(i-1, pt, h1pm2->GetMean());
    vgm[n]->SetPointError(i-1, 0., h1pm2->GetMeanError());

    // Full range
    gd0->SetPoint(i-1, pt, h1pd->GetMean());
    gd0->SetPointError(i-1, 0., h1pd->GetMeanError());
    gm0->SetPoint(i-1, pt, h1pm->GetMean());
    gm0->SetPointError(i-1, 0., h1pm->GetMeanError());
    //delete c1;
  } // for i

  TH1D *h2 = tdrHist("h2","Parameter",0,2.,"p_{T} (GeV)",3,50);
  lumi_13TeV = "UL17+UL18 ZeroBias";
  TCanvas *c2 = tdrCanvas("c2",h2,4,11,kSquare);
  gPad->SetLogx();
  //for (int i = 0; i != f1->GetNpar(); ++i) {
  for (int i = 0; i != vgd.size()-1; ++i) {
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

  int n = vgd.size()-1;//f1->GetNpar();
  vgd[n]->SetMarkerStyle(kFullCircle);
  vgd[n]->SetMarkerColor(kGray+1);
  vgd[n]->SetLineColor(kGray+1);
  vgd[n]->Draw("SAMEPz");
  vgm[n]->SetMarkerStyle(kOpenCircle);
  vgm[n]->SetMarkerColor(kGray+1);
  vgm[n]->SetLineColor(kGray+1);
  vgm[n]->Draw("SAMEPz");

  fp5->Draw("SAME");

  //delete c2;
  TGraphErrors *gd = (TGraphErrors*)vgd[n]->Clone("gd");
  TGraphErrors *gm = (TGraphErrors*)vgm[n]->Clone("gm");

  TH1D *h3u = tdrHist("h3","Response",0.78,1.22,"p_{T} (GeV)",3,500);
  TH1D *h3d = tdrHist("h3","Data/MC",0.93,1.08,"p_{T} (GeV)",3,500);
  lumi_13TeV = "UL17+UL18 ZeroBias";
  TCanvas *c3 = tdrDiCanvas("c3",h3u,h3d,4,11);

  c3->cd(1);
  gPad->SetLogx();
  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(3,1,500,1);
  l->DrawLine(10,0.85,10,1.00);
  l->DrawLine(50,0.95,50,1.05);
  l->DrawLine(250,1.00,250,1.15);
  tdrDraw(gm0,"Pz",kOpenSquare,kGray+1);
  tdrDraw(gd0,"Pz",kFullSquare,kGray+1);
  tdrDraw(gm,"Pz",kOpenCircle,kBlack);
  tdrDraw(gd,"Pz",kFullCircle,kBlack);

  c3->cd(2);
  gPad->SetLogx();
  l->DrawLine(3,1,500,1);
  l->DrawLine(50,0.95,50,1.05);
  l->DrawLine(10,0.95,10,1.00);
  l->DrawLine(250,1.00,250,1.05);
  l->SetLineStyle(kDotted);
  l->DrawLine(3,0.97,500,0.97);
  l->DrawLine(3,1.03,500,1.03);
  TGraphErrors *gr0 = tools::ratioGraphs(gd0,gm0);
  tdrDraw(gr0,"Pz",kOpenSquare,kGray+1);
  TGraphErrors *gr = tools::ratioGraphs(gd,gm);
  tdrDraw(gr,"Pz",kOpenCircle,kBlack);


  c3->cd(1);

  for (int i = 0; i != gd->GetN(); ++i) {
    gd->SetPointError(i, 0., sqrt(pow(gd->GetEY()[i],2)+pow(0.003,2)));
  }
  //int nd = gd->GetN();
  int nd = 0;
  TGraphErrors *gdx = new TGraphErrors(0);
  gdx->SetPoint(nd, 50, 1.0);
  gdx->SetPointError(nd, 0., 0.003);
  gdx->SetPoint(nd+2, 100, 1.06);
  gdx->SetPointError(nd+2, 0., 0.003);
  gdx->SetPoint(nd+1, 250, 1.08);
  gdx->SetPointError(nd+1, 0., 0.003);
  //gd->SetPoint(nd, 250, 1.15);
  //gd->SetPointError(nd, 0., 0.005);
  for (int i = 0; i != gdx->GetN(); ++i) {
    int n = gd->GetN();
    gd->SetPoint(n, gdx->GetX()[i], gdx->GetY()[i]);
    gd->SetPointError(n, gdx->GetEX()[i], gdx->GetEY()[i]);
  }

  for (int i = 0; i != gm->GetN(); ++i) {
    gm->SetPointError(i, 0., sqrt(pow(gm->GetEY()[i],2)+pow(0.003,2)));
  }
  //int nm = gm->GetN();
  int nm = 0;
  TGraphErrors *gmx = new TGraphErrors(0);
  gmx->SetPoint(nm, 50, 1.0);
  gmx->SetPointError(nm, 0., 0.003);
  gmx->SetPoint(nm+2, 100, 1.03);
  gmx->SetPointError(nm+2, 0., 0.003);
  gmx->SetPoint(nm+1, 250, 1.04);
  gmx->SetPointError(nm+1, 0., 0.003);
  //gm->SetPoint(nm, 250, 1.12);
  //gm->SetPointError(nm, 0., 0.005);
  for (int i = 0; i != gmx->GetN(); ++i) {
    int n = gm->GetN();
    gm->SetPoint(n, gmx->GetX()[i], gmx->GetY()[i]);
    gm->SetPointError(n, gmx->GetEX()[i], gmx->GetEY()[i]);
  }

  tdrDraw(gmx,"Pz",kOpenCircle,kBlue);
  tdrDraw(gdx,"Pz",kFullCircle,kRed);
  
  //TF1 *fd = new TF1("fd","[0]+[1]*pow(x,[2])+[3]*pow(x,[4])",3,50);
  //fd->SetParameters(1,-16,-1.7,7.4,-1);
  //fd->FixParameter(3,-1.7);
  /*
  TF1 *fd = new TF1("fd","[0]+log(x/50.)*([1]+log(x/50.)*([2]+log(x/50.)*[3]))"
		    "+ [4]/x + [5]/(x*x)",3,50);
  fd->SetParameters(1,-0.02,-0.01,-0.001, 0.1,0.1);
  TF1 *fm = new TF1("fm","[0]+log(x/50.)*([1]+log(x/50.)*([2]+log(x/50.)*[3]))"
		    "+ [4]/x + [5]/(x*x)",3,50);
  fm->SetParameters(1,-0.02,-0.01,-0.001, 0.1,0.1);
  fm->SetLineColor(kBlue);
  */
  TF1 *fd = new TF1("fd","1.15+[0]*pow(x,-2)+[1]*pow(x,-1.5)+[2]*pow(x,-1)"
		    "+[3]*pow(x,-0.75)+[4]*pow(x,-0.5)+[5]*pow(x,-0.25)",3,500);
  fd->SetParameters(1,0.5,0.1,-0.1,-0.05,-0.01);
  TF1 *fm = new TF1("fm","1.12+[0]*pow(x,-2)+[1]*pow(x,-1.5)+[2]*pow(x,-1)"
		    "+[3]*pow(x,-0.75)+[4]*pow(x,-0.5)+[5]*pow(x,-0.25)",3,500);
  fm->SetParameters(1,0.5,0.1,-0.1,-0.05,-0.01);
  fm->SetLineColor(kBlue);

  gm->Fit(fm,"RN");
  gd->Fit(fd,"RN");

  fm->Draw("SAME");
  fd->Draw("SAME");

  c3->cd(2);

  TGraph *gf = new TGraph(0);
  for (double pt = 3.; pt < 500; pt *= 1.05) {
    int n = gf->GetN();
    gf->SetPoint(n, pt, fd->Eval(pt) / fm->Eval(pt));
  }
  tdrDraw(gf,"C",kNone,kRed,kSolid,-1);

  //TF1 *fr0 = new TF1("fr0","[0]+[1]*0.5*(1+TMath::Erf((x-[2])/[3]))",10,40);
  TF1 *fr0 = new TF1("fr0","[0]+[1]*0.5*(1+TMath::Erf((x-[2])/[3]))"
  		     "+[4]/x+[5]/(x*x)",3.9,40);
  fr0->SetParameters(0.975,0.054,18.0,4.99,0.03,-0.01);
  fr0->SetLineColor(kGreen+2);
  fr0->SetLineWidth(2);
  gr0->Fit(fr0,"RN");
  fr0->Draw("SAME");

  TF1 *fr0x = (TF1*)fr0->Clone("fr0x");
  fr0x->SetLineStyle(kDashed);
  fr0x->SetRange(3,500);
  fr0x->Draw("SAME");

  c3->cd(1);
  
  TLegend *leg3 = tdrLeg(0.40,0.90-0.05*4,0.65,0.90);
  leg3->AddEntry(gd0,"Data Profile Mean","PLE");
  leg3->AddEntry(gm0,"MC Profile Mean","PLE");
  leg3->AddEntry(gd,"Data Hist Mean","PLE");
  leg3->AddEntry(gm,"MC Hist Mean","PLE");


  //delete c0;
  //delete c1s;
  delete c2;
} // drawPFhadronLoop()
