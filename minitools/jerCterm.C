// Purpose: derive JER C-term scale factor from modified hotzone maps
#include "TFile.h"
#include "TProfile2D.h"
#include "TRandom3.h"
#include "TLatex.h"
#include "TMath.h"
#include "TGraphErrors.h"

#include <fstream>

#include "../tdrstyle_mod15.C"

bool doJetVeto = true;
bool doColdVeto = true;
bool doZoom = true;
bool plotBinPDF = false;
bool plotEtaPDF = true;

TH1D* jerCterms(string trg, string mod, string iov, string data);
void drawJERSF();

void jerCterm() {
  //jerCterms("jt0","asymm","ABCD","DATA");
  //jerCterms("jt0","mpf","ABCD","DATA");
  /*
  jerCterms("jt40","asymm");
  */
  //jerCterms("jt200","asymm","ABCD","DATA");
  //jerCterms("jt200","mpf");

  //jerCterms("jt500","asymm","ABCD","DATA");
  //jerCterms("jt500","asymm","ABCD","MC");
  //jerCterms("jt500","mpf","ABCD","DATA");
  //jerCterms("jt500","mpf","ABCD","MC");
  //jerCterms("jt500","asymm","ABCD","HT");
  //jerCterms("jt500","mpf","ABCD","HT");
  /*
  jerCterms("jt500","asymm","A","DATA");
  jerCterms("jt500","asymm","B","DATA");
  jerCterms("jt500","asymm","C","DATA");
  jerCterms("jt500","asymm","D","DATA");
  */


  /*
  jerCterms("jt0","asymm","ABCD","DATA");
  jerCterms("jt40","asymm","ABCD","DATA");
  jerCterms("jt80","asymm","ABCD","DATA");
  jerCterms("jt140","asymm","ABCD","DATA");
  jerCterms("jt200","asymm","ABCD","DATA");
  jerCterms("jt260","asymm","ABCD","DATA");
  jerCterms("jt320","asymm","ABCD","DATA");
  jerCterms("jt400","asymm","ABCD","DATA");
  jerCterms("jt450","asymm","ABCD","DATA");
  jerCterms("jt500","asymm","ABCD","DATA");

  jerCterms("jt0","asymm","ABCD","MC");
  jerCterms("jt40","asymm","ABCD","MC");
  jerCterms("jt80","asymm","ABCD","MC");
  jerCterms("jt140","asymm","ABCD","MC");
  jerCterms("jt200","asymm","ABCD","MC");
  jerCterms("jt260","asymm","ABCD","MC");
  jerCterms("jt320","asymm","ABCD","MC");
  jerCterms("jt400","asymm","ABCD","MC");
  jerCterms("jt450","asymm","ABCD","MC");
  jerCterms("jt500","asymm","ABCD","MC");
  */
  if (false) {
    
    setTDRStyle();
    TDirectory *curdir = gDirectory;

    TH1D *h = new TH1D("hs",";Jet #eta;Tower-level RMS (%);",60,0,5.191);
    h->GetYaxis()->SetRangeUser(0+1e-5,8-1e-5);
    lumi_13TeV = Form("%s 2018%s %s - RMS of %s means",
		      "All trigs","ABCD","","DB");
    TCanvas *cs = tdrCanvas("cs",h,4,11,kSquare);

    // Systematic scan of all triggers vs eta
    const string vtrg[] = {/*"jt0",*/"jt40","jt80","jt140","jt200","jt260",
			   "jt320","jt400","jt450","jt500"};
    map<string, int> color;
    color["jt500"] = kBlack;
    color["jt450"] = kBlue+3;
    color["jt400"] = kBlue+2;
    color["jt320"] = kBlue+1;
    color["jt260"] = kBlue;
    color["jt200"] = kBlue-4;
    color["jt140"] = kBlue-7;
    color["jt80"] = kBlue-9;
    color["jt40"] = kBlue-10;
    color["jt0"] =  kGreen+2;
    const int ntrg = sizeof(vtrg)/sizeof(vtrg[0]);

    TLegend *legm = tdrLeg(0.18,0.78-0.035*(ntrg+2),0.38,0.78);
    legm->SetTextSize(0.035);
    legm->SetHeader("MC");
    TLegend *legd = tdrLeg(0.23,0.78-0.035*(ntrg+2),0.43,0.78);
    legd->SetTextSize(0.035);
    legd->SetHeader("  Data");

    TH1D *hd0(0), *hm0(0), *hm0b(0), *hsumw(0), *hsumwm(0);
    for (int i = 0; i != ntrg; ++i) {
      TH1D *hd = jerCterms(vtrg[i], "asymm", "ABCD", "DATA");
      TH1D *hm = jerCterms(vtrg[i], "asymm", "ABCD", "MC");

      // Initialize combined histograms
      if (!hd0) {
	hd0 = (TH1D*)hd->Clone("hd0"); hd0->Reset();
	hm0 = (TH1D*)hm->Clone("hm0"); hm0->Reset();
	hm0b = (TH1D*)hm->Clone("hm0b"); hm0b->Reset();
	hsumw = (TH1D*)hd->Clone("hsumw"); hsumw->Reset();
	hsumwm = (TH1D*)hm->Clone("hsumwm"); hsumwm->Reset();
      }
      
      // Add weighted data
      for (int ieta = 1; ieta != hd0->GetNbinsX()+1; ++ieta) {
	int ipt(0); sscanf(vtrg[i].c_str(),"jt%d",&ipt);
	double eta = hd->GetBinLowEdge(ieta);
	double ejet = ipt * cosh(eta);
	// Approximate fraction of JER measurement coming from C term
	double cd = hd->GetBinContent(ieta);
	double kw = cd / sqrt(100.*100./ipt + cd*cd);
	double cm = hm->GetBinContent(ieta);
	double kwm = cm / sqrt(100.*100./ipt + cm*cm);
	// Otherwise use statistical uncertainty as basis for weight
	double errd = hd->GetBinError(ieta);
	double w = (errd!=0 ? pow(kw,2)*1./pow(errd,2) : 0);
	double errm = hm->GetBinError(ieta);
	double wm = (errm!=0 ? pow(kwm,2)*1./pow(errm,2) : 0);
	// Add to previous bins
	double w0 = hsumw->GetBinContent(ieta);
	double sumw = w0 + w;
	if (sumw==0) continue;
	if (hd->GetBinError(ieta)==0) continue;
	if (hm->GetBinError(ieta)==0) continue;
	hsumw->SetBinContent(ieta, sumw);
	//double cd = hd->GetBinContent(ieta);
	//double errd = hd->GetBinError(ieta);
	hd0->SetBinContent(ieta, (hd0->GetBinContent(ieta)*w0 + cd*w) / sumw);
	hd0->SetBinError(ieta, sqrt(pow(hd0->GetBinError(ieta)*w0,2) + 
				    pow(errd*w,2)) / sumw);
	//double cm = hm->GetBinContent(ieta);
	//double errm = hm->GetBinError(ieta);
	hm0->SetBinContent(ieta, (hm0->GetBinContent(ieta)*w0 + cm*w) / sumw);
	hm0->SetBinError(ieta, sqrt(pow(hm0->GetBinError(ieta)*w0,2) + 
				    pow(errm*w,2)) / sumw);
	//
	double w0m = hsumwm->GetBinContent(ieta);
	double sumwm = w0m + wm;
	hsumwm->SetBinContent(ieta, sumwm);
	hm0b->SetBinContent(ieta, (hm0b->GetBinContent(ieta)*w0m+cm*wm)/sumwm);
	hm0b->SetBinError(ieta, sqrt(pow(hm0b->GetBinError(ieta)*w0m,2) + 
				     pow(errm*wm,2)) / sumwm);
      }

      cs->cd();
      tdrDraw(hd,"Pz",kFullCircle,color[vtrg[i]]);//ntrg-i);
      tdrDraw(hm,"Pz",kOpenCircle,color[vtrg[i]]);//ntrg-i);

      legd->AddEntry(hd,vtrg[i].c_str(),"PLE");
      legm->AddEntry(hm," ","PLE");
    } // for i

    tdrDraw(hd0,"Pz",kFullStar,kRed);
    tdrDraw(hm0,"Pz",kOpenStar,kRed);
    tdrDraw(hm0b,"Pz",kOpenStar,kRed-9); hm0b->SetMarkerSize(0.7);

    legd->AddEntry(hd0,"Wgtd avg.","PLE");
    legm->AddEntry(hm0," ","PLE");    

    cs->SaveAs(Form("pdf/jerCterm/jerCterm_SmearVsEta_%s_%s_%s%s.pdf",
		    "AllTrigs","asymm","ABCD",""));

    TFile *fout = new TFile("rootfiles/jerCterm.root","RECREATE");
    hd0->SetFillStyle(kNone); hd0->SetMarkerStyle(kFullCircle);
    hd0->Write("jerc_rms_data");
    hm0->SetFillStyle(kNone); hm0->SetMarkerStyle(kOpenCircle);
    hm0->Write("jerc_rms_mc");
    hm0b->SetFillStyle(kNone); hm0->SetMarkerStyle(kOpenDiamond);
    hm0b->Write("jerc_rms_mc_v2");
    fout->Close();
  } // eta scan

  cout << "Calling drawJERSF" << endl;
  drawJERSF();
} // jerCterm

TH1D* jerCterms(string trg, string mod, string iov, string data) {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  // Plotting settings
  map<string,map<string, double> > range;
  range["jt0"]["mpf"] = (doZoom ? 0.22 : 0.50);
  range["jt40"]["asymm"] = 0.22;
  range["jt200"]["mpf"] = (doZoom ? 0.08 : 0.22);
  range["jt500"]["mpf"] = (doZoom ? 0.08 : 0.22);

  //range["jt0"]["asymm"] = (doZoom ? 0.22 : 0.44);
  range["jt0"]["asymm"] = (doZoom ? 0.08 : 0.22);
  range["jt40"]["asymm"]  = (doZoom ? 0.08 : 0.22);
  range["jt80"]["asymm"]  = (doZoom ? 0.08 : 0.22);
  range["jt140"]["asymm"] = (doZoom ? 0.08 : 0.22);
  range["jt200"]["asymm"] = (doZoom ? 0.08 : 0.22);
  range["jt260"]["asymm"] = (doZoom ? 0.08 : 0.22);
  range["jt320"]["asymm"] = (doZoom ? 0.08 : 0.22);
  range["jt400"]["asymm"] = (doZoom ? 0.08 : 0.22);
  range["jt450"]["asymm"] = (doZoom ? 0.08 : 0.22);
  range["jt500"]["asymm"] = (doZoom ? 0.08 : 0.22);
  double ymax = range[trg][mod];

  // Open input file
  //TFile *f = new TFile("rootfiles/output-DATA-1-Ctest.root","READ");
  //TFile *f = new TFile("rootfiles/UL2018_JERCterm/D/output-DATA-1.root","READ");
  //TFile *f = new TFile("rootfiles/UL2018_JERCterm/ABCD/output-DATA-1.root","READ");
  //TFile *f = new TFile("rootfiles/UL2018_JERCterm/ABCD/output-MC-1.root","READ");
  const char *ci = iov.c_str();
  const char *cd = data.c_str();
  TFile *f = new TFile(Form("rootfiles/UL2018_JERCterm/%s/output-%s-1.root",
			    ci,cd),"READ");
  assert(f && !f->IsZombie());
  curdir->cd();

  // Retrieve input 2D profile
  const char *ct = trg.c_str();
  const char *cm = mod.c_str();
  string label = (mod == "mpf" ? "MPF" : "DB");
  const char *cl = label.c_str();
  TProfile2D *p2d = (TProfile2D*)f->Get(Form("FullEta_Reco/%s/p2dj%s",ct,cm));
  assert(p2d);
  TH2D *h2d = p2d->ProjectionXY(Form("h2dj%s_%s_%s%s",cm,ct,ci,cd));

  // Because h2d is R = (0.5*(p-t)) / (0.5*(p+t)), we have p/t = (1+R)/(1-R)
  // This is approximately p/t = 1 + 2*R
  // Scale h2d by x2 so that range is visually about same as for final p/t
  h2d->Scale(2.);

  if (doJetVeto) {

    // Jet veto maps
    TFile *fjv = new TFile("../JECDatabase/jet_veto_maps/Summer19UL18_V1/hotjets-UL18.root","READ");
    assert(fjv && !fjv->IsZombie());
    curdir->cd();

    TH2D *h2jv = (TH2D*)fjv->Get("h2hot_ul18_plus_hem1516_and_hbp2m1");
    assert(h2jv);

    for (int i = 1; i != h2d->GetNbinsX()+1; ++i) {
      for (int j = 1; j != h2d->GetNbinsY()+1; ++j) {
	double eta = h2d->GetXaxis()->GetBinCenter(i);
	double phi = h2d->GetYaxis()->GetBinCenter(j);
	int iv = h2jv->GetXaxis()->FindBin(eta);
	int jv = h2jv->GetYaxis()->FindBin(phi);
	if (h2jv->GetBinContent(iv,jv)>0) {
	  h2d->SetBinContent(i,j,0.);
	  h2d->SetBinError(i,j,0.);
	}
      } // for j
    } // for i

    fjv->Close();
  } // doJetVeto

  if (doColdVeto) {

    // Jet veto maps
    TFile *fjv = new TFile("rootfiles/coldjets-18runABCD.root","READ");
    assert(fjv && !fjv->IsZombie());
    curdir->cd();

    TH2D *h2jv = (TH2D*)fjv->Get("all/h2hole");
    assert(h2jv);

    for (int i = 1; i != h2d->GetNbinsX()+1; ++i) {
      for (int j = 1; j != h2d->GetNbinsY()+1; ++j) {
	double eta = h2d->GetXaxis()->GetBinCenter(i);
	double phi = h2d->GetYaxis()->GetBinCenter(j);
	int iv = h2jv->GetXaxis()->FindBin(eta);
	int jv = h2jv->GetYaxis()->FindBin(phi);
	if (h2jv->GetBinContent(iv,jv)>0) {
	  h2d->SetBinContent(i,j,0.);
	  h2d->SetBinError(i,j,0.);
	}
      } // for j
    } // for i

    fjv->Close();
  } // doJetVeto

  // Setup plotting for 2D profile
  TH1D *h0 = new TH1D(Form("h0%s_%s_%s%s",cm,ct,ci,cd),";#eta;#phi",
		      120,-5.236,5.236);
  h0->GetYaxis()->SetRangeUser(-TMath::Pi(),TMath::Pi());
  h2d->GetZaxis()->SetRangeUser(-ymax,ymax);
  lumi_13TeV = Form("%s 2018%s %s - 2D profile of %s",ct,ci,cd,cl);
  TCanvas *c0 = tdrCanvas(Form("c0%s_%s_%s%s",cm,ct,ci,cd),h0,4,11,kSquare);
  h2d->Draw("SAME COLZ");
  gPad->SetRightMargin(0.15);

  TLine *l = new TLine();
  l->SetLineStyle(kDotted);
  l->DrawLine(-1.31,-TMath::Pi(),-1.31,+TMath::Pi());
  l->DrawLine(  0.0,-TMath::Pi(),  0.0,+TMath::Pi());
  l->DrawLine(+1.31,-TMath::Pi(),+1.31,+TMath::Pi());
  
  gPad->RedrawAxis();


  // Setup plotting for eta strip scan
  TH1D *h = new TH1D(Form("h%s_%s_%s%s",cm,ct,ci,cd),
		     Form(";#phi;%s mean for probe vs tag",cl),
		     72,-TMath::Pi(),TMath::Pi());
  h->GetYaxis()->SetRangeUser(-ymax,+1.4*ymax);
  lumi_13TeV = Form("%s 2018%s %s - #eta scan of %s",ct,ci,cd,cl);
  TCanvas *c1 = tdrCanvas(Form("c1%s_%s_%s%s",cm,ct,ci,cd),h,4,11,kSquare);
  
  // Loop over barrel eta bins
  int ieta1 = h2d->GetXaxis()->FindBin(-1.3+1e-4);
  int ieta2 = h2d->GetXaxis()->FindBin(+1.3-1e-4);
  TH1D *hrms = new TH1D(Form("hrms%s_%s_%s%s",cm,ct,ci,cd),
			";Mean of probe vs tag;Towers",100,-ymax,ymax);
  TH1D *hrmsh = new TH1D(Form("hrmsh%s_%s_%s%s",cm,ct,ci,cd),
			 ";Mean of probe vs tag;Halves",100,-ymax,ymax);
  TH1D *hrmss = new TH1D(Form("hrmss%s_%s_%s%s",cm,ct,ci,cd),
			 ";Mean of probe vs tag;Strips",100,-ymax,ymax);
  TH1D *hexp = new TH1D(Form("hexp%s_%s_%s%s",cm,ct,ci,cd),
			";Expected means distribution;Towers",100,-ymax,ymax);
  TH1D *hexph = new TH1D(Form("hexph%s_%s_%s%s",cm,ct,ci,cd),
			 ";Expected means distribution;Halves",100,-ymax,ymax);
  TH1D *hexps = new TH1D(Form("hexps%s_%s_%s%s",cm,ct,ci,cd),
			 ";Expected means distribution;Strips",100,-ymax,ymax);
  TRandom3 rnd;

  // JER SFC vs eta. Bins from Summer20UL16_JRV3_MC_SF_AK4PFchs.txt
  //const double vx[] = {0, 1.3, 1.740, 1.930, 2.043,
  const double vx[] = {0, 0.522, 0.783, 1.131, 1.305, 1.740, 1.930, 2.043,
		       2.322, 2.500, 2.650, 2.853, 2.964, 3.139, 5.191};
  const int nx = sizeof(vx)/sizeof(vx[0]) - 1;
  TH2D *h2rms = new TH2D(Form("h2rms_%s_%s_%s%s",cm,ct,ci,cd),
			 ";#eta;Mean of probe vs tag;",
			 nx,vx,150,-ymax*100*1.5,ymax*100*1.5);
  TH2D *h2exp = new TH2D(Form("h2exp_%s_%s_%s%s",cm,ct,ci,cd),
			 ";#eta;Expected means distribution;",
			 nx,vx,150,-ymax*100*1.5,ymax*100*1.5);

  // Project narrow bins
  //for (int ieta = ieta1; ieta != ieta2+1; ++ieta) {
  for (int ieta = 1; ieta != h2d->GetNbinsX()+1; ++ieta) {

    // Project eta bin to phi
    TH1D *h1 = h2d->ProjectionY(Form("h1%s_%s_%s_%s_%d",cm,ct,ci,cd,ieta),
				ieta,ieta);
  
    // Modify bin contents from approximate 1+R to exact (1+R/2)/(1-R/2)
    for (int j = 1; j != h1->GetNbinsX()+1; ++j) {
      double r = h1->GetBinContent(j);
      h1->SetBinContent(j, (1+0.5*r)/(1-0.5*r) - 1);
    }

    double abseta = fabs(h2d->GetXaxis()->GetBinCenter(ieta));
    bool isb = (ieta>=ieta1 && ieta<ieta2+1);
    if (isb) tdrDraw(h1,"HISTE",kNone,ieta-ieta1+1,kSolid,-1,kNone);

    // Fill scatter of means (non-empty bins only for now)
    for (int j = 1; j != h1->GetNbinsX()+1; ++j) {
      if (h1->GetBinError(j)!=0) {
	double mean = h1->GetBinContent(j);
	double err = h1->GetBinError(j);
	if (isb) hrms->Fill(mean);
	h2rms->Fill(abseta,mean*100);
	// Simulated mean=0 scatter with x100 sampling
	for (int k = 0; k != 100; ++k) {
	  if (isb) hexp->Fill(rnd.Gaus(0., err));
	  h2exp->Fill(abseta, rnd.Gaus(0., err)*100);
	}
      }
    }

    delete h1;
  } // for ieta

  // Project wide bin on top
  if (true) {
    // Profile eta bins to phi
    TH1D *h0 = h2d->ProjectionY(Form("h0%s_%s_%s%s_2",cm,ct,ci,cd),ieta1,ieta2);
    h0->Scale(1./(ieta2-ieta1+1));

    TH1D *h1 = h2d->ProjectionY(Form("h1x%s_%s_%s_%s",cm,ct,ci,cd),ieta1,ieta2);
    h1->Clear();

    int ieta0m = h2d->GetXaxis()->FindBin(0.0-1e-4);
    TH1D *h1m = h2d->ProjectionY(Form("h1m%s_%s_%s%s",cm,ct,ci,cd),
				 ieta1,ieta0m);
    h1m->Clear();
    int ieta0p = h2d->GetXaxis()->FindBin(0.0+1e-4);
    TH1D *h1p = h2d->ProjectionY(Form("h1p%s_%s_%s%s",cm,ct,ci,cd),
				 ieta0p,ieta2);
    h1p->Clear();
  
    // Do this by hand, due to e.g. empty towers
    for (int iphi = 1; iphi != h1->GetNbinsX()+1; ++iphi) {
      TH1D *h1x = h2d->ProjectionX(Form("h1x%s_%s_%s%s_%d",cm,ct,ci,cd,iphi),
				   iphi,iphi);
      double ysum(0), ey2sum(0), wsum(0);
      double ysumm(0), ey2summ(0), wsumm(0);
      double ysump(0), ey2sump(0), wsump(0);
      for (int ieta = ieta1; ieta != ieta2+1; ++ieta) {
	if (h1x->GetBinError(ieta)!=0) {
	  ysum += h1x->GetBinContent(ieta);
	  ey2sum += pow(h1x->GetBinError(ieta),2);
	  wsum += 1;
	  if (ieta<=ieta0m) {
	    ysumm += h1x->GetBinContent(ieta);
	    ey2summ += pow(h1x->GetBinError(ieta),2);
	    wsumm += 1;
	  }
	  if (ieta>=ieta0p) {
	    ysump += h1x->GetBinContent(ieta);
	    ey2sump += pow(h1x->GetBinError(ieta),2);
	    wsump += 1;
	  }
	}
      } // for ieta
      h1->SetBinContent(iphi, wsum!=0 ? ysum / wsum : 0);
      h1->SetBinError(iphi, wsum!=0 ? sqrt(ey2sum) / wsum : 0);
      h1m->SetBinContent(iphi, wsumm!=0 ? ysumm / wsumm : 0);
      h1m->SetBinError(iphi, wsumm!=0 ? sqrt(ey2summ) / wsumm : 0);
      h1p->SetBinContent(iphi, wsump!=0 ? ysump / wsump : 0);
      h1p->SetBinError(iphi, wsump!=0 ? sqrt(ey2sump) / wsump : 0);
    } // for iphi

    // Modify bin contents from approximate 1+R to exact (1+R/2)/(1-R/2)
    for (int j = 1; j != h1->GetNbinsX()+1; ++j) {
      double r0 = h0->GetBinContent(j);
      h0->SetBinContent(j, (1+0.5*r0)/(1-0.5*r0) - 1);
      double r = h1->GetBinContent(j);
      h1->SetBinContent(j, (1+0.5*r)/(1-0.5*r) - 1);
      double rm = h1m->GetBinContent(j);
      h1m->SetBinContent(j, (1+0.5*rm)/(1-0.5*rm) - 1);
      double rp = h1p->GetBinContent(j);
      h1p->SetBinContent(j, (1+0.5*rp)/(1-0.5*rp) - 1);
    } // for j
    
    //tdrDraw(h0,"HISTE",kOpenCircle,kBlack,kSolid,-1,kNone);
    //h0->SetMarkerSize(0.5);
    tdrDraw(h1,"HISTE",kFullCircle,kBlack,kSolid,-1,kNone);
    h1->SetMarkerSize(0.5);
    tdrDraw(h1m,"HISTE",kOpenCircle,kBlue,kSolid,-1,kNone);
    h1m->SetMarkerSize(0.5);
    tdrDraw(h1p,"HISTE",kOpenCircle,kRed,kSolid,-1,kNone);
    h1p->SetMarkerSize(0.5);

    // Fill scatter of means (non-empty bins only for now)
    for (int j = 1; j != h1->GetNbinsX()+1; ++j) {
      if (h1->GetBinError(j)!=0) {
	double mean = h1->GetBinContent(j);
	double err = h1->GetBinError(j);
	hrmss->Fill(mean);
	// Simulated mean=0 scatter with x100 sampling
	for (int k = 0; k != 100; ++k) {
	  hexps->Fill(rnd.Gaus(0., err));
	}
      }
      if (h1m->GetBinError(j)!=0) {
	double mean = h1m->GetBinContent(j);
	double err = h1m->GetBinError(j);
	hrmsh->Fill(mean);
	// Simulated mean=0 scatter with x100 sampling
	for (int k = 0; k != 100; ++k) {
	  hexph->Fill(rnd.Gaus(0., err));
	}
      }
      if (h1p->GetBinError(j)!=0) {
	double mean = h1p->GetBinContent(j);
	double err = h1p->GetBinError(j);
	hrmsh->Fill(mean);
	// Simulated mean=0 scatter with x100 sampling
	for (int k = 0; k != 100; ++k) {
	  hexph->Fill(rnd.Gaus(0., err));
	}
      }
    } // for j

    //delete h1;
    //delete h1m;
    //delete h1p;
  } // Wide bin

  gPad->RedrawAxis();

  // Setup plotting for RMS
  TH1D *h2 = new TH1D(Form("h2%s_%s_%s%s",cm,ct,ci,cd),
		      Form(";%s mean for probe vs tag;"
			   "Fraction of towers/strips",cl),
		      100,-ymax,ymax);
  h2->GetYaxis()->SetRangeUser(0,0.15);
  lumi_13TeV = Form("%s 2018%s %s - RMS of %s means",ct,ci,cd,cl);
  TCanvas *c2 = tdrCanvas(Form("c2%s_%s_%s%s",cm,ct,ci,cd),h2,4,11,kSquare);

  hexp->Scale(1./hexp->Integral());
  tdrDraw(hexp,"HIST",kNone,kGray,kSolid,-1,1001,kGray);

  hexph->Scale(1./hexph->Integral());
  tdrDraw(hexph,"HIST",kNone,kGray+1,kSolid,-1,1001,kGray+1);

  hexps->Scale(1./hexps->Integral());
  tdrDraw(hexps,"HIST",kNone,kGray+2,kSolid,-1,1001,kGray+2);

  hrmss->Scale(1./hrmss->Integral());
  tdrDraw(hrmss,"HIST",kNone,kRed,kSolid,-1,1001,kRed-9);
  hrmss->SetFillColorAlpha(kRed-9,0.3);

  hrmsh->Scale(1./hrmsh->Integral());
  tdrDraw(hrmsh,"HIST",kNone,kGreen+2,kSolid,-1,1001,kGreen-9);
  hrmsh->SetFillColorAlpha(kGreen-9,0.3);

  hrms->Scale(1./hrms->Integral());
  tdrDraw(hrms,"HIST",kNone,kBlue,kSolid,-1,1001,kBlue-9);
  hrms->SetFillColorAlpha(kBlue-9,0.3);

  TLatex *tex = new TLatex();
  tex->SetTextSize(0.035);
  tex->SetNDC();
  
  // resolution_Rap9_MC_UL18V2V3_D.pdf
  // C = 0.0368 for Gaus fit, but 3 TeV RMS is about 20% worse
  double c = 0.044; // 0.0368*1.2
  double jert = sqrt(max(pow(hrms->GetRMS(),2)  - pow(hexp->GetRMS(),2),0.));
  double jerh = sqrt(max(pow(hrmsh->GetRMS(),2) - pow(hexph->GetRMS(),2),0.));
  double jers = sqrt(max(pow(hrmss->GetRMS(),2) - pow(hexps->GetRMS(),2),0.));
  double sfct = sqrt(c*c + jert*jert) / c;
  double sfch = sqrt(c*c + jerh*jerh) / c;
  double sfcs = sqrt(c*c + jers*jers) / c;

  tex->SetTextColor(kBlue);
  tex->DrawLatex(0.63,0.85,Form("RMS(tower)=%1.2f%%",hrms->GetRMS()*100.));
  tex->DrawLatex(0.63,0.81,Form("Exp(tower)=%1.2f%%",hexp->GetRMS()*100.));
  tex->DrawLatex(0.63,0.77,Form("SFC(tower)=%1.3f",sfct));
  tex->DrawLatex(0.17,0.73,Form("JERC(tower)=%1.2f%%",jert*100));
  tex->SetTextColor(kGreen+2);
  tex->DrawLatex(0.63,0.73,Form("RMS(halves)=%1.2f%%",hrmsh->GetRMS()*100.));
  tex->DrawLatex(0.63,0.69,Form("Exp(halves)=%1.2f%%",hexph->GetRMS()*100.));
  tex->DrawLatex(0.63,0.65,Form("SFC(halves)=%1.3f",sfch));
  tex->DrawLatex(0.17,0.69,Form("JERC(halves)=%1.2f%%",jerh*100));
  tex->SetTextColor(kRed);
  tex->DrawLatex(0.63,0.61,Form("RMS(strip)=%1.2f%%",hrmss->GetRMS()*100.));
  tex->DrawLatex(0.63,0.57,Form("Exp(strip)=%1.2f%%",hexps->GetRMS()*100.));
  tex->DrawLatex(0.63,0.53,Form("SFC(strip)=%1.3f",sfcs));
  tex->DrawLatex(0.17,0.65,Form("JERC(strip)=%1.2f%%",jers*100));

  gPad->RedrawAxis();

  if (plotBinPDF) {
    c0->SaveAs(Form("pdf/jerCterm/jerCterm_2D_%s_%s_%s%s.pdf",ct,cm,ci,cd));
    c1->SaveAs(Form("pdf/jerCterm/jerCterm_1D_%s_%s_%s%s.pdf",ct,cm,ci,cd));
    c2->SaveAs(Form("pdf/jerCterm/jerCterm_SF_%s_%s_%s%s.pdf",ct,cm,ci,cd));
  }


  // Plot DB/MPF scatter vs eta
  TH1D *h3 = new TH1D(Form("h3%s_%s_%s%s",cm,ct,ci,cd),
		      Form(";#eta;%s mean for probe vs tag (%%);",cl),
		      nx,vx);//,100,-ymax,ymax);
  h3->GetYaxis()->SetRangeUser(-ymax*100*1.5,ymax*100*1.5);
  lumi_13TeV = Form("%s 2018%s %s - RMS of %s means",ct,ci,cd,cl);
  TCanvas *c3 = tdrCanvas(Form("c3%s_%s_%s%s",cm,ct,ci,cd),h3,4,11,kSquare);

  h2rms->Scale(nx/h2rms->Integral());
  h2rms->GetZaxis()->SetRangeUser(0,0.15);
  h2rms->Draw("SAMECOLZ");
  gPad->SetRightMargin(0.15);
  gPad->RedrawAxis();

  if (plotBinPDF) {
    c3->SaveAs(Form("pdf/jerCterm/jerCterm_1DvsEta_%s_%s_%s%s.pdf",
		    ct,cm,ci,cd));
  }

  curdir->cd();

  // Determine SF vs eta
  TH1D *h4 = new TH1D(Form("h4%s_%s_%s%s",cm,ct,ci,cd),";#eta;RMS (%);",nx,vx);
  h4->GetYaxis()->SetRangeUser(0,8);//4);
  lumi_13TeV = Form("%s 2018%s %s - RMS of %s means",ct,ci,cd,cl);
  TCanvas *c4 = tdrCanvas(Form("c4_%s_%s_%s%s",cm,ct,ci,cd),h4,4,11,kSquare);

  TH1D *h1rms = new TH1D(Form("h1rms_%s_%s_%s_%s",cm,ct,ci,cd),
			 ";#eta;RMS (%)",nx,vx);
  for (int i = 1; i != h1rms->GetNbinsX()+1; ++i) {
    TH1D *hrms = h2rms->ProjectionY("tmp1",i,i);
    TH1D *hexp = h2exp->ProjectionY("tmp2",i,i);

    double rms = hrms->GetRMS();
    double erms = hrms->GetRMSError();
    double exp = hexp->GetRMS();
    double eexp = hexp->GetRMSError();

    double smear = sqrt(max(0., rms*rms - exp*exp));
    double esmear = sqrt(max(0., pow(rms+erms,2) - pow(exp-eexp,2))) - smear;

    if (esmear>0 && smear>0) {
      h1rms->SetBinContent(i, smear);
      h1rms->SetBinError(i, esmear);
    }

    delete hrms;
    delete hexp;
  } // for i

  tdrDraw(h1rms,"Pz",kFullCircle,kBlack);
  
  if (plotEtaPDF) {
    c4->SaveAs(Form("pdf/jerCterm/jerCterm_SmearVsEta_%s_%s_%s%s.pdf",
		    ct,cm,ci,cd));
  }

  f->Close();
  curdir->cd();
  
  return h1rms;
} // jerCterms


void drawJERSF() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  ifstream fd("../JRDatabase/textFiles/Summer19UL18_JRV2_MC/Summer19UL18_JRV2_MC_SF_AK4PFchs.txt");
  //ifstream fd("../JERCProtoLab/Summer19UL18/JER_SF/Summer19UL18_JRV1_MC_SF_AK4PFchs.txt"); // same as above?
  assert(fd.is_open());

  ifstream fm("../JRDatabase/textFiles/Summer19UL18_JRV2_MC/Summer19UL18_JRV2_MC_PtResolution_AK4PFchs.txt");
  assert(fm.is_open());

  TFile *fc = new TFile("rootfiles/jerCterm.root","READ");
  assert(fc && !fc->IsZombie());

  TFile *fn = new TFile("rootfiles/noiseTerm2018.root","READ");
  assert(fn && !fn->IsZombie());


  curdir->cd();

  // Read in C-term extra RMS
  TH1D *hcd = (TH1D*)fc->Get("jerc_rms_data"); assert(hcd);
  TH1D *hcm = (TH1D*)fc->Get("jerc_rms_mc"); assert(hcm);

  // Read in N-term RMS
  TGraphAsymmErrors *gnd = (TGraphAsymmErrors*)fn->Get("Data/RMS"); assert(gnd);
  TGraphAsymmErrors *gnm = (TGraphAsymmErrors*)fn->Get("MC/RMS"); assert(gnm);

  // Read in dijet SF
  char header[512];
  fd.getline(header,512);
  cout << "*"<<header<<"*"<< endl;
  double etamin, etamax, sf, sfdw, sfup;
  int ncol;
  TGraphErrors *gd = new TGraphErrors(0);
  while (fd >> etamin >> etamax >> ncol >> sf >> sfdw >> sfup) {
    assert(ncol==3);
    if (etamin>=0) { // symmetric
      int n = gd->GetN();
      gd->SetPoint(n, 0.5*(etamin+etamax), sf);
      gd->SetPointError(n, 0.5*(etamax-etamin), 0.5*(sfup-sfdw));
    }
  }
  // Read in MC truth C-term
  double rhomin, rhomax, jn, js, jc, jd;
  int ptmin, ptmax;
  fm.getline(header,512);
  cout << "*"<<header<<"*"<< endl;
  TGraphErrors *gjc = new TGraphErrors(0);
  TH1D *hjc = (TH1D*)hcd->Clone("hjc"); hjc->Reset();
  TH1D *hcx = (TH1D*)hcd->Clone("hcx"); hcx->Reset();
  while (fm >> etamin >> etamax >> rhomin >> rhomax >> ncol >>
	 ptmin >>  ptmax >> jn >> js >> jc >> jd) {
    assert(ncol==6);
    if (etamin>=0 && rhomin>13 && rhomin<14) { // symmetric(?)
      int n = gjc->GetN();
      gjc->SetPoint(n, 0.5*(etamin+etamax), 1+jc);
      gjc->SetPointError(n, 0.5*(etamax-etamin), 0.);
      int j = hjc->FindBin(0.5*(etamin+etamax));
      assert(hjc->GetBinContent(j)==0);
      double jcmin = 0;//0.04;
      jc = max(jcmin, jc);
      hjc->SetBinContent(j, jc);
      hcx->SetBinContent(j, 1+jc);
    }
  }
  // patch empty hc bin
  if (hjc->GetBinContent(hjc->FindBin(2.5))==0) {
    hjc->SetBinContent(hjc->FindBin(2.5),hjc->GetBinContent(hjc->FindBin(2.7)));
    hcx->SetBinContent(hcx->FindBin(2.5),hcx->GetBinContent(hcx->FindBin(2.7)));
  }


  // Re-interpret extra C-term RMS as MC C-term SF
  TH1D *hc1 = (TH1D*)hcd->Clone("hc1"); hc1->Reset();
  TH1D *hc2 = (TH1D*)hcd->Clone("hc2"); hc2->Reset();
  for (int i = 1; i != hc1->GetNbinsX()+1; ++i) {
    double cd = 0.01*hcd->GetBinContent(i);
    double ed = 0.01*hcd->GetBinError(i);
    double cm = 0.01*hcm->GetBinContent(i);
    double em = 0.01*hcm->GetBinError(i);
    double c0 = hjc->GetBinContent(i);
    double sfc1 = sqrt(cd*cd + c0*c0) / sqrt(cm*cm + c0*c0);
    double sfu1 = sqrt(pow(cd+ed,2) + c0*c0) / sqrt(pow(cm-em,2) + c0*c0);
    double sfc2 = sqrt(max(0.,cd*cd - cm*cm + c0*c0)) / c0;
    hc1->SetBinContent(i, sfc1);
    hc1->SetBinError(i, sfu1-sfc1);
    hc2->SetBinContent(i, sfc2);
  }
  
  // Turn data and MC into SF for N-term
  assert(gnd->GetN()==gnm->GetN());
  TGraphErrors *gn = new TGraphErrors(gnd->GetN());
  for (int i = 0; i != gn->GetN(); ++i) {
    gn->SetPoint(i, gnd->GetX()[i], gnd->GetY()[i] / gnm->GetY()[i]);
  }


  TH1D *h = tdrHist("hd","JER SF",0.8,1.5,"#eta_{jet}",0,5.191);
  lumi_13TeV = "UL2018";
  TCanvas *c1 = tdrCanvas("cd",h,4,11,kSquare);
  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(0,1,5.191,1);

  tdrDraw(gn,"Pz",kFullSquare,kOrange+2);
  tdrDraw(gd,"Pz",kFullCircle,kGreen+2);
  tdrDraw(hc1,"Pz",kFullDiamond,kBlue);
  tdrDraw(hc2,"Pz",kOpenDiamond,kBlue);
  tdrDraw(gjc,"Pz",kOpenStar,kRed-9);
  tdrDraw(hcx,"Pz",kOpenStar,kRed);

  TLegend *leg = tdrLeg(0.62,0.90-4*0.05,0.82,0.90);
  leg->SetTextSize(0.040);
  leg->AddEntry(gd,"S-term SF (DJ)","PLE");
  leg->AddEntry(hc1,"C-term SF (2D)","PLE");
  leg->AddEntry(hcx,"1+C (MC)","PLE");
  leg->AddEntry(gn,"N-term SF (RC)","PLE");
  
  c1->SaveAs("pdf/jerCterm/jerCterm_drawJERSF.pdf");
} // drawJERSF
