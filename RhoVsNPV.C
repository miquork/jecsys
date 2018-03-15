#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMultiGraph.h"

#include "tdrstyle_mod15.C"

void cleanGraph(TGraphErrors *g) {
  for (int i = g->GetN(); i != -1; --i) {
    if (g->GetEY()[i]==0 || (g->GetEY()[i]>1 || g->GetEX()[i]>1))
      g->RemovePoint(i);
  }
}

// map TProfiles
TGraphErrors* makeGraph(TProfile *p1, TProfile *p2,
			double dnpv=0, double due=0) {

  int nbinsx = max(p1->GetNbinsX(),p2->GetNbinsX());
  TGraphErrors *g = new TGraphErrors(0);
  for (int i = 1; i != nbinsx+1; ++i) {
    if (p1->GetBinCenter(i)>4 && p1->GetBinCenter(i)<33) {
      int n = g->GetN();
      g->SetPoint(n, p1->GetBinContent(i)+dnpv, p2->GetBinContent(i)+due);
      g->SetPointError(n, p1->GetBinError(i), p2->GetBinError(i));
    }
  }  

  return g;
}

// map TGraphErrors
TGraphErrors* makeGraph(TGraphErrors *g1, TGraphErrors *g2,
			double dnpv=0, double due=0) {

  TGraphErrors *g = new TGraphErrors(0);
  for (int i = 1; i != g1->GetN(); ++i) {
    for (int j = 1; j != g2->GetN(); ++j) {
    
      if (fabs(g1->GetX()[i]-g2->GetX()[j])<0.5) {
	double mu = 0.5*(g1->GetX()[i]+g2->GetX()[j]);
	double npv = g1->GetY()[i];
	double rho = g2->GetY()[j];
	if (mu>4 && mu<33) {
	  int n = g->GetN();
	  g->SetPoint(n, npv+dnpv, rho+due);
	  g->SetPointError(n, g1->GetEY()[i], g2->GetEY()[j]);
	}
      }
    }
  }

  return g;
}

void scaleGraph(TGraphErrors *g, TF1 *f) {
  for (int i = 0; i != g->GetN(); ++i) {
    g->SetPoint(i, g->GetX()[i], g->GetY()[i]/f->Eval(g->GetX()[i]));
    g->SetPointError(i, g->GetEX()[i], g->GetEY()[i]/f->Eval(g->GetX()[i]));
  }
}
void shiftGraph(TGraphErrors *g, TF1 *f) {
  for (int i = 0; i != g->GetN(); ++i) {
    g->SetPoint(i, g->GetX()[i], g->GetY()[i] - f->Eval(g->GetX()[i]));
    g->SetPointError(i, g->GetEX()[i], g->GetEY()[i]);
  }
}
void shiftGraph(TGraphErrors *g, double dy) {
  for (int i = 0; i != g->GetN(); ++i) {
    g->SetPoint(i, g->GetX()[i], g->GetY()[i] - dy);
    g->SetPointError(i, g->GetEX()[i], g->GetEY()[i]);
  }
}

// Plot <Rho> vs <NPV>, both from bins of TruePU
void RhoVsNPV() {

  TDirectory *curdir = gDirectory;
  setTDRStyle();

  double rhoDJ = 1.528;
  double rhoDJMC = rhoDJ;
  double rhoSJ = 1.;
  double rhoSJMC = rhoSJ;
  double rhoSJNR = rhoSJ;
  double rhoZmm = 1.2;
  double rhoZmmMC = rhoZmm;
  double rhoZmmHW = rhoZmm;
  double rhoGJ = 1.2;
  double rhoGJMC = rhoGJ;
  double rhoTT = 1.2;
  double rhoTTMC = rhoTT;
  double rhoZB = 0.;
  double rhoZBMC = rhoZB;
  double rhoZBNR = rhoZB;
  //double vtxeff = 0.7050; // Run I
  //double vtxeff2 = -0.001334; // Run I
  double vtxeff = 0.764614; // Run II fG
  double vtxeff2 = -0.0030028; // Run II fG
  double rholin = 0.573921; // Run II fG
  double rhoquad = -0.000112284; // Run II fG
  double nfakeDJ = 0.3;

  lumi_13TeV = "2016fG, 8.0 fb^{-1}";
  //lumi_13TeV = "2016fG(H), 8.0 (16.8) fb^{-1}";
  //lumi_13TeV = "2016fGH, 16.8 fb^{-1}";
  extraText = "Private Work";

  TH1D *h3 = new TH1D("h3",";#mu;#LTN_{PV,good}-N_{HS}#GT / #epsilon_{vtx}(#mu)",50,0,50);
  h3->SetMaximum(3.00);
  h3->SetMinimum(0.50);
  TCanvas *c3 = tdrCanvas("c3",h3,4,11,kSquare);

  TF1 *fs = new TF1("fz","([0]+[1]*x)*x",0,50);
  fs->SetParameters(vtxeff,vtxeff2);


  TH1D *h4 = new TH1D("h4",";#mu;#LT#rho-#rho_{UE}#GT / #rho_{PU}(#mu)",50,0,50);
  h4->SetMaximum(3.00);
  h4->SetMinimum(0.50);
  TCanvas *c4 = tdrCanvas("c4",h4,4,11,kSquare);

  TF1 *fpu = new TF1("fpu","([0]+[1]*x)*x",0,50);
  fpu->SetParameters(rholin,rhoquad);


  map<string, TGraphErrors* > gds;
  map<string, TGraphErrors* > gms; // RD MC
  map<string, TGraphErrors* > gns; // non-RD MC
  if (true) { // Z+jet, updated to 2016 Legacy re-reco

    TFile *fp = new TFile("rootfiles/combination_ZJet_Zmm_G_2017-12-08_wide.root","READ");
    //TFile *fp = new TFile("rootfiles/combination_ZJet_Zmm_H_2017-12-08_wide.root","READ");
    assert(fp && !fp->IsZombie());

    TGraphErrors *gmn = (TGraphErrors*)fp->Get("MC_npv_vs_npumean_CHS_a20_eta_00_13_L1L2L3");
    assert(gmn);
    TGraphErrors *gmr = (TGraphErrors*)fp->Get("MC_rho_vs_npumean_CHS_a20_eta_00_13_L1L2L3");
    assert(gmr);
    TGraphErrors *gdn = (TGraphErrors*)fp->Get("Data_npv_vs_npumean_CHS_a20_eta_00_13_L1L2L3");
    assert(gdn);
    TGraphErrors *gdr = (TGraphErrors*)fp->Get("Data_rho_vs_npumean_CHS_a20_eta_00_13_L1L2L3");
    assert(gdr);

    TGraphErrors *gm = makeGraph(gmn, gmr, -1, 0);
    TGraphErrors *gd = makeGraph(gdn, gdr, -1, 0);

    cleanGraph(gd);
    cleanGraph(gm);

    gds["zmm"] = gd;
    gms["zmm"] = gm;

    c3->cd();
    TGraphErrors *gd1 = (TGraphErrors*)gdn->Clone("gd1z"); 
    cleanGraph(gd1); shiftGraph(gd1,1.); scaleGraph(gd1,fs);
    tdrDraw(gd1,"Pz",kFullCircle,kRed);
    TGraphErrors *gm1 = (TGraphErrors*)gmn->Clone("gm1z"); 
    cleanGraph(gm1); shiftGraph(gm1,1.); scaleGraph(gm1,fs);
    tdrDraw(gm1,"Pz",kOpenCircle,kRed);

    c4->cd();
    TGraphErrors *gd2 = (TGraphErrors*)gdr->Clone("gd2z"); 
    cleanGraph(gd2); shiftGraph(gd2,1.4); scaleGraph(gd2,fpu);
    tdrDraw(gd2,"Pz",kFullCircle,kRed);
    TGraphErrors *gm2 = (TGraphErrors*)gmr->Clone("gm2z"); 
    cleanGraph(gm2); shiftGraph(gm2,1.4); scaleGraph(gm2,fpu);
    tdrDraw(gm2,"Pz",kOpenCircle,kRed);
  } // Z+jet

  if (true) { // gamma+jet, missing for 2016 Legacy re-reco

    //TFile *f = new TFile("files/gjet_rho_vs_mu_info-nobadruns.root","READ");
    //TFile *f = new TFile("rootfiles/gjet_rho_vs_mu_2016BCDEFGH.root","READ");
    TFile *fd = new TFile("rootfiles/gjet_rho_vs_mu_2016fG.root","READ");
    assert(fd && !fd->IsZombie());
    TFile *fm = new TFile("rootfiles/gjet_rho_vs_mu_2016MC.root","READ");
    assert(fm && !fm->IsZombie());
    //assert(f->cd("histos"));
    //TDirectory *fp = gDirectory;
    // ptgt180, ptgt70, allpt
    //TProfile *pp1 = (TProfile*)fp->Get("ptgt180_npv_vs_mu_MCgjetANDqcd");
    TProfile *pp1 = (TProfile*)fm->Get("NpvGood_vs_mu_eta0013_ptPhot_40_50");
    assert(pp1);
    //TProfile *pp2 = (TProfile*)fp->Get("ptgt180_rho_vs_mu_MCgjetANDqcd");
    TProfile *pp2 = (TProfile*)fm->Get("rho_vs_mu_eta0013_ptPhot_40_50");
    assert(pp2);
    //TProfile *ppd1 = (TProfile*)fp->Get("ptgt180_npv_vs_mu_DATA");
    TProfile *ppd1 = (TProfile*)fd->Get("NpvGood_vs_mu_eta0013_ptPhot_40_50");
    //175_230");//"40_50");
    assert(ppd1);
    //TProfile *ppd2 = (TProfile*)fp->Get("ptgt180_rho_vs_mu_DATA");
    TProfile *ppd2 = (TProfile*)fd->Get("rho_vs_mu_eta0013_ptPhot_40_50");
    //175_230");//40_50");
    assert(ppd2);

    TGraphErrors *gm = makeGraph(pp1, pp2, -1, 0);
    TGraphErrors *gd = makeGraph(ppd1, ppd2, -1, 0);

    cleanGraph(gd);
    cleanGraph(gm);

    gds["gj"] = gd;
    gms["gj"] = gm;
    
    c3->cd();
    TGraphErrors *gd1 = new TGraphErrors(ppd1);
    cleanGraph(gd1); shiftGraph(gd1,1.); scaleGraph(gd1,fs);
    tdrDraw(gd1,"Pz",kFullDiamond,kBlue);
    TGraphErrors *gm1 = new TGraphErrors(pp1);
    cleanGraph(gm1); shiftGraph(gm1,1.); scaleGraph(gm1,fs);
    tdrDraw(gm1,"Pz",kOpenDiamond,kBlue);

    c4->cd();
    TGraphErrors *gd2 = new TGraphErrors(ppd2);
    cleanGraph(gd2); shiftGraph(gd2,1.4); scaleGraph(gd2,fpu);
    tdrDraw(gd2,"Pz",kFullDiamond,kBlue);
    TGraphErrors *gm2 = new TGraphErrors(pp2);
    cleanGraph(gm2); shiftGraph(gm2,1.4); scaleGraph(gm2,fpu);
    tdrDraw(gm2,"Pz",kOpenDiamond,kBlue);
  } // gamma+jet

  if (false) { // ttbar

    //TFile *f = new TFile("files/ttbarue_obspu.root","READ");
    TFile *f = new TFile("files/ttbarue_truepu.root","READ");
    assert(f && !f->IsZombie());

    TProfile *pd1 = (TProfile*)f->Get("Profile_nvx_mu_dilep_data");
    assert(pd1);
    TProfile *pd2 = (TProfile*)f->Get("Profile_rho_mu_dilep_data");
    assert(pd2);    

    TProfile *pm1 = (TProfile*)f->Get("Profile_nvx_mu_dilep_mc");
    assert(pm1);
    TProfile *pm2 = (TProfile*)f->Get("Profile_rho_mu_dilep_mc");
    assert(pm2);

    TProfile *pn1 = (TProfile*)f->Get("Profile_nvx_mu_dilep_mc");
    assert(pn1);
    TProfile *pn2 = (TProfile*)f->Get("Profile_rho_mu_dilep_mc");
    assert(pn2);

    // Already vs <Npv-NHS> ?
    TGraphErrors *gd = makeGraph(pd1, pd2, 0, 0);
    TGraphErrors *gm = makeGraph(pm1, pm2, 0, 0);
    TGraphErrors *gn = makeGraph(pn1, pn2, 0, 0);

    cleanGraph(gd);
    cleanGraph(gm);
    cleanGraph(gn);

    gds["tt"] = gd;
    gms["tt"] = gm;
    gns["tt"] = gn;

    TGraphErrors *gd1 = new TGraphErrors(pd1); cleanGraph(gd1); scaleGraph(gd1,fs);
    tdrDraw(gd1,"P",kFullTriangleUp,kMagenta+2);
    TGraphErrors *gn1 = new TGraphErrors(pn1); cleanGraph(gn1); scaleGraph(gn1,fs);
    tdrDraw(gn1,"P",kOpenTriangleUp,kMagenta+2);
  } // ttbar


  const int ntrg = 7;
  const char *trg[ntrg] = {"jt40","jt80","jt140","jt200","jt260",
			   "jt320","jt400"};
  if (false) { // dijet

    for (int itrg = 0; itrg != ntrg; ++itrg) {

      const char *ct = trg[itrg];
      TFile *fd = new TFile("files/output-DATA-1.root","READ");
      assert(fd && !fd->IsZombie());
      TProfile *pd1 = (TProfile*)fd->Get(Form("Standard/Eta_0.0-1.3/%s"
					      "/pnpvvstrpu",ct));
      assert(pd1);
      TProfile *pd2 = (TProfile*)fd->Get(Form("Standard/Eta_0.0-1.3/%s"
					      "/prhovstrpu",ct));
      assert(pd2);    
      
      TFile *fm = new TFile("files/output-MC-1.root","READ");
      assert(fm && !fm->IsZombie());
      TProfile *pm1 = (TProfile*)fm->Get(Form("Standard/Eta_0.0-1.3/%s"
					      "/pnpvvstrpu",ct));
      assert(pm1);
      TProfile *pm2 = (TProfile*)fm->Get(Form("Standard/Eta_0.0-1.3/%s"
					      "/prhovstrpu",ct));
      assert(pm2);

      TFile *fh = new TFile("files/output-HW-1.root","READ");
      assert(fh && !fh->IsZombie());
      TProfile *ph1 = (TProfile*)fh->Get(Form("Standard/Eta_0.0-1.3/%s"
					      "/pnpvvstrpu",ct));
      assert(ph1);
      TProfile *ph2 = (TProfile*)fh->Get(Form("Standard/Eta_0.0-1.3/%s"
					      "/prhovstrpu",ct));
      assert(ph2);
      
      TGraphErrors *gd = makeGraph(pd1, pd2, -1, 0);
      TGraphErrors *gm = makeGraph(pm1, pm2, -1, 0);
      TGraphErrors *gh = makeGraph(ph1, ph2, -1, 0);
      
      cleanGraph(gd);
      cleanGraph(gm);
      cleanGraph(gh);

      gds[Form("dj_%s",ct)] = gd;
      gms[Form("dj_%s",ct)] = gm;
      gns[Form("dj_%s",ct)] = gh;
     
      if (string(ct)=="jt400") {

	TF1 *fs = new TF1(Form("fs1_%s",ct),"x/(x-[2])*([0]+[1]*x)*x",0,50);
	fs->SetParameters(vtxeff,vtxeff2,1+nfakeDJ);
	TGraphErrors *gd1 = new TGraphErrors(pd1); cleanGraph(gd1); scaleGraph(gd1,fs);
	tdrDraw(gd1,"P",kFullCircle,kBlack);
	TGraphErrors *gm1 = new TGraphErrors(pm1); cleanGraph(gm1); scaleGraph(gm1,fs);
	tdrDraw(gm1,"P",kOpenCircle,kBlack);
      }
    } // for itrg

    gds["dj"] = gds["dj_jt400"]; assert(gds["dj"]);
    gms["dj"] = gms["dj_jt400"]; assert(gms["dj"]);

    gds["sj"] = gds["dj_jt40"]; assert(gds["sj"]);
    gms["sj"] = gms["dj_jt40"]; assert(gms["sj"]);
    gns["sj"] = gns["dj_jt40"]; assert(gns["sj"]);
  } // dijet

  if (true) { // ZB soft jet (sj), updated for 2016 Legacy re-reco

    //TFile *fd = new TFile("rootfiles/zerobiasL1DataGH.root","READ");
    TFile *fd = new TFile("rootfiles/zerobiasL1DataG.root","READ");
    assert(fd && !fd->IsZombie());
    TProfile *pd1 = (TProfile*)fd->Get("pnpv30");
    assert(pd1);
    TProfile *pd2 = (TProfile*)fd->Get("prho30");
    assert(pd2);    
      
    TFile *fm = new TFile("rootfiles/zerobiasL1MCMC.root","READ");
    assert(fm && !fm->IsZombie());
    TProfile *pm1 = (TProfile*)fm->Get("pnpv30");
    assert(pm1);
    TProfile *pm2 = (TProfile*)fm->Get("prho30");
    assert(pm2);

    TGraphErrors *gd = makeGraph(pd1, pd2, -1, 0);
    TGraphErrors *gm = makeGraph(pm1, pm2, -1, 0);
      
    cleanGraph(gd);
    cleanGraph(gm);

    gds["sj"] = gd;
    gms["sj"] = gm;

    c3->cd();
    TGraphErrors *gd1 = new TGraphErrors(pd1);
    cleanGraph(gd1); shiftGraph(gd1,1.); scaleGraph(gd1,fs);
    tdrDraw(gd1,"Pz",kFullSquare,kBlack);
    TGraphErrors *gm1 = new TGraphErrors(pm1);
    cleanGraph(gm1); shiftGraph(gm1,1.); scaleGraph(gm1,fs);
    tdrDraw(gm1,"Pz",kOpenSquare,kBlack);

    c4->cd();
    TGraphErrors *gd2 = new TGraphErrors(pd2);
    cleanGraph(gd2); shiftGraph(gd2,2.0); scaleGraph(gd2,fpu);
    tdrDraw(gd2,"Pz",kFullSquare,kBlack);
    TGraphErrors *gm2 = new TGraphErrors(pm2);
    cleanGraph(gm2); shiftGraph(gm2,2.0); scaleGraph(gm2,fpu);
    tdrDraw(gm2,"Pz",kOpenSquare,kBlack);
  } // ZB soft jets

  if (true) { // ZB, updated for 2016 Legacy re-reco

    const char *p = " rootfiles/Legacy_zerobias2016_Aug072017_"
      "SingleNeutrinoRunIISumer16/";
    TFile *fd = new TFile(Form("%sLegacy_FG_R4.root",p),"READ");
    //TFile *fd = new TFile(Form("%sLegacy_H_R4.root",p),"READ");
    assert(fd && !fd->IsZombie());
    TProfile *pd1 = (TProfile*)fd->Get("p_nPV_nPU");
    assert(pd1);
    TProfile *pd2 = (TProfile*)fd->Get("p_rho_nPU");
    assert(pd2);    

    TFile *fm = new TFile(Form("%sSingleNeutrino_MC_R4.root",p),"READ");
    assert(fm && !fm->IsZombie());
    TProfile *pm1 = (TProfile*)fm->Get("p_nPV_nPU");
    assert(pm1);
    TProfile *pm2 = (TProfile*)fm->Get("p_rho_nPU");
    assert(pm2);

    TGraphErrors *gd = makeGraph(pd1, pd2, +0, 0);
    TGraphErrors *gm = makeGraph(pm1, pm2, +0, 0);

    cleanGraph(gd);
    cleanGraph(gm);

    gds["zb"] = gd;
    gms["zb"] = gm;

    TF1 *fsx = new TF1("fsx","[0]+([1]+[2]*x)*x",5,30);
    fsx->FixParameter(0,0);
    pd1->Fit(fsx,"QRN");
    cout << "fsx: " << fsx->GetParameter(0) << " + "
	 << fsx->GetParameter(1) << "*x + "
	 << fsx->GetParameter(2) << "*x*x" << endl;
    cout << "fsx_err: " << fsx->GetParError(0) << " + "
	 << fsx->GetParError(1) << "*x + "
	 << fsx->GetParError(2) << "*x*x" << endl;


    c3->cd();
    TGraphErrors *gd1 = new TGraphErrors(pd1);
    cleanGraph(gd1); scaleGraph(gd1,fs);
    tdrDraw(gd1,"Pz",kFullDotMedium,kGray+2);
    TGraphErrors *gm1 = new TGraphErrors(pm1);
    cleanGraph(gm1); scaleGraph(gm1,fs);
    tdrDraw(gm1,"Pz",kFullDotSmall,kGray+2);

    c4->cd();
    TGraphErrors *gd2 = new TGraphErrors(pd2);
    cleanGraph(gd2); scaleGraph(gd2,fpu);
    tdrDraw(gd2,"Pz",kFullDotMedium,kGray+2);
    TGraphErrors *gm2 = new TGraphErrors(pm2);
    cleanGraph(gm2); scaleGraph(gm2,fpu);
    tdrDraw(gm2,"Pz",kFullDotSmall,kGray+2);
  } // ZB

  c3->cd();
  TLegend *leg3a = tdrLeg(0.49,0.67,0.75,0.92);
  leg3a->SetHeader("MC");
  leg3a->AddEntry(gms["sj"]," ","PL");
  //leg3a->AddEntry(gms["dj"]," ","PL");
  leg3a->AddEntry(gms["gj"]," ","PL");
  leg3a->AddEntry(gms["zmm"]," ","PL");
  //leg3a->AddEntry(gms["tt"]," ","PL");
  leg3a->AddEntry(gms["zb"]," ","PL");
  
  TLegend *leg3b = tdrLeg(0.56,0.67,0.82,0.92);
  leg3b->SetHeader("Data");
  //leg3b->AddEntry(gds["sj"],"Jet (p_{T} > 74 GeV)","PL"); // jt40
  leg3b->AddEntry(gds["sj"],"Jet (p_{T} > 30 GeV)","PL"); // ZB jets
  //leg3b->AddEntry(gds["dj"],"Jet (p_{T} > 967 GeV)","PL"); // jt400
  //leg3b->AddEntry(gds["gj"],"#gamma+jet (p_{T} > 180 GeV)","PL");
  leg3b->AddEntry(gds["gj"],"#gamma+jet (p_{T} > 40 GeV)","PL");
  leg3b->AddEntry(gds["zmm"],"Z+jet","PL");
  //leg3b->AddEntry(gds["tt"],"t#bar{t} (non-RD)","PL");
  leg3b->AddEntry(gds["zb"],"Zero Bias","PL");

  c4->cd();
  TLegend *leg4a = tdrLeg(0.49,0.67,0.75,0.92);
  leg4a->SetHeader("MC");
  leg4a->AddEntry(gms["sj"]," ","PL");
  leg4a->AddEntry(gms["gj"]," ","PL");
  leg4a->AddEntry(gms["zmm"]," ","PL");
  leg4a->AddEntry(gms["zb"]," ","PL");
  
  TLegend *leg4b = tdrLeg(0.56,0.67,0.82,0.92);
  leg4b->SetHeader("Data");
  leg4b->AddEntry(gds["sj"],"Jet (p_{T} > 30 GeV)","PL"); // ZB jets
  leg4b->AddEntry(gds["gj"],"#gamma+jet (p_{T} > 40 GeV)","PL");
  leg4b->AddEntry(gds["zmm"],"Z+jet","PL");
  leg4b->AddEntry(gds["zb"],"Zero Bias","PL");
  
  c3->SaveAs("pdf/RhoVsNPV_NpvVsMu.pdf");
  c4->SaveAs("pdf/RhoVsNPV_RhoVsMu.pdf");

  curdir->cd();

  TH1D *h2 = new TH1D("h2",";#LTN_{PV,good} - N_{HS}#GT;"
		      "#LT#rho#GT - ZB (GeV)",25,0,25);
  h2->SetMaximum(2.8);
  h2->SetMinimum(-0.1);//1.1);
  
  TH1D *h = new TH1D("h",";#LTN_{PV,good}#GT;#LT#rho#GT (GeV)",25,0,25);
  h->SetMaximum(25.);
  h->SetMinimum(0.);

  TCanvas *c1 = tdrDiCanvas("c1",h,h2,4,0);
  h->GetYaxis()->SetTitleOffset(1.10);
  h2->GetXaxis()->SetTitleOffset(0.85);

  c1->cd(1);
  gPad->Update();

  TF1 *f1 = new TF1("f1","[0] + [1]*x + [2]*x*x",0,30);
  TF1 *f2 = new TF1("f2","[0] + [1]*x + [2]*x*x",0,30);
  TF1 *f3 = new TF1("f2","[0] + [1]*x + [2]*x*x",0,30);
  TF1 *fue = new TF1("fue","[0]+[1]*x",0,30);
  TF1 *fuemc = new TF1("fuemc","[0]+[1]*x",0,30);
  fuemc->SetLineStyle(kDashed);

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.030);

  vector<string> v;
  v.push_back("sj");
  //v.push_back("dj");
  v.push_back("zmm");
  v.push_back("gj");
  //v.push_back("tt");
  v.push_back("zb");

  /*
  int colors[] = {kBlack, kRed,kRed+2, kBlue, kGreen+2, kGray+2};
  int mcmarker[] = {kOpenSquare, kOpenCircle,kOpenCross, kOpenDiamond,
		    kOpenDiamond, kFullDotSmall};
  int dtmarker[] = {kFullSquare, kFullCircle,kFullTriangleUp, kFullDiamond,
		    kFullDiamond, kFullDotMedium};
  */
  int colors[] = {kBlack, kRed, kBlue, kGray+2};
  int mcmarker[] = {kOpenSquare, kOpenCircle, kOpenDiamond, kFullDotSmall};
  int dtmarker[] = {kFullSquare, kFullCircle, kFullDiamond, kFullDotMedium};

  gms["zb"]->Fit(f2,"QRN");
  f2->SetLineWidth(2);//4);
  f2->SetLineColor(kGray);
  f2->Draw("SAME");

  gds["zb"]->Fit(f1,"QRN");
  f1->SetLineWidth(4);
  f1->SetLineColor(kGray);
  f1->Draw("SAME");

  for (unsigned int i = 0; i != v.size(); ++i) {

    const char *c = v[i].c_str();

    TGraphErrors *gd = gds[c];
    tdrDraw(gd,"Pz",dtmarker[i],colors[i],kSolid,colors[i]);

    tex->SetTextColor(gd->GetMarkerColor());

    TGraphErrors *gm = gms[c];
    if (gm && v[i]!="tt") {
      gm->SetMarkerStyle(mcmarker[i]);
      gm->SetMarkerColor(colors[i]);
      gm->SetLineColor(colors[i]);
    }
    tdrDraw(gm,"Pz",mcmarker[i],colors[i],kSolid,colors[i]);

    c1->cd(2);
    TGraphErrors *gdr = (TGraphErrors*)gd->Clone();
    shiftGraph(gdr, f1);

    // Clean bigger error bars
    for (int i = gdr->GetN()-1; i != -1; --i) {
      if (gdr->GetEY()[i]>0.2) gdr->RemovePoint(i);
    }
    gdr->Draw("SAMEPz");

    TGraphErrors *gmr = (TGraphErrors*)gm->Clone();
    shiftGraph(gmr, f2);
    if (v[i]=="sj" || v[i]=="zmm" || v[i]=="gj") {
      
      // Clean bigger error bars
      for (int i = gmr->GetN()-1; i != -1; --i) {
	if (gmr->GetEY()[i]>0.2) gmr->RemovePoint(i);
      }
      tdrDraw(gmr,"Pz",mcmarker[i],colors[i],kSolid,colors[i]);
      gmr->SetMarkerSize(0.7);
    }

    gdr->Fit(fue,"QRN");
    TF1 *fuec = (TF1*)fue->Clone(Form("fue_%s",c));
    fuec->SetLineColor(colors[i]);
    fuec->Draw("SAME");

    gmr->Fit(fuemc,"QRN");
    TF1 *fuemcc = (TF1*)fuemc->Clone(Form("fuemc_%s",c));
    fuemcc->SetLineColor(colors[i]);
    if (v[i]!="zb") fuemc->Draw("SAME");

    if (v[i]=="sj") {
      rhoSJ = fue->Eval(14.);
      rhoSJMC = fuemc->Eval(14.);
      rhoSJNR = fuemcc->Eval(14.);
    }
    if (v[i]=="dj") {
      rhoDJ = fue->Eval(14.);
      rhoDJMC = fuemc->Eval(14.);
    }
    if (v[i]=="zmm") {
      rhoZmm = fue->Eval(14.);
      rhoZmmMC = fuemc->Eval(14.);
    }
    if (v[i]=="gj") {
      rhoGJ = fue->Eval(14.);
      rhoGJMC = fuemc->Eval(14.);
    }
    if (v[i]=="tt") {
      rhoTT = fue->Eval(14.);
      rhoTTMC = fuemc->Eval(14.);
    }
    if (v[i]=="zb") {
      rhoZB = f1->Eval(14.)/20.;
      rhoZBMC = f2->Eval(14.)/20.;
      rhoZBNR = f3->Eval(14.)/20.;
    }

    //c1->cd();
    c1->cd(1);
  }

  TLegend *leg2 = tdrLeg(0.19,0.57,0.45,0.87);
  //leg2->AddEntry(gds["sj"],"Jet (p_{T}>74 GeV)","PL"); // jt40
  leg2->AddEntry(gds["sj"],"Jet (p_{T}>30 GeV)","PL"); // ZB jets
  //leg2->AddEntry(gds["gj"],"#gamma+jet (p_{T}>180 GeV)","PL");
  leg2->AddEntry(gds["gj"],"#gamma+jet (p_{T}>40 GeV)","PL");
  leg2->AddEntry(gds["zmm"],"Z+jet (p_{T}>30 GeV)","PL");
  //leg2->AddEntry(gds["tt"],"t#bar{t}","PL");
  leg2->AddEntry(gds["zb"],"Zero Bias","PL");

  tex->SetTextColor(kBlack);
  tex->SetTextSize(0.045);
  tex->DrawLatex(0.20,0.50,"|#eta| < 1.3");
  tex->SetTextSize(0.035);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.60,0.27,Form("#rho_{UE}(Jet) = %1.2f GeV",
  				rhoSJ));
  tex->SetTextColor(kBlue);
  tex->DrawLatex(0.60,0.22,Form("#rho_{UE}(#gamma+jet) = %1.2f GeV",
  				rhoGJ));
  tex->SetTextColor(kRed);
  tex->DrawLatex(0.60,0.17,Form("#rho_{UE}(Z+jet) = %1.2f GeV",
				rhoZmm));
  //tex->SetTextColor(kGreen+2);
  //tex->DrawLatex(0.60,0.12,Form("#rho_{UE}(t#bar{t}) = %1.2f GeV",
  //				rhoTT));
  tex->SetTextColor(kGray+2);
  tex->DrawLatex(0.60,0.07,Form("#rho_{UE}(ZB) = %1.2f GeV",
				rhoZB));

  c1->cd(2);

  TLegend *leg2d = tdrLeg(0.19,0.78,0.93,0.90);
  leg2d->SetNColumns(3);
  leg2d->SetTextSize(0.045*2.);
  leg2d->AddEntry(gms["sj"],"Jet MC","PL");
  leg2d->AddEntry(gms["gj"],"#gamma+jet MC","PL");
  leg2d->AddEntry(gms["zmm"],"Z+jet MC","PL");

  curdir->cd();

  c1->Update();
  c1->SaveAs("pdf/RhovsNPV_Legacy2016G.pdf");
  //c1->SaveAs("pdf/RhovsNPV_Legacy2016GH.pdf");

  // Dijet study
  if (false) {

    TH1D *h22 = new TH1D("h22",";#LTN_{PV,good} - N_{HS}#GT;"
			"#LT#rho#GT- ZB (GeV)",30,0,30);
    h22->SetMaximum(3.5);
    h22->SetMinimum(-0.5);
    TCanvas *c2 = tdrDiCanvas("c2",h,h22,2,0); 
    h2->GetXaxis()->SetTitleOffset(0.85);
   
    int jcol[ntrg] = {kBlack, kRed, kOrange+2, kGreen+2, kCyan+2, kBlue,
		      kMagenta+2};
    
    c2->cd(1);

    TF1 *f12 = new TF1("f12","[0] + [1]*x + [2]*x*x",0,30);
    TF1 *fue = new TF1("fue2","[0]+[1]*x",0,30);
    gds["zb"]->Fit(f1,"QRN");
    f1->SetLineWidth(4);
    f1->SetLineColor(kGray);
    f1->Draw("SAME");
    tdrDraw((TGraphErrors*)gds["zb"]->Clone(), "P", kOpenCircle, kBlack);

    // Treat one ZB interaction as "hard scatter"
    TGraphErrors *gzb = (TGraphErrors*)gds["zb"]->Clone("zbue");
    for (int i = 0; i != gzb->GetN(); ++i) {
      gzb->SetPoint(i, gzb->GetX()[i]-0.70, gzb->GetY()[i]);
    }

    tdrDraw(gzb,"P",kFullCircle,kGray+1);
    TF1 *fuezb = (TF1*)fue->Clone("fue_zb");
    fuezb->SetLineColor(kGray+1);
    fuezb->Draw("SAME");

    c2->cd(2);

    TGraphErrors *gdrzb = (TGraphErrors*)gzb->Clone();
    shiftGraph(gdrzb, f1);
    gdrzb->Draw("SAMEP");

    gdrzb->Fit(fue,"QRN");
    TF1 *fuzb = (TF1*)fue->Clone("fue_zb");
    fuzb->SetLineColor(kGray+1);
    fuzb->Draw("SAME");

    for (int itrg = 0; itrg != ntrg; ++itrg) {
      
      c2->cd(1);
      
      const char *ct = Form("dj_%s",trg[itrg]);
      TGraphErrors *gd = gds[ct];
      tdrDraw(gd,"P",kFullCircle,jcol[itrg]);
      
      c2->cd(2);
      
      TGraphErrors *gdr = (TGraphErrors*)gd->Clone();
      shiftGraph(gdr, f1);
      gdr->Draw("SAMEP");

      gdr->Fit(fue,"QRN");
      TF1 *fuec = (TF1*)fue->Clone(Form("fue_%s",ct));
      fuec->SetLineColor(jcol[itrg]);
      fuec->Draw("SAME");
    } // for itrg

    //delete c2;
  } // jet trigger comparisons
  
}
