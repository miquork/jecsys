// Purpose: Study JER SF with RC noise, dijet S and 2D C term measurements
#include "TFile.h"
#include "TGraphAsymmErrors.h"

#include "../tdrstyle_mod15.C"
#include "../tools.C"

#include <fstream>

void cleanGraph(TGraph *g, double xmin=0, double xmax=0) {
  for (int i = g->GetN()-1; i != -1; --i) {
    if (g->GetY()[i]==0 ||
	(xmax!=0 && g->GetX()[i]>xmax) ||
	(xmin!=0 && g->GetX()[i]<xmin)) g->RemovePoint(i);
  }
}

struct JER {
  double eta, deta;
  double n, s, c;
  double kn, ks, kc;
  double sfmin, sfmid, sfdif, sfmax;
  double sf8, sf80, sf3500, sfxmin, sfxmax;
};

JER jerSFs(int ieta);

TCanvas *_c0(0);
void jerSF() {

  // 5x4 of 1x1.25
  _c0 = new TCanvas("c0","c0",1250,1000);
  _c0->Divide(5,4,0,0);
  
  const int neta = 19;//14;
  TGraphErrors *ge = new TGraphErrors(neta);
  TGraphErrors *g8 = new TGraphErrors(neta);
  TGraphErrors *g80 = new TGraphErrors(neta);
  TGraphErrors *gxmin = new TGraphErrors(neta);
  TGraphErrors *gxmax = new TGraphErrors(neta);
  TGraphErrors *gc = new TGraphErrors(neta);
  TGraphErrors *gs = new TGraphErrors(neta);
  TGraphErrors *gn = new TGraphErrors(neta);

  // Plots vs pT in eta bins
  vector<JER> vjer(neta);
  for (int ieta = 0; ieta != neta; ++ieta) {
    JER jer = jerSFs(ieta);
    vjer[ieta] = jer;
    //cout << "jer.eta="<<jer.eta<<" jer.sfmid="<<jer.sfmid
    //	 << " jer.sf80="<<jer.sf80<<endl;
    ge->SetPoint(ieta, jer.eta, jer.sfmid);
    ge->SetPointError(ieta, jer.deta, jer.sfdif);
    g8->SetPoint(ieta, jer.eta, jer.sf8);
    g8->SetPointError(ieta, jer.deta, 0.);
    g80->SetPoint(ieta, jer.eta, jer.sf80);
    g80->SetPointError(ieta, jer.deta, 0.);
    gxmin->SetPoint(ieta, jer.eta, jer.sfxmin);
    gxmin->SetPointError(ieta, jer.deta, 0.);
    gxmax->SetPoint(ieta, jer.eta, jer.sfxmax);
    gxmax->SetPointError(ieta, jer.deta, 0.);
    gc->SetPoint(ieta, jer.eta, jer.kc);
    gc->SetPointError(ieta, jer.deta, 0.);
    gs->SetPoint(ieta, jer.eta, jer.ks);
    gs->SetPointError(ieta, jer.deta, 0.);
    gn->SetPoint(ieta, jer.eta, jer.kn);
    gn->SetPointError(ieta, jer.deta, 0.);
  } // for ieta

  _c0->SaveAs("pdf/jerSF/jerSF_5x4.pdf");

  // Load official JER SF as reference
  //TFile *fj = new TFile("rootfiles/jerCterm_v5.root","READ");
  TFile *fj = new TFile("rootfiles/jerCterm.root","READ");
  assert(fj && !fj->IsZombie());
  TGraphErrors *gsf = (TGraphErrors*)fj->Get("JER_SF_UL18_V2"); assert(gsf);
  TGraphErrors *gn0 = (TGraphErrors*)fj->Get("JER_N_RC"); assert(gn0);

  // Draw summary plots vs eta, overlaid width dijet SF
  TH1D *h = tdrHist("h","JER SF",0.85,1.75,"#eta_{jet}",0,5.2);
  lumi_13TeV = "UL2018, 59.9 fb^{-1}";
  TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);

  TLine *l = new TLine();
  tdrDraw(ge,"E2",kNone,kYellow+3,kSolid,-1,1001,kYellow+1);
  l->SetLineStyle(kDashed); l->DrawLine(0,1.0,5.2,1.0);
  l->SetLineStyle(kDotted); l->DrawLine(0,1.1,5.2,1.1);

  tdrDraw(gxmin,"Pz",kOpenDiamond,kRed,kSolid);
  tdrDraw(g8,"Pz",kNone,kRed,kSolid);
  tdrDraw(g80,"Pz",kFullCircle,kBlack,kSolid);
  tdrDraw(gxmax,"Pz",kOpenDiamond,kBlue,kSolid);
  tdrDraw(gsf,"Pz",kOpenCircle,kGreen+2);

  TLegend *leg = tdrLeg(0.63,0.90-5*0.05,0.83,0.90);
  leg->AddEntry(gsf,"UL18_V2","PLE");
  leg->AddEntry(g80,"p_{T} = 80 GeV","PL");
  leg->AddEntry(gxmin,"p_{T} = 8 GeV","PL");
  leg->AddEntry(gxmax,"p_{T} = E/cosh(#eta)","PL");
  leg->AddEntry(ge,"SF range","F");

  gPad->RedrawAxis();
  c1->SaveAs("pdf/jerSF/jerSF_vsEta.pdf");

  TLegend *leg2 = tdrLeg(0.20,0.70-3*0.05,0.40,0.70);
  leg2->AddEntry(gn,"RC noise SF","PL");
  leg2->AddEntry(gs,"Dijet S SF","PL");
  leg2->AddEntry(gc,"2D const. SF","PL");

  tdrDraw(gn0,"Pz",kOpenSquare,kOrange+2); gn0->SetMarkerSize(0.5);
  tdrDraw(gn,"Pz",kFullSquare,kOrange+2); gn->SetMarkerSize(0.5);
  tdrDraw(gs,"Pz",kFullSquare,kBlack); gs->SetMarkerSize(0.5);
  tdrDraw(gc,"Pz",kFullSquare,kBlue); gc->SetMarkerSize(0.5);

  gPad->RedrawAxis();
  c1->SaveAs("pdf/jerSF/jerSF_vsEta_wSF.pdf");

  // Produce output text file (FactorizedJetCorrecter Style
  ofstream txt("pdf/jerSF/Summer19UL18_JRV3_MC_SF_AK4PFchs.txt");
  //txt << "{1 JetEta 0 None ScaleFactor}" << endl;
  txt << "{1 JetEta 2 JetPt Rho  "
      << "sqrt(([0]*[0]*[3]*[3]*y/[6])/(x*x)+[1]*[1]*[4]*[4]/x+[2]*[2]*[5]*[5])"
      << "/sqrt(([0]*[0]*y/[6])/(x*x)+[1]*[1]/x+[2]*[2])"
      << " Correction L2Relative}" << endl; // test: runs
      //<< " Correction JERScaleFactor}" << endl; // test: fails
    //    << " Resolution}" << endl; // ScaleFactor does not run, nor does this
  //<< " ScaleFactor}" << endl; // doesn't run
  double rho = 20.85; // UL18 Z+jet vs pT at [60,130] GeV, |eta|<1.3
  for (int i = neta-1; i != -1; --i) {
    JER &jer = vjer[i];
    txt << Form("%5.3f %6.3f %2d %2d %4d %2d %2d %5.2f %5.3f %6.4f "
		" %5.3f %5.3f %5.3f %5.2f\n",
		-jer.eta-jer.deta, -jer.eta+jer.deta, 11, 8, 4000, 0, 70,
		jer.n, jer.s, jer.c, jer.kn, jer.ks, jer.kc, rho);
  } // for i in -neta
  for (int i = 0; i != neta; ++i) {
    JER &jer = vjer[i];
    txt << Form("%5.3f %6.3f %2d %2d %4d %2d %2d %5.2f %5.3f %6.4f "
		" %5.3f %5.3f %5.3f %5.2f\n",
		+jer.eta-jer.deta, +jer.eta+jer.deta, 11, 8, 4000, 0, 70,
		jer.n, jer.s, jer.c, jer.kn, jer.ks, jer.kc, rho);
  } // for i in -neta
  
} // jerSF

JER jerSFs(int ieta) {

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  
  TFile *fn = new TFile("rootfiles/jerCombo/RC.root","READ");
  assert(fn && !fn->IsZombie());

  TFile *fd = new TFile("rootfiles/jerCombo/dijet.root","READ");
  assert(fd && !fd->IsZombie());

  TFile *fs = new TFile("rootfiles/jerZjet.root","READ");
  assert(fs && !fs->IsZombie());

  //TFile *fc = new TFile("rootfiles/jerCombo/Cterm.root","READ");
  TFile *fc = new TFile("rootfiles/jerCombo/Cterm_v2.root","READ");
  //TFile *fc = new TFile("rootfiles/jerCombo/Cterm_v3.root","READ");
  assert(fc && !fc->IsZombie());

  curdir->cd();

  // try to map ieta to jeta for dijet
  const double xc[] = {0, 0.261, 0.522, 0.783, 1.044, 1.305,
		       1.566, 1.74, 1.93, 2.043, 2.172, 2.322, 2.5,
		       2.65, 2.853, 2.964, 3.139, 3.489, 3.839, 5.19};
  const int nxc = sizeof(xc)/sizeof(xc[0])-1;
  const double xd[] = {0, 0.522, 0.783, 1.044, 1.305,
		       1.74, 1.93, 2.043, 2.322, 2.5,
		       2.65, 2.853, 2.964, 3.139, 5.19};
  const int nxd = sizeof(xd)/sizeof(xd[0])-1;
  TH1D *hxc = new TH1D(Form("hxc_%d",ieta),"",nxc,xc);
  TH1D *hxd = new TH1D(Form("hxd_%d",ieta),"",nxd,xd);

  // Set ks, cs by hand
  const int neta = 19;
  double acs[neta];
  // ieta12 (2.5-2.65) behaving very strange, dipping in JER
  acs[0]=0.90; acs[1]=0.90; acs[2]=0.90; acs[3]=0.90; acs[4]=0.98;
  acs[5]=1.02; acs[6]=1.02; acs[7]=1.00; acs[8]=0.98; acs[9]=0.95;
  acs[10]=0.93; acs[11]=0.92; acs[12]=0.70; acs[13]=0.95;
  acs[14]=0.95; acs[15]=1.20; acs[16]=0.90; acs[17]=0.80; acs[18]=0.70;
  double aks[neta];
  aks[0]=1.14; aks[1]=1.14; aks[2]=1.15; aks[3]=1.18; aks[4]=1.19;
  //aks[0]=1.17; aks[1]=1.17; aks[2]=1.17; aks[3]=1.18; aks[4]=1.19;
  aks[5]=1.19; aks[6]=1.22; aks[7]=1.20; aks[8]=1.20; aks[9]=1.25;
  aks[10]=1.27; aks[11]=1.24; aks[12]=1.40; aks[13]=1.40;
  aks[14]=2.0; aks[15]=1.65; aks[16]=1.20; aks[17]=1.25; aks[18]=1.30;
  double acc[neta];
  acc[0]=0.040; acc[1]=0.040; acc[2]=0.043; acc[3]=0.045; acc[4]=0.057;
  acc[5]=0.055; acc[6]=0.055; acc[7]=0.045; acc[8]=0.045; acc[9]=0.047;
  acc[10]=0.049; acc[11]=0.043; acc[12]=0.040; acc[13]=0.040;
  acc[14]=0.075; acc[15]=0.117; acc[16]=0.075; acc[17]=0.085; acc[18]=0.085;
  /*
  acs[0]=0.90; acs[1]=0.92; acs[2]=0.93; acs[3]=1.00; acs[4]=1.07;
  acs[5]=0.97; acs[6]=0.97; acs[7]=0.91; acs[8]=0.92; acs[9]=0.80;
  acs[10]=1.10; acs[11]=1.00; acs[12]=1.00; acs[13]=1.00;
  acs[14]=1.00; acs[15]=1.00; acs[16]=1.00; acs[17]=1.00; acs[18]=1.00;
  double aks[neta];
  aks[0]=1.14; aks[1]=1.14; aks[2]=1.15; aks[3]=1.15; aks[4]=1.18;
  aks[5]=1.17; aks[6]=1.19; aks[7]=1.27; aks[8]=1.21; aks[9]=1.32;
  aks[10]=1.30; aks[11]=1.00; aks[12]=1.00; aks[13]=1.00;
  aks[14]=1.00; aks[15]=1.00; aks[16]=1.00; aks[17]=1.00; aks[18]=1.00;
  double acc[neta];
  acc[0]=0.040; acc[1]=0.043; acc[2]=0.044; acc[3]=0.057; acc[4]=0.053;
  acc[5]=0.044; acc[6]=0.044; acc[7]=0.050; acc[8]=0.043; acc[9]=0.035;
  acc[10]=0.035; acc[11]=0.04; acc[12]=0.04; acc[13]=0.04;
  acc[14]=0.035; acc[15]=0.04; acc[16]=0.04; acc[17]=0.04; acc[18]=0.04;
  */

  // Load data

  //TGraphAsymmErrors *gnd(0), *gnm(0);
  TGraphAsymmErrors *gnd = (TGraphAsymmErrors*)fn->Get("Data/RMS"); assert(gnd);
  TGraphAsymmErrors *gnm = (TGraphAsymmErrors*)fn->Get("MC/RMS"); assert(gnm);

  //TH1F *hdd = (TH1F*)fd->Get(Form("data_JER_standard_SM%d",ieta+1)); assert(hdd);
  //TH1F *hdm = (TH1F*)fd->Get(Form("MC_JER_standard_SM%d",ieta+1)); assert(hdm);
  TH1F *hdd(0), *hdm(0);
  int jeta = hxd->FindBin(hxc->GetBinCenter(ieta+1));
  /*
  if (ieta==0) {
    hdd = (TH1F*)fd->Get(Form("data_JER_standard_SM%d",ieta+1)); assert(hdd);
    hdm = (TH1F*)fd->Get(Form("MC_JER_standard_SM%d",ieta+1)); assert(hdm);
  }
  else {
    hdd = (TH1F*)fd->Get(Form("data_JER_standard_FE%d",ieta)); assert(hdd);
    hdm = (TH1F*)fd->Get(Form("MC_JER_standard_FE%d",ieta)); assert(hdm);
  }
  */
  if (jeta==1) {
    hdd = (TH1F*)fd->Get(Form("data_JER_standard_SM%d",jeta)); assert(hdd);
    hdm = (TH1F*)fd->Get(Form("MC_JER_standard_SM%d",jeta)); assert(hdm);
  }
  else {
    hdd = (TH1F*)fd->Get(Form("data_JER_standard_FE%d",jeta-1)); assert(hdd);
    hdm = (TH1F*)fd->Get(Form("MC_JER_standard_FE%d",jeta-1)); assert(hdm);
  }
  TGraphErrors *gdd = new TGraphErrors(hdd); cleanGraph(gdd);
  TGraphErrors *gdm = new TGraphErrors(hdm); cleanGraph(gdm);

  TH1D *hsd(0), *hsm(0);
  TGraphErrors *gsd(0), *gsm(0);
  if (ieta<5) {
    hsd = (TH1D*)fs->Get("zjet_jer_sn_data"); assert(hsd);
    hsm = (TH1D*)fs->Get("zjet_jer_sn_mc");   assert(hsm);
    gsd = new TGraphErrors(hsd); cleanGraph(gsd,15.,130.);
    gsm = new TGraphErrors(hsm); cleanGraph(gsm,15.,130.);
  }

  //TH1D *hcd = (TH1D*)fc->Get("jerc_rms_data"); assert(hcd); // Cterm
  //TH1D *hcm = (TH1D*)fc->Get("jerc_rms_mc"); assert(hcm); // Cterm
  TH1D *hcd = (TH1D*)fc->Get("jerc_rms_data_v3"); assert(hcd); // Cterm_v2
  TH1D *hcm = (TH1D*)fc->Get("jerc_rms_mc_v3"); assert(hcm); // Cterm_v2
  TGraphErrors *gcd = new TGraphErrors(hcd); //cleanGraph(gcd);
  TGraphErrors *gcm = new TGraphErrors(hcm); //cleanGraph(gcm);

  double eta1 = hcd->GetBinLowEdge(ieta+1);
  double eta2 = hcd->GetBinLowEdge(ieta+2);

  // Map eta-binned results to pT
  // For ptn: N-dominated region, from where pT*(1+2sigma)~15 GeV
  // Now JER~0.45 at 8 GeV, so 8*(1+2*0.45)=15.2
  double ptn = 8.;
  double ks = aks[ieta];//1.15;
  double cs = acs[ieta];//0.90;
  TGraphErrors *gndx = new TGraphErrors(1);
  gndx->SetPoint(0, ptn, tools::oplus(gnd->GetY()[ieta]/ptn, ks*cs/sqrt(ptn)));
  gndx->SetPointError(0, 0., gnd->GetEYhigh()[ieta]/ptn);
  TGraphErrors *gnmx = new TGraphErrors(1);
  gnmx->SetPoint(0, ptn, tools::oplus(gnm->GetY()[ieta]/ptn, cs/sqrt(ptn)));
  gnmx->SetPointError(0, 0., gnm->GetEYhigh()[ieta]/ptn);

  // For ptn: N-dominated region, from where pT*(1+2sigma)~15 GeV
  // Now JER~0.45 at 8 GeV, so 8*(1+2*0.45)=15.2
  double ptc = 3500.;
  double kc0 = 1.00;
  double cc = acc[ieta];
  TGraphErrors *gcdx = new TGraphErrors(1);
  gcdx->SetPoint(0, ptc, tools::oplus(gcd->GetY()[ieta]*0.01, kc0*cc));//2D(+)MC
  gcdx->SetPoint(0, ptc, tools::oplus(gcdx->GetY()[0], ks*cs/sqrt(ptc)));
  gcdx->SetPointError(0, 0., gcd->GetEY()[ieta]*0.01);
  TGraphErrors *gcmx = new TGraphErrors(1);
  gcmx->SetPoint(0, ptc, tools::oplus(gcm->GetY()[ieta]*0.01, cc)); // 2D(+)MC
  gcmx->SetPoint(0, ptc, tools::oplus(gcmx->GetY()[0], cs/sqrt(ptc)));
  gcmx->SetPointError(0, 0., gcm->GetEY()[ieta]*0.01);

  // Calculate ratios

  TH1F *hdr = (TH1F*)hdd->Clone(Form("hdr_%d",ieta));
  hdr->Divide(hdm);
  TGraphErrors *gdr = new TGraphErrors(hdr);
  TGraphErrors *gnrx = tools::ratioGraphs(gndx,gnmx);
  TGraphErrors *gcrx = tools::ratioGraphs(gcdx,gcmx);

  TH1D *hsr(0);
  TGraphErrors *gsr(0);
  if (hsd && hsm) {
    hsr = (TH1D*)hsd->Clone(Form("hsr_%d",ieta));
    hsr->Divide(hsm);
    gsr = new TGraphErrors(hsr); cleanGraph(gsr,15.,130.);
  }

  // Define resolution functions
  double xmin = 7;//6.5;
  double xmax = 4000;
  TF1 *fnd = new TF1(Form("fnd_%d",ieta),"[0]/x",xmin,xmax);
  fnd->SetParameter(0,gnd->GetY()[ieta]);
  TF1 *fnm = new TF1(Form("fnm_%d",ieta),"[0]/x",xmin,xmax);
  fnm->SetParameter(0,gnm->GetY()[ieta]);
  TF1 *fsd = new TF1(Form("fsd_%d",ieta),"[0]/sqrt(x)",xmin,xmax);
  fsd->SetParameter(0,ks*cs);
  TF1 *fsm = new TF1(Form("fsm_%d",ieta),"[0]/sqrt(x)",xmin,xmax);
  fsm->SetParameter(0,cs);
  TF1 *fcd = new TF1(Form("fcd_%d",ieta),"[0]",xmin,xmax);
  fcd->SetParameter(0,tools::oplus(gcd->GetY()[ieta]*0.01, kc0*cc));
  TF1 *fcm = new TF1(Form("fcm_%d",ieta),"[0]",xmin,xmax);
  fcm->SetParameter(0,tools::oplus(gcm->GetY()[ieta]*0.01, cc));

  TF1 *fjd = new TF1(Form("fjd_%d",ieta),
		     "sqrt([0]*[0]/(x*x)+[1]*[1]/x+[2]*[2])",xmin,xmax);
  fjd->SetParameters(fnd->GetParameter(0),fsd->GetParameter(0),
		     fcd->GetParameter(0));
  TF1 *fjm = new TF1(Form("fjm_%d",ieta),
		     "sqrt([0]*[0]/(x*x)+[1]*[1]/x+[2]*[2])",xmin,xmax);
  fjm->SetParameters(fnm->GetParameter(0),fsm->GetParameter(0),
		     fcm->GetParameter(0));
  TF1 *fjr = new TF1(Form("fjr_%d",ieta),
		     "sqrt([0]*[0]/(x*x)+[1]*[1]/x+[2]*[2])/"
		     "sqrt([3]*[3]/(x*x)+[4]*[4]/x+[5]*[5])",xmin,xmax);
  fjr->SetParameters(fnd->GetParameter(0),fsd->GetParameter(0),
		     fcd->GetParameter(0),
		     fnm->GetParameter(0),fsm->GetParameter(0),
		     fcm->GetParameter(0));
  // Get minimum and maximum in the physical range
  fjr->SetRange(ptn,min(ptc,6500./cosh(eta1)));
  double sfmax = fjr->GetMaximum();
  double sfmin = fjr->GetMinimum();
  double sfmid = 0.5*(sfmax+sfmin);
  double sfdif = 0.5*(sfmax-sfmin);
  fjr->SetRange(xmin,xmax);

  TH1F *hddr = (TH1F*)hdd->Clone(Form("hddr_%d",ieta)); hddr->Divide(fjd);
  TH1F *hdmr = (TH1F*)hdm->Clone(Form("hdmr_%d",ieta)); hdmr->Divide(fjm);

  TH1D *hsdr(0), *hsmr(0);
  if (hsd && hsm) {
    hsdr = (TH1D*)hsd->Clone(Form("hsdr_%d",ieta)); hsdr->Divide(fjd);
    hsmr = (TH1D*)hsm->Clone(Form("hsmr_%d",ieta)); hsmr->Divide(fjm);
  }

  TH1D *hu = tdrHist(Form("hu_%d",ieta),"JER",0.,0.65,"p_{T} (GeV)",xmin,xmax);
  TH1D *hd = tdrHist(Form("hd_%d",ieta),"Data/MC",0.95,1.25,"p_{T} (GeV)",xmin,xmax);
  if (ieta+1>=11) hd->GetYaxis()->SetRangeUser(0.85,1.75);//0.8,1.5);
  hd->GetXaxis()->SetMoreLogLabels(kFALSE);
  lumi_13TeV = "2018";
  extraText = "Private";
  TCanvas *c1 = tdrDiCanvas(Form("c1_%d",ieta),hu,hd,4,11);

  c1->cd(1);
  gPad->SetLogx();

  fnd->SetLineColor(kRed);
  fnd->Draw("SAME");
  fnm->SetLineColor(kRed);
  fnm->SetLineStyle(kDotted);
  fnm->Draw("SAME");

  fsd->SetLineColor(kGreen+2);
  fsd->Draw("SAME");
  fsm->SetLineColor(kGreen+2);
  fsm->SetLineStyle(kDotted);
  fsm->Draw("SAME");

  fcd->SetLineColor(kBlue);
  fcd->Draw("SAME");
  fcm->SetLineColor(kBlue);
  fcm->SetLineStyle(kDotted);
  fcm->Draw("SAME");

  fjd->SetLineColor(kBlack);
  fjd->Draw("SAME");
  fjm->SetLineColor(kBlack);
  fjm->SetLineStyle(kDotted);
  fjm->Draw("SAME");

  tdrDraw(gdd,"Pz",kFullCircle,kBlack);
  tdrDraw(gdm,"Pz",kOpenCircle,kBlack);

  tdrDraw(gndx,"Pz",kFullCircle,kRed);
  tdrDraw(gnmx,"Pz",kOpenCircle,kRed);

  if (gsd && gsm) {
    tdrDraw(gsd,"Pz",kFullDiamond,kGreen+2); gsd->SetMarkerSize(1.5);
    tdrDraw(gsm,"Pz",kOpenDiamond,kGreen+2); gsm->SetMarkerSize(1.5);
  }

  tdrDraw(gcdx,"Pz",kFullCircle,kBlue);
  tdrDraw(gcmx,"Pz",kOpenCircle,kBlue);

  TLegend *leg = tdrLeg(0.58,0.90-(gsd ? 5 : 4)*0.05,0.78,0.90);
  leg->AddEntry(gdd,"Dijet data","PLE");
  leg->AddEntry(gdm,"Dijet MC","PLE");
  leg->AddEntry(gndx,"RC noise #oplusS","PLE");
  if (gsd) leg->AddEntry(gsd,"Z+jet stoch. #oplusN","PLE");
  leg->AddEntry(gcdx,"2D const. #oplusC_{0}#oplusS","PLE");
  TLegend *leg1 = tdrLeg(0.32,0.90-10*0.05,0.52,0.90-5*0.05);
  leg1->AddEntry(fjd,Form("JER data (SF=[%1.3f,%1.3f])",sfmin,sfmax),"L");
  double n = fjm->GetParameter(0);
  double s = fjm->GetParameter(1);
  double c = fjm->GetParameter(2);
  leg1->AddEntry(fjm,Form("MC (N=%1.1f, S=%1.2f, C=%1.3f)",n,s,c),
		 //fjm->GetParameter(0),fjm->GetParameter(1),
		 //fjm->GetParameter(2)),
		 "L");
  double kn = fjd->GetParameter(0)/fjm->GetParameter(0);
  leg1->AddEntry(fnm,Form("Noise (SF=%1.3f)",kn),
		 //fjd->GetParameter(0)/fjm->GetParameter(0)),
		 "L");
  double ksb = fjd->GetParameter(1)/fjm->GetParameter(1);
  assert(ks==ksb);
  leg1->AddEntry(fsm,Form("Stochastic (SF=%1.3f)",ks),
		 //fjd->GetParameter(1)/fjm->GetParameter(1)),
		 "L");
  double kc = fjd->GetParameter(2)/fjm->GetParameter(2);
  leg1->AddEntry(fcm,Form("Constant (SF=%1.3f)",kc),
		 //fjd->GetParameter(2)/fjm->GetParameter(2)),
		 "L");
  TLegend *leg1b = tdrLeg(0.32,0.90-10*0.05+0.012,0.52,0.90-7*0.05+0.012);
  leg1b->AddEntry(fnd," ","L");
  leg1b->AddEntry(fsd," ","L");
  leg1b->AddEntry(fcd," ","L");

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  if (ieta==0) tex->DrawLatex(0.32,0.80,"|#eta|<0.261, SM1");
  else tex->DrawLatex(0.32,0.80,Form("%1.3f#leq|#eta|<%1.3f",eta1,eta2));

  c1->cd(2);
  gPad->SetLogx();

  // Calculate error band for JER SF using max/min
  TGraphErrors *ge = new TGraphErrors(100);
  TGraphErrors *ge2 = new TGraphErrors(100);
  for (int i = 0; i != ge->GetN(); ++i) {
    double x = xmin*pow(xmax/xmin,i/(ge->GetN()-1));
    ge->SetPoint(i, x, sfmid);//0.5*(sfmax+sfmin));
    ge->SetPointError(i, 0., sfdif);//0.5*(sfmax-sfmin));
    double xmax2 = min(ptc,6500./cosh(eta1));
    double xmin2 = 8.;
    double x2 = xmin2*pow(xmax2/xmin2,i/(ge->GetN()-1));
    ge2->SetPoint(i, x2, sfmid);//0.5*(sfmax+sfmin));
    ge2->SetPointError(i, 0., sfdif);//0.5*(sfmax-sfmin));
  }
  tdrDraw(ge,"E3",kNone,kYellow,kSolid,-1,1001,kYellow-9);
  tdrDraw(ge2,"E3",kNone,kYellow+1,kSolid,-1,1001,kYellow+1);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(xmin,1,xmax,1);
  l->SetLineStyle(kDotted);
  l->DrawLine(xmin,sfmax,xmax,sfmax);
  l->DrawLine(xmin,sfmin,xmax,sfmin);

  fjr->SetLineColor(kBlack);
  fjr->Draw("SAME");

 
  tdrDraw(gdr,"Pz",kFullCircle,kBlack);
  tdrDraw(gnrx,"Pz",kFullCircle,kRed);
  if (gsr) { tdrDraw(gsr,"Pz",kFullDiamond,kGreen+2); gsr->SetMarkerSize(1.5); }
  tdrDraw(gcrx,"Pz",kFullCircle,kBlue);


  tdrDraw(hddr,"Pz",kFullCircle,kGray+1); hddr->SetMarkerSize(0.5);
  tdrDraw(hdmr,"Pz",kOpenCircle,kGray+1); hdmr->SetMarkerSize(0.5);

  if (hsdr && hsmr) {
    tdrDraw(hsdr,"Pz",kFullDiamond,kGreen-9); hsdr->SetMarkerSize(0.75);
    tdrDraw(hsmr,"Pz",kOpenDiamond,kGreen-9); hsmr->SetMarkerSize(0.75);
  }
  
  TLegend *leg2 = tdrLeg(0.30,0.35,0.50,0.55);
  leg2->SetTextColor(kGray+1);
  leg2->AddEntry(hddr,"Dijet Data / JER","PLE");
  leg2->AddEntry(hdmr,"Dijet MC / JER","PLE");
  
  gPad->RedrawAxis();

  if (jeta==1)
    //c1->SaveAs(Form("pdf/jerSF/jerSF_SM%d.pdf",ieta+1));
    c1->SaveAs(Form("pdf/jerSF/jerSF_SM%d_ieta%d.pdf",jeta,ieta));
  else
    //c1->SaveAs(Form("pdf/jerSF/jerSF_FE%d.pdf",ieta+1));
    c1->SaveAs(Form("pdf/jerSF/jerSF_FE%d_ieta%d.pdf",jeta-1,ieta));


  // Summarize all results
  _c0->cd(ieta+1);
  gPad->SetLogx();
  if (ieta%5==0) {
    gPad->SetLeftMargin(0.25);
    hd->GetYaxis()->SetTitleOffset(0.9);
  }
  if (ieta/5==0) {
    gPad->SetTopMargin(0.10);
  }
  if (ieta/5==3) {
    gPad->SetBottomMargin(0.30);
    hd->GetXaxis()->SetTitleOffset(0.9);
  }

  hd->Draw();
  tdrDraw(ge,"E3",kNone,kYellow,kSolid,-1,1001,kYellow-9);
  tdrDraw(ge2,"E3",kNone,kYellow+1,kSolid,-1,1001,kYellow+1);

  l->SetLineStyle(kDashed);
  l->DrawLine(xmin,1,xmax,1);
  l->SetLineStyle(kDotted);
  l->DrawLine(xmin,sfmax,xmax,sfmax);
  l->DrawLine(xmin,sfmin,xmax,sfmin);

  fjr->Draw("SAME");
  tdrDraw(gdr,"Pz",kFullCircle,kBlack);
  tdrDraw(gnrx,"Pz",kFullCircle,kRed);
  if (gsr) { tdrDraw(gsr,"Pz",kFullDiamond,kGreen+2); gsr->SetMarkerSize(1.5); }
  tdrDraw(gcrx,"Pz",kFullCircle,kBlue);

  tex->SetTextSize(tex->GetTextSize()*2.0);
  tex->DrawLatex(ieta%5==0 ? 0.30 : 0.05, ieta/5==0 ? 0.70 : 0.80,
		 Form("%s%d_%d", jeta==1 ? "SM" : "FE",
		      jeta==1 ? jeta : jeta-1, ieta));
  double x1 = (ieta%5==0 ? 0.30 : 0.05);
  double y1 = (ieta/5==0 ? 0.80 : 0.90);
  tex->SetTextColor(kGreen+2);
  if (ieta==0) tex->DrawLatex(x1,y1,Form("|#eta|<%1.3f",
					 hcd->GetBinLowEdge(ieta+2)));
  else tex->DrawLatex(x1,y1,Form("%1.3f#leq|#eta|<%1.3f",
				 hcd->GetBinLowEdge(ieta+1),
				 hcd->GetBinLowEdge(ieta+2)));

  gPad->RedrawAxis();

  // Draw legend pad
  //if (ieta==18) {
  if (ieta==0) {
    _c0->cd(20);
    gPad->SetLogx();
    gPad->SetBottomMargin(0.30);
    hd->GetXaxis()->SetTitleOffset(0.9);
    hd->Draw();
    TLegend *leg = tdrLeg(0.05,0.95-0.10*(gsr ? 6: 5),0.55,0.95);
    leg->SetTextSize(leg->GetTextSize()*2.0);
    leg->AddEntry(gdr,"Dijet (S #oplus N #oplus C)","PLE");
    leg->AddEntry(gnrx,"RC noise #oplus S","PLE");
    if (gsr) leg->AddEntry(gsr,"Z+jet stoch. #oplus N","PLE");
    leg->AddEntry(gcrx,"2D const. #oplus C_{0} #oplus S","PLE");
    leg->AddEntry(fjr,"JER SF vs p_{T}","L");
    leg->AddEntry(ge2,"Range of JER SF","F");
    //leg->AddEntry(ge,"JER SF vs p_{T}","FL");
  }

  // Return detailed information for analysis vs eta
  //struct JER {
  //double n, s, c;
  //double kn, ks, kc;
  //double sfmin, sfmid, sfdif, sfmax;
  //double sf8, sf80, sf3500;
  //};
  JER jer;
  jer.eta = hcd->GetBinCenter(ieta+1);
  jer.deta = 0.5*hcd->GetBinWidth(ieta+1);
  jer.n  =  n; jer.s  =  s; jer.c  =  c;
  jer.kn = kn; jer.ks = ks; jer.kc = kc;
  jer.sfmin = sfmin; jer.sfmid = sfmid; jer.sfdif = sfdif; jer.sfmax = sfmax;
  jer.sf8    = fjr->Eval(8.);
  jer.sf80   = fjr->Eval(80.);
  jer.sf3500 = fjr->Eval(3500.);
  jer.sfxmin  = fjr->Eval(ptn);
  jer.sfxmax  = fjr->Eval(min(ptc,6500./cosh(eta1)));

  return jer;
} // jerSF
