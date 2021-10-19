// Purpose: estimate HCAL, ECAL and tracker variations with toyPF results
//          fit these for input into global fitter
#include "TFile.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TMultiGraph.h"

#include "../tdrstyle_mod15.C"

const double ptmin = 15;
const double ptmax = 4500;
const double kd = 0.22;//0.15;
const double km = 0.12;//0.22;//0.12;//0.075;

const bool trackerOnly = false;//true;
const bool useZjet = true; // else dijet

void varPlotsResp();
void varPlotsComp(string sf);
void cleanGraph(TGraphErrors *g);

TCanvas *_c1(0);
TLegend *_leg(0);
void varPlots() {

  // Summary plot of all composition changes overlaid
  setTDRStyle();
  TDirectory *curdir = gDirectory;
  TH1D *h = tdrHist("hCompAll","PF composition changes (10^{-2})",
		    //-1,1.5,"p_{T} (GeV)",ptmin,ptmax);
		    -1.2+1e-4,1.3-1e-4,"p_{T} (GeV)",ptmin,ptmax);
  lumi_13TeV = "toyPF";
  TCanvas *c1 = tdrCanvas("c1CompAll",h,4,11,kSquare);
  gPad->SetLogx();
  TLegend *leg = tdrLeg(0.60,0.69,0.90,0.87);
  TLine *l = new TLine(); l->SetLineStyle(kDashed);
  l->DrawLine(ptmin,0,ptmax,0);
  _c1 = c1;
  _leg = leg;
  curdir->cd();

  varPlotsResp();

  varPlotsComp("chf");
  varPlotsComp("gammaf");
  varPlotsComp("nhf");
  varPlotsComp("ef");
  //varPlotsComp("muf");
  
  c1->SaveAs(Form("pdf/varPlotsComp_all_%s.pdf",useZjet ? "zjet" : "dijet"));
}

void varPlotsResp() {
  setTDRStyle();

  TDirectory *curdir = gDirectory;  
  TFile *f3(0);
  //TFile *f3 = new TFile("rootfiles/varPlots_Mikael_5M_20200514.root","READ");
  if (useZjet) f3 = new TFile("rootfiles/varPlots_Mikael_5M_20200604_Zjet.root","READ"); // same as above 20200514 (should be)
  //TFile *f3 = new TFile("rootfiles/varPlots_Mikael_1M_20200604_dijet.root","READ");
  else f3 = new TFile("rootfiles/varPlots_Mikael_5M_20200615_dijet.root","READ"); 
  assert(f3 && !f3->IsZombie());

  // For data tracking efficiency
  TFile *fd1 = new TFile("rootfiles/output-DATA-2b-UL17V4_BCDEF_sub2p7.root","READ");
  assert(fd1 && !fd1->IsZombie());
  //TFile *fd2 = new TFile("rootfiles/output-DATA-2b-UL17V4_BCDEF.root","READ");
  //assert(fd2 && !fd2->IsZombie());
  TFile *fd3 = new TFile("rootfiles/output-DATA-2b-UL17V4_BCDEF_plus2p7.root","READ");
  assert(fd3 && !fd3->IsZombie());

  // For MC tracking efficiency
  TFile *fm1 = new TFile("rootfiles/output-MCNU-2b-UL17V4_BCDEF_sub2p7.root","READ");
  assert(fm1 && !fm1->IsZombie());
  TFile *fm3 = new TFile("rootfiles/output-MCNU-2b-UL17V4_BCDEF_plus2p7.root","READ");
  assert(fm3 && !fm3->IsZombie());

  // For Gen (MC truth) tracking efficiency
  TFile *fg1 = new TFile("rootfiles/output-MC-1-UL17V4_BCDEF_sub2p7.root","READ");
  assert(fg1 && !fg1->IsZombie());
  TFile *fg3 = new TFile("rootfiles/output-MC-1-UL17V4_BCDEF_plus2p7.root","READ");
  assert(fg3 && !fg3->IsZombie());

  // For MC MinBias cross section biases (80 mb vs 69.2 mb)
  TFile *fm69 = new TFile("rootfiles/output-MC-1-UL17V4_BCDEF.root","READ");
  assert(fm69 && !fm69->IsZombie());
  TFile *fm80 = new TFile("rootfiles/output-MC80-1-UL17V4_BCDEF.root","READ");
  assert(fm80 && !fm80->IsZombie());

  // ToyPF or fragmentation (Herwig7 vs Pythia8 M2)
  /*
  TFile *ftp = new TFile("rootfiles/CMSJES_P8_Zjet_5000000_290420.root","READ");
  assert(ftp && !ftp->IsZombie());
  TFile *fth = new TFile("rootfiles/CMSJES_H7_Zjet_1000000.root","READ");
  assert(fth && !fth->IsZombie());
  */
  TFile *ftpz = new TFile("rootfiles/P8H7DijetZjet_1M_080620/CMSJES_P8_Zjet_1000000.root","READ");
  assert(ftpz && !ftpz->IsZombie());
  TFile *fthz = new TFile("rootfiles/P8H7DijetZjet_1M_080620/CMSJES_H7_Zjet_1000000.root","READ");
  assert(fthz && !fthz->IsZombie());
  //
  TFile *ftpd = new TFile("rootfiles/P8H7DijetZjet_1M_080620/CMSJES_P8_dijet_1000000.root","READ");
  assert(ftpd && !ftpd->IsZombie());
  TFile *fthd = new TFile("rootfiles/P8H7DijetZjet_1M_080620/CMSJES_H7_dijet_1000000.root","READ");
  assert(fthd && !fthd->IsZombie());

  TFile *f8cp5 = new TFile("rootfiles/output-P8CP5-1-EOY17_DE.root","READ");
  assert(f8cp5 && !f8cp5->IsZombie());
  TFile *f8m1 = new TFile("rootfiles/output-P8M1-1-EOY17_DE.root","READ");
  assert(f8m1 && !f8m1->IsZombie());
  TFile *fhw = new TFile("rootfiles/output-HW-1-EOY17_DE.root","READ");
  assert(fhw && !fhw->IsZombie());

  curdir->cd();

  TH1D *hcp3 = (TH1D*)f3->Get("h_Rjet_Cp3"); assert(hcp3);
  TH1D *hcm3 = (TH1D*)f3->Get("h_Rjet_Cm3"); assert(hcm3);
  TH1D *hhp3 = (TH1D*)f3->Get("h_Rjet_HadHCALp3"); assert(hcp3);
  TH1D *hhm3 = (TH1D*)f3->Get("h_Rjet_HadHCALm3"); assert(hcm3);
  TH1D *hep3 = (TH1D*)f3->Get("h_Rjet_HadECALp3"); assert(hcp3);
  TH1D *hem3 = (TH1D*)f3->Get("h_Rjet_HadECALm3"); assert(hcm3);

  TH1D *ht3 = (TH1D*)f3->Get("h_Rjet_Trkm3"); assert(ht3);
  TH1D *hp = (TH1D*)f3->Get("h_Rjet_Photonm3"); assert(hp);
  curdir->cd();

  // Difference in RunBCDEF high and low tracking efficiency regions
  TProfile *pd1 = (TProfile*)fd1->Get("Standard/Eta_0.0-1.3/ppt_probepertag");
  assert(pd1);
  //TProfile *pd2 = (TProfile*)fd2->Get("Standard/Eta_0.0-1.3/ppt_probepertag");
  //assert(pd2);
  TProfile *pd3 = (TProfile*)fd3->Get("Standard/Eta_0.0-1.3/ppt_probepertag");
  assert(pd3);

  TProfile *pm1 = (TProfile*)fm1->Get("Standard/Eta_0.0-1.3/ppt_probepertag");
  assert(pm1);
  TProfile *pm3 = (TProfile*)fm3->Get("Standard/Eta_0.0-1.3/ppt_probepertag");
  assert(pm3);

  // options: p2r_guw (no weights), p2r_g (weights); guw smoother, consistent
  TProfile *pg1 = (TProfile*)fg1->Get("Standard/Eta_0.0-1.3/mc/p2r_guw");
  assert(pg1);
  TProfile *pg3 = (TProfile*)fg3->Get("Standard/Eta_0.0-1.3/mc/p2r_guw");
  assert(pg3);

  // p2r_guw does not show effect (PU weighed in) so have to use pr2_g
  TProfile *pm69 = (TProfile*)fm69->Get("Standard/Eta_0.0-1.3/mc/p2r_g");
  assert(pm69);
  TProfile *pm80 = (TProfile*)fm80->Get("Standard/Eta_0.0-1.3/mc/p2r_g");
  assert(pm80);

  // toyPF Herwig7 H7-UE-MMHT vs Pythia8 M2
  TProfile *ptpz = (TProfile*)ftpz->Get("prMPF");
  assert(ptpz);
  TProfile *pthz = (TProfile*)fthz->Get("prMPF");
  assert(pthz);
  TProfile *ptpd = (TProfile*)ftpd->Get("prMPF");
  assert(ptpd);
  TProfile *pthd = (TProfile*)fthd->Get("prMPF");
  assert(pthd);

  // options: p2r_guw and p2r_g (guw is smoother)
  TProfile *p8cp5 = (TProfile*)f8cp5->Get("Standard/Eta_0.0-1.3/mc/p2r_guw");
  assert(p8cp5);
  TProfile *p8m1 = (TProfile*)f8m1->Get("Standard/Eta_0.0-1.3/mc/p2r_guw");
  assert(p8m1);
  TProfile *phw = (TProfile*)fhw->Get("Standard/Eta_0.0-1.3/mc/p2r_guw");
  assert(phw);

  TH1D *hd = pd1->ProjectionX("hd_ppt_probpertag");
  for (int i = 1; i != hd->GetNbinsX()+1; ++i) {
    double y1 = pd1->GetBinContent(i); // phi<2.7
    double ey1 = pd1->GetBinError(i);
    //double y2 = pd2->GetBinContent(i); // all
    //double ey2 = pd2->GetBinError(i);
    double y3 = pd3->GetBinContent(i); // phi>2.7
    double ey3 = pd3->GetBinError(i);
    //double k = 0.15;//0.2;
    hd->SetBinContent(i, kd*100.*(y3-y1));
    hd->SetBinError(i, kd*100.*ey3);
  }

  TH1D *hm = pm1->ProjectionX("hm_ppt_probpertag");
  for (int i = 1; i != hm->GetNbinsX()+1; ++i) {
    double y1 = pm1->GetBinContent(i); // phi<2.7
    double ey1 = pm1->GetBinError(i);
    double y3 = pm3->GetBinContent(i); // phi>2.7
    double ey3 = pm3->GetBinError(i);
    //double k = 0.075;//0.15;//4*0.075;//0.15;//0.2;
    hm->SetBinContent(i, km*100.*(y3-y1));
    hm->SetBinError(i, km*100.*ey3);
  }

  TH1D *hg = pg1->ProjectionX("hg_p2r_guw");
  TGraphErrors *gg = new TGraphErrors(0);
  for (int i = 1; i != hm->GetNbinsX()+1; ++i) {
    double y1 = pg1->GetBinContent(i); // phi<2.7
    double ey1 = pg1->GetBinError(i);
    double y3 = pg3->GetBinContent(i); // phi>2.7
    double ey3 = pg3->GetBinError(i);
    hg->SetBinContent(i, km*100.*(y3-y1));
    hg->SetBinError(i, km*100.*ey3);
    // Shift truth result vs pT,tag
    double pt = pg1->GetBinCenter(i);
    double y = pm1->GetBinContent(pm1->FindBin(pt));
    if (y3!=0 && ey3 !=0 && y!=0) {
      int n = gg->GetN();
      gg->SetPoint(n, pt/y, km*100*(y3-y1));
      gg->SetPointError(n, 0.5*hg->GetBinWidth(i), km*100*ey3);
    }
  } // for i

  TH1D *hm80 = pm80->ProjectionX("hm80_p2r_g");
  for (int i = 1; i != hm80->GetNbinsX()+1; ++i) {
    double y80 = pm80->GetBinContent(i);
    double ey80 = pm80->GetBinError(i);
    double y69 = pm69->GetBinContent(i); // phi>2.7
    double ey69 = pm69->GetBinError(i);
    hm80->SetBinContent(i, 100.*(y69-y80));
    double cf = 0.5; // correlation factor
    hm80->SetBinError(i, cf*100.*ey80);
  }

  TH1D *htfz = ptpz->ProjectionX("htfz");
  for (int i = 1; i != htfz->GetNbinsX()+1; ++i) {
    double y1 = ptpz->GetBinContent(i);
    double ey1 = ptpz->GetBinError(i);
    assert(ptpz->GetBinLowEdge(i)==pthz->GetBinLowEdge(i));
    double y3 = pthz->GetBinContent(i);
    double ey3 = pthz->GetBinError(i);
    htfz->SetBinContent(i, 100.*(y3-y1));
    htfz->SetBinError(i, 100.*ey3);
  } // for i
  TH1D *htfd = ptpd->ProjectionX("htfz");
  for (int i = 1; i != htfd->GetNbinsX()+1; ++i) {
    double y1 = ptpd->GetBinContent(i);
    double ey1 = ptpd->GetBinError(i);
    assert(ptpd->GetBinLowEdge(i)==pthd->GetBinLowEdge(i));
    double y3 = pthd->GetBinContent(i);
    double ey3 = pthd->GetBinError(i);
    htfd->SetBinContent(i, 100.*(y3-y1));
    htfd->SetBinError(i, 100.*ey3);
  } // for i

  TH1D *hp8m1 = p8m1->ProjectionX("hm_p2r_g");
  TH1D *hw = phw->ProjectionX("hm_p2r_g");
  for (int i = 1; i != hd->GetNbinsX()+1; ++i) {
    double y5 = p8cp5->GetBinContent(i); // P8CP5
    double ey5 = p8cp5->GetBinError(i);
    double y1 = p8m1->GetBinContent(i); // P8M1
    double ey1 = p8m1->GetBinError(i);
    double y3 = phw->GetBinContent(i); // Herwig++
    double ey3 = phw->GetBinError(i);
    hp8m1->SetBinContent(i, 100.*(y1-y5));
    hp8m1->SetBinError(i, 100.*sqrt(ey1*ey1+ey5*ey5));
    hw->SetBinContent(i, 100.*(y3-y5));
    hw->SetBinError(i, 100.*sqrt(ey3*ey3+ey5*ey5));
  }


  ht3 = (TH1D*)ht3->Clone("ht3");
  hp = (TH1D*)hp->Clone("hp");

  TH1D *h1 = (TH1D*)hcp3->Clone("h1"); h1->Divide(h1);
  for (int i = 1; i != h1->GetNbinsX()+1; ++i) h1->SetBinError(i, 0);
  TH1D *h1p = (TH1D*)hp->Clone("h1p"); h1p->Divide(h1p);
  for (int i = 1; i != h1p->GetNbinsX()+1; ++i) h1p->SetBinError(i, 0);

  hcp3->Add(h1,-1);
  hcp3->Scale(100.);
  hcm3->Add(h1,-1);
  hcm3->Scale(100.);

  hhp3->Add(h1,-1);
  hhp3->Scale(100.);
  hhm3->Add(h1,-1);
  hhm3->Scale(100.);
  hep3->Add(h1,-1);
  hep3->Scale(100.);
  hem3->Add(h1,-1);
  hem3->Scale(100.);

  ht3->Add(h1p,-1);
  ht3->Scale(100.);
  hp->Add(h1p,-1);
  hp->Scale(100.);

  TH1D *hcmm3 = (TH1D*)hcm3->Clone("hcmm3");
  hcmm3->Scale(-1);

  // Set limit to 1.7% at 3 TeV and 0.3% at 15 GeV
  // (75% of hadrons, with 50% depositing 100% to HCAL, 50% depositing 50%)
  // (low pT has O(10%) HCAL energy)
  // (so 75%*75*3%=1.7%, 10%*3%=0.3%)
  int j15 = hcp3->FindBin(15.);
  hcp3->SetBinContent(j15, 0.25);

  // Log-ling interpolated SPR from -3% to +3%
  TH1D *hcx = (TH1D*)hcp3->Clone("hcx"); hcx->Reset();
  TH1D *hhx = (TH1D*)hhp3->Clone("hhx"); hhx->Reset();
  for (int i = 1; i != hcx->GetNbinsX()+1; ++i) {
    double pt = hcp3->GetBinCenter(i);
    double p = hcp3->GetBinContent(i);
    if (pt>2000) p = min(p, -hcm3->GetBinContent(i));
    // Log-lin interpolation from -1 at 15 GeV to +1 at 2884 GeV (0 at 208 GeV)
    double w = -1 + log(pt/15.)/log(208./15.);
    double x = log(pt/208.);
    hcx->SetBinContent(i, w*p);
    hcx->SetBinError(i, hcp3->GetBinError(i));
    //
    assert(hhx->GetBinLowEdge(i)==hcx->GetBinLowEdge(i));
    hhx->SetBinContent(i, w*hhp3->GetBinContent(i));
    hhx->SetBinError(i, hhp3->GetBinError(i));
  }

  /*
  // Calculate shapes shifted to zero at 208 GeV
  // should also probably scale their magnitude for better comparison
  int j208 = hcp3->FindBin(208.);
  double dc = hcp3->GetBinContent(j208)-hcx->GetBinContent(j208);
  double dt = ht->GetBinContent(j208);
  double dp = hp->GetBinContent(j208);
  TH1D *hcxp = (TH1D*)hcp3->Clone("hcxp");
  TH1D *hcp03 = (TH1D*)hcp3->Clone("hcp03");
  TH1D *ht0 = (TH1D*)ht->Clone("ht0");
  TH1D *hp0 = (TH1D*)hp->Clone("hp0");
  for (int i = 1; i != hcxp->GetNbinsX()+1; ++i) {
    hcxp->SetBinContent(i, hcx->GetBinContent(i)+dc);
    hcp03->SetBinContent(i, hcp03->GetBinContent(i)-dc);
    ht0->SetBinContent(i, ht0->GetBinContent(i)-dt);
    hp0->SetBinContent(i, hp0->GetBinContent(i)-dp);
  }
  */

  TH1D *h = tdrHist("h","Response change (%)",-2.,3.,
		    "p_{T} (GeV)",ptmin,ptmax);
  if (trackerOnly) { h->GetYaxis()->SetRangeUser(-1.2,1.3); }
  lumi_13TeV = "toyPF";
  TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);
  gPad->SetLogx();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(ptmin,0.,ptmax,0.);

  // Run I SPR for comparison
  TF1 *fhb = new TF1("fhb","100.*(max(0.,[0]+[1]*pow(x,[2]))-1)",10,ptmax);
  fhb->SetParameters(1.03091e+00, -5.11540e-02, -1.54227e-01); // SPRH
  fhb->SetLineColor(kMagenta+2);
  fhb->SetLineStyle(kDashed);
  fhb->SetLineWidth(3);
  //fhb->Draw("SAME");

  if (trackerOnly) {
    tdrDraw(ht3,"P",kFullSquare,kBlack); //ht3->SetMarkerSize(0.7);
    tdrDraw(hd,"P",kOpenSquare,kRed); //hd->SetMarkerSize(0.5);
    tdrDraw(hm,"P",kOpenSquare,kBlue); //hm->SetMarkerSize(0.5);
    tdrDraw(gg,"Pz",kFullDiamond,kGreen+2); gg->SetMarkerSize(1.3);//0.8);
  }
  else {
    //tdrDraw(hhp3,"P",kOpenTriangleUp,kRed);
    //tdrDraw(hhm3,"P",kOpenTriangleDown,kBlue);
    //tdrDraw(hhx,"P",kOpenDiamond,kBlack);
    tdrDraw(hhp3,"P",kFullCircle,kRed);
    tdrDraw(hhm3,"P",kFullCircle,kBlue);
    tdrDraw(hhx,"P",kFullCircle,kBlack);
    //tdrDraw(hep3,"P",kOpenTriangleDown,kBlue);//kRed);
    tdrDraw(hem3,"P",kOpenTriangleUp,kBlue);
    
    //tdrDraw(hcmm3,"P",kOpenCircle,kBlue);
    //tdrDraw(hcp3,"P",kFullCircle,kRed);
    //tdrDraw(hcm3,"P",kFullCircle,kBlue);
    //tdrDraw(hcx,"P",kFullCircle,kBlack);
    
    tdrDraw(ht3,"P",kFullSquare,kGreen+2); ht3->SetMarkerSize(0.7);
    tdrDraw(hp,"P",kFullDiamond,kCyan+2);
    
    
    //hd->GetXaxis()->SetRangeUser(max(49.,ptmin),ptmax);
    tdrDraw(hd,"P",kOpenSquare,kGreen+3); hd->SetMarkerSize(0.5);
    //hm->GetXaxis()->SetRangeUser(max(49.,ptmin),ptmax);
    tdrDraw(hm,"P",kOpenSquare,kGreen+4); hm->SetMarkerSize(0.5);
    //tdrDraw(hg,"P",kOpenDiamond,kGreen+4); hg->SetMarkerSize(0.8);
    tdrDraw(gg,"Pz",kFullDiamond,kGreen+4); gg->SetMarkerSize(0.8);
    //tdrDraw(hm1,"P",kOpenDiamond,kMagenta+1);
    tdrDraw(hw,"P",kFullDiamond,kMagenta+1);
    hw->GetXaxis()->SetRangeUser(max(25.,ptmin),min(4500.,ptmax));

    tdrDraw(htfz,"P",kFullDiamond,kMagenta+2); htfz->SetMarkerSize(1.5);
    tdrDraw(htfd,"P",kOpenDiamond,kMagenta+2); htfd->SetMarkerSize(1.5);
    tdrDraw(hm80,"Pz",kFullStar,kOrange+1);
  }

  TF1 *fcp3 = new TF1("fcp3","[0]+[1]*pow(x/[2],[3])/(1+pow(x/[2],[3]))*"
		     "(1-pow(x/[2],-[3]))",ptmin,ptmax);
  fcp3->SetParameters(+0.3,+1.5,208.,1);
  hcp3->Fit(fcp3,"RN");
  fcp3->SetLineColor(kRed);
  fcp3->SetLineWidth(3);
  fcp3->SetLineStyle(kDotted);
  //fcp3->Draw("SAME");

  TF1 *fcm3 = new TF1("fcm3","[0]+[1]*pow(x/[2],[3])/(1+pow(x/[2],[3]))*"
		     "(1-pow(x/[2],-[3]))",ptmin,ptmax);
  fcm3->SetParameters(-0.3,-1.5,208.,1);
  hcm3->Fit(fcm3,"QRN");
  fcm3->SetLineColor(kBlue);
  fcm3->SetLineStyle(kDotted);
  fcm3->SetLineWidth(3);
  //fcm3->Draw("SAME");

  TF1 *fhp3 = new TF1("fhp3","[0]+[1]*pow(x/[2],[3])/(1+pow(x/[2],[3]))*"
		      "(1-pow(x/[2],-[3]))",ptmin,ptmax);
  fhp3->SetParameters(+0.3,+1.5,208.,1);
  hhp3->Fit(fhp3,"QRN");
  fhp3->SetLineColor(kRed);
  fhp3->SetLineStyle(kDotted);//kDashed);
  fhp3->SetLineWidth(2);
  if (!trackerOnly) fhp3->Draw("SAME");

  TF1 *fhm3 = new TF1("fhm3","[0]+[1]*pow(x/[2],[3])/(1+pow(x/[2],[3]))*"
		      "(1-pow(x/[2],-[3]))",ptmin,ptmax);
  fhm3->SetParameters(-0.3,-1.5,208.,1);
  hhm3->Fit(fhm3,"QRN");
  fhm3->SetLineColor(kBlue);
  fhm3->SetLineStyle(kDotted);
  fhm3->SetLineWidth(2);
  if (!trackerOnly) fhm3->Draw("SAME");

  TF1 *fep3 = new TF1("fep3","[0]+[1]*pow(x/[2],[3])/(1+pow(x/[2],[3]))*"
		      "(1-pow(x/[2],-[3]))",ptmin,ptmax);
  fep3->SetParameters(-0.1,-0.5,208.,1);
  hep3->Fit(fep3,"QRN");
  fep3->SetLineColor(kBlue);
  fep3->SetLineStyle(kDashed);
  fep3->SetLineWidth(2);
  //fep3->Draw("SAME");

  TF1 *fem3 = new TF1("fem3","[0]+[1]*pow(x/[2],[3])/(1+pow(x/[2],[3]))*"
		      "(1-pow(x/[2],-[3]))",ptmin,ptmax);
  fem3->SetParameters(-0.1,-0.5,208.,1);
  hem3->Fit(fem3,"QRN");
  fem3->SetLineColor(kBlue);
  fem3->SetLineStyle(kDashed);
  fem3->SetLineWidth(2);
  if (!trackerOnly) fem3->Draw("SAME");

  TF1 *fcx = new TF1("fcx","[0]+[1]*pow(x/[2],[3])/(1+pow(x/[2],[3]))*"
		     "(1-pow(x/[2],-[3]))",ptmin,ptmax);
  fcx->SetParameters(-0.3,2.4,208.,1);
  hcx->Fit(fcx,"QRN");
  fcx->SetLineColor(kBlack);
  //fcx->Draw("SAME");

  TF1 *fhx = new TF1("fhx","[0]+[1]*pow(x/[2],[3])/(1+pow(x/[2],[3]))*"
		     "(1-pow(x/[2],-[3]))",ptmin,ptmax);
  fhx->SetParameters(-0.3,2.4,208.,1);
  hhx->Fit(fhx,"QRN");
  fhx->SetLineColor(kBlack);
  if (!trackerOnly) fhx->Draw("SAME");

  /*
  TF1 *fc = new TF1("fc","[0]+[1]*pow(x/[2],[3])/(1+pow(x/[2],[3]))*"
		    "(1-pow(x/[2],-[3]))",ptmin,ptmax);
  fc->SetParameters(fcp3->GetParameter(0)-fcp3->Eval(208.),
		    fcp3->GetParameter(1), fcp3->GetParameter(2),
		    fcp3->GetParameter(3));
  hcp03->Fit(fc,"QRN");
  fc->SetLineColor(kOrange+2);
  //fc->Draw("SAME");
  */

  TF1 *ft3 = new TF1("ft3","[0]+[1]*pow(x/208.,[2])",ptmin,ptmax);
  ft3->SetParameters(0,-0.15,-0.3);
  ht3->Fit(ft3,"QRN");
  ft3->SetLineColor(kGreen+2);
  ft3->Draw("SAME");

  // Limit fd fit range to pT>49 GeV as there is some discontinuity that is not
  // consistent with toyPF expectations. Something with Zero Bias trigger?
  //TF1 *fd = new TF1("fd","[0]+[1]*pow(x/208.,[2])+[3]/x",max(49.,ptmin),ptmax);
  TF1 *fd = new TF1("fd","[0]+[1]*pow(x/208.,[2])+[3]/x",ptmin,ptmax);
  fd->SetParameters(ft3->GetParameter(0),ft3->GetParameter(1),
		    ft3->GetParameter(2),0.5);
  fd->FixParameter(2,ft3->GetParameter(2));
  hd->Fit(fd,"QRN");
  fd->SetLineColor(kGreen+3);
  fd->SetLineWidth(2);
  fd->SetLineStyle(kDashed);
  fd->SetRange(ptmin,ptmax);
  fd->Draw("SAME");

  //TF1 *fm = new TF1("fm","[0]+[1]*pow(x/208.,[2])+[3]/x",max(49.,ptmin),ptmax);
  TF1 *fm = new TF1("fm","[0]+[1]*pow(x/208.,[2])+[3]/x",ptmin,ptmax);
  fm->SetParameters(ft3->GetParameter(0),ft3->GetParameter(1),
		    ft3->GetParameter(2),0.5);
  fm->FixParameter(2,ft3->GetParameter(2));
  hm->Fit(fm,"QRN");
  fm->SetLineColor(kGreen+4);
  fm->SetLineWidth(2);
  fm->SetLineStyle(kDashed);
  fm->SetRange(ptmin,ptmax);
  fm->Draw("SAME");

  // Combined fit to data and MC simultaneously
  // Same slope+exponent, different offsets
  TGraphErrors *gmg = new TGraphErrors(0);
  for (int i = hm->GetNbinsX(); i != -1; --i) {
    double pt = hm->GetBinCenter(i);
    if (pt>ptmin && pt<ptmax && hm->GetBinError(i)!=0) {
      int n = gmg->GetN();
      gmg->SetPoint(n, -pt, hm->GetBinContent(i));
      gmg->SetPointError(n, 0.25*hm->GetBinWidth(i), hm->GetBinError(i));
    }
  }
  for (int i = 1; i != hd->GetNbinsX()+1; ++i) {
    double pt = hd->GetBinCenter(i);
    if (pt>ptmin && pt<ptmax && hd->GetBinError(i)!=0) {
      int n = gmg->GetN();
      //gmg->SetPoint(n, pt, hd->GetBinContent(i));
      //gmg->SetPointError(n, 0.25*hd->GetBinWidth(i), hd->GetBinError(i));
    }
  }
  //tdrDraw(gmg,"Pz",kFullSquare,kBlack);
  /*
  TF1 *fmg = new TF1("fmg",
		     "(x>0)*([0]+[1]*pow(x/208.,[2])+[3]/x)"
		      "+(x<0)*([4]+[1]*pow(abs(x)/208.,[2])+[5]/abs(x))",
		     -ptmax,ptmax);
  fmg->SetParameters(ft3->GetParameter(0),ft3->GetParameter(1),
		     ft3->GetParameter(2),0.5,
		     ft3->GetParameter(0),0.5);
  //fmg->SetParameters(fd->GetParameter(0),fd->GetParameter(1),
  //		     fd->GetParameter(2),fd->GetParameter(3),
  //		     fm->GetParameter(0), fm->GetParameter(3));
  gmg->Fit(fmg,"QRN");
  fmg->SetLineColor(kGreen+3);
  fmg->SetLineWidth(2);
  fmg->SetLineStyle(kSolid);//kDotted);
  fmg->SetRange(ptmin,ptmax);
  //fmg->Draw("SAME");

  // Draw data side of the combined fit
  TF1 *fmgd = new TF1("fmgd","[0]+[1]*pow(x/208.,[2])+[3]/x",ptmin,ptmax);
  fmgd->SetParameters(fmg->GetParameter(0),fmg->GetParameter(1),
  		      fmg->GetParameter(2),fmg->GetParameter(3));
  fmgd->SetLineColor(kGreen+3);
  fmgd->SetLineWidth(2);
  fmgd->SetLineStyle(kSolid);
  fmgd->Draw("SAME");

  // Draw also MC side of the fit (pars, 4,1,2,5)
  TF1 *fmgm = new TF1("fmgm","[0]+[1]*pow(x/208.,[2])+[3]/x",ptmin,ptmax);
  //fmgm->SetParameters(fm->GetParameter(0),fm->GetParameter(1),
  //		      fm->GetParameter(2),fm->GetParameter(3));
  fmgm->SetParameters(fmg->GetParameter(4),fmg->GetParameter(1),
  		      fmg->GetParameter(2),fmg->GetParameter(5));
  fmgm->SetLineColor(kGreen+4);
  fmgm->SetLineWidth(2);
  fmgm->SetLineStyle(kSolid);
  fmgm->Draw("SAME");
  */

  TF1 *ftfz = new TF1("ftfz","[0]+[1]*max(log(x/[2]),0.)",ptmin,ptmax);
  ftfz->SetParameters(0.2,0.6,350);
  //ftfz->FixParameter(2,540);
  htfz->Fit(ftfz,"QRN");
  ftfz->SetLineColor(kMagenta+2);
  ftfz->SetLineWidth(2);
  ftfz->SetLineStyle(kDashed);
  if (!trackerOnly) ftfz->Draw("SAME");
  //
  TF1 *ftfd = new TF1("ftfd","[0]+[1]*max(log(x/[2]),0.)",ptmin,ptmax);
  ftfd->SetParameters(0.2,0.6,350);
  //ftfd->FixParameter(2,540);
  htfd->Fit(ftfd,"QRN");
  ftfd->SetLineColor(kMagenta+2);
  ftfd->SetLineWidth(2);
  ftfd->SetLineStyle(kDashed);
  if (!trackerOnly) ftfd->Draw("SAME");

  TF1 *fp8m1 = new TF1("fp8m1","[0]+[1]*pow(x/208.,[2])+[3]/x",25,ptmax);
  fp8m1->SetParameters(1.5,-0.5,+0.3,0);
  hp8m1->Fit(fp8m1,"QRN");
  fp8m1->SetLineColor(kMagenta+1);
  fp8m1->SetLineWidth(2);
  fp8m1->SetLineStyle(kDashed);
  //fp8m1->Draw("SAME");

  TF1 *f1m80 = new TF1("fm80","[0]+[1]*pow(x/208.,[2])+[3]*exp(-[4]*x)",
		       ptmin,ptmax);
  f1m80->SetParameters(0.1,-0.1,-0.3,0.,1.);
  //f1m80->FixParameter(2,ft3->GetParameter(2));
  hm80->Fit(f1m80,"QRN");
  f1m80->SetLineColor(kOrange+2);
  f1m80->SetLineWidth(2);
  f1m80->SetLineStyle(kDotted);
  f1m80->SetRange(ptmin,ptmax);
  if (!trackerOnly) f1m80->Draw("SAME");

  TF1 *fw = new TF1("fhw","[0]+[1]*(1+(pow(x/[2],[3])-1)/(pow(x/[2],[3])+1))+"
		    "[4]/x+[5]*log(x)/x",25,ptmax);
  fw->SetParameters(0.9526,-0.3883,1285,2.46,18.1,-2.062); // FullMC
  hw->Fit(fw,"QRN");
  fw->SetLineColor(kMagenta+1);
  fw->SetLineWidth(2);
  fw->SetLineStyle(kDashDotted);
  if (!trackerOnly) fw->Draw("SAME");

  TF1 *fp = new TF1("fp","[0]",ptmin,ptmax);
  fp->SetParameters(0);
  hp->Fit(fp,"QRN");
  fp->SetLineColor(kCyan+2);
  if (!trackerOnly) fp->Draw("SAME");


  if (trackerOnly) {
    //TLegend *leg1 = tdrLeg(0.38,0.675,0.68,0.855);
    TLegend *leg1 = tdrLeg(0.40,0.70,0.70,0.90);
    leg1->AddEntry(ht3,"ToyPF -3%","PL");
    leg1->AddEntry(hd,Form("Data #times%1.3f",kd),"PLE");
    leg1->AddEntry(hm,Form("MC reco #times%1.3f",km),"PLE");
    leg1->AddEntry(gg,Form("MC gen #times%1.3f",km),"PLE");
  }
  else {
    TLegend *leg1 = tdrLeg(0.38,0.63,0.68,0.90);
    leg1->SetTextSize(0.040);
    //leg1->AddEntry(hcp3,"SPR +3%","PL");
    leg1->AddEntry(hhp3,"HCAL Had. +3%","PL");
    leg1->AddEntry(hhx,"HCAL Had. #pm3%","PL");
    leg1->AddEntry(hhm3,"HCAL Had. -3%","PL");
    //leg1->AddEntry(hep3,"ECAL Had. +3%","PL");
    //leg1->AddEntry(hw,"HS1 vs CP5","PLE");
    leg1->AddEntry(hem3,"ECAL Had. -3%","PL");
    leg1->AddEntry(ht3,"Tracking -3%","PL");
    leg1->AddEntry(hp,"Photons -3%","PL");
    //leg1->AddEntry(hem3,"ECAL Had. -3%","PL");
    //leg1->AddEntry(hcm3,"SPR -3%","PL");

    TLegend *leg2 = tdrLeg(0.17,0.150,0.47,0.27);
    //TLegend *leg2 = tdrLeg(0.17,0.19,0.47,0.27);
    leg2->SetTextSize(0.040);
    //leg2->AddEntry(hcx,"log-lin -3% to +3%","PL");
    leg2->AddEntry(hw,"HS1 vs CP5","PLE");
    //leg2->AddEntry(hd,"RunBCDEF Trk#times0.15","PLE");
    leg2->AddEntry(hd,Form("RunBCDEF Trk#times%1.3f",kd),"PLE");
    //leg2->AddEntry(hm,"MC Trk#times0.15","PLE");
    //leg2->AddEntry(hm,"MC Trk#times0.075","PLE");
    leg2->AddEntry(hm,Form("MC Trk#times%1.3f",km),"PLE");
    //leg2->AddEntry(hm,"MC Trk#times0.30","PLE");
  }

  cout << endl;
  cout << "  // Fits from minitools/varPlots.C" << endl;

  cout << "  // SPR -3% to +3% cross variation (ECAL+HCAL)" << endl;
  cout << Form("  if (!fcx) fcx = new TF1(\"fcx\",\"%s\",15,4500);\n",
	       fcx->GetExpFormula().Data());
  cout << Form("  fcx->SetParameters(%1.4g, %1.4g, %1.4g, %1.4g); // toyPF\n\n",
	       fcx->GetParameter(0), fcx->GetParameter(1),
	       fcx->GetParameter(2), fcx->GetParameter(3));

  cout << "  // SPRH -3% to +3% cross variation (HCAL only)" << endl;
  cout << Form("  if (!fhx) fhx = new TF1(\"fhx\",\"%s\",15,4500);\n",
	       fhx->GetExpFormula().Data());
  cout << Form("  fhx->SetParameters(%1.4g, %1.4g, %1.4g, %1.4g); // toyPF\n\n",
	       fhx->GetParameter(0), fhx->GetParameter(1),
	       fhx->GetParameter(2), fhx->GetParameter(3));
  /*
  cout << "  // SPR +3% variation" << endl;
  cout << Form("  if (!fch) fch = new TF1(\"fch\",\"%s\",15,4500);\n",
	       fc->GetExpFormula().Data());
  cout << Form("  fch->SetParameters(%1.4g, %1.4g, %1.4g, %1.4g); // toyPF\n\n",
	       fc->GetParameter(0), fc->GetParameter(1),
	       fc->GetParameter(2), fc->GetParameter(3));
  */
  cout << "  // SPRH -3% variation" << endl;
  cout << Form("  if (!fhh) fhh = new TF1(\"fhh\",\"%s\",15,4500);\n",
	       fhm3->GetExpFormula().Data());
  cout << Form("  fhh->SetParameters(%1.4g, %1.4g, %1.4g, %1.4g); // toyPF\n\n",
	       fhm3->GetParameter(0), fhm3->GetParameter(1),
	       fhm3->GetParameter(2), fhm3->GetParameter(3));

  cout << "  // SPRE -3% variation" << endl;
  cout << Form("  if (!feh) feh = new TF1(\"feh\",\"%s\",15,4500);\n",
	       fem3->GetExpFormula().Data());
  cout << Form("  feh->SetParameters(%1.4g, %1.4g, %1.4g, %1.4g); // toyPF\n\n",
	       fem3->GetParameter(0), fem3->GetParameter(1),
	       fem3->GetParameter(2), fem3->GetParameter(3));

  cout << "  // Tracking -3% variation in toyPF" << endl;
  cout << Form("  if (!ft) ft = new TF1(\"ft\",\"%s\",15,4500);\n",
	       ft3->GetExpFormula().Data());
  cout << Form("  ft->SetParameters(%1.4g, %1.4g, %1.4g); // toyPF\n\n",
	       ft3->GetParameter(0),ft3->GetParameter(1),ft3->GetParameter(2));

  cout << "  // Tracking '-1%' variation in UL2017 data" << endl;
  cout << Form("  if (!ftd) ftd = new TF1(\"ftd\",\"%s\",15,4500);\n",
	       fd->GetExpFormula().Data());
  cout << Form("  ftd->SetParameters(%1.4g, %1.4g, %1.4g, %1.4g); // Data\n\n",
	       fd->GetParameter(0), fd->GetParameter(1), fd->GetParameter(2),
	       fd->GetParameter(3));

  cout << "  // Tracking '-1%' variation in UL2017 MC" << endl;
  cout << Form("  if (!ftm) ftm = new TF1(\"ftm\",\"%s\",15,4500);\n",
	       fm->GetExpFormula().Data());
  cout << Form("  ftm->SetParameters(%1.4g, %1.4g, %1.4g, %1.4g); // MC\n\n",
	       fm->GetParameter(0), fm->GetParameter(1), fm->GetParameter(2),
	       fm->GetParameter(3));

  cout << "  // Photon -3% variation" << endl;
  cout << Form("  if (!fp) fp = new TF1(\"fp\",\"%s\",15,4500);\n",
	       fp->GetExpFormula().Data());
  cout << Form("  fp->SetParameter(0,%1.4g); // toyPF\n\n",
	       fp->GetParameter(0));

  /*
  cout << "  // P8M1 vs P8CP5" << endl;
  cout << Form("  if (!fm) fm = new TF1(\"fm\",\"%s\",15,4500);\n",
	       fp8m1->GetExpFormula().Data());
  cout << Form("  fm->SetParameter(0,%1.4g); // FullMC\n\n",
	       f8m1->GetParameter(0));
  */

  cout << "  // sigmaMB 69.2 mb to 80 mb variation" << endl;
  cout << Form("  if (!fm80) fm80 = new TF1(\"fm80\",\"%s\",15,4500);\n",
	       f1m80->GetExpFormula().Data());
  cout << Form("  fm80->SetParameters(%1.4g,%1.4g,%1.4g, %1.4g,%1.4g); // toyPF\n\n",
	       f1m80->GetParameter(0),f1m80->GetParameter(1),
	       f1m80->GetParameter(2),
	       f1m80->GetParameter(3),f1m80->GetParameter(4));

  cout << "  // ToyPF H7 MMHT vs P8 M2 (Z+jet)" << endl;
  cout << Form("  if (!ftfz) ftfz = new TF1(\"ftfz\",\"%s\",15,4500);\n",
	       ftfz->GetExpFormula().Data());
  cout << Form("  ftfz->SetParameters(%1.4g, %1.4g, %1.4g); // toyPF\n\n",
	       ftfz->GetParameter(0), ftfz->GetParameter(1),
	       ftfz->GetParameter(2));
  //
  cout << "  // ToyPF H7 MMHT vs P8 M2 (Dijet)" << endl;
  cout << Form("  if (!ftfd) ftfd = new TF1(\"ftfd\",\"%s\",15,4500);\n",
	       ftfd->GetExpFormula().Data());
  cout << Form("  ftfd->SetParameters(%1.4g, %1.4g, %1.4g); // toyPF\n\n",
	       ftfd->GetParameter(0), ftfd->GetParameter(1),
	       ftfd->GetParameter(2));

  cout << "  // H++ vs P8CP5" << endl;
  cout << Form("  if (!fhw) fhw = new TF1(\"fhw\",\"%s\",15,4500);\n",
	       fw->GetExpFormula().Data());
  cout << Form("  fhw->SetParameters(%1.4g,%1.4g,%1.4g,%1.4g,%1.4g,%1.4g); // FullMC\n\n",
	       fw->GetParameter(0),fw->GetParameter(1),fw->GetParameter(2),
	       fw->GetParameter(3),fw->GetParameter(4),fw->GetParameter(5));

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  tex->DrawLatex(0.77,0.87,"|#eta| < 1.3");

  if (trackerOnly)
    c1->SaveAs(Form("pdf/varPlotsResp_trackerOnly_%s.pdf",
		    useZjet ? "zjet" : "dijet"));
  else
    c1->SaveAs(Form("pdf/varPlotsResp_%s.pdf",useZjet ? "zjet" : "dijet"));

  // Cumulative sum
  // Need to double check the factor 3 for vp[6] in global fit
  //double vp[7] = {1,1,1,1,1,0,0};
  double vp[7] = {1.682,1.779,1.516,1.322,-1.144,0.015,-1.122}; // posterior
  TH1D *ha = (TH1D*)hd->Clone("haResp"); ha->Reset();
  for (int i = 1; i != ha->GetNbinsX()+1; ++i) {
    double pt = ha->GetBinCenter(i);
    ha->SetBinContent(i, fd->Eval(pt)  * vp[0]
		      + fp->Eval(pt)   * vp[1] 
		      + fhx->Eval(pt)  * vp[2]
		      + fhm3->Eval(pt) * vp[3]
		      + fem3->Eval(pt) * vp[4]
		      + fw->Eval(pt)   * vp[5]
		      + 3*(fd->Eval(pt)-fm->Eval(pt)) * vp[6]);
    ha->SetBinError(i, (fd->Eval(pt)-fm->Eval(pt)) * 1);
  }

  TH1D *h2 = tdrHist("h2","Response change (%)",-2.,3.,
		     "p_{T} (GeV)",ptmin,ptmax);
  if (trackerOnly) { h2->GetYaxis()->SetRangeUser(-3.0,2.0); }//-1.2,1.3); }
  TCanvas *c2 = tdrCanvas("c2",h2,4,11,kSquare);
  gPad->SetLogx();

  l->SetLineStyle(kDashed);
  l->DrawLine(ptmin,0.,ptmax,0.);

  tdrDraw(ha,"PE3",kFullCircle);
  tdrDraw(ha,"PHIST",kFullCircle);

  tex->DrawLatex(0.30,0.75,Form("p_{0}(%5.2f), p_{1}(%5.2f), p_{2}(%5.2f)",
				vp[0],vp[1],vp[2]));
  tex->DrawLatex(0.30,0.70,Form("p_{3}(%5.2f), p_{4}(%5.2f), p_{5}(%5.2f)",
				vp[3],vp[4],vp[5]));
  tex->DrawLatex(0.30,0.65,Form("p_{7}(%5.2f)",
				vp[6]));

  gPad->RedrawAxis();
  c2->SaveAs(Form("pdf/varPlotsResp_sum_%s.pdf",useZjet ? "zjet" : "dijet"));

} // varPlotsResp

void scaleStat(TH1D *h, double k) {
  for (int i = 1; i != h->GetNbinsX()+1; ++i)
    h->SetBinError(i, h->GetBinError(i)*k);
} // scaleStat

//TGraphErrors *_g(0);
//Double_t _fc(Double_t *x, Double_t *p) {
//assert(_g);
//return _g->Eval(*x);
//}

void varPlotsComp(string sf) {

  const char *ctp = "tp";//"tp" for tag-and-probe or "" for direct match
  const char *ctp2 = "tp"; // for EOY17 MC compositions
  setTDRStyle();

  TDirectory *curdir = gDirectory;  
  TFile *f(0);
  //TFile *f = new TFile("rootfiles/varPlots_Mikael_5M_20200514.root","READ");
  if (useZjet) f = new TFile("rootfiles/varPlots_Mikael_5M_20200604_Zjet.root","READ"); // same as above 20200514 (should be)
  //TFile *f = new TFile("rootfiles/varPlots_Mikael_1M_20200604_Dijet.root","READ");
  else f = new TFile("rootfiles/varPlots_Mikael_5M_20200615_dijet.root","READ");
  assert(f && !f->IsZombie());

  TFile *fd1 = new TFile("rootfiles/output-DATA-2b-UL17V4_BCDEF_sub2p7.root","READ");
  assert(fd1 && !fd1->IsZombie());
  //TFile *fd2 = new TFile("rootfiles/output-DATA-2b-UL17V4_BCDEF.root","READ");
  //assert(fd2 && !fd2->IsZombie());
  TFile *fd3 = new TFile("rootfiles/output-DATA-2b-UL17V4_BCDEF_plus2p7.root","READ");
  assert(fd3 && !fd3->IsZombie());
  
  TFile *fm1 = new TFile("rootfiles/output-MCNU-2b-UL17V4_BCDEF_sub2p7.root","READ");
  TFile *fm3 = new TFile("rootfiles/output-MCNU-2b-UL17V4_BCDEF_plus2p7.root","READ");

  // For MC MinBias cross section biases (80 mb vs 69.2 mb)
  TFile *fm69 = new TFile("rootfiles/output-MCNU-2b-UL17V4_BCDEF.root","READ");
  assert(fm69 && !fm69->IsZombie());
  TFile *fm80 = new TFile("rootfiles/output-MC80NU-2b-UL17V4_BCDEF.root","READ");
  assert(fm80 && !fm80->IsZombie());

  TFile *ftpz = new TFile("rootfiles/P8H7DijetZjet_1M_080620/CMSJES_P8_Zjet_1000000.root","READ");
  assert(ftpz && !ftpz->IsZombie());
  TFile *fthz = new TFile("rootfiles/P8H7DijetZjet_1M_080620/CMSJES_H7_Zjet_1000000.root","READ");
  assert(fthz && !fthz->IsZombie());
  //
  TFile *ftpd = new TFile("rootfiles/P8H7DijetZjet_1M_080620/CMSJES_P8_dijet_1000000.root","READ");
  assert(ftpd && !ftpd->IsZombie());
  TFile *fthd = new TFile("rootfiles/P8H7DijetZjet_1M_080620/CMSJES_H7_dijet_1000000.root","READ");
  assert(fthd && !fthd->IsZombie());

  TFile *fp8cp5 = new TFile("rootfiles/output-P8CP5-2b-EOY17_DE.root","READ");
  assert(fp8cp5 && !fp8cp5->IsZombie());
  TFile *fp8m1 = new TFile("rootfiles/output-P8M1-2b-EOY17_DE.root","READ");
  assert(fp8m1 && !fp8m1->IsZombie());
  TFile *fhw = new TFile("rootfiles/output-HW-2b-EOY17_DE.root","READ");
  assert(fhw && !fhw->IsZombie());

  curdir->cd();

  //string sf = "nhf";
  //string sf2 = (sf=="gammaf" ? "nef" : sf);
  string sf2 = sf;
  if (sf=="gammaf") sf2 = "nef";
  if (sf=="ef") sf2 = "cef";
  //if (sf=="muf") { sf2 = "cef"; } // PATCH
  const char *cf = sf.c_str();
  const char *cf2 = sf2.c_str();
  TH1D *hcp = (TH1D*)f->Get(Form("h%s_Cp3",cf)); assert(hcp);
  TH1D *hcm = (TH1D*)f->Get(Form("h%s_Cm3",cf)); assert(hcm);
  TH1D *hhp = (TH1D*)f->Get(Form("h%s_HadHCALp3",cf)); assert(hhp);
  TH1D *hhm = (TH1D*)f->Get(Form("h%s_HadHCALm3",cf)); assert(hhm);
  TH1D *hep = (TH1D*)f->Get(Form("h%s_HadECALp3",cf)); assert(hep);
  TH1D *hem = (TH1D*)f->Get(Form("h%s_HadECALm3",cf)); assert(hem);
  TH1D *ht1 = (TH1D*)f->Get(Form("h%s_Trkm1",cf)); assert(ht1);
  TH1D *ht3 = (TH1D*)f->Get(Form("h%s_Trkm3",cf)); assert(ht3);
  //TH1D *hp = (TH1D*)f->Get(Form("h%s_ECALm3",cf)); assert(hp); // was wrong!
  TH1D *hp = (TH1D*)f->Get(Form("h%s_Photonm3",cf)); assert(hp); // 20200604

  // Difference in RunF high and low tracking efficiency regions
  TProfile *pd1 = (TProfile*)fd1->Get(Form("Standard/Eta_0.0-1.3/p%s%s",cf2,ctp));
  assert(pd1);
  //TProfile *pd2 = (TProfile*)fd2->Get(Form("Standard/Eta_0.0-1.3/p%s%s",cf2,ctp));
  //assert(pd2);
  TProfile *pd3 = (TProfile*)fd3->Get(Form("Standard/Eta_0.0-1.3/p%s%s",cf2,ctp));
  assert(pd3);

  TProfile *pm1 = (TProfile*)fm1->Get(Form("Standard/Eta_0.0-1.3/p%s%s",cf2,ctp));
  assert(pm1);
  TProfile *pm3 = (TProfile*)fm3->Get(Form("Standard/Eta_0.0-1.3/p%s%s",cf2,ctp));
  assert(pm3);

  TProfile *pm69 = (TProfile*)fm69->Get(Form("Standard/Eta_0.0-1.3/p%s%s",cf2,ctp));
  assert(pm69);
  TProfile *pm80 = (TProfile*)fm80->Get(Form("Standard/Eta_0.0-1.3/p%s%s",cf2,ctp));
  assert(pm80);

  // toyPF Herwig7 H7-UE-MMHT vs Pythia8 M2
  TProfile *ptpz = (TProfile*)ftpz->Get(Form("pr%s",cf));
  assert(ptpz);
  TProfile *pthz = (TProfile*)fthz->Get(Form("pr%s",cf));
  assert(pthz);
  TProfile *ptpd = (TProfile*)ftpd->Get(Form("pr%s",cf));
  assert(ptpd);
  TProfile *pthd = (TProfile*)fthd->Get(Form("pr%s",cf));
  assert(pthd);

  TProfile *p8cp5 = (TProfile*)fp8cp5->Get(Form("Standard/Eta_0.0-1.3/p%s%s",cf2,ctp2));
  assert(p8cp5);
  TProfile *p8m1 = (TProfile*)fp8m1->Get(Form("Standard/Eta_0.0-1.3/p%s%s",cf2,ctp2));
  assert(p8m1);
  TProfile *phw = (TProfile*)fhw->Get(Form("Standard/Eta_0.0-1.3/p%s%s",cf2,ctp2));
  assert(phw);

  TH1D *hd = pd1->ProjectionX(Form("hd_%s",cf));
  for (int i = 1; i != hd->GetNbinsX()+1; ++i) {
    double y1 = pd1->GetBinContent(i); // phi<2.7
    double ey1 = pd1->GetBinError(i);
    //double y2 = pd2->GetBinContent(i); // all
    //double ey2 = pd2->GetBinError(i);
    double y3 = pd3->GetBinContent(i); // phi>2.7
    double ey3 = pd3->GetBinError(i);
    //double k = 0.15;//0.2;
    hd->SetBinContent(i, kd*100.*(y3-y1));
    hd->SetBinError(i, kd*100.*ey3);
  }

  TH1D *hm = pm1->ProjectionX(Form("hm_%s",cf));
  for (int i = 1; i != hd->GetNbinsX()+1; ++i) {
    double y1 = pm1->GetBinContent(i); // phi<2.7
    double ey1 = pm1->GetBinError(i);
    double y3 = pm3->GetBinContent(i); // phi>2.7
    double ey3 = pm3->GetBinError(i);
    //double k = 0.075;//0.15;//0.2;
    hm->SetBinContent(i, km*100.*(y3-y1));
    hm->SetBinError(i, km*100.*ey3);
  }

  TH1D *htfz = ptpz->ProjectionX("htfz");
  for (int i = 1; i != htfz->GetNbinsX()+1; ++i) {
    double y1 = ptpz->GetBinContent(i);
    double ey1 = ptpz->GetBinError(i);
    assert(ptpz->GetBinLowEdge(i)==pthz->GetBinLowEdge(i));
    double y3 = pthz->GetBinContent(i);
    double ey3 = pthz->GetBinError(i);
    htfz->SetBinContent(i, 100.*(y3-y1));
    htfz->SetBinError(i, 100.*ey3);
  } // for i
  TH1D *htfd = ptpd->ProjectionX("htfz");
  for (int i = 1; i != htfd->GetNbinsX()+1; ++i) {
    double y1 = ptpd->GetBinContent(i);
    double ey1 = ptpd->GetBinError(i);
    assert(ptpd->GetBinLowEdge(i)==pthd->GetBinLowEdge(i));
    double y3 = pthd->GetBinContent(i);
    double ey3 = pthd->GetBinError(i);
    htfd->SetBinContent(i, 100.*(y3-y1));
    htfd->SetBinError(i, 100.*ey3);
  } // for i

  TH1D *hm80 = pm80->ProjectionX(Form("hm80_%s",cf));
  for (int i = 1; i != hd->GetNbinsX()+1; ++i) {
    double y80 = pm80->GetBinContent(i); // sigmaMB=80 mb
    double ey80 = pm80->GetBinError(i);
    double y69 = pm69->GetBinContent(i); // sigmaMB=69.2 mb
    double ey69 = pm69->GetBinError(i);
    hm80->SetBinContent(i, 100.*(y69-y80));
    double cf = 0.5; // correlation factor to reduce stat uncertainty
    hm80->SetBinError(i, cf*100.*ey80);
  }

  // Statistical combination of data and MC
  TH1D *hmg = (TH1D*)hm->Clone(Form("hmg_%s",cf));
  hmg->Reset();
  for (int i = 1; i != hmg->GetNbinsX()+1; ++i) {
    double y1 = hd->GetBinContent(i);
    double ey1 = hd->GetBinError(i);
    double y2 = hm->GetBinContent(i);
    double ey2 = hm->GetBinError(i);
    double w1(0), w2(0);
    if (ey1!=0 && ey2!=0) {
      w1 = ey2*ey2/(ey1*ey1+ey2*ey2);
      w2 = ey1*ey1/(ey1*ey1+ey2*ey2);
    }
    else if (ey1!=0) { w1=1; w2=0; }
    else if (ey2!=0) { w1=0; w2=1; }
    double y = w1*y1+w2*y2;
    double ey = sqrt(w1*w1*ey1*ey1+w2*w2*ey2*ey2);
    hmg->SetBinContent(i, y);
    hmg->SetBinError(i, ey);
  }

  TH1D *h8m1 = p8m1->ProjectionX(Form("hm_%s",cf));
  TH1D *hw = phw->ProjectionX(Form("hw_%s",cf));
  for (int i = 1; i != hm->GetNbinsX()+1; ++i) {
    double y5 = p8cp5->GetBinContent(i); // P8CP5
    double ey5 = p8cp5->GetBinError(i);
    double y1 = p8m1->GetBinContent(i); // P8M1
    double ey1 = p8m1->GetBinError(i);
    double y3 = phw->GetBinContent(i); // H++
    double ey3 = phw->GetBinError(i);
    h8m1->SetBinContent(i, 100.*(y1-y5));
    h8m1->SetBinError(i, 100.*sqrt(ey1*ey1+ey5*ey5));
    hw->SetBinContent(i, 100.*(y3-y5));
    hw->SetBinError(i, 100.*sqrt(ey3*ey3+ey5*ey5));
  }

  // For 20200514 to turn to percentage
  hcp->Scale(100);
  hcm->Scale(100);
  hhp->Scale(100);
  hhm->Scale(100);
  hep->Scale(100);
  hem->Scale(100);
  ht1->Scale(100);
  ht3->Scale(100);
  hp->Scale(100);

  scaleStat(hcp,0.1);
  scaleStat(hcm,0.1);
  scaleStat(hhp,0.1);
  scaleStat(hhm,0.1);
  scaleStat(hep,0.1);
  scaleStat(hem,0.1);
  scaleStat(ht1,0.1);
  scaleStat(ht3,0.1);
  scaleStat(hp,0.1);

  // Interpolate SPR
  TH1D *hcx = (TH1D*)hcp->Clone(Form("hcx_%s",cf2)); hcx->Reset();
  TH1D *hhx = (TH1D*)hhp->Clone(Form("hhx_%s",cf2)); hhx->Reset();
  for (int i = 1; i != hcx->GetNbinsX()+1; ++i) {
    double pt = hcp->GetBinCenter(i);
    double p = hcp->GetBinContent(i);
    double w = -1 + log(pt/15.)/log(208./15.);
    double x = log(pt/208.);
    hcx->SetBinContent(i, w*p);
    hcx->SetBinError(i, hcp->GetBinError(i));
    //
    assert(hhx->GetBinLowEdge(i)==hcx->GetBinLowEdge(i));
    hhx->SetBinContent(i, w*hhp->GetBinContent(i));
    hhx->SetBinError(i, hhp->GetBinError(i));
  }

  TH1D *h = tdrHist(Form("hComp_%s",cf2),
		    Form("Composition change for %s (10^{-2})",cf2),
		    -0.7,1.3,"p_{T} (GeV)",ptmin,ptmax);
  if (trackerOnly) { h->GetYaxis()->SetRangeUser(-1.2,1.3); }
  if (string(ctp)=="")   h->SetXTitle("p_{T,probe} (GeV)");
  if (string(ctp)=="tp") h->SetXTitle("p_{T,tag} (GeV)");
  lumi_13TeV = "toyPF";
  TCanvas *c1 = tdrCanvas(Form("c1Comp_%s",cf2),h,4,11,kSquare);
  gPad->SetLogx();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(ptmin,0.,ptmax,0.);

  if (trackerOnly) {
    tdrDraw(ht1,"Pz",kFullSquare,kBlack);    //ht1->SetMarkerSize(0.7);
    tdrDraw(hd,"Pz",kOpenSquare,kRed); //hd->SetMarkerSize(0.5);
    tdrDraw(hm,"Pz",kOpenSquare,kBlue); //hm->SetMarkerSize(0.5);
    tdrDraw(hmg,"Pz",kFullSquare,kGreen+2); hmg->SetMarkerSize(0.8);//0.5);
  }
  else {
    //tdrDraw(hcp,"Pz",kFullCircle,kRed);
    //tdrDraw(hcm,"Pz",kFullCircle,kBlue);
    //tdrDraw(hcx,"Pz",kFullCircle,kBlack);
    //tdrDraw(hhp,"Pz",kFullCircle,kRed);
    tdrDraw(hhx,"Pz",kFullDiamond,kBlack);
    tdrDraw(hhm,"Pz",kFullCircle,kBlue);
    //tdrDraw(hh,"Pz",kOpenTriangleUp,kRed);   hh->SetMarkerSize(0.7);
    tdrDraw(hem,"Pz",kOpenTriangleUp,kBlue);  hem->SetMarkerSize(0.7);
    tdrDraw(ht1,"Pz",kFullSquare,kGreen+2);    ht1->SetMarkerSize(0.7);
    tdrDraw(hp,"Pz",kFullDiamond,kCyan+2); 
    
    tdrDraw(hd,"Pz",kOpenSquare,kGreen+3); hd->SetMarkerSize(0.5);
    tdrDraw(hm,"Pz",kOpenSquare,kGreen+4); hm->SetMarkerSize(0.5);
    tdrDraw(hmg,"Pz",kFullSquare,kGreen+4); hmg->SetMarkerSize(0.5);
    //tdrDraw(h8m1,"Pz",kOpenDiamond,kMagenta+1);
    tdrDraw(hw,"Pz",kFullDiamond,kMagenta+1); //hw->SetMarkerSize(0.7);

    tdrDraw(htfz,"P",kFullDiamond,kMagenta+2); htfz->SetMarkerSize(1.5);
    tdrDraw(htfd,"P",kOpenDiamond,kMagenta+2); htfd->SetMarkerSize(1.5);
    tdrDraw(hm80,"Pz",kFullStar,kOrange+1);

    hw->GetXaxis()->SetRangeUser(max(25.,ptmin),min(2000.,ptmax));
    
    //tdrDraw(ha,"Pz",kFullDiamond,kOrange+2);
  }

  if (trackerOnly) {
    TLegend *leg = tdrLeg(0.40,0.70,0.70,0.90);
    leg->AddEntry(ht1,"ToyPF -1%","PLE");
    leg->AddEntry(hd,Form("Data #times%1.3f",kd),"PLE");
    leg->AddEntry(hm,Form("MC reco #times%1.3f",km),"PLE");
    leg->AddEntry(hmg,"Data+MC combo","PLE");
  }
  else {
    //TLegend *leg = tdrLeg(0.62,0.615,0.90,0.90);
    TLegend *leg = tdrLeg(0.62,0.66,0.90,0.90);
    //leg->AddEntry(hhp,"HCAL +3%","PLE");
    leg->AddEntry(hhx,"HCAL #pm3%","PLE");
    leg->AddEntry(hhm,"HCAL -3%","PLE");
    leg->AddEntry(hem,"ECAL-h -3%","PLE");
    leg->AddEntry(ht1,"Tracking -1%","PLE");
    leg->AddEntry(hp,"Photons -3%","PLE");
    
    //TLegend *leg2 = tdrLeg(0.20,0.625,0.48,0.795);
    TLegend *leg2 = tdrLeg(0.20,0.650,0.48,0.795);
    leg2->AddEntry(hw,"HS1 vs CP5","PLE");
    //leg2->AddEntry(hd,"BCDEF Trk#times0.15","PLE");
    leg2->AddEntry(hd,Form("BCDEF Trk#times%1.3f",kd),"PLE");
    //leg2->AddEntry(hm,"MC Trk#times0.15","PLE");
    //leg2->AddEntry(hm,"MC Trk#times0.075","PLE");
    leg2->AddEntry(hm,Form("MC Trk#times%1.3f",km),"PLE");
  }

  // Fit variations
  TF1 *ft3 = new TF1(Form("ft3_%s",cf2),
		     "[0]+[1]*(1+(pow(x/[2],[3])-1)/(pow(x/[2],[3])+1))"
		     "+[4]*pow(x,-0.3)",
		     ptmin,ptmax);
  ft3->SetLineColor(kBlue);
  if (sf2=="chf") ft3->SetParameters(-1.5,+0.2,1000.,0.);//1.3,1,0.);
  if (sf2=="nhf") ft3->SetParameters(1.0,-0.3,1000.,1.3);
  if (sf2=="nef") ft3->SetParameters(0.4,+0.1,1000.,1.3);
  ht3->Fit(ft3,"QRN");
  //ft3->Draw("SAME");

  TF1 *ft1 = new TF1(Form("ft3_%s",cf2),
		     "[0]+[1]*(1+(pow(x/[2],[3])-1)/(pow(x/[2],[3])+1))"
		     "+[4]*pow(x,-0.3)",
		     ptmin,ptmax);
  ft1->SetParameters(ft3->GetParameter(0)/3.,ft3->GetParameter(1)/3.,
		     ft3->GetParameter(2),ft3->GetParameter(3),
		     ft3->GetParameter(4)/3.);
  ft1->SetLineColor(kGreen+2);
  ft1->Draw("SAME");

  TF1 *fd = new TF1(Form("fd_%s",cf2),
		     "[0]+[1]*(1+(pow(x/[2],[3])-1)/(pow(x/[2],[3])+1))"
		     "+[4]*pow(x,+0.3)+[5]/x",
		     ptmin,ptmax);
  fd->SetParameters(ft3->GetParameter(0)/3.,ft3->GetParameter(1)/3.,
  		    ft3->GetParameter(2),ft3->GetParameter(3),
  		    ft3->GetParameter(4)/3.,0);
  hd->Fit(fd,"QRN");
  fd->SetLineColor(kGreen+3);
  fd->SetLineWidth(2);
  fd->SetLineStyle(kDotted);//kDashed);
  fd->Draw("SAME");

  TF1 *fm = new TF1(Form("fm_%s",cf2),
		     "[0]+[1]*(1+(pow(x/[2],[3])-1)/(pow(x/[2],[3])+1))"
		     "+[4]*pow(x,+0.3)+[5]/x",
		     ptmin,ptmax);
  fm->SetParameters(fd->GetParameter(0),fd->GetParameter(1),
  		    fd->GetParameter(2),fd->GetParameter(3),
  		    fd->GetParameter(4),fd->GetParameter(5));
  hm->Fit(fm,"QRN");
  fm->SetLineColor(kGreen+4);
  fm->SetLineWidth(2);
  fm->SetLineStyle(kDotted);//kDashed);
  fm->Draw("SAME");

  // Fit data and MC together, as their shapes are very similar
  TGraphErrors *gd = new TGraphErrors(hd);
  TGraphErrors *gm = new TGraphErrors(hm);
  cleanGraph(gd);
  cleanGraph(gm);
  TMultiGraph *tmg = new TMultiGraph();
  tmg->Add(gd);
  tmg->Add(gm);
  TF1 *ftmg = new TF1(Form("ftmg_%s",cf2),
		      "[0]+[1]*(1+(pow(x/[2],[3])-1)/(pow(x/[2],[3])+1))"
		      //"+[4]*pow(x,+0.3)+[5]/x",
		      "+[4]*pow(x,[5])+[6]/x",
		      ptmin,ptmax);
  ftmg->SetParameters(fm->GetParameter(0),fm->GetParameter(1),
		      fm->GetParameter(2),fm->GetParameter(3),
		      fm->GetParameter(4),+0.3,fm->GetParameter(5));
  if (sf2=="chf") {
    ftmg->FixParameter(5,+0.3);
  }
  if (sf2=="nhf") {
    ftmg->SetParameters(0, -0.2,4000,1, 0.04, 0.3, 0);
    ftmg->FixParameter(5,+0.3);
    ftmg->FixParameter(2,4000);
  }
  if (sf2=="nef") {
    //ftmg->SetParameters(-3.304,4.75,2727,0.08363,-0.1305,-1.28); // 236.6/115
    //ftmg->SetParameters(0.1, -0.2,2000,1, 0.05,-0.2, -0.8);
    //ftmg->SetParameters(0.2257,-0.168,4000,1,0.01264,+0.3,-3.072);
    //ftmg->SetParameters(0.007747,-2.772,5000,0.9828,0.0315,0.5413,-0.7262);
    //ftmg->FixParameter(2,5000);
    ftmg->SetParameters(0.1143, -3.469, 5000, 1.026, 0.01281, 0.6582, -1.613); // 167.0/115
    ftmg->FixParameter(5,-0.3);
    //ftmg->FixParameter(2,4000);//3000);//165.0; 4000);//166.3 //5000); //167.0
  }
  tmg->Fit(ftmg,"QRN");
  ftmg->SetLineColor(kGreen+4);
  ftmg->SetLineWidth(2);
  ftmg->SetLineStyle(kSolid);
  ftmg->Draw("SAME");

  TF1 *f8m1 = new TF1(Form("f8m1_%s",cf2),
		    "[0]+[1]*pow(x,[2])+[3]/x",
		     ptmin,ptmax);
  if (sf2=="chf") f8m1->SetParameters(-0.2,0.,0.3,0);
  if (sf2=="nhf") f8m1->SetParameters(0.14,4.,0.3,0);
  if (sf2=="nef") f8m1->SetParameters(0,-5,0.3,0);
  h8m1->Fit(f8m1,"QRN");
  f8m1->SetLineColor(kMagenta+1);
  f8m1->SetLineWidth(2);
  f8m1->SetLineStyle(kDashed);
  //f8m1->Draw("SAME");

  TF1 *f1m80 = new TF1(Form("fm80_%s",cf2),
		       "[0]+[1]*pow(x,[2])",
		       max(49.,ptmin),ptmax);
  f1m80->SetParameters(0.2,-0.1,-0.3);
  hm80->Fit(f1m80,"QRN");
  f1m80->SetLineColor(kOrange+2);
  f1m80->SetLineWidth(2);
  f1m80->SetLineStyle(kDotted);
  f1m80->SetRange(ptmin,ptmax);
  f1m80->Draw("SAME");

  TF1 *fw = new TF1(Form("fhw_%s",cf2),
		    "[0]+[1]*pow(x,[2])+[3]/x",
		    //"[0]+[1]*(1+(pow(x/[2],[3])-1)/(pow(x/[2],[3])+1))"
		    //"+[4]*pow(x,[5])",
		    max(25.,ptmin),min(2000.,ptmax));
  if (sf2=="chf") {
    fw->SetParameters(-0.2, 0.,1,0);
    fw->FixParameter(3,0); 
  }
  if (sf2=="nhf") fw->SetParameters(0.14,4.,0.3,0);
  if (sf2=="nef") fw->SetParameters(0,-5,0.3,0);
  hw->Fit(fw,"QRN");
  fw->SetLineColor(kMagenta+1);
  fw->SetLineWidth(3);
  fw->SetLineStyle(kDashed);//kDashDotted);
  if (!trackerOnly) fw->Draw("SAME");

  TF1 *fp = new TF1(Form("fp_%s",cf2),
		    "[0]+[1]*(1+(pow(x/[2],[3])-1)/(pow(x/[2],[3])+1))"
		    "+[4]*pow(x,[5])",
		    ptmin,ptmax);
  fp->SetLineColor(kCyan+2);
  if (sf2!="chf") fp->FixParameter(1,0);
  if (sf2=="chf") fp->SetParameters(0.4,0.0,1000,1.3,-0.1,0.3);
  if (sf2=="nhf") fp->SetParameters(0.1,0.,1000.,1.3,+0.3,0.3);
  if (sf2=="nef") fp->SetParameters(-0.2,0.0,1000.,1.3,-0.3);
  hp->Fit(fp,"QRN");
  if (!trackerOnly) fp->Draw("SAME");

  TF1 *fcx = new TF1(Form("fcx_%s",cf2),
		     "[0]+[1]*(1+(pow(x/[2],[3])-1)/(pow(x/[2],[3])+1))"
		     "+[4]*pow(x,[5])",
		     ptmin,ptmax);
  fcx->SetLineColor(kBlack);
  if (sf2=="chf") fcx->SetParameters(-0.1,0,400.,1.3,2,-0.3);
  if (sf2=="nhf") fcx->SetParameters(-0.1,0,400.,1.3,2,+0.3);
  if (sf2=="nef") fcx->SetParameters(0.0,0,1000.,1.3,-2,+0.3);
  hcx->Fit(fcx,"QRN");
  //fcx->Draw("SAME");

  TF1 *fhx = new TF1(Form("fhx_%s",cf2),
		     "[0]+[1]*(1+(pow(x/[2],[3])-1)/(pow(x/[2],[3])+1))"
		     "+[4]*pow(x,[5])",
		     ptmin,ptmax);
  fhx->SetLineColor(kBlack);
  if (sf2=="chf") fhx->SetParameters(-0.1,0,400.,1.3,2,-0.3);
  if (sf2=="nhf") fhx->SetParameters(-0.1,0,400.,1.3,2,+0.3);
  if (sf2=="nef") fhx->SetParameters(0.0,0,1000.,1.3,-2,+0.3);
  hhx->Fit(fhx,"QRN");
  if (!trackerOnly) fhx->Draw("SAME");

  TF1 *fhp = new TF1(Form("fhp_%s",cf2),
		     "[0]+[1]*(1+(pow(x/[2],[3])-1)/(pow(x/[2],[3])+1))"
		     "+[4]*pow(x,[5])",
		     ptmin,ptmax);
  fhp->SetLineColor(kRed);
  if (sf2=="chf") fhp->SetParameters(-0.2,0.1,1000.,1.3,0.,0.3);
  if (sf2=="nhf") fhp->SetParameters(0.14,0,1000,1.3,4.,0.3);
  if (sf2=="nef") fhp->SetParameters(0,0,1000,1.3,-5,0.3);
  hhp->Fit(fhp,"QRN");
  //fhp->Draw("SAME");

  TF1 *fhm = new TF1(Form("fhm_%s",cf2),
		     "[0]+[1]*(1+(pow(x/[2],[3])-1)/(pow(x/[2],[3])+1))"
		     "+[4]*pow(x,[5])",
		     ptmin,ptmax);
  fhm->SetLineColor(kBlue);
  if (sf2=="chf") fhm->SetParameters(+0.2,-0.1,1000.,1.3,0.,0.3);
  if (sf2=="nhf") fhm->SetParameters(-0.14,0,1000,1.3,+4.,0.3);
  if (sf2=="nef") fhm->SetParameters(0,0,1000,1.3,+5,0.3);
  hhm->Fit(fhm,"QRN");
  if (!trackerOnly) fhm->Draw("SAME");

  TF1 *fep = new TF1(Form("fep_%s",cf2),
		     "[0]+[1]*(1+(pow(x/[2],[3])-1)/(pow(x/[2],[3])+1))"
		     "+[4]*pow(x,[5])",
		     ptmin,ptmax);
  fep->SetLineColor(kBlue);
  fep->FixParameter(1,0);
  if (sf2=="chf") fep->SetParameters(0.05,0.,1000.,1.3,-0.5,0.3);
  if (sf2=="nhf") fep->SetParameters(0.0,0,1000,1.3,4.,0.3);
  if (sf2=="nef") fep->SetParameters(0,0,1000,1.3,-4.,0.3);
  hep->Fit(fep,"QRN");
  //fep->Draw("SAME");

  TF1 *fem = new TF1(Form("fem_%s",cf2),
		     "[0]+[1]*(1+(pow(x/[2],[3])-1)/(pow(x/[2],[3])+1))"
		     "+[4]*pow(x,[5])",
		     ptmin,ptmax);
  fem->SetLineColor(kBlue);
  fem->FixParameter(1,0);
  if (sf2=="chf") fem->SetParameters(-0.05,0.,1000.,1.3,+0.5,0.3);
  if (sf2=="nhf") fem->SetParameters(0.0,0,1000,1.3,-4.,0.3);
  if (sf2=="nef") fem->SetParameters(0,0,1000,1.3,+4.,0.3);
  hem->Fit(fem,"QRN");
  if (!trackerOnly) fem->Draw("SAME");

  vector<TF1*> vf1;
  vf1.push_back(ft3);
  vf1.push_back(fp);
  vf1.push_back(fcx);
  vf1.push_back(fhx);
  vf1.push_back(fhm);
  vf1.push_back(fem);
  //vf1.push_back(fd);
  //vf1.push_back(fm);
  vf1.push_back(ftmg);
  vf1.push_back(f1m80);
  vf1.push_back(fw);

  cout << endl;
  cout << " // Fits from minitools/varPlots.C" << endl;
  for (unsigned int i = 0; i != vf1.size(); ++i) {
    TF1 *f1 = vf1[i];
    cout << Form("  TF1 *%s = new TF1(\"%s\",\"%s\",%1.0f,%1.0f);\n",
		 f1->GetName(),f1->GetName(),f1->GetExpFormula().Data(),
		 ptmin,ptmax);
    cout << Form("  %s->SetParameters(",f1->GetName());
    for (int j = 0; j != f1->GetNpar()-1; ++j)
      cout << Form("%1.4g, ",f1->GetParameter(j));
    cout << Form("%1.4g); // %1.1f/%d\n",
		 f1->GetParameter(f1->GetNpar()-1),
		 f1->GetChisquare(), f1->GetNDF());
  }
  cout << endl;


  if (trackerOnly)
    c1->SaveAs(Form("pdf/varPlotsComp_trackerOnly_%s_%s.pdf",cf2,
		    useZjet ? "zjet" : "dijet"));
  else
    c1->SaveAs(Form("pdf/varPlotsComp_%s_%s.pdf",cf2,
		    useZjet ? "zjet" : "dijet"));

  assert(_c1);
  assert(_leg);
  _c1->cd();

  map<string,int> color;
  color["chf"] = kRed;
  color["nhf"] = kGreen+2;
  color["nef"] = kBlue;
  color["cef"] = kCyan+1;

  map<string,int> marker;
  marker["chf"] = kFullCircle;
  marker["nhf"] = kFullDiamond;
  marker["nef"] = kFullSquare;
  marker["cef"] = kFullDiamond;

  //TH1D *ha2 = (TH1D*)ht1->Clone(Form("ha2_%s",cf2));
  //ha2->Add(ha2,hp,1,1.0);
  //ha2->Add(ha2,hcx,1,1.5);
  //ha2->Add(ha2,hhm,1,-2.0);
  //ha2->Add(ha2,hem,1,1.0);
  //
  TH1D *ha2 = (TH1D*)hd->Clone(Form("ha2_%s",cf2));
  TH1D *ha3 = (TH1D*)hd->Clone(Form("ha3_%s",cf2));
  ha2->Reset();
  ha3->Reset();
  // {Tr, Gam, HX, HH, HE, HW, TD
  //double vp[6] = {1.807, 1.699, 1.476, 1.366, -1.277, -0.206};
  //double vp[6] = {1.738, 1.669, 1.480, 1.401, -1.423, -0.094};
  //double vp[6] = {1, 1., 1., 1., 1, 0}; // prior
  double vp[6] = {1.682, 1.779, 1.516, 1.322, -1.144, 0.015}; // posterior
  for (int i = 1; i != hd->GetNbinsX()+1; ++i) {
    double pt = ha2->GetBinCenter(i);
    // NB: hd should be mg or statistical average of data and MC
    ha2->SetBinContent(i, 0
		       //+ hd ->GetBinContent( hd->FindBin(pt)) * vp[0]
		       + hmg->GetBinContent( hd->FindBin(pt)) * vp[0]
		       + hp ->GetBinContent( hp->FindBin(pt)) * vp[1]
		       + hhx->GetBinContent(hhx->FindBin(pt)) * vp[2]
		       + hhm->GetBinContent(hhm->FindBin(pt)) * vp[3]
		       + hem->GetBinContent(hem->FindBin(pt)) * vp[4] 
		       + hw ->GetBinContent( hw->FindBin(pt)) * vp[5]); 
    ha2->SetBinError(i, 0
		     //+ hd ->GetBinError( hd->FindBin(pt)) * vp[0]
		     + hmg->GetBinError( hd->FindBin(pt)) * vp[0]
		     + hp ->GetBinError( hp->FindBin(pt)) * vp[1]
		     + hhx->GetBinError(hhx->FindBin(pt)) * vp[2]
		     + hhm->GetBinError(hhm->FindBin(pt)) * vp[3]
		     + hem->GetBinError(hem->FindBin(pt)) * vp[4] 
		     + hw ->GetBinError( hw->FindBin(pt)) * vp[5]); 
    ha3->SetBinContent(i, 0
		       + ftmg->Eval(pt) * vp[0]
		       + fp  ->Eval(pt) * vp[1]
		       + fhx ->Eval(pt) * vp[2]
		       + fhm ->Eval(pt) * vp[3]
		       + fem ->Eval(pt) * vp[4] 
		       + fw  ->Eval(pt) * vp[5]); 
    ha3->SetBinError(i, 0.5*ha2->GetBinError(i));
  }
  /*
  ha2->Add(hd,  1.807);
  ha2->Add(hp,  1.699);
  ha2->Add(hhx, 1.476);
  ha2->Add(hhm, 1.366);
  ha2->Add(hem,-1.227);
  ha2->Add(hw, -0.206);
  */

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.035);
  //tex->DrawLatex(0.20,0.17,"P(-1.0%),CX(+4.5%),H(-4.5%),E(-3.0%)");

  if (sf2=="chf") {
    tex->DrawLatex(0.30,0.22,Form("p_{0}(%5.2f), p_{1}(%5.2f), p_{2}(%5.2f)",
				  vp[0],vp[1],vp[2]));
    tex->DrawLatex(0.30,0.17,Form("p_{3}(%5.2f), p_{4}(%5.2f), p_{7}(%5.2f)",
				  vp[3],vp[4],vp[5]));
  }

  ha2->SetMarkerSize(1.0);
  tdrDraw(ha3,"E3",marker[sf2],color[sf2],kSolid,-1,1001,color[sf2]-9);
  ha3->SetFillColorAlpha(color[sf2]-9,0.7);
  tdrDraw(ha2,"Pz",marker[sf2],color[sf2],kSolid);
  _leg->AddEntry(ha2,cf2,"PLE");
} // varPlotsComp

void cleanGraph(TGraphErrors *g) {
  for (int i = g->GetN()-1; i != -1; --i) {
    if (g->GetEY()[i]==0) g->RemovePoint(i);
  }
} // cleanGraph
