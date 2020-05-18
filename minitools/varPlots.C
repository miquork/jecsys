// Purpose: estimate HCAL, ECAL and tracker variations with toyPF results
//          fit these for input into global fitter
#include "TFile.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TProfile.h"

#include "../tdrstyle_mod15.C"

const double ptmin = 15;
const double ptmax = 4500;

void varPlotsResp();
void varPlotsComp(string sf);

TCanvas *_c1(0);
TLegend *_leg(0);
void varPlots() {

  // Summary plot of all composition changes overlaid
  setTDRStyle();
  TDirectory *curdir = gDirectory;
  TH1D *h = tdrHist("hCompAll","PF composition changes (10^{-2})",
		    -1,1.5,"p_{T} (GeV)",ptmin,ptmax);
  lumi_13TeV = "toyPF";
  TCanvas *c1 = tdrCanvas("c1CompAll",h,4,11,kSquare);
  gPad->SetLogx();
  TLegend *leg = tdrLeg(0.60,0.69,0.90,0.87);
  TLine *l = new TLine(); l->SetLineStyle(kDashed);
  l->DrawLine(ptmin,0,ptmax,0);
  _c1 = c1;
  _leg = leg;
  curdir->cd();

  //varPlotsResp();
  varPlotsComp("nhf");
  varPlotsComp("gammaf");
  varPlotsComp("chf");

  c1->SaveAs("pdf/varPlotsComp_all.pdf");
}

void varPlotsResp() {
  setTDRStyle();

  TDirectory *curdir = gDirectory;  
  //TFile *f = new TFile("rootfiles/varPlots_Mikael_5M_20200203.root","READ");
  //assert(f && !f->IsZombie());
  //TFile *f2 = new TFile("rootfiles/varPlots_Mikael_5M_20200425.root","READ");
  //assert(f2 && !f2->IsZombie());
  //TFile *f3 = new TFile("rootfiles/varHCAL_Mikael_5M_20200429.root","READ");
  //TFile *f3 = new TFile("rootfiles/varPlots_Mikael_5M_20200511.root","READ");
  TFile *f3 = new TFile("rootfiles/varPlots_Mikael_5M_20200514.root","READ");
  assert(f3 && !f3->IsZombie());
  curdir->cd();

  //TH1D *hcp3 = (TH1D*)f2->Get("h_Rjet_Cp3"); assert(hcp3);
  //TH1D *hcm3 = (TH1D*)f2->Get("h_Rjet_Cm3"); assert(hcm3);
  //TH1D *hcp3 = (TH1D*)f3->Get("h_Rjet_HCALp3"); assert(hcp3);
  //TH1D *hcm3 = (TH1D*)f3->Get("h_Rjet_HCALm3"); assert(hcm3);
  TH1D *hcp3 = (TH1D*)f3->Get("h_Rjet_Cp3"); assert(hcp3);
  TH1D *hcm3 = (TH1D*)f3->Get("h_Rjet_Cm3"); assert(hcm3);
  TH1D *hhp3 = (TH1D*)f3->Get("h_Rjet_HadHCALp3"); assert(hcp3);
  TH1D *hhm3 = (TH1D*)f3->Get("h_Rjet_HadHCALm3"); assert(hcm3);
  TH1D *hep3 = (TH1D*)f3->Get("h_Rjet_HadECALp3"); assert(hcp3);
  TH1D *hem3 = (TH1D*)f3->Get("h_Rjet_HadECALm3"); assert(hcm3);

  //TH1D *ht = (TH1D*)f->Get("h_Rjet_Trk"); assert(ht);
  //TH1D *hp = (TH1D*)f->Get("h_Rjet_Photon"); assert(hp);
  TH1D *ht = (TH1D*)f3->Get("h_Rjet_Trkm3"); assert(ht);
  //TH1D *hp = (TH1D*)f3->Get("h_Rjet_ECALm3"); assert(hp); // 20200511
  TH1D *hp = (TH1D*)f3->Get("h_Rjet_Photonm3"); assert(hp); // 20200514
  curdir->cd();

  ht = (TH1D*)ht->Clone("ht");
  hp = (TH1D*)hp->Clone("hp");
  TH1D *hpe = (TH1D*)hp->Clone("hpe");

  TH1D *h1 = (TH1D*)hcp3->Clone("h1"); h1->Divide(h1);
  for (int i = 1; i != h1->GetNbinsX()+1; ++i) h1->SetBinError(i, 0);
  TH1D *h1p = (TH1D*)hp->Clone("h1p"); h1p->Divide(h1p);
  for (int i = 1; i != h1p->GetNbinsX()+1; ++i) h1p->SetBinError(i, 0);

  hcp3->Add(h1,-1);
  hcp3->Scale(100.);
  //hcp3->Scale(sqrt(2.));
  hcm3->Add(h1,-1);
  hcm3->Scale(100.);
  //hcm3->Scale(sqrt(2.));

  hhp3->Add(h1,-1);
  hhp3->Scale(100.);
  hhm3->Add(h1,-1);
  hhm3->Scale(100.);
  hep3->Add(h1,-1);
  hep3->Scale(100.);
  hem3->Add(h1,-1);
  hem3->Scale(100.);

  ht->Add(h1p,-1);
  ht->Scale(100.);
  hp->Add(h1p,-1);
  hp->Scale(100.);
  //hp->Add(hem3,-1); // Photon only, 20200511
  hpe->Add(h1p,-1);
  hpe->Scale(100.); // Photon+SPRE

  TH1D *hcmm3 = (TH1D*)hcm3->Clone("hcmm3");
  hcmm3->Scale(-1);

  // Set limit to 1.7% at 3 TeV and 0.3% at 15 GeV
  // (75% of hadrons, with 50% depositing 100% to HCAL, 50% depositing 50%)
  // (low pT has O(10%) HCAL energy)
  // (so 75%*75*3%=1.7%, 10%*3%=0.3%)
  int j15 = hcp3->FindBin(15.);
  hcp3->SetBinContent(j15, 0.25);
  //hcp3->SetBinContent(j15, 0.3);
  /*
  int j3000 = hcp3->FindBin(3000.);
  hcp3->SetBinContent(j3000, 2.25);
  */

  TH1D *hcx = (TH1D*)hcp3->Clone("hcx"); hcx->Reset();
  TH1D *hcx2 = (TH1D*)hcp3->Clone("hcx2"); hcx2->Reset();
  for (int i = 1; i != hcx->GetNbinsX()+1; ++i) {
    double pt = hcp3->GetBinCenter(i);
    double p = hcp3->GetBinContent(i);
    if (pt>2000) p = min(p, -hcm3->GetBinContent(i));
    //if (pt>3000) p = 1.8;
    //double m = hcm3->GetBinContent(i); // = -p
    // Log-lin interpolation from -1 at 15 GeV to +1 at 2884 GeV (0 at 208 GeV)
    //double w = max(-1.,min(+1., -1 + log(pt/15.)/log(208./15.)));
    double w = -1 + log(pt/15.)/log(208./15.);
    double x = log(pt/208.);
    // 1.34 is exponent of x from fit of fc3
    // it also produces limits of +/-1 for x of -3 (10 GeV) and +3 (4.2 TeV)
    double w2 = 1.34*x/(1+fabs(x));
    //double w2 = 1.433*x/(1+fabs(x));
    hcx->SetBinContent(i, w*p);
    hcx->SetBinError(i, hcp3->GetBinError(i));
    //if (pt>3000) hcx->SetBinError(i, 0.2);
    hcx2->SetBinContent(i, w2*p);
    hcx2->SetBinError(i, hcp3->GetBinError(i));
  }

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
  //TH1D *h = new TH1D("h",";p_{T} (GeV);Change (%)",2485,15,ptmax);
  //h->GetXaxis()->SetMoreLogLabels();
  //h->GetXaxis()->SetNoExponent();

  TH1D *h = tdrHist("h","Response change (%)",-2.,3.,
		    "p_{T} (GeV)",ptmin,ptmax);
  lumi_13TeV = "toyPF";
  TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);
  gPad->SetLogx();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(ptmin,0.,ptmax,0.);

  //TF1 *fhb = new TF1("fhb","max(0.,[0]+[1]*pow(x,[2]))",10,ptmax);
  //TF1 *fhb = new TF1("fhb","100.*(max(1.,[0]+[1]*pow(x,[2]))-1)",10,ptmax);
  // The max(0.,) should have been max(1,), but the former was used up to now...
  TF1 *fhb = new TF1("fhb","100.*(max(0.,[0]+[1]*pow(x,[2]))-1)",10,ptmax);
  fhb->SetParameters(1.03091e+00, -5.11540e-02, -1.54227e-01); // SPRH
  fhb->SetLineColor(kMagenta+2);
  fhb->SetLineStyle(kDashed);
  fhb->SetLineWidth(3);
  fhb->Draw("SAME");

  tdrDraw(hhp3,"P",kOpenTriangleUp,kRed);  ht->SetMarkerSize(0.7);
  //tdrDraw(hhm3,"P",kOpenTriangleDown,kBlue);
  //tdrDraw(hep3,"P",kOpenTriangleDown,kRed);
  tdrDraw(hem3,"P",kOpenTriangleUp,kBlue);  ht->SetMarkerSize(0.7);

  tdrDraw(hcmm3,"P",kOpenCircle,kBlue);
  tdrDraw(hcp3,"P",kFullCircle,kRed);
  tdrDraw(hcm3,"P",kFullCircle,kBlue);

  //tdrDraw(ht0,"P",kOpenSquare,kGreen+2);
  //tdrDraw(hp0,"P",kOpenDiamond,kCyan+2);
  tdrDraw(ht,"P",kFullSquare,kGreen+2); ht->SetMarkerSize(0.7);
  //tdrDraw(hpe,"P",kOpenDiamond,kCyan+2);
  tdrDraw(hp,"P",kFullDiamond,kCyan+2);

  tdrDraw(hcx,"P",kFullCircle,kBlack);
  //tdrDraw(hcx2,"P",kFullCircle,kGray+2); hcx->SetMarkerSize(0.7);
  //tdrDraw(hcxp,"P",kOpenCircle,kBlack);
  //tdrDraw(hcp03,"P",kOpenCircle,kOrange+2);

  // Min corresponds to 10% (NHF) times 3% uncertainty = 0.3%
  // Max corresponds to 75% (hadrons) * 80% (HCAL fraction) times 3% = 1.8%,
  // where HCAL fraction at high pT is 35%*100% + 65%*68% = 79% ~ 80%
  // (pfhadrons_frac_AK.pdf and pfhadrons_ECAL_10_04.pdf)
  /*
  TF1 *fc3 = new TF1("fc3",Form("max(%1.1f,min(%1.1f,"
				"[0]+[1]*pow(x/208.,[2])+[3]*log(x)/x))",
				0.3,1.7),
		     15,ptmax);
  fc3->SetParameters(0.3,0.2,1,0);
  */
  /*
  TF1 *fc3 = new TF1("fc3","[0]+log(x/208.)*([1]+pow(log(x/208.),2)*("
		     "[2],ptmin,ptmax);
  fc3->SetParameters(0.1,0.01,0.003,-0.001,+0.0003);
  */
  TF1 *fc3 = new TF1("fc3","[0]+[1]*pow(x/[2],[3])/(1+pow(x/[2],[3]))*"
		     "(1-pow(x/[2],-[3]))",ptmin,ptmax);
  fc3->SetParameters(0.3,1.5,208.,1);
  hcp3->Fit(fc3,"RN");
  fc3->SetLineColor(kRed);
  fc3->SetLineWidth(3);
  fc3->SetLineStyle(kDotted);
  fc3->Draw("SAME");

  /*
  TF1 *fm3 = new TF1("fn3",Form("min(%1.1f,max(%1.1f,"
				"[0]+[1]*pow(x/208.,[2])+[3]*log(x)/x))",
				-0.3,-1.7),
		     ptmin,ptmax);
  fm3->SetParameters(-0.3,-0.2,1,0);
  */
  TF1 *fm3 = new TF1("fm3","[0]+[1]*pow(x/[2],[3])/(1+pow(x/[2],[3]))*"
		     "(1-pow(x/[2],-[3]))",ptmin,ptmax);
  fm3->SetParameters(-0.3,-1.5,208.,1);
  hcm3->Fit(fm3,"QRN");
  fm3->SetLineColor(kBlue);
  fm3->SetLineStyle(kDotted);
  fm3->SetLineWidth(3);
  fm3->Draw("SAME");

  TF1 *fh3 = new TF1("fh3","[0]+[1]*pow(x/[2],[3])/(1+pow(x/[2],[3]))*"
		     "(1-pow(x/[2],-[3]))",ptmin,ptmax);
  fh3->SetParameters(+0.3,+1.5,208.,1);
  hhp3->Fit(fh3,"QRN");
  fh3->SetLineColor(kRed);
  fh3->SetLineStyle(kDashed);
  fh3->SetLineWidth(2);
  fh3->Draw("SAME");

  TF1 *fe3 = new TF1("fe3","[0]+[1]*pow(x/[2],[3])/(1+pow(x/[2],[3]))*"
		     "(1-pow(x/[2],-[3]))",ptmin,ptmax);
  fe3->SetParameters(-0.1,-0.5,208.,1);
  hem3->Fit(fe3,"QRN");
  fe3->SetLineColor(kBlue);
  fe3->SetLineStyle(kDashed);
  fe3->SetLineWidth(2);
  fe3->Draw("SAME");

  /*
  TF1 *fx = new TF1("fx","max(-0.3,min(1.7,[0]+[1]*pow(x/208.,[2])))",ptmin,ptmax);
  fx->SetParameters(-0.2,0.2,1);
  */
  TF1 *fx = new TF1("fx","[0]+[1]*pow(x/[2],[3])/(1+pow(x/[2],[3]))*"
		    "(1-pow(x/[2],-[3]))",ptmin,ptmax);
  fx->SetParameters(-0.3,2.4,208.,1);
  hcx->Fit(fx,"QRN");
  fx->SetLineColor(kBlack);
  fx->Draw("SAME");

  /*
  TF1 *fc = new TF1("fc",Form("max(%1.1f,min(%1.1f,"
			      "[0]+[1]*pow(x/208.,[2])+[3]*log(x)/x))",
			      0.3-dc,1.7-dc),
		    ptmin,ptmax);
  fc->SetParameters(-0.2,0.2,1,0);
  */
  TF1 *fc = new TF1("fc","[0]+[1]*pow(x/[2],[3])/(1+pow(x/[2],[3]))*"
		    "(1-pow(x/[2],-[3]))",ptmin,ptmax);
  //fc->SetParameters(-0.7,1.5,208.,1);
  fc->SetParameters(fc3->GetParameter(0)-fc3->Eval(208.),
		    fc3->GetParameter(1), fc3->GetParameter(2),
		    fc3->GetParameter(3));
  hcp03->Fit(fc,"QRN");
  fc->SetLineColor(kOrange+2);
  //fc->Draw("SAME");

  TF1 *ft = new TF1("ft","[0]+[1]*pow(x/208.,[2])",ptmin,ptmax);
  ft->SetParameters(0,-0.15,-0.3);
  ht->Fit(ft,"QRN");
  ft->SetLineColor(kGreen+2);
  ft->Draw("SAME");

  TF1 *fp = new TF1("fp","[0]",ptmin,ptmax);
  fp->SetParameters(0);
  hp->Fit(fp,"QRN");
  fp->SetLineColor(kCyan+2);
  fp->Draw("SAME");

  TF1 *fpe = new TF1("fpe","[0]+[1]*pow(x/208.,[2])",ptmin,ptmax);
  fpe->SetParameters(0,-0.15,-0.3);
  hpe->Fit(fpe,"QRN");
  fpe->SetLineColor(kCyan+2);
  //fpe->Draw("SAME");

  //TLegend *leg1 = tdrLeg(0.40,0.65,0.70,0.90);
  //TLegend *leg1 = tdrLeg(0.38,0.65,0.68,0.90);
  TLegend *leg1 = tdrLeg(0.38,0.63,0.68,0.90);
  leg1->SetTextSize(0.040);
  //leg1->AddEntry(fhb,"Orig. SPRH +3%","L");
  //leg1->AddEntry(fhb,"Orig. HCAL +3%","L");
  //leg1->AddEntry(hcp3,"SPRH +3%","PL");
  leg1->AddEntry(hcp3,"SPR +3%","PL");
  leg1->AddEntry(hhp3,"HCAL Had. +3%","PL");
  //leg1->AddEntry(hcp3,"SPR HCAL +3%","PL");
  //leg1->AddEntry(hcx,"log-lin -3% to +3%","PL");
  //leg1->AddEntry(hcp03,Form("+3%% #Delta -%1.1f%%",dc),"PL");
  //leg1->AddEntry(ht,"Tracking -1%","PL");
  leg1->AddEntry(ht,"Tracking -3%","PL");
  //leg1->AddEntry(hp,"ECAL -1%","PL");
  //leg1->AddEntry(hp,"Photons -1%","PL");
  //leg1->AddEntry(hp,"ECAL -3%","PL");
  leg1->AddEntry(hp,"Photons -3%","PL");
  //leg1->AddEntry(hcm3,"SPRH -3%","PL");
  //leg1->AddEntry(hcm3,"SPR HCAL -3%","PL");
  leg1->AddEntry(hem3,"ECAL Had. -3%","PL");
  leg1->AddEntry(hcm3,"SPR -3%","PL");


  //TLegend *leg2 = tdrLeg(0.17,0.15,0.47,0.33);
  TLegend *leg2 = tdrLeg(0.17,0.150,0.47,0.305);
  leg2->SetTextSize(0.040);
  leg2->SetHeader("SPR modified");
  leg2->AddEntry(hcx,"log-lin -3% to +3%","PL");
  leg2->AddEntry(fhb,"Run I SPRH","L");
  //leg2->AddEntry(hcp03,Form("+3%% #Delta -%1.1f%%",dc),"PL");

  cout << endl;
  cout << "  // Fits from minitools/varPlots.C" << endl;

  cout << "  // SPR -3% to +3% cross variation" << endl;
  cout << Form("  if (!fxh) fxh = new TF1(\"fxh\",\"%s\",15,4500);\n",
	       fx->GetExpFormula().Data());
  cout << Form("  fxh->SetParameters(%1.4g, %1.4g, %1.4g, %1.4g); // toyPF\n\n",
	       fx->GetParameter(0), fx->GetParameter(1),
	       fx->GetParameter(2), fx->GetParameter(3));

  cout << "  // SPR +3% variation" << endl;
  cout << Form("  if (!fch) fch = new TF1(\"fch\",\"%s\",15,4500);\n",
	       fc->GetExpFormula().Data());
  cout << Form("  fch->SetParameters(%1.4g, %1.4g, %1.4g, %1.4g); // toyPF\n\n",
	       fc->GetParameter(0), fc->GetParameter(1),
	       fc->GetParameter(2), fc->GetParameter(3));

  cout << "  // SPRH +3% variation" << endl;
  cout << Form("  if (!fhh) fhh = new TF1(\"fhh\",\"%s\",15,4500);\n",
	       fh3->GetExpFormula().Data());
  cout << Form("  fhh->SetParameters(%1.4g, %1.4g, %1.4g, %1.4g); // toyPF\n\n",
	       fh3->GetParameter(0), fh3->GetParameter(1),
	       fh3->GetParameter(2), fh3->GetParameter(3));

  cout << "  // SPRE -3% variation" << endl;
  cout << Form("  if (!feh) feh = new TF1(\"feh\",\"%s\",15,4500);\n",
	       fe3->GetExpFormula().Data());
  cout << Form("  feh->SetParameters(%1.4g, %1.4g, %1.4g, %1.4g); // toyPF\n\n",
	       fe3->GetParameter(0), fe3->GetParameter(1),
	       fe3->GetParameter(2), fe3->GetParameter(3));

  cout << "  // Tracking -3% variation" << endl;
  cout << Form("  if (!ft) ft = new TF1(\"ft\",\"%s\",15,4500);\n",
	       ft->GetExpFormula().Data());
  cout << Form("  ft->SetParameters(%1.4g, %1.4g, %1.4g); // toyPF\n\n",
	       ft->GetParameter(0), ft->GetParameter(1), ft->GetParameter(2));

  cout << "  // Photon -3% variation" << endl;
  cout << Form("  if (!fp) fp = new TF1(\"fp\",\"%s\",15,4500);\n",
	       fp->GetExpFormula().Data());
  cout << Form("  fp->SetParameter(0,%1.4g); // toyPF\n\n",
	       fp->GetParameter(0));

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  tex->DrawLatex(0.77,0.87,"|#eta| < 1.3");

  c1->SaveAs("pdf/varPlotsResp.pdf");
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

  setTDRStyle();

  TDirectory *curdir = gDirectory;  
  //TFile *f = new TFile("rootfiles/varPlots_Mikael_5M_20200511.root","READ");
  TFile *f = new TFile("rootfiles/varPlots_Mikael_5M_20200514.root","READ");
  assert(f && !f->IsZombie());

  //TFile *f2 = new TFile("rootfiles/output-MCNU-2b-UL17V4_BCDEF.root","READ");
  //assert(f2 && !f2->IsZombie());
  //f2->cd("Standard/Eta_0.0-1.3");
  //TDirectory *d = gDirectory;

  curdir->cd();

  //string sf = "nhf";
  string sf2 = (sf=="gammaf" ? "nef" : sf);
  const char *cf = sf.c_str();
  const char *cf2 = sf2.c_str();
  TH1D *hcp = (TH1D*)f->Get(Form("h%s_Cp3",cf)); assert(hcp);
  TH1D *hcm = (TH1D*)f->Get(Form("h%s_Cm3",cf)); assert(hcm);
  TH1D *hh = (TH1D*)f->Get(Form("h%s_HadHCALp3",cf)); assert(hh);
  TH1D *he = (TH1D*)f->Get(Form("h%s_HadECALm3",cf)); assert(he);
  //TH1D *ht = (TH1D*)f->Get(Form("h%s_Trkm3",cf)); assert(ht);
  TH1D *ht1 = (TH1D*)f->Get(Form("h%s_Trkm1",cf)); assert(ht1);
  TH1D *ht3 = (TH1D*)f->Get(Form("h%s_Trkm3",cf)); assert(ht3);
  TH1D *hp = (TH1D*)f->Get(Form("h%s_ECALm3",cf)); assert(hp);

  //TProfile *p = (TProfile*)d->Get(Form("p%stp",cf2)); assert(p);
  //TGraphErrors *g = new TGraphErrors(p->ProjectionX(Form("p%s",cf2)));
  //_g = g;
  //TF1 *fc = new TF1("fc",_fc,ptmin,ptmax,0);

  /*
  // For 20200511 to turn ratio to diff
  TH1D *h1 = (TH1D*)hcp->Clone("h1Comp"); h1->Divide(h1);
  for (int i = 1; i != h1->GetNbinsX()+1; ++i) h1->SetBinError(i, 0);

  hcp->Add(hcp,h1,100,-100);
  hcm->Add(hcm,h1,100,-100);
  hh->Add(hh,h1,100,-100);
  he->Add(he,h1,100,-100);
  //ht->Add(ht,h1,100,-100);
  ht1->Add(ht1,h1,100,-100);
  ht3->Add(ht3,h1,100,-100);
  hp->Add(hp,h1,100,-100);
  hp->Add(he,-1);
  */

  /*
  // For 20200511 to turn ratio to diff
  hcp->Multiply(fc);
  hcm->Multiply(fc);
  hh->Multiply(fc);
  he->Multiply(fc);
  ht1->Multiply(fc);
  ht3->Multiply(fc);
  hp->Multiply(fc);
  */

  // For 20200514 to turn to percentage
  hcp->Scale(100);
  hcm->Scale(100);
  hh->Scale(100);
  he->Scale(100);
  ht1->Scale(100);
  ht3->Scale(100);
  hp->Scale(100);

  scaleStat(hcp,0.1);
  scaleStat(hcm,0.1);
  scaleStat(hh,0.1);
  scaleStat(he,0.1);
  //scaleStat(ht,0.1);
  scaleStat(ht1,0.1);
  scaleStat(ht3,0.1);
  scaleStat(hp,0.1);

  // Reproduce combination
  TH1D *hcx = (TH1D*)hcp->Clone(Form("hcx_%s",cf2)); hcx->Reset();
  for (int i = 1; i != hcx->GetNbinsX()+1; ++i) {
    double pt = hcp->GetBinCenter(i);
    double p = hcp->GetBinContent(i);
    double w = -1 + log(pt/15.)/log(208./15.);
    double x = log(pt/208.);
    hcx->SetBinContent(i, w*p);
    hcx->SetBinError(i, hcp->GetBinError(i));
  }

  TH1D *ha = (TH1D*)ht1->Clone(Form("ha_%s",cf2));
  ha->Add(ha,hp,1,1.);
  ha->Add(ha,hcx,1,1.5);
  //ha->Add(ha,hcm,1,1);
  ha->Add(ha,hh,1,-1.5);
  //ha->Add(ha,ht1,1,0.5);

  TH1D *h = tdrHist("hComp",
		    //Form("Composition change for %s (%%)",cf),
		    Form("Composition change for %s (10^{-2})",cf2),
		    -2.,3.,"p_{T} (GeV)",ptmin,ptmax);
  lumi_13TeV = "toyPF";
  TCanvas *c1 = tdrCanvas("c1Comp",h,4,11,kSquare);
  gPad->SetLogx();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(ptmin,0.,ptmax,0.);

  tdrDraw(hcp,"Pz",kFullCircle,kRed);
  tdrDraw(hcm,"Pz",kFullCircle,kBlue);
  tdrDraw(hcx,"Pz",kFullCircle,kBlack);
  tdrDraw(hh,"Pz",kOpenTriangleUp,kRed);   hh->SetMarkerSize(0.7);
  tdrDraw(he,"Pz",kOpenTriangleUp,kBlue);  he->SetMarkerSize(0.7);
  //tdrDraw(ht,"Pz",kFullSquare,kBlue);      ht->SetMarkerSize(0.7);
  tdrDraw(ht3,"Pz",kFullSquare,kBlue);      ht1->SetMarkerSize(0.7);
  tdrDraw(ht1,"Pz",kFullSquare,kBlue-9);    ht3->SetMarkerSize(0.7);
  tdrDraw(hp,"Pz",kFullDiamond,kCyan+2); 

  //tdrDraw(ha,"Pz",kFullDiamond,kOrange+2);

  TLegend *leg = tdrLeg(0.62,0.615,0.90,0.90);
  //TLegend *leg = tdrLeg(0.62,0.66,0.90,0.90);
  leg->AddEntry(ht3,"Tracker -3%","PLE");
  leg->AddEntry(ht1,"Tracker -1%","PLE");
  leg->AddEntry(hh,"HCAL-h +3%","PLE");
  leg->AddEntry(hcp,"SPR +3%","PLE");
  leg->AddEntry(hcm,"SPR -3%","PLE");
  leg->AddEntry(hp,"Photons -3%","PLE");
  leg->AddEntry(he,"ECAL-h -3%","PLE");

  //TLegend *leg2 = tdrLeg(0.20,0.57,0.48,0.66);
  TLegend *leg2 = tdrLeg(0.20,0.66,0.48,0.75);
  leg2->AddEntry(hcx,"SPR #pm3%","PLE");
  //leg2->AddEntry(ha,"T1+P3+H3+SPR#pm3%","PLE");

  c1->SaveAs(Form("pdf/varPlotsComp_%s.pdf",cf2));

  assert(_c1);
  assert(_leg);
  _c1->cd();

  map<string,int> color;
  color["chf"] = kRed;
  color["nhf"] = kGreen+2;
  color["nef"] = kBlue;

  map<string,int> marker;
  marker["chf"] = kFullCircle;
  marker["nhf"] = kFullDiamond;
  marker["nef"] = kFullSquare;

  TH1D *ha2 = (TH1D*)ht1->Clone(Form("ha2_%s",cf2));
  ha2->Add(ha2,hp,1,1.0);
  ha2->Add(ha2,hcx,1,1.5);
  //ha2->Add(ha2,hcm,1,1);
  ha2->Add(ha2,hh,1,-2.0);
  ha2->Add(ha2,he,1,1.0);
  //ha2->Add(ha2,ht1,1,0.5);

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.035);
  tex->DrawLatex(0.20,0.17,"P(-1.0%),CX(+4.5%),H(-4.5%),E(-3.0%)");

  ha2->SetMarkerSize(1.0);
  tdrDraw(ha2,"Pz",marker[sf2],color[sf2]);
  _leg->AddEntry(ha2,cf2,"PLE");
} // varPlotsComp
