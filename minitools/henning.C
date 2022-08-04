// Purpose: analyse Henning's CHS/PUPPI JMENANO files
// NB: currently consuming lots of memory and very slow to finish up
#include "TFile.h"
#include "TString.h"
#include "TF1.h"
#include "TMath.h"
#include "TH2D.h"

#include "TGraphErrors.h"
#include "TLine.h"
#include "TLatex.h"

#include "../tdrstyle_mod15.C"
#include "../tools.C"

#include <vector>

bool debug = false;
bool interpolateMedian = true;
double Median(const TH1D * h1) { 

   int n = h1->GetXaxis()->GetNbins();  
   std::vector<double>  x(n);
   h1->GetXaxis()->GetCenter( &x[0] );
   const double * y = h1->GetArray(); 
   if (interpolateMedian) {
     // exclude underflow/overflows from bin content array y
     double xmid = TMath::Median(n, &x[0], &y[1]);
     int ix = h1->GetXaxis()->FindBin(xmid);
     double nlow = h1->Integral(1,ix-1);
     double nhigh = h1->Integral(1,ix);
     //double nhigh = nlow+h1->Integral(ix,ix);
     double ntot = h1->Integral();
     //double ntot = h1->Integral(0,h1->GetNbinsX()+1);
     //double ntot = nhigh+h1->Integral(ix+1,h1->GetNbins()+1);
     if (!(nlow<=0.5*ntot) || !(nhigh>=0.5*ntot)) {
       cout << endl << "nlow="<<nlow<<", nhigh="<<nhigh<<", ntot="<<ntot
	    << ", 0.5*ntot="<<0.5*ntot << endl << flush;
       cout << "0.5*ntot-nlow="<<(0.5*ntot-nlow)
	    << ", nhigh-0.5*ntot="<<(nhigh-0.5*ntot) << endl << flush;
     }
     //assert(nlow<=0.5*ntot);
     assert(0.5*ntot-nlow>=-2e-15);
     //assert(nhigh>=0.5*ntot);
     assert(nhigh-0.5*ntot>=-1e-15);
     double w = (nhigh>nlow ? (0.5*ntot-nlow) / (nhigh-nlow) : 0);
     return (h1->GetBinLowEdge(ix) + w * h1->GetBinWidth(ix));
   }
   else { // bin edges only
     // exclude underflow/overflows from bin content array y
     return TMath::Median(n, &x[0], &y[1]); 
   }
} // Median

// Calculate (95%) confidence interval around median
// => could replace with TH1D::GetQuantiles or TMath::Quantiles?
//    Do they interpolate within bins?
double Confidence(const TH1D *h1, double median, double confLevel = 0.95) {

  int ix = h1->GetXaxis()->FindBin(median);
  int ixlow = ix;
  int ixhigh = ix;
  int nb = h1->GetNbinsX();
  double ntot = h1->Integral();
  double nsum = h1->GetBinContent(ix);
  double width = h1->GetBinWidth(ix);
  if (ntot==0) return 0;
  while (nsum < confLevel * ntot) {
    double nlow = (ixlow>0 ? h1->GetBinContent(ixlow-1) : 0);
    double nhigh = (ixhigh<nb ? h1->GetBinContent(ixhigh+1) : 0);
    if (nsum+max(nlow,nhigh) < confLevel * ntot) {
      if (nlow>=nhigh && ixlow>0) {
	nsum += nlow; --ixlow; width += h1->GetBinWidth(ixlow);
      }
      else if (ixhigh<nb) {
	nsum += nhigh; ++ixhigh; width += h1->GetBinWidth(ixhigh); 
      }
      else
	assert(false);
    }
    else {
      if (nlow>nhigh) {
	width += h1->GetBinWidth(ixlow-1) * (confLevel * ntot - nsum) / nlow;
      }
      else {
	width += h1->GetBinWidth(ixhigh+1) * (confLevel * ntot - nsum) / nhigh;
      }
      nsum = ntot; // finish loop
    }
  } // while
  return width;
} // Confidence

bool isCorr = false;
double ptmin = 10.;
void henning() {

  double ptrec = (isCorr ? ptmin : ptmin * 0.8);

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  //TFile *fp = new TFile("rootfiles/Henning_CHSPUPPIcorrected_ptresponse_export.root","READ");
  TFile *fp(0);
  if (isCorr)
    fp =new TFile("rootfiles/Henning_corrected_ptresponse_export.root","READ");
  else
    fp = new TFile("rootfiles/Henning_ptresponse_export.root","READ");
  assert(fp && !fp->IsZombie());
  TFile *fg = new TFile("rootfiles/Henning_GenJetConts_export.root","READ");
  assert(fg && !fg->IsZombie());

  double ptbins_corr[] = 
    //{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 20,
    {11, 12, 13, 14, 15, 17, 20,
     23, 27, 30, 35, 40, 45, 57, 72, 90, 120, 
     150, 200, 300, 400, 550, 750, 1000, 1500, 2000, 2500, 3000,
     3500, 4000, 4500, 5000, 10000};
  const int npt_corr = sizeof(ptbins_corr)/sizeof(ptbins_corr[0])-1;
  /*
  double ptbins_uncorr[] = 
    {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 20,
     23, 27, 30, 35, 40, 45, 57, 72, 90, 120, 
     150, 200, 300, 400, 550, 750, 1000, 1500, 2000, 2500, 3000,
     3500, 4000, 4500, 5000, 10000};
  */
  double ptbins_uncorr[] = 
    {8, 9, 10, 11, 12, 13, 14, 15, 17, 20,
     23, 27, 30, 35, 40, 45, 57, 72, 90, 120, 
     150, 200, 300, 400, 550, 750, 1000, 1500, 2000, 2500, 3000,
     3500, 4000, 4500, 5000, 10000};
  const int npt_uncorr = sizeof(ptbins_uncorr)/sizeof(ptbins_uncorr[0])-1;
  double *ptbins = (isCorr ? ptbins_corr : ptbins_uncorr);
  const int npt = (isCorr ? npt_corr : npt_uncorr);
  double etabins[] =
    {-4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, -3.489,
     -3.314, -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172, -2.043,
     -1.93, -1.83, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218,
     -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, -0.435, 
     -0.348, -0.261, -0.174, -0.087, 0, 0.087, 0.174, 0.261, 0.348, 0.435,
     0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305,
     1.392, 1.479, 1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5,
     2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191,
     4.363, 4.538, 4.716, 4.889}; //]#5.191 ]
    //{0, 0.087, 0.174, 0.261, 0.348, 0.435};
  //{0.261, 0.348};
    //{0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305};
    //{0.609, 0.696};
    //{1.392, 1.479, 1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5};
    //{2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191};
  //{3.664, 3.839};
     //4.363, 4.538, 4.716, 4.889}; //]#5.191 ]
  //double etabins[] = {0, 0.087};
  const int neta = sizeof(etabins)/sizeof(etabins[0])-1;

  TLine *l = new TLine();
  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);

  // Record results over all eta and pT in TH2Ds
  TH2D *h2eff_any = new TH2D("h2eff_any",";p_{T} (GeV);#eta",
			     npt,ptbins,neta,etabins);
  TH2D *h2eff_pt8 = new TH2D("h2eff_pt8",";p_{T} (GeV);#eta",
			     npt,ptbins,neta,etabins);
  TH2D *h2eff_core = new TH2D("h2eff_core",";p_{T} (GeV);#eta",
			     npt,ptbins,neta,etabins);

  TH2D *h2jes_med = new TH2D("h2jes_med",";p_{T} (GeV);#eta",
			     npt,ptbins,neta,etabins);
  TH2D *h2jes_mean = new TH2D("h2jes_mean",";p_{T} (GeV);#eta",
			      npt,ptbins,neta,etabins);
  TH2D *h2jes_gaus = new TH2D("h2jes_gaus",";p_{T} (GeV);#eta",
			      npt,ptbins,neta,etabins);

  TH2D *h2jer_qrt = new TH2D("h2jer_qrt",";p_{T} (GeV);#eta",
			     npt,ptbins,neta,etabins);
  TH2D *h2jer_rms = new TH2D("h2jer_rms",";p_{T} (GeV);#eta",
			     npt,ptbins,neta,etabins);
  TH2D *h2jer_gaus = new TH2D("h2jer_gaus",";p_{T} (GeV);#eta",
			      npt,ptbins,neta,etabins);

  for (int ieta = 0; ieta != neta; ++ieta) {

    cout << endl << ":" << flush;

    TCanvas *c0lin = new TCanvas(Form("c0lin_%d",ieta),"c0lin",1600,1200);
    int ny = int(sqrt(npt/1.5)+0.5);
    int nx = int(ceil(double(npt)/ny));
    c0lin->Divide(nx,ny,0,0);
    TH1D *h0 = tdrHist(Form("h0_%d",ieta),"dN/dR/N",0.,5.,
		       Form("RECO / GEN (%s)", isCorr ? "corr." : "raw"),0,2);

    TCanvas *c0log = new TCanvas(Form("c0log_%d",ieta),"c0log",1600,1200);
    c0log->Divide(nx,ny,0,0);
    TH1D *h0log = tdrHist(Form("h0log_%d",ieta),"dN/dR/N",1e-4,2e1,
			  Form("RECO / GEN (%s)",isCorr ? "corr." : "raw"),0,2);
    
    lumi_13TeV = "Run 2 Legacy";
    TH1D *h = tdrHist(Form("h_%d",ieta),"dN/dR/N",0.,2.,
		      Form("RECO / GEN (%s)", isCorr ? "corr." : "raw"),0,3);
    TCanvas *c1 = tdrCanvas(Form("c1_%d",ieta),h,4,11,kSquare);

    TH1D *h2up = tdrHist(Form("h2up_%d",ieta),
			 Form("RECO / GEN (%s)", isCorr ? "corr." : "raw"),
			 0.65,1.15,"p_{T} (GeV)",10,4500);
    TH1D *h2dw = tdrHist(Form("h2dw_%d",ieta),"Ratio",0.98,1.02,
			 "p_{T} (GeV)",10,4500);
    TCanvas *c2 = tdrDiCanvas(Form("c2_%d",ieta),h2up,h2dw,4,11);
    
    TH1D *h3up = tdrHist(Form("h3up_%d",ieta),
			 Form("RECO / GEN (%s)", isCorr ? "corr." : "raw"),
			 0.,0.5,"p_{T} (GeV)",10,4500);
    TH1D *h3dw = tdrHist(Form("h3dw_%d",ieta),"Ratio",0.8,1.4,
			 "p_{T} (GeV)",10,4500);
    TCanvas *c3 = tdrDiCanvas(Form("c3_%d",ieta),h3up,h3dw,4,11);
    
    TH1D *h4up = tdrHist(Form("h4up_%d",ieta),
			 Form("RECO / GEN (%s)", isCorr ? "corr." : "raw"),
			 0.65,1.15,"p_{T} (GeV)",10,4500);
    TH1D *h4dw = tdrHist(Form("h4dw_%d",ieta),"Ratio",0.98,1.02,
			 "p_{T} (GeV)",10,4500);
    TCanvas *c4 = tdrDiCanvas(Form("c4_%d",ieta),h4up,h4dw,4,11);

    if (debug) cout << "c1s-" << flush;
    
    TGraphErrors *gmean = new TGraphErrors(0);
    TGraphErrors *gmeanp = new TGraphErrors(0);
    TGraphErrors *gmed = new TGraphErrors(0);
    TGraphErrors *gmedp = new TGraphErrors(0);
    TGraphErrors *gmu = new TGraphErrors(0);
    TGraphErrors *gmup = new TGraphErrors(0);
    TGraphErrors *gmun = new TGraphErrors(0);
    TGraphErrors *gmupn = new TGraphErrors(0);
    TGraphErrors *gmun2 = new TGraphErrors(0);
    TGraphErrors *gmupn2 = new TGraphErrors(0);
    TGraphErrors *gnorm0 = new TGraphErrors(0);
    TGraphErrors *gnorm1 = new TGraphErrors(0);
    TGraphErrors *gnorm = new TGraphErrors(0);
    TGraphErrors *gnormp = new TGraphErrors(0);
    
    TGraphErrors *grms = new TGraphErrors(0);
    TGraphErrors *gconf68 = new TGraphErrors(0);
    TGraphErrors *gconf90 = new TGraphErrors(0);
    TGraphErrors *gconf95 = new TGraphErrors(0);
    TGraphErrors *gconf98 = new TGraphErrors(0);
    TGraphErrors *gconf99 = new TGraphErrors(0);
    TGraphErrors *gsigma = new TGraphErrors(0);
    TGraphErrors *gsigmap = new TGraphErrors(0);
    TGraphErrors *gsigman = new TGraphErrors(0);
    
    grms->SetName("grms");
    gconf68->SetName("gconf68");
    gconf90->SetName("gconf90");
    gconf95->SetName("gconf95");
    gconf98->SetName("gconf98");
    gconf99->SetName("gconf99");
    gsigma->SetName("gsigma");
    gsigmap->SetName("gsigmap");
    gsigman->SetName("gsigman");
    
    if (debug) cout << "grs-" << flush;

    double eta1 = etabins[ieta];
    double eta2 = etabins[ieta+1];
    TString teta = Form("GenJetCounts___eta_%06.3f_to_%06.3f",eta1,eta2);
    teta.ReplaceAll(".","_");
    TH1D *hg = (TH1D*)fg->Get(teta.Data());
    if (!hg) {
      cout << teta << endl << flush;
      assert(hg);
      continue;
    }
    
    for (int ipt = 0; ipt != npt; ++ipt) {
      
      cout << "." << flush;
      double pt1 = ptbins[ipt];
      double pt2 = ptbins[ipt+1];
      double pt = 0.5*(pt1+pt2);

      // Skip bins with part in unphysical regions
      //if (pt2*cosh(eta2)>6500.) continue;
      // Skip bins with all in unphysical regions
      //if (pt1*cosh(eta1)>6500.) continue;

      //TString tpt = Form("CHSPUPPIcorrected_ptresponse__"
      //TString tpt = Form("corrected_ptresponse__"
      TString tpt = Form("%sptresponse__"
			 "pT_%06.0f_to_%06.0f__"
			 "eta_%06.3f_to_%06.3f",
			 isCorr ? "corrected_" : "",
			 pt1,pt2,eta1,eta2);
      tpt.ReplaceAll(".","_");
      TH1D *hp = (TH1D*)fp->Get(tpt.Data());
      if (!hp) {
	cout << tpt << endl << flush;
	assert(hp);
	continue;
      }

      double ng = hg->GetBinContent(hg->FindBin(pt));
      double nw = hp->GetBinWidth(hp->FindBin(1.0));
      hp->Scale(ng!=0 ? 1./ng/nw : 1.);
      double err, err2;
      double np = hp->Integral() * nw;
      hp->IntegralAndError(-1,-1,err); err *= nw;
      if (np==0) continue;

      // Zero out pTreco<10 GeV
      double rms = hp->GetRMS();
      double mean = hp->GetMean();
      double xmin0 = ptrec/pt1;
      double xmin = xmin0;
      TH1D *hp2 = (TH1D*)hp->Clone(Form("hp2_%d_%d",ieta,ipt));
      for (int ix = 0; ix != hp->GetNbinsX()+1; ++ix) {
	if (hp2->GetBinLowEdge(ix)<xmin0) {
	  hp2->SetBinContent(ix,0);
	  hp2->SetBinError(ix,0);
	  if (hp2->GetBinLowEdge(ix+1)>xmin0)
	    xmin = hp2->GetBinLowEdge(ix+1);
	}
      }

      double np2 = hp2->Integral() * nw;
      hp2->IntegralAndError(-1,-1,err2); err2 *= nw;
      if (np2==0) continue;

      double krms = 2.0;
      double ksigma = 2.0;
      double mus = mean;
      double sig = rms;
      double norm = np;//norm

      if (debug) cout << "hp-" << flush;

      // Normalized Gaussian, full range
      TF1 *f1 = new TF1(Form("f1_%d_%d",ieta,ipt),
			"TMath::Gaus(x,[0],[1],1)*[2]",0.,5);
      f1->SetLineColor(kMagenta+1);
      f1->SetRange(max(0.,mean-krms*rms),mean+krms*rms);
      f1->SetParameters(mean,rms,1);
      f1->FixParameter(2,norm);
      f1->SetParLimits(1,0.01,max(0.1,1.0*mean));
      hp->Fit(f1,"QRN");

      // Normalized Gaussian, partial range
      TF1 *f1p = new TF1(Form("f1p_%d_%d",ieta,ipt),
			 "TMath::Gaus(x,[0],[1],1)*[2]",0.,5);
      f1p->SetLineColor(kRed);
      f1p->SetRange(max(xmin,mean-krms*rms),mean+krms*rms);
      f1p->SetParameters(mean,rms,1);
      f1p->FixParameter(2,norm);
      hp2->Fit(f1p,"QRN");

      // Re-normalized Gaussian, full range
      TF1 *f1n = new TF1(Form("f1n_%d_%d",ieta,ipt),
			 "TMath::Gaus(x,[0],[1],1)*[2]",0.,5);
      f1n->SetLineColor(kGreen+2);
      f1n->FixParameter(0,f1->GetParameter(0));
      f1n->FixParameter(1,f1->GetParameter(1));
      f1n->SetParLimits(2,0.5,1.5);
      f1n->SetRange(f1->GetXmin(), f1->GetXmax());
      f1n->SetParameter(2,1);
      hp->Fit(f1n,"QRN");

      // Re-normalized Gaussian, partial range
      TF1 *f1pn = new TF1(Form("f1pn_%d_%d",ieta,ipt),
			  "TMath::Gaus(x,[0],[1],1)*[2]",0.,5);
      f1pn->SetLineColor(kOrange+1);
      f1pn->FixParameter(0,f1->GetParameter(0));
      f1pn->FixParameter(1,f1->GetParameter(1));
      f1pn->SetRange(f1p->GetXmin(), f1p->GetXmax());
      f1pn->SetParameter(2,1);
      hp2->Fit(f1pn,"QRN");

      if (debug) cout << "f1s-" << flush;

      // Un-normalized Gaussian fit to full range
      TF1 *f0 = new TF1(Form("f0_%d_%d",ieta,ipt),
			//"gaus",max(0.,mu-ksigma*sigma),mu+ksigma*sigma);
			"gaus",max(0.,mus-ksigma*sig),mus+ksigma*sig);
      f0->SetLineColor(kCyan+1);
      f0->SetParameters(1,f1->GetParameter(0),f1->GetParameter(1));
      hp->Fit(f0,"QRN");

      // Un-normalized Gaussian fit to partial range
      TF1 *f0p = new TF1(Form("f0p_%d_%d",ieta,ipt),
			 //"gaus",xmin,mu+ksigma*sigma);
			 "gaus",xmin,mus+ksigma*sig);
      f0p->SetLineColor(kBlue);
      f0p->SetParameters(1,f1p->GetParameter(0),f1p->GetParameter(1));
      hp2->Fit(f0p,"QRN");

      if (debug) cout << "f0s-" << flush;

      double eff(0.), eff8(0.), effc(0.);
      double med(0.), /*mean(0.),*/ mu(0.), qrt(0.), /*rms(0.),*/ sigma(0.);
      double efferr(0.), eff8err(0.), effcerr(0.);
      double meanerr(0.), muerr(0.), rmserr(0.), sigmaerr(0.);
      if (hp->GetMean()!=0) {
	int jpt = gmean->GetN();
	//mean = hp->GetMean();
	//meanerr = hp->GetMeanError();
	gmean->SetPoint(jpt, pt, mean);
	gmean->SetPointError(jpt, 0, meanerr);
	gmeanp->SetPoint(jpt, pt, hp2->GetMean());
	gmeanp->SetPointError(jpt, 0, hp2->GetMeanError());
	//
	med = Median(hp);
	gmed->SetPoint(jpt, pt, med);
	gmedp->SetPoint(jpt, pt, Median(hp2));
	//
	//mu = f0->GetParameter(1);
	//muerr = f0->GetParError(1);
	//gmu->SetPoint(jpt, pt, mu);
	//gmu->SetPointError(jpt, 0, muerr);
	gmu->SetPoint(jpt, pt, f0->GetParameter(1));
	gmu->SetPointError(jpt, 0, f0->GetParError(1));
	gmup->SetPoint(jpt, pt, f0p->GetParameter(1));
	gmup->SetPointError(jpt, 0, f0p->GetParError(1));
	
	mu = f1->GetParameter(0);
	muerr = f1->GetParError(0);
	sigma = f1->GetParameter(1)/f1->GetParameter(0);
	sigmaerr = f1->GetParError(1)/f1->GetParameter(0);

	gmun->SetPoint(jpt, pt, f1->GetParameter(0));
	gmun->SetPointError(jpt, 0, f1->GetParError(0));
	gmupn->SetPoint(jpt, pt, f1p->GetParameter(0));
	gmupn->SetPointError(jpt, 0, f1p->GetParError(0));
	gmun2->SetPoint(jpt, pt, f1n->GetParameter(0));
	gmun2->SetPointError(jpt, 0, f1n->GetParError(0));
	gmupn2->SetPoint(jpt, pt, f1pn->GetParameter(0));
	gmupn2->SetPointError(jpt, 0, f1pn->GetParError(0));
	
	if (debug) cout << "med-" << flush;

	eff = np;
	efferr = err;
	eff8 = np2;
	eff8err = err2;
	effc = f1n->GetParameter(2);
	effcerr = f1n->GetParError(2);
	gnorm0->SetPoint(jpt, pt, eff);
	gnorm1->SetPoint(jpt, pt, eff8);
	gnorm->SetPoint(jpt, pt, effc);
	gnormp->SetPoint(jpt, pt, f1pn->GetParameter(2));
	
	rms = hp->GetRMS() / hp->GetMean();
	rmserr = hp->GetRMSError() / hp->GetMean();
	grms->SetPoint(jpt, pt, rms);
	grms->SetPointError(jpt, 0., rmserr);
	
	// 99% confidence is [-2.5sigma,+2.5sigma] so divide by 2*2.576
	gconf99->SetPoint(jpt, pt, Confidence(hp, med, 0.99)/(2*2.576*med));
	if (debug) cout << "qrt99-" << flush;
	// 98% confidence is [-2.3sigma,+2.3sigma] so divide by 2*2.326
	gconf98->SetPoint(jpt, pt, Confidence(hp, med, 0.98)/(2*2.326*med));
	// 95% confidence is [-2sigma,+2sigma] so divide by 2*1.960
	gconf95->SetPoint(jpt, pt, Confidence(hp, med, 0.95)/(2*1.960*med));
	// 90% confidence is [-1.6*sigma,+1.6sigma] so divide by 2*1.645
	gconf90->SetPoint(jpt, pt, Confidence(hp, med, 0.90)/(2*1.645*med));
	// 68% confidence is [-sigma,+sigma] so divide by 2*0.9945
	qrt = Confidence(hp, med, 0.68)/(2*0.9945*med);
	gconf68->SetPoint(jpt, pt, qrt);
	
	if (debug) cout << "qrt-" << flush;

	//sigma = f0->GetParameter(2)/f0->GetParameter(1);
	//sigmaerr = f0->GetParError(2)/f0->GetParameter(1);
	//gsigma->SetPoint(jpt, pt, sigma);
	//gsigma->SetPointError(jpt, 0, sigmaerr);
	gsigma->SetPoint(jpt, pt, f0->GetParameter(2)/f0->GetParameter(1));
	gsigma->SetPointError(jpt,0,f0->GetParError(2)/f0->GetParameter(1));
	gsigmap->SetPoint(jpt, pt, f0p->GetParameter(2)/f0p->GetParameter(1));
	gsigmap->SetPointError(jpt,0,f0p->GetParError(2)/f0p->GetParameter(1));
	gsigman->SetPoint(jpt, pt, f1p->GetParameter(1)/f1p->GetParameter(0));
	gsigman->SetPointError(jpt,0,f1p->GetParError(1)/f1p->GetParameter(2));
      }

      // Draw plots
      if (fabs(pt1-11)<0.5) {
	c1->cd();
	TF1 *f1p_c = (TF1*)f1p->Clone("f1p_c");
	//f1p_c->SetRange(mu-ksigma*sigma,mu+ksigma*sigma);
	f1p_c->SetRange(max(0.,mus-ksigma*sig),mus+ksigma*sig);
	f1p_c->SetLineStyle(kDotted);
	TF1 *f0p_c = (TF1*)f0p->Clone("f0p_c");
	//f0p_c->SetRange(mu-ksigma*sigma,mu+ksigma*sigma);
	f0p_c->SetRange(max(0.,mus-ksigma*sig),mus+ksigma*sig);
	f0p_c->SetLineStyle(kDotted);
	
	tdrDraw(hp,"H",kNone,kGray+1,kSolid,-1,1001,kGray);
	tdrDraw(hp2,"H",kNone,kYellow+2,kSolid,-1,1001,kYellow+1);
	f1->Draw("SAME");
	f1p_c->Draw("SAME");
	f1p->Draw("SAME");
	f1n->Draw("SAME");
	f1pn->Draw("SAME");
	f0->Draw("SAME");
	f0p_c->Draw("SAME");
	f0p->Draw("SAME");

	tex->SetTextAlign(31); // Left
	tex->DrawLatex(0.90,0.85,Form("%1.0f#leqp_{T}<%1.0f GeV",pt1,pt2));
	tex->DrawLatex(0.90,0.79,Form("%1.3f#leq#eta<%1.3f",eta1,eta2));
	
	gPad->RedrawAxis();

	if (debug) cout << "c1-" << flush;
      }

      c0lin->cd(ipt+1);
      TH1D *h0i = (TH1D*)h0->DrawClone(Form("h0_%d_%d",ieta+1,ipt+1));
      tex->SetTextSize(0.045*2); tex->SetTextAlign(31); // Left
      tex->DrawLatex(0.95,0.90,Form("%1.0f#leqp_{T}<%1.0f GeV",pt1,pt2));
      if (ipt==0)
	tex->DrawLatex(0.95,0.80,Form("%1.3f#leq#eta<%1.3f",eta1,eta2));
      double maxy[] = {2.0,3.0,6.0,10.0,10.0};
      assert(sizeof(maxy)/sizeof(maxy[0])>=ny);
      h0i->SetMaximum(maxy[ipt/nx]);
      tdrDraw(hp,"H",kNone,kGray+1,kSolid,-1,1001,kGray);
      tdrDraw(hp2,"H",kNone,kYellow+2,kSolid,-1,1001,kYellow+1);
      f1->Draw("SAME");
      gPad->RedrawAxis();

      c0log->cd(ipt+1);
      gPad->SetLogy();
      TH1D *h0logi = (TH1D*)h0log->DrawClone(Form("h0log_%d_%d",ieta+1,ipt+1));
      tex->DrawLatex(0.95,0.90,Form("%1.0f#leqp_{T}<%1.0f GeV",pt1,pt2));
      if (ipt==0)
	tex->DrawLatex(0.95,0.80,Form("%1.3f#leq#eta<%1.3f",eta1,eta2));
      tex->SetTextSize(0.045); tex->SetTextAlign(11); // Right
      tdrDraw(hp,"H",kNone,kGray+1,kSolid,-1,1001,kGray);
      tdrDraw(hp2,"H",kNone,kYellow+2,kSolid,-1,1001,kYellow+1);
      f1->Draw("SAME");
      TF1 *f1f = (TF1*)f1->Clone(Form("f1f_%d_%d",ieta,ipt));
      f1f->SetLineStyle(kDotted);
      f1f->SetRange(0,3);
      f1f->Draw("SAME");
      gPad->RedrawAxis();

      if (debug) cout << "c0-" << flush;

      h2eff_any->SetBinContent(ipt+1, ieta+1, eff);
      h2eff_any->SetBinError(ipt+1, ieta+1, efferr);
      h2eff_pt8->SetBinContent(ipt+1, ieta+1, eff8);
      h2eff_pt8->SetBinError(ipt+1, ieta+1, eff8err);
      h2eff_core->SetBinContent(ipt+1, ieta+1, eff>0 ? effc / eff : 0);
      h2eff_core->SetBinError(ipt+1, ieta+1, eff>0 ? effcerr / eff : 0);

      h2jes_gaus->SetBinContent(ipt+1, ieta+1, mu);
      h2jes_gaus->SetBinError(ipt+1, ieta+1, muerr);
      h2jer_gaus->SetBinContent(ipt+1, ieta+1, sigma);
      h2jer_gaus->SetBinError(ipt+1, ieta+1, sigmaerr);

      h2jes_med->SetBinContent(ipt+1, ieta+1, med);
      h2jes_med->SetBinError(ipt+1, ieta+1, meanerr);
      h2jer_qrt->SetBinContent(ipt+1, ieta+1, qrt);
      h2jer_qrt->SetBinError(ipt+1, ieta+1, rmserr);

      h2jes_mean->SetBinContent(ipt+1, ieta+1, mean);
      h2jes_mean->SetBinError(ipt+1, ieta+1, meanerr);
      h2jer_rms->SetBinContent(ipt+1, ieta+1, rms);
      h2jer_rms->SetBinError(ipt+1, ieta+1, rmserr);
      
      if (debug) cout << "h2-" << endl << flush;
    } // for ipt

    //if (ieta!=neta/2) continue;
    
    // Efficiencies
    c4->cd(1);
    gPad->SetLogx();
    tdrDraw(gnorm0,"Lz",kNone,kBlack); // full histogram
    tdrDraw(gnorm1,"Lz",kNone,kGray+1); // partial histogram
    tdrDraw(gnorm,"Lz",kNone,kGreen+2); // norm of renormalized Gaus, fill range
    //tdrDraw(gnormp,"Lz",kNone,kCyan+2); // norm of ren. Gaus, partial range

    TLegend *leg1b = tdrLeg(0.40,0.90-0.05*3,0.70,0.90);
    leg1b->AddEntry(gnorm0,"Any pT","PL");
    leg1b->AddEntry(gnorm1,Form("pT>%1.1f",ptrec),"PL");
    leg1b->AddEntry(gnorm,"Core (any pT)","PL");
    //leg1b->AddEntry(gnormp,Form("Core (pT>%1.1f)",ptrec),"PL");

    c4->cd(2);
    gPad->SetLogx();

    // JES
    c2->cd(1);
    gPad->SetLogx();
    l->SetLineStyle(kDashed);
    l->DrawLine(10,1,4500,1);
    tex->DrawLatex(0.66,0.86,Form("%1.3f<|#eta|<%1.3f",eta1,eta2));

    tdrDraw(gmed,"Lz",kNone,kGray+3); // mean down to zero
    tdrDraw(gmedp,"Lz",kNone,kYellow+3); // mean at pTreco>10
    tdrDraw(gmean,"Lz",kNone,kGray+2); // mean down to zero
    tdrDraw(gmeanp,"Lz",kNone,kYellow+2); // mean at pTreco>10
    gmed->SetLineWidth(2);
    gmedp->SetLineWidth(2);
    tdrDraw(gmu,"Lz",kNone,kCyan+1); // free Gaus in full range
    tdrDraw(gmup,"Lz",kNone,kBlue); // free Gaus in partial range
    //tdrDraw(gmun,"Lz",kNone,kMagenta+1); // normalized Gaus in full range
    tdrDraw(gmupn,"Lz",kNone,kRed); // normalized Gaus in partial range
    gPad->RedrawAxis();

    TLegend *leg1 = tdrLeg(0.40,0.05,0.70,0.05+0.05*8);
    leg1->AddEntry(gmed,"Median (any pT)","PL");
    leg1->AddEntry(gmedp,Form("Median (pT>%1.1f)",ptrec),"PL");
    leg1->AddEntry(gmean,"Mean (any pT)","PLE");
    leg1->AddEntry(gmeanp,Form("Mean (pT>%1.1f)",ptrec),"PLE");
    leg1->AddEntry(gmu,"Gaus (any pT)","PLE");
    leg1->AddEntry(gmup,Form("Gaus (pT>%1.1f)",ptrec),"PLE");
    //leg1->AddEntry(gmun,"NGaus (any pT)","PLE");
    leg1->AddEntry(gmupn,Form("NGaus (pT>%1.1f)",ptrec),"PLE");

    c2->cd(2);
    gPad->SetLogx();
    //l->DrawLine(10,1,4500,1);
    TGraphErrors *grmean = tools::ratioGraphs(gmean,gmed);
    TGraphErrors *grmeanp = tools::ratioGraphs(gmeanp,gmed);
    TGraphErrors *grmed = tools::ratioGraphs(gmed,gmed);
    TGraphErrors *grmedp = tools::ratioGraphs(gmedp,gmed);
    TGraphErrors *grmu = tools::ratioGraphs(gmu,gmed);
    TGraphErrors *grmup = tools::ratioGraphs(gmup,gmed);
    TGraphErrors *grmun = tools::ratioGraphs(gmun,gmed);
    TGraphErrors *grmupn = tools::ratioGraphs(gmupn,gmed);
    tdrDraw(grmed,"L0",kNone,kGray+3); // median at pTreco>0
    tdrDraw(grmedp,"L0",kNone,kYellow+3); // median at pTreco>10
    tdrDraw(grmean,"L",kNone,kGray+2); // mean at pTreco>0
    tdrDraw(grmeanp,"L",kNone,kYellow+2); // mean at pTreco>10
    grmed->SetLineWidth(2);
    grmedp->SetLineWidth(2);
    tdrDraw(grmu,"L",kNone,kCyan+1); // free Gaus, full range
    tdrDraw(grmup,"L",kNone,kBlue); // free Gaus, partial range
    //tdrDraw(grmun,"L",kNone,kMagenta+1); // normalized Gaus, full range
    tdrDraw(grmupn,"L",kNone,kRed); // normalized Gaus, partial range
    gPad->RedrawAxis();


    c3->cd(1);
    gPad->SetLogx();
    tex->DrawLatex(0.66,0.86,Form("%1.3f<|#eta|<%1.3f",eta1,eta2));

    tdrDraw(gconf68,"L0",kNone,kGray+1);
    tdrDraw(gconf90,"L0",kNone,kGray+2);
    //tdrDraw(gconf95,"L0",kNone,kGray+3);
    tdrDraw(gconf98,"L0",kNone,kGray+3);
    //tdrDraw(gconf99,"L0",kNone,kGray+3);
    tdrDraw(grms,"L",kNone,kGray+2);
    tdrDraw(gsigma,"L",kNone,kCyan+1);
    tdrDraw(gsigmap,"L",kNone,kBlue);
    tdrDraw(gsigman,"L",kNone,kRed);

    c3->cd(2);
    gPad->SetLogx();

    TGraphErrors *grconf68 = tools::ratioGraphs(gconf68,gconf68);
    TGraphErrors *grconf90 = tools::ratioGraphs(gconf90,gconf68);
    TGraphErrors *grconf98 = tools::ratioGraphs(gconf98,gconf68);
    TGraphErrors *grrms = tools::ratioGraphs(grms,gconf68);
    TGraphErrors *grsigma = tools::ratioGraphs(gsigma,gconf68);
    TGraphErrors *grsigmap = tools::ratioGraphs(gsigmap,gconf68);
    TGraphErrors *grsigman = tools::ratioGraphs(gsigman,gconf68);

    tdrDraw(grconf68,"L0",kNone,kGray+1);
    tdrDraw(grconf90,"L0",kNone,kGray+2);
    tdrDraw(grconf98,"L0",kNone,kGray+3);
    tdrDraw(grrms,"L",kNone,kGray+2);
    tdrDraw(grsigma,"L",kNone,kCyan+1);
    tdrDraw(grsigmap,"L",kNone,kBlue);
    tdrDraw(grsigman,"L",kNone,kRed);

    cout << endl << flush;
    if (fabs(eta1-0)<0.1) {
      c0lin->SaveAs("pdf/henning/henning_linfits.pdf");
      c0log->SaveAs("pdf/henning/henning_logfits.pdf");
      c1->SaveAs("pdf/henning/henning_resp.pdf");
      c4->SaveAs("pdf/henning/henning_effvspT.pdf");
      c2->SaveAs("pdf/henning/henning_jesvspT.pdf");
      c3->SaveAs("pdf/henning/henning_jervspT.pdf");
    }
    string sy = Form("eta_%04.0f_%05.0f",1000*eta1,1000*eta2);
    const char *cy = sy.c_str();
    c0lin->SaveAs(Form("pdf/henning/bins/henning_linfits_%s.pdf",cy));
    c0log->SaveAs(Form("pdf/henning/bins/henning_logfits_%s.pdf",cy));
    c1->SaveAs(Form("pdf/henning/bins/henning_resp_%s.pdf",cy));
    c4->SaveAs(Form("pdf/henning/bins/henning_effvspT_%s.pdf",cy));
    c2->SaveAs(Form("pdf/henning/bins/henning_jesvspT_%s.pdf",cy));
    c3->SaveAs(Form("pdf/henning/bins/henning_jervspT_%s.pdf",cy));
  } // for ieta

  TFile *fout = new TFile("rootfiles/henning.root","RECREATE");
  h2eff_any->Write("h2eff_any");
  h2eff_pt8->Write("h2eff_pt8");
  h2eff_core->Write("h2eff_core");

  h2jes_med->Write("h2jes_med");
  h2jes_mean->Write("h2jes_mean");
  h2jes_gaus->Write("h2jes_gaus");

  h2jer_qrt->Write("h2jer_qrt");
  h2jer_rms->Write("h2jer_rms");
  h2jer_gaus->Write("h2jer_gaus");

  fout->Write();
  fout->Close();
  curdir->cd();
  
} // henning
