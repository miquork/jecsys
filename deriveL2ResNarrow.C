// Attempt to derive L2Res from single jet triggers
// as a by-product on inclusive jet cross section measurement
// Inputs from Hannu Siikonen
#include "TFile.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TMatrixD.h"
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TMath.h"

#include "tdrstyle_mod15.C"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

struct l2fit {

  double chi2;
  int ndf;
  TF1 *f1;
  TMatrixD *emat;
  
  l2fit(double chi2_ = 0, int ndf_ = 0, TF1 *f1_ = 0, TMatrixD *emat_ = 0) :
    chi2(chi2_), ndf(ndf_), f1(f1_), emat(emat_) { };
};

const double _xmin(30.), _xmax(3500.);

FactorizedJetCorrector *_jer_dt(0), *_jer_mc(0);
Double_t fJER(Double_t *x, Double_t *p) {
  // Note: still need to implement quadratic averaging over
  //       eta and rho bins for |eta|<1.3 and all rho
  //       Where is the sqrt(2) compared to sigma(asymmetry)?

  double pt = x[0];
  double eta = p[0];
  bool ismc = (p[1]==0);
  double rho = 25.*1.0;//p[1];

  if (!_jer_dt) {
    const char *sd = "CondFormats/JetMETObjects/data";
    const char *st = "Spring16_25nsV6_DATA";
    const char *sa = "PtResolution_AK4PFchs";
    const char *s = Form("%s/%s_%s.txt",sd,st,sa); cout << s << endl << flush;

    JetCorrectorParameters *pjer = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*pjer);
    _jer_dt = new FactorizedJetCorrector(v);
  }
  if (!_jer_mc) {
    const char *sd = "CondFormats/JetMETObjects/data";
    const char *st = "Spring16_25nsV6_MC";
    const char *sa = "PtResolution_AK4PFchs";
    const char *s = Form("%s/%s_%s.txt",sd,st,sa); cout << s << endl << flush;

    JetCorrectorParameters *pjer = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*pjer);
    _jer_mc = new FactorizedJetCorrector(v);
  }

  FactorizedJetCorrector *_jer = (ismc ? _jer_mc : _jer_dt);
  _jer->setJetPt(pt);
  _jer->setJetEta(eta);
  _jer->setRho(rho);
  
  double jer = _jer->getCorrection();
  
  return jer / sqrt(2.);
  //return jer;
} // fJER

double _eta(2.0);
Double_t fXsec(Double_t *x, Double_t *p) {

  double pt = x[0];
  double xsec = p[0]*pow(pt,p[1])*pow(1-2*pt*cosh(_eta)/13000.,p[2]);

  return xsec;
}

TF1 *_xsec(0);
Double_t deltaPtKernel2(Double_t *x, Double_t *p) {

  double dpt = x[0];
  double pTave = p[0];
  double sigma1 = p[1];
  double sigma2 = p[2];
  double pTgen = p[3];
  bool norm = (p[4]==0);

  if (!_xsec) {
    _xsec = new TF1("xsec",fXsec,_xmin,_xmax,3);
    //_xsec->SetParameters(1.03331e+09, -4.27019e+00, 6.92540e+00); // pT,ave
    _xsec->SetParameters(1.42389e+09, -4.29955e+00, 6.22351e+00); // pT,tag
  }

  // dpt = (pt2-pT1)/(2*pTave)
  double s1 = TMath::Gaus(pTave-dpt*pTave-pTgen, 0, sigma1*pTave, kTRUE);
  double s2 = TMath::Gaus(pTave+dpt*pTave-pTgen, 0, sigma2*pTave, kTRUE);
  double xsec = _xsec->Eval(pTgen);
  // Could the cross section fall steeper vs jet2 pT?
  //xsec /= (1 - 2.*pTgen*cosh(_eta)/13000.,6.9254);
  //xsec *= (1 - 2.*pTge*cosh(_eta)/13000.,6.9254);

  return (norm ? s1*s2*xsec : dpt*s1*s2*xsec);
  //return (norm ? s1*s2*xsec : fabs(dpt)*s1*s2*xsec); // does show large effect
} // deltaPtKernel2


TF1 *_fdpt2(0); // Integral over dpt
Double_t deltaPtKernel1(Double_t *x, Double_t *p) {

  double pTgen = x[0];
  double pTave = p[0];
  double sigma1 = p[1];
  double sigma2 = p[2];
  double norm = p[3];

  if (!_fdpt2)
    _fdpt2 = new TF1("fdpt2",deltaPtKernel2,_xmin,_xmax,5);

  double ksigma = 2.5;//3.5;
  double sigma = sqrt(sigma1*sigma1 + sigma2*sigma2);
  _fdpt2->SetParameters(pTave, sigma1, sigma2, pTgen, norm);
  double Xdn = _fdpt2->Integral(-ksigma*sigma, +ksigma*sigma, 1e-4);

  return Xdn; // dpt*dn or just dn
} // deltaPtKernel1

// Calculate deltaPt given pT,ave spectrum and JER
// Assume for simplicity pTgen1=pTgen=pTgen => pTave,gen=pTgen
// deltaPt = (pTreco2 - pTreco1) / 2
// double dpt = \integral_pTgen \integral_dpt dpt
//             * JER(pTave+dpt-pTgen; sigma2) * JER(pTave-dpt-pTgen; sigma1)
//             * f(pTgen) d(dpt) d(pTgen)
TF1 *_fdpt1(0); // Integral over pTgen
Double_t deltaPt(Double_t *x, Double_t *p) {

  double pTave = x[0];
  double eta1 = p[0];
  double eta2 = p[1];
  double isdt = p[2]; // 0==mc, 1==data
  double p1[2] = {eta1,isdt};
  double p2[2] = {eta2,isdt};

  double ksigma = 2.5;//3.5;
  double sigma1 = sqrt(2.)*fJER(&pTave, &p1[0]);
  double sigma2 = sqrt(2.)*fJER(&pTave, &p2[0]);

  if (!_fdpt1)
    _fdpt1 = new TF1("fdpt1",deltaPtKernel1,_xmin,_xmax,4);

  double sigma = sqrt(sigma1*sigma1 + sigma2*sigma2);
  _fdpt1->SetParameters(pTave, sigma1, sigma2, 1);
  double dptdn = _fdpt1->Integral(pTave*(1-ksigma*sigma),
				  pTave*(1+ksigma*sigma),1e-4);
  _fdpt1->SetParameters(pTave, sigma1, sigma2, 0);
  double dn = _fdpt1->Integral(pTave*(1-ksigma*sigma),
			       pTave*(1+ksigma*sigma),1e-4);

  return (dptdn / dn);
} // deltaPt

double getErr(double x, TF1 *f1, TMatrixD *emat) {

  int n = f1->GetNpar();
  assert(emat->GetNcols()==n);
  assert(emat->GetNrows()==n);
  vector<double> df(n);
  for (int i = 0; i != n; ++i) {
    double p = f1->GetParameter(i);
    double dp = f1->GetParError(i);
    assert(dp!=0);
    double y1 = f1->Eval(x);
    f1->SetParameter(i, p+0.01*dp);
    double y2 = f1->Eval(x);
    f1->SetParameter(i, p);
    df[i] = (y2-y1)/(0.01*dp);
  }

  double err2(0);
  for (int i = 0; i != n; ++i) {
    for (int j = 0; j != n; ++j) {
      err2 += df[i]*df[j]*(*emat)[i][j];
    }
  } // for i

  return sqrt(err2);
}

// Derive truncated RMS (and it's uncertainty) from histogram
// by finding smallest range that contains frac events
//double truncRMS(TH1D *h, double *err, double frac = 0.985) {
double truncRMS(TH1D *h, double frac = 0.985) {
  //double truncRMS(TH1D *h, double frac = 0.975) {
  
  //return h->GetRMS(); // TMP speedup
 
  TH1D *htmp = (TH1D*)h->Clone("htmp");
  //htmp->Clear();

  //double sum = h->Integral();
  double sum = h->Integral(1,h->GetNbinsX()); // avoid overflows
  int imean = h->FindBin(h->GetMean());
  int imin(-1), jmin(h->GetNbinsX()+1);
  double dxmin(h->GetNbinsX()+1);
  for (int i = 1; i != imean; ++i) {
    //if (h->Integral(i, h->GetNbinsX()+1) < frac) continue;
    for (int j = h->GetNbinsX(); j != imean; --j) {
      int dx = j-i;
      //if (h->Integral(i,j) / sum < frac) continue;
      if (h->Integral(i,j) / sum >= frac && dx < dxmin) {
	dxmin = dx;
	imin = i;
	jmin = j;
      }  
    } // for j
  } // for i

  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    if (i<imin || i > jmin) {
      htmp->SetBinContent(i, 0);
      htmp->SetBinError(i, 0);
    }
  } // for i

  double rms = htmp->GetRMS();
  //*err = htmp->GetRMS()/sqrt(h->GetEffectiveEntries());
  delete htmp;
  
  return rms; 
} // truncRMS


l2fit deriveL2ResNarrow_x(double y1=2.0, double y2=2.5,
		    double alphamax=0.2, double xmax=2000.);

TH1D *_htrg(0);
string stp = "";//"tp";
const char *ctp = stp.c_str();//"tp";
void deriveL2ResNarrow() {

  // For spectrum reconstruction
  const int ntrg = 9;   // 1st 74 (jet40) 
  //double trigpt[ntrg+1] = {0, 84, 114, 196, 272, 330, 395, 468, 548,3500};
  double trigpt[ntrg+1] = {0, 84, 114, 174, 245, 330, 395, 468, 548,3500};
  // Effective prescales for RunFLateG_rereco from Engin root tuple
  //152728, 152655, 152261, 151443, 149893, 146976, 141570, 129504, 20007, 19543.2, 18588.5, 16868.2, 1986.39, 1902.99, 1758.02, 564.405, 511.402, 66.7384, 61.396, 22.0983, 20.2887, 7.31431, 6.60524,
  //double trigpre[ntrg] = {152728, 129504, 20007, 1986, 564.4, 66.74, 22.10, 7.314, 1};
  //double trigpre[ntrg] = {152728, 50000, 20007, 1986, 564.4, 66.74, 22.10, 7.314, 1};
  double trigpre[ntrg] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
  _htrg = new TH1D("htrg",";p_{T};prescale;",ntrg,trigpt);
  for (int i = 1; i != _htrg->GetNbinsX()+1; ++i) {
    _htrg->SetBinContent(i, trigpre[i-1]);
  }

  //double etas[] = {0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0,
  //	   3.2, 4.7};
  //double xmax[] = {2000, 2000, 2000, 1500, 1200, 800, 500, 350};
  double etas[] = {0,0.261,0.522,0.783,1.044,1.305,
		   1.479,1.653,1.93,2.172,2.322,2.5,
		   2.65,2.853,2.964,3.139,3.489,3.839,
		   5.191};
  double xmax[] = {2000, 2000, 2000, 2000, 2000, 1780,
		   1600, 1400, 1100, 1000, 900, 800,
		   680, 630, 420, 290, 190};
  // Debug one bin first
//   double etas[] = {2.5,3.0};//2.65};
//   double xmax[] = {800};
  const int neta = sizeof(etas)/sizeof(etas[0]) - 1;
  //double alphamax = 0.3;
  //double alphamax = 0.2;
  //double alphamax = 0.1;
  //double alphamax = 0.05;
  //double alphamax = 0.025;
  //double alphas[] = {0.025, 0.05, 0.1, 0.2, 0.3};
  //double alphas[] = {0.2};
  //double alphas[] = {0.05};
  double alphas[] = {0.10};
  const int nalpha = sizeof(alphas)/sizeof(alphas[0]);

  double pts[] = {60, 120, 240, 480};
  int colorpt[] = {kBlack, kRed, kGreen+2, kBlue};
  int markerpt[] = {0, kFullSquare, 0, 0};
  const int npt = sizeof(pts)/sizeof(pts[0]);
  vector<TH1D*> hpts(npt);
  for (int ipt = 0; ipt != npt; ++ipt) {
    hpts[ipt] = new TH1D(Form("hpt_%1.0f",pts[ipt]),
			 ";|#eta|;Relative correction",
			 neta, &etas[0]);
  }

  for (int ialpha = 0; ialpha != nalpha; ++ialpha) {

    const double alphamax = alphas[ialpha];

    vector<l2fit> l2s(neta);
    TGraphErrors *gchi2 = new TGraphErrors(neta);
    for (int ieta = 0; ieta != neta; ++ieta) {
      
      double eta = 0.5*(etas[ieta] + etas[ieta+1]);
      double deta = 0.5*(etas[ieta+1] - etas[ieta]);
      l2s[ieta] = deriveL2ResNarrow_x(etas[ieta], etas[ieta+1],
				      alphamax, xmax[ieta]);
      
      gchi2->SetPoint(ieta, eta, l2s[ieta].chi2 / l2s[ieta].ndf);
      gchi2->SetPointError(ieta, deta, 1. / sqrt(l2s[ieta].ndf));
      
      for (int ipt = 0; ipt != npt; ++ipt) {
	double pt = pts[ipt];
	double B = l2s[ieta].f1->Eval(pt);
	double Rrel = (1 + B) / (1 - B);
	if (pt*cosh(eta) < 6500.) {
	  hpts[ipt]->SetBinContent(ieta+1, 1./Rrel);
	  if (markerpt[ipt]!=0) {
	    double err = getErr(pt, l2s[ieta].f1, l2s[ieta].emat);
	    hpts[ipt]->SetBinError(ieta+1, err);
	  }
	}
      }
    } // for ieta
    
    TH1D *h1 = new TH1D("h1",";#eta;#chi^{2} / NDF",neta,&etas[0]);
    h1->SetMinimum(0.);
    h1->SetMaximum(3);
    
    TCanvas *c1 = tdrCanvas("c1",h1,4,11,kSquare);
    
    TLine *l1 = new TLine();
    l1->SetLineStyle(kDashed);
    l1->DrawLine(etas[0],1,etas[neta],1);
    
    tdrDraw(gchi2,"Pz",kFullCircle);
    
    c1->SaveAs(Form("pdf/deriveL2ResNarrow_chi2_amax%1.0f%s.pdf",
		    alphamax*100,ctp));
    
    
    
    TH1D *h2 = new TH1D("h2",";#eta;Relative correction",neta,&etas[0]);
    h2->SetMinimum(0.8);//0.98 - 0.03 - 0.05);
    h2->SetMaximum(1.3);//1.19-1e-4 - 0.03 + 0.05);
    
    TCanvas *c2 = tdrCanvas("c2",h2,4,11,kSquare);
    
    TLine *l2 = new TLine();
    l2->SetLineStyle(kDashed);
    l2->DrawLine(etas[0],1,etas[neta],1);
    
    TLegend *leg2 = tdrLeg(0.17,0.48,0.47,0.78);
    leg2->SetHeader("p_{T}(jet)");
    
    for (int ipt = 0; ipt != npt; ++ipt) {
      tdrDraw(hpts[ipt], markerpt[ipt]!=0 ? "HP][" : "H][", markerpt[ipt],
	      colorpt[ipt],kSolid,colorpt[ipt],kNone,0);
      hpts[ipt]->SetLineWidth(2);
      leg2->AddEntry(hpts[ipt],Form("%1.0f GeV",pts[ipt]),
		     markerpt[ipt]!=0 ? "LP" : "L");
    }
    
    gPad->RedrawAxis();
    c2->SaveAs(Form("pdf/deriveL2ResNarrow_Fig15_amax%1.0f%s.pdf",
		    alphamax*100,ctp));

    if (ialpha!=nalpha-1) { // leave last ones open
      delete h1;
      delete h2;
    }
  } // for ialpha
}


l2fit deriveL2ResNarrow_x(double y1, double y2, double alphamax, double xmax) {

  TDirectory *curdir = gDirectory;
  TFile *fd = new TFile("rootfiles/output_RunGV6x/output-DATA-2b.root","READ");
  //TFile *fd = new TFile("rootfiles/output_RunGV6/output-DATA-2b.root","READ");
  //TFile *fd = new TFile("rootfiles/output_RunGV6/output-DATA-2b-noL2Res.root",
  //		"READ");
  //TFile *fd = new TFile("rootfiles/output_RunGV6/output-DATA-2b.root","READ");
  //TFile *fd = new TFile("rootfiles/output_RunGV5/output-DATA-2b.root","READ");
  //TFile *fd = new TFile("rootfiles/output_RunGV3/output-DATA-2b.root","READ");
  //TFile *fd = new TFile("rootfiles/output_RunG_old/output-DATA-2b.root","READ");
  assert(fd && !fd->IsZombie());
  TFile *fm = new TFile("rootfiles/output_RunGV6x/output-MC-2b.root","READ");
  //TFile *fm = new TFile("rootfiles/output_RunGV6/output-MC-2b.root","READ");
  //TFile *fm = new TFile("rootfiles/output_RunGV5/output-MC-2b.root","READ");
  //TFile *fm = new TFile("rootfiles/output_RunGV3/output-MC-2b.root","READ");
  //TFile *fm = new TFile("rootfiles/output_RunG_old/output-MC-2b.root","READ");
  assert(fm && !fm->IsZombie());

  TFile *fz = new TFile("rootfiles/L3res_mumu_Run2016FlateG.root","READ");
  assert(fz && !fz->IsZombie());
  curdir->cd();

// See: (pt-alpha-asymm/mpf)
//  KEY: TH3D      hdjasymm;1    => pT,ave
//  KEY: TH3D      hdjasymmtp;1  => pT,tag
//  KEY: TH3D      hdjmpf;1
//  KEY: TH3D      hdjmpftp;1

// Note: mpf here means the same as balance/B.

  //const char *cd = "Standard";
  //const char *ce0 = "Eta_0.0-1.3"; // Does this make a difference? Event counts?
  const char *cd = "FullEta";
  const char *ce0 = "";
  //string se = Form("Eta_%1.1f-%1.1f",y1,y2);
  string se = Form("Eta_%1.3f-%1.3f",y1,y2);
  const char *ce = se.c_str();
  string se2 = Form("Eta_%02.0f-%02.0f",10*y1,10*y2);
  const char *ce2 = se2.c_str();
  double xmin(30.);//60.);//, xmax(2000.);
  double fitxmin(64.);
  const char *ca = "_a01";//Form("_a%");

  setTDRStyle();

  // Tag-and-probe with pT,tag binning: ctp = "tp"
  // Regular asymmetry with pT,ave binning: ctp = ""
  //TH3D *h3da = (TH3D*)fd->Get(Form("%s/%s/hdjasymm%s",cd,ce,ctp)); assert(h3da);
  //TH3D *h3db = (TH3D*)fd->Get(Form("%s/%s/hdjmpf%s",cd,ce,ctp));   assert(h3db);
  //TH3D *h3ma = (TH3D*)fm->Get(Form("%s/%s/hdjasymm%s",cd,ce,ctp)); assert(h3ma);
  //TH3D *h3mb = (TH3D*)fm->Get(Form("%s/%s/hdjmpf%s",cd,ce,ctp));   assert(h3mb);

  TH3D *h3da = (TH3D*)fd->Get(Form("%s/%s/hdjasymm%s%s",cd,ce0,ctp,ca));
  TH3D *h3db = (TH3D*)fd->Get(Form("%s/%s/hdjmpf%s%s",cd,ce0,ctp,ca));
  TH3D *h3ma = (TH3D*)fm->Get(Form("%s/%s/hdjasymm%s%s",cd,ce0,ctp,ca));
  TH3D *h3mb = (TH3D*)fm->Get(Form("%s/%s/hdjmpf%s%s",cd,ce0,ctp,ca));
  assert(h3da);
  assert(h3db);
  assert(h3ma);
  assert(h3mb);


  // First, look at mean and pT,ave spectrum
  //int ymin = h3da->GetYaxis()->FindBin(0.);
  //int ymax = h3da->GetYaxis()->FindBin(alphamax-1e-4);
  int ymin = h3da->GetYaxis()->FindBin(y1);
  int ymax = h3da->GetYaxis()->FindBin(y2-1e-4);
  h3da->Sumw2();
  h3da->GetYaxis()->SetRange(ymin,ymax);
  TH2D *h2da = (TH2D*)h3da->Project3D("zxe"); // X=pT,ave; Y=Asymmetry
  h2da->Sumw2();
  h2da->SetName("h2da");
  TProfile *p1da = h2da->ProfileX("p1da",1,-1,"e");
  if (stp=="tp") p1da->Scale(0.5);
  TH1D *hnda = h2da->ProjectionX("hnda",0,-1,"e");
  //
  h3ma->Sumw2();
  h3ma->GetYaxis()->SetRange(ymin,ymax);
  TH2D *h2ma = (TH2D*)h3ma->Project3D("zxe"); // X=pT,ave; Y=Asymmetry
  h2ma->Sumw2();
  h2ma->SetName("h2ma");
  TProfile *p1ma = h2ma->ProfileX("p1ma",1,-1,"e");
  if (stp=="tp") p1ma->Scale(0.5);
  TH1D *hnma = h2ma->ProjectionX("hnma",0,-1,"e");

  h3db->Sumw2();
  h3db->GetYaxis()->SetRange(ymin,ymax);
  TH2D *h2db = (TH2D*)h3db->Project3D("zxe"); // X=pT,ave; Y=Asymmetry
  h2db->Sumw2();
  h2db->SetName("h2db");
  TProfile *p1db = h2db->ProfileX("p1db",1,-1,"e");
  p1db->Scale(0.5);
  //
  h3mb->Sumw2();
  h3mb->GetYaxis()->SetRange(ymin,ymax);
  TH2D *h2mb = (TH2D*)h3mb->Project3D("zxe"); // X=pT,ave; Y=Asymmetry
  h2mb->Sumw2();
  h2mb->SetName("h2mb");
  TProfile *p1mb = h2mb->ProfileX("p1mb",1,-1,"e");
  p1mb->Scale(0.5);

  TH1D *h1ra = p1da->ProjectionX("h1ra","e");
  h1ra->Add(p1ma,-1);
  TH1D *h1rb = p1db->ProjectionX("h1rb","e");
  h1rb->Add(p1mb,-1);

  // Second, look at RMS / Gaussian sigma
  TF1 *fga = new TF1("fga","gaus",-1,1);
  TF1 *fgb = new TF1("fgb","gaus",-1,1);
  // fit range [mean-ksigma*sigma, mean+ksigma*sigma]
  const double ksigma = 1.5; 

  TH1D *h1da = p1da->ProjectionX("h1da","e"); h1da->Clear();
  TH1D *h1ma = p1ma->ProjectionX("h1ma","e"); h1ma->Clear();
  TH1D *h1db = p1db->ProjectionX("h1db","e"); h1db->Clear();
  TH1D *h1mb = p1mb->ProjectionX("h1mb","e"); h1mb->Clear();
  for (int i = 1; i != h1da->GetNbinsX()+1; ++i) {
    TH1D *htmpa = h2da->ProjectionY("htmpa",i,i,"e");
    TH1D *htmpb = h2db->ProjectionY("htmpb",i,i,"e");
    double rmsa = truncRMS(htmpa);
    double rmsb = truncRMS(htmpb);
    fga->SetRange(htmpa->GetMean() - rmsa*ksigma,
		  htmpa->GetMean() + rmsa*ksigma);
    fgb->SetRange(htmpb->GetMean() - rmsb*ksigma,
		  htmpb->GetMean() + rmsb*ksigma);
    fga->SetParameters(htmpa->Integral(),htmpa->GetMean(),htmpa->GetRMS());
    fgb->SetParameters(htmpb->Integral(),htmpb->GetMean(),htmpb->GetRMS());
    htmpa->Fit(fga,"QRN");
    htmpb->Fit(fgb,"QRN");
    if (htmpa->Integral()!=0.) {
      h1da->SetBinContent(i, (stp=="tp" ? 0.5 : 1)*rmsa);
      h1da->SetBinError(i, (stp=="tp" ? 0.5 : 1)*htmpa->GetRMS()/sqrt(htmpa->GetEffectiveEntries()));
      h1db->SetBinContent(i, 0.5*rmsb);
      h1db->SetBinError(i, 0.5*htmpb->GetRMS()/sqrt(htmpb->GetEffectiveEntries()));
    }
    delete htmpa;
    delete htmpb;
    // Use data fit as starting value for MC
    htmpa = h2ma->ProjectionY("htmpa",i,i,"e");
    htmpb = h2mb->ProjectionY("htmpb",i,i,"e");
    rmsa = truncRMS(htmpa);
    rmsb = truncRMS(htmpb);
    fga->SetRange(htmpa->GetMean() - rmsa*ksigma,
		  htmpa->GetMean() + rmsa*ksigma);
    fgb->SetRange(htmpb->GetMean() - rmsb*ksigma,
		  htmpb->GetMean() + rmsb*ksigma);
    htmpa->Fit(fga,"QRN");
    htmpb->Fit(fgb,"QRN");
    if (htmpa->Integral()!=0.) {
      h1ma->SetBinContent(i, (stp=="tp" ? 0.5 : 1)*rmsa);
      h1ma->SetBinError(i, (stp=="tp" ? 0.5 : 1)*htmpa->GetRMS()/sqrt(htmpa->GetEffectiveEntries())); 
      h1mb->SetBinContent(i, 0.5*rmsb);
      h1mb->SetBinError(i, 0.5*htmpb->GetRMS()/sqrt(htmpb->GetEffectiveEntries())); 
    }
    delete htmpa;
    delete htmpb;
  }

  curdir->cd();

  TH1D *hup = new TH1D("hup",";p_{T,ave} (GeV);#LTAsymmetry#GT",
		       //1940,60,2000);
		       1970,30,2000);
  hup->SetMaximum(+0.205);//+0.155);
  hup->SetMinimum(-0.155);//-0.105);
  TH1D *hdw = new TH1D("hdw",";p_{T,ave} (GeV);Data-MC",
		       //1940,60,2000);
		       1930,30,2000);
  hdw->SetMaximum(+0.105);//+0.055);
  hdw->SetMinimum(-0.105);//-0.055);
  hdw->GetXaxis()->SetMoreLogLabels();
  hdw->GetXaxis()->SetNoExponent();
  if (stp=="tp") hdw->SetXTitle("p_{T,tag}");


  lumi_13TeV = "RunG";
  TCanvas *c2 = tdrDiCanvas("c2",hup,hdw,4,11);

  c2->cd(1);
  gPad->SetLogx();
  p1da->GetXaxis()->SetRangeUser(xmin,xmax);
  p1ma->GetXaxis()->SetRangeUser(xmin,xmax);
  p1db->GetXaxis()->SetRangeUser(xmin,xmax);
  p1mb->GetXaxis()->SetRangeUser(xmin,xmax);
  tdrDraw(p1db,"P",kFullSquare,kRed);
  tdrDraw(p1mb,"P",kOpenSquare,kRed);
  tdrDraw(p1da,"P",kFullCircle,kBlue);
  tdrDraw(p1ma,"P",kOpenCircle,kBlue);
  
  TF1 *f1da = new TF1("f1da","[0]+[1]*log(x)+[2]*log(x)*log(x)",fitxmin,xmax);
  f1da->SetParameters(0,0.01,0.0001);
  p1da->Fit(f1da,"QRN");
  f1da->SetLineColor(kBlue);
  f1da->Draw("SAME");
  TF1 *f1ma = new TF1("f1ma","[0]+[1]*log(x)+[2]*log(x)*log(x)",fitxmin,xmax);
  f1ma->SetParameters(0,0.01,0.0001);
  p1ma->Fit(f1ma,"QRN");
  f1ma->SetLineColor(kBlue);
  f1ma->SetLineStyle(kDashed);
  f1ma->Draw("SAME");
  //
  TF1 *f1db = new TF1("f1db","[0]+[1]*log(x)+[2]*log(x)*log(x)",fitxmin,xmax);
  f1db->SetParameters(0,0.01,0.0001);
  p1db->Fit(f1db,"QRN");
  f1db->SetLineColor(kRed);
  f1db->Draw("SAME");
  TF1 *f1mb = new TF1("f1mb","[0]+[1]*log(x)+[2]*log(x)*log(x)",fitxmin,xmax);
  f1mb->SetParameters(0,0.01,0.0001);
  p1mb->Fit(f1mb,"QRN");
  f1mb->SetLineColor(kRed);
  f1mb->SetLineStyle(kDashed);
  f1mb->Draw("SAME");

  // Check estimated resolution bias
  TF1 *fddpt = new TF1("fddpt",deltaPt,fitxmin,xmax,3);
  fddpt->SetParameters(0.65,0.5*(y1+y2), 1);//0.0, 2.5, 1);
  fddpt->SetLineStyle(kDotted);
  fddpt->SetLineColor(kGreen+2);
  fddpt->SetLineWidth(2);
  fddpt->Draw("SAME");
  TF1 *fmdpt = new TF1("fmdpt",deltaPt,fitxmin,xmax,3);
  fmdpt->SetParameters(0.65, 0.5*(y1+y2), 0);//0.0, 2.5, 0);
  fmdpt->SetLineStyle(kDotted);
  fmdpt->SetLineColor(kMagenta+1);
  fmdpt->SetLineWidth(2);
  fmdpt->Draw("SAME");
  
  TLatex *tex = new TLatex();
  tex->SetTextSize(0.045);
  tex->SetNDC();
  tex->DrawLatex(0.40,0.85,ce);
  tex->DrawLatex(0.40,0.75,Form("#alpha < %1.2f",alphamax));
  if (stp=="tp")  tex->DrawLatex(0.60,0.75,"(tp)");

  TLegend *leg2up = tdrLeg(0.7,0.66,0.9,0.88);
  leg2up->AddEntry(p1db,"#LTB#GT_{data}","PL");
  leg2up->AddEntry(p1mb,"#LTB#GT_{MC}","PL");
  leg2up->AddEntry(p1da,"#LTA#GT_{data}","PL");
  leg2up->AddEntry(p1ma,"#LTA#GT_{MC}","PL");


  c2->cd(2);
  gPad->SetLogx();
  h1ra->GetXaxis()->SetRangeUser(xmin,xmax);
  h1rb->GetXaxis()->SetRangeUser(xmin,xmax);
  tdrDraw(h1rb,"P",kFullSquare,kRed);
  tdrDraw(h1ra,"P",kFullCircle,kBlue);

  TF1 *f2a = new TF1("f2a","[0]+[1]*log(x)+[2]*log(x)*log(x)",fitxmin,xmax);
  f2a->SetParameters(0,0.01,0.0001);
  h1ra->Fit(f2a,"QRN");
  f2a->SetLineColor(kBlue);
  f2a->Draw("SAME");
  TF1 *f2b = new TF1("f2b","[0]+[1]*log(x)+[2]*log(x)*log(x)",fitxmin,xmax);
  f2b->SetParameters(0,0.01,0.0001);
  h1rb->Fit(f2b,"QRN");
  f2b->SetLineColor(kRed);
  f2b->Draw("SAME");

  TMatrixD *emat = new TMatrixD(f2b->GetNpar(), f2b->GetNpar());
  gMinuit->mnemat(emat->GetMatrixArray(), f2b->GetNpar());

  TLegend *leg2dw = tdrLeg(0.6,0.8,0.9,0.9);
  leg2dw->SetTextSize(2.*0.045);
  leg2dw->SetNColumns(2);
  leg2dw->AddEntry(h1rb,"#LTB#GT","PL");
  leg2dw->AddEntry(h1ra,"#LTA#GT","PL");
  
  // Add Z+jet
  string sz0 = Form("Ratio_MPF_CHS_a30_eta_%04.0f_%04.0f_L1L2L3",
		   1000.*0,1000.*1.3);
  cout << sz0 << endl << flush;
  TGraphErrors *gz0 = (TGraphErrors*)fz->Get(sz0.c_str());
  assert(gz0);
  string sz = Form("Ratio_MPF_CHS_a30_eta_%04.0f_%04.0f_L1L2L3",
		   1000.*y1,1000.*y2);
  cout << sz << endl << flush;
  TGraphErrors *gz = (TGraphErrors*)fz->Get(sz.c_str());
  assert(gz);

  for (int i = 0; i != gz->GetN(); ++i) {

    double x = gz->GetX()[i];
    double y = gz->GetY()[i];

    int k(-1);
    for (int j = 0; j != gz->GetN(); ++j) {
      if (fabs(gz->GetX()[j] - x) < 0.03*x) {
	assert(k==-1);
	k = j;
      }
    } // for j
    assert(k!=-1);
    double y0 = gz0->GetY()[i];

    gz->SetPoint(i, x, y / y0 - 1);
  } // for i

  TF1 *f2az = new TF1("f2az","[0]+[1]*log(x)+[2]*log(x)*log(x)",fitxmin,xmax);
  f2az->SetParameters(f2a->GetParameter(0),f2a->GetParameter(1),
		      f2a->GetParameter(2));
  f2az->SetLineStyle(kDotted);
  f2az->SetLineColor(kBlue);
  f2az->Draw("SAME");
  TF1 *f2bz = new TF1("f2bz","[0]+[1]*log(x)+[2]*log(x)*log(x)",fitxmin,xmax);
  f2bz->SetParameters(f2b->GetParameter(0),f2b->GetParameter(1),
		      f2b->GetParameter(2));
  f2bz->SetLineStyle(kDotted);
  f2bz->SetLineColor(kRed);
  f2bz->Draw("SAME");

  tdrDraw(gz,"Pz",kOpenSquare,kGreen+2);


  c2->SaveAs(Form("pdf/deriveL2ResNarrow_%s_amax%1.0f%s.pdf",
		  ce2,alphamax*100,ctp));


  TH1D *hup3 = new TH1D("hup3",";p_{T,ave} (GeV);#sigma(asymmetry)",
			1940,60,2000);
  hup3->SetMaximum(0.3);//0.5);
  hup3->SetMinimum(0.+1e-4);
  TH1D *hdw3 = new TH1D("hdw3",";p_{T,ave} (GeV);Data/MC",
		       1940,60,2000);
  hdw3->SetMaximum(1+0.5-1e-4);
  hdw3->SetMinimum(1-0.3+1e-4);
  hdw3->GetXaxis()->SetMoreLogLabels();
  hdw3->GetXaxis()->SetNoExponent();
  if (stp=="tp") hdw3->SetXTitle("p_{T,tag}");

  TCanvas *c3 = tdrDiCanvas("c3",hup3,hdw3,4,11);

  TF1 *fsdb = new TF1("fsdb","sqrt([0]*[0]+[1]*[1]/x+[2]*[2]/(x*x))",xmin,xmax);
  fsdb->SetParameters(0.038, 0.9, 10.4);
  h1db->Fit(fsdb,"RN");
  TF1 *fsmb = new TF1("fsmb","sqrt([0]*[0]+[1]*[1]/x+[2]*[2]/(x*x))",xmin,xmax);
  fsmb->SetParameters(0.028, 0.9, 10.4);
  h1mb->Fit(fsmb,"QRN");

  TF1 *fsda = new TF1("fsda","sqrt([0]*[0]+[1]*[1]/x+[2]*[2]/(x*x)+"
		      "[3]*[3]*pow(x,2*[4]))",xmin,xmax);
  fsda->SetParameters(0.038, 0.9, 4.1, 0.46*alphamax,-0.29);
  fsda->FixParameter(0, fsdb->GetParameter(0));
  fsda->FixParameter(1, fsdb->GetParameter(1));
  //fsda->FixParameter(2, fsdb->GetParameter(2));
  // https://twiki.cern.ch/twiki/pub/CMSPublic/MultipleConeSizes14/jerplots_ParFits_Combined.pdf (muxA = 25*pi*0.5^2 = 12.5 for 2016, N~2.5);
  double parN = sqrt(-1*fabs(-1) + 0.75*0.75*25*TMath::Pi()*0.4*0.4);
  fsda->FixParameter(2, parN);
  h1da->Fit(fsda,"RN");
  TF1 *fsma = new TF1("fsma","sqrt([0]*[0]+[1]*[1]/x+[2]*[2]/(x*x)+"
		      "[3]*[3]*pow(x,2*[4]))",xmin,xmax);
  //"[3]*[3]*exp(x*2*[4]))",xmin,xmax);
  fsma->SetParameters(0.038, 0.9, 4.1, 0.46*alphamax,-0.29);
  fsma->FixParameter(0, fsmb->GetParameter(0));
  fsma->FixParameter(1, fsmb->GetParameter(1));
  //fsma->FixParameter(2, fsmb->GetParameter(2));
  fsma->FixParameter(2, parN);
  h1ma->Fit(fsma,"QRN");

  // Note: PLI model (fit vs pTave and alphamax) is not yet ideal
  //       It's not fully describing the low pTave end of asymmetry
  //       However, could be that alphamax >= 24 GeV/pTave is kicking in
  //       for low values of pTave and alphamax, confusing the picture


  c3->cd(1);
  gPad->SetLogx();
  h1da->GetXaxis()->SetRangeUser(xmin,xmax);
  h1ma->GetXaxis()->SetRangeUser(xmin,xmax);
  h1db->GetXaxis()->SetRangeUser(xmin,xmax);
  h1mb->GetXaxis()->SetRangeUser(xmin,xmax);

  TF1 *fspli = new TF1("fspli","sqrt([0]*[0]*pow(x,2*[1]))",xmin,xmax);
  fspli->SetParameters(fsda->GetParameter(3),fsda->GetParameter(4));
  TF1 *fspu = new TF1("fspu","sqrt([0]*[0]/(x*x) - [1]*[1]/(x*x))",xmin,xmax);
  fspu->SetParameters(fsdb->GetParameter(2),fsda->GetParameter(2));
  TF1 *fsjer = new TF1("fsjer","sqrt([0]*[0]+[1]*[1]/x+[2]*[2]/(x*x))",
		       xmin,xmax);
  fsjer->SetParameters(fsdb->GetParameter(0),fsdb->GetParameter(1),
		       fsda->GetParameter(2));
  TF1 *fjer = new TF1("fjer",fJER,xmin,xmax,2);
  fjer->SetParameters(0.65, 0);//, 25*1.0); // eta,rho

  fjer->SetLineStyle(kDotted);
  fjer->SetLineColor(kBlack);
  fjer->SetLineWidth(2);
  fjer->Draw("SAME");

  fspli->SetLineStyle(kDashed);
  fspli->SetLineColor(kCyan+2);
  fspli->Draw("SAME");
  fsjer->SetLineStyle(kDotted);
  fsjer->SetLineColor(kGreen+2);
  fsjer->SetLineWidth(2);
  fsjer->Draw("SAME");
  fspu->SetLineStyle(kDashDotted);
  fspu->SetLineColor(kOrange+2);
  fspu->Draw("SAME");

  if (stp=="tp") 
    tex->DrawLatex(0.50,0.50,Form("PLI: %1.2f #times p_{T,tag}^{%1.3f}",
				  fspli->GetParameter(0),
				  fspli->GetParameter(1)));
  else
    tex->DrawLatex(0.50,0.50,Form("PLI: %1.2f #times p_{T,ave}^{%1.3f}",
				  fspli->GetParameter(0),
				  fspli->GetParameter(1)));

  fsda->SetLineColor(kBlue);
  fsda->Draw("SAME");
  fsma->SetLineColor(kBlue);
  fsma->SetLineStyle(kDashed);
  fsma->Draw("SAME");
  //
  fsdb->Draw("SAME");
  fsmb->SetLineStyle(kDashed);
  fsmb->Draw("SAME");

  tdrDraw(h1db,"P",kFullSquare,kRed);
  tdrDraw(h1mb,"P",kOpenSquare,kRed);
  tdrDraw(h1da,"P",kFullCircle,kBlue);
  tdrDraw(h1ma,"P",kOpenCircle,kBlue);

  tex->DrawLatex(0.40,0.85,ce);
  tex->DrawLatex(0.40,0.75,Form("#alpha < %1.2f",alphamax));
  if (stp=="tp") tex->DrawLatex(0.60,0.75,"(tp)");

  TLegend *leg3up = tdrLeg(0.7,0.66,0.9,0.88);
  leg3up->AddEntry(h1db,"#LTB#GT_{data}","PL");
  leg3up->AddEntry(h1mb,"#LTB#GT_{MC}","PL");
  leg3up->AddEntry(h1da,"#LTA#GT_{data}","PL");
  leg3up->AddEntry(h1ma,"#LTA#GT_{MC}","PL");


  c3->cd(2);
  gPad->SetLogx();

  TH1D *hsra = (TH1D*)h1da->Clone("hsra");
  hsra->Divide(h1ma);
  TH1D *hsrb = (TH1D*)h1db->Clone("hsrb");
  hsrb->Divide(h1mb);

  TF1 *fsrb = new TF1("fsrb","sqrt([0]*[0]+[1]*[1]/x+[2]*[2]/(x*x))/"
		      "sqrt([3]*[3]+[4]*[4]/x+[5]*[5]/(x*x))",xmin,xmax);
  fsrb->SetParameters(fsdb->GetParameter(0), fsdb->GetParameter(1), 
		      fsdb->GetParameter(2), 
		      fsmb->GetParameter(0), fsmb->GetParameter(1), 
		      fsmb->GetParameter(2));
  TF1 *fsra = new TF1("fsra","sqrt([0]*[0]+[1]*[1]/x+[2]*[2]/(x*x))/"
		      "sqrt([3]*[3]+[4]*[4]/x+[5]*[5]/(x*x))",xmin,xmax);
  fsra->SetParameters(fsda->GetParameter(0), fsda->GetParameter(1), 
		      fsda->GetParameter(2), 
		      fsma->GetParameter(0), fsma->GetParameter(1), 
		      fsma->GetParameter(2));

  fsrb->Draw("SAME");
  fsra->SetLineColor(kBlue);
  fsra->Draw("SAME");
  tdrDraw(hsrb,"P",kFullSquare,kRed);
  tdrDraw(hsra,"P",kFullCircle,kBlue);

  TLegend *leg3dw = tdrLeg(0.6,0.8,0.9,0.9);
  leg3dw->SetTextSize(2.*0.045);
  leg3dw->SetNColumns(2);
  leg3dw->AddEntry(hsrb,"#sigma(B)","PL");
  leg3dw->AddEntry(hsra,"#sigma(A)","PL");

  c3->SaveAs(Form("pdf/deriveL2ResNarrow_RMS_%s_amax%1.0f%s.pdf",
		  ce2,alphamax*100,ctp));


  TH1D *hup4 = new TH1D("hup4",";p_{T,ave} (GeV);N_{entries}",// / pb^{-1}",
			1940,60,2000);
  const double ntot = 2.*36e6;
  hup4->SetMaximum(1e4 * ntot);
  hup4->SetMinimum(1e-8 * ntot); // was 1e-6
  TH1D *hdw4 = new TH1D("hdw4",";p_{T,ave} (GeV);Data/MC",
			1940,60,2000);
  hdw4->SetMaximum(1+0.8-1e-4);
  hdw4->SetMinimum(1-0.5+1e-4);
  hdw4->GetXaxis()->SetMoreLogLabels();
  hdw4->GetXaxis()->SetNoExponent();
  if (stp=="tp") hdw4->SetXTitle("p_{T,tag}");

  TCanvas *c4 = tdrDiCanvas("c4",hup4,hdw4,4,11);

  TF1 *fpt = new TF1("fpt","[0]*pow(x,[1])*pow(1-2.*x/13000.*cosh([3]),[2])",
		     xmin,xmax);
  fpt->SetParameters(1e16,-5,10,0);
  fpt->FixParameter(3,0.);

  c4->cd(1);
  gPad->SetLogx();
  gPad->SetLogy();

  hnma->GetXaxis()->SetRangeUser(xmin,xmax);
  hnda->GetXaxis()->SetRangeUser(xmin,xmax);

  //hnda->Scale(1./(2*36e6));
  hnma->Scale(ntot/20e6); // pThat bins require /20e6

  TH1D *hnda2 = (TH1D*)hnda->Clone("hnda2");
  for (int i = 1; i != hnda->GetNbinsX()+1; ++i) {
    int j = _htrg->FindBin(hnda->GetBinCenter(i));
    double pre = _htrg->GetBinContent(j);
    hnda2->SetBinContent(i, hnda->GetBinContent(i)*pre);
    hnda2->SetBinError(i, hnda->GetBinError(i)*pre);
  }
  hnda2->Fit(fpt,"QRN");

  fpt->Draw("SAME");
  tdrDraw(hnma,"P",kOpenCircle);
  tdrDraw(hnda2,"P",kFullCircle,kRed);
  tdrDraw(hnda,"P",kFullCircle);

  tex->DrawLatex(0.40,0.85,ce);
  tex->DrawLatex(0.40,0.75,Form("#alpha < %1.2f",alphamax));
  if (stp=="tp") tex->DrawLatex(0.60,0.75,"(tp)");

  TLegend *leg4up = tdrLeg(0.7,0.66,0.9,0.88);
  leg4up->AddEntry(hnda,"Data","P");
  leg4up->AddEntry(hnda2,"Data#timespre.","P");
  leg4up->AddEntry(hnma,"MC","P");

  c4->cd(2);
  gPad->SetLogx();

  TH1D *hnra = (TH1D*)hnda->Clone("hnra");
  //hnra->Divide(hnma);
  hnra->Divide(hnda,hnma,1,1,"e");
  TH1D *hnra2 = (TH1D*)hnda2->Clone("hnra2");
  //hnra2->Divide(hnma);
  hnra2->Divide(hnda2,hnma,1,1,"e");

  tdrDraw(hnra2,"P",kFullCircle,kRed);
  tdrDraw(hnra,"P",kFullCircle);

  c4->SaveAs(Form("pdf/deriveL2ResNarrow_Nev_%s_amax%1.0f%s.pdf",
		  ce2,alphamax*100,ctp));



  TH1D *hup5 = new TH1D("hup5",";p_{T,ave} (GeV);Cross section (pb/GeV)",
			1940,60,2000);
  const double lumi = 7.8e6; // RunG
  hup5->SetMaximum(1e4);
  hup5->SetMinimum(1e-6*0.1);
  TH1D *hdw5 = new TH1D("hdw5",";p_{T,ave} (GeV);Data/MC",
			1940,60,2000);
  hdw5->SetMaximum(1+0.8-1e-4);
  hdw5->SetMinimum(1-0.5+1e-4);
  hdw5->GetXaxis()->SetMoreLogLabels();
  hdw5->GetXaxis()->SetNoExponent();
  if (stp=="tp") hdw5->SetXTitle("p_{T,tag}");

  TCanvas *c5 = tdrDiCanvas("c5",hup5,hdw5,4,11);

  _eta = y1;//2.0;
  //TF1 *fxs = new TF1("fxs","[0]*pow(x,[1])*pow(1-2.*x/13000.*cosh([3]),[2])",
  //	     xmin,xmax);
  //fxs->SetParameters(1e11,-5,10,_eta);
  //fxs->FixParameter(3,_eta);
  TF1 *fxs = new TF1("fxs",fXsec,xmin,xmax,3);
  fxs->SetParameters(1e11,-5,10);

  c5->cd(1);
  gPad->SetLogx();
  gPad->SetLogy();

  TH1D *hxda = (TH1D*)hnda2->Clone("hxda");
  TH1D *hxma = (TH1D*)hnma->Clone("hxma");
  for (int i = 1; i != hnda2->GetNbinsX()+1; ++i) {
    double dpt = hnda2->GetBinWidth(i);
    hxda->SetBinContent(i, hnda2->GetBinContent(i)/dpt/lumi);
    hxda->SetBinError(i, hnda2->GetBinError(i)/dpt/lumi);
    hxma->SetBinContent(i, hnma->GetBinContent(i)/dpt/lumi);
    hxma->SetBinError(i, hnma->GetBinError(i)/dpt/lumi);
  }
  hxda->Fit(fxs,"RN");

  fxs->Draw("SAME");
  tdrDraw(hxma,"P",kOpenCircle);
  tdrDraw(hxda,"P",kFullCircle);

  tex->DrawLatex(0.40,0.85,ce);
  tex->DrawLatex(0.40,0.75,Form("#alpha < %1.2f",alphamax));
  if (stp=="tp") tex->DrawLatex(0.60,0.75,"(tp)");

  TLegend *leg5up = tdrLeg(0.7,0.73,0.9,0.88);
  leg5up->AddEntry(hxda,"Data","P");
  leg5up->AddEntry(hxma,"MC","P");

  c5->cd(2);
  gPad->SetLogx();

  TH1D *hxra = (TH1D*)hxda->Clone("hxra");
  hxra->Divide(hxda,hxma,1,1,"e");
  TH1D *hxra2 = (TH1D*)hxda->Clone("hxra2");
  //hxra2->Divide(fxs);
  // Ensure that the data is divived by the integral over the bin
  // Checked, and is not the same as direct division
  for (int i = 1; i != hxra2->GetNbinsX()+1; ++i) {
    double pt1 = hxra2->GetBinLowEdge(i);
    double pt2 = hxra2->GetBinLowEdge(i+1);
    double val = fxs->Integral(pt1, pt2) / (pt2 - pt1);
    hxra2->SetBinContent(i, hxra2->GetBinContent(i) / val);
    hxra2->SetBinError(i, hxra2->GetBinError(i) / val);
  }

  tdrDraw(hxra2,"P",kOpenCircle);
  tdrDraw(hxra,"P",kFullCircle);

  TLegend *leg5dw = tdrLeg(0.6,0.8,0.9,0.9);
  leg5dw->SetTextSize(2.*0.045);
  leg5dw->SetNColumns(2);
  leg5dw->AddEntry(hxra,"Data/MC","PL");
  leg5dw->AddEntry(hxra2,"Data/Fit","PL");

  c5->SaveAs(Form("pdf/deriveL2ResNarrow_Xsec_%s_amax%1.0f%s.pdf",
		  ce2,alphamax*100,ctp));

  return l2fit(f2b->GetChisquare(), f2b->GetNDF(), f2b, emat);

} // deriveL2ResNarrow_x
