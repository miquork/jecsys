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

// Ugly global variables...
//const char *ct = "G";
//const char *ct = "BCDEFGH";
const char *ct = "H";

struct l2fit {

  double chi2a, chi2b;
  int ndfa, ndfb;
  TF1 *f1a, *f1b;
  TMatrixD *emata, *ematb;
  
  l2fit(double chi2a_ = 0, int ndfa_ = 0, TF1 *f1a_ = 0, TMatrixD *emata_ = 0, 
	double chi2b_ = 0, int ndfb_ = 0, TF1 *f1b_ = 0, TMatrixD *ematb_ = 0) :
    chi2a(chi2a_), ndfa(ndfa_), f1a(f1a_), emata(emata_),
    chi2b(chi2b_), ndfb(ndfb_), f1b(f1b_), ematb(ematb_) { };
};

const double _xmin(30.), _xmax(3500.);

FactorizedJetCorrector *_jer_dt(0), *_jer_mc(0);
Double_t fJER(Double_t *x, Double_t *p) {
  // Note: still need to implement quadratic averaging over
  //       eta and rho bins for |eta|<1.3 and all rho
  //       Where is the sqrt(2) compared to sigma(asymmetry)?

  double pt = x[0];
  double eta = p[0];
  double eta0 = 0.65;
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

  _jer->setJetPt(pt);
  _jer->setJetEta(eta0);
  _jer->setRho(rho);
  double jer0 = _jer->getCorrection();
  
  //return jer / sqrt(2.);
  return sqrt(jer*jer +jer0*jer0) / 2.;
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

// Clean out bins with only 1 entry from data-MC pair
void cleanPair(TH1D *h1, TH1D* h2) {
  assert(h1);
  assert(h2);
  assert(h1->GetNbinsX()==h2->GetNbinsX());
  for (int i = 1; i != h1->GetNbinsX()+1; ++i) {
    if (h1->GetBinError(i)==0 || h2->GetBinError(i)==0) {
      h1->SetBinContent(i, 0);
      h1->SetBinError(i, 0);
      h2->SetBinContent(i, 0);
      h2->SetBinError(i, 0);
    }
  } // for i
} // cleanPair

l2fit deriveL2ResNarrow_x(double y1=2.0, double y2=2.5,
		    double alphamax=0.2, double xmax=2000.);

TH1D *_htrg(0);
// Main method (pTave, tag pT, probe pT)
string stp = "";//"tp";
const char *ctp = stp.c_str();//"tp";
// Cross-check method 1 (tag pT)
string stp2 = "tp";
const char *ctp2 = stp2.c_str();
// Cross-check method 2 (probe pT)
string stp3 = "pt";
const char *ctp3 = stp3.c_str();

void deriveL2ResNarrow() {

  // For spectrum reconstruction
  const int ntrg = 9;   // 1st 74 (jet40) 
  //double trigpt[ntrg+1] = {0, 84, 114, 196, 272, 330, 395, 468, 548,3500};
  double trigpt[ntrg+1] = {0, 84, 114, 174, 245, 330, 395, 468, 548,3500};
  // Effective prescales for RunFLateG_rereco from Engin root tuple
  //152728, 152655, 152261, 151443, 149893, 146976, 141570, 129504, 20007, 19543.2, 18588.5, 16868.2, 1986.39, 1902.99, 1758.02, 564.405, 511.402, 66.7384, 61.396, 22.0983, 20.2887, 7.31431, 6.60524,
  //double trigpre[ntrg] = {152728, 129504, 20007, 1986, 564.4, 66.74, 22.10, 7.314, 1};
  double trigpre[ntrg] = {152728, 50000, 20007, 1986, 564.4, 66.74, 22.10, 7.314, 1};
  //double trigpre[ntrg] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
  _htrg = new TH1D("htrg",";p_{T};prescale;",ntrg,trigpt);
  for (int i = 1; i != _htrg->GetNbinsX()+1; ++i) {
    _htrg->SetBinContent(i, trigpre[i-1]);
  }

  //double etas[] = {0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0,
  //	   3.2, 4.7};
  //double xmax[] = {2000, 2000, 2000, 1500, 1200, 800, 500, 350};
  // double etas[] = {0,0.261,0.522,0.783,1.044,1.305,
  // 		   1.479,1.653,1.93,2.172,2.322,2.5,
  // 		   2.65,2.853,2.964,3.139,3.489,3.839,
  // 		   5.191};
  // double xmax[] = {2000, 2000, 2000, 2000, 2000,
  // 		   1780, 1600, 1400, 1100, 1000, 900,
  // 		   800, 680, 630, 420, 290, 190,
  // 		   133};
  double etas[] = {0,0.261,0.522,0.783, 1.044,1.305,1.479,
		   1.653,1.93,2.172, 2.322,2.5,2.65,
		   2.853,2.964,3.139, 3.489,3.839,5.191};
  double xmax[] = {2000,2000,2000, 1784,1684,1497, // min 1%, 5%/sqrt(N)->25
		   1410,1248,1032, 905,790,686, // same min 1%, N>25
		   638,507,468, 395,300,220}; // eyeballing limits, low stat

  // SMP-J binning:
  //28, 32, 37, 43, 49, 56, 64, 74, 84,
  //97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
  // 507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
  // 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000,

  // Debug one bin first
//   double etas[] = {2.5,3.0};//2.65};
//   double xmax[] = {800};
  const int neta = sizeof(etas)/sizeof(etas[0]) - 1;
  assert(sizeof(xmax)/sizeof(xmax[0]) == neta);
  //double alphamax = 0.3;
  //double alphamax = 0.2;
  //double alphamax = 0.1;
  //double alphamax = 0.05;
  //double alphamax = 0.025;
  //double alphas[] = {0.025, 0.05, 0.1, 0.2, 0.3};
  //double alphas[] = {0.2};
  //double alphas[] = {0.05};
  //double alphas[] = {0.10};
  //double alphas[] = {0.20};
  //double alphas[] = {0.30};
  double alphas[] = {0.1, 0.15, 0.20, 0.30};
  const int nalpha = sizeof(alphas)/sizeof(alphas[0]);

  double pts[] = {60, 120, 240, 480};
  int colorpt[] = {kBlack, kRed, kGreen+2, kBlue};
  int markerpt[] = {0, kFullSquare, 0, 0};
  const int npt = sizeof(pts)/sizeof(pts[0]);
  vector<TH1D*> hptas(npt), hptbs(npt);
  for (int ipt = 0; ipt != npt; ++ipt) {
    hptas[ipt] = new TH1D(Form("hpt_%1.0f",pts[ipt]),
			  ";|#eta|;Relative correction (DB)",
			  neta, &etas[0]);
    hptbs[ipt] = new TH1D(Form("hpt_%1.0f",pts[ipt]),
			  ";|#eta|;Relative correction (MPF)",
			  neta, &etas[0]);
  }

  for (int ialpha = 0; ialpha != nalpha; ++ialpha) {

    const double alphamax = alphas[ialpha];

    vector<l2fit> l2s(neta);
    double maxchi2 = 5;
    TGraphErrors *gchi2a = new TGraphErrors(neta);
    TGraphErrors *gchi2b = new TGraphErrors(neta);
    for (int ieta = 0; ieta != neta; ++ieta) {
      
      double eta = 0.5*(etas[ieta] + etas[ieta+1]);
      double deta = 0.5*(etas[ieta+1] - etas[ieta]);
      l2s[ieta] = deriveL2ResNarrow_x(etas[ieta], etas[ieta+1],
				      alphamax, xmax[ieta]);
      
      gchi2a->SetPoint(ieta, eta, min(l2s[ieta].chi2a/l2s[ieta].ndfa, maxchi2));
      gchi2a->SetPointError(ieta, deta, 1. / sqrt(l2s[ieta].ndfa));

      gchi2b->SetPoint(ieta, eta, min(l2s[ieta].chi2b/l2s[ieta].ndfb, maxchi2));
      gchi2b->SetPointError(ieta, deta, 1. / sqrt(l2s[ieta].ndfb));
      
      for (int ipt = 0; ipt != npt; ++ipt) {
	double pt = pts[ipt];

	double A = l2s[ieta].f1a->Eval(pt);
	double Rrel_A = (1 + A) / (1 - A);
	double B = l2s[ieta].f1b->Eval(pt);
	double Rrel_B = (1 + B) / (1 - B);
	if (pt*cosh(eta) < 6500.) {

	  hptas[ipt]->SetBinContent(ieta+1, 1./Rrel_A);
	  if (markerpt[ipt]!=0) {
	    double err = getErr(pt, l2s[ieta].f1a, l2s[ieta].emata);
	    hptas[ipt]->SetBinError(ieta+1, err);
	  }

	  hptbs[ipt]->SetBinContent(ieta+1, 1./Rrel_B);
	  if (markerpt[ipt]!=0) {
	    double err = getErr(pt, l2s[ieta].f1b, l2s[ieta].ematb);
	    hptbs[ipt]->SetBinError(ieta+1, err);
	  }

	}
      }
    } // for ieta
    
    TH1D *h1 = new TH1D("h1",";#eta;#chi^{2} / NDF",neta,&etas[0]);
    h1->SetMinimum(0.);
    //h1->SetMaximum(3);
    h1->SetMaximum(maxchi2+1e-4);//5);
    
    TCanvas *c1 = tdrCanvas("c1",h1,4,11,kSquare);
    
    TLine *l1 = new TLine();
    l1->SetLineStyle(kDashed);
    l1->DrawLine(etas[0],1,etas[neta],1);
    
    tdrDraw(gchi2a,"Pz",kFullCircle);//kFullCircle);
    tdrDraw(gchi2b,"Pz",kOpenSquare);//kOpenCircle);

    TLegend *leg = tdrLeg(0.60,0.78,0.80,0.90);
    //leg->AddEntry(gchi2,"Log-lin","PL");
    //leg->AddEntry(gchi2b,"Quad-log","PL");
    leg->AddEntry(gchi2a,"DB quad-log","PL");
    leg->AddEntry(gchi2b,"MPF quad-log","PL");
    
    c1->SaveAs(Form("pdf/l2res/%s/deriveL2ResNarrow_chi2_amax%1.0f%s.pdf",
		    ct,alphamax*100,ctp));
    
    
    
    TH1D *h2a = new TH1D("h2",";#eta;Relative correction (DB)",neta,&etas[0]);
    h2a->SetMinimum(0.8);//0.98 - 0.03 - 0.05);
    h2a->SetMaximum(1.3);//1.19-1e-4 - 0.03 + 0.05);
    
    TCanvas *c2a = tdrCanvas("c2a",h2a,4,11,kSquare);
    
    TLine *l2a = new TLine();
    l2a->SetLineStyle(kDashed);
    l2a->DrawLine(etas[0],1,etas[neta],1);
    
    TLegend *leg2a = tdrLeg(0.17,0.48,0.47,0.78);
    leg2a->SetHeader("p_{T}(jet)");
    
    for (int ipt = 0; ipt != npt; ++ipt) {
      tdrDraw(hptas[ipt], markerpt[ipt]!=0 ? "HP][" : "H][", markerpt[ipt],
	      colorpt[ipt],kSolid,colorpt[ipt],kNone,0);
      hptas[ipt]->SetLineWidth(2);
      leg2a->AddEntry(hptas[ipt],Form("%1.0f GeV",pts[ipt]),
		      markerpt[ipt]!=0 ? "LP" : "L");
    }
    
    gPad->RedrawAxis();
    c2a->SaveAs(Form("pdf/l2res/%s/deriveL2ResNarrow_Fig15a_amax%1.0f%s.pdf",
		     ct,alphamax*100,ctp));

    TH1D *h2b = new TH1D("h2",";#eta;Relative correction (MPF)",neta,&etas[0]);
    h2b->SetMinimum(0.8);//0.98 - 0.03 - 0.05);
    h2b->SetMaximum(1.3);//1.19-1e-4 - 0.03 + 0.05);
    
    TCanvas *c2b = tdrCanvas("c2b",h2b,4,11,kSquare);
    
    TLine *l2b = new TLine();
    l2b->SetLineStyle(kDashed);
    l2b->DrawLine(etas[0],1,etas[neta],1);
    
    TLegend *leg2b = tdrLeg(0.17,0.48,0.47,0.78);
    leg2b->SetHeader("p_{T}(jet)");
    
    for (int ipt = 0; ipt != npt; ++ipt) {
      tdrDraw(hptbs[ipt], markerpt[ipt]!=0 ? "HP][" : "H][", markerpt[ipt],
	      colorpt[ipt],kSolid,colorpt[ipt],kNone,0);
      hptbs[ipt]->SetLineWidth(2);
      leg2b->AddEntry(hptbs[ipt],Form("%1.0f GeV",pts[ipt]),
		      markerpt[ipt]!=0 ? "LP" : "L");
    }
    
    gPad->RedrawAxis();
    c2b->SaveAs(Form("pdf/l2res/%s/deriveL2ResNarrow_Fig15b_amax%1.0f%s.pdf",
		     ct,alphamax*100,ctp));



    if (ialpha!=nalpha-1) { // leave last ones open
      delete h1;
      delete h2a;
      delete h2b;
    }
  } // for ialpha
}


l2fit deriveL2ResNarrow_x(double y1, double y2, double alphamax, double xmax) {

  TDirectory *curdir = gDirectory;
  //TFile *fd = new TFile(Form("rootfiles/exclusion_cmsweek/mikko/%s/output-DATA-2b.root",ct),"READ");
  TFile *fd = new TFile(Form("rootfiles/output_Run%se/mikko/output-DATA-2b.root",ct),"READ");
  //TFile *fd = new TFile("rootfiles/exclusion_cmsweek/mikko/GH/output-DATA-2b.root","READ");
  //TFile *fd = new TFile("rootfiles/exclusion_cmsweek/mikko/G/output-DATA-2b.root","READ");
  //TFile *fd = new TFile("rootfiles/exclusiontests/mikkofG/output-DATA-2b.root","READ");
  //TFile *fd = new TFile("rootfiles/exclusiontests/normalfG/output-DATA-2b.root","READ");
  //TFile *fd = new TFile("rootfiles/exclusiontests/normalH/output-DATA-2b.root","READ");
  //TFile *fd = new TFile("rootfiles/output_RunFlateG_Feb03V0/nol2l3res/output-DATA-2b.root","READ");
  //TFile *fd = new TFile("rootfiles/output_RunGV6x/output-DATA-2b.root","READ");
  //TFile *fd = new TFile("rootfiles/output_RunGV6/output-DATA-2b.root","READ");
  //TFile *fd = new TFile("rootfiles/output_RunGV6/output-DATA-2b-noL2Res.root",
  //		"READ");
  //TFile *fd = new TFile("rootfiles/output_RunGV6/output-DATA-2b.root","READ");
  //TFile *fd = new TFile("rootfiles/output_RunGV5/output-DATA-2b.root","READ");
  //TFile *fd = new TFile("rootfiles/output_RunGV3/output-DATA-2b.root","READ");
  //TFile *fd = new TFile("rootfiles/output_RunG_old/output-DATA-2b.root","READ");
  assert(fd && !fd->IsZombie());
  TFile *fm = new TFile(Form("rootfiles/exclusion_cmsweek/mc/%s/output-MC-2b.root",ct),"READ");
  //TFile *fm = new TFile("rootfiles/exclusion_cmsweek/mc/GH/output-MC-2b.root","READ");
  //TFile *fm = new TFile("rootfiles/exclusion_cmsweek/mc/G/output-MC-2b.root","READ");
  //TFile *fm = new TFile("rootfiles/exclusiontests/normalMC_fG/output-MC-2b.root","READ");
  //TFile *fm = new TFile("rootfiles/output_RunFlateG_Feb03V0/output-MC-2b.root","READ");
  //TFile *fm = new TFile("rootfiles/output_RunGV6x/output-MC-2b.root","READ");
  //TFile *fm = new TFile("rootfiles/output_RunGV6/output-MC-2b.root","READ");
  //TFile *fm = new TFile("rootfiles/output_RunGV5/output-MC-2b.root","READ");
  //TFile *fm = new TFile("rootfiles/output_RunGV3/output-MC-2b.root","READ");
  //TFile *fm = new TFile("rootfiles/output_RunG_old/output-MC-2b.root","READ");
  assert(fm && !fm->IsZombie());

  //TFile *fz = new TFile("rootfiles/L3res_mumu_Run2016FlateG.root","READ");
  //TFile *fz = new TFile("rootfiles/L3res_mumu_inclusive.root","READ");
  TFile *fz = new TFile(Form("rootfiles/L3res_mumu_Run2016%s.root",ct),"READ");
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
  //const char *cd = "FullEta";
  const char *cd = "FullEta_Reco";
  const char *ce0 = "";
  //string se = Form("Eta_%1.1f-%1.1f",y1,y2);
  string se = Form("Eta_%1.3f-%1.3f",y1,y2);
  const char *ce = se.c_str();
  string se2 = Form("Eta_%02.0f-%02.0f",10*y1,10*y2);
  const char *ce2 = se2.c_str();
  double xmin(30.);//60.);//, xmax(2000.);
  double fitxmin(64.);
  //const char *ca = "_a01";//Form("_a%");
  const char *ca = TString(Form("_a%1.3g",alphamax)).ReplaceAll(".","").Data();

  setTDRStyle();

  // Tag-and-probe with pT,tag binning: ctp = "tp"
  // Regular asymmetry with pT,ave binning: ctp = ""
  TH3D *h3da = (TH3D*)fd->Get(Form("%s/%s/hdjasymm%s%s",cd,ce0,ctp,ca));
  TH3D *h3db = (TH3D*)fd->Get(Form("%s/%s/hdjmpf%s%s",cd,ce0,ctp,ca));
  TH3D *h3ma = (TH3D*)fm->Get(Form("%s/%s/hdjasymm%s%s",cd,ce0,ctp,ca));
  TH3D *h3mb = (TH3D*)fm->Get(Form("%s/%s/hdjmpf%s%s",cd,ce0,ctp,ca));
  assert(h3da);
  assert(h3db);
  assert(h3ma);
  assert(h3mb);

  TH3D *h3dav2 = (TH3D*)fd->Get(Form("%s/%s/hdjasymm%s%s",cd,ce0,ctp2,ca));
  TH3D *h3dbv2 = (TH3D*)fd->Get(Form("%s/%s/hdjmpf%s%s",cd,ce0,ctp2,ca));
  TH3D *h3mav2 = (TH3D*)fm->Get(Form("%s/%s/hdjasymm%s%s",cd,ce0,ctp2,ca));
  TH3D *h3mbv2 = (TH3D*)fm->Get(Form("%s/%s/hdjmpf%s%s",cd,ce0,ctp2,ca));
  assert(h3dav2);
  assert(h3dbv2);
  assert(h3mav2);
  assert(h3mbv2);

  TH3D *h3dav3 = (TH3D*)fd->Get(Form("%s/%s/hdjasymm%s%s",cd,ce0,ctp3,ca));
  TH3D *h3dbv3 = (TH3D*)fd->Get(Form("%s/%s/hdjmpf%s%s",cd,ce0,ctp3,ca));
  TH3D *h3mav3 = (TH3D*)fm->Get(Form("%s/%s/hdjasymm%s%s",cd,ce0,ctp3,ca));
  TH3D *h3mbv3 = (TH3D*)fm->Get(Form("%s/%s/hdjmpf%s%s",cd,ce0,ctp3,ca));
  assert(h3dav3);
  assert(h3dbv3);
  assert(h3mav3);
  assert(h3mbv3);

  // https://github.com/miquork/jetphys/blob/chs/fillHistos.C#L1941-L1944
  // A = (P-T)/(T+P) => (1+A)*T = (1-A)*P => P/T = (1+A)/(1-A)
  // B = (P-T)/((P+T)/2) => (1+B/2)*T = (1-B/2)*P => P/T = (1-B/2)/(1-B/2)
  // => Scale B by 0.5 to get JEC 8 TeV paper definition, use A formula later
  //
  // B*pTave = METvec . unitVecTag = (-Pvec - Tvec) . unitTvec
  //         = -|P|*cos(pi) - |T|*cos(0) = P-T
  // for DeltaPhi=pi, but otherwise P*(1-delta) - T

  // First, look at mean and pT,ave spectrum
  int ymin = h3da->GetYaxis()->FindBin(y1);
  int ymax = h3da->GetYaxis()->FindBin(y2-1e-4);
  h3da->Sumw2();
  h3da->GetYaxis()->SetRange(ymin,ymax);
  TH2D *h2da = (TH2D*)h3da->Project3D("zxe"); // X=pT,ave; Y=Asymmetry
  h2da->Sumw2();
  h2da->SetName("h2da");
  TProfile *p1da = h2da->ProfileX("p1da",1,-1,"e");
  //if (stp=="tp") p1da->Scale(0.5);
  //p1da->Scale(0.5);
  TH1D *h1da = p1da->ProjectionX("h1da", "e");
  TH1D *hnda = h2da->ProjectionX("hnda",0,-1,"e");
  //
  h3ma->Sumw2();
  h3ma->GetYaxis()->SetRange(ymin,ymax);
  TH2D *h2ma = (TH2D*)h3ma->Project3D("zxe"); // X=pT,ave; Y=Asymmetry
  h2ma->Sumw2();
  h2ma->SetName("h2ma");
  TProfile *p1ma = h2ma->ProfileX("p1ma",1,-1,"e");
  //if (stp=="tp") p1ma->Scale(0.5);
  //p1ma->Scale(0.5);
  TH1D *h1ma = p1ma->ProjectionX("h1ma", "e");
  TH1D *hnma = h2ma->ProjectionX("hnma",0,-1,"e");

  h3db->Sumw2();
  h3db->GetYaxis()->SetRange(ymin,ymax);
  TH2D *h2db = (TH2D*)h3db->Project3D("zxe"); // X=pT,ave; Y=Asymmetry
  h2db->Sumw2();
  h2db->SetName("h2db");
  TProfile *p1db = h2db->ProfileX("p1db",1,-1,"e");
  //p1db->Scale(0.5); // Fixed for cmsweek
  TH1D *h1db = p1db->ProjectionX("h1db", "e");
  //
  h3mb->Sumw2();
  h3mb->GetYaxis()->SetRange(ymin,ymax);
  TH2D *h2mb = (TH2D*)h3mb->Project3D("zxe"); // X=pT,ave; Y=Asymmetry
  h2mb->Sumw2();
  h2mb->SetName("h2mb");
  TProfile *p1mb = h2mb->ProfileX("p1mb",1,-1,"e");
  //p1mb->Scale(0.5);  // Fixed for cmsweek
  TH1D *h1mb = p1mb->ProjectionX("h1dmb", "e");

  TH1D *h1ra = p1da->ProjectionX("h1ra","e");
  cleanPair(h1da, h1ma);
  h1ra->Add(h1da, h1ma, 1, -1);
  TH1D *h1rb = p1db->ProjectionX("h1rb","e");
  cleanPair(h1db, h1mb);
  h1rb->Add(h1db, h1mb, 1, -1);

  // Repeat for v2 DB and MPF
  h3dav2->GetYaxis()->SetRange(ymin,ymax);
  TH2D *h2dav2 = (TH2D*)h3dav2->Project3D("zxe"); // X=pT,ave; Y=Asymmetry
  h2dav2->Sumw2();
  h2dav2->SetName("h2dav2");
  TProfile *p1dav2 = h2dav2->ProfileX("p1dav2",1,-1,"e");
  //if (stp2=="tp") p1dav2->Scale(0.5);
  //p1dav2->Scale(0.5);
  TH1D *h1dav2 = p1dav2->ProjectionX("h1dav2", "e");
  //
  h3mav2->Sumw2();
  h3mav2->GetYaxis()->SetRange(ymin,ymax);
  TH2D *h2mav2 = (TH2D*)h3mav2->Project3D("zxe"); // X=pT,ave; Y=Asymmetry
  h2mav2->Sumw2();
  h2mav2->SetName("h2mav2");
  TProfile *p1mav2 = h2mav2->ProfileX("p1mav2",1,-1,"e");
  //if (stp2=="tp") p1mav2->Scale(0.5);
  //p1mav2->Scale(0.5);
  TH1D *h1mav2 = p1mav2->ProjectionX("h1mav2", "e");
  
  h3dbv2->Sumw2();
  h3dbv2->GetYaxis()->SetRange(ymin,ymax);
  TH2D *h2dbv2 = (TH2D*)h3dbv2->Project3D("zxe"); // X=pT,ave; Y=Asymmetry
  h2dbv2->Sumw2();
  h2dbv2->SetName("h2dbv2");
  TProfile *p1dbv2 = h2dbv2->ProfileX("p1dbv2",1,-1,"e");
  //p1dbv2->Scale(0.5); // Fixed for cmsweek
  TH1D *h1dbv2 = p1dbv2->ProjectionX("h1dbv2", "e");
  //
  h3mbv2->Sumw2();
  h3mbv2->GetYaxis()->SetRange(ymin,ymax);
  TH2D *h2mbv2 = (TH2D*)h3mbv2->Project3D("zxe"); // X=pT,ave; Y=Asymmetry
  h2mbv2->Sumw2();
  h2mbv2->SetName("h2mbv2");
  TProfile *p1mbv2 = h2mbv2->ProfileX("p1mbv2",1,-1,"e");
  //p1mbv2->Scale(0.5); // Fixed for cmsweek
  TH1D *h1mbv2 = p1mbv2->ProjectionX("h1mbv2", "e");

  TH1D *h1rav2 = p1dav2->ProjectionX("h1rav2","e");
  cleanPair(h1dav2, h1mav2);
  h1rav2->Add(h1dav2, h1mav2, 1, -1);
  TH1D *h1rbv2 = p1dbv2->ProjectionX("h1rbv2","e");
  cleanPair(h1dbv2, h1mbv2);
  h1rbv2->Add(h1dbv2, h1mbv2, 1, -1);


  // Repeat for v3 DB and MPF
  h3dav3->GetYaxis()->SetRange(ymin,ymax);
  TH2D *h2dav3 = (TH2D*)h3dav3->Project3D("zxe"); // X=pT,ave; Y=Asymmetry
  h2dav3->Sumw2();
  h2dav3->SetName("h2dav3");
  TProfile *p1dav3 = h2dav3->ProfileX("p1dav3",1,-1,"e");
  //if (stp3=="tp") p1dav3->Scale(0.5);
  //p1dav3->Scale(0.5);
  TH1D *h1dav3 = p1dav3->ProjectionX("h1dav3", "e");
  //
  h3mav3->Sumw2();
  h3mav3->GetYaxis()->SetRange(ymin,ymax);
  TH2D *h2mav3 = (TH2D*)h3mav3->Project3D("zxe"); // X=pT,ave; Y=Asymmetry
  h2mav3->Sumw2();
  h2mav3->SetName("h2mav3");
  TProfile *p1mav3 = h2mav3->ProfileX("p1mav3",1,-1,"e");
  //if (stp3=="tp") p1mav3->Scale(0.5);
  //p1mav3->Scale(0.5);
  TH1D *h1mav3 = p1mav3->ProjectionX("h1mav3", "e");
  
  h3dbv3->Sumw2();
  h3dbv3->GetYaxis()->SetRange(ymin,ymax);
  TH2D *h2dbv3 = (TH2D*)h3dbv3->Project3D("zxe"); // X=pT,ave; Y=Asymmetry
  h2dbv3->Sumw2();
  h2dbv3->SetName("h2dbv3");
  TProfile *p1dbv3 = h2dbv3->ProfileX("p1dbv3",1,-1,"e");
  //p1dbv3->Scale(0.5); // Fixed for cmsweek
  TH1D *h1dbv3 = p1dbv3->ProjectionX("h1dbv3", "e");
  //
  h3mbv3->Sumw2();
  h3mbv3->GetYaxis()->SetRange(ymin,ymax);
  TH2D *h2mbv3 = (TH2D*)h3mbv3->Project3D("zxe"); // X=pT,ave; Y=Asymmetry
  h2mbv3->Sumw2();
  h2mbv3->SetName("h2mbv3");
  TProfile *p1mbv3 = h2mbv3->ProfileX("p1mbv3",1,-1,"e");
  //p1mbv3->Scale(0.5); // Fixed for cmsweek
  TH1D *h1mbv3 = p1mbv3->ProjectionX("h1mbv3", "e");

  TH1D *h1rav3 = p1dav3->ProjectionX("h1rav3","e");
  cleanPair(h1dav3, h1mav3);
  h1rav3->Add(h1dav3, h1mav3, 1, -1);
  TH1D *h1rbv3 = p1dbv3->ProjectionX("h1rbv3","e");
  cleanPair(h1dbv3, h1mbv3);
  h1rbv3->Add(h1dbv3, h1mbv3, 1, -1);



  // Second, look at RMS / Gaussian sigma
  TF1 *fga = new TF1("fga","gaus",-1,1);
  TF1 *fgb = new TF1("fgb","gaus",-1,1);
  // fit range [mean-ksigma*sigma, mean+ksigma*sigma]
  const double ksigma = 1.5; 

  TH1D *h1das = p1da->ProjectionX("h1das","e"); h1das->Clear();
  TH1D *h1mas = p1ma->ProjectionX("h1mas","e"); h1mas->Clear();
  TH1D *h1dbs = p1db->ProjectionX("h1dbs","e"); h1dbs->Clear();
  TH1D *h1mbs = p1mb->ProjectionX("h1mbs","e"); h1mbs->Clear();
  for (int i = 1; i != h1das->GetNbinsX()+1; ++i) {
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
      h1das->SetBinContent(i, rmsa);
      h1das->SetBinError(i, htmpa->GetRMS()/sqrt(htmpa->GetEffectiveEntries()));
      h1dbs->SetBinContent(i, rmsb);
      h1dbs->SetBinError(i, htmpb->GetRMS()/sqrt(htmpb->GetEffectiveEntries()));
      //h1das->SetBinContent(i, (stp=="tp" ? 0.5 : 1)*rmsa);
      //h1das->SetBinError(i, (stp=="tp" ? 0.5 : 1)*htmpa->GetRMS()/sqrt(htmpa->GetEffectiveEntries()));
      //h1dbs->SetBinContent(i, 0.5*rmsb);
      //h1dbs->SetBinError(i, 0.5*htmpb->GetRMS()/sqrt(htmpb->GetEffectiveEntries()));
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
      h1mas->SetBinContent(i, rmsa);
      h1mas->SetBinError(i, htmpa->GetRMS()/sqrt(htmpa->GetEffectiveEntries())); 
      h1mbs->SetBinContent(i, rmsb);
      h1mbs->SetBinError(i, htmpb->GetRMS()/sqrt(htmpb->GetEffectiveEntries())); 
      //h1mas->SetBinContent(i, (stp=="tp" ? 0.5 : 1)*rmsa);
      //h1mas->SetBinError(i, (stp=="tp" ? 0.5 : 1)*htmpa->GetRMS()/sqrt(htmpa->GetEffectiveEntries())); 
      //h1mbs->SetBinContent(i, 0.5*rmsb);
      //h1mbs->SetBinError(i, 0.5*htmpb->GetRMS()/sqrt(htmpb->GetEffectiveEntries())); 
    }
    delete htmpa;
    delete htmpb;
  }

  curdir->cd();

  TH1D *hup = new TH1D("hup",";p_{T,ave} (GeV);#LTAsymmetry#GT",
		       1970,30,2000);
  hup->SetMaximum(+0.205);
  hup->SetMinimum(-0.155);
  TH1D *hdw = new TH1D("hdw",";p_{T,ave} (GeV);Data-MC",
		       1930,30,2000);
  hdw->SetMaximum(+0.105);
  hdw->SetMinimum(-0.105);
  hdw->GetXaxis()->SetMoreLogLabels();
  hdw->GetXaxis()->SetNoExponent();
  if (stp=="tp") hdw->SetXTitle("p_{T,tag}");


  //lumi_13TeV = "RunG";
  lumi_13TeV = Form("Run%s 03FebV3, 7.5 fb^{-1}",ct); // string() if autoclean
  //lumi_13TeV = "RunH 03FebV0, 8.8 fb^{-1}";
  TCanvas *c2 = tdrDiCanvas("c2",hup,hdw,4,11);

  c2->cd(1);
  gPad->SetLogx();

  // error bands from pT,tag, pT,probe binning
  TH1D *h1dae = p1dav2->ProjectionX("h1dae");
  TH1D *h1dbe = p1dbv2->ProjectionX("h1dbe");
  for (int i = 1; i != h1dae->GetNbinsX()+1; ++i) {
    if (p1dav2->GetBinContent(i)!=0 && p1dav2->GetBinError(i)!=0) {

      double aup = max(h1dav3->GetBinContent(i),
		       h1mav3->GetBinContent(i));
      double adw = min(h1dav2->GetBinContent(i),
		       h1mav2->GetBinContent(i));
      double bup = max(h1dbv3->GetBinContent(i),
		       h1mbv3->GetBinContent(i));
      double bdw = min(h1dbv2->GetBinContent(i),
		       h1mbv2->GetBinContent(i));

      h1dae->SetBinContent(i, 0.5*(aup+adw));
      h1dae->SetBinError(i, 0.5*(aup-adw));

      h1dbe->SetBinContent(i, 0.5*(bup+bdw));
      h1dbe->SetBinError(i, 0.5*(bup-bdw));
    }
  }
  h1dae->GetXaxis()->SetRangeUser(xmin,xmax);
  h1dbe->GetXaxis()->SetRangeUser(xmin,xmax);
  tdrDraw(h1dae,"E3",kNone,kBlue-9,kSolid,-1,1001,kBlue-9);
  h1dae->SetFillColorAlpha(kBlue-9, 0.35); // 35% transparent
  tdrDraw(h1dbe,"E3",kNone,kRed-9,kSolid,-1,1001,kRed-9);
  h1dbe->SetFillColorAlpha(kRed-9, 0.35); // 35% transparent

  // pT,tag
  h1dav2->GetXaxis()->SetRangeUser(xmin,xmax);
  h1mav2->GetXaxis()->SetRangeUser(xmin,xmax);
  h1dbv2->GetXaxis()->SetRangeUser(xmin,xmax);
  h1mbv2->GetXaxis()->SetRangeUser(xmin,xmax);
  tdrDraw(h1dav2,"HISTPL",kFullDiamond,kBlue-9,kSolid,-1,kNone);
  tdrDraw(h1mav2,"HISTPL",kOpenDiamond,kBlue-9,kSolid,-1,kNone);
  tdrDraw(h1dbv2,"HISTPL",kFullDiamond,kRed-9,kSolid,-1,kNone);
  tdrDraw(h1mbv2,"HISTPL",kOpenDiamond,kRed-9,kSolid,-1,kNone);
  h1dav2->SetMarkerSize(0.7);
  h1mav2->SetMarkerSize(0.7);
  h1dbv2->SetMarkerSize(0.7);
  h1mbv2->SetMarkerSize(0.7);
  // pT,probe
  h1dav3->GetXaxis()->SetRangeUser(xmin,xmax);
  h1mav3->GetXaxis()->SetRangeUser(xmin,xmax);
  h1dbv3->GetXaxis()->SetRangeUser(xmin,xmax);
  h1mbv3->GetXaxis()->SetRangeUser(xmin,xmax);
  tdrDraw(h1dav3,"HISTPL",kFullDiamond,kAzure-9,kSolid,-1,kNone);
  tdrDraw(h1mav3,"HISTPL",kOpenDiamond,kAzure-9,kSolid,-1,kNone);
  tdrDraw(h1dbv3,"HISTPL",kFullDiamond,kOrange-9,kSolid,-1,kNone);
  tdrDraw(h1mbv3,"HISTPL",kOpenDiamond,kOrange-9,kSolid,-1,kNone);
  h1dav3->SetMarkerSize(0.7);
  h1mav3->SetMarkerSize(0.7);
  h1dbv3->SetMarkerSize(0.7);
  h1mbv3->SetMarkerSize(0.7);
  // pT,ave
  h1da->GetXaxis()->SetRangeUser(xmin,xmax);
  h1ma->GetXaxis()->SetRangeUser(xmin,xmax);
  h1db->GetXaxis()->SetRangeUser(xmin,xmax);
  h1mb->GetXaxis()->SetRangeUser(xmin,xmax);
  tdrDraw(h1db,"P",kFullSquare,kRed);
  tdrDraw(h1mb,"P",kOpenSquare,kRed);
  tdrDraw(h1da,"P",kFullCircle,kBlue);
  tdrDraw(h1ma,"P",kOpenCircle,kBlue);

  // Add Z+jet
  vector<string> zs;
  zs.push_back("Data_MPF");
  zs.push_back("MC_MPF");
  zs.push_back("Ratio_MPF");
  zs.push_back("Data_PtBal");
  zs.push_back("MC_PtBal");
  zs.push_back("Ratio_PtBal");
  map<string, TGraphErrors*> gzs;

  for (int iz = 0; iz != zs.size(); ++iz) {

    string sz0 = Form("%s_CHS_a%d_eta_%04.0f_%04.0f_L1L2L3",
		      zs[iz].c_str(), int(alphamax*100+0.5), 1000.*0,1000.*1.3);
    cout << sz0 << endl << flush;
    TGraphErrors *gz0 = (TGraphErrors*)fz->Get(sz0.c_str());
    assert(gz0);

    string sz = Form("%s_CHS_a%d_eta_%04.0f_%04.0f_L1L2L3",
		     zs[iz].c_str(), int(alphamax*100+0.5), 1000.*y1,1000.*y2);
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

    gzs[zs[iz]] = gz;
  } // for iz

  // Z+jet
  tdrDraw(gzs["Data_MPF"],"Pz",kFullSquare,kGreen+3);
  gzs["Data_MPF"]->SetMarkerSize(0.7);
  tdrDraw(gzs["MC_MPF"],"Pz",kOpenSquare,kGreen+3);
  gzs["MC_MPF"]->SetMarkerSize(0.7);
  tdrDraw(gzs["Data_PtBal"],"Pz",kFullCircle,kGreen+2);
  gzs["Data_PtBal"]->SetMarkerSize(0.7);
  tdrDraw(gzs["MC_PtBal"],"Pz",kOpenCircle,kGreen+2);
  gzs["MC_PtBal"]->SetMarkerSize(0.7);
  
  TF1 *f1da = new TF1("f1da","[0]+[1]*log(x)+[2]*log(x)*log(x)",fitxmin,xmax);
  f1da->SetParameters(0,0.01,0.0001);
  h1da->Fit(f1da,"QRN");
  f1da->SetLineColor(kBlue);
  f1da->Draw("SAME");
  TF1 *f1ma = new TF1("f1ma","[0]+[1]*log(x)+[2]*log(x)*log(x)",fitxmin,xmax);
  f1ma->SetParameters(0,0.01,0.0001);
  h1ma->Fit(f1ma,"QRN");
  f1ma->SetLineColor(kBlue);
  f1ma->SetLineStyle(kDashed);
  f1ma->Draw("SAME");
  //
  TF1 *f1db = new TF1("f1db","[0]+[1]*log(x)+[2]*log(x)*log(x)",fitxmin,xmax);
  f1db->SetParameters(0,0.01,0.0001);
  h1db->Fit(f1db,"QRN");
  f1db->SetLineColor(kRed);
  f1db->Draw("SAME");
  TF1 *f1mb = new TF1("f1mb","[0]+[1]*log(x)+[2]*log(x)*log(x)",fitxmin,xmax);
  f1mb->SetParameters(0,0.01,0.0001);
  h1mb->Fit(f1mb,"QRN");
  f1mb->SetLineColor(kRed);
  f1mb->SetLineStyle(kDashed);
  f1mb->Draw("SAME");

  // Check estimated resolution bias
  //TF1 *fddpt = new TF1("fddpt",deltaPt,fitxmin,xmax,3);
  //fddpt->SetParameters(0.65,0.5*(y1+y2), 1);//0.0, 2.5, 1);
  //fddpt->SetLineStyle(kDotted);
  //fddpt->SetLineColor(kGreen+2);
  //fddpt->SetLineWidth(2);
  //fddpt->Draw("SAME");
  //TF1 *fmdpt = new TF1("fmdpt",deltaPt,fitxmin,xmax,3);
  //fmdpt->SetParameters(0.65, 0.5*(y1+y2), 0);//0.0, 2.5, 0);
  //fmdpt->SetLineStyle(kDotted);
  //fmdpt->SetLineColor(kMagenta+1);
  //fmdpt->SetLineWidth(2);
  //fmdpt->Draw("SAME");
  
  TLatex *tex = new TLatex();
  tex->SetTextSize(0.045);
  tex->SetNDC();
  tex->DrawLatex(0.40,0.85,ce);
  tex->DrawLatex(0.40,0.75,Form("#alpha < %1.2f",alphamax));
  if (stp=="tp")  tex->DrawLatex(0.60,0.75,"(tp)");

  TLegend *leg2up = tdrLeg(0.7,0.66,0.9,0.88);
  leg2up->AddEntry(h1db,"#LTB#GT_{data}","PL");
  leg2up->AddEntry(h1mb,"#LTB#GT_{MC}","PL");
  leg2up->AddEntry(h1da,"#LTA#GT_{data}","PL");
  leg2up->AddEntry(h1ma,"#LTA#GT_{MC}","PL");

  TLegend *leg2up2 = tdrLeg(0.7,0.06,0.9,0.28);
  leg2up2->AddEntry(h1dbv2,"#LTB#GT_{tag}","PL");
  leg2up2->AddEntry(h1dbv3,"#LTB#GT_{probe}","PL");
  leg2up2->AddEntry(h1dav2,"#LTA#GT_{tag}","PL");
  leg2up2->AddEntry(h1dav3,"#LTA#GT_{probe}","PL");


  c2->cd(2);
  gPad->SetLogx();

  // error bands from pT,tag, pT,probe binning
  TH1D *h1rae = (TH1D*) h1rav2->Clone("h1rae");
  TH1D *h1rbe = (TH1D*) h1rbv2->Clone("h1rbe");
  for (int i = 1; i != h1rae->GetNbinsX()+1; ++i) {
    if (h1rav2->GetBinContent(i)!=0 && h1rav2->GetBinError(i)!=0) {

      h1rae->SetBinContent(i, 0.5*(h1rav2->GetBinContent(i)
				   +h1rav3->GetBinContent(i)));
      h1rae->SetBinError(i, 0.5*fabs(h1rav2->GetBinContent(i)
				     -h1rav3->GetBinContent(i)));

      h1rbe->SetBinContent(i, 0.5*(h1rbv2->GetBinContent(i)
				   +h1rbv3->GetBinContent(i)));
      h1rbe->SetBinError(i, 0.5*fabs(h1rbv2->GetBinContent(i)
				     -h1rbv3->GetBinContent(i)));
    }
  }
  h1rae->GetXaxis()->SetRangeUser(xmin,xmax);
  h1rbe->GetXaxis()->SetRangeUser(xmin,xmax);
  tdrDraw(h1rae,"E3",kNone,kBlue-9,kSolid,-1,1001,kBlue-9);
  h1rae->SetFillColorAlpha(kBlue-9, 0.35); // 35% transparent
  tdrDraw(h1rbe,"E3",kNone,kRed-9,kSolid,-1,1001,kRed-9);
  h1rbe->SetFillColorAlpha(kRed-9, 0.35); // 35% transparent

  // pT,tag
  h1rav2->GetXaxis()->SetRangeUser(xmin,xmax);
  h1rbv2->GetXaxis()->SetRangeUser(xmin,xmax);
  tdrDraw(h1rav2,"HISTPL",kFullDiamond,kBlue-9,kSolid,-1,kNone);
  tdrDraw(h1rbv2,"HISTPL",kFullDiamond,kRed-9,kSolid,-1,kNone);
  h1rav2->SetMarkerSize(0.7);
  h1rbv2->SetMarkerSize(0.7);
  // pT,probe
  h1rav3->GetXaxis()->SetRangeUser(xmin,xmax);
  h1rbv3->GetXaxis()->SetRangeUser(xmin,xmax);
  tdrDraw(h1rav3,"HISTPL",kFullDiamond,kAzure-9,kSolid,-1,kNone);
  tdrDraw(h1rbv3,"HISTPL",kFullDiamond,kOrange-9,kSolid,-1,kNone);
  h1rav3->SetMarkerSize(0.7);
  h1rbv3->SetMarkerSize(0.7);
  // pT,ave binning on top
  h1ra->GetXaxis()->SetRangeUser(xmin,xmax);
  h1rb->GetXaxis()->SetRangeUser(xmin,xmax);
  tdrDraw(h1rb,"P",kFullSquare,kRed);
  tdrDraw(h1ra,"P",kFullCircle,kBlue);

  // Z+jet
  tdrDraw(gzs["Ratio_PtBal"],"Pz",kFullCircle,kGreen+2);
  gzs["Ratio_PtBal"]->SetMarkerSize(0.7);
  tdrDraw(gzs["Ratio_MPF"],"Pz",kOpenSquare,kGreen+3);
  gzs["Ratio_MPF"]->SetMarkerSize(0.7);
  
  TF1 *f1a = new TF1("f1a","[0]+[1]*log(x)",fitxmin,xmax);
  f1a->SetParameters(0,0.01);
  h1ra->Fit(f1a,"QRN");
  f1a->SetLineColor(kBlue);
  f1a->Draw("SAME");

  // Get error matrix for DB (until MPF fixed in SMP-J tuples)
  TMatrixD *emat1a = new TMatrixD(f1a->GetNpar(), f1a->GetNpar());
  gMinuit->mnemat(emat1a->GetMatrixArray(), f1a->GetNpar());

  TF1 *f1b = new TF1("f1b","[0]+[1]*log(x)",fitxmin,xmax);
  f1b->SetParameters(0,0.01);
  h1rb->Fit(f1b,"QRN");
  f1b->SetLineColor(kRed);
  f1b->Draw("SAME");

  // Get matrix for MPF
  TMatrixD *emat1b = new TMatrixD(f1b->GetNpar(), f1b->GetNpar());
  gMinuit->mnemat(emat1b->GetMatrixArray(), f1b->GetNpar());

  TF1 *f2a = new TF1("f2a","[0]+[1]*log(x)+[2]*log(x)*log(x)",fitxmin,xmax);
  f2a->SetParameters(0,0.01,0.0001);
  h1ra->Fit(f2a,"QRN");
  f2a->SetLineColor(kBlue);
  f2a->SetLineStyle(kDashed);
  f2a->Draw("SAME");

  // Get error matrix for DB (until MPF fixed in SMP-J tuples)
  TMatrixD *emat2a = new TMatrixD(f2a->GetNpar(), f2a->GetNpar());
  gMinuit->mnemat(emat2a->GetMatrixArray(), f2a->GetNpar());

  TF1 *f2b = new TF1("f2b","[0]+[1]*log(x)+[2]*log(x)*log(x)",fitxmin,xmax);
  f2b->SetParameters(0,0.01,0.0001);
  h1rb->Fit(f2b,"QRN");
  f2b->SetLineColor(kRed);
  f2b->SetLineStyle(kDashed);
  f2b->Draw("SAME");

  // Get matrix for MPF
  TMatrixD *emat2b = new TMatrixD(f2b->GetNpar(), f2b->GetNpar());
  gMinuit->mnemat(emat2b->GetMatrixArray(), f2b->GetNpar());

  //TMatrixD *emat = new TMatrixD(f2b->GetNpar(), f2b->GetNpar());
  //gMinuit->mnemat(emat->GetMatrixArray(), f2b->GetNpar());

  TLegend *leg2dw = tdrLeg(0.6,0.8,0.9,0.9);
  leg2dw->SetTextSize(2.*0.045);
  leg2dw->SetNColumns(2);
  leg2dw->AddEntry(h1rb,"#LTB#GT","PL");
  leg2dw->AddEntry(h1ra,"#LTA#GT","PL");
  
  TF1 *f1az = new TF1("f1az","[0]+[1]*log(x)",fitxmin,xmax);
  f1az->SetParameters(f1a->GetParameter(0),f1a->GetParameter(1));
  f1az->SetLineStyle(kDotted);
  f1az->SetLineColor(kBlue);
  f1az->Draw("SAME");
  TF1 *f1bz = new TF1("f1bz","[0]+[1]*log(x)",fitxmin,xmax);
  f1bz->SetParameters(f1b->GetParameter(0),f1b->GetParameter(1));
  f1bz->SetLineStyle(kDotted);
  f1bz->SetLineColor(kRed);
  f1bz->Draw("SAME");

  TF1 *f2az = new TF1("f2az","[0]+[1]*log(x)+[2]*log(x)*log(x)",fitxmin,xmax);
  f2az->SetParameters(f2a->GetParameter(0),f2a->GetParameter(1),
		      f2a->GetParameter(2));
  f2az->SetLineStyle(kDotted);
  f2az->SetLineColor(kBlue);
  f2az->SetLineStyle(kDashed);
  f2az->Draw("SAME");
  TF1 *f2bz = new TF1("f2bz","[0]+[1]*log(x)+[2]*log(x)*log(x)",fitxmin,xmax);
  f2bz->SetParameters(f2b->GetParameter(0),f2b->GetParameter(1),
		      f2b->GetParameter(2));
  f2bz->SetLineStyle(kDotted);
  f2bz->SetLineColor(kRed);
  f2bz->SetLineStyle(kDashed);
  f2bz->Draw("SAME");

  //tdrDraw(gz,"Pz",kOpenSquare,kGreen+2);

  TLine *l = new TLine();
  l->SetLineStyle(kDotted);
  l->DrawLine(xmin,0,xmax,0);

  gPad->RedrawAxis();

  c2->SaveAs(Form("pdf/l2res/%s/deriveL2ResNarrow_%s_amax%1.0f%s.pdf",
		  ct,ce2,alphamax*100,ctp));


  TH1D *hup3 = new TH1D("hup3",";p_{T,ave} (GeV);#sigma(asymmetry)",
			1970,30,2000);
  hup3->SetMaximum(0.3);//0.5);
  hup3->SetMinimum(0.+1e-4);
  TH1D *hdw3 = new TH1D("hdw3",";p_{T,ave} (GeV);Data/MC",
		       1970,30,2000);
  hdw3->SetMaximum(1+0.5-1e-4);
  hdw3->SetMinimum(1-0.3+1e-4);
  hdw3->GetXaxis()->SetMoreLogLabels();
  hdw3->GetXaxis()->SetNoExponent();
  if (stp=="tp") hdw3->SetXTitle("p_{T,tag}");

  TCanvas *c3 = tdrDiCanvas("c3",hup3,hdw3,4,11);

  TF1 *fsdb = new TF1("fsdb","sqrt([0]*[0]+[1]*[1]/x+[2]*[2]/(x*x))",xmin,xmax);
  fsdb->SetParameters(0.038, 0.9, 10.4);
  h1dbs->Fit(fsdb,"RN");
  TF1 *fsmb = new TF1("fsmb","sqrt([0]*[0]+[1]*[1]/x+[2]*[2]/(x*x))",xmin,xmax);
  fsmb->SetParameters(0.028, 0.9, 10.4);
  h1mbs->Fit(fsmb,"QRN");

  TF1 *fsda = new TF1("fsda","sqrt([0]*[0]+[1]*[1]/x+[2]*[2]/(x*x)+"
		      "[3]*[3]*pow(x,2*[4]))",xmin,xmax);
  fsda->SetParameters(0.038, 0.9, 4.1, 0.46*alphamax,-0.29);
  fsda->FixParameter(0, fsdb->GetParameter(0));
  fsda->FixParameter(1, fsdb->GetParameter(1));
  //fsda->FixParameter(2, fsdb->GetParameter(2));
  // https://twiki.cern.ch/twiki/pub/CMSPublic/MultipleConeSizes14/jerplots_ParFits_Combined.pdf (muxA = 25*pi*0.5^2 = 12.5 for 2016, N~2.5);
  double parN = sqrt(-1*fabs(-1) + 0.75*0.75*25*TMath::Pi()*0.4*0.4);
  fsda->FixParameter(2, parN);
  h1das->Fit(fsda,"RN");
  TF1 *fsma = new TF1("fsma","sqrt([0]*[0]+[1]*[1]/x+[2]*[2]/(x*x)+"
		      "[3]*[3]*pow(x,2*[4]))",xmin,xmax);
  //"[3]*[3]*exp(x*2*[4]))",xmin,xmax);
  fsma->SetParameters(0.038, 0.9, 4.1, 0.46*alphamax,-0.29);
  fsma->FixParameter(0, fsmb->GetParameter(0));
  fsma->FixParameter(1, fsmb->GetParameter(1));
  //fsma->FixParameter(2, fsmb->GetParameter(2));
  fsma->FixParameter(2, parN);
  h1mas->Fit(fsma,"QRN");

  // Note: PLI model (fit vs pTave and alphamax) is not yet ideal
  //       It's not fully describing the low pTave end of asymmetry
  //       However, could be that alphamax >= 24 GeV/pTave is kicking in
  //       for low values of pTave and alphamax, confusing the picture


  c3->cd(1);
  gPad->SetLogx();
  h1das->GetXaxis()->SetRangeUser(xmin,xmax);
  h1mas->GetXaxis()->SetRangeUser(xmin,xmax);
  h1dbs->GetXaxis()->SetRangeUser(xmin,xmax);
  h1mbs->GetXaxis()->SetRangeUser(xmin,xmax);

  TF1 *fspli = new TF1("fspli","sqrt([0]*[0]*pow(x,2*[1]))",xmin,xmax);
  fspli->SetParameters(fsda->GetParameter(3),fsda->GetParameter(4));
  TF1 *fspu = new TF1("fspu","sqrt([0]*[0]/(x*x) - [1]*[1]/(x*x))",xmin,xmax);
  fspu->SetParameters(fsdb->GetParameter(2),fsda->GetParameter(2));
  TF1 *fsjer = new TF1("fsjer","sqrt([0]*[0]+[1]*[1]/x+[2]*[2]/(x*x))",
		       xmin,xmax);
  fsjer->SetParameters(fsdb->GetParameter(0),fsdb->GetParameter(1),
		       fsda->GetParameter(2));
  TF1 *fjer = new TF1("fjer",fJER,xmin,xmax,2);
  fjer->SetParameters(0.5*(y1+y2), 15.41);//0.65, 0);//, 25*1.0); // eta,rho

  fjer->SetLineStyle(kDotted);
  fjer->SetLineColor(kBlack);
  fjer->SetLineWidth(3);//2);
  fjer->Draw("SAME");

  fspli->SetLineStyle(kDashed);
  fspli->SetLineColor(kCyan+2);
  //fspli->Draw("SAME");
  fsjer->SetLineStyle(kDotted);
  fsjer->SetLineColor(kGreen+2);
  fsjer->SetLineWidth(2);
  //fsjer->Draw("SAME");
  fspu->SetLineStyle(kDashDotted);
  fspu->SetLineColor(kOrange+2);
  //fspu->Draw("SAME");

  //if (stp=="tp") 
  //tex->DrawLatex(0.50,0.50,Form("PLI: %1.2f #times p_{T,tag}^{%1.3f}",
  //			  fspli->GetParameter(0),
  //			  fspli->GetParameter(1)));
  //else
  //tex->DrawLatex(0.50,0.50,Form("PLI: %1.2f #times p_{T,ave}^{%1.3f}",
  //			  fspli->GetParameter(0),
  //			  fspli->GetParameter(1)));

  fsda->SetLineColor(kBlue);
  fsda->Draw("SAME");
  fsma->SetLineColor(kBlue);
  fsma->SetLineStyle(kDashed);
  fsma->Draw("SAME");
  //
  fsdb->Draw("SAME");
  fsmb->SetLineStyle(kDashed);
  fsmb->Draw("SAME");

  tdrDraw(h1dbs,"P",kFullSquare,kRed);
  tdrDraw(h1mbs,"P",kOpenSquare,kRed);
  tdrDraw(h1das,"P",kFullCircle,kBlue);
  tdrDraw(h1mas,"P",kOpenCircle,kBlue);

  tex->DrawLatex(0.40,0.85,ce);
  tex->DrawLatex(0.40,0.75,Form("#alpha < %1.2f",alphamax));
  if (stp=="tp") tex->DrawLatex(0.60,0.75,"(tp)");

  TLegend *leg3up = tdrLeg(0.7,0.60,0.9,0.88);
  leg3up->AddEntry(h1dbs,"#LTB#GT_{data}","PL");
  leg3up->AddEntry(h1mbs,"#LTB#GT_{MC}","PL");
  leg3up->AddEntry(h1das,"#LTA#GT_{data}","PL");
  leg3up->AddEntry(h1mas,"#LTA#GT_{MC}","PL");
  //leg3up->AddEntry(fjer,"JER/#sqrt{2}","L");
  leg3up->AddEntry(fjer,"#sigma_{0}#oplus#sigma_{1}/2","L");


  c3->cd(2);
  gPad->SetLogx();

  TH1D *hsra = (TH1D*)h1das->Clone("hsra");
  hsra->Divide(h1mas);
  TH1D *hsrb = (TH1D*)h1dbs->Clone("hsrb");
  hsrb->Divide(h1mbs);

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

  c3->SaveAs(Form("pdf/l2res/%s/deriveL2ResNarrow_RMS_%s_amax%1.0f%s.pdf",
		  ct,ce2,alphamax*100,ctp));


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
  hnma->Scale(1e-9);

  TH1D *hnda2 = (TH1D*)hnda->Clone("hnda2");
  for (int i = 1; i != hnda->GetNbinsX()+1; ++i) {
    int j = _htrg->FindBin(hnda->GetBinCenter(i));
    double pre = _htrg->GetBinContent(j);
    hnda2->SetBinContent(i, hnda->GetBinContent(i)/pre);
    hnda2->SetBinError(i, hnda->GetBinError(i)/pre);
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

  c4->SaveAs(Form("pdf/l2res/%s/deriveL2ResNarrow_Nev_%s_amax%1.0f%s.pdf",
		  ct,ce2,alphamax*100,ctp));



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

  c5->SaveAs(Form("pdf/l2res/%s/deriveL2ResNarrow_Xsec_%s_amax%1.0f%s.pdf",
		  ct,ce2,alphamax*100,ctp));

  //return l2fit(f2b->GetChisquare(), f2b->GetNDF(), f2b, emat); // MPF
  // Use pT balance until MET fixed in SMP-J tuples
  //return l2fit(f2a->GetChisquare(), f2a->GetNDF(), f2a, emat); // pTbal
  return l2fit(f2a->GetChisquare(), f1a->GetNDF(), f2a, emat2a,
	       f2b->GetChisquare(), f2b->GetNDF(), f2b, emat2b);

} // deriveL2ResNarrow_x
