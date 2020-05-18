
// Purpose: Estimate effective bias in JEC (deltaJEC) from inclusive jet
//          cross section relative to published 2015 data (50 ns collisions).
//          Useful for estimating time stability of JEC throughout Run 2.
// run with 'root -l minitools/mk_drawDeltaJEC.C'
#include "TFile.h"
#include "TF1.h"
#include "TLatex.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TF2.h"
#include "TGraphErrors.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "ptresolution.h"
#include "../tdrstyle_mod15.C"

const double _prec = 1e-5; // precision of xsec and JEC inversion

const bool correctECALprefire17H = false;//true;
const bool correctECALprefire = false;//true;
const bool unfold17H = true;
const bool unfoldData = false;
const bool getDetData = false; // detector level data
const bool getFwd = true; // forward unfolded instead of Dagostini unfolded
//const bool correctFilterEff = false;

bool plotMC = true; // 0.0-1.3 missing still
bool plot2015 = false;
bool plotVs17UL = false; // set to true for 17UL

// "Rebin(2)"
const int nptb = 51+10+6;
const float ptbins[nptb+1] =
  {15, 18, 21, 24, 28, 32, 37, 43, 49, 56,
   64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362,
   395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032,
   1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116,
   2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832,
   4037, 4252, 4477, 4713, 4961, 5220};

//TH1F *rebinXsec(TH1F *h) {
TH1D *rebinXsec(TH1D *h) {

  //TH1F *h2 = new TH1F(Form("%s_rb2",h->GetName()),"",nptb,ptbins);
  TH1D *h2 = new TH1D(Form("%s_rb2",h->GetName()),"",nptb,ptbins);
  for (int i = 1; i != h2->GetNbinsX()+1; ++i) {
    float pt1 = h2->GetBinLowEdge(i);
    float pt2 = h2->GetBinLowEdge(i+1);
    int j1 = h->FindBin(pt1);
    int j2 = h->FindBin(pt2)-1;
    if (!(h2->GetBinLowEdge(i)==h->GetBinLowEdge(j1))) {
      cout << h2->GetBinLowEdge(i) << ", " << h->GetBinLowEdge(j1)
	   << endl << flush;
      assert(false);
    }
    if (!(h2->GetBinLowEdge(i+1)==h->GetBinLowEdge(j2+1))) {
      cout << h2->GetBinLowEdge(i+1) << ", " << h->GetBinLowEdge(j2+1)
	   << endl << flush;
      assert(false);
    }

    double sumx(0);
    double errx2(0);
    for (int j = j1; j != j2+1; ++j) {
      sumx += h->GetBinContent(j) * h->GetBinWidth(j);
      if (h->GetBinContent(j)!=0)
	errx2 += pow(h->GetBinError(j) * h->GetBinWidth(j),2);
    }
    h2->SetBinContent(i, sumx / h2->GetBinWidth(i));
    h2->GetBinError(i, sqrt(errx2) / h2->GetBinWidth(i));
  }

  return h2;
}

// Return results of interpolated histogram
// For interpolation, first turn it into log-log graph
TH1D *_href(0);
TGraphErrors *_gref(0);
Double_t fRef(Double_t *x, Double_t *p) {
  
  double pt = x[0];
  double eta = p[0];
  //double reset = p[1];

  //if (pt<10 || pt*cosh(eta)>6500.) return 0;
  if (pt*cosh(eta)>6500.*0.95) return 0;

  // New function
  assert(_href);
  if (_gref==0) { // create new reference graph

    _gref = new TGraphErrors(0);

    TF1 *f1 = new TF1("fRefTmp","pow(x,-5)*pow(1-2*x*cosh([0])/13000.,10)",
    		      5,6500./cosh(eta));
    f1->SetParameter(0,eta);

    for (int i = 1; i != _href->GetNbinsX()+1; ++i) {

      if (_href->GetBinContent(i)>0) {

	double xmin = _href->GetBinLowEdge(i);
	double xmax = _href->GetBinLowEdge(i+1);
	double yint = f1->Integral(xmin,xmax) / (xmax-xmin);
	double x = f1->GetX(yint, xmin, xmax, _prec);
	//double x = _href->GetBinCenter(i);
	double y = _href->GetBinContent(i);
	double ey = _href->GetBinError(i);

	int n = _gref->GetN();
	_gref->SetPoint(n, log(x), log(y));
	_gref->SetPointError(n, 0, log(y+ey)-log(y));
      }
    } // for i

    //delete f1;
  } // create new reference graph
  assert(_gref);

  double logxsec = _gref->Eval(log(pt));

  return exp(logxsec);
} // fRef


// Helper functions to find JEC for corrected pt
void setEtaPtRho(FactorizedJetCorrector *jec, double eta, double pt,
		 double rho){

  assert(jec);
  jec->setJetEta(eta);
  jec->setJetPt(pt);
  jec->setRho(rho);
  jec->setJetA(0.50265);

  return;
}

FactorizedJetCorrector *_thejec(0);
TF1 *fCorrPt(0);
Double_t funcCorrPt(Double_t *x, Double_t *p) {
  
  double eta = p[0];
  double pt = x[0];
  double rho = p[1];
  setEtaPtRho(_thejec, eta, pt, rho);

  return (_thejec->getCorrection() * pt);
}

const bool _useptgen = true;
const double _rhoDE = 19.21; // EOY17 DE jt450
double getJEC(FactorizedJetCorrector *jec, double eta, double pt, double rho) {

  setEtaPtRho(jec, eta, pt, rho);

  // if using pTgen, need to iterate to solve ptreco
  if (_useptgen) {

    double ptgen = pt;
    _thejec = jec;
    if (!fCorrPt) fCorrPt = new TF1("fCorrPt",funcCorrPt,5,6500,2);
    fCorrPt->SetParameters(eta, rho);
    // Find ptreco that gives pTreco*JEC = pTgen
    double ptreco = fCorrPt->GetX(ptgen,5,6500);

    setEtaPtRho(jec, eta, ptreco, rho);
  }

  return (jec->getCorrection());
} // getEtaPtE

const double ptmin(15);//114;//64;
const double ptmax(3500.);
const double emax(6500.);
const double xsecmin(1.0001e-7);
const double xsecmax(0.9999e11);
void drawDeltaJEC(string sy = "0.0-0.5", string sdir = "") {

  // string to char* for easier use with TForm later
  const char *cy = sy.c_str();

  // pdflatex cannot handle periods inside file name so remove them
  string sy2 = string(TString(cy).ReplaceAll(".",""));
  const char *cy2 = sy2.c_str();

  setTDRStyle();
  
  TDirectory *curdir = gDirectory;

  // Common format ROOT tuple from Engin Eren
  // https://gitlab.cern.ch/CMS-SMP-J/InclusiveJetsLegacy/blob/master/
  // => common2016.root (JEC V8?)
  // => common2016_October2018_V17.root (corrected for lumi?)
  // https://gitlab.cern.ch/lamartik/combinationfiles
  // => common17_V11.root
  // => common2018_V7.root
  // => common2018_V10.root
  //TFile *fin1 = new TFile("rootfiles/common2018_V7.root","READ");
  //TFile *fin1 = new TFile("rootfiles/common2018_V10.root","READ");
  //TFile *fin1 = new TFile("rootfiles/common2018_V13h.root","READ");
  //TFile *fin1 = new TFile("rootfiles/common2018_V10_hotzone.root","READ");
  //TFile *fin1 = new TFile("rootfiles/common2018_V10_hotzone-2.root","READ");
  //TFile *fin1 = new TFile("rootfiles/common2018_V18_unfolding.root","READ");
  //TFile *fin1 = new TFile("rootfiles/common2018_V18_unfolding.root","READ");
  //TFile *fin1 = new TFile("rootfiles/common2018_V19.root","READ");
  TFile *fin1 = new TFile("rootfiles/common2018_V19-5.root","READ");
  assert(fin1 && !fin1->IsZombie());
  
  //TFile *fin2 = new TFile("rootfiles/common2016_LegacyIOVs_v3.root","READ");
  //TFile *fin2 = new TFile("rootfiles/common2017_V11.root","READ");
  //TFile *fin2 = new TFile("rootfiles/common2017_V32.root","READ");
  //TFile *fin2 = new TFile("rootfiles/common2017_V32_tight.root","READ");
  //TFile *fin2 = new TFile("rootfiles/common2017_V32_hotzone.root","READ");
  //TFile *fin2 = new TFile("rootfiles/common2017_V32_hotzone-3.root","READ");
  //TFile *fin2 = new TFile("rootfiles/common2017_V32_hotzone_unfolding.root","READ");
  TFile *fin2 = new TFile("rootfiles/common2017_V32_hotzone_unfolding-3.root","READ");
  assert(fin2 && !fin2->IsZombie());

  //TFile *fin2ul = new TFile("rootfiles/commonUL2017_V1.root","READ");
  //TFile *fin2ul = new TFile("rootfiles/commonUL2017_V1-2.root","READ");
  //TFile *fin2ul = new TFile("rootfiles/commonUL2017_V1-3.root","READ");
  //TFile *fin2ul = new TFile("rootfiles/commonUL2017_V1-4.root","READ");
  //TFile *fin2ul = new TFile("rootfiles/commonUL2017_V2.root","READ");
  //TFile *fin2ul = new TFile("rootfiles/commonUL2017_V4.root","READ");
  //TFile *fin2ul = new TFile("rootfiles/commonUL2017_V4_V2M4res.root","READ");
  //TFile *fin2ul = new TFile("rootfiles/commonUL2017_V4_V2M4res_hotzone.root","READ");
  //TFile *fin2ul = new TFile("rootfiles/commonUL2017_V4_V2M4res_hotzone_scaled.root","READ");
  TFile *fin2ul = new TFile("rootfiles/commonUL2017_V4_V2M4res_hotzone_scaled2p.root","READ");
  assert(fin2ul && !fin2ul->IsZombie());

  //TFile *finh = new TFile("rootfiles/outHdata-Hdata-lowpu.root","READ");
  //TFile *finh = new TFile("rootfiles/outHdata-Hdata-lowpu-2.root","READ");
  //TFile *finh = new TFile("rootfiles/outH_l2l3data-Hdata-lowpu.root","READ");
  //TFile *finh = new TFile("rootfiles/outH_recorrected_l2l3data-Hdata-lowpu.root","READ");
  TFile *finh = new TFile("rootfiles/outH_defaultcorr-Hdata-lowpu.root","READ");
  assert(finh && !finh->IsZombie());

  TFile *finh2 = new TFile("rootfiles/outH_Fall17_v32-Hdata-lowpu.root","READ");
  assert(finh2 && !finh2->IsZombie());

  TFile *finh3 = new TFile("rootfiles/outH_recorrected_l2l3data-Hdata-lowpu.root","READ");
  assert(finh3 && !finh3->IsZombie());

  //TFile *fin3 = new TFile("rootfiles/common2016_October2018_V17.root","READ");
  //TFile *fin3 = new TFile("rootfiles/common2016_V11.root","READ");
  //TFile *fin3 = new TFile("rootfiles/common2016_V11_hotzone.root","READ");
  //TFile *fin3 = new TFile("rootfiles/common2016_V11_hotzone-2.root","READ");
  //TFile *fin3 = new TFile("rootfiles/common2016_V11_hotzone_unfolding.root","READ");
  TFile *fin3 = new TFile("rootfiles/common2016_V11_hotzone_unfolding-2.root","READ");
  assert(fin3 && !fin3->IsZombie());

  // Cross check of Suman's R-scan results
  //TFile *finx = new TFile("rootfiles/Differential_CrossSection_AK4_GH_Folded_Suman.root","READ");
  //TFile *finx = new TFile("rootfiles/Differential_CrossSection_AK4_EF_Folded_Suman.root","READ");
  //TFile *finx = new TFile("rootfiles/Differential_CrossSection_AK4_BCD_Folded_Suman.root","READ");
  TFile *finx = (getDetData ?
		 new TFile("rootfiles/Differential_CrossSection_AK4_BCDEFGH_Folded_Suman.root","READ") : // new
		 new TFile("rootfiles/Differential_CrossSection_AK4_BCDEFGH_UnFolded_Suman.root","READ")); // new file
  //TFile *finx = new TFile("rootfiles/Differential_CrossSection_AK4_GH_UnFolded_Suman.root","READ");
  //TFile *finx = new TFile("rootfiles/Differential_CrossSection_AK4_EF_UnFolded_Suman.root","READ");
  //TFile *finx = new TFile("rootfiles/Differential_CrossSection_AK4_BCD_UnFolded_Suman.root","READ");
  //TFile *finx = new TFile("rootfiles/Differential_CrossSection_AK4_BCDEFGH_UnFolded_Suman.root","READ"); // new file
  //TFile *finx = new TFile("rootfiles/Differential_CrossSection_AK4_Suman.root","READ"); // old file
  assert(finx && !finx->IsZombie());

  // MET filter efficiencies from Patrick
  //TFile *finmet1 = new TFile("rootfiles/Pythia16Flat_forMikko.root","READ");
  //assert(finmet1 && !finmet1->IsZombie());

  TFile *fu = new TFile("rootfiles/unfold.root","READ");
  assert(fu && !fu->IsZombie());

  TFile *frun2 = new TFile("rootfiles/drawDeltaJEC_Run2.root","READ");
  assert(frun2 && !frun2->IsZombie());

  // Current uncertainty
  const char *p = "CondFormats/JetMETObjects/data/";
  //const char *t18 = "Autumn18_RunD_V8_DATA";
  //const char *t18 = "Autumn18_RunD_V8_DATA";
  //const char *t18 = "Autumn18_V13temp_DATA";
  //const char *t18 = "Autumn18_V14_DATA";
  //const char *t18 = "Autumn18_V18_DATA";
  const char *t18 = "Autumn18_V19_MC";
  const char *a18 = "AK4";
  const char *s18 = Form("%s%s_Uncertainty_%sPFchs.txt",p,t18,a18);
  cout<<"**"<<s18<<endl<<flush;
  JetCorrectionUncertainty *unc18 = new JetCorrectionUncertainty(s18);

  const char *ss1 = Form("%s%s_UncertaintySources_%sPFchs.txt",p,t18,a18);
  const char *sss1 = "RelativeBal";//Sample";
  cout<<"**"<<ss1<<":"<<sss1<<endl<<flush;
  JetCorrectorParameters *ps1 = new JetCorrectorParameters(ss1,sss1);
  JetCorrectionUncertainty *uncs1 = new JetCorrectionUncertainty(*ps1);

  const char *ss2 = Form("%s%s_UncertaintySources_%sPFchs.txt",p,t18,a18);
  const char *sss2a = "RelativePtEC1";
  cout<<"**"<<ss2<<":"<<sss2a<<endl<<flush;
  JetCorrectorParameters *ps2a = new JetCorrectorParameters(ss2,sss2a);
  JetCorrectionUncertainty *uncs2a = new JetCorrectionUncertainty(*ps2a);
  const char *sss2b = "RelativePtEC2";
  cout<<"**"<<ss2<<":"<<sss2b<<endl<<flush;
  JetCorrectorParameters *ps2b = new JetCorrectorParameters(ss2,sss2b);
  JetCorrectionUncertainty *uncs2b = new JetCorrectionUncertainty(*ps2b);

  // 2017 uncertainty
  //const char *t17 = "Fall17_17Nov2017B_V11_DATA";
  const char *t17 = "Fall17_17Nov2017B_V32_DATA";
  const char *a17 = "AK4";
  const char *s17 = Form("%s%s_Uncertainty_%sPFchs.txt",p,t17,a17);
  cout<<"**"<<s17<<endl<<flush;
  JetCorrectionUncertainty *unc17 = new JetCorrectionUncertainty(s17);

  // 2016 uncertainty
  //const char *t16 = "Summer16_07Aug2017GH_V17_DATA";
  const char *t16 = "Summer16_07Aug2017GH_V11_DATA";
  const char *a16 = "AK4";
  const char *s16 = Form("%s%s_Uncertainty_%sPFchs.txt",p,t16,a16);
  cout<<"**"<<s16<<endl<<flush;
  JetCorrectionUncertainty *unc16 = new JetCorrectionUncertainty(s16);

  // 2015 uncertainty (50 ns)
  const char *t15 = "Summer15_50nsV5_DATA";
  const char *a15 = "AK4";
  const char *s15 = Form("%s%s_Uncertainty_%sPFchs.txt",p,t15,a15);
  cout<<"**"<<s15<<endl<<flush;
  JetCorrectionUncertainty *unc15 = new JetCorrectionUncertainty(s15);

  // 2012 uncertainty (8 TeV)
  const char *t12 = "Winter14_V8_DATA";
  const char *a12 = "AK5";
  const char *s12 = Form("%s%s_Uncertainty_%sPFchs.txt",p,t12,a12);
  cout<<"**"<<s12<<endl<<flush;
  JetCorrectionUncertainty *unc12 = new JetCorrectionUncertainty(s12);

  // Difference of L2 in SimpleL1 and ComplexL1 variants
  const char *tl2 = "Summer19UL17_V1";
  const char *ss = Form("%s%s_SimpleL1_MC_L2Relative_AK4PFchs.txt",p,tl2);
  cout << ss << endl;
  JetCorrectorParameters *par_l2S = new JetCorrectorParameters(ss);
  vector<JetCorrectorParameters> vparS;
  vparS.push_back(*par_l2S);
  FactorizedJetCorrector *jecS = new FactorizedJetCorrector(vparS);
  const char *sc = Form("%s%s_ComplexL1_MC_L2Relative_AK4PFchs.txt",p,tl2);
  cout << sc << endl;
  JetCorrectorParameters *par_l2C = new JetCorrectorParameters(sc);
  vector<JetCorrectorParameters> vparC;
  vparC.push_back(*par_l2C);
  FactorizedJetCorrector *jecC = new FactorizedJetCorrector(vparC);

  // List of eras to plot
  //const int neras = 3;//5;//12;//11;
  //string eras[neras] =
  // Directory to save plots (useful for IOV breakup)
  // const char *cdir = "";
  // Eras to be plotted
  string eras_run2[] =
    {"16LM",
     //"16BCD","16EF","16GH",
     //"17B","17C","17DE","17F",
     //"18A","18B","18C",
     //"18D",
     "17LM","18LM","Run2",
     "17UL"};
     //"16SC"};
  const int neras_run2 = sizeof(eras_run2)/sizeof(eras_run2[0]);
  //const char *cdir = "2016";
  string eras_2016[] = {"16BCD","16EF","16GH", "16LM"};//, "18LM"};
  const int neras_2016 = sizeof(eras_2016)/sizeof(eras_2016[0]);
  //const char *cdir = "2017";
  //string eras_2017[] = {"17B","17C","17DE","17F", "17LM"};//,"18LM"};
  string eras_2017[] = {"17B","17C","17DE","17F", "17LM"};//, "17H"};
  const int neras_2017 = sizeof(eras_2017)/sizeof(eras_2017[0]);
  //string eras_2017ul[] = {"17ULB","17ULC","17ULD","17ULE","17ULF","17UL"};
  string eras_2017ul[] = {"17UL","17ULB","17ULC","17ULD","17ULE","17ULF"};
  const int neras_2017ul = sizeof(eras_2017ul)/sizeof(eras_2017ul[0]);
  //string eras_2017h[] = {"17F","17ULF","17H","17H2","17H3","17LM","17UL"};
  string eras_2017h[] = {"17F","17ULF","17H","17H2","17H3"};
  const int neras_2017h = sizeof(eras_2017h)/sizeof(eras_2017h[0]);
  //const char *cdir = "2018";
  string eras_2018[] = {"18A","18B","18C","18D", "18LM"};//,"18LM"};
  const int neras_2018 = sizeof(eras_2018)/sizeof(eras_2018[0]);
    
  // Map eras to legend labels
  map<string,const char*> label;
  label["16BCD"] = "BCD";
  label["16EF"] = "EF";
  label["16GH"] = "GH";
  label["16LM"] = "2016";
  label["17B"] = "B";
  label["17C"] = "C";
  label["17DE"] = "DE";
  label["17F"] = "F";
  //label["17H"] = "H (low PU)";
  label["17LM"] = "2017";
  label["17ULB"] = "B (UL)";
  label["17ULC"] = "C (UL)";
  label["17ULD"] = "D (UL)";
  label["17ULE"] = "E (UL)";
  label["17ULF"] = "F (UL)";
  label["17H"] = "H (def)";//"H (low PU)";
  label["17H2"] = "H (V32)";//"H (low PU)";
  label["17H3"] = "H (UL17)";//"H (low PU)";
  label["17UL"] = "2017UL";
  label["18A"] = "A";
  label["18B"] = "B";
  label["18C"] = "C";
  label["18D"] = "D";
  label["18LM"] = "2018";
  label["16SC"] = "JME-19-003";
  label["Run2"] = "Run2";

  const char *cdir = sdir.c_str();
  string *eras(0);
  int neras(0);
  if (sdir=="")     { eras = eras_run2; neras = neras_run2; }
  if (sdir=="2016") { eras = eras_2016; neras = neras_2016; }
  if (sdir=="2017") { eras = eras_2017; neras = neras_2017; }
  if (sdir=="2017UL") { eras = eras_2017ul; neras = neras_2017ul; }
  if (sdir=="2017H") { eras = eras_2017h; neras = neras_2017h; }
  if (sdir=="2018") { eras = eras_2018; neras = neras_2018; }
  plotVs17UL = false;
  if (sdir=="17UL") { eras = eras_2017ul; neras = neras_2017ul; plotVs17UL = true; }

  // Mapping used later
  int irund(-1), irunsc(-1);//, irunbcd(-1);
  for (int i = 0; i != neras; ++i) {
    if (eras[i]=="18D") irund = i;
    if (eras[i]=="16SC") irunsc = i;
    //if (eras[i]=="16BCD") irunbcd = i;
    //if (eras[i]=="16LM" && irunbcd==-1) irunbcd = i;
    if (eras[i]=="18LM" && irund==-1) irund = i;
  }
  //assert(irund!=-1);
  //assert(irunbcd!=-1);

  // Updated mappings
  map<string, const char*> meras;
  meras["16BCD"] = "2016_BCD";
  meras["16EF"] = "2016_EF";
  meras["16GH"] = "2016_GH";
  meras["17B"] = "2017_B";
  meras["17C"] = "2017_C";
  meras["17DE"] = "2017_DE";
  meras["17F"] = "2017_F";
  meras["17ULB"] = "2017_B";
  meras["17ULC"] = "2017_C";
  meras["17ULD"] = "2017_D";
  meras["17ULE"] = "2017_E";
  meras["17ULF"] = "2017_F";
  meras["17H"] = "2017_H";
  meras["17H2"] = "2017_H";
  meras["17H3"] = "2017_H";
  meras["18A"] = "2018_A";
  meras["18B"] = "2018_B";
  meras["18C"] = "2018_C";
  meras["18D"] = "2018_D";
  meras["16SC"] = "2016";
  meras["16LM"] = "2016_all";
  meras["17LM"] = "2017_all";
  meras["17UL"] = "2017_all";
  meras["18LM"] = "2018_all";
  meras["Run2"] = "Run2";
  map<string, TFile*> mfins;
  mfins["16BCD"] = mfins["16EF"] = mfins["16GH"] = fin3;
  mfins["17B"] = mfins["17C"] = mfins["17DE"] = mfins["17F"] = fin2;
  mfins["17ULB"] = mfins["17ULC"] = mfins["17ULD"] = mfins["17ULE"] = mfins["17ULF"] = mfins["17UL"] = fin2ul;
  mfins["17H"] = finh;
  mfins["17H2"] = finh2;
  mfins["17H3"] = finh3;
  mfins["18A"] = mfins["18B"] = mfins["18C"] = mfins["18D"] = fin1;
  mfins["16SC"] = finx;
  mfins["16LM"] = fin3;
  mfins["17LM"] = fin2;
  mfins["18LM"] = fin1;
  mfins["Run2"] = 0;
  map<string, double> mlumi;
  mlumi["16BCD"] = mlumi["16EF"] = mlumi["16GH"] = 1;
  mlumi["17B"] = mlumi["17C"] = mlumi["17DE"] = mlumi["17F"] = 1;
  mlumi["17ULB"] = mlumi["17ULC"] = mlumi["17ULD"] = mlumi["17ULE"] = mlumi["17ULF"] = mlumi["17UL"] = 1;
  mlumi["17H"] = 1;//1e3;//1e-6;
  mlumi["17H2"] = 1;//1e3;//1e-6;
  mlumi["17H3"] = 1;//1e3;//1e-6;
  mlumi["18A"] = mlumi["18B"] = mlumi["18C"] = mlumi["18D"] = 1;
  mlumi["16SC"] = 1000;
  mlumi["16LM"] = mlumi["17LM"] = mlumi["18LM"] = mlumi["Run2"] = 1;
  map<string, int> mmarker;
  mmarker["16BCD"] = kOpenSquare;//kFullSquare;
  mmarker["16EF"] = kOpenCircle;//kOpenSquare;
  mmarker["16GH"] = kOpenDiamond;//kFullSquare;
  mmarker["17B"] = kOpenSquare;//kOpenCircle;
  mmarker["17C"] = kOpenCircle;
  mmarker["17DE"] = kOpenDiamond;//kFullCircle;
  mmarker["17F"] = kOpenStar;//kOpenCircle;
  mmarker["17ULB"] = kOpenSquare;
  mmarker["17ULC"] = kOpenCircle;
  mmarker["17ULD"] = kOpenCross;
  mmarker["17ULE"] = kOpenDiamond;
  mmarker["17ULF"] = kOpenStar;
  mmarker["17H"] = kFullCircle;//kFullStar;
  mmarker["17H2"] = kOpenCircle;//kOpenStar;
  mmarker["17H3"] = kFullDiamond;//kOpenStar;
  mmarker["18A"] = kOpenSquare;//kFullDiamond;
  mmarker["18B"] = kOpenCircle;//kOpenDiamond;
  mmarker["18C"] = kOpenDiamond;
  mmarker["18D"] = kOpenStar;//kFullDiamond;
  mmarker["16SC"] = kFullStar;
  mmarker["16LM"] = kFullSquare;
  mmarker["17LM"] = kFullCircle;
  mmarker["17UL"] = kOpenCircle;//kFullCross;
  mmarker["18LM"] = kFullDiamond;
  mmarker["Run2"] = kFullStar;
  map<string, int> mcolor;
  mcolor["16BCD"] = kMagenta+1;
  mcolor["16EF"] = kMagenta+2;//kCyan+2;
  mcolor["16GH"] = kMagenta+3;//kGray+2;
  mcolor["17B"] = kCyan+3;//kBlue;
  mcolor["17C"] = kCyan+2;
  mcolor["17DE"] = kBlue+0;//1;//kGreen+2;
  mcolor["17F"] = kBlue+0;//1;//kGreen+2;
  mcolor["17ULB"] = kCyan+3;
  mcolor["17ULC"] = kCyan+2;
  mcolor["17ULD"] = kBlue+2;
  mcolor["17ULE"] = kBlue+2;
  mcolor["17ULF"] = kBlue+2;
  mcolor["17H"] = kBlack;
  mcolor["17H2"] = kBlack;
  mcolor["17H3"] = kGray+3;
  mcolor["18A"] = kYellow+3;
  mcolor["18B"] = kOrange+2;
  mcolor["18C"] = kRed+1;
  mcolor["18D"] = kGray+2;//kBlack;
  mcolor["16SC"] = kRed+2;//kBlack;
  mcolor["16LM"] = kMagenta+1;
  mcolor["17LM"] = kBlue+1;
  mcolor["17UL"] = kBlue+2;
  mcolor["18LM"] = kBlack;
  mcolor["Run2"] = kRed;

  map<string, double> lumi;
  lumi["17ULBCDEF"] = 41.3;//41.5;
  lumi["17ULB"] = 4.8;
  lumi["17ULC"] = 9.6;
  lumi["17ULD"] = 4.2;
  lumi["17ULE"] = 9.3;
  lumi["17ULF"] = 13.4;

  /*
  const int nera = 11;
  TFile *fins[nera] =
    {fin1, fin1, fin1, fin1,
     fin2, fin2, fin2, fin2, //fin2,
     fin3, fin3, fin3};
  string eras[nera] =
    //{"A","B","C", "D",
    {"2018_A","2018_B","2018_C", "2018_D",
     //"2017_B","2017_C","2017_D","2017_E","2017_F",
     "2017_B","2017_C","2017_DE","2017_F",
     //"BCD2016","EF2016","GH2016"};
     //"2016_B","2016_Fe","2016_G"};
     "2016_BCD","2016_EF","2016_GH"};
  // luminosity re-normalization, if any needed
  double k17 = 1e5;
  double k16 = 0.5;
  double lumi[nera] =
    //{14.0/59.9, 7.1/59.9, 6.9/59.9, 31.9/59.9,
    {1, 1, 1, 1,
     //4800, 9600, 4200, 9300, 13400,
     1, 1, 1, 1, //4800, 9600, 13500, 13400,
     1, 1, 1};
     //4.8*k17/41.4, 9.6*k17/41.4, 4.2*k17/41.4, 9.3*k17/41.4, 13.4*k17/41.4};
  //lumimap["A"] = "Run2018A 14.0 fb^{-1}"; //PdmV Analysis TWiki
  //lumimap["B"] = "Run2018B 7.1 fb^{-1}"; //PdmV Analysis TWiki
  //lumimap["C"] = "Run2018C 6.9 fb^{-1}"; //PdmV Analysis TWiki
  //lumimap["D"] = "Run2018D 31.9 fb^{-1}"; //PdmV Analysis TWiki
  //lumimap["ABC"] = "Run2018ABC 28.0 fb^{-1}"; //PdmV Analysis TWiki
  //lumimap["ABCD"] = "Run2018ABCD 59.9 fb^{-1}"; //PdmV Analysis TWiki

  int color[nera] =
    {kRed+2, kOrange+2, kBlue+1, kBlack,
     kRed+2, kOrange+2, kBlue+1, kBlack,
     kMagenta+2, kCyan+2, kGray+2};
  int marker[nera] =
    {kFullCircle, kFullDiamond, kFullStar, kFullSquare,
     kOpenCircle, kOpenDiamond, kOpenStar, kOpenSquare,
     kFullCircle, kFullDiamond, kFullStar};
  const char* label[nera] =
    //{"RunA","RunB","RunC","RunD"};
    {"18A","18B","18C","18D",
     //"17B","17C","17D","17E","17F",
     "17B","17C","17DE","17F",
     "16BCD","16EF","16GH"};
*/


  // Settings for the spectrum fit and JEC-equivalent plot y-axis range
  double eta(0.), ymin(-4), ymax(+6);
  /*
  if (sy=="0.0-0.5") { eta = 0;   ymin = -6; ymax = +9; }
  if (sy=="0.5-1.0") { eta = 0.5; ymin = -6; ymax = +9; }
  if (sy=="1.0-1.5") { eta = 1.0; ymin = -6; ymax = +9; }
  if (sy=="1.5-2.0") { eta = 1.5; ymin = -10; ymax = +15; }
  if (sy=="2.0-2.5") { eta = 2.0; ymin = -10; ymax = +15; }
  if (sy=="2.5-3.0") { eta = 2.5; ymin = -25; ymax = +45; }
  if (sy=="3.2-4.7") { eta = 3.3; ymin = -25; ymax = +45; }
  */
  if (sy=="0.0-1.3") { eta = 0;   ymin = -6; ymax = +6; }
  if (sy=="0.0-0.5") { eta = 0;   ymin = -6; ymax = +6; }
  if (sy=="0.5-1.0") { eta = 0.5; ymin = -6; ymax = +6; }
  if (sy=="1.0-1.5") { eta = 1.0; ymin = -6; ymax = +6; }
  if (sy=="1.5-2.0") { eta = 1.5; ymin = -9; ymax = +9; }
  if (sy=="2.0-2.5") { eta = 2.0; ymin = -9; ymax = +9; }
  if (sy=="2.5-3.0") { eta = 2.5; ymin = -18; ymax = +18; }
  // HF likely has slight rapidity bias affecting minimum |y|
  // This is very important at high pT near kinematic limit
  // Or rapidity smearing, high pT JER tails above kinematic limit
  if (sy=="3.2-4.7") { eta = 3.3; ymin = -12; ymax = +12; }
  // Drop out 16SC for higher eta
  int nera = (eta>1.9 && irunsc!=-1 ? neras-1 : neras);

  ///////////////////////////////////////////////////////////////////

  // Load reference MC
  TH1D *hmc(0);
  TH1D *hmc2(0);
  TH1D *hmc3(0);
  if (plotMC) {
    //TFile *fmc = new TFile("../jecsys2018/rootfiles/common2016_V11_hotzone-2.root","READ");
    TFile *fmc = new TFile("rootfiles/common2018_V19.root","READ");
    //TFile *fmc = new TFile("rootfiles/common2018_V19-5.root","READ");
    //TFile *fmc = new TFile("../jecsys2018/rootfiles/output-MC-1-Fall18V8-D.root","READ");
    assert(fmc && !fmc->IsZombie());
    //TFile *fmc2 = new TFile("rootfiles/P8_dijet_5000000_ptw.root","READ");
    TFile *fmc2 = new TFile("rootfiles/pThat10_P8_dijet_5M_pt.root","READ");
    assert(fmc2 && !fmc2->IsZombie());
    //TFile *fmc3 = new TFile("rootfiles/pThat5_P8_dijet_5000000_pt.root","READ");
    TFile *fmc3 = new TFile("rootfiles/pThat5_P8_dijet_5M_pt.root","READ");
    assert(fmc3 && !fmc3->IsZombie());
    //hmc = (TH1D*)fmc->Get(Form("ak4/Eta_%s/hpt_MCP8M1_2016_GH_gen",cy));
    hmc = (TH1D*)fmc->Get(Form("ak4/Eta_%s/hpt_P8CP5_2018_A_gen",cy)); //V19
    if (!hmc && eta==0)
      hmc = (TH1D*)fmc->Get(Form("ak4/Eta_%s/hpt_P8CP5_2018_A_gen","0.5-1.0"));
    assert(hmc);
    //hmc = (TH1D*)fmc->Get(Form("Standard/Eta_%s/mc/hpt_g0",cy));
    hmc2 = (TH1D*)fmc2->Get(Form("Eta_%s",cy)); // ptw
    assert(hmc2);
    hmc3 = (TH1D*)fmc3->Get(Form("Eta_%s",cy)); // ptw
    assert(hmc3);

    hmc->Scale(5000.); // A_gen
    //hmc->Scale(5000/0.8,"width"); // hpt_g0
    //hmc->SetMarkerStyle(kOpenDiamond);
    //hmc->SetMarkerColor(kGreen+3);
    //hmc->SetLineColor(kGreen+3);

    //hmc2 = (TH1D*)hmc2->Clone(Form("hmc_%s",cy));
    //hmc->Scale(eta<3.2 ? 0.8 : 0.6);
    //hmc2->Scale(0.7e4,"width"); // random guess
    double sumwp8m2pt10 = 985699;
    //double xsecp8cp5pt10 = 1820000000*pow(10./15.,-5)*5./15.+1820000000+138900000+19100000+2735000+467500+117400+7753+642.1+185.9+32.05+9.365+0.8398+0.1124+0.006752+0.0001626;
    //hmc2->Scale(xsecp8cp5pt10/sumwp8m2pt10,"width");
    double xsecp8m2pt10 = 8.4271780e+09;
    //  = (8.42987e+09+8.43093e+09+8.42949e+09+8.42071e+09+8.42489e+09)/5.
    double kp8m2pt10 = 0.8;
    hmc2->Scale(kp8m2pt10*xsecp8m2pt10/sumwp8m2pt10,"width");
    if (eta>=3.2) hmc2->Scale(1./3.);
    if (sy=="0.0-1.3") hmc2->Scale(1./2.6);

    double sumwp8m2pt5 = 5.87183e+06;
    double xsecp8m2pt5 = 7.2921120e+10;
    // =(7.28676e+10+7.29845e+10+7.29394e+10+7.29014e+10+7.29127e+10)/5.
    double kp8m2pt5 = 0.8;
    hmc3->Scale(kp8m2pt5*xsecp8m2pt5/sumwp8m2pt5,"width");
    if (eta>=3.2) hmc3->Scale(1./3.);
    if (sy=="0.0-1.3") hmc3->Scale(1./2.6);

    //if (eta>=3.2) hmc2->Scale(1./3.); // x3 wider rapidity bin
    //hmc2->SetMarkerStyle(kOpenDiamond);
    //hmc2->SetMarkerColor(kBlue);
    //hmc2->SetLineColor(kBlue);

    // Apply ad-hoc MET filter efficiency correction to MC
    // Also remove low pT (<30-50 GeV) biased by pThat>15 GeV
    // (until somehow fixed)
    /*
    if (correctFilterEff) {
      
      TH2D *h2b = (TH2D*)finmet1->Get("before"); assert(h2b);
      TH2D *h2a = (TH2D*)finmet1->Get("after"); assert(h2a);
      int iy = int((eta+0.25)/0.5)+1;
      TH1D *hb = h2b->ProjectionX(Form("hbmc_%d",iy),iy,iy);
      TH1D *ha = h2a->ProjectionX(Form("hamc_%d",iy),iy,iy);
      TH1D *hr = (TH1D*)ha->Clone(Form("hrmetmc_%d",iy));
      hr->Divide(hb);
      
      for (int i = 1; i != hmc->GetNbinsX()+1; ++i) {
	if (hmc->GetBinContent(i)!=0) {
	  double pt = hmc->GetBinCenter(i);
	  double rmet = hr->GetBinContent(hr->FindBin(pt));
	  if (eta<3.2) { // only available for barrel
	    hmc->SetBinContent(i, rmet ? hmc->GetBinContent(i)/rmet : 0);
	    hmc->SetBinError(i, rmet ? hmc->GetBinError(i)/rmet : 0);
	  }
	  if (pt<30.) {
	    hmc->SetBinContent(i, 0);
	    hmc->SetBinError(i, 0);
	  }
	}
      }
    } // ad-hoc MET filter
    */
  } // load ref MC

  // Load reference data
  TH1D *href(0);
  TF1 *fref = new TF1("fref",
		      Form("[0]*pow(x,[1]+[2]*log(x))"
			   "*pow(1-2.*x*cosh(%1.2f)/13000.,[3])",eta),
		      0.85*ptmin,emax/cosh(eta));
  fref->SetParameters(1e12,-5,0.001,10);
  if (true) {
    //href = (TH1D*)fin2ul->Get(Form("ak4/Eta_%s/hpt_data_2017_all_ptcl_fwd",cy));
    //href = (TH1D*)fin2->Get(Form("ak4/Eta_%s/hpt_data_2017_all_ptcl_fwd",cy));
    href = (TH1D*)frun2->Get(Form("hrun2_%s",cy));
    //if (!href && sy=="0.0-1.3") { // PATCH
    if (plotVs17UL && sdir == "17UL") {
      //TH1D *href1 = (TH1D*)frun2->Get("hrun2_0.0-0.5"); assert(href1);
      //TH1D *href2 = (TH1D*)frun2->Get("hrun2_0.5-1.0"); assert(href2);
      //TH1D *href3 = (TH1D*)frun2->Get("hrun2_1.0-1.5"); assert(href3);
      //href = (TH1D*)href1->Clone("href");
      //href->Add(href2);
      //href->Add(href3);
      //href->Scale(1./3);
      //href = (TH1D*)fin2ul->Get("ak4/Eta_0.0-1.3/hpt_data_2017_all_ptcl_fwd");
      href = (TH1D*)fin2ul->Get(Form("ak4/Eta_%s/hpt_data_2017_all_%s",cy,
				     getDetData ?
				     (unfoldData ? "ptcl_fwd" : "det") :
				     (getFwd ? "ptcl_fwd" : "ptcl_dag")));
      assert(href);
      href = (TH1D*)href->Clone(Form("href_%s",cy));
      // Patch reference for own forward smearing
      if (getDetData && unfoldData) {
	TFile *f = new TFile("pdf/drawDeltaJEC_UL17patch.root","READ");
	//if (f && !f->IsZombie()) {
	assert(f && !f->IsZombie());
	TH1D *hs = (TH1D*)f->Get("jetRatio_BCDEF"); assert(hs);
	// Divide by hand not to change href uncertainties
	for (int i = 1; i != href->GetNbinsX()+1; ++i) {
	  href->SetBinContent(i, href->GetBinContent(i)*hs->GetBinContent(i));
	}
	//}
      }
    }
    assert(href);
    href->Fit(fref,"QRN");
  }
  //fref->SetRange(ptmin*0.85,emax/cosh(eta));

  // Load new results
  //TH1F *h1s[nera];
  //map<string, TH1F*> mh1s;
  TH1D *h1s[nera];
  map<string, TH1D*> mh1s;
  for (int iera = 0; iera != nera; ++iera) {

    //const char *cera = eras[iera].c_str(); 
    string se = eras[iera];
    const char *cera = meras[se]; 

    // Combine all of Run 2 data after efficiency corrections and unfolding
    if (se=="Run2") {

      double lum16 = 1;
      double lum17 = 1;
      double lum18 = 1;
      double lumtot = lum16 + lum17 + lum18;
      //TH1F *h16 = mh1s["16LM"]; assert(h16);
      //TH1F *h17 = mh1s["17LM"]; assert(h17);
      //TH1F *h18 = mh1s["18LM"]; assert(h18);
      TH1D *h16 = mh1s["16LM"]; assert(h16);
      TH1D *h17 = mh1s["17LM"]; assert(h17);
      TH1D *h18 = mh1s["18LM"]; assert(h18);
      TH1D *hrun2 = (TH1D*)h16->Clone(Form("hrun2_%s",cy));
      hrun2->Reset();
      hrun2->Add(h16,lum16);
      hrun2->Add(h17,lum17);
      hrun2->Add(h18,lum18);
      hrun2->Scale(1./lumtot);
      //h1s[iera] = (TH1F*)hrun2;
      //mh1s[se] = (TH1F*)hrun2;
      h1s[iera] = (TH1D*)hrun2;
      mh1s[se] = (TH1D*)hrun2;

      // Store results to a file for later reference
      TDirectory *curdir = gDirectory;
      TFile *fout = new TFile("pdf/drawDeltaJEC_Run2.root","UPDATE");
      hrun2->Write(hrun2->GetName(), TObject::kOverwrite);
      fout->Close();
      curdir->cd();

      continue; // don't process corrections etc. again
    }

    //string hname = Form("ak4/Eta_%s/hpt_data_%s_det",cy,cera);
    string hname = (getDetData ? Form("ak4/Eta_%s/hpt_data_%s_det",cy,cera) :
		    getFwd ? Form("ak4/Eta_%s/hpt_data_%s_ptcl_fwd",cy,cera) :
		    Form("ak4/Eta_%s/hpt_data_%s_ptcl_dag",cy,cera));
    //TH1F *hera = (TH1F*)fins[iera]->Get(hname.c_str());
    //TH1F *hera = (TH1F*)mfins[se]->Get(hname.c_str());
    TH1D *hera = (TH1D*)mfins[se]->Get(hname.c_str());

    // 2017H only has det and different binning
    if (!hera && (se=="17H" || se=="17H2" || se=="17H3")) {
      string hname = Form("ak4/Eta_%s/hpt_data_%s_det",cy,cera);
      //hera = (TH1F*)mfins[se]->Get(hname.c_str());
      hera = (TH1D*)mfins[se]->Get(hname.c_str());

      if (hera) {
	// PATCH missing pT bin width normalization
	//for (int i = 1; i != hera->GetNbinsX()+1; ++i) {
	//hera->SetBinContent(i, hera->GetBinContent(i)/hera->GetBinWidth(i));
	//hera->SetBinError(i, hera->GetBinError(i)/hera->GetBinWidth(i));
	//}
	// PATCH different pT binning, scale by 0.5
	//TH1F *htmp = (TH1F*)mh1s["17B"]; assert(htmp);
	//htmp = (TH1F*)htmp->Clone("17H");
	//TH1D *htmp = (TH1D*)mh1s["17B"]; assert(htmp);
	TH1D *htmp = (TH1D*)mh1s["17F"]; assert(htmp);
	htmp = (TH1D*)htmp->Clone(se.c_str());
	htmp->Reset();
	for (int i = 1; i != hera->GetNbinsX()+1; ++i) {
	  int j = htmp->FindBin(hera->GetBinCenter(i));
	  //htmp->SetBinContent(j, 0.5*hera->GetBinContent(i));
	  //htmp->SetBinError(j, 0.5*hera->GetBinError(i));
	  //htmp->SetBinContent(j, 0.5*hera->GetBinContent(i));
	  //htmp->SetBinError(j, 0.5*hera->GetBinError(i));
	  htmp->SetBinContent(j, hera->GetBinContent(i));
	  htmp->SetBinError(j, hera->GetBinError(i));
	}
	hera = htmp;
      }

    }

    // Engin's file has different naming scheme
    // Same with Suman's file
    if (!hera) {
      int iy = int((eta+0.25)/0.5)+1;
      hname = Form("ak4/y_%s/hptData_%s_detector_%dbin",cy,cera,iy);
      if (iy>4)
	//hera = (TH1F*)mfins[se]->Get("ak4/y_1.5-2.0/hptData_2016_detector_4bin");
	hera = (TH1D*)mfins[se]->Get("ak4/y_1.5-2.0/hptData_2016_detector_4bin");
      else
	//hera = (TH1F*)mfins[se]->Get(hname.c_str());
	hera = (TH1D*)mfins[se]->Get(hname.c_str());
    }
    if (!hera) cout << "Histogram " << hname << " not found!" << endl << flush;
    assert(hera);

    // Clone histogram so not overwriting later
    hera = (TH1D*)hera->Clone(Form("hera_%s",se.c_str()));
				   
    // patch/fix a strange normalization bug for 16BCD, 16EF and 16GH
    //string se = label[iera];
    //string se = eras[iera];
    //if (se=="16BCD"||se=="16EF"||se=="16GH") {//||se=="17DE") {
    //for (int i = 1; i != hera->GetNbinsX()+1; ++i) {
    //double k = int((eta+0.25)/0.5)+1;
    //hera->SetBinContent(i, hera->GetBinContent(i)*k);
    //hera->SetBinError(i, hera->GetBinError(i)*k);
    //}
    //}

    // Correct ECAL prefire for 2016 and 2017 data at 2<|eta|<3
    if (correctECALprefire ||
	(correctECALprefire17H && (se=="17H" || se=="17H2" || se=="17H3"))) {

      jer_iov run(run2018);
      //if (se=="16BCD"||se=="16EF"||se=="16GH")         run = run2016;
      //if (se=="17B"||se=="17C"||se=="17DE"||se=="17F") run = run2017;
      if (se=="16BCD") run = run2016bcd;
      if (se=="16EF")  run = run2016ef;
      if (se=="16GH")  run = run2016gh;
      if (se=="16LM")  run = run2016;
      if (se=="17B")   run = run2017b;
      if (se=="17C")   run = run2017c;
      if (se=="17DE")  run = run2017de;
      if (se=="17F")   run = run2017f;
      if (se=="17H")   run = run2016bcd;//run2017b;//run2017f; // LowPU
      if (se=="17H2")  run = run2016bcd;//run2017b;//run2017f; // LowPU
      if (se=="17H3")  run = run2016bcd;//run2017b;//run2017f; // LowPU
      if (se=="17LM")  run = run2017;
      if (se=="17ULB")  run = run2017b;
      if (se=="17ULC")  run = run2017c;
      if (se=="17ULD")  run = run2017de;
      if (se=="17ELE")  run = run2017de;
      if (se=="17ULF")  run = run2017f;
      if (se=="17UL")   run = run2017de;//run = run2017;
      if (se=="18A")   run = run2018abc;
      if (se=="18B")   run = run2018abc;
      if (se=="18C")   run = run2018abc;
      if (se=="18D")   run = run2018d;
      if (se=="18LM")  run = run2018;

      for (int i = 1; i != hera->GetNbinsX()+1; ++i) {

	double pt = hera->GetBinCenter(i);
	double ineff = ecalprefire(pt, eta+0.1, run);
	double eff = 1 - ineff;
	double corr = 1./eff;
	hera->SetBinContent(i, hera->GetBinContent(i)*corr);
	hera->SetBinError(i, hera->GetBinError(i)*corr);
      }
    } // correctECALprefire

    // Apply ad-hoc unfolding to det-level data
    //if (true && se!="16SC") {
    //if (true) {
    //if (unfoldData) {
    if (unfoldData || (unfold17H && (se=="17H" || se=="17H2" || se=="17H3"))) {
      //TF1 *f1 = (TF1*)fu->Get(Form("fr_%s",cy)); assert(f1);
      //TH1D *hr = (TH1D*)fu->Get(Form("hr_%s",cy)); assert(hr);
      int year(0);
      if (TString(se.c_str()).Contains("16")) year = 2016;
      if (TString(se.c_str()).Contains("17")) year = 2017; // except...
      if (TString(se.c_str()).Contains("17H")) year = 2016; // ...if 2017H
      //if (TString(se.c_str()).Contains("18")) year = 2018;
      if (se=="18A") year = 2019;//2018;//abc;
      if (se=="18B") year = 2019;//2018;//abc;
      if (se=="18C") year = 2019;//2018;//abc;
      if (se=="18D") year = 2020;//2019;//2018d;
      if (se=="18LM") year = 2018;//abc;
      // separate eras for each UL
      if (se=="17UL") year = 2021;
      if (se=="17ULB") year = 2022;
      if (se=="17ULC") year = 2023;
      if (se=="17ULD") year = 2024;
      if (se=="17ULE") year = 2025;
      if (se=="17ULF") year = 2026;

      // special case for 2018 2.5-3.0 bin
      //if (TString(se.c_str()).Contains("18") && sy=="2.5-3.0") year = 2016;
      // PATCH special case for 2018 V13h
      //if (TString(se.c_str()).Contains("18") && sy=="2.5-3.0") year = 2016;
      // PATCH Suman's data
      //if (TString(se.c_str()).Contains("16SC")) year = 2016;
      //if (se=="18D" && sy=="2.5-3.0") year = 2016;
      TH1D *hr = (TH1D*)fu->Get(Form("hr_%s_%d",cy,year)); assert(hr);
      for (int i = 1; i != hera->GetNbinsX()+1; ++i) {
	if (hera->GetBinContent(i)!=0) {
	  double pt = hera->GetBinCenter(i);
	  //double rdu = f1->Eval(pt);
	  double rdu = hr->GetBinContent(hr->FindBin(pt));
	  hera->SetBinContent(i, rdu ? hera->GetBinContent(i)/rdu : 0);
	  hera->SetBinError(i, rdu ? hera->GetBinError(i)/rdu : 0);
	}
      }
    } // ad-hoc unfolding

    // Apply ad-hoc filter efficiency correction to data
    /*
    if (correctFilterEff && fabs(eta)<2.0 && se!="16SC" &&
	TString(se.c_str()).Contains("16")) {
      
      TH2D *h2b = (TH2D*)finmet1->Get("before"); assert(h2b);
      TH2D *h2a = (TH2D*)finmet1->Get("after"); assert(h2a);
      int iy = int((eta+0.25)/0.5)+1;
      TH1D *hb = h2b->ProjectionX(Form("hb_%d",iy),iy,iy);
      TH1D *ha = h2a->ProjectionX(Form("ha_%d",iy),iy,iy);
      TH1D *hr = (TH1D*)ha->Clone(Form("hrmet_%d",iy));
      hr->Divide(hb);
      
      for (int i = 1; i != hera->GetNbinsX()+1; ++i) {
	if (hera->GetBinContent(i)!=0) {
	  double pt = hera->GetBinCenter(i);
	  double rmet = hr->GetBinContent(hr->FindBin(pt));
	  hera->SetBinContent(i, rmet ? hera->GetBinContent(i)/rmet : 0);
	  hera->SetBinError(i, rmet ? hera->GetBinError(i)/rmet : 0);
	}
      }
    } // ad-hoc MET filter
    */

    //hera->Scale(1./lumi[iera]);
    hera->Scale(1./mlumi[se]);
    h1s[iera] = hera;
    mh1s[se] = hera;
  } // for iera

  // Fit new results with powerlaw
  TF1 *f1s[nera];
  for (int iera = 0; iera != nera; ++iera) {
    TF1 *f1 = new TF1(Form("f1_%d",iera),
		      Form("[0]*pow(x,[1]+[2]*log(x))"
			   "*pow(1-2.*x*cosh(%1.2f)/13000.,[3])",eta),
		      ptmin,emax/cosh(eta));
    f1->SetParameters(1e12,-5,0.001,10);
    f1->SetNpx(2640.*2); // this limits precision of DeltaJEC
    h1s[iera]->Fit(f1,"RNI");
    //f1->SetLineColor(color[iera]);
    f1->SetLineColor(mcolor[eras[iera]]);
    f1s[iera] = f1;
  } // for iera


  ///////////////////////////////////////////////////////////////////

  // Unfolded reference data from SMP-15-007 stored in HEPDATA
  //TH1F *h50ns4(0), *h50ns7(0);
  TH1D *h50ns4(0), *h50ns7(0);
  if (plot2015) {
    TFile *fhepd = new TFile("rootfiles/HEPData-ins1459051-v1-root.root",
			     "READ");
    assert(fhepd && !fhepd->IsZombie());
    curdir->cd();

    // HEPData (1-7 AK7, 8-14 AK4; both unfolded)
    map<string, string> mhep7, mhep4;
    mhep7["0.0-0.5"] = "Table 1";
    mhep7["0.5-1.0"] = "Table 2";
    mhep7["1.0-1.5"] = "Table 3";
    mhep7["1.5-2.0"] = "Table 4";
    mhep7["2.0-2.5"] = "Table 5";
    mhep7["2.5-3.0"] = "Table 6";
    mhep7["3.2-4.7"] = "Table 7";
    //
    mhep4["0.0-0.5"] = "Table 8";
    mhep4["0.5-1.0"] = "Table 9";
    mhep4["1.0-1.5"] = "Table 10";
    mhep4["1.5-2.0"] = "Table 11";
    mhep4["2.0-2.5"] = "Table 12";
    mhep4["2.5-3.0"] = "Table 13";
    mhep4["3.2-4.7"] = "Table 14";
    //
    const char *ct7 = mhep7[sy].c_str();
    //TH1F *hhepdy7 = (TH1F*)fhepd->Get(Form("%s/Hist1D_y1",ct7));
    TH1D *hhepdy7 = (TH1D*)fhepd->Get(Form("%s/Hist1D_y1",ct7));
    assert(hhepdy7);
    //TH1F *hhepde7 = (TH1F*)fhepd->Get(Form("%s/Hist1D_y1_e1",ct7));
    TH1D *hhepde7 = (TH1D*)fhepd->Get(Form("%s/Hist1D_y1_e1",ct7));
    assert(hhepde7);
    for (int i = 1; i != hhepdy7->GetNbinsX()+1; ++i) {
      hhepdy7->SetBinError(i, hhepde7->GetBinContent(i));
    }
    h50ns7 = hhepdy7;
    //
    const char *ct4 = mhep4[sy].c_str();
    //TH1F *hhepdy4 = (TH1F*)fhepd->Get(Form("%s/Hist1D_y1",ct4));
    TH1D *hhepdy4 = (TH1D*)fhepd->Get(Form("%s/Hist1D_y1",ct4));
    assert(hhepdy4);
    //TH1F *hhepde4 = (TH1F*)fhepd->Get(Form("%s/Hist1D_y1_e1",ct4));
    TH1D *hhepde4 = (TH1D*)fhepd->Get(Form("%s/Hist1D_y1_e1",ct4));
    assert(hhepde4);
    for (int i = 1; i != hhepdy4->GetNbinsX()+1; ++i) {
      hhepdy4->SetBinError(i, hhepde4->GetBinContent(i));
    }
    //h50ns4 = hhepdy4;
    //h50ns4 = (TH1F*)href;
    h50ns4 = (TH1D*)href;
    
    // Copy AK4 over to new binning, and set 30<pT<40 GeV from Run18D (was 18A) 
    // Also set/replace emax>2000 from Run18D to constrain very high pT
    // Constraints set at 30% level (xerr)
    // => replace Run2018D with Run2016BCD at low pT to reduce ZB bias
    /*
    h50ns4 = new TH1F("h50ns4",";p_{T} (GeV);Cross section",nptb,ptbins);
    for (int i = 1; i != h50ns4->GetNbinsX()+1; ++i) {
      int x = h50ns4->GetBinCenter(i);
      int j = hhepdy4->FindBin(x);
      //int irund = 3; // Index of 2018D
      //int irunbcd = 8; // Index of 2016BCD
      //int irund = 1;//3;//10; // Index of 2018D
      //int irunbcd = 0; // Index of 2016BCD
      //assert(string(label[irunbcd])=="16BCD");
      //assert(string(label[irund])=="18D");
      //assert(eras[irunbcd]=="16BCD" || eras[irunbcd]=="16LM");
      assert(irund!=-1);
      assert(eras[irund]=="18D" || eras[irund]=="18LM" || eras[irund]=="Run2");
      const double emax = 2000;
      //const double xerr = 0.05;//0.3; // 50%
      double xerrhi = 0.01;
      if (eta==2.0 || eta==2.5) xerrhi = 0.1;
      double xerrlo = 0.10;
      //int k1 = h1s[irunbcd]->FindBin(x);
      int k1 = h1s[irund]->FindBin(x);
      int k2 = h1s[irund]->FindBin(x);
      if (fabs(eta)>=3.2) { // Use 2016 for HF reference => 2018
	double pt = h50ns4->GetBinCenter(i);
	double xerr = 0.15;//(pt<300 ? 0.3 : 0.2);//xerrlo;
	int irun =  irund;//irunbcd;
	int k = k1;
	h50ns4->SetBinContent(i, h1s[irun]->GetBinContent(k));
	h50ns4->SetBinError(i, sqrt(pow(h1s[irun]->GetBinError(k),2)
				    +pow(h1s[irun]->GetBinContent(k)*xerr,2)));
	if (pt<45) {
	  h50ns4->SetBinContent(i, 0);
	  h50ns4->SetBinError(i, 0);
	}
      }
      else if (hhepdy4->GetBinContent(j)!=0 && (x*cosh(eta)<emax)) {
	h50ns4->SetBinContent(i, hhepdy4->GetBinContent(j));
	h50ns4->SetBinError(i, hhepdy4->GetBinError(j));
      }
      else if ((x>=30 && x<=40) || (x*cosh(eta)>=emax)) {
	//int irun =  (x<=40 ? irunbcd : irund);
	int irun =  (x<=40 ? irund : irund);
	int k = (x<=40 ? k1 : k2);
	double xerr = (x<=40 ? xerrlo : xerrhi);
	h50ns4->SetBinContent(i, h1s[irun]->GetBinContent(k));
	h50ns4->SetBinError(i, sqrt(pow(h1s[irun]->GetBinError(k),2)
				    +pow(h1s[irun]->GetBinContent(k)*xerr,2)));
      }
    }
    */
  } // plot2015
  assert(h50ns4 || !plot2015);
  assert(h50ns7 || !plot2015);


  // Fit old AK7 results with powerlaw
  TF1 *f50ns7 = new TF1("f50ns7",
			Form("[0]*pow(x,[1]+[2]*log(x))"
			     "*pow(1-2.*x*cosh(%1.2f)/13000.,[3])",eta),
			ptmin,emax/cosh(eta));
  f50ns7->SetParameters(1e12,-5,0.001,10);
  f50ns7->SetNpx(2640.*2); // this limits precision of DeltaJEC
  if (plot2015) h50ns7->Fit(f50ns7,"RNI");
  f50ns7->SetLineColor(kGreen+2);

  // Fit old AK4 results with powerlaw
  TF1 *f50ns4 = new TF1("f50ns4",
			Form("[0]*pow(x,[1]+[2]*log(x))"
			     "*pow(1-2.*x*cosh(%1.2f)/13000.,[3])",eta),
			ptmin,emax/cosh(eta));
  f50ns4->SetParameters(1e12,-5,0.001,10);
  //TF1 *f50ns4 = new TF1("f50ns4",
  //		Form("[0]*pow(x,[1]+[2]*log(x)+[4]*log(x)*log(x))"
  //		     "*pow(1-2.*x*cosh(%1.2f)/13000.,[3])",eta),
  //		ptmin,emax/cosh(eta));
  //f50ns4->SetParameters(1e12,-5,0.001,10,0.01);
  f50ns4->SetNpx(2640.*2); // this limits precision of DeltaJEC
  if (plot2015) h50ns4->Fit(f50ns4,"RNI");
  f50ns4->SetLineColor(kGreen+2);

  // Fit gen MC with powerlaw
  assert(hmc || !plotMC);
  assert(hmc2 || !plotMC);
  assert(hmc3 || !plotMC);
  TF1 *fmc = new TF1("fmc",
		     Form("[0]*pow(x,[1]+[2]*log(x))"
			  "*pow(1-2.*x*cosh(%1.2f)/13000.,[3])",eta),
			ptmin,emax/cosh(eta));
  fmc->SetParameters(1e12,-5,0.001,10);
  fmc->SetNpx(2640.*2); // this limits precision of DeltaJEC
  //hmc->Fit(fmc,"RNI");
  //if (plotMC) hmc2->Fit(fmc,"RNI");
  if (plotMC) hmc3->Fit(fmc,"RNI");
  fmc->SetLineColor(kGreen+3);

  //TF1 *fref = f50ns7; // compare to 50 ns data for AK7
  //TF1 *fref = f50ns4; // compare to 50 ns data for AK4 (previous default)
  //TF1 *fref = fmc; // compare to gen MC

  // Calculate ratio of fits
  TF1 *f1rs[nera];
  for (int iera = 0; iera != nera; ++iera) {

    //fref = f50ns4;

    TF1 *f1 = f1s[iera];
    TF1 *f1r = new TF1(Form("f1r_%d",iera),
		       Form("[0]*pow(x,[1]+[2]*log(x))"
			    "*pow(1-2.*x*cosh(%1.2f)/13000.,[3]) / "
			    "([4]*pow(x,[5]+[6]*log(x))"
			    "*pow(1-2.*x*cosh(%1.2f)/13000.,[7]))",
			    eta,eta),
		       ptmin,emax/cosh(eta));
    f1r->SetParameters(f1->GetParameter(0),f1->GetParameter(1),
		       f1->GetParameter(2),f1->GetParameter(3),
		       fref->GetParameter(0),fref->GetParameter(1),
		       fref->GetParameter(2),fref->GetParameter(3));
    //f1r->SetLineColor(color[iera]);
    f1r->SetLineColor(mcolor[eras[iera]]);
    f1rs[iera] = f1r;
  } // for iera

  // Calculate ratio of fit (nominally should be 1 when 50 ns is fref)
  TF1 *f50nsr = new TF1("f50nsr",
			Form("[0]*pow(x,[1]+[2]*log(x))"
		       "*pow(1-2.*x*cosh(%1.2f)/13000.,[3]) / "
			     "([4]*pow(x,[5]+[6]*log(x))"
			     "*pow(1-2.*x*cosh(%1.2f)/13000.,[7]))",
			     eta,eta),
			ptmin,emax/cosh(eta));
  f50nsr->SetParameters(f50ns4->GetParameter(0),f50ns4->GetParameter(1),
			f50ns4->GetParameter(2),f50ns4->GetParameter(3),
			fref->GetParameter(0),fref->GetParameter(1),
			fref->GetParameter(2),fref->GetParameter(3));
  f50nsr->SetLineColor(kGreen+2);



  TF1 *fmcr = new TF1("fmcr",
		      Form("[0]*pow(x,[1]+[2]*log(x))"
			   "*pow(1-2.*x*cosh(%1.2f)/13000.,[3]) / "
			   "([4]*pow(x,[5]+[6]*log(x))"
			   "*pow(1-2.*x*cosh(%1.2f)/13000.,[7]))",
			   eta,eta),
		      ptmin,emax/cosh(eta));
  fmcr->SetParameters(fmc->GetParameter(0),fmc->GetParameter(1),
		      fmc->GetParameter(2),fmc->GetParameter(3),
		      fref->GetParameter(0),fref->GetParameter(1),
		      fref->GetParameter(2),fref->GetParameter(3));
  fmc->SetLineColor(kGreen+3);

  // Divide new data by reference (was fit now Run2 data)
  //TH1F *hrs[nera];
  TH1D *hrs[nera];
  for (int iera = 0; iera != nera; ++iera) {

    //fref = f50ns4;
    //TH1F *hr = (TH1F*)h1s[iera]->Clone(Form("hr_%d",iera));
    TH1D *hr = (TH1D*)h1s[iera]->Clone(Form("hr_%d",iera));
    //hr->Divide(fref);
    if (eras[iera]!="16SC") {
      if (sdir=="17UL") {
	for (int i = 0; i != hr->GetNbinsX()+1; ++i) {
	  // c = a/(a+b) => (dc/c)^2 = (1-c)^2 * ( (da/a)^2 + (db/b)^2 )
	  // t = (a+b) => dt^2 = da^2 + db^2 => db^2 = dt^2 - da^2 
	  string se = eras[iera];
	  double ref = href->GetBinContent(i);
	  double eref = href->GetBinError(i);
	  double enom = hr->GetBinError(i);
	  double da = enom / ref;
	  double db = sqrt(fabs(eref*eref - enom*enom)) / ref;
	  double c = lumi[eras[iera]]/lumi["17ULBCDEF"];
	  //double dc = c*(1-c)*sqrt(da*da + db*db);
	  double dc = (1-c)*sqrt(da*da + db*db);
	  hr->SetBinContent(i, hr->GetBinContent(i) / ref);
	  hr->SetBinError(i, dc);
	}
      } // "UL17"
      else
	hr->Divide(href);
    }
    else { // Suman had different bin range
      for (int i = 1; i != hr->GetNbinsX()+1; ++i) {
	int j = href->FindBin(hr->GetBinCenter(i));
	double ref = href->GetBinContent(j);
	hr->SetBinContent(i, hr->GetBinContent(i) / ref);
	hr->SetBinError(i, hr->GetBinError(i) / ref);
      }
    }

    hrs[iera] = hr;
  } // for iera

  // Divide new data by ABCD for direct time stability
  /*
  TH1F *h1rs[nera];
  for (int iera = 1; iera != nera; ++iera) {
    TH1F *h1r = (TH1F*)h1s[iera]->Clone(Form("h1r_%d",iera));
    h1r->Divide(h1s[0]);
    h1rs[iera] = h1r;
  } // for iera
  */

  // Ratio of 50 ns data to reference (50 ns) fit
  //TH1F *hr50ns4 = (TH1F*)h50ns4->Clone("hr50ns4");
  TH1D *hr50ns4 = (plot2015 ? (TH1D*)h50ns4->Clone("hr50ns4") : 0);
  if (plot2015) hr50ns4->Divide(fref);//f50ns4);//fref);

  //TH1F *hr50ns7 = (TH1F*)h50ns7->Clone("hr50ns7");
  TH1D *hr50ns7 = (plot2015 ? (TH1D*)h50ns7->Clone("hr50ns7") : 0);
  if (plot2015) hr50ns7->Divide(fref);//f50ns7);//fref);

  TH1D *hrmc = (plotMC ? (TH1D*)hmc->Clone("hrmc") : 0);
  //hrmc->Divide(fref);
  if (plotMC) hrmc->Divide(href);

  TH1D *hrmc2 = (plotMC ? (TH1D*)hmc2->Clone("hrmc2") : 0);
  //hrmc2->Divide(fref);
  if (plotMC) hrmc2->Divide(href);

  TH1D *hrmc3 = (plotMC ? (TH1D*)hmc3->Clone("hrmc3") : 0);
  //hrmc3->Divide(fref);
  if (plotMC) hrmc3->Divide(href);
  
  TH1D *h = new TH1D("h",";Jet p_{T} (GeV);"
		     "d#sigma^{2} / dy dp_{T} (pb / GeV)",
		     int(ptmax-ptmin),ptmin,ptmax);
  h->SetMinimum(xsecmin);//1e-6 *1.0001);
  h->SetMaximum(xsecmax);//1e+8 *0.9999);

  //TH1D *h2 = new TH1D("h2",";Jet p_{T} (GeV);Ratio to 50 ns"
  //TH1D *h2 = new TH1D("h2",";Jet p_{T} (GeV);Ratio to 2017UL",
  //TH1D *h2 = new TH1D("h2",";Jet p_{T} (GeV);Ratio to 2017",
  TH1D *h2 = new TH1D("h2",";Jet p_{T} (GeV);Ratio to Run2",
		      int(ptmax-ptmin),ptmin,ptmax);
  if (plotVs17UL) h2->SetYTitle("Ratio to 2017UL");
  h2->SetMinimum(eta==2.5 ? 0.0 : 0.5);//0.);
  h2->SetMaximum(eta==2.5 ? 2.5 : 1.5);//fabs(eta)>=2.5 ? 5.0 : 2.0);
  //if (fabs(eta)>=3.0) h2->GetYaxis()->SetTitle("Ratio to 2016");
  h2->GetXaxis()->SetNoExponent();
  h2->GetXaxis()->SetMoreLogLabels();

  //lumi_13TeV = "Run2018ABCD 59.9 fb^{-1}";
  //lumi_13TeV = "Run2017BCDEF 41.4 fb^{-1}";
  //lumi_13TeV = "Run2016BCDEFGH 36.5 fb^{-1}";
  //lumi_13TeV = "Run2015 71 pb^{-1}";
  lumi_13TeV = "Run2 137.8 fb^{-1}";
  if (sdir=="2016") lumi_13TeV = "Run2016 X fb^{-1}";
  if (sdir=="2017") lumi_13TeV = "Run2017 Y fb^{-1}";
  if (sdir=="2018") lumi_13TeV = "Run2018 Z fb^{-1}";
  TCanvas *c1 = tdrDiCanvas("c1",h,h2,4,11);

  c1->cd(1);
  //tdrDraw(h50ns7,"Pz",kOpenSquare,kGreen+2,kSolid,kGreen+2);
  for (int iera = 0; iera != nera; ++iera) {
    //tdrDraw(h1s[iera],"Pz",marker[iera],color[iera],kSolid,color[iera]);
    string se = eras[iera];
    tdrDraw(h1s[iera],"Pz",mmarker[se],mcolor[se],kSolid,mcolor[se]);
    if (mmarker[se]==kFullDiamond || mmarker[se]==kOpenDiamond)
      h1s[iera]->SetMarkerSize(1.5);
  }
  //tdrDraw(h50ns4,"Pz",kOpenDiamond,kGreen+2,kSolid,kGreen+2);
  if (plotMC) {
    tdrDraw(hmc3,"Pz",kOpenDiamond,kRed+2,kSolid,kRed+2);
    hmc3->SetMarkerSize(0.7);
    tdrDraw(hmc2,"Pz",kOpenDiamond,kOrange+2,kSolid,kOrange+2);
    hmc2->SetMarkerSize(0.7);
    tdrDraw(hmc,"Pz",kOpenDiamond,kGreen+3,kSolid,kGreen+3);
    hmc->SetMarkerSize(0.7);
  }
  gPad->SetLogx();
  gPad->SetLogy();
  
  //TLegend *leg = tdrLeg(0.50,0.90-0.06*nera,0.80,0.90);
  //for (int iera = 0; iera != nera; ++iera) {
  //leg->AddEntry(h1s[iera],label[iera],"PL");
  //}
  const int nmax = 6;
  const int nx = 3;//2;//(eta<3.0 ? 2 : 1);
  TLegend *leg1 = tdrLeg(0.50,0.90-0.06*min(nmax-1,nera+nx),0.8,0.9);
  TLegend *leg2 = tdrLeg(0.70,0.90-0.06*max(2,min(nmax,nera+nx-nmax+1)),1.0,0.9);
  for (int iera = 0; iera != nera; ++iera) {
    //if (iera<nmax)  leg1->AddEntry(h1s[iera],label[iera],"PL");
    //if (iera>=nmax) leg2->AddEntry(h1s[iera],label[iera],"PL");
    string se = eras[iera];
    if (iera<nmax)  leg1->AddEntry(h1s[iera],label[se],"PL");
    if (iera>=nmax) leg2->AddEntry(h1s[iera],label[se],"PL");
  }
  if (fabs(eta)<3.0) {
    //if (nera<nmax)  leg1->AddEntry(h50ns4,"2015 pub.","PL");
    //if (nera>=nmax) leg2->AddEntry(h50ns4,"2015 pub.","PL");
  }
  if (plotMC) {
    if (nera<nmax-nx-1)  leg1->AddEntry(hmc,"MC","PL");
    if (nera>=nmax-nx-1) leg2->AddEntry(hmc,"MC","PL");
    if (nera<nmax-nx-1)  leg1->AddEntry(hmc2,"P8M2-10","PL");
    if (nera>=nmax-nx-1) leg2->AddEntry(hmc2,"P8M2-10","PL");
    if (nera<nmax-nx-1)  leg1->AddEntry(hmc3,"P8M2-5","PL");
    if (nera>=nmax-nx-1) leg2->AddEntry(hmc3,"P8M2-5","PL");
  }
  //leg->AddEntry(h50ns4,"74X 50 ns AK4","PL");
  //leg->AddEntry(h50ns7,"74X 50 ns AK7","PL");

  c1->cd(2);
  //tdrDraw(hr50ns7,"Pz",kOpenSquare,kGreen+2,kSolid,kGreen+2);
  for (int iera = 0; iera != nera; ++iera) {
    //tdrDraw(hrs[iera],"Pz",marker[iera],color[iera],kSolid,color[iera]);
    string se = eras[iera];
    tdrDraw(hrs[iera],"Pz",mmarker[se],mcolor[se],kSolid,mcolor[se]);
    if (mmarker[se]==kFullDiamond || mmarker[se]==kOpenDiamond)
      hrs[iera]->SetMarkerSize(1.5);
    hrs[iera]->GetXaxis()->SetRangeUser(ptmin,ptmax);
  }
  //tdrDraw(hr50ns4,"Pz",kOpenDiamond,kGreen+2,kSolid,kGreen+2);
  if (plotMC) {
    tdrDraw(hrmc3,"Pz",kOpenDiamond,kRed+2,kSolid,kRed+2);
    hrmc3->SetMarkerSize(0.7);
    tdrDraw(hrmc2,"Pz",kOpenDiamond,kOrange+2,kSolid,kOrange+2);
    hrmc2->SetMarkerSize(0.7);
    tdrDraw(hrmc,"Pz",kOpenDiamond,kGreen+3,kSolid,kGreen+3);
    hrmc->SetMarkerSize(0.7);
  }
  gPad->SetLogx();

  c1->cd(1);
  f50ns4->SetRange(ptmin,0.9*emax/cosh(eta));//3500);
  //f50ns4->DrawClone("SAME");
  //f50ns7->SetRange(ptmin,3500);
  //f50ns7->Draw("SAME");
  // Draw new fits as well (or don't)
  for (int iera = 0; iera != nera; ++iera) {
    //f1s[iera]->Draw("SAME");
  }

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  tex->SetTextColor(kGreen+2);
  //tex->DrawLatex(0.20,0.10,Form("fit(50ns) AK7: #chi^{2}/NDF = %1.1f/%d",
  //			f50ns7->GetChisquare(),f50ns7->GetNDF()));
  //if (fabs(eta)>=3.0) 
  //tex->DrawLatex(0.20,0.05,Form("fit(2016) AK4: #chi^{2}/NDF = %1.1f/%d",
  //			  f50ns4->GetChisquare(),f50ns4->GetNDF()));
  //else
  //tex->DrawLatex(0.20,0.05,Form("fit(50ns) AK4: #chi^{2}/NDF = %1.1f/%d",
  //			  f50ns4->GetChisquare(),f50ns4->GetNDF()));
  //tex->DrawLatex(0.20,0.05,Form("fit(MC) AK4: #chi^{2}/NDF = %1.1f/%d",
  //			fmc->GetChisquare(),fmc->GetNDF()));
						
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.20,0.20,Form("y #in %s",cy)); 

  c1->cd(2);
  f50nsr->Draw("SAME");
  //fmcr->Draw("SAME");
  // Draw ratio of fits as well (or don't)
  for (int iera = 0; iera != nera; ++iera) {
    //f1rs[iera]->Draw("SAME");
  }

  c1->SaveAs(Form("pdf/%s/drawDeltaJEC_jetPt_%s.pdf",cdir,cy2));


  //TH1D *h3 = new TH1D("h3",";Jet p_{T} (GeV);#DeltaJEC-equivalent vs 17UL (%)",
  //TH1D *h3 = new TH1D("h3",";Jet p_{T} (GeV);#DeltaJEC-equivalent vs 17 (%)",
  TH1D *h3 = new TH1D("h3",";Jet p_{T} (GeV);#DeltaJEC-equivalent vs Run2 (%)",
		      int(ptmax-ptmin),ptmin,ptmax);
  if (plotVs17UL) h3->SetYTitle("#DeltaJEC-equivalent vs 17UL (%)");
  h3->SetMinimum(ymin);
  h3->SetMaximum(ymax);
  h3->GetXaxis()->SetMoreLogLabels();
  h3->GetXaxis()->SetNoExponent();
  
  // Calculate reference uncertainty
  /*
  TH1F *hunc18 = (TH1F*)h50ns4->Clone("hunc18");
  TH1F *hunc17 = (TH1F*)h50ns4->Clone("hunc17");
  TH1F *hunc16 = (TH1F*)h50ns4->Clone("hunc16");
  TH1F *huncs1 = (TH1F*)h50ns4->Clone("huncs1");
  TH1F *huncs2 = (TH1F*)h50ns4->Clone("huncs2");
  TH1F *hunc15 = (TH1F*)h50ns4->Clone("hunc15");
  TH1F *hunc12 = (TH1F*)h50ns4->Clone("hunc12");
  */
  TH1D *hunc18 = (TH1D*)href->Clone("hunc18");
  TH1D *hunc17 = (TH1D*)href->Clone("hunc17");
  TH1D *hunc16 = (TH1D*)href->Clone("hunc16");
  TH1D *huncs1 = (TH1D*)href->Clone("huncs1");
  TH1D *huncs2 = (TH1D*)href->Clone("huncs2");
  TH1D *hunc15 = (TH1D*)href->Clone("hunc15");
  TH1D *hunc12 = (TH1D*)href->Clone("hunc12");
  TH1D *huncCS = (TH1D*)href->Clone("huncCS");
  for (int i = 1; i != hunc18->GetNbinsX()+1; ++i) {
    
    double pt = hunc18->GetBinCenter(i);

    unc18->setJetEta(eta+0.25);
    unc18->setJetPt(pt);
    hunc18->SetBinContent(i, 0);
    hunc18->SetBinError(i, 100.*unc18->getUncertainty(true));

    unc17->setJetEta(eta+0.25);
    unc17->setJetPt(pt);
    hunc17->SetBinContent(i, 0);
    hunc17->SetBinError(i, 100.*unc17->getUncertainty(true));

    unc16->setJetEta(eta+0.25);
    unc16->setJetPt(pt);
    hunc16->SetBinContent(i, 0);
    hunc16->SetBinError(i, 100.*unc16->getUncertainty(true));

    unc15->setJetEta(eta+0.25);
    unc15->setJetPt(pt);
    hunc15->SetBinContent(i, 0);
    hunc15->SetBinError(i, 100.*unc15->getUncertainty(true));

    unc12->setJetEta(eta+0.25);
    unc12->setJetPt(pt);
    hunc12->SetBinContent(i, 0);
    hunc12->SetBinError(i, 100.*unc12->getUncertainty(true));

    //double kC = 1.01;
    double kC = 1.01+0.23/40.-0.23/pt;
    jecC->setJetPt(pt);
    jecC->setJetEta(eta);
    //double cC = jecC->getCorrection();
    double cC = getJEC(jecC, eta, pt, _rhoDE);
    jecS->setJetPt(pt);
    jecS->setJetEta(eta);
    //double cS = jecS->getCorrection();
    double cS = getJEC(jecS, eta, pt, _rhoDE);
    huncCS->SetBinContent(i, 0);
    huncCS->SetBinError(i, 0.5*100.*fabs(kC*cC-cS));

    uncs1->setJetEta(eta+0.25);
    uncs1->setJetPt(pt);
    huncs1->SetBinContent(i, 100.*uncs1->getUncertainty(true));
    huncs1->SetBinError(i, 0.);

    uncs2a->setJetEta(eta+0.25);
    uncs2a->setJetPt(pt);
    uncs2b->setJetEta(eta+0.25);
    uncs2b->setJetPt(pt);
    huncs2->SetBinContent(i, 100.*2.*uncs2a->getUncertainty(true) +
			  100.*2.*uncs2b->getUncertainty(true));
    huncs2->SetBinError(i, 0.);
    
  }

  // Use Run2 histogram 'href' as reference fit 'fref'
  _href = href;
  TF1 *fref2 = new TF1("fref2",fRef,10.,6500./cosh(eta),1);
  fref2->SetParameter(0, eta);

  // Calculate JEC equivalent shift
  //TH1F *hds[nera];
  //TH1F *huncl = (TH1F*)h50ns4->Clone("huncl");
  TH1D *hds[nera];
  TH1D *huncl = (TH1D*)href->Clone("huncl");
  for (int iera = 0; iera != nera; ++iera) {
    
    //fref = f50ns4;
    //TH1F *hr = hrs[iera];
    //TH1F *hd = (TH1F*)hr->Clone(Form("hd_%d",iera));
    TH1D *hr = hrs[iera];
    TH1D *hd = (TH1D*)hr->Clone(Form("hd_%d",iera));
    for (int i = 1; i != hr->GetNbinsX()+1; ++i) {

      if (hr->GetBinContent(i)==0) continue;
      //double x = hr->GetBinCenter(i);
      double x = hr->GetBinLowEdge(i);
      double xmin = max(10.,0.5*x);
      double xmax = min(1.5*x, 6500./cosh(eta)-1.);
      if (x<15) continue;
      
      double y = fref2->Eval(x);
      double y2 = y * hr->GetBinContent(i);
      //double x2 = fref2->GetX(y2,0.85*x,1.15*x,1e-3);
      //double x2 = fref2->GetX(y2,0.85*x,1.15*x);
      //double x2 = fref2->GetX(y2,0.5*x,xmax,_prec);//1e-3);
      double x2 = fref2->GetX(y2,xmin,xmax,_prec);//1e-3);
      hd->SetBinContent(i, 100.*(x/x2-1));

      // Luminosity uncertainty from ABCD
      if (iera==0) {
	double y3 = y / 1.027; // lumi uncertainty for 2016, after 3.3% patch
	//double x3 = fref2->GetX(y3,0.85*x,1.15*x,1e-3);
	//double x3 = fref2->GetX(y3,0.5*x,xmax,_prec);//1e-3);
	double x3 = fref2->GetX(y3,xmin,xmax,_prec);//1e-3);
	int j = huncl->FindBin(hd->GetBinCenter(i));
	huncl->SetBinContent(j, 0);
	huncl->SetBinError(j, 100.*fabs(x/x3-1));
      }

      double y2_up = y * (hr->GetBinContent(i) + hr->GetBinError(i));
      //double x2_dw = fref2->GetX(y2_up,0.85*x,1.15*x);
      //double x2_dw = fref2->GetX(y2_up,0.5*x,xmax,_prec);//1e-3);
      double x2_dw = fref2->GetX(y2_up,xmin,xmax,_prec);//1e-3);
      double y2_dw = y * (hr->GetBinContent(i) - hr->GetBinError(i));
      //double x2_up = fref2->GetX(y2_dw,0.85*x,1.15*x);
      //double x2_up = fref2->GetX(y2_dw,0.5*x,xmax,_prec);//x1e-3);
      double x2_up = fref2->GetX(y2_dw,xmin,xmax,_prec);//x1e-3);
      hd->SetBinError(i, 100.*sqrt(pow((x2_up-x2_dw)/x2,2) + pow(0.001,2)));
    } // for i
    hds[iera] = hd;

  } // for iera


  // Calculate cross section uncertainty
  TH1D *huncx = (TH1D*)href->Clone("huncx");
  if (true) {
    
    huncx->Reset();
    for (int i = 1; i != huncx->GetNbinsX()+1; ++i) {
      if (href->GetBinContent(i)!=0) {
	double pt = huncx->GetBinCenter(i);
	//if (pt<10 || pt*cosh(eta)>6500*0.93) continue;
	if (pt<10 || pt*cosh(eta)>6500*0.90) continue;
	double ejes = 0.01*hunc12->GetBinError(hunc12->FindBin(pt));
	double y = fref2->Eval(pt);
	double yup = fref2->Eval(max(pt*(1-ejes),5.));
	double ydw = fref2->Eval(min(pt*(1+ejes),0.99*6500./cosh(eta)));
	double ymid = 0.5*(yup+ydw);
	double ey = 0.5*fabs(yup-ydw);
	huncx->SetBinContent(i, ymid/ymid);
	huncx->SetBinError(i, ey/ymid);
      }
    } // for i

  } // huncx


  

  /////////////////////
  // Data/Run2 plots //
  /////////////////////

  //TH1D *h1b = new TH1D("h1b",";Jet p_{T} (GeV);Ratio to 2017UL",
  //TH1D *h1b = new TH1D("h1b",";Jet p_{T} (GeV);Ratio to 2017",
  TH1D *h1b = new TH1D("h1b",";Jet p_{T} (GeV);Ratio to Run2",
		       int(ptmax-ptmin),ptmin,ptmax);
  if (plotVs17UL) h1b->SetYTitle("Ratio to 2017UL");
  h1b->SetMinimum(eta==2.5 ? 0.0 : 0.5);
  h1b->SetMaximum(eta==2.5 ? 2.5 : 1.5);
  h1b->GetXaxis()->SetNoExponent();
  h1b->GetXaxis()->SetMoreLogLabels();

  TCanvas *c1b = tdrCanvas("c1b",h1b,4,11,kSquare);

  tdrDraw(huncx,"E3", kNone, kBlack,kSolid,-1,1001,kRed-9); // 2012
  huncx->SetFillColorAlpha(kRed-9,0.70);  

  for (int iera = 0; iera != nera; ++iera) {
    string se = eras[iera];
    tdrDraw(hrs[iera],"Pz",mmarker[se],mcolor[se],kSolid,mcolor[se]);
    if (mmarker[se]==kFullDiamond || mmarker[se]==kOpenDiamond)
      hrs[iera]->SetMarkerSize(1.0);//1.5);
    else
      hrs[iera]->SetMarkerSize(0.7);

    hrs[iera]->GetXaxis()->SetRangeUser(ptmin,ptmax);
    
    // PATCH UL17
    if (sdir=="17UL" && se=="17UL" && getDetData && unfoldData && false) {
      TFile *fout = new TFile("pdf/drawDeltaJEC_UL17patch.root","RECREATE");
      hrs[iera]->Write("jetRatio_BCDEF",TObject::kOverwrite);
      fout->Close();
    }
  }
  if (plotMC) {
    tdrDraw(hrmc3,"Pz",kOpenDiamond,kRed+2,kSolid,kRed+2);
    hrmc3->SetMarkerSize(0.7);
    tdrDraw(hrmc2,"Pz",kOpenDiamond,kOrange+2,kSolid,kOrange+2);
    hrmc2->SetMarkerSize(0.7);
    tdrDraw(hrmc,"Pz",kOpenDiamond,kGreen+3,kSolid,kGreen+3);
    hrmc->SetMarkerSize(0.7);
  }
  gPad->SetLogx();

  TLegend *leg1b = tdrLeg(0.40,0.90-0.05*min(nmax-1,nera+nx),0.70,0.9);
  TLegend *leg2b = tdrLeg(0.60,0.90-0.05*max(2,min(nmax,nera+nx-nmax+1)),0.90,0.9);
  for (int iera = 0; iera != nera; ++iera) {
    string se = eras[iera];
    if (iera<nmax)  leg1b->AddEntry(h1s[iera],label[se],"PL");
    if (iera>=nmax) leg2b->AddEntry(h1s[iera],label[se],"PL");
  }
  if (plotMC) {
    if (nera<nmax-nx-1)  leg1b->AddEntry(hmc,"MC","PL");
    if (nera>=nmax-nx-1) leg2b->AddEntry(hmc,"MC","PL");
    if (nera<nmax-nx-1)  leg1b->AddEntry(hmc2,"P8M2-10","PL");
    if (nera>=nmax-nx-1) leg2b->AddEntry(hmc2,"P8M2-10","PL");
    if (nera<nmax-nx-1)  leg1b->AddEntry(hmc3,"P8M2-5","PL");
    if (nera>=nmax-nx-1) leg2b->AddEntry(hmc3,"P8M2-5","PL");
  }

  TLegend *leg3b = tdrLeg(0.57,0.15,0.87,0.19);
  leg3b->AddEntry(huncx,"2012 Win14_V8","F");

  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.20,0.20,Form("|y| #in %s",cy)); 

  c1b->SaveAs(Form("pdf/%s/drawDeltaJEC_xsecRatio_%s.pdf",cdir,cy2));


  ////////////////////
  // DeltaJEC plots //
  ////////////////////

  TCanvas *c2 = tdrCanvas("c2",h3,4,11,kSquare);

  //tdrDraw(hunc15,"E3", kNone, kBlack,kSolid,-1,1001,kOrange-9); // 2015
  //
  if (sdir=="" || sdir=="2018")
    tdrDraw(hunc18,"E3", kSolid, kBlack,kSolid,-1,1001,kYellow+1); // 2018
  hunc18->SetFillColorAlpha(kYellow+1,0.70);
  if (sdir=="" || sdir=="2017" || sdir=="2017UL" || sdir=="17UL")
    tdrDraw(hunc17,"E3", kSolid, kBlack,kSolid,-1,1001,kCyan-6); // 2017
  hunc17->SetFillColorAlpha(kCyan-6,0.70);
  //tdrDraw(hunc16,"E3", kSolid, kBlack,kSolid,-1,1001,kGreen-6); // 2016
  //hunc16->SetFillColorAlpha(kGreen-6,0.70);
  if (sdir=="" || sdir=="2016")
    tdrDraw(hunc16,"E3", kSolid, kBlack,kSolid,-1,1001,kViolet-8); // 2016
  hunc16->SetFillColorAlpha(kViolet-8,0.70);

  if (sdir=="2017UL" && false) {
    tdrDraw(huncCS,"E3", kSolid, kBlack,kSolid,-1,1001,kCyan-6); // 2017UL
    huncCS->SetFillColorAlpha(kCyan-6,0.70);
    TF1 *fCS0 = new TF1("fCS0","100.*[0]/x",15,3500);
    //fCS0->SetParameter(0,1.72396*0.138); // p3 fit uncertainty for BCDEF
    fCS0->SetParameter(0,1.72396*0.138); // p3 fit uncertainty for BCDEF
    fCS0->SetLineWidth(3);
    fCS0->SetLineColor(kBlack);
    fCS0->Draw("SAME");
    TF1 *fCS1 = new TF1("fCS1","100.*[0]/x",15,3500);
    fCS1->SetParameter(0,0.025*15.);
    fCS1->SetLineWidth(2);
    fCS1->Draw("SAME");
    TF1 *fCS2 = new TF1("fCS2","100.*[0]/(x*log(x))",15,3500);
    fCS2->SetParameter(0,0.025*15.*log(15.));
    fCS2->SetLineColor(kBlue);
    fCS2->Draw("SAME");
  }

  //tdrDraw(hunc15,"E3", kNone, kBlack,kSolid,-1,1001,kOrange-9); // 2015
  //hunc15->SetFillColorAlpha(kOrange-9,0.70);
  tdrDraw(hunc12,"E3", kNone, kBlack,kSolid,-1,1001,kRed-9); // 2012
  hunc12->SetFillColorAlpha(kRed-9,0.70);
  //tdrDraw(huncl,"E3", kNone, kBlack,kSolid,-1,1001,kBlue-9);
  //tdrDraw(hunc,"E3", kSolid, kBlack,kSolid,-1,3001,kYellow+1);
  //tdrDraw(huncs1,"HIST][",kSolid,kBlack,kSolid,-1,kNone,kOrange-9); //RelPtBal
  //tdrDraw(huncs2,"HIST][",kSolid,kBlue,kSolid,-1,kNone,kOrange-9); //ReltPtEC1+2
  
  if (sdir=="" || sdir=="2018")
    tdrDraw((TH1D*)hunc18->Clone("hunc18b"),"E3",kNone,0,kSolid,kYellow+2,kNone,kYellow+2);
  if (sdir=="" || sdir=="2017" || sdir=="2017UL" ||  sdir=="17UL")
    tdrDraw((TH1D*)hunc17->Clone("hunc17b"),"E3",kNone,0,kSolid,kCyan+2,kNone,kCyan+2);
  if (sdir=="" || sdir=="2016")
    tdrDraw((TH1D*)hunc16->Clone("hunc16b"),"E3",kNone,0,kSolid,kViolet+2,kNone,kViolet+2);

  for  (int iera = 0; iera != nera; ++iera) {
    //tdrDraw(hds[iera],"Pz",marker[iera],color[iera],kSolid,color[iera]);
    string se = eras[iera];
    tdrDraw(hds[iera],"Pz",mmarker[se],mcolor[se],kSolid,mcolor[se]);
    if (mmarker[se]==kFullDiamond || mmarker[se]==kOpenDiamond)
      //hds[iera]->SetMarkerSize(1.5);
      hds[iera]->SetMarkerSize(1.0);
    else 
      hds[iera]->SetMarkerSize(0.7);
    hds[iera]->GetXaxis()->SetRangeUser(ptmin,ptmax);

  }
  gPad->SetLogx();

  //TLegend *legd = tdrLeg(0.50,0.90-0.04*(nera+3),0.80,0.90);
  //for  (int iera = 0; iera != nera; ++iera) {
  //legd->AddEntry(hrs[iera],label[iera],"PL");
  //}
  const int nmaxd = 6;
  TLegend *legd1 = tdrLeg(0.40,0.9-0.04*min(nmaxd,nera),0.70,0.9);
  TLegend *legd2 = tdrLeg(0.60,0.9-0.04*max(0,min(nmaxd,nera-nmaxd)),0.90,0.9);
  for  (int iera = 0; iera != nera; ++iera) {
    //if (iera<nmaxd)  legd1->AddEntry(hrs[iera],label[iera],"PL");
    //if (iera>=nmaxd) legd2->AddEntry(hrs[iera],label[iera],"PL");
    string se = eras[iera];
    if (iera<nmaxd)  legd1->AddEntry(hrs[iera],label[se],"PL");
    if (iera>=nmaxd) legd2->AddEntry(hrs[iera],label[se],"PL");
  }
  
  //TLegend *legu = tdrLeg(0.50,0.15,0.80,0.23);
  const int nmaxu1 = 3;
  const int nmaxu2 = 1;//3;
  TLegend *legu2 = tdrLeg(0.57,0.15,0.87,0.15+0.04*nmaxu2);
  TLegend *legu1 = tdrLeg(0.17,0.15,0.47,0.15+0.04*nmaxu1);
  //legu1->AddEntry(hunc18,"2018 Aut18_V10+unc","F");
  //legu1->AddEntry(hunc18,"Aut18_V10+V13tmp","F");
  //legu1->AddEntry(hunc18,"Aut18_V10+V14unc","F");
  if (sdir=="" || sdir=="2018")
    //legu1->AddEntry(hunc18,"2018 Aut18_V18","F");
    legu1->AddEntry(hunc18,"2018 Aut18_V19","F");
    //legu1->AddEntry(hunc18,"Aut18_V10+V18unc","F");
  //legu1->AddEntry(hunc18,"Aut18_V13h+V14unc","F");
  //legu1->AddEntry(hunc18,"2018 Aut18_V8","F");
  //legu1->AddEntry(hunc17,"2017 17Nov_V11","F");
  if (sdir=="" || sdir=="2017" || sdir=="2017UL" || sdir=="17UL")
    legu1->AddEntry(hunc17,"2017 17Nov_V32","F");
  //legu1->AddEntry(hunc16,"2016 07Aug_V17","F");
  if (sdir=="" || sdir=="2016")
    legu1->AddEntry(hunc16,"2016 07Aug_V11","F");
  //legu2->AddEntry(hunc15,"2015 Sum16_50nsV5","F");
  legu2->AddEntry(hunc12,"2012 Win14_V8","F");
  //legd->AddEntry(huncl,"2.6% lum.","F");
  
  /*
  tex->SetTextColor(kBlack); tex->SetTextSize(0.040);
  tex->DrawLatex(0.19,0.75,Form("Anti-k_{T} R=0.4"));
  tex->DrawLatex(0.19,0.71,Form("PF+CHS unf."));
  tex->DrawLatex(0.19,0.67,Form("|y|#in %s",cy));
  */
  tex->SetTextColor(kBlack); tex->SetTextSize(0.045);
  tex->DrawLatex(0.65,0.88,Form("|y|#in %s",cy));
  tex->DrawLatex(0.65,0.83,Form("Anti-k_{T} R=0.4"));
  tex->DrawLatex(0.65,0.78,Form("PF+CHS unf."));
  tex->SetTextSize(0.025);
  //tex->DrawLatex(0.19,0.64,Form("vs HEPData-ins1459051-v1"));
  tex->SetTextSize(0.040);
  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(ptmin,0,ptmax,0);
  l->SetLineStyle(kDotted);
  l->DrawLine(ptmin,+1.5,ptmax,+1.5);
  l->DrawLine(ptmin,-1.5,ptmax,-1.5);

  // Save UL17 ratio for global fit input
  if (sdir=="17UL" && sy=="0.0-1.3") {
    TFile *fout = new TFile("pdf/drawDeltaJEC_17UL.root",
			    getDetData && !unfoldData ? "RECREATE" : "UPDATE");
    fout->cd();
    for (int iera = 0; iera != nera; ++iera) {
      string se = eras[iera];
      TH1D *h = (TH1D*)hds[iera]->Clone(Form("jet_Run%s",se.c_str()));
      for (int i = 1; i != h->GetNbinsX()+1; ++i) {
	h->SetBinContent(i, 1 + 0.01*h->GetBinContent(i));
	h->SetBinError(i, 0.01*h->GetBinError(i));
      }
      h->Write(Form("jet_Run%s_%s",se.c_str(),
		    getDetData ? (unfoldData ? "fwd3" : "det") :
		    (getFwd ? "fwd" : "dag")), TObject::kOverwrite);
    } // for iera
    fout->Close();
    curdir->cd();    
  } // sdir

  c2->RedrawAxis();
  c2->SaveAs(Form("pdf/%s/drawDeltaJEC_%s.pdf",cdir,cy2));
}


// Ansatz Kernel
//const double emax = 6500.;
int cnt_a = 0;
const int nk = 3; // number of kernel parameters (excluding pt, eta)
bool _dojes = true;
Double_t smearedAnsatzKernel(Double_t *x, Double_t *p) {

  if (++cnt_a%1000000==0) {
    cout << "+" << flush;
  }

  const double pt = x[0]; // true pT
  const double ptmeas = p[0]; // measured pT
  const double eta = p[4]; // rapidity

  double res = ptresolution(pt, eta+1e-3) * pt;
  double jes = (_dojes ? ptresponse(pt, eta+1e-3) : 1);
  const double s = TMath::Gaus(ptmeas, jes*pt, res, kTRUE);
  const double f = p[1] * pow(pt, p[2])
    * pow(1 - pt*cosh(eta) / emax, p[3]);

  return (f * s);
}

// Smeared Ansatzz
double _epsilon = 1e-12;
TF1 *_kernel = 0; // global variable, not pretty but works
Double_t smearedAnsatz(Double_t *x, Double_t *p) {

  const double pt = x[0];
  const double eta = p[3];

  if (!_kernel) _kernel = new TF1("_kernel", smearedAnsatzKernel,
				  1., emax/cosh(eta), nk+2);

  double res = ptresolution(pt, eta+1e-3) * pt;
  const double sigma = max(0.10, min(res/pt, 0.30));
  double ptmin = pt / (1. + 4.*sigma); // xmin*(1+4*sigma)=x
  ptmin = max(1.,ptmin); // safety check
  double ptmax = pt / (1. - 3.*sigma); // xmax*(1-3*sigma)=x
  //cout << Form("1pt %10.5f sigma %10.5f ptmin %10.5f ptmax %10.5f eta %10.5f",pt, sigma, ptmin, ptmax, eta) << endl << flush;
  ptmax = min(emax/cosh(eta), ptmax); // safety check
  //cout << Form("2pt %10.5f sigma %10.5f ptmin %10.5f ptmax %10.5f eta %10.5f",pt, sigma, ptmin, ptmax, eta) << endl << flush;

  const double par[nk+2] = {pt, p[0], p[1], p[2], p[3]};
  _kernel->SetParameters(&par[0]);

  // Set pT bin limits needed in smearing matrix generation
  //if (p[5]>0 && p[5]<emax/cosh(eta)) ptmin = p[5];
  //if (p[6]>0 && p[6]<emax/cosh(eta)) ptmax = p[6];

  return ( _kernel->Integral(ptmin, ptmax, _epsilon) );
}

// Tool to estimate unfolding corrections
void unfold(string sy = "0.0-1.3", int ieta = 0) {


  //TFile *f = new TFile("rootfiles/common2018_V7.root","READ");
  //TFile *f = new TFile("rootfiles/common2018_V10.root","READ");
  //TFile *f = new TFile("rootfiles/common2018_V10_hotzone-2.root","READ");
  //TFile *f = new TFile("rootfiles/common2018_V13h.root","READ");
  //TFile *f = new TFile("rootfiles/common2016_LegacyIOVs_v3.root","READ");
  //TFile *f = new TFile("rootfiles/common2016_October2018_V17.root","READ");
  //TFile *f = new TFile("rootfiles/common2018_V19-5.root","READ");
  //TFile *f = new TFile("rootfiles/common2018_V19-5.root","READ");
  TFile *f = new TFile("rootfiles/commonUL2017_V4_V2M4res_hotzone_scaled2p.root","READ");
  assert(f && !f->IsZombie());
  TFile *fout = new TFile("rootfiles/unfold.root",
			  ieta==0 ? "RECREATE" : "UPDATE");

  const char *cy = sy.c_str();
  //TH1D *hd = (TH1D*)f->Get(Form("ak4/y_%s/hptData_full2016_detector_%dbin",
  //cy,ieta));
  //TH1D *hd = (TH1D*)f->Get(Form("ak4/Eta_%s/hpt_data_2018_D_det",cy));
  TH1D *hd = (TH1D*)f->Get(Form("ak4/Eta_%s/hpt_data_2017_all_det",cy));

  assert(hd);
  //TH1D *hu = (TH1D*)f->Get(Form("ak4/y_%s/hptData_full2016_particle_%dbin",
  //				cy,ieta));
  //TH1D *hu = (TH1D*)f->Get(Form("ak4/Eta_%s/hpt_data_2018_D_det",cy));
  //TH1D *hu = (TH1D*)f->Get(Form("ak4/Eta_%s/hpt_data_2018_all_det",cy));
  TH1D *hu = (TH1D*)f->Get(Form("ak4/Eta_%s/hpt_data_2017_all_det",cy));
  assert(hu);

  /*
  //TF1 *f1 = new TF1(Form("fr_%d",ieta),"[0]+[1]*log(x)+[2]*log(x)*log(x)",
  TF1 *f1 = new TF1(Form("fr_%s",cy),"[0]+[1]*log(x)+[2]*log(x)*log(x)",
		    114,3000);
  f1->SetParameters(1,0.01,-0.001);
  hr->Fit(f1,"QRN");
  f1->Write();
  */

  // Values of NLO fits done by hand
  const int neta = 7;
  const double p[neta][nk+1] =
    {{1,-4.8332,12.101,0.0},
     {1,-4.8219,11.4926,0.5},
     {1,-4.84815,9.46945,1.0},
     {1,-4.77036,8.65309,1.5},
     {1,-4.6652,8.09749,2},
     {1,-4.97904,6.4373,2.5},
     {1,-5.72102,4.34217,3.2}};
  int i = max(0,ieta-1);
  double eta = p[i][nk];
  double maxpt = emax/cosh(eta);
  // For |eta|<1.3
  double maxpt2 = emax/cosh(0.5);
  double maxpt3 = emax/cosh(1.0);

  _ismcjer = false;
  _usejme = false;
  _jer_iov = run2016;
  
  // Initial fit of the NLO curve to a histogram
  TF1 *fus = new TF1(Form("fus%s",cy),
		     "[0]*pow(x,[1])*pow(1-x*cosh([3])/6500.,[2])",5,maxpt);
  fus->SetParameters(p[i][0], p[i][1], p[i][2], p[i][3]);

  // For |eta|<1.3
  TF1 *fus2 = new TF1(Form("fus2%s",cy),
		      "[0]*pow(x,[1])*pow(1-x*cosh([3])/6500.,[2])",5,maxpt2);
  fus2->SetParameters(p[1][0], p[1][1], p[1][2], p[1][3]);
  TF1 *fus3 = new TF1(Form("fus3%s",cy),
		      "[0]*pow(x,[1])*pow(1-x*cosh([3])/6500.,[2])",5,maxpt3);
  fus3->SetParameters(p[2][0], p[2][1], p[2][2], p[2][3]);
  
  // Smeared spectrum
  TF1 *fs = new TF1(Form("fs%s",cy),smearedAnsatz,5.,maxpt,nk+1);
  fs->SetParameters(fus->GetParameter(0), fus->GetParameter(1),
		    fus->GetParameter(2), fus->GetParameter(3));
  // For |eta|<1.3
  TF1 *fs2 = new TF1(Form("fs2%s",cy),smearedAnsatz,5.,maxpt2,nk+1);
  fs2->SetParameters(fus2->GetParameter(0), fus2->GetParameter(1),
		     fus2->GetParameter(2), fus2->GetParameter(3));
  TF1 *fs3 = new TF1(Form("fs3%s",cy),smearedAnsatz,5.,maxpt3,nk+1);
  fs3->SetParameters(fus3->GetParameter(0), fus3->GetParameter(1),
		     fus3->GetParameter(2), fus3->GetParameter(3));

  fout->cd();

  const int niov = 12;//6;//5;
  jer_iov iovs[niov] = {run1, run2016, run2017, run2018, run2018abc, run2018d,
			ul17, ul17b, ul17c, ul17d, ul17e, ul17f};
  for (int iov = 0; iov != niov; ++iov) {

    _jer_iov = iovs[iov]; // use by ptresolution in fs and fus
    // For UL17, scale also low pT JES time dependence
    if (_jer_iov==ul17 || _jer_iov==ul17b || _jer_iov==ul17c || 
	_jer_iov==ul17d || _jer_iov==ul17e || _jer_iov==ul17f) _dojes = true;
    else _dojes = false;

    TH1D *hr = (TH1D*)hu->Clone(Form("hr_%s_%d",cy,2015+iov));
    for (int i = 1; i != hr->GetNbinsX()+1; ++i) {
      int i1 = hd->FindBin(hr->GetBinLowEdge(i));
      int i2 = hd->FindBin(hr->GetBinLowEdge(i+1)-0.5);
      double yd = (hd->GetBinContent(i1)*hd->GetBinWidth(i1) +
		   hd->GetBinContent(i2)*hd->GetBinWidth(i2)) /
	(hd->GetBinWidth(i1) + hd->GetBinWidth(i2));
      double yu = hu->GetBinContent(i);
      double pt = hr->GetBinCenter(i);
      if (yu!=0 && yd!=0 && pt*cosh(eta)<emax) {
	
	double nd = fs->Eval(pt);
	double nu = fus->Eval(pt);
	hr->SetBinContent(i, nd/nu);//yd/yu);
	hr->SetBinError(i, 0);//hu->GetBinError(i)/hu->GetBinContent(i)
	//* hr->GetBinContent(i));
	if (ieta==0) {
	  double nd2 = (pt*cosh(0.5)<emax ? fs2->Eval(pt) : 0);
	  double nu2 = (pt*cosh(0.5)<emax ? fus2->Eval(pt) : 0);
	  double nd3 = (pt*cosh(1.0)<emax ? fs3->Eval(pt) : 0);
	  double nu3 = (pt*cosh(1.0)<emax ? fus3->Eval(pt) : 0);
	  hr->SetBinContent(i, (nd+nd2+0.6*nd3)/(nu+nu2+0.6*nd3));
	  hr->SetBinError(i, 0);//hu->GetBinError(i)/hu->GetBinContent(i)
	}
      }
      else {
	hr->SetBinContent(i, 0);
	hr->SetBinError(i, 0);
      }
    }
    //hr->Write(hr->GetName(),TObject::kOverwrite);
  } // for iov

  fout->Write();
  fout->Close();

} // unfold
