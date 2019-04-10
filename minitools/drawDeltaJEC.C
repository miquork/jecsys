
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

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "ptresolution.h"
#include "../tdrstyle_mod15.C"

// "Rebin(2)"
const int nptb = 51+10+6;
const float ptbins[nptb+1] =
  {15, 18, 21, 24, 28, 32, 37, 43, 49, 56,
   64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362,
   395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032,
   1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116,
   2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832,
   4037, 4252, 4477, 4713, 4961, 5220};

TH1F *rebinXsec(TH1F *h) {

  TH1F *h2 = new TH1F(Form("%s_rb2",h->GetName()),"",nptb,ptbins);
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


const bool correctECALprefire = true;
const double ptmin(15);//114;//64;
const double ptmax(3500.);
const double emax(6500.);
const double xsecmin(1.0001e-7);
const double xsecmax(0.9999e11);
void drawDeltaJEC(string sy = "0.0-0.5") {

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
  TFile *fin1 = new TFile("rootfiles/common2018_V7.root","READ");
  assert(fin1 && !fin1->IsZombie());
  
  //TFile *fin2 = new TFile("rootfiles/common2016_LegacyIOVs_v3.root","READ");
  //TFile *fin2 = new TFile("rootfiles/common2017_V11.root","READ");
  TFile *fin2 = new TFile("rootfiles/common2017_V32.root","READ");
  assert(fin2 && !fin2->IsZombie());

  //TFile *fin3 = new TFile("rootfiles/common2016_October2018_V17.root","READ");
  TFile *fin3 = new TFile("rootfiles/common2016_V11.root","READ");
  assert(fin3 && !fin3->IsZombie());

  TFile *fu = new TFile("rootfiles/unfold.root","READ");
  assert(fu && !fu->IsZombie());

  // Current uncertainty
  const char *p = "CondFormats/JetMETObjects/data/";
  const char *t18 = "Autumn18_RunD_V8_DATA";
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
    {1, 1, 6.9/3.0, 1,
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
     kRed+2, kOrange+2, kBlue+1, /*kGreen+2,*/ kBlack,
     kMagenta+2, kCyan+2, kGray+2};
  int marker[nera] =
    {kFullCircle, kFullDiamond, kFullStar, kFullSquare,
     kOpenCircle, kOpenDiamond, kOpenStar, kOpenSquare, /*kOpenTriangleDown,*/
     kFullCircle, kFullDiamond, kFullStar};//kFullSquare};
  const char* label[nera] =
    //{"RunA","RunB","RunC","RunD"};
    {"18A","18B","18C","18D",
     //"17B","17C","17D","17E","17F",
     "17B","17C","17DE","17F",
     "16BCD","16EF","16GH"};

  // Settings for the spectrum fit and JEC-equivalent plot y-axis range
  double eta(0.), ymin(-4), ymax(+6);
  if (sy=="0.0-0.5") { eta = 0;   ymin = -6; ymax = +9; }
  if (sy=="0.5-1.0") { eta = 0.5; ymin = -6; ymax = +9; }
  if (sy=="1.0-1.5") { eta = 1.0; ymin = -6; ymax = +9; }
  if (sy=="1.5-2.0") { eta = 1.5; ymin = -10; ymax = +15; }
  if (sy=="2.0-2.5") { eta = 2.0; ymin = -10; ymax = +15; }
  if (sy=="2.5-3.0") { eta = 2.5; ymin = -25; ymax = +45; }
  if (sy=="3.2-4.7") { eta = 3.2; ymin = -25; ymax = +45; }

  ///////////////////////////////////////////////////////////////////

  // Load new results
  TH1F *h1s[nera];
  for (int iera = 0; iera != nera; ++iera) {

    const char *cera = eras[iera].c_str(); 
    string hname = Form("ak4/Eta_%s/hpt_data_%s_det",cy,cera);
    TH1F *hera = (TH1F*)fins[iera]->Get(hname.c_str());

    // Engin's file has different naming scheme
    if (!hera) {
      int iy = int((eta+0.25)/0.5)+1;
      hname = Form("ak4/y_%s/hptData_%s_detector_%dbin",cy,cera,iy);
      hera = (TH1F*)fins[iera]->Get(hname.c_str());
    }
    if (!hera) cout << "Histogram " << hname << " not found!" << endl << flush;
    assert(hera);

    // fix a strange normalization bug for 16BCD, 16EF and 16GH
    string se = label[iera];
    if (se=="16BCD"||se=="16EF"||se=="16GH"||se=="17DE") {
      for (int i = 1; i != hera->GetNbinsX()+1; ++i) {
	double k = int((eta+0.25)/0.5)+1;
	hera->SetBinContent(i, hera->GetBinContent(i)*k);
	hera->SetBinError(i, hera->GetBinError(i)*k);
      }
    }

    // Correct ECAL prefire for 2016 and 2017 data at 2<|eta|<3
    if (correctECALprefire) {

      jer_iov run(run2018);
      //if (se=="16BCD"||se=="16EF"||se=="16GH")         run = run2016;
      //if (se=="17B"||se=="17C"||se=="17DE"||se=="17F") run = run2017;
      if (se=="16BCD") run = run2016bcd;
      if (se=="16EF")  run = run2016ef;
      if (se=="16GH")  run = run2016gh;
      if (se=="17B")   run = run2017b;
      if (se=="17C")   run = run2017c;
      if (se=="17DE")  run = run2017de;
      if (se=="17F")   run = run2017f;

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
    if (true) {
      //TF1 *f1 = (TF1*)fu->Get(Form("fr_%s",cy)); assert(f1);
      //TH1D *hr = (TH1D*)fu->Get(Form("hr_%s",cy)); assert(hr);
      int year(0);
      if (TString(se.c_str()).Contains("16")) year = 2016;
      if (TString(se.c_str()).Contains("17")) year = 2017;
      if (TString(se.c_str()).Contains("18")) year = 2018;
      // special case for 2018 2.5-3.0 bin
      if (TString(se.c_str()).Contains("18") && sy=="2.5-3.0") year = 2016;
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

    hera->Scale(1./lumi[iera]);
    h1s[iera] = hera;
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
    f1->SetLineColor(color[iera]);
    f1s[iera] = f1;
  } // for iera

  ///////////////////////////////////////////////////////////////////

  // Unfolded reference data from SMP-15-007 stored in HEPDATA
  TH1F *h50ns4(0), *h50ns7(0);
  {
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
    TH1F *hhepdy7 = (TH1F*)fhepd->Get(Form("%s/Hist1D_y1",ct7));
    assert(hhepdy7);
    TH1F *hhepde7 = (TH1F*)fhepd->Get(Form("%s/Hist1D_y1_e1",ct7));
    assert(hhepde7);
    for (int i = 1; i != hhepdy7->GetNbinsX()+1; ++i) {
      hhepdy7->SetBinError(i, hhepde7->GetBinContent(i));
    }
    h50ns7 = hhepdy7;
    //
    const char *ct4 = mhep4[sy].c_str();
    TH1F *hhepdy4 = (TH1F*)fhepd->Get(Form("%s/Hist1D_y1",ct4));
    assert(hhepdy4);
    TH1F *hhepde4 = (TH1F*)fhepd->Get(Form("%s/Hist1D_y1_e1",ct4));
    assert(hhepde4);
    for (int i = 1; i != hhepdy4->GetNbinsX()+1; ++i) {
      hhepdy4->SetBinError(i, hhepde4->GetBinContent(i));
    }
    //h50ns4 = hhepdy4;
    
    // Copy AK4 over to new binning, and set 30<pT<40 GeV from Run18D (was 18A) 
    // Also set/replace emax>2000 from Run18D to constrain very high pT
    // Constraints set at 30% level (xerr)
    // => replace Run2018D with Run2016BCD at low pT to reduce ZB bias
    h50ns4 = new TH1F("h50ns4",";p_{T} (GeV);Cross section",nptb,ptbins);
    for (int i = 1; i != h50ns4->GetNbinsX()+1; ++i) {
      int x = h50ns4->GetBinCenter(i);
      int j = hhepdy4->FindBin(x);
      int irund = 3; // Index of 2018D
      int irunbcd = 8; // Index of 2016BCD
      const double emax = 2000;
      const double xerr = 0.3; // 50%
      int k1 = h1s[irunbcd]->FindBin(x);
      int k2 = h1s[irund]->FindBin(x);
      if (hhepdy4->GetBinContent(j)!=0 && (x*cosh(eta)<emax)) {
	h50ns4->SetBinContent(i, hhepdy4->GetBinContent(j));
	h50ns4->SetBinError(i, hhepdy4->GetBinError(j));
      }
      else if ((x>=30 && x<=40) || (x*cosh(eta)>=emax)) {
	int irun =  (x<=40 ? irunbcd : irund);
	int k = (x<=40 ? k1 : k2);
	h50ns4->SetBinContent(i, h1s[irun]->GetBinContent(k));
	h50ns4->SetBinError(i, sqrt(pow(h1s[irun]->GetBinError(k),2)
				    +pow(h1s[irun]->GetBinContent(k)*xerr,2)));
      }
    }

  }
  assert(h50ns4);
  assert(h50ns7);
  
  // Fit old AK7 results with powerlaw
  TF1 *f50ns7 = new TF1("f50ns7",
			Form("[0]*pow(x,[1]+[2]*log(x))"
			     "*pow(1-2.*x*cosh(%1.2f)/13000.,[3])",eta),
			ptmin,emax/cosh(eta));
  f50ns7->SetParameters(1e12,-5,0.001,10);
  f50ns7->SetNpx(2640.*2); // this limits precision of DeltaJEC
  h50ns7->Fit(f50ns7,"RNI");
  f50ns7->SetLineColor(kGreen+2);

  // Fit old AK4 results with powerlaw
  TF1 *f50ns4 = new TF1("f50ns4",
			Form("[0]*pow(x,[1]+[2]*log(x))"
			     "*pow(1-2.*x*cosh(%1.2f)/13000.,[3])",eta),
			ptmin,emax/cosh(eta));
  f50ns4->SetParameters(1e12,-5,0.001,10);
  f50ns4->SetNpx(2640.*2); // this limits precision of DeltaJEC
  h50ns4->Fit(f50ns4,"RNI");
  f50ns4->SetLineColor(kGreen+2);

  //TF1 *fref = f50ns7; // compare to 50 ns data for AK7
  TF1 *fref = f50ns4; // compare to 50 ns data for AK4

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
    f1r->SetLineColor(color[iera]);
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

  // Divide new data by reference fit
  TH1F *hrs[nera];
  for (int iera = 0; iera != nera; ++iera) {

    //fref = f50ns4;
    TH1F *hr = (TH1F*)h1s[iera]->Clone(Form("hr_%d",iera));
    hr->Divide(fref);
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
  TH1F *hr50ns4 = (TH1F*)h50ns4->Clone("hr50ns4");
  hr50ns4->Divide(fref);//f50ns4);//fref);

  TH1F *hr50ns7 = (TH1F*)h50ns7->Clone("hr50ns7");
  hr50ns7->Divide(fref);//f50ns7);//fref);

  
  TH1D *h = new TH1D("h",";Jet p_{T} (GeV);"
		     "d#sigma^{2} / dy dp_{T} (pb / GeV)",
		     int(ptmax-ptmin),ptmin,ptmax);
  h->SetMinimum(xsecmin);//1e-6 *1.0001);
  h->SetMaximum(xsecmax);//1e+8 *0.9999);

  TH1D *h2 = new TH1D("h2",";Jet p_{T} (GeV);Ratio to 50 ns",
		      int(ptmax-ptmin),ptmin,ptmax);
  h2->SetMinimum(0.);
  h2->SetMaximum(2.0);
  h2->GetXaxis()->SetNoExponent();
  h2->GetXaxis()->SetMoreLogLabels();

  //lumi_13TeV = "Run2018ABCD 59.9 fb^{-1}";
  //lumi_13TeV = "Run2017BCDEF 41.4 fb^{-1}";
  //lumi_13TeV = "Run2016BCDEFGH 36.5 fb^{-1}";
  //lumi_13TeV = "Run2015 71 pb^{-1}";
  lumi_13TeV = "Run2 137.8 fb^{-1}";
  TCanvas *c1 = tdrDiCanvas("c1",h,h2,4,11);

  c1->cd(1);
  //tdrDraw(h50ns7,"Pz",kOpenSquare,kGreen+2,kSolid,kGreen+2);
  for (int iera = 0; iera != nera; ++iera) {
    tdrDraw(h1s[iera],"Pz",marker[iera],color[iera],kSolid,color[iera]);
  }
  tdrDraw(h50ns4,"Pz",kOpenDiamond,kGreen+2,kSolid,kGreen+2);
  gPad->SetLogx();
  gPad->SetLogy();
  
  //TLegend *leg = tdrLeg(0.50,0.90-0.06*nera,0.80,0.90);
  //for (int iera = 0; iera != nera; ++iera) {
  //leg->AddEntry(h1s[iera],label[iera],"PL");
  //}
  const int nmax = 6;
  TLegend *leg1 = tdrLeg(0.50,0.90-0.06*min(nmax,nera),0.80,0.90);
  TLegend *leg2 = tdrLeg(0.65,0.90-0.06*max(1,min(nmax,nera-nmax)+1),0.95,0.90);
  for (int iera = 0; iera != nera; ++iera) {
    if (iera<nmax)  leg1->AddEntry(h1s[iera],label[iera],"PL");
    if (iera>=nmax) leg2->AddEntry(h1s[iera],label[iera],"PL");
  }
  if (nera<nmax)  leg1->AddEntry(h50ns4,"2015 pub.","PL");
  if (nera>=nmax) leg2->AddEntry(h50ns4,"2015 pub.","PL");
  //leg->AddEntry(h50ns4,"74X 50 ns AK4","PL");
  //leg->AddEntry(h50ns7,"74X 50 ns AK7","PL");

  c1->cd(2);
  //tdrDraw(hr50ns7,"Pz",kOpenSquare,kGreen+2,kSolid,kGreen+2);
  for (int iera = 0; iera != nera; ++iera) {
    tdrDraw(hrs[iera],"Pz",marker[iera],color[iera],kSolid,color[iera]);
    hrs[iera]->GetXaxis()->SetRangeUser(ptmin,ptmax);
  }
  tdrDraw(hr50ns4,"Pz",kOpenDiamond,kGreen+2,kSolid,kGreen+2);
  gPad->SetLogx();

  c1->cd(1);
  f50ns4->SetRange(ptmin,0.9*emax/cosh(eta));//3500);
  f50ns4->DrawClone("SAME");
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
  tex->DrawLatex(0.20,0.05,Form("fit(50ns) AK4: #chi^{2}/NDF = %1.1f/%d",
				f50ns4->GetChisquare(),f50ns4->GetNDF()));
						
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.20,0.20,Form("y #in %s",cy)); 

  c1->cd(2);
  f50nsr->Draw("SAME");
  // Draw ratio of fits as well (or don't)
  for (int iera = 0; iera != nera; ++iera) {
    //f1rs[iera]->Draw("SAME");
  }

  c1->SaveAs(Form("pdf/drawDeltaJEC_jetPt_%s.pdf",cy2));


  TH1D *h3 = new TH1D("h3",";Jet p_{T} (GeV);#DeltaJEC-equivalent (%)",
		      int(ptmax-ptmin),ptmin,ptmax);
  h3->SetMinimum(ymin);
  h3->SetMaximum(ymax);
  h3->GetXaxis()->SetMoreLogLabels();
  h3->GetXaxis()->SetNoExponent();
  
  // Calculate reference uncertainty
  TH1F *hunc18 = (TH1F*)h50ns4->Clone("hunc18");
  TH1F *hunc17 = (TH1F*)h50ns4->Clone("hunc17");
  TH1F *hunc16 = (TH1F*)h50ns4->Clone("hunc16");
  TH1F *huncs1 = (TH1F*)h50ns4->Clone("huncs1");
  TH1F *huncs2 = (TH1F*)h50ns4->Clone("huncs2");
  TH1F *hunc15 = (TH1F*)h50ns4->Clone("hunc15");
  TH1F *hunc12 = (TH1F*)h50ns4->Clone("hunc12");
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

  TH1F *hds[nera];
  TH1F *huncl = (TH1F*)h50ns4->Clone("huncl");
  for (int iera = 0; iera != nera; ++iera) {
    
    //fref = f50ns4;
    TH1F *hr = hrs[iera];
    TH1F *hd = (TH1F*)hr->Clone(Form("hd_%d",iera));
    for (int i = 1; i != hr->GetNbinsX()+1; ++i) {

      if (hr->GetBinContent(i)==0) continue;
      double x = hr->GetBinCenter(i);
      double y = fref->Eval(x);
      double y2 = y * hr->GetBinContent(i);
      double x2 = fref->GetX(y2,0.85*x,1.15*x,1e-3);
      hd->SetBinContent(i, 100.*(x/x2-1));

      // Luminosity uncertainty from ABCD
      if (iera==0) {
	double y3 = y / 1.027; // lumi uncertainty for 2016, after 3.3% patch
	double x3 = fref->GetX(y3,0.85*x,1.15*x,1e-3);
	int j = huncl->FindBin(hd->GetBinCenter(i));
	huncl->SetBinContent(j, 0);
	huncl->SetBinError(j, 100.*fabs(x/x3-1));
      }

      double y2_up = y * (hr->GetBinContent(i) + hr->GetBinError(i));
      double x2_dw = fref->GetX(y2_up,0.85*x,1.15*x);
      double y2_dw = y * (hr->GetBinContent(i) - hr->GetBinError(i));
      double x2_up = fref->GetX(y2_dw,0.85*x,1.15*x);
      hd->SetBinError(i, 100.*sqrt(pow((x2_up-x2_dw)/x2,2) + pow(0.001,2)));
    } // for i
    hds[iera] = hd;

  } // for iera


  TCanvas *c2 = tdrCanvas("c2",h3,4,11,kSquare);

  tdrDraw(hunc18,"E3", kSolid, kBlack,kSolid,-1,1001,kYellow+1); // 2018
  tdrDraw(hunc17,"E3", kSolid, kBlack,kSolid,-1,1001,kCyan-6); // 2017
  hunc17->SetFillColorAlpha(kCyan-6,0.70);
  //tdrDraw(hunc16,"E3", kSolid, kBlack,kSolid,-1,1001,kGreen-6); // 2016
  //hunc16->SetFillColorAlpha(kGreen-6,0.70);
  tdrDraw(hunc16,"E3", kSolid, kBlack,kSolid,-1,1001,kViolet-8); // 2016
  hunc16->SetFillColorAlpha(kViolet-8,0.70);

  tdrDraw(hunc15,"E3", kNone, kBlack,kSolid,-1,1001,kOrange-9); // 2015
  hunc15->SetFillColorAlpha(kOrange-9,0.70);
  tdrDraw(hunc12,"E3", kNone, kBlack,kSolid,-1,1001,kRed-9); // 2012
  hunc12->SetFillColorAlpha(kRed-9,0.70);
  //tdrDraw(huncl,"E3", kNone, kBlack,kSolid,-1,1001,kBlue-9);
  //tdrDraw(hunc,"E3", kSolid, kBlack,kSolid,-1,3001,kYellow+1);
  tdrDraw(huncs1,"HIST][",kSolid,kBlack,kSolid,-1,kNone,kOrange-9); //RelPtBal
  tdrDraw(huncs2,"HIST][",kSolid,kBlue,kSolid,-1,kNone,kOrange-9); //ReltPtEC1+2
  
  for  (int iera = 0; iera != nera; ++iera) {
    tdrDraw(hds[iera],"Pz",marker[iera],color[iera],kSolid,color[iera]);
    hds[iera]->GetXaxis()->SetRangeUser(ptmin,ptmax);
  }
  gPad->SetLogx();

  //TLegend *legd = tdrLeg(0.50,0.90-0.04*(nera+3),0.80,0.90);
  //for  (int iera = 0; iera != nera; ++iera) {
  //legd->AddEntry(hrs[iera],label[iera],"PL");
  //}
  const int nmaxd = 6;
  TLegend *legd1 = tdrLeg(0.50,0.9-0.04*min(nmaxd,nera),0.80,0.9);
  TLegend *legd2 = tdrLeg(0.65,0.9-0.04*max(0,min(nmaxd,nera-nmaxd)),0.95,0.9);
  for  (int iera = 0; iera != nera; ++iera) {
    if (iera<nmaxd)  legd1->AddEntry(hrs[iera],label[iera],"PL");
    if (iera>=nmaxd) legd2->AddEntry(hrs[iera],label[iera],"PL");
  }
  
  //TLegend *legu = tdrLeg(0.50,0.15,0.80,0.23);
  const int nmaxu1 = 3;
  const int nmaxu2 = 3;
  TLegend *legu1 = tdrLeg(0.20,0.15,0.50,0.15+0.04*nmaxu1);
  TLegend *legu2 = tdrLeg(0.60,0.15,0.90,0.15+0.04*nmaxu2);
  legu1->AddEntry(hunc18,"2018 Aut18_V8","F");
  //legu1->AddEntry(hunc17,"2017 17Nov_V11","F");
  legu1->AddEntry(hunc17,"2017 17Nov_V32","F");
  //legu1->AddEntry(hunc16,"2016 07Aug_V17","F");
  legu1->AddEntry(hunc16,"2016 07Aug_V11","F");
  legu2->AddEntry(hunc15,"2015 Sum16_50nsV5","F");
  legu2->AddEntry(hunc12,"2012 Win14_V8","F");
  //legd->AddEntry(huncl,"2.6% lum.","F");
  
  tex->SetTextColor(kBlack); tex->SetTextSize(0.040);
  tex->DrawLatex(0.19,0.75,Form("Anti-k_{T} R=0.4"));
  tex->DrawLatex(0.19,0.71,Form("PF+CHS unf."));
  tex->DrawLatex(0.19,0.67,Form("|y|#in %s",cy));
  tex->SetTextSize(0.025);
  tex->DrawLatex(0.19,0.64,Form("vs HEPData-ins1459051-v1"));
  tex->SetTextSize(0.040);
  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(ptmin,0,ptmax,0);
  l->SetLineStyle(kDotted);
  l->DrawLine(ptmin,+1.5,ptmax,+1.5);
  l->DrawLine(ptmin,-1.5,ptmax,-1.5);

  c2->RedrawAxis();
  c2->SaveAs(Form("pdf/drawDeltaJEC_%s.pdf",cy2));
}


// Ansatz Kernel
//const double emax = 6500.;
int cnt_a = 0;
const int nk = 3; // number of kernel parameters (excluding pt, eta)
Double_t smearedAnsatzKernel(Double_t *x, Double_t *p) {

  if (++cnt_a%1000000==0) {
    cout << "+" << flush;
  }

  const double pt = x[0]; // true pT
  const double ptmeas = p[0]; // measured pT
  const double eta = p[4]; // rapidity

  double res = ptresolution(pt, eta+1e-3) * pt;
  const double s = TMath::Gaus(ptmeas, pt, res, kTRUE);
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
void unfold(string sy = "0.0-0.5", int ieta = 1) {


  TFile *f = new TFile("rootfiles/common2018_V7.root","READ");
  //TFile *f = new TFile("rootfiles/common2016_LegacyIOVs_v3.root","READ");
  //TFile *f = new TFile("rootfiles/common2016_October2018_V17.root","READ");
  assert(f && !f->IsZombie());
  TFile *fout = new TFile("rootfiles/unfold.root",
			  ieta==1 ? "RECREATE" : "UPDATE");

  const char *cy = sy.c_str();
  //TH1D *hd = (TH1D*)f->Get(Form("ak4/y_%s/hptData_full2016_detector_%dbin",
  //cy,ieta));
  TH1D *hd = (TH1D*)f->Get(Form("ak4/Eta_%s/hpt_data_2018_D_det",cy));

  assert(hd);
  //TH1D *hu = (TH1D*)f->Get(Form("ak4/y_%s/hptData_full2016_particle_%dbin",
  //				cy,ieta));
  TH1D *hu = (TH1D*)f->Get(Form("ak4/Eta_%s/hpt_data_2018_D_det",cy));
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
  int i = ieta-1;
  double eta = p[i][nk];
  double maxpt = emax/cosh(eta);

  _ismcjer = false;
  _usejme = false;
  _jer_iov = run2016;
  
  // Initial fit of the NLO curve to a histogram
  TF1 *fus = new TF1(Form("fus%s",cy),
		     "[0]*pow(x,[1])"
		     "*pow(1-x*cosh([3])/6500.,[2])",
		     5,maxpt);
  fus->SetParameters(p[i][0], p[i][1], p[i][2], p[i][3]);
  
  // Smeared spectrum
  TF1 *fs = new TF1(Form("fs%s",cy),smearedAnsatz,5.,maxpt,nk+1);
  fs->SetParameters(fus->GetParameter(0), fus->GetParameter(1),
		    fus->GetParameter(2), fus->GetParameter(3));
  

  fout->cd();

  const int niov = 4;
  jer_iov iovs[niov] = {run1, run2016, run2017, run2018};
  for (int iov = 0; iov != niov; ++iov) {

    _jer_iov = iovs[iov]; // use by ptresolution in fs and fus
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
