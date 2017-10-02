#include "TCanvas.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TLine.h"

#include "TLatex.h"
//#include "tdrstyle_mod12.C"
#include "tdrstyle_mod14.C"
//#include "settings12.h"

#include <fstream>
#include <map>

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

// symmetric if both below false
// positive true overrides negative
bool _usenegative = false;
bool _usepositive = false;

bool _useptgen = true; // iterate to produce JEC vs pTgen
bool _dothree  = false;//true; // compare three JECs instead of just two
bool _paper    = true; // graphical settings for the paper (e.g. y-axis range)

const double _mu = 24.68;//12.8;//20;//19.83; // 20/fb at 8 TeV (htrpu)
//const double _lumi = 19800.;
bool _pdf = true; // save .pdf
bool _C   = false;//true; // save .C
bool _mc(false);
string _alg("");

// Mapping from mu to NPV and rho
// Derived from inclusive jet data, |eta|<1.3, jt320: prhovstrpu, pnpvvstrpu
// TF1 *f1 = new TF1("f1","[0]+[1]*x+[2]*x*x",7,33);
// for (int i=0; i!=f1->GetNpar(); ++i) cout<<Form("%1.4g, ",f1->GetParameter(i)); cout<<endl;
//
// prhovstrpu->Fit(f1,"R")
double getRho(double mu) {

  double p[3] = {1.009, 0.5515, 0.0003597}; // DATA

  return (p[0]+p[1]*mu+p[2]*mu*mu);

}
// pnpvvstrpu->Fit(f1,"R")
double getNPV(double mu) {

  double p[3] = {1.032, 0.7039, -0.001319}; // DATA

  return (p[0]+p[1]*mu+p[2]*mu*mu);
}

void setEtaPtE(FactorizedJetCorrector *jec, double eta, double pt, double e,
	       int mu) {

  assert(jec);
  
  int npv = getNPV(mu);
  double rho = getRho(mu);

  jec->setJetEta(eta);
  jec->setJetPt(pt);
  // L1Offset
  jec->setJetE(e);
  jec->setNPV(npv);
  // L1FastJet
  bool is5 = (_alg=="AK5PF"||_alg=="AK5PFchs"||_alg=="AK5CALO");
  bool is4 = (_alg=="AK4PF"||_alg=="AK4PFchs"||_alg=="AK4CALO");
  double jeta = TMath::Pi()*(is5 ? 0.5*0.5 : (is4 ? 0.4*0.4 : 0.7*0.7));
  jec->setJetA(jeta);
  jec->setRho(rho);

  return;
}

FactorizedJetCorrector *_thejec(0);
TF1 *fCorrPt(0);
Double_t funcCorrPt(Double_t *x, Double_t *p) {
  
  double eta = p[0];
  double pt = x[0];
  double e = pt * cosh(eta);
  double mu = p[1];
  setEtaPtE(_thejec, eta, pt, e, mu);

  return (_thejec->getCorrection() * pt);
}

double getEtaPtE(FactorizedJetCorrector *jec, double eta, double pt, double e,
		 int mu = _mu) {

  setEtaPtE(jec, eta, pt, e, mu);

  // if using pTgen, need to iterate to solve ptreco
  if (_useptgen) {

    double ptgen = pt;
    _thejec = jec;
    fCorrPt->SetParameters(eta, mu);
    // Find ptreco that gives pTreco*JEC = pTgen
    double ptreco = fCorrPt->GetX(ptgen,5,6500);

    setEtaPtE(jec, eta, ptreco, e, mu);
  }

  return (jec->getCorrection());
} // getEtaPtE

double getEtaPtUncert(JetCorrectionUncertainty *unc,
		      FactorizedJetCorrector *jec,
		      double eta, double pt, double mu = _mu) {

  assert(unc);
  
  unc->setJetEta(eta);
  unc->setJetPt(pt);

  // if not using pTgen, need to solve it
  if (!_useptgen) {

    assert(jec);
    double ptreco = pt;
    _thejec = jec;
    fCorrPt->SetParameters(eta, mu);
    double ptgen = fCorrPt->Eval(ptreco);

    unc->setJetEta(eta);
    unc->setJetPt(ptgen);
  }

  return (unc->getUncertainty(true));
} // getEtaPtUncert


void compareJECversions(string algo="AK4PFchs",
			bool l1=true, bool l2l3=true, bool res=true,
			string type="DATA") {

  setTDRStyle();
  writeExtraText = false; // for JEC paper CWR

  assert(type=="DATA" || type=="MC");
  const bool mc = (type=="MC");
  _mc = mc;
  _alg = algo;
  assert(!mc || !res);
  const char *a = algo.c_str();
  const char *str;

  // Initialize function used to invert JEC
  fCorrPt = new TF1("fCorrPt",funcCorrPt,1,6500,2);
  //bool dothree = (_dothree && !(l1 && TString(a).Contains("chs")));
  bool dothree = _dothree;

  // Legends
  const char *cm0 = "DATA";
  const char *cm = type.c_str();
  string sgen = (_useptgen ? "corr" : "raw");
  const char *cgen = sgen.c_str();

  // 2015 JEC, 76X
  //string sid2 = (_mc ? "Summer15_25nsV6_MC" : "Summer15_25nsV6_DATA");
  //string sid2 = (_mc ? "Fall15_25nsV1_MC" : "Fall15_25nsV1_DATA");
  //string sid2 = (_mc ? "Spring16_25nsV8p2_MC" : "Spring16_25nsV8p2_DATA");
  //string sid2 = (_mc ? "Spring16_23Sep2016GHV1_MC" : "Spring16_23Sep2016GHV1_DATA");
  //string sid2 = (_mc ? "Spring16_23Sep2016BCDV1_MC" : "Spring16_23Sep2016BCDV1_DATA");
  //string sid2 = (_mc ? "Summer16_23Sep2016BCDV1_MC" : "Summer16_23Sep2016BCDV1_DATA");
  //string sid2 = (_mc ? "Summer16_23Sep2016BCDV1_MC" : "Summer16_23Sep2016BCDV1_DATA");
  //string sid2 = (_mc ? "Summer16_23Sep2016GV2_MC" : "Summer16_23Sep2016GV2_DATA");
  //string sid2 = (_mc ? "Summer16_23Sep2016V2_MC" : "Summer16_23Sep2016GV3_DATA");
  //string sid2 = (_mc ? "Summer16_23Sep2016V3_MC" : "Summer16_23Sep2016HV3_DATA");
  string sid2 = (_mc ? "Summer16_03Feb2017_V1_MC" : "Summer16_03Feb2017H_V3_DATA");
  const char *cid2 = sid2.c_str();
  const char *a2 = a;
  //const char *s2 = "1.3 fb^{-1} (13 TeV)";
  //const char *s2s = "2012";
  //const char *s2 = "2.1 fb^{-1} (13 TeV)";
  //const char *s2 = "76Xv1 (13 TeV)";
  //const char *s2s = "76X";
  //const char *s2 = "80Xv8 G";// (13 TeV)";
  //const char *s2 = "80XreV1 GH";// (13 TeV)";
  //const char *s2s = "reGH";
  //const char *s2 = "80XreV1 BCD";
  //const char *s2s = "reV1 BCD";
  //const char *s2 = "80Xre Sum16";
  //const char *s2s = "Sum16";
  //const char *s2 = "80X Sum16 BCD";
  //const char *s2s = "BCD";
  //const char *s2 = "80X Sum16 G";
  //const char *s2 = "Summer16GV3";
  //const char *s2s = "Sum16V3";
  const char *s2 = "03Feb2017H_V3";
  const char *s2s = "03Feb17HV3";
  // PATCH 2012 with clones
  //if (algo=="AK4PF") a2 = "AK5PF";
  //if (algo=="AK4PFchs") a2 = "AK5PFchs";

  // 2012 JEC
  //string sid1 = (_mc ? "Winter14_V8_MC" : "Winter14_V8_DATA");
  // 74X JEC
  //string sid1 = (_mc ? "Summer15_25nsV7_MC" : "Summer15_25nsV7_DATA");
  //string sid1 = (_mc ? "Spring16_25nsV8BCD_MC" : "Spring16_25nsV8BCD_DATA");
  //string sid1 = (_mc ? "Spring16_23Sep2016BCDV1_MC" : "Spring16_23Sep2016BCDV1_DATA");
  //string sid1 = (_mc ? "Summer16_23Sep2016EFV1_MC" : "Summer16_23Sep2016EFV1_DATA");
  //string sid1 = (_mc ? "Summer16_23Sep2016EFV2_MC" : "Summer16_23Sep2016EFV2_DATA");
  //string sid1 = (_mc ? "Spring16_23Sep2016V1_MC" : "Spring16_23Sep2016GHV1_DATA");
  //string sid1 = (_mc ? "Spring16_23Sep2016V1_MC" : "Spring16_23Sep2016GV2_DATA");
  //string sid1 = (_mc ? "Summer16_23Sep2016V3_MC" : "Summer16_23Sep2016HV3_DATA");
  //string sid1 = (_mc ? "Summer16_03Feb2017_V3_MC" : "Summer16_03Feb2017G_V3_DATA");
  string sid1 = (_mc ? "Summer16_23Sep2016_V3_MC" : "Summer16_23Sep2016HV3_DATA");
  const char *cid1 = sid1.c_str();
  //const char *a1 = "AK5PFchs";//a;
  const char *a1 = "AK4PFchs";//a;
  //const char *s1 = "23Sep2016HV3";
  //const char *s1s = "23Sep16HV3";
  //const char *s1 = "03Feb2017G_V3";
  //const char *s1s = "03Feb17G_V3";
  const char *s1 = "23Sep2016H_V3";
  const char *s1s = "23Sep16H_V3";
  //const char *s1 = "20 fb^{-1} (8 TeV)";
  //const char *s1s = "2012";
  //const char *s1 = "1.3 fb^{-1} (13 TeV)";
  //const char *s1 = "74Xv7 (13 TeV)";
  //const char *s1s = "74X";
  //const char *s1 = "80XprV8 BCD";// (13 TeV)";
  //const char *s1 = "80Xre Spr16";// (13 TeV)";
  //const char *s1 = "80XreV1 BCD";// (13 TeV)";
  //const char *s1s = "prV8 BCD";
  //const char *s1s = "Spr16";
  //const char *s1 = "80X Sum16 EFearly";
  //const char *s1s = "EF";
  //const char *s1 = "80X Spr16 GH";
  //const char *s1s = "Spr16";
  // PATCH 2012 with clones
  //if (algo=="AK4PF") a1 = "AK5PF";
  //if (algo=="AK4PFchs") a1 = "AK5PFchs";

  // 2012 JEC
  //string sid3 = (_mc ? "Winter14_V8_MC" : "Winter14_V8_DATA");
  string sid3 = (_mc ? "Summer16_03Feb2017_V1_MC" : "Summer16_03Feb2017H_V3_DATA");
  //string sid3 = (_mc ? "Spring16_25nsV8E_MC" : "Spring16_25nsV8E_DATA");
  //string sid3 = (_mc ? "Spring16_25nsV6_MC" : "Spring16_25nsV6_DATA");
  //string sid3 = (_mc ? "Spring16_25nsV8F_MC" : "Spring16_25nsV8F_DATA");
  //string sid3 = (_mc ? "Spring16_25nsV8E_MC" : "Spring16_25nsV8E_DATA");
  //string sid3 = (_mc ? "Summer15_50nsV4_MC" : "Summer15_50nsV4_DATA");
  //string sid3 = (_mc ? "Spring16_23Sep2016EV1_MC" : "Spring16_23Sep2016EV1_DATA");
  //string sid3 = (_mc ? "Summer16_23Sep2016GV1_MC" : "Summer16_23Sep2016GV1_DATA");
  //string sid3 = (_mc ? "Summer16_23Sep2016GV2_MC" : "Summer16_23Sep2016GV2_DATA");
  //string sid3 = (_mc ? "Summer16_23Sep2016V2_MC" : "Summer16_23Sep2016BCDV3_DATA");
  const char *cid3 = sid3.c_str();
  const char *a3 = a;
  //const char *a3 = "AK5PFchs"; // for Winter14
  //const char *s3 = "20 fb^{-1} (8 TeV)"; // for Winter14
  //const char *s3 = "R=0.5, 20 fb^{-1} (8 TeV)"; // for Winter14
  //const char *s3 = "R=0.5, Winter14_V8"; // for Winter14
  const char *s3 = "03Feb2017H_V3";
  const char *s3s = "03Feb17H_V3";
  //const char *s3s = "2012";
  //const char *s3s = "74X";
  //const char *s3s = "Run I"; // for Winter14
  //const char *s3 = "80Xv8 E";// (13 TeV)";
  //const char *s3s = "E";
  //const char *s3 = "80Xv6 BCD";// (13 TeV)";
  //const char *s3s = "V6";
  //const char *s3 = "80Xv8 F";// (13 TeV)";
  //const char *s3s = "F";
  //const char *s3 = "50ns v4";
  //const char *s3s = "50ns";
  //const char *s3 = "80XreV1 E";
  //const char *s3s = "reE";
  //const char *s3 = "80Xv8 E";// (13 TeV)";
  //const char *s3s = "E";
  //const char *s3 = "80X Sum16 FlateG";
  //const char *s3s = "G";
  //const char *s3 = "80X Sum16 BCD";
  //const char *s3s = "BCD";
  //if (algo=="AK4PF") a3 = "AK5PF";
  //if (algo=="AK4PFchs") a3 = "AK5PFchs";

  // 2011 JEC
  //string sid3 = "GR_R_42_V23";
  //const char *cid3 = sid3.c_str();
  //const char *a3 = "AK5PFchs";//a;
  //const char *s3 = "5 fb^{-1} (7 TeV)";
  //const char *s3s = "2011";

  // 2010 JEC
  //string sid3 = "START38_V13";
  //const char *cid3 = sid3.c_str();
  //const char *a3 = a;
  //const char *s3 = "36 pb^{-1} (7 TeV)";
  //const char *s3s = "2010";
  // PATCH 2010 with clones AK7PF/PFchs
  //if (algo=="AK7PF") a3 = "AK5PF";
  //if (algo=="AK7PFchs") a3 = "AK5PFchs";


  str=Form("CondFormats/JetMETObjects/data/%s_L1FastJet_%s.txt",cid1,a1);
  //str=Form("CondFormats/JetMETObjects/data/%s_L1RC_%s.txt",cid1,a1);
  cout << str << endl << flush;
  JetCorrectorParameters *JetCorPar1L1 = new JetCorrectorParameters(str);
  str=Form("CondFormats/JetMETObjects/data/%s_L2Relative_%s.txt",cid1,a1);
  cout << str << endl << flush;
  JetCorrectorParameters *JetCorPar1L2 = new JetCorrectorParameters(str);
  str=Form("CondFormats/JetMETObjects/data/%s_L3Absolute_%s.txt",cid1,a1);
  cout << str << endl << flush;
  JetCorrectorParameters *JetCorPar1L3 = new JetCorrectorParameters(str);
  str=Form("CondFormats/JetMETObjects/data/%s_L2L3Residual_%s.txt",cid1,a1);
  if (!mc) cout << str << endl << flush;
  JetCorrectorParameters *JetCorPar1 = (mc ? 0 : new JetCorrectorParameters(str));
  str=Form("CondFormats/JetMETObjects/data/%s_Uncertainty_%s.txt",cid1,a1);
  cout << str << endl << flush;
  JetCorrectionUncertainty *jecUnc1 = new JetCorrectionUncertainty(str);

  str=Form("CondFormats/JetMETObjects/data/%s_L1FastJet_%s.txt",cid2,a2);
  //str=Form("CondFormats/JetMETObjects/data/%s_L1RC_%s.txt",cid2,a2);
  cout << str << endl << flush;
  JetCorrectorParameters *JetCorPar2L1 = new JetCorrectorParameters(str);
  str=Form("CondFormats/JetMETObjects/data/%s_L2Relative_%s.txt",cid2,a2);
  cout << str << endl << flush;
  JetCorrectorParameters *JetCorPar2L2 = new JetCorrectorParameters(str);
  str=Form("CondFormats/JetMETObjects/data/%s_L3Absolute_%s.txt",cid2,a2);
  cout << str << endl << flush;
  JetCorrectorParameters *JetCorPar2L3 = new JetCorrectorParameters(str);
  str=Form("CondFormats/JetMETObjects/data/%s_L2L3Residual_%s.txt",cid2,a2);
  if (!mc) cout << str << endl << flush;
  JetCorrectorParameters *JetCorPar2 = (mc ? 0 : new JetCorrectorParameters(str));
  str=Form("CondFormats/JetMETObjects/data/%s_Uncertainty_%s.txt",cid2,a2);
  cout << str << endl << flush;
  JetCorrectionUncertainty *jecUnc2 = new JetCorrectionUncertainty(str);

  vector<JetCorrectorParameters> vParam1;
  if (l1)   vParam1.push_back(*JetCorPar1L1);
  if (l2l3) vParam1.push_back(*JetCorPar1L2);
  if (l2l3) vParam1.push_back(*JetCorPar1L3);
  if (res && !mc && JetCorPar1)  vParam1.push_back(*JetCorPar1);
  vector<JetCorrectorParameters> vParam2;
  if (l1)   vParam2.push_back(*JetCorPar2L1);
  if (l2l3) vParam2.push_back(*JetCorPar2L2);
  if (l2l3) vParam2.push_back(*JetCorPar2L3);
  if (res && !mc && JetCorPar2)  vParam2.push_back(*JetCorPar2);

  FactorizedJetCorrector *JEC1 = new FactorizedJetCorrector(vParam1);
  FactorizedJetCorrector *JEC2 = new FactorizedJetCorrector(vParam2);

  FactorizedJetCorrector *JEC3(0);
  JetCorrectionUncertainty *jecUnc3(0);
  if (dothree) {
    // Note the reversed naming scheme in 2010
    /*
    str=Form("CondFormats/JetMETObjects/data/%s_%s_L1FastJet.txt",cid3,a3);
    cout << str << endl << flush;
    JetCorrectorParameters *JetCorPar3L1 = new JetCorrectorParameters(str);
    str=Form("CondFormats/JetMETObjects/data/%s_%s_L2Relative.txt",cid3,a3);
    cout << str << endl << flush;
    JetCorrectorParameters *JetCorPar3L2 = new JetCorrectorParameters(str);
    str=Form("CondFormats/JetMETObjects/data/%s_%s_L3Absolute.txt",cid3,a3);
    cout << str << endl << flush;
    JetCorrectorParameters *JetCorPar3L3 = new JetCorrectorParameters(str);
    str=Form("CondFormats/JetMETObjects/data/%s_%s_L2L3Residual.txt",cid3,a3);
    if (!mc) cout << str << endl << flush;
    JetCorrectorParameters *JetCorPar3 = (mc ? 0 : new JetCorrectorParameters(str));
    str=Form("CondFormats/JetMETObjects/data/%s_%s_Uncertainty.txt",cid3,a3);
    cout << str << endl << flush;
    jecUnc3 = new JetCorrectionUncertainty(str);
    */

    str=Form("CondFormats/JetMETObjects/data/%s_L1FastJet_%s.txt",cid3,a3);
    //str=Form("CondFormats/JetMETObjects/data/%s_RC_%s.txt",cid3,a3); // L1RC
    cout << str << endl << flush;
    JetCorrectorParameters *JetCorPar3L1 = new JetCorrectorParameters(str);
    str=Form("CondFormats/JetMETObjects/data/%s_L2Relative_%s.txt",cid3,a3);
    cout << str << endl << flush;
    JetCorrectorParameters *JetCorPar3L2 = new JetCorrectorParameters(str);
    str=Form("CondFormats/JetMETObjects/data/%s_L3Absolute_%s.txt",cid3,a3);
    cout << str << endl << flush;
    JetCorrectorParameters *JetCorPar3L3 = new JetCorrectorParameters(str);
    str=Form("CondFormats/JetMETObjects/data/%s_L2L3Residual_%s.txt",cid3,a3);
    if (!mc) cout << str << endl << flush;
    JetCorrectorParameters *JetCorPar3 = (mc ? 0 : new JetCorrectorParameters(str));
    str=Form("CondFormats/JetMETObjects/data/%s_Uncertainty_%s.txt",cid3,a3);
    cout << str << endl << flush;
    jecUnc3 = new JetCorrectionUncertainty(str);

    vector<JetCorrectorParameters> vParam3;
    if (l1)   vParam3.push_back(*JetCorPar3L1);
    if (l2l3) vParam3.push_back(*JetCorPar3L2);
    if (l2l3) vParam3.push_back(*JetCorPar3L3);
    if (res && !mc && JetCorPar3)  vParam3.push_back(*JetCorPar3);
    JEC3 = new FactorizedJetCorrector(vParam3);

    assert(JEC3);
    assert(jecUnc3);
  }

  //_JEC1 = JEC1;

  TCanvas *c0 = new TCanvas(Form("c0_%s",a),Form("c0_%s",a),600,600);

  TCanvas *c2 = new TCanvas(Form("c2_%s",a),Form("c2_%s",a),600,600);

  TH1D *h = new TH1D(Form("h_%s",a),Form(";|#eta|;%s L2L3 residual",a),
		     50,0,5);
  if (_usenegative) h->GetXaxis()->SetTitle("-|#eta|");
  if (_usepositive) h->GetXaxis()->SetTitle("+|#eta|");
  const char *cl1 = (l1 ? "L1" : "");
  const char *cl2l3 = (l2l3 ? "L2L3" : "");
  const char *cpl = (res&&(l1||l2l3) ? "+" : "");
  const char *cplus = (res&&(l1||l2l3) ? "Plus" : "");
  const char *cres = (res ? "L2L3res" : "");
  
  // Create suitable binning
  const double x_pt[] =
    {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
     97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 362, 430,
     507, 592, 686, 790, 905, 1032, 1172, 1327, 1497, 1684,
     //_paper ? 1999. : 1890.,
     2000, 2238, 2500, 2787, 3103, 3450};
  const int ndiv_pt = sizeof(x_pt)/sizeof(x_pt[0])-1;// - (_paper ? 3 : 0);
  TH1D *hpt = new TH1D(Form("hpt_%s",a),
		       Form(";p_{T,%s} (GeV);%s L2L3 residual",cgen,a),
		       ndiv_pt, x_pt);
  hpt->GetYaxis()->SetTitle(Form("%s%s%s%s",cl1,cl2l3,cpl,cres));
  hpt->GetXaxis()->SetMoreLogLabels();
  hpt->GetXaxis()->SetNoExponent();
  hpt->SetMinimum(0.);
  hpt->SetMaximum(2.);

  h->GetYaxis()->SetTitle(Form("%s%s%s%s",cl1,cl2l3,cpl,cres));
  if (_paper) {
    if (l1 && !l2l3 && !res) h->SetYTitle("Pileup offset correction");
    //if (l1 && !l2l3 && !res) h->SetYTitle("Random cone offset correction"); // L1RC
    if (!l1 && l2l3 && !res) h->SetYTitle("Simulated response correction");
    if (!l1 && !l2l3 && res) h->SetYTitle("Residual response correction");
    //
    if (l1 && !l2l3 && !res) hpt->SetYTitle("Pileup offset correction");
    //if (l1 && !l2l3 && !res) hpt->SetYTitle("Random cone offset correction"); // L1RC
    if (!l1 && l2l3 && !res) hpt->SetYTitle("Simulated response correction");
    if (!l1 && !l2l3 && res) hpt->SetYTitle("Residual response correction");
  }
  h->SetMinimum(0.3);
  h->SetMaximum(2.0);
  if (_paper) {
    if (l1 && !l2l3 && !res) h->GetYaxis()->SetRangeUser(0.65,1.25);
    if (!l1 && l2l3 && !res) h->GetYaxis()->SetRangeUser(0.85,1.8);
    if (!l1 && !l2l3 && res) h->GetYaxis()->SetRangeUser(0.85,1.45);
    //
    if (l1 && !l2l3 && !res) hpt->GetYaxis()->SetRangeUser(0.5,1.4);
    if (!l1 && l2l3 && !res) hpt->GetYaxis()->SetRangeUser(0.85,1.6);
    //if (!l1 && !l2l3 && res) hpt->GetYaxis()->SetRangeUser(0.85,1.45);
    //if (!l1 && !l2l3 && res) hpt->GetYaxis()->SetRangeUser(0.80,1.20);
    //if (!l1 && !l2l3 && res) hpt->GetYaxis()->SetRangeUser(0.96,1.06);
    //if (!l1 && !l2l3 && res) hpt->GetYaxis()->SetRangeUser(0.95,1.15);
    //if (!l1 && !l2l3 && res) hpt->GetYaxis()->SetRangeUser(0.95,1.06);
if (!l1 && !l2l3 && res) hpt->GetYaxis()->SetRangeUser(0.96,1.07);
    //
    if (l1 && !l2l3 && !res) hpt->GetXaxis()->SetRangeUser(10,1999);
    if (!l1 && l2l3 && !res) hpt->GetXaxis()->SetRangeUser(10,1999);
    //if (!l1 && !l2l3 && res) hpt->GetXaxis()->SetRangeUser(10,1999);
    //if (!l1 && !l2l3 && res) hpt->GetXaxis()->SetRangeUser(10,3500);
    if (!l1 && !l2l3 && res) hpt->GetXaxis()->SetRangeUser(50,3500);
  }

  //lumi_7TeV  = (dothree ? "36 pb^{-1} + 4.9 fb^{-1}" : "4.9 fb^{-1}");
  //lumi_13TeV  = "19.8 fb^{-1} (8 TeV) + 1.3--2.1 fb^{-1}";
  //lumi_13TeV  = "27 fb^{-1} (13 TeV)";
  //lumi_13TeV  = "2016 re-reco 36.5 fb^{-1}";
  lumi_13TeV  = "2016 re-mAOD 36.5 fb^{-1}";

  TH1D *h1a = (TH1D*)h->Clone(Form("h1a_%s",a));
  TCanvas *c1a = tdrCanvas(Form("c1a_%s",a),h1a,4,11,kSquare);
  TH1D *h1b = (TH1D*)h->Clone(Form("h1b_%s",a));
  TCanvas *c1b = tdrCanvas(Form("c1b_%s",a),h1b,4,11,kSquare);
  TH1D *h1c = (TH1D*)h->Clone(Form("h1c_%s",a));
  TCanvas *c1c = tdrCanvas(Form("c1c_%s",a),h1c,4,11,kSquare);
  TH1D *h1e = (TH1D*)h->Clone(Form("h1e_%s",a));
  TCanvas *c1e = tdrCanvas(Form("c1e_%s",a),h1e,4,11,kSquare);

  TCanvas *c1d = tdrCanvas(Form("c1d_%s",a),hpt,4,11,kSquare);
  if (_paper) hpt->SetTitleOffset(0.97); // comma otherwise cut off

  TGraph *g1a = new TGraph(0);
  TGraph *g1b = new TGraph(0);
  TGraph *g1c = new TGraph(0);
  TGraph *g1d = new TGraph(0);
  TGraph *g1e = new TGraph(0);
  //
  TGraph *g2a = new TGraph(0);
  TGraph *g2b = new TGraph(0);
  TGraph *g2c = new TGraph(0);
  TGraph *g2d = new TGraph(0);
  TGraph *g2e = new TGraph(0);
  //
  TGraph *g3a = new TGraph(0);
  TGraph *g3b = new TGraph(0);
  TGraph *g3c = new TGraph(0);
  TGraph *g3d = new TGraph(0);
  TGraph *g3e = new TGraph(0);
  //
  TGraph *g21a = new TGraph(0);
  TGraph *g21b = new TGraph(0);
  TGraph *g21c = new TGraph(0);
  TGraph *g21e = new TGraph(0);

  TGraphErrors *g1a_e = new TGraphErrors(0);
  TGraphErrors *g1b_e = new TGraphErrors(0);
  TGraphErrors *g1c_e = new TGraphErrors(0);
  TGraphErrors *g1d_e = new TGraphErrors(0);
  TGraphErrors *g1e_e = new TGraphErrors(0);
  TGraph *g1a_pl = new TGraph(0);
  TGraph *g1a_mn = new TGraph(0);
  TGraph *g1b_pl = new TGraph(0);
  TGraph *g1b_mn = new TGraph(0);
  TGraph *g1c_pl = new TGraph(0);
  TGraph *g1c_mn = new TGraph(0);
  TGraph *g1d_pl = new TGraph(0);
  TGraph *g1d_mn = new TGraph(0);
  TGraph *g1e_pl = new TGraph(0);
  TGraph *g1e_mn = new TGraph(0);

  TGraphErrors *g2a_e = new TGraphErrors(0);
  TGraphErrors *g2b_e = new TGraphErrors(0);
  TGraphErrors *g2c_e = new TGraphErrors(0);
  TGraphErrors *g2d_e = new TGraphErrors(0);
  TGraphErrors *g2e_e = new TGraphErrors(0);
  TGraph *g2a_pl = new TGraph(0);
  TGraph *g2a_mn = new TGraph(0);
  TGraph *g2b_pl = new TGraph(0);
  TGraph *g2b_mn = new TGraph(0);
  TGraph *g2c_pl = new TGraph(0);
  TGraph *g2c_mn = new TGraph(0);
  TGraph *g2d_pl = new TGraph(0);
  TGraph *g2d_mn = new TGraph(0);
  TGraph *g2e_pl = new TGraph(0);
  TGraph *g2e_mn = new TGraph(0);

  TGraphErrors *g3a_e = new TGraphErrors(0);
  TGraphErrors *g3b_e = new TGraphErrors(0);
  TGraphErrors *g3c_e = new TGraphErrors(0);
  TGraphErrors *g3d_e = new TGraphErrors(0);
  TGraphErrors *g3e_e = new TGraphErrors(0);
  TGraph *g3a_pl = new TGraph(0);
  TGraph *g3a_mn = new TGraph(0);
  TGraph *g3b_pl = new TGraph(0);
  TGraph *g3b_mn = new TGraph(0);
  TGraph *g3c_pl = new TGraph(0);
  TGraph *g3c_mn = new TGraph(0);
  TGraph *g3d_pl = new TGraph(0);
  TGraph *g3d_mn = new TGraph(0);
  TGraph *g3e_pl = new TGraph(0);
  TGraph *g3e_mn = new TGraph(0);

  
  // Different projections in one loop
  /*
  for (int icase = 0; icase != 3; ++icase) {

    TGraph *g1 = new TGraph(0);
    TGraph *g2 = new TGraph(0);
    TGraph *g21 = new TGraph(0);

    TGraphErrors *g1_e = new TGraphErrors(0);
    TGraph *g1_pl = new TGraph(0);
    TGraph *g1_mn = new TGraph(0);
    TGraphErrors *g2_e = new TGraphErrors(0);
    TGraph *g2_pl = new TGraph(0);
    TGraph *g2_mn = new TGraph(0);

    for (int i = 1; i != h->GetNbinsX()+1; ++i) {
      double eta = h->GetBinCenter(i);
    } // for i

    TGraphErrors *g1_e = new TGraphErrors(0);
    TGraph *g1_pl = new TGraph(0);
    TGraph *g1_mn = new TGraph(0);
    TGraphErrors *g2_e = new TGraphErrors(0);
    TGraph *g2_pl = new TGraph(0);
    TGraph *g2_mn = new TGraph(0);
  } // icase
  */

  const int npt = 6;
  //double ptbins[npt] = {30, 40, 50, 80, 140,500};
  double ptbins[npt] = {30, 60, 120, 240, 480, 960};
  TGraphErrors *g21s[npt];
  for (int i = 0; i != npt; ++i) {
    g21s[i] = new TGraphErrors(0);
  }
  
  for (int i = 1; i != h->GetNbinsX()+1; ++i) {

    double eta = h->GetBinCenter(i);
    if (fabs(eta)>4.7) continue;

    // ***** Pt = 30, 50, 80, 120, 200, 500 *****
    {
      for (int j = 0; j != npt; ++j) {

	TGraphErrors *g21 = g21s[j];
	double pt = ptbins[j];
	double energy = pt*cosh(eta);

	if (energy < 6500.) {
	  // Asymmetric corrections now
	  double y1 = 0.5*(getEtaPtE(JEC1, +eta, pt, energy)
			    + getEtaPtE(JEC1, -eta, pt, energy));
	  double y2 = 0.5*(getEtaPtE(JEC2, +eta, pt, energy)
			    + getEtaPtE(JEC2, -eta, pt, energy));

	  g21->SetPoint(g21->GetN(), eta, y2/y1);
	} // energy < 6500
      } // for j
    } // pt bins
 
    // ***** Pt = 30 
    {
      double pt = 30.;
      double energy = pt*cosh(eta);
      
      if (energy < 6500.) {
	// Asymmetric corrections now
	double y1 = 0.5*(getEtaPtE(JEC1, +eta, pt, energy)
			  + getEtaPtE(JEC1, -eta, pt, energy));
	double y2 = 0.5*(getEtaPtE(JEC2, +eta, pt, energy)
			  + getEtaPtE(JEC2, -eta, pt, energy));
	double y3(0);
	if (dothree) y3 = 0.5*(getEtaPtE(JEC3, +eta, pt, energy)
				+ getEtaPtE(JEC3, -eta, pt, energy));
	// negative side
	if (_usenegative) {
	  y1 = getEtaPtE(JEC1, -eta, pt, energy);
	  y2 = getEtaPtE(JEC2, -eta, pt, energy);
	  y3 = (dothree ? getEtaPtE(JEC3, -eta, pt, energy) : 0);
	}
	// positive side
	if (_usepositive) {
	  y1 = getEtaPtE(JEC1, +eta, pt, energy);
	  y2 = getEtaPtE(JEC2, +eta, pt, energy);
	  y3 = (dothree ? getEtaPtE(JEC3, +eta, pt, energy) : 0);
	}
	double e1 = getEtaPtUncert(jecUnc1, JEC1, eta, pt);
	double e2 = getEtaPtUncert(jecUnc2, JEC2, eta, pt);
	double e3 = (dothree ? getEtaPtUncert(jecUnc3, JEC3, eta, pt) : 0);
	
	g1a->SetPoint(g1a->GetN(), eta, y1);
	g2a->SetPoint(g2a->GetN(), eta, y2);
	g3a->SetPoint(g3a->GetN(), eta, y3);
	g21a->SetPoint(g21a->GetN(),eta, y2/y1);
	//
	g1a_pl->SetPoint(g1a_pl->GetN(), eta, y1*(1+e1));
	g1a_mn->SetPoint(g1a_mn->GetN(), eta, y1*(1-e1));
	g1a_e->SetPoint(i-1, eta, y1);
	g1a_e->SetPointError(i-1, 0., y1*e1);
	//
	g2a_pl->SetPoint(g2a_pl->GetN(), eta, y2*(1+e2));
	g2a_mn->SetPoint(g2a_mn->GetN(), eta, y2*(1-e2));
	g2a_e->SetPoint(i-1, eta, y2);
	g2a_e->SetPointError(i-1, 0., y2*e2);
	//
	g3a_pl->SetPoint(g3a_pl->GetN(), eta, y3*(1+e3));
	g3a_mn->SetPoint(g3a_mn->GetN(), eta, y3*(1-e3));
	g3a_e->SetPoint(i-1, eta, y3);
	g3a_e->SetPointError(i-1, 0., y3*e3);
      }
    }

    // ***** Pt = 100 
    {
      double pt = 100.;
      double energy = pt*cosh(eta);
      
      if (energy < 6500.) {
	// Asymmetric corrections now
	double y1 = 0.5*(getEtaPtE(JEC1, +eta, pt, energy)
			 + getEtaPtE(JEC1, -eta, pt, energy));
	double y2 = 0.5*(getEtaPtE(JEC2, +eta, pt, energy)
			 + getEtaPtE(JEC2, -eta, pt, energy));
	double y3(0);
	if (dothree) y3 = 0.5*(getEtaPtE(JEC3, +eta, pt, energy)
			       + getEtaPtE(JEC3, -eta, pt, energy));
	// negative side
	if (_usenegative) {
	  y1 = getEtaPtE(JEC1, -eta, pt, energy);
	  y2 = getEtaPtE(JEC2, -eta, pt, energy);
	  y3 = (dothree ? getEtaPtE(JEC3, -eta, pt, energy) : 0);
	}
	// positive side
	if (_usepositive) {
	  y1 = getEtaPtE(JEC1, +eta, pt, energy);
	  y2 = getEtaPtE(JEC2, +eta, pt, energy);
	  y3 = (dothree ? getEtaPtE(JEC3, +eta, pt, energy) : 0);
	}
	double e1 = getEtaPtUncert(jecUnc1, JEC1, eta, pt);
	double e2 = getEtaPtUncert(jecUnc2, JEC2, eta, pt);
	double e3 = (dothree ? getEtaPtUncert(jecUnc3, JEC3, eta, pt) : 0);
	
	g1b->SetPoint(g1b->GetN(), eta, y1);
	g2b->SetPoint(g2b->GetN(), eta, y2);
	g3b->SetPoint(g3b->GetN(), eta, y3);
	g21b->SetPoint(g21b->GetN(),eta, y2/y1);
	//
	g1b_pl->SetPoint(g1b_pl->GetN(), eta, y1*(1+e1));
	g1b_mn->SetPoint(g1b_mn->GetN(), eta, y1*(1-e1));
	g1b_e->SetPoint(i-1, eta, y1);
	g1b_e->SetPointError(i-1, 0., y1*e1);
	//
	g2b_pl->SetPoint(g2b_pl->GetN(), eta, y2*(1+e2));
	g2b_mn->SetPoint(g2b_mn->GetN(), eta, y2*(1-e2));
	g2b_e->SetPoint(i-1, eta, y2);
	g2b_e->SetPointError(i-1, 0., y2*e2);
	//
	g3b_pl->SetPoint(g3b_pl->GetN(), eta, y3*(1+e3));
	g3b_mn->SetPoint(g3b_mn->GetN(), eta, y3*(1-e3));
	g3b_e->SetPoint(i-1, eta, y3);
	g3b_e->SetPointError(i-1, 0., y3*e3);
      }
    }

    // ***** Pt = 1000
    // ***** Pt = 600 (extends better up to interesting |eta|~3)
    {
      //double pt = 1000.;
      double pt = 600.;
      double energy = pt*cosh(eta);
      
      if (energy < 6500.) {
	// Asymmetric corrections now
	double y1 = 0.5*(getEtaPtE(JEC1, +eta, pt, energy)
			 + getEtaPtE(JEC1, -eta, pt, energy));
	double y2 = 0.5*(getEtaPtE(JEC2, +eta, pt, energy)
			 + getEtaPtE(JEC2, -eta, pt, energy));
	double y3(0);
	if (dothree) y3 = 0.5*(getEtaPtE(JEC3, +eta, pt, energy)
			       + getEtaPtE(JEC3, -eta, pt, energy));
	// negative side
	if (_usenegative) {
	  y1 = getEtaPtE(JEC1, -eta, pt, energy);
	  y2 = getEtaPtE(JEC2, -eta, pt, energy);
	  y3 = (dothree ? getEtaPtE(JEC3, -eta, pt, energy) : 0);
	}
	// positive side
	if (_usepositive) {
	  y1 = getEtaPtE(JEC1, +eta, pt, energy);
	  y2 = getEtaPtE(JEC2, +eta, pt, energy);
	  y3 = (dothree ? getEtaPtE(JEC3, +eta, pt, energy) : 0);
	}
	double e1 = getEtaPtUncert(jecUnc1, JEC1, eta, pt);
	double e2 = getEtaPtUncert(jecUnc2, JEC2, eta, pt);
	double e3 = (dothree ? getEtaPtUncert(jecUnc3, JEC3, eta, pt) : 0);
	
	g1c->SetPoint(g1c->GetN(), eta, y1);
	g2c->SetPoint(g2c->GetN(), eta, y2);
	g3c->SetPoint(g3c->GetN(), eta, y3);
	g21c->SetPoint(g21c->GetN(),eta, y2/y1);
	//
	g1c_pl->SetPoint(g1c_pl->GetN(), eta, y1*(1+e1));
	g1c_mn->SetPoint(g1c_mn->GetN(), eta, y1*(1-e1));
	g1c_e->SetPoint(i-1, eta, y1);
	g1c_e->SetPointError(i-1, 0., y1*e1);
	//
	g2c_pl->SetPoint(g2c_pl->GetN(), eta, y2*(1+e2));
	g2c_mn->SetPoint(g2c_mn->GetN(), eta, y2*(1-e2));
	g2c_e->SetPoint(i-1, eta, y2);
	g2c_e->SetPointError(i-1, 0., y2*e2);
	//
	g3c_pl->SetPoint(g3c_pl->GetN(), eta, y3*(1+e3));
	g3c_mn->SetPoint(g3c_mn->GetN(), eta, y3*(1-e3));
	g3c_e->SetPoint(i-1, eta, y3);
	g3c_e->SetPointError(i-1, 0., y3*e3);
      }
    }

    // ***** E = 1000 
    {
      double energy = 1000.;
      double pt = energy/cosh(eta);
     
      if (pt > 10.) {
	// Asymmetric corrections now
	double y1 = 0.5*(getEtaPtE(JEC1, +eta, pt, energy)
			  + getEtaPtE(JEC1, -eta, pt, energy));
	double y2 = 0.5*(getEtaPtE(JEC2, +eta, pt, energy)
			  + getEtaPtE(JEC2, -eta, pt, energy));
	double y3(0);
	if (dothree) y3 = 0.5*(getEtaPtE(JEC3, +eta, pt, energy)
				+ getEtaPtE(JEC3, -eta, pt, energy));
	// negative side
	if (_usenegative) {
	  y1 = getEtaPtE(JEC1, -eta, pt, energy);
	  y2 = getEtaPtE(JEC2, -eta, pt, energy);
	  y3 = (dothree ? getEtaPtE(JEC3, -eta, pt, energy) : 0);
	}
	// positive side
	if (_usepositive) {
	  y1 = getEtaPtE(JEC1, +eta, pt, energy);
	  y2 = getEtaPtE(JEC2, +eta, pt, energy);
	  y3 = (dothree ? getEtaPtE(JEC3, +eta, pt, energy) : 0);
	}
	double e1 = getEtaPtUncert(jecUnc1, JEC1, eta, pt);
	double e2 = getEtaPtUncert(jecUnc2, JEC2, eta, pt);
	double e3 = (dothree ? getEtaPtUncert(jecUnc3, JEC3, eta, pt) : 0);
	
	g1e->SetPoint(g1e->GetN(), eta, y1);
	g2e->SetPoint(g2e->GetN(), eta, y2);
	g3e->SetPoint(g3e->GetN(), eta, y3);
	g21e->SetPoint(g21e->GetN(),eta, y2/y1);
	//
	g1e_pl->SetPoint(g1e_pl->GetN(), eta, y1*(1+e1));
	g1e_mn->SetPoint(g1e_mn->GetN(), eta, y1*(1-e1));
	g1e_e->SetPoint(i-1, eta, y1);
	g1e_e->SetPointError(i-1, 0., y1*e1);
	//
	g2e_pl->SetPoint(g2e_pl->GetN(), eta, y2*(1+e2));
	g2e_mn->SetPoint(g2e_mn->GetN(), eta, y2*(1-e2));
	g2e_e->SetPoint(i-1, eta, y2);
	g2e_e->SetPointError(i-1, 0., y2*e2);
	//
	g3e_pl->SetPoint(g3e_pl->GetN(), eta, y3*(1+e3));
	g3e_mn->SetPoint(g3e_mn->GetN(), eta, y3*(1-e3));
	g3e_e->SetPoint(i-1, eta, y3);
	g3e_e->SetPointError(i-1, 0., y3*e3);
      }
    }
  } // for i

  for (int i = 1; i != hpt->GetNbinsX()+1; ++i) {

    // ***** Eta = 0  => |eta|<1.3
    const double etas[] = {0, 0.261, 0.522, 0.783, 1.044, 1.305};
    const int neta = sizeof(etas)/sizeof(etas[0])-1;
    double sumy1(0), sumy2(0), sumy3(0), sume1(0), sume2(0), sume3(0);
    int nsum(0);

    double pt = hpt->GetBinCenter(i);
    for (int j = 0; j != neta; ++j) {
      
      //double eta = 0.;
      double eta = 0.5*(etas[j]+etas[j+1]);
      double energy = pt*cosh(eta);
      
      if (pt>10 && energy < 6500.) {
	// Asymmetric corrections now
	double y1 = 0.5*(getEtaPtE(JEC1, +eta, pt, energy)
			  + getEtaPtE(JEC1, -eta, pt, energy));
	double y2 = 0.5*(getEtaPtE(JEC2, +eta, pt, energy)
			  + getEtaPtE(JEC2, -eta, pt, energy));
	double y3(0);
	if (dothree) y3 = 0.5*(getEtaPtE(JEC3, +eta, pt, energy)
				+ getEtaPtE(JEC3, -eta, pt, energy));
	// negative side
	if (_usenegative) {
	  y1 = getEtaPtE(JEC1, -eta, pt, energy);
	  y2 = getEtaPtE(JEC2, -eta, pt, energy);
	  y3 = (dothree ? getEtaPtE(JEC3, -eta, pt, energy) : 0);
	}
	// positive side
	if (_usepositive) {
	  y1 = getEtaPtE(JEC1, +eta, pt, energy);
	  y2 = getEtaPtE(JEC2, +eta, pt, energy);
	  y3 = (dothree ? getEtaPtE(JEC3, +eta, pt, energy) : 0);
	}
	double e1 = getEtaPtUncert(jecUnc1, JEC1, eta, pt);
	double e2 = getEtaPtUncert(jecUnc2, JEC2, eta, pt);
	double e3 = (dothree ? getEtaPtUncert(jecUnc3, JEC3, eta, pt) : 0);
	
	sumy1 = (sumy1*nsum + y1) / (1.+nsum);
	sumy2 = (sumy2*nsum + y2) / (1.+nsum);
	sumy3 = (sumy3*nsum + y3) / (1.+nsum);
	sume1 = sqrt(sume1*sume1*nsum + e1*e1) / sqrt(1.+nsum);
	sume2 = sqrt(sume2*sume2*nsum + e2*e2) / sqrt(1.+nsum);
	sume3 = sqrt(sume3*sume3*nsum + e3*e3) / sqrt(1.+nsum);
	++nsum;
      } // pt>ptmin && e<emax
    } // for j

    g1d->SetPoint(g1d->GetN(), pt, sumy1);
    g2d->SetPoint(g2d->GetN(), pt, sumy2);
    g3d->SetPoint(g3d->GetN(), pt, sumy3);
    //
    g1d_pl->SetPoint(g1d_pl->GetN(), pt, sumy1*(1+sume1));
    g1d_mn->SetPoint(g1d_mn->GetN(), pt, sumy1*(1-sume1));
    g1d_e->SetPoint(i-1, pt, sumy1);
    g1d_e->SetPointError(i-1, 0., sumy1*sume1);
    //
    g2d_pl->SetPoint(g2d_pl->GetN(), pt, sumy2*(1+sume2));
    g2d_mn->SetPoint(g2d_mn->GetN(), pt, sumy2*(1-sume2));
    g2d_e->SetPoint(i-1, pt, sumy2);
    g2d_e->SetPointError(i-1, 0., sumy2*sume2);
    //
    g3d_pl->SetPoint(g3d_pl->GetN(), pt, sumy3*(1+sume3));
    g3d_mn->SetPoint(g3d_mn->GetN(), pt, sumy3*(1-sume3));
    g3d_e->SetPoint(i-1, pt, sumy3);
    g3d_e->SetPointError(i-1, 0., sumy3*sume3);
    //}
    //} // *** Eta = 0 => |eta|<1.3
  } // for i

  // Generic legend
  //TLegend *leg = new TLegend(0.20,dothree ? 0.70 : 0.75,0.40,0.85,"","brNDC");
  //leg->SetTextSize(0.045);
  //leg->SetBorderSize(0);
  //leg->SetFillStyle(kNone);
  //leg->AddEntry(g2a,s2,"LPF");
  //leg->AddEntry(g1a,s1,"LPF");
  //if (dothree) leg->AddEntry(g3a,s3,"LPF");
  //leg->Draw();

  // For legends
  g3a->SetFillStyle(3003);
  g3a->SetFillColor(kGreen+2);
  g3b->SetFillStyle(3003);
  g3b->SetFillColor(kGreen+2);
  g3c->SetFillStyle(3003);
  g3c->SetFillColor(kGreen+2);
  g3d->SetFillStyle(3003);
  g3d->SetFillColor(kGreen+2);
  g3e->SetFillStyle(3003);
  g3e->SetFillColor(kGreen+2);

  g1a->SetFillStyle(3003);
  g1a->SetFillColor(kBlue);
  g1b->SetFillStyle(3003);
  g1b->SetFillColor(kBlue);
  g1c->SetFillStyle(3003);
  g1c->SetFillColor(kBlue);
  g1d->SetFillStyle(3003);
  g1d->SetFillColor(kBlue);
  g1e->SetFillStyle(3003);
  g1e->SetFillColor(kBlue);

  g2a->SetFillStyle(3003);
  g2a->SetFillColor(kRed);
  g2b->SetFillStyle(3003);
  g2b->SetFillColor(kRed);
  g2c->SetFillStyle(3003);
  g2c->SetFillColor(kRed);
  g2d->SetFillStyle(3003);
  g2d->SetFillColor(kRed);
  g2e->SetFillStyle(3003);
  g2e->SetFillColor(kRed);
  
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);

  map<string,const char*> texmap;
  texmap["AK5PF"] = "R = 0.5, PF";
  texmap["AK5PFchs"] = "R = 0.5, PF+CHS";
  texmap["AK7PF"] = "R = 0.7, PF";
  texmap["AK7PFchs"] = "R = 0.7, PF+CHS";


  // ***** Pt = 30
  {
    c1a->cd();
    //h->DrawClone("AXIS");

    if (dothree) {
      g3a_e->SetFillStyle(3003);
      g3a_e->SetFillColor(kGreen+2);
      g3a_e->Draw("SAME E3");
      g3a_pl->SetLineColor(kGreen-9);
      g3a_pl->SetLineStyle(kSolid);//kDotted);
      g3a_pl->Draw("SAMEL");
      g3a_mn->SetLineColor(kGreen-9);
      g3a_mn->SetLineStyle(kSolid);//kDotted);
      g3a_mn->Draw("SAMEL");
    }

    g1a_e->SetFillStyle(3003);
    g1a_e->SetFillColor(kBlue);
    g1a_e->Draw("SAME E3");
    g1a_pl->SetLineColor(kBlue-9);
    g1a_pl->SetLineStyle(kSolid);//kDotted);
    g1a_pl->Draw("SAMEL");
    g1a_mn->SetLineColor(kBlue-9);
    g1a_mn->SetLineStyle(kSolid);//kDotted);
    g1a_mn->Draw("SAMEL");

    g2a_e->SetFillStyle(3003);
    g2a_e->SetFillColor(kRed);
    g2a_e->Draw("SAME E3");
    g2a_pl->SetLineColor(kRed-9);
    g2a_pl->SetLineStyle(kSolid);//kDotted);
    g2a_pl->Draw("SAMEL");
    g2a_mn->SetLineColor(kRed-9);
    g2a_mn->SetLineStyle(kSolid);//kDotted);
    g2a_mn->Draw("SAMEL");
        
    if (dothree) {
      g3a->SetMarkerStyle(kOpenSquare);
      g3a->SetMarkerColor(kGreen+2);
      g3a->SetLineColor(kGreen+2);
      g3a->Draw("SAMEPL");
    }

    g1a->SetMarkerStyle(kFullSquare);
    g1a->SetMarkerColor(kBlue);
    g1a->SetLineColor(kBlue);
    g1a->Draw("SAMEPL");

    g2a->SetMarkerStyle(kFullCircle);
    g2a->SetMarkerColor(kRed);
    g2a->SetLineColor(kRed);
    g2a->Draw("SAMEPL");

    //tex->DrawLatex(0.20,0.88,Form("p_{T,%s} = 30 GeV%s, %s",cgen,
    //				  l1||l2l3 ? Form(", #LT#mu#GT = %1.1f",_mu) : "",
    //				  cm));
    //tex->DrawLatex(0.65,0.80,a);
    //leg->Draw();

    tex->DrawLatex(0.19,0.81,Form("p_{T,%s} = 30 GeV",cgen));
    if (l1) tex->DrawLatex(0.19,0.74,Form("#LT#mu#GT = %1.1f",_mu));
    if (l1) tex->DrawLatex(0.19,0.68,Form("R = %1.1f, PF+chs",0.4));

    //TLegend *leg1a = tdrLeg(0.60, dothree ? 0.66 : 0.72, 0.90, 0.90);
    //TLegend *leg1a = tdrLeg(0.57, dothree ? 0.66 : 0.72, 0.87, 0.90);
    //TLegend *leg1a = tdrLeg(0.47, dothree ? 0.66 : 0.72, 0.87, 0.90);
    TLegend *leg1a = tdrLeg(0.47, dothree ? 0.72 : 0.78, 0.87, 0.93);
    leg1a->SetHeader(texmap[a]);
    leg1a->AddEntry(g2a,s2,"LPF");
    leg1a->AddEntry(g1a,s1,"LPF");
    if (dothree) leg1a->AddEntry(g3a,s3,"LPF");

    //if (!mc) cmsPrel(_lumi);
    //if (mc)  cmsPrel(0);
    gPad->RedrawAxis();
  }

  // ***** Pt = 100
  {
    c1b->cd();
    //h->DrawClone("AXIS");

    if (dothree) {
      g3b_e->SetFillStyle(3003);
      g3b_e->SetFillColor(kGreen+2);
      g3b_e->Draw("SAME E3");
      g3b_pl->SetLineColor(kGreen-9);
      g3b_pl->SetLineStyle(kSolid);//kDotted);
      g3b_pl->Draw("SAMEL");
      g3b_mn->SetLineColor(kGreen-9);
      g3b_mn->SetLineStyle(kSolid);//kDotted);
      g3b_mn->Draw("SAMEL");    
    }

    g1b_e->SetFillStyle(3003);
    g1b_e->SetFillColor(kBlue);
    g1b_e->Draw("SAME E3");
    g1b_pl->SetLineColor(kBlue-9);
    g1b_pl->SetLineStyle(kSolid);//kDotted);
    g1b_pl->Draw("SAMEL");
    g1b_mn->SetLineColor(kBlue-9);
    g1b_mn->SetLineStyle(kSolid);//kDotted);
    g1b_mn->Draw("SAMEL");    

    g2b_e->SetFillStyle(3003);
    g2b_e->SetFillColor(kRed);
    g2b_e->Draw("SAME E3");
    g2b_pl->SetLineColor(kRed-9);
    g2b_pl->SetLineStyle(kSolid);//kDotted);
    g2b_pl->Draw("SAMEL");
    g2b_mn->SetLineColor(kRed-9);
    g2b_mn->SetLineStyle(kSolid);//kDotted);
    g2b_mn->Draw("SAMEL");    

    if (dothree) {
      g3b->SetMarkerStyle(kOpenSquare);
      g3b->SetMarkerColor(kGreen+2);
      g3b->SetLineColor(kGreen+2);
      g3b->Draw("SAMEPL");
    }

    g1b->SetMarkerStyle(kFullSquare);
    g1b->SetMarkerColor(kBlue);
    g1b->SetLineColor(kBlue);
    g1b->Draw("SAMEPL");
    
    g2b->SetMarkerStyle(kFullCircle);
    g2b->SetMarkerColor(kRed);
    g2b->SetLineColor(kRed);
    g2b->Draw("SAMEPL");
    
    //tex->DrawLatex(0.20,0.88,Form("p_{T,%s} = 100 GeV%s, %s",cgen,
    //				  l1||l2l3 ? Form(", #LT#mu#GT = %1.1f",_mu) : "",
    //				  cm));
    //tex->DrawLatex(0.65,0.80,a);
    //leg->Draw();

    tex->DrawLatex(0.19,0.75,Form("p_{T,%s} = 100 GeV",cgen));
    if (l1) tex->DrawLatex(0.19,0.68,Form("#LT#mu#GT = %1.1f",_mu));

    //TLegend *leg1b = tdrLeg(0.60, dothree ? 0.66 : 0.73, 0.90, 0.90);
    //TLegend *leg1b = tdrLeg(0.57, dothree ? 0.66 : 0.73, 0.87, 0.90);
    TLegend *leg1b = tdrLeg(0.47, dothree ? 0.66 : 0.73, 0.87, 0.90);
    leg1b->SetHeader(texmap[a]);
    leg1b->AddEntry(g2b,s2,"LPF");
    leg1b->AddEntry(g1b,s1,"LPF");
    if (dothree) leg1b->AddEntry(g3b,s3,"LPF");

    //if (!mc) cmsPrel(_lumi);
    //if (mc)  cmsPrel(0);
    gPad->RedrawAxis();
  }

  // ***** Pt = 1000
  {
    c1c->cd();

    if (dothree) {
      g3c_e->SetFillStyle(3003);
      g3c_e->SetFillColor(kGreen+2);
      g3c_e->Draw("SAME E3");
      g3c_pl->SetLineColor(kGreen-9);
      g3c_pl->SetLineStyle(kSolid);
      g3c_pl->Draw("SAMEL");
      g3c_mn->SetLineColor(kGreen-9);
      g3c_mn->SetLineStyle(kSolid);
      g3c_mn->Draw("SAMEL");    
    }

    g1c_e->SetFillStyle(3003);
    g1c_e->SetFillColor(kBlue);
    g1c_e->Draw("SAME E3");
    g1c_pl->SetLineColor(kBlue-9);
    g1c_pl->SetLineStyle(kSolid);
    g1c_pl->Draw("SAMEL");
    g1c_mn->SetLineColor(kBlue-9);
    g1c_mn->SetLineStyle(kSolid);
    g1c_mn->Draw("SAMEL");    

    g2c_e->SetFillStyle(3003);
    g2c_e->SetFillColor(kRed);
    g2c_e->Draw("SAME E3");
    g2c_pl->SetLineColor(kRed-9);
    g2c_pl->SetLineStyle(kSolid);
    g2c_pl->Draw("SAMEL");
    g2c_mn->SetLineColor(kRed-9);
    g2c_mn->SetLineStyle(kSolid);
    g2c_mn->Draw("SAMEL");    

    if (dothree) {
      g3c->SetMarkerStyle(kOpenSquare);
      g3c->SetMarkerColor(kGreen+2);
      g3c->SetLineColor(kGreen+2);
      g3c->Draw("SAMEPL");
    }

    g1c->SetMarkerStyle(kFullSquare);
    g1c->SetMarkerColor(kBlue);
    g1c->SetLineColor(kBlue);
    g1c->Draw("SAMEPL");
    
    g2c->SetMarkerStyle(kFullCircle);
    g2c->SetMarkerColor(kRed);
    g2c->SetLineColor(kRed);
    g2c->Draw("SAMEPL");
    
    //tex->DrawLatex(0.19,0.75,Form("p_{T,%s} = 1000 GeV",cgen));
    tex->DrawLatex(0.19,0.75,Form("p_{T,%s} = 600 GeV",cgen));
    if (l1) tex->DrawLatex(0.19,0.68,Form("#LT#mu#GT = %1.1f",_mu));

    //TLegend *leg1c = tdrLeg(0.60, dothree ? 0.66 : 0.73, 0.90, 0.90);
    //TLegend *leg1c = tdrLeg(0.57, dothree ? 0.66 : 0.73, 0.87, 0.90);
    TLegend *leg1c = tdrLeg(0.47, dothree ? 0.66 : 0.73, 0.87, 0.90);
    leg1c->SetHeader(texmap[a]);
    leg1c->AddEntry(g2c,s2,"LPF");
    leg1c->AddEntry(g1c,s1,"LPF");
    if (dothree) leg1c->AddEntry(g3c,s3,"LPF");

    gPad->RedrawAxis();
  }

  // ***** E = 1000
  {
    c1e->cd();
    //h->DrawClone("AXIS");

    if (dothree) {
      g3e_e->SetFillStyle(3003);
      g3e_e->SetFillColor(kGreen+2);
      g3e_e->Draw("SAME E3");
      g3e_pl->SetLineColor(kGreen-9);
      g3e_pl->SetLineStyle(kSolid);//kDotted);
      g3e_pl->Draw("SAMEL");
      g3e_mn->SetLineColor(kGreen-9);
      g3e_mn->SetLineStyle(kSolid);//kDotted);
      g3e_mn->Draw("SAMEL");
    }

    g1e_e->SetFillStyle(3003);
    g1e_e->SetFillColor(kBlue);
    g1e_e->Draw("SAME E3");
    g1e_pl->SetLineColor(kBlue-9);
    g1e_pl->SetLineStyle(kSolid);//kDotted);
    g1e_pl->Draw("SAMEL");
    g1e_mn->SetLineColor(kBlue-9);
    g1e_mn->SetLineStyle(kSolid);//kDotted);
    g1e_mn->Draw("SAMEL");

    g2e_e->SetFillStyle(3003);
    g2e_e->SetFillColor(kRed);
    g2e_e->Draw("SAME E3");
    g2e_pl->SetLineColor(kRed-9);
    g2e_pl->SetLineStyle(kSolid);//kDotted);
    g2e_pl->Draw("SAMEL");
    g2e_mn->SetLineColor(kRed-9);
    g2e_mn->SetLineStyle(kSolid);//kDotted);
    g2e_mn->Draw("SAMEL");
    
    if (dothree) {
      g3e->SetMarkerStyle(kOpenSquare);
      g3e->SetMarkerColor(kGreen+2);
      g3e->SetLineColor(kGreen+2);
      g3e->Draw("SAMEPL");
    }

    g1e->SetMarkerStyle(kFullSquare);
    g1e->SetMarkerColor(kBlue);
    g1e->SetLineColor(kBlue);
    g1e->Draw("SAMEPL");
    
    g2e->SetMarkerStyle(kFullCircle);
    g2e->SetMarkerColor(kRed);
    g2e->SetLineColor(kRed);
    g2e->Draw("SAMEPL");
    
    //tex->DrawLatex(0.20,0.88,Form("E_{%s} = 1000 GeV%s, %s",cgen,
    //				  l1||l2l3 ? Form(", #LT#mu#GT = %1.1f",_mu) : "",
    //				  cm));
    //tex->DrawLatex(0.65,0.80,a);
    //leg->Draw();

    tex->DrawLatex(0.19,0.75,Form("E_{%s} = 1000 GeV",cgen));
    if (l1) tex->DrawLatex(0.19,0.68,Form("#LT#mu#GT = %1.1f",_mu));

    //TLegend *leg1e = tdrLeg(0.60, dothree ? 0.66 : 0.73, 0.90, 0.90);
    //TLegend *leg1e = tdrLeg(0.57, dothree ? 0.66 : 0.73, 0.87, 0.90);
    TLegend *leg1e = tdrLeg(0.47, dothree ? 0.66 : 0.73, 0.87, 0.90);
    leg1e->SetHeader(texmap[a]);
    leg1e->AddEntry(g2e,s2,"LPF");
    leg1e->AddEntry(g1e,s1,"LPF");
    if (dothree) leg1e->AddEntry(g3e,s3,"LPF");

    //if (!mc) cmsPrel(_lumi);
    //if (mc)  cmsPrel(0);
    gPad->RedrawAxis();
  }

  // ***** Eta = 0
  {
    c1d->cd();
    c1d->SetLogx();
    
    if (dothree) {
      g3d_e->SetFillStyle(3003);
      g3d_e->SetFillColor(kGreen+2);
      g3d_e->Draw("SAME E3");
      g3d_pl->SetLineColor(kGreen-9);
      g3d_pl->SetLineStyle(kSolid);
      g3d_pl->Draw("SAMEL");
      g3d_mn->SetLineColor(kGreen-9);
      g3d_mn->SetLineStyle(kSolid);
      g3d_mn->Draw("SAMEL");
    }

    g1d_e->SetFillStyle(3003);
    g1d_e->SetFillColor(kBlue);
    g1d_e->Draw("SAME E3");
    g1d_pl->SetLineColor(kBlue-9);
    g1d_pl->SetLineStyle(kSolid);
    g1d_pl->Draw("SAMEL");
    g1d_mn->SetLineColor(kBlue-9);
    g1d_mn->SetLineStyle(kSolid);
    g1d_mn->Draw("SAMEL");

    g2d_e->SetFillStyle(3003);
    g2d_e->SetFillColor(kRed);
    g2d_e->Draw("SAME E3");
    g2d_pl->SetLineColor(kRed-9);
    g2d_pl->SetLineStyle(kSolid);
    g2d_pl->Draw("SAMEL");
    g2d_mn->SetLineColor(kRed-9);
    g2d_mn->SetLineStyle(kSolid);
    g2d_mn->Draw("SAMEL");
        
    if (dothree) {
      g3d->SetMarkerStyle(kOpenSquare);
      g3d->SetMarkerColor(kGreen+2);
      g3d->SetLineColor(kGreen+2);
      g3d->Draw("SAMEPL");
    }

    g1d->SetMarkerStyle(kFullSquare);
    g1d->SetMarkerColor(kBlue);
    g1d->SetLineColor(kBlue);
    g1d->Draw("SAMEPL");

    g2d->SetMarkerStyle(kFullCircle);
    g2d->SetMarkerColor(kRed);
    g2d->SetLineColor(kRed);
    g2d->Draw("SAMEPL");

    //tex->DrawLatex(0.19,0.75,"|#eta| = 0");
    tex->DrawLatex(0.19,0.75,"|#eta| < 1.3");
    if (l1) tex->DrawLatex(0.19,0.68,Form("#LT#mu#GT = %1.1f",_mu));

    //TLegend *leg1d = tdrLeg(0.60, dothree ? 0.66 : 0.72, 0.90, 0.90);
    //TLegend *leg1d = tdrLeg(0.57, dothree ? 0.66 : 0.72, 0.87, 0.90);
    //TLegend *leg1d = tdrLeg(0.57, dothree ? 0.71 : 0.77, 0.87, 0.95);
    TLegend *leg1d = tdrLeg(0.47, dothree ? 0.71 : 0.77, 0.87, 0.95);
    leg1d->SetHeader(texmap[a]);
    leg1d->AddEntry(g2d,s2,"LPF");
    leg1d->AddEntry(g1d,s1,"LPF");
    if (dothree) leg1d->AddEntry(g3d,s3,"LPF");

    gPad->RedrawAxis();
  }

  string ctype = string(cl1)+string(cl2l3)+string(cplus)+string(cres);
  const char *cs = ctype.c_str();
  if (_pdf) {
    c1a->SaveAs(Form("pdf/compareJECversions_%s_%s_%s_Pt030.pdf",a,cm,cs));
    c1b->SaveAs(Form("pdf/compareJECversions_%s_%s_%s_Pt100.pdf",a,cm,cs));
    //c1c->SaveAs(Form("pdf/compareJECversions_%s_%s_%s_Pt1000.pdf",a,cm,cs));
    c1c->SaveAs(Form("pdf/compareJECversions_%s_%s_%s_Pt600.pdf",a,cm,cs));
    c1d->SaveAs(Form("pdf/compareJECversions_%s_%s_%s_Eta00.pdf",a,cm,cs));
    c1e->SaveAs(Form("pdf/compareJECversions_%s_%s_%s_Q1000.pdf",a,cm,cs));
  }
  if (_C) {
    c1a->SaveAs(Form("pdfC/compareJECversions_%s_%s_%s_Pt030.C",a,cm,cs));
    c1b->SaveAs(Form("pdfC/compareJECversions_%s_%s_%s_Pt100.C",a,cm,cs));
    //c1c->SaveAs(Form("pdfC/compareJECversions_%s_%s_%s_Pt1000.C",a,cm,cs));
    c1c->SaveAs(Form("pdfC/compareJECversions_%s_%s_%s_Pt600.C",a,cm,cs));
    c1d->SaveAs(Form("pdfC/compareJECversions_%s_%s_%s_Eta00.C",a,cm,cs));
    c1e->SaveAs(Form("pdfC/compareJECversions_%s_%s_%s_Q1000.C",a,cm,cs));
  }

  // ***** Multiple pT bins for ratio only
  {
    int colors[] = {kBlack, kBlue, kCyan+2, kGreen+2, kOrange+2, kRed};
    int styles[] = {kSolid, kDashed, kDotted, kDashDotted, kDashed, kSolid};
    c0->cd();
    h->SetMinimum(0.85);//0.90);
    h->SetMaximum(1.25);//1.20);
    TH1D *h0 = (TH1D*)h->DrawClone("AXIS");
    h0->GetYaxis()->SetTitle(Form("%s%s%s%s (%s / %s)",cl1,cl2l3,cpl,cres,
				  s2s, s1s));

    TLine *l = new TLine();
    l->SetLineStyle(kDashed);
    l->DrawLine(0,1,5,1);
    l->DrawLine(0,1.05,5,1.05);
    l->SetLineStyle(kDotted);
    l->DrawLine(0,1.02,5,1.02);
    l->DrawLine(0,0.98,5,0.98);

    TLegend *leg = new TLegend(0.20,0.62,0.60,0.92,"","brNDC");
    leg->SetTextSize(0.045);
    leg->SetBorderSize(0);
    leg->SetFillStyle(kNone);
    leg->Draw();

    for (int i = 0; i != npt; ++i) {
      TGraphErrors *gr = g21s[i];
      gr->SetLineColor(colors[i]);//i+1);
      gr->SetLineStyle(styles[i]);//i+1);
      gr->SetLineWidth(3);
      gr->Draw("SAME L");

      //leg->AddEntry(gr,Form("%s / %s (p_{T,%s}=%1.0f GeV)",
      //s2s, s1s, cgen,ptbins[i]),"L");
      leg->AddEntry(gr,Form("p_{T,%s} = %1.0f GeV",
			    cgen,ptbins[i]),"L");
    } // for i
    tex->DrawLatex(0.70,0.85,a);
    //cmsPrel(_lumi);
  }
  c0->SaveAs(Form("pdf/compareJECversions_%s_%s_%s_%sover%s.pdf",a,cm,cs,s2s,s1s));

  
  {// Ratio plots
    c2->cd();

    h->SetMinimum(l1||l2l3 ? 0.80 : 0.95);
    h->SetMaximum(l1||l2l3 ? 1.20 : 1.12);
    if (res && (algo=="CALO" || algo=="JPT")) h->SetMinimum(0.80);
    h->DrawClone("AXIS");

    TLine *l = new TLine();
    l->SetLineStyle(kDashed);
    l->DrawLine(0,1,5,1);
    l->SetLineColor(kGreen+2);
    l->SetLineStyle(kDotted);
    l->DrawLine(0,1.02,5,1.02);
    l->DrawLine(0,0.98,5,0.98);
    l->SetLineColor(kBlue);
    l->SetLineStyle(kDotted);
    l->DrawLine(0,1.05,5,1.05);
    l->DrawLine(0,0.95,5,0.95);

    g21a->SetMarkerStyle(kOpenCircle);
    g21a->SetMarkerColor(kBlue);
    g21a->SetLineColor(kBlue);
    g21a->SetLineStyle(kDashed);
    g21a->SetLineWidth(3);
    g21a->Draw("SAMEL");

    g21c->SetMarkerStyle(kOpenSquare);
    g21c->SetMarkerColor(kRed);
    g21c->SetLineColor(kRed);
    g21c->SetLineStyle(kDotted);
    g21c->SetLineWidth(3);
    g21c->Draw("SAMEL");

    g21b->SetMarkerStyle(kFullCircle);
    g21b->SetMarkerColor(kGreen+2);
    g21b->SetLineColor(kGreen+2);
    g21b->SetLineStyle(kSolid);
    g21b->SetLineWidth(3);
    g21b->Draw("SAMEL");
    
    TLegend *leg = new TLegend(0.20,0.77,0.60,0.92,"","brNDC");
    leg->SetTextSize(0.045);
    leg->SetBorderSize(0);
    leg->SetFillStyle(kNone);
    leg->AddEntry(g21a,Form("%s / %s (p_{T,%s}=30 GeV)",s2s,s1s,cgen),"LP");
    leg->AddEntry(g21b,Form("%s / %s (p_{T,%s}=100 GeV)",s2s,s1s,cgen),"LP");
    leg->AddEntry(g21c,Form("%s / %s (p_{T,%s}=600 GeV)",s2s,s1s,cgen),"LP");
    //leg->AddEntry(g21c,Form("%s / %s (p_{T,%s}=1000 GeV)",s2s,s1s,cgen),"LP");
    //leg->AddEntry(g21c,Form("%s / %s (E_{%s}=1000 GeV)",s2s,s1s,cgen),"LP");
    leg->Draw();

    //if (!mc) cmsPrel(_lumi);
    //if (mc)  cmsPrel(0);
    gPad->RedrawAxis();
    
    if(_pdf) c2->SaveAs(Form("pdf/compareJECversions_%s_%s_%s_Ratios.pdf",a,cm,cs));
  } // Ratio plots
} // compareJECversions


