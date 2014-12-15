#include "JECUncertainty.hpp"
#include "TMath.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"
#include "Math/BrentRootFinder.h"
//#include "Math/RootFinderAlgorithms.h"

#include <cmath>
#include <map>

using namespace std;

// Turn debug mode on if the code fails with an exception from JEC packages
// This is most likely a missing/misnamed file given to JetCorrectorParameters
// The last file printed out before the crash in debug mode is usually the fault
bool debug = false;


JECUncertainty::JECUncertainty(const jec::JetAlgo& algo, 
			       const jec::DataType& type, 
			       const jec::ErrorTypes& errType,
			       const double mu) :
  _algo(algo), _type(type), _errType(errType), _mu(mu)
{

  _fjes = 0; _emat = 0; _fhb = _fl1 = 0;
  _fl1ref = _fl1up = _fl1dw = 0;

  _algo = algo;
  _calo = (_algo==jec::AK5CALO || _algo==jec::AK7CALO);
  _jpt = (_algo==jec::AK5JPT);
  _pfchs = (_algo==jec::AK5PFchs || _algo==jec::AK7PFchs);
  _pflow = (_algo==jec::AK5PF || _algo==jec::AK7PF || _pfchs);
  _ideal = false;//(_algo==IDEAL);
  _trkbase = (_pflow || _jpt);

  _ajet = 1;
  if (algo==jec::AK7PF || algo==jec::AK7PFchs || algo==jec::AK7CALO)
    _ajet = pow(0.7/0.5,2);

  _InitL1();
  _InitJEC();
  _InitL2Res();
  _InitL3Res();


}


// Uncertainty is returned as _relative_ uncertainty
// Keep systematics signed for correlations
double JECUncertainty::Uncert(const double pTprime, const double eta) {

  double err2 = 0;
  double eta2 = min(max(eta,-5.190),5.190); // fix for drawing macro

  // 2013 systematics
  double errAbs(0), errRel(0), errPileUp(0), errFlavor(0), errTime(0);
  if (_errType & jec::kAbsolute) {
    errAbs = _Absolute(pTprime);
    err2 += errAbs * errAbs;
  }
  if (_errType & jec::kRelative) {
    errRel = _Relative(pTprime, eta2);
    err2 += errRel * errRel;
  }
  if (_errType & (jec::kPileUp | jec::kPileUpMuZero | jec::kPileUpEnvelope)) {
    errPileUp = _PileUp(pTprime, eta2);
    err2 += errPileUp * errPileUp;
  }
  if (_errType & jec::kFlavorMask) {
    errFlavor = _Flavor(pTprime, eta);
    err2 += errFlavor*errFlavor;
  }
  if (_errType & (jec::kTime | jec::kTimePtMask)) {
    errTime = _Time(pTprime, eta2);
    err2 += errTime * errTime;
  }
  
  double err = sqrt(err2);

  // if requesting single source, return signed for sign-changing cases
  if (!(_errType & ~jec::kAbsoluteFrag)) return errAbs;
  if (!(_errType & ~jec::kAbsoluteSPRE)) return errAbs;
  if (!(_errType & ~jec::kAbsoluteSPRH)) return errAbs;
  if (!(_errType & ~jec::kRelativePtBB))  return errRel;
  if (!(_errType & ~jec::kRelativePtEC1)) return errRel;
  if (!(_errType & ~jec::kRelativePtEC2)) return errRel;
  if (!(_errType & ~jec::kRelativePtHF))  return errRel;
  if (!(_errType & ~jec::kRelativePt))    return errRel; // EXTRA
  if (!(_errType & ~jec::kPileUpDataMC)) return errPileUp;
  if (!(_errType & ~jec::kPileUpPtBB))   return errPileUp;
  if (!(_errType & ~jec::kPileUpPtEC1))  return errPileUp;
  if (!(_errType & ~jec::kPileUpPtEC2))  return errPileUp;
  if (!(_errType & ~jec::kPileUpPtHF))   return errPileUp;
  if (!(_errType & ~jec::kPileUpPt))     return errPileUp; // EXTRA
  if (!(_errType & ~jec::kPileUpMuZero)) return errPileUp; // OPT
  if (!(_errType & ~jec::kTimePtRunA)) return errTime;
  if (!(_errType & ~jec::kTimePtRunB)) return errTime;
  if (!(_errType & ~jec::kTimePtRunC)) return errTime;
  if (!(_errType & ~jec::kTimePtRunD)) return errTime;
  //
  if (!(_errType & ~jec::kFlavorMask))   return errFlavor;

  return err;
} // Uncert

void JECUncertainty::_InitL1() {

  // RandomCone (V0) files from Ia Iashvili by e-mail (DropBox link)
  // On 17 May 2014, at 16:05
  // Re: Summer14 combination files with RD MC

  // L1FastJet (V1) files for data from Ia Iashvili by e-mail (tar file)
  // On 24 Mar 2014, at 16:52
  // Re: L2 Residuals Corrections - this time packed and ready

  // DataMcSF files from Ia Iasvili by e-mail (DropBox link)
  // On 21 May 2014, at 16:15
  // Re: New MC-based L1
  // https://www.dropbox.com/s/dm1ndhdmx9zdr72/DataMcSF.tar.gz

  map<jec::JetAlgo, const char*> names;
  names[jec::AK5PF] = "AK5PF";
  names[jec::AK5PFchs] = "AK5PFchs";
  names[jec::AK5CALO] = "AK5Calo";
  names[jec::AK7PF] = "AK7PF";
  names[jec::AK7PFchs] = "AK7PFchs";
  names[jec::AK7CALO] = "AK7Calo";
  const char *a = names[_algo];
  string directory = "CondFormats/JetMETObjects/data/";
  const char *d = directory.c_str();

  // For PileUpPt in DATA
  {
    const char *s = Form("%sWinter14_V0_DATA_L1FastJetPU_%s_pt.txt",d,a);
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l1);
    _jecL1DTflat = new FactorizedJetCorrector(v);
  }
  {
    const char *s = Form("%sWinter14_V0_MC_L1FastJetPU_%s_pt.txt",d,a);
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l1);
    _jecL1MCflat = new FactorizedJetCorrector(v);
  }
  // For PileUpPt in MC
  {
    const char *s = Form("%sWinter14_V1_DATA_L1FastJet_%s.txt",d,a);
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l1);
    _jecL1DTpt = new FactorizedJetCorrector(v);
  }
  {
    const char *s = Form("%sWinter14_V1_MC_L1FastJet_%s.txt",d,a);
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l1);
    _jecL1MCpt = new FactorizedJetCorrector(v);
  }
  // For PileUpDataMC
  {
    const char *s = Form("%sWinter14_DataMcSF_L1FastJetPU_%s.txt",d,a);
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l1);
    _jecL1sf = new FactorizedJetCorrector(v);
  }

  // For PileUpPtRef
  {
    const char *a = "AK5PFchs";
    const char *s = Form("%sWinter14_V0_DATA_L1FastJetPU_%s_pt.txt",d,a);
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l1);
    _jecL1DTflat_ak5pfchs = new FactorizedJetCorrector(v);
  }
  {
    const char *a = "AK5PFchs";
    const char *s = Form("%sWinter14_V1_DATA_L1FastJet_%s.txt",d,a);
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l1);
    _jecL1DTpt_ak5pfchs = new FactorizedJetCorrector(v);
  }

} // InitL1


void JECUncertainty::_InitJEC() {

  // JEC files collected by Alexx Perloff to web directory
  // http://people.physics.tamu.edu/aperloff/CMS_JEC/index.php?path=Winter_14%2FWinter14_V4_txts/
  // V5 has L2L3Residuals changed on top of these (V4 below is really new V5)

  map<jec::JetAlgo, const char*> names;
  names[jec::AK5PF] = "AK5PF";
  names[jec::AK5PFchs] = "AK5PFchs";
  names[jec::AK5CALO] = "AK5Calo";
  names[jec::AK7PF] = "AK7PF";
  names[jec::AK7PFchs] = "AK7PFchs";
  names[jec::AK7CALO] = "AK7Calo";  
  const char *a = names[_algo];
  string directory = "CondFormats/JetMETObjects/data/";
  const char *d = directory.c_str();

  const char *s;
  s = Form("%sWinter14_V4_DATA_L1FastJet_%s.txt",d,a);
  if (debug) cout << s << endl << flush;
  JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
  s = Form("%sWinter14_V4_DATA_L2Relative_%s.txt",d,a);
  if (debug) cout << s << endl << flush;
  JetCorrectorParameters *l2 = new JetCorrectorParameters(s);
  s = Form("%sWinter14_V4_DATA_L3Absolute_%s.txt",d,a);
  if (debug) cout << s << endl << flush;
  JetCorrectorParameters *l3 = new JetCorrectorParameters(s);
  // Only one L3Residual derived for now (although we will later clone this)
  s = Form("%sWinter14_V4_DATA_L2L3Residual_AK5PFchs.txt",d);
  if (debug) cout << s << endl << flush;
  JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);

  vector<JetCorrectorParameters> v;
  v.push_back(*l1);
  v.push_back(*l2);
  v.push_back(*l3);
  v.push_back(*l2l3res);
  _jecDefault = new FactorizedJetCorrector(v);
  _jec = _jecDefault;

  // Another version using Random Cone offset (L1 V0)
  s = Form("%sWinter14_V0_DATA_L1FastJetPU_%s_pt.txt",d,a);
  if (debug) cout << s << endl << flush;
  JetCorrectorParameters *l1v0 = new JetCorrectorParameters(s);

  vector<JetCorrectorParameters> v0;
  v0.push_back(*l1v0);
  v0.push_back(*l2);
  v0.push_back(*l3);
  v0.push_back(*l2l3res);
  _jecWithL1V0 = new FactorizedJetCorrector(v0);

} // InitJEC

void JECUncertainty::_InitL2Res() {

  // LOGLIN/FLAT + JERup/JERdown + STAT for AK5PFchs, AK5PF and AK7PF
  // On 26 Sep 2014, at 16:25, Rathjens, Denis
  // Re: winter14 correction update
  // => Winter14_V5_uncertaintyFilesL2res.tar.gz (directory 'summary')
  // Fixed empty trailing lines and wrong column numbers in the files by hand

  map<jec::JetAlgo, const char*> names;
  names[jec::AK5PF] = "AK5PF";
  names[jec::AK5PFchs] = "AK5PFchs";
  names[jec::AK5CALO] = "AK5PF"; // Replace "AK5Calo";
  names[jec::AK7PF] = "AK7PF";
  names[jec::AK7PFchs] = "AK7PF"; // Replace "AK7PFchs";
  names[jec::AK7CALO] = "AK7PF"; // Replace "AK5Calo";
  string directory = "CondFormats/JetMETObjects/data/";
  const char *d = directory.c_str();
  const char *a = names[_algo];

  const char *s;
  // For RelativePt
  {
    s = Form("%sWinter14_V5_DATA_L2L3Residual_%s.txt.FLAT",d,a);
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l2l3res);
    _jecL2ResFlat = new FactorizedJetCorrector(v);
  }
  {
    s = Form("%sWinter14_V5_DATA_L2L3Residual_%s.txt.LOGLIN",d,a);
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l2l3res);
    _jecL2ResPt = new FactorizedJetCorrector(v);
  }
  // For RelativeJER
  {
    s = Form("%sWinter14_V5_DATA_L2L3Residual_%s.txt.JERup",d,a);
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l2l3res);
    _jecL2jerup = new FactorizedJetCorrector(v);
  }
  {
    s = Form("%sWinter14_V5_DATA_L2L3Residual_%s.txt.JERdown",d,a);
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l2l3res);
    _jecL2jerdw = new FactorizedJetCorrector(v);
  }
  // For RelativeStat
  {
    s = Form("%sWinter14_V5_DATA_L2L3Residual_%s.txt.STAT",d,a);
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l2l3res);
    _jecL2stat = new FactorizedJetCorrector(v);
  }
  
} // InitL2Res

void JECUncertainty::_InitL3Res() {

  // V3PT
  const int n = 2;
  const double pars[n] =
    {0.9773, -0.0442};
  const double emata[n][n] =
    {{3.508e-06,  1.033e-07},
     {1.033e-07,   0.000231}};

  // Sub-optimal to copy every time, but this is the most flexible interface
  if (!_emat) {
    _emat = new TMatrixD(n, n);
    for (int i = 0; i != n; ++i) {
      for (int j = 0; j != n; ++j) {
	(*_emat)[i][j] = emata[i][j];
      }
    }
  }
  if (!_fjes) {
    _fjes = new TF1("fjes",_jesfit,10.,4000.,n);
    for (int i = 0; i != n; ++i) {
      _fjes->SetParameter(i, pars[i]);
    }
  }

} // _InitL3Res


// Solve pTraw from pTprime = pTraw / R(pTraw) using Brent's method
// We want to provide JEC uncertainties vs pTprime, not pTraw, but JEC
// is only available as a function of pTraw
double JECUncertainty::_Rjet(double pTprime, double eta,
			     double ajet = -1, double mu = -1,
			     FactorizedJetCorrector *jec = 0) {
  
  if (ajet<0) ajet = 0.785*_ajet;
  if (mu<0) mu = 11.85; 
  double npv = _NpvFromMu(mu);
  double rho = _RhoFromMu(mu);
  if (!jec) jec = _jec;

  ResponseFunc f(pTprime,jec,npv,eta,rho,ajet);
  
  //ROOT::Math::Roots::Brent brf;
  ROOT::Math::BrentRootFinder brf;
  //brf.SetLogScan(true);
  // Set parameters of the method
  brf.SetFunction(f,std::max(2.0,0.25*pTprime),std::min(4*pTprime,3500.0));
  bool found_root = brf.Solve(50,1e-4,1e-5);
  double pTraw = brf.Root();
  double rjet = pTraw /pTprime;
  jec->setJetPt(pTraw);
  jec->setJetEta(eta);
  jec->setJetA(ajet);
  jec->setRho(rho);
  double corr = jec->getCorrection();
  if(std::abs((pTraw * corr - pTprime)/pTprime) > 0.001) {
    std::cout << "NPV:" << npv << "   Brent: status:" << brf.Status() << " " << brf.Iterations() << '\n';
    std::cout << "R: pTprime:" << pTprime << "  eta:" << eta << " pTraw:" << pTraw << " corr:" << corr 
              << " pTcor:" << corr * pTraw << " rjet*pTprime:" << rjet * pTprime << '\n'; 
  }
  assert(std::abs((pTraw * corr - pTprime)/pTprime) < 0.001);  
  assert(found_root);
  assert(brf.Status() == 0);
  /*
  if((eta > 2.6) && (eta < 2.8)) {
    _jec->setJetE(pTprime*cosh(eta)); 
    _jec->setJetEta(eta); _jec->setNPV(_npv);
    _jec->setJetPt(pTprime);
    std::cout << "NPV:" << _npv << " cor:" << _jec->getCorrection() << "   Brent: status:" << brf.Status() << " " << brf.Iterations() << '\n';
    std::cout << "R: pTprime:" << pTprime << "  eta:" << eta << " pTraw:" << pTraw << " corr:" << corr 
              << " pTcor:" << corr * pTraw << " rjet*pTprime:" << _rjet * pTprime << '\n'; 
  }
  */
  return rjet;
} // _Rjet


// Combine pT dependent absolute scale uncertainties
double JECUncertainty::_Absolute(const double pt) const {

  double stat          = (_errType & jec::kAbsoluteStat          ? _AbsoluteStat(pt)          : 0.);
  double scale         = (_errType & jec::kAbsoluteScale         ? _AbsoluteScale()         : 0.);
  double FlMap         = (_errType & jec::kAbsoluteFlavorMapping ? _AbsoluteFlavorMapping() : 0.); //for backward-compatibility/historical reasons
  double MPFBias       = (_errType & jec::kAbsoluteMPFBias       ? _AbsoluteMPFBias()       : 0.);
  double spr           = (_errType & jec::kAbsoluteSPR           ? _AbsoluteSPR(pt)         : 0.);
  double frag          = (_errType & jec::kAbsoluteFrag          ? _AbsoluteFrag(pt)        : 0.);

  // signed sources
  if (!(_errType & ~jec::kAbsoluteStat)) return stat;
  if (!(_errType & ~jec::kAbsoluteSPRE)) return spr;
  if (!(_errType & ~jec::kAbsoluteSPRH)) return spr;
  if (!(_errType & ~jec::kAbsoluteFrag)) return frag;

  return sqrt(stat*stat + scale*scale + FlMap*FlMap + MPFBias*MPFBias + spr*spr + frag*frag);
} // Absolute


// pT dependent fit of L3Res
//TF1 *fhb(0), *fl1(0);
//Double_t JECUncertainty::_jesfit(Double_t *x, Double_t *p) {
Double_t _jesfit(Double_t *x, Double_t *p) {

  double pt = x[0];

  // 2012 RDMC
  if (!_fhb) _fhb = new TF1("fhb","max(0,[0]+[1]*pow(x,[2]))",10,3500);
  _fhb->SetParameters(1.03091e+00, -5.11540e-02, -1.54227e-01); // SPRH
  if (!_fl1) _fl1 = new TF1("fl1","1+([0]+[1]*log(x))/x",30,2000);
  _fl1->SetParameters(-2.36997, 0.413917);
  // p[0]: overall scale shift, p[1]: HCAL shift in % (full band +3%)
  // p[2]: fraction of PileUpPtBB uncertainty
  double ptref = 208;//225.;
  double jes =  (p[0] + p[1]/3.*100*(_fhb->Eval(pt)-_fhb->Eval(ptref))
		 //+ p[2]*(fl1->Eval(pt)-fl1->Eval(ptref)));
		 + -0.090*(_fl1->Eval(pt)-_fl1->Eval(ptref)));

  return jes;
} // jesfit
// Make sure this matches above, it is used in TimePt
//Double_t JECUncertainty::_jeshb(double pt, double hb) {
Double_t _jeshb(double pt, double hb) {

  if (!_fhb) _fhb = new TF1("fhb","max(0,[0]+[1]*pow(x,[2]))",10,3500);
  _fhb->SetParameters(1.03091e+00, -5.11540e-02, -1.54227e-01); // SPRH
  double hb0 = -0.0442; // V3PT
  double jes = hb/3.*100*(_fhb->Eval(pt)-1) - hb0/3.*100*(_fhb->Eval(pt)-1);

  return jes;
}

// Fit uncertainty
double JECUncertainty::_jesfitunc(double x, TF1 *f, TMatrixD *emat) const {

  assert(f);
  assert(emat);
  int n = f->GetNpar();
  vector<double> df(n);
  for (int i = 0; i != n; ++i) {
    Double_t p = f->GetParameter(i);
    Double_t dp = 0.1*sqrt((*emat)[i][i]);
    f->SetParameter(i, p + dp);
    double fup = f->Eval(x);
    f->SetParameter(i, p - dp);
    double fdw = f->Eval(x);
    f->SetParameter(i, p);
    df[i] = (dp ? (fup - fdw) / (2. * dp) : 0);
  }
  double sumerr2 = 0;
  for (int i = 0; i != n; ++i) {
    for (int j = 0; j != n; ++j) {
      sumerr2 += (*emat)[i][j] * df[i] * df[j];
    }
  }
  
  return sqrt(sumerr2);
}

// Continuation of 2011 ATLAS/CMS JES correlation discussions
// Absolute scale split-up: Statistics
double JECUncertainty::_AbsoluteStat(double pTprime) const {

  // Uncertainty band from the pT-dependent global fit with SPRH
  // Fit includes Zmm, Zee and gamma scales as floating nuisance parameters,
  // so subtract fit scale uncertainty in quadrature to avoid double-counting
  double AbsStat = _jesfitunc(pTprime, _fjes, _emat);
  double AbsScale = _AbsoluteScale();
  double AbsStatSys = sqrt(max(AbsStat*AbsStat - AbsScale*AbsScale, 0.));

  // Should we add sign to the "statistical" uncerainty source
  // (this is really the fit uncertainty for the SPRH component)

  return AbsStatSys;
}


// Continuation of 2011 ATLAS/CMS JES correlation discussions
// Absolute scale split-up: Scale
double JECUncertainty::_AbsoluteScale() const {

  // - Global fit is given Zmm, Zee and gamma scale uncertainties, as well
  //   as extrapolation systematics for MPF (was mpfbias) and pT balance
  // - Extract the constant part from the error matrix as AbsoluteScale. This
  //   is achieved by shifting the reference pT=208, which decouples the
  //   constant term from slope parameters (correlations falls from ~90% to 9%)
  double AbsScaleSys = 0.0019; // for pTref=208 GeV

  return AbsScaleSys;
}


// Continuation of 2011 ATLAS/CMS JES correlation discussions
// Absolute scale split-up: MPF bias
double JECUncertainty::_AbsoluteMPFBias() const {

  // - global fit is done to MPF and pT balance simultaneously so previous
  //   mpfbias uncertainty from FSR/ISR is included in fit uncertainty => drop
  // - beambias and neutrinos shared by MPF and pT balance, and not part
  //   of global fit => keep (but update estimates?)

  // Estimate beam bias from Pythia/Herwig difference
  // of genMET with (a) all particles, (b) all particles at |eta|<5,
  // (c) no neutrinos and no |eta|>5 particles
  // (not done properly, these are just rough estimates)
  // Also, would need to propagate the pT dependence in the global fit
  double beambias = 0.002; // (a)-(b)

  // Primary estimate of neutrino bias is excess in PF electron and
  // muon fraction in dijet composition studies, which show 0.1% extra
  // in both, flat in pT. These should be associated with an equal
  // amount of neutrino energy on average, i.e. 0.2% total.
  // In addition to this, there is evidence of large mismodeling of
  // heavy flavor production from gluon splitting, let's say 100% for
  // gluon-split fraction of 50% in MC, so 50% overall.
  // MC has about 4% b and 10% c, with 25% of each decaying semileptonically
  // and semileptonic decays carrying about 12% energy in neutrinos for b's
  // (how much for c's?). Therefore
  // 4% b's (x) 50% GS (x) 25% of semileptonic (x) 12% of nu_E = 0.06%
  // 10% c's (x) 50% GS (x) 25% semileptonic (x) 6%(??) of E_nu = 0.15%
  // => 0.21% obtained this way is very compatible
  double neutrinos = 0.002; // (b)-(c)
  
  double err2 = beambias*beambias + neutrinos*neutrinos;

  double AbsMPFBiasSys = sqrt(err2);

  return AbsMPFBiasSys;
}

// Continuation of 2011 ATLAS/CMS JES correlation discussions
// Absolute scale split-up: Flavor mapping
double JECUncertainty::_AbsoluteFlavorMapping() const {

  // Flavor uncertainties now fully propagate each flavor through
  // Z/gamma+jet and dijet balancing steps so this original
  // flavor mapping uncertainty becomes zero

  double AbsFlavorMappingSys = 0.; // for "20% glue" reference point at pT=200 GeV

  return AbsFlavorMappingSys;
}


// High pT systematics from Herwig/Pythia ratio extrapolation
double JECUncertainty::_AbsoluteFrag(const double pTprime) const {

  double dr(0);
  if (_errType & jec::kAbsoluteFrag) {

    TF1 *f1 = new TF1("f1","[0]+[1]*log10(0.01*x)+[2]/x+[3]/(x*x)",10,3000);
    if (_calo) {      
      // drawFragFlavor[4]
      double p[4] = {1.0098, -0.0073, 1.788, -19.95}; // CALO
      f1->SetParameters(p[0],p[1],p[2],p[3]);
    }
    if (_pflow) {
      double etaref = 0;
      return _FlavorMixed(pTprime, etaref, "20% glue");
    }

    double r = f1->Eval(pTprime);
    double pTref = 208; // 2012 RDMC V3PT: effective <pT> for global fit
    double rref = f1->Eval(pTref);
    dr = r/rref-1;
    delete f1;
  }

  if (_ideal) dr *= 0.5; // Use Pythia/Herwig mean for extrapolation
  
  return dr;
} // AbsoluteFrag


// Single pion response uncertainty from propagating +/-3% SPR into JEC
// using FastSim. Results have been fitted with power law functions
double JECUncertainty::_AbsoluteSPR(const double pTprime) const {

  // New fits from Juska on Nov 26, 2012, 5:26 pm
  double errSPR(0), errSPRE(0), errSPRH(0);
  double difSPR(0), difSPRE(0), difSPRH(0);
  double refSPR(0), refSPRE(0), refSPRH(0);
  double refpt = 208.; // 2012 RDMC V3PT

  TF1 *f = new TF1("fmore","max(0,[0]+[1]*pow(x,[2]))",10,3500);

  //if (_errType & jec::kAbsoluteSPR) { // this is now obsolete
  //f->SetParameters(1.02829e+00, -6.22540e-02, -2.67123e-01);
  //if (_calo) f->SetParameters(1.02084e+00, 3.23027e-02, -9.46202e-01);
  //errSPR = f->Eval(pTprime);
  //refSPR = f->Eval(refpt);
  //}

  // Fix 2014-05-21: multiply errSPRE and errSRPH residuals by sqrt(2)
  // to keep errSPRE_3%(oplus)errSPRH_3% ~ errSPR_3%
  if (_errType & jec::kAbsoluteSPRE) { // SPR in ECAL
    f->SetParameters(1.00567e+00, -3.04275e-02, -6.75493e-01);
    if (_calo) f->SetParameters(1.00166e+00, 1.57065e-02, -2.06585e-01);
    difSPRE = sqrt(2.)*(f->Eval(pTprime)-1) + 1;
    refSPRE = sqrt(2.)*(f->Eval(refpt)-1) + 1;
  }
  if (_errType & jec::kAbsoluteSPRH) { // SPR in HCAL
    f->SetParameters(1.03091e+00, -5.11540e-02, -1.54227e-01);
    if (_calo) f->SetParameters(1.02246e+00, -1.55689e-02, -1.17219e-01);
    //errSPRH = sqrt(2.)*(f->Eval(pTprime)-1) + 1;
    //refSPRH = sqrt(2.)*(f->Eval(refpt)-1) + 1;
    // - SPRH is obtained from the global fit as -0.0442 +/- 0.0152,
    //   thus uncertainty can be scaled from 3% down to 1.52%
    // - 4.4% from fit is ok, since a-priori should have had sqrt(2)*3%*=4.2%
    difSPRH = 0.0152/0.03 * (f->Eval(pTprime)-1) + 1; // V3PT
    refSPRH = 0.0152/0.03 * (f->Eval(refpt)-1) + 1; // V3PT
  }

  // replace directly done SPR with pieces broken up into ECAL and HCAL
  errSPRE = (difSPRE-refSPRE);
  errSPRH = (difSPRH-refSPRH);
  errSPR = sqrt(pow(errSPRE,2) + pow(errSPRH,2));

  // signed sources
  if (!(_errType & ~jec::kAbsoluteSPRE)) return errSPRE;
  if (!(_errType & ~jec::kAbsoluteSPRH)) return errSPRH;

  return errSPR;
} // AbsoluteSPR


// Combine relative uncertainties
double JECUncertainty::_Relative(const double pTprime,
				 const double eta) const {
  
  double sjer = (_errType & jec::kRelativeJER ? _RelativeJER(pTprime, eta) : 0.);
  double sfsr = (_errType & jec::kRelativeFSR ? _RelativeFSR(eta) : 0.);
  double stat = (_errType & jec::kRelativeStat ? _RelativeStat(pTprime, eta) : 0.);
  double spt = (_errType & jec::kRelativePt ? _RelativePt(pTprime, eta) : 0.);

  // signed sources
  if (!(_errType & ~jec::kRelativePtBB)) return spt;
  if (!(_errType & ~jec::kRelativePtEC1)) return spt;
  if (!(_errType & ~jec::kRelativePtEC2)) return spt;
  if (!(_errType & ~jec::kRelativePtHF))  return spt;
  if (!(_errType & ~jec::kRelativePt))  return spt; // XTRA

  return sqrt(sjer*sjer + sfsr*sfsr + stat*stat + spt*spt);
}


// Relative scale uncertainty vs eta from JER bias
double JECUncertainty::_RelativeJER(const double pTprime,
				    const double eta) const {

  _jecL2jerup->setJetEta(eta);
  _jecL2jerup->setJetPt(pTprime);
  double up = _jecL2jerup->getCorrection();
  _jecL2jerdw->setJetEta(eta);
  _jecL2jerdw->setJetPt(pTprime);
  double dw = _jecL2jerdw->getCorrection();

  // Uncertainty is half of up and down, so full difference to mean
  double err = 0.5 * (up - dw);

  double x = fabs(eta);
  if (x<1.5) return 0; // Assuming BB negligible for now
  if ( (x>=1.5 && x<2.5 && (_errType & jec::kRelativeJEREC1)) ||
       (x>=2.5 && x<3.0 && (_errType & jec::kRelativeJEREC2)) ||
       (x>=3.0 && x<5.2 && (_errType & jec::kRelativeJERHF)) )
    return err;

  return 0;
} // RelativeJER


// Relative scale uncertainty vs eta from soft radiation (FSR)
double JECUncertainty::_RelativeFSR(const double eta) const {

  // New estimate based on MPF (SJ, DJ) closure plots from Henning by e-mail
  // Subject: 	Re: updated systematics for 53X
  // Date: 	April 25, 2013 1:12:53 PM GMT+03:00
  
  TF1 f("f","[0]*(pow(cosh(min(3.2,x)),2)-1)",0,5.2);
  f.SetParameter(0,1.5/(pow(cosh(3.2),2)-1));

  return 0.01*f.Eval(fabs(eta));
} // RelativeFSR


// Statistical uncertainty in L2Res (symmetrized, wide bins)
double JECUncertainty::_RelativeStat(const double pTprime,
				     const double eta) const {

  _jecL2stat->setJetEta(eta);
  _jecL2stat->setJetPt(pTprime); // pT doesn't matter
  double err = _jecL2stat->getCorrection();

  double x = fabs(eta);
  if ( (x>=2.5 && x<3.0 && (_errType & jec::kRelativeStatEC2)) ||
       (x>=3.0 && x<5.5 && (_errType & jec::kRelativeStatHF )) )
    return err;

  return 0;
} // RelativeStat

// Uncertainty in L2Res pT dependence: log-linear fit vs constant fit
// To-do: solve pTprime with Brent's method to compare at same pTprime?
double JECUncertainty::_RelativePt(const double pTprime,
				   const double eta) const {

  // limit pt to accessible range
  const double ptmin = 10.;
  const double emax = 4000;
  double pt = max(ptmin, min(pTprime, emax/cosh(eta)));

  _jecL2ResFlat->setJetPt(pt);
  _jecL2ResFlat->setJetEta(eta);
  double corrflat = _jecL2ResFlat->getCorrection();
  double ptraw =  pt / corrflat; // close enough in first approx
  _jecL2ResPt->setJetPt(ptraw);
  _jecL2ResPt->setJetEta(eta);
  double corrpt = _jecL2ResPt->getCorrection();

  // Use 50% of the slope as uncertainty consistently everywhere
  // We do correct for it for a reason, so 100% seems too conservative
  double kfactor = 0.5;
  double err = kfactor * (corrflat / corrpt - 1); 

  double x = fabs(eta);
  if ((x>=0.0 && x<1.3 && _errType & jec::kRelativePtBB) ||
      (x>=1.3 && x<2.5 && _errType & jec::kRelativePtEC1) ||
      (x>=2.5 && x<3.0 && _errType & jec::kRelativePtEC2) ||
      (x>=3.0 && x<5.5 && _errType & jec::kRelativePtHF))
    return err;

  return 0;
} // RelativePt


// Combine pileup uncertainty sources
double JECUncertainty::_PileUp(const double pTprime, const double eta) {
  
  //double pT = _Rjet(pTprime, eta) * pTprime;
  //if((pT < 10) && (eta > 2.5)) std::cout <<  "_PileUp: Pt: prime: " << pTprime << ", " << pT << '\n';

  double smc =   (_errType & jec::kPileUpDataMC ? _PileUpDataMC(pTprime, eta) : 0);
  double spt =   (_errType & (jec::kPileUpPt | jec::kPileUpMuZero) ?
		  _PileUpPt(pTprime, eta) : 0.);
  double senv = (_errType & jec::kPileUpEnvelope ?
		 _PileUpEnvelope(pTprime, eta) : 0);
  double err = sqrt(smc*smc + spt*spt + senv*senv);

  // signed sources
  if (!(_errType & ~jec::kPileUpDataMC)) return smc;
  if (!(_errType & ~jec::kPileUpPtRef)) return spt;
  if (!(_errType & ~jec::kPileUpPtBB)) return spt;
  if (!(_errType & ~jec::kPileUpPtEC1)) return spt;
  if (!(_errType & ~jec::kPileUpPtEC2)) return spt;
  if (!(_errType & ~jec::kPileUpPtHF)) return spt;
  if (!(_errType & ~jec::kPileUpPtEta)) return spt; // EXTRA
  if (!(_errType & ~jec::kPileUpMuZero)) return spt; // OPTIONAL
  if (!(_errType & ~jec::kPileUpEnvelope)) return senv; // OPTIONAL

  return err;
} // PileUp


// Pile-up uncertainty from data / MC difference
double JECUncertainty::_PileUpDataMC(const double pTprime, const double eta) {
  
  // Winter14 Data/MC uncertainty from scale factor variation vs rho
  // https://indico.cern.ch/event/308741/contribution/4/material/slides/1.pdf
  // pages 7-8 (JERC talk by Ia Iashvili)
  double rhomin = 7;
  double rhomax = 14;
  double rhoavg = 10;
  double sfmin = _L1SF(pTprime, eta, rhomin);
  double sfmax = _L1SF(pTprime, eta, rhomax);
  double sfavg = _L1SF(pTprime, eta, rhoavg);
  
  double pTraw = _Rjet(pTprime, eta) * pTprime;
  double l1 = _L1Data(pTraw, eta);
  
  double sys = fabs(l1-1) * max(fabs(sfmin-sfavg), fabs(sfmax-sfavg));
			     
  return sys;
} // PileUpDataMC


// Pile-up uncertainty from pT dependence
// Implemented as difference between MC truth (V1) and Random Cone (V0),
// with log-linear pT dependence corrected for as with L2Residuals
// BB uncertainty spans the whole detector due to dijet balance
double JECUncertainty::_PileUpPt(const double pTprime, const double eta) {

  // Limit eta to [-5,5] because V0 files don't go further out
  // Even closer in, because bias goes nuts in the last bin out
  double maxeta = 4.5;
  double etax = max(-maxeta,min(maxeta,eta));
  double x = fabs(etax);

  FactorizedJetCorrector *_l1flat = (_errType==jec::kData ?
				     _jecL1DTflat : _jecL1MCflat);
  FactorizedJetCorrector *_l1pt = (_errType==jec::kData ?
				   _jecL1DTpt : _jecL1MCpt);

  double sysref(0), syseta(0), syszero(0), sys(0);
  // Absolute scale offset fit with log(x)/x for AK5PFchs, which gives the
  // reference scale uncertainty. Assumption is that this uncertainty is only
  // coming from the neutral part of the jet core so that the uncertainty is
  // the same for all the algorithms and cone sizes (or should we scale?)

  // Only do this once since it's very time-consuming
  // Reference L1 residual is from AK5PFchs L3Residual fit
  if (!_fl1ref) {
    
    FactorizedJetCorrector *_l1flatref = _jecL1DTflat_ak5pfchs;
    FactorizedJetCorrector *_l1ptref = _jecL1DTpt_ak5pfchs;

    // Shape effectively used in L3Residual fit for low pT
    // ([2] was fixed in the end to a reasonable value, but shape was tested)
    _fl1ref = new TF1("fl1ref","[0]*([1]+[2]*log(x))/x",30,1000);
    _fl1ref->SetParameters(1,-2.36997, 0.413917); 
    _fl1ref->FixParameter(1,-2.36997); // shape fixed in global fit
    _fl1ref->FixParameter(2, 0.413917); // shape fixed in global fit
    _fl1up = (TF1*)_fl1ref->Clone("fl1up");
    _fl1dw = (TF1*)_fl1ref->Clone("fl1dw");
    
    const double x_pt[] =
      {28, 32, 37, 43, 49, 56, 64, 74, 84,
       97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 362, 430,
       507, 592, 686, 790, 905, 1032};
    const int ndiv_pt = sizeof(x_pt)/sizeof(x_pt[0])-1;
    
    double x_eta[] = {-1.4,-1.2,-1.0, -0.8,-0.6,-0.4,-0.2,0.,
		      0, 0.2,0.4,0.6,0.8,1.0, 1.2,1.4};
    const int ndiv_eta = sizeof(x_eta)/sizeof(x_eta[0])-1;
    
    // Average offset uncertainty over barrel (|eta|<1.3),
    // then fit with log(x)/x to extra effective reference uncertainty
    // (since L3Residual was studied including log(x)/x in the fit)
    // Difference to log(x) will be absorbed to PileUpPtBB (syseta) 
    //
    // What fraction of the full Flat-Pt difference do we use?
    // Global fit gave -35 +/- 25% for the PileUpPtBB as fit parameter
    // Use max of 10%(algo), 60% (algo) vs 35%(ak5pfchs)) as systematic
    // Difference between algo and AK5PFchs should result in larger systematic
    // than 25% variation for AK5PFchs itself
    const double k0(0.35), kdw(0.10), kup(0.60);
    TGraph *gr = new TGraph(0);
    TGraph *gup = new TGraph(0);
    TGraph *gdw = new TGraph(0);
    for (int ipt = 0; ipt != ndiv_pt; ++ipt) {
      
      double pt = 0.5*(x_pt[ipt] + x_pt[ipt+1]);
      
      double sumw(0), sumsysr(0), sumsysup(0), sumsysdw(0);
      for (int ieta = 0; ieta != ndiv_eta; ++ieta) {
	
	double etab = 0.5*(x_eta[ieta]+x_eta[ieta+1]);
	double l1fr = _Rjet(pt, etab, -1, -1, _l1flatref);
	double l1pr = _Rjet(pt, etab, -1, -1, _l1ptref);
	double sysr = k0 * (l1fr / l1pr - 1);
	sumsysr += sysr;
	sumw    += 1;
	double l1f = _Rjet(pt, etab, -1, -1, _l1flat);
	double l1p = _Rjet(pt, etab, -1, -1, _l1pt);
	double sysup = kup * (l1f / l1p - 1);
	sumsysup += sysup;
	double sysdw = kdw * (l1f / l1p - 1);
	sumsysdw += sysdw;
      } // for ieta
      double sysr = sumsysr / sumw;
      gr->SetPoint(ipt, pt, sysr);
      double sysup = sumsysup / sumw;
      gup->SetPoint(ipt, pt, sysup);
      double sysdw = sumsysdw / sumw;
      gdw->SetPoint(ipt, pt, sysdw);
    } // for ipt
    
    gr->Fit(_fl1ref, "QRN");
    gup->Fit(_fl1up, "QRN");
    gdw->Fit(_fl1dw, "QRN");
    delete gr;
    delete gup;
    delete gdw;
  } // !_fl1ref

  if (_errType & jec::kPileUpPtRef) {
    sysref = absmax(_fl1up->Eval(pTprime) - _fl1ref->Eval(pTprime),
		    _fl1dw->Eval(pTprime) - _fl1ref->Eval(pTprime));
    //delete f1;
  } // sysref

  // Residual relative PU offset after applying dijet balance comes from
  // different shapes 'a+(b+c*log(x))/x' for L1 and 'a+b*log(x)' for L2Residual
  // Also subtract the residual offset from barrel, assuming the sign
  // and magnitude are correlated between different eta regions
  // (use +60% consistently for variation and barrel average)
  if ( (x>=0.0 && x<1.3 && _errType & jec::kPileUpPtBB) ||
       (x>=1.3 && x<2.5 && _errType & jec::kPileUpPtEC1) ||
       (x>=2.5 && x<3.0 && _errType & jec::kPileUpPtEC2) ||
       (x>=3.0 && x<5.5 && _errType & jec::kPileUpPtHF) ||
       _errType & jec::kPileUpMuZero ) {

    // kfactor gives the maximal size of the effect,
    // which we estimate from the global fit range [10%,60%]
    // by just taking 60% (assuming this gives largest variation)
    // Even if this is not true for every single case, keep 60%
    // to have the signs of the variations correlated
    double kfactor = 1;
    if (x<1.3)           kfactor = 0.60;
    if (x>=1.3 && x<2.5) kfactor = 0.60;
    if (x>=2.5 && x<3.0) kfactor = 0.60;
    if (x>=3.0)          kfactor = 0.60;

    const double x_pt[] =
      {28, 32, 37, 43, 49, 56, 64, 74, 84,
       97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 362, 430,
       507, 592, 686, 790, 905, 1032};
    const int ndiv_pt = sizeof(x_pt)/sizeof(x_pt[0])-1;

    assert(_fl1ref);
    TGraph *g = new TGraph(0);
    for (int ipt = 0; ipt != ndiv_pt; ++ipt) {

      double pt = 0.5*(x_pt[ipt] + x_pt[ipt+1]);
      double l1f = _Rjet(pt, etax, -1, -1, _l1flat);
      double l1p = _Rjet(pt, etax, -1, -1, _l1pt);
      //double sysb = kfactor * (l1f / l1p - 1) - _fl1ref->Eval(pTprime);
      double sysb = kfactor * (l1f / l1p - 1) - _fl1up->Eval(pTprime);
      g->SetPoint(ipt, pt, sysb);
    } // for ipt
    
    // Shape used in L2Residual fit
    TF1 *f1 = new TF1("f1", "[0]+[1]*log(x)",
		      _ajet>1 ? 71 : 60, 2000./cosh(etax));
    // AK5PFchs was used for L3Residual fit so effectively goes down to 30 GeV
    if (_algo==jec::AK5PFchs) f1->SetRange(30, 2000./cosh(etax));
    g->Fit(f1,"QRN");

    if ( (x>=0.0 && x<1.3 && _errType & jec::kPileUpPtBB) ||
	 (x>=1.3 && x<2.5 && _errType & jec::kPileUpPtEC1) ||
	 (x>=2.5 && x<3.0 && _errType & jec::kPileUpPtEC2) ||
	 (x>=3.0 && x<5.5 && _errType & jec::kPileUpPtHF) ) {

      double l1f = _Rjet(pTprime, etax, -1, -1, _l1flat);
      double l1p = _Rjet(pTprime, etax, -1, -1, _l1pt);
      syseta = kfactor * (l1f / l1p - 1) - _fl1ref->Eval(pTprime)
	- f1->Eval(pTprime);
    }
    
    if (_errType & jec::kPileUpMuZero) {
      syszero = _fl1ref->Eval(pTprime) + f1->Eval(pTprime);
    }
    
    delete g;
    delete f1;
  } // syseta

  sys = sqrt(sysref*sysref + syseta*syseta + syszero*syszero);

  // For single sources
  if (!(_errType & ~jec::kPileUpPtRef)) return sysref;
  if (!(_errType & ~jec::kPileUpPt)) return syseta;
  if (!(_errType & ~jec::kPileUpPtBB)) return syseta;
  if (!(_errType & ~jec::kPileUpPtEC1)) return syseta;
  if (!(_errType & ~jec::kPileUpPtEC2)) return syseta;
  if (!(_errType & ~jec::kPileUpPtHF)) return syseta;
  if (!(_errType & ~jec::kPileUpMuZero)) return syszero;

  return sys;
} // _PileUpPt


// (Obsolete version of) Pile-up uncertainty from pT dependence
// Implemented as difference between MC truth (V1) and Random Cone (V0),
//double JECUncertainty::_PileUpPt(const double pTprime, const double eta) {
double JECUncertainty::_PileUpEnvelope(const double pTprime, const double eta) {

  // Limit eta to [-5,5] because V0 files don't go further out
  // Even closer in, because bias goes nuts in the last bin out
  double maxeta = 4.5;
  double etax = max(-maxeta,min(maxeta,eta));

  // - PileUpPtBB is included in the global fit, with result -0.3489 +/- 0.2550,
  //   thus the |eta|<1.3 uncertainty can be reduced to 25.5%
  // - It seems justified that the uncertainty for the other regions should be
  //   comparable to bias limit for BB, i.e. smaller than 35% + 25% = 60%
  // - In the future, one could also consider how dijet balance (+Z/gamma+jet)
  //   constrains PileUpPt uncertainty in other eta regions

  double kfactor = 0.60;//1;
  //double x = fabs(etax);
  //if (x<1.3)           kfactor = 0.25;//0.2550;
  //if (x>=1.3 && x<2.5) kfactor = 0.30;
  //if (x>=2.5 && x<3.0) kfactor = 0.60;
  //if (x>=3.0)          kfactor = 1.00;

  _jec = _jecDefault;
  double pT_p = _Rjet(pTprime, etax) * pTprime;  
  _jec = _jecWithL1V0;
  double pT_v0p = _Rjet(pTprime, etax) * pTprime;
  _jec = _jecDefault;
  //
  double l1n_p = _L1DataFlat(pT_v0p, etax); // V0Data
  double l1s_p = _L1Data(pT_p, etax); // V1Data = V1MC * V0Data / V0MC
  double sys_p = kfactor * (l1n_p / l1s_p - 1);
  //
  double l1nb_p = _L1MCFlat(pT_v0p, etax);
  double l1sb_p = _L1MC(pT_p, etax);
  double sysb_p = kfactor * (l1nb_p / l1sb_p - 1);

  _jec = _jecDefault;
  double pT_m = _Rjet(pTprime, -etax) * pTprime;
  _jec = _jecWithL1V0;
  double pT_v0m = _Rjet(pTprime, -etax) * pTprime;
  _jec = _jecDefault;    
  //
  double l1n_m = _L1DataFlat(pT_v0m, -etax);
  double l1s_m = _L1Data(pT_m, -etax);
  double sys_m = kfactor * (l1n_m / l1s_m - 1);
  //
  double l1nb_m = _L1MCFlat(pT_v0m, -etax);
  double l1sb_m = _L1MC(pT_m, -etax);
  double sysb_m = kfactor * (l1nb_m / l1sb_m - 1);
    
  // Symmetrize
  double sys = 0.5*(sys_p + sys_m);
  double sysb = 0.5*(sysb_p + sysb_m);

  if(sys > 1) {
    std::cout << "strange uncertainty value: << " << sys*_ajet << " _PileUpBias:" << pT_p << "," << eta << " l1:" << l1n_p << "," << l1s_p << '\n';
  }

  // Check that bias consistent for data and MC files
  if (fabs(sys-sysb)>0.001) {
    if (debug) {
      cout << "In PileUpBias(pTprime="<<pTprime<<",eta="<<eta<<"): ";
      cout << "sys="<<sys<<" sysb="<<sysb<<" diff="<<fabs(sys-sysb)<<endl<<flush;
    }
    // Fail if difference outside comfort zone;
    //assert(fabs(sys-sysb)<0.002);
  }

  //double x = fabs(etax);
  //if ((x>=0.0 && x<1.3 && _errType & jec::kPileUpPtBB) ||
  //  (x>=1.3 && x<2.5 && _errType & jec::kPileUpPtEC1) ||
  //  (x>=2.5 && x<3.0 && _errType & jec::kPileUpPtEC2) ||
  //  (x>=3.0 && x<5.5 && _errType & jec::kPileUpPtHF))
  return sys * _ajet;// _ajet: jet area/ 0.5^2
  
  //return 0;
} // _PileUpPt

// Combine jet flavor uncertainties
double JECUncertainty::_Flavor(double pTprime, double eta) const {

  double err2(0), errFlavor(0);
  if (_errType & jec::kFlavorQCD) {
    assert( !(_errType & (jec::kFlavorMask & ~jec::kFlavorQCD)) ); 
    errFlavor = _FlavorMixed(pTprime,eta,"QCD");
    err2 += errFlavor*errFlavor;
  }
  if (_errType & jec::kFlavorZJet) {
    assert( !(_errType & (jec::kFlavorMask & ~jec::kFlavorZJet)) ); 
    errFlavor = _FlavorMixed(pTprime,eta,"Z+jet");
    err2 += errFlavor*errFlavor;
  }
  if (_errType & jec::kFlavorPhotonJet) {
    assert( !(_errType & (jec::kFlavorMask & ~jec::kFlavorPhotonJet)) ); 
    errFlavor = _FlavorMixed(pTprime,eta,"photon+jet");
    err2 += errFlavor*errFlavor;
  }
  if (_errType & jec::kFlavorPureQuark) {
    assert( !(_errType & (jec::kFlavorMask & ~jec::kFlavorPureQuark)) ); 
    errFlavor = _FlavorMixed(pTprime,eta,"quark");
    err2 += errFlavor*errFlavor;
  }
  if (_errType & jec::kFlavorPureGluon) {
    assert( !(_errType & (jec::kFlavorMask & ~jec::kFlavorPureGluon)) ); 
    errFlavor = _FlavorMixed(pTprime,eta,"gluon");
    err2 += errFlavor*errFlavor;
  }
  if (_errType & jec::kFlavorPureCharm) {
    assert( !(_errType & (jec::kFlavorMask & ~jec::kFlavorPureCharm)) ); 
    errFlavor = _FlavorMixed(pTprime,eta,"charm");
    err2 += errFlavor*errFlavor;
  }
  if (_errType & jec::kFlavorPureBottom) {
    assert( !(_errType & (jec::kFlavorMask & ~jec::kFlavorPureBottom)) ); 
    errFlavor = _FlavorMixed(pTprime,eta,"bottom");
    err2 += errFlavor*errFlavor;
  }

  return errFlavor;
} // _Flavor


// Flavor systematic for predefined mixture
// relative to "20% glue" (barrel) and "QCD" (forward)
double JECUncertainty::_FlavorMixed(double pTprime, double eta,
					     string smix) const {

  // Calculate the "20% glue" reference point at 200 GeV in the barrel
  // (if applying pT-dependent L3Residual, could use pTprime instead)
  double refb(0);
  {
    double ptref = 200;
    double etaref = 0;

    double fL = _FlavorFraction(ptref, etaref, 0, 1);
    double fG = _FlavorFraction(ptref, etaref, 1, 1);
    double fC = _FlavorFraction(ptref, etaref, 2, 1);
    double fB = _FlavorFraction(ptref, etaref, 3, 1);
  
    refb = _FlavorMix(ptref, etaref, fL, fG, fC, fB);
  }

  // Calculate the QCD reference point at pTprime in the barrel
  // (if not applying pT-dependent L2Residual, need effective pTref instead)
  double qcdb(0);
  {
    double ptref = pTprime;
    double etaref = 0;
    
    double fL = _FlavorFraction(ptref, etaref, 0, 0);
    double fG = _FlavorFraction(ptref, etaref, 1, 0);
    double fC = _FlavorFraction(ptref, etaref, 2, 0);
    double fB = _FlavorFraction(ptref, etaref, 3, 0);
  
    qcdb = _FlavorMix(ptref, etaref, fL, fG, fC, fB);
  }

  // Calculate the QCD reference point at pTprime in the forward region
  // (if not applying pT-dependent L2Residual, need effective pTref instead)
  double qcdf(0);
  {
    double ptref = pTprime;
    
    double fL = _FlavorFraction(ptref, eta, 0, 0);
    double fG = _FlavorFraction(ptref, eta, 1, 0);
    double fC = _FlavorFraction(ptref, eta, 2, 0);
    double fB = _FlavorFraction(ptref, eta, 3, 0);
  
    qcdf = _FlavorMix(ptref, eta, fL, fG, fC, fB);
  }

  // Consider how calibration is propagated from Z+jet at 200 GeV
  // to QCD jet in barrel to QCD jet in the forward region
  double ref = refb + (qcdf-qcdb);

  // Calculate fractions and call _FlavorMix
  if (smix == "20% glue" ||  smix == "QCD" ||
      smix == "Z+jet" || smix == "photon+jet") {

    double pt = pTprime;
    if (smix=="20% glue")
      pt = 200;// Z+jet effective mean pT

    int iflavor(-1);
    if (smix=="QCD")        iflavor = 0;
    if (smix=="20% glue")   iflavor = 1; // Z/gamma+jet
    if (smix=="Z+jet")      iflavor = 3; // Zee+jet //2; // Zmm+jet
    if (smix=="photon+jet") iflavor = 4;

    double fL = _FlavorFraction(pt, eta, 0, iflavor);
    double fG = _FlavorFraction(pt, eta, 1, iflavor);
    double fC = _FlavorFraction(pt, eta, 2, iflavor);
    double fB = _FlavorFraction(pt, eta, 3, iflavor);
  
    return (_FlavorMix(pTprime, eta, fL, fG, fC, fB) - ref);
  }
  else if (smix=="quark")  {
    return _FlavorMix(pTprime,eta, 1, 0, 0, 0) - ref;
  }
  else if (smix=="gluon")  {
    return _FlavorMix(pTprime,eta, 0, 1, 0, 0) - ref;
  }
  else if (smix=="charm")  {
    return _FlavorMix(pTprime,eta, 0, 0, 1, 0) - ref;
  }
  else if (smix=="bottom") {
    return _FlavorMix(pTprime,eta, 0, 0, 0, 1) - ref;
  }
  else {
    assert(false);
  }

} // _Flavor


// Flavor systematic for any given mixture
double JECUncertainty::_FlavorMix(double pTprime, double eta,
					   double fl, double fg,
					   double fc, double fb) const {

  assert(fabs(fl+fg+fc+fb-1)<0.0001);
  assert(fl>=0 && fl<=1);
  assert(fg>=0 && fg<=1);
  assert(fc>=0 && fc<=1);
  assert(fb>=0 && fb<=1);

  double uL = _FlavorResponse(pTprime, eta, 0);
  double uG = _FlavorResponse(pTprime, eta, 1);
  double uC = _FlavorResponse(pTprime, eta, 2);
  double uB = _FlavorResponse(pTprime, eta, 3);

  double diff = fl*uL + fg*uG + fc*uC + fb*uB;

  return diff;
} // _FlavorMix


// Flavor response difference (Pythia/Herwig) for any given flavor
// Determined from QCD dijet MC
double JECUncertainty::_FlavorResponse(double pt, double eta,
				       int iflavor) const{

 // flavors: 0/uds, 1/gluon, 2/charm, 3/bottom (, 4/unmatched->1/gluon)
  // samples; 0/dijet, 1/Z/gamma+jet, 2/Zmm+jet, 3/Zee+jet 4/gamma+jet
  
  // drawFragFlavor::drawPureFlavor()
  // Sample = Dijet (eta,flavor(phys),par), alpha=0.2 (June 27)
  static double pFlavor[5][4][4] =
    {{{ 1.0030, -0.0005, -0.0245,  0.0000},
      { 1.0122, -0.0009,  0.5333,  0.0000},
      { 1.0034, -0.0005, -0.0245,  0.0000},
      { 1.0005, -0.0005, -0.0245,  0.0000}},
     {{ 1.0007,  0.0011,  0.1224,  0.0000},
      { 1.0124, -0.0043,  0.5622,  0.0000},
      { 1.0014,  0.0011,  0.1224,  0.0000},
      { 1.0006,  0.0011,  0.1224,  0.0000}},
     {{ 0.9939,  0.0094,  0.6444,  0.0000},
      { 1.0092,  0.0013,  0.6444,  0.0000},
      { 0.9944,  0.0094,  0.6444,  0.0000},
      { 0.9906,  0.0094,  0.6444,  0.0000}},
     {{ 0.9832,  0.0387,  1.4110,  0.0000},
      { 1.0028,  0.0269,  1.4110,  0.0000},
      { 0.9832,  0.0387,  1.4110,  0.0000},
      { 0.9832,  0.0387,  1.4110,  0.0000}},
     {{ 0.9970,  0.0129,  0.3621,  0.0000},
      { 1.0155,  0.0081,  0.3621,  0.0000},
      { 0.9970,  0.0129,  0.3621,  0.0000},
      { 0.9970,  0.0129,  0.3621,  0.0000}}};

  double *p;
  int ieta = 0;
  if (fabs(eta)<1.3) ieta = 0;
  if (fabs(eta)>=1.3 && fabs(eta)<2.5) ieta = 1;
  if (fabs(eta)>=2.5 && fabs(eta)<3.0) ieta = 2;
  if (fabs(eta)>=3.0 && fabs(eta)<3.2) ieta = 3;
  if (fabs(eta)>=3.2 && fabs(eta)<5.2) ieta = 4;
  p = &pFlavor[ieta][iflavor][0];

  // Ensure reliable pT range
  double ptmin = 30;
  double ptmax = 2000;
  double emax = 3000;
  ptmax = min(ptmax, emax/cosh(eta));
  double x = max(ptmin, min(ptmax, pt));

  double f = p[0]+p[1]*log10(0.01*x)+p[2]/x+p[3]/(x*x);

  return f;
} // _FlavorResponse

// Fraction of different flavors in various samples using Pythia6
// This uses the physics definition (algorithmic does not work for Herwig)
double JECUncertainty::_FlavorFraction(double pt, double eta,
				       int iflavor, int isample) const{

  // flavors: 0/uds, 1/gluon, 2/charm, 3/bottom
  // samples; 0/dijet, 1/z/gamma+jet, 2/Zmm+jet 3/Zee+jet 4/gamma+jet

  // drawFragFlavor::drawFractions()
  // Sample = Dijet (eta,flavor*2,par), alpha=0.2 (June 27)
  static double pDijet[5][8][4] =
    {{{ 0.2681,  0.1416,  0.5286, -0.1943},
      { 0.6716, -0.1239, -0.4865,  0.1786},
      { 0.0376, -0.0088, -0.0346,  0.0158},
      { 0.0222, -0.0028, -0.0228,  0.0095},
      { 0.3332,  0.1158,  0.2623, -0.1127},
      { 0.5703, -0.1271, -0.2496,  0.1280},
      { 0.0661, -0.0001,  0.0111, -0.0187},
      { 0.0300,  0.0124, -0.0251,  0.0042}},
     {{ 0.3500,  0.2959,  0.6164, -0.4113},
      { 0.5987, -0.2676, -0.5571,  0.3616},
      { 0.0331, -0.0303, -0.0053,  0.0095},
      { 0.0182, -0.0100, -0.0193,  0.0155},
      { 0.3902,  0.1816,  0.3468, -0.2645},
      { 0.5243, -0.1677, -0.3028,  0.2260},
      { 0.0587, -0.0148, -0.0093,  0.0153},
      { 0.0266, -0.0038, -0.0152,  0.0066}},
     {{ 0.5699,  0.4292, -0.0625,  0.0010},
      { 0.4017, -0.4050,  0.0563,  0.0020},
      { 0.0181, -0.0143, -0.0110,  0.0010},
      { 0.0091, -0.0050,  0.0052,  0.0005},
      { 0.4581,  0.3083, -0.0340,  0.0010},
      { 0.4729, -0.2925,  0.0459,  0.0020},
      { 0.0488, -0.0160,  0.0059,  0.0010},
      { 0.0201, -0.0074, -0.0046,  0.0005}},
     {{ 0.6010,  0.4530, -0.0100,  0.0010},
      { 0.3698, -0.4079, -0.0200,  0.0020},
      { 0.0194, -0.0342, -0.0100,  0.0010},
      { 0.0059,  0.0041, -0.0050,  0.0005},
      { 0.5452,  0.2126, -0.0100,  0.0010},
      { 0.4038, -0.1705, -0.0200,  0.0020},
      { 0.0359, -0.0048, -0.0100,  0.0010},
      { 0.0129,  0.0037, -0.0050,  0.0005}},
     {{ 0.7012,  0.3214, -0.0100,  0.0010},
      { 0.2766, -0.2810, -0.0200,  0.0020},
      { 0.0124, -0.0232, -0.0100,  0.0010},
      { 0.0046, -0.0078, -0.0050,  0.0005},
      { 0.5706,  0.1454, -0.0100,  0.0010},
      { 0.3856, -0.1318, -0.0200,  0.0020},
      { 0.0348, -0.0146, -0.0100,  0.0010},
      { 0.0078,  0.0159, -0.0050,  0.0005}}};

  // drawFragFlavor::drawFractions()
  // Sample = ZmmJet (eta,flavor*2,par), alpha=0.2 (June 27)
  static double pZmmJet[5][8][4] =
    {{{ 0.6841,  0.2204, -0.2097,  0.0935},
      { 0.2065, -0.1531,  0.3741, -0.2253},
      { 0.0625, -0.0482, -0.1222,  0.1786},
      { 0.0451, -0.0465, -0.0728,  0.0792},
      { 0.5458,  0.1337, -0.2113,  0.1023},
      { 0.3352, -0.1096,  0.3097, -0.1567},
      { 0.0717, -0.0132, -0.0842,  0.0875},
      { 0.0457, -0.0413, -0.0635,  0.0816}},
     {{ 0.6627,  0.2196, -0.1827,  0.0500},
      { 0.2470, -0.1321,  0.3700, -0.1297},
      { 0.0495, -0.0345, -0.0787, -0.1273},
      { 0.0397, -0.0542, -0.1389,  0.1126},
      { 0.5273,  0.1654, -0.0820, -0.1585},
      { 0.3693, -0.1174,  0.2863,  0.1857},
      { 0.0627, -0.0047, -0.1610, -0.3629},
      { 0.0399, -0.0452, -0.1224,  0.1183}},
     {{ 0.5722,  0.3495,  0.9643,  0.0010},
      { 0.3587, -0.2969, -0.5397,  0.0020},
      { 0.0384, -0.0869, -0.4412,  0.0010},
      { 0.0217,  0.0783,  0.2421,  0.0005},
      { 0.4660,  0.2667,  0.9236,  0.0010},
      { 0.4551, -0.2763, -0.6320,  0.0020},
      { 0.0464, -0.0685, -0.4831,  0.0010},
      { 0.0197,  0.0619,  0.2568,  0.0005}},
     {{ 0.6258,  0.2493, -0.0100,  0.0010},
      { 0.3371,  0.0600, -0.0200,  0.0020},
      { 0.0653,  0.0020, -0.0100,  0.0010},
      { 0.0657,  0.0930, -0.0050,  0.0005},
      { 0.5002,  0.1766, -0.0100,  0.0010},
      { 0.4643,  0.4190, -0.0200,  0.0020},
      { 0.0685, -0.0989, -0.0100,  0.0010},
      { 0.0592,  0.0493, -0.0050,  0.0005}},
     {{ 0.5921,  0.6777, -0.0100,  0.0010},
      { 0.3615, -0.6090, -0.0200,  0.0020},
      { 0.0293,  0.0192, -0.0100,  0.0010},
      { 0.0263,  0.0197, -0.0050,  0.0005},
      { 0.4623,  0.3278, -0.0100,  0.0010},
      { 0.4692, -0.3279, -0.0200,  0.0020},
      { 0.0476,  0.0528, -0.0100,  0.0010},
      { 0.0315,  0.0395, -0.0050,  0.0005}}};

  // drawFragFlavor::drawFractions()
  // Sample = ZeeJet (eta,flavor*2,par), alpha=0.3, phys only (June27)
  static double pZeeJet[5][8][4] =
    {{{ 0.6796,  0.1055, -0.4136,  0.2034},
      { 0.2117, -0.0107,  0.5236, -0.5488},
      { 0.0625, -0.0542, -0.1059,  0.0681},
      { 0.0477, -0.0354, -0.0936,  0.0553},
      { 0.6796,  0.1055, -0.4136,  0.2034},
      { 0.2117, -0.0107,  0.5236, -0.5488},
      { 0.0625, -0.0542, -0.1059,  0.0681},
      { 0.0477, -0.0354, -0.0936,  0.0553}},
     {{ 0.6719,  0.2185, -0.4875, -0.2509},
      { 0.2472, -0.0818,  0.4388, -0.4202},
      { 0.0453, -0.0804, -0.0552,  0.1770},
      { 0.0364, -0.0623, -0.0395,  0.1994},
      { 0.6712,  0.2174, -0.4643, -0.2005},
      { 0.2472, -0.0818,  0.4388, -0.4202},
      { 0.0453, -0.0804, -0.0552,  0.1770},
      { 0.0364, -0.0623, -0.0395,  0.1994}},
     {{ 0.6208,  0.2249, -0.3039,  0.0010},
      { 0.3204, -0.1559,  0.3343,  0.0020},
      { 0.0316, -0.0643, -0.0650,  0.0010},
      { 0.0198, -0.0327,  0.0113,  0.0005},
      { 0.6214,  0.2330, -0.2884,  0.0010},
      { 0.3202, -0.1599,  0.3259,  0.0020},
      { 0.0316, -0.0643, -0.0650,  0.0010},
      { 0.0198, -0.0327,  0.0113,  0.0005}},
     {{ 0.5998,  0.3859, -0.0100,  0.0010},
      { 0.3531, -0.3259, -0.0200,  0.0020},
      { 0.0042, -0.0954, -0.0100,  0.0010},
      { 0.0152, -0.0537, -0.0050,  0.0005},
      { 0.5998,  0.3859, -0.0100,  0.0010},
      { 0.3531, -0.3259, -0.0200,  0.0020},
      { 0.0042, -0.0954, -0.0100,  0.0010},
      { 0.0152, -0.0537, -0.0050,  0.0005}},
     {{ 0.5242,  0.3670, -0.0100,  0.0010},
      { 0.4170, -0.3651, -0.0200,  0.0020},
      { 0.0106, -0.0748, -0.0100,  0.0010},
      { 0.0121, -0.0375, -0.0050,  0.0005},
      { 0.5242,  0.3670, -0.0100,  0.0010},
      { 0.4170, -0.3651, -0.0200,  0.0020},
      { 0.0106, -0.0748, -0.0100,  0.0010},
      { 0.0121, -0.0375, -0.0050,  0.0005}}};

  // drawFragFlavor::drawFractions()
  // Sample = GJet (eta,flavor*2,par), alpha=0.2 (June 27)
  static double pGJet[5][8][4] =
    {{{ 0.7059,  0.0849,  0.0320, -0.1264},
      { 0.1157,  0.1101,  0.1160, -0.0568},
      { 0.1570, -0.1833, -0.0971,  0.1372},
      { 0.0210, -0.0163, -0.0285,  0.0252},
      { 0.5780,  0.0398,  0.0492, -0.1030},
      { 0.2330,  0.1178,  0.0399, -0.0543},
      { 0.1658, -0.1509, -0.0692,  0.1294},
      { 0.0231, -0.0076, -0.0173,  0.0258}},
     {{ 0.6814,  0.0423, -0.0287, -0.1011},
      { 0.1799,  0.1947,  0.0344, -0.0278},
      { 0.1197, -0.2006, -0.0071,  0.1058},
      { 0.0169, -0.0221, -0.0218,  0.0320},
      { 0.5661,  0.0293,  0.0326, -0.1311},
      { 0.2798,  0.1653, -0.0962,  0.0428},
      { 0.1331, -0.1589, -0.0014,  0.1120},
      { 0.0203, -0.0138, -0.0140,  0.0418}},
     {{ 0.5629, -0.1085, -0.2172,  0.0010},
      { 0.3523,  0.2662,  0.0244,  0.0020},
      { 0.0841, -0.2336,  0.0449,  0.0010},
      { 0.0084, -0.0137, -0.0090,  0.0005},
      { 0.4811, -0.1264,  0.1273,  0.0010},
      { 0.4135,  0.2715, -0.6560,  0.0020},
      { 0.1048, -0.1644,  0.0237,  0.0010},
      { 0.0134, -0.0048, -0.0321,  0.0005}},
     {{ 0.5580, -0.3131, -0.0100,  0.0010},
      { 0.2486, -0.2173, -0.0200,  0.0020},
      { 0.1101, -0.0154, -0.0100,  0.0010},
      { 0.0525, -0.0012, -0.0050,  0.0005},
      { 0.6201,  0.1566, -0.0100,  0.0010},
      { 0.5144,  1.2762, -0.0200,  0.0020},
      { 0.0824,  0.0543, -0.0100,  0.0010},
      { 0.0525, -0.0012, -0.0050,  0.0005}},
     {{ 0.3955, -0.0753, -0.0100,  0.0010},
      { 0.8985,  1.9568, -0.0200,  0.0020},
      { 0.1028,  0.0031, -0.0100,  0.0010},
      { 0.0386,  0.0338, -0.0050,  0.0005},
      { 0.5116,  0.3740, -0.0100,  0.0010},
      { 0.7278,  1.4165, -0.0200,  0.0020},
      { 0.0334, -0.1712, -0.0100,  0.0010},
      { 0.0047, -0.0514, -0.0050,  0.0005}}};

  // drawFragFlavor::drawFractions()
  // Sample = ZJet (eta,flavor*2,par), alpha=0.2 (Zee 0.3 + missing algo; June 27)
  static double pZJet[5][8][4] =
    {{{ 0.6821,  0.2125, -0.2349,  0.0574},
      { 0.2097, -0.1333,  0.3774, -0.1934},
      { 0.0643, -0.0253, -0.0623,  0.0328},
      { 0.0478, -0.0423, -0.0936,  0.0951},
      { 0.5422,  0.1602, -0.1219, -0.0055},
      { 0.3401, -0.1352,  0.1903, -0.0593},
      { 0.0718,  0.0140, -0.0238,  0.0012},
      { 0.0477, -0.0339, -0.0663,  0.0794}},
     {{ 0.6630,  0.2181, -0.3145,  0.0245},
      { 0.2491, -0.0936,  0.4686, -0.2604},
      { 0.0506, -0.0406, -0.0480,  0.0426},
      { 0.0395, -0.0448, -0.0678,  0.0804},
      { 0.5254,  0.1642, -0.0737, -0.1236},
      { 0.3747, -0.1048,  0.1828, -0.0703},
      { 0.0614,  0.0032,  0.0016, -0.0075},
      { 0.0413, -0.0302, -0.0885,  0.1181}},
     {{ 0.6066,  0.3226,  0.0440,  0.0010},
      { 0.3441, -0.1151,  0.2897,  0.0020},
      { 0.0354, -0.0547, -0.0653,  0.0010},
      { 0.0265, -0.0028,  0.0442,  0.0005},
      { 0.4673,  0.3349,  1.1103,  0.0010},
      { 0.4642, -0.1988, -0.4877,  0.0020},
      { 0.0515,  0.0175, -0.2655,  0.0010},
      { 0.0291,  0.0191,  0.0561,  0.0005}},
     {{ 0.6392,  0.4737, -0.0100,  0.0010},
      { 0.3840, -0.2258, -0.0200,  0.0020},
      { 0.0354, -0.0268, -0.0100,  0.0010},
      { 0.0241, -0.0419, -0.0050,  0.0005},
      { 0.5289,  0.3133, -0.0100,  0.0010},
      { 0.5056,  0.6265, -0.0200,  0.0020},
      { 0.0577, -0.1386, -0.0100,  0.0010},
      { 0.0582,  0.0430, -0.0050,  0.0005}},
     {{ 0.5674,  0.4754, -0.0100,  0.0010},
      { 0.4476, -0.2899, -0.0200,  0.0020},
      { 0.0207, -0.0494, -0.0100,  0.0010},
      { 0.0191, -0.0226, -0.0050,  0.0005},
      { 0.4823,  0.3806, -0.0100,  0.0010},
      { 0.5002, -0.2172, -0.0200,  0.0020},
      { 0.0522,  0.0742, -0.0100,  0.0010},
      { 0.0397,  0.0747, -0.0050,  0.0005}}};

  double *p(0);
  int ieta = 0;
  if (fabs(eta)<1.3) ieta = 0;
  if (fabs(eta)>=1.3 && fabs(eta)<2.5) ieta = 1;
  if (fabs(eta)>=2.5 && fabs(eta)<3.0) ieta = 2;
  if (fabs(eta)>=3.0 && fabs(eta)<3.2) ieta = 3;
  if (fabs(eta)>=3.2 && fabs(eta)<5.2) ieta = 4;
  if (isample==0) p = &pDijet[ieta][iflavor][0];
  if (isample==1) p = &pZJet[ieta][iflavor][0];
  if (isample==2) p = &pZmmJet[ieta][iflavor][0];
  if (isample==3) p = &pZeeJet[ieta][iflavor][0];
  if (isample==4) p = &pGJet[ieta][iflavor][0];

  // Ensure reliable pT range
  double ptmin = 30;
  if (isample==0) ptmin = 50;
  double ptmax = 2000;
  if (isample==1) ptmax = 800;
  if (isample==2) ptmax = 600;
  if (isample==3) ptmax = 600;
  if (isample==4) ptmax = 800;
  double emax = 3000;
  if (isample==1) emax = 2000;  
  if (isample==2) emax = 1500;
  if (isample==3) emax = 1500;
  if (isample==4) emax = 2000;  
  ptmax = min(ptmax, emax/cosh(eta));
  double x = max(ptmin, min(ptmax, pt));

  double f = p[0] + p[1]*log10(0.01*x) + p[2]*pow(log10(0.01*x),2)
    + p[3]*pow(log10(0.01*x),3);
  f = max(0., min(1., f));

  // Calculate sum of all flavors to get correct normalization
  double sumf(0);
  const int nf = 4;
  for (int jf = 0; jf != nf; ++jf) {
    
    double *pi(0);
    const int kf = iflavor/nf;
    if (isample==0) pi = &pDijet[ieta][nf*kf+jf][0];
    if (isample==1) pi = &pZJet[ieta][nf*kf+jf][0];
    if (isample==2) pi = &pZmmJet[ieta][nf*kf+jf][0];
    if (isample==3) pi = &pZeeJet[ieta][nf*kf+jf][0];
    if (isample==4) pi = &pGJet[ieta][nf*kf+jf][0];
    double fi = pi[0] + pi[1]*log10(0.01*x) + pi[2]*pow(log10(0.01*x),2)
      + pi[3]*pow(log10(0.01*x),3);
    fi = max(0., min(1., fi));
    sumf += fi;
  } // for jf
  
  // Normalize sum of flavor fractions to 1, if too large
  if (sumf>1)
    f /= sumf;
  // Add missing fraction (undefined, usually) to gluons
  if (sumf<1 && (iflavor==1 || iflavor==5))
    f += (1-sumf);

  return f;
} // _FlavorFraction

// Time-dependence uncertainties from L2 and L3
double JECUncertainty::_Time(const double pt, const double eta) const {

  // Optional epoch time uncertainties
  double spt(0);
  if (_errType & jec::kTimePtRunA) {
    assert( !(_errType & (jec::kTimePtMask & ~jec::kTimePtRunA)) ); 
    spt = _TimePt(pt,1);
  }
  if (_errType & jec::kTimePtRunB) {
    assert( !(_errType & (jec::kTimePtMask & ~jec::kTimePtRunB)) ); 
    spt = _TimePt(pt,2);
  }
  if (_errType & jec::kTimePtRunC) {
    assert( !(_errType & (jec::kTimePtMask & ~jec::kTimePtRunC)) ); 
    spt = _TimePt(pt,3);
  }
  if (_errType & jec::kTimePtRunD) {
    assert( !(_errType & (jec::kTimePtMask & ~jec::kTimePtRunD)) ); 
    spt = _TimePt(pt,4);
  }
  // Normal time uncertainties
  if (_errType & jec::kTimePt) {
    assert( !(_errType & (jec::kTimePtMask & ~jec::kTimePt)) ); 
    spt  =   _TimePt(pt);
  }
  double seta =   (_errType & jec::kTimeEta ? _TimeEta(eta) : 0);

  double err = sqrt(seta*seta + spt*spt);

  // signed sources
  if (!(_errType & ~jec::kTimeEta)) return seta;
  if (!(_errType & ~jec::kTimePt))  return spt;
  if (!(_errType & ~jec::kTimePtRunA))  return spt;
  if (!(_errType & ~jec::kTimePtRunB))  return spt;
  if (!(_errType & ~jec::kTimePtRunC))  return spt;
  if (!(_errType & ~jec::kTimePtRunD))  return spt;

  return err;
}

// Time-dependence uncertainty from L2
double JECUncertainty::_TimeEta(const double eta) const {

  // Time dependence systematics from Henning Kirschenmann
  // Subject: 	Re: updated systematics for 53X
  // Date: 	April 23, 2013 1:17:12 PM GMT+03:00
  // Arrays (rounded to 0.01%) produced from histogram
  // Output53Z2FullReRecoJECBin2012_11TimeDependence.root:/AbsMPFVsRunNumber10/ResolutionPlots_MPFResponseVsRunNumber_DeviationsOfRatioVsBinVar40_AbsMPFVsRunNumber40_DataDev

  const int neta = 7;
  double etabins[neta+1] =
    {0, 0.783, 1.305, 1.93, 2.5, 2.964, 3.2, 5.5};//5.191};
  double errs_mpf[neta] =
    {0.03, 0.06, 0.36, 0.51, 0.70, 1.10, 0.88};

  double err(0.1);
  for (int i = 0; i != neta; ++i) {
    if (fabs(eta)>=etabins[i] && fabs(eta)<etabins[i+1])
      err = 0.01*errs_mpf[i];
  }

  // Add 3% in quadrature for Calo and JPT to cover for
  // Calo/PF ratio vs eta after L1L2L3Res
  if (_calo||_jpt) err = sqrt(err*err + 0.03*0.03);

  return err;
} // Time

// Time-dependence uncertainty from L3
// based on isolated charged hadron E/p in barrel from Matthieu
// On 23 Jul 2014, at 14:59, Matthieu Marionneau
// Re: HCAL raddam studies: where are we?
double JECUncertainty::_TimePt(const double pt, int epoch) const {
  // epochs: 0 (RMS all), 1 (runA/all), 2 (runB/all), 3 (runC/all), 4 (runD/all)
  assert(epoch>=0 && epoch<=4);

  const int nepoch = 4;
  // Values eye-balled from Eop_BH_3_10.pdf
  double pt3[nepoch]  = {0.667, 0.661, 0.655, 0.651};
  double pt10[nepoch] = {0.815, 0.797, 0.782, 0.778};
  double lum[nepoch] = {851.3, 4411.7, 7031.7, 7185.7};

  double sumlum(0);
  for (int i = 0; i != nepoch; ++i) sumlum += lum[i];
  if (debug) cout << "Total lumi: " <<  sumlum/1000. << " fb-1" <<  endl;

  if (debug) cout << "Epoch relative weights:" << endl;
  double w[nepoch];
  double scale3(0), scale10(0);
  for (int i = 0; i != nepoch; ++i) {
    w[i] = lum[i] / sumlum;
    scale3  += w[i] * pt3[i];
    scale10 += w[i] * pt10[i];
    if (debug) cout << Form("%5.1f%% ",w[i]*100);
  }
  if (debug) cout << endl;

  if (debug) cout << "Epoch relative jes: " << endl;
  double r3[nepoch];
  double r10[nepoch];
  double rms2(0);
  for (int i = 0; i != nepoch; ++i) {
    r3[i]  = pt3[i] / scale3;
    r10[i] = pt10[i] / scale10;
    rms2 += w[i] * pow(r10[i]-1,2);
  }
  if (debug) {
    for (int i = 0; i != nepoch; ++i) cout << Form("%6.3f ",r3[i]); cout<<endl;
    for (int i = 0; i != nepoch; ++i) cout << Form("%6.3f ",r10[i]); cout<<endl;
  }

  double hb = (100-4.42)/100.;//V3PT
  if (debug) cout << Form("Epoch absolute jes (Full 2012 is %1.3f):\n",hb);
  double jes3[nepoch];
  double jes10[nepoch];
  for (int i = 0; i != nepoch; ++i) {
    jes3[i]  = r3[i] * hb;
    jes10[i] = r10[i] * hb;
  }
  if (debug) {
    for (int i = 0; i != nepoch; ++i) cout<<Form("%6.3f ",jes3[i]); cout<<endl;
    for (int i = 0; i != nepoch; ++i) cout<<Form("%6.3f ",jes10[i]); cout<<endl;
    cout << "http://cms-physics.web.cern.ch/cms-physics/public/JME-10-008-pas.pdf Figure 4:" << endl;
    cout << "Barrel(pT=3 GeV)~0.98; Barrel(pT=10 GeV)~1.00" << endl;
  }

  double rms = sqrt(rms2);
  if (debug) cout << Form("RMS: %1.2f%%",rms*100.) << endl;
  
  double hbnew = (hb-1) + (epoch==0 ? rms : r10[epoch-1]-1);
  if (debug) cout << Form("hbnew=%1.3f",hbnew) << endl;
  double err = _jeshb(pt, hbnew);

  return err;
}

// Random Cone offset for data
double JECUncertainty::_L1DataFlat(const double pT, const double eta) {

  //L1Offset VO Data
  assert(_jecL1DTflat);
  _jecL1DTflat->setRho(_RhoFromMu(_mu));
  _jecL1DTflat->setJetA(TMath::Pi()*0.5*0.5*_ajet);
  _jecL1DTflat->setJetEta(eta);
  _jecL1DTflat->setJetPt(pT);
  return ( _jecL1DTflat->getCorrection() );
}

// Random Cone offset for MC
double JECUncertainty::_L1MCFlat(const double pT, const double eta) {

  //L1Offset VO MC
  assert(_jecL1MCflat);
  _jecL1MCflat->setRho(_RhoFromMu(_mu));
  _jecL1MCflat->setJetA(TMath::Pi()*0.5*0.5*_ajet);
  _jecL1MCflat->setJetEta(eta);
  _jecL1MCflat->setJetPt(pT);
  return ( _jecL1MCflat->getCorrection() );
}

// Scaled offset for data (MC truth * RandomConeData / RandomConeMC)
double JECUncertainty::_L1Data(const double pT, const double eta) {

  //L1Offset V5 Data
  assert(_jecL1DTpt);
  _jecL1DTpt->setRho(_RhoFromMu(_mu));
  _jecL1DTpt->setJetA(TMath::Pi()*0.5*0.5*_ajet);
  _jecL1DTpt->setJetEta(eta);
  _jecL1DTpt->setJetPt(pT);
  return ( _jecL1DTpt->getCorrection() );
}

// Scaled offset for MC (MC truth)
double JECUncertainty::_L1MC(const double pT, const double eta) {

  assert(_jecL1MCpt);
  _jecL1MCpt->setRho(_RhoFromMu(_mu));
  _jecL1MCpt->setJetA(TMath::Pi()*0.5*0.5*_ajet);
  _jecL1MCpt->setJetEta(eta);
  _jecL1MCpt->setJetPt(pT);
  return ( _jecL1MCpt->getCorrection() );
}

// Data/MC scale factor as function of rho
double JECUncertainty::_L1SF(const double pT, const double eta,
			     const double rho) {

  assert(_jecL1sf);
  _jecL1sf->setRho(rho);
  _jecL1sf->setJetA(TMath::Pi()*0.5*0.5*_ajet);
  _jecL1sf->setJetEta(eta);
  _jecL1sf->setJetPt(pT);
  return ( _jecL1sf->getCorrection() );
}

// Move to mu-based mapping, which is better for comparing
// different PU scenarios as it considers both IT and OOT PU,
// plus we have a number directly comparable to ATLAS
double JECUncertainty::_RhoFromMu(double mu) {
  // Eta_0.0-1.3, jt320
  return (1.01272 + 0.551183*mu + 0.000362936*mu*mu);
}
double JECUncertainty::_NpvFromMu(double mu) {
  // Eta_0.0-1.3, jt400
  return (0.851334 + 0.722608*mu - 0.00184534*mu*mu);
}
