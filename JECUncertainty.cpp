#include "JECUncertainty.hpp"
#include "TMath.h"
#include "TFile.h"
#include "TF1.h"
#include "Math/BrentRootFinder.h"
#include "Math/RootFinderAlgorithms.h"

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
				      const double npv) :
  _algo(algo), _type(type), _errType(errType), _npv(npv)
{

  this->SetJetAlgo(algo); // initialize tables for jetalg
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
  if (_errType & jec::kPileUp) {
    errPileUp = _PileUp(pTprime, eta2);
    err2 += errPileUp * errPileUp;
  }
  //if (_errType & jec::kFlavor) {
  if (_errType & jec::kFlavorMask) {
    errFlavor = _Flavor(pTprime, eta);
    err2 += errFlavor*errFlavor;
  }
  if (_errType & jec::kTime) {
    errTime = _Time(eta2);
    err2 += errTime * errTime;
  }
  
  double err = sqrt(err2);

  // if requesting single source, return signed for sign-changing cases
  if (!(_errType & ~jec::kAbsoluteFrag)) return errAbs;
  if (!(_errType & ~jec::kAbsoluteSPRE)) return errAbs;
  if (!(_errType & ~jec::kAbsoluteSPRH)) return errAbs;
  if (!(_errType & ~jec::kRelativePtBB)) return errRel;
  if (!(_errType & ~jec::kRelativePtEC1)) return errRel;
  if (!(_errType & ~jec::kRelativePtEC2)) return errRel;
  if (!(_errType & ~jec::kRelativePtHF))  return errRel;
  if (!(_errType & ~jec::kRelativePt))    return errRel; // EXTRA
  if (!(_errType & ~jec::kPileUpDataMC)) return errPileUp;
  if (!(_errType & ~jec::kPileUpBias))   return errPileUp;
  if (!(_errType & ~jec::kPileUpPtBB))   return errPileUp;
  if (!(_errType & ~jec::kPileUpPtEC))   return errPileUp;
  if (!(_errType & ~jec::kPileUpPtHF))   return errPileUp;
  if (!(_errType & ~jec::kPileUpPt))     return errPileUp; // EXTRA
  //
  if (!(_errType & ~jec::kFlavorMask))   return errFlavor;

  return err;
}


void JECUncertainty::SetJetAlgo(const jec::JetAlgo& algo) {

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

}


void JECUncertainty::_InitL1() {

  // MC truth L1 (MC_V1) and L2L3 files from Ricardo Euseby by e-mailed weblink
  // Subject: 	Re: New MC-based L1
  // Date: 	May 9, 2013 10:50:54 PM GMT+03:00
  // http://people.physics.tamu.edu/eusebi/jec/53X/SQLfiles/Summer13/

  // RandomCone (V0) and scaled L1 (DATA_V1) files from Ia Iashvili by e-mail
  // Subject: 	Re: New MC-based L1
  // Date: 	May 14, 2013 8:54:27 PM GMT+03:00

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

  // L1Offset systematics
  {
    const char *s = Form("%sSummer13_V0_DATA_L1Offset_%s_pt.txt",d,a);
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l1);
    _jecL1nominal = new FactorizedJetCorrector(v);
  }
  {
    const char *s = Form("%sSummer13_V0_MC_L1Offset_%s_pt.txt",d,a);
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l1);
    _jecL1MCnominal = new FactorizedJetCorrector(v);
  }
  {
    const char *s = Form("%sSummer13_V1_DATA_L1FastJet_%s.txt",d,a);
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l1);
    _jecL1scaled = new FactorizedJetCorrector(v);
  }
  {
    const char *s = Form("%sSummer13_V1_MC_L1FastJet_%s.txt",d,a);
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l1);
    _jecL1pt = new FactorizedJetCorrector(v);
  }
} // InitL1


void JECUncertainty::_InitJEC() {

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
  s = Form("%sSummer13_V1_DATA_L1FastJet_%s.txt",d,a);
  if (debug) cout << s << endl << flush;
  JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
  s = Form("%sSummer13_V1_DATA_L2Relative_%s.txt",d,a);
  if (debug) cout << s << endl << flush;
  JetCorrectorParameters *l2 = new JetCorrectorParameters(s);
  s = Form("%sSummer13_V1_DATA_L3Absolute_%s.txt",d,a);
  if (debug) cout << s << endl << flush;
  JetCorrectorParameters *l3 = new JetCorrectorParameters(s);
  s = Form("%sSummer13_V1_DATA_L2L3Residual_AK5PF.txt.V4handMadeUpdatedL3",d);
  if (_algo==jec::AK5CALO || _algo==jec::AK7CALO) {
    s = Form("%sWinter12_V1_DATA_L2L3Residual_%s.txt.kFSRone",d,a);
  }  
  if (debug) cout << s << endl << flush;
  JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);

  vector<JetCorrectorParameters> v;
  v.push_back(*l1);
  v.push_back(*l2);
  v.push_back(*l3);
  v.push_back(*l2l3res);
  _jecDefault = new FactorizedJetCorrector(v);
  _jec = _jecDefault;

  // Another versions using Random Cone offset (L1 V0)
  s = Form("%sSummer13_V0_DATA_L1Offset_%s_pt.txt",d,a);
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

  // JER, PT, statistical uncertainty from Denis Rathjens
  // Subject: 	Re: New SQL Files
  // Date: 	May 21, 2013 9:02:44 PM GMT+03:00

  // (Winter12(13) files for PF,PFCHS,CALO,JPT from Denis Rathjens by e-mail)
  // (Subject: 	Re: l2res in rereco)
  // (Date: 	April 22, 2013 5:04:16 PM GMT+03:00)

  // Fixed PTDEP file from Denis Rathjens by e-mail
  // Subject: 	Re: l2res in rereco
  // Date: 	April 23, 2013 5:57:54 PM GMT+03:00

  // JER variations from Denis Ratjhens
  // Subject: 	Re: updated systematics for 53X
  // Date: 	April 23, 2013 3:53:01 PM GMT+03:00

  map<jec::JetAlgo, const char*> names;
  names[jec::AK5PF] = "AK5PF";
  names[jec::AK5PFchs] = "AK5PFchs";
  names[jec::AK5CALO] = "AK5Calo";
  // No L2Res systematics files for AK7, use AK5
  names[jec::AK7PF] = "AK5PF";
  names[jec::AK7PFchs] = "AK5PFchs";
  names[jec::AK7CALO] = "AK5Calo";

  string directory = "CondFormats/JetMETObjects/data/";
  const char *d = directory.c_str();

  const char *s;
  {
    s = Form("%sSummer13_V1_DATA_L2L3Residual_AK5PF.txt.PTDEPENDENCE",d);
    if (_algo==jec::AK5CALO || _algo==jec::AK7CALO)
      s = Form("%sSummer13_V1_DATA_L2L3Residual_AK5CALO.txt.PTDEPENDENCE",d);
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l2l3res);
    _jecL2ResFlat = new FactorizedJetCorrector(v);
  }
  {
    s = Form("%sSummer13_V1_DATA_L2L3Residual_AK5PF.txt.FineBinnedResidual",d);
    if (_algo==jec::AK5CALO || _algo==jec::AK7CALO)
      s = Form("%sSummer13_V1_DATA_L2L3Residual_AK5CALO.txt.FineBinnedResidual",d);
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l2l3res);
    _jecL2ResPt = new FactorizedJetCorrector(v);
  }
  {
    s = Form("%sSummer13_V1_DATA_L2L3Residual_AK5PF.txt.StatisticalUncertainties",d);
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l2l3res);
    _jecL2stat = new FactorizedJetCorrector(v);
  }

  {
    s = Form("%sSummer13_V1_DATA_L2L3Residual_AK5PF.txt.JERvariationUp",d);
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l2l3res);
    _jecL2jerup = new FactorizedJetCorrector(v);
  }
  {
    s = Form("%sSummer13_V1_DATA_L2L3Residual_AK5PF.txt.JERvariationDown",d);
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l2l3res);
    _jecL2jerdw = new FactorizedJetCorrector(v);
  }

} // InitL2Res


void JECUncertainty::SetNPV(const double npv) {
  _npv = npv;
}


void JECUncertainty::SetErrType(const jec::ErrorTypes& errType) {
  _errType = errType;
}


jec::ErrorTypes JECUncertainty::GetErrType() {
  return _errType;
}


// Solve pTraw from pTprime = pTraw / R(pTraw) using Brent's method
// We want to provide JEC uncertainties vs pTprime, not pTraw, but JEC
// is only available as a function of pTraw
double JECUncertainty::_Rjet(const double pTprime, const double eta) {
  
  _pTprime = pTprime;
  _eta = eta;
  
  ResponseFunc f(pTprime,_jec,_npv,_eta,_Rho(_npv),TMath::Pi()*0.5*0.5*_ajet);
  
  ROOT::Math::Roots::Brent brf;
  //brf.SetLogScan(true);
  // Set parameters of the method
  brf.SetFunction(f,std::max(2.0,0.25*pTprime),std::min(4*pTprime,3500.0));
  bool found_root = brf.Solve(50,1e-4,1e-5);
  double pTraw = brf.Root();
  _rjet = pTraw /pTprime;
  _jec->setJetE(pTraw*cosh(eta)); 
  _jec->setJetEta(eta); _jec->setNPV(_npv);
  _jec->setRho(_Rho(_npv)); _jec->setJetA(TMath::Pi()*0.5*0.5*_ajet);
  _jec->setJetPt(pTraw);
  double corr = _jec->getCorrection();
  if(std::abs((pTraw * corr - pTprime)/pTprime) > 0.001) {
    std::cout << "NPV:" << _npv << "   Brent: status:" << brf.Status() << " " << brf.Iterations() << '\n';
    std::cout << "R: pTprime:" << pTprime << "  eta:" << eta << " pTraw:" << pTraw << " corr:" << corr 
              << " pTcor:" << corr * pTraw << " rjet*pTprime:" << _rjet * pTprime << '\n'; 
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
  return _rjet;
} // _Rjet


// Combine pT dependent absolute scale uncertainties
double JECUncertainty::_Absolute(const double pt) const {

  double abs = (_errType & jec::kAbsoluteScale ? _AbsoluteScale() : 0.);
  double spr = (_errType & jec::kAbsoluteSPR ? _AbsoluteSPR(pt) : 0.);
  double frag = (_errType & jec::kAbsoluteFrag ? _AbsoluteFrag(pt) : 0.);

  // signed sources
  if (!(_errType & ~jec::kAbsoluteSPRE)) return spr;
  if (!(_errType & ~jec::kAbsoluteSPRH)) return spr;

  return sqrt(abs*abs + spr*spr + frag*frag);
} // Absolute


// Absolute scale from photon+jet and Z+jet
double JECUncertainty::_AbsoluteScale() const {

  // Combined uncertainty for Zmumu reference scale,
  // MPF method bias, and flavor mapping to QCD mixture

  // https://indico.cern.ch/getFile.py/access?contribId=4&resId=0&materialId=slides&confId=216765 
  // also here: http://www-ekp.physik.uni-karlsruhe.de/~dhaitz/plots_archive/2012_11_13/extrapolation/ratio_extrapolation_alpha_0_30__AK5PFCHSL1L2L3Res.png
  double zscale = 0.001; // data/MC difference, s11
  double mpfbias = 0.002; // alpha=0.1->0.0 difference, s10
  double stat = 0.002; // for alpha<0.30 fit, s10
  
  // This needs to be updated
  //double flavormap = 0.005; // Pythia/Herwig difference
  double flavormap = 0.; // for "20% glue" reference point at pT=200 GeV
  // Estimate also these from Pythia/Herwig difference
  // of genMET with (a) all particles, (b) all particles at |eta|<5,
  // (c) no neutrinos and no |eta|>5 particles
  // (not done properly, these are just rough estimates)
  double beambias = 0.002; // (a)-(b)
  double neutrinos = 0.002; // (b)-(c)
  
  double err2(0);
  err2 += zscale*zscale + mpfbias*mpfbias + stat*stat;
  err2 += flavormap*flavormap;
  err2 += beambias*beambias + neutrinos*neutrinos;

  return sqrt(err2);
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
    double pTref = 200.; // effective <pT> for Zgamma+jet data
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
  double refSPR(0), refSPRE(0), refSPRH(0);
  double refpt = 200.;

  TF1 *f = new TF1("fmore","max(0,[0]+[1]*pow(x,[2]))",10,3500);
  if (_errType & jec::kAbsoluteSPR) { // this is now obsolete
    f->SetParameters(1.02829e+00, -6.22540e-02, -2.67123e-01);
    if (_calo) f->SetParameters(1.02084e+00, 3.23027e-02, -9.46202e-01);
    errSPR = f->Eval(pTprime);
    refSPR = f->Eval(refpt);
  }
  if (_errType & jec::kAbsoluteSPRE) { // SPR in ECAL
    f->SetParameters(1.00567e+00, -3.04275e-02, -6.75493e-01);
    if (_calo) f->SetParameters(1.00166e+00, 1.57065e-02, -2.06585e-01);
    errSPRE = f->Eval(pTprime);
    refSPRE = f->Eval(refpt);
  }
  if (_errType & jec::kAbsoluteSPRH) { // SPR in HCAL
    f->SetParameters(1.03091e+00, -5.11540e-02, -1.54227e-01);
    if (_calo) f->SetParameters(1.02246e+00, -1.55689e-02, -1.17219e-01);
    errSPRH = f->Eval(pTprime);
    refSPRH = f->Eval(refpt);
  }

  // replace directly done SPR with pieces broken up into ECAL and HCAL
  errSPR = sqrt(pow(errSPRE-refSPRE,2) + pow(errSPRH-refSPRH,2));

  // signed sources
  if (!(_errType & ~jec::kAbsoluteSPRE)) return (errSPRE-refSPRE);
  if (!(_errType & ~jec::kAbsoluteSPRH)) return (errSPRH-refSPRH);

  return errSPR;
} // AbsoluteSPR


// Combine relative uncertainties
double JECUncertainty::_Relative(const double pTprime,
					const double eta) const {
  
  double sjer = (_errType & jec::kRelativeJER ? _RelativeJER(eta) : 0.);
  double sfsr = (_errType & jec::kRelativeFSR ? _RelativeFSR(eta) : 0.);
  double stat = (_errType & jec::kRelativeStat ? _RelativeStat(eta) : 0.);
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
double JECUncertainty::_RelativeJER(const double eta) const {

  // JER variations from Denis Rathjens, Nov 28, 2012, 10:22:27pm
  // Re: time and pT dependence systematics

  double x = fabs(eta);
  _jecL2jerup->setJetEta(x);
  _jecL2jerup->setJetPt(100.); // pT doesn't matter
  double up = _jecL2jerup->getCorrection();
  _jecL2jerdw->setJetEta(x);
  _jecL2jerdw->setJetPt(100.); // pT doesn't matter
  double dw = _jecL2jerdw->getCorrection();
  double err = 0.5*(up-dw);

  if (x<1.5) return 0;
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
double JECUncertainty::_RelativeStat(const double eta) const {

  double x = fabs(eta);

  // Statistical uncertainties from Denis Rathjens by e-mail
  // Subject: 	Re: l2res in rereco
  // Date: 	April 22, 2013 5:04:16 PM GMT+03:00

  _jecL2stat->setJetEta(x);
  _jecL2stat->setJetPt(100.); // pT doesn't matter
  double err = _jecL2stat->getCorrection();
  if ( (x>=2.5 && x<3. && (_errType & jec::kRelativeStatEC2)) ||
       (x>=3. && (_errType & jec::kRelativeStatHF)) )
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

  double abseta = fabs(eta);
  _jecL2ResFlat->setJetPt(pt);
  _jecL2ResFlat->setJetEta(abseta);
  double corrflat_p = _jecL2ResFlat->getCorrection();
  double ptraw_p =  pt / corrflat_p; // close enough in first approx
  _jecL2ResPt->setJetPt(ptraw_p);
  _jecL2ResPt->setJetEta(abseta);
  double corrpt_p = _jecL2ResPt->getCorrection();

  // Symmetrize systematics
  _jecL2ResFlat->setJetPt(pt);
  _jecL2ResFlat->setJetEta(-abseta);
  double corrflat_m = _jecL2ResFlat->getCorrection();
  double ptraw_m =  pt / corrflat_m; // close enough in first approx
  _jecL2ResPt->setJetPt(ptraw_m);
  _jecL2ResPt->setJetEta(-abseta);
  double corrpt_m = _jecL2ResPt->getCorrection();

  double corrpt = 0.5*(corrpt_m + corrpt_p);
  double corrflat = 0.5*(corrflat_m + corrflat_p);


  double err = (corrpt/corrflat - 1);
  double kfactor_hf = 0.50; // HF slope about same as slope stat
  if (abseta>=2.964) err *= kfactor_hf;
  double kfactor_ec = 0.50; // EC slope corrected
  if (abseta>=1.305 && abseta<2.964) err *= kfactor_ec;
  double x = fabs(eta);

  // Note 2014-Jan-10: was kRelativePtBB missing from published total?
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
  double sbias = (_errType & jec::kPileUpBias ? _PileUpBias(pTprime, eta) : 0.);
  double spt =   (_errType & jec::kPileUpPt ? _PileUpPt(pTprime, eta) : 0.);
  double err = sqrt(smc*smc + sbias*sbias + spt*spt);

  // signed sources
  if (!(_errType & ~jec::kPileUpDataMC)) return smc;
  if (!(_errType & ~jec::kPileUpBias)) return sbias;

  if (!(_errType & ~jec::kPileUpPtBB)) return (eta>=0.0 && eta<1.3) ? spt : 0;
  if (!(_errType & ~jec::kPileUpPtEC)) return (eta>=1.3 && eta<3.0) ? spt : 0;
  if (!(_errType & ~jec::kPileUpPtHF)) return (eta>=3.0 && eta<5.5) ? spt : 0;

  if (!(_errType & ~jec::kPileUpPt)) return spt; // EXTRA

  return err;
} // PileUp


// Pile-up uncertainty from data / MC difference
double JECUncertainty::_PileUpDataMC(const double pTprime, const double eta) {

  // Use 20% as uncertainty since we now have separate data and MC corrections
  // This should cover for uncertainty in the L1 scale factor
  double kfactor = 0.20;
  double pT_p = _Rjet(pTprime, eta) * pTprime;
  double l1d_p = _L1Data(pT_p, eta); // V1MC * V0Data / V0MC
  double l1m_p = _L1MC(pT_p, eta); // V1MC
  double sys_p = kfactor * (l1m_p / l1d_p - 1);
  //
  double l1db_p = _L1DataRaw(pT_p, eta);
  double l1mb_p = _L1MCRaw(pT_p, eta);
  double sysb_p = kfactor * (l1mb_p / l1db_p - 1);

  double pT_m = _Rjet(pTprime, -eta) * pTprime;
  double l1d_m = _L1Data(pT_m, -eta); // V1MC * V0Data / V0MC
  double l1m_m = _L1MC(pT_m, -eta); // V1MC
  double sys_m = kfactor * (l1m_m / l1d_m - 1);
  //
  double l1db_m = _L1DataRaw(pT_m, -eta);
  double l1mb_m = _L1MCRaw(pT_m, -eta);
  double sysb_m = kfactor * (l1mb_m / l1db_m - 1);


  // Symmetrize uncertainty
  double sys = 0.5*(sys_p + sys_m);
  double sysb = 0.5*(sysb_p + sysb_m);

  // Check that Data/MC difference consistent for scaled and nominal files
  if (fabs(sys-sysb)>0.001) {
    if (debug) {
      cout << "In PileUpDataMC(pt="<<pT_p<<",eta="<<eta<<"): ";
      cout << "sys="<<sys<<" sysb="<<sysb<<" diff="<<fabs(sys-sysb)
	   << " ratio="<<(sys-sysb)/sys<<endl<<flush;
    }
    // Fail if difference outside comfort zone;
    //assert(fabs(sys-sysb)<0.002 || fabs(sys-sysb)<0.20*fabs(sys));
  }

  return sys * _ajet;
} // PileUpDataMC


// Pile-up uncertainty from poorly understood MC bias (obsolete)
// Update: now though to come from jet pT dependence => in PileUpPt
double JECUncertainty::_PileUpBias(const double pTprime, const double eta) {
  return 0;
}


// Pile-up uncertainty from pT dependence
// Implemented as difference between MC truth (V1) and Random Cone (V0)
double JECUncertainty::_PileUpPt(const double pTprime, const double eta) {

  // Limit eta to [-5,5] because V0 files don't go further out
  // Even closer in, because bias goes nuts in the last bin out
  double maxeta = 4.5;
  double etax = max(-maxeta,min(maxeta,eta));

  // Given the detailed correction, kfactor of 50% should be justified already
  // ...Nope, couldn't be confirmed to better than 100% with Z+jet
  // Need to study PU profile bias vs NPV more (do vs instlum per LS instead?)
  double kfactor = 1.0;
  _jec = _jecDefault;
  double pT_p = _Rjet(pTprime, etax) * pTprime;  
  _jec = _jecWithL1V0;
  double pT_v0p = _Rjet(pTprime, etax) * pTprime;
  _jec = _jecDefault;
  //
  double l1n_p = _L1DataRaw(pT_v0p, etax); // V0Data
  double l1s_p = _L1Data(pT_p, etax); // V1Data = V1MC * V0Data / V0MC
  double sys_p = kfactor * (l1n_p / l1s_p - 1);
  //
  double l1nb_p = _L1MCRaw(pT_v0p, etax);
  double l1sb_p = _L1MC(pT_p, etax);
  double sysb_p = kfactor * (l1nb_p / l1sb_p - 1);

  _jec = _jecDefault;
  double pT_m = _Rjet(pTprime, -etax) * pTprime;
  _jec = _jecWithL1V0;
  double pT_v0m = _Rjet(pTprime, -etax) * pTprime;
  _jec = _jecDefault;    
  //
  double l1n_m = _L1DataRaw(pT_v0m, -etax);
  double l1s_m = _L1Data(pT_m, -etax);
  double sys_m = kfactor * (l1n_m / l1s_m - 1);
  //
  double l1nb_m = _L1MCRaw(pT_v0m, -etax);
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

  return sys * _ajet;// _ajet: jet area/ 0.5^2
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

  // This is MC truth flavor variation relative to QCD
  // It is not applied at the moment
  if (_errType & jec::kFlavorMC) {

    // Update for 38X, Mar14: Base new estimates on Pythia/Herwig
    // difference of flavor mapping relative to reference L2L3
    // Both agree rather well for gluons, but not quarks
    double dr = 0.;
    if (_calo) {
      
      //drawFragFlavor[4] = // uds (CALO)  
      double pl[4] = {0.9857, 0.0113, -0.938, 0.00};
      //drawFragFlavor[4] = // gluon (CALO)
      double pg[4] = {1.0081, -0.0025, 0.090, 0.00};

      TF1 *f1 = new TF1("f1","[0]+[1]*log10(0.01*x)+[2]/x+[3]/(x*x)",10,3000);
      f1->SetParameters(pl[0], pl[1], pl[2], pl[3]);
      double fl = fabs(1-f1->Eval(pTprime));
      f1->SetParameters(pg[0], pg[1], pg[2], pg[3]);
      double fg = fabs(1-f1->Eval(pTprime));
      delete f1;

      dr = max(fl,fg);
    }
    if (_pflow) {
      //drawFragFlavor[4] = // uds  
      double pl[4] = {0.9904, 0.0071, -0.202, 0.00};
      //drawFragFlavor[4] = // gluon  
      double pg[4] = {1.0035, 0.0029, 0.049, 0.00};
      //drawFragFlavor[4] = // bottom  
      double pb[4] = {0.9878, 0.0104, -0.202, 0.00};
      //drawFragFlavor[4] = // charm  
      double pc[4] = {0.9916, 0.0071, -0.202, 0.00};

      TF1 *f1 = new TF1("f1","[0]+[1]*log10(0.01*x)+[2]/x+[3]/(x*x)",10,3000);
      f1->SetParameters(pl[0], pl[1], pl[2], pl[3]);
      double fl = fabs(1-f1->Eval(pTprime));
      f1->SetParameters(pg[0], pg[1], pg[2], pg[3]);
      double fg = fabs(1-f1->Eval(pTprime));
      f1->SetParameters(pb[0], pb[1], pb[2], pb[3]);
      double fb = fabs(1-f1->Eval(pTprime));
      f1->SetParameters(pc[0], pc[1], pc[2], pc[3]);
      double fc = fabs(1-f1->Eval(pTprime));
      delete f1;

      dr = max(max(max(fl,fg),fb),fc);
    }

    if (_ideal) dr *= 0.5; // Pythia/Herwig mean for reference

    err2 += dr * dr;
    errFlavor = sqrt(err2);
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
  
    //return (_FlavorMix(pt, eta, fL, fG, fC, fB) - ref);
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


// Time-dependence uncertainty from L2
double JECUncertainty::_Time(const double eta) const {

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


// Random Cone offset for data
double JECUncertainty::_L1DataRaw(const double pT, const double eta) {

  //L1Offset VO Data
  assert(_jecL1nominal);
  _jecL1nominal->setRho(_Rho(_npv));
  _jecL1nominal->setJetA(TMath::Pi()*0.5*0.5*_ajet);
  _jecL1nominal->setNPV(_npv);
  _jecL1nominal->setJetEta(eta);
  _jecL1nominal->setJetE(pT*cosh(eta));
  _jecL1nominal->setJetPt(pT);
  return ( _jecL1nominal->getCorrection() );
}

// Random Cone offset for MC
double JECUncertainty::_L1MCRaw(const double pT, const double eta) {

 //L1Offset VO MC
  assert(_jecL1MCnominal);
  _jecL1MCnominal->setNPV(_npv);
  _jecL1MCnominal->setRho(_Rho(_npv));
  _jecL1MCnominal->setJetA(TMath::Pi()*0.5*0.5*_ajet);
  _jecL1MCnominal->setJetEta(eta);
  _jecL1MCnominal->setJetE(pT*cosh(eta));
  _jecL1MCnominal->setJetPt(pT);
  return ( _jecL1MCnominal->getCorrection() );
}

// Scaled offset for data (MC truth * RandomConeData / RandomConeMC)
double JECUncertainty::_L1Data(const double pT, const double eta) {

  //L1Offset V5 Data
  assert(_jecL1scaled);
  _jecL1scaled->setNPV(_npv);
  _jecL1scaled->setRho(_Rho(_npv));
  _jecL1scaled->setJetA(TMath::Pi()*0.5*0.5*_ajet);
  _jecL1scaled->setJetEta(eta);
  _jecL1scaled->setJetE(pT*cosh(eta));
  _jecL1scaled->setJetPt(pT);
  return ( _jecL1scaled->getCorrection() );
}

// Scaled offset for MC (MC truth)
double JECUncertainty::_L1MC(const double pT, const double eta) {

  assert(_jecL1pt);
  _jecL1pt->setNPV(_npv);
  _jecL1pt->setRho(_Rho(_npv));
  _jecL1pt->setJetA(TMath::Pi()*0.5*0.5*_ajet);
  _jecL1pt->setJetEta(eta);
  _jecL1pt->setJetE(pT*cosh(eta));
  _jecL1pt->setJetPt(pT);
  return ( _jecL1pt->getCorrection() );
}


// Offset density in events (needed to map average NPV to Rho for L1FastJet)
double JECUncertainty::_Rho(const double npvmean) {

  // For data 2012-10-24 (9.2/fb 53X Fall12 V1; range 0.5-40.5)
  //const float _uepf = 1.068; // from 2010 studies
  const float _ootpf = 1.606; // UE+OOT GeV/A 2012
  const float _itpf1 = 0.745; // IT p1 GeV/A 2012
  const float _itpf2 = +0.0016; // IT p2 GeV/A 2012
  double rho = _ootpf + (npvmean-1) * (_itpf1 + (npvmean-1) * _itpf2);

  return rho;
}
