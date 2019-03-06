#include "JECUncertainty.hpp"
#include "TMath.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"
#include "Math/BrentRootFinder.h"
//#include "Math/RootFinderAlgorithms.h"

#include <cmath>
#include <map>

// Autumn18 updates
// update MC-corrections --> Autumn18_V5_MC
// update data corrections --> take RunA for data corrections Autumn18_RunA_V5, but Autumn18_RunA_V5M2_DATA_L2L3Residual_AK4PFchs.txt for L2L3 residual (global fit result + dijet L2Res)
// [RelativeSample]:  RunD results for Dijet+barrel global fit vs. zjet only (no L2Res from dijets before)
// [AbsoluteSample]: use to illustrate ad-hoc difference between Friday's guesstimate and ZJet-only result for RunD


// Fall17_17Nov2017_V31 uncertainty files
// Reference for V31 "definition": https://www.dropbox.com/s/5up6c8g1vrau95v/20181025_JERC_2017_17NovV31_DATA.pdf?dl=0
// [Update central values]
// MC: Keep V24
// Data: Update to V31
// BCDEF CUSTOMHYBRID (const in 2.65-2.964; 3.839-5.191) Fall17_17Nov2017BCDEF_181022_MPF_CUSTOMHYBRID_181022_V31_L2L3Residual_pythia8_AK4PFchs.txt assembled from
// RunBCDEF_17Nov17_2017_METCorrectionCorrected3_woL1BXcleaning_JERnominalV2/output/Fall17_17Nov2017_MPF_FLAT_L2Residual_pythia8_AK4PFchs.txt + LOGLIN
// https://indico.cern.ch/event/767001/#49-l2res-with-jer-sfs-runbcdef
// and global fit  0.977 0.032 0.0
// https://www.dropbox.com/s/ee0kpg8zhohxlj2/jecslides_20181019_Fall17V28_good_p2.pdf?dl=0
// [TimeDependence]:
// - update above --> V31
// [Relative]
// BCDEF for flat vs. loglin MPF vs. ptbal
// [RelativeSample]
// - closure of Z+jet/gamma+jet on V31 vs. 1 (i.e. full correction in CUSTOMHYBRID style)
// [RelativeStat]
// - recycled
// [RelativeFSR]
// - recycled
// [RelativeJER]
// - recycled (from V27)
// [RelativePrefire]
// - sizeable effect at medium pt, L1BXClean vs. noL1BXClean not a good proxy
// MISSING:
// - more proper Prefire uncetainty (if any)
// - multijets
// - update to Flavor uncertainties (a la 2016)


// Fall17_17Nov2017_V27 uncertainty files
// [Update central values]
// Fall17_17Nov2017_V24_MC_L1RC and L1FastJet (_jecL1MCflat vs. _jecL1MCpt
// _jec updated to V27 Data
// Central L2L3Res values for V27 are from derivedWithDiJetTrigger_withoutL1SeedCleaning
// https://indico.cern.ch/event/759977/#30-update-of-l2residuals-for-a
// [Time Dependence]:
// - take BCDEF from https://indico.cern.ch/event/759977/contributions/3156284/attachments/1723118/2793691/L2ResCorrection_RunBCDEFcombined_SiTrg_withoutJERSF_withoutL1BXclean_V23_noClosure.tar
// - resuse old L3Fit 2-par [0.9787 0.0656]) which may be very outdated by now !?
// - use V27 for the IOVs;
// [Relative]
// - use https://indico.cern.ch/event/759977/contributions/3156284/attachments/1723118/2793691/L2ResCorrection_RunBCDEFcombined_SiTrg_withoutJERSF_withoutL1BXclean_V23_noClosure.tar
// -- for: flat vs. loglin; MPF vs. ptbal
// [RelativeSample]:
// - use Const/2par results of V27 closure results from Z+jet vs. dijet; Const for |eta|>=2.65
// [RelativeStat]
// - recycled
// [RelativeFSR]
// - recycled
// [RelativeJER]
// - taken from closure: https://indico.cern.ch/event/759977/contributions/3157465/attachments/1724807/2792149/L2ResCorrection_RunBCDEFcombined_SiTrg_JERSF2016_withAndWithoutL1BXclean_V27_closure.tar
// - use up/down from "without L1BXcleaning"
// - change paths and use updated files in CondFormats/JetMETObjects/data
// NEW: [RelativePrefire]
//  - taken from closure: https://indico.cern.ch/event/759977/contributions/3157465/attachments/1724807/2792149/L2ResCorrection_RunBCDEFcombined_SiTrg_JERSF2016_withAndWithoutL1BXclean_V27_closure.tar
//  - use full difference between with/without L1BXPrefiring JERnominal



using namespace std;

// Future improvements:
// 1) Add L2L3(+Res) to InitL1, so that JEC inversion really gets pTprime
//    Currently the L1 uncertainty could be evaluated at slightly wrong point,
//    but difference is small because nom and denom have the same

// Turn debug mode on if the code fails with an exception from JEC packages
// This is most likely a missing/misnamed file given to JetCorrectorParameters
// The last file printed out before the crash in debug mode is usually the fault
bool debug = false;//true;//false;


JECUncertainty::JECUncertainty(const jec::JetAlgo& algo, 
			       const jec::DataType& type, 
			       const jec::ErrorTypes& errType,
			       const double mu) :
  _algo(algo), _type(type), _errType(errType), _mu(mu)
{

  _fjes = 0; _emat = 0; _fhb = _fl1 = 0;
  //_fl3ref = _fl3up = _fl3dw = _fl2up = 0;
  _fl3 =_fl2 = 0;
  _hl3 = _hl3ref = _hl2ref = 0;

  _algo = algo;
  _calo = (_algo==jec::AK4CALO || _algo==jec::AK5CALO ||
	   _algo==jec::AK7CALO || _algo==jec::AK8CALO);
  _jpt = (_algo==jec::AK4JPT);
  _pfchs = (_algo==jec::AK4PFchs || _algo==jec::AK5PFchs ||
	    _algo==jec::AK7PFchs || _algo==jec::AK8PFchs);
  _pflow = (_algo==jec::AK4PF || _algo==jec::AK5PF ||
	    _algo==jec::AK7PF || _algo==jec::AK8PF || _pfchs);
  _ideal = false;//(_algo==IDEAL);
  _trkbase = (_pflow || _jpt);

  _ajet = 1;
  if (algo==jec::AK4PF || algo==jec::AK4PFchs || algo==jec::AK4CALO)
    _ajet = pow(0.4/0.5,2);
  if (algo==jec::AK7PF || algo==jec::AK7PFchs || algo==jec::AK7CALO)
    _ajet = pow(0.7/0.5,2);
  if (algo==jec::AK8PF || algo==jec::AK8PFchs || algo==jec::AK8CALO)
    _ajet = pow(0.8/0.5,2);

  _InitL1();
  _InitJEC();
  _InitL2Res();
  _InitL3Res();
  
  // initialize Run I uncertainty
  {
    const char *cd = "CondFormats/JetMETObjects/data";
    const char *s = Form("%s/Winter14_V8_DATA_UncertaintySources_AK5PFchs.txt",cd); // V8 Run I (official file, but same as above)
    const char *s2 = "TotalNoFlavorNoTime";
    //cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *pref = new JetCorrectorParameters(s,s2);
    _uncRunI = new JetCorrectionUncertainty(*pref);
  }

}


// Uncertainty is returned as _relative_ uncertainty
// Keep systematics signed for correlations
double JECUncertainty::Uncert(const double pTprime, const double eta) {

  if (debug) cout << "JECUncertainty::Uncert("<<pTprime
		  << ", "<<eta<<")" << endl << flush;
  
  double err2 = 0;
  double eta2 = min(max(eta,-5.190),5.190); // fix for drawing macro

  // 2013 systematics
  double errAbs(0), errRel(0), errPileUp(0), errFlavor(0), errTime(0);
  double errRunI(0);
  if (_errType & jec::kAbsolute) {
    if (debug) cout << "_Absolute" << endl << flush;
    errAbs = _Absolute(pTprime);
    err2 += errAbs * errAbs;
  }
  if (_errType & jec::kRelative) {
    if (debug) cout << "_Relative" << endl << flush;
    errRel = _Relative(pTprime, eta2);
    err2 += errRel * errRel;
  }
  if (_errType & (jec::kPileUp | jec::kPileUpMuZero | jec::kPileUpEnvelope)) {
    if (debug) cout << "_PileUp" << endl << flush;
    errPileUp = _PileUp(pTprime, eta2);
    err2 += errPileUp * errPileUp;
  }
  if (_errType & jec::kFlavorMask) {
    if (debug) cout << "_Flavor" << endl << flush;
    errFlavor = _Flavor(pTprime, eta);
    err2 += errFlavor*errFlavor;
  }
  if (_errType & (jec::kTime | jec::kTimePtEtaMask)) {
    if (debug) cout << "_Time" << endl << flush;
    errTime = _Time(pTprime, eta2);
    err2 += errTime * errTime;
  }
  if (_errType & jec::kRunI) {
    _uncRunI->setJetPt(pTprime);
    _uncRunI->setJetEta(eta);
    errRunI = _uncRunI->getUncertainty(true);
    err2 += errRunI*errRunI;
  }
  
  double err = sqrt(err2);

  // if requesting single source, return signed for sign-changing cases
  if (!(_errType & ~jec::kAbsoluteFrag))   return errAbs;
  if (!(_errType & ~jec::kAbsoluteSPRE))   return errAbs;
  if (!(_errType & ~jec::kAbsoluteSPRH))   return errAbs;
  if (!(_errType & ~jec::kAbsoluteSample)) return errAbs;
  if (!(_errType & ~jec::kRelativeFSR))   return errRel;
  if (!(_errType & ~jec::kRelativePtBB))  return errRel;
  if (!(_errType & ~jec::kRelativePtEC1)) return errRel;
  if (!(_errType & ~jec::kRelativePtEC2)) return errRel;
  if (!(_errType & ~jec::kRelativePtHF))  return errRel;
  if (!(_errType & ~jec::kRelativePt))    return errRel; // EXTRA
  if (!(_errType & ~jec::kRelativeBal))    return errRel;  // Run II
  if (!(_errType & ~jec::kRelativeSample)) return errRel;  // Run II
  if (!(_errType & ~jec::kPileUpDataMC)) return errPileUp;
  if (!(_errType & ~jec::kPileUpPtRef))  return errPileUp;
  if (!(_errType & ~jec::kPileUpPtBB))   return errPileUp;
  if (!(_errType & ~jec::kPileUpPtEC1))  return errPileUp;
  if (!(_errType & ~jec::kPileUpPtEC2))  return errPileUp;
  if (!(_errType & ~jec::kPileUpPtHF))   return errPileUp;
  if (!(_errType & ~jec::kPileUpPt))     return errPileUp; // EXTRA
  if (!(_errType & ~jec::kPileUpPtEta))  return errPileUp; // EXTRA
  if (!(_errType & ~jec::kPileUpMuZero)) return errPileUp; // OPT
  if (!(_errType & ~jec::kPileUpEnvelope)) return errPileUp; // OPT
  if (!(_errType & ~jec::kTimeRunB)) return errTime;
  if (!(_errType & ~jec::kTimeRunC)) return errTime;
  if (!(_errType & ~jec::kTimeRunDE)) return errTime;
  //  if (!(_errType & ~jec::kTimeRunE)) return errTime;
  if (!(_errType & ~jec::kTimeRunF)) return errTime;
  if (!(_errType & ~jec::kTimePtEta))     return errTime;
  //
  if (!(_errType & ~jec::kFlavorMask))   return errFlavor;

  return err;
} // Uncert

void JECUncertainty::_InitL1() {

  // Inputs taken:
  // - regular L1FastJet files
  // - random cone (RC) files
  //   (used for type-I MET and scaling L1FastJet_DATA)
  // Both types are now part of the JEC database file by default

  map<jec::JetAlgo, const char*> names;
  names[jec::AK4PF] = "AK4PF";
  names[jec::AK4PFchs] = "AK4PFchs";
  names[jec::AK4PFpuppi] = "AK4PFpuppi";
  names[jec::AK4CALO] = "AK4Calo";
  names[jec::AK8PF] = "AK8PF";
  names[jec::AK8PFchs] = "AK8PFchs";
  names[jec::AK8PFpuppi] = "AK8PFpuppi";
  //names[jec::AK8CALO] = "AK8Calo";
  const char *a = names[_algo];
  string directory = "CondFormats/JetMETObjects/data/";
  const char *d = directory.c_str();

  {
    //const char *s = Form("%sSummer16_03Feb2017_V1_MC_L1RC_%s.txt",d,a);
    const char *s = Form("%sAutumn18_V5_MC_L1RC_%s.txt",d,a);
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l1);
    _jecL1MCflat = new FactorizedJetCorrector(v);
  }
  {
    //const char *s = Form("%sSummer16_03Feb2017_V1_MC_L1FastJet_%s.txt",d,a);
    const char *s = Form("%sAutumn18_V5_MC_L1FastJet_%s.txt",d,a);
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l1);
    _jecL1MCpt = new FactorizedJetCorrector(v);
  }

  // For PileUpPtRef (L3Res only using AK4PFchs)
  {
    const char *a = "AK4PFchs"; // !! L3Res only for this
    //const char *s = Form("%sSummer16_03Feb2017_V1_MC_L1RC_%s.txt",d,a);
    const char *s = Form("%sAutumn18_V5_MC_L1RC_%s.txt",d,a);

    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l1);
    _jecL1MCflat_ak4pfchs = new FactorizedJetCorrector(v);
  }
  {
    const char *a = "AK4PFchs"; // !! L3Res only for this
    //const char *s = Form("%sSummer16_03Feb2017_V1_MC_L1FastJet_%s.txt",d,a);
    const char *s = Form("%sAutumn18_V5_MC_L1FastJet_%s.txt",d,a);
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l1);
    _jecL1MCpt_ak4pfchs = new FactorizedJetCorrector(v);
  }

} // InitL1


void JECUncertainty::_InitJEC() {

  // Inputs
  // - regular L1L2L3+Res files
  // All files are stored in JEC database:
  // https://github.com/cms-jet/JECDatabase/tree/master/tarballs)

  map<jec::JetAlgo, const char*> names;
  names[jec::AK4PF] = "AK4PF";
  names[jec::AK4PFchs] = "AK4PFchs";
  names[jec::AK4PFpuppi] = "AK4PFpuppi";
  names[jec::AK4CALO] = "AK4Calo";
  names[jec::AK8PF] = "AK8PF";
  names[jec::AK8PFchs] = "AK8PFchs";
  names[jec::AK8PFpuppi] = "AK8PFpuppi";
  //names[jec::AK8CALO] = "AK8Calo";
  const char *a = names[_algo];
  string directory = "CondFormats/JetMETObjects/data/";
  const char *d = directory.c_str();

  const char *s;
  s = Form("%sAutumn18_RunA_V5_DATA_L1FastJet_%s.txt",d,a);
  if (debug) cout << s << endl << flush;
  JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
  s = Form("%sAutumn18_RunA_V5_DATA_L2Relative_%s.txt",d,a);
  if (debug) cout << s << endl << flush;
  JetCorrectorParameters *l2 = new JetCorrectorParameters(s);
  s = Form("%sAutumn18_RunA_V5_DATA_L3Absolute_%s.txt",d,a);
  if (debug) cout << s << endl << flush;
  JetCorrectorParameters *l3 = new JetCorrectorParameters(s);
  // Only one L3Residual derived for now for AK4PFchs
  // (although we clone this later on)
  s = Form("%sAutumn18_RunA_V5M2_DATA_L2L3Residual_%s.txt",d,a);
  if (debug) cout << s << endl << flush;
  JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);

  vector<JetCorrectorParameters> v;
  v.push_back(*l1);
  v.push_back(*l2);
  v.push_back(*l3);
  v.push_back(*l2l3res);
  _jecDefault = new FactorizedJetCorrector(v);
  _jec = _jecDefault;

  { // RunB 
    s = Form("%sFall17_17Nov2017B_V31_DATA_L2L3Residual_%s.txt",d,a);
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);
    
    vector<JetCorrectorParameters> v;
    v.push_back(*l1);
    v.push_back(*l2);
    v.push_back(*l3);
    v.push_back(*l2l3res);
    _jecB = new FactorizedJetCorrector(v);
  }
  { // RunC
    s = Form("%sFall17_17Nov2017C_V31_DATA_L2L3Residual_%s.txt",d,a);
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);
    
    vector<JetCorrectorParameters> v;
    v.push_back(*l1);
    v.push_back(*l2);
    v.push_back(*l3);
    v.push_back(*l2l3res);
    _jecC = new FactorizedJetCorrector(v);
  }
  { // RunDE
    s = Form("%sFall17_17Nov2017DE_V31_DATA_L2L3Residual_%s.txt",d,a);
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);
    
    vector<JetCorrectorParameters> v;
    v.push_back(*l1);
    v.push_back(*l2);
    v.push_back(*l3);
    v.push_back(*l2l3res);
    _jecDE = new FactorizedJetCorrector(v);
  }
//  { // RunE
//    s = Form("%sFall17_17Nov2017E_V31_DATA_L2L3Residual_%s.txt",d,a);
//    if (debug) cout << s << endl << flush;
//    JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);
//    
//    vector<JetCorrectorParameters> v;
//    v.push_back(*l1);
//    v.push_back(*l2);
//    v.push_back(*l3);
//    v.push_back(*l2l3res);
//    _jecE = new FactorizedJetCorrector(v);
//  }
  { // RunF
    s = Form("%sFall17_17Nov2017F_V31_DATA_L2L3Residual_%s.txt",d,a);
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);
    
    vector<JetCorrectorParameters> v;
    v.push_back(*l1);
    v.push_back(*l2);
    v.push_back(*l3);
    v.push_back(*l2l3res);
    _jecF = new FactorizedJetCorrector(v);
  }
  { // RunBCDEF (all 2017)
    s = Form("%sFall17_17Nov2017BCDEF_181022_MPF_CUSTOMHYBRID_181022_V31_L2L3Residual_pythia8_%s.txt",d,a); // see above
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);
    
    vector<JetCorrectorParameters> v;
    v.push_back(*l1);
    v.push_back(*l2);
    v.push_back(*l3);
    v.push_back(*l2l3res);
    _jecBCDEF = new FactorizedJetCorrector(v);
  }

  // Special JEC with RC correction
  /*
  {
    //s = Form("%sFall15_25nsV1_DATA_L1RC_%s.txt",d,a); // 76X
    s = Form("%sSpring16_25nsV3_DATA_L1RC_%s.txt",d,a); // 80X V3
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l1rc = new JetCorrectorParameters(s);
    
    vector<JetCorrectorParameters> rc;
    rc.push_back(*l1rc);
    rc.push_back(*l2);
    rc.push_back(*l3);
    rc.push_back(*l2l3res);
    _jecWithL1RC = new FactorizedJetCorrector(rc);
  }
  */

} // InitJEC

void JECUncertainty::_InitL2Res() {

  // Inputs (+ updated, - recycled):
  // + Pythia MPF loglin (=default)
  // + Pythia MPF flat (vs MPF loglin: RelativePt)
  // + Pythia MPF loglin (vs pT loglin: RelativeBal)
  // + L2L3Res closure test for Z+jet vs. dijet on top of V27 (RelativeSample)
  // + Pythia JERup vs JERdw (RelativeJER; MPF loglin)
  // + NEW: L1BX filter with vs. without (MPF loglin) (RelativePrefire)
  // - Herwig MPF loglin (vs Pythia MPF loglin: RelativeFSR)  ; recycled
  // - Pythia MPF loglin .STAT (RelativeStat)                 ; recycled
  // These come in variants AK4PFchs, AK4PFpuppi, AK8PFchs, AK8puppi

  map<jec::JetAlgo, const char*> names;
  names[jec::AK4PF] = "AK4PFchs"; // Replace "AK4PF";
  names[jec::AK4PFchs] = "AK4PFchs";
  names[jec::AK4PFpuppi] = "AK4PFpuppi";
  names[jec::AK4CALO] = "AK4PFchs"; // Replace "AK4Calo";
  names[jec::AK8PF] = "AK4PFchs"; // NB: AK4 // Replace "AK4PF";
  names[jec::AK8PFchs] = "AK4PFchs"; // NB: AK4
  names[jec::AK8PFpuppi] = "AK8PFpuppi";
  //names[jec::AK8CALO] = "AK8KFchs"; // Replace "AK8Calo";
  string directory = "CondFormats/JetMETObjects/data/";
  const char *d = directory.c_str();
  const char *a = names[_algo];

  const char *s, *s2;
  // For RelativePt (flat vs loglin)
  {
    s = Form("%sAutumn18_RunD_V5_20190304_DATA_MPF_FLAT_L2Residual_%s.txt",d,a);
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l2l3res);
    _jecL2ResFlat = new FactorizedJetCorrector(v);
  }
  {
    s = Form("%sAutumn18_RunD_V5_20190304_DATA_MPF_LOGLIN_L2Residual_%s.txt",d,a);
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l2l3res);
    _jecL2ResPt = new FactorizedJetCorrector(v);
  }

  // For RelativeBal (MPF vs pTbal); separate source for Sum16
  {
    s = Form("%sAutumn18_RunD_V5_20190304_DATA_MPF_LOGLIN_L2Residual_%s.txt",d,a);
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l2l3res);
    _jecL2ResMPF = new FactorizedJetCorrector(v);
  }
  {
    s = Form("%sAutumn18_RunD_V5_20190304_DATA_pT_LOGLIN_L2Residual_%s.txt",d,a);
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l2l3res);
    _jecL2ResBal = new FactorizedJetCorrector(v);
  }

  // For RelativeSample (dijet vs Z/gam+jet)
  {
    s = Form("%sV31_CUSTOMHYBRID_BCDEF_CollectL2Output_gam_zll_PtBalMPF_Fall17_17Nov2017BCDEF_VXXX_DATA_L2L3Residual_%s.txt",d,a); // Closure result after V31 custom hybrid
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l2l3res);
    _jecL2L3ResCustomHybridZGamJet = new FactorizedJetCorrector(v);
  }

  {
    s = Form("%sAutumn18_RunD_V5M_DATA_L2L3Residual_%s.txt",d,a); // RunD result for Dijet+guesstimate for barrel global fit (possibly accounting for lack of gamma+jet, pulling down at high pt)
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l2l3res);
    _jecL2L3Res2ParDiJetGuesstimate = new FactorizedJetCorrector(v);
  }
  {
    s = Form("%sAutumn18_RunD_V5M2_DATA_L2L3Residual_%s.txt",d,a); // RunD results for Dijet+barrel global fit vs. zjet only (no L2Res from dijets before)
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l2l3res);
    _jecL2L3Res2ParDiJet = new FactorizedJetCorrector(v);
  }
  {
    s = Form("%sAutumn18_RunD_V5MZJetOnly2Par_DATA_L2L3Residual_%s.txt",d,a); // RunD results for Dijet+barrel global fit vs. zjet only (no L2Res from dijets before)
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l2l3res);
    _jecL2L3Res2ParZGamJet = new FactorizedJetCorrector(v);
  }
  {
    s = Form("%sAutumn18_RunD_V5M2_FLAT_DATA_L2L3Residual_%s.txt",d,a); // RunD results for Dijet+barrel global fit vs. zjet only (no L2Res from dijets before)
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l2l3res);
    _jecL2L3ResConstDiJet = new FactorizedJetCorrector(v);
  }
  {
    s = Form("%sAutumn18_RunD_V5MZJetOnly1Par_DATA_L2L3Residual_%s.txt",d,a); // RunD results for Dijet+barrel global fit vs. zjet only (no L2Res from dijets before)
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l2l3res);
    _jecL2L3ResConstZGamJet = new FactorizedJetCorrector(v);
  }
  //{
  //s = Form("%s/Summer16_03Feb2017V7BCDEFGH_GlobalFit_ZJet_GammaJet_DiJet_2_remove_lowpt_ZJet_L2L3Residual_pythia8_%s.txt",d,a);
  //if (debug) cout << s << endl << flush;
  //JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);
  //vector<JetCorrectorParameters> v;
  //v.push_back(*l2l3res);
  //_jecL2L3ResZGamDiJet = new FactorizedJetCorrector(v);
  //}

  {
    s = Form("%sRunBCDEF_17Nov17_2017_withoutL1BXcleaning_JERnominal_output_Fall17_17Nov2017_MPF_LOGLIN_L2Residual_pythia8_%s.txt",d,a); // Closure result after V27
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l2l3res);
    _jecL2ResPrefireFilterNo = new FactorizedJetCorrector(v);
  }

  {
    s = Form("%sRunBCDEF_17Nov17_2017_withL1BXcleaning_JERnominal_output_Fall17_17Nov2017_MPF_LOGLIN_L2Residual_pythia8_%s.txt",d,a); // Closure result after V27
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l2l3res);
    _jecL2ResPrefireFilterYes = new FactorizedJetCorrector(v);
  }
  
  
  // RelativeFSR (Pythia vs Herwig) ; recycled
  {
    s = Form("%sSpring16_25ns_MPF_LOGLIN_L2Residual_pythia8_v4_%s.txt",d,a); // 80XV8 -- no update for Sum16V2, Sum16V3, 03FebV3, 03FebV7, Fall17, Fall17V27
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l2l3res);
    _jecL2ResPY = new FactorizedJetCorrector(v);
  }
  {
    s = Form("%sSpring16_25ns_MPF_LOGLIN_L2Residual_herwigpp_v4_%s.txt",d,a); // 80XV8 -- no update for Sum16V2, Sum16V3, 03FebV3, 03FebV7, Fall17, Fall17V27
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l2l3res);
    _jecL2ResHW = new FactorizedJetCorrector(v);
  }

  // For RelativeStat (true statistical only) ; recycled
  {
    s = Form("%sSummer16_03Feb2017BCDEFGHV7_MPF_LOGLIN_L2Residual_pythia8_%s.txt.STAT",d,a); // no update for Fall17V27
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l2l3res);
    _jecL2ResStat = new FactorizedJetCorrector(v);
  }

  // For RelativeJER
  {
    //s = Form("%sSummer16_03Feb2017BCDEFGHV7_MPF_LOGLIN_JERup_L2Residual_pythia8_%s.txt",d,a);
    s = Form("%sRunBCDEF_17Nov17_2017_withoutL1BXcleaning_JERup_output_Fall17_17Nov2017_MPF_LOGLIN_L2Residual_pythia8_%s.txt",d,a);
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l2l3res);
    _jecL2jerup = new FactorizedJetCorrector(v);
  }
  {
    //s = Form("%sSummer16_03Feb2017BCDEFGHV7_MPF_LOGLIN_JERdown_L2Residual_pythia8_%s.txt",d,a);
    s = Form("%sRunBCDEF_17Nov17_2017_withoutL1BXcleaning_JERdown_output_Fall17_17Nov2017_MPF_LOGLIN_L2Residual_pythia8_%s.txt",d,a);
    if (debug) cout << s << endl << flush;
    JetCorrectorParameters *l2l3res = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l2l3res);
    _jecL2jerdw = new FactorizedJetCorrector(v);
  }
  
} // InitL2Res

void JECUncertainty::_InitL3Res() {

  // 03FebV3 BCDEFGH fit 0.9911 +/- 0.0017, 0.0437 +/- 0.0029,
  // p0 and p1 correlation 0.10, seems reasonable for pTref=208 GeV
  /*
  const int n = 2;
  const double pars[n] =
    { 0.9911, +0.0437};
  const double emata[n][n] =
    {{2.819e-06,  4.872e-07},
     {4.872e-07,  8.488e-06}};
  */

  // Autumn18_RunA_V5 2-par fit with Zll + gamjet, chi2/NDF = 20.6466 / 29
  // Fit parameters:
  // 1:    0.9932 +/- 0.0027,    2:    0.0659 +/- 0.0112,
  // p0-p1 correlation 0.17
  const int n = 2;
  const double pars[n] =
    { 0.9932, 0.0659};
  const double emata[n][n] =
    {{ 7.361e-06,  5.126e-06},
     { 5.126e-06,  0.0001256}};

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
    _fjes = new TF1("fjes",_jesfit,10.,6500.,n);
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
  if (mu<0) mu = _mu;
  double npv = _NpvFromMu(mu);
  double rho = _RhoFromMu(mu);
  if (!jec) jec = _jec;

  ResponseFunc f(pTprime,jec,npv,eta,rho,ajet);
  
  //ROOT::Math::Roots::Brent brf;
  ROOT::Math::BrentRootFinder brf;
  //brf.SetLogScan(true);
  // Set parameters of the method
  //brf.SetFunction(f,std::max(2.0,0.25*pTprime),std::min(4*pTprime,4000.0));
  //bool found_root = brf.Solve(50,1e-4,1e-5);
 // Enlarge upper range and improve precision to cover also AK8
  //brf.SetFunction(f,std::max(2.0,0.25*pTprime),std::min(10*pTprime,4000.0));
  brf.SetFunction(f,std::max(2.0,0.25*pTprime),std::min(10*pTprime,6500*1.2));
  bool found_root = brf.Solve(50,1e-4,1e-5);
  double pTraw = brf.Root();
  double rjet = pTraw /pTprime;
  assert(jec);
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
double JECUncertainty::_Absolute(const double pt){

  double stat          = (_errType & jec::kAbsoluteStat          ? _AbsoluteStat(pt)          : 0.);
  double scale         = (_errType & jec::kAbsoluteScale         ? _AbsoluteScale()         : 0.);
  double sample        = (_errType & jec::kAbsoluteSample        ? _AbsoluteSample(pt)      : 0.);
  double FlMap         = (_errType & jec::kAbsoluteFlavorMapping ? _AbsoluteFlavorMapping() : 0.); //for backward-compatibility/historical reasons
  double MPFBias       = (_errType & jec::kAbsoluteMPFBias       ? _AbsoluteMPFBias()       : 0.);
  double spr           = (_errType & jec::kAbsoluteSPR           ? _AbsoluteSPR(pt)         : 0.);
  double frag          = (_errType & jec::kAbsoluteFrag          ? _AbsoluteFrag(pt)        : 0.);

  // signed sources
  if (!(_errType & ~jec::kAbsoluteStat)) return stat;
  if (!(_errType & ~jec::kAbsoluteSPRE)) return spr;
  if (!(_errType & ~jec::kAbsoluteSPRH)) return spr;
  if (!(_errType & ~jec::kAbsoluteFrag)) return frag;
  if (!(_errType & ~jec::kAbsoluteSample)) return sample;

  return sqrt(stat*stat + scale*scale + FlMap*FlMap + MPFBias*MPFBias + spr*spr + frag*frag + sample*sample);
} // Absolute


// pT dependent fit of L3Res
TF1 *fhb(0);
Double_t _jesfit(Double_t *x, Double_t *p) {

  double pt = x[0];

  if (!fhb) fhb = new TF1("fhb","max(0.,[0]+[1]*pow(x,[2]))",10,6500);
  fhb->SetParameters(1.03091e+00, -5.11540e-02, -1.54227e-01); // SPRH
  // p[0]: overall scale shift, p[1]: HCAL shift in % (full band +3%)
  double ptref = 208; // pT that minimizes correlation in p[0] and p[1]
  double jes = (p[0] + p[1]/3.*100*(fhb->Eval(pt)-fhb->Eval(ptref)));

  return jes;
} // jesfit

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
  // Same for SPRH, in which case AbsStatSys is probably very close to zero
  double AbsStat = _jesfitunc(pTprime, _fjes, _emat);
  double AbsScale = _AbsoluteScale();
  double AbsSPRH = _AbsoluteSPRH(pTprime);
  double AbsStatSys = sqrt(max(AbsStat*AbsStat - AbsScale*AbsScale
  			       - AbsSPRH*AbsSPRH, 0.));

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

//  //double AbsScaleSys = 0.0017; // 03FebV3 (mu=e=0.2%, g=2.0%(!!), EMFoot=0.5%)
  double AbsScaleSys = 0.0010; // 17NovV6 (Zll=0.05%, g=0.5%, EMFoot=0.5%)
  return AbsScaleSys;
}

double JECUncertainty::_AbsoluteSample(const double pTprime) {
  //Temporary hack for illustration: Difference between guesstimate and global fit with only Z+Jet

  // Take Z/gam+jet-only closure fit as extra source
  // Nominal correction from dijet, so dijet should close
  assert(_jecL2L3Res2ParDiJet);
  assert(_jecL2L3Res2ParDiJetGuesstimate);
  
  double rfit = _Rjet(pTprime, 0., -1, -1, _jecL2L3Res2ParDiJet);
  double rg = _Rjet(pTprime, 0., -1, -1, _jecL2L3Res2ParDiJetGuesstimate);
  
  double kfactor = 1;
  double err = kfactor*(rg / rfit - 1);
  //std::cout << "_AbsoluteSample: " << err << std::endl;
  return err;
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

    TF1 *f1 = new TF1("f1","[0]+[1]*TMath::Log10(0.01*x)+[2]/x+[3]/(x*x)",10,3000);
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

  if (_errType & jec::kAbsoluteSPRE) { // SPR in ECAL
    errSPRE = _AbsoluteSPRE(pTprime);
  }
  if (_errType & jec::kAbsoluteSPRH) { // SPR in HCAL
    errSPRH = _AbsoluteSPRH(pTprime);
  }

  // replace directly done SPR with pieces broken up into ECAL and HCAL
  errSPR = sqrt(pow(errSPRE,2) + pow(errSPRH,2));

  // signed sources
  if (!(_errType & ~jec::kAbsoluteSPRE)) return errSPRE;
  if (!(_errType & ~jec::kAbsoluteSPRH)) return errSPRH;

  return errSPR;
} // AbsoluteSPR

// Single pion response in HCAL
double JECUncertainty::_AbsoluteSPRH(const double pTprime) const {

  // New fits from Juska on Nov 26, 2012, 5:26 pm
  double errSPRH(0), difSPRH(0), refSPRH(0);
  double refpt = 208.; // 2012 RDMC V3PT

  TF1 *f = new TF1("fsprh","TMath::Max(0.,[0]+[1]*pow(x,[2]))",10,6500);

  f->SetParameters(1.03091e+00, -5.11540e-02, -1.54227e-01);
  if (_calo) f->SetParameters(1.02246e+00, -1.55689e-02, -1.17219e-01);
  // 2016BCDEFGH fit
  //difSPRH = 0.0029/0.03 * (f->Eval(pTprime)-1) + 1; // 03FebV3
  //refSPRH = 0.0029/0.03 * (f->Eval(refpt)-1) + 1; // 03FebV3
  // 2017BCDEF fit
  difSPRH = 0.0084/0.03 * (f->Eval(pTprime)-1) + 1; // 17NovV6
  refSPRH = 0.0084/0.03 * (f->Eval(refpt)-1) + 1; // 17NovV6
  //
  errSPRH = (difSPRH-refSPRH);

  // NB: returns signed systematic
  return errSPRH;
} // AbsoluteSPRH

// Single pion response in ECAL
double JECUncertainty::_AbsoluteSPRE(const double pTprime) const {

  // New fits from Juska on Nov 26, 2012, 5:26 pm
  double errSPRE(0), difSPRE(0), refSPRE(0);
  double refpt = 208.; // 2012 RDMC V3PT

  TF1 *f = new TF1("fspre","TMath::Max(0.,[0]+[1]*pow(x,[2]))",10,3500);


  // Fix 2014-05-21: multiply errSPRE and errSRPH residuals by sqrt(2)
  // to keep errSPRE_3%(oplus)errSPRH_3% ~ errSPR_3%
  f->SetParameters(1.00567e+00, -3.04275e-02, -6.75493e-01);
  if (_calo) f->SetParameters(1.00166e+00, 1.57065e-02, -2.06585e-01);
  // 80X V4M3: increase pT dependence uncertainty 50% over Run I
  //difSPRE = 2.0*sqrt(2.)*(f->Eval(pTprime)-1) + 1;
  //refSPRE = 2.0*sqrt(2.)*(f->Eval(refpt)-1) + 1;
  // 80XV8: back to Run I a priori uncertainty
  difSPRE = sqrt(2.)*(f->Eval(pTprime)-1) + 1;
  refSPRE = sqrt(2.)*(f->Eval(refpt)-1) + 1;
  errSPRE = (difSPRE-refSPRE);

  // NB: returns signed systematic
  return errSPRE;
} // AbsoluteSPR

// Combine relative uncertainties
double JECUncertainty::_Relative(const double pTprime,
				 const double eta) {
  
  double sjer = (_errType & jec::kRelativeJER ? _RelativeJER(pTprime, eta) : 0.);
  double sfsr = (_errType & jec::kRelativeFSR ? _RelativeFSR(pTprime, eta) : 0.);
  double stat = (_errType & jec::kRelativeStat ? _RelativeStat(pTprime, eta) : 0.);
  double spt =  (_errType & jec::kRelativePt ? _RelativePt(pTprime, eta) : 0.);
  double sbal = (_errType & jec::kRelativeBal ? _RelativeBal(pTprime, eta) : 0.);
  double spf =  0.;//(_errType & jec::kRelativePrefire ? _RelativePrefire(pTprime, eta) : 0.);
  double ssam = (_errType & jec::kRelativeSample ? _RelativeSample(pTprime, eta) : 0.);

  // signed sources
  if (!(_errType & ~jec::kRelativeFSR)) return sfsr;
  if (!(_errType & ~jec::kRelativePtBB)) return spt;
  if (!(_errType & ~jec::kRelativePtEC1)) return spt;
  if (!(_errType & ~jec::kRelativePtEC2)) return spt;
  if (!(_errType & ~jec::kRelativePtHF))  return spt;
  if (!(_errType & ~jec::kRelativePt))  return spt; // XTRA
  if (!(_errType & ~jec::kRelativeBal))  return sbal; // Sum16
  if (!(_errType & ~jec::kRelativePrefire))  return spf; // Sum16
  if (!(_errType & ~jec::kRelativeSample))  return ssam; // 03FebV7

  return sqrt(sjer*sjer + sfsr*sfsr + stat*stat + spt*spt
	      + sbal*sbal + spf*spf + ssam*ssam);
}


// Relative scale uncertainty vs eta from JER bias
double JECUncertainty::_RelativeJER(const double pTprime,
				    const double eta) {

  assert(_jecL2jerup);
  double rup = _Rjet(pTprime, eta, -1, -1, _jecL2jerup);

  assert(_jecL2jerdw);
  double rdw = _Rjet(pTprime, eta, -1, -1, _jecL2jerdw);


  // Uncertainty is half of up and down, so full difference to mean
  double rmean = 0.5 * (rup + rdw);
  double err = 0.5 * (rup - rdw) / rmean;

  double x = fabs(eta);
  if (x<1.5) return 0; // Assuming BB negligible for now
  if ( (x>=1.5 && x<2.5 && (_errType & jec::kRelativeJEREC1)) ||
       (x>=2.5 && x<3.0 && (_errType & jec::kRelativeJEREC2)) ||
       (x>=3.0 && x<5.2 && (_errType & jec::kRelativeJERHF)) )
    return err;

  return 0;
} // RelativeJER


// Helper function for _RelativeFSR
TF1 *fkFSR(0);
Double_t _kFSR(Double_t *xx, Double_t *p) {

  double x = *xx;

  const int np = 3;
  double p0 = p[0];
  double p1 = p[1];
  double p2 = p[2];

  double val = p0 + p1*cosh(x) / (1 + p2*cosh(x));

  double df[np] = {1, cosh(x) / (1 + p2*cosh(x)),
		   -p1 * pow(cosh(x) / (1 + p2*cosh(x)), 2)}; 

  double em[np][np] = {{p[3], (p[4]), (p[6])},
		       {p[4],  p[5],  (p[7])},
		       {p[6],  p[7],   p[8]}};

  double sigma = p[9];

  double err2(0);
  for (int i = 0; i != np; ++i) {
    for (int j = 0; j != np; ++j) {
      err2 += df[i]*df[j]*em[i][j];
    }
  }
  double err = sqrt(err2);
  
  return (val + sigma * err);
}

// Relative scale uncertainty vs eta from soft radiation (FSR)
// Estimated from difference between MPF and pT balance methods
double JECUncertainty::_RelativeFSR(const double pTprime, const double eta) {

  assert(_jecL2ResPY);
  double rpy = _Rjet(pTprime, eta, -1, -1, _jecL2ResPY);
  assert(_jecL2ResHW);
  double rhw = _Rjet(pTprime, eta, -1, -1, _jecL2ResHW);
  
  double kfactor = 1.0;
  double diff = kfactor * (rhw / rpy - 1);

  return diff;
} // RelativeFSR


// Statistical uncertainty in L2Res (symmetrized, wide bins)
double JECUncertainty::_RelativeStat(const double pTprime,
				     const double eta) const {

  double err(0);
  double x = fabs(eta);
  if ( (x>=1.3 && x<2.5 && (_errType & jec::kRelativeStatEC)) ||
       (x>=2.5 && x<3.0 && (_errType & jec::kRelativeStatEC)) ||
       (x>=3.0 && x<5.5 && (_errType & jec::kRelativeStatHF )) ) {
    
    _jecL2ResStat->setJetEta(eta);
    _jecL2ResStat->setJetPt(pTprime); // pT doesn't matter
    double stat = _jecL2ResStat->getCorrection();

    err = stat;
  }
  
  if ( (_errType & jec::kRelativeStatFSR ) ) {

    // On 13 Jan 2015, at 14:46, from Denis Rathjens
    // Re: Out of Office AutoReply: Pseudo Pt plot
    if (!fkFSR) fkFSR = new TF1("fkFSR",_kFSR,0,5.2,10);
    fkFSR->SetParameters(1.00187, -0.00207088, 0.276750,
			 1.46992e-06,
			 -2.15806e-06, 3.3642e-06,
			 0.000333186, -0.000551704, 0.0984257,
			 0);
    fkFSR->SetParameter(9, +1);
    double stat_up = fkFSR->Eval(eta);
    fkFSR->SetParameter(9, -1);
    double stat_dw = fkFSR->Eval(eta);
    double stat = 0.5 * ( stat_up - stat_dw );
    err = sqrt(err*err + stat*stat);
  }
  
  return err;
} // RelativeStat

// Uncertainty in L2Res pT dependence: log-linear fit vs constant fit
// To-do: solve pTprime with Brent's method to compare at same pTprime?
double JECUncertainty::_RelativePt(const double pTprime,
				   const double eta) {

  // limit pt to accessible range
  const double ptmin = 10.;
  const double emax = 6500;
  double pt = max(ptmin, min(pTprime, emax/cosh(eta)));

  assert(_jecL2ResFlat);
  double rflat = _Rjet(pTprime, eta, -1, -1, _jecL2ResFlat);

  assert(_jecL2ResPt);
  double rpt = _Rjet(pTprime, eta, -1, -1, _jecL2ResPt);

  // Use 50% of the slope as uncertainty consistently everywhere
  // We do correct for it for a reason, so 100% seems too conservative
  double kfactor = 0.5;
  double err = kfactor * (rflat / rpt - 1);

  double x = fabs(eta);
  if ((x>=0.0 && x<1.3 && _errType & jec::kRelativePtBB) ||
      (x>=1.3 && x<2.5 && _errType & jec::kRelativePtEC1) ||
      (x>=2.5 && x<3.0 && _errType & jec::kRelativePtEC2) ||
      (x>=3.0 && x<5.5 && _errType & jec::kRelativePtHF))
    return err;

  return 0;
} // RelativePt

double JECUncertainty::_RelativeBal(const double pTprime,
				    const double eta) {

  // limit pt to accessible range
  const double ptmin = 10.;
  const double emax = 6500;
  double pt = max(ptmin, min(pTprime, emax/cosh(eta)));

  // Take full MPF vs pTbal difference as extra source
  // These two should agree, but don't, yet
  assert(_jecL2ResMPF);
  double rmpf = _Rjet(pTprime, eta, -1, -1, _jecL2ResMPF);
  assert(_jecL2ResBal);
  double rbal = _Rjet(pTprime, eta, -1, -1, _jecL2ResBal);
  
  double kfactor = 1;
  double err = kfactor*(rbal / rmpf - 1);

  return err;
} // RelativeBal

double JECUncertainty::_RelativePrefire(const double pTprime,
				    const double eta) {

  // limit pt to accessible range
  const double ptmin = 10.;
  const double emax = 6500;
  double pt = max(ptmin, min(pTprime, emax/cosh(eta)));

  // Take full with/without L1BXcleaning MPF loglin difference as extra source
  // These two should agree, but don't, yet
  assert(_jecL2ResPrefireFilterYes);
  double rYes = _Rjet(pTprime, eta, -1, -1, _jecL2ResPrefireFilterYes);
  assert(_jecL2ResPrefireFilterNo);
  double rNo = _Rjet(pTprime, eta, -1, -1, _jecL2ResPrefireFilterNo);

  //  std::cout << "rYes: " << rYes << " rNo: " << rNo << std::endl;
  double kfactor = 1;
  double err = kfactor*(rYes / rNo - 1);

  return err;
} // RelativePrefire

double JECUncertainty::_RelativeSample(const double pTprime,
				       const double eta) {

  // limit pt to accessible range
  const double ptmin = 10.;
  const double emax = 6500;
  double pt = max(ptmin, min(pTprime, emax/cosh(eta)));

  // Take Z/gam+jet-only closure fit as extra source
  // Nominal correction from dijet, so dijet should close
  assert(_jecL2L3Res2ParDiJet);
  double rdj = _Rjet(pt, eta, -1, -1, _jecL2L3Res2ParDiJet);
  //  assert(_jecL2L3ResCustomHybridZGamJet);
  //  double rzg = _Rjet(pt, eta, -1, -1, _jecL2L3ResCustomHybridZGamJet);
  
  assert(_jecL2L3Res2ParZGamJet);
  double rzg = _Rjet(pt, eta, -1, -1, _jecL2L3Res2ParZGamJet);

  double kfactor = 1;
  double loglinerror = kfactor*(rzg / rdj - 1);
    
  //outside of eta 2.65, take const fit difference for now
  if(fabs(eta)>2.649){
    assert(_jecL2L3ResConstDiJet);
    rdj = _Rjet(pt, eta, -1, -1, _jecL2L3ResConstDiJet);
    assert(_jecL2L3ResConstZGamJet);
    rzg = _Rjet(pt, eta, -1, -1, _jecL2L3ResConstZGamJet);
  }
  //  double kfactor = 1;
  double err = kfactor*(rzg / rdj - 1);
  //double err = kfactor*(rzg - 1);
  
  //  return fabs(err)>fabs(loglinerror) ? err : loglinerror;
  return err; //return const difference above 2.65
} // RelativeSample

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
  
  if (debug) cout << "_PileUpDataMC" << endl << flush;

  // Guesstimate of data/MC scale factor uncertainty
  //double errsf = 0.05; // 03FebV3
  double errsf = 0.10; // 17NovV6

  double l1 = _Rjet(pTprime, eta, -1, -1, _jecL1MCpt);
  double sys = fabs(l1-1) * errsf;
  		     
  return sys;
} // PileUpDataMC


// PileUpPt re-implemented for Run II
// Simplified a bit by using up variation in MC only
// no separate reference algo for L2res and L3res (AK4PFchs both),
// and with new L2res, L3res parameterizations
double JECUncertainty::_PileUpPt(const double pTprime, const double eta) {

  if (debug) cout << "_PileUpPt" << endl << flush;

  double x = fabs(eta);
  double etax = max(-4.5,min(4.5,eta));

  FactorizedJetCorrector *_l1flat = _jecL1MCflat; assert(_l1flat);
  FactorizedJetCorrector *_l1pt = _jecL1MCpt;     assert(_l1pt);

  FactorizedJetCorrector *_l1flatref = _jecL1MCflat_ak4pfchs; 
  FactorizedJetCorrector *_l1ptref = _jecL1MCpt_ak4pfchs;
  assert(_l1flatref);
  assert(_l1ptref);

  double sysref(0), syseta(0), syszero(0), sys(0);

  if (!_fl3 && !_hl3 && !_hl3ref) {

    // Shape used in L3Res fit (17NovV6: 2-par fit)
    _fl3 = new TF1("fl3","[0]+[1]*max(0.,[2]+[3]*pow(x,[4]))",30,1400);
    _fl3->SetParameters(1,0, 0.03091e+00,-5.11540e-02,-1.54227e-01); // SPRH
    _fl3->FixParameter(2, 0.03091e+00); // shape fixed in global fit
    _fl3->FixParameter(3, -5.11540e-02); // shape fixed in global fit
    _fl3->FixParameter(4, -1.54227e-01); // shape fixed in global fit

    const double x_pt[] =
      {10, 12, 15, 18, 21, 24, 
       28, 32, 37, 43, 49, 56, 64, 74, 84,
       97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 362, 430,
       507, 592, 686, 790, 905, 1032};
    const int ndiv_pt = sizeof(x_pt)/sizeof(x_pt[0])-1;

    _hl3 = new TH1D("hl3",";p_{T} (GeV);Offset error",ndiv_pt,x_pt);
    _hl3ref = new TH1D("hl3ref",";p_{T} (GeV);Offset error",ndiv_pt,x_pt);

    double x_eta[] = //{0,0};
      {-1.4,-1.2,-1.0, -0.8,-0.6,-0.4,-0.2,0.,
       0, 0.2,0.4,0.6,0.8,1.0, 1.2,1.4};
    const int ndiv_eta = sizeof(x_eta)/sizeof(x_eta[0])-1;
    
    // Average offset uncertainty over barrel (|eta|<1.3)
    // Ways to estimate (not yet done):
    // 1) fit (L1Data-L1RC) difference in L3Res global fit (-100 +/- 40% or so)
    // 2) use ZeroBias jet pT>30 GeV stability vs mu (100%)
    // 3) use Z-jet pT (MPF) difference vs mu at low pT (??)
    const double kfactor(1);
    for (int ipt = 0; ipt != ndiv_pt; ++ipt) {
      
      double pt = 0.5*(x_pt[ipt] + x_pt[ipt+1]);
      
      double sumw(0), sumsys(0), sumsysref(0);
      for (int ieta = 0; ieta != ndiv_eta; ++ieta) {
	
	double etab = 0.5*(x_eta[ieta]+x_eta[ieta+1]);
	double l1f = _Rjet(pt, etab, -1, -1, _l1flat);
	double l1p = _Rjet(pt, etab, -1, -1, _l1pt);
	double sys = kfactor * (l1f / l1p - 1);
	sumsys += sys;
	sumw   += 1;

	// L3res is really only fitted to AK4PFchs (=ref)
	// Fix 2015-10-02: use AK4 area (0.503)
	double l1fref = _Rjet(pt, etab, 0.503, -1, _l1flatref);
	double l1pref = _Rjet(pt, etab, 0.503, -1, _l1ptref);
	double sysref = kfactor * (l1fref / l1pref - 1);
	sumsysref += sysref;
      } // for ieta
      double sys = sumsys / sumw;
      double sysref = sumsysref / sumw;

      _hl3->SetBinContent(_hl3->FindBin(pt), sys);
      _hl3ref->SetBinContent(_hl3ref->FindBin(pt), sysref);
    } // for ipt
    
    _hl3ref->Fit(_fl3, "QRN");
  } // !_fl3

  // Residual offset uncertainty remaining after L3Res
  // Double-counting a bit, because also included in L2(L3)Res part
  // Maybe relevant, if reference algos for dijet and Z+jet balance differ
  // (e.g. AK8chs vs AK4chs only)
  if (_errType & jec::kPileUpPtRef) {
    // Don't recalculate eta average, but just interpolate in pT
    double pt = max(10.,min(1032.,pTprime));
    sysref = _hl3->Interpolate(pt) - _fl3->Eval(pt);
  } // sysref

  // Residual relative PU offset after applying dijet balance comes from
  // fitting offset residual with L2res shape, for AK4PFchs (=ref)
  // (now a constant, was [0]+[1]*log(x) in Run I)
  // The residual offset from barrel cancels out some of the residual
  // offset in the forward regions in the pTprobe/pTtag ratio
  // (use +100% consistently for forward variation and barrel average)
  if ( (x>=0.0 && x<1.3 && _errType & jec::kPileUpPtBB) ||
       (x>=1.3 && x<2.5 && _errType & jec::kPileUpPtEC1) ||
       (x>=2.5 && x<3.0 && _errType & jec::kPileUpPtEC2) ||
       (x>=3.0 && x<5.5 && _errType & jec::kPileUpPtHF) ||
       _errType & jec::kPileUpMuZero ) {

    double kfactor = 1;

    const double x_pt[] =
      {28, 32, 37, 43, 49, 56, 64, 74, 84,
       97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 362, 430,
       507, 592, 686, 790, 905, 1032};
    const int ndiv_pt = sizeof(x_pt)/sizeof(x_pt[0])-1;

    if (!_hl2ref) _hl2ref = new TH1D("hl2ref",";p_{T} (GeV); Offset residual",
				     ndiv_pt, x_pt);

    for (int ipt = 0; ipt != ndiv_pt; ++ipt) {

      double pt = 0.5*(x_pt[ipt] + x_pt[ipt+1]);
      // Fix 2015-10-02: use AK4 area (0.503)
      double l1fref = _Rjet(pt, etax, 0.503, -1, _l1flatref);
      double l1pref = _Rjet(pt, etax, 0.503, -1, _l1ptref);
      // probe vs tag, so barrel offset residual cancels out
      double sysrefb = kfactor*(l1fref/l1pref-1) - _hl3ref->Interpolate(pt);
      _hl2ref->SetBinContent(_hl2ref->FindBin(pt), sysrefb);
    } // for ipt
    
    TF1 *fl2 = new TF1("fl2", "[0]+[1]*log(x)", 62, 3250/cosh(etax));
    _hl2ref->Fit(fl2,"QRN");

    // Residual offset uncertainty remaining after L2Res
    if ( (x>=0.0 && x<1.3 && _errType & jec::kPileUpPtBB) ||
	 (x>=1.3 && x<2.5 && _errType & jec::kPileUpPtEC1) ||
	 (x>=2.5 && x<3.0 && _errType & jec::kPileUpPtEC2) ||
	 (x>=3.0 && x<5.5 && _errType & jec::kPileUpPtHF) ) {

      double l1f = _Rjet(pTprime, etax, -1, -1, _l1flat);
      double l1p = _Rjet(pTprime, etax, -1, -1, _l1pt);
      //double sysb = kfactor * (l1f / l1p - 1) - _hl3->Interpolate(pTprime);
      // probe corrected for offset absorbed in L3Res and L2Res for AK4PFchs
      double sysb = kfactor * (l1f / l1p - 1) - _fl3->Eval(pTprime);
      syseta = sysb - fl2->Eval(pTprime); // Run II
    }

    // Residual offset absorbed into L3Res and biasing <mu>=0
    if (_errType & jec::kPileUpMuZero) {
      syszero = fl2->Eval(pTprime) + _fl3->Eval(pTprime); // Run II
    }
    
    delete fl2;
  } // syseta

  sys = sqrt(sysref*sysref + syseta*syseta + syszero*syszero);

  // For single sources
  if (!(_errType & ~jec::kPileUpPtRef)) return sysref;
  if (!(_errType & ~jec::kPileUpPt)) return syseta;
  if (!(_errType & ~jec::kPileUpPtBB)) return syseta;
  if (!(_errType & ~jec::kPileUpPtEC1)) return syseta;
  if (!(_errType & ~jec::kPileUpPtEC2)) return syseta;
  if (!(_errType & ~jec::kPileUpPtHF)) return syseta;
  if (!(_errType & ~jec::kPileUpPtEta)) return syseta;
  if (!(_errType & ~jec::kPileUpMuZero)) return syszero;

  return sys;
} // _PileUpPt

// (Obsolete version of) Pile-up uncertainty from pT dependence
// Implemented as difference between MC truth (V1) and Random Cone (V0),
double JECUncertainty::_PileUpEnvelope(const double pTprime, const double eta) {

  if (debug) cout << "_PileUpEnvelope" << endl << flush;

  // Limit eta to [-5,5] because V0 files don't go further out
  // Even closer in, because bias goes nuts in the last bin out
  double maxeta = 4.5;
  double etax = max(-maxeta,min(maxeta,eta));

  double kfactor = 1;

  FactorizedJetCorrector *_l1flat = _jecL1MCflat;
  FactorizedJetCorrector *_l1pt = _jecL1MCpt;
  //
  double l1f = _Rjet(pTprime, etax, -1, -1, _l1flat);
  double l1p = _Rjet(pTprime, etax, -1, -1, _l1pt);
  double sys = kfactor * (l1f / l1p - 1);

  return sys;
} // _PileUpEnvelope

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
double JECUncertainty::_Time(const double pt, const double eta) {

  // Optional epoch time uncertainties
  double spt(0);
  if (_errType & jec::kTimeRunB) {
    assert( !(_errType & (jec::kTimePtEtaMask & ~jec::kTimeRunB)) ); 
    spt = _TimePtEta(pt,eta,1);
  }
  if (_errType & jec::kTimeRunC) {
    assert( !(_errType & (jec::kTimePtEtaMask & ~jec::kTimeRunC)) ); 
    spt = _TimePtEta(pt,eta,2);
  }
  if (_errType & jec::kTimeRunDE) {
    assert( !(_errType & (jec::kTimePtEtaMask & ~jec::kTimeRunDE)) ); 
    spt = _TimePtEta(pt,eta,3);
  }
//  if (_errType & jec::kTimeRunE) {
//    assert( !(_errType & (jec::kTimePtEtaMask & ~jec::kTimeRunE)) ); 
//    spt = _TimePtEta(pt,eta,4);
//  }
  if (_errType & jec::kTimeRunF) {
    assert( !(_errType & (jec::kTimePtEtaMask & ~jec::kTimeRunF)) ); 
    spt = _TimePtEta(pt,eta,4);
  }
  // Normal time uncertainties
  if (_errType & jec::kTimePtEta) {
    assert( !(_errType & (jec::kTimePtEtaMask & ~jec::kTimePtEta)) ); 
    spt  =   _TimePtEta(pt,eta,0);
  }

  double err = spt;

  // signed sources
  if (!(_errType & ~jec::kTimePtEta))  return spt;
  if (!(_errType & ~jec::kTimeRunB))  return spt;
  if (!(_errType & ~jec::kTimeRunC))  return spt;
  if (!(_errType & ~jec::kTimeRunDE))  return spt;
  //  if (!(_errType & ~jec::kTimeRunE))  return spt;
  if (!(_errType & ~jec::kTimeRunF))  return spt;

  return err;
}

// Time-dependence uncertainty from L2L3Res
// Time uncertainty is difference between directly fitted BCDEF
// and weighted average of the four IOVs
double JECUncertainty::_TimePtEta(const double pt, const double eta,
				  int epoch) {

  // epochs: 0 (runsBCDEF), 1 (runB), 2 (runC), 3 (runDE), 4 (runF) //for V31 ff.
  assert(epoch>=0 && epoch<=4);

  const int nepoch = 4;
  double lum[nepoch] = {4.8, 9.6, 13.5, 13.4}; // B,C,DE,f; tot 41.4(3)/fb

  double sumlum(0);
  for (int i = 0; i != nepoch; ++i) sumlum += lum[i];
  if (debug) cout << "Total lumi: " <<  sumlum << " fb-1" <<  endl;

  double w[nepoch];
  for (int i = 0; i != nepoch; ++i) w[i] = lum[i] / sumlum;

  assert(_jecB);
  assert(_jecC);
  assert(_jecDE);
  assert(_jecF);
  assert(_jecBCDEF);
  double jecB = _Rjet(pt, eta, -1, -1, _jecB);
  double jecC  = _Rjet(pt, eta, -1, -1, _jecC);
  double jecDE   = _Rjet(pt, eta, -1, -1, _jecDE);
  double jecF   = _Rjet(pt, eta, -1, -1, _jecF);
  double jecs[nepoch] = {jecB, jecC, jecDE, jecF};

  // BCDEF from weighted average of IOVs
  double jecSum(0);
  for (int i = 0; i != nepoch; ++i) {
    jecSum += w[i] * jecs[i];
  }

  // Directly fitted BCDEF
  double jecAll = _Rjet(pt, eta, -1, -1, _jecBCDEF);

  // Indirect epoch from direct BCDEFGH minus other epochs, All~Sum
  // Sum = w1*jes1 + w2*jes2 + ... => jes1 ~ (All - w2*jes2 - ...)/w1
  //double sumOth(0);
  //for (int i = 0; i != nepoch; ++i) {
  //if (i != epoch-1) sumOth += w[i] * jecs[i];
  //}
  //double jecAlt = (jecAll - sumOth) / w[epoch-1];
  //
  // double err(0);
  //if (epoch==0) err = jecAll - jecSum;
  //if (epoch!=0) err = jecAlt - jecs[epoch-1];
  // => didn't work too well, all periods on same side
  //
  // Try instead All = Sum + Delta, such that
  // Delta = k * sqrt(sum_j w_i * (jes_i-Sum)^2),
  // solve for k and use delta_i = k * (jes_i-Sum)
  double sumDelta2(0);
  for (int i = 0; i != nepoch; ++i) {
    if (i != epoch-1) sumDelta2 += w[i] * pow(jecs[i] - jecSum, 2);
  }
  double Delta = jecAll - jecSum;
  double kErr = Delta / sqrt(sumDelta2);

  double err(0);
  if (epoch==0) err = jecAll - jecSum;
  if (epoch!=0) err = kErr * (jecs[epoch-1] - jecSum);

  return err;
}

// Move to mu-based mapping, which is better for comparing
// different PU scenarios as it considers both IT and OOT PU,
// plus we have a number directly comparable to ATLAS
double JECUncertainty::_RhoFromMu(double mu) {
  // Eta_0.0-1.3, jt320
  //return (1.01272 + 0.551183*mu + 0.000362936*mu*mu); Run I?
  // RunG data, Eta_0.0-1.3, jt320
  return (1.7964 + 0.565301*mu + -2.22401e-04*mu*mu); // 80XV8
}
double JECUncertainty::_NpvFromMu(double mu) {
  // Eta_0.0-1.3, jt400
  //return (0.851334 + 0.722608*mu - 0.00184534*mu*mu);
  // RunG data, Eta_0.0-1.3, jt320
  return (1.92228 + 0.666238*mu + -5.35377e-04*mu*mu); // 80XV8
}
