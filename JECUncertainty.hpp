#ifndef __JECUNCERTAINTY__
#define __JECUNCERTAINTY__

//
// Purpose: JEC uncertainty sources for correlation analysis
//
#define STANDALONE
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "ErrorTypes.hpp"
#include "JetDefs.hpp"

#include "Math/IFunction.h"
// ROOT (root.cern.ch) modules
#include "TMatrixD.h"
#include "TF1.h"
#include "TH1D.h"

#include <iostream>
#include <string>

namespace jec {
  
  enum JetAlgo  {AK4PF, AK4PFchs, AK4PFpuppi, AK4CALO, AK4JPT,
		 AK5PF, AK5PFchs, AK5CALO, AK5JPT,
		 AK7PF, AK7PFchs, AK7CALO, AK7JPT,
		 AK8PF, AK8PFchs, AK8PFpuppi, AK8CALO, AK8JPT};
  enum DataType {DATA, MC, PY, HW};
}

// Figure out a better way instead of global variables. Static?
TF1 *_fhb(0);
TF1 *_fl1(0);
Double_t _jesfit(Double_t *x, Double_t *p);
//Double_t _jeshb(double pt, double hb);

inline double absmax(double a, double b) {
  return (fabs(a) > fabs(b) ? a : b);
}

class JECUncertainty {
public:
  
  JECUncertainty(const jec::JetAlgo& algo = jec::AK4PF,
		 const jec::DataType& type = jec::DATA,
		 const jec::ErrorTypes& errType = jec::kData,
		 const double mu = 24.68);//19.81);
  ~JECUncertainty(){};
  
  double Uncert(const double pTprime, const double eta);
  // double Rjet(const double pTprime, const double eta); // add this?

  private:

  // Jet response
  void _InitL1();
  void _InitJEC();
  void _InitL2Res();
  void _InitL3Res();
  double _Rjet(const double pTprime, const double eta,
	       const double ajet, const double mu,
	       FactorizedJetCorrector *jec);

  // Statistical and systematic uncertainties
  double _Absolute(const double pTprime) const;
  double _AbsoluteStat(const double pTprime) const;
  double _AbsoluteScale() const;
  double _AbsoluteMPFBias() const;
  double _AbsoluteFlavorMapping() const;
  double _AbsoluteFrag(const double pTprime) const;
  double _AbsoluteSPR(const double pTprime) const;
  double _AbsoluteSPRH(const double pTprime) const;
  double _AbsoluteSPRE(const double pTprime) const;
  //
  double _Relative(double pTprime, double eta);
  double _RelativeJER(double pTprime, double eta);
  double _RelativeFSR(double pTprime, double eta);
  double _RelativeStat(double pTprime, double eta) const;
  double _RelativePt(double pTprime, double eta);
  double _RelativeBal(double pTprime, double eta);
  double _RelativeSample(double pTprime, double eta);
  //
  double _PileUp(double pTprime, double eta);
  double _PileUpDataMC(double pTprime, double eta);
  double _PileUpPt(double pTprime, double eta);
  double _PileUpEnvelope(double pTprime, double eta);
  //
  double _Flavor(double pTprime, double eta) const;
  double _FlavorMixed(double pTprime, double eta, std::string smix) const;
  double _FlavorMix(double pTprime, double eta, double fl, double fg, 
		    double fc, double fb) const;
  double _FlavorResponse(double pTprime, double eta, int iflavor) const;
  double _FlavorFraction(double pTprime, double eta,
			 int iflavor, int isample) const;
  //
  double _Time(double pTprime, double eta);
  //double _TimeEta(const double eta);
  double _TimePtEta(const double pt, const double eta, int epoch=0);
  
  // pieces of L1Offset
  double _L1MCFlat(double pTraw, double eta);
  double _L1DataFlat(double pTraw, double eta);
  double _L1Data(double pTraw, double eta);
  double _L1MC(double pTraw, double eta);
  double _L1SF(double pTraw, double eta, double rho);

  double _RhoFromMu(double mu);
  double _NpvFromMu(double mu);

  // helpers for AbsoluteStat
  TF1 *_fjes;
  TMatrixD *_emat;
  double _jesfitunc(double x, TF1 *f, TMatrixD *emat) const;

  // helpers for calculating PileUpPt systematics
  TF1 *_fl3, *_fl2;
  TH1D *_hl3, *_hl3ref, *_hl2ref;

private:
  jec::JetAlgo _algo;
  jec::DataType _type;
  bool _calo;
  bool _jpt;
  bool _pfchs;
  bool _pflow;
  bool _ideal;
  bool _trkbase;
  jec::ErrorTypes _errType;
  bool _isMC;

  // JEC variations
  FactorizedJetCorrector *_jec;
  FactorizedJetCorrector *_jecDefault;
  //
  FactorizedJetCorrector *_jecBCD;
  FactorizedJetCorrector *_jecEF;
  FactorizedJetCorrector *_jecG;
  FactorizedJetCorrector *_jecH;
  FactorizedJetCorrector *_jecBCDEFGH;
  //
  //FactorizedJetCorrector *_jecWithL1RC;
  //FactorizedJetCorrector *_jecL1DTflat;
  //FactorizedJetCorrector *_jecL1DTpt;
  FactorizedJetCorrector *_jecL1MCflat;
  FactorizedJetCorrector *_jecL1MCpt;
  //FactorizedJetCorrector *_jecL1sf;
  FactorizedJetCorrector *_jecL1MCflat_ak4pfchs;
  FactorizedJetCorrector *_jecL1MCpt_ak4pfchs;
  //
  FactorizedJetCorrector *_jecL2ResFlat;
  FactorizedJetCorrector *_jecL2ResPt;
  FactorizedJetCorrector *_jecL2ResMPF; // 80XV8
  FactorizedJetCorrector *_jecL2ResBal; // 80XV8
  FactorizedJetCorrector *_jecL2ResPY;
  FactorizedJetCorrector *_jecL2ResHW;
  FactorizedJetCorrector *_jecL2L3ResDiJet; // 03FebV7
  FactorizedJetCorrector *_jecL2L3ResZGamJet; // 03FebV7
  FactorizedJetCorrector *_jecL2L3ResZGamDiJet; // 03FebV7
  //JetCorrectionUncertainty *_uncL2ResStat;
  FactorizedJetCorrector *_jecL2ResStat;
  //
  FactorizedJetCorrector *_jecL2jerup;
  FactorizedJetCorrector *_jecL2jerdw;
  //FactorizedJetCorrector *_jecL2stat;

  JetCorrectionUncertainty *_uncRunI;

  double _mu;

  // scale factor for AK7 offset (jet area R=0.7/R=0.5)
  double _ajet;
};

  class ResponseFunc : public ROOT::Math::IBaseFunctionOneDim
  {
  public:
    ResponseFunc(double pTprime, FactorizedJetCorrector *jec, int npv, double eta, double rho, double jeta) :
      ROOT::Math::IBaseFunctionOneDim(),_pTprime(pTprime),_jec(jec),_npv(npv),_eta(eta),_rho(rho),_jeta(jeta) {}
    
    double DoEval(double pTraw) const {
      _jec->setJetPt(pTraw);
      _jec->setJetEta(_eta); 
      _jec->setRho(_rho);
      _jec->setJetA(_jeta);
      double cor = _jec->getCorrection();
      //std::cout << "in Brent: pTraw:" << pTraw << "  cor:" << cor << "  dist:" << pTraw * cor - _pTprime << '\n';
      return pTraw * cor - _pTprime; 
    }
    
    ROOT::Math::IBaseFunctionOneDim* Clone() const {
      return new ResponseFunc(_pTprime,_jec,_npv,_eta,_rho,_jeta);
    }
  private:
    double _pTprime;
    FactorizedJetCorrector *_jec;
    int _npv; // obsolete
    double _eta;
    double _rho;
    double _jeta;
  };


#endif /* __JECUNCERTAINTY__ */
