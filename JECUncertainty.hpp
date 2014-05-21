// Adapted from D0 Experiment jetcorr/jetcorr/RjetCorr.hpp

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

#include <iostream>
#include <string>

namespace jec {
  
  enum JetAlgo  {AK5PF, AK5PFchs, AK5JPT, AK5CALO,
		 AK7PF, AK7PFchs, AK7CALO};
  enum DataType {DATA, MC, PY, HW};
}

Double_t _jesfit(Double_t *x, Double_t *p);

class JECUncertainty {
public:
  
  JECUncertainty(const jec::JetAlgo& algo = jec::AK5PF,
			const jec::DataType& type = jec::DATA,
			const jec::ErrorTypes& errType = jec::kData,
			const double npv = 20);
  ~JECUncertainty(){};
  
  double Uncert(const double pTprime, const double eta);
  
  void SetJetAlgo(const jec::JetAlgo& algo);
  //void SetDataType(const jec::DataType& type = jec::DATA);
  void SetNPV(const double npv);
  void SetErrType(const jec::ErrorTypes& errType);

  jec::ErrorTypes GetErrType();

  // Half-internal calculations
  // Keep accessible for plotting scripts

  // Jet response
  void _InitL1();
  void _InitJEC();
  void _InitL2Res();
  void _InitL3Res();
  double _Rjet(const double pTprime, const double eta);

  // Statistical and systematic uncertainties
  //double _SystErr(const double pTprime, const double eta);// const;

  double _Absolute(const double pTprime) const;
  double _AbsoluteStat(const double pTprime) const;
  double _AbsoluteScale() const;
  double _AbsoluteMPFBias() const;
  double _AbsoluteFlavorMapping() const;
  double _AbsoluteFrag(const double pTprime) const;
  double _AbsoluteSPR(const double pTprime) const;
  //
  double _Relative(double pTprime, double eta) const;
  double _RelativeJER(double eta) const;
  double _RelativeFSR(double eta) const;
  double _RelativeStat(double eta) const;
  double _RelativePt(double pTprime, double eta) const;
  //
  double _PileUp(double pTprime, double eta);
  double _PileUpDataMC(double pTprime, double eta);
  double _PileUpBias(double pTprime, double eta);
  double _PileUpPt(double pTprime, double eta);
  //
  double _Flavor(double pTprime, double eta) const;
  double _FlavorMixed(double pTprime, double eta, std::string smix) const;
  double _FlavorMix(double pTprime, double eta, double fl, double fg, 
		    double fc, double fb) const;
  double _FlavorResponse(double pTprime, double eta, int iflavor) const;
  double _FlavorFraction(double pTprime, double eta,
			 int iflavor, int isample) const;
  //
  double _Time(double eta) const;
  
  // pieces of L1Offset
  double _L1MCFlat(double pTraw, double eta);
  double _L1DataFlat(double pTraw, double eta);
  double _L1Data(double pTraw, double eta);
  double _L1MC(double pTraw, double eta);
  double _L1SF(double pTraw, double eta, double rho);
  double _Rho(double npvmean);

  // helpers for AbsoluteStat
  TF1 *_fjes;
  TMatrixD *_emat;
  double _jesfitunc(double x, TF1 *f, TMatrixD *emat) const;

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
  FactorizedJetCorrector *_jecWithL1V0;
  FactorizedJetCorrector *_jecL1DTflat;
  FactorizedJetCorrector *_jecL1DTpt;
  FactorizedJetCorrector *_jecL1MCflat;
  FactorizedJetCorrector *_jecL1pt;
  FactorizedJetCorrector *_jecL1sf;
  FactorizedJetCorrector *_jecL2ResFlat;
  FactorizedJetCorrector *_jecL2ResPt;
  FactorizedJetCorrector *_jecL2jerup;
  FactorizedJetCorrector *_jecL2jerdw;
  FactorizedJetCorrector *_jecL2stat;
  // Time dependence histograms
  //TH1D *_hMay10, *_hV4, *_hAug5, *_hV6, *_hV1;

  double _npv;
  double _pTprime, _eta, _rjet;

  // scale factor for AK7 offset (jet area R=0.7/R=0.5)
  double _ajet;

  class ResponseFunc : public ROOT::Math::IBaseFunctionOneDim
  {
  public:
    ResponseFunc(double pTprime, FactorizedJetCorrector *jec, int npv, double eta, double rho, double jeta) :
      ROOT::Math::IBaseFunctionOneDim(),_pTprime(pTprime),_jec(jec),_npv(npv),_eta(eta),_rho(rho),_jeta(jeta) {}
    
    double DoEval(double pTraw) const {
      _jec->setJetPt(pTraw);
      _jec->setJetE(pTraw*cosh(_eta)); 
      _jec->setJetEta(_eta); 
      _jec->setNPV(_npv);
      _jec->setRho(_rho); _jec->setJetA(_jeta);
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
    int _npv;
    double _eta;
    double _rho;
    double _jeta;
  };
};

#endif /* __JECUNCERTAINTY_SUMMER12__ */
