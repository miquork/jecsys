#ifndef __FLAVOR__
#define __FLAVOR__

// Purpose: class to facilitate retrieving flavor responses and fractions
#include "TF1.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include <string>
#include <map>

using namespace std;

class Flavor {

 public:
  Flavor(){};
  ~Flavor(){};

  // Reproduce style of Robin Aggleton's plots to check code
  bool doRobin = false;

  // Flavor response for arbitrary mix (unitarized) over eta
  double getResp(double pt, double absetamin, double absetamax,
		 double f[5], double w=0, string s="");
  // Flavor response for pre-defined mixture averaged over eta
  double getResp(double pt, double absetamin, double absetamax,
		 string smix, double w=0, string s="");
  // Flavor response for pre-defined mixture or pure flavor
  double getResp(double pt, double eta,
		 string smix, double w=0, string s="");
  // Flavor response for arbitrary mix (unitarized); w=-1:data, w=0:P8, w=1:HW
  double getResp(double pt, double eta, double f[5],
		 //double fud, double fs, double fc, double fb,
		 double w=0, string s="");

  // Flavor response uncertainty (Herwig-Pythia) for pre-defined mixture
  //double getUnc(double pt, double eta, string smix);
  // Flavor response uncertainty (Herwig-Pythia) for arbitrary mix (unitarized)
  //double getUnc(double pt, double eta,
  //		double fud, double fs, double fc, double fb);

  // Flavor fractions for pre-defined mixture
  //void getFracs(string smix, double &fud, double &fs, double &fc, double &fb);
  void getFracs(double pt, double eta, string smix, double (&f)[5]);

  // Load JEC objects into _mjec map
  std::map<std::string, FactorizedJetCorrector*> _mjec;
  std::map<std::string, FactorizedJetCorrector*> _mhw7;
  void loadJEC(std::string smc);

  // Load residual responses as TF1 into _mjec map
  std::map<std::string, TF1*> _mres;
  void loadRES();
};

#endif //__FLAVOR__
