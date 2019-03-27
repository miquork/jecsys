#ifndef __ptresolution_h__
#define __ptresolution_h__

#include "TMath.h"

#include <string>

#include "JetMETCorrections/Modules/interface/JetResolution.h"

using std::string;

// Switch MC truth or data resolutions
bool _ismcjer = true;
// Use official JME resolutions
bool _usejme = true;
//bool _jerkscale = true;//false;
// Define cone size (default Run I AK5)
//bool _ak7 = false;
//bool _ak4 = true;
bool _run2012 = false;
bool _run2018 = true;

const int _nres = 7;

// Define as global variables to make code faster
JME::JetResolution *_jer(0);
JME::JetResolutionScaleFactor *_jer_sf(0);

// From Sanmay by email 30 Mar 2015
// (direct recommendations from different eta bins => not ideal?)
const double kpar2012[_nres][2] = {
  {1.079, 0.026},
  {1.099, 0.028},
  {1.121, 0.029},
  {1.208, 0.039},
  {1.208, 0.053},
  {1.395, 0.062},
  {1.395, 0.062}, // copy
};

// Placeholder
const double kpar2018[_nres][2] = {
  {1.079, 0.026},
  {1.099, 0.028},
  {1.121, 0.029},
  {1.208, 0.039},
  {1.208, 0.053},
  {1.395, 0.062},
  {1.395, 0.062}, // copy
};


// Resolutions with Pythia Z2* + ak5ak7resolution12.C
// (produced on iMac desktop, with ROOT 5.30/00, iterating with AK5+2sigma)
// On Sep 15, 2014, using Winter14_V5 private pre-version (root tuples v12)
// Fit of JER for R=0.5, 8 TeV, 53X
const double vpar2012[_nres][3] =
  {{3.13, 0.897, 0.0337},  // y 0.0-0.5, chi2 21.6/33
   {3.58, 0.868, 0.0376},  // y 0.5-1.0, chi2 12.9/33
   {4.78, 0.885, 0.0438},  // y 1.0-1.5, chi2 26.5/33
   {5.36, 0.745, 0.0265},  // y 1.5-2.0, chi2 14.5/32
   {4.45, 0.680, 0.0253},  // y 2.0-2.5, chi2 10.5/25
   {2.86, 0.874, 0.0000}, // y 2.5-3.0, chi2 10.8/19
   {2.86, 0.874, 0.0000}}; // copy

const double vpar2018[_nres][3] = 
// Fit of JER for R=0.4, 13 TeV file Fall18V8-D
  {{-2.67, 0.991, 0.0330},  // y 0.0-0.5, chi2 33.5/35
   {-2.64, 1.009, 0.0374},  // y 0.5-1.0, chi2 34.4/35
   {-3.01, 1.155, 0.0514},  // y 1.0-1.5, chi2 31.3/35
   {1.08, 1.035, 0.0452},  // y 1.5-2.0, chi2 114.7/35
   {1.78, 1.054, 0.0419},  // y 2.0-2.5, chi2 72.8/35
   {1.01, 1.237, 0.0000}, // y 2.5-3.0, chi2 130.6/33
   {1.01, 1.237, 0.0000}}; // copy

double ptresolution(double pt, double eta) {

  //int ieta = min(_nres-1, int(fabs(eta) / 0.5));
  //int iy = min(5, ieta);
  //int iy = min(_nres-1, int(fabs(eta) / 0.5 + 0.5));
  int iy = min(_nres-1, int(fabs(eta) / 0.5));
  double res = 0;

  // Own parameterized JER
  if (!_usejme) {
    if (_run2012) {
      double kn = pow(0.4/0.5,2); //  scale noise term down by jet area
      res = sqrt(pow(kn*vpar2012[iy][0]/pt,2) + pow(vpar2012[iy][1],2)/pt + 
		 pow(vpar2012[iy][2],2));
      if (!_ismcjer) res *= kpar2012[iy][0];
    }
    if (_run2018) {
      res = sqrt(pow(vpar2018[iy][0]/pt,2) + pow(vpar2018[iy][1],2)/pt + 
		 pow(vpar2018[iy][2],2));
      if (!_ismcjer) res *= kpar2018[iy][0];
    }
  }

  // Official JME resolutions
  if (_usejme) {
    // Example code from:
    // https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_25/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h

    if (_run2018 && !_jer && !_jer_sf) {
      string resolutionFile = "../JRDatabase/textFiles/Autumn18_V1_MC/"
	"Autumn18_V1_MC_PtResolution_AK4PFchs.txt";
      string scaleFactorFile = "../JRDatabase/textFiles/Autumn18_V1_MC/"
	"Autumn18_V1_MC_SF_AK4PFchs.txt";
      string weightFile = "rootfiles/jerweights.root";

      _jer = new JME::JetResolution(resolutionFile);
      _jer_sf = new JME::JetResolutionScaleFactor(scaleFactorFile);
    }

    assert(_jer);
    assert(_jer_sf);
    

    //double sumres2(0);
    //const int neta = 1;
    //const double veta[neta+1] =
    //{eta,eta};
    //const int nrho = 1;
    //const double vrho[nrho+1] =
    //{rho,rho}

    double rho = 19.31; // Data-Fall17V8-D, jt500 hrho
    double jet_resolution = _jer->getResolution({{JME::Binning::JetPt, pt}, {JME::Binning::JetEta, eta}, {JME::Binning::Rho, rho}});
    double jer_sf = _jer_sf->getScaleFactor({{JME::Binning::JetEta, eta}}, Variation::NOMINAL);//m_systematic_variation);

    res = jet_resolution;
    if (!_ismcjer) res *= jer_sf;
  }

  return res;
}

#endif // __ptresolution_h__
