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
enum jer_iov { none, run1, run2016, run2017, run2018};
jer_iov _jer_iov = none;

const int _nres = 8;

// Define as global variables to make code faster
double _rho = 19.31; // Data-Fall17V8-D, jt500 hrho
JME::JetResolution *_jer(0);
JME::JetResolutionScaleFactor *_jer_sf(0);

// From Sanmay by email 30 Mar 2015
// (direct recommendations from different eta bins => not ideal?)
/*
const double kpar2012[_nres][2] = {
  {1.079, 0.026},
  {1.099, 0.028},
  {1.121, 0.029},
  {1.208, 0.039},
  {1.208, 0.053},
  {1.395, 0.062},
  {1.395, 0.062}, // copy
  {1.395, 0.062}, // copy
  };*/
const double kpar2012[_nres][2] = {
// JER SF produced with minitools/resolution.C:redoJER("Run1")
 {1.079, 0.026}, // 0.0-0.5
 {1.098, 0.028}, // 0.5-1.0
 {1.114, 0.029}, // 1.0-1.5
 {1.163, 0.037}, // 1.5-2.0
 {1.222, 0.051}, // 2.0-2.5
 {1.283, 0.062}, // 2.5-3.0
 {1.324, 0.088}, // 3.0-3.2
 {1.056, 0.191}, // 3.2-4.7
};


// ptresolution.h has Data/MC ratio produced with ak7ak5ratio.C:redoJER()

// Placeholder
const double kpar2016[_nres][2] = {
// JER SF produced with minitools/resolution.C:redoJER("Run2016")
 {1.159, 0.064}, // 0.0-0.5
 {1.173, 0.064}, // 0.5-1.0
 {1.145, 0.089}, // 1.0-1.5
 {1.119, 0.104}, // 1.5-2.0
 {1.197, 0.157}, // 2.0-2.5
 {1.421, 0.204}, // 2.5-3.0
 {1.188, 0.128}, // 3.0-3.2
 {1.192, 0.149}, // 3.2-4.7
};

// Placeholder
const double kpar2017[_nres][2] = {
// JER SF produced with minitools/resolution.C:redoJER("Run2017")
 {1.143, 0.022}, // 0.0-0.5
 {1.145, 0.046}, // 0.5-1.0
 {1.115, 0.114}, // 1.0-1.5
 {1.153, 0.134}, // 1.5-2.0
 {1.305, 0.173}, // 2.0-2.5
 {2.011, 0.511}, // 2.5-3.0
 {1.254, 0.114}, // 3.0-3.2
 {1.154, 0.152}, // 3.2-4.7
};

// Placeholder
const double kpar2018[_nres][2] = {
// JER SF produced with minitools/resolution.C:redoJER("Run2018")
 {1.150, 0.043}, // 0.0-0.5
 {1.121, 0.066}, // 0.5-1.0
 {1.114, 0.127}, // 1.0-1.5
 {1.125, 0.192}, // 1.5-2.0
 {1.196, 0.222}, // 2.0-2.5
 {2.020, 0.579}, // 2.5-3.0
 {1.206, 0.194}, // 3.0-3.2
 {1.082, 0.198}, // 3.2-4.7
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
   {2.86, 0.874, 0.0000}, // copy
   {2.86, 0.874, 0.0000}}; // copy

//////////////////////////////////////////////////////////
// NB: resolution.C uses JER to set Gaussian fit range so
//     iterating on the fit may change the fitted values.
//     Here we have started from 2012 JER and iterated once,
//     only updating fits that had improved chi2.
/////////////////////////////////////////////////////////
const double vpar2016[_nres][3] = 
// Fit of JER for R=0.4, 13 TeV file MC-Legacy16-GH
  {{-2.37, 0.948, 0.0350},  // y 0.0-0.5, chi2 68.7/35
   {-2.62, 0.982, 0.0365},  // y 0.5-1.0, chi2 109.9/35
   {-2.39, 1.086, 0.0485},  // y 1.0-1.5, chi2 111.8/35
   {-1.20, 1.070, 0.0311},  // y 1.5-2.0, chi2 92.5/35
   {1.58, 1.003, 0.0248},  // y 2.0-2.5, chi2 31.8/35
   {3.33, 0.964, 0.0242},  // y 2.5-3.0, chi2 53.6/33
   {3.86, 0.829, 0.1233},  // y 3.0-3.2, chi2 33.8/33
   {2.16, 0.779, 0.0892}}; // y 3.2-4.7, chi2 50.5/33

const double vpar2017[_nres][3] = 
// Fit of JER for R=0.4, 13 TeV file P8CP5-17nov17-DE
  {{-2.73, 0.997, 0.0344},  // y 0.0-0.5, chi2 25.0/35
   {-2.73, 1.007, 0.0396},  // y 0.5-1.0, chi2 23.3/35
   {-2.98, 1.159, 0.0527},  // y 1.0-1.5, chi2 24.5/35
   {-1.16, 1.134, 0.0283},  // y 1.5-2.0, chi2 18.6/35
   {2.40, 1.082, 0.0210},  // y 2.0-2.5, chi2 26.9/35
   {5.24, 1.037, -0.0000},  // y 2.5-3.0, chi2 43.8/33
   {-2.57, 1.379, 0.1106},  // y 3.0-3.2, chi2 166.3/33
   {2.83, 0.772, 0.0967}}; // y 3.2-4.7, chi2 8.8/33

const double vpar2018[_nres][3] = 
// Fit of JER for R=0.4, 13 TeV file MC-Fall18V8-D
  {{-2.67, 0.991, 0.0330},  // y 0.0-0.5, chi2 33.5/35
   {-2.76, 1.009, 0.0374},  // y 0.5-1.0, chi2 24.0/35
   {-3.04, 1.156, 0.0516},  // y 1.0-1.5, chi2 41.9/35
   {1.56, 1.027, 0.0445},  // y 1.5-2.0, chi2 70.7/35
   {1.92, 1.052, 0.0422},  // y 2.0-2.5, chi2 39.6/35
   {3.49, 1.209, 0.0073},  // y 2.5-3.0, chi2 29.8/33
   {1.84, 0.965, 0.1333},  // y 3.0-3.2, chi2 615.9/33 (pretty bad)
   {2.24, 0.750, 0.0901}}; // y 3.2-4.7, chi2 76.2/33

double ptresolution(double pt, double eta) {

  //int ieta = min(_nres-1, int(fabs(eta) / 0.5));
  //int iy = min(5, ieta);
  //int iy = min(_nres-1, int(fabs(eta) / 0.5 + 0.5));
  int iy = min(_nres-1, int(fabs(eta) / 0.5));
  if (fabs(eta)>3.2) iy = _nres-1; // 3.2-4.7 bin instead of 3.5-4.7
  double res = 0;

  // Own parameterized JER
  if (!_usejme) {
    if (_jer_iov==none) {
      assert(false); // not initialized
    }
    if (_jer_iov==run1) {
      double kn = pow(0.4/0.5,2); //  scale noise term down by jet area
      res = sqrt(pow(kn*vpar2012[iy][0]/pt,2) + pow(vpar2012[iy][1],2)/pt + 
		 pow(vpar2012[iy][2],2));
      if (!_ismcjer) res *= kpar2012[iy][0];
    }
    if (_jer_iov==run2016) {
      res = sqrt(vpar2016[iy][0]*fabs(vpar2016[iy][0])/(pt*pt) +
		 pow(vpar2016[iy][1],2)/pt + 
		 pow(vpar2016[iy][2],2));
      if (!_ismcjer) res *= kpar2016[iy][0];
    }
    if (_jer_iov==run2017) {
      res = sqrt(vpar2017[iy][0]*fabs(vpar2017[iy][0])/(pt*pt) +
		 pow(vpar2017[iy][1],2)/pt + 
		 pow(vpar2017[iy][2],2));
      if (!_ismcjer) res *= kpar2017[iy][0];
    }
    if (_jer_iov==run2018) {
      res = sqrt(vpar2018[iy][0]*fabs(vpar2018[iy][0])/(pt*pt) +
		 pow(vpar2018[iy][1],2)/pt + 
		 pow(vpar2018[iy][2],2));
      if (!_ismcjer) res *= kpar2018[iy][0];
    }
  }

  // Official JME resolutions
  if (_usejme) {
    // Example code from:
    // https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_25/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h

    if (_jer_iov==none) {
      assert(false); // not initialized
    }
    if (_jer_iov==run1 && !_jer && !_jer_sf) {
      assert(false); // Run 1 not supported
    }
    if (_jer_iov==run2016 && !_jer && !_jer_sf) {
      string resolutionFile = "../JRDatabase/textFiles/Summer16_25nsV1_MC/"
	"Summer16_25nsV1_MC_PtResolution_AK4PFchs.txt";
      string scaleFactorFile = "../JRDatabase/textFiles/Summer16_25nsV1_MC/"
	"Summer16_25nsV1_MC_SF_AK4PFchs.txt";
      string weightFile = "rootfiles/jerweights.root";

      _jer = new JME::JetResolution(resolutionFile);
      _jer_sf = new JME::JetResolutionScaleFactor(scaleFactorFile);
    }
    if (_jer_iov==run2017 && !_jer && !_jer_sf) {
      string resolutionFile = "../JRDatabase/textFiles/Fall17_V3_MC/"
	"Fall17_V3_MC_PtResolution_AK4PFchs.txt";
      string scaleFactorFile = "../JRDatabase/textFiles/Fall17_V3_MC/"
	"Fall17_V3_MC_SF_AK4PFchs.txt";
      string weightFile = "rootfiles/jerweights.root";

      _jer = new JME::JetResolution(resolutionFile);
      _jer_sf = new JME::JetResolutionScaleFactor(scaleFactorFile);
    }
    if (_jer_iov==run2018 && !_jer && !_jer_sf) {
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

    //double rho = 19.31; // Data-Fall17V8-D, jt500 hrho
    double jet_resolution = _jer->getResolution({{JME::Binning::JetPt, pt}, {JME::Binning::JetEta, eta}, {JME::Binning::Rho, _rho}});
    double jer_sf = _jer_sf->getScaleFactor({{JME::Binning::JetEta, eta}}, Variation::NOMINAL);//m_systematic_variation);

    res = jet_resolution;
    if (!_ismcjer) res *= jer_sf;
  }

  return res;
}

#endif // __ptresolution_h__
