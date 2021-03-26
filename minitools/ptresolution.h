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
enum jer_iov { none, run1, run2016, run2017, run2018, run2018abc, run2018d,
	       run2016bcd, run2016ef, run2016gh,
	       run2017b, run2017c, run2017d, run2017e, run2017f, run2017de,
	       ul17, ul17b, ul17c, ul17d, ul17e, ul17f,
	       ul18, ul18a, ul18b, ul18c, ul18d,
               ul16apv, ul16bcd, ul16ef, ul16gh};
jer_iov _jer_iov = none;

const int _nres = 8;

// Define as global variables to make code faster
double _rho = 19.31; // Data-Fall17V8-D, jt500 hrho
JME::JetResolution *_jer(0);
JME::JetResolutionScaleFactor *_jer_sf(0);

bool _usenewsf = false;
TF1 *_newsf(0);

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
/*
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
*/
const double kpar2018[_nres][2] = {
// JER SF produced with minitools/resolution.C:redoJER("Run2018")
 {1.155, 0.031}, // 0.0-0.5
 {1.128, 0.045}, // 0.5-1.0
 {1.101, 0.066}, // 1.0-1.5
 {1.101, 0.070}, // 1.5-2.0
 {1.191, 0.083}, // 2.0-2.5
 {2.002, 0.479}, // 2.5-3.0
 {1.174, 0.039}, // 3.0-3.2
 {1.092, 0.049}, // 3.2-4.7
};
const double kpar2018abc[_nres][2] = {
// JER SF produced with minitools/resolution.C:redoJER("Run2018ABC")
 {1.168, 0.048}, // 0.0-0.5
 {1.130, 0.043}, // 0.5-1.0
 {1.093, 0.071}, // 1.0-1.5
 {1.096, 0.075}, // 1.5-2.0
 {1.170, 0.090}, // 2.0-2.5
 {1.663, 0.335}, // 2.5-3.0
 {1.206, 0.048}, // 3.0-3.2
 {1.073, 0.068}, // 3.2-4.7
};
const double kpar2018d[_nres][2] = {
// JER SF produced with minitools/resolution.C:redoJER("Run2018D")
 {1.159, 0.039}, // 0.0-0.5
 {1.140, 0.055}, // 0.5-1.0
 {1.118, 0.070}, // 1.0-1.5
 {1.100, 0.073}, // 1.5-2.0
 {1.229, 0.103}, // 2.0-2.5
 {2.295, 0.573}, // 2.5-3.0
 {1.163, 0.054}, // 3.0-3.2
 {1.078, 0.054}, // 3.2-4.7
};

// JER SF produced with minitools/resolution.C:redoJER("RunUL17")
const double kparul17[_nres+1][2] = {
 {1.108, 0.056}, // 0.0-0.5
 {1.112, 0.027}, // 0.5-1.0
 {1.149, 0.053}, // 1.0-1.5
 {1.145, 0.074}, // 1.5-2.0
 {1.132, 0.087}, // 2.0-2.5
 {1.257, 0.128}, // 2.5-3.0
 {1.238, 0.056}, // 3.0-3.2 (fixed)
 {1.038, 0.067}, // 3.2-4.7
 {1.111, 0.042}, // 0.0-1.3 (new)
};

const double kparul18[_nres+1][2] = {
// JER SF produced with minitools/resolution.C:redoJER("RunUL18")
 {1.144, 0.032}, //0.010}, // 0.0-0.5
 {1.151, 0.034}, // 0.5-1.0
 {1.146, 0.040}, // 1.0-1.5
 {1.159, 0.032}, // 1.5-2.0
 {1.144, 0.087}, // 2.0-2.5
 {1.271, 0.089}, // 2.5-3.0
 {1.238, 0.072}, // 3.0-3.2
 {1.037, 0.157}, // 3.2-4.7
 {1.145, 0.032}, // 0.0-1.3
};
// NB: 0.0-0.5 uncertainty suspiciously low compared to UL17 and other bins
// => changed to 0.0-1.3 value by hand
// NB2: all JER SF ~1.15 at |y|<2.0. Could be Pythia8 modelling (vs Herwig++)

const double kparul16bcd[_nres+1][2] = {
// JER SF produced with minitools/resolution.C:redoJER("RunUL16BCD")
 {1.091, 0.023}, // 0.0-0.5
 {1.097, 0.019}, // 0.5-1.0
 {1.069, 0.030}, // 1.0-1.5
 {1.035, 0.030}, // 1.5-2.0
 {1.013, 0.040}, // 2.0-2.5
 {1.151, 0.083}, // 2.5-3.0
 {1.173, 0.036}, // 3.0-3.2
 {1.006, 0.046}, // 3.2-4.7
 {1.089, 0.023}, // 0.0-1.3
};

const double kparul16ef[_nres+1][2] = {
// JER SF produced with minitools/resolution.C:redoJER("RunUL16EF")
 {1.091, 0.023}, // 0.0-0.5
 {1.097, 0.019}, // 0.5-1.0
 {1.069, 0.030}, // 1.0-1.5
 {1.035, 0.030}, // 1.5-2.0
 {1.013, 0.040}, // 2.0-2.5
 {1.150, 0.083}, // 2.5-3.0
 {1.173, 0.036}, // 3.0-3.2
 {1.006, 0.046}, // 3.2-4.7
 {1.089, 0.023}, // 0.0-1.3
};

const double kparul16apv[_nres+1][2] = {
// JER SF produced with minitools/resolution.C:redoJER("RunUL16BCDEF")
 {1.091, 0.023}, // 0.0-0.5
 {1.097, 0.019}, // 0.5-1.0
 {1.069, 0.030}, // 1.0-1.5
 {1.035, 0.030}, // 1.5-2.0
 {1.013, 0.040}, // 2.0-2.5
 {1.150, 0.083}, // 2.5-3.0
 {1.173, 0.036}, // 3.0-3.2
 {1.006, 0.046}, // 3.2-4.7
 {1.089, 0.023}, // 0.0-1.3
};

const double kparul16gh[_nres+1][2] = {
// JER SF produced with minitools/resolution.C:redoJER("RunUL16GH")
 {1.099, 0.013}, // 0.0-0.5
 {1.112, 0.029}, // 0.5-1.0
 {1.088, 0.055}, // 1.0-1.5
 {1.064, 0.043}, // 1.5-2.0
 {1.043, 0.046}, // 2.0-2.5
 {1.171, 0.139}, // 2.5-3.0
 {1.146, 0.034}, // 3.0-3.2
 {1.067, 0.045}, // 3.2-4.7
 {1.102, 0.033}, // 0.0-1.3
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

// Redid UL17 fits a second time to fix kinematic cutoff bias in HF
// (Gaussian fit range startig a bit higher for 3.2-4.7)

const double vparul17[_nres+1][3] =
// Fit of JER for R=0.4, 13 TeV file MC-UL17V4_BCDEF
  {{-2.64, 1.012, 0.0352},  // y 0.0-0.5, chi2 46.0/55
   {-2.42, 1.022, 0.0403},  // y 0.5-1.0, chi2 36.5/54
   {-2.77, 1.212, 0.0518},  // y 1.0-1.5, chi2 84.7/49
   {1.07, 1.158, 0.0339},  // y 1.5-2.0, chi2 21.6/41
   {2.15, 1.138, 0.0323},  // y 2.0-2.5, chi2 14.6/33
   {4.35, 1.169, 0.0115},  // y 2.5-3.0, chi2 11.1/26
   {4.78, 0.871, 0.1405},  // y 3.0-3.2, chi2 31.8/20
   {3.44, 0.659, 0.1033},  // y 3.2-4.7, chi2 28.5/18
   {-2.95, 1.080, 0.0372}}; // y 0.0-1.3, chi2 40.7/54

const double vparul17b[_nres+1][3] =
// Fit of JER for R=0.4, 13 TeV file MC-UL17V4_B
  {{-2.84, 1.004, 0.0353},  // y 0.0-0.5, chi2 56.2/55
   {-2.56, 1.012, 0.0403},  // y 0.5-1.0, chi2 36.0/54
   {-2.86, 1.196, 0.0519},  // y 1.0-1.5, chi2 78.1/49
   {0.68, 1.140, 0.0341},  // y 1.5-2.0, chi2 23.5/41
   {2.17, 1.102, 0.0330},  // y 2.0-2.5, chi2 17.6/33
   {4.13, 1.148, 0.0176},  // y 2.5-3.0, chi2 10.1/26
   {3.20, 1.073, 0.1279},  // y 3.0-3.2, chi2 88.9/20
   {2.98, 0.733, 0.0991},  // y 3.2-4.7, chi2 46.5/18
   {-3.06, 1.070, 0.0373}}; // y 0.0-1.3, chi2 53.1/54

const double vparul17c[_nres+1][3] =
// Fit of JER for R=0.4, 13 TeV file MC-UL17V4_C
  {{-2.85, 1.003, 0.0353},  // y 0.0-0.5, chi2 56.5/55
   {-2.58, 1.011, 0.0404},  // y 0.5-1.0, chi2 36.2/54
   {-2.88, 1.195, 0.0519},  // y 1.0-1.5, chi2 76.0/49
   {0.55, 1.139, 0.0341},  // y 1.5-2.0, chi2 23.1/41
   {2.17, 1.098, 0.0332},  // y 2.0-2.5, chi2 17.7/33
   {4.02, 1.154, 0.0161},  // y 2.5-3.0, chi2 9.9/26
   {3.65, 1.012, 0.1313},  // y 3.0-3.2, chi2 69.8/20
   {2.30, 0.832, 0.0932},  // y 3.2-4.7, chi2 112.1/18
   {-3.08, 1.069, 0.0373}}; // y 0.0-1.3, chi2 53.9/54

const double vparul17d[_nres+1][3] =
// Fit of JER for R=0.4, 13 TeV file MC-UL17V4_D
  {{-2.85, 1.002, 0.0353},  // y 0.0-0.5, chi2 57.3/55
   {-2.58, 1.012, 0.0403},  // y 0.5-1.0, chi2 36.4/54
   {-2.86, 1.193, 0.0520},  // y 1.0-1.5, chi2 70.5/49
   {0.63, 1.138, 0.0341},  // y 1.5-2.0, chi2 24.7/41
   {2.22, 1.095, 0.0334},  // y 2.0-2.5, chi2 16.8/33
   {4.10, 1.143, 0.0187},  // y 2.5-3.0, chi2 11.3/26
   {3.04, 1.076, 0.1281},  // y 3.0-3.2, chi2 112.4/20
   {3.02, 0.724, 0.0994},  // y 3.2-4.7, chi2 42.2/18
   {-3.07, 1.068, 0.0373}}; // y 0.0-1.3, chi2 52.3/54

const double vparul17e[_nres+1][3] =
// Fit of JER for R=0.4, 13 TeV file MC-UL17V4_F
  {{-2.44, 1.020, 0.0351},  // y 0.0-0.5, chi2 41.8/55
   {-2.24, 1.035, 0.0402},  // y 0.5-1.0, chi2 39.1/54
   {-2.50, 1.226, 0.0517},  // y 1.0-1.5, chi2 80.4/49
   {1.45, 1.181, 0.0338},  // y 1.5-2.0, chi2 21.6/41
   {2.18, 1.179, 0.0316},  // y 2.0-2.5, chi2 14.5/33
   {4.70, 1.181, 0.0001},  // y 2.5-3.0, chi2 14.3/26
   {4.66, 0.975, 0.1341},  // y 3.0-3.2, chi2 42.3/20
   {3.82, 0.576, 0.1055},  // y 3.2-4.7, chi2 38.4/18
   {-2.74, 1.090, 0.0371}}; // y 0.0-1.3, chi2 35.2/54

const double vparul17f[_nres+1][3] =
// Fit of JER for R=0.4, 13 TeV file MC-UL17V4_F
  {{-2.31, 1.024, 0.0350},  // y 0.0-0.5, chi2 33.9/55
   {-2.33, 1.035, 0.0402},  // y 0.5-1.0, chi2 38.5/54
   {-2.61, 1.225, 0.0517},  // y 1.0-1.5, chi2 81.7/49
   {1.20, 1.183, 0.0333},  // y 1.5-2.0, chi2 28.1/41
   {2.18, 1.176, 0.0309},  // y 2.0-2.5, chi2 20.8/33
   {4.76, 1.173, -0.0000},  // y 2.5-3.0, chi2 12.4/26
   {4.90, 0.898, 0.1412},  // y 3.0-3.2, chi2 38.0/20
   {4.57, -0.000, 0.1106}, // y 3.2-4.7, chi2 112.8/18
   {-2.74, 1.090, 0.0371}}; // y 0.0-1.3, chi2 35.2/54 (new)

// UL18 fits iterated once, even though many fits got slightly worse chi2
// Starting point was UL17 JER

const double vparul18[_nres+1][3] =
// Fit of JER for R=0.4, 13 TeV file MC-UL18V2V3_ABCD
  {{-2.03, 1.006, 0.0352},  // y 0.0-0.5, chi2 38.5/55
   {-2.05, 1.033, 0.0395},  // y 0.5-1.0, chi2 53.7/54
   {-2.07, 1.213, 0.0512},  // y 1.0-1.5, chi2 129.3/49
   {3.03, 1.081, 0.0425},  // y 1.5-2.0, chi2 24.3/41
   {3.68, 1.052, 0.0395},  // y 2.0-2.5, chi2 18.7/33
   {5.74, 1.106, 0.0276},  // y 2.5-3.0, chi2 9.6/26
   {5.01, 1.132, 0.1391},  // y 3.0-3.2, chi2 30.3/20
   {3.34, 0.820, 0.0965},  // y 3.2-4.7, chi2 37.5/18
   {-2.53, 1.083, 0.0368}}; // y 0.0-1.3, chi2 90.9/54

const double vparul18a[_nres+1][3] =
// Fit of JER for R=0.4, 13 TeV file MC-UL18V2V3_A
  {{-1.83, 1.007, 0.0352},  // y 0.0-0.5, chi2 39.3/55
   {-1.80, 1.035, 0.0395},  // y 0.5-1.0, chi2 48.0/54
   {-1.86, 1.218, 0.0512},  // y 1.0-1.5, chi2 133.5/49
   {3.19, 1.090, 0.0425},  // y 1.5-2.0, chi2 24.3/41
   {3.88, 1.063, 0.0398},  // y 2.0-2.5, chi2 19.2/33
   {5.94, 1.107, 0.0274},  // y 2.5-3.0, chi2 9.8/26
   {5.14, 1.145, 0.1390},  // y 3.0-3.2, chi2 27.7/20
   {3.17, 0.871, 0.0936},  // y 3.2-4.7, chi2 83.3/18
   {-2.37, 1.085, 0.0368}}; // y 0.0-1.3, chi2 86.2/54

const double vparul18b[_nres+1][3] =
// Fit of JER for R=0.4, 13 TeV file MC-UL18V2V3_B
  {{-2.01, 1.006, 0.0352},  // y 0.0-0.5, chi2 38.4/55
   {-2.03, 1.033, 0.0395},  // y 0.5-1.0, chi2 53.8/54
   {-2.05, 1.213, 0.0512},  // y 1.0-1.5, chi2 129.7/49
   {3.04, 1.081, 0.0425},  // y 1.5-2.0, chi2 24.3/41
   {3.71, 1.052, 0.0395},  // y 2.0-2.5, chi2 18.7/33
   {5.75, 1.107, 0.0275},  // y 2.5-3.0, chi2 9.8/26
   {5.01, 1.132, 0.1392},  // y 3.0-3.2, chi2 30.1/20
   {3.34, 0.823, 0.0964},  // y 3.2-4.7, chi2 40.0/18
   {-2.51, 1.082, 0.0368}}; // y 0.0-1.3, chi2 90.5/54

const double vparul18c[_nres+1][3] =
// Fit of JER for R=0.4, 13 TeV file MC-UL18V2V3_C
  {{-1.93, 1.005, 0.0353},  // y 0.0-0.5, chi2 39.2/55
   {-1.87, 1.032, 0.0395},  // y 0.5-1.0, chi2 51.0/54
   {-1.91, 1.211, 0.0516},  // y 1.0-1.5, chi2 121.5/49
   {2.91, 1.094, 0.0426},  // y 1.5-2.0, chi2 20.8/41
   {3.66, 1.068, 0.0394},  // y 2.0-2.5, chi2 23.5/33
   {5.72, 1.133, 0.0264},  // y 2.5-3.0, chi2 8.4/26
   {5.61, 1.112, 0.1445},  // y 3.0-3.2, chi2 13.4/20
   {3.00, 0.890, 0.0923},  // y 3.2-4.7, chi2 88.4/18
   {-2.43, 1.083, 0.0369}}; // y 0.0-1.3, chi2 91.2/54

const double vparul18d[_nres+1][3] =
// Fit of JER for R=0.4, 13 TeV file MC-UL18V2V3_D
  {{-2.13, 1.005, 0.0352},  // y 0.0-0.5, chi2 37.5/55
   {-2.14, 1.031, 0.0395},  // y 0.5-1.0, chi2 52.9/54
   {-2.20, 1.211, 0.0515},  // y 1.0-1.5, chi2 123.8/49
   {2.73, 1.088, 0.0426},  // y 1.5-2.0, chi2 20.6/41
   {3.41, 1.064, 0.0389},  // y 2.0-2.5, chi2 24.0/33
   {5.47, 1.139, 0.0260},  // y 2.5-3.0, chi2 8.4/26
   {5.41, 1.103, 0.1448},  // y 3.0-3.2, chi2 13.1/20
   {3.23, 0.828, 0.0963},  // y 3.2-4.7, chi2 36.3/18
   {-2.60, 1.082, 0.0368}}; // y 0.0-1.3, chi2 93.3/54

const double vparul16bcd[_nres+1][3] =
// Fit of JER for R=0.4, 13 TeV file MC-UL16V3V1_BCD
  {{-2.77, 1.000, 0.0356},  // y 0.0-0.5, chi2 57.8/55
   {-2.56, 1.007, 0.0405},  // y 0.5-1.0, chi2 32.7/54
   {-1.57, 1.131, 0.0535},  // y 1.0-1.5, chi2 62.0/49
   {1.84, 1.161, 0.0320},  // y 1.5-2.0, chi2 28.3/41
   {2.70, 1.057, 0.0327},  // y 2.0-2.5, chi2 29.6/33
   {2.93, 1.116, 0.0262},  // y 2.5-3.0, chi2 51.9/26
   {4.11, 0.514, 0.1336},  // y 3.0-3.2, chi2 256.9/20
   {3.28, 0.649, 0.0977},  // y 3.2-4.7, chi2 42.0/18
   {-3.01, 1.061, 0.0374}}; // y 0.0-1.3, chi2 64.9/54

const double vparul16ef[_nres+1][3] =
// Fit of JER for R=0.4, 13 TeV file MC-UL16V3V1_EF
  {{-2.43, 1.042, 0.0356},  // y 0.0-0.5, chi2 49.0/55
   {-2.11, 1.046, 0.0407},  // y 0.5-1.0, chi2 39.2/54
   {0.91, 1.170, 0.0546},  // y 1.0-1.5, chi2 88.5/49
   {2.31, 1.209, 0.0334},  // y 1.5-2.0, chi2 33.4/41
   {3.31, 1.097, 0.0368},  // y 2.0-2.5, chi2 44.6/33
   {3.96, 1.140, 0.0264},  // y 2.5-3.0, chi2 26.4/26
   {5.49, 0.000, 0.1311},  // y 3.0-3.2, chi2 197.6/20
   {3.13, 0.772, 0.0927},  // y 3.2-4.7, chi2 76.4/18
   {-2.63, 1.103, 0.0375}}; // y 0.0-1.3, chi2 67.7/54

const double vparul16apv[_nres+1][3] =
// Fit of JER for R=0.4, 13 TeV file MC-UL16V3V1_BCDEF
  {{-2.66, 1.019, 0.0356},  // y 0.0-0.5, chi2 54.5/55
   {-2.36, 1.025, 0.0406},  // y 0.5-1.0, chi2 33.6/54
   {-1.11, 1.150, 0.0540},  // y 1.0-1.5, chi2 74.4/49
   {2.03, 1.180, 0.0330},  // y 1.5-2.0, chi2 26.8/41
   {2.91, 1.074, 0.0349},  // y 2.0-2.5, chi2 33.1/33
   {3.47, 1.133, 0.0261},  // y 2.5-3.0, chi2 33.2/26
   {5.42, 0.000, 0.1348},  // y 3.0-3.2, chi2 122.0/20
   {3.03, 0.746, 0.0948},  // y 3.2-4.7, chi2 43.2/18
   {-2.85, 1.080, 0.0374}}; // y 0.0-1.3, chi2 66.6/54

const double vparul16gh[_nres+1][3] =
// Fit of JER for R=0.4, 13 TeV file MC-UL16V2V1_GH
  {{-2.08, 1.000, 0.0354},  // y 0.0-0.5, chi2 52.8/55
   {-1.36, 1.004, 0.0406},  // y 0.5-1.0, chi2 35.6/54
   {1.44, 1.160, 0.0529},  // y 1.0-1.5, chi2 72.3/49
   {2.64, 1.199, 0.0334},  // y 1.5-2.0, chi2 32.5/41
   {3.51, 1.114, 0.0351},  // y 2.0-2.5, chi2 46.4/33
   {4.03, 1.071, 0.0325},  // y 2.5-3.0, chi2 163.7/26
   {4.83, -0.000, 0.1419},  // y 3.0-3.2, chi2 202.0/20
   {2.57, 0.885, 0.0827},  // y 3.2-4.7, chi2 202.1/18
   {-2.21, 1.065, 0.0372}}; // y 0.0-1.3, chi2 74.8/54

const double vjesul17[_nres+1][3] =
// Fit of JES for R=0.4, 13 TeV file MC-UL17V4_BCDEF
  {{0.0417, 2.708, 0.3844},  // y 0.0-0.5, chi2 31.4/6
   {0.0519, 2.708, 0.3844},  // y 0.5-1.0, chi2 13.3/6
   {0.0708, 2.708, 0.3851},  // y 1.0-1.5, chi2 2.2/6
   {0.0726, 2.708, 0.4110},  // y 1.5-2.0, chi2 2.7/6
   {0.0787, 2.708, 0.3978},  // y 2.0-2.5, chi2 6.1/6
   {0.1848, 2.708, 0.3350},  // y 2.5-3.0, chi2 12.1/6
   {0.2564, 2.708, 0.4779},  // y 3.0-3.2, chi2 26.3/6
   {0.1285, 2.708, 0.4126},  // y 3.2-4.7, chi2 57.5/6
   {0.0472, 2.708, 0.3934}}; // y 0.0-1.3, chi2 22.6/6

const double vjesul17b[_nres+1][3] =
// Fit of JES for R=0.4, 13 TeV file MC-UL17V4_B
  {{0.0383, 2.708, 0.3593},  // y 0.0-0.5, chi2 23.1/6
   {0.0479, 2.708, 0.3533},  // y 0.5-1.0, chi2 9.7/6
   {0.0601, 2.708, 0.3616},  // y 1.0-1.5, chi2 3.1/6
   {0.0591, 2.708, 0.3686},  // y 1.5-2.0, chi2 4.0/6
   {0.0692, 2.708, 0.3446},  // y 2.0-2.5, chi2 8.1/6
   {0.1709, 2.708, 0.3206},  // y 2.5-3.0, chi2 17.0/6
   {0.2335, 2.708, 0.4743},  // y 3.0-3.2, chi2 24.8/6
   {0.0871, 2.708, 0.4502},  // y 3.2-4.7, chi2 52.7/6
   {0.0439, 2.708, 0.3616}}; // y 0.0-1.3, chi2 16.0/6

const double vjesul17c[_nres+1][3] =
// Fit of JES for R=0.4, 13 TeV file MC-UL17V4_C
  {{0.0372, 2.708, 0.3603},  // y 0.0-0.5, chi2 23.3/6
   {0.0470, 2.708, 0.3510},  // y 0.5-1.0, chi2 9.9/6
   {0.0573, 2.708, 0.3567},  // y 1.0-1.5, chi2 3.8/6
   {0.0563, 2.708, 0.3576},  // y 1.5-2.0, chi2 6.7/6
   {0.0680, 2.708, 0.3333},  // y 2.0-2.5, chi2 11.2/6
   {0.1673, 2.708, 0.3183},  // y 2.5-3.0, chi2 18.7/6
   {0.2248, 2.708, 0.4840},  // y 3.0-3.2, chi2 26.4/6
   {0.0795, 2.708, 0.4669},  // y 3.2-4.7, chi2 50.1/6
   {0.0431, 2.708, 0.3584}}; // y 0.0-1.3, chi2 14.2/6

const double vjesul17d[_nres+1][3] =
// Fit of JES for R=0.4, 13 TeV file MC-UL17V4_D
  {{0.0376, 2.708, 0.3580},  // y 0.0-0.5, chi2 20.2/6
   {0.0477, 2.708, 0.3508},  // y 0.5-1.0, chi2 8.1/6
   {0.0599, 2.708, 0.3567},  // y 1.0-1.5, chi2 3.0/6
   {0.0575, 2.708, 0.3577},  // y 1.5-2.0, chi2 5.0/6
   {0.0685, 2.708, 0.3360},  // y 2.0-2.5, chi2 9.2/6
   {0.1691, 2.708, 0.3167},  // y 2.5-3.0, chi2 19.0/6
   {0.2256, 2.708, 0.4778},  // y 3.0-3.2, chi2 23.5/6
   {0.0830, 2.708, 0.4532},  // y 3.2-4.7, chi2 50.9/6
   {0.0440, 2.708, 0.3572}}; // y 0.0-1.3, chi2 13.1/6

const double vjesul17e[_nres+1][3] =
// Fit of JES for R=0.4, 13 TeV file MC-UL17V4_E
  {{0.0459, 2.708, 0.3891},  // y 0.0-0.5, chi2 25.8/6
   {0.0571, 2.708, 0.3880},  // y 0.5-1.0, chi2 12.8/6
   {0.0800, 2.708, 0.3946},  // y 1.0-1.5, chi2 1.5/6
   {0.0854, 2.708, 0.4255},  // y 1.5-2.0, chi2 2.2/6
   {0.0886, 2.708, 0.4247},  // y 2.0-2.5, chi2 6.0/6
   {0.1975, 2.708, 0.3410},  // y 2.5-3.0, chi2 8.5/6
   {0.2663, 2.708, 0.4767},  // y 3.0-3.2, chi2 25.0/6
   {0.1624, 2.708, 0.3862},  // y 3.2-4.7, chi2 59.5/6
   {0.0522, 2.708, 0.3998}}; // y 0.0-1.3, chi2 23.7/6

const double vjesul17f[_nres+1][3] =
// Fit of JES for R=0.4, 13 TeV file MC-UL17V4_F
  {{0.0524, 2.708, 0.3932},  // y 0.0-0.5, chi2 16.8/6
   {0.0613, 2.708, 0.4083},  // y 0.5-1.0, chi2 11.6/6
   {0.0909, 2.708, 0.4004},  // y 1.0-1.5, chi2 1.5/6
   {0.0991, 2.708, 0.4380},  // y 1.5-2.0, chi2 2.7/6
   {0.1001, 2.708, 0.4384},  // y 2.0-2.5, chi2 3.6/6
   {0.2077, 2.708, 0.3496},  // y 2.5-3.0, chi2 5.1/6
   {0.2814, 2.708, 0.4813},  // y 3.0-3.2, chi2 25.4/6
   {0.1812, 2.708, 0.3959},  // y 3.2-4.7, chi2 56.6/6
   {0.0578, 2.708, 0.4131}}; // y 0.0-1.3, chi2 22.8/6

const double vjesul18[_nres+1][3] =
// Fit of JES for R=0.4, 13 TeV file MC-UL18V2V3_ABCD
  {{-0.0098, 2.708, -0.6398},  // y 0.0-0.5, chi2 14.8/6
   {-0.1260, 2.708, -0.1465},  // y 0.5-1.0, chi2 18.7/6
   {-0.0151, 2.708, -0.4094},  // y 1.0-1.5, chi2 5.0/6
   {-0.2272, 2.708, 0.0464},  // y 1.5-2.0, chi2 19.0/6
   {-0.0196, 2.708, 0.5397},  // y 2.0-2.5, chi2 8.0/6
   {-0.0206, 2.708, 0.7535},  // y 2.5-3.0, chi2 26.9/6
   {0.0357, 2.708, 0.7886},  // y 3.0-3.2, chi2 8.6/6
   {-0.0500, 2.708, 0.2247},  // y 3.2-4.7, chi2 30.0/6
   {-0.0090, 2.708, -0.5652}}; // y 0.0-1.3, chi2 13.7/6

const double vjesul18a[_nres+1][3] =
// Fit of JES for R=0.4, 13 TeV file MC-UL18V2V3_A
  {{-0.0118, 2.708, -0.5764},  // y 0.0-0.5, chi2 13.2/6
   {-0.0107, 2.708, -0.5594},  // y 0.5-1.0, chi2 8.8/6
   {-0.0180, 2.708, 0.3913},  // y 1.0-1.5, chi2 5.6/6
   {-0.0128, 2.708, -0.4021},  // y 1.5-2.0, chi2 15.6/6
   {-0.0230, 2.708, -0.4616},  // y 2.0-2.5, chi2 9.7/6
   {-0.0194, 2.708, 0.7128},  // y 2.5-3.0, chi2 27.0/6
   {0.0344, 2.708, 0.7835},  // y 3.0-3.2, chi2 6.5/6
   {-0.0460, 2.708, 0.2407},  // y 3.2-4.7, chi2 28.1/6
   {-0.0130, 2.708, 0.4726}}; // y 0.0-1.3, chi2 16.5/6

const double vjesul18b[_nres+1][3] =
// Fit of JES for R=0.4, 13 TeV file MC-UL18V2V3_B
  {{-0.0099, 2.708, -0.6382},  // y 0.0-0.5, chi2 14.5/6
   {-0.1248, 2.708, -0.1485},  // y 0.5-1.0, chi2 18.6/6
   {-0.0156, 2.708, 0.4040},  // y 1.0-1.5, chi2 5.0/6
   {-0.2222, 2.708, 0.0465},  // y 1.5-2.0, chi2 19.4/6
   {-0.0199, 2.708, -0.5302},  // y 2.0-2.5, chi2 8.1/6
   {-0.0202, 2.708, 0.7593},  // y 2.5-3.0, chi2 27.2/6
   {0.0367, 2.708, 0.7775},  // y 3.0-3.2, chi2 8.6/6
   {-0.0495, 2.708, 0.2268},  // y 3.2-4.7, chi2 29.8/6
   {-0.0095, 2.708, 0.5493}}; // y 0.0-1.3, chi2 16.2/6

const double vjesul18c[_nres+1][3] =
// Fit of JES for R=0.4, 13 TeV file MC-UL18V2V3_C
  {{-0.0103, 2.708, 0.5792},  // y 0.0-0.5, chi2 12.4/6
   {-0.0962, 2.708, -0.1511},  // y 0.5-1.0, chi2 16.2/6
   {0.0022, 2.708, 0.0094},  // y 1.0-1.5, chi2 20.1/6
   {-0.0106, 2.708, -0.4534},  // y 1.5-2.0, chi2 10.4/6
   {-0.0210, 2.708, 0.5652},  // y 2.0-2.5, chi2 9.1/6
   {-0.0259, 2.708, 0.9157},  // y 2.5-3.0, chi2 12.3/6
   {0.0891, 2.708, 0.3182},  // y 3.0-3.2, chi2 25.3/6
   {-0.0407, 2.708, 0.2180},  // y 3.2-4.7, chi2 36.8/6
   {-0.0097, 2.708, 0.5094}}; // y 0.0-1.3, chi2 16.2/6

const double vjesul18d[_nres+1][3] =
// Fit of JES for R=0.4, 13 TeV file MC-UL18V2V3_D
  {{-0.0076, 2.708, -0.6849},  // y 0.0-0.5, chi2 15.9/6
   {-0.0063, 2.708, 0.7091},  // y 0.5-1.0, chi2 11.2/6
   {-0.0136, 2.708, -0.4184},  // y 1.0-1.5, chi2 3.4/6
   {-0.4531, 2.708, 0.0454},  // y 1.5-2.0, chi2 20.5/6
   {-0.0185, 2.708, 0.7128},  // y 2.0-2.5, chi2 9.7/6
   {-0.0272, 2.708, 0.9444},  // y 2.5-3.0, chi2 13.6/6
   {0.0826, 2.708, 0.3502},  // y 3.0-3.2, chi2 34.0/6
   {-0.0506, 2.708, 0.1898},  // y 3.2-4.7, chi2 39.9/6
   {-0.0066, 2.708, -0.6222}}; // y 0.0-1.3, chi2 15.2/6

const double vjesul16bcd[_nres+1][3] =
// Fit of JES for R=0.4, 13 TeV file MC-UL16V3V1_BCD
  {{-0.0127, 2.708, 1.1655},  // y 0.0-0.5, chi2 4.0/6
   {-0.0064, 2.708, 1.4175},  // y 0.5-1.0, chi2 8.9/6
   {-0.0098, 2.708, 0.5450},  // y 1.0-1.5, chi2 5.6/6
   {-0.0256, 2.708, 0.3660},  // y 1.5-2.0, chi2 7.2/6
   {-0.0232, 2.708, 0.6011},  // y 2.0-2.5, chi2 10.9/6
   {-0.0233, 2.708, 0.4010},  // y 2.5-3.0, chi2 30.8/6
   {2.0378, 2.708, -0.0438},  // y 3.0-3.2, chi2 139.0/6
   {-0.0569, 2.708, 0.1399},  // y 3.2-4.7, chi2 31.2/6
   {-0.0088, 2.708, 1.0792}}; // y 0.0-1.3, chi2 16.2/6

const double vjesul16ef[_nres+1][3] =
// Fit of JES for R=0.4, 13 TeV file MC-UL16V3V1_EF
//  {{-0.0105, 2.708, 1046.3662},  // y 0.0-0.5, chi2 16.1/6 (bad?)
//   {-0.0049, 2.708, 715.0220},  // y 0.5-1.0, chi2 11.7/6 (bad?)
  {{-0.0130, 2.708, 1.2000},  // y 0.0-0.5, chi2 23.2/6 (redone)
   {-0.0061, 2.708, 1.2000},  // y 0.5-1.0, chi2 13.2/6 (redone)
   {-0.0092, 2.708, 1.0555},  // y 1.0-1.5, chi2 20.3/6
   {-0.0203, 2.708, 0.6954},  // y 1.5-2.0, chi2 29.4/6
   {-0.0336, 2.708, 0.6085},  // y 2.0-2.5, chi2 12.4/6
   {-11.2324, 2.708, 0.0070},  // y 2.5-3.0, chi2 77.4/6
   {0.2079, 2.708, 0.3865},  // y 3.0-3.2, chi2 108.6/6
   //{0.0093, 2.708, 1026.2131},  // y 3.2-4.7, chi2 32.3/6 (bad?)
   //{-0.0070, 2.708, 1897.2814}}; // y 0.0-1.3, chi2 20.6/6 (bad?)
   {0.0097, 2.708, 1.2000},  // y 3.2-4.7, chi2 30.9/6 (redone)
   {-0.0086, 2.708, 1.2000}}; // y 0.0-1.3, chi2 26.2/6 (redone)

const double vjesul16apv[_nres+1][3] =
// Fit of JES for R=0.4, 13 TeV file MC-UL16V3V1_BCDEF
  {{-0.0108, 2.708, 2.3927},  // y 0.0-0.5, chi2 8.2/6
   {-0.0060, 2.708, 1.5840},  // y 0.5-1.0, chi2 6.3/6
   {-0.0097, 2.708, 0.6702},  // y 1.0-1.5, chi2 8.9/6
   {-0.0221, 2.708, 0.4745},  // y 1.5-2.0, chi2 11.9/6
   {-0.0268, 2.708, 0.6092},  // y 2.0-2.5, chi2 12.4/6
   {-0.0069, 2.708, 0.3911},  // y 2.5-3.0, chi2 41.9/6
   //{0.0060, 2.708, 1704.5604},  // y 3.2-4.7, chi2 21.6/6 (bad?)
   //{0.0065, 2.708, 1535.2922},  // y 3.2-4.7, chi2 27.6/6 (bad?)
   //{-0.0079, 2.708, 1.6986}}; // y 0.0-1.3, chi2 11.2/6
   {0.0200, 2.708, 1.2000},  // y 3.0-3.2, chi2 213.1/6 (redone)
   {-0.0200, 2.708, 0.1330},  // y 3.2-4.7, chi2 48.4/6 (redone)
   {-0.0085, 2.708, 1.2000}}; // y 0.0-1.3, chi2 17.7/6 (redone)

const double vjesul16gh[_nres+1][3] =
// Fit of JES for R=0.4, 13 TeV file MC-UL16V2V1_GH
  {{-0.0154, 2.708, -0.7738},  // y 0.0-0.5, chi2 10.7/6
   {-0.0091, 2.708, 1.4061},  // y 0.5-1.0, chi2 16.5/6
   {-0.0154, 2.708, 0.7654},  // y 1.0-1.5, chi2 7.5/6
   {-0.0314, 2.708, 0.5431},  // y 1.5-2.0, chi2 6.9/6
   {-0.0515, 2.708, 0.4929},  // y 2.0-2.5, chi2 8.1/6
   {-0.0060, 2.708, 0.3588},  // y 2.5-3.0, chi2 158.0/6
   {0.2630, 2.708, -0.3512},  // y 3.0-3.2, chi2 100.4/6
   {0.0728, 2.708, 0.1893},  // y 3.2-4.7, chi2 49.0/6
   {-0.0114, 2.708, 0.9130}}; // y 0.0-1.3, chi2 22.9/6


double ptresolution(double pt, double eta, double eta2=0) {

  if (eta>3.2 && eta<3.5) eta = 3.55; // 3.2-4.7 bin instead of 3.5-4.7

  //int ieta = min(_nres-1, int(fabs(eta) / 0.5));
  //int iy = min(5, ieta);
  //int iy = min(_nres-1, int(fabs(eta) / 0.5 + 0.5));
  int iy = min(_nres-1, int(fabs(eta) / 0.5));
  //if (fabs(eta)>3.2) iy = _nres-1; // 3.2-4.7 bin instead of 3.5-4.7
  if (eta<0.5 && eta2>1.0) int iy = _nres; // 0.0-1.3 (UL17 && UL18)
  // Only supported for UL17 for now
  if (iy==_nres) {
    assert(_jer_iov==ul17 || _jer_iov==ul17b || _jer_iov==ul17c || 
	   _jer_iov==ul17d || _jer_iov==ul17e || _jer_iov==ul17f ||
	   _jer_iov==ul18 || _jer_iov==ul18a || _jer_iov==ul18b || 
	   _jer_iov==ul18c || _jer_iov==ul18d ||
	   _jer_iov==ul16apv || _jer_iov==ul16bcd || _jer_iov==ul16ef ||
	   _jer_iov==ul16gh);
  }
  double res = 0;

  // New scaling for 2018 JER SF
  if (!_newsf) {
    _newsf = new TF1("SFnew3","sqrt(pow([0],2)/(x*x)+pow(sqrt([3])*[1],2)/x+pow([3]*[2],2)) / sqrt(pow([3]*[0],2)/(x*x)+pow([3]*[1],2)/x+pow([3]*[2],2)) * [3]",10,3000);
    _newsf->SetParameters(5*2, 0, 0.05*sqrt(2), 2.020);
  } // newsf

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
      double jer_sf = kpar2018[iy][0];
      if (_usenewsf) {
	_newsf->SetParameter(3, jer_sf); jer_sf = _newsf->Eval(pt);
      }
      if (!_ismcjer) res *= jer_sf;
    }
    if (_jer_iov==run2018abc) {
      res = sqrt(vpar2018[iy][0]*fabs(vpar2018[iy][0])/(pt*pt) +
		 pow(vpar2018[iy][1],2)/pt + 
		 pow(vpar2018[iy][2],2));
      double jer_sf = kpar2018abc[iy][0];
      if (_usenewsf) {
	_newsf->SetParameter(3, jer_sf); jer_sf = _newsf->Eval(pt);
      }
      if (!_ismcjer) res *= jer_sf;
    }
    if (_jer_iov==run2018d) {
      res = sqrt(vpar2018[iy][0]*fabs(vpar2018[iy][0])/(pt*pt) +
		 pow(vpar2018[iy][1],2)/pt + 
		 pow(vpar2018[iy][2],2));
      double jer_sf = kpar2018d[iy][0];
      if (_usenewsf) {
	_newsf->SetParameter(3, jer_sf); jer_sf = _newsf->Eval(pt);
      }
      if (!_ismcjer) res *= jer_sf;
    }
    // UL17
    if (_jer_iov==ul17) {
      res = sqrt(vparul17[iy][0]*fabs(vparul17[iy][0])/(pt*pt) +
		 pow(vparul17[iy][1],2)/pt + 
		 pow(vparul17[iy][2],2));
      if (!_ismcjer) res *= kparul17[iy][0];
    }
    if (_jer_iov==ul17b) {
      res = sqrt(vparul17b[iy][0]*fabs(vparul17b[iy][0])/(pt*pt) +
		 pow(vparul17b[iy][1],2)/pt + 
		 pow(vparul17b[iy][2],2));
      if (!_ismcjer) res *= kparul17[iy][0];
    }
    if (_jer_iov==ul17c) {
      res = sqrt(vparul17c[iy][0]*fabs(vparul17c[iy][0])/(pt*pt) +
		 pow(vparul17c[iy][1],2)/pt + 
		 pow(vparul17c[iy][2],2));
      if (!_ismcjer) res *= kparul17[iy][0];
    }
    if (_jer_iov==ul17d) {
      res = sqrt(vparul17d[iy][0]*fabs(vparul17d[iy][0])/(pt*pt) +
		 pow(vparul17d[iy][1],2)/pt + 
		 pow(vparul17d[iy][2],2));
      if (!_ismcjer) res *= kparul17[iy][0];
    }
    if (_jer_iov==ul17e) {
      res = sqrt(vparul17e[iy][0]*fabs(vparul17e[iy][0])/(pt*pt) +
		 pow(vparul17e[iy][1],2)/pt + 
		 pow(vparul17e[iy][2],2));
      if (!_ismcjer) res *= kparul17[iy][0];
    }
    if (_jer_iov==ul17f) {
      res = sqrt(vparul17f[iy][0]*fabs(vparul17f[iy][0])/(pt*pt) +
		 pow(vparul17f[iy][1],2)/pt + 
		 pow(vparul17f[iy][2],2));
      if (!_ismcjer) res *= kparul17[iy][0];
    }
    // UL18
    if (_jer_iov==ul18) {
      res = sqrt(vparul18[iy][0]*fabs(vparul18[iy][0])/(pt*pt) +
		 pow(vparul18[iy][1],2)/pt + 
		 pow(vparul18[iy][2],2));
      if (!_ismcjer) res *= kparul18[iy][0];
    }
    if (_jer_iov==ul18a) {
      res = sqrt(vparul18a[iy][0]*fabs(vparul18a[iy][0])/(pt*pt) +
		 pow(vparul18a[iy][1],2)/pt + 
		 pow(vparul18a[iy][2],2));
      if (!_ismcjer) res *= kparul18[iy][0];
    }
    if (_jer_iov==ul18b) {
      res = sqrt(vparul18b[iy][0]*fabs(vparul18b[iy][0])/(pt*pt) +
		 pow(vparul18b[iy][1],2)/pt + 
		 pow(vparul18b[iy][2],2));
      if (!_ismcjer) res *= kparul18[iy][0];
    }
    if (_jer_iov==ul18c) {
      res = sqrt(vparul18c[iy][0]*fabs(vparul18c[iy][0])/(pt*pt) +
		 pow(vparul18c[iy][1],2)/pt + 
		 pow(vparul18c[iy][2],2));
      if (!_ismcjer) res *= kparul18[iy][0];
    }
    if (_jer_iov==ul18d) {
      res = sqrt(vparul18d[iy][0]*fabs(vparul18d[iy][0])/(pt*pt) +
		 pow(vparul18d[iy][1],2)/pt + 
		 pow(vparul18d[iy][2],2));
      if (!_ismcjer) res *= kparul18[iy][0];
    }

    // UL16
    if (_jer_iov==ul16bcd) {
      res = sqrt(vparul16bcd[iy][0]*fabs(vparul16bcd[iy][0])/(pt*pt) +
		 pow(vparul16bcd[iy][1],2)/pt + 
		 pow(vparul16bcd[iy][2],2));
      if (!_ismcjer) res *= kparul16bcd[iy][0];
    }
    if (_jer_iov==ul16ef) {
      res = sqrt(vparul16ef[iy][0]*fabs(vparul16ef[iy][0])/(pt*pt) +
		 pow(vparul16ef[iy][1],2)/pt + 
		 pow(vparul16ef[iy][2],2));
      if (!_ismcjer) res *= kparul16ef[iy][0];
    }
    if (_jer_iov==ul16apv) {
      res = sqrt(vparul16apv[iy][0]*fabs(vparul16apv[iy][0])/(pt*pt) +
		 pow(vparul16apv[iy][1],2)/pt + 
		 pow(vparul16apv[iy][2],2));
      if (!_ismcjer) res *= kparul16apv[iy][0];
    } // placeholder
    if (_jer_iov==ul16gh) {
      res = sqrt(vparul16gh[iy][0]*fabs(vparul16gh[iy][0])/(pt*pt) +
		 pow(vparul16gh[iy][1],2)/pt + 
		 pow(vparul16gh[iy][2],2));
      if (!_ismcjer) res *= kparul16gh[iy][0];
    }

  } // !_usejme

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
      string resolutionFile = "../JRDatabase/textFiles/Autumn18_V4_MC/"
      "Autumn18_V4_MC_PtResolution_AK4PFchs.txt";
      string scaleFactorFile = "../JRDatabase/textFiles/Autumn18_V4_MC/"
	"Autumn18_V4_MC_SF_AK4PFchs.txt";
      string weightFile = "rootfiles/jerweights.root";

      _jer = new JME::JetResolution(resolutionFile);
      _jer_sf = new JME::JetResolutionScaleFactor(scaleFactorFile);
    }
    if (_jer_iov==run2018abc && !_jer && !_jer_sf) {
      string resolutionFile = "../JRDatabase/textFiles/Autumn18_RunABC_V4_MC/"
      "Autumn18_RunABC_V4_MC_PtResolution_AK4PFchs.txt";
      string scaleFactorFile = "../JRDatabase/textFiles/Autumn18_RunABC_V4_MC/"
	"Autumn18_RunABC_V4_MC_SF_AK4PFchs.txt";
      string weightFile = "rootfiles/jerweights.root";

      _jer = new JME::JetResolution(resolutionFile);
      _jer_sf = new JME::JetResolutionScaleFactor(scaleFactorFile);
    }
    if (_jer_iov==run2018d && !_jer && !_jer_sf) {
      string resolutionFile = "../JRDatabase/textFiles/Autumn18_RunD_V4_MC/"
      "Autumn18_RunD_V4_MC_PtResolution_AK4PFchs.txt";
      string scaleFactorFile = "../JRDatabase/textFiles/Autumn18_RunD_V4_MC/"
	"Autumn18_RunD_V4_MC_SF_AK4PFchs.txt";
      string weightFile = "rootfiles/jerweights.root";

      _jer = new JME::JetResolution(resolutionFile);
      _jer_sf = new JME::JetResolutionScaleFactor(scaleFactorFile);
    }
    if ((_jer_iov==ul17 || _jer_iov==ul17b || _jer_iov==ul17c
	 || _jer_iov==ul17d || _jer_iov==ul17e || _jer_iov==ul17f)
	&& !_jer && !_jer_sf) {
      string resolutionFile = "../JRDatabase/textFiles/Summer19UL17_JRV2_MC/"
	"Summer19UL17_JRV2_MC_PtResolution_AK4PFchs.txt";
      string scaleFactorFile = "../JRDatabase/textFiles/Summer19UL17_JRV2_MC/"
	"Summer19UL17_JRV2_MC_SF_AK4PFchs.txt";
      string weightFile = "rootfiles/jerweights.root";

      _jer = new JME::JetResolution(resolutionFile);
      _jer_sf = new JME::JetResolutionScaleFactor(scaleFactorFile);
    }
    if ((_jer_iov==ul18 || _jer_iov==ul18a || _jer_iov==ul18b
	 || _jer_iov==ul18c || _jer_iov==ul18d)
	&& !_jer && !_jer_sf) {
      string resolutionFile = "../JRDatabase/textFiles/Summer19UL18_JRV2_MC/"
	"Summer19UL18_JRV2_MC_PtResolution_AK4PFchs.txt";
      string scaleFactorFile = "../JRDatabase/textFiles/Summer19UL18_JRV2_MC/"
	"Summer19UL18_JRV2_MC_SF_AK4PFchs.txt";
      string weightFile = "rootfiles/jerweights.root";

      _jer = new JME::JetResolution(resolutionFile);
      _jer_sf = new JME::JetResolutionScaleFactor(scaleFactorFile);
    }
    if ((_jer_iov==ul16apv || _jer_iov==ul16bcd || _jer_iov==ul16ef)
	&& !_jer && !_jer_sf) {
      string resolutionFile = "../JRDatabase/textFiles/Summer20UL16_JRV1_MC/"
	"Summer20UL16_JRV1_MC_PtResolution_AK4PFchs.txt"; // placeholder
      //string scaleFactorFile = "../JRDatabase/textFiles/Summer19UL16_JRV1_MC/"
      //"Summer19UL16_JRV1_MC_SF_AK4PFchs.txt";
      string scaleFactorFile = "../JERCProtoLab/Summer19UL16/JER_SF/Summer19UL16APV_JRV1_MC_SF_AK4PFchs.txt"; // APV
      string weightFile = "rootfiles/jerweights.root";

      _jer = new JME::JetResolution(resolutionFile);
      _jer_sf = new JME::JetResolutionScaleFactor(scaleFactorFile);
    }
    if (_jer_iov==ul16gh && !_jer && !_jer_sf) {
      string resolutionFile = "../JRDatabase/textFiles/Summer20UL16_JRV1_MC/"
	"Summer20UL16_JRV1_MC_PtResolution_AK4PFchs.txt"; // placeholder
      //string scaleFactorFile = "../JRDatabase/textFiles/Summer19UL16_JRV1_MC/"
      //"Summer19UL16_JRV1_MC_SF_AK4PFchs.txt";
      string scaleFactorFile = "../JERCProtoLab/Summer19UL16/JER_SF/Summer19UL16_JRV1_MC_SF_AK4PFchs.txt"; // Non-APV
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
    
    if (_jer_iov==run2018 || _jer_iov==run2018abc || _jer_iov==run2018d) {
      if (_usenewsf) {
	_newsf->SetParameter(3, jer_sf); jer_sf = _newsf->Eval(pt);
      }
    }

    res = jet_resolution;
    if (!_ismcjer) res *= jer_sf;
  }

  return res;
} // ptresolution

// return bias in low pT JES
double ptresponse(double pt, double eta, double eta2=0) {

  if (eta>3.2 && eta<3.5) eta = 3.55; // 3.2-4.7 bin instead of 3.5-4.7
  int iy = min(_nres-1, int(fabs(eta) / 0.5));
  //if (fabs(eta)>3.2) iy = _nres-1; // 3.2-4.7 bin instead of 3.5-4.7
  if (eta<0.5 && eta2>1.0) int iy = _nres; // 0.0-1.3 (UL17)
  // Only supported for UL17 for now
  if (iy==_nres) {
    assert(_jer_iov==ul17 || _jer_iov==ul17b || _jer_iov==ul17c || 
	   _jer_iov==ul17d || _jer_iov==ul17e || _jer_iov==ul17f);
  }

  const double *p;
  p = &vjesul17[iy][0];
  if (_jer_iov==ul17)  p = &vjesul17[iy][0];
  if (_jer_iov==ul17b) p = &vjesul17b[iy][0];
  if (_jer_iov==ul17c) p = &vjesul17c[iy][0];
  if (_jer_iov==ul17d) p = &vjesul17d[iy][0];
  if (_jer_iov==ul17e) p = &vjesul17e[iy][0];
  if (_jer_iov==ul17f) p = &vjesul17f[iy][0];

  double jes = 1+p[0]*exp(-0.5*pow(log(pt)-p[1],2)/(p[2]*p[2]));

  p = &vjesul17[iy][0];
  double jes0 = 1+p[0]*exp(-0.5*pow(log(pt)-p[1],2)/(p[2]*p[2]));

  return (jes / jes0);
} // ptresponse


// retun ECAL prefire fraction
//enum ecal_iov {run2016, run2017, run2018};
TF1 *_fecalpf(0);
double ecalprefire(double pt, double eta, jer_iov run) {

  if (!_fecalpf) {
    _fecalpf = new TF1("feff3","[0]*0.5*(1+erf((x-[1])/([2]*sqrt(x))))",
		       55,4000./cosh(eta));
    _fecalpf->SetParameters(0,0,1);
  }

  // fits done with jecsys/minitools/resolution.C:redoECALprefire()
  const int neta = 2;
  const int npar = 3;
  const double pars16[neta][npar] = {
    //{0.205, 324.0, 31.4}, // 2016BtoH, eta 2.0-2.5
    {0.148, 205.0, 17.1}, // 2016BtoH, eta 2.0-2.5
    //{0.588, 227.7, 12.3}, // 2016BtoH, eta 2.5-3.0
    {0.582, 225.5, 12.1}, // 2016BtoH, eta 2.5-3.0
  };
  const double pars16bcd[neta][npar] = {
  //{0.091, 233.4, 22.9}, // 2016BCD, eta 2.0-2.5
    {0.082, 205.0, 17.1}, // 2016BCD, eta 2.0-2.5
    //{0.522, 261.5, 15.7}, // 2016BCD, eta 2.5-3.0
    {0.447, 225.8, 12.1}, // 2016BCD, eta 2.5-3.0
  };
  const double pars16ef[neta][npar] = {
    //{0.120, 154.8, 8.6}, // 2016EF, eta 2.0-2.5
    {0.156, 205.0, 17.1}, // 2016EF, eta 2.0-2.5
    //{0.534, 237.7, 11.8}, // 2016EF, eta 2.5-3.0
    {0.541, 240.3, 12.1}, // 2016EF, eta 2.5-3.0
  };
  const double pars16gh[neta][npar] = {
    //{0.162, 204.2, 17.1}, // 2016FGH, eta 2.0-2.5
    {0.163, 205.0, 17.1}, // 2016FGH, eta 2.0-2.5
    //{0.699, 196.4, 12.1}, // 2016FGH, eta 2.5-3.0
    {0.699, 196.3, 12.1}, // 2016FGH, eta 2.5-3.0
  };
  //
  const double pars17[neta][npar] = {
    //{0.521, 543.5, 52.0}, // 2017BtoF, eta 2.0-2.5
    {0.243, 175.0, 13.4}, // 2017BtoF, eta 2.0-2.5
    //{0.837, 191.6, 11.3}, // 2017BtoF, eta 2.5-3.0
    {0.780, 179.5, 9.4}, // 2017BtoF, eta 2.5-3.0

  };
  const double pars17b[neta][npar] = {
    //{0.110, 164.7, 11.0}, // 2017B, eta 2.0-2.5
    {0.117, 175.0, 13.4}, // 2017B, eta 2.0-2.5
    //{0.770, 219.0, 13.2}, // 2017B, eta 2.5-3.0
    {0.664, 191.2, 9.4}, // 2017B, eta 2.5-3.0
  };
  const double pars17c[neta][npar] = {
    //{0.391, 556.8, 49.4}, // 2017C, eta 2.0-2.5
    {0.172, 175.0, 13.4}, // 2017C, eta 2.0-2.5
    //{0.729, 166.2, 9.3}, // 2017C, eta 2.5-3.0
    {0.733, 166.9, 9.4}, // 2017C, eta 2.5-3.0
  };
  const double pars17de[neta][npar] = {
    //{0.230, 172.2, 13.4}, // 2017DE, eta 2.0-2.5
    {0.232, 175.0, 13.4}, // 2017DE, eta 2.0-2.5
    //{0.767, 160.1, 9.4}, // 2017DE, eta 2.5-3.0
    {0.768, 160.3, 9.4}, // 2017DE, eta 2.5-3.0
  };
  const double pars17f[neta][npar] = {
    //{0.615, 442.8, 41.7}, // 2017F, eta 2.0-2.5
    {0.316, 175.0, 13.4}, // 2017F, eta 2.0-2.5
    //{0.850, 164.1, 9.5}, // 2017F, eta 2.5-3.0
    {0.846, 163.4, 9.4}, // 2017F, eta 2.5-3.0
  };

// resolution.C:redoECALprefire(2,2017BtoF)
  const double parsul17[neta][npar] = {
  {0.243, 175.0, 13.4}, // 2017BtoF, eta 2.0-2.5 (72.1/5)
  {0.780, 179.5, 9.4}, // 2017BtoF, eta 2.5-3.0 (17.3/4)
  };
  const double parsul17b[neta][npar] = {
  {0.116, 175.0, 13.4}, // 2017B, eta 2.0-2.5 (2.5/5)
  {0.647, 186.2, 9.4}, // 2017B, eta 2.5-3.0 (5.9/4)
  };
  const double parsul17c[neta][npar] = {
  {0.174, 175.0, 13.4}, // 2017C, eta 2.0-2.5 (3.1/5)
  {0.735, 168.3, 9.4}, // 2017C, eta 2.5-3.0 (1.2/4)
  };
  const double parsul17d[neta][npar] = {
  {0.201, 175.0, 13.4}, // 2017D, eta 2.0-2.5 (4.1/5)
  {0.746, 161.7, 9.4}, // 2017D, eta 2.5-3.0 (4.3/4)
  };
  const double parsul17e[neta][npar] = {
  {0.246, 175.0, 13.4}, // 2017E, eta 2.0-2.5 (2.7/5)
  {0.779, 160.3, 9.4}, // 2017E, eta 2.5-3.0 (1.9/4)
  };
  const double parsul17f[neta][npar] = {
  {0.321, 175.0, 13.4}, // 2017F, eta 2.0-2.5 (13.4/5)
  {0.847, 163.9, 9.4}, // 2017F, eta 2.5-3.0 (1.8/4)
  };
  
  if (fabs(eta)<2.0 || fabs(eta)>3.0) return 0;
  if (run==run2018) return 0;
  if (run==run2018d) return 0;
  if (run==run2018abc) return 0;

  int iy = (fabs(eta)<2.5 ? 0 : 1);
  const double (*pars)[neta][npar](0);// = (run==run2016 ? pars16 : pars17);
  if (run==run2016)    pars = &pars16;
  if (run==run2016bcd) pars = &pars16bcd;
  if (run==run2016ef)  pars = &pars16ef;
  if (run==run2016gh)  pars = &pars16gh;
  if (run==run2017)    pars = &pars17;
  if (run==run2017b)   pars = &pars17b;
  if (run==run2017c)   pars = &pars17c;
  if (run==run2017de)  pars = &pars17de;
  if (run==run2017f)   pars = &pars17f;
  assert(pars);
  _fecalpf->SetParameters((*pars)[iy][0],(*pars)[iy][1],(*pars)[iy][2]);

  return (_fecalpf->Eval(pt));
}

#endif // __ptresolution_h__
