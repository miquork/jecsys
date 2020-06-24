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
	       ul17, ul17b, ul17c, ul17d, ul17e, ul17f};
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

const double vparul17[_nres+1][3] =
// Fit of JER for R=0.4, 13 TeV file MC-UL17V4_BCDEF
  {{-2.54, 1.016, 0.0350},  // y 0.0-0.5, chi2 36.5/55
   {-2.49, 1.023, 0.0402},  // y 0.5-1.0, chi2 35.7/54
   {-2.81, 1.212, 0.0518},  // y 1.0-1.5, chi2 84.2/49
   {0.70, 1.164, 0.0332},  // y 1.5-2.0, chi2 25.8/41
   {2.19, 1.134, 0.0320},  // y 2.0-2.5, chi2 18.8/33
   {4.37, 1.165, 0.0101},  // y 2.5-3.0, chi2 11.8/26
   {5.03, 0.791, 0.1470},  // y 3.0-3.2, chi2 37.6/20
   {4.52, -0.000, 0.1114}, // y 3.2-4.7, chi2 88.1/18
   {-2.95, 1.080, 0.0372}}; // y 0.0-1.3, chi2 40.7/54 (new)

const double vparul17b[_nres+1][3] =
// Fit of JER for R=0.4, 13 TeV file MC-UL17V4_B
  {{-2.68, 1.009, 0.0351},  // y 0.0-0.5, chi2 42.6/55
   {-2.61, 1.014, 0.0403},  // y 0.5-1.0, chi2 35.9/54
   {-2.89, 1.196, 0.0519},  // y 1.0-1.5, chi2 77.9/49
   {-0.32, 1.145, 0.0334},  // y 1.5-2.0, chi2 24.7/41
   {2.23, 1.097, 0.0333},  // y 2.0-2.5, chi2 18.3/33
   {4.17, 1.144, 0.0169},  // y 2.5-3.0, chi2 11.8/26
   {4.98, 0.734, 0.1477},  // y 3.0-3.2, chi2 22.9/20
   {4.48, 0.000, 0.1099}, // y 3.2-4.7, chi2 99.3/18
   {-3.06, 1.070, 0.0373}}; // y 0.0-1.3, chi2 53.1/54 (new)

const double vparul17c[_nres+1][3] =
// Fit of JER for R=0.4, 13 TeV file MC-UL17V4_C
  {{-2.70, 1.008, 0.0351},  // y 0.0-0.5, chi2 42.5/55
   {-2.63, 1.013, 0.0403},  // y 0.5-1.0, chi2 36.6/54
   {-2.90, 1.195, 0.0519},  // y 1.0-1.5, chi2 75.8/49
   {-0.46, 1.143, 0.0334},  // y 1.5-2.0, chi2 27.4/41
   {2.21, 1.094, 0.0333},  // y 2.0-2.5, chi2 18.1/33
   {4.05, 1.151, 0.0150},  // y 2.5-3.0, chi2 12.3/26
   {5.09, 0.698, 0.1490},  // y 3.0-3.2, chi2 24.7/20
   {4.48, 0.000, 0.1100}, // y 3.2-4.7, chi2 89.7/18
   {-3.08, 1.069, 0.0373}}; // y 0.0-1.3, chi2 53.9/54 (new)

const double vparul17d[_nres+1][3] =
// Fit of JER for R=0.4, 13 TeV file MC-UL17V4_D
  {{-2.69, 1.007, 0.0351},  // y 0.0-0.5, chi2 43.5/55
   {-2.62, 1.013, 0.0402},  // y 0.5-1.0, chi2 35.9/54
   {-2.88, 1.193, 0.0520},  // y 1.0-1.5, chi2 69.9/49
   {-0.31, 1.142, 0.0335},  // y 1.5-2.0, chi2 25.7/41
   {2.26, 1.091, 0.0336},  // y 2.0-2.5, chi2 16.4/33
   {4.13, 1.141, 0.0177},  // y 2.5-3.0, chi2 13.3/26
   {5.16, 0.673, 0.1488},  // y 3.0-3.2, chi2 20.8/20
   {4.46, -0.000, 0.1098}, // y 3.2-4.7, chi2 103.1/18
   {-3.07, 1.068, 0.0373}}; // y 0.0-1.3, chi2 52.3/54 (new)

const double vparul17e[_nres+1][3] =
// Fit of JER for R=0.4, 13 TeV file MC-UL17V4_E
  {{-2.45, 1.019, 0.0350},  // y 0.0-0.5, chi2 32.2/55
   {-2.42, 1.026, 0.0402},  // y 0.5-1.0, chi2 38.9/54
   {-2.75, 1.218, 0.0517},  // y 1.0-1.5, chi2 84.8/49
   {1.00, 1.168, 0.0333},  // y 1.5-2.0, chi2 28.9/41
   {2.26, 1.146, 0.0314},  // y 2.0-2.5, chi2 19.4/33
   {4.50, 1.170, 0.0066},  // y 2.5-3.0, chi2 12.0/26
   {4.82, 0.855, 0.1449},  // y 3.0-3.2, chi2 37.6/20
   {4.54, 0.000, 0.1112}, // y 3.2-4.7, chi2 117.6/18
   {-2.87, 1.084, 0.0372}}; // y 0.0-1.3, chi2 39.2/54 (new)

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


const double vjesul17[_nres+1][3] =
// Fit of JES for R=0.4, 13 TeV file MC-UL17V4_BCDEF
  {{0.0417, 2.708, 0.3844},  // y 0.0-0.5, chi2 31.4/6
   {0.0519, 2.708, 0.3844},  // y 0.5-1.0, chi2 13.3/6
   {0.0708, 2.708, 0.3851},  // y 1.0-1.5, chi2 2.2/6
   {0.0726, 2.708, 0.4110},  // y 1.5-2.0, chi2 2.7/6
   {0.0787, 2.708, 0.3978},  // y 2.0-2.5, chi2 6.1/6
   {0.1848, 2.708, 0.3350},  // y 2.5-3.0, chi2 12.1/6
   {0.2564, 2.708, 0.4779},  // y 3.0-3.2, chi2 26.3/6
   {0.1503, 2.708, 0.3853}, // y 3.2-4.7, chi2 54.2/6
   {0.0472, 2.708, 0.3934}}; // y 0.0-1.3, chi2 22.6/6 (new)

const double vjesul17b[_nres+1][3] =
// Fit of JES for R=0.4, 13 TeV file MC-UL17V4_B
  {{0.0411, 2.708, 0.3460},  // y 0.0-0.5, chi2 14.8/6
   {0.0479, 2.708, 0.3533},  // y 0.5-1.0, chi2 9.7/6
   {0.0601, 2.708, 0.3616},  // y 1.0-1.5, chi2 3.1/6
   {0.0591, 2.708, 0.3686},  // y 1.5-2.0, chi2 4.0/6
   {0.0692, 2.708, 0.3446},  // y 2.0-2.5, chi2 8.1/6
   {0.1709, 2.708, 0.3206},  // y 2.5-3.0, chi2 17.0/6
   {0.2335, 2.708, 0.4743},  // y 3.0-3.2, chi2 24.8/6
   {0.1150, 2.708, 0.3914}, // y 3.2-4.7, chi2 55.3/6
   {0.0439, 2.708, 0.3616}}; // y 0.0-1.3, chi2 16.0/6 (new)

const double vjesul17c[_nres+1][3] =
// Fit of JES for R=0.4, 13 TeV file MC-UL17V4_C
  {{0.0397, 2.708, 0.3479},  // y 0.0-0.5, chi2 14.7/6
   {0.0470, 2.708, 0.3510},  // y 0.5-1.0, chi2 9.9/6
   {0.0573, 2.708, 0.3567},  // y 1.0-1.5, chi2 3.8/6
   {0.0563, 2.708, 0.3576},  // y 1.5-2.0, chi2 6.7/6
   {0.0680, 2.708, 0.3333},  // y 2.0-2.5, chi2 11.2/6
   {0.1673, 2.708, 0.3183},  // y 2.5-3.0, chi2 18.7/6
   {0.2248, 2.708, 0.4840},  // y 3.0-3.2, chi2 26.4/6
   {0.1051, 2.708, 0.4045}, // y 3.2-4.7, chi2 54.7/6
   {0.0431, 2.708, 0.3584}}; // y 0.0-1.3, chi2 14.2/6 (new)

const double vjesul17d[_nres+1][3] =
// Fit of JES for R=0.4, 13 TeV file MC-UL17V4_D
  {{0.0403, 2.708, 0.3440},  // y 0.0-0.5, chi2 12.9/6
   {0.0477, 2.708, 0.3508},  // y 0.5-1.0, chi2 8.1/6
   {0.0599, 2.708, 0.3567},  // y 1.0-1.5, chi2 3.0/6
   {0.0575, 2.708, 0.3577},  // y 1.5-2.0, chi2 5.0/6
   {0.0685, 2.708, 0.3360},  // y 2.0-2.5, chi2 9.2/6
   {0.1691, 2.708, 0.3167},  // y 2.5-3.0, chi2 19.0/6
   {0.2256, 2.708, 0.4778},  // y 3.0-3.2, chi2 23.5/6
   {0.1107, 2.708, 0.3906}, // y 3.2-4.7, chi2 56.0/6
   {0.0440, 2.708, 0.3572}}; // y 0.0-1.3, chi2 13.1/6 (new)

const double vjesul17e[_nres+1][3] =
// Fit of JES for R=0.4, 13 TeV file MC-UL17V4_E
  {{0.0477, 2.708, 0.3823},  // y 0.0-0.5, chi2 17.0/6
   {0.0571, 2.708, 0.3880},  // y 0.5-1.0, chi2 12.8/6
   {0.0800, 2.708, 0.3946},  // y 1.0-1.5, chi2 1.5/6
   {0.0854, 2.708, 0.4255},  // y 1.5-2.0, chi2 2.2/6
   {0.0886, 2.708, 0.4247},  // y 2.0-2.5, chi2 6.0/6
   {0.1975, 2.708, 0.3410},  // y 2.5-3.0, chi2 8.5/6
   {0.2663, 2.708, 0.4767},  // y 3.0-3.2, chi2 25.0/6
   {0.1765, 2.708, 0.3755}, // y 3.2-4.7, chi2 52.4/6
   {0.0522, 2.708, 0.3998}}; // y 0.0-1.3, chi2 23.7/6 (new)

const double vjesul17f[_nres+1][3] =
// Fit of JES for R=0.4, 13 TeV file MC-UL17V4_F
  {{0.0529, 2.708, 0.3929},  // y 0.0-0.5, chi2 13.4/6
   {0.0613, 2.708, 0.4083},  // y 0.5-1.0, chi2 11.6/6
   {0.0909, 2.708, 0.4004},  // y 1.0-1.5, chi2 1.5/6
   {0.0991, 2.708, 0.4380},  // y 1.5-2.0, chi2 2.7/6
   {0.1001, 2.708, 0.4384},  // y 2.0-2.5, chi2 3.6/6
   {0.2077, 2.708, 0.3496},  // y 2.5-3.0, chi2 5.1/6
   {0.2814, 2.708, 0.4813},  // y 3.0-3.2, chi2 25.4/6
   {0.1965, 2.708, 0.3849}, // y 3.2-4.7, chi2 51.0/6
   {0.0578, 2.708, 0.4131}}; // y 0.0-1.3, chi2 22.8/6 (new)


double ptresolution(double pt, double eta, double eta2=0) {

  if (eta>3.2 && eta<3.5) eta = 3.55; // 3.2-4.7 bin instead of 3.5-4.7

  //int ieta = min(_nres-1, int(fabs(eta) / 0.5));
  //int iy = min(5, ieta);
  //int iy = min(_nres-1, int(fabs(eta) / 0.5 + 0.5));
  int iy = min(_nres-1, int(fabs(eta) / 0.5));
  //if (fabs(eta)>3.2) iy = _nres-1; // 3.2-4.7 bin instead of 3.5-4.7
  if (eta<0.5 && eta2>1.0) int iy = _nres; // 0.0-1.3 (UL17)
  // Only supported for UL17 for now
  if (iy==_nres) {
    assert(_jer_iov==ul17 || _jer_iov==ul17b || _jer_iov==ul17c || 
	   _jer_iov==ul17d || _jer_iov==ul17e || _jer_iov==ul17f);
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
    //
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
