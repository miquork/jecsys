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
bool _ak7 = false;
bool _ak4 = true;

// This is needed by drawRunHistos to unfold JER from rate vs NPV
//const double _spu = 1.5; // sigma_PU/PV = 1.5 GeV
//const double _spu = 2.6; // sigma_PU/PV = 1.5 GeV
//const double _npv0 = 13.45;
//double _npv = 13.45;

/*
// Data / MC ratio for JER, produced with ak7ak5ratio.C:redoJER() (2012/08/29)
const double kpar[6][2] = {
  {1.052, 0.063}, // 0.0-0.5
  {1.057, 0.057}, // 0.5-1.0
  {1.088, 0.063}, // 1.0-1.5
  {1.117, 0.082}, // 1.5-2.0
  {1.189, 0.134}, // 2.0-2.5
  {1.288, 0.200} // 2.5-3.
};
*/
/*
// From Sanmay's talk (Aug 6)
const double kpar[6][2] = {
  {1.12, 0.063},
  {1.11, 0.057},
  {1.11, 0.063},
  {1.22, 0.082},
  {1.23, 0.134},
  {1.1,  0.200}
};
*/

// From Sanmay by e-mail (May 8, 2014)
/*
const double kpar[6][2] = {
  {1.12, 0.063},
  {1.11, 0.057},
  {1.11, 0.063},
  {1.22, 0.082},
  {1.23, 0.134},
  {1.1,  0.200}
  //{1.21,  0.200}
};
*/
// Data / MC ratio for JER, produced with ak5ak7ratio12.C:redoJER() (2014/09/15)
// Based on averaging recommended numbers from 2012
/*
const double kpar[6][2] = {
 {1.079, 0.026}, // 0.0-0.5
 {1.099, 0.028}, // 0.5-1.0
 {1.116, 0.029}, // 1.0-1.5
 {1.169, 0.039}, // 1.5-2.0
 {1.227, 0.053}, // 2.0-2.5
 {1.296, 0.062}, // 2.5-3.0
};
*/
// From Sanmay by email 30 Mar 2015
// (direct recommendations from different eta bins => not ideal?)
const double kpar[6][2] = {
  {1.079, 0.026},
  {1.099, 0.028},
  {1.121, 0.029},
  {1.208, 0.039},
  {1.208, 0.053},
  {1.395, 0.062},
};

// Resolutions with Pythia Z2* + ak5ak7resolution12.C
// (produced on iMac desktop, with ROOT 5.30/00, iterating with AK5+2sigma)
// On Sep 15, 2014, using Winter14_V5 private pre-version (root tuples v12)
// Fit of JER for R=0.5, 8 TeV, 53X
const double vpar5[6][3] =
  {{3.13, 0.897, 0.0337},  // y 0.0-0.5, chi2 21.6/33
   {3.58, 0.868, 0.0376},  // y 0.5-1.0, chi2 12.9/33
   {4.78, 0.885, 0.0438},  // y 1.0-1.5, chi2 26.5/33
   {5.36, 0.745, 0.0265},  // y 1.5-2.0, chi2 14.5/32
   {4.45, 0.680, 0.0253},  // y 2.0-2.5, chi2 10.5/25
   {2.86, 0.874, 0.0000}}; // y 2.5-3.0, chi2 10.8/19
/*
  // values from Sanmay by email 30 Mar 2015
  // duplicate AK7 just to check right MC truth is picked up
  {{5.79356, 0.984061, 0.0290218},
   {6.10575, 0.952320, 0.0328014},
   {6.34643, 0.998053, 0.0370081},
   {6.49438, 0.866373, 0.00}, //8.97208e-07},
   {6.06794, 0.734516, 0.00}, //1.31400e-05},
   {4.57993, 0.853656, 0.00}};//1.30360e-06}};
*/
// Fit of JER for R=0.7, 8 TeV, 53X
const double vpar7[6][3] =
  // values from Sanmay by email 30 Mar 2015
  {{5.79356, 0.984061, 0.0290218},
   {6.10575, 0.952320, 0.0328014},
   {6.34643, 0.998053, 0.0370081},
   {6.49438, 0.866373, 0.00}, //8.97208e-07},
   {6.06794, 0.734516, 0.00}, //1.31400e-05},
   {4.57993, 0.853656, 0.00}};//1.30360e-06}};
  /*
  // values from Sanmays AN2012_223_v19.pdf Table 5 page 20.
  {{6.251, 1.062, 0.031},  // 0-0.5
   {6.710, 1.047, 0.036},  // 0.5-1
   {7.114, 1.119, 0.041},  // 1-1.5
   {7.845, 1.047, 0.000},  // 1.5-2
   {7.330, 0.887, 0.000},  // 2-2.5
   {6.389, 1.191, 0.000}}; // 2.5-3
  */
  /*
  // values from Sanmays AN2012_223_v19.pdf Table 5 page 20.
  {{6.251/kpar[0][0], 1.062/kpar[0][0], 0.031/kpar[0][0]},  // 0-0.5
   {6.710/kpar[1][0], 1.047/kpar[1][0], 0.036/kpar[1][0]},  // 0.5-1
   {7.114/kpar[2][0], 1.119/kpar[2][0], 0.041/kpar[2][0]},  // 1-1.5
   {7.845/kpar[3][0], 1.047/kpar[3][0], 0.000/kpar[3][0]},  // 1.5-2
   {7.330/kpar[4][0], 0.887/kpar[4][0], 0.000/kpar[4][0]},  // 2-2.5
   {6.389/kpar[5][0], 1.191/kpar[5][0], 0.000/kpar[5][0]}}; // 2.5-3
  */
  /*
  {{5.47, 0.896, 0.0329},  // y 0.0-0.5, chi2 17.8/33
   {5.67, 0.870, 0.0365},  // y 0.5-1.0, chi2 14.2/33
   {6.72, 0.848, 0.0435},  // y 1.0-1.5, chi2 29.7/33
   {6.83, 0.726, 0.0245},  // y 1.5-2.0, chi2 14.2/32
   {5.92, 0.660, 0.0235},  // y 2.0-2.5, chi2 11.0/25
   {3.92, 0.862, 0.0000}}; // y 2.5-3.0, chi2 29.5/19
  */

// Resolutions with Pythia Z2 + ak5ak7resolution12.C
// (produced on iMac desktop, with ROOT 5.30/00, iterating with AK5+2sigma)
// Fit of JER for R=0.5, 8 TeV, 53X (2014/08/11)
/*
  {{3.20, 0.907, 0.0343},  // y 0.0-0.5, chi2 16.7/33
   {3.58, 0.882, 0.0374},  // y 0.5-1.0, chi2 17.0/33
   {4.58, 0.913, 0.0433},  // y 1.0-1.5, chi2 22.8/33
   {5.49, 0.739, 0.0256},  // y 1.5-2.0, chi2 18.6/32
   {4.73, 0.678, 0.0234},  // y 2.0-2.5, chi2 19.2/25
   {1.32, 0.916, 0.0000}}; // y 2.5-3.0, chi2 7.0/19
*/
// Fit of JER for R=0.5, 8 TeV, 53X (20/fb, April 10)
/*
  {{1.88, 0.933, 0.0304},  // y 0.0-0.5, chi2 18.8/33
   {2.44, 0.921, 0.0339},  // y 0.5-1.0, chi2 40.4/33
   {2.23, 1.008, 0.0356},  // y 1.0-1.5, chi2 61.1/33
   {4.12, 0.832, 0.0179},  // y 1.5-2.0, chi2 23.5/32
   {4.79, 0.644, 0.0250},  // y 2.0-2.5, chi2 16.2/25
   {3.77, 0.819, 0.0090}}; // y 2.5-3.0, chi2 9.6/19
*/
// Fit of JER for 8 TeV, 9/fb, 53X V5 (2012/12/11)
/*
  {{2.70, 0.888, 0.0336},  // y 0.0-0.5, chi2 24.9/35
   {2.95, 0.883, 0.0365},  // y 0.5-1.0, chi2 17.2/34
   {4.49, 0.886, 0.0434},  // y 1.0-1.5, chi2 18.1/31
   {4.96, 0.757, 0.0246},  // y 1.5-2.0, chi2 10.6/27
   {4.56, 0.650, 0.0254},  // y 2.0-2.5, chi2 9.5/21
   {3.96, 0.812, 0.0138}}; // y 2.5-3.0, chi2 3.0/16
*/
/*
  // 8 TeV, 9/fb, 53X Fall12 V1 (2012/10/25)
  {{2.63, 0.885, 0.0338},  // y 0.0-0.5, chi2 26.7/35
   {2.91, 0.879, 0.0366},  // y 0.5-1.0, chi2 17.3/34
   {4.49, 0.880, 0.0437},  // y 1.0-1.5, chi2 18.6/31
   {4.92, 0.754, 0.0250},  // y 1.5-2.0, chi2 13.2/27
   {4.42, 0.649, 0.0254},  // y 2.0-2.5, chi2 10.0/21
   {3.63, 0.825, 0.0091}}; // y 2.5-3.0, chi2 5.2/16
*/
  /*
  // 7 TeV, 52X GRV23 (2012/08/29)
  {{0.83, 0.845, 0.0369},  // y 0.0-0.5, chi2 29.6/35
   {0.98, 0.858, 0.0379},  // y 0.5-1.0, chi2 24.1/34
   {2.81, 0.894, 0.0435},  // y 1.0-1.5, chi2 19.7/31
   {4.14, 0.732, 0.0248},  // y 1.5-2.0, chi2 15.7/27
   {3.67, 0.623, 0.0247},  // y 2.0-2.5, chi2 8.5/21
   {3.40, 0.782, 0.0109}}; // y 2.5-3.0, chi2 12.9/16
  */
// Fit of JER for AK7
// Fit of JER for R=0.7, 8 TeV, 53X  (2014/08/11)
/*
  {{5.82, 0.893, 0.0337},  // y 0.0-0.5, chi2 15.5/33
   {5.78, 0.879, 0.0364},  // y 0.5-1.0, chi2 15.8/33
   {6.72, 0.869, 0.0433},  // y 1.0-1.5, chi2 24.3/33
   {7.14, 0.702, 0.0249},  // y 1.5-2.0, chi2 21.5/32
   {5.70, 0.675, 0.0221},  // y 2.0-2.5, chi2 6.9/25
   {3.93, 0.834, -0.0000}}; // y 2.5-3.0, chi2 2.6/19
*/
  // Parameter values from Sanmay by e-mail (May 8, 2014)
  /* // older version?
  {{6.13, 0.950, 0.0308},
   {6.82, 0.896, 0.0351},
   {7.32, 0.940, 0.0393},
   {7.73, 0.737, 0.0218568},
   {8.00, 0.590, 0.0231838},
   {6.85, 0.821, 5.21952e-08}};
  */
/* // newer version?
  {{6.13049, 0.949756, 0.0308104},
   {6.81515, 0.896017, 0.0350882},
   {7.32160, 0.939579, 0.0392782},
   {7.73126, 0.736956, 0.0218568},
   {8.00339, 0.589606, 0.0231838},
   {6.84775, 0.821282, 5.21952e-08}};
*/
/*
// Fit of JER for R=0.7, 8 TeV, 53X (20/fb, April 10)
  {{4.56, 0.942, 0.0295},  // y 0.0-0.5, chi2 20.8/33
   {5.07, 0.919, 0.0329},  // y 0.5-1.0, chi2 39.4/33
   {5.47, 0.961, 0.0372},  // y 1.0-1.5, chi2 29.8/33
   {5.81, 0.822, 0.0139},  // y 1.5-2.0, chi2 28.4/32
   {6.26, 0.666, 0.0210},  // y 2.0-2.5, chi2 16.2/25
   {5.33, 0.845, 0.0000}}; // y 2.5-3.0, chi2 11.6/19
*/
/*
// Parameter values from Sanmay (talk on Aug 6)
  {{5.507, 0.942, 0.029},
   {6.058, 0.893, 0.033},
   {5.695, 0.999, 0.033},
   {6.434, 0.810, 0.013},
   {7.474, 0.672, 0.016},
   {7.396, 0.557, 0.034}};
*/
/*
// Fit of JER for 8 TeV, 53X (Sanmay's JEC)
  {{4.35, 0.931, 0.0304},  // y 0.0-0.5, chi2 25.3/35
   {4.63, 0.927, 0.0325},  // y 0.5-1.0, chi2 44.5/35
   {5.35, 0.951, 0.0379},  // y 1.0-1.5, chi2 27.2/32
   {5.52, 0.816, 0.0151},  // y 1.5-2.0, chi2 26.5/28
   {5.94, 0.665, 0.0213},  // y 2.0-2.5, chi2 15.4/23
   {5.41, 0.825, 0.0000}}; // y 2.5-3.0, chi2 11.6/16
*/
/*
// Fit of JER for 8 TeV, 53X
  {{4.23, 0.940, 0.0299},  // y 0.0-0.5, chi2 22.4/35
   {4.62, 0.919, 0.0331},  // y 0.5-1.0, chi2 48.5/35
   {4.95, 0.965, 0.0374},  // y 1.0-1.5, chi2 34.8/32
   {5.69, 0.802, 0.0159},  // y 1.5-2.0, chi2 22.1/28
   {5.06, 0.765, 0.0069},  // y 2.0-2.5, chi2 28.0/23
   {4.76, 0.811, 0.0000}}; // y 2.5-3.0, chi2 32.5/16
*/
/*
  // Fit of JER for 8 TeV, 9/fb, 53X V5 (2012/12/11)
  {{5.18, 0.892, 0.0326},  // y 0.0-0.5, chi2 20.4/35
   {5.55, 0.867, 0.0360},  // y 0.5-1.0, chi2 23.0/34
   {6.43, 0.857, 0.0427},  // y 1.0-1.5, chi2 26.4/31
   {6.41, 0.742, 0.0228},  // y 1.5-2.0, chi2 13.4/27
   {6.49, 0.653, 0.0196},  // y 2.0-2.5, chi2 11.2/21
   {5.99, 0.810, -0.0000}}; // y 2.5-3.0, chi2 7.7/16
*/
/*
  // 7 TeV, 52X GRV23 (2012/08/29)
  {{3.51, 0.826, 0.0364},  // y 0.0-0.5, chi2 32.8/35
   {3.23, 0.851, 0.0367},  // y 0.5-1.0, chi2 23.9/34
   {4.36, 0.871, 0.0415},  // y 1.0-1.5, chi2 30.5/31
   {5.22, 0.713, 0.0229},  // y 1.5-2.0, chi2 16.9/27
   {5.07, 0.610, 0.0207},  // y 2.0-2.5, chi2 8.8/21
   {5.64, 0.735, 0.0000}}; // y 2.5-3.0, chi2 24.9/16
*/

// Hauke's resolutions
//const int _nres = 7;
//const int _npar = 4;
//float _pres_ak5calo[_nres][_npar] =
 // From CMSSW Spring10_PtResolution_AK5Calo
// {{/*0 0.5*/ 5.0913, 0.511619, 0, 0.325477},
//  {/*0.5 1*/ 4.93688, 0.543349, 0, 0.310638},
//  {/*1 1.5*/ 4.95729, 0.574202, 0, 0.32344},
//  {/*1.5 2*/ 3.37322, 0.96262, 0, 0.0860736},
//  {/*2 2.5*/ 4.06368, 0.351561, 0, 0.356167},
//  {/*2.5 3*/ 2.59292, 0.357679, 0, 0.321092},
//  {/*2.5 3*/ 2.59292, 0.357679, 0, 0.321092}}; // copy

// Hauke's numbers, Feb 11; confirmed Feb25 (in Fall10_PtResolution_AK5*.txt)
////float _pconst_ak5calo[_nres] = {0.033, 0.051, 0.039, 0.067, 0.045, 0.071,0.086};
// mk_resolution.C->resolution_datamc.pdf values on Jan 31st for time-dep JEC
//float _pscale_ak5pf[_nres] = {1.068, 1.061, 1.049, 1.168, 1.155, 1.193};//, 1.072}; // used for 7 TeV, official combination

double ptresolution(double pt, double eta) {

  //int ieta = min(_nres-1, int(fabs(eta) / 0.5));
  //int iy = min(5, ieta);
  int iy = min(5, int(fabs(eta) / 0.5 + 0.5));
  double res = 0;
  if (_ak7) {
    assert(!_ak4);
    res = sqrt(pow(vpar7[iy][0]/pt,2) + pow(vpar7[iy][1],2)/pt + 
	       pow(vpar7[iy][2],2));
  }
  else if (_ak4) {
    double kn = pow(0.4/0.5,2); //  scale noise term down by jet area
    res = sqrt(pow(kn*vpar5[iy][0]/pt,2) + pow(vpar5[iy][1],2)/pt + 
	       pow(vpar5[iy][2],2));
  }
  else
    res = sqrt(pow(vpar5[iy][0]/pt,2) + pow(vpar5[iy][1],2)/pt + 
	       pow(vpar5[iy][2],2));
  if (!_ismcjer) res *= kpar[iy][0];

  // replace with JME resolutions
  if (_usejme) {
    // Example code from:
    // https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_25/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h

    string resolutionFile = "../JRDatabase/textFiles/Autumn18_V1_MC/"
      "Autumn18_V1_MC_PtResolution_AK4PFchs.txt";
    string scaleFactorFile = "../JRDatabase/textFiles/Autumn18_V1_MC/"
      "Autumn18_V1_MC_SF_AK4PFchs.txt";
    
    //m_resolution_from_file.reset(new JME::JetREsolution(resolutionFile));
    //m_scale_factor_from_file.reset(new JME::JetResolutionScaleFactor(scaleFactorFile));

    JME::JetResolution *resolution = new JME::JetResolution(resolutionFile);
    JME::JetResolutionScaleFactor *resolution_sf = new JME::JetResolutionScaleFactor(scaleFactorFile);
  
    double rho = 18.;
    double jet_resolution = resolution->getResolution({{JME::Binning::JetPt, pt}, {JME::Binning::JetEta, eta}, {JME::Binning::Rho, rho}});
    double jer_sf = resolution_sf->getScaleFactor({{JME::Binning::JetEta, eta}}, Variation::NOMINAL);//m_systematic_variation);

    res = jet_resolution;
    if (!_ismcjer) res *= jer_sf;
  }

  // extra smearing from PU
  //double cpu2 = (_npv - _npv0) * pow(_spu/pt,2);
  //if (cpu2<-res*res) cpu2 = -0.75*res*res;
  //res = sqrt(res*res + cpu2);

  return res;
}

/*
double ptresolution_calo(double pt, double eta) {
  
  int ieta = min(_nres-1, int(fabs(eta) / 0.5));
  double N = _pres_ak5calo[ieta][0];
  double S = _pres_ak5calo[ieta][1];
  double C = _pres_ak5calo[ieta][2];
  double m = _pres_ak5calo[ieta][3];
  double c = _pconst_ak5calo[ieta];
  double k = _pscale_ak5pf[ieta]; // !!! Using PF scaling
  if (_ismcjer || !_jerkscale) k = 1.;
  if (_ismcjer ||  _jerkscale) c = 0.;

  double res = k*sqrt(TMath::Sign(N*N,N)/(pt*pt) + S*S*pow(pt,m-1) + C*C + c*c);

  return res;
}
*/

#endif // __ptresolution_h__
