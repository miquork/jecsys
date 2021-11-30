#ifndef globalFitSettings
#define globalFitSettings
#include <array>
#include <string>
#include <map>

using std::array;
using std::string;
using std::map;

// Input data
struct fitData {
  string name;
  string type;
  TGraphErrors *input;
  TGraphErrors *input2;
  TGraphErrors *output;
  TGraphErrors *output2;
};

// Systematics for input data
struct fitSyst {
  int idx;
  string name;
  string appliesTo;
  TH1D *hist;
};

// Shapes possible for JES and PF composition
struct fitShape {
  int idx;
  string name;
  string appliesTo;
  TF1 *func;
};

// Gaussian prior for fit parameters
bool _gf_penalizeFitPars = true;

// Listing of all available data sets
const unsigned int ndt = 2;
const array<array<string,2>,ndt> _gf_datasets = {{
    //{"mpfchs1_gamjet_a100","Rjet"}
    {"hdm_mpfchs1_gamjet","Rjet"},
    {"hdm_mpfchs1_zjet","Rjet"}
  }};

// Use only these data sets (empty for all)
const array<string,2> _gf_datasets_whitelist = {
  "hdm_mpfchs1_gamjet",
  "hdm_mpfchs1_zjet",
  //"mpfchs1_gamjet_a100",
};


// Listing of all available uncertainty sources
// {name, appliesTo, histogram name}
const int nsrc = 3;
const array<array<string,3>,nsrc> _gf_sources = {{
    //{"gamscale","mpfchs1_gamjet_a100","zee_gamjet_eig0"},
    {"FSRuncl","mpfchs1_gamjet_a100","hkfsr3_mpfchs1_gamjet_mpfu1"},
    //{"gamscale0","hdm_mpfchs1_gamjet","zee_gamjet_eig0"},
    //{"gamscale1","hdm_mpfchs1_gamjet","zee_gamjet_eig1"},
    //{"gamscale2","hdm_mpfchs1_gamjet","zee_gamjet_eig2"},
    {"FSRuncl","hdm_mpfchs1_gamjet","hkfsr3_mpfchs1_gamjet_mpfu1"},
    //{"FSRumc","hdm_mpfchs1_gamjet","hkfsr3_mpfchs1_gamjet_mpfu2"},
    //{"FSRnjet","hdm_mpfchs1_gamjet","hkfsr3_mpfchs1_gamjet_mpfn1"}
    {"FSRuncl","hdm_mpfchs1_zjet","hkfsr3_mpfchs1_zjet_mpfu1"},
  }};


// Listing of all available JEC fit shapes
// {name, appliesTo, function string}
const int nshp = 2;
const array<array<string,3>, nshp> _gf_shapes = {{
    {"gamscale","Rjet","0.02"},
    {"eescale","Rjet","0.01*log10(x/200.)"},
  }};

#endif
