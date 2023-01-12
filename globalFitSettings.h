#ifndef globalFitSettings
#define globalFitSettings
#include <array>
#include <string>
#include <map>

using std::array;
using std::string;
using std::map;

// Global minimum (stat.) uncertainty for input data
double globalErrMin = 0.002559; // chi2/NDF=101.0/101

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
  bool ispos;
  TF1 *func;
};

// Gaussian prior for fit parameters
bool _gf_penalizeFitPars = true;

// Permit baseline shift for JES with L2L3Res inputs. Needed for PFcomp fit
bool _gf_useJESref = true;//false;

// Listing of all available 'datasets'
// {name, type, name2}
// 'name' must match graph name in rootfiles/jecdata[X].root:ratio/eta_[Y]/
// 'type' decides which of 'shapes' below is applied, e.g. 'Rjet' or 'chf'
// 'name2' is secondary input, so far only used for 'multijet' modes
// How to use: comment out datasets that are not needed to switch them off
const unsigned int ndt = 44;
const array<array<string,3>,ndt> _gf_datasets = {{
    //{"mpfchs1_gamjet_a100", "Rjet",""},
    //{"mpfchs1_zjet_a100", "Rjet",""},
    //{"hdm_mpfchs1_gamjet",  "Rjet",""},
    {"hdm_mpfchs1_multijet","Rjet","crecoil_multijet_a30"},
    {"hdm_gamjet",  "Rjet",""},
    {"hdm_mpfchs1_zjet",    "Rjet",""},
    {"hdm_mpfchs1_zlljet",    "Rjet",""},
    //{"mpfchs1_hadw_a30",    "Rjet",""},
    //{"hadw_a30",    "Rjet",""},
    {"hdm_hadw",    "Rjet",""},
    {"hdm_incjet",  "Rjet",""},
    //
    {"chf_pfjet_a30",  "chf",""},
    {"nhf_pfjet_a30",  "nhf",""},
    {"nef_pfjet_a30",  "nef",""},
    {"cef_pfjet_a30",  "cef",""},
    {"muf_pfjet_a30",  "muf",""},
    {"chf_gamjet_a100",  "chf",""},
    {"nhf_gamjet_a100",  "nhf",""},
    {"nef_gamjet_a100",  "nef",""},
    {"cef_gamjet_a100",  "cef",""},
    {"muf_gamjet_a100",  "muf",""},
    {"chf_zjet_a100",  "chf",""},
    {"nhf_zjet_a100",  "nhf",""},
    {"nef_zjet_a100",  "nef",""},
    {"cef_zjet_a100",  "cef",""},
    {"muf_zjet_a100",  "muf",""},
    {"chf_zlljet_a100",  "chf",""},
    {"nhf_zlljet_a100",  "nhf",""},
    {"nef_zlljet_a100",  "nef",""},
    {"cef_zlljet_a100",  "cef",""},
    {"muf_zlljet_a100",  "muf",""},
    //{"cef_zmmjet_a100",  "cef",""},
    //{"muf_zmmjet_a100",  "muf",""},
    //
    {"chf_zjet_ren",  "chf",""},
    {"nhf_zjet_ren",  "nhf",""},
    {"nef_zjet_ren",  "nef",""},
    {"chf_zlljet_ren",  "chf",""},
    {"nhf_zlljet_ren",  "nhf",""},
    {"nef_zlljet_ren",  "nef",""},
    {"chf_pfjet_ren",  "chf",""},
    {"nhf_pfjet_ren",  "nhf",""},
    {"nef_pfjet_ren",  "nef",""},
    {"chf_gamjet_ren",  "chf",""},
    {"nhf_gamjet_ren",  "nhf",""},
    {"nef_gamjet_ren",  "nef",""},
    //
    {"chf_cmb_ren",  "chf",""},
    {"nhf_cmb_ren",  "nhf",""},
    {"nef_cmb_ren",  "nef",""},
    //
    {"hdm_cmb",  "Rjet",""},
    {"hdm_cmb_mj",  "Rjet",""},
  }};

// Use only these data sets (empty for all)
// {name}
// 'name' should match the dataset name in the list above
// How to use: uncomment individual datasets to use only those
const array<string,29> _gf_datasets_whitelist = {
  //"mpfchs1_gamjet_a100",
  //"mpfchs1_zjet_a100",
  //"hdm_mpfchs1_gamjet",
  /*
  "hdm_gamjet",
  "hdm_mpfchs1_zjet",
  "hdm_mpfchs1_zlljet",
  //"mpfchs1_hadw_a30",
  "hdm_hadw",
  "hdm_mpfchs1_multijet",
  */
  //"hdm_mpfchs1_zlljet",
  /*
  "chf_gamjet_a100",
  "nhf_gamjet_a100",
  "nef_gamjet_a100",
  "cef_gamjet_a100",
  "muf_gamjet_a100",
  "chf_pfjet_a30",
  "nhf_pfjet_a30",
  "nef_pfjet_a30",
  "cef_pfjet_a30",
  "muf_pfjet_a30",
  "chf_zjet_a100",
  "nhf_zjet_a100",
  "nef_zjet_a100",
  "cef_zjet_a100",
  "muf_zjet_a100",
  "chf_zlljet_a100",
  "nhf_zlljet_a100",
  "nef_zlljet_a100",
  "cef_zlljet_a100",
  "muf_zlljet_a100",
  */
  /*
  "chf_pfjet_ren",
  "nhf_pfjet_ren",
  "nef_pfjet_ren",
  "cef_pfjet_a30",
  "muf_pfjet_a30",
  "chf_gamjet_ren",
  "nhf_gamjet_ren",
  "nef_gamjet_ren",
  "cef_gamjet_a100",
  "muf_gamjet_a100",

  "chf_zjet_ren",
  "nhf_zjet_ren",
  "nef_zjet_ren",
  "cef_zjet_a100",
  "muf_zjet_a100",

  "chf_zlljet_ren",
  "nhf_zlljet_ren",
  "nef_zlljet_ren",
  "cef_zlljet_a100",
  "muf_zlljet_a100",
  */
  "chf_cmb_ren",
  "nhf_cmb_ren",
  "nef_cmb_ren",
  //"cef_pfjet_a100",
  //"muf_pfjet_a100",

  //"hdm_incjet",
  //"hdm_cmb", // 2017H
  "hdm_cmb_mj", // Not 2017H

  //"cef_zmmjet_a100",
  //"muf_zmmjet_a100",
};

// Listing one-sided (positive-definite) sources
const int npos = 17;
const array<string,npos> _gf_posdef =
  {"tv4","tv402","tv404","tv405","tv410","tv416","tv420","tv430",
   "tv3n1","tv300pn",
   //"hhp3",
   "hhred103","hhred100","hhred097", "hhblue103","hhblue100","hhblue097"};

// Listing source limits
// (How to best implement this?)

// Listing of all available uncertainty 'sources'
// {name, appliesTo, histname}
// 'name' can repeat, e.g. 'FSRuncl' for multiple inputs, and uses same nuisance
// 'appliesTo' must match the dataset name above
// 'histname' must match histogram name in rootfiles/jecdata[X].root:[Y]/sys/
// How to use: list all 'sources'. Ones not matching active inputs are ignored
// NB: For one-sided (positive) sources, add them to _gf_posdef list
const int nsrc = 13;
const array<array<string,3>,nsrc> _gf_sources = {{
    /*
    {"gamscale0","hdm_gamjet","zee_gamjet_eig0"},
    {"gamscale1","hdm_gamjet","zee_gamjet_eig1"},
    {"gamscale2","hdm_gamjet","zee_gamjet_eig2"},

    {"FSR","hdm_hadw","hadw_ptave_fsr"},
    */
    /*
    {"FSRunclDT","hdm_gamjet","hkfsr3_mpfchs1_gamjet_mpfu1"},
    {"FSRunclDT","hdm_mpfchs1_zjet","hkfsr3_mpfchs1_zjet_mpfu1"},
    {"FSRunclDT","hdm_mpfchs1_zlljet","hkfsr3_mpfchs1_zlljet_mpfu1"},
    //{"FSRunclDT","hdm_mpfchs1_multijet","hkfsr3_mpfchs1_multijet_mpfu1"},

    {"FSRunclMC","hdm_gamjet","hkfsr3_mpfchs1_gamjet_mpfu2"},
    {"FSRunclMC","hdm_mpfchs1_zjet","hkfsr3_mpfchs1_zjet_mpfu2"},
    {"FSRunclMC","hdm_mpfchs1_zlljet","hkfsr3_mpfchs1_zlljet_mpfu2"},
    //{"FSRunclMC","hdm_mpfchs1_multijet","hkfsr3_mpfchs1_multijet_mpfu2"},
    */

    //{"FSRunclDTMJ","hdm_mpfchs1_multijet","hkfsr3_mpfchs1_multijet_mpfu1"},
    //{"FSRunclMCMJ","hdm_mpfchs1_multijet","hkfsr3_mpfchs1_multijet_mpfu2"},
  }};


// Listing of all available JEC fit 'shapes'
// {name, appliesTo, funcstring}
// 'name' can repeat e.g. for 'Rjet' and 'chf', and uses same fit parameter
// 'appliesTo' must match one of the 'type' in dataset listing
// 'funcstring' is a TF1-style string for a fixed functions. Use 'x' for pT
// How this works: each unique 'name' is assigned a fit parameter that
// multiplies 'funcstring' when applied to dataset with same 'appliesTo' type
//
// Fits are produced with minitools/fullSimShapes.C with doProduction=true
// Copy them from pdf/fullSimShapes/txt and /prod to /handpicked
// Each plot/txt pair is identified by unique MD5 checksum calculated from
// concatenated Rjet,chf,nef,nhf functions strings
const int nshp = 112;
const array<array<string,3>, nshp> _gf_shapes = {{
    //{"p0","Rjet","0.0075"},
    //{"p1","Rjet","1/x"},
    //{"p2","Rjet","0.01*x/1000."},
    //{"p0","nef","0.0050"},
    //{"p3","Rjet","0.01*pow(x/1000.,2)"},
    //{"hadscaledummy","Rjet","0.01*log10(x/200.)"},    
    //{"hl1diff","Rjet","0+100.*(0.130492/x+0.00261293*log(x)/x)"},
    //{"hl1diff","chf","0"},
    //{"hl1diff","nhf","0"},
    //{"hl1diff","nef","0"},
    //
    //{"tm3","Rjet","0.4804-(5+0.003886)*pow(x,-(0.3+-0.008537))"},
    //{"tm3","chf","-(6+-0.1424)+(7+-0.6886)*pow(x,-(0.1+-0.02696))+(6e-4+0.0001251)*pow(x,1+-0.04998)"},
    //{"tm3","nef","(3+4.378)-(5+3.762)*pow(x,-(0.1+-0.03012))-(1e-2+0.001722)*pow(x,1+-0.3309)"},
    //{"tm3","nhf","(1+-0.7921)-(1+-1.827)*pow(x,-(0.1+0.2451))-(1e-4+-0.01532)*pow(x,1+-0.5104)"},
    // "tm3" checksum: e22a45499ef35f55ae41acd5d5155392
    //{"tv2c","Rjet","2.953-(1+7.396)*pow(max(60.,x),-(0.3+-0.2205))+12.02/log(max(60.,x))"},
    //{"tv2c","chf","1.822-(0.2+0.3817)*pow(max(60.,x),0.3+0.056)+(2e-2+0.01495)*pow(max(60.,x),1+-0.3286)"},
    //{"tv2c","nef","-2.446+(0.2+0.9376)*pow(max(60.,x),0.3+-0.08704)-(1e-2+-8.191e-05)*pow(max(60.,x),1+-0.2494)"},
    //{"tv2c","nhf","-0.794+(0.1+0.296)*pow(max(60.,x),0.3+-0.1164)"},
    // "tv2c" checksum: e5e2f58673dc24541ecd61b63a3c0206
    //{"tv3b","Rjet","3.826-(1.0+5.546)*pow(x,-(0.3+-0.2528))+4.993/log(x)"},
    //{"tv3b","chf","-0.1952+log(max(60.,x))*(0.1042+log(max(60.,x))*(-0.00815+log(max(60.,x))*(-0.004963+log(max(60.,x))*(-0.0008004+log(max(60.,x))*(-3.058e-05+log(max(60.,x))*2.125e-05)))))"},
    //{"tv3b","nef","12.92+log(max(60.,x))*(-5.953+log(max(60.,x))*(0.432+log(max(60.,x))*(0.1087+log(max(60.,x))*(-0.00318+log(max(60.,x))*(-0.002581+log(max(60.,x))*0.0001877)))))"},
    //{"tv3b","nhf","-13.69+log(max(60.,x))*(6.225+log(max(60.,x))*(-0.4248+log(max(60.,x))*(-0.1153+log(max(60.,x))*(0.003793+log(max(60.,x))*(0.002942+log(max(60.,x))*-0.0002327)))))"},
    // "tv3b" checksum: dce2a5d0406417661c36474156d2955e

    //{"hm3","Rjet","-2.031+log(max(30.,x))*(1.086+log(max(30.,x))*(-0.1244+log(max(30.,x))*(-0.02109+log(max(30.,x))*(0.001102+log(max(30.,x))*(0.0005024+log(max(30.,x))*-3.976e-05)))))"},
    //{"hm3","chf","2.123+log(max(45.,x))*(-1.032+log(max(45.,x))*(0.09046+log(max(45.,x))*(0.02061+log(max(45.,x))*(-0.0006251+log(max(45.,x))*(-0.0004836+log(max(45.,x))*3.337e-05)))))"},
    ///{"hm3","nef","-11.17+log(max(45.,x))*(5.404+log(max(45.,x))*(-0.4343+log(max(45.,x))*(-0.09588+log(max(45.,x))*(0.003923+log(max(45.,x))*(0.002342+log(max(45.,x))*-0.0001815)))))"},
    //{"hm3","nhf","8.609+log(max(45.,x))*(-4.098+log(max(45.,x))*(0.2895+log(max(45.,x))*(0.07626+log(max(45.,x))*(-0.002491+log(max(45.,x))*(-0.001932+log(max(45.,x))*0.0001493)))))"},
    // "hm3" checksum: 34a9d830647004bce3c067d3051e49a4

    //{"hg3","Rjet","-9.593+log(x)*(5.365+log(x)*(-0.5851+log(x)*(-0.0916+log(x)*(0.006879+log(x)*(0.00247+log(x)*-0.0002162)))))"},
    //{"hg3","chf","4.932+log(x)*(-2.826+log(x)*(0.3241+log(x)*(0.04885+log(x)*(-0.004092+log(x)*(-0.001363+log(x)*0.0001255)))))"},
    //{"hg3","nef","-1.993+log(x)*(1.246+log(x)*(-0.1619+log(x)*(-0.02078+log(x)*(0.002244+log(x)*(0.0006096+log(x)*-6.341e-05)))))"},
    //{"hg3","nhf","-2.119+log(x)*(1.13+log(x)*(-0.1092+log(x)*(-0.02241+log(x)*(0.001192+log(x)*(0.0006236+log(x)*-4.888e-05)))))"},
    // "hg3" checksum: 4f67c3d26422ded6ad4913a4acccc671
    //{"hr3","Rjet","-11.47+log(max(60.,x))*(5.175+log(max(60.,x))*(-0.3588+log(max(60.,x))*(-0.0907+log(max(60.,x))*(0.002202+log(max(60.,x))*(0.002044+log(max(60.,x))*-0.0001393)))))"},
    //{"hr3","chf","9.162+log(max(60.,x))*(-4.171+log(max(60.,x))*(0.2949+log(max(60.,x))*(0.07514+log(max(60.,x))*(-0.002166+log(max(60.,x))*(-0.001776+log(max(60.,x))*0.0001296)))))"},
    //{"hr3","nef","-11.49+log(max(60.,x))*(5.444+log(max(60.,x))*(-0.4254+log(max(60.,x))*(-0.09711+log(max(60.,x))*(0.003954+log(max(60.,x))*(0.002406+log(max(60.,x))*-0.0001913)))))"},
    //{"hr3","nhf","3.501+log(max(60.,x))*(-1.777+log(max(60.,x))*(0.1526+log(max(60.,x))*(0.03359+log(max(60.,x))*(-0.001877+log(max(60.,x))*(-0.0009396+log(max(60.,x))*8.466e-05)))))"},
    // "hr3" checksum: 52dde108aadd6f6212c7788ea94eef14
    //{"hb1","Rjet","-5.919+log(max(60.,x))*(3.485+log(max(60.,x))*(-0.3927+log(max(60.,x))*(-0.06647+log(max(60.,x))*(0.004971+log(max(60.,x))*(0.001908+log(max(60.,x))*-0.0001726)))))"},
    //{"hb1","chf","1.881+log(max(60.,x))*(-1.296+log(max(60.,x))*(0.1736+log(max(60.,x))*(0.02616+log(max(60.,x))*(-0.002446+log(max(60.,x))*(-0.0008044+log(max(60.,x))*7.711e-05)))))"},
    //{"hb1","nef","-1.578+log(max(60.,x))*(0.9897+log(max(60.,x))*(-0.1209+log(max(60.,x))*(-0.01783+log(max(60.,x))*(0.001579+log(max(60.,x))*(0.0004951+log(max(60.,x))*-4.711e-05)))))"},
    //{"hb1","nhf","0.9254+log(max(60.,x))*(-0.2228+log(max(60.,x))*(-0.02179+log(max(60.,x))*(0.001546+log(max(60.,x))*(0.0006709+log(max(60.,x))*(6.646e-05+log(max(60.,x))*-1.192e-05)))))"},
    // "hb1" checksum: 4b5bbd4c5e7ad7566afede7502149be1

    //{"em3","Rjet","5.4+log(max(60.,x))*(-2.137+log(max(60.,x))*(0.0747+log(max(60.,x))*(0.02768+log(max(60.,x))*(-7.826e-05+log(max(60.,x))*(-0.0004993+log(max(60.,x))*3.088e-05)))))"},
    //{"em3","chf","1.233+log(max(60.,x))*(-0.6129+log(max(60.,x))*(0.07748+log(max(60.,x))*(0.01447+log(max(60.,x))*(-0.0008236+log(max(60.,x))*(-0.0004073+log(max(60.,x))*3.229e-05)))))"},
    //{"em3","nef","0.6707+log(max(60.,x))*(-0.2559+log(max(60.,x))*(-0.02444+log(max(60.,x))*(0.0009631+log(max(60.,x))*(0.0005403+log(max(60.,x))*(5.475e-05+log(max(60.,x))*-8.592e-06)))))"},
    //{"em3","nhf","0.00299+log(max(60.,x))*(0.0287+log(max(60.,x))*(0.001202+log(max(60.,x))*(-0.0005139+log(max(60.,x))*(-0.0001163+log(max(60.,x))*(-6.21e-06+log(max(60.,x))*3.392e-06)))))"},
    // "em3" checksum: aacef53137acaa484ece7c34b7731dad  (old v11)
    /*
    {"em3","Rjet","-0.7339+log(max(60.,x))*(-0.01379+log(max(60.,x))*(0.0005509+log(max(60.,x))*(-0.0006042+log(max(60.,x))*(-0.000185+log(max(60.,x))*(-1.734e-05+log(max(60.,x))*4.528e-06)))))"},
    {"em3","chf","0.6556+log(max(60.,x))*(-0.1551+log(max(60.,x))*(0.0182+log(max(60.,x))*(0.003993+log(max(60.,x))*(-0.0002847+log(max(60.,x))*(-0.0001214+log(max(60.,x))*1.098e-05)))))"},
    {"em3","nef","0.1762+log(max(60.,x))*(-0.242+log(max(60.,x))*(0.0004845+log(max(60.,x))*(0.004878+log(max(60.,x))*(0.0002149+log(max(60.,x))*(-0.000119+log(max(60.,x))*6.146e-06)))))"},
    {"em3","nhf","-0.5815+log(max(60.,x))*(0.2909+log(max(60.,x))*(-0.01477+log(max(60.,x))*(-0.00633+log(max(60.,x))*(7.409e-05+log(max(60.,x))*(0.0001685+log(max(60.,x))*-1.181e-05)))))"},
    // "em3" checksum: b00d510a6c444e71353f5b2bd28d99d2 (new v15)
    */
    {"em3","Rjet","-0.3222+log(x)*(0.1973+log(x)*(-0.0774+log(x)*(-0.006124+log(x)*(0.001035+log(x)*(0.0002058+log(x)*-2.089e-05)))))"},
    {"em3","chf","2.59+log(x)*(-1.188+log(x)*(0.1324+log(x)*(0.01962+log(x)*(-0.001522+log(x)*(-0.0004592+log(x)*3.848e-05)))))"},
    {"em3","nef","-2.715+log(x)*(1.225+log(x)*(-0.1397+log(x)*(-0.02023+log(x)*(0.001699+log(x)*(0.0004898+log(x)*-4.382e-05)))))"},
    {"em3","nhf","2.168+log(x)*(-1.761+log(x)*(0.4903+log(x)*(-0.0334+log(x)*(-0.007356+log(x)*(0.001323+log(x)*-5.843e-05)))))"},
    // "em3" checksum: c8cb4a54ecee28c98c7fe286c0beafa9

    /*
    {"pm3","Rjet","-0.5754+log(x)*(-0.08557+log(x)*0.007061)"},
    {"pm3","chf","0.05363+log(x)*(0.1341+log(x)*-0.01254)"},
    {"pm3","nef","-0.2875+log(x)*(-0.06739+log(x)*0.005384)"},
    {"pm3","nhf","0.2565+log(x)*(-0.07617+log(x)*0.008004)"},
    // "pm3" checksum: 2c3e797437198f24f27b96b299ca072c (toyPF only!)
    */
    /*
    {"tv4","Rjet","-1.744+log(max(20.,x))*(-0.1024+log(max(20.,x))*(0.05999+log(max(20.,x))*(0.01157+log(max(20.,x))*(-0.00277+log(max(20.,x))*(8.454e-05+log(max(20.,x))*5.22e-06)))))"},
    {"tv4","chf","-1.45+log(max(20.,x))*(0.06786+log(max(20.,x))*(0.08777+log(max(20.,x))*(-0.001014+log(max(20.,x))*(-0.01193+log(max(20.,x))*(0.002165+log(max(20.,x))*-0.000104)))))"},
    {"tv4","nef","2.276+log(max(20.,x))*(-0.3972+log(max(20.,x))*(-0.226+log(max(20.,x))*(0.01108+log(max(20.,x))*(0.02858+log(max(20.,x))*(-0.005714+log(max(20.,x))*0.0003052)))))"},
    {"tv4","nhf","-0.8347+log(max(20.,x))*(0.3361+log(max(20.,x))*(0.1384+log(max(20.,x))*(-0.0109+log(max(20.,x))*(-0.01658+log(max(20.,x))*(0.003571+log(max(20.,x))*-0.0002039)))))"},
    // "tv4" checksum: e7fa22497bb42f656cbdeaeef502c1ae
    */
    //{"tv402","Rjet","1.909+log(max(20.,x))*(-0.3779+log(max(20.,x))*(-0.1912+log(max(20.,x))*(0.01457+log(max(20.,x))*(0.01258+log(max(20.,x))*(-0.002258+log(max(20.,x))*0.000108)))))"},//old
    //{"tv402","chf","0.5855+log(max(20.,x))*(0.05936+log(max(20.,x))*(-0.02901+log(max(20.,x))*(-0.01958+log(max(20.,x))*(-0.003297+log(max(20.,x))*(0.00149+log(max(20.,x))*-9.738e-05)))))"},//old
    //{"tv402","nef","1.47+log(max(20.,x))*(-0.4284+log(max(20.,x))*(-0.1878+log(max(20.,x))*(0.01976+log(max(20.,x))*(0.02608+log(max(20.,x))*(-0.005632+log(max(20.,x))*0.0003119)))))"},//old
    //{"tv402","nhf","-1.98+log(max(20.,x))*(0.347+log(max(20.,x))*(0.2087+log(max(20.,x))*(0.0005292+log(max(20.,x))*(-0.02154+log(max(20.,x))*(0.003855+log(max(20.,x))*-0.0001971)))))"},//old
    // "tv402" checksum: d56e060d3b86fe97caa42cc49235db86 //old
    /*
    {"tv402","Rjet","3.958+log(max(20.,x))*(-1.89+log(max(20.,x))*(0.1272+log(max(20.,x))*(0.02788+log(max(20.,x))*(-0.0007582+log(max(20.,x))*(-0.0005579+log(max(20.,x))*3.858e-05)))))"},
    {"tv402","chf","-0.1554+log(max(20.,x))*(0.6983+log(max(20.,x))*(-0.1954+log(max(20.,x))*(-0.01841+log(max(20.,x))*(0.002612+log(max(20.,x))*(0.0006059+log(max(20.,x))*-5.758e-05)))))"},
    {"tv402","nef","-11.5+log(max(20.,x))*(10.24+log(max(20.,x))*(-3.201+log(max(20.,x))*(0.3069+log(max(20.,x))*(0.04157+log(max(20.,x))*(-0.01009+log(max(20.,x))*0.0005254)))))"},
    {"tv402","nhf","11.4+log(max(20.,x))*(-10.2+log(max(20.,x))*(2.98+log(max(20.,x))*(-0.2092+log(max(20.,x))*(-0.04651+log(max(20.,x))*(0.008715+log(max(20.,x))*-0.0004105)))))"},
    // "tv402" checksum: 46dd5517ed3970927e2cf61d0debf267
    */
    {"tv402","Rjet","-3.512+log(x)*(3.975+log(x)*(-1.433+log(x)*(0.1443+log(x)*(0.01588+log(x)*(-0.003781+log(x)*0.0001845)))))"},
    {"tv402","chf","-14.72+log(x)*(12.93+log(x)*(-3.671+log(x)*(0.2582+log(x)*(0.04383+log(x)*(-0.007821+log(x)*0.0003428)))))"},
    {"tv402","nef","2.082+log(x)*(-0.8289+log(x)*(-0.2106+log(x)*(0.1027+log(x)*(0.00176+log(x)*(-0.002916+log(x)*0.0002057)))))"},
    {"tv402","nhf","14.14+log(x)*(-13.34+log(x)*(4.169+log(x)*(-0.3458+log(x)*(-0.06092+log(x)*(0.01287+log(x)*-0.0006435)))))"},
    // "tv402" checksum: 6b76d4d2922fa185b8fe54ae5a5a56a5

    /*
    {"tv404","Rjet","4.181+log(max(40.,x))*(-0.892+log(max(40.,x))*(-0.2645+log(max(40.,x))*(0.02442+log(max(40.,x))*(0.01801+log(max(40.,x))*(-0.003297+log(max(40.,x))*0.0001588)))))"},
    {"tv404","chf","1.446+log(max(40.,x))*(0.05596+log(max(40.,x))*(-0.06008+log(max(40.,x))*(-0.02352+log(max(40.,x))*(-0.002307+log(max(40.,x))*(0.001514+log(max(40.,x))*-0.0001039)))))"},
    {"tv404","nef","2.543+log(max(40.,x))*(-0.9361+log(max(40.,x))*(-0.1988+log(max(40.,x))*(0.04779+log(max(40.,x))*(0.02337+log(max(40.,x))*(-0.005779+log(max(40.,x))*0.000332)))))"},
    {"tv404","nhf","-3.898+log(max(40.,x))*(0.854+log(max(40.,x))*(0.2537+log(max(40.,x))*(-0.02336+log(max(40.,x))*(-0.02055+log(max(40.,x))*(0.004147+log(max(40.,x))*-0.0002217)))))"},
    // "tv404" checksum: 529835cc15ef6137843f7dade538cbe7

    //{"tv405","Rjet","5.259+log(max(55.,x))*(-1.148+log(max(55.,x))*(-0.2753+log(max(55.,x))*(0.0303+log(max(55.,x))*(0.01758+log(max(55.,x))*(-0.003261+log(max(55.,x))*0.0001557)))))"},
    //{"tv405","chf","2.461+log(max(55.,x))*(0.02618+log(max(55.,x))*(-0.09818+log(max(55.,x))*(-0.0296+log(max(55.,x))*(-0.001592+log(max(55.,x))*(0.001828+log(max(55.,x))*-0.0001361)))))"},
    //{"tv405","nef","-1.845+log(max(55.,x))*(-0.01624+log(max(55.,x))*(0.07542+log(max(55.,x))*(0.02256+log(max(55.,x))*(0.001107+log(max(55.,x))*(-0.001437+log(max(55.,x))*0.000108)))))"},
    //{"tv405","nhf","-7.361+log(max(55.,x))*(1.759+log(max(55.,x))*(0.3868+log(max(55.,x))*(-0.05384+log(max(55.,x))*(-0.02855+log(max(55.,x))*(0.00607+log(max(55.,x))*-0.0003244)))))"},
    // "tv405" checksum: 8fce5c498456fe02c0363856463570fc

    {"tv410","Rjet","1.008+log(max(80.,x))*(-0.02373+log(max(80.,x))*(-0.0399+log(max(80.,x))*(-0.008647+log(max(80.,x))*(0.0001365+log(max(80.,x))*(0.0005318+log(max(80.,x))*-4.521e-05)))))"},
    {"tv410","chf","5.1+log(max(130.,x))*(-0.1557+log(max(130.,x))*(-0.16+log(max(130.,x))*(-0.02868+log(max(130.,x))*(0.0005512+log(max(130.,x))*(0.001442+log(max(130.,x))*-0.0001122)))))"},
    {"tv410","nef","-4.621+log(max(130.,x))*(0.1698+log(max(130.,x))*(0.152+log(max(130.,x))*(0.02635+log(max(130.,x))*(-0.0008613+log(max(130.,x))*(-0.001435+log(max(130.,x))*0.0001173)))))"},
    {"tv410","nhf","-0.7254+log(max(130.,x))*(0.003569+log(max(130.,x))*(0.01724+log(max(130.,x))*(0.003567+log(max(130.,x))*(0.0001682+log(max(130.,x))*(-9.359e-05+log(max(130.,x))*3.722e-06)))))"},
    // "tv410" checksum: fe44ab7f7e70aeeeccaccfdb009e14b5

    {"tv416","Rjet","2.377+log(max(210.,x))*(-0.157+log(max(210.,x))*(-0.07272+log(max(210.,x))*(-0.008669+log(max(210.,x))*(0.0009315+log(max(210.,x))*(0.0005091+log(max(210.,x))*-4.791e-05)))))"},
    {"tv416","chf","5.568+log(max(210.,x))*(-0.2539+log(max(210.,x))*(-0.1519+log(max(210.,x))*(-0.02099+log(max(210.,x))*(0.001111+log(max(210.,x))*(0.0009707+log(max(210.,x))*-7.739e-05)))))"},
    {"tv416","nef","-5.902+log(max(210.,x))*(0.3252+log(max(210.,x))*(0.1717+log(max(210.,x))*(0.02233+log(max(210.,x))*(-0.001709+log(max(210.,x))*(-0.001169+log(max(210.,x))*0.0001005)))))"},
    {"tv416","nhf","-0.2774+log(max(210.,x))*(-0.02241+log(max(210.,x))*(0.001053+log(max(210.,x))*(0.0009887+log(max(210.,x))*(0.0002543+log(max(210.,x))*(3.266e-05+log(max(210.,x))*-5.764e-06)))))"},
    // "tv416" checksum: d0312fcdecc345f5cb895597dcb31e15

    {"tv420","Rjet","2.257+log(max(260.,x))*(-0.1588+log(max(260.,x))*(-0.06471+log(max(260.,x))*(-0.006918+log(max(260.,x))*(0.0008498+log(max(260.,x))*(0.0003971+log(max(260.,x))*-3.745e-05)))))"},
    {"tv420","chf","3.781+log(max(260.,x))*(-0.153+log(max(260.,x))*(-0.09059+log(max(260.,x))*(-0.01219+log(max(260.,x))*(0.0004525+log(max(260.,x))*(0.000463+log(max(260.,x))*-3.137e-05)))))"},
    {"tv420","nef","-4.801+log(max(260.,x))*(0.2746+log(max(260.,x))*(0.1287+log(max(260.,x))*(0.01533+log(max(260.,x))*(-0.001261+log(max(260.,x))*(-0.0007476+log(max(260.,x))*6.217e-05)))))"},
    {"tv420","nhf","-0.2152+log(max(260.,x))*(-0.01666+log(max(260.,x))*(0.0006959+log(max(260.,x))*(0.0006495+log(max(260.,x))*(0.0001601+log(max(260.,x))*(2.012e-05+log(max(260.,x))*-2.888e-06)))))"},
    // "tv420" checksum: 558d56f719bfc29cfc8312ee172b2df5

    {"tv430","Rjet","1.429+log(max(400.,x))*(-0.1116+log(max(400.,x))*(-0.03549+log(max(400.,x))*(-0.002987+log(max(400.,x))*(0.000466+log(max(400.,x))*(0.0001632+log(max(400.,x))*-1.524e-05)))))"},
    {"tv430","chf","0.8448+log(max(400.,x))*(0.0439+log(max(400.,x))*(-0.006351+log(max(400.,x))*(-0.002596+log(max(400.,x))*(-0.0004886+log(max(400.,x))*(-3.931e-05+log(max(400.,x))*1.323e-05)))))"},
    {"tv430","nef","-0.7171+log(max(400.,x))*(-0.03446+log(max(400.,x))*(0.006104+log(max(400.,x))*(0.002317+log(max(400.,x))*(0.0004232+log(max(400.,x))*(3.117e-05+log(max(400.,x))*-1.254e-05)))))"},
    {"tv430","nhf","-0.1236+log(max(400.,x))*(-0.008847+log(max(400.,x))*(0.0002758+log(max(400.,x))*(0.0002674+log(max(400.,x))*(6.114e-05+log(max(400.,x))*(7.497e-06+log(max(400.,x))*-5.947e-07)))))"},
    // "tv430" checksum: 6530c2d9c8da33e62d283e3ab8005c9b
    */
    /*
    {"tv3n1","Rjet","-3.701+log(max(20.,x))*(0.2984+log(max(20.,x))*(0.2522+log(max(20.,x))*(-0.004803+log(max(20.,x))*(-0.01514+log(max(20.,x))*(0.002354+log(max(20.,x))*-0.000105)))))"},
    {"tv3n1","chf","-2.666+log(max(20.,x))*(0.162+log(max(20.,x))*(0.1868+log(max(20.,x))*(0.01373+log(max(20.,x))*(-0.01726+log(max(20.,x))*(0.002511+log(max(20.,x))*-0.0001114)))))"},
    {"tv3n1","nef","1.055+log(max(20.,x))*(-0.02634+log(max(20.,x))*(-0.06437+log(max(20.,x))*(-0.007024+log(max(20.,x))*(0.00537+log(max(20.,x))*(-0.0006669+log(max(20.,x))*2.568e-05)))))"},
    {"tv3n1","nhf","1.798+log(max(20.,x))*(-0.1881+log(max(20.,x))*(-0.1403+log(max(20.,x))*(-0.00308+log(max(20.,x))*(0.01302+log(max(20.,x))*(-0.002139+log(max(20.,x))*0.0001031)))))"},
    // "tv3n1" checksum: 7ca4e89f1263f751c45b03474cf1f81d
    */
    {"tv3n1","Rjet","17.03+log(x)*(-14.94+log(x)*(4.142+log(x)*(-0.3059+log(x)*(-0.047+log(x)*(0.008868+log(x)*-0.0003968)))))"},
    {"tv3n1","chf","17.76+log(x)*(-16.17+log(x)*(4.619+log(x)*(-0.3464+log(x)*(-0.05633+log(x)*(0.01074+log(x)*-0.0004876)))))"},
    {"tv3n1","nef","-8.682+log(x)*(7.874+log(x)*(-2.232+log(x)*(0.164+log(x)*(0.02771+log(x)*(-0.00521+log(x)*0.0002358)))))"},
    {"tv3n1","nhf","-9.61+log(x)*(8.698+log(x)*(-2.455+log(x)*(0.1624+log(x)*(0.03711+log(x)*(-0.006591+log(x)*0.0002973)))))"},
    // "tv3n1" checksum: 348f764eaddd02005d2c66cabf4fc863

    /*
    {"tv300pn","Rjet","0.3318+log(max(60.,x))*(0.0003536+log(max(60.,x))*(-0.0135+log(max(60.,x))*(-0.003663+log(max(60.,x))*(-8.071e-05+log(max(60.,x))*(0.0002107+log(max(60.,x))*-1.596e-05)))))"},
    {"tv300pn","chf","-3.905+log(max(60.,x))*(1.074+log(max(60.,x))*(0.2158+log(max(60.,x))*(-0.03556+log(max(60.,x))*(-0.01723+log(max(60.,x))*(0.003427+log(max(60.,x))*-0.0001636)))))"},
    {"tv300pn","nef","6.304+log(max(60.,x))*(-1.732+log(max(60.,x))*(-0.3375+log(max(60.,x))*(0.06036+log(max(60.,x))*(0.02726+log(max(60.,x))*(-0.005913+log(max(60.,x))*0.0003072)))))"},
    {"tv300pn","nhf","0.08994+log(max(60.,x))*(-0.00834+log(max(60.,x))*(-0.006685+log(max(60.,x))*(-0.00142+log(max(60.,x))*(0.0001282+log(max(60.,x))*(0.0001672+log(max(60.,x))*-1.586e-05)))))"},
    // "tv300pn" checksum: ec49e9738a835f17217d3eb1bad31a12
    */
    {"tv300pn","Rjet","-0.8581+log(max(30.,x))*(0.5263+log(max(30.,x))*(-0.06701+log(max(30.,x))*(-0.008783+log(max(30.,x))*(0.0008284+log(max(30.,x))*(0.0002308+log(max(30.,x))*-2.113e-05)))))"},
    {"tv300pn","chf","6.831+log(max(30.,x))*(-6.417+log(max(30.,x))*(2.039+log(max(30.,x))*(-0.181+log(max(30.,x))*(-0.02781+log(max(30.,x))*(0.005889+log(max(30.,x))*-0.0002787)))))"},
    {"tv300pn","nef","-14.02+log(max(30.,x))*(12.55+log(max(30.,x))*(-3.805+log(max(30.,x))*(0.3161+log(max(30.,x))*(0.05375+log(max(30.,x))*(-0.01122+log(max(30.,x))*0.0005451)))))"},
    {"tv300pn","nhf","6.388+log(max(30.,x))*(-5.397+log(max(30.,x))*(1.524+log(max(30.,x))*(-0.1035+log(max(30.,x))*(-0.026+log(max(30.,x))*(0.00501+log(max(30.,x))*-0.0002474)))))"},
    // "tv300pn" checksum: a25af6468140b2a27a00c0fc1bab83df

    /*
    {"hhm3","Rjet","-1.903+log(x)*(1.073+log(x)*(-0.1386+log(x)*(-0.01856+log(x)*(0.001308+log(x)*(0.0004344+log(x)*-3.608e-05)))))"},
    {"hhm3","chf","1.164+log(x)*(-0.6436+log(x)*(0.07755+log(x)*(0.01185+log(x)*(-0.0006605+log(x)*(-0.000255+log(x)*1.741e-05)))))"},
    {"hhm3","nef","9.419+log(x)*(-9.221+log(x)*(3.072+log(x)*(-0.2995+log(x)*(-0.04007+log(x)*(0.009592+log(x)*-0.0004909)))))"},
    {"hhm3","nhf","-9.835+log(x)*(9.161+log(x)*(-2.909+log(x)*(0.2538+log(x)*(0.04149+log(x)*(-0.009093+log(x)*0.0004571)))))"},
    // "hhm3" checksum: a3cae68333a5673779cc657f95c00917
    */
    /*
    {"hhp3","Rjet","1.482+log(x)*(-0.8367+log(x)*(0.1146+log(x)*(0.0144+log(x)*(-0.00103+log(x)*(-0.0003238+log(x)*2.623e-05)))))"},
    {"hhp3","chf","-1.205+log(x)*(0.6638+log(x)*(-0.07904+log(x)*(-0.01214+log(x)*(0.0006711+log(x)*(0.0002605+log(x)*-1.765e-05)))))"},
    {"hhp3","nef","-8.908+log(x)*(8.713+log(x)*(-2.901+log(x)*(0.2831+log(x)*(0.03767+log(x)*(-0.009033+log(x)*0.0004623)))))"},
    {"hhp3","nhf","9.373+log(x)*(-8.649+log(x)*(2.721+log(x)*(-0.2341+log(x)*(-0.03887+log(x)*(0.008437+log(x)*-0.0004227)))))"},
    // "hhp3" checksum: 2d61256b68385e77d19346cad7dff15
    */
    {"hhp3","Rjet","2.817+log(x)*(-2.491+log(x)*(0.698+log(x)*(-0.04045+log(x)*(-0.008571+log(x)*(0.001366+log(x)*-5.705e-05)))))"},
    {"hhp3","chf","1.363+log(x)*(-1.602+log(x)*(0.5867+log(x)*(-0.06816+log(x)*(-0.00732+log(x)*(0.002+log(x)*-0.0001039)))))"},
    {"hhp3","nef","-7.619+log(x)*(7.721+log(x)*(-2.651+log(x)*(0.2725+log(x)*(0.03317+log(x)*(-0.00842+log(x)*0.000441)))))"},
    {"hhp3","nhf","7.178+log(x)*(-6.88+log(x)*(2.246+log(x)*(-0.1991+log(x)*(-0.03391+log(x)*(0.007537+log(x)*-0.0003856)))))"},
    // "hhp3" checksum: 57569f81de14365a30d012311898ae2a

    /*
    {"hhred103","Rjet","8.837+log(x)*(-8.669+log(x)*(2.824+log(x)*(-0.2632+log(x)*(-0.0366+log(x)*(0.008155+log(x)*-0.000392)))))"},
    {"hhred103","chf","-6.194+log(x)*(6.147+log(x)*(-2.033+log(x)*(0.1931+log(x)*(0.02723+log(x)*(-0.00626+log(x)*0.0003123)))))"},
    {"hhred103","nef","7.881+log(x)*(-7.771+log(x)*(2.609+log(x)*(-0.2595+log(x)*(-0.03347+log(x)*(0.008263+log(x)*-0.0004328)))))"},
    {"hhred103","nhf","0.9881+log(x)*(-0.6753+log(x)*(0.08679+log(x)*(0.01227+log(x)*(-0.001412+log(x)*(-0.0004022+log(x)*4.383e-05)))))"},
    // "hhred103" checksum: 907a81747461b1cf65c14dad5979198c

    {"hhred100","Rjet","8.122+log(x)*(-8.077+log(x)*(2.674+log(x)*(-0.2577+log(x)*(-0.0343+log(x)*(0.007968+log(x)*-0.0003958)))))"},
    {"hhred100","chf","-5.229+log(x)*(5.244+log(x)*(-1.751+log(x)*(0.169+log(x)*(0.02334+log(x)*(-0.005438+log(x)*0.0002729)))))"},
    {"hhred100","nef","7.64+log(x)*(-7.502+log(x)*(2.505+log(x)*(-0.2463+log(x)*(-0.03231+log(x)*(0.007848+log(x)*-0.0004064)))))"},
    {"hhred100","nhf","0.763+log(x)*(-0.5192+log(x)*(0.06223+log(x)*(0.009601+log(x)*(-0.0009671+log(x)*(-0.0003035+log(x)*3.088e-05)))))"},
    // "hhred100" checksum: 69199af5667cfa77a5ce5e81a2ce3005

    {"hhred097","Rjet","3.955+log(x)*(-4.299+log(x)*(1.536+log(x)*(-0.166+log(x)*(-0.01905+log(x)*(0.004961+log(x)*-0.0002624)))))"},
    {"hhred097","chf","0.9744+log(x)*(-0.5195+log(x)*(0.05018+log(x)*(0.009843+log(x)*(-0.0002913+log(x)*(-0.0001925+log(x)*1.092e-05)))))"},
    {"hhred097","nef","8.258+log(x)*(-8.024+log(x)*(2.641+log(x)*(-0.2521+log(x)*(-0.03449+log(x)*(0.008077+log(x)*-0.0004084)))))"},
    {"hhred097","nhf","-5.417+log(x)*(5.002+log(x)*(-1.575+log(x)*(0.1357+log(x)*(0.02288+log(x)*(-0.004978+log(x)*0.0002491)))))"},
    // "hhred097" checksum: 6ef5d3eb909154c914c4486400e3a441
    */
    /*
    {"hhblue103","Rjet","7.659+log(x)*(-8.383+log(x)*(3.072+log(x)*(-0.3504+log(x)*(-0.0365+log(x)*(0.01052+log(x)*-0.0005788)))))"},
    {"hhblue103","chf","-4.358+log(x)*(4.718+log(x)*(-1.722+log(x)*(0.1967+log(x)*(0.02094+log(x)*(-0.006097+log(x)*0.000342)))))"},
    {"hhblue103","nef","-1.369+log(x)*(0.9264+log(x)*(-0.1342+log(x)*(-0.01511+log(x)*(0.001976+log(x)*(0.0004621+log(x)*-5.121e-05)))))"},
    {"hhblue103","nhf","3.779+log(x)*(-3.806+log(x)*(1.257+log(x)*(-0.1162+log(x)*(-0.01886+log(x)*(0.004334+log(x)*-0.0002215)))))"},
    // "hhblue103" checksum: c29962961dc22e7ee360df6163004f8b
    */
    {"hhblue103","Rjet","-52.44+log(max(55.,x))*(36.13+log(max(55.,x))*(-7.941+log(max(55.,x))*(0.2816+log(max(55.,x))*(0.116+log(max(55.,x))*(-0.01448+log(max(55.,x))*0.0004894)))))"},
    {"hhblue103","chf","6.354+log(max(55.,x))*(-3.394+log(max(55.,x))*(0.3485+log(max(55.,x))*(0.06261+log(max(55.,x))*(-0.004791+log(max(55.,x))*(-0.001639+log(max(55.,x))*0.0001499)))))"},
    {"hhblue103","nef","-0.5177+log(max(55.,x))*(0.4554+log(max(55.,x))*(-0.07566+log(max(55.,x))*(-0.008602+log(max(55.,x))*(0.001275+log(max(55.,x))*(0.0002701+log(max(55.,x))*-3.214e-05)))))"},
    {"hhblue103","nhf","-5.249+log(max(55.,x))*(2.559+log(max(55.,x))*(-0.1941+log(max(55.,x))*(-0.05668+log(max(55.,x))*(0.002404+log(max(55.,x))*(0.001512+log(max(55.,x))*-0.0001229)))))"},
    // "hhblue103" checksum: b4332baa28b900ace3c09361eee1850e

    /*
    {"hhblue100","Rjet","2.58+log(x)*(-3.599+log(x)*(1.556+log(x)*(-0.2132+log(x)*(-0.01701+log(x)*(0.006162+log(x)*-0.0003663)))))"},
    {"hhblue100","chf","1.957+log(x)*(-1.221+log(x)*(0.1683+log(x)*(0.02064+log(x)*(-0.002401+log(x)*(-0.000606+log(x)*6.107e-05)))))"},
    {"hhblue100","nef","-1.127+log(x)*(0.7446+log(x)*(-0.1012+log(x)*(-0.0125+log(x)*(0.001363+log(x)*(0.0003537+log(x)*-3.549e-05)))))"},
    {"hhblue100","nhf","-1.047+log(x)*(0.592+log(x)*(-0.07548+log(x)*(-0.0108+log(x)*(0.001101+log(x)*(0.0003304+log(x)*-3.178e-05)))))"},
    // "hhblue100" checksum: f6b72f796aa50e8d4afb41193e3d403e

    {"hhblue097","Rjet","-4.789+log(x)*(2.964+log(x)*(-0.4036+log(x)*(-0.04671+log(x)*(0.005588+log(x)*(0.001409+log(x)*-0.0001475)))))"},
    {"hhblue097","chf","1.64+log(x)*(-1.035+log(x)*(0.1484+log(x)*(0.01714+log(x)*(-0.002177+log(x)*(-0.0005157+log(x)*5.458e-05)))))"},
    {"hhblue097","nef","-0.9217+log(x)*(0.5897+log(x)*(-0.07294+log(x)*(-0.01026+log(x)*(0.0008351+log(x)*(0.0002602+log(x)*-2.197e-05)))))"},
    {"hhblue097","nhf","-1.054+log(x)*(0.6373+log(x)*(-0.09353+log(x)*(-0.0111+log(x)*(0.001566+log(x)*(0.000381+log(x)*-4.376e-05)))))"},
    // "hhblue097" checksum: 890b7aa27992fbd7ceab251cbcfffc23
    */

  }};

#endif
