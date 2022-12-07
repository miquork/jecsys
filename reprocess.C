// File: reprocess.C
// Created by Mikko Voutilainen, on Sep 6th, 2012
// Updated on Jan 19, 2015 (Winter14_V8 for the paper)
// Purpose: Combine graphs from difference channels for simpler JEC analysis
//           Macro examplePlot() shows how to create plots with "CMS JEC style"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TF1.h"
#include "TMultiGraph.h"
#include "TLatex.h"
#include "TMinuit.h"
#include "TMatrixD.h"
#include "TMath.h"
#include "TLine.h"
#include "TProfile.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "tools.h"
#include "tdrstyle_mod15.C"

#include <string>
#include <iostream>
#include <vector>
#include <map>

using namespace std;

// rho used to calculate l1bias
// These values from Z(mm)+jet UL17BCDEF |eta|<1.3 at close to 202.5 GeV
// Maximilian: [MC,Data]_Rho_CHS_a30_eta_00_13_L1L2Res at 198 GeV (alpha<1.0)
// Sami:       h_Zpt_rho_alpha100_eta00_13 at 208 GeV (alpha<0.3)
const double gRhoDT = 21.39; // +/-0.04 (21.39+/-0.05) (18.41+/-0.04)
const double gRhoMC = 20.60; // +/-0.03 (20.60+/-0.09) (20.67+/-0.03)
const double gRho = 21.00; // +/-0.40 Average of gRhoDT and gRhoMC
const bool _dcsonly = false;
const bool rp_debug = true; // verbose messages

// Appy mass corrections to Z+jet (and gamma+jet)
// Update 20200325: start moving these to globalFitSyst.C as syst. eigenvectors
// Need to correct Zee and Zmm here to combine them, but photon could be later
bool useFixedFit = true; // with minitools/drawZmass.C (def:true)
double fitUncMin = 0.00000; // Some bug if unc<0?
bool correctZmmMass = true; // pT with mumu mass (def:true)
bool correctZeeMass = true; // pT with ee mass (def:true)
bool correctZMass = false; // DP_2021 // pT with Run 2 Z+jet mass (def:true)
bool correctGamMass = true; //!!UL17_V3 false, pT with ee mass at 2*pT
bool correctGamScale = false; // additional to GamMass, with value below
double valueGamScale = 1;//1.01; //1.;
bool correctUncert = false;  // ll mass uncertainty => globalFitSyst.C

// PS weight variations, 0,1,2,3 and "" for reference
string sPSWgtZ = ""; // use 'PSWeight[X]/' variant of Z+jet MC (def:"") [PSWeight0 has very poor effective statistics]
string sPSWgtG = ""; // PS weights '_ps[X]' for gamma+jet (def:"")
string sPSWgtW = ""; // PS weights '_ps[X]' for hadW (def:"")

// Which binning to use for Z+jet ("zpt" (default),"jetpt","ptave")
string zjetMode = "zpt";

// Which binning to use for multijets ("leading", "recoil" or "ptave" (default))
string multijetMode = "ptave"; // check also multijetModeS in globalFitSyst.C
bool correctMultijetLeading = false;//true; // correct for difference to JER-hybrid (UL17 true, UL18 false)
bool correctMultijetRecoilScale = true; // recoil JES relative to leading JES (due to higher gluon fraction) for data/MC ratio
double multijetRecoilScale = 1.000; // Extra data/MC difference in recoil JES
bool useMultijetRecoilScaleFunction = true; // Implement function
bool halveMultijetRecoilScaleFunction = false; // Halve effect for UL16
bool snipeMultijetBins = false; // remove individual outliers from multijets
bool snipeIncjetBins = false; // remove individual outliers from inclusive jets
bool snipeZlljetBins = false; // remove individual outliers from Zll+jet

// Which binning to use for PF jets ("tp" (default),"pt", "both", "none" (off))
string pfMode = "tp"; // "both" not supported at the moment (default: "tp")

bool confirmWarnings=false;//true; //if active, deficiencies in input files are patched after confirmation via pushing "any key to continue"
bool confirmedGamJet=false; // to avoid too many confirmations
//need to adjust corlevel in multijet.C as well!!! Not synced automatically right now

//"L1L2": "MCTruth corrections" applied, in reality L1Data/MC,L2
//"L1L2Res": MCTruth + L2Res applied
//"L1L2L3Res": MCTruth + L2L3Res applied; reference JES central values set to 1.0 (affects plotting as well)
//string CorLevel = "L1L2";    // Rarely used
//string CorLevel = "L1L2Res"; // Input to (L2)L3Res
string CorLevel = "L1L2L3Res"; // Closure test for L2L3Res


// Settings for cleaned up global fit
/////////////////////////////////////

// Zll+jet
double fzeeptmin(15.);   // Zee+jet pTmin
double fzeeptmax(800.);  // Zee+jet pTmax
double fzmmptmin(15.);   // Zmm+jet pTmin
double fzmmptmax(800.);  // Zmm+jet pTmax
double fzllmpfptmin(15); // Zll+jet (combined) MPF pTmin
double fzllmpfptmax(800);// Zll+jet (combined) MPF pTmax
double fzllbalptmin(15); // Zll+jet (combined) DB pTmin
double fzllbalptmax(800);// Zll+jet (combined) DB pTmax
double fzllptmin(15.);   // Zll+jet (combined) pTmin

// Z+jet (Sami)
double fzptmin(15.);   // Z+jet (Sami) pTmin
double fzptmax(800.);  // Z+jet (Sami) pTmax
double fzmpfptmin(15); // Z+jet (Sami) MPF pTmin
double fzmpfptmax(800);// Z+jet (Sami) MPF pTmax
double fzbalptmin(15); // Z+jet (Sami) DB pTmin
double fzbalptmax(800);// Z+jet (Sami) DB pTmax
double fzbptmax(300.); // Z+b (Sami) pTmax

// PF composition (for plots, tighter range in global fit)
double fpfjetptmin(15.);   // Dijet PF composition pTmin
double fpfjetptmax(2116.); // Dijet PF composition pTmax
double fzllpfzptmin(15);   // Z+jet PF composition pTmin
double fzllpfzptmax(230);  // Z+jet PF composition pTmax
double fzpfzptmin(15);   // Z+jet PF composition pTmin
double fzpfzptmax(230);  // Z+jet PF composition pTmax
double fppfgptmin(40);   // G+jet PF composition pTmin
double fppfgptmax(800);  // G+jet PF composition pTmax

// Inclusive jets
double fincjetptmin(21);    // Inclusive jets pTmin
double fincjetptmax(2116.); // Inclusive jets pTmin

// Hadronic W>qq'
double fhadwptamin(35);  // W>qq' pTave pTmin
double fhadwptamax(200); // W>qq' pTave pTmax 
double fhadwptbmin(40);  // W>qq' pTboth pTmin
double fhadwptbmax(175); // W>qq' pTboth pTmin

// Multijet
double fmultijetptmin(114);     // Multijet pTmin
double fmultijetptmax(2640);    // Multijet pTmax
double fmultijetptmax2(1890);   // Multijet pTmax2 (for low stats IOVs)
double fmultijetmjbptmin(114);  // Multijet MJB pTmin
double fmultijetmjbptmax(2640); // Multijet MJB pTmax
double fmultijetmpfptmin(114);  // Multijet MPF pTmin

//for fine etabins deactivate ptbal
double fdijetmpfptmin(30.); // Dijet MPF pTmin
double fdijetbalptmin(30.); // Dijet DB pTmin
double fdijetptmax(1500.);  // Dijet pT max

// Photon+jet
double fpmpfptmin(230.);  // photon+jet MPF pTmin
double fpmpfptmax(1500.); // photon+jet MPF pTmin
double fpbalptmin(230.);  // photon+jet DB pTmin
double fpbalptmax(700.);  // photon+jet DB pTmax

//minimum event counts
const double neventsmin = 20.;


// Helper functions to handle JEC
FactorizedJetCorrector *_thejec(0);
TF1 *fCorrPt(0);
FactorizedJetCorrector *getFJC(string l1, string l2="", string res="",
			       string path="");
void setEtaPtRho(FactorizedJetCorrector *jec, double eta,double pt,double rho);
Double_t funcCorrPt(Double_t *x, Double_t *p);
double getJEC(FactorizedJetCorrector *jec, double eta, double ptcorr,
	      double rho = gRho);


// put all the different methods in a single file for easy access for everybody
void reprocess(string epoch="") {

  // Set TDR style to have correct graphical setttings when storing graphs
  setTDRStyle();
 
  TDirectory *curdir = gDirectory;

  const char *cep = epoch.c_str();
  TFile *fout = new TFile(Form("rootfiles/jecdata%s.root", cep), "RECREATE");

  ////////////////////////
  // Define input files //
  ////////////////////////

  bool isUL18 = (epoch=="2018ABCD" || epoch=="2018A" || 
		 epoch=="2018B" || epoch=="2018C" || epoch=="2018D");
  bool isUL17 = (epoch=="2017BCDEF" || epoch=="2017B" || epoch=="2017C" ||
		 epoch=="2017D" || epoch=="2017E" || epoch=="2017F");
  bool isUL16 = (epoch=="2016BCDEF" || epoch=="2016BCD" || 
		 epoch=="2016EF" || epoch=="2016GH");
  bool isAPV = (epoch=="2016BCDEF" || epoch=="2016BCD" || 
		epoch=="2016EF");
  bool isRun2 = (epoch=="Run2Test");
  bool isLowPU = (epoch=="2017H");

  string sd;
  string sd18 = "../JERCProtoLab/Summer19UL18/L3Residual_";
  string sd17 = "../JERCProtoLab/Summer19UL17/L3Residual_";
  string sd16 = "../JERCProtoLab/Summer19UL16/L3Residual_";
  if (isUL18) sd = sd18;
  if (isUL17) sd = sd17;
  if (isUL16) sd = sd16; 
  const char *cd = sd.c_str();

  // Overall switches
  //if (isUL16 || isUL17) { CorLevel = "L1L2L3Res"; } // Closure

  // Multijet switches
  if (isUL16 || isUL17 || isUL18 || isRun2) { 
    fmultijetptmin = 114;//800;
    fmultijetptmax = (isRun2 ? 2500 : 2116);

    // Turn MJB "off" to avoid double counting HDM
    fmultijetmjbptmin = 114;//153;
    fmultijetmjbptmax = (isRun2 ? 2500 : 2116);

    // Recoil ranges are {0,97},{97,196},{196,220},{220,272},{272,395},{395,468},{468,592},{592,686},{686,790},{790,6500}
  }
  //if (isUL18 && CorLevel=="L1L2Res") {
  //fmultijetptmin = 330;
  //}
  if (isRun2) { // DP_2021
    fmultijetptmin = 114;
    fmultijetptmax = 2500;

    // Turn MJB "off" to avoid double counting HDM (or not)
    fmultijetmjbptmin = 114;
    fmultijetmjbptmax = 2500;
  }

  // Photon+jet switches
  if (isUL16 || isUL17 || isUL18 || isRun2) { 
    correctGamMass = true;
    correctGamScale = false;
    // Correct photon+jet scale variation within UL16. Average is still ~1.000
    //if (isUL16 &&  isAPV) { correctGamScale = true; valueGamScale = 0.995; } // latest,DP_2021?
    //if (isUL16 && !isAPV) { correctGamScale = true; valueGamScale = 1.005; } // latest,DP_2021?
    //if (false && (isUL17 || isUL18)) {
    //correctGamScale = true;
    //valueGamScale = 1.01; // ad-hoc to match Z+jet
    //if (isUL18) valueGamScale = 1;//1.005;
    //}
    // 230-300 GeV still systematically low in UL18 and others?
    fpmpfptmin = 40; // DP_2021
    //fpmpfptmin = 60; // latest; 40;
    fpmpfptmax = 1200;
    fpbalptmin = 40; // DP_2021
    //fpbalptmin = 60;//latest 40;
    fpbalptmax = 1200;
  } 
  //if (isRun2) {
  //correctGamMass = true;
  //correctGamScale = false;
  //fpmpfptmin = 60;//40;
  //fpmpfptmax = 1200;
  //fpbalptmin = 60;//40;
  //fpbalptmax = 1200;
  //}

  // Z+jet switches
  if (isUL16 && false) {
    // Z+jet pT<35 GeV biased?
    fzllmpfptmin = 30;
    if (!isAPV) fzllmpfptmin = 35;
    fzllmpfptmax = 230;
    fzllbalptmin = 30;
    if (!isAPV) fzllbalptmin = 35;
    fzllbalptmax = 175;

    // For UH same settings
    //fzmpfptmin = 30;
    //if (!isAPV) fzmpfptmin = 30;
    //fzmpfptmax = 230;
    //fzbalptmin = 30;
    //if (!isAPV) fzbalptmin = 35;
    //fzbalptmax = 175;
  }
  if ((isUL17 || isUL18) && false ) {
    fzllmpfptmin = 30;
    fzllmpfptmax = 230;
    if (isUL18) fzllmpfptmax = 300;
    fzllbalptmin = 30;
    fzllbalptmax = 230;
    if (isUL18) fzllbalptmax = 300;

    // For UH same settings
    //fzmpfptmin = 30;
    //fzmpfptmax = 230;
    //if (isUL18) fzmpfptmax = 300;
    //fzbalptmin = 30;
    //fzbalptmax = 230;
    //if (isUL18) fzbalptmax = 300;
  }
  if (isUL18 && CorLevel=="L1L2Res") {
    fzllbalptmin = 35;
    fzllbalptmax = 400;
  }
  if (isRun2) {
    fzllmpfptmin = 20;
    fzllmpfptmax = 800;
    fzllbalptmin = 20;
    fzllbalptmax = 800;
  }


  // W>qq' switches
  if (isUL16 && CorLevel=="L1L2Res") {
    fhadwptamin = 70;
    fhadwptamax = 175;
    fhadwptbmin = 60;
    fhadwptbmax = 70;
  }
  if ((isUL16 || isUL17 || isUL18) && CorLevel=="L1L2L3Res") {
    fhadwptamin = 70;
    fhadwptamax = 175;
    fhadwptbmin = 70;
    fhadwptbmax = 80;
  }
  if (isRun2) {
    fhadwptamin = 40;
    fhadwptamax = 150;//200;
    //fhadwptbmin = 70; // DP_2021
    fhadwptbmin = 80; // DP_2021
    fhadwptbmax = 80;
  }

  // Inclusive jet switches
  //if (isUL16) {
  //fincjetptmin = 21;
    //if (epoch=="2016EF") fincjetptmin = 64;
  //}
  if (isLowPU) {
    fincjetptmin = 21;
    fincjetptmax = 790;
  }


  ////////////////////////
  // Dijets             //
  ////////////////////////

  map<string,const char*> fdj_files;
  fdj_files["B"] = "B";
  fdj_files["C"] = "C";
  fdj_files["D"] = "D";
  fdj_files["E"] = "E";
  fdj_files["F"] = "F";
  fdj_files["BCDEF"] = "BCDEF";
  TFile *fdj(0), *fdj2(0);
  if (CorLevel=="L1L2Res") {
    fdj = new TFile(Form("rootfiles/L2Res_GlobalFitInput_JERSF_V2_Fall17_17Nov2017_V31_DATA/Run%s_17Nov17_2017_JERnominalV2/output/JEC_L2_Dijet_AK4PFchs_pythia8_eta_0_13.root","BCDEF"),"READ"); // regular
    assert(fdj && !fdj->IsZombie());
    fdj2 = new TFile(Form("rootfiles/L2Res_GlobalFitInput_JERSF_V2_Fall17_17Nov2017_V31_DATA/Run%s_17Nov17_2017_JERnominalV2/output/JEC_L2_Dijet_AK4PFchs_pythia8.root","BCDEF"),"READ"); // narrow bin
    assert(fdj2 && !fdj2->IsZombie());
  }

  if(CorLevel=="L1L2L3Res"){
    if (isUL18) {
      // skip for now
    }
    else if (isUL17) { // UL17
      // skip for now
    }
    else if (isUL16) {
      // skip or now
    }
    else if (isRun2) {
      // skip for now
    }
    else if (isLowPU) {
      // skip for now
    }
    else
      assert(false);
  }

  if(CorLevel!="L1L2"&&CorLevel!="L1L2L3Res"){
    cout << Form("Dijet files are only available for Autumn18-V17 (L1L2+L2L3Res) narrow/wide bins. L1L2L3Res used for L1L2Res as placeholder. Please confirm by pushing any key.") << endl;
    if(confirmWarnings)cin.ignore();
  }

  
  //////////////////////////
  // Multijets            //
  //////////////////////////

  map<string,const char*> fm_files;
  fm_files["2018A"] = "A";
  fm_files["2018B"] = "B";
  fm_files["2018C"] = "C";
  fm_files["2018D"] = "D";
  fm_files["2018ABCD"] = "ABCD";
  fm_files["2017BCDEF"] = "BCDEF";
  fm_files["2016BCD"] = "BCD";
  fm_files["2016EF"] = "EF";
  fm_files["2016GH"] = "FGH";
  fm_files["2016BCDEF"] = "BCDEF";
  fm_files["2016BCDEFGH"] = "BCDEFGH";
  TFile *fmj(0);
  if (CorLevel=="L1L2Res") {
    if (isUL18) {
      fmj = new TFile(Form("rootfiles/multijet_Rebin2_20201207_UL2018%s_jecV5_jerV2.root",fm_files[epoch]),"READ");
    }
    else if (isUL17) {
      fmj = new TFile(Form("rootfiles/multijet_Rebin2_20201209_UL2017%s_jecV6_jerV3.root",fm_files[epoch]),"READ");
    }
    else if (isUL16) {
      fmj = new TFile(Form("../JERCProtoLab/Summer19UL16/L3Residual_multijet/multijet_UL2016%s_jecV5_jerV2.root",fm_files[epoch]),"READ");
    }
    else if (isRun2) {
      assert(false);
    }
    else
      assert(false);
  }

  if(CorLevel=="L1L2L3Res") {
    if (isUL18) {
      fmj = new TFile(Form("rootfiles/multijet_Rebin2_20201207_UL2018%s_jecV5_jerV2.root",fm_files[epoch]),"READ");
    }
    else if (isUL17) {
      fmj = new TFile(Form("rootfiles/multijet_Rebin2_20201209_UL2017%s_jecV6_jerV3.root",fm_files[epoch]),"READ");
    }
    else if (isUL16) {
      //fmj = new TFile(Form("../JERCProtoLab/Summer19UL16/L3Residual_multijet/multijet_UL2016%s_jecV5_jerV2_closure.root",fm_files[epoch]),"READ");
      fmj = new TFile(Form("../JERCProtoLab/Summer19UL16/L3Residual_multijet/multijet_UL2016%s_jecV7_jerV2_closure.root",fm_files[epoch]),"READ");
    }
    else if (isRun2) {
      // Inputs created with recombine.C
      fmj = new TFile("rootfiles/jecdataRun2TestData.root","READ");
    }
    else if (isLowPU) {
      // Skip for now
    }
    else
      assert(false);
  }
  assert((fmj && !fmj->IsZombie()) || isLowPU);


  ////////////////////////
  // Inclusive jets     //
  ////////////////////////

  map<string,const char*> fij_eras;
  fij_eras["2018A"] = "A";
  fij_eras["2018B"] = "B";
  fij_eras["2018C"] = "C";
  fij_eras["2018D"] = "D";
  //fij_eras["2018ABCD"] = "";
  fij_eras["2018ABCD"] = "UL18"; // vsUL5x
  //fij_eras["2017BCDEF"] = "";
  fij_eras["2017BCDEF"] = "UL17"; // vsUL5X
  fij_eras["2016BCD"] = "BCD";
  fij_eras["2016EF"] = "EF";
  //fij_eras["2016BCDEF"] = "BCDEF";
  fij_eras["2016BCDEF"] = "UL16APV"; // vsUL5X
  //fij_eras["2016GH"] = "GH";
  fij_eras["2016GH"] = "UL16GH"; // vsUL5X
  //fij_eras["Run2Test"] = "UL4X";
  fij_eras["Run2Test"] = "UL5X"; // vsUL5X
  fij_eras["2017H"] = "UL17H";
  TFile *fij(0);
  if (CorLevel=="L1L2Res" || CorLevel=="L1L2L3Res") {
    /*
    if (isUL18) {
      fij = new TFile("rootfiles/drawDeltaJEC_18UL_JECV3.root","READ");
    }
    else if (isUL17) {
      fij = new TFile("rootfiles/drawDeltaJEC_17UL_V2M4res_cp2_all_v4.root","READ");
    }
    else if (isUL16) {
      fij = new TFile("../JERCProtoLab/Summer19UL16/L3Residual_incljet/drawDeltaJEC_16ULJECV7_vsBCDEF.root","READ"); // was V5_vsBCDandEF
    }
    else if (isRun2) {
      fij = new TFile("rootfiles/drawDeltaJEC_Run2_vsUL4X.root","READ");
    }
    else if (isLowPU) {
      fij = new TFile("rootfiles/drawDeltaJEC_Run2_vsUL4X.root","READ");
    }
    */
    if (isUL16 || isUL17 || isUL18 || isRun2 || isLowPU) {
      fij = new TFile("rootfiles/drawDeltaJEC_Run2_vsUL5X.root","READ");
    }
    else
      assert(false);
  }
  assert((fij && !fij->IsZombie()));// || isLowPU);


  ////////////////////////////
  // PF composition (dijet) //
  ////////////////////////////

  map<string,const char*> fpf_files;
  fpf_files["2018A"] = "A";
  fpf_files["2018B"] = "B";
  fpf_files["2018C"] = "C";
  fpf_files["2018D"] = "D";
  fpf_files["2018ABCD"] = "ABCD";
  fpf_files["2017BCDEF"] = "BCDEF";
  fpf_files["2016BCD"] = "BCD";
  fpf_files["2016EF"] = "EF";
  fpf_files["2016GH"] = "GH";
  fpf_files["2016BCDEF"] = "BCDEF";
  TFile *fpfdt(0), *fpfmc(0);
  if (CorLevel=="L1L2Res" || CorLevel=="L1L2L3Res") {
    if (isUL18) {
      //fpfdt = new TFile(Form("rootfiles/output-DATA-2b-UL18V2V3-%s.root",
      //		     fpf_files[epoch]),"READ"); // DP_2021
      //fpfmc = new TFile(Form("rootfiles/output-MCNU-2b-UL18V2V3-%s.root",
      //		     fpf_files[epoch]),"READ"); // DP_2021
      /*
      fpfdt = new TFile(Form("rootfiles/output-DATA-2b-UL18V5V2_%s.root",
			     fpf_files[epoch]),"READ");
      fpfmc = new TFile(Form("rootfiles/output-MC-2b-UL18V5V2_%s.root",
			     fpf_files[epoch]),"READ");
      */
      /*
      fpfdt = new TFile("rootfiles/output-DATA-2b-UL2018-PFtest.root","READ");
      fpfmc = new TFile("rootfiles/output-MC-2b-UL2018-PFtest.root","READ");
      */
      fpfdt = new TFile("rootfiles/output-DATA-2b-UL18V5V2_ABCD-pThat.root","READ"); // latest
      fpfmc = new TFile("rootfiles/output-MC-2b-UL18V5V2_ABCD-pThat.root","READ"); // latest
    }
    else if (isUL17) {
      //fpfdt = new TFile(Form("rootfiles/output-DATA-2b-UL17V4_%s.root",
      //		     fpf_files[epoch]),"READ"); // DP_2021
      //fpfmc = new TFile(Form("rootfiles/output-MCNU-2b-UL17V4_%s.root",
      //		     fpf_files[epoch]),"READ"); // DP_2021
      fpfdt = new TFile(Form("rootfiles/output-DATA-2b-UL17V5V3_%s.root",
			     fpf_files[epoch]),"READ"); // latest
      fpfmc = new TFile(Form("rootfiles/output-MC-2b-UL17V5V3_%s.root",
			     fpf_files[epoch]),"READ"); // latest
    }
    else if (isUL16 && !isAPV) { // old MC-UL16V2V1, DATA-UL16V5V2
      //fpfdt = new TFile(Form("rootfiles/output-DATA-2b-UL16V5V2_%s.root",
      //                     fpf_files[epoch]),"READ"); // DP_2021
      //fpfmc = new TFile(Form("rootfiles/output-MCNU-2b-UL16V2V1_%s.root",
      //                     fpf_files[epoch]),"READ"); // DP_2021
      fpfdt = new TFile(Form("rootfiles/output-DATA-2b-UL16V7V3_%s.root",
			     fpf_files[epoch]),"READ"); // latest
      /*
      fpfmc = new TFile(Form("rootfiles/output-MCNU-2b-UL16V7V3_%s.root",
			     fpf_files[epoch]),"READ");
      */
      fpfmc = new TFile(Form("rootfiles/output-MC-2b-UL16V7V3_%s.root",
			     fpf_files[epoch]),"READ"); // latest
    }
    else if (isUL16 && isAPV) { // old MC-UL16V3V1, DATA-UL16V5V2
      //fpfdt = new TFile(Form("rootfiles/output-DATA-2b-UL16V5V2_%s.root",
      //                     fpf_files[epoch]),"READ"); // DP_2021
      //fpfmc = new TFile(Form("rootfiles/output-MCNU-2b-UL16V3V1_%s.root",
      //                     fpf_files[epoch]),"READ"); // DP_2021
      fpfdt = new TFile(Form("rootfiles/output-DATA-2b-UL16V7V3_%s.root",
			     fpf_files[epoch]),"READ"); // latest
      /*
      fpfmc = new TFile(Form("rootfiles/output-MCNU-2b-UL16V7V3_%s.root",
			     fpf_files[epoch]),"READ");
      */
      fpfmc = new TFile(Form("rootfiles/output-MC-2b-UL16V7V3_%s.root",
			     fpf_files[epoch]),"READ"); // latest
    }
    else if (isRun2) {
      fpfdt = new TFile("rootfiles/jecdataRun2TestData.root","READ");
      fpfmc = new TFile("rootfiles/jecdataRun2TestData.root","READ");
    }
    else if (isLowPU) {
      //fpfdt = new TFile("rootfiles/output-DATA-2b-SevgiMar1.root","READ");
      //fpfmc = new TFile("rootfiles/output-MC-2b-SevgiMar1.root","READ");
      // With new pThat bins 5-10 and 10-15 added
      fpfdt = new TFile("rootfiles/output-DATA-2b-Sevgi-pfcomp_lowpu.root","READ");
      fpfmc = new TFile("rootfiles/output-MC-2b-Sevgi-pfcomp_lowpu.root","READ");
    }
    else
      assert(false);
  }
  assert((fpfdt && !fpfdt->IsZombie()) || pfMode=="none");// || isLowPU);
  assert((fpfmc && !fpfmc->IsZombie()) || pfMode=="none");// || isLowPU);


  ////////////////////////////
  // Hadronic W (W>qq')     //
  ////////////////////////////

  map<string,const char*> fw_files;
  fw_files["2018A"] = "A";
  fw_files["2018B"] = "B";
  fw_files["2018C"] = "C";
  fw_files["2018D"] = "D";
  fw_files["2018ABCD"] = "ABCD";
  fw_files["2017BCDEF"] = "BCDEF";
  fw_files["2016BCD"] = "BCD";
  fw_files["2016EF"] = "EF";
  fw_files["2016BCDEF"] = "APV";
  fw_files["2016GH"] = "GH";
  TFile *fw(0), *fw2(0);
  if (CorLevel=="L1L2Res") {
    if (isUL18) {
      fw = new TFile(Form("rootfiles/hadwUL18%s.root",fw_files[epoch]),"READ");
    }
    else if (isUL17) {
      fw = new TFile("rootfiles/hadW.root","READ"); // no eras, yet
    }
    else if (isUL16 && !isAPV) {
      fw = new TFile("../JERCProtoLab/Summer19UL16/L3Residual_W/hadW16V2_Glu_v3.root","READ");
    }
    else if (isUL16 && isAPV) {
      fw = new TFile(Form("../JERCProtoLab/Summer19UL16/L3Residual_W/hadW16%sV3_Glu.root",fw_files[epoch]),"READ");
    }
    else if (isRun2) {
      assert(false);
    }
    else
      assert(false);
  }

  if (CorLevel=="L1L2L3Res") {
    // NB: newer files use |eta|<2.5 instead of |eta|<1.3
    if (isUL18) {
      //fw = new TFile("rootfiles/hadW18V5_MPDGcorrNoW.root","READ"); // fp02
      //fw = new TFile("rootfiles/hadW1718V5_MPDGcorrNoW.root","READ"); // fp02
      fw = new TFile("rootfiles/hadW18V5_JEC.root","READ"); // fp001, |eta|<1.3 // DP_2021
      //fw = new TFile(Form("rootfiles/hadW18V5_JECv2%s.root",sPSWgtW.c_str()),
      //	     "READ"); // fp001, |eta|<1.3, // latest
    }
    else if (isUL17) {
      //fw = new TFile("rootfiles/hadWUL17V5_MPDGcorrNoW.root","READ"); // fitprob001
      //fw = new TFile("rootfiles/hadW17V5_MPDGcorrNoW.root","READ"); // fp02
      //fw = new TFile("rootfiles/hadW17V5_Glu.root","READ"); // fp001
      fw = new TFile("rootfiles/hadW17V5_JEC.root","READ"); // fp001, |eta|<1.3 // DP_2021
      //fw = new TFile(Form("rootfiles/hadW17V5_JECv2%s.root",sPSWgtW.c_str()),
      //	     "READ"); // fp001, |eta|<1.3 // latest
    }
    else if (isUL16) {
      //fw = new TFile(Form("rootfiles/hadW16%sV5_Glu.root",fw_files[epoch]),
      fw = new TFile(Form("rootfiles/hadW16%sV7_JEC.root",fw_files[epoch]), "READ"); // DP_2021
      //fw = new TFile(Form("rootfiles/hadW16%sV7_JECv2%s.root",fw_files[epoch],
      //  sPSWgtW.c_str()),"READ");
    }
    else if (isRun2) {
      fw = new TFile("rootfiles/jecdataRun2TestData.root");
      //fw = new TFile(Form("rootfiles/jecdataRun2TestData%s.root",
      //		  sPSWgtW.c_str()),"READ");
    }
    else if (isLowPU) {
      // skip for now
    }
    else
      assert(false);
  }

  assert((fw && !fw->IsZombie()) || isLowPU);


  ////////////////////////////
  // Photon+jet             //
  ////////////////////////////

  map<string,const char*> fp_files;
  fp_files["2018A"] = "A";
  fp_files["2018B"] = "B";
  fp_files["2018C"] = "C";
  fp_files["2018D"] = "D";
  fp_files["2018ABCD"] = "ABCD";
  fp_files["2017BCDEF"] = "BCDEF";
  fp_files["2016BCD"] = "BCD";
  fp_files["2016EF"] = "EF";
  fp_files["2016GH"] = "FGH";
  fp_files["2016BCDEF"] = "BCDEF";
  TFile *fp(0);
  if (CorLevel=="L1L2Res") {
    if (isUL18) {
      fp = new TFile(Form("%sgamma/Gjet_combinationfile_only_L2Res_%s_only_L2Res_JEC-v4_Data-v2_MC.root",cd,fp_files[epoch]),"READ");
    }
    else if (isUL17) {
      fp = new TFile(Form("rootfiles/2020-04-02/SimpleL1_only_L2Res/Gjet_combinationfile_SimpleL1_only_L2Res_%s_SimpleL1_only_L2Res.root",fp_files[epoch]),"READ");
    }
    else if (isUL16 && !isAPV) {
      fp = new TFile(Form("../JERCProtoLab/Summer19UL16/L3Residual_gamma/V1/UL16NonAPVGjet_combinationfile_L2L3Res_%s_L2L3Res_V1Bis.root",fp_files[epoch]),"READ");
    }
    else if (isUL16 && isAPV) {
      if (epoch=="2016BCDEF") 
	fp = new TFile(Form("../JERCProtoLab/Summer19UL16/L3Residual_gamma/V1/UL16APVGjet_combinationfile_L2L3Res_%s_L2L3Res_V1.root",fp_files[epoch]),"READ");
      else if (epoch=="2016EF" || epoch=="2016BCD") 
	fp = new TFile(Form("../JERCProtoLab/Summer19UL16/L3Residual_gamma/V1/UL16APVSplitGjet_combinationfile_L2L3Res_%s_L2L3Res_V1Bis.root",fp_files[epoch]),"READ");
      else
	assert(false);
    }
    else if (isRun2) {
      assert(false);
    }
    else
      assert(false);
  } // L2L3Res

  if(CorLevel=="L1L2L3Res"){
    if (isUL18) {
      //fp = new TFile(Form("%sgamma/Gjet_combinationfile_L2L3Res_%s_L2L3Res_JEC-v4_Data-v2_MC.root",cd,fp_files[epoch]),"READ");
      //fp = new TFile("../gamjet/files/GamHistosRatio_2018ABCD_P8_v18.root"); // DP_2021
      //fp = new TFile("../gamjet/files/GamHistosRatio_2018ABCD_P8_v19.root");
      //fp = new TFile("../gamjet/files/GamHistosRatio_2018ABCD_P8QCD_v19.root");
      fp = new TFile("../gamjet/files/GamHistosRatio_2018ABCD_P8QCD_v20.root"); // latest
    }
    else if (isUL17) {
      //fp = new TFile(Form("rootfiles/Gjet_combinationfile_17_Nov_V31_L2L3res_%s.root", fp_files[epoch]),"READ");
      //fp = new TFile("../gamjet/files/GamHistosRatio_2017BCDEF_P8_v18.root"); // DP_2021
      //fp = new TFile("../gamjet/files/GamHistosRatio_2017BCDEF_P8_v19.root");
      //fp = new TFile("../gamjet/files/GamHistosRatio_2017BCDEF_P8QCD_v19.root");
      fp = new TFile("../gamjet/files/GamHistosRatio_2017BCDEF_P8QCD_v20.root"); // latest
    }
    else if (isUL16 && !isAPV) {
      //fp = new TFile(Form("../JERCProtoLab/Summer19UL16/L3Residual_gamma/V1Closure/UL16NonAPVGjet_combinationfile_L2L3Res_%s_L2L3Res_V1Closure.root",fp_files[epoch]),"READ");
      //fp = new TFile(Form("../JERCProtoLab/Summer19UL16/L3Residual_gamma/V7Closure/ClosureV7UL16NONAPVFGHjet_combinationfile_L2L3Res_%s_L2L3Res_V2.root",fp_files[epoch]),"READ");
      //fp = new TFile("../gamjet/files/GamHistosRatio_2016FGH_P8_v18.root"); // DP_2021
      //fp = new TFile("../gamjet/files/GamHistosRatio_2016FGH_P8_v19.root");
      //fp = new TFile("../gamjet/files/GamHistosRatio_2016FGH_P8QCD_v19.root");
      fp = new TFile("../gamjet/files/GamHistosRatio_2016FGH_P8QCD_v20.root"); // latest
    }
    else if (isUL16 && isAPV) {
      if (epoch=="2016BCDEF") 
	//fp = new TFile(Form("../JERCProtoLab/Summer19UL16/L3Residual_gamma/V1Closure/UL16APVGjet_combinationfile_L2L3Res_%s_L2L3Res_V1Closure.root",fp_files[epoch]),"READ");
	//fp = new TFile(Form("../JERCProtoLab/Summer19UL16/L3Residual_gamma/V7Closure/ClosureV7UL16APVBCDEFjet_combinationfile_L2L3Res_%s_L2L3Res_V2.root",fp_files[epoch]),"READ");
	//fp = new TFile("../gamjet/files/GamHistosRatio_2016BCDEF_P8APV_v18.root"); // DP_2021
	//fp = new TFile("../gamjet/files/GamHistosRatio_2016BCDEF_P8APV_v19.root");
	//fp = new TFile("../gamjet/files/GamHistosRatio_2016BCDEF_P8QCDAPV_v19.root");
	fp = new TFile("../gamjet/files/GamHistosRatio_2016BCDEF_P8QCDAPV_v20.root"); // latest
      else if (epoch=="2016EF" || epoch=="2016BCD")  {
	fp = new TFile(Form("../JERCProtoLab/Summer19UL16/L3Residual_gamma/V1Closure/UL16APVSplitGjet_combinationfile_L2L3Res_%s_L2L3Res_V1Closure.root",fp_files[epoch]),"READ"); // NB: old!
	assert(false); // don't trip on old files
      }
      else
	assert(false);
    }
    else if (isRun2) {
      fp = new TFile("rootfiles/jecdataRun2TestData.root","RAED");
    }
    else if (isLowPU) {
      // skip for now
    }
    else
      assert(false);
  } // L1L2L3Res

  assert((fp && !fp->IsZombie()) || isLowPU);


  ////////////////////////////
  // Z(ee,mumu)+jet / KIT   //
  ////////////////////////////

  map<string,const char*> fz_files;
  fz_files["2018A"] = "A";
  fz_files["2018B"] = "B";
  fz_files["2018C"] = "C";
  fz_files["2018D"] = "D";
  fz_files["2018ABCD"] = "ABCD";
  fz_files["2017BCDEF"] = "BCDEF";
  fz_files["2016BCD"] = "preVFPBCD";
  fz_files["2016EF"] = "preVFPEFearly";
  fz_files["2016BCDEF"] = "preVFPBCDEFearly";
  fz_files["2016GH"] = "postVFPFlateGH";
  TFile *fzmm(0), *fzee(0);
  if (CorLevel=="L1L2Res") {
    if (isUL18) {
      // With UL18V4, pT bin split
      fzmm = new TFile(Form("%sZ/JEC_Combination_Zmm/splitZPtBin70/ZJetCombination_Zmm_DYJets_Madgraph_12Nov2019_Summer19UL18_JECV4_L1L2Res.root",cd),"READ");
      fzee = new TFile(Form("%sZ/JEC_Combination_Zee/splitZPtBin70/ZJetCombination_Zee_DYJets_Madgraph_12Nov2019_Summer19UL18_JECV4_L1L2Res.root",cd),"READ");
    }
    else if (isUL17) {
      fzmm = new TFile("../JERCProtoLab/Summer19UL17/L3Residual_Z/JEC_Combination_Zmm/ZJetCombination_Zmm_DYJets_Madgraph_09Aug2019_Summer19UL17_JECV5_L1L2Res.root","READ");
      fzee = new TFile("../JERCProtoLab/Summer19UL17/L3Residual_Z/JEC_Combination_Zee/ZJetCombination_Zee_DYJets_Madgraph_09Aug2019_Summer19UL17_JECV5_L1L2Res.root","READ");
    }
    else if (isUL16 && !isAPV) {
      fzmm = new TFile("../JERCProtoLab/Summer19UL16/L3Residual_Z/nonAPV/JEC_Combination_Zmm/ZJetCombination_Zmm_DYJets_amc_21Feb2020_Summer19UL16_L1L2Res.root","READ");
      fzee = new TFile("../JERCProtoLab/Summer19UL16/L3Residual_Z/nonAPV/JEC_Combination_Zee/ZJetCombination_Zee_DYJets_amc_21Feb2020_Summer19UL16_V2_L1L2Res.root","READ"); // |eta|<1.3 for electrons
      assert(false); // make sure not going here when testing V7 closure
    }
    else if (isUL16 && isAPV) {
      fzmm = new TFile("../JERCProtoLab/Summer19UL16/L3Residual_Z/APV/JEC_Combination_Zmm/ZJetCombination_Zmm_DYJets_amc_21Feb2020_Summer19UL16_V3_L1L2Res.root","READ");
      fzee = new TFile("../JERCProtoLab/Summer19UL16/L3Residual_Z/APV/JEC_Combination_Zee/electrons_in_barrel/ZJetCombination_Zee_DYJets_amc_21Feb2020_Summer19UL16_V3_L1L2Res.root","READ");
      assert(false); // make sure not going here when testing V7 closure
    }
    else if (isRun2) {
      assert(false);
    }
    else
      assert(false);
  }

  if(CorLevel=="L1L2L3Res"){
    if (isUL18) {
      fzmm = new TFile(Form("%sZ/JEC_Combination_Zmm/splitZPtBin70/ZJetCombination_Zmm_DYJets_Madgraph_12Nov2019_Summer19UL18_JECV4_L1L2L3Res.root",cd),"READ");
      fzee = new TFile(Form("%sZ/JEC_Combination_Zee/splitZPtBin70/ZJetCombination_Zee_DYJets_Madgraph_12Nov2019_Summer19UL18_JECV4_L1L2L3Res.root",cd),"READ");
    }
    else if (isUL17) {
      fzmm = new TFile("../JERCProtoLab/Summer19UL17/L3Residual_Z/JEC_Combination_Zmm/ZJetCombination_Zmm_DYJets_Madgraph_09Aug2019_Summer19UL17_JECV6_L1L2L3Res.root","READ"); // was V5
      fzee = new TFile("../JERCProtoLab/Summer19UL17/L3Residual_Z/JEC_Combination_Zee/ZJetCombination_Zee_DYJets_Madgraph_09Aug2019_Summer19UL17_JECV6_L1L2L3Res.root","READ"); // was V5
    }
    else if (isUL16 && !isAPV) {
      fzmm = new TFile("../JERCProtoLab/Summer19UL16/L3Residual_Z/nonAPV/JEC_Combination_Zmm/ZJetCombination_Zmm_DYJets_amc_21Feb2020_Summer19UL16_V7_L1L2L3Res.root","READ"); // was V5old
      fzee = new TFile("../JERCProtoLab/Summer19UL16/L3Residual_Z/nonAPV/JEC_Combination_Zee/ZJetCombination_Zee_DYJets_amc_21Feb2020_Summer19UL16_V7_L1L2L3Res.root","READ"); // was V5old (|eta|<1.3 for electrons?)
    }
    else if (isUL16 && isAPV) {
      fzmm = new TFile("../JERCProtoLab/Summer19UL16/L3Residual_Z/APV/JEC_Combination_Zmm/ZJetCombination_Zmm_DYJets_amc_21Feb2020_Summer19UL16_V7_L1L2L3Res.root","READ"); // was V5old
      fzee = new TFile("../JERCProtoLab/Summer19UL16/L3Residual_Z/APV/JEC_Combination_Zee/ZJetCombination_Zee_DYJets_amc_21Feb2020_Summer19UL16_V7_L1L2L3Res.root","READ"); // was V5old
    }
    else if (isRun2) {
      fzmm = new TFile("rootfiles/jecdataRun2TestData.root","READ");
      fzee = new TFile("rootfiles/jecdataRun2TestData.root","READ");
    }
    else if (isLowPU) {
      // Use 2017F as placeholder for 2017H (hard to switch of Z mass)
      fzmm = new TFile("../JERCProtoLab/Summer19UL17/L3Residual_Z/JEC_Combination_Zmm/ZJetCombination_Zmm_DYJets_Madgraph_09Aug2019_Summer19UL17_JECV6_L1L2L3Res.root","READ"); // place-holder
      fzee = new TFile("../JERCProtoLab/Summer19UL17/L3Residual_Z/JEC_Combination_Zee/ZJetCombination_Zee_DYJets_Madgraph_09Aug2019_Summer19UL17_JECV6_L1L2L3Res.root","READ"); // place-holder
    }
    else
      assert(false);
  }
  assert(fzmm && !fzmm->IsZombie());
  assert(fzee && !fzee->IsZombie());

  if (isUL18) {
    fzmm->cd(Form("Run2018%s",fz_files[epoch])); fzmm = (TFile*)gDirectory;
    fzee->cd(Form("Run2018%s",fz_files[epoch])); fzee = (TFile*)gDirectory;
  }
  else if (isUL17) {
    fzmm->cd(Form("Run2017%s",fz_files[epoch])); fzmm = (TFile*)gDirectory;
    fzee->cd(Form("Run2017%s",fz_files[epoch])); fzee = (TFile*)gDirectory;
  }
  else if (isUL16) {
    fzmm->cd(Form("Run2016%s",fz_files[epoch])); fzmm = (TFile*)gDirectory;
    fzee->cd(Form("Run2016%s",fz_files[epoch])); fzee = (TFile*)gDirectory;
    /*
    if (epoch=="2016GH") {
      fzmm->cd("Run2016postVFPFlateGH"); fzmm = (TFile*)gDirectory;
      fzee->cd("Run2016postVFPFlateGH"); fzee = (TFile*)gDirectory;
    }
    else if (epoch=="2016BCDEF") {
      fzmm->cd("Run2016preVFPBCDEFearly"); fzmm = (TFile*)gDirectory;
      fzee->cd("Run2016preVFPBCDEFearly"); fzee = (TFile*)gDirectory;
    }
    else if (epoch=="2016BCD") {
      fzmm->cd("Run2016preVFPBCD"); fzmm = (TFile*)gDirectory;
      fzee->cd("Run2016preVFPBCD"); fzee = (TFile*)gDirectory;
    }
    else if (epoch=="2016EF") {
      fzmm->cd("Run2016preVFPEFearly"); fzmm = (TFile*)gDirectory;
      fzee->cd("Run2016preVFPEFearly"); fzee = (TFile*)gDirectory;
    }
    else
      assert(false);
    */
  }
  else if (isRun2) {
    // do nothing extra
  }
  else if (isLowPU) {
    // Use 2017F as placeholder for 2017H (hard to switch of Z mass)
    fzmm->cd(Form("Run2017%s","F")); fzmm = (TFile*)gDirectory;
    fzee->cd(Form("Run2017%s","F")); fzee = (TFile*)gDirectory;
  }
  else
    assert(false);


  ////////////////////////////
  // Z+jet / Helsinki       //
  ////////////////////////////

  TFile *fz(0), *fzjes(0);
  TH1D *hzjes(0);
  if (CorLevel=="L1L2Res") {
    if (isUL18) {
      fz = new TFile("rootfiles/jme_bplusZ_merged_v34.root","READ"); // not L2L3Res, but has latest bells and whistles
      fzjes = new TFile("rootfiles/jecdataBCDEF_V4.root","READ");
      hzjes = (TH1D*)fzjes->Get("ratio/eta00-13/herr_ref");
    }
    else if (isUL17) {
      fz = new TFile("rootfiles/jme_bplusZ_merged_v35.root","READ"); // UL18 placeholder
      // JES correction for zjet
      fzjes = new TFile("rootfiles/jecdataBCDEF_V4.root","READ");
      hzjes = (TH1D*)fzjes->Get("ratio/eta00-13/herr_ref");
    }
    else if (isUL16) {
      // skip for now, use UL17 as place holder
      // UL17 L1L2L3Res as better placeholder (contains Z mass)
      fz = new TFile("rootfiles/jme_bplusZ_merged_v28_2017emu_wTTJets.root","READ");
      fzjes = new TFile("rootfiles/jecdataBCDEF_V4.root","READ");
      hzjes = (TH1D*)fzjes->Get("ratio/eta00-13/herr_ref");
    }
    else if (isRun2) {
      assert(false);
    }
    else
      assert(false);
  }    

  if(CorLevel=="L1L2L3Res"){
    if (isUL18) {
      //fz = new TFile("rootfiles/jme_bplusZ_merged_v35.root","READ"); // test (17+18 e+mu)
      fz = new TFile("rootfiles/jme_bplusZ_merged_v35_2018_emu_wTTJets.root","READ"); // botched PSWeights // DP_2021
      //fz = new TFile("rootfiles/jme_bplusZ_merged_v38_muon2018.root","READ"); // only a100
      //fz = new TFile("rootfiles/jme_bplusZ_merged_v43_muon2018.root","READ"); // only a100,30
      //fz = new TFile("rootfiles/jme_bplusZ_merged_v44_muon2018.root","READ"); // only a100,30 // latest
    }
    else if (isUL17) {
      //fz = new TFile("rootfiles/jme_bplusZ_merged_v28_2017emu_wTTJets.root","READ");
      fz = new TFile("rootfiles/jme_bplusZ_merged_v35_2017_emu_wTTJets.root","READ");
    }
    else if (isUL16 && !isAPV) {
      //fz = new TFile("rootfiles/jme_bplusZ_merged_v36_2016fh_emu_wTTJets.root","READ"); // was v35; v36 buggy, used BCDEF JEC for MC
      //fz = new TFile("rootfiles/jme_bplusZ_merged_vX_2016FH_mu.root","READ"); // was v35
      fz = new TFile("rootfiles/jme_bplusZ_merged_v37_2016FH_mu.root","READ"); // was v35
      //fzjes = new TFile("rootfiles/jecdataBCDEF_V4.root","READ");
      //hzjes = (TH1D*)fzjes->Get("ratio/eta00-13/herr_ref");
    }
    else if (isUL16 && isAPV) {
      fz = new TFile("rootfiles/jme_bplusZ_merged_v36_2016bf_emu_wTTJets.root","READ"); // was v35
    }
    else if (isRun2) {
      fz = new TFile("rootfiles/jecdataRun2TestData.root","READ");
    }
    else if (epoch=="2017H") { // Low PU dataset
      //fz = new TFile("rootfiles/jme_bplusZ_merged_v37_muon2017H.root","READ");
      fz = new TFile("rootfiles/jme_bplusZ_merged_v38_muon2017H.root","READ");
    }
    else
      assert(false);
  }

  assert(fz && !fz->IsZombie());
  //assert(fzjes && !fzjes->IsZombie());
  //assert(hzjes);


  TFile *fmz = fz;
  TFile *fmzee = fzee;
  TFile *fmzmm = fzmm;
  assert(fmz && !fmz->IsZombie());
  assert(fmzee && !fmzee->IsZombie());
  assert(fmzmm && !fmzmm->IsZombie());

  string sr = (epoch=="L4" ? "eta_00_24" : "eta_00_13");
  const char *cr = sr.c_str();
  const char *cl = CorLevel.c_str();
  TH1D *hmzee = (TH1D*)fmzee->Get(Form("Ratio_ZMass_CHS_a30_%s_%s",cr,cl));
  TH1D *hmzmm = (TH1D*)fmzmm->Get(Form("Ratio_ZMass_CHS_a30_%s_%s",cr,cl));
  TH1D *hmz(0);
  if (isRun2) {
    hmzee = (TH1D*)fmzee->Get("ratio/eta00-13/mass_zeejet_a30");
    hmzmm = (TH1D*)fmzmm->Get("ratio/eta00-13/mass_zmmjet_a30");
    hmz = (TH1D*)fmz->Get("ratio/eta00-13/mass_zjet_a30");
  }
  assert(hmzee);
  assert(hmzmm);
  assert(hmz || !isRun2);

  TH2D *hmz_dt2 = (TH2D*)fmz->Get(Form("data/%s/h_Zpt_mZ_alpha30",cr));
  if (!hmz_dt2) hmz_dt2 = (TH2D*)fmz->Get(Form("data/%s/h_Zpt_mZ_alpha100",cr));
  if (!hmz_dt2) hmz_dt2 = (TH2D*)fmz->Get(Form("data/%s/h_Zpt_mZ_alpha100_eta_00_13",cr));
  assert(hmz_dt2 || isRun2);
  TH1D *hmz_dt = (hmz_dt2 ? hmz_dt2->ProfileX()->ProjectionX("hmz_dt") : 0);
  TH1D *hmzee_dt = (TH1D*)fmzee->Get(Form("Data_ZMass_CHS_a30_%s_%s",cr,cl));
  TH1D *hmzmm_dt = (TH1D*)fmzmm->Get(Form("Data_ZMass_CHS_a30_%s_%s",cr,cl));
  if (isRun2) {
    hmzee_dt = (TH1D*)fmzee->Get("data/eta00-13/mass_zeejet_a30");
    hmzmm_dt = (TH1D*)fmzmm->Get("data/eta00-13/mass_zmmjet_a30");
    hmz_dt = (TH1D*)fmz->Get("data/eta00-13/mass_zjet_a30");
  }
  assert(hmzee_dt);
  assert(hmzmm_dt);
  assert(hmz_dt);

  TH2D *hmz_mc2 = (TH2D*)fmz->Get(Form("mc/%s/h_Zpt_mZ_alpha30",cr));
  if (!hmz_mc2) hmz_mc2 = (TH2D*)fmz->Get(Form("mc/%s/h_Zpt_mZ_alpha100",cr));
  if (!hmz_mc2) hmz_mc2 = (TH2D*)fmz->Get(Form("mc/%s/h_Zpt_mZ_alpha100_eta_00_13",cr));
  assert(hmz_mc2 || isRun2);
  TH1D *hmz_mc = (hmz_mc2 ? hmz_mc2->ProfileX()->ProjectionX("hmz_mc") : 0);
  TH1D *hmzee_mc = (TH1D*)fmzee->Get(Form("MC_ZMass_CHS_a30_%s_%s",cr,cl));
  TH1D *hmzmm_mc = (TH1D*)fmzmm->Get(Form("MC_ZMass_CHS_a30_%s_%s",cr,cl));
  if (isRun2) {
    hmzee_mc = (TH1D*)fmzee->Get("mc/eta00-13/mass_zeejet_a30");
    hmzmm_mc = (TH1D*)fmzmm->Get("mc/eta00-13/mass_zmmjet_a30");
    hmz_mc = (TH1D*)fmz->Get("mc/eta00-13/mass_zjet_a30");
  }
  assert(hmzee_mc);
  assert(hmzmm_mc);
  assert(hmz_mc);

  TH1D *hmzee1 = (TH1D*)fmzee->Get(Form("Ratio_ZMass_CHS_a100_%s_%s",cr,cl));
  TH1D *hmzmm1 = (TH1D*)fmzmm->Get(Form("Ratio_ZMass_CHS_a100_%s_%s",cr,cl));
  TH1D *hmz1(0);
  if (isRun2) {
    hmzee1 = (TH1D*)fmzee->Get("ratio/eta00-13/mass_zeejet_a100");
    hmzmm1 = (TH1D*)fmzmm->Get("ratio/eta00-13/mass_zmmjet_a100");
    hmz1 = (TH1D*)fmz->Get("ratio/eta00-13/mass_zjet_a100");
  }
  assert(hmzee1);
  assert(hmzmm1);
  assert(hmz1 || !isRun2);

  TH2D *hmz1_dt2 = (TH2D*)fmz->Get(Form("data/%s/h_Zpt_mZ_alpha100",cr));
  if (!hmz1_dt2) hmz1_dt2 = (TH2D*)fmz->Get(Form("data/%s/h_Zpt_mZ_alpha100_eta_00_13",cr));
  assert(hmz1_dt2 || isRun2);
  TH1D *hmz1_dt = (hmz1_dt2 ? hmz1_dt2->ProfileX()->ProjectionX("hmz1_dt") : 0);
  TH1D *hmzee1_dt =(TH1D*)fmzee->Get(Form("Data_ZMass_CHS_a100_%s_%s",cr,cl));
  TH1D *hmzmm1_dt =(TH1D*)fmzmm->Get(Form("Data_ZMass_CHS_a100_%s_%s",cr,cl));
  if (isRun2) {
    hmzee1_dt = (TH1D*)fmzee->Get("data/eta00-13/mass_zeejet_a100");
    hmzmm1_dt = (TH1D*)fmzmm->Get("data/eta00-13/mass_zmmjet_a100");
    hmz1_dt = (TH1D*)fmz->Get("data/eta00-13/mass_zjet_a100");
  }
  assert(hmzee1_dt);
  assert(hmzmm1_dt);
  assert(hmz1_dt);

  TH2D *hmz1_mc2 = (TH2D*)fmz->Get(Form("mc/%s/h_Zpt_mZ_alpha100",cr));
  if (!hmz1_mc2) hmz1_mc2 = (TH2D*)fmz->Get(Form("mc/%s/h_Zpt_mZ_alpha100_eta_00_13",cr));
  assert(hmz1_mc2 || isRun2);
  TH1D *hmz1_mc = (hmz1_mc2 ? hmz1_mc2->ProfileX()->ProjectionX("hmz1_mc") : 0);
  TH1D *hmzee1_mc = (TH1D*)fmzee->Get(Form("MC_ZMass_CHS_a100_%s_%s",cr,cl));
  TH1D *hmzmm1_mc = (TH1D*)fmzmm->Get(Form("MC_ZMass_CHS_a100_%s_%s",cr,cl));
  if (isRun2) {
    hmzee1_mc = (TH1D*)fmzee->Get("mc/eta00-13/mass_zeejet_a100");
    hmzmm1_mc = (TH1D*)fmzmm->Get("mc/eta00-13/mass_zmmjet_a100");
    hmz1_mc = (TH1D*)fmz->Get("mc/eta00-13/mass_zjet_a100");
  }
  assert(hmzee1_mc);
  assert(hmzmm1_mc);
  assert(hmz1_mc);

  if (!isRun2) {
    hmz = (TH1D*)hmz_dt->Clone("hmz");
    hmz->Divide(hmz_mc);
    hmz1 = (TH1D*)hmz1_dt->Clone("hmz1");
    hmz1->Divide(hmz1_mc);
  }
  assert(hmz);
  assert(hmz1);


  // Extra functions for gamma+jet, modified from Zee+jet
  TF1 *f1mgam = new TF1("f1mgam","[0]+[1]*log(0.01*x)+[2]*pow(log(0.01*x),2)",
			fpmpfptmin, fpmpfptmax);
  TF1 *f1egam = new TF1("f1egam","sqrt([0]+pow(log(0.01*x),2)*[1]"
                        "+pow(log(0.01*x),4)*[2]"
                        "+2*log(0.01*x)*[3]+2*pow(log(0.01*x),2)*[4]"
                        "+2*pow(log(0.01*x),3)*[5])",
			fpmpfptmin, fpmpfptmax);

  // \BEGIN copy-paste from minitools/drawZmass.C

  // Smoothen mass corrections
  TF1 *f1mzee = new TF1("f1mzee","[0]+[1]*log(0.01*x)+[2]*pow(log(0.01*x),2)",
			fzeeptmin, fzeeptmax);
  TF1 *f1ezee = new TF1("f1ezee","sqrt([0]+pow(log(0.01*x),2)*[1]"
                        "+pow(log(0.01*x),4)*[2]"
                        "+2*log(0.01*x)*[3]+2*pow(log(0.01*x),2)*[4]"
                        "+2*pow(log(0.01*x),3)*[5])",
			fzeeptmin, fzeeptmax);
  if (correctZeeMass || correctGamMass) {
    if (useFixedFit) {

      if (isUL18) {
	// UL18 Run2018ABCD fit with minitools/drawZmass.C
	//f1mzee->SetParameters(1.00153, 0.00214, -0.00012);
	//f1ezee->SetParameters( +5.5e-08,  +1.5e-07, +3.01e-07,
	//		       +2.35e-08, -7.73e-08, -4.34e-08);
	f1mzee->SetParameters(1.00132, 0.00186, 0.00061);
	f1ezee->SetParameters(+1.31e-08, +2.53e-08, +9.92e-09,
			      +8.94e-09, -1.28e-09, +1.08e-08);
	// UL18 Run2018ABCD Zee from above
	f1mgam->SetParameters(1.00132, 0.00186, 0.00061);
	f1egam->SetParameters(+1.31e-08, +2.53e-08, +9.92e-09,
			      +8.94e-09, -1.28e-09, +1.08e-08);
      }
      else if (isUL17) {
	//assert(false); // CHECK
	// UL17 RunBCDEF fit with minitools/drawZmass.C
	//f1mzee->SetParameters(0.99780, 0.00225, 0.00031);
	//f1ezee->SetParameters( +6.1e-08, +1.66e-07, +3.27e-07,
	//		       +2.56e-08,  -8.5e-08, -5.39e-08);
	// UL17 Run2017BCDEF fit with minitools/drawZmass.C
	f1mzee->SetParameters(1.00080, 0.00215, 0.00048);
	f1ezee->SetParameters(+2.58e-08, +4.18e-08,  +1.3e-08,
			      +1.84e-08, -1.88e-09, +1.38e-08);
	// UL17 RunBCDEF Zee from above
	f1mgam->SetParameters(1.00080, 0.00215, 0.00048);
	f1egam->SetParameters(+2.58e-08, +4.18e-08,  +1.3e-08,
			      +1.84e-08, -1.88e-09, +1.38e-08);
      }
      else if (isUL16 && !isAPV) {
	// V2 with |eta|<1.3 elextrons and no p2 constraint
	// UL16 2016GH fit with minitools/drawZmass.C
	f1mzee->SetParameters(1.00175, 0.00349, 0.00159);
	f1ezee->SetParameters(+1.09e-07, +2.61e-07, +1.37e-07,
			      +5.34e-08, -2.42e-08, +1.48e-07);
	// UL16 Run2016GH from above
	f1mgam->SetParameters(1.00175, 0.00349, 0.00159);
	f1egam->SetParameters(+1.09e-07, +2.61e-07, +1.37e-07,
			      +5.34e-08, -2.42e-08, +1.48e-07);
      }
      else if (isUL16 && isAPV) {
	// UL16 2016BCDEF fit with minitools/drawZmass.C
	f1mzee->SetParameters(1.00122, 0.00256, 0.00098);
	f1ezee->SetParameters(+1.05e-07, +2.65e-07, +1.33e-07,
			      +5.45e-08, -2.02e-08,  +1.5e-07);
	// UL16 Run2016BCDEF from above Zee
	f1mgam->SetParameters(1.00122, 0.00256, 0.00098);
	f1egam->SetParameters(+1.05e-07, +2.65e-07, +1.33e-07,
			      +5.45e-08, -2.02e-08,  +1.5e-07);
      }
      else if (isRun2) { // 137.6/fb

	// Run2Test fit with minitools/drawZmass.C
	f1mzee->SetParameters(1.00101, 0.00206, 0.00055);
	f1ezee->SetParameters(+7.14e-09, +1.54e-08, +6.75e-09,
			      +4.57e-09, -8.73e-10, +7.45e-09);
	
	// Run2Test Zee from above
	f1mgam->SetParameters(1.00101, 0.00206, 0.00055);
	f1egam->SetParameters(+7.14e-09, +1.54e-08, +6.75e-09,
			      +4.57e-09, -8.73e-10, +7.45e-09);
      }
      else if (isLowPU) {
	// Use 2017BCDEF as placeholder
	// UL17 Run2017BCDEF fit with minitools/drawZmass.C
	f1mzee->SetParameters(1.00080, 0.00215, 0.00048);
	f1ezee->SetParameters(+2.58e-08, +4.18e-08,  +1.3e-08,
			      +1.84e-08, -1.88e-09, +1.38e-08);
	// UL17 RunBCDEF Zee from above
	f1mgam->SetParameters(1.00080, 0.00215, 0.00048);
	f1egam->SetParameters(+2.58e-08, +4.18e-08,  +1.3e-08,
			      +1.84e-08, -1.88e-09, +1.38e-08);
      }
      else
	assert(false);
      // No correction
      //f1mgam->SetParameters(1,0,0);
    }
    else
      hmzee->Fit(f1mzee);

  }

  TF1 *f1mz = new TF1("f1mz","[0]+[1]*log(0.01*x)+[2]*pow(log(0.01*x),2)",
		      fzptmin, fzptmax);
  TF1 *f1ez = new TF1("f1ez","sqrt([0]+pow(log(0.01*x),2)*[1]"
		      "+pow(log(0.01*x),4)*[2]"
		      "+2*log(0.01*x)*[3]+2*pow(log(0.01*x),2)*[4]"
		      "+2*pow(log(0.01*x),3)*[5])",
		      fzptmin, fzptmax);
  if (correctZMass) {
      // Run2Test fit with minitools/drawZmass.C
    f1mz->SetParameters(0.99875, 0.00118, 0.00059);
    f1ez->SetParameters(+2.51e-09, +8.73e-09, +3.99e-09,
			+2.15e-09, +7.07e-11, +4.95e-09);
  }

  TF1 *f1mzmm = new TF1("f1mzmm","[0]+[1]*log(0.01*x)+[2]*pow(log(0.01*x),2)",
			fzmmptmin, fzmmptmax);
  TF1 *f1ezmm = new TF1("f1ezmm","sqrt([0]+pow(log(0.01*x),2)*[1]"
                        "+pow(log(0.01*x),4)*[2]"
                        "+2*log(0.01*x)*[3]+2*pow(log(0.01*x),2)*[4]"
                        "+2*pow(log(0.01*x),3)*[5])",
			fzeeptmin, fzeeptmax);
  if (correctZmmMass) {
    if (useFixedFit) {
      
      if (isUL18) {
	// UL18 ABCD fit with minitools/drawZmass.C
	//f1mzmm->SetParameters(0.99839, 0.00000, 0.00000);
	//f1ezmm->SetParameters(+1.64e-08,        +0,        +0,
	//		             +0,        +0,        +0);	
	f1mzmm->SetParameters(0.99796, 0.00000, 0.00000);
	f1ezmm->SetParameters(+2.36e-09,        +0,        +0,
			      +0,        +0,        +0);
      }
      else if (isUL17) {
	//assert(false); // CHECK
	// UL17 BCDEF fit with minitools/drawZmass.C
	//f1mzmm->SetParameters(0.99821, 0.00000, 0.00000);
	//f1ezmm->SetParameters(+3.43e-08,        +0,        +0,
	//		             +0,        +0,        +0);
	// UL17 2017BCDEF fit with minitools/drawZmass.C
	f1mzmm->SetParameters(0.99889, 0.00000, 0.00000);
	f1ezmm->SetParameters(+4.82e-09,        +0,        +0,
			      +0,        +0,        +0);
      }
      else if (isUL16 && !isAPV) {
	// UL16 2016GH fit with minitools/drawZmass.C
	f1mzmm->SetParameters(1.00017, 0.00000, 0.00000);
	f1ezmm->SetParameters(+1.21e-08,        +0,        +0,
			      +0,        +0,        +0);
      }
      else if (isUL16 && isAPV) {
	// UL16 2016BCDEF fit with minitools/drawZmass.C
	f1mzmm->SetParameters(1.00053, 0.00000, 0.00000);
	f1ezmm->SetParameters(+1.17e-08,        +0,        +0,
			      +0,        +0,        +0);
      }
      else if (isRun2) {

	// Run2Test fit with minitools/drawZmass.C
	f1mzmm->SetParameters(0.99888, 0.00000, 0.00000);
	f1ezmm->SetParameters(+1.31e-09,        +0,        +0,
			      +0,        +0,        +0);
      }
      else if (isLowPU) {
	// Use 2017BCDEF as placeholder
	// UL17 2017BCDEF fit with minitools/drawZmass.C
	f1mzmm->SetParameters(0.99889, 0.00000, 0.00000);
	f1ezmm->SetParameters(+4.82e-09,        +0,        +0,
			      +0,        +0,        +0);
      }
      else
	assert(false);
    }
    else
      hmzmm->Fit(f1mzmm);
  }

  // \END copy-paste from minitools/drawZmass.C


  TF1 *f2 = new TF1("f2","[0]*pow(x,2)+[1]",30,3000);
  f2->SetParameters(0,0);

  if (correctMultijetLeading && multijetMode=="leading") {
    f2->SetParameters(1.105e-07, 0.03063); // V4

    // Fits done in minitools/systMultijet.C on UL17BCDEF
    //f2 Chi2/NDF = 15.2 / 19
    f2->SetParameters(1.462e-07, -0.001596); // 20200419-SimpleL1 MPF
  }

  // MultijetRecoilScale from minitools/drawFlavor.C
  TF1 *f1mjb = new TF1("f1mjb","[p0]+[p1]*pow(x,[p2])",114,2500);
  //f1mjb->SetParameters(1.001, 8.262, -1.315); // chi2/NDF=0.0/2383
  // New Run2Test update for Flavor.C
  f1mjb->SetParameters(1.002, 99.07, -1.97); // chi2/NDF=2718.9/2383 // DP_2021
  // f1zmjb update to account for gJES bias in parameterized Z+jet
  //f1mjb->SetParameters(1.002, 21.69, -1.705); // chi2/NDF=1187.1/2383 // latest


  // Link to Z+jet 2D distribution for JEC calculations
  // This is used for correctly averaging JEC and its uncertainty
  // for the wide eta bins used in global fit combinations
  // => need to update for Run 2
  // On 28 Feb 2014, at 14:00, from Dominik Haitz
  // Re: pT~200 GeV
  TFile *feta = new TFile("rootfiles/th2d_jeteta_zpt.root","READ");
  assert(feta && !feta->IsZombie());

  TH2D *h2eta = (TH2D*)feta->Get("data"); assert(h2eta);

  // New 2D pT-eta distribution down to 5 GeV based on P8 dijet sample
  TFile *feta2 = new TFile("rootfiles/P8_dijet_45M_TH2D_correct_weighting.root","READ");
  assert(feta2 && !feta2->IsZombie());
  TH2D *h2pteta = (TH2D*)feta2->Get("pTjet_etajet"); assert(h2pteta);
  //TFile *feta2 = new TFile("rootfiles/RunIISummer20UL16NanoAODAPVv2_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_LeadingJet2DDistribution.root","READ");
  //assert(feta2 && !feta2->IsZombie());
  //TH2D *h2pteta = (TH2D*)feta2->Get("Dist"); assert(h2pteta);


  // Store pointers to all files in a map for easy handling later
  map<string, TFile*> files;
  files["pfjet_data"] = fpfdt;
  files["pfjet_mc"] = fpfmc;
  files["dijet"] = fdj;
  files["incjet"] = fij;
  files["hadw"] = fw;
  files["multijet"] = fmj;
  files["gamjet"] = fp;
  files["zeejet"] = fzee;
  files["zmmjet"] = fzmm;
  files["zjet"] = fz;
  // Flavor response variants (Z+jet)
  files["zi"] = fz;
  files["zb"] = fz;
  files["zc"] = fz;
  files["zq"] = fz;
  files["zg"] = fz;
  files["zn"] = fz;
  files["zii"] = fz;
  files["zib"] = fz;
  files["zic"] = fz;
  files["ziq"] = fz;
  files["ziu"] = fz; //new
  files["zis"] = fz; //new
  files["zig"] = fz;
  files["zin"] = fz;
  files["zbi"] = fz;
  files["zbb"] = fz;
  files["zbc"] = fz;
  files["zbq"] = fz;
  files["zbg"] = fz;
  files["zbn"] = fz;
  files["zci"] = fz;
  files["zcb"] = fz;
  files["zcc"] = fz;
  files["zcq"] = fz;
  files["zcg"] = fz;
  files["zcn"] = fz;
  files["zqi"] = fz;
  files["zqb"] = fz;
  files["zqc"] = fz;
  files["zqq"] = fz;
  files["zqg"] = fz;
  files["zqn"] = fz;
  files["zgi"] = fz;
  files["zgb"] = fz;
  files["zgc"] = fz;
  files["zgq"] = fz;
  files["zgg"] = fz;
  files["zgn"] = fz;
  files["zni"] = fz;
  files["znb"] = fz;
  files["znc"] = fz;
  files["znq"] = fz;
  files["zng"] = fz;
  files["znn"] = fz;
  /*
  // Flavor response variants (photon+jet)
  files["gi"] = fp;
  files["gb"] = fp;
  files["gc"] = fp;
  files["gq"] = fp;
  files["gg"] = fp;
  files["gn"] = fp;
  files["gii"] = fp;
  files["gib"] = fp;
  files["gic"] = fp;
  files["giq"] = fp;
  files["gig"] = fp;
  files["gin"] = fp;
  files["gbi"] = fp;
  files["gbb"] = fp;
  files["gbc"] = fp;
  files["gbq"] = fp;
  files["gbg"] = fp;
  files["gbn"] = fp;
  files["gci"] = fp;
  files["gcb"] = fp;
  files["gcc"] = fp;
  files["gcq"] = fp;
  files["gcg"] = fp;
  files["gcn"] = fp;
  files["gqi"] = fp;
  files["gqb"] = fp;
  files["gqc"] = fp;
  files["gqq"] = fp;
  files["gqg"] = fp;
  files["gqn"] = fp;
  files["ggi"] = fp;
  files["ggb"] = fp;
  files["ggc"] = fp;
  files["ggq"] = fp;
  files["ggg"] = fp;
  files["ggn"] = fp;
  files["gni"] = fp;
  files["gnb"] = fp;
  files["gnc"] = fp;
  files["gnq"] = fp;
  files["gng"] = fp;
  files["gnn"] = fp;
  */

  // Map variable names used in different files to common scheme
  map<string, map<string, const char*> > rename;

  // Note: non-CHS results are not available for all samples, so not used
  rename["dijet"]["mpfchs"] = "mpfchs";
  rename["dijet"]["mpfchs1"] = "mpfchs";
  rename["dijet"]["ptchs"] = "ptchs";

  // Minsuk in Run II
  rename["multijet"]["ratio"] = "Data";// => PATCH
  rename["multijet"]["data"] = "Data";
  rename["multijet"]["mc"] = "MC";
  rename["multijet"]["ratiocrecoil"] = "CRecoil_L1L2Res";
  rename["multijet"]["datacrecoil"] = "CRecoil_L1L2Res";
  rename["multijet"]["mccrecoil"] = "CRecoil_MG";//P8";
  // Switch to leading jet pT binning instead or recoil pT binning
  if (multijetMode=="leading") {
    if (isUL18) {
      rename["multijet"]["ratiocrecoil"] = "CRecoil_leading_L1L2res";
      rename["multijet"]["datacrecoil"] = "CRecoil_leading_L1L2res";
      rename["multijet"]["mccrecoil"] = "CRecoil_leading_MG";
      //
      rename["multijet"]["ratiompfchs1"] = "MPF_leading_L1L2res";
      rename["multijet"]["ratioptchs"] = "MJB_leading_L1L2res";
      rename["multijet"]["datampfchs1"] = "MPF_leading_L1L2res";
      rename["multijet"]["dataptchs"] = "MJB_leading_L1L2res";
      rename["multijet"]["mcmpfchs1"] = "MPF_leading_MG";
      rename["multijet"]["mcptchs"] = "MJB_leading_MG";
      rename["multijet"]["mccounts"] = "h2_MPF_leading_MG";

      if (CorLevel=="L1L2L3Res") {
	// Minsuk v2
	rename["multijet"]["ratiocrecoil"] = "CRecoil_leading_L1L2L3Res";
	rename["multijet"]["datacrecoil"] = "CRecoil_leading_L1L2L3Res";
	rename["multijet"]["ratiompfchs1"] = "MPF_leading_L1L2L3Res";
	rename["multijet"]["ratioptchs"] = "MPF1_leading_L1L2L3Res"; //MJB
	rename["multijet"]["datampfchs1"] = "MPF_leading_L1L2L3Res";
	rename["multijet"]["dataptchs"] = "MPF1_leading_L1L2L3Res"; //MJB
	rename["multijet"]["ratiocounts"] = "h2_MPF_leading_L1L2L3Res";
	rename["multijet"]["datacounts"] = "h2_MPF_leading_L1L2L3Res";
	rename["multijet"]["mcptchs"] = "MPF1_leading_MG"; //MJB

	// new stuff v2
	rename["multijet"]["ratiompf1"] = "MPF1_leading_L1L2L3Res";
	rename["multijet"]["ratiompfn"] = "MPFn_leading_L1L2L3Res";
	rename["multijet"]["ratiompfu"] = "MPFu_leading_L1L2L3Res";//v3
	rename["multijet"]["datampf1"] = "MPF1_leading_L1L2L3Res";
	rename["multijet"]["datampfn"] = "MPFn_leading_L1L2L3Res";
	rename["multijet"]["datampfu"] = "MPFu_leading_L1L2L3Res"; //v3
	rename["multijet"]["mcmpf1"] = "MPF1_leading_MG";
	rename["multijet"]["mcmpfn"] = "MPFn_leading_MG";
	rename["multijet"]["mcmpfu"] = "MPFu_leading_MG"; // v3
      }
    }
    else if (isUL17) {
      rename["multijet"]["ratiocrecoil"] = "CRecoil_leading_L1L2Res";
      rename["multijet"]["datacrecoil"] = "CRecoil_leading_L1L2Res";
      rename["multijet"]["mccrecoil"] = "CRecoil_leading_MG";
      //
      rename["multijet"]["ratiompfchs1"] = "MPF_leading_L1L2Res";
      rename["multijet"]["ratioptchs"] = "MJB_leading_L1L2Res";
      rename["multijet"]["datampfchs1"] = "MPF_leading_L1L2Res";
      rename["multijet"]["dataptchs"] = "MJB_leading_L1L2Res";
      rename["multijet"]["mcmpfchs1"] = "MPF_leading_MG";
      rename["multijet"]["mcptchs"] = "MJB_leading_MG";
    }
    else
      assert(false);
  }
  else if (multijetMode=="recoil") {
    if (isUL18) {
      rename["multijet"]["ratiompfchs1"] = "MPF_recoil_L1L2res";
      rename["multijet"]["ratioptchs"] = "MJB_recoil_L1L2res";
      rename["multijet"]["ratiocrecoil"] = "CRecoil_L1L2res";
      rename["multijet"]["datampfchs1"] = "MPF_recoil_L1L2res";
      rename["multijet"]["dataptchs"] = "MJB_recoil_L1L2res";
      rename["multijet"]["datacrecoil"] = "CRecoil_L1L2res";
      rename["multijet"]["mcmpfchs1"] = "MPF_recoil_MG";
      rename["multijet"]["mcptchs"] = "MJB_recoil_MG";
      rename["multijet"]["mccrecoil"] = "CRecoil_MG";
      rename["multijet"]["mccounts"] = "h2_MPF_leading_MG";

      if (CorLevel=="L1L2L3Res") {

	rename["multijet"]["ratiompfchs1"] = "MPF_recoil_L1L2L3Res";
	rename["multijet"]["ratioptchs"] = "MPF1_recoil_L1L2L3Res"; //MJB
	rename["multijet"]["ratiocrecoil"] = "CRecoil_L1L2L3Res";
	rename["multijet"]["datampfchs1"] = "MPF_recoil_L1L2L3Res";
	rename["multijet"]["dataptchs"] = "MPF1_recoil_L1L2L3Res"; //MJB
	rename["multijet"]["datacrecoil"] = "CRecoil_L1L2L3Res";
	rename["multijet"]["mcmpfchs1"] = "MPF_recoil_MG";
	rename["multijet"]["mcptchs"] = "MPF1_recoil_MG"; //MJB
	rename["multijet"]["mccrecoil"] = "CRecoil_MG";

	rename["multijet"]["ratiocounts"] = "h2_MPF_leading_L1L2L3Res";
	rename["multijet"]["datacounts"] = "h2_MPF_leading_L1L2L3Res";

	// new stuff v3
	rename["multijet"]["ratiompf1"] = "MPF1_recoil_L1L2L3Res";
	rename["multijet"]["ratiompfn"] = "MPFn_recoil_L1L2L3Res";
	rename["multijet"]["ratiompfu"] = "MPFu_recoil_L1L2L3Res";
	rename["multijet"]["datampf1"] = "MPF1_recoil_L1L2L3Res";
	rename["multijet"]["datampfn"] = "MPFn_recoil_L1L2L3Res";
	rename["multijet"]["datampfu"] = "MPFu_recoil_L1L2L3Res";
	rename["multijet"]["mcmpf1"] = "MPF1_recoil_MG";
	rename["multijet"]["mcmpfn"] = "MPFn_recoil_MG";
	rename["multijet"]["mcmpfu"] = "MPFu_recoil_MG";
      }

    }
    else if (isUL17) {
      rename["multijet"]["ratiompfchs1"] = "MPF_recoil_L1L2Res";
      rename["multijet"]["ratioptchs"] = "MJB_recoil_L1L2Res";
      rename["multijet"]["datampfchs1"] = "MPF_recoil_L1L2Res";
      rename["multijet"]["dataptchs"] = "MJB_recoil_L1L2Res";
      rename["multijet"]["mcmpfchs1"] = "MPF_recoil_MG";
      rename["multijet"]["mcptchs"] = "MJB_recoil_MG";
    }
    else
      assert(false);
  }
  else if (multijetMode=="ptave") {
    //if (isUL18 || isUL16) {
    if (isUL16) {
      rename["multijet"]["ratiocrecoil"] = "CRecoil_ptave_L1L2res";
      rename["multijet"]["datacrecoil"] = "CRecoil_ptave_L1L2res";
      rename["multijet"]["mccrecoil"] = "CRecoil_ptave_MG";
      //
      rename["multijet"]["ratiompfchs1"] = "MPF_ptave_L1L2res";
      rename["multijet"]["ratioptchs"] = "MPF1_ptave_L1L2res";
      rename["multijet"]["datampfchs1"] = "MPF_ptave_L1L2res";
      rename["multijet"]["dataptchs"] = "MPF1_ptave_L1L2res";
      rename["multijet"]["mcmpfchs1"] = "MPF_ptave_MG";
      rename["multijet"]["mcptchs"] = "MPF1_ptave_MG";

      rename["multijet"]["ratiocounts"] = "Stats_ptave_L1L2res";
      rename["multijet"]["datacounts"] = "Stats_ptave_L1L2res";
      rename["multijet"]["mccounts"] = "Stats_ptave_MG";

      if (CorLevel=="L1L2Res") {
      
	// new stuff v3+UL16
	rename["multijet"]["ratiompf1"] = "MPF1_ptave_L1L2res";
	rename["multijet"]["ratiompfn"] = "MPFn_ptave_L1L2res";
	rename["multijet"]["ratiompfu"] = "MPFu_ptave_L1L2res";
	rename["multijet"]["ratiorho"] = "prho_ptave_L1L2res";
	rename["multijet"]["datampf1"] = "MPF1_ptave_L1L2res";
	rename["multijet"]["datampfn"] = "MPFn_ptave_L1L2res";
	rename["multijet"]["datampfu"] = "MPFu_ptave_L1L2res";
	rename["multijet"]["datarho"] = "prho_ptave_L1L2res";
	rename["multijet"]["mcmpf1"] = "MPF1_ptave_MG";
	rename["multijet"]["mcmpfn"] = "MPFn_ptave_MG";
	rename["multijet"]["mcmpfu"] = "MPFu_ptave_MG";
	rename["multijet"]["mcrho"] = "prho_ptave_MG";
      }

      if (CorLevel=="L1L2L3Res") { 
	// Minsuk v2
	rename["multijet"]["ratiocrecoil"] = "CRecoil_ptave_L1L2L3res";
	rename["multijet"]["ratiompfchs1"] = "MPF_ptave_L1L2L3res";
	rename["multijet"]["ratioptchs"] = "MPF1_ptave_L1L2L3res";//MJB
	rename["multijet"]["datacrecoil"] = "CRecoil_ptave_L1L2L3res";
	rename["multijet"]["datampfchs1"] = "MPF_ptave_L1L2L3res";
	rename["multijet"]["dataptchs"] = "MPF1_ptave_L1L2L3res";//MJB
	rename["multijet"]["mcptchs"] = "MPF1_ptave_MG";//MJB
	
	// new stuff v3
	rename["multijet"]["ratiompf1"] = "MPF1_ptave_L1L2L3res";
	rename["multijet"]["ratiompfn"] = "MPFn_ptave_L1L2L3res";
	rename["multijet"]["ratiompfu"] = "MPFu_ptave_L1L2L3res";
	rename["multijet"]["ratiorho"] = "prho_ptave_L1L2L3res";
	rename["multijet"]["datampf1"] = "MPF1_ptave_L1L2L3res";
	rename["multijet"]["datampfn"] = "MPFn_ptave_L1L2L3res";
	rename["multijet"]["datampfu"] = "MPFu_ptave_L1L2L3res";
	rename["multijet"]["datarho"] = "prho_ptave_L1L2L3res";
	rename["multijet"]["mcmpf1"] = "MPF1_ptave_MG";
	rename["multijet"]["mcmpfn"] = "MPFn_ptave_MG";
	rename["multijet"]["mcmpfu"] = "MPFu_ptave_MG";
	rename["multijet"]["mcrho"] = "prho_ptave_MG";
      }
    }
    else if (isUL17 || isUL18) {
      rename["multijet"]["ratiocrecoil"] = "CRecoil_ptave_L1L2Res";
      rename["multijet"]["datacrecoil"] = "CRecoil_ptave_L1L2Res";
      rename["multijet"]["mccrecoil"] = "CRecoil_ptave_MG";
      //
      rename["multijet"]["ratiompfchs1"] = "MPF_ptave_L1L2Res";
      rename["multijet"]["ratioptchs"] = "MJB_ptave_L1L2Res";
      rename["multijet"]["datampfchs1"] = "MPF_ptave_L1L2Res";
      rename["multijet"]["dataptchs"] = "MJB_ptave_L1L2Res";
      rename["multijet"]["mcmpfchs1"] = "MPF_ptave_MG";
      rename["multijet"]["mcptchs"] = "MJB_ptave_MG";
      rename["multijet"]["ratiocounts"] = "h2_MPF0_ptave_L1L2Res";
      rename["multijet"]["mccounts"] = "h2_MPF0_ptave_MG";
      rename["multijet"]["datacounts"] = "h2_MPF0_ptave_L1L2Res";

      rename["multijet"]["ratiompf1"] = "MPF1_ptave_L1L2Res";
      rename["multijet"]["ratiompfn"] = "MPFn_ptave_L1L2Res";
      rename["multijet"]["ratiompfu"] = "MPFu_ptave_L1L2Res";
      rename["multijet"]["datampf1"] = "MPF1_ptave_L1L2Res";
      rename["multijet"]["datampfn"] = "MPFn_ptave_L1L2Res";
      rename["multijet"]["datampfu"] = "MPFu_ptave_L1L2Res";
      rename["multijet"]["mcmpf1"] = "MPF1_ptave_MG";
      rename["multijet"]["mcmpfn"] = "MPFn_ptave_MG";
      rename["multijet"]["mcmpfu"] = "MPFu_ptave_MG";

      if (CorLevel=="L1L2L3Res") {
	rename["multijet"]["ratiocrecoil"] = "CRecoil_ptave_L1L2L3Res";
	rename["multijet"]["ratiompfchs1"] = "MPF_ptave_L1L2L3Res";
	rename["multijet"]["ratioptchs"] = "MJB_ptave_L1L2L3Res";//MJB,MPF1
	rename["multijet"]["datacrecoil"] = "CRecoil_ptave_L1L2L3Res";
	rename["multijet"]["datampfchs1"] = "MPF_ptave_L1L2L3Res";
	rename["multijet"]["dataptchs"] = "MJB_ptave_L1L2L3Res";//MJB,MPF1
	rename["multijet"]["mcptchs"] = "MJB_ptave_MG";//MJB,MPF1
	
	// Use MPF0 or MPFx? We don't have MPFx. Probably either is fine
	rename["multijet"]["ratiocounts"] = "h2_MPF0_ptave_L1L2L3Res";
	rename["multijet"]["datacounts"] = "h2_MPF0_ptave_L1L2L3Res";

	rename["multijet"]["ratiompf1"] = "MPF1_ptave_L1L2L3Res";
	rename["multijet"]["ratiompfn"] = "MPFn_ptave_L1L2L3Res";
	rename["multijet"]["ratiompfu"] = "MPFu_ptave_L1L2L3Res";
	rename["multijet"]["datampf1"] = "MPF1_ptave_L1L2L3Res";
	rename["multijet"]["datampfn"] = "MPFn_ptave_L1L2L3Res";
	rename["multijet"]["datampfu"] = "MPFu_ptave_L1L2L3Res";
	rename["multijet"]["mcmpf1"] = "MPF1_ptave_MG";
	rename["multijet"]["mcmpfn"] = "MPFn_ptave_MG";
	rename["multijet"]["mcmpfu"] = "MPFu_ptave_MG";
      }
    }
    else if (isRun2) {
      // nothing extra
    }
    else if (isLowPU) {
      // nothing extra
    }
    else
      assert(false);
  }
  else
    assert(false);

  /*
  if (isUL18) {
    rename["incjet"]["mpfchs1"] = "jet_Run18UL%s_det";
    rename["incjet"]["ptchs"] = "jet_Run18UL%s_det"; // copy
  }
  else if (isUL17) {
    rename["incjet"]["mpfchs1"] = "jet_Run17UL%s_fwd3"; // JER V2 + JES fix
    rename["incjet"]["ptchs"] = "jet_Run17UL%s_fwd3"; // copy
  }
  else if (isUL16) {
    rename["incjet"]["mpfchs1"] = "jet_Run16UL%s_fwd";
    rename["incjet"]["ptchs"] = "jet_Run16UL%s_fwd";

  }
  else if (isRun2) {
    rename["incjet"]["mpfchs1"] = "jet_Run%s_fwd";
    rename["incjet"]["ptchs"] = "jet_Run%s_fwd";
  }
  else if (isLowPU) {
    rename["incjet"]["mpfchs1"] = "jet_Run%s_fwd";
    rename["incjet"]["ptchs"] = "jet_Run%s_fwd";
  }
  */
  if (isUL16 || isUL17 || isUL18 || isRun2 || isLowPU) {
    rename["incjet"]["mpfchs1"] = "jet_Run%s_fwd";
    rename["incjet"]["ptchs"] = "jet_Run%s_fwd";
  }
  else
    assert(false);

  rename["hadw"]["mpfchs1"] = "mass_ptave";
  rename["hadw"]["ptchs"] = "mass_ptboth";
  rename["hadw"]["counts"] = "nevents_ptave";
  rename["hadw"]["rho"] = "rho_ptave";

  rename["gamjet"]["ratio"] = "";
  rename["gamjet"]["data"] = "_DATA"; 
  rename["gamjet"]["mc"] = "_MC"; 
  rename["gamjet"]["mpfchs"] = "resp_MPFchs";
  rename["gamjet"]["mpfchs1"] = "resp_MPFchs"; 
  //rename["gamjet"]["ptchs"] = "resp_PtBalchs";  // v18 DP_2021
  rename["gamjet"]["ptchs"] = "resp_DBchs"; // v19 latest
  rename["gamjet"]["counts"] = "RawNEvents_data_vs_pt";
  //
  rename["gamjet"]["mpf1"] = "resp_MPFR1chs";
  rename["gamjet"]["mpfn"] = "resp_MPFRnchs";
  rename["gamjet"]["mpfu"] = "resp_MpfRuchs";
  rename["gamjet"]["rho"] = "resp_Rho_CHS"; // added in v20
  //rename["gamjet"]["rho"] = "prhovspt"; // v19

  rename["zeejet"]["ratio"] = "Ratio";
  rename["zeejet"]["data"] = "Data";
  rename["zeejet"]["mc"] = "MC";
  rename["zeejet"]["mpfchs"] = "MPF-notypeI_CHS";
  rename["zeejet"]["mpfchs1"] = "MPF_CHS";
  //rename["zeejet"]["ptchs"] = (isUL18 ? "MPF-leading-jet_CHS" : "PtBal_CHS");
  rename["zeejet"]["ptchs"] = "MPF-leading-jet_CHS";
  rename["zeejet"]["counts"] = "RawNEvents_CHS";
  rename["zeejet"]["chf"] = "PFjetCHF_CHS";
  rename["zeejet"]["nef"] = "PFjetPF_CHS";
  rename["zeejet"]["nhf"] = "PFjetNHF_CHS";
  rename["zeejet"]["cef"] = "PFjetEF_CHS";
  rename["zeejet"]["muf"] = "PFjetMF_CHS";
  //
  rename["zeejet"]["mpf1"] = "MPF-leading-jet_CHS";
  rename["zeejet"]["mpfn"] = "MPF-subleading-jets_CHS";
  rename["zeejet"]["mpfu"] = "MPF-unclustered_CHS";
  rename["zeejet"]["rho"] = "Rho_CHS";

  rename["zmmjet"]["ratio"] = "Ratio";
  rename["zmmjet"]["data"] = "Data";
  rename["zmmjet"]["mc"] = "MC";
  rename["zmmjet"]["mpfchs"] = "MPF-notypeI_CHS";
  rename["zmmjet"]["mpfchs1"] = "MPF_CHS"; 
  //rename["zmmjet"]["ptchs"] = (isUL18 ? "MPF-leading-jet_CHS" : "PtBal_CHS"); 
  rename["zmmjet"]["ptchs"] = "MPF-leading-jet_CHS";
  rename["zmmjet"]["counts"] = "RawNEvents_CHS";
  rename["zmmjet"]["chf"] = "PFjetCHF_CHS";
  rename["zmmjet"]["nef"] = "PFjetPF_CHS";
  rename["zmmjet"]["nhf"] = "PFjetNHF_CHS";
  rename["zmmjet"]["cef"] = "PFjetEF_CHS";
  rename["zmmjet"]["muf"] = "PFjetMF_CHS";
  //
  rename["zmmjet"]["mpf1"] = "MPF-leading-jet_CHS";
  rename["zmmjet"]["mpfn"] = "MPF-subleading-jets_CHS";
  rename["zmmjet"]["mpfu"] = "MPF-unclustered_CHS";
  rename["zmmjet"]["rho"] = "Rho_CHS";

  // Results from Sami's Z+b analysis
  rename["zjet"]["ratio"] = "data"; // missing => PATCH
  rename["zjet"]["data"] = "data";
  rename["zjet"]["mc"] = "mc";
  bool isNewName = (isUL16 && !isAPV);// DP_2021 ((isUL16 && !isAPV) || isLowPU || isUL18);
  rename["zjet"]["mpfchs"] = (isNewName ? "rmpf" : "mpfchs");
  rename["zjet"]["mpfchs1"] = (isNewName ? "rmpf" : "mpfchs"); 
  rename["zjet"]["ptchs"] = (isNewName ? "rmpfjet1" : "mpfjet1");//"ptchs"; 
  rename["zjet"]["counts"] = (isNewName ? "statistics_rmpf" : "statistics_mpfchs");
  rename["zjet"]["chf"] = "chHEF";
  rename["zjet"]["nef"] = "neEmEF";
  rename["zjet"]["nhf"] = "neHEF";
  rename["zjet"]["cef"] = "chEmEF";
  rename["zjet"]["muf"] = "muEF";
  //
  rename["zjet"]["mpf1"] = (isNewName ? "rmpfjet1" : "mpfjet1");
  rename["zjet"]["mpfn"] = (isNewName ? "rmpfjetn" : "mpfjetn");
  rename["zjet"]["mpfu"] = (isNewName ? "rmpfuncl" : "mpfuncl");
  rename["zjet"]["rho"] = "rho";

  // Flavor response variants (Z+jet)
  rename["zi"]["tag"] = "";
  rename["zb"]["tag"] = "_btagdeepbtight";
  rename["zc"]["tag"] = "_btagdeepctight";
  rename["zq"]["tag"] = "_quarktag";
  rename["zg"]["tag"] = "_gluontag";
  rename["zn"]["tag"] = "_notag";
  rename["zi"]["parent"] = "zjet";
  rename["zb"]["parent"] = "zjet";
  rename["zc"]["parent"] = "zjet";
  rename["zq"]["parent"] = "zjet";
  rename["zg"]["parent"] = "zjet";
  rename["zn"]["parent"] = "zjet";
  rename["zi"]["gen"] = "";
  rename["zb"]["gen"] = "";//_genb";
  rename["zc"]["gen"] = "";//_genc";
  rename["zq"]["gen"] = "";//_genuds";
  rename["zg"]["gen"] = "";//_geng";
  rename["zn"]["gen"] = "";//_geng";
  //
  rename["zii"]["tag"] = "";
  rename["zib"]["tag"] = "";
  rename["zic"]["tag"] = "";
  rename["ziq"]["tag"] = "";
  rename["ziu"]["tag"] = "";//_quarktag"; //new (but botched in v37)
  rename["zis"]["tag"] = "";//_quarktag"; //new (but botched in v37)
  rename["zig"]["tag"] = "";
  rename["zin"]["tag"] = "";
  rename["zii"]["parent"] = "zjet";
  rename["zib"]["parent"] = "zjet";
  rename["zic"]["parent"] = "zjet";
  rename["ziq"]["parent"] = "zjet";
  rename["ziu"]["parent"] = "zjet"; //new
  rename["zis"]["parent"] = "zjet"; //new
  rename["zig"]["parent"] = "zjet";
  rename["zin"]["parent"] = "zjet";
  rename["zii"]["gen"] = "";
  rename["zib"]["gen"] = "_genb";
  rename["zic"]["gen"] = "_genc";
  rename["ziq"]["gen"] = "_genuds";
  rename["ziu"]["gen"] = "_genud"; //new
  rename["zis"]["gen"] = "_gens"; //new
  rename["zig"]["gen"] = "_geng";
  rename["zin"]["gen"] = "_unclassified";
  //
  rename["zbi"]["tag"] = "_btagdeepbtight";
  rename["zbb"]["tag"] = "_btagdeepbtight";
  rename["zbc"]["tag"] = "_btagdeepbtight";
  rename["zbq"]["tag"] = "_btagdeepbtight";
  rename["zbg"]["tag"] = "_btagdeepbtight";
  rename["zbn"]["tag"] = "_btagdeepbtight";
  rename["zbi"]["parent"] = "zjet";
  rename["zbb"]["parent"] = "zjet";
  rename["zbc"]["parent"] = "zjet";
  rename["zbq"]["parent"] = "zjet";
  rename["zbg"]["parent"] = "zjet";
  rename["zbn"]["parent"] = "zjet";
  rename["zbi"]["gen"] = "";
  rename["zbb"]["gen"] = "_genb";
  rename["zbc"]["gen"] = "_genc";
  rename["zbq"]["gen"] = "_genuds";
  rename["zbg"]["gen"] = "_geng";
  rename["zbn"]["gen"] = "_unclassified";
  //
  rename["zci"]["tag"] = "_btagdeepctight";
  rename["zcb"]["tag"] = "_btagdeepctight";
  rename["zcc"]["tag"] = "_btagdeepctight";
  rename["zcq"]["tag"] = "_btagdeepctight";
  rename["zcg"]["tag"] = "_btagdeepctight";
  rename["zcn"]["tag"] = "_btagdeepctight";
  rename["zci"]["parent"] = "zjet";
  rename["zcb"]["parent"] = "zjet";
  rename["zcc"]["parent"] = "zjet";
  rename["zcq"]["parent"] = "zjet";
  rename["zcg"]["parent"] = "zjet";
  rename["zcn"]["parent"] = "zjet";
  rename["zci"]["gen"] = "";
  rename["zcb"]["gen"] = "_genb";
  rename["zcc"]["gen"] = "_genc";
  rename["zcq"]["gen"] = "_genuds";
  rename["zcg"]["gen"] = "_geng";
  rename["zcn"]["gen"] = "_unclassified";
  //
  rename["zqi"]["tag"] = "_quarktag";
  rename["zqb"]["tag"] = "_quarktag";
  rename["zqc"]["tag"] = "_quarktag";
  rename["zqq"]["tag"] = "_quarktag";
  rename["zqg"]["tag"] = "_quarktag";
  rename["zqn"]["tag"] = "_quarktag";
  rename["zqi"]["parent"] = "zjet";
  rename["zqb"]["parent"] = "zjet";
  rename["zqc"]["parent"] = "zjet";
  rename["zqq"]["parent"] = "zjet";
  rename["zqg"]["parent"] = "zjet";
  rename["zqn"]["parent"] = "zjet";
  rename["zqi"]["gen"] = "";
  rename["zqb"]["gen"] = "_genb";
  rename["zqc"]["gen"] = "_genc";
  rename["zqq"]["gen"] = "_genuds";
  rename["zqg"]["gen"] = "_geng";
  rename["zqn"]["gen"] = "_unclassified";
  //
  rename["zgi"]["tag"] = "_gluontag";
  rename["zgb"]["tag"] = "_gluontag";
  rename["zgc"]["tag"] = "_gluontag";
  rename["zgq"]["tag"] = "_gluontag";
  rename["zgg"]["tag"] = "_gluontag";
  rename["zgn"]["tag"] = "_gluontag";
  rename["zgi"]["parent"] = "zjet";
  rename["zgb"]["parent"] = "zjet";
  rename["zgc"]["parent"] = "zjet";
  rename["zgq"]["parent"] = "zjet";
  rename["zgg"]["parent"] = "zjet";
  rename["zgn"]["parent"] = "zjet";
  rename["zgi"]["gen"] = "";
  rename["zgb"]["gen"] = "_genb";
  rename["zgc"]["gen"] = "_genc";
  rename["zgq"]["gen"] = "_genuds";
  rename["zgg"]["gen"] = "_geng";
  rename["zgn"]["gen"] = "_unclassified";
  //
  rename["zni"]["tag"] = "_notag";
  rename["znb"]["tag"] = "_notag";
  rename["znc"]["tag"] = "_notag";
  rename["znq"]["tag"] = "_notag";
  rename["zng"]["tag"] = "_notag";
  rename["znn"]["tag"] = "_notag";
  rename["zni"]["parent"] = "zjet";
  rename["znb"]["parent"] = "zjet";
  rename["znc"]["parent"] = "zjet";
  rename["znq"]["parent"] = "zjet";
  rename["zng"]["parent"] = "zjet";
  rename["znn"]["parent"] = "zjet";
  rename["zni"]["gen"] = "";
  rename["znb"]["gen"] = "_genb";
  rename["znc"]["gen"] = "_genc";
  rename["znq"]["gen"] = "_genuds";
  rename["zng"]["gen"] = "_geng";
  rename["znn"]["gen"] = "_unclassified";

  // Photon+jet flavor
  rename["gi"]["parent"] = rename["gb"]["parent"] = rename["gc"]["parent"] =
  rename["gq"]["parent"] = rename["gg"]["parent"] = rename["gn"]["parent"] =
  rename["gii"]["parent"] = rename["gib"]["parent"] = rename["gic"]["parent"] =
  rename["giq"]["parent"] = rename["gig"]["parent"] = rename["gin"]["parent"] =
  rename["gbi"]["parent"] = rename["gbb"]["parent"] = rename["gbc"]["parent"] =
  rename["gbq"]["parent"] = rename["gbg"]["parent"] = rename["gbn"]["parent"] =
  rename["gci"]["parent"] = rename["gcb"]["parent"] = rename["gcc"]["parent"] =
  rename["gcq"]["parent"] = rename["gcg"]["parent"] = rename["gcn"]["parent"] =
  rename["gqi"]["parent"] = rename["gqb"]["parent"] = rename["gqc"]["parent"] =
  rename["gqq"]["parent"] = rename["gqg"]["parent"] = rename["gqn"]["parent"] =
  rename["ggi"]["parent"] = rename["ggb"]["parent"] = rename["ggc"]["parent"] =
  rename["ggq"]["parent"] = rename["ggg"]["parent"] = rename["ggn"]["parent"] =
  rename["gni"]["parent"] = rename["gnb"]["parent"] = rename["gnc"]["parent"] =
  rename["gnq"]["parent"] = rename["gng"]["parent"] = rename["gnn"]["parent"] =
    "gamjet";

  // color and style codes
  map<string, map<string, int> > style;

  style["pfjet"]["chf"] = kFullCircle;
  style["pfjet"]["nhf"] = kFullDiamond;
  style["pfjet"]["nef"] = kFullSquare;
  style["pfjet"]["cef"] = kFullDiamond;
  style["pfjet"]["muf"] = kFullDiamond;
  style["pfjet"]["puf"] = kFullDiamond;
  style["pfjet_mc"]["chf"] = kOpenCircle;
  style["pfjet_mc"]["nhf"] = kOpenDiamond;
  style["pfjet_mc"]["nef"] = kOpenSquare;
  style["pfjet_mc"]["cef"] = kOpenDiamond;
  style["pfjet_mc"]["muf"] = kOpenDiamond;
  style["pfjet_mc"]["puf"] = kOpenDiamond;
  style["dijet"]["mpfchs1"] = kFullDiamond;
  style["dijet"]["ptchs"] = kOpenDiamond;
  style["incjet"]["mpfchs1"] = kFullDiamond;
  style["incjet"]["ptchs"] = kOpenDiamond;
  style["hadw"]["mpfchs1"] = kFullCircle;
  style["hadw"]["ptchs"] = kOpenCircle;
  style["multijet"]["mpfchs1"] = kFullTriangleUp;
  style["multijet"]["ptchs"] = kOpenTriangleUp;
  style["multijet"]["mpf1"] = kFullTriangleUp;
  style["multijet"]["mpfn"] = kFullTriangleUp;
  style["multijet"]["mpfu"] = kFullTriangleDown;
  style["multijet"]["rho"] = kFullTriangleDown;
  style["multijet_mc"]["mpf1"] = kOpenTriangleUp;
  style["multijet_mc"]["mpfn"] = kOpenTriangleUp;
  style["multijet_mc"]["mpfu"] = kOpenTriangleDown;
  style["multijet_mc"]["rho"] = kOpenTriangleDown;
  style["gamjet"]["mpfchs1"] = kFullSquare;
  style["gamjet"]["ptchs"] = kOpenSquare;
  style["gamjet"]["mpf1"] = kFullTriangleUp;
  style["gamjet"]["mpfn"] = kFullTriangleUp;
  style["gamjet"]["mpfu"] = kFullTriangleDown;
  style["gamjet"]["rho"] = kFullTriangleDown;
  style["gamjet_mc"]["mpf1"] = kOpenTriangleUp;
  style["gamjet_mc"]["mpfn"] = kOpenTriangleUp;
  style["gamjet_mc"]["mpfu"] = kOpenTriangleDown;
  style["gamjet_mc"]["rho"] = kOpenTriangleDown;
  style["zeejet"]["mpfchs1"] = kFullCircle;
  style["zeejet"]["ptchs"] = kOpenCircle;
  style["zeejet"]["mpf1"] = kFullTriangleUp;
  style["zeejet"]["mpfn"] = kFullTriangleUp;
  style["zeejet"]["mpfu"] = kFullTriangleDown;
  style["zeejet"]["rho"] = kFullTriangleDown;
  style["zeejet_mc"]["mpf1"] = kOpenTriangleUp;
  style["zeejet_mc"]["mpfn"] = kOpenTriangleUp;
  style["zeejet_mc"]["mpfu"] = kOpenTriangleDown;
  style["zeejet_mc"]["rho"] = kOpenTriangleDown;
  style["zmmjet"]["mpfchs1"] = kFullStar;
  style["zmmjet"]["ptchs"] = kOpenStar;
  style["zmmjet"]["mpf1"] = kFullTriangleUp;
  style["zmmjet"]["mpfn"] = kFullTriangleUp;
  style["zmmjet"]["mpfu"] = kFullTriangleDown;
  style["zmmjet"]["rho"] = kFullTriangleDown;
  style["zmmjet_mc"]["mpf1"] = kOpenTriangleUp;
  style["zmmjet_mc"]["mpfn"] = kOpenTriangleUp;
  style["zmmjet_mc"]["mpfu"] = kOpenTriangleDown;
  style["zmmjet_mc"]["rho"] = kOpenTriangleDown;
  style["zlljet"]["mpfchs1"] = kFullDiamond;
  style["zlljet"]["ptchs"] = kOpenDiamond;
  style["zlljet"]["mpf1"] = kFullTriangleUp;
  style["zlljet"]["mpfn"] = kFullTriangleUp;
  style["zlljet"]["mpfu"] = kFullTriangleDown;
  style["zlljet"]["rho"] = kFullTriangleDown;
  style["zlljet_mc"]["mpf1"] = kOpenTriangleUp;
  style["zlljet_mc"]["mpfn"] = kOpenTriangleUp;
  style["zlljet_mc"]["mpfu"] = kOpenTriangleDown;
  style["zlljet_mc"]["rho"] = kOpenTriangleDown;
  style["zjet"]["mpfchs1"] = kFullDiamond;
  style["zjet"]["ptchs"] = kOpenDiamond;
  style["zlljet"]["ptchs"] = kOpenDiamond;
  style["zlljet"]["chf"] = kFullCircle;
  style["zlljet"]["nhf"] = kFullDiamond;
  style["zlljet"]["nef"] = kFullSquare;
  style["zlljet"]["cef"] = kFullDiamond;
  style["zlljet"]["muf"] = kFullDiamond;
  style["zlljet_mc"]["chf"] = kOpenCircle;
  style["zlljet_mc"]["nhf"] = kOpenDiamond;
  style["zlljet_mc"]["nef"] = kOpenSquare;
  style["zlljet_mc"]["cef"] = kOpenDiamond;
  style["zlljet_mc"]["muf"] = kOpenDiamond;
  style["zlljet"]["mpf1"] = kFullTriangleUp;
  style["zlljet"]["mpfn"] = kFullTriangleUp;
  style["zlljet"]["mpfu"] = kFullTriangleDown;
  style["zlljet"]["rho"] = kFullTriangleDown;
  style["zlljet_mc"]["mpf1"] = kOpenTriangleUp;
  style["zlljet_mc"]["mpfn"] = kOpenTriangleUp;
  style["zlljet_mc"]["mpfu"] = kOpenTriangleDown;
  style["zlljet_mc"]["rho"] = kOpenTriangleDown;
  style["zjet"]["chf"] = kFullCircle;
  style["zjet"]["nhf"] = kFullDiamond;
  style["zjet"]["nef"] = kFullSquare;
  style["zjet"]["cef"] = kFullDiamond;
  style["zjet"]["muf"] = kFullDiamond;
  style["zjet_mc"]["chf"] = kOpenCircle;
  style["zjet_mc"]["nhf"] = kOpenDiamond;
  style["zjet_mc"]["nef"] = kOpenSquare;
  style["zjet_mc"]["cef"] = kOpenDiamond;
  style["zjet_mc"]["muf"] = kOpenDiamond;
  style["zjet"]["mpf1"] = kFullTriangleUp;
  style["zjet"]["mpfn"] = kFullTriangleUp;
  style["zjet"]["mpfu"] = kFullTriangleDown;
  style["zjet"]["rho"] = kFullTriangleDown;
  style["zjet_mc"]["mpf1"] = kOpenTriangleUp;
  style["zjet_mc"]["mpfn"] = kOpenTriangleUp;
  style["zjet_mc"]["mpfu"] = kOpenTriangleDown;
  style["zjet_mc"]["rho"] = kOpenTriangleDown;
  style["gamjet"]["chf"] = kFullCircle;
  style["gamjet"]["nhf"] = kFullDiamond;
  style["gamjet"]["nef"] = kFullSquare;
  style["gamjet"]["cef"] = kFullDiamond;
  style["gamjet"]["muf"] = kFullDiamond;
  style["gamjet_mc"]["chf"] = kOpenCircle;
  style["gamjet_mc"]["nhf"] = kOpenDiamond;
  style["gamjet_mc"]["nef"] = kOpenSquare;
  style["gamjet_mc"]["cef"] = kOpenDiamond;
  style["gamjet_mc"]["muf"] = kOpenDiamond;

  map<string, int> color;
  color["pfjet"] = kBlack;
  color["pfjet_chf"] = kRed;
  color["pfjet_nhf"] = kGreen+2;
  color["pfjet_nef"] = kBlue;
  color["pfjet_cef"] = kCyan+1;
  color["pfjet_muf"] = kMagenta+1;
  color["dijet"] = kBlack;
  color["incjet"] = kOrange+2;
  color["hadw"] = kGreen+2;
  color["multijet"] = kBlack;
  color["multijet_mpf1"] = kRed;
  color["multijet_mpfn"] = kGreen+2;
  color["multijet_mpfu"] = kBlue;
  color["multijet_rho"] = kBlack;
  color["gamjet"] = kBlue;
  color["gamjet_mpf1"] = kRed;
  color["gamjet_mpfn"] = kGreen+2;
  color["gamjet_mpfu"] = kBlue;
  color["gamjet_rho"] = kBlack;
  color["gamjet_chf"] = kRed;
  color["gamjet_nhf"] = kGreen+2;
  color["gamjet_nef"] = kBlue;
  color["gamjet_cef"] = kCyan+1;
  color["gamjet_muf"] = kMagenta+1;
  color["zeejet"] = kGreen+2;
  color["zeejet_mpf1"] = kRed;
  color["zeejet_mpfn"] = kGreen+2;
  color["zeejet_mpfu"] = kBlue;
  color["zeejet_rho"] = kBlack;
  color["zeejet"] = kMagenta+2;
  color["zmmjet_mpf1"] = kRed;
  color["zmmjet_mpfn"] = kGreen+2;
  color["zmmjet_mpfu"] = kBlue;
  color["zmmjet_rho"] = kBlack;
  color["zlljet"] = kRed;
  color["zlljet_mpf1"] = kRed;
  color["zlljet_mpfn"] = kGreen+2;
  color["zlljet_mpfu"] = kBlue;
  color["zlljet_rho"] = kBlack;
  color["zlljet_chf"] = kRed;
  color["zlljet_nhf"] = kGreen+2;
  color["zlljet_nef"] = kBlue;
  color["zlljet_cef"] = kCyan+1;
  color["zlljet_muf"] = kMagenta+1;
  color["zjet"] = kRed+1;
  color["zjet_muf"] = kMagenta+1;
  color["zjet_mpf1"] = kRed;
  color["zjet_mpfn"] = kGreen+2;
  color["zjet_mpfu"] = kBlue;
  color["zjet_rho"] = kBlack;
  color["zjet_chf"] = kRed;
  color["zjet_nhf"] = kGreen+2;
  color["zjet_nef"] = kBlue;
  color["zjet_cef"] = kCyan+1;
  color["zjet_muf"] = kMagenta+1;
  color["zi"] = kBlack;
  color["zb"] = kRed;
  color["zc"] = kGreen+2;
  color["zq"] = kBlue;
  color["zg"] = kOrange+1;
  color["zn"] = kGray;


  map<double, double> size;
  size[0.10] = 0.6;
  size[0.15] = 0.8;
  size[0.20] = 1.0;
  size[0.30] = 1.2;
  size[0.50] = 1.2; // zjet only
  size[1.00] = 1.2; // zjet only


  // Select which subsets of data to process:
  // datamc x method x sample x etabin x alphacut
  //    3   x    3   x   4    x   6-8  x    4      = 864-1152 graphs
  vector<string> dirs;
  dirs.push_back("mc");
  dirs.push_back("data");
  dirs.push_back("ratio");

  vector<string> types;
  types.push_back("counts");
  types.push_back("crecoil");
  //types.push_back("mpfchs"); // raw MET
  types.push_back("mpfchs1"); // Type-I MET
  types.push_back("ptchs");
  // for pfjet only (activate puf, cef, muf later?)
  types.push_back("chf");
  types.push_back("nef");
  types.push_back("nhf");
  types.push_back("cef");
  types.push_back("muf");
  types.push_back("puf");
  types.push_back("mpf1");
  types.push_back("mpfn");
  types.push_back("mpfu");
  types.push_back("rho");
  // <pT,reco> and <pT,gen> vs ref pT (MC only)
  types.push_back("rjet");
  types.push_back("gjet");

  vector<string> sets;
  //sets.push_back("dijet");
  sets.push_back("incjet");
  sets.push_back("multijet");
  sets.push_back("gamjet");
  sets.push_back("zeejet");
  sets.push_back("zmmjet");
  sets.push_back("zlljet");
  sets.push_back("zjet");
  sets.push_back("hadw");
  sets.push_back("pfjet");
  //
  sets.push_back("zi");
  sets.push_back("zb");
  sets.push_back("zc");
  sets.push_back("zq");
  sets.push_back("zg");
  sets.push_back("zn");
  //
  sets.push_back("zii");
  sets.push_back("zib");
  sets.push_back("zic");
  sets.push_back("ziq");
  if ((isUL16 && !isAPV) && false) {
    sets.push_back("ziu"); //new (in v37, but not yet working)
    sets.push_back("zis"); //new (in v37, but not yet working)
  }
  sets.push_back("zig");
  sets.push_back("zin");
  //
  sets.push_back("zbi");
  sets.push_back("zbb");
  sets.push_back("zbc");
  sets.push_back("zbq");
  sets.push_back("zbg");
  sets.push_back("zbn");
  //
  sets.push_back("zci");
  sets.push_back("zcb");
  sets.push_back("zcc");
  sets.push_back("zcq");
  sets.push_back("zcg");
  sets.push_back("zcn");
  //
  sets.push_back("zqi");
  sets.push_back("zqb");
  sets.push_back("zqc");
  sets.push_back("zqq");
  sets.push_back("zqg");
  sets.push_back("zqn");
  //
  sets.push_back("zgi");
  sets.push_back("zgb");
  sets.push_back("zgc");
  sets.push_back("zgq");
  sets.push_back("zgg");
  sets.push_back("zgn");
  //
  sets.push_back("zni");
  sets.push_back("znb");
  sets.push_back("znc");
  sets.push_back("znq");
  sets.push_back("zng");
  sets.push_back("znn");
  //
  sets.push_back("gi");
  sets.push_back("gb");
  sets.push_back("gc");
  sets.push_back("gq");
  sets.push_back("gg");
  sets.push_back("gn");
  //
  sets.push_back("gii");
  sets.push_back("gib");
  sets.push_back("gic");
  sets.push_back("giq");
  sets.push_back("gig");
  sets.push_back("gin");
  //
  sets.push_back("gbi");
  sets.push_back("gbb");
  sets.push_back("gbc");
  sets.push_back("gbq");
  sets.push_back("gbg");
  sets.push_back("gbn");
  //
  sets.push_back("gci");
  sets.push_back("gcb");
  sets.push_back("gcc");
  sets.push_back("gcq");
  sets.push_back("gcg");
  sets.push_back("gcn");
  //
  sets.push_back("gqi");
  sets.push_back("gqb");
  sets.push_back("gqc");
  sets.push_back("gqq");
  sets.push_back("gqg");
  sets.push_back("gqn");
  //
  sets.push_back("ggi");
  sets.push_back("ggb");
  sets.push_back("ggc");
  sets.push_back("ggq");
  sets.push_back("ggg");
  sets.push_back("ggn");
  //
  sets.push_back("gni");
  sets.push_back("gnb");
  sets.push_back("gnc");
  sets.push_back("gnq");
  sets.push_back("gng");
  sets.push_back("gnn");

  if (isLowPU) {
    sets.clear();
    sets.push_back("incjet");
    sets.push_back("zjet");
    sets.push_back("pfjet"); //Mar1

    sets.push_back("zi");
    sets.push_back("zb");
    sets.push_back("zc");
    sets.push_back("zq");
    sets.push_back("zg");
    sets.push_back("zn");
    //
    sets.push_back("zii");
    sets.push_back("zib");
    sets.push_back("zic");
    sets.push_back("ziq");
    sets.push_back("ziu"); //new in v38
    sets.push_back("zis"); //new in v38
    sets.push_back("zig");
    sets.push_back("zin");
    //
    sets.push_back("zbi");
    sets.push_back("zbb");
    sets.push_back("zbc");
    sets.push_back("zbq");
    sets.push_back("zbg");
    sets.push_back("zbn");
    //
    sets.push_back("zci");
    sets.push_back("zcb");
    sets.push_back("zcc");
    sets.push_back("zcq");
    sets.push_back("zcg");
    sets.push_back("zcn");
    //
    sets.push_back("zqi");
    sets.push_back("zqb");
    sets.push_back("zqc");
    sets.push_back("zqq");
    sets.push_back("zqg");
    sets.push_back("zqn");
    //
    sets.push_back("zgi");
    sets.push_back("zgb");
    sets.push_back("zgc");
    sets.push_back("zgq");
    sets.push_back("zgg");
    sets.push_back("zgn");
    //
    sets.push_back("zni");
    sets.push_back("znb");
    sets.push_back("znc");
    sets.push_back("znq");
    sets.push_back("zng");
    sets.push_back("znn");
  }

  vector<pair<double,double> > etas;
  etas.push_back(make_pair<double,double>(0,1.305));
  etas.push_back(make_pair<double,double>(0,2.500));
  // Narrow eta bins for L2Res
  
  /*
  etas.push_back(make_pair<double,double>(0.000,0.261)); 
  etas.push_back(make_pair<double,double>(0.261,0.522)); 
  // PATCH!! V15 gamjet combines 0-0.5 bins
  etas.push_back(make_pair<double,double>(0.522,0.783)); 
  etas.push_back(make_pair<double,double>(0.783,1.044)); 
  etas.push_back(make_pair<double,double>(1.044,1.305)); 
  // PATCH!! V15 gamjet has 0.8-1.1? 
  etas.push_back(make_pair<double,double>(1.305,1.479)); 
  etas.push_back(make_pair<double,double>(1.479,1.653)); 
  // PATCH!! V15 gamjet has 1.3-1.7? 
  etas.push_back(make_pair<double,double>(1.653,1.930)); 
  etas.push_back(make_pair<double,double>(1.930,2.172)); 
  etas.push_back(make_pair<double,double>(2.172,2.322)); 
  etas.push_back(make_pair<double,double>(2.322,2.500)); 
  etas.push_back(make_pair<double,double>(2.500,2.650)); 
  etas.push_back(make_pair<double,double>(2.650,2.853)); 
  etas.push_back(make_pair<double,double>(2.853,2.964)); 
  etas.push_back(make_pair<double,double>(2.964,3.139)); 
  etas.push_back(make_pair<double,double>(3.139,3.489)); 
  etas.push_back(make_pair<double,double>(3.489,3.839)); 
  etas.push_back(make_pair<double,double>(3.839,5.191));

  // Wide eta bins for L2L3Res closure
  etas.push_back(make_pair<double,double>(1.305,1.93));
  etas.push_back(make_pair<double,double>(1.93,2.5));
  etas.push_back(make_pair<double,double>(2.5,2.964));
  etas.push_back(make_pair<double,double>(2.964,3.2));
  etas.push_back(make_pair<double,double>(3.2,5.191));
  */


  vector<double> alphas;
  alphas.push_back(0.10);
  alphas.push_back(0.15);
  alphas.push_back(0.20); //  => patch because of a20 missing for Z+jet
  alphas.push_back(0.30);
  alphas.push_back(0.50); // for zjet only
  alphas.push_back(1.00); // for zjet only

  if (isLowPU) {
    alphas.clear();
    alphas.push_back(0.30); // for pfjet
    alphas.push_back(1.00); // for zjet
  }

  ///////////////////////////////////////////
  // Rename selected graphs and store them //
  ///////////////////////////////////////////

  map<string, map<string, map<string, map<int, map<int, TGraphErrors*> > > > > grs;
  map<string, map<string, map<int, map<int, TH1D*> > > > counts;

  // Loop over data, MC, ratio
  for (unsigned int idir = 0; idir != dirs.size(); ++idir) {

    string d = dirs[idir];
    const char *dd = d.c_str();

    fout->mkdir(dd);
    assert(fout->cd(dd));
    //TDirectory *dout0 = gDirectory; // broke in ROOT 6.04/08
    TDirectory *dout0 = fout->GetDirectory(dd); assert(dout0);

    // Loop over eta bins
    for (unsigned int ieta = 0; ieta != etas.size(); ++ieta) {
      
      double eta1 = etas[ieta].first;
      double eta2 = etas[ieta].second;
      const char *dd0 = Form("eta%02.0f-%02.0f",eta1*10.,eta2*10.);
      cout << dd0 << endl << flush;
      
      dout0->mkdir(dd0);
      assert(dout0->cd(dd0));
      //TDirectory *dout = gDirectory; // broke in ROOT 6.04/08
      TDirectory *dout = dout0->GetDirectory(dd0); assert(dout);
      dout->mkdir("orig");
      dout->cd();

      // Save background histogram for easy drawing from root file
      TH1D *h = new TH1D("h",Form(";p_{T} (GeV);Response (%s)",dd),
			 int(2000-10),10,2000);
      h->GetXaxis()->SetMoreLogLabels();
      h->GetXaxis()->SetNoExponent();
      h->SetMinimum(0.50);
      h->SetMaximum(1.50);
      h->Write();
      
      // Loop over methods (MPF, pT balance)
      for (unsigned int itype = 0; itype != types.size(); ++itype) {
	
	string t = types[itype];
	const char* tt = t.c_str();

	// Loop over samples (Z+jet, gamma+jet, dijet)
	for (unsigned int iset = 0; iset != sets.size(); ++iset) {
	
	  string s = sets[iset];
	  const char* ss = s.c_str();
	  string sp = (rename[s]["parent"]!=0 ? rename[s]["parent"] : "");
	  
	  TFile *f = files[s];
	  if (!f && sp=="gamjet") f = files[sp];
	  // pfjet has two files
	  if (!f && s=="pfjet") f = files[s+"_"+d];
	  if (s=="pfjet" && pfMode=="none") continue;

	  // Take pT and MPF from different files for gamma+jet (or not)
	  //if (s=="gamjet" && f==0) {
	  //if (t=="mpfchs" || t=="mpfchs1") f = fp;//fp1;
	  //if (t=="ptchs") f = fp;//fp2;
	  //}
	  if (!(f || s=="zlljet" || (s=="pfjet" && d=="ratio"))) 
	    cout << "File missing for " << s << endl << flush;
	  assert(f || s=="zlljet" || (s=="pfjet" && d=="ratio"));

	  if (eta1==0 && eta2==2.5 &&
	      !(s=="zjet" || sp=="zjet" || s=="whad"))
	    continue;

	  // C_recoil is only meant for multijets
	  if (t=="crecoil" && s!="multijet") continue;

	  // MPF decomposition only for Z+jet for now (+multijet, +photon+jet)
	  if ((t=="mpf1" || t=="mpfn" || t=="mpfu" || t=="rho") &&
	      !(((s=="zeejet" || s=="zmmjet" || s=="zlljet" || s=="zjet" ||
		  sp=="zjet" || 
		  s=="multijet") && isUL18) ||
		((s=="zeejet" || s=="zmmjet" || s=="zlljet" || s=="zjet" ||
		  s=="multijet") && (isUL17 || isUL16 || isRun2)) ||
		(s=="gamjet" && (epoch=="2016EF" || epoch=="2016GH" ||
				 epoch=="2016BCD" || epoch=="2016BCDEF" ||
				 epoch=="2018ABCD" || epoch=="2017BCDEF" ||
				 epoch=="Run2Test")
		 ) || //&& t!="rho") ||
		(sp=="gamjet") ||
		(s=="zjet" || sp=="zjet") ||
		(s=="hadw" && t=="rho") ||
		(s=="pfjet" && t=="rho")
		))
	    continue;
	  //if (t=="rho" && s=="zjet") continue; // still missing for zjet
	  if (t=="rho" && sp=="zjet") continue;
	  //if (t=="rho" && s=="gamjet") continue; // DP_2021
	  if (t=="rho" && sp=="gamjet") continue;
	  if (t=="rho" && s=="multijet"  // and for multijet
	      && !(CorLevel=="L1L2L3Res" && isUL16)) continue;
	  if (t=="rho" && s=="pfjet" && isLowPU) continue;

	  // fractions are only for pfjet, and only fractions are for pfjet
	  // (now adding fractions to zjet)
	  bool isfrac = (t=="chf"||t=="nef"||t=="nhf"||
			 t=="cef"||t=="muf"||t=="puf");
	  if (isfrac && s!="pfjet" && s!="zjet" && s!="gamjet" && // v19 latest
	      //if (isfrac && s!="pfjet" && s!="zjet" && // DP_2021
	      s!="zeejet" && s!="zmmjet" && s!="zlljet") continue;
	  if (!(isfrac || t=="rho") && s=="pfjet") continue;
	  if (t=="puf" && s!="pfjet") continue;

	  bool isflavor = 
	    (s=="zi"  || s=="zb"  || s=="zc"  || s=="zq"  || s=="zg"||s=="zn")||
	    (s=="gi"  || s=="gb"  || s=="gc"  || s=="gq"  || s=="gg"||s=="gn");
	  bool isflavormc = 
	    (s=="zii" || s=="zbi" || s=="zci" || s=="zqi"||s=="zgi"||s=="zni"||
	     s=="zib" || s=="zbb" || s=="zcb" || s=="zqb"||s=="zgb"||s=="znb"||
	     s=="zic" || s=="zbc" || s=="zcc" || s=="zqc"||s=="zgc"||s=="znc"||
	     s=="ziq" || s=="zbq" || s=="zcq" || s=="zqq"||s=="zgq"||s=="znq"||
	     s=="zig" || s=="zbg" || s=="zcg" || s=="zqg"||s=="zgg"||s=="zng"||
	     s=="zin" || s=="zbn" || s=="zcn" || s=="zqn"||s=="zgn"||s=="znn"||
	     s=="ziu" || s=="zis")||
	    (s=="gii" || s=="gbi" || s=="gci" || s=="gqi"||s=="ggi"||s=="gni"||
	     s=="gib" || s=="gbb" || s=="gcb" || s=="gqb"||s=="ggb"||s=="gnb"||
	     s=="gic" || s=="gbc" || s=="gcc" || s=="gqc"||s=="ggc"||s=="gnc"||
	     s=="giq" || s=="gbq" || s=="gcq" || s=="gqq"||s=="ggq"||s=="gnq"||
	     s=="gig" || s=="gbg" || s=="gcg" || s=="gqg"||s=="ggg"||s=="gng"||
	     s=="gin" || s=="gbn" || s=="gcn" || s=="gqn"||s=="ggn"||s=="gnn");
	  if (isflavormc && d!="mc") continue;

	  // true responses
	  if ((t=="rjet" || t=="gjet") && sp!="gamjet") continue;

	  // Some inputs missing for UL16 and UL17
	  //if ((isUL16 || isUL17 || isRun2) && 
	  //  (s=="zn" || s=="zin" || s=="zni" || s=="znn" ||
	  //   s=="zbn" || s=="zcn" || s=="zqn" || s=="zgn" || 
	  //   s=="znb" || s=="znc" || s=="znq" || s=="zng"))// || 
	  //  continue;

	  // MPF composition
	  bool ismpfc = (t=="mpf1" || t=="mpfn" || t=="mpfu" || t=="rho");

	  for (unsigned int ialpha = 0; ialpha != alphas.size(); ++ialpha) {

	    double alpha = alphas[ialpha];

	    eta1 = etas[ieta].first; // reset to avoid trouble with below
	    eta2 = etas[ieta].second; // reset to avoid trouble with below
	    
	    // Take wide and narrow bins from different file for dijet and Z+jet
	    bool narrowBin = ((fabs(eta2-eta1)<0.4 && 
	    		       !(fabs(eta1-3.0)<0.05 && fabs(eta2-3.2)<0.05)) ||
	    		      fabs(eta1)>3.8);
	    if (narrowBin) {
	    if (s=="dijet")  f = fdj2;
	    }
	    assert(f || s=="zlljet" || (s=="pfjet"&&d=="ratio"));

	    // Only store a30,a100 for zb,zc,zq,zg for now to not bloat the file
	    //if ((isflavor||isflavormc) && !(alpha==0.30 || alpha==1.00))
	    if ((isflavor||isflavormc) && alpha!=1.00)
	      continue;

	    if (alpha!=1.00 && isUL18 && isNewName && (s=="zjet"||sp=="zjet"))
	      continue;
	    if (alpha!=1.00 && isUL16 && !isAPV && (s=="zjet" || sp=="zjet"))
	      continue;
	    if (alpha!=1.00 && isLowPU && (s=="zjet" || sp=="zjet"))
	      continue;

	    if (alpha!=1.00 && sp=="gamjet") continue;
	    if (alpha>0.35 && (s!="zjet" && sp!="zjet" && s!="zeejet" && s!="zmmjet" && s!="zlljet" && s!="gamjet" && sp!="gamjet" && s!="incjet")) continue;
	    //if (alpha>0.35 && s=="gamjet" && isUL17) continue;
	    //if (alpha>0.35 && s=="gamjet" && isRun2) continue;

            if (t=="counts" && s!="zmmjet" && s!="zeejet" && s!="gamjet" &&
		s!="zjet" && sp!="zjet" && sp!="gamjet" &&
		s!="hadw" && s!="multijet")
              continue; // counts available only for z+jet and gamjet, so far

	    if (s=="pfjet" && !(fabs(eta1-0)<0.1 && fabs(eta2-1.3)<0.1
				&& fabs(alpha-0.30)<0.01))
	      continue; // barrel only, no specific alpha
	    if (s=="incjet" && !(fabs(eta1-0)<0.1 && fabs(eta2-1.3)<0.1
				 && fabs(alpha-1.00)<0.01))
	      //&& fabs(alpha-0.30)<0.01)) // && t=="mpfchs1"))
	      continue; // barrel only, no specific pt or alpha
	    if (s=="hadw" && !(fabs(eta1-0)<0.1 && fabs(eta2-1.3)<0.1
			       && fabs(alpha-0.30)<0.01)) // && t=="mpfchs1"))
	      continue; // barrel only, no specific pt or alpha
	    if (s=="multijet" && (!(fabs(eta1-0)<0.1 && fabs(eta2-1.3)<0.1)
				  || fabs(alpha-0.10)<0.01))
	      continue; // only barrel for multijet balance, pT=(15),(20),30
	    if (s=="gamjet"  && fabs(eta1-3.2)<0.1) { eta1=3.0; eta2=3.2; }

	    //if ((s=="ziu"||s=="zis") && t!="counts" && !isLowPU)
	    //continue; // v37 botch

	    int ieta1 = int(10.*eta1+0.5);
	    int ieta2 = int(10.*eta2+0.5);

	    // Reconstruct naming scheme used in each of the files
	    // If non-conventional naming schemes, patch here
	    const char *c(0);
	    if (s=="pfjet") {
	      if (t=="rho")
		c = Form("Standard/Eta_%02.1f-%02.1f/p%s",
			 eta1, eta2, tt);
	      else
		c = Form("Standard/Eta_%02.1f-%02.1f/p%s%s",
			 eta1, eta2, tt, pfMode=="tp" ? "tp" : ""); 
	    }
	    if (s=="dijet") {
	      c = Form("%s/eta_%02.0f_%02.0f/%s_%s_a%1.0f",
                       dd, 10.001*eta1, 10.001*eta2, rename[s][t], ss, 100.*alpha); 
	    } // dijet
	    if (s=="incjet") {
	      //if (isUL18) {
	      //c = Form(rename[s][t],fij_eras[epoch]);
	      //if (epoch=="2018ABCD" && CorLevel=="L1L2Res")
	      //  c = Form(rename[s][t],"");
	      //}
	      //else 
	      if (isUL17 || isUL16 || isUL18 || isRun2 || isLowPU) {
		c = Form(rename[s][t],fij_eras[epoch]);
	      }
	      else
		assert(false);
	    }
	    if (s=="hadw") {
	      if (CorLevel=="L1L2L3Res") {
		//if (isUL18 || isUL17)
		//if (isUL18)
		//c = Form("%s_%s_%s_fitprob02_L1L2L3",dd,rename[s][t],ss);
		//else 
		if (isUL16 || isUL17 || isUL18)
		  c = Form("%s_%s_%s_fitprob001_L1L2L3",dd,rename[s][t],ss);
		else
		  assert(isRun2);
	      }
	      if (CorLevel=="L1L2Res") {
		if (isUL18 || isUL17)
		  c = Form("%s_%s_%s_fitprob02_L1L2L3",dd,rename[s][t],ss);
		else if (isUL16)
		  c = Form("%s_%s_%s_fitprob001_L1L2L3",dd,rename[s][t],ss);
		else
		  assert(false);
	      }
	    }
	    if (s=="multijet") {
	      c = Form("%s/Pt%1.0f/%s",rename[s][d], 100.*alpha,
		       rename[s][(d+t)]);
	    } // multijet
	    if (s=="gamjet" && t=="counts") {
	      c = Form("%s%s_a%1.0f_eta%02.0f_%02.0f_%s",
		       rename[s]["mpfchs1"],
		       d=="ratio" ? rename[s]["data"] : rename[s][d],
		       100.*alpha, 10.*eta1, 10.*eta2, rename[s][t]);
	    }
	    //else if (s=="gamjet" && t=="rho") {
	    //c = Form("%s_%s",rename[s][t],dd);
	    //}
	    else if (s=="gamjet" && (ismpfc||t=="mpfchs1"||t=="ptchs") &&
		     d=="ratio") { // patch missing
	      c = Form("%s%s_a%1.0f_eta%02.0f_%02.0f",
		       rename[s][t], rename[s]["data"], // ratio->data
		       100.*alpha, 10.*eta1, 10.*eta2);
	    } // gamjet
	    else if (s=="gamjet" && isfrac) { // new PF composition
	      c = Form("pf/p%s_%s",tt,dd);
	    }
	    else if (s=="gamjet") {
	      c = Form("%s%s_a%1.0f_eta%02.0f_%02.0f%s",
		       rename[s][t], rename[s][d],
		       100.*alpha, 10.*eta1, 10.*eta2,
		       (d!="data" ? sPSWgtG.c_str() : ""));
	    } // gamjet
	    if (sp=="gamjet") {
	      c = Form("flavor/%s_%s_%s",tt,ss,dd);
	    }
	    if (s=="zmmjet" || s=="zeejet") {
	      if (CorLevel=="L1L2Res") {
		c = Form("%s_%s_a%1.0f_eta_%02.0f_%02.0f_L1L2Res",
			 rename[s][d], rename[s][t], 100.*alpha,
			 10.01*eta1, 10.01*eta2);
		assert(false); // make sure not going here for V7 closure
	      }
	      else if (CorLevel=="L1L2L3Res") {
		c = Form("%s_%s_a%1.0f_eta_%02.0f_%02.0f_L1L2L3Res",
			 rename[s][d], rename[s][t], 100.*alpha,
			 10.01*eta1, 10.01*eta2);
	      }
	      else
		assert(false);
	    } // Zll+jet
	    if (s=="zjet") {
	      if (zjetMode=="zpt")
		c = Form("%s/eta_%02.0f_%02.0f/%s%s_zmmjet_a%1.0f",
			 rename[s][d],10*eta1,10*eta2,
			 ((d=="data"||d=="ratio"||t=="counts") ? "" : 
			  sPSWgtZ.c_str()),
			 rename[s][t],100.*alpha);
	      else if (zjetMode=="jetpt")
		c = Form("%s/eta_%02.0f_%02.0f/%s%s_jetpt_a%1.0f",
			 rename[s][d],10*eta1,10*eta2,
			 ((d=="data"||d=="ratio"||t=="counts") ? "" : 
			  sPSWgtZ.c_str()),
			 rename[s][t],100.*alpha);
	      else if (zjetMode=="ptave")
		c = Form("%s/eta_%02.0f_%02.0f/%s%s_ptave_a%1.0f",
			 rename[s][d],10*eta1,10*eta2,
			 ((d=="data"||d=="ratio"||t=="counts") ? "" : 
			  sPSWgtZ.c_str()),
			 rename[s][t],100.*alpha);
	      else
		assert(false);
	      //
	      if (isfrac || t=="rho") {
		if (t=="rho")
		  c = Form("%s/eta_%02.0f_%02.0f/h_Zpt_%s_alpha%1.0f",rename[s][d],10*eta1,10*eta2,rename[s][t],100.*alpha); // v44
		else if (zjetMode=="zpt")
		  c = Form("%s/eta_%02.0f_%02.0f/h_Zpt_%s_alpha%1.0f",rename[s][d],10*eta1,10*eta2,rename[s][t],100.*alpha);
		//c = Form("%s/eta_%02.0f_%02.0f/h_Zpt_%s_alpha%1.0f_eta_%02.0f_%02.0f",rename[s][d],10*eta1,10*eta2,rename[s][t],100.*alpha,10*eta1,10*eta2); // v43
		else if (zjetMode=="jetpt")
		  c = Form("%s/eta_%02.0f_%02.0f/h_JetPt_%s_alpha%1.0f",rename[s][d],10*eta1,10*eta2,rename[s][t],100.*alpha);
		//c = Form("%s/eta_%02.0f_%02.0f/h_JetPt_%s_alpha%1.0f_eta_%02.0f_%02.0f",rename[s][d],10*eta1,10*eta2,rename[s][t],100.*alpha,10*eta1,10*eta2); // v43
		else if (zjetMode=="ptave")
		  c = Form("%s/eta_%02.0f_%02.0f/h_PtAve_%s_alpha%1.0f",rename[s][d],10*eta1,10*eta2,rename[s][t],100.*alpha);
		//assert(false);
		else
		  assert(false);
	      }
	    } // zjet
	    if (sp=="zjet") {
	      c = Form("%s/eta_%02.0f_%02.0f%s/%s%s_zmmjet%s_a%1.0f",
	    	       rename["zjet"][d],10*eta1,10*eta2,
		       (d=="mc" ? rename[s]["gen"] : ""),
		       ((d=="data"||d=="ratio"||t=="counts") ? "" : 
			sPSWgtZ.c_str()),
		       rename["zjet"][t],rename[s]["tag"],100.*alpha);
	      if (isfrac) {
		assert(false);
	      }
	    } // Z+jet (Z+b)

	    // Run 2: overwrite any previous object name c, except incjet
	    // Inputs created with recombine.C
	    if (isRun2 && s!="incjet")
	      c = Form("%s/eta%02d-%02d/orig/%s_%s_a%1.0f",
		       dd,ieta1,ieta2,tt,ss,100.*alpha);
	    if (isRun2 && t=="counts")
	      c = Form("%s/eta%02d-%02d/%s_%s_a%1.0f",
		       dd,ieta1,ieta2,tt,ss,100.*alpha);


	    TObject *obj = ((s=="zlljet"||(s=="pfjet"&&d=="ratio")
			     //||(s=="zjet"&&d=="ratio")
			     ) ? 0 :
			    f->Get(c));
	    if (!obj && !(s=="zlljet" || (s=="pfjet"&&d=="ratio")
			  //|| (s=="zjet"&&d=="ratio")
			  )) {
	      cout << "Graph " << c << " not found for "
		   << s << " " << t << " " << d << "!" <<endl << flush;
	      cout << "File: " << f->GetName() << endl << flush;
	      cout << "Eta " << eta1 << " - " << eta2 <<  endl << flush;
	      if (narrowBin) cout << "Narrow bin" << endl << flush;
	      else           cout << "Wide bin" << endl << flush;
	    }
	    
	    // Merge Zee+jet and Zmumu+jet
	    if (s=="zlljet") {

	      TGraphErrors *gee =  grs[d][t]["zeejet"][ieta][ialpha];
	      TGraphErrors *gmm =  grs[d][t]["zmmjet"][ieta][ialpha];
	      assert(gee);
	      assert(gmm);

	      TGraphErrors *g = tools::mergeGraphs(gee, gmm);
	      
	      // clean out low pT points
	      for (int i = g->GetN()-1; i != -1; --i) {
		if (g->GetX()[i]<fzllptmin) g->RemovePoint(i);
	      }

	      obj = (TObject*)g;
	    }
	    // Calculate data-MC diff
	    if (s=="pfjet" && d=="ratio") {
	      TGraphErrors *gd = grs["data"][t][s][ieta][ialpha];
	      TGraphErrors *gm = grs["mc"][t][s][ieta][ialpha];
	      assert(gd);
	      assert(gm);

	      TGraphErrors *g = tools::diffGraphs(gd, gm);
	      obj = (TObject*)g;
	    }
	    // Calculate data/MC ratio
	    if ((s=="zjet" || sp=="zjet") &&
		d=="ratio" && (t=="mpfchs1"||t=="ptchs")) {
	      TGraphErrors *gd = grs["data"][t][s][ieta][ialpha];
	      TGraphErrors *gm = grs["mc"][t][s][ieta][ialpha];
	      assert(gd);
	      assert(gm);

	      TGraphErrors *g = tools::ratioGraphs(gd, gm);
	      obj = (TObject*)g;
	    }
	    // Calculate data/MC ratio or difference
	    if (s=="gamjet" && d=="ratio" &&
		(ismpfc || t=="mpfchs1" | t=="ptchs")) {
	      TGraphErrors *gd = grs["data"][t][s][ieta][ialpha];
	      TGraphErrors *gm = grs["mc"][t][s][ieta][ialpha];
	      assert(gd);
	      assert(gm);

	      if (t=="mpf1" || t=="mpfchs1" || t=="ptchs") {
		TGraphErrors *g = tools::ratioGraphs(gd, gm);
		obj = (TObject*)g;
	      }
	      else {
		TGraphErrors *g = tools::diffGraphs(gd, gm);
		obj = (TObject*)g;
	      }
	    }
	    // Calculate data/MC ratio or difference for gamjet-derivates
	    if (sp=="gamjet" && d=="ratio" && (ismpfc || t=="mpfchs1")) {
	      TGraphErrors *gd = grs["data"][t][s][ieta][ialpha];
	      TGraphErrors *gm = grs["mc"][t][s][ieta][ialpha];
	      assert(gd);
	      assert(gm);

	      if (t=="mpf1" || t=="mpfchs1") {
		TGraphErrors *g = tools::ratioGraphs(gd, gm);
		obj = (TObject*)g;
	      }
	      else {
		TGraphErrors *g = tools::diffGraphs(gd, gm);
		obj = (TObject*)g;
	      }
	    }
	    // Calculate data-MC diff
	    // (for zlljet, replace original ratio)
	    if ((s=="zjet" && d=="ratio" && isfrac) ||
		((s=="zeejet" || s=="zmmjet" || s=="zlljet") &&
		 d=="ratio" && isfrac)) {
	      TGraphErrors *gd = grs["data"][t][s][ieta][ialpha];
	      TGraphErrors *gm = grs["mc"][t][s][ieta][ialpha];
	      assert(gd);
	      assert(gm);

	      TGraphErrors *g = tools::diffGraphs(gd, gm);
	      obj = (TObject*)g;
	    }

	    if (!obj) {
	      cout << "Missing " << c << endl << flush;
	    }
	    assert(obj);

	    // write out counts to jecdata.root (as TH1F)
            if (t=="counts" && (s=="zmmjet" || s=="zeejet" || s=="gamjet" ||
				s=="zjet" || sp=="zjet" || s=="hadw" ||
				s=="multijet" || sp=="gamjet") ){
	      // Multijet counts are retrieved from TH2D
	      if (s=="multijet" && !(isUL16 || isRun2)) {
		assert(obj->InheritsFrom("TH2D"));
		dout->cd();
		obj = ((TH2D*)obj)->RebinX(2);
		obj = ((TH2D*)obj)->ProjectionX(Form("%s_%s_a%1.0f",tt,ss,
						     100.*alphas[ialpha]));
	      }
	      
              assert(obj->InheritsFrom("TH1D") ||obj->InheritsFrom("TH1F"));
	      dout->cd();
	      TH1D *h = (TH1D*)obj;
	      // 20200326 Below doesn't work for some reason?
	      //if (obj->InheritsFrom("TH1F")) h = (TH1D*)obj;
	      //else if (obj->InheritsFrom("TH1D")){
	      //TH1D *temp = (TH1D*) obj;
	      //temp->Copy(*h);// copy contents to different type histogram
              //}
              h->SetName(Form("%s_%s_a%1.0f",tt,ss,100.*alphas[ialpha]));
              //if (obj->InheritsFrom("TH1F"))       cout << Form("%s_%s_a%1.0f is TH1Float and has %f entries, %f effective entries, and integral %f",tt,ss,100.*alphas[ialpha],h->GetEntries(),h->GetEffectiveEntries(),h->Integral()) << endl << flush;
              //else if (obj->InheritsFrom("TH1D"))       cout << Form("%s_%s_a%1.0f is TH1Double and has %f entries, %f effective entries, and integral %f",tt,ss,100.*alphas[ialpha],h->GetEntries(),h->GetEffectiveEntries(),h->Integral()) << endl << flush;

              h->Write();
              counts[d][s][ieta][ialpha] = h;
              continue; // counts reliably available only for z+jet, so far. For gamma+jet MC events contain weights
            }

	    // If data stored in TH2D instead of TGrapherrors, patch here
	    // Patch for zjet that has PF fractions in TH2D (v24)
	    if (obj->InheritsFrom("TH2")) {
	      obj = new TGraphErrors(((TH2D*)obj)->ProfileX()->ProjectionX());
	    }

	    // If data stored in TH1D instead of TGraphErrors, patch here
	    // Patch for dijet file that has TH1D's instead of graphs
	    if (obj->InheritsFrom("TH1")) {
	      obj = new TGraphErrors((TH1D*)obj);
	    }

	    // If data stored in TProfile instead of TGraphErrors, patch here
	    // Patch for pfjet file that has TProfiles's instead of graphs
	    // Patch for zjet file that has TProfiles instead of graphs
	    if (obj->InheritsFrom("TProfile")) {
	      //obj = new TGraphErrors((TProfile*)obj);
	      obj = new TGraphErrors(((TProfile*)obj)->ProjectionX());
	    }

	    assert(obj->InheritsFrom("TGraphErrors"));
	    TGraphErrors *g = (TGraphErrors*)obj;
	    
	    // Clean out empty points from TH1D->TGraphErrors conversion
	    for (int i = g->GetN()-1; i != -1; --i) {
              assert(i<=g->GetN()-1);
	      // Clean out spurious empty pooints
	      if (g->GetY()[i]==0 && g->GetEY()[i]==0) g->RemovePoint(i);
	    } // for i

	    // Make a copy of raw data before cleaning cuts and corrections
	    // for documenting it
	    TGraphErrors *g_orig = (TGraphErrors*)g->Clone();
	    g_orig->SetName(Form("%s_%s_a%1.0f",tt,ss,100.*alphas[ialpha]));
	    
	    // Select stable range and good statistics points for global fit
	    for (int i = g->GetN()-1; i != -1; --i) {
              assert(i<=g->GetN()-1);
              // remove points with large error (more than 0.2 right now)
              if (g->GetEY()[i]>0.2 && !isflavormc)  g->RemovePoint(i);
	      // Clean out point outside good ranges
	      else if ((s=="gamjet"||sp=="gamjet") && (t=="mpfchs1"||ismpfc) &&
		       (g->GetX()[i]<fpmpfptmin || g->GetX()[i]>fpmpfptmax))
		g->RemovePoint(i);
	      else if ((s=="gamjet"||sp=="gamjet") && t=="ptchs" &&
		       (g->GetX()[i]<fpbalptmin || g->GetX()[i]>fpbalptmax))
		g->RemovePoint(i);
	      else if ((s=="gamjet") && isfrac &&
		       (g->GetX()[i]<fppfgptmin || g->GetX()[i]>fppfgptmax))
		g->RemovePoint(i);
	      else if (s=="dijet" && t=="mpfchs1" &&
		       (g->GetX()[i]<fdijetmpfptmin || g->GetX()[i]>fdijetptmax))
		g->RemovePoint(i);
	      else if (s=="dijet" && t=="ptchs" &&
		       (g->GetX()[i]<fdijetbalptmin || g->GetX()[i]>fdijetptmax))
		g->RemovePoint(i);
	      else if (s=="pfjet" &&
		       (g->GetX()[i]<fpfjetptmin || g->GetX()[i]>fpfjetptmax))
		g->RemovePoint(i);
	      else if (s=="incjet" &&
		       (g->GetX()[i]<fincjetptmin || g->GetX()[i]>fincjetptmax))
		g->RemovePoint(i);
	      else if (s=="hadw" && t=="mpfchs1" &&
		       (g->GetX()[i]<fhadwptamin || g->GetX()[i]>fhadwptamax))
		g->RemovePoint(i);
	      else if (s=="hadw" && t=="ptchs" &&
		       (g->GetX()[i]<fhadwptbmin || g->GetX()[i]>fhadwptbmax))
		g->RemovePoint(i);
	      else if (s=="multijet" &&
		       (g->GetX()[i]<fmultijetptmin ||
			g->GetX()[i]>fmultijetptmax ||
			((epoch!="2018ABCD" && epoch!="2017BCDEF" &&
			  epoch!="2016BCDEF" && epoch!="2016GH") &&
			 g->GetX()[i]>fmultijetptmax2) ||
			(t=="ptchs" && (g->GetX()[i]<fmultijetmjbptmin ||
					g->GetX()[i]>fmultijetmjbptmax)) ||
			(t=="mpfchs1" && g->GetX()[i]<fmultijetmpfptmin)))
		g->RemovePoint(i);
	      else if (s=="zeejet" && 
		       (g->GetX()[i]<fzeeptmin || g->GetX()[i]>fzeeptmax))
		g->RemovePoint(i);
	      else if (s=="zmmjet" && 
		       (g->GetX()[i]<fzmmptmin || g->GetX()[i]>fzmmptmax))
		g->RemovePoint(i);
	      else if (s=="zjet" &&
		       (g->GetX()[i]<fzptmin || g->GetX()[i]>fzptmax))
		g->RemovePoint(i);
	      else if (sp=="zjet" && 
		       (g->GetX()[i]<fzptmin || g->GetX()[i]>fzbptmax))
		g->RemovePoint(i);
	      else if ((s=="zmmjet" || s=="zeejet") && t=="mpfchs1" &&
		       (g->GetX()[i]<fzllmpfptmin || g->GetX()[i]>fzllmpfptmax))
		g->RemovePoint(i);
	      else if (s=="zjet" && t=="mpfchs1" &&
		       (g->GetX()[i]<fzmpfptmin || g->GetX()[i]>fzmpfptmax))
		g->RemovePoint(i);
	      else if ((s=="zmmjet" || s=="zeejet") && t=="ptchs" &&
		       (g->GetX()[i]<fzllbalptmin || g->GetX()[i]>fzllbalptmax))
		g->RemovePoint(i);
	      else if (s=="zjet" && t=="ptchs" &&
		       (g->GetX()[i]<fzbalptmin || g->GetX()[i]>fzbalptmax))
		g->RemovePoint(i);
	      else if ((s=="zmmjet" || s=="zeejet") && isfrac &&
		       (g->GetX()[i]<fzllpfzptmin || g->GetX()[i]>fzllpfzptmax))
		g->RemovePoint(i);
	      else if ((s=="zjet") && isfrac &&
		       (g->GetX()[i]<fzpfzptmin || g->GetX()[i]>fzpfzptmax))
		g->RemovePoint(i);
              else if (s=="zmmjet" || s=="zeejet" || s=="gamjet" || s=="zjet"){ // patch: clean away points with low statistics based on event counts histograms, currently Z+jet
                assert(counts[d][s][ieta][ialpha]);
                double pt = g->GetX()[i];
                int ipt = counts[d][s][ieta][ialpha]->FindBin(pt);
		double nentries = counts[d][s][ieta][ialpha]->GetBinContent(ipt);
                if (d=="ratio"){
                  assert(counts["mc"][s][ieta][ialpha]);
                  assert(counts["data"][s][ieta][ialpha]);
		  int ipt_dt = counts["data"][s][ieta][ialpha]->FindBin(pt);
		  int ipt_mc = counts["mc"][s][ieta][ialpha]->FindBin(pt);

                  if(counts["mc"][s][ieta][ialpha]->GetBinContent(ipt_mc) < neventsmin || counts["data"][s][ieta][ialpha]->GetBinContent(ipt_dt) < neventsmin){
                    //if (rp_debug) cout << ipt << " pt " <<pt << " g->GetX()[i] " << g->GetX()[i] << " nentries MC "<< counts["mc"][s][ieta][ialpha]->GetBinContent(ipt) << " nentries data " <<  counts["data"][s][ieta][ialpha]->GetBinContent(ipt)<<  " y: " << g->GetY()[i] << endl;
		    if (rp_debug) cout << "ratio " << s << " " << t << " ieta " << ieta << " ialpha " << ialpha << " ipt " << ipt << " pt " << pt << " REMOVE POINT... nentries MC "<< counts["mc"][s][ieta][ialpha]->GetBinContent(ipt_mc) << " data " <<  counts["data"][s][ieta][ialpha]->GetBinContent(ipt_dt) << " y: " << g->GetY()[i] << "+/-" << g->GetEY()[i] << endl;
                    g->RemovePoint(i);
                  }
                }
                else if (nentries < neventsmin){
                  //if (rp_debug) cout << ipt << " pt " <<pt << " g->GetX()[i] " << g->GetX()[i] << " nentries "<< nentries <<  " y: " << g->GetY()[i] << endl;
		  if (rp_debug) cout << d << " " << s << " " << t << " ieta " << ieta << " ialpha " << ialpha << " ipt " << ipt << " pt " <<pt << " REMOVE POINT... nentries "<< nentries <<  " y: " << g->GetY()[i] << "+/-" << g->GetEY()[i] << endl;
                  g->RemovePoint(i);
                }
              } // zmm/zee

              //else if (i==0)continue;
              //remove points where difference of central value between points is larger than 5*sigma(higher pt point) to remove spurious high pt gamma+jet points with presumably low stats
              //else if (i>0 && i==g->GetN()-1 && fabs(g->GetY()[i]-g->GetY()[i-1])>g->GetEY()[i]*5.) g->RemovePoint(i);
	    } // for i

	    // patch MC/data to data/MC for dijet samples
	    if (s=="dijet" && d=="ratio") {
	      for (int i = 0; i != g->GetN(); ++i) {
	    	double x = g->GetX()[i];
	    	double ex = g->GetEX()[i];
	    	double y = g->GetY()[i];
	    	double ey = g->GetEY()[i];
	    	assert(y!=0);
	    	g->SetPoint(i, x, y!=0 ? 1./y : 0);
	    	g->SetPointError(i, ex, y!=0 ? ey/(y*y) : 0);
	      }
	    }

	    // patch L1L2L3Res to L1L2Res for "zjet" data
	    // (ratio calculated on the fly from data so patched also)
	    if (false && s=="zjet" && d=="data" && (t=="mpfchs1" || t=="ptchs")) {
	      //if (s=="zjet" && d=="data" && (t=="ptchs")) {
	      assert(hzjes);
	      for (int i = 0; i != g->GetN(); ++i) {
		double x = g->GetX()[i];
		double jes = hzjes->GetBinContent(hzjes->FindBin(x));
		g->SetPoint(i, x, g->GetY()[i]*jes);
	      } // for i
	    } // patch zjet
	    
	    // PATCH hadw ratio uncertainties to 2 or sqrt(2) larger
	    // Not sure if justified, but scatter looks too large atm
	    if (s=="hadw" && d=="ratio" && (t=="mpfchs1" || t=="ptchs")
		&& !isUL16) {
	      double k(1);
	      if (t=="mpfchs1") k = 2.;
	      if (t=="ptchs")  k = sqrt(2.);
	      for (int i = 0; i != g->GetN(); ++i) {
		g->SetPointError(i, g->GetEX()[i], g->GetEY()[i]*k);
	      }
	    }

	    // patch Z+jet pT ratio uncertainty (80X-590/pb)
	    if ((s=="zeejet" || s=="zmmjet") && d=="ratio") {

	      TGraphErrors *gfixd = grs["data"][t][s][ieta][ialpha];
	      TGraphErrors *gfixm = grs["mc"][t][s][ieta][ialpha];
              assert(gfixd);
              assert(gfixm);

	      for (int i = 0; i != g->GetN(); ++i) {

		double x = g->GetX()[i];
		double ex = g->GetEX()[i];

		double exd(0), eyd(0), eym(0);
		for (int jd = 0; jd != gfixd->GetN(); ++jd) {

		  if (fabs(x-gfixd->GetX()[jd])<ex) {
		    //assert(eyd==0);
		    exd = gfixd->GetEX()[jd];
		    eyd = gfixd->GetEY()[jd];
		  }
		} // for jd
		//assert(eyd!=0);
		for (int jm = 0; jm != gfixm->GetN(); ++jm) {

		  if (fabs(x-gfixm->GetX()[jm])<ex) {
		    //assert(eym==0);
		    eym = gfixm->GetEY()[jm];
		  }
		} // for jm
		//assert(eym!=0);

		g->SetPoint(i, x, g->GetY()[i]);
		g->SetPointError(i, exd, sqrt(eyd*eyd+eym*eym));

	      } // for i
	    } // patch ratio uncertainty

	    // Patch missing multijet ratio
	    if (s=="multijet" && d=="ratio") {

	      TGraphErrors *gfixd = grs["data"][t][s][ieta][ialpha];
	      TGraphErrors *gfixm = grs["mc"][t][s][ieta][ialpha];
              assert(gfixd);
              assert(gfixm);

	      for (int i = 0; i != g->GetN(); ++i) {

		double x = g->GetX()[i];
		double ex = 0.05*x;//g->GetEX()[i];

		double xd(0), yd(0), exd(0), eyd(0);
		double xm(0), ym(0), exm(0), eym(0);
		for (int jd = 0; jd != gfixd->GetN(); ++jd) {

		  if (fabs(x-gfixd->GetX()[jd])<ex) {
		    //assert(eyd==0);
		    xd = gfixd->GetX()[jd];
		    yd = gfixd->GetY()[jd];
		    exd = gfixd->GetEX()[jd];
		    eyd = gfixd->GetEY()[jd];
		  }
		} // for jd
		//assert(eyd!=0);
		for (int jm = 0; jm != gfixm->GetN(); ++jm) {

		  if (fabs(x-gfixm->GetX()[jm])<ex) {
		    //assert(eym==0);
		    xm = gfixm->GetX()[jm];
		    ym = gfixm->GetY()[jm];
		    exm = gfixm->GetEX()[jm];
		    eym = gfixm->GetEY()[jm];
		  }
		} // for jm
		//assert(eym!=0);

		if (t=="crecoil") {
		  g->SetPoint(i, 0.5*(xd+xm), yd!=0 && ym!=0 ? 0.5*(yd+ym) : 0);
		  g->SetPointError(i, 0.5*fabs(xd-xm), yd!=0 && ym!=0 ?
				   sqrt(0.5*eyd*eyd + 0.5*eym*eym
					+ pow(yd-ym,2)) : 0);
		}
		else { // mpfchs1, ptchs
		  g->SetPoint(i, 0.5*(xd+xm), ym ? yd / ym : 0.);
		  g->SetPointError(i, 0.5*fabs(xd-xm),
				   yd!=0 && ym!=0 ?
				   yd / ym * sqrt(pow(eyd/yd,2) + pow(eym/ym,2))
				   : 0.);
		}
	      }

	      for (int i = g->GetN()-1; i != -1; --i) {
		// Clean out spurious empty points
		if (g->GetY()[i]==0 && g->GetEY()[i]==0) g->RemovePoint(i);
	      }
	    } // patch missing multijet ratio

	    // Mass corrections for Zmm
	    // Limit uncertainty at high pT to 0.5% (about twice correction)
 	    if (correctZmmMass && s=="zmmjet" && (d=="data" || d=="ratio") &&
		(t=="mpfchs1" || t=="mpf1" || t=="ptchs")) {
 	      for (int i = 0; i != g->GetN(); ++i) {
		double pt = g->GetX()[i];
 		double ipt = hmzmm->FindBin(pt);
 		double ek = min(0.005,hmzmm->GetBinError(ipt));
		if (useFixedFit) ek = max(fitUncMin,f1ezmm->Eval(pt));
		double k = f1mzmm->Eval(pt);
		//double ek = femzmm->Eval(pt);
 		g->SetPoint(i, g->GetX()[i], g->GetY()[i]*k);
		if (correctUncert)
		  g->SetPointError(i, g->GetEX()[i], 
				   sqrt(pow(g->GetEY()[i]*k,2) + ek*ek));
 	      }
 	    }
	    // Mass corrections for Z+jet
 	    if (correctZMass && s=="zjet" && (d=="data" || d=="ratio") &&
		(t=="mpfchs1" || t=="mpf1" || t=="ptchs")) {
 	      for (int i = 0; i != g->GetN(); ++i) {
		double pt = g->GetX()[i];
		double ek = max(fitUncMin,f1ez->Eval(pt));
		double k = f1mz->Eval(pt);
 		g->SetPoint(i, g->GetX()[i], g->GetY()[i]*k);
		if (correctUncert)
		  g->SetPointError(i, g->GetEX()[i], 
				   sqrt(pow(g->GetEY()[i]*k,2) + ek*ek));
 	      }
 	    }
	    // Mass corrections for Zee
	    // Limit uncertainty at high pT to 0.5% (about half correction)
	    if (correctZeeMass && s=="zeejet" && (d=="data" || d=="ratio") &&
		(t=="mpfchs1" || t=="mpf1" || t=="ptchs")) {
	      for (int i = 0; i != g->GetN(); ++i) {
		double pt = g->GetX()[i];
		double ipt = hmzee->FindBin(pt);
		double ek = min(0.005,hmzee->GetBinError(ipt));
		if (useFixedFit) ek = max(fitUncMin,f1ezee->Eval(pt));
		double k = f1mzee->Eval(pt);
		//double ek = femzee->Eval(pt);
		g->SetPoint(i, g->GetX()[i], g->GetY()[i]*k);
		if (correctUncert)
		  g->SetPointError(i, g->GetEX()[i], 
				   sqrt(pow(g->GetEY()[i]*k,2) + ek*ek));
	      }
	    }
	    // Zee mass corrections for photon+jet
	    // Note fix (post Legacy2016): pT,gamma=pT,lepton=pT,Z/2
	    // Limit uncertainty at high pT to 0.5% (about half correction)
	      /*
	    if (correctGamMass && s=="gamjet" && (d=="data" || d=="ratio") &&
		(t=="mpfchs1" || t=="ptchs" || t=="mpf1" ||
		 t=="mpfn" || t=="mpfu") // latest
	      //(t=="mpfchs1" || t=="ptchs" || t=="mpf1") // DP_2021 (bug!)
		&& !(d=="ratio" && t=="mpf1")) {
	      // NB: mpf1 ratio calculated on the fly => "data" propagates to it
	      */
	    if (correctGamMass && s=="gamjet" && d=="data" &&
		(t=="mpfchs1"||t=="ptchs"||t=="mpf1"||t=="mpfn"||t=="mpfu")) {
	      // NB: ratios calculated on the fly => "data" propagates to it

	      for (int i = 0; i != g->GetN(); ++i) {
		//assert(!correctGamScale); // UL17
		double ptgam = g->GetX()[i];
		double pt = 2*ptgam;
		//double ipt = hmzee->FindBin(pt);
		double ipt = min(hmzee->FindBin(fzmmptmax), hmzee->FindBin(pt));
		double ek = min(0.005,hmzee->GetBinError(ipt));
		if (useFixedFit) ek = max(fitUncMin,f1egam->Eval(pt));
		double k = f1mgam->Eval(pt); // vs pTZ
		g->SetPoint(i, g->GetX()[i], g->GetY()[i]*k);
		if (correctUncert)
		  g->SetPointError(i, g->GetEX()[i], 
				   sqrt(pow(g->GetEY()[i]*k,2) + ek*ek));
	      }
	    } // correctGamMass

	    // Photon scale correction (from drawGamVsZmm)
	    // No separate unceratinty added, is already in global fit
	    //if (correctGamScale && s=="gamjet" && (d=="data" || d=="ratio") &&
	    if (correctGamScale && s=="gamjet" && d=="data" &&
		(t=="mpfchs1"||t=="ptchs"||t=="mpf1"||t=="mpfn"||t=="mpfu")) {
	      //&& !(d=="ratio" && t=="mpf1")) {
	      // NB: mpf1 ratio calculated on the fly => "data" propagates to it

	      //assert(!correctGamMass); // UL17
	      for (int i = 0; i != g->GetN(); ++i) {
		g->SetPoint(i, g->GetX()[i], g->GetY()[i]*valueGamScale);
	      }
	    } // correctGamScale


	    // Multijet JER-hybrid corrections
	    if (correctMultijetLeading && multijetMode=="leading" &&
		s=="multijet" && (d=="data" || d=="ratio")) {

	      for (int i = 0; i != g->GetN(); ++i) {
		double pt = g->GetX()[i];
		double k0 = f2->Eval(pt); // MJBlead / JER-hybrid - 1 (%)
		double kHybridOverLead = 1./(1+0.01*k0);
		g->SetPoint(i, pt, g->GetY()[i] * kHybridOverLead);
	      }
	    }
	    // Multijet recoil scale correction
	    if (correctMultijetRecoilScale && s=="multijet" &&
		d=="data" && // ratio calculated on the fly, so propagates there
		(t=="mpfchs1" || t=="ptchs" || t=="mpf1")) {
	      for (int i = 0; i != g->GetN(); ++i) {
		double pt = g->GetX()[i];
		if (useMultijetRecoilScaleFunction) {
		  double kScale = f1mjb->Eval(pt);
		  if (halveMultijetRecoilScaleFunction) 
		    kScale = 0.5*(kScale+1);
		  g->SetPoint(i, pt, g->GetY()[i]/kScale * multijetRecoilScale);
		}
		else 
		  g->SetPoint(i, pt, g->GetY()[i] * multijetRecoilScale);
	      }
	    }
	    // Remove targeted single outliers from multijet results
	    if (snipeMultijetBins && s=="multijet") {
	      for (int i = g->GetN()-1; i != -1; --i) {
		double pt = g->GetX()[i];
		//if (epoch=="2016BCDEF" && fabs(pt-1400)<100) {
		if (isUL16 && isAPV && fabs(pt-1400)<100) {
		  g->RemovePoint(i);
		}
	      }
	    }
	    // Remove targeted single outliers from inclusive jet results
	    if (snipeIncjetBins && s=="incjet") {
	      for (int i = g->GetN()-1; i != -1; --i) {
		double pt = g->GetX()[i];
		if (isUL16 && isAPV && fabs(pt-69)<5) {
		  g->RemovePoint(i);
		}
	      }
	    }
	    // Remove targeted single outliers from Zll+jet results
	    if (snipeZlljetBins && s=="zlljet") {
	      for (int i = g->GetN()-1; i != -1; --i) {
		double pt = g->GetX()[i];
		//if (isUL16 && isAPV && epoch=="2016BCD" && fabs(pt-60)<10) {
		if ((isUL16 && isAPV && fabs(pt-65)<5) ||
		    (isUL18 && fabs(pt-37.5)<2.5)) {
		  g->RemovePoint(i);
		}
	      }
	    }
	    

	    dout->cd("orig");
	    g_orig->Write();

	    dout->cd();

	    // Set uniform naming scheme and graphical style
	    g->SetName(Form("%s_%s_a%1.0f",tt,ss,100.*alphas[ialpha]));
	    g->UseCurrentStyle(); // Basic TDR style
	    g->SetMarkerStyle(style[s][t]);
	    g->SetMarkerColor(color[s]);
	    g->SetLineColor(color[s]);
	    if (s=="pfjet" || (s=="zjet" && isfrac) ||
		(s=="zlljet" && isfrac) || (s=="gamjet" && isfrac) ||
		ismpfc) {
	      g->SetMarkerStyle((d=="mc" && style[s+"_"+d][t]!=0) ?
				style[s+"_"+d][t] : style[s][t]);
	      g->SetMarkerColor(color[s+"_"+t]!=0 ? color[s+"_"+t] : color[s]);
	      g->SetLineColor(color[s+"_"+t]!=0 ? color[s+"_"+t] : color[s]);
	    }
	    if (sp=="zjet" || sp=="gamjet") {
	      g->SetMarkerStyle(style[sp][t]);
	      g->SetMarkerColor(color[s]!=0 ? color[s] : color[sp]);
	      g->SetLineColor(color[s]!=0 ? color[s] : color[sp]);
	    }
	    g->SetMarkerSize(size[alpha]);
	    g->SetDrawOption("SAMEP");

            //if (  s=="zmmjet" ||  s=="zeejet")cout << Form("Writing out %s_%s_a%1.0f_eta%02.0f-%02.0f",tt,ss,100.*alphas[ialpha] ,eta1*10.,eta2*10.) << endl;
	    g->Write();

	    grs[d][t][s][ieta][ialpha] = g;

	  } // for ialpha
	} // for iset
      } // for itype
    } // for ieta
  } // for itier

  // Save mass histos as well
  cout << "Save mass histos" << endl << flush;
  fout->cd("ratio/eta00-13");
  hmzee->Write("mass_zeejet_a30");
  hmzmm->Write("mass_zmmjet_a30");
  hmz->Write("mass_zjet_a30"); // v28
  hmzee1->Write("mass_zeejet_a100");
  hmzmm1->Write("mass_zmmjet_a100");
  hmz1->Write("mass_zjet_a100"); // v28

  fout->cd("data/eta00-13");
  hmzee_dt->Write("mass_zeejet_a30");
  hmzmm_dt->Write("mass_zmmjet_a30");
  hmz_dt->Write("mass_zjet_a30");
  hmzee1_dt->Write("mass_zeejet_a100");
  hmzmm1_dt->Write("mass_zmmjet_a100");
  hmz1_dt->Write("mass_zjet_a100");

  fout->cd("mc/eta00-13");
  hmzee_mc->Write("mass_zeejet_a30");
  hmzmm_mc->Write("mass_zmmjet_a30");
  hmz_mc->Write("mass_zjet_a30");
  hmzee1_mc->Write("mass_zeejet_a100");
  hmzmm1_mc->Write("mass_zmmjet_a100");
  hmz1_mc->Write("mass_zjet_a100");

  if (fdj) fdj->Close();
  if (fp) fp->Close();
  if (fzee) fzee->Close();
  if (fzmm) fzmm->Close();
  if (fz) fz->Close();
  curdir->cd();


  /////////////////////////////////////////////////////
  // Calculate JEC central values and uncertainties, //
  // and store them for easy visualialization later  //
  /////////////////////////////////////////////////////

  if (rp_debug) cout << "Instantiating JECs..." << endl << flush;

  map<string,const char*> mera;
  mera["2018ABCD"] = "ABCD";
  mera["2018A"] = "A";
  mera["2018B"] = "B";
  mera["2018C"] = "C";
  mera["2018D"] = "D";
  mera["2017BCDEF"] = "BCDEF";
  mera["2016BCDEF"] = "BCDEF";
  mera["2016GH"] = "FGH";
  mera["2016BCD"] = "BCD";
  mera["2016EF"] = "EF";
  mera["2016BCDEF"] = "BCDEF";
  mera["2016BCDEFGH"] = "BCDEFGH";

  // Calculate L2L3Res with JEC uncertainty
  {

    const char *s, *s2;
    // Usual directory for text files
    const char *cd = "CondFormats/JetMETObjects/data";
    // Directory for development version
    const char *cdx = "../JERCProtoLab/Summer19UL16/global_fit/V5M1";
    const char *ce = epoch.c_str();
    
    // New JEC for plotting on the back and mcjec for softrad3.C
    // ** Also used as reference for CorLevel=="L1L2L3Res" **
    // So be careful when running minitools/createL2L3Res.C
    FactorizedJetCorrector *mcjec(0);
    FactorizedJetCorrector *jec(0), *jeca(0);
    FactorizedJetCorrector *jecb(0), *jecc(0), *jecd(0), *jece(0), *jecf(0);
    FactorizedJetCorrector *jecg(0), *jech(0);
    double jecwa(0), jecwb(0), jecwc(0), jecwd(0), jecwe(0), jecwf(0);
    double jecwg(0), jecwh(0);
    // all of Run 2 in a more systematic fashion
    vector< pair<double,FactorizedJetCorrector*> > vjec;
    {
      if (rp_debug) cout << "Individual JECs..." << endl << flush;

      if (isUL18) {
	jec = getFJC("","",Form("Summer19UL18_Run%s_V5_DATA_L2L3Residual",
				epoch=="2018ABCD" ? "A" : mera[epoch]));
	mcjec = getFJC("",Form("Summer19UL18_Run%s_V5_DATA_L2Relative",
			       epoch=="2018ABCD" ? "A" : mera[epoch]));
      }
      else if (isUL17) {
	jec =getFJC("","",Form("Summer19UL17_Run%s_V5_DATA_L2L3Residual",
			       epoch=="2017BCDEF" ? "C" : mera[epoch]));
	mcjec = getFJC("",Form("Summer19UL17_Run%s_V5_DATA_L2Relative",
			       epoch=="2017BCDEF" ? "C" : mera[epoch]));
      }
      // Add MC JEC for Runcl for UL16
      else if (isUL16 && !isAPV) {
	jec = getFJC("",Form("Summer19UL16_Run%s_V7_DATA_L2L3Residual",mera[epoch]),"",cd); // was V5
	//jec = getFJC("",Form("Summer19UL16_Run%s_V5M1_DATA_L2L3Residual",mera[epoch]),"",cdx);
	mcjec = getFJC("",Form("Summer19UL16_Run%s_V7_DATA_L2Relative",mera[epoch]),"",cd); // was V5
      }
      else if (isUL16 && isAPV) {
	jec = getFJC("",Form("Summer19UL16APV_Run%s_V7_DATA_L2L3Residual",mera[epoch]),"",cd); // was V5
	//jec = getFJC("",Form("Summer19UL16_Run%s_V5M1_DATA_L2L3Residual",mera[epoch]),"",cdx);
	mcjec = getFJC("",Form("Summer19UL16APV_Run%s_V7_DATA_L2Relative",mera[epoch]),"",cd); // was V5
      }
      else if (isRun2) { // clone 2018A
	//jec = getFJC("","","Summer19UL18_RunA_V5_DATA_L2L3Residual");
	jec = getFJC("","","mergeL2L3ResTextFiles_Run2Test","textFiles");
	mcjec = getFJC("","Summer19UL18_RunA_V5_DATA_L2Relative");
      }
      else if (isLowPU) {
	jec =getFJC("","",Form("Summer19UL17_Run%s_V5_DATA_L2L3Residual","F"));
	mcjec = getFJC("",Form("Summer19UL17_Run%s_V5_DATA_L2Relative","F"));
      }
      else
	assert(false);


      if (rp_debug) cout << "Combined JECs..." << endl << flush;
      
      if (epoch=="2018ABCD") {
	
	jeca =getFJC("","","Summer19UL18_RunA_V5_DATA_L2L3Residual");
	jecb =getFJC("","","Summer19UL18_RunB_V5_DATA_L2L3Residual");
	jecc =getFJC("","","Summer19UL18_RunC_V5_DATA_L2L3Residual");
	jecd =getFJC("","","Summer19UL18_RunD_V5_DATA_L2L3Residual");

	double lumtot = 14.0+7.1+6.9+31.9; //59.9/fb
	jecwa = 14.0/lumtot;
	jecwb = 7.1/lumtot;
	jecwc = 6.9/lumtot;
	jecwd = 31.9/lumtot;
      }
      else if (epoch=="2017BCDEF") {

	jecb =getFJC("","","Summer19UL17_RunB_V5_DATA_L2L3Residual");
	jecc =getFJC("","","Summer19UL17_RunC_V5_DATA_L2L3Residual");
	jecd =getFJC("","","Summer19UL17_RunD_V5_DATA_L2L3Residual");
	jece =getFJC("","","Summer19UL17_RunE_V5_DATA_L2L3Residual");
	jecf =getFJC("","","Summer19UL17_RunF_V5_DATA_L2L3Residual");
	
	// luminosity from Hugues Lattaud, photon+jet 11 June 2018
	double lumtot = 4.8+9.6+4.2+9.3+13.4; // 41.3/fb
	jecwb = 4.8/lumtot;
	jecwc = 9.6/lumtot;
	jecwd = 4.2/lumtot;
	jecwe = 9.3/lumtot;
	jecwf = 13.4/lumtot;
      }
      else if (epoch=="2016BCDEF") {

	// were V5
	jecb = getFJC("","","Summer19UL16APV_RunBCD_V7_DATA_L2L3Residual");
	jece = getFJC("","","Summer19UL16APV_RunEF_V7_DATA_L2L3Residual");
	//jecb = getFJC("","","Summer19UL16_RunBCD_V5M1_DATA_L2L3Residual",cdx);
	//jece = getFJC("","","Summer19UL16_RunEF_V5M1_DATA_L2L3Residual",cdx);
	
	// luminosity from Hugues Lattaud, photon+jet 11 June 2018
	double lumtot = 5.9+2.6+4.4 +4.0+3.2; // 20.1 (vs 19.8?)
	jecwb = (5.9+2.6+4.4)/lumtot; // 12.9
	jecwe = (4.0+3.2)/lumtot; // 7.2 (vs 6.8?)
      }
      else if (epoch=="2016GH") {

	jecg = getFJC("","","Summer19UL16_RunFGH_V7_DATA_L2L3Residual");
	jecwg = 1;
      }
      else if (isRun2) {

#define PAIR(a,b) (make_pair<double,FactorizedJetCorrector*>((a),getFJC("","",(b))))

	vjec.push_back(PAIR(12.9,"Summer19UL16APV_RunBCD_V7_DATA_L2L3Residual"));
	vjec.push_back(PAIR( 6.8,"Summer19UL16APV_RunEF_V7_DATA_L2L3Residual"));
	vjec.push_back(PAIR(16.8,"Summer19UL16_RunFGH_V7_DATA_L2L3Residual"));

	vjec.push_back(PAIR( 4.8,"Summer19UL17_RunB_V5_DATA_L2L3Residual"));
	vjec.push_back(PAIR( 9.6,"Summer19UL17_RunC_V5_DATA_L2L3Residual"));
	vjec.push_back(PAIR( 4.2,"Summer19UL17_RunD_V5_DATA_L2L3Residual"));
	vjec.push_back(PAIR( 9.3,"Summer19UL17_RunE_V5_DATA_L2L3Residual"));
	vjec.push_back(PAIR(13.4,"Summer19UL17_RunF_V5_DATA_L2L3Residual"));

	vjec.push_back(PAIR(14.0,"Summer19UL18_RunA_V5_DATA_L2L3Residual"));
	vjec.push_back(PAIR( 7.1,"Summer19UL18_RunB_V5_DATA_L2L3Residual"));
	vjec.push_back(PAIR( 6.9,"Summer19UL18_RunC_V5_DATA_L2L3Residual"));
	vjec.push_back(PAIR(31.9,"Summer19UL18_RunD_V5_DATA_L2L3Residual"));
      }
      // Alternative combined JEC
      //else if (epoch=="2016BCDEF") {
      //
      //jecb =getFJC("","","Summer19UL16_RunBCDEF_V5_DATA_L2L3Residual",
      //	      cd);
      //jecwb = 1;
      //}
    }


    if (rp_debug) cout << "Loading reference JECs..." << endl << flush;
    
    FactorizedJetCorrector *jecrun1, *jecold, *jecl1flat,*jecl1rcdt, *jecl1pt,
      *jecl1mc, *jecl1rcmc, *jecl1dt, *jecl1c, *jecl1s, *jecl1rc,
      *jecl2c, *jecl2s;

    // Reference Run I JEC for plotting on the back
    jecrun1 =   getFJC("","","Winter14_V8_DATA_L2L3Residual_AK5PFchs"); 
    
    // Store old JEC for undoing it in global fit
    // NB: global fit ideally needs a temporary pT-flat L3Res as input
    // But even with this pT-dependent L2Res can cause problems
    if (isUL18) {
      //jecold = getFJC("","",Form("Summer19UL18_Run%s_V5_DATA_L2L3Residual",
      //			 epoch=="2018ABCD" ? "C" : mera[epoch]));
      jecold = getFJC("","",Form("Summer19UL18_Run%s_V3M1_DATA_L2L3Residual",
				 epoch=="2018ABCD" ? "C" : mera[epoch]));
    }
    else if (isUL17) {
      //jecold = getFJC("","",Form("Summer19UL17_Run%s_V5_DATA_L2L3Residual",
      //epoch=="2017BCDEF" ? "C" : mera[epoch]));
      jecold = getFJC("","",Form("Summer19UL17_Run%s_V2M5_SimpleL1_DATA_L2L3Residual",mera[epoch]));
    }
    else if (isUL16 && !isAPV) {
      if (CorLevel=="L1L2L3Res")
	jecold = getFJC("","",Form("Summer19UL16_Run%s_V7_DATA_L2L3Residual",
				   mera[epoch])); // was V5
      else
	assert(false);
    }
    else if (isUL16 && isAPV) {
      if (CorLevel=="L1L2L3Res")
	jecold = getFJC("","",Form("Summer19UL16APV_Run%s_V7_DATA_L2L3Residual",
				   mera[epoch])); // was V5
      else
	assert(false);
    }
    else if (isRun2) { // use UL18D
      //jecold = getFJC("","","Summer19UL18_RunD_V5_DATA_L2L3Residual");
      jecold = getFJC("","","mergeL2L3ResTextFiles_Run2Test","textFiles");
    }
    else if (isLowPU) {
      jecold = getFJC("","",Form("Summer19UL17_Run%s_V2M5_SimpleL1_DATA_L2L3Residual","F"));
    }
    else
      assert(false);


    if (rp_debug) cout << "Loading L1 JECs..." << endl << flush;
 
    // Difference between pT-dependent and flat L1
    if (isUL17 || isUL18 || isRun2 || isLowPU) {
      jecl1flat = getFJC("Summer19UL17_RunC_V1_SimpleL1_DATA_L1RC");
      jecl1rcdt = getFJC("Summer19UL17_RunC_V1_SimpleL1_DATA_L1RC");
      jecl1pt =   getFJC("Summer19UL17_RunC_V1_SimpleL1_DATA_L1FastJet");
      jecl1mc =   getFJC("Summer19UL17_V1_SimpleL1_MC_L1FastJet");
      jecl1rcmc = getFJC("Summer19UL17_V1_SimpleL1_MC_L1RC");
      jecl1dt =   getFJC("Summer19UL17_RunC_V1_SimpleL1_DATA_L1FastJet");
    }
    else if (isUL16 && !isAPV) {
      // were V5 
      jecl1flat = getFJC("Summer19UL16_RunFGH_V7_DATA_L1RC");
      jecl1rcdt = getFJC("Summer19UL16_RunFGH_V7_DATA_L1RC");
      jecl1pt =   getFJC("Summer19UL16_RunFGH_V7_DATA_L1FastJet");
      jecl1mc =   getFJC("Summer19UL16_RunFGH_V7_DATA_L1RC");
      jecl1rcmc = getFJC("Summer19UL16_RunFGH_V7_DATA_L1RC");
      jecl1dt =   getFJC("Summer19UL16_RunFGH_V7_DATA_L1FastJet");
    }
    else if (isUL16 && isAPV) {
      // were V5
      jecl1flat = getFJC("Summer19UL16APV_RunBCDEF_V7_DATA_L1RC");
      jecl1rcdt = getFJC("Summer19UL16APV_RunBCDEF_V7_DATA_L1RC");
      jecl1pt =   getFJC("Summer19UL16APV_RunBCDEF_V7_DATA_L1FastJet");
      jecl1mc =   getFJC("Summer19UL16APV_RunBCDEF_V7_DATA_L1RC");
      jecl1rcmc = getFJC("Summer19UL16APV_RunBCDEF_V7_DATA_L1RC");
      jecl1dt =   getFJC("Summer19UL16APV_RunBCDEF_V7_DATA_L1FastJet");
    }
    else 
      assert(false);

    // With L2Relative included
    if (isUL17 || isUL18 || isRun2 || isLowPU) {
      jecl1c =    getFJC("Summer19UL17_V1_ComplexL1_MC_L1FastJet",
			 // This is on purpose SimpleL1_L2 for ComplexL1_L1
			 // Need L1C L2 on same footing with L1S and L1RC
			 "Summer19UL17_V1_SimpleL1_MC_L2Relative");
      jecl1s =    getFJC("Summer19UL17_V1_SimpleL1_MC_L1FastJet",
			 "Summer19UL17_V1_SimpleL1_MC_L2Relative");
      jecl1rc =   getFJC("Summer19UL17_V1_SimpleL1_MC_L1RC",
			 "Summer19UL17_V1_SimpleL1_MC_L2Relative");
      jecl2c =    getFJC("","Summer19UL17_V1_ComplexL1_MC_L2Relative");
      jecl2s =    getFJC("","Summer19UL17_V1_SimpleL1_MC_L2Relative");
    }
    else if (isUL16 && !isAPV) {
      // were V5
      jecl1c =    getFJC("Summer19UL16_V7_MC_L1FastJet",
			 "Summer19UL16_V7_MC_L2Relative");
      jecl1s =    getFJC("Summer19UL16_V7_MC_L1FastJet",
			 "Summer19UL16_V7_MC_L2Relative");
      jecl1rc =   getFJC("Summer19UL16_V7_MC_L1RC",
			 "Summer19UL16_V7_MC_L2Relative");
      jecl2c =    getFJC("","Summer19UL16_V7_MC_L2Relative");
      jecl2s =    getFJC("","Summer19UL16_V7_MC_L2Relative");
    }
    else if (isUL16 && isAPV) {
      // were V5
      jecl1c =    getFJC("Summer19UL16APV_V7_MC_L1FastJet",
			 "Summer19UL16APV_V7_MC_L2Relative");
      jecl1s =    getFJC("Summer19UL16APV_V7_MC_L1FastJet",
			 "Summer19UL16APV_V7_MC_L2Relative");
      jecl1rc =   getFJC("Summer19UL16APV_V7_MC_L1RC",
			 "Summer19UL16APV_V7_MC_L2Relative");
      jecl2c =    getFJC("","Summer19UL16APV_V7_MC_L2Relative");
      jecl2s =    getFJC("","Summer19UL16APV_V7_MC_L2Relative");
    }
    else
      assert(false);


    if (rp_debug) cout << "Loading L1 uncertainty JECs..." << endl << flush;

    // New additions (RC: mean vs median; L1: Simple vs SemiSimple)
    FactorizedJetCorrector *jrcdnom(0),*jrcmnom(0),*jrcdalt(0),*jrcmalt(0); // 1
    FactorizedJetCorrector *jl1dnomsf(0), *jl1daltsf(0), *jl1mnomsf(0); // 2a
    FactorizedJetCorrector *jl1dnom(0),*jl1mnom(0),*jl1dalt(0),*jl1malt(0);// 2b

    // Use consistent BCDEF files instead of the official files
    if (isUL17 || isUL18 || isUL16 || isRun2 || isLowPU) { // Add UL16, although uncert. may be old
      jrcdnom = getFJC("UL17_RunBCDEF_V1_DATA_L1RC"); // (1)
      jrcmnom = getFJC("UL17_RunBCDEF_V1_MC_L1RC"); // (1)
      jrcdalt = getFJC("UL17_MED_RunBCDEF_V1_DATA_L1RC"); // (1)
      jrcmalt = getFJC("UL17_MED_RunBCDEF_V1_MC_L1RC"); // (1)

      // Use place-holder for the missing custom data files
      jl1dnomsf = getFJC("UL17_RunBCDEF_L1Simple_L1FastJet_AK4PFchs"); // (2a)
      jl1daltsf = getFJC("UL17_RunBCDEF_L1Simple_L1FastJet_AK4PFMEDchs"); // (2a)
      jl1mnomsf = getFJC("Summer19UL17_V1_SimpleL1_MC_L1FastJet"); //xtra

      // Shell magic to create new new data file from old data and MC files:
      //   cp -pi [Summer19UL17...MC...] mc.txt
      //   cp -pi [Summer19UL17...DATA...] dt.txt
      //   cat dt.txt | awk '{print "     "$12"    "$13}' > sf.txt
      //   paste mc.txt sf.txt > dt1.txt
      //   sed -e "s/^M//" dt1.txt > dt2.txt
      // (above, use Ctrl+V Ctrl+M to produce the carriage return ^M)
      //   sed -e "s/8       0/10       0/" dt2.txt > dt3.txt
      // ( sed -e "s/9       0/11       0/" dt2.txt > dt3.txt )
      //   head -n1 dt.txt > dt4.txt
      // ( head -n1 mc.txt > dt4.txt )
      // ( and by hand modificatinos for L1SemiSimple )
      //   tail -n82 dt3.txt >> dt4.txt

      // List of data files still needed (RunE->RunBCDEF, SimpleL1->SemiSimpleL1)
      jl1dnom = getFJC("UL17_RunBCDEF_L1Simple_L1FastJet", //(2b)
		       "Summer19UL17_V1_SimpleL1_MC_L2Relative"); // (2b)
      jl1mnom = getFJC("Summer19UL17_V1_SimpleL1_MC_L1FastJet", // (2b)
		       "Summer19UL17_V1_SimpleL1_MC_L2Relative"); // (2b)
      jl1dalt = getFJC("UL17_RunBCDEF_L1SemiSimple_L1FastJet", //(2b)
		       "ParallelMCL1_L2Relative_AK4PFchs_L1SemiSimple"
		       "_L2L3Splines_ptclip8"); // (2b)
      jl1malt = getFJC("ParallelMCL1_L1FastJet_AK4PFchs_L1SemiSimple", // (2b)
		       "ParallelMCL1_L2Relative_AK4PFchs_L1SemiSimple"
		       "_L2L3Splines_ptclip8"); // (2b)
    }
    

    if (rp_debug) cout << "Loading uncertainty sources..." << endl << flush;

    // Run I uncertainty => 80XV6 uncertainty
    s = Form("%s/Winter14_V8_DATA_UncertaintySources_AK5PFchs.txt",cd); // V8 Run I
    s2 = "SubTotalAbsolute";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_ref1 = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc_ref1 = new JetCorrectionUncertainty(*p_ref1);

    // Total uncertainty, excluding Flavor and Time
    if (isUL18) {
      s = Form("%s/Summer19UL18_RunA_V5_DATA_UncertaintySources_AK4PFchs.txt",cd);
    }
    else if (isUL17) {
      s = Form("%s/Summer19UL17_RunC_V5_DATA_UncertaintySources_AK4PFchs.txt",cd);
    }
    else if (isUL16 && !isAPV) {
      s = Form("%s/Summer19UL16_RunFGH_V7_DATA_UncertaintySources_AK4PFchs.txt",cd); // was V5
    }
    else if (isUL16 && isAPV) {
      s = Form("%s/Summer19UL16APV_RunBCDEF_V7_DATA_UncertaintySources_AK4PFchs.txt",cd); // was V5
    }
    else if (isRun2) { // use UL18A
      s = Form("%s/Summer19UL18_RunA_V5_DATA_UncertaintySources_AK4PFchs.txt",cd); 
   }
    else if (isLowPU) {
     s = Form("%s/Summer19UL17_RunF_V5_DATA_UncertaintySources_AK4PFchs.txt",cd);
    }
    else
      assert(false);

    s2 = "TotalNoFlavorNoTime";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_unc = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p_unc);
    //
    s2 = "SubTotalAbsolute";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_ref = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc_ref = new JetCorrectionUncertainty(*p_ref);
    //
    s2 = "SubTotalPt";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_pt = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc_pt = new JetCorrectionUncertainty(*p_pt);
    //
    s2 = "SinglePionHCAL";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_hcal = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc_hcal = new JetCorrectionUncertainty(*p_hcal);
    //
    s2 = "SinglePionECAL";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_ecal = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc_ecal = new JetCorrectionUncertainty(*p_ecal);
    //
    s2 = "FlavorPureGluon";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_glu = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc_glu = new JetCorrectionUncertainty(*p_glu);
    //
    s2 = "FlavorPureQuark";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_uds = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc_uds = new JetCorrectionUncertainty(*p_uds);
    //
    s2 = "FlavorPureCharm";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_cha = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc_cha = new JetCorrectionUncertainty(*p_cha);
    //
    s2 = "FlavorPureBottom";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_bot = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc_bot = new JetCorrectionUncertainty(*p_bot);
    //
    s2 = "FlavorZJet";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_zjt = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc_zjt = new JetCorrectionUncertainty(*p_zjt);
    //
    s2 = "FlavorPhotonJet";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_gjt = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc_gjt = new JetCorrectionUncertainty(*p_gjt);
    //
    s2 = "FlavorQCD";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_qcd = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc_qcd = new JetCorrectionUncertainty(*p_qcd);
    //
    s2 = "SubTotalPileUp";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_pu = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc_pu = new JetCorrectionUncertainty(*p_pu);
    //
    s2 = "TotalNoFlavor";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_noflv = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc_noflv = new JetCorrectionUncertainty(*p_noflv);

    // Loop over eta bins, but do JEC for data/MC ratio only
    for (unsigned int ieta = 0; ieta != etas.size(); ++ieta) {
      
      assert(fout->cd("ratio"));
      TDirectory *dout1 = fout->GetDirectory("ratio"); assert(dout1);
      double eta1 = etas[ieta].first; double eta2 = etas[ieta].second;
      const char *dd1 = Form("eta%02.0f-%02.0f",eta1*10.,eta2*10.);
      assert(dout1->cd(dd1));
      TDirectory *dout2 = dout1->GetDirectory(dd1); assert(dout2);
      dout2->cd();
      cout << "Eta bin:" << eta1 <<"-"<< eta2 << endl;

      const double ptbins[] = {15, 16, 18, 20, 22, 25,
			       30, 35, 40, 50, 60, 70, 85, 100, 125, 155, 180,
			       210, 250, 300, 350, 400, 500, 600, 800, 1000,
			       1200, 1500,
			       1800, 2100, 2400, 2700, 3000, 3300, 3600,
			       3900, 4200, 4500, 4501};
      const int npt = sizeof(ptbins)/sizeof(ptbins[0])-1;

      // Uncertainty bands
      TH1D *herr = new TH1D("herr",";p_{T} (GeV);JEC uncertainty;",
			    npt, &ptbins[0]);
      TH1D *herr_l2l3res = new TH1D("herr_l2l3res",";p_{T} (GeV);L2L3Res;",
				npt, &ptbins[0]);
      TH1D *herr_ref = new TH1D("herr_ref",";p_{T} (GeV);TotalNoFlavorNoTime;",
				npt, &ptbins[0]);
      TH1D *herr_spr = new TH1D("herr_spr",";p_{T} (GeV);SPR uncertainty;",
				npt, &ptbins[0]);
      TH1D *herr_glu = new TH1D("herr_glu",";p_{T} (GeV);Gluon uncertainty;",
				npt, &ptbins[0]);
      TH1D *herr_uds = new TH1D("herr_uds",";p_{T} (GeV);Quark uncertainty;",
				npt, &ptbins[0]);
      TH1D *herr_cha = new TH1D("herr_cha",";p_{T} (GeV);Charm uncertainty;",
				npt, &ptbins[0]);
      TH1D *herr_bot = new TH1D("herr_bot",";p_{T} (GeV);Bottom uncertainty;",
				npt, &ptbins[0]);
      TH1D *herr_zjt = new TH1D("herr_zjt",";p_{T} (GeV);Z+jet uncertainty;",
				npt, &ptbins[0]);
      TH1D *herr_gjt = new TH1D("herr_gjt",";p_{T} (GeV);#gamma+jet uncert.;",
				npt, &ptbins[0]);
      TH1D *herr_qcd = new TH1D("herr_qcd",";p_{T} (GeV);QCD uncertainty;",
				npt, &ptbins[0]);
      TH1D *herr_pu = new TH1D("herr_pu",";p_{T} (GeV);PU uncertainty;",
				npt, &ptbins[0]);
      TH1D *herr_noflv = new TH1D("herr_noflv",";p_{T} (GeV);TotalNoFlavor;",
				npt, &ptbins[0]);
      TH1D *herr_mpf = new TH1D("herr_mpf",";p_{T} (GeV);MPF uncertainty;",
				npt, &ptbins[0]);
      TH1D *herr_pt = new TH1D("herr_pt",";p_{T} (GeV);p_{T} uncertainty;",
			       npt, &ptbins[0]);

      // Run I JEC central value
      TH1D *hrun1 = new TH1D("hrun1",";p_{T} (GeV);Run I JES;",
			     npt, &ptbins[0]);

      // JEC central value
      TH1D *hjes = new TH1D("hjes",";p_{T} (GeV);Old JES;",npt, &ptbins[0]);

      // PileUpPtEnvelope, effectively (flat L1RC vs pT-dependent L1FastJet)
      TH1D *hl1bias = new TH1D("hl1bias",";p_{T} (GeV);"
			       "JEC L1 p_{T} / JEC L1 flat;",
			       npt, &ptbins[0]);
      // Data/MC difference in L1FastJet
      TH1D *hl1diff = new TH1D("hl1diff",";p_{T} (GeV);"
			       "JEC L1 Data / JEC L1 MC;",
			       npt, &ptbins[0]);
      // Impact from 1 GeV difference in input rho (e.g. due to UE)
      TH1D *hl1drho = new TH1D("hl1drho",";p_{T} (GeV);"
			       "JEC L1 Data (#rho) / JEC L1 Data (#rho+1);",
			       npt, &ptbins[0]);
      //
      TH1D *hl1dtomc = new TH1D("hl1dtomc",";p_{T} (GeV);"
				"JEC L1 Simple data / JEC L1Simple MC;",
				npt, &ptbins[0]);
      TH1D *hrcdtomc = new TH1D("hrcdtomc",";p_{T} (GeV);"
				"JEC L1RC data / JEC L1RC MC;",
				npt, &ptbins[0]);
      TH1D *hl1cos = new TH1D("hl1cos",";p_{T} (GeV);"
			      "JEC L1Complex MC / JEC L1Simple MC;",
			      npt, &ptbins[0]);
      TH1D *hl1rcos = new TH1D("hl1rcos",";p_{T} (GeV);"
			       "JEC L1RC MC / JEC L1Simple MC;",
			       npt, &ptbins[0]);
      TH1D *hl2cos = new TH1D("hl2cos",";p_{T} (GeV);"
			      "JEC L2 of L1Complex MC / L1Simple MC;",
			      npt, &ptbins[0]);

  // Subject: Re-designing L1 systematics for UL2017
  // Date: 19 May 2020
  //
  //0) Break L1RC data/MC uncertainty into parts affecting (1) MPF measurement and/or (2b) MPF and pTbal, and (2a,c) into parts absorbed in (L2)L3Res
  //
  //1) MPF: use ratio of L1RC(data/MC) with nominal method (mean) vs alternative method (median) as a nuisance parameter in global fit. The idea here is that we don’t know which estimator is yet optimal for removing bias from MET when filling the hole in PU left by correcting the jet. Doesn’t affect pTbal methods. I’ve in the past included L1 uncertainties as nuisances instead of components in global fit, so this is not entirely new.
  //
  //2) (L2)L3Res: evaluate uncertainty from (2a) uncertainty in data/MC SF applied on nominal method (L1Simple from NoPU), and (2b) uncertainty from alternative method (L1SemiSimple from EpsilonPU)
  //
  //	2a) data/MC SF uncertainty as difference between nominal (mean) SF and SF from alternative (median) method. About same as past PileUpDataMC, but now gets explicitly absorbed in L3Res fit.
  //
  //	2b) L1Method uncertainty as difference between data/MC ratio with nominal (SimpleL1) and alternative (SemiSimpleL1) applied at <rho>. This is about the same as past PileUpEnvelope, which was also absorbed in L3Res fit.
  //
  //	2b) L1Method uncertainty as difference between data/MC ratio with nominal (SimpleL1) and alternative (SemiSimpleL1) applied at <rho>*(1+/-20%) relative to <rho>. The idea here being that data/MC ratio evaluated at <rho> in step (2b) will be absorbed in L3Res, but variations around it will not. The +/-20% covers the range of PU variations in single jet triggers and different eras (2017F is +20%, 2017BCD is -20% relative to 2017BCDEF).
  //
  // Later addition:
  // 2c) => 3a), 3b), for rho+/-20% variations for SF and SemiSimple

      // New additions
      TH1D *hrcmpf = new TH1D("hsys1_rcmpf",";p_{T} (GeV);"
			      "(RC_{dt,alt}/RC_{MC,alt}) / "
			      "(RC_{dt,nom}/RC_{MC,nom})",
			      npt, &ptbins[0]);
      TH1D *hrcmpfdt = new TH1D("hrcmpfdt",";p_{T} (GeV);"
				"RC_{dt,alt} / RC_{dt,nom}",
				npt, &ptbins[0]);
      TH1D *hrcmpfmc = new TH1D("hrcmpfmc",";p_{T} (GeV);"
				"RC_{MC,alt} / RC_{MC,nom}",
				npt, &ptbins[0]);
      TH1D *hrcdt = new TH1D("hrcdt",";p_{T} (GeV);"
			     "RC_{dt,nom}",
			     npt, &ptbins[0]);
      TH1D *hrcmc = new TH1D("hrcmc",";p_{T} (GeV);"
			     "RC_{MC,nom}",
			     npt, &ptbins[0]);
      //
      TH1D *hl1sf = new TH1D("hsys2a_l1sf",";p_{T} (GeV);"
			     "L1_{dt,altsf}/L1_{dt,nomsf}",
			     npt, &ptbins[0]);
      TH1D *hl1sfnom = new TH1D("hl1sfnom",";p_{T} (GeV);"
				"SF_{dt,nom} = L1_{dt,nomsf}/L1_{mc,nom}",
				npt, &ptbins[0]);
      TH1D *hl1sfalt = new TH1D("hl1sfalt",";p_{T} (GeV);"
				"SF_{dt,alt} = L1_{dt,altsf}/L1_{mc,nom}",
				npt, &ptbins[0]);
      TH1D *hl1offdt = new TH1D("hl1offdt",";p_{T} (GeV);"
				"L1offdt = p_{T,corr}*(L1_{dt,nomsf}-1)",
				npt, &ptbins[0]);
      TH1D *hl1dt = new TH1D("hl1dt",";p_{T} (GeV);"
			     //"L1_{dt,nom} = p_{T,corr}/(p_{T,corr}+L1offdt)",
			     "L1_{dt,nomsf}",
			     npt, &ptbins[0]);
      //
      TH1D *hl1 = new TH1D("hsys2b_l1",";p_{T} (GeV);"
			   "(L1_{dt,alt}/L1_{mc,alt}) /"
			   "(L1_{dt,nom}/L1_{mc,nom})",
			   npt, &ptbins[0]);
      TH1D *hl1d = new TH1D("hl1d",";p_{T} (GeV);"
			    "(L1_{dt,alt}/L1_{dt,nom})",
			    npt, &ptbins[0]);
      TH1D *hl1m = new TH1D("hl1m",";p_{T} (GeV);"
			    "(L1_{mc,alt}/L1_{mc,nom})",
			    npt, &ptbins[0]);
      TH1D *hl1dnom = new TH1D("hl1dnom",";p_{T} (GeV);"
			       "L1_{dt,nom})",
			       npt, &ptbins[0]);
      TH1D *hl1mnom = new TH1D("hl1mnom",";p_{T} (GeV);"
			       "L1_{mc,nom})",
			       npt, &ptbins[0]);
      //
      TH1D *hl1sfr = new TH1D("hsys3a_l1sfr",";p_{T} (GeV);"
			       "[L1_{dt,altsf}/L1_{dt,nomsf}](#LT#rho#GT+20%) /"
			       "[L1_{dt,altsf}/L1_{dt,nomsf}](#LT#rho#GT) ",
			       npt, &ptbins[0]);
      TH1D *hl1sf20 = new TH1D("hl1sfp20",";p_{T} (GeV);"
			       "[L1_{dt,altsf}/L1_{dt,nomsf}](#LT#rho#GT+20%)",
			       npt, &ptbins[0]);
      TH1D *hl1dt20 = new TH1D("hl1dt20",";p_{T} (GeV);"
			       "(p_{T,corr}+L1offdt(#LT#rho#GT+20%)) /"
			       "(p_{T,corr}+L1offdt(#LT#rho#GT))",
			       npt, &ptbins[0]);
      //
      TH1D *hl1r = new TH1D("hsys3b_l1r",";p_{T} (GeV);"
			    "[(L1_{dt,alt}/L1_{mc,alt})/"
			    "(L1_{dt,nom}/L1_{mc,nom})] "
			    "(#LT#rho#GT+20%) / (#LT#rho#GT)",
			    npt, &ptbins[0]);
      TH1D *hl120 = new TH1D("hl120",";p_{T} (GeV);"
			     "[(L1_{dt,alt}/L1_{mc,alt})](#LT#rho#GT+20%) /"
			     "[(L1_{dt,nom}/L1_{mc,nom})](#LT#rho#GT)",
			     npt, &ptbins[0]);
      TH1D *hl1d20 = new TH1D("hl1d20",";p_{T} (GeV);"
			      "[(L1_{dt,alt}/L1_{dt,nom})](#LT#rho#GT+20%)",
			      npt, &ptbins[0]);

      if(rp_debug) cout << "calculate reference R_uncl (hruncl)" << endl<<flush;
      TH1D *hruncl = new TH1D("hruncl",";mode;Runcl",4,0.5,4.5);
      if (CorLevel=="L1L2L3Res") { // R_uncl
	double ptref = 15.;
	int jref = h2pteta->GetXaxis()->FindBin(ptref);
	jref = max(1,min(jref, h2pteta->GetXaxis()->GetNbins()));
	TH1D *hpteta = h2pteta->ProjectionY("hpteta",jref,jref);

	double suml2l3(0), suml2l3res(0), sumw5(0), sumw8(0);
	for (int ieta = 1; ieta != hpteta->GetNbinsX()+1; ++ieta) {

	  double abseta = hpteta->GetBinCenter(ieta);
	  double w = hpteta->GetBinContent(ieta);
	  sumw8 += w;
	  if (abseta<5.191) {
	    sumw5 += w;
	    double l2l3p = 1./getJEC(mcjec,+abseta,ptref);
	    double l2l3m = 1./getJEC(mcjec,-abseta,ptref);
	    double l2l3 = 0.5*(l2l3p+l2l3m);
	    suml2l3 += w*l2l3;
	    double l2l3resp = 1./getJEC(jec,+abseta,ptref);
	    double l2l3resm = 1./getJEC(jec,-abseta,ptref);
	    double l2l3res = 0.5*(l2l3resp+l2l3resm);
	    suml2l3res += w*l2l3res;
	  }
	} // for ieta
	hruncl->SetBinContent(1, sumw5/sumw8);
	hruncl->GetXaxis()->SetBinLabel(1,"|#eta|<5.2");
	hruncl->SetBinContent(2, suml2l3/sumw5);
	hruncl->GetXaxis()->SetBinLabel(2,"L2L3");
	hruncl->SetBinContent(3, suml2l3res/sumw5);
	hruncl->GetXaxis()->SetBinLabel(3,"L2L3Res");
	hruncl->SetBinContent(4, (suml2l3/sumw5)*(suml2l3res/sumw5));
	hruncl->GetXaxis()->SetBinLabel(4,"DATA");
      } // R_uncl

      if(rp_debug) cout << "create reference JES bands" << endl << flush;
      for (int ipt = 1; ipt != herr->GetNbinsX()+1; ++ipt) {

	double pt = herr->GetBinCenter(ipt);

	//int jpt = h2eta->GetXaxis()->FindBin(pt);
	int jpt = h2pteta->GetXaxis()->FindBin(pt);
	// avoid underflow bin for pT<30 GeV
	//jpt = max(1,min(jpt, h2eta->GetXaxis()->GetNbins()));
	jpt = max(1,min(jpt, h2pteta->GetXaxis()->GetNbins()));

	// Take Z+jet eta,pT distribution for correctly averaging JEC
	//TH1D *heta = h2eta->ProjectionY(Form("heta_%d",ipt), jpt, jpt);
	TH1D *heta = h2pteta->ProjectionY(Form("heta_%d",ipt), jpt, jpt);
	const int ieta2 = heta->FindBin(eta2);
	const int ieta1 = heta->FindBin(eta1);
	const int intw = heta->Integral(ieta1,ieta2-1);
	const int neta = 2*(ieta2 - ieta1);

	// Loop over fine eta bins to average JEC and uncertainties
	double sumval(0), sumvall2l3res(0), sumerr2(0), sumw(0);
	double sumrun1(0), sumjes(0);
	double sumvall1flat(0), sumvall1pt(0), sumvall1mc(0), sumvall1rho(0);
	double sumvall1dt(0), sumvall1c(0), sumvall1s(0), sumvall1rc(0);
	double sumvall1rcdt(0), sumvall1rcmc(0);
	double sumvall2c(0), sumvall2s(0);
	double sumerr2_pt(0), sumerr2_hcal(0), sumerr2_ecal(0), sumerr2_pu(0);
	double sumerr2_glu(0), sumerr2_uds(0), sumerr2_cha(0), sumerr2_bot(0);
	double sumerr2_zjt(0), sumerr2_gjt(0), sumerr2_qcd(0);
	double sumerr2_noflv(0), sumerr2_ref(0), sumerr2_ref1(0);

	// New additions
	double sumvalrcdnom(0),sumvalrcmnom(0), sumvalrcdalt(0),sumvalrcmalt(0);
	double sumvall1dnomsf(0), sumvall1daltsf(0), sumvall1mnomsf(0);
	double sumvall1dnom(0),sumvall1mnom(0), sumvall1dalt(0),sumvall1malt(0);
	double sumvall1dnomsf20(0), sumvall1daltsf20(0);
	double sumvall1dnom20(0),sumvall1mnom20(0);
	double sumvall1dalt20(0),sumvall1malt20(0);


	for (int jeta = 0; jeta != neta; ++jeta) {

	  assert(eta2 > eta1);
	  assert(eta1 >= 0);

	  // average over both plus and minus sides
	  int keta = ieta1 + (jeta % (neta/2));
	  double eta = (jeta<neta/2 ? -1 : +1) * heta->GetBinCenter(keta);
	  double w = (intw ? heta->GetBinContent(keta) / (2*intw) : 1./neta);
	  if (!(w<=1)) {
	    cout << "pt =" << pt << " w="<<w << endl << flush;
	    assert(w<=1);
	  }

	  // JEC central values first
	  double val = 1./getJEC(jec,eta,pt);

	  // Special combined JECs
	  if (epoch=="2018ABCD") {
	    assert(jeca); assert(jecb);assert(jecc);assert(jecd);
	    assert(fabs(jecwa+jecwb+jecwc+jecwd-1)<1e-4);
	    
	    double vala = 1./getJEC(jeca,eta,pt);
	    double valb = 1./getJEC(jecb,eta,pt);
	    double valc = 1./getJEC(jecc,eta,pt);
	    double vald = 1./getJEC(jecd,eta,pt);

	    val = jecwa*vala +jecwb*valb +jecwc*valc +jecwd*vald;
	  }
	  else if (epoch=="2017BCDEF") {
	    assert(jecb);assert(jecc);assert(jecd);assert(jece);assert(jecf);
	    assert(fabs(jecwb+jecwc+jecwd+jecwe+jecwf-1)<1e-4);
	    
	    double valb = 1./getJEC(jecb,eta,pt);
	    double valc = 1./getJEC(jecc,eta,pt);
	    double vald = 1./getJEC(jecd,eta,pt);
	    double vale = 1./getJEC(jece,eta,pt);
	    double valf = 1./getJEC(jecf,eta,pt);

	    val = jecwb*valb +jecwc*valc +jecwd*vald +jecwe*vale +jecwf*valf;
	  }
	  else if (epoch=="2016BCDEF") {
	    assert(jecb); assert(jece);
	    assert(fabs(jecwb+jecwe-1)<1e-4);
	    
	    double valb = 1./getJEC(jecb,eta,pt);
	    double vale = 1./getJEC(jece,eta,pt);

	    val = jecwb*valb +jecwe*vale;
	    //val = jecwb*valb; // alternative version
	  }
	  else if (epoch=="2016GH") {
	    assert(jecg);
	    
	    double valg = 1./getJEC(jecg,eta,pt);
	    assert(fabs(jecwg-1)<1e-4);

	    val = jecwg*valg;
	  }
	  else if (isRun2) {

	    // combination now done in minitools/mergeL2L3ResTextFiles.C
	    double val = 1./getJEC(jec,eta,pt);
	    //double sum(0), sumw(0);
	    //for (unsigned int i = 0; i != vjec.size(); ++i) {
	    //double w = vjec[i].first;
	    //sum += w/getJEC(vjec[i].second,eta,pt);
	    //sumw += w;
	    //} // for i in vjec
	    //val = (sum / sumw);
	  }

	  // reference JEC
	  double jesrun1 = 1./getJEC(jecrun1,eta,pt);

	  // old JEC
	  double jes = 1./getJEC(jecold,eta,pt);

          //if(rp_debug) cout << "ABC/ABCD special treatment done" << endl;
	  double vall2l3res = val; // For closure test 
	  //double vall2l3res = jes; // For closure test  (V5M2)
          //if(CorLevel=="L1L2L3Res") val = 1; //to get proper bands during closure test
	  if(CorLevel=="L1L2L3Res") val /= jes; //to get proper bands during closure test
	  if(CorLevel=="L1L2L3Res") jesrun1 /= jes;  //to get proper bands during closure test
	  if(CorLevel=="L1L2L3Res") jes /= jes;  //to get proper bands during closure test

	  sumvall2l3res += w*vall2l3res;
	  sumrun1 += w*jesrun1;
	  sumval += w*val;
	  sumjes += w*jes;
	  sumw += w; // sum weights only once
	  
	  // Various JEC for systematics
	  sumvall1flat += w/     getJEC(jecl1flat, eta,pt);
	  sumvall1rcdt += w/     getJEC(jecl1rcdt, eta,pt);
	  sumvall1pt   += w/     getJEC(jecl1pt,   eta,pt);
	  sumvall1mc   += w/     getJEC(jecl1mc,   eta,pt);
	  sumvall1rcmc += w/     getJEC(jecl1rcmc, eta,pt);
	  sumvall1rho  += w/     getJEC(jecl1pt,   eta,pt,gRho+1);
	  sumvall1dt   += w/     getJEC(jecl1dt,   eta,pt);
	  sumvall1c    += w/     getJEC(jecl1c,    eta,pt);
	  sumvall1s    += w/     getJEC(jecl1s,    eta,pt);
	  sumvall1rc   += w/max( getJEC(jecl1rc,   eta,pt),0.0001);
	  sumvall2c    += w/     getJEC(jecl2c,    eta,pt);
	  sumvall2s    += w/     getJEC(jecl2s,    eta,pt);

	  // New additions
	  sumvalrcdnom += w/     getJEC(jrcdnom,    eta,pt,gRhoDT);
	  sumvalrcmnom += w/     getJEC(jrcmnom,    eta,pt,gRhoMC);
	  sumvalrcdalt += w/     getJEC(jrcdalt,    eta,pt,gRhoDT);
	  sumvalrcmalt += w/     getJEC(jrcmalt,    eta,pt,gRhoMC);
	  //
	  sumvall1dnomsf += w/   getJEC(jl1dnomsf, eta,pt,gRhoDT);
	  sumvall1daltsf += w/   getJEC(jl1daltsf, eta,pt,gRhoDT);
	  sumvall1mnomsf += w/   getJEC(jl1mnomsf, eta,pt,gRhoMC);
	  //
	  sumvall1dnom += w/     getJEC(jl1dnom, eta,pt,gRhoDT);
	  sumvall1mnom += w/     getJEC(jl1mnom, eta,pt,gRhoMC);
	  sumvall1dalt += w/     getJEC(jl1dalt, eta,pt,gRhoDT);
	  sumvall1malt += w/     getJEC(jl1malt, eta,pt,gRhoMC);
	  //
	  sumvall1dnomsf20 += w/ getJEC(jl1dnomsf, eta,pt,gRhoDT*1.2);
	  sumvall1daltsf20 += w/ getJEC(jl1daltsf, eta,pt,gRhoDT*1.2);
	  //
	  sumvall1dnom20 += w/   getJEC(jl1dnom, eta,pt,gRhoDT*1.2);
	  sumvall1mnom20 += w/   getJEC(jl1mnom, eta,pt,gRhoMC*1.2);
	  sumvall1dalt20 += w/   getJEC(jl1dalt, eta,pt,gRhoDT*1.2);
	  sumvall1malt20 += w/   getJEC(jl1malt, eta,pt,gRhoMC*1.2);

	  // JEC uncertainties
	  unc->setJetEta(eta);
	  unc->setJetPt(pt);
	  double err = unc->getUncertainty(true);
	  sumerr2 += w*err*err; // sum squared weights only once

	  unc_ref->setJetEta(eta);
	  unc_ref->setJetPt(pt);
	  double err_ref = unc_ref->getUncertainty(true);
	  sumerr2_ref += w*err_ref*err_ref;

	  unc_ref1->setJetEta(eta);
	  unc_ref1->setJetPt(pt);
	  double err_ref1 = unc_ref1->getUncertainty(true);
	  sumerr2_ref1 += w*err_ref1*err_ref1;

	  unc_pt->setJetEta(eta);
	  unc_pt->setJetPt(pt);
	  double err_pt = unc_pt->getUncertainty(true);
	  sumerr2_pt += w*err_pt*err_pt;

	  unc_hcal->setJetEta(eta);
	  unc_hcal->setJetPt(pt);
	  double err_hcal = unc_hcal->getUncertainty(true);
	  sumerr2_hcal += w*err_hcal*err_hcal;

	  unc_ecal->setJetEta(eta);
	  unc_ecal->setJetPt(pt);
	  double err_ecal = unc_ecal->getUncertainty(true);
	  sumerr2_ecal += w*err_ecal*err_ecal;

	  unc_glu->setJetEta(eta);
	  unc_glu->setJetPt(pt);
	  double err_glu = unc_glu->getUncertainty(true);
	  sumerr2_glu += w*err_glu*err_glu;

	  unc_uds->setJetEta(eta);
	  unc_uds->setJetPt(pt);
	  double err_uds = unc_uds->getUncertainty(true);
	  sumerr2_uds += w*err_uds*err_uds;

	  unc_cha->setJetEta(eta);
	  unc_cha->setJetPt(pt);
	  double err_cha = unc_cha->getUncertainty(true);
	  sumerr2_cha += w*err_cha*err_cha;

	  unc_bot->setJetEta(eta);
	  unc_bot->setJetPt(pt);
	  double err_bot = unc_bot->getUncertainty(true);
	  sumerr2_bot += w*err_bot*err_bot;

	  unc_zjt->setJetEta(eta);
	  unc_zjt->setJetPt(pt);
	  double err_zjt = unc_zjt->getUncertainty(true);
	  sumerr2_zjt += w*err_zjt*err_zjt;

	  unc_gjt->setJetEta(eta);
	  unc_gjt->setJetPt(pt);
	  double err_gjt = unc_gjt->getUncertainty(true);
	  sumerr2_gjt += w*err_gjt*err_gjt;

	  unc_qcd->setJetEta(eta);
	  unc_qcd->setJetPt(pt);
	  double err_qcd = unc_qcd->getUncertainty(true);
	  sumerr2_qcd += w*err_qcd*err_qcd;

	  unc_pu->setJetEta(eta);
	  unc_pu->setJetPt(pt);
	  double err_pu = unc_pu->getUncertainty(true);
	  sumerr2_pu += w*err_pu*err_pu;

	  unc_noflv->setJetEta(eta);
	  unc_noflv->setJetPt(pt);
	  double err_noflv = unc_noflv->getUncertainty(true);
	  sumerr2_noflv += w*err_noflv*err_noflv;
	} // for jeta

	// normalize by total weight for correct average
	double vall2l3res = sumvall2l3res / sumw;
	double val = sumval / sumw;
        //if(CorLevel=="L1L2L3Res") val = 1.0; //to get proper bands during closure test -- should already be 1.0 ... TESTING

	// normalize uncertainties (quadratic instead of linear addition)
	double err = sqrt(sumerr2 / sumw);
	double err_ref = sqrt(sumerr2_ref / sumw);
	double err_ref1 = sqrt(sumerr2_ref1 / sumw);
	double err_pt = sqrt(sumerr2_pt / sumw);
	double err_hcal = sqrt(sumerr2_hcal / sumw);
	double err_ecal = sqrt(sumerr2_ecal / sumw);
	double err_glu = sqrt(sumerr2_glu / sumw);
	double err_uds = sqrt(sumerr2_uds / sumw);
	double err_cha = sqrt(sumerr2_cha / sumw);
	double err_bot = sqrt(sumerr2_bot / sumw);
	double err_zjt = sqrt(sumerr2_zjt / sumw);
	double err_gjt = sqrt(sumerr2_gjt / sumw);
	double err_qcd = sqrt(sumerr2_qcd / sumw);
	double err_pu = sqrt(sumerr2_pu / sumw);
	double err_noflv = sqrt(sumerr2_noflv / sumw);

	// center uncertainties around JEC central value
	herr->SetBinContent(ipt, val);
	herr_l2l3res->SetBinContent(ipt, vall2l3res);
	herr_ref->SetBinContent(ipt, val);
	herr_spr->SetBinContent(ipt, val);
	herr_pu->SetBinContent(ipt, 1);
	herr_noflv->SetBinContent(ipt, val);
	herr_mpf->SetBinContent(ipt, val);
	herr_pt->SetBinContent(ipt, val);

	herr->SetBinError(ipt, val*err);
	herr_l2l3res->SetBinError(ipt, vall2l3res*err_ref);
	herr_ref->SetBinError(ipt, val*err_ref);
	herr_spr->SetBinError(ipt,val*sqrt(err_hcal*err_hcal
					   + err_ecal*err_ecal));
	herr_glu->SetBinError(ipt, val*err_glu);
	herr_uds->SetBinError(ipt, val*err_uds);
	herr_cha->SetBinError(ipt, val*err_cha);
	herr_bot->SetBinError(ipt, val*err_bot);
	herr_zjt->SetBinError(ipt, val*err_zjt);
	herr_gjt->SetBinError(ipt, val*err_gjt);
	herr_qcd->SetBinError(ipt, val*err_qcd);
	herr_pu->SetBinError(ipt, val*err_pu);
	herr_noflv->SetBinError(ipt, val*err_noflv);
	herr_mpf->SetBinError(ipt, val*err_pt);
	herr_pt->SetBinError(ipt, val*sqrt(err_pt*err_pt+err_pu*err_pu));

	double run1 = (sumrun1 / sumw);
        //if(CorLevel=="L1L2L3Res") run1 = 1.0;  //to get proper bands during closure test
	hrun1->SetBinContent(ipt, run1);
	hrun1->SetBinError(ipt, run1*err_ref1);

	double jes = (sumjes / sumw);
        //if(CorLevel=="L1L2L3Res") jes = 1.0;  //to get proper bands during closure test
	hjes->SetBinContent(ipt, jes);
	double l1pt = (sumvall1pt / sumw);
	double l1mc = (sumvall1mc / sumw);
	double l1rho = (sumvall1rho / sumw);
	double l1flat = (sumvall1flat / sumw);
	double l1bias = l1pt / l1flat;
	double l1diff = l1pt / l1mc;
	double l1drho = l1pt / l1rho;
	hl1bias->SetBinContent(ipt, l1bias);
	hl1bias->SetBinError(ipt, 0.5*fabs(1-l1bias));
	hl1diff->SetBinContent(ipt, l1diff);
	hl1diff->SetBinError(ipt, 0.5*fabs(1-l1diff));
	hl1drho->SetBinContent(ipt, l1drho);
	hl1drho->SetBinError(ipt, 0.5*fabs(1-l1drho));
	// New ratios: data/MC(simple), complex/simple, RC/simple
	// Do these for JEC=1/JES, so easier to track impact on L3Res
	// (how much JES would have been different with another JEC)
	double l1dt = (sumvall1dt / sumw);
	double l1c = (sumvall1c / sumw);
	double l1s = (sumvall1s / sumw);
	double l1rc = (sumvall1rc / sumw);
	double l1rcdt = (sumvall1rcdt / sumw);
	double l1rcmc = (sumvall1rcmc / sumw);
	hl1dtomc->SetBinContent(ipt, (1./l1dt) / (1./l1mc));
	hrcdtomc->SetBinContent(ipt, (1./l1rcdt) / (1./l1rcmc));
	hl1cos->SetBinContent(ipt, (1./l1c) / (1./l1s));
	hl1rcos->SetBinContent(ipt, (1./l1rc) / (1./l1s));
	double l2c = (sumvall2c / sumw);
	double l2s = (sumvall2s / sumw);
	hl2cos->SetBinContent(ipt, (1./l2c) / (1./l2s));

	// New additions
	double rcdnom = (sumvalrcdnom / sumw);
	double rcmnom = (sumvalrcmnom / sumw);
	double rcdalt = (sumvalrcdalt / sumw);
	double rcmalt = (sumvalrcmalt / sumw);
	hrcmpf->SetBinContent(ipt, (rcdalt/rcmalt) / (rcdnom/rcmnom)); //(1)
	hrcmpfdt->SetBinContent(ipt, rcdalt/rcdnom); //xtra
	hrcmpfmc->SetBinContent(ipt, rcmalt/rcmnom); //xtra
	hrcdt->SetBinContent(ipt, rcdnom); //xtra
	hrcmc->SetBinContent(ipt, rcmnom); //xtra
	//
	// JES*pTcorr = pTcorr + offset => offset = pTcorr*(JES-1)
	double l1dnomsf = (sumvall1dnomsf / sumw);
	double l1daltsf = (sumvall1daltsf / sumw);
	double l1mnomsf = (sumvall1mnomsf / sumw);
	double offl1dnomsf = pt*(l1dnomsf-1);
	hl1sf->SetBinContent(ipt, l1daltsf/l1dnomsf); //(2a)
	hl1sfnom->SetBinContent(ipt, (l1dnomsf-1)/(l1mnomsf-1)); //xtra
	hl1sfalt->SetBinContent(ipt, (l1daltsf-1)/(l1mnomsf-1)); //xtra
	hl1offdt->SetBinContent(ipt, offl1dnomsf); //xtra
	hl1dt->SetBinContent(ipt, l1dnomsf); //xtra
	//
	double l1dnom = (sumvall1dnom / sumw);
	double l1mnom = (sumvall1mnom / sumw);
	double l1dalt = (sumvall1dalt / sumw);
	double l1malt = (sumvall1malt / sumw);
	hl1->SetBinContent(ipt, (l1dalt/l1malt)/(l1dnom/l1mnom)); //(2b)
	hl1d->SetBinContent(ipt, l1dalt/l1dnom); //xtra
	hl1m->SetBinContent(ipt, l1malt/l1mnom); //xtra
	hl1dnom->SetBinContent(ipt, l1dnom); //xtra
	hl1mnom->SetBinContent(ipt, l1mnom); //xtra
	//
	double l1dnomsf20 = (sumvall1dnomsf20 / sumw);
	double l1daltsf20 = (sumvall1daltsf20 / sumw);
	hl1sfr->SetBinContent(ipt, (l1daltsf20/l1dnomsf20)
			      / (l1daltsf/l1dnomsf) ); //(3a)
	hl1sf20->SetBinContent(ipt, l1daltsf20/l1dnomsf20); //xtra
	hl1dt20->SetBinContent(ipt, l1dnomsf20); //xtra
	//
	double l1dnom20 = (sumvall1dnom20 / sumw);
	double l1mnom20 = (sumvall1mnom20 / sumw);
	double l1dalt20 = (sumvall1dalt20 / sumw);
	double l1malt20 = (sumvall1malt20 / sumw);
	hl1r->SetBinContent(ipt, ((l1dalt20/l1malt20)/(l1dnom20/l1mnom20)) /
			    ((l1dalt/l1malt)/(l1dnom/l1mnom))); //(3b)
	hl120->SetBinContent(ipt, (l1dalt20/l1malt20)/(l1dnom20/l1mnom20)); //xtra
	hl1d20->SetBinContent(ipt, l1dalt20/l1dnom20); //xtra
      } // ipt
      if(rp_debug) cout << "done creating reference JES bands" << endl;


      dout2->cd();

      herr->SetMarkerSize(0);
      herr->SetFillStyle(1001);
      herr->SetFillColorAlpha(kYellow+1,0.5);
      herr->Write();

      herr_l2l3res->SetMarkerSize(0);
      herr_l2l3res->SetFillStyle(1001);
      herr_l2l3res->SetFillColorAlpha(kYellow+1,0.5);
      herr_l2l3res->Write();

      herr_ref->SetMarkerSize(0);
      herr_ref->SetFillStyle(1001);
      herr_ref->SetFillColorAlpha(kYellow+1,0.5);
      herr_ref->Write();

      herr_spr->SetMarkerSize(0);
      herr_spr->SetFillStyle(1001);
      herr_spr->SetFillColorAlpha(kYellow,0.5);
      herr_spr->Write();

      herr_glu->SetMarkerSize(0);
      herr_glu->SetFillStyle(1001);
      herr_glu->SetFillColorAlpha(kYellow+1,0.5);
      herr_glu->Write();

      herr_uds->SetMarkerSize(0);
      herr_uds->SetFillStyle(1001);
      herr_uds->SetFillColorAlpha(kYellow+1,0.5);
      herr_uds->Write();

      herr_cha->SetMarkerSize(0);
      herr_cha->SetFillStyle(1001);
      herr_cha->SetFillColorAlpha(kYellow+1,0.5);
      herr_cha->Write();

      herr_bot->SetMarkerSize(0);
      herr_bot->SetFillStyle(1001);
      herr_bot->SetFillColorAlpha(kYellow+1,0.5);
      herr_bot->Write();

      herr_zjt->SetMarkerSize(0);
      herr_zjt->SetFillStyle(1001);
      herr_zjt->SetFillColorAlpha(kYellow+1,0.5);
      herr_zjt->Write();

      herr_gjt->SetMarkerSize(0);
      herr_gjt->SetFillStyle(1001);
      herr_gjt->SetFillColorAlpha(kYellow+1,0.5);
      herr_gjt->Write();

      herr_qcd->SetMarkerSize(0);
      herr_qcd->SetFillStyle(1001);
      herr_qcd->SetFillColorAlpha(kYellow+1,0.5);
      herr_qcd->Write();

      herr_pu->SetMarkerSize(0);
      herr_pu->SetFillStyle(1001);
      herr_pu->SetFillColorAlpha(kYellow+1,0.5);
      herr_pu->Write();

      herr_noflv->SetMarkerSize(0);
      herr_noflv->SetFillStyle(1001);
      herr_noflv->SetFillColorAlpha(kYellow+1,0.5);
      herr_noflv->Write();

      herr_mpf->SetMarkerSize(0);
      herr_mpf->SetFillStyle(1001);
      herr_mpf->SetFillColorAlpha(kYellow+1,0.5);
      herr_mpf->Write();

      herr_pt->SetMarkerSize(0);
      herr_pt->SetFillStyle(1001);
      herr_pt->SetFillColorAlpha(kYellow+1,0.5);
      herr_pt->Write();

      hrun1->SetMarkerSize(0);
      hrun1->SetFillStyle(1001);
      hrun1->SetFillColorAlpha(kCyan+1,0.5);
      hrun1->Write();

      hjes->Write();
      hl1bias->Write();
      hl1diff->Write();
      hl1drho->Write();
      //
      hl1dtomc->Write();
      hrcdtomc->Write();
      hl1cos->Write();
      hl1rcos->Write();
      hl2cos->Write();

      // New additions
      hrcmpf->Write(); // (1)
      hrcmpfdt->Write(); //xtra
      hrcmpfmc->Write(); //xtra
      hrcdt->Write(); //xtra
      hrcmc->Write(); //xtra
      //
      hl1sf->Write(); // (2a)
      hl1sfnom->Write(); //xtra
      hl1sfalt->Write(); //xtra
      hl1offdt->Write(); //xtra
      hl1dt->Write(); //xtra
      //
      hl1->Write(); // (2b)
      hl1d->Write(); //xtra
      hl1m->Write(); //xtra
      hl1dnom->Write(); //xtra
      hl1mnom->Write(); //xtra
      //
      hl1sfr->Write(); // (3a)
      hl1sf20->Write(); //xtra
      hl1dt20->Write(); //xtra
      //
      hl1r->Write(); // (3b)
      hl120->Write(); //xtra
      hl1d20->Write(); //xtra
      //
      hruncl->Write();
    } // ieta
  } // JEC+sys

  fout->Close();

} // reprocess

// Helper function to retrieve FactorizedJetCorrector
FactorizedJetCorrector *getFJC(string l1, string l2, string res, string path) {

  // Set default jet algo
  if (l1!="" && !(TString(l1.c_str()).Contains("_AK")))
    l1 += "_AK4PFchs";
  if (l2!="" && !(TString(l2.c_str()).Contains("_AK")))
    l2 += "_AK4PFchs";
  if (res!="" && !(TString(res.c_str()).Contains("_AK")))
    res += "_AK4PFchs";

  //string path = "CondFormats/JetMETObjects/data";
  if (path=="") path = "CondFormats/JetMETObjects/data";
  const char *cd = path.c_str();
  const char *cl1 = l1.c_str();
  const char *cl2 = l2.c_str();
  const char *cres = res.c_str();
  string s("");

  vector<JetCorrectorParameters> v;
  if (l1!=""){
    s = Form("%s/%s.txt",cd,cl1);
    cout << s << endl << flush;
    JetCorrectorParameters *pl1 = new JetCorrectorParameters(s);
    v.push_back(*pl1);
  }
  if (l2!="") {
    s = Form("%s/%s.txt",cd,cl2);
    cout << s << endl << flush;
    JetCorrectorParameters *pl2 = new JetCorrectorParameters(s);
    v.push_back(*pl2);
  }
  if (res!="") {
    s = Form("%s/%s.txt",cd,cres);
    cout << s << endl << flush;
    JetCorrectorParameters *pres = new JetCorrectorParameters(s);
    v.push_back(*pres);
  }
  FactorizedJetCorrector *jec = new FactorizedJetCorrector(v);

  return jec;
} // getJFC

// Invert JEC to get differences as a function of pTcorr
// Helper functions to find JEC for corrected pt
void setEtaPtRho(FactorizedJetCorrector *jec, double eta, double pt,
                 double rho){

  assert(jec);
  jec->setJetEta(eta);
  jec->setJetPt(pt);
  jec->setRho(rho);
  jec->setJetA(0.50265);

  return;
}

Double_t funcCorrPt(Double_t *x, Double_t *p) {
  
  double eta = p[0];
  double pt = x[0];
  double rho = p[1];
  setEtaPtRho(_thejec, eta, pt, rho);

  return (_thejec->getCorrection() * pt);
}

//const double _rhoDE = 19.21; // EOY17 DE jt450
double getJEC(FactorizedJetCorrector *jec, double eta, double ptcorr,
	      double rho) {

  // iterate to solve ptreco for given ptcorr
  _thejec = jec;
  if (!fCorrPt) fCorrPt = new TF1("fCorrPt",funcCorrPt,5,6500,2);
  fCorrPt->SetParameters(eta, rho);
  // Find ptreco that gives pTreco*JEC = pTcorr
  double ptreco = fCorrPt->GetX(ptcorr,5,6500);

  setEtaPtRho(jec, eta, ptreco, rho);

  return (jec->getCorrection());
} // getEtaPtE
