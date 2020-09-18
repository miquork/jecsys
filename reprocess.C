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
//const double gRho = 15.36; // 2017-06-02 for 2016 data
//const double gRho = 19.21; // EOY17 DE jt450
// These values from Z(mm)+jet UL17BCDEF |eta|<1.3 at close to 202.5 GeV
// Maximilian: [MC,Data]_Rho_CHS_a30_eta_00_13_L1L2Res at 198 GeV (alpha<1.0)
// Sami:       h_Zpt_rho_alpha100_eta00_13 at 208 GeV (alpha<0.3)
const double gRhoDT = 21.39; // +/-0.04 (21.39+/-0.05) (18.41+/-0.04)
const double gRhoMC = 20.60; // +/-0.03 (20.60+/-0.09) (20.67+/-0.03)
const double gRho = 21.00; // +/-0.40 Average of gRhoDT and gRhoMC
const bool _dcsonly = false;
const bool rp_debug = true; // verbose messages

// Appy mass corrections to Z+jet
// Update 20200325: start moving these to globalFitSyst.C as syst. eigenvectors
// Need to correct Zee and Zmm here to combine them, but photon could be later
bool useFixedFit = true; // with minitools/drawZmass.C (def:true)
double fitUncMin = 0.00000; // Some bug if unc<0?
bool correctZmmMass = true; // pT with mumu mass (def:true)
bool correctZeeMass = true; // pT with ee mass (def:true)
bool correctGamMass = true; //!!UL17_V3 false, pT with ee mass at 2*pT
bool correctUncert = false;  // ll mass uncertainty => globalFitSyst.C
//
// Which binning to use for multijets ("leading", "recoil" or "ptave")
string multijetMode = "ptave"; // check also multijetModeS in globalFitSyst.C
bool correctMultijetLeading = false;//true; // correct for diference to JER-hybrid (UL17 true, UL18 false)
bool patchMultijetPtAveMCMJBPt20 = false;
//bool patchMultijetStatOverPt = 1684; // RMSfit/sqrt(N); zero for none
//
// Which binning to use for PF composition ("tp","pt", "both" or "none" (off))
string pfMode = "tp"; // "both" not supported at the moment

bool confirmWarnings=false;//true; //if active, deficiencies in input files are patched after confirmation via pushing "any key to continue"
bool confirmedGamJet=false; // to avoid too many confirmations
//need to adjust corlevel in multijet.C as well!!! Not synced automatically right now
string CorLevel="L1L2Res"; // same for gamma+jet and Z+jet on RunD
//string CorLevel="L1L2L3Res"; // For Fall17 only had L1L
//string CorLevel="L1L2";
//"L1L2": "MCTruth corrections" applied, in reality L1Data/MC,L2
//"L1L2Res": MCTruth + L2Res applied
//"L1L2L3Res": MCTruth + L2L3Res applied; reference JES central values set to 1.0 (affects pliotting as well)

bool correctGamScale = true; 
double valueGamScale = 1.01;//1.;



// Settings for cleaned up global fit
/////////////////////////////////////
// Minimum pTcut for gamma+jet
// 85: 141.1, 105: 124.7/72, 135: 122.4/70
//double fpmpfptmin(175);//105);//85);//30);//100.); // photon+jet MPF
//double fpbalptmin(175);//105);//30);//105);//85);//30);//100.); // photon+jet pTbal
double fzeeptmin(15.);   // Zee+jet both methods
double fzmmptmin(15.);//30.);   // Zmm+jet both methods
double fzptmin(15.);//40.);//30.);  // Z+jet both methods
// Additional cuts to Z+jet MPF / balance methods
//double fzmpfptmin(30.);   // Z+jet MPF
//double fzbalptmin(30.);//100);//30.);   // Z+jet pTbal

// multijet minimum and maximum pT
// (high for now until FSR working for low pT)
// 49: 228.6/78, 64: 155.7/76, 84: 135.8/74, 114: 124.7/72, 153: 120.7/70
// 20200330: 49:238.3/84, 64:164.6/82, 84:143.5/80, 114:131.4/78,
// ...20200330: 153:125.4/76 196:108.5/74 245:107.7/72
// ...20200419 ptave: 300:
double fpfjetptmin(15.);
double fpfjetptmax(2116.);
double fincjetptmin(21);//18)15.); // UL17 21
double fincjetptmax(2116.);//4037.);
double fhadwptamin(35);//35.);//30.); // UL17 30.
double fhadwptamax(200);//200.); // UL17 200
double fhadwptbmin(35);//40.);//30.); // UL17 30.
double fhadwptbmax(175);//200.); // UL17 200
// 200.5/141 => 188.0/139 => 172.9/135 => 175.0/137 => 207.5/141 => 201.4/141
// 174.1/138
double fmultijetptmin(507);//300);//153);//114); // 114 in UL17, 430 UL18
double fmultijetptmax(2640);//2116); // 2116 in UL17
double fmultijetptmax2(1890);//2366);//1890); // 1890 in UL17
double fmultijetmjbptmin(1032); // avoid FSR bias in UL18 with 1032
// 124.5/113 => 776.2/149 => 164.8/125 => 137.3/117 => 105.7/111 => 117.1/119

//for fine etabins deactivate ptbal
double fdijetmpfptmin(30);
double fdijetbalptmin(30.);
double fdijetptmax(1500.);

/*
// Settings for maximally inclusive glofal fit
/////////////////////////////////////////////
// (use for checking data/data, MC/MC ratios before fit)
// Minimum pTcut for gamma+jet
double fpmpfptmin(30.);//100. // photon+jet MPF
double fpbalptmin(30.);//100; // photon+jet pTbal
double fzeeptmin(30.); // Zee+jet
double fzmmptmin(30.); // Zmm+jet
// Additional cuts to Z+jet MPF / balance methods
double fzmpfptmin(30.); // Z+jet MPF
double fzbalptmin(30.); // Z+jet pTbal
*/

// Maximum pTcut for samples (to avoid bins with too large uncertainty)
double fpmpfptmax(1500.); // photon+jet MPF
double fpbalptmax(700);//1500);//700.);  // photon+jet pTbal
double fzeeptmax(700.);   // Zee+jet
double fzmmptmax(700.);   // Zmm+jet
double fzptmax(700.);   // Z+jet
// Additional cuts to Z+jet MPF / balance methods
//double fzmpfptmax(700);//500.);  // Z+jet MPF
//double fzbalptmax(500);//300.);  // Z+jet pTbal

// Regular settings for fits with Z+jet
double fpmpfptmin(230.); double fpbalptmin(230.); double fzllmpfptmin(15);/*40);30);*/ double fzmpfptmin(15);/*40);30);*/ double fzllbalptmin(80);/*15);*//*80);*//*105.);80.);*/ double fzbalptmin(105);/*15);*//*105.);*//*30);*/ double fzllmpfptmax(700); double fzmpfptmax(700); double fzllbalptmax(500); double fzbalptmax(500); // UL17
//double fpmpfptmin(100.); double fpbalptmin(100.); double fzllmpfptmin(40); double fzllbalptmin(100); double fzllmpfptmax(500); double fzllbalptmax(300); // EOY17 settings (incl. effective 40 GeV from bug and gamma+jet min pT)

// Special setting for "multijet only" fit 
//double fzllmpfptmin(130); double fzllbalptmin(130); double fzllmpfptmax(175); double fzllbalptmax(175); // option B
//double fzllmpfptmin(175); double fzllbalptmin(175); double fzllmpfptmax(230); double fzllbalptmax(230); // option A



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

  // Monday 19 Feb 2018
  //https://indico.cern.ch/event/706518/#2-l2res-with-dijet-2017-data
  //https://indico.cern.ch/event/706518/contributions/2899252/attachments/1602656/2549097/L2Res_Fall17_17Nov_V5V6_dijet_nominal.tar.gz
  map<string,const char*> fdj_files;
  fdj_files["B"] = "B";
  fdj_files["C"] = "C";
  fdj_files["D"] = "D";
  fdj_files["E"] = "E";
  fdj_files["F"] = "F";
  fdj_files["BCDEF"] = "BCDEF";
  //  TFile *fdj = new TFile(Form("rootfiles/L2Res_Fall17_17Nov_V5V6/Run%s/JEC_L2_Dijet_AK4PFchs_pythia8.root",fdj_files[epoch]),"READ");
  //patch... using closure files...
  //TFile *fdj = new TFile(Form("rootfiles/L2Res_GlobalFitInput_JERSF_V2_Fall17_17Nov2017_V31_DATA/Run%s_17Nov17_2017_JERnominalV2/output/JEC_L2_Dijet_AK4PFchs_pythia8_eta_0_13.root",fdj_files[epoch]),"READ"); //only RunBCDEF, nominal JERSF_V2
  TFile *fdj = new TFile(Form("rootfiles/L2Res_GlobalFitInput_JERSF_V2_Fall17_17Nov2017_V31_DATA/Run%s_17Nov17_2017_JERnominalV2/output/JEC_L2_Dijet_AK4PFchs_pythia8_eta_0_13.root","BCDEF"),"READ"); //UL2017-v1
    assert(fdj && !fdj->IsZombie());
    //TFile *fdj2 = new TFile(Form("rootfiles/L2Res_GlobalFitInput_JERSF_V2_Fall17_17Nov2017_V31_DATA/Run%s_17Nov17_2017_JERnominalV2/output/JEC_L2_Dijet_AK4PFchs_pythia8.root",fdj_files[epoch]),"READ"); //only RunBCDEF, nominal JERSF_V2
    TFile *fdj2 = new TFile(Form("rootfiles/L2Res_GlobalFitInput_JERSF_V2_Fall17_17Nov2017_V31_DATA/Run%s_17Nov17_2017_JERnominalV2/output/JEC_L2_Dijet_AK4PFchs_pythia8.root","BCDEF"),"READ"); // UL2017-v1
    assert(fdj2 && !fdj2->IsZombie());

  // Anastasia Karavdina, 2016 Legacy re-reco (8 Dec 2017) :
  // https://indico.cern.ch/event/682570/
  //TFile *fdj = new TFile("rootfiles/JEC_L2_Dijet_AK4PFchs_pythia8_07122017hcalCleaning_wideEtabins.root"); // BCDEFGH only?
  //assert(fdj && !fdj->IsZombie());

  //TFile *fdj2 =0;
  //TFile *fdj2 = new TFile("rootfiles/JEC_L2_Dijet_AK4PFchs_pythia8_07122017hcalCleaning.root"); // narrow bins to complement above wide ones
  //assert(fdj2 && !fdj2->IsZombie());

  if(CorLevel=="L1L2L3Res"){
    //Monday 22 Oct 2018, Fall17_17Nov2017_V31_DATA
    //  https://indico.cern.ch/event/767001/#51-closure-test-for-fall17_17n
    //    rootfiles/L2Res_GlobalFitInput_JERSF_V2_RunBCDEF_Fall17_17Nov2017_V31_DATA/RunBCDEF_17Nov17_2017_JERnominalV2/output/
    fdj = new TFile(Form("rootfiles/L2Res_GlobalFitInput_JERSF_V2_Fall17_17Nov2017_V31_DATA/Run%s_17Nov17_2017_JERnominalV2/output/JEC_L2_Dijet_AK4PFchs_pythia8_eta_0_13.root",fdj_files[epoch]),"READ"); //only RunBCDEF, nominal JERSF_V2
    assert(fdj && !fdj->IsZombie());
    fdj2 = new TFile(Form("rootfiles/L2Res_GlobalFitInput_JERSF_V2_Fall17_17Nov2017_V31_DATA/Run%s_17Nov17_2017_JERnominalV2/output/JEC_L2_Dijet_AK4PFchs_pythia8.root",fdj_files[epoch]),"READ"); //only RunBCDEF, nominal JERSF_V2
    assert(fdj2 && !fdj2->IsZombie());
  }




  if(CorLevel!="L1L2"&&CorLevel!="L1L2L3Res"){
    cout << Form("Dijet files are only available for Autumn18-V17 (L1L2+L2L3Res) narrow/wide bins. L1L2L3Res used for L1L2Res as placeholder. Please confirm by pushing any key.") << endl;
    if(confirmWarnings)cin.ignore();
  }


  // Andrey Popov, April, 2017 (Feb03_L2ResV2)
  // https://indico.cern.ch/event/634367/
  // last traditional style input
  map<string,const char*> fm_files;
  //dummy files, temp
//   fm_files["A"] = "0428_Run2016All"; // also update multijet.C
//   fm_files["B"] = "0428_Run2016All"; // also update multijet.C
//   fm_files["C"] = "0428_Run2016All"; // also update multijet.C
//   fm_files["D"] = "0428_Run2016All"; // also update multijet.C
//   fm_files["ABC"] = "0428_Run2016All"; // also update multijet.C
//   fm_files["ABCD"] = "0428_Run2016All"; // also update multijet.C
//   TFile *fmj = new TFile(Form("rootfiles/multijet_2017%s.root",
//   		      fm_files[epoch]),"READ");

  // Andrey Popov, March 19, 2018
  // https://indico.cern.ch/event/713034/#4-residuals-with-multijet-2016
  // `NEW MULTIJET INPUT`
  //fm_files["BCD"] = "BCD";
  //fm_files["EF"] = "EFearly";
  //fm_files["G"] = "FlateGH"; //X
  //fm_files["H"] = "FlateGH"; //X
  //fm_files["GH"] = "FlateGH";
  //fm_files["BCDEFGH"] = "All";
  //TFile *fmj = new TFile(Form("rootfiles/multijet_180319_2016%s.root",
  //		      fm_files[epoch]),"READ");


  
  // add Minsuk's new multijet files in reprocess.C and in multijet.C
  // they are found in rootfiles/multijet_20190911_JEC_Autunm18_V17_JER_Autumn18_V7
  //  
  //Subject: 	RE: Multijet combination file format
  //Date: 	Wed, 11 Sep 2019 15:53:41 +0200
  fm_files["B"] = "B"; // also update multijet.C
  fm_files["C"] = "C"; // also update multijet.C
  fm_files["BC"] = "BC";
  fm_files["D"] = "D"; // also update multijet.C
  fm_files["E"] = "E"; // also update multijet.C
  fm_files["DE"] = "DE"; // also update multijet.C
  fm_files["F"] = "F"; // also update multijet.C
  fm_files["BCDEF"] = "BCDEF"; // also update multijet.C
  fm_files["2018ABCD"] = "ABCD";
  fm_files["2018A"] = "A";
  fm_files["2018B"] = "B";
  fm_files["2018C"] = "C";
  fm_files["2018D"] = "D";
  TFile *fmj(0);
  bool isUL18 = (epoch=="2018ABCD" || epoch=="2018A" || 
		 epoch=="2018B" || epoch=="2018C" || epoch=="2018D");
  //if (epoch=="2018ABCD" ||
  //epoch=="2018A" || epoch=="2018B" || epoch=="2018C" || epoch=="2018D") {
  if (isUL18) {
    //fmj = new TFile(Form("rootfiles/multijet_Rebin2_20200526_UL2017%s_SimpleL1V4JRV2_jecV4_jerV2.root","BCDEF"),"READ"); // placeholder
    fmj = new TFile(Form("rootfiles/multijet_UL2018%s_jecDTV3MCV2_jerUL17V2-2.root",fm_files[epoch]),"READ");
  }
  else { // UL17
    //fmj = new TFile(Form("rootfiles/multijet_Rebin2_20200526_UL2017%s_SimpleL1V4JRV2Sigma80MB_jecV4_jerV2.root",fm_files[epoch]),"READ"); // 80 mb - NG
    fmj = new TFile(Form("rootfiles/multijet_Rebin2_20200526_UL2017%s_SimpleL1V4JRV2_jecV4_jerV2.root",fm_files[epoch]),"READ"); // 69.2 mb
  }
  assert(fmj && !fmj->IsZombie());
  //TFile *fmj =0;

  //TFile *fij = new TFile("rootfiles/drawDeltaJEC_17UL_V2M4res_cp2_all_v2.root","READ");
  //TFile *fij = new TFile("rootfiles/drawDeltaJEC_17UL_V2M4res_cp2_all_v3.root","READ"); // Update JER SF V2 + MC truth JER per IOV
  map<string,const char*> fij_eras;
  fij_eras["2018A"] = "A";
  fij_eras["2018B"] = "B";
  fij_eras["2018C"] = "C";
  fij_eras["2018D"] = "D";
  fij_eras["2018ABCD"] = "";
  TFile *fij = (isUL18 ? 
		new TFile("rootfiles/drawDeltaJEC_18UL_JECV3.root","READ") :
		new TFile("rootfiles/drawDeltaJEC_17UL_V2M4res_cp2_all_v4.root","READ")); // Update JER SF V2 + MC truth JER per IOV (V2M5 input)
  assert(fij && !fij->IsZombie());

  map<string,const char*> fpf_files;
  fpf_files["2018A"] = "A";
  fpf_files["2018B"] = "B";
  fpf_files["2018C"] = "C";
  fpf_files["2018D"] = "D";
  fpf_files["2018ABCD"] = "ABCD";
  TFile *fpfdt(0), *fpfmc(0);
  if (isUL18) {
    fpfdt = new TFile(Form("rootfiles/output-DATA-2b-UL18V2V3-%s.root",
			   fpf_files[epoch]),"READ");
    fpfmc = new TFile(Form("rootfiles/output-MCNU-2b-UL18V2V3-%s.root",
			   fpf_files[epoch]),"READ");
  }
  else { // UL17
    fpfdt = new TFile(Form("rootfiles/output-DATA-2b-UL17V4_%s.root",
			   epoch.c_str()),"READ");
    // MC80NU 80mb chi2=189.2/137, MCNU 69.2mb=218.137/137
    //TFile *fpfmc = new TFile(Form("rootfiles/output-MC80NU-2b-UL17V4_%s.root",
    fpfmc = new TFile(Form("rootfiles/output-MCNU-2b-UL17V4_%s.root",
			   epoch.c_str()),"READ");
  }

  assert((fpfdt && !fpfdt->IsZombie()) || pfMode=="none");
  assert((fpfmc && !fpfmc->IsZombie()) || pfMode=="none");

  map<string,const char*> fw_files;
  fw_files["2018A"] = "A";
  fw_files["2018B"] = "B";
  fw_files["2018C"] = "C";
  fw_files["2018D"] = "D";
  fw_files["2018ABCD"] = "ABCD";
  TFile *fw = new TFile(Form("rootfiles/hadW%s%s.root",
			     isUL18 ? "UL18" : "",
			     isUL18 ? fw_files[epoch] : ""),
			"READ");
  //	       new TFile("rootfiles/hadW.root","READ"));
  assert(fw && !fw->IsZombie());

  map<string,const char*> fp_files;
  fp_files["B"] = "B";
  fp_files["C"] = "C";
  fp_files["BC"] = "BC";
  fp_files["D"] = "D";
  fp_files["E"] = "E";
  fp_files["DE"] = "DE";
  fp_files["F"] = "F";
  fp_files["BCDEF"] = "BCDEF";
  fp_files["2018A"] = "A";
  fp_files["2018B"] = "B";
  fp_files["2018C"] = "C";
  fp_files["2018D"] = "D";
  fp_files["2018ABCD"] = "ABCD";
  TFile *fp(0);
  if (isUL18) {
    fp = new TFile(Form("../JERCProtoLab/Summer19UL18/L3Residual_gamma/Gjet_combinationfile_only_L2Res_%s_only_L2Res.root",fp_files[epoch]),"READ");
  }
  else { // UL17
    fp = new TFile(Form("rootfiles/2020-04-02/SimpleL1_only_L2Res/Gjet_combinationfile_SimpleL1_only_L2Res_%s_SimpleL1_only_L2Res.root",fp_files[epoch]),"READ");
  //fp = new TFile(Form("rootfiles/2020-04-02/ComplexL1_only_L2Res/Gjet_combinationfile_ComplexL1_only_L2Res_%s_ComplexL1_only_L2Res.root",fp_files[epoch]),"READ");
  }

  if(CorLevel=="L1L2L3Res"){
    assert(false); // re-enable later
    // Hugues Lattaud, 2017 V31 closure, Nov 12, 2018
    // https://indico.cern.ch/event/772590/#6-l3res-gammajet-closure
    fp_files["B"] = "B";
    fp_files["C"] = "C";
    fp_files["D"] = "D";
    fp_files["E"] = "E";
    fp_files["F"] = "F";
    fp_files["BCDEF"] = "BCDEF";
    fp = new TFile(Form("rootfiles/Gjet_combinationfile_17_Nov_V31_L2L3res_%s.root", fp_files[epoch]),"READ");
    // can not update, yet, due to incompatible binning
    //    fp = new TFile(Form("rootfiles/Gjet_combinationfile_17_Nov_V24_L2L3res_%s.root", fp_files[epoch]),"READ"); //residual residual input https://indico.cern.ch/event/758051/#21-closure-test-with-gammajets
  }

  assert(fp && !fp->IsZombie());

  map<string,const char*> fz_files;
  fz_files["B"] = "B";
  fz_files["C"] = "C";
  fz_files["D"] = "D";
  fz_files["E"] = "E";
  fz_files["F"] = "F";
  fz_files["BCDEF"] = "BCDEF";
  fz_files["2018A"] = "A";
  fz_files["2018B"] = "B";
  fz_files["2018C"] = "C";
  fz_files["2018D"] = "D";
  fz_files["2018ABCD"] = "ABCD";
  TFile *fzmm(0), *fzee(0);
  if (isUL18) {
    fzmm = new TFile("../JERCProtoLab/Summer19UL18/L3Residual_Z/JEC_Combination_Zmm/ZJetCombination_Zmm_DYJets_Madgraph_12Nov2019_Summer19UL18_JECV3_L1L2Res.root","READ"); // UL18-v1
    fzee = new TFile("../JERCProtoLab/Summer19UL18/L3Residual_Z/JEC_Combination_Zee/ZJetCombination_Zee_DYJets_Madgraph_12Nov2019_Summer19UL18_JECV3_L1L2Res.root","READ"); // UL18-v1
  }
  else { // UL17
    fzmm = new TFile("../JERCProtoLab/Summer19UL17/L3Residual_Z/JEC_Combination_Zmm/ZJetCombination_Zmm_DYJets_Madgraph_09Aug2019_Summer19UL17_JECV4_L1L2Res.root","READ"); // V4 as above, V5 new and worse
    fzee = new TFile("../JERCProtoLab/Summer19UL17/L3Residual_Z/JEC_Combination_Zee/ZJetCombination_Zee_DYJets_Madgraph_09Aug2019_Summer19UL17_JECV4_L1L2Res.root","READ"); // V4 as above, V5 new and worse
  }
  assert(fzmm && !fzmm->IsZombie());
  assert(fzee && !fzee->IsZombie());

  //TFile *fz = new TFile("rootfiles/jme_bplusZ_merged_vX.root","READ");
  //TFile *fz = new TFile("rootfiles/jme_bplusZ_merged_v22.root","READ");
  //TFile *fz = new TFile("rootfiles/jme_bplusZ_merged_vRawChsMET.root","READ");
  //TFile *fz = new TFile("rootfiles/jme_bplusZ_merged_v24.root","READ");
  //TFile *fz = new TFile("rootfiles/jme_bplusZ_merged_v25.root","READ");
  TFile *fz = new TFile("rootfiles/jme_bplusZ_merged_v26_DYJets.root","READ");
  assert(fz && !fz->IsZombie());
  // JES correction for zjet
  TFile *fzjes = new TFile("rootfiles/jecdataBCDEF_V4.root","READ");
  assert(fzjes && !fzjes->IsZombie());
  TH1D *hzjes = (TH1D*)fzjes->Get("ratio/eta00-13/herr_ref"); assert(hzjes);

  // UL2017-v1
  if (isUL18) {
    fzmm->cd(Form("Run2018%s",fz_files[epoch])); fzmm = (TFile*)gDirectory;
    fzee->cd(Form("Run2018%s",fz_files[epoch])); fzee = (TFile*)gDirectory;
  }
  else { // UL17
    fzmm->cd(Form("Run2017%s",fz_files[epoch])); fzmm = (TFile*)gDirectory;
    fzee->cd(Form("Run2017%s",fz_files[epoch])); fzee = (TFile*)gDirectory;
  }

  if(CorLevel=="L1L2L3Res"){
    // Monday 22 Oct V31 closure
    // https://indico.cern.ch/event/767001/#52-closure-test-for-fall17_17n
    fzmm = new TFile(Form("rootfiles/zjet_combination_Fall17_JECV31_Zmm_%s_2018-10-26.root",fz_files[epoch]),"READ"); // https://indico.cern.ch/event/759977/#35-closure-test-for-fall17_17n
    fzee = new TFile(Form("rootfiles/zjet_combination_Fall17_JECV31_Zee_%s_2018-10-26.root",fz_files[epoch]),"READ"); // for October 2018 release
 }

  // Link to Z mass files (same as above now) and histograms
  // This is used for scaling Z mass back to 1 for Zee and Zmm
  TFile *fmz = fz;
  TFile *fmzee = fzee;
  TFile *fmzmm = fzmm;
  assert(fmzee && !fmzee->IsZombie());
  assert(fmzmm && !fmzmm->IsZombie());

  string sr = (epoch=="L4" ? "eta_00_24" : "eta_00_13");
  const char *cr = sr.c_str();
  /*
  TH1D *hmzee = (TH1D*)fmzee->Get(Form("Ratio_ZMass_CHS_a30_%s_L1L2L3",cr));
  TH1D *hmzmm = (TH1D*)fmzmm->Get(Form("Ratio_ZMass_CHS_a30_%s_L1L2L3",cr));
  assert(hmzee);
  assert(hmzmm);
  TH1D *hmzee_dt = (TH1D*)fmzee->Get(Form("Data_ZMass_CHS_a30_%s_L1L2L3",cr));
  TH1D *hmzmm_dt = (TH1D*)fmzmm->Get(Form("Data_ZMass_CHS_a30_%s_L1L2L3",cr));
  assert(hmzee_dt);
  assert(hmzmm_dt);
  TH1D *hmzee_mc = (TH1D*)fmzee->Get(Form("MC_ZMass_CHS_a30_%s_L1L2L3",cr));
  TH1D *hmzmm_mc = (TH1D*)fmzmm->Get(Form("MC_ZMass_CHS_a30_%s_L1L2L3",cr));
  assert(hmzee_mc);
  assert(hmzmm_mc);
  */
  TH1D *hmzee = (TH1D*)fmzee->Get(Form("Ratio_ZMass_CHS_a30_%s_L1L2Res",cr));
  TH1D *hmzmm = (TH1D*)fmzmm->Get(Form("Ratio_ZMass_CHS_a30_%s_L1L2Res",cr));
  assert(hmzee);
  assert(hmzmm);
  TH1D *hmzee_dt = (TH1D*)fmzee->Get(Form("Data_ZMass_CHS_a30_%s_L1L2Res",cr));
  TH1D *hmzmm_dt = (TH1D*)fmzmm->Get(Form("Data_ZMass_CHS_a30_%s_L1L2Res",cr));
  assert(hmzee_dt);
  assert(hmzmm_dt);
  TH1D *hmzee_mc = (TH1D*)fmzee->Get(Form("MC_ZMass_CHS_a30_%s_L1L2Res",cr));
  TH1D *hmzmm_mc = (TH1D*)fmzmm->Get(Form("MC_ZMass_CHS_a30_%s_L1L2Res",cr));
  assert(hmzee_mc);
  assert(hmzmm_mc);

  TH1D *hmzee1 = (TH1D*)fmzee->Get(Form("Ratio_ZMass_CHS_a100_%s_L1L2Res",cr));
  TH1D *hmzmm1 = (TH1D*)fmzmm->Get(Form("Ratio_ZMass_CHS_a100_%s_L1L2Res",cr));
  assert(hmzee1);
  assert(hmzmm1);
  TH1D *hmzee1_dt =(TH1D*)fmzee->Get(Form("Data_ZMass_CHS_a100_%s_L1L2Res",cr));
  TH1D *hmzmm1_dt =(TH1D*)fmzmm->Get(Form("Data_ZMass_CHS_a100_%s_L1L2Res",cr));
  assert(hmzee1_dt);
  assert(hmzmm1_dt);
  TH1D *hmzee1_mc = (TH1D*)fmzee->Get(Form("MC_ZMass_CHS_a100_%s_L1L2Res",cr));
  TH1D *hmzmm1_mc = (TH1D*)fmzmm->Get(Form("MC_ZMass_CHS_a100_%s_L1L2Res",cr));
  assert(hmzee1_mc);
  assert(hmzmm1_mc);

  // v24 has _eta_00_13, v25 not
  // v26 missing Z mass again
  /*
  TH2D *h2mz = (TH2D*)fmz->Get("data/eta_00_13/h_Zpt_mZ_alpha30");
  TH2D *h2mz_dt = (TH2D*)fmz->Get("data/eta_00_13/h_Zpt_mZ_alpha30");
  TH2D *h2mz_mc = (TH2D*)fmz->Get("mc/eta_00_13/h_Zpt_mZ_alpha30");
  assert(h2mz);
  assert(h2mz_dt);
  assert(h2mz_mc);
  TH1D *hmz_dt = h2mz_dt->ProfileX("h2mz_dt")->ProjectionX("hmz_dt");
  TH1D *hmz_mc = h2mz_mc->ProfileX("h2mz_mc")->ProjectionX("hmz_mc");
  TH1D *hmz = h2mz->ProfileX("h2mz")->ProjectionX("hmz");
  hmz->Divide(hmz_mc);
  */

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
	f1mzee->SetParameters(1.00153, 0.00214, -0.00012);
	f1ezee->SetParameters( +5.5e-08,  +1.5e-07, +3.01e-07,
			       +2.35e-08, -7.73e-08, -4.34e-08);
	// UL18 Run2018ABCD Zee from above
	f1mgam->SetParameters(1.00153, 0.00214, -0.00012);
	f1egam->SetParameters( +5.5e-08,  +1.5e-07, +3.01e-07,
			       +2.35e-08, -7.73e-08, -4.34e-08);
      }
      else {
	assert(false); // CHECK
	// UL17 RunBCDEF fit with minitools/drawZmass.C
	f1mzee->SetParameters(0.99780, 0.00225, 0.00031);
	f1ezee->SetParameters( +6.1e-08, +1.66e-07, +3.27e-07,
			       +2.56e-08,  -8.5e-08, -5.39e-08);
	// UL17 RunBCDEF Zee from above
	f1mgam->SetParameters(0.99780, 0.00225, 0.00031);
	f1egam->SetParameters( +6.1e-08, +1.66e-07, +3.27e-07,
			       +2.56e-08,  -8.5e-08, -5.39e-08);
      }
      // No correction
      //f1mgam->SetParameters(1,0,0);
    }
    else
      hmzee->Fit(f1mzee);

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
	f1mzmm->SetParameters(0.99839, 0.00000, 0.00000);
	f1ezmm->SetParameters(+1.64e-08,        +0,        +0,
			             +0,        +0,        +0);	
      }
      else {
	assert(false); // CHECK
	// UL17 BCDEF fit with minitools/drawZmass.C
	f1mzmm->SetParameters(0.99821, 0.00000, 0.00000);
	f1ezmm->SetParameters(+3.43e-08,        +0,        +0,
			             +0,        +0,        +0);
      }
    }
    else
      hmzmm->Fit(f1mzmm);
  }


  // \END copy-paste from minitools/drawZmass.C


  //TF1 *f2 = new TF1("f2","[p0]*pow(x,2)+[p1]",30,3000); // bug in V4?
  TF1 *f2 = new TF1("f2","[0]*pow(x,2)+[1]",30,3000);
  f2->SetParameters(0,0);

  if (correctMultijetLeading && multijetMode=="leading") {
    f2->SetParameters(1.105e-07, 0.03063); // V4

    // Fits done in minitools/systMultijet.C on UL17BCDEF
    //f2 Chi2/NDF = 15.2 / 19
    //TF1 *f2 = new TF1("f2","[p0]*pow(x,2)+[p1]",30,3000);
    f2->SetParameters(1.462e-07, -0.001596); // 20200419-SimpleL1 MPF
  }

  // Link to Z+jet 2D distribution for JEC calculations
  // This is used for correctly averaging JEC and its uncertainty
  // for the wide eta bins used in global fit combinations
  // => need to update for Run 2
  // On 28 Feb 2014, at 14:00, from Dominik Haitz
  // Re: pT~200 GeV
  TFile *feta = new TFile("rootfiles/th2d_jeteta_zpt.root","READ");
  assert(feta && !feta->IsZombie());

  TH2D *h2eta = (TH2D*)feta->Get("data"); assert(h2eta);


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


  // Map variable names used in different files to common scheme
  map<string, map<string, const char*> > rename;

  //rename["pfjet"]["mpfchs1"] = "tp";
  //rename["pfjet"]["ptchs"] = ""; // "pt"

  // Note: non-CHS results are not available for all samples, so not used
  rename["dijet"]["mpfchs"] = "mpfchs";
  rename["dijet"]["mpfchs1"] = "mpfchs";
  rename["dijet"]["ptchs"] = "ptchs";

  // Andrey and Run I
//   rename["multijet"]["ratio"] = "Data";//"Ratio"; => PATCH
//   rename["multijet"]["data"] = "Data";
//   rename["multijet"]["mc"] = "MC";
//   rename["multijet"]["crecoil"] = "CRecoil";
//   rename["multijet"]["mpfchs"] = "MPF";
//   rename["multijet"]["mpfchs1"] = "MPF";
//   rename["multijet"]["ptchs"] = "MJB";

// Minsuk in Run II
  rename["multijet"]["ratio"] = "Data";// => PATCH
  rename["multijet"]["data"] = "Data";
  rename["multijet"]["mc"] = "MC";
  rename["multijet"]["ratiocrecoil"] = "CRecoil_L1L2Res";
  rename["multijet"]["datacrecoil"] = "CRecoil_L1L2Res";
  rename["multijet"]["mccrecoil"] = "CRecoil_MG";//P8";
  // Switch to leading jet pT binning instead or recoil pT binning
  if (multijetMode=="leading") {
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
    }
    else {
      rename["multijet"]["ratiompfchs1"] = "MPF_recoil_L1L2Res";
      rename["multijet"]["ratioptchs"] = "MJB_recoil_L1L2Res";
      rename["multijet"]["datampfchs1"] = "MPF_recoil_L1L2Res";
      rename["multijet"]["dataptchs"] = "MJB_recoil_L1L2Res";
      rename["multijet"]["mcmpfchs1"] = "MPF_recoil_MG";
      rename["multijet"]["mcptchs"] = "MJB_recoil_MG";
    }
  }
  else if (multijetMode=="ptave") {
    if (isUL18) {
      rename["multijet"]["ratiocrecoil"] = "CRecoil_ptave_L1L2res";
      rename["multijet"]["datacrecoil"] = "CRecoil_ptave_L1L2res";
      rename["multijet"]["mccrecoil"] = "CRecoil_ptave_MG";
      //
      rename["multijet"]["ratiompfchs1"] = "MPF_ptave_L1L2res_metType2";
      rename["multijet"]["ratioptchs"] = "MJB_ptave_L1L2res";
      rename["multijet"]["datampfchs1"] = "MPF_ptave_L1L2res_metType2";
      rename["multijet"]["dataptchs"] = "MJB_ptave_L1L2res";
      rename["multijet"]["mcmpfchs1"] = "MPF_ptave_MG_metType2";
      rename["multijet"]["mcptchs"] = "MJB_ptave_MG";
    }
    else {
      //rename["multijet"]["ratiocrecoil"] = "CRecoil_ptave";
      //rename["multijet"]["datacrecoil"] = "CRecoil_ptave";
      rename["multijet"]["ratiocrecoil"] = "CRecoil_ptave_L1L2Res";
      rename["multijet"]["datacrecoil"] = "CRecoil_ptave_L1L2Res";
      rename["multijet"]["mccrecoil"] = "CRecoil_ptave_MG";
      //
      //rename["multijet"]["ratiompfchs1"] = "MPF_ptave_L1L2Res";
      rename["multijet"]["ratiompfchs1"] = "MPF_ptave_L1L2Res_metType2";
      rename["multijet"]["ratioptchs"] = "MJB_ptave_L1L2Res";
      //rename["multijet"]["datampfchs1"] = "MPF_ptave_L1L2Res";
      rename["multijet"]["datampfchs1"] = "MPF_ptave_L1L2Res_metType2";
      rename["multijet"]["dataptchs"] = "MJB_ptave_L1L2Res";
      //rename["multijet"]["mcmpfchs1"] = "MPF_ptave_MG";
      rename["multijet"]["mcmpfchs1"] = "MPF_ptave_MG_metType2";
      rename["multijet"]["mcptchs"] = "MJB_ptave_MG";
    }
  }
  else
    assert(false);

  if (isUL18) {
    rename["incjet"]["mpfchs1"] = "jet_Run18UL%s_det";
    rename["incjet"]["ptchs"] = "jet_Run18UL%s_det"; // copy
  }
  else {
    //rename["incjet"]["ratio"] = "Run17UL";
    //rename["incjet"]["mpfchs1"] = "jet_Run17UL%s_fwd";
    //rename["incjet"]["mpfchs1"] = "jet_Run17UL%s_fwd2"; // JER SF V2 + JER per IOV
    rename["incjet"]["mpfchs1"] = "jet_Run17UL%s_fwd3"; // JER V2 + JES fix
    //rename["incjet"]["ptchs"] = "jet_Run17UL%s_dag";
    //rename["incjet"]["ptchs"] = "jet_Run17UL%s_fwd2"; // copy
    rename["incjet"]["ptchs"] = "jet_Run17UL%s_fwd3"; // copy
  }

  rename["hadw"]["mpfchs1"] = "mass_ptave";//ptboth";
  rename["hadw"]["ptchs"] = "mass_ptboth";//ptave";
  rename["hadw"]["counts"] = "nevents_ptave";//ptboth";

  rename["gamjet"]["ratio"] = "";
  rename["gamjet"]["data"] = "_DATA"; 
  rename["gamjet"]["mc"] = "_MC"; 
  rename["gamjet"]["mpfchs"] = "resp_MPFchs";
  rename["gamjet"]["mpfchs1"] = "resp_MPFchs"; 
  rename["gamjet"]["ptchs"] = "resp_PtBalchs"; 
  rename["gamjet"]["counts"] = "RawNEvents_data_vs_pt";

  rename["zeejet"]["ratio"] = "Ratio";
  rename["zeejet"]["data"] = "Data";
  rename["zeejet"]["mc"] = "MC";
  rename["zeejet"]["mpfchs"] = "MPF-notypeI_CHS";
  rename["zeejet"]["mpfchs1"] = "MPF_CHS"; 
  rename["zeejet"]["ptchs"] = "PtBal_CHS"; 
  rename["zeejet"]["counts"] = "RawNEvents_CHS";
  rename["zeejet"]["chf"] = "PFjetCHF_CHS";
  rename["zeejet"]["nef"] = "PFjetPF_CHS";
  rename["zeejet"]["nhf"] = "PFjetNHF_CHS";
  rename["zeejet"]["cef"] = "PFjetEF_CHS";
  rename["zeejet"]["muf"] = "PFjetMF_CHS";

  rename["zmmjet"]["ratio"] = "Ratio";
  rename["zmmjet"]["data"] = "Data";
  rename["zmmjet"]["mc"] = "MC";
  rename["zmmjet"]["mpfchs"] = "MPF-notypeI_CHS";
  //rename["zmmjet"]["mpfchs1"] = "MPF-notypeI_CHS"; //TEMP TEMP TEMP!!!
  rename["zmmjet"]["mpfchs1"] = "MPF_CHS"; 
  rename["zmmjet"]["ptchs"] = "PtBal_CHS"; 
  rename["zmmjet"]["counts"] = "RawNEvents_CHS";
  rename["zmmjet"]["chf"] = "PFjetCHF_CHS";
  rename["zmmjet"]["nef"] = "PFjetPF_CHS";
  rename["zmmjet"]["nhf"] = "PFjetNHF_CHS";
  rename["zmmjet"]["cef"] = "PFjetEF_CHS";
  rename["zmmjet"]["muf"] = "PFjetMF_CHS";

  // Results from Sami's Z+b analysis
  rename["zjet"]["ratio"] = "data"; // missing => PATCH
  rename["zjet"]["data"] = "data";
  rename["zjet"]["mc"] = "mc";
  rename["zjet"]["mpfchs"] = "mpfchs";
  rename["zjet"]["mpfchs1"] = "mpfchs"; 
  rename["zjet"]["ptchs"] = "ptchs"; 
  rename["zjet"]["counts"] = "statistics_mpfchs";
  rename["zjet"]["chf"] = "chHEF";
  rename["zjet"]["nef"] = "neEmEF";
  rename["zjet"]["nhf"] = "neHEF";
  rename["zjet"]["cef"] = "chEmEF";
  rename["zjet"]["muf"] = "muEF";


  // color and style codes
  map<string, map<string, int> > style;

  style["pfjet"]["chf"] = kFullCircle;
  style["pfjet"]["nhf"] = kFullDiamond;
  style["pfjet"]["nef"] = kFullSquare;
  style["pfjet"]["cef"] = kFullDiamond;
  style["pfjet"]["muf"] = kFullDiamond;
  style["pfjet_mc"]["chf"] = kOpenCircle;
  style["pfjet_mc"]["nhf"] = kOpenDiamond;
  style["pfjet_mc"]["nef"] = kOpenSquare;
  style["pfjet_mc"]["cef"] = kOpenDiamond;
  style["pfjet_mc"]["muf"] = kOpenDiamond;
  style["dijet"]["mpfchs1"] = kFullDiamond;//kFullCircle;
  style["dijet"]["ptchs"] = kOpenDiamond;//kOpenSquare;
  style["incjet"]["mpfchs1"] = kFullDiamond;
  style["incjet"]["ptchs"] = kOpenDiamond;
  style["hadw"]["mpfchs1"] = kFullCircle;
  style["hadw"]["ptchs"] = kOpenCircle;
  style["multijet"]["mpfchs1"] = kFullTriangleUp;
  style["multijet"]["ptchs"] = kOpenTriangleUp;
  style["gamjet"]["mpfchs1"] = kFullSquare;
  style["gamjet"]["ptchs"] = kOpenSquare;
  style["zeejet"]["mpfchs1"] = kFullCircle;
  style["zeejet"]["ptchs"] = kOpenCircle;
  style["zmmjet"]["mpfchs1"] = kFullStar;
  style["zmmjet"]["ptchs"] = kOpenStar;
  style["zlljet"]["mpfchs1"] = kFullDiamond;
  style["zlljet"]["ptchs"] = kOpenDiamond;
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

  map<string, int> color;
  color["pfjet"] = kBlack;
  color["pfjet_chf"] = kRed;
  color["pfjet_nhf"] = kGreen+2;
  color["pfjet_nef"] = kBlue;
  color["pfjet_cef"] = kCyan+1;
  color["pfjet_muf"] = kMagenta+1;
  color["dijet"] = kBlack;
  color["incjet"] = kOrange+2;//kBlack;
  color["hadw"] = kGreen+2;
  color["multijet"] = kBlack;
  color["gamjet"] = kBlue;
  color["zeejet"] = kGreen+2;
  color["zmmjet"] = kRed;
  color["zlljet"] = kMagenta+2;
  color["zlljet_chf"] = kRed;
  color["zlljet_nhf"] = kGreen+2;
  color["zlljet_nef"] = kBlue;
  color["zlljet_cef"] = kCyan+1;
  color["zlljet_muf"] = kMagenta+1;
  color["zjet"] = kRed+1;
  color["zjet_chf"] = kRed;
  color["zjet_nhf"] = kGreen+2;
  color["zjet_nef"] = kBlue;
  color["zjet_cef"] = kCyan+1;
  color["zjet_muf"] = kMagenta+1;

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
  //types.push_back("puf");

  vector<string> sets;
  sets.push_back("dijet");
  sets.push_back("incjet");
  sets.push_back("multijet");
  sets.push_back("gamjet");
  sets.push_back("zeejet");
  sets.push_back("zmmjet");
  sets.push_back("zlljet");
  sets.push_back("zjet");
  sets.push_back("hadw");
  sets.push_back("pfjet");

  vector<pair<double,double> > etas;
  // reference region |eta|<1.3
  if (epoch!="L4") etas.push_back(make_pair<double,double>(0,1.305));
  if (epoch=="L4") etas.push_back(make_pair<double,double>(0,2.4));
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
	  
	  TFile *f = files[s];
	  // pfjet has two files
	  if (!f && s=="pfjet") f = files[s+"_"+d];
	  if (s=="pfjet" && pfMode=="none") continue;

	  // Take pT and MPF from different files for gamma+jet (or not)
	  //if (s=="gamjet" && f==0) {
	  //if (t=="mpfchs" || t=="mpfchs1") f = fp;//fp1;
	  //if (t=="ptchs") f = fp;//fp2;
	  //}
	  assert(f || s=="zlljet" || (s=="pfjet" && d=="ratio"));

	  // C_recoil is only meant for multijets
	  if (t=="crecoil" && s!="multijet") continue;

	  // fractions are only for pfjet, and only fractions are for pfjet
	  // (now adding fractions to zjet)
	  bool isfrac = (t=="chf"||t=="nef"||t=="nhf"||
			 t=="cef"||t=="muf"||t=="puf");
	  if (isfrac && s!="pfjet" && //s!="zjet" && // v26
	      s!="zeejet" && s!="zmmjet" && s!="zlljet") continue; // v25
	  //if (isfrac && s!="pfjet" && s!="zjet") continue; // v25
	  //if (isfrac && s!="pfjet") continue; // v24
	  if (!isfrac && s=="pfjet") continue;

	  for (unsigned int ialpha = 0; ialpha != alphas.size(); ++ialpha) {

	    double alpha = alphas[ialpha];

	    eta1 = etas[ieta].first; // reset to avoid trouble with below
	    eta2 = etas[ieta].second; // reset to avoid trouble with below
	    
	    // Take wide and narrow bins from different file for dijet and Z+jet
	    /*
	      bool narrowBin = ((fabs(eta2-eta1)<0.4 && fabs(eta2)<3.1) ||
	      (fabs(eta1-3.0)<0.1 && fabs(eta2-3.1)<0.1) ||
	      (fabs(eta1-3.139)<0.1 && fabs(eta2-3.489)<0.1) ||
	      fabs(eta1)>3.4);
	    */
	    bool narrowBin = ((fabs(eta2-eta1)<0.4 && 
	    		       !(fabs(eta1-3.0)<0.05 && fabs(eta2-3.2)<0.05)) ||
	    		      fabs(eta1)>3.8);
	    if (narrowBin) {
	    if (s=="dijet")  f = fdj2;
	    //if (s=="zeejet") f = fzee2; // GH only
	    //if (s=="zmmjet") f = fzmm2; // GH only
	    }
	    assert(f || s=="zlljet" || (s=="pfjet"&&d=="ratio"));

	    //if (alpha>0.35 && s!="zjet") continue;
	    if (alpha>0.35 && (s!="zjet" && s!="zeejet" && s!="zmmjet" && s!="zlljet")) continue;

            if (t=="counts" && s!="zmmjet" && s!="zeejet" && s!="gamjet"
		&& s!="zjet" && s!="hadw")
              continue; // counts available only for z+jet and gamjet, so far

	    if (s=="pfjet" && !(fabs(eta1-0)<0.1 && fabs(eta2-1.3)<0.1
				&& fabs(alpha-0.30)<0.01))
	      continue; // barrel only, no specific alpha
	    if (s=="incjet" && !(fabs(eta1-0)<0.1 && fabs(eta2-1.3)<0.1
				 && fabs(alpha-0.30)<0.01)) // && t=="mpfchs1"))
	      continue; // barrel only, no specific pt or alpha
	    if (s=="hadw" && !(fabs(eta1-0)<0.1 && fabs(eta2-1.3)<0.1
			       && fabs(alpha-0.30)<0.01)) // && t=="mpfchs1"))
	      continue; // barrel only, no specific pt or alpha
	    if (s=="multijet" && (!((epoch!="L4" &&
				     fabs(eta1-0)<0.1 && fabs(eta2-1.3)<0.1) ||
				    (epoch=="L4" &&
				     fabs(eta1-0)<0.1 && fabs(eta2-2.4)<0.1))
				  //|| fabs(alpha-0.10)<0.01 
				  //|| (fabs(alpha-0.15)<0.01 && isUL18)
				  //|| (fabs(alpha-0.20)<0.01 && isUL18)
				  || fabs(alpha-0.10)<0.01))
	      continue; // only barrel for multijet balance, pT=(15),(20),30
	    if (s=="gamjet"  && fabs(eta1-3.2)<0.1) { eta1=3.0; eta2=3.2; }

	    // Reconstruct naming scheme used in each of the files
	    // If non-conventional naming schemes, patch here
	    const char *c(0);
	    if (s=="pfjet") {
	      c = Form("Standard/Eta_%02.1f-%02.1f/p%s%s",
                       eta1, eta2, tt, pfMode=="tp" ? "tp" : ""); 
	    }
	    if (s=="dijet") {
	      c = Form("%s/eta_%02.0f_%02.0f/%s_%s_a%1.0f",
                       dd, 10.001*eta1, 10.001*eta2, rename[s][t], ss, 100.*alpha); 
	    } // dijet
	    if (s=="incjet") {
	      //c = Form("jet_Run17UL%s",epoch=="BCDEF" ? "" : epoch.c_str());
	      if (isUL18)
		c = Form(rename[s][t],fij_eras[epoch]);
	      else
		c = Form(rename[s][t],epoch=="BCDEF" ? "" : epoch.c_str());
	    }
	    if (s=="hadw") {
	      c = Form("%s_%s_%s_fitprob02_L1L2L3",dd,rename[s][t],ss);
	    }
	    if (s=="multijet") {
	      //c = Form("%s/Pt%1.0f/%s", rename[s][d], 100.*alpha, rename[s][t]);
	      //if (isUL18)
	      //c = Form("%s/%s", rename[s][d], rename[s][(d+t)]);
	      //else 
	      c = Form("%s/Pt%1.0f/%s", rename[s][d], 100.*alpha, rename[s][(d+t)]);
	    } // multijet
	    
	    //if (s=="gamjet" && t=="counts" && d=="ratio")continue;
	    if (s=="gamjet" && t=="counts") {
	      c = Form("%s%s_a%1.0f_eta%02.0f_%02.0f_%s",
		       rename[s]["mpfchs1"],
		       d=="ratio" ? rename[s]["data"] : rename[s][d],
		       100.*alpha, 10.*eta1, 10.*eta2, rename[s][t]);
	    }
	    else if (s=="gamjet") {
	      c = Form("%s%s_a%1.0f_eta%02.0f_%02.0f",
		       rename[s][t], rename[s][d],
		       100.*alpha, 10.*eta1, 10.*eta2);
	    } // gamjet
	    if (s=="zmmjet" || s=="zeejet") {
	      //c = Form("%s_%s_a%1.0f_eta_%02.0f_%02.0f_L1L2L3",
	      c = Form("%s_%s_a%1.0f_eta_%02.0f_%02.0f_L1L2Res",
	    	       rename[s][d], rename[s][t], 100.*alpha,
	    	       10.01*eta1, 10.01*eta2);
	      if(CorLevel=="L1L2L3Res")
		c = Form("%s_%s_a%1.0f_eta_%02.0f_%02.0f_L1L2L3Res",
			 rename[s][d], rename[s][t], 100.*alpha,
			 10.01*eta1, 10.01*eta2);
	    } // Zll+jet
	    if (s=="zjet") {
	      c = Form("%s/eta_%02.0f_%02.0f/%s_zmmjet_a%1.0f",
	    	       rename[s][d],10*eta1,10*eta2,rename[s][t],100.*alpha);
	      if (isfrac) {
		//c = Form("%s/eta_%02.0f_%02.0f/h_Zpt_%s_alpha%1.0f_eta_%02.0f_%02.0f",rename[s][d],10*eta1,10*eta2,rename[s][t],100.*alpha,10*eta1,10*eta2); // v24
		c = Form("%s/eta_%02.0f_%02.0f/h_Zpt_%s_alpha%1.0f",rename[s][d],10*eta1,10*eta2,rename[s][t],100.*alpha);
	      }
	    }
	    //assert(c || s=="zlljet");

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
	    if (s=="zjet" && d=="ratio" && (t=="mpfchs1"||t=="ptchs")) {
	      TGraphErrors *gd = grs["data"][t][s][ieta][ialpha];
	      TGraphErrors *gm = grs["mc"][t][s][ieta][ialpha];
	      assert(gd);
	      assert(gm);

	      TGraphErrors *g = tools::ratioGraphs(gd, gm);
	      obj = (TObject*)g;
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

	    assert(obj);

            if (t=="counts" && (s=="zmmjet" || s=="zeejet" || s=="gamjet" || s=="zjet" || s=="hadw") ){ // write out counts to jecdata.root (as TH1F)
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

	    // If data stored in TH2D instead of TGrapherrors, patch her
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
              if (g->GetEY()[i]>0.2)  g->RemovePoint(i);
	      // Clean out point outside good ranges
	      else if (s=="gamjet" && t=="mpfchs1" &&
		       (g->GetX()[i]<fpmpfptmin || g->GetX()[i]>fpmpfptmax))
		g->RemovePoint(i);
	      else if (s=="gamjet" && t=="ptchs" &&
		       (g->GetX()[i]<fpbalptmin || g->GetX()[i]>fpbalptmax))
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
			((epoch!="BCDEF" && epoch!="2018ABCD") &&
			 g->GetX()[i]>fmultijetptmax2) ||
			(t=="ptchs" && g->GetX()[i]<fmultijetmjbptmin)))
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
	    if (s=="zjet" && d=="data" && (t=="mpfchs1" || t=="ptchs")) {
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
	    if (s=="hadw" && d=="ratio" && (t=="mpfchs1" || t=="ptchs")) {
	      double k(1);
	      if (t=="mpfchs1") k = 2.;
	      if (t=="ptchs")  k = sqrt(2.);
	      for (int i = 0; i != g->GetN(); ++i) {
		g->SetPointError(i, g->GetEX()[i], g->GetEY()[i]*k);
	      }
	    }

	    // PATCH new Jan15 EGamma corrections for Legacy 2016
	    /*
	    if (s=="gamjet" && (d=="data" || d=="ratio")) {

	      for (int i = 0; i != g->GetN(); ++i) {

		double x = g->GetX()[i];
		double y = g->GetY()[i];
		// Scale pT>450 GeV up by 0.8% based on minitools/drawGamVsGam.C
		if (x>400 && x<500) {
		  g->SetPoint(i, x, y*1.0040);//1.0045);
		}
		if (x>500) {
		  g->SetPoint(i, x, y*1.0080);//1.0090);
		}
	      }
	    } // patch EGamma corrections
	    */

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

	    // patch botched grapin in 20200419-SimpleL1
	    // by cloning alpha15
	    if (patchMultijetPtAveMCMJBPt20 && multijetMode=="ptave" &&
		s=="multijet" && d=="mc" && t=="ptchs" && alpha==0.20) {
	      
	      // The 1.01 shift is close but not enough at low pT
	      //for (int i = 0; i != g->GetN(); ++i) {
	      //assert(g->GetY()[i]<0.05);
	      //g->SetPoint(i, g->GetX()[i], g->GetY()[i]+1.01);
	      //} // for i
	      TGraphErrors *gfixm = grs["mc"][t][s][ieta][1]; assert(gfixm);
	      g = (TGraphErrors*)gfixm->Clone("ptchs_multijet_a20");
	    } // patch botched bin
	    //if (patchMultiJetStatOverPt>0 && s=="multijet") {
	      // approximate RMS with minitools/systMultijet.C:doStat
	      //double rms = (multijetMode=="ptave" ? 0.11 : 0.09);
	      //for (int i = 0; i != g->GetN(); ++i) {
	    //double pt = g->GetX();
	    //	if (pt > patchMultiJetStatOverPt) {
	    //	  g->SetPointError(i, g->GetEX()[i], rms/sqrt(n));
	    //	}
	    //}
	    //} // patchMultiJetStatOverPt

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
		  //if (isUL18) { // PATCH upside down C_recoil
		  //g->SetPoint(i, 0.5*(xd+xm), yd*ym!=0 ? 1./g->GetY()[i] : 0);
		  //}
		}
		else { // mpfchs1, ptchs
		  g->SetPoint(i, 0.5*(xd+xm), ym ? yd / ym : 0.);
		  //g->SetPoint(i, 0.5*(xd+xm), ym ? 0.975 * yd / ym : 0.); // !!PATCH!!
		  //g->SetPoint(i, 0.5*(xd+xm), yd ? ym / yd : 0.); // !!PATCH!!g => did't work as well, although post-fit better
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
 	    if (correctZmmMass && s=="zmmjet" && (d=="data" || d=="ratio")) {
 	      for (int i = 0; i != g->GetN(); ++i) {
		double pt = g->GetX()[i];
 		double ipt = hmzmm->FindBin(pt);
		//double k = max(0.99,min(1.01,hmzmm->GetBinContent(ipt)));
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
	    // Mass corrections for Zee
	    // Limit uncertainty at high pT to 0.5% (about half correction)
	    if (correctZeeMass && s=="zeejet" && (d=="data" || d=="ratio")) {
	      for (int i = 0; i != g->GetN(); ++i) {
		double pt = g->GetX()[i];
		double ipt = hmzee->FindBin(pt);
		//double k = max(0.99,min(1.01,hmzee->GetBinContent(ipt)));
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
	    if (correctGamMass && s=="gamjet" && (d=="data" || d=="ratio")) {
	      for (int i = 0; i != g->GetN(); ++i) {
		//assert(!correctGamScale); // UL17
		double ptgam = g->GetX()[i];
		double pt = 2*ptgam; // fix post Legacy2016
		//double ipt = hmzee->FindBin(pt);
		double ipt = min(hmzee->FindBin(fzmmptmax), hmzee->FindBin(pt));
		double ek = min(0.005,hmzee->GetBinError(ipt));
		//if (useFixedFit) ek = max(fitUncMin,f1ezee->Eval(pt));
		//double k = f1mzee->Eval(pt);
		if (useFixedFit) ek = max(fitUncMin,f1egam->Eval(pt));
		double k = f1mgam->Eval(pt); // vs pTZ
		g->SetPoint(i, g->GetX()[i], g->GetY()[i]*k);
		if (correctUncert)
		  g->SetPointError(i, g->GetEX()[i], 
				   sqrt(pow(g->GetEY()[i]*k,2) + ek*ek));
	      }
	    }

	    // Photon scale correction (from drawGamVsZmm)
	    // No separate unceratinty added, is already in global fit
	    if (correctGamScale && s=="gamjet" && (d=="data" || d=="ratio")) {
	      //assert(!correctGamMass); // UL17
	      for (int i = 0; i != g->GetN(); ++i) {
		g->SetPoint(i, g->GetX()[i], g->GetY()[i]*valueGamScale);
	      }
	    }


	    // Multijet JER-hybrid corrections
	    // 196.5/78 => 190.2/78; flatter at high pT, less MPF bias
	    if (correctMultijetLeading && multijetMode=="leading" &&
		s=="multijet" && (d=="data" || d=="ratio")) {

	      for (int i = 0; i != g->GetN(); ++i) {
		double pt = g->GetX()[i];
		double k0 = f2->Eval(pt); // MJBlead / JER-hybrid - 1 (%)
		double kHybridOverLead = 1./(1+0.01*k0);
		g->SetPoint(i, pt, g->GetY()[i] * kHybridOverLead);
	      }
	    }

	    dout->cd("orig");
	    g_orig->Write();

	    dout->cd();

	    // Set uniform naming scheme and graphical style
	    g->SetName(Form("%s_%s_a%1.0f",tt,ss,100.*alphas[ialpha]));
	    g->UseCurrentStyle(); // Basic TDR style
	    //g->SetMarkerStyle(style[t]);
	    g->SetMarkerStyle(style[s][t]);
	    g->SetMarkerColor(color[s]);
	    g->SetLineColor(color[s]);
	    if (s=="pfjet" || (s=="zjet" && isfrac) ||
		(s=="zlljet" && isfrac)) {
	      g->SetMarkerStyle(d=="mc" ? style[s+"_"+d][t] : style[s][t]);
	      g->SetMarkerColor(color[s+"_"+t]);
	      g->SetLineColor(color[s+"_"+t]);
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
  fout->cd("ratio/eta00-13");
  hmzee->Write("mass_zeejet_a30");
  hmzmm->Write("mass_zmmjet_a30");
  //hmz->Write("mass_zjet_a30"); // v26
  hmzee1->Write("mass_zeejet_a100");
  hmzmm1->Write("mass_zmmjet_a100");
  //hmz1->Write("mass_zjet_a100");

  fout->cd("data/eta00-13");
  hmzee_dt->Write("mass_zeejet_a30");
  hmzmm_dt->Write("mass_zmmjet_a30");
  //hmz_dt->Write("mass_zjet_a30");
  hmzee1_dt->Write("mass_zeejet_a100");
  hmzmm1_dt->Write("mass_zmmjet_a100");
  //hmz1_dt->Write("mass_zjet_a100");

  fout->cd("mc/eta00-13");
  hmzee_mc->Write("mass_zeejet_a30");
  hmzmm_mc->Write("mass_zmmjet_a30");
  //hmz_mc->Write("mass_zjet_a30");
  hmzee1_mc->Write("mass_zeejet_a100");
  hmzmm1_mc->Write("mass_zmmjet_a100");
  //hmz1_mc->Write("mass_zjet_a100");

  fdj->Close();
  fp->Close(); // single file
  //{ fp1->Close(); fp2->Close(); } // two files
  fzee->Close();
  fzmm->Close();
  fz->Close();
  curdir->cd();


  /////////////////////////////////////////////////////
  // Calculate JEC central values and uncertainties, //
  // and store them for easy visualialization later  //
  /////////////////////////////////////////////////////

  map<string,const char*> mera;
  mera["2018ABCD"] = "ABCD";
  mera["2018A"] = "A";
  mera["2018B"] = "B";
  mera["2018C"] = "C";
  mera["2018D"] = "D";

  // Calculate L2L3Res with JEC uncertainty
  {

    const char *s, *s2;
    const char *cd = "CondFormats/JetMETObjects/data";
    const char *cpl = "../JERCProtoLab/Summer19UL18/global_fit/V3A2M1J2";
    const char *ce = epoch.c_str();
    
    // New JEC for plotting on the back
    FactorizedJetCorrector *jec(0), *jeca(0);
    FactorizedJetCorrector *jecb(0), *jecc(0), *jecd(0), *jece(0), *jecf(0);
    double jecwb(0), jecwc(0), jecwd(0), jecwe(0), jecwf(0);       // for BCDEF
    double jecwa(0); // for UL18
    {
      /*
      jec = getFJC("","",
		   Form("Summer19UL17_Run%s_V2M5_SimpleL1_DATA_L2L3Residual",
			epoch=="BCDEF" ? "B" : epoch.c_str()));
      */
      jec = (isUL18 ?
	     //getFJC("","",Form("Summer19UL18_Run%s_V3M1_DATA_L2L3Residual",
	     getFJC("","",Form("Summer19UL18_Run%s_V3A2M1J2_DATA_L2L3Residual",
			       mera[epoch]),cpl) :
	     getFJC("","",Form("Summer19UL17_Run%s_V4_DATA_L2L3Residual",
			       epoch=="BCDEF" ? "B" : epoch.c_str())));

      if (epoch=="2018ABCD") {
	
	jeca =getFJC("","","Summer19UL18_RunA_V3M1_DATA_L2L3Residual");
	jecb =getFJC("","","Summer19UL18_RunB_V3M1_DATA_L2L3Residual");
	jecc =getFJC("","","Summer19UL18_RunC_V3M1_DATA_L2L3Residual");
	jecd =getFJC("","","Summer19UL18_RunD_V3M1_DATA_L2L3Residual");
	jecf =0;
	
	//lumtot = 9.6+4.2+9.3; // UL17 just CDE, drop B,F
	double lumtot = 14.0+7.1+6.9+31.9; //59.9/fb
	jecwa = 14.0/lumtot;
	jecwb = 7.1/lumtot;
	jecwc = 6.9/lumtot;
	jecwd = 31.9/lumtot;
	jecwe = 0;
	jecwf = 0;
      }
      if (epoch=="BCDEF") {
	/*
	jecb =getFJC("","","Summer19UL17_RunB_V2M5_SimpleL1_DATA_L2L3Residual");
	jecc =getFJC("","","Summer19UL17_RunC_V2M5_SimpleL1_DATA_L2L3Residual");
	jecd =getFJC("","","Summer19UL17_RunD_V2M5_SimpleL1_DATA_L2L3Residual");
	jece =getFJC("","","Summer19UL17_RunE_V2M5_SimpleL1_DATA_L2L3Residual");
	jecf =getFJC("","","Summer19UL17_RunF_V2M5_SimpleL1_DATA_L2L3Residual");
	*/
	jeca =0;
	jecb =getFJC("","","Summer19UL17_RunB_V4_DATA_L2L3Residual");
	jecc =getFJC("","","Summer19UL17_RunC_V4_DATA_L2L3Residual");
	jecd =getFJC("","","Summer19UL17_RunD_V4_DATA_L2L3Residual");
	jece =getFJC("","","Summer19UL17_RunE_V4_DATA_L2L3Residual");
	jecf =getFJC("","","Summer19UL17_RunF_V4_DATA_L2L3Residual");
	
	// luminosity from Hugues Lattaud, photon+jet 11 June 2018
	double lumtot = 4.8+9.6+4.2+9.3+13.4; // 41.3/fb
	jecwa = 0;
	jecwb = 4.8/lumtot;
	jecwc = 9.6/lumtot;
	jecwd = 4.2/lumtot;
	jecwe = 9.3/lumtot;
	jecwf = 13.4/lumtot;
      }
    }
    //}

    FactorizedJetCorrector *jecrun1, *jecold, *jecl1flat,*jecl1rcdt, *jecl1pt,
      *jecl1mc, *jecl1rcmc, *jecl1dt, *jecl1c, *jecl1s, *jecl1rc,
      *jecl2c, *jecl2s;

    // Reference Run I JEC for plotting on the back
    jecrun1 =   getFJC("","","Winter14_V8_DATA_L2L3Residual_AK5PFchs"); 
    // Store old JEC for undoing it in global fit
    // NB: global fit ideally needs a temporary pT-flat L3Res as input
    // But even with this pT-dependent L2Res can cause problems
    if (isUL18) {
      jecold =    getFJC("","",Form("Summer19UL18_Run%s_V3_DATA_L2Residual",
				    epoch=="2018ABCD" ? "C" : mera[epoch]));
    }
    else {
      jecold =    getFJC("","","Summer19UL17_RunC_V1_SimpleL1_DATA_L2Residual");
    }
 
    // Difference between pT-dependent and flat L1
    jecl1flat = getFJC("Summer19UL17_RunC_V1_SimpleL1_DATA_L1RC");
    jecl1rcdt = getFJC("Summer19UL17_RunC_V1_SimpleL1_DATA_L1RC");
    jecl1pt =   getFJC("Summer19UL17_RunC_V1_SimpleL1_DATA_L1FastJet");
    jecl1mc =   getFJC("Summer19UL17_V1_SimpleL1_MC_L1FastJet");
    jecl1rcmc = getFJC("Summer19UL17_V1_SimpleL1_MC_L1RC");
    jecl1dt =   getFJC("Summer19UL17_RunC_V1_SimpleL1_DATA_L1FastJet");

    // With L2Relative included
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

    // New additions (RC: mean vs median; L1: Simple vs SemiSimple)
    FactorizedJetCorrector *jrcdnom, *jrcmnom, *jrcdalt, *jrcmalt; // 1
    FactorizedJetCorrector *jl1dnomsf, *jl1daltsf, *jl1mnomsf; // 2a
    FactorizedJetCorrector *jl1dnom, *jl1mnom, *jl1dalt, *jl1malt; // 2b,2c

    //jrcdnom = getFJC("Summer19UL17_RunE_V1_SimpleL1_DATA_L1RC"); // (1)
    //jrcmnom = getFJC("Summer19UL17_V1_SimpleL1_MC_L1RC"); // (1)
    // Use consistent BCDEF files instead of the official files above
    jrcdnom = getFJC("UL17_RunBCDEF_V1_DATA_L1RC"); // (1)
    jrcmnom = getFJC("UL17_RunBCDEF_V1_MC_L1RC"); // (1)
    jrcdalt = getFJC("UL17_MED_RunBCDEF_V1_DATA_L1RC"); // (1)
    jrcmalt = getFJC("UL17_MED_RunBCDEF_V1_MC_L1RC"); // (1)

    // Use place-holder for the missing custom data files
    //jl1dnomsf = getFJC("Summer19UL17_RunE_V1_SimpleL1_DATA_L1FastJet"); // (2a) -- placeholder for the two below
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
  //jl1dnom = getFJC("Summer19UL17_RunE_V1_SimpleL1_DATA_L1FastJet", // (2b)--PH
    jl1dnom = getFJC("UL17_RunBCDEF_L1Simple_L1FastJet", //(2b)
		     "Summer19UL17_V1_SimpleL1_MC_L2Relative"); // (2b)
    jl1mnom = getFJC("Summer19UL17_V1_SimpleL1_MC_L1FastJet", // (2b)
		     "Summer19UL17_V1_SimpleL1_MC_L2Relative"); // (2b)
    //jl1dalt = getFJC("Summer19UL17_RunBCDEF_V1_SemiSimpleL1_DATA_L1FastJet",
  //jl1dalt = getFJC("ParallelMCL1_L1FastJet_AK4PFchs_L1SemiSimple",//(2b)--PH
  //jl1dalt = getFJC("ParallelMCL1_RunE_L1SemiSimple_DATA_L1FastJet",//(2b)--PH
    jl1dalt = getFJC("UL17_RunBCDEF_L1SemiSimple_L1FastJet", //(2b)
		     "ParallelMCL1_L2Relative_AK4PFchs_L1SemiSimple"
		     "_L2L3Splines_ptclip8"); // (2b)
    jl1malt = getFJC("ParallelMCL1_L1FastJet_AK4PFchs_L1SemiSimple", // (2b)
		     "ParallelMCL1_L2Relative_AK4PFchs_L1SemiSimple"
		     "_L2L3Splines_ptclip8"); // (2b)
    

    // Run I uncertainty => 80XV6 uncertainty
    s = Form("%s/Winter14_V8_DATA_UncertaintySources_AK5PFchs.txt",cd); // V8 Run I
    //s2 = "TotalNoFlavorNoTime";
    s2 = "SubTotalAbsolute"; // 07AugV4
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_ref1 = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc_ref1 = new JetCorrectionUncertainty(*p_ref1);

    // Total uncertainty, excluding Flavor and Time
    s =Form("%s/Fall17_17Nov2017B_V32_DATA_UncertaintySources_AK4PFchs.txt",cd);
    s2 = "TotalNoFlavorNoTime";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_unc = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p_unc);

    // Partial uncertainties
    s =Form("%s/Fall17_17Nov2017B_V32_DATA_UncertaintySources_AK4PFchs.txt",cd);
    //s2 = "TotalNoFlavorNoTime";
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
      //assert(gDirectory->cd(dd1)); // broke in ROOT 6.04/08
      assert(dout1->cd(dd1));
      TDirectory *dout2 = dout1->GetDirectory(dd1); assert(dout2);
      dout2->cd();
      cout << "Eta bin:" << eta1 <<"-"<< eta2 << endl;

      
      //const double ptbins[] = {15, 16, 17, 18, 19, 20, 22, 24, 26, 28,
			       //30, 40, 50, 60, 75, 100, 125, 155, 180,
      //const double ptbins[] = {29,30, 40, 50, 60, 75, 100, 125, 155, 180,
      const double ptbins[] = {15, 16, 18, 20, 22, 25,
			       30, 35, 40, 50, 60, 75, 100, 125, 155, 180,
			       210, 250, 300, 350, 400, 500, 600, 800, 1000,
			       1200, 1500, //1600, 1601};
			       //1800, 2100, 2500, 3000, 3500, 3501};
			       1800, 2100, 2400, 2700, 3000, 3300, 3600,
			       3900, 4200, 4500, 4501};
      const int npt = sizeof(ptbins)/sizeof(ptbins[0])-1;

      // Uncertainty bands
      TH1D *herr = new TH1D("herr",";p_{T} (GeV);JEC uncertainty;",
			    npt, &ptbins[0]);
      TH1D *herr_ref = new TH1D("herr_ref",";p_{T} (GeV);TotalNoFlavorNoTime;",
				npt, &ptbins[0]);
      TH1D *herr_spr = new TH1D("herr_spr",";p_{T} (GeV);SPR uncertainty;",
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

      if(rp_debug) cout << "create reference JES bands" << endl;
      for (int ipt = 1; ipt != herr->GetNbinsX()+1; ++ipt) {

	double pt = herr->GetBinCenter(ipt);

	int jpt = h2eta->GetXaxis()->FindBin(pt);
	//jpt = max(0,min(jpt, h2eta->GetXaxis()->GetNbins()));
	// 20200622: avoid underflow bin for pT<30 GeV
	jpt = max(1,min(jpt, h2eta->GetXaxis()->GetNbins()));

	// Take Z+jet eta,pT distribution for correctly averaging JEC
	TH1D *heta = h2eta->ProjectionY(Form("heta_%d",ipt), jpt, jpt);
	const int ieta2 = heta->FindBin(eta2);
	const int ieta1 = heta->FindBin(eta1);
	const int intw = heta->Integral(ieta1,ieta2-1);
	const int neta = 2*(ieta2 - ieta1);

	// Loop over fine eta bins to average JEC and uncertainties
	double sumval(0), sumerr2(0), sumw(0);
	double sumrun1(0), sumjes(0);
	double sumvall1flat(0), sumvall1pt(0), sumvall1mc(0), sumvall1rho(0);
	double sumvall1dt(0), sumvall1c(0), sumvall1s(0), sumvall1rc(0);
	double sumvall1rcdt(0), sumvall1rcmc(0);
	double sumvall2c(0), sumvall2s(0);
	double sumerr2_pt(0), sumerr2_hcal(0), sumerr2_ecal(0), sumerr2_pu(0);
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
	  double val = (epoch=="L4" ? 1 : 1./getJEC(jec,eta,pt));

	  if (epoch=="2018ABCD") {
	    assert(jeca); assert(jecb);assert(jecc);assert(jecd);
	    assert(fabs(jecwa+jecwb+jecwc+jecwd-1)<1e-4);
	    
	    double vala = 1./getJEC(jeca,eta,pt);
	    double valb = 1./getJEC(jecb,eta,pt);
	    double valc = 1./getJEC(jecc,eta,pt);
	    double vald = 1./getJEC(jecd,eta,pt);

	    val = jecwa*vala +jecwb*valb +jecwc*valc +jecwd*vald;
	  }
	  if (epoch=="BCDEF") {
	    assert(jecb);assert(jecc);assert(jecd);assert(jecf);
	    assert(fabs(jecwb+jecwc+jecwd+jecwe+jecwf-1)<1e-4);
	    
	    double valb = 1./getJEC(jecb,eta,pt);
	    double valc = 1./getJEC(jecc,eta,pt);
	    double vald = 1./getJEC(jecd,eta,pt);
	    double vale = 1./getJEC(jece,eta,pt);
	    double valf = 1./getJEC(jecf,eta,pt);

	    val = jecwb*valb +jecwc*valc +jecwd*vald +jecwe*vale +jecwf*valf;
	  }

          //if(rp_debug) cout << "ABC/ABCD special treatment done" << endl;
          if(CorLevel=="L1L2L3Res")val = 1.0; //to get proper bands during closure test
	  if(CorLevel=="L1L2L3Res") val = 1;
	  sumval += w*val;
	  sumw += w; // sum weights only once
	  
	  // reference JEC
	  double jesrun1 = (epoch=="L4" ? 1 : 1./getJEC(jecrun1,eta,pt));
	  sumrun1 += w*jesrun1;

	  // old JEC
	  double jes = (epoch=="L4" ? 1 : 1./getJEC(jecold,eta,pt));
	  sumjes += w*jes;

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
	double val = sumval / sumw;
        if(CorLevel=="L1L2L3Res")val = 1.0; //to get proper bands during closure test -- should already be 1.0 ... TESTING

	// normalize uncertainties (quadratic instead of linear addition)
	double err = sqrt(sumerr2 / sumw);
	double err_ref = sqrt(sumerr2_ref / sumw);
	double err_ref1 = sqrt(sumerr2_ref1 / sumw);
	double err_pt = sqrt(sumerr2_pt / sumw);
	double err_hcal = sqrt(sumerr2_hcal / sumw);
	double err_ecal = sqrt(sumerr2_ecal / sumw);
	double err_pu = sqrt(sumerr2_pu / sumw);
	double err_noflv = sqrt(sumerr2_noflv / sumw);

	// center uncertainties around JEC central value
	herr->SetBinContent(ipt, val);
	herr_ref->SetBinContent(ipt, val);
	herr_spr->SetBinContent(ipt, val);
	herr_pu->SetBinContent(ipt, 1);
	herr_noflv->SetBinContent(ipt, val);
	herr_mpf->SetBinContent(ipt, val);
	herr_pt->SetBinContent(ipt, val);

	herr->SetBinError(ipt, val*err);
	herr_ref->SetBinError(ipt, val*err_ref);
	herr_spr->SetBinError(ipt,val*sqrt(err_hcal*err_hcal
					   + err_ecal*err_ecal));
	herr_pu->SetBinError(ipt, val*err_pu);
	herr_noflv->SetBinError(ipt, val*err_noflv);
	herr_mpf->SetBinError(ipt, val*err_pt);
	herr_pt->SetBinError(ipt, val*sqrt(err_pt*err_pt+err_pu*err_pu));

	double run1 = (sumrun1 / sumw);
        if(CorLevel=="L1L2L3Res") run1 = 1.0;  //to get proper bands during closure test
	hrun1->SetBinContent(ipt, run1);
	hrun1->SetBinError(ipt, run1*err_ref1);//run1*err);

	double jes = (sumjes / sumw);
        if(CorLevel=="L1L2L3Res") jes = 1.0;  //to get proper bands during closure test
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
	/*
	double l1dnomsf = (sumvall1dnomsf / sumw);
	double offl1 = pt * (l1dnomsf - 1);
	double offrcdnom = pt * (rcdnom - 1);
	double offrcdalt = pt * (rcdalt - 1);
	double offrcmnom = pt * (rcmnom - 1);
	double offrcmalt = pt * (rcmalt - 1);
	double sfnom = offrcdnom / offrcmnom;
	double sfalt = offrcdalt / offrcmalt;
	hl1sf->SetBinContent(ipt, (pt+offl1*sfalt-offl1*sfnom)/pt); //(2a)
	hl1sfnom->SetBinContent(ipt, sfnom); //xtra
	hl1sfalt->SetBinContent(ipt, sfalt); //xtra
	hl1offdt->SetBinContent(ipt, offl1); //xtra
	hl1dt->SetBinContent(ipt, pt/(pt+offl1)); //xtra
	*/
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
	//hl1dalt->SetBinContent(ipt, l1dalt); //xtra
	//hl1malt->SetBinContent(ipt, l1malt); //xtra
	//
	double l1dnomsf20 = (sumvall1dnomsf20 / sumw);
	double l1daltsf20 = (sumvall1daltsf20 / sumw);
	hl1sfr->SetBinContent(ipt, (l1daltsf20/l1dnomsf20)
			      / (l1daltsf/l1dnomsf) ); //(3a)
	hl1sf20->SetBinContent(ipt, l1daltsf20/l1dnomsf20); //xtra
	hl1dt20->SetBinContent(ipt, l1dnomsf20); //xtra
	/*
	double offl120 = pt * (l1dnomsf20 - 1);
	hl1sfr->SetBinContent(ipt, ((pt+offl120*sfalt-offl120*sfnom)/pt)
			      / ((pt+offl1*sfalt-offl1*sfnom)/pt)); //(3a)
	hl1sf20->SetBinContent(ipt, (pt+offl120*sfalt-offl120*sfnom)/pt); //xtra
	hl1dt20->SetBinContent(ipt, (pt+offl120)/(pt+offl1)); //xtra
	*/
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

      herr_ref->SetMarkerSize(0);
      herr_ref->SetFillStyle(1001);
      herr_ref->SetFillColorAlpha(kYellow+1,0.5);
      herr_ref->Write();

      herr_spr->SetMarkerSize(0);
      herr_spr->SetFillStyle(1001);
      herr_spr->SetFillColorAlpha(kYellow,0.5);
      herr_spr->Write();

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
    JetCorrectorParameters *pl2 = new JetCorrectorParameters(s);
    cout << s << endl << flush;
    v.push_back(*pl2);
  }
  if (res!="") {
    s = Form("%s/%s.txt",cd,cres);
    JetCorrectorParameters *pres = new JetCorrectorParameters(s);
    cout << s << endl << flush;
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
