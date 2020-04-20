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

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "tools.h"
#include "tdrstyle_mod14.C"

#include <string>
#include <iostream>
#include <vector>
#include <map>

using namespace std;

// rho used to calculate l1bias
//const double gRho = 15.36; // 2017-06-02 for 2016 data
const double gRho = 19.21; // EOY17 DE jt450
const bool _dcsonly = false;
const bool rp_debug = true; // verbose messages

// Appy mass corrections to Z+jet
// Update 20200325: start moving these to globalFitSyst.C as syst. eigenvectors
// Need to correct Zee and Zmm here to combine them, but photon could be later
bool useFixedFit = true; // with minitools/drawZmass.C
double fitUncMin = 0.00000; // Some bug if unc<0?
bool correctZmmMass = true; // pT with mumu mass
bool correctZeeMass = true; // pT with ee mass
bool correctGamMass = true; //!!UL17_V3 false, pT with ee mass at 2*pT
bool correctUncert = false;  // ll mass uncertainty => globalFitSyst.C
//
// Which binning to use for multijets ("leading" or "recoil")
string multijetMode = "leading"; // check also multijetModeS in globalFitSyst.C
bool correctMultijetLeading = true; // correct for difference to JER-hybrid

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
double fzeeptmin(30.);   // Zee+jet both methods
double fzmmptmin(30.);//30.);   // Zmm+jet both methods
// Additional cuts to Z+jet MPF / balance methods
//double fzmpfptmin(30.);   // Z+jet MPF
//double fzbalptmin(30.);//100);//30.);   // Z+jet pTbal

// multijet minimum and maximum pT
// (high for now until FSR working for low pT)
// 49: 228.6/78, 64: 155.7/76, 84: 135.8/74, 114: 124.7/72, 153: 120.7/70
// 20200330: 49:238.3/84, 64:164.6/82, 84:143.5/80, 114:131.4/78,
// ...20200330: 153:125.4/76 196:108.5/74 245:107.7/72
double fmultijetptmin(114);//114);//49);//64);//84);114);//153);
// 2640 good for BCDEF, but destabilizes some IOVs?
double fmultijetptmax(2116);//2640.);//1890.);//1800.);//3000.);


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
// Additional cuts to Z+jet MPF / balance methods
//double fzmpfptmax(700);//500.);  // Z+jet MPF
//double fzbalptmax(500);//300.);  // Z+jet pTbal

// Regular settings for fits with Z+jet
double fpmpfptmin(230.); double fpbalptmin(230.); double fzmpfptmin(30); double fzbalptmin(80.);/*30);*/ double fzmpfptmax(700); double fzbalptmax(500); // UL17
//double fpmpfptmin(100.); double fpbalptmin(100.); double fzmpfptmin(40); double fzbalptmin(100); double fzmpfptmax(500); double fzbalptmax(300); // EOY17 settings (incl. effective 40 GeV from bug and gamma+jet min pT)

// Special setting for "multijet only" fit 
//double fzmpfptmin(130); double fzbalptmin(130); double fzmpfptmax(175); double fzbalptmax(175); // option B
//double fzmpfptmin(175); double fzbalptmin(175); double fzmpfptmax(230); double fzbalptmax(230); // option A



//minimum event counts
const double neventsmin = 20.;

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

FactorizedJetCorrector *_thejec(0);
TF1 *fCorrPt(0);
Double_t funcCorrPt(Double_t *x, Double_t *p) {
  
  double eta = p[0];
  double pt = x[0];
  double rho = p[1];
  setEtaPtRho(_thejec, eta, pt, rho);

  return (_thejec->getCorrection() * pt);
}

//const double _rhoDE = 19.21; // EOY17 DE jt450
double getJEC(FactorizedJetCorrector *jec, double eta, double ptcorr,
	      double rho = gRho) {


  // terate to solve ptreco for given ptcorr
  _thejec = jec;
  if (!fCorrPt) fCorrPt = new TF1("fCorrPt",funcCorrPt,5,6500,2);
  fCorrPt->SetParameters(eta, rho);
  // Find ptreco that gives pTreco*JEC = pTcorr
  double ptreco = fCorrPt->GetX(ptcorr,5,6500);

  setEtaPtRho(jec, eta, ptreco, rho);

  return (jec->getCorrection());
} // getEtaPtE



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
  //TFile *fmj = new TFile(Form("rootfiles/multijet_20190911_JEC_Autunm18_V17_JER_Autumn18_V7/multijet_20190911_Run2018%s_P8CP5_jecV17_jerV7.root",fm_files[epoch]),"READ"); // LO Pythia8 off by 2.5% on multijet scale
  //TFile *fmj = new TFile(Form("rootfiles/multijet_20190911_JEC_Autunm18_V17_JER_Autumn18_V7/multijet_20190911_Run2018%s_MGP8CP5_jecV17_jerV4.root",fm_files[epoch]),"READ"); // MadGraph much better match to data than LO P8 (just not JER V4)
  //TFile *fmj = new TFile(Form("rootfiles/multijet_20190911_JEC_Autunm18_V17_JER_Autumn18_V7/multijet_20190912_Run2018%s_MC_jecV17_jerV7.root",fm_files[epoch]),"READ"); // All MC in one file (JERV4+ABC only for MG)
  //
  //TFile *fmj = new TFile(Form("rootfiles/multijet_20190911_JEC_Autunm18_V17_JER_Autumn18_V7/multijet_Rebin2_20190920_Run2018%s_jecV17_jerV7.root",fm_files[epoch]),"READ"); // All MC in one file (JERV4+ABC only for MG) EOY2017
  //if (!fmj || fmj->IsZombie()) 
  //fmj = new TFile(Form("rootfiles/multijet_20190911_JEC_Autunm18_V17_JER_Autumn18_V7/multijet_20190912_Run2018%s_MC_jecV17_jerV7.root",fm_files[epoch]),"READ"); // All MC in one file (JERV4+ABC only for MG) EOY2017
  //TFile *fmj = new TFile(Form("rootfiles/multijet_Rebin2_20200312_UL2017%s_jecV1_jerV3.root",fm_files[epoch]),"READ"); // UL2017
  //TFile *fmj = new TFile(Form("rootfiles/multijet_Rebin2_20200318_UL2017%s_jecV1_jerV3.root",fm_files[epoch]),"READ"); // UL2017
  //TFile *fmj = new TFile(Form("rootfiles/multijet_Rebin2_20200320_UL2017%s_jecV1_jerV3.root",fm_files[epoch]),"READ"); // UL2017
  //TFile *fmj = new TFile(Form("rootfiles/multijet_Rebin2_20200330_UL2017%s_jecV1_jerV3.root",fm_files[epoch]),"READ"); // UL2017 (TProfile means?)
  //TFile *fmj = new TFile(Form("rootfiles/multijet_Rebin2_20200404_UL2017%s_ComplexL1_jecV1_jerV1.root",fm_files[epoch]),"READ"); // UL2017 (Complex+new JER etc.)
  //TFile *fmj = new TFile(Form("rootfiles/multijet_Rebin2_20200404_UL2017%s_SimpleL1_jecV1_jerV3.root",fm_files[epoch]),"READ"); // UL2017 (Simple+old JER)
  //TFile *fmj = new TFile(Form("rootfiles/multijet_Rebin2_20200404_UL2017%s_ComplexL1_jecV1_jerV3.root",fm_files[epoch]),"READ"); // UL2017 (Complex+new JER, partially complete)
  TFile *fmj = new TFile(Form("rootfiles/multijet_Rebin2_20200419_UL2017%s_SimpleL1_jecV2_jerV1.root",fm_files[epoch]),"READ"); // UL2017_V2
  assert(fmj && !fmj->IsZombie());
  //TFile *fmj =0;

  // Hugues Lattaud, 2017 V27 inputs with L2Res (not incl. JER SF as in V28)
  // https://indico.cern.ch/event/765393/#47-l3res-gammajets-with-fall17
  //map<string,const char*> fp_files;
  //fp_files["B"] = "B_B";
  //fp_files["C"] = "C_C";
  //fp_files["D"] = "D_D";
  //fp_files["E"] = "E_E";
  //fp_files["F"] = "F_F";
  //fp_files["BCDEF"] = "BCDEF";
  //TFile *fp = new TFile(Form("rootfiles/Gjet_combinationfile_17_Nov_V27_L2res_%s.root", fp_files[epoch]),"READ");

  map<string,const char*> fp_files;
  fp_files["B"] = "B";
  fp_files["C"] = "C";//"E";
  fp_files["BC"] = "BC";
  fp_files["D"] = "D";
  fp_files["E"] = "E";
  fp_files["DE"] = "DE";
  fp_files["F"] = "F";
  fp_files["BCDEF"] = "BCDEF";//"E";
  TFile *fp = new TFile(Form("rootfiles/2020-04-02/SimpleL1_only_L2Res/Gjet_combinationfile_SimpleL1_only_L2Res_%s_SimpleL1_only_L2Res.root",fp_files[epoch]),"READ");
  //TFile *fp = new TFile(Form("rootfiles/2020-04-02/ComplexL1_only_L2Res/Gjet_combinationfile_ComplexL1_only_L2Res_%s_ComplexL1_only_L2Res.root",fp_files[epoch]),"READ");

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

  // Daniel Savoiu, 2017 V28 inputs with L2Res (incl. JER SF on top of V27)
  // https://indico.cern.ch/event/763555/#31-l3res-and-closure-test-resu
  map<string,const char*> fz_files;
  fz_files["B"] = "B";
  fz_files["C"] = "C";
  fz_files["D"] = "D";
  fz_files["E"] = "E";
  fz_files["F"] = "F";
  fz_files["BCDEF"] = "BCDEF"; // 20200317 update
  //fz_files["BCDEF"] = "B"; // UL2017-v1
  //TFile *fzmm = new TFile(Form("rootfiles/zjet_combination_Fall17_JECV28_Zmm_%s_2018-10-08.root",fz_files[epoch]),"READ"); // EOY2017
  //TFile *fzee = new TFile(Form("rootfiles/zjet_combination_Fall17_JECV28_Zee_%s_2018-10-08.root",fz_files[epoch]),"READ"); // EOY2017
  //TFile *fzmm = new TFile("rootfiles/ZJetCombination_Zmm_09Aug2019_Summer19UL17_JECV1_SimpleL1.root","READ"); // UL2017-v1
  //TFile *fzee = new TFile("rootfiles/ZJetCombination_Zee_09Aug2019_Summer19UL17_JECV1_SimpleL1.root","READ"); // UL2017-v1
  //TFile *fzmm = new TFile("rootfiles/ZJetCombination_Zmm_09Aug2019_Summer19UL17_JECV1_ComplexL1.root","READ"); // UL2017-v1b
  //TFile *fzee = new TFile("rootfiles/ZJetCombination_Zee_09Aug2019_Summer19UL17_JECV1_ComplexL1.root","READ"); // UL2017-v1b
  //TFile *fzmm = new TFile("rootfiles/ZJetCombination_Zmm_DYJets_Madgraph_09Aug2019_Summer19UL17_JECV1_ComplexL1.root","READ"); // UL2017-v1b
  //TFile *fzee = new TFile("rootfiles/ZJetCombination_Zee_DYJets_Madgraph_09Aug2019_Summer19UL17_JECV1_ComplexL1.root","READ"); // UL2017-v1b
  TFile *fzmm = new TFile("rootfiles/ZJetCombination_Zmm_DYJets_Madgraph_09Aug2019_Summer19UL17_JECV1_SimpleL1.root","READ"); // UL2017-v1b
  TFile *fzee = new TFile("rootfiles/ZJetCombination_Zee_DYJets_Madgraph_09Aug2019_Summer19UL17_JECV1_SimpleL1.root","READ"); // UL2017-v1b
  assert(fzmm && !fzmm->IsZombie());
  assert(fzee && !fzee->IsZombie());

  // UL2017-v1
  fzmm->cd(Form("Run2017%s",fz_files[epoch])); fzmm = (TFile*)gDirectory;
  fzee->cd(Form("Run2017%s",fz_files[epoch])); fzee = (TFile*)gDirectory;


  if(CorLevel=="L1L2L3Res"){
    // Monday 22 Oct V31 closure
    // https://indico.cern.ch/event/767001/#52-closure-test-for-fall17_17n
    fzmm = new TFile(Form("rootfiles/zjet_combination_Fall17_JECV31_Zmm_%s_2018-10-26.root",fz_files[epoch]),"READ"); // https://indico.cern.ch/event/759977/#35-closure-test-for-fall17_17n
    fzee = new TFile(Form("rootfiles/zjet_combination_Fall17_JECV31_Zee_%s_2018-10-26.root",fz_files[epoch]),"READ"); // for October 2018 release
 }

  // Link to Z mass files (same as above now) and histograms
  // This is used for scaling Z mass back to 1 for Zee and Zmm
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
      // 17Nov17_V10 BCDEF (EOY2017)
      //f1mzee->SetParameters(1.00246, 0.00214, 0.00116);
      //f1ezee->SetParameters(+3.73e-08, +1.17e-07, +2.02e-07,
      //                     +1.5e-08, -5.07e-08, -6.39e-08);
      // UL17 RunB fit with minitools/drawZmass.C (3 free)
      //f1mzee->SetParameters(0.99904, 0.00170, 0.00156);
      //f1ezee->SetParameters(+1.94e-07, +5.65e-07, +9.65e-07,
      //                    +8.95e-08, -2.53e-07,  -2.2e-07);
      // UL17 RunB fit with minitools/drawZmass.C (p1,p2 fixed to EOY2017)
      //f1mzee->SetParameters(0.99917, 0.00214, 0.00116);
      //f1ezee->SetParameters(+1.26e-07,        +0,        +0,
      //                           +0,        +0,        +0);
      // UL17 RunC fit with minitools/drawZmass.C (p1,p2 fixed to EOY2017)
      //f1mzee->SetParameters(0.99772, 0.00214, 0.00116);
      //f1mzee->SetParameters(0.997557, 0.00214, 0.00116); // p0 from B+C+D+E+F
      //f1ezee->SetParameters( +6.8e-08,        +0,        +0,
      //                           +0,        +0,        +0);
      // UL17 RunBCDEF fit with minitools/drawZmass.C
      f1mzee->SetParameters(0.99780, 0.00225, 0.00031);
      f1ezee->SetParameters( +6.1e-08, +1.66e-07, +3.27e-07,
                            +2.56e-08,  -8.5e-08, -5.39e-08);

      // 17Nov17_V10 BCDEF (EOY2017)
      //f1mgam->SetParameters(1.00246, 0.00214, 0.00116);
      //f1egam->SetParameters(+3.73e-08, +1.17e-07, +2.02e-07,
      //		    +1.5e-08, -5.07e-08, -6.39e-08);
      // UL17 RunBCDEF Zee from above
      f1mgam->SetParameters(0.99780, 0.00225, 0.00031);
      f1egam->SetParameters( +6.1e-08, +1.66e-07, +3.27e-07,
                            +2.56e-08,  -8.5e-08, -5.39e-08);
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
      // 17Nov17_V10 BCDEF (EOY2017)
      //f1mzmm->SetParameters(0.99854, 0.00000, 0.00000);
      //f1ezmm->SetParameters(+1.53e-08,        +0,        +0,
      //                           +0,        +0,        +0);
      // UL17 C fit with minitools/drawZmass.C
      //f1mzmm->SetParameters(0.99829, 0.00000, 0.00000);
      //f1mzmm->SetParameters(0.998235, 0.00000, 0.00000); // p0 from B+C+D+E+F
      // UL17 BCDEF fit with minitools/drawZmass.C
      f1mzmm->SetParameters(0.99821, 0.00000, 0.00000);
      f1ezmm->SetParameters(+3.43e-08,        +0,        +0,
                                   +0,        +0,        +0);
    }
    else
      hmzmm->Fit(f1mzmm);
  }


  // \END copy-paste from minitools/drawZmass.C


  TF1 *f2 = new TF1("f2","[p0]*pow(x,2)+[p1]",30,3000);
  f2->SetParameters(0,0);

  if (correctMultijetLeading && multijetMode=="leading") {
    f2->SetParameters(1.105e-07, 0.03063);
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
  files["dijet"] = fdj;
  files["multijet"] = fmj;
  files["gamjet"] = fp;
  files["zeejet"] = fzee;
  files["zmmjet"] = fzmm;


  // Map variable names used in different files to common scheme
  map<string, map<string, const char*> > rename;

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
  rename["multijet"]["ratiocrecoil"] = "CRecoil";
  rename["multijet"]["datacrecoil"] = "CRecoil";
  rename["multijet"]["mccrecoil"] = "CRecoil_MG";//P8";
  // Switch to leading jet pT binning instead or recoil pT binning
  if (multijetMode=="leading") {
    rename["multijet"]["ratiompfchs1"] = "MPF_leading_L1L2Res";
    rename["multijet"]["ratioptchs"] = "MJB_leading_L1L2Res";
    rename["multijet"]["datampfchs1"] = "MPF_leading_L1L2Res";
    rename["multijet"]["dataptchs"] = "MJB_leading_L1L2Res";
    rename["multijet"]["mcmpfchs1"] = "MPF_leading_MG";
    rename["multijet"]["mcptchs"] = "MJB_leading_MG";
  }
  else if (multijetMode=="recoil") {
    rename["multijet"]["ratiompfchs1"] = "MPF_recoil_L1L2Res";
    rename["multijet"]["ratioptchs"] = "MJB_recoil_L1L2Res";
    rename["multijet"]["datampfchs1"] = "MPF_recoil_L1L2Res";
    rename["multijet"]["dataptchs"] = "MJB_recoil_L1L2Res";
    rename["multijet"]["mcmpfchs1"] = "MPF_recoil_MG";
    rename["multijet"]["mcptchs"] = "MJB_recoil_MG";
  }
  else
    assert(false);


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

  rename["zmmjet"]["ratio"] = "Ratio";
  rename["zmmjet"]["data"] = "Data";
  rename["zmmjet"]["mc"] = "MC";
  rename["zmmjet"]["mpfchs"] = "MPF-notypeI_CHS";
  rename["zmmjet"]["mpfchs1"] = "MPF_CHS"; 
  rename["zmmjet"]["ptchs"] = "PtBal_CHS"; 
  rename["zmmjet"]["counts"] = "RawNEvents_CHS";
  

  // color and style codes
  map<string, map<string, int> > style;

  style["dijet"]["mpfchs1"] = kFullDiamond;//kFullCircle;
  style["dijet"]["ptchs"] = kOpenDiamond;//kOpenSquare;
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

  map<string, int> color;
  color["dijet"] = kBlack;
  color["multijet"] = kBlack;
  color["gamjet"] = kBlue;
  color["zeejet"] = kGreen+2;
  color["zmmjet"] = kRed;
  color["zlljet"] = kMagenta+2;

  map<double, double> size;
  size[0.10] = 0.6;
  size[0.15] = 0.8;
  size[0.20] = 1.0;
  size[0.30] = 1.2;


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

  vector<string> sets;
  sets.push_back("dijet");
  sets.push_back("multijet");
  sets.push_back("gamjet");
  sets.push_back("zeejet");
  sets.push_back("zmmjet");
  sets.push_back("zlljet");

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

	  // Take pT and MPF from different files for gamma+jet (or not)
	  //if (s=="gamjet" && f==0) {
	  //if (t=="mpfchs" || t=="mpfchs1") f = fp;//fp1;
	  //if (t=="ptchs") f = fp;//fp2;
	  //}
	  assert(f || s=="zlljet");

	  // C_recoil is only meant for multijets
	  if (t=="crecoil" && s!="multijet") continue;

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
	    assert(f || s=="zlljet");


            if (t=="counts" && s!="zmmjet" && s!="zeejet" && s!="gamjet")
              continue; // counts available only for z+jet and gamjet, so far

	    if (s=="multijet" && (!((epoch!="L4" &&
				     fabs(eta1-0)<0.1 && fabs(eta2-1.3)<0.1) ||
				    (epoch=="L4" &&
				     fabs(eta1-0)<0.1 && fabs(eta2-2.4)<0.1))
				  //|| fabs(alpha-0.10)<0.01 
				  // || fabs(alpha-0.15)<0.01
				  //|| fabs(alpha-0.20)<0.01))
				  || fabs(alpha-0.10)<0.01))
	      continue; // only barrel for multijet balance, pT=(15),(20),30
	    if (s=="gamjet"  && fabs(eta1-3.2)<0.1) { eta1=3.0; eta2=3.2; }

	    // Reconstruct naming scheme used in each of the files
	    // If non-conventional naming schemes, patch here
	    const char *c(0);
	    if (s=="dijet") {
	      c = Form("%s/eta_%02.0f_%02.0f/%s_%s_a%1.0f",
                       dd, 10.001*eta1, 10.001*eta2, rename[s][t], ss, 100.*alpha); 
	    } // dijet
	    if (s=="multijet") {
	      //c = Form("%s/Pt%1.0f/%s", rename[s][d], 100.*alpha, rename[s][t]);
	      c = Form("%s/Pt%1.0f/%s", rename[s][d], 100.*alpha, rename[s][(d+t)]);
	    } // multijet
	    
	    //if (s=="gamjet" && t=="counts" && d=="ratio")continue;
	    if (s=="gamjet" && t=="counts") {
	      c = Form("%s%s_a%1.0f_eta%02.0f_%02.0f_%s",
		       rename[s]["mpfchs1"], d=="ratio" ? rename[s]["data"] :rename[s][d],
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
	    } // Z+jet
	    assert(c || s=="zlljet");

	    TObject *obj = (s=="zlljet" ? 0 : f->Get(c));
	    if (!obj && s!="zlljet") {
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
	    
	    assert(obj);


            if (t=="counts" && (s=="zmmjet" || s=="zeejet" || s=="gamjet") ){ // write out counts to jecdata.root (as TH1F)
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

	    // If data stored in TH1D instead of TGraphErrors, patch here
	    // Patch for dijet file that has TH1D's instead of graphs
	    if (obj->InheritsFrom("TH1")) {
	      obj = new TGraphErrors((TH1D*)obj);
	    }

	    assert(obj->InheritsFrom("TGraphErrors"));
	    TGraphErrors *g = (TGraphErrors*)obj;

	    // Clean out empty points from TH1D->TGraphErrors conversion
	    for (int i = g->GetN()-1; i != -1; --i) {
              assert(i<=g->GetN()-1);
	      // Clean out spurious empty pooints
	      if (g->GetY()[i]==0 && g->GetEY()[i]==0) g->RemovePoint(i);
              // remove points with large error (more than 0.2 right now)
              else if (g->GetEY()[i]>0.2)  g->RemovePoint(i);
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
	      else if (s=="multijet" &&
		       (g->GetX()[i]<fmultijetptmin || g->GetX()[i]>fmultijetptmax))
		g->RemovePoint(i);
	      else if (s=="zeejet" && 
		       (g->GetX()[i]<fzeeptmin || g->GetX()[i]>fzeeptmax))
		g->RemovePoint(i);
	      else if (s=="zmmjet" && 
		       (g->GetX()[i]<fzmmptmin || g->GetX()[i]>fzmmptmax))
		g->RemovePoint(i);
	      else if ((s=="zmmjet" || s=="zeejet") && t=="mpfchs1" &&
		       (g->GetX()[i]<fzmpfptmin || g->GetX()[i]>fzmpfptmax))
		g->RemovePoint(i);
	      else if ((s=="zmmjet" || s=="zeejet") && t=="ptchs" &&
		       (g->GetX()[i]<fzbalptmin || g->GetX()[i]>fzbalptmax))
		g->RemovePoint(i);
              else if (s=="zmmjet" || s=="zeejet" || s=="gamjet"){ // patch: clean away points with low statistics based on event counts histograms, currently Z+jet
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

	    dout->cd();

	    // Set uniform naming scheme and graphical style
	    g->SetName(Form("%s_%s_a%1.0f",tt,ss,100.*alphas[ialpha]));
	    g->UseCurrentStyle(); // Basic TDR style
	    //g->SetMarkerStyle(style[t]);
	    g->SetMarkerStyle(style[s][t]);
	    g->SetMarkerColor(color[s]);
	    g->SetLineColor(color[s]);
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
  fout->cd("data/eta00-13");
  hmzee_dt->Write("mass_zeejet_a30");
  hmzmm_dt->Write("mass_zmmjet_a30");
  fout->cd("mc/eta00-13");
  hmzee_mc->Write("mass_zeejet_a30");
  hmzmm_mc->Write("mass_zmmjet_a30");

  fdj->Close();
  fp->Close(); // single file
  //{ fp1->Close(); fp2->Close(); } // two files
  fzee->Close();
  fzmm->Close();
  curdir->cd();


  /////////////////////////////////////////////////////
  // Calculate JEC central values and uncertainties, //
  // and store them for easy visualialization later  //
  /////////////////////////////////////////////////////

  // Calculate L2L3Res with JEC uncertainty
  {

    const char *s, *s2;
    const char *cd = "CondFormats/JetMETObjects/data";
    const char *ce = epoch.c_str();
    
    // New JEC for plotting on the back
    FactorizedJetCorrector *jec;
    FactorizedJetCorrector *jecc(0), *jecd(0), *jece(0), *jecf(0); // for BCDEF
    double jecwb(1), jecwc(0), jecwd(0), jecwe(0), jecwf(0);       // for BCDEF
    /*
    FactorizedJetCorrector *jecc(0), *jecde(0), *jecf(0); // for BCDEF
    double jecwb(1), jecwc(0), jecwde(0), jecwf(0);       // for BCDEF
    */
    {
      //s = Form("%s/Fall17_17Nov2017%s_V32_DATA_L2L3Residual_AK4PFchs.txt",cd,
      //       epoch=="BCDEF" ? "B" : 
      //       (epoch=="D"||epoch=="E") ? "DE" : epoch.c_str());
      //s = Form("%s/Summer19UL17_Run%s_V1_SimpleL1_DATA_L2L3Residual_AK4PFchs.txt",cd,epoch=="BCDEF"||epoch=="DE" ? "E" : epoch.c_str());
      s = Form("textFiles/UL17V2M3-L2L3Res+JERSF/Summer19UL17_Run%s_V2M3_SimpleL1_DATA_L2L3Residual_AK4PFchs.txt",epoch=="BCDEF" ? "B" : epoch.c_str());

      cout << s << endl;
      JetCorrectorParameters *par_l2l3res = new JetCorrectorParameters(s);
      vector<JetCorrectorParameters> vpar;
      vpar.push_back(*par_l2l3res);
      jec = new FactorizedJetCorrector(vpar);

      if (epoch=="BCDEF") {

	// luminosity from Hugues Lattaud, photon+jet 11 June 2018
	double lumtot = 4.8+9.6+4.2+9.3+13.4; // 41.3/fb
	jecwb = 4.8/lumtot;

	//s=Form("%s/Fall17_17Nov2017C_V32_DATA_L2L3Residual_AK4PFchs.txt",cd);
	//s = Form("%s/Summer19UL17_RunC_V1_SimpleL1_DATA_L2L3Residual_AK4PFchs.txt",cd);
	s = "textFiles/UL17V2M3-L2L3Res+JERSF/Summer19UL17_RunC_V2M3_SimpleL1_DATA_L2L3Residual_AK4PFchs.txt";
	cout << s << endl;
	JetCorrectorParameters *par_c = new JetCorrectorParameters(s);
	vector<JetCorrectorParameters> vpar_c;
	vpar_c.push_back(*par_c);
	jecc = new FactorizedJetCorrector(vpar_c);
	jecwc = 9.6/lumtot;


	//s=Form("%s/Fall17_17Nov2017D_V32_DATA_L2L3Residual_AK4PFchs.txt",cd);
	//s = Form("%s/Summer19UL17_RunD_V1_SimpleL1_DATA_L2L3Residual_AK4PFchs.txt",cd);	
	s = "textFiles/UL17V2M3-L2L3Res+JERSF/Summer19UL17_RunD_V2M3_SimpleL1_DATA_L2L3Residual_AK4PFchs.txt";
	cout << s << endl;
	JetCorrectorParameters *par_d = new JetCorrectorParameters(s);
	vector<JetCorrectorParameters> vpar_d;
	vpar_d.push_back(*par_d);
	jecd = new FactorizedJetCorrector(vpar_d);
	jecwd = 4.2/lumtot;

	//s=Form("%s/Fall17_17Nov2017E_V32_DATA_L2L3Residual_AK4PFchs.txt",cd);
	//s = Form("%s/Summer19UL17_RunE_V1_SimpleL1_DATA_L2L3Residual_AK4PFchs.txt",cd);	
	s = "textFiles/UL17V2M3-L2L3Res+JERSF/Summer19UL17_RunE_V2M3_SimpleL1_DATA_L2L3Residual_AK4PFchs.txt";
	cout << s << endl;
	JetCorrectorParameters *par_e = new JetCorrectorParameters(s);
	vector<JetCorrectorParameters> vpar_e;
	vpar_e.push_back(*par_e);
	jece = new FactorizedJetCorrector(vpar_e);
	jecwe = 9.3/lumtot;
	
	/*
	s=Form("%s/Fall17_17Nov2017DE_V32_DATA_L2L3Residual_AK4PFchs.txt",cd);
	cout << s << endl;
	JetCorrectorParameters *par_de = new JetCorrectorParameters(s);
	vector<JetCorrectorParameters> vpar_de;
	vpar_de.push_back(*par_de);
	jecde = new FactorizedJetCorrector(vpar_de);
	jecwde = (4.2+9.3)/lumtot;
	*/

	//s=Form("%s/Fall17_17Nov2017F_V10_DATA_L2L3Residual_AK4PFchs.txt",cd);
	//s=Form("%s/Fall17_17Nov2017F_V32_DATA_L2L3Residual_AK4PFchs.txt",cd);
	//s = Form("%s/Summer19UL17_RunF_V1_SimpleL1_DATA_L2L3Residual_AK4PFchs.txt",cd);	
	s = "textFiles/UL17V2M3-L2L3Res+JERSF/Summer19UL17_RunF_V2M3_SimpleL1_DATA_L2L3Residual_AK4PFchs.txt";
	cout << s << endl;
	JetCorrectorParameters *par_f = new JetCorrectorParameters(s);
	vector<JetCorrectorParameters> vpar_f;
	vpar_f.push_back(*par_f);
	jecf = new FactorizedJetCorrector(vpar_f);
	jecwf = 13.4/lumtot;
      }
    }

    // Reference Run I JEC for plotting on the back
    FactorizedJetCorrector *jecrun1;
    {
      s = Form("%s/Winter14_V8_DATA_L2L3Residual_AK5PFchs.txt",cd);
      cout << s << endl;
      JetCorrectorParameters *par_l2l3res = new JetCorrectorParameters(s);
      vector<JetCorrectorParameters> vpar;
      vpar.push_back(*par_l2l3res);
      jecrun1 = new FactorizedJetCorrector(vpar);
    }

    // Store old JEC for undoing it in global fit
    // NB: global fit ideally needs a temporary pT-flat L3Res as input
    // But even with this pT-dependent L2Res can cause problems
    FactorizedJetCorrector *jecold;
    {
      //s = Form("%s/Fall17_17Nov2017%s_V32_DATA_L2L3Residual_AK4PFchs.txt",cd,
      //epoch=="BCDEF"||epoch=="L4" ? "B" :
      //(epoch=="D"||epoch=="E") ? "DE" : epoch.c_str()); // Fall17
      //s = Form("%s/Summer19UL17_Run%s_V1_ComplexL1_DATA_L2Residual_AK4PFchs.txt", cd, "C");
      s = Form("%s/Summer19UL17_Run%s_V1_SimpleL1_DATA_L2Residual_AK4PFchs.txt",cd, "C");
      cout << s << endl;
      JetCorrectorParameters *par_old = new JetCorrectorParameters(s);
      vector<JetCorrectorParameters> v;
      v.push_back(*par_old);
      jecold = new FactorizedJetCorrector(v);
    }

    // Difference between pT-dependent and flat L1
    FactorizedJetCorrector *jecl1flat;
    {
      //s = Form("%s/Fall17_17Nov2017%s_V32_DATA_L1RC_AK4PFchs.txt",cd,epoch=="BCDEF"||epoch=="L4" ? "B" : (epoch=="D"||epoch=="E") ? "DE" : epoch.c_str()); // Fall17
      //s = Form("%s/Summer19UL17_Run%s_V1_ComplexL1_DATA_L1RC_AK4PFchs.txt",
      s = Form("%s/Summer19UL17_Run%s_V1_SimpleL1_DATA_L1RC_AK4PFchs.txt",
	       cd,"C");
      cout << s << endl << flush;
      JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
      vector<JetCorrectorParameters> v;
      v.push_back(*l1);
      jecl1flat = new FactorizedJetCorrector(v);
    }
    FactorizedJetCorrector *jecl1pt;
    {
      //s = Form("%s/Fall17_17Nov2017%s_V32_DATA_L1FastJet_AK4PFchs.txt",cd,epoch=="BCDEF"||epoch=="L4" ? "B" : (epoch=="D"||epoch=="E") ? "DE" : epoch.c_str()); // Fall17
      //s = Form("%s/Summer19UL17_Run%s_V1_ComplexL1_DATA_L1FastJet_AK4PFchs.txt",cd,"C");
      s = Form("%s/Summer19UL17_Run%s_V1_SimpleL1_DATA_L1FastJet_AK4PFchs.txt",
	       cd,"C");
      cout << s << endl << flush;
      JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
      vector<JetCorrectorParameters> v;
      v.push_back(*l1);
      jecl1pt = new FactorizedJetCorrector(v);
    }
    FactorizedJetCorrector *jecl1mc;
    {
      //s = Form("%s/Fall17_17Nov2017_V32_MC_L1FastJet_AK4PFchs.txt",cd); // Fall17
      //s = Form("%s/Summer19UL17_V1_ComplexL1_MC_L1FastJet_AK4PFchs.txt",cd);
      s = Form("%s/Summer19UL17_V1_SimpleL1_MC_L1FastJet_AK4PFchs.txt",cd);
      cout << s << endl << flush;
      JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
      vector<JetCorrectorParameters> v;
      v.push_back(*l1);
      jecl1mc = new FactorizedJetCorrector(v);
    }

    FactorizedJetCorrector *jecl1dt;
    {
      s = Form("%s/Summer19UL17_Run%s_V1_SimpleL1_DATA_L1FastJet_AK4PFchs.txt",
	       cd,"C");
      cout << s << endl << flush;
      JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
      vector<JetCorrectorParameters> v;
      v.push_back(*l1);
      jecl1dt = new FactorizedJetCorrector(v);
    }
    FactorizedJetCorrector *jecl1c;
    {
      s = Form("%s/Summer19UL17_V1_ComplexL1_MC_L1FastJet_AK4PFchs.txt",cd);
      cout << s << endl << flush;
      JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
      //s = Form("%s/Summer19UL17_V1_ComplexL1_MC_L2Relative_AK4PFchs.txt",cd);
      // This is on purpose SimpleL1_L2 for ComplexL1_L1
      // Need L1C L2 on same footing with L1S and L1RC
      s = Form("%s/Summer19UL17_V1_SimpleL1_MC_L2Relative_AK4PFchs.txt",cd);
      cout << s << endl << flush;
      JetCorrectorParameters *l2 = new JetCorrectorParameters(s);
      vector<JetCorrectorParameters> v;
      v.push_back(*l1);
      v.push_back(*l2);
      jecl1c = new FactorizedJetCorrector(v);
    }
    FactorizedJetCorrector *jecl1s;
    {
      s = Form("%s/Summer19UL17_V1_SimpleL1_MC_L1FastJet_AK4PFchs.txt",cd);
      cout << s << endl << flush;
      JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
      s = Form("%s/Summer19UL17_V1_SimpleL1_MC_L2Relative_AK4PFchs.txt",cd);
      cout << s << endl << flush;
      JetCorrectorParameters *l2 = new JetCorrectorParameters(s);
      vector<JetCorrectorParameters> v;
      v.push_back(*l1);
      v.push_back(*l2);
      jecl1s = new FactorizedJetCorrector(v);
    }
    FactorizedJetCorrector *jecl1rc;
    {
      s = Form("%s/Summer19UL17_V1_SimpleL1_MC_L1RC_AK4PFchs.txt",cd);
      cout << s << endl << flush;
      JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
      s = Form("%s/Summer19UL17_V1_SimpleL1_MC_L2Relative_AK4PFchs.txt",cd);
      cout << s << endl << flush;
      JetCorrectorParameters *l2 = new JetCorrectorParameters(s);
      vector<JetCorrectorParameters> v;
      v.push_back(*l1);
      v.push_back(*l2);
      jecl1rc = new FactorizedJetCorrector(v);
    }
    FactorizedJetCorrector *jecl2c;
    {
      s = Form("%s/Summer19UL17_V1_ComplexL1_MC_L2Relative_AK4PFchs.txt",cd);
      cout << s << endl << flush;
      JetCorrectorParameters *l2 = new JetCorrectorParameters(s);
      vector<JetCorrectorParameters> v;
      v.push_back(*l2);
      jecl2c = new FactorizedJetCorrector(v);
    }
    FactorizedJetCorrector *jecl2s;
    {
      s = Form("%s/Summer19UL17_V1_SimpleL1_MC_L2Relative_AK4PFchs.txt",cd);
      cout << s << endl << flush;
      JetCorrectorParameters *l2 = new JetCorrectorParameters(s);
      vector<JetCorrectorParameters> v;
      v.push_back(*l2);
      jecl2s = new FactorizedJetCorrector(v);
    }


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

      
      const double ptbins[] = {15, 16, 17, 18, 19, 20, 22, 24, 26, 28,
			       30, 40, 50, 60, 75, 100, 125, 155, 180,
      //const double ptbins[] = {29,30, 40, 50, 60, 75, 100, 125, 155, 180,
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
      TH1D *hl1cos = new TH1D("hl1cos",";p_{T} (GeV);"
			      "JEC L1Complex MC / JEC L1Simple MC;",
			      npt, &ptbins[0]);
      TH1D *hl1rcos = new TH1D("hl1rcos",";p_{T} (GeV);"
			       "JEC L1RC MC / JEC L1Simple MC;",
			       npt, &ptbins[0]);
      TH1D *hl2cos = new TH1D("hl2cos",";p_{T} (GeV);"
			      "JEC L2 of L1Complex MC / L1Simple MC;",
			      npt, &ptbins[0]);

      if(rp_debug) cout << "create reference JES bands" << endl;
      for (int ipt = 1; ipt != herr->GetNbinsX()+1; ++ipt) {

	double pt = herr->GetBinCenter(ipt);

	int jpt = h2eta->GetXaxis()->FindBin(pt);
	jpt = max(0,min(jpt, h2eta->GetXaxis()->GetNbins()));

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
	double sumvall2c(0), sumvall2s(0);
	double sumerr2_pt(0), sumerr2_hcal(0), sumerr2_ecal(0), sumerr2_pu(0);
	double sumerr2_noflv(0), sumerr2_ref(0), sumerr2_ref1(0);
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
	  //jec->setJetEta(eta);
	  //jec->setJetPt(pt);
	  //double val = (epoch=="L4" ? 1 : 1./jec->getCorrection());
	  double val = (epoch=="L4" ? 1 : 1./getJEC(jec,eta,pt));

	  if (epoch=="BCDEF") {

	    assert(jecc); assert(jecd); assert(jece); assert(jecf);
	    assert(fabs(jecwb+jecwc+jecwd+jecwe+jecwf-1)<1e-4);
	    /*
	    assert(jecc); assert(jecde); assert(jecf);
	    assert(fabs(jecwb+jecwc+jecwde+jecwf-1)<1e-4);
	    */
	    double valb = val;

	    //jecc->setJetEta(eta);
	    //jecc->setJetPt(pt);
	    //double valc = 1./jecc->getCorrection();
	    double valc = 1./getJEC(jecc,eta,pt);

	    //jecd->setJetEta(eta);
	    //jecd->setJetPt(pt);
	    //double vald = 1./jecd->getCorrection();
	    double vald = 1./getJEC(jecd,eta,pt);

	    //jece->setJetEta(eta);
	    //jece->setJetPt(pt);
	    //double vale = 1./jece->getCorrection();
	    double vale = 1./getJEC(jece,eta,pt);
	    /*
	    jecde->setJetEta(eta);
	    jecde->setJetPt(pt);
	    double valde = 1./jecde->getCorrection();
	    */
	    //jecf->setJetEta(eta);
	    //jecf->setJetPt(pt);
	    //double valf = 1./jecf->getCorrection();
	    double valf = 1./getJEC(jecf,eta,pt);

	    val = jecwb*valb +jecwc*valc +jecwd*vald +jecwe*vale +jecwf*valf;
	    /*
	    val = jecwb*valb +jecwc*valc +jecwde*valde +jecwf*valf;
	    */
	  }
          //if(rp_debug) cout << "ABC/ABCD special treatment done" << endl;
          if(CorLevel=="L1L2L3Res")val = 1.0; //to get proper bands during closure test
	  if(CorLevel=="L1L2L3Res") val = 1;
	  sumval += w*val;
	  sumw += w; // sum weights only once
	  
	  // reference JEC
	  //jecrun1->setJetEta(eta);
	  //jecrun1->setJetPt(pt);
	  //double jesrun1 = (epoch=="L4" ? 1 : 1./jecrun1->getCorrection());
	  double jesrun1 = (epoch=="L4" ? 1 : 1./getJEC(jecrun1,eta,pt));
	  sumrun1 += w*jesrun1;

	  // old JEC
	  //jecold->setJetEta(eta);
	  //jecold->setJetPt(pt);
	  //double jes = (epoch=="L4" ? 1 : 1./jecold->getCorrection());
	  double jes = (epoch=="L4" ? 1 : 1./getJEC(jecold,eta,pt));
	  sumjes += w*jes;

	  //jecl1flat->setJetEta(eta);
	  //jecl1flat->setJetPt(pt);
	  //jecl1flat->setRho(gRho);
	  //jecl1flat->setJetA(TMath::Pi()*0.4*0.4);
	  //double vall1flat = 1./jecl1flat->getCorrection();
	  double vall1flat = 1./getJEC(jecl1flat,eta,pt);
	  sumvall1flat += w*vall1flat;
	  //
	  //jecl1pt->setJetEta(eta);
	  //jecl1pt->setJetPt(pt);
	  //jecl1pt->setRho(gRho);
	  //jecl1pt->setJetA(TMath::Pi()*0.4*0.4);
	  //double vall1pt = 1./jecl1pt->getCorrection();
	  double vall1pt = 1./getJEC(jecl1pt,eta,pt);
	  sumvall1pt += w*vall1pt;
	  //
	  //jecl1mc->setJetEta(eta);
	  //jecl1mc->setJetPt(pt);
	  //jecl1mc->setRho(gRho);
	  //jecl1mc->setJetA(TMath::Pi()*0.4*0.4);
	  //double vall1mc = 1./jecl1mc->getCorrection();
	  double vall1mc = 1./getJEC(jecl1mc,eta,pt);
	  sumvall1mc += w*vall1mc;
	  //
	  //jecl1pt->setJetEta(eta);
	  //jecl1pt->setJetPt(pt);
	  //jecl1pt->setRho(gRho+1);
	  //jecl1pt->setJetA(TMath::Pi()*0.4*0.4);
	  //double vall1rho = 1./jecl1pt->getCorrection();
	  double vall1rho = 1./getJEC(jecl1pt,eta,pt,gRho+1);
	  sumvall1rho += w*vall1rho;
	  //
	  //
	  //jecl1dt->setJetEta(eta);
	  //jecl1dt->setJetPt(pt);
	  //jecl1dt->setRho(gRho);
	  //jecl1dt->setJetA(TMath::Pi()*0.4*0.4);
	  //double vall1dt = 1./jecl1dt->getCorrection();
	  double vall1dt = 1./getJEC(jecl1dt,eta,pt);
	  sumvall1dt += w*vall1dt;
	  //
	  //jecl1c->setJetEta(eta);
	  //jecl1c->setJetPt(pt);
	  //jecl1c->setRho(gRho);
	  //jecl1c->setJetA(TMath::Pi()*0.4*0.4);
	  //double vall1c = 1./jecl1c->getCorrection();
	  double vall1c = 1./getJEC(jecl1c,eta,pt);
	  sumvall1c += w*vall1c;
	  //
	  //jecl1s->setJetEta(eta);
	  //jecl1s->setJetPt(pt);
	  //jecl1s->setRho(gRho);
	  //jecl1s->setJetA(TMath::Pi()*0.4*0.4);
	  //double vall1s = 1./jecl1s->getCorrection();
	  double vall1s = 1./getJEC(jecl1s,eta,pt);
	  sumvall1s += w*vall1s;
	  //
	  //jecl1rc->setJetEta(eta);
	  //jecl1rc->setJetPt(pt);
	  //jecl1rc->setRho(gRho);
	  //jecl1rc->setJetA(TMath::Pi()*0.4*0.4);
	  //double vall1rc = 1./jecl1rc->getCorrection();
	  double vall1rc = 1./max(0.0001,getJEC(jecl1rc,eta,pt));
	  sumvall1rc += w*vall1rc;
	  //
	  //jecl2c->setJetEta(eta);
	  //jecl2c->setJetPt(pt);
	  //double vall2c = 1./jecl2c->getCorrection();
	  double vall2c = 1./getJEC(jecl2c,eta,pt);
	  sumvall2c += w*vall2c;
	  //
	  //jecl2s->setJetEta(eta);
	  //jecl2s->setJetPt(pt);
	  //double vall2s = 1./jecl2s->getCorrection();
	  double vall2s = 1./getJEC(jecl2s,eta,pt);
	  sumvall2s += w*vall2s;


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
	hl1dtomc->SetBinContent(ipt, (1./l1dt) / (1./l1mc));
	hl1cos->SetBinContent(ipt, (1./l1c) / (1./l1s));
	hl1rcos->SetBinContent(ipt, (1./l1rc) / (1./l1s));
	double l2c = (sumvall2c / sumw);
	double l2s = (sumvall2s / sumw);
	hl2cos->SetBinContent(ipt, (1./l2c) / (1./l2s));
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
      hl1cos->Write();
      hl1rcos->Write();
      hl2cos->Write();
    } // ieta
  } // JEC+sys

  fout->Close();

} // reprocess



