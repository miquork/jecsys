// File: reprocess.C
// Created by Mikko Voutilainen, on Sep 6th, 2012
// Updated on Jan 19, 2015 (Winter14_V8 for the paper)
// Purpose: Combine graphs from difference channels for simpler JEC analysis
//           Macro examplePlot() shows how to create plots with "CMS JEC style"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
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

#include "tdrstyle_mod14.C"

#include <string>
#include <iostream>
#include <vector>
#include <map>

using namespace std;

// rho used to calculate l1bias
const double gRho = 7.1;
const bool _dcsonly = false;

// Minimum pTcut for gamma+jet (raise to 60 GeV for MadGraph samples
double fpptmin(40.);   //76XV1M1:140
double fzeeptmin(30.); //76XV1M1:50
double fzmmptmin(30.); //76XV1M1:50

// put all the different methods in a single file for easy access for everybody
void reprocess() {

  // Set TDR style to have correct graphical setttings when storing graphs
  setTDRStyle();
 
  TDirectory *curdir = gDirectory;

  TFile *fout = new TFile("rootfiles/jecdata.root","RECREATE");

  ////////////////////////
  // Define input files //
  ////////////////////////

  // On 22 Dec 2014, at 14:14, from Rathjens, Denis
  // Re: Out of Office AutoReply: Pseudo Pt plot
  //TFile *fdj = new TFile("rootfiles/Winter14V6.root","READ"); // newL1V6
  //
  // On 31 Jul 2015, at 12:43, Marc Stoever
  // Re: Combination of Run II data sets
  //TFile *fdj = new TFile("rootfiles/combination_DJ_2015-07-31.root","READ");
  // On 05 Aug 2015, at 13:29, Marc Stoever
  // Re: Combination of Run II data sets
  //TFile *fdj = new TFile("rootfiles/combination_DJ_2015-08-05.root","READ");//***
  // On 06 Aug 2015, at 21:22, Marc Stoever
  // Re: Next JERC meeting : Thursday August 6th 18:00
  //TFile *fdj = new TFile("rootfiles/combination_DJ_2015-08-06.root","READ");//badeta
  // On 07 Aug 2015, at 16:41, Marc Stoever
  // Re: Next JERC meeting : Thursday August 6th 18:00
  //TFile *fdj = new TFile("rootfiles/combination_DJ_2015-08-07.root","READ"); //bad!!
  //
  // On 13 Aug 2015, at 19:21, Marc Stoever
  // Re: Dijet combination file
  //TFile *fdj = new TFile("rootfiles/combination_DJ_2015-08-13.root","READ");
  //
  // On 24 Aug 2015, at 15:26, Marc Stoever
  // Re: Dijet combination file
  //TFile *fdj = new TFile("rootfiles/combination_DJ_2015-08-15.root","READ");
  //
  // On 09 Sep 2015, at 19:12, Marc Stoever
  // Re: JSON files: Collisions15 with magnet at 3.8T
  //TFile *fdj = new TFile("rootfiles/combination_DJ_2015-09-09-25ns.root","READ");
  //
  // On 30 Sep 2015, at 18:28, Marc Stoever
  // Re: dijet l2res on 2015D
  //TFile *fdj = new TFile("rootfiles/combination_DJ_2015-09-30-25ns.root","READ");
  //
  // On 20 Oct 2015, at 12:46, Marc Stoever
  // Re: JEC update with 1/fb -- orange alert?
  //TFile *fdj = new TFile(_dcsonly ?
  //			 "rootfiles/CombiFile_Dijet_DConly_974pb.root" :
  //			 "rootfiles/CombiFile_Dijet_goldenJSON_547pb.root",
  //
  // On 27 Oct 2015, at 16:04, Marc Stoever
  // Re: Update with October 19 json file
  TFile *fdj = new TFile("rootfiles/CombiFile_Dijet_Cert_246908-258750.root",
  		 "READ");
  //
  // Anastasia:
  //
  //TFile *fdj = new TFile("rootfiles/JECcombifile_DiJet_1stAttempt.root",
  //	 "READ");

  assert(fdj && !fdj->IsZombie());
  
  // On 14 Jan 2015, at 20:36, from Viola Sordini
  // Re: news from g+j
  // (Renamed 14Jan15 to 17Jan15)
  //TFile *fp = new TFile("rootfiles/gamma_jet_response_all_alphas_with_raw_footprint_regression_17Jan15.root","READ"); // footprint+regression
  //
  // On 03 Aug 2015, at 20:00, Federico Preiato
  // Re: Combination of Run II data sets
  //TFile *fp = new TFile("rootfiles/GammaJet_response_AllMethods_alphacut020_Aug032015.root","READ");
  // On 04 Aug 2015, at 16:26, Federico Preiato
  // Re: Combination of Run II data sets
  //TFile *fp = new TFile("rootfiles/GammaJet_AlphaCut010_020_030_Aug_04_2015.root","READ");
  // On 12 Aug 2015, at 20:42, Federico Preiato
  // Re: Closure Test Residual
  //TFile *fp = new TFile("rootfiles/GammaJet_AlphaCut010_015_020_030_Aug_12_2015.root","READ");
  //
  // On 11 Sep 2015, at 19:18, Federico Preiato
  // Re: File for global fit
  //TFile *fp = new TFile("rootfiles/GammaJet_AllAlphaCut_Sept_11_2015.root","READ");
  //
  // On 26 Sep 2015, at 01:03, Federico Preiato
  // Re: Golden JSON out for Run2015D JEC
  //TFile *fp = new TFile("rootfiles/PhotonJet_Run2015D_DCSONLY_Sept_25_2015.root","READ");
  //
  // On 28 Sep 2015, at 12:46, Federico Preiato
  // Re: Golden JSON out for Run2015D JEC
  //TFile *fp = new TFile("rootfiles/PhotonJet_Run2015D_CertJson_Sept_27_2015.root","READ");
  //
  // On 29 Sep 2015, at 18:53, Federico Preiato
  // Re: New file for global fit
  //TFile *fp = new TFile("rootfiles/PhotonJet_Run2015D_CertJson_Sept_29_2015.root","READ");
  //
  // On 18 Oct 2015, at 15:08, Federico Preiato
  // Re: File for combination
  // On 21 Oct 2015, at 12:57, Federico Preiato
  // Re: a15 in DCSOnly missing?
  //TFile *fp = new TFile(_dcsonly ? 
  //			//"rootfiles/GammaJet_2015_10_18_DCSOnly_980pb.root" :
  //			"rootfiles/GammaJet_2015_10_18_DCSOnly_980pb_Fixed.root":
  //			"rootfiles/GammaJet_2015_10_18_CertifiedJson_552pb.root",
  //
  // On 29 Oct 2015, at 15:24, Federico Preiato
  // Re: a15 in DCSOnly missing?
  //TFile *fp = new TFile("rootfiles/GammaJet_28-10-2015.root",
  //			"READ");
  //
  // On 18 Nov 2015, at 20:13, Federico Preiato
  // Re: New file
  //TFile *fp = new TFile("rootfiles/GammaJet_ppCollision2015_17-11-2015.root",
  //			"READ"); // 2.1 fb-1
  //assert(fp && !fp->IsZombie());
  //
  // On 05 Jan 2016, at 19:20, Federico Preiato
  // Re: New results 
  //TFile *fp = new TFile("rootfiles/GammaJet_Data_vs_GJet_Pythia8_04Jan2016.root","READ"); // 2.1 fb-1 + L2Res + MadGraph
  //TFile *fp = new TFile("rootfiles/GammaJet_Data_vs_GJet_Madgraph_04Jan2016.root","READ"); // 2.1 fb-1 + L2Res + MadGraph
  //assert(fp && !fp->IsZombie());
  //
  // On 12 Jan 2016, at 18:47, Federico Preiato
  // Re: New files
  //TFile *fp = new TFile("rootfiles/GammaJet_Data_vs_GJet_Pythia_2016-01-11.root","READ"); // 129.8/71, 0.982 vs 0.0768
  //TFile *fp = new TFile("rootfiles/GammaJet_Data_vs_GJet_plus_QCD_Pythia_2016-01-11.root","READ"); // 97.1/71, 0.9768 vs 0.0476
  //TFile *fp = new TFile("rootfiles/GammaJet_Data_vs_GJet_plus_QCD_Pythia_2016-01-11.root","READ"); fpptmin = 140; // 68.0/63, 0.9773 vs 0.0249 (best!) => 74X V7
  //TFile *fp = new TFile("rootfiles/GammaJet_Data_vs_GJet_Madgraph_2016-01-11.root","READ"); fpptmin = 60; // 116.3/69, 0.983 vs 0.0789
  //TFile *fp = new TFile("rootfiles/GammaJet_Data_vs_GJet_plus_QCD_Madgraph_2016-01-11.root","READ"); fpptmin = 60; // 116.0/69, 0.9800 vs 0.0771
  //TFile *fp = new TFile("rootfiles/GammaJet_Data_vs_GJet_plus_QCD_Madgraph_2016-01-11.root","READ"); fpptmin = 180; // 74.5/61, 0.9772 vs 0.0245
  //
  // On 15 Jan 2016, at 19:12, Federico Preiato
  // Re: New files
  //TFile *fp = new TFile("rootfiles/GammaJet_Data_vs_GJet_Madgraph_plus_QCD_Pythia_2016-01-11.root","READ"); fpptmin = 60; // 97.4/69, 0.9752 vs 0.0474
  //TFile *fp = new TFile("rootfiles/GammaJet_Data_vs_GJet_Madgraph_plus_QCD_Pythia_2016-01-11.root","READ"); fpptmin = 140; // 70.2/63, 0.9764 vs 0.0260 (2nd)
  // => should still fix pT binning of FSR+ISR corrections, add trk ineff shape
  //
  // On 21 Jan 2016, at 22:48, Federico 
  // Re: global fit, 76X
  //TFile *fp = new TFile("rootfiles/GammaJet_Data_vs_GJet_plus_QCD_Pythia_2016-01-11.root","READ"); fpptmin = 140; // 68.0/63, 0.9773 vs 0.0249 (best!)=>74X V7
  //TFile *fp = new TFile("rootfiles/GammaJet_Data_vs_GJet_Pythia_76X_JECv6_Jan_20_2016.root","READ"); //fpptmin = 140; // 76X
  //
  // On 16 Feb 2016, at 03:25, Federico Preiato
  // Re: global fit, 76X
  //TFile *fp = new TFile("rootfiles/GammaJet_Data_vs_GJetPythia_Rel76X_JECFallV1_Feb_16_2016.root","READ"); // 76X
  //
  // On 18 Feb 2016, at 17:55, Federico Preiato 
  // Re: Combination Binning [Re: Data set comparisons]
  TFile *fp = new TFile("rootfiles/GammaJet_Data_vs_GJetPythia_Rel76X_JECFallV1_Feb_18_2016_NewBins.root","READ"); // 76X V2
  assert(fp && !fp->IsZombie());


  // On 08 Jan 2015, at 14:21, Dominik Haitz
  // Re: Plot updates for the paper
  //TFile *fzee = new TFile("rootfiles/Zee_2015-01-08.root","READ"); // newL1V6
  //assert(fzee && !fzee->IsZombie());

  // On 08 Jan 2015, at 14:21, Dominik Haitz
  // Re: Plot updates for the paper
  //TFile *fzmm = new TFile("rootfiles/Zmm_2015-01-08.root","READ"); // newL1V6
  //
  // On 31 Jul 2015, at 16:23, Fischer, Max (SCC)
  // Re: Combination of Run II data sets
  //TFile *fzmm = new TFile("rootfiles/combination_Zjet_2015-07-31.root","READ"); // eta bins missing
  // On 04 Aug 2015, at 18:49, Fischer, Max (SCC)
  // Re: Combination of Run II data sets
  //TFile *fzmm = new TFile("rootfiles/combination_ZJet_2015-08-04.root","READ"); // buggy MC?
  // On 05 Aug 2015, at 12:14, Dominik Wilhelm Haitz
  // Re: Combination of Run II data sets
  //TFile *fzmm = new TFile("rootfiles/combination_ZJet_2015-08-05.root","READ");
  //
  // On 11 Sep 2015, at 16:58, Fischer, Max 
  // Re: Z+jet for 25 ns
  //TFile *fzmm = new TFile("rootfiles/combination_ZJet_2015-09-11.root","READ");
  //
  // On 19 Oct 2015, at 09:16, Fischer, Max 
  // Re: JEC update with 1/fb -- orange alert?
  //TFile *fzmm = new TFile(_dcsonly ? 
  //			  "rootfiles/combination_ZJet_DCSOnly_2015-10-18.root" :
  //			  "rootfiles/combination_ZJet_2015-10-18.root",
  //
  // On 29 Oct 2015, at 20:06, Fischer, Max
  // Re: [ekp-excalibur] Update with October 19 json file
  //TFile *fzmm = new TFile(_dcsonly ?
  //			  "rootfiles/combination_ZJet_2015-10-29_DCSOnly2015D.root" :
  //			  "rootfiles/combination_ZJet_2015-10-29_golden2015D.root",
  //			  "READ");
  //
  // On 19 Nov 2015, at 11:16, Fischer, Max
  // Re: [ekp-excalibur] Update with November 13 json file
  //TFile *fzmm = new TFile("rootfiles/combination_ZJet_2015-11-19.root","READ");
  //
  // On 01 Dec 2015, at 20:27, Fischer, Max
  // Re: [ekp-excalibur] Update with November 13 json file
  //TFile *fzmm = new TFile("rootfiles/combination_ZJet_2015-12-01.root","READ");
  //
  // On 07 Dec 2015, at 14:27, Fischer, Max (SCC)
  // Re: Zee + Zmm Combination
  //TFile *fzmm = new TFile("rootfiles/combination_ZmmJet_2015-12-07.root","READ");
  //TFile *fzee = new TFile("rootfiles/combination_ZeeJet_2015-12-07.root","READ");
  //
  // On 10 Dec 2015, at 20:51, Fischer, Max 
  // Re: Zee + Zmm Combination
  //TFile *fzmm = new TFile("rootfiles/combination_ZmmJet_2015-12-10.root","READ"); // => 74X V7
  //TFile *fzee = new TFile("rootfiles/combination_ZeeJet_2015-12-10.root","READ"); // => 74X V7
  //
  // On 27 Jan 2016, at 18:20, Fischer, Max (SCC)
  // Re: global fit, 76X
  //TFile *fzmm = new TFile("rootfiles/combination_ZmmJet_2016-01-27.root","READ");
  //TFile *fzee = new TFile("rootfiles/combination_ZeeJet_2016-01-27.root","READ");
  //
  // On 01 Feb 2016, at 22:06, Fischer, Max
  // Re: Corrupt files (was Re: global fit, 76X)
  //TFile *fzmm = new TFile("rootfiles/combination_ZJet_Zmm201602012016-02-01.root","READ");
  //TFile *fzee = new TFile("rootfiles/combination_ZJet_Zee201602012016-02-01.root","READ");
  //
  // On 15 Feb 2016, at 16:45, Fischer, Max (SCC)
  // Re: Corrupt files (was Re: global fit, 76X)
  //TFile *fzmm = new TFile("rootfiles/combination_ZJet_Zmm201602152016-02-15.root","READ");
  //TFile *fzee = new TFile("rootfiles/combination_ZJet_Zee201602152016-02-15.root","READ");
  //
  // Re: Corrupt files (was Re: global fit, 76X)
  // On 17 Feb 2016, at 11:47, Fischer, Max (SCC)
  //TFile *fzmm = new TFile("rootfiles/combination_ZJet_Zmm201602162016-02-16.root","READ");
  //TFile *fzee = new TFile("rootfiles/combination_ZJet_Zee201602162016-02-16.root","READ");
  //
  // Re: Corrupt files (was Re: global fit, 76X)
  // On 18 Feb 2016, at 12:36, Fischer, Max (SCC)
  TFile *fzmm = new TFile("rootfiles/combination_ZJet_Zmm201602182016-02-18.root","READ"); // 76X V2
  TFile *fzee = new TFile("rootfiles/combination_ZJet_Zee201602182016-02-18.root","READ"); // 76X V2


  assert(fzmm && !fzmm->IsZombie());
  assert(fzee && !fzee->IsZombie());
  

  // On 28 Feb 2014, at 14:00, from Dominik Haitz
  // Re: pT~200 GeV
  TFile *feta = new TFile("rootfiles/th2d_jeteta_zpt.root","READ");
  assert(feta && !feta->IsZombie());

  // On 27 Jan 2015, at 14:04, Dominik Haitz
  // Re: 53X electron momentum corrections
  //TFile *fm = new TFile("rootfiles/Zmm-Zee.root","READ");
  //assert(fm && !fm->IsZombie());

  // On 27 Jan 2015, at 16:31, Rajdeep Mohan Chatterjee  
  // Re: Updated L2Res for new L1 - please rerun
  TFile *fm = new TFile("rootfiles/Zee_JEC_V6_ZMassVsZPt.root");
  assert(fm && !fm->IsZombie());
  

  // On 29 Jan 2015, at 14:28, Dominik Haitz
  // Re: 53X electron momentum corrections
  //TFile *fmzee = new TFile("rootfiles/Zee-massvspt.root","READ"); // 74X V7
  TFile *fmzee = fzee; // 76X V2
  assert(fmzee && !fmzee->IsZombie());
  //TFile *fmzmm = new TFile("rootfiles/Zmm-massvspt.root","READ"); // 74X V7
  TFile *fmzmm = fzmm; // 76X V2
  assert(fmzmm && !fmzmm->IsZombie());

  // This is used for correctly averaging JEC and its uncertainty
  // for the wide eta bins used in global fit combinations
  TH2D *h2eta = (TH2D*)feta->Get("data"); assert(h2eta);

  // This is used for scaling Z mass back to 1 for Zee and Zmm
  //TProfile *pzm = (TProfile*)fm->Get(""); assert(pzm);
  //TH1D *hmzee = (TH1D*)fmzee->Get("zmass_zpt_Zee_ratio"); assert(hmzee); // 74X V7
  //TH1D *hmzmm = (TH1D*)fmzmm->Get("zmass_zpt_Zmm_ratio"); assert(hmzmm); // 74X V7
  TH1D *hmzee = (TH1D*)fmzee->Get("Ratio_ZMass_CHS_a30_eta_00_13_L1L2L3"); // 76X V2
  TH1D *hmzmm = (TH1D*)fmzmm->Get("Ratio_ZMass_CHS_a30_eta_00_13_L1L2L3"); // 76X V2
  assert(hmzee);
  assert(hmzmm);

  map<string, TFile*> files;
  files["dijet"] = fdj;
  files["gamjet"] = fp;
  files["zeejet"] = fzee;
  files["zmmjet"] = fzmm;


  // Map variable names used in different files to common scheme
  map<string, map<string, const char*> > rename;

  // Note: non-CHS results are not available for all samples, so not used
  //rename["dijet"]["ratio"] = "ratio"; 
  //rename["dijet"]["data"] = "data"; 
  //rename["dijet"]["mc"] = "MC"; 
  //rename["dijet"]["mpfchs"] = "MPFchs_VsPtAve"; 
  //rename["dijet"]["mpfchs1"] = "MPFT1chs_VsPtAve";
  //rename["dijet"]["ptchs"] = "PtBalchs";
  rename["dijet"]["mpfchs"] = "mpfchs";
  rename["dijet"]["mpfchs1"] = "mpfchs";
  rename["dijet"]["ptchs"] = "ptchs";

  rename["gamjet"]["ratio"] = "";
  rename["gamjet"]["data"] = "_DATA"; 
  rename["gamjet"]["mc"] = "_MC"; 
  rename["gamjet"]["mpfchs"] = "resp_MPFchs";
  rename["gamjet"]["mpfchs1"] = "resp_MPFchs"; 
  rename["gamjet"]["ptchs"] = "resp_PtBalchs"; 

  /*
  rename["zeejet"]["ratio"] = "ratio"; 
  rename["zeejet"]["data"] = "Data_Ee_Corr"; 
  rename["zeejet"]["mc"] = "MC_Ee_Corr"; 
  rename["zeejet"]["mpfchs"] = "MPF-notypeI_CHS";
  rename["zeejet"]["mpfchs1"] = "MPF_CHS";
  rename["zeejet"]["ptchs"] = "PtBal_CHS"; 
  */

  rename["zeejet"]["ratio"] = "Ratio";
  rename["zeejet"]["data"] = "Data";
  rename["zeejet"]["mc"] = "MC";
  rename["zeejet"]["mpfchs"] = "MPF-notypeI_CHS";
  rename["zeejet"]["mpfchs1"] = "MPF_CHS"; 
  rename["zeejet"]["ptchs"] = "PtBal_CHS"; 

  rename["zmmjet"]["ratio"] = "Ratio";//"ratio"; 
  rename["zmmjet"]["data"] = "Data";
  rename["zmmjet"]["mc"] = "MC";
  rename["zmmjet"]["mpfchs"] = "MPF-notypeI_CHS";
  rename["zmmjet"]["mpfchs1"] = "MPF_CHS"; 
  rename["zmmjet"]["ptchs"] = "PtBal_CHS"; 


  

  // color and style codes
  //map<string, int> style;
  map<string, map<string, int> > style;
  //style["mpfchs"] = kFullStar;
  //style["mpf"] = kFullDiamond;
  //style["mpfchs1"] = kFullCircle;
  //style["mpf1"] = kFullSquare;
  //style["ptchs"] = kOpenCircle;
  //style["pt"] = kOpenSquare;
  //style["mpfchs1"] = kFullCircle;
  //style["ptchs"] = kOpenSquare;

  style["dijet"]["mpfchs1"] = kFullCircle;
  style["dijet"]["ptchs"] = kOpenSquare;
  style["gamjet"]["mpfchs1"] = kFullSquare;
  style["gamjet"]["ptchs"] = kOpenSquare;
  style["zeejet"]["mpfchs1"] = kFullCircle;
  style["zeejet"]["ptchs"] = kOpenCircle;
  style["zmmjet"]["mpfchs1"] = kFullStar;//kFullDiamond;
  style["zmmjet"]["ptchs"] = kOpenStar;//kOpenDiamond;

  map<string, int> color;
  color["dijet"] = kBlack;
  color["gamjet"] = kBlue;
  color["zeejet"] = kGreen+2;
  color["zmmjet"] = kRed;

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
  //types.push_back("mpfchs");
  types.push_back("mpfchs1");
  types.push_back("ptchs");

  vector<string> sets;
  sets.push_back("dijet");
  sets.push_back("gamjet");
  sets.push_back("zeejet");
  sets.push_back("zmmjet");

  vector<pair<double,double> > etas;
  // reference region |eta|<1.3
  etas.push_back(make_pair<double,double>(0,1.305));
  // bins for coarse eta-dependent corrections
  etas.push_back(make_pair<double,double>(0,0.783)); // do these exist?
  etas.push_back(make_pair<double,double>(0.783,1.305)); // do these exist?
  etas.push_back(make_pair<double,double>(1.305,1.93));
  etas.push_back(make_pair<double,double>(1.93,2.5));
  etas.push_back(make_pair<double,double>(2.5,2.964));
  etas.push_back(make_pair<double,double>(2.964,3.2));
  etas.push_back(make_pair<double,double>(3.2,5.191));

  vector<double> alphas;
  alphas.push_back(0.10);
  alphas.push_back(0.15);
  alphas.push_back(0.20);
  alphas.push_back(0.30);

  ///////////////////////////////////////////
  // Rename selected graphs and store them //
  ///////////////////////////////////////////

  map<string, map<string, map<string, map<int, map<int, TGraphErrors*> > > > > grs;

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
	  if (s=="gamjet" && f==0) {
	    if (t=="mpfchs" || t=="mpfchs1") f = fp;//fp1;
	    if (t=="ptchs") f = fp;//fp2;
	  }
	  assert(f);

	  for (unsigned int ialpha = 0; ialpha != alphas.size(); ++ialpha) {

	    double alpha = alphas[ialpha];

	    eta1 = etas[ieta].first; // reset to avoid trouble with below
	    eta2 = etas[ieta].second; // reset to avoid trouble with below

	    // If samples missing break-up into e.g. [3,3.2] and [3.2,5.2] bins
	    // or merged [0,1.3] bin, patch here
	    //if (s=="zeejet" && eta1==0 && fabs(eta2-1.3)<0.1) { eta2=1.305; }
	    //if (s=="zmmjet" && fabs(eta1-3.0)<0.1) { eta2=3.1; }
	    //if (s=="zmmjet" && fabs(eta1-3.2)<0.1) { eta1=3.1; }
	    if (s=="dijet"  && fabs(eta1-3.2)<0.1) { eta2=5.0; } // 74X
	    //if (s=="gamjet"  && fabs(alpha-0.3)<0.1) { alpha=0.2; }
	    //if (s=="gamjet"  && fabs(eta1-0.0)<0.1) { eta2=1.3; } // ??
	    //if (s=="gamjet"  && fabs(eta1-0.8)<0.1) { eta1=0.0; } // ??
	    //if (s=="gamjet"  && fabs(eta1-1.3)<0.1) { eta2=2.0; }
	    //if (s=="gamjet"  && fabs(eta1-1.9)<0.1) { eta1=2.0; eta2=3.0; }
	    //if (s=="gamjet"  && fabs(eta1-2.5)<0.1) { eta1=2.0; eta2=3.0; }
	    //if (s=="gamjet"  && fabs(eta1-3.0)<0.1) { eta2=5.2; }
	    //if (s=="gamjet"  && fabs(eta1-3.2)<0.1) { eta1=3.0; }

	    // If some subsets of data missing, patch (skip) here
	    // gamjet missing non-CHS graphs
	    //if (s=="gamjet" && (t=="mpf" || t=="mpf1" || t=="pt"))
	    //continue;
	    //if (s=="dijet" && (fabs(eta1-0.783)<0.1 || fabs(eta2-0.783)<0.1))
	    //continue;

	    // Reconstruct naming scheme used in each of the files
	    // If non-conventional naming schemes, patch here
	    const char *c(0);
	    if (s=="dijet") {
	      //c = Form("%s_a%1.0f_eta%1.0f_%1.0f_%s",
	      //       rename[s][t], 100.*alpha, 10.*eta1, 10.*eta2,
	      //       rename[s][d]);
	      c = Form("%s/eta%02.0f-%02.0f/%s_%s_a%1.0f", // 74X
	             dd, 10.*eta1, 10.*eta2, rename[s][t], ss, 100.*alpha);
	      //c = Form("%s/eta%04.0f-%04.0f/%s_%s_a%1.0f", // 76X
		//       dd, 1000.*eta1, 1000.*eta2, rename[s][t], ss, 100.*alpha);
	    } // dijet
	    if (s=="gamjet") {
	      //if (t=="mpfchs") // special *_raw_* naming for non-type-I MPF
		//c = Form("%s%s_raw_a%1.0f_eta%02.0f_%02.0f",
	      // rename[s][t], rename[s][d],
	      //	 100.*alpha, 10.*eta1, 10.*eta2);
	      //else
	      c = Form("%s%s_a%1.0f_eta%02.0f_%02.0f",
		       rename[s][t], rename[s][d],
		       100.*alpha, 10.*eta1, 10.*eta2);
	    } // gamjet
	    //if (s=="zeejet" ) {
	    //c = Form("%s_%s_a%1.0f_eta%02.0f_%02.0f_L1L2L3Res",
	    //	       rename[s][d], rename[s][t], 100.*alpha,
	    //	       10.*eta1, 10.*eta2);
	    //} // zeejet
	    if (s=="zmmjet" || s=="zeejet") {
	      //c = Form("%s_%s_a%1.0f_eta%02.0f_%02.0f_L1L2L3Res",
	      c = Form("%s_%s_a%1.0f_eta_%02.0f_%02.0f_L1L2L3",
		       rename[s][d], rename[s][t], 100.*alpha,
		       10.*eta1, 10.*eta2);
		       
	    } // zmmjet
	    assert(c);

	    TObject *obj = f->Get(c);
	    if (!obj) {
	      cout << "Graph " << c << " not found for "
		   << s << " " << t << "!" <<endl << flush;
	      cout << "File: " << f->GetName() << endl << flush;
	    }
	    
	    assert(obj);

	    // If data stored in TH1D instead of TGraphErrors, patch here
	    // Patch for dijet file that has TH1D's instead of graphs
	    if (obj->InheritsFrom("TH1")) {
	      obj = new TGraphErrors((TH1D*)obj);
	    }

	    assert(obj->InheritsFrom("TGraphErrors"));
	    TGraphErrors *g = (TGraphErrors*)obj;

	    // Clean out empty points from TH1D->TGraphErrors conversion
	    for (int i = g->GetN()-1; i != -1; --i) {
	      if (g->GetY()[i]==0 && g->GetEY()[i]==0) g->RemovePoint(i);
	      // clean out weird points from Sebastien's data
	      //else if (g->GetX()[i]>4000 || g->GetX()[i]<10) g->RemovePoint(i);
	      // Remove too large uncertainties
	      //else if (g->GetEY()[i]>0.05 && g->GetX()[i]>1000.) g->RemovePoint(i); // 74X
	      //else if (g->GetEY()[i]>0.03 && g->GetX()[i]>1000.) g->RemovePoint(i); // 76x
	      else if (g->GetEY()[i]>0.03 && g->GetX()[i]>400.) g->RemovePoint(i); // 76x
	      // Z+jet points with too wide bin at high pT
	      //else if ((s=="zmmjet" || s=="zeejet") &&
	      else if ((s=="zeejet" && t=="ptchs") &&
		       (g->GetX()[i]>300 && g->GetEX()[i]>100 &&
			fabs(eta1-0.0)<0.1))
		g->RemovePoint(i);
	      else if (s=="gamjet" && g->GetX()[i]<fpptmin)
		g->RemovePoint(i);
	      else if (s=="zeejet" && g->GetX()[i]<fzeeptmin)
		g->RemovePoint(i);
	      else if (s=="zmmjet" && g->GetX()[i]<fzmmptmin)
		g->RemovePoint(i);
	      // Remove bad point with too small uncertainty
	      //else if (s=="dijet" &&
	      //       ((fabs(eta1-1.3)<0.1 && fabs(g->GetX()[i]-200)<10) ||
	      //	(fabs(eta1-1.9)<0.1 && fabs(g->GetX()[i]-200)<10) ||
	      //	(fabs(eta1-2.5)<0.1 && fabs(g->GetX()[i]-200)<10)))
	      //g->RemovePoint(i);
	    } // for i

	    // patch Z+jet pT center and uncertainty
	    if ((s=="zeejet" || s=="zmmjet") && d=="ratio") {

	      TGraphErrors *gfix = grs["mc"][t][s][ieta][ialpha];
              assert(gfix);
	      for (int i = 0; i != g->GetN(); ++i) {
		double x = g->GetX()[i];
		double ex = g->GetEX()[i];
		for (int j = 0; j != gfix->GetN(); ++j) {
		  double xx = gfix->GetX()[j];
		  double exx = gfix->GetEX()[j];
		  if (fabs(x-xx)<ex) {
		    g->SetPoint(i, xx, g->GetY()[i]);
		    g->SetPointError(i, exx, g->GetEY()[i]);
		  }
		} // for j
	      } // for i
	    } // patch <pT>

	    // chi2/NDF=126.0/71 before Zmm+Zee scales (gam pT>40)
	    // chi2/NDF=79.7/63 before Zmm+Zee scales (gam pT>140)
	    // chi2/NDF=126.8/71 with Zee scale (no Zmm scale, gam pT>40)
	    if (false && s=="zeejet" && (d=="data" || d=="ratio")) {
	      for (int i = 0; i != g->GetN(); ++i) {
		double ipt = hmzee->FindBin(g->GetX()[i]);
		double k = hmzee->GetBinContent(ipt);
		double ek = hmzee->GetBinError(ipt);
		g->SetPoint(i, g->GetX()[i], g->GetY()[i]*k);
 		g->SetPointError(i, g->GetEX()[i], 
				 sqrt(pow(g->GetEY()[i]*k,2) + ek*ek));
	      }
	    }
	    // chi2/NDF=122.7/71 with Zmm scale (no Zee scale, gam pT>40)
	    // chi2/NDF=123.4/71 with Zmm+Zee scales (gam pT>40)
	    // chi2/NDF=77.4/63 with Zmm scale (no Zee scale, gam pT>140)
	    // chi2/NDF=78.6/63 with Zmm+Zee scale (gam pT>140)
 	    if (false && s=="zmmjet" && (d=="data" || d=="ratio")) {
 	      for (int i = 0; i != g->GetN(); ++i) {
 		double ipt = hmzmm->FindBin(g->GetX()[i]);
 		double k = hmzmm->GetBinContent(ipt);
 		double ek = hmzmm->GetBinError(ipt);
 		g->SetPoint(i, g->GetX()[i], g->GetY()[i]*k);
 		g->SetPointError(i, g->GetEX()[i], 
				 sqrt(pow(g->GetEY()[i]*k,2) + ek*ek));
 	      }
 	    }

	    // TEST 2015-01-29: remove first two points at pT<60 GeV
	    // from pTbal for: (i) photon+jet (-2.8/2), (ii) Zmm+jet (-9.2/3),
	    // (iii) Zee+jet (-7.3 / 3)
	    // The best we have with Z mass fixes plus cutting out pTbal<60 GeV
	    // is 87.5 / 84, compared to 111.5 / 92 otherwise
	    if (false && (s=="gamjet" || s=="zmmjet" || s=="zeejet")
		&& t=="ptchs") {
	      for (int i = g->GetN()-1; i != -1; --i) {
		if (g->GetX()[i]<60.) g->RemovePoint(i);
	      }
	    }

	    dout->cd();

	    // Set uniform naming scheme and graphical style
	    //g->SetName(Form("%s_%s_a%1.0f",tt,ss,100.*alpha));
	    g->SetName(Form("%s_%s_a%1.0f",tt,ss,100.*alphas[ialpha]));
	    g->UseCurrentStyle(); // Basic TDR style
	    //g->SetMarkerStyle(style[t]);
	    g->SetMarkerStyle(style[s][t]);
	    g->SetMarkerColor(color[s]);
	    g->SetLineColor(color[s]);
	    g->SetMarkerSize(size[alpha]);
	    g->SetDrawOption("SAMEP");

	    g->Write();

	    grs[d][t][s][ieta][ialpha] = g;

	  } // for ialpha
	} // for iset
      } // for itype
    } // for ieta
  } // for itier

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
    
    // New JEC for plotting on the back
    FactorizedJetCorrector *jec;
    {
      //s = Form("%s/Winter14_V8_DATA_L2L3Residual_AK5PFchs.txt",cd);
      //s = Form("%s/Summer15_50ns_DATA_L2L3Residual_AK4PFchs.txt",cd);
      //s = Form("%s/Summer15_50nsV1_DATA_L2L3Residual_AK4PFchs.txt",cd);
      //s = Form("%s/Summer15_50nsV2_DATA_L2L3Residual_AK4PFchs.txt",cd);
      //s = Form("%s/Summer15_50nsV2_DATA_L2L3Residual_AK4PFchs.txt.LOGLIN",cd);
      //s = Form("%s/Summer15_50nsV3_DATA_L2L3Residual_AK4PFchs.txt",cd);
      //s = Form("%s/Summer15_50nsV3M1_DATA_L2L3Residual_AK4PFchs.txt",cd);
      //s = Form("%s/Summer15_50nsV3M2_DATA_L2L3Residual_AK4PFchs.txt",cd);
      //s = Form("%s/Summer15_50nsV4_DATA_L2L3Residual_AK4PFchs.txt",cd); // V4
      //s = Form("%s/Summer15_50nsV3M2_DATA_L2L3Residual_AK4PFchs.txt.LOGLIN",cd);
      //s = Form("%s/Summer15_25nsV2H1M1_DATA_L2L3Residual_AK4PFchs.txt",cd);
      //s = Form("%s/Summer15_25nsV3_DATA_L2L3Residual_AK4PFchs.txt",cd);
      //s = Form("%s/Summer15_25nsV3M2_DATA_L2L3Residual_AK4PFchs.txt",cd);
      //s = Form("%s/Summer15_25nsV3M3_DATA_L2L3Residual_AK4PFchs.txt",cd);
      //s = Form("%s/Summer15_25nsV6M1_DATA_L2L3Residual_AK4PFchs.txt",cd);
      //s = Form("%s/Summer15_25nsV6M2_DATA_L2L3Residual_AK4PFchs.txt",cd);
      //s = Form("%s/Summer15_25nsV6M3_DATA_L2L3Residual_AK4PFchs.txt",cd);
      //s = Form("%s/Summer15_25nsV6_DATA_L2L3Residual_AK4PFchs.txt",cd);
      //s = Form("%s/Summer15_25nsV7M1_DATA_L2L3Residual_AK4PFchs.txt",cd); // 74X V7
      //s = Form("%s/Fall15_25nsV1M1_DATA_L2L3Residual_AK4PFchs.txt",cd);
      s = Form("%s/Fall15_25nsV1M2_DATA_L2L3Residual_AK4PFchs.txt",cd); // 76X V2
      cout << s << endl;
      JetCorrectorParameters *par_l2l3res = new JetCorrectorParameters(s);
      vector<JetCorrectorParameters> vpar;
      vpar.push_back(*par_l2l3res);
      jec = new FactorizedJetCorrector(vpar);
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
    // But event with this pT-dependent L2Res can cause problems
    FactorizedJetCorrector *jecold;
    {
      //s = Form("%s/Winter14_V6_DATA_L2L3Residual_AK5PFchs.txt",cd); // 74X V7
      s = Form("%s/Fall15_25ns_COMB_LOGLIN_L2Residual_v2_AK4PFchs_nokFSR.txt",cd); // 76X V2
      cout << s << endl;
      JetCorrectorParameters *par_old = new JetCorrectorParameters(s);
      vector<JetCorrectorParameters> v;
      v.push_back(*par_old);
      jecold = new FactorizedJetCorrector(v);
    }

    // Difference between pT-dependent and flat L1
    FactorizedJetCorrector *jecl1flat;
    {
      //s = Form("%s/Winter14_V0_DATA_L1FastJetPU_AK5PFchs_pt.txt",cd); // 74X V7
      s = Form("%s/Fall15_25nsV1_DATA_L1RC_AK4PFchs.txt",cd); // 76X V2
      cout << s << endl << flush;
      JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
      vector<JetCorrectorParameters> v;
      v.push_back(*l1);
      jecl1flat = new FactorizedJetCorrector(v);
    }
    FactorizedJetCorrector *jecl1pt;
    {
      //s = Form("%s/Winter14_V1_DATA_L1FastJet_AK5PFchs.txt",cd); // GT
      //s = Form("%s/Winter14_V6_DATA_L1FastJet_AK5PFchs.txt",cd); // newL1V6
      //s = Form("%s/Winter14_V8_DATA_L1FastJet_AK5PFchs.txt",cd); // 74X V7
      s = Form("%s/Fall15_25nsV1_DATA_L1FastJet_AK4PFchs.txt",cd); // 76X V2
      cout << s << endl << flush;
      JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
      vector<JetCorrectorParameters> v;
      v.push_back(*l1);
      jecl1pt = new FactorizedJetCorrector(v);
    }

    // Run I uncertainty
    s = Form("%s/Winter14_V10M_DATA_UncertaintySources_AK5PFchs.txt",cd); // V8
    s2 = "TotalNoFlavorNoTime";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_ref1 = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc_ref1 = new JetCorrectionUncertainty(*p_ref1);

    // Total uncertainty, excluding Flavor and Time
    //s = Form("%s/Winter14_V10M_DATA_Uncertainty_AK5PFchs.txt",cd); // V8 Run I
    //s = Form("txt/Summer15_50nsV2_DATA_UncertaintySources_AK4PFchs.txt",cd); // Run II
    //s = Form("txt/Summer15_50nsV3_DATA_UncertaintySources_AK4PFchs.txt",cd);
    //s = Form("txt/Summer15_50nsV3M2_DATA_UncertaintySources_AK4PFchs.txt",cd); // Run II
    //s = Form("%s/Summer15_50nsV4_DATA_UncertaintySources_AK4PFchs.txt",cd); // Run II
    //s = Form("txt/Summer15_50nsV4M1_DATA_UncertaintySources_AK4PFchs.txt",cd); // Run II
    //s = Form("txt/Summer15_25nsV3M3_DATA_UncertaintySources_AK4PFchs.txt",cd); // Run II
    //s = Form("%s/Summer15_25nsV3M3_DATA_UncertaintySources_AK4PFchs.txt",cd); // Run II
    //s = Form("%s/Summer15_25nsV6M1_DATA_UncertaintySources_AK4PFchs.txt",cd);
    //s = Form("%s/Summer15_25nsV6M2_DATA_UncertaintySources_AK4PFchs.txt",cd);
    //s = Form("%s/Summer15_25nsV6M3_DATA_UncertaintySources_AK4PFchs.txt",cd);
    //s = Form("%s/Summer15_25nsV6_DATA_UncertaintySources_AK4PFchs.txt",cd);
    //s = Form("%s/Summer15_25nsV7M1_DATA_UncertaintySources_AK4PFchs.txt",cd); // 74X V7
    s = Form("%s/Fall15_25nsV1M2_DATA_UncertaintySources_AK4PFchs.txt",cd); // 76X V2
    s2 = "TotalNoFlavorNoTime";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_unc = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p_unc);

    // Partial uncertainties
    //s = Form("%s/Winter14_V10M_DATA_UncertaintySources_AK5PFchs.txt",cd); // V8
    //s = Form("txt/Summer15_50nsV2_DATA_UncertaintySources_AK4PFchs.txt",cd); // Run II
    //s = Form("txt/Summer15_50nsV3_DATA_UncertaintySources_AK4PFchs.txt",cd);
    //s = Form("txt/Summer15_50nsV3M2_DATA_UncertaintySources_AK4PFchs.txt",cd); // Run II
    //s = Form("%s/Summer15_50nsV4_DATA_UncertaintySources_AK4PFchs.txt",cd); // Run II
    //s = Form("txt/Summer15_50nsV4M1_DATA_UncertaintySources_AK4PFchs.txt",cd); // Run II
    //s = Form("txt/Summer15_25nsV3M3_DATA_UncertaintySources_AK4PFchs.txt",cd); // Run II
    //s = Form("%s/Summer15_25nsV3M3_DATA_UncertaintySources_AK4PFchs.txt",cd); // Run II
    //s = Form("%s/Summer15_25nsV6M1_DATA_UncertaintySources_AK4PFchs.txt",cd);
    //s = Form("%s/Summer15_25nsV6M2_DATA_UncertaintySources_AK4PFchs.txt",cd);
    //s = Form("%s/Summer15_25nsV6M3_DATA_UncertaintySources_AK4PFchs.txt",cd);
    //s = Form("%s/Summer15_25nsV6_DATA_UncertaintySources_AK4PFchs.txt",cd);
    //s = Form("%s/Summer15_25nsV7M1_DATA_UncertaintySources_AK4PFchs.txt",cd);
    //s = Form("%s/Summer15_25nsV7M1_DATA_UncertaintySources_AK4PFchs.txt",cd); // 74X V7
    s = Form("%s/Fall15_25nsV1M2_DATA_UncertaintySources_AK4PFchs.txt",cd); // 76X V2
    s2 = "TotalNoFlavorNoTime";
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

      
      const double ptbins[] = {29,30, 40, 50, 60, 75, 100, 125, 155, 180,
			       210, 250, 300, 350, 400, 500, 600, 800, 1000,
			       1200, 1500, //1600, 1601};
			       1800, 2100, 2500, 3000, 3500, 3501};
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

      // PileUpPtEnvelope, effectively
      TH1D *hl1bias = new TH1D("hl1bias",";p_{T} (GeV);"
			       "JEC L1 p_{T} / JEC L1 flat;",
			       npt, &ptbins[0]);

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
	double sumrun1(0), sumjes(0), sumvall1flat(0), sumvall1pt(0);
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
	  jec->setJetEta(eta);
	  jec->setJetPt(pt);
	  double val = 1./jec->getCorrection();
	  sumval += w*val;
	  sumw += w; // sum weights only once

	  // reference JEC
	  jecrun1->setJetEta(eta);
	  jecrun1->setJetPt(pt);
	  double jesrun1 = 1./jecrun1->getCorrection();
	  sumrun1 += w*jesrun1;

	  // old JEC
	  jecold->setJetEta(eta);
	  jecold->setJetPt(pt);
	  //double jes = 1.; // 74X V7
	  double jes = 1./jecold->getCorrection(); // 76X V2
	  sumjes += w*jes;

	  jecl1flat->setJetEta(eta);
	  jecl1flat->setJetPt(pt);
	  //jecl1flat->setNPV(int(14.4));
	  jecl1flat->setRho(gRho);
	  jecl1flat->setJetA(TMath::Pi()*0.4*0.4);//0.5*0.5);
	  double vall1flat = 1./jecl1flat->getCorrection();
	  sumvall1flat += w*vall1flat;
	  //
	  jecl1pt->setJetEta(eta);
	  jecl1pt->setJetPt(pt);
	  //jecl1pt->setNPV(int(14.4));
	  jecl1pt->setRho(gRho);
	  jecl1pt->setJetA(TMath::Pi()*0.4*0.4);//0.5*0.5);
	  double vall1pt = 1./jecl1pt->getCorrection();
	  sumvall1pt += w*vall1pt;

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
	hrun1->SetBinContent(ipt, run1);
	hrun1->SetBinError(ipt, run1*err_ref1);//run1*err);

	double jes = (sumjes / sumw);
	hjes->SetBinContent(ipt, jes);
	double l1pt = (sumvall1pt / sumw);
	double l1flat = (sumvall1flat / sumw);
	double l1bias = l1pt / l1flat;
	hl1bias->SetBinContent(ipt, l1bias);
	hl1bias->SetBinError(ipt, 0.5*fabs(1-l1bias));
      } // ipt

      dout2->cd();

      herr->SetMarkerSize(0);
      herr->SetFillStyle(1001);
      herr->SetFillColor(kYellow+1);
      herr->Write();

      herr_ref->SetMarkerSize(0);
      herr_ref->SetFillStyle(1001);
      herr_ref->SetFillColor(kYellow+1);
      herr_ref->Write();

      herr_spr->SetMarkerSize(0);
      herr_spr->SetFillStyle(1001);
      herr_spr->SetFillColor(kYellow);
      herr_spr->Write();

      herr_pu->SetMarkerSize(0);
      herr_pu->SetFillStyle(1001);
      herr_pu->SetFillColor(kYellow+1);
      herr_pu->Write();

      herr_noflv->SetMarkerSize(0);
      herr_noflv->SetFillStyle(1001);
      herr_noflv->SetFillColor(kYellow+1);
      herr_noflv->Write();

      herr_mpf->SetMarkerSize(0);
      herr_mpf->SetFillStyle(1001);
      herr_mpf->SetFillColor(kYellow+1);
      herr_mpf->Write();

      herr_pt->SetMarkerSize(0);
      herr_pt->SetFillStyle(1001);
      herr_pt->SetFillColor(kYellow+1);
      herr_pt->Write();

      hrun1->SetMarkerSize(0);
      hrun1->SetFillStyle(1001);
      hrun1->SetFillColor(kCyan+1);
      hrun1->Write();

      hjes->Write();
      hl1bias->Write();
    } // ieta
  } // JEC+sys

  fout->Close();

} // reprocess
