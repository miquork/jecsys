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
const double gRho = 15.36; // 2017-06-02 for 2016 data
const bool _dcsonly = false;
bool doFall17Closure = false;//true;

// Appy mass corrections to Z+jet
bool useFixedFit = true; // with minitools/drawZmass.C
double fitUncMin = 0.00000; // Some bug if unc<0?
bool correctZmmMass = true; // pT with mumu mass
bool correctZeeMass = true; // pT with ee mass
bool correctGamMass = true; // pT with ee mass at 2*pT
bool correctUncert = true;  // ll mass uncertainty
//
bool correctGamScale = false;     // separate fixed scale with value below
double valueGamScale = 1.;//1./0.989;  // drawGamVsZmm BCDEFGH

// Minimum pTcut for gamma+jet
double fpmpfptmin(100.); // photon+jet MPF
double fpbalptmin(100.); // photon+jet pTbal
double fzeeptmin(30.);   // Zee+jet both methods
double fzmmptmin(30.);   // Zmm+jet both methods
// Additional cuts to Z+jet MPF / balance methods
double fzmpfptmin(30.);   // Z+jet MPF
double fzbalptmin(100);//30.);   // Z+jet pTbal

//for fine etabins deactivate ptbal
double fdijetmpfptmin(30);
double fdijetbalptmin(30.);
double fdijetptmax(1500.);

// Maximum pTcut for samples (to avoid bins with too large uncertainty)
double fpmpfptmax(1500.);//2500.); // photon+jet MPF
double fpbalptmax(700.);//2500.);//700.);  // photon+jet pTbal
double fzeeptmax(700.);   // Zee+jet
double fzmmptmax(700.);   // Zmm+jet
// Additional cuts to Z+jet MPF / balance methods
double fzmpfptmax(500.);//700.);  // Z+jet MPF
double fzbalptmax(300.);//700.);  // Z+je pTbal

//minimum event counts
const double neventsmin = 20.;

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
  TFile *fdj = new TFile(Form("rootfiles/L2Res_Fall17_17Nov_V5V6/Run%s/JEC_L2_Dijet_AK4PFchs_pythia8.root",fdj_files[epoch]),"READ");

  // Anastasia Karavdina, 2016 Legacy re-reco (8 Dec 2017) :
  // https://indico.cern.ch/event/682570/
  //TFile *fdj = new TFile("rootfiles/JEC_L2_Dijet_AK4PFchs_pythia8_07122017hcalCleaning_wideEtabins.root"); // BCDEFGH only?
  //assert(fdj && !fdj->IsZombie());

  //TFile *fdj2 = new TFile("rootfiles/JEC_L2_Dijet_AK4PFchs_pythia8_07122017hcalCleaning.root"); // narrow bins to complement above wide ones
  //assert(fdj2 && !fdj2->IsZombie());

  // Andrey Popov, April 28, 2017 (Feb03_L2ResV2)
  // https://indico.cern.ch/event/634367/
  map<string,const char*> fm_files;
  fm_files["BCD"] = "0428_Run2016BCD"; // also update multijet.C
  fm_files["EF"] = "0428_Run2016EFearly"; // also update multijet.C
  fm_files["G"] = "0428_Run2016FlateG"; // also update multijet.C
  fm_files["H"] = "0428_Run2016H"; // also update multijet.C
  fm_files["GH"] = "0428_Run2016H"; // duplicate H
  fm_files["BCDEFGH"] = "0428_Run2016All"; // also update multijet.C
  //TFile *fmj = new TFile(Form("rootfiles/multijet_2017%s.root",
  //		      fm_files[epoch]),"READ");

  //assert(fmj && !fmj->IsZombie());
  
  // Hugues Lattaud, 2017 V27 inputs with L2Res (not incl. JER SF as in V28)
  // https://indico.cern.ch/event/765393/#47-l3res-gammajets-with-fall17
  map<string,const char*> fp_files;
  fp_files["B"] = "B_B";
  fp_files["C"] = "C_C";
  fp_files["D"] = "D_D";
  fp_files["E"] = "E_E";
  fp_files["F"] = "F_F";
  fp_files["BCDEF"] = "BCDEF";
  TFile *fp = new TFile(Form("rootfiles/Gjet_combinationfile_17_Nov_V27_L2res_%s.root", fp_files[epoch]),"READ");

  if (doFall17Closure){
    // Hugues Lattaud, 2017 V24 closure, Sep 17, 2018
    // https://indico.cern.ch/event/758051/
    fp_files["B"] = "B";
    fp_files["C"] = "C";
    fp_files["D"] = "D";
    fp_files["E"] = "E";
    fp_files["F"] = "F";
    fp_files["BCDEF"] = "BCDEF";
    fp = new TFile(Form("rootfiles/Gjet_combinationfile_17_Nov_V24_L2L3Res_%s.root", fp_files[epoch]),"READ");
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
  fz_files["BCDEF"] = "BCDEF";
  TFile *fzmm = new TFile(Form("rootfiles/zjet_combination_Fall17_JECV28_Zmm_%s_2018-10-08.root",fz_files[epoch]),"READ");
  TFile *fzee = new TFile(Form("rootfiles/zjet_combination_Fall17_JECV28_Zee_%s_2018-10-08.root",fz_files[epoch]),"READ");
  assert(fzmm && !fzmm->IsZombie());
  assert(fzee && !fzee->IsZombie());

  if (doFall17Closure){
    // Daniel Savoiu, 2017 V24 closure, 17 Sep 2018
    // https://indico.cern.ch/event/758051/
    fzmm = new TFile(Form("rootfiles/zjet_combination_Fall17_JECV24_Zmm_%s_2018-09-18.root",fz_files[epoch]),"READ");
    fzee = new TFile(Form("rootfiles/zjet_combination_Fall17_JECV24_Zee_%s_2018-09-18.root",fz_files[epoch]),"READ");
    assert(fzmm && !fzmm->IsZombie());
    assert(fzee && !fzee->IsZombie());
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
      // BCDEFGH fit with minitools/drawZmass.C
      //f1mzee->SetParameters(1.01732, -0.00909, 0.00116);
      //f1ezee->SetParameters(7.54e-05, 1.41e-05, 1.63e-07,
      //		    -3.26e-05, 3.47e-06, -1.51e-06);
      // Fall17 BCDEF
      //f1mzee->SetParameters(1.02106, -0.01001, 0.00130);
      //f1ezee->SetParameters(1.02e-04, 1.83e-05, 2.01e-07,
      //		    -4.31e-05, 4.48e-06, -1.91e-06);
      // 17Nov17_V10 BCDEF
      f1mzee->SetParameters(1.00246, 0.00214, 0.00116);
      f1ezee->SetParameters(+3.73e-08, +1.17e-07, +2.02e-07,
                             +1.5e-08, -5.07e-08, -6.39e-08);
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
      // BCDEFGH fit with minitools/drawZmass.C
      //f1mzmm->SetParameters(0.99855, 0., 0.);
      //f1ezmm->SetParameters(pow(0.00010,2),0,0, 0,0,0);
      // Fall17 BCDEF
      //f1mzmm->SetParameters(0.99837, 0., 0.);
      //f1ezmm->SetParameters(pow(0.00019,2),0,0, 0,0,0);
      // Fall17 BCDEF (full Zmm statistics)
      //f1mzmm->SetParameters(0.998440, 0., 0.);
      //f1ezmm->SetParameters(pow(0.00011,2),0,0, 0,0,0);
      // 17Nov17_V10 BCDEF
      f1mzmm->SetParameters(0.99854, 0.00000, 0.00000);
      f1ezmm->SetParameters(+1.53e-08,        +0,        +0,
                                   +0,        +0,        +0);
    }
    else
      hmzmm->Fit(f1mzmm);
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
  //files["multijet"] = fmj;
  files["gamjet"] = fp;
  files["zeejet"] = fzee;
  files["zmmjet"] = fzmm;


  // Map variable names used in different files to common scheme
  map<string, map<string, const char*> > rename;

  // Note: non-CHS results are not available for all samples, so not used
  rename["dijet"]["mpfchs"] = "mpfchs";
  rename["dijet"]["mpfchs1"] = "mpfchs";
  rename["dijet"]["ptchs"] = "ptchs";

  rename["multijet"]["ratio"] = "Data";//"Ratio"; => PATCH
  rename["multijet"]["data"] = "Data";
  rename["multijet"]["mc"] = "MC";
  rename["multijet"]["crecoil"] = "CRecoil";
  rename["multijet"]["mpfchs"] = "MPF";
  rename["multijet"]["mpfchs1"] = "MPF";
  rename["multijet"]["ptchs"] = "MJB";

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
  //sets.push_back("dijet");
  //sets.push_back("multijet");
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
  etas.push_back(make_pair<double,double>(0.522,0.783)); 
  etas.push_back(make_pair<double,double>(0.783,1.044)); 
  etas.push_back(make_pair<double,double>(1.044,1.305)); 
  etas.push_back(make_pair<double,double>(1.305,1.479)); 
  etas.push_back(make_pair<double,double>(1.479,1.653)); 
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
  */
  // Wide eta bins for L2L3Res closure
  etas.push_back(make_pair<double,double>(1.305,1.93));
  etas.push_back(make_pair<double,double>(1.93,2.5));
  etas.push_back(make_pair<double,double>(2.5,2.964));
  etas.push_back(make_pair<double,double>(2.964,3.2));
  etas.push_back(make_pair<double,double>(3.2,5.191));



  vector<double> alphas;
  alphas.push_back(0.10);
  alphas.push_back(0.15);
  alphas.push_back(0.20); //  => patch because of a20 missing for Z+jet
  alphas.push_back(0.30);

  ///////////////////////////////////////////
  // Rename selected graphs and store them //
  ///////////////////////////////////////////

  map<string, map<string, map<string, map<int, map<int, TGraphErrors*> > > > > grs;
  map<string, map<string, map<int, map<int, TH1F*> > > > counts;

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
	    //if (narrowBin) {
	    //  if (s=="dijet")  f = fdj2;  // BCDEFGH only
	    //  if (s=="zeejet") f = fzee2; // GH only
	    //  if (s=="zmmjet") f = fzmm2; // GH only
	    //}
	    assert(f || s=="zlljet");


            if (t=="counts" && s!="zmmjet" && s!="zeejet" && s!="gamjet")
              continue; // counts available only for z+jet and gamjet, so far

	    if (s=="multijet" && (!((epoch!="L4" &&
				     fabs(eta1-0)<0.1 && fabs(eta2-1.3)<0.1) ||
				    (epoch=="L4" &&
				     fabs(eta1-0)<0.1 && fabs(eta2-2.4)<0.1))
				  || fabs(alpha-0.10)<0.01))
	      continue; // only barrel for multijet balance, pT=15,20,30
	    if (s=="gamjet"  && fabs(eta1-3.2)<0.1) { eta1=3.0; eta2=3.2; }

	    // Reconstruct naming scheme used in each of the files
	    // If non-conventional naming schemes, patch here
	    const char *c(0);
	    if (s=="dijet") {
	      c = Form("%s/eta_%02.0f_%02.0f/%s_%s_a%1.0f",
                       dd, 10.001*eta1, 10.001*eta2, rename[s][t], ss, 100.*alpha); 
	    } // dijet
	    if (s=="multijet") {
	      c = Form("%s/Pt%1.0f/%s", rename[s][d], 100.*alpha, rename[s][t]);
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
	      if (doFall17Closure)
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


            if (t=="counts" && (s=="zmmjet" || s=="zeejet" || s=="gamjet") ){ // write out counts to jecdata.root (as TH1D)
              assert(obj->InheritsFrom("TH1D")||obj->InheritsFrom("TH1F"));
	      TH1F *g;
	      if (obj->InheritsFrom("TH1F")) g = (TH1F*)obj;
	      else if (obj->InheritsFrom("TH1D")) g->Copy( *(TH1D*)obj );
	      //              TH1D *g = (TH1D*)obj;
              g->SetName(Form("%s_%s_a%1.0f",tt,ss,100.*alphas[ialpha]));
              g->UseCurrentStyle(); // Basic TDR style
              g->SetMarkerStyle(style[s][t]);
              g->SetMarkerColor(color[s]);
              g->SetLineColor(color[s]);
              g->SetMarkerSize(size[alpha]);
              g->SetDrawOption("SAMEP");
              g->Write();
              counts[d][s][ieta][ialpha] = g;
              continue; // counts available only for z+jet, so far
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
                double ipt = counts[d][s][ieta][ialpha]->FindBin(pt);
		double nentries = counts[d][s][ieta][ialpha]->GetBinContent(ipt);
                if (d=="ratio"){
                  assert(counts["mc"][s][ieta][ialpha]);
                  assert(counts["data"][s][ieta][ialpha]);
                  if(counts["mc"][s][ieta][ialpha]->GetBinContent(ipt) < neventsmin || counts["data"][s][ieta][ialpha]->GetBinContent(ipt) < neventsmin){
                    cout << ipt << " pt " <<pt << " g->GetX()[i] " << g->GetX()[i] << " nentries MC "<< counts["mc"][s][ieta][ialpha]->GetBinContent(ipt) << " nentries data " <<  counts["data"][s][ieta][ialpha]->GetBinContent(ipt)<<  " y: " << g->GetY()[i] << endl;
                    g->RemovePoint(i);
                  }
                }
                else if (nentries < neventsmin){
                  cout << ipt << " pt " <<pt << " g->GetX()[i] " << g->GetX()[i] << " nentries "<< nentries <<  " y: " << g->GetY()[i] << endl;
                  g->RemovePoint(i);
                }
              } // zmm/zee

              //else if (i==0)continue;
              //remove points where difference of central value between points is larger than 5*sigma(higher pt point) to remove spurious high pt gamma+jet points with presumably low stats
              //else if (i>0 && i==g->GetN()-1 && fabs(g->GetY()[i]-g->GetY()[i-1])>g->GetEY()[i]*5.) g->RemovePoint(i);
	    } // for i

	    // // patch MC/data to data/MC for dijet samples
	    // if (s=="dijet" && d=="ratio") {
	    //   for (int i = 0; i != g->GetN(); ++i) {
	    // 	double x = g->GetX()[i];
	    // 	double ex = g->GetEX()[i];
	    // 	double y = g->GetY()[i];
	    // 	double ey = g->GetEY()[i];
	    // 	assert(y!=0);
	    // 	g->SetPoint(i, x, y!=0 ? 1./y : 0);
	    // 	g->SetPointError(i, ex, y!=0 ? ey/(y*y) : 0);
	    //   }
	    // }

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
		assert(!correctGamScale);
		double ptgam = g->GetX()[i];
		double pt = 2*ptgam; // fix post Legacy2016
		//double ipt = hmzee->FindBin(pt);
		double ipt = min(hmzee->FindBin(fzmmptmax), hmzee->FindBin(pt));
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
	    // Photon scale correction (from drawGamVsZmm)
	    // No separate unceratinty added, is already in global fit
	    if (correctGamScale && s=="gamjet" && (d=="data" || d=="ratio")) {
	      assert(!correctGamMass);
	      for (int i = 0; i != g->GetN(); ++i) {
		g->SetPoint(i, g->GetX()[i], g->GetY()[i]*valueGamScale);
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

  //fdj->Close();
  //fp->Close(); // single file
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
    {
      s = Form("%s/Fall17_17Nov2017%s_V28_DATA_L2L3Residual_AK4PFchs.txt",cd,
	       epoch=="BCDEF" ? "B" : epoch.c_str());
      cout << s << endl;
      JetCorrectorParameters *par_l2l3res = new JetCorrectorParameters(s);
      vector<JetCorrectorParameters> vpar;
      vpar.push_back(*par_l2l3res);
      jec = new FactorizedJetCorrector(vpar);

      if (epoch=="BCDEF") {

	// luminosity from Hugues Lattaud, photon+jet 11 June 2018
	double lumtot = 4.8+9.6+4.2+9.3+13.4; // 41.3/fb
	jecwb = 4.8/lumtot;

	s=Form("%s/Fall17_17Nov2017C_V28_DATA_L2L3Residual_AK4PFchs.txt",cd);
	cout << s << endl;
	JetCorrectorParameters *par_c = new JetCorrectorParameters(s);
	vector<JetCorrectorParameters> vpar_c;
	vpar_c.push_back(*par_c);
	jecc = new FactorizedJetCorrector(vpar_c);
	jecwc = 9.6/lumtot;

	s=Form("%s/Fall17_17Nov2017D_V28_DATA_L2L3Residual_AK4PFchs.txt",cd);
	cout << s << endl;
	JetCorrectorParameters *par_d = new JetCorrectorParameters(s);
	vector<JetCorrectorParameters> vpar_d;
	vpar_d.push_back(*par_d);
	jecd = new FactorizedJetCorrector(vpar_d);
	jecwd = 4.2/lumtot;

	s=Form("%s/Fall17_17Nov2017E_V28_DATA_L2L3Residual_AK4PFchs.txt",cd);
	cout << s << endl;
	JetCorrectorParameters *par_e = new JetCorrectorParameters(s);
	vector<JetCorrectorParameters> vpar_e;
	vpar_e.push_back(*par_e);
	jece = new FactorizedJetCorrector(vpar_e);
	jecwe = 9.3/lumtot;

	s=Form("%s/Fall17_17Nov2017F_V10_DATA_L2L3Residual_AK4PFchs.txt",cd);
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
      s = Form("%s/Fall17_17Nov2017%s_V28_DATA_L2L3Residual_AK4PFchs.txt",cd,
	       epoch=="BCDEF"||epoch=="L4" ? "B" : epoch.c_str()); // Fall17
      cout << s << endl;
      JetCorrectorParameters *par_old = new JetCorrectorParameters(s);
      vector<JetCorrectorParameters> v;
      v.push_back(*par_old);
      jecold = new FactorizedJetCorrector(v);
    }

    // Difference between pT-dependent and flat L1
    FactorizedJetCorrector *jecl1flat;
    {
      s = Form("%s/Fall17_17Nov2017%s_V28_DATA_L1RC_AK4PFchs.txt",cd,epoch=="BCDEF"||epoch=="L4" ? "B" : epoch.c_str()); // Fall17
      cout << s << endl << flush;
      JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
      vector<JetCorrectorParameters> v;
      v.push_back(*l1);
      jecl1flat = new FactorizedJetCorrector(v);
    }
    FactorizedJetCorrector *jecl1pt;
    {
      s = Form("%s/Fall17_17Nov2017%s_V28_DATA_L1FastJet_AK4PFchs.txt",cd,epoch=="BCDEF"||epoch=="L4" ? "B" : epoch.c_str()); // Fall17
      cout << s << endl << flush;
      JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
      vector<JetCorrectorParameters> v;
      v.push_back(*l1);
      jecl1pt = new FactorizedJetCorrector(v);
    }

    // Run I uncertainty => 80XV6 uncertainty
    s = Form("%s/Winter14_V8_DATA_UncertaintySources_AK5PFchs.txt",cd); // V8 Run I
    //s2 = "TotalNoFlavorNoTime";
    s2 = "SubTotalAbsolute"; // 07AugV4
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_ref1 = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc_ref1 = new JetCorrectionUncertainty(*p_ref1);

    // Total uncertainty, excluding Flavor and Time
    s =Form("%s/Fall17_17Nov2017B_V28_DATA_UncertaintySources_AK4PFchs.txt",cd);
    s2 = "TotalNoFlavorNoTime";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_unc = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p_unc);

    // Partial uncertainties
    s =Form("%s/Fall17_17Nov2017B_V28_DATA_UncertaintySources_AK4PFchs.txt",cd);
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
	  double val = (epoch=="L4" ? 1 : 1./jec->getCorrection());

	  if (epoch=="BCDEF") {
	    assert(jecc); assert(jecd); assert(jece); assert(jecf);
	    assert(fabs(jecwb+jecwc+jecwd+jecwe+jecwf-1)<1e-4);

	    double valb = val;

	    jecc->setJetEta(eta);
	    jecc->setJetPt(pt);
	    double valc = 1./jecc->getCorrection();

	    jecd->setJetEta(eta);
	    jecd->setJetPt(pt);
	    double vald = 1./jecd->getCorrection();

	    jece->setJetEta(eta);
	    jece->setJetPt(pt);
	    double vale = 1./jece->getCorrection();

	    jecf->setJetEta(eta);
	    jecf->setJetPt(pt);
	    double valf = 1./jecf->getCorrection();

	    val = jecwb*valb +jecwc*valc +jecwd*vald +jecwe*vale +jecwf*valf;
	  }

	  if (doFall17Closure) val = 1;
	  sumval += w*val;
	  sumw += w; // sum weights only once
	  
	  // reference JEC
	  jecrun1->setJetEta(eta);
	  jecrun1->setJetPt(pt);
	  double jesrun1 = (epoch=="L4" ? 1 : 1./jecrun1->getCorrection());
	  sumrun1 += w*jesrun1;

	  // old JEC
	  jecold->setJetEta(eta);
	  jecold->setJetPt(pt);
	  double jes = (epoch=="L4" ? 1 : 1./jecold->getCorrection());
	  sumjes += w*jes;

	  jecl1flat->setJetEta(eta);
	  jecl1flat->setJetPt(pt);
	  jecl1flat->setRho(gRho);
	  jecl1flat->setJetA(TMath::Pi()*0.4*0.4);
	  double vall1flat = 1./jecl1flat->getCorrection();
	  sumvall1flat += w*vall1flat;
	  //
	  jecl1pt->setJetEta(eta);
	  jecl1pt->setJetPt(pt);
	  jecl1pt->setRho(gRho);
	  jecl1pt->setJetA(TMath::Pi()*0.4*0.4);
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
	if (doFall17Closure) jes = 1;
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



