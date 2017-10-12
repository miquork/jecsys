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
const double gRho = 15.36; // 2017-06-02 for 2016 data
const bool _dcsonly = false;

// Minimum pTcut for gamma+jet
double fpmpfptmin(30);   // photon+jet (pTbal issue in 80X-Sum16)
double fpbalptmin(30);   // photon+jet (pTbal issue in 80X-Sum16)
double fzeeptmin(30.); // Zee+jet
double fzmmptmin(30.); // Zmm+jet
// Additional cuts to Z+jet MPF / balance methods
double fzmpfptmin(30.); // Z+jet MPF
double fzbalptmin(30.);

//for fine etabins deactivate ptbal
double fdijetmpfptmin(30);
double fdijetbalptmin(30.);
double fdijetptmax(1500.);

// Maximum pTcut for samples (to avoid bins with too large uncertainty)
double fpptmax(1500.);
double fzeeptmax(1000.);
double fzmmptmax(1000.);

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


  // Fine eta-bins
  // https://indico.cern.ch/event/651009/contributions/2658485/
  // 03Feb2017_L2Res_InputGlobalFit_12July2017.tar

  map<string,const char*> fdj_files;
  fdj_files["BCD"] = "BCD";
  fdj_files["EF"] = "EFearly";
  fdj_files["G"] = "FlateG";
  fdj_files["H"] = "H";
  fdj_files["BCDEFGH"] = "BCDEFGH";
  //  TFile *fdj = new TFile(Form("rootfiles/2017_10_L2ResGlobalFit/JEC_L2_Dijet_%s_2017-09-22_AK4PFchs_pythia8.root",fdj_files[epoch]),"READ");

  //    Subject: Re: BCDEFGH combination of combination files for global fit with 03Feb_03
  //Date: Mon, 9 Oct 2017 15:13:18 +0200
  //From: Multhaup, Jens <jens.multhaup@desy.de>
  TFile *fdj = new TFile(Form("rootfiles/2017_10_L2ResGlobalFit/L2Res_Dijet_09Oct2017_InputGlobalFit/Run%s/output/JEC_L2_Dijet_AK4PFchs_pythia8.root",fdj_files[epoch]),"READ");

  assert(fdj && !fdj->IsZombie());

  // Andrey Popov, Apri 28, 2017 (Feb03_L2ResV2)
  // https://indico.cern.ch/event/634367/
  map<string,const char*> fm_files;
  fm_files["BCD"] = "0428_Run2016BCD"; // also update multijet.C
  fm_files["EF"] = "0428_Run2016EFearly"; // also update multijet.C
  fm_files["G"] = "0428_Run2016FlateG"; // also update multijet.C
  fm_files["H"] = "0428_Run2016H"; // also update multijet.C
  fm_files["BCDEFGH"] = "0428_Run2016All"; // also update multijet.C
  TFile *fmj = new TFile(Form("rootfiles/multijet_2017%s.root",
  			      fm_files[epoch]),"READ");

  assert(fmj && !fmj->IsZombie());
  
  // Hugues Lattaud, May 6, 2017 (03Feb_L2ResV2 + TESTB)
  // https://indico.cern.ch/event/636628/
  map<string,const char*> fp_files;
  fp_files["BCD"] = "runBCD_TESTB";
  fp_files["EF"] = "runEFearly_TESTB";
  fp_files["G"] = "runFlateG_TESTB";
  fp_files["H"] = "runH_TESTB";
  //  fp_files["BCDEFGH"] = "BCDEFG"; // May 19 update (H included?)
  fp_files["BCDEFGH"] = "BCDEFGH"; ///afs/cern.ch/user/h/hlattaud/public/Combinationfile_for_Henning/GammaJet_combinationfile_finebinning_BCDEFGH_03Feb17_V3_noresidual.root
  //  TFile *fp = new TFile(Form("rootfiles/Gjet_%s_03feb2017_V2.root",
  //			     fp_files[epoch]),"READ");
  // /afs/cern.ch/user/h/hlattaud/public/Combinationfile_for_Henning/GammaJet_combinationfile_finebinning_BCDEFGH_03Feb17_V3_noresidual.root
  //Subject: Re: Combination files (no residual corrections applied) Summer16_03Feb2017H_V3; merged BCDEFGH
  //Date: Mon, 9 Oct 2017 14:35:04 +0200
  //From: Hugues LATTAUD <lattaud@ipnl.in2p3.fr>
  // /afs/cern.ch/user/h/hlattaud/public/Combinationfile_for_Henning/Combination_file_gammaplusjet_BCDEFGH_03FeB17_nores_fixed_etabinning.root (fine/wide/all bins)
  TFile *fp = new TFile(Form("rootfiles/2017_10_L2ResGlobalFit/Combination_file_gammaplusjet_%s_03FeB17_nores_fixed_etabinning.root",
  			     fp_files[epoch]),"READ");

  assert(fp && !fp->IsZombie());

  // Subject: 	Re: [urgent] BCDEFGH combination of combination files for global fit with 03Feb_03
  // Date: 	Mon, 25 Sep 2017 15:16:04 +0200
  // JECv7 applied, no residual but MC-based corrections)
  // https://indico.cern.ch/event/665707/
  // Fine eta bins for global fit

  map<string,const char*> fz_files;
  fz_files["BCD"] = "BCD";
  fz_files["EF"] = "EF";
  fz_files["G"] = "G";
  fz_files["H"] = "H";
  fz_files["BCDEFGH"] = "BCDEFGH"; 
  TFile *fzmm = new TFile(Form("rootfiles/2017_10_L2ResGlobalFit/combination_ZJet_Zmm_%s"
                               "_2017-10-10.root",
                               //"_2017-09-25.root",
			       //"_2017-10-04.root",
   			       fz_files[epoch]),"READ");
  TFile *fzee = new TFile(Form("rootfiles/2017_10_L2ResGlobalFit/combination_ZJet_Zee_%s"
                               "_2017-10-10.root",
                               //"_2017-09-25.root",
			       //"_2017-10-04.root",
  			       fz_files[epoch]),"READ");
  //
  assert(fzmm && !fzmm->IsZombie());
  assert(fzee && !fzee->IsZombie());

  // Link to Z mass files (same as above now) and histograms
  // This is used for scaling Z mass back to 1 for Zee and Zmm
  TFile *fmzee = fzee;
  TFile *fmzmm = fzmm;
  assert(fmzee && !fmzee->IsZombie());
  assert(fmzmm && !fmzmm->IsZombie());

  string sr = (epoch=="L4" ? "eta_00_24" : "eta_00_13");
  const char *cr = sr.c_str();
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

  rename["zeejet"]["ratio"] = "Ratio";
  rename["zeejet"]["data"] = "Data";
  rename["zeejet"]["mc"] = "MC";
  rename["zeejet"]["mpfchs"] = "MPF-notypeI_CHS";
  rename["zeejet"]["mpfchs1"] = "MPF_CHS"; 
  rename["zeejet"]["ptchs"] = "PtBal_CHS"; 

  rename["zmmjet"]["ratio"] = "Ratio";
  rename["zmmjet"]["data"] = "Data";
  rename["zmmjet"]["mc"] = "MC";
  rename["zmmjet"]["mpfchs"] = "MPF-notypeI_CHS";
  rename["zmmjet"]["mpfchs1"] = "MPF_CHS"; 
  rename["zmmjet"]["ptchs"] = "PtBal_CHS"; 
  

  // color and style codes
  map<string, map<string, int> > style;

  style["dijet"]["mpfchs1"] = kFullCircle;
  style["dijet"]["ptchs"] = kOpenSquare;
  style["multijet"]["mpfchs1"] = kFullTriangleUp;
  style["multijet"]["ptchs"] = kOpenTriangleUp;
  style["gamjet"]["mpfchs1"] = kFullSquare;
  style["gamjet"]["ptchs"] = kOpenSquare;
  style["zeejet"]["mpfchs1"] = kFullCircle;
  style["zeejet"]["ptchs"] = kOpenCircle;
  style["zmmjet"]["mpfchs1"] = kFullStar;
  style["zmmjet"]["ptchs"] = kOpenStar;

  map<string, int> color;
  color["dijet"] = kBlack;
  color["multijet"] = kBlack;
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

  vector<pair<double,double> > etas;
  // reference region |eta|<1.3
  if (epoch!="L4") etas.push_back(make_pair<double,double>(0,1.305));
  if (epoch=="L4") etas.push_back(make_pair<double,double>(0,2.4));
  // Narrow eta bins for L2Res
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
  // Wide eta bins for L2L3Res closure
//  etas.push_back(make_pair<double,double>(1.305,1.93));
//  etas.push_back(make_pair<double,double>(1.93,2.5));
//  etas.push_back(make_pair<double,double>(2.5,2.964));
//  etas.push_back(make_pair<double,double>(2.964,3.2));
//  etas.push_back(make_pair<double,double>(3.2,5.191));



  vector<double> alphas;
  alphas.push_back(0.10);
  alphas.push_back(0.15);
  alphas.push_back(0.20); //  => patch because of a20 missing for Z+jet
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
	  if (s=="gamjet" && f==0) {
	    if (t=="mpfchs" || t=="mpfchs1") f = fp;//fp1;
	    if (t=="ptchs") f = fp;//fp2;
	  }
	  assert(f);

	  // C_recoil is only meant for multijets
	  if (t=="crecoil" && s!="multijet") continue;

	  for (unsigned int ialpha = 0; ialpha != alphas.size(); ++ialpha) {

	    double alpha = alphas[ialpha];

	    eta1 = etas[ieta].first; // reset to avoid trouble with below
	    eta2 = etas[ieta].second; // reset to avoid trouble with below

	    //continue; // only barrel for multijet balance, pT=10,20,30
	    if (s=="multijet" && (!((epoch!="L4" &&
				     fabs(eta1-0)<0.1 && fabs(eta2-1.3)<0.1) ||
				    (epoch=="L4" &&
				     fabs(eta1-0)<0.1 && fabs(eta2-2.4)<0.1))
				  || fabs(alpha-0.10)<0.01))
	      continue; // only barrel for multijet balance, pT=15,20,30
	    //if (s=="gamjet"  && fabs(eta1-0.0)<0.1) { eta2=1.3; } // ??
	    if (s=="gamjet"  && fabs(eta1-3.2)<0.1) { eta1=3.0; eta2=3.2; }

            //	    // Patch missing Z+jet a20 in Sep25_2017 files
            //            if ((s=="zmmjet" || s=="zeejet") && (alpha==0.20))// => patch because of a20 missing for Z+jet
            //              alpha = 0.15;

	    // Reconstruct naming scheme used in each of the files
	    // If non-conventional naming schemes, patch here
	    const char *c(0);
	    if (s=="dijet") {
	      c = Form("%s/eta_%02.0f_%02.0f/%s_%s_a%1.0f", // Sep2017 eta fine bins
                       dd, 10.001*eta1, 10.001*eta2, rename[s][t], ss, 100.*alpha); 
//              if((eta1==0&&       eta2==1.305)|| //add extra_wide suffix for wide bins
//                 (eta1==0&&       eta2==2.4  )||
//                 (eta1==1.305&&   eta2==1.93 )||
//                 (eta1==1.93&&    eta2==2.5  )||
//                 (eta1==2.5&&     eta2==2.964)||
//                 (eta1==2.964&&   eta2==3.2  )||
//                 (eta1==3.2&&     eta2==5.191)) c = Form("%s_wide/eta_%02.0f_%02.0f/%s_%s_a%1.0f", 
//                                                       dd, 10.001*eta1, 10.001*eta2, rename[s][t], ss, 100.*alpha); 
	    } // dijet
	    if (s=="multijet") {
	      c = Form("%s/Pt%1.0f/%s", rename[s][d], 100.*alpha, rename[s][t]);
	    } // multijet
	    if (s=="gamjet") {
	      c = Form("%s%s_a%1.0f_eta%02.0f_%02.0f",
		       rename[s][t], rename[s][d],
                       //		       100.*alpha, 10.001*eta1, 10.001*eta2);
		       100.*alpha, 10.*eta1, 10.*eta2);
//              if(eta1==0.783 && eta2==1.044)c = Form("%s%s_a%1.0f_eta07_10", //finebinning_BCDEFGH_03Feb17_V3_noresidual.root patch
//                                                     rename[s][t], rename[s][d],
//                                                     100.*alpha);
//              if(eta1==1.930 && eta2==2.172)c = Form("%s%s_a%1.0f_eta19_21", //finebinning_BCDEFGH_03Feb17_V3_noresidual.root patch
//                                                     rename[s][t], rename[s][d],
//                                                     100.*alpha);
//              if(eta1==2.172 && eta2==2.322)c = Form("%s%s_a%1.0f_eta21_23", //finebinning_BCDEFGH_03Feb17_V3_noresidual.root patch
//                                                     rename[s][t], rename[s][d],
//                                                     100.*alpha);
//              if(eta1==2.172 && eta2==2.322)c = Form("%s%s_a%1.0f_eta21_23", //finebinning_BCDEFGH_03Feb17_V3_noresidual.root patch
//                                                     rename[s][t], rename[s][d],
//                                                     100.*alpha);
	    } // gamjet
	    if (s=="zmmjet" || s=="zeejet") {
	      c = Form("%s_%s_a%1.0f_eta_%02.0f_%02.0f_L1L2L3",
	    	       rename[s][d], rename[s][t], 100.*alpha,
	    	       10.01*eta1, 10.01*eta2); // TB // changed to 10.01 to avoid strange rounding bug with eta 2.65 bin edge
	      //	      cout << c << endl << flush;
	    } // Z+jet (TB)
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
	      // Clean out spurious empty pooints
	      if (g->GetY()[i]==0 && g->GetEY()[i]==0) g->RemovePoint(i);
	      // Clean out point outside good ranges
	      //else if (s=="gamjet" && 
	      //       (g->GetX()[i]<fpptmin || g->GetX()[i]>fpptmax))
	      //g->RemovePoint(i);
	      else if (s=="gamjet" && t=="mpfchs1" &&
		       (g->GetX()[i]<fpmpfptmin || g->GetX()[i]>fpptmax))
		g->RemovePoint(i);
	      else if (s=="gamjet" && t=="ptchs" &&
		       (g->GetX()[i]<fpbalptmin || g->GetX()[i]>fpptmax))
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
		       (g->GetX()[i]<fzmpfptmin))
		g->RemovePoint(i);
	      else if ((s=="zmmjet" || s=="zeejet") && t=="ptchs" &&
		       (g->GetX()[i]<fzbalptmin))
		g->RemovePoint(i);
	      // Remove bad point from zmmjet MPF and pT balance at 600
	      //else if (s=="zmmjet" && t=="mpfchs1" && d=="ratio" &&
	      //else if (s=="zmmjet" && d=="ratio" &&
	      //     fabs(g->GetX()[i]-600)<50)
	      //g->RemovePoint(i);
	      //else if (s=="zmmjet" && t=="ptchs" && d=="ratio" &&
	      //       fabs(g->GetX()[i])>500)
	      //g->RemovePoint(i);
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

	    // patch Z+jet pT ratio uncertainty (80X-590/pb)
	    if ((s=="zeejet" || s=="zmmjet") && d=="ratio") { // TB
	      //if (false && (s=="zeejet" || s=="zmmjet") && d=="ratio") { // RS

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

	    // patch Z+jet pT center and uncertainty (76X)
	    /*
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
	    */

	    // Turn on mass corrections on for Zee
	    // (1: add mass, 2: limit to +/-1%, 3: Zmm unc 0.5->0.2%)
	    // Before, NDF=75: 129.8+0.1, 211.2-0.4, 136.5-0.4, 106.8-0.2
	    // After1, NDF=75: 131.2-0.3, 182.9-0.3, 127.1-0.3, 107.9-0.1
	    // After2, NDF=75: 124.5-0.2, 182.2-0.2, 125.2-0.3, 104.9-0.1
	    // After3: NDF=75: 124.6-0.4, 184.1-0.6, 126.5-0.6, 105.3-0.3
	    //
	    // Turn off for 80X BCDEF + GH: only flat -0.15% shift
	    if (false && s=="zeejet" && (d=="data" || d=="ratio")) {
	      for (int i = 0; i != g->GetN(); ++i) {
		double ipt = hmzee->FindBin(g->GetX()[i]);
		double k = max(0.99,min(1.01,hmzee->GetBinContent(ipt)));
		double ek = hmzee->GetBinError(ipt);
		g->SetPoint(i, g->GetX()[i], g->GetY()[i]*k);
 		g->SetPointError(i, g->GetEX()[i], 
				 sqrt(pow(g->GetEY()[i]*k,2) + ek*ek));
	      }
	    }
	    // Turn on mass corrections for Zmm
	    // Before, NDF=75: 129.8-0.2, 211.2+0.6, 136.5+0.5, 106.8+0.4
	    // After1, NDF=75: 131.2+0.2, 182.9+0.6, 127.1+0.5, 107.9+0.3
	    // After2, NDF=75: 124.5+0.2, 182.2+0.6, 125.2+0.5, 104.9+0.3
	    // After3: NDF=75: 124.6+0.1, 184.1+0.4, 126.5+0.3, 105.3+0.2
	    //
	    // Turn off for 80X BCDEF + GH: only flat -0.15% shift
 	    if (false && s=="zmmjet" && (d=="data" || d=="ratio")) {
 	      for (int i = 0; i != g->GetN(); ++i) {
 		double ipt = hmzmm->FindBin(g->GetX()[i]);
		double k = max(0.99,min(1.01,hmzmm->GetBinContent(ipt)));
 		double ek = hmzmm->GetBinError(ipt);
 		g->SetPoint(i, g->GetX()[i], g->GetY()[i]*k);
 		g->SetPointError(i, g->GetEX()[i], 
				 sqrt(pow(g->GetEY()[i]*k,2) + ek*ek));
 	      }
 	    }

	    /*
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
	    */

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
    {
      //s = Form("%s/Fall15_25nsV1M2_DATA_L2L3Residual_AK4PFchs.txt",cd); // 76X V2
      //s = Form("%s/Spring16_25nsV3M1_DATA_L2L3Residual_AK4PFchs.txt",cd); // 80X V3M1
      //s = Form("%s/Spring16_25nsV3M1std_DATA_L2L3Residual_AK4PFchs.txt",cd); // 80X V3M1std
      //s = Form("%s/Spring16_25nsV3_DATA_L2L3Residual_AK4PFchs.txt",cd); // 80X V3
      //s = Form("%s/Spring16_25nsV4M1_DATA_L2L3Residual_AK4PFchs.txt",cd); // 80X V4
      //s = Form("%s/Spring16_25nsV4M1_DATA_L2L3Residual_AK4PFchs.txt",cd); // 80X V4
      //s = Form("%s/Spring16_25nsV7M1%s_DATA_L2L3Residual_AK4PFchs.txt",cd,ce); // 80X V7
      //s = Form("%s/Spring16_25nsV7M2%s_DATA_L2L3Residual_AK4PFchs.txt",cd,ce); // 80X V7G special
      //s = Form("%s/Spring16_25nsV8%s_DATA_L2L3Residual_AK4PFchs.txt",cd,epoch=="G"||epoch=="H"||epoch=="GH" ? "p2" : ce); // 80XV8
      //s = Form("%s/Spring16_23Sep2016%sV1_DATA_L2L3Residual_AK4PFchs.txt",cd,epoch=="G"||epoch=="H"||epoch=="GH" ? "GH" : ce); // 80XreV1
      //s = Form("%s/Spring16_23Sep2016%sV1_DATA_L2L3Residual_AK4PFchs.txt",cd,epoch=="G"||epoch=="H"||epoch=="GH"||epoch=="BCDEFGH" ? "GH" : (epoch=="EF" ? "E" : (epoch=="BCDEF" ? "BCD" : ce))); // 80XreV1+Sum16 pre
      //s = Form("%s/Summer16_23Sep2016%sV1_DATA_L2Residual_AK4PFchs.txt",cd,epoch=="GH" ? "G" : ce); // 80XreV1+Sum16
      //s = Form("%s/Summer16_23Sep2016%sV2_DATA_L2L3Residual_AK4PFchs.txt",cd,epoch=="GH"||epoch=="L4" ? "G" : (epoch=="E"||epoch=="F" ? "EF" : (epoch=="BCDEF" ? "BCD" : ce))); // Sum16
      //s = Form("%s/Summer16_23Sep2016%sV3_DATA_L2L3Residual_AK4PFchs.txt",cd,epoch=="GH"||epoch=="L4" ? "G" : (epoch=="E"||epoch=="F" ? "EF" : (epoch=="BCDEF" ? "BCD" : ce))); // Sum16
      //s = Form("%s/Summer16_03Feb2017%s_V3_DATA_L2L3Residual_AK4PFchs.txt",cd,epoch=="GH"||epoch=="L4" ? "G" : (epoch=="E"||epoch=="F" ? "EF" : (epoch=="BCDEF" ||epoch=="BCDEFGH" ? "BCDEFGH" : ce))); // 03Feb
      //s = Form("%s/Summer16_03Feb2017_10Oct2017_MPF_Hybrid_%s_L2L3Residual_pythia8_AK4PFchs.txt",cd,epoch=="GH"||epoch=="L4" ? "G" : (epoch=="E"||epoch=="F" ? "EF" : (epoch=="BCDEF" ||epoch=="BCDEFGH" ? "BCDEFGH" : ce))); // 03Feb 10Oct2017 with L3 from V3
      //s = Form("%s/Summer16_03Feb2017_10Oct2017_Val_pT_Hybrid_%s_L2L3Residual_pythia8_AK4PFchs.txt",cd,epoch=="GH"||epoch=="L4" ? "G" : (epoch=="E"||epoch=="F" ? "EF" : (epoch=="BCDEF" ||epoch=="BCDEFGH" ? "BCDEFGH" : ce))); // 03Feb 10Oct2017 with L3 from V3
      s = Form("%s/Summer16_03Feb2017_10Oct2017_Val_MPF_LOGLIN_%s_L2L3Residual_pythia8_AK4PFchs.txt",cd,epoch=="GH"||epoch=="L4" ? "G" : (epoch=="E"||epoch=="F" ? "EF" : (epoch=="BCDEF" ||epoch=="BCDEFGH" ? "BCDEFGH" : ce))); // 03Feb 10Oct2017 with L3 from V3
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
      //s = Form("%s/Spring16_25nsV4M1_DATA_L2L3Residual_AK4PFchs.txt",cd); // 80X V4
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
      //s = Form("%s/Winter14_V6_DATA_L2L3Residual_AK5PFchs.txt",cd); // 74X V7
      //s = Form("%s/Fall15_25ns_COMB_LOGLIN_L2Residual_v2_AK4PFchs_nokFSR.txt",cd); // 76X V2
      //s = Form("%s/Spring16_25ns_MPF_LOGLIN_L2Residual_pythia8_v3_AK4PFchs.txt",cd); // 80X V3M1
      //s = Form("%s/Spring16_25nsV7%s_DATA_L2Residual_AK4PFchs.txt",cd,ce); // 80X V7
      //s = Form("%s/Spring16_25nsV8%s_DATA_L2Residual_AK4PFchs.txt",cd,epoch=="G"||epoch=="H"||epoch=="GH" ? "p2" : ce); // 80XV8
      //s = Form("%s/Spring16_23Sep2016%sV1_DATA_L2Residual_AK4PFchs.txt",cd,epoch=="G"||epoch=="H"||epoch=="GH" ? "GH" : ce); // 80XreV1
      //s = Form("%s/Summer16_23Sep2016%sV1_DATA_L2Residual_AK4PFchs.txt",cd,epoch=="GH"||epoch=="BCDEFGH"||epoch=="L4" ? "G" : (epoch=="BCDEF" ? "BCD" : ce)); // 80XreV1+Sum16
      //s = Form("%s/Summer16_03Feb2017%s_V2_DATA_L2Residual_AK4PFchs.txt",cd,epoch=="GH"||epoch=="BCDEFGH"||epoch=="L4" ? "G" : (epoch=="BCDEF" ? "BCD" : ce));
      //s = Form("%s/Summer16_03Feb2017%s_V3_DATA_L2Residual_AK4PFchs.txt",cd,epoch=="GH"||epoch=="L4" ? "G" : (epoch=="BCDEF" ? "BCD" : ce)); // 03Feb
      s = Form("%s/Summer16_03Feb2017_10Oct2017_Val_MPF_LOGLIN_%s_L2Residual_pythia8_AK4PFchs.txt",cd,epoch=="GH"||epoch=="L4" ? "G" : (epoch=="E"||epoch=="F" ? "EF" : (epoch=="BCDEF" ||epoch=="BCDEFGH" ? "BCDEFGH" : ce))); // 03Feb 10Oct2017
      cout << s << endl;
      JetCorrectorParameters *par_old = new JetCorrectorParameters(s);
      vector<JetCorrectorParameters> v;
      v.push_back(*par_old);
      jecold = new FactorizedJetCorrector(v);
    }

    // Difference between pT-dependent and flat L1
    FactorizedJetCorrector *jecl1flat;
    {
      //s = Form("%s/Fall15_25nsV1_DATA_L1RC_AK4PFchs.txt",cd); // 76X V2
      //s = Form("%s/Spring16_25nsV2_DATA_L1RC_AK4PFchs.txt",cd); // 80X V3M1
      //s = Form("%s/Spring16_25nsV3_DATA_L1RC_AK4PFchs.txt",cd); // 80X V3
      //s = Form("%s/Spring16_25nsV7%s_DATA_L1RC_AK4PFchs.txt",cd,ce); // 80X VV7
      //s = Form("%s/Spring16_25nsV8%s_DATA_L1RC_AK4PFchs.txt",cd,epoch=="G"||epoch=="H"||epoch=="GH" ? "p2" : ce); // 80XV8
      //s = Form("%s/Spring16_23Sep2016%sV1_DATA_L1RC_AK4PFchs.txt",cd,epoch=="G"||epoch=="H"||epoch=="GH" ? "GH" : ce); // 80XreV1
      //s = Form("%s/Summer16_23Sep2016%sV1_DATA_L1RC_AK4PFchs.txt",cd,epoch=="GH"||epoch=="BCDEFGH"||epoch=="L4" ? "G" : (epoch=="BCDEF" ? "BCD" : ce)); // 80XreV1+Sum16
      s = Form("%s/Summer16_03Feb2017%s_V3_DATA_L1RC_AK4PFchs.txt",cd,epoch=="GH"||epoch=="BCDEFGH"||epoch=="L4" ? "G" : (epoch=="BCDEF" ? "BCD" : ce));
      cout << s << endl << flush;
      JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
      vector<JetCorrectorParameters> v;
      v.push_back(*l1);
      jecl1flat = new FactorizedJetCorrector(v);
    }
    FactorizedJetCorrector *jecl1pt;
    {
      //s = Form("%s/Fall15_25nsV1_DATA_L1FastJet_AK4PFchs.txt",cd); // 76X V2
      //s = Form("%s/Spring16_25nsV2_DATA_L1FastJet_AK4PFchs.txt",cd); // 80X V3M1
      //s = Form("%s/Spring16_25nsV3_DATA_L1FastJet_AK4PFchs.txt",cd); // 80X V3
      //s = Form("%s/Spring16_25nsV7%s_DATA_L1FastJet_AK4PFchs.txt",cd,ce); // 80X V7
      //s = Form("%s/Spring16_25nsV8%s_DATA_L1FastJet_AK4PFchs.txt",cd,epoch=="G"||epoch=="H"||epoch=="GH" ? "p2" : ce); // 80XV8
      //s = Form("%s/Spring16_23Sep2016%sV1_DATA_L1FastJet_AK4PFchs.txt",cd,epoch=="G"||epoch=="H"||epoch=="GH" ? "GH" : ce); // 80XreV1
      //s = Form("%s/Summer16_23Sep2016%sV1_DATA_L1FastJet_AK4PFchs.txt",cd,epoch=="GH"||epoch=="BCDEFGH"||epoch=="L4" ? "G" : (epoch=="BCDEF" ? "BCD" : ce)); // 80XreV1+Sum16
      s = Form("%s/Summer16_03Feb2017%s_V3_DATA_L1FastJet_AK4PFchs.txt",cd,epoch=="GH"||epoch=="BCDEFGH"||epoch=="L4" ? "G" : (epoch=="BCDEF" ? "BCD" : ce)); 
      cout << s << endl << flush;
      JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
      vector<JetCorrectorParameters> v;
      v.push_back(*l1);
      jecl1pt = new FactorizedJetCorrector(v);
    }

    // Run I uncertainty => 80XV6 uncertainty
    //s = Form("%s/Winter14_V10M_DATA_UncertaintySources_AK5PFchs.txt",cd); // V8 Run I
    s = Form("%s/Winter14_V8_DATA_UncertaintySources_AK5PFchs.txt",cd); // V8 Run I (official file, but same as above)
    //s = Form("%s/Spring16_25nsV6_DATA_UncertaintySources_AK4PFchs.txt",cd); // 80X V6
    //s = Form("%s/Spring16_25nsV8M1_DATA_UncertaintySources_AK4PFchs.txt",cd); // current (p2 bogged, use M1)
    s2 = "TotalNoFlavorNoTime";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_ref1 = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc_ref1 = new JetCorrectionUncertainty(*p_ref1);

    // Total uncertainty, excluding Flavor and Time
    //s = Form("%s/Fall15_25nsV1M2_DATA_UncertaintySources_AK4PFchs.txt",cd); // 76X V2
    //s = Form("%s/Spring16_25nsV4M3_DATA_UncertaintySources_AK4PFchs.txt",cd); // 80X V4
    //s = Form("%s/Winter14_V10M_DATA_UncertaintySources_AK5PFchs.txt",cd); // 80X V7 special
    //s = Form("%s/Spring16_25nsV8%s_DATA_UncertaintySources_AK4PFchs.txt",cd,epoch=="G" ? "M1" : ce); // current (p2 bogged, use M1)
    //s = Form("%s/Spring16_25nsV8M1_DATA_UncertaintySources_AK4PFchs.txt",cd); // current (M1=G, should add Time for others)
    //s = Form("%s/Spring16_25nsV10p2_DATA_UncertaintySources_AK4PFchs.txt",cd);
    //s = Form("%s/Spring16_23Sep2016GHV1_DATA_UncertaintySources_AK4PFchs.txt",cd);
    //s = Form("%s/Summer16_23Sep2016V3_DATA_UncertaintySources_AK4PFchs.txt",cd); // Sum16
    s = Form("%s/Summer16_03Feb2017_V3_DATA_UncertaintySources_AK4PFchs.txt",cd); // Sum16
    s2 = "TotalNoFlavorNoTime";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_unc = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p_unc);

    // Partial uncertainties
    //s = Form("%s/Fall15_25nsV1M2_DATA_UncertaintySources_AK4PFchs.txt",cd); // 76X V2
    //s = Form("%s/Spring16_25nsV4M3_DATA_UncertaintySources_AK4PFchs.txt",cd); // 80X V3
    //s = Form("%s/Winter14_V10M_DATA_UncertaintySources_AK5PFchs.txt",cd); // 80X V7 special
    //s = Form("%s/Spring16_25nsV8%s_DATA_UncertaintySources_AK4PFchs.txt",cd,epoch=="G" ? "M1" : ce); // current (p2 bogged, use M1)
    //s = Form("%s/Spring16_25nsV8M1_DATA_UncertaintySources_AK4PFchs.txt",cd); // current (M1=G, should add Time for others)
    //s = Form("%s/Spring16_25nsV10p2_DATA_UncertaintySources_AK4PFchs.txt",cd);
    //s = Form("%s/Spring16_23Sep2016GHV1_DATA_UncertaintySources_AK4PFchs.txt",cd);
    //s = Form("%s/Summer16_23Sep2016V3_DATA_UncertaintySources_AK4PFchs.txt",cd);
    s = Form("%s/Summer16_03Feb2017_V3_DATA_UncertaintySources_AK4PFchs.txt",cd);
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
	  double val = (epoch=="L4" ? 1 : 1./jec->getCorrection());
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
	  //double jes = 1.; // 74X V7
	  double jes = (epoch=="L4" ? 1 : 1./jecold->getCorrection()); // 76X V2
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



