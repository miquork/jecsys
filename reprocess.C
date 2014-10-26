// File: reprocess.C
// Created by Mikko Voutilainen, on Sep 6th, 2012
// Updated on Oct 25, 2014 (cleanup for 8 TeV JEC paper)
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

// put all the different methods in a single file for easy access for everybody
void reprocess() {

  // Set TDR style to have correct graphical setttings when storing graphs
  setTDRStyle();

  TDirectory *curdir = gDirectory;

  TFile *fout = new TFile("rootfiles/jecdata.root","RECREATE");

  ////////////////////////
  // Define input files //
  ////////////////////////

  // On 30 Jul 2014, at 14:51, from Denis Rathjens
  // Re: Summer14 combination files with RD MC
  TFile *fdj = new TFile("rootfiles/Winter13V4_DJ.root","READ");
  assert(fdj && !fdj->IsZombie());
  
  // On 23 Jun 2014, at 22:48, from Viola Sordini
  // Re: typeI fix - errata ?
  TFile *fp = new TFile("rootfiles/gamma_jet_response_all_alphas_with_raw_footprint_and_typeIfix_NEW.root","READ");

  // On 30 Jul 2014, at 11:47, from Rajdeep Mohan Chatterjee
  // Re: L1Fix for type-I MET
  TFile *fzee = new TFile("rootfiles/Zee_JEC_Summary_Summer14_RDMC_L1Fix.root",
			  "READ");
  assert(fzee && !fzee->IsZombie());

  // On 20 Jun 2014, at 12:29, from Dominik Haitz
  // Re: type-I fix
  TFile *fzmm = new TFile("rootfiles/2014-06-20_Zmm_typeI-fixed.root","READ");
  assert(fzmm && !fzmm->IsZombie());

  // On 28 Feb 2014, at 14:00, from Dominik Haitz
  // Re: pT~200 GeV
  TFile *feta = new TFile("rootfiles/th2d_jeteta_zpt.root","READ");
  assert(feta && !feta->IsZombie());

  // This is used for correctly averaging JEC and its uncertainty
  // for the wide eta bins used in global fit combinations
  TH2D *h2eta = (TH2D*)feta->Get("data"); assert(h2eta);

  map<string, TFile*> files;
  files["dijet"] = fdj;
  files["gamjet"] = fp;
  files["zeejet"] = fzee;
  files["zmmjet"] = fzmm;


  // Map variable names used in different files to common scheme
  map<string, map<string, const char*> > rename;

  // Note: non-CHS results are not available for all samples, so not used
  rename["dijet"]["ratio"] = "ratio"; 
  rename["dijet"]["data"] = "data"; 
  rename["dijet"]["mc"] = "MC"; 
  rename["dijet"]["mpfchs"] = "MPFchs_VsPtAve"; 
  rename["dijet"]["mpfchs1"] = "MPFT1chs_VsPtAve";
  rename["dijet"]["ptchs"] = "PtBalchs";

  rename["gamjet"]["ratio"] = "";
  rename["gamjet"]["data"] = "_DATA"; 
  rename["gamjet"]["mc"] = "_MC"; 
  rename["gamjet"]["mpfchs"] = "MPFchs";
  rename["gamjet"]["mpfchs1"] = "MPFchs"; 
  rename["gamjet"]["ptchs"] = "PtBalchs"; 

  rename["zeejet"]["ratio"] = "DataByMC"; 
  rename["zeejet"]["data"] = "Data"; 
  rename["zeejet"]["mc"] = "MC"; 
  rename["zeejet"]["mpfchs"] = "MPFchs";
  rename["zeejet"]["mpfchs1"] = "MPFchs_Fix";
  rename["zeejet"]["ptchs"] = "PtBalchs"; 

  rename["zmmjet"]["ratio"] = "ratio"; 
  rename["zmmjet"]["data"] = "data";
  rename["zmmjet"]["mc"] = "MC";
  rename["zmmjet"]["mpfchs"] = "MPF-notypeI_CHS";
  rename["zmmjet"]["mpfchs1"] = "MPF_CHS"; 
  rename["zmmjet"]["ptchs"] = "PtBal_CHS"; 
  

  // color and style codes
  map<string, int> style;
  style["mpfchs"] = kFullStar;
  style["mpf"] = kFullDiamond;
  style["mpfchs1"] = kFullCircle;
  style["mpf1"] = kFullSquare;
  style["ptchs"] = kOpenCircle;
  style["pt"] = kOpenSquare;

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
  dirs.push_back("ratio");
  dirs.push_back("data");
  dirs.push_back("mc");

  vector<string> types;
  types.push_back("mpfchs");
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

  // Loop over data, MC, ratio
  for (unsigned int idir = 0; idir != dirs.size(); ++idir) {

    string d = dirs[idir];
    const char *dd = d.c_str();

    fout->mkdir(dd);
    assert(fout->cd(dd));
    TDirectory *dout0 = gDirectory;

    // Loop over eta bins
    for (unsigned int ieta = 0; ieta != etas.size(); ++ieta) {
      
      double eta1 = etas[ieta].first;
      double eta2 = etas[ieta].second;
      dout0->mkdir(Form("eta%02.0f-%02.0f",eta1*10.,eta2*10.));
      assert(dout0->cd(Form("eta%02.0f-%02.0f",eta1*10.,eta2*10.)));
      TDirectory *dout = gDirectory;

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

	  // Take pT and MPF from different files for gamma+jet
	  //if (s=="gamjet" && f==0) {
	  //if (t=="mpfchs" || t=="mpfchs1") f = fp1;
	  //if (t=="ptchs") f = fp2;
	  //}
	  //assert(f);

	  for (unsigned int ialpha = 0; ialpha != alphas.size(); ++ialpha) {

	    double alpha = alphas[ialpha];

	    eta1 = etas[ieta].first; // reset to avoid trouble with below
	    eta2 = etas[ieta].second; // reset to avoid trouble with below

	    // If samples missing break-up into e.g. [3,3.2] and [3.2,5.2] bins
	    // or merged [0,1.3] bin, patch here
	    //if (s=="zeejet" && eta1==0 && fabs(eta2-1.3)<0.1) { eta2=1.305; }

	    // If some subsets of data missing, skip here
	    // gamjet missing non-CHS graphs
	    //if (s=="gamjet" && (t=="mpf" || t=="mpf1" || t=="pt"))
	    //continue;
	    if (s=="dijet" && (fabs(eta1-0.783)<0.1 || fabs(eta2-0.783)<0.1))
		continue;

	    // Reconstruct naming scheme used in each of the files
	    // If non-conventional naming schemes, patch here
	    const char *c(0);
	    if (s=="dijet") {
	      c = Form("%s_a%1.0f_eta%1.0f_%1.0f_%s",
		       rename[s][t], 100.*alpha, 10.*eta1, 10.*eta2,
		       rename[s][d]);
	    } // dijet
	    if (s=="gamjet") {
	      if (t=="mpfchs") // special *_raw_* naming for non-type-I MPF
		c = Form("%s%s_raw_a%1.0f_eta%02.0f_%02.0f",
			 rename[s][t], rename[s][d],
			 100.*alpha, 10.*eta1, 10.*eta2);
	      else
		c = Form("%s%s_a%1.0f_eta%02.0f_%02.0f",
			 rename[s][t], rename[s][d],
			 100.*alpha, 10.*eta1, 10.*eta2);
	    } // gamjet
	    if (s=="zeejet") {
	      c = Form("%s_%s_a%1.0f_eta%1.0f_%1.0f",
		       rename[s][d], rename[s][t],
		       100.*alpha, 1000.*eta1, 1000.*eta2);
	    } // zeejet
	    if (s=="zmmjet") {
	      c = Form("%s_%s_a%1.0f_eta%02.0f_%02.0f_L1L2L3Res",
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
	      else if (g->GetX()[i]>4000 || g->GetX()[i]<10) g->RemovePoint(i);
	    }

	    dout->cd();

	    // Set uniform naming scheme and graphical style
	    g->SetName(Form("%s_%s_a%1.0f",tt,ss,100.*alpha));
	    g->UseCurrentStyle(); // Basic TDR style
	    g->SetMarkerStyle(style[t]);
	    g->SetMarkerColor(color[s]);
	    g->SetLineColor(color[s]);
	    g->SetMarkerSize(size[alpha]);
	    g->SetDrawOption("SAMEP");

	    g->Write();

	  } // for ialpha
	} // for iset
      } // for itype
    } // for ieta
  } // for itier

  fdj->Close();
  fp->Close();
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
      s = Form("%s/Winter14_V5_DATA_L2L3Residual_AK5PFchs.txt",cd);
      cout << s << endl;
      JetCorrectorParameters *par_l2l3res = new JetCorrectorParameters(s);
      vector<JetCorrectorParameters> vpar;
      vpar.push_back(*par_l2l3res);
      jec = new FactorizedJetCorrector(vpar);
    }

    // Store old JEC for undoing it in global fit
    FactorizedJetCorrector *jecold;
    {
      s = Form("%s/Winter14_V1_DATA_L2L3Residual_AK5PFchs.txt",cd);
      cout << s << endl;
      JetCorrectorParameters *par_old = new JetCorrectorParameters(s);
      vector<JetCorrectorParameters> v;
      v.push_back(*par_old);
      jecold = new FactorizedJetCorrector(v);
    }

    // Difference between pT-dependent and flat L1
    FactorizedJetCorrector *jecl1flat;
    {
      s = Form("%s/Winter14_V0_DATA_L1FastJetPU_AK5PFchs_pt.txt",cd);
      cout << s << endl << flush;
      JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
      vector<JetCorrectorParameters> v;
      v.push_back(*l1);
      jecl1flat = new FactorizedJetCorrector(v);
    }
    FactorizedJetCorrector *jecl1pt;
    {
      s = Form("%s/Winter14_V1_DATA_L1FastJet_AK5PFchs.txt",cd);
      cout << s << endl << flush;
      JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
      vector<JetCorrectorParameters> v;
      v.push_back(*l1);
      jecl1pt = new FactorizedJetCorrector(v);
    }

    // Total uncertainty, excluding Flavor and Time
    s = Form("%s/Winter14_V5_DATA_Uncertainty_AK5PFchs.txt",cd);
    cout << s << endl << flush;
    JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(s);

    // Partial uncertainties
    s = Form("%s/Winter14_V5_DATA_UncertaintySources_AK5PFchs.txt",cd);
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
      double eta1 = etas[ieta].first; double eta2 = etas[ieta].second;
      assert(gDirectory->cd(Form("eta%02.0f-%02.0f",eta1*10.,eta2*10.)));
      cout << "Eta bin:" << eta1 <<"-"<< eta2 << endl;

      
      const double ptbins[] = {29,30, 40, 50, 60, 75, 100, 125, 155, 180,
			       210, 250, 300, 350, 400, 500, 600, 800, 1000,
			       1200, 1500, 1600,1601};
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
	double sumjes(0), sumvall1flat(0), sumvall1pt(0);
	double sumerr2_pt(0), sumerr2_hcal(0), sumerr2_ecal(0), sumerr2_pu(0);
	double sumerr2_noflv(0), sumerr2_ref(0);
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

	  jecold->setJetEta(eta);
	  jecold->setJetPt(pt);
	  double jes = 1./jecold->getCorrection();
	  sumjes += w*jes;

	  jecl1flat->setJetEta(eta);
	  jecl1flat->setJetPt(pt);
	  //jecl1flat->setNPV(int(14.4));
	  jecl1flat->setRho(12.0);
	  jecl1flat->setJetA(TMath::Pi()*0.5*0.5);
	  double vall1flat = 1./jecl1flat->getCorrection();
	  sumvall1flat += w*vall1flat;
	  //
	  jecl1pt->setJetEta(eta);
	  jecl1pt->setJetPt(pt);
	  //jecl1pt->setNPV(int(14.4));
	  jecl1pt->setRho(12.0);
	  jecl1pt->setJetA(TMath::Pi()*0.5*0.5);
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

	double jes = (sumjes / sumw);
	hjes->SetBinContent(ipt, jes);
	double l1pt = (sumvall1pt / sumw);
	double l1flat = (sumvall1flat / sumw);
	double l1bias = l1pt / l1flat;
	hl1bias->SetBinContent(ipt, l1bias);
	hl1bias->SetBinError(ipt, 0.5*fabs(1-l1bias));
      } // ipt

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

      hjes->Write();
      hl1bias->Write();
    } // ieta
  } // JEC+sys

  fout->Close();

} // reprocess
