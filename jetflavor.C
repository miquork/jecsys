// File: flavor.C
// Created by Mikko Voutilainen, on May 20th, 2018
// Purpose: add jet flavor fractions into the global fit combination file
//          These are MC-only and independent of IOV, so ok to add once for all
//          Especially in the current initial stages some hacking is needed
//          so better not clutter the reprocess.C
//          The fractions go into the 'flavor' subdirectory in each folder

// Check results vs Fig.28 in https://arxiv.org/pdf/1607.03663.pdf
// Should have fractions shifted by ratio of sqrt(s), but otherwise match 'Ph'
// May need to check alpha dependence, if alpha<0.2 vs alpha<0.3 matters

// For dijet, consider adding separate fractions for tag jet,
// then study if tag jet fractions vary with probe jet eta
// (need to extract tag from dijet_XX samples by summing)

#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
//#include "TProfile.h"
//#include "TCanvas.h"
//#include "TLegend.h"
#include "TF1.h"
//#include "TMultiGraph.h"
//#include "TLatex.h"
//#include "TMinuit.h"
//#include "TMatrixD.h"
//#include "TVectorD.h"
//#include "TMath.h"

//#include "../tools.h"
#include "tdrstyle_mod14.C"

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>

using namespace std;

// Loop over all eta bins
//void jetflavorBin(double etamin=0.0, double etamax=1.3, string epoch="BCDEFGH",
void jetflavorBin(double etamin, double etamax, string epoch,
		  map<string, map<string, vector<double> > > &pars);

void jetflavor() {
  
  cout << "Starting jetflavor()" << endl << flush;

  // wider bins
  const int npar = 4;
  const int neta = 6;
  double etabins[neta+1] = {0, 13, 19, 25, 30, 32, 52};
  
  map<int, map<string, map<string, vector<double> > > > pars;
  for (int i = 0; i != neta; ++i) {
    
    //pars[int(etabins[i]+0.5)] = 
      jetflavorBin(0.1*etabins[i],0.1*etabins[i+1],
		   "BCDEFGH", pars[i]);
  }

  cout << "Ready to write results out" << endl << flush;
  
  // Store results in a text file
  //ofstream fout("CondFormats/JetMETObjects/data/Summer16_07Aug2017_V10_Flavor_Pythia8_MC_Fractions_AK4PFchs.txt");
  // NB: name must contain "L5Flavor" for the JEC code to work
  //     (or possibly tag would need to be different?)
  ofstream fout("txt/Summer16_07Aug2017_V10_Flavor_Fractions_MC_L5Flavor_AK4PFchs.txt",ios::out);
  
  const int nsamples = 3;
  string samples[nsamples] = {"zlljet","gamjet","dijet"};
  const int nflavors = 6;
  string flavors[nflavors] = {"ud","s","c","b","g","oth"};

  map<string, double> ptmin;
  ptmin["dijet"]  = 20;
  ptmin["zlljet"] = 30;
  ptmin["gamjet"] = 40;
  map<string, double> emax;
  emax["dijet"]  = 2000;
  emax["zlljet"] = 2000;
  emax["gamjet"] = 3000;

  for (int isample = 0; isample != nsamples; ++isample) {
    
    const char *cs = samples[isample].c_str();
    
    for (int iflv = 0; iflv != nflavors; ++iflv) {
      
      const char *cf = flavors[iflv].c_str();
      cout << "Sample " << cs << " flavor " << cf << endl << flush;

      fout << Form("[%s_%s]",cs,cf) << endl;
      fout << "{1 JetEta 1 JetPt "
	"[0]+[1]*log10(0.01*x)+[2]*pow(log10(0.01*x),2)+"
	"[3]*pow(log10(0.01*x),3) Correction L5Flavor}" << endl;
      
      // Negative side eta side
      for (int ieta = neta-1; ieta != -1; --ieta) {
	
	double eta1 = 0.1*etabins[ieta];
	double eta2 = 0.1*etabins[ieta+1];
	double xmin = ptmin[cs];
	double xmax = emax[cs]/cosh(eta1);
	fout << Form("  %4.1f  %4.1f  6 %2d %4d  ",
		     -eta2, -eta1, int(xmin+0.5), int(xmax+0.5));
		     //-eta2, -eta1, 10, int(3000./cosh(eta1)));
	//vector<double> &p = pars[int(eta1+0.5)][cs][cf];
	//vector<double> &p = pars[ieta][cs][cf];
	double *p = &pars[ieta][cs][cf][0];
	//for (unsigned int ipar = 0; ipar != p.size(); ++ipar) {
	for (unsigned int ipar = 0; ipar != npar; ++ipar) {
	  fout << Form("%11.5g ", p[ipar]);
	} // for ipar
	fout << endl;
      } // for ieta

      // Positive eta side
      for (int ieta = 0; ieta != neta; ++ieta) {
	
	double eta1 = 0.1*etabins[ieta];
	double eta2 = 0.1*etabins[ieta+1];
	double xmin = ptmin[cs];
	double xmax = emax[cs]/cosh(eta1);
	fout << Form("  %4.1f  %4.1f  6 %2d %4d  ",
		     eta1, eta2, int(xmin+0.5), int(xmax+0.5));
	//eta1, eta2, 10, int(3000./cosh(eta1)));
	//vector<double> &p = pars[int(eta1+0.5)][cs][cf];
	//vector<double> &p = pars[ieta][cs][cf];
	double *p = &pars[ieta][cs][cf][0];
	//for (unsigned int ipar = 0; ipar != p.size(); ++ipar) {
	for (unsigned int ipar = 0; ipar != npar; ++ipar) {
	  fout << Form("%11.5g ", p[ipar]);
	} // for ipar
	fout << endl;
      } // for ieta
      //fout << endl;
    } // for iflv
  } // for isample
  
} // jetflavors

//void jetflavor(double etamin=0.0, double etamax=1.3, string epoch="BCDEFGH") {

void jetflavorBin(double etamin, double etamax, string epoch,
		  map<string, map<string, vector<double> > > &pars) {

  // Flavor fit function and storage for parameters
  TF1 *f1 = new TF1(Form("f1_%1.2f-%1.1f",etamin,etamax),
		    "[0] + [1]*log10(0.01*x) + [2]*pow(log10(0.01*x),2) +"
		    "[3]*pow(log10(0.01*x),3)",10,6500);
  //map<string, map<string, vector<double> > > pars;


  cout << "Calling jetflavorBin("<<etamin<<","<<etamax<<");\n"<<flush;

  setTDRStyle();
  //writeExtraText = false; // for JEC paper CWR

  TDirectory *curdir = gDirectory;

  // Open jecdata.root produced by reprocess.C
  TFile *finout = new TFile(Form("rootfiles/jecdata%s.root",epoch.c_str()),
			    "UPDATE");
  assert(finout && !finout->IsZombie());
  
  TFile *fp = new TFile("rootfiles/Gjet_combinationfile_07Aug17_L2res_V6_BCDEFGH_EGM3.root","READ");
  assert(fp && !fp->IsZombie());

  TFile *fd = new TFile("rootfiles/JEC_L2_Dijet_AK4PFchs_pythia8_07122017hcalCleaning_wideEtabins.root","READ");
  assert(fd && !fd->IsZombie());

  // https://indico.cern.ch/event/715277/#7-residulas-with-zjets-status
  TFile *fz = new TFile("rootfiles/flavor_fractions_Fall17_JECV4_Zmm_DYNJ_Madgraph_2018-03-14.root","READ");
  assert(fz && !fz->IsZombie());

  cout << "Files open" << endl << flush;

  const int ndirs = 1;//3;
  const char* dirs[ndirs] = {"ratio"};//{"data", "mc", "ratio"};
  const int nsamples = 3;//5;
  const char* samples[nsamples] = {"gamjet","dijet", "zlljet"};//,"zeejet", "zmmjet", "zlljet", "dijet"};
  //const int nalphas = 4;
  //const int alphas[nalphas] = {30, 20, 15, 10};

  // If use uds, keep it last so can add up ud and s as needed
  const int nflavors = 7;
  const char *flavors[nflavors] = {"ud","s","c","b","g","oth","uds"};

  map<string, map<string, string> > rename;
  rename["gamjet"]["oth"] = "Undefined";
  rename["gamjet"]["uds"] = "uds";
  rename["gamjet"]["ud"] = "ud";
  rename["gamjet"]["s"] = "s";
  rename["gamjet"]["c"] = "c";
  rename["gamjet"]["b"] = "b";
  rename["gamjet"]["g"] = "glu";
  rename["dijet"]["oth"] = "undefined";
  rename["dijet"]["uds"] = "uds";
  rename["dijet"]["ud"] = "ud";
  rename["dijet"]["s"] = "s";
  rename["dijet"]["c"] = "c";
  rename["dijet"]["b"] = "b";
  rename["dijet"]["g"] = "glu"; // or gluExt
  rename["zlljet"]["oth"] = "Undefined";
  rename["zlljet"]["uds"] = "uds";
  rename["zlljet"]["ud"] = "ud";
  rename["zlljet"]["s"] = "s";
  rename["zlljet"]["c"] = "c";
  rename["zlljet"]["b"] = "b";
  rename["zlljet"]["g"] = "glu";

  // color and style codes
  map<string, map<string, int> > marker;
  //marker["dijet"]["ud"] = kFullCircle;
  //marker["gamjet"]["ud"] = kFullSquare;
  //marker["zlljet"]["ud"] = kFullDiamond;
  marker["dijet"]["oth"] = kFullDotLarge;
  marker["dijet"]["uds"] = kFullDotLarge;
  marker["dijet"]["ud"] = kFullDotLarge;
  marker["dijet"]["s"] = kFullDotLarge;
  marker["dijet"]["c"] = kFullDotLarge;
  marker["dijet"]["b"] = kFullDotLarge;
  marker["dijet"]["g"] = kFullDotLarge;
  marker["gamjet"]["oth"] = kFullSquare;
  marker["gamjet"]["uds"] = kFullSquare;
  marker["gamjet"]["ud"] = kFullSquare;
  marker["gamjet"]["s"] = kFullSquare;
  marker["gamjet"]["c"] = kFullSquare;
  marker["gamjet"]["b"] = kFullSquare;
  marker["gamjet"]["g"] = kFullSquare;
  marker["zlljet"]["oth"] = kFullDiamond;
  marker["zlljet"]["uds"] = kFullDiamond;
  marker["zlljet"]["ud"] = kFullDiamond;
  marker["zlljet"]["s"] = kFullDiamond;
  marker["zlljet"]["c"] = kFullDiamond;
  marker["zlljet"]["b"] = kFullDiamond;
  marker["zlljet"]["g"] = kFullDiamond;


  map<string, map<string, int> > color;
  //color["dijet"]["ud"] = kBlack;
  //color["gamjet"]["ud"] = kBlue;
  //color["zlljet"]["ud"] = kMagenta+2;
  color["dijet"]["oth"] = kGray;
  color["dijet"]["uds"] = kMagenta+2;
  color["dijet"]["ud"] = kBlue-9;
  color["dijet"]["s"] = kCyan+2;
  color["dijet"]["c"] = kGreen+2;
  color["dijet"]["b"] = kRed;
  color["dijet"]["g"] = kBlue;//kBlack;
  color["gamjet"]["oth"] = kGray;
  color["gamjet"]["uds"] = kMagenta+2;
  color["gamjet"]["ud"] = kBlue-9;
  color["gamjet"]["s"] = kCyan+2;
  color["gamjet"]["c"] = kGreen+2;
  color["gamjet"]["b"] = kRed;
  color["gamjet"]["g"] = kBlue;//kBlack;
  color["zlljet"]["oth"] = kGray;
  color["zlljet"]["uds"] = kMagenta+2;
  color["zlljet"]["ud"] = kBlue-9;
  color["zlljet"]["s"] = kCyan+2;
  color["zlljet"]["c"] = kGreen+2;
  color["zlljet"]["b"] = kRed;
  color["zlljet"]["g"] = kBlue;//kBlack;

  string sbin = Form("eta%02.0f-%02.0f",10*etamin,10*etamax);
  const char* bin = sbin.c_str();

 
  for (int idir = 0; idir != ndirs; ++idir) {

    cout << "  " <<  dirs[idir] << " - " << sbin << endl << flush;

    finout->cd();
    assert(finout->cd(dirs[idir]));
    TDirectory *din1 = finout->GetDirectory(dirs[idir]); assert(din1);

    assert(din1->cd(sbin.c_str()));
    TDirectory *din2 = din1->GetDirectory(sbin.c_str()); assert(din2);

    din2->cd();
    if (!din2->GetDirectory("flv"))
      din2->mkdir("flv");

    assert(din2->cd("flv"));
    TDirectory *dout = din2->GetDirectory("flv"); assert(dout);
    
    dout->cd();

    for (int isample = 0; isample != nsamples; ++isample) {
      
      string ss = samples[isample];

      cout << "    " <<  ss << endl << flush;

      if (ss=="gamjet") {

	assert(fp->cd("MC"));
	TDirectory *dp1 = fp->GetDirectory("MC"); assert(dp1);	

	// wider bins
	const int neta = 7;
	double etabins[neta+1] = {0, 8, 13, 19, 25, 30, 32, 52};
	// widest bins (0013 not available)
	//const int neta = 6;
	//double etabins[neta+1] = {0, 13, 19, 25, 30, 32, 52};

	TH1D *hp = new TH1D("hp","",neta,etabins);
	int ieta = hp->FindBin(0.5*(10*etamin+10*etamax));
	double eta1 = hp->GetBinLowEdge(ieta);
	double eta2 = hp->GetBinLowEdge(ieta+1);
	string sp = Form("eta%02d%02d",int(eta1+0.5),int(eta2+0.5));

	assert(dp1->cd(sp.c_str()));
	TDirectory *dp2 = dp1->GetDirectory(sp.c_str()); assert(dp2);

	TH1D *hsum(0);
	TH1F *hud(0), *hs(0);
	for (int iflv = 0; iflv != nflavors; ++iflv) {

	  string sf = flavors[iflv];
	  string sh = Form("flavFraction_%s_a30",rename[ss][sf].c_str());
	  TH1F *hf = (TH1F*)dp2->Get(sh.c_str());

	  cout << "      " <<  sf << endl << flush;

	  // Patch missing uds
	  // (binomial uncertainty can still be a big off)
	  if (sf=="uds" && !hf) {
	    assert(hud);
	    assert(hs);
	    hf = (TH1F*)hud->Clone("uds");
	    hf->Add(hs);
	  }
	  if (sf=="ud") hud = hf;
	  if (sf=="s")  hs  = hf;

	  if (!hf) cout << sh << endl << flush;
	  assert(hf);

	  dout->cd();
	  TH1D *h = new TH1D();
	  hf->Copy(*h);
	  sh = Form("%s_%s",ss.c_str(),sf.c_str());

	  // Fit shape (must fit before h written and deleted)
	  f1->SetParameters(0.5,0,0,0);
	  h->Fit(f1,"QRN");
	  cout << "    h# " << h << " hf# " << hf << endl;
	  pars[ss][sf].resize(f1->GetNpar());
	  for (int i = 0; i != f1->GetNpar(); ++i) {
	    pars[ss][sf][i] = f1->GetParameter(i);
	    cout << f1->GetParameter(i) << " ";
	  } // for i
	  cout << endl;
	  //f1->Write(Form("f1_%s_%s",ss.c_str(),sf.c_str()),TObject::kOverwrite);
	  h->SetNameTitle(sh.c_str(), sh.c_str());
	  h->SetLineColor(color[ss][sf]);
	  h->SetMarkerColor(color[ss][sf]);
	  h->SetMarkerStyle(marker[ss][sf]);
	  h->Write(sh.c_str(),TObject::kOverwrite);

	  if (sf!="uds") {
	    if (!hsum) hsum = (TH1D*)h->Clone("gamjet_sum");
	    else       hsum->Add(h);
	  }

	} // for iflv

	dout->cd();
	hsum->SetNameTitle("gamjet_sum","gamjet_sum");
	hsum->Write("gamjet_sum",TObject::kOverwrite);
      } // gamjet

      if (ss=="dijet") {

	assert(fd->cd("mc"));
	TDirectory *dd1 = fd->GetDirectory("mc"); assert(dd1);

	const int neta = 6;
	double etabins[neta+1] = {0, 13, 19, 25, 30, 32, 52};

	TH1D *hd = new TH1D("hd","",neta,etabins);
	int ieta = hd->FindBin(0.5*(10*etamin+10*etamax));
	double eta1 = hd->GetBinLowEdge(ieta);
	double eta2 = hd->GetBinLowEdge(ieta+1);
	string sd = Form("eta_%02d_%02d",int(eta1+0.5),int(eta2+0.5));

	assert(dd1->cd(sd.c_str()));
	TDirectory *dd2 = dd1->GetDirectory(sd.c_str()); assert(dd2);

	TH1D *hsum(0);
	for (int iflv = 0; iflv != nflavors; ++iflv) {

	  string sf = flavors[iflv];
	  string sh = Form("dijet_flavFraction_probejet_%s_a30",
			   rename[ss][sf].c_str());

	  TGraphErrors *gf = (TGraphErrors*)dd2->Get(sh.c_str());
	  if (!gf) cout << sh << endl;
	  assert(gf);
	  TH1D *hf = (TH1D*)dd2->Get("mc_RawNEvents_a30"); assert(hf);
	  //hf->Clear();
	  for (int i = 1; i != hf->GetNbinsX()+1; ++i) {
	    hf->SetBinContent(i,0);
	    hf->SetBinError(i,0);
	  }
	  
	  dout->cd();
	  TH1D *h = (TH1D*)hf->Clone("hf");
	  for (int i = 0; i != gf->GetN(); ++i) {
	    int ipt = h->FindBin(gf->GetX()[i]);
	    h->SetBinContent(ipt,gf->GetY()[i]);
	    h->SetBinError(ipt,gf->GetEY()[i]);
	  }
	  sh = Form("%s_%s",ss.c_str(),sf.c_str());

	  // Fit shape
	  f1->SetParameters(0.5,0,0,0);
	  h->Fit(f1,"QRN");
	  pars[ss][sf].resize(f1->GetNpar());
	  for (int i = 0; i != f1->GetNpar(); ++i) {
	    pars[ss][sf][i] = f1->GetParameter(i);
	  } // for i
	  //f1->Write(Form("f1_%s_%s",ss.c_str(),sf.c_str()),TObject::kOverwrite);

	  h->SetNameTitle(sh.c_str(), sh.c_str());
	  h->SetLineColor(color[ss][sf]);
	  h->SetMarkerColor(color[ss][sf]);
	  h->SetMarkerStyle(marker[ss][sf]);
	  h->Write(sh.c_str(),TObject::kOverwrite);

	  if (sf!="uds") { // don't double count ud and s
	    if (!hsum) hsum = (TH1D*)h->Clone("dijet_sum");
	    else       hsum->Add(h);
	  }

	}

	dout->cd();
	hsum->SetNameTitle("dijet_sum","dijet_sum");
	hsum->Write("dijet_sum",TObject::kOverwrite);
      } // dijet

      if (ss=="zlljet") {

	// fine bins
	//const int neta = 18;
	//double etabins[neta+1] = {0, 3, 5, 8, 10, 13, 15, 17, 19, 22, 23, 25,
	//			  27, 29, 30, 31, 35, 38, 52};
	// wider bins
	//const int neta = 7;
	//double etabins[neta+1] = {0, 8, 13, 19, 25, 30, 32, 52};
	// widest bins
	const int neta = 6;
	double etabins[neta+1] = {0, 13, 19, 25, 30, 32, 52};

	TH1D *hz = new TH1D("hz","",neta,etabins);
	int ieta = hz->FindBin(0.5*(10*etamin+10*etamax));
	double eta1 = hz->GetBinLowEdge(ieta);
	double eta2 = hz->GetBinLowEdge(ieta+1);
	string sz = Form("eta_%02d_%02d",int(eta1+0.5),int(eta2+0.5));

	TH1D *hsum(0);
	for (int iflv = 0; iflv != nflavors; ++iflv) {

	  string sf = flavors[iflv];
	  string sh = Form("MC_%s_CHS_a30_%s",rename[ss][sf].c_str(),
			   sz.c_str());
	  TH1F *hf = (TH1F*)fz->Get(sh.c_str());
	  if (!hf) cout << sh << endl;
	  assert(hf);

	  dout->cd();
	  //TH1D *h = (TH1D*)hf->Clone("hf");
	  //TH1D *h = new TH1D();
	  //hf->Copy(*h);

	  TH1D *hfn = (TH1D*)fz->Get(Form("MC_RawNEvents_CHS_a30_%s",
					  sz.c_str())); assert(hfn);
	  dout->cd();
	  // Copy NEvents for getting error bars and style ok,
	  // and fix bin uncertainties to binomial
	  TH1D *h = (TH1D*)hfn->Clone("hf");
	  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
	    int ipt = hf->FindBin(h->GetBinCenter(i));
	    double f = hf->GetBinContent(i);
	    double n = hfn->GetBinContent(i);
	    h->SetBinContent(ipt, f*n);
	    h->SetBinError(ipt, sqrt(f)*hfn->GetBinError(i));
	    //h->SetBinContent(ipt, f);
	    //h->SetBinError(ipt, f*(1-f)*hfn->GetBinError(i)
	    //		   / hfn->GetBinContent(i));
	  }
	  h->Divide(h,hfn);

	  sh = Form("%s_%s",ss.c_str(),sf.c_str());

	  // Fit shape
	  f1->SetParameters(0.5,0,0,0);
	  h->Fit(f1,"QRN");
	  pars[ss][sf].resize(f1->GetNpar());
	  for (int i = 0; i != f1->GetNpar(); ++i) {
	    pars[ss][sf][i] = f1->GetParameter(i);
	  } // for i
	  //f1->Write(Form("f1_%s_%s",ss.c_str(),sf.c_str()),TObject::kOverwrite);

	  h->SetNameTitle(sh.c_str(), sh.c_str());
	  h->SetLineColor(color[ss][sf]);
	  h->SetMarkerColor(color[ss][sf]);
	  h->SetMarkerStyle(marker[ss][sf]);
	  h->Write(sh.c_str(),TObject::kOverwrite);

	  if (sf!="uds") { // don't double count ud and s
	    if (!hsum) hsum = (TH1D*)h->Clone("zlljet_sum");
	    else       hsum->Add(h);
	  }

	}

	dout->cd();
	hsum->SetNameTitle("zlljet_sum","zlljet_sum");
	hsum->Write("zlljet_sum",TObject::kOverwrite);
      }
      
    } // for isample
  } // for idir
  
  finout->Close();
  fp->Close();
  fd->Close();
  fz->Close();
  
  //return pars;
}
