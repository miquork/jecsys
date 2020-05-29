// File: softrad.C
// Created by Mikko Voutilainen, on June 5th, 2014
// Updated on Oct 25, 2014 (cleanup for 8 TeV JEC paper)
// Purpose: Use JEC combination file to derive FSR+ISR corrections,
//          including uncertainty eigenvectors for global refit
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TF1.h"
#include "TMultiGraph.h"
#include "TLatex.h"
#include "TMinuit.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMath.h"
#include "TLine.h"

//#include "../tools.h"
#include "tdrstyle_mod15.C"

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>

const double _lumi = 2065.;//19800.;
// UL17: adjust threhsolds bit for FSR fits
const double ptrec_zlljet = 8.9;//5.;
const double ptrec_zjet = 15;
const double ptrec_gjet = 8.9;//15.;

bool debug = false;
bool dodijet = false;
bool domultijet = false; // =dodijet at |eta|<1.3
bool dropZee = true; // clean up plots a bit (leave Zll)
bool dropZmm = true; // clean up plots a bit  (leave Zll)

bool verbose = true;

using namespace std;

TF1 *_sr_fitError_func(0);
TMatrixD *_sr_fitError_emat(0);
Double_t sr_fitError(Double_t *xx, Double_t *p);

// Soft radiation corrections for L3Res
void softrad(double etamin=0.0, double etamax=1.3, bool dodijet=false,
	     string epoch="") {

  if (dodijet && fabs(etamin-0.0)<0.1 && fabs(etamax-1.3)<0.1)
    domultijet = true;
  
  cout << "Calling softrad("<<etamin<<","<<etamax<<","<<dodijet<<");\n"<<flush;
  const char *cep = epoch.c_str();
  cout << "For epoch " << epoch << endl;

  setTDRStyle();
  //writeExtraText = false; // for JEC paper CWR

  TDirectory *curdir = gDirectory;

  // Open jecdata.root produced by reprocess.C
  TFile *finout = new TFile(Form("rootfiles/jecdata%s.root",epoch.c_str()),
			    "UPDATE");
  assert(finout && !finout->IsZombie());
  curdir->cd();

  const int ndirs = 3;
  const char* dirs[ndirs] = {"data", "mc", "ratio"};
  const int nmethods = 2;
  const char* methods[nmethods] = {"mpfchs1", "ptchs"};
  //  const int nsamples = (dodijet ? 4 : 3);
  //const char* samples[4] = {"zeejet", "zmmjet", "zlljet", "dijet"};
  //const int nsamples = (dodijet ? 5 : 4);
  const int nsamples = (dodijet ? 6 : 5);
  //const char* samples[5] = {"gamjet","zeejet", "zmmjet", "zlljet",
  //const char* samples[5] = {"gamjet","zeejet", "zmmjet", "zjet",
  //			    domultijet ? "multijet" : "dijet"};
  const char* samples[6] = {"gamjet","zeejet", "zmmjet", "zlljet", "zjet",
			    domultijet ? "multijet" : "dijet"};
  //const int nsamples = (dodijet ? 4 : 3);
  //const char* samples[4] = {"gamjet","zeejet", "zmmjet", "dijet"};
  //const int nsamples = (dodijet ? 3 : 2);
  //const char* samples[3] = {"gamjet", "zmmjet", "dijet"};
  //const int nsamples = (dodijet ? 3 : 2);
  //const char* samples[3] = {"gamjet", "zlljet", "dijet"};
  const int idj = (dodijet ? nsamples-1 : -1);
  string sbin = Form("eta%02.0f-%02.0f",10*etamin,10*etamax);
  const char* bin = sbin.c_str();
  const int nalphas = 4;
  const int alphas[nalphas] = {30, 20, 15, 10};

  // Zll+jet bins
  const double ptbins1[] =
    {30, 40, 50, 60, 85, 105, 130, 175, 230, 300, 400, 500, 700, 1000, 1500};
  const int npt1 = sizeof(ptbins1)/sizeof(ptbins1[0])-1;
  double ptmax1 = ptbins1[npt1];
  TH1D *hpt1 = new TH1D("hpt1","",npt1,&ptbins1[0]);
  TProfile *ppt1 = new TProfile("ppt1","",npt1,&ptbins1[0]);

  // Z+jet bins
  const double ptbins1b[] =
    {30, 40, 50, 60, 70,85, 105, 130, 175, 230, 300, 400, 500, 700, 1000, 1500};
  const int npt1b = sizeof(ptbins1b)/sizeof(ptbins1b[0])-1;
  double ptmax1b = ptbins1b[npt1];
  TH1D *hpt1b = new TH1D("hpt1b","",npt1b,&ptbins1b[0]);
  TProfile *ppt1b = new TProfile("ppt1b","",npt1b,&ptbins1b[0]);


  // gamma+jet bins
  const double ptbins2[] = {40, 50, 60, 85, 105, 130, 175, 230, 300, 400,
  			    500, 700, 1000, 1500};
  //const double ptbins2[] = {40, 50, 60, 85, 105, 130, 175, 230, 300, 400, 500, 700, 1000, 3000}; // more correct, but too wide last bin
  const int npt2 = sizeof(ptbins2)/sizeof(ptbins2[0])-1;
  cout << Form("npt2: %i", npt2) << endl << flush;
  double ptmax2 = ptbins2[npt2];
  TH1D *hpt2 = new TH1D("hpt2","",npt2,&ptbins2[0]);
  TProfile *ppt2 = new TProfile("ppt2","",npt2,&ptbins2[0]);

  // dijet bins
  const double ptbins4[] = {20, 62, 107, 175, 242, 310, 379, 467,
			    628, 839, 1121, 1497, 2000};
  const int npt4 = sizeof(ptbins4)/sizeof(ptbins4[0])-1;
  TH1D *hpt4 = new TH1D("hpt4","",npt4,&ptbins4[0]);
  TProfile *ppt4 = new TProfile("ppt4","",npt4,&ptbins4[0]);

  // multijet bins
  const double ptbins5[] = {10, 15, 21, 28, 37, 49, 64, 84, 114, 153, 196, 245, 300, 362, 430, 507, 592, 686, 790, 905, 1032, 1172, 1327, 1497, 1684, 1890, 2116, 2366, 2640, 2941, 3273, 5220};
  const int npt5 = sizeof(ptbins5)/sizeof(ptbins5[0])-1;
  double ptmax5 = ptbins5[npt5];
  TH1D *hpt5 = new TH1D("hpt5","",npt5,&ptbins5[0]);
  TProfile *ppt5 = new TProfile("ppt5","",npt5,&ptbins5[0]);

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);

  map<string,const char*> texlabel;
  texlabel["gamjet"] = "#gamma+jet";
  texlabel["zeejet"] = "Z(#rightarrowee)+jet";// TB";
  //texlabel["zeejet"] = "Z(#rightarrowee)+jet RS";
  texlabel["zmmjet"] = "Z(#rightarrow#mu#mu)+jet";// TB";
  //texlabel["zmmjet"] = "Z(#rightarrow#mu#mu)+jet RS";
  texlabel["zlljet"] = "Z(#rightarrowl^{+}l^{-})+jet";
  texlabel["zjet"] = "Z+jet";
  texlabel["dijet"] = "Dijet";
  texlabel["multijet"] = "Multijet";
  texlabel["ptchs"] = "p_{T} balance (CHS)";
  texlabel["mpfchs"] = "MPF raw (CHS)";
  texlabel["mpfchs1"] = "MPF type-I (CHS)";

  // overlay of various alpha values
  TCanvas *c1 = new TCanvas("c1","c1",ndirs*400,nmethods*400);
  c1->Divide(ndirs,nmethods);

  TH1D *h1 = new TH1D("h1",";p_{T} (GeV);Response",2610,30,2640);
  h1->GetXaxis()->SetNoExponent();
  h1->GetXaxis()->SetMoreLogLabels();

  // extrapolation vs alpha for each pT bin
  vector<TCanvas*> c2s(ndirs*nmethods);
  for (unsigned int icanvas = 0; icanvas != c2s.size(); ++icanvas) {
    TCanvas *c2 = new TCanvas(Form("c2_%d",icanvas),Form("c2_%d",icanvas),
			      1200,1200);
    c2->Divide(3,3);
    c2s[icanvas] = c2;
  }

  TH1D *h2 = new TH1D("h2",";#alpha;Response",10,0.,0.4);
  h2->SetMaximum(1.10);//1.08);
  h2->SetMinimum(0.85);//0.88);

  // krad corrections
  TCanvas *c3 = new TCanvas("c3","c3",ndirs*400,nmethods*400);
  c3->Divide(ndirs,nmethods);

  TH1D *h3 = new TH1D("h3",";p_{T,ref} (GeV);FSR sensitivity: -dR/d#alpha [%]",
		      2610,30,2640);
  h3->GetXaxis()->SetNoExponent();
  h3->GetXaxis()->SetMoreLogLabels();

  cout << "Reading in data" << endl << flush;
  // Read in plots vs pT (and alpha)
  map<string, map<string, map<string, map<int, TGraphErrors*> > > > gemap;
  map<string, map<string, map<string, map<int, TGraphErrors*> > > > gamap;
  for (int idir = 0; idir != ndirs; ++idir) {
    for (int  imethod = 0; imethod != nmethods; ++imethod) {
      for (int  isample = 0; isample != nsamples; ++isample) {

        finout->cd();
        assert(finout->cd(dirs[idir]));
        TDirectory *din1 = finout->GetDirectory(dirs[idir]); assert(din1);
        assert(din1->cd(bin));
        TDirectory *d = din1->GetDirectory(bin); assert(d);
        d->cd();

        const char *cd = dirs[idir];
        const char *cm = methods[imethod];
        const char *cs = samples[isample];
	string ss = samples[isample];

	for (int  ialpha = 0; ialpha != nalphas; ++ialpha) {
	  const int a = alphas[ialpha];
	  // Get graph made vs pT
	  string s = Form("%s/%s/%s_%s_a%d",dirs[idir],bin,cm,cs,a);
	  TGraphErrors *g = (TGraphErrors*)finout->Get(s.c_str());
	  if (!g) cout << "Missing " << s << endl << flush;
	  if (!g && ss=="multijet" && a==10) g = new TGraphErrors(0);
	  assert(g);

	  // Clean out empty points
	  // as well as trigger-biased ones for dijets
	  // as well as weird gamma+jet high pT point
	  // as well as half-empty Z+jet high pT points
	  for (int i = g->GetN()-1; i != -1; --i) {
            //cout << g->GetX()[i] << " ";
	    if (g->GetY()[i]==0 || g->GetEY()[i]==0 ||
		(ss=="dijet" && g->GetX()[i]<70.)  ||
		(ss=="multijet" && g->GetX()[i]>ptmax5)  ||
		//(ss=="multijet" && g->GetX()[i]<49.)  ||
		//(ss=="multijet" && g->GetX()[i]<84.)  ||
		(ss=="gamjet" && g->GetX()[i]>ptmax2) ||
		(ss=="zeejet" && g->GetX()[i]>ptmax1) ||
		(ss=="zmmjet" && g->GetX()[i]>ptmax1) ||
		(ss=="zlljet" && g->GetX()[i]>ptmax1) ||
		(ss=="zjet" && g->GetX()[i]>ptmax1b))
	      g->RemovePoint(i);
	  }

	  gemap[cd][cm][cs][a] = g;
	  
	  // Sort points into new graphs vs alpha
	  //TH1D *hpt = (isample==0 ? hpt2 : hpt1);
	  //TProfile *ppt = (isample==0 ? ppt2 : ppt1);
	  TH1D *hpt = (ss=="gamjet" ? hpt2 : 
		       (ss=="zjet" ? hpt1b : hpt1));
	  TProfile *ppt = (ss=="gamjet" ? ppt2 : 
			   (ss=="zjet" ? ppt1b : ppt1));
	  //if (isample==3) { hpt = hpt4; ppt = ppt4; } // pas-v6
	  //if (isample==idj) { hpt = hpt4; ppt = ppt4; } // pas-v6
	  if (ss=="dijet") { hpt = hpt4; ppt = ppt4; } // pas-v6
	  if (ss=="multijet") { hpt = hpt5; ppt = ppt5; }
	  for (int i = 0; i != g->GetN(); ++i) {
	    
	    double pt = g->GetX()[i];
	    ppt->Fill(pt, pt);
            //if(verbose)cout << Form("Fill %s %s ppt bin %i with pt %6.2f", samples[isample], methods[imethod],i,pt) <<endl << flush; //this also limits the range of kfsr-plotting
	    int ipt = int(hpt->GetBinLowEdge(hpt->FindBin(pt))+0.5);
	    //int ipt = int(pt+0.5);
	    TGraphErrors *ga = gamap[cd][cm][cs][ipt];
	    if (!ga) {
	      ga = new TGraphErrors(0);
	      ga->SetMarkerStyle(g->GetMarkerStyle());
	      ga->SetMarkerColor(g->GetMarkerColor());
	      ga->SetLineColor(g->GetLineColor());
	      gamap[cd][cm][cs][ipt] = ga;
	    }
	    int n = ga->GetN();
	    ga->SetPoint(n, 0.01*a, g->GetY()[i]);
	    ga->SetPointError(n, 0, g->GetEY()[i]);
	  } // for i 

	} // for ialpha
        
        map<int, TGraphErrors*> &gam = gamap[cd][cm][cs];
        map<int, TGraphErrors*>::iterator itpt;
        for (itpt = gam.begin(); itpt != gam.end(); ++itpt) {
	  int ipt = itpt->first;
	  int jpt = hpt1->FindBin(ipt);
	  
	  TGraphErrors *ga = itpt->second; assert(ga);

          //patch here
          
	  // Clean out  points at low alpha that have lower uncertainty than next one
	  for (int i = ga->GetN()-1; i != -1; --i) {
            if (debug) cout << ga->GetX()[i] << " " << ga->GetY()[i] << " i " << i <<"  ga->GetN() " << ga->GetN() << endl;
            if (i>0 && ga->GetX()[i]!=ga->GetX()[i-1] && ga->GetEY()[i]<0.5*ga->GetEY()[i-1]){
              cout << ga->GetX()[i] << " " << ga->GetY()[i]  << " remove point because ga->GetEY()[i]<0.5*ga->GetEY()[i-1] " << ga->GetEY()[i] << " < 0.5* " << ga->GetEY()[i-1] << endl;
              ga->RemovePoint(i);
            }
	  }
          //if (ga->GetN()<3){
          //  ga->Set(0);
          //  cout << "removing graph" <<endl;
	  //}

        }
      } // for isample
    } // for imethod
  } // for idir

  cout << "Drawing plots vs pT for each alpha" << endl << flush;

  // 2x6 plots
  for (int idir = 0; idir != ndirs; ++idir) {
    for (int  imethod = 0; imethod != nmethods; ++imethod) {

      const char *cd = dirs[idir];
      const char *cm = methods[imethod];

      int ipad = ndirs*imethod + idir + 1; assert(ipad<=6);
      c1->cd(ipad);
      gPad->SetLogx();
      h1->SetMaximum(1.15);//idir<2 ? 1.15 : 1.08);
      h1->SetMinimum(0.80);//idir<2 ? 0.85 : 0.93);
      h1->SetYTitle(Form("Response (%s)",cd));
      h1->DrawClone("AXIS");
      tex->DrawLatex(0.20,0.85,texlabel[cm]);
      //if (epoch!="L4") tex->DrawLatex(0.20,0.80,"|#eta| < 1.3, #alpha=0.1--0.3");
      //if (epoch=="L4") tex->DrawLatex(0.20,0.80,"|#eta| < 2.4, #alpha=0.1--0.3");
      if (etamin==0)
	tex->DrawLatex(0.20,0.80,Form("|#eta| < %1.1f, #alpha=0.1--0.3",etamax));
      else
	tex->DrawLatex(0.20,0.80,Form("%1.1f < |#eta| < %1.1f, #alpha=0.1--0.3",
				      etamin, etamax));
      TLegend *leg = tdrLeg(0.60,0.75,0.90,0.90);

      for (int  isample = 0; isample != nsamples; ++isample) {
	for (int  ialpha = 0; ialpha != nalphas; ++ialpha) {

	  const char *cs = samples[isample];
	  string ss = cs;
	  if (etamin==0 && ss=="dijet") continue;

	  const int a = alphas[ialpha];
	  TGraphErrors *g = gemap[cd][cm][cs][a]; assert(g);

	  // Clean out points with very large uncertainty for plot readability
	  //for (int i = g->GetN()-1; i != -1; --i) {
	  //if (g->GetEY()[i]>0.02) g->RemovePoint(i);
	  //}
	  
	  if (!((dropZee && ss=="zeejet") || (dropZmm && ss=="zmmjet"))) {
	    g->Draw("SAME Pz");
	    
	    if (ialpha==0) leg->AddEntry(g,texlabel[cs],"P");
	  }
	}
      } // for isample

      // Individual plots for JEC paper
      if ( true ) { // paper

	TH1D *h = new TH1D(Form("h_5%s_%s",cd,cm),
			   Form(";p_{T} (GeV);Response (%s)",cd),
			   //1470,30,1500);
			   1270,30,1300); // Run I paper
	h->GetXaxis()->SetMoreLogLabels();
	h->GetXaxis()->SetNoExponent();
	h->SetMinimum(0.88);//0.83);//0.88);
	h->SetMaximum(1.13);//1.13);

	// better range for pT balance plots
	if (string(cm)=="ptchs" && string(cd)!="ratio") {
	  h->SetMinimum(0.83);
	  h->SetMaximum(1.08);
	}
	
	//map<string, const char*> lumimap;
	//lumimap["BCD"] = "Run2016BCD Legacy, 12.9 fb^{-1}";
	//lumimap["E"] = "Run2016E Legacy, 4.0 fb^{-1}";
	//lumimap["F"] = "Run2016F Legacy, 2.8 fb^{-1}";//3.1 fb^{-1}";
	//lumimap["G"] = "Run2016fG Legacy, 8.0 fb^{-1}";
	//lumimap["H"] = "Run2016H Legacy, 8.8 fb^{-1}";
	//lumimap["GH"] = "Run2016fGH Legacy, 16.8 fb^{-1}";
	//lumimap["BCDEF"] = "Run2016BCDEF Legacy, 19.7 fb^{-1}";
	//lumimap["BCDEFGH"] = "Run2016BCDEFGH Legacy, 36.5 fb^{-1}";
	//lumimap["EF"] = "Run2016EF Legacy, 6.8 fb^{-1}";
	//lumimap["L4"] = "Run2016BCDEFGH closure, 36.5 fb^{-1}";
        map<string, const char*> lumimap;
	/*
        lumimap["A"] = "Run2018A 14.0 fb^{-1}"; //PdmV Analysis TWiki
        lumimap["B"] = "Run2018B 7.1 fb^{-1}"; //PdmV Analysis TWiki
        lumimap["C"] = "Run2018C 6.9 fb^{-1}"; //PdmV Analysis TWiki
        lumimap["D"] = "Run2018D 31.9 fb^{-1}"; //PdmV Analysis TWiki
        lumimap["ABC"] = "Run2018ABC 28.0 fb^{-1}"; //PdmV Analysis TWiki
        lumimap["ABCD"] = "Run2018ABCD 59.9 fb^{-1}"; //PdmV Analysis TWiki
	*/
	lumimap["BCDEF"] = "2017, 41.5 fb^{-1}"; // for DP note
	lumimap["B"] = "Run2017B, 4.8 fb^{-1}";
	lumimap["C"] = "Run2017C, 9.6 fb^{-1}";
	lumimap["D"] = "Run2017D, 4.2 fb^{-1}";
	lumimap["E"] = "Run2017E, 9.3 fb^{-1}";
	lumimap["F"] = "Run2017F, 13.4 fb^{-1}";
	lumi_13TeV = lumimap[epoch];

	TCanvas *c0 = tdrCanvas(Form("c0_%s_%s",cm,cd), h, 4, 11, true);
	c0->SetLogx();
	

	TLegend *leg = tdrLeg(0.55,0.68,0.85,0.83);
	tex->DrawLatex(0.55,0.85,texlabel[cm]);
	//if (epoch!="L4") tex->DrawLatex(0.55,0.18,"|#eta| < 1.3, #alpha < 0.3");
	//if (epoch=="L4") tex->DrawLatex(0.55,0.18,"|#eta| < 2.4, #alpha < 0.3");
	if (etamin==0)
	  tex->DrawLatex(0.55,0.18,Form("|#eta| < %1.1f, #alpha < 0.3",etamax));
	else
	  tex->DrawLatex(0.55,0.18,Form("%1.1f < |#eta| < %1.1f, #alpha < 0.3",
					etamin, etamax));
	//tex->DrawLatex(0.55,0.18,"Anti-k_{T} R=0.5");

	// Loop over Z+jet and gamma+jet (only, no dijet/multijet)
	//for (int  isample = 0; isample != min(3,nsamples); ++isample) {
	//for (int  isample = 0; isample != min(idj,nsamples); ++isample) {
	for (int  isample = 0; isample != nsamples; ++isample) {
	  
	  const char *cs = samples[isample];
	  string ss = cs;
	  if (ss=="dijet") continue;
	  TGraphErrors *g = gemap[cd][cm][cs][30];
	  if (!g) cout << cd <<"_"<< cm <<"_"<< cs <<"_30"<< endl << flush;
	  assert(g);
	  if (!((dropZee && ss=="zeejet") || (dropZmm && ss=="zmmjet"))) {
	    g->Draw("SAME Pz");
	  
	    leg->AddEntry(g,texlabel[cs],"P");
	  }
	} // for isample

	if (etamin==0 && (fabs(etamax-1.3)<0.1 || fabs(etamax-2.4)<0.1)) {
	  c0->SaveAs(Form("pdf/%s/paper_softrad_%s_%s_vspt.pdf",cep,cd,cm));
	  //c0->SaveAs(Form("pdfC/paper_softrad_%s_%s_vspt.C",cd,cm));
	}
	else {
	  c0->SaveAs(Form("pdf/%s/an_softrad_%s_%s_eta%02.0f-%02.0f_vspt.pdf",
			  cep,cd,cm,10*etamin,10*etamax));
	}
      } // paper
    } // for imethod
  } // for idir
  
  c1->cd(0);
  //cmsPrel(_lumi, true);
  //CMS_lumi(c1, 2, 33);
  c1->SaveAs(Form("pdf/%s/softrad_2x6_vspt_eta%02.0f-%02.0f.pdf",cep,10*etamin,10*etamax));


  cout << "Drawing plots vs alpha for each pT" << endl << flush;
  cout << "...and fitting slope vs alpha" << endl << flush;

  map<string, map<string, map<string, TGraphErrors* > > > gkmap;

  // 2x6 plots
  for (int idir = 0; idir != ndirs; ++idir) {
    for (int  imethod = 0; imethod != nmethods; ++imethod) {
      
      int icanvas = nmethods*imethod + idir; assert(icanvas<=6);
      TCanvas *c2 = c2s[icanvas]; assert(c2);

      const char *cd = dirs[idir];
      string dd = cd;
      const char *cm = methods[imethod];
      string mm = cm;

      const int npads = 9;
      for (int ipad = 0; ipad != npads; ++ipad) {
	c2->cd(ipad+1);
	h2->SetYTitle(Form("Response (%s)",cd));
	h2->DrawClone("AXIS");
	tex->DrawLatex(0.20,0.85,texlabel[cm]);
	//tex->DrawLatex(0.20,0.80,epoch=="L4" ? "|#eta| < 2.4" : "|#eta| < 1.3");
	if (etamin==0)
	  tex->DrawLatex(0.20,0.80,Form("|#eta| < %1.1f",etamax));
	else
	  tex->DrawLatex(0.20,0.80,Form("%1.1f < |#eta| < %1.1f",
					etamin,etamax));
	tex->DrawLatex(0.20,0.75,Form("%1.0f < p_{T} < %1.0f GeV",
				      hpt1->GetBinLowEdge(ipad+1),
				      hpt1->GetBinLowEdge(ipad+2)));
	TLegend *leg = tdrLeg(0.65,0.75,0.90,0.90);
	leg->AddEntry(gemap[cd][cm]["gamjet"][30], texlabel["gamjet"], "P");
	leg->AddEntry(gemap[cd][cm]["zeejet"][30], texlabel["zeejet"], "P");
	leg->AddEntry(gemap[cd][cm]["zmmjet"][30], texlabel["zmmjet"], "P");
	leg->AddEntry(gemap[cd][cm]["zlljet"][30], texlabel["zlljet"], "P");
	leg->AddEntry(gemap[cd][cm]["zjet"][30],   texlabel["zjet"], "P");
	if (domultijet)
	  leg->AddEntry(gemap[cd][cm]["multijet"][30],texlabel["multijet"],"P");
	else if (dodijet)
	  leg->AddEntry(gemap[cd][cm]["dijet"][30], texlabel["dijet"], "P");
      }

      for (int  isample = 0; isample != nsamples; ++isample) {

	const char *cs = samples[isample];
	string ss = cs;

	map<int, TGraphErrors*> &gam = gamap[cd][cm][cs];
	map<int, TGraphErrors*>::iterator itpt;
	for (itpt = gam.begin(); itpt != gam.end(); ++itpt) {

	  int ipt = itpt->first;
	  int jpt = hpt1->FindBin(ipt);
	  //double ptmin = hpt1->GetBinLowEdge(jpt); // 2017-06-13
	  //if (jpt>npads) continue;
	  //assert(jpt<=npads);
	  //c2->cd(jpt);
	  c2->cd(min(npads,jpt));
	  
	  TGraphErrors *ga = itpt->second; assert(ga);
	  
	  if (!((dropZee && ss=="zeejet") || (dropZmm && ss=="zmmjet")))
	    ga->Draw("SAME Pz");

	  // Fit slope
	  TF1 *f1 = new TF1(Form("f1_%s_%s_%s_%d",cd,cm,cs,ipt),
			    "(x<1)*([0]+[1]*x) + (x>1 && x<2)*[0] +"
			    "(x>2)*[1]",-1,1);
	  f1->SetLineColor(ga->GetLineColor());
	  f1->SetParameters(1,0);
	  double minalpha = (isample==0 ? ptrec_gjet/ipt :
			     (ss=="zjet" ? ptrec_zjet/ipt : ptrec_zlljet/ipt));
	  if (ss=="multijet") minalpha=0.145; // 0.15 included
	  // Constrain slope to within reasonable values
	  // in the absence of sufficient data using priors
	  // apply at max(ptrec_gjet,ptrec_zjet)/0.15 = 59.3 GeV ~ 60 GeV
	  if (ipt<60) { // use priors
	    int n = ga->GetN();
	    bool iszjet=(ss=="zeejet"||ss=="zmmjet"||ss=="zlljet"||ss=="zjet");
	    if (mm=="ptchs") {
	      // For pT balance in UL17, looking at higher pT:
	      // Best fit for Z+jet data is -30%, MC -25%, ratio -7%
	      // Best fit for gamma+jet data and MC is -25%, ratio -3%
	      // Data,MC variation 5%, ratio variation 2%(gamma+jet),7%(Z+jet)
	      // Use two sigma of above for constraints
	      if (dd=="data"  && iszjet) ga->SetPoint(n,   2.5, -0.30);
	      if (dd=="data"  && iszjet) ga->SetPointError(n,0,  0.05*2);
	      if (dd=="mc"    && iszjet) ga->SetPoint(n,   2.5, -0.25);
	      if (dd=="mc"    && iszjet) ga->SetPointError(n,0,  0.05*2);
	      if (dd=="ratio" && iszjet) ga->SetPoint(n,   2.5, -0.07);
	      if (dd=="ratio" && iszjet) ga->SetPointError(n,0,  0.07*2);
	      //
	      if (dd=="data"  && ss=="gamjet") ga->SetPoint(n,   2.5, -0.25);
	      if (dd=="data"  && ss=="gamjet") ga->SetPointError(n,0,  0.05*2);
	      if (dd=="mc"    && ss=="gamjet") ga->SetPoint(n,   2.5, -0.25);
	      if (dd=="mc"    && ss=="gamjet") ga->SetPointError(n,0,  0.05*2);
	      if (dd=="ratio" && ss=="gamjet") ga->SetPoint(n,   2.5, -0.03);
	      if (dd=="ratio" && ss=="gamjet") ga->SetPointError(n,0,  0.02*2);
	    }
	    if (mm=="mpfchs1") { // MPF
	      // For MPF, expectation is no slope
	      // Maximal slope would be approximately
	      // (<vecpT2>/alpha ~ 25% from pT balance) times
	      // (response difference between pT1 and vecpT2~10%)
	      // => 0.25*0.10 = 2.5%
	      // For data/MC, estimate uncertainty as half of this
	      // => 1.25%
	      // UL17 constraint, no slope in MPF, two sigma of above
	      if (iszjet || ss=="gamjet") {
		ga->SetPoint(n, 2.5, dd=="mc" ? 0. : 0.00);
		if (dd=="data")  ga->SetPointError(n, 0, 0.050);
		if (dd=="mc")    ga->SetPointError(n, 0, 0.05);
		if (dd=="ratio") ga->SetPointError(n, 0, 0.025);
	      }
	    } // MPF
	  } // use priors

	  // Create graph to store results
	  TGraphErrors *gk = gkmap[cd][cm][cs];
	  if (!gk) {
	    gk = new TGraphErrors(0);
	    gk->SetMarkerStyle(ga->GetMarkerStyle());
	    gk->SetMarkerColor(ga->GetMarkerColor());
	    gk->SetLineColor(ga->GetLineColor());
	    gkmap[cd][cm][cs] = gk;
	  }

	  if (ga->GetN()>2) {

	    f1->SetRange(minalpha, 3.);
	    ga->Fit(f1,"QRN");

	    if (f1->GetNDF()>=0) {
	      if (!((dropZee && ss=="zeejet") || (dropZmm && ss=="zmmjet")))
		f1->DrawClone("SAME");
	      f1->SetRange(0,0.4);
	      f1->SetLineStyle(kDashed);
	      if (!((dropZee && ss=="zeejet") || (dropZmm && ss=="zmmjet")))
		f1->DrawClone("SAME");

	      // Store results
	      TGraphErrors *gk = gkmap[cd][cm][cs];
	      int n = gk->GetN();

	      TProfile *ppt(0);// = (isample==0 ? ppt2 : ppt1);
	      //if (isample==3) { ppt = ppt4; } // pas-v6
	      //if (isample==idj) { ppt = ppt4; } // pas-v6
	      if (ss=="gamjet") { ppt = ppt2; }
	      if (ss=="zlljet")   { ppt = ppt1; }
	      if (ss=="zjet")     { ppt = ppt1b; }
	      if (ss=="zeejet")   { ppt = ppt1; }
	      if (ss=="zmmjet")   { ppt = ppt1; }
	      if (ss=="dijet")    { ppt = ppt4; }
	      if (ss=="multijet") { ppt = ppt5; }
	      assert(ppt);
	      double pt = ppt->GetBinContent(ppt->FindBin(ipt));
	      gk->SetPoint(n, pt, f1->GetParameter(1));
	      gk->SetPointError(n, 0, f1->GetParError(1));
	    } // f1->GetNDF()>=0
	    else {
	      cout << " f1->GetNDF()<0 in "<<cd<<" "<<cs<<" "<<cm<<endl<<flush; 
	    }
	  } // ga->GetN()>2
	  else {
	    cout << " ga->GetN()<=2 for "<<cd<<" "<<cs<<" "<<cm<<endl<<flush; 
	  }
	} // for itpt
	
      } // for isample
      
      c2->SaveAs(Form("pdf/%s/softrad_3x3_%s_%s_vsalpha_eta%02.0f-%02.0f.pdf",cep,cd,cm,10*etamin,10*etamax));
      
    } // for method
  } // for dir


  cout << "Drawing plots of kFSR vs pT" << endl;

  // 2x6 plots
  for (int idir = 0; idir != ndirs; ++idir) {
    for (int  imethod = 0; imethod != nmethods; ++imethod) {

      const char *cd = dirs[idir];
      const char *cm = methods[imethod];
      string sm(cm);

      TMultiGraph *mgk = new TMultiGraph();

      int ipad = ndirs*imethod + idir + 1; assert(ipad<=6);
      c3->cd(ipad);
      gPad->SetLogx();
      //h3->SetMaximum(imethod==0 ? 0.05 : (idir!=2 ? 0.1 : 0.25));
      //h3->SetMinimum(imethod==0 ? -0.05 : (idir!=2 ? -0.4 : -0.25));
      h3->SetMaximum(imethod==0 ? 0.10 : (idir!=2 ? 0.45 : 0.25));
      h3->SetMinimum(imethod==0 ? -0.10 : (idir!=2 ? -0.55 : -0.25));
      h3->SetYTitle(Form("k_{FSR} = dR/d#alpha (%s)",cd));
      h3->DrawClone("AXIS");
      TLine *l = new TLine();
      l->SetLineStyle(kDashed);
      l->DrawLine(h3->GetXaxis()->GetXmin(),0,h3->GetXaxis()->GetXmax(),0);
      tex->DrawLatex(0.20,0.85,texlabel[cm]);
      //if (epoch!="L4") tex->DrawLatex(0.20,0.80,"|#eta| < 1.3");
      //if (epoch=="L4") tex->DrawLatex(0.20,0.80,"|#eta| < 2.4");
      if (etamin==0)
	tex->DrawLatex(0.20,0.80,Form("|#eta| < %1.1f",etamax));
      else
	tex->DrawLatex(0.20,0.80,Form("%1.1f < |#eta| < %1.1f",
				      etamin,etamax));
      tex->DrawLatex(0.20,0.75,epoch.c_str());
      TLegend *leg = tdrLeg(0.60,0.75,0.90,0.90);

      for (int  isample = 0; isample != nsamples; ++isample) {

	const char *cs = samples[isample];
	string ss(cs);
	TGraphErrors *gk = gkmap[cd][cm][cs];
	if (!gk) cout << cd << " " << cm << " " << cs
		      << "eta_"<<etamin<<"_"<<etamax << endl << flush;
	assert(gk);
	
	if (!((dropZee && ss=="zeejet") || (dropZmm && ss=="zmmjet")))
	  leg->AddEntry(gk,texlabel[cs],"P");

	// Fit each sample separately for pT balance
	if (true) {

	  TF1 *fk(0);
	  // For UL17_V2
	  if ((ss=="zlljet" || ss=="zeejet" || ss=="zmmjet")||ss=="zjet") {
	    // Log-lin works well, except maybe for zlljet pTbal in data
	    // Want to reduce freedom a bit at the edges of phase space
	    if (sm=="ptchs")
	      fk = new TF1(Form("fk_%s_%s_%s",cd,cm,cs),
			   "[0]+[1]*log(0.01*x)",30,3000);
			   //"[0]+[1]*log(0.01*x)+[2]*pow(log(0.01*x),2)",30,3000);
			   //"[0]+[1]/log(x/0.218)",30,3000);
	    if (sm=="mpfchs1")
	      // Similar to ptchs, but now scaled to pt>15 GeV part only
	      fk = new TF1(Form("fk_%s_%s_%s",cd,cm,cs),
			   "([0]+[1]*log(0.01*x))*(15./x)/0.3",30,3000);
			   //"([0]+[1]*log(0.01*x)+[2]*pow(log(0.01*x),2))*(15./x)/0.3",30,3000);
			   //"([0]+[1]/log(x/0.218))*(15./x)/0.3",30,3000);
	  }
	  if (ss=="gamjet") {
	    
	    if (sm=="ptchs")
	      // Log-lin works well, especially for short span at pT>230 GeV
	      fk = new TF1(Form("fk_%s_%s_%s",cd,cm,cs),
			   "[0]+[1]*log(0.01*x)",30,3000);
	    if (sm=="mpfchs1")
	      // Similar to ptchs, but now scaled to pt>15 GeV part only
	      fk = new TF1(Form("fk_%s_%s_%s",cd,cm,cs),
			   "([0]+[1]*log(0.01*x))*(15./x)/0.3",30,3000);
	  }
	  if (ss=="multijet" || ss=="dijet") {
	    // Log-lin ok for high pT, but low pT needs extra shaping
	    //fk = new TF1(Form("fk_%s_%s_%s",cd,cm,cs),
	    //		 "[0]+[1]*log(0.01*x)+[2]/(log(x/0.218) * x)",30,3000);
			 //"[0]+[1]*log(0.01*x)+[2]*pow(log(0.01*x),2)"
			 //"+ [3] / (log(x/0.218) * x)",30,3000);
	    fk = new TF1(Form("fk_%s_%s_%s",cd,cm,cs),
			 "[0]+[1]/(log(x/0.218) * x)",30,3000);
	  }
	  assert(fk);
	  /*
	  if (ss=="multijet") {// || sm=="mpfchs1") {
	    // MJB lead: proportional to alpha_s times 15 GeV / pT
	    // MPF lead: proportional to MJB times unclustered pT response diff.
	    // both simply to 1/(log(x/LambdaQCD)*x) shape at high x
	    // but, code does not like a single parameter, need at least two?
	    // alpha_s + quad-log shape?
	    // globalFitL3Res.C may also be confused by other than 3 pars
	    fk = new TF1(Form("fk_%s_%s_%s",cd,cm,cs),
			 //"[0] / (log(x/0.218) * x) + [1]/x + [2]*log(0.01*x)",
			 //"[0] / (log(x/0.218) * x) + [1] + [2]*log(0.01*x)",
			 "[0]+[1]*log(0.01*x)+[2]*pow(log(0.01*x),2)"
			 "+ [3] / (log(x/0.218) * x)",
			 30,3000);
	    //if (sm=="ptchs")   fk->SetParameters(45,-0.05,+0.01,0);
	    //if (sm=="mpfchs1") fk->SetParameters(45*0.2,-0.05,+0.01,0);
	    //if (sm=="ptchs")   fk->SetParameters(45,+0.01,0);
	    //if (sm=="mpfchs1") fk->SetParameters(45*0.2,+0.01,0);
	    if (sm=="ptchs")   fk->SetParameters(-0.05,+0.01,0,45);
	    if (sm=="mpfchs1") fk->SetParameters(-0.05,+0.01,0,45*0.2);
	  }
	  else {
	    // Z+jet still better with quad-log
	    fk = new TF1(Form("fk_%s_%s_%s",cd,cm,cs),
			 "[0]+[1]*log(0.01*x)+[2]*pow(log(0.01*x),2)",
			 30,1500);
	    fk->SetParameters(-0.05,+0.01,0);
	    //fk->FixParameter(2,0.); // post-Moriond19, reduce quad-log to log0-lin => UL17 again quad-log
	  }
	  */

	  double dy = (isample!=0 && dropZee && dropZmm ? 0.09 : 0);

	  fk->SetLineColor(gk->GetLineColor());
	  gk->Fit(fk, "QRN");
	  tex->SetTextColor(fk->GetLineColor());
	  if (!((dropZee && ss=="zeejet") || (dropZmm && ss=="zmmjet")))
	    tex->DrawLatex(0.55,0.27-0.045*isample+dy,
			   Form("#chi^{2}/NDF = %1.1f / %d",
				fk->GetChisquare(), fk->GetNDF()));
	  tex->SetTextColor(kBlack);

	  // Error band
	  const int n = fk->GetNpar();
	  TMatrixD emat(n,n);
	  gMinuit->mnemat(emat.GetMatrixArray(), n);
	  TF1 *fke = new TF1(Form("fk_%s_%s_%s",cd,cm,cs),
			     sr_fitError, 30, 2640, 1);
	  _sr_fitError_func = fk;
	  _sr_fitError_emat = &emat;

	  fke->SetLineStyle(kSolid);
	  fke->SetLineColor(ss=="multijet" ? kGray+1 : fk->GetLineColor()-10);
	  fke->SetParameter(0,-1);
	  if (!((dropZee && ss=="zeejet") || (dropZmm && ss=="zmmjet")))
	    fke->DrawClone("SAME");
	  fke->SetParameter(0,+1);
	  if (!((dropZee && ss=="zeejet") || (dropZmm && ss=="zmmjet")))
	    fke->DrawClone("SAME");

	  if (!((dropZee && ss=="zeejet") || (dropZmm && ss=="zmmjet"))) {
	    fk->DrawClone("SAME");
	    gk->DrawClone("SAME Pz");
	  }


          ofstream txtFSRDiJet(Form("txt2/GlobalFitOutput_FSRFit_%s_%s_%s.txt",cd, cs,cm),ios_base::app);
          if(etamin==0.&&etamax==0.261)txtFSRDiJet << "{ 1 JetEta 1 JetPt [0]+[1]*log(0.01*x)+[2]*pow(log(0.01*x),2) Correction L2Relative}";
          if(!(etamin==0.&&etamax==1.3))txtFSRDiJet      << Form("\n%7.4f   %7.4f  5 10 6500 %7.4f  %7.4f  %7.4f ", etamin, etamax, fk->GetParameter(0), fk->GetParameter(1), fk->GetParameter(2) );

          

	  // Store soft radiation corrections in fsr subdirectory
	  assert(finout->cd(cd));
	  TDirectory *dout1 = finout->GetDirectory(cd); assert(dout1);
	  //assert(gDirectory->cd(bin)); // broke in ROOT 6.04/08
	  //if (!gDirectory->FindObject("fsr")) gDirectory->mkdir("fsr");
	  //assert(gDirectory->cd("fsr"));
	  assert(dout1->cd(bin));
	  TDirectory *dout2 = dout1->GetDirectory(bin); assert(dout2);
	  if (!dout2->FindObject("fsr")) dout2->mkdir("fsr");
	  assert(dout2->cd("fsr"));
	  TDirectory *dout3 = dout2->GetDirectory("fsr"); assert(dout3);
	  dout3->cd();

	  TH1D *hk = (TH1D*)(isample==0 ? hpt2->Clone() : 
			     (ss=="zjet" ? hpt1b->Clone() : hpt1->Clone()));
	  TProfile *ppt = (isample==0 ? ppt2 :
			   (ss=="zjet" ? ppt1b : ppt1));
          //if (isample==idj) { hk = (TH1D*)hpt4->Clone(); ppt = ppt4; } 
	  if (ss=="dijet") { hk = (TH1D*)hpt4->Clone(); ppt = ppt4; } 
          if (ss=="multijet") { hk = (TH1D*)hpt5->Clone(); ppt = ppt5; } 

	  hk->SetName(Form("hkfsr_%s_%s",cm,cs));
          
          cout << "processing isample " << isample  << " " << samples[isample]<< endl << flush;

	  // for (int i = 1; i != hk->GetNbinsX()+1; ++i) {
	  //   double pt = ppt->GetBinContent(i);
	  //   if (pt>30 && pt<1500) {
	  //     hk->SetBinContent(i, fk->Eval(pt));
          //     if(verbose)cout << "pt " << pt << " fk->eval(pt) " << fk->Eval(pt) << endl << flush;

	  //     hk->SetBinError(i, fabs(fke->Eval(pt)-fk->Eval(pt)));
	  //   }
	  //   else {
          //     if(verbose)cout << "Out of range...: pt " << pt << " fk->eval(pt) " << fk->Eval(pt) << endl << flush;
	  //     hk->SetBinContent(i, 0);
	  //     hk->SetBinError(i, 0);
	  //   }
	  // }
	  
	  // Set non-zero values only for points in graphs
	  hk->Reset();
	  for (int i = 0; i != gk->GetN(); ++i) {

	    double pt = gk->GetX()[i];
	    int ipt = hk->FindBin(pt);
	    pt = hk->GetBinCenter(ipt); // IN UL2017 to reproduce EOY2017
	    hk->SetBinContent(ipt, fk->Eval(pt));
	    hk->SetBinError(ipt, fabs(fke->Eval(pt)-fk->Eval(pt)));
	  }

	  hk->Write(hk->GetName(), TObject::kOverwrite);

	  // Factorize error matrix into eigenvectors
	  // Remember: A = Q*Lambda*Q^-1, where
	  // A is emat, Q is eigmat, and Lambda is a diagonal matrix with
	  // eigenvalues from eigvec on the diagonal. For eigenmatrix
	  // Q^-1 = Q^T, i.e. inverse matrix is the original transposed
	  TVectorD eigvec(n);
	  TMatrixD eigmat = emat.EigenVectors(eigvec);

	  // Eigenvectors are the columns and sum of eigenvectors squared
	  // equals original uncertainty. Calculate histograms from the
	  // eigenvectors and store them
	  TF1 *fkeig = (TF1*)fk->Clone(Form("%s_eig",fk->GetName()));
	  fkeig->SetLineStyle(kDotted);
	  for (int ieig = 0; ieig != n; ++ieig) {

	    // Eigenvector functions
	    for (int i = 0; i != n; ++i) {
	      fkeig->SetParameter(i, fk->GetParameter(i)
				  + eigmat[i][ieig] * sqrt(eigvec[ieig]));
	    }
	    // Draw eigenvectors (plots gets busy, so turned off)
	    // fkeig->DrawClone("SAMEL");

	    // Eigenvector histograms evaluated at bin mean pT
	    TH1D *hke = (TH1D*)hk->Clone(Form("%s_eig%d",hk->GetName(),ieig));
	    hke->Reset();

	    for (int i = 0; i != gk->GetN(); ++i) {

	      double pt = gk->GetX()[i];
	      int ipt = hke->FindBin(pt);
	      // Need to store central value as well, because
	      // uncertainty sources are signed
	      hke->SetBinContent(ipt, fkeig->Eval(pt)-fk->Eval(pt));
	      hke->SetBinError(ipt, fabs(fkeig->Eval(pt)-fk->Eval(pt)));
	    }
	    hke->Write(hke->GetName(), TObject::kOverwrite);
	  }

	  cout << "." << flush;
	} // if tree
      } // for isample
    } // for imethod
  } // for idir
  
  c3->cd(0);
  //cmsPrel(_lumi, true);
  CMS_lumi(c3, 2, 33);
  c3->SaveAs(Form("pdf/%s/softrad_2x6_kfsr_eta%02.0f-%02.0f.pdf",cep,10*etamin,10*etamax));

  ///////////////////////////////////////////////////////////
  // Recreate MPF bias and uncertainties for data/MC ratio:
  // 0) Central value: pTbal factored into jetn and uncl
  //    times expected data/MC response difference for these
  //    plus jetn and uncl data/MC difference times MC response
  // 1) Source 1: Uncl(MC) x DeltaRuncl(data/MC)
  // 2) Source 2: JetN(MC) x DeltaRjetn(data/MC) 
  // 3) Source 3: DeltaUncl(data/MC) x Runlc(MC)
  // 4) Source 4: DeltaJetN(data/MC) x Rjetn(MC)
  // Cross-terms omitted for now
  ////////////////////////////////////////

  finout->Close();
  //fout->Close();
  curdir->cd();
} // softrad

Double_t sr_fitError(Double_t *xx, Double_t *p) {

  assert(_sr_fitError_func);
  assert(_sr_fitError_emat);
  double x = *xx;
  double k = p[0];
  TF1 *f = _sr_fitError_func;
  int n = f->GetNpar();
  TMatrixD &emat = (*_sr_fitError_emat);
  assert(emat.GetNrows()==n);
  assert(emat.GetNcols()==n);
  
  vector<double> df(n);
  for (int i = 0; i != n; ++i) {

    double p = f->GetParameter(i);
    double dp = 0.1*f->GetParError(i);
    f->SetParameter(i, p+dp);
    double fup = f->Eval(x);
    f->SetParameter(i, p-dp);
    double fdw = f->Eval(x);
    f->SetParameter(i, p);
    df[i] = (dp ? (fup - fdw) / (2.*dp) : 0);
  }

  double sumerr2(0);
  for (int i = 0; i != n; ++i) {
    for (int j = 0; j != n; ++j) {
      sumerr2 += emat[i][j]*df[i]*df[j];
    }
  }

  double err = sqrt(sumerr2);

  return (f->Eval(x) + k*err);
}
