// Purpose: Perform Run 2 Legacy global fit
//
// Pre-requisites:
// - softrad3.C : produce FSR+ISR corrections for HDM (MPF+DB)
// - globalFitSyst.C : produce uncertainty shapes
// - globalFitSettings.h : inputs definitions and configurable settings
// Post-processing:
// - globalFitPulls.C : plot pull distributions
// - minitools/createL2L3ResTextFile.C : produce simplified (L2)L3Res text file
// [- minitools/mergerL2L3ResTextFiles.C : combine L3Res with L2Res]
//
// Authors: Mikko Voutilainen, Henning Kirschenmann
//
// Notes: enable systematic source offsetting?
#include "TFile.h"
#include "TF1.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TMinuit.h"
#include "TFitter.h"
#include "TLine.h"
#include "TProfile.h"

#include <string>
#include <vector>
#include <map>
#include <iostream>

#include "tools.C"
#include "tdrstyle_mod15.C"
#include "globalFitSettings.h"

//using namespace globalFitRun2;
using namespace std;

// Helper functions to draw fit uncertainty band for arbitrary TF1 
Double_t fitError(Double_t *x, Double_t *p);
Double_t jesFit(Double_t *x, Double_t *p);
void jesFitter(Int_t& npar, Double_t* grad, Double_t& chi2, Double_t* par,
	       Int_t flag);
void cleanGraph(TGraphErrors *g);
void globalFitEtaBin(double etamin, double etamax);

// Define global variables used in fitError
TF1 *_fitError_func(0);                 // fitError uncertainty function and
TMatrixD *_fitError_emat(0);            // error matrix

// Define global variables used in jesFit and jesFitter
// The structs are defined in globalFitSettings.h
vector<fitData> _vdt;                   // jesFitter input data,
map<string, vector<fitSyst> > _msrc;    // data->sources,
map<string, vector<fitShape> > _mshape; // obs->shapes,
string _obs;                            // data type switch,
TH1D *_hjesref(0);                      // and reference JES
int cnt(0), Nk(0);
TF1 *_jesFit(0);                        // JES fit used in jesFitter


// Call global fit for each eta bin separately
void globalFitRun2() {

  globalFitEtaBin(0.0, 1.3);
} // globalFitRun2

void globalFitEtaBin(double etamin, double etamax) {

  // Set fancy plotting style (CMS TDR style)
  setTDRStyle();

  // Keep track of current working directory so plots don't disappear
  TDirectory *curdir = gDirectory;
  

  // 1. Load input data sets
  ///////////////////////////

  // Open file created by minitools/runAllIOVs.py and recombine.C
  TFile *f = new TFile("rootfiles/jecdataRun2Test.root","READ");
  assert(f && !f->IsZombie());

  // Prepare links to all relevant subdirectories
  f->cd("ratio");
  TDirectory *dratio = gDirectory;
  dratio->cd(Form("eta%02.0f-%02.0f", 10.*etamin, 10*etamax));
  TDirectory *deta = gDirectory;
  deta->cd(Form("sys"));
  TDirectory *dsys = gDirectory;
  deta->cd(Form("fsr"));
  TDirectory *dfsr = gDirectory;
  
  // Load reference JES and uncertainty
  _hjesref = (TH1D*)deta->Get("herr_l2l3res"); assert(_hjesref);
  TH1D *herr = (TH1D*)deta->Get("herr"); assert(herr);

  // Set whitelist for quickly selecting only subset of datasets
  set<string> whitelist;
  for (unsigned int i = 0; i != _gf_datasets_whitelist.size(); ++i) {
    if (_gf_datasets_whitelist[i]!="")
      whitelist.insert(_gf_datasets_whitelist[i]);
  }

  // Create listing of all active datasets
  set<string> datasets;
  for (unsigned int i = 0; i != _gf_datasets.size(); ++i) {

    // Read in information from globalFitSettings.h
    const char *name  = _gf_datasets[i][0].c_str();
    const char *type  = _gf_datasets[i][1].c_str();
    const char *name2 = _gf_datasets[i][2].c_str();

    // Check if input dataset is whitelisted
    if (!whitelist.empty() && whitelist.find(name)==whitelist.end()) continue;
    if (string(name)=="") continue; // missing elements in dataset array
    
    // Retrieve graph for input data
    TGraphErrors *g = (TGraphErrors*)deta->Get(name);
    if (!g) cout << "Input " << name << " not found." << endl << flush;
    assert(g);

    // Patch HDM input from TH1D to graph [temporary] => fix in softrad3.C
    if (TString(name).Contains("hdm_") && !TString(name).Contains("hadw")) {
      TH1D *h = (TH1D*)deta->Get(name);
      assert(h);
      g = new TGraphErrors(h);
      cleanGraph(g);
    }

    // Create fitData object
    fitData data;
    data.name = name;
    data.type = type;
    data.input = (TGraphErrors*)g->Clone(Form("%s_in",name));
    data.output = (TGraphErrors*)g->Clone(Form("%s_out",name));

    // Multijet special
    if (string(name2)!="") {
      TGraphErrors *g2 = (TGraphErrors*)deta->Get(name2);
      assert(g2);
      
      data.input2 = (TGraphErrors*)g2->Clone(Form("%s_in2",name2));
      data.output2 = (TGraphErrors*)g2->Clone(Form("%s_out2",name2));
    }

    // Save fitData for list and vector
    datasets.insert(name);
    _vdt.push_back(data);
  } // for i in datasets


  // 2. Load input uncertainty sources
  ////////////////////////////////////

  // Create listing and mapping of active uncertainty sources
  set<string> sources;
  map<string, int> msrcidx;
  for (unsigned int i = 0; i != _gf_sources.size(); ++i) {

    // Read in information from globalFitSettings.h
    const char *name      = _gf_sources[i][0].c_str();
    const char *appliesTo = _gf_sources[i][1].c_str();
    const char *histname  = _gf_sources[i][2].c_str();

    // Activate only sources that apply to input data sets
    if (datasets.find(appliesTo)==datasets.end()) continue;

    // Use running indexing for active sources
    if (sources.find(name)==sources.end()) {
      msrcidx[name] = sources.size();
    }

    // Retrieve histogram for shape
    TH1D *h = (TH1D*)dsys->Get(histname);
    if (!h) h = (TH1D*)dfsr->Get(histname);
    assert(h);
  
    // Create fitSyst object
    fitSyst syst;
    syst.idx = msrcidx[name];
    syst.name = name;
    syst.appliesTo = appliesTo;
    syst.hist = (TH1D*)h->Clone(Form("%s_sys",histname));

    // Save fitSyst to list and map
    sources.insert(name);
    _msrc[syst.appliesTo].push_back(syst);
  } // for i in sources  


  // 3. Load input fit shapes
  ///////////////////////////

  // Create listing and mapping of active response/composition shapes
  set<string> shapes;
  map<string, int> mshapeidx;
  for (unsigned int i = 0; i != _gf_shapes.size(); ++i) {
    
    string name       = _gf_shapes[i][0];
    string appliesTo  = _gf_shapes[i][1];
    string funcstring = _gf_shapes[i][2];

    if (name=="") continue; // ignore exta empty (commented out) elements

    // Use running indexing for active shapes
    if (shapes.find(name)==shapes.end())
      mshapeidx[name] = shapes.size();

    // Create fitShape object
    fitShape shape;
    shape.idx = mshapeidx[name];
    shape.name = name;
    shape.appliesTo = appliesTo;
    shape.func = new TF1(Form("f1_%s_%s",name.c_str(),appliesTo.c_str()),
			 funcstring.c_str(),15.,6500.);

    // Save fitShapes to list and map
    shapes.insert(name);
    _mshape[appliesTo].push_back(shape);
  } // for i in shapes

  // Create function to plot for JES (or composition)
  double minpt = 15.;
  double maxpt = 2000;//1500.;
  int njesFit = shapes.size();
  _jesFit = new TF1("jesFit",jesFit,minpt,maxpt,njesFit);


  // 4. Perform chi2 fit
  ///////////////////////

  //const int npar = _jesFit->GetNpar();
  const int npar = shapes.size();
  const int nsrc = sources.size();//_msrc->size();
  Int_t ntot = npar+nsrc;

  cout << endl;
  cout << "Global fit has " << npar << " fit parameters and "
       << nsrc << " nuisance parameters." << endl;
  if (_gf_penalizeFitPars)
    cout << "Fit parameters have Gaussian prior" << endl;
  cout << endl;

  //for (int i = 0; i != npar; ++i) a[i] = _jesFit->GetParameter(i);

  // Setup global chi2 fit (jesFitter is our function)
  TFitter *fitter = new TFitter(ntot);
  fitter->SetFCN(jesFitter);

  // Set parameters
  vector<double> a(ntot, 0);
  vector<string> parnames(ntot); // empty for now, to fill with shapes/sources
  for (int i = 0; i != ntot; ++i)
    fitter->SetParameter(i, parnames[i].c_str(), a[i], (i<npar ? 0.01 : 1),
			 -100, 100);

  // Run fitter (multiple times if needed)
  const int nfit = 1;
  cnt = 0;
  for (int i = 0; i != nfit; ++i)
    fitter->ExecuteCommand("MINI", 0, 0);
  TMatrixD emat(ntot, ntot);
  gMinuit->mnemat(emat.GetMatrixArray(), ntot);

  // Retrieve the chi2 the hard way
  Double_t tmp_par[ntot], tmp_err[ntot];
  TVectorD vpar(ntot);
  TVectorD verr(ntot);
  Double_t chi2_gbl(0), chi2_src(0), chi2_par(0), chi2_data(0);
  //vector<double> vchi2_data(gs2.size(),0);
  //vector<int> vndata(gs2.size(),0);
  int npar_true(0), nsrc_true(0), ndt(0);
  Double_t grad[ntot];
  Int_t flag = 1;
  //TH1D *hsrc = new TH1D("hsrc",";Nuisance parameter;",12,-3,3);

  for (int i = 0; i != ntot; ++i) {
    tmp_par[i] = fitter->GetParameter(i);
    tmp_err[i] = fitter->GetParError(i);
    vpar[i] = fitter->GetParameter(i);
    verr[i] = fitter->GetParError(i);
  }
  jesFitter(ntot, grad, chi2_gbl, tmp_par, flag);

  //cout << "List of fitted sources that count (nonzero):" << endl;
  for (int i = 0; i != ntot; ++i) {
    if (fabs(tmp_par[i])!=0 || fabs(tmp_err[i]-1)>1e-2) {
      if (i < npar) ++npar_true;
      else          ++nsrc_true;
      //hsrc->Fill(tmp_par[i]);
    }
    if (i < npar) chi2_par += pow(tmp_par[i],2);
    else          chi2_src += pow(tmp_par[i],2);
  }
  for (unsigned int i = 0; i != _vdt.size(); ++i) {
    TGraphErrors *gout = _vdt[i].output;
    _obs = _vdt[i].type; // for _jesFit
    for (int j = 0; j != gout->GetN(); ++j) {
      double x = gout->GetX()[j];
      double y = gout->GetY()[j];
      double ey = gout->GetEY()[j];
      chi2_data += pow((y - _jesFit->Eval(x)) / ey, 2);
      //vchi2_data[i] += pow((y-jesfit->Eval(x))/ey,2);
      //++vndata[i];
      ++ndt;
    }
  }
  
  cout << endl;
  cout << Form("Used %d data points, %d fit parameters and %d nuisances.\n"
	       "Data chi2/NDF = %1.1f / %d [%1.0f,%1.0f]\n"
	       "Nuisance chi2/Nsrc = %1.1f / %d\n"
	       "Parameter chi2/Npar = %1.1f / %d\n"
	       "Total chi2/NDF = %1.1f / %d\n",
	       ndt, npar, nsrc_true, chi2_data, ndt - npar,
	       _jesFit->GetXmin(), _jesFit->GetXmax(),
	       chi2_src, nsrc_true,
	       chi2_par, npar_true,
	       chi2_gbl, ndt - npar);
  cout << endl;
  vector<fitShape> &v = _mshape["Rjet"];
  assert(int(v.size())==njesFit);
  for (unsigned int i = 0; i != v.size(); ++i) {
    cout << Form("  %5s : %+5.2f +/- %5.2f", v[i].name.c_str(),
		 vpar[v[i].idx], verr[v[i].idx]) << endl;
		 //_jesFit->GetParameter(v[i].idx),
		 //_jesFit->GetParError(v[i].idx)) << endl;
  } // for i in _mshape

  // 5. Draw results
  //////////////////
  bool drawResults = true;
  if (drawResults) {

    // Graphical settings in globalFitStyles.h
    #include "globalFitStyles.h"
    curdir->cd();

    // Create canvas
    lumi_13TeV = "globalFitRun2.C(\"Run2Test\")";
    TH1D *h = tdrHist("h","JES",0.982+1e-5,1.025-1e-5); // ratio (hdm)
    TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);
    gPad->SetLogx();

    TLine *l = new TLine();
    TLegend *leg = tdrLeg(0.60,0.90,0.80,0.90);
    TLegend *leg2 = tdrLeg(0.45,0.15,0.65,0.30);

    // Draw fit on the back
    _jesFit->SetRange(15.,3500.); // nice range
    _jesFit->SetNpx(3485.); // dense binning for log scale
    _obs = "Rjet";
    TGraph *gr = new TGraph(_jesFit); // use graph to keep '_obs' setting
    TGraphErrors *gre = new TGraphErrors(gr->GetN());
    _fitError_func = _jesFit;
    _fitError_emat = &emat;
    double k = 1;
    for (int i = 0; i != gr->GetN(); ++i) {
      double pt = gr->GetX()[i];
      gre->SetPoint(i, pt, gr->GetY()[i]);
      gre->SetPointError(i, 0., fitError(&pt, &k) - gr->GetY()[i]);
    }
    tdrDraw(herr,"E3",kFullCircle,kCyan+2,kSolid,-1,1001,kCyan+1);
    tdrDraw(gre,"E3",kNone,kBlack,kSolid,-1,1001,kYellow+1);
    tdrDraw(gr,"Lz",kNone,kBlack);
    l->SetLineStyle(kDashed);
    l->DrawLine(15,1,3500,1);

    leg2->AddEntry(l,"Run 2 avg.","L");
    leg2->AddEntry(herr,"Total unc.","F");
    leg2->AddEntry(gre,"Fit unc.","FL");

    // Separate canvas for CHF, NHF, NEF
    TH1D *hc = tdrHist("h","PF composition (0.01)",-2e-2+1e-7,2e-2-1e-7);
    TCanvas *c1c = tdrCanvas("c1c",hc,4,11,kSquare);
    gPad->SetLogx();
    l->DrawLine(15,0,3500,0);

    // Separate canvas for CEF, MUF
    TH1D *hl = tdrHist("hl","PF composition (0.01)",-2e-3+1e-7,2.5e-3-1e-7);
    TCanvas *c1l = tdrCanvas("c1l",hl,4,11,kSquare);
    gPad->SetLogx();
    l->DrawLine(15,0,3500,0);

    // Sanity check PF composition sums
    TGraphErrors *gpfjet(0);
    TGraphErrors *gzjet(0);
    TGraphErrors *gzljet(0);
    TGraphErrors *gzmjet(0);
    TGraphErrors *ggjet(0);

    for (unsigned int i = 0; i != _vdt.size(); ++i) {
      
      string name = _vdt[i].name;
      string type = _vdt[i].type;

      if (type=="Rjet") {
	c1->cd();
	//leg->AddEntry(_vdt[i].input,_gf_label[name],"PLE");
	leg->AddEntry(_vdt[i].output,_gf_label[name],"PLE");
	leg->SetY1NDC(leg->GetY1NDC()-0.05);
      }
      else if (type=="cef" || type=="muf") c1l->cd();
      else c1c->cd();
      
      // Default settings
      if (_gf_color[name]==0)  _gf_color[name] = kBlack;
      if (_gf_marker[name]==0) _gf_marker[name] = kFullCircle;
      if (_gf_label[name]==0)  _gf_label[name] = name.c_str();
      if (_gf_size[name]==0)   _gf_size[name] = 1.0;
      
      //tdrDraw(_vdt[i].input,"Pz",kOpenSquare,mcolor[name]);
      //tdrDraw(_vdt[i].output,"Pz",kFullCircle,mcolor[name]);
      tdrDraw(_vdt[i].output,"Pz",_gf_marker[name],_gf_color[name]);
      if (name=="hdm_mpfchs1_multijet")
	tdrDraw(_vdt[i].input,"Pz",kOpenTriangleUp,_gf_color[name]);

      _vdt[i].output->SetMarkerSize(_gf_size[name]);
      _vdt[i].input->SetMarkerSize(_gf_size[name]);
      
      /*
      // Add up all the fractions to check they sum up to unity
      if (TString(name.c_str()).Contains("pfjet") && type!="Rjet") {
	if (!gpfjet) gpfjet = (TGraphErrors*)_vdt[i].input->Clone("gpfjet");
	else gpfjet = tools::diffGraphs(gpfjet,_vdt[i].input,1,-1);
      }
      if (TString(name.c_str()).Contains("zjet") && type!="Rjet") {
	if (!gzjet) gzjet = (TGraphErrors*)_vdt[i].input->Clone("gzjet");
	else gzjet = tools::diffGraphs(gzjet,_vdt[i].input,1,-1);
      }
      if (TString(name.c_str()).Contains("zlljet") && type!="Rjet") {
	if (!gzljet) gzljet = (TGraphErrors*)_vdt[i].input->Clone("gzljet");
	else gzljet = tools::diffGraphs(gzljet,_vdt[i].input,1,-1);
      }
      if (TString(name.c_str()).Contains("zmmjet") && type!="Rjet") {
	if (!gzmjet) gzmjet = (TGraphErrors*)_vdt[i].input->Clone("gzmjet");
	else gzmjet = tools::diffGraphs(gzmjet,_vdt[i].input,1,-1);
      }
      if (TString(name.c_str()).Contains("gamjet") && type!="Rjet") {
	if (!ggjet) ggjet = (TGraphErrors*)_vdt[i].input->Clone("ggjet");
	else ggjet = tools::diffGraphs(ggjet,_vdt[i].input,1,-1);
      }
      */
    } // for i in _vdt   

    c1l->cd();
    _obs = "cef";
    TGraph *gcef = new TGraph(_jesFit);  // use graph to keep '_obs' setting
    tdrDraw(gcef,"Lz",kNone,kCyan+2,kSolid);

    _obs = "muf";
    TGraph *gmuf = new TGraph(_jesFit);  // use graph to keep '_obs' setting
    tdrDraw(gmuf,"Lz",kNone,kMagenta+2,kSolid);

    c1c->cd();
    /*
    tdrDraw(gzmjet,"PLz",kFullSquare,kGray); gzmjet->SetMarkerSize(0.8);
    tdrDraw(gzljet,"PLz",kFullSquare,kGray);
    tdrDraw(gzjet,"PLz",kFullSquare,kGray+1);
    tdrDraw(ggjet,"PLz",kFullCircle,kGray+2);
    tdrDraw(gpfjet,"PLz",kFullDiamond,kBlack);
    */
    _obs = "chf";
    TGraph *gchf = new TGraph(_jesFit);  // use graph to keep '_obs' setting
    TGraphErrors *gchfe = new TGraphErrors(gchf->GetN());
    for (int i = 0; i != gchf->GetN(); ++i) {
      double pt = gr->GetX()[i];
      gchfe->SetPoint(i, pt, gchf->GetY()[i]);
      gchfe->SetPointError(i, 0., fitError(&pt, &k) - gchf->GetY()[i]);
    }
    tdrDraw(gchfe,"E3",kNone,kRed+1,kSolid,-1,1001,kRed-9);
    tdrDraw(gchf,"Lz",kNone,kRed,kSolid);

    _obs = "nhf";
    TGraph *gnhf = new TGraph(_jesFit);  // use graph to keep '_obs' setting
    TGraphErrors *gnhfe = new TGraphErrors(gnhf->GetN());
    for (int i = 0; i != gnhf->GetN(); ++i) {
      double pt = gr->GetX()[i];
      gnhfe->SetPoint(i, pt, gnhf->GetY()[i]);
      gnhfe->SetPointError(i, 0., fitError(&pt, &k) - gnhf->GetY()[i]);
    }
    tdrDraw(gnhfe,"E3",kNone,kGreen+3,kSolid,-1,1001,kGreen-9);
    tdrDraw(gnhf,"Lz",kNone,kGreen+2,kSolid);

    _obs = "nef";
    TGraph *gnef = new TGraph(_jesFit);  // use graph to keep '_obs' setting
    TGraphErrors *gnefe = new TGraphErrors(gnef->GetN());
    for (int i = 0; i != gnef->GetN(); ++i) {
      double pt = gr->GetX()[i];
      gnefe->SetPoint(i, pt, gnef->GetY()[i]);
      gnefe->SetPointError(i, 0., fitError(&pt, &k) - gnef->GetY()[i]);
    }
    tdrDraw(gnefe,"E3",kNone,kBlue+1,kSolid,-1,1001,kBlue-9);
    tdrDraw(gnef,"Lz",kNone,kBlue,kSolid);

    _obs = "Rjet";
    TGraph *grjt = new TGraph(_jesFit);  // use graph to keep '_obs' setting
    tdrDraw(grjt,"Lz",kNone,kBlack,kSolid);


    // test case: gamma+jet true response in MC
    if (false) {
      c1->cd();
      TFile *fg = new TFile("../gamjet/files/GamHistosFill_mc_2018P8.root",
			    "READ");
      assert(fg && !fg->IsZombie());
      TProfile *p = (TProfile*)fg->Get("control/prgen"); assert(p);
      //tdrDraw(p,"HIST",kNone,kCyan+2,kSolid,-1,kNone);
      TProfile *pr = (TProfile*)fg->Get("control/prjet"); assert(pr);
      TProfile *pg = (TProfile*)fg->Get("control/pgjet"); assert(pg);
      TH1D *hp = pr->ProjectionX("hp");
      hp->Divide(pg);
      tdrDraw(hp,"HIST",kNone,kCyan+2,kSolid,-1,kNone);
      //tdrDraw(p,"HIST",kNone,kCyan+3,kDashed,-1,kNone);
    }
    // test case: Z+jet true response in MC
    if (false) {
      //TFile *fz = new TFile("rootfiles/jme_bplusZ_merged_v35_2018_emu_wTTJets.root","READ");
      TFile *fz = new TFile("rootfiles/jme_bplusZ_merged_v37_2016FH_mu.root","READ");
      assert(fz && !fz->IsZombie());
      // rz_zmmjet_a100 : (Zpt,genZpt/Zpt) => Z resp. log-lin slope +1% to -2%
      // rpt_zmmjet_a100 : (Zpt,RpT) => pTreco/pTZ vs pTZ
      // rbal_zmmjet_a100 : (Zpt,ljet_pt/Zpt) => pTrgen/pTZ vs pTZ
      // *rgenjet1_zmmjet_a100 : (Zpt,...) => pTgen*cos(dphi)/pTZ vs pTZ
      // *rmpfjet1_zmmjet_a100 : (Zpt,RMPFjet1) => pTreco*cos(dphi)/pTZ vs pTZ
      TGraphErrors *gr2 = (TGraphErrors*)fz->Get("mc/eta_00_13/rmpfjet1_zmmjet_a100"); assert(gr2);
      TGraphErrors *gg2 = (TGraphErrors*)fz->Get("mc/eta_00_13/rgenjet1_zmmjet_a100"); assert(gg2);
      TGraphErrors *g2 = tools::ratioGraphs(gr2,gg2);
      TGraphErrors *gr = (TGraphErrors*)fz->Get("mc/eta_00_13/rpt_zmmjet_a100"); assert(gr);
      TGraphErrors *gg = (TGraphErrors*)fz->Get("mc/eta_00_13/rbal_zmmjet_a100"); assert(gg);
      TGraphErrors *g = tools::ratioGraphs(gr,gg);
      tdrDraw(g,"Lz",kNone,kRed+2,kSolid,-1,kNone);
      tdrDraw(g2,"Lz",kNone,kRed+2,kDashed,-1,kNone);
      // rz may be relevant, HDM is above MC truth at high pT
      TGraphErrors *gz = (TGraphErrors*)fz->Get("mc/eta_00_13/rz_zmmjet_a100");
      assert(gz);
      TGraphErrors *gz1 = tools::ratioGraphs(g,gz);
      tdrDraw(gz1,"Lz",kNone,kRed+2,kDotted,-1,kNone);
      // Much too large fix at high pT. Different in 2016FH? Something else?
    }

    c1->cd(); gPad->RedrawAxis();
    c1c->cd(); gPad->RedrawAxis();
    c1l->cd(); gPad->RedrawAxis();
   
    c1->SaveAs("pdf/globalFitRun2/Run2Test_rjet.pdf");
    c1c->SaveAs("pdf/globalFitRun2/Run2Test_pf.pdf");
    c1l->SaveAs("pdf/globalFitRun2/Run2Test_mu.pdf");
  }
} // globalFitEtaBin


// Generic fit function for JES and composition
Double_t jesFit(Double_t *x, Double_t *p) {

  // Choose JES (var~1) or PF composition (var~0) as baseline
  double var = (_obs=="Rjet" ? 1. : 0.);
  double pt = x[0];

  // Permit baseline shift for JES with L2L3Res inputs. Needed for PFcomp fit
  double jesref = 1;
  if (_gf_useJESref && _obs=="Rjet") {
    assert(_hjesref);
    jesref = _hjesref->Interpolate(pt);
  }

  // Load available shapes for this observable
  vector<fitShape> &v = _mshape[_obs];
  for (unsigned int i = 0; i != v.size(); ++i) {

    // Calculate variation and add it to total
    int idx = v[i].idx;   assert(idx>=0);
    TF1 *f1 = v[i].func;  assert(f1);
    var += p[idx] * f1->Eval(pt) * 0.01; // fullSimShapes in %'s
  } // for i in mshape

  return (var / jesref);
} // jesFit


// General chi2 fit function
void jesFitter(Int_t& npar, Double_t* grad, Double_t& chi2, Double_t* par,
	       Int_t flag) {

  // Basic checks
  assert(_jesFit);

  // Parametes for nuisances (sources)
  Double_t *ps = &par[_jesFit->GetNpar()];
  int ns = npar - _jesFit->GetNpar();

  if (flag) {
    
    // do the following calculation:
    // chi2 = sum_i( (x_i+sum_s(a_s y_si) -fit)^2/sigma_i^2) + sum_s(a_s^2)
    chi2 = 0;
    Nk = 0;

    // Loop over input data (graphs x bins)
    // - for each point, add up source eigenvectors x nuisance parameters
    // - then calculate chi2 adding up residuals + nuisance parameters
    for (unsigned int ig = 0; ig != _vdt.size(); ++ig) {

      string name         = _vdt[ig].name;
      string type         = _vdt[ig].type;
      TGraphErrors *gin   = _vdt[ig].input;  assert(gin);
      TGraphErrors *gin2  = _vdt[ig].input2;
      TGraphErrors *gout  = _vdt[ig].output; assert(gout);
      TGraphErrors *gout2 = _vdt[ig].output2;
      
      for (int i = 0; i != gin->GetN(); ++i) {

	// Retrieve central value and uncertainty for this point
	double pt = gin->GetX()[i];
	double data = gin->GetY()[i];
	double sigma = gin->GetEY()[i];

	//if (pt < ptminJesFit) continue;

	// Calculate fit value at this point
	for (int ipar = 0; ipar != _jesFit->GetNpar(); ++ipar) {
	  _jesFit->SetParameter(ipar, par[ipar]);
	}
	_obs = type;
	double fit = _jesFit->EvalPar(&pt,par);

	// For multijet balancing, multiply data by reference JES
	if (TString(name.c_str()).Contains("multijet")) {
	  assert(gin2);
	  double ptref = gin2->GetY()[i] * pt;
	  double fitRef = _jesFit->EvalPar(&ptref,par);
	  data *= fitRef;
	  sigma *= fitRef;
	} // multijet

	// Calculate total shift caused by all nuisance parameters
	double shifts = 0;
	vector<fitSyst> &_vsrc = _msrc[name];
	for (unsigned int is = 0; is != _vsrc.size(); ++is) {

	  TH1D *hsrc = _vsrc[is].hist; assert(hsrc);
	  int idx = _vsrc[is].idx;

	  int ipt = hsrc->FindBin(pt);
	  shifts += ps[idx] * hsrc->GetBinContent(ipt);
	}
	
	// Add chi2 from residual
	double chi = (data + shifts - fit) / sigma;
	chi2 += chi * chi;
	++Nk;

	// Store shifted data
	assert(gout->GetN()==gin->GetN() && gout->GetX()[i]==pt);
	gout->SetPoint(i, pt, data + shifts);

	// For multijets, store also downward extrapolation
	if (TString(name.c_str()).Contains("multijet")) {
	  double jes = _jesFit->EvalPar(&pt, par);
	  double ptref = pt * gin2->GetY()[i];
	  double jesref = _jesFit->EvalPar(&ptref, par);
	  // MJB = jes / jesref
	  // data = MJB*jesref => "jesref" = jes / MJB = jes * jesref/data
	  // NB: double check the logic here
	  gout2->SetPoint(i, ptref, jes * jesref / (data + shifts));
	  gout2->SetPointError(i, gout->GetEX()[i], gout2->GetEY()[i]);
	} // multijet
      } // for ipt
    } // for ig
  
    // Add chi2 from nuisance parameters
    for (int is = 0; is != ns; ++is) {
      
      chi2 += ps[is]*ps[is];
    } // for ipar

    // Add penalty for fit parameters (Bayesian prior, essentially)
    if (_gf_penalizeFitPars) {
      for (int ipar = 0; ipar != _jesFit->GetNpar(); ++ipar) {
	chi2 += par[ipar] * par[ipar];
	++Nk;
      }
    } // penalizeFitPars
    
    // Give some feedback on progress in case loop gets stuck
    if ((++cnt)%100==0) cout << "." << flush;
  } // if flag
  else {
    if (grad) {}; // suppress warning;
    return;
  }

} // jesFitter


// Calculate fit uncertainty for an arbitrary function                          
// using standard error propagation with differentials                          
Double_t fitError(Double_t *xx, Double_t *pp) {

  // Sanity checks on inputs                                                    
  assert(_fitError_func);
  assert(_fitError_emat);
  double x = *xx;
  double k = pp[0];
  TF1 *f = _fitError_func;
  int n = f->GetNpar();
  TMatrixD &emat = (*_fitError_emat);
  assert(emat.GetNrows()==n);
  assert(emat.GetNcols()==n);

  // Partial derivatives as differentials with 10% step size                    
  vector<double> df(n);
  for (int i = 0; i != n; ++i) {

    double p = f->GetParameter(i);
    double dp = 0.1*sqrt(emat[i][i]); // step size 10% of uncertainty           
    f->SetParameter(i, p+dp);
    double fup = f->Eval(x);
    f->SetParameter(i, p-dp);
    double fdw = f->Eval(x);
    f->SetParameter(i, p);

    // Calculate partial derivative as a simple differential                    
    df[i] = (dp ? (fup - fdw) / (2.*dp) : 0);
  }

  // Perform standard error propagation                                         
  double sumerr2(0);
  for (int i = 0; i != n; ++i) {
    for (int j = 0; j != n; ++j) {
      sumerr2 += emat[i][j]*df[i]*df[j];
    }
  }

  double err = sqrt(sumerr2);

  return (f->Eval(x) + k*err);
}

void cleanGraph(TGraphErrors *g) {
  for (int i = g->GetN()-1; i != -1; --i) {
    if (g->GetY()[i]==0 && g->GetEY()[i]==0)
      g->RemovePoint(i);
  }
} // cleanGraph
