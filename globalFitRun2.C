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

#include <string>
#include <vector>
#include <map>
#include <iostream>

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
string _obs;                            // and data type switch
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
  
  set<string> whitelist;
  for (unsigned int i = 0; i != _gf_datasets_whitelist.size(); ++i) {
    whitelist.insert(_gf_datasets_whitelist[i]);
  }

  // Create listing of all active datasets
  set<string> datasets;
  for (unsigned int i = 0; i != _gf_datasets.size(); ++i) {

    // Read in information from globalFitSettings.h
    const char *name = _gf_datasets[i][0].c_str();
    const char *type = _gf_datasets[i][1].c_str();

    // Check if input data is whitelisted
    if (!whitelist.empty() && whitelist.find(name)==whitelist.end()) continue;

    // Retrieve graph for input data
    TGraphErrors *g = (TGraphErrors*)deta->Get(name);
    assert(g);

    // Patch HDM input from TH1D to graph [temporary] => fix in softrad3.C
    if (TString(name).Contains("hdm_")) {
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
    data.input2 = 0;
    data.output = (TGraphErrors*)g->Clone(Form("%s_out",name));
    data.output2 = 0;

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
  double maxpt = 1500.;
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

  vector<double> a(ntot, 0);
  //for (int i = 0; i != npar; ++i) a[i] = _jesFit->GetParameter(i);

  // Setup global chi2 fit (jesFitter is our function)
  TFitter *fitter = new TFitter(ntot);
  fitter->SetFCN(jesFitter);

  // Set parameters
  vector<string> parnames(ntot);
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
  Double_t chi2_gbl(0), chi2_src(0), chi2_data(0);
  //vector<double> vchi2_data(gs2.size(),0);
  //vector<int> vndata(gs2.size(),0);
  int nsrc_true(0), ndt(0);
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
  for (int i = npar; i != ntot; ++i) {
    if (fabs(tmp_par[i])!=0 || fabs(tmp_err[i]-1)>1e-2) {
      ++nsrc_true;
      //hsrc->Fill(tmp_par[i]);
    }
    chi2_src += pow(tmp_par[i],2);
  }
  for (unsigned int i = 0; i != _vdt.size(); ++i) {
    TGraphErrors *gout = _vdt[i].output;
    _obs = _vdt[i].type; // for jesfit
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
	       "Data chi2/NDF = %1.1f / %d\n"
	       "Nuisance chi2/Nsrc = %1.1f / %d\n"
	       "Total chi2/NDF = %1.1f / %d\n",
	       ndt, npar, nsrc_true, chi2_data, ndt - npar, chi2_src, nsrc_true,
	       chi2_gbl, ndt - npar);
  cout << endl;


  // 5. Draw results
  //////////////////
  bool drawResults = true;
  if (drawResults) {

    // Graphical settings
    map<string, int> mcolor;
    mcolor["hdm_mpfchs1_gamjet"] = kBlue;
    mcolor["hdm_mpfchs1_zjet"]   = kRed;


    curdir->cd();
    TH1D *h = tdrHist("h","JES",0.95,1.05);
    lumi_13TeV = "Run2Test";
    TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);
    gPad->SetLogx();

    for (int i = 0; i != _vdt.size(); ++i) {
      
      string name = _vdt[i].name;
      tdrDraw(_vdt[i].input,"Pz",kOpenSquare,mcolor[name]);
      tdrDraw(_vdt[i].output,"Pz",kFullCircle,mcolor[name]);
    }      

    _obs = "Rjet";
    _jesFit->Draw("SAME");
  }
} // globalFitEtaBin


// Generic fit function for JES and composition
Double_t jesFit(Double_t *x, Double_t *p) {

  // Choose JES (var~1) or PF composition (var~0) as baseline
  double var = (_obs=="Rjet" ? 1. : 0.);
  double pt = x[0];

  // Load available shapes for this observable
  vector<fitShape> &v = _mshape[_obs];
  for (unsigned int i = 0; i != v.size(); ++i) {

    // Calculate variation and add it to total
    int idx = v[i].idx;   assert(idx>=0);
    TF1 *f1 = v[i].func;  assert(f1);
    var += p[idx] * f1->Eval(pt);
  } // for i in mshape

  return var;
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
