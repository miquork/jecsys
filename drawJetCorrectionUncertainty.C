// Adapted from D0 Experiment jetcorr/macros/RjetUncertainty.C
//
// Draw the uncertainty plots for L3 absolute response correction
// The uncertainties are stored in L3Corr.cpp/hpp

#include "TCanvas.h"
#include "TGraph.h"
#include "TH1.h"
#include "TLine.h"
#include "TPave.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TMath.h"

#include "settings.h" // _lumi, _pdf, _eps
#include "JECUncertainty.hpp"
// => Load .cpp in mk_drawJetCorrectionUncertainty.C 

//#include "tdrstyle_mod.C"
#include "tdrstyle_mod14.C"

#include <string>
#include <fstream>
#include <map>

using namespace std;

// Don't plot individual bins, just keep 4x2
bool _minimal = false;//true;
bool _fourbytwo = false;
bool _twobythree = false;//true; // for paper
bool _extra = false; // single pion plots
bool _paper = true; //  for paper

// Plot uncertainty (true) or source (false)
bool _absUncert = true;//true;//false
// NB: All source files are currently printed together with AK4PFchs uncertainty
bool _doTXT = true; // create uncertainty and source text files
bool _didTXT = false;

// List of (hard-coded) default parameters
jec::JetAlgo  d_algo = jec::AK4PFchs; // Replaced in function call
const jec::DataType d_type = jec::DATA; // Uncertainties for data (or data/MC)
//const double d_mu = 24.68; // 80XV8 RunG eta0.0-1.3 jet450
// from https://hypernews.cern.ch/HyperNews/CMS/get/AUX/2018/02/18/13:45:27-50622-pileup-2018-02-16.gif :
const double d_mu = 32.5; // 2017 golden 
const bool d_mpf = true; // L2L3Res uncertainties for MPF method
map<jec::JetAlgo,const char*> *_algnames;

// global histogram counter to avoid memory leaks from duplicate names
int icnt(0);

TCanvas *_canvas(0); // global canvas
int _icanvas = 0;

const int kBeige = 41;
const int kDarkGray = kGray+1;
const int kRightHatch = 3654;
const int kLeftHatch = 3645;
const int kVertHatch = 3699;
const int kHorixHatch = 3600;
const int kCrossHatch = 3644;

// Helper routine to do much the same as TLegend
void DrawLine(TGraph *g, double xmin, double ymin, const char *str) {

  TLine* l1 = new TLine();
  l1->SetLineStyle(g->GetLineStyle());
  l1->SetLineColor(g->GetLineColor());
  l1->SetLineWidth(g->GetLineWidth());
  l1->DrawLineNDC(xmin, ymin+0.01, xmin+0.07, ymin+0.01);
  TLatex* tStat = new TLatex(xmin+0.09, ymin, str);
  tStat->SetNDC();
  tStat->SetTextSize(0.045);
  tStat->Draw();
}
void DrawLine(TH1 *h, double xmin, double ymin, const char *str) {

  TLine* l1 = new TLine();
  l1->SetLineStyle(h->GetLineStyle());
  l1->SetLineColor(h->GetLineColor());
  l1->SetLineWidth(h->GetLineWidth());
  l1->DrawLineNDC(xmin, ymin+0.01, xmin+0.07, ymin+0.01);
  
  if (string(str)!="") {
    TLatex* tStat = new TLatex(xmin+0.09, ymin, str);
    tStat->SetNDC();
    tStat->SetTextSize(0.045);
    tStat->Draw();
  }
}
void DrawFill(TH1 *h, double xmin, double ymin, const char *str) {

  TPave *b1 = new TPave(xmin, ymin+0.01-0.025, xmin+0.07, ymin+0.01+0.025,
			1, "NDC");
  b1->SetFillStyle(h->GetFillStyle());
  b1->SetFillColor(h->GetFillColor());
  b1->SetLineStyle(1);
  b1->SetLineColor(h->GetLineColor());
  b1->Draw("L");
  
  if (string(str)!="") {
    TLatex* tTot = new TLatex(xmin+0.09, ymin, str);
    tTot->SetNDC();
    tTot->SetTextSize(0.045);
    tTot->Draw();
  }
}

struct uncert {
  string name;
  string title;
  jec::ErrorTypes type;
  string method; // mpf, pt or default
  string jetalgo; // ak5pf, ak5calo, ak5pf etc.
  double mu; // average <mu>
  int lcolor;
  int lstyle;
  int lwid;
  int mcolor;
  int mstyle;
  int fcolor;
  int fstyle;
  string option;
  //const char* opt; // fails in ROOT6
  uncert(string name_, string title_, jec::ErrorTypes type_,
	 string method_, string jetalgo_, double mu_,
	 int lcolor_, int lstyle_, int lwid_ = 1,
	 int mcolor_ = kBlack, int mstyle_ = kNone,
	 int fcolor_ = kBlack, int fstyle_ = kNone,
	 string option_ = "L") :
    name(name_), title(title_), type(type_),
    method(method_), jetalgo(jetalgo_), mu(mu_),
    lcolor(lcolor_), lstyle(lstyle_), lwid(lwid_),
    mcolor(mcolor_), mstyle(mstyle_),
    fcolor(fcolor_), fstyle(fstyle_),
    option(option_) {

    //opt = option.c_str();
    assert(method=="default" || method=="mpf" || method=="pt");
    assert(mu>=0 || mu==-1);
  }
};

void plotUncertainty(vector<uncert> const& sys,
		     unsigned int nsys1, unsigned int nsys2,
		     jec::JetAlgo jetAlg,
		     string name, string ytitle,
		     string label1, string label2,
		     double emax, double ptmin,//);//, bool plotLog);
		     string type="fixPt", double typevar=0.);

void drawJetCorrectionUncertainty(string algo = "AK4PFchs",
				  bool doTXT = _doTXT,
				  bool minimal = _minimal) {

  writeExtraText = true;//false; // for JEC paper CWR
  
  _doTXT = doTXT;
  _didTXT = false;
  _minimal = minimal;
  if (algo=="AK4PF") d_algo = jec::AK4PF;
  if (algo=="AK4PFchs") d_algo = jec::AK4PFchs;
  if (algo=="AK7PF") d_algo = jec::AK7PF;
  if (algo=="AK7PFchs") d_algo = jec::AK7PFchs;
  if (algo=="AK8PF") d_algo = jec::AK8PF;
  if (algo=="AK8PFchs") d_algo = jec::AK8PFchs;
  if (algo=="AK4CALO") d_algo = jec::AK4CALO;
  if (algo=="AK7CALO") d_algo = jec::AK7CALO;

  cout << "drawJetCorrectionUncertainty" << endl << flush;

  vector<uncert> sy;
  sy.push_back(uncert("tot", "Total uncertainty", jec::kData,
		      "default", "default", -1, // defaults
		      kBlack, kSolid, 1, // line
		      kBlack, kNone, // marker
		      kDarkGray, 1001, "LF")); // fill
  sy.push_back(uncert("notime", "Excl. flavor, time", jec::kDataNoFlavorNoTime,
		      "default", "default", -1, // defaults
		      kOrange+2, kSolid, 1, // line
		      kBlack, kNone, // marker
  		      kOrange, 1001, "LF")); // fill
  sy.push_back(uncert("runI", "Run I", jec::kRunI,
		      "default", "default", -1, // defaults
		      kRed-9, kSolid, 1, // line
		      kBlack, kNone, // marker
  		      kRed-9, 1001, "LF")); // fill
  sy.push_back(uncert("absolute", "Absolute scale", jec::kAbsolute,
  		      "default", "default", -1, // defaults
  		      kRed, kSolid, 1, // line
 		      kRed, kOpenCircle, // marker
  		      kNone, kNone, "LP")); // fill
  sy.push_back(uncert("relative", "Relative scale",
		      jec::kRelative & ~jec::kRelativeBal & ~jec::kRelativeSample,
		      "default", "default", -1, // defaults
		      kBlack, kSolid, 1, // line
		      kBlack, kFullTriangleDown, // marker
		      kNone, kNone, "LP")); // fill
  sy.push_back(uncert("pileup", Form("Pileup (#LT#mu#GT=%1.0f)",d_mu), jec::kPileUp,
		      "default", "default", -1, // defaults
		      kBlue, kNone, 1, // line
		      kBlue, kFullSquare, // marker
		      kNone, kNone, "LP")); // fill

  sy.push_back(uncert("relativebugs", "Method & sample",jec::kRelativeBal|jec::kRelativeSample,
		      "default", "default", -1, // defaults
		      kBlack, kSolid, 1, // line
		      kBlack, kOpenTriangleDown, // marker
		      kNone, kNone, "LP")); // fill

  sy.push_back(uncert("flavjec", "Jet flavor (QCD)", jec::kFlavorQCD,
		      "default", "default", -1, // defaults
		      kGreen+2, kSolid, 1, // line
		      kGreen+2, kOpenSquare, // marker
		      kNone, kNone, "LP")); // fill
  sy.push_back(uncert("time", "Time stability", jec::kTime,
  		      "default", "default", -1, // defaults
  		      kMagenta+2, kNone, 1, // line
  		      kMagenta+2, kOpenTriangleDown, // marker
  		      kNone, kNone, "LP")); // fill

  vector<uncert> sym;
  sym.push_back(uncert("tot", "Total uncertainty", jec::kMC,
		       "default", "default", -1, // defaults
		       kBlack, kSolid, 1, // line
		       kBlack, kNone, // marker
		       kDarkGray, 1001, "LF")); // fill
  sym.push_back(uncert("absolute", "Absolute scale", jec::kAbsolute,
		       "default", "default", -1, // defaults
		       kYellow+3, kSolid, 1, // line
		       kBlack, kNone, // marker
		       kYellow, 1001, "LF")); // fill
  sym.push_back(uncert("relative", "Relative scale",jec::kRelative,
		       "default", "default", -1, // defaults
		       kBlack, kSolid, 1, // line
		       kBlack, kFullTriangleDown, // marker
		       kNone, kNone, "LP")); // fill
  sym.push_back(uncert("highpt", "Extrapolation", jec::kAbsolutePt,//kPtExtra,
		       "default", "default", -1, // defaults
		       kRed, kSolid, 1, // line
		       kRed, kOpenCircle, // marker
		       kNone, kNone, "LP")); // fill
  sym.push_back(uncert("pileup", Form("Pile-up, #LT#mu#GT=%1.0f",d_mu),
		       jec::kPileUpDataMC,
		       "default", "default", -1, // defaults
		       kBlue, kNone, 1, // line
		       kBlue, kFullSquare, // marker
		       kNone, kNone, "LP")); // fill
  sym.push_back(uncert("flavjec", "Jet flavor (QCD)", jec::kFlavorQCD,
		       "default", "default", -1, // defaults
		       kGreen+2, kSolid, 1, // line
		       kGreen+2, kOpenSquare, // marker
		       kNone, kNone, "LP")); // fill
  sym.push_back(uncert("time", "Time stability", jec::kTime,
		       "default", "default", -1, // defaults
		       kMagenta+2, kNone, 1, // line
		       kMagenta+2, kOpenTriangleDown, // marker
		       kNone, kNone, "LP")); // fill


  vector<uncert> sypu;
  sypu.push_back(uncert("pileup", "SubTotalPileUp",
			jec::kPileUp,
			"default", "default", -1, // defaults
			kBlue, kNone, 1, // line
			kBlue, kFullSquare, // marker
			kBlue-9, 1001, "LFP")); // fill
  sypu.push_back(uncert("pupt", "PileUpPtEta",
			jec::kPileUpPt & ~jec::kPileUpPtRef,
			"default", "default", -1, // defaults
			kYellow+3, kNone, 1, // line
			kBlack, kNone, // marker
			kYellow, 1001, "LF")); // fill
  sypu.push_back(uncert("puref", "PileUpPtRef",
			jec::kPileUpPtRef,
			"default", "default", -1, // defaults
			kRed, kNone, 1, // line
			kRed, kOpenCircle, // marker
			kNone, kNone, "LP")); // fill
  sypu.push_back(uncert("pudatamc", "PileUpDataMC",
			jec::kPileUpDataMC,
			"default", "default", -1, // defaults
			kGreen+2, kNone, 1, // line
			kGreen+2, kFullTriangleDown, // marker
			kNone, kNone, "LP")); // fill
  sypu.push_back(uncert("puzero", "PileUpMuZero (opt)",
			jec::kPileUpMuZero,
			"default", "default", -1, // defaults
			kGray+2, kNone, 1, // line
			kGray+2, kOpenDiamond, // marker
			kNone, kNone, "LP")); // fill
  // Remove from paper plots after pas-v6
  sypu.push_back(uncert("puenvelope", "PU envelope (opt)",
  		jec::kPileUpEnvelope,
  		"default", "default", -1, // defaults
  		kGray+1, kSolid, 1, // line
  		kGray+1, kNone, // marker
  		kNone, kNone, "LP")); // fill

  vector<uncert> syrel;
  syrel.push_back(uncert("relative", "SubTotalRelative",
			 jec::kRelative,
			 "default", "default", -1, // defaults
			 kBlack, kNone, 1, // line
			 kBlack, kFullTriangleDown, // marker
			 kBlue-9, 1001, "LFP")); // fill
  syrel.push_back(uncert("relpt", "RelativePt",
			 jec::kRelativePt,
			 "default", "default", -1, // defaults
			 kYellow+3, kNone, 1, // line
			 kBlack, kNone, // marker
			 kYellow, 1001, "LF")); // fill
  syrel.push_back(uncert("relbal", "RelativeBal",
			 jec::kRelativeBal,
			 "default", "default", -1, // defaults
			 kRed+1, kNone, 1, // line
			 kRed+1, kFullTriangleDown, // marker
			 kBlack, kNone, "LP")); // fill
  syrel.push_back(uncert("relsam", "RelativeSample",
			 jec::kRelativeSample,
			 "default", "default", -1, // defaults
			 kBlue+1, kNone, 1, // line
			 kBlue+1, kFullTriangleUp, // marker
			 kBlack, kNone, "LP")); // fill
  syrel.push_back(uncert("reljer", "RelativeJER",
			 jec::kRelativeJER,
			 "default", "default", -1, // defaults
			 kRed, kNone, 1, // line
			 kRed, kOpenCircle, // marker
			 kBlack, kNone, "LP")); // fill
  syrel.push_back(uncert("relfsr", "RelativeFSR",
			 jec::kRelativeFSR,
			 "default", "default", -1, // defaults
			 kBlue, kNone, 1, // line
			 kBlue, kFullSquare, // marker
			 kBlack, kNone, "LP")); // fill
  syrel.push_back(uncert("relstat", "RelativeStat",
			 jec::kRelativeStat,
			 "default", "default", -1, // defaults
			 kGreen+2, kNone, 1, // line
			 kGreen+2, kOpenSquare, // marker
			 kBlack, kNone, "LP")); // fill

  vector<uncert> sypt;
  sypt.push_back(uncert("abspt", "SubTotalAbsolute",
			jec::kAbsolute,
			"default", "default", -1, // defaults
			kRed, kNone, 1, // line
			kBlack, kNone, // marker
			kRed-9, 1001, "LFP")); // fill
  sypt.push_back(uncert("absscale", "AbsoluteScale", jec::kAbsoluteScale,
			"default", "default", -1, // defaults
			kYellow+3, kSolid, 1, // line
			kBlack, kNone, // marker
			kYellow, 1001, "LF")); // fill
  sypt.push_back(uncert("absstat", "AbsoluteStat", jec::kAbsoluteStat,
			"default", "default", -1, // defaults
			kGray+2, kSolid, 1, // line
			kGray+2, kNone, // marker
			kGray, 1001, "LF")); // fill
  //
  sypt.push_back(uncert("absfrag", "Fragmentation",
			jec::kAbsoluteFrag,
			"default", "default", -1, // defaults
			kRed, kNone, 1, // line
			kRed, kOpenCircle, // marker
			kNone, kNone, "LP")); // fill
  sypt.push_back(uncert("absspre", "SinglePionECAL",
			jec::kAbsoluteSPRE,
			"default", "default", -1, // defaults
			kBlack, kNone, 1, // line
			kBlack, kFullTriangleDown, // marker
			kNone, kNone, "LP")); // fill
  sypt.push_back(uncert("abssprh", "SinglePionHCAL",
			jec::kAbsoluteSPRH,
			"default", "default", -1, // defaults
			kBlue, kNone, 1, // line
			kBlue, kFullSquare, // marker
			kNone, kNone, "LP")); // fill
  sypt.push_back(uncert("absstat", "MPFBias",
			jec::kAbsoluteMPFBias,
  			"default", "default", -1, // defaults
  			kGreen+2, kNone, 1, // line
  			kGreen+2, kOpenSquare, // marker
  			kNone, kNone, "LP")); // fill

  vector<uncert> syf;
  syf.push_back(uncert("flavor_gluon", "Gluons",
                       jec::kFlavorPureGluon,
                      "default", "default", -1, // defaults
                      kYellow+3, kSolid, 1, // line
                      kBlack, kNone, // marker
                      kYellow, 1001, "LF")); // fill
  syf.push_back(uncert("flavor_qcd", "QCD mixture",
                       jec::kFlavorQCD,
                      "default", "default", -1, // defaults
                      kGreen+2, kSolid, 1, // line
                      kGreen+2, kOpenSquare, // marker
                      kNone, kNone, "LP")); // fill
  syf.push_back(uncert("flavor_zjet", "Z+jet mixture",
                       jec::kFlavorZJet,
                       "default", "default", -1, // defaults
                       kRed, kNone, 1, // line
                       kRed, kOpenCircle, // marker
                       kNone, kNone, "LP")); // fill
  syf.push_back(uncert("flavor_zjet", "#gamma+jet mixture",
                       jec::kFlavorPhotonJet,
                       "default", "default", -1, // defaults
                       kMagenta+1, kNone, 1, // line
                       kMagenta+1, kOpenDiamond, // marker
                       kNone, kNone, "LP")); // fill
  syf.push_back(uncert("flavor", "Light quarks",
                       jec::kFlavorPureQuark,
                       "default", "default", -1, // defaults
                       kBlue, kNone, 1, // line
                       kBlue, kFullTriangleUp, // marker
                       kNone, kNone, "LP")); // fill
  syf.push_back(uncert("flavor", "Bottom quarks",
                       jec::kFlavorPureBottom,
                       "default", "default", -1, // defaults
                       kBlack, kNone, 1, // line
                       kBlack, kFullTriangleDown, // marker
                       kNone, kNone, "LP")); // fill

  vector<uncert> syt;
  syt.push_back(uncert("time_eta", "TimePtEta",
		       jec::kTimePtEta,
		       "default", "default", -1, // defaults
		       kYellow+3, kSolid, 1, // line
		       kBlack, kNone, // marker
		       kYellow, 1001, "LF")); // fill
  syt.push_back(uncert("time_runa", "TimeRunB",
                       jec::kTimeRunB,
                       "default", "default", -1, // defaults
                       kRed, kNone, 1, // line
                       kRed, kOpenCircle, // marker
                       kNone, kNone, "LP")); // fill
  syt.push_back(uncert("flavor_runb", "TimeRunC",
                       jec::kTimeRunC,
                       "default", "default", -1, // defaults
                       kMagenta+1, kNone, 1, // line
                       kMagenta+1, kOpenDiamond, // marker
                       kNone, kNone, "LP")); // fill
  syt.push_back(uncert("time_runc", "TimeRunD",
                       jec::kTimeRunD,
                       "default", "default", -1, // defaults
                       kBlue, kNone, 1, // line
                       kBlue, kFullTriangleUp, // marker
                       kNone, kNone, "LP")); // fill
  syt.push_back(uncert("time_runc", "TimeRunE",
                       jec::kTimeRunE,
                       "default", "default", -1, // defaults
                       kBlue+2, kNone, 1, // line
                       kBlue+2, kOpenTriangleUp, // marker
                       kNone, kNone, "LP")); // fill
  syt.push_back(uncert("time_rund", "TimeRunF",
                       jec::kTimeRunF,
                       "default", "default", -1, // defaults
                       kBlack, kNone, 1, // line
                       kBlack, kFullTriangleDown, // marker
                       kNone, kNone, "LP")); // fill

  vector<uncert> sys;
  sys.push_back(uncert("singlepion", "SinglePionHCAL",
                       jec::kAbsoluteSPRH,
		       "default", "default", -1, // defaults
		       kYellow+3, kSolid, 1, // line
		       kBlack, kNone, // marker
		       kYellow, 1001, "LF")); // fill
  sys.push_back(uncert("pileup", "PileUpEnvelope",
                       jec::kPileUpEnvelope,
		       "default", "default", -1, // defaults
		       kOrange+2, kDotted, 2, // line
		       kNone, kNone, // marker
		       kNone, kNone, "L")); // fill
  sys.push_back(uncert("relativept", "RelativePt",
                       jec::kRelativePt,
		       "default", "default", -1, // defaults
		       kBlack, kSolid, 1, // line
		       kBlack, kFullCircle, // marker
		       kNone, kNone, "LP")); // fill

  vector<uncert> syCorrGroups;
  syCorrGroups.push_back(uncert("tot", "Total uncertainty", jec::kData,
		      "default", "default", -1, // defaults
		      kBlack, kSolid, 1, // line
		      kBlack, kNone, // marker
		      kDarkGray, 1001, "LF")); // fill
  syCorrGroups.push_back(uncert("uncorr", "Uncorrelated group", jec::kCorrelationGroupUncorrelated,
		      "default", "default", -1, // defaults
		      kYellow+3, kSolid, 1, // line
		      kBlack, kNone, // marker
		      kYellow, 1001, "LF")); // fill
  syCorrGroups.push_back(uncert("insitu", "MPF/in-situ group",jec::kCorrelationGroupMPFInSitu,
		      "default", "default", -1, // defaults
		      kBlack, kSolid, 1, // line
		      kBlack, kFullTriangleDown, // marker
		      kNone, kNone, "LP")); // fill
  syCorrGroups.push_back(uncert("flavorgroup", "Flavor group", jec::kCorrelationGroupFlavor,
		      "default", "default", -1, // defaults
		      kRed, kSolid, 1, // line
		      kRed, kOpenCircle, // marker
		      kNone, kNone, "LP")); // fill
  syCorrGroups.push_back(uncert("intercalibration", "Intercalibration group", jec::kCorrelationGroupIntercalibration,
		      "default", "default", -1, // defaults
		      kBlue, kNone, 1, // line
		      kBlue, kFullSquare, // marker
		      kNone, kNone, "LP")); // fill
  syCorrGroups.push_back(uncert("bjes", "b-JES group", jec::kCorrelationGroupbJES,
		      "default", "default", -1, // defaults
		      kGreen+2, kSolid, 1, // line
		      kGreen+2, kOpenSquare, // marker
		      kNone, kNone, "LP")); // fill



  jec::JetAlgo jetAlg = d_algo;

  double r(0);
  if (jetAlg==jec::AK8PF||jetAlg==jec::AK8PFchs||jetAlg==jec::AK8CALO) r = 0.8;
  if (jetAlg==jec::AK7PF||jetAlg==jec::AK7PFchs||jetAlg==jec::AK7CALO) r = 0.7;
  if (jetAlg==jec::AK5PF||jetAlg==jec::AK5PFchs||jetAlg==jec::AK5CALO) r = 0.5;
  if (jetAlg==jec::AK4PF||jetAlg==jec::AK4PFchs||jetAlg==jec::AK4CALO) r = 0.4;
  string sa = "";
  if (jetAlg==jec::AK4CALO||jetAlg==jec::AK5CALO||
      jetAlg==jec::AK7CALO||jetAlg==jec::AK8CALO) sa = "Calo";
  if (jetAlg==jec::AK4PFchs||jetAlg==jec::AK5PFchs||
      jetAlg==jec::AK7PFchs||jetAlg==jec::AK8PFchs) sa = "PF+CHS";
  if (jetAlg==jec::AK4PF||jetAlg==jec::AK5PF||
      jetAlg==jec::AK7PF||jetAlg==jec::AK8PF) sa = "PF";
  string ss = Form("R=%1.1f %s", r, sa.c_str());
  const char *s = ss.c_str();

  const char *cu = (_absUncert ? "JECUncert" : "JECSource");

  map<jec::JetAlgo, const char*> names;
  _algnames = &names;
  names[jec::AK4PF] = "AK4PF";
  names[jec::AK4PFchs] = "AK4PFchs";
  names[jec::AK5PF] = "AK5PF";
  names[jec::AK5PFchs] = "AK5PFchs";
  names[jec::AK7PF] = "AK7PF";
  names[jec::AK7PFchs] = "AK7PFchs";
  names[jec::AK8PF] = "AK8PF";
  names[jec::AK8PFchs] = "AK8PFchs";
  names[jec::AK4CALO] = "AK4CALO";
  names[jec::AK5CALO] = "AK5CALO";
  names[jec::AK7CALO] = "AK7CALO";
  names[jec::AK8CALO] = "AK8CALO";

  string ssd = Form("%s_DATA_Summary_%s", cu, names[jetAlg]);
  const char *sd = ssd.c_str();

  TCanvas *c(0);
   if (_fourbytwo) {
     c  = new TCanvas("c4x2","c4x2",2400,1200);
     c->Divide(4,2);
   }
   if (_twobythree) {
     c  = new TCanvas("c2x3","c2x3",1200,1800);
     c->Divide(2,3);
   }
   //assert(c);
  _canvas = c;
  _icanvas = 1;

  bool minimaltmp = _minimal;
  if (algo=="AK4PFchs") _minimal = false;

  double ym =  (jetAlg==jec::AK4PF || jetAlg==jec::AK4PFchs ?  18. : 24.);
  double ym0 =  (jetAlg==jec::AK4PF || jetAlg==jec::AK4PFchs ?  7. : 10.);

  // Data uncertainty
  // vs pT
  if (_paper && _absUncert) {
  plotUncertainty(sy, 0, sy.size(), jetAlg, Form("%s_Eta00",sd),
		  //"JEC uncertainty", s, "|#eta_{jet}|=0", ym,10,"fixEta",0.);
		  "JEC uncertainty",s, "|#eta_{jet}| = 0", ym0,10,"fixEta",0.);
  plotUncertainty(sy, 0, sy.size(), jetAlg, Form("%s_Eta27",sd),
		  "JEC uncertainty",s, "|#eta_{jet}| = 2.7",ym,10,"fixEta",2.7);
  plotUncertainty(sy, 0, sy.size(), jetAlg, Form("%s_Eta35",sd),
		  "JEC uncertainty",s, "|#eta_{jet}| = 3.5",ym,10,"fixEta",3.5);
  //
  plotUncertainty(sy, 0, sy.size(), jetAlg, Form("%s_Pt30",sd),
		  "JEC uncertainty", s, "p_{T} = 30 GeV", ym,10,"fixPt",30.);
  plotUncertainty(sy, 0, sy.size(), jetAlg, Form("%s_Pt100",sd),
		  //"JEC uncertainty", s, "p_{T}=100 GeV", ym,10,"fixPt",100.);
		  "JEC uncertainty", s, "p_{T} = 100 GeV", 9,10,"fixPt",100.);
  plotUncertainty(sy, 0, sy.size(), jetAlg, Form("%s_Pt1000",sd),
		  "JEC uncertainty", s, "p_{T} = 1000 GeV",ym0,10,"fixPt",1000);
  }
  if (_fourbytwo) {

  plotUncertainty(sy, 0, sy.size(), jetAlg, Form("%s_Eta00",sd),
		  "JEC uncertainty", s, "|#eta_{jet}|=0", ym,10,"fixEta",0.);
  plotUncertainty(sy, 0, sy.size(), jetAlg, Form("%s_Eta20",sd),
		  "JEC uncertainty", s, "|#eta_{jet}|=2.0", ym,10,"fixEta",2.0);
  plotUncertainty(sy, 0, sy.size(), jetAlg, Form("%s_Eta27",sd),
		  "JEC uncertainty", s, "|#eta_{jet}|=2.7", ym,10,"fixEta",2.7);
  plotUncertainty(sy, 0, sy.size(), jetAlg, Form("%s_Eta42",sd),
		  "JEC uncertainty", s, "|#eta_{jet}|=4.2", ym,10,"fixEta",4.2);
  // vs eta
  plotUncertainty(sy, 0, sy.size(), jetAlg, Form("%s_Pt30",sd),
		  "JEC uncertainty", s, "p_{T}=30 GeV", ym,10,"fixPt",30.);
  plotUncertainty(sy, 0, sy.size(), jetAlg, Form("%s_Pt100",sd),
		  "JEC uncertainty", s, "p_{T}=100 GeV", ym,10,"fixPt",100.);
  plotUncertainty(sy, 0, sy.size(), jetAlg, Form("%s_E500",sd),
		  "JEC uncertainty", s, "E=500 GeV", ym,10,"fixE",500.);
  plotUncertainty(sy, 0, sy.size(), jetAlg, Form("%s_E2000",sd),
		  "JEC uncertainty", s, "E=2000 GeV", ym,10,"fixE",2000.);
  plotUncertainty(sy, 0, sy.size(), jetAlg, Form("%s_E1000",sd),
		  "JEC uncertainty", s, "E=1000 GeV", ym,10,"fixE",1000.);//last
  }
  if (_twobythree) {

  plotUncertainty(sy, 0, sy.size(), jetAlg, Form("%s_Eta00",sd),
		  "JEC uncertainty", s, "|#eta_{jet}|=0", ym,10,"fixEta",0.);
  plotUncertainty(sy, 0, sy.size(), jetAlg, Form("%s_Pt30",sd),
		  "JEC uncertainty", s, "p_{T}=30 GeV", ym,10,"fixPt",30.);
  //
  plotUncertainty(sy, 0, sy.size(), jetAlg, Form("%s_Eta27",sd),
		  "JEC uncertainty", s, "|#eta_{jet}|=2.7", ym,10,"fixEta",2.7);
  plotUncertainty(sy, 0, sy.size(), jetAlg, Form("%s_Pt100",sd),
		  "JEC uncertainty", s, "p_{T}=100 GeV", ym,10,"fixPt",100.);
  //
  plotUncertainty(sy, 0, sy.size(), jetAlg, Form("%s_Eta35",sd),
		  "JEC uncertainty", s, "|#eta_{jet}|=3.5", ym,10,"fixEta",3.5);
  plotUncertainty(sy, 0, sy.size(), jetAlg, Form("%s_Pt1000",sd),
		  "JEC uncertainty", s, "p_{T}=1000 GeV", ym,10,"fixPt",1000.);
  }

  if (_canvas) _canvas->SaveAs(Form("pdf/%s.pdf",sd));
  _icanvas = 1;

  _minimal = minimaltmp;

  string ssm = Form("%s_MC_Summary_%s",cu,names[jetAlg]);
  const char *sm = ssm.c_str();

  // Data/MC (MC) uncertainty
  // vs pT
  if (_fourbytwo) {

  plotUncertainty(sym, 0, sym.size(), jetAlg, Form("%s_Eta00",sm),
		  "JEC uncertainty", s,
		  "|#eta_{jet}|=0", ym,10,"fixEta",0.);
  plotUncertainty(sym, 0, sym.size(), jetAlg, Form("%s_Eta20",sm),
		  "JEC uncertainty", s,
		  "|#eta_{jet}|=2.0", ym,10,"fixEta",2.0);
  plotUncertainty(sym, 0, sym.size(), jetAlg, Form("%s_Eta27",sm),
		  "JEC uncertainty", s,
		  "|#eta_{jet}|=2.7", ym,10,"fixEta",2.7);
  plotUncertainty(sym, 0, sym.size(), jetAlg, Form("%s_Eta42",sm),
		  "JEC uncertainty", s,
		  "|#eta_{jet}|=4.2", ym,10,"fixEta",4.2);
  // vs eta
  plotUncertainty(sym, 0, sym.size(), jetAlg, Form("%s_Pt30",sm),
		  "JEC uncertainty", s,
		  "p_{T}=30 GeV", ym,10,"fixPt",30.);
  plotUncertainty(sym, 0, sym.size(), jetAlg, Form("%s_Pt100",sm),
		  "JEC uncertainty", s,
		  "p_{T}=100 GeV", ym,10,"fixPt",100.);
  plotUncertainty(sym, 0, sym.size(), jetAlg, Form("%s_E500",sm),
		  "JEC uncertainty", s,
		  "E=500 GeV", ym,10,"fixE",500.);
  plotUncertainty(sym, 0, sym.size(), jetAlg, Form("%s_E2000",sm),
		  "JEC uncertainty", s,
		  "E=2000 GeV", ym,10,"fixE",2000.);
  plotUncertainty(sym, 0, sym.size(), jetAlg, Form("%s_E1000",sm),
		  "JEC uncertainty", s,
		  "E=1000 GeV", ym,10,"fixE",1000.); // moved last
  }
  if (_twobythree) {

  plotUncertainty(sym, 0, sym.size(), jetAlg, Form("%s_Eta00",sm),
		  "JEC uncertainty", s,
		  "|#eta_{jet}|=0", ym,10,"fixEta",0.);
  plotUncertainty(sym, 0, sym.size(), jetAlg, Form("%s_Pt30",sm),
		  "JEC uncertainty", s,
		  "p_{T}=30 GeV", ym,10,"fixPt",30.);
  //
  plotUncertainty(sym, 0, sym.size(), jetAlg, Form("%s_Eta27",sm),
		  "JEC uncertainty", s,
		  "|#eta_{jet}|=2.7", ym,10,"fixEta",2.7);
  plotUncertainty(sym, 0, sym.size(), jetAlg, Form("%s_Pt100",sm),
		  "JEC uncertainty", s,
		  "p_{T}=100 GeV", ym,10,"fixPt",100.);
  //
  plotUncertainty(sym, 0, sym.size(), jetAlg, Form("%s_Eta35",sm),
		  "JEC uncertainty", s,
		  "|#eta_{jet}|=3.5", ym,10,"fixEta",3.5);
  plotUncertainty(sym, 0, sym.size(), jetAlg, Form("%s_Pt500",sm),
		  "JEC uncertainty", s,
		  "p_{T}=500 GeV", ym,10,"fixPt",500.);
  }

  if (_canvas) _canvas->SaveAs(Form("pdf/%s.pdf",sm));
  _icanvas = 1;

  string sspu = Form("%s_PileUp_%s",cu,names[jetAlg]);
  const char *spu = sspu.c_str();

  double ymaxpu = 18;
  double ymaxpu0 = 10;

  // PU uncertainty
  // vs pT
  if (_paper && _absUncert) {
  plotUncertainty(sypu, 0, sypu.size(), jetAlg, Form("%s_Eta00",spu),
		  "JEC uncertainty", s,
		  "|#eta_{jet}| = 0", ymaxpu0,10,"fixEta",0.);
  plotUncertainty(sypu, 0, sypu.size(), jetAlg, Form("%s_Eta27",spu),
		  "JEC uncertainty", s,
		  "|#eta_{jet}| = 2.7", ymaxpu,10,"fixEta",2.7);
  plotUncertainty(sypu, 0, sypu.size(), jetAlg, Form("%s_Eta35",spu),
		  "JEC uncertainty", s,
		  "|#eta_{jet}| = 3.5", ymaxpu,10,"fixEta",3.5);
  plotUncertainty(sypu, 0, sypu.size(), jetAlg, Form("%s_Pt30",spu),
		  "JEC uncertainty", s,
		  "p_{T} = 30 GeV", ymaxpu,10,"fixPt",30.);
  }
  if (_fourbytwo) {

  plotUncertainty(sypu, 0, sypu.size(), jetAlg, Form("%s_Eta00",spu),
		  "JEC uncertainty", s,
		  "|#eta_{jet}|=0", ymaxpu,10,"fixEta",0.);
  plotUncertainty(sypu, 0, sypu.size(), jetAlg, Form("%s_Eta20",spu),
		  "JEC uncertainty", s,
		  "|#eta_{jet}|=2.0", ymaxpu,10,"fixEta",2.0);
  plotUncertainty(sypu, 0, sypu.size(), jetAlg, Form("%s_Eta27",spu),
		  "JEC uncertainty", s,
		  "|#eta_{jet}|=2.7", ymaxpu,10,"fixEta",2.7);
  plotUncertainty(sypu, 0, sypu.size(), jetAlg, Form("%s_Eta42",spu),
		  "JEC uncertainty", s,
		  "|#eta_{jet}|=4.2", ymaxpu,10,"fixEta",4.2);
  // vs eta
  plotUncertainty(sypu, 0, sypu.size(), jetAlg, Form("%s_Pt30",spu),
		  "JEC uncertainty", s,
		  "p_{T}=30 GeV", ymaxpu,10,"fixPt",30.);
  plotUncertainty(sypu, 0, sypu.size(), jetAlg, Form("%s_Pt100",spu),
		  "JEC uncertainty", s,
		  "p_{T}=100 GeV", ymaxpu,10,"fixPt",100.);
  plotUncertainty(sypu, 0, sypu.size(), jetAlg, Form("%s_E500",spu),
		  "JEC uncertainty", s,
		  "E=500 GeV", ymaxpu,10,"fixE",500.);
  plotUncertainty(sypu, 0, sypu.size(), jetAlg, Form("%s_E2000",spu),
		  "JEC uncertainty", s,
		  "E=2000 GeV", ymaxpu,10,"fixE",2000.);
  plotUncertainty(sypu, 0, sypu.size(), jetAlg, Form("%s_E1000",spu),
		  "JEC uncertainty", s,
		  "E=1000 GeV", ymaxpu,10,"fixE",1000.); // moved last
  }
  if (_twobythree) {
    
  plotUncertainty(sypu, 0, sypu.size(), jetAlg, Form("%s_Eta00",spu),
		  "JEC uncertainty", s,
		  "|#eta_{jet}|=0", ymaxpu,10,"fixEta",0.);
  plotUncertainty(sypu, 0, sypu.size(), jetAlg, Form("%s_Pt30",spu),
		  "JEC uncertainty", s,
		  "p_{T}=30 GeV", ymaxpu,10,"fixPt",30.);
  //
  plotUncertainty(sypu, 0, sypu.size(), jetAlg, Form("%s_Eta27",spu),
		  "JEC uncertainty", s,
		  "|#eta_{jet}|=2.7", ymaxpu,10,"fixEta",2.7);
  plotUncertainty(sypu, 0, sypu.size(), jetAlg, Form("%s_Pt100",spu),
		  "JEC uncertainty", s,
		  "p_{T}=100 GeV", ymaxpu,10,"fixPt",100.);
  //
  plotUncertainty(sypu, 0, sypu.size(), jetAlg, Form("%s_Eta35",spu),
		  "JEC uncertainty", s,
		  "|#eta_{jet}|=3.5", ymaxpu,10,"fixEta",3.5);
  plotUncertainty(sypu, 0, sypu.size(), jetAlg, Form("%s_Pt500",spu),
		  "JEC uncertainty", s,
		  "p_{T}=500 GeV", ymaxpu,10,"fixPt",500.);
  }

  if (_canvas) _canvas->SaveAs(Form("pdf/%s.pdf",spu));
  _icanvas = 1;

  string ssrel = Form("%s_Relative_%s",cu,names[jetAlg]);
  const char *srel = ssrel.c_str();

  double ymaxrel = 12;//18.;

  // Relative scale uncertainty
  // vs pT
  if (_paper && _absUncert) {

  plotUncertainty(syrel, 0, syrel.size(), jetAlg, Form("%s_Eta00",srel),
		  "JEC uncertainty", s,
		  "|#eta_{jet}| = 0.0", ymaxrel/3,10,"fixEta",0.0);
  plotUncertainty(syrel, 0, syrel.size(), jetAlg, Form("%s_Eta17",srel),
		  "JEC uncertainty", s,
		  "|#eta_{jet}| = 1.7", ymaxrel/3,10,"fixEta",1.7);
  plotUncertainty(syrel, 0, syrel.size(), jetAlg, Form("%s_Eta27",srel),
		  "JEC uncertainty", s,
		  "|#eta_{jet}| = 2.7", ymaxrel,10,"fixEta",2.7);

  plotUncertainty(syrel, 0, syrel.size(), jetAlg, Form("%s_Pt30",srel),
		  "JEC uncertainty", s,
		  "p_{T} = 30 GeV", ymaxrel,10,"fixPt",30.);
  plotUncertainty(syrel, 0, syrel.size(), jetAlg, Form("%s_Pt100",srel),
		  "JEC uncertainty", s,
		  "p_{T} = 100 GeV", ymaxrel,10,"fixPt",100.);
  plotUncertainty(syrel, 0, syrel.size(), jetAlg, Form("%s_Pt500",srel),
		  "JEC uncertainty", s,
		  "p_{T} = 500 GeV", ymaxrel,10,"fixPt",500.);
  }
  if (_fourbytwo) {

  plotUncertainty(syrel, 0, syrel.size(), jetAlg, Form("%s_Eta00",srel),
		  "JEC uncertainty", s,
		  "|#eta_{jet}|=0", ymaxrel,10,"fixEta",0.);
  plotUncertainty(syrel, 0, syrel.size(), jetAlg, Form("%s_Eta20",srel),
		  "JEC uncertainty", s,
		  "|#eta_{jet}|=2.0", ymaxrel,10,"fixEta",2.0);
  plotUncertainty(syrel, 0, syrel.size(), jetAlg, Form("%s_Eta27",srel),
		  "JEC uncertainty", s,
		  "|#eta_{jet}|=2.7", ymaxrel,10,"fixEta",2.7);
  plotUncertainty(syrel, 0, syrel.size(), jetAlg, Form("%s_Eta42",srel),
		  "JEC uncertainty", s,
		  "|#eta_{jet}|=4.2", ymaxrel,10,"fixEta",4.2);
  // vs eta
  plotUncertainty(syrel, 0, syrel.size(), jetAlg, Form("%s_Pt30",srel),
		  "JEC uncertainty", s,
		  "p_{T}=30 GeV", ymaxrel,10,"fixPt",30.);
  plotUncertainty(syrel, 0, syrel.size(), jetAlg, Form("%s_Pt100",srel),
		  "JEC uncertainty", s,
		  "p_{T}=100 GeV", ymaxrel,10,"fixPt",100.);
  plotUncertainty(syrel, 0, syrel.size(), jetAlg, Form("%s_E500",srel),
		  "JEC uncertainty", s,
		  "E=500 GeV", ymaxrel,10,"fixE",500.);
  plotUncertainty(syrel, 0, syrel.size(), jetAlg, Form("%s_E2000",srel),
		  "JEC uncertainty", s,
		  "E=2000 GeV", ymaxrel,10,"fixE",2000.);
  plotUncertainty(syrel, 0, syrel.size(), jetAlg, Form("%s_E1000",srel),
		  "JEC uncertainty", s,
		  "E=1000 GeV", ymaxrel,10,"fixE",1000.); // moved last
  }
  if (_twobythree) {

  plotUncertainty(syrel, 0, syrel.size(), jetAlg, Form("%s_Eta00",srel),
		  "JEC uncertainty", s,
		  "|#eta_{jet}|=0", ymaxrel,10,"fixEta",0.);
  plotUncertainty(syrel, 0, syrel.size(), jetAlg, Form("%s_Pt30",srel),
		  "JEC uncertainty", s,
		  "p_{T}=30 GeV", ymaxrel,10,"fixPt",30.);
  //
  plotUncertainty(syrel, 0, syrel.size(), jetAlg, Form("%s_Eta27",srel),
		  "JEC uncertainty", s,
		  "|#eta_{jet}|=2.7", ymaxrel,10,"fixEta",2.7);
  plotUncertainty(syrel, 0, syrel.size(), jetAlg, Form("%s_Pt100",srel),
		  "JEC uncertainty", s,
		  "p_{T}=100 GeV", ymaxrel,10,"fixPt",100.);
  //
  plotUncertainty(syrel, 0, syrel.size(), jetAlg, Form("%s_Eta35",srel),
		  "JEC uncertainty", s,
		  "|#eta_{jet}|=3.5", ymaxrel,10,"fixEta",3.5);
  plotUncertainty(syrel, 0, syrel.size(), jetAlg, Form("%s_Pt500",srel),
		  "JEC uncertainty", s,
		  "p_{T}=500 GeV", ymaxrel,10,"fixPt",500.);
  }

  if (_canvas) _canvas->SaveAs(Form("pdf/%s.pdf",srel));
  _icanvas = 1;


  string ssCorrGroups = Form("%s_CorrelationGroups_%s",cu,names[jetAlg]);
  const char *sCorrGroups = ssCorrGroups.c_str();

  // Correlation groups
  // vs pT
  if (_fourbytwo) {

  plotUncertainty(syCorrGroups, 0, syCorrGroups.size(), jetAlg, Form("%s_Eta00",sCorrGroups),
		  "JEC uncertainty", "R=0.4 PF", "|#eta_{jet}|=0",
		  10,10,"fixEta",0.);
  plotUncertainty(syCorrGroups, 0, syCorrGroups.size(), jetAlg, Form("%s_Eta20",sCorrGroups),
		  "JEC uncertainty", "R=0.4 PF", "|#eta_{jet}|=2.0",
		  10,10,"fixEta",2.0);
  plotUncertainty(syCorrGroups, 0, syCorrGroups.size(), jetAlg, Form("%s_Eta27",sCorrGroups),
		  "JEC uncertainty", "R=0.4 PF", "|#eta_{jet}|=2.7",
		  10,10,"fixEta",2.7);
  plotUncertainty(syCorrGroups, 0, syCorrGroups.size(), jetAlg, Form("%s_Eta42",sCorrGroups),
		  "JEC uncertainty", "R=0.4 PF", "|#eta_{jet}|=4.2",
		  10,10,"fixEta",4.2);
  // vs eta
  plotUncertainty(syCorrGroups, 0, syCorrGroups.size(), jetAlg, Form("%s_Pt30",sCorrGroups),
		  "JEC uncertainty", "R=0.4 PF", "p_{T}=30 GeV",
		  10,10,"fixPt",30.);
  plotUncertainty(syCorrGroups, 0, syCorrGroups.size(), jetAlg, Form("%s_Pt100",sCorrGroups),
		  "JEC uncertainty", "R=0.4 PF", "p_{T}=100 GeV",
		  10,10,"fixPt",100.);
  plotUncertainty(syCorrGroups, 0, syCorrGroups.size(), jetAlg, Form("%s_E500",sCorrGroups),
		  "JEC uncertainty", "R=0.4 PF", "E=500 GeV",
		  10,10,"fixPt",500.);
  plotUncertainty(syCorrGroups, 0, syCorrGroups.size(), jetAlg, Form("%s_E2000",sCorrGroups),
		  "JEC uncertainty", "R=0.4 PF", "E=2000 GeV",
		  10,10,"fixE",500.);
  plotUncertainty(syCorrGroups, 0, syCorrGroups.size(), jetAlg, Form("%s_E1000",sCorrGroups),
		  "JEC uncertainty", "R=0.4 PF", "E=1000 GeV",
		  10,10,"fixE",2000.);
  }
  if (_twobythree) {

  plotUncertainty(syCorrGroups, 0, syCorrGroups.size(), jetAlg, Form("%s_Eta00",sCorrGroups),
		  "JEC uncertainty", "R=0.4 PF", "|#eta_{jet}|=0",
		  10,10,"fixEta",0.);
  plotUncertainty(syCorrGroups, 0, syCorrGroups.size(), jetAlg, Form("%s_Pt30",sCorrGroups),
		  "JEC uncertainty", "R=0.4 PF", "p_{T}=30 GeV",
		  10,10,"fixPt",30.);
  //
  plotUncertainty(syCorrGroups, 0, syCorrGroups.size(), jetAlg, Form("%s_Eta27",sCorrGroups),
		  "JEC uncertainty", "R=0.4 PF", "|#eta_{jet}|=2.7",
		  10,10,"fixEta",2.7);
  plotUncertainty(syCorrGroups, 0, syCorrGroups.size(), jetAlg, Form("%s_Pt100",sCorrGroups),
		  "JEC uncertainty", "R=0.4 PF", "p_{T}=100 GeV",
		  10,10,"fixPt",100.);
  //
  plotUncertainty(syCorrGroups, 0, syCorrGroups.size(), jetAlg, Form("%s_Eta35",sCorrGroups),
		  "JEC uncertainty", "R=0.4 PF", "|#eta_{jet}|=3.5",
		  10,10,"fixEta",3.5);
  plotUncertainty(syCorrGroups, 0, syCorrGroups.size(), jetAlg, Form("%s_Pt500",sCorrGroups),
		  "JEC uncertainty", "R=0.4 PF", "p_{T}=500 GeV",
		  10,10,"fixPt",500.);
  }

  if (_canvas) _canvas->SaveAs(Form("pdf/%s.pdf",sCorrGroups));
  _icanvas = 1;


  bool absUncertTmp = _absUncert;
  _absUncert = false;
  string sspt = Form("%s_AbsolutePt_%s",cu,names[jetAlg]);
  const char *spt = sspt.c_str();

  double ymaxpt = 2;//5; // 2

  // Absolute scale pT dependent uncertainty
  // vs pT
  if (_paper && !_absUncert) {
  plotUncertainty(sypt, 0, sypt.size(), jetAlg, Form("%s_Eta00",spt),
		  "JEC uncertainty", s,
		  "|#eta_{jet}| = 0", ymaxpt,10,"fixEta",0.);
  }
  if (_fourbytwo) {

  plotUncertainty(sypt, 0, sypt.size(), jetAlg, Form("%s_Eta00",spt),
		  "JEC uncertainty", s,
		  "|#eta_{jet}|=0", ymaxpt,10,"fixEta",0.);
  plotUncertainty(sypt, 0, sypt.size(), jetAlg, Form("%s_Eta20",spt),
		  "JEC uncertainty", s,
		  "|#eta_{jet}|=2.0", ymaxpt,10,"fixEta",2.0);
  plotUncertainty(sypt, 0, sypt.size(), jetAlg, Form("%s_Eta27",spt),
		  "JEC uncertainty", s,
		  "|#eta_{jet}|=2.7", ymaxpt,10,"fixEta",2.7);
  plotUncertainty(sypt, 0, sypt.size(), jetAlg, Form("%s_Eta42",spt),
		  "JEC uncertainty", s,
		  "|#eta_{jet}|=4.2", ymaxpt,10,"fixEta",4.2);
  // vs eta
  plotUncertainty(sypt, 0, sypt.size(), jetAlg, Form("%s_Pt30",spt),
		  "JEC uncertainty", s,
		  "p_{T}=30 GeV", ymaxpt,10,"fixPt",30.);
  plotUncertainty(sypt, 0, sypt.size(), jetAlg, Form("%s_Pt100",spt),
		  "JEC uncertainty", s,
		  "p_{T}=100 GeV", ymaxpt,10,"fixPt",100.);
  plotUncertainty(sypt, 0, sypt.size(), jetAlg, Form("%s_E500",spt),
		  "JEC uncertainty", s,
		  "E=500 GeV", ymaxpt,10,"fixE",500.);
  plotUncertainty(sypt, 0, sypt.size(), jetAlg, Form("%s_E2000",spt),
		  "JEC uncertainty", s,
		  "E=2000 GeV", ymaxpt,10,"fixE",2000.);
  plotUncertainty(sypt, 0, sypt.size(), jetAlg, Form("%s_E1000",spt),
		  "JEC uncertainty", s,
		  "E=1000 GeV", ymaxpt,10,"fixE",1000.); // moved last
  }
  if (_twobythree) {

  plotUncertainty(sypt, 0, sypt.size(), jetAlg, Form("%s_Eta00",spt),
		  "JEC uncertainty", s,
		  "|#eta_{jet}|=0", ymaxpt,10,"fixEta",0.);
  plotUncertainty(sypt, 0, sypt.size(), jetAlg, Form("%s_Pt30",spt),
		  "JEC uncertainty", s,
		  "p_{T}=30 GeV", ymaxpt,10,"fixPt",30.);
  //
  plotUncertainty(sypt, 0, sypt.size(), jetAlg, Form("%s_Eta27",spt),
		  "JEC uncertainty", s,
		  "|#eta_{jet}|=2.7", ymaxpt,10,"fixEta",2.7);
  plotUncertainty(sypt, 0, sypt.size(), jetAlg, Form("%s_Pt100",spt),
		  "JEC uncertainty", s,
		  "p_{T}=100 GeV", ymaxpt,10,"fixPt",100.);
  //
  plotUncertainty(sypt, 0, sypt.size(), jetAlg, Form("%s_Eta35",spt),
		  "JEC uncertainty", s,
		  "|#eta_{jet}|=3.5", ymaxpt,10,"fixEta",3.5);
  plotUncertainty(sypt, 0, sypt.size(), jetAlg, Form("%s_Pt500",spt),
		  "JEC uncertainty", s,
		  "p_{T}=500 GeV", ymaxpt,10,"fixPt",500.);
  }
  //
  _absUncert = absUncertTmp;

  if (_canvas) _canvas->SaveAs(Form("pdf/%s.pdf",spt));
  _icanvas = 1;

  string ssf = Form("%s_Flavor_%s",cu,names[jetAlg]);
  const char *sf = ssf.c_str();

  minimaltmp = _minimal;
  if (algo=="AK4PFchs") _minimal = false;

  // Flavor uncertainty
  // vs pT
  if (_paper && !_absUncert) {
  plotUncertainty(syf, 0, syf.size(), jetAlg, Form("%s_Eta00",sf),
                  "JEC uncertainty", s,
                  "|#eta_{jet}| = 0", 5,10,"fixEta",0.);
  plotUncertainty(syf, 0, syf.size(), jetAlg, Form("%s_Eta27",sf),
                  "JEC uncertainty", s,
                  "|#eta_{jet}| = 2.7", 5,10,"fixEta",2.7);
  plotUncertainty(syf, 0, syf.size(), jetAlg, Form("%s_Pt30",sf),
                  "JEC uncertainty", s,
                  "p_{T} = 30 GeV", 5,10,"fixPt",30.);
  //
  plotUncertainty(syf, 0, syf.size(), jetAlg, Form("%s_Pt100",sf),
                  "JEC uncertainty", s,
                  "p_{T} = 100 GeV", 5,10,"fixPt",100.);
  }
  if (_fourbytwo) {

  plotUncertainty(syf, 0, syf.size(), jetAlg, Form("%s_Eta00",sf),
                  "JEC uncertainty", s,
                  "|#eta_{jet}|=0", 5,10,"fixEta",0.);
  plotUncertainty(syf, 0, syf.size(), jetAlg, Form("%s_Eta20",sf),
                  "JEC uncertainty", s,
                  "|#eta_{jet}|=2.0", 5,10,"fixEta",2.0);
  plotUncertainty(syf, 0, syf.size(), jetAlg, Form("%s_Eta27",sf),
                  "JEC uncertainty", s,
                  "|#eta_{jet}|=2.7", 5,10,"fixEta",2.7);
  plotUncertainty(syf, 0, syf.size(), jetAlg, Form("%s_Eta42",sf),
                  "JEC uncertainty", s,
                  "|#eta_{jet}|=4.2", 5,10,"fixEta",4.2);
  // vs eta
  plotUncertainty(syf, 0, syf.size(), jetAlg, Form("%s_Pt30",sf),
                  "JEC uncertainty", s,
                  "p_{T}=30 GeV", 5,10,"fixPt",30.);
  plotUncertainty(syf, 0, syf.size(), jetAlg, Form("%s_Pt100",sf),
                  "JEC uncertainty", s,
                  "p_{T}=100 GeV", 5,10,"fixPt",100.);
  plotUncertainty(syf, 0, syf.size(), jetAlg, Form("%s_E500",sf),
                  "JEC uncertainty", s,
                  "E=500 GeV", 5,10,"fixE",500.);
  plotUncertainty(syf, 0, syf.size(), jetAlg, Form("%s_E2000",sf),
                  "JEC uncertainty", s,
                  "E=2000 GeV", 5,10,"fixE",2000.);
  plotUncertainty(syf, 0, syf.size(), jetAlg, Form("%s_E1000",sf),
                  "JEC uncertainty", s,
                  "E=1000 GeV", 5,10,"fixE",1000.); // moved last
  }
  if (_twobythree) {

  plotUncertainty(syf, 0, syf.size(), jetAlg, Form("%s_Eta00",sf),
                  "JEC uncertainty", s,
                  "|#eta_{jet}|=0", 5,10,"fixEta",0.);
  plotUncertainty(syf, 0, syf.size(), jetAlg, Form("%s_Pt30",sf),
                  "JEC uncertainty", s,
                  "p_{T}=30 GeV", 5,10,"fixPt",30.);
  //
  plotUncertainty(syf, 0, syf.size(), jetAlg, Form("%s_Eta27",sf),
                  "JEC uncertainty", s,
                  "|#eta_{jet}|=2.7", 5,10,"fixEta",2.7);
  plotUncertainty(syf, 0, syf.size(), jetAlg, Form("%s_Pt100",sf),
                  "JEC uncertainty", s,
                  "p_{T}=100 GeV", 5,10,"fixPt",100.);
  //
  plotUncertainty(syf, 0, syf.size(), jetAlg, Form("%s_Eta35",sf),
                  "JEC uncertainty", s,
                  "|#eta_{jet}|=3.5", 5,10,"fixEta",3.5);
  plotUncertainty(syf, 0, syf.size(), jetAlg, Form("%s_Pt500",sf),
                  "JEC uncertainty", s,
                  "p_{T}=500 GeV", 5,10,"fixPt",500.);
  }

  if (_canvas) _canvas->SaveAs(Form("pdf/%s.pdf",sf));
  _icanvas = 1;

  _minimal = minimaltmp;

  absUncertTmp = _absUncert;
  _absUncert = false;
  string sst = Form("%s_Time_%s",cu,names[jetAlg]);
  const char *st = sst.c_str();

  double ymaxt = 3.2;//6;//3.2;

  // Drop TimeEta for |eta|=0 plot
  //vector<uncert> syt0(syt.begin()+1,syt.end());

  // Time uncertainty
  // vs pT
  if (_paper && !_absUncert) {
  plotUncertainty(syt, 0, syt.size(), jetAlg, Form("%s_Eta00",st),
                  "JEC uncertainty", s,
                  "|#eta_{jet}| = 0", ymaxt,10,"fixEta",0.);
  plotUncertainty(syt, 0, syt.size(), jetAlg, Form("%s_Pt30",st),
                  "JEC uncertainty", s,
                  "p_{T} = 30 GeV", ymaxt,10,"fixPt",30.);
  plotUncertainty(syt, 0, syt.size(), jetAlg, Form("%s_Pt100",st),
                  "JEC uncertainty", s,
                  "p_{T} = 100 GeV", ymaxt,10,"fixPt",100.);
  plotUncertainty(syt, 0, syt.size(), jetAlg, Form("%s_Pt500",st),
                  "JEC uncertainty", s,
                  "p_{T} = 500 GeV", ymaxt,10,"fixPt",500.);
  }
  if (_fourbytwo) {

  plotUncertainty(syt, 0, syt.size(), jetAlg, Form("%s_Eta00",st),
                  "JEC uncertainty", s,
                  "|#eta_{jet}|=0", ymaxt,10,"fixEta",0.);
  plotUncertainty(syt, 0, syt.size(), jetAlg, Form("%s_Eta20",st),
                  "JEC uncertainty", s,
                  "|#eta_{jet}|=2.0", ymaxt,10,"fixEta",2.0);
  plotUncertainty(syt, 0, syt.size(), jetAlg, Form("%s_Eta27",st),
                  "JEC uncertainty", s,
                  "|#eta_{jet}|=2.7", ymaxt,10,"fixEta",2.7);
  plotUncertainty(syt, 0, syt.size(), jetAlg, Form("%s_Eta42",st),
                  "JEC uncertainty", s,
                  "|#eta_{jet}|=4.2", ymaxt,10,"fixEta",4.2);
  // vs eta
  plotUncertainty(syt, 0, syt.size(), jetAlg, Form("%s_Pt30",st),
                  "JEC uncertainty", s,
                  "p_{T}=30 GeV", ymaxt,10,"fixPt",30.);
  plotUncertainty(syt, 0, syt.size(), jetAlg, Form("%s_Pt100",st),
                  "JEC uncertainty", s,
                  "p_{T}=100 GeV", ymaxt,10,"fixPt",100.);
  plotUncertainty(syt, 0, syt.size(), jetAlg, Form("%s_E500",st),
                  "JEC uncertainty", s,
                  "E=500 GeV", ymaxt,10,"fixE",500.);
  plotUncertainty(syt, 0, syt.size(), jetAlg, Form("%s_E2000",st),
                  "JEC uncertainty", s,
                  "E=2000 GeV", ymaxt,10,"fixE",2000.);
  plotUncertainty(syt, 0, syt.size(), jetAlg, Form("%s_E1000",st),
                  "JEC uncertainty", s,
                  "E=1000 GeV", ymaxt,10,"fixE",1000.); // moved last
  }
  if (_twobythree) {
    
  plotUncertainty(syt, 0, syt.size(), jetAlg, Form("%s_Eta00",st),
                  "JEC uncertainty", s,
                  "|#eta_{jet}|=0", ymaxt,10,"fixEta",0.);
  plotUncertainty(syt, 0, syt.size(), jetAlg, Form("%s_Pt30",st),
                  "JEC uncertainty", s,
                  "p_{T}=30 GeV", ymaxt,10,"fixPt",30.);
  //
  plotUncertainty(syt, 0, syt.size(), jetAlg, Form("%s_Eta27",st),
                  "JEC uncertainty", s,
                  "|#eta_{jet}|=2.7", ymaxt,10,"fixEta",2.7);
  plotUncertainty(syt, 0, syt.size(), jetAlg, Form("%s_Pt100",st),
                  "JEC uncertainty", s,
                  "p_{T}=100 GeV", ymaxt,10,"fixPt",100.);
  //
  plotUncertainty(syt, 0, syt.size(), jetAlg, Form("%s_Eta35",st),
                  "JEC uncertainty", s,
                  "|#eta_{jet}|=3.5", ymaxt,10,"fixEta",3.5);
  plotUncertainty(syt, 0, syt.size(), jetAlg, Form("%s_Pt500",st),
                  "JEC uncertainty", s,
                  "p_{T}=500 GeV", ymaxt,10,"fixPt",500.);
  }
  //
  _absUncert = absUncertTmp;

  if (_canvas) _canvas->SaveAs(Form("pdf/%s.pdf",st));
  _icanvas = 1;

  absUncertTmp = _absUncert;
  _absUncert = false;
  string sss = Form("%s_SinglePion_%s",cu,names[jetAlg]);
  const char *css = sss.c_str();

  // Single pion versus eta dependence
  // vs pT
  if (_fourbytwo) {

  plotUncertainty(sys, 0, sys.size(), jetAlg, Form("%s_Eta00",css),
                  "JEC uncertainty", s,
                  "|#eta_{jet}|=0", -5,10,"fixEta",0.);
  plotUncertainty(sys, 0, sys.size(), jetAlg, Form("%s_Eta03",css),
                  "JEC uncertainty", s,
                  "|#eta_{jet}|=0.25", -5,10,"fixEta",0.25);
  plotUncertainty(sys, 0, sys.size(), jetAlg, Form("%s_Eta08",css),
                  "JEC uncertainty", s,
                  "|#eta_{jet}|=0.75", -5,10,"fixEta",0.75);
  plotUncertainty(sys, 0, sys.size(), jetAlg, Form("%s_Eta13",css),
                  "JEC uncertainty", s,
                  "|#eta_{jet}|=1.25", -5,10,"fixEta",1.25);
  plotUncertainty(sys, 0, sys.size(), jetAlg, Form("%s_Eta18",css),
                  "JEC uncertainty", s,
                  "|#eta_{jet}|=1.75", -5,10,"fixEta",1.75);
  plotUncertainty(sys, 0, sys.size(), jetAlg, Form("%s_Eta23",css),
                  "JEC uncertainty", s,
                  "|#eta_{jet}|=2.25", -5,10,"fixEta",2.25);
  plotUncertainty(sys, 0, sys.size(), jetAlg, Form("%s_Pt100",css),
                  "JEC uncertainty", s,
                  "p_{T}=100 GeV", -5,10,"fixPt",100.);
  plotUncertainty(sys, 0, sys.size(), jetAlg, Form("%s_E2000",css),
                "JEC uncertainty", s,
                "E=2000 GeV", -5,10,"fixE",2000.);
  }
  if (_twobythree) {

  plotUncertainty(sys, 0, sys.size(), jetAlg, Form("%s_Eta00",css),
                  "JEC uncertainty", s,
                  "|#eta_{jet}|=0", -5,10,"fixEta",0.);
  plotUncertainty(sys, 0, sys.size(), jetAlg, Form("%s_Eta03",css),
                  "JEC uncertainty", s,
                  "|#eta_{jet}|=0.25", -5,10,"fixEta",0.25);
  plotUncertainty(sys, 0, sys.size(), jetAlg, Form("%s_Eta08",css),
                  "JEC uncertainty", s,
                  "|#eta_{jet}|=0.75", -5,10,"fixEta",0.75);
  plotUncertainty(sys, 0, sys.size(), jetAlg, Form("%s_Eta13",css),
                  "JEC uncertainty", s,
                  "|#eta_{jet}|=1.25", -5,10,"fixEta",1.25);
  plotUncertainty(sys, 0, sys.size(), jetAlg, Form("%s_Eta18",css),
                  "JEC uncertainty", s,
                  "|#eta_{jet}|=1.75", -5,10,"fixEta",1.75);
  plotUncertainty(sys, 0, sys.size(), jetAlg, Form("%s_Eta23",css),
                  "JEC uncertainty", s,
                  "|#eta_{jet}|=2.25", -5,10,"fixEta",2.25);
  }

  // Additional plots for JEC plots
  if (_extra) {

    plotUncertainty(sys, 0, sys.size(), jetAlg, Form("%s_Pt070",css),
		    "JEC uncertainty", s,
		    "p_{T}=70 GeV", -5,10,"fixPt",70.);
    plotUncertainty(sys, 0, sys.size(), jetAlg, Form("%s_Pt100",css),
		    "JEC uncertainty", s,
		    "p_{T}=100 GeV", -5,10,"fixPt",100.);
    plotUncertainty(sys, 0, sys.size(), jetAlg, Form("%s_Pt600",css),
		    "JEC uncertainty", s,
		    "p_{T}=600 GeV", -5,10,"fixPt",600.);
  }

  _absUncert = absUncertTmp;

  if (_canvas) _canvas->SaveAs(Form("pdf/%s.pdf",css));
  _icanvas = 0;
  if (_canvas) delete _canvas;

} // L3Uncertainty_new

void plotUncertainty(vector<uncert> const& sys,
		     unsigned int nsys1, unsigned int nsys2,
		     jec::JetAlgo jetAlg,
		     string name, string ytitle,
		     string label1, string label2,
		     double emax, double ptmin,//) {
		     string type, double typevar) {

  cout << "plotUncertainty("<<name<<")"<<endl<<flush;

  // Create suitable binning
  const double x_pt[] =
    {8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
     97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 362, 430,
     507, 592, 686, 790, 905, 1032, 1172, 1327, 1497, 1684, 1890, //1999};
     //2000, 2238, 2500, 2787, 3103, 3450,
     2116, 2366, 2640, 2941, 3273, 3637, 
     4037, 4477, 4961, 5492, 6076, 7000};
  const int ndiv_pt = sizeof(x_pt)/sizeof(x_pt[0])-1;
  const double x_eta[] =
    {-5.4,-5.0,-4.4,-4,-3.5,-3,-2.8,-2.6,-2.4,-2.2,-2.0,
     -1.8,-1.6,-1.4,-1.2,-1.0, -0.8,-0.6,-0.4,-0.2,0.,
     0.2,0.4,0.6,0.8,1.0,1.2,1.4,
     1.6,1.8,2.0,2.2,2.4,2.6,2.8,3,3.5,4,4.4,5.0,5.4};
  const int ndiv_eta = sizeof(x_eta)/sizeof(x_eta[0])-1;     

  // Re-determine bin edges based on maxe=4000.
  int jx(0), ndiv_new(0);
  if (type=="fixEta") {
    double maxe = 6500;//4000;
    double minpt = 10.;
    double maxpt = maxe/cosh(typevar);
    int imin(0), imax(ndiv_pt+1);
    for (int i = 0; i != ndiv_pt; ++i) {
      if (x_pt[i]<=minpt) imin = i;
      if (x_pt[i]>maxpt && imax==ndiv_pt+1) imax = i+1; // is ok?
    }
    jx = imin;
    ndiv_new = imax-imin-1; // is ok?
  }
  if (type=="fixPt") {
    double maxe = 6500;//4000;
    double maxeta = TMath::ACosH(maxe/typevar);
    int imin(0), imax(ndiv_eta+1);
    for (int i = 0; i != ndiv_eta; ++i) {
      if (x_eta[i]<-maxeta) imin = i;
      if (x_eta[i]>maxeta && imax==ndiv_eta+1) imax = i+1; // is ok?
    }
    jx = imin;
    ndiv_new = imax-imin-1; // is ok?
  }
  if (type=="fixE") {
    double minpt = 10;
    double maxeta = TMath::ACosH(typevar/minpt);
    int imin(0), imax(ndiv_eta+1);
    for (int i = 0; i != ndiv_eta; ++i) {
      if (x_eta[i]<-maxeta) imin = i;
      if (x_eta[i]>maxeta && imax==ndiv_eta+1) imax = i+1; // is ok?
    }
    jx = imin;
    ndiv_new = imax-imin-1; // is ok?
  }

  //cout << "Got here 1" << endl << flush;
  
  const double *x(0), *x0(0);
  int ndiv0(0);
  assert(type=="fixEta" || type=="fixPt" || type=="fixE");
  double eta(0), pt(0), e(0), ndiv(0);
  if (type=="fixEta") { 
    eta = typevar;
    x = &x_pt[jx]; ndiv = ndiv_new;//-1;
    x0 = &x_pt[1]; ndiv0 = ndiv_pt-1;
  }
  if (type=="fixPt") {
    pt = typevar;
    x = &x_eta[jx]; ndiv = ndiv_new;
    x0 = &x_eta[0]; ndiv0 = ndiv_eta;
  }
  if (type=="fixE") {
    e = typevar;
    x = &x_eta[jx]; ndiv = ndiv_new;
    x0 = &x_eta[0]; ndiv0 = ndiv_eta;
  }

  //cout << "Got here 2" << endl << flush;
  
  TVirtualPad *c2(0);
  if (_canvas) {
    _canvas->cd(_icanvas++);
    c2 = gPad;
  }

  if (c2) { if (type=="fixEta") c2->SetLogx(); }

  const char *cy = ytitle.c_str();
  TH1D *h0 = new TH1D(Form("h0_%s_%d",name.c_str(),++icnt),
		      Form(";p_{T} (GeV);%s (%%)",cy),
		      ndiv0, &x0[0]);
  if (type=="fixEta") h0->GetXaxis()->SetTitle("p_{T} (GeV)");
  if (type=="fixPt") h0->GetXaxis()->SetTitle("#eta_{jet}");
  if (type=="fixE") h0->GetXaxis()->SetTitle("#eta_{jet}");
  h0->SetMinimum(_absUncert ? 0 : -0.75*fabs(emax));
  h0->SetMaximum(_absUncert ? fabs(emax) : 1.5*fabs(emax));
  h0->GetXaxis()->SetMoreLogLabels();
  h0->GetXaxis()->SetNoExponent();
  h0->GetYaxis()->SetTitleOffset(1.0);

  if (_paper) h0->GetXaxis()->SetRangeUser(10,3500);
  //lumi_13TeV = "Run2016BCDEFGH re-reco, 36.5 fb^{-1}"; // Sum16
  lumi_13TeV = "2017BCDEF re-reco, 41.4 fb^{-1}"; // Fall17
  TCanvas *c1 = tdrCanvas(Form("c1_%s",name.c_str()),h0,4,11,kSquare);
  if (type=="fixEta") c1->SetLogx();

  if (c2) { c2->cd(); h0->DrawClone("AXIS"); c1->cd(); }


  //cout << "Got here 3" << endl << flush;

  TLegend *leg1 = new TLegend(0.18,0.90-0.05*nsys1,0.38,0.90,"","brNDC");
  leg1->SetFillStyle(kNone);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.045);
  const double tx = (type=="fixEta" ? 0.50 : 0.50);
  TLegend *leg2 = new TLegend(tx,0.90-0.05*nsys2,tx+0.20,0.90,"","brNDC");
  leg2->SetFillStyle(kNone);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.045);
  TLatex *tex1 = new TLatex(0.185, 0.73, label1.c_str());
  tex1->SetNDC();
  tex1->SetTextSize(0.045);
  TLatex *tex2 = new TLatex(0.185, 0.68, label2.c_str());
  tex2->SetNDC();
  tex2->SetTextSize(0.045);

  //cout << "Got here 4" << endl << flush;

  // Create and draw uncertainties
  const unsigned int nsys = nsys1 + nsys2;
  assert(sys.size()>=nsys);
  //cout << "sys.size(): " << sys.size() << " nsys: " << nsys << endl << flush;
  for (unsigned int isys = 0; isys != nsys; ++isys) {

    //cout << "...isys: " << isys << endl << flush;
    uncert const& u = sys[isys];
    double mu = (u.mu==-1 ? d_mu : u.mu);
    JECUncertainty rjet(jetAlg, jec::DATA, u.type, mu);
    TGraph *g = new TGraph(0);
    g->SetName(Form("L3_%s",u.name.c_str()));
    TGraph *gup = new TGraph(0);
    gup->SetName(Form("L3_%s_up",u.name.c_str()));
    TGraph *gdw = new TGraph(0);
    gdw->SetName(Form("L3_%s_dw",u.name.c_str()));
    TH1D *h = new TH1D(Form("L3_%s_%s_%d",name.c_str(),u.name.c_str(),++icnt),
		       "",ndiv,&x[0]);

    //cout << "Loop bins" << endl << flush;
    for (int ix = 0; ix != ndiv; ++ix) {
      //cout << "." << flush;
      double var = 0.5*(x[ix]+x[ix+1]);
      double pt(0), eta(0);
      if (type=="fixEta") { pt = var; eta = typevar; }
      if (type=="fixPt") { eta = var; pt = typevar; }
      if (type=="fixE") { eta = var; pt = typevar/cosh(eta); }
      if (pt > ptmin) { // DP note
	double err(0);
	double r = 1;//rjet.Rjet(pt, eta, err);
	err = rjet.Uncert(pt, eta);
	if (_absUncert) err = fabs(err);
	int n = g->GetN();
	g->SetPoint(n, var, 100. * err / r);
	gup->SetPoint(n, var, +100. * err / r);
	gdw->SetPoint(n, var, -100. * err / r);
	h->SetBinError(ix+1, 100. * err / r);
	h->SetBinContent(ix+1, 0);
      }
    } // for ix

    //cout << "u.opt="<<u.opt<<" / u.option="<<u.option << endl << flush;
    assert(u.option=="L" || u.option=="LF" || u.option=="LP" ||
	   u.option=="LFP");
    TString *opt2 = new TString(u.option.c_str());
    if (opt2->Contains("F")) {
      //cout << "x" << flush;
      h->SetMarkerStyle(kNone);
      h->SetLineStyle(kNone);
      h->SetFillStyle(u.fstyle);
      h->SetFillColor(u.fcolor);
      TH1D *hc = (TH1D*)h->DrawClone("SAME E3");
      // Make Run I transparent
      if (u.title=="Run I") {
	hc->SetFillColorAlpha(u.fcolor, 0.70); // 70% transparent (def 35%)
      }
      if (c2) { c2->cd(); h->DrawClone("SAME E3"); c1->cd(); }
      opt2->ReplaceAll("F","");
    }
    g->SetLineColor(u.lcolor);
    g->SetLineStyle(u.lstyle);
    g->SetLineWidth(u.lwid);
    g->SetMarkerColor(u.mcolor);
    g->SetMarkerStyle(u.mstyle);
    g->SetFillColor(u.fcolor);
    g->SetFillStyle(u.fstyle);
    g->DrawClone(*opt2);
    if (c2) { c2->cd(); g->DrawClone(*opt2); c1->cd(); }

    if (isys<nsys1) leg1->AddEntry(g, u.title.c_str(), u.option.c_str());
    else            leg2->AddEntry(g, u.title.c_str(), u.option.c_str());

  } // for isys

  //cout << "Got here 5" << endl << flush;

  // Hide low pT stuff
  if (type=="fixEta") {
    TBox *box = new TBox(h0->GetXaxis()->GetXmin(), h0->GetMinimum(),
			 ptmin, h0->GetMaximum());
    box->SetFillStyle(1001);
    box->SetFillColor(kWhite);
    box->Draw();
    if (c2) { c2->cd(); box->DrawClone(); c1->cd(); }
  }

  leg1->Draw();
  leg2->Draw();
  tex1->Draw();  
  tex2->Draw();  
  if (c2) { c2->cd(); leg1->DrawClone(); leg2->DrawClone(); tex1->DrawClone(); tex2->DrawClone(); c1->cd(); }

  gPad->RedrawAxis();
  if (c2) {
    c2->cd();
    gPad->RedrawAxis();
    c1->cd();
  }

  if (!_minimal) {
    if (_eps) c1->SaveAs(("eps/"+name+".eps").c_str());
    if (_pdf) c1->SaveAs(("pdf/"+name+".pdf").c_str());
  }    

  //cout << "Got here 6" << endl << flush;

  if (_doTXT && !_didTXT) {

    JECUncertainty rjets(d_algo, jec::DATA, jec::kData, d_mu);
    //ofstream fouts(Form("txt/Summer16_03Feb2017_V9_DATA_Uncertainty_%s.txt",
    ofstream fouts(Form("txt/Fall17_07Nov2017_V6M_DATA_Uncertainty_%s.txt",
			(*_algnames)[d_algo]), ios::out);
    fouts << "{1 JetEta 1 JetPt \"\" Correction Uncertainty}" << endl;

    for (int ieta = 0; ieta != ndiv_eta; ++ieta) {
      
      double etamin = x_eta[ieta];
      double etamax = x_eta[ieta+1];
      double eta = 0.5*(etamin+etamax);

      fouts << Form("%1.1f %1.1f %d ",etamin,etamax,ndiv_pt*3);
      
      for (int ipt = 0; ipt != ndiv_pt; ++ipt) {
	
	double pt = 0.5*(x_pt[ipt]+x_pt[ipt+1]);

	double errs = rjets.Uncert(pt, eta);
	fouts << Form("%1.1f %1.4f %1.4f ", pt, errs, errs);
      } // for ipt

      fouts << endl;
    } // for ieta

    //ofstream fout(Form("txt/Summer16_03Feb2017_V9_DATA_UncertaintySources_%s.txt",
    ofstream fout(Form("txt/Fall17_07Nov2017_V6M_DATA_UncertaintySources_%s.txt",
		       (*_algnames)[d_algo]), ios::out);
    //fout << Form("#Uncertainty sources for Summer16_03Feb2017_V9_DATA_%s",
    fout << Form("#Uncertainty sources for Fall17_07Nov2017_V6M_DATA_%s",
		 (*_algnames)[d_algo]) << endl;
    cout << "Storing uncertainties to: "
      //<< Form("txt/Summer16_03Feb2017_V9_DATA_UncertaintySources_%s.txt",
	 << Form("txt/Fall17_07Nov2017_V6M_DATA_UncertaintySources_%s.txt",
		 (*_algnames)[d_algo]) << endl;

    jec::ErrorTypes vsrc[] =
      {jec::kAbsoluteStat, jec::kAbsoluteScale, jec::kAbsoluteFlavorMapping, jec::kAbsoluteMPFBias,
       jec::kAbsoluteFrag,
       jec::kAbsoluteSPRE, jec::kAbsoluteSPRH,
       jec::kFlavorQCD,
       jec::kTimePtEta,
       jec::kRelativeJEREC1, jec::kRelativeJEREC2, jec::kRelativeJERHF,
       jec::kRelativePtBB,
       jec::kRelativePtEC1, jec::kRelativePtEC2, jec::kRelativePtHF,
       jec::kRelativeBal, jec::kRelativeSample,
       jec::kRelativeFSR, jec::kRelativeStatFSR, 
       jec::kRelativeStatEC, jec::kRelativeStatHF,
       jec::kPileUpDataMC, jec::kPileUpPtRef,
       jec::kPileUpPtBB, jec::kPileUpPtEC1, jec::kPileUpPtEC2, jec::kPileUpPtHF,
       jec::kPileUpMuZero,
       jec::kPileUpEnvelope,
       jec::kPileUp, jec::kRelative, jec::kAbsolutePt, jec::kAbsoluteFlat,
       jec::kAbsolute,
       jec::kMC,
       jec::kData,
       jec::kDataNoFlavor, jec::kDataNoTime, jec::kDataNoFlavorNoTime,
       jec::kFlavorZJet, jec::kFlavorPhotonJet,
       jec::kFlavorPureGluon, jec::kFlavorPureQuark, 
       jec::kFlavorPureCharm, jec::kFlavorPureBottom,
       jec::kTimeRunB, jec::kTimeRunC, jec::kTimeRunD, jec::kTimeRunE,
       jec::kTimeRunF,
       jec::kCorrelationGroupMPFInSitu, jec::kCorrelationGroupIntercalibration, jec::kCorrelationGroupbJES,
       jec::kCorrelationGroupFlavor, jec::kCorrelationGroupUncorrelated
      };

    const int nsrc = sizeof(vsrc)/sizeof(vsrc[0]);
    map<jec::ErrorTypes, string> srcname;
    srcname[jec::kAbsoluteStat] =           "AbsoluteStat";    // L3Res fit
    srcname[jec::kAbsoluteScale] =          "AbsoluteScale";   // Zll scale
    srcname[jec::kAbsoluteFlavorMapping] =  "AbsoluteFlavMap"; // obsolete
    srcname[jec::kAbsoluteMPFBias] =        "AbsoluteMPFBias"; // old, 0.2%
    srcname[jec::kAbsoluteFrag] =           "Fragmentation";   // update?
    srcname[jec::kAbsoluteSPRE] =           "SinglePionECAL";  // +/-3%
    srcname[jec::kAbsoluteSPRH] =           "SinglePionHCAL";  // L3Res fit
    srcname[jec::kRelativeJEREC1] =  "RelativeJEREC1";
    srcname[jec::kRelativeJEREC2] =  "RelativeJEREC2";
    srcname[jec::kRelativeJERHF] =   "RelativeJERHF";
    srcname[jec::kRelativePtBB] =    "RelativePtBB";   // flat vs log-lin
    srcname[jec::kRelativePtEC1] =   "RelativePtEC1";
    srcname[jec::kRelativePtEC2] =   "RelativePtEC2";
    srcname[jec::kRelativePtHF] =    "RelativePtHF";
    srcname[jec::kRelativeBal] =     "RelativeBal";    // pT balance vs MPf
    srcname[jec::kRelativeSample] =  "RelativeSample"; // dijet vs Z+jet
    srcname[jec::kRelativeFSR] =     "RelativeFSR";    // Pythia vs Herwig
    srcname[jec::kRelativeStatFSR] = "RelativeStatFSR";
    srcname[jec::kRelativeStatEC] =  "RelativeStatEC";
    srcname[jec::kRelativeStatHF] =  "RelativeStatHF";
    srcname[jec::kPileUpDataMC] =   "PileUpDataMC";
    srcname[jec::kPileUpPtRef] =    "PileUpPtRef";
    srcname[jec::kPileUpPtBB] =     "PileUpPtBB";
    srcname[jec::kPileUpPtEC1] =    "PileUpPtEC1";
    srcname[jec::kPileUpPtEC2] =    "PileUpPtEC2";
    srcname[jec::kPileUpPtHF] =     "PileUpPtHF";
    srcname[jec::kPileUpMuZero] =   "PileUpMuZero";
    srcname[jec::kPileUpEnvelope] = "PileUpEnvelope"; // L1Data-L1RC
    srcname[jec::kFlavorQCD] = "FlavorQCD";
    srcname[jec::kTimePtEta] = "TimePtEta";
    srcname[jec::kPileUp] = "SubTotalPileUp";
    srcname[jec::kRelative] = "SubTotalRelative";
    srcname[jec::kAbsolutePt] = "SubTotalPt";
    srcname[jec::kAbsoluteFlat] = "SubTotalScale";
    srcname[jec::kAbsolute] = "SubTotalAbsolute";
    srcname[jec::kMC] = "SubTotalMC";
    srcname[jec::kData] = "Total";
    srcname[jec::kDataNoFlavor] = "TotalNoFlavor";
    srcname[jec::kDataNoTime] = "TotalNoTime";
    srcname[jec::kDataNoFlavorNoTime] = "TotalNoFlavorNoTime";
    srcname[jec::kFlavorZJet] = "FlavorZJet";
    srcname[jec::kFlavorPhotonJet] = "FlavorPhotonJet";
    srcname[jec::kFlavorPureGluon] = "FlavorPureGluon";
    srcname[jec::kFlavorPureQuark] = "FlavorPureQuark";
    srcname[jec::kFlavorPureBottom] = "FlavorPureBottom";
    srcname[jec::kFlavorPureCharm] = "FlavorPureCharm";
    srcname[jec::kTimeRunB] =  "TimeRunB"; 
    srcname[jec::kTimeRunC] =  "TimeRunC";
    srcname[jec::kTimeRunD] = "TimeRunD";
    srcname[jec::kTimeRunE] = "TimeRunE";
    srcname[jec::kTimeRunF] =  "TimeRunF";
    srcname[jec::kCorrelationGroupMPFInSitu] = "CorrelationGroupMPFInSitu";
    srcname[jec::kCorrelationGroupFlavor] = "CorrelationGroupFlavor";
    srcname[jec::kCorrelationGroupIntercalibration] = "CorrelationGroupIntercalibration";
    srcname[jec::kCorrelationGroupbJES] = "CorrelationGroupbJES";
    srcname[jec::kCorrelationGroupUncorrelated] = "CorrelationGroupUncorrelated";

    if (!(srcname.size()==nsrc)) {
      cout << "srcname.size()="<<srcname.size()
	   << " nsrc="<<nsrc<<endl<<flush;
    }
    assert(srcname.size()==nsrc); // check that '<' defined for ErrorTypes

    for (int isrc = 0; isrc != nsrc; ++isrc) {

      jec::ErrorTypes &src = vsrc[isrc];
      std::cout << srcname[src] << "\", \"" << std::endl;

      JECUncertainty rjet(d_algo, jec::DATA, src, d_mu);
      fout << "["<<srcname[src]<<"]" << endl;
      fout << "{1 JetEta 1 JetPt \"\" Correction JECSource}" << endl;

      for (int ieta = 0; ieta != ndiv_eta; ++ieta) {

	double etamin = x_eta[ieta];
	double etamax = x_eta[ieta+1];
	double eta = 0.5*(etamin+etamax);
	fout << Form("%1.1f %1.1f %d ",etamin,etamax,ndiv_pt*3);

	for (int ipt = 0; ipt != ndiv_pt; ++ipt) {
	
	  double pt = 0.5*(x_pt[ipt]+x_pt[ipt+1]);

	  double err(0);
	  double r = 1.;//rjet.Rjet(pt, eta, err);
	  err = rjet.Uncert(pt, eta);
	  err /= r; // relative uncertainty
	  fout << Form("%1.1f %1.4f %1.4f ", pt, err, err);
	} // for ipt

	fout << endl;
      } // for ieta
    } // for isrc
    
    _didTXT = true;
  } // print uncertainty sources

} // plotUncertainty

