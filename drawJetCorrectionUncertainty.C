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

#include "tdrstyle_mod.C"

#include <string>
#include <fstream>

using namespace std;

// Don't plot individual bins, just keep 4x2
bool _minimal = true;

// Print uncertainty (true) or source (false)
bool _absUncert = true;
// NB: All source files are currently printed together with AK5PF uncertainty

// List of (hard-coded) default parameters
jec::JetAlgo  d_algo = jec::AK5PFchs; // Replaced in function call
const jec::DataType d_type = jec::DATA; // Uncertainties for data (or data/MC)
const double d_npv = 14; // Average pile-up for L1 uncertainties
const bool d_mpf = true; // L2L3Res uncertainties for MPF method

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

  /*
    TLine *l1 = new TLine();
    l1->SetLineStyle(1);
    l1->SetLineWidth(25);
    l1->SetLineColor(h->GetFillColor());
    l1->DrawLineNDC(xmin, ymin+0.01, xmin+0.07, ymin+0.01);
  */
  
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
  double npv; // average pv
  int lcolor;
  int lstyle;
  int lwid;
  int mcolor;
  int mstyle;
  int fcolor;
  int fstyle;
  string option;
  const char* opt;
  uncert(string name_, string title_, jec::ErrorTypes type_,
	 string method_, string jetalgo_, double npv_,
	 int lcolor_, int lstyle_, int lwid_ = 1,
	 int mcolor_ = kBlack, int mstyle_ = kNone,
	 int fcolor_ = kBlack, int fstyle_ = kNone,
	 string option_ = "L") :
    name(name_), title(title_), type(type_),
    method(method_), jetalgo(jetalgo_), npv(npv_),
    lcolor(lcolor_), lstyle(lstyle_), lwid(lwid_),
    mcolor(mcolor_), mstyle(mstyle_),
    fcolor(fcolor_), fstyle(fstyle_),
    option(option_) {

    opt = option.c_str();
    assert(method=="default" || method=="mpf" || method=="pt");
    assert(npv>=0 || npv==-1);
  }
};

void plotUncertainty(vector<uncert> const& sys,
		     unsigned int nsys1, unsigned int nsys2,
		     jec::JetAlgo jetAlg,
		     string name, string ytitle,
		     string label1, string label2,
		     double emax, double ptmin,//);//, bool plotLog);
		     string type="fixPt", double typevar=0.);

void drawJetCorrectionUncertainty(string algo = "AK5PF") {
  
  if (algo=="AK5PF") d_algo = jec::AK5PF;
  if (algo=="AK5PFchs") d_algo = jec::AK5PFchs;
  if (algo=="AK7PF") d_algo = jec::AK7PF;
  if (algo=="AK7PFchs") d_algo = jec::AK7PFchs;
  if (algo=="AK5CALO") d_algo = jec::AK5CALO;
  if (algo=="AK7CALO") d_algo = jec::AK7CALO;

  cout << "drawJetCorrectionUncertainty" << endl << flush;

  vector<uncert> sy;
  sy.push_back(uncert("tot", "Total uncertainty", jec::kData,
		      "default", "default", -1, // defaults
		      kBlack, kSolid, 1, // line
		      kBlack, kNone, // marker
		      kDarkGray, 1001, "LF")); // fill
  sy.push_back(uncert("absolute", "Absolute scale", jec::kAbsoluteScale,
		      "default", "default", -1, // defaults
		      kYellow+3, kSolid, 1, // line
		      kBlack, kNone, // marker
		      kYellow, 1001, "LF")); // fill
  sy.push_back(uncert("relative", "Relative scale",jec::kRelative,
		      "default", "default", -1, // defaults
		      kBlack, kSolid, 1, // line
		      kBlack, kFullTriangleDown, // marker
		      kNone, kNone, "LP")); // fill
  sy.push_back(uncert("highpt", "Extrapolation", jec::kAbsolutePt,//kPtExtra,
		      "default", "default", -1, // defaults
		      kRed, kSolid, 1, // line
		      kRed, kOpenCircle, // marker
		      kNone, kNone, "LP")); // fill
  sy.push_back(uncert("pileup", Form("Pile-up, NPV=%1.0f",d_npv), jec::kPileUp,
		      "default", "default", -1, // defaults
		      kBlue, kNone, 1, // line
		      kBlue, kFullSquare, // marker
		      kNone, kNone, "LP")); // fill
  //sy.push_back(uncert("flavjec", "Jet flavor", jec::kFlavor,
  sy.push_back(uncert("flavjec", "Jet flavor (QCD)", jec::kFlavorQCD,
		      "default", "default", -1, // defaults
		      kGreen+2, kSolid, 1, // line
		      kGreen+2, kOpenSquare, // marker
		      kNone, kNone, "LP")); // fill
  /*
  sy.push_back(uncert("residual", "Residuals", jec::kResidual,
		      "default", "default", -1, // defaults
		      kBlack, kNone, 1, // line
		      kBlack, kOpenTriangleUp, // marker
		      kNone, kNone, "LP")); // fill
  */
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
  sym.push_back(uncert("absolute", "Absolute scale", jec::kAbsoluteScale,
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
  sym.push_back(uncert("pileup", Form("Pile-up, NPV=%1.0f",d_npv),
		       jec::kPileUpDataMC,
		       "default", "default", -1, // defaults
		       kBlue, kNone, 1, // line
		       kBlue, kFullSquare, // marker
		       kNone, kNone, "LP")); // fill
  //sym.push_back(uncert("flavjec", "Jet flavor", jec::kFlavorMC,
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
  sypu.push_back(uncert("pupt", "PileUpPt",
			jec::kPileUpPtBB | jec::kPileUpPtEC | jec::kPileUpPtHF,
			"default", "default", -1, // defaults
			kYellow+3, kNone, 1, // line
			kBlack, kNone, // marker
			kYellow, 1001, "LF")); // fill
  sypu.push_back(uncert("pubias", "PileUpBias (off)",
			jec::kPileUpBias,
			"default", "default", -1, // defaults
			kRed, kNone, 1, // line
			kRed, kOpenCircle, // marker
			kNone, kNone, "LP")); // fill
  sypu.push_back(uncert("pudatamc", "PileUpDataMC",
			jec::kPileUpDataMC,
			"default", "default", -1, // defaults
			kBlack, kNone, 1, // line
			kBlack, kFullTriangleDown, // marker
			kNone, kNone, "LP")); // fill

  vector<uncert> syrel;
  syrel.push_back(uncert("relative", "SubTotalRelative",
			 jec::kRelative,
			 "default", "default", -1, // defaults
			 kBlack, kNone, 1, // line
			 kBlack, kFullTriangleDown, // marker
			 kBlue-9, 1001, "LFP")); // fill
  syrel.push_back(uncert("relpt", "RelativePt",
			 jec::kRelativePtBB | jec::kRelativePtEC1 |
			 jec::kRelativePtEC2 | jec::kRelativePtHF,
			 "default", "default", -1, // defaults
			 kYellow+3, kNone, 1, // line
			 kBlack, kNone, // marker
			 kYellow, 1001, "LF")); // fill
  syrel.push_back(uncert("reljer", "RelativeJER",
			 jec::kRelativeJEREC1 | jec::kRelativeJEREC2 |
			 jec::kRelativeJERHF,
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
			 jec::kRelativeStatEC2 | jec::kRelativeStatHF,
			 "default", "default", -1, // defaults
			 kGreen+2, kNone, 1, // line
			 kGreen+2, kOpenSquare, // marker
			 kBlack, kNone, "LP")); // fill
  /*
  syrel.push_back(uncert("relsample", "RelativeSample (off)",
			 jec::kRelativeSample,
			 "default", "default", -1, // defaults
			 kMagenta+2, kNone, 1, // line
			 kMagenta+2, kOpenTriangleDown, // marker
			 kBlack, kNone, "LP")); // fill
  */

  vector<uncert> sypt;
  sypt.push_back(uncert("abspt", "SubTotalPt",
			jec::kAbsolutePt,//kPtExtra,
			"default", "default", -1, // defaults
			kRed, kNone, 1, // line
			kRed, kOpenCircle, // marker
			kRed-9, 1001, "LFP")); // fill
  sypt.push_back(uncert("absfrag", "HighPtExtra",
			jec::kAbsoluteFrag,
			"default", "default", -1, // defaults
			kYellow+3, kNone, 1, // line
			kBlack, kNone, // marker
			kYellow, 1001, "LF")); // fill
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
  sypt.push_back(uncert("absecal", "ECAL (off)",
			jec::kAbsoluteECAL,
			"default", "default", -1, // defaults
			kGreen+2, kNone, 1, // line
			kGreen+2, kOpenSquare, // marker
			kNone, kNone, "LP")); // fill
  sypt.push_back(uncert("abstrack", "Tracker (off)",
			jec::kAbsoluteTrack,
			"default", "default", -1, // defaults
			kMagenta+2, kNone, 1, // line
			kMagenta+2, kOpenTriangleDown, // marker
			kNone, kNone, "LP")); // fill

  vector<uncert> syf;
  syf.push_back(uncert("flavor_gluon", "Gluon",
                       jec::kFlavorPureGluon,
                      "default", "default", -1, // defaults
                      kYellow+3, kSolid, 1, // line
                      kBlack, kNone, // marker
                      kYellow, 1001, "LF")); // fill
  syf.push_back(uncert("flavor_qcd", "QCD Mixture",
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
  syf.push_back(uncert("flavor", "Bottom",
                       jec::kFlavorPureBottom,
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
  sys.push_back(uncert("relativept", "RelativePt",
                       jec::kRelativePt,
		       "default", "default", -1, // defaults
		       kBlack, kSolid, 1, // line
		       kBlack, kFullCircle, // marker
		       kNone, kNone, "LP")); // fill

  vector<uncert> syCorrGroups;
  syCorrGroups.push_back(uncert("tot", "Total uncertainty", jec::kMC,
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

  double r = (jetAlg==jec::AK7PF||jetAlg==jec::AK7PFchs||jetAlg==jec::AK7CALO ?
	      0.7 : 0.5);
  string sa = ((jetAlg==jec::AK5CALO||jetAlg==jec::AK7CALO) ?
	       "Calo" : ((jetAlg==jec::AK5PFchs||jetAlg==jec::AK7PFchs)
			 ? "PFchs" : "PF"));
  string ss = Form("Anti-k_{T} R=%1.1f %s", r, sa.c_str());
  const char *s = ss.c_str();

  const char *cu = (_absUncert ? "JECUncert" : "JECSource");

  map<jec::JetAlgo, const char*> names;
  names[jec::AK5PF] = "AK5PF";
  names[jec::AK5PFchs] = "AK5PFchs";
  names[jec::AK7PF] = "AK7PF";
  names[jec::AK7PFchs] = "AK7PFchs";
  names[jec::AK5CALO] = "AK5CALO";
  names[jec::AK7CALO] = "AK7CALO";

  string ssd = Form("%s_DATA_Summary_%s", cu, names[jetAlg]);
  //d_algo==jec::AK5PF ? "AK5PF" :
  //(d_algo==jec::AK5PFchs ? "AK5PFchs" : "OTHER"));
  const char *sd = ssd.c_str();

  TCanvas *c = new TCanvas("c4x2","c4x2",2400,1200);
  c->Divide(4,2);
  _canvas = c;
  _icanvas = 1;

  bool minimaltmp = _minimal;
  //if (algo=="AK5PFchs") _minimal = false;
  if (algo=="AK5PF") _minimal = false;

  // Data uncertainty
  // vs pT
  plotUncertainty(sy, 0, sy.size(), jetAlg, Form("%s_Eta00",sd),
		  "JEC uncertainty", s, "|#eta_{jet}|=0", 10,10,"fixEta",0.);
  plotUncertainty(sy, 0, sy.size(), jetAlg, Form("%s_Eta20",sd),
		  "JEC uncertainty", s, "|#eta_{jet}|=2.0", 10,10,"fixEta",2.0);
  plotUncertainty(sy, 0, sy.size(), jetAlg, Form("%s_Eta27",sd),
		  "JEC uncertainty", s, "|#eta_{jet}|=2.7", 10,10,"fixEta",2.7);
  plotUncertainty(sy, 0, sy.size(), jetAlg, Form("%s_Eta42",sd),
		  "JEC uncertainty", s, "|#eta_{jet}|=4.2", 10,10,"fixEta",4.2);
  // vs eta
  plotUncertainty(sy, 0, sy.size(), jetAlg, Form("%s_Pt30",sd),
		  "JEC uncertainty", s, "p_{T}=30 GeV", 10,10,"fixPt",30.);
  plotUncertainty(sy, 0, sy.size(), jetAlg, Form("%s_Pt100",sd),
		  "JEC uncertainty", s, "p_{T}=100 GeV", 10,10,"fixPt",100.);
  plotUncertainty(sy, 0, sy.size(), jetAlg, Form("%s_E500",sd),
		  "JEC uncertainty", s, "E=500 GeV", 10,10,"fixE",500.);
  plotUncertainty(sy, 0, sy.size(), jetAlg, Form("%s_E2000",sd),
		  "JEC uncertainty", s, "E=2000 GeV", 10,10,"fixE",2000.);
  plotUncertainty(sy, 0, sy.size(), jetAlg, Form("%s_E1000",sd),
		  "JEC uncertainty", s, "E=1000 GeV", 10,10,"fixE",1000.);//last

  _canvas->SaveAs(Form("pdf/%s.pdf",sd));
  _icanvas = 1;

  _minimal = minimaltmp;

  string ssm = Form("%s_MC_Summary_%s",cu,names[jetAlg]);
  //		    d_algo==jec::AK5PF ? "AK5PF" :
  //		    (d_algo==jec::AK5PFchs ? "AK5PFchs" : "OTHER"));
  const char *sm = ssm.c_str();

  // Data/MC (MC) uncertainty
  // vs pT
  plotUncertainty(sym, 0, sym.size(), jetAlg, Form("%s_Eta00",sm),
		  "JEC uncertainty (Data/MC)", s,
		  "|#eta_{jet}|=0", 10,10,"fixEta",0.);
  plotUncertainty(sym, 0, sym.size(), jetAlg, Form("%s_Eta20",sm),
		  "JEC uncertainty (Data/MC)", s,
		  "|#eta_{jet}|=2.0", 10,10,"fixEta",2.0);
  plotUncertainty(sym, 0, sym.size(), jetAlg, Form("%s_Eta27",sm),
		  "JEC uncertainty (Data/MC)", s,
		  "|#eta_{jet}|=2.7", 10,10,"fixEta",2.7);
  plotUncertainty(sym, 0, sym.size(), jetAlg, Form("%s_Eta42",sm),
		  "JEC uncertainty (Data/MC)", s,
		  "|#eta_{jet}|=4.2", 10,10,"fixEta",4.2);
  // vs eta
  plotUncertainty(sym, 0, sym.size(), jetAlg, Form("%s_Pt30",sm),
		  "JEC uncertainty (Data/MC)", s,
		  "p_{T}=30 GeV", 10,10,"fixPt",30.);
  plotUncertainty(sym, 0, sym.size(), jetAlg, Form("%s_Pt100",sm),
		  "JEC uncertainty (Data/MC)", s,
		  "p_{T}=100 GeV", 10,10,"fixPt",100.);
  plotUncertainty(sym, 0, sym.size(), jetAlg, Form("%s_E500",sm),
		  "JEC uncertainty (Data/MC)", s,
		  "E=500 GeV", 10,10,"fixE",500.);
  plotUncertainty(sym, 0, sym.size(), jetAlg, Form("%s_E2000",sm),
		  "JEC uncertainty (Data/MC)", s,
		  "E=2000 GeV", 10,10,"fixE",2000.);
  plotUncertainty(sym, 0, sym.size(), jetAlg, Form("%s_E1000",sm),
		  "JEC uncertainty (Data/MC)", s,
		  "E=1000 GeV", 10,10,"fixE",1000.); // moved last

  _canvas->SaveAs(Form("pdf/%s.pdf",sm));
  _icanvas = 1;

  string sspu = Form("%s_PileUp_%s",cu,names[jetAlg]);
  //		    d_algo==jec::AK5PF ? "AK5PF" :
  //		    (d_algo==jec::AK5PFchs ? "AK5PFchs" : "OTHER"));
  const char *spu = sspu.c_str();

  // PU uncertainty
  // vs pT
  plotUncertainty(sypu, 0, sypu.size(), jetAlg, Form("%s_Eta00",spu),
		  "JEC uncertainty (Data/MC)", s,
		  "|#eta_{jet}|=0", 10,10,"fixEta",0.);
  plotUncertainty(sypu, 0, sypu.size(), jetAlg, Form("%s_Eta20",spu),
		  "JEC uncertainty (Data/MC)", s,
		  "|#eta_{jet}|=2.0", 10,10,"fixEta",2.0);
  plotUncertainty(sypu, 0, sypu.size(), jetAlg, Form("%s_Eta27",spu),
		  "JEC uncertainty (Data/MC)", s,
		  "|#eta_{jet}|=2.7", 10,10,"fixEta",2.7);
  plotUncertainty(sypu, 0, sypu.size(), jetAlg, Form("%s_Eta42",spu),
		  "JEC uncertainty (Data/MC)", s,
		  "|#eta_{jet}|=4.2", 10,10,"fixEta",4.2);
  // vs eta
  plotUncertainty(sypu, 0, sypu.size(), jetAlg, Form("%s_Pt30",spu),
		  "JEC uncertainty (Data/MC)", s,
		  "p_{T}=30 GeV", 10,10,"fixPt",30.);
  plotUncertainty(sypu, 0, sypu.size(), jetAlg, Form("%s_Pt100",spu),
		  "JEC uncertainty (Data/MC)", s,
		  "p_{T}=100 GeV", 10,10,"fixPt",100.);
  plotUncertainty(sypu, 0, sypu.size(), jetAlg, Form("%s_E500",spu),
		  "JEC uncertainty (Data/MC)", s,
		  "E=500 GeV", 10,10,"fixE",500.);
  plotUncertainty(sypu, 0, sypu.size(), jetAlg, Form("%s_E2000",spu),
		  "JEC uncertainty (Data/MC)", s,
		  "E=2000 GeV", 10,10,"fixE",2000.);
  plotUncertainty(sypu, 0, sypu.size(), jetAlg, Form("%s_E1000",spu),
		  "JEC uncertainty (Data/MC)", s,
		  "E=1000 GeV", 10,10,"fixE",1000.); // moved last

  _canvas->SaveAs(Form("pdf/%s.pdf",spu));
  _icanvas = 1;

  string ssrel = Form("%s_Relative_%s",cu,names[jetAlg]);
  //		    d_algo==jec::AK5PF ? "AK5PF" :
  //		    (d_algo==jec::AK5PFchs ? "AK5PFchs" : "OTHER"));
  const char *srel = ssrel.c_str();

  // Relative scale uncertainty
  // vs pT
  plotUncertainty(syrel, 0, syrel.size(), jetAlg, Form("%s_Eta00",srel),
		  "JEC uncertainty (Data/MC)", s,
		  "|#eta_{jet}|=0", 10,10,"fixEta",0.);
  plotUncertainty(syrel, 0, syrel.size(), jetAlg, Form("%s_Eta20",srel),
		  "JEC uncertainty (Data/MC)", s,
		  "|#eta_{jet}|=2.0", 10,10,"fixEta",2.0);
  plotUncertainty(syrel, 0, syrel.size(), jetAlg, Form("%s_Eta27",srel),
		  "JEC uncertainty (Data/MC)", s,
		  "|#eta_{jet}|=2.7", 10,10,"fixEta",2.7);
  plotUncertainty(syrel, 0, syrel.size(), jetAlg, Form("%s_Eta42",srel),
		  "JEC uncertainty (Data/MC)", s,
		  "|#eta_{jet}|=4.2", 10,10,"fixEta",4.2);
  // vs eta
  plotUncertainty(syrel, 0, syrel.size(), jetAlg, Form("%s_Pt30",srel),
		  "JEC uncertainty (Data/MC)", s,
		  "p_{T}=30 GeV", 10,10,"fixPt",30.);
  plotUncertainty(syrel, 0, syrel.size(), jetAlg, Form("%s_Pt100",srel),
		  "JEC uncertainty (Data/MC)", s,
		  "p_{T}=100 GeV", 10,10,"fixPt",100.);
  plotUncertainty(syrel, 0, syrel.size(), jetAlg, Form("%s_E500",srel),
		  "JEC uncertainty (Data/MC)", s,
		  "E=500 GeV", 10,10,"fixE",500.);
  plotUncertainty(syrel, 0, syrel.size(), jetAlg, Form("%s_E2000",srel),
		  "JEC uncertainty (Data/MC)", s,
		  "E=2000 GeV", 10,10,"fixE",2000.);
  plotUncertainty(syrel, 0, syrel.size(), jetAlg, Form("%s_E1000",srel),
		  "JEC uncertainty (Data/MC)", s,
		  "E=1000 GeV", 10,10,"fixE",1000.); // moved last

  _canvas->SaveAs(Form("pdf/%s.pdf",srel));
  _icanvas = 1;


  string ssCorrGroups = Form("%s_CorrelationGroups_%s",cu,names[jetAlg]);
  //		    d_algo==jec::AK5PF ? "AK5PF" :
  //		    (d_algo==jec::AK5PFchs ? "AK5PFchs" : "OTHER"));
  const char *sCorrGroups = ssCorrGroups.c_str();

  // Correlation groups
  // vs pT
  plotUncertainty(syCorrGroups, 0, syCorrGroups.size(), jetAlg, Form("%s_Eta00",sCorrGroups),
		  "JEC uncertainty", "Anti-k_{T} R=0.5 PF", "|#eta_{jet}|=0",
		  10,10,"fixEta",0.);
  plotUncertainty(syCorrGroups, 0, syCorrGroups.size(), jetAlg, Form("%s_Eta20",srel),
		  "JEC uncertainty", "Anti-k_{T} R=0.5 PF", "|#eta_{jet}|=2.0",
		  10,10,"fixEta",2.0);
  plotUncertainty(syCorrGroups, 0, syCorrGroups.size(), jetAlg, Form("%s_Eta27",srel),
		  "JEC uncertainty", "Anti-k_{T} R=0.5 PF", "|#eta_{jet}|=2.7",
		  10,10,"fixEta",2.7);
  plotUncertainty(syCorrGroups, 0, syCorrGroups.size(), jetAlg, Form("%s_Eta42",srel),
		  "JEC uncertainty", "Anti-k_{T} R=0.5 PF", "|#eta_{jet}|=4.2",
		  10,10,"fixEta",4.2);
  // vs eta
  plotUncertainty(syCorrGroups, 0, syCorrGroups.size(), jetAlg, Form("%s_Pt30",srel),
		  "JEC uncertainty", "Anti-k_{T} R=0.5 PF", "p_{T}=30 GeV",
		  10,10,"fixPt",30.);
  plotUncertainty(syCorrGroups, 0, syCorrGroups.size(), jetAlg, Form("%s_Pt100",srel),
		  "JEC uncertainty", "Anti-k_{T} R=0.5 PF", "p_{T}=100 GeV",
		  10,10,"fixPt",100.);
  plotUncertainty(syCorrGroups, 0, syCorrGroups.size(), jetAlg, Form("%s_E500",srel),
		  "JEC uncertainty", "Anti-k_{T} R=0.5 PF", "E=500 GeV",
		  10,10,"fixPt",500.);
  plotUncertainty(syCorrGroups, 0, syCorrGroups.size(), jetAlg, Form("%s_E2000",srel),
		  "JEC uncertainty", "Anti-k_{T} R=0.5 PF", "E=2000 GeV",
		  10,10,"fixE",500.);
  plotUncertainty(syCorrGroups, 0, syCorrGroups.size(), jetAlg, Form("%s_E1000",srel),
		  "JEC uncertainty", "Anti-k_{T} R=0.5 PF", "E=1000 GeV",
		  10,10,"fixE",2000.);

  _canvas->SaveAs(Form("pdf/%s.pdf",sCorrGroups));
  _icanvas = 1;



  string sspt = Form("%s_AbsolutePt_%s",cu,names[jetAlg]);
  //		    d_algo==jec::AK5PF ? "AK5PF" :
  //		    (d_algo==jec::AK5PFchs ? "AK5PFchs" : "OTHER"));
  const char *spt = sspt.c_str();

  // Absolute scale pT dependent uncertainty
  // vs pT
  plotUncertainty(sypt, 0, sypt.size(), jetAlg, Form("%s_Eta00",spt),
		  "JEC uncertainty (Data/MC)", s,
		  "|#eta_{jet}|=0", 10,10,"fixEta",0.);
  plotUncertainty(sypt, 0, sypt.size(), jetAlg, Form("%s_Eta20",spt),
		  "JEC uncertainty (Data/MC)", s,
		  "|#eta_{jet}|=2.0", 10,10,"fixEta",2.0);
  plotUncertainty(sypt, 0, sypt.size(), jetAlg, Form("%s_Eta27",spt),
		  "JEC uncertainty (Data/MC)", s,
		  "|#eta_{jet}|=2.7", 10,10,"fixEta",2.7);
  plotUncertainty(sypt, 0, sypt.size(), jetAlg, Form("%s_Eta42",spt),
		  "JEC uncertainty (Data/MC)", s,
		  "|#eta_{jet}|=4.2", 10,10,"fixEta",4.2);
  // vs eta
  plotUncertainty(sypt, 0, sypt.size(), jetAlg, Form("%s_Pt30",spt),
		  "JEC uncertainty (Data/MC)", s,
		  "p_{T}=30 GeV", 10,10,"fixPt",30.);
  plotUncertainty(sypt, 0, sypt.size(), jetAlg, Form("%s_Pt100",spt),
		  "JEC uncertainty (Data/MC)", s,
		  "p_{T}=100 GeV", 10,10,"fixPt",100.);
  plotUncertainty(sypt, 0, sypt.size(), jetAlg, Form("%s_E500",spt),
		  "JEC uncertainty (Data/MC)", s,
		  "E=500 GeV", 10,10,"fixE",500.);
  plotUncertainty(sypt, 0, sypt.size(), jetAlg, Form("%s_E2000",spt),
		  "JEC uncertainty (Data/MC)", s,
		  "E=2000 GeV", 10,10,"fixE",2000.);
  plotUncertainty(sypt, 0, sypt.size(), jetAlg, Form("%s_E1000",spt),
		  "JEC uncertainty (Data/MC)", s,
		  "E=1000 GeV", 10,10,"fixE",1000.); // moved last

  _canvas->SaveAs(Form("pdf/%s.pdf",spt));
  _icanvas = 1;

  string ssf = Form("%s_Flavor_%s",cu,names[jetAlg]);
  const char *sf = ssf.c_str();

  minimaltmp = _minimal;
  if (algo=="AK5PF") _minimal = false;

  // Flavor uncertainty
  // vs pT
  plotUncertainty(syf, 0, syf.size(), jetAlg, Form("%s_Eta00",sf),
                  "JEC uncertainty (Data/MC)", s,
                  "|#eta_{jet}|=0", 10,10,"fixEta",0.);
  plotUncertainty(syf, 0, syf.size(), jetAlg, Form("%s_Eta20",sf),
                  "JEC uncertainty (Data/MC)", s,
                  "|#eta_{jet}|=2.0", 10,10,"fixEta",2.0);
  plotUncertainty(syf, 0, syf.size(), jetAlg, Form("%s_Eta27",sf),
                  "JEC uncertainty (Data/MC)", s,
                  "|#eta_{jet}|=2.7", 10,10,"fixEta",2.7);
  plotUncertainty(syf, 0, syf.size(), jetAlg, Form("%s_Eta42",sf),
                  "JEC uncertainty (Data/MC)", s,
                  "|#eta_{jet}|=4.2", 10,10,"fixEta",4.2);
  // vs eta
  plotUncertainty(syf, 0, syf.size(), jetAlg, Form("%s_Pt30",sf),
                  "JEC uncertainty (Data/MC)", s,
                  "p_{T}=30 GeV", 10,10,"fixPt",30.);
  plotUncertainty(syf, 0, syf.size(), jetAlg, Form("%s_Pt100",sf),
                  "JEC uncertainty (Data/MC)", s,
                  "p_{T}=100 GeV", 10,10,"fixPt",100.);
  plotUncertainty(syf, 0, syf.size(), jetAlg, Form("%s_E500",sf),
                  "JEC uncertainty (Data/MC)", s,
                  "E=500 GeV", 10,10,"fixE",500.);
  plotUncertainty(syf, 0, syf.size(), jetAlg, Form("%s_E2000",sf),
                  "JEC uncertainty (Data/MC)", s,
                  "E=2000 GeV", 10,10,"fixE",2000.);
  plotUncertainty(syf, 0, syf.size(), jetAlg, Form("%s_E1000",sf),
                  "JEC uncertainty (Data/MC)", s,
                  "E=1000 GeV", 10,10,"fixE",1000.); // moved last

  _canvas->SaveAs(Form("pdf/%s.pdf",sf));
  _icanvas = 1;

  _minimal = minimaltmp;

  bool absUncertTmp = _absUncert;
  _absUncert = false;
  string sss = Form("%s_SinglePion_%s",cu,names[jetAlg]);
  const char *css = sss.c_str();

  // Single pion versus eta dependence
  // vs pT
  plotUncertainty(sys, 0, sys.size(), jetAlg, Form("%s_Eta00",css),
                  "JEC uncertainty (Data/MC)", s,
                  "|#eta_{jet}|=0", -5,10,"fixEta",0.);
  plotUncertainty(sys, 0, sys.size(), jetAlg, Form("%s_Eta00",css),
                  "JEC uncertainty (Data/MC)", s,
                  "|#eta_{jet}|=0.25", -5,10,"fixEta",0.25);
  plotUncertainty(sys, 0, sys.size(), jetAlg, Form("%s_Eta00",css),
                  "JEC uncertainty (Data/MC)", s,
                  "|#eta_{jet}|=0.75", -5,10,"fixEta",0.75);
  plotUncertainty(sys, 0, sys.size(), jetAlg, Form("%s_Eta00",css),
                  "JEC uncertainty (Data/MC)", s,
                  "|#eta_{jet}|=1.25", -5,10,"fixEta",1.25);
  plotUncertainty(sys, 0, sys.size(), jetAlg, Form("%s_Eta00",css),
                  "JEC uncertainty (Data/MC)", s,
                  "|#eta_{jet}|=1.75", -5,10,"fixEta",1.75);
  plotUncertainty(sys, 0, sys.size(), jetAlg, Form("%s_Eta20",css),
                  "JEC uncertainty (Data/MC)", s,
                  "|#eta_{jet}|=2.25", -5,10,"fixEta",2.25);
  //plotUncertainty(sys, 0, sys.size(), jetAlg, Form("%s_Eta27",css),
  //              "JEC uncertainty (Data/MC)", s,
  //              "|#eta_{jet}|=2.7", -5,10,"fixEta",2.7);
  //plotUncertainty(sys, 0, sys.size(), jetAlg, Form("%s_Eta42",css),
  //              "JEC uncertainty (Data/MC)", s,
  //              "|#eta_{jet}|=4.2", -5,10,"fixEta",4.2);
  // vs eta
  //plotUncertainty(sys, 0, sys.size(), jetAlg, Form("%s_Pt30",css),
  //            "JEC uncertainty (Data/MC)", s,
  //            "p_{T}=30 GeV", -5,10,"fixPt",30.);
  plotUncertainty(sys, 0, sys.size(), jetAlg, Form("%s_Pt100",css),
                  "JEC uncertainty (Data/MC)", s,
                  "p_{T}=100 GeV", -5,10,"fixPt",100.);
  //plotUncertainty(sys, 0, sys.size(), jetAlg, Form("%s_E500",css),
  //              "JEC uncertainty (Data/MC)", s,
  //              "E=500 GeV", -5,10,"fixE",500.);
  plotUncertainty(sys, 0, sys.size(), jetAlg, Form("%s_E2000",css),
                "JEC uncertainty (Data/MC)", s,
                "E=2000 GeV", -5,10,"fixE",2000.);
  //plotUncertainty(sys, 0, sys.size(), jetAlg, Form("%s_E1000",css),
  //              "JEC uncertainty (Data/MC)", s,
  //              "E=1000 GeV", -5,10,"fixE",1000.); // moved last
  _absUncert = absUncertTmp;

  _canvas->SaveAs(Form("pdf/%s.pdf",css));
  _icanvas = 0;
  delete _canvas;

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
     2000, 2238, 2500, 2787, 3103, 3450};
  const int ndiv_pt = sizeof(x_pt)/sizeof(x_pt[0])-1;
  const double x_eta[] =
    /*
    {-5.2,-5.199,-4.4,-4,-3.5,-3,-2.8,-2.6,-2.4,-2.2,-2.0,
     -1.8,-1.6,-1.4,-1.2,-1.0, -0.8,-0.6,-0.4,-0.2,0.,
     0.2,0.4,0.6,0.8,1.0,1.2,1.4,
     1.6,1.8,2.0,2.2,2.4,2.6,2.8,3,3.5,4,4.4,5.199,5.2};
    */
  // Up to Summer13
    {-5.4,-5.0,-4.4,-4,-3.5,-3,-2.8,-2.6,-2.4,-2.2,-2.0,
     -1.8,-1.6,-1.4,-1.2,-1.0, -0.8,-0.6,-0.4,-0.2,0.,
     0.2,0.4,0.6,0.8,1.0,1.2,1.4,
     1.6,1.8,2.0,2.2,2.4,2.6,2.8,3,3.5,4,4.4,5.0,5.4};
  /*
    {-5.4,-5.2,-5.0,-4.8,-4.6,-4.4,-4.2,-4,-3.8,
     -3.6,-3.4,-3.2, -3,-2.8,-2.6,-2.4,-2.2,-2.0,
     -1.8,-1.6,-1.4,-1.2,-1.0, -0.8,-0.6,-0.4,-0.2,0.,
     0.2,0.4,0.6,0.8,1.0,1.2,1.4, 1.6,1.8,
     2.0,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,
     3.8,4,4.2,4.4,4.6,4.8,5.0,5.2,5.4};
  */
  const int ndiv_eta = sizeof(x_eta)/sizeof(x_eta[0])-1;     

  // Re-determine bin edges based on maxe=4000.
  int jx(0), ndiv_new(0);
  if (type=="fixEta") {
    double maxe = 4000;
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
    double maxe = 4000;
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
  int ndiv0;
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

  TCanvas *c1 = new TCanvas(Form("c1_%s",name.c_str()),"c1",600,600);
  if (type=="fixEta") c1->SetLogx();
  if (c2) { if (type=="fixEta") c2->SetLogx(); }

  const char *cy = ytitle.c_str();
  TH1D *h0 = new TH1D(Form("h0_%s_%d",name.c_str(),++icnt),
		      Form(";p_{T} (GeV);%s [%%]",cy),
		      ndiv0, &x0[0]);//ndiv-2,&x[1]);
  if (type=="fixEta") h0->GetXaxis()->SetTitle("p_{T} (GeV)");
  if (type=="fixPt") h0->GetXaxis()->SetTitle("#eta_{jet}");
  if (type=="fixE") h0->GetXaxis()->SetTitle("#eta_{jet}");
  //h0->SetMinimum(emax<0 ? emax : 0.);
  h0->SetMinimum(_absUncert ? 0 : -0.5*fabs(emax));
  h0->SetMaximum(fabs(emax));
  h0->GetXaxis()->SetMoreLogLabels();
  h0->GetXaxis()->SetNoExponent();
  h0->GetYaxis()->SetTitleOffset(1.0);
  h0->Draw("AXIS");
  if (c2) { c2->cd(); h0->DrawClone("AXIS"); c1->cd(); }

  //cout << "Got here 3" << endl << flush;

  TLegend *leg1 = new TLegend(0.18,0.93-0.05*nsys1,0.38,0.93,"","brNDC");
  leg1->SetFillStyle(kNone);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.045);
  const double tx = (type=="fixEta" ? 0.54 : 0.40);
  TLegend *leg2 = new TLegend(tx,0.93-0.05*nsys2,tx+0.20,0.93,"","brNDC");
  leg2->SetFillStyle(kNone);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.045);
  TLatex *tex1 = new TLatex(tx, 0.85-0.05*max(nsys1,nsys2),label1.c_str());
  tex1->SetNDC();
  tex1->SetTextSize(0.045);
  TLatex *tex2 = new TLatex(tx, 0.80-0.05*max(nsys1,nsys2),label2.c_str());
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
    //bool mpf = (u.method=="default" ? d_mpf : u.method=="mpf");
    double npv = (u.npv==-1 ? d_npv : u.npv);
    //L3Corr rjet(jetAlg, u.type, d_id, npu, mpf);
    JECUncertainty rjet(jetAlg, jec::DATA, u.type, npv);
    TGraph *g = new TGraph(0);//ndiv);
    g->SetName(Form("L3_%s",u.name.c_str()));
    TGraph *gup = new TGraph(0);//ndiv);
    gup->SetName(Form("L3_%s_up",u.name.c_str()));
    TGraph *gdw = new TGraph(0);//ndiv);
    gdw->SetName(Form("L3_%s_dw",u.name.c_str()));
    TH1D *h = new TH1D(Form("L3_%s_%s_%d",name.c_str(),u.name.c_str(),++icnt),
		       "",ndiv,&x[0]);

    for (int ix = 0; ix != ndiv; ++ix) {
      //cout << "." << flush;
      double var = 0.5*(x[ix]+x[ix+1]);
      double pt(0), eta(0);
      if (type=="fixEta") { pt = var; eta = typevar; }
      if (type=="fixPt") { eta = var; pt = typevar; }
      if (type=="fixE") { eta = var; pt = typevar/cosh(eta); }
      //if (pt*cosh(eta) < 3500. && pt > ptmin) {
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
      }
    } // for ix

    TString opt2 = u.opt;
    if (opt2.Contains("F")) {
      //cout << "x" << flush;
      h->SetMarkerStyle(kNone);
      h->SetLineStyle(kNone);
      h->SetFillStyle(u.fstyle);
      h->SetFillColor(u.fcolor);
      h->Draw("SAME E3");
      if (c2) { c2->cd(); h->DrawClone("SAME E3"); c1->cd(); }
      opt2.ReplaceAll("F","");
    }
    g->SetLineColor(u.lcolor);
    g->SetLineStyle(u.lstyle);
    g->SetLineWidth(u.lwid);
    g->SetMarkerColor(u.mcolor);
    g->SetMarkerStyle(u.mstyle);
    g->SetFillColor(u.fcolor);
    g->SetFillStyle(u.fstyle);
    g->Draw(opt2);
    if (c2) { c2->cd(); g->DrawClone(opt2); c1->cd(); }

    if (isys<nsys1) leg1->AddEntry(g, u.title.c_str(), u.opt);
    else            leg2->AddEntry(g, u.title.c_str(), u.opt);


    // Print uncertainty into a text file for Kostas
    // Remember the following:
    // 1) no signs in front of the numbers
    // 2) we need the original numbers, not the ones used for presentation purposes (e.g. 0.04 NOT 4%).
    // 3) the grid has the form (pt, uncertainty, uncertainty) 
    /*
    if (u.name=="tot") {
      ofstream fout(Form("txt/%s.txt",name.c_str()),ios::out);
      fout << "{1 JetEta 1 JetPt \"\" Correction L2Relative}" << endl;
      for (int i = 0; i != g->GetN(); ++i) {
	//fout << Form("%1.1f +%1.3g -%1.3g ", g->GetX()[i],
	//	       g->GetY()[i], g->GetY()[i]);
	fout << Form("%1.1f %1.3g %1.3g ", g->GetX()[i],
		     0.01*g->GetY()[i], 0.01*g->GetY()[i]);
      }
      fout << endl;
    }
    */
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
  if (TString(name.c_str()).Contains("JECUncert_Flavor"))
    cmsPrel(0);
  else
  if (name=="JECUncert_Offset_PFAK5" ||
      name=="JECUncert_Offset_CALOAK5" ||
      name=="JECUncert_MPF" ||
      name=="JECUncert_HighPt_PFAK5" ||
      name=="JECUncert_HighPt_JPTAK5" ||
      name=="JECUncert_HighPt_CALOAK5" ||
      name=="JECUncert_PFAK5_summary" ||
      name=="JECUncert_JPTAK5_summary" ||
      name=="JECUncert_AK5_summary")
    cmsFinal(_lumi);
  else
    cmsPrel(_lumi);//2.9);//1.2);//60);
  gPad->RedrawAxis();
  if (c2) { c2->cd(); cmsPrel(_lumi); gPad->RedrawAxis(); c1->cd(); }
  //h0->Draw("SAME AXIS");

  if (!_minimal) {
    if (_eps) c1->SaveAs(("eps/"+name+".eps").c_str());
    if (_pdf) c1->SaveAs(("pdf/"+name+".pdf").c_str());
  }    

  //cout << "Got here 6" << endl << flush;

  if (name=="JECUncert_DATA_Summary_AK5PF_Eta00") {
    //if (name=="JECUncert_DATA_AK5PFchs_Eta00") {

    JECUncertainty rjet5p(jec::AK5PF, jec::DATA, jec::kData, d_npv);
    JECUncertainty rjet7p(jec::AK7PF, jec::DATA, jec::kData, d_npv);
    JECUncertainty rjet5s(jec::AK5PFchs, jec::DATA, jec::kData, d_npv);
    JECUncertainty rjet7s(jec::AK7PFchs, jec::DATA, jec::kData, d_npv);
    JECUncertainty rjet5c(jec::AK5CALO, jec::DATA, jec::kData, d_npv);
    JECUncertainty rjet7c(jec::AK5CALO, jec::DATA, jec::kData, d_npv);

    ofstream fout5p("txt/Summer13_V3_DATA_Uncertainty_AK5PF.txt",ios::out);
    fout5p << "{1 JetEta 1 JetPt \"\" Correction Uncertainty}" << endl;
    ofstream fout5s("txt/Summer13_V3_DATA_Uncertainty_AK5PFchs.txt",ios::out);
    fout5s << "{1 JetEta 1 JetPt \"\" Correction Uncertainty}" << endl;
    ofstream fout5c("txt/Summer13_V3_DATA_Uncertainty_AK5Calo.txt",ios::out);
    fout5c << "{1 JetEta 1 JetPt \"\" Correction Uncertainty}" << endl;
    ofstream fout7p("txt/Summer13_V3_DATA_Uncertainty_AK7PF.txt",ios::out);
    fout7p << "{1 JetEta 1 JetPt \"\" Correction Uncertainty}" << endl;
    ofstream fout7s("txt/Summer13_V3_DATA_Uncertainty_AK7PFchs.txt",ios::out);
    fout7s << "{1 JetEta 1 JetPt \"\" Correction Uncertainty}" << endl;
    ofstream fout7c("txt/Summer13_V3_DATA_Uncertainty_AK7Calo.txt",ios::out);
    fout7c << "{1 JetEta 1 JetPt \"\" Correction Uncertainty}" << endl;

    JECUncertainty rjet5px(jec::AK5PF, jec::DATA, jec::kMC, d_npv);
    JECUncertainty rjet7px(jec::AK7PF, jec::DATA, jec::kMC, d_npv);
    JECUncertainty rjet5sx(jec::AK5PFchs, jec::DATA, jec::kMC, d_npv);
    JECUncertainty rjet7sx(jec::AK7PFchs, jec::DATA, jec::kMC, d_npv);
    JECUncertainty rjet5cx(jec::AK5CALO, jec::DATA, jec::kMC, d_npv);
    JECUncertainty rjet7cx(jec::AK7CALO, jec::DATA, jec::kMC, d_npv);

    ofstream fout5px("txt/Summer13_V3_MC_Uncertainty_AK5PF.txt",ios::out);
    fout5px << "{1 JetEta 1 JetPt \"\" Correction Uncertainty}" << endl;
    ofstream fout5sx("txt/Summer13_V3_MC_Uncertainty_AK5PFchs.txt",ios::out);
    fout5sx << "{1 JetEta 1 JetPt \"\" Correction Uncertainty}" << endl;
    ofstream fout5cx("txt/Summer13_V3_MC_Uncertainty_AK5Calo.txt",ios::out);
    fout5cx << "{1 JetEta 1 JetPt \"\" Correction Uncertainty}" << endl;
    ofstream fout7px("txt/Summer13_V3_MC_Uncertainty_AK7PF.txt",ios::out);
    fout7px << "{1 JetEta 1 JetPt \"\" Correction Uncertainty}" << endl;
    ofstream fout7sx("txt/Summer13_V3_MC_Uncertainty_AK7PFchs.txt",ios::out);
    fout7sx << "{1 JetEta 1 JetPt \"\" Correction Uncertainty}" << endl;
    ofstream fout7cx("txt/Summer13_V3_MC_Uncertainty_AK7Calo.txt",ios::out);
    fout7cx << "{1 JetEta 1 JetPt \"\" Correction Uncertainty}" << endl;

    for (int ieta = 0; ieta != ndiv_eta; ++ieta) {
      
      double etamin = x_eta[ieta];
      double etamax = x_eta[ieta+1];
      double eta = 0.5*(etamin+etamax);
      
      fout5p << Form("%1.1f %1.1f %d ",etamin,etamax,ndiv_pt*3);
      fout5s << Form("%1.1f %1.1f %d ",etamin,etamax,ndiv_pt*3);
      fout5c << Form("%1.1f %1.1f %d ",etamin,etamax,ndiv_pt*3);
      fout7p << Form("%1.1f %1.1f %d ",etamin,etamax,ndiv_pt*3);
      fout7s << Form("%1.1f %1.1f %d ",etamin,etamax,ndiv_pt*3);
      fout7c << Form("%1.1f %1.1f %d ",etamin,etamax,ndiv_pt*3);
      //
      fout5px << Form("%1.1f %1.1f %d ",etamin,etamax,ndiv_pt*3);
      fout5sx << Form("%1.1f %1.1f %d ",etamin,etamax,ndiv_pt*3);
      fout5cx << Form("%1.1f %1.1f %d ",etamin,etamax,ndiv_pt*3);
      fout7px << Form("%1.1f %1.1f %d ",etamin,etamax,ndiv_pt*3);
      fout7sx << Form("%1.1f %1.1f %d ",etamin,etamax,ndiv_pt*3);
      fout7cx << Form("%1.1f %1.1f %d ",etamin,etamax,ndiv_pt*3);
      
      for (int ipt = 0; ipt != ndiv_pt; ++ipt) {
	
	double pt = 0.5*(x_pt[ipt]+x_pt[ipt+1]);

	{
	  double err5(0);
	  double r5 = 1;//rjet5.Rjet(pt, eta, err5);
	  err5 = rjet5p.Uncert(pt, eta);
	  err5 /= r5; // relative uncertainty
	  fout5p << Form("%1.1f %1.4f %1.4f ", pt, err5, err5);
	  //double err5c = sqrt(err5*err5 + 0.015*0.015);
	  //fout5c << Form("%1.1f %1.4f %1.4f ", pt, err5c, err5c);
	  double err5s = rjet5s.Uncert(pt, eta);
	  fout5s << Form("%1.1f %1.4f %1.4f ", pt, err5s, err5s);
	  double err5c = rjet5c.Uncert(pt, eta);
	  fout5c << Form("%1.1f %1.4f %1.4f ", pt, err5c, err5c);
	  
	  double err7(0);
	  double r7 = 1.;//rjet7.Rjet(pt, eta, err7);
	  err7 = rjet7p.Uncert(pt, eta);
	  err7 /= r7; // relative uncertainty
	  fout7p << Form("%1.1f %1.4f %1.4f ", pt, err7, err7);
	  //double err7c = sqrt(err7*err7 + 0.015*0.015);
	  //fout7c << Form("%1.1f %1.4f %1.4f ", pt, err7c, err7c);
	  double err7s = rjet7s.Uncert(pt, eta);
	  fout7s << Form("%1.1f %1.4f %1.4f ", pt, err7s, err7s);
	  double err7c = rjet7c.Uncert(pt, eta);
	  fout7c << Form("%1.1f %1.4f %1.4f ", pt, err7c, err7c);
	}
	{
	  double err5x(0);
	  double r5x = 1.;//rjet5x.Rjet(pt, eta, err5x);
	  err5x = rjet5px.Uncert(pt, eta);
	  err5x /= r5x; // relative uncertainty
	  fout5px << Form("%1.1f %1.4f %1.4f ", pt, err5x, err5x);
	  //double err5cx = sqrt(err5x*err5x + 0.015*0.015);
	  //fout5cx << Form("%1.1f %1.4f %1.4f ", pt, err5cx, err5cx);	  
	  double err5sx = rjet5sx.Uncert(pt, eta);
	  fout5sx << Form("%1.1f %1.4f %1.4f ", pt, err5sx, err5sx);
	  double err5cx = rjet5cx.Uncert(pt, eta);
	  fout5cx << Form("%1.1f %1.4f %1.4f ", pt, err5cx, err5cx);
	  
	  double err7x(0);
	  double r7x = 1.;//rjet7x.Rjet(pt, eta, err7x);
	  err7x = rjet7px.Uncert(pt, eta);
	  err7x /= r7x; // relative uncertainty
	  fout7px << Form("%1.1f %1.4f %1.4f ", pt, err7x, err7x);
	  //double err7cx = sqrt(err7x*err7x + 0.015*0.015);
	  //fout7cx << Form("%1.1f %1.4f %1.4f ", pt, err7cx, err7cx);
	  double err7sx = rjet7sx.Uncert(pt, eta);
	  fout7sx << Form("%1.1f %1.4f %1.4f ", pt, err7sx, err7sx);
	  double err7cx = rjet7cx.Uncert(pt, eta);
	  fout7cx << Form("%1.1f %1.4f %1.4f ", pt, err7cx, err7cx);
	}
      } // for ipt
      fout5p << endl;
      fout5s << endl;
      fout5c << endl;
      fout7p << endl;
      fout7s << endl;
      fout7c << endl;
      //
      fout5px << endl;
      fout5sx << endl;
      fout5cx << endl;
      fout7px << endl;
      fout7sx << endl;
      fout7cx << endl;
    } // for ieta
  } // print uncertainty

  if (name=="JECUncert_DATA_Summary_AK5PF_Eta00") {
    //if (name=="JECUncert_DATA_AK5PFchs_Eta00") {
    
    // Note: AK5PFchs is CHS, AK7PF is non-CHS (AK7PFchs on Jan 25)
    ofstream fout5("txt/Summer13_V3_DATA_UncertaintySources_AK5PFchs.txt",ios::out);
    fout5 << "#Uncertainty sources for Summer13_V3_DATA_AK5PFchs" << endl;
    cout << "Storing uncertainties to: "
	 << "txt/Summer13_V3_DATA_UncertaintySources_AK5PFchs.txt" << endl;
    ofstream fout5x("txt/Summer13_V3_DATA_UncertaintySources_AK5PF.txt",ios::out);
    fout5x << "#Uncertainty sources for Summer13_V3_DATA_AK5PF" << endl;
    cout << "Storing uncertainties to: "
	 << "txt/Summer13_V3_DATA_UncertaintySources_AK5PF.txt" << endl;
    //
    ofstream fout7("txt/Summer13_V3_DATA_UncertaintySources_AK7PFchs.txt",
		   ios::out);
    fout7 << "#Uncertainty sources for Summer13_V3_DATA_AK7PFchs" << endl;
    cout << "Storing uncertainties to: "
	 << "txt/Summer13_V3_DATA_UncertaintySources_AK7PFchs.txt" << endl;
    ofstream fout7x("txt/Summer13_V3_DATA_UncertaintySources_AK7PF.txt",
		    ios::out);
    fout7x << "#Uncertainty sources for Summer13_V3_DATA_AK7PF" << endl;
    cout << "Storing uncertainties to: "
	 << "txt/Summer13_V3_DATA_UncertaintySources_AK7PF.txt" << endl;

    jec::ErrorTypes vsrc[] =
      //{jec::kAbsolute, jec::kRelative, jec::kPtExtra};
      {jec::kAbsoluteStat, jec::kAbsoluteScale, jec::kAbsoluteFlavorMapping, jec::kAbsoluteMPFBias,
       jec::kAbsoluteFrag, /*jec::kAbsoluteSPR,*/
       jec::kAbsoluteSPRE, jec::kAbsoluteSPRH,
       /*jec::kAbsoluteECAL, jec::kAbsoluteTrack,*/
       /*jec::kFlavorMC,*/ jec::kFlavorQCD,/*new*/ jec::kTime,
       jec::kRelativeJEREC1, jec::kRelativeJEREC2, jec::kRelativeJERHF,
       jec::kRelativePtBB, /*new*/
       jec::kRelativePtEC1, jec::kRelativePtEC2, jec::kRelativePtHF,
       jec::kRelativeFSR,/*new*/ jec::kRelativeStatEC2, jec::kRelativeStatHF,
       /*jec::kRelativeSample,*/
       jec::kPileUpDataMC, /*jec::kPileUpOOT,*/ /*jec::kPileUpPt,*/
       jec::kPileUpPtBB, jec::kPileUpPtEC, jec::kPileUpPtHF,
       jec::kPileUpBias,
       /*jec::kPileUpJetRate,*/
       jec::kPileUp, jec::kRelative, jec::kAbsolutePt,//jec::kPtExtra,
       jec::kMC,
       jec::kData,//Total};
       jec::kDataNoFlavor, jec::kFlavorZJet, jec::kFlavorPhotonJet,
       jec::kFlavorPureGluon, jec::kFlavorPureQuark, 
       jec::kFlavorPureCharm, jec::kFlavorPureBottom,
       jec::kCorrelationGroupMPFInSitu, jec::kCorrelationGroupIntercalibration, jec::kCorrelationGroupbJES,
       jec::kCorrelationGroupFlavor, jec::kCorrelationGroupUncorrelated
      };

    const int nsrc = sizeof(vsrc)/sizeof(vsrc[0]);
    map<jec::ErrorTypes, string> srcname;
    //srcname[jec::kPtExtra] = "PtExtra";
    //srcname[jec::kAbsoluteScale] = "Absolute";
    srcname[jec::kAbsoluteStat] = "AbsoluteStat";
    srcname[jec::kAbsoluteScale] = "AbsoluteScale";
    srcname[jec::kAbsoluteFlavorMapping] = "AbsoluteFlavMap";
    srcname[jec::kAbsoluteMPFBias] = "AbsoluteMPFBias";
    //srcname[jec::kRelative] = "Relative";
    srcname[jec::kRelativeJEREC1] = "RelativeJEREC1";
    srcname[jec::kRelativeJEREC2] = "RelativeJEREC2";
    srcname[jec::kRelativeJERHF] = "RelativeJERHF";
    srcname[jec::kRelativePtBB] = "RelativePtBB"; // new in Summer13_V1
    srcname[jec::kRelativePtEC1] = "RelativePtEC1";
    srcname[jec::kRelativePtEC2] = "RelativePtEC2";
    srcname[jec::kRelativePtHF] = "RelativePtHF";
    srcname[jec::kRelativeFSR] = "RelativeFSR"; // new in Summer13_V1
    srcname[jec::kRelativeStatEC2] = "RelativeStatEC2";
    srcname[jec::kRelativeStatHF] = "RelativeStatHF";
    //srcname[jec::kRelativeSample] = "RelativeSample";
    srcname[jec::kAbsoluteFrag] = "HighPtExtra"; // update uncertainty?
    //srcname[jec::kAbsoluteSPR] = "SinglePion"; // shape ok?
    srcname[jec::kAbsoluteSPRE] = "SinglePionECAL";
    srcname[jec::kAbsoluteSPRH] = "SinglePionHCAL";
    //srcname[jec::kAbsoluteECAL] = "SinglePion";
    //srcname[jec::kAbsoluteTrack] = "SinglePion";
    //srcname[jec::kPileUp] = "PileUp";
    srcname[jec::kPileUpDataMC] = "PileUpDataMC";
    //srcname[jec::kPileUpOOT] = "PileUpOOT";
    //srcname[jec::kPileUpPt] = "PileUpPt";
    srcname[jec::kPileUpPtBB] = "PileUpPtBB";
    srcname[jec::kPileUpPtEC] = "PileUpPtEC";
    srcname[jec::kPileUpPtHF] = "PileUpPtHF";
    srcname[jec::kPileUpBias] = "PileUpBias";
    //srcname[jec::kPileUpJetRate] = "PileUpJetRate";
    //srcname[jec::kFlavorMC] = "Flavor";
    srcname[jec::kFlavorQCD] = "FlavorQCD";
    srcname[jec::kTime] = "Time";
    srcname[jec::kPileUp] = "SubTotalPileUp";
    srcname[jec::kRelative] = "SubTotalRelative";
    //srcname[jec::kPtExtra] = "SubTotalPt";
    srcname[jec::kAbsolutePt] = "SubTotalPt";
    srcname[jec::kMC] = "SubTotalMC";
    srcname[jec::kData] = "Total";
    srcname[jec::kDataNoFlavor] = "TotalNoFlavor";
    srcname[jec::kFlavorZJet] = "FlavorZJet";
    srcname[jec::kFlavorPhotonJet] = "FlavorPhotonJet";
    srcname[jec::kFlavorPureGluon] = "FlavorPureGluon";
    srcname[jec::kFlavorPureQuark] = "FlavorPureQuark";
    srcname[jec::kFlavorPureBottom] = "FlavorPureBottom";
    srcname[jec::kFlavorPureCharm] = "FlavorPureCharm";
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
      JECUncertainty rjet5(jec::AK5PFchs, jec::DATA, src, d_npv);
      fout5 << "["<<srcname[src]<<"]" << endl;
      fout5 << "{1 JetEta 1 JetPt \"\" Correction JECSource}" << endl;
      JECUncertainty rjet5x(jec::AK5PF, jec::DATA, src, d_npv);
      fout5x << "["<<srcname[src]<<"]" << endl;
      fout5x << "{1 JetEta 1 JetPt \"\" Correction JECSource}" << endl;
      JECUncertainty rjet7(jec::AK7PFchs, jec::DATA, src, d_npv);
      fout7 << "["<<srcname[src]<<"]" << endl;
      fout7 << "{1 JetEta 1 JetPt \"\" Correction JECSource}" << endl;
      JECUncertainty rjet7x(jec::AK7PF, jec::DATA, src, d_npv);
      fout7x << "["<<srcname[src]<<"]" << endl;
      fout7x << "{1 JetEta 1 JetPt \"\" Correction JECSource}" << endl;

      for (int ieta = 0; ieta != ndiv_eta; ++ieta) {

	double etamin = x_eta[ieta];
	double etamax = x_eta[ieta+1];
	double eta = 0.5*(etamin+etamax);
	fout5 << Form("%1.1f %1.1f %d ",etamin,etamax,ndiv_pt*3);
	fout5x << Form("%1.1f %1.1f %d ",etamin,etamax,ndiv_pt*3);
	fout7 << Form("%1.1f %1.1f %d ",etamin,etamax,ndiv_pt*3);
	fout7x << Form("%1.1f %1.1f %d ",etamin,etamax,ndiv_pt*3);

	for (int ipt = 0; ipt != ndiv_pt; ++ipt) {
	
	  double pt = 0.5*(x_pt[ipt]+x_pt[ipt+1]);
	  double err5(0);
	  double r5 = 1.;//rjet5.Rjet(pt, eta, err5);
	  err5 = rjet5.Uncert(pt, eta);
	  err5 /= r5; // relative uncertainty
	  fout5 << Form("%1.1f %1.4f %1.4f ", pt, err5, err5);
	  //
	  double err5x(0);
	  double r5x = 1.;//rjet5x.Rjet(pt, eta, err5x);
	  err5x = rjet5x.Uncert(pt, eta);
	  err5x /= r5x; // relative uncertainty
	  fout5x << Form("%1.1f %1.4f %1.4f ", pt, err5x, err5x);
	  //
	  double err7(0);
	  double r7 = 1.;//rjet7.Rjet(pt, eta, err7);
	  err7 = rjet7.Uncert(pt, eta);
	  err7 /= r7; // relative uncertainty
	  fout7 << Form("%1.1f %1.4f %1.4f ", pt, err7, err7);
	  //
	  double err7x(0);
	  double r7x = 1.;//rjet7x.Rjet(pt, eta, err7x);
	  err7x = rjet7x.Uncert(pt, eta);
	  err7x /= r7x; // relative uncertainty
	  fout7x << Form("%1.1f %1.4f %1.4f ", pt, err7x, err7x);
	} // for ipt
	fout5 << endl;
	fout5x << endl;
	fout7 << endl;
	fout7x << endl;
      } // for ieta
    } // for isrc
    
    
  } // print uncertainty sources

} // plotUncertainty

