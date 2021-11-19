// Purpose: Fit response and composition variation shapes from FullSim.
//          This is evolution of minitools/varPlots.C, dropping old baggage
//          Comparisons are done to old shapes from toyPF
//          (use ToyPF as placeholder for missing plots)
#include "TFile.h"
#include "TF1.h"

#include "../tdrstyle_mod15.C"

#include <vector>
#include <string>

// Settings
bool patchECALm3 = false;//true;
bool plotToyPF = false;//true;
string CorLevel = "L1L2L3";
#include "../globalFitL3Res.C" // toyPF parameterizations


// Helper functions
void scaleError(TH1D *h, double scale); // scale histogram error
void toPercentage(TH1D *h); // change relative response to percentage change
void fullSimShape(string mode);

// Main function call
void fullSimShapes() {

  fullSimShape("Rjet");
  fullSimShape("chf");
  fullSimShape("nhf");
  fullSimShape("nef");
  /*
  fullSimShape("hp3");
  fullSimShape("hc3");
  fullSimShape("hm3");
  fullSimShape("em3");
  fullSimShape("pm3");
  fullSimShape("tm3");
  */
}

void fullSimShape(string mode) {
  
  setTDRStyle();
  TDirectory *curdir = gDirectory;

  // Open input file
  //TFile *f = new TFile("rootfiles/FullSim_100k_variations_v3.root","READ");
  //TFile *f = new TFile("rootfiles/FullSim_100k_variations_v4.root","READ");
  //TFile *f = new TFile("rootfiles/FullSim_100k_variations_v5.root","READ");
  //TFile *f = new TFile("rootfiles/FullSim_100k_variations_v6.root","READ");
  TFile *f = new TFile("rootfiles/FullSim_100k_variations_v7.root","READ");
  assert(f && !f->IsZombie());

  // Open input file for toyPF placeholders
  TFile *ftoy = new TFile("rootfiles/varPlots_Mikael_5M_20200604_Zjet.root","READ");
  assert(ftoy && !ftoy->IsZombie());

  curdir->cd();

  // Normal mode is to loop over systematics (sysMode==true)
  // Alternative will loop over observables
  bool sysMode = (mode=="Rjet" || mode=="chf" || mode=="nhf" || mode=="nef");

  // List variations to be plotted
  vector<string> vars;
  if (sysMode) {
    vars.push_back("hp3"); // HCAL +3%
    //vars.push_back("hx3"); // HCAL +/-3%
    //vars.push_back("hc1"); // HCAL Custom #1
    vars.push_back("hr3"); // HCAL Red
    vars.push_back("hg3"); // HCAL Green
    vars.push_back("hm3"); // HCAL -3%
    //vars.push_back("em3"); // ECAL Had. -3%
    //vars.push_back("pm3"); // Photons -3%
    //vars.push_back("tm3"); // Tracking -3%
  }   
  else {
    vars.push_back("Rjet");
    vars.push_back("chf");
    vars.push_back("nhf");
    vars.push_back("nef");
  }

  // Map systematics/observables to a function (shape) to be fitted
  // Code in initial guess with parameter representing difference to it
  string slogpol2 = "[0]+log(x)*([1]+log(x)*[2])";
  const char *clogpol2 = slogpol2.c_str();
  string slogpol4 = "[0]+log(x)*([1]+log(x)*([2]+log(x)*([3]+log(x)*[4])))";
  const char *clogpol4 = slogpol4.c_str();
  //string slogpol5 = "[0]+log(x)*([1]+log(x)*([2]+log(x)*([3]+log(x)*"
  //"([4]+log(x)*[5]))))";
  //const char *clogpol5 = slogpol5.c_str();
  string slogpol6 = "[0]+log(x)*([1]+log(x)*([2]+log(x)*([3]+log(x)*"
    "([4]+log(x)*([5]+log(x)*[6])))))";
  const char *clogpol6 = slogpol6.c_str();
  map<string,const char*> func;
  // systematics
  //func["hp3"] = "max(+(1.6+[0])-(1.0+[1])*pow(x,-(0.3+[2])),+abs(0.15+[3]))";
  //func["hp3"] = "sqrt(pow(max((1.6+[0])-(1.0+[1])*pow(x,-(0.3+[2])),0.),2)"
  //"+pow(0.15+[3],2))";
  func["hp3"] = clogpol4;
  //func["hx3"] = "(-1 + log(x/15.)/log(208./15.))*"
  //"max(+(1.6+[0])-(1.0+[1])*pow(x,-(0.3+[2])),+abs(0.15+[3]))";
  func["hx3"] = "sqrt(pow(max((1.6+[0])-(1.0+[1])*pow(x,-(0.3+[2])),0.),2)"
    "+pow(0.15+[3],2)) * (-1 + log(x/15.)/log(208./15.))";
  //func["hm3"] = "min(-(1.6+[0])+(1.0+[1])*pow(x,-(0.3+[2])),-abs(0.15+[3]))";
  func["hc1"] = "[0]+log(x)*([1]+log(x)*([2]+log(x)*([3]+log(x)*[4])))";
  //func["hc3"] = "[0]+log(x)*([1]+log(x)*([2]+log(x)*([3]+log(x)*[4])))";
  func["hg3"] = clogpol6;
  func["hr3"] = clogpol6;
  //func["hm3"] = "-sqrt(pow(max((+1.6+[0])-(1.0+[1])*pow(x,-(0.3+[2])),0.),2)"
  //"+pow(0.15+[3],2))";
  func["hm3"] = clogpol4;
  //func["em3"] = "min(-(0.6+[0])+(0.4+[1])*pow(x,-(0.3+[2])),-abs(0.06+[3]))";
  func["em3"] = clogpol6;
  //func["em3"] = "-sqrt(pow(max((+0.6+[0])-(0.4+[1])*pow(x,-(0.3+[2])),0.),2)"
  //"+pow(0.06+[3],2))";
  //func["pm3"] = "-0.75+[0]";
  func["pm3"] = clogpol2;
  func["tm3"] = "[0]-(1.0+[1])*pow(x,-(0.3+[2]))";
  // observables
  func["Rjet"] = "[0]";
  func["chf"] = "[0]";
  func["nhf"] = "[0]";
  func["nef"] = "[0]";

  map<string, map<string, const char*> > funcs;
  funcs["chf"]["hp3"] = clogpol4;
  funcs["chf"]["hr3"] = clogpol4;
  funcs["chf"]["hg3"] = clogpol6;
  funcs["chf"]["hm3"] = clogpol4;
  funcs["chf"]["em3"] = clogpol4;
  funcs["chf"]["pm3"] = clogpol2;
  funcs["chf"]["tm3"] = clogpol4;
  //
  funcs["nhf"]["hp3"] = clogpol6;
  funcs["nhf"]["hr3"] = clogpol6;
  funcs["nhf"]["hg3"] = clogpol6;
  funcs["nhf"]["hm3"] = clogpol6;
  funcs["nhf"]["em3"] = clogpol4;
  funcs["nhf"]["pm3"] = clogpol2;
  funcs["nhf"]["tm3"] = clogpol4;
  //
  funcs["nef"]["hp3"] = clogpol6;
  funcs["nef"]["hr3"] = clogpol6;
  funcs["nef"]["hg3"] = clogpol6;
  funcs["nef"]["hm3"] = clogpol6;
  funcs["nef"]["em3"] = clogpol4;
  funcs["nef"]["pm3"] = clogpol2;
  funcs["nef"]["tm3"] = clogpol4;

  // Map old toyPF parameterizations
  map<string, map<string, TF1*> > toyf;
  double x[1], p[9];
  //jesFit(x,p); // initialize funcs
  //setToyShapeFuncs(); // initialize funcs
  setFullShapeFuncs(); // initialize funcs
  toyf["Rjet"]["hp3"] = 0;//fhx; assert(fhx);
  toyf["Rjet"]["hx3"] = fhx; assert(fhx);
  toyf["Rjet"]["hc1"] = fhx; assert(fhx);
  toyf["Rjet"]["hr3"] = fhx; assert(fhx);
  toyf["Rjet"]["hg3"] = fhx; assert(fhx);
  toyf["Rjet"]["hm3"] = fhh; assert(fhh);
  toyf["Rjet"]["em3"] = feh; assert(feh);
  toyf["Rjet"]["pm3"] = fp;  assert(fp);
  toyf["Rjet"]["tm3"] = ftd; assert(ftd);

  toyf["chf"]["hp3"] = 0;
  toyf["chf"]["hx3"] = _mpf["chf"][2];
  toyf["chf"]["hc1"] = _mpf["chf"][2];
  toyf["chf"]["hr3"] = _mpf["chf"][2];
  toyf["chf"]["hg3"] = _mpf["chf"][2];
  toyf["chf"]["hm3"] = _mpf["chf"][3];
  toyf["chf"]["em3"] = _mpf["chf"][4];
  toyf["chf"]["pm3"] = _mpf["chf"][1];
  toyf["chf"]["tm3"] = _mpf["chf"][0]; 
  //
  toyf["nhf"]["hp3"] = 0;
  toyf["nhf"]["hx3"] = _mpf["nhf"][2];
  toyf["nhf"]["hc1"] = _mpf["nhf"][2];
  toyf["nhf"]["hr3"] = _mpf["nhf"][2];
  toyf["nhf"]["hg3"] = _mpf["nhf"][2];
  toyf["nhf"]["hm3"] = _mpf["nhf"][3];
  toyf["nhf"]["em3"] = _mpf["nhf"][4];
  toyf["nhf"]["pm3"] = _mpf["nhf"][1];
  toyf["nhf"]["tm3"] = _mpf["nhf"][0]; 
  //
  toyf["nef"]["hp3"] = 0;
  toyf["nef"]["hx3"] = _mpf["nef"][2];
  toyf["nef"]["hc1"] = _mpf["nef"][2];
  toyf["nef"]["hr3"] = _mpf["nef"][2];
  toyf["nef"]["hg3"] = _mpf["nef"][2];
  toyf["nef"]["hm3"] = _mpf["nef"][3];
  toyf["nef"]["em3"] = _mpf["nef"][4];
  toyf["nef"]["pm3"] = _mpf["nef"][1];
  toyf["nef"]["tm3"] = _mpf["nef"][0]; 

  // Map new fullSimShapes
  map<string, map<string, TF1*> > fits;

  // Map systematics/observables to histogram names
  map<string,const char*> name;
  // systematics
  name["hp3"] = "HadHCALp3";
  name["hx3"] = "HadHCALm3";
  name["hc1"] = "customHCALgreen";
  name["hg3"] = "customHCALgreenNptcl";
  name["hr3"] = "customHCALredNptcl";
  name["hm3"] = "HadHCALm3";
  //name["em3"] = "HadECALm3"; // toyPF placeholder
  name["em3"] = "ECALm3"; // toyPF placeholder
  name["pm3"] = "Photonm3";  // toyPF placeholder
  name["tm3"] = "Trkm3";
  // observables
  name["Rjet"] = "Rjet";
  name["chf"] = "chf";
  name["nhf"] = "nhf";
  name["nef"] = "gammaf";

  // Maps variations/observables to labels
  map<string,const char*> label;
  // systematics
  label["hp3"] = "HCAL Had. +3%";
  label["hx3"] = "HCAL Had. #pm3%";
  label["hc1"] = "HCAL Custom #1";
  label["hr3"] = "HCAL Red";
  label["hg3"] = "HCAL Green";
  label["hm3"] = "HCAL Had. -3%";
  label["em3"] = "ECAL Had. -3% (hybrid)";//(toyPF)";
  label["pm3"] = "Photons -3% (toyPF)";
  //label["tm3"] = "Tracking -1%";//patched from -3%";
  label["tm3"] = "Tracking -3%";
  // observables
  label["Rjet"] = "Jet response (R_{jet})";
  label["chf"] = "Charged hadrons (CHF)";
  label["nhf"] = "Neutral hadrons (NHF)";
  label["nef"] = "Photons (NEF)";

  // Map systematics/observables to marker style
  map<string, int> marker;
  // systematics
  marker["hp3"] = kFullCircle;
  marker["hx3"] = kFullCircle;
  marker["hc1"] = kFullStar;
  marker["hr3"] = kFullStar;
  marker["hg3"] = kFullStar;
  marker["hm3"] = kFullCircle;
  marker["em3"] = kOpenTriangleUp;
  marker["pm3"] = kFullDiamond;
  marker["tm3"] = kFullSquare;
  // observables
  marker["Rjet"] = kFullCircle;
  marker["chf"] = kFullCircle;
  marker["nhf"] = kFullDiamond;
  marker["nef"] = kFullSquare;

  // Map variations/observables to marker color
  map<string, int> color;
  // systematics
  color["hp3"] = kRed+1;
  color["hx3"] = kBlack;
  color["hc1"] = kBlack;
  color["hr3"] = kRed;
  color["hg3"] = kGreen+2;
  color["hm3"] = kBlue+1;
  color["em3"] = kBlue-9;
  color["pm3"] = kCyan+2;
  color["tm3"] = kBlack;//kGreen+2;
  // observables
  color["Rjet"] = kBlack;
  color["chf"] = kRed;
  color["nhf"] = kGreen+2;
  color["nef"] = kBlue;
  
  // Load histograms
  map<string, TH1D*> hist;
  for (int i = 0; i != vars.size(); ++i) {

    string sv = vars[i];
    const char *cobs = (sysMode ? name[mode] : name[sv]);
    const char *csys = (sysMode ? name[sv] : name[mode]);
    string obs = cobs;
    string sys = csys;

    TH1D *hv = (TH1D*)f->Get(Form("%s_%s",cobs,csys));
    if (hv) {
      // Clone to avoid changing original
      hv = (TH1D*)hv->Clone(Form("%s_%s_%s",mode.c_str(),cobs,csys));

      scaleError(hv,0.1); // Patch FullSim uncertainty
      //if (obs=="Rjet") toPercentage(hv);
      //else {
      //hv->Scale(100.); scaleError(hv,0.1);
      //}
    }
    else { // toyPF placeholders
      hv = (TH1D*)ftoy->Get(Form("h_%s_%s",cobs,csys));
      if (!hv) hv = (TH1D*)ftoy->Get(Form("h%s_%s",cobs,csys));
      if (!hv) cout << "Hist " << sv << " not found!" << endl << flush;
      assert(hv);

      // Clone to avoid changing original
      hv = (TH1D*)hv->Clone(Form("%s_%s_%s",mode.c_str(),cobs,csys));
    }

    if (obs=="Rjet") toPercentage(hv);
    else {
      hv->Scale(100.); scaleError(hv,0.1);
      if (obs=="nef") scaleError(hv,0.1);
    }

    // Reduce tracking effect on composition to fit on plots
    //if (sys=="Trkm3") {
    //hv->Scale(1./3);
    //}
    hist[sv] = (TH1D*)hv;
  } // for i in vars

  // Patch ECALm3 to HadECALm3 by subtracting Photonm3
  if (patchECALm3 && sysMode) {
    TH1D *he = hist["em3"]; assert(he);
    TH1D *hp = hist["pm3"]; assert(hp);
    TH1D *h = (TH1D*)he->Clone("hem3");
    for (int i = 1; i != h->GetNbinsX()+1; ++i) {
      //h->SetBinContent(i, min(he->GetBinContent(i)-hp->GetBinContent(i),0.));
      h->SetBinContent(i, (he->GetBinContent(i)-hp->GetBinContent(i)));
      h->SetBinError(i, sqrt(pow(he->GetBinError(i),2) +
			     pow(hp->GetBinError(i),2)));
    } // for i
    hist["em3"] = h;
  } // patchECALm3

  // Calculate HCAL cross
  if (true && hist["hp3"]) {
    // Interpolate SPR
    TH1D *hp3 = hist["hp3"];
    TH1D *hx3 = (TH1D*)hp3->Clone("hx3");

    // Log-lin interpolated SPR from -3% to +3%
    for (int i = 1; i != hp3->GetNbinsX()+1; ++i) {
      double pt = hp3->GetBinCenter(i);
      double p = hp3->GetBinContent(i);
      // Log-lin interpolation from -1 @15 GeV to +1 @2884 GeV (0 @208 GeV)
      double w = -1 + log(pt/15.)/log(208./15.);
      double x = log(pt/208.);
      hx3->SetBinContent(i, w*hp3->GetBinContent(i));
      hx3->SetBinError(i, hp3->GetBinError(i));
    } // for i
    hist["hx3"] = hx3;
  }

  // Setup plotting
  double maxy = (mode=="Rjet" ? 3 : 3-1e-5);//2.5-1e-5);
  double miny = (mode=="Rjet" ? -2 : -2+1e-5);
  const char *title = (mode=="Rjet" ? "Response change (%)" :
		       sysMode ? "PF composition change (10^{-2})" :
		       "PF changes (% or 10^{-2})");
  TH1D *h = tdrHist(Form("h_%s",mode.c_str()), title,miny,maxy);
  lumi_13TeV = "FullSim (+toyPF placeholders)";
  TCanvas *c1 = tdrCanvas(Form("c1_%s",mode.c_str()),h,4,11,kSquare);
  gPad->SetLogx();
  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(h->GetXaxis()->GetXmin(),0,h->GetXaxis()->GetXmax(),0);
  
  // Setup legends
  //TLegend *leg = tdrLeg(0.38,0.91-0.045*vars.size(),0.68,0.91);
  TLegend *leg = tdrLeg(0.38,0.91-0.040*vars.size(),0.68,0.91);
  leg->SetTextSize(0.040);
  TLatex *tex = new TLatex();
  tex->SetTextSize(0.045);
  tex->SetNDC();
  //tex->DrawLatex(0.78,0.87,"|#eta| < 1.3");
  if (mode=="Rjet") 
    tex->DrawLatex(0.19,0.75,"|#eta| < 1.3");
  else
    tex->DrawLatex(0.19,0.75,mode.c_str());

  // Plot histograms
  for (int i = 0; i != vars.size(); ++i) {

    string sv = vars[i];
    const char *cobs = (sysMode ? name[mode] : name[sv]);
    const char *csys = (sysMode ? name[sv] : name[mode]);
    string obs = cobs;
    string sys = csys;

    tdrDraw(hist[sv], "Pz", marker[sv], color[sv], kSolid, -1, kNone);
    leg->AddEntry(hist[sv], label[sv], "PL");

    if (mode=="Rjet") {
    assert(func[sv]!="");
      TF1 *f1 = new TF1(Form("f1_%s_%s",cobs,csys),func[sv],15.,3500.);
      f1->SetParameters(0,0,0,0);
      if (sv=="hc3") hist[sv]->Fit(f1,"QRNW");
      else hist[sv]->Fit(f1,"QRN");
      f1->SetLineColor(color[sv]);
      f1->SetLineWidth(2);
      f1->Draw("SAME");

      fits[mode][sv] = f1;
    }
    else if (funcs[mode][sv]!=0) {
      TF1 *f1 = new TF1(Form("f1_%s_%s_%s",mode.c_str(),cobs,csys),
			funcs[mode][sv],15.,3500.);
      if (sv=="tm3"&&mode=="chf") f1->SetRange(25.,3500.);
      f1->SetParameters(0,0,0,0);
      //hist[sv]->Fit(f1,"QRN");
      if (sv=="hc3") hist[sv]->Fit(f1,"QRNW");
      else hist[sv]->Fit(f1,"QRN");
      f1->SetLineColor(color[sv]);
      f1->SetLineWidth(2);
      f1->Draw("SAME");

      fits[mode][sv] = f1;
    }

    if (plotToyPF && toyf[mode][sv]!=0) {
      TF1 *f1 = toyf[mode][sv];
      f1->SetLineColor(color[sv]);
      f1->SetLineWidth(3);
      f1->SetLineStyle(kDashed);
      f1->Draw("SAME");
    }
  } // for i

  //c1->SaveAs(Form("pdf/fullSimShapes/fullSimShapes_%s.pdf",mode.c_str()));
  c1->SaveAs(Form("pdf/fullSimShapes/fullSimShapes_%s_%s.pdf",
		  sysMode ? "sys" : "var", mode.c_str()));

  // Print fullSimShapes to be used in globalFitL3Res.C
  if (sysMode) {
    cout << "    // Fits from minitools/fullSimShapes.C" << endl;
    cout << "    //////////////////////////////////////" << endl;
    cout << endl;

    if (mode=="Rjet") cout << "    // Jet response (Rjet)\n";
    if (mode=="chf")  cout << "    // Charged hadron fraction (CHF)\n";
    if (mode=="nhf")  cout << "    // Neutral hadron fraction (NHF)\n";
    if (mode=="nef")  cout << "    // Photon fraction (NEF)\n";
    cout << "   // Fits from minitools/fullSimShapes.C" << endl;

    for (int i = 0; i != vars.size(); ++i) {
      string sv = vars[i];
      const char *cv = sv.c_str();
      const char *cm = mode.c_str();
      TF1 *f1 = fits[mode][sv]; //assert(f1);
      if (!f1) continue;

// [p0]+[p1]*(1+(pow(x/[p2],[p3])-1)/(pow(x/[p2],[p3])+1))+[p4]*pow(x,[p5])
      cout << Form("    TF1 *f%s_%s = new TF1(\"f%s_%s\",\"%s\",15,4500);",
		   cv,cm,cv,cm,f1->GetExpFormula().Data()) << endl;
      cout << Form("    f%s_%s->SetParameters(",cv,cm);
      for (int j = 0; j != f1->GetNpar(); ++j) {
	cout << Form("%s%1.4g",j==0 ? "" : ",",f1->GetParameter(j));
      } // for j in GetNpar
      cout << Form("); // %1.1f/%d\n",f1->GetChisquare(),f1->GetNDF());
    } // for in in vars
    cout << endl;
  } // if sysMode

} // fullSimShapes


// Scale histogram error
void scaleError(TH1D *h, double scale) {
  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    h->SetBinError(i, h->GetBinError(i)*scale);
  }
} // scaleError

// Change relative response to percentage change
void toPercentage(TH1D *h) {
  
  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    h->SetBinContent(i, 100.*(h->GetBinContent(i)-1));
    h->SetBinError(i, 100.*h->GetBinError(i));
  }
}
