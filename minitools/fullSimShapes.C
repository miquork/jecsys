// Purpose: Fit response and composition variation shapes from FullSim.
//          This is evolution of minitools/varPlots.C, dropping old baggage
//          Comparisons are done to old shapes from toyPF
//          (use ToyPF as placeholder for missing plots)
#include "TFile.h"
#include "TF1.h"
#include "TSystem.h"
#include "TGraphErrors.h"

#include "../tdrstyle_mod15.C"
#include "../tools.C"

#include <fstream>
#include <vector>
#include <string>

// Settings
bool patchECALm3 = false;//true;
bool plotToyPF = true;
bool doProduction = true; // produce lines for globalFitSettings.h
string CorLevel = "L1L2L3";
#include "../globalFitL3Res.C" // toyPF parameterizations


// Helper functions
void scaleError(TH1D *h, double scale); // scale histogram error
void toPercentage(TH1D *h); // change relative response to percentage change
void scaleGraph(TGraph *g, double scale);
void fullSimShape(string mode);
void fullSimFlavor(string obs, string var);

// Main function call
void fullSimShapes() {

  /*
  fullSimShape("Rjet");
  fullSimShape("chf");
  fullSimShape("nhf");
  fullSimShape("nef");
  */
  /*
  Fullsimshape("hp3");
  fullSimShape("hc3");
  fullSimShape("hm3");
  fullSimShape("em3");
  fullSimShape("pm3");
  fullSimShape("tm3");
  */
  /*
  fullSimShape("tm3");
  fullSimShape("tv2");
  //fullSimShape("tv3");
  fullSimShape("tv2b");
  fullSimShape("tv2c");
  fullSimShape("tv3b");

  fullSimShape("cmb");
  */

  //fullSimShape("tm3");
  //fullSimShape("tv2c");
  //fullSimShape("tv3b");
  //fullSimShape("hg3");
  //fullSimShape("hr3");
  //fullSimShape("pm3");
  //fullSimShape("hm3");
  //fullSimShape("em3");
  //fullSimShape("hb3");
  //fullSimShape("hb1");

  //fullSimShape("tv4");
  fullSimShape("tv402");
  //fullSimShape("tv404");
  //fullSimShape("tv410");
  //fullSimShape("tv430");
  //fullSimShape("tv405");
  //fullSimShape("tv416");
  //fullSimShape("tv420");

  //fullSimShape("tv3n1");
  //fullSimShape("tv300pn");

  //fullSimShape("hhm3");
  //fullSimShape("hhp3");
  //fullSimShape("hhred103");
  //fullSimShape("hhred100");
  //fullSimShape("hhred097");
  //fullSimShape("hhblue103");
  //fullSimShape("hhblue100");
  //fullSimShape("hhblue097");

  //fullSimShape("em3");
  //fullSimShape("pm3");

  //fullSimFlavor("chf","Trkv3b");
  //fullSimFlavor("Rjet","Trkv3b");
  //fullSimFlavor("nhf","customHCALgreenNptcl");
  //fullSimFlavor("nef","customHCALgreenNptcl");
  //fullSimFlavor("Rjet","customHCALgreenNptcl");
}

void fullSimShape(string mode) {
  
  setTDRStyle();
  TDirectory *curdir = gDirectory;

  // Open input file
  //TFile *f = new TFile("rootfiles/FullSim_100k_variations_v3.root","READ");
  //TFile *f = new TFile("rootfiles/FullSim_100k_variations_v4.root","READ");
  //TFile *f = new TFile("rootfiles/FullSim_100k_variations_v5.root","READ");
  //TFile *f = new TFile("rootfiles/FullSim_100k_variations_v6.root","READ");
  //TFile *f = new TFile("rootfiles/FullSim_100k_variations_v7.root","READ");
  //TFile *f = new TFile("rootfiles/FullSim_100k_variations_v8.root","READ");
  //TFile *f = new TFile("rootfiles/FullSim_100k_variations_v9.root","READ");
  //TFile *f = new TFile("rootfiles/FullSim_100k_variations_v10.root","READ");
  //TFile *f = new TFile("rootfiles/FullSim_100k_variations_v11.root","READ");
  //TFile *f = new TFile("rootfiles/FullSim_100k_variations_v12.root","READ");
  TFile *f = new TFile("rootfiles/FullSim_100k_variations_v15.root","READ");
  assert(f && !f->IsZombie());

  // Open input file for toyPF placeholders
  //TFile *ftoy = new TFile("rootfiles/varPlots_Mikael_5M_20200604_Zjet.root","READ");
  //assert(ftoy && !ftoy->IsZombie());

  curdir->cd();

  // Normal mode is to loop over systematics (sysMode==true)
  // Alternative will loop over observables
  bool sysMode = (mode=="Rjet" || mode=="chf" || mode=="nhf" || mode=="nef");

  // List variations to be plotted
  vector<string> vars;
  vector<pair<string,double> > cmbs;
  if (mode=="cmb") { // cmbMode
    vars.push_back("Rjet");
    vars.push_back("chf");
    vars.push_back("nhf");
    vars.push_back("nef");
    /*
    cmbs.push_back(make_pair<string,double>("hg3",1.0)); // HCAL Green
    cmbs.push_back(make_pair<string,double>("hm3",0.5)); // HCAL -3%
    //cmbs.push_back(make_pair<string,double>("em3",0.4)); // ECAL -3%
    cmbs.push_back(make_pair<string,double>("em3",0.2)); // ECAL -3%
    //cmbs.push_back(make_pair<string,int>("pm3",0.1)); // Photons -3%
    //cmbs.push_back(make_pair<string,double>("tm3",0.7)); // Tracking -3%
    //cmbs.push_back(make_pair<string,double>("tv2",1.0)); // Tracks -1% to -3%
    cmbs.push_back(make_pair<string,double>("tv2b",1.0)); // Tracks -1% to -3%
    */
    // (x*3.0+y*3)=3; (x*0.5+y*3)=1 => x*2.5=2, y=1-0.5*x => x=0.8, y=0.2
    //cmbs.push_back(make_pair<string,double>("tm3",0.2)); // Trk -3%
    //cmbs.push_back(make_pair<string,double>("tv2c",0.8)); // Trk -0.5% to -3%
    //cmbs.push_back(make_pair<string,double>("tv2b",-1.0)); // Trk -1.0% to -3%
    // (x*3.0+y*3)=3; (x*1.0+y*3)=0 => x*2.0=3, y=-x/3 => x=1.5, y=-0.5
    //cmbs.push_back(make_pair<string,double>("tv2b",+1.5)); // Trk -1.0% to -3%
    //cmbs.push_back(make_pair<string,double>("tm3",-0.5)); // Trk -3%
    //cmbs.push_back(make_pair<string,double>("tv2",+1.5)); // Trk -1.0% to -3%
    //cmbs.push_back(make_pair<string,double>("tm3",-0.5)); // Trk -3%
    // (x*3.0+y*3)=3; (x*0.5+y*3)=0 => x*2.5=3, y=-0.5*x/3 => x=1.2, y=-0.2
    //cmbs.push_back(make_pair<string,double>("tv2c",+1.2)); // Trk -1.0% to -3%
    //cmbs.push_back(make_pair<string,double>("tm3",-0.2)); // Trk -3%
    cmbs.push_back(make_pair<string,double>("tv2",+1.0)); // Trk -1.0% to -3%
    cmbs.push_back(make_pair<string,double>("tv2b",-1.0)); // Trk -1.0% to -3%
  }
  else if (sysMode) { // sysMode
    /*
    vars.push_back("hp3"); // HCAL +3%
    vars.push_back("hx3"); // HCAL +/-3%
    //vars.push_back("hc1"); // HCAL Custom #1
    vars.push_back("hr3"); // HCAL Red
    vars.push_back("hg3"); // HCAL Green
    vars.push_back("hb3"); // HCAL Blue097
    vars.push_back("hb1"); // HCAL Blue1
    vars.push_back("hm3"); // HCAL -3%
    */
    //vars.push_back("em3"); // ECAL Had. -3%
    //vars.push_back("pm3"); // Photons -3%
    /*
    vars.push_back("tm3"); // Tracking -3%
    vars.push_back("tv2"); // Tracking v2
    vars.push_back("tv2b"); // Tracking v2b
    vars.push_back("tv2c"); // Tracking v2c
    //vars.push_back("tv3"); // Tracking v3
    vars.push_back("tv3b"); // Tracking v3b
    */
    /*
    // Tracking family
    vars.push_back("tv4"); // Tracking -3%
    vars.push_back("tv402"); // Tracking -3% at ntrk>=2 
    vars.push_back("tv404"); // Tracking -3% at ntrk>=4 
    vars.push_back("tv405"); // Tracking -3% at ntrk>=5 
    vars.push_back("tv410"); // Tracking -3% at ntrk>=10 
    vars.push_back("tv416"); // Tracking -3% at ntrk>=16 
    vars.push_back("tv420"); // Tracking -3% at ntrk>=20 
    vars.push_back("tv430"); // Tracking -3% at ntrk>=30 
    */
    // Second tracking family
    vars.push_back("tv4"); // Tracking -3%
    vars.push_back("tv3n1"); // Tracking -3% at ntrk==1 
    vars.push_back("tv402"); // Tracking -3% at ntrk>=2
    vars.push_back("tv300pn"); // Tracking 0.999^{n-1} = 1.000*0.999^{n-1} (3c)
    vars.push_back("tv301pn"); // Tracking 0.999^n = 0.999*0.999^{n-1} (3d)
    vars.push_back("tv310pn"); // Tracking 0.990*0.999^{n-1} (3b)

    /*
    // HCAL families
    vars.push_back("hhm3"); // HCAL -3%
    vars.push_back("hhp3"); // HCAL +3%
    vars.push_back("hhred103"); // HCAL red fit +3%
    vars.push_back("hhred100"); // HCAL red fit 0%
    vars.push_back("hhred097"); // HCAL red fit -3%
    vars.push_back("hhblue103"); // HCAL blue fit +3% (green)
    vars.push_back("hhblue100"); // HCAL blue fit 0%
    vars.push_back("hhblue097"); // HCAL blue fit -3%
    */
  }   
  else { // obsMode
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
  string slogpol6_20 = "[0]+log(max(20.,x))*([1]+log(max(20.,x))*([2]+log(max(20.,x))*([3]+log(max(20.,x))*([4]+log(max(20.,x))*([5]+log(max(20.,x))*[6])))))";
  const char *clogpol6_20 = slogpol6_20.c_str();
  string slogpol6_30 = "[0]+log(max(30.,x))*([1]+log(max(30.,x))*([2]+log(max(30.,x))*([3]+log(max(30.,x))*([4]+log(max(30.,x))*([5]+log(max(30.,x))*[6])))))";
  const char *clogpol6_30 = slogpol6_30.c_str();
  string slogpol6_40 = "[0]+log(max(40.,x))*([1]+log(max(40.,x))*([2]+log(max(40.,x))*([3]+log(max(40.,x))*([4]+log(max(40.,x))*([5]+log(max(40.,x))*[6])))))";
  const char *clogpol6_40 = slogpol6_40.c_str();
  string slogpol6_45 = "[0]+log(max(45.,x))*([1]+log(max(45.,x))*([2]+log(max(45.,x))*([3]+log(max(45.,x))*([4]+log(max(45.,x))*([5]+log(max(45.,x))*[6])))))";
  const char *clogpol6_45 = slogpol6_45.c_str();
  string slogpol6_50 = "[0]+log(max(50.,x))*([1]+log(max(50.,x))*([2]+log(max(50.,x))*([3]+log(max(50.,x))*([4]+log(max(50.,x))*([5]+log(max(50.,x))*[6])))))";
  const char *clogpol6_50 = slogpol6_50.c_str();
  string slogpol6_55 = "[0]+log(max(55.,x))*([1]+log(max(55.,x))*([2]+log(max(55.,x))*([3]+log(max(55.,x))*([4]+log(max(55.,x))*([5]+log(max(55.,x))*[6])))))";
  const char *clogpol6_55 = slogpol6_55.c_str();
  string slogpol6_60 = "[0]+log(max(60.,x))*([1]+log(max(60.,x))*([2]+log(max(60.,x))*([3]+log(max(60.,x))*([4]+log(max(60.,x))*([5]+log(max(60.,x))*[6])))))";
  const char *clogpol6_60 = slogpol6_60.c_str();
  string slogpol6_65 = "[0]+log(max(65.,x))*([1]+log(max(65.,x))*([2]+log(max(65.,x))*([3]+log(max(65.,x))*([4]+log(max(65.,x))*([5]+log(max(65.,x))*[6])))))";
  const char *clogpol6_65 = slogpol6_65.c_str();
  string slogpol6_70 = "[0]+log(max(70.,x))*([1]+log(max(70.,x))*([2]+log(max(70.,x))*([3]+log(max(70.,x))*([4]+log(max(70.,x))*([5]+log(max(70.,x))*[6])))))";
  const char *clogpol6_70 = slogpol6_70.c_str();
  string slogpol6_80 = "[0]+log(max(80.,x))*([1]+log(max(80.,x))*([2]+log(max(80.,x))*([3]+log(max(80.,x))*([4]+log(max(80.,x))*([5]+log(max(80.,x))*[6])))))";
  const char *clogpol6_80 = slogpol6_80.c_str();
  string slogpol6_100 = "[0]+log(max(100.,x))*([1]+log(max(100.,x))*([2]+log(max(100.,x))*([3]+log(max(100.,x))*([4]+log(max(100.,x))*([5]+log(max(100.,x))*[6])))))";
  const char *clogpol6_100 = slogpol6_100.c_str();
  string slogpol6_130 = "[0]+log(max(130.,x))*([1]+log(max(130.,x))*([2]+log(max(130.,x))*([3]+log(max(130.,x))*([4]+log(max(130.,x))*([5]+log(max(130.,x))*[6])))))";
  const char *clogpol6_130 = slogpol6_130.c_str();
  string slogpol6_210 = "[0]+log(max(210.,x))*([1]+log(max(210.,x))*([2]+log(max(210.,x))*([3]+log(max(210.,x))*([4]+log(max(210.,x))*([5]+log(max(210.,x))*[6])))))";
  const char *clogpol6_210 = slogpol6_210.c_str();
  string slogpol6_260 = "[0]+log(max(260.,x))*([1]+log(max(260.,x))*([2]+log(max(260.,x))*([3]+log(max(260.,x))*([4]+log(max(260.,x))*([5]+log(max(260.,x))*[6])))))";
  const char *clogpol6_260 = slogpol6_260.c_str();
  string slogpol6_400 = "[0]+log(max(400.,x))*([1]+log(max(400.,x))*([2]+log(max(400.,x))*([3]+log(max(400.,x))*([4]+log(max(400.,x))*([5]+log(max(400.,x))*[6])))))";
  const char *clogpol6_400 = slogpol6_400.c_str();
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
  func["hr3"] = clogpol6_60;
  func["hb3"] = clogpol6_60;
  func["hb1"] = clogpol6_60;
  //func["hm3"] = "-sqrt(pow(max((+1.6+[0])-(1.0+[1])*pow(x,-(0.3+[2])),0.),2)"
  //"+pow(0.15+[3],2))";
  func["hm3"] = clogpol6_30;
  //func["em3"] = "min(-(0.6+[0])+(0.4+[1])*pow(x,-(0.3+[2])),-abs(0.06+[3]))"
  func["em3"] = clogpol6_60;
  //func["em3"] = "-sqrt(pow(max((+0.6+[0])-(0.4+[1])*pow(x,-(0.3+[2])),0.),2)"
  //"+pow(0.06+[3],2))";
  //func["pm3"] = "-0.75+[0]";
  func["pm3"] = clogpol2;
  //func["tm3"] = "[0]-(1.0+[1])*pow(x,-(0.3+[2]))";
  //func["tv2"] = clogpol4;
  func["tv2"] = "[0]-(1.0+[1])*pow(x,-(0.3+[2]))+[3]/log(x)";
  //func["tv3"] = clogpol4;
  func["tv3"] = "[0]-(1.0+[1])*pow(x,-(0.3+[2]))+[3]/log(x)";
  //func["tv3b"] = clogpol4;
  func["tv3b"] = "[0]-(1.0+[1])*pow(x,-(0.3+[2]))+[3]/log(x)";
  //func["tv2b"] = clogpol4; // [25,3500 GeV]
  func["tv2b"] = "[0]-(1.0+[1])*pow(x,-(0.3+[2]))+[3]/log(x)";
  //func["tv2c"] = clogpol4;
  //func["tv2c"] = "[0]-(1.0+[1])*pow(x,-(0.3+[2]))+[3]/log(x)";
  //
  func["tv4"] = clogpol6_20;
  func["tv402"] = clogpol6_20;
  func["tv404"] = clogpol6_40;//30;
  func["tv405"] = clogpol6_55;
  func["tv410"] = clogpol6_80;//130;
  func["tv416"] = clogpol6_210;
  func["tv420"] = clogpol6_260;
  func["tv430"] = clogpol6_400;
  //
  func["tv3"] = clogpol6_20;
  func["tv3n1"] = clogpol6_20;
  func["tv300pn"] = clogpol6_60;
  func["tv301pn"] = clogpol6_40;
  func["tv310pn"] = clogpol6_20;
  //
  func["hhp3"] = clogpol6;
  func["hhm3"] = clogpol6;
  func["hhred103"] = clogpol6;
  func["hhred100"] = clogpol6;
  func["hhred097"] = clogpol6;
  func["hhblue103"] = clogpol6;
  func["hhblue100"] = clogpol6;
  func["hhblue097"] = clogpol6;
  //
  func["cmb"] = clogpol6;
  // observables
  func["Rjet"] = "[0]";
  func["chf"] = "[0]";
  func["nhf"] = "[0]";
  func["nef"] = "[0]";

  map<string, map<string, const char*> > funcs;
  funcs["cmb"]["Rjet"] = clogpol6;
  //
  funcs["chf"]["hp3"] = clogpol4;
  funcs["chf"]["hr3"] = clogpol6_60;
  funcs["chf"]["hg3"] = clogpol6;
  funcs["chf"]["hb3"] = clogpol6_60;
  funcs["chf"]["hb1"] = clogpol6_60;
  funcs["chf"]["hm3"] = clogpol6_45;
  funcs["chf"]["em3"] = clogpol6_60;
  funcs["chf"]["pm3"] = clogpol2;
  //funcs["chf"]["tm3"] = clogpol4;
  funcs["chf"]["tv2"] = clogpol6;
  funcs["chf"]["tv3"] = clogpol6;
  funcs["chf"]["tv3b"] = clogpol6_60;
  funcs["chf"]["tv2b"] = clogpol6;
  //funcs["chf"]["tv2c"] = clogpol6;
  //
  funcs["chf"]["tv4"] = clogpol6_20;
  funcs["chf"]["tv402"] = clogpol6_20;
  funcs["chf"]["tv404"] = clogpol6_40;//30;
  funcs["chf"]["tv405"] = clogpol6_55;
  funcs["chf"]["tv410"] = clogpol6_130;
  funcs["chf"]["tv416"] = clogpol6_210;
  funcs["chf"]["tv420"] = clogpol6_260;
  funcs["chf"]["tv430"] = clogpol6_400;
  //
  funcs["chf"]["tv3"] = clogpol6_20;
  funcs["chf"]["tv3n1"] = clogpol6_20;
  funcs["chf"]["tv300pn"] = clogpol6_60;
  funcs["chf"]["tv301pn"] = clogpol6_40;
  funcs["chf"]["tv310pn"] = clogpol6_20;
  //
  funcs["chf"]["hhp3"] = clogpol6;
  funcs["chf"]["hhm3"] = clogpol6;
  funcs["chf"]["hhred103"] = clogpol6;
  funcs["chf"]["hhred100"] = clogpol6;
  funcs["chf"]["hhred097"] = clogpol6;
  funcs["chf"]["hhblue103"] = clogpol6;
  funcs["chf"]["hhblue100"] = clogpol6;
  funcs["chf"]["hhblue097"] = clogpol6;
  //
  funcs["cmb"]["chf"] = clogpol6;
  //
  funcs["nhf"]["hp3"] = clogpol6;
  funcs["nhf"]["hr3"] = clogpol6_60;
  funcs["nhf"]["hg3"] = clogpol6;
  funcs["nhf"]["hb3"] = clogpol6_60;
  funcs["nhf"]["hb1"] = clogpol6_60;
  funcs["nhf"]["hm3"] = clogpol6_45;
  funcs["nhf"]["em3"] = clogpol6_60;
  funcs["nhf"]["pm3"] = clogpol2;
  //funcs["nhf"]["tm3"] = clogpol4;
  funcs["nhf"]["tv2"] = clogpol6;
  funcs["nhf"]["tv3"] = clogpol6;
  funcs["nhf"]["tv3b"] = clogpol6_60;
  funcs["nhf"]["tv2b"] = clogpol6;
  //funcs["nhf"]["tv2c"] = clogpol6;
  //
  funcs["nhf"]["tv4"] = clogpol6_20;
  funcs["nhf"]["tv402"] = clogpol6_20;
  funcs["nhf"]["tv404"] = clogpol6_40;//30;
  funcs["nhf"]["tv405"] = clogpol6_55;
  funcs["nhf"]["tv410"] = clogpol6_130;
  funcs["nhf"]["tv416"] = clogpol6_210;
  funcs["nhf"]["tv420"] = clogpol6_260;
  funcs["nhf"]["tv430"] = clogpol6_400;
  //
  funcs["nhf"]["tv3"] = clogpol6_20;
  funcs["nhf"]["tv3n1"] = clogpol6_20;
  funcs["nhf"]["tv300pn"] = clogpol6_60;
  funcs["nhf"]["tv301pn"] = clogpol6_40;
  funcs["nhf"]["tv310pn"] = clogpol6_20;
  //
  funcs["nhf"]["hhp3"] = clogpol6;
  funcs["nhf"]["hhm3"] = clogpol6;
  funcs["nhf"]["hhred103"] = clogpol6;
  funcs["nhf"]["hhred100"] = clogpol6;
  funcs["nhf"]["hhred097"] = clogpol6;
  funcs["nhf"]["hhblue103"] = clogpol6;
  funcs["nhf"]["hhblue100"] = clogpol6;
  funcs["nhf"]["hhblue097"] = clogpol6;
  //
  funcs["cmb"]["nhf"] = clogpol6;
  //
  funcs["nef"]["hp3"] = clogpol6;
  funcs["nef"]["hr3"] = clogpol6_60;
  funcs["nef"]["hg3"] = clogpol6;
  funcs["nef"]["hb3"] = clogpol6_60;
  funcs["nef"]["hb1"] = clogpol6_60;
  funcs["nef"]["hm3"] = clogpol6_45;
  funcs["nef"]["em3"] = clogpol6_60;
  funcs["nef"]["pm3"] = clogpol2;
  //funcs["nef"]["tm3"] = clogpol4;
  funcs["nef"]["tv2"] = clogpol6;
  funcs["nef"]["tv3"] = clogpol6;
  funcs["nef"]["tv3b"] = clogpol6_60;
  funcs["nef"]["tv2b"] = clogpol6;
  //funcs["nef"]["tv2c"] = clogpol6;
  //
  funcs["nef"]["tv4"] = clogpol6_20;
  funcs["nef"]["tv402"] = clogpol6_20;
  funcs["nef"]["tv404"] = clogpol6_40;//30;
  funcs["nef"]["tv405"] = clogpol6_55;
  funcs["nef"]["tv410"] = clogpol6_130;
  funcs["nef"]["tv416"] = clogpol6_210;
  funcs["nef"]["tv420"] = clogpol6_260;
  funcs["nef"]["tv430"] = clogpol6_400;
  //
  funcs["nef"]["tv3"] = clogpol6_20;
  funcs["nef"]["tv3n1"] = clogpol6_20;
  funcs["nef"]["tv300pn"] = clogpol6_60;
  funcs["nef"]["tv301pn"] = clogpol6_40;
  funcs["nef"]["tv310pn"] = clogpol6_20;
  //
  funcs["nef"]["hhp3"] = clogpol6;
  funcs["nef"]["hhm3"] = clogpol6;
  funcs["nef"]["hhred103"] = clogpol6;
  funcs["nef"]["hhred100"] = clogpol6;
  funcs["nef"]["hhred097"] = clogpol6;
  funcs["nef"]["hhblue103"] = clogpol6;
  funcs["nef"]["hhblue100"] = clogpol6;
  funcs["nef"]["hhblue097"] = clogpol6;
  //
  funcs["cmb"]["nef"] = clogpol6;

  // Hand-optimized sources
  func["tm3"] = "[0]-(5+[1])*pow(x,-(0.3+[2]))";
  funcs["chf"]["tm3"] = "-(6+[0])+(7+[1])*pow(x,-(0.1+[2]))+(6e-4+[3])*pow(x,1+[4])";
  funcs["nef"]["tm3"] = "(3+[0])-(5+[1])*pow(x,-(0.1+[2]))-(1e-2+[3])*pow(x,1+[4])";
  funcs["nhf"]["tm3"] = "(1+[0])-(1+[1])*pow(x,-(0.1+[2]))-(1e-4+[3])*pow(x,1+[4])";
  //
  func["tv2c"] = "[0]-(1+[1])*pow(max(60.,x),-(0.3+[2]))+[3]/log(max(60.,x))";
  funcs["chf"]["tv2c"] = "[0]-(0.2+[1])*pow(max(60.,x),0.3+[2])+(2e-2+[3])*pow(max(60.,x),1+[4])";
    //"-(6+[0])+(7+[1])*pow(x,-(0.1+[2]))+(6e-4+[3])*pow(x,1+[4])-(5+[5])*pow(x,-1+[6])";
  funcs["nef"]["tv2c"] = "[0]+(0.2+[1])*pow(max(60.,x),0.3+[2])-(1e-2+[3])*pow(max(60.,x),1+[4])";
  funcs["nhf"]["tv2c"] = "[0]+(0.1+[1])*pow(max(60.,x),0.3+[2])";
    //"(1+[0])-(1+[1])*pow(x,-(0.1+[2]))-(1e-4+[3])*pow(x,1+[4])+(6+[5])*pow(x,-1+[6])";
  //
  //func["tv3b"] = "[0]-(1+[1])*pow(max(60.,x),-(0.3+[2]))+[3]/log(max(60.,x))";
  //funcs["chf"]["tv3b"] = "[0]-(0.2+[1])*pow(max(60.,x),0.3+[2])+(2e-2+[3])*pow(max(60.,x),1+[4])";
  //funcs["nef"]["tv3b"] = "[0]+(0.2+[1])*pow(max(60.,x),0.3+[2])-(1e-2+[3])*pow(max(60.,x),1+[4])";
  //funcs["nhf"]["tv3b"] = "[0]+(0.1+[1])*pow(max(60.,x),0.3+[2])";

  // Map old toyPF parameterizations
  map<string, map<string, TF1*> > toyf;
  double x[1], p[9];
  //jesFit(x,p); // initialize funcs
  setToyShapeFuncs(); // initialize funcs
  //setFullShapeFuncs(); // initialize funcs
  toyf["Rjet"]["hp3"] = 0;//fhx; assert(fhx);
  toyf["Rjet"]["hx3"] = fhx; assert(fhx);
  toyf["Rjet"]["hc1"] = fhx; assert(fhx);
  toyf["Rjet"]["hr3"] = fhx; assert(fhx);
  toyf["Rjet"]["hg3"] = fhx; assert(fhx);
  toyf["Rjet"]["hb3"] = fhx; assert(fhx);
  toyf["Rjet"]["hb1"] = fhx; assert(fhx);
  toyf["Rjet"]["hm3"] = fhh; assert(fhh);
  toyf["Rjet"]["em3"] = feh; assert(feh);
  toyf["Rjet"]["pm3"] = fp;  assert(fp);
  toyf["Rjet"]["tm3"] = ftd; assert(ftd);
  toyf["Rjet"]["tv2"] = ftd; assert(ftd);
  toyf["Rjet"]["tv3"] = ftd; assert(ftd);
  toyf["Rjet"]["tv3b"] = ftd; assert(ftd);
  toyf["Rjet"]["tv2b"] = ftd; assert(ftd);
  toyf["Rjet"]["tv2c"] = ftd; assert(ftd);
  toyf["Rjet"]["cmb"] = 0;

  toyf["chf"]["hp3"] = 0;
  toyf["chf"]["hx3"] = _mpf["chf"][2];
  toyf["chf"]["hc1"] = _mpf["chf"][2];
  toyf["chf"]["hr3"] = _mpf["chf"][2];
  toyf["chf"]["hg3"] = _mpf["chf"][2];
  toyf["chf"]["hb3"] = _mpf["chf"][2];
  toyf["chf"]["hb1"] = _mpf["chf"][2];
  toyf["chf"]["hm3"] = _mpf["chf"][3];
  toyf["chf"]["em3"] = _mpf["chf"][4];
  toyf["chf"]["pm3"] = _mpf["chf"][1];
  toyf["chf"]["tm3"] = _mpf["chf"][0]; 
  toyf["chf"]["tv2"] = _mpf["chf"][0]; 
  toyf["chf"]["tv3"] = _mpf["chf"][0]; 
  toyf["chf"]["tv3b"] = _mpf["chf"][0]; 
  toyf["chf"]["tv2b"] = _mpf["chf"][0]; 
  toyf["chf"]["tv2c"] = _mpf["chf"][0]; 
  toyf["chf"]["cmb"] = 0;
  //
  toyf["nhf"]["hp3"] = 0;
  toyf["nhf"]["hx3"] = _mpf["nhf"][2];
  toyf["nhf"]["hc1"] = _mpf["nhf"][2];
  toyf["nhf"]["hr3"] = _mpf["nhf"][2];
  toyf["nhf"]["hg3"] = _mpf["nhf"][2];
  toyf["nhf"]["hb3"] = _mpf["nhf"][2];
  toyf["nhf"]["hb1"] = _mpf["nhf"][2];
  toyf["nhf"]["hm3"] = _mpf["nhf"][3];
  toyf["nhf"]["em3"] = _mpf["nhf"][4];
  toyf["nhf"]["pm3"] = _mpf["nhf"][1];
  toyf["nhf"]["tm3"] = _mpf["nhf"][0]; 
  toyf["nhf"]["tv2"] = _mpf["nhf"][0]; 
  toyf["nhf"]["tv3"] = _mpf["nhf"][0]; 
  toyf["nhf"]["tv3b"] = _mpf["nhf"][0]; 
  toyf["nhf"]["tv2b"] = _mpf["nhf"][0]; 
  toyf["nhf"]["tv2c"] = _mpf["nhf"][0]; 
  toyf["nhf"]["cmb"] = 0;
  //
  toyf["nef"]["hp3"] = 0;
  toyf["nef"]["hx3"] = _mpf["nef"][2];
  toyf["nef"]["hc1"] = _mpf["nef"][2];
  toyf["nef"]["hr3"] = _mpf["nef"][2];
  toyf["nef"]["hg3"] = _mpf["nef"][2];
  toyf["nef"]["hb3"] = _mpf["nef"][2];
  toyf["nef"]["hb1"] = _mpf["nef"][2];
  toyf["nef"]["hm3"] = _mpf["nef"][3];
  toyf["nef"]["em3"] = _mpf["nef"][4];
  toyf["nef"]["pm3"] = _mpf["nef"][1];
  toyf["nef"]["tm3"] = _mpf["nef"][0]; 
  toyf["nef"]["tv2"] = _mpf["nef"][0]; 
  toyf["nef"]["tv3"] = _mpf["nef"][0]; 
  toyf["nef"]["tv3b"] = _mpf["nef"][0]; 
  toyf["nef"]["tv2b"] = _mpf["nef"][0]; 
  toyf["nef"]["tv2c"] = _mpf["nef"][0]; 
  toyf["nef"]["cmb"] = 0;

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
  name["hb3"] = "customHCALblue097";
  name["hb1"] = "customHCALblue1";
  name["hm3"] = "HadHCALm3";
  //name["em3"] = "HadECALm3"; // toyPF placeholder
  name["em3"] = "ECALm3"; // toyPF placeholder
  name["pm3"] = "Photonm3";  // toyPF placeholder
  name["tm3"] = "Trkm3";
  name["tv2"] = "Trkv2";
  name["tv2b"] = "Trkv2b";
  name["tv2c"] = "Trkv2c";
  name["tv3"] = "Trkv3";
  name["tv3b"] = "Trkv3b";
  //
  name["tv4"] = "Trkm3";
  name["tv402"] = "Trkv4a";
  name["tv410"] = "Trkv4b";
  name["tv404"] = "Trkv4d";
  name["tv430"] = "Trkv4e";
  name["tv405"] = "Trkv4f";
  name["tv416"] = "Trkv4g";
  name["tv420"] = "Trkv4h";
  //
  name["tv3"] = "Trkm3";
  name["tv3n1"] = "TrkNtrk1";
  name["tv300pn"] = "Trkv3c";
  name["tv301pn"] = "Trkv3d";
  name["tv310pn"] = "Trkv3b";
  //
  name["hhp3"] = "HadHCALp3";
  name["hhm3"] = "HadHCALm3";
  name["hhred103"] = "customHCALredNptcl";
  name["hhred100"] = "customHCALred1";
  name["hhred097"] = "customHCALred097";
  name["hhblue103"] = "customHCALgreenNptcl";
  name["hhblue100"] = "customHCALblue1";
  name["hhblue097"] = "customHCALblue097";
  //
  name["cmb"] = "Combined";
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
  label["hb3"] = "HCAL Blue -3%";
  label["hb1"] = "HCAL Blue 0%";
  label["hm3"] = "HCAL Had. -3%";
  label["em3"] = "ECAL Had. -3%";// (hybrid)";//(toyPF)";
  label["pm3"] = "Photons -3% (toyPF)";
  //label["tm3"] = "Tracking -1%";//patched from -3%";
  label["tm3"] = "Tracking -3%";
  label["tv2"] = "Tracking -1% to -3% (n>1)";
  label["tv3"] = "Tracking (-1%)^{n}";
  label["tv3b"] = "Tracking (-1%)#times(0.1%)^{n}";
  label["tv2b"] = "Tracking -1% to -3% (n>10)";
  label["tv2c"] = "Tracking -0.5% to -3% (n>10)";
  //
  label["tv4"] = "Tracking -3%";
  label["tv402"] = "Tracking -3% (n#geq2)";
  label["tv410"] = "Tracking -3% (n#geq10)";
  label["tv404"] = "Tracking -3% (n#geq4)";
  label["tv430"] = "Tracking -3% (n#geq30)";
  label["tv405"] = "Tracking -3% (n#geq5)";
  label["tv416"] = "Tracking -3% (n#geq16)";
  label["tv420"] = "Tracking -3% (n#geq20)";
  //
  label["tv3"] = "Tracking -3%";
  label["tv3n1"] = "Tracking -3% (n=1)";
  label["tv300pn"] = "Tracking -0.1%#times(n-1)";
  label["tv301pn"] = "Tracking -0.1%-0.1%#times(n-1)";
  label["tv310pn"] = "Tracking -1.0%-0.1%#times(n-1)";
  //
  label["hhp3"] = "HCAL Had. +3%";
  label["hhm3"] = "HCAL Had. -3%";
  label["hhred103"] = "HCAL Red +3%";
  label["hhred100"] = "HCAL Red 0%";
  label["hhred097"] = "HCAL Red -3%";
  label["hhblue103"] = "HCAL Blue +3% (\"Green\")";
  label["hhblue100"] = "HCAL Blue 0%";
  label["hhblue097"] = "HCAL Blue -3%";
  //
  label["cmb"] = "Combined";
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
  marker["hb3"] = kFullStar;
  marker["hb1"] = kFullStar;
  marker["hm3"] = kFullCircle;
  marker["em3"] = kOpenTriangleUp;
  marker["pm3"] = kFullDiamond;
  marker["tm3"] = kFullSquare;
  marker["tv2"] = kOpenSquare;
  marker["tv3"] = kOpenSquare;
  marker["tv3b"] = kOpenSquare;
  marker["tv2b"] = kOpenSquare;
  marker["tv2c"] = kOpenSquare;
  //
  marker["tv4"] = kFullSquare;
  marker["tv402"] = kFullSquare;
  marker["tv410"] = kFullSquare;
  marker["tv404"] = kOpenSquare;
  marker["tv430"] = kFullSquare;
  marker["tv405"] = kOpenDiamond;
  marker["tv416"] = kOpenSquare;
  marker["tv420"] = kOpenDiamond;
  //
  marker["tv3"] = kFullSquare;
  marker["tv3n1"] = kFullSquare;
  marker["tv300pn"] = kFullCircle;
  marker["tv301pn"] = kOpenCircle;
  marker["tv310pn"] = kFullDiamond;
  //
  marker["hhp3"] = kFullSquare;
  marker["hhm3"] = kFullSquare;
  marker["hhred103"] = kFullCircle;
  marker["hhred100"] = kOpenCircle;
  marker["hhred097"] = kFullCircle;
  marker["hhblue103"] = kFullDiamond;
  marker["hhblue100"] = kOpenDiamond;
  marker["hhblue097"] = kFullDiamond;
  //
  marker["cmb"] = kFullCircle;
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
  color["hb3"] = kBlue+1;
  color["hb1"] = kBlue;
  color["hm3"] = kBlue+1;
  color["em3"] = kBlue-9;
  color["pm3"] = kCyan+2;
  color["tm3"] = kBlack;//kGreen+2;
  color["tv2"] = kCyan+1;
  color["tv3"] = kMagenta+1;
  color["tv3b"] = kMagenta+3;
  color["tv2b"] = kCyan+2;
  color["tv2c"] = kCyan+4;
  //
  color["tv4"] = kRed+2;
  color["tv402"] = kBlue-6;//kRed-9;
  color["tv410"] = kRed-7;
  color["tv404"] = kRed-8;
  color["tv430"] = kRed-6;
  color["tv405"] = kBlue-9;
  color["tv416"] = kBlue-8;
  color["tv420"] = kBlue-7;
  //
  color["tv3"] = kRed+2;
  color["tv3n1"] = kRed-9;
  color["tv300pn"] = kRed-7;
  color["tv301pn"] = kRed-8;
  color["tv310pn"] = kRed-6;
  //
  color["hhp3"] = kRed+2;
  color["hhm3"] = kBlue+2;
  color["hhred103"] = kRed-9;
  color["hhred100"] = kRed;
  color["hhred097"] = kRed-8;
  color["hhblue103"] = kBlue-9;
  color["hhblue100"] = kBlue;
  color["hhblue097"] = kBlue-8;
  //
  color["cmb"] = kBlack;
  // observables
  color["Rjet"] = kBlack;
  color["chf"] = kRed;
  color["nhf"] = kGreen+2;
  color["nef"] = kBlue;
  
  // Load histograms
  map<string, TH1D*> hist;
  for (unsigned int i = 0; i != vars.size(); ++i) {

    string sv = vars[i];
    const char *cobs = (sysMode ? name[mode] : name[sv]);
    const char *csys = (sysMode ? name[sv] : name[mode]);
    string obs = cobs;
    string sys = csys;

    TH1D *hv = (TH1D*)f->Get(Form("%s_%s",cobs,csys));
    if (mode=="cmb") { // Combination mode

      assert(cmbs.size()!=0);
      hv = (TH1D*)f->Get(Form("%s_%s",cobs,name[cmbs[0].first]));
      assert(hv);
      // Clone to avoid changing original
      hv = (TH1D*)hv->Clone(Form("%s_%s_%s",mode.c_str(),cobs,csys));
      hv->Reset();

      for (unsigned int j = 0; j != cmbs.size(); ++j) {
	const char *csys2 = name[cmbs[j].first];
	TH1D *hv2 = (TH1D*)f->Get(Form("%s_%s",cobs,csys2));
	hv2 = (TH1D*)hv2->Clone(Form("%s_%s_%s_2",mode.c_str(),cobs,csys));
	assert(hv2);
	scaleError(hv2,0.1); // Patch FullSim uncertainty
	if (obs=="Rjet") toPercentage(hv2);
	hv->Add(hv2,cmbs[j].second);
      } // for j
    }
    else if (hv) {
      // Clone to avoid changing original
      hv = (TH1D*)hv->Clone(Form("%s_%s_%s",mode.c_str(),cobs,csys));

      //scaleError(hv,0.1); // Patch FullSim uncertainty
      if (obs=="Rjet") scaleError(hv,0.1); // Patch FullSim uncertainty
      //if (obs=="Rjet") toPercentage(hv);
      //else {
      //hv->Scale(100.); scaleError(hv,0.1);
      //}
    }
    else { // toyPF placeholders
      assert(false);
      //hv = (TH1D*)ftoy->Get(Form("h_%s_%s",cobs,csys));
      //if (!hv) hv = (TH1D*)ftoy->Get(Form("h%s_%s",cobs,csys));
      //if (!hv) cout << "Hist " << sv << " not found!" << endl << flush;
      //assert(hv);

      // Clone to avoid changing original
      //hv = (TH1D*)hv->Clone(Form("%s_%s_%s",mode.c_str(),cobs,csys));
    }
    
    if (obs=="Rjet") {
      if (mode!="cmb") toPercentage(hv);
    }
    else {
      hv->Scale(100.); scaleError(hv,0.1);
      if (obs=="nef") scaleError(hv,0.1);
    }

    // Set 15-20 GeV to zero
    string sv2 = (sysMode ? sv : mode);    
    if (TString(sv2.c_str()).Contains("tv4") ||
	TString(sv2.c_str()).Contains("tv3")) {
      if (sv2!="tv402") hv->SetBinContent(1,0.);
      bool isn1 = (sv2=="tv4" || sv2=="tv3n1" ||
		   sv2=="tv301pn" || sv2=="tv310pn");
      double k = 1;
      if (sv2=="tv301pn") k = 1./30.;
      if (sv2=="tv310pn") k = 1./3.;
      if (isn1 && obs=="Rjet") hv->SetBinContent(1,-1.4*k);
      if (isn1 && obs=="chf")  hv->SetBinContent(1,-1.0*k);
      if (isn1 && obs=="nhf")  hv->SetBinContent(1,+0.5*k);
      if (isn1 && obs=="gammaf")  hv->SetBinContent(1,+0.5*k);
      // Special treatment of tv402 to improve 15-100 GeV for global fit
      if (sv2=="tv402") {
	//if (obs=="Rjet") hv->SetBinContent(1,-0.3);
	//if (obs=="chf")  hv->SetBinContent(1,-0.2);
	hv->SetBinContent(1,0.);
	if (obs=="nhf")  { //hv->SetBinContent(1,+0.1);
	  int j600 = hv->FindBin(600.);
	  hv->SetBinContent(j600, 0); hv->SetBinError(j600, 0.);
	  // Ensure ECAL>HCAL>0 at pT<100 GeV, fix points with anomalous Rjet
	  hv->SetBinContent(3, hv->GetBinContent(3)+0.30+0.05); // 30
	  hv->SetBinError(3, 0.05); // 30
	  hv->SetBinContent(4, hv->GetBinContent(4)-0.15+0.05); // 42
	  hv->SetBinContent(5, hv->GetBinContent(5)-0.30); // 55
	  //hv->SetBinError(5, 0.1); // 55
	  hv->SetBinContent(7, hv->GetBinContent(7)-0.10); // 100
	  hv->SetBinContent(8, hv->GetBinContent(8)-0.10); // 150
	}
	if (obs=="gammaf")  { //hv->SetBinContent(1,+0.1);
	  // Ensure ECAL>HCAL>0 at pT<100 GeV, fix points with anomalous Rjet
	  hv->SetBinContent(3, hv->GetBinContent(3)-0.30+0.03); // 30
	  hv->SetBinError(3, 0.05); // 30
	  hv->SetBinContent(4, hv->GetBinContent(4)+0.15+0.00); // 42
	  hv->SetBinContent(5, hv->GetBinContent(5)+0.30); // 55
	  //hv->SetBinError(5, 0.05); // 55
	  hv->SetBinContent(7, hv->GetBinContent(7)+0.10); // 100
	  hv->SetBinContent(8, hv->GetBinContent(8)+0.10); // 150
	}
	//if (obs=="gammaf")  hv->SetBinContent(1,+0.1);
      }
      if (sv2!="tv402") hv->SetBinError(1,1e-4);
      if (sv2=="tv402") hv->SetBinError(1,0.05);
    }
    // Constraint 15-25 GeV for HCAL impacts
    if (TString(sv2.c_str()).Contains("hh")) {
      hv->SetBinContent(1,0.);
      // Boundaries
      if (sv2=="hhp3" && obs=="Rjet") hv->SetBinContent(1,+0.25);
      if (sv2=="hhp3" && obs=="chf")  hv->SetBinContent(1,-0.15);
      if (sv2=="hhp3" && obs=="nhf")  hv->SetBinContent(1,+0.20);
      if (sv2=="hhp3" && obs=="gammaf")  hv->SetBinContent(1,-0.05);
      if (sv2=="hhm3" && obs=="Rjet") hv->SetBinContent(1,-0.25);
      if (sv2=="hhm3" && obs=="chf")  hv->SetBinContent(1,+0.15);
      if (sv2=="hhm3" && obs=="nhf")  hv->SetBinContent(1,-0.20);
      if (sv2=="hhm3" && obs=="gammaf")  hv->SetBinContent(1,+0.05);
      // Reds
      if (TString(sv2.c_str()).Contains("hhred")) {
	if (obs=="Rjet") hv->SetBinContent(1,-0.125);
	if (obs=="chf")  hv->SetBinContent(1,+0.075);
	if (obs=="nhf")  hv->SetBinContent(1,-0.10);
	if (obs=="gammaf")  hv->SetBinContent(1,+0.025);
      }
      // Blues
      if (TString(sv2.c_str()).Contains("hhblue")) {
	if (obs=="Rjet") hv->SetBinContent(1,-0.125);
	if (obs=="chf")  hv->SetBinContent(1,+0.075);
	if (obs=="nhf")  hv->SetBinContent(1,-0.10);
	if (obs=="gammaf")  hv->SetBinContent(1,+0.025);
      }
      hv->SetBinError(1,1e-2);
      hv->SetBinContent(2,hv->GetBinContent(1));
      hv->SetBinError(2,2e-2);
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
  /*
  double maxy = (mode=="cmb" ? 2.499 : (mode=="Rjet" ? 3 : 3-1e-5));//2.5-1e-5);
  double miny = (mode=="cmb" ? -1.999 : (mode=="Rjet" ? -2 : -2+1e-5));
  */
  double maxy(3-1e-5), miny(-2+1e-5);
  if (mode=="tv3") { maxy = 10; miny = -10; }
  const char *title = (mode=="Rjet" ? "Response change (%)" :
		       sysMode ? "PF composition change (10^{-2})" :
		       "PF changes (% or 10^{-2})");
  TH1D *h = tdrHist(Form("h_%s",mode.c_str()), title,miny,maxy);
  lumi_13TeV = "FullSim";// (+toyPF placeholders)";
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
  TH1D *hsum(0); // for varMode
  for (unsigned int i = 0; i != vars.size(); ++i) {

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
      //if (sv=="tv2b") f1->SetRange(25.,3500.);
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
      if (sv=="Rjet"&&mode=="cmb") f1->SetRange(25.,3500.);
      f1->SetParameters(0,0,0,0);
      //hist[sv]->Fit(f1,"QRN");
      if (sv=="hc3") hist[sv]->Fit(f1,"QRNW");
      else hist[sv]->Fit(f1,"QRN");
      f1->SetLineColor(color[sv]);
      f1->SetLineWidth(2);
      f1->Draw("SAME");

      fits[mode][sv] = f1;
    }
    else if (funcs[sv][mode]!=0 || (sv=="Rjet" && func[mode]!=0)) {
      TF1 *f1 = new TF1(Form("f1_%s_%s_%s",mode.c_str(),cobs,csys),
			sv=="Rjet" ? func[mode] : funcs[sv][mode],15.,3500.);
      f1->SetParameters(0,0,0,0);
      if (mode=="tm3") f1->SetRange(20,2500);
      if (mode=="tv2c") f1->SetRange(15,2500);
      hist[sv]->Fit(f1,"QRN");
      f1->SetLineColor(color[sv]);
      f1->SetLineWidth(2);
      f1->Draw("SAME");

      if (!hsum) {
	hsum = (TH1D*)hist[sv]->Clone(Form("hsum_%s",mode.c_str()));
	hsum->Reset();
      }
      if (hsum && sv!="Rjet") {
	for (int i = 1; i != hsum->GetNbinsX()+1; ++i) {
	  hsum->SetBinContent(i, hsum->GetBinContent(i) + 
			      f1->Eval(hsum->GetBinCenter(i)));
	} // for i
      }

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

  // Check that cumulative sum of PF fits is small (<0.1%) to avoid biases
  if (hsum) {
    tdrDraw(hsum,"HIST",kNone,kBlack,kDotted,-1,kNone);
  }
  // Print functions for use in globalFitSettings.h
  if (hsum && doProduction) {

    ofstream fout(Form("pdf/fullSimShapes/txt/%s.txt",mode.c_str()));
    cout << "Copy-paste these into globalFitSettings.h:_gf_shapes array"<<endl;
    cout << "=========================================================="<<endl;

    TString md5;
    typedef map<string, map<string, TF1*> >::iterator IT;
    typedef map<string, TF1*>::iterator JT;
    for (IT it = fits.begin(); it != fits.end(); ++it) {
      for (JT jt = it->second.begin(); jt != it->second.end(); ++jt) {

	//cout << "Test fits[" << it->first << "][" << jt->first << "]" << endl;
	//cout << jt->second->GetExpFormula().Data() << endl;
	TF1 *f1 = jt->second;
	TString t = f1->GetExpFormula();
	for (int i = 0; i != f1->GetNpar(); ++i) {
	  t.ReplaceAll(Form("[p%d]",i),Form("%1.4g",f1->GetParameter(i)));
	} // for i in f1
	//cout << t.Data() << "(checksum " << t.MD5() << ")" << endl;

	cout << Form("    {\"%s\",\"%s\",\"%s\"},\n",
		     it->first.c_str(), jt->first.c_str(), t.Data());
	fout << Form("    {\"%s\",\"%s\",\"%s\"},\n",
		     it->first.c_str(), jt->first.c_str(), t.Data());
	md5 += t;
      } // for jt in fits.second
    } // for it in fits
    //cout << md5 << endl;
    cout << "    // \""<<mode<<"\" checksum: " << md5.MD5() << endl;
    fout << "    // \""<<mode<<"\" checksum: " << md5.MD5() << endl;
    cout << "=========================================================="<<endl;
    c1->SaveAs(Form("pdf/fullSimShapes/prod/fullSimShapes_%s_%s.pdf",
		    mode.c_str(), md5.MD5().Data()));
    fout.close();
    gSystem->Exec(Form("cp pdf/fullSimShapes/txt/%s.txt pdf/fullSimShapes/txt/%s_%s.txt",mode.c_str(),mode.c_str(),md5.MD5().Data()));
  } // doProduction


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

    for (unsigned int i = 0; i != vars.size(); ++i) {
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


// Compare composition and response flavor differences to split sources
void fullSimFlavor(string obs, string var) {

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  const char *co = obs.c_str();
  const char *cv = var.c_str();
  string obs2 = (obs=="nef" ? "gammaf" : obs);
  const char *co2 = obs2.c_str();

  //TFile *f = new TFile("rootfiles/FullSim_100k_variations_v13.root","READ");
  //TFile *f = new TFile("rootfiles/FullSim_100k_variations_v14.root","READ");
  TFile *f = new TFile("rootfiles/FullSim_100k_variations_v15.root","READ");
  assert(f && !f->IsZombie());

  TFile *f2 = new TFile("rootfiles/jecdataRun2Test.root","READ");
  assert(f2 && !f2->IsZombie());
  f2->cd("ratio");
  gDirectory->cd("eta00-13");
  TDirectory *d2 = gDirectory;

  curdir->cd();


  // 4a=402, 4d=404, 4f=405, 4b=410
  // trkEffv3b: trkEff = 0.99*0.999^(Ntrk-1)
  //const char *label = "N_{trk}#geq2"; // 4a
  //const char *label = "N_{trk}#geq5"; // 4f
  //const char *label = "N_{trk}#geq10"; // 4b
  //const char *label = "#epsilon(trk) -3%";
  string slabel;
  if (var=="Trkv3b") slabel = "0.99#times0.999^{N_{trk}-1}";
  if (var=="customHCALgreenNptcl") slabel = "HCAL \"Green\"";
  const char *label = slabel.c_str();
  //const char *label = "0.999^{N_{trk}-1}";//Trkv3b"; // 3c
  //const char *label = "0.99^{N_{trk}-1}";//Trkv3";
  double ymin(-2.0+1e-5), ymax(+1.0-1e-5); // chf
  //double ymin(-0.4), ymax(+2.6); // nhf
  //TH1D *ha = (TH1D*)f->Get("chf_Trkv4b"); assert(ha);
  //TH1D *hq = (TH1D*)f->Get("chfud_Trkv4b"); assert(hq);
  //TH1D *hg = (TH1D*)f->Get("chfg_Trkv4b"); assert(hg);
  //TH1D *ha = (TH1D*)f->Get("chf_Trkm3"); assert(ha);
  //TH1D *hq = (TH1D*)f->Get("chfud_Trkm3"); assert(hq);
  //TH1D *hg = (TH1D*)f->Get("chfg_Trkm3"); assert(hg);
  //TH1D *ha = (TH1D*)f->Get("chf_Trkv3b"); assert(ha);
  //TH1D *hq = (TH1D*)f->Get("chfud_Trkv3b"); assert(hq);
  //TH1D *hg = (TH1D*)f->Get("chfg_Trkv3b"); assert(hg);

  TH1D *ha = (TH1D*)f->Get(Form("%s_%s",co2,cv)); assert(ha);
  TH1D *hq = (TH1D*)f->Get(Form("%sud_%s",co2,cv)); assert(hq);
  TH1D *hg = (TH1D*)f->Get(Form("%sg_%s",co2,cv)); assert(hg);

  //TH1D *hc = (TH1D*)d2->Get("chf_cmb_ren_stat"); assert(hc);
  //TGraphErrors *gd = (TGraphErrors*)d2->Get("chf_pfjet_a30"); assert(gd);
  //TGraphErrors *gz = (TGraphErrors*)d2->Get("chf_zjet_a100"); assert(gz);
  //TGraphErrors *gp = (TGraphErrors*)d2->Get("chf_gamjet_ren"); assert(gp);
  //TH1D *hc = (TH1D*)d2->Get("nhf_cmb_ren_stat"); assert(hc);
  //TGraphErrors *gd = (TGraphErrors*)d2->Get("nhf_pfjet_a30"); assert(gd);
  //TGraphErrors *gz = (TGraphErrors*)d2->Get("nhf_zjet_a100"); assert(gz);
  //TGraphErrors *gp = (TGraphErrors*)d2->Get("nhf_gamjet_ren"); assert(gp);

  TH1D *hc(0);
  TGraphErrors *gd(0), *gz(0), *gp(0);
  if (obs=="Rjet") {
    hc = (TH1D*)d2->Get("hdm_cmb_mj"); assert(hc);
    TFile *f = new TFile("rootfiles/Zflavor_datamc_Eta13_Run2Test.root","READ");
    assert(f && !f->IsZombie());
    gd = (TGraphErrors*)f->Get("grg"); assert(gd);
    gz = (TGraphErrors*)f->Get("grq"); assert(gz);
    gp = (TGraphErrors*)f->Get("grc"); assert(gp);
  }
  else {
    hc = (TH1D*)d2->Get(Form("%s_cmb_ren_stat",co)); assert(hc);
    gd = (TGraphErrors*)d2->Get(Form("%s_pfjet_a30",co)); assert(gd);
    gz = (TGraphErrors*)d2->Get(Form("%s_zjet_a100",co)); assert(gz);
    gp = (TGraphErrors*)d2->Get(Form("%s_gamjet_ren",co)); assert(gp);
  }

  /*
  //const char *label = "HCAL -3%";
  const char *label = "HCAL Red +3%";
  double ymin(-0.4), ymax(+2.6);
  //TH1D *ha = (TH1D*)f->Get("nhf_HadHCALp3"); assert(ha);
  //TH1D *hq = (TH1D*)f->Get("nhfud_HadHCALp3"); assert(hq);
  //TH1D *hg = (TH1D*)f->Get("nhfg_HadHCALp3"); assert(hg);
  TH1D *ha = (TH1D*)f->Get("nhf_customHCALredNptcl"); assert(ha);
  TH1D *hq = (TH1D*)f->Get("nhfud_customHCALredNptcl"); assert(hq);
  TH1D *hg = (TH1D*)f->Get("nhfg_customHCALredNptcl"); assert(hg);

  TH1D *hc = (TH1D*)d2->Get("nhf_cmb_ren_stat"); assert(hc);
  TGraphErrors *gd = (TGraphErrors*)d2->Get("nhf_pfjet_a30"); assert(gd);
  TGraphErrors *gz = (TGraphErrors*)d2->Get("nhf_zjet_a100"); assert(gz);
  TGraphErrors *gp = (TGraphErrors*)d2->Get("nhf_gamjet_ren"); assert(gp);
  */

  TH1D *hgq = (TH1D*)ha->Clone("hgq"); hgq->Reset();
  hgq->Add(hg,hq,100,-100);

  TH1D *hr = (TH1D*)ha->Clone("hr"); hr->Reset();
  TH1D *hr2 = (TH1D*)ha->Clone("hr2"); hr2->Reset();
  TH1D *hd = (TH1D*)ha->Clone("hd"); hd->Reset();
  TH1D *hz = (TH1D*)ha->Clone("hz"); hz->Reset();
  TH1D *hp = (TH1D*)ha->Clone("hp"); hp->Reset();
  double k = (obs=="Rjet" ? 1 : 100.);
  for (int i = 0; i != gd->GetN(); ++i) {
    int j = hd->FindBin(gd->GetX()[i]);
    hd->SetBinContent(j, k*gd->GetY()[i]);
    hd->SetBinError(j, k*gd->GetEY()[i]);
  }
  for (int i = 0; i != gz->GetN(); ++i) {
    int j = hd->FindBin(gz->GetX()[i]);
    hz->SetBinContent(j, k*gz->GetY()[i]);
    hz->SetBinError(j, k*gz->GetEY()[i]);
  }
  for (int i = gp->GetN()-1; i != -1; --i) {
    if (gp->GetX()[i]<200.) {
      gp->RemovePoint(i);
    }
    else {
      int j = hd->FindBin(gp->GetX()[i]);
      hp->SetBinContent(j, k*gp->GetY()[i]);
      hp->SetBinError(j, k*gp->GetEY()[i]);
    }
  }
  hr->Add(hd,hz,1,-1);
  hr2->Add(hd,hp,1,-1);
  // Zero out non-overlapping range
  for (int i = 1; i != hr->GetNbinsX()+1; ++i) {
    if (hd->GetBinError(i)==0 || hz->GetBinError(i)==0 ||
	hr->GetBinCenter(i)<30.) {
      hr->SetBinContent(i, 0.);
      hr->SetBinError(i, 0.);
    }
    if (hd->GetBinError(i)==0 || hp->GetBinError(i)==0 ||
	hr->GetBinCenter(i)<30.) {
      hr2->SetBinContent(i, 0.);
      hr2->SetBinError(i, 0.);
    }
  }

  lumi_13TeV = "Run 2 Legacy, 138 fb^{-1}";
  extraText = "Private";
  TH1D *h1u = tdrHist("h1u","Data - MC (%)",ymin,ymax);//-2+1e-5,+1-1e-5);
  TH1D *h1d = tdrHist("h1d","G - UD (%)",-0.5,+0.5);
  if (obs=="Rjet") {
    h1d->GetYaxis()->SetRangeUser(-1.6,+0.8);
  }
  if (var=="customHCALgreenNptcl" && obs=="nhf") {
    h1u->GetYaxis()->SetRangeUser(-0.4+1e-5,+2.6-1e-5);
  }
  TCanvas *c1 = tdrDiCanvas("c1",h1u,h1d,4,11);

  c1->cd(1);
  gPad->SetLogx();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(15,0,3500,0);

  if (obs=="Rjet") {
    hq->Add(ha,-1); hq->Divide(ha); hq->Scale(100.);
    hg->Add(ha,-1); hg->Divide(ha); hg->Scale(100.);
  }
  else {
    ha->Scale(100.); hq->Scale(100.); hg->Scale(100.);
  }
  tdrDraw(hg,"HIST",kNone,kRed,kSolid,-1,kNone);
  tdrDraw(hq,"HIST",kNone,kBlue,kSolid,-1,kNone);
  if (obs!="Rjet") tdrDraw(ha,"HIST",kNone,kBlack,kSolid,-1,kNone);

  if (obs!="Rjet") {
    hc->Scale(100.); scaleGraph(gd,100); scaleGraph(gz,100); scaleGraph(gp,100);
  }
  tdrDraw(hc,"Pz",kFullDiamond,kBlack);
  tdrDraw(gd,"Pz",kFullSquare,kRed);
  tdrDraw(gz,"Pz",kFullCircle,kBlue);
  tdrDraw(gp,"Pz",kOpenCircle,kBlue);

  if (obs=="Rjet") {
    TLegend *leg = tdrLeg(0.60,0.90-4*0.05,0.80,0.90);
    leg->AddEntry(gz,"Z+jet (q-jet)","PLE");
    leg->AddEntry(gd,"Z+jet (g-jet)","PLE");
    leg->AddEntry(hq,"Quark MC","L");
    leg->AddEntry(hg,"Gluon MC","L");
  }
  else {
    TLegend *leg = tdrLeg(0.60,0.90-7*0.05,0.80,0.90);
    leg->AddEntry(gz,"Z+jet (q-rich)","PLE");
    leg->AddEntry(gp,"#gamma+jet","PLE");
    leg->AddEntry(hc,Form("Combined %s",co),"PLE");
    //leg->AddEntry(hc,"Combined NHF","PLE");
    leg->AddEntry(gd,"Dijet (g-rich)","PLE");
    leg->AddEntry(hq,"Quark MC","L");
    leg->AddEntry(ha,"Combined MC","L");
    leg->AddEntry(hg,"Gluon MC","L");
  }

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  tex->DrawLatex(0.34,0.85,label);

  c1->cd(2);
  gPad->SetLogx();

  l->DrawLine(15,0,3500,0);

  tdrDraw(hgq,"HIST",kNone,kRed,kSolid,-1,kNone);
  tdrDraw(hr,"Pz",kFullCircle,kRed,kSolid,-1,kNone);
  tdrDraw(hr2,"Pz",kOpenCircle,kRed,kSolid,-1,kNone);

  c1->SaveAs(Form("pdf/fullSimShapes/fullSimShapes_fullSimFlavor_%s_%s.pdf",
		  co,cv));
} // fullSimFlavor


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
} // toPercentage

void scaleGraph(TGraph *g, double scale) {
  for (int i = 0; i != g->GetN(); ++i) {
    g->SetPoint(i, g->GetX()[i], g->GetY()[i]*scale);
    if (g->InheritsFrom("TGraphErrors")) {
      ((TGraphErrors*)g)->SetPointError(i, g->GetEX()[i], g->GetEY()[i]*scale);
    }
  } // for i
} // scaleGraph
