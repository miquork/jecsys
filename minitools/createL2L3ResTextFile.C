#include "TString.h"

#include <iostream>
#include <fstream>
#include <map>

using namespace std;

const bool debug = false;

void createL2L3ResTextFiles(string set="BCDEF SimpleL1");

void createL2L3ResTextFile() {
  createL2L3ResTextFiles("BCDEF SimpleL1");
  createL2L3ResTextFiles("B SimpleL1");
  createL2L3ResTextFiles("C SimpleL1");
  createL2L3ResTextFiles("D SimpleL1");
  createL2L3ResTextFiles("E SimpleL1");
  createL2L3ResTextFiles("F SimpleL1");
  //
  createL2L3ResTextFiles("BCDEF ComplexL1");
  createL2L3ResTextFiles("B ComplexL1");
  createL2L3ResTextFiles("C ComplexL1");
  createL2L3ResTextFiles("D ComplexL1");
  createL2L3ResTextFiles("E ComplexL1");
  createL2L3ResTextFiles("F ComplexL1");
}

void createL2L3ResTextFiles(string set) {

  if (debug) cout << "Warning: sscanf only works correctly when code is compiled (.C+)\n";

  ////////////////////////////////////////////////////
  // Copy these values by hand from                 //
  // textFiles/globalFitL3Res.txt produced by       //
  // globalFitL3Res.C                               //
  //////////e/////////////////////////////////////////
  const int np = 3;
  double p[np] = {0,0,0};
  if (set=="BCDEF SimpleL1") { p[0] = 0.9809, p[1] = 0.04856, p[2] = -0.586; }
  if (set=="B SimpleL1") { p[0] = 0.9870, p[1] = 0.08822, p[2] = -0.586; }
  if (set=="C SimpleL1") { p[0] = 0.9851, p[1] = 0.04675, p[2] = -0.586; }
  if (set=="D SimpleL1") { p[0] = 0.9836, p[1] = 0.04383, p[2] = -0.586; }
  if (set=="E SimpleL1") { p[0] = 0.9811, p[1] = 0.04080, p[2] = -0.586; }
  if (set=="F SimpleL1") { p[0] = 0.9739, p[1] = 0.04800, p[2] = -0.586; }
  // Complex cloned from Simple above for now
  if (set=="BCDEF ComplexL1") { p[0] = 0.9809, p[1] = 0.04856, p[2] = -0.586; }
  if (set=="B ComplexL1") { p[0] = 0.9870, p[1] = 0.08822, p[2] = -0.586; }
  if (set=="C ComplexL1") { p[0] = 0.9851, p[1] = 0.04675, p[2] = -0.586; }
  if (set=="D ComplexL1") { p[0] = 0.9836, p[1] = 0.04383, p[2] = -0.586; }
  if (set=="E ComplexL1") { p[0] = 0.9811, p[1] = 0.04080, p[2] = -0.586; }
  if (set=="F ComplexL1") { p[0] = 0.9739, p[1] = 0.04800, p[2] = -0.586; }
  assert(!(p[0]==0 && p[1]==0 && p[2]==0));
  
  /////////////////////////////////////////////////////////////
  // Generate input and output file names semi-automatically  //
  //////////////////////////////////////////////////////////////
  char run[512], l1[512];
  sscanf(set.c_str(),"%s %s",run,l1);
  cout << "Processing " << run << "_" << l1 << endl;
  
  ifstream fin(Form("textFiles/UL17V2-L2Res+JERSF/%s/Run%s/Summer19UL17_V1_%s_MPF_LOGLIN_L2Residual_pythia8_AK4PFchs.txt",l1,run,l1));
  ofstream fout(Form("textFiles/UL17V2-L2L3Res+JERSF/Summer19UL17_Run%s_V2M2_%s_DATA_L2L3Residual_AK4PFchs.txt",run,l1));


  /////////////////////////////////
  /////////////////////////////////
  // Rest is automatic
  ////////////////////////////////
  ////////////////////////////////

  string header;
  getline(fin, header);
  if (debug) cout << header << endl;
  header = "{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[8]*1.72396*(1./208.-1./x))) Correction L2Relative}";
  if (debug) cout << header << endl;
  fout << header << endl;
  
  string line;
  double etamin, etamax;
  int npar, xmin, xmax, ptmin0, ptmax1;
  double p2, p3, p4, p5;
  double p6, p7, p8;
  int cnt(0); int cntmax(0);
  while (getline(fin,line)) {
    if (cnt<cntmax && debug) cout << line << endl;
    assert(sscanf(line.c_str(),"%lf %lf  %d  %d %d  %d %d  %lf %lf %lf %lf"
		  "  %lf %lf %lf",
		  &etamin, &etamax, &npar, &xmin, &xmax, &ptmin0, &ptmax1,
		  &p2, &p3, &p4, &p5,  &p6, &p7, &p8)==14);
    if (cnt<cntmax && debug)
      cout << Form("  %9.6f %9.6f   %d   %d %d   %d   %d   %8.6f %8.6f"
		   "   %8.6f %8.6f   %5.3f %5.3f %5.3f",
		   etamin, etamax, npar, xmin, xmax, ptmin0, ptmax1,
		   p2, p3, p4, p5,  p[0], p[1], p[2]) << endl;
      fout << Form("  %9.6f %9.6f   %d   %d %d   %d   %d   %8.6f %8.6f"
		   "   %8.6f %8.6f   %5.3f %5.3f %5.3f",
		   etamin, etamax, npar, xmin, xmax, ptmin0, ptmax1,
		   p2, p3, p4, p5,  p[0], p[1], p[2]) << endl;
    ++cnt;
  }

} // creataL2L3ResTextFile
