#include "TString.h"
#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"

#include "../tdrstyle_mod15.C"

#include <iostream>
#include <fstream>
#include <map>

using namespace std;

const bool debug = false;

void createL2L3ResTextFiles(string set="BCDEF SimpleL1");

TLegend *_leg(0);
void createL2L3ResTextFile() {
  /*
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
  */

  setTDRStyle();

  double ptmin = 15;
  double ptmax = 4500;
  TH1D *h = tdrHist("h","Absolute response at |#eta| < 1.3",
		    0.97,1.01,"p_{T} (GeV)",ptmin,ptmax);
  lumi_13TeV = "UL2017";
  TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);
  c1->SetLeftMargin(0.17);
  c1->SetRightMargin(0.03);
  h->SetTitleOffset(1.5,"Y");
  gPad->SetLogx();

  _leg = tdrLeg(0.45,0.60,0.75,0.90);

  createL2L3ResTextFiles("BCDEF");
  createL2L3ResTextFiles("B");
  createL2L3ResTextFiles("C");
  createL2L3ResTextFiles("D");
  createL2L3ResTextFiles("E");
  createL2L3ResTextFiles("F");
  
  c1->SaveAs("pdf/createL2L3ResTextFile.pdf");
}

void createL2L3ResTextFiles(string set) {

  if (debug) cout << "Warning: sscanf only works correctly when code is compiled (.C+)\n";

  ////////////////////////////////////////////////////
  // Copy these values by hand from                 //
  // textFiles/globalFitL3Res.txt produced by       //
  // globalFitL3Res.C                               //
  //////////e/////////////////////////////////////////
  /*
  const int np = 3;
  double p[np] = {0,0,0};

  // V2M3
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
  */
  // V2M4
  /*
  if (set=="BCDEF SimpleL1") { p[0] = 0.9811, p[1] = 1.0868, p[2] = -0.1813; }
  if (set=="B SimpleL1") { p[0] = 0.9811, p[1] = 1.0868, p[2] = -0.1813; }
  if (set=="C SimpleL1") { p[0] = 0.9811, p[1] = 1.0868, p[2] = -0.1813; }
  if (set=="D SimpleL1") { p[0] = 0.9811, p[1] = 1.0868, p[2] = -0.1813; }
  if (set=="E SimpleL1") { p[0] = 0.9811, p[1] = 1.0868, p[2] = -0.1813; }
  if (set=="F SimpleL1") { p[0] = 0.9811, p[1] = 1.0868, p[2] = -0.1813; }
  // Complex cloned from Simple above for now
  if (set=="BCDEF ComplexL1") { p[0] = 0.9811, p[1] = 1.0868, p[2] = -0.1813; }
  if (set=="B ComplexL1") { p[0] = 0.9811, p[1] = 1.0868, p[2] = -0.1813; }
  if (set=="C ComplexL1") { p[0] = 0.9811, p[1] = 1.0868, p[2] = -0.1813; }
  if (set=="D ComplexL1") { p[0] = 0.9811, p[1] = 1.0868, p[2] = -0.1813; }
  if (set=="E ComplexL1") { p[0] = 0.9811, p[1] = 1.0868, p[2] = -0.1813; }
  if (set=="F ComplexL1") { p[0] = 0.9811, p[1] = 1.0868, p[2] = -0.1813; }
  assert(!(p[0]==0 && p[1]==0 && p[2]==0));
  */

  // For V2M5, simplify complex sum into an effective formula
  // Need good starting values and/or a few iterations to converge
  // Fit done to hjesfit from each IOV
  TDirectory *curdir = gDirectory;
  TFile *f = new TFile(Form("rootfiles/jecdata%s.root",set.c_str()),"READ");
  assert(f && !f->IsZombie());
  TH1D *h = (TH1D*)f->Get("ratio/eta00-13/sys/hjesfit"); assert(h);
  curdir->cd();

  TF1 *f1 = new TF1(Form("f1_%s",set.c_str()),"[0]+[1]/x+[2]*log(x)/x+[3]*(pow(x/[4],[5])-1)/(pow(x/[4],[5])+1)+[6]*pow(x,-0.3051)",15,4500);
  f1->SetParameters(0.98, 0.1,0.01, 0.01,500.,1.3, 0.001);

  map<string,int> color;
  color["BCDEF"] = kYellow+2;
  color["B"] = kBlue;
  color["C"] = kGreen+2;
  color["D"] = kOrange+2;
  color["E"] = kMagenta+2;
  color["F"] = kRed;

  h->Fit(f1,"QRN");
  h->Fit(f1,"QRNM");
  h->Fit(f1,"QRNM");
  tdrDraw(h,"LE3",kNone,color[set],kSolid,-1,1001,color[set]-9);
  h->SetFillColorAlpha(color[set]-9,0.7);
  f1->SetLineColor(color[set]);
  f1->Draw("SAME");

  _leg->AddEntry(h,set.c_str(),"FL");
  
  const int np = 7;
  double p[np];
  for (int i = 0; i != np; ++i) {
    p[i] = f1->GetParameter(i);
  } // for i


  /////////////////////////////////////////////////////////////
  // Generate input and output file names semi-automatically  //
  //////////////////////////////////////////////////////////////
  //char run[512], l1[512];
  //sscanf(set.c_str(),"%s %s",run,l1);
  //cout << "Processing " << run << "_" << l1 << endl;
  const char *run = set.c_str();
  const char *l1 = "SimpleL1";
  cout << "Processing " << set << endl;

  ifstream fin(Form("textFiles/UL17V2-L2Res+JERSF/%s/Run%s/Summer19UL17_V1_%s_MPF_LOGLIN_L2Residual_pythia8_AK4PFchs.txt",l1,run,l1));
  //ofstream fout(Form("textFiles/UL17V2-L2L3Res+JERSF/Summer19UL17_Run%s_V2M2_%s_DATA_L2L3Residual_AK4PFchs.txt",run,l1));
  //ofstream fout(Form("textFiles/UL17V2MX-L2L3Res+JERSF/Summer19UL17_Run%s_V2M3_%s_DATA_L2L3Residual_AK4PFchs.txt",run,l1));
  //ofstream fout(Form("textFiles/UL17V2MX-L2L3Res+JERSF/Summer19UL17_Run%s_V2M4_%s_DATA_L2L3Residual_AK4PFchs.txt",run,l1));
  ofstream fout(Form("textFiles/UL17V2MX-L2L3Res+JERSF/Summer19UL17_Run%s_V2M5_%s_DATA_L2L3Residual_AK4PFchs.txt",run,l1));

  // 2018: https://twiki.cern.ch/twiki/bin/view/CMSPublic/PixelOfflinePlotsOctober2018

  // Print out response function to copy below by hand
  /*
  // V2M4
  TF1 *fx = new TF1("fx","[p0]+[p1]*pow(x/[p2],[p3])/(1+pow(x/[p2],[p3]))*(1-pow(x/[p2],-[p3]))",15,4500);
  fx->SetParameters(1.184, 1.428, 1402, 1.225); // toyPF
  if (set=="BCDEF SimpleL1")
    cout << Form("1./([6]+[7]*(%1.5f+%1.5f*TMath::Power(x/%1.0f.,%1.3f)"
		 "/(1+TMath::Power(x/%1.0f.,%1.3f))"
		 "*(1-TMath::Power(x/%1.0f.,-%1.3f)))",
		 0.01*fx->GetParameter(0), 0.01*fx->GetParameter(1),
		 fx->GetParameter(2), fx->GetParameter(3),
		 fx->GetParameter(2), fx->GetParameter(3),
		 fx->GetParameter(2), fx->GetParameter(3));

  // Print out L1 function to copy below by hand
  TF1 *fl1 = new TF1("fl1","1-([0]+[1]*log(x)+[2]*pow(log(x),2))/x",10,3500);
  fl1->SetParameters(0.350077, 0.553560, -0.0527681);
  if (set=="BCDEF SimpleL1")
    cout << Form("+[8]*(%1.6f-(%1.6f+TMath::Log(x)*(%1.6f%+1.7f*TMath::Log(x)))"
		 "/x))",
		 1-fl1->Eval(208.),
		 fl1->GetParameter(0),
		 fl1->GetParameter(1),
		 fl1->GetParameter(2)) << endl;
  */

  /////////////////////////////////
  /////////////////////////////////
  // Rest is automatic
  ////////////////////////////////
  ////////////////////////////////

  string header;
  getline(fin, header);
  if (debug) cout << header << endl;
  //header = "{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[8]*1.72396*(1./208.-1./x))) Correction L2Relative}";
  //header = "{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[8]*(0.00866-(0.35008+TMath::Log(x)*(0.55356-0.05277*TMath::Log(x)))/x))) Correction L2Relative}";
  //header = "{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]*(0.01184+0.01428*TMath::Power(x/1402.,1.225)/(1+TMath::Power(x/1402.,1.225))*(1-TMath::Power(x/1402.,-1.225)))+[8]*(0.008661-(0.350077+TMath::Log(x)*(0.553560-0.0527681*TMath::Log(x)))/x))) Correction L2Relative}"; // V2M4
  header = "{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]/x+[8]*log(x)/x+[9]*(pow(x/[10],[11])-1)/(pow(x/[10],[11])+1)+[12]*pow(x,-0.3051))) Correction L2Relative}"; // V2M5
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
    /*
    // V2M4
    if (cnt<cntmax && debug)
      cout << Form("  %9.6f %9.6f   %d   %d %d   %d   %d   %8.6f %8.6f"
		   "   %8.6f %8.6f   %5.3f %5.4f %5.3f",
		   etamin, etamax, npar, xmin, xmax, ptmin0, ptmax1,
		   p2, p3, p4, p5,  p[0], p[1], p[2]) << endl;
      fout << Form("  %9.6f %9.6f   %d   %d %d   %d   %d   %8.6f %8.6f"
		   "   %8.6f %8.6f   %5.4f %5.4f %5.4f",
		   etamin, etamax, npar, xmin, xmax, ptmin0, ptmax1,
		   p2, p3, p4, p5,  p[0], p[1], p[2]) << endl;
    */
    // V2M5
    int nparnew = 15;
    if (cnt<cntmax && debug)
      cout << Form("  %9.6f %9.6f   %d   %d %d   %d   %d   %8.6f %8.6f"
		   "   %8.6f %8.6f   "
		   "%5.4f %5.4f %5.5f %5.5f %5.1f %5.4f %5.5f",
		   etamin, etamax, nparnew, xmin, xmax, ptmin0, ptmax1,
		   p2, p3, p4, p5,
		   p[0], p[1], p[2], p[3], p[4], p[5], p[6]) << endl;
      fout << Form("  %9.6f %9.6f   %d   %d %d   %d   %d   %8.6f %8.6f"
		   "   %8.6f %8.6f   "
		   "%5.4f %5.4f %5.5f %5.5f %5.1f %5.4f %5.5f",
		   etamin, etamax, nparnew, xmin, xmax, ptmin0, ptmax1,
		   p2, p3, p4, p5,
		   p[0], p[1], p[2], p[3], p[4], p[5], p[6]) << endl;
    ++cnt;
  }

} // creataL2L3ResTextFile
