// Purpose: Create L2L3Res text file with simple parameterization
//          Takes as input previous L2Res and complex 9p global JES fit
//          Outputs same L2Res and "simple" 7p fit
//          ("simple" as in removing main parameter degeneracies)
//          For merging IOVs together at text file level, use
//          minitools/mergeL2L3ResTextFiles.C
#include "TString.h"
#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TLine.h"

#include "../tdrstyle_mod15.C"

#include <iostream>
#include <fstream>
#include <map>

using namespace std;

const bool debug = true;

void createL2L3ResTextFiles(string set="2016GH");

TLegend *_leg(0);
void createL2L3ResTextFile() {

  setTDRStyle();

  double ptmin = 15;
  double ptmax = 4500;
  //TH1D *h = tdrHist("h","Absolute response at |#eta| < 1.3",
  //		    0.97,1.01,"p_{T} (GeV)",ptmin,ptmax);
  //TH1D *h = tdrHist("h","Absolute response at |#eta| < 1.3",
  //		    0.982,1.022,"p_{T} (GeV)",ptmin,ptmax);
  //TH1D *h = tdrHist("h","Absolute response at |#eta| < 1.3",
  //		    0.97,1.02,"p_{T} (GeV)",ptmin,ptmax);
  TH1D *h = tdrHist("h","Absolute response at |#eta| < 1.3",
		    //0.960,1.025,"p_{T} (GeV)",ptmin,ptmax); // V5M2
		    //0.960,1.025,"p_{T} (GeV)",ptmin,ptmax); // V5M3
		    0.94+1e-4,1.06-1e-4,"p_{T} (GeV)",ptmin,ptmax); // low PU
  		    //0.970,1.035,"p_{T} (GeV)",ptmin,ptmax);
  //lumi_13TeV = "UL2017";
  //lumi_13TeV = "UL2016GH";
  //lumi_13TeV = "UL2016, 36.5 fb^{-1}";
  lumi_13TeV = "Run2, 137.9 fb^{-1}";
  TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);
  c1->SetLeftMargin(0.17);
  c1->SetRightMargin(0.03);
  h->SetTitleOffset(1.5,"Y");
  gPad->SetLogx();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(ptmin,1,ptmax,1);

  _leg = tdrLeg(0.45,0.90,0.75,0.90);

  /*
  createL2L3ResTextFiles("2018ABCD");
  createL2L3ResTextFiles("2018A");
  createL2L3ResTextFiles("2018B");
  createL2L3ResTextFiles("2018C");
  createL2L3ResTextFiles("2018D");

  c1->Update();
  c1->SaveAs("pdf/createL2L3ResTextFile_UL18JECV3.pdf");
  */

  /*
  createL2L3ResTextFiles("Run2Test");
  createL2L3ResTextFiles("2016BCDEF");
  createL2L3ResTextFiles("2016GH");
  createL2L3ResTextFiles("2017BCDEF");
  createL2L3ResTextFiles("2017H"); // low PU
  createL2L3ResTextFiles("2018ABCD");

  c1->Update();
  c1->SaveAs("pdf/createL2L3ResTextFile_Run2_JECV7.pdf");
  */

  createL2L3ResTextFiles("2017BCDEF");
  createL2L3ResTextFiles("2017B");
  createL2L3ResTextFiles("2017C");
  createL2L3ResTextFiles("2017D");
  createL2L3ResTextFiles("2017E");
  createL2L3ResTextFiles("2017F");
  createL2L3ResTextFiles("2017H");

  c1->Update();
  c1->SaveAs("pdf/createL2L3ResTextFile_2017H_JECV7.pdf");

  /*
  createL2L3ResTextFiles("2016BCDEF");
  createL2L3ResTextFiles("2016BCD");
  createL2L3ResTextFiles("2016EF");
  createL2L3ResTextFiles("2016GH");
  //createL2L3ResTextFiles("2017BCDEF");
  //createL2L3ResTextFiles("2018ABCD");

  c1->Update();
  //c1->SaveAs("pdf/createL2L3ResTextFile_UL16GHJECV2.pdf");
  //c1->SaveAs("pdf/createL2L3ResTextFile_UL16_BtoF_JECV3_GH_JECV2.pdf");
  //c1->SaveAs("pdf/createL2L3ResTextFile_UL16_BtoH_JECV2.pdf");
  //c1->SaveAs("pdf/createL2L3ResTextFile_UL16_BtoH_JECV3M2.pdf");
  //c1->SaveAs("pdf/createL2L3ResTextFile_Run2_JECV3M2.pdf");
  //c1->SaveAs("pdf/createL2L3ResTextFile_Run2_UL16JECV5M1.pdf");
  //c1->SaveAs("pdf/createL2L3ResTextFile_UL16_BtoH_JECV5M2.pdf");
  //c1->SaveAs("pdf/createL2L3ResTextFile_UL16_BtoH_JECV5M3.pdf");
  c1->SaveAs("pdf/createL2L3ResTextFile_UL16_BtoH_JECV5M4.pdf");
  */
}

void createL2L3ResTextFiles(string set) {

  //if (debug) 
  cout << "Warning: sscanf only works correctly when code is compiled (.C+)\n";

  cout << "Processing " << set << endl << flush;

  // For V2M5, simplify complex sum into an effective formula
  // Need good starting values and/or a few iterations to converge
  // Fit done to hjesfit from each IOV
  TDirectory *curdir = gDirectory;
  TFile *f = new TFile(Form("rootfiles/jecdata%s.root",set.c_str()),"READ");

  assert(f && !f->IsZombie());
  //TH1D *h = (TH1D*)f->Get("ratio/eta00-13/sys/hjesfit"); assert(h);
  TH1D *h(0);
  h = (TH1D*)f->Get("ratio/eta00-13/sys/hjesfit2");
  if (!h) h = (TH1D*)f->Get("ratio/eta00-13/sys/hjesfit");
  assert(h);
  curdir->cd();

  TF1 *f1 = new TF1(Form("f1_%s",set.c_str()),"[0]+[1]/x+[2]*log(x)/x+[3]*(pow(x/[4],[5])-1)/(pow(x/[4],[5])+1)+[6]*pow(x,-0.3051)",15,4500);
  f1->SetParameters(0.98, 0.1,0.01, 0.01,500.,1.3, 0.001);
  if (set=="2017H") 
    f1->SetParameters(0.99, 1.5,0.01, 0.01,1000.,1.3, 0.001);

  map<string,int> color;
  color["2017BCDEF"] = kGreen+2;
  color["BCDEF"] = kYellow+2;
  color["B"] = kBlue;
  color["C"] = kGreen+2;
  color["D"] = kOrange+2;
  color["E"] = kMagenta+2;
  color["F"] = kRed;
  color["2017H"] = kCyan+2; // low PU
  //
  //color["2018ABCD"] = kYellow+2;
  color["2018ABCD"] = kBlue;
  color["2018A"] = kCyan+2;//kBlue;
  color["2018B"] = kGreen+2;
  color["2018C"] = kOrange+2;
  color["2018D"] = kRed;
  //
  color["2016BCDEF"] = kYellow+2;//kCyan+2;
  color["2016BCD"] = kBlue;
  color["2016EF"] = kGreen+2;
  color["2016GH"] = kRed;
  // 
  color["Run2Test"] = kMagenta+2;

  h->Fit(f1,"QRN");
  h->Fit(f1,"QRNM");
  h->Fit(f1,"QRNM");
  tdrDraw(h,"LE3",kNone,color[set],kSolid,-1,1001,color[set]-9);
  h->SetFillColorAlpha(color[set]-9,0.7);
  f1->SetLineColor(color[set]);
  f1->Draw("SAME");

  _leg->SetY1(_leg->GetY1()-0.05);
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
  //cout << "Processing " << set << endl;

  bool isUL18 = (set=="2018ABCD" || set=="2018A" || set=="2018B" || 
		 set=="2018C" || set=="2018D");

  bool isUL17 = (set=="2017BCDEF" || set=="2017B" || set=="2017C" || 
		 set=="2017D" || set=="2017E" || set=="2017F");
  bool isLowPU = (set=="2017H");

  bool isUL16 = (set=="2016BCDEFGH" || set=="2016BCDEF" || set=="2016BCD" ||
		 set=="2016EF" || set=="2016GH");
  bool isRun2 = (set=="Run2Test");

  string sin, sout;
  if (isUL18) {
    map<string,const char*> mera;
    mera["2018ABCD"] = "ABCD";
    mera["2018A"] = "A";
    mera["2018B"] = "B";
    mera["2018C"] = "C";
    mera["2018D"] = "D";
    //sin = Form("CondFormats/JetMETObjects/data/Summer19UL18_Run%s_V3_DATA_L2Residual_AK4PFchs.txt",set=="2018ABCD" ? "C" : mera[set]);
    //sin = Form("../JERCProtoLab/Summer19UL18/L2Residual/V2/Run%s/Summer19UL18_V2_MPF_LOGLIN_L2Residual_pythia8_AK4PFchs.txt",mera[set]); // official for V5
    sin = Form("../JECDatabase/textFiles/Summer19UL18_Run%s_V5_DATA/Summer19UL18_Run%s_V5_DATA_L2Residual_AK4PFchs.txt",mera[set],mera[set]);
    if (set=="2018ABCD") // minitools/mk_mergeL2L3ResTextFiles.C:("2018")
      sin = "textFiles/mergeL2ResTextFiles_2018.txt";
    //
    //sout = Form("textFiles/UL18V3-L2L3Res/Summer19UL18_Run%s_V3M1_DATA_L2L3Residual_AK4PFchs.txt",mera[set]);
    //sout = Form("../JERCProtoLab/Summer19UL18/global_fit/V3A2M1J2/Summer19UL18_Run%s_V3A2M1J2_DATA_L2L3Residual_AK4PFchs.txt",mera[set]);
    //sout = Form("../JERCProtoLab/Summer19UL16/global_fit/V3M2/Summer19UL18_Run%s_V3M2_DATA_L2L3Residual_AK4PFchs.txt",mera[set]);
    //sout = Form("../JERCProtoLab/Summer19UL16/global_fit/V5M1/Summer19UL18_Run%s_V5M1_DATA_L2L3Residual_AK4PFchs.txt",mera[set]);
    //sout = Form("../JERCProtoLab/Summer19UL16/global_fit/V5M2/Summer19UL18_Run%s_V5M2_DATA_L2L3Residual_AK4PFchs.txt",mera[set]);
    //sout = Form("../JERCProtoLab/Summer19UL16/global_fit/V5M4/Summer19UL18_Run%s_V5M4_DATA_L2L3Residual_AK4PFchs.txt",mera[set]);
    sout = Form("../JERCProtoLab/Summer20UL/global_fit/Summer19UL18_Run%s_VX_DATA_L2L3Residual_AK4PFchs.txt",mera[set]);
  }
  else if (isUL16) {
    map<string,const char*> mera;
    mera["2016BCDEFGH"] = "BCDEFGH";
    mera["2016BCDEF"] = "BCDEF";
    mera["2016BCD"] = "BCD";
    mera["2016EF"] = "EF";
    mera["2016GH"] = "FGH";
    if (set=="2016GH") {
      //sin = Form("../JECDatabase/textFiles/Summer19UL16_RunFGH_V2_DATA/Summer19UL16_Run%s_V2_DATA_L2Residual_AK4PFchs.txt",mera[set]);
      sin = Form("../JECDatabase/textFiles/Summer19UL16_RunFGH_V5_DATA/Summer19UL16_Run%s_V5_DATA_L2Residual_AK4PFchs.txt",mera[set]);
      //sout = Form("../JERCProtoLab/Summer19UL16/global_fit/V2M1/Summer19UL16_Run%s_V2M1_DATA_L2L3Residual_AK4PFchs.txt",mera[set]);
      //sout = Form("../JERCProtoLab/Summer19UL16/global_fit/V3M1/Summer19UL16_Run%s_V3M1_DATA_L2L3Residual_AK4PFchs.txt",mera[set]);
      //sout = Form("../JERCProtoLab/Summer19UL16/global_fit/V3M2/Summer19UL16_Run%s_V3M2_DATA_L2L3Residual_AK4PFchs.txt",mera[set]);
      //sout = Form("../JERCProtoLab/Summer19UL16/global_fit/V5M1/Summer19UL16_Run%s_V5M1_DATA_L2L3Residual_AK4PFchs.txt",mera[set]);
      //sout = Form("../JERCProtoLab/Summer19UL16/global_fit/V5M2/Summer19UL16_Run%s_V5M2_DATA_L2L3Residual_AK4PFchs.txt",mera[set]);
      //sout = Form("../JERCProtoLab/Summer19UL16/global_fit/V5M4/Summer19UL16_Run%s_V5M4_DATA_L2L3Residual_AK4PFchs.txt",mera[set]);
      sout = Form("../JERCProtoLab/Summer20UL/global_fit/Summer19UL16_Run%s_VX_DATA_L2L3Residual_AK4PFchs.txt",mera[set]);
    }
    else if (set=="2016BCDEF" || set=="2016BCD" || set=="2016EF") {
      //sin = Form("../JECDatabase/textFiles/Summer19UL16APV_Run%s_V3_DATA/Summer19UL16APV_Run%s_V3_DATA_L2Residual_AK4PFchs.txt",mera[set],mera[set]);
      sin = Form("../JECDatabase/textFiles/Summer19UL16APV_Run%s_V5_DATA/Summer19UL16APV_Run%s_V5_DATA_L2Residual_AK4PFchs.txt",mera[set],mera[set]);
      //
      //sout = Form("../JERCProtoLab/Summer19UL16/global_fit/V3M1/Summer19UL16_Run%s_V3M1_DATA_L2L3Residual_AK4PFchs.txt",mera[set]);
      //sout = Form("../JERCProtoLab/Summer19UL16/global_fit/V3M2/Summer19UL16_Run%s_V3M2_DATA_L2L3Residual_AK4PFchs.txt",mera[set]);
      //sout = Form("../JERCProtoLab/Summer19UL16/global_fit/V5M1/Summer19UL16_Run%s_V5M1_DATA_L2L3Residual_AK4PFchs.txt",mera[set]);
      //sout = Form("../JERCProtoLab/Summer19UL16/global_fit/V5M2/Summer19UL16_Run%s_V5M2_DATA_L2L3Residual_AK4PFchs.txt",mera[set]);
      //sout = Form("../JERCProtoLab/Summer19UL16/global_fit/V5M4/Summer19UL16_Run%s_V5M4_DATA_L2L3Residual_AK4PFchs.txt",mera[set]);
      sout = Form("../JERCProtoLab/Summer20UL/global_fit/Summer19UL16APV_Run%s_VX_DATA_L2L3Residual_AK4PFchs.txt",mera[set]);
    }
    else
      assert(false);
  }
  else if (isUL17) {
    map<string,const char*> mera;
    mera["2017BCDEF"] = "BCDEF";
    mera["2017B"] = "B";
    mera["2017C"] = "C";
    mera["2017D"] = "D";
    mera["2017E"] = "E";
    mera["2017F"] = "F";
    sin = Form("../JECDatabase/textFiles/Summer19UL17_Run%s_V6_DATA/Summer19UL17_Run%s_V6_DATA_L2Residual_AK4PFchs.txt",mera[set],mera[set]);
    if (set=="2017BCDEF") // minitools/mk_mergeL2L3ResTextFiles.C:("2017")
      sin = "textFiles/mergeL2ResTextFiles_2017.txt";
    //sout = Form("../JERCProtoLab/Summer19UL16/global_fit/V3M2/Summer19UL17_Run%s_V3M2_DATA_L2L3Residual_AK4PFchs.txt",mera[set]);
    //sout = Form("../JERCProtoLab/Summer19UL16/global_fit/V5M1/Summer19UL17_Run%s_V5M1_DATA_L2L3Residual_AK4PFchs.txt",mera[set]);
    //sout = Form("../JERCProtoLab/Summer19UL16/global_fit/V5M2/Summer19UL17_Run%s_V5M2_DATA_L2L3Residual_AK4PFchs.txt",mera[set]);
    //sout = Form("../JERCProtoLab/Summer19UL16/global_fit/V5M4/Summer19UL17_Run%s_V5M4_DATA_L2L3Residual_AK4PFchs.txt",mera[set]);
    sout = Form("../JERCProtoLab/Summer20UL/global_fit/Summer19UL17_Run%s_VX_DATA_L2L3Residual_AK4PFchs.txt",mera[set]);
  }
  else if (isLowPU) {
    sin = "CondFormats/JetMETObjects/data/Summer19UL17_V1_SimpleL1_MPF_LOGLIN_L2Residual_pythia8_AK4PFchs.txt";
    //sout = "../JERCProtoLab/Summer19UL17/global_fit/Summer19UL17_RunH_V1M1_DATA_L2L3Residual_AK4PFchs.txt";
    sout = "../JERCProtoLab/Summer20UL/global_fit/Summer19UL17_RunH_VX_DATA_L2L3Residual_AK4PFchs.txt";
  }
  else if (isRun2) {
    sin = "textFiles/mergeL2ResTextFiles_Run2.txt";
    sout = "../JERCProtoLab/Summer20UL/global_fit/Summer19UL_Run2_VX_DATA_L2L3Residual_AK4PFchs.txt";
  }
  else {
    assert(false);
    //sin = Form("textFiles/UL17V2-L2Res+JERSF/%s/Run%s/Summer19UL17_V1_%s_MPF_LOGLIN_L2Residual_pythia8_AK4PFchs.txt",l1,run,l1);
    //sout = Form("textFiles/UL17V2MX-L2L3Res+JERSF/Summer19UL17_Run%s_V2M5_%s_DATA_L2L3Residual_AK4PFchs.txt",run,l1);
  }

  ifstream fin(sin.c_str());
  assert(fin.is_open());
  ofstream fout(sout.c_str());

  // 2018: https://twiki.cern.ch/twiki/bin/view/CMSPublic/PixelOfflinePlotsOctober2018

  string header;
  getline(fin, header);
  if (debug) cout << "Old L2L3Residual header:" << endl;
  if (debug) cout << header << endl;

  header = "{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]/x+[8]*log(x)/x+[9]*(pow(x/[10],[11])-1)/(pow(x/[10],[11])+1)+[12]*pow(x,-0.3051))) Correction L2Relative}";
  if (debug) cout << "New L2L3Residual header:" << endl;
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
    assert(!(etamin==0 && etamax==0));
    if (fabs(etamin)<0.01 && fabs(etamax)<0.01) {
      cout << "sscanf failed!" << endl << flush;
      exit(1);
    }
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
