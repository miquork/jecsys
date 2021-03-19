// Purpose: calibrate Zee and Zmm masses for JEC
//          use minitools/drawZeeVsZmm.C to cross-check parameterization
// run with 'root -l -b -q minitools/drawZmass.C+g'
#include "../tdrstyle_mod15.C"
#include "../tools.C"

#include "TFile.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TLine.h"
#include "TMatrixD.h"

// Plot Z/gamma+jet MPF statistical uncertaity band as reference
bool plotZmmStat = true;
bool plotZeeStat = true;
bool plotGamStat = true;
bool useAlpha100 = true;//false;//true;

bool ispr = true; // PR=fancy plot
const double ptzmax = 1000;
void cleanGraph(TGraphErrors *g, double xmax) {
  for (int i = g->GetN(); i != -1; --i) {
    if (g->GetX()[i]>xmax) g->RemovePoint(i);
  } // for i
} // cleanGraph


void drawZmasses(string run="ABC") {

  setTDRStyle();

  bool isUL18 = (run=="2018ABCD" || run=="2018A" || run=="2018B" ||
		 run=="2018C" || run=="2018D");

  bool isUL17 = (run=="2017BCDEF");

  bool isUL16 = (run=="2016GH");

  TFile *f = new TFile(Form("rootfiles/jecdata%s.root",run.c_str()),"READ");
  assert(f && !f->IsZombie());

  TGraphErrors *gzmmd = (TGraphErrors*)f->Get("data/eta00-13/mass_zmmjet_a30");
  assert(gzmmd);
  TGraphErrors *gzmmm = (TGraphErrors*)f->Get("mc/eta00-13/mass_zmmjet_a30");
  assert(gzmmm);
  //TGraphErrors *gzmmr = (TGraphErrors*)f->Get("ratio/eta00-13/mass_zmmjet_a30");
  //assert(gzmmr);

  TGraphErrors *gzeed = (TGraphErrors*)f->Get("data/eta00-13/mass_zeejet_a30");
  assert(gzeed);
  TGraphErrors *gzeem = (TGraphErrors*)f->Get("mc/eta00-13/mass_zeejet_a30");
  assert(gzeem);
  //TGraphErrors *gzeer = (TGraphErrors*)f->Get("ratio/eta00-13/mass_zeejet_a30");
  //assert(gzeer);

  // MPFchs for statistical uncertainty band
  TGraphErrors *gzmmr = (TGraphErrors*)f->Get("ratio/eta00-13/mpfchs1_zmmjet_a30");
  //assert(gzmmr);
  TGraphErrors *gzeer = (TGraphErrors*)f->Get("ratio/eta00-13/mpfchs1_zeejet_a30");
  //assert(gzeer);
  TGraphErrors *ggamr = (TGraphErrors*)f->Get("ratio/eta00-13/mpfchs1_gamjet_a30");
  //assert(ggamr);

  // Take alpha<1.0 instead of alpha<0.3
  double ptmin = 30.;
  if (useAlpha100) {
    ptmin = 15;

    gzmmd = (TGraphErrors*)f->Get("data/eta00-13/mass_zmmjet_a100");
    gzmmm = (TGraphErrors*)f->Get("mc/eta00-13/mass_zmmjet_a100");

    gzeed = (TGraphErrors*)f->Get("data/eta00-13/mass_zeejet_a100");
    gzeem = (TGraphErrors*)f->Get("mc/eta00-13/mass_zeejet_a100");

    gzmmr = (TGraphErrors*)f->Get("ratio/eta00-13/mpfchs1_zmmjet_a100");
    gzeer = (TGraphErrors*)f->Get("ratio/eta00-13/mpfchs1_zeejet_a100");
    ggamr = (TGraphErrors*)f->Get("ratio/eta00-13/mpfchs1_gamjet_a30");
  }
  assert(gzmmd);
  assert(gzmmm);
  assert(gzeed);
  assert(gzeem);


  cleanGraph(gzmmd,ptzmax);
  cleanGraph(gzmmm,ptzmax);
  //cleanGraph(gzmmr,ptzmax);
  cleanGraph(gzeed,ptzmax);
  cleanGraph(gzeem,ptzmax);

  double ptmax = 1200.*2.;
  //TH1D *hup = new TH1D("hup",";p_{T,Z} (GeV);m_{Z} (GeV)",670,30,ptmax);
  TH1D *hup = new TH1D("hup",";p_{T,Z} (GeV);m_{Z} (GeV)",670,15,ptmax);
  hup->SetMinimum(89.8);//90.0);//80);//90.0);
  hup->SetMaximum(92.8);//93.0);//110);//93.0);
  if (isUL16) {
    hup->SetMinimum(89.8-0.25);
    hup->SetMaximum(92.8+0.5);
  }
  hup->GetXaxis()->SetMoreLogLabels();
  hup->GetXaxis()->SetNoExponent();

  //TH1D *hdw = new TH1D("hdw",";p_{T,Z} (GeV);Data / MC - 1 (%)",670,30,ptmax);
  TH1D *hdw = new TH1D("hdw",";p_{T,Z} (GeV);Data / MC - 1 (%)",670,15,ptmax);
  hdw->SetMinimum((0.990+1e-5-1)*100);
  hdw->SetMaximum((1.025+1e-5-1)*100);
  if (isUL16) {
    hdw->SetMinimum((0.985+1e-5-1)*100);
    hdw->SetMaximum((1.035+1e-5-1)*100);
  }
  hdw->GetXaxis()->SetMoreLogLabels();
  hdw->GetXaxis()->SetNoExponent();

  //extraText = "Private Work";
  map<string, const char*> lumimap;
  //lumimap["A"] = "Run2018A 14.0 fb^{-1}"; //PdmV Analysis TWiki
  //lumimap["B"] = "Run2018B 7.1 fb^{-1}"; //PdmV Analysis TWiki
  //lumimap["C"] = "Run2018C 6.9 fb^{-1}"; //PdmV Analysis TWiki
  //lumimap["D"] = "Run2018D 31.9 fb^{-1}"; //PdmV Analysis TWiki
  //lumimap["ABC"] = "Run2018ABC 28.0 fb^{-1}"; //PdmV Analysis TWiki
  //lumimap["ABCD"] = "Run2018ABCD 59.9 fb^{-1}"; //PdmV Analysis TWiki
  lumimap["BCDEF"] = "2017, 41.5 fb^{-1}"; // for DP note
  lumimap["B"] = "Run2017B, 4.8 fb^{-1}";
  lumimap["C"] = "Run2017C, 9.6 fb^{-1}";
  lumimap["D"] = "Run2017D, 4.2 fb^{-1}";
  lumimap["E"] = "Run2017E, 9.3 fb^{-1}";
  lumimap["F"] = "Run2017F, 13.4 fb^{-1}";
  lumimap["2018ABCD"] = "2018, 59.9 fb^{-1}"; // placeholder
  lumimap["2018A"] = "Run2018A, 14.0 fb^{-1}";
  lumimap["2018B"] = "Run2018B, 7.1 fb^{-1}";
  lumimap["2018C"] = "Run2018C, 6.9 fb^{-1}";
  lumimap["2018D"] = "Run2018D, 31.9 fb^{-1}";
  lumimap["2016BCD"] = "Run2016BCD, 12.9 fb^{-1}";
  lumimap["2016EF"] = "Run2016EF, 6.8 fb^{-1}";
  lumimap["2016GH"] = "Run2016GH, 16.8 fb^{-1}";
  lumimap["2016BCDEFGH"] = "Run2016BCDEFGH, 36.5 fb^{-1}";
  lumi_13TeV = lumimap[run];
  TCanvas *c1 = tdrDiCanvas("c1",hdw,hup,4,11);

  c1->cd(2);
  gPad->SetLogx();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  double mzpdg = 91.2;
  //l->DrawLine(30,mzpdg,ptmax,mzpdg);
  l->DrawLine(15,mzpdg,ptmax,mzpdg);

  tdrDraw(gzeed,"Pz",kFullSquare,kGreen+2);
  tdrDraw(gzeem,"Pz",kOpenSquare,kGreen+2);

  tdrDraw(gzmmd,"Pz",kFullCircle,kRed);
  tdrDraw(gzmmm,"Pz",kOpenCircle,kRed);


  c1->cd(1);
  gPad->SetLogx();

  // Plot originals on the back to see the default fits
  //tdrDraw(gzeer,"Pz",kFullSquare,kGreen+2);
  //tdrDraw(gzmmr,"Pz",kFullCircle,kRed);

  // Plot the old fit used in EOY2017 reprocess.C (f1mzee)
  TF1 *f1eoy17ee = new TF1("f1eoy17ee","(([0]+[1]*log(0.01*x)"
			   "+[2]*pow(log(0.01*x),2))"
			   "-1)*100", 15, ptmax);
  //f1eoy17ee->SetParameters(1.00298, 0.00260, 0.00026); // 2018?
  f1eoy17ee->SetParameters(1.00246, 0.00214, 0.00116); // EOY2017
  f1eoy17ee->SetLineColor(kGray);
  f1eoy17ee->SetLineWidth(2);
  if (isUL17) f1eoy17ee->Draw("SAME");
  //
  TF1 *f1eoy17mm = new TF1("f1eoy17","([0]-1)*100", 15, ptmax);
  f1eoy17mm->SetParameter(0, 0.99854); // EOY2017
  f1eoy17mm->SetLineColor(kGray);
  f1eoy17mm->SetLineWidth(2);
  if (isUL17) f1eoy17mm->Draw("SAME");

  // Plot the old fit used in EOY2018
  TF1 *f1eoy18ee = new TF1("f1eoy18ee","(([0]+[1]*log(0.01*x)"
			   "+[2]*pow(log(0.01*x),2))"
			   "-1)*100", 15, ptmax);
  f1eoy18ee->SetParameters(1.00298, 0.00260, 0.00026); // EOY2018 ABC
  f1eoy18ee->SetLineColor(kGray);
  f1eoy18ee->SetLineWidth(2);
  if (isUL18) f1eoy18ee->Draw("SAME");
  //
  TF1 *f1eoy18mm = new TF1("f1eoy18mm","([0]-1)*100", 15, ptmax);
  f1eoy18mm->SetParameters(0.99851, 0.00000, 0.00000); // EOY18 ABC
  f1eoy18mm->SetLineColor(kGray);
  f1eoy18mm->SetLineWidth(2);
  if (isUL18) f1eoy18mm->Draw("SAME");

  // Plot the old fit used in EOY2016
  TF1 *f1eoy16ee = new TF1("f1eoy16ee","(([0]+[1]*log(0.01*x)"
			   "+[2]*pow(log(0.01*x),2))"
			   "-1)*100", 15, ptmax);
  f1eoy16ee->SetParameters(1.00025, 0.00092, -0.00001); // BCDEFGH (EOY16)
  f1eoy16ee->SetLineColor(kGray);
  f1eoy16ee->SetLineWidth(2);
  if (isUL16) f1eoy16ee->Draw("SAME");
  //
  TF1 *f1eoy16mm = new TF1("f1eoy16mm","([0]-1)*100", 15, ptmax);
  f1eoy16mm->SetParameters(0.99865, 0.00000, 0.00000); // BCDEFGH (EOY16)
  f1eoy16mm->SetLineColor(kGray);
  f1eoy16mm->SetLineWidth(2);
  if (isUL16) f1eoy16mm->Draw("SAME");

  // Plot the new fit used in UL17 reprocess.C (f1mzee)
  double kee = -(4.8*0.329 + 9.6*0.474 + 4.2*0.513 + 9.3*0.517 + 13.4*0.534)
    / (4.8+9.6+4.2+9.3+13.4)*0.01 + 1.00246;
  double kmm = -(4.8*0.181 + 9.6*0.171 + 4.2*0.213 + 9.3*0.190 + 13.4*0.158)
    / (4.8+9.6+4.2+9.3+13.4)*0.01 + 1;
  cout << "kee(UL17 BCDEF) = " << kee << endl;
  cout << "kmm(UL17 BCDEF) = " << kmm << endl;
  TF1 *f1ul17ee = new TF1("f1ul17ee","(([0]+[1]*log(0.01*x)"
			  "+[2]*pow(log(0.01*x),2))"
			  "-1)*100", 15, ptmax);
  //f1ul17ee->SetParameters(1.00298, 0.00260, 0.00026); // UL17B
  //f1ul17ee->SetParameters(0.997557, 0.00214, 0.00116); // UL17B+C+D+E+F
  f1ul17ee->SetParameters(0.99780, 0.00225, 0.00031); // UL17BCDEF-v1
  if (isUL18 || isUL16) {
    f1ul17ee->SetLineColor(kGray);
    f1ul17ee->SetLineStyle(kDashed);
  }
  else
    f1ul17ee->SetLineColor(kBlack);
  f1ul17ee->SetLineWidth(2);
  f1ul17ee->Draw("SAME");
  //
  TF1 *f1ul17mm = new TF1("f1ul17mm","([0]-1)*100", 15, ptmax);
  //f1ul17mm->SetParameter(0, 0.998235); // UL17B+C+D+E+F
  f1ul17mm->SetParameter(0, 0.99821); // UL17BCDEF
  if (isUL18 || isUL16)  {
    f1ul17mm->SetLineColor(kGray);
    f1ul17mm->SetLineStyle(kDashed);
  }
  else 
    f1ul17mm->SetLineColor(kBlack);
  f1ul17mm->SetLineWidth(2);
  f1ul17mm->Draw("SAME");

  // Plot the new fit used in UL18 reprocess.C (f1mzee)
  TF1 *f1ul18ee = new TF1("f1ul18ee","(([0]+[1]*log(0.01*x)"
			  "+[2]*pow(log(0.01*x),2))"
			  "-1)*100", 15, ptmax);
  //f1ul18ee->SetParameters(1.00143, 0.00225, 0.00031); // UL18ABCD-v0
  f1ul18ee->SetParameters(1.00153, 0.00214, -0.00012);
  f1ul18ee->SetLineColor(kBlack);
  f1ul18ee->SetLineWidth(2);
  if (isUL16) {
    f1ul18ee->SetLineStyle(kDotted);
    f1ul18ee->SetLineColor(kGray);
  }
  if (isUL18 || isUL16) f1ul18ee->Draw("SAME");
  //
  TF1 *f1ul18mm = new TF1("f1ul18mm","([0]-1)*100", 15, ptmax);
  f1ul18mm->SetParameters(0.99839, 0.00000, 0.00000); // UL18ABCD
  f1ul18mm->SetLineColor(kBlack);
  f1ul18mm->SetLineWidth(2);
  if (isUL16) {
    f1ul18mm->SetLineStyle(kDotted);
    f1ul18mm->SetLineColor(kGray);
  }
  if (isUL18 || isUL16) f1ul18mm->Draw("SAME");

  // Plot the new fit used in UL16 reprocess.C (f1mzee)
  TF1 *f1ul16ee = new TF1("f1ul16ee","(([0]+[1]*log(0.01*x)"
			  "+[2]*pow(log(0.01*x),2))"
			  "-1)*100", 15, ptmax);
  //f1ul16ee->SetParameters(0.99780, 0.00225, 0.00031); // UL17 placeholder
  //f1ul16ee->SetParameters(1.00110, 0.00225, 0.00031); // UL16GH-UL17 slope
  //f1ul16ee->SetParameters(1.00128, 0.00287, 0.00062); // UL16GH+p2 penalty
  f1ul16ee->SetParameters(1.00175, 0.00349, 0.00159); // UL16GH+eta13+free p2
  f1ul16ee->SetLineColor(kBlack);
  f1ul16ee->SetLineWidth(2);
  if (isUL16) f1ul16ee->Draw("SAME");
  //
  TF1 *f1ul16mm = new TF1("f1ul16mm","([0]-1)*100", 15, ptmax);
  //f1ul16mm->SetParameters(0.99821, 0.00000, 0.00000); // UL17 placeholder
  f1ul16mm->SetParameters(1.00017, 0.00000, 0.00000); // UL16GH
  f1ul16mm->SetLineColor(kBlack);
  f1ul16mm->SetLineWidth(2);
  if (isUL16) f1ul16mm->Draw("SAME");

  l->DrawLine(15,0,ptmax,0);

  // Recreate ratios, because GetListOfFuctions segfaults
  // Also have small shift in pT,Z values within bins
  TGraphErrors *gzeer2 = tools::ratioGraphs(gzeed,gzeem); 
  TGraphErrors *gzmmr2 = tools::ratioGraphs(gzmmd,gzmmm);

  // Increase uncertainty for Zmm mZ kink at pT,Z~70 GeV
  int ik = 3;
  double ek = 0.0000; // no extra (UL17)
  // same as ee * mumu(100)/ee(100):
  //double ek = gzeer2->GetEY()[3] * gzmmr2->GetEY()[4]/gzeer2->GetEY()[4]; // EOY2017
  gzmmr2->SetPointError(ik, gzmmr2->GetEX()[ik],
			sqrt(pow(gzmmr2->GetEY()[ik],2)+pow(ek,2)));

  // Change result to (data/MC-1)*100
  for (int i = 0; i != gzeer2->GetN(); ++i) {
    gzeer2->SetPoint(i, gzeer2->GetX()[i], (gzeer2->GetY()[i]-1)*100);
    gzeer2->SetPointError(i, gzeer2->GetEX()[i], (gzeer2->GetEY()[i])*100);
  }
  for (int i = 0; i != gzmmr2->GetN(); ++i) {
    gzmmr2->SetPoint(i, gzmmr2->GetX()[i], (gzmmr2->GetY()[i]-1)*100);
    gzmmr2->SetPointError(i, gzmmr2->GetEX()[i], (gzmmr2->GetEY()[i])*100);
  }

  // Remove points below ptmin
  for (int i = gzeer2->GetN()-1; i>-1; --i) {
    if (gzeer2->GetX()[i]<ptmin) gzeer2->RemovePoint(i);
  }
  for (int i = gzmmr2->GetN()-1; i>-1; --i) {
    if (gzmmr2->GetX()[i]<ptmin) gzmmr2->RemovePoint(i);
  }

  //TF1 *f1mzee = new TF1("f1zee","([0]+[1]*log(0.01*x)+[2]*pow(log(0.01*x),2)-1)*100",30,ptmax);
  // Add possibility to penalize quadratic term
  TF1 *f1mzee = new TF1("f1zee","(x>15)*(([0]+[1]*log(0.01*x)+[2]*pow(log(0.01*x),2)-1)*100)+(x<15)*[2]",-15,ptmax);
  if (true && !isUL16) { // Penalty for quadratic term
    int n = gzeer2->GetN();
    gzeer2->SetPoint(n, 10, 0.0005);
    //if (isUL16)
    //gzeer2->SetPointError(n, 0, 0.00030); // 1sigma of p2 in a100 fit
    //else 
    gzeer2->SetPointError(n, 0, 0.00014);//0.00021); // 1sigma of p2 in a100 fit
  }
  //f1mzee->FixParameter(2,0);
  f1mzee->SetParameters(1,0.001,0.0001);
  //f1mzee->FixParameter(1, 0.00214); // EOY2017BCDEF
  //f1mzee->FixParameter(2, 0.00116); // EOY2017BCDEF
  if (run!="BCDEF" && run!="2018ABCD" && run!="2016GH") {
    f1mzee->FixParameter(1, 0.00225); // UL2017BCDEF-v1
    f1mzee->FixParameter(2, 0.00031); // UL2017BCDEF-v1
  }
  f1mzee->SetLineColor(kGreen+2);
  gzeer2->Fit(f1mzee,"QRN");

  // Put Zee+jet and gamma+jet statistical uncertainties around the fit
  TGraphErrors *gstatze(0);
  if (gzeer && plotZeeStat) {
    gstatze = (TGraphErrors*)gzeer->Clone();
    for (int i = 0; i != gstatze->GetN(); ++i) {
      gstatze->SetPointError(i,0,100.*gstatze->GetEY()[i]/gstatze->GetY()[i]);
      gstatze->SetPoint(i,gstatze->GetX()[i],f1mzee->Eval(gstatze->GetX()[i]));
    }
    gstatze->SetLineColor(kGreen+2);
    gstatze->SetFillStyle(1001);
    gstatze->SetFillColorAlpha(kGreen-9,0.3);
    gstatze->Draw("E3");
  }
  TGraphErrors *gstatg(0);
  if (ggamr && plotGamStat) {
    gstatg = (TGraphErrors*)ggamr->Clone();
    for (int i = 0; i != gstatg->GetN(); ++i) {
      gstatg->SetPointError(i, 0, 100.*gstatg->GetEY()[i]/gstatg->GetY()[i]);
      gstatg->SetPoint(i, 2*gstatg->GetX()[i], f1mzee->Eval(2*gstatg->GetX()[i]));
    }
    gstatg->SetLineColor(kBlue);
    gstatg->SetFillStyle(1001);
    gstatg->SetFillColorAlpha(kBlue-9,0.3);
    gstatg->Draw("E3");
  }
  
  f1mzee->Draw("SAME");

  TMatrixD emat(f1mzee->GetNpar(),f1mzee->GetNpar());
  gMinuit->mnemat(&emat[0][0],f1mzee->GetNpar());

  // d/dp0=1, d/dp1=log(x), d/dp2=log(x)^2
  // eps =   (d/dp0)^2 m_00 +   (d/dp1)^2 m_11   +   (d/dp2)^2 m_22
  //     + 2*(d/dp0/p1)m_01 + 2*(d/dp0/p2)m_02 + 2*(d/dp1/p2)m_12
  TF1 *f1e = new TF1("f1e","([0]+[1]*log(0.01*x)+[2]*pow(log(0.01*x),2)"
		     "+ [3]*sqrt([4] + pow(log(0.01*x),2)*[5]"
		     "+ pow(log(0.01*x),4)*[6]"
		     "+ 2*log(0.01*x)*[7]+2*pow(log(0.01*x),2)*[8]"
		     "+ 2*pow(log(0.01*x),3)*[9])"
		     " - 1)*100",15,ptmax);
  f1e->SetParameters(f1mzee->GetParameter(0),f1mzee->GetParameter(1),
		     f1mzee->GetParameter(2), +1,
		     emat[0][0], emat[1][1], emat[2][2],
		     emat[0][1], emat[0][2], emat[1][2]);

  f1e->SetLineColor(kGreen-8);
  f1e->DrawClone("SAME");
  f1e->SetParameter(3, -1);
  f1e->DrawClone("SAME");

  // Print out corrections for reprocess.C
  cout <<
    "  // \\BEGIN copy-paste from minitools/drawZmass.C\n"
    "\n"
    "  // Smoothen mass corrections\n"
    "  TF1 *f1mzee = new TF1(\"f1mzee\",\"[0]+[1]*log(0.01*x)+[2]*pow(log(0.01*x),2)\",\n"
    "			     fzeeptmin, fzeeptmax);\n"
    "  TF1 *f1ezee = new TF1(\"f1ezee\",\"sqrt([0]+pow(log(0.01*x),2)*[1]\"\n"
    "                        \"+pow(log(0.01*x),4)*[2]\"\n"
    "                        \"+2*log(0.01*x)*[3]+2*pow(log(0.01*x),2)*[4]\"\n"
    "                        \"+2*pow(log(0.01*x),3)*[5])\",\n"
    "			     fzeeptmin, fzeeptmax);\n"
    "  if (correctZeeMass || correctGamMass) {\n"
    "    if (useFixedFit) {\n";
  //"      // ABCD fit with minitools/drawZmass.C\n";
  if (isUL17)
    cout << Form("      // UL17 Run%s fit with minitools/drawZmass.C\n",
		 run.c_str());
  if (isUL16)
    cout << Form("      // UL16 Run%s fit with minitools/drawZmass.C\n",
		 run.c_str());
  cout << Form(""
	       "      f1mzee->SetParameters(%1.5f, %1.5f, %1.5f);\n"
	       "      f1ezee->SetParameters(%+9.3g, %+9.3g, %+9.3g,\n"
	       "			    %+9.3g, %+9.3g, %+9.3g);\n"
	       "    }\n"
	       "    else\n"
	       "      hmzee->Fit(f1mzee);\n"
	       "\n"
	       "  }\n",
	       f1mzee->GetParameter(0),f1mzee->GetParameter(1),
	       f1mzee->GetParameter(2),
	       emat[0][0], emat[1][1], emat[2][2],
	       emat[0][1], emat[0][2], emat[1][2]);


  //TF1 *f1mzee0 = new TF1("f1zee0","(1+[0]*log(x)+[1]*log(x)*log(x)-1)*100",
  //			 15,ptmax);
  //f1mzee0->SetParameters(0.001,0.0001);
  //f1mzee0->SetLineColor(kGreen+3);
  //gzeer2->Fit(f1mzee0,"QRN");

  TF1 *f1mzmm = new TF1("f1zmm","([0]+[1]*log(0.01*x)+[2]*pow(log(0.01*x),2)-1)*100",
			15,ptmax);
  f1mzmm->SetParameters(1,0,0);
  f1mzmm->FixParameter(2,0);
  f1mzmm->FixParameter(1,0);
  f1mzmm->SetLineColor(kRed);
  gzmmr2->Fit(f1mzmm,"QRN");

  TGraphErrors *gstatzm(0);
  if (gzmmr && plotZmmStat) {
    gstatzm = (TGraphErrors*)gzmmr->Clone();
    for (int i = 0; i != gstatzm->GetN(); ++i) {
      gstatzm->SetPointError(i,0,100.*gstatzm->GetEY()[i]/gstatzm->GetY()[i]);
      gstatzm->SetPoint(i,gstatzm->GetX()[i],f1mzmm->Eval(gstatzm->GetX()[i]));
    }
    gstatzm->SetLineColor(kRed);
    gstatzm->SetFillStyle(1001);
    gstatzm->SetFillColorAlpha(kRed-9,0.3);
    gstatzm->Draw("E3");
  }

  f1mzmm->Draw("SAME");

  TMatrixD emat2(f1mzmm->GetNpar(),f1mzmm->GetNpar());
  gMinuit->mnemat(&emat2[0][0],f1mzmm->GetNpar());
  
  f1e->SetLineStyle(kSolid);
  f1e->SetLineColor(kRed-9);
  f1e->SetParameters(f1mzmm->GetParameter(0),0,0,+1,
		     pow(f1mzmm->GetParError(0),2),0,0, 0,0,0);
  f1e->DrawClone("SAME");
  f1e->SetParameter(3,-1);
  f1e->DrawClone("SAME");

  // Print out corrections for reprocess.C
  cout <<
    "\n"
    "  TF1 *f1mzmm = new TF1(\"f1mzmm\",\"[0]+[1]*log(0.01*x)+[2]*pow(log(0.01*x),2)\",\n"
    "			     fzmmptmin, fzmmptmax);\n"
    "  TF1 *f1ezmm = new TF1(\"f1ezmm\",\"sqrt([0]+pow(log(0.01*x),2)*[1]\"\n"
    "                        \"+pow(log(0.01*x),4)*[2]\"\n"
    "                        \"+2*log(0.01*x)*[3]+2*pow(log(0.01*x),2)*[4]\"\n"
    "                        \"+2*pow(log(0.01*x),3)*[5])\"\n"
    "			     ,fzeeptmin, fzeeptmax);\n"
    "  if (correctZmmMass) {\n"
    "    if (useFixedFit) {\n";
    //"      // ABC fit with minitools/drawZmass.C\n";
  if (isUL17) 
    cout << Form("      // UL17 %s fit with minitools/drawZmass.C\n",
		 run.c_str());
  if (isUL16) 
    cout << Form("      // UL16 %s fit with minitools/drawZmass.C\n",
		 run.c_str());
  cout << Form(""
	       "      f1mzmm->SetParameters(%1.5f, %1.5f, %1.5f);\n"
	       "      f1ezmm->SetParameters(%+9.3g, %+9.3g, %+9.3g,\n"
	       "			    %+9.3g, %+9.3g, %+9.3g);\n"
	       "    }\n"
	       "    else\n"
	       "      hmzmm->Fit(f1mzmm);\n"
	       "\n"
	       "  }\n"
               "\n"
               "  // \\END copy-paste from minitools/drawZmass.C\n",
	       f1mzmm->GetParameter(0),f1mzmm->GetParameter(1),
	       f1mzmm->GetParameter(2),
	       emat2[0][0], emat2[1][1], emat2[2][2],
	       emat2[0][1], emat2[0][2], emat2[1][2]);

  // Draw points last so they stay on top of fits
  tdrDraw(gzeer2,"Pz",kFullSquare,kGreen+2); gzeer2->SetMarkerSize(0.5);
  tdrDraw(gzmmr2,"Pz",kFullCircle,kRed);     gzmmr2->SetMarkerSize(0.5);

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  if (useAlpha100) {
    tex->SetTextSize(0.055);
    tex->DrawLatex(0.35,0.84,"#alpha<1.0");
    tex->SetTextSize(0.045);
  }

  TLegend *leg = tdrLeg(0.20,0.50,0.40,0.74);
  leg->AddEntry(gzeed,"Z(#rightarrow e^{+}e^{-}) data","PL");
  leg->AddEntry(gzmmd,"Z(#rightarrow #mu^{+}#mu^{-}) data","PL");
  leg->AddEntry(gzeem,"Z(#rightarrow e^{+}e^{-}) MC","PL");
  leg->AddEntry(gzmmm,"Z(#rightarrow #mu^{+}#mu^{-}) MC","PL");

  if (plotZmmStat || plotZeeStat || plotGamStat) {
    
    int n(0);
    if (plotZmmStat) ++n;
    if (plotZeeStat) ++n;
    if (plotGamStat) ++n;
    if (run=="BCDEF" || run=="2018ABCD" || run=="2016GH") {
      TLegend *leg2 = tdrLeg(0.50,0.72,0.80,0.72+0.06*n);
      if (plotZeeStat) leg2->AddEntry(gstatzm,"Zee+jet MPF stat.","F");
      if (plotZmmStat) leg2->AddEntry(gstatze,"Z#mu#mu+jet MPF stat.","F");
      if (plotGamStat) leg2->AddEntry(gstatg, "#gamma+jet stat. @2#timesp_{T,Z}","F");
      if (run=="2018ABCD") {
	// Add UL17, EOY18, UL18 lines
	leg2->SetY1(leg2->GetY1()-3*0.06);
	leg2->AddEntry(f1ul17ee,"UL17 Zee","L");
	leg2->AddEntry(f1eoy18ee,"EOY18 Zee","L");
	leg2->AddEntry(f1ul18ee,"UL18 Zee","L");
      }
      if (run=="2016GH") {
	leg2->SetY1(leg2->GetY1()-4*0.055);
	leg2->AddEntry(f1ul16mm,"UL16","L");
	leg2->AddEntry(f1eoy16mm,"EYO16","L");
	leg2->AddEntry(f1ul17mm,"UL17","L");
	leg2->AddEntry(f1ul18mm,"UL18","L");
      }
    }
  }

  // Don't print all the internal information to the PR plots
  if (!ispr) {

    tex->SetTextColor(kGreen+2);
    tex->DrawLatex(0.20,0.40,Form("#chi^{2} / NDF = %1.1f / %d (p3)",
				  f1mzee->GetChisquare(), f1mzee->GetNDF()));
    tex->SetTextColor(kGreen+2);
    tex->DrawLatex(0.20,0.70,Form("p_{0} = %1.5f #pm %1.5f",
				  f1mzee->GetParameter(0),
				  f1mzee->GetParError(0)));
    tex->DrawLatex(0.20,0.65,Form("p_{1} = %1.5f #pm %1.5f",
				  f1mzee->GetParameter(1),
				  f1mzee->GetParError(1)));
    tex->DrawLatex(0.20,0.60,Form("p_{2} = %1.5f #pm %1.5f",
				  f1mzee->GetParameter(2),
				  f1mzee->GetParError(2)));
    
    tex->SetTextSize(0.035);
    tex->DrawLatex(0.20,0.54,Form("%10.3g %10.3g %10.3g",
				  emat[0][0], emat[0][1], emat[0][2]));
    tex->DrawLatex(0.20,0.50,Form("%10.3g %10.3g %10.3g",
				  emat[1][0], emat[1][1], emat[1][2]));
    tex->DrawLatex(0.20,0.46,Form("%10.3g %10.3g %10.3g",
				  emat[2][0], emat[2][1], emat[2][2]));
    tex->SetTextSize(0.045);
    
    tex->SetTextColor(kRed);
    tex->DrawLatex(0.20,0.05,Form("#chi^{2} / NDF = %1.1f / %d (p1)",
				  f1mzmm->GetChisquare(), f1mzmm->GetNDF()));
    tex->DrawLatex(0.20,0.11,Form("p_{0} = %1.5f #pm %1.5f",
				  f1mzmm->GetParameter(0),
				  f1mzmm->GetParError(0)));
  } // !ispr
  if (ispr) {

    tex->SetTextColor(kGreen+2);
    tex->DrawLatex(0.2,0.45,Form("#chi^{2} / NDF = %1.1f / %d (e^{+}e^{-})",
				  f1mzee->GetChisquare(), f1mzee->GetNDF()));
    if (run!="BCDEF" && run!="2018ABCD" && run!="2016GH")
      tex->DrawLatex(0.4,0.85,Form("#Delta^{2}M(e^{+}e^{-})"
                                   " = %+1.3f #pm %1.3f%%",
                                   100.*(f1mzee->GetParameter(0)-
					 f1ul17ee->GetParameter(0)),
					 //f1eoy17ee->GetParameter(0)),
                                   100.*f1mzee->GetParError(0)));

    tex->SetTextColor(kRed);
    tex->DrawLatex(0.2,0.05,Form("#chi^{2} / NDF = %1.1f / %d (#mu^{+}#mu^{-})",
				 f1mzmm->GetChisquare(), f1mzmm->GetNDF()));
    tex->DrawLatex(0.2,0.11,Form("#DeltaM(#mu^{+}#mu^{-}) = %1.3f #pm %1.3f%%",
				 (f1mzmm->GetParameter(0)-1)*100,
				 f1mzmm->GetParError(0)*100));
    if (run!="BCDEF" && run!="2018ABCD" && run!="2016GH")
      tex->DrawLatex(0.4,0.80,Form("#Delta^{2}M(#mu^{+}#mu^{-})"
                                   " = %+1.3f #pm %1.3f%%",
                                   100.*(f1mzmm->GetParameter(0)-
					 f1ul17mm->GetParameter(0)),
					 //f1eoy17mm->GetParameter(0)),
                                   100.*f1mzmm->GetParError(0)));
  }


  c1->SaveAs(Form("pdf/drawZmass_Run%s.pdf",run.c_str()));
} // drawZmasses

void drawZmass() {

  /*
  drawZmasses("B");
  drawZmasses("C");
  drawZmasses("D");
  drawZmasses("E");
  drawZmasses("F");
  drawZmasses("BCDEF");
  */
  /*
  drawZmasses("2018ABCD");
  drawZmasses("2018A");
  drawZmasses("2018B");
  drawZmasses("2018C");
  drawZmasses("2018D");
  */
  //drawZmasses("2017BCDEF");

  drawZmasses("2016GH");
} // drawZmass
