// Purpose: draw flavor response results from Andris to understand details
//          take DY-MG-Py u quark response as new baseline
#include "TFile.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TSpline.h"
#include "TMultiGraph.h"

#include "../tdrstyle_mod22.C"

const double minerr = 0.001*0.5;
double oplus(double a, double b) {
  return sqrt(a*a+b*b);
}
TGraphErrors* cleanRange(TGraphErrors *g, double x1, double x2) {
  g = (TGraphErrors*)g->Clone(Form("%s_clean",g->GetName()));
  for (int i = g->GetN()-1; i != -1; --i) {
    if (g->GetX()[i]<x1 || g->GetX()[i]>x2) g->RemovePoint(i);
  }
  return g;
} // cleanRange
TGraphErrors *merge(TGraphErrors *g1, TGraphErrors *g2) {
  TGraphErrors *g = new TGraphErrors(0);
  for (int i = 0; i != g1->GetN(); ++i) {
    int n = g->GetN();
    g->SetPoint(n, g1->GetX()[i], g1->GetY()[i]);
    //g->SetPointError(n, g1->GetEX()[i], g1->GetEY()[i]);
    g->SetPointError(n, g1->GetEX()[i], oplus(g1->GetEY()[i],minerr));
  }
  for (int i = 0; i != g2->GetN(); ++i) {
    int n = g->GetN();
    g->SetPoint(n, g2->GetX()[i], g2->GetY()[i]);
    //g->SetPointError(n, g2->GetEX()[i], g2->GetEY()[i]);
    g->SetPointError(n, g2->GetEX()[i], oplus(g2->GetEY()[i],minerr));
  }
  return g;
} // merge
TGraphErrors *divide(TGraphErrors *g, TSpline3 *s3) {
  g = (TGraphErrors*)g->Clone(Form("%s_divide",g->GetName()));
  for (int i = 0; i != g->GetN(); ++i) {
    double x = g->GetX()[i];
    g->SetPoint(i, x, g->GetY()[i]/s3->Eval(x));
    g->SetPointError(i, g->GetEX()[i], g->GetEY()[i]/s3->Eval(x));
  }
  return g;
} // divide
void print(TF1 *f1, string s="") {
  cout << Form("  %s->SetParameters(",f1->GetName());
  for (int i = 0; i != f1->GetNpar(); ++i) {
    cout << Form("%s%1.4g",i==0?"":" ,",f1->GetParameter(i));
  }
  cout << ");" << (s!="" ? " // " : "") << s << endl;
} // print

void drawAndrisFlavor(string smc="Pythia") {

  // Set graphical style
  setTDRStyle();
  TDirectory *curdir = gDirectory;

  // Open input file
  TFile *f = new TFile("rootfiles/andris_response_fit_results_root/response_fit_results_L5.root","READ");
  assert(f && !f->IsZombie());
  curdir->cd();

  // Read in reference response curves
  TGraphErrors *gr1(0), *gr2(0), *gr(0), *grc1(0), *grc2(0);
  if (smc=="Herwig") {
    gr1 = (TGraphErrors*)f->Get("DY-MG-Her/response_ud_eta0p0to1p305");assert(gr1);
    gr2 = (TGraphErrors*)f->Get("QCD-MG-Her/response_ud_eta0p0to1p305");assert(gr2);
  }
  else {
    gr1 = (TGraphErrors*)f->Get("DY-MG-Py/response_ud_eta0p0to1p305");assert(gr1);
    gr2 = (TGraphErrors*)f->Get("QCD-MG-Py/response_ud_eta0p0to1p305");assert(gr2);
  }
    
  // Remove biased range and large error bars
  grc1 = cleanRange(gr1,15.,60.);//400.);
  grc2 = cleanRange(gr2,60.,3500.);
  gr = merge(grc1,grc2);

  // Draw background canvas
  TH1D *h = tdrHist("h","Response for ud-jet",0.99,1.09);
  lumi_13TeV = "Run2 Legacy MC";
  extraText = "Private";
  TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);

  TLine *l = new TLine();
  l->SetLineColor(kGray+1);
  l->SetLineStyle(kDashed);
  l->DrawLine(15,1,3500,1);
  
  tdrDraw(gr1,"Pz",kFullCircle,kRed);
  tdrDraw(gr2,"Pz",kFullCircle,kBlue);
  tdrDraw(gr,"Pz",kFullCircle,kBlack); gr->SetMarkerSize(0.6);

  // Fit reference, or create a spline out of it, so we can renormalize
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr1);
  mg->Add(gr2);

  // Fit smooth response shape as well as alternative to splines
  TF1 *f1 = new TF1("f1","[0]+[1]*pow(x,[2])"
		    "+[3]*pow(x/1000.,[4])"
		    "+[5]/x",15,3500);
  f1->SetParameters(0.8657 ,-0.8603 ,-0.2572 ,0.2917 ,-0.1611 ,0.8741); // ud
  f1->SetLineColor(kRed-9);
  f1->DrawClone("SAME");
  //
  gr->Fit(f1,"RN");
  f1->SetLineColor(kRed);
  f1->DrawClone("SAME");
  print(f1,"ud");
  
  TSpline3 *s3z = new TSpline3("s3z",gr1);
  s3z->SetLineColor(kRed);
  s3z->Draw("SAME");
  
  TSpline3 *s3q = new TSpline3("s3q",gr2);
  s3q->SetLineColor(kBlue);
  s3q->Draw("SAME");

  TSpline3 *s3 = new TSpline3("s3",gr);
  s3->Draw("SAME");
  
  TLegend *leg1 = tdrLeg(0.60,0.85-0.05*3,0.85,0.85);
  leg1->AddEntry(gr1,"DY-MG-Her","PLE");
  leg1->AddEntry(gr2,"QCD-MG-Her","PLE");
  leg1->AddEntry(gr,"Combination","PLE");

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.040);
  if (smc=="Herwig")
    tex->DrawLatex(0.35,0.86,"MadGraph + Herwig HS1 (DY+QCD)");
  else
    tex->DrawLatex(0.35,0.86,"MadGraph + Pythia CP5 (DY+QCD)");
  tex->DrawLatex(0.35,0.81,"|#eta| < 1.3");
  
  gPad->SetLogx();
  if (smc=="Herwig")
    c1->SaveAs("pdf/drawAndrisFlavor_ref_Herwig.pdf");
  else
    c1->SaveAs("pdf/drawAndrisFlavor_ref_Pythia.pdf");

  // Start adding ratios of various responses
  // Up quark
  TGraphErrors *gu1(0), *gu2(0), *guc1(0), *guc2(0), *gu(0);
  if (smc=="Herwig") {
    gu1 = (TGraphErrors*)f->Get("DY-MG-Her/response_u_eta0p0to1p305");assert(gu1);
    gu2 = (TGraphErrors*)f->Get("QCD-MG-Her/response_u_eta0p0to1p305");assert(gu2);
  }
  else {
    gu1 = (TGraphErrors*)f->Get("DY-MG-Py/response_u_eta0p0to1p305");assert(gu1);
    gu2 = (TGraphErrors*)f->Get("QCD-MG-Py/response_u_eta0p0to1p305");assert(gu2);
  }
  guc1 = cleanRange(gu1,15.,60.);
  guc2 = cleanRange(gu2,60.,3500.);
  gu = merge(guc1,guc2);

  gu1 = divide(gu1,s3z);//DY
  gu2 = divide(gu2,s3q);//QCD
  gu = divide(gu,s3);
  
  // Down quark
  TGraphErrors *gd1(0), *gd2(0), *gdc1(0), *gdc2(0), *gd(0);
  if (smc=="Herwig") {
    gd1 = (TGraphErrors*)f->Get("DY-MG-Her/response_d_eta0p0to1p305");assert(gd1);
    gd2 = (TGraphErrors*)f->Get("QCD-MG-Her/response_d_eta0p0to1p305");assert(gd2);
  }
  else {
    gd1 = (TGraphErrors*)f->Get("DY-MG-Py/response_d_eta0p0to1p305");assert(gd1);
    gd2 = (TGraphErrors*)f->Get("QCD-MG-Py/response_d_eta0p0to1p305");assert(gd2);
  }
  gdc1 = cleanRange(gd1,15.,60.);
  gdc2 = cleanRange(gd2,60.,3500.);
  gd = merge(gdc1,gdc2);

  gd1 = divide(gd1,s3z);//DY
  gd2 = divide(gd2,s3q);//QCD
  gd = divide(gd,s3);

  // Strange quark
  TGraphErrors *gs1(0), *gs2(0), *gsc1(0), *gsc2(0), *gs(0);
  if (smc=="Herwig") {
    gs1 = (TGraphErrors*)f->Get("DY-MG-Her/response_s_eta0p0to1p305");assert(gs1);
    gs2 = (TGraphErrors*)f->Get("QCD-MG-Her/response_s_eta0p0to1p305");assert(gs2);
  }
  else {
    gs1 = (TGraphErrors*)f->Get("DY-MG-Py/response_s_eta0p0to1p305");assert(gs1);
    gs2 = (TGraphErrors*)f->Get("QCD-MG-Py/response_s_eta0p0to1p305");assert(gs2);
  }
    
  gsc1 = cleanRange(gs1,17.,60.);
  gsc2 = cleanRange(gs2,60.,3500.);
  gs = merge(gsc1,gsc2);

  gs1 = divide(gs1,s3z);//DY
  gs2 = divide(gs2,s3q);//QCD
  gs = divide(gs,s3);

  // Charm quark
  TGraphErrors *gc1(0), *gc2(0), *gcc1(0), *gcc2(0), *gc(0);
  if (smc=="Herwig") {
    gc1 = (TGraphErrors*)f->Get("DY-MG-Her/response_c_eta0p0to1p305");assert(gc1);
    gc2 = (TGraphErrors*)f->Get("QCD-MG-Her/response_c_eta0p0to1p305");assert(gc2);
  }
  else {
    gc1 = (TGraphErrors*)f->Get("DY-MG-Py/response_c_eta0p0to1p305");assert(gc1);
    gc2 = (TGraphErrors*)f->Get("QCD-MG-Py/response_c_eta0p0to1p305");assert(gc2);
  }
  //gcc1 = cleanRange(gc1,15.,60.);
  //gcc2 = cleanRange(gc2,60.,3500.);
  gcc1 = cleanRange(gc1,15.,15.);
  gcc2 = cleanRange(gc2,20.,3500.);
  gc = merge(gcc1,gcc2);

  gc1 = divide(gc1,s3z);//DY
  gc2 = divide(gc2,s3q);//QCD
  //gc = divide(gc,s3);
  gc = divide(gc,s3q);//QCD

  // Bottom quark
  TGraphErrors *gb1(0), *gb2(0), *gbc1(0), *gbc2(0), *gb(0);
  if (smc=="Herwig") {
    gb1 = (TGraphErrors*)f->Get("DY-MG-Her/response_b_eta0p0to1p305");assert(gb1);
    gb2 = (TGraphErrors*)f->Get("QCD-MG-Her/response_b_eta0p0to1p305");assert(gb2);
  }
  else {
    gb1 = (TGraphErrors*)f->Get("DY-MG-Py/response_b_eta0p0to1p305");assert(gb1);
    gb2 = (TGraphErrors*)f->Get("QCD-MG-Py/response_b_eta0p0to1p305");assert(gb2);
  }
  //gbc1 = cleanRange(gb1,15.,60.);
  //gbc2 = cleanRange(gb2,60.,3500.);
  gbc1 = cleanRange(gb1,15.,15.);
  gbc2 = cleanRange(gb2,25.,3500.);
  gb = merge(gbc1,gbc2);

  gb1 = divide(gb1,s3z);//DY
  gb2 = divide(gb2,s3q);//QCD
  //gb = divide(gb,s3);
  gb = divide(gb,s3q);//QCD

  // Gluon
  TGraphErrors *gg1(0), *gg2(0), *ggc1(0), *ggc2(0), *gg(0);
  if (smc=="Herwig") {
    gg1 = (TGraphErrors*)f->Get("DY-MG-Her/response_g_eta0p0to1p305");assert(gg1);
    gg2 = (TGraphErrors*)f->Get("QCD-MG-Her/response_g_eta0p0to1p305");assert(gg2);
  }
  else {
    gg1 = (TGraphErrors*)f->Get("DY-MG-Py/response_g_eta0p0to1p305");assert(gg1);
    gg2 = (TGraphErrors*)f->Get("QCD-MG-Py/response_g_eta0p0to1p305");assert(gg2);
  }
  ggc1 = cleanRange(gg1,15.,60.);
  ggc2 = cleanRange(gg2,60.,3500.);
  gg = merge(ggc1,ggc2);

  gg1 = divide(gg1,s3);
  gg2 = divide(gg2,s3q);//QCD
  gg = divide(gg,s3);

  TF1 *fr = new TF1("fr","1+[0]*pow(x,[1])"
		    "+[2]/x"
		    "+[3]*log(x)/x"
		    "+[4]*pow(x,fabs([5]))"
		    ,15,3500.);
  fr->SetParameters(-1,-0.3,0.01,0.0001, 1e-6,1);
  fr->FixParameter(4,0);

  TF1 *fr2 = new TF1("fr2","1+[0]*pow(x,[1])"
		     "+[2]*TMath::Gaus(log(x),[3],[4])"
		     ,15,3500.);
  fr2->SetParameters(-1,-0.3, 
		     -0.025,log(150.),0.3,
		     1e-6,1);
  
  // Draw background canvas
  TH1D *h2 = tdrHist("h2","Response / ud-jet",0.92,1.02);
  TCanvas *c2 = tdrCanvas("c2",h2,4,11,kSquare);

  l->DrawLine(15,1,3500,1);
  tdrDraw(gu,"Pz",kFullCircle,kBlack); gu->SetMarkerSize(0.6);

  //tdrDraw(gd1,"Pz",kOpenCircle,kMagenta);
  //tdrDraw(gd2,"Pz",kOpenCircle,kMagenta+1);
  tdrDraw(gd,"Pz",kOpenCircle,kMagenta+1); gd->SetMarkerSize(0.6);

  //tdrDraw(gs1,"Pz",kFullSquare,kOrange);
  //tdrDraw(gs2,"Pz",kFullSquare,kOrange+1);
  tdrDraw(gs,"Pz",kFullSquare,kOrange+1); gs->SetMarkerSize(0.6);

  //tdrDraw(gc1,"Pz",kOpenSquare,kGreen);
  //tdrDraw(gc2,"Pz",kOpenSquare,kGreen+1);
  tdrDraw(gc,"Pz",kOpenSquare,kGreen+2); gc->SetMarkerSize(0.6);

  //tdrDraw(gb1,"Pz",kFullTriangleUp,kRed-9);
  //tdrDraw(gb2,"Pz",kFullTriangleUp,kRed-8);
  tdrDraw(gb,"Pz",kFullTriangleUp,kRed); gb->SetMarkerSize(0.6);

  //tdrDraw(gg1,"Pz",kOpenTriangleUp,kBlue-9);
  //tdrDraw(gg2,"Pz",kOpenTriangleUp,kBlue-8);
  tdrDraw(gg,"Pz",kOpenTriangleUp,kBlue); gg->SetMarkerSize(0.6);

  //fr->SetParameters(-0.9431 ,-0.7933 ,0.9809 ,0.1929 ,2e-05 ,0.3963); // gd (u-jet)
  fr->SetParameters(-0.9367 ,-0.8106 ,0.8523 ,0.2358 ,2e-05 ,0.2892); // gd
  fr->SetLineColor(kMagenta-9);
  fr->DrawClone("SAME");
  //
  gd->Fit(fr,"QRN");
  fr->SetLineColor(kMagenta+1);
  fr->DrawClone("SAME");
  print(fr,"gd");

  //fr->SetParameters(-0.1172 ,-0.8942 ,0.1084 ,0.01739 ,2e-05 ,1.119e-05); // gu (u-jet)
  fr->SetParameters(-0.1399 ,-0.9999 ,-0.06641 ,0.1023 ,2e-05 ,1.007e-06); // gu
  fr->SetLineColor(kGray+1);
  fr->DrawClone("SAME");
  //
  gu->Fit(fr,"QRN");
  fr->SetLineColor(kBlack);
  fr->DrawClone("SAME");
  print(fr,"gu");

  //fr->SetParameters(-0.03003 ,-0.03497 ,0.5745 ,-0.345 ,2e-05 ,0.8265); // gc (u-jet)
  fr->SetParameters(-0.02159 ,0.01307 ,0.7188 ,-0.3794 ,2e-05 ,0.8394); // gc
  fr->SetLineColor(kGreen-9);
  fr->DrawClone("SAME");
  //
  gc->Fit(fr,"QRN");
  fr->SetLineColor(kGreen+2);
  fr->DrawClone("SAME");
  print(fr,"gc");
  
  //fr->SetParameters(-0.00426 ,0.2732 ,2.86 ,-1.379 ,2e-05 ,0.9152); // gb (u-jet)
  fr->SetParameters(-0.003196 ,0.3194 ,2.886 ,-1.356 ,2e-05 ,0.9292); // gb
  fr->SetLineColor(kRed-9);
  fr->DrawClone("SAME");
  //
  gb->Fit(fr,"QRN");
  fr->SetLineColor(kRed);
  fr->DrawClone("SAME");
  print(fr,"gb");
  
  //fr->SetParameters(-0.01607 ,0.104 ,0.9676 ,-0.6181 ,2e-05 ,0.8872); // gg (u-jet)
  fr->SetParameters(-0.01472 ,0.1171 ,0.8334 ,-0.5488 ,2e-05 ,0.8902); // gg
  fr->SetLineColor(kBlue-9);
  fr->DrawClone("SAME");
  //
  gg->Fit(fr,"QRN");
  fr->SetLineColor(kBlue);
  fr->DrawClone("SAME");
  print(fr,"gg");


  fr2->SetParameters(-0.03628 ,-0.1735 ,-0.01336 ,4.918 ,-0.9191); // gs
  fr2->SetLineColor(kOrange-9);
  fr2->DrawClone("SAME");
  //
  gs->Fit(fr2,"QRN");
  fr2->SetLineColor(kOrange+1);
  fr2->DrawClone("SAME");
  print(fr2,"gs");
  //
  fr2->SetParameter(2,0);
  fr2->SetLineStyle(kDotted);
  fr2->DrawClone("SAME");
  
  TLegend *leg2 = tdrLeg(0.70,0.17,0.90,0.17+0.04*8);
  leg2->AddEntry(gu,"u-jet","PLE");
  leg2->AddEntry(gd,"d-jet","PLE");
  leg2->AddEntry(gs,"s-jet","PLE");
  leg2->AddEntry(gc,"c-jet","PLE");
  leg2->AddEntry(gb,"b-jet","PLE");
  leg2->AddEntry(gg,"g-jet","PLE");

  if (smc=="Herwig")
    tex->DrawLatex(0.35,0.86,"MadGraph + Herwig HS1 (DY+QCD)");
  else
    tex->DrawLatex(0.35,0.86,"MadGraph + Pythia CP5 (DY+QCD)");
  tex->DrawLatex(0.35,0.81,"|#eta| < 1.3");
    
  
  gPad->SetLogx();
  gPad->RedrawAxis();
  if (smc=="Herwig")
    c2->SaveAs("pdf/drawAndrisFlavor_flavors_Herwig.pdf");
  else
    c2->SaveAs("pdf/drawAndrisFlavor_flavors_Pythia.pdf");
  
} // void drawAndrisFlavor
