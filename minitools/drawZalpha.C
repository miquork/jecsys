// Purpose: Study parton correction using ratio of pTbal and MPF in Z+jet.
//          Ideally constrain alpha_S,FSR in MC to reduce uncertainty for m_top
// Run with 'root -l minitools/drawZalpha.C'
#include "TFile.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TLine.h"
//#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TFitResult.h"

#include <map>
#include <string>

#include "../tdrstyle_mod15.C"
#include "../tools.C"


bool correctUnclustMPF = true; // using fCMPF
bool correctGluonMPF = true;   // using fCMPF
bool correctPtBalUE = false;//true;
bool correctQuarkResp = false;//true;

bool drawPS = false;//true;//false;
bool drawData = true;//false;

// Replace FSR fit by pTbal/MPF(alpha<0.3) ratio
const bool patchLowPtUnclus = true;
const double ptpatch = 105;//130;//105;
const double ptthr = 17.5;//18; // min pT for alpha_max
const double pttype1 = 15.; // min pT for type-I MET
// ue_dg: use search

const double etamax = 2.5;//1.3;//2.5
const double amin = 0.;
const double amax = 0.35;//0.275;//0.35

//pair<double,double> drawZalphas(double ptref);
struct FSR {
  double data;
  double data_err;
  double data_chi2;
  double data_mpf;
  double data_mpferr;
  double mc;
  double mc_err;
  double mc_chi2;
  double mc_mpf;
  double mc_mpferr;
  double ps0;
  double ps1;
  double ps2;
  double ps3;
  double ps0_chi2;
  double ps1_chi2;
  double ps2_chi2;
  double ps3_chi2;
};
//pair<double,double> drawZalphas(double ptref);
FSR drawZalphas(double ptref);

TGraphErrors *clipGraph(TGraphErrors *gref, double a, double b) {
  TGraphErrors *g = (TGraphErrors*)gref->Clone();
  for (int i = g->GetN()-1; i != -1; --i) {
    if (g->GetX()[i]<a || g->GetX()[i]>b) g->RemovePoint(i);
  } // for i
  return g;
} // clipGraph

TGraphErrors *fitErrGraph(TF1 *f1, TMatrixDSym &emat, TGraphErrors *gref) {

  assert(f1);
  assert(gref);
  //assert(emat);

  TGraphErrors *g = new TGraphErrors(gref->GetN()+2);
  for (int i = 0; i != g->GetN(); ++i) {
    double x(0);
    if (i==0) x = f1->GetXmin();
    else if (i==g->GetN()-1) x = f1->GetXmax();
    else x = g->GetX()[i-1];
    double y = f1->Eval(x);
    double ey2(0);

    // partial derivatives as differentials
    double dy[f1->GetNpar()];
    for (int j = 0; j != f1->GetNpar(); ++j) {
      double dp = 0.1*sqrt(emat[j][j]);
      double p = f1->GetParameter(j);
      f1->SetParameter(j, p+dp);
      double yup = f1->Eval(x);
      f1->SetParameter(j, p-dp);
      double ydw = f1->Eval(x);
      f1->SetParameter(j, p);
      dy[j] = (dp!=0 ? (yup-ydw)/(2*dp) : 0);
    }

    for (int j = 0; j != f1->GetNpar(); ++j) {
      for (int k = 0; k != f1->GetNpar(); ++k) {
	ey2 += emat[j][k]*dy[j]*dy[k];
      }
    }
    double ey = sqrt(ey2);
    g->SetPoint(i, x, y);
    g->SetPointError(i, 0, ey);
  } // for i
  
  return g;
} // fitErrGraph

// Correct MPF for unclustered pT
//TGraphErrors *
//void fCMPF(TGraphErrors *g, TGraphErrors *gu, TGraphErrors *gn) {  
void fCMPF(TGraphErrors *g, double ptz,
	   TGraphErrors *g1, TGraphErrors *gn, TGraphErrors *gu,
	   bool don, bool dou) {  

  assert(g);
  assert(g1 || !don);
  assert(gn || !don);
  assert(gu || !dou);
  if (don) assert(g->GetN()==g1->GetN());
  if (don) assert(g->GetN()==gn->GetN());
  if (dou) assert(g->GetN()==gu->GetN());

  for (int i = 0; i != g->GetN(); ++i) {

    assert(!g1 || fabs(g->GetX()[i]-g1->GetX()[i])<0.01);
    assert(!gn || fabs(g->GetX()[i]-gn->GetX()[i])<0.01);
    assert(!gu || fabs(g->GetX()[i]-gu->GetX()[i])<0.01);
    double x = g->GetX()[i];
    double R = g->GetY()[i];
    double eR = g->GetEY()[i];
    double P1 = (g1 ? g1->GetY()[i] : 0);
    double eP1 = (g1 ? g->GetEY()[i] : 0);
    double Pn = (gn ? gn->GetY()[i] : 0);
    double ePn = (gn ? gn->GetEY()[i] : 0);
    double Pu = (gu ? gu->GetY()[i] : 0);
    double ePu = (gu ? gu->GetEY()[i] : 0);
    const double R1 = 1.0;
    const double ru = 0.65;
    //const double rn = 0.92;//0.95; // do later vs pT
    // NB: use R_g at low pT ([15,alpha*pTz]) relative to R_incl at pTz
    //const double ptue = 1; // GeV
    //P1 -= ptue/ptz;
    //Pu += ptue/ptz;
    //double Rc = R / ( 1 + (rn-1) * Pn/rn + (ru-1) * Pu/ru);
    //double eRc = sqrt(pow(eR*Rc/R,2) + pow(ePn*(1-rn)*Rc*Rc/R,2) + pow(ePu*(1-1/ru)*Rc*Rc/R,2));
    //
    // Method on 2020/03/02
    //double Rc = R / ( 1 + (rn-1) * (1-P1/R1) + (1-rn/ru) * Pu/R1);
    //double eRc = sqrt(pow(eR*Rc/R,2) + pow(eP1*(rn-1)*Rc*Rc/R,2) + pow(ePu*(1-rn/ru)/R1*Rc*Rc/R,2));
    //
    // Function outlined for JERC talk on 2020/03/02
    double alpha = 0.3;
    double ptmin = 15;
    double ptg = 0.5*(ptmin + max(alpha*ptz,ptmin));
    double rzjet = 1 + 0.01*( -1.52 - 0.37*pow(0.01*ptz,-2.089) );
    double rgjet = 1 + 0.01*( -2.58 - 2.46*pow(0.01*ptg,-1.157) );
    double rn = rgjet / rzjet;
    double fn = Pn / (rn*R1);
    double fu = Pu / (ru*R1);
    double Rc = R / ( 1 + (rn-1)*fn + (ru-1)*fu);
    double eRc = sqrt(pow(eR*Rc/R,2) + pow(eP1*(rn-1)*Rc*Rc/R,2) + pow(ePu*(1-rn/ru)/R1*Rc*Rc/R,2)); // update
    //
    // from estimateUnclust.C:
    //double P1 = rdj1->GetY()[i]-ptue/ptz;
    //double Pu = rdunc->GetY()[i]+ptue/ptz;
    // double kbias = (1 + (rn-1)*(1-P1/R1) + (1-rn/ru)*Pu/R1);
    // double Rcorr = R / kbias;
    g->SetPoint(i, x, Rc);
    g->SetPointError(i, 0, eRc);
  } // for i

} // fCMPF


TH1D *_hgf(0);
TF1 *_fas(0);
void drawZalpha() {

  // arXiv:1604.08082 Eq.(5.6)
  // Wikipedia: Lambda_QCD = 218 +/- 24 MeV
  TF1 *fas = new TF1("fas","4*TMath::Pi()/[0]*(1./log(x*x/([1]*[1])) + [1]*[1]/([1]*[1] - x*x))",30,3000);
  //fas->SetParameters(-0.038,0.218);
  //fas->SetParameters(-0.041,0.218); // before UE change Zr->Dg
  //fas->SetParameters(-0.0385,0.218); // after UE change Zr->Dg
  const double alpharef = 0.3;
  fas->SetParameters(-0.035/alpharef,0.218); // rmpfjet1

  _fas = fas;

  double vpt[] = {35,45,55,65,78,95,118,152,202,265, 350,450,600,850};
  const int npt = sizeof(vpt)/sizeof(vpt[0]);
  TGraphErrors *gd = new TGraphErrors(npt);
  TGraphErrors *gd2 = new TGraphErrors(npt);
  TGraphErrors *gm = new TGraphErrors(npt);
  TGraphErrors *gm2 = new TGraphErrors(npt);
  TGraphErrors *gps0 = new TGraphErrors(npt);
  TGraphErrors *gps1 = new TGraphErrors(npt);
  TGraphErrors *gps2 = new TGraphErrors(npt);
  TGraphErrors *gps3 = new TGraphErrors(npt);
  //
  TGraph *gchi2dt = new TGraph(npt);
  TGraph *gchi2mc = new TGraph(npt);
  TGraph *gchi2ps0 = new TGraph(npt);
  TGraph *gchi2ps1 = new TGraph(npt);
  TGraph *gchi2ps2 = new TGraph(npt);
  TGraph *gchi2ps3 = new TGraph(npt);
  //
  TGraphErrors *gmpfdt = new TGraphErrors(npt);
  TGraphErrors *gmpfmc = new TGraphErrors(npt);
  TGraphErrors *gmpfr = new TGraphErrors(npt);

  for (int i = 0; i != npt; ++i) {
    //pair<double,double> fsr = 
    FSR fsr = drawZalphas(vpt[i]);
    double pt = vpt[i];
    gd->SetPoint(i, pt, fsr.data);//fsr.first);
    gd->SetPointError(i, 0, fsr.data_err);//fsr.second);
    gchi2dt->SetPoint(i, pt, fsr.data_chi2);
    gmpfdt->SetPoint(i, pt, fsr.data_mpf);
    gmpfdt->SetPointError(i, 0, fsr.data_mpferr);
    gm->SetPoint(i, pt, fsr.mc);//second);
    gm->SetPointError(i, 0, fsr.mc_err);//second);
    gchi2mc->SetPoint(i, pt, fsr.mc_chi2);
    gmpfmc->SetPoint(i, pt, fsr.mc_mpf);
    gmpfmc->SetPointError(i, 0, fsr.mc_mpferr);
    //
    gmpfr->SetPoint(i, pt, 100*((1+fsr.data_mpf/100) / (1+fsr.mc_mpf/100) - 1));
    gmpfr->SetPointError(i, 0, sqrt(pow(fsr.data_mpferr,2)
				    + pow(fsr.mc_mpferr,2)));
    //
    gchi2ps0->SetPoint(i, pt, fsr.ps0_chi2);
    gchi2ps1->SetPoint(i, pt, fsr.ps1_chi2);
    gchi2ps2->SetPoint(i, pt, fsr.ps2_chi2);
    gchi2ps3->SetPoint(i, pt, fsr.ps3_chi2);

    /*
    // Map points to pure quark sample
    double fg = _hgf->GetBinContent(_hgf->FindBin(vpt[i]));
    //double CX = fg*2. + (1-fg)*4./3.; // bug: C_A=3, not 2! 
    double CX = fg*3. + (1-fg)*4./3.; 
    double ke = 1;//2
    gd2->SetPoint(i, pt, fsr.data * pow(4./3. / CX, ke));
    gd2->SetPointError(i, 0, fsr.data_err * pow(4./3. / CX, ke));
    gm2->SetPoint(i, pt, fsr.mc * pow(4./3. / CX, ke));
    gm2->SetPointError(i, 0, fsr.mc_err * pow(4./3. / CX, ke));
    //
    double dfsr = -4;
    gps0->SetPoint(i, pt, fsr.ps0 * pow(4./3. / CX, ke) + dfsr);
    gps1->SetPoint(i, pt, fsr.ps1 * pow(4./3. / CX, ke) + dfsr);
    gps2->SetPoint(i, pt, fsr.ps2 * pow(4./3. / CX, ke) + dfsr);
    gps3->SetPoint(i, pt, fsr.ps3 * pow(4./3. / CX, ke) + dfsr);
    */
    // Drop mapping
    gd2->SetPoint(i, pt, fsr.data);
    gd2->SetPointError(i, 0, fsr.data_err);
    gm2->SetPoint(i, pt, fsr.mc);
    gm2->SetPointError(i, 0, fsr.mc_err);
    //
    gps0->SetPoint(i, pt, fsr.ps0);
    gps1->SetPoint(i, pt, fsr.ps1);
    gps2->SetPoint(i, pt, fsr.ps2);
    gps3->SetPoint(i, pt, fsr.ps3);
  }

  //TCanvas *c2 = new TCanvas("c2","c2",600,600);
  TH1D *h2 = new TH1D("h2",";p_{T,Z} (GeV);FSR+ISR (%)",970,30,1000);
  h2->GetXaxis()->SetMoreLogLabels();
  h2->GetXaxis()->SetNoExponent();
  h2->SetMaximum(0);
  h2->SetMinimum(-15);//-45);

  TH1D *h2dw = new TH1D("h2dw",";p_{T,Z} (GeV);Data-MC (%)",970,30,1000);
  h2dw->GetXaxis()->SetMoreLogLabels();
  h2dw->GetXaxis()->SetNoExponent();
  h2dw->SetMaximum(+1.5);//+4.0);
  h2dw->SetMinimum(-1.5);//-4.0);

  TH1D *h2b = new TH1D("h2b",";p_{T,Z} (GeV);FSR+ISR (%)",970,30,1000);
  h2b->GetXaxis()->SetMoreLogLabels();
  h2b->GetXaxis()->SetNoExponent();
  h2b->SetMaximum(0);
  h2b->SetMinimum(-15);//-45);

  TH1D *h2bdw = new TH1D("h2bdw",";p_{T,Z} (GeV);Diff. (%)",970,30,1000);
  h2bdw->GetXaxis()->SetMoreLogLabels();
  h2bdw->GetXaxis()->SetNoExponent();
  h2bdw->SetMaximum(+1.5);//+4.0);
  h2bdw->SetMinimum(-1.5);//-4.0);

  lumi_13TeV = "Run 2, 136.5 fb^{-1}";
  //TCanvas *c2b = tdrCanvas("c2b",h2b,4,11,kSquare);
  TCanvas *c2b = tdrDiCanvas("c2b",h2b,h2bdw,4,11);

  c2b->cd(1);
  gPad->SetLogx();

  TGraphErrors *gps = new TGraphErrors(gps0->GetN());
  TGraphErrors *gpsd = new TGraphErrors(gps0->GetN());
  TGraphErrors *gpsd0 = new TGraphErrors(gps0->GetN());
  TGraphErrors *gpsd1 = new TGraphErrors(gps0->GetN());
  TGraphErrors *gpsd2 = new TGraphErrors(gps0->GetN());
  TGraphErrors *gpsd3 = new TGraphErrors(gps0->GetN());

  for (int i = 0; i != gps->GetN(); ++i) {

    double x = gps0->GetX()[i];
    double ex = gps0->GetEX()[i];
    double yref = 0.5*(gps0->GetY()[i] + gps2->GetY()[i]);
    double eyref = 0.5*(gps0->GetEY()[i] + gps2->GetEY()[i]);
    gps->SetPoint(i, x, yref);
    gps->SetPointError(i, ex, eyref);
    //
    gpsd->SetPoint(i, x, gps->GetY()[i]-yref);
    gpsd->SetPointError(i, ex, 0.);
    gpsd0->SetPoint(i, x, gps0->GetY()[i]-yref);
    gpsd0->SetPointError(i, ex, sqrt(pow(gps0->GetEY()[i],2)+eyref*eyref));
    gpsd1->SetPoint(i, x, gps1->GetY()[i]-yref);
    gpsd1->SetPointError(i, ex, sqrt(pow(gps1->GetEY()[i],2)+eyref*eyref));
    gpsd2->SetPoint(i, x, gps2->GetY()[i]-yref);
    gpsd2->SetPointError(i, ex, sqrt(pow(gps2->GetEY()[i],2)+eyref*eyref));
    gpsd3->SetPoint(i, x, gps3->GetY()[i]-yref);
    gpsd3->SetPointError(i, ex, sqrt(pow(gps3->GetEY()[i],2)+eyref*eyref));
  }

  fas->Draw("SAME");
  tdrDraw(gps,"Pz",kFullSquare,kBlue);
  tdrDraw(gps0,"Pz",kOpenTriangleDown,kRed+0);
  tdrDraw(gps1,"Pz",kFullTriangleDown,kRed+1);
  tdrDraw(gps2,"Pz",kOpenTriangleUp,kRed+2);
  tdrDraw(gps3,"Pz",kFullTriangleUp,kRed+3);

  TLegend *leg2b = tdrLeg(0.50,0.60,0.70,0.90);
  leg2b->AddEntry(gps,"Reference","P");
  leg2b->AddEntry(gps0,"ISR=0.5, FSR=1","P");
  leg2b->AddEntry(gps1,"ISR=1,   FSR=0.5","P");
  leg2b->AddEntry(gps2,"ISR=2,   FSR=1","P");
  leg2b->AddEntry(gps3,"ISR=1,   FSR=2","P");

  leg2b->AddEntry(fas,"#alpha_{s} p_{T} dependence","L");

  c2b->cd(2);
  gPad->SetLogx();

  tdrDraw(gpsd,"Pz",kFullSquare,kBlue);
  tdrDraw(gpsd0,"Pz",kOpenTriangleDown,kRed+0);
  tdrDraw(gpsd1,"Pz",kFullTriangleDown,kRed+1);
  tdrDraw(gpsd2,"Pz",kOpenTriangleUp,kRed+2);
  tdrDraw(gpsd3,"Pz",kFullTriangleUp,kRed+3);

  c2b->SaveAs("pdf/drawZalpha_Run2_PSWeights_FSRvsPt.pdf");


  //TCanvas *c2 = tdrCanvas("c2",h2,4,11,kSquare);
  TCanvas *c2 = tdrDiCanvas("c2",h2,h2dw,4,11);

  c2->cd(1);
  gPad->SetLogx();
  fas->Draw("SAME");

  tdrDraw(gd,"Pz",kFullCircle,kGray+1); gd->SetMarkerSize(0.5);
  tdrDraw(gm,"Pz",kOpenCircle,kGray+1); gm->SetMarkerSize(0.5);
  tdrDraw(gd2,"Pz",kFullSquare);        gd2->SetMarkerSize(0.5);
  tdrDraw(gm2,"Pz",kOpenSquare);        gm2->SetMarkerSize(0.5);
  // Draw PS variations
  if (drawPS) {
    if (drawData) {
      tdrDraw(gps0,"Pz",kOpenDiamond,kRed+0);  gps0->SetMarkerSize(0.5);
      tdrDraw(gps1,"Pz",kOpenDiamond,kRed+1);  gps1->SetMarkerSize(0.5);
      tdrDraw(gps2,"Pz",kOpenDiamond,kRed+2);  gps2->SetMarkerSize(0.5);
      tdrDraw(gps3,"Pz",kOpenDiamond,kRed+3);  gps3->SetMarkerSize(0.5);
    }
    //else {
    //c2b->cd();
    //}
  }
  //c2->SetLogx();

  //TLegend *leg = tdrLeg(0.20,0.60,0.40,0.85);
  TLegend *leg2 = tdrLeg(0.50,0.65,0.70,0.90);
  //leg2->AddEntry(gd2,"Data (q equiv.)","P");
  //leg2->AddEntry(gm2,"MC (q equiv.)","P");
  leg2->AddEntry(gd2,"Data","P");
  leg2->AddEntry(gm2,"MC","P");
  leg2->AddEntry(fas,"#alpha_{s} p_{T} dependence","L");
  //leg2->AddEntry(gd,"Data (Z+jet mix)","P");
  //leg2->AddEntry(gm,"MC (Z+jet mix)","P");

  c2->cd(2);
  gPad->SetLogx();

  TGraphErrors *gd2d = new TGraphErrors(gd2->GetN());
  TGraphErrors *gm2d = new TGraphErrors(gm2->GetN());

  for (int i = 0; i != gd2d->GetN(); ++i) {

    double x = gm2->GetX()[i];
    double ex = gm2->GetEX()[i];
    double yref = gm2->GetY()[i];
    double eyref = gm2->GetEY()[i];
    gd2d->SetPoint(i, x, gd2->GetY()[i]-yref);
    gd2d->SetPointError(i, ex, sqrt(pow(gd2->GetEY()[i],2)+eyref*eyref));
    gm2d->SetPoint(i, x, gm2->GetY()[i]-yref);
    gm2d->SetPointError(i, ex, 0);//sqrt(pow(gm2->GetEY()[i],2)+eyref*eyref));
  }

  tdrDraw(gd2d,"Pz",kFullSquare);        gd2d->SetMarkerSize(0.5);
  tdrDraw(gm2d,"Pz",kOpenSquare);        gm2d->SetMarkerSize(0.5);

  c2->SaveAs("pdf/drawZalpha_Run2_FSRvsPt.pdf");


  TH1D *h3 = new TH1D("h3",";p_{T,Z} (GeV);Fit #chi^{2}/NDF",970,30,1000);
  h3->GetXaxis()->SetMoreLogLabels();
  h3->GetXaxis()->SetNoExponent();
  h3->SetMaximum(20);//100);
  h3->SetMinimum(0);

  TCanvas *c3 = tdrCanvas("c3",h3,4,11,kSquare);
  c3->SetLogx();
  
  tdrDraw(gchi2dt,"P",kFullCircle);
  tdrDraw(gchi2mc,"P",kOpenCircle);
  // Draw parton shower variations
  if (drawPS) {
    tdrDraw(gchi2ps0,"P",kOpenDiamond,kRed+0); gchi2ps0->SetMarkerSize(0.5);
    tdrDraw(gchi2ps1,"P",kOpenDiamond,kRed+1); gchi2ps1->SetMarkerSize(0.5);
    tdrDraw(gchi2ps2,"P",kOpenDiamond,kRed+2); gchi2ps2->SetMarkerSize(0.5);
    tdrDraw(gchi2ps3,"P",kOpenDiamond,kRed+3); gchi2ps3->SetMarkerSize(0.5);
  }

  TLegend *leg3 = tdrLeg(0.50,0.70,0.70,0.80);
  leg3->AddEntry(gchi2dt,"Data","P");
  leg3->AddEntry(gchi2mc,"MC","P");

  c3->SaveAs("pdf/drawZalpha_Run2_Chi2vsPt.pdf");


  TH1D *h4 = new TH1D("h4",";p_{T,Z} (GeV);Fit MPF (%)",970,30,1000);
  h4->GetXaxis()->SetMoreLogLabels();
  h4->GetXaxis()->SetNoExponent();
  h4->SetMaximum(+1.0-1e-4);
  h4->SetMinimum(-2.0+1e-4);

  TH1D *h4d = new TH1D("h4",";p_{T,Z} (GeV);Data/MC-1 (%)",970,30,1000);
  h4d->GetXaxis()->SetMoreLogLabels();
  h4d->GetXaxis()->SetNoExponent();
  h4d->SetMaximum(+0.5+1e-3);//+0.5
  h4d->SetMinimum(-0.5-1e-3);//-0.5

  //TCanvas *c4 = tdrCanvas("c4",h4,4,11,kSquare);
  //c4->SetLogx();
  TCanvas *c4 = tdrDiCanvas("c4",h4,h4d,4,11);
  
  c4->cd(1);
  gPad->SetLogx();
  tdrDraw(gmpfdt,"Pz",kFullCircle); gmpfdt->SetMarkerSize(0.5);
  tdrDraw(gmpfmc,"Pz",kOpenCircle); gmpfmc->SetMarkerSize(0.5);

  TLegend *leg4 = tdrLeg(0.50,0.70,0.70,0.80);
  leg4->AddEntry(gmpfdt,"Data","P");
  leg4->AddEntry(gmpfmc,"MC","P");

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  if (etamax==2.5)
    tex->DrawLatex(0.50,0.85,"|#eta| < 2.5, R = 0.4 PF+CHS");
  if (etamax==1.3)
    tex->DrawLatex(0.50,0.85,"|#eta| < 1.3, R = 0.4 PF+CHS");

  c4->cd(2);
  gPad->SetLogx();

  tdrDraw(gmpfr,"Pz",kFullCircle); gmpfr->SetMarkerSize(0.5);

  TLine *l = new TLine();
  l->SetLineStyle(kDotted);
  l->DrawLine(30,0,1000,0);

  c4->SaveAs("pdf/drawZalpha_Run2_MPFvsPt.pdf");

  // Next steps: add R_MPF_jet1, R_MPF_jetn, R_MPF_unclus to corrections
  // Figure out how to decouple unclustered correction and combined FSR fit
  // (or this is natural by using R_MPF_unclus?)
  // (R_MPF -> R_MPF + R_unclus/Rsoft vs R_MPF -> R_MPF + FSR*(1-Rsoft))
  // (where FSR = 1 - R_MPF_jet1 - R_MPF_jetn + UEcorr or from pTbal slope)
  // TO-DO: correct bias from muon JER (and electron/muon scales?)
  // TO-DO: figure out Fit MPF bias in MC
  // TO-DO: uncertainty estimates (FSR, Rsoft, unclust etc.)
  // TO-DO: extend model to alpha>=0.5 to get pT<60 GeV back in
} // drawZalpha

//void drawZalpha() {
//pair<double,double> drawZalphas(double ptref) {
FSR drawZalphas(double ptref) {
  //double ptref = 100.;
  //double ptref = 200.;
  //double ptref = 250.;

  setTDRStyle();

  //TFile *f = new TFile("rootfiles/jecdataABC.root","READ");
  //TFile *f = new TFile("rootfiles/jecdataD.root","READ");
  //Tfile *f = new TFile("rootfiles/jecdataABCD.root","READ");
  //TFile *f = new TFile("rootfiles/jme_bplusZ_Run2_v18.root","READ"); // Samimode
  //TFile *f = new TFile("rootfiles/jme_bplusZ_Run2_v20_noTTBar.root","READ");
  //TFile *f = new TFile("rootfiles/jme_bplusZ_Run2_v20.root","READ");
  //TFile *f = new TFile("rootfiles/jme_bplusZ_Run2_v23.root","READ");
  //TFile *f = new TFile("rootfiles/jme_bplusZ_Run2_v25.root","READ");
  //TFile *f = new TFile("rootfiles/jme_bplusZ_Run2_v26.root","READ");
  //TFile *f = new TFile("rootfiles/jme_bplusZ_Run2_Ver1.root","READ");
  //TFile *f = new TFile("rootfiles/jme_bplusZ_Run2_Ver2.root","READ");
  //TFile *f = new TFile("rootfiles/jme_bplusZ_Run2_Ver3.root","READ");
  //TFile *f = new TFile("rootfiles/jme_bplusZ_Run2_Ver4.root","READ");
  //TFile *f = new TFile("rootfiles/jme_bplusZ_Run2_Ver5.root","READ");
  TFile *f = new TFile("rootfiles/jme_bplusZ_Run2_Ver6.root","READ");
  assert(f && !f->IsZombie());

  /*
  string sdp = "data/eta00-13/ptchs_zlljet_a";
  string sdm = "data/eta00-13/mpfchs1_zlljet_a";
  string smp = "mc/eta00-13/ptchs_zlljet_a";
  string smm = "mc/eta00-13/mpfchs1_zlljet_a";
  string sst = "";
  string sgf = "";
  if (true) { // Samimode
    if (etamax==2.5) {
      sdp = "data/eta_00_25/ptchs_zmmjet_a";
      sdm = "data/eta_00_25/mpfchs_zmmjet_a";
      smp = "mc/eta_00_25/ptchs_zmmjet_a";
      smm = "mc/eta_00_25/mpfchs_zmmjet_a";
      sst = "mc/eta_00_25/statistics_zmmjet_a";
      sgf = "mc/eta_00_25/gpurity_zmmjet_a";
    }
    if (etamax==1.3) {
      //sdp = "data/eta_00_13/ptchs_zmmjet_a";
      sdp = "data/eta_00_13/rmpfjet1_zmmjet_a";//v25
      //sdm = "data/eta_00_13/mpfchs_zmmjet_a";
      sdm = "data/eta_00_13/rmpf_zmmjet_a";//v25
      //smp = "mc/eta_00_13/ptchs_zmmjet_a";
      smp = "mc/eta_00_13/rmpfjet1_zmmjet_a";//v25
      //smm = "mc/eta_00_13/mpfchs_zmmjet_a";
      smm = "mc/eta_00_13/rmpf_zmmjet_a";//v25
      sst = "mc/eta_00_13/statistics_zmmjet_a";
      sgf = "mc/eta_00_13/gpurity_zmmjet_a";
    }
  }
  const char *cdp = sdp.c_str();
  const char *cdm = sdm.c_str();
  const char *cmp = smp.c_str();
  const char *cmm = smm.c_str();
  const char *cst = sst.c_str();
  const char *cgf = sgf.c_str();
  */

  // alpha bins
  const int na = 10;
  //int va[na] = {100,80,60,50,40,30,25,20,15,10};
  const int va[na] = {10,15,20,25,30,40,50,60,80,100};
  // data type
  const int nd = 6;//2;
  //const string vd[nd] = {"data","mc"}; // "ratio"
  const string vd[nd] = {"data","mc","ps0","ps1","ps2","ps3"}; // "ratio"
  // histogram / graph type
  const int nh = 6;//5;//4;
  //const string vh[nh] = {"mpfchs","ptchs","statistics","gpurity"};
  const string vh[nh] = {"stat","mpf","pt","gf","unc","jn"};
  // shorter names to longer ones
  map<string,string> mstol;
  //mstol["mpf"] = "mpfchs";
  mstol["mpf"] = "rmpf";//v25
  //mstol["pt"] = "ptchs";
  mstol["pt"] = "rmpfjet1";//v25
  mstol["unc"] = "rmpfuncl";//v25
  mstol["jn"] = "rmpfjetn";//v25
  mstol["stat"] = "statistics";
  mstol["gf"] = "gpurity";
  mstol["ps0"] = "PSWeight0";
  mstol["ps1"] = "PSWeight1";
  mstol["ps2"] = "PSWeight2";
  mstol["ps3"] = "PSWeight3";
  // maps of histograms and graphs
  map<string, map<string, map<int, TH1D* > > > mh;
  map<string, map<string, map<int, TGraphErrors* > > > mg;
  map<string, map<string, TGraphErrors* > > mga;
  TGraphErrors *gd(0), *gm(0), *gdm(0), *gmm(0);
  TGraphErrors *gpps0(0), *gpps1(0), *gpps2(0), *gpps3(0);
  TGraphErrors *gmps0(0), *gmps1(0), *gmps2(0), *gmps3(0);

  // Retrieve all the relevant information from file
  string seta = (etamax==1.3 ? "00_13" : "00_25");
  const char *ceta = seta.c_str();
  int ipt(-1);

  for (int id = 0; id != nd; ++id) {
    for (int ih = 0; ih != nh; ++ih) {

      string sd = vd[id]; const char *cd = sd.c_str();
      string sh = vh[ih]; const char *ch = mstol[sh].c_str();

      TGraphErrors *g(0); TH1D *h(0);
      TGraphErrors *ga = new TGraphErrors(0);
      for (int ia = 0; ia != na; ++ia) {
	
	// Construct histogram/graph name
	int a = va[ia];
	string shname = Form("%s/eta_%s/%s_zmmjet_a%d",cd,ceta,ch,a);
	if (sd=="ps0"||sd=="ps1"||sd=="ps2"||sd=="ps3") {
	  const char *cd = mstol[sd].c_str();
	  shname = Form("mc/eta_%s/%s/%s_zmmjet_a%d",ceta,cd,ch,a);	  
	}

	// histograms and graphs for data and mc
	if (sh=="stat") {
	  h = (TH1D*)f->Get(shname.c_str());
	  if (!h) cout << "Missing " << shname << endl << flush;
	  assert(h);
	  
	  mh[sd][sh][a] = h;
	} // histos

	if (sh=="mpf" || sh=="pt" || sh=="unc" || sh=="jn" ||
	    (sd=="mc" && sh=="gf")) {
	  TGraphErrors *g = (TGraphErrors*)f->Get(shname.c_str());
	  if (!g) cout << "Missing " << shname << endl << flush;
	  assert(g);
	  
	  mg[sd][sh][a] = g;
	 
	  // store purity for alpha<0.3 in a histogram for easy retrieval
	  if (sd=="mc" && sh=="gf" && a==30 && g && !_hgf) {
	    
	    TH1D *href = mh["data"]["stat"][30]; assert(href);
	    _hgf = (TH1D*)href->Clone("_hgf");
	    for (int i = 0; i != g->GetN(); ++i) {
	      int ipt = _hgf->FindBin(g->GetX()[i]);
	      _hgf->SetBinContent(i, g->GetY()[i]);
	      _hgf->SetBinError(i, g->GetEY()[i]);
	    } // for i
	  }

	  // map points vs alpha for given pt
	  if (ipt==-1 && g) {
	    double dpt(1000); //int ipt(-1);
	    for (int i = 0; i != g->GetN(); ++i) {
	      if (fabs(ptref-g->GetX()[i])<dpt) {
		dpt = fabs(ptref-g->GetX()[i]);
		ipt = i;
	      }
	    } // for i
	  }

	  // Check that pT bin is always the same for ipt
	  assert(g->GetN()>=ipt);
	  TH1D *hdst = mh["data"]["stat"][10]; //assert(hdst);
	  if (hdst) {
	    int ibin = hdst->FindBin(ptref);
	    assert(hdst->FindBin(g->GetX()[ipt])==ibin);
	  }

	  tools::SetPoint(ga, ia, 0.01*a, g->GetY()[ipt], 0, g->GetEY()[ipt]);
	} // graphs
      } // for ia

      if (sd=="data" && sh=="pt")  gd = ga;
      if (sd=="data" && sh=="mpf") gdm = ga;
      if (sd=="mc"   && sh=="pt")  gm = ga;
      if (sd=="mc"   && sh=="mpf") gmm = ga;
      //
      if (sd=="ps0"   && sh=="pt")  gpps0 = ga;
      if (sd=="ps1"   && sh=="pt")  gpps1 = ga;
      if (sd=="ps2"   && sh=="pt")  gpps2 = ga;
      if (sd=="ps3"   && sh=="pt")  gpps3 = ga;
      if (sd=="ps0"   && sh=="mpf")  gmps0 = ga;
      if (sd=="ps1"   && sh=="mpf")  gmps1 = ga;
      if (sd=="ps2"   && sh=="mpf")  gmps2 = ga;
      if (sd=="ps3"   && sh=="mpf")  gmps3 = ga;

      mga[sd][sh] = ga;
    } // for ih
  } // for id
  assert(ipt>=0);
  assert(gd);
  assert(gm);
  assert(gdm);
  assert(gmm);
  //
  assert(gpps0);
  assert(gpps1);
  assert(gpps2);
  assert(gpps3);
  assert(gmps0);
  assert(gmps1);
  assert(gmps2);
  assert(gmps3);

  assert(mg["data"]["pt"][30]);
  assert(mg["mc"]["pt"][30]);
  double pt = mg["data"]["pt"][30]->GetX()[ipt];
  double pt_mc = mg["mc"]["pt"][30]->GetX()[ipt];
  TH1D *hdst = mh["data"]["stat"][30]; assert(hdst);
  int ibin = hdst->FindBin(pt);
  assert(hdst->FindBin(pt_mc)==ibin);
  double minpt = hdst->GetBinLowEdge(ibin);
  double maxpt = hdst->GetBinLowEdge(ibin+1);

  if ((correctUnclustMPF || correctGluonMPF) && true) {

    bool a = correctGluonMPF;
    bool b = correctUnclustMPF;

    fCMPF(gdm, pt, mga["data"]["pt"],mga["data"]["jn"],mga["data"]["unc"],a,b);
    fCMPF(gmm, pt, mga["mc"]["pt"],  mga["mc"]["jn"],  mga["mc"]["unc"],  a,b);

    // gmps needs extra code to patch
    fCMPF(gmps0, pt, mga["ps0"]["pt"],mga["ps0"]["jn"],mga["ps0"]["unc"],a,b);
    fCMPF(gmps1, pt, mga["ps1"]["pt"],mga["ps1"]["jn"],mga["ps1"]["unc"],a,b);
    fCMPF(gmps2, pt, mga["ps2"]["pt"],mga["ps2"]["jn"],mga["ps2"]["unc"],a,b);
    fCMPF(gmps3, pt, mga["ps3"]["pt"],mga["ps3"]["jn"],mga["ps3"]["unc"],a,b);
  }

  // Scale to %
  // and shift all by small factor to get MC MPF at 1
  // and remove alpha<15/pT,Z
  double sf = 1.000;//0.994;
  for (int i = 0; i != gd->GetN(); ++i) {
    tools::SetPoint(gd, i, gd->GetX()[i], 100.*(sf*gd->GetY()[i]-1),
		    0, gd->GetEY()[i]*100);
    tools::SetPoint(gm, i, gm->GetX()[i], 100.*(sf*gm->GetY()[i]-1),
		    0, gm->GetEY()[i]*100);
    tools::SetPoint(gdm, i, gdm->GetX()[i], 100.*(sf*gdm->GetY()[i]-1),
		    0, gdm->GetEY()[i]*100);
    tools::SetPoint(gmm, i, gmm->GetX()[i], 100.*(sf*gmm->GetY()[i]-1),
		    0, gmm->GetEY()[i]*100);
    //
    tools::SetPoint(gpps0, i, gpps0->GetX()[i], 100.*(sf*gpps0->GetY()[i]-1),
		    0, gpps0->GetEY()[i]*100);
    tools::SetPoint(gpps1, i, gpps1->GetX()[i], 100.*(sf*gpps1->GetY()[i]-1),
		    0, gpps1->GetEY()[i]*100);
    tools::SetPoint(gpps2, i, gpps2->GetX()[i], 100.*(sf*gpps2->GetY()[i]-1),
		    0, gpps2->GetEY()[i]*100);
    tools::SetPoint(gpps3, i, gpps3->GetX()[i], 100.*(sf*gpps3->GetY()[i]-1),
		    0, gpps3->GetEY()[i]*100);
    tools::SetPoint(gmps0, i, gmps0->GetX()[i], 100.*(sf*gmps0->GetY()[i]-1),
		    0, gmps0->GetEY()[i]*100);
    tools::SetPoint(gmps1, i, gmps1->GetX()[i], 100.*(sf*gmps1->GetY()[i]-1),
		    0, gmps1->GetEY()[i]*100);
    tools::SetPoint(gmps2, i, gmps2->GetX()[i], 100.*(sf*gmps2->GetY()[i]-1),
		    0, gmps2->GetEY()[i]*100);
    tools::SetPoint(gmps3, i, gmps3->GetX()[i], 100.*(sf*gmps3->GetY()[i]-1),
		    0, gmps3->GetEY()[i]*100);
  }

  // Remove points in ineffective range
  for (int i = gd->GetN()-1; i!=-1; --i) {

    //if (gd->GetX()[i]<15./pt) {
    //if (gd->GetX()[i]<15./maxpt) {
    //if (gd->GetX()[i]<20./maxpt) {
    if (gd->GetX()[i]<ptthr/minpt) {
      gd->RemovePoint(i);
      gdm->RemovePoint(i);
    }
  }
  for (int i = gm->GetN()-1; i!=-1; --i) {
    if (gm->GetX()[i]<ptthr/minpt) {
      gm->RemovePoint(i);
      gmm->RemovePoint(i);
      //
      gpps0->RemovePoint(i);
      gpps1->RemovePoint(i);
      gpps2->RemovePoint(i);
      gpps3->RemovePoint(i);
      gmps0->RemovePoint(i);
      gmps1->RemovePoint(i);
      gmps2->RemovePoint(i);
      gmps3->RemovePoint(i);
    }
  }

  double ptmax = 1200.*2.;
  //TH1D *hup0 = new TH1D("hup0",";#alpha_{max};#LTp_{T,jet}#GT / p_{T,Z+#LTUE#GT}^{ } - 1 (%)",10,0,1.05);
  TH1D *hup0 = new TH1D("hup0",";#alpha_{max};#LTp_{T,jet}#GT / p_{T,Z}^{ } - 1 (%)",10,0,1.05);
  hup0->SetMinimum(-17);
  hup0->SetMaximum(+5);
  hup0->GetXaxis()->SetMoreLogLabels();
  hup0->GetXaxis()->SetNoExponent();

  TH1D *hdw0 = new TH1D("hdw0",";#alpha_{max};MC/ref-1 (%)",10,0,1.05);
  hdw0->SetMinimum(-1.0);
  hdw0->SetMaximum(+2.0);
  hdw0->GetXaxis()->SetMoreLogLabels();
  hdw0->GetXaxis()->SetNoExponent();

  TH1D *hup = new TH1D("hup",";#alpha_{max};#LTp_{T,jet}#GT / p_{T,Z}^{ } - 1 (%)",
		       10,0,1.05);//0.52);//0.34);
  hup->SetMinimum(-17);//-22);
  hup->SetMaximum(+5);

  TH1D *hdw = new TH1D("hdw",";#alpha_{max};Data/MC-1 (%)",10,0,1.05);//0.52);//0.34);
  hdw->SetMinimum(-1.0);//-0.9);//-0.5);//-1.0);//-3.5);
  hdw->SetMaximum(+2.0);//+0.7);//+1.5);//+1.0);
  hdw->GetXaxis()->SetMoreLogLabels();
  hdw->GetXaxis()->SetNoExponent();

  //lumi_13TeV = "Run2016BCDEFGH 36.5 fb^{-1}";
  //lumi_13TeV = "Run2016 36.5 fb^{-1}";
  //lumi_13TeV = "Run2018ABCD 59.9 fb^{-1}";
  lumi_13TeV = "Run2 136.5 fb^{-1}";
  extraText = "Private Work";
  //TCanvas *c0 = tdrCanvas("c0",hup0,4,11,kSquare);
  TCanvas *c0 = tdrDiCanvas("c0",hup0,hdw0,4,11);
  hup0->GetYaxis()->SetTitleOffset(1.00);

  TGraphErrors *gpps = new TGraphErrors(gpps0->GetN());
  TGraphErrors *gppsd = new TGraphErrors(gpps0->GetN());
  TGraphErrors *gppsd0 = new TGraphErrors(gpps0->GetN());
  TGraphErrors *gppsd1 = new TGraphErrors(gpps0->GetN());
  TGraphErrors *gppsd2 = new TGraphErrors(gpps0->GetN());
  TGraphErrors *gppsd3 = new TGraphErrors(gpps0->GetN());
  TGraphErrors *gmps = new TGraphErrors(gmps0->GetN());
  TGraphErrors *gmpsd = new TGraphErrors(gmps0->GetN());
  TGraphErrors *gmpsd0 = new TGraphErrors(gmps0->GetN());
  TGraphErrors *gmpsd1 = new TGraphErrors(gmps0->GetN());
  TGraphErrors *gmpsd2 = new TGraphErrors(gmps0->GetN());
  TGraphErrors *gmpsd3 = new TGraphErrors(gmps0->GetN());

  for (int i = 0; i != gpps->GetN(); ++i) {
    double x = gpps0->GetX()[i];
    double ex = gpps0->GetEX()[i];
    double y = 0.5*(gpps0->GetY()[i] + gpps2->GetY()[i]);
    double ey = 0.5*(gpps0->GetEY()[i] + gpps2->GetEY()[i]);
    double ke = 0.1; // approximate for highly correlated statistics
    gpps->SetPoint(i, x, y);
    gpps->SetPointError(i, ex, ey);
    gppsd->SetPoint(i, x, gpps->GetY()[i]-y);
    gppsd->SetPointError(i, ex, 0);
    gppsd0->SetPoint(i, x, gpps0->GetY()[i]-y);
    gppsd0->SetPointError(i, ex, sqrt(pow(gpps0->GetEY()[i],2)+ey*ey)*ke);
    gppsd1->SetPoint(i, x, gpps1->GetY()[i]-y);
    gppsd1->SetPointError(i, ex, sqrt(pow(gpps1->GetEY()[i],2)+ey*ey)*ke);
    gppsd2->SetPoint(i, x, gpps2->GetY()[i]-y);
    gppsd2->SetPointError(i, ex, sqrt(pow(gpps2->GetEY()[i],2)+ey*ey)*ke);
    gppsd3->SetPoint(i, x, gpps3->GetY()[i]-y);
    gppsd3->SetPointError(i, ex, sqrt(pow(gpps3->GetEY()[i],2)+ey*ey)*ke);
    //
    double z = 0.5*(gmps0->GetY()[i] + gmps2->GetY()[i]);
    double ez = 0.5*(gmps0->GetEY()[i] + gmps2->GetEY()[i]);
    gmps->SetPoint(i, x, z);
    gmps->SetPointError(i, ex, ez);
    gmpsd->SetPoint(i, x, gmps->GetY()[i]-z);
    gmpsd->SetPointError(i, ex, 0);
    gmpsd0->SetPoint(i, x, gmps0->GetY()[i]-z);
    gmpsd0->SetPointError(i, ex, sqrt(pow(gmps0->GetEY()[i],2)+ez*ez)*ke);
    gmpsd1->SetPoint(i, x, gmps1->GetY()[i]-z);
    gmpsd1->SetPointError(i, ex, sqrt(pow(gmps1->GetEY()[i],2)+ez*ez)*ke);
    gmpsd2->SetPoint(i, x, gmps2->GetY()[i]-z);
    gmpsd2->SetPointError(i, ex, sqrt(pow(gmps2->GetEY()[i],2)+ez*ez)*ke);
    gmpsd3->SetPoint(i, x, gmps3->GetY()[i]-z);
    gmpsd3->SetPointError(i, ex, sqrt(pow(gmps3->GetEY()[i],2)+ez*ez)*ke);
  }

  /*
  c0->cd(1);

  tdrDraw(gpps,"Pz",kFullSquare,kRed);
  tdrDraw(gpps0,"Pz",kOpenTriangleDown,kRed);
  tdrDraw(gpps1,"Pz",kFullTriangleDown,kRed+1);
  tdrDraw(gpps2,"Pz",kOpenTriangleUp,kRed+2);
  tdrDraw(gpps3,"Pz",kFullTriangleUp,kRed+3);
  tdrDraw(gmps,"Pz",kFullSquare,kBlue);
  tdrDraw(gmps0,"Pz",kOpenTriangleDown,kBlue); 
  tdrDraw(gmps1,"Pz",kFullTriangleDown,kBlue+1);
  tdrDraw(gmps2,"Pz",kOpenTriangleUp,kBlue+2);
  tdrDraw(gmps3,"Pz",kFullTriangleUp,kBlue+3);

  gpps->SetMarkerSize(0.5);
  gpps0->SetMarkerSize(0.5);
  gpps1->SetMarkerSize(0.5);
  gpps2->SetMarkerSize(0.5);
  gpps3->SetMarkerSize(0.5);
  gmps->SetMarkerSize(0.4);
  gmps0->SetMarkerSize(0.4);
  gmps1->SetMarkerSize(0.4);
  gmps2->SetMarkerSize(0.4);
  gmps3->SetMarkerSize(0.4);
  */

  c0->cd(2);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(0,0,1,0);

  //tdrDraw(gppsd,"Pz",kFullSquare,kRed);
  tdrDraw(gppsd0,"Pz",kOpenTriangleDown,kRed);
  tdrDraw(gppsd1,"Pz",kFullTriangleDown,kRed+1);
  tdrDraw(gppsd2,"Pz",kOpenTriangleUp,kRed+2);
  tdrDraw(gppsd3,"Pz",kFullTriangleUp,kRed+3);
  //tdrDraw(gmpsd,"Pz",kFullSquare,kBlue);
  tdrDraw(gmpsd0,"Pz",kOpenTriangleDown,kBlue); 
  tdrDraw(gmpsd1,"Pz",kFullTriangleDown,kBlue+1);
  tdrDraw(gmpsd2,"Pz",kOpenTriangleUp,kBlue+2);
  tdrDraw(gmpsd3,"Pz",kFullTriangleUp,kBlue+3);

  gppsd->SetMarkerSize(0.5);
  gppsd0->SetMarkerSize(0.5);
  gppsd1->SetMarkerSize(0.5);
  gppsd2->SetMarkerSize(0.5);
  gppsd3->SetMarkerSize(0.5);
  gmpsd->SetMarkerSize(0.4);
  gmpsd0->SetMarkerSize(0.4);
  gmpsd1->SetMarkerSize(0.4);
  gmpsd2->SetMarkerSize(0.4);
  gmpsd3->SetMarkerSize(0.4);


  TCanvas *c1 = tdrDiCanvas("c1",hup,hdw,4,11);
  hup->GetYaxis()->SetTitleOffset(1.00);

  c1->cd(1);
  tdrDraw(gdm,"Pz",kFullSquare,kBlue); gdm->SetMarkerSize(0.5);
  tdrDraw(gmm,"Pz",kOpenSquare,kBlue); gmm->SetMarkerSize(0.5);
  tdrDraw(gd,"Pz",kFullCircle,kRed);   gd->SetMarkerSize(0.5);
  tdrDraw(gm,"Pz",kOpenCircle,kRed);   gm->SetMarkerSize(0.5);

  // Draw also PS
  if (drawPS) {
    if (drawData) {
      tdrDraw(gpps0,"Pz",kOpenDiamond,kRed); gpps0->SetMarkerSize(0.5);
      tdrDraw(gpps1,"Pz",kOpenDiamond,kRed+1); gpps1->SetMarkerSize(0.5);
      tdrDraw(gpps2,"Pz",kOpenDiamond,kRed+2); gpps2->SetMarkerSize(0.5);
      tdrDraw(gpps3,"Pz",kOpenDiamond,kRed+3); gpps3->SetMarkerSize(0.5);
      tdrDraw(gmps0,"Pz",kOpenDiamond,kBlue); gmps0->SetMarkerSize(0.5);
      tdrDraw(gmps1,"Pz",kOpenDiamond,kBlue+1); gmps1->SetMarkerSize(0.5);
      tdrDraw(gmps2,"Pz",kOpenDiamond,kBlue+2); gmps2->SetMarkerSize(0.5);
      tdrDraw(gmps3,"Pz",kOpenDiamond,kBlue+3); gmps3->SetMarkerSize(0.5);
    }
    //else {
      //c0->cd();
      //c1->cd(1);
    //}
    
  }

  // Should correct MPF for lower recoil response:
  // type-I not applied for pT,corr<15 GeV, response between 0.65-0.9
  // it is however applied for 2nd jet, which reduces impact (mostly 3rd on)
  // Would need:
  // 1) vector projection of leading jet
  // 2) vector projection of subleading jet
  // 3) vector projection of unclustered energy (-MET - pT,lead - pT,jetsum)
  // or maybe just the last one?
  // Expectation is small slope down for MPF, as observed. However, different
  // for data and MC mainly because of L2Res not applied to unclustered pT?

  // Maybe should also correct pT balance for know UE effects?
  // pTgen ~ pT,Z + pT,UE, where pT,UE at particle level for *dijets*
  // This is because R_MC absorbs large R_jet-R_UE difference into a small
  // change of R_jet. So correcting pT,jet (>>pT,UE) by slightly smaller R_jet
  // accounts also for the low response of R_UE.
  // For Z+jet, L1FastJet puts back in the dijet UE at reco level, then
  // corrects the whole jet by the smaller R_jet to bring dijet UE to gen level

  TF1 *fdm = new TF1("fdm","[0]+[1]*x",amin,amax);
  fdm->FixParameter(1,0); // no slope for MPF
  fdm->SetParameters(0, 0);
  fdm->SetLineColor(kBlue-0);
  fdm->SetLineStyle(kSolid);
  fdm->SetLineWidth(2);
  gdm->Fit(fdm,"QRN");
  fdm->Draw("SAME");

  TF1 *fmm = new TF1("fmm","[0]+[1]*x",0,amax);
  fmm->FixParameter(1,0); // no slope for MPF
  fmm->SetParameters(0, 0);
  gmm->Fit(fmm,"QRN");
  fmm->SetLineColor(kBlue);//-9);
  fmm->SetLineStyle(kDashed);
  fmm->SetLineWidth(2);
  fmm->Draw("SAME");

  // Scale data by MPF data/MC ratio
  if (false) {
    double  k = (fmm->GetParameter(0)/100.+1) / (fdm->GetParameter(0)/100.+1);
    for (int i = 0; i != gd->GetN(); ++i) {
      gd->SetPoint(i, gd->GetX()[i], ((gd->GetY()[i]/100.+1)*k-1)*100.);
      gdm->SetPoint(i, gdm->GetX()[i], ((gdm->GetY()[i]/100.+1)*k-1)*100.);
    }
    gdm->Fit(fdm,"QRN");
  }

  TF1 *fm = new TF1("fm","[0]+[1]*x",amin,amax);
  fm->SetParameters(0, _fas->Eval(pt));
  //TF1 *fm = new TF1("fm","[0]+[1]*(x-[0]/[2])",amin,amax);
  //fm->FixParameter(2,pt);
  fm->SetLineColor(kRed);//-9);
  fm->SetLineStyle(kDashed);
  //double ue = TMath::Pi()*0.4*0.4*2.00;//1.52*log(13)/log(8); // Z+jet UE
  //
  //double uepred = TMath::Pi()*0.4*0.4*(1.30/0.6 + // Z+jet UE over R_UE
    //			       0.7/0.9); // dijet UE excess over R_jet
  double uepred = TMath::Pi()*0.4*0.4*2.0/0.6; // dijet UE at ptcl (*1/R_UE)
  //
  //double ue = TMath::Pi()*0.4*0.4*2.00*log(13)/log(8); // dijet UE
  //double ue = TMath::Pi()*0.4*0.4*2.00; // dijet UE??
  //fm->FixParameter(0, ue/pt * 100.);
  //fm->FixParameter(0, -0.5); // PATCH
  //fm->SetParameter(0, ue/pt * 100.);
  gm->Fit(fm,"QRN");
  /*
  double uemc = fm->GetParameter(0)*pt/100.;
  double duemc = fm->GetParError(0)*pt/100.;
  //gm->SetPoint(gm->GetN(), 0, ue/pt * 100.); // UE/pT,Z-1
  gm->SetPoint(gm->GetN(), uemc/pt, uemc/pt * 100.); // UE/pT,Z
  //gm->SetPoint(gm->GetN(), 0, -0.5); // PATCH
  //gmm->SetPoint(gmm->GetN(), 0, 1.50*ue/pt * 100.); // R_UE*UE/pT,Z-1
  fm->Draw("SAME");
  */

  TF1 *fd = new TF1("fd","[0]+[1]*x",amin,amax);
  fd->SetParameters(0, _fas->Eval(pt));
  //TF1 *fd = new TF1("fd","[0]+[1]*(x-[0]/[2])",amin,amax);
  //fd->FixParameter(2,pt);
  //fd->FixParameter(0, ue/pt * 100.);
  //fd->FixParameter(0, -0.5); // PATCH
  fd->SetLineColor(kRed);//-9);
  //gd->SetPoint(gd->GetN(),0,fd->GetParameter(0)); // UE/pT,Z-1 (before fit)
  fd->SetLineStyle(kSolid);
  gd->Fit(fd,"QRN");
  /*
  double uedt = fd->GetParameter(0)*pt/100.;
  double duedt = fd->GetParError(0)*pt/100.;
  //gd->SetPoint(gd->GetN(),0,fd->GetParameter(0)); // UE/pT,Z-1 (after fit)
  gd->SetPoint(gd->GetN(), uedt/pt, uedt/pt * 100.);
  //gd->SetPoint(gd->GetN(),0,-0.5); // PATCH
  fd->Draw("SAME");
  */

  // Correct MPF for unclustered energy at pT,corr<15 GeV by dPt*(1-Rsoft)
  // (Combined FSR+ISR fit quite sensitive to both Rsoft and data/MC SF)
  // (Need to also look at chi2/NDF and other proxies)
  if (correctUnclustMPF && true) { // most code moved earlier
    gdm->Fit(fdm,"QRN");
    gmm->Fit(fmm,"QRN");
  }
  if (correctUnclustMPF && false) { // old method, not using unclustered pT

    int i3(-1);
    for (int i = 0; i != gdm->GetN() && i3<0; ++i) {
      if (fabs(gdm->GetX()[i]-0.3)<0.01) i3 = i;
    }
    //assert(i3>=0 || pt<50);
    assert(i3>=0 || pt<60);
    if (i3<0) i3=0;

    for (int i = 0; i != gdm->GetN(); ++i) {

      // Fit MPF at low pT very sensitive to Rsoft (and data/MC)
      double Rsoft = 0.7;//0.5;
      double Rsoftm = Rsoft;
      double Rsoftd = Rsoft*1.05;

      double fg = _hgf->GetBinContent(_hgf->FindBin(pt));
      double CX = fg*2. + (1-fg)*4./3.; 
      double ck = pow(CX / (4./3.), 2);

      double a = gmm->GetX()[i];
      double k = (a<0.5 ? 1 : 0.5); // suppression on unclustered pT
      //double pttype1 = 15.;
      double dptd = -k*fd->GetParameter(1)*pttype1/pt;
      if (pt < ptpatch && patchLowPtUnclus) {
	// Would need to correct this for Unclust and UE:
	//dptd = k*(gd->GetY()[i3]-gdm->GetY()[i3])/0.3*ptthr/pt;
	//dptd = -k * _fas->Eval(pt) * ck * ptthr/pt;
	const double alpharef = 0.3;
	dptd = -k * _fas->Eval(pt)/alpharef * ck * ptthr/pt;
      }
      gdm->SetPoint(i, gdm->GetX()[i], gdm->GetY()[i] + dptd*(1-Rsoftd));

      double dptm = -k*fm->GetParameter(1)*pttype1/pt;
      if (pt < ptpatch && patchLowPtUnclus) {
	// Would need to correct this for Unclust and UE:
	//dptm = k*(gm->GetY()[i3]-gmm->GetY()[i3])/0.3*ptthr/pt;
	//dptm = -k * _fas->Eval(pt) * ck * pttype1/pt;
	const double alpharef = 0.3;
	dptm = -k * _fas->Eval(pt)/alpharef * ck * pttype1/pt;
      }
      gmm->SetPoint(i, gmm->GetX()[i], gmm->GetY()[i] + dptm*(1-Rsoftm));
      // use same dpt for PS for 1st iteration
      gmps->SetPoint(i, gmps->GetX()[i], gmps->GetY()[i] + dptm*(1-Rsoftm));
      gmps0->SetPoint(i, gmps0->GetX()[i], gmps0->GetY()[i] + dptm*(1-Rsoftm));
      gmps1->SetPoint(i, gmps1->GetX()[i], gmps1->GetY()[i] + dptm*(1-Rsoftm));
      gmps2->SetPoint(i, gmps2->GetX()[i], gmps2->GetY()[i] + dptm*(1-Rsoftm));
      gmps3->SetPoint(i, gmps3->GetX()[i], gmps3->GetY()[i] + dptm*(1-Rsoftm));
    }
    gdm->Fit(fdm,"QRN");
    gmm->Fit(fmm,"QRN");
  } // correctUnclustMPF
  
  // Correct MPF for gluon jet (recoil) vs quark jet (lead jet) response
  if (correctGluonMPF && false) {

    // Gluon fraction for leading jet vs recoil
    //double gfraclead = 0.20;
    //double gfracrecoil = 0.80;
    //double dgf = gfracrecoil-gfraclead;
    // Quark-rich response (at pt) over gluon-rich response (at 15-pt*alpha)
    // Effective factor includes gluon fraction (20% lead, 805 for recoil?)
    double qovg = 1.03;

    for (int i = 0; i != gdm->GetN(); ++i) {

      //double ptthr = 15.; // pttype1
      double dptd = -fd->GetParameter(1)*(gd->GetX()[i]-pttype1/pt);
      gdm->SetPoint(i, gdm->GetX()[i], gdm->GetY()[i] + dptd*(qovg-1));

      double dptm = -fm->GetParameter(1)*(gm->GetX()[i]-pttype1/pt);
      gmm->SetPoint(i, gmm->GetX()[i], gmm->GetY()[i] + dptm*(qovg-1));
      gmps->SetPoint(i, gmps->GetX()[i], gmps->GetY()[i] + dptm*(qovg-1));
      gmps0->SetPoint(i, gmps0->GetX()[i], gmps0->GetY()[i] + dptm*(qovg-1));
      gmps1->SetPoint(i, gmps1->GetX()[i], gmps1->GetY()[i] + dptm*(qovg-1));
      gmps2->SetPoint(i, gmps2->GetX()[i], gmps2->GetY()[i] + dptm*(qovg-1));
      gmps3->SetPoint(i, gmps3->GetX()[i], gmps3->GetY()[i] + dptm*(qovg-1));
    }
    gdm->Fit(fdm,"QRN");
    gmm->Fit(fmm,"QRN");
  } // correctGluonMPF

  // CorrectPtBal for UE
  if (correctPtBalUE) {

    // UE at gen and FullSim levels from Minsuk (2017DE CP5)
    double ue_zg = 2.75;
    double ue_zr = 1.8; // Z+jet RECO CP5
    double ue_zd = 1.4; // Z+jet DATA (M1)
    double ue_dg = 3.4; // dijet GEN (tune?)
    double ue_dr = 2.2; // dijet RECO (tune?);
    double Rjet = 0.95; // jet response estimate
    double ajet = TMath::Pi()*0.4*0.4;

    double qovg = 1.03;

    for (int i = 0; i != gd->GetN(); ++i) {

      //double dptd = ue_zd*ajet/Rjet/pt*100.;
      double dptd = ue_dg*ajet/pt*100.;
      //double dptd = ue_dr/Rjet*ajet/pt*100.;
      gd->SetPoint(i, gd->GetX()[i], gd->GetY()[i] - dptd);

      //double dptm = ue_zr*ajet/Rjet/pt*100.;
      double dptm = ue_dg*ajet/pt*100.;
      //double dptm = ue_dr/Rjet*ajet/pt*100.;
      gm->SetPoint(i, gm->GetX()[i], gm->GetY()[i] - dptm);
      gpps->SetPoint(i, gpps->GetX()[i], gpps->GetY()[i] - dptm);
      gpps0->SetPoint(i, gpps0->GetX()[i], gpps0->GetY()[i] - dptm);
      gpps1->SetPoint(i, gpps1->GetX()[i], gpps1->GetY()[i] - dptm);
      gpps2->SetPoint(i, gpps2->GetX()[i], gpps2->GetY()[i] - dptm);
      gpps3->SetPoint(i, gpps3->GetX()[i], gpps3->GetY()[i] - dptm);
    }
    gd->Fit(fd,"QRN");
    gm->Fit(fm,"QRN");
    //hup->GetYaxis()->SetTitle("#LTp_{T,jet}#GT / (p_{T,Z}^{ } + "
    //		      "#LTp_{T,UE}#GT) - 1 (%)");
    hup->GetYaxis()->SetTitle("#LTp_{T,jet}#GT / p_{T,Z+#LTUE#GT}^{ } - 1 (%)");
  } // correctPtBalUE

  // Correct quark-rich (80%) Z+jet response to nominal
  if (correctQuarkResp) {

    // ud and g responses from Mikael
    double Rq = 1.005;
    double Rg = 0.990;
    double Rc = 1.000;
    double Rs = 0.990;
    double Rz = 0.6*Rq+0.1*Rs+0.2*Rg+0.1*Rc;
    double Rd = 0.3*Rq+0.6*Rg+0.1*Rc;
    
    double qovg = 1.03;

    for (int i = 0; i != gd->GetN(); ++i) {

      double dptd = (Rz-Rd)*100;
      gd->SetPoint(i, gd->GetX()[i], gd->GetY()[i] - dptd);
      gdm->SetPoint(i, gdm->GetX()[i], gdm->GetY()[i] - dptd);

      double dptm = (Rz-Rd)*100;
      gm->SetPoint(i, gm->GetX()[i], gm->GetY()[i] - dptm);
      gmm->SetPoint(i, gmm->GetX()[i], gmm->GetY()[i] - dptm);
      //
      gpps->SetPoint(i, gpps->GetX()[i], gpps->GetY()[i] - dptm);
      gpps0->SetPoint(i, gpps0->GetX()[i], gpps0->GetY()[i] - dptm);
      gpps1->SetPoint(i, gpps1->GetX()[i], gpps1->GetY()[i] - dptm);
      gpps2->SetPoint(i, gpps2->GetX()[i], gpps2->GetY()[i] - dptm);
      gpps3->SetPoint(i, gpps3->GetX()[i], gpps3->GetY()[i] - dptm);
      gmps->SetPoint(i, gmps->GetX()[i], gmps->GetY()[i] - dptm);
      gmps0->SetPoint(i, gmps0->GetX()[i], gmps0->GetY()[i] - dptm);
      gmps1->SetPoint(i, gmps1->GetX()[i], gmps1->GetY()[i] - dptm);
      gmps2->SetPoint(i, gmps2->GetX()[i], gmps2->GetY()[i] - dptm);
      gmps3->SetPoint(i, gmps3->GetX()[i], gmps3->GetY()[i] - dptm);
    }
    gd->Fit(fd,"QRN");
    gdm->Fit(fdm,"QRN");
    gm->Fit(fm,"QRN");
    gmm->Fit(fmm,"QRN");
  } // correctQuarkResp

  // Copy remaining ones for combined MPF and pT balance fit
  TGraphErrors *mgd = new TGraphErrors(gd->GetN()+gdm->GetN());
  TGraphErrors *mgm = new TGraphErrors(gm->GetN()+gmm->GetN());
  TGraphErrors *mgps = new TGraphErrors(gpps->GetN()+gmps->GetN());
  TGraphErrors *mgps0 = new TGraphErrors(gpps0->GetN()+gmps0->GetN());
  TGraphErrors *mgps1 = new TGraphErrors(gpps1->GetN()+gmps1->GetN());
  TGraphErrors *mgps2 = new TGraphErrors(gpps2->GetN()+gmps2->GetN());
  TGraphErrors *mgps3 = new TGraphErrors(gpps3->GetN()+gmps3->GetN());
  for (int i = gdm->GetN()-1; i != -1; --i) {
    mgd->SetPoint(gdm->GetN()-1-i, -gdm->GetX()[i], gdm->GetY()[i]);
    mgd->SetPointError(gdm->GetN()-1-i, gdm->GetEX()[i], gdm->GetEY()[i]);
    mgm->SetPoint(gmm->GetN()-1-i, -gmm->GetX()[i], gmm->GetY()[i]);
    mgm->SetPointError(gmm->GetN()-1-i, gmm->GetEX()[i], gmm->GetEY()[i]);
    //
    mgps->SetPoint(gmps->GetN()-1-i, -gmps->GetX()[i], gmps->GetY()[i]);
    mgps->SetPointError(gmps->GetN()-1-i, gmps->GetEX()[i], gmps->GetEY()[i]);
    mgps0->SetPoint(gmps0->GetN()-1-i, -gmps0->GetX()[i], gmps0->GetY()[i]);
    mgps0->SetPointError(gmps0->GetN()-1-i, gmps0->GetEX()[i], gmps0->GetEY()[i]);
    mgps1->SetPoint(gmps1->GetN()-1-i, -gmps1->GetX()[i], gmps1->GetY()[i]);
    mgps1->SetPointError(gmps1->GetN()-1-i, gmps1->GetEX()[i], gmps1->GetEY()[i]);
    mgps2->SetPoint(gmps2->GetN()-1-i, -gmps2->GetX()[i], gmps2->GetY()[i]);
    mgps2->SetPointError(gmps2->GetN()-1-i, gmps2->GetEX()[i], gmps2->GetEY()[i]);
    mgps3->SetPoint(gmps3->GetN()-1-i, -gmps3->GetX()[i], gmps3->GetY()[i]);
    mgps3->SetPointError(gmps3->GetN()-1-i, gmps3->GetEX()[i], gmps3->GetEY()[i]);
  } // for i
  for (int i = 0; i != gd->GetN(); ++i) {
    mgd->SetPoint(gdm->GetN()+i, gd->GetX()[i], gd->GetY()[i]);
    mgd->SetPointError(gdm->GetN()+i, gd->GetEX()[i], gd->GetEY()[i]);
    mgm->SetPoint(gmm->GetN()+i, gm->GetX()[i], gm->GetY()[i]);
    mgm->SetPointError(gmm->GetN()+i, gm->GetEX()[i], gm->GetEY()[i]);
    //
    mgps->SetPoint(gmps->GetN()+i, gpps->GetX()[i], gpps->GetY()[i]);
    mgps->SetPointError(gmps->GetN()+i, gpps->GetEX()[i], gpps->GetEY()[i]);
    mgps0->SetPoint(gmps0->GetN()+i, gpps0->GetX()[i], gpps0->GetY()[i]);
    mgps0->SetPointError(gmps0->GetN()+i, gpps0->GetEX()[i], gpps0->GetEY()[i]);
    mgps1->SetPoint(gmps1->GetN()+i, gpps1->GetX()[i], gpps1->GetY()[i]);
    mgps1->SetPointError(gmps1->GetN()+i, gpps1->GetEX()[i], gpps1->GetEY()[i]);
    mgps2->SetPoint(gmps2->GetN()+i, gpps2->GetX()[i], gpps2->GetY()[i]);
    mgps2->SetPointError(gmps2->GetN()+i, gpps2->GetEX()[i], gpps2->GetEY()[i]);
    mgps3->SetPoint(gmps3->GetN()+i, gpps3->GetX()[i], gpps3->GetY()[i]);
    mgps3->SetPointError(gmps3->GetN()+i, gpps3->GetEX()[i], gpps3->GetEY()[i]);
  }
  //tdrDraw(mgd,"Pz",kFullDiamond,kBlack);
  //tdrDraw(mgm,"Pz",kOpenDiamond,kBlack);

  double uemc = fm->GetParameter(0)*pt/100.;
  double duemc = fm->GetParError(0)*pt/100.;
  //gm->SetPoint(gm->GetN(), uemc/pt, uemc/pt * 100.);
  fm->Draw("SAME");
  //
  double uedt = fd->GetParameter(0)*pt/100.;
  double duedt = fd->GetParError(0)*pt/100.;
  //gd->SetPoint(gd->GetN(), uedt/pt, uedt/pt * 100.);
  fd->Draw("SAME");


  // Combined fit for data
  TF1 *fmgd = new TF1("fmgd","[0]+[1]*x*(x>0)",-amax,+amax);
  fmgd->SetParameters(0,1);
  fmgd->SetLineColor(kBlack);
  fmgd->SetLineStyle(kSolid);
  mgd->Fit(fmgd,"QRN");
  fmgd->Draw("SAME");
  TF1 *fmgdb = (TF1*)fmgd->Clone("fmgdb");
  fmgdb->SetParameter(1,0);
  fmgdb->Draw("SAME");

  // Combined fit for MC
  TF1 *fmgm = new TF1("fmgm","[0]+[1]*x*(x>0)",-amax,+amax);
  fmgm->SetParameters(0,1);
  fmgm->SetLineColor(kBlack);
  fmgm->SetLineStyle(kDashed);
  mgm->Fit(fmgm,"QRN");
  fmgm->Draw("SAME");
  TF1 *fmgmb = (TF1*)fmgm->Clone("fmgmb");
  fmgmb->SetParameter(1,0);
  fmgmb->Draw("SAME");
  //
  TF1 *fmgps = (TF1*)fmgm->Clone("fmgps");
  fmgps->SetLineStyle(kDotted);
  TF1 *fmgps0 = (TF1*)fmgm->Clone("fmgps0");
  //fmgps0->SetLineStyle(kDotted);
  TF1 *fmgps1 = (TF1*)fmgps0->Clone("fmgps1");
  TF1 *fmgps2 = (TF1*)fmgps0->Clone("fmgps2");
  TF1 *fmgps3 = (TF1*)fmgps0->Clone("fmgps3");
  mgps->Fit(fmgps,"QRN");
  mgps0->Fit(fmgps0,"QRN");
  mgps1->Fit(fmgps1,"QRN");
  mgps2->Fit(fmgps2,"QRN");
  mgps3->Fit(fmgps3,"QRN");
  if (drawPS) {
    c0->cd();
    fmgps->Draw("SAME"); fmgps->SetLineColor(kRed);
    fmgps0->Draw("SAME"); fmgps0->SetLineColor(kRed);
    fmgps1->Draw("SAME"); fmgps1->SetLineColor(kRed+1);
    fmgps2->Draw("SAME"); fmgps2->SetLineColor(kRed+2);
    fmgps3->Draw("SAME"); fmgps3->SetLineColor(kRed+3);
    c1->cd(1);
  }
  TF1 *fmgpsb = (TF1*)fmgps->Clone("fmgpsb");
  TF1 *fmgps0b = (TF1*)fmgps0->Clone("fmgps0b");
  TF1 *fmgps1b = (TF1*)fmgps1->Clone("fmgps1b");
  TF1 *fmgps2b = (TF1*)fmgps2->Clone("fmgps2b");
  TF1 *fmgps3b = (TF1*)fmgps3->Clone("fmgps3b");
  fmgpsb->SetParameter(1,0);
  fmgps0b->SetParameter(1,0);
  fmgps1b->SetParameter(1,0);
  fmgps2b->SetParameter(1,0);
  fmgps3b->SetParameter(1,0);
  if (drawPS) {
    c0->cd();
    fmgpsb->Draw("SAME");
    fmgps0b->Draw("SAME");
    fmgps1b->Draw("SAME");
    fmgps2b->Draw("SAME");
    fmgps3b->Draw("SAME");
    c1->cd(1);
  }

  TLatex *tex = new TLatex();
  tex->SetTextSize(0.045);
  tex->SetNDC();

  if (etamax==2.5)
    //tex->DrawLatex(0.20,0.11,"|#eta| < 2.5, R = 0.4 PF+CHS");
    tex->DrawLatex(0.50,0.80,"|#eta| < 2.5, R = 0.4 PF+CHS");
  if (etamax==1.3)
    //tex->DrawLatex(0.20,0.11,"|#eta| < 1.3, R = 0.4 PF+CHS");
    tex->DrawLatex(0.50,0.80,"|#eta| < 1.3, R = 0.4 PF+CHS");
  tex->DrawLatex(0.20,0.05,Form("%1.0f < p_{T,Z} < %1.0f GeV, #LTp_{T,Z}#GT = %1.0f GeV",minpt,maxpt,pt));

  /*
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.50,0.85,Form("#chi^{2} / NDF = %1.1f / %d (%1.1f / %d)",
				fmgd->GetChisquare(), fmgd->GetNDF(),
				fmgm->GetChisquare(), fmgm->GetNDF()));

  tex->SetTextColor(kBlue);//-9);
  tex->DrawLatex(0.50,0.80,Form("#chi^{2} / NDF = %1.1f / %d (%1.1f / %d)",
				fdm->GetChisquare(), fdm->GetNDF(),
				fmm->GetChisquare(), fmm->GetNDF()));

  tex->SetTextColor(kRed);//-9);
  tex->DrawLatex(0.20,0.17,Form("#chi^{2} / NDF = %1.1f / %d (%1.1f / %d)",
				fd->GetChisquare(), fd->GetNDF(),
				fm->GetChisquare(), fm->GetNDF()));
  tex->SetTextSize(0.035);
  tex->DrawLatex(0.20,0.22,Form("k_{FSR} = %1.1f#pm%1.1f%% (%1.1f#pm%1.1f%%)",
				-fd->GetParameter(1),
				fd->GetParError(1),
				-fm->GetParameter(1),
				fm->GetParError(1)));
  tex->DrawLatex(0.20,0.26,Form("#rho_{UE} = %1.2f#pm%1.2f "
				"(%1.2f#pm%1.2f) GeV",
				uedt/(TMath::Pi()*0.4*0.4),
				duedt/(TMath::Pi()*0.4*0.4),
				uemc/(TMath::Pi()*0.4*0.4),
				duemc/(TMath::Pi()*0.4*0.4)));
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.58,0.22,Form("k_{FSR} = %1.1f#pm%1.1f%% (%1.1f#pm%1.1f%%)",
				-fmgd->GetParameter(1),
				fmgd->GetParError(1),
				-fmgm->GetParameter(1),
				fmgm->GetParError(1)));
  */
  gPad->RedrawAxis();

  /*  
  TLegend *leg = tdrLeg(0.65,0.45,0.85,0.70);
  leg->SetHeader("Z(#rightarrowl^{+}l^{-})+jet");
  leg->AddEntry(gdm,"Data MPF","PL");
  leg->AddEntry(gmm,"MC MPF","PL");
  leg->AddEntry(gd,"Data p_{T}-bal.","PL");
  leg->AddEntry(gm,"MC p_{T}-bal.","PL");
  */
  TLegend *leg = tdrLeg(0.45,0.25,0.85,0.60);
  leg->SetHeader("Z(#rightarrowl^{+}l^{-})+jet  (#chi^{2}/NDF)");
  leg->AddEntry(gdm,Form("Data MPF (%1.1f)",fdm->GetChisquare()/fdm->GetNDF()),
		"PL");
  leg->AddEntry(gmm,Form("MC MPF   (%1.1f)",fmm->GetChisquare()/fmm->GetNDF()),
		"PL");
  leg->AddEntry(gd,Form("Data p_{T}-bal. (%1.1f)",
			fd->GetChisquare()/fd->GetNDF()),"PL");
  leg->AddEntry(gm,Form("MC p_{T}-bal.   (%1.1f)",
			fm->GetChisquare()/fm->GetNDF()),"PL");
  //
  leg->AddEntry(fmgd,Form("Data both (%1.1f)",
			  fmgd->GetChisquare()/fmgd->GetNDF()),"L");
  leg->AddEntry(fmgm,Form("MC both   (%1.1f)",
			  fmgm->GetChisquare()/fmgm->GetNDF()),"L");

  c1->cd(2);

  TGraphErrors *grm = (TGraphErrors*)gdm->Clone();
  TGraphErrors *grmps = (TGraphErrors*)gmps->Clone();
  TGraphErrors *grmps0 = (TGraphErrors*)gmps0->Clone();
  TGraphErrors *grmps1 = (TGraphErrors*)gmps1->Clone();
  TGraphErrors *grmps2 = (TGraphErrors*)gmps2->Clone();
  TGraphErrors *grmps3 = (TGraphErrors*)gmps3->Clone();
  for (int i = 0; i != grm->GetN(); ++i) {
    double y = gmm->GetY()[i]/100.+1;
    double ey2 = pow(gmm->GetEY()[i],2);
    double x = gmm->GetX()[i];
    tools::SetPoint(grm, i, x, 100*((gdm->GetY()[i]/100.+1)/y-1),
		    0, sqrt(pow(gdm->GetEY()[i],2)+ey2));
    //
    tools::SetPoint(grmps, i, x, 100*((gmps->GetY()[i]/100.+1)/y-1),
		    0, sqrt(fabs(pow(gmps->GetEY()[i],2)-ey2)));
    tools::SetPoint(grmps0, i, x, 100*((gmps0->GetY()[i]/100.+1)/y-1),
		    0, sqrt(fabs(pow(gmps0->GetEY()[i],2)-ey2)));
    tools::SetPoint(grmps1, i, x, 100*((gmps1->GetY()[i]/100.+1)/y-1),
		    0, sqrt(fabs(pow(gmps1->GetEY()[i],2)-ey2)));
    tools::SetPoint(grmps2, i, x, 100*((gmps2->GetY()[i]/100.+1)/y-1),
		    0, sqrt(fabs(pow(gmps2->GetEY()[i],2)-ey2)));
    tools::SetPoint(grmps3, i, x, 100*((gmps3->GetY()[i]/100.+1)/y-1),
		    0, sqrt(fabs(pow(gmps3->GetEY()[i],2)-ey2)));
  }

  TGraphErrors *gr = (TGraphErrors*)gd->Clone();
  TGraphErrors *grpps = (TGraphErrors*)gpps->Clone();
  TGraphErrors *grpps0 = (TGraphErrors*)gpps0->Clone();
  TGraphErrors *grpps1 = (TGraphErrors*)gpps1->Clone();
  TGraphErrors *grpps2 = (TGraphErrors*)gpps2->Clone();
  TGraphErrors *grpps3 = (TGraphErrors*)gpps3->Clone();
  for (int i = 0; i != gr->GetN(); ++i) {
    double y = gm->GetY()[i]/100.+1;
    double ey2 = pow(gm->GetEY()[i],2);
    double x = gm->GetX()[i];
    tools::SetPoint(gr, i, x, 100*((gd->GetY()[i]/100.+1)/y-1),
		    0, sqrt(pow(gd->GetEY()[i],2)+ey2));
    //
    tools::SetPoint(grpps, i, x, 100*((gpps->GetY()[i]/100.+1)/y-1),
		    0, sqrt(fabs(pow(gpps->GetEY()[i],2)-ey2)));
    tools::SetPoint(grpps0, i, x, 100*((gpps0->GetY()[i]/100.+1)/y-1),
		    0, sqrt(fabs(pow(gpps0->GetEY()[i],2)-ey2)));
    tools::SetPoint(grpps1, i, x, 100*((gpps1->GetY()[i]/100.+1)/y-1),
		    0, sqrt(fabs(pow(gpps1->GetEY()[i],2)-ey2)));
    tools::SetPoint(grpps2, i, x, 100*((gpps2->GetY()[i]/100.+1)/y-1),
		    0, sqrt(fabs(pow(gpps2->GetEY()[i],2)-ey2)));
    tools::SetPoint(grpps3, i, x, 100*((gpps3->GetY()[i]/100.+1)/y-1),
		    0, sqrt(fabs(pow(gpps3->GetEY()[i],2)-ey2)));
  }

  TGraphErrors *mgr = new TGraphErrors(gr->GetN()+grm->GetN());
  for (int i = grm->GetN()-1; i != -1; --i) {
    mgr->SetPoint(grm->GetN()-1-i, -grm->GetX()[i], grm->GetY()[i]);
    mgr->SetPointError(grm->GetN()-1-i, grm->GetEX()[i], grm->GetEY()[i]);
  } // for i
  for (int i = 0; i != gr->GetN(); ++i) {
    mgr->SetPoint(grm->GetN()+i, gr->GetX()[i], gr->GetY()[i]);
    mgr->SetPointError(grm->GetN()+i, gr->GetEX()[i], gr->GetEY()[i]);
  }

  tdrDraw(grm,"Pz",kFullSquare,kBlue); grm->SetMarkerSize(0.5);
  tdrDraw(gr,"Pz",kFullCircle,kRed);   gr->SetMarkerSize(0.5);
  //tdrDraw(mgr,"Pz",kFullDiamond,kBlack);
  //
  // Parton shower variations
  if (drawPS) {
    if (drawData) {
      tdrDraw(grmps0,"Pz",kOpenDiamond,kBlue);
      tdrDraw(grmps1,"Pz",kOpenDiamond,kBlue+1);
      tdrDraw(grmps2,"Pz",kOpenDiamond,kBlue+2);
      tdrDraw(grmps3,"Pz",kOpenDiamond,kBlue+3);
      tdrDraw(grpps0,"Pz",kOpenDiamond,kRed);
      tdrDraw(grpps1,"Pz",kOpenDiamond,kRed+1);
      tdrDraw(grpps2,"Pz",kOpenDiamond,kRed+2);
      tdrDraw(grpps3,"Pz",kOpenDiamond,kRed+3);
    }
  }

  TF1 *frm = new TF1("frm","((([0]+[1]*x)/100.+1)"
  		     " / (([2]+[3]*x)/100.+1) - 1)*100",amin,amax);
  frm->SetParameters(fdm->GetParameter(0), fdm->GetParameter(1),
  		     fmm->GetParameter(0), fmm->GetParameter(1));
  //TF1 *frm = new TF1("frm","((([0]+[1]*(x-[0]/[4]))/100.+1)"
  //		   " / (([2]+[3]*(x-[0]/4))/100.+1) - 1)*100",0,0.35);
  //frm->SetParameters(fdm->GetParameter(0), fdm->GetParameter(1),
  //		     fmm->GetParameter(0), fmm->GetParameter(1),
  //		     pt);
  frm->SetLineColor(kBlue);//-9);
  frm->SetLineStyle(kDashed);
  frm->Draw("SAME");

  TF1 *frmn = new TF1("frmn","[0]+[1]*x",amin,amax);
  frmn->FixParameter(1,0); // no slope for MPF
  grm->Fit(frmn,"QRN");
  frmn->SetLineColor(kBlue);
  frmn->SetLineStyle(kSolid);
  frmn->Draw("SAME");

  TF1 *fr = new TF1("fr","((([0]+[1]*x)/100.+1)"
  		    " / (([2]+[3]*x)/100.+1) - 1)*100",0,0.35);
  fr->SetParameters(fd->GetParameter(0), fd->GetParameter(1),
  		    fm->GetParameter(0), fm->GetParameter(1));
  //TF1 *fr = new TF1("fr","((([0]+[1]*(x-[0]/[4]))/100.+1)"
  //		    " / (([2]+[3]*(x-[0]/[4]))/100.+1) - 1)*100",amin,amax);
  //fr->SetParameters(fd->GetParameter(0), fd->GetParameter(1),
  //		    fm->GetParameter(0), fm->GetParameter(1),
  //		    pt);
  fr->SetLineColor(kRed);//-9);
  fr->SetLineStyle(kDashed);
  fr->Draw("SAME");

  TF1 *frn = new TF1("frn","[0]+[1]*x",amin,amax);
  gr->Fit(frn,"QRN");
  frn->SetLineColor(kRed);
  frn->SetLineStyle(kSolid);
  frn->Draw("SAME");

  // Combined fit
  TF1 *frr = new TF1("frr","[0]+[1]*x*(x>0)",-amax,+amax);
  frr->SetParameters(0,1);
  frr->SetLineColor(kBlack);
  frr->SetLineStyle(kDashed);
  mgr->Fit(frr,"QRN");
  frr->Draw("SAME");
  TF1 *frrb = (TF1*)frr->Clone("frrb");
  frrb->SetParameter(1,0);
  frrb->Draw("SAME");

  gPad->RedrawAxis();

  c1->Update();
  //c1->SaveAs("pdf/drawZalpha_2018.pdf");
  c1->SaveAs(Form("pdf/drawZalpha_Run2_%1.0f-%1.0f.pdf",minpt,maxpt));
  

  c0->cd(1);

  TF1 *fpps = new TF1("fpps","[0]+[1]*x",amin,amax);
  // ROOT 6.18/04 does not seem to support simultaneous range and option "S"
  // nor does it return a fit result or gMinuit without option "S" :(
  // so, need to clip graph first to intended range
  // => ah, was rather fitting empty graphs
  TGraphErrors *gpps_clip = clipGraph(gpps,amin,amax);
  if (gpps_clip->GetN()>0) {

    TFitResultPtr frpps = gpps_clip->Fit(fpps,"SRQN");//"QRN");
    fpps->SetLineColor(kRed);
    fpps->SetLineStyle(kSolid);

    //TMatrixD ematpps(2,2); assert(gMinuit);
    //gMinuit->mnemat(ematpps.GetMatrixArray(),2);
    TMatrixDSym ematpps = frpps->GetCovarianceMatrix();
    TGraphErrors *gppse = fitErrGraph(fpps,ematpps,gpps);
    tdrDraw(gppse,"E3",kNone,kRed-9,kSolid,-1,1001,kRed-9);
    fpps->Draw("SAME");
  }

  TF1 *fmps = new TF1("fmps","[0]+[1]*x",amin,amax);
  fmps->FixParameter(1,0);
  gmps->Fit(fmps,"QRN");
  fmps->SetLineColor(kBlue);
  fmps->SetLineStyle(kSolid);
  fmps->Draw("SAME");

  // Combined fit
  TF1 *fps = new TF1("fps","[0]+[1]*x*(x>0)",-amax,+amax);
  fps->SetParameters(0,1);
  fps->SetLineColor(kBlack);
  fps->SetLineStyle(kDashed);
  mgps->Fit(fps,"QRN");
  fps->Draw("SAME");
  TF1 *fpsb = (TF1*)fps->Clone("fpsb");
  fpsb->SetParameter(1,0);
  fpsb->Draw("SAME");


  tdrDraw(gpps,"Pz",kFullSquare,kRed);
  tdrDraw(gpps0,"Pz",kOpenTriangleDown,kRed);
  tdrDraw(gpps1,"Pz",kFullTriangleDown,kRed+1);
  tdrDraw(gpps2,"Pz",kOpenTriangleUp,kRed+2);
  tdrDraw(gpps3,"Pz",kFullTriangleUp,kRed+3);
  tdrDraw(gmps,"Pz",kFullSquare,kBlue);
  tdrDraw(gmps0,"Pz",kOpenTriangleDown,kBlue); 
  tdrDraw(gmps1,"Pz",kFullTriangleDown,kBlue+1);
  tdrDraw(gmps2,"Pz",kOpenTriangleUp,kBlue+2);
  tdrDraw(gmps3,"Pz",kFullTriangleUp,kBlue+3);

  gpps->SetMarkerSize(0.5);
  gpps0->SetMarkerSize(0.5);
  gpps1->SetMarkerSize(0.5);
  gpps2->SetMarkerSize(0.5);
  gpps3->SetMarkerSize(0.5);
  gmps->SetMarkerSize(0.4);
  gmps0->SetMarkerSize(0.4);
  gmps1->SetMarkerSize(0.4);
  gmps2->SetMarkerSize(0.4);
  gmps3->SetMarkerSize(0.4);

  c0->SaveAs(Form("pdf/drawZalpha_Run2_PSWeight_%1.0f-%1.0f.pdf",minpt,maxpt));

  //return make_pair<double,double>(fmgm->GetParameter(1),fmgm->GetParError(1));
  //return make_pair<double,double>(fmgd->GetParameter(1),fmgd->GetParError(1));
  //return make_pair<double,double>(fmgd->GetParameter(1),fmgm->GetParameter(1));

  FSR fsr;
  const double alpharef = 0.30;
  fsr.data = fmgd->GetParameter(1) * alpharef;
  fsr.data_err = fmgd->GetParError(1) * alpharef;
  fsr.data_chi2 = fmgd->GetChisquare() / fmgd->GetNDF();
  fsr.data_mpf = fmgd->GetParameter(0);
  fsr.data_mpferr = fmgd->GetParError(0);
  fsr.mc = fmgm->GetParameter(1) * alpharef;
  fsr.mc_err = fmgm->GetParError(1) * alpharef;
  fsr.mc_chi2 = fmgm->GetChisquare() / fmgm->GetNDF();
  fsr.mc_mpf = fmgm->GetParameter(0);
  fsr.mc_mpferr = fmgm->GetParError(0);
  //
  fsr.ps0 = fmgps0->GetParameter(1) * alpharef;
  fsr.ps1 = fmgps1->GetParameter(1) * alpharef;
  fsr.ps2 = fmgps2->GetParameter(1) * alpharef;
  fsr.ps3 = fmgps3->GetParameter(1) * alpharef;
  fsr.ps0_chi2 = fmgps0->GetChisquare() / fmgps0->GetNDF();
  fsr.ps1_chi2 = fmgps1->GetChisquare() / fmgps1->GetNDF();
  fsr.ps2_chi2 = fmgps2->GetChisquare() / fmgps2->GetNDF();
  fsr.ps3_chi2 = fmgps3->GetChisquare() / fmgps3->GetNDF();

  return fsr;
} // drawZalphas
