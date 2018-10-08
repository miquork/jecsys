#include "TFile.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMultiGraph.h"

#include "../tools.C"

#include <string>
using namespace std;

void compareIOV() {

  //const int niov = 6;
  //string iovs[niov] = {"BCDEF","B","C","D","E","F"};
  //int colors[niov] = {kGray+2,kBlack, kBlue, kGreen+2, kOrange+2, kRed};

  const int niov = 4;
  string iovs[niov] = {"BCDEF","D","E","F"};
  int colors[niov] = {kGray+2, kBlue, kGreen+2, kRed};

  TH1D *h = new TH1D("h",";p_{T} (GeV);R_{data}",100,30,250);//3500);
  h->SetMinimum(0.85);//0.65);//0.7);
  h->SetMaximum(1.15);//1.00);//1.45);
  h->Draw();

  TGraphErrors *g10(0), *g20(0);
  for (int i = 0; i != niov; ++i) {

    //TFile *f = new TFile(Form("rootfiles/zjet_combination_Fall17_JECV5_Zmm_%s_2018-02-24.root",iovs[iov]),"READ");
    TFile *f = new TFile(Form("../rootfiles/jecdata%s.root",iovs[i].c_str()),"READ");
    assert(f && !f->IsZombie());
    
    //TGraphErrors *g = (TGraphErrors*)f->Get("data/eta29-30/ptchs_zmmjet_a30");
    //TGraphErrors *g1 = (TGraphErrors*)f->Get("data/eta26-29/ptchs_zmmjet_a30");
    //TGraphErrors *g1 = (TGraphErrors*)f->Get("data/eta26-29/ptchs_zeejet_a30");
    TGraphErrors *g1 = (TGraphErrors*)f->Get("data/eta26-29/ptchs_gamjet_a30");
    assert(g1);
    if (g10==0) g10 = (TGraphErrors*)g1->Clone("g10");
    g1 = tools::ratioGraphs(g1,g10);

    g1->SetLineColor(colors[i]);
    g1->SetMarkerColor(colors[i]);
    if (g1->GetN()>0) g1->Draw("SAMEPz");

    //TGraphErrors *g2 = (TGraphErrors*)f->Get("data/eta26-29/mpfchs1_zmmjet_a30");
    //TGraphErrors *g2 = (TGraphErrors*)f->Get("data/eta26-29/mpfchs1_zeejet_a30");
    TGraphErrors *g2 = (TGraphErrors*)f->Get("data/eta26-29/mpfchs1_gamjet_a30");
    assert(g2);
    if (g20==0) g20 = (TGraphErrors*)g2->Clone("g10");
    g2 = tools::ratioGraphs(g2,g20);

    g2->SetLineColor(colors[i]);
    g2->SetMarkerColor(colors[i]);
    if (g2->GetN()>0) g2->Draw("SAMEPz");

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(g1);
    mg->Add(g2);

    TF1 *f1 = new TF1(Form("f1_%d",i),"[0]+[1]*log(x)",30,250);
    f1->SetParameters(1,0.);
    mg->Fit(f1,"QRN");
    f1->SetLineColor(colors[i]);
    f1->Draw("SAME");

    // ECAL 0.17 -> 0.10
    // jet energy ~50% of ECAL
    // (0.17-0.10)/0.10*0.5=35%
    // s5: data 0.89 vs MC 0.85, so (0.15-0.11)/0.11=36%
    // 35%*36%=13%
    // low pT: photons in ECAL (25%/75%), hadrons in HCAL (37.5%/75%)
    // high pT: photons 25% and hadrons 37.5% in ECAL, hadrons 37.5% in HCAL
  }  

}
