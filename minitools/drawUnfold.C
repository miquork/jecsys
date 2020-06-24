// Purpose: Draw comparison of (un)folding for inclusive jets
//          with JES and JER per IOV
#include "TFile.h"

#include "../tdrstyle_mod15.C"

#include <map>
#include <string>
using namespace std;

bool _useZoom = true;
bool _useIncJet = true;

void drawUnfold() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  //TFile *f = new TFile("rootfiles/unfold.root","READ");
  //TFile *f = new TFile("rootfiles/unfold_jes_mc.root","READ");
  //TFile *f = new TFile("rootfiles/unfold_jes_v3.root","READ");
  TFile *f = new TFile("rootfiles/commonUL2017_V4_V2M5res_hotzone_jes-2.root","READ");
  assert(f && !f->IsZombie());

  //TFile *f2 = new TFile("rootfiles/unfold.root","READ");
  //TFile *f2 = new TFile("rootfiles/unfold_nojes_mc.root","READ");
  //TFile *f2 = new TFile("rootfiles/unfold_nojes_v3.root","READ");
  TFile *f2 = new TFile("rootfiles/commonUL2017_V4_V2M5res_hotzone-2.root","READ");
  assert(f2 && !f2->IsZombie());

  //string aiov[] = {"BCDEF","B","C","D","E","F"};
  //int aidxp[] = {0,1,2,3,4,5};
  //string aiov[] = {"BCDEF","F","E","D","C","B"};
  //string aiov[] = {"UL17","F","E","D","C","B"};
  string aiov[] = {"all","F","E","D","C","B"};
  int aidx[] = {0,5,4,3,2,1};
  int niov = sizeof(aiov)/sizeof(aiov[0]);

  map<string,int> marker;
  //marker["BCDEF"] = kFullSquare;
  //marker["UL17"] = kFullSquare;
  marker["all"] = kFullSquare;
  marker["B"] = kFullCircle;
  marker["C"] = kFullDiamond;
  marker["D"] = kFullSquare;
  marker["E"] = kFullCross;
  marker["F"] = kFullCrossX;

  map<string,int> color;
  //color["BCDEF"] = kBlack;
  //color["UL17"] = kBlack;
  color["all"] = kBlack;
  color["B"] = kBlue;
  color["C"] = kRed;
  color["D"] = kGreen+2;
  color["E"] = kMagenta+1;
  color["F"] = kOrange+2;

  double ptmin(15),ptmax(4500);
  TH1D *h = tdrHist("h","JER smearing",
		    1.03,1.35,"p_{T} (GeV)",ptmin,ptmax);
  if (_useZoom) {
    h->GetXaxis()->SetRangeUser(15,100);
    h->GetYaxis()->SetRangeUser(1.1,1.35);
  }
  lumi_13TeV = "UL2017";
  TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);
  gPad->SetLogx();

  double dx = 0.10;
  TLegend *leg = tdrLeg(0.60+dx,0.57,0.90+dx,0.87);
  //leg->SetHeader("+JES");
  leg->SetHeader(" #Delta#mu");
  TLegend *leg2 = tdrLeg(0.52+dx,0.57,0.82+dx,0.87);
  //leg2->SetHeader("JER");
  leg2->SetHeader("#Delta#sigma +");
  TLine *l = new TLine(); l->SetLineStyle(kDashed);
  l->DrawLine(ptmin,0,ptmax,0);

  for (int iov = 0; iov != niov; ++iov) {
    
    string s = aiov[iov];
    const char *c = s.c_str();
    int j = aidx[iov];
    TH1D *hr(0), *hr2(0);
    if (_useIncJet) {
      TH1D *hu1 = (TH1D*)f->Get(Form("ak4/Eta_0.0-1.3/"
				     "hpt_data_2017_%s_ptcl_fwd",c));
      TH1D *hu2 = (TH1D*)f2->Get(Form("ak4/Eta_0.0-1.3/"
				      "hpt_data_2017_%s_ptcl_fwd",c));
      TH1D *h1 = (TH1D*)f->Get(Form("ak4/Eta_0.0-1.3/hpt_data_2017_%s_det",c));
      TH1D *h2 = (TH1D*)f2->Get(Form("ak4/Eta_0.0-1.3/hpt_data_2017_%s_det",c));
      hr = (TH1D*)h1->Clone(Form("hr_%d",j));
      hr2 = (TH1D*)h2->Clone(Form("hr_%d",j));
      hr->Divide(hu1);
      hr2->Divide(hu2);
    }
    else {
      hr = (TH1D*)f->Get(Form("hr_0.0-1.3_%d",2021+j));
      hr2 = (TH1D*)f2->Get(Form("hr_0.0-1.3_%d",2021+j));
    }
    assert(hr);
    assert(hr2);
    tdrDraw(hr,"HIST",marker[s],color[s],kSolid,-1,kNone);
    tdrDraw(hr,"HISTP",marker[s],color[s],kSolid,-1,kNone);
    hr2->SetMarkerSize(0.8);
    tdrDraw(hr2,"HIST",marker[s],color[s],kDashed,-1,kNone);
    tdrDraw(hr2,"HISTP",marker[s],color[s],kDashed,-1,kNone);
    leg->AddEntry(hr,s.c_str(),"PL");
    leg2->AddEntry(hr2," ","PL");
  } // for iov

  c1->SaveAs("pdf/drawUnfold.pdf");
} // drawUnfold
