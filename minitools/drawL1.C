// Created by Mikko Voutilainen, on Oct 11th, 2018
// Purpose is to raw:
// 1) rho vs mu for data and MC (+ratio), determine k0 = (<rho>_data-rho0)/(<rho>_MC-rho0)
// 2) off(<rho>_MC)/(<rho>_MC-rho0) vs rc(<rho>_MC)/(<rho>_MC-rho0) (PF and PFchs)
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TF1.h"
#include "TMultiGraph.h"
#include "TLatex.h"
#include "TMinuit.h"
#include "TMatrixD.h"
#include "TMath.h"
#include "TLine.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "tools.h"
#include "tdrstyle_mod14.C"

#include <string>

using namespace std;

FactorizedJetCorrector *getJEC(string version) {

  string s;
  const char *cd = "CondFormats/JetMETObjects/data";
  const char *gt = "Summer16_07Aug2017";

  FactorizedJetCorrector *jecl1;
  s = Form("%s/%s%s.txt",cd,gt,version.c_str());
  cout << s << endl << flush;
  JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
  vector<JetCorrectorParameters> v;
  v.push_back(*l1);
  jecl1 = new FactorizedJetCorrector(v);

  return jecl1;
}

void drawL1(string run = "GH", string mode = "DataMC") {

  assert(mode=="DataMC" || mode=="CHSPF" || mode=="RCL1");
  setTDRStyle();

  double gRho = 20.;
  double gRho0 = 1.519;
  double gPt = 30.;
  double gJetArea = TMath::Pi()*0.4*0.4;

  /*
  string s;
  const char *cd = "CondFormats/JetMETObjects/data";
    
  FactorizedJetCorrector *jecl1mcpf;
  {
    s = Form("%s/Summer16_07Aug2017_V15_MC_L1FastJet_AK4PF.txt",cd);
    cout << s << endl << flush;
    JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l1);
    jecl1mcpf = new FactorizedJetCorrector(v);
  }

  FactorizedJetCorrector *jecl1mcchs;
  {
    s = Form("%s/Summer16_07Aug2017_V15_MC_L1FastJet_AK4PFchs.txt",cd);
    cout << s << endl << flush;
    JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
    vector<JetCorrectorParameters> v;
    v.push_back(*l1);
    jecl1mcchs = new FactorizedJetCorrector(v);
  }
  */
  FactorizedJetCorrector *jecrcmcpf = getJEC("_V15_MC_L1RC_AK4PF");
  FactorizedJetCorrector *jecrcmcchs = getJEC("_V15_MC_L1RC_AK4PFchs");
  FactorizedJetCorrector *jecl1mcpf = getJEC("_V15_MC_L1FastJet_AK4PF");
  FactorizedJetCorrector *jecl1mcchs = getJEC("_V15_MC_L1FastJet_AK4PFchs");
  FactorizedJetCorrector *jecrcdtpf = getJEC("GH_V18_DATA_L1RC_AK4PF");
  FactorizedJetCorrector *jecrcdtchs = getJEC("GH_V18_DATA_L1RC_AK4PFchs");
  FactorizedJetCorrector *jecl1dtpf = getJEC("GH_V18_DATA_L1FastJet_AK4PF");
  FactorizedJetCorrector *jecl1dtchs = getJEC("GH_V18_DATA_L1FastJet_AK4PFchs");

  TH1D *h = new TH1D("h",";#eta;off(#LT#rho#GT) / A#times(#LT#rho#GT-#rho_{0})",
		     94,-4.7,4.7);//
  //104,-5.2,5.2);
  TH1D *hrcmcpf = (TH1D*)h->Clone("hrcmcpf");
  TH1D *hrcmcchs = (TH1D*)h->Clone("hrcmcchs");
  TH1D *hl1mcpf = (TH1D*)h->Clone("hl1mcpf");
  TH1D *hl1mcchs = (TH1D*)h->Clone("hl1mcchs");
  TH1D *hrcdtpf = (TH1D*)h->Clone("hrcdtpf");
  TH1D *hrcdtchs = (TH1D*)h->Clone("hcdtchs");
  TH1D *hl1dtpf = (TH1D*)h->Clone("hl1dtpf");
  TH1D *hl1dtchs = (TH1D*)h->Clone("hl1dtchs");

  //TH1D *hl1mcpf = new TH1D("hl1mcpf",";#eta;off(#LT#rho#GT) / A#times(#LT#rho#GT-#rho_{0})",
  //		   104,-5.2,5.2);
  //TH1D *hl1mcch = new TH1D("hl1mch",";#eta;off(#LT#rho#GT) / A#times(#LT#rho#GT-#rho_{0})",
  //		   104,-5.2,5.2);

  const int njec = 8;
  FactorizedJetCorrector *vjec[njec] = {jecrcmcpf, jecrcmcchs, jecl1mcpf, jecl1mcchs,
					jecrcdtpf, jecrcdtchs, jecl1dtpf, jecl1dtchs};
  TH1D *vh[njec] = {hrcmcpf, hrcmcchs, hl1mcpf, hl1mcchs,
		    hrcdtpf, hrcdtchs, hl1dtpf, hl1dtchs};
  int colors[njec] = {kCyan+2, kCyan+1, kBlue, kRed, 
		      kMagenta+2, kMagenta+1, kGreen+2, kBlack};
  int style[njec] = {kSolid, kSolid, kSolid, kSolid, kSolid, kSolid, kSolid, kSolid};
  //double drho[njec] = {gRho0, gRho0, 0, 0, gRho0, gRho0, 0, 0};

  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    
    double eta = h->GetBinCenter(i);
    for (int j = 0; j != njec; ++j) {

      FactorizedJetCorrector *jec = vjec[j];
      jec->setJetEta(eta);
      jec->setJetPt(gPt);
      jec->setRho(gRho);//-drho[j]);
      jec->setJetA(gJetArea);
      double c = jec->getCorrection();
      // c = 1-off/pt => off = (1-c)*pt
      double off = (1-c)*gPt;
    
      TH1D *h = vh[j];
      h->SetBinContent(i, off / (gJetArea*(gRho-gRho0)));
    }
  }


  double ptmax = 1200.*2.;
  TH1D *hup = new TH1D("hup",";#eta;offset(#LT#rho#GT) / A#times(#LT#rho#GT-#rho_{0})",
		       94,-4.7,4.7);//104,-5.2,5.2);
  hup->SetMinimum(0.1+1e-4);
  hup->SetMaximum(1.7-1e-4);

  TH1D *hdw = new TH1D("hdw",";#eta;Data / MC",94,-4.7,4.7);//104,-5.2,5.2);
  hdw->SetMinimum(0.9+1e-4);
  hdw->SetMaximum(1.2-1e-4);

  map<string, const char*> lumimap;
  lumimap["BCD"] = "Run2016BCD Legacy, 12.9 fb^{-1}";
  lumimap["EF"] = "Run2016EF Legacy, 6.8 fb^{-1}";
  lumimap["FG"] = "Run2016fG Legacy, 8.0 fb^{-1}";
  lumimap["H"] = "Run2016H Legacy, 8.8 fb^{-1}";
  lumimap["GH"] = "Run2016fGH Legacy V18, 16.8 fb^{-1}";
  lumimap["BCDEFGH"] = "Run2016BCDEFGH Legacy, 36.5 fb^{-1}";
  lumi_13TeV = lumimap[run];

  TCanvas *c1 = tdrDiCanvas("c1",hup,hdw,4,11);

  c1->cd(1);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kBlack);//kBlue);
  //l->DrawLine(-5.2,1,5.2,1);
  //l->DrawLine(-4.7,1,4.7,1);

  double ymc = hl1mcpf->Integral()/hl1mcpf->GetNbinsX();
  double yrc = hrcmcpf->Integral()/hrcmcpf->GetNbinsX();

  /*
  hrcmcpf->Scale(ymc/yrc);
  hrcmcchs->Scale(ymc/yrc);
  hrcdtpf->Scale(ymc/yrc);
  hrcdtchs->Scale(ymc/yrc);
  */

  l->SetLineStyle(kDotted);
  //l->DrawLine(-5.2,ymc,5.2,ymc);
  //l->DrawLine(-4.7,ymc,4.7,ymc);

  /*
  tdrDraw(hrcmcpf,"HIST",0,kCyan+2,kSolid,-1,0,0);
  tdrDraw(hrcmcchs,"HIST",0,kCyan+1,kSolid,-1,0,0);
  tdrDraw(hl1mcpf,"HIST",0,kBlue,kSolid,-1,0,0);
  tdrDraw(hl1mcchs,"HIST",0,kRed,kSolid,-1,0,0);

  tdrDraw(hrcdtpf,"HIST",0,kMagenta+2,kSolid,-1,0,0);
  tdrDraw(hrcdtchs,"HIST",0,kMagenta+1,kSolid,-1,0,0);
  tdrDraw(hl1dtpf,"HIST",0,kGreen+2,kSolid,-1,0,0);
  tdrDraw(hl1dtchs,"HIST",0,kBlack,kSolid,-1,0,0);
  */

  // MC red, data black, RC blue
  tdrDraw(hrcmcchs,"HIST",0,kRed,kDashed,-1,0,0);
  tdrDraw(hl1mcchs,"HIST",0,kOrange+2,kDashed,-1,0,0);

  tdrDraw(hrcdtchs,"HIST",0,kRed,kSolid,-1,0,0);
  tdrDraw(hl1dtchs,"HIST",0,kOrange+2,kSolid,-1,0,0);

  tdrDraw(hl1mcpf,"HIST",0,kCyan+2,kDashed,-1,0,0);
  tdrDraw(hl1dtpf,"HIST",0,kCyan+2,kSolid,-1,0,0);

  tdrDraw(hrcmcpf,"HIST",0,kBlue,kDashed,-1,0,0);
  tdrDraw(hrcdtpf,"HIST",0,kBlue,kSolid,-1,0,0);

  hrcdtchs->SetLineWidth(3);
  hl1dtchs->SetLineWidth(3);

  TLegend *leg1up = tdrLeg(0.43,0.66,0.63,0.76);
  leg1up->AddEntry(hrcdtpf,"PF RC DT","L");
  leg1up->AddEntry(hrcmcpf,"PF RC MC","L");

  TLegend *leg1m1 = tdrLeg(0.45,0.43,0.65,0.53);
  leg1m1->AddEntry(hl1dtpf,"PF L1 DT","L");
  leg1m1->AddEntry(hl1mcpf,"PF L1 MC","L");

  TLegend *leg1m2 = tdrLeg(0.43,0.30,0.63,0.40);
  leg1m2->AddEntry(hrcdtchs,"CHS RC DT","L");
  leg1m2->AddEntry(hrcmcchs,"CHS RC MC","L");

  TLegend *leg1dw = tdrLeg(0.43,0.10,0.63,0.20);
  leg1dw->AddEntry(hl1dtchs,"CHS L1 DT","L");
  leg1dw->AddEntry(hl1mcchs,"CHS L1 MC","L");

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);
  tex->SetTextColor(kGray+2);
  tex->DrawLatex(0.40,0.85,Form("p_{T}=%1.0f GeV, #LT#rho#GT = %1.1f GeV",gPt,gRho));
  tex->DrawLatex(0.40,0.80,Form("Anti-k_{T} R=0.4, #rho_{0}= %1.3f GeV",gRho0));

  gPad->RedrawAxis();


  c1->cd(2);

  l->SetLineStyle(kDashed);
  l->SetLineColor(kBlack);
  l->DrawLine(-4.7,1,4.7,1);

  TH1D *hrcpf = (TH1D*)hrcdtpf->Clone("hrcpf");
  TH1D *hrcchs = (TH1D*)hrcdtchs->Clone("hrcchs");
  TH1D *hl1pf = (TH1D*)hl1dtpf->Clone("hl1pf");
  TH1D *hl1chs = (TH1D*)hl1dtchs->Clone("hl1chs");

  hrcpf->Divide(hrcmcpf);
  hrcchs->Divide(hrcmcchs);
  hl1pf->Divide(hl1mcpf);
  hl1chs->Divide(hl1mcchs);

  TH1D *hrcmc = (TH1D*)hrcmcchs->Clone("hrcmc");
  TH1D *hrcdt = (TH1D*)hrcdtchs->Clone("hrcdt");
  TH1D *hl1mc = (TH1D*)hl1mcchs->Clone("hl1mc");
  TH1D *hl1dt = (TH1D*)hl1dtchs->Clone("hl1cdt");

  hrcmc->Divide(hrcmcpf);
  hrcdt->Divide(hrcdtpf);
  hl1mc->Divide(hl1mcpf);
  hl1dt->Divide(hl1dtpf);

  TH1D *hmcpf = (TH1D*)hrcmcpf->Clone("hmcpf");
  TH1D *hmcchs = (TH1D*)hrcmcchs->Clone("hmcchs");
  TH1D *hdtpf = (TH1D*)hrcdtpf->Clone("hdtpf");
  TH1D *hdtchs = (TH1D*)hrcdtchs->Clone("hdtchs");

  hmcpf->Divide(hl1mcpf);
  hmcchs->Divide(hl1mcchs);
  hdtpf->Divide(hl1dtpf);
  hdtchs->Divide(hl1dtchs);


  /*
  tdrDraw(hrcpf,"HIST",0,kBlue,kSolid,-1,0,0);
  tdrDraw(hrcchs,"HIST",0,kRed,kSolid,-1,0,0);
  tdrDraw(hl1pf,"HIST",0,kCyan+2,kSolid,-1,0,0);
  tdrDraw(hl1chs,"HIST",0,kMagenta+2,kSolid,-1,0,0);
  */
  
  if (mode=="DataMC") {

    hdw->SetYTitle("Data / MC");

    // L1CHS ratio red, L1PF ratio blue, RCCHS magenta, RCPF blue
    tdrDraw(hrcpf,"HIST",0,kBlue,kSolid,-1,0,0);
    tdrDraw(hrcchs,"HIST",0,kRed,kSolid,-1,0,0);
    tdrDraw(hl1pf,"HIST",0,kCyan+2,kSolid,-1,0,0);
    tdrDraw(hl1chs,"HIST",0,kOrange+2,kSolid,-1,0,0);
    
    hrcchs->SetLineWidth(3);
    hl1chs->SetLineWidth(3);
    
    TLegend *leg2 = tdrLeg(0.17,0.50,0.32,0.90);
    leg2->SetTextSize(0.045*2);
    leg2->AddEntry(hrcchs,"CHS RC","L");
    leg2->AddEntry(hl1chs,"CHS L1","L");
    leg2->AddEntry(hrcpf,"PF RC","L");
    leg2->AddEntry(hl1pf,"PF L1","L");
    
    c1->SaveAs(Form("pdf/drawL1_DataMC.pdf"));
  }
  if (mode=="CHSPF") {

    hdw->SetYTitle("CHS / PF");
    hdw->GetYaxis()->SetRangeUser(0.4+1e-4,0.55-1e-4);

    // L1CHS ratio red, L1PF ratio blue, RCCHS magenta, RCPF blue
    tdrDraw(hrcdt,"HIST",0,kRed,kSolid,-1,0,0); 
    tdrDraw(hrcmc,"HIST",0,kRed,kDotted,-1,0,0);
    tdrDraw(hl1dt,"HIST",0,kOrange+2,kSolid,-1,0,0);
    tdrDraw(hl1mc,"HIST",0,kOrange+2,kDotted,-1,0,0);
    
    hrcdt->SetLineWidth(3);
    hl1dt->SetLineWidth(3);
    hrcmc->SetLineWidth(2);
    hl1mc->SetLineWidth(2);
    
    TLegend *leg2 = tdrLeg(0.17,0.50,0.32,0.90);
    leg2->SetTextSize(0.045*2);
    leg2->AddEntry(hrcdt,"RC DT","L");
    leg2->AddEntry(hrcmc,"RC MC","L");
    leg2->AddEntry(hl1dt,"L1 DT","L");
    leg2->AddEntry(hl1mc,"L1 MC","L");
    
    c1->SaveAs(Form("pdf/drawL1_CHSPF.pdf"));
  }
  if (mode=="RCL1") {

    hdw->SetYTitle("RC / L1");
    hdw->GetYaxis()->SetRangeUser(0.9+1e-4,2.5);//0.55-1e-4);

    // L1CHS ratio red, L1PF ratio blue, RCCHS magenta, RCPF blue
    tdrDraw(hdtpf,"HIST",0,kBlue,kSolid,-1,0,0); 
    tdrDraw(hdtchs,"HIST",0,kRed,kSolid,-1,0,0);
    tdrDraw(hmcpf,"HIST",0,kBlue,kDotted,-1,0,0);
    tdrDraw(hmcchs,"HIST",0,kRed,kDotted,-1,0,0);
    
    hdtpf->SetLineWidth(3);
    hdtchs->SetLineWidth(3);
    hmcpf->SetLineWidth(2);
    hmcchs->SetLineWidth(2);
    
    TLegend *leg2 = tdrLeg(0.17,0.50,0.32,0.90);
    leg2->SetTextSize(0.045*2);
    leg2->AddEntry(hdtpf,"PF DT","L");
    leg2->AddEntry(hdtchs,"CHS DT","L");
    leg2->AddEntry(hmcpf,"PF MC","L");
    leg2->AddEntry(hmcchs,"CHS MC","L");
    
    c1->SaveAs(Form("pdf/drawL1_RCL1.pdf"));
  }

}
