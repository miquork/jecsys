// Purpose: Plot comparison of KIT and UH Z+jet results.
//          Discussion on Mattermost/CMS-JME-POG/Sync Task
//          Input files to JERCProtoLab with _SyncTask_[version]_ in the name
//          Success is declared when results match *exactly*
//          Things I think we should check vs pT,Z:
//          event counts, MPF, DB, mpf1, mpfn, mpfu, mZ.
#include "TFile.h"
#include "TLine.h"
#include "TLatex.h"
#include "TGraphErrors.h"

#include "../tdrstyle_mod15.C"

void drawSyncSet(string obs, string year, string lep, string data);
void drawSyncTask() {

  //drawSyncSet("counts","UL16nonAPV","mm","data");
  //drawSyncSet("mpf","UL16nonAPV","mm","data");

  drawSyncSet("counts","UL17B","mm","data");
  drawSyncSet("mpf","UL17B","mm","data");

} // drawSyncTask

void drawSyncSet(string obs, string year, string lep, string data) {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  // Hard-coded list of KIT file locations
  map<string, map<string, const char*> > mfk;
  mfk["UL16nonAPV"]["mm"] =
    "../JERCProtoLab/Summer19UL16/L3Residual_Z/nonAPV/JEC_Combination_Zmm/"
    "ZJetCombination_Zmm_DYJets_amc_21Feb2020_Summer19UL16_V7_L1L2L3Res.root";
  mfk["UL17B"]["mm"] =
    //"../JERCProtoLab/Summer19UL17/L3Residual_Z/JEC_Combination_Zmm/"
    //"ZJetCombination_Zmm_DYJets_Madgraph_09Aug2019_Summer19UL17_JECV6_L1L2L3Res.root";
    "rootfiles/ZJetCombination_Zmm_DYJets_amcatnlo_MiniAODv2_Summer19UL17_JECV6_L1L2L3Res.root";

  // Hard-coded list of UH file locations
  map<string, map<string, const char*> > mfh;
  mfh["UL16nonAPV"]["mm"] =
    "rootfiles/jme_bplusZ_merged_v37_2016FH_mu.root";
  mfh["UL17B"]["mm"] =
    "rootfiles/jme_bplusZ_merged_vX_2017B.root"; // mu, no MC

  // Open files
  TFile *fk = new TFile(mfk[year][lep]);
  assert(fk && !fk->IsZombie());
  TFile *fh = new TFile(mfh[year][lep]);
  assert(fh && !fh->IsZombie());

  // Load event counts to get binning
  TH1D *hk(0);
  if (year=="UL16nonAPV") hk = (TH1D*)fk->Get("Run2016postVFPFlateGH/Data_RawNEvents_CHS_a100_eta_00_13_L1L2L3Res");
  if (year=="UL17B") hk = (TH1D*)fk->Get("Run2017B/Data_RawNEvents_CHS_a100_eta_00_13_L1L2L3Res");
  assert(hk);
  TH1D *hh = (TH1D*)fh->Get("data/eta_00_13/statistics_rmpf_zmmjet_a100");
  assert(hh);
  
  // Load observables
  TGraphErrors *gk(0);
  if (year=="UL16nonAPV") gk = (TGraphErrors*)fk->Get("Run2016postVFPFlateGH/Data_MPF_CHS_a100_eta_00_13_L1L2L3Res");
  if (year=="UL17B") gk = (TGraphErrors*)fk->Get("Run2017B/Data_MPF_CHS_a100_eta_00_13_L1L2L3Res");
  assert(gk);
  TGraphErrors *gh = (TGraphErrors*)fh->Get("data/eta_00_13/rmpf_zmmjet_a100");
  assert(gh);

  // Setup plotting, initially for event counts
  TH1D *hu(0), *hd(0);
  if (obs=="counts") {
    hu = tdrHist("hu","Events per bin",0.5,2e6,"p_{T,Z} (GeV)",15,1500);
    hd = tdrHist("hd","UH/KIT-1 (%)",-20,+30,"p_{T,Z} (GeV)",15,1500);
  }
  if (obs=="mpf") {
      hu = tdrHist("hu","MPF",0.95,1.25,"p_{T,Z} (GeV)",15,1500);
      hd = tdrHist("hd","UH/KIT-1 (%)",-1.5,+1.5,"p_{T,Z} (GeV)",15,1500);
  }
  lumi_13TeV = "2016GH #mu#mu data"; // GENERALIZE
  TCanvas *c1 = tdrDiCanvas("c1",hu,hd,4,11);

  c1->cd(1);
  gPad->SetLogx();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(15,1,1500,1);

  TLegend *leg = tdrLeg(0.68,0.90-3*0.06,0.98,0.90);
  leg->SetHeader("L1L2L3Res");
  if (obs=="counts") {
    gPad->SetLogy();
    tdrDraw(hk,"HIST",kNone,kBlue,kSolid,kBlue,1001,kBlue-9);
    tdrDraw(hh,"Pz",kFullCircle,kRed,kSolid,kRed,kNone);
    leg->AddEntry(hh,"UH Z+jet","PLE");
    leg->AddEntry(hk,"KIT Z+jet","F");
  }
  else {
    tdrDraw(gk,"PLz",kFullSquare,kBlue,kSolid,kBlue,kNone);
    tdrDraw(gh,"Pz",kFullCircle,kRed,kSolid,kRed,kNone);
    leg->AddEntry(gh,"UH Z+jet","PE");
    leg->AddEntry(gk,"KIT Z+jet","PLE");
  }

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  tex->DrawLatex(0.42,0.855,"|#eta|<1.3, #alpha<1.0");

  gPad->RedrawAxis();

  c1->cd(2);
  gPad->SetLogx();
  l->DrawLine(15,0,1500,0);

  //TH1D *hr = (TH1D*)hh->Clone("hr");
  //hr->Add(hk,-1);
  //hr->Divide(hk);
  //hr->Scale(100.);

  // Divide by hand due to different number of bins
  TH1D *hr = (TH1D*)hk->Clone("hr"); hr->Clear();
  if (obs=="counts") {
    for (int i = 1; i != hr->GetNbinsX()+1; ++i) {
      int j = hh->GetXaxis()->FindBin(hr->GetBinCenter(i));
      hr->SetBinContent(i, 100.*(hh->GetBinContent(j)/hk->GetBinContent(i)-1));
      hr->SetBinError(i, 100.*hk->GetBinError(i)/hk->GetBinContent(i));
    }
  }
  else {
    for (int i = 0; i != gk->GetN()+1; ++i) {
      int j = hr->GetXaxis()->FindBin(gk->GetX()[i]);
      hr->SetBinContent(j, gk->GetY()[i]);
      hr->SetBinError(j, 100. * gk->GetEY()[i] / gk->GetY()[i]);
    } // for i
    for (int i = 0; i != gh->GetN()+1; ++i) {
      int j = hr->GetXaxis()->FindBin(gh->GetX()[i]);
      hr->SetBinContent(j, 100.*(gh->GetY()[i]/hr->GetBinContent(j)-1));
    } // for i
  }

  tdrDraw(hr,"Pz",kFullCircle,kRed);

  c1->SaveAs(Form("pdf/drawSyncTask/drawSyncTask_%s_%s_%s_%s.pdf",
		  obs.c_str(), year.c_str(), lep.c_str(), data.c_str()));
} // void drawSyncSet
