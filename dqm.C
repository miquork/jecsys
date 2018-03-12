// Purpose: plot PF composition distributions for multiple variables
// - (CHF, NHF, NEF) x (data, MC) x (direct, TP) for each pT, eta bin
// - same for multiplicities
#include "TFile.h"
#include "TH1D.h"

#include "tdrstyle_mod15.C"

void dqm() {

  TDirectory *curdir = gDirectory;
  setTDRStyle();

  const char *crun = "GH";
  TFile *fm = new TFile(Form("rootfiles/exclusion_cmsweek/mc/"
			     "%s/output-MC-1.root",crun));
  assert(fm && !fm->IsZombie());

  TFile *fd = new TFile(Form("rootfiles/exclusion_cmsweek/normal/"
			     "%s/output-DATA-1.root",crun));
  assert(fd && !fd->IsZombie());

  curdir->cd();

  //////////////////////
  // Plotter settings //
  //////////////////////

  // Choose fractions to be plotted
  const int nfrac = 3;//4;
  string fracs[nfrac] = {"chf", "nhf", "nef"};// "cef"};

  // Define colors and markers for each fraction
  map<string,int> colors;
  colors["chf"] = kRed;
  colors["nhf"] = kGreen+2;
  colors["nef"] = kBlue;
  colors["cef"] = kCyan+2;
  map<string, int> markers;
  markers["chf"] = kFullCircle;
  markers["nhf"] = kFullDiamond;
  markers["nef"] = kFullSquare;
  markers["cef"] = kFullStar;
  map<string, int> dpmarkers;
  dpmarkers["chf"] = kOpenCircle;
  dpmarkers["nhf"] = kOpenDiamond;
  dpmarkers["nef"] = kOpenSquare;
  dpmarkers["cef"] = kOpenStar;
  map<string, const char*> namesbb;
  namesbb["chf"] = "h^{#pm}";
  namesbb["nhf"] = "K^{0}_{had}";
  namesbb["nef"] = "#pi^{0}+K^{0}_{EM}";
  namesbb["cef"] = "e^{#pm}";
  map<string, const char*> namesec;
  namesec["chf"] = "h^{#pm}_{trk}";
  namesec["nhf"] = "h_{had}";
  namesec["nef"] = "#pi^{0}+h_{EM}";
  namesec["cef"] = "e^{#pm}";
  map<string,double> scalesbb;
  scalesbb["chf"] = 1.;
  scalesbb["nhf"] = 0.10;
  scalesbb["nef"] = 1.;
  scalesbb["cef"] = 0.01;
  map<string,double> scalesec;
  scalesec["chf"] = 0.05;
  scalesec["nhf"] = 1.0;
  scalesec["nef"] = 1.;
  scalesec["cef"] = 0.01;

  //////////////////////////////////////////////


  // Default bin, add loop later
  const char *cd = "Standard";
  //const char *ceta = "Eta_0.0-0.5";
  //const char *ctrg = "jt450";

  const int neta = 6;
  const char* etas[neta] = {"Eta_0.0-0.5","Eta_0.5-1.0","Eta_1.0-1.5",
			    "Eta_1.5-2.0","Eta_2.0-2.5","Eta_2.5-3.0"};
  //                        "Eta_3.2-4.7"};
  const int ntrg = 9;
  const char* trgs[ntrg] = {"jt40","jt60","jt80","jt140","jt200",
			    "jt260","jt320","jt400","jt450"};

  vector<string> filenames;
  for (int ieta = 0; ieta != neta; ++ieta) {
    for (int itrg = 0; itrg != ntrg; ++itrg) {


      const char *ceta = etas[ieta];
      const char *ctrg = trgs[itrg];


  // Find out the jet pT range for this trigger
  TH1D *hselpt = (TH1D*)fm->Get(Form("%s/%s/%s/hselpt",cd,ceta,ctrg));
  assert(hselpt);
  double minpt(0), maxpt(0);
  for (int i = 2; i != hselpt->GetNbinsX()+1; ++i) {
    if (hselpt->GetBinContent(i-1)==0 && hselpt->GetBinContent(i)!=0)
      if (minpt==0) minpt = hselpt->GetBinLowEdge(i);
    if (hselpt->GetBinContent(i)!=0 && hselpt->GetBinContent(i+1)==0)
      maxpt = hselpt->GetBinLowEdge(i+1);
  }
  float y1(0),y2(0);
  assert(sscanf(ceta,"Eta_%f-%f",&y1,&y2)==2);
  double y = 0.5*(y1+y2);

  // For loop, delete old histograms
  TObject *tob;
  tob = gROOT->FindObject("h"); if (tob) delete tob;
  tob = gROOT->FindObject("h2"); if (tob) delete tob;
  tob = gROOT->FindObject("c1"); if (tob) delete tob;

  // Create background histogram and canvas
  TH1D *h = new TH1D("h",";PF component fraction;dN/df (%)",101,0,1.01);
  h->SetMaximum(3.5);
  h->SetMinimum(0.+1e-4); // for tdrDiCanvas
  TH1D *h2 = new TH1D("h2",";PF component fraction;Data/MC",101,0,1.01);
  h2->SetMaximum(2.0-1e-5);//1.75-1e-5);
  h2->SetMinimum(0.5+1e-5);//0.75);

  lumi_13TeV = "Run2106GH re-mAOD 16.8 fb^{-1}";
  //TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);
  TCanvas *c1 = tdrDiCanvas("c1",h,h2,4,11);
		
  c1->cd(1);

  // Add basic labelling
  TLatex *tex = new TLatex();
  tex->SetTextSize(0.040);
  tex->SetNDC();
  if (y1==0) tex->DrawLatex(0.40,0.87,Form("|y|<%1.1f (%s)",y2,ctrg));
  if (y1!=0) tex->DrawLatex(0.40,0.87,Form("%1.1f<|y|<%1.1f (%s)",
					   y1,y2,ctrg));
  //tex->DrawLatex(0.40,0.82,"395<p_{T}<468 GeV");
  tex->DrawLatex(0.40,0.82,Form("%1.0f<p_{T}<%1.0f GeV",minpt,maxpt));
  
  TLegend *leg = tdrLeg(0.70,0.89-5*0.05,0.95,0.89);

  for (int i = 0; i != nfrac; ++i) {

    const char *cf = fracs[i].c_str();
    // Least biased tag-and-probe (TP) fractions sensitive to JEC
    TH1D *hd = (TH1D*)fd->Get(Form("%s/%s/%s/h%stp",cd,ceta,ctrg,cf));
    assert(hd);
    TH1D *hm = (TH1D*)fm->Get(Form("%s/%s/%s/h%stp",cd,ceta,ctrg,cf));
    assert(hm);

    // More biased direct probe (DP) that are sensitive to JER also
    TH1D *hd2 = (TH1D*)fd->Get(Form("%s/%s/%s/h%s",cd,ceta,ctrg,cf));
    assert(hd2);
    TH1D *hm2 = (TH1D*)fm->Get(Form("%s/%s/%s/h%s",cd,ceta,ctrg,cf));
    assert(hm2);

    // Normalize data and MC to unit area
    double scale = (y<2.5 ? scalesbb[cf] : scalesec[cf]);
    hd->Scale(100.*scale/hd->Integral());
    hm->Scale(100.*scale/hm->Integral());
    hd2->Scale(100.*scale/hd2->Integral());
    hm2->Scale(100.*scale/hm2->Integral());
    
    // calculate TP-DP difference for hatching
    for (int j = 1; j != hm2->GetNbinsX()+1; ++j) {
      double mean = 0.5*(hm->GetBinContent(j)+hm2->GetBinContent(j));
      double diff = 0.5*fabs(hm->GetBinContent(j)-hm2->GetBinContent(j));
      hm2->SetBinContent(j, mean);
      hm2->SetBinError(j, diff);
    }

    c1->cd(1);
    
    // Regular TP fractions
    int color = colors[cf];
    int marker = markers[cf];
    tdrDraw(hm,"HIST",kNone,0,kSolid,color,1001,color-9);
    hm->SetFillColorAlpha(color-9, 0.35); // 35% transparent
    tdrDraw(hd,"P",marker,color,kSolid,color,1001,color-9);
    if (marker!=kFullDiamond && marker!=kOpenDiamond) hd->SetMarkerSize(0.5);

    // Additional DP fractions
    int marker2 = dpmarkers[cf];
    tdrDraw(hm2,"E2",kNone,0,kDotted,color,1001,color-9);
    hm2->SetFillColorAlpha(color-9, 0.70); // 70% transparent
    tdrDraw(hd2,"P",marker2,color,kSolid,color,kNone,0);
    if (marker!=kFullDiamond && marker!=kOpenDiamond) hd2->SetMarkerSize(0.5);
        
    const char *name = (y<2.5 ? namesbb[cf] : namesec[cf]);
    if (scale==1) leg->AddEntry(hd,name,"FPL");
    if (scale!=1) leg->AddEntry(hd,Form("%s #times %1.2f",name,scale),"FPL");

    c1->cd(2);

    TLine *l = new TLine();
    l->SetLineStyle(kDotted);
    l->DrawLine(0,1,1,1);

    int nrb = 5;
    TH1D *hr = (TH1D*)hd->Clone(Form("r%s",cf));
    TH1D *hmrb = (TH1D*)hm->Clone(Form("mrb%s",cf));
    hr->Rebin(nrb);
    hmrb->Rebin(nrb);
    hr->Divide(hmrb);
    tdrDraw(hr,"r",marker,color,kSolid,color,1001,color-9);
    hr->SetMarkerSize(hr->GetMarkerSize()*2.);

    TH1D *hr2 = (TH1D*)hd2->Clone(Form("r%s2",cf));
    TH1D *hmrb2 = (TH1D*)hm2->Clone(Form("mrb%s2",cf));
    hr2->Rebin(nrb);
    hmrb2->Rebin(nrb);
    hr2->Divide(hmrb2);
    tdrDraw(hr2,"r",marker2,color,kSolid,color,1001,color-9);
    hr2->SetMarkerSize(hr2->GetMarkerSize()*2.);
  } // for i

  // Add labelling for TP and DP methods
  TH1D *hl1 = (TH1D*)h->Clone("hl1");
  hl1->SetMarkerStyle(kFullCircle);
  leg->AddEntry(hl1,"tag&probe","P");
  TH1D *hl2 = (TH1D*)h->Clone("hl2");
  hl2->SetMarkerStyle(kOpenCircle);
  leg->AddEntry(hl2,"probe only","P");
  
  gPad->RedrawAxis();

  string sfile = Form("dqm_y%02d-%02d_%s.pdf",
		      int(10.*y1),int(10*y2),ctrg);

  c1->cd(0);
  c1->SaveAs(Form("pdf/dqm/%s/%s",crun,sfile.c_str()));

  filenames.push_back(Form("\\includegraphics[width=0.30\\textwidth]{%s}",
			   sfile.c_str()));

    } // itrg
  } // ieta

  for (unsigned int i = 0; i != filenames.size(); ++i) {
    cout << filenames[i] << endl;
  }
} 
