// Purpose: Plot comparison of multijet flavor for leading jet and recoil
//          Use pT,ave binning and JetFlavor, which contains heavy hadrons as bc
//          (that is g>bbbar and g>cccbar are classified as b,c and not g)
//          and has no 'none' flavored jets.
//          Parameterize these fractions to be used in Flavor.C/h
// Credits to Minsuk Kim for producing the inputs
//
// Run with 'root -l -b -q minitools/drawMultijetFlavor.C+g'
#include "TFile.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TF1.h"
#include "TMultiGraph.h"

#include "../tdrstyle_mod15.C"

// NB: HF fractions (b,c) look consistent between recoil and lead,
//     and also s fraction is very similar (recoil<leading) => use multigraphs?
// Q: why is s so similar to b fraction?

void cleanGraph(TGraph *g) {
  for (int i = g->GetN()-1; i != -1; --i) {
    if (g->GetY()[i]==0) g->RemovePoint(i);
  }
}

void drawMultijetFlavor() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TFile *f = new TFile("rootfiles/multijet_flavour_v2.root","READ");
  assert(f && !f->IsZombie());
  f->cd("MC/Pt30/flavour");
  TDirectory *d = gDirectory;

  curdir->cd();

  const int nf = 5;
  string flavor[nf] = {"ud","s","c","b","g"};
  int color[nf] = {kBlue,kCyan+2,kGreen+2,kRed,kOrange+2};
  int fmarker[nf] = {kFullSquare,kFullStar,kFullDiamond,kFullDiamond,
		     kFullCircle};
  int omarker[nf] = {kOpenSquare,kOpenStar,kOpenDiamond,kOpenDiamond,
		     kOpenCircle};
  const int nt = 2;
  string type[nt] = {"recoiljet","leadingjet"};
  
  lumi_13TeV = "Multijet UL17+UL18 MC";
  TH1D *h = tdrHist("h","Flavor fraction",0,1.1,
		    "p_{T,ave} (GeV)",40,2000);
  TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);
  gPad->SetLogx();

  TLegend *leg1 = tdrLeg(0.37,0.90-5*0.045,0.57,0.90);
  leg1->SetTextSize(0.04);
  leg1->SetHeader("Recoil");
  TLegend *leg2 = tdrLeg(0.60,0.90-5*0.045,0.80,0.90);
  leg2->SetTextSize(0.04);
  leg2->SetHeader("Leading jet");

  cout << " // Parameters from minitools/drawMultijetFlavor.C" << endl;

  double p1s(0), p1b(0); // fix slopes for b,s
  for (int i = 0; i != nt; ++i) {

    if (i==0) cout << "  double pmjrecoil[nf][3] = {" << endl;
    if (i==1) cout << "  double pmjlead[nf][3] = {" << endl;

    for (int j = 0; j != nf; ++j) {
      
      const char *ct = type[i].c_str();
      const char *cf = flavor[j].c_str();

      TGraphErrors *g = (TGraphErrors*)d->Get(Form("ptave_%s_%s_ul17",ct,cf));
      assert(g);
      TGraphErrors *g2 = (TGraphErrors*)d->Get(Form("ptave_%s_%s_ul18",ct,cf));
      assert(g2);
      
      tdrDraw(g,"Pz",i==0 ? fmarker[j] : omarker[j], color[j]);
      //tdrDraw(g2,"Pz",i==0 ? fmarker[j] : omarker[j], color[j]);

      // Remove empty points to avoid bias in fits
      cleanGraph(g);
      cleanGraph(g2);

      g2->SetMarkerSize(0.7);

      TMultiGraph *mg = new TMultiGraph();
      mg->Add(g);
      mg->Add(g2);

      TF1 *f1 = new TF1(Form("f1_%s_%s",ct,cf),"[0]+[1]*pow(x,[2])",
			//40,2000);
			40,1650.);
      //if (j==2) f1->SetRange(40,1450.); // charm drop-off at high pT,ave
      if (j==1||j==2||j==3) f1->SetRange(40,1450.); // charm drop-off at high pT,ave
      //if (i==0 && (j==0||j==4))
      //f1->SetRange(40,1650.); // ud drop off at high pT,ave
      f1->SetLineColor(color[j]);
      f1->SetLineStyle(i==0 ? kSolid : kDashed);
      if (j==4) f1->SetParameters(0.84,-0.028,0.3333); //gluons
      if (j==0) f1->SetParameters(0.04,+0.021,0.3333); //ud
      if (j==2) f1->SetParameters(0.04,+0.007,0.3333); //c
      //if (j==1) f1->SetParameters(0.04,+0.000,0.3333); //s
      //if (j==3) f1->SetParameters(0.04,+0.000,0.3333); //b
      if (j==1) f1->SetParameters(0.06,-0.003,0.3333); //s
      if (j==3) f1->SetParameters(0.02,+0.003,0.3333); //b
      // NB: initial values sum up to exactly 1.000 everywhere
      // Would need multifit to obey the same boundary after fit
      if (!(j==0||j==4)) f1->FixParameter(2,0.3333);
      if (i==1 && (j==1)) f1->FixParameter(1, p1s);
      if (i==1 && (j==3)) f1->FixParameter(1, p1b);
      //if (j!=2) 
      //g->Fit(f1,"QRN");
      mg->Fit(f1,"QRN");
      f1->Draw("SAME");

      if (i==0 && (j==1)) p1s = f1->GetParameter(1);
      if (i==0 && (j==3)) p1b = f1->GetParameter(1);

      if (i==0) leg1->AddEntry(g,Form("(%1.1f/%d)",
				      f1->GetChisquare(),f1->GetNDF()),"PLE");
      if (i==1) leg2->AddEntry(g,Form("%s-jet (%1.1f/%d)",cf,
				      f1->GetChisquare(),f1->GetNDF()),"PLE");

      if (j==0) cout << "   // " << f1->GetExpFormula() << endl;

      cout << Form("    {%1.5g, %1.5g, %1.5g}",
		   f1->GetParameter(0),f1->GetParameter(1),f1->GetParameter(2));
      cout << Form("%s // (%s %1.1f/%d)",(j==4 ? "};" : ","),
		   cf,f1->GetChisquare(), f1->GetNDF()) << endl;
    } // for j
  } // for i
  
  c1->SaveAs("pdf/drawMultijetFlavor.pdf");
} // drawMultijetFlavor
