// Purpose: Draw new JES uncertainty sources from global fit
//          See also code to plot correlations globalFitPulls.C
//          Code to prooduce uncertainty files: drawJetCorrectionUncertainty.C
//          Code to plot individual sources: test/mk_drawUncertainty.C
//          Code to redo plots from source .txt files: test/...
#include "TFile.h"
#include "TSystem.h"

#include "../tdrstyle_mod15.C"

void drawJESunc() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  //TFile *f = new TFile("rootfiles/jecdata2018ABCD.root");
  //TFile *f = new TFile("rootfiles/jecdataBCDEF.root");
  TFile *f = new TFile("rootfiles/jecdata2016GH.root");
  assert(f && !f->IsZombie());
  
  curdir->cd();

  TH1D *hjes = (TH1D*)f->Get("ratio/eta00-13/sys/hjesfit"); assert(hjes);

  int neig(0);
  TH1D *heig(0);
  vector<TH1D*> vhe;
  while ((heig=(TH1D*)f->Get(Form("ratio/eta00-13/sys/hjesfit_eig%d",neig)))) {
    vhe.push_back(heig);
    ++neig;
  }
  assert(vhe.size()!=0);

  int npar(0);
  TH1D *hpar(0);
  vector<TH1D*> vhp;
  while ((hpar=(TH1D*)f->Get(Form("ratio/eta00-13/sys/hjesfit_par%d",npar)))) {
    vhp.push_back(hpar);
    ++npar;
  }
  assert(vhp.size()==vhe.size());

  // Turn hjesfit into uncertainty band in %'s around zero
  TH1D *herr = (TH1D*)hjes->Clone("herr");
  for (int i = 1; i != hjes->GetNbinsX()+1; ++i) {
    herr->SetBinError(i, 100.*hjes->GetBinError(i));///hjes->GetBinContent(i));
    herr->SetBinContent(i, 0);
    hjes->SetBinContent(i, 100.*(hjes->GetBinContent(i)-1));
    hjes->SetBinError(i, 100.*hjes->GetBinError(i));
  }

  double ptmin = 15;
  double ptmax = 4000;//4500;
  TH1D *h = tdrHist("h","Absolute JES uncertainty (%)",
		    -0.8,1.0,"p_{T} (GeV)",ptmin,ptmax);
  lumi_13TeV = "UL16GH";
  //lumi_13TeV = "UL17";
  //lumi_13TeV = "UL18";
  TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);
  gPad->SetLogx();

  // (1+6+1)*0.035=0.28
  TLegend *leg = tdrLeg(0.60,0.61,0.90,0.89);
  leg->SetTextSize(0.035);
  TLine *l = new TLine(); l->SetLineStyle(kDashed);

  tdrDraw(herr,"E3",kNone,kYellow+2,kSolid,-1,1001,kYellow+1);
  l->DrawLine(ptmin,0,ptmax,0);

  leg->AddEntry(herr,"Total","F");

  int color[] = {kBlack,kRed,kGreen+2,kBlue,kOrange+2,kMagenta+2,kCyan+2,
		 kViolet,kPink,kGray,kGray,kGray};

  TH1D *sumhe = (TH1D*)vhe[0]->Clone("sumhe"); sumhe->Reset();
  for (int i = 0; i != vhe.size(); ++i) {

    TH1D *he = vhe[i];
    he->Scale(100.);
    tdrDraw(he,"L",kNone,color[i],kSolid,-1,kNone,-1);

    //if (i<6) {
    if (i<7) {
      for (int j = 1; j != sumhe->GetNbinsX()+1; ++j) {
	sumhe->SetBinContent(j, sumhe->GetBinContent(j) +
			     pow(he->GetBinContent(j),2));
      } // for j
      leg->AddEntry(he,Form("Source %d",i+1),"L");
    }
  } // for i

  for (int j = 1; j != sumhe->GetNbinsX()+1; ++j) {
    sumhe->SetBinContent(j, sqrt(sumhe->GetBinContent(j)));
  } // for j

  tdrDraw(sumhe,"L",kNone,kYellow+3,kSolid,-1,kNone,-1);
  sumhe->SetLineWidth(2);
  //leg->AddEntry(sumhe,"Top 6","L");

  gPad->RedrawAxis();
  //c1->SaveAs("pdf/drawJESunc_err_UL18.pdf");
  //c1->SaveAs("pdf/drawJESunc_err_UL17.pdf");
  c1->SaveAs("pdf/drawJESunc_err_UL16GH.pdf");

  /////////////////////////////////////////////////////////

  TH1D *h2 = tdrHist("h2","JES components",
		     //0.97,1.02,
		     //-3,+2,
		     //-2.5,+3, // UL17
		     -3.5,+3.5, // UL16
		     "p_{T} (GeV)",ptmin,ptmax);
  TCanvas *c2 = tdrCanvas("c2",h2,4,0,kSquare);
  gPad->SetLogx();
  
  //TLegend *leg2 = tdrLeg(0.30,0.61,0.90,0.89);
  // 0.04*(1+8+1)=0.40, divided in two
  TLegend *leg2 = tdrLeg(0.25,0.70,0.85,0.90);
  leg2->SetNColumns(2);
  leg2->SetTextSize(0.035);

  tdrDraw(hjes,"E3",kNone,kYellow+2,kSolid,-1,1001,kYellow+1);
  l->DrawLine(ptmin,0,ptmax,0);

  leg2->AddEntry(hjes,"Total JES","FL");

  string names[] = {"Tracker","Photon","Hadron","HCAL","ECAL",
		    "Herwig","L1RC","Trk data","#sigma_{MB}"};
  
  TH1D *sumhp = (TH1D*)vhp[0]->Clone("sumhp"); sumhp->Reset();
  for (int i = 0; i != vhp.size(); ++i) {

    TH1D *hp = vhp[i];
    // Turn into percentage difference for easier plottig
    for (int j = 1; j != hp->GetNbinsX()+1; ++j) {
      hp->SetBinContent(j, 100.*(hp->GetBinContent(j)-1));
      hp->SetBinError(j, 100.*hp->GetBinError(j));
    }

    tdrDraw(hp,"LE3",kNone,color[i],kSolid,-1,1001,
	    color[i]==kBlack? kGray : color[i]-9);
    hp->SetFillColorAlpha(color[i]==kBlack ? kGray : color[i]-9, 0.7);

    TGraph *gk = new TGraph(hp);
    gk->SetLineColor(color[i]);
    gk->Draw("SAMEL");

    //if (i<6) {
    for (int j = 1; j != sumhp->GetNbinsX()+1; ++j) {
      sumhp->SetBinContent(j, sumhp->GetBinContent(j) +
			   hp->GetBinContent(j));
      //(hp->GetBinContent(j)-1));

      /*
      // Cumulative sums
      for (int k = 0; k != vhs.size(); ++k) {
	int k2 = nhs[k];
	if (vhs[k2]==0) {
	  vhs[k2] = (TH1D*)sumhp->Clone(Form("sumhp%d",k));
	  vhs[k2]->Reset();
	} 
	TH1D *hs = vhs[k2];
	const double r = 0.5; // reduce size of uncertainty band
	//if (i<=k2) {
	
	  hs->SetBinContent(j, hs->GetBinContent(j) + hp->GetBinContent(j));
	  hs->SetBinError(j, min(max(hs->GetBinError(j),r*hp->GetBinError(j)),
				 r*hjes->GetBinError(j)));
	}
      } // for k
      */
    } // for j
    //leg2->AddEntry(hp,Form("Parameter %d",i+1),"FL");
    leg2->AddEntry(hp,Form("%s (p%d)",names[i].c_str(),i),"FL");
    //}

  } // for i

  //for (int j = 1; j != sumhp->GetNbinsX()+1; ++j) {
  //sumhp->SetBinContent(j, sumhp->GetBinContent(j));
  //} // for j

  tdrDraw(sumhp,"L",kNone,kYellow+3,kSolid,-1,kNone,-1);
  sumhp->SetLineWidth(2);
  leg2->AddEntry(sumhp,"Sum of pars","L");

  gPad->RedrawAxis();
  //c2->SaveAs("pdf/drawJESunc_par_UL18.pdf");
  //c2->SaveAs("pdf/drawJESunc_par_UL17.pdf");
  c2->SaveAs("pdf/drawJESunc_par_UL16GH.pdf");

  ///////////////////////////////////////////////////////////
  gSystem->Exec("rm pdf/drawJESunc_par2.gif");
  vector<TH1D*> vhs(vhp.size()); // cumulative shifts
  // ordering of cumulative sum
  //int nhs[] = {1,0,2,7,3,4,6,5};
  //int nhs[] = {1,0,7,2,3,4,6,5}; // v1
  //int nhs[] = {1,2,0,7,3,4,6,5};
  //int nhs[] = {1,2,7,0,3,4,6,5};
  //int nhs[] = {5, 1,0, 7, 2,3,4, 6}; // v2
  int nhs[] = {5, 1,2,0, 7, 3,4, 6, 8}; // v2
  for (int k = 0; k != vhs.size(); ++k) {

    TH1D *hs = (TH1D*)sumhp->Clone(Form("sumhp%d",k));
    hs->Reset();

    for (int i = 0; i != k+1; ++i) {
      TH1D *hp = vhp[nhs[i]];
      
      for (int j = 1; j != hs->GetNbinsX()+1; ++j) {
	hs->SetBinContent(j, hs->GetBinContent(j)+hp->GetBinContent(j));
	//const double r = 0.5; // reduce size of uncertainty band
	//hs->SetBinError(j, min(max(hs->GetBinError(j),r*hp->GetBinError(j)),
	//		       r*hjes->GetBinError(j)));
	hs->SetBinError(j, hjes->GetBinError(j)*0.5*(1-0.1*k));
      } // for j
    } // for i
    vhs[k] = hs;
  } // for k

  TH1D *h3 = tdrHist("h3","Cumulative JES components",
		     //-3.5,+2, // UL17
		     -3.5,+2, // UL16
		     "p_{T} (GeV)",ptmin,ptmax);
  TCanvas *c3 = tdrCanvas("c3",h3,4,0,kSquare);
  gPad->SetLogx();

  TLegend *leg3 = tdrLeg(0.25,0.70,0.85,0.90);
  leg3->SetNColumns(2);
  leg3->SetTextSize(0.035);
  
  tdrDraw(hjes,"E3",kNone,kYellow+2,kSolid,-1,1001,kYellow+1);
  l->DrawLine(ptmin,0,ptmax,0);
  leg3->AddEntry(hjes,"Total JES","FL");

  for (int i = 0; i != vhs.size(); ++i) {

    int k2 = nhs[i];
    TH1D *hs = vhs[i];
    tdrDraw(hs,"LE3",kNone,color[k2],kSolid,-1,1001,
	    color[k2]==kBlack? kGray : color[k2]-9);
    hs->SetFillColorAlpha(color[k2]==kBlack ? kGray : color[k2]-9, 0.7);

    TGraph *gs = new TGraph(hs);
    gs->SetLineColor(color[k2]);
    gs->Draw("SAMEL");

    leg3->AddEntry(hs,Form("%s (p%d)",names[k2].c_str(),k2),"FL");

    c3->SaveAs("pdf/drawJESunc_par2.gif+50");
  }
  //c3->SaveAs("pdf/drawJESunc_par2_UL18.pdf");
  //c3->SaveAs("pdf/drawJESunc_par2_UL17.pdf");
  c3->SaveAs("pdf/drawJESunc_par2_UL16GH.pdf");

} // drawJESunc
