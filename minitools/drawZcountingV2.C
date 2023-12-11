// Purpose: Draw jet counting vs Z counting per fill, then vs lumi
// Author:  mikko (dot) voutilainen (at) cern (dot) ch
// Data:    2023-08-15
#include "TFile.h"
#include "TTree.h"
#include "TProfile.h"
#include "TF1.h"

#include "../tdrstyle_mod22.C"
double minLumi = 500;//1000 pb^{-1}

void drawZcountingV2() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TFile *f = new TFile("rootfiles/zcountingV2.root","READ");
  assert(f && !f->IsZombie());
  //TFile *f2 = new TFile("rootfiles/zcountingV2_ecalPrefire_v1.root","READ");
  TFile *f2 = new TFile("rootfiles/zcountingV2_emuPrefires_v2.root","READ");
  assert(f2 && !f2->IsZombie());

  TTree *t = (TTree*)f->Get("counts");
  TTree *t2 = (TTree*)f2->Get("counts");
  curdir->cd();

  t->Draw("jetcount/reclumj:run>>pj(60000,270000,330000)","(reclumj>0)*reclumj","prof");
  t2->Draw("jetcount/reclumj:run>>pjj(60000,270000,330000)","(reclumj>0)*reclumj","prof");
  t2->Draw("zcount/reclumz:run>>pz(60000,270000,330000)","(reclumz>0)*reclumz","prof");
  t2->Draw("run>>hzlum(60000,270000,330000)","reclumz","");
  
  TProfile *pj = (TProfile*)gROOT->FindObject("pj");
  TProfile *pjj = (TProfile*)gROOT->FindObject("pjj");
  TProfile *pz = (TProfile*)gROOT->FindObject("pz");
  TH1D *hzlum = (TProfile*)gROOT->FindObject("hzlum");

  TH1D *hj = pj->ProjectionX("hj");
  TH1D *hjj = pjj->ProjectionX("hjj");
  TH1D *hz = pz->ProjectionX("hz");
  
  TF1 *f1 = new TF1("f1","[0]",0,500000);
  hj->Fit(f1,"QRNW");
  hj->Scale(1./f1->GetParameter(0));
  hjj->Fit(f1,"QRNW");
  hjj->Scale(1./f1->GetParameter(0));
  hz->Fit(f1,"QRNW");
  hz->Scale(1./f1->GetParameter(0));

  // Create helper histogram for checking IOV boundaries
  //double viov[] = {270000, 290000, 310000, 330000};
  //double viov[] = {273158,278777, 278769,284044, 297050,306460, 315257,325172};
  double viov[] = {273158,278777, 278778,284044, 297050,306460, 315257,325172};
  const int niov = sizeof(viov)/sizeof(viov[0])-1;
  TH1D *hiov = new TH1D("hiov","",niov,viov);

  // From David Walter:
  // https://github.com/CMS-LUMI-POG/ZCounting/blob/Run2/ZHarvester/Plotting/plot_jet_comparison.py
  map<int,string> rundict;
  rundict[273158] = "16BCD";
  rundict[276831] = "EF";//"16EF";
  rundict[278820] = "16GH";
  rundict[297020] = "17B";          
  rundict[299337] = "C";//17C";          
  rundict[302030] = "D";//"17D";          
  rundict[303435] = "E";//"17E";          
  rundict[304911] = "F";//"17F";          
  rundict[315252] = "18A";  
  rundict[316998] = "B";//"18B";  
  rundict[319313] = "C";//"18C";  
  rundict[320394] = "D";//"18D";
  rundict[325172] = "end"; // the end
  
  // Determine new binning to have >1/fb per bin
  // TBD: check if IOV boundary to not go across those
  vector<int> ix;
  vector<double> x;
  for (int i = 1; i != hzlum->GetNbinsX()+1; ++i) {

    // Sum up active runs
    if (hzlum->GetBinContent(i)>0) {

      int thisRun = hzlum->GetBinLowEdge(i);
      int nextRun = hzlum->GetBinLowEdge(i+1);
      
      // Set starting point to first active run
      if (x.size()==0) {
	x.push_back(thisRun);
	ix.push_back(i);
      }

      // Check if latest bins add up to >1/fb, or hitting IOV boundary
      bool isBoundary = (hiov->FindBin(thisRun) != hiov->FindBin(nextRun));
      int ix1 = ix[ix.size()-1];
      double sumLumi = hzlum->Integral(ix1,i); // in /pb
      if (sumLumi>minLumi || isBoundary) {

	// If boundary and not full 1/fb, replace previous upper edge
	if (isBoundary && sumLumi<minLumi) {
	  x[x.size()-1] = nextRun;
	  ix[ix.size()-1] = i+1;
	}
	// Otherwise just start next bin
	else {
	  x.push_back(nextRun);
	  ix.push_back(i+1);
	}
      }
    }
  } // hzlum
  cout << "x bins size " << x.size() << endl << flush;
  cout << "ix bins size " << ix.size() << endl << flush;
  
  // Fill integrated luminosity in each bins
  TH1D *hsumlum = new TH1D("hsumlum","",x.size()-1,&x[0]);
  for (int i = 1; i != hzlum->GetNbinsX()+1; ++i) {
    hsumlum->Fill(hzlum->GetBinCenter(i), hzlum->GetBinContent(i)*0.001);
  }
  // Fill cumulative luminosity per bin, and create new binning
  TH1D *hcumlum = new TH1D("hcumlum","",x.size()-1,&x[0]);
  vector<double> vlum; vlum.push_back(0);
  for (int i = 1; i != hcumlum->GetNbinsX()+1; ++i) {
    hcumlum->SetBinContent(i, hcumlum->GetBinContent(i-1) +
			   hsumlum->GetBinContent(i));
    vlum.push_back(hcumlum->GetBinContent(i));
  }

  // New histograms of jet/Z xsec vs cumulative luminosity;
  TH1D *hz2 = new TH1D("hz2",";Cumulative luminosity (fb^{-1});"
		       "#sigma / #LT#sigma#GT",vlum.size()-1,&vlum[0]);
  TH1D *hj2 = new TH1D("hj2",";Cumulative luminosity (fb^{-1});"
		       "#sigma / #LT#sigma#GT",vlum.size()-1,&vlum[0]);
  TH1D *hjj2 = new TH1D("hjh2",";Cumulative luminosity (fb^{-1});"
			"#sigma / #LT#sigma#GT",vlum.size()-1,&vlum[0]);
  for (int i = 1; i != hz->GetNbinsX()+1; ++i) {
    int run = hz->GetBinLowEdge(i);
    double lum = hzlum->GetBinContent(i)*0.001;
    int lumbin = hsumlum->FindBin(run);
    double sumlum = hsumlum->GetBinContent(lumbin);
    double cumlum = hcumlum->GetBinContent(lumbin);
    double frac = (sumlum ? lum / sumlum : 0.);
    hz2->Fill(cumlum, hz->GetBinContent(i)*frac);
    hz2->SetBinError(lumbin, 0.);

    hj2->Fill(cumlum, hj->GetBinContent(i)*frac);
    hj2->SetBinError(lumbin, 0.);
    hjj2->Fill(cumlum, hjj->GetBinContent(i)*frac);
    hjj2->SetBinError(lumbin, 0.);
  }
  
  //double lum16g = hcumlum->GetBinContent(hcumlum->FindBin(viov[1]));
  //double lum17 = hcumlum->GetBinContent(hcumlum->FindBin(viov[3]));
  //double lum18 = hcumlum->GetBinContent(hcumlum->FindBin(viov[5]));

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  TLatex *tex = new TLatex();
  //tex->SetNDC(); 
  tex->SetTextSize(0.045);

  double totlum = hz2->GetBinLowEdge(hz2->GetNbinsX()+1);;
  lumi_13TeV = Form("Run 2, %1.1f fb^{-1}",totlum);
  extraText = "Work-in-progress";
		    
  TH1D *h1 = tdrHist("h1","#sigma / #LT#sigma#GT",0.90,1.10,
		     "Cumulative luminosity (fb^{-1})",0,135);
  TCanvas *c1 = tdrCanvas("c1",h1,4,11,kRectangular);

  // Draw grid of IOVs
  l->SetLineStyle(kSolid);
  l->DrawLine(0,1,135,1);
  for (map<int,string>::const_iterator it = rundict.begin();
       it != rundict.end(); ++it) {
    int run = (it->first);
    string siov = (it->second);
    if (siov=="end") break;
    double lum = hcumlum->GetBinContent(hcumlum->FindBin(run));
    if (siov=="16BCD" || siov=="16GH" || siov=="17B" || siov=="18A") {
      l->SetLineStyle(kDashed);
      l->DrawLine(lum,0.90,lum,1.10);
      tex->SetTextSize(0.045);
      tex->DrawLatex(lum+1,0.91,siov.c_str());
    }
    else {
      l->SetLineStyle(kDotted);
      l->DrawLine(lum,0.90,lum,1.10);
      tex->SetTextSize(0.035);
      tex->DrawLatex(lum+0.5,0.922,siov.c_str());
    }
  }

  for (int i = 0; i != niov; ++i) {
    if (i%2==1) continue;
    double run1 = viov[i];
    double run2 = viov[i+1];
    double lum1 = hcumlum->GetBinContent(hcumlum->FindBin(run1));
    double lum2 = max(hcumlum->GetBinContent(hcumlum->FindBin(run2)),
		      hcumlum->GetBinContent(hcumlum->FindBin(run2)-1));
    cout << "  run1 / lum1 : " << run1 << " / " << lum1
	 << ", run2 / lum2 : " << run2 << " / " << lum2 << endl << flush;
    f1->SetRange(lum1,lum2);

    hz2->Fit(f1,"QRNW");
    double p0z = f1->GetParameter(0);
    l->SetLineStyle(kSolid);
    l->SetLineColor(kGreen+2);
    l->DrawLine(lum1,p0z,lum2,p0z);

    hjj2->Fit(f1,"QRNW");
    double p0 = f1->GetParameter(0);
    l->SetLineColor(kRed);
    l->DrawLine(lum1,p0,lum2,p0);

    hj2->Fit(f1,"QRNW");
    double p0j = f1->GetParameter(0);
    l->SetLineStyle(kDotted);
    //l->DrawLine(lum1,p0j,lum2,p0j);
  }

  tdrDraw(hz2,"Pz",kFullSquare,kGreen+2);
  //tdrDraw(hj2,"Pz",kOpenCircle,kRed);
  tdrDraw(hjj2,"Pz",kFullCircle,kRed);

  hz2->SetMarkerSize(0.8);
  hj2->SetMarkerSize(0.8);
  hjj2->SetMarkerSize(0.8);

  TLegend *leg = tdrLeg(0.60,0.90,0.80,0.90);
  leg->AddEntry(hz2,"Z boson","PLE");
  leg->AddEntry(hjj2,"Jet","PLE");
  //leg->AddEntry(hj2,"Jet (w/o ECAL prefire)","PLE");
  leg->SetY1(leg->GetY2()-0.05*leg->GetNRows());
  
  tex->SetNDC();
  tex->SetTextSize(0.040);
  tex->DrawLatex(0.62,leg->GetY1()-0.05,
		 "|#eta_{jet}|<1.3, p_{T,jet}>592 GeV");
  tex->SetTextSize(0.035);
  tex->DrawLatex(0.60,0.30,"ECAL prefire 16+17 (per year)");
  tex->DrawLatex(0.60,0.26,"Jet veto map (per year)");

  gPad->RedrawAxis();
  c1->SaveAs("pdf/drawZcountingV2.pdf");


  // Ratio of jet and Z counting
  TH1D *h2 = tdrHist("h2","#sigma_{jet} / #LT#sigma_{jet}#GT)"
		     " / (#sigma_{Z} / #LT#sigma_{Z}#GT)",0.90,1.10,
		     "Cumulative luminosity (fb^{-1})",0,135);
  TCanvas *c2 = tdrCanvas("c2",h2,4,11,kRectangular);

  TH1D *hr2 = (TH1D*)hjj2->Clone("hr2");
  hr2->Divide(hz2);

  // Draw grid of IOVs
  l->SetLineColor(kBlack);
  l->SetLineStyle(kSolid);
  l->DrawLine(0,1,135,1);
  tex->SetNDC(kFALSE);
  for (map<int,string>::const_iterator it = rundict.begin();
       it != rundict.end(); ++it) {
    int run = (it->first);
    string siov = (it->second);
    if (siov=="end") break;
    double lum = hcumlum->GetBinContent(hcumlum->FindBin(run));
    l->SetLineColor(kBlack);
    if (siov=="16BCD" || siov=="16GH" || siov=="17B" || siov=="18A") {
      l->SetLineStyle(kDashed);
      l->DrawLine(lum,0.90,lum,1.10);
      tex->SetTextSize(0.045);
      tex->DrawLatex(lum+1,0.91,siov.c_str());
    }
    else {
      l->SetLineStyle(kDotted);
      l->DrawLine(lum,0.90,lum,1.10);
      tex->SetTextSize(0.035);
      tex->DrawLatex(lum+0.5,0.922,siov.c_str());
    }

    map<int,string>::const_iterator jt = it; ++jt; // next run
    double run2 = jt->first; // std::next in C++0x?
    double lum2 = max(hcumlum->GetBinContent(hcumlum->FindBin(run2)),
		      hcumlum->GetBinContent(hcumlum->FindBin(run2)-1));
    cout << "  run / lum : " << run << " / " << lum
    	 << ", run2 / lum2 : " << run2 << " / " << lum2 << endl << flush;
    f1->SetRange(lum,lum2);

    hr2->Fit(f1,"QRNW");
    double p0 = f1->GetParameter(0);
    l->SetLineStyle(kDotted);
    l->SetLineColor(kRed);
    l->SetLineWidth(2);
    l->DrawLine(lum,p0,lum2,p0);
    l->SetLineWidth(1);
  }
  
  tdrDraw(hr2,"Pz",kFullCircle,kRed);
  hr2->SetMarkerSize(0.8);

  TLegend *leg2 = tdrLeg(0.60,0.90,0.80,0.90);
  leg2->AddEntry(hr2,"Jet / Z boson","PLE");
  leg2->SetY1(leg2->GetY2()-0.05*leg2->GetNRows());
  
  tex->SetNDC();
  tex->SetTextSize(0.040);
  tex->DrawLatex(0.62,leg->GetY1()-0.05,
		 "|#eta_{jet}|<1.3, p_{T,jet}>592 GeV");
  tex->SetTextSize(0.035);
  tex->DrawLatex(0.60,0.30,"ECAL prefire 16+17 (per year)");
  tex->DrawLatex(0.60,0.26,"Jet veto map (per year)");

  for (int i = 0; i != niov; ++i) {
    if (i%2==1) continue;
    double run1 = viov[i];
    double run2 = viov[i+1];
    double lum1 = hcumlum->GetBinContent(hcumlum->FindBin(run1));
    double lum2 = max(hcumlum->GetBinContent(hcumlum->FindBin(run2)),
		      hcumlum->GetBinContent(hcumlum->FindBin(run2)-1));
    //cout << "  run1 / lum1 : " << run1 << " / " << lum1
    //	 << ", run2 / lum2 : " << run2 << " / " << lum2 << endl << flush;
    f1->SetRange(lum1,lum2);

    hr2->Fit(f1,"QRNW");
    double p0 = f1->GetParameter(0);
    l->SetLineStyle(kSolid);
    l->SetLineColor(kRed);
    l->DrawLine(lum1,p0,lum2,p0);
  }

  gPad->RedrawAxis();
  c2->SaveAs("pdf/drawZcountingV2_ratio.pdf");

  
} // drawZcountingV2
