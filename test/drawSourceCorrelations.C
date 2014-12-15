#include "TH2D.h" 
#include "TLine.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "tdrstyle_mod14.C"

void drawSourceCorrelations() {

  setTDRStyle();

  // Autogenerate list of sources from UncertaintySources.txt with this macro:
  // grep '\[' CondFormats/JetMETObjects/data/Winter14_V5_DATA_UncertaintySources_AK5PFchs.txt | sed 's/\]/"\);/' | sed 's/\[/s.push_back\("/'
  // The ones that are included in Total (kData) are listed in JetDefs.hpp,
  // with bit-name mapping found in drawJetCorrectionUncertainty.cpp

  // 8 TeV Winter14_V5 sources
  vector<string> s;
  s.push_back("AbsoluteStat");
  s.push_back("AbsoluteScale");
  s.push_back("AbsoluteFlavMap");
  s.push_back("AbsoluteMPFBias");
  //s.push_back("HighPtExtra"); // 8 TeV GR
  s.push_back("Fragmentation"); // new patched 8 TeV
  s.push_back("SinglePionECAL");
  s.push_back("SinglePionHCAL");
  s.push_back("FlavorQCD");
  //s.push_back("TimeEta");
  //s.push_back("TimePt");
  s.push_back("RelativeJEREC1");
  s.push_back("RelativeJEREC2");
  s.push_back("RelativeJERHF");
  s.push_back("RelativePtBB");
  s.push_back("RelativePtEC1");
  s.push_back("RelativePtEC2");
  s.push_back("RelativePtHF");
  s.push_back("RelativeFSR");
  s.push_back("RelativeStatEC2");
  s.push_back("RelativeStatHF");
  s.push_back("PileUpDataMC");
  s.push_back("PileUpPtRef");
  s.push_back("PileUpPtBB");
  s.push_back("PileUpPtEC1");
  s.push_back("PileUpPtEC2");
  s.push_back("PileUpPtHF");
  //s.push_back("PileUpMuZero");
  //s.push_back("PileUpEnvelope");
  //s.push_back("SubTotalPileUp");
  //s.push_back("SubTotalRelative");
  //s.push_back("SubTotalPt");
  //s.push_back("SubTotalScale");
  //s.push_back("SubTotalMC");
  //s.push_back("Total");
  //s.push_back("TotalNoFlavor");
  //s.push_back("TotalNoTime");
  //s.push_back("TotalNoFlavorNoTime");
  //s.push_back("FlavorZJet");
  //s.push_back("FlavorPhotonJet");
  //s.push_back("FlavorPureGluon");
  //s.push_back("FlavorPureQuark");
  //s.push_back("FlavorPureCharm");
  //s.push_back("FlavorPureBottom");
  //s.push_back("TimeRunA");
  //s.push_back("TimeRunB");
  //s.push_back("TimeRunC");
  //s.push_back("TimeRunD");
  //s.push_back("CorrelationGroupMPFInSitu");
  //s.push_back("CorrelationGroupIntercalibration");
  //s.push_back("CorrelationGroupbJES");
  //s.push_back("CorrelationGroupFlavor");
  //s.push_back("CorrelationGroupUncorrelated");

  // 7 TeV JEC11_V12 sources
  /*
  vector<string> s;
  s.push_back("Absolute");
  s.push_back("HighPtExtra");
  s.push_back("SinglePion");
  s.push_back("Flavor");
  //s.push_back("Time");
  s.push_back("RelativeJEREC1");
  s.push_back("RelativeJEREC2");
  s.push_back("RelativeJERHF");
  s.push_back("RelativeFSR");
  s.push_back("RelativeStatEC2");
  s.push_back("RelativeStatHF");
  s.push_back("PileUpDataMC");
  s.push_back("PileUpOOT");
  s.push_back("PileUpPt");
  s.push_back("PileUpBias");
  s.push_back("PileUpJetRate");
  //s.push_back("SubTotalPileUp");
  //s.push_back("SubTotalRelative");
  //s.push_back("SubTotalPt");
  //s.push_back("SubTotalDataMC");
  //s.push_back("Total");
  */

  // define this bit to test that quadratic sum equals total
  string sref = "Total";

  double etabins[] = {0,1.3,2.5,3.0,4.7};
  const int neta = sizeof(etabins)/sizeof(etabins[0])-1;
  double ptbins[] =
    {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
     97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 362, 430,
     507, 592, 686, 790, 905, 1032, 1172, 1327, 1497, 1684, 1890, //1999};
     2000, 2238, 2500}; // pT=2504 observed limit with 20/fb at 8 TeV
  //, 2787, 3103, 3450};
  const int npt = sizeof(ptbins)/sizeof(ptbins[0])-1;

  const int netapt = neta*npt;
  double etaptbins[netapt+1];
  for (int ieta = 0; ieta != neta; ++ieta) {
    for (int ipt = 0; ipt != npt; ++ipt) {
      double x = pow(2500/10,ieta)*(ptbins[ipt]);
      etaptbins[ieta*npt+ipt] = x;
    } // for npt
  } // for ieta
  etaptbins[netapt] = pow(2500/10,neta)*10;

  TH2D *h2 = new TH2D("h2",";p_{T} (GeV);p_{T} (GeV)",
		      npt, ptbins, npt, ptbins);
  TH2D *h2x = new TH2D("h2x",";250^{i_{#eta}}p_{T} (GeV)"
		       ";250^{i_{#eta}}p_{T} (GeV)",
		       netapt, etaptbins, netapt, etaptbins);

  const double zmin = -0.4;//-0.25;
  // Set empty cells to white (below lower limit of zmin+eps)
  for (int i = 1; i != h2x->GetNbinsX()+1; ++i) {
    for (int j = 1; j != h2x->GetNbinsY()+1; ++j) {
      h2x->SetBinContent(i, j, zmin);
    }
  }

  // Create uncertainty sources
  vector<JetCorrectionUncertainty*> uncs(s.size());
  for (unsigned int k = 0; k != s.size(); ++k) {

    // GT 8 TeV
    //const char *sf = "CondFormats/JetMETObjects/data/"
    ////"Winter14_V5_DATA_UncertaintySources_AK7PF.txt"; // 8 TeV
    //"Winter14_V5_DATA_UncertaintySources_AK5PFchs.txt"; // 8 TeV
    // New patched 8 TeV (fixed fragmentation source sign)
    const char *sf = "../txt/"
    "Winter14_V8M_DATA_UncertaintySources_AK5PFchs.txt"; // 8 TeV
    // Older 7 TeV
    //const char *sf = "../../CondFormats/JetMETObjects/data/"
    //"JEC11_V12_AK7PF_UncertaintySources.txt"; // 7 TeV
    const char *src = s[k].c_str();
    JetCorrectorParameters *p = new JetCorrectorParameters(sf, src);
    JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
    uncs[k] = unc;
  }

  // Correlation plot between (eta,pT) pairs
  for (int ii = 0; ii != neta; ++ii) {
  for (int jj = 0; jj != neta; ++jj) {

    double eta1 = etabins[ii];
    double eta2 = etabins[jj];

  for (int i = 1; i != h2->GetNbinsX()+1; ++i) {
    for (int j = 1; j != h2->GetNbinsY()+1; ++j) {
      
      double pt1 = h2->GetXaxis()->GetBinCenter(i);
      double pt2 = h2->GetYaxis()->GetBinCenter(j);

      // Calculate correlation according to Sec. 8 Eq.(46) in JEC paper v3
      double sumski2(0), sumskj2(0), sumskiskj(0);
      for (unsigned int k = 0; k != s.size(); ++k) {
	
	JetCorrectionUncertainty *unc = uncs[k]; assert(unc);

	unc->setJetPt(pt1);
	unc->setJetEta(eta1);
	double ski = unc->getUncertainty(true);

	unc->setJetPt(pt2);
	unc->setJetEta(eta2);
	double skj = unc->getUncertainty(true);

	sumski2 += ski*ski;
	sumskj2 += skj*skj;
	sumskiskj += ski*skj;
      }
      double si = sqrt(sumski2);
      double sj = sqrt(sumskj2);
      double rhoij = sumskiskj / (si*sj);
      
      if (pt1*cosh(eta1)<4000 && pt2*cosh(eta2)<4000) {
	
	// |eta|<1.3
	if (ii==0 && jj==0) {
	  if (!(pt1<400 && pt2>900)) // white background for legend
	    h2->SetBinContent(i, j, rhoij);
	  else
	    h2->SetBinContent(i, j, zmin);
	}
	// all etas
	if (!(ii==jj && pt1<150 && pt2>600)) // white background for labels
	  h2x->SetBinContent(ii*npt+i, jj*npt+j, rhoij);
	else
	  h2x->SetBinContent(i, j, zmin);
      }
    } // for j
  } // for i

  } // for jj
  } // for ii

  // Sigle eta bin
  {
    TH1D *h1 = new TH1D("h1",";p_{T} (GeV);p_{T} (GeV)",npt,ptbins);
    h1->SetMinimum(ptbins[0]);
    h1->SetMaximum(ptbins[npt]);
    //h1->GetXaxis()->SetMoreLogLabels();
    h1->GetXaxis()->SetNoExponent();
    //h1->GetYaxis()->SetMoreLogLabels();
    h1->GetYaxis()->SetNoExponent();
    

    TCanvas *c1 = tdrCanvas("c1",h1,2,0,kSquare); // 8 TeV
    //TCanvas *c1 = tdrCanvas("c1",h1,1,0,kSquare); // 7 TeV
    gPad->SetLogx();
    gPad->SetLogy();
    
    gStyle->SetPalette(1);
    gStyle->SetOptStat(0);
    
    h2->Draw("COLZ SAME");
    //h2->GetZaxis()->SetRangeUser(0,1-1e-4);
    h2->GetZaxis()->SetRangeUser(zmin+1e-4,1-1e-4);

    gPad->SetRightMargin(0.12);
    gPad->SetLeftMargin(0.17);
    gPad->RedrawAxis();
    gPad->Update();

    TLine *l = new TLine();
    l->DrawLine(10,10,2500,2500);
    l->SetLineStyle(kDashed);
    l->DrawLine(1000,10,1000,1000);
    l->DrawLine(1000,1000,2500,1000);
    l->DrawLine(200,10,200,200);
    l->DrawLine(200,200,2500,200);
    l->DrawLine(40,10,40,40);
    l->DrawLine(40,40,2500,40);

    TLatex *tex = new TLatex();
    tex->SetNDC(); tex->SetTextSize(0.045);
    tex->DrawLatex(0.20,0.87,"Anti-k_{T} R=0.5, PF+CHS");
    //tex->DrawLatex(0.20,0.80,"TotalNoTime");
    tex->DrawLatex(0.20,0.80,"Total excl. time");

    tex->SetTextFont(42);
    tex->SetTextSize(0.05);
    tex->DrawLatex(0.10,0.09,"10");

    c1->SaveAs("pdf/drawSourceCorrelations_8TeV_Eta13.pdf");
    //c1->SaveAs("pdf/drawSourceCorrelations_7TeV.pdf");
  }

  // Multiple eta bins
  {
    TH1D *h1x = new TH1D("h1x",";250^{i#eta} #times p_{T} (GeV)"
			 ";250^{i#eta} #times p_{T} (GeV)",
			 netapt,etaptbins);
    h1x->SetMinimum(etaptbins[0]);
    h1x->SetMaximum(etaptbins[netapt]);
    //h1x->GetXaxis()->SetMoreLogLabels();
    //h1x->GetXaxis()->SetNoExponent();
    //h1x->GetYaxis()->SetMoreLogLabels();
    //h1x->GetYaxis()->SetNoExponent();

    TCanvas *c1x = tdrCanvas("c1x",h1x,2,0,kSquare); // 8 TeV cross-eta
    gPad->SetLogx();
    gPad->SetLogy();
    
    gStyle->SetPalette(1);
    gStyle->SetOptStat(0);
    
    h2x->Draw("COLZ SAME");
    //h2x->GetZaxis()->SetRangeUser(0,1-1e-4);
    h2x->GetZaxis()->SetRangeUser(zmin+1e-4,1-1e-4);

    gPad->SetBottomMargin(0.15);
    h1x->GetXaxis()->SetTitleOffset(1.2);
    
    gPad->SetRightMargin(0.12);
    gPad->SetLeftMargin(0.17);
    gPad->RedrawAxis();
    gPad->Update();
    
    TLine *l = new TLine();
    l->DrawLine(10,10,2500*pow(250,3),2500*pow(250,3));
    l->SetLineStyle(kDashed);
    l->DrawLine(1000,10,1000,1000);
    l->DrawLine(1000,1000,2500*pow(250,3),1000);
    l->DrawLine(200,10,200,200);
    l->DrawLine(200,200,2500*pow(250,3),200);
    l->DrawLine(40,10,40,40);
    l->DrawLine(40,40,2500*pow(250,3),40);
    
    TLatex *tex = new TLatex();
    tex->SetNDC(); tex->SetTextSize(0.045);
    tex->DrawLatex(0.20,0.88,"Anti-k_{T} R=0.5, PF+CHS");//, TotalNoTime");
    //tex->DrawLatex(0.20,0.80,"TotalNoTime");
    tex->DrawLatex(0.71,0.88,"HF");
    tex->DrawLatex(0.530,0.69,"EC2");
    tex->DrawLatex(0.355,0.50,"EC1");
    tex->DrawLatex(0.190,0.30,"BB");

    tex->SetTextFont(42);
    tex->SetTextSize(0.05);
    tex->DrawLatex(0.10,0.09,"10");
    
    c1x->SaveAs("pdf/drawSourceCorrelations_8TeV.pdf");
  }
}
