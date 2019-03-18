// Purpose: Draw weight of each (eta,rho) bin for each y bin
//          (as used in inclusive jet analysis).
//          Store results in TH2Ds and TH3D in jerweights2018.root
// Run with 'root -l minitools/drawJERWeights.C'

#include "TFile.h"
#include "TH2D.h"
#include "TH3D.h"
#include <iostream>

#include "../tdrstyle_mod15.C"

using namespace std;

void drawJERweights() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TFile *f = new TFile("rootfiles/output-DATA-1-Fall18V8-D.root","READ");
  assert(f && !f->IsZombie());

  TFile *fout = new TFile("rootfiles/jecweights_2018.root","RECREATE");
  curdir->cd();

  //string path = "Standard/Eta_0.0-0.5/jt500";
  //string path = "Standard/Eta_0.5-1.0/jt500";
  //string path = "Standard/Eta_1.0-1.5/jt500";
  //string path = "Standard/Eta_1.5-2.0/jt500";
  //string path = "Standard/Eta_2.0-2.5/jt500";
  //string path = "Standard/Eta_2.5-3.0/jt500";
  //string path = "Standard/Eta_3.0-3.2/jt500";
  //string path = "Standard/Eta_3.2-4.7/jt500";

  const int ny = 8;
  double vy[ny+1] =
    {0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.2, 4.7};
  const int neta = 13 ;
  double veta[neta+1] =
    {0, 0.5, 0.8, 1.1, 1.3, 1.7, 1.9, 2.1, 2.3, 2.5, 2.8, 3, 3.2, 4.7};
  // must shuffle this around to match heta bin edges and avoid double-counting
  double veta2[neta+1] =
    {0, 0.522, 0.783, 1.131, 1.305, 1.653, 1.930, 2.172, 2.322, 2.500,
     2.853, 2.964, 3.314, 4.716};
  const int nrho = 7;
  double vrho[nrho+1] =
  //{0, 6.37, 12.4, 18.42, 24.45, 30.47, 36.49, 42.52};
  // must round this to match hrho bin edges and avoid double-counting
    {0, 6.5, 12.5, 18.5, 24.5, 30.5, 36.5, 42.5};

  TH3D *h3 = new TH3D("h3",";y_{jet};#eta;#rho (GeV);",
		      ny,vy, neta,veta, nrho,vrho);

  for (int iy = 0; iy != ny; ++iy) {

    double ymin(vy[iy]), ymax(vy[iy+1]); int trg(500);
    string path = Form("Standard/Eta_%1.1f-%1.1f/jt%d",ymin,ymax,trg);
    //float ymin(0), ymax(0); int trg(0);
    //assert(sscanf(cp,"Standard/Eta_%f-%f/jt%d",&ymin,&ymax,&trg)==3);
    cout << "ymin="<<ymin<<" ymax="<<ymax<<" trg="<<trg<<endl;
    cout << "path " << path << endl;

    const char *cp = path.c_str();
    TH1D *hrho = (TH1D*)f->Get(Form("%s/hrho",cp));
    assert(hrho);
    TH1D *heta = (TH1D*)f->Get(Form("%s/heta",cp));
    assert(heta);
    
    TH2D *h2 = new TH2D(Form("h2_%d",iy),";#eta;#rho (GeV);",
			neta,veta,nrho,vrho);

    double sumfij(0); double eps(1e-4);
    for (int ieta = 1; ieta != h2->GetNbinsX()+1; ++ieta) {
      for (int irho = 1; irho != h2->GetNbinsY()+1; ++irho) {
	
	// patch eta bin edges in heta
	double etamin = veta2[ieta-1];//h2->GetXaxis()->GetBinLowEdge(ieta);
	double etamax = veta2[ieta];//h2->GetXaxis()->GetBinLowEdge(ieta+1);
	
	// patch eta bin edges in heta
	//if (fabs(ymax-0.5)<eps) {
	//if (fabs(etamin-0.522)<eps) etamin=0.522;
	//if (fabs(etamax-0.522)<eps) etamax=0.522;
	//}
	if (fabs(ymax-1.0)<eps) {
	  if (fabs(etamin-0.522)<eps) etamin=0.435;
	  if (fabs(etamax-0.522)<eps) etamax=0.435;
	}
	if (fabs(ymax-3.0)<eps) {
	  if (fabs(etamin-2.964)<eps) etamin=3.139;
	  if (fabs(etamax-2.964)<eps) etamax=3.139;
	}
	if (fabs(ymax-4.7)<eps) {
	  if (fabs(etamin-3.314)<eps) etamin=3.139;
	  if (fabs(etamax-3.314)<eps) etamax=3.139;
	}

	double fieta = (heta->Integral(heta->FindBin(-etamax),
				       heta->FindBin(-etamin-eps)) +
			heta->Integral(heta->FindBin(etamin),
				       heta->FindBin(etamax-eps))) /
	  heta->Integral();
	
	double rhomin = h2->GetYaxis()->GetBinLowEdge(irho);
	double rhomax = h2->GetYaxis()->GetBinLowEdge(irho+1);
	
	double firho = (hrho->Integral(hrho->FindBin(rhomin),
				       hrho->FindBin(rhomax-eps))) / 
	  hrho->Integral();
	
	double fij = fieta * firho;
	sumfij += fij;
	h2->SetBinContent(ieta, irho, fij);
	h3->SetBinContent(iy, ieta, irho, fij);
      } // for irho
    } // for ieta
    
    TCanvas *c1 = new TCanvas(Form("c1_%d",iy),Form("c1_%d",iy),600,600);
    h2->Draw("COLZ");
    gPad->SetRightMargin(0.15);
    cout << " sum_ij (f_ij) = " << sumfij << endl;
    cout << " integral TH2D = " << h2->Integral() << endl;
    gPad->Update();

    // Pause for each bin (does not draw canvas?)
    //cout << "Press any key to continue." << endl;
    //cin.ignore();
    fout->cd();
    h2->Write();
    curdir->cd();
  } // for iy

  cout << " integral T32D = " << h3->Integral() << endl;
  gPad->Update();

  fout->cd();
  h3->Write();
  curdir->cd();

  fout->Write();
  fout->Close();
}
