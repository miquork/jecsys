#define PFhadronLoop_cxx
#include "PFhadronLoop.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "TProfile.h"
#include "TH2D.h"

#include <iostream>

#include "../tdrstyle_mod15.C"

using namespace std;

void PFhadronLoop::Loop()
{
//   In a ROOT session, you can do:
//      root> .L PFhadronLoop.C
//      root> PFhadronLoop t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   setTDRStyle();
   TDirectory *curdir = gDirectory;

   // Jet veto maps
   //TFile *fjv = new TFile("../JECDatabase/jet_veto_maps/Summer19UL16_V0/hotjets-UL16.root","READ");
   TFile *fjv = new TFile("../JECDatabase/jet_veto_maps/Summer19UL18_V1/hotjets-UL18.root","READ");
   assert(fjv && !fjv->IsZombie());
   
   //TH2D *h2jv = (TH2D*)fjv->Get("h2hot_ul16_plus_hbm2_hbp12_qie11");
   TH2D *h2jv = (TH2D*)fjv->Get("h2hot_ul18_plus_hem1516_and_hbp2m1");
   assert(h2jv);
   //TH2D *h2mc = (TH2D*)fjv->Get("h2hot_mc");
   //assert(h2mc);
   //h2jv->Add(h2mc);
   

   // Open output file for profiles
   TFile *fout = new TFile("rootfiles/Hadrons.root","UPDATE");
   //fout->mkdir("MC16APV");
   //fout->cd("MC16APV");
   //fout->mkdir("MC18");
   //fout->cd("MC18");
   fout->mkdir("DT18");
   fout->cd("DT18");
   TDirectory *d = gDirectory;

   // PtTrk binning
   const double ax[] = {3,3.5,4,4.5,5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,
			11.5,12,13,14,15,16,18,20,22,25,30,35,40,50};
   const int nx = sizeof(ax)/sizeof(ax[0])-1;

   // Maps of tracks
   TH2D *h2b = (TH2D*)h2jv->Clone("h2b"); h2b->Reset(); // b=before
   TH2D *h2a = (TH2D*)h2jv->Clone("h2a"); h2a->Reset(); // a=after

   // Summed calorimeter energies
   // NB: Ef_* seems to be internally rounded to about full 1%
   // => to avoid stripes, use bin width 1% and edges as 0.5% (middle full 1%)
   TH2D *h2p =new TH2D("h2p",";p_{T} (GeV);"
		       "(E_{ECAL} + E_{HCAL}) / E_{track}",nx,ax,
		       250,-0.005,2.495);
   TProfile *p = new TProfile("p",";p_{T} (GeV);"
			      "(E_{ECAL} + E_{HCAL}) / E_{track}",nx,ax);
   TProfile *pr = new TProfile("pr",";p_{T} (GeV);"
			       "(E_{ECAL,raw} + E_{HCAL,raw}) / E_{track}",
			       nx,ax);
   TProfile *px = new TProfile("px",";p_{T} (GeV);"
			       "(E_{ECAL,raw} + E_{HCAL,corr}) / E_{track}",
			       nx,ax);
   
   // Individual calorimeter energies
   TProfile *pe = new TProfile("pe",";p_{T} (GeV);E_{ECAL,raw} / E_{track}",
			       nx,ax);
   TProfile *ph = new TProfile("ph",";p_{T} (GeV);E_{HCAL,raw} / E_{track}",
			       nx,ax);
   TProfile *pre = new TProfile("pre",";p_{T} (GeV);E_{ECAL,raw} / E_{track}",
				nx,ax);
   TProfile *prh = new TProfile("prh",";p_{T} (GeV);E_{HCAL,raw} / E_{track}",
				nx,ax);

   // PF candidate vs track
   TProfile *pt = new TProfile("pt",";p_{T} (GeV);"
			       "E_{cand} / E_{track}",nx,ax);

   // Special selections
   /*
   TProfile *ps = (TProfile*)p->Clone("ps");
   TProfile *psr = (TProfile*)pr->Clone("psr");
   TProfile *psx = (TProfile*)px->Clone("psx");
   TProfile *pst = (TProfile*)pt->Clone("pst");
   */

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (jentry%100000==0) cout << "." << flush;
      
      // Track distribution before cuts
      h2b->Fill(Eta, Phi);

      // Jet veto
      int iv = h2jv->GetXaxis()->FindBin(Eta);
      int jv = h2jv->GetYaxis()->FindBin(Phi);
      bool isveto = (h2jv->GetBinContent(iv,jv)>0);
      
      // Correction for cases where Pt is PtTrk scaled to ECAL+HCAL
      double c = (PtTrk>0 ? Pt/PtTrk : 1);
      // For PtTrk>20 GeV in MC18:
      // 1) Only removes 2/3 of the spike at 1 for ECAL+HCAL?
      // 2) But also increases peak at 0.1 by same amount; x3 in relative terms
      // For PtTrk>20 GeV in MC16APV:
      // 1) Removes 90% of the spike at 1 for ECAL+HCAL
      // 2) Increases peak at 0.1 by same amount; +25%(?) in relative terms
      // But ok, lot coming from MC16 hot region. Not such behavior in MC18
      // MC18 has 7k/2.3M tracks with PtTrk>20 GeV, while MC16APV has 45k/4.9M

      // Pt!=PtTrk when E_track>>E_calo, probably to avoid high pt fake tracks
      bool isGood = (PtTrk>0 && Pt==PtTrk);

      // Sum of calorimeter energies vs PtTrk
      double fc = c*(Ef_ECAL+Ef_HCAL);

      // Regular seletion with jet veto
      if (fabs(Eta)<1.3 && isGood && !isveto) {

	// Track distribution after cuts
	h2a->Fill(Eta, Phi);

	// Summed calorimeter energies
	h2p->Fill(PtTrk,  Ef_ECAL+Ef_HCAL);
	p->Fill(PtTrk,  Ef_ECAL+Ef_HCAL);
	pr->Fill(PtTrk, Ef_ECALRaw+Ef_HCALRaw);
	px->Fill(PtTrk, Ef_ECALRaw+Ef_HCAL);
	
	// Individual calorimeter energies
	pe->Fill(PtTrk, Ef_ECAL);
	ph->Fill(PtTrk, Ef_HCAL);
	pre->Fill(PtTrk, Ef_ECALRaw);
	prh->Fill(PtTrk, Ef_HCALRaw);

	// PF candidate vs track
	pt->Fill(PtTrk, Pt / PtTrk);
      } // good barrel track with jet veto

      // Special selection with good matching track and PF candidate
      //if (fabs(Eta)<1.3 && PtTrk>0 && !isveto && fabs(Pt/PtTrk-1)<0.05) {

      /*
      // Special selection using PF candidate Pt
      if (fabs(Eta)<1.3 && isGood && !isveto) { // && fabs(fc-1)<0.7) {

	// Summed calorimeter energies
	ps->Fill(Pt,  c*(Ef_ECAL+Ef_HCAL));
	psr->Fill(Pt, c*(Ef_ECALRaw+Ef_HCALRaw));
	psx->Fill(Pt, c*(Ef_ECALRaw+Ef_HCAL));

	// PF candidate vs track
	pst->Fill(Pt, Pt / PtTrk);
      } // good barrel track without jet veto
      */
   } // for jentry in nentries
   cout << endl << "Loop finished, writing results" << endl << flush;

   fout->Write(0,TObject::kOverwrite);
   fout->Close();
} // PFhadronLoop::Loop()
