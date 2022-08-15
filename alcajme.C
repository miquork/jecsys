#define alcajme_cxx
#include "alcajme.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TLorentzVector.h"
#include "TProfile.h"
#include "TProfile2D.h"


#include <iostream>

double DELTAPHI(double a, double b) {
  double phi1 = max(a,b);
  double phi2 = min(a,b);
  double d = phi1-phi2;
  if (d>TMath::Pi()) d -= TMath::TwoPi();
  return fabs(d);
}

void alcajme::Loop()
{
//   In a ROOT session, you can do:
//      root> .L alcajme.C
//      root> alcajme t
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

   bool cloneTree = false;
   TTree *fChain2(0);
   TFile *fout2(0);
   TDirectory *curdir = gDirectory;
   if (cloneTree) {
     cout << "Cloning TTree..." << flush;
     fout2 = new TFile("../data/AlCaRaw/skimmed_ntuple.root","RECREATE");
     fChain2 = fChain->CloneTree(1000000);
     fout2->Write();
     fout2->Close();
     curdir->cd();
     cout << "Ready." << endl << flush;
   }

   TLorentzVector p4, p4mht, p4mht2, p4mhtc, p4mhtc3, p4t, p4p, p4b, p4bx;
   TLorentzVector p4m0, p4m2, p4mn, p4mu;
   TFile *fout = new TFile("rootfiles/alcajme_out.root","RECREATE");

   // Regular L2Relative and L2Res eta binning
   double etabins[] =
     {-5.191,
      -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, -3.489,
      -3.314, -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172, -2.043,
      -1.93, -1.83, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218,
      -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, -0.435, 
      -0.348, -0.261, -0.174, -0.087, 0, 0.087, 0.174, 0.261, 0.348, 0.435,
      0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305,
      1.392, 1.479, 1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5,
      2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191,
      4.363, 4.538, 4.716, 4.889, 5.191};
   const int neta = sizeof(etabins)/sizeof(etabins[0])-1;
   const int ny = 600;
   double vy[ny+1];
   for (int i = 0; i != ny+1; ++i) vy[i] = -2. + (+4.+2.)/ny*i;

   // L2Res profiles for HDM method
   TProfile *pm0aeta = new TProfile("pm0aeta","PtAve;#eta;MPF",neta,etabins);
   TProfile *pm2aeta = new TProfile("pm2aeta","PtAve;#eta;MPF2",neta,etabins);
   TProfile *pmnaeta = new TProfile("pmnaeta","PtAve;#eta;MPFn",neta,etabins);
   TProfile *pmuaeta = new TProfile("pmuaeta","PtAve;#eta;MPFu",neta,etabins);

   TProfile *pm0teta = new TProfile("pm0teta","PtTag;#eta;MPF",neta,etabins);
   TProfile *pm2teta = new TProfile("pm2teta","PtTag;#eta;MPF2",neta,etabins);
   TProfile *pmnteta = new TProfile("pmnteta","PtTag;#eta;MPFn",neta,etabins);
   TProfile *pmuteta = new TProfile("pmuteta","PtTag;#eta;MPFu",neta,etabins);

   TProfile *pm0peta = new TProfile("pm0peta","PtProbe;#eta;MPF",neta,etabins);
   TProfile *pm2peta = new TProfile("pm2peta","PtProbe;#eta;MPF2",neta,etabins);
   TProfile *pmnpeta = new TProfile("pmnpeta","PtProbe;#eta;MPFn",neta,etabins);
   TProfile *pmupeta = new TProfile("pmupeta","PtProbe;#eta;MPFu",neta,etabins);

   // 2D distributions for MPF-MPFX method
   TH2D *h2m0aeta = new TH2D("h2m0aeta","PtAve;#eta;MPF",neta,etabins,ny,vy);
   TH2D *h2m0xaeta = new TH2D("h2m0xaeta","PtAve;#eta;MPFX",neta,etabins,ny,vy);
   TH2D *h2m0teta = new TH2D("h2m0teta","PtTag;#eta;MPF",neta,etabins,ny,vy);
   TH2D *h2m0xteta= new TH2D("h2m0xteta","PtTag;#eta;MPFX",neta,etabins,ny,vy);
   TH2D *h2m0peta = new TH2D("h2m0peta","PtProbe;#eta;MPF",neta,etabins,ny,vy);
   TH2D *h2m0xpeta=new TH2D("h2m0xpeta","PtProbe;#eta;MPFX",neta,etabins,ny,vy);
   //
   TH2D *h2m2aeta = new TH2D("h2m2aeta","PtAve;#eta;MPF2",neta,etabins,ny,vy);
   TH2D *h2m2xaeta= new TH2D("h2m2xaeta","PtAve;#eta;MPF2X",neta,etabins,ny,vy);
   TH2D *h2m2teta = new TH2D("h2m2teta","PtTag;#eta;MPF2",neta,etabins,ny,vy);
   TH2D *h2m2xteta= new TH2D("h2m2xteta","PtTag;#eta;MPF2X",neta,etabins,ny,vy);
   TH2D *h2m2peta = new TH2D("h2m2peta","PtProb;#eta;MPF2",neta,etabins,ny,vy);
   TH2D *h2m2xpeta=new TH2D("h2m2xpeta","PtProb;#eta;MPF2X",neta,etabins,ny,vy);


   // 2D eta-phi histograms and profiles for jet veto maps
   TH2D *h2etaphi = new TH2D("hetaphi","All;#eta;#phi",neta,etabins,
			     72,-TMath::Pi(),TMath::Pi());
   TH2D *h2etaphi40 = new TH2D("hetaphi40","All;#eta;#phi",neta,etabins,
			       72,-TMath::Pi(),TMath::Pi());
   TProfile2D *p2etaphia = new TProfile2D("p2etaphia",
					 "Tag-and-probe PtAve;#eta;#phi;MPF2",
					  neta,etabins,
					  72,-TMath::Pi(),TMath::Pi());
   TProfile2D *p2etaphit = new TProfile2D("p2etaphit",
					 "Tag-and-probe PtTag;#eta;#phi;MPF2",
					  neta,etabins,
					  72,-TMath::Pi(),TMath::Pi());
   TProfile2D *p2etaphip = new TProfile2D("p2etaphip",
					 "Tag-and-probe PtProbe;#eta;#phi;MPF2",
					  neta,etabins,
					  72,-TMath::Pi(),TMath::Pi());
   
   // Controls
   TH1D *hpt = new TH1D("hpt","hpt",500,0,500);
   TH1D *hpt30 = new TH1D("hpt30","hpt30",500,0,500);
   TH1D *heta = new TH1D("heta","heta",104,-5.2,5.2);
   TH1D *heta40 = new TH1D("heta40","heta40",104,-5.2,5.2);
   TH1D *hnjet = new TH1D("hnjet","hnjet",500,0,500);
   TH1D *hnjet30 = new TH1D("hnjet30","hnjet30",500,0,500);

   TH1D *hpta = new TH1D("hpta","hpta",500,0,500);
   TH1D *hptt = new TH1D("hptt","hpta",500,0,500);
   TH1D *hptp = new TH1D("hptp","hpta",500,0,500);

   TH1D *hav = new TH1D("hav","Asymmetry",400,-2,2);
   TH1D *hat = new TH1D("hat","Pt,probe / Pt,tag",400,0,4);
   TH1D *hap = new TH1D("hap","Pt,tag / Pt,probe",400,0,4);

   TH1D *havb = new TH1D("havb","Asymmetry",400,-2,2);
   TH1D *hatb = new TH1D("hatb","Pt,probe / Pt,tag",400,0,4);
   TH1D *hapb = new TH1D("hapb","Pt,tag / Pt,probe",400,0,4);
   
   TH1D *hmv = new TH1D("hmv","MPF PtAve",800,-2,6);
   TH1D *hmt = new TH1D("hmt","MPF PtTag",800,-2,6);
   TH1D *hmp = new TH1D("hmp","MPF PtProbe",800,-2,6);

   TH1D *hmva = new TH1D("hmva","MPFA PtAve",800,-2,6);
   TH1D *hmvc = new TH1D("hmvc","MPFC PtAve",800,-2,6);
   TH1D *hmvc3 = new TH1D("hmvc3","MPFC3 PtAve",800,-2,6);
   TH1D *hmv2 = new TH1D("hmv2","MPF2 PtAve",800,-2,6);
   //TH1D *hmv2 = new TH1D("hmv2","MPF2 PtAve",400,-2,2);


   //Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nentries = fChain->GetEntries();
   cout << "Loaded " << nentries << " entries" << endl;
   //nentries = min((int)nentries, (int)1000000);
   cout << "Processing " << nentries << " entries" << endl;
   //TTimer t;
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      //if (jentry%100==0) cout << "." << flush;
      //if (jentry%10000==0) cout << "." << flush;
      if (jentry%100000==0) cout << "." << flush;
      if (jentry%5000000==0) cout << "n="<<jentry<<endl<<flush;
      //if (jentry%1000000==0) cout << timer

      vector<float> *vpt = hltPFJetsCorrectedMatchedToCaloJets10_pt;
      vector<float> *veta = hltPFJetsCorrectedMatchedToCaloJets10_eta;
      vector<float> *vphi = hltPFJetsCorrectedMatchedToCaloJets10_phi;
      vector<float> *vmass = hltPFJetsCorrectedMatchedToCaloJets10_mass;

      int njet = vpt->size();
      int njet3 = 0;
      int njetn = 0;

      p4mht.SetPtEtaPhiM(0,0,0,0);
      p4mht2.SetPtEtaPhiM(0,0,0,0);
      p4mhtc.SetPtEtaPhiM(0,0,0,0);
      p4mhtc3.SetPtEtaPhiM(0,0,0,0);
      p4m0.SetPtEtaPhiM(0,0,0,0);
      p4m2.SetPtEtaPhiM(0,0,0,0);
      p4mn.SetPtEtaPhiM(0,0,0,0);
      p4mu.SetPtEtaPhiM(0,0,0,0);
      for (int i = 0; i != njet; ++i) {
	
	p4.SetPtEtaPhiM((*vpt)[i], (*veta)[i], (*vphi)[i], (*vmass)[i]);
	
	//hpt->Fill((*vpt)[i]);
	hpt->Fill(p4.Pt());
	if (fabs(p4.Eta())<3.0)
	  hpt30->Fill(p4.Pt());
	heta->Fill(p4.Eta());
	h2etaphi->Fill(p4.Eta(),p4.Phi());
	if (p4.Pt()>40.) {
	  heta40->Fill(p4.Eta());
	  h2etaphi40->Fill(p4.Eta(),p4.Phi());
	}
	p4mht -= p4;
	if (i<2) p4mht2 -= p4;
	if (fabs(p4.Eta())<3.0 || i<2) {
	  p4mhtc -= p4;
	  if (njet3<3) {
	    p4mhtc3 -= p4;
	  }
	  ++njet3;
	}
	
	// L2Res HDM
	p4m0 -= p4;
	if (i<2) { // leading jet
	  p4m2 -= p4;
	}
	else if (fabs(p4.Eta())<3.0 && njetn<3) { // soft jets
	  p4mn -= p4;
	  ++njetn;
	}
	else { // other "unclustered"
	  p4mu -= p4;
	}
      } //  for i in njet
      hnjet->Fill(njet);
      hnjet30->Fill(njet3);

      // dijet selection
      if (njet>=2) {

	// both leading jets act as tag and probe in turn
	for (int itag = 0; itag != 2; ++itag) {

	  int iprobe = (itag == 0 ? 1 : 0);
	  p4t.SetPtEtaPhiM((*vpt)[itag], (*veta)[itag], (*vphi)[itag],
			   (*vmass)[itag]);
	  p4p.SetPtEtaPhiM((*vpt)[iprobe], (*veta)[iprobe], (*vphi)[iprobe],
			   (*vmass)[iprobe]);
	  // bisector axis
	  p4b.SetPtEtaPhiM(0,0,0,0);
	  p4b += p4t;
	  p4b -= p4p;
	  p4b.SetPtEtaPhiM(p4b.Pt(),0.,p4b.Phi(),0.);
	  //p4b *= 1./fabs(p4b.Mag());
	  p4b *= 1./p4b.Pt();
	  p4bx.SetPtEtaPhiM(p4b.Pt(),0.,p4b.Phi()+0.5*TMath::Pi(),0.);

	  double ptave = 0.5*(p4t.Pt()+p4p.Pt());
	  double asymm = (p4p.Pt() - p4t.Pt()) / ptave;
	  //double deltaphi = fabs(p4t.DeltaPhi(p4p));
	  double deltaphi = DELTAPHI(p4t.Phi(),p4p.Phi());
	  //double mpf = 1 + p4p.Dot(p4mht)/ptave;
	  double mpf = 1 + (p4mht.Vect().Dot(p4b.Vect()))/ptave;
	  double mpf2 = 1 + (p4mht2.Vect().Dot(p4b.Vect()))/ptave;
	  double mpfc = 1 + (p4mhtc.Vect().Dot(p4b.Vect()))/ptave;
	  double mpfc3 = 1 + (p4mhtc3.Vect().Dot(p4b.Vect()))/ptave;

	  double m0 = 1 + (p4m0.Vect().Dot(p4b.Vect()))/ptave;
	  double m2 = 1 + (p4m2.Vect().Dot(p4b.Vect()))/ptave;
	  double mn = 0 + (p4mn.Vect().Dot(p4b.Vect()))/ptave;
	  double mu = 0 + (p4mu.Vect().Dot(p4b.Vect()))/ptave;

	  double m0x = 1 + (p4m0.Vect().Dot(p4bx.Vect()))/ptave;
	  double m2x = 1 + (p4m2.Vect().Dot(p4bx.Vect()))/ptave;

	  if (fabs(p4t.Eta())<1.3 && deltaphi>2.7) {
	    double eta = p4p.Eta();
	    if (ptave >= 40 && ptave <50) {
	      pm0aeta->Fill(eta, m0);
	      pm2aeta->Fill(eta, m2);
	      pmnaeta->Fill(eta, mn);
	      pmuaeta->Fill(eta, mu);

	      h2m0aeta->Fill(eta, m0);
	      h2m0xaeta->Fill(eta, m0x);
	      h2m2aeta->Fill(eta, m2);
	      h2m2xaeta->Fill(eta, m2x);

	      p2etaphia->Fill(eta, p4p.Phi(), asymm);
	    }
	    double pttag = p4t.Pt();
	    if (pttag >= 40 && pttag <50) {
	      pm0teta->Fill(eta, m0);
	      pm2teta->Fill(eta, m2);
	      pmnteta->Fill(eta, mn);
	      pmuteta->Fill(eta, mu);

	      h2m0teta->Fill(eta, m0);
	      h2m0xteta->Fill(eta, m0x);
	      h2m2teta->Fill(eta, m2);
	      h2m2xteta->Fill(eta, m2x);

	      p2etaphit->Fill(eta, p4p.Phi(), p4p.Pt()/pttag);
	    }
	    double ptprobe = p4p.Pt();
	    if (ptprobe >= 40 && ptprobe <50) {
	      pm0peta->Fill(eta, m0);
	      pm2peta->Fill(eta, m2);
	      pmnpeta->Fill(eta, mn);
	      pmupeta->Fill(eta, mu);

	      h2m0peta->Fill(eta, m0);
	      h2m0xpeta->Fill(eta, m0x);
	      h2m2peta->Fill(eta, m2);
	      h2m2xpeta->Fill(eta, m2x);

	      p2etaphip->Fill(eta, p4p.Phi(), p4t.Pt()/ptprobe);
	    }
	  }

	  if (fabs(p4t.Eta())<1.3 && deltaphi>2.7) {

	    hpta->Fill(ptave);
	    hptt->Fill(p4t.Pt());
	    hptp->Fill(p4p.Pt());
	    
	    if (ptave>=40) {
	      hav->Fill(asymm);
	      if (fabs(p4p.Eta())<1.3) {
		havb->Fill(asymm);
	      }

	      /*
	      cout << "itag="<<itag<<" hav="<<asymm<<" hmv2="<<mpf2
		   << " ptt="<<p4t.Pt()<<" ptp="<<p4p.Pt()<<endl
		   << " phit="<<p4t.Phi()<<" phip="<<p4p.Phi()
		   << " ptave="<<ptave<<endl
		   << " phib="<<p4b.Phi()
		   << " dphi="<<deltaphi 
		   << " magb="<<p4b.Mag()<<endl
		   << " ptmht2="<<p4mht2.Pt()
		   << " phimht2="<<p4mht2.Phi()
		   << endl;
	      */
	      hmv->Fill(mpf);
	      hmva->Fill(mpf);
	      hmvc->Fill(mpfc);
	      hmvc3->Fill(mpfc3);
	      hmv2->Fill(mpf2);
	      //cout << p4mht.Pt() << " / " << ptave << " vs a=" << mpf << endl;
	      //cout << p4mhtc.Pt() << " / " << ptave << " vs c=" << mpfc << endl;
	      //cout << p4mhtc3.Pt() << " / "<< ptave << " vs c3="<< mpfc3<< endl;
	      //if (fabs(p4p.Eta())<1.3)
	      //hmvb->Fill(mpf);
	    }
	    if (p4t.Pt()>=40) {
	      hat->Fill(p4p.Pt() / p4t.Pt());
	      if (fabs(p4p.Eta())<1.3) {
		hatb->Fill(p4p.Pt() / p4t.Pt());
	      }
	      hmt->Fill(mpf);
	    }
	    if (p4p.Pt()>=40) {
	      hap->Fill(p4t.Pt() / p4p.Pt());
	      if (fabs(p4p.Eta())<1.3) {
		hapb->Fill(p4t.Pt() / p4p.Pt());
	      }
	      hmp->Fill(mpf);
	    }

	  } // tag and deltaphi
	} // for itag
      } // dijet


      //if (cloneTree) fChain2->Fill();
   } // for jentry

   //TFile *fout = new TFile("rootfiles/alcajme_out.root","RECREATE");
   //hpt->Write("hpt",TObject::kOverwrite);
   fout->Write();
   fout->Close();

   //if (cloneTree) {
   //fChain2->Write();
   //fout2->Write();
   //fout2->Close();
   //}

   /*
   //hpt->Draw();
   //hpta->Draw();
   //hptt->Draw();
   //hptp->Draw();
   //heta->Draw();
   hav->Rebin(10);
   //hav->Draw();
   hmv->Rebin(10);
   hmv->Draw();
   hat->Rebin(10);
   //hat->Draw();
   hap->Rebin(10);
   //hap->Draw();

   // Barrel
   havb->Rebin(10);
   //havb->Draw();
   hatb->Rebin(10);
   //hatb->Draw();
   hapb->Rebin(10);
   //hapb->Draw();
   */
}
