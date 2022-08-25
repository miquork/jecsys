#define alcajme_cxx
#include "alcajme.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TLorentzVector.h"
#include "TProfile.h"
#include "TProfile2D.h"


#include <iostream>
#include <fstream>
#include <set>
#include <map>

double DELTAPHI(double a, double b) {
  double phi1 = max(a,b);
  double phi2 = min(a,b);
  double d = phi1-phi2;
  if (d>TMath::Pi()) d -= TMath::TwoPi();
  return fabs(d);
}

class basicHistos {
public:
  double ptmin;
  double ptmax;

  TProfile *pm0ab;
  TProfile *pm2ab;
  TProfile *pmnab;
  TProfile *pmuab;
  TProfile *pmoab;

  TProfile *pm0tb;
  TProfile *pm2tb;
  TProfile *pmntb;
  TProfile *pmutb;
  TProfile *pmotb;

  TProfile *pm0pb;
  TProfile *pm2pb;
  TProfile *pmnpb;
  TProfile *pmupb;
  TProfile *pmopb;
 
  TProfile *pm0ac;
  TProfile *pm2ac;
  TProfile *pmnac;
  TProfile *pmuac;
  TProfile *pmoac;

  TProfile *pm0tc;
  TProfile *pm2tc;
  TProfile *pmntc;
  TProfile *pmutc;
  TProfile *pmotc;

  TProfile *pm0af;
  TProfile *pm2af;
  TProfile *pmnaf;
  TProfile *pmuaf;
  TProfile *pmoaf;

  TProfile *pm0pf;
  TProfile *pm2pf;
  TProfile *pmnpf;
  TProfile *pmupf;
  TProfile *pmopf;

  // 2D distributions for MPF-MPFX method
  TH2D *h2m0ab;
  TH2D *h2m0abx;
  TH2D *h2m0tb;
  TH2D *h2m0tbx;
  TH2D *h2m0pb;
  TH2D *h2m0pbx;
   //
  TH2D *h2m2ab;
  TH2D *h2m2abx;
  TH2D *h2m2tb;
  TH2D *h2m2tbx;
  TH2D *h2m2pb;
  TH2D *h2m2pbx;

  TH2D *h2m0tc;
  TH2D *h2m0tcx;
  TH2D *h2m0pf;
  TH2D *h2m0pfx;
   //
  TH2D *h2m2tc;
  TH2D *h2m2tcx;
  TH2D *h2m2pf;
  TH2D *h2m2pfx;

  TProfile2D *p2etaphia, *p2etaphit, *p2etaphip;  
};

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

   fChain->SetBranchStatus("*",0);  // disable all branches

   fChain->SetBranchStatus("run",1);
   fChain->SetBranchStatus("luminosityBlock",1);
   fChain->SetBranchStatus("event",1);

   fChain->SetBranchStatus("hltPFJetsCorrectedMatchedToCaloJets10_pt",1);
   fChain->SetBranchStatus("hltPFJetsCorrectedMatchedToCaloJets10_eta",1);
   fChain->SetBranchStatus("hltPFJetsCorrectedMatchedToCaloJets10_phi",1);
   fChain->SetBranchStatus("hltPFJetsCorrectedMatchedToCaloJets10_mass",1);

   fChain->SetBranchStatus("hltPFJetsCorrectedMatchedToCaloJets10AK8_pt",1);

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

   bool doJSON = true;
   if (doJSON) {
     if (!LoadJSON()) {
       cout << "Issues loading the JSON file; aborting..." << endl << flush;
       exit(0);
     }
   }
   int _nbadevts_json(0);
   bool doAK4 = true;
   int _nbadevts_ak4(0);

   TLorentzVector p4dj;
   TLorentzVector p4, p4s, p4mht, p4mht2, p4mhtc, p4mhtc3, p4t, p4p;
   TLorentzVector p4b, p4bx, p4c, p4cx, p4f, p4fx;
   TLorentzVector p4m0, p4m2, p4mn, p4mu, p4mo;
   TFile *fout = new TFile("rootfiles/alcajme_out.root","RECREATE");

   // Regular L2Relative and L2Res eta binning
   double vx[] =
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
   const int nx = sizeof(vx)/sizeof(vx[0])-1;
   const int ny = 800;
   double vy[ny+1];
   for (int i = 0; i != ny+1; ++i) vy[i] = -3. + (+5.+3.)/ny*i;
   const int nz = 400;
   double vz[nz+1];
   for (int i = 0; i != nz+1; ++i) vz[i] = -1. + (+3.+1.)/nz*i;

   // L2Res profiles for HDM method
   // coding: m0=MPF, m2=DB, mn=n-jet, mu=uncl. (observable)
   //         a=PtAve, t=PtTag, p=PtProbe (binning)
   //         b=bisector, c=central(tag), f=forward(probe) (projection axis)
   TProfile *pm0ab = new TProfile("pm0ab","PtAve Bisector;#eta;MPF",nx,vx);
   TProfile *pm2ab = new TProfile("pm2ab","PtAve Bisector;#eta;MPF2",nx,vx);
   TProfile *pmnab = new TProfile("pmnab","PtAve Bisector;#eta;MPFn",nx,vx);
   TProfile *pmuab = new TProfile("pmuab","PtAve Bisector;#eta;MPFu",nx,vx);
   TProfile *pmoab = new TProfile("pmoab","PtAve Bisector;#eta;MPFo",nx,vx);

   TProfile *pm0tb = new TProfile("pm0tb","PtTag Bisector;#eta;MPF",nx,vx);
   TProfile *pm2tb = new TProfile("pm2tb","PtTag Bisector;#eta;MPF2",nx,vx);
   TProfile *pmntb = new TProfile("pmntb","PtTag Bisector;#eta;MPFn",nx,vx);
   TProfile *pmutb = new TProfile("pmutb","PtTag Bisector;#eta;MPFu",nx,vx);
   TProfile *pmotb = new TProfile("pmotb","PtTag Bisector;#eta;MPFo",nx,vx);

   TProfile *pm0pb = new TProfile("pm0pb","PtProbe Bisector;#eta;MPF",nx,vx);
   TProfile *pm2pb = new TProfile("pm2pb","PtProbe Bisector;#eta;MPF2",nx,vx);
   TProfile *pmnpb = new TProfile("pmnpb","PtProbe Bisector;#eta;MPFn",nx,vx);
   TProfile *pmupb = new TProfile("pmupb","PtProbe Bisector;#eta;MPFu",nx,vx);
   TProfile *pmopb = new TProfile("pmopb","PtProbe Bisector;#eta;MPFo",nx,vx);
 
   TProfile *pm0ac = new TProfile("pm0ac","PtAve Central;#eta;MPF",nx,vx);
   TProfile *pm2ac = new TProfile("pm2ac","PtAve Central;#eta;MPF2",nx,vx);
   TProfile *pmnac = new TProfile("pmnac","PtAve Central;#eta;MPFn",nx,vx);
   TProfile *pmuac = new TProfile("pmuac","PtAve Central;#eta;MPFu",nx,vx);
   TProfile *pmoac = new TProfile("pmoac","PtAve Central;#eta;MPFo",nx,vx);

   TProfile *pm0tc = new TProfile("pm0tc","PtTag Central;#eta;MPF",nx,vx);
   TProfile *pm2tc = new TProfile("pm2tc","PtTag Central;#eta;MPF2",nx,vx);
   TProfile *pmntc = new TProfile("pmntc","PtTag Central;#eta;MPFn",nx,vx);
   TProfile *pmutc = new TProfile("pmutc","PtTag Central;#eta;MPFu",nx,vx);
   TProfile *pmotc = new TProfile("pmotc","PtTag Central;#eta;MPFo",nx,vx);

   TProfile *pm0af = new TProfile("pm0af","PtAve Forward;#eta;MPF",nx,vx);
   TProfile *pm2af = new TProfile("pm2af","PtAve Forward;#eta;MPF2",nx,vx);
   TProfile *pmnaf = new TProfile("pmnaf","PtAve Forward;#eta;MPFn",nx,vx);
   TProfile *pmuaf = new TProfile("pmuaf","PtAve Forward;#eta;MPFu",nx,vx);
   TProfile *pmoaf = new TProfile("pmoaf","PtAve Forward;#eta;MPFo",nx,vx);

   TProfile *pm0pf = new TProfile("pm0pf","PtProbe Forward;#eta;MPF",nx,vx);
   TProfile *pm2pf = new TProfile("pm2pf","PtProbe Forward;#eta;MPF2",nx,vx);
   TProfile *pmnpf = new TProfile("pmnpf","PtProbe Forward;#eta;MPFn",nx,vx);
   TProfile *pmupf = new TProfile("pmupf","PtProbe Forward;#eta;MPFu",nx,vx);
   TProfile *pmopf = new TProfile("pmopf","PtProbe Forward;#eta;MPFo",nx,vx);

   // 2D distributions for MPF-MPFX method
   TH2D *h2m0ab = new TH2D("h2m0ab","PtAve Bisector;#eta;MPF",nx,vx,ny,vy);
   TH2D *h2m0abx = new TH2D("h2m0abx","PtAve Bisector;#eta;MPFX",nx,vx,ny,vy);
   TH2D *h2m0tb = new TH2D("h2m0tb","PtTag Bisector;#eta;MPF",nx,vx,ny,vy);
   TH2D *h2m0tbx = new TH2D("h2m0tbx","PtTag Bisector;#eta;MPFX",nx,vx,ny,vy);
   TH2D *h2m0pb  = new TH2D("h2m0pb","PtProbe Bisector;#eta;MPF",nx,vx,ny,vy);
   TH2D *h2m0pbx = new TH2D("h2m0pbx","PtProbe Bisector;#eta;MPFX",nx,vx,ny,vy);
   //
   TH2D *h2m2ab = new TH2D("h2m2ab","PtAve Bisector;#eta;MPF2",nx,vx,nz,vz);
   TH2D *h2m2abx = new TH2D("h2m2abx","PtAve Bisector;#eta;MPF2X",nx,vx,nz,vz);
   TH2D *h2m2tb = new TH2D("h2m2tb","PtTag Bisector;#eta;MPF2",nx,vx,nz,vz);
   TH2D *h2m2tbx = new TH2D("h2m2tbx","PtTag Bisector;#eta;MPF2X",nx,vx,nz,vz);
   TH2D *h2m2pb = new TH2D("h2m2pb","PtProb Bisector;#eta;MPF2",nx,vx,nz,vz);
   TH2D *h2m2pbx =new TH2D("h2m2pbx","PtProb Bisector;#eta;MPF2X",nx,vx,nz,vz);

   TH2D *h2m0tc = new TH2D("h2m0tc","PtTag Central;#eta;MPF",nx,vx,ny,vy);
   TH2D *h2m0tcx = new TH2D("h2m0tcx","PtTag Central;#eta;MPFX",nx,vx,ny,vy);
   TH2D *h2m0pf  = new TH2D("h2m0pf","PtProbe Forward;#eta;MPF",nx,vx,ny,vy);
   TH2D *h2m0pfx = new TH2D("h2m0pfx","PtProbe Forward;#eta;MPFX",nx,vx,ny,vy);
   //
   TH2D *h2m2tc = new TH2D("h2m2tc","PtTag Central;#eta;MPF2",nx,vx,nz,vz);
   TH2D *h2m2tcx = new TH2D("h2m2tcx","PtTag Central;#eta;MPF2X",nx,vx,nz,vz);
   TH2D *h2m2pf = new TH2D("h2m2pf","PtProb Forward;#eta;MPF2",nx,vx,nz,vz);
   TH2D *h2m2pfx =new TH2D("h2m2pfx","PtProb Forward;#eta;MPF2X",nx,vx,nz,vz);


   // 2D eta-phi histograms and profiles for jet veto maps
   TH2D *h2etaphi = new TH2D("h2etaphi","All;#eta;#phi",nx,vx,
			     72,-TMath::Pi(),TMath::Pi());
   TH2D *h2etaphi40 = new TH2D("h2etaphi40","All;#eta;#phi",nx,vx,
			       72,-TMath::Pi(),TMath::Pi());
   TProfile2D *p2etaphia = new TProfile2D("p2etaphia",
					 "Tag-and-probe PtAve;#eta;#phi;MPF2",
					  nx,vx,
					  72,-TMath::Pi(),TMath::Pi());
   TProfile2D *p2etaphit = new TProfile2D("p2etaphit",
					 "Tag-and-probe PtTag;#eta;#phi;MPF2",
					  nx,vx,
					  72,-TMath::Pi(),TMath::Pi());
   TProfile2D *p2etaphip = new TProfile2D("p2etaphip",
					 "Tag-and-probe PtProbe;#eta;#phi;MPF2",
					  nx,vx,
					  72,-TMath::Pi(),TMath::Pi());
   
   // Controls
   TH1D *hpt = new TH1D("hpt","hpt",500,0,500);
   TH1D *hpt30 = new TH1D("hpt30","hpt30",500,0,500);
   TH1D *heta = new TH1D("heta","heta",104,-5.2,5.2);
   TH1D *heta40 = new TH1D("heta40","heta40",104,-5.2,5.2);
   TH1D *hnjet = new TH1D("hnjet","hnjet",500,0,500);
   TH1D *hnjet30 = new TH1D("hnjet30","hnjet30",500,0,500);
   TH1D *hndup = new TH1D("hndup","hndup",100,0,100);

   TH1D *hnjetak8 = new TH1D("hnjetak8","hnjetak8",500,0,500);
   TH1D *hpt0withak8 = new TH1D("hpt0withak8","hpt0withak8",500,0,500);
   TH1D *hpt0noak8 = new TH1D("hpt0noak8","hptno0ak8",500,0,500);

   TH1D *hpta = new TH1D("hpta","PtAve",500,0,500);
   TH1D *hptt = new TH1D("hptt","PtTag",500,0,500);
   TH1D *hptp = new TH1D("hptp","PtProbe",500,0,500);

   const double twopi = TMath::TwoPi();
   TH1D *hdphi = new TH1D("hdphi","#Delta#phi",126,-twopi,+twopi);
   TH1D *hdphib = new TH1D("hdphib","#Delta#phi Barrel",126,-twopi,+twopi);
   TH2D *h2dphi = new TH2D("h2dphi","#Delta#phi vs #eta;#eta;#Delta#phi",
			   nx,vx,126,-TMath::TwoPi(),+TMath::TwoPi());
   TH1D *hdr = new TH1D("hdr","#DeltaR",220,0,11);//126,0,twopi);
   TH1D *hdrb = new TH1D("hdrb","#DeltaR Barrel",220,0,11);//126,0,twopi);
   TH1D *hdri = new TH1D("hdri","#DeltaR Inclusive",220,0,11);//126,0,twopi);
   TH1D *hdrib = new TH1D("hdrib","#DeltaR Inclusive Barrel",220,0,11);
   TH1D *hdribb = new TH1D("hdribb","#DeltaR Inclusive Barrel^{2}",220,0,11);
   TH1D *hdri2 = new TH1D("hdri2","#DeltaR Leads",220,0,11);//126,0,twopi);
   TH1D *hdri2b = new TH1D("hdri2b","#DeltaR Leads Barrel",220,0,11);
   TH1D *hdri2bb = new TH1D("hdri2bb","#DeltaR Leads Barrel^{2}",220,0,11);
   
   /*
   TH2D *h2etaphidr0 = new TH2D("h2etaphidr0","#DeltaR=0;#eta;#phi",nx,vx,
				72,-TMath::Pi(),TMath::Pi());
   TH2D *h2ijdr0 = new TH2D("hijdr0","#DeltaR=0;i;j",50,0,50,50,0,50);
   TH1D *havdr0 = new TH1D("havddr0","Asymmetry #DeltaR=0",400,-2,2);
   TH1D *hpt1dr0 = new TH1D("hpt1dr0","Pt1 #DeltaR=0",500,0,500);
   TH1D *hpt2dr0 = new TH1D("hpt2dr0","Pt2 #DeltaR=0",500,0,500);
   */

   TH1D *hav = new TH1D("hav","Asymmetry",400,-2,2);
   TH1D *hat = new TH1D("hat","Pt,probe / Pt,tag",400,0,4);
   TH1D *hap = new TH1D("hap","Pt,tag / Pt,probe",400,0,4);

   TH1D *havb = new TH1D("havb","Asymmetry Barrel",400,-2,2);
   TH1D *hatb = new TH1D("hatb","Pt,probe / Pt,tag Barrel",400,0,4);
   TH1D *hapb = new TH1D("hapb","Pt,tag / Pt,probe Barrel",400,0,4);
   
   TH1D *hmv = new TH1D("hmv","MPF PtAve",800,-2,6);
   TH1D *hmt = new TH1D("hmt","MPF PtTag",800,-2,6);
   TH1D *hmp = new TH1D("hmp","MPF PtProbe",800,-2,6);

   TH1D *hmva = new TH1D("hmva","MPFA PtAve",800,-2,6);
   TH1D *hmvc = new TH1D("hmvc","MPFC PtAve",800,-2,6);
   TH1D *hmvc3 = new TH1D("hmvc3","MPFC3 PtAve",800,-2,6);
   TH1D *hmv2 = new TH1D("hmv2","MPF2 PtAve",800,-2,6);
   //TH1D *hmv2 = new TH1D("hmv2","MPF2 PtAve",400,-2,2);

   // Extra fun
   // Run vs BX => need BX
   TH1D *hmjji = new TH1D("hmjji","Dijet mass, inclusive",1000,0,1000);
   TH1D *hmjji13 = new TH1D("hmjji13","Dijet mass, |#DeltaEta|<1.3, incl.",
			    1000,0,1000);
   TH1D *hmjj1 = new TH1D("hmjj1","Dijet mass, 1-lead",1000,0,1000);
   TH1D *hmjj113 = new TH1D("hmjj113","Dijet mass, |#DeltaEta|<1.3, 1-lead",
			    1000,0,1000);
   TH1D *hmjj2 = new TH1D("hmjj2","Dijet mass, 2-lead",1000,0,1000);
   TH1D *hmjj213 = new TH1D("hmjj213","Dijet mass, |#DeltaEta|<1.3, 2-lead",
			    1000,0,1000);
   TH1D *hmjj = new TH1D("hmjj","Dijet mass, tag-and-probe",1000,0,1000);
   TH1D *hmjj13 = new TH1D("hmjj13","Dijet mass, |#DeltaEta|<1.3, TnP",
			   1000,0,1000);

   // Listing of runs and LS
   int nrun(0), nls(0), nevt(0);
   map<int, map<int, int> > mrunls;

   // Copy L2Res histograms for multiple pT bins
   const double vpt[] = {40, 50, 60, 70, 85, 100, 125, 155, 180, 210, 250, 300};
   const int npt = sizeof(vpt)/sizeof(vpt[0])-1;
   TH1D *hptbins = new TH1D("hptbins",";p_{T} (GeV);N_{events}",npt,vpt);
   vector<basicHistos*> vh(npt);
   const bool doPtBins = true;
   if (doPtBins) {
     
     for (int ipt = 0; ipt != npt; ++ipt) {

       fout->mkdir(Form("Pt_%d_%d",int(vpt[ipt]),int(vpt[ipt+1])));
       fout->cd(Form("Pt_%d_%d",int(vpt[ipt]),int(vpt[ipt+1])));
       basicHistos *h = new basicHistos();
       vh[ipt] = h;
       h->ptmin = vpt[ipt];
       h->ptmax = vpt[ipt+1];

       //if (ptave >= 40 && ptave <50) {
       h->pm0ab = (TProfile*)pm0ab->Clone();
       h->pm2ab = (TProfile*)pm2ab->Clone();
       h->pmnab = (TProfile*)pmnab->Clone();
       h->pmuab = (TProfile*)pmuab->Clone();
       h->pmoab = (TProfile*)pmoab->Clone();
     
       h->pm0ac = (TProfile*)pm0ac->Clone();
       h->pm2ac = (TProfile*)pm2ac->Clone();
       h->pmnac = (TProfile*)pmnac->Clone();
       h->pmuac = (TProfile*)pmuac->Clone();
       h->pmoac = (TProfile*)pmoac->Clone();
       
       h->pm0af = (TProfile*)pm0af->Clone();
       h->pm2af = (TProfile*)pm2af->Clone();
       h->pmnaf = (TProfile*)pmnaf->Clone();
       h->pmuaf = (TProfile*)pmuaf->Clone();
       h->pmoaf = (TProfile*)pmoaf->Clone();
       
       h->h2m0ab = (TH2D*)h2m0ab->Clone();
       h->h2m0abx = (TH2D*)h2m0abx->Clone();
       h->h2m2ab = (TH2D*)h2m2ab->Clone();
       h->h2m2abx = (TH2D*)h2m2abx->Clone();
       
       h->p2etaphia = (TProfile2D*)p2etaphia->Clone();
       //}
       //if (pttag >= 40 && pttag <50) {
       h->pm0tb = (TProfile*)pm0tb->Clone();
       h->pm2tb = (TProfile*)pm2tb->Clone();
       h->pmntb = (TProfile*)pmntb->Clone();
       h->pmutb = (TProfile*)pmutb->Clone();
       h->pmotb = (TProfile*)pmotb->Clone();
       
       h->pm0tc = (TProfile*)pm0tc->Clone();
       h->pm2tc = (TProfile*)pm2tc->Clone();
       h->pmntc = (TProfile*)pmntc->Clone();
       h->pmutc = (TProfile*)pmutc->Clone();
       h->pmotc = (TProfile*)pmotc->Clone();
       
       h->h2m0tb = (TH2D*)h2m0tb->Clone();
       h->h2m0tbx = (TH2D*)h2m0tbx->Clone();
       h->h2m2tb = (TH2D*)h2m2tb->Clone();
       h->h2m2tbx = (TH2D*)h2m2tbx->Clone();
       
       h->h2m0tc = (TH2D*)h2m0tc->Clone();
       h->h2m0tcx = (TH2D*)h2m0tcx->Clone();
       h->h2m2tc = (TH2D*)h2m2tc->Clone();
       h->h2m2tcx = (TH2D*)h2m2tcx->Clone();
       
       h->p2etaphit = (TProfile2D*)p2etaphit->Clone();
       //}
       //if (ptprobe >= 40 && ptprobe <50) {
       h->pm0pb = (TProfile*)pm0pb->Clone();
       h->pm2pb = (TProfile*)pm2pb->Clone();
       h->pmnpb = (TProfile*)pmnpb->Clone();
       h->pmupb = (TProfile*)pmupb->Clone();
       h->pmopb = (TProfile*)pmopb->Clone();
       
       h->pm0pf = (TProfile*)pm0pf->Clone(); // fix v9
       h->pm2pf = (TProfile*)pm2pf->Clone();
       h->pmnpf = (TProfile*)pmnpf->Clone();
       h->pmupf = (TProfile*)pmupf->Clone();
       h->pmopf = (TProfile*)pmopf->Clone();
       
       h->h2m0pb = (TH2D*)h2m0pb->Clone();
       h->h2m0pbx = (TH2D*)h2m0pbx->Clone();
       h->h2m2pb = (TH2D*)h2m2pb->Clone();
       h->h2m2pbx = (TH2D*)h2m2pbx->Clone();
       
       h->h2m0pf = (TH2D*)h2m0pf->Clone();
       h->h2m0pfx = (TH2D*)h2m0pfx->Clone();
       h->h2m2pf = (TH2D*)h2m2pf->Clone();
       h->h2m2pfx = (TH2D*)h2m2pfx->Clone();
       
       h->p2etaphip = (TProfile2D*)p2etaphip->Clone();
     } // for ipt
   } // doPtBins

   //Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nentries = fChain->GetEntries();
   cout << "Loaded " << nentries << " entries" << endl;
   //nentries = min((int)nentries, (int)1000000);
   //nentries = min((int)nentries, (int)100000);
   cout << "Processing " << nentries << " entries" << endl;
   //TTimer t;
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;

      //if (jentry%100==0) cout << "." << flush;
      //if (jentry%10000==0) cout << "." << flush;
      if (jentry%100000==0) cout << "." << flush;
      if (jentry%5000000==0) cout << "n="<<jentry<<endl<<flush;
      //if (jentry%1000000==0) cout << timer
     
      if (doJSON) {
	
	//b_branchname->GetEntry(ientry); //read only this branch
	b_run->GetEntry(ientry); //read only this branch
	b_luminosityBlock->GetEntry(ientry); //read only this branch
	
	// Does the run/LS pass the latest JSON selection?
	if (_json[run][luminosityBlock]==0) {
	  //_badjson.insert(pair<int, int>(run, lbn));
	  ++_nbadevts_json;
	  //return false;
	  continue;
	}
      } // doJSON

      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (mrunls.find(run)==mrunls.end()) ++nrun;
      if (mrunls[run].find(luminosityBlock)==mrunls[run].end()) ++nls;
      ++nevt;
      mrunls[run][luminosityBlock] = 1;
      
      vector<float> *vpt = hltPFJetsCorrectedMatchedToCaloJets10_pt;
      vector<float> *veta = hltPFJetsCorrectedMatchedToCaloJets10_eta;
      vector<float> *vphi = hltPFJetsCorrectedMatchedToCaloJets10_phi;
      vector<float> *vmass = hltPFJetsCorrectedMatchedToCaloJets10_mass;

      vector<float> *vptak8 = hltPFJetsCorrectedMatchedToCaloJets10AK8_pt;
      hnjetak8->Fill(vptak8->size());
      double ptak4 = (vpt->size()!=0 ? (*vpt)[0] : 0);
      if (vptak8->size()!=0) hpt0withak8->Fill(ptak4);
      if (vptak8->size()==0) hpt0noak8->Fill(ptak4);

      if (doAK4) {
	if (vptak8->size()!=0) {
	  ++_nbadevts_ak4;
	  continue;
	}
      }

      // Duplicate finding. Single HLT jet matching multiple L1 jets?
      int njet0 = vpt->size();
      set<int> duplicates; // set to work with triplicates also
      for (int i = 0; i != njet0; ++i) {
	for (int j = i+1; j != njet0; ++j) {
	  if ((*vpt)[i]==(*vpt)[j] &&
	      (*veta)[i]==(*veta)[j] &&
	      (*vphi)[i]==(*vphi)[j]) {
	    duplicates.insert(j);
	  }
	} // for j
      } // for i
      // Remove duplicates
      /*
      int ndup(0);
      typedef set<int>::const_iterator IT;
      for (IT it = duplicates.begin(); it != duplicates.end(); ++it) {
	vpt->erase(vpt->begin()+*it-ndup);
	veta->erase(veta->begin()+*it-ndup);
	vphi->erase(vphi->begin()+*it-ndup);
	vmass->erase(vmass->begin()+*it-ndup);
	++ndup;
      }
      */
      // Hopefully faster way to do this (above increased processing >x4)
      // (minimum amount of copying and resizing)
      int ndup = duplicates.size();
      vector<int> vi(njet0-ndup); int idx(0);
      for (int i = 0; i != njet0; ++i) {
	if (duplicates.find(i)==duplicates.end()) {
	  vi[idx] = i; ++idx;
	}
      }
      for (int i = 0; i != njet0-ndup; ++i) {
	if (i!=vi[i]) {
	  (*vpt)[i] = (*vpt)[vi[i]];
	  (*veta)[i] = (*veta)[vi[i]];
	  (*vphi)[i] = (*vphi)[vi[i]];
	  (*vmass)[i] = (*vmass)[vi[i]];
	}
      } // for i
      vpt->resize(njet0-ndup);
      veta->resize(njet0-ndup);
      vphi->resize(njet0-ndup);
      vmass->resize(njet0-ndup);
      
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
      p4mo.SetPtEtaPhiM(0,0,0,0);
      for (int i = 0; i != njet; ++i) {
	
	p4.SetPtEtaPhiM((*vpt)[i], (*veta)[i], (*vphi)[i], (*vmass)[i]);
	
	for (int j = i+1; j != njet; ++j) {
	  p4s.SetPtEtaPhiM((*vpt)[j], (*veta)[j], (*vphi)[j], (*vmass)[j]);
	  double dr = p4.DeltaR(p4s);
	  hdri->Fill(dr);

	  // Dijet mass
	  p4dj = p4; p4dj += p4s;
	  double mjj = p4dj.M();
	  //double mjj = (p4+p4s).M();
	  double deta = fabs(p4.Eta()-p4s.Eta());
	  hmjji->Fill(mjj);
	  if (deta<1.3) hmjji13->Fill(mjj);
	  if (i<2) hmjj1->Fill(mjj);
	  if (i<2 && deta<1.3) hmjj113->Fill(mjj);
	  /*
	  if (dr==0) {
	    h2etaphidr0->Fill(p4.Eta(),p4.Phi());
	    h2ijdr0->Fill(i,j);
	    hpt1dr0->Fill(p4.Pt());
	    hpt2dr0->Fill(p4s.Pt());
	    havdr0->Fill((p4.Pt()-p4s.Pt())/(0.5*(p4.Pt()+p4s.Pt())));
	  }
	  */
	  if (i<2 && j<2) hdri2->Fill(dr);
	  if (fabs(p4.Eta())<1.3 || fabs(p4s.Eta())<1.3) {
	    hdrib->Fill(dr);
	    if (i<2 && j<2) hdri2b->Fill(dr);
	    if (fabs(p4.Eta())<1.3 && fabs(p4s.Eta())<1.3) {
	      hdribb->Fill(dr);
	      if (i<2 && j<2) hdri2bb->Fill(dr);
	    }
	  }
	} // for j

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
	  p4mo -= p4;
	  ++njetn;
	}
	else { // other "unclustered"
	  p4mu -= p4;
	  p4mo -= p4;
	}
      } //  for i in njet
      hnjet->Fill(njet);
      hnjet30->Fill(njet3);
      hndup->Fill(ndup);

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
	  p4b *= 1./p4b.Pt();
	  p4bx.SetPtEtaPhiM(p4b.Pt(),0.,p4b.Phi()+0.5*TMath::Pi(),0.);

	  double eta = p4p.Eta();
	  double pttag = p4t.Pt();
	  double ptprobe = p4p.Pt();
	  double ptave = 0.5*(pttag+ptprobe);
	  double asymm = (ptprobe - pttag) / ptave;
	  // Maximum asymmetry of 2/3 corresponds to x2 ratio of tag and probe
	  // Permit ~0.7 extra scaling to allow for HF L3Res
	  const double maxa = 10;//0.963;//2./3;
	  //double deltaphi = fabs(p4t.DeltaPhi(p4p));
	  double deltaphi = DELTAPHI(p4t.Phi(),p4p.Phi());
	  double dr = p4t.DeltaR(p4p);
	  double mpf = 1 + (p4mht.Vect().Dot(p4b.Vect()))/ptave;
	  double mpf2 = 1 + (p4mht2.Vect().Dot(p4b.Vect()))/ptave;
	  double mpfc = 1 + (p4mhtc.Vect().Dot(p4b.Vect()))/ptave;
	  double mpfc3 = 1 + (p4mhtc3.Vect().Dot(p4b.Vect()))/ptave;

	  // bisector axis
	  double m0b = 1 + (p4m0.Vect().Dot(p4b.Vect()))/ptave;
	  double m2b = 1 + (p4m2.Vect().Dot(p4b.Vect()))/ptave;
	  double mnb = 0 + (p4mn.Vect().Dot(p4b.Vect()))/ptave;
	  double mub = 0 + (p4mu.Vect().Dot(p4b.Vect()))/ptave;
	  double mob = 0 + (p4mo.Vect().Dot(p4b.Vect()))/ptave;

	  double m0bx = 1 + (p4m0.Vect().Dot(p4bx.Vect()))/ptave;
	  double m2bx = 1 + (p4m2.Vect().Dot(p4bx.Vect()))/ptave;

	  // central axis + pttag norm&binning
	  p4c.SetPtEtaPhiM(0,0,0,0);
	  p4c += p4t;
	  p4c.SetPtEtaPhiM(p4c.Pt(),0.,p4c.Phi(),0.);
	  p4c *= 1./p4c.Pt();
	  p4cx.SetPtEtaPhiM(p4c.Pt(),0.,p4c.Phi()+0.5*TMath::Pi(),0.);

	  double m0c = 1 + (p4m0.Vect().Dot(p4c.Vect()))/pttag;
	  double m2c = 1 + (p4m2.Vect().Dot(p4c.Vect()))/pttag;
	  double mnc = 0 + (p4mn.Vect().Dot(p4c.Vect()))/pttag;
	  double muc = 0 + (p4mu.Vect().Dot(p4c.Vect()))/pttag;
	  double moc = 0 + (p4mo.Vect().Dot(p4c.Vect()))/pttag;

	  double m0cx = 1 + (p4m0.Vect().Dot(p4cx.Vect()))/pttag;
	  double m2cx = 1 + (p4m2.Vect().Dot(p4cx.Vect()))/pttag;

	  // forward axis ("backwards" to tag hemisphere) + ptprobe norm&binning
	  p4f.SetPtEtaPhiM(0,0,0,0);
	  p4f -= p4p;
	  p4f.SetPtEtaPhiM(p4f.Pt(),0.,p4f.Phi(),0.);
	  p4f *= 1./p4f.Pt();
	  p4fx.SetPtEtaPhiM(p4f.Pt(),0.,p4f.Phi()+0.5*TMath::Pi(),0.);

	  double m0f = 1 + (p4m0.Vect().Dot(p4f.Vect()))/ptprobe;
	  double m2f = 1 + (p4m2.Vect().Dot(p4f.Vect()))/ptprobe;
	  double mnf = 0 + (p4mn.Vect().Dot(p4f.Vect()))/ptprobe;
	  double muf = 0 + (p4mu.Vect().Dot(p4f.Vect()))/ptprobe;
	  double mof = 0 + (p4mo.Vect().Dot(p4f.Vect()))/ptprobe;

	  double m0fx = 1 + (p4m0.Vect().Dot(p4fx.Vect()))/ptprobe;
	  double m2fx = 1 + (p4m2.Vect().Dot(p4fx.Vect()))/ptprobe;

	  // Dijet mass
	  p4dj = p4t; p4dj += p4p;
	  double mjj = p4dj.M();
	  //double mjj = (p4p+p4t).M();
	  double deta = fabs(p4p.Eta()-p4t.Eta());
	  hmjj2->Fill(mjj);
	  if (deta<1.3) hmjj213->Fill(mjj);

	  basicHistos *h(0);
	  if (fabs(p4t.Eta())<1.3 && deltaphi>2.7 && fabs(asymm)<maxa) {
	    //double eta = p4p.Eta();
	    //if (ptave >= 40 && ptave <50) {
	    int iptave = hptbins->FindBin(ptave)-1;
	    if (iptave>=0 && iptave<npt) {
	      h = vh[iptave];
	      h->pm0ab->Fill(eta, m0b);
	      h->pm2ab->Fill(eta, m2b);
	      h->pmnab->Fill(eta, mnb);
	      h->pmuab->Fill(eta, mub);
	      h->pmoab->Fill(eta, mob);

	      h->pm0ac->Fill(eta, m0c);
	      h->pm2ac->Fill(eta, m2c);
	      h->pmnac->Fill(eta, mnc);
	      h->pmuac->Fill(eta, muc);
	      h->pmoac->Fill(eta, moc);

	      h->pm0af->Fill(eta, m0f);
	      h->pm2af->Fill(eta, m2f);
	      h->pmnaf->Fill(eta, mnf);
	      h->pmuaf->Fill(eta, muf);
	      h->pmoaf->Fill(eta, mof);

	      h->h2m0ab->Fill(eta, m0b);
	      h->h2m0abx->Fill(eta, m0bx);
	      h->h2m2ab->Fill(eta, m2b);
	      h->h2m2abx->Fill(eta, m2bx);

	      h->p2etaphia->Fill(eta, p4p.Phi(), asymm);
	    }
	    //double pttag = p4t.Pt();
	    //if (pttag >= 40 && pttag <50) {
	    int ipttag = hptbins->FindBin(pttag)-1;
	    if (ipttag>=0 && ipttag<npt) {
	      h = vh[ipttag];
	      h->pm0tb->Fill(eta, m0b);
	      h->pm2tb->Fill(eta, m2b);
	      h->pmntb->Fill(eta, mnb);
	      h->pmutb->Fill(eta, mub);
	      h->pmotb->Fill(eta, mob);

	      h->pm0tc->Fill(eta, m0c);
	      h->pm2tc->Fill(eta, m2c);
	      h->pmntc->Fill(eta, mnc);
	      h->pmutc->Fill(eta, muc);
	      h->pmotc->Fill(eta, moc);

	      h->h2m0tb->Fill(eta, m0b);
	      h->h2m0tbx->Fill(eta, m0bx);
	      h->h2m2tb->Fill(eta, m2b);
	      h->h2m2tbx->Fill(eta, m2bx);

	      h->h2m0tc->Fill(eta, m0c);
	      h->h2m0tcx->Fill(eta, m0cx);
	      h->h2m2tc->Fill(eta, m2c);
	      h->h2m2tcx->Fill(eta, m2cx);

	      h->p2etaphit->Fill(eta, p4p.Phi(), p4p.Pt()/pttag);
	    }
	    //double ptprobe = p4p.Pt();
	    //if (ptprobe >= 40 && ptprobe <50) {
	    int iptprobe = hptbins->FindBin(ptprobe)-1;
	    if (iptprobe>=0 && iptprobe<npt) {
	      h = vh[iptprobe];
	      h->pm0pb->Fill(eta, m0b);
	      h->pm2pb->Fill(eta, m2b);
	      h->pmnpb->Fill(eta, mnb);
	      h->pmupb->Fill(eta, mub);
	      h->pmopb->Fill(eta, mob);

	      h->pm0pf->Fill(eta, m0f);
	      h->pm2pf->Fill(eta, m2f);
	      h->pmnpf->Fill(eta, mnf);
	      h->pmupf->Fill(eta, muf);
	      h->pmopf->Fill(eta, mof);

	      h->h2m0pb->Fill(eta, m0b);
	      h->h2m0pbx->Fill(eta, m0bx);
	      h->h2m2pb->Fill(eta, m2b);
	      h->h2m2pbx->Fill(eta, m2bx);

	      h->h2m0pf->Fill(eta, m0f);
	      h->h2m0pfx->Fill(eta, m0fx);
	      h->h2m2pf->Fill(eta, m2f);
	      h->h2m2pfx->Fill(eta, m2fx);

	      h->p2etaphip->Fill(eta, p4p.Phi(), p4t.Pt()/ptprobe);
	    }
	  }

	  if (fabs(p4t.Eta()<1.3) && fabs(asymm)<maxa) {

	    if (ptave>=40) {
	      hdphi->Fill(deltaphi);
	      h2dphi->Fill(p4p.Eta(),deltaphi);
	      hdr->Fill(dr);

	      // Dijet mass
	      hmjj->Fill(mjj);
	      if (deta<1.3) hmjj13->Fill(mjj);

	      if (fabs(p4p.Eta())<1.3) {
		hdphib->Fill(deltaphi);
		hdrb->Fill(dr);
	      }
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

   if (doJSON)
     cout << Form("Found %d bad events according to new JSON (events cut)",_nbadevts_json) << endl;
   if (doAK4)
     cout << Form("Found %d bad events according to new AK8 vector (events cut)",_nbadevts_ak4) << endl;

   cout << "Processed " << nrun << " runs, "
	<< nls << " luminosity blocks and " << nevt << " events" << endl;
   cout << "Saving these to file rootfiles/alcajme.json for brilcalc" << endl;

   ofstream fjson("rootfiles/alcajme.json");
   fjson << "{" << endl;
   for (map<int, map<int,int> >::iterator it = mrunls.begin();
	it != mrunls.end(); ++it) {

     int run = it->first;
     int ls_prev(0), ls(0), ok(0);
     bool firstblock(true), newblock(true);

     if (it!=mrunls.begin()) fjson << "," << endl;
     fjson << "  \"" << run << "\": [";
     for (map<int,int>::iterator jt = it->second.begin();
	  jt != it->second.end(); ++jt) {

       ok = jt->second;
       if (ok) {
	 ls_prev = ls;
	 ls = jt->first;
       }
       else
	 continue;

       newblock = (firstblock || (ls - ls_prev > 1));
       if (newblock && ok) {
	 if (firstblock) {
	   fjson << "[" << ls << ", ";
	   firstblock = false;
	 }
	 else {
	   fjson << ls_prev << "], [" << ls << ", ";
	 }
       }
     } // for jt
     fjson << ls << "]]";
   } // for it
   fjson << "}" << endl;
   
}


bool alcajme::LoadJSON()
{
  // Golden 1.44/fb
  //string json = "rootfiles/Cert_Collisions2022_355100_356615_Golden.json";
  // Golden JSON, 4.86/fb
  string json = "rootfiles/Cert_Collisions2022_355100_357550_Golden..json";
  // DCSOnly 6.40/fb
  //string json = "rootfiles/Cert_Collisions2022_355100_357550_13p6TeV_DCSOnly_TkPx.json";
  bool debug = false;
  cout << "Processing LoadJSON() with " + json + " ..." << endl;
  ifstream file(json, ios::in);
  if (!file.is_open()) return false;
  char c;
  string s, s2;
  char s1[256];
  int rn(0), ls1(0), ls2(0), nrun(0), nls(0);
  file.get(c);
  if (c!='{') return false;
  while (file >> s and sscanf(s.c_str(),"\"%d\":",&rn)==1) {
    if (debug) cout << Form("\"%d\": ",rn);
    
    while (file.get(c) and c==' ') {};
    if (debug) { cout << Form("%c",c) << flush; assert(c=='['); }
    ++nrun;

    bool endrun = false;
    while (!endrun and file >> s >> s2 and sscanf((s+s2).c_str(),"[%d,%d]%s",&ls1,&ls2,s1)==3) {
      //if (debug) cout << Form("[%d,%d,%s]",ls1,ls2,s1);
      if (debug) cout << Form("[%d, %d]%s ",ls1,ls2,s1);

      for (int ls = ls1; ls != ls2+1; ++ls) {
        _json[rn][ls] = 1;
        ++nls;
      }

      s2 = s1;
      endrun = (s2=="]," || s2=="]}");
      if (debug and !endrun and s2!=",") { cout << string("s1: ")+s2 << endl << flush; assert(s2==","); }
    } // while ls
    if (debug) cout << endl;

    if (s2=="]}") continue;
    else if (debug and s2!="],") cout << string("s2: ")+s2 << endl << flush;
    assert(s2=="],");
  } // while run
  //if (s2!="]}") { PrintInfo(string("s3: ")+s2,true); return false; }
  if (s2!="]}") { cout <<  string("s3: ")+s2 << endl; return false; }

  cout << string("Called LoadJSON() with ") + json + ":" << endl;
  cout << Form("Loaded %d good runs and %d good lumi sections",nrun,nls) << endl;
  return true;
} // LoadJSON
