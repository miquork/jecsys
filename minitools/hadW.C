#define hadW_cxx
#include "hadW.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "TProfile.h"
#include "TF1.h"
#include "TLine.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TBox.h"
#include "TGraphErrors.h"
#include "TH2D.h"
#include "TProfile2D.h"

// For JEC
#include "../CondFormats/JetMETObjects/src/Utilities.cc"
#include "../CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "../CondFormats/JetMETObjects/interface/SimpleJetCorrector.h"
#include "../CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
// For JEC uncertainty
//#include "CondFormats/JetMETObjects/interface/SimpleJetCorrectionUncertainty.h"
//#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "../tdrstyle_mod15.C"

#include <vector>

void hadW::Loop()
{
//   In a ROOT session, you can do:
//      root> .L hadW.C
//      root> hadW t
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

   TDirectory *curdir = gDirectory;
   //TFile *fout = new TFile(Form("rootfiles/hadW%s_v2.root",_s.c_str()),
   TFile *fout = new TFile(Form("rootfiles/hadW%s.root",_s.c_str()),
			   "RECREATE");

   const int nbins = 200;
   const double xmin = 30;
   const double xmax = 200;
   const int nx1 = (xmax-xmin)/5.;
   const int nx2 = (xmax-xmin)/10.;
   //const double xb[] = {30,35,40,45,50,60,80,110,150,200};
   //const double xb[] = {30,40,50,62,80,110,150,200}; // UL17?
   const double xb[] = {30,35,40,45,50,60,70,85,105,130,175,230};
   // UL18 Zll+jet as refernce
   // {15, 20, 25, 30, 35, 40, 45, 50, 60, 85, 105, 130, 175, 230,}
   // Z+jet bins (70 extra vs Zll) as reference
   // {15,20,25,30,35,40,50,60,70,85,105,130,175,230,300,400,500,700, 1000, 1500};
   //const double xb[] = {30,42,50,62,80,110,150,200}; // => 42-50 weird
   const int nxb = sizeof(xb)/sizeof(xb[0])-1;
   // For 2D calibration in pT-eta (from deriveL2ResNarrow.C)
   const double yb[] = {0,0.261,0.522,0.783, 1.044,1.305,1.479,
			1.653,1.93,2.172, 2.322,2.5,5.2};//,2.65};
   //2.853,2.964,3.139, 3.489,3.839,5.191};
   const double nyb = sizeof(yb)/sizeof(yb[0])-1;
   // For quick mapping to indeces
   TH1D *hxb = new TH1D("hxb",";p_{T} (GeV)",nxb,xb);
   TH1D *hyb = new TH1D("hyb",";#eta",nyb,yb);
   TH2D *hxbyb = new TH2D("hxbyb",";p_{PT} (GeV);#eta",nxb,xb,nyb,yb);

   double mwmin = 55;//65;
   double mwmax = 105;
   TH1D *hmw = new TH1D("hmw",";m_{W}",nbins,mwmin,mwmax);
   TH1D *hmwveto1b = new TH1D("hmwveto1b",";m_{W}",nbins,mwmin,mwmax);
   TH1D *hmwveto2b = new TH1D("hmwveto2b",";m_{W}",nbins,mwmin,mwmax);
   TH1D *hmwveto1e = new TH1D("hmwveto1e",";m_{W}",nbins,mwmin,mwmax);
   TH1D *hmwveto2e = new TH1D("hmwveto2e",";m_{W}",nbins,mwmin,mwmax);
   //
   TH1D *hmm = new TH1D("hmm",";m_{W,0}",nbins,mwmin,mwmax);
   TH1D *hmmveto1b = new TH1D("hmmveto1b",";m_{W,0}",nbins,mwmin,mwmax);
   TH1D *hmmveto2b = new TH1D("hmmveto2b",";m_{W,0}",nbins,mwmin,mwmax);
   TH1D *hmmveto1e = new TH1D("hmmveto1e",";m_{W,0}",nbins,mwmin,mwmax);
   TH1D *hmmveto2e = new TH1D("hmmveto2e",";m_{W,0}",nbins,mwmin,mwmax);
   //
   TH1D *hpt = new TH1D("hpt",";p_{T,jet}",40,0,200);
   TH1D *hptave = new TH1D("hptave",";p_{T,jet}",40,0,200);
   TH1D *hptmin = new TH1D("hptmin",";p_{T,jet,min}",40,0,200);
   TH2D *hpt2 = new TH2D("hpt2",";p_{T,jet1};p_{T,jet2}",40,0,200,40,0,200);
   TH2D *hpt213 = new TH2D("hpt213",";p_{T,jet1};p_{T,jet2}",40,0,200,40,0,200);
   //
   TH1D *hmw30 = new TH1D("hmw30",";m_{W};N",nbins,mwmin,mwmax);
   TH1D *hmw40 = new TH1D("hmw40",";m_{W};N",nbins,mwmin,mwmax);
   TH1D *hmw50 = new TH1D("hmw50",";m_{W};N",nbins,mwmin,mwmax);
   TH1D *hmw60 = new TH1D("hmw60",";m_{W};N",nbins,mwmin,mwmax);
   //
   TH1D *hmw1330 = new TH1D("hmw1330",";m_{W};N",nbins,mwmin,mwmax);
   TH1D *hmw1340 = new TH1D("hmw1340",";m_{W};N",nbins,mwmin,mwmax);
   TH1D *hmw1350 = new TH1D("hmw1350",";m_{W};N",nbins,mwmin,mwmax);
   TH1D *hmw1360 = new TH1D("hmw1360",";m_{W};N",nbins,mwmin,mwmax);
   //
   TH1D *hmm30 = new TH1D("hmm30",";m_{W,0};N",nbins,mwmin,mwmax);
   TH1D *hmm40 = new TH1D("hmm40",";m_{W};N",nbins,mwmin,mwmax);
   TH1D *hmm50 = new TH1D("hmm50",";m_{W};N",nbins,mwmin,mwmax);
   TH1D *hmm60 = new TH1D("hmm60",";m_{W};N",nbins,mwmin,mwmax);
   //
   TH1D *hmm1330 = new TH1D("hmm1330",";m_{W,0};N",nbins,mwmin,mwmax);
   TH1D *hmm1340 = new TH1D("hmm1340",";m_{W,0};N",nbins,mwmin,mwmax);
   TH1D *hmm1350 = new TH1D("hmm1350",";m_{W,0};N",nbins,mwmin,mwmax);
   TH1D *hmm1360 = new TH1D("hmm1360",";m_{W,0};N",nbins,mwmin,mwmax);
   //
   TH1D *hptjet13 = new TH1D("hptjet13",";p_{T,jet}",nx1,xmin,xmax);
   TH1D *hptave13 = new TH1D("hptave13",";p_{T,ave}",nx1,xmin,xmax);
   TH1D *hptave13b = new TH1D("hptave13b",";p_{T,ave}",nxb,xb);
   TH1D *hptmin13 = new TH1D("hptmin13",";p_{T,min}",nx1,xmin,xmax);
   TH1D *hptthr13 = new TH1D("hpttrh13",";p_{T,thr};N",nx1,xmin,xmax);
   TH1D *hptboth13 = new TH1D("hptboth13",";p_{T,both}",nx2,xmin,xmax);
   TH1D *hptboth13b = new TH1D("hptboth13b",";p_{T,both}",nxb,xb);
   TH1D *hptboth13f = new TH1D("hptboth13f",";p_{T,both}",170,30,200);
   TProfile *pptboth13b = new TProfile("pptboth13b",";p_{T,both};p_{T,both}",
				       nxb,xb);
   TProfile *pptave13b = new TProfile("pptave13b",";p_{T,ave};p_{T,ave}",
				      nxb,xb);
   //
   TProfile *pmw13ptjet = new TProfile("pmw13ptjet",";p_{T,jet};"
				       "#LTm_{W}#GT (GeV)",nx2,xmin,xmax);
   TProfile *pmw13ptmin = new TProfile("pmw13ptmin",";p_{T,min};"
				       "#LTm_{W}#GT (GeV)",nx2,xmin,xmax);
   TProfile *pmw13ptave = new TProfile("pmw13ptave",";p_{T,ave};"
				       "#LTm_{W}#GT (GeV)",nx2,xmin,xmax);
   TProfile *pmw13bptave = new TProfile("pmw13bptave",";p_{T,ave};"
					"#LTm_{W,0}#GT (GeV)",nxb,xb);
   TProfile *pmw13ptthr = new TProfile("pmw13ptthr",";p_{T,thr};"
				       "#LTm_{W}#GT (GeV)",nx2,xmin,xmax);
   TProfile *pmw13ptboth = new TProfile("pmw13ptboth",";p_{T,both};"
					"#LTm_{W}#GT (GeV)",nx2,xmin,xmax);
   TProfile *pmw13bptboth = new TProfile("pmw13bptboth",";p_{T,both};"
					 "#LTm_{W,0}#GT (GeV)",nxb,xb);
   //
   TProfile *pmm13ptjet = new TProfile("pmm13ptjet",";p_{T,jet};"
				       "#LTm_{W,0}#GT (GeV)",nx2,xmin,xmax);
   TProfile *pmm13ptmin = new TProfile("pmm13ptmin",";p_{T,min};"
				       "#LTm_{W,0}#GT (GeV)",nx2,xmin,xmax);
   TProfile *pmm13ptave = new TProfile("pmm13ptave",";p_{T,ave};"
				       "#LTm_{W,0}#GT (GeV)",nx2,xmin,xmax);
   TProfile *pmm13bptave = new TProfile("pmm13bptave",";p_{T,ave};"
					"#LTm_{W,0}#GT (GeV)",nxb,xb);
   TProfile *pmm13bptave_noL3Res = new TProfile("pmm13bptave_noL3Res",
						";p_{T,ave};"
						"#LTm_{W,0,noL3Res}#GT (GeV)",
						nxb,xb);
   TProfile *pmm13ptthr = new TProfile("pmm13ptthr",";p_{T,thr};"
				       "#LTm_{W,0}#GT (GeV)",nx2,xmin,xmax);
   TProfile *pmm13ptboth = new TProfile("pmm13ptboth",";p_{T,both};"
					"#LTm_{W,0}#GT (GeV)",nx2,xmin,xmax);
   TProfile *pmm13bptboth = new TProfile("pmm13bptboth",";p_{T,both};"
					"#LTm_{W,0}#GT (GeV)",nxb,xb);
   TProfile *pmm13bptboth_noL3Res = new TProfile("pmm13bptboth_noL3Res",
						 ";p_{T,both};"
						 "#LTm_{W,0,noL3Res}#GT (GeV)",
						 nxb,xb);
   //
   TProfile *pmm13bptboth_fp0 = new TProfile("pmm13bptboth_fp0",";p_{T,both};"
					     "#LTm_{W,0}#GT (GeV) (fp>0.0)",
					     nxb,xb);
   TProfile *pmm13bptboth_fp1 = new TProfile("pmm13bptboth_fp1",";p_{T,both};"
					     "#LTm_{W,0}#GT (GeV) (fp>0.1)",
					     nxb,xb);
   TProfile *pmm13bptboth_fp2 = new TProfile("pmm13bptboth_fp2",";p_{T,both};"
					     "#LTm_{W,0}#GT (GeV) (fp>0.2)",
					     nxb,xb);
   TProfile *pmm13bptboth_fp3 = new TProfile("pmm13bptboth_fp3",";p_{T,both};"
					     "#LTm_{W,0}#GT (GeV) (fp>0.3)",
					     nxb,xb);
   TProfile *pmm13bptboth_fp4 = new TProfile("pmm13bptboth_fp4",";p_{T,both};"
					     "#LTm_{W,0}#GT (GeV) (fp>0.4)",
					     nxb,xb);
   TProfile *pmm13bptboth_fp5 = new TProfile("pmm13bptboth_fp5",";p_{T,both};"
					     "#LTm_{W,0}#GT (GeV) (fp>0.5)",
					     nxb,xb);
   TProfile *pmm13bptboth_fp6 = new TProfile("pmm13bptboth_fp6",";p_{T,both};"
					     "#LTm_{W,0}#GT (GeV) (fp>0.6)",
					     nxb,xb);
   TProfile *pmm13bptboth_fp7 = new TProfile("pmm13bptboth_fp7",";p_{T,both};"
					     "#LTm_{W,0}#GT (GeV) (fp>0.7)",
					     nxb,xb);
   TProfile *pmm13bptboth_fp8 = new TProfile("pmm13bptboth_fp8",";p_{T,both};"
					     "#LTm_{W,0}#GT (GeV) (fp>0.8)",
					     nxb,xb);
   TProfile *pmm13bptboth_fp9 = new TProfile("pmm13bptboth_fp9",";p_{T,both};"
					     "#LTm_{W,0}#GT (GeV) (fp>0.9)",
					     nxb,xb);
   //
   //
   TProfile *pmm13bptave_fp0 = new TProfile("pmm13bptave_fp0",";p_{T,ave};"
					     "#LTm_{W,0}#GT (GeV) (fp>0.0)",
					     nxb,xb);
   TProfile *pmm13bptave_fp1 = new TProfile("pmm13bptave_fp1",";p_{T,ave};"
					     "#LTm_{W,0}#GT (GeV) (fp>0.1)",
					     nxb,xb);
   TProfile *pmm13bptave_fp2 = new TProfile("pmm13bptave_fp2",";p_{T,ave};"
					     "#LTm_{W,0}#GT (GeV) (fp>0.2)",
					     nxb,xb);
   TProfile *pmm13bptave_fp3 = new TProfile("pmm13bptave_fp3",";p_{T,ave};"
					     "#LTm_{W,0}#GT (GeV) (fp>0.3)",
					     nxb,xb);
   TProfile *pmm13bptave_fp4 = new TProfile("pmm13bptave_fp4",";p_{T,ave};"
					     "#LTm_{W,0}#GT (GeV) (fp>0.4)",
					     nxb,xb);
   TProfile *pmm13bptave_fp5 = new TProfile("pmm13bptave_fp5",";p_{T,ave};"
					     "#LTm_{W,0}#GT (GeV) (fp>0.5)",
					     nxb,xb);
   TProfile *pmm13bptave_fp6 = new TProfile("pmm13bptave_fp6",";p_{T,ave};"
					     "#LTm_{W,0}#GT (GeV) (fp>0.6)",
					     nxb,xb);
   TProfile *pmm13bptave_fp7 = new TProfile("pmm13bptave_fp7",";p_{T,ave};"
					     "#LTm_{W,0}#GT (GeV) (fp>0.7)",
					     nxb,xb);
   TProfile *pmm13bptave_fp8 = new TProfile("pmm13bptave_fp8",";p_{T,ave};"
					     "#LTm_{W,0}#GT (GeV) (fp>0.8)",
					     nxb,xb);
   TProfile *pmm13bptave_fp9 = new TProfile("pmm13bptave_fp9",";p_{T,ave};"
					     "#LTm_{W,0}#GT (GeV) (fp>0.9)",
					     nxb,xb);
   //
   TH2D *hptetaboth13 = new TH2D("hptetaboth13",";p_{T,both};#eta_{jet}",
				 //24,30,150,26,-1.3,1.3);
				 //12,30,150,
				 nx2,xmin,xmax,26,-1.3,1.3);
   TH2D *hptetaboth13b = new TH2D("hptetaboth13b",";p_{T,both};#eta_{jet}",
				  nxb,xb,26,-1.3,1.3);


   // 2D maps
   //const int n2b = (nxb+2)*(nyb+2);
   //const int n2b = nxb*nyb; // excluding over/underflow bins (assert)
   const int n2b = (nxb+1)*nyb; // keeping overflow bin for pt only
   TH2D *h2 = new TH2D("h2",";bin i (fwd);bin j (cnt)",n2b,0,n2b,n2b,0,n2b);
   TProfile2D *p2 = new TProfile2D("p2",";bin i (fwd);bin j (cnt);"
				   "m_{W,0}/m_{W,PDG}",n2b,0,n2b,n2b,0,n2b);
   TProfile2D *p2m2 = new TProfile2D("p2m2",";bin i (fwd);bin j (cnt);"
				     "(m_{W,0}/m_{W,PDG})^{2}",
				     n2b,0,n2b,n2b,0,n2b);

   curdir->cd();

   Long64_t nentries = fChain->GetEntriesFast();

   /*
   fChain->SetBranchStatus("*",0);  // disable all branches
   fChain->SetBranchStatus("recoWMass",1);  // activate branchname
   fChain->SetBranchStatus("fitProb",1);  // activate branchname
   fChain->SetBranchStatus("pt1",1);  // activate branchname
   fChain->SetBranchStatus("pt2",1);  // activate branchname
   fChain->SetBranchStatus("eta1",1);  // activate branchname
   fChain->SetBranchStatus("eta2",1);  // activate branchname
   */

   TBox HBPw89(0,2.705,1.4835,3.1416);
   TBox HEP17(1.31,-0.5236,2.96,-0.8727);

   TLorentzVector j1, j2;
   
   const char *cd = "rootfiles/flavor";
   const char *cfud = "Autumn18_V3_MC_Pythia8_ud_L2Relative_AK4PFchs";
   const char *cfs = "Autumn18_V3_MC_Pythia8_s_L2Relative_AK4PFchs";
   const char *cfc = "Autumn18_V3_MC_Pythia8_c_L2Relative_AK4PFchs";
   const char *cf0 = "Autumn18_V3_MC_Pythia8_all_L2Relative_AK4PFchs";
   const char *cdref = "CondFormats/JetMETObjects/data";
   //const char *cref = "Summer19UL17_RunBCDEF_V2M5_L2L3Residual_AK4PFchs";
   const char *cref = "Summer19UL18_Run%s_V3_DATA_L2Residual_AK4PFchs";
   string sfud = Form("%s/%s.txt",cd,cfud);
   cout << sfud << endl;
   string sfs = Form("%s/%s.txt",cd,cfs);
   cout << sfs << endl;
   string sfc = Form("%s/%s.txt",cd,cfc);
   cout << sfc << endl;
   string s0 = Form("%s/%s.txt",cd,cf0);
   cout << s0 << endl;
   //string sref = Form("%s/%s.txt",cdref,cref);
   //cout << sref << endl;
   vector<JetCorrectorParameters> vfud;
   JetCorrectorParameters *pfud = new JetCorrectorParameters(sfud);
   vfud.push_back(*pfud);
   FactorizedJetCorrector *jecfud = new FactorizedJetCorrector(vfud);
   vector<JetCorrectorParameters> vfs;
   JetCorrectorParameters *pfs = new JetCorrectorParameters(sfs);
   vfs.push_back(*pfs);
   FactorizedJetCorrector *jecfs = new FactorizedJetCorrector(vfs);
   vector<JetCorrectorParameters> vfc;
   JetCorrectorParameters *pfc = new JetCorrectorParameters(sfc);
   vfc.push_back(*pfc);
   FactorizedJetCorrector *jecfc = new FactorizedJetCorrector(vfc);
   //
   vector<JetCorrectorParameters> v0;
   JetCorrectorParameters *pfl0 = new JetCorrectorParameters(s0);
   v0.push_back(*pfl0);
   FactorizedJetCorrector *jec0 = new FactorizedJetCorrector(v0);
   //
   //const int nrun = 5;
   //const char* runs[nrun] = {"B","C","D","E","F"};
   //double runlums[nrun] = {4.8,9.6,4.2,9.3,13.4};
   const int nrun = 4;
   const char* runs[nrun] = {"A","B","C","D"};
   double runlums[nrun] = {14.0, 7.1, 6.9, 31.9};
   vector<FactorizedJetCorrector*> jecrefs(nrun);
   for (int irun = 0; irun != nrun; ++irun) {
     string srun = Form(cref,runs[irun]);
     const char *cref = srun.c_str();
     string sref = Form("%s/%s.txt",cdref,cref);
     cout << sref << endl;
     vector<JetCorrectorParameters> vref;
     JetCorrectorParameters *pref = new JetCorrectorParameters(sref);
     vref.push_back(*pref);
     FactorizedJetCorrector *jecref = new FactorizedJetCorrector(vref);
     jecrefs[irun] = jecref;
   }

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (jentry%100000==00) cout << "." << flush;

      double recoWmassOrig = recoWMass;
      double recoWmass = recoWMass;
      // Alternative with massless jets
      double m0 = 0;
      // or with quark masses (mostly charm mass)
      double mq = 0;//(0.5*0.03+0.25*0.095+0.25*1.275);
      // or with QCD jet mass (guesstimate sqrt(pt)) => too high, too uncertain
      double mj1 = sqrt(pt1);
      double mj2 = sqrt(pt2);
      j1.SetPtEtaPhiM(pt1,eta1,phi1,mq);
      j2.SetPtEtaPhiM(pt2,eta2,phi2,mq);
      double recoWmass0 = (j1+j2).M();

      // Recalculate JEC for light quarks
      jecfud->setJetPt(pt1);
      jecfud->setJetEta(eta1);
      double jecfud1 = jecfud->getCorrection();
      jecfs->setJetPt(pt1);
      jecfs->setJetEta(eta1);
      double jecfs1 = jecfs->getCorrection();
      jecfc->setJetPt(pt1);
      jecfc->setJetEta(eta1);
      double jecfc1 = jecfc->getCorrection();
      double jecfl1 = (0.50*jecfud1+0.25*jecfs1+0.25*jecfc1);
      jec0->setJetPt(pt1);
      jec0->setJetEta(eta1);
      double jec01 = jec0->getCorrection();
      //
      jecfud->setJetPt(pt2);
      jecfud->setJetEta(eta2);
      double jecfud2 = jecfud->getCorrection();
      jecfs->setJetPt(pt2);
      jecfs->setJetEta(eta2);
      double jecfs2 = jecfs->getCorrection();
      jecfc->setJetPt(pt1);
      jecfc->setJetEta(eta1);
      double jecfc2 = jecfc->getCorrection();
      double jecfl2 = (0.50*jecfud2+0.25*jecfs2+0.25*jecfc2);
      jec0->setJetPt(pt2);
      jec0->setJetEta(eta2);
      double jec02 = jec0->getCorrection();
      j1.SetPtEtaPhiM(jecfl1/jec01*pt1,eta1,phi1,mq);
      j2.SetPtEtaPhiM(jecfl2/jec02*pt2,eta2,phi2,mq);
      double recoWmassQ = (j1+j2).M();

      // Also account for neutrino pT lost by charm quarks,
      // which is ~6% of charm quark pT for semileptonic decays (2*16.0%)?
      double jesnu = 1-0.25*0.32*0.06;
      j1.SetPtEtaPhiM(jecfl1/jec01*pt1/jesnu,eta1,phi1,mq);
      j2.SetPtEtaPhiM(jecfl2/jec02*pt2/jesnu,eta2,phi2,mq);
      double recoWmassQNU = (j1+j2).M();

      // Calculate recoWMassQNU without L3Res
      double jesref1(0), jesref2(0), lumsum(0);
      for (int irun = 0; irun != nrun; ++irun) {
	double lumi = runlums[irun];
	lumsum += lumi;
	FactorizedJetCorrector *jecref = jecrefs[irun];
	jecref->setJetPt(pt1);
	jecref->setJetEta(eta1);
	jesref1 += (_s=="MC17" ? 1 : 1./jecref->getCorrection()) * lumi;
	jecref->setJetPt(pt2);
	jecref->setJetEta(eta2);
	jesref2 += (_s=="MC17" ? 1 : 1./jecref->getCorrection()) * lumi;
      }
      jesref1 /= lumsum;
      jesref2 /= lumsum;
      j1.SetPtEtaPhiM(jecfl1/jec01*pt1/jesnu*jesref1,eta1,phi1,mq);
      j2.SetPtEtaPhiM(jecfl2/jec02*pt2/jesnu*jesref2,eta2,phi2,mq);
      double recoWmassUnc = (j1+j2).M();

      // Replacements
      recoWmass = recoWmass0;
      //recoWmass0 = recoWmassQ;
      recoWmass0 = recoWmassQNU;

      if (fitProb>0.2) {
	
	hmw->Fill(recoWmass);
	hmm->Fill(recoWmass0);
	hpt->Fill(pt1);
	hpt->Fill(pt2);
	hpt2->Fill(pt1,pt2);
	if (fabs(eta1)<1.3 && fabs(eta2)<1.3) hpt213->Fill(pt1,pt2);
	hptave->Fill(0.5*(pt1+pt2));
	hptmin->Fill(min(pt1,pt2));
	
	// Check HBPw8/9 veto region
	if ((eta1>HBPw89.GetX1() && eta1<HBPw89.GetX2() &&
	     phi1>HBPw89.GetY1() && phi1<HBPw89.GetY2()) ||
	    (eta2>HBPw89.GetX1() && eta2<HBPw89.GetX2() &&
	     phi2>HBPw89.GetY1() && phi2<HBPw89.GetY2())) {
	  hmwveto1b->Fill(recoWmass);
	  hmmveto1b->Fill(recoWmass0);
	}
	if ((eta1>HBPw89.GetX1() && eta1<HBPw89.GetX2() &&
	     phi1>HBPw89.GetY1() && phi1<HBPw89.GetY2()) &&
	    (eta2>HBPw89.GetX1() && eta2<HBPw89.GetX2() &&
	     phi2>HBPw89.GetY1() && phi2<HBPw89.GetY2())) {
	  hmwveto2b->Fill(recoWmass);
	  hmmveto2b->Fill(recoWmass0);
	}
	// Check HEP17 veto region
	if ((eta1>HEP17.GetX1() && eta1<HEP17.GetX2() &&
	     phi1>HEP17.GetY1() && phi1<HEP17.GetY2()) ||
	    (eta2>HEP17.GetX1() && eta2<HEP17.GetX2() &&
	     phi2>HEP17.GetY1() && phi2<HEP17.GetY2())) {
	  hmwveto1e->Fill(recoWmass);
	  hmmveto1e->Fill(recoWmass0);
	}
	if ((eta1>HEP17.GetX1() && eta1<HEP17.GetX2() &&
	     phi1>HEP17.GetY1() && phi1<HEP17.GetY2()) &&
	    (eta2>HEP17.GetX1() && eta2<HEP17.GetX2() &&
	     phi2>HEP17.GetY1() && phi2<HEP17.GetY2())) {
	  hmwveto2e->Fill(recoWmass);
	  hmmveto2e->Fill(recoWmass0);
	}

	if (pt1>30 && pt2>30) hmw30->Fill(recoWmass);
	if (pt1>40 && pt2>40) hmw40->Fill(recoWmass);
	if (pt1>50 && pt2>50) hmw50->Fill(recoWmass);
	if (pt1>60 && pt2>60) hmw60->Fill(recoWmass);
	//
	if (pt1>30 && pt2>30) hmm30->Fill(recoWmass0);
	if (pt1>40 && pt2>40) hmm40->Fill(recoWmass0);
	if (pt1>50 && pt2>50) hmm50->Fill(recoWmass0);
	if (pt1>60 && pt2>60) hmm60->Fill(recoWmass0);

	bool goodw = (recoWmass>65 && recoWmass<105);
	bool goodw0 = (recoWmass0>60 && recoWmass0<100);
	if (fabs(eta1)<1.3 && fabs(eta2)<1.3) {
	  if (pt1>30 && pt2>30) hmw1330->Fill(recoWmass);
	  if (pt1>40 && pt2>40) hmw1340->Fill(recoWmass);
	  if (pt1>50 && pt2>50) hmw1350->Fill(recoWmass);
	  if (pt1>60 && pt2>60) hmw1360->Fill(recoWmass);
	  //
	  if (pt1>30 && pt2>30) hmm1330->Fill(recoWmass0);
	  if (pt1>40 && pt2>40) hmm1340->Fill(recoWmass0);
	  if (pt1>50 && pt2>50) hmm1350->Fill(recoWmass0);
	  if (pt1>60 && pt2>60) hmm1360->Fill(recoWmass0);

	  if (goodw || goodw0) {
	    hptjet13->Fill(pt1);
	    hptjet13->Fill(pt2);
	    hptave13->Fill(0.5*(pt1+pt2));
	    hptave13b->Fill(0.5*(pt1+pt2));
	    pptave13b->Fill(0.5*(pt1+pt2),0.5*(pt1+pt2));
	    hptmin13->Fill(min(pt1,pt2));
	    //
	    if (goodw) pmw13ptjet->Fill(pt1,recoWmass);
	    if (goodw) pmw13ptjet->Fill(pt2,recoWmass);
	    if (goodw) pmw13ptave->Fill(0.5*(pt1+pt2),recoWmass);
	    if (goodw) pmw13bptave->Fill(0.5*(pt1+pt2),recoWmass);
	    if (goodw) pmw13ptmin->Fill(min(pt1,pt2),recoWmass);
	    //
	    if (goodw0) pmm13ptjet->Fill(pt1,recoWmass0);
	    if (goodw0) pmm13ptjet->Fill(pt2,recoWmass0);
	    if (goodw0) pmm13ptave->Fill(0.5*(pt1+pt2),recoWmass0);
	    if (goodw0) pmm13bptave->Fill(0.5*(pt1+pt2),recoWmass0);
	    if (goodw0) pmm13bptave_noL3Res->Fill(0.5*(pt1+pt2),recoWmassUnc);
	    if (goodw0) pmm13ptmin->Fill(min(pt1,pt2),recoWmass0);

	    for (int i = 1; i != hptboth13->GetNbinsX()+1; ++i) {
	      double ptmin = hptboth13->GetBinLowEdge(i);
	      double ptmax = hptboth13->GetBinLowEdge(i+1);
	      double ptmid = 0.5*(ptmin+ptmax);
	      if (pt1>ptmin && pt2>ptmin && pt1<ptmax && pt2<ptmax) {
		hptboth13->Fill(ptmid);
		if (goodw) pmw13ptboth->Fill(ptmid,recoWmass);
		if (goodw0) pmm13ptboth->Fill(ptmid,recoWmass0);
		if (goodw0) hptetaboth13->Fill(ptmid, eta1);
		if (goodw0) hptetaboth13->Fill(ptmid, eta2);
	      }
	    } // for i
	    for (int i = 1; i != hptboth13b->GetNbinsX()+1; ++i) {
	      double ptmin = hptboth13b->GetBinLowEdge(i);
	      double ptmax = hptboth13b->GetBinLowEdge(i+1);
	      double ptmid = 0.5*(ptmin+ptmax);
	      double ptave = 0.5*(pt1+pt2);
	      if (pt1>ptmin && pt2>ptmin && pt1<ptmax && pt2<ptmax) {
		hptboth13b->Fill(ptmid);
		hptboth13f->Fill(ptave);
		pptboth13b->Fill(ptmid,ptave);
		if (goodw) pmw13bptboth->Fill(ptmid,recoWmass);
		if (goodw0) pmm13bptboth->Fill(ptmid,recoWmass0);
		if (goodw0) pmm13bptboth_noL3Res->Fill(ptmid,recoWmassUnc);
		hptetaboth13b->Fill(ptmid, eta1);
		hptetaboth13b->Fill(ptmid, eta2);
	      }
	    } // for i
	    for (int i = 1; i != hptthr13->GetNbinsX()+1; ++i) {
	      double ptmin = hptthr13->GetBinCenter(i);
	      double ptmax = xmax;//150.;
	      if (pt1>ptmin && pt2>ptmin &&
		  pt1<ptmax && pt2<ptmax) {
		hptthr13->Fill(ptmin);
		if (goodw) pmw13ptthr->Fill(ptmin,recoWmass);
		if (goodw0) pmm13ptthr->Fill(ptmin,recoWmass0);
	      }
	    } // for i
	  } // Wmass
	} // barrel

	// 2D results for fitProb>0.2
	int ix1 = hxb->FindBin(pt1);
	int ix2 = hxb->FindBin(pt2);
	int iy1 = hyb->FindBin(fabs(eta1));
	int iy2 = hyb->FindBin(fabs(eta2));
	// Remember range 0=overflow, 1,...,n, n+1 = overflow => n+2 bins
	//int ib1 = (nxb+2)*iy1 + ix1;
	//int ib2 = (nxb+2)*iy2 + ix2;
	// Keep overflow for pT
	assert(ix1>0); //assert(ix1-1<nxb);
	assert(ix2>0); //assert(ix2-1<nxb);
	assert(iy1>0); assert(iy1-1<nyb);
	assert(iy2>0); assert(iy2-1<nyb);
	int ib1 = (nxb+1)*(iy1-1) + (ix1-1);
	int ib2 = (nxb+1)*(iy2-1) + (ix2-1);
	// Symmetric filling of jets
	const double mwpdg = 80.4;
	//if (goodw0 && fitProb>0.2) {
	//h2->Fill(ib1, ib2);
	//h2->Fill(ib2, ib1);
	//p2->Fill(ib1, ib2, recoWmass0 / mwpdg);
	//p2->Fill(ib2, ib1, recoWmass0 / mwpdg);
	//}
	// Rapidity ordering of jets
	int ifwd = (fabs(eta1)>=fabs(eta2) ? ib1 : ib2);
	int icnt = (fabs(eta1)>=fabs(eta2) ? ib2 : ib1);
	if (goodw0 && fitProb>0.2) {
	  h2->Fill(ifwd, icnt);
	  p2->Fill(ifwd, icnt, recoWmass0 / mwpdg);
	  p2m2->Fill(ifwd, icnt, pow(recoWmass0 / mwpdg,2));
	}
      } // fitProb>0.2
   

      // Fit probability variations
      if (fabs(eta1)<1.3 && fabs(eta2)<1.3) {
	for (int i = 1; i != pmm13bptboth_fp1->GetNbinsX()+1; ++i) {
	  double ptmin = pmm13bptboth_fp1->GetBinLowEdge(i);
	  double ptmax = pmm13bptboth_fp1->GetBinLowEdge(i+1);
	  double ptmid = 0.5*(ptmin+ptmax);
	  bool goodw0 = (recoWmass0>60 && recoWmass0<100);
	  if (pt1>ptmin && pt2>ptmin && pt1<ptmax && pt2<ptmax) {
	    if (fitProb>0.0 && goodw0) 
	      pmm13bptboth_fp0->Fill(ptmid,recoWmass0);
	    if (fitProb>0.1 && goodw0) 
	      pmm13bptboth_fp1->Fill(ptmid,recoWmass0);
	    if (fitProb>0.2 && goodw0) 
	      pmm13bptboth_fp2->Fill(ptmid,recoWmass0);
	    if (fitProb>0.3 && goodw0) 
	      pmm13bptboth_fp3->Fill(ptmid,recoWmass0);
	    if (fitProb>0.4 && goodw0)
	      pmm13bptboth_fp4->Fill(ptmid,recoWmass0);
	    if (fitProb>0.5 && goodw0)
	      pmm13bptboth_fp5->Fill(ptmid,recoWmass0);
	    if (fitProb>0.6 && goodw0)
	      pmm13bptboth_fp6->Fill(ptmid,recoWmass0);
	    if (fitProb>0.7 && goodw0)
	      pmm13bptboth_fp7->Fill(ptmid,recoWmass0);
	    if (fitProb>0.8 && goodw0)
	      pmm13bptboth_fp8->Fill(ptmid,recoWmass0);
	    if (fitProb>0.9 && goodw0)
	      pmm13bptboth_fp9->Fill(ptmid,recoWmass0);
	  } // ptboth
	  //
	  double ptave = 0.5*(pt1+pt2);
	  if (fitProb>0.0 && goodw0) 
	    pmm13bptave_fp0->Fill(ptave,recoWmass0);
	  if (fitProb>0.1 && goodw0) 
	    pmm13bptave_fp1->Fill(ptave,recoWmass0);
	  if (fitProb>0.2 && goodw0) 
	    pmm13bptave_fp2->Fill(ptave,recoWmass0);
	  if (fitProb>0.3 && goodw0) 
	    pmm13bptave_fp3->Fill(ptave,recoWmass0);
	  if (fitProb>0.4 && goodw0)
	    pmm13bptave_fp4->Fill(ptave,recoWmass0);
	  if (fitProb>0.5 && goodw0)
	    pmm13bptave_fp5->Fill(ptave,recoWmass0);
	  if (fitProb>0.6 && goodw0)
	    pmm13bptave_fp6->Fill(ptave,recoWmass0);
	  if (fitProb>0.7 && goodw0)
	    pmm13bptave_fp7->Fill(ptave,recoWmass0);
	  if (fitProb>0.8 && goodw0)
	    pmm13bptave_fp8->Fill(ptave,recoWmass0);
	  if (fitProb>0.9 && goodw0)
	    pmm13bptave_fp9->Fill(ptave,recoWmass0);
	} // for i
      } // fitProb>0.1



   } // for jentry
   
   // Post-processing
   TH1D *ah[] = {hmw, hmw30, hmw40, hmw50, hmw60,
		 hmm, hmm30, hmm40, hmm50, hmm60,
		 hptthr13, 
		 hmw1330, hmw1340, hmw1350, hmw1360,
		 hmm1330, hmm1340, hmm1350, hmm1360,
		 hmwveto1b,hmwveto2b,hmwveto1e,hmwveto2e,
		 hmmveto1b,hmmveto2b,hmmveto1e,hmmveto2e};
   const int nh = sizeof(ah)/sizeof(ah[0]);
   for (int i = 0; i != nh; ++i) {
     TH1D *h = ah[i];
     h->Scale(1./h->Integral());
   }

   // Save histograms
   fout->cd();
   h2->Write("h2",TObject::kOverwrite);
   p2->Write("p2",TObject::kOverwrite);
   p2m2->Write("p2m2",TObject::kOverwrite);
   hxb->Write("hxb",TObject::kOverwrite);
   hyb->Write("hyb",TObject::kOverwrite);
   hxbyb->Write("hxbyb",TObject::kOverwrite);
   //
   hmw->Write("hmw",TObject::kOverwrite);
   hmwveto1b->Write("hmwveto1b",TObject::kOverwrite);
   hmwveto2b->Write("hmwveto2b",TObject::kOverwrite);
   hmwveto1e->Write("hmwveto1e",TObject::kOverwrite);
   hmwveto2e->Write("hmwveto2e",TObject::kOverwrite);
   //
   hmm->Write("hmm",TObject::kOverwrite);
   hmmveto1b->Write("hmmveto1b",TObject::kOverwrite);
   hmmveto2b->Write("hmmveto2b",TObject::kOverwrite);
   hmmveto1e->Write("hmmveto1e",TObject::kOverwrite);
   hmmveto2e->Write("hmmveto2e",TObject::kOverwrite);
   //
   hpt->Write("hpt",TObject::kOverwrite);
   hptave->Write("hptave",TObject::kOverwrite);
   hptmin->Write("hptmin",TObject::kOverwrite);
   hpt2->Write("hpt2",TObject::kOverwrite);
   hpt213->Write("hpt213",TObject::kOverwrite);
   //
   hmw30->Write("hmw30",TObject::kOverwrite);
   hmw40->Write("hmw40",TObject::kOverwrite);
   hmw50->Write("hmw50",TObject::kOverwrite);
   hmw60->Write("hmw60",TObject::kOverwrite);
   hmw1330->Write("hmw1330",TObject::kOverwrite);
   hmw1340->Write("hmw1340",TObject::kOverwrite);
   hmw1350->Write("hmw1350",TObject::kOverwrite);
   hmw1360->Write("hmw1360",TObject::kOverwrite);
   //
   hmm30->Write("hmm30",TObject::kOverwrite);
   hmm40->Write("hmm40",TObject::kOverwrite);
   hmm50->Write("hmm50",TObject::kOverwrite);
   hmm60->Write("hmm60",TObject::kOverwrite);
   hmm1330->Write("hmm1330",TObject::kOverwrite);
   hmm1340->Write("hmm1340",TObject::kOverwrite);
   hmm1350->Write("hmm1350",TObject::kOverwrite);
   hmm1360->Write("hmm1360",TObject::kOverwrite);
   //
   hptjet13->Write("hptjet13",TObject::kOverwrite);
   hptmin13->Write("hptmin13",TObject::kOverwrite);
   hptave13->Write("hptave13",TObject::kOverwrite);
   hptave13b->Write("hptave13b",TObject::kOverwrite);
   hptthr13->Write("hptthr13",TObject::kOverwrite);
   hptboth13->Write("hptboth13",TObject::kOverwrite);
   hptboth13b->Write("hptboth13b",TObject::kOverwrite);
   hptboth13f->Write("hptboth13f",TObject::kOverwrite);
   pptboth13b->Write("pptboth13b",TObject::kOverwrite);
   pptave13b->Write("pptave13b",TObject::kOverwrite);
   //
   pmw13ptjet->Write("pmw13ptjet",TObject::kOverwrite);
   pmw13ptmin->Write("pmw13ptmin",TObject::kOverwrite);
   pmw13ptave->Write("pmw13ptave",TObject::kOverwrite);
   pmw13bptave->Write("pmw13bptave",TObject::kOverwrite);
   pmw13ptthr->Write("pmw13ptthr",TObject::kOverwrite);
   pmw13ptboth->Write("pmw13ptboth",TObject::kOverwrite);
   pmw13bptboth->Write("pmw13bptboth",TObject::kOverwrite);
   //
   pmm13ptjet->Write("pmm13ptjet",TObject::kOverwrite);
   pmm13ptmin->Write("pmm13ptmin",TObject::kOverwrite);
   pmm13ptave->Write("pmm13ptave",TObject::kOverwrite);
   pmm13bptave->Write("pmm13bptave",TObject::kOverwrite);
   pmm13bptave_noL3Res->Write("pmm13bptave_noL3Res",TObject::kOverwrite);
   pmm13ptthr->Write("pmm13ptthr",TObject::kOverwrite);
   pmm13ptboth->Write("pmm13ptboth",TObject::kOverwrite);
   pmm13bptboth->Write("pmm13bptboth",TObject::kOverwrite);
   pmm13bptboth_noL3Res->Write("pmm13bptboth_noL3Res",TObject::kOverwrite);
   //
   pmm13bptboth_fp0->Write("pmm13bptboth_fp0",TObject::kOverwrite);
   pmm13bptboth_fp1->Write("pmm13bptboth_fp1",TObject::kOverwrite);
   pmm13bptboth_fp2->Write("pmm13bptboth_fp2",TObject::kOverwrite);
   pmm13bptboth_fp3->Write("pmm13bptboth_fp3",TObject::kOverwrite);
   pmm13bptboth_fp4->Write("pmm13bptboth_fp4",TObject::kOverwrite);
   pmm13bptboth_fp5->Write("pmm13bptboth_fp5",TObject::kOverwrite);
   pmm13bptboth_fp6->Write("pmm13bptboth_fp6",TObject::kOverwrite);
   pmm13bptboth_fp7->Write("pmm13bptboth_fp7",TObject::kOverwrite);
   pmm13bptboth_fp8->Write("pmm13bptboth_fp8",TObject::kOverwrite);
   pmm13bptboth_fp9->Write("pmm13bptboth_fp9",TObject::kOverwrite);
   //
   pmm13bptave_fp0->Write("pmm13bptave_fp0",TObject::kOverwrite);
   pmm13bptave_fp1->Write("pmm13bptave_fp1",TObject::kOverwrite);
   pmm13bptave_fp2->Write("pmm13bptave_fp2",TObject::kOverwrite);
   pmm13bptave_fp3->Write("pmm13bptave_fp3",TObject::kOverwrite);
   pmm13bptave_fp4->Write("pmm13bptave_fp4",TObject::kOverwrite);
   pmm13bptave_fp5->Write("pmm13bptave_fp5",TObject::kOverwrite);
   pmm13bptave_fp6->Write("pmm13bptave_fp6",TObject::kOverwrite);
   pmm13bptave_fp7->Write("pmm13bptave_fp7",TObject::kOverwrite);
   pmm13bptave_fp8->Write("pmm13bptave_fp8",TObject::kOverwrite);
   pmm13bptave_fp9->Write("pmm13bptave_fp9",TObject::kOverwrite);
   //
   hptetaboth13->Write("hptetaboth13",TObject::kOverwrite);
   hptetaboth13b->Write("hptetaboth13b",TObject::kOverwrite);

   fout->Close();
   curdir->cd();
} // Loop

void hadW::Draw() {//string sdt, string smc) {

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  //TFile *fdt = new TFile("rootfiles/hadWUL17.root","READ");
  TFile *fdt = new TFile("rootfiles/hadWUL18.root","READ");
  assert(fdt && !fdt->IsZombie());
  //TFile *fmc = new TFile("rootfiles/hadWMC17.root","READ");
  TFile *fmc = new TFile("rootfiles/hadWMC18.root","READ");
  assert(fdt && !fdt->IsZombie());
  //TFile *fout = new TFile("rootfiles/hadWUL17.root","RECREATE");
  TFile *fout = new TFile("rootfiles/hadWUL18.root","UPDATE");
  curdir->cd();

  /*
  TH1D *hdt = new TH1D("hdt","Min p_{T,scale};m_{W}",4,25,65);
  TH1D *hmc = new TH1D("hmc","Min p_{T,scale};m_{W}",4,25,65);
  TH1D *hr = new TH1D("hr","Min p_{T,scale};Data/MC-1 (%)",4,25,65);
  for (int i = 1; i != hdt->GetNbinsX()+1; ++i) {
    int ipt = int(hdt->GetBinCenter(i)+0.5);
    TH1D *hd = (TH1D*)fdt->Get(Form("hmw13%d",ipt)); assert(hd);
    //if (!hd) continue;
    hdt->SetBinContent(i, hd->GetMean());
    hdt->SetBinError(i, hd->GetMeanError());
    TH1D *hm = (TH1D*)fmc->Get(Form("hmw13%d",ipt)); assert(hm);
    hmc->SetBinContent(i, hm->GetMean());
    hmc->SetBinError(i, hm->GetMeanError());
    //
    hr->SetBinContent(i, 100.*(hd->GetMean()/hm->GetMean()-1));
    hr->SetBinError(i, 100.*sqrt(pow(hd->GetMeanError()/hd->GetMean(),2)+
				 pow(hm->GetMeanError()/hm->GetMean(),2)));
  }
  */

  TProfile *ptdt = (TProfile*)fdt->Get("pmw13ptthr"); assert(ptdt);
  TProfile *ptmc = (TProfile*)fmc->Get("pmw13ptthr"); assert(ptmc);
  TH1D *ptr = ptdt->ProjectionX("ptr");
  TF1 *f1 = new TF1("f1","1",0,7500);
  ptr->Divide(ptmc);
  ptr->Add(f1,-1);
  ptr->Scale(100.);

  TProfile *pbdt = (TProfile*)fdt->Get("pmw13bptboth"); assert(pbdt);
  TProfile *pbmc = (TProfile*)fmc->Get("pmw13bptboth"); assert(pbmc);
  TH1D *pbr = pbdt->ProjectionX("pbr");
  pbr->Divide(pbmc);
  pbr->Add(f1,-1);
  pbr->Scale(100.);

  TProfile *pt0dt = (TProfile*)fdt->Get("pmm13ptthr"); assert(pt0dt);
  TProfile *pt0mc = (TProfile*)fmc->Get("pmm13ptthr"); assert(pt0mc);
  TH1D *pt0r = pt0dt->ProjectionX("pt0r");
  pt0r->Divide(pt0mc);
  pt0r->Add(f1,-1);
  pt0r->Scale(100.);

  TProfile *pa0dt = (TProfile*)fdt->Get("pmm13bptave"); assert(pa0dt);
  TProfile *pa0mc = (TProfile*)fmc->Get("pmm13bptave"); assert(pa0mc);
  TH1D *pa0r = pa0dt->ProjectionX("pa0r");
  pa0r->Divide(pa0mc);
  pa0r->Add(f1,-1);
  pa0r->Scale(100.);

  TProfile *pb0dt = (TProfile*)fdt->Get("pmm13bptboth"); assert(pb0dt);
  TProfile *pb0mc = (TProfile*)fmc->Get("pmm13bptboth"); assert(pb0mc);
  TH1D *pb0r = pb0dt->ProjectionX("pb0r");
  pb0r->Divide(pb0mc);
  pb0r->Add(f1,-1);
  pb0r->Scale(100.);

  TProfile *pdt = (TProfile*)fdt->Get("pmm13bptboth_noL3Res"); assert(pdt);
  TProfile *pmc = (TProfile*)fmc->Get("pmm13bptboth_noL3Res"); assert(pmc);
  TH1D *pr = pdt->ProjectionX("pr");
  pr->Divide(pmc);
  pr->Add(f1,-1);
  pr->Scale(100.);

  TProfile *padt = (TProfile*)fdt->Get("pmm13bptave_noL3Res"); assert(padt);
  TProfile *pamc = (TProfile*)fmc->Get("pmm13bptave_noL3Res"); assert(pamc);
  TH1D *par = padt->ProjectionX("par");
  par->Divide(pmc);
  par->Add(f1,-1);
  par->Scale(100.);

  TH1D *hdt = (TH1D*)fdt->Get("hptboth13b"); assert(hdt);
  TH1D *hmc = (TH1D*)fmc->Get("hptboth13b"); assert(hmc);
  TH1D *hr = (TH1D*)hdt->Clone("hptboth13br");
  hr->Divide(hmc);

  TH1D *hadt = (TH1D*)fdt->Get("hptave13b"); assert(hadt);
  TH1D *hamc = (TH1D*)fmc->Get("hptave13b"); assert(hamc);
  TH1D *har = (TH1D*)hadt->Clone("hptave13br");
  har->Divide(hamc);

  //TH1D *hup = tdrHist("hup","#LTm_{W}#GT (GeV)",74,102,"p_{T,ref} (GeV)",25,205);
  //TH1D *hdw = tdrHist("hdw","Data/MC-1 (%)",-1.5,+1.5,"p_{T,ref} (GeV)",25,205);
  TH1D *hup =tdrHist("hup","#LTm_{W}#GT (GeV)",73,90,"p_{T,jets} (GeV)",30,200);
  TH1D *hdw =tdrHist("hdw","Data/MC-1 (%)",-1.5,+1.5,"p_{T,jets} (GeV)",30,200);

  //lumi_13TeV = "UL17 ttbar lepton+jet hadW";
  lumi_13TeV = "UL18 ttbar lepton+jet hadW";
  TCanvas *c1 = tdrDiCanvas("c1",hup,hdw,4,11);
  
  TLine *l = new TLine();
  l->SetLineStyle(kDashed);

  c1->cd(1);
  gPad->SetLogx();

  l->DrawLine(30,80.4,200,80.4);
  tdrDraw(ptdt,"Pz",kFullCircle,kRed);
  tdrDraw(ptmc,"Pz",kOpenCircle,kRed);
  tdrDraw(pbdt,"Pz",kFullCircle,kBlue);
  tdrDraw(pbmc,"Pz",kOpenCircle,kBlue);
  //
  tdrDraw(pt0dt,"Pz",kFullCircle,kMagenta+2);
  tdrDraw(pt0mc,"Pz",kOpenCircle,kMagenta+2);
  tdrDraw(pb0dt,"Pz",kFullCircle,kCyan+2);
  tdrDraw(pb0mc,"Pz",kOpenCircle,kCyan+2);

  /*
  TLegend *leg1 = tdrLeg(0.40,0.70,0.60,0.90);
  leg1->SetHeader("m_{j}#neq0");
  leg1->AddEntry(ptdt," ","PLE");
  leg1->AddEntry(ptmc," ","PLE");
  leg1->AddEntry(pbdt," ","PLE");
  leg1->AddEntry(pbmc," ","PLE");

  TLegend *leg2 = tdrLeg(0.48,0.70,0.68,0.90);
  leg2->SetHeader("=0");
  leg2->AddEntry(pt0dt,"Data ref<p_{T}<200 GeV","PLE");
  leg2->AddEntry(pt0mc,"MC   ref<p_{T}<200 GeV","PLE");
  leg2->AddEntry(pb0dt,"Data ref<p_{T}<ref+1","PLE");
  leg2->AddEntry(pb0mc,"MC   ref<p_{T}<ref+1","PLE");
  */
  TLatex *tex = new TLatex();
  tex->SetTextSize(0.045); tex->SetNDC();
  tex->DrawLatex(0.72,0.855,"|#eta_{j}| < 1.3");
  tex->DrawLatex(0.19,0.71,"fitProb>0.2");

  TLegend *leg1 = tdrLeg(0.39,0.64,0.59,0.89);
  leg1->SetHeader("DT");
  leg1->AddEntry(ptdt," ","PLE");
  leg1->AddEntry(pbdt," ","PLE");
  leg1->AddEntry(pt0dt," ","PLE");
  leg1->AddEntry(pb0dt," ","PLE");

  TLegend *leg2 = tdrLeg(0.45,0.64,0.65,0.89);
  leg2->SetHeader("MC");
  /*
  leg2->AddEntry(ptmc,"ref<p_{T,j}<200 GeV, m_{j}#neq0","PLE");
  leg2->AddEntry(pbmc,"ref<p_{T,j}<ref+1, m_{j}#neq0","PLE");
  leg2->AddEntry(pt0mc,"ref<p_{T,j}<200 GeV, m_{j}=0","PLE");
  leg2->AddEntry(pb0mc,"ref<p_{T,j}<ref+1, m_{j}=0","PLE");
  */
  leg2->AddEntry(ptmc,"p_{T,j}<200 GeV, m_{j}=0#otimesR_{dj}","PLE");
  leg2->AddEntry(pbmc,"m_{j}=0 & R_{dj}","PLE");
  leg2->AddEntry(pt0mc,"p_{T,j}<200 GeV, m_{j}=0#otimesR_{q}#timesR_{#nu}","PLE");
  leg2->AddEntry(pb0mc,"m_{j}=0 & R_{udsc}#timesR_{#nu}","PLE");

  c1->cd(2);
  gPad->SetLogx();

  l->DrawLine(30,0,200,0);
  tdrDraw(ptr,"Pz",kFullCircle,kRed);
  tdrDraw(pbr,"Pz",kFullCircle,kBlue);
  //
  tdrDraw(pt0r,"Pz",kFullCircle,kMagenta+2);
  tdrDraw(pb0r,"Pz",kFullCircle,kCyan+2);

  c1->SaveAs("pdf/hadW.pdf");


  //
  TProfile *pptdt = (TProfile*)fdt->Get("pptboth13b"); assert(pptdt);
  TProfile *pptmc = (TProfile*)fmc->Get("pptboth13b"); assert(pptmc);
  TProfile *pptadt = (TProfile*)fdt->Get("pptave13b"); assert(pptadt);
  TProfile *pptamc = (TProfile*)fmc->Get("pptave13b"); assert(pptamc);

  TGraphErrors *gb0dt = new TGraphErrors(0);
  TGraphErrors *gb0mc = new TGraphErrors(0);
  TGraphErrors *gb0r = new TGraphErrors(0);
  TGraphErrors *ga0dt = new TGraphErrors(0);
  TGraphErrors *ga0mc = new TGraphErrors(0);
  TGraphErrors *ga0r = new TGraphErrors(0);
  TGraphErrors *gdt = new TGraphErrors(0);
  TGraphErrors *gmc = new TGraphErrors(0);
  TGraphErrors *gr = new TGraphErrors(0);
  TGraphErrors *gadt = new TGraphErrors(0);
  TGraphErrors *gamc = new TGraphErrors(0);
  TGraphErrors *gar = new TGraphErrors(0);
  for (int i = 1; i != pptdt->GetNbinsX()+1; ++i) {
    double ptdt = pptdt->GetBinContent(i);
    double ptmc = pptmc->GetBinContent(i);
    int j = pb0dt->FindBin(ptdt);
    gb0dt->SetPoint(i-1, ptdt, pb0dt->GetBinContent(j));
    gb0dt->SetPointError(i-1, pptdt->GetBinError(i), pb0dt->GetBinError(j));
    gb0mc->SetPoint(i-1, ptmc, pb0mc->GetBinContent(j));
    gb0mc->SetPointError(i-1, pptmc->GetBinError(i), pb0mc->GetBinError(j));
    gb0r->SetPoint(i-1, 0.5*(ptdt+ptmc),
		   100*(pb0dt->GetBinContent(j) / pb0mc->GetBinContent(j)-1));
    gb0r->SetPointError(i-1, 0.5*(pptdt->GetBinError(i)+pptmc->GetBinError(i)),
			sqrt(pow(pb0dt->GetBinError(j)/pb0dt->GetBinContent(j),2)+pow(pb0mc->GetBinError(j)/pb0mc->GetBinContent(j),2)) *
			pb0dt->GetBinContent(j) / pb0mc->GetBinContent(j)*100);
    //
    double ptadt = pptadt->GetBinContent(i);
    double ptamc = pptamc->GetBinContent(i);
    assert(pa0dt->FindBin(ptadt)==j);
    ga0dt->SetPoint(i-1, ptadt, pa0dt->GetBinContent(j));
    ga0dt->SetPointError(i-1, pptadt->GetBinError(i), pa0dt->GetBinError(j));
    ga0mc->SetPoint(i-1, ptamc, pa0mc->GetBinContent(j));
    ga0mc->SetPointError(i-1, pptamc->GetBinError(i), pa0mc->GetBinError(j));
    ga0r->SetPoint(i-1, 0.5*(ptadt+ptamc),
		   100*(pa0dt->GetBinContent(j) / pa0mc->GetBinContent(j)-1));
    ga0r->SetPointError(i-1,0.5*(pptadt->GetBinError(i)+pptamc->GetBinError(i)),
			sqrt(pow(pa0dt->GetBinError(j)/pa0dt->GetBinContent(j),2)+pow(pa0mc->GetBinError(j)/pa0mc->GetBinContent(j),2)) *
			pa0dt->GetBinContent(j) / pa0mc->GetBinContent(j)*100);
    //
    double mw = 80.4;
    gdt->SetPoint(i-1, ptdt, pdt->GetBinContent(j)/mw);
    gdt->SetPointError(i-1, pptdt->GetBinError(i), pdt->GetBinError(j)/mw);
    gmc->SetPoint(i-1, ptmc, pmc->GetBinContent(j)/mw);
    gmc->SetPointError(i-1, pptmc->GetBinError(i), pmc->GetBinError(j)/mw);
    gr->SetPoint(i-1, 0.5*(ptdt+ptmc),
		 pdt->GetBinContent(j) / pmc->GetBinContent(j));
    //100*(pdt->GetBinContent(j) / pmc->GetBinContent(j)-1));
    //gr->SetPointError(i-1, 0.5*(pptdt->GetBinError(i)+pptmc->GetBinError(i)),
    //		      sqrt(pow(pdt->GetBinError(j)/pdt->GetBinContent(j),2)+pow(pmc->GetBinError(j)/pmc->GetBinContent(j),2)) *
    //		      pdt->GetBinContent(j) / pmc->GetBinContent(j)*100);
    gr->SetPointError(i-1, 0.5*(pptdt->GetBinError(i)+pptmc->GetBinError(i)),
		      sqrt(pow(pdt->GetBinError(j)/pdt->GetBinContent(j),2)+pow(pmc->GetBinError(j)/pmc->GetBinContent(j),2)) *
    		      pdt->GetBinContent(j) / pmc->GetBinContent(j));
    //
    gadt->SetPoint(i-1, ptadt, padt->GetBinContent(j)/mw);
    gadt->SetPointError(i-1, pptadt->GetBinError(i), padt->GetBinError(j)/mw);
    gamc->SetPoint(i-1, ptamc, pamc->GetBinContent(j)/mw);
    gamc->SetPointError(i-1, pptamc->GetBinError(i), pamc->GetBinError(j)/mw);
    gar->SetPoint(i-1, 0.5*(ptadt+ptamc),
		  padt->GetBinContent(j) / pamc->GetBinContent(j));
    gar->SetPointError(i-1, 0.5*(pptadt->GetBinError(i)+pptamc->GetBinError(i)),
		       sqrt(pow(padt->GetBinError(j)/padt->GetBinContent(j),2)+pow(pamc->GetBinError(j)/pamc->GetBinContent(j),2)) *
		       padt->GetBinContent(j) / pamc->GetBinContent(j));
  } // for i

  //TH1D *hup2 = tdrHist("hup2","#LTm_{W}#GT (GeV)",72.9+1e-5,83.9-1e-5,
  //TH1D *hup2 = tdrHist("hup2","#LTm_{W}#GT (GeV)",79.9+1e-5,80.9-1e-5,
  //TH1D *hup2 = tdrHist("hup2","#LTm_{W}#GT (GeV)",79.5+1e-5,81.5-1e-5,
  TH1D *hup2 = tdrHist("hup2","#LTm_{W}#GT (GeV)",79.0+1e-5,82.0-1e-5,
		       "p_{T,jets} (GeV)",30,200);
  TH1D *hdw2 = tdrHist("hdw2","Data/MC-1 (%)",-1.6,+0.8,
		       "p_{T,jets} (GeV)",30,200);

  TCanvas *c2 = tdrDiCanvas("c2",hup2,hdw2,4,11);
  
  c2->cd(1);
  gPad->SetLogx();

  l->DrawLine(30,80.4,200,80.4);
  //tdrDraw(pt0dt,"Pz",kFullCircle,kMagenta+2);
  //tdrDraw(pt0mc,"Pz",kOpenCircle,kMagenta+2);
  tdrDraw(pb0dt,"Pz",kFullCircle,kCyan+2);
  tdrDraw(pb0mc,"Pz",kOpenCircle,kCyan+2);
  tdrDraw(gb0dt,"Pz",kFullCircle,kCyan+3);
  tdrDraw(gb0mc,"Pz",kOpenCircle,kCyan+3);
  tdrDraw(pa0dt,"Pz",kFullDiamond,kMagenta+1);
  tdrDraw(pa0mc,"Pz",kOpenDiamond,kMagenta+2);
  pa0dt->SetMarkerSize(1.5);
  pa0mc->SetMarkerSize(1.5);
  // Check that L3Res was correctly removed
  tdrDraw(gdt,"Pz",kFullCircle,kBlack);
  tdrDraw(gmc,"Pz",kOpenCircle,kBlack);

  gb0dt->GetXaxis()->SetTitle("p_{T,jets} (GeV)");
  gb0mc->GetXaxis()->SetTitle("p_{T,jets} (GeV)");
  gdt->GetXaxis()->SetTitle("p_{T,jets} (GeV)");
  gmc->GetXaxis()->SetTitle("p_{T,jets} (GeV)");
  //
  gb0dt->GetYaxis()->SetTitle("m_{W,0} (GeV)");
  gb0mc->GetYaxis()->SetTitle("m_{W,0} (GeV)");
  gdt->GetYaxis()->SetTitle("m_{W,0}/m_{W,PDG} (GeV)");
  gmc->GetYaxis()->SetTitle("m_{W,0}/m_{W,PDG} (GeV)");

  tex->DrawLatex(0.40,0.84,"fitProb>0.2");
  tex->DrawLatex(0.40,0.78,"|#eta_{j}| < 1.3");
  tex->DrawLatex(0.30,0.72,"60<m_{W}<100 GeV");

  TLegend *leg12 = tdrLeg(0.62,0.72,0.82,0.87);
  leg12->SetHeader("DT");
  //leg12->AddEntry(pt0dt," ","PLE");
  leg12->AddEntry(pb0dt," ","PLE");
  leg12->AddEntry(pa0dt," ","PLE");

  TLegend *leg22 = tdrLeg(0.68,0.72,0.88,0.87);
  leg22->SetHeader("MC");
  //leg22->AddEntry(pt0mc,"ref<p_{T,j}<200 GeV, m_{j}=0","PLE");
  //leg22->AddEntry(pb0mc,"ref<p_{T,j}<ref+1, m_{j}=0","PLE");
  //leg22->AddEntry(pb0mc," m_{j}=0, R_{dj}","PLE");
  leg22->AddEntry(pb0mc," m_{j}=0#otimesR_{q}#timesR_{#nu}","PLE");
  leg22->AddEntry(pa0mc," vs p_{T,ave}","PLE");

  c2->cd(2);
  gPad->SetLogx();

  l->DrawLine(30,0,200,0);
  //tdrDraw(pt0r,"Pz",kFullCircle,kMagenta+2);
  tdrDraw(pb0r,"Pz",kFullCircle,kCyan+2);
  tdrDraw(gb0r,"Pz",kFullCircle,kCyan+3);
  // Check that L3Res was correctly removed
  tdrDraw(gr,"Pz",kFullCircle,kBlack);
  tdrDraw(pa0r,"Pz",kFullDiamond,kMagenta+1);
  tdrDraw(ga0r,"Pz",kFullDiamond,kMagenta+2);
  pa0r->SetMarkerSize(1.5);

  gb0r->GetXaxis()->SetTitle("p_{T,jets} (GeV)");
  gr->GetXaxis()->SetTitle("p_{T,jets} (GeV)");
  //
  gb0r->GetYaxis()->SetTitle("m_{W,0,data}/m_{W,0,MC}-1 (%)");
  gr->GetYaxis()->SetTitle("m_{W,0,data}/m_{W,0,MC}");

  c2->SaveAs("pdf/hadW0.pdf");
  
  fout->cd();
  hdt->Write("data_nevents_ptboth_hadw_fitprob02_L1L2L3");
  hmc->Write("mc_nevents_ptboth_hadw_fitprob02_L1L2L3");
  hr->Write("ratio_nevents_ptboth_hadw_fitprob02_L1L2L3");
  gb0dt->Write("data_mass_ptboth_hadw_fitprob02_L1L2L3Res",TObject::kOverwrite);
  gb0mc->Write("mc_mass_ptboth_hadw_fitprob02_L1L2L3Res",TObject::kOverwrite);
  gb0r->Write("ratio_mass_ptboth_hadw_fitprob02_L1L2L3Res",TObject::kOverwrite);
  gdt->Write("data_mass_ptboth_hadw_fitprob02_L1L2L3",TObject::kOverwrite);
  gmc->Write("mc_mass_ptboth_hadw_fitprob02_L1L2L3",TObject::kOverwrite);
  gr->Write("ratio_mass_ptboth_hadw_fitprob02_L1L2L3",TObject::kOverwrite);
  //
  hadt->Write("data_nevents_ptave_hadw_fitprob02_L1L2L3");
  hamc->Write("mc_nevents_ptave_hadw_fitprob02_L1L2L3");
  har->Write("ratio_nevents_ptave_hadw_fitprob02_L1L2L3");
  ga0dt->Write("data_mass_ptave_hadw_fitprob02_L1L2L3Res",TObject::kOverwrite);
  ga0mc->Write("mc_mass_ptave_hadw_fitprob02_L1L2L3Res",TObject::kOverwrite);
  ga0r->Write("ratio_mass_ptave_hadw_fitprob02_L1L2L3Res",TObject::kOverwrite);
  gadt->Write("data_mass_ptave_hadw_fitprob02_L1L2L3",TObject::kOverwrite);
  gamc->Write("mc_mass_ptave_hadw_fitprob02_L1L2L3",TObject::kOverwrite);
  gar->Write("ratio_mass_ptave_hadw_fitprob02_L1L2L3",TObject::kOverwrite);
  fout->Close();

  // Draw m=0 case only as well? Proper <pT,j> centering
  // Add flavor corrections. Ad-hoc mass correction for PU as well?
  // Matrix case (TH2D of #jets and their mW in neighbouring bins)?

} // draw

// Draw a scan of fitProbability cuts
void hadW::DrawFP(string spt) {

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  //TFile *fdt = new TFile("rootfiles/hadWUL17.root","READ");
  TFile *fdt = new TFile("rootfiles/hadWUL18.root","READ");
  assert(fdt && !fdt->IsZombie());
  //TFile *fmc = new TFile("rootfiles/hadWMC17.root","READ");
  TFile *fmc = new TFile("rootfiles/hadWMC18.root","READ");
  assert(fdt && !fdt->IsZombie());
  curdir->cd();

  const char *cpt = spt.c_str();

  const int nfp = 10;
  TProfile *pdt[nfp], *pmc[nfp];
  TH1D *hr[nfp];
  for (int ifp = 0; ifp != nfp; ++ifp) {
    TProfile *pd = (TProfile*)fdt->Get(Form("pmm13b%s_fp%d",cpt,ifp));
    assert(pd);
    TProfile *pm = (TProfile*)fmc->Get(Form("pmm13b%s_fp%d",cpt,ifp));
    assert(pm);
    pdt[ifp] = pd;
    pmc[ifp] = pm;
    TH1D *h = pd->ProjectionX(Form("hr%d",ifp));
    h->Divide(pm);
    for (int i = 1; i != h->GetNbinsX()+1; ++i) {
      h->SetBinContent(i, 100.*(h->GetBinContent(i)-1));
      h->SetBinError(i, 100.*h->GetBinError(i));
    } // for i
    hr[ifp] = h;
  } // for ifp
  TH1D *hsys = (TH1D*)hr[2]->Clone("hrsys");
  TH1D *hse = (TH1D*)hr[2]->Clone("hrse");
  for (int i = 1; i != hsys->GetNbinsX()+1; ++i) {
    hsys->SetBinContent(i, 0);
    hsys->SetBinError(i, fabs(hr[9]->GetBinContent(i)-hr[2]->GetBinContent(i))
		      * (1.0-0.2)/(0.9-0.2));
    hse->SetBinContent(i, fabs(hr[9]->GetBinContent(i)-hr[2]->GetBinContent(i))
		       * (1.0-0.2)/(0.9-0.2));
    hse->SetBinError(i, hr[9]->GetBinError(i));
  } // for i

  int color[nfp] = {kBlue, kCyan+1, kBlack, kGreen+2, kYellow+1,
		    kYellow+2, kOrange+1, kOrange+3, kRed, kRed+2};

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  TLatex *tex = new TLatex();
  tex->SetTextSize(0.045); tex->SetNDC();

  TH1D *hup = tdrHist("hupfp","#LTm_{W,0}#GT (GeV)",78.1+1e-5,83.1-1e-5,
		      "p_{T,jets} (GeV)",30,200);
  TH1D *hdw = tdrHist("hdwfp","Data/MC-1 (%)",-0.9,+0.9,//-0.8,+0.4,
		      "p_{T,jets} (GeV)",30,200);
  if (spt=="ptboth") hdw->SetXTitle("p_{T,jets} (GeV)");
  if (spt=="ptave")  hdw->SetXTitle("p_{T,ave} (GeV)");

  TCanvas *c1 = tdrDiCanvas("c1fp",hup,hdw,4,11);
  
  c1->cd(1);
  gPad->SetLogx();

  l->DrawLine(30,80.4,200,80.4);
  tex->DrawLatex(0.40,0.84,"|#eta_{jets}| < 1.3");
  //tex->DrawLatex(0.40,0.84,"|#eta_{jets}| < 2.5");
  tex->DrawLatex(0.40,0.78,"60<m_{W}<100 GeV");
  tex->DrawLatex(0.40,0.72,"m_{j}=0#otimesR_{q}#timesR_{#nu}");

  TLegend *leg1 = tdrLeg(0.67,0.59,0.87,0.89);
  leg1->SetHeader("DT");
  leg1->SetTextSize(0.035);

  TLegend *leg2 = tdrLeg(0.73,0.59,0.93,0.89);
  leg2->SetHeader("MC  fitProb");
  leg2->SetTextSize(0.035);

  for (int ifp = 0; ifp != nfp; ++ifp) {
    tdrDraw(pdt[ifp],"Pz",kFullCircle,color[ifp]);
    tdrDraw(pmc[ifp],"Pz",kOpenSquare,color[ifp]);

    leg1->AddEntry(pdt[ifp]," ","PLE");
    leg2->AddEntry(pmc[ifp],Form("  > %1.1f",0.1*ifp),"PLE");
  } // for ifp
  gPad->RedrawAxis();

  c1->cd(2);
  gPad->SetLogx();

  l->DrawLine(30,0,200,0);
  
  tdrDraw(hsys,"E3",kNone,kYellow,kSolid,kYellow+1,1001,kYellow);
  tdrDraw(hse,"Pz",kFullDiamond,kYellow+3);
  TF1 *fse = new TF1("fse","[0]+[1]/(log(x/0.218) * x)",30,200);
  hse->Fit(fse,"RN");
  // Is limit at infinity consistent with 0 as expected? If yes, refit
  if (fabs(fse->GetParameter(0))<2.*fse->GetParError(0)) {
    fse->FixParameter(0,0.);
    hse->Fit(fse,"RN");
  }
  fse->SetLineColor(kYellow+3);
  fse->Draw("SAME");

  TF1 *fse2 = new TF1("fse2","[0]+fabs([1])/(log(x/0.218))",40,200);
  hse->Fit(fse2,"RN");
  fse2->SetLineColor(kOrange+2);
  fse2->Draw("SAME");

  cout << Form("  // Fit from minitools/hadW.C::DrawFP(\"%s\")",cpt)<< endl;
  cout << Form("  TF1 *fhadw_%s = new TF1(\"fhadw_%s\","
	       "\"[0]+[1]/(log(x/0.218) * x)\","
	       "30,200);",cpt,cpt) << endl;
  cout << Form("  fhadw_%s->SetParameters(%1.5g, %1.5g);",
	       cpt,fse->GetParameter(0), fse->GetParameter(1)) << endl;
  cout << Form("  TF1 *fhadw2_%s = new TF1(\"fhadw2_%s\","
	       "\"[0]+fabs([1])/(log(x/0.218))\","
	       "30,200);",cpt,cpt) << endl;
  cout << Form("  fhadw2_%s->SetParameters(%1.5g, %1.5g);",
	       cpt,fse2->GetParameter(0), fse2->GetParameter(1)) << endl;

  for (int ifp = 0; ifp != nfp; ++ifp) {
    tdrDraw(hr[ifp],"Pz",kFullCircle,color[ifp]);
  } // for ifp
  gPad->RedrawAxis();

  c1->SaveAs(Form("pdf/hadW_drawFP_%s.pdf",cpt));

} // drawFP

// Analyze and draw 2D response
void hadW::Draw2D(string set) {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  //TFile *f = new TFile("rootfiles/hadWMC17_v2.root","READ");
  //TFile *f = new TFile(Form("rootfiles/hadW%s_v2.root",set.c_str()),"READ");
  TFile *f = new TFile(Form("rootfiles/hadW%s.root",set.c_str()),"READ");
  assert(f && !f->IsZombie());

  //TFile *fout = new TFile("rootfiles/hadW.root","UPDATE");
  TFile *fout = new TFile("rootfiles/hadWUL18.root","UPDATE");
  assert(fout && !fout->IsZombie());
  curdir->cd();

  TH2D *h2 = (TH2D*)f->Get("h2"); assert(h2);
  TProfile2D *p2 = (TProfile2D*)f->Get("p2"); assert(p2);
  TProfile2D *p2m2 = (TProfile2D*)f->Get("p2m2"); assert(p2m2);

  // Later take hxbyb directly from file
   const double xb[] = {30,40,50,62,80,110,150,200};
   const int nxb = sizeof(xb)/sizeof(xb[0])-1;
   const double yb[] = {0,0.261,0.522,0.783, 1.044,1.305,1.479,
			1.653,1.93,2.172, 2.322,2.5};
   const double nyb = sizeof(yb)/sizeof(yb[0])-1;
   TH2D *h2k = new TH2D("h2k",";p_{T} (GeV);#eta",nxb,xb,nyb,yb);

  // Normalize counts to forward jet
  TH2D *h2n = (TH2D*)h2->Clone("h2n");
  for (int i = 1; i != h2->GetNbinsX()+1; ++i) {
    double ni = h2->Integral(i,i,1,h2->GetNbinsY());
    for (int j = 1; j != h2->GetNbinsY()+1; ++j) {
      if (ni!=0) {
	h2n->SetBinContent(i, j, h2->GetBinContent(i,j)/ni);
	h2n->SetBinError(i, j, h2->GetBinError(i,j)/ni);
      }
    } // for j
  } // for i

  // Multiply by mass squared
  TH2D *h2x = (TH2D*)h2n->Clone("h2x");
  for (int i = 1; i != h2->GetNbinsX()+1; ++i) {
    for (int j = 1; j != h2->GetNbinsY()+1; ++j) {
      if (h2->GetBinContent(i,j)!=0) {
	h2x->SetBinContent(i, j, h2n->GetBinContent(i,j)
			   * p2m2->GetBinContent(i,j));
	h2x->SetBinError(i, j, sqrt(pow(h2n->GetBinError(i,j)/
					h2->GetBinContent(i,j),2)+
				    pow(p2m2->GetBinError(i,j)/
					p2m2->GetBinContent(i,j),2)) *
			 h2x->GetBinContent(i,j));
      }
    } // for j
  } // for i

  TCanvas *c1 = new TCanvas("c1_h2","c1_h2",600,600);
  gPad->SetRightMargin(0.15);
  h2->Draw("COLZ");
  TCanvas *c2 = new TCanvas("c2_h2n","c2_h2n",600,600);
  gPad->SetRightMargin(0.15);
  h2n->Draw("COLZ");
  TCanvas *c3 = new TCanvas("c3_p2","c3_p2",600,600);
  gPad->SetRightMargin(0.15);
  p2->Draw("COLZ");
  p2->GetZaxis()->SetRangeUser(0.85,1.15);
  TCanvas *c4 = new TCanvas("c4_p2m2","c4_p2m2",600,600);
  gPad->SetRightMargin(0.15);
  p2m2->Draw("COLZ");
  p2m2->GetZaxis()->SetRangeUser(0.85,1.15);
  TCanvas *c5 = new TCanvas("c5_h2x","c5_h2x",600,600);
  gPad->SetRightMargin(0.15);
  h2x->Draw("COLZ");

  // Invariant mass:
  // https://en.wikipedia.org/wiki/Invariant_mass
  // M^2 = 2pT1*pT2*(cosh(eta1-eta2)-cos(phi1-phi2))

  // All pairs sum up to mW => N variables, 1 equation => not enough
  // mW = (sum_ij N_ij * m_ij * sqrt(k_i * k_j)) / (sum_kl N_kl)
  // All pairs for any bin also sum up to mW => N vars, N eqs => ok?
  // mW(i) = (sum_j N_ij * m_ij * sqrt(k_i * k_j)) / (sum_k N_ik)
  // => mW(i) = sqrt(k_i) * (sum_j f_ij * m_ij * sqrt(k_j)) 
  // Tricky... maybe start with mW(i)^2 to simplify?
  // mW(i)^2 = k_i * (sum_j f_ij * m_ij^2 * k_j)
  // => mW(i)^2 - k_i*( sum_j!=i f_ij*m_ij^2*k_j) - k_i^2*f_ii*m_ii^2 = 0
  // => k_i^2*(f_ii*m_ii^2) + k_i*( sum_j!=i f_ij*m_ij^2*k_j) - mW(i)^2 = 0
  // This is quadratic equation for k_i, given k_j. Could solve it iteratively?
  // k_i = (-b +/- sqrt(b^2 - 4*a*c)) / (2*a)
  // Or just set mW(i)^2 = 1 and thus have
  // I[n,n] = X[n,n] * A[n,n] * X[n,n], where X[n,n] = I[n,n] * x[n,1]
  // Hmm, original quadratic equation of N variables can be decomposed into
  // a product of N polynomials of first degree, i.e. linear eqs. 

  // Solve for k_i iteratively
  vector<double> vk(h2x->GetNbinsX(),1.); // initial guess 1 for each
  double sumk2(h2x->GetNbinsX());
  double eps = 1e-6;
  int it(0);
  const int maxit(100);
  while (sqrt(sumk2) > eps * h2x->GetNbinsX() && it<maxit) {

    vector<double> vkref = vk;
    sumk2 = 0;
    for (int i = 1; i != h2x->GetNbinsX()+1; ++i) {
      //double fii = h2n->GetBinContent(i,i);
      //double mii = p2->GetBinContent(i,i);
      //double a = fii*mii;
      double a(h2x->GetBinContent(i,i));
      double b(0);
      double c(-1);
      for (int j = 1; j != h2x->GetNbinsY()+1; ++j) {
	if (j!=i) {
	  b += h2x->GetBinContent(i,j) * vkref[j-1];
	}
      } // for j
      // quadratic equation, positive solution
      if (a==0) { vk[i-1] = (b!=0 ? 1./b : 1.); }
      else { vk[i-1] = (-b + sqrt(b*b - 4*a*c))/(2*a); }
      
      //sumk2 += vk[i-1]*vk[i-1];
      sumk2 += pow(a*vkref[i-1]*vkref[i-1] + b*vkref[i-1] - 1,2)
	+ pow(vk[i-1]-vkref[i-1],2);
    } // for i
    cout << Form("it=%d, eps=%1.3g",it,sqrt(sumk2)/h2x->GetNbinsX())
	 << endl << flush;
    ++it;
  } // while sumk2 > eps

  // Plot k_i results in 2D pt-eta plane
  // Take hxbyb from file for binning, clone to h2k
  //vector<double> xb(h2->GetNbinsX()+1);
  //vector<double> yb(h2->GetNbinsY()+1);
  //for (unsigned int i = 0; i != xb.size(); ++i) {
  //xb[i] = h2->GetXaxis()->GetBinLowEdge(i+1);
  //}
  //for (unsigned int i = 0; i != yb.size(); ++i) {
  //yb[i] = h2->GetYaxis()->GetBinLowEdge(i+1);
  //}
  //TH2D *h2k = new TH2D("h2k",";p_{T} (GeV);#eta",
  //xb.size()-1,xb,yb.size()-1,yb);
  for (unsigned int i = 0; i != vk.size(); ++i) {
    int ix = i % (h2k->GetNbinsX()+1) + 1;
    int iy = i / (h2k->GetNbinsX()+1) + 1;
    h2k->SetBinContent(ix, iy, vk[i]);
  }

  TCanvas *c6 = new TCanvas("c6_h2k","c6_h2k",600,600);
  gPad->SetRightMargin(0.15);
  h2k->Draw("COLZ");
  h2k->GetZaxis()->SetRangeUser(0.95,1.05);
  
  c6->SaveAs(Form("pdf/hadW_draw2D_2D_%s.pdf",set.c_str()));

  // Save results to output file
  fout->cd();
  h2k->Write(Form("h2k_%s",set.c_str()),TObject::kOverwrite);
  fout->Write();
  curdir->cd();

  // If data and MC both available, also plot and store their ratio
  TH2D *h2kd = (TH2D*)fout->Get("h2k_UL17");
  TH2D *h2km = (TH2D*)fout->Get("h2k_MC17");
  if (h2kd && h2km) {
    TH2D *h2kr = (TH2D*)h2kd->Clone("h2k_ratio");
    h2kr->Divide(h2km);

    TCanvas *c7 = new TCanvas("c7_h2k","c7_h2kr",600,600);
    gPad->SetRightMargin(0.15);
    h2kr->Draw("COLZ");
    h2kr->GetZaxis()->SetRangeUser(0.98,1.02);
    
    c7->SaveAs(Form("pdf/hadW_draw2D_2D_%s.pdf","ratio"));

    fout->cd();
    h2kr->Write("h2k_ratio",TObject::kOverwrite);
    curdir->cd();
  } // h2kd && h2km

  fout->Close();
} // Draw2D()

/*
void hadW::drawVeto() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  TFile *fdt = new TFile("rootfiles/hadWUL17.root","READ");
  assert(fdt && !fdt->IsZombie());
  TFile *fmc = new TFile("rootfiles/hadWMC17.root","READ");
  assert(fdt && !fdt->IsZombie());
  curdir->cd();

  TH1D *hmw = (TH1D*)f->Get("hmw"); assert(hmw);
  TH1D *h1b = (TH1D*)f->Get("hmwveto1b"); assert(h1b);
  TH1D *h2b = (TH1D*)f->Get("hmwveto3b"); assert(h2b);
  TH1D *h1e = (TH1D*)f->Get("hmwveto1e"); assert(h1e);
  TH1D *h2e = (TH1D*)f->Get("hmwveto2r"); assert(h2e);

  

} // drawVeto
*/
