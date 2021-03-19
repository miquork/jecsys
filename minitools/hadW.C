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

#include "parsePileUpJSON.C"

#include "../tdrstyle_mod15.C"

#include <vector>
#include <algorithm> // for std::min

//const double fitProbRef = 0.2;// Z+b studies //0.01;//0.2;
const double fitProbRef = 0.01; // JEC L2Res studies
//const double etaMax = 2.5; // Z+b studies
const double etaMax = 1.3; // JEC studies
const double minLepPt = 35;//34;//30; // electron min pT 34, muon 26/29 GeV
const bool redoJEC = false;
const bool correctJetMassData = false; // data/MC ratio for mass
const bool usePtWeight = true; // reweigh top pT in MC

using namespace std;

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
   //TFile *fout = new TFile(Form("rootfiles/hadW%s_MPDGcorrNoW.root",_s.c_str()),
   //TFile *fout = new TFile(Form("rootfiles/hadW%s_EMUF.root",_s.c_str()),
   TFile *fout = new TFile(Form("rootfiles/hadW%s_Glu.root",_s.c_str()),
			   "RECREATE");
   fout->mkdir("jet");
   fout->cd("jet");
   TDirectory *doutjet = gDirectory;
   curdir->cd();

   //parsePileUpJSON("rootfiles/pileup_ASCII_2018.txt");
   parsePileUpJSON("rootfiles/pileup_ASCII_2017-2018.txt");

   bool isMC = (TString(_s.c_str()).Contains("MC"));
   bool isDT = (TString(_s.c_str()).Contains("UL"));
   assert(isMC||isDT);
   bool is18 = (TString(_s.c_str()).Contains("18"));
   bool is17 = (TString(_s.c_str()).Contains("17"));
   bool is16 = (TString(_s.c_str()).Contains("16"));
   //assert((is18 || is17) && (!is18 || !is17)); // XOR
   assert(is18 || is17 || is16);
   if ( isMC) cout << "Looping over MC" << endl;
   if ( isDT) cout << "Looping over data" << endl;

   // Fit of reco mass/pT ratio vs gen from rootfiles/drawWjetMass.C
   TF1 *fmj = new TF1("fmj","1+[1]*log(x/[0])",30,230);
   fmj->SetParameters(62.2, 0.09327); // ptave = 66.5 GeV

   // Fit of pTgen/pTave ratio from rootfiles/drawWjetMass.C
   TF1 *fpt = new TF1("fpt","[0]+[1]*exp(-[2]*x)",30,230);
   fpt->SetParameters(0.958, 6.664, 0.1424); // chi2/NDF=276.3/8

   // Fit of mjet(data)/mjet(mc)-1 (%) from rootfiles/drawWjetMass.C
   TF1 *fmd = new TF1("fmd","[0]",30,230);
   fmd->SetParameter(0, 3.294); // chi2/NDF=13.1/10

   const int nbins = 200;
   const double xmin = 30;
   const double xmax = 200;
   const int nx1 = (xmax-xmin)/5.;
   const int nx2 = (xmax-xmin)/10.;
   const double xb[] = {30,35,40,45,50,60,70,85,105,130,175,230};
   const int nxb = sizeof(xb)/sizeof(xb[0])-1;
   // For QGL studies to match dijet
   const double xd[] = {28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245};
     //{15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395};
   const int nxd = sizeof(xd)/sizeof(xd[0])-1;

   // For 2D calibration in pT-eta (from deriveL2ResNarrow.C)
   const double yb[] = {0,0.261,0.522,0.783, 1.044,1.305,1.479,
			1.653,1.93,2.172, 2.322,2.5,5.2};//,2.65};
   //2.853,2.964,3.139, 3.489,3.839,5.191};
   const double nyb = sizeof(yb)/sizeof(yb[0])-1;

   // For quick mapping to indeces
   TH1D *hxb = new TH1D("hxb",";p_{T} (GeV)",nxb,xb);
   TH1D *hyb = new TH1D("hyb",";#eta",nyb,yb);
   TH2D *hxbyb = new TH2D("hxbyb",";p_{PT} (GeV);#eta",nxb,xb,nyb,yb);

   double mwmin = 55;
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
   TProfile *pm2m13bptave = new TProfile("pm2m13bptave",";p_{T,ave};"
					 "#LTm_{W,0}^{2}#GT (GeV)",nxb,xb);
   // Flavor composition for W>qq'
   TProfile *pfuds = new TProfile("pfuds",";p_{T,ave};Frac. of uds",nxb,xb);
   TProfile *pfud = new TProfile("pfud",";p_{T,ave};Frac. of ud",nxb,xb);
   TProfile *pfu = new TProfile("pfu",";p_{T,ave};Frac. of u",nxb,xb);
   TProfile *pfd = new TProfile("pfd",";p_{T,ave};Frac. of d",nxb,xb);
   TProfile *pfs = new TProfile("pfs",";p_{T,ave};Frac. of s",nxb,xb);
   TProfile *pfc = new TProfile("pfc",";p_{T,ave};Frac. of c",nxb,xb);
   TProfile *pfb = new TProfile("pfb",";p_{T,ave};Frac. of b",nxb,xb);
   TProfile *pfg = new TProfile("pfg",";p_{T,ave};Frac. of g",nxb,xb);
   TProfile *pfo = new TProfile("pfo",";p_{T,ave};Frac. of o",nxb,xb);
   //  
   TProfile *pm1 = new TProfile("pm1",";p_{T,ave};#LTm_{W,0}^{1}#GT (GeV)",nxb,xb);
   TProfile *pm2 = new TProfile("pm2",";p_{T,ave};#LTm_{W,0}^{2}#GT (GeV)",nxb,xb);
   TProfile *pw1 = new TProfile("pw1",";p_{T,ave};#LTm_{W,gen0}^{1}#GT (GeV)",nxb,xb);
   TProfile *pw2 = new TProfile("pw2",";p_{T,ave};#LTm_{W,gen0}^{2}#GT (GeV)",nxb,xb);
   TProfile *pww1 = new TProfile("pww1",";p_{T,ave};#LTm_{W,gen0}^{1}#GT (GeV)",nxb,xb);
   TProfile *pww2 = new TProfile("pww2",";p_{T,ave};#LTm_{W,gen0}^{2}#GT (GeV)",nxb,xb);
   TProfile *p1 = new TProfile("p1",";p_{T,gen};(p_{T,reco}/p_{T,gen})^{1}",nxb,xb);
   TProfile *p2 = new TProfile("p2",";p_{T,gen};(p_{T,reco}/p_{T,gen})^{2}",nxb,xb);
   TProfile *pj1 = new TProfile("pj1",";p_{T,ave};(p_{T,reco}/p_{T,gen})^{1}",nxb,xb);
   TProfile *pj2 = new TProfile("pj2",";p_{T,ave};(p_{T,reco}/p_{T,gen})^{2}",nxb,xb);
   TProfile *pd1 = new TProfile("pd1",";p_{T,ave};((p_{T,1}-p_{T,2})/p_{T,ave})^{1}",nxb,xb);
   TProfile *pd2 = new TProfile("pd2",";p_{T,ave};((p_{T,1}-p_{T,2})/p_{T,ave})^{2}",nxb,xb);
   TProfile *pg1 = new TProfile("pg1",";p_{T,ave};((p_{T,gen1}-p_{T,gen2})/p_{T,ave})^{1}",nxb,xb);
   TProfile *pg2 = new TProfile("pg2",";p_{T,ave};((p_{T,gen1}-p_{T,gen2})/p_{T,ave})^{2}",nxb,xb);
   //
   TProfile *pmrg1 = new TProfile("pmrg1",";p_{T,gen};(m_{reco}/m_{gen})^{1}",nxb,xb);
   TProfile *pmrg2 = new TProfile("pmrg2",";p_{T,gen};(m_{reco}/m_{gen})^{2}",nxb,xb);
   TProfile *pmra1 = new TProfile("pmra1",";p_{T,ave};(m_{reco}/m_{gen})^{1}",nxb,xb);
   TProfile *pmra2 = new TProfile("pmra2",";p_{T,ave};(m_{reco}/m_{gen})^{2}",nxb,xb);
   //
   TProfile *pjm1 = new TProfile("pjm1",";p_{T,ave};(m_{reco}/p_{T,reco}^{CF})^{1}",nxb,xb);
   TProfile *pjm2 = new TProfile("pjm2",";p_{T,ave};(m_{reco}/p_{T,reco}^{CF})^{1}",nxb,xb);
   TProfile *pgm1 = new TProfile("pgm1",";p_{T,ave};(m_{gen}/p_{T,gen}^{CF})^{1}",nxb,xb);
   TProfile *pgm2 = new TProfile("pgm2",";p_{T,ave};(m_{gen}/p_{T,gen}^{CF})^{1}",nxb,xb);
   //
   TProfile *pmjet = new TProfile("pmjet",";p_{T,ave};m_{jet}/p_{T,ave}",nxb,xb);
   TProfile *pavecf = new TProfile("pavecf",";p_{T,ave};p_{T,ave}^{cf}",nxb,xb);
   TProfile *pmreco = new TProfile("pmreco",";p_{T,ave};m_{reco}/p_{T,ave}",nxb,xb);
   TProfile *pmgen = new TProfile("pmgen",";p_{T,ave};m_{gen}/p_{T,ave}",nxb,xb);
   TProfile *pavegen = new TProfile("pavegen",";p_{T,ave};p_{T,ave,gen}/p_{T,ave}",nxb,xb);
   TProfile *pavereco = new TProfile("pavereco",";p_{T,ave};p_{T,ave}",nxb,xb);
   TProfile *pavegencf = new TProfile("pavegencf",";p_{T,ave};p_{T,ave,gen}^{CF}/p_{T,ave}",nxb,xb);
   TProfile *paverecocf = new TProfile("paverecocf",";p_{T,ave};p_{T,ave,reco}^{CF}/p_{T,ave}",nxb,xb);
   //
   TH1D *hrbq = new TH1D("hrbq",";R_{bq}",500,0,5);
   TH2D *h2bq = new TH2D("h2bq",";p_{T,q,ave};p_{T,b,ave}",nxb,xb,nxb,xb);
   TH2D *h2rbqq = new TH2D("h2rbqq",";p_{T,q,ave};R_{bq}",nxb,xb,250,0,5);
   TProfile *prbqq = new TProfile("prbqq",";p_{T,q,ave};R_{bq}",nxb,xb);
   TH2D *h2rbqb = new TH2D("h2rbqb",";p_{T,b,ave};1./R_{bq}",nxb,xb,250,0,5);
   TProfile *prbqb = new TProfile("prbqb",";p_{T,b,ave};1./R_{bq}",nxb,xb);
   TH2D *h2rbqa = new TH2D("h2rbqa",";p_{T,b+q,ave};1+#DeltaR_{bq}",nxb,xb,250,-2,3);
   TProfile *prbqa = new TProfile("prbqa",";p_{T,b+q,ave};1+#DeltaR_{bq}",nxb,xb);
   TH2D *h2rbqc = new TH2D("h2rbqc",";p_{T,b+1.2*q,ave};1+#DeltaR_{bqc}",nxb,xb,250,-2,3);
   TProfile *prbqc = new TProfile("prbqc",";p_{T,b+1.2*q,ave};1+#DeltaR_{bqc}",nxb,xb);
   TProfile *pb1 = new TProfile("pb1",";p_{T,b,gen};R_{b}^{1}",nxb,xb);
   TProfile *pb2 = new TProfile("pb2",";p_{T,b,gen};R_{b}^{2}",nxb,xb);
   TProfile *pbq1 = new TProfile("pbq1",";p_{T,q,ave};R_{b}^{1}",nxb,xb);
   TProfile *pbq2 = new TProfile("pbq2",";p_{T,q,ave};R_{b}^{2}",nxb,xb);
   TProfile *pbb1 = new TProfile("pbb1",";p_{T,b,ave};R_{b}^{1}",nxb,xb);
   TProfile *pbb2 = new TProfile("pbb2",";p_{T,b,ave};R_{b}^{2}",nxb,xb);
   TProfile *pba1 = new TProfile("pba1",";p_{T,b+q,ave};R_{b}^{1}",nxb,xb);
   TProfile *pba2 = new TProfile("pba2",";p_{T,b+q,ave};R_{b}^{2}",nxb,xb);
   //
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
   //
   TProfile *pmm13bptboth_fp01 = new TProfile("pmm13bptboth_fp01",";p_{T,both};"
					      "#LTm_{W,0}#GT (GeV) (fp>0.01)",
					      nxb,xb);
   TProfile *pmm13bptboth_fp02 = new TProfile("pmm13bptboth_fp02",";p_{T,both};"
					      "#LTm_{W,0}#GT (GeV) (fp>0.02)",
					      nxb,xb);
   TProfile *pmm13bptboth_fp05 = new TProfile("pmm13bptboth_fp05",";p_{T,both};"
					      "#LTm_{W,0}#GT (GeV) (fp>0.05)",
					      nxb,xb);
   //
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
   //
   TProfile *pmm13bptave_fp01 = new TProfile("pmm13bptave_fp01",";p_{T,ave};"
					     "#LTm_{W,0}#GT (GeV) (fp>0.01)",
					     nxb,xb);
   TProfile *pmm13bptave_fp02 = new TProfile("pmm13bptave_fp02",";p_{T,ave};"
					     "#LTm_{W,0}#GT (GeV) (fp>0.02)",
					     nxb,xb);
   TProfile *pmm13bptave_fp05 = new TProfile("pmm13bptave_fp05",";p_{T,ave};"
					     "#LTm_{W,0}#GT (GeV) (fp>0.05)",
					     nxb,xb);
   //
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

   // TOP
   TH1D *h1lep = new TH1D("h1ptlep",";p_{T,lep} (GeV);",nxb,xb);
   TH1D *h1pt1 = new TH1D("h1pttt1",";p_{T,t1} (GeV);",nxb,xb);
   TH1D *h1mt = new TH1D("h1mt",";m_{t1,2} (GeV);",95,125,220);
   TH1D *h1mt1 = new TH1D("h1mt1",";m_{t1} (Gev);",95,125,220);
   TH1D *h1mt2 = new TH1D("h1mt2",";m_{t2} (GeV)",95,125,220);
   TH1D *h1mtmin = new TH1D("h1mtmin",";m_{t,max} (GeV)",95,125,220);
   TH1D *h1mtmax = new TH1D("h1mtmax",";m_{t,min} (GeV)",95,125,220);
   TH1D *h1rb = new TH1D("h1rb",";p_{t,gen-b}/p_{t,reco-b};",200,0.0,2.0);
   TH1D *h1gb = new TH1D("h1gb",";p_{t,reco-b}/p_{t,gen-b};",200,0.0,2.0);
   //
   TH2D *h2mt = new TH2D("h2mt",";p_{T,b1};m_{t1,2}",nxb,xb,95,125,220);
   TH2D *h2mt1 = new TH2D("h2mt1",";p_{T,b1};m_{t1}",nxb,xb,95,125,220);
   TH2D *h2mt2 = new TH2D("h2mt2",";p_{T,b1};m_{t2}",nxb,xb,95,125,220);
   TH2D *h2mtmin = new TH2D("h2mtmin",";p_{T,b1};m_{t,min}",nxb,xb,95,125,220);
   TH2D *h2mtmax = new TH2D("h2mtmax",";p_{T,b1};m_{t,max}",nxb,xb,95,125,220);
   TH2D *h2rb = new TH2D("h2rb",";p_{t,reco-b} (GeV);p_{t,gen-b}/p_{t,reco-b}",
			 nxb,xb,200,0.0,2.0);
   TH2D *h2gb = new TH2D("h2gb",";p_{T,gen-b} (GeV);p_{t,reco-b}/p_{t,gen-b};",
			 nxb,xb,200,0.0,2.0);
   //
   TProfile *pmt = new TProfile("pmt",";p_{T,b};m_{t}",nxb,xb);
   TProfile *pmtup = new TProfile("pmtup",";p_{T,b,up};m_{t,up}",nxb,xb);
   TProfile *pmtdw = new TProfile("pmtdw",";p_{T,b,dw};m_{t,dw}",nxb,xb);
   TProfile *pmta = new TProfile("pmta",";p_{T,bq,ave};m_{t}",nxb,xb);
   TProfile *pmt1 = new TProfile("pmt1",";p_{T,b1};m_{t1}",nxb,xb);
   TProfile *pmt2 = new TProfile("pmt2",";p_{T,b1};m_{t2}",nxb,xb);
   TProfile *pmtmin = new TProfile("pmtmin",";p_{T,b1};m_{t,min}",nxb,xb);
   TProfile *pmtmax = new TProfile("pmtmax",";p_{T,b1};m_{t,max}",nxb,xb);
   //{ {272760,276819},{276820,278801},{278802,284045} }, // BCD, EF(early), (Flate)GH
    //{ {297020,299329},{299330,302029},{302030,303434},{303435,304899},{304900,306460} }, // B,C,D,E,F: DE fused at 303434},{303435,
    //{ {306926,307082} }, // H
    //{ {315000,316999},{317000,319319},{319320,320399},{320400,326000} } // A,B, C,D
   double iovs[] = {1,
		    272760, 276820, 278802, 
		    297020, 299330, 302030, 303435, 304900,
		    306926,
		    315000, 317000, 319320, 320400,326000+1};
   const int niov = sizeof(iovs)/sizeof(iovs[0])-1;
   TProfile *pmwiov = new TProfile("pmwiov",";Run;#LTm_{W}#GT (GeV)",niov,iovs);
   TProfile *pbqiov = new TProfile("pbqiov",";Run;#LTR_{bq}#GT",niov,iovs);
   TProfile *plbiov = new TProfile("plbiov",";Run;#LTm_{lb}#GT (GeV)",niov,iovs);
   TProfile *pmtiov = new TProfile("pmtiov",";Run;#LTm_{t}#GT (GeV)",niov,iovs);
   TProfile *pmt1iov = new TProfile("pmt1iov",";Run;#LTm_{t,1}#GT (GeV)",niov,iovs);
   TProfile *ppt1iov = new TProfile("pptt1iov",";Run;#LTp_{T,t,1}#GT (GeV)",niov,iovs);
   TProfile *plepiov = new TProfile("plepiov",";Run;#LTp_{T,lep}#GT (GeV)",niov,iovs);

   // UE studies
   TProfile *pmuiov = new TProfile("pmuiov",";Run;#LT#mu#GT",niov,iovs);
   TProfile *pnpviov = new TProfile("pnpviov",";Run;#LTN_{PV}#GT",niov,iovs);
   TProfile *prhoiov = new TProfile("prhoiov",";Run;#LT#rho#GT (GeV)",niov,iovs);
   //
   TProfile *pmuvsmu = new TProfile("pmuvsmu",";#mu;#LT#mu#GT",60,0,60);
   TProfile *pnpvvsmu = new TProfile("pnpvvsmu",";#mu;#LTN_{PV}#GT",60,0,60);
   TProfile *prhovsmu = new TProfile("prhovsmu",";#mu;#LT#rho#GT (GeV)",60,0,60);

   // bJES
   TProfile *prb = new TProfile("prb",";p_{T,reco-b} (GeV);"
				"p_{T,gen-b}/p_{T,reco-b}",nxb,xb);
   TProfile *pgb = new TProfile("pgb",";p_{T,gen-b} (GeV);"
				"p_{T,reco-b}/p_{T,gen-b}",nxb,xb);
   TProfile *prb0 = new TProfile("prb0",";p_{T,reco-b,orig} (GeV);"
				 "p_{T,gen-b}/p_{T,reco-b,orig}",nxb,xb);
   TProfile *pgb0 = new TProfile("pgb0",";p_{T,gen-b} (GeV);"
				 "p_{T,reco-b,orig}/p_{T,gen-b}",nxb,xb);
   // for b-jet drawWjetMass.C::drawBjetMass mass corrections
   TProfile *pmbjet = new TProfile("pmbjet",";p_{T,b};m_{b,jet}/p_{T,b}",nxb,xb);
   TProfile *pmbreco = new TProfile("pmbreco",";p_{T,b};m_{b}/p_{T,b}",nxb,xb);
   TProfile *pmbgen = new TProfile("pmbgen",";p_{T,b};m_{b,gen}/p_{T,b}",nxb,xb);
   TProfile *pbreco = new TProfile("pbreco",";p_{T,b};p_{T,b}",nxb,xb);
   TProfile *pbgen = new TProfile("pbgen",";p_{T,b};p_{T,b,gen}/p_{T,b}",nxb,xb);

   // ATLAS-style plots
   // (https://arxiv.org/abs/1503.05472, Fig.4)
   TH1D *hatmw = new TH1D("atlas_mw",";m_{W}^{reco} (GeV)",55,55,110);
   TH1D *hatbq = new TH1D("atlas_rbq",";R_{bq}^{reco}",90,0.3,3);
   TH1D *hatmt = new TH1D("atlas_mt",";m_{top}^{reco} (GeV)",90,130,220);
   TH1D *hatlb = new TH1D("atlas_mlb",";m_{lb}^{reco} (GeV)",45,35,170);

   // Gap fractions, pT5 / ISR
   TH1D *hatlbmet = new TH1D("atlas_mlbmet",";m_{lbmet}^{reco} (GeV)",
			     285,35,320);
   //TProfile *pgap = new TProfile("pgap",";p_{T,gap} (GeV);Gap fraction",nxb,xb);
   TProfile *pgap = new TProfile("pgap",";p_{T,gap} (GeV);Gap fraction",
				 200,30,230);
   TH1D *hatpt5 = new TH1D("hatpt5",";p_{T,j5} (GeV)",nxb,xb);
   TH1D *hatpttt = new TH1D("hatpttt",";p_{T,tt} (GeV)",nxb,xb);
   TH1D *hatpt5ott = new TH1D("hatpt5ott",";p_{T,j5}/p_{T,tt}",100,0,4);
   TH2D *h2atpt5vtt = new TH2D("h2atpt5vtt",";p_{t,tt};p_{T,5}",nxb,xb,nxb,xb);
   TProfile *patpt5ott = new TProfile("patpt5ott",";p_{T,tt};p_{T,j5}/p_{T,tt}",
				      nxb,xb);

   // Check impact of event weights on mW and Rbq. Do they explain 0.3% in Rbq?
   TH1D *hatmwsw = new TH1D("atlas_mw_signal",";m_{W}^{reco} (GeV)",55,55,110);
   TH1D *hatmws1 = new TH1D("atlas_mw_signal_noweight",";m_{W}^{reco} (GeV)",
			    55,55,110);
   TH1D *hatbqsw = new TH1D("atlas_rbq_signal",";R_{bq}^{reco}",90,0.3,3);
   TH1D *hatbqs1 = new TH1D("atlas_rbq_signal_noweight",";R_{bq}^{reco}",
			    90,0.3,3);

   // Check impact of QGL>0.5 cut on mW vs pT,ave and pT,jet
   TH1D *hatmwqq = new TH1D("atlas_mwqq",";m_{W,qq}^{reco} (GeV)",55,55,110);
   TH1D *hatmwqg = new TH1D("atlas_mwqg",";m_{W,qg}^{reco} (GeV)",55,55,110);
   TH1D *hatmwgg = new TH1D("atlas_mwgg",";m_{W,gg}^{reco} (GeV)",55,55,110);
   TProfile *patmwa =   new TProfile("patmwa",  ";p_{T,ave};m_{W}",nxb,xb);
   TProfile *patmwqqa = new TProfile("patmwqqa",";p_{T,ave,qq};m_{W}",nxb,xb);
   TProfile *patmwqga = new TProfile("patmwqga",";p_{T,ave,qg};m_{W}",nxb,xb);
   TProfile *patmwgga = new TProfile("patmwgga",";p_{T,ave,gg};m_{W}",nxb,xb);
   TProfile *patmwb =  new TProfile("patmwb", ";p_{T,ave};m_{W}",nxb,xb);
   TProfile *patmwqb = new TProfile("patmwqb",";p_{T,ave,q};m_{W}",nxb,xb);
   TProfile *patmwgb = new TProfile("patmwgb",";p_{T,ave,g};m_{W}",nxb,xb);
   
   // Checks on lepton fraction for semileptonic BR and neutrinoes
   TProfile *pqelfa = new TProfile("pqelfa",";p_{T,q} (GeV);Incl. Elf",nxb,xb);
   TProfile *pqelfb = new TProfile("pqelfb",";p_{T,q} (GeV);Excl. Elf",nxb,xb);
   TProfile *pqelfn = new TProfile("pqelfn",";p_{T,q} (GeV);ElN frac.",nxb,xb);
   TProfile *pqmufa = new TProfile("pqmufa",";p_{T,q} (GeV);Incl. Muf",nxb,xb);
   TProfile *pqmufb = new TProfile("pqmufb",";p_{T,q} (GeV);Excl. Muf",nxb,xb);
   TProfile *pqmufn = new TProfile("pqmufn",";p_{T,q} (GeV);MuN frac.",nxb,xb);
   //
   TProfile *pbelfa = new TProfile("pbelfa",";p_{T,b} (GeV);Incl. Elf",nxb,xb);
   TProfile *pbelfb = new TProfile("pbelfb",";p_{T,b} (GeV);Excl. Elf",nxb,xb);
   TProfile *pbelfn = new TProfile("pbelfn",";p_{T,b} (GeV);ElN frac.",nxb,xb);
   TProfile *pbmufa = new TProfile("pbmufa",";p_{T,b} (GeV);Incl. Muf",nxb,xb);
   TProfile *pbmufb = new TProfile("pbmufb",";p_{T,b} (GeV);Excl. Muf",nxb,xb);
   TProfile *pbmufn = new TProfile("pbmufn",";p_{T,b} (GeV);MuN frac.",nxb,xb);

   // QGL shape for quarks
   TH1D *hqglq = new TH1D("hqglq",";QGL (q);N_{jets}",100,0,1);
   TH1D *hqglb = new TH1D("hqglb",";QGL (b);N_{jets}",100,0,1);
   TH2D *h2qglqa = new TH2D("h2qglqa",";QGL (q);pT_{T,q-ave};",100,0,1,nxd,xd);
   TH2D *h2qglqb = new TH2D("h2qglqb",";QGL (q);pT_{T,q};",100,0,1,nxd,xd);
   TH2D *h2qglqb_g = new TH2D("h2qglqb_g",";QGL (q);pT_{T,q};",100,0,1,nxd,xd);
   TH2D *h2qglqb_q = new TH2D("h2qglqb_q",";QGL (q);pT_{T,q};",100,0,1,nxd,xd);
   TH2D *h2qglqb_u = new TH2D("h2qglqb_u",";QGL (q);pT_{T,u};",100,0,1,nxd,xd);
   TH2D *h2qglqb_d = new TH2D("h2qglqb_d",";QGL (q);pT_{T,d};",100,0,1,nxd,xd);
   TH2D *h2qglqb_s = new TH2D("h2qglqb_s",";QGL (q);pT_{T,s};",100,0,1,nxd,xd);
   TH2D *h2qglqb_c = new TH2D("h2qglqb_c",";QGL (q);pT_{T,c};",100,0,1,nxd,xd);
   TH2D *h2qglqb_b = new TH2D("h2qglqb_b",";QGL (q);pT_{T,b};",100,0,1,nxd,xd);
   TH2D *h2qglqb_o = new TH2D("h2qglqb_o",";QGL (q);pT_{T,o};",100,0,1,nxd,xd);
   TH2D *h2qglba = new TH2D("h2qglba",";QGL (q);pT_{T,b-ave};",100,0,1,nxd,xd);
   TH2D *h2qglbb = new TH2D("h2qglbb",";QGL (q);pT_{T,b};",100,0,1,nxd,xd);
   TH2D *h2qglbb_b = new TH2D("h2qglbb_b",";QGL (q);pT_{T,b};",100,0,1,nxd,xd);
   // ...and gluons
   TH1D *hqglg = new TH1D("hqglg",";QGL (g);N_{jets}",100,0,1);
   TH2D *h2qglg = new TH2D("h2qglg",";QGL (g);pT_{T,g};",100,0,1,nxd,xd);
   TH2D *h2qglg_g = new TH2D("h2qglg_g",";QGL (g);pT_{T,g};",100,0,1,nxd,xd);
   TH2D *h2qglg_q = new TH2D("h2qglg_q",";QGL (g);pT_{T,q};",100,0,1,nxd,xd);
   TH2D *h2qglg_u = new TH2D("h2qglg_u",";QGL (g);pT_{T,u};",100,0,1,nxd,xd);
   TH2D *h2qglg_d = new TH2D("h2qglg_d",";QGL (g);pT_{T,d};",100,0,1,nxd,xd);
   TH2D *h2qglg_s = new TH2D("h2qglg_s",";QGL (g);pT_{T,s};",100,0,1,nxd,xd);
   TH2D *h2qglg_c = new TH2D("h2qglg_c",";QGL (g);pT_{T,c};",100,0,1,nxd,xd);
   TH2D *h2qglg_b = new TH2D("h2qglg_b",";QGL (g);pT_{T,b};",100,0,1,nxd,xd);
   TH2D *h2qglg_o = new TH2D("h2qglg_o",";QGL (g);pT_{T,o};",100,0,1,nxd,xd);

   // QGL efficiencies and fractions (reco flavor - truth flavor)
   TProfile *peqq = new TProfile("peqq",";p_{T,qq};Eff. QGL>0.5",nxd,xd);
   TProfile *pequ = new TProfile("pequ",";p_{T,qu};Eff. QGL>0.5",nxd,xd);
   TProfile *peqd = new TProfile("peqd",";p_{T,qd};Eff. QGL>0.5",nxd,xd);
   TProfile *peqs = new TProfile("peqs",";p_{T,qs};Eff. QGL>0.5",nxd,xd);
   TProfile *peqc = new TProfile("peqc",";p_{T,qc};Eff. QGL>0.5",nxd,xd);
   TProfile *peqb = new TProfile("peqb",";p_{T,qb};Eff. QGL>0.5",nxd,xd);
   TProfile *peqg = new TProfile("peqg",";p_{T,qg}<Eff. QGL>0.5",nxd,xd);
   TProfile *peqo = new TProfile("peqo",";p_{T,qo}<Eff. QGL>0.5",nxd,xd);

   TProfile *pebb = new TProfile("pebb",";p_{T,bb};Eff. QGL>0.5",nxd,xd);
   TProfile *pebq = new TProfile("pebq",";p_{T,bq};Eff. QGL>0.5",nxd,xd);
   TProfile *pebg = new TProfile("pebg",";p_{T,bg};Eff. QGL>0.5",nxd,xd);
   TProfile *pebo = new TProfile("pebo",";p_{T,bo};Eff. QGL>0.5",nxd,xd);

   TProfile *pegg = new TProfile("pegg",";p_{T,gg};Eff. QGL>0.5",nxd,xd);
   TProfile *pegq = new TProfile("pegq",";p_{T,gq};Eff. QGL>0.5",nxd,xd);
   TProfile *pegu = new TProfile("pegu",";p_{T,gq};Eff. QGL>0.5",nxd,xd);
   TProfile *pegd = new TProfile("pegd",";p_{T,gq};Eff. QGL>0.5",nxd,xd);
   TProfile *pegs = new TProfile("pegs",";p_{T,gq};Eff. QGL>0.5",nxd,xd);
   TProfile *pegc = new TProfile("pegc",";p_{T,gq};Eff. QGL>0.5",nxd,xd);
   TProfile *pegb = new TProfile("pegb",";p_{T,gq};Eff. QGL>0.5",nxd,xd);
   TProfile *pego = new TProfile("pego",";p_{T,go};Eff. QGL>0.5",nxd,xd);

   TProfile *pfqq = new TProfile("pfqq",";p_{T,q};Fraction of q",nxd,xd);
   TProfile *pfqu = new TProfile("pfqu",";p_{T,q};Fraction of u",nxd,xd);
   TProfile *pfqd = new TProfile("pfqd",";p_{T,q};Fraction of d",nxd,xd);
   TProfile *pfqs = new TProfile("pfqs",";p_{T,q};Fraction of s",nxd,xd);
   TProfile *pfqc = new TProfile("pfqc",";p_{T,q};Fraction of c",nxd,xd);
   TProfile *pfqb = new TProfile("pfqb",";p_{T,q};Fraction of b",nxd,xd);
   TProfile *pfqg = new TProfile("pfqg",";p_{T,q};Fraction of g",nxd,xd);
   TProfile *pfqo = new TProfile("pfqo",";p_{T,q};Fraction of o",nxd,xd);

   TProfile *pfbb = new TProfile("pfbb",";p_{T,b};Fraction of b",nxd,xd);
   TProfile *pfbq = new TProfile("pfbq",";p_{T,b};Fraction of q",nxd,xd);
   TProfile *pfbg = new TProfile("pfbg",";p_{T,b};Fraction of g",nxd,xd);
   TProfile *pfbo = new TProfile("pfbo",";p_{T,b};Fraction of o",nxd,xd);

   TProfile *pfgg = new TProfile("pfgg",";p_{T,g};Fraction of g",nxd,xd);
   TProfile *pfgq = new TProfile("pfgq",";p_{T,g};Fraction of q",nxd,xd);
   TProfile *pfgu = new TProfile("pfgu",";p_{T,g};Fraction of u",nxd,xd);
   TProfile *pfgd = new TProfile("pfgd",";p_{T,g};Fraction of d",nxd,xd);
   TProfile *pfgs = new TProfile("pfgs",";p_{T,g};Fraction of s",nxd,xd);
   TProfile *pfgc = new TProfile("pfgc",";p_{T,g};Fraction of c",nxd,xd);
   TProfile *pfgb = new TProfile("pfgb",";p_{T,g};Fraction of b",nxd,xd);
   TProfile *pfgo = new TProfile("pfgo",";p_{T,g};Fraction of o",nxd,xd);

   // 2D maps
   //const int n2b = (nxb+2)*(nyb+2);
   //const int n2b = nxb*nyb; // excluding over/underflow bins (assert)
   const int n2b = (nxb+1)*nyb; // keeping overflow bin for pt only
   TH2D *h2 = new TH2D("h2",";bin i (fwd);bin j (cnt)",n2b,0,n2b,n2b,0,n2b);
   TProfile2D *p2m1 = new TProfile2D("p2",";bin i (fwd);bin j (cnt);"
				     "(m_{W,0}/m_{W,PDG})^{2}",
				     n2b,0,n2b,n2b,0,n2b);
   TProfile2D *p2m2 = new TProfile2D("p2m2",";bin i (fwd);bin j (cnt);"
				     "(m_{W,0}/m_{W,PDG})^{2}",
				     n2b,0,n2b,n2b,0,n2b);

   // tag pairs
   TH2D *h2LHGO = new TH2D("h2LHGO",";flav1;flav2",12,-5.5,6.5,12,-5.5,6.5);
   TH2D *h2L = new TH2D("h2L",";flav1;flav2",12,-5.5,6.5,12,-5.5,6.5);
   TH2D *h2H = new TH2D("h2H",";flav1;flav2",12,-5.5,6.5,12,-5.5,6.5);
   TH2D *h2G = new TH2D("h2G",";flav1;flav2",12,-5.5,6.5,12,-5.5,6.5);
   TH2D *h2O = new TH2D("h2O",";flav1;flav2",12,-5.5,6.5,12,-5.5,6.5);
   TH1D *hmABCD = new TH1D("hmABCD",";m_{W} (GeV)",50,60,110);
   TH1D *hmA = new TH1D("hmA",";m_{W} (GeV)",50,60,110);
   TH1D *hmB = new TH1D("hmB",";m_{W} (GeV)",50,60,110);
   TH1D *hmC = new TH1D("hmC",";m_{W} (GeV)",50,60,110);
   TH1D *hmD = new TH1D("hmD",";m_{W} (GeV)",50,60,110);
   TH1D *hmLHGO = new TH1D("hmLHGO",";m_{W} (GeV)",50,60,110);
   TH1D *hmL = new TH1D("hmL",";m_{W} (GeV)",50,60,110);
   TH1D *hmH = new TH1D("hmH",";m_{W} (GeV)",50,60,110);
   TH1D *hmG = new TH1D("hmG",";m_{W} (GeV)",50,60,110);
   TH1D *hmO = new TH1D("hmO",";m_{W} (GeV)",50,60,110);
   TH1D *hnABCD = new TH1D("hnABCD",";tagLHGO",5,-0.5,4.5);
   TH1D *hnA = new TH1D("hnA",";tagLHGO",5,-0.5,4.5);
   TH1D *hnB = new TH1D("hnB",";tagLHGO",5,-0.5,4.5);
   TH1D *hnC = new TH1D("hnC",";tagLHGO",5,-0.5,4.5);
   TH1D *hnD = new TH1D("hnD",";tagLHGO",5,-0.5,4.5);
   TH1D *hnLHGO = new TH1D("hnLHGO",";tagABCD",5,-0.5,4.5);
   TH1D *hnL = new TH1D("hnL",";tagABCD",5,-0.5,4.5);
   TH1D *hnH = new TH1D("hnH",";tagABCD",5,-0.5,4.5);
   TH1D *hnG = new TH1D("hnG",";tagABCD",5,-0.5,4.5);
   TH1D *hnO = new TH1D("hnO",";tagABCD",5,-0.5,4.5);
   TH2D *h2n = new TH2D("h2n",";tagABCD;tagLHGO",5,-0.5,4.5,5,-0.5,4.5);
   TProfile2D *p2m = new TProfile2D("p2m",";tagABCD;tagLHGO",
				    5,-0.5,4.5,5,-0.5,4.5);

   TH1D *hl = new TH1D("hl",";p_{T,l} (GeV)",200,0,200);
   TH1D *hnu = new TH1D("hnu",";p_{T,nu} (GeV)",200,0,200);
   TH1D *hlnu = new TH1D("hlnu",";p_{T,l,nu} (GeV)",200,0,200);
   TH1D *hu = new TH1D("hu",";p_{T,u} (GeV)",200,0,200);
   TH1D *hd = new TH1D("hd",";p_{T,d} (GeV)",200,0,200);
   TH1D *hud = new TH1D("hud",";p_{T,u,d} (GeV)",200,0,200);
   TH1D *hc = new TH1D("hc",";p_{T,c} (GeV)",200,0,200);
   TH1D *hs = new TH1D("hs",";p_{T,s} (GeV)",200,0,200);
   TH1D *hcs = new TH1D("hcs",";p_{T,c,s} (GeV)",200,0,200);
   TProfile *plnu =new TProfile("plnu",";l, nu or l,nu;p_{T} (GeV)",3,-0.5,2.5);
   TProfile *pud =new TProfile("pud",";u, d or u,d;p_{T} (GeV)",3,-0.5,2.5);
   TProfile *pcs =new TProfile("pcs",";c, s or c,s;p_{T} (GeV)",3,-0.5,2.5);

   TProfile *gud =new TProfile("gud",";u, d or u,d;p_{T,gen} (GeV)",3,-0.5,2.5);
   TProfile *gcs =new TProfile("gcs",";c, s or c,s;p_{T,gen} (GeV)",3,-0.5,2.5);
   
   // PS weight variations, all 46-1=45
   TProfile *psmw = new TProfile("psmw",";iPS;mw",45,0.5,45.5);
   TProfile *psmt = new TProfile("psmt",";iPS;mt",45,0.5,45.5);
   TProfile *psmwb = new TProfile("psmwb",";iPS;mwb",45,0.5,45.5);
   TProfile *psmwB = new TProfile("psmwbk",";iPS;mwB",45,0.5,45.5);
   TProfile *psmlb = new TProfile("psmlb",";iPS;mlb",45,0.5,45.5);
   TProfile *psrbq = new TProfile("psrbq",";iPS;rbq",45,0.5,45.5);
   TProfile *psrbb = new TProfile("psrbb",";iPS;rbb",45,0.5,45.5);
   TProfile *psrqq = new TProfile("psrqq",";iPS;rqq",45,0.5,45.5);
   TProfile *psptb = new TProfile("psptb",";iPS;ptb",45,0.5,45.5);
   TProfile *psptB = new TProfile("psptbk",";iPS;ptB",45,0.5,45.5);
   TProfile *psptqq = new TProfile("psptqq",";iPS;ptqq",45,0.5,45.5);
   TProfile *psptw = new TProfile("psptw",";iPS;ptw",45,0.5,45.5);
   TProfile *psptt = new TProfile("psptt",";iPS;ptt",45,0.5,45.5);
   TProfile *psrb = new TProfile("psrb",";iPS;Rb",45,0.5,45.5);
   TProfile *psrq = new TProfile("psrq",";iPS;Rq",45,0.5,45.5);
   //
   TProfile *psgmw = new TProfile("psgmw",";iPS;gmw",45,0.5,45.5);
   TProfile *psgmt = new TProfile("psgmt",";iPS;gmt",45,0.5,45.5);
   TProfile *psgmwb = new TProfile("psgmwb",";iPS;gmwb",45,0.5,45.5);
   TProfile *psgmwB = new TProfile("psgmwbk",";iPS;gmwB",45,0.5,45.5);
   TProfile *psgmlb = new TProfile("psgmlb",";iPS;gmlb",45,0.5,45.5);
   TProfile *psgrbq = new TProfile("psgrbq",";iPS;grbq",45,0.5,45.5);
   TProfile *psgrbb = new TProfile("psgrbb",";iPS;grbb",45,0.5,45.5);
   TProfile *psgrqq = new TProfile("psgrqq",";iPS;grqq",45,0.5,45.5);
   TProfile *psgptb = new TProfile("psgptb",";iPS;gptb",45,0.5,45.5);
   TProfile *psgptB = new TProfile("psgptbk",";iPS;gptB",45,0.5,45.5);
   TProfile *psgptqq = new TProfile("psgptqq",";iPS;gptqq",45,0.5,45.5);
   TProfile *psgptw = new TProfile("psgptw",";iPS;gptw",45,0.5,45.5);
   TProfile *psgptt = new TProfile("psgptt",";iPS;gptt",45,0.5,45.5);

   // Decay products' radius <0.4 at pT,t>860, but plot until 1000 GeV
   TH1D *htopptorig = new TH1D("htopptorig",";top p_{T,orig} (GeV)",100,0,1000);
   TH1D *htopptgen = new TH1D("htoppttgen",";top p_{T,gen} (GeV)",100,0,1000);
   TH1D *htoppt = new TH1D("htopptt",";top p_{T} (GeV)",100,0,1000);
   TH1D *htopptw0 = new TH1D("htoppttw0",";top p_{T} (GeV)",100,0,1000);
   TH1D *htoppt_fsr = new TH1D("htoppt_fsr",";top p_{T} (GeV)",100,0,1000);
   TH1D *htoppt_q2qg = new TH1D("htoppt_q2qg",";top p_{T} (GeV)",100,0,1000);
   TH1D *htoppt_x2xg = new TH1D("htoppt_x2xg",";top p_{T} (GeV)",100,0,1000);
   TH1D *htoppt_isr = new TH1D("htoppt_isr",";top p_{T} (GeV)",100,0,1000);
   
   // MPF components: 0=MHT(nu), 1=lep, 2=b2, 3=b1, 4=q1+q2 (, 5=t1)
   TProfile *pmpf = new TProfile("pmpf",";MPF component",6,-0.5,5.5);

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
   const char *cfb = "Autumn18_V3_MC_Pythia8_b_L2Relative_AK4PFchs";
   const char *cf0 = "Autumn18_V3_MC_Pythia8_all_L2Relative_AK4PFchs";
   const char *cdref = "CondFormats/JetMETObjects/data";
   //const char *cref = "Summer19UL17_Run%s_V2M5_SimpleL1_Data_L2L3Residual_AK4PFchs";//UL17??
   string sref, sref17, sref18, sref16;
   //sref17 = "Summer19UL17_Run%s_V5_Data_L2Residual_AK4PFchs";//UL17??
   //sref17 = "Summer19UL17_Run%s_V5_Data_L2L3Residual_AK4PFchs";//UL17 v2??
   //sref17 = "Summer19UL17_RunBCDEF_V2M5_L2L3Residual_AK4PFchs";//UL17 v2?
   sref18 = "Summer19UL18_Run%s_V3_Data_L2Residual_AK4PFchs";//UL18
   sref17 = "Summer19UL17_Run%s_V2M5_SimpleL1_Data_L2L3Residual_AK4PFchs";//UL17??
   sref16 = "Summer19UL16_RunFGH_V2_DATA_L2Residual_AK4PFchs";
   if (is18) sref = sref18;
   if (is17) sref = sref17;
   if (is16) sref = sref16;
   const char *cref = sref.c_str();
   //const char *cref = "Summer19UL18_Run%s_V3_DATA_L2Residual_AK4PFchs";//UL18
   //const char *cref = "Summer19UL18_Run%s_V4_DATA_L2L3Residual_AK4PFchs";
   //const char *cnew = "Summer19UL18_Run%s_V4_DATA_L2L3Residual_AK4PFchs";//UL18
   string snew, snew17, snew18, snew16;
   //snew18 = "Summer19UL18_Run%s_V4_DATA_L2L3Residual_AK4PFchs";//UL18
   snew18 = "Summer19UL18_Run%s_V5_DATA_L2L3Residual_AK4PFchs";//UL18
   snew17 = "Summer19UL17_Run%s_V5_DATA_L2L3Residual_AK4PFchs";//UL17
   snew16 = "Summer19UL16_RunFGH_V2_DATA_L2Residual_AK4PFchs";
   if (is18) snew = snew18;
   if (is17) snew = snew17;
   if (is16) snew = snew16;
   const char *cnew = snew.c_str();
   string sfud = Form("%s/%s.txt",cd,cfud);
   cout << sfud << endl;
   string sfs = Form("%s/%s.txt",cd,cfs);
   cout << sfs << endl;
   string sfc = Form("%s/%s.txt",cd,cfc);
   cout << sfc << endl;
   string sfb = Form("%s/%s.txt",cd,cfb);
   cout << sfb << endl;
   string s0 = Form("%s/%s.txt",cd,cf0);
   cout << s0 << endl;
   //string sref = Form("%s/%s.txt",cdref,cref);
   //cout << sref << endl;
   vector<JetCorrectorParameters> vfud;
   JetCorrectorParameters *parud = new JetCorrectorParameters(sfud);
   vfud.push_back(*parud);
   FactorizedJetCorrector *jecfud = new FactorizedJetCorrector(vfud);
   vector<JetCorrectorParameters> vfs;
   JetCorrectorParameters *pars = new JetCorrectorParameters(sfs);
   vfs.push_back(*pars);
   FactorizedJetCorrector *jecfs = new FactorizedJetCorrector(vfs);
   vector<JetCorrectorParameters> vfc;
   JetCorrectorParameters *parc = new JetCorrectorParameters(sfc);
   vfc.push_back(*parc);
   FactorizedJetCorrector *jecfc = new FactorizedJetCorrector(vfc);
   vector<JetCorrectorParameters> vfb;
   JetCorrectorParameters *parb = new JetCorrectorParameters(sfb);
   vfb.push_back(*parb);
   FactorizedJetCorrector *jecfb = new FactorizedJetCorrector(vfb);
   //
   vector<JetCorrectorParameters> v0;
   JetCorrectorParameters *pfl0 = new JetCorrectorParameters(s0);
   v0.push_back(*pfl0);
   FactorizedJetCorrector *jec0 = new FactorizedJetCorrector(v0);
   //


   /*
   const int nrunmax = 5;
   const int nrun17 = 5;
   const char* runs17[nrunmax] = {"B","C","D","E","F"};
   double runlums17[nrunmax] = {4.8,9.6,4.2,9.3,13.4};
   const int nrun18 = 4;
   const char* runs18[nrunmax] = {"A","B","C","D", ""}; //UL18
   double runlums18[nrunmax] = {14.0, 7.1, 6.9, 31.9, 0}; //UL18
   int nrun(is17 ? nrun17 : (is18 ? nrun18 : 0));
   const char* *runs = (is17 ? runs17 : (is18 ? runs18 : 0));
   double *runlums = (is17 ? runlums17 : (is18 ? runlums18 : 0));
   */
   const int nrun = pmtiov->GetNbinsX();
   const char *runs[] = {"MC","BCD","EF","GH", "B","C","D","E","F", "H",
			 "A","B","C","D"};
   assert(nrun==(sizeof(runs)/sizeof(runs[0])));
   double runlums[] = {150, 1,1,1, 4.8,9.6,4.2,9.3,13.4, 1,
		       14.0,7.1,6.9,31.9}; //UL18

   vector<FactorizedJetCorrector*> jecrefs(nrun);
   vector<FactorizedJetCorrector*> jecnews(nrun);
   //for (int irun = 0; irun != nrun; ++irun) {
   for (int irun = 1; irun != nrun; ++irun) {

     if (is16 && (irun<1||irun>3)) continue;
     if (is17 && (irun<4||irun>8)) continue;
     if (is18 && (irun<10||irun>13)) continue;

     string srun = Form(cref,runs[irun]); // UL18,UL17
     const char *cref = srun.c_str(); // UL18,UL17
     string sref = Form("%s/%s.txt",cdref,cref);
     cout << sref << endl;
     vector<JetCorrectorParameters> vref;
     JetCorrectorParameters *pref = new JetCorrectorParameters(sref);
     vref.push_back(*pref);
     FactorizedJetCorrector *jecref = new FactorizedJetCorrector(vref);
     jecrefs[irun] = jecref;
     //
     string srun2 = Form(cnew,runs[irun]);
     const char *cnew = srun2.c_str();
     string snew = Form("%s/%s.txt",cdref,cnew);
     cout << snew << endl;
     vector<JetCorrectorParameters> vnew;
     JetCorrectorParameters *pnew = new JetCorrectorParameters(snew);
     vnew.push_back(*pnew);
     FactorizedJetCorrector *jecnew = new FactorizedJetCorrector(vnew);
     jecnews[irun] = jecnew;
   }

   TFile *fjv(0);
   if (is18) fjv = new TFile("../JECDatabase/jet_veto_maps/Summer19UL18_V1/hotjets-UL18.root","READ");
   if (is17) fjv = new TFile("../JECDatabase/jet_veto_maps/Summer19UL17_V2/hotjets-UL17_v2.root","READ");
   if (is16) fjv = new TFile("../JECDatabase/jet_veto_maps/Summer19UL16_V0/hotjets-UL16.root","READ");
   assert(fjv && !fjv->IsZombie());
   TH2D *h2jv(0);
   if (is18) h2jv = (TH2D*)fjv->Get("h2hot_ul18_plus_hem1516_and_hbp2m1");
   if (is17) h2jv = (TH2D*)fjv->Get("h2hot_ul17_plus_hep17_plus_hbpw89");
   if (is16) h2jv = (TH2D*)fjv->Get("h2hot_ul16_plus_hbm2_hbp12_qie11");
   assert(h2jv);
   int igood(0);
   int ibad(0);

   // To-do:
   // Separate out different flavor combinations
   // TCut cs("(flav1==4&&flav2==-3)||(flav1==3&&flav2==-4)||(flav1==-4&&flav2==3)||(flav1==-3&&flav2==4)")
   // TCut ud("(flav1==1&&flav2==-2)||(flav1==2&&flav2==-1)||(flav1==-2&&flav2==1)||(flav1==-1&&flav2==2)")
   // TCut cd("(flav1==4&&flav2==-1)||(flav1==1&&flav2==-4)||(flav1==-4&&flav2==1)||(flav1==-1&&flav2==4)")
   // TCut us("(flav1==3&&flav2==-2)||(flav1==2&&flav2==-3)||(flav1==-3&&flav2==2)||(flav1==-2&&flav2==3)")
   // TCut hit = (cs||ud||cd||us);
   // TCut miss = !hit;
   // TCut g("flav1==21||flav2==2");
   
   int cnt(0);
   map<int, map<int, map<int, int> > > dubs;
   
   // Data/MC ratio from minitools/hadW_toppt.C
   TF1 *f1pt = new TF1("f1pt","[p0]+[p1]*pow(x,[p2])",1,400);
   f1pt->SetParameters(1.249,-0.0208,0.5182);

   // Data/MC ratio from minitools/hadW_toppt.C (iteration 2)
   // Chisquare / NDF = 59.1 / 37 (1.598)
   TF1 *f1pt2 = new TF1("f1pt2","[p0]+[p1]*pow(x,[p2])",1,400);
   f1pt2->SetParameters(2.186,-1.22,-0.005904);


   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (jentry%100000==00) cout << "." << flush;

      // Clean out duplicates for 2018 data
      // MC18 also has duplicates and quadruplicates, plus unique event number
      // UL17 quadruplicates have 1 filled, 3 empty
      // and for now UL17 is missing event number branch
      if (is18 || is17 || is16) {
	if (dubs[run][lumiBlock][event]==0) {
	  dubs[run][lumiBlock][event] = 1;
	}
	else {
	  if (++cnt<10) cout << "Duplicate run "<<run<<" ls "<<lumiBlock
			     << " event "<<event<<". Skipping" << endl<< flush;
	  continue;
	}
      }
      // UL17 duplicates and quadruplicates have bpt1 empty
      //if (is17) {
      //if (bpt1==0) continue;
      //}

      // Homogenize lepton pT cuts
      if (lep_pt < minLepPt) continue;

      // Jet veto filters
      //assert(false); // Add jet veto filters
      if (true) { // jet veto
	int i1 = h2jv->GetXaxis()->FindBin(eta1);
	int j1 = h2jv->GetYaxis()->FindBin(phi1);
	int i2 = h2jv->GetXaxis()->FindBin(eta2);
	int j2 = h2jv->GetYaxis()->FindBin(phi2);
	if (h2jv->GetBinContent(i1,j1)>0 || h2jv->GetBinContent(i2,j2)>0) {
	  ++ibad;
	  continue;
	}
	else
	  ++igood;
      } // jet veto


      // Make backup copy of original pt and mass
      double recoWmass_orig = recoWMass;
      double pt1_orig = pt1;
      double pt2_orig = pt2;
      double m1_orig = m1;
      double m2_orig = m2;
      j1.SetPtEtaPhiM(pt1_orig,eta1,phi1,m1_orig);
      j2.SetPtEtaPhiM(pt2_orig,eta2,phi2,m2_orig);
      double recoWmass_origM = (j1+j2).M();      

      double bpt1_orig = bpt1;
      double bpt2_orig = bpt2;
      double bm1_orig = bm1;
      double bm2_orig = bm2;

      // Update measured masses and pt's
      pt1 = pt1_orig * fpt->Eval(pt1_orig);
      pt2 = pt2_orig * fpt->Eval(pt2_orig);
      m1 = m1_orig / fmj->Eval(pt1_orig) * fpt->Eval(pt1_orig);
      m2 = m2_orig / fmj->Eval(pt2_orig) * fpt->Eval(pt1_orig);
      if (!isMC && correctJetMassData) {
	m1 *= 1./(1+0.01*fmd->Eval(pt1_orig));
	m2 *= 1./(1+0.01*fmd->Eval(pt2_orig));
      }

      //double recoWmass = recoWMass;
      // Alternative with massless jets
      //double m0 = 0;
      // or with quark masses (mostly charm mass)
      //double mq = 0;//(0.5*0.03+0.25*0.095+0.25*1.275);
      // with reco/gen corrected mass
      //double mj1 = m1/fmj->Eval(pt1);
      //double mj2 = m2/fmj->Eval(pt2);
      // or with QCD jet mass (guesstimate sqrt(pt)) => too high, too uncertain
      //double mj1 = sqrt(pt1);
      //double mj2 = sqrt(pt2);
      //j1.SetPtEtaPhiM(pt1,eta1,phi1,m1);
      //j2.SetPtEtaPhiM(pt2,eta2,phi2,m2);
      //double recoWmass0 = (j1+j2).M();

      // Almost genWmass (keep reco jet directions): worse than reco
      //j1.SetPtEtaPhiM(gen_pt1,eta1,phi1,mq);
      //j2.SetPtEtaPhiM(gen_pt2,eta2,phi2,mq);
      // genWmass0 (McorrNoW: change mq=0 to gen_m1, gen_m2)
      j1.SetPtEtaPhiM(gen_pt1,gen_eta1,gen_phi1,gen_m1);
      j2.SetPtEtaPhiM(gen_pt2,gen_eta2,gen_phi2,gen_m2);
      double genWmass0 = (j1+j2).M();
      double genWpt0 = (j1+j2).Pt();

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

      // Recalculate JEC for b jets
      jecfb->setJetPt(bpt1);
      jecfb->setJetEta(beta1);
      double jecfb1 = jecfb->getCorrection();
      jec0->setJetPt(bpt1);
      jec0->setJetEta(beta1);
      double jecb01 = jec0->getCorrection();
      //
      jecfb->setJetPt(bpt1);
      jecfb->setJetEta(beta1);
      double jecfb2 = jecfb->getCorrection();
      jec0->setJetPt(bpt2);
      jec0->setJetEta(beta2);
      double jecb02 = jec0->getCorrection();

      // Calculate effective L2Res or L2L3Res to undo it
      // Calculate new L2L3Res to reapply it
      double jesref1(0), jesref2(0), lumsum(0);
      double jesnew1(0), jesnew2(0);
      double jesbref1(0), jesbref2(0);
      double jesbnew1(0), jesbnew2(0);
      int irunref = pmtiov->FindBin(run+0.5)-1;
      for (int irun = 0; irun != nrun; ++irun) {

	if (isMC || !redoJEC) { // No recalculation needed for MC
	  jesref1 = jesref2 = jesnew1 = jesnew2 = 1;
	  jesbref1 = jesbref2 = jesbnew1 = jesbnew2 = 1;
	  lumsum = 1;
	  continue;
	}
	if (irun!=irunref) continue; // Don't average, use IOV JEC

	double lumi = runlums[irun];
	lumsum += lumi;
	FactorizedJetCorrector *jecref = jecrefs[irun];
	jecref->setJetPt(pt1);
	jecref->setJetEta(eta1);
	jesref1 += (isMC ? 1 : 1./jecref->getCorrection()) * lumi;
	jecref->setJetPt(pt2);
	jecref->setJetEta(eta2);
	jesref2 += (isMC ? 1 : 1./jecref->getCorrection()) * lumi;
	//
	jecref->setJetPt(bpt1);
	jecref->setJetEta(beta1);
	jesbref1 += (isMC ? 1 : 1./jecref->getCorrection()) * lumi;
	jecref->setJetPt(bpt2);
	jecref->setJetEta(beta2);
	jesbref2 += (isMC ? 1 : 1./jecref->getCorrection()) * lumi;
	//
	FactorizedJetCorrector *jecnew = jecnews[irun];
	jecnew->setJetPt(pt1);
	jecnew->setJetEta(eta1);
	jesnew1 += (isMC ? 1 : 1./jecnew->getCorrection()) * lumi;
	jecnew->setJetPt(pt2);
	jecnew->setJetEta(eta2);
	jesnew2 += (isMC ? 1 : 1./jecnew->getCorrection()) * lumi;
	//
	jecnew->setJetPt(bpt1);
	jecnew->setJetEta(beta1);
	jesbnew1 += (isMC ? 1 : 1./jecnew->getCorrection()) * lumi;
	jecnew->setJetPt(bpt2);
	jecnew->setJetEta(beta2);
	jesbnew2 += (isMC ? 1 : 1./jecnew->getCorrection()) * lumi;
      }
      jesref1 /= lumsum;
      jesref2 /= lumsum;
      jesnew1 /= lumsum;
      jesnew2 /= lumsum;
      //
      jesbref1 /= lumsum;
      jesbref2 /= lumsum;
      jesbnew1 /= lumsum;
      jesbnew2 /= lumsum;

      // W mass with new L3L2Res
      double k021 = jesref1/jesnew1;
      double k022 = jesref2/jesnew2;
      j1.SetPtEtaPhiM(pt1*k021,eta1,phi1,m1*k021);
      j2.SetPtEtaPhiM(pt2*k022,eta2,phi2,m2*k022);
      double recoWmass02 = (j1+j2).M();
      
      // W mass with flavor corrections and quark (mq=0) mass
      double kQ1 = jecfl1/jec01 * k021;
      double kQ2 = jecfl2/jec02 * k022;
      j1.SetPtEtaPhiM(pt1*kQ1,eta1,phi1,m1*kQ1);
      j2.SetPtEtaPhiM(pt2*kQ2,eta2,phi2,m2*kQ2);
      double recoWmassQ2 = (j1+j2).M();

      // Also account for neutrino pT lost by charm quarks,
      // which is ~6% of charm quark pT for semileptonic decays (2*16.0%)?
      double jesnu = 1-0.25*0.32*0.06;
      double kQNU1 = 1./jesnu * kQ1;
      double kQNU2 = 1./jesnu * kQ2;
      j1.SetPtEtaPhiM(pt1*kQNU1,eta1,phi1,m1*kQNU1);
      j2.SetPtEtaPhiM(pt2*kQNU2,eta2,phi2,m2*kQNU2);
      double recoWmassQNU2 = (j1+j2).M();

      // Calculate recoWMassQNU without L2Res (or L2L3Res)
      double kUnc1 = kQNU1 * jesnew1;
      double kUnc2 = kQNU2 * jesnew2;
      j1.SetPtEtaPhiM(pt1*kUnc1,eta1,phi1,m1*kUnc2);
      j2.SetPtEtaPhiM(pt2*kUnc2,eta2,phi2,m2*kUnc2);
      //double recoWmassUnc = (j1+j2).M(); // TMP PATCH TMP

      // Calculate recoWMassQNU with new L2L3Res
      // (had a bug for j2, jesnew2->jesnew1, fixed 20201027)
      //double kQNU21 = kQNU1; 
      //double kQNU22 = kQNU2; 
      //j1.SetPtEtaPhiM(pt1*kQNU21,eta1,phi1,m1*kQNU21);
      //j2.SetPtEtaPhiM(pt2*kQNU22,eta2,phi2,m2*kQNU22);
      //double recoWmassQNU2 = (j1+j2).M();

      // Replacements
      //
      double recoWmassUnc = recoWmassQNU2; // TMP PATCH TMP V4 closure
      //
      double recoWmass = recoWmass02; // V4 closure
      double recoWmass0 = recoWmassQNU2; // fitProb studies - V4

      // First stop, W mass(es)
      j1.SetPtEtaPhiM(pt1*kQNU1,eta1,phi1,m1*kQNU1);
      j2.SetPtEtaPhiM(pt2*kQNU2,eta2,phi2,m2*kQNU2);
      TLorentzVector qq, W, b1, b2, t1, t2, Wb, WB, nu;
      TLorentzVector t1up, t1dw, t2up, t2dw;

      // W boson without (0D) and with (1D) correction to PDG mass
      qq = (j1+j2); // raw hadronic W, mass not fixed to PDG value (0D)
      W = (qq * (80.4 / qq.M())); // Fix mass to PDG (1D)

      // Next, b jets with corrections
      double kb21 = (jesbref1/jesbnew1);
      double kb22 = (jesbref2/jesbnew2);
      double kB1 = (jecfb1/jecb01) * kb21;
      double kB2 = (jecfb2/jecb02) * kb22;
      double jesbnu = 1-1.00*0.32*0.15; // 32% semileptonic with 15% nu
      double kBNU1 = (1./jesbnu) * kB1;
      double kBNU2 = (1./jesbnu) * kB2;
      b1.SetPtEtaPhiM(bpt1*kBNU1,beta1,bphi1,bm1*kBNU1);
      b2.SetPtEtaPhiM(bpt2*kBNU2,beta2,bphi2,bm2*kBNU2);

      // Leptonic mass observable m_lb
      TLorentzVector l,lb2,met,lb2met;
      l.SetPtEtaPhiM(lep_pt,lep_eta,lep_phi,lep_m);
      b2.SetPtEtaPhiM(kBNU2*bpt2,beta2,bphi2,kBNU2*bm2);
      lb2 = (l+b2);
      met.SetPxPyPzE(METx,METy,0.,sqrt(METx*METx+METy*METy));
      lb2met = (lb2+met);

      // Two alternative pairings of b and hadronic W
      // (first one is now sorted to be hadronic top candidate)
      t1 = (b1+qq); // (0D)
      t2 = (b2+qq);
      Wb = (W + b1); // (1D)
      WB = (W + b1 * (80.4 / qq.M())); // (2D)
      nu = (qq+b1+b2+l); // nu=MHT
 
      // bJES variations
      double djes = 0.005;
      b1.SetPtEtaPhiM(bpt1*kBNU1*(1+djes),beta1,bphi1,bm1*kBNU1*(1+djes));
      b2.SetPtEtaPhiM(bpt2*kBNU2*(1+djes),beta2,bphi2,bm2*kBNU2*(1+djes));
      t1up = (b1+qq);
      t2up = (b2+qq);
      b1.SetPtEtaPhiM(bpt1*kBNU1*(1-djes),beta1,bphi1,bm1*kBNU1*(1-djes));
      b2.SetPtEtaPhiM(bpt2*kBNU2*(1-djes),beta2,bphi2,bm2*kBNU2*(1-djes));
      t1dw = (b1+qq);
      t2dw = (b2+qq);

      // Bunch of gen observables
      TLorentzVector gq1, gq2, gb1, gb2;
      gq1.SetPtEtaPhiM(gen_pt1, gen_eta1, gen_phi1, gen_m1);
      gq2.SetPtEtaPhiM(gen_pt2, gen_eta2, gen_phi2, gen_m2);
      gb1.SetPtEtaPhiM(gen_bpt1, gen_beta1, gen_bphi1, gen_bm1);
      gb2.SetPtEtaPhiM(gen_bpt2, gen_beta2, gen_bphi2, gen_bm2);
      bool goodgen = (gq1.Pt()>0 && gq2.Pt()>0 &&
		      gb1.Pt()>0 && gb2.Pt()>0);
      TLorentzVector gw, gW, gt1, gWb, gWB, glb2, nil(0,0,0,0);
      gw = (gq1+gq2);
      double gmw = gw.M();
      gW = (gmw!=0 ? gw * (80.4 / gmw) : nil);
      gt1 = (gw + gb1);
      gWb = (gW + gb1);
      gWB = (gmw!=0 ? gW + gb1 * (80.4 / gmw) : nil);
      glb2 = (l + gb2);
      
      // Weight events with fitProb (as done in mt analysis)
      //double w = fitProb;
      //double w = (isMC ? (*PSWgts)[1] : 1.);
      // should use weights, but results unstable
      double w = (isMC ? weight : 1.);
      //double w = (isMC ? min((double)weight,4.0e8) : 1.);
      double w_nopt = w; // no pT weights

      if (usePtWeight && isMC) {
	//double pt = gt1.Pt();
	//w *= (pt>1 ? f1pt->Eval(pt)*f1pt2->Eval(pt) : 1.);
	w *= tptweight;
      }

      // Should this also get corrections?
      double ptave = 0.5*(pt1+pt2);

      if (fitProb>fitProbRef) {
	
	hmw->Fill(recoWmass,w);
	hmm->Fill(recoWmass0,w);
	hpt->Fill(pt1,w);
	hpt->Fill(pt2,w);
	hpt2->Fill(pt1,pt2,w);
	if (fabs(eta1)<1.3 && fabs(eta2)<1.3) hpt213->Fill(pt1,pt2,w);
	hptave->Fill(ptave,w);
	hptmin->Fill(min(pt1,pt2),w);
	
	// Check HBPw8/9 veto region
	if ((eta1>HBPw89.GetX1() && eta1<HBPw89.GetX2() &&
	     phi1>HBPw89.GetY1() && phi1<HBPw89.GetY2()) ||
	    (eta2>HBPw89.GetX1() && eta2<HBPw89.GetX2() &&
	     phi2>HBPw89.GetY1() && phi2<HBPw89.GetY2())) {
	  hmwveto1b->Fill(recoWmass,w);
	  hmmveto1b->Fill(recoWmass0,w);
	}
	if ((eta1>HBPw89.GetX1() && eta1<HBPw89.GetX2() &&
	     phi1>HBPw89.GetY1() && phi1<HBPw89.GetY2()) &&
	    (eta2>HBPw89.GetX1() && eta2<HBPw89.GetX2() &&
	     phi2>HBPw89.GetY1() && phi2<HBPw89.GetY2())) {
	  hmwveto2b->Fill(recoWmass,w);
	  hmmveto2b->Fill(recoWmass0,w);
	}
	// Check HEP17 veto region
	if ((eta1>HEP17.GetX1() && eta1<HEP17.GetX2() &&
	     phi1>HEP17.GetY1() && phi1<HEP17.GetY2()) ||
	    (eta2>HEP17.GetX1() && eta2<HEP17.GetX2() &&
	     phi2>HEP17.GetY1() && phi2<HEP17.GetY2())) {
	  hmwveto1e->Fill(recoWmass,w);
	  hmmveto1e->Fill(recoWmass0,w);
	}
	if ((eta1>HEP17.GetX1() && eta1<HEP17.GetX2() &&
	     phi1>HEP17.GetY1() && phi1<HEP17.GetY2()) &&
	    (eta2>HEP17.GetX1() && eta2<HEP17.GetX2() &&
	     phi2>HEP17.GetY1() && phi2<HEP17.GetY2())) {
	  hmwveto2e->Fill(recoWmass,w);
	  hmmveto2e->Fill(recoWmass0,w);
	}

	if (pt1>30 && pt2>30) hmw30->Fill(recoWmass,w);
	if (pt1>40 && pt2>40) hmw40->Fill(recoWmass,w);
	if (pt1>50 && pt2>50) hmw50->Fill(recoWmass,w);
	if (pt1>60 && pt2>60) hmw60->Fill(recoWmass,w);
	//
	if (pt1>30 && pt2>30) hmm30->Fill(recoWmass0,w);
	if (pt1>40 && pt2>40) hmm40->Fill(recoWmass0,w);
	if (pt1>50 && pt2>50) hmm50->Fill(recoWmass0,w);
	if (pt1>60 && pt2>60) hmm60->Fill(recoWmass0,w);

	bool goodwa = (recoWmass0>55 && recoWmass0<110);
	bool goodw = (recoWmass>65 && recoWmass<105);
	bool goodw0 = (recoWmass0>60 && recoWmass0<100);
	bool goodgw0 = (genWmass0>60 && genWmass0<100);
	// do proper deltaPhi later
	bool gooddr1 = (pow(eta1-gen_eta1,2)+pow(phi1-gen_phi1,2) < 0.2*0.2);
	bool gooddr2 = (pow(eta2-gen_eta2,2)+pow(phi2-gen_phi2,2) < 0.2*0.2);

	if (fabs(eta1)<etaMax && fabs(eta2)<etaMax) {

	  if (pt1>30 && pt2>30) hmw1330->Fill(recoWmass,w);
	  if (pt1>40 && pt2>40) hmw1340->Fill(recoWmass,w);
	  if (pt1>50 && pt2>50) hmw1350->Fill(recoWmass,w);
	  if (pt1>60 && pt2>60) hmw1360->Fill(recoWmass,w);
	  //
	  if (pt1>30 && pt2>30) hmm1330->Fill(recoWmass0,w);
	  if (pt1>40 && pt2>40) hmm1340->Fill(recoWmass0,w);
	  if (pt1>50 && pt2>50) hmm1350->Fill(recoWmass0,w);
	  if (pt1>60 && pt2>60) hmm1360->Fill(recoWmass0,w);

	  if (goodw || goodw0 || goodwa) {
	    hptjet13->Fill(pt1,w);
	    hptjet13->Fill(pt2,w);
	    hptave13->Fill(ptave,w);
	    hptave13b->Fill(ptave,w);
	    pptave13b->Fill(ptave,ptave,w);
	    hptmin13->Fill(min(pt1,pt2),w);
	    //
	    if (goodw) pmw13ptjet->Fill(pt1,recoWmass,w);
	    if (goodw) pmw13ptjet->Fill(pt2,recoWmass,w);
	    if (goodw) pmw13ptave->Fill(ptave,recoWmass,w);
	    if (goodw) pmw13bptave->Fill(ptave,recoWmass,w);
	    if (goodw) pmw13ptmin->Fill(min(pt1,pt2),recoWmass,w);
	    //
	    if (goodw0) pmm13ptjet->Fill(pt1,recoWmass0,w);
	    if (goodw0) pmm13ptjet->Fill(pt2,recoWmass0,w);
	    if (goodw0) pmm13ptave->Fill(ptave,recoWmass0,w);
	    if (goodw0) pmm13bptave->Fill(ptave,recoWmass0,w);
	    if (goodw0) pm2m13bptave->Fill(ptave,pow(recoWmass0,2),w);
	    if (goodw0) pmm13bptave_noL3Res->Fill(ptave,recoWmassUnc,w);
	    if (goodw0) pmm13ptmin->Fill(min(pt1,pt2),recoWmass0,w);

	    // Fill flavor composition compatible with pmm13bptave
	    // i.e. recoWmass0 (recoWmassQNU2) and pt1,pt2
	    if (goodw0) {
	      // Calculate pT-weighted W>qq' flavors
	      double w1 = pt1/(2*ptave);
	      int id1 = abs(flav1);
	      int id2 = abs(flav2);
	      pfuds->Fill(ptave, w1*(id1==1||id1==2||id1==3 ? 1 : 0) +
			  (1-w1)   *(id2==1||id2==2||id2==3 ? 1 : 0), w);
	      pfud->Fill(ptave, w1*(id1==1||id1==2 ? 1 : 0) +
			 (1-w1)   *(id2==1||id2==2 ? 1 : 0), w);
	      pfu->Fill(ptave, w1*(id1==2 ? 1 : 0) +
			(1-w1)   *(id2==2 ? 1 : 0), w);
	      pfd->Fill(ptave, w1*(id1==1 ? 1 : 0) +
			(1-w1)   *(id2==1 ? 1 : 0), w);
	      pfs->Fill(ptave, w1*(id1==3 ? 1 : 0) +
			(1-w1)   *(id2==3 ? 1 : 0), w);
	      pfc->Fill(ptave, w1*(id1==4 ? 1 : 0) +
			(1-w1)   *(id2==4 ? 1 : 0), w);
	      pfb->Fill(ptave, w1*(id1==5 ? 1 : 0) +
			(1-w1)   *(id2==5 ? 1 : 0), w);
	      pfg->Fill(ptave, w1*(id1==21 ? 1 : 0) +
			(1-w1)   *(id2==21 ? 1 : 0), w);
	      pfo->Fill(ptave, w1*(id1==0 ? 1 : 0) +
			(1-w1)   *(id2==0 ? 1 : 0), w);
	    }

	    double CF = 2./3.;
	    if (goodw0) pmjet->Fill(ptave,m1/ptave,w);
	    if (goodw0) pmjet->Fill(ptave,m2/ptave,w);
	    if (goodw0) pavecf->Fill(ptave,CF*pow(ptave,CF)/ptave,w);

	    if (goodw0 &&
		gen_pt1>0 && fabs(pt1/gen_pt1-1)<0.5 && gooddr1 &&
		gen_pt2>0 && fabs(pt2/gen_pt2-1)<0.5 && gooddr2) {
	      pm1->Fill(ptave,recoWmass0,w);
	      pm2->Fill(ptave,pow(recoWmass0,2),w);
	      pw1->Fill(ptave,genWmass0,w);
	      pw2->Fill(ptave,pow(genWmass0,2),w);
	      if (goodgw0) {
		pww1->Fill(ptave,genWmass0,w);
		pww2->Fill(ptave,pow(genWmass0,2),w);
	      }
	      p1->Fill(gen_pt1,pt1/gen_pt1,w);
	      p1->Fill(gen_pt2,pt2/gen_pt2,w);
	      p2->Fill(gen_pt1,pow(pt1/gen_pt1,2),w);
	      p2->Fill(gen_pt2,pow(pt2/gen_pt2,2),w);
	      pj1->Fill(ptave,pt1/gen_pt1,w);
	      pj1->Fill(ptave,pt2/gen_pt2,w);
	      pj2->Fill(ptave,pow(pt1/gen_pt1,2),w);
	      pj2->Fill(ptave,pow(pt2/gen_pt2,2),w);
	      pd1->Fill(ptave,(pt1-pt2)/ptave,w);
	      pd2->Fill(ptave,pow((pt1-pt2)/ptave,2),w);
	      pg1->Fill(ptave,(gen_pt1-gen_pt2)/ptave,w);
	      pg2->Fill(ptave,pow((gen_pt1-gen_pt2)/ptave,2),w);
	      //
	      if (gen_m1>0.01) pmrg1->Fill(gen_pt1, m1/gen_m1, w);
	      if (gen_m2>0.01) pmrg1->Fill(gen_pt2, m2/gen_m2, w);
	      if (gen_m1>0.01) pmrg2->Fill(gen_pt1, pow(m1/gen_m1,2), w);
	      if (gen_m2>0.01) pmrg2->Fill(gen_pt2, pow(m2/gen_m2,2), w);
	      if (gen_m1>0.01) pmra1->Fill(ptave, m1/gen_m1, w);
	      if (gen_m2>0.01) pmra1->Fill(ptave, m2/gen_m2, w);
	      if (gen_m1>0.01) pmra2->Fill(ptave, pow(m1/gen_m1,2), w);
	      if (gen_m2>0.01) pmra2->Fill(ptave, pow(m2/gen_m2,2), w);
	      //
	      pjm1->Fill(ptave,m1/(CF*pow(pt1,CF)),w);
	      pjm1->Fill(ptave,m2/(CF*pow(pt2,CF)),w);
	      pjm2->Fill(ptave,pow(m1/(CF*pow(pt1,CF)),2),w);
	      pjm2->Fill(ptave,pow(m2/(CF*pow(pt2,CF)),2),w);
	      pgm1->Fill(ptave,gen_m1/(CF*pow(gen_pt1,CF)),w);
	      pgm1->Fill(ptave,gen_m2/(CF*pow(gen_pt2,CF)),w);
	      pgm2->Fill(ptave,pow(gen_m1/(CF*pow(gen_pt1,CF)),2),w);
	      pgm2->Fill(ptave,pow(gen_m2/(CF*pow(gen_pt2,CF)),2),w);
	      //
	      pmreco->Fill(ptave,m1/ptave,w);
	      pmreco->Fill(ptave,m2/ptave,w);
	      pmgen->Fill(ptave,gen_m1/ptave,w);
	      pmgen->Fill(ptave,gen_m2/ptave,w);
	      pavegen->Fill(ptave,0.5*(gen_pt1+gen_pt2)/ptave,w);
	      pavereco->Fill(ptave,ptave,w);
	      pavegencf->Fill(ptave,CF*pow(0.5*(gen_pt1+gen_pt2),CF)/ptave,w);
	      paverecocf->Fill(ptave,CF*pow(ptave,CF)/ptave,w);
	    } // gooddr1, gooddr2

	    // TOP
	    double ptt1 = t1.Pt();
	    double mt1 = t1.M();
	    double mt2 = t2.M();
	    double mwb = Wb.M();
	    double mwB = WB.M();
	    double mt1up = t1up.M();
	    double mt1dw = t1dw.M();
	    double mt2up = t2up.M();
	    double mt2dw = t2dw.M();
	    double mtmin = min(mt1,mt2);
	    double mtmax = max(mt1,mt2);
	    double ptbqave1 = 0.5*(bpt1 + 0.5*(pt1+pt2));
	    double ptbqave2 = 0.5*(bpt2 + 0.5*(pt1+pt2));
	    bool goodmt1 = (mt1>125 && mt1<220);
	    bool goodmt2 = (mt2>125 && mt2<220);
	    bool goodmt1a = (mt1>130 && mt1<220);
	    bool goodmt2a = (mt2>130 && mt2<220);
	    bool goodmt1up = (mt1up>125 && mt1up<220);
	    bool goodmt1dw = (mt1dw>125 && mt1dw<220);
	    bool goodmt2up = (mt2up>125 && mt2up<220);
	    bool goodmt2dw = (mt2dw>125 && mt2dw<220);
	    bool goodbdr1 = (pow(beta1-gen_beta1,2)+pow(bphi1-gen_bphi1,2) < 0.2*0.2);
	    bool goodbdr2 = (pow(beta2-gen_beta2,2)+pow(bphi2-gen_bphi2,2) < 0.2*0.2);
	    if (isDT) TruePU = getTruePU(run,lumiBlock);
	    if (goodw0) pmuiov->Fill(run+0.5,TruePU,w);
	    if (goodw0) pnpviov->Fill(run+0.5,NPrVtx,w);
	    if (goodw0) prhoiov->Fill(run+0.5,pfRho,w);
	    //
	    if (goodw0) pmuvsmu->Fill(TruePU,TruePU,w);
	    if (goodw0) pnpvvsmu->Fill(TruePU,NPrVtx,w);
	    if (goodw0) prhovsmu->Fill(TruePU,pfRho,w);
	    
	    if (goodw0) pmwiov->Fill(run+0.5,recoWmass0,w);

	    // TOP1
	    if (fabs(beta1)<etaMax && goodw0) {

	      if (goodmt1) pmtiov->Fill(run+0.5,mt1,w);
	      if (goodmt1) pmt1iov->Fill(run+0.5,mt1,w);
	      if (goodmt1) ppt1iov->Fill(run+0.5,ptt1,w);
	      if (goodmt1) plepiov->Fill(run+0.5,lep_pt,w);
	      h1mt->Fill(     mt1,w);
	      h2mt->Fill(bpt1,mt1,w);
	      if (goodmt1)   pmt  ->Fill(bpt1,mt1,w);
	      if (goodmt1up) pmtup->Fill(bpt1,mt1up,w);
	      if (goodmt1dw) pmtdw->Fill(bpt1,mt1dw,w);
	      if (goodmt1)   pmta  ->Fill(ptbqave1,mt1,w);
	      //
	      if (goodmt1) h1pt1->Fill(ptt1,w);
	      if (goodmt1) h1lep->Fill(lep_pt,w);
	      h1mt1->Fill(     mt1,w);
	      h2mt1->Fill(bpt1,mt1,w);
	      if (goodmt1) pmt1 ->Fill(bpt1,mt1,w);
	      if (mt1<mt2) {
		h1mtmin->Fill(     mt1,w);
		h2mtmin->Fill(bpt1,mt1,w);
		if (goodmt1) pmtmin ->Fill(bpt1,mt1,w);
	      }
	      if (mt1>mt2) {
		h1mtmax->Fill(     mt1,w);
		h2mtmax->Fill(bpt1,mt1,w);
		if (goodmt1) pmtmax ->Fill(bpt1,mt1,w);
	      }
	      // bJES1
	      if (goodmt1 && gen_bpt1>0) {
		h1rb->Fill(gen_bpt1/(kBNU1*bpt1), w);
		h2rb->Fill(bpt1, gen_bpt1/(kBNU1*bpt1), w);
		prb ->Fill(bpt1, gen_bpt1/(kBNU1*bpt1), w);
		prb0->Fill(bpt1, gen_bpt1/bpt1, w);
		//
		h1gb->Fill((kBNU1*bpt1)/gen_bpt1, w);
		h2gb->Fill(gen_bpt1, (kBNU1*bpt1)/gen_bpt1, w);
		pgb ->Fill(gen_bpt1, (kBNU1*bpt1)/gen_bpt1, w);
		pgb0->Fill(gen_bpt1, bpt1/gen_bpt1, w);
		if (goodbdr1) {
		  pmbreco->Fill(bpt1,bm1/bpt1,w);
		  pmbgen->Fill(bpt1,gen_bm1/bpt1,w);
		  pbgen->Fill(bpt1,gen_bpt1/bpt1,w);
		  pbreco->Fill(bpt1,bpt1,w);
		}
	      }
	      if (goodmt1) {
		pmbjet->Fill(bpt1,bm1/bpt1,w);
	      }
	    } // TOP1
	    // TOP2
	    if (fabs(beta2)<etaMax && goodw0) {

	      if (goodmt2) pmtiov->Fill(run+0.5,mt2,w);
	      h1mt->Fill(     mt2,w);
	      h2mt->Fill(bpt2,mt2,w);
	      if (goodmt2)   pmt  ->Fill(bpt2,mt2,w);
	      if (goodmt2up) pmtup->Fill(bpt2,mt2up,w);
	      if (goodmt2dw) pmtdw->Fill(bpt2,mt2dw,w);
	      if (goodmt2)   pmta  ->Fill(ptbqave2,mt2,w);
	      //
	      h1mt2->Fill(     mt2,w);
	      h2mt2->Fill(bpt2,mt2,w);
	      if (goodmt2) pmt2 ->Fill(bpt2,mt2,w);
	      if (mt2<mt1) {
		h1mtmin->Fill(     mt2,w);
		h2mtmin->Fill(bpt2,mt2,w);
		if (goodmt2) pmtmin ->Fill(bpt2,mt2,w);
	      }
	      if (mt2>mt1) {
		h1mtmax->Fill(     mt2,w);
		h2mtmax->Fill(bpt2,mt2,w);
		if (goodmt2) pmtmax ->Fill(bpt2,mt2,w);
	      }
	      // bJES
	      if (goodmt2 && gen_bpt2>0) {
		h1rb->Fill(gen_bpt2/(kBNU2*bpt2), w);
		h2rb->Fill(bpt2, gen_bpt2/(kBNU2*bpt2), w);
		prb ->Fill(bpt2, gen_bpt2/(kBNU2*bpt2), w);
		//
		h1gb->Fill((kBNU2*bpt2)/gen_bpt2, w);
		h2gb->Fill(gen_bpt2, (kBNU2*bpt2)/gen_bpt2, w);
		pgb ->Fill(gen_bpt2, (kBNU2*bpt2)/gen_bpt2, w);
		pmbreco->Fill(bpt2,bm2/bpt2,w);
		if (goodbdr2) {
		  pmbreco->Fill(bpt2,bm2/bpt2,w);
		  pmbgen->Fill(bpt2,gen_bm2/bpt2,w);
		  pbgen->Fill(bpt2,gen_bpt2/bpt2,w);
		  pbreco->Fill(bpt2,bpt2,w);
		}
	      }
	      if (goodmt2) {
		pmbjet->Fill(bpt2,bm2/bpt2,w);
	      }
	    } // TOP2

	    // bJES Rbq attempts
	    if (fabs(beta1)<etaMax && fabs(beta2)<etaMax && goodw0) {

	      double rbq = (bpt1+bpt2)/(pt1+pt2);
	      double bptave = 0.5*(bpt1+bpt2);
	      double aptave = 0.5*(ptave + bptave);
	      double cptave = 0.5*(1.2*ptave + bptave);
	      double drbq = ((bpt1+bpt2) - (pt1+pt2)) / (4. * aptave);
	      double drbqc = ((bpt1+bpt2) - 1.2*(pt1+pt2)) / (4. * cptave);
	      hrbq->Fill(rbq,w);
	      h2bq->Fill(ptave,bptave,w); 
	      h2rbqq->Fill(ptave,rbq,w); 
	      prbqq->Fill(ptave,rbq,w); 
	      h2rbqb->Fill(bptave,1./rbq,w); 
	      prbqb->Fill(bptave,1./rbq,w); 
	      h2rbqa->Fill(aptave,1+drbq,w); 
	      prbqa->Fill(aptave,1+drbq,w); 
	      h2rbqc->Fill(cptave,1+drbqc,w); 
	      prbqc->Fill(cptave,1+drbqc,w); 

	      if (gen_bpt1>0 && fabs(bpt1/gen_bpt1-1)<0.5 &&
		  gen_bpt2>0 && fabs(bpt2/gen_bpt2-1)<0.5) {
		pb1->Fill(gen_bpt1,bpt1/gen_bpt1,w);
		pb1->Fill(gen_bpt2,bpt2/gen_bpt2,w);
		pb2->Fill(gen_bpt1,pow(bpt1/gen_bpt1,2),w);
		pb2->Fill(gen_bpt2,pow(bpt2/gen_bpt2,2),w);
		pbq1->Fill(ptave,bpt1/gen_bpt1,w);
		pbq1->Fill(ptave,bpt2/gen_bpt2,w);
		pbq2->Fill(ptave,pow(bpt1/gen_bpt1,2),w);
		pbq2->Fill(ptave,pow(bpt2/gen_bpt2,2),w);
		pbb1->Fill(bptave,bpt1/gen_bpt1,w);
		pbb1->Fill(bptave,bpt2/gen_bpt2,w);
		pbb2->Fill(bptave,pow(bpt1/gen_bpt1,2),w);
		pbb2->Fill(bptave,pow(bpt2/gen_bpt2,2),w);
		pba1->Fill(aptave,bpt1/gen_bpt1,w);
		pba1->Fill(aptave,bpt2/gen_bpt2,w);
		pba2->Fill(aptave,pow(bpt1/gen_bpt1,2),w);
		pba2->Fill(aptave,pow(bpt2/gen_bpt2,2),w);
	      }
	      
	    } // Rbq

	    // ATLAS-style plots. Tighter eta cuts, though,
	    // and corrections for jets
	    double rbq = (kBNU1*bpt1+kBNU2*bpt2)/(kQNU1*pt1+kQNU2*pt2);
	    double rbb = (kBNU1*bpt1+kBNU2*bpt2);
	    double rqq = (kQNU1*pt1+kQNU2*pt2);
	    double rb1 = (gen_bpt1!=0 ? kBNU1*bpt1 / gen_bpt1 : 0);
	    double rb2 = (gen_bpt2!=0 ? kBNU2*bpt2 / gen_bpt2 : 0);
	    double rq1 = (gen_pt1!=0 ? kQNU1*pt1 / gen_pt1 : 0);
	    double rq2 = (gen_pt2!=0 ? kQNU2*pt2 / gen_pt2 : 0);
	    double ptqq = qq.Pt();
	    double ptb = kBNU1*bpt1;
	    double ptB = ptb * (80.4 / qq.M());
	    double ptw = W.Pt();
	    bool goodrbqa = (rbq>0.3 && rbq<3.0);
	    double mlb2 = lb2.M();
	    bool goodmlb2a = (mlb2>30 && mlb2<170);
	    double mlb2met = lb2met.M();

	    // pT5 and ISR variables
	    double pttt = (t1+lb2met).Pt();

	    double gmt1 = gt1.M();
	    double gmwb = gWb.M();
	    double gmwB = gWB.M();
	    double grbq = (gq1.Pt()+gq2.Pt()>0 ? 
			   (gb1.Pt()+gb2.Pt())/(gq1.Pt()+gq2.Pt()) : 0);
	    double grbb = (gb1.Pt()+gb2.Pt());
	    double grqq = (gq1.Pt()+gq2.Pt());
	    double gptb = gb1.Pt();
	    double gptB = gb1.Pt() * (80.4 / gmw);
	    double gptqq = gw.Pt();
	    double gptw = gW.Pt();
	    double gmlb2 = glb2.M();

	    if (fabs(eta1)<etaMax && fabs(eta2)<etaMax &&
		fabs(beta1)<etaMax && fabs(beta2)<etaMax &&
		goodwa && goodrbqa && goodmt1a && goodmlb2a) {

	      // ATLAS 3D variables
	      hatmw->Fill(recoWmass0,w);
	      hatbq->Fill(rbq,w);
	      hatmt->Fill(mt1,w);
	      hatlb->Fill(mlb2,w);
	      
	      // Gap fraction, pT5 absolute and relative pT
	      hatlbmet->Fill(mlb2met,w);
	      for (int i = 0; i != pgap->GetNbinsX()+1; ++i) {
		double ptgap = pgap->GetBinLowEdge(i);
		pgap->Fill(ptgap+1e-3, gpt1<ptgap ? 1 : 0, w);
	      }
	      hatpt5->Fill(gpt1,w);
	      hatpttt->Fill(pttt,w);
	      hatpt5ott->Fill(pttt!=0 ? gpt1/pttt : 0,w);
	      h2atpt5vtt->Fill(pttt, gpt1, w);
	      if (gpt1>0) patpt5ott->Fill(pttt, pttt!=0 ? gpt1/pttt : 0., w);

	      // Weight variants, SIdx==0 means signal
	      // Could also require right combotype? ComboType==1
	      if (SIdx==0) hatmwsw->Fill(recoWmass0,w);
	      if (SIdx==0) hatmws1->Fill(recoWmass0,1);
	      if (SIdx==0) hatbqsw->Fill(rbq,w);
	      if (SIdx==0) hatbqs1->Fill(rbq,1);

	      // QGL>0.5 variants
	      bool isq1 = (qgl1>=0.5);
	      bool isq2 = (qgl2>=0.5);
	      bool isg1 = (qgl1>=0 && qgl1<0.5);
	      bool isg2 = (qgl2>=0 && qgl2<0.5);
	      bool isqq = (isq1 && isq2);
	      bool isqg = (isq1||isq2) && (isg1||isg2);
	      bool isgg = (isg1 && isg2);
	      if (isqq)    hatmwqq->Fill(recoWmass0,w);
	      if (isqg) hatmwqg->Fill(recoWmass0,w);
	      if (isgg)    hatmwgg->Fill(recoWmass0,w);
	      patmwa->Fill(0.5*(kQNU1*pt1+kQNU2*pt2), recoWmass0,w);
	      if (isqq) patmwqqa->Fill(0.5*(kQNU1*pt1+kQNU2*pt2), recoWmass0,w);
	      if (isqg) patmwqga->Fill(0.5*(kQNU1*pt1+kQNU2*pt2), recoWmass0,w);
	      if (isgg) patmwgga->Fill(0.5*(kQNU1*pt1+kQNU2*pt2), recoWmass0,w);
	      patmwb->Fill(kQNU1*pt1, recoWmass0,w);
	      patmwb->Fill(kQNU2*pt2, recoWmass0,w);
	      if (isq1) patmwqb->Fill(kQNU1*pt1, recoWmass0,w);
	      if (isq2) patmwqb->Fill(kQNU2*pt2, recoWmass0,w);
	      if (isg1) patmwgb->Fill(kQNU1*pt1, recoWmass0,w);
	      if (isg2) patmwgb->Fill(kQNU2*pt2, recoWmass0,w);

	      // Semileptonic BRs
	      pqelfa->Fill(kQNU1*pt1, elf1, w);
	      pqelfa->Fill(kQNU2*pt2, elf2, w);
	      if (elf1>0) pqelfb->Fill(kQNU1*pt1, elf1, w);
	      if (elf2>0) pqelfb->Fill(kQNU2*pt2, elf2, w);
	      pqelfn->Fill(kQNU1*pt1, elf1>0 ? 1 : 0, w);
	      pqelfn->Fill(kQNU2*pt2, elf2>0 ? 1 : 0, w);
	      //
	      pqmufa->Fill(kQNU1*pt1, muf1, w);
	      pqmufa->Fill(kQNU2*pt2, muf2, w);
	      if (muf1>0) pqmufb->Fill(kQNU1*pt1, muf1, w);
	      if (muf2>0) pqmufb->Fill(kQNU2*pt2, muf2, w);
	      pqmufn->Fill(kQNU1*pt1, muf1>0 ? 1 : 0, w);
	      pqmufn->Fill(kQNU2*pt2, muf2>0 ? 1 : 0, w);

	      pbelfa->Fill(kBNU1*bpt1, belf1, w);
	      pbelfa->Fill(kBNU2*bpt2, belf2, w);
	      if (belf1>0) pbelfb->Fill(kBNU1*bpt1, belf1, w);
	      if (belf2>0) pbelfb->Fill(kBNU2*bpt2, belf2, w);
	      pbelfn->Fill(kBNU1*bpt1, belf1>0 ? 1 : 0, w);
	      pbelfn->Fill(kBNU2*bpt2, belf2>0 ? 1 : 0, w);
	      //
	      pbmufa->Fill(kBNU1*bpt1, bmuf1, w);
	      pbmufa->Fill(kBNU2*bpt2, bmuf2, w);
	      if (bmuf1>0) pbmufb->Fill(kBNU1*bpt1, bmuf1, w);
	      if (bmuf2>0) pbmufb->Fill(kBNU2*bpt2, bmuf2, w);
	      pbmufn->Fill(kBNU1*bpt1, bmuf1>0 ? 1 : 0, w);
	      pbmufn->Fill(kBNU2*bpt2, bmuf2>0 ? 1 : 0, w);

	      // QGL shape for quarks
	      hqglq->Fill(qgl1, w);
	      hqglq->Fill(qgl2, w);
	      hqglb->Fill(bqgl1, w);
	      hqglb->Fill(bqgl2, w);
	      h2qglqa->Fill(qgl1, 0.5*(kQNU1*pt1+kQNU2*pt2), w);
	      h2qglqa->Fill(qgl2, 0.5*(kQNU1*pt1+kQNU2*pt2), w);
	      h2qglqb->Fill(qgl1, kQNU1*pt1, w);
	      h2qglqb->Fill(qgl2, kQNU2*pt2, w);
	      if (abs(flav1)<5 && flav1!=0)
		                 h2qglqb_q->Fill(qgl1, kQNU1*pt1, w);
	      if (abs(flav2)<5 && flav2!=0)
		                 h2qglqb_q->Fill(qgl2, kQNU2*pt2, w);
	      if (flav1==21)     h2qglqb_g->Fill(qgl1, kQNU1*pt1, w);
	      if (flav2==21)     h2qglqb_g->Fill(qgl2, kQNU2*pt2, w);
	      if (abs(flav1)==1) h2qglqb_d->Fill(qgl1, kQNU1*pt1, w);
	      if (abs(flav2)==1) h2qglqb_d->Fill(qgl2, kQNU2*pt2, w);
	      if (abs(flav1)==2) h2qglqb_u->Fill(qgl1, kQNU1*pt1, w);
	      if (abs(flav2)==2) h2qglqb_u->Fill(qgl2, kQNU2*pt2, w);
	      if (abs(flav1)==3) h2qglqb_s->Fill(qgl1, kQNU1*pt1, w);
	      if (abs(flav2)==3) h2qglqb_s->Fill(qgl2, kQNU2*pt2, w);
	      if (abs(flav1)==4) h2qglqb_c->Fill(qgl1, kQNU1*pt1, w);
	      if (abs(flav2)==4) h2qglqb_c->Fill(qgl2, kQNU2*pt2, w);
	      if (abs(flav1)==5) h2qglqb_b->Fill(qgl1, kQNU1*pt1, w);
	      if (abs(flav2)==5) h2qglqb_b->Fill(qgl2, kQNU2*pt2, w);
	      if (flav1==0)      h2qglqb_o->Fill(qgl1, kQNU1*pt1, w);
	      if (flav2==0)      h2qglqb_o->Fill(qgl2, kQNU2*pt2, w);
	      h2qglba->Fill(bqgl1, 0.5*(kBNU1*bpt1+kBNU2*bpt2), w);
	      h2qglba->Fill(bqgl2, 0.5*(kBNU1*bpt1+kBNU2*bpt2), w);
	      h2qglbb->Fill(bqgl1, kBNU1*bpt1, w);
	      h2qglbb->Fill(bqgl2, kBNU2*bpt2, w);
	      if (abs(bflav1)==5) h2qglbb_b->Fill(bqgl1, kBNU1*bpt1, w);
	      if (abs(bflav2)==5) h2qglbb_b->Fill(bqgl2, kBNU2*bpt2, w);
	      // ...and gluons
	      if (gpt1>0) {
		hqglg->Fill(gqgl1, w);
		h2qglg->Fill(gqgl1, gpt1, w);
		if (gflav1==21)     h2qglg_g->Fill(gqgl1, gpt1, w);
		if (abs(gflav1)<5 && gflav1!=0)
		                    h2qglg_q->Fill(gqgl1, gpt1, w);
		if (abs(gflav1)==1) h2qglg_u->Fill(gqgl1, gpt1, w);
		if (abs(gflav1)==2) h2qglg_d->Fill(gqgl1, gpt1, w);
		if (abs(gflav1)==3) h2qglg_s->Fill(gqgl1, gpt1, w);
		if (abs(gflav1)==4) h2qglg_c->Fill(gqgl1, gpt1, w);
		if (abs(gflav1)==5) h2qglg_b->Fill(gqgl1, gpt1, w);
		if (gflav1==0)      h2qglg_o->Fill(gqgl1, gpt1, w);
	      }

	      // Quark / gluon fractions
	      if (abs(flav1)<5 && flav1!=0)
		                 peqq->Fill(kQNU1*pt1, qgl1>0.5 ? 1 : 0, w);
	      if (abs(flav2)<5 && flav2!=0)
		                 peqq->Fill(kQNU2*pt2, qgl2>0.5 ? 1 : 0, w);
	      if (abs(flav1)==1) peqd->Fill(kQNU1*pt1, qgl1>0.5 ? 1 : 0, w);
	      if (abs(flav2)==1) peqd->Fill(kQNU2*pt2, qgl2>0.5 ? 1 : 0, w);
	      if (abs(flav1)==2) pequ->Fill(kQNU1*pt1, qgl1>0.5 ? 1 : 0, w);
	      if (abs(flav2)==2) pequ->Fill(kQNU2*pt2, qgl2>0.5 ? 1 : 0, w);
	      if (abs(flav1)==3) peqs->Fill(kQNU1*pt1, qgl1>0.5 ? 1 : 0, w);
	      if (abs(flav2)==3) peqs->Fill(kQNU2*pt2, qgl2>0.5 ? 1 : 0, w);
	      if (abs(flav1)==4) peqc->Fill(kQNU1*pt1, qgl1>0.5 ? 1 : 0, w);
	      if (abs(flav2)==4) peqc->Fill(kQNU2*pt2, qgl2>0.5 ? 1 : 0, w);
	      if (abs(flav1)==5) peqb->Fill(kQNU1*pt1, qgl1>0.5 ? 1 : 0, w);
	      if (abs(flav2)==5) peqb->Fill(kQNU2*pt2, qgl2>0.5 ? 1 : 0, w);
	      if (flav1==21)     peqg->Fill(kQNU1*pt1, qgl1>0.5 ? 1 : 0, w);
	      if (flav2==21)     peqg->Fill(kQNU2*pt2, qgl2>0.5 ? 1 : 0, w);
	      if (flav1==0)      peqo->Fill(kQNU1*pt1, qgl1>0.5 ? 1 : 0, w);
	      if (flav2==0)      peqo->Fill(kQNU2*pt2, qgl2>0.5 ? 1 : 0, w);

	      if (abs(bflav1)==5)  pebb->Fill(kBNU1*bpt1, bqgl1>0.5 ? 1 : 0, w);
	      if (abs(bflav2)==5)  pebb->Fill(kBNU2*bpt2, bqgl2>0.5 ? 1 : 0, w);
	      if (abs(bflav1)<5 && bflav1!=0)
		                   pebq->Fill(kBNU1*bpt1, bqgl1>0.5 ? 1 : 0, w);
	      if (abs(bflav2)<5 && bflav2!=0)
		                   pebq->Fill(kBNU2*bpt2, bqgl2>0.5 ? 1 : 0, w);
	      if (abs(bflav1)==21) pebg->Fill(kBNU1*bpt1, bqgl1>0.5 ? 1 : 0, w);
	      if (abs(bflav2)==21) pebg->Fill(kBNU2*bpt2, bqgl2>0.5 ? 1 : 0, w);
	      if (bflav1==0)       pebo->Fill(kBNU1*bpt1, bqgl1>0.5 ? 1 : 0, w);
	      if (bflav2==0)       pebo->Fill(kBNU2*bpt2, bqgl2>0.5 ? 1 : 0, w);

	      if (gpt1>0) {
		if (gflav1==21)     pegg->Fill(gpt1, gqgl1>0.5 ? 1 : 0, w);
		if (abs(gflav1)<5 && gflav1!=0)
		                    pegq->Fill(gpt1, gqgl1>0.5 ? 1 : 0, w);
		if (abs(gflav1)==1) pegd->Fill(gpt1, gqgl1>0.5 ? 1 : 0, w);
		if (abs(gflav1)==2) pegu->Fill(gpt1, gqgl1>0.5 ? 1 : 0, w);
		if (abs(gflav1)==3) pegs->Fill(gpt1, gqgl1>0.5 ? 1 : 0, w);
		if (abs(gflav1)==4) pegc->Fill(gpt1, gqgl1>0.5 ? 1 : 0, w);
		if (abs(gflav1)==5) pegb->Fill(gpt1, gqgl1>0.5 ? 1 : 0, w);
		if (gflav1==0)      pego->Fill(gpt1, gqgl1>0.5 ? 1 : 0, w);
	      }		

	      pfqq->Fill(kQNU1*pt1, abs(flav1)<5 && flav1!=0 ? 1 : 0, w);
	      pfqq->Fill(kQNU2*pt2, abs(flav2)<5 && flav2!=0 ? 1 : 0, w);
	      pfqd->Fill(kQNU1*pt1, abs(flav1)==1 ? 1 : 0, w);
	      pfqd->Fill(kQNU2*pt2, abs(flav2)==1 ? 1 : 0, w);
	      pfqu->Fill(kQNU1*pt1, abs(flav1)==2 ? 1 : 0, w);
	      pfqu->Fill(kQNU2*pt2, abs(flav2)==2 ? 1 : 0, w);
	      pfqs->Fill(kQNU1*pt1, abs(flav1)==3 ? 1 : 0, w);
	      pfqs->Fill(kQNU2*pt2, abs(flav2)==3 ? 1 : 0, w);
	      pfqc->Fill(kQNU1*pt1, abs(flav1)==4 ? 1 : 0, w);
	      pfqc->Fill(kQNU2*pt2, abs(flav2)==4 ? 1 : 0, w);
	      pfqb->Fill(kQNU1*pt1, abs(flav1)==5 ? 1 : 0, w);
	      pfqb->Fill(kQNU2*pt2, abs(flav2)==5 ? 1 : 0, w);
	      pfqg->Fill(kQNU1*pt1, flav1==21 ? 1 : 0, w);
	      pfqg->Fill(kQNU2*pt2, flav2==21 ? 1 : 0, w);
	      pfqo->Fill(kQNU1*pt1, flav1==0 ? 1 : 0, w);
	      pfqo->Fill(kQNU2*pt2, flav2==0 ? 1 : 0, w);

	      pfbb->Fill(kBNU1*bpt1, abs(bflav1)==5 ? 1 : 0, w);
	      pfbb->Fill(kBNU2*bpt2, abs(bflav2)==5 ? 1 : 0, w);
	      pfbq->Fill(kBNU1*bpt1, abs(bflav1)<5 && bflav1!=0 ? 1 : 0, w);
	      pfbq->Fill(kBNU2*bpt2, abs(bflav2)<5 && bflav2!=0 ? 1 : 0, w);
	      pfbg->Fill(kBNU1*bpt1, bflav1==21 ? 1 : 0, w);
	      pfbg->Fill(kBNU2*bpt2, bflav2==21 ? 1 : 0, w);
	      pfbo->Fill(kBNU1*bpt1, bflav1==0 ? 1 : 0, w);
	      pfbo->Fill(kBNU2*bpt2, bflav2==0 ? 1 : 0, w);

	      if (gpt1>0) {
		pfgg->Fill(gpt1, gflav1==21 ? 1 : 0, w);
		pfgq->Fill(gpt1, abs(gflav1)<5 && gflav1!=0 ? 1 : 0, w);
		pfgd->Fill(gpt1, abs(gflav1)==1 ? 1 : 0, w);
		pfgu->Fill(gpt1, abs(gflav1)==2 ? 1 : 0, w);
		pfgs->Fill(gpt1, abs(gflav1)==3 ? 1 : 0, w);
		pfgc->Fill(gpt1, abs(gflav1)==4 ? 1 : 0, w);
		pfgb->Fill(gpt1, abs(gflav1)==5 ? 1 : 0, w);
		pfgo->Fill(gpt1, gflav1==0 ? 1 : 0, w);
	      }

	      // Top pT
	      htopptorig->Fill(hadTopPt,w);
	      htopptgen->Fill(gt1.Pt(),w);
	      htoppt->Fill(ptt1,w);
	      htopptw0->Fill(ptt1,w_nopt);
	      if (isMC) {
		htoppt_fsr->Fill(ptt1, (*PSWgts)[4]/(*PSWgts)[0]);
		htoppt_q2qg->Fill(ptt1, (*PSWgts)[12]/(*PSWgts)[0]);
		htoppt_x2xg->Fill(ptt1, (*PSWgts)[14]/(*PSWgts)[0]);
		htoppt_isr->Fill(ptt1, (*PSWgts)[26]/(*PSWgts)[0]);
	      }		

	      // MPF components: 0=MHT(nu), 1=lep, 2=b2, 3=b1, 4=q1+q2 (, 5=t1)
	      // bins 6,-0.5,5.5);
	      //TLorentzVector ref = (nu * (1./nu.Mag()));
	      TLorentzVector ref;
	      ref.SetPtEtaPhiM(nu.Pt(),0,nu.Phi(),0);
	      ref *= (ref.Pt()!=0 ? 1./nu.Pt() : 1);
	      pmpf->Fill(0., nu.Dot(ref),w);
	      pmpf->Fill(1., l.Dot(ref),w);
	      pmpf->Fill(2., b2.Dot(ref),w);
	      pmpf->Fill(3., b1.Dot(ref),w);
	      pmpf->Fill(4., qq.Dot(ref),w);
	      pmpf->Fill(5., t1.Dot(ref),w);
	      /*
	      TVector3 ref(nu.Px(),nu.Py(),0);
	      ref *= (nu.Pt()!= 0 ? 1./nu.Pt() : 1);
	      pmpf->Fill(0., ref * nu,w);
	      pmpf->Fill(1., ref * l,w);
	      pmpf->Fill(2., ref * b2,w);
	      pmpf->Fill(3., ref * b1,w);
	      pmpf->Fill(4., ref * qq,w);
	      pmpf->Fill(5., ref * t1,w);
	      */
	      pbqiov->Fill(run+0.5,rbq,w);
	      plbiov->Fill(run+0.5,mlb2,w);

	      // Fill PS weights variations with ATLAS cuts
	      // all jets in barrel, good mw, mt, mlb plus fitProb>0.2
	      if (isMC) {
		const int nps = 46;
		assert(PSWgts->size()==nps);
		for (int ips = 1; ips != nps; ++ips) {

		  //double wps = (*PSWgts)[ips]/(*PSWgts)[0];
		  double wps = (*PSWgts)[ips]/(*PSWgts)[0] * w;

		  // Reco observables
		  psmw->Fill(ips, recoWmass0, wps);
		  psmt->Fill(ips, mt1, wps);
		  psmwb->Fill(ips, mwb, wps);
		  psmwB->Fill(ips, mwB, wps);
		  psmlb->Fill(ips, mlb2, wps);
		  psrbq->Fill(ips, rbq, wps);
		  psrbb->Fill(ips, rbb, wps);
		  psrqq->Fill(ips, rqq, wps);
		  psptb->Fill(ips, ptb, wps);
		  psptB->Fill(ips, ptB, wps);
		  psptqq->Fill(ips, ptqq, wps);
		  psptw->Fill(ips, ptw, wps);
		  psptt->Fill(ips, ptt1, wps);
		  if (fabs(rb1-1)<0.5) psrb->Fill(ips, rb1, wps);
		  if (fabs(rb2-1)<0.5) psrb->Fill(ips, rb2, wps);
		  if (fabs(rq1-1)<0.5) psrq->Fill(ips, rq1, wps);
		  if (fabs(rq2-1)<0.5) psrq->Fill(ips, rq2, wps);

		  // Gen observables to factor out changes in JES
		  if (goodgen) {
		    psgmw->Fill(ips, gmw, wps);
		    psgmt->Fill(ips, gmt1, wps);
		    psgmwb->Fill(ips, gmwb, wps);
		    psgmwB->Fill(ips, gmwB, wps);
		    psgmlb->Fill(ips, gmlb2, wps);
		    psgrbq->Fill(ips, grbq, wps);
		    psgrbb->Fill(ips, grbb, wps);
		    psgrqq->Fill(ips, grqq, wps);
		    psgptb->Fill(ips, gptb, wps);
		    psgptB->Fill(ips, gptB, wps);
		    psgptqq->Fill(ips, gptqq, wps);
		    psgptw->Fill(ips, gptw, wps);
		    psgptt->Fill(ips, gt1.Pt(), wps);
		  }
		} // for ips
	      }
	    } // ATLAS cuts

	    // cs(+cd) vs ud(+us) scale

	    // true flavor pairs
	    bool cs = ((flav1==+4&&flav2==-3) || (flav1==-3&&flav2==+4) || 
		       (flav1==-4&&flav2==+3) || (flav1==+3&&flav2==-4)); 
	    bool cd = ((flav1==+4&&flav2==-1) || (flav1==-1&&flav2==+4) || 
		       (flav1==-4&&flav2==+1) || (flav1==+1&&flav2==-4)); 
	    bool ud = ((flav1==+2&&flav2==-1) || (flav1==-1&&flav2==+2) || 
		       (flav1==-2&&flav2==+1) || (flav1==+1&&flav2==-2)); 
	    bool us = ((flav1==+2&&flav2==-3) || (flav1==-3&&flav2==+2) || 
		       (flav1==-2&&flav2==+3) || (flav1==+3&&flav2==-2));
	    //bool gx = (flav1==21||flav2==22); // bug found 20201218
	    bool gx = (flav1==21||flav2==21);
	    bool ot = !(cs||cd||ud||us||gx);

	    // truth tags
	    bool isL = (ud||us);
	    bool isH = (cs||cd);
	    bool isG = gx;
	    bool isO = ot;
	    
	    int tagLHGO(0);
	    if (isL) tagLHGO = 1;
	    if (isH) tagLHGO = 2;
	    if (isG) tagLHGO = 3;
	    if (isO) tagLHGO = 4;

	    // tight ctag
	    // https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL18
	    /*
	    bool ct1 = (ctag1/(ctag1+udstag1+gtag1)>0.282 &&
			ctag1/(ctag1+btag1+bleptag1)>0.267); // ctag1
	    bool ct2 = (ctag2/(ctag2+udstag2+gtag2)>0.282 &&
			ctag2/(ctag2+btag2+bleptag2)>0.267); //ctag2
	    */
	    // Only tag vs light, because W>qq' has very few b
	    bool ct1 = (ctag1/(ctag1+udstag1+gtag1)>0.282); // ctag1
	    bool ct2 = (ctag2/(ctag2+udstag2+gtag2)>0.282); // ctag1
	    bool ct = (ct1||ct2); // ctag

	    // medium ctag
	    // https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL18
	    bool cm1 = (ctag1/(ctag1+udstag1+gtag1)>0.099 &&
			ctag1/(ctag1+btag1+bleptag1)>0.325); // ctag1
	    bool cm2 = (ctag2/(ctag2+udstag2+gtag2)>0.099 &&
			ctag2/(ctag2+btag2+bleptag2)>0.325); //ctag2
	    bool cm = (cm1||cm2); // ctag

	    // gluon tag
	    bool cg1 = (gtag1>0.8);
	    bool cg2 = (gtag2>0.8);
	    bool cg = (cg1||cg2);

	    // reco tags
	    bool isD = cg;
	    bool isC = ct && !cg;
	    bool isB = cm && !ct && !cg;
	    bool isA = !cm && !ct && !cg;
	    int tagABCD(0);
	    if (isA) tagABCD = 1;
	    if (isB) tagABCD = 2;
	    if (isC) tagABCD = 3;
	    if (isD) tagABCD = 4;

	    // flavor correlationsx
	    h2LHGO->Fill(min(6,flav1),min(6,flav2),w);
	    if (isL) h2L->Fill(min(6,flav1),min(6,flav2),w);
	    if (isH) h2H->Fill(min(6,flav1),min(6,flav2),w);
	    if (isG) h2G->Fill(min(6,flav1),min(6,flav2),w);
	    if (isO) h2O->Fill(min(6,flav1),min(6,flav2),w);

	    // W mass distribution for flavor pairs
	    hmLHGO->Fill(recoWmass0,w);
	    if (isL) hmL->Fill(recoWmass0,w);
	    if (isH) hmH->Fill(recoWmass0,w);
	    if (isG) hmH->Fill(recoWmass0,w);
	    if (isO) hmO->Fill(recoWmass0,w);

	    // W mass distribution for tag regions
	    hmABCD->Fill(recoWmass0,w);
	    if (isA) hmA->Fill(recoWmass0,w);
	    if (isB) hmB->Fill(recoWmass0,w);
	    if (isC) hmC->Fill(recoWmass0,w);
	    if (isD) hmD->Fill(recoWmass0,w);

	    // Number of tags in each ABCD region
	    hnLHGO->Fill(tagABCD,w); hnLHGO->Fill(0.,w);
	    if (isL) { hnL->Fill(tagABCD,w); hnL->Fill(0.,w); }
	    if (isH) { hnH->Fill(tagABCD,w); hnH->Fill(0.,w); }
	    if (isG) { hnG->Fill(tagABCD,w); hnG->Fill(0.,w); }
	    if (isO) { hnO->Fill(tagABCD,w); hnO->Fill(0.,w); }

	    // Number of jets in each LHGO  region
	    hnABCD->Fill(tagLHGO,w); hnABCD->Fill(0.,w);
	    if (isA) { hnA->Fill(tagLHGO,w); hnA->Fill(0.,w); }
	    if (isB) { hnB->Fill(tagLHGO,w); hnB->Fill(0.,w); }
	    if (isC) { hnC->Fill(tagLHGO,w); hnC->Fill(0.,w); }
	    if (isD) { hnD->Fill(tagLHGO,w); hnD->Fill(0.,w); }

	    // 2D counts
	    h2n->Fill(tagABCD,tagLHGO,w);
	    h2n->Fill(tagABCD,0.,w);
	    h2n->Fill(0.,tagLHGO,w);
	    h2n->Fill(0.,0.,w);

	    // 2D mass profiles
	    p2m->Fill(tagABCD,tagLHGO,recoWmass0,w);
	    p2m->Fill(tagABCD,0.,recoWmass0,w);
	    p2m->Fill(0.,tagLHGO,recoWmass0,w);
	    p2m->Fill(0.,0.,recoWmass0,w);

	    // s in A vs c in H region
	    if ((ct1 && !cg2) && !(ct2||cm2||cg2)) {
	      hc->Fill(kQNU1*pt1,w);
	      hs->Fill(kQNU2*pt2,w);
	      hcs->Fill(kQNU1*pt1,w);
	      hcs->Fill(kQNU2*pt2,w);
	      pcs->Fill(0.,kQNU1*pt1,w);
	      pcs->Fill(1.,kQNU2*pt2,w);
	      pcs->Fill(2.,kQNU1*pt1,w);
	      pcs->Fill(2.,kQNU2*pt2,w);
	    }
	    if ((ct2 && !cg2) && !(ct1||cm1||cg1)) {
	      hc->Fill(kQNU2*pt2,w);
	      hs->Fill(kQNU1*pt1,w);
	      hcs->Fill(kQNU2*pt2,w);
	      hcs->Fill(kQNU1*pt1,w);
	      pcs->Fill(0.,kQNU2*pt2,w);
	      pcs->Fill(1.,kQNU1*pt1,w);
	      pcs->Fill(2.,kQNU2*pt2,w);
	      pcs->Fill(2.,kQNU1*pt1,w);
	    }
	    if (((ct1 && !cg2) && !(ct2||cm2||cg2)) ||
		((ct2 && !cg2) && !(ct1||cm1||cg1))) {
	      if (isH && cs) {
		if (abs(flav1)==4) gcs->Fill(0., gen_pt1,w);
		if (abs(flav2)==4) gcs->Fill(0., gen_pt2,w);
		if (abs(flav1)==3) gcs->Fill(1., gen_pt1,w);
		if (abs(flav2)==3) gcs->Fill(1., gen_pt2,w);
		gcs->Fill(2., gen_pt1,w);
		gcs->Fill(2., gen_pt2,w);
	      }
	    }
	    // u,d both in A
	    if (!(ct1||cm1||cg1) && !(ct2 || cm2 || cg2)) {
	      hu->Fill(kQNU1*pt1,w);
	      hud->Fill(kQNU1*pt1,w);
	      pud->Fill(0.,kQNU1*pt1,w);
	      hd->Fill(kQNU2*pt2,w);
	      hud->Fill(kQNU2*pt2,w);
	      pud->Fill(1.,kQNU2*pt2,w);
	      pud->Fill(2.,kQNU1*pt1,w);
	      pud->Fill(2.,kQNU2*pt2,w);
	      if (isL && ud) {
		if (abs(flav1)==2) gud->Fill(0., gen_pt1,w);
		if (abs(flav2)==2) gud->Fill(0., gen_pt2,w);
		if (abs(flav1)==1) gud->Fill(1., gen_pt1,w);
		if (abs(flav2)==1) gud->Fill(1., gen_pt2,w);
		gud->Fill(2., gen_pt1,w);
		gud->Fill(2., gen_pt2,w);
	      }
	    }
	    hl->Fill(lep_pt,w);
	    plnu->Fill(0., lep_pt,w);
	    //double mht = (qq+b1+b2+l).Pt();
	    double mht = nu.Pt();
	    hnu->Fill(mht,w);
	    hlnu->Fill(lep_pt,w);
	    hlnu->Fill(mht,w);
	    plnu->Fill(1., mht,w);
	    plnu->Fill(2., lep_pt,w);
	    plnu->Fill(2., mht,w);

	    // Response of tags in each ABCD region

	    for (int i = 1; i != hptboth13->GetNbinsX()+1; ++i) {
	      double ptmin = hptboth13->GetBinLowEdge(i);
	      double ptmax = hptboth13->GetBinLowEdge(i+1);
	      double ptmid = 0.5*(ptmin+ptmax);
	      if (pt1>ptmin && pt2>ptmin && pt1<ptmax && pt2<ptmax) {
		hptboth13->Fill(ptmid,w);
		if (goodw) pmw13ptboth->Fill(ptmid,recoWmass,w);
		if (goodw0) pmm13ptboth->Fill(ptmid,recoWmass0,w);
		if (goodw0) hptetaboth13->Fill(ptmid, eta1, w);
		if (goodw0) hptetaboth13->Fill(ptmid, eta2, w);
	      }
	    } // for i
	    for (int i = 1; i != hptboth13b->GetNbinsX()+1; ++i) {
	      double ptmin = hptboth13b->GetBinLowEdge(i);
	      double ptmax = hptboth13b->GetBinLowEdge(i+1);
	      double ptmid = 0.5*(ptmin+ptmax);
	      //double ptave = 0.5*(pt1+pt2);
	      if (pt1>ptmin && pt2>ptmin && pt1<ptmax && pt2<ptmax) {
		hptboth13b->Fill(ptmid,w);
		hptboth13f->Fill(ptave,w);
		pptboth13b->Fill(ptmid,ptave,w);
		if (goodw) pmw13bptboth->Fill(ptmid,recoWmass,w);
		if (goodw0) pmm13bptboth->Fill(ptmid,recoWmass0,w);
		if (goodw0) pmm13bptboth_noL3Res->Fill(ptmid,recoWmassUnc,w);
		hptetaboth13b->Fill(ptmid, eta1,w);
		hptetaboth13b->Fill(ptmid, eta2,w);
	      }
	    } // for i
	    for (int i = 1; i != hptthr13->GetNbinsX()+1; ++i) {
	      double ptmin = hptthr13->GetBinCenter(i);
	      double ptmax = xmax;//150.;
	      if (pt1>ptmin && pt2>ptmin &&
		  pt1<ptmax && pt2<ptmax) {
		hptthr13->Fill(ptmin,w);
		if (goodw) pmw13ptthr->Fill(ptmin,recoWmass,w);
		if (goodw0) pmm13ptthr->Fill(ptmin,recoWmass0,w);
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
	//if (goodw0 && fitProb>0.2) {
	if (goodw0 && fitProb>fitProbRef) {
	  h2->Fill(ifwd, icnt, w);
	  p2m1->Fill(ifwd, icnt, recoWmass0 / mwpdg, w);
	  p2m2->Fill(ifwd, icnt, pow(recoWmass0 / mwpdg,2), w);
	}
      } // fitProb>0.2
   

      // Fit probability variations
      if (fabs(eta1)<etaMax && fabs(eta2)<etaMax) {
	for (int i = 1; i != pmm13bptboth_fp1->GetNbinsX()+1; ++i) {
	  double ptmin = pmm13bptboth_fp1->GetBinLowEdge(i);
	  double ptmax = pmm13bptboth_fp1->GetBinLowEdge(i+1);
	  double ptmid = 0.5*(ptmin+ptmax);
	  bool goodw0 = (recoWmass0>60 && recoWmass0<100);
	  if (pt1>ptmin && pt2>ptmin && pt1<ptmax && pt2<ptmax) {
	    if (fitProb>0.0 && goodw0) 
	      pmm13bptboth_fp0->Fill(ptmid,recoWmass0,w);
	    //
	    if (fitProb>0.01 && goodw0) 
	      pmm13bptboth_fp01->Fill(ptmid,recoWmass0,w);
	    if (fitProb>0.02 && goodw0) 
	      pmm13bptboth_fp02->Fill(ptmid,recoWmass0,w);
	    if (fitProb>0.05 && goodw0) 
	      pmm13bptboth_fp05->Fill(ptmid,recoWmass0,w);
	    //
	    if (fitProb>0.1 && goodw0) 
	      pmm13bptboth_fp1->Fill(ptmid,recoWmass0,w);
	    if (fitProb>0.2 && goodw0) 
	      pmm13bptboth_fp2->Fill(ptmid,recoWmass0,w);
	    if (fitProb>0.3 && goodw0) 
	      pmm13bptboth_fp3->Fill(ptmid,recoWmass0,w);
	    if (fitProb>0.4 && goodw0)
	      pmm13bptboth_fp4->Fill(ptmid,recoWmass0,w);
	    if (fitProb>0.5 && goodw0)
	      pmm13bptboth_fp5->Fill(ptmid,recoWmass0,w);
	    if (fitProb>0.6 && goodw0)
	      pmm13bptboth_fp6->Fill(ptmid,recoWmass0,w);
	    if (fitProb>0.7 && goodw0)
	      pmm13bptboth_fp7->Fill(ptmid,recoWmass0,w);
	    if (fitProb>0.8 && goodw0)
	      pmm13bptboth_fp8->Fill(ptmid,recoWmass0,w);
	    if (fitProb>0.9 && goodw0)
	      pmm13bptboth_fp9->Fill(ptmid,recoWmass0,w);
	  } // ptboth
	  //
	  //double ptave = 0.5*(pt1+pt2);
	  if (fitProb>0.0 && goodw0) 
	    pmm13bptave_fp0->Fill(ptave,recoWmass0,w);
	  //
	  if (fitProb>0.01 && goodw0) 
	    pmm13bptave_fp01->Fill(ptave,recoWmass0,w);
	  if (fitProb>0.02 && goodw0) 
	    pmm13bptave_fp02->Fill(ptave,recoWmass0,w);
	  if (fitProb>0.05 && goodw0) 
	    pmm13bptave_fp05->Fill(ptave,recoWmass0,w);
	  //
	  if (fitProb>0.1 && goodw0) 
	    pmm13bptave_fp1->Fill(ptave,recoWmass0,w);
	  if (fitProb>0.2 && goodw0) 
	    pmm13bptave_fp2->Fill(ptave,recoWmass0,w);
	  if (fitProb>0.3 && goodw0) 
	    pmm13bptave_fp3->Fill(ptave,recoWmass0,w);
	  if (fitProb>0.4 && goodw0)
	    pmm13bptave_fp4->Fill(ptave,recoWmass0,w);
	  if (fitProb>0.5 && goodw0)
	    pmm13bptave_fp5->Fill(ptave,recoWmass0,w);
	  if (fitProb>0.6 && goodw0)
	    pmm13bptave_fp6->Fill(ptave,recoWmass0,w);
	  if (fitProb>0.7 && goodw0)
	    pmm13bptave_fp7->Fill(ptave,recoWmass0,w);
	  if (fitProb>0.8 && goodw0)
	    pmm13bptave_fp8->Fill(ptave,recoWmass0,w);
	  if (fitProb>0.9 && goodw0)
	    pmm13bptave_fp9->Fill(ptave,recoWmass0,w);
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
   pm2m13bptave->Write("pm2m13bptave",TObject::kOverwrite);
   pmm13bptave_noL3Res->Write("pmm13bptave_noL3Res",TObject::kOverwrite);
   pmm13ptthr->Write("pmm13ptthr",TObject::kOverwrite);
   pmm13ptboth->Write("pmm13ptboth",TObject::kOverwrite);
   pmm13bptboth->Write("pmm13bptboth",TObject::kOverwrite);
   pmm13bptboth_noL3Res->Write("pmm13bptboth_noL3Res",TObject::kOverwrite);
   //
   pfuds->Write("pfuds",TObject::kOverwrite);
   pfud->Write("pfud",TObject::kOverwrite);
   pfu->Write("pfu",TObject::kOverwrite);
   pfd->Write("pfd",TObject::kOverwrite);
   pfs->Write("pfs",TObject::kOverwrite);
   pfc->Write("pfc",TObject::kOverwrite);
   pfb->Write("pfb",TObject::kOverwrite);
   pfg->Write("pfg",TObject::kOverwrite);
   pfo->Write("pfo",TObject::kOverwrite);
   //
   pmm13bptboth_fp0->Write("pmm13bptboth_fp0",TObject::kOverwrite);
   //
   pmm13bptboth_fp01->Write("pmm13bptboth_fp01",TObject::kOverwrite);
   pmm13bptboth_fp02->Write("pmm13bptboth_fp02",TObject::kOverwrite);
   pmm13bptboth_fp05->Write("pmm13bptboth_fp05",TObject::kOverwrite);
   //
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
   //
   pmm13bptave_fp01->Write("pmm13bptave_fp01",TObject::kOverwrite);
   pmm13bptave_fp02->Write("pmm13bptave_fp02",TObject::kOverwrite);
   pmm13bptave_fp05->Write("pmm13bptave_fp05",TObject::kOverwrite);
   //
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
   //
   pm1->Write("pm1",TObject::kOverwrite);
   pm2->Write("pm2",TObject::kOverwrite);
   pw1->Write("pw1",TObject::kOverwrite);
   pw2->Write("pw2",TObject::kOverwrite);
   pww1->Write("pww1",TObject::kOverwrite);
   pww2->Write("pww2",TObject::kOverwrite);
   p1->Write("p1",TObject::kOverwrite);
   p2->Write("p2",TObject::kOverwrite);
   pj1->Write("pj1",TObject::kOverwrite);
   pj2->Write("pj2",TObject::kOverwrite);
   pd1->Write("pd1",TObject::kOverwrite);
   pd2->Write("pd2",TObject::kOverwrite);
   pg1->Write("pg1",TObject::kOverwrite);
   pg2->Write("pg2",TObject::kOverwrite);
   //
   pmrg1->Write("pmrg1",TObject::kOverwrite);
   pmrg2->Write("pmrg2",TObject::kOverwrite);
   pmra1->Write("pmra1",TObject::kOverwrite);
   pmra2->Write("pmra2",TObject::kOverwrite);
   pjm1->Write("pjm1",TObject::kOverwrite);
   pjm2->Write("pjm2",TObject::kOverwrite);
   pgm1->Write("pgm1",TObject::kOverwrite);
   pgm2->Write("pgm2",TObject::kOverwrite);
   //
   pmjet->Write("pmjet",TObject::kOverwrite);
   pavecf->Write("pavecf",TObject::kOverwrite);
   pmreco->Write("pmreco",TObject::kOverwrite);
   pmgen->Write("pmgen",TObject::kOverwrite);
   pavegen->Write("pavegen",TObject::kOverwrite);
   pavereco->Write("pavereco",TObject::kOverwrite);
   pavegencf->Write("pavegencf",TObject::kOverwrite);
   paverecocf->Write("paverecocf",TObject::kOverwrite);

   hrbq->Write("hrbq",TObject::kOverwrite);
   h2bq->Write("h2bq",TObject::kOverwrite);
   h2rbqq->Write("h2rbqq",TObject::kOverwrite);
   prbqq->Write("prbqq",TObject::kOverwrite);
   h2rbqb->Write("h2rbqb",TObject::kOverwrite);
   prbqb->Write("prbqb",TObject::kOverwrite);
   h2rbqa->Write("h2rbqa",TObject::kOverwrite);
   prbqa->Write("prbqa",TObject::kOverwrite);
   h2rbqc->Write("h2rbqc",TObject::kOverwrite);
   prbqc->Write("prbqc",TObject::kOverwrite);
   pb1->Write("pb1",TObject::kOverwrite);
   pb2->Write("pb2",TObject::kOverwrite);
   pbq1->Write("pbq1",TObject::kOverwrite);
   pbq2->Write("pbq2",TObject::kOverwrite);
   pbb1->Write("pbb1",TObject::kOverwrite);
   pbb2->Write("pbb2",TObject::kOverwrite);
   pba1->Write("pba1",TObject::kOverwrite);
   pba2->Write("pba2",TObject::kOverwrite);
   //
   h1lep->Write("h1ptlep",TObject::kOverwrite);
   h1pt1->Write("h1ptt1",TObject::kOverwrite);
   h1mt->Write("h1mt",TObject::kOverwrite);
   h1mt1->Write("h1mt1",TObject::kOverwrite);
   h1mt2->Write("h1mt2",TObject::kOverwrite);
   h1mtmin->Write("h1mtmin",TObject::kOverwrite);
   h1mtmax->Write("h1mtmax",TObject::kOverwrite);
   h2mt->Write("h2mt",TObject::kOverwrite);
   h2mt1->Write("h2mt1",TObject::kOverwrite);
   h2mt2->Write("h2mt2",TObject::kOverwrite);
   h2mtmin->Write("h2mtmin",TObject::kOverwrite);
   h2mtmax->Write("h2mtmax",TObject::kOverwrite);
   pmwiov->Write("pmwiov",TObject::kOverwrite);
   pbqiov->Write("prbqiov",TObject::kOverwrite);
   plbiov->Write("pmlbiov",TObject::kOverwrite);
   ppt1iov->Write("pptt1iov",TObject::kOverwrite);
   plepiov->Write("plepiov",TObject::kOverwrite);
   pmtiov->Write("pmtiov",TObject::kOverwrite);
   pmt1iov->Write("pmt1iov",TObject::kOverwrite);
   pmt->Write("pmt",TObject::kOverwrite);
   pmtup->Write("pmtup",TObject::kOverwrite);
   pmtdw->Write("pmtdw",TObject::kOverwrite);
   pmta->Write("pmta",TObject::kOverwrite);
   pmt1->Write("pmt1",TObject::kOverwrite);
   pmt2->Write("pmt2",TObject::kOverwrite);
   pmtmin->Write("pmtmin",TObject::kOverwrite);
   pmtmax->Write("pmtmax",TObject::kOverwrite);
   //
   pmuiov->Write("pmuiov",TObject::kOverwrite);
   pnpviov->Write("pnpviov",TObject::kOverwrite);
   prhoiov->Write("prhoiov",TObject::kOverwrite);
   pmuvsmu->Write("pmuvsmu",TObject::kOverwrite);
   pnpvvsmu->Write("pnpvvsmu",TObject::kOverwrite);
   prhovsmu->Write("prhovsmu",TObject::kOverwrite);
   //
   h1rb->Write("h1rb",TObject::kOverwrite);
   h1gb->Write("h1gb",TObject::kOverwrite);
   h2rb->Write("h2rb",TObject::kOverwrite);
   h2gb->Write("h2gb",TObject::kOverwrite);
   prb->Write("prb",TObject::kOverwrite);
   prb0->Write("prb0",TObject::kOverwrite);
   pgb->Write("pgb",TObject::kOverwrite);
   pgb0->Write("pgb0",TObject::kOverwrite);
   //
   pmbjet->Write("pmbjet",TObject::kOverwrite);
   pmbreco->Write("pmbreco",TObject::kOverwrite);
   pmbgen->Write("pmbgen",TObject::kOverwrite);
   pbgen->Write("pbgen",TObject::kOverwrite);
   pbreco->Write("pbreco",TObject::kOverwrite);
   //
   hatmw->Write("atlas_mw",TObject::kOverwrite);
   hatbq->Write("atlas_rbq",TObject::kOverwrite);
   hatmt->Write("atlas_mt",TObject::kOverwrite);
   hatlb->Write("atlas_mlb",TObject::kOverwrite);
   //
   hatlbmet->Write("atlas_mlbmet",TObject::kOverwrite);
   //hatlbmetz->Write("atlas_mlbmetz",TObject::kOverwrite);
   pgap->Write("pgap",TObject::kOverwrite);
   hatpt5->Write("hatpt5",TObject::kOverwrite);
   hatpttt->Write("hatpttt",TObject::kOverwrite);
   hatpt5ott->Write("hatpt5ott",TObject::kOverwrite);
   h2atpt5vtt->Write("h2atpt5vtt",TObject::kOverwrite);
   patpt5ott->Write("patpt5ott",TObject::kOverwrite);
   //
   hatmwsw->Write("atlas_mw_signal",TObject::kOverwrite);
   hatmws1->Write("atlas_mw_signal_noweight",TObject::kOverwrite);
   hatbqsw->Write("atlas_rbq_signal",TObject::kOverwrite);
   hatbqs1->Write("atlas_rbq_signal_noweight",TObject::kOverwrite);
   //
   hatmwqq->Write("atlas_mw_qq",TObject::kOverwrite);
   hatmwqg->Write("atlas_mw_qg",TObject::kOverwrite);
   hatmwgg->Write("atlas_mw_gg",TObject::kOverwrite);
   patmwa->Write("patmwa",TObject::kOverwrite);
   patmwqqa->Write("patmwqqa",TObject::kOverwrite);
   patmwqga->Write("patmwqga",TObject::kOverwrite);
   patmwgga->Write("patmwgga",TObject::kOverwrite);
   patmwb->Write("patmwb",TObject::kOverwrite);
   patmwqb->Write("patmwqb",TObject::kOverwrite);
   patmwgb->Write("patmwgb",TObject::kOverwrite);
   //
   pqelfa->Write("pqelfa",TObject::kOverwrite);
   pqelfb->Write("pqelfb",TObject::kOverwrite);
   pqelfn->Write("pqelfn",TObject::kOverwrite);
   pqmufa->Write("pqmufa",TObject::kOverwrite);
   pqmufb->Write("pqmufb",TObject::kOverwrite);
   pqmufn->Write("pqmufn",TObject::kOverwrite);
   //
   pbelfa->Write("pbelfa",TObject::kOverwrite);
   pbelfb->Write("pbelfb",TObject::kOverwrite);
   pbelfn->Write("pbelfn",TObject::kOverwrite);
   pbmufa->Write("pbmufa",TObject::kOverwrite);
   pbmufb->Write("pbmufb",TObject::kOverwrite);
   pbmufn->Write("pbmufn",TObject::kOverwrite);

   hqglq->Write("hqglq",TObject::kOverwrite);
   hqglb->Write("hqglb",TObject::kOverwrite);
   h2qglqa->Write("h2qglqa",TObject::kOverwrite);
   h2qglqb->Write("h2qglqb",TObject::kOverwrite);
   h2qglqb_q->Write("h2qglqb_q",TObject::kOverwrite);
   h2qglqb_g->Write("h2qglqb_g",TObject::kOverwrite);
   h2qglqb_u->Write("h2qglqb_u",TObject::kOverwrite);
   h2qglqb_d->Write("h2qglqb_d",TObject::kOverwrite);
   h2qglqb_s->Write("h2qglqb_s",TObject::kOverwrite);
   h2qglqb_c->Write("h2qglqb_c",TObject::kOverwrite);
   h2qglqb_b->Write("h2qglqb_b",TObject::kOverwrite);
   h2qglqb_o->Write("h2qglqb_o",TObject::kOverwrite);
   h2qglba->Write("h2qglba",TObject::kOverwrite);
   h2qglbb->Write("h2qglbb",TObject::kOverwrite);
   h2qglbb_b->Write("h2qglbb_b",TObject::kOverwrite);

   hqglg->Write("hqglg",TObject::kOverwrite);
   h2qglg->Write("h2qglg",TObject::kOverwrite);
   h2qglg_q->Write("h2qglg_q",TObject::kOverwrite);
   h2qglg_u->Write("h2qglg_u",TObject::kOverwrite);
   h2qglg_d->Write("h2qglg_d",TObject::kOverwrite);
   h2qglg_s->Write("h2qglg_s",TObject::kOverwrite);
   h2qglg_c->Write("h2qglg_c",TObject::kOverwrite);
   h2qglg_b->Write("h2qglg_b",TObject::kOverwrite);
   h2qglg_g->Write("h2qglg_g",TObject::kOverwrite);
   h2qglg_o->Write("h2qglg_o",TObject::kOverwrite);

   peqq->Write("peqq",TObject::kOverwrite);
   pequ->Write("pequ",TObject::kOverwrite);
   peqd->Write("peqd",TObject::kOverwrite);
   peqs->Write("peqs",TObject::kOverwrite);
   peqc->Write("peqc",TObject::kOverwrite);
   peqb->Write("peqb",TObject::kOverwrite);
   peqg->Write("peqg",TObject::kOverwrite);
   peqo->Write("peqo",TObject::kOverwrite);

   pebb->Write("pebb",TObject::kOverwrite);
   pebq->Write("pebq",TObject::kOverwrite);
   pebg->Write("pebg",TObject::kOverwrite);
   pebo->Write("pebo",TObject::kOverwrite);

   pegg->Write("pegg",TObject::kOverwrite);
   pegq->Write("pegq",TObject::kOverwrite);
   pegu->Write("pegu",TObject::kOverwrite);
   pegd->Write("pegd",TObject::kOverwrite);
   pegs->Write("pegs",TObject::kOverwrite);
   pegc->Write("pegc",TObject::kOverwrite);
   pegb->Write("pegb",TObject::kOverwrite);
   pego->Write("pego",TObject::kOverwrite);

   pfqq->Write("pfqq",TObject::kOverwrite);
   pfqu->Write("pfqu",TObject::kOverwrite);
   pfqd->Write("pfqd",TObject::kOverwrite);
   pfqs->Write("pfqs",TObject::kOverwrite);
   pfqc->Write("pfqc",TObject::kOverwrite);
   pfqb->Write("pfqb",TObject::kOverwrite);
   pfqg->Write("pfqg",TObject::kOverwrite);
   pfqo->Write("pfqo",TObject::kOverwrite);

   pfbb->Write("pfbb",TObject::kOverwrite);
   pfbq->Write("pfbq",TObject::kOverwrite);
   pfbg->Write("pfbg",TObject::kOverwrite);
   pfbo->Write("pfbo",TObject::kOverwrite);

   pfgg->Write("pfgg",TObject::kOverwrite);
   pfgq->Write("pfgq",TObject::kOverwrite);
   pfgu->Write("pfgu",TObject::kOverwrite);
   pfgd->Write("pfgd",TObject::kOverwrite);
   pfgs->Write("pfgs",TObject::kOverwrite);
   pfgc->Write("pfgc",TObject::kOverwrite);
   pfgb->Write("pfgb",TObject::kOverwrite);
   pfgo->Write("pfgo",TObject::kOverwrite);

   h2LHGO->Write("h2LHGO",TObject::kOverwrite);
   h2L->Write("h2L",TObject::kOverwrite);
   h2H->Write("h2H",TObject::kOverwrite);
   h2G->Write("h2G",TObject::kOverwrite);
   h2O->Write("h2O",TObject::kOverwrite);

   hmABCD->Write("hmABCD",TObject::kOverwrite);
   hmA->Write("hmA",TObject::kOverwrite);
   hmB->Write("hmB",TObject::kOverwrite);
   hmC->Write("hmC",TObject::kOverwrite);
   hmD->Write("hmD",TObject::kOverwrite);

   hmLHGO->Write("hmLHGO",TObject::kOverwrite);
   hmL->Write("hmL",TObject::kOverwrite);
   hmH->Write("hmH",TObject::kOverwrite);
   hmG->Write("hmG",TObject::kOverwrite);
   hmO->Write("hmO",TObject::kOverwrite);

   hnABCD->Write("hnABCD",TObject::kOverwrite);
   hnA->Write("hnA",TObject::kOverwrite);
   hnB->Write("hnB",TObject::kOverwrite);
   hnC->Write("hnC",TObject::kOverwrite);
   hnD->Write("hnD",TObject::kOverwrite);

   hnLHGO->Write("hnLHGO",TObject::kOverwrite);
   hnL->Write("hnL",TObject::kOverwrite);
   hnH->Write("hnH",TObject::kOverwrite);
   hnG->Write("hnG",TObject::kOverwrite);
   hnO->Write("hnO",TObject::kOverwrite);
   
   h2n->Write("h2n",TObject::kOverwrite);
   p2m->Write("p2m",TObject::kOverwrite);

   hl->Write("hl",TObject::kOverwrite);
   hnu->Write("hnu",TObject::kOverwrite);
   hlnu->Write("hlnu",TObject::kOverwrite);
   hu->Write("hu",TObject::kOverwrite);
   hd->Write("hd",TObject::kOverwrite);
   hud->Write("hud",TObject::kOverwrite);
   hs->Write("hs",TObject::kOverwrite);
   hc->Write("hc",TObject::kOverwrite);
   hcs->Write("hcs",TObject::kOverwrite);
   plnu->Write("plnu",TObject::kOverwrite);
   pud->Write("pud",TObject::kOverwrite);
   pcs->Write("pcs",TObject::kOverwrite);

   psmw->Write("psmw",TObject::kOverwrite);
   psmt->Write("psmt",TObject::kOverwrite);
   psmwb->Write("psmwb",TObject::kOverwrite);
   psmwB->Write("psmwbk",TObject::kOverwrite);
   psmlb->Write("psmlb",TObject::kOverwrite);
   psrbq->Write("psrbq",TObject::kOverwrite);
   psrbb->Write("psrbb",TObject::kOverwrite);
   psrqq->Write("psrqq",TObject::kOverwrite);
   psptb->Write("psptb",TObject::kOverwrite);
   psptB->Write("psptbk",TObject::kOverwrite);
   psptqq->Write("psptqq",TObject::kOverwrite);
   psptw->Write("psptw",TObject::kOverwrite);
   psptt->Write("psptt",TObject::kOverwrite);
   psrb->Write("psrb",TObject::kOverwrite);
   psrq->Write("psrq",TObject::kOverwrite);
   //
   psgmw->Write("psgmw",TObject::kOverwrite);
   psgmt->Write("psgmt",TObject::kOverwrite);
   psgmwb->Write("psgmwb",TObject::kOverwrite);
   psgmwB->Write("psgmwbk",TObject::kOverwrite);
   psgmlb->Write("psgmlb",TObject::kOverwrite);
   psgrbq->Write("psgrbq",TObject::kOverwrite);
   psgrbb->Write("psgrbb",TObject::kOverwrite);
   psgrqq->Write("psgrqq",TObject::kOverwrite);
   psgptb->Write("psgptb",TObject::kOverwrite);
   psgptB->Write("psgptbk",TObject::kOverwrite);
   psgptqq->Write("psgptqq",TObject::kOverwrite);
   psgptw->Write("psgptw",TObject::kOverwrite);
   psgptt->Write("psgptt",TObject::kOverwrite);

   htopptorig->Write("htopptorig",TObject::kOverwrite);
   htopptgen->Write("htopptgen",TObject::kOverwrite);
   htoppt->Write("htoppt",TObject::kOverwrite);
   htopptw0->Write("htopptw0",TObject::kOverwrite);
   htoppt_fsr->Write("htoppt_fsr",TObject::kOverwrite);
   htoppt_q2qg->Write("htoppt_q2qg",TObject::kOverwrite);
   htoppt_x2xg->Write("htoppt_x2xg",TObject::kOverwrite);
   htoppt_isr->Write("htoppt_isr",TObject::kOverwrite);

   pmpf->Write("pmpf",TObject::kOverwrite);

   fout->Close();
   curdir->cd();

   cout << "Processed " << igood << " good events, "
	<< ibad << " vetoed ones" << endl;
   cout << "Veto inefficiency " << 100.*ibad/max(igood+ibad,1) << "%" << endl;
   cout << "Removed " << cnt << " duplicates" << endl;
} // Loop

void hadW::Draw(string s) {//string sdt, string smc) {

  const char *c = s.c_str();
   bool is18 = (TString(s.c_str()).Contains("18"));
   bool is17 = (TString(s.c_str()).Contains("17"));
   bool is16 = (TString(s.c_str()).Contains("16"));
   //assert((is18 || is17) && (!is18 || !is17)); // XOR
   assert(is18 || is17 || is16);
   if ( is18) cout << "Drawing UL18" << endl;
   if ( is17) cout << "Drawing UL17" << endl;
   if ( is16) cout << "Drawing UL16" << endl;

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  //TFile *fdt = new TFile("rootfiles/hadWUL17.root","READ");
  //TFile *fdt = new TFile("rootfiles/hadWUL18.root","READ");
  //TFile *fdt = new TFile("rootfiles/hadWUL18V4_MPDGcorrNoW.root","READ");
  //TFile *fdt = new TFile("rootfiles/hadWUL17V5_MPDGcorrNoW.root","READ");
  //TFile *fdt = new TFile(Form("rootfiles/hadWUL%s_MPDGcorrNoW.root",c),"READ");
  //TFile *fdt = new TFile(Form("rootfiles/hadWUL%s_EMUF.root",c),"READ");
  TFile *fdt = new TFile(Form("rootfiles/hadWUL%s_Glu.root",c),"READ");
  assert(fdt && !fdt->IsZombie());
  //TFile *fmc = new TFile("rootfiles/hadWMC17.root","READ");
  //TFile *fmc = new TFile("rootfiles/hadWMC18.root","READ");
  //TFile *fmc = new TFile("rootfiles/hadWMC18V4_MPDGcorrNoW.root","READ");
  //TFile *fmc = new TFile("rootfiles/hadWMC17V5_MPDGcorrNoW.root","READ");
  //TFile *fmc = new TFile(Form("rootfiles/hadWMC%s_MPDGcorrNoW.root",c),"READ");
  //TFile *fmc = new TFile(Form("rootfiles/hadWMC%s_EMUF.root",c),"READ");
  TFile *fmc = new TFile(Form("rootfiles/hadWMC%s_Glu.root",c),"READ");
  assert(fdt && !fdt->IsZombie());
  //TFile *fout = new TFile("rootfiles/hadWUL17.root","RECREATE");
  //TFile *fout = new TFile("rootfiles/hadWUL18.root","UPDATE");
  //TFile *fout = new TFile("rootfiles/hadWUL18V4_MPDGcorrNoW.root","UPDATE");
  //TFile *fout = new TFile("rootfiles/hadWUL17V5_MPDGcorrNoW.root","UPDATE");
  //TFile *fout = new TFile(Form("rootfiles/hadW%s_MPDGcorrNoW.root",c),"UPDATE");
  //TFile *fout = new TFile(Form("rootfiles/hadW%s_MPDGcorrNoW.root",c),"RECREATE");
  //TFile *fout = new TFile(Form("rootfiles/hadW%s_EMUF.root",c),"RECREATE");
  TFile *fout = new TFile(Form("rootfiles/hadW%s_Glu.root",c),"RECREATE");
  curdir->cd();

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

  TH1D *hup =tdrHist("hup","#LTm_{W}#GT (GeV)",73,90,"p_{T,jets} (GeV)",30,200);
  TH1D *hdw =tdrHist("hdw","Data/MC-1 (%)",-1.5,+1.5,"p_{T,jets} (GeV)",30,200);

  if (is18) lumi_13TeV = "UL18 ttbar lepton+jet hadW";
  if (is17) lumi_13TeV = "UL17 ttbar lepton+jet hadW";
  if (is16) lumi_13TeV = "UL16 ttbar lepton+jet hadW";
  if (is17&&is18) lumi_13TeV = "UL17+18 ttbar lepton+jet hadW";
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

  TLatex *tex = new TLatex();
  tex->SetTextSize(0.045); tex->SetNDC();
  tex->DrawLatex(0.72,0.855,Form("|#eta_{j}| < %1.1f",etaMax));
  tex->DrawLatex(0.19,0.71,Form("fitProb>%1.2g",fitProbRef));

  TLegend *leg1 = tdrLeg(0.39,0.64,0.59,0.89);
  leg1->SetHeader("DT");
  leg1->AddEntry(ptdt," ","PLE");
  leg1->AddEntry(pbdt," ","PLE");
  leg1->AddEntry(pt0dt," ","PLE");
  leg1->AddEntry(pb0dt," ","PLE");

  TLegend *leg2 = tdrLeg(0.45,0.64,0.65,0.89);
  leg2->SetHeader("MC");
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

  //c1->SaveAs("pdf/hadW_MPDGcorrNoW.pdf");
  //c1->SaveAs("pdf/hadW_UL17V5_MPDGcorrNoW.pdf");
  c1->SaveAs(Form("pdf/hadW_%s_MPDGcorrNoW.pdf",c));

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

  TH1D *hup2 = tdrHist("hup2","#LTm_{W}#GT (GeV)",75.0+1e-5,90.0-1e-5,
		       "p_{T,jets} (GeV)",30,200);
  //TH1D *hdw2 = tdrHist("hdw2","Data/MC-1 (%)",-1.2,+1.1,
  TH1D *hdw2 = tdrHist("hdw2","Data/MC-1 (%)",-3.5,+2.0,
		       "p_{T,jets} (GeV)",30,200);

  TCanvas *c2 = tdrDiCanvas("c2",hup2,hdw2,4,11);
  
  c2->cd(1);
  gPad->SetLogx();

  l->DrawLine(30,80.4,200,80.4);
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

  tex->DrawLatex(0.40,0.84,Form("fitProb>%1.2g",fitProbRef));
  tex->DrawLatex(0.40,0.78,Form("|#eta_{j}| < %1.1f",etaMax));
  tex->DrawLatex(0.30,0.72,"60<m_{W}<100 GeV");

  TLegend *leg12 = tdrLeg(0.62,0.72,0.82,0.87);
  leg12->SetHeader("DT");
  leg12->AddEntry(pb0dt," ","PLE");
  leg12->AddEntry(pa0dt," ","PLE");

  TLegend *leg22 = tdrLeg(0.68,0.72,0.88,0.87);
  leg22->SetHeader("MC");
  leg22->AddEntry(pb0mc," vs p_{T,both}","PLE");
  leg22->AddEntry(pa0mc," vs p_{T,ave}","PLE");

  c2->cd(2);
  gPad->SetLogx();

  l->DrawLine(30,0,200,0);
  //tdrDraw(pt0r,"Pz",kFullCircle,kMagenta+2);
  tdrDraw(pb0r,"Pz",kFullCircle,kCyan+2);
  tdrDraw(gb0r,"Pz",kFullCircle,kCyan+3);
  // Check that L3Res was correctly removed
  //tdrDraw(gr,"Pz",kFullCircle,kBlack);
  tdrDraw(pa0r,"Pz",kFullDiamond,kMagenta+1);
  tdrDraw(ga0r,"Pz",kFullDiamond,kMagenta+2);
  pa0r->SetMarkerSize(1.5);

  gb0r->GetXaxis()->SetTitle("p_{T,jets} (GeV)");
  gr->GetXaxis()->SetTitle("p_{T,jets} (GeV)");
  //
  gb0r->GetYaxis()->SetTitle("m_{W,0,data}/m_{W,0,MC}-1 (%)");
  gr->GetYaxis()->SetTitle("m_{W,0,data}/m_{W,0,MC}");

  //c2->SaveAs("pdf/hadW0_MPDGcorrNoW.pdf");
  //c2->SaveAs("pdf/hadW0_UL17V5_MPDGcorrNoW.pdf");
  c2->SaveAs(Form("pdf/hadW0_%s_MPDGcorrNoW.pdf",c));
  
  fout->cd();
  if (fabs(fitProbRef-0.2)<1e-3) {
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
  }
  if (fabs(fitProbRef-0.01)<1e-3) {
  hdt->Write("data_nevents_ptboth_hadw_fitprob001_L1L2L3");
  hmc->Write("mc_nevents_ptboth_hadw_fitprob001_L1L2L3");
  hr->Write("ratio_nevents_ptboth_hadw_fitprob001_L1L2L3");
  gb0dt->Write("data_mass_ptboth_hadw_fitprob001_L1L2L3Res",TObject::kOverwrite);
  gb0mc->Write("mc_mass_ptboth_hadw_fitprob001_L1L2L3Res",TObject::kOverwrite);
  gb0r->Write("ratio_mass_ptboth_hadw_fitprob001_L1L2L3Res",TObject::kOverwrite);
  gdt->Write("data_mass_ptboth_hadw_fitprob001_L1L2L3",TObject::kOverwrite);
  gmc->Write("mc_mass_ptboth_hadw_fitprob001_L1L2L3",TObject::kOverwrite);
  gr->Write("ratio_mass_ptboth_hadw_fitprob001_L1L2L3",TObject::kOverwrite);
  //
  hadt->Write("data_nevents_ptave_hadw_fitprob001_L1L2L3");
  hamc->Write("mc_nevents_ptave_hadw_fitprob001_L1L2L3");
  har->Write("ratio_nevents_ptave_hadw_fitprob001_L1L2L3");
  ga0dt->Write("data_mass_ptave_hadw_fitprob001_L1L2L3Res",TObject::kOverwrite);
  ga0mc->Write("mc_mass_ptave_hadw_fitprob001_L1L2L3Res",TObject::kOverwrite);
  ga0r->Write("ratio_mass_ptave_hadw_fitprob001_L1L2L3Res",TObject::kOverwrite);
  gadt->Write("data_mass_ptave_hadw_fitprob001_L1L2L3",TObject::kOverwrite);
  gamc->Write("mc_mass_ptave_hadw_fitprob001_L1L2L3",TObject::kOverwrite);
  gar->Write("ratio_mass_ptave_hadw_fitprob001_L1L2L3",TObject::kOverwrite);
  }
  fout->Close();

  // Draw m=0 case only as well? Proper <pT,j> centering
  // Add flavor corrections. Ad-hoc mass correction for PU as well?
  // Matrix case (TH2D of #jets and their mW in neighbouring bins)?

} // draw

// Draw a scan of fitProbability cuts
void hadW::DrawFP(string spt, string mode) {

  const char *cm = mode.c_str();

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  //TFile *fdt = new TFile("rootfiles/hadWUL17.root","READ");
  //TFile *fdt = new TFile("rootfiles/hadWUL18.root","READ");
  //TFile *fdt = new TFile("rootfiles/hadWUL18V4_MPDGcorrNoW.root","READ");
  //TFile *fdt = new TFile("rootfiles/hadWUL17V5_MPDGcorrNoW.root","READ");
  //TFile *fdt = new TFile("rootfiles/hadWUL17V5_MPDGcorrNoW.root","READ");
  //TFile *fdt = new TFile(Form("rootfiles/hadWUL%s_MPDGcorrNoW.root",cm),"READ");
  //TFile *fdt = new TFile(Form("rootfiles/hadWUL%s_EMUF.root",cm),"READ");
  TFile *fdt = new TFile(Form("rootfiles/hadWUL%s_Glu.root",cm),"READ");
  assert(fdt && !fdt->IsZombie());
  //TFile *fmc = new TFile("rootfiles/hadWMC17.root","READ");
  //TFile *fmc = new TFile("rootfiles/hadWMC18.root","READ");
  //TFile *fmc = new TFile("rootfiles/hadWMC18V4_MPDGcorrNoW.root","READ");
  //TFile *fmc = new TFile("rootfiles/hadWMC17V5_MPDGcorrNoW.root","READ");
  //TFile *fmc = new TFile(Form("rootfiles/hadWMC%s_MPDGcorrNoW.root",cm),"READ");
  //TFile *fmc = new TFile(Form("rootfiles/hadWMC%s_EMUF.root",cm),"READ");
  TFile *fmc = new TFile(Form("rootfiles/hadWMC%s_Glu.root",cm),"READ");
  assert(fdt && !fdt->IsZombie());
  curdir->cd();

  const char *cpt = spt.c_str();

  const int nfp = 11;
  TProfile *pdt[nfp], *pmc[nfp];
  TH1D *hr[nfp];
  for (int ifp = 0; ifp != nfp; ++ifp) {
    TProfile *pd = (TProfile*)fdt->Get(Form("pmm13b%s_fp%d",cpt,ifp));
    if (ifp==10)
      pd = (TProfile*)fdt->Get(Form("pmm13b%s_fp%02d",cpt,1));
    assert(pd);
    TProfile *pm = (TProfile*)fmc->Get(Form("pmm13b%s_fp%d",cpt,ifp));
    if (ifp==10)
      pm = (TProfile*)fmc->Get(Form("pmm13b%s_fp%02d",cpt,1));
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

  int color[nfp] = {kBlue, kCyan+1, /*kBlack,*/kGreen, kGreen+2, kYellow+1,
		    kYellow+2, kOrange+1, kOrange+3, kRed, kRed+2, kBlack};

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  TLatex *tex = new TLatex();
  tex->SetTextSize(0.045); tex->SetNDC();

  TH1D *hup = tdrHist("hupfp","#LTm_{W,0}#GT (GeV)",75.1+1e-5,90.1-1e-5,
		      "p_{T,jets} (GeV)",30,200);
  double djec = 0.3; // New L3Res
  //TH1D *hdw = tdrHist("hdwfp","Data/MC-1 (%)",-1.9+djec,+0.9+djec,
  TH1D *hdw = tdrHist("hdwfp","Data/MC-1 (%)",-3.5,2.0,
		      "p_{T,jets} (GeV)",30,200);
  if (spt=="ptboth") hdw->SetXTitle("p_{T,jets} (GeV)");
  if (spt=="ptave")  hdw->SetXTitle("p_{T,ave} (GeV)");

  TCanvas *c1 = tdrDiCanvas("c1fp",hup,hdw,4,11);
  
  c1->cd(1);
  gPad->SetLogx();

  l->DrawLine(30,80.4,200,80.4);
  tex->DrawLatex(0.40,0.84,Form("|#eta_{jets}| < %1.1f",etaMax));
  tex->DrawLatex(0.40,0.78,"60<m_{W}<100 GeV");

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

  tdrDraw(hsys,"E3",kNone,kYellow,kSolid,kYellow+1,1001,kYellow);
  l->DrawLine(30,0,200,0);
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

  c1->SaveAs(Form("pdf/hadW_drawFP_%s_MPDGcorrNoW_%s.pdf",cpt,cm));

} // drawFP


void hadW::DrawRMS(string s) {

  const char *c = s.c_str();

  // m = sqrt(pT1*pT2)*c(angular terms)
  // => sigma_m/m = 1/2*sqrt(sigma_pT1^2/pT1^2 + sigma_pT2^2/pT2^2) by error propagation
  
  setTDRStyle();
  TDirectory *curdir = gDirectory;
  //TFile *fdt = new TFile("rootfiles/hadWUL18V4_MPDGcorrNoW.root","READ");
  //TFile *fdt = new TFile("rootfiles/hadWUL17V5_MPDGcorrNoW.root","READ");
  //TFile *fdt = new TFile(Form("rootfiles/hadWUL%s_MPDGcorrNoW.root",c),"READ");
  //TFile *fdt = new TFile(Form("rootfiles/hadWUL%s_EMUF.root",c),"READ");
  TFile *fdt = new TFile(Form("rootfiles/hadWUL%s_Glu.root",c),"READ");
  assert(fdt && !fdt->IsZombie());
  //TFile *fmc = new TFile("rootfiles/hadWMC18V4_MPDGcorrNoW.root","READ");
  //TFile *fmc = new TFile("rootfiles/hadWMC17V5_MPDGcorrNoW.root","READ");
  //TFile *fmc = new TFile(Form("rootfiles/hadWMC%s_MPDGcorrNoW.root",c),"READ");
  //TFile *fmc = new TFile(Form("rootfiles/hadWMC%s_EMUF.root",c),"READ");
  TFile *fmc = new TFile(Form("rootfiles/hadWMC%s_Glu.root",c),"READ");
  assert(fdt && !fdt->IsZombie());
  curdir->cd();

  TProfile *pm1dt = (TProfile*)fdt->Get("pmm13bptave"); assert(pm1dt);
  TProfile *pm2dt = (TProfile*)fdt->Get("pm2m13bptave"); assert(pm2dt);
  TProfile *pm1mc = (TProfile*)fmc->Get("pmm13bptave"); assert(pm1mc);
  TProfile *pm2mc = (TProfile*)fmc->Get("pm2m13bptave"); assert(pm2mc);
  TProfile *pm1g = (TProfile*)fmc->Get("pm1"); assert(pm1g);
  TProfile *pm2g = (TProfile*)fmc->Get("pm2"); assert(pm2g);
  TProfile *pw1 = (TProfile*)fmc->Get("pw1"); assert(pw1);
  TProfile *pw2 = (TProfile*)fmc->Get("pw2"); assert(pw2);
  TProfile *pww1 = (TProfile*)fmc->Get("pww1"); assert(pww1);
  TProfile *pww2 = (TProfile*)fmc->Get("pww2"); assert(pww2);
  TProfile *pj1 = (TProfile*)fmc->Get("pj1"); assert(pj1);
  TProfile *pj2 = (TProfile*)fmc->Get("pj2"); assert(pj2);


  TH1D *hrmsdt = pm1dt->ProjectionX("hrmsdt");
  TH1D *hrmsmc = pm1mc->ProjectionX("hrmsmc");
  TH1D *hrmsg = pm1g->ProjectionX("hrmsg");
  TH1D *hrmsw = pw1->ProjectionX("hrmsw");
  TH1D *hrmsww = pww1->ProjectionX("hrmsww");
  TH1D *hrmsj = pj1->ProjectionX("hrmsj");
  TH1D *hrms = pm1dt->ProjectionX("hrms");
  for (int i = 1; i != hrmsdt->GetNbinsX()+1; ++i) {

    hrmsdt->SetBinContent(i, 100.*sqrt(max(pm2dt->GetBinContent(i) -
					   pow(pm1dt->GetBinContent(i),2),0.))
			  / 80.4);
    hrmsmc->SetBinContent(i, 100.*sqrt(max(pm2mc->GetBinContent(i) -
					   pow(pm1mc->GetBinContent(i),2),0.))
			  / 80.4);
    hrmsg->SetBinContent(i, 100.*sqrt(max(pm2g->GetBinContent(i) -
					  pow(pm1g->GetBinContent(i),2),0.))
			 / 80.4);
    hrmsw->SetBinContent(i, 100.*sqrt(max(pw2->GetBinContent(i) -
					  pow(pw1->GetBinContent(i),2),0.))
			 / 80.4);
    hrmsww->SetBinContent(i, 100.*sqrt(max(pww2->GetBinContent(i) -
					  pow(pww1->GetBinContent(i),2),0.))
			 / 80.4);
    hrmsj->SetBinContent(i, 100.*sqrt(max(pj2->GetBinContent(i) -
					  pow(pj1->GetBinContent(i),2),0.)));

    hrmsdt->SetBinError(i, hrmsdt->GetBinContent(i) *
			pm1dt->GetBinError(i)
			/ (hrmsdt->GetBinContent(i) * 80.4 / 100.));
    hrmsmc->SetBinError(i, hrmsmc->GetBinContent(i) *
			pm1mc->GetBinError(i)
			/ (hrmsmc->GetBinContent(i) * 80.4 / 100.));
    hrmsg->SetBinError(i, hrmsg->GetBinContent(i) *
		       pm1g->GetBinError(i)
		       / (hrmsg->GetBinContent(i) * 80.4 / 100.));
    hrmsw->SetBinError(i, hrmsw->GetBinContent(i) *
		       pw1->GetBinError(i)
		       / (hrmsw->GetBinContent(i) * 80.4 / 100.));
    hrmsww->SetBinError(i, hrmsww->GetBinContent(i) *
			pww1->GetBinError(i)
			/ (hrmsww->GetBinContent(i) * 80.4 / 100.));
    hrmsj->SetBinError(i, hrmsj->GetBinContent(i) *
		       pj1->GetBinError(i)
		       / (hrmsj->GetBinContent(i) / 100.));

    hrms->SetBinContent(i, 100.*(hrmsdt->GetBinContent(i) /
				 hrmsmc->GetBinContent(i) - 1));

    hrms->SetBinError(i, 100.*sqrt(pow(hrmsdt->GetBinError(i)
				       / hrmsdt->GetBinContent(i),2) +
				   pow(hrmsdt->GetBinError(i)
				       / hrmsdt->GetBinContent(i),2)));
  } // for i
  
  // M^2 = 2*pT1*pT2*(cosh(eta1-eta1) - cos(phi1-phi2))
  // pTave=(pT1+pT2)/2, pTdiff=(pT1-pT2)/2
  // pT1=pTave+pTdiff, pT2=pTave-pTdiff
  // pT1*pT2 = pTave^2*(1 + (pTdiff/pTave)^2)
  // => M ~ 2*pTave*sqrt(1+(pTdiff/pTave)^2)*sqrt((Cdeta-Cdphi)/2)

  TH1D *hup = tdrHist("huprms","RMS(m_{W,0} / 80.4 GeV) (%)",5.0,10.5,
		      //4.5,8.5,
		      "p_{T,ave} (GeV)",30,200);
  TH1D *hdw = tdrHist("hdwrms","Data/MC-1 (%)",-1,5,
		      "p_{T,ave} (GeV)",30,200);

  TCanvas *c1 = tdrDiCanvas("c1rms",hup,hdw,4,11);

  c1->cd(1);
  gPad->SetLogx();
  hrmsj->Scale(0.5);
  tdrDraw(hrmsj,"Pz",kOpenDiamond,kBlack);
  tdrDraw(hrmsw,"Pz",kFullDiamond,kGray+1);
  tdrDraw(hrmsww,"Pz",kFullDiamond,kBlack);
  tdrDraw(hrmsg,"Pz",kOpenSquare,kBlack);
  tdrDraw(hrmsmc,"Pz",kOpenCircle,kBlack);
  tdrDraw(hrmsdt,"Pz",kFullCircle,kBlack);
  
  TLegend *leg = tdrLeg(0.5,0.49,0.8,0.85);
  leg->AddEntry(hrmsdt,"Data","PLE");
  leg->AddEntry(hrmsmc,"MC","PLE");
  leg->AddEntry(hrmsg,"MC (gen-match #times 2)","PLE");
  leg->AddEntry(hrmsw,"Gen","PLE");
  leg->AddEntry(hrmsww,"Gen (80-100 GeV)","PLE");
  leg->AddEntry(hrmsj,"JER #times 0.5","PLE");

  c1->cd(2);
  gPad->SetLogx();
  tdrDraw(hrms,"Pz",kFullCircle,kBlack);

  //c1->SaveAs("pdf/hadw_drawRMS_ptave_MPDGcorrNoW.pdf");
  //c1->SaveAs("pdf/hadw_drawRMS_ptave_UL17V5_MPDGcorrNoW.pdf");
  c1->SaveAs(Form("pdf/hadw_drawRMS_ptave_%s_MPDGcorrNoW.pdf",c));
} // DrawRMS


// Analyze and draw 2D response
void hadW::Draw2D(string set) {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  //TFile *f = new TFile("rootfiles/hadWMC17_v2.root","READ");
  //TFile *f = new TFile(Form("rootfiles/hadW%s_v2.root",set.c_str()),"READ");
  //TFile *f = new TFile(Form("rootfiles/hadW%s.root",set.c_str()),"READ");
  //TFile *f = new TFile(Form("rootfiles/hadW%s_EMUF2.root",set.c_str()),"READ");
  TFile *f = new TFile(Form("rootfiles/hadW%s_Glu2.root",set.c_str()),"READ");
  assert(f && !f->IsZombie());

  //TFile *fout = new TFile("rootfiles/hadW.root","UPDATE");
  //TFile *fout = new TFile("rootfiles/hadWUL18.root","UPDATE");
  //TFile *fout = new TFile("rootfiles/hadWUL18V4_MPDGcorrNoW.root","UPDATE");
  //TFile *fout = new TFile("rootfiles/hadWUL17V5_MPDGcorrNoW.root","UPDATE");
  //TFile *fout = new TFile(Form("rootfiles/hadW%s_MPDGcorrNoW.root",set.c_str()),
  //TFile *fout = new TFile(Form("rootfiles/hadW%s_EMUF.root",set.c_str()),
  TFile *fout = new TFile(Form("rootfiles/hadW%s_Glu.root",set.c_str()),
			  "UPDATE");
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
    
    c7->SaveAs(Form("pdf/hadW_draw2D_2D_%s_UL17V5.pdf","ratio"));

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
