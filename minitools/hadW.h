//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jun 26 15:44:54 2020 by ROOT version 6.18/04
// from TChain tree/
//////////////////////////////////////////////////////////

// Updated with:
//////////////////////////////////////////////////////////                      
// This class has been automatically generated on                               
// Tue Dec  1 22:45:54 2020 by ROOT version 6.18/04                             
// from TTree tree/tree                                                         
// found on file: WmassUL18.root                                                
//////////////////////////////////////////////////////////   

#ifndef hadW_h
#define hadW_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <string>
#include <iostream>

using namespace std;

// Header file for the classes stored in the TTree if any.

class hadW {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   Float_t         pt1;
   Float_t         pt2;
   Float_t         eta1;
   Float_t         eta2;
   Float_t         phi1;
   Float_t         phi2;
   Float_t         m1;
   Float_t         m2;
   Float_t         sc1;
   Float_t         sc2;
   Float_t         muf1;
   Float_t         muf2;
   Float_t         elf1;
   Float_t         elf2;
   Float_t         nhf1;
   Float_t         nhf2;
   Float_t         qgl1;
   Float_t         qgl2;
   Float_t         btag1;
   Float_t         btag2;
   Float_t         bleptag1;
   Float_t         bleptag2;
   Float_t         ctag1;
   Float_t         ctag2;
   Float_t         udstag1;
   Float_t         udstag2;
   Float_t         gtag1;
   Float_t         gtag2;
   Int_t           flav1;
   Int_t           flav2;
   /*
   Int_t           ntrk1;
   Int_t           ntrk2;
   Int_t           nputrk1;
   Int_t           nputrk2;
   Float_t         recoWEta;
   Float_t         recoWPhi;
   */
   Float_t         recoWMass;
   Float_t         gen_pt1;
   Float_t         gen_pt2;
   Float_t         gen_eta1;
   Float_t         gen_eta2;
   Float_t         gen_phi1;
   Float_t         gen_phi2;
   Float_t         gen_m1;
   Float_t         gen_m2;
   Float_t         bpt1;
   Float_t         bpt2;
   Float_t         beta1;
   Float_t         beta2;
   Float_t         bphi1;
   Float_t         bphi2;
   Float_t         bm1;
   Float_t         bm2;
   Float_t         bsc1;
   Float_t         bsc2;
   Float_t         bmuf1;
   Float_t         bmuf2;
   Float_t         belf1;
   Float_t         belf2;
   //Float_t         bnhf1;
   //Float_t         bnhf2;
   Float_t         bqgl1;
   Float_t         bqgl2;
   /*
   Float_t         bbtag1;
   Float_t         bbtag2;
   Float_t         bbleptag1;
   Float_t         bbleptag2;
   Float_t         bctag1;
   Float_t         bctag2;
   Float_t         budstag1;
   Float_t         budstag2;
   Float_t         bgtag1;
   Float_t         bgtag2;
   */
   Int_t           bflav1;
   Int_t           bflav2;
   /*
   Int_t           bntrk1;
   Int_t           bntrk2;
   Int_t           bnputrk1;
   Int_t           bnputrk2;
   */
   Float_t         gen_bpt1;
   Float_t         gen_bpt2;
   Float_t         gen_beta1;
   Float_t         gen_beta2;
   Float_t         gen_bphi1;
   Float_t         gen_bphi2;
   Float_t         gen_bm1;
   Float_t         gen_bm2;

   // New for testing =>
   Int_t           gen_bLeadId1;
   Int_t           gen_bLeadId2;
   Int_t           gen_bFlags1;
   Int_t           gen_bFlags2;
   Float_t         gen_bXB1;
   Float_t         gen_bXB2;
   Float_t         gen_bNuCorr1;
   Float_t         gen_bNuCorr2;
   // <= New for testing

   Float_t         gpt1;
   Float_t         geta1;
   /*
   Float_t         gphi1;
   Float_t         gm1;
   Float_t         gsc1;
   Float_t         gmuf1;
   Float_t         gelf1;
   Float_t         gnhf1;
   */
   Float_t         gqgl1;
   /*
   Float_t         gbtag1;
   Float_t         gbleptag1;
   Float_t         gctag1;
   Float_t         gudstag1;
   Float_t         ggtag1;
   */
   Int_t           gflav1;
   /*
   Int_t           gntrk1;
   Int_t           gnputrk1;
   Float_t         gen_gpt1;
   Float_t         gen_geta1;
   Float_t         gen_gphi1;
   Float_t         gen_gm1;
   */

   Float_t         lep_pt;
   Float_t         lep_eta;
   Float_t         lep_phi;
   Float_t         lep_m;
   Float_t         hadTopPt;
   Float_t         hadTopPtG;
   Float_t         lepTopPtG;
   Float_t         fitProb;
   Float_t         weight;
   Float_t         tptweight;
   Int_t           run;
   Int_t           lumiBlock;
   //Int_t           event;
   ULong_t           event;
   Float_t         pfRho;
   Int_t           TruePU;
   Int_t           NPrVtx;
   Int_t           NPrVtxGood;
   Float_t         METx;
   Float_t         METy;
   Int_t           SIdx;
   Int_t           ComboType;
   vector<float>  *PSWgts;

   // List of branches
   TBranch        *b_pt1;   //!
   TBranch        *b_pt2;   //!
   TBranch        *b_eta1;   //!
   TBranch        *b_eta2;   //!
   TBranch        *b_phi1;   //!
   TBranch        *b_phi2;   //!
   TBranch        *b_m1;   //!
   TBranch        *b_m2;   //!
   TBranch        *b_sc1;   //!
   TBranch        *b_sc2;   //!
   TBranch        *b_muf1;   //!
   TBranch        *b_muf2;   //!
   TBranch        *b_elf1;   //!
   TBranch        *b_elf2;   //!
   TBranch        *b_nhf1;   //!
   TBranch        *b_nhf2;   //!
   TBranch        *b_qgl1;   //!
   TBranch        *b_qgl2;   //!
   TBranch        *b_btag1;   //!
   TBranch        *b_btag2;   //!
   TBranch        *b_bleptag1;   //!
   TBranch        *b_bleptag2;   //!
   TBranch        *b_ctag1;   //!
   TBranch        *b_ctag2;   //!
   TBranch        *b_udstag1;   //!
   TBranch        *b_udstag2;   //!
   TBranch        *b_gtag1;   //!
   TBranch        *b_gtag2;   //!
   TBranch        *b_flav1;   //!
   TBranch        *b_flav2;   //!
   TBranch        *b_ntrk1;   //!
   TBranch        *b_ntrk2;   //!
   TBranch        *b_nputrk1;   //!
   TBranch        *b_nputrk2;   //!
   TBranch        *b_recoWEta;   //!
   TBranch        *b_recoWPhi;   //!
   TBranch        *b_recoWMass;   //!
   TBranch        *b_gen_pt1;   //!
   TBranch        *b_gen_pt2;   //!
   TBranch        *b_gen_eta1;   //!
   TBranch        *b_gen_eta2;   //!
   TBranch        *b_gen_phi1;   //!
   TBranch        *b_gen_phi2;   //!
   TBranch        *b_gen_m1;   //!
   TBranch        *b_gen_m2;   //!
   TBranch        *b_bpt1;   //!
   TBranch        *b_bpt2;   //!
   TBranch        *b_beta1;   //!
   TBranch        *b_beta2;   //!
   TBranch        *b_bphi1;   //!
   TBranch        *b_bphi2;   //!
   TBranch        *b_bm1;   //!
   TBranch        *b_bm2;   //!
   TBranch        *b_bsc1;   //!
   TBranch        *b_bsc2;   //!
   TBranch        *b_bmuf1;   //!
   TBranch        *b_bmuf2;   //!
   TBranch        *b_belf1;   //!
   TBranch        *b_belf2;   //!
   TBranch        *b_bnhf1;   //!
   TBranch        *b_bnhf2;   //!
   TBranch        *b_bqgl1;   //!
   TBranch        *b_bqgl2;   //!
   TBranch        *b_bbtag1;   //!
   TBranch        *b_bbtag2;   //!
   TBranch        *b_bbleptag1;   //!
   TBranch        *b_bbleptag2;   //!
   TBranch        *b_bctag1;   //!
   TBranch        *b_bctag2;   //!
   TBranch        *b_budstag1;   //!
   TBranch        *b_budstag2;   //!
   TBranch        *b_bgtag1;   //!
   TBranch        *b_bgtag2;   //!
   TBranch        *b_bflav1;   //!
   TBranch        *b_bflav2;   //!
   TBranch        *b_bntrk1;   //!
   TBranch        *b_bntrk2;   //!
   TBranch        *b_bnputrk1;   //!
   TBranch        *b_bnputrk2;   //!
   TBranch        *b_gen_bpt1;   //!
   TBranch        *b_gen_bpt2;   //!
   TBranch        *b_gen_beta1;   //!
   TBranch        *b_gen_beta2;   //!
   TBranch        *b_gen_bphi1;   //!
   TBranch        *b_gen_bphi2;   //!
   TBranch        *b_gen_bm1;   //!
   TBranch        *b_gen_bm2;   //!
   // New for testing =>
   TBranch        *b_gen_leadBId1;   //!                                        
   TBranch        *b_gen_leadBId2;   //!                                        
   TBranch        *b_gen_bFlags1;   //!
   TBranch        *b_gen_bFlags2;   //!
   TBranch        *b_gen_bXB1;   //!                                            
   TBranch        *b_gen_bXB2;   //!                                            
   TBranch        *b_gen_bnucorr1;   //!                                        
   TBranch        *b_gen_bnucorr2;   //!     
   // <= New for testing
   TBranch        *b_gpt1;   //!
   TBranch        *b_geta1;   //!
   TBranch        *b_gphi1;   //!
   TBranch        *b_gm1;   //!
   TBranch        *b_gsc1;   //!
   TBranch        *b_gmuf1;   //!
   TBranch        *b_gelf1;   //!
   TBranch        *b_gnhf1;   //!
   TBranch        *b_gqgl1;   //!
   TBranch        *b_gbtag1;   //!
   TBranch        *b_gbleptag1;   //!
   TBranch        *b_gctag1;   //!
   TBranch        *b_gudstag1;   //!
   TBranch        *b_ggtag1;   //!
   TBranch        *b_gflav1;   //!
   TBranch        *b_gntrk1;   //!
   TBranch        *b_gnputrk1;   //!
   TBranch        *b_gen_gpt1;   //!
   TBranch        *b_gen_geta1;   //!
   TBranch        *b_gen_gphi1;   //!
   TBranch        *b_gen_gm1;   //!
   TBranch        *b_lep_pt;   //!
   TBranch        *b_lep_eta;   //!
   TBranch        *b_lep_phi;   //!
   TBranch        *b_lep_m;   //!
   TBranch        *b_hadTopPt;   //!
   TBranch        *b_hadTopPtGen;   //!
   TBranch        *b_lepTopPtGen;   //!
   TBranch        *b_fitProb;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_tptweight;   //!
   TBranch        *b_run;   //!
   TBranch        *b_lumiBlock;   //!
   TBranch        *b_event;   //!
   TBranch        *b_pfRho;   //!
   TBranch        *b_TruePU;   //!
   TBranch        *b_NPrVtx;   //!
   TBranch        *b_NPrVtxGood;   //!
   TBranch        *b_METx;   //!
   TBranch        *b_METy;   //!
   TBranch        *b_ComboType;   //!
   TBranch        *b_SampleIdx;   //!
   TBranch        *b_PSWgts;   //!

   std::string _s;
   Int_t runMC;
   Double_t runMCWeight;
   
   hadW(TTree *tree=0, std::string s="MC17");
   virtual ~hadW();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   void Draw(string mode);//string smc, string sdt);
   void DrawFP(string spt, string mode);
   void DrawRMS(string mode);
   void Draw2D(string set);
};

#endif

#ifdef hadW_cxx
hadW::hadW(TTree *tree, string s) : fChain(0), _s(s)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("tree",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("tree","");
      chain->Add("rootfiles/HadW/WmassMC17.root/tree");
      chain->Add("rootfiles/HadW/WmassUL17.root/tree");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

hadW::~hadW()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t hadW::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t hadW::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void hadW::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   PSWgts = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("pt1", &pt1, &b_pt1);
   fChain->SetBranchAddress("pt2", &pt2, &b_pt2);
   fChain->SetBranchAddress("eta1", &eta1, &b_eta1);
   fChain->SetBranchAddress("eta2", &eta2, &b_eta2);
   fChain->SetBranchAddress("phi1", &phi1, &b_phi1);
   fChain->SetBranchAddress("phi2", &phi2, &b_phi2);
   fChain->SetBranchAddress("m1", &m1, &b_m1);
   fChain->SetBranchAddress("m2", &m2, &b_m2);
   fChain->SetBranchAddress("sc1", &sc1, &b_sc1);
   fChain->SetBranchAddress("sc2", &sc2, &b_sc2);
   fChain->SetBranchAddress("muf1", &muf1, &b_muf1);
   fChain->SetBranchAddress("muf2", &muf2, &b_muf2);
   fChain->SetBranchAddress("elf1", &elf1, &b_elf1);
   fChain->SetBranchAddress("elf2", &elf2, &b_elf2);
   fChain->SetBranchAddress("nhf1", &nhf1, &b_nhf1);
   fChain->SetBranchAddress("nhf2", &nhf2, &b_nhf2);
   fChain->SetBranchAddress("qgl1", &qgl1, &b_qgl1);
   fChain->SetBranchAddress("qgl2", &qgl2, &b_qgl2);
   fChain->SetBranchAddress("btag1", &btag1, &b_btag1);
   fChain->SetBranchAddress("btag2", &btag2, &b_btag2);
   fChain->SetBranchAddress("bleptag1", &bleptag1, &b_bleptag1);
   fChain->SetBranchAddress("bleptag2", &bleptag2, &b_bleptag2);
   fChain->SetBranchAddress("ctag1", &ctag1, &b_ctag1);
   fChain->SetBranchAddress("ctag2", &ctag2, &b_ctag2);
   fChain->SetBranchAddress("udstag1", &udstag1, &b_udstag1);
   fChain->SetBranchAddress("udstag2", &udstag2, &b_udstag2);
   fChain->SetBranchAddress("gtag1", &gtag1, &b_gtag1);
   fChain->SetBranchAddress("gtag2", &gtag2, &b_gtag2);
   fChain->SetBranchAddress("flav1", &flav1, &b_flav1);
   fChain->SetBranchAddress("flav2", &flav2, &b_flav2);
   /*
   fChain->SetBranchAddress("ntrk1", &ntrk1, &b_ntrk1);
   fChain->SetBranchAddress("ntrk2", &ntrk2, &b_ntrk2);
   fChain->SetBranchAddress("nputrk1", &nputrk1, &b_nputrk1);
   fChain->SetBranchAddress("nputrk2", &nputrk2, &b_nputrk2);
   fChain->SetBranchAddress("recoWEta", &recoWEta, &b_recoWEta);
   fChain->SetBranchAddress("recoWPhi", &recoWPhi, &b_recoWPhi);
   */
   fChain->SetBranchAddress("recoWMass", &recoWMass, &b_recoWMass);
   fChain->SetBranchAddress("gen_pt1", &gen_pt1, &b_gen_pt1);
   fChain->SetBranchAddress("gen_pt2", &gen_pt2, &b_gen_pt2);
   fChain->SetBranchAddress("gen_eta1", &gen_eta1, &b_gen_eta1);
   fChain->SetBranchAddress("gen_eta2", &gen_eta2, &b_gen_eta2);
   fChain->SetBranchAddress("gen_phi1", &gen_phi1, &b_gen_phi1);
   fChain->SetBranchAddress("gen_phi2", &gen_phi2, &b_gen_phi2);
   fChain->SetBranchAddress("gen_m1", &gen_m1, &b_gen_m1);
   fChain->SetBranchAddress("gen_m2", &gen_m2, &b_gen_m2);
   fChain->SetBranchAddress("bpt1", &bpt1, &b_bpt1);
   fChain->SetBranchAddress("bpt2", &bpt2, &b_bpt2);
   fChain->SetBranchAddress("beta1", &beta1, &b_beta1);
   fChain->SetBranchAddress("beta2", &beta2, &b_beta2);
   fChain->SetBranchAddress("bphi1", &bphi1, &b_bphi1);
   fChain->SetBranchAddress("bphi2", &bphi2, &b_bphi2);
   fChain->SetBranchAddress("bm1", &bm1, &b_bm1);
   fChain->SetBranchAddress("bm2", &bm2, &b_bm2);
   fChain->SetBranchAddress("bsc1", &bsc1, &b_bsc1);
   fChain->SetBranchAddress("bsc2", &bsc2, &b_bsc2);
   fChain->SetBranchAddress("bmuf1", &bmuf1, &b_bmuf1);
   fChain->SetBranchAddress("bmuf2", &bmuf2, &b_bmuf2);
   fChain->SetBranchAddress("belf1", &belf1, &b_belf1);
   fChain->SetBranchAddress("belf2", &belf2, &b_belf2);
   //fChain->SetBranchAddress("bnhf1", &bnhf1, &b_bnhf1);
   //fChain->SetBranchAddress("bnhf2", &bnhf2, &b_bnhf2);
   fChain->SetBranchAddress("bqgl1", &bqgl1, &b_bqgl1);
   fChain->SetBranchAddress("bqgl2", &bqgl2, &b_bqgl2);
   /*
   fChain->SetBranchAddress("bbtag1", &bbtag1, &b_bbtag1);
   fChain->SetBranchAddress("bbtag2", &bbtag2, &b_bbtag2);
   fChain->SetBranchAddress("bbleptag1", &bbleptag1, &b_bbleptag1);
   fChain->SetBranchAddress("bbleptag2", &bbleptag2, &b_bbleptag2);
   fChain->SetBranchAddress("bctag1", &bctag1, &b_bctag1);
   fChain->SetBranchAddress("bctag2", &bctag2, &b_bctag2);
   fChain->SetBranchAddress("budstag1", &budstag1, &b_budstag1);
   fChain->SetBranchAddress("budstag2", &budstag2, &b_budstag2);
   fChain->SetBranchAddress("bgtag1", &bgtag1, &b_bgtag1);
   fChain->SetBranchAddress("bgtag2", &bgtag2, &b_bgtag2);
   */
   fChain->SetBranchAddress("bflav1", &bflav1, &b_bflav1);
   fChain->SetBranchAddress("bflav2", &bflav2, &b_bflav2);
   /*
   fChain->SetBranchAddress("bntrk1", &bntrk1, &b_bntrk1);
   fChain->SetBranchAddress("bntrk2", &bntrk2, &b_bntrk2);
   fChain->SetBranchAddress("bnputrk1", &bnputrk1, &b_bnputrk1);
   fChain->SetBranchAddress("bnputrk2", &bnputrk2, &b_bnputrk2);
   */
   fChain->SetBranchAddress("gen_bpt1", &gen_bpt1, &b_gen_bpt1);
   fChain->SetBranchAddress("gen_bpt2", &gen_bpt2, &b_gen_bpt2);
   fChain->SetBranchAddress("gen_beta1", &gen_beta1, &b_gen_beta1);
   fChain->SetBranchAddress("gen_beta2", &gen_beta2, &b_gen_beta2);
   fChain->SetBranchAddress("gen_bphi1", &gen_bphi1, &b_gen_bphi1);
   fChain->SetBranchAddress("gen_bphi2", &gen_bphi2, &b_gen_bphi2);
   fChain->SetBranchAddress("gen_bm1", &gen_bm1, &b_gen_bm1);
   fChain->SetBranchAddress("gen_bm2", &gen_bm2, &b_gen_bm2);

   // New for testing =>
   fChain->SetBranchAddress("gen_bLeadId1", &gen_bLeadId1, &b_gen_leadBId1);
   fChain->SetBranchAddress("gen_bLeadId2", &gen_bLeadId2, &b_gen_leadBId2);
   fChain->SetBranchAddress("gen_bFlags1", &gen_bFlags1, &b_gen_bFlags1);
   fChain->SetBranchAddress("gen_bFlags2", &gen_bFlags2, &b_gen_bFlags2);
   fChain->SetBranchAddress("gen_bXB1", &gen_bXB1, &b_gen_bXB1);
   fChain->SetBranchAddress("gen_bXB2", &gen_bXB2, &b_gen_bXB2);
   fChain->SetBranchAddress("gen_bNuCorr1", &gen_bNuCorr1, &b_gen_bnucorr1);
   fChain->SetBranchAddress("gen_bNuCorr2", &gen_bNuCorr2, &b_gen_bnucorr2);
   // <= New for testing

   fChain->SetBranchAddress("gpt1", &gpt1, &b_gpt1);
   fChain->SetBranchAddress("geta1", &geta1, &b_geta1);
   /*
   fChain->SetBranchAddress("gphi1", &gphi1, &b_gphi1);
   fChain->SetBranchAddress("gm1", &gm1, &b_gm1);
   fChain->SetBranchAddress("gsc1", &gsc1, &b_gsc1);
   fChain->SetBranchAddress("gmuf1", &gmuf1, &b_gmuf1);
   fChain->SetBranchAddress("gelf1", &gelf1, &b_gelf1);
   fChain->SetBranchAddress("gnhf1", &gnhf1, &b_gnhf1);
   */
   fChain->SetBranchAddress("gqgl1", &gqgl1, &b_gqgl1);
   /*
   fChain->SetBranchAddress("gbtag1", &gbtag1, &b_gbtag1);
   fChain->SetBranchAddress("gbleptag1", &gbleptag1, &b_gbleptag1);
   fChain->SetBranchAddress("gctag1", &gctag1, &b_gctag1);
   fChain->SetBranchAddress("gudstag1", &gudstag1, &b_gudstag1);
   fChain->SetBranchAddress("ggtag1", &ggtag1, &b_ggtag1);
   */
   fChain->SetBranchAddress("gflav1", &gflav1, &b_gflav1);
   /*
   fChain->SetBranchAddress("gntrk1", &gntrk1, &b_gntrk1);
   fChain->SetBranchAddress("gnputrk1", &gnputrk1, &b_gnputrk1);
   fChain->SetBranchAddress("gen_gpt1", &gen_gpt1, &b_gen_gpt1);
   fChain->SetBranchAddress("gen_geta1", &gen_geta1, &b_gen_geta1);
   fChain->SetBranchAddress("gen_gphi1", &gen_gphi1, &b_gen_gphi1);
   fChain->SetBranchAddress("gen_gm1", &gen_gm1, &b_gen_gm1);
   */

   fChain->SetBranchAddress("lep_pt", &lep_pt, &b_lep_pt);
   fChain->SetBranchAddress("lep_eta", &lep_eta, &b_lep_eta);
   fChain->SetBranchAddress("lep_phi", &lep_phi, &b_lep_phi);
   fChain->SetBranchAddress("lep_m", &lep_m, &b_lep_m);
   fChain->SetBranchAddress("hadTopPt", &hadTopPt, &b_hadTopPt);
   fChain->SetBranchAddress("hadTopPtG", &hadTopPtG, &b_hadTopPtGen);
   fChain->SetBranchAddress("lepTopPtG", &lepTopPtG, &b_lepTopPtGen);
   fChain->SetBranchAddress("fitProb", &fitProb, &b_fitProb);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("tptweight", &tptweight, &b_tptweight);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("pfRho", &pfRho, &b_pfRho);
   fChain->SetBranchAddress("TruePU", &TruePU, &b_TruePU);
   fChain->SetBranchAddress("NPrVtx", &NPrVtx, &b_NPrVtx);
   fChain->SetBranchAddress("NPrVtxGood", &NPrVtxGood, &b_NPrVtxGood);
   fChain->SetBranchAddress("METx", &METx, &b_METx);
   fChain->SetBranchAddress("METy", &METy, &b_METy);
   fChain->SetBranchAddress("SIdx", &SIdx, &b_SampleIdx);
   fChain->SetBranchAddress("ComboType", &ComboType, &b_ComboType);
   fChain->SetBranchAddress("PSWgts", &PSWgts, &b_PSWgts);

   Notify();
}

Bool_t hadW::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.
  
  cout << "Notify: new file #" << fCurrent << endl << flush;
  if (fChain->GetCurrentFile()) {
    map<string, int> runMap;
    runMap["rootfiles/HadW/UL16_L3Res_V7/Muo16APV_MC.root"] = 274000;
    runMap["rootfiles/HadW/UL16_L3Res_V7/Muo16_MC.root"] =    282000;
    runMap["rootfiles/HadW/UL17V5_WithGlu/Muo17_MC.root"] =   298000;
    runMap["rootfiles/HadW/UL17V5_WithGlu/Ele17_MC.root"] =   300000;
    runMap["rootfiles/HadW/UL18V5_WithGlu/Muo18_MC.root"] =   316000;
    runMap["rootfiles/HadW/UL18V5_WithGlu/Ele18_MC.root"] =   318000;
    //
    double w16a = 19.7;
    double w16g = 16.8;
    double w17  = 41.5;
    double w18  = 59.9;
    double sumw = w16a + w16g + w17 + w18;
    map<string, double> weightMap;
    weightMap["rootfiles/HadW/UL16_L3Res_V7/Muo16APV_MC.root"] = w16a / sumw;
    weightMap["rootfiles/HadW/UL16_L3Res_V7/Muo16_MC.root"] = w16g / sumw;
    weightMap["rootfiles/HadW/UL17V5_WithGlu/Muo17_MC.root"] = 0.5 * w17 / sumw;
    weightMap["rootfiles/HadW/UL17V5_WithGlu/Ele17_MC.root"] = 0.5 * w17 / sumw;
    weightMap["rootfiles/HadW/UL18V5_WithGlu/Muo18_MC.root"] = 0.5 * w18 / sumw;
    weightMap["rootfiles/HadW/UL18V5_WithGlu/Ele18_MC.root"] = 0.5 * w18 / sumw;

    runMC = runMap[fChain->GetCurrentFile()->GetName()];
    runMCWeight = weightMap[fChain->GetCurrentFile()->GetName()];
    if (runMCWeight==0) runMCWeight = 1;
    cout << "(" << fChain->GetCurrentFile()->GetName()
	 << " for run " << runMC
	 << " weight " << runMCWeight << ")" << endl << flush;
  }

   return kTRUE;
}

void hadW::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t hadW::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef hadW_cxx
