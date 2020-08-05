//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jun 26 15:44:54 2020 by ROOT version 6.18/04
// from TChain tree/
//////////////////////////////////////////////////////////

#ifndef hadW_h
#define hadW_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <string>

using namespace std;

// Header file for the classes stored in the TTree if any.

class hadW {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         pt1;
   Float_t         pt2;
   Float_t         eta1;
   Float_t         eta2;
   Float_t         phi1;
   Float_t         phi2;
   Float_t         recoWMass;
   Float_t         fitProb;
   Float_t         pfRho;

   // List of branches
   TBranch        *b_pt1;   //!
   TBranch        *b_pt2;   //!
   TBranch        *b_eta1;   //!
   TBranch        *b_eta2;   //!
   TBranch        *b_phi1;   //!
   TBranch        *b_phi2;   //!
   TBranch        *b_recoWMass;   //!
   TBranch        *b_fitProb;   //!
   TBranch        *b_pfRho;   //!

   std::string _s;

   hadW(TTree *tree=0, std::string s="MC17");
   virtual ~hadW();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   void Draw();//string smc, string sdt);
   void DrawFP(string spt);
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
   fChain->SetBranchAddress("recoWMass", &recoWMass, &b_recoWMass);
   fChain->SetBranchAddress("fitProb", &fitProb, &b_fitProb);
   fChain->SetBranchAddress("pfRho", &pfRho, &b_pfRho);
   Notify();
}

Bool_t hadW::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

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
