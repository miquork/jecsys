//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jun  2 08:56:45 2021 by ROOT version 6.18/04
// from TChain ak4/ProcessedTree/
//////////////////////////////////////////////////////////

#ifndef PFhadronLoop_h
#define PFhadronLoop_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class PFhadronLoop {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         prescale;
   Float_t         PtTrk;
   Float_t         Pt;
   Float_t         Eta;
   Float_t         Phi;
   Float_t         E;
   Float_t         Ef_ECALRaw;
   Float_t         Ef_ECAL;
   Float_t         Ef_HCALRaw;
   Float_t         Ef_HCAL;

   // List of branches
   TBranch        *b_prescale;   //!
   TBranch        *b_PtTrk;   //!
   TBranch        *b_Pt;   //!
   TBranch        *b_Eta;   //!
   TBranch        *b_Phi;   //!
   TBranch        *b_E;   //!
   TBranch        *b_Ef_ECALRaw;   //!
   TBranch        *b_Ef_ECAL;   //!
   TBranch        *b_Ef_HCALRaw;   //!
   TBranch        *b_Ef_HCAL;   //!

   PFhadronLoop(TTree *tree=0);
   virtual ~PFhadronLoop();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef PFhadronLoop_cxx
PFhadronLoop::PFhadronLoop(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("rootfiles/Hadrons_UL16/data/ZeroBias_UL16RunH.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("rootfiles/Hadrons_UL16/data/ZeroBias_UL16RunH.root");
      }
      f->GetObject("ak4/ProcessedTree",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("ak4/ProcessedTree","");
      chain->Add("rootfiles/Hadrons/Had18MC.root/ak4/ProcessedTree");
      chain->Add("rootfiles/Hadrons/Had18MC_HEM.root/ak4/ProcessedTree");
      chain->Add("rootfiles/Hadrons/Had17MC.root/ak4/ProcessedTree");
      chain->Add("rootfiles/Hadrons_UL16/mc/SingleNeutrino_UL16APV.root/ak4/ProcessedTree");
      chain->Add("rootfiles/Hadrons_UL16/mc/SingleNeutrino_UL16nonAPV.root/ak4/ProcessedTree");
      chain->Add("rootfiles/Hadrons/Had18A.root/ak4/ProcessedTree");
      chain->Add("rootfiles/Hadrons/Had18B.root/ak4/ProcessedTree");
      chain->Add("rootfiles/Hadrons/Had18C.root/ak4/ProcessedTree");
      chain->Add("rootfiles/Hadrons/Had18D.root/ak4/ProcessedTree");
      chain->Add("rootfiles/Hadrons/Had17B.root/ak4/ProcessedTree");
      chain->Add("rootfiles/Hadrons/Had17C.root/ak4/ProcessedTree");
      chain->Add("rootfiles/Hadrons/Had17D.root/ak4/ProcessedTree");
      chain->Add("rootfiles/Hadrons/Had17E.root/ak4/ProcessedTree");
      chain->Add("rootfiles/Hadrons/Had17F.root/ak4/ProcessedTree");
      chain->Add("rootfiles/Hadrons_UL16/data/ZeroBias_UL16RunB.root/ak4/ProcessedTree");
      chain->Add("rootfiles/Hadrons_UL16/data/ZeroBias_UL16RunC.root/ak4/ProcessedTree");
      chain->Add("rootfiles/Hadrons_UL16/data/ZeroBias_UL16RunD.root/ak4/ProcessedTree");
      chain->Add("rootfiles/Hadrons_UL16/data/ZeroBias_UL16RunE.root/ak4/ProcessedTree");
      chain->Add("rootfiles/Hadrons_UL16/data/ZeroBias_UL16RunFe.root/ak4/ProcessedTree");
      chain->Add("rootfiles/Hadrons_UL16/data/ZeroBias_UL16RunFl.root/ak4/ProcessedTree");
      chain->Add("rootfiles/Hadrons_UL16/data/ZeroBias_UL16RunG.root/ak4/ProcessedTree");
      chain->Add("rootfiles/Hadrons_UL16/data/ZeroBias_UL16RunH.root/ak4/ProcessedTree");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

PFhadronLoop::~PFhadronLoop()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t PFhadronLoop::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t PFhadronLoop::LoadTree(Long64_t entry)
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

void PFhadronLoop::Init(TTree *tree)
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

   fChain->SetBranchAddress("prescale", &prescale, &b_prescale);
   fChain->SetBranchAddress("PtTrk", &PtTrk, &b_PtTrk);
   fChain->SetBranchAddress("Pt", &Pt, &b_Pt);
   fChain->SetBranchAddress("Eta", &Eta, &b_Eta);
   fChain->SetBranchAddress("Phi", &Phi, &b_Phi);
   fChain->SetBranchAddress("E", &E, &b_E);
   fChain->SetBranchAddress("Ef_ECALRaw", &Ef_ECALRaw, &b_Ef_ECALRaw);
   fChain->SetBranchAddress("Ef_ECAL", &Ef_ECAL, &b_Ef_ECAL);
   fChain->SetBranchAddress("Ef_HCALRaw", &Ef_HCALRaw, &b_Ef_HCALRaw);
   fChain->SetBranchAddress("Ef_HCAL", &Ef_HCAL, &b_Ef_HCAL);
   Notify();
}

Bool_t PFhadronLoop::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void PFhadronLoop::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t PFhadronLoop::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef PFhadronLoop_cxx
