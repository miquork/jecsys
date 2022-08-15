//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Aug 10 13:24:52 2022 by ROOT version 6.24/04
// from TTree Events/Events
// found on file: rootfiles/tst_1K_runIII_AlCaRaw.root
//////////////////////////////////////////////////////////

#ifndef alcajme_h
#define alcajme_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"

class alcajme {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          luminosityBlock;
   ULong64_t       event;
   Float_t         HepMCGenEvent_scale;
   Float_t         GenEventInfo_qScale;
   Int_t           pileupInfo_BX0_numPUInteractions;
   Float_t         pileupInfo_BX0_numTrueInteractions;
   Float_t         pileupInfo_BX0_max_pT_hats;
   UInt_t          pileupInfo_BX0_n_pThat000to020;
   UInt_t          pileupInfo_BX0_n_pThat020to030;
   UInt_t          pileupInfo_BX0_n_pThat030to050;
   UInt_t          pileupInfo_BX0_n_pThat050to080;
   UInt_t          pileupInfo_BX0_n_pThat080to120;
   UInt_t          pileupInfo_BX0_n_pThat120to170;
   UInt_t          pileupInfo_BX0_n_pThat170to300;
   UInt_t          pileupInfo_BX0_n_pThat300to470;
   UInt_t          pileupInfo_BX0_n_pThat470to600;
   UInt_t          pileupInfo_BX0_n_pThat600toInf;
   vector<float>   *hltPFJetsCorrectedMatchedToCaloJets10_pt;
   vector<float>   *hltPFJetsCorrectedMatchedToCaloJets10_eta;
   vector<float>   *hltPFJetsCorrectedMatchedToCaloJets10_phi;
   vector<float>   *hltPFJetsCorrectedMatchedToCaloJets10_mass;
   vector<float>   *hltPFJetsCorrectedMatchedToCaloJets10_jesc;
   vector<float>   *hltPFJetsCorrectedMatchedToCaloJets10_jetArea;
   vector<unsigned int> *hltPFJetsCorrectedMatchedToCaloJets10_numberOfDaughters;
   vector<float>   *hltPFJetsCorrectedMatchedToCaloJets10_chargedHadronEnergyFraction;
   vector<float>   *hltPFJetsCorrectedMatchedToCaloJets10_neutralHadronEnergyFraction;
   vector<float>   *hltPFJetsCorrectedMatchedToCaloJets10_electronEnergyFraction;
   vector<float>   *hltPFJetsCorrectedMatchedToCaloJets10_photonEnergyFraction;
   vector<float>   *hltPFJetsCorrectedMatchedToCaloJets10_muonEnergyFraction;
   vector<int>     *hltPFJetsCorrectedMatchedToCaloJets10_chargedHadronMultiplicity;
   vector<int>     *hltPFJetsCorrectedMatchedToCaloJets10_neutralHadronMultiplicity;
   vector<int>     *hltPFJetsCorrectedMatchedToCaloJets10_electronMultiplicity;
   vector<int>     *hltPFJetsCorrectedMatchedToCaloJets10_photonMultiplicity;
   vector<int>     *hltPFJetsCorrectedMatchedToCaloJets10_muonMultiplicity;
   vector<float>   *hltPFJetsCorrectedMatchedToCaloJets10AK8_pt;
   vector<float>   *hltPFJetsCorrectedMatchedToCaloJets10AK8_eta;
   vector<float>   *hltPFJetsCorrectedMatchedToCaloJets10AK8_phi;
   vector<float>   *hltPFJetsCorrectedMatchedToCaloJets10AK8_mass;
   vector<float>   *hltPFJetsCorrectedMatchedToCaloJets10AK8_jesc;
   vector<float>   *hltPFJetsCorrectedMatchedToCaloJets10AK8_jetArea;
   vector<unsigned int> *hltPFJetsCorrectedMatchedToCaloJets10AK8_numberOfDaughters;
   vector<float>   *hltPFJetsCorrectedMatchedToCaloJets10AK8_chargedHadronEnergyFraction;
   vector<float>   *hltPFJetsCorrectedMatchedToCaloJets10AK8_neutralHadronEnergyFraction;
   vector<float>   *hltPFJetsCorrectedMatchedToCaloJets10AK8_electronEnergyFraction;
   vector<float>   *hltPFJetsCorrectedMatchedToCaloJets10AK8_photonEnergyFraction;
   vector<float>   *hltPFJetsCorrectedMatchedToCaloJets10AK8_muonEnergyFraction;
   vector<int>     *hltPFJetsCorrectedMatchedToCaloJets10AK8_chargedHadronMultiplicity;
   vector<int>     *hltPFJetsCorrectedMatchedToCaloJets10AK8_neutralHadronMultiplicity;
   vector<int>     *hltPFJetsCorrectedMatchedToCaloJets10AK8_electronMultiplicity;
   vector<int>     *hltPFJetsCorrectedMatchedToCaloJets10AK8_photonMultiplicity;
   vector<int>     *hltPFJetsCorrectedMatchedToCaloJets10AK8_muonMultiplicity;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_event;   //!
   TBranch        *b_HepMCGenEvent_scale;   //!
   TBranch        *b_GenEventInfo_qScale;   //!
   TBranch        *b_pileupInfo_BX0_numPUInteractions;   //!
   TBranch        *b_pileupInfo_BX0_numTrueInteractions;   //!
   TBranch        *b_pileupInfo_BX0_max_pT_hats;   //!
   TBranch        *b_pileupInfo_BX0_n_pThat000to020;   //!
   TBranch        *b_pileupInfo_BX0_n_pThat020to030;   //!
   TBranch        *b_pileupInfo_BX0_n_pThat030to050;   //!
   TBranch        *b_pileupInfo_BX0_n_pThat050to080;   //!
   TBranch        *b_pileupInfo_BX0_n_pThat080to120;   //!
   TBranch        *b_pileupInfo_BX0_n_pThat120to170;   //!
   TBranch        *b_pileupInfo_BX0_n_pThat170to300;   //!
   TBranch        *b_pileupInfo_BX0_n_pThat300to470;   //!
   TBranch        *b_pileupInfo_BX0_n_pThat470to600;   //!
   TBranch        *b_pileupInfo_BX0_n_pThat600toInf;   //!
   TBranch        *b_hltPFJetsCorrectedMatchedToCaloJets10_pt;   //!
   TBranch        *b_hltPFJetsCorrectedMatchedToCaloJets10_eta;   //!
   TBranch        *b_hltPFJetsCorrectedMatchedToCaloJets10_phi;   //!
   TBranch        *b_hltPFJetsCorrectedMatchedToCaloJets10_mass;   //!
   TBranch        *b_hltPFJetsCorrectedMatchedToCaloJets10_jesc;   //!
   TBranch        *b_hltPFJetsCorrectedMatchedToCaloJets10_jetArea;   //!
   TBranch        *b_hltPFJetsCorrectedMatchedToCaloJets10_numberOfDaughters;   //!
   TBranch        *b_hltPFJetsCorrectedMatchedToCaloJets10_chargedHadronEnergyFraction;   //!
   TBranch        *b_hltPFJetsCorrectedMatchedToCaloJets10_neutralHadronEnergyFraction;   //!
   TBranch        *b_hltPFJetsCorrectedMatchedToCaloJets10_electronEnergyFraction;   //!
   TBranch        *b_hltPFJetsCorrectedMatchedToCaloJets10_photonEnergyFraction;   //!
   TBranch        *b_hltPFJetsCorrectedMatchedToCaloJets10_muonEnergyFraction;   //!
   TBranch        *b_hltPFJetsCorrectedMatchedToCaloJets10_chargedHadronMultiplicity;   //!
   TBranch        *b_hltPFJetsCorrectedMatchedToCaloJets10_neutralHadronMultiplicity;   //!
   TBranch        *b_hltPFJetsCorrectedMatchedToCaloJets10_electronMultiplicity;   //!
   TBranch        *b_hltPFJetsCorrectedMatchedToCaloJets10_photonMultiplicity;   //!
   TBranch        *b_hltPFJetsCorrectedMatchedToCaloJets10_muonMultiplicity;   //!
   TBranch        *b_hltPFJetsCorrectedMatchedToCaloJets10AK8_pt;   //!
   TBranch        *b_hltPFJetsCorrectedMatchedToCaloJets10AK8_eta;   //!
   TBranch        *b_hltPFJetsCorrectedMatchedToCaloJets10AK8_phi;   //!
   TBranch        *b_hltPFJetsCorrectedMatchedToCaloJets10AK8_mass;   //!
   TBranch        *b_hltPFJetsCorrectedMatchedToCaloJets10AK8_jesc;   //!
   TBranch        *b_hltPFJetsCorrectedMatchedToCaloJets10AK8_jetArea;   //!
   TBranch        *b_hltPFJetsCorrectedMatchedToCaloJets10AK8_numberOfDaughters;   //!
   TBranch        *b_hltPFJetsCorrectedMatchedToCaloJets10AK8_chargedHadronEnergyFraction;   //!
   TBranch        *b_hltPFJetsCorrectedMatchedToCaloJets10AK8_neutralHadronEnergyFraction;   //!
   TBranch        *b_hltPFJetsCorrectedMatchedToCaloJets10AK8_electronEnergyFraction;   //!
   TBranch        *b_hltPFJetsCorrectedMatchedToCaloJets10AK8_photonEnergyFraction;   //!
   TBranch        *b_hltPFJetsCorrectedMatchedToCaloJets10AK8_muonEnergyFraction;   //!
   TBranch        *b_hltPFJetsCorrectedMatchedToCaloJets10AK8_chargedHadronMultiplicity;   //!
   TBranch        *b_hltPFJetsCorrectedMatchedToCaloJets10AK8_neutralHadronMultiplicity;   //!
   TBranch        *b_hltPFJetsCorrectedMatchedToCaloJets10AK8_electronMultiplicity;   //!
   TBranch        *b_hltPFJetsCorrectedMatchedToCaloJets10AK8_photonMultiplicity;   //!
   TBranch        *b_hltPFJetsCorrectedMatchedToCaloJets10AK8_muonMultiplicity;   //!

   alcajme(TTree *tree=0);
   virtual ~alcajme();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef alcajme_cxx
alcajme::alcajme(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  const bool only1M = false;
   if (tree == 0) {
     //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("rootfiles/tst_1K_runIII_AlCaRaw.root");
     //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../data/AlCaRaw/hadded_ntuple.root");
     //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../data/AlCaRaw/skimmed_ntuple.root");
     TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(only1M ? 
"../data/AlCaRaw/skimmed_ntuple.root" :
"../data/AlCaRaw/hadded_ntuple.root");
      if (!f || !f->IsOpen()) {
	//f = new TFile("rootfiles/tst_1K_runIII_AlCaRaw.root");
	//f = new TFile("../data/AlCaRaw/hadded_ntuple.root");
	//f = new TFile("../data/AlCaRaw/skimmed_ntuple.root");
	f = new TFile(only1M ? "../data/AlCaRaw/skimmed_ntuple.root":
		      "../data/AlCaRaw/hadded_ntuple.root");
      }
      //TDirectory * dir = (TDirectory*)f->Get("rootfiles/tst_1K_runIII_AlCaRaw.root:/JMETriggerNTuple");
      //TDirectory * dir = (TDirectory*)f->Get("../data/AlCaRaw/hadded_ntuple.root:/JMETriggerNTuple");
      if (only1M) f->GetObject("Events",tree); // skimmed_ntuple.root
      else {
	TDirectory * dir = (TDirectory*)f->Get("../data/AlCaRaw/hadded_ntuple.root:/JMETriggerNTuple");
	dir->GetObject("Events",tree);
      }

   }
   Init(tree);
}

alcajme::~alcajme()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t alcajme::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t alcajme::LoadTree(Long64_t entry)
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

void alcajme::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   hltPFJetsCorrectedMatchedToCaloJets10_pt = 0;
   hltPFJetsCorrectedMatchedToCaloJets10_eta = 0;
   hltPFJetsCorrectedMatchedToCaloJets10_phi = 0;
   hltPFJetsCorrectedMatchedToCaloJets10_mass = 0;
   hltPFJetsCorrectedMatchedToCaloJets10_jesc = 0;
   hltPFJetsCorrectedMatchedToCaloJets10_jetArea = 0;
   hltPFJetsCorrectedMatchedToCaloJets10_numberOfDaughters = 0;
   hltPFJetsCorrectedMatchedToCaloJets10_chargedHadronEnergyFraction = 0;
   hltPFJetsCorrectedMatchedToCaloJets10_neutralHadronEnergyFraction = 0;
   hltPFJetsCorrectedMatchedToCaloJets10_electronEnergyFraction = 0;
   hltPFJetsCorrectedMatchedToCaloJets10_photonEnergyFraction = 0;
   hltPFJetsCorrectedMatchedToCaloJets10_muonEnergyFraction = 0;
   hltPFJetsCorrectedMatchedToCaloJets10_chargedHadronMultiplicity = 0;
   hltPFJetsCorrectedMatchedToCaloJets10_neutralHadronMultiplicity = 0;
   hltPFJetsCorrectedMatchedToCaloJets10_electronMultiplicity = 0;
   hltPFJetsCorrectedMatchedToCaloJets10_photonMultiplicity = 0;
   hltPFJetsCorrectedMatchedToCaloJets10_muonMultiplicity = 0;
   hltPFJetsCorrectedMatchedToCaloJets10AK8_pt = 0;
   hltPFJetsCorrectedMatchedToCaloJets10AK8_eta = 0;
   hltPFJetsCorrectedMatchedToCaloJets10AK8_phi = 0;
   hltPFJetsCorrectedMatchedToCaloJets10AK8_mass = 0;
   hltPFJetsCorrectedMatchedToCaloJets10AK8_jesc = 0;
   hltPFJetsCorrectedMatchedToCaloJets10AK8_jetArea = 0;
   hltPFJetsCorrectedMatchedToCaloJets10AK8_numberOfDaughters = 0;
   hltPFJetsCorrectedMatchedToCaloJets10AK8_chargedHadronEnergyFraction = 0;
   hltPFJetsCorrectedMatchedToCaloJets10AK8_neutralHadronEnergyFraction = 0;
   hltPFJetsCorrectedMatchedToCaloJets10AK8_electronEnergyFraction = 0;
   hltPFJetsCorrectedMatchedToCaloJets10AK8_photonEnergyFraction = 0;
   hltPFJetsCorrectedMatchedToCaloJets10AK8_muonEnergyFraction = 0;
   hltPFJetsCorrectedMatchedToCaloJets10AK8_chargedHadronMultiplicity = 0;
   hltPFJetsCorrectedMatchedToCaloJets10AK8_neutralHadronMultiplicity = 0;
   hltPFJetsCorrectedMatchedToCaloJets10AK8_electronMultiplicity = 0;
   hltPFJetsCorrectedMatchedToCaloJets10AK8_photonMultiplicity = 0;
   hltPFJetsCorrectedMatchedToCaloJets10AK8_muonMultiplicity = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("HepMCGenEvent_scale", &HepMCGenEvent_scale, &b_HepMCGenEvent_scale);
   fChain->SetBranchAddress("GenEventInfo_qScale", &GenEventInfo_qScale, &b_GenEventInfo_qScale);
   fChain->SetBranchAddress("pileupInfo_BX0_numPUInteractions", &pileupInfo_BX0_numPUInteractions, &b_pileupInfo_BX0_numPUInteractions);
   fChain->SetBranchAddress("pileupInfo_BX0_numTrueInteractions", &pileupInfo_BX0_numTrueInteractions, &b_pileupInfo_BX0_numTrueInteractions);
   fChain->SetBranchAddress("pileupInfo_BX0_max_pT_hats", &pileupInfo_BX0_max_pT_hats, &b_pileupInfo_BX0_max_pT_hats);
   fChain->SetBranchAddress("pileupInfo_BX0_n_pThat000to020", &pileupInfo_BX0_n_pThat000to020, &b_pileupInfo_BX0_n_pThat000to020);
   fChain->SetBranchAddress("pileupInfo_BX0_n_pThat020to030", &pileupInfo_BX0_n_pThat020to030, &b_pileupInfo_BX0_n_pThat020to030);
   fChain->SetBranchAddress("pileupInfo_BX0_n_pThat030to050", &pileupInfo_BX0_n_pThat030to050, &b_pileupInfo_BX0_n_pThat030to050);
   fChain->SetBranchAddress("pileupInfo_BX0_n_pThat050to080", &pileupInfo_BX0_n_pThat050to080, &b_pileupInfo_BX0_n_pThat050to080);
   fChain->SetBranchAddress("pileupInfo_BX0_n_pThat080to120", &pileupInfo_BX0_n_pThat080to120, &b_pileupInfo_BX0_n_pThat080to120);
   fChain->SetBranchAddress("pileupInfo_BX0_n_pThat120to170", &pileupInfo_BX0_n_pThat120to170, &b_pileupInfo_BX0_n_pThat120to170);
   fChain->SetBranchAddress("pileupInfo_BX0_n_pThat170to300", &pileupInfo_BX0_n_pThat170to300, &b_pileupInfo_BX0_n_pThat170to300);
   fChain->SetBranchAddress("pileupInfo_BX0_n_pThat300to470", &pileupInfo_BX0_n_pThat300to470, &b_pileupInfo_BX0_n_pThat300to470);
   fChain->SetBranchAddress("pileupInfo_BX0_n_pThat470to600", &pileupInfo_BX0_n_pThat470to600, &b_pileupInfo_BX0_n_pThat470to600);
   fChain->SetBranchAddress("pileupInfo_BX0_n_pThat600toInf", &pileupInfo_BX0_n_pThat600toInf, &b_pileupInfo_BX0_n_pThat600toInf);
   fChain->SetBranchAddress("hltPFJetsCorrectedMatchedToCaloJets10_pt", &hltPFJetsCorrectedMatchedToCaloJets10_pt, &b_hltPFJetsCorrectedMatchedToCaloJets10_pt);
   fChain->SetBranchAddress("hltPFJetsCorrectedMatchedToCaloJets10_eta", &hltPFJetsCorrectedMatchedToCaloJets10_eta, &b_hltPFJetsCorrectedMatchedToCaloJets10_eta);
   fChain->SetBranchAddress("hltPFJetsCorrectedMatchedToCaloJets10_phi", &hltPFJetsCorrectedMatchedToCaloJets10_phi, &b_hltPFJetsCorrectedMatchedToCaloJets10_phi);
   fChain->SetBranchAddress("hltPFJetsCorrectedMatchedToCaloJets10_mass", &hltPFJetsCorrectedMatchedToCaloJets10_mass, &b_hltPFJetsCorrectedMatchedToCaloJets10_mass);
   fChain->SetBranchAddress("hltPFJetsCorrectedMatchedToCaloJets10_jesc", &hltPFJetsCorrectedMatchedToCaloJets10_jesc, &b_hltPFJetsCorrectedMatchedToCaloJets10_jesc);
   fChain->SetBranchAddress("hltPFJetsCorrectedMatchedToCaloJets10_jetArea", &hltPFJetsCorrectedMatchedToCaloJets10_jetArea, &b_hltPFJetsCorrectedMatchedToCaloJets10_jetArea);
   fChain->SetBranchAddress("hltPFJetsCorrectedMatchedToCaloJets10_numberOfDaughters", &hltPFJetsCorrectedMatchedToCaloJets10_numberOfDaughters, &b_hltPFJetsCorrectedMatchedToCaloJets10_numberOfDaughters);
   fChain->SetBranchAddress("hltPFJetsCorrectedMatchedToCaloJets10_chargedHadronEnergyFraction", &hltPFJetsCorrectedMatchedToCaloJets10_chargedHadronEnergyFraction, &b_hltPFJetsCorrectedMatchedToCaloJets10_chargedHadronEnergyFraction);
   fChain->SetBranchAddress("hltPFJetsCorrectedMatchedToCaloJets10_neutralHadronEnergyFraction", &hltPFJetsCorrectedMatchedToCaloJets10_neutralHadronEnergyFraction, &b_hltPFJetsCorrectedMatchedToCaloJets10_neutralHadronEnergyFraction);
   fChain->SetBranchAddress("hltPFJetsCorrectedMatchedToCaloJets10_electronEnergyFraction", &hltPFJetsCorrectedMatchedToCaloJets10_electronEnergyFraction, &b_hltPFJetsCorrectedMatchedToCaloJets10_electronEnergyFraction);
   fChain->SetBranchAddress("hltPFJetsCorrectedMatchedToCaloJets10_photonEnergyFraction", &hltPFJetsCorrectedMatchedToCaloJets10_photonEnergyFraction, &b_hltPFJetsCorrectedMatchedToCaloJets10_photonEnergyFraction);
   fChain->SetBranchAddress("hltPFJetsCorrectedMatchedToCaloJets10_muonEnergyFraction", &hltPFJetsCorrectedMatchedToCaloJets10_muonEnergyFraction, &b_hltPFJetsCorrectedMatchedToCaloJets10_muonEnergyFraction);
   fChain->SetBranchAddress("hltPFJetsCorrectedMatchedToCaloJets10_chargedHadronMultiplicity", &hltPFJetsCorrectedMatchedToCaloJets10_chargedHadronMultiplicity, &b_hltPFJetsCorrectedMatchedToCaloJets10_chargedHadronMultiplicity);
   fChain->SetBranchAddress("hltPFJetsCorrectedMatchedToCaloJets10_neutralHadronMultiplicity", &hltPFJetsCorrectedMatchedToCaloJets10_neutralHadronMultiplicity, &b_hltPFJetsCorrectedMatchedToCaloJets10_neutralHadronMultiplicity);
   fChain->SetBranchAddress("hltPFJetsCorrectedMatchedToCaloJets10_electronMultiplicity", &hltPFJetsCorrectedMatchedToCaloJets10_electronMultiplicity, &b_hltPFJetsCorrectedMatchedToCaloJets10_electronMultiplicity);
   fChain->SetBranchAddress("hltPFJetsCorrectedMatchedToCaloJets10_photonMultiplicity", &hltPFJetsCorrectedMatchedToCaloJets10_photonMultiplicity, &b_hltPFJetsCorrectedMatchedToCaloJets10_photonMultiplicity);
   fChain->SetBranchAddress("hltPFJetsCorrectedMatchedToCaloJets10_muonMultiplicity", &hltPFJetsCorrectedMatchedToCaloJets10_muonMultiplicity, &b_hltPFJetsCorrectedMatchedToCaloJets10_muonMultiplicity);
   fChain->SetBranchAddress("hltPFJetsCorrectedMatchedToCaloJets10AK8_pt", &hltPFJetsCorrectedMatchedToCaloJets10AK8_pt, &b_hltPFJetsCorrectedMatchedToCaloJets10AK8_pt);
   fChain->SetBranchAddress("hltPFJetsCorrectedMatchedToCaloJets10AK8_eta", &hltPFJetsCorrectedMatchedToCaloJets10AK8_eta, &b_hltPFJetsCorrectedMatchedToCaloJets10AK8_eta);
   fChain->SetBranchAddress("hltPFJetsCorrectedMatchedToCaloJets10AK8_phi", &hltPFJetsCorrectedMatchedToCaloJets10AK8_phi, &b_hltPFJetsCorrectedMatchedToCaloJets10AK8_phi);
   fChain->SetBranchAddress("hltPFJetsCorrectedMatchedToCaloJets10AK8_mass", &hltPFJetsCorrectedMatchedToCaloJets10AK8_mass, &b_hltPFJetsCorrectedMatchedToCaloJets10AK8_mass);
   fChain->SetBranchAddress("hltPFJetsCorrectedMatchedToCaloJets10AK8_jesc", &hltPFJetsCorrectedMatchedToCaloJets10AK8_jesc, &b_hltPFJetsCorrectedMatchedToCaloJets10AK8_jesc);
   fChain->SetBranchAddress("hltPFJetsCorrectedMatchedToCaloJets10AK8_jetArea", &hltPFJetsCorrectedMatchedToCaloJets10AK8_jetArea, &b_hltPFJetsCorrectedMatchedToCaloJets10AK8_jetArea);
   fChain->SetBranchAddress("hltPFJetsCorrectedMatchedToCaloJets10AK8_numberOfDaughters", &hltPFJetsCorrectedMatchedToCaloJets10AK8_numberOfDaughters, &b_hltPFJetsCorrectedMatchedToCaloJets10AK8_numberOfDaughters);
   fChain->SetBranchAddress("hltPFJetsCorrectedMatchedToCaloJets10AK8_chargedHadronEnergyFraction", &hltPFJetsCorrectedMatchedToCaloJets10AK8_chargedHadronEnergyFraction, &b_hltPFJetsCorrectedMatchedToCaloJets10AK8_chargedHadronEnergyFraction);
   fChain->SetBranchAddress("hltPFJetsCorrectedMatchedToCaloJets10AK8_neutralHadronEnergyFraction", &hltPFJetsCorrectedMatchedToCaloJets10AK8_neutralHadronEnergyFraction, &b_hltPFJetsCorrectedMatchedToCaloJets10AK8_neutralHadronEnergyFraction);
   fChain->SetBranchAddress("hltPFJetsCorrectedMatchedToCaloJets10AK8_electronEnergyFraction", &hltPFJetsCorrectedMatchedToCaloJets10AK8_electronEnergyFraction, &b_hltPFJetsCorrectedMatchedToCaloJets10AK8_electronEnergyFraction);
   fChain->SetBranchAddress("hltPFJetsCorrectedMatchedToCaloJets10AK8_photonEnergyFraction", &hltPFJetsCorrectedMatchedToCaloJets10AK8_photonEnergyFraction, &b_hltPFJetsCorrectedMatchedToCaloJets10AK8_photonEnergyFraction);
   fChain->SetBranchAddress("hltPFJetsCorrectedMatchedToCaloJets10AK8_muonEnergyFraction", &hltPFJetsCorrectedMatchedToCaloJets10AK8_muonEnergyFraction, &b_hltPFJetsCorrectedMatchedToCaloJets10AK8_muonEnergyFraction);
   fChain->SetBranchAddress("hltPFJetsCorrectedMatchedToCaloJets10AK8_chargedHadronMultiplicity", &hltPFJetsCorrectedMatchedToCaloJets10AK8_chargedHadronMultiplicity, &b_hltPFJetsCorrectedMatchedToCaloJets10AK8_chargedHadronMultiplicity);
   fChain->SetBranchAddress("hltPFJetsCorrectedMatchedToCaloJets10AK8_neutralHadronMultiplicity", &hltPFJetsCorrectedMatchedToCaloJets10AK8_neutralHadronMultiplicity, &b_hltPFJetsCorrectedMatchedToCaloJets10AK8_neutralHadronMultiplicity);
   fChain->SetBranchAddress("hltPFJetsCorrectedMatchedToCaloJets10AK8_electronMultiplicity", &hltPFJetsCorrectedMatchedToCaloJets10AK8_electronMultiplicity, &b_hltPFJetsCorrectedMatchedToCaloJets10AK8_electronMultiplicity);
   fChain->SetBranchAddress("hltPFJetsCorrectedMatchedToCaloJets10AK8_photonMultiplicity", &hltPFJetsCorrectedMatchedToCaloJets10AK8_photonMultiplicity, &b_hltPFJetsCorrectedMatchedToCaloJets10AK8_photonMultiplicity);
   fChain->SetBranchAddress("hltPFJetsCorrectedMatchedToCaloJets10AK8_muonMultiplicity", &hltPFJetsCorrectedMatchedToCaloJets10AK8_muonMultiplicity, &b_hltPFJetsCorrectedMatchedToCaloJets10AK8_muonMultiplicity);
   Notify();
}

Bool_t alcajme::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void alcajme::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t alcajme::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef alcajme_cxx
