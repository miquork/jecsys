//#include "minitools/hadW.h"

// For JEC
#include "../CondFormats/JetMETObjects/src/Utilities.cc"
#include "../CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "../CondFormats/JetMETObjects/interface/SimpleJetCorrector.h"
#include "../CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
// For JEC uncertainty
//#include "CondFormats/JetMETObjects/interface/SimpleJetCorrectionUncertainty.h"
//#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+)
R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+)
//R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+)
//R__LOAD_LIBRARY(CondFormats/JetMETObjects/src/JetResolutionObject.cc+)
//R__LOAD_LIBRARY(JetMETCorrections/Modules/src/JetResolution.cc+)

R__LOAD_LIBRARY(minitools/hadW.C+g)

void mk_hadW() {
  
  TChain *cmc = new TChain("tree","tree");
  cmc->AddFile("rootfiles/HadW/WmassMC17.root");
  TChain *cdt = new TChain("tree","tree");
  cdt->AddFile("rootfiles/HadW/WmassUL17.root");

  hadW wmc(cmc,"MC17");
  wmc.Loop();

  hadW wdt(cdt,"UL17");
  wdt.Loop();

  wmc.Draw();//"MC17","UL17");
  wmc.DrawFP("ptboth");
  wmc.DrawFP("ptave");

}
