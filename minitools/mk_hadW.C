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
  //cmc->AddFile("rootfiles/HadW/UL17/WmassMC17.root");
  cmc->AddFile("rootfiles/HadW/UL18/WmassMC18.root");
  TChain *cdt = new TChain("tree","tree");
  //cdt->AddFile("rootfiles/HadW/UL17/WmassUL17.root");
  //cdt->AddFile("rootfiles/HadW/UL18/WmassUL18.root");
  cdt->AddFile("rootfiles/HadW/UL18/AEl.root");
  cdt->AddFile("rootfiles/HadW/UL18/AMu.root");
  
  //hadW wmc(cmc,"MC17");
  hadW wmc(cmc,"MC18");
  //wmc.Loop();

  //hadW wdt(cdt,"UL17");
  hadW wdt(cdt,"UL18");
  wdt.Loop();

  wmc.Draw();//"MC17","UL17");
  wmc.DrawFP("ptboth");
  wmc.DrawFP("ptave");

  /*
  wmc.Draw2D("MC17");
  wmc.Draw2D("UL17");
  */
}
