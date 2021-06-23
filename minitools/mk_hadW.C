//#include "minitools/hadW.h"
// NB: if fails to compile, clean out .so .d .pcm also from JetMETObjects/src,
// then run mk_reprocess.C first

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

void mk_hadW(string mode = "18V5") {
//void mk_hadW(string mode = "1718V5") {
//void mk_hadW(string mode = "161718V5") {
//void mk_hadW(string mode = "16V2") { // L2Res
//void mk_hadW(string mode = "16GHV5") { // L2L3Res
//void mk_hadW(string mode = "16APVV5") {
//void mk_hadW(string mode = "16BCDV5") {
//void mk_hadW(string mode = "16EFV5") {
//void mk_hadW(string mode = "17V5") {
//void mk_hadW(string mode = "16APVV7") {
//void mk_hadW(string mode = "16BCDV7") {
//void mk_hadW(string mode = "16EFV7") {
//void mk_hadW(string mode = "16GHV7") {
//void mk_hadW(string mode = "17V5New") {
//void mk_hadW(string mode = "161718V7") {

  if (!(mode=="16APVV5" || mode=="16BCDV5" || mode=="16EFV5" || mode=="16V2" ||
	mode=="16GHV5" || mode=="17V5" || mode=="18V5" ||
	mode=="16APVV7" || mode=="16BCDV7"||mode=="16EFV7" || mode=="16GHV7" ||
	mode=="1718V5" || mode=="161718V5" || mode=="161718V7" ||
	mode=="17V5New")) return;
  
  // Caveat for 1718V5: veto regions may not work properly (only one used),
  // same with JEC reapplication

  TChain *cmc = new TChain("tree","tree");
  if (mode=="16V2") {
    cmc->AddFile("rootfiles/HadW/UL16GH/Muo16_MC.root"); // L2Res
    //cmc->AddFile("rootfiles/HadW/UL16GH/Ele16_MC.root");
  }
  if (mode=="16GHV5" || mode=="161718V5") {
    cmc->AddFile("rootfiles/HadW/UL16_L3Res_V5/Muo16_MC.root"); // L2L3Res
  }
  if (mode=="16GHV7" || mode=="161718V7") {
    cmc->AddFile("rootfiles/HadW/UL16_L3Res_V7/Muo16_MC.root"); // L2L3Res
  }
  if (mode=="16APVV5" || mode=="161718V5") {
    cmc->AddFile("rootfiles/HadW/UL16_L3Res_V5/Muo16APV_MC.root"); // L2L3Res
    //cmc->AddFile("rootfiles/HadW/UL16BCDEF/Muo16APV_MC.root"); // L2Res
    //cmc->AddFile("rootfiles/HadW/UL16GH/Ele16_MC.root");
  }
  if (mode=="16APVV7" || mode=="161718V7") {
    cmc->AddFile("rootfiles/HadW/UL16_L3Res_V7/Muo16APV_MC.root"); // L2L3Res
  }
  if (mode=="16BCDV5") {
    cmc->AddFile("rootfiles/HadW/UL16_L3Res_V5/Muo16BCD_MC.root"); // L2L3Res
    //cmc->AddFile("rootfiles/HadW/UL16BCDEF/Muo16BCD_MC.root"); // L2Res
  }
  if (mode=="16BCDV7") {
    cmc->AddFile("rootfiles/HadW/UL16_L3Res_V7/Muo16BCD_MC.root"); // L2L3Res
  }
  if (mode=="16EFV5") {
    cmc->AddFile("rootfiles/HadW/UL16_L3Res_V5/Muo16EF_MC.root"); // L2L3Res
    //cmc->AddFile("rootfiles/HadW/UL16BCDEF/Muo16EF_MC.root"); // L2Res
  }
  if (mode=="16EFV7") {
    cmc->AddFile("rootfiles/HadW/UL16_L3Res_V7/Muo16EF_MC.root"); // L2L3Res
  }
  //if (mode=="17V5") cmc->AddFile("rootfiles/HadW/UL17/WmassMC17.root");
  if (mode=="17V5" || mode=="1718V5" || mode=="161718V5" || mode=="161718V7") {
    //cmc->AddFile("rootfiles/HadW/UL17V5/WMass_Muo17_MadGraph.root");
    //cmc->AddFile("rootfiles/HadW/UL17V5/WMass_Ele17_MadGraph.root");
    //cmc->AddFile("rootfiles/HadW/UL17V5_NoMu/WMass_Muo17_PowHeg.root");
    //cmc->AddFile("rootfiles/HadW/UL17V5_NoMu/WMass_Ele17_PowHeg.root");
    //cmc->AddFile("rootfiles/HadW/UL17V5_WithMu/Muo17_MC.root");
    //cmc->AddFile("rootfiles/HadW/UL17V5_WithMu/Ele17_MC.root");

    //cmc->AddFile("rootfiles/HadW/UL17V5_WithEMu/Muo17_MC.root");
    //cmc->AddFile("rootfiles/HadW/UL17V5_WithEMu/Ele17_MC.root");
    cmc->AddFile("rootfiles/HadW/UL17V5_WithGlu/Muo17_MC.root");
    cmc->AddFile("rootfiles/HadW/UL17V5_WithGlu/Ele17_MC.root");
  }
  if (mode=="17V5New") {
    cmc->AddFile("rootfiles/HadW/UL17V5_WithGlu/NewTestingMuo17_MC.root");
  }
  //if (mode=="18V5") cmc->AddFile("rootfiles/HadW/UL18/WmassMC18.root");
  if (mode=="18V5" || mode=="1718V5" || mode=="161718V5" || mode=="161718V7") {
    //cmc->AddFile("rootfiles/HadW/UL18V5/WMass_Muo18_MadGraph.root");
    //cmc->AddFile("rootfiles/HadW/UL18V5/WMass_Ele18_MadGraph.root");
    //cmc->AddFile("rootfiles/HadW/UL18V5_NoMu/WMass_Muo18_PowHeg.root");
    //cmc->AddFile("rootfiles/HadW/UL18V5_NoMu/WMass_Ele18_PowHeg.root");
    //cmc->AddFile("rootfiles/HadW/UL18V5_WithMu/Muo18_MC.root");
    //cmc->AddFile("rootfiles/HadW/UL18V5_WithMu/Ele18_MC.root");

    //cmc->AddFile("rootfiles/HadW/UL18V5_WithEMu/Muo18_MC.root");
    //cmc->AddFile("rootfiles/HadW/UL18V5_WithEMu/Ele18_MC.root");
    cmc->AddFile("rootfiles/HadW/UL18V5_WithGlu/Muo18_MC.root");
    cmc->AddFile("rootfiles/HadW/UL18V5_WithGlu/Ele18_MC.root");
  }
  TChain *cdt = new TChain("tree","tree");
  if (mode=="16V2") {
    //cdt->AddFile("rootfiles/HadW/UL16GH/Muo16_Run2016.root");
    cdt->AddFile("rootfiles/HadW/UL16GH/Muo16_DATA.root"); // L2Res
  }
  if (mode=="16GHV5" || mode=="161718V5") {
    cdt->AddFile("rootfiles/HadW/UL16_L3Res_V5/Muo16_DATA.root"); // L2L3Res
  }
  if (mode=="16GHV7" || mode=="161718V7") {
    cdt->AddFile("rootfiles/HadW/UL16_L3Res_V7/Muo16_DATA.root"); // L2L3Res
  }
  if (mode=="16APVV5" || mode=="161718V5") {
    //cdt->AddFile("rootfiles/HadW/UL16BCDEF/Muo16APV_DATA.root"); // L2Res
    cdt->AddFile("rootfiles/HadW/UL16_L3Res_V5/Muo16APV_DATA.root"); // L2L3Res
  }
  if (mode=="16APVV7" || mode=="161718V7") {
    cdt->AddFile("rootfiles/HadW/UL16_L3Res_V7/Muo16APV_DATA.root"); // L2L3Res
  }
  if (mode=="16BCDV5") {
    //cdt->AddFile("rootfiles/HadW/UL16BCDEF/Muo16BCD_DATA.root"); // L2Res
    cdt->AddFile("rootfiles/HadW/UL16_L3Res_V5/Muo16BCD_DATA.root"); // L2L3Res
  }
  if (mode=="16BCDV7") {
    cdt->AddFile("rootfiles/HadW/UL16_L3Res_V7/Muo16BCD_DATA.root"); // L2L3Res
  }
  if (mode=="16EFV5") {
    //cdt->AddFile("rootfiles/HadW/UL16BCDEF/Muo16EF_DATA.root"); // L2Res
    cdt->AddFile("rootfiles/HadW/UL16_L3Res_V5/Muo16EF_DATA.root"); // L2L3Res
  }
  if (mode=="16EFV7") {
    cdt->AddFile("rootfiles/HadW/UL16_L3Res_V7/Muo16EF_DATA.root"); // L2L3Res
  }
  //if (mode=="17V5") cdt->AddFile("rootfiles/HadW/UL17/WmassUL17.root");
  if (mode=="17V5" || mode=="1718V5" || mode=="161718V5" || mode=="17V5New" ||
      mode=="161718V7") {
    /*
    cdt->AddFile("rootfiles/HadW/UL17V5_NoMu/WMass_Muo17_DTB.root");
    cdt->AddFile("rootfiles/HadW/UL17V5_NoMu/WMass_Muo17_DTC.root");
    cdt->AddFile("rootfiles/HadW/UL17V5_NoMu/WMass_Muo17_DTD.root");
    cdt->AddFile("rootfiles/HadW/UL17V5_NoMu/WMass_Muo17_DTE.root");
    cdt->AddFile("rootfiles/HadW/UL17V5_NoMu/WMass_Muo17_DTF.root");
    cdt->AddFile("rootfiles/HadW/UL17V5_NoMu/WMass_Ele17_DTB.root");
    cdt->AddFile("rootfiles/HadW/UL17V5_NoMu/WMass_Ele17_DTC.root");
    cdt->AddFile("rootfiles/HadW/UL17V5_NoMu/WMass_Ele17_DTD.root");
    cdt->AddFile("rootfiles/HadW/UL17V5_NoMu/WMass_Ele17_DTE.root");
    cdt->AddFile("rootfiles/HadW/UL17V5_NoMu/WMass_Ele17_DTF.root");
    */
    /*
    cdt->AddFile("rootfiles/HadW/UL17V5_WithMu/Muo17_Run2017B.root");
    cdt->AddFile("rootfiles/HadW/UL17V5_WithMu/Muo17_Run2017C.root");
    cdt->AddFile("rootfiles/HadW/UL17V5_WithMu/Muo17_Run2017D.root");
    cdt->AddFile("rootfiles/HadW/UL17V5_WithMu/Muo17_Run2017E.root");
    cdt->AddFile("rootfiles/HadW/UL17V5_WithMu/Muo17_Run2017F.root");
    */
    /*
    cdt->AddFile("rootfiles/HadW/UL17V5_WithEMu/Muo17_Run2017B.root");
    cdt->AddFile("rootfiles/HadW/UL17V5_WithEMu/Muo17_Run2017C.root");
    cdt->AddFile("rootfiles/HadW/UL17V5_WithEMu/Muo17_Run2017D.root");
    cdt->AddFile("rootfiles/HadW/UL17V5_WithEMu/Muo17_Run2017E.root");
    cdt->AddFile("rootfiles/HadW/UL17V5_WithEMu/Muo17_Run2017F.root");
    cdt->AddFile("rootfiles/HadW/UL17V5_WithEMu/Ele17_Run2017B.root");
    cdt->AddFile("rootfiles/HadW/UL17V5_WithEMu/Ele17_Run2017C.root");
    cdt->AddFile("rootfiles/HadW/UL17V5_WithEMu/Ele17_Run2017D.root");
    cdt->AddFile("rootfiles/HadW/UL17V5_WithEMu/Ele17_Run2017E.root");
    cdt->AddFile("rootfiles/HadW/UL17V5_WithEMu/Ele17_Run2017F.root");
    */
    cdt->AddFile("rootfiles/HadW/UL17V5_WithGlu/Muo17_Run2017B.root");
    cdt->AddFile("rootfiles/HadW/UL17V5_WithGlu/Muo17_Run2017C.root");
    cdt->AddFile("rootfiles/HadW/UL17V5_WithGlu/Muo17_Run2017D.root");
    cdt->AddFile("rootfiles/HadW/UL17V5_WithGlu/Muo17_Run2017E.root");
    cdt->AddFile("rootfiles/HadW/UL17V5_WithGlu/Muo17_Run2017F.root");
    cdt->AddFile("rootfiles/HadW/UL17V5_WithGlu/Ele17_Run2017B.root");
    cdt->AddFile("rootfiles/HadW/UL17V5_WithGlu/Ele17_Run2017C.root");
    cdt->AddFile("rootfiles/HadW/UL17V5_WithGlu/Ele17_Run2017D.root");
    cdt->AddFile("rootfiles/HadW/UL17V5_WithGlu/Ele17_Run2017E.root");
    cdt->AddFile("rootfiles/HadW/UL17V5_WithGlu/Ele17_Run2017F.root");
  }
  //if (mode=="18V5") cdt->AddFile("rootfiles/HadW/UL18/WmassUL18.root");
  if (mode=="18V5" || mode=="1718V5" || mode=="161718V5" || mode=="161718V7") {
    /*
    cdt->AddFile("rootfiles/HadW/UL18V5_NoMu/WMass_Muo18_DTA.root");
    cdt->AddFile("rootfiles/HadW/UL18V5_NoMu/WMass_Muo18_DTB.root");
    cdt->AddFile("rootfiles/HadW/UL18V5_NoMu/WMass_Muo18_DTC.root");
    cdt->AddFile("rootfiles/HadW/UL18V5_NoMu/WMass_Muo18_DTD.root");
    cdt->AddFile("rootfiles/HadW/UL18V5_NoMu/WMass_Ele18_DTA.root");
    cdt->AddFile("rootfiles/HadW/UL18V5_NoMu/WMass_Ele18_DTB.root");
    cdt->AddFile("rootfiles/HadW/UL18V5_NoMu/WMass_Ele18_DTC.root");
    cdt->AddFile("rootfiles/HadW/UL18V5_NoMu/WMass_Ele18_DTD.root");
    */
    /*
    cdt->AddFile("rootfiles/HadW/UL18V5_WithMu/Muo18_Run2018A.root");
    cdt->AddFile("rootfiles/HadW/UL18V5_WithMu/Muo18_Run2018B.root");
    cdt->AddFile("rootfiles/HadW/UL18V5_WithMu/Muo18_Run2018C.root");
    cdt->AddFile("rootfiles/HadW/UL18V5_WithMu/Muo18_Run2018D.root");
    */
    /*
    cdt->AddFile("rootfiles/HadW/UL18V5_WithEMu/Muo18_Run2018A.root");
    cdt->AddFile("rootfiles/HadW/UL18V5_WithEMu/Muo18_Run2018B.root");
    cdt->AddFile("rootfiles/HadW/UL18V5_WithEMu/Muo18_Run2018C.root");
    cdt->AddFile("rootfiles/HadW/UL18V5_WithEMu/Muo18_Run2018D.root");
    cdt->AddFile("rootfiles/HadW/UL18V5_WithEMu/Ele18_Run2018A.root");
    cdt->AddFile("rootfiles/HadW/UL18V5_WithEMu/Ele18_Run2018B.root");
    cdt->AddFile("rootfiles/HadW/UL18V5_WithEMu/Ele18_Run2018C.root");
    cdt->AddFile("rootfiles/HadW/UL18V5_WithEMu/Ele18_Run2018D.root");
    */
    cdt->AddFile("rootfiles/HadW/UL18V5_WithGlu/Muo18_Run2018A.root");
    cdt->AddFile("rootfiles/HadW/UL18V5_WithGlu/Muo18_Run2018B.root");
    cdt->AddFile("rootfiles/HadW/UL18V5_WithGlu/Muo18_Run2018C.root");
    cdt->AddFile("rootfiles/HadW/UL18V5_WithGlu/Muo18_Run2018D.root");
    cdt->AddFile("rootfiles/HadW/UL18V5_WithGlu/Ele18_Run2018A.root");
    cdt->AddFile("rootfiles/HadW/UL18V5_WithGlu/Ele18_Run2018B.root");
    cdt->AddFile("rootfiles/HadW/UL18V5_WithGlu/Ele18_Run2018C.root");
    cdt->AddFile("rootfiles/HadW/UL18V5_WithGlu/Ele18_Run2018D.root");
  }

  hadW wmc(cmc,"MC"+mode);
  wmc.Loop();

  hadW wdt(cdt,"UL"+mode);
  wdt.Loop();

  wmc.Draw(mode);
  //wmc.DrawFP("ptboth",mode);
  wmc.DrawFP("ptave",mode);
  wmc.DrawRMS(mode);

  /*
  wmc.Draw2D("MC17");
  wmc.Draw2D("UL17");
  */
  
  // Plotting scripts using these outputs:
  // minitools/drawATLASmt.C - plot ATLAS-style mw, rbq, mt and mlb
  // minitools/drawWjetMass.C - plot jet mass/pT and mass reco/gen ratios
  // minitools/drawWb.C - plot mbqq variants vs b-jet pT
  // add tools to plot time stability

  // drawRbq
}
