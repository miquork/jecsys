// Adapted from D0 Experiment jetcorr/macros/mk_RjetUncertainty.C

{
  //gSystem->SetIncludePath("-Iinclude -I.");
  //cout << "Include path: " << gSystem->GetIncludePath() << endl;
  gROOT->ProcessLine(".L tdrstyle_mod.C");

  // For JEC central value
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+");
  // For JEC uncertainty
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+");

  gROOT->ProcessLine(".L ErrorTypes.cpp+");
  gROOT->ProcessLine(".L JECUncertainty.cpp+");

  gROOT->ProcessLine(".exception");
  gROOT->ProcessLine(".L drawJetCorrectionUncertainty.C+");
  
  setTDRStyle();

  // Single source test
  //drawJetCorrectionUncertainty("AK5PFchs"); // no source files
  // Print out source files (only setup for AK5PF + true pair)
  //drawJetCorrectionUncertainty("AK5PF",true); // also source files

  // Whole shebang

  // Switch JECUncert to JECSource with 'bool _useAbsUncert = false'
  drawJetCorrectionUncertainty("AK5PF",false); // no source files (quick)
  //drawJetCorrectionUncertainty("AK5PF,true"); // also source files
  drawJetCorrectionUncertainty("AK5PFchs");
  drawJetCorrectionUncertainty("AK7PF");
  drawJetCorrectionUncertainty("AK7PFchs");
  //drawJetCorrectionUncertainty("AK5CALO");
  //drawJetCorrectionUncertainty("AK7CALO");


  // Fixes for pas-v6
  //drawJetCorrectionUncertainty("AK5PFchs",false);
}
