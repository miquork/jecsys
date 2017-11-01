// Adapted from D0 Experiment jetcorr/macros/mk_RjetUncertainty.C

{
  //gSystem->SetIncludePath("-Iinclude -I.");
  //cout << "Include path: " << gSystem->GetIncludePath() << endl;
  //gROOT->ProcessLine(".L tdrstyle_mod.C");

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

  //gROOT->ProcessLine(".exception");
  gROOT->ProcessLine(".L drawJetCorrectionUncertainty.C+");
  
  setTDRStyle();

  // NB: remember to delete old plots before pdflatex pdf/sysplots.tex
  // NB2: JECUncert vs JECSource with 'bool _useAbsUncert'

  // Single source test
  //drawJetCorrectionUncertainty("AK4PFchs",false); //no source files, minimal
  // Print out source files (only setup for AK4PFchs + true pair)
  drawJetCorrectionUncertainty("AK4PFchs",true); // also source files
  //drawJetCorrectionUncertainty("AK8PFchs",true); // also source files

  // Whole shebang

  // Switch JECUncert to JECSource with 'bool _useAbsUncert = false'
  /*
  drawJetCorrectionUncertainty("AK4PFchs",true); // also source files
  drawJetCorrectionUncertainty("AK4PF",true); // also source files
  drawJetCorrectionUncertainty("AK8PFchs",true); // also source files
  drawJetCorrectionUncertainty("AK8PF",true); // also source files
  */
  //drawJetCorrectionUncertainty("AK4PFchs",false); // no source files
  //drawJetCorrectionUncertainty("AK4PF");
  //drawJetCorrectionUncertainty("AK8PFchs",true);
  //drawJetCorrectionUncertainty("AK8PF");
  //drawJetCorrectionUncertainty("AK4CALO");
  //drawJetCorrectionUncertainty("AK7CALO");


  // Final paper plots, only (run with both _absUncert=true, false)
  //drawJetCorrectionUncertainty("AK4PFchs",false);
}
