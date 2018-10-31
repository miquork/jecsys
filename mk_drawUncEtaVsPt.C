{
  // make sure asserts are run
  #undef NDEBUG
  // does not seem to be enough, also need to combile with +g

  string currentWorkingDir = gSystem->pwd();
  cout <<currentWorkingDir.c_str() <<endl;
  gSystem->AddIncludePath(Form("-I%s",currentWorkingDir.c_str()));



  // For JEC residual (and pile-up)
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
  gROOT->ProcessLine(".L CondFormats//JetMETObjects/src/JetCorrectorParameters.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+");
  //
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+");


  gROOT->ProcessLine(".L drawUncEtaVsPt.C+g");

  drawUncEtaVsPt("txt/Summer16_03Feb2017_V9_DATA_UncertaintySources_AK4PFchs.txt");
}
