{
  // make sure asserts are run                                                  
  #undef NDEBUG
  // does not seem to be enough, also need to combile with +g                   

  // For JEC residual (and pile-up)                                             
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
  gROOT->ProcessLine(".L CondFormats//JetMETObjects/src/JetCorrectorParameters.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+");
  //                                                                            
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+");


  // Compile with +g to make sure asserts are run                               
  gROOT->ProcessLine(".L minitools/mergeL2L3ResTextFiles.C+g");

  //mergeL2L3ResTextFiles("2016APV");
  //mergeL2L3ResTextFiles("2016GH");
  //mergeL2L3ResTextFiles("2017",0.01);
  //mergeL2L3ResTextFiles("2018",0.002);
  mergeL2L3ResTextFiles("Run2",0.002);
}
