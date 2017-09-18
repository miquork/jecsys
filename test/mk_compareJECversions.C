{
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+");

  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+");

  gROOT->ProcessLine(".L test/compareJECversions.C+");
  gROOT->ProcessLine(".exception");
  
  // Test Run2015D V3M2 L2L3Res
  //  compareJECversions("AK4PFchs",false,false,true,"DATA"); // Fig.24topleft (L2L3Res)
  //compareJECversions("AK4PFchs",true,true,true,"DATA"); // L1L2L3+Res
  //  compareJECversions("AK4PFchs",true,true,false,"DATA"); // L1L2L3
  //  compareJECversions("AK4PFchs",false,true,false,"DATA"); // L2L3
  compareJECversions("AK4PFchs",true,false,false,"DATA"); // L1
  compareJECversions("AK4PFchs",true,false,false,"MC"); // L1
  /*


  compareJECversions("AK4PFchs",false,false,true,"DATA"); // Fig.24topleft
  //compareJECversions("AK4PFchs",false,true,false,"MC"); // L2L3
  //compareJECversions("AK4PFchs",true,false,false,"MC"); // L1 MC

  /*
  compareJECversions("AK4PFchs",true,false,false,"MC"); // L1 MC
  compareJECversions("AK4PFchs",true,false,false,"DATA"); // L1
  compareJECversions("AK4PFchs",false,true,false,"DATA"); // L2L3
  compareJECversions("AK4PFchs",false,false,true,"DATA"); // Res

  compareJECversions("AK4PFchs",true,true,false,"DATA"); // L1L2L3
  compareJECversions("AK4PFchs",true,true,true,"DATA"); // L1L2L3+Res
  compareJECversions("AK4PFchs",false,true,true,"DATA"); // L2L3+Res
  */
  // Patch for pas-v6
  //compareJECversions("AK5PF",true,false,false); // Fig.9left(top)
  //compareJECversions("AK5PFchs",true,false,false); // Fig.9left(bottom)
  //compareJECversions("AK5PFchs",false,true,false); // Fig.13topleft
  //compareJECversions("AK5PFchs",false,false,true); // Fig.24topleft

  /*
  compareJECversions("AK5PFchs",true,true,true); // L1+L2L3+Res
  compareJECversions("AK5PF",   true,true,true); // L1+L2L3+Res
  compareJECversions("AK7PFchs",true,true,true); // L1+L2L3+Res
  compareJECversions("AK7PF",   true,true,true); // L1+L2L3+Res

  compareJECversions("AK5PFchs",false,true,true); // L2L3+Res
  compareJECversions("AK5PF",   false,true,true); // L2L3+Res
  compareJECversions("AK7PFchs",false,true,true); // L2L3+Res
  compareJECversions("AK7PF",   false,true,true); // L2L3+Res

  compareJECversions("AK5PFchs",false,false,true); // Res
  compareJECversions("AK5PF",   false,false,true); // Res
  compareJECversions("AK7PFchs",false,false,true); // Res
  compareJECversions("AK7PF",   false,false,true); // Res

  compareJECversions("AK5PFchs",false,true,false); // L2L3
  compareJECversions("AK5PF",   false,true,false); // L2L3
  compareJECversions("AK7PFchs",false,true,false); // L2L3
  compareJECversions("AK7PF",   false,true,false); // L2L3

  compareJECversions("AK5PFchs",true,false,false); // L1
  compareJECversions("AK5PF",   true,false,false); // L1
  compareJECversions("AK7PFchs",true,false,false); // L1
  compareJECversions("AK7PF",   true,false,false); // L1


  compareJECversions("AK5PF",true,false,false); // L1
  compareJECversions("AK5PFchs",true,false,false); // L1
  compareJECversions("AK5PFchs",false,true,false); // L2L3
  compareJECversions("AK5PFchs",false,false,true); // Res
  compareJECversions("AK5PFchs",true,true,true); // All
  */
}
