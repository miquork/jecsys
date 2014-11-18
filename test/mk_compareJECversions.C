{
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+");

  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+");

  gROOT->ProcessLine(".L compareJECversions.C+");
  gROOT->ProcessLine(".exception");

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
  */

  compareJECversions("AK5PF",true,false,false); // L1
  compareJECversions("AK5PFchs",true,false,false); // L1
  compareJECversions("AK5PFchs",false,true,false); // L2L3
  compareJECversions("AK5PFchs",false,false,true); // Res
}
