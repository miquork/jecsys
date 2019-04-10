{

  // For JEC (for uncertainty)
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+");

  // For JEC uncertainty
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+");

  // For JER (used in ptresolution.C)
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetResolutionObject.cc+g");
  gROOT->ProcessLine(".L JetMETCorrections/Modules/src/JetResolution.cc+g");

  //gROOT->ProcessLine(".L minitools/ptresolution.h+g");
  gROOT->ProcessLine(".L minitools/drawDeltaJEC.C+g");

  /*
  unfold("0.0-0.5",1);
  unfold("0.5-1.0",2);
  unfold("1.0-1.5",3);
  unfold("1.5-2.0",4);
  unfold("2.0-2.5",5);
  unfold("2.5-3.0",6);
  unfold("3.2-4.7",7);
  */

  drawDeltaJEC("0.0-0.5");
  drawDeltaJEC("0.5-1.0");
  drawDeltaJEC("1.0-1.5");
  drawDeltaJEC("1.5-2.0");
  drawDeltaJEC("2.0-2.5");
  drawDeltaJEC("2.5-3.0");
  drawDeltaJEC("3.2-4.7");

}
