{
  // JEC central value
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+");
  // JEC uncertainties
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+");

  // drawing macro
  gROOT->ProcessLine(".L drawSourceCorrelations.C+");

  // run
  drawSourceCorrelations("Total");
  drawSourceCorrelations("TotalNoFlavorNoTime");
  drawSourceCorrelations("TotalNoTime","AK7PF");
}
