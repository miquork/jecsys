{

  // For JEC (for uncertainty)
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+");

  // For JEC uncertainty
  //gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+");
  //gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+");

  gROOT->ProcessLine(".L minitools/estimateMinJetPt.C+g");
  
  estimateMinJetPt();
}
