{

  // Official JME resolutions
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+g");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetResolutionObject.cc+g");
  gROOT->ProcessLine(".L JetMETCorrections/Modules/src/JetResolution.cc+g");

  gROOT->ProcessLine(".L minitools/tools.C+g");
  gROOT->ProcessLine(".L minitools/resolution.C+g");
  
  resolution();
}
