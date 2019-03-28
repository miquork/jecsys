{

  // Official JME resolutions
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+g");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetResolutionObject.cc+g");
  gROOT->ProcessLine(".L JetMETCorrections/Modules/src/JetResolution.cc+g");

  gROOT->ProcessLine(".L minitools/tools.C+g");
  gROOT->ProcessLine(".L minitools/resolution.C+g");
  
  //resolution();
  //resolution("MC","Fall18V8-D");
  //resolution("P8CP5","17nov17-DE");
  //resolution("MC","Legacy16-GH");
  //resolution("HW","Legacy16-GH");

  redoJER("Run2018");
  redoJER("Run2017");
  redoJER("Run2016");
  redoJER("Run1");
}
