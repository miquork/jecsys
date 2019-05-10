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

  //redoJER("Run2018");
  //redoJER("Run2018ABC");
  //redoJER("Run2018D");
  //redoJER("Run2017");
  //redoJER("Run2016");
  //redoJER("Run1");

  //redoECALprefire(2.0,run2016);
  //redoECALprefire(2.5,run2016);
  //redoECALprefire(2.0,run2016bcd);
  //redoECALprefire(2.5,run2016bcd);
  //redoECALprefire(2.0,run2016ef);
  //redoECALprefire(2.5,run2016ef);
  //redoECALprefire(2.0,run2016gh);
  //redoECALprefire(2.5,run2016gh);

  //redoECALprefire(2.0,run2017);
  //redoECALprefire(2.5,run2017);
  //redoECALprefire(2.0,run2017b);
  //redoECALprefire(2.5,run2017b);
  //redoECALprefire(2.0,run2017c);
  //redoECALprefire(2.5,run2017c);
  //redoECALprefire(2.0,run2017de);
  //redoECALprefire(2.5,run2017de);
  //redoECALprefire(2.0,run2017f);
  //redoECALprefire(2.5,run2017f);

}
