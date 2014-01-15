{
  // Example code for accessing JEC uncertainties standalone
  // Download code with 'git clone http://www.github.com/miquork/jecsys'
  // Execute with 'root -l -b -q mk_testSources.C'
  
  // Compile stand-alone JEC libraries included in the package
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+");

  // Compile test code
  gROOT->ProcessLine(".L testSources.C+");

  // Test quadratic sum of sources in file1 against total in file2 at pT, eta
  testSources("txt/Summer13_V2_DATA_UncertaintySources_AK5PF.txt",
	      "txt/Summer13_V2_DATA_Uncertainty_AK5PF.txt", 50., 0.0);
  // NB1: Above example files are provided in the package
  // NB2: More can be produced with mk_drawJetCorrectionUncertainty.C

}
