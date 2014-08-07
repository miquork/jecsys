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
  testSources("txt/Winter14_V5_DATA_UncertaintySources_AK5PF.txt",
	      "txt/Winter14_V5_DATA_Uncertainty_AK5PF.txt", 49.25, 2.0,"Winter14_V5");
  // NB1: Above example files are provided in the package
  // NB2: More can be produced with mk_drawJetCorrectionUncertainty.C
  testSourcesLoopVars("txt/Winter14_V5_DATA_UncertaintySources_AK5PF.txt",
		      "txt/Winter14_V5_DATA_Uncertainty_AK5PF.txt","Winter14_V5");//check that sum of sources equals total in file1 and file2 for a number of eta/pt-combinations

//   testCommonSources  ("txt/Summer13_V5_DATA_UncertaintySources_AK5PF.txt",
// 		      "txt/Summer13_V2_DATA_UncertaintySources_AK5PF.txt", 49.25, 2.0,"Summer13_V2V5CommonSources"); // check individual common sources for backward-compatibility
//   testSourcesLoopVars("txt/Winter14_V5_DATA_UncertaintySources_AK5PF.txt",
// 		      "txt/Winter14_V5_DATA_Uncertainty_AK5PF.txt","Winter14_V5"
// 		      ); // check for backward-compatibility
  testSourcesLoopVars("txt/Winter14_V5_DATA_UncertaintySources_AK5PF.txt",
		      "","Winter14_V5_SubTotalPt"); // SubTotalPt
  testSourcesLoopVars("txt/Winter14_V5_DATA_UncertaintySources_AK5PF.txt",
		      "","Winter14_V5_SubTotalRelative"); // SubTotalRelative
  testSourcesLoopVars("txt/Winter14_V5_DATA_UncertaintySources_AK5PF.txt",
		      "","Winter14_V5_SubTotalPileUp"); // SubTotalPileUp
  testSourcesLoopVars("txt/Winter14_V5_DATA_UncertaintySources_AK5PF.txt",
		      "","Winter14_V5_CorrGroups"); // CorrGroups
  testSourcesLoopVars("txt/Winter14_V5_DATA_UncertaintySources_AK5PF.txt",
		      "txt/Winter14_V5_DATA_Uncertainty_AK5PF.txt","Winter14_V5");//check that sum of sources equals total in file1 and file2 for a number of eta/pt-combinations


  cout << "NB: Only basic tets implemented yet, skipping rest..." << endl;
  exit();



  if(nFailedTests)std::cout <<  "WARNING: There have been " << nFailedTests << " failed tests "  << std::endl;
}
