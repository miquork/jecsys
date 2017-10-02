{
  // make sure asserts are run
  #undef NDEBUG
  // does not seem to be enough, also need to combile with +g

  // For JEC residual (and pile-up)
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
  gROOT->ProcessLine(".L CondFormats//JetMETObjects/src/JetCorrectorParameters.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+");
  //
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+");

  // Compile with +g to make sure asserts are run
  gROOT->ProcessLine(".L tools.C+g");
  gROOT->ProcessLine(".L reprocess.C+g");
  gROOT->ProcessLine(".L softrad.C+g");
  gROOT->ProcessLine(".L multijet.C+g");
  gROOT->ProcessLine(".L globalFitL3Res.C+g");

  // Merge inputs from separate groups
  // NB: this does not need to be run, if the merged inputs
  //     are already available in 'rootfiles/jecdata.root'
  string epoch = "H";
  //"BCD";//"E";//"F";//"EF";//"G";//"H";//"GH";//"BCDEF";//BCDEFGH
  //"L4"; // for BCDEFGH closure test including |eta|<2.4

  reprocess(epoch); // Switched off for JetMET100


//  // Calculate soft radiation (ISR+FSR) corrections
//  // and uncertainty eigenvectors for global fit
  softrad(0.0,epoch=="L4" ? 2.4 : 1.3,true,epoch);
//  //  softrad(0.0,0.8,true,epoch); // missing dijet
//  //  softrad(0.8,1.3,true,epoch); // missing dijet
//  softrad(1.3,1.9,true,epoch);
//  softrad(1.9,2.5,true,epoch);
//  softrad(2.5,3.0,true,epoch);
//  softrad(3.0,3.2,true,epoch);
//  softrad(3.2,5.2,true,epoch);
//  softrad(0.0,1.3,true,epoch);
  softrad(0.0,0.3,true,epoch);
  softrad(0.3,0.5,true,epoch);
  softrad(0.5,0.8,true,epoch);
  softrad(0.8,1.0,true,epoch);
  softrad(1.0,1.3,true,epoch);
  softrad(1.3,1.5,true,epoch);
  softrad(1.5,1.7,true,epoch);
  softrad(1.7,1.9,true,epoch);
  softrad(1.9,2.2,true,epoch);
  softrad(2.2,2.3,true,epoch);
  softrad(2.3,2.5,true,epoch);
  softrad(2.5,2.6,true,epoch);
  softrad(2.6,2.9,true,epoch);
  softrad(2.9,3.0,true,epoch);
  softrad(3.0,3.1,true,epoch);
  softrad(3.1,3.5,true,epoch);
  softrad(3.5,3.8,true,epoch);
  softrad(3.8,5.2,true,epoch);

  globalFitL3Res(0.0,epoch=="L4" ? 2.4 : 1.3, epoch);
  globalFitL3Res(0.0,0.3,epoch);
  globalFitL3Res(0.3,0.5,epoch);
  globalFitL3Res(0.5,0.8,epoch);
  globalFitL3Res(0.8,1.0,epoch);
  globalFitL3Res(1.0,1.3,epoch);
  globalFitL3Res(1.3,1.5,epoch);
  globalFitL3Res(1.5,1.7,epoch);
  globalFitL3Res(1.7,1.9,epoch);
  globalFitL3Res(1.9,2.2,epoch);
  globalFitL3Res(2.2,2.3,epoch);
  globalFitL3Res(2.3,2.5,epoch);
  globalFitL3Res(2.5,2.6,epoch);
  globalFitL3Res(2.6,2.9,epoch);
  globalFitL3Res(2.9,3.0,epoch);
  globalFitL3Res(3.0,3.1,epoch);
  globalFitL3Res(3.1,3.5,epoch);
  globalFitL3Res(3.5,3.8,epoch);
  globalFitL3Res(3.8,5.2,epoch);

  
//
//  //softrad(0.0,1.3,true,epoch); // redo for plots
//  softrad(0.0,epoch=="L4" ? 2.4 : 1.3,true,epoch); // redo for plots
//  // Run multijet analysis to store information for later global fit
//  // => multijet central values now old, but FSR still needed
//  multijet(false,epoch);
//  multijet(true,epoch);
//
//  // Perform final global fit (goes into GT)
//  //globalFitL3Res(0.0,1.3,epoch); // L3Res
//  // These are just checks for now:
//  //  globalFitL3Res(0.0,0.8,epoch); // coarse L2Res, missing dijet
//  //  globalFitL3Res(0.8,1.3,epoch); // coarse L2Res, missing dijet
//  globalFitL3Res(1.3,1.9,epoch); // coarse L2Res
//  globalFitL3Res(1.9,2.5,epoch); // coarse L2Res
//  globalFitL3Res(2.5,3.0,epoch); // coarse L2Res
//  globalFitL3Res(3.0,3.2,epoch); // coarse L2Res
//  globalFitL3Res(3.2,5.2,epoch); // coarse L2Res
//  // Redo for plots and fit results for L3Res
//  // (above eta bins may overwrite some plots)
//
//  globalFitL3Res(0.0,epoch=="L4" ? 2.4 : 1.3, epoch);
//
//  //globalFitL3Res(1.9,2.5,epoch);
//  //globalFitL3Res(2.5,3.0,epoch);

}
