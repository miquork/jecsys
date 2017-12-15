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
  //"BCD", "EF", "G", "H", "BCDEFGH", "L4" (closure for |eta|<2.4)
  // BCD 47->46.9, EF 48.4->47.9, G 33.5->33.5, H 50.5

  reprocess(epoch); // Switched off for JetMET100

  //softrad(0.0,epoch=="L4" ? 2.4 : 1.3,true,epoch); // redo for plots
  softrad(0.0,epoch=="L4" ? 2.4 : 1.3,false,epoch); // without dijets
  // Run multijet analysis to store information for later global fit
  // => multijet central values now old, but FSR still needed
  multijet(false,epoch);
  multijet(true,epoch);
  // Perform final global fit (goes into GT)
  globalFitL3Res(0.0,epoch=="L4" ? 2.4 : 1.3, epoch); // L3Res

  //now do narrow bins for L2Res
  // Calculate soft radiation (ISR+FSR) corrections
  // and uncertainty eigenvectors for global fit
  // softrad(0.000,0.261, true, epoch); 
  // softrad(0.261,0.522, true, epoch); 
  // softrad(0.522,0.783, true, epoch); 
  // softrad(0.783,1.044, true, epoch); 
  // softrad(1.044,1.305, true, epoch); 
  // softrad(1.305,1.479, true, epoch); 
  // softrad(1.479,1.653, true, epoch); 
  // softrad(1.653,1.930, true, epoch); 
  // softrad(1.930,2.172, true, epoch); 
  // softrad(2.172,2.322, true, epoch); 
  // softrad(2.322,2.500, true, epoch); 
  // softrad(2.500,2.650, true, epoch); 
  // softrad(2.650,2.853, true, epoch); 
  // softrad(2.853,2.964, true, epoch); 
  // softrad(2.964,3.139, true, epoch); 
  // softrad(3.139,3.489, true, epoch); 
  // softrad(3.489,3.839, true, epoch); 
  // softrad(3.839,5.191, true, epoch);

  // globalFitL3Res(0.000,0.261, epoch); 
  // globalFitL3Res(0.261,0.522, epoch); 
  // globalFitL3Res(0.522,0.783, epoch); 
  // globalFitL3Res(0.783,1.044, epoch); 
  // globalFitL3Res(1.044,1.305, epoch); 
  // globalFitL3Res(1.305,1.479, epoch); 
  // globalFitL3Res(1.479,1.653, epoch); 
  // globalFitL3Res(1.653,1.930, epoch); 
  // globalFitL3Res(1.930,2.172, epoch); 
  // globalFitL3Res(2.172,2.322, epoch); 
  // globalFitL3Res(2.322,2.500, epoch); 
  // globalFitL3Res(2.500,2.650, epoch); 
  // globalFitL3Res(2.650,2.853, epoch); 
  // globalFitL3Res(2.853,2.964, epoch); 
  // globalFitL3Res(2.964,3.139, epoch); 
  // globalFitL3Res(3.139,3.489, epoch); 
  // globalFitL3Res(3.489,3.839, epoch); 
  // globalFitL3Res(3.839,5.191, epoch);

  //wide eta bins
   // softrad(0.0,0.8,true,epoch); // missing dijet
   // softrad(0.8,1.3,true,epoch); // missing dijet
   softrad(1.3,1.9,true,epoch);
   softrad(1.9,2.5,true,epoch);
   softrad(2.5,3.0,true,epoch);
   softrad(3.0,3.2,true,epoch);
   softrad(3.2,5.2,true,epoch);
   softrad(0.0,1.3,true,epoch);
   //These are just checks for now:
   // globalFitL3Res(0.0,0.8,epoch); // coarse L2Res, missing dijet
   // globalFitL3Res(0.8,1.3,epoch); // coarse L2Res, missing dijet
   globalFitL3Res(1.3,1.9,epoch); // coarse L2Res
   globalFitL3Res(1.9,2.5,epoch); // coarse L2Res
   globalFitL3Res(2.5,3.0,epoch); // coarse L2Res
   globalFitL3Res(3.0,3.2,epoch); // coarse L2Res
   globalFitL3Res(3.2,5.2,epoch); // coarse L2Res



}
