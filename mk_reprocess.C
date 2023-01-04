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

  //if you have a symbolic link in jecsys to jec-fit-prototype, you can just do
  //string currentWorkingDir = gSystem->pwd();
  //cout <<currentWorkingDir.c_str() <<endl;
  //gSystem->AddIncludePath(Form("-I%s/jec-fit-prototype/include",currentWorkingDir.c_str()));
  //gSystem->Load(Form("%s/jec-fit-prototype/lib/libjecfit",currentWorkingDir.c_str()));

  
  // Compile with +g to make sure asserts are run
  gROOT->ProcessLine(".L tools.C+g");
  gROOT->ProcessLine(".L Flavor.C+g");
  gROOT->ProcessLine(".L reprocess.C+g");
  //gROOT->ProcessLine(".L softrad.C+g");
  gROOT->ProcessLine(".L softrad3.C+g");
  gROOT->ProcessLine(".L globalFitSyst.C+g");
  gROOT->ProcessLine(".L globalFitRenormPF.C+g");
  gROOT->ProcessLine(".L globalFitL3Res.C+g");

  // Merge inputs from separate groups
  // NB: this does not need to be run, if the merged inputs
  //     are already available in 'rootfiles/jecdata.root'
  string epoch = "Run2Test";
  #ifdef epochname
  std::cout << epoch.c_str()<< std::endl;
  std::cout << inputepoch.c_str()<< std::endl;
  epoch = inputepoch;
  #endif

  // Read in files from different groups and merge them in jecdata[epoch].root
  reprocess(epoch); // Comment out if using archived jecdata[epoch].root

  // Alpha extrapolation: use alpha<0.1,0.15,0.20,0.30 to derive FSR correction
  //softrad(0.0,1.3,true,epoch); // redo for plots

  // HDM method: use HT decomposition (lead, soft jets, unclustered) for FSR
  softrad3(0.0,1.3,true,epoch); // 3-point FSR
  softrad3(0.0,2.5,true,epoch); // 3-point FSR

  // Produce central systematic uncertainties for globalFitL3Res
  globalFitSyst(epoch);     // also for globalFitRun2.C
  globalFitRenormPF(epoch); // for globalFitRun2.C

  // Run global fit
  /////////////////

  // Reference IOVs: don't use inclusive jets
  if (epoch=="2018ABCD" || epoch=="2017BCDEF" || epoch=="2016BCDEF" ||
      epoch=="2016GH" || 
      epoch=="Run2Test") {
    //globalFitL3Res(0.0,1.3, epoch, "MJDJ_gam_zll_hadw", "PtBalMPF");
    globalFitL3Res(0.0,1.3, epoch, "MJDJ_gam_z_hadw", "PtBalMPF"); // DP_2021
  }
  // Other non-reference IOVs: add inclusive jets
  else if (epoch=="2016BCD" || epoch=="2016EF" || //epoch=="2016GH" ||
	   epoch=="2017B" || epoch=="2017C" || epoch=="2017D" ||
	   epoch=="2017E" || epoch=="2017F" ||
	   epoch=="2018A" || epoch=="2018B" || epoch=="2018C" ||
	   epoch=="2018D")
    globalFitL3Res(0.0,1.3, epoch, "MJDJ_inc_gam_zll_hadw", "PtBalMPF");
    //globalFitL3Res(0.0,1.3, epoch, "MJDJ_inc_gam_z_hadw", "PtBalMPF");
  // Non-supported options, bail out
  else if (epoch=="2017H") {
    globalFitL3Res(0.0,1.3, epoch, "inc_z", "PtBalMPF");
  }
  else
    assert(false);
  
}
