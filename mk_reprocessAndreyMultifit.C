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


  //Quick recipe to make the new multijet work:
  //setup the jecsys as usual. You are already on the right branch if you can read this (2016LegacyWithNewMultijet for integration)
  //-make sure you have a working ROOT version with C++14 support enabled (see below for instructions)
  //-fetch Andrey's package (you can do this directly in the jecsys folder):
  /*
    git clone git@github.com:IPNL-CMS/jec-fit-prototype.git
    mkdir jec-fit-prototype/build
    cd jec-fit-prototype/build
    cmake ..
    make
    cd ..
    #reading root file directly from https only works if ROOT was compiled with -Dssl=ON (see below), otherwise test with local file
    #inputdir="https://aapopov.web.cern.ch/aapopov/jec_inputs/prototype"
    #bin/fit --photonjet-run1 $inputdir/photonjet_Run1.root --multijet-binnedsum $inputdir/multijet_BinnedSum.root --zjet-run1 $inputdir/Zjet_Run1.root
    #run locally instead
    bin/fit --zjet-run1 ../rootfiles/Zjet_Run1.root --multijet-binnedsum ../rootfiles/multijet_BinnedSum.root
   */
  //-go back to jecsys
  //run root -b -q 'mk_reprocessAndreyMultifit.C'
  //if there are any error messages, you should consider to turn the '+' to '++' at least temporary to make sure everything is recompiled with the fixed ROOT version
  

  //~Universal recipe to get a compatible root version compiled (standard ROOT releases don't seem to have C++14 activated by default, yet, at least not for Mac and Ubuntu16:)
  /*
  git clone http://root.cern.ch/git/root.git
  mkdir root_build
  cd root_build
  cmake  -Dcxx14=ON -Dminuit2=ON ../root
  //precompiled Mac version with also -Dssl=ON available from https://cernbox.cern.ch/index.php/s/ZG7V7L3FivXyUWS
  //standard Mac environment does not have openssl anymore, can do 
  //brew install openssl 
  //and then do 
  //export OPENSSL_ROOT_DIR=/usr/local/Cellar/openssl/1.0.2l
  //before compiling root for a "one-time-fix"
  cmake --build .
  */ 

  //to read in Andrey's global fit library [for multijet inclusion] - edit path by hand for now (needs to be absolute?)
  //  gSystem->AddIncludePath("-I/Users/kirschen/cernbox/JERC/jec-fit-prototype/include");
  //  gSystem->Load("/Users/kirschen/cernbox/JERC/jec-fit-prototype/lib/libjecfit.so");

  //if you have a symbolic link in jecsys to jec-fit-prototype, you can just do
  string currentWorkingDir = gSystem->pwd();
  cout <<currentWorkingDir.c_str() <<endl;
  gSystem->AddIncludePath(Form("-I%s/jec-fit-prototype/include",currentWorkingDir.c_str()));
  gSystem->Load(Form("%s/jec-fit-prototype/lib/libjecfit",currentWorkingDir.c_str()));

  // Compile with +g to make sure asserts are run
  gROOT->ProcessLine(".L tools.C+g");
  gROOT->ProcessLine(".L reprocess.C+g");
  gROOT->ProcessLine(".L softrad.C+g");
  //gROOT->ProcessLine(".L multijet.C+g");
  gROOT->ProcessLine(".L globalFitL3Res.C+g");

  // Merge inputs from separate groups
  // NB: this does not need to be run, if the merged inputs
  //     are already available in 'rootfiles/jecdata.root'
  string epoch = "GH";//"BCDEFGH";//"BCDEFGH";
  #ifdef epochname
  std::cout << epoch.c_str()<< std::endl;
  std::cout << inputepoch.c_str()<< std::endl;
  epoch = inputepoch;
  #endif
  // 0.8->0.9% 52.6->54.7, 53.5->60.8, 39.4->54.4, 71.8->72.5
  //"BCD", "EF", "G", "H", "BCDEFGH", "L4" (closure for |eta|<2.4)
  // BCD 47->46.9, EF 48.4->47.9, G 33.5->33.5, H 50.5

  reprocess(epoch); // Switched off for JetMET100
//
//  //softrad(0.0,epoch=="L4" ? 2.4 : 1.3,true,epoch); // redo for plots
  softrad(0.0,epoch=="L4" ? 2.4 : 1.3,false,epoch); // without dijets
//  // Run multijet analysis to store information for later global fit
//  // => multijet central values now old, but FSR still needed
//  multijet(false,epoch);
//  multijet(true,epoch);
//  // Perform final global fit (goes into GT)
  globalFitL3Res(0.0,epoch=="L4" ? 2.4 : 1.3, epoch); // L3Res

  //now do narrow bins for L2Res
  // Calculate soft radiation (ISR+FSR) corrections
  // and uncertainty eigenvectors for global fit

//  Bool_t dodijetsoftrad=true;
//
//  softrad(0.000,0.261, dodijetsoftrad, epoch); 
//  softrad(0.261,0.522, dodijetsoftrad, epoch); 
//  softrad(0.522,0.783, dodijetsoftrad, epoch); 
//  softrad(0.783,1.044, dodijetsoftrad, epoch); 
//  softrad(1.044,1.305, dodijetsoftrad, epoch); 
//  softrad(1.305,1.479, dodijetsoftrad, epoch); 
//  softrad(1.479,1.653, dodijetsoftrad, epoch); 
//  softrad(1.653,1.930, dodijetsoftrad, epoch); 
//  softrad(1.930,2.172, dodijetsoftrad, epoch); 
//  softrad(2.172,2.322, dodijetsoftrad, epoch); 
//  softrad(2.322,2.500, dodijetsoftrad, epoch); 
//  softrad(2.500,2.650, dodijetsoftrad, epoch); 
//  softrad(2.650,2.853, dodijetsoftrad, epoch); 
//  softrad(2.853,2.964, dodijetsoftrad, epoch); 
//  softrad(2.964,3.139, dodijetsoftrad, epoch); 
//  softrad(3.139,3.489, dodijetsoftrad, epoch); 
//  softrad(3.489,3.839, dodijetsoftrad, epoch); 
//  softrad(3.839,5.191, dodijetsoftrad, epoch);
//
//  globalFitL3Res(0.000,0.261, epoch); 
//  globalFitL3Res(0.261,0.522, epoch); 
//  globalFitL3Res(0.522,0.783, epoch); 
//  globalFitL3Res(0.783,1.044, epoch); 
//  globalFitL3Res(1.044,1.305, epoch); 
//  globalFitL3Res(1.305,1.479, epoch); 
//  globalFitL3Res(1.479,1.653, epoch); 
//  globalFitL3Res(1.653,1.930, epoch); 
//  globalFitL3Res(1.930,2.172, epoch); 
//  globalFitL3Res(2.172,2.322, epoch); 
//  globalFitL3Res(2.322,2.500, epoch); 
//  globalFitL3Res(2.500,2.650, epoch); 
//  globalFitL3Res(2.650,2.853, epoch); 
//  globalFitL3Res(2.853,2.964, epoch); 
//  globalFitL3Res(2.964,3.139, epoch); 
//  globalFitL3Res(3.139,3.489, epoch); 
//  globalFitL3Res(3.489,3.839, epoch); 
//  globalFitL3Res(3.839,5.191, epoch);
//
//  //produce summary pdf with all plots according to era
//  gSystem->Exec(Form("pdflatex '\\def\\RunPeriod{pdf/%s}\\input{pdf/jecslides_FineEta_2016Legacy.tex}'", epoch.c_str()));
//  gSystem->Exec("mkdir CollectL2Output");
//  gSystem->Exec(Form("mv jecslides_FineEta_2016Legacy.tex.pdf CollectL2Output/jecslides_FineEta_2016_%s.pdf", epoch.c_str()));
//  gSystem->Exec(Form("./minitools/convertGlobalFitOutputToStandardTxt.sh txt2/GlobalFitOutput_L2L3Residuals.txt  CollectL2Output/Summer16Legacy%s_VXXX_DATA_L2L3Residual_AK4PFchs.txt", epoch.c_str()));
//  gSystem->Exec("rm txt2/*");
//  
//  //wide eta bins
//   // softrad(0.0,0.8,dodijetsoftrad,epoch); // missing dijet
//   // softrad(0.8,1.3,dodijetsoftrad,epoch); // missing dijet
//   softrad(1.3,1.9,dodijetsoftrad,epoch);
//   softrad(1.9,2.5,dodijetsoftrad,epoch);
//   softrad(2.5,3.0,dodijetsoftrad,epoch);
//   softrad(3.0,3.2,dodijetsoftrad,epoch);
//   softrad(3.2,5.2,dodijetsoftrad,epoch);
//   softrad(0.0,1.3,dodijetsoftrad,epoch);
//   //These are just checks for now:
//   // globalFitL3Res(0.0,0.8,epoch); // coarse L2Res, missing dijet
//   // globalFitL3Res(0.8,1.3,epoch); // coarse L2Res, missing dijet
//   globalFitL3Res(1.3,1.9,epoch); // coarse L2Res
//   globalFitL3Res(1.9,2.5,epoch); // coarse L2Res
//   globalFitL3Res(2.5,3.0,epoch); // coarse L2Res
//   globalFitL3Res(3.0,3.2,epoch); // coarse L2Res
//   globalFitL3Res(3.2,5.2,epoch); // coarse L2Res
//
//   // Repeat to see parameters
//  globalFitL3Res(0.0,epoch=="L4" ? 2.4 : 1.3, epoch); // L3Res

}
