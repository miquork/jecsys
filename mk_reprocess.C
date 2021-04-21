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
  string currentWorkingDir = gSystem->pwd();
  cout <<currentWorkingDir.c_str() <<endl;
  gSystem->AddIncludePath(Form("-I%s/jec-fit-prototype/include",currentWorkingDir.c_str()));
  gSystem->Load(Form("%s/jec-fit-prototype/lib/libjecfit",currentWorkingDir.c_str()));

  
  // Compile with +g to make sure asserts are run
  gROOT->ProcessLine(".L tools.C+g");
  gROOT->ProcessLine(".L reprocess.C+g");
  gROOT->ProcessLine(".L softrad.C+g");
  gROOT->ProcessLine(".L softrad3.C+g");
  //gROOT->ProcessLine(".L multijet.C+g"); // obsolete since 20200318
  gROOT->ProcessLine(".L globalFitSyst.C+g");
  gROOT->ProcessLine(".L globalFitL3Res.C+g");

  // Merge inputs from separate groups
  // NB: this does not need to be run, if the merged inputs
  //     are already available in 'rootfiles/jecdata.root'
  string epoch = "A";//"BCDEFGH";
  #ifdef epochname
  std::cout << epoch.c_str()<< std::endl;
  std::cout << inputepoch.c_str()<< std::endl;
  epoch = inputepoch;
  #endif
  // 0.8->0.9% 52.6->54.7, 53.5->60.8, 39.4->54.4, 71.8->72.5
  //"BCD", "EF", "G", "H", "BCDEFGH", "L4" (closure for |eta|<2.4)
  // BCD 47->46.9, EF 48.4->47.9, G 33.5->33.5, H 50.5

  reprocess(epoch); // Switched off for JetMET100

  softrad(0.0,epoch=="L4" ? 2.4 : 1.3,true,epoch); // redo for plots
  //softrad3(0.0,epoch=="L4" ? 2.4 : 1.3,true,epoch); // 3-point FSR
  softrad3(0.0,1.3,true,epoch); // 3-point FSR
  softrad3(0.0,2.5,true,epoch); // 3-point FSR
  
  globalFitSyst(epoch);

  //  softrad(0.0,epoch=="L4" ? 2.4 : 1.3,false,epoch); // without dijets
  // Run multijet analysis to store information for later global fit
  // => multijet central values now old, but FSR still needed
  // 20200319: multijet had FSR=0, sys obsolete; do FSR now in softrad
  //multijet(false,epoch);
  //multijet(true,epoch);
  // Perform final global fit (goes into GT) - 2018v0
  //globalFitL3Res(0.0,epoch=="L4" ? 2.4 : 1.3, epoch, epoch=="B"||epoch=="ABC" ? "zee_zmm" : "gam_zee_zmm", "PtBalMPF");
  // Perform final global fit (goes into GT) - 2018 V5M
  //globalFitL3Res(0.0,epoch=="L4" ? 2.4 : 1.3, epoch, epoch=="B"||epoch=="ABC"||epoch=="D"||epoch=="ABCD" ? "zll" : "gam_zll", "PtBalMPF");
  //globalFitL3Res(0.0,epoch=="L4" ? 2.4 : 1.3, epoch, epoch=="ABC"||epoch=="ABCD" ? "zll" : "gam_zll", "PtBalMPF");
  //globalFitL3Res(0.0,epoch=="L4" ? 2.4 : 1.3, epoch, "gam_zll", "PtBalMPF");
  //globalFitL3Res(0.0,epoch=="L4" ? 2.4 : 1.3, epoch, "MJDJ_gam_zll", "PtBalMPF");
  if (epoch=="2018ABCD" || epoch=="2017BCDEF" || epoch=="2016BCDEF" || epoch=="2016GH") {
    //globalFitL3Res(0.0,epoch=="L4" ? 2.4 : 1.3, epoch, "gam_zll", "PtBalMPF");
    //globalFitL3Res(0.0,epoch=="L4" ? 2.4 : 1.3, epoch, "gam_zll_hadw", "PtBalMPF");
  //globalFitL3Res(0.0,epoch=="L4" ? 2.4 : 1.3, epoch, "MJDJ_gam_zll", "PtBalMPF");
    if (epoch=="2018ABCD")
      globalFitL3Res(0.0,epoch=="L4" ? 2.4 : 1.3, epoch, "MJDJ_gam_zll_hadw", "PtBalMPF"); // V4 (KIT)
    //globalFitL3Res(0.0,epoch=="L4" ? 2.4 : 1.3, epoch, "MJDJ_gam_z_hadw", "PtBalMPF"); // V4 (UH)
    if (epoch=="2017BCDEF")
      globalFitL3Res(0.0,epoch=="L4" ? 2.4 : 1.3, epoch, "MJDJ_gam_zll_hadw", "PtBalMPF"); // V5 (KIT)
      //globalFitL3Res(0.0,epoch=="L4" ? 2.4 : 1.3, epoch, "MJDJ_gam_z_hadw", "PtBalMPF"); // V4 (UH)
    if (epoch=="2016BCDEF")
      //globalFitL3Res(0.0,1.3, epoch, "MJDJ_inc_gam_zll_hadw", "PtBalMPF");
      globalFitL3Res(0.0,1.3, epoch, "MJDJ_gam_zll_hadw", "PtBalMPF");
      //globalFitL3Res(0.0,1.3, epoch, "MJDJ_gam_zmm", "PtBalMPF");
      //globalFitL3Res(0.0,1.3, epoch, "MJDJ_gam_zll", "PtBalMPF");
    if (epoch=="2016GH")
      //globalFitL3Res(0.0,1.3, epoch, "MJDJ_gam_zll_hadw", "MPF"); // not working
      //globalFitL3Res(0.0,1.3, epoch, "MJDJ_gam_zll_hadw", "PtBalMPF");
      globalFitL3Res(0.0,1.3, epoch, "MJDJ_inc_gam_zll_hadw", "PtBalMPF");
      //globalFitL3Res(0.0,1.3, epoch, "MJDJ_gam_zll", "PtBalMPF");
    //globalFitL3Res(0.0,1.3, epoch, "MJDJ_zll", "PtBalMPF");
  }
  else if (epoch=="2016BCD" || epoch=="2016EF")
    //globalFitL3Res(0.0,1.3, epoch, "MJDJ_gam_zll_hadw", "PtBalMPF");
    globalFitL3Res(0.0,1.3, epoch, "MJDJ_inc_gam_zll_hadw", "PtBalMPF");
  else if (epoch=="BCDEF")
    //globalFitL3Res(0.0,epoch=="L4" ? 2.4 : 1.3, epoch, "MJDJ_gam_zll", "PtBalMPF");
    //globalFitL3Res(0.0,epoch=="L4" ? 2.4 : 1.3, epoch, "MJDJ_gam_z", "PtBalMPF");
    globalFitL3Res(0.0,epoch=="L4" ? 2.4 : 1.3, epoch, "MJDJ_gam_zll_hadw", "PtBalMPF"); // hadW
    //globalFitL3Res(0.0,epoch=="L4" ? 2.4 : 1.3, epoch, "MJDJ_gam_z_hadw", "PtBalMPF"); // hadW
  else if (epoch=="2018A" || epoch=="2018B" ||
	   epoch=="2018C" || epoch=="2018D")
    //globalFitL3Res(0.0,epoch=="L4" ? 2.4 : 1.3, epoch, "MJDJ_gam_zll", "PtBalMPF");
    //globalFitL3Res(0.0,epoch=="L4" ? 2.4 : 1.3, epoch, "MJDJ_gam_zll_hadw", "PtBalMPF");
    globalFitL3Res(0.0,epoch=="L4" ? 2.4 : 1.3, epoch, "MJDJ_inc_gam_zll_hadw", "PtBalMPF");
  else
    globalFitL3Res(0.0,epoch=="L4" ? 2.4 : 1.3, epoch, "MJDJ_inc_gam_zll", "PtBalMPF");
  //globalFitL3Res(0.0,epoch=="L4" ? 2.4 : 1.3, epoch, "MJDJ_zll", "PtBalMPF");
  //globalFitL3Res(0.0,epoch=="L4" ? 2.4 : 1.3, epoch, "zll", "PtBalMPF");

  
////  
////  //now do narrow bins for L2Res
////  // Calculate soft radiation (ISR+FSR) corrections
////  // and uncertainty eigenvectors for global fit
////
//  Bool_t dodijetsoftrad=true;
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
//
//  //std::vector <string> sampleconfigs = {epoch=="B"||epoch=="ABC"||epoch=="ABCD" ? "zee_zmm" : "gam_zee_zmm"};//"MJDJ_gam_zll","DJ","gam","zll","gam_zll"};
//  //  std::vector <string> sampleconfigs = {epoch=="ABC"||epoch=="ABCD" ? "zll" : "gam_zll"};
//  //  std::vector <string> sampleconfigs = {epoch=="D"||epoch=="ABCD" ? "zll" : "gam_zll"};
//  //  std::vector <string> sampleconfigs = {"gam_zll"};
//  //std::vector <string> sampleconfigs = {"zmm"};//"MJDJ_gam_zll","DJ","gam","zll","gam_zll"};
//  std::vector <string> sampleconfigs = {"gam_zll","MJDJ_gam_zll","DJ"};//,"gam","zll","gam_zll"};
//  //std::vector <string> sampleconfigs = {"MJDJ_gam_zll"};//,"DJ"};//,"gam","zll","gam_zll"};
//  std::vector <string> methodconfigs = {  "PtBalMPF"};//,"PtBal","MPF"};
//
//    for(auto s : sampleconfigs){
//    for(auto m : methodconfigs){
//      globalFitL3Res(0.000,0.261, epoch, s, m); //default: Standard_MJDJ_gam_zee_zmm; PtBalMPF
//      globalFitL3Res(0.261,0.522, epoch, s, m); 
//      globalFitL3Res(0.522,0.783, epoch, s, m); 
//      globalFitL3Res(0.783,1.044, epoch, s, m); 
//      globalFitL3Res(1.044,1.305, epoch, s, m); 
//      globalFitL3Res(1.305,1.479, epoch, s, m); 
//      globalFitL3Res(1.479,1.653, epoch, s, m); 
//      globalFitL3Res(1.653,1.930, epoch, s, m); 
//      globalFitL3Res(1.930,2.172, epoch, s, m); 
//      globalFitL3Res(2.172,2.322, epoch, s, m); 
//      globalFitL3Res(2.322,2.500, epoch, s, m); 
//      globalFitL3Res(2.500,2.650, epoch, s, m); 
//      globalFitL3Res(2.650,2.853, epoch, s, m); 
//      globalFitL3Res(2.853,2.964, epoch, s, m); 
//      globalFitL3Res(2.964,3.139, epoch, s, m); 
//      globalFitL3Res(3.139,3.489, epoch, s, m); 
//      globalFitL3Res(3.489,3.839, epoch, s, m); 
//      globalFitL3Res(3.839,5.191, epoch, s, m);
//      
//      //produce summary pdf with all plots according to era
//      gSystem->Exec(Form("pdflatex '\\def\\RunPeriod{pdf/%s}\\input{pdf/jecslides_FineEta_2016Legacy.tex}'", epoch.c_str()));
//      string FolderName = Form("CollectL2Output_%s_%s",s.c_str(),m.c_str());
//      gSystem->Exec(Form("mkdir %s",FolderName.c_str()));
//      gSystem->Exec(Form("mv jecslides_FineEta_2016Legacy.pdf %s/jecslides_FineEta_Autumn18_V16_%s.pdf", FolderName.c_str(), epoch.c_str()));
//      gSystem->Exec(Form("./minitools/convertGlobalFitOutputToStandardTxt.sh txt2/GlobalFitOutput_L2L3Residuals.txt  %s/Summer16Legacy%s_VXXX_DATA_L2L3Residual_AK4PFchs.txt", FolderName.c_str(), epoch.c_str()));
//      gSystem->Exec(Form("./minitools/convertGlobalFitOutputToStandardTxt.sh txt2/GlobalFitOutput_L2L3Residuals_Chi2OverNDF.txt  %s/Summer16Legacy%s_VXXX_DATA_L2L3Residual_Chi2OverNDF_AK4PFchs.txt", FolderName.c_str(), epoch.c_str()));
//      gSystem->Exec(Form("cp  %s/Summer16Legacy%s_VXXX_DATA_L2L3Residual_AK4PFchs.txt %s_Summer16Legacy%s_VXXX_DATA_L2L3Residual_AK4PFchs.txt", FolderName.c_str(), epoch.c_str(), FolderName.c_str(), epoch.c_str()));
//      gSystem->Exec("rm txt2/*");
//  
//      
//    }
//  }
//
//


  
//  //wide eta bins
//   // softrad(0.0,0.8,true,epoch); // missing dijet
//   // softrad(0.8,1.3,true,epoch); // missing dijet
//   softrad(1.3,1.9,true,epoch);
//   softrad(1.9,2.5,true,epoch);
//   softrad(2.5,3.0,true,epoch);
//   softrad(3.0,3.2,true,epoch);
//   softrad(3.2,5.2,true,epoch);
//   softrad(0.0,1.3,true,epoch);
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
//  //  globalFitL3Res(0.0,epoch=="L4" ? 2.4 : 1.3, epoch); // L3Res
  
}
