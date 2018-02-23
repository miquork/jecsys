void mk_reprocess_RunEpoch(std::string inepoch=""){
  #define epochname
  gROOT->ProcessLine(Form("string inputepoch(\"%s\");",inepoch.c_str()));
  gROOT->ProcessLine(".x mk_reprocess.C");
}
