{
  //gROOT->ProcessLine(".L alcajme.C"); // 16.24.04
  gROOT->ProcessLine(".L alcajme.C+g"); // 16.18.04
  gROOT->ProcessLine("alcajme t;");
  t.Loop();
}
