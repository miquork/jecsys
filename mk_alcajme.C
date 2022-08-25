{
  //gROOT->ProcessLine(".L alcajme.C"); // 16.24.04
  gROOT->ProcessLine(".L alcajme.C+g"); // 16.18.04
  //gROOT->ProcessLine("alcajme t;");
  gROOT->ProcessLine("alcajme t");
  t.Loop();

  // README / CHANGELOG
  // v10 - 4.86/fb golden JSON (up to 357550 vs 357112, 09:00->09:44 => 44 min)
  // Found 11097277 bad events according to new JSON (events cut)
  // Found 21218842 bad events according to new AK8 vector (events cut)
  // Processed 36 runs, 7805 luminosity blocks and 39415771 events 
  // 356482 -> 357112, 1.79/fb out of 4.86/fb
  //
  // v9 - 1.44/fb golden JSON (up to 356616, 08:27->08:33 => 6 min)
  // v8 - veto nAK8!=0 to get only AK4 path, doAK4 (08:53->09:54 => 1h 1 min)
  // v7 - add ptbins as separate folders (23:16->01:17 => 2h 1min)
  // v6 - add h2dphi and hmjj histograms (multiple variants)
  // processing time 14:45-12:36 => 2h 9min
}
