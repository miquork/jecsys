{
  gROOT->ProcessLine(".L compareJECdata.C+g");
  //gROOT->ProcessLine(".L compareJECdata.C");

  string dir = "mc";
  string ch = "gamjet";
  string var = "mpfchs1";
  //
  string alpha = "a30";
  string eta = "eta00-13";
  string file1 = "rootfiles/jecdataG_80XSum16V5.root";
  string label1 = "Gv5";
  string file2 = "rootfiles/jecdataG_80XSum16V4.root";
  string label2 = "Gv4";

  //compareJECdata("mc", "gamjet", "mpfchs1",
  //		 alpha,eta,file1,label1,file2,label2);

  //compareJECdata("mc", "zmmjet", "mpfchs1");
  //compareJECdata("data", "zmmjet", "mpfchs1");
  // => Zmm+jet now completely different for MPF, but for both data and MC
  //compareJECdata("mc", "zeejet", "mpfchs1");
  //compareJECdata("data", "zeejet", "mpfchs1");
  // => Zee+jet shows similar behavior. The |eta|<2.4 issue here?
  //compareJECdata("mc", "zmmjet", "ptchs");
  //compareJECdata("data", "zmmjet", "ptchs");
  // => Maybe bit less effect for pTbal, but less similar for data and MC

  //compareJECdata("mc", "gamjet", "mpfchs1"); // => stable
  //compareJECdata("data", "gamjet", "mpfchs1"); // => different slopes
  //compareJECdata("ratio", "gamjet", "mpfchs1"); // => different slopes
  //compareJECdata("mc", "gamjet", "ptchs"); // => stable, except pT<60 GeV
  //compareJECdata("data", "gamjet", "ptchs"); // => different slopes
  compareJECdata("ratio", "gamjet", "ptchs");
  // // => G+jet MC is stable (except pT<70 GeV for ptchs), but data changes slope
  // // => Need to fix the a10,a15,a20 being same as a30 issue, though

  // compareJECdata("mc", "zmmjet", "mpfchs1"); // => stable
  // compareJECdata("data", "zmmjet", "mpfchs1"); // => different slopes
  //compareJECdata("ratio", "zmmjet", "mpfchs1"); // => different slopes
  // compareJECdata("mc", "zmmjet", "ptchs"); // => stable, except pT<60 GeV
  // compareJECdata("data", "zmmjet", "ptchs"); // => different slopes
  // compareJECdata("ratio", "zmmjet", "ptchs");

  // compareJECdata("mc", "zeejet", "mpfchs1"); // => stable
  // compareJECdata("data", "zeejet", "mpfchs1"); // => different slopes
  //compareJECdata("ratio", "zeejet", "mpfchs1"); // => different slopes
  // compareJECdata("mc", "zeejet", "ptchs"); // => stable, except pT<60 GeV
  // compareJECdata("data", "zeejet", "ptchs"); // => different slopes
  // compareJECdata("ratio", "zeejet", "ptchs");

  // compareJECdata("mc", "multijet", "mpfchs1"); // => V4 dip at 650 GeV
  // compareJECdata("mc", "multijet", "ptchs"); // => V4 dip at 650 GeV
  // compareJECdata("data", "multijet", "mpfchs1"); // => V5 -0.4% @ <500-700 GeV
  // compareJECdata("data", "multijet", "ptchs"); // => V5 -0.6% @ <500-700 GeV
  // compareJECdata("ratio", "multijet", "mpfchs1"); // => crash?
  // compareJECdata("ratio", "multijet", "ptchs"); // => crash?
  // => Multijet has substantial changes in both data and MC
  // => Net effect should be significant drop for data/MC at low pT
}
