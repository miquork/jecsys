{

  // For JEC (for uncertainty)
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+");

  // For JEC uncertainty
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+");

  // For JER (used in ptresolution.C)
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetResolutionObject.cc+g");
  gROOT->ProcessLine(".L JetMETCorrections/Modules/src/JetResolution.cc+g");

  //gROOT->ProcessLine(".L minitools/ptresolution.h+g");
  gROOT->ProcessLine(".L minitools/drawDeltaJEC.C+g");

  const bool y = true;
  const bool n = false;
  bool fwd = y; // forward smearing for unfolding
  bool use13 = y;
  bool doAll = n;//y;//n;
  bool run2 = n;//y;
  bool ul2017 = n;//y;
  bool ul17x = n;//y;
  bool ul18x = y;
  bool eoy2016 = n;
  bool eoy2017 = n;
  bool eoy2018 = n;

  if (fwd) {
    if (use13) unfold("0.0-1.3",0);
    if (doAll) {
      unfold("0.0-0.5",1);
      unfold("0.5-1.0",2);
      unfold("1.0-1.5",3);
      unfold("1.5-2.0",4);
      unfold("2.0-2.5",5);
      unfold("2.5-3.0",6);
      unfold("3.2-4.7",7);
    }
  }
  if (run2) {
    if (use13) drawDeltaJEC("0.0-1.3");
    if (doAll) {
      drawDeltaJEC("0.0-0.5");
      drawDeltaJEC("0.5-1.0");
      drawDeltaJEC("1.0-1.5");
      drawDeltaJEC("1.5-2.0");
      drawDeltaJEC("2.0-2.5");
      drawDeltaJEC("2.5-3.0");
      drawDeltaJEC("3.2-4.7");
    }
  }
  if (ul2017) {
    if (use13)
      drawDeltaJEC("0.0-1.3","2017UL");
    if (doAll) {
      drawDeltaJEC("0.0-0.5","2017UL");
      drawDeltaJEC("0.5-1.0","2017UL");
      drawDeltaJEC("1.0-1.5","2017UL");
      drawDeltaJEC("1.5-2.0","2017UL");
      drawDeltaJEC("2.0-2.5","2017UL");
      drawDeltaJEC("2.5-3.0","2017UL");
      drawDeltaJEC("3.2-4.7","2017UL");
    }
  }
  if (ul17x) {
    if (use13) 
      drawDeltaJEC("0.0-1.3","17UL");
    if (doAll) {
      drawDeltaJEC("0.0-0.5","17UL");
      drawDeltaJEC("0.5-1.0","17UL");
      drawDeltaJEC("1.0-1.5","17UL");
      drawDeltaJEC("1.5-2.0","17UL");
      drawDeltaJEC("2.0-2.5","17UL");
      drawDeltaJEC("2.5-3.0","17UL");
      drawDeltaJEC("3.2-4.7","17UL");
    }
  }
  if (ul18x) {
    if (use13) 
      drawDeltaJEC("0.0-1.3","18UL");
    if (doAll) {
      drawDeltaJEC("0.0-0.5","18UL");
      drawDeltaJEC("0.5-1.0","18UL");
      drawDeltaJEC("1.0-1.5","18UL");
      drawDeltaJEC("1.5-2.0","18UL");
      drawDeltaJEC("2.0-2.5","18UL");
      drawDeltaJEC("2.5-3.0","18UL");
      drawDeltaJEC("3.2-4.7","18UL");
    }
  }
  /*
  //drawDeltaJEC("0.0-1.3","2017H");
  drawDeltaJEC("0.0-0.5","2017H");
  drawDeltaJEC("0.5-1.0","2017H");
  drawDeltaJEC("1.0-1.5","2017H");
  drawDeltaJEC("1.5-2.0","2017H");
  drawDeltaJEC("2.0-2.5","2017H");
  drawDeltaJEC("2.5-3.0","2017H");
  drawDeltaJEC("3.2-4.7","2017H");
  */
  if (eoy2018) {
    if (use13)
      drawDeltaJEC("0.0-1.3","2018");
    if (doAll) {
      drawDeltaJEC("0.0-0.5","2018");
      drawDeltaJEC("0.5-1.0","2018");
      drawDeltaJEC("1.0-1.5","2018");
      drawDeltaJEC("1.5-2.0","2018");
      drawDeltaJEC("2.0-2.5","2018");
      drawDeltaJEC("2.5-3.0","2018");  
      drawDeltaJEC("3.2-4.7","2018");
    }
  }
  if (eoy2017) {
    if (use13)
      drawDeltaJEC("0.0-1.3","2017");
    if (doAll) {
      drawDeltaJEC("0.0-0.5","2017");
      drawDeltaJEC("0.5-1.0","2017");
      drawDeltaJEC("1.0-1.5","2017");
      drawDeltaJEC("1.5-2.0","2017");
      drawDeltaJEC("2.0-2.5","2017");
      drawDeltaJEC("2.5-3.0","2017");
      drawDeltaJEC("3.2-4.7","2017");
    }
  }
  if (eoy2016) {
    if (use13)
      drawDeltaJEC("0.0-1.3","2016");
    if (doAll) {
      drawDeltaJEC("0.0-0.5","2016");
      drawDeltaJEC("0.5-1.0","2016");
      drawDeltaJEC("1.0-1.5","2016");
      drawDeltaJEC("1.5-2.0","2016");
      drawDeltaJEC("2.0-2.5","2016");
      drawDeltaJEC("2.5-3.0","2016");
      drawDeltaJEC("3.2-4.7","2016");
    }
  }
}
