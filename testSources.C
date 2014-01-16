#include "TFile.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include <iostream>
#include <vector>
#include <string>

using namespace std;

bool testSources(string file1, string file2,
		 double jetpt=50., double jeteta=2.4) {
  
  // Instantiate uncertainty sources
  const char* srcnames[] =
    {"Absolute", "HighPtExtra", "SinglePionECAL", "SinglePionHCAL",
     "FlavorQCD", "Time",
     "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF",
     "RelativePtBB","RelativePtEC1", "RelativePtEC2", "RelativePtHF",
     "RelativeFSR", "RelativeStatEC2", "RelativeStatHF",
     "PileUpDataMC",
     "PileUpPtBB", "PileUpPtEC", "PileUpPtHF",
     "Total"};
  // Subtotals and flavor uncertainties
  //"PileUpBias"};
  //"SubTotalPileUp","SubTotalRelative","SubTotalPt","SubTotalMC",
  //"Total","TotalNoFlavor",
  //"FlavorZJet","FlavorPhotonJet","FlavorPureGluon","FlavorPureQuark",
  //"FlavorPureCharm","FlavorPureBottom"};

  // Can use one of these to test subtotals. In this case, comment out above
  //const char* srcnames[] =
  //{"HighPtExtra", "SinglePionECAL", "SinglePionHCAL","SubTotalPt"};
  //{"RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF", "RelativePtBB","RelativePtEC1", "RelativePtEC2", "RelativePtHF", "RelativeFSR", "RelativeStatEC2", "RelativeStatHF", "SubTotalRelative"};
  //{"PileUpDataMC","PileUpPtBB", "PileUpPtEC", "PileUpPtHF","SubTotalPileUp"};
  const int nsrc = sizeof(srcnames)/sizeof(char*)-1;
  std::vector<JetCorrectionUncertainty*> vsrc(nsrc);
  
  // Load individual sources from file1
  cout <<"\nLoading individual sources ("<<nsrc<<") from file1..." << endl;
  for (int isrc = 0; isrc < nsrc; isrc++) {
    
    const char *name = srcnames[isrc];
    cout << Form("=> %s:%s",file1.c_str(),name) << endl << flush;
    JetCorrectorParameters *p = new JetCorrectorParameters(file1.c_str(), name);
    JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
    vsrc[isrc] = unc;
  } // for isrc
  
  // Load total uncertainty source from file1
  cout << "\nLoading Total source from file1..." << endl;
  cout << Form("=> %s:%s",file1.c_str(),srcnames[nsrc]) << endl << flush;
  JetCorrectorParameters *p = new JetCorrectorParameters(file1.c_str(),srcnames[nsrc]);
  JetCorrectionUncertainty *total = new JetCorrectionUncertainty(*p);
  
  // Load total uncertainty from file2
  cout << "\nLoading total uncertainty from file2..." << endl;
  cout << "=> " << file2 << endl << flush;
  JetCorrectorParameters *p2 = new JetCorrectorParameters(file2.c_str());
  JetCorrectionUncertainty *total2 = new JetCorrectionUncertainty(*p2);
  
  // Calculate uncertainty per source and as a total
  double sum2_up(0), sum2_dw(0), sum2(0);
  for (int isrc = 0; isrc != nsrc; isrc++) {
    
    JetCorrectionUncertainty *unc = vsrc[isrc];
    unc->setJetPt(jetpt);
    unc->setJetEta(jeteta);
    double sup = unc->getUncertainty(true); // up variation
    unc->setJetPt(jetpt);
    unc->setJetEta(jeteta);
    double sdw = unc->getUncertainty(false); // down variation
    
    sum2_up += pow(max(sup,sdw),2);
    sum2_dw += pow(min(sup,sdw),2);
    sum2 += pow(max(fabs(sup),fabs(sdw)),2);
  } // for isrc

  // Check that uncertainties are symmetric
  assert(fabs(sum2_up-sum2_dw)<5e-4);
  assert(fabs(sum2-fabs(sum2_up))<5e-4);
  
  total->setJetPt(jetpt);
  total->setJetEta(jeteta);
  double uncert = total->getUncertainty(true);
  
  total2->setJetPt(jetpt);
  total2->setJetEta(jeteta);
  double uncert2 = total2->getUncertainty(true);
  
  bool fail = false;

  // Check that quadratic sum of sources equals total uncertainty
  cout << "\nTesting sum versus file1 total...";
  if (fail = (fabs(uncert - sqrt(sum2)) > 5e-4)) {
    cout << "FAIL" << endl;
    cout << "Error: Quadratic sum of sources does not agree with total" << endl
	 << "       Difference is "<<uncert<<"(total) - "<<sqrt(sum2)
	 << "(sum) = " << uncert - sqrt(sum2) << endl;
    //assert(fabs(uncert - sqrt(sum2_up)) < 1e-3);
  }
  else
    cout << "ok" << endl;
  
  cout << "Testing sum versus file2 total...";
  if (fail = (fabs(uncert2 - sqrt(sum2)) > 5e-4)) {
    cout << "FAIL" << endl;
    cout << "Error: Quadratic sum of sources does not agree with total" << endl
	 << "       Difference is "<<uncert2<<"(total) - "<<sqrt(sum2)
	 << "(sum) = " << uncert2 - sqrt(sum2) << endl;
    //assert(fabs(uncert2 - sqrt(sum2_up)) < 1e-3);
  }
  else
    cout << "ok" << endl;

  cout << "Test result: " << (fail ? "FAIL" : "PASS") << endl;
  cout << "Tested eta="<<jeteta<<" pt="<<jetpt<<endl;

  return (!fail);
} // testSources
