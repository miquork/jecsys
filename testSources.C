#include "TFile.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include <iostream>
#include <vector>
#include <string>

using namespace std;

// Turn debug on to produce more verbose output for debugging (e.g. listing individual sources)
bool debug = false;
bool verbose = true;

unsigned int nFailedTests = 0;

const char* srcnames_Summer13_V2[] =
  {"Absolute", "HighPtExtra", "SinglePionECAL", "SinglePionHCAL",
   "FlavorQCD", "Time",
   "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF",
   "RelativePtBB","RelativePtEC1", "RelativePtEC2", "RelativePtHF",
   "RelativeFSR", "RelativeStatEC2", "RelativeStatHF",
   "PileUpDataMC",
   "PileUpPtBB", "PileUpPtEC", "PileUpPtHF",
   "Total"};

const int nsrc_Summer13_V2 = sizeof(srcnames_Summer13_V2)/sizeof(char*)-1;

const char* srcnames_Summer13_V5[] =
  {"AbsoluteStat", "AbsoluteScale", "AbsoluteFlavMap", 
   "AbsoluteMPFBias", "HighPtExtra", "SinglePionECAL", "SinglePionHCAL",
   "FlavorQCD", "Time",
   "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF",
   "RelativePtBB","RelativePtEC1", "RelativePtEC2", "RelativePtHF",
   "RelativeFSR", "RelativeStatEC2", "RelativeStatHF",
   "PileUpDataMC",
   "PileUpPtBB", "PileUpPtEC", "PileUpPtHF",
   "Total"};

const int nsrc_Summer13_V5 = sizeof(srcnames_Summer13_V5)/sizeof(char*)-1;

const char* srcnames_Summer13_V2V5CommonSources[] =
  {"HighPtExtra", "SinglePionECAL", "SinglePionHCAL",
   "FlavorQCD", "Time",
   "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF",
   "RelativePtBB","RelativePtEC1", "RelativePtEC2", "RelativePtHF",
   "RelativeFSR", "RelativeStatEC2", "RelativeStatHF",
   "PileUpDataMC",
   "PileUpPtBB", "PileUpPtEC", "PileUpPtHF", "PileUpBias",
   "SubTotalPileUp", "SubTotalRelative", "SubTotalPt", "SubTotalMC",
   "Total", "TotalNoFlavor", "FlavorZJet", "FlavorPhotonJet", "FlavorPureGluon", 
   "FlavorPureQuark", "FlavorPureCharm", "FlavorPureBottom"};

const int nsrc_Summer13_V2V5CommonSources = sizeof(srcnames_Summer13_V2V5CommonSources)/sizeof(char*)-1;

const char* srcnames_Summer13_V5_SubTotalPt[] =
  {"HighPtExtra", "SinglePionECAL", "SinglePionHCAL","SubTotalPt"}; 
const int nsrc_Summer13_V5_SubTotalPt = sizeof(srcnames_Summer13_V5_SubTotalPt)/sizeof(char*)-1;

const char* srcnames_Summer13_V5_SubTotalRelative[] =
  {"RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF", "RelativePtBB","RelativePtEC1", "RelativePtEC2", "RelativePtHF", "RelativeFSR", "RelativeStatEC2", "RelativeStatHF", "SubTotalRelative"};
const int nsrc_Summer13_V5_SubTotalRelative = sizeof(srcnames_Summer13_V5_SubTotalRelative)/sizeof(char*)-1;

const char* srcnames_Summer13_V5_SubTotalPileUp[] =
  {"PileUpDataMC","PileUpPtBB", "PileUpPtEC", "PileUpPtHF","SubTotalPileUp"};
const int nsrc_Summer13_V5_SubTotalPileUp = sizeof(srcnames_Summer13_V5_SubTotalPileUp)/sizeof(char*)-1;

const char* srcnames_Summer13_V5_CorrGroups[] =
  {"CorrelationGroupMPFInSitu", "CorrelationGroupIntercalibration",
   "CorrelationGroupFlavor", "CorrelationGroupUncorrelated","Total"};
const int nsrc_Summer13_V5_CorrGroups = sizeof(srcnames_Summer13_V5_CorrGroups)/sizeof(char*)-1;

const char* srcnames_Winter14_V5[] =
  {"AbsoluteStat", "AbsoluteScale", "AbsoluteFlavMap", 
   "AbsoluteMPFBias", "HighPtExtra", "SinglePionECAL", "SinglePionHCAL",
   "FlavorQCD", "Time",
   "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF",
   "RelativePtBB","RelativePtEC1", "RelativePtEC2", "RelativePtHF",
   "RelativeFSR", "RelativeStatEC2", "RelativeStatHF",
   "PileUpDataMC",
   "PileUpPtBB", "PileUpPtEC1", "PileUpPtEC2", "PileUpPtHF",
   "Total"};

const int nsrc_Winter14_V5 = sizeof(srcnames_Winter14_V5)/sizeof(char*)-1;

const char* srcnames_Winter14_V5_SubTotalPt[] =
  {"HighPtExtra", "SinglePionECAL", "SinglePionHCAL","SubTotalPt"}; 
const int nsrc_Winter14_V5_SubTotalPt = sizeof(srcnames_Winter14_V5_SubTotalPt)/sizeof(char*)-1;

const char* srcnames_Winter14_V5_SubTotalRelative[] =
  {"RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF", "RelativePtBB","RelativePtEC1", "RelativePtEC2", "RelativePtHF", "RelativeFSR", "RelativeStatEC2", "RelativeStatHF", "SubTotalRelative"};
const int nsrc_Winter14_V5_SubTotalRelative = sizeof(srcnames_Winter14_V5_SubTotalRelative)/sizeof(char*)-1;

const char* srcnames_Winter14_V5_SubTotalPileUp[] =
  {"PileUpDataMC","PileUpPtBB", "PileUpPtEC1", "PileUpPtEC2", "PileUpPtHF",
   "SubTotalPileUp"};
const int nsrc_Winter14_V5_SubTotalPileUp = sizeof(srcnames_Winter14_V5_SubTotalPileUp)/sizeof(char*)-1;

const char* srcnames_Winter14_V5_CorrGroups[] =
  {"CorrelationGroupMPFInSitu", "CorrelationGroupIntercalibration",
   "CorrelationGroupFlavor", "CorrelationGroupUncorrelated","Total"};
const int nsrc_Winter14_V5_CorrGroups = sizeof(srcnames_Winter14_V5_CorrGroups)/sizeof(char*)-1;


//check all sources listed in src_selection to agree in file1/file2
bool testCommonSources(string file1, string file2,
			 double jetpt=50., double jeteta=2.4, std::string src_selection="Summer13_V2V5CommonSources") {
  if(debug){
    std::cout<< "WARNING: debug/verbose mode activated" << std::endl;
    verbose = true;
  }
  
  const char** srcnames;
  int nsrc;

  // Instantiate uncertainty sources
  if(src_selection == "Summer13_V2V5CommonSources"){
    srcnames = srcnames_Summer13_V2V5CommonSources;
    nsrc=nsrc_Summer13_V2V5CommonSources;
  }

  std::vector<JetCorrectionUncertainty*> vsrcfile1(nsrc);
  std::vector<JetCorrectionUncertainty*> vsrcfile2(nsrc);
  
  // Load individual sources from file1
  if(verbose)cout <<"\nLoading individual sources ("<<nsrc<<") from file1 and file2..." << endl;
  for (int isrc = 0; isrc < nsrc; isrc++) {
    
    const char *name = srcnames[isrc];
    if(verbose)cout << Form("=> %s:%s",file1.c_str(),name) << endl << flush;
    if(verbose)cout << Form("=> %s:%s",file2.c_str(),name) << endl << flush;
    JetCorrectorParameters *p1 = new JetCorrectorParameters(file1.c_str(), name);
    JetCorrectionUncertainty *unc1 = new JetCorrectionUncertainty(*p1);
    JetCorrectorParameters *p2 = new JetCorrectorParameters(file2.c_str(), name);
    JetCorrectionUncertainty *unc2 = new JetCorrectionUncertainty(*p2);
    vsrcfile1[isrc] = unc1;
    vsrcfile2[isrc] = unc2;
  } // for isrc
  
  bool fail = false;
  // Calculate uncertainty per source 
  if(verbose)cout << "Testing each \"common\" source for file1 and file2...";
  for (int isrc = 0; isrc != nsrc; isrc++) {
    
    JetCorrectionUncertainty *unc1 = vsrcfile1[isrc];
    unc1->setJetPt(jetpt);
    unc1->setJetEta(jeteta);
    double sup1 = unc1->getUncertainty(true); // up variation
    JetCorrectionUncertainty *unc2 = vsrcfile2[isrc];
    unc2->setJetPt(jetpt);
    unc2->setJetEta(jeteta);
    double sup2 = unc2->getUncertainty(true); // up variation
    if(debug)printf("%20s: file1 - %10.7f; file2 - %10.7f \n",srcnames[isrc],sup1, sup2); 
    // Check that each "common" source gives same result in file1 and file2
    if (fail = (sup1!=sup2)) {
      cout << "FAIL" << endl;
      cout << "Error: Sources do not agree across files" << endl
	   << "       Difference is "<<sup1<<"(file1) - "<<sup2
	 << "(file2) = " << sup1 - sup2 << endl;
    }
    else
      if(verbose)cout << "ok" << endl;

  } // for isrc

  if(verbose)cout << "Test result: " << (fail ? "FAIL" : "PASS") << endl;
  if(verbose)cout << "Tested eta="<<jeteta<<" pt="<<jetpt<<endl;

  return (!fail);
}

bool testSources(string file1, string file2,
			 double jetpt=50., double jeteta=2.4, std::string src_selection="Summer13_V2") {
  if(debug){
    std::cout<< "WARNING: debug/verbose mode activated" << std::endl;
    verbose = true;
  }
  
  const char** srcnames;
  int nsrc;

  // Instantiate uncertainty sources
  if(src_selection == "Summer13_V2"){
    srcnames = srcnames_Summer13_V2;
    nsrc=nsrc_Summer13_V2;
  }
  else if(src_selection == "Summer13_V5"){
    srcnames = srcnames_Summer13_V5;
    nsrc=nsrc_Summer13_V5;
  }
  else if(src_selection == "Summer13_V5_SubTotalPt"){
    srcnames = srcnames_Summer13_V5_SubTotalPt;
    nsrc=nsrc_Summer13_V5_SubTotalPt;
  }
  else if(src_selection == "Summer13_V5_SubTotalRelative"){
    srcnames = srcnames_Summer13_V5_SubTotalRelative;
    nsrc=nsrc_Summer13_V5_SubTotalRelative;
  }
  else if(src_selection == "Summer13_V5_SubTotalPileUp"){
    srcnames = srcnames_Summer13_V5_SubTotalPileUp;
    nsrc=nsrc_Summer13_V5_SubTotalPileUp;
  }
  else if(src_selection == "Summer13_V5_CorrGroups"){
    srcnames = srcnames_Summer13_V5_CorrGroups;
    nsrc=nsrc_Summer13_V5_CorrGroups;
  }
  else if(src_selection == "Winter14_V5"){
    srcnames = srcnames_Winter14_V5;
    nsrc=nsrc_Winter14_V5;
  }
  else if(src_selection == "Winter14_V5_SubTotalPt"){
    srcnames = srcnames_Winter14_V5_SubTotalPt;
    nsrc=nsrc_Winter14_V5_SubTotalPt;
  }
  else if(src_selection == "Winter14_V5_SubTotalRelative"){
    srcnames = srcnames_Winter14_V5_SubTotalRelative;
    nsrc=nsrc_Winter14_V5_SubTotalRelative;
  }
  else if(src_selection == "Winter14_V5_SubTotalPileUp"){
    srcnames = srcnames_Winter14_V5_SubTotalPileUp;
    nsrc=nsrc_Winter14_V5_SubTotalPileUp;
  }
  else if(src_selection == "Winter14_V5_CorrGroups"){
    srcnames = srcnames_Winter14_V5_CorrGroups;
    nsrc=nsrc_Winter14_V5_CorrGroups;
  }

  std::vector<JetCorrectionUncertainty*> vsrc(nsrc);
  
  // Load individual sources from file1
  if(verbose)cout <<"\nLoading individual sources ("<<nsrc<<") from file1..." << endl;
  for (int isrc = 0; isrc < nsrc; isrc++) {
    
    const char *name = srcnames[isrc];
    if(verbose)cout << Form("=> %s:%s",file1.c_str(),name) << endl << flush;
    JetCorrectorParameters *p = new JetCorrectorParameters(file1.c_str(), name);
    JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
    vsrc[isrc] = unc;
  } // for isrc
  
  // Load total uncertainty source from file1
  if(verbose)cout << "\nLoading Total source from file1..." << endl;
  if(verbose)cout << Form("=> %s:%s",file1.c_str(),srcnames[nsrc]) << endl << flush;
  JetCorrectorParameters *p = new JetCorrectorParameters(file1.c_str(),srcnames[nsrc]);
  JetCorrectionUncertainty *total = new JetCorrectionUncertainty(*p);
  
  JetCorrectorParameters *p2=0;
  JetCorrectionUncertainty *total2=0;
  if(file2!=""){
    // Load total uncertainty from file2
    if(verbose)cout << "\nLoading total uncertainty from file2..." << endl;
    if(verbose)cout << "=> " << file2 << endl << flush;
    p2 = new JetCorrectorParameters(file2.c_str());
    total2 = new JetCorrectionUncertainty(*p2);
  }
  
  // Calculate uncertainty per source and as a total
  double sum2_up(0), sum2_dw(0), sum2(0);
  for (int isrc = 0; isrc != nsrc; isrc++) {
    
    JetCorrectionUncertainty *unc = vsrc[isrc];
    unc->setJetPt(jetpt);
    unc->setJetEta(jeteta);
    double sup = unc->getUncertainty(true); // up variation
    if(debug)printf("%20s: %10.7f \n",srcnames[isrc],sup); 
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
  
  double uncert2=-1;
  if(file2!=""){
    total2->setJetPt(jetpt);
    total2->setJetEta(jeteta);
    uncert2 = total2->getUncertainty(true);
  }
  
  bool fail = false;

  // Check that quadratic sum of sources equals total uncertainty
  if(verbose)cout << "\nTesting sum versus file1 total...";
  if (fail = (fabs(uncert - sqrt(sum2)) > 5e-4)) {
    cout << "FAIL" << endl;
    cout << "Error: Quadratic sum of sources does not agree with total" << endl
	 << "       Difference is "<<uncert<<"(total) - "<<sqrt(sum2)
	 << "(sum) = " << uncert - sqrt(sum2) << endl;
    //assert(fabs(uncert - sqrt(sum2_up)) < 1e-3);
  }
  else
    if(verbose)cout << "ok" << endl;
  
  if(file2!=""){
    if(verbose)cout << "Testing sum versus file2 total...";
    if (fail = (fabs(uncert2 - sqrt(sum2)) > 5e-4)) {
      cout << "FAIL" << endl;
      cout << "Error: Quadratic sum of sources does not agree with total" << endl
	   << "       Difference is "<<uncert2<<"(total) - "<<sqrt(sum2)
	   << "(sum) = " << uncert2 - sqrt(sum2) << endl;
      //assert(fabs(uncert2 - sqrt(sum2_up)) < 1e-3);
    }
    else
      if(verbose)cout << "ok" << endl;
  }

  if(verbose)cout << "Test result: " << (fail ? "FAIL" : "PASS") << endl;
  if(verbose)cout << "Tested eta="<<jeteta<<" pt="<<jetpt<<endl;

  if(fail)nFailedTests++;

  return (!fail);
} // testSources



//helper function to run cross-checks for a number of important eta/pt-combinations
bool testSourcesLoopVars(string file1, string file2, std::string src_selection="Summer13_V2"){
  double etas[] = {0, 2.0, 2.7, 4.0, 5.2};
  const int neta = sizeof(etas)/sizeof(etas[0]);
  
  double pts[] = {30, 100, 500, 1000};
  const int npt = sizeof(pts)/sizeof(pts[0]);

  bool tempverbose= verbose;
  verbose=false;
  for (int ipt = 0; ipt != npt; ++ipt) {
    for (int ieta = 0; ieta != neta; ++ieta) {
      testSources(file1, file2,	pts[ipt], etas[ieta], src_selection);
    }
  }
  verbose = tempverbose;

  std::cout << "Tested " << file1.c_str() << (file2!="" ? " and " : "") << file2.c_str() << " for " << src_selection << std::endl;
  std::cout << "Test result: " << (nFailedTests ? "FAIL" : "PASS") << endl;
  return nFailedTests==0 ? true : false;
}

