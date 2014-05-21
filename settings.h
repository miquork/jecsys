#ifndef __settings_h__
#define __settings_h__
#include <string>

//const double _lumi = 4961.;//4942.9; // pixelLumi
const double _lumi = 19700;//18893.;//19466.;//9213.;//4961.;//5003.;//2930.;//1616.;//570.35;//366.221;//355.3;//21.3; // 2012 lumiCalc2
//const double _lumi = 4665;//4832;//4958;//4666;//3172;//2177.;//1100.;//715.;//492;//36;

// Algorithm to use ("AK7" or "AK5")
std::string _algo = "AK7";
// do additional AK5/AK7 histograms
const bool _ak5ak7 = (_algo=="AK7" && true);
// also run CaloJets (needed for trigger emulation, TP+trigger)
const bool dofriends = true;
// Apply "CHS" through betaStar
const bool _doCHS = false;
// Veto jets near ECAL boundaries in JetID
const bool _doECALveto = false;//true;
// Fill histograms separately for each five eras
const bool _doEras = false;//true;
const bool _useIOV = false;//true; // use IOV's for data
// Process pThatbins instead of flat sample
const bool _pthatbins = false;
// Correct for trigger efficiency based on MC
const bool _dotrigeff = false;//true;
// Correct for time-dependence (prescales) in data
const bool _dotimedep = false;//true;
// For creating smearing matrix
const bool _doMatrix = false;

// which MC version
//const bool _spring10 = true;//false;
//const bool _spring10 = false;
//const bool _fall10 = true; //mc-fall10-v2,data-36pbinv-nov4-cmssw387p2-fall10jec
// patch MinBias in normalizeHistos.C
//const bool _fixmb = false;
// reapply json selection based on the latest one (check lumicalc if false!)
const bool _dojson = true;//false;//true;
// check that run / lumi section was listed in the .csv file
const bool _dolumcheck = false;
// veto bad run list
const bool _dorunveto = false;
// check for duplicates (warning: takes a lot of memory!)
const bool _checkduplicates = true;//false;
// only load selected branches
const bool _quick = true;
// print out debugging information
const bool _debug = false;
// use parameterized tag rates instead of binned
//const bool _usetagfit = true;
//const bool _usetagfit0 = false;//true;
// Center uncertainties around ansatz (true) or data (false)
const bool _centerOnAnsatz = false;
const bool _centerOnTheory = true;
//const bool _centerBOnAnsatz = true;
// Plot Pythia for final PRL results
const bool _plotPythia = false;
// Minimum and maximum pT range to be plotted and fitted
const double _recopt = 21.;//15.;//25;//15;//12;//10.;
const double fitptmin = 57;//15;//12;//10.;
//const double fitbptmin = 57;//15;//12;//10.;
const double xminpas = 74;//24.;
const double xmin = 20.;//10.;
const double xmax = 2500.;//1500;//900;//500;
// Minimum and maximum pT range for b-tagging to be fitted
//const double bxmin = 30.;
//const double bxmax = 1000.;//500.;//400;//220.;//160;

// Write histograms to PAS file
const bool _pas = true;
// Draw b-jets against MC@NLO instead of reco MC
//const bool _mcnlo = true;
// Draw againts HERAPDF1.5 instead of PDF4LHC
const bool _herapdf = false;

// Produce extra file types
const bool _eps = false;//true;
const bool _pdf = true;
const bool _gif = false;
const bool _png = false;

#endif // __settings_h__
