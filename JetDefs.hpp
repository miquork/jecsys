// Adapted from D0 Experiment jetcorr/jetcorr/JetDefs.hpp

#ifndef INC_JETDEFS
#define INC_JETDEFS
///////////////////////////////////////////////////////////////////////////////
// $Id: JetDefs.hpp,v 1.2 2008/12/19 17:16:05 voutila Exp $
// 
// File: JetDefs.hpp
//
// Purpose:  various definitions and parameters for jet corrections package
//
// Created:  Nov-3-20017   Mikko Voutilainen (adapted from D0 Experiment)
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "ErrorTypes.hpp"

namespace jec {

  const ErrorTypes kNone                = ErrorTypes(0L, 0L);
  // Pile-up systematics, bits 1-5
  const ErrorTypes kPileUpDataMC        = ErrorTypes(0L, 1L << 0);
  const ErrorTypes kPileUpBias          = ErrorTypes(0L, 1L << 1); //zero
  const ErrorTypes kPileUpPtBB          = ErrorTypes(0L, 1L << 2);
  const ErrorTypes kPileUpPtEC          = ErrorTypes(0L, 1L << 3);
  const ErrorTypes kPileUpPtHF          = ErrorTypes(0L, 1L << 4);
  // Relative correction systematics, bits 6-16
  const ErrorTypes kRelativeJEREC1      = ErrorTypes(0L, 1L << 6);
  const ErrorTypes kRelativeJEREC2      = ErrorTypes(0L, 1L << 7);
  const ErrorTypes kRelativeJERHF       = ErrorTypes(0L, 1L << 8);
  const ErrorTypes kRelativeFSR         = ErrorTypes(0L, 1L << 9);
  const ErrorTypes kRelativeStatEC2     = ErrorTypes(0L, 1L << 10);
  const ErrorTypes kRelativeStatHF      = ErrorTypes(0L, 1L << 11);
  //const ErrorTypes kRelativeSample      = ErrorTypes(0L, 1L << 12);
  const ErrorTypes kRelativePtBB        = ErrorTypes(0L, 1L << 13);
  const ErrorTypes kRelativePtEC1       = ErrorTypes(0L, 1L << 14);
  const ErrorTypes kRelativePtEC2       = ErrorTypes(0L, 1L << 15);
  const ErrorTypes kRelativePtHF        = ErrorTypes(0L, 1L << 16);
  // Absolute scale (pT dependence) systematics, bits 17-27
  const ErrorTypes kAbsoluteScale       = ErrorTypes(0L, 1L << 17);
  const ErrorTypes kAbsoluteSPRE        = ErrorTypes(0L, 1L << 18);
  const ErrorTypes kAbsoluteSPRH        = ErrorTypes(0L, 1L << 19);
  const ErrorTypes kAbsoluteECAL        = ErrorTypes(0L, 1L << 20); //off
  const ErrorTypes kAbsoluteTrack       = ErrorTypes(0L, 1L << 21); //off
  const ErrorTypes kAbsoluteFrag        = ErrorTypes(0L, 1L << 22);
  // Flavor systematics for L5(residual), bits 26-35
  // => Obsoleted
  // optional Flavors systematics for pure flavors or mixtures, bits 36-45
  const ErrorTypes kFlavorQCD           = ErrorTypes(0L, 1L << 36); //default
  const ErrorTypes kFlavorZJet          = ErrorTypes(0L, 1L << 37); //opt
  const ErrorTypes kFlavorPhotonJet     = ErrorTypes(0L, 1L << 38); //opt
  const ErrorTypes kFlavorPureGluon     = ErrorTypes(0L, 1L << 39); //opt
  const ErrorTypes kFlavorPureQuark     = ErrorTypes(0L, 1L << 40); //opt
  const ErrorTypes kFlavorPureCharm     = ErrorTypes(0L, 1L << 41); //opt
  const ErrorTypes kFlavorPureBottom    = ErrorTypes(0L, 1L << 42); //opt
  // time dependence, bits 46-50
  const ErrorTypes kTime                = ErrorTypes(0L, 1L << 46);

  // Add this to single sources (e.g. kPileUpDataMC) to get unsigned uncertainty
  const ErrorTypes kDoUnsigned          = ErrorTypes(0L, 1L << 51);

  // Combinations of bits
  const ErrorTypes kPileUpPt            = kPileUpPtBB | kPileUpPtEC | kPileUpPtHF;
  const ErrorTypes kRelativeJER         = kRelativeJEREC1 | kRelativeJEREC2 | kRelativeJERHF;
  const ErrorTypes kRelativePt          = kRelativePtBB | kRelativePtEC1 | kRelativePtEC2 | kRelativePtHF;
  const ErrorTypes kRelativeStat        = kRelativeStatEC2 | kRelativeStatHF;
  const ErrorTypes kAbsoluteSPR         = kAbsoluteSPRE | kAbsoluteSPRH;

  // SubTotalPileUp, SubTotalRelative, [SubTotalAbsolute], SubTotalPt
  const ErrorTypes kPileUp              = kPileUpDataMC | kPileUpPt;
  const ErrorTypes kRelative            = kRelativeJER | kRelativeFSR | kRelativeStat | kRelativePt;
  const ErrorTypes kAbsolute            = kAbsoluteScale | kAbsoluteFrag | kAbsoluteSPR;
  const ErrorTypes kAbsolutePt          = kAbsoluteFrag | kAbsoluteSPR;

  // Test mask: only one of these should be on at a time
  const ErrorTypes kFlavorMask          = kFlavorQCD | kFlavorZJet | kFlavorPhotonJet | kFlavorPureQuark | kFlavorPureGluon | kFlavorPureCharm | kFlavorPureBottom;

  // Total uncertainty bits
  const ErrorTypes kMC = kPileUpDataMC | kRelative | kAbsolute | kFlavorQCD | kTime; // for Data/MC comparisons (excludes kPileUpPt)
  const ErrorTypes kData = kPileUp | kRelative | kAbsolute | kFlavorQCD | kTime; // for analyses with only data

  // SubTotalNoFlavor
  const ErrorTypes kDataNoFlavor = kData & ~kFlavorQCD;
} // namespace jec

#endif // INC_JETDEFS
