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

  const ErrorTypes kNone                    = ErrorTypes(0L, 0L);
  // Pile-up systematics, bits 1-5        
  const ErrorTypes kPileUpDataMC            = ErrorTypes(0L, 1L << 0);
  //const ErrorTypes kPileUpBias              = ErrorTypes(0L, 1L << 1); //zero
  const ErrorTypes kPileUpPtRef             = ErrorTypes(0L, 1L << 1);
  const ErrorTypes kPileUpPtBB              = ErrorTypes(0L, 1L << 2);
  const ErrorTypes kPileUpPtEC1             = ErrorTypes(0L, 1L << 3);
  const ErrorTypes kPileUpPtEC2             = ErrorTypes(0L, 1L << 4);
  const ErrorTypes kPileUpPtHF              = ErrorTypes(0L, 1L << 5);
  // Relative correction systematics, bits 6-16
  const ErrorTypes kRelativeJEREC1          = ErrorTypes(0L, 1L << 6);
  const ErrorTypes kRelativeJEREC2          = ErrorTypes(0L, 1L << 7);
  const ErrorTypes kRelativeJERHF           = ErrorTypes(0L, 1L << 8);
  const ErrorTypes kRelativeFSR             = ErrorTypes(0L, 1L << 9);
  const ErrorTypes kRelativeStatEC          = ErrorTypes(0L, 1L << 10);
  const ErrorTypes kRelativeStatHF          = ErrorTypes(0L, 1L << 11);
  const ErrorTypes kRelativeStatFSR         = ErrorTypes(0L, 1L << 12);
  //const ErrorTypes kRelativeSample          = ErrorTypes(0L, 1L << 12);
  const ErrorTypes kRelativePtBB            = ErrorTypes(0L, 1L << 13);
  const ErrorTypes kRelativePtEC1           = ErrorTypes(0L, 1L << 14);
  const ErrorTypes kRelativePtEC2           = ErrorTypes(0L, 1L << 15);
  const ErrorTypes kRelativePtHF            = ErrorTypes(0L, 1L << 16);
  const ErrorTypes kRelativeBal             = ErrorTypes(0L, 1L << 26); //Sum16 
  // Absolute scale (pT dependence) systematics, bits 17-27
  const ErrorTypes kAbsoluteScale           = ErrorTypes(0L, 1L << 17);
  const ErrorTypes kAbsoluteSPRE            = ErrorTypes(0L, 1L << 18);
  const ErrorTypes kAbsoluteSPRH            = ErrorTypes(0L, 1L << 19);
  const ErrorTypes kAbsoluteECAL            = ErrorTypes(0L, 1L << 20); //off
  const ErrorTypes kAbsoluteTrack           = ErrorTypes(0L, 1L << 21); //off
  const ErrorTypes kAbsoluteFrag            = ErrorTypes(0L, 1L << 22);
  const ErrorTypes kAbsoluteStat            = ErrorTypes(0L, 1L << 23); //splitting for correlation groups
  const ErrorTypes kAbsoluteMPFBias         = ErrorTypes(0L, 1L << 24);	//splitting for correlation groups
  const ErrorTypes kAbsoluteFlavorMapping   = ErrorTypes(0L, 1L << 25); //splitting for correlation groups

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
  // time dependence, bits 46-51
  //const ErrorTypes kTimeEta             = ErrorTypes(0L, 1L << 46);
  const ErrorTypes kTimePtEta           = ErrorTypes(0L, 1L << 47);
  // optional time bits for individual epochs (not included in total), 48-51
  const ErrorTypes kTimeRunBCD          = ErrorTypes(0L, 1L << 48); //opt
  //const ErrorTypes kTimeRunE            = ErrorTypes(0L, 1L << 49); //opt
  //const ErrorTypes kTimeRunF            = ErrorTypes(0L, 1L << 50); //opt
  const ErrorTypes kTimeRunEF           = ErrorTypes(0L, 1L << 49); //opt
  //const ErrorTypes kTimeRunGH           = ErrorTypes(0L, 1L << 51); //opt
  const ErrorTypes kTimeRunG            = ErrorTypes(0L, 1L << 50); //opt
  const ErrorTypes kTimeRunH            = ErrorTypes(0L, 1L << 51); //opt
  // optional PU term for <mu>=0 sample (bias from fitting L2Res with <mu>=20)
  const ErrorTypes kPileUpMuZero        = ErrorTypes(0L, 1L << 52); //opt
  const ErrorTypes kPileUpEnvelope      = ErrorTypes(0L, 1L << 53); //xtra

  const ErrorTypes kRunI = ErrorTypes(0L, 1L << 54); // reference

  // Add this to single sources (e.g. kPileUpDataMC) to get unsigned uncertainty
  //const ErrorTypes kDoUnsigned          = ErrorTypes(0L, 1L << 54);

  //  // Extra bits to to be able to write-out extra entries for correlation groups, bits 55-58 
  //  // (even though they are only composed of a single source)
  //  //TOPLHCWG CMS/ATLAS JEC correlation groups
  const ErrorTypes kCorrelationGroupMPFInSitu            = ErrorTypes(1L << 55, 0L) | kAbsoluteMPFBias;
  const ErrorTypes kCorrelationGroupFlavor               = ErrorTypes(1L << 56, 0L) | kAbsoluteFlavorMapping | kFlavorQCD ;
  const ErrorTypes kCorrelationGroupIntercalibration     = ErrorTypes(1L << 57, 0L) | kRelativeFSR;
  const ErrorTypes kCorrelationGroupbJES                 = ErrorTypes(1L << 58, 0L) | kFlavorQCD;


  // Combinations of bits
  const ErrorTypes kTime                = kTimePtEta;//kTimeEta | kTimePt;
  const ErrorTypes kPileUpPtEta         = kPileUpPtBB | kPileUpPtEC1 | kPileUpPtEC2 | kPileUpPtHF;
  const ErrorTypes kPileUpPt            = kPileUpPtRef | kPileUpPtEta;
  const ErrorTypes kRelativeJER         = kRelativeJEREC1 | kRelativeJEREC2 | kRelativeJERHF;
  const ErrorTypes kRelativePt          = kRelativePtBB | kRelativePtEC1 | kRelativePtEC2 | kRelativePtHF;
  const ErrorTypes kRelativeStat        = kRelativeStatEC | kRelativeStatHF | kRelativeStatFSR;
  const ErrorTypes kAbsoluteSPR         = kAbsoluteSPRE | kAbsoluteSPRH;

  // SubTotalPileUp, SubTotalRelative, [SubTotalAbsolute], SubTotalPt
  const ErrorTypes kPileUp              = kPileUpDataMC | kPileUpPt;
  const ErrorTypes kRelative            = kRelativeJER | kRelativeFSR | kRelativeStat | kRelativePt | kRelativeBal;
  const ErrorTypes kAbsolutePt          = kAbsoluteFrag | kAbsoluteSPR;
  const ErrorTypes kAbsoluteFlat        = kAbsoluteStat | kAbsoluteMPFBias | kAbsoluteFlavorMapping | kAbsoluteScale;
  const ErrorTypes kAbsolute            = kAbsoluteFlat | kAbsolutePt;

  // Test mask: only one of these should be on at a time
  const ErrorTypes kFlavorMask          = kFlavorQCD | kFlavorZJet | kFlavorPhotonJet | kFlavorPureQuark | kFlavorPureGluon | kFlavorPureCharm | kFlavorPureBottom;
  //const ErrorTypes kTimePtMask          = kTimePt | kTimePtRunBCD | kTimePtRunE | kTimePtRunF | kTimePtRunGH;
  //const ErrorTypes kTimePtEtaMask          = kTimePtEta | kTimeRunBCD | kTimeRunE | kTimeRunF | kTimeRunGH;
  const ErrorTypes kTimePtEtaMask          = kTimePtEta | kTimeRunBCD | kTimeRunEF | kTimeRunG | kTimeRunH;

  // Total uncertainty bits
  const ErrorTypes kMC = kPileUpDataMC | kRelative | kAbsolute | kFlavorQCD | kTime; // for Data/MC comparisons (excludes kPileUpPt)
  const ErrorTypes kData = kPileUp | kRelative | kAbsolute | kFlavorQCD | kTime; // for analyses with only data

  // SubTotalNoFlavor => for mixing in flavor separately
  const ErrorTypes kDataNoFlavor = kData & ~kFlavorQCD;
  // SubTotalNoFlavor => for QCD / inclusive jets
  const ErrorTypes kDataNoTime = kData & ~kTime;
  // SubTotalNoFlavorNoTime => for top mass
  const ErrorTypes kDataNoFlavorNoTime = kData & ~kFlavorQCD & ~kTime;

  //TOPLHC CMS/ATLAS JEC correlation groups
  const ErrorTypes kCorrelationGroupPartiallyCorrelated          = kFlavorQCD | kRelativeFSR | kAbsoluteFlavorMapping | kAbsoluteMPFBias;
  const ErrorTypes kCorrelationGroupUncorrelated                 = ~kCorrelationGroupPartiallyCorrelated & kData;


} // namespace jec

#endif // INC_JETDEFS
