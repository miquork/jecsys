# jecsys
CMS jet energy scale systematics

This 'jecsys' package is meant for deriving and testing CMS jet energy scale uncertainties.
If you are a USER (rather than EXPERT) of JEC uncertainties, please start with the 'test' subdirectory.

There you will find:
- mk_drawUncertainty.C (using drawUncertainty.C) to draw an uncertainty source of your choice
- mk_drawSourceCorrelations.C to draw JEC correlations across pT and eta
- mk_drawCMSResponse.C to draw the JEC as a function of eta for various jet pT choices

The UncertaintySources text files, as well as the central values of JEC are available at
  https://github.com/cms-jet/JECDatabase/tree/master/tarballs
  
If you are interested in the input data used for JEC global fit and uncertainties, you can look into
  rootfiles/jecdata[X].root
were [X] would vary depending on dataset (e.g. BCD, EF, G and H for 2016 data). However, these files are
only reliable when used together with a certain tag to identify the JEC used, e.g. 'Summer16_23SepV4'.

For the real EXPERTS (i.e. if you get EPR from JetMET work):
- mk_reprocess.C (reprocess.C, softrad.C, globalFitL3Res.C) to produce jecdata.root and perform L3Res global fit
- mk_drawJetCorrectionUncertainty.C (JECUncertainty.*, ErrorTypes.*, JetDefs.hpp) to derive UncertaintySources
- mk_testSources.C to test UncertaintySources
