// Purpose: plot time stability of physics observables from data
// Current short list: ptlep, ptt1; mw, rbq, mt, mlb; mu(truePU), NPV, rho
// Uses output from minitools/mk_hadW.C
// run with 'root -l -b -q minitools/drawTimeStability.C+g'
#include "TFile.h"
#include "TProfile.h"
#include "TF1.h"
#include "TGraphErrors.h"

#include "../tdrstyle_mod15.C"

void drawTimeStability() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;
//const char *cm = mode.c_str();

  TFile *fd = new TFile(Form("rootfiles/hadWUL%s_MPDGcorrNoW.root",cm),"READ");
  assert(fd && !fd->IsZombie());

  TFile *fm = new TFile(Form("rootfiles/hadWMC%s_MPDGcorrNoW.root",cm),"READ");
  assert(fm && !fm->IsZombie());

  curdir->cd();
} // drawTimeStability
