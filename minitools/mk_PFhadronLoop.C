// Purpose: Run studies on isolated PF charged hadrons for calibration
//          Histogrammin/profiling and plotting are separated
// Run with 'root -l -b -q minitools/mk_PFhadronLoop.C

R__LOAD_LIBRARY(minitools/PFhadronLoop.C+g)


void mk_PFhadronLoop() {

  TChain *c = new TChain("ak4/ProcessedTree");
  //c->AddFile("rootfiles/Hadrons_UL16/mc/SingleNeutrino_UL16APV.root");
  //c->AddFile("rootfiles/Hadrons_UL16/mc/SingleNeutrino_UL16nonAPV.root");

  // Hadrons_18MC.root by hand
  //c->AddFile("rootfiles/Hadrons/Had18MC.root");

  // Hadrons_18DT.root by hand
  c->AddFile("rootfiles/Hadrons/Had18A.root");
  c->AddFile("rootfiles/Hadrons/Had18B.root");
  c->AddFile("rootfiles/Hadrons/Had18C.root");
  c->AddFile("rootfiles/Hadrons/Had18D.root");

  PFhadronLoop p(c);
  p.Loop();
} // mk_PFhadronLoop
