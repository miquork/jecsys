// Purpose: Run studies on isolated PF charged hadrons for calibration
//          Histogrammin/profiling and plotting are separated
// Run with 'root -l -b -q minitools/mk_PFhadronLoop.C

R__LOAD_LIBRARY(minitools/PFhadronLoop.C+g)


void mk_PFhadronLoop() {

  TChain *c = new TChain("ak4/ProcessedTree");
  //c->AddFile("rootfiles/Hadrons_UL16/mc/SingleNeutrino_UL16APV.root");
  //c->AddFile("rootfiles/Hadrons_UL16/mc/SingleNeutrino_UL16nonAPV.root");

  /*
  // Hadrons_161718MC.root by hand
  //c->AddFile("rootfiles/Hadrons_UL16/mc/SingleNeutrino_UL16APV.root");
  c->AddFile("rootfiles/Hadrons_UL16/mc/SingleNeutrino_UL16nonAPV.root");
  c->AddFile("rootfiles/Hadrons/Had17MC.root");
  c->AddFile("rootfiles/Hadrons/Had18MC.root");
  c->AddFile("rootfiles/Hadrons/Had18MC_HEM.root");
  */

  // Hadrons_161718.root, Hadrons_16GH1718.root by hand
  //c->AddFile("rootfiles/Hadrons_UL16/data/ZeroBias_UL16RunB.root");
  //c->AddFile("rootfiles/Hadrons_UL16/data/ZeroBias_UL16RunC.root");
  //c->AddFile("rootfiles/Hadrons_UL16/data/ZeroBias_UL16RunD.root");
  //c->AddFile("rootfiles/Hadrons_UL16/data/ZeroBias_UL16RunE.root");
  //c->AddFile("rootfiles/Hadrons_UL16/data/ZeroBias_UL16RunFe.root");
  c->AddFile("rootfiles/Hadrons_UL16/data/ZeroBias_UL16RunFl.root");
  c->AddFile("rootfiles/Hadrons_UL16/data/ZeroBias_UL16RunG.root");
  c->AddFile("rootfiles/Hadrons_UL16/data/ZeroBias_UL16RunH.root");

  // Hadrons_1718.root by hand
  c->AddFile("rootfiles/Hadrons/Had17B.root");
  c->AddFile("rootfiles/Hadrons/Had17C.root");
  c->AddFile("rootfiles/Hadrons/Had17D.root");
  c->AddFile("rootfiles/Hadrons/Had17E.root");
  c->AddFile("rootfiles/Hadrons/Had17F.root");

  // Hadrons_18DT.root by hand
  c->AddFile("rootfiles/Hadrons/Had18A.root");
  c->AddFile("rootfiles/Hadrons/Had18B.root");
  c->AddFile("rootfiles/Hadrons/Had18C.root");
  c->AddFile("rootfiles/Hadrons/Had18D.root");


  PFhadronLoop p(c);
  p.Loop();
} // mk_PFhadronLoop
