// Purpose: Run studies on isolated PF charged hadrons for calibration
//          Histogrammin/profiling and plotting are separated
// Run with 'root -l -b -q minitools/mk_PFhadronLoop.C

R__LOAD_LIBRARY(minitools/PFhadronLoop.C+g)
const string mode = "DT1718";
//const string mode = "MC1718";
//const string mode = "DT18";
//const string mode = "MC18";
//const string mode = "DT17";
//const string mode = "MC17";
//const string mode = "DT16APV";
//const string mode = "MC16APV";
//const string mode = "DT16GH";
//const string mode = "MC16GH";

//const string mode = "DT18A";
//const string mode = "MC18A";
//const string mode = "DT18B";
//const string mode = "MC18B";
//const string mode = "DT18C";
//const string mode = "MC18C";
//const string mode = "DT18D";
//const string mode = "MC18D";

void mk_PFhadronLoop() {

  TChain *c = new TChain("ak4/ProcessedTree");
  if (mode=="MC16APV") {
    c->AddFile("rootfiles/Hadrons_UL16/mc/SingleNeutrino_UL16APV.root");
  }
  if (mode=="MC16GH") {
    c->AddFile("rootfiles/Hadrons_UL16/mc/SingleNeutrino_UL16nonAPV.root");
  }

  // Hadrons_161718MC.root by hand
  //c->AddFile("rootfiles/Hadrons_UL16/mc/SingleNeutrino_UL16APV.root");
  //c->AddFile("rootfiles/Hadrons_UL16/mc/SingleNeutrino_UL16nonAPV.root");
  if (mode=="MC1718" || mode=="MC17") {
    c->AddFile("rootfiles/Hadrons/Had17MC.root");
  }
  if (mode=="MC1718" || mode=="MC18" ||
      mode=="MC18A" || mode=="MC18B" || mode=="MC18C" || mode=="MC18D") { 
    c->AddFile("rootfiles/Hadrons/Had18MC.root");
    c->AddFile("rootfiles/Hadrons/Had18MC_HEM.root");
  }

  // Hadrons_161718.root, Hadrons_16GH1718.root by hand
  if (mode=="DT16APV") {
    c->AddFile("rootfiles/Hadrons_UL16/data/ZeroBias_UL16RunB.root");
    c->AddFile("rootfiles/Hadrons_UL16/data/ZeroBias_UL16RunC.root");
    c->AddFile("rootfiles/Hadrons_UL16/data/ZeroBias_UL16RunD.root");
    c->AddFile("rootfiles/Hadrons_UL16/data/ZeroBias_UL16RunE.root");
    c->AddFile("rootfiles/Hadrons_UL16/data/ZeroBias_UL16RunFe.root");
  }
  if (mode=="DT16GH") {
    c->AddFile("rootfiles/Hadrons_UL16/data/ZeroBias_UL16RunFl.root");
    c->AddFile("rootfiles/Hadrons_UL16/data/ZeroBias_UL16RunG.root");
    c->AddFile("rootfiles/Hadrons_UL16/data/ZeroBias_UL16RunH.root");
  }

  // Hadrons_1718.root by hand
  if (mode=="DT1718" || mode=="DT17") {
    c->AddFile("rootfiles/Hadrons/Had17B.root");
    c->AddFile("rootfiles/Hadrons/Had17C.root");
    c->AddFile("rootfiles/Hadrons/Had17D.root");
    c->AddFile("rootfiles/Hadrons/Had17E.root");
    c->AddFile("rootfiles/Hadrons/Had17F.root");
  }

  // Hadrons_18DT.root by hand
  if (mode=="DT1718" || mode=="DT18") {
    c->AddFile("rootfiles/Hadrons/Had18A.root");
    c->AddFile("rootfiles/Hadrons/Had18B.root");
    c->AddFile("rootfiles/Hadrons/Had18C.root");
    c->AddFile("rootfiles/Hadrons/Had18D.root");
  }
  if (mode=="DT18A") c->AddFile("rootfiles/Hadrons/Had18A.root");
  if (mode=="DT18B") c->AddFile("rootfiles/Hadrons/Had18B.root");
  if (mode=="DT18C") c->AddFile("rootfiles/Hadrons/Had18C.root");
  if (mode=="DT18D") c->AddFile("rootfiles/Hadrons/Had18D.root");

  PFhadronLoop p(c,mode);
  p.Loop();
} // mk_PFhadronLoop
