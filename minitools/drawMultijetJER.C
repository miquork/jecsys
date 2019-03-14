// Purpose: draw plots demonstratig JER and unclustered energy bias in multijets
//          using dijet and inclusive jet events
// run with  'root -l minitools/drawMultijetJER.C'
#include "../tdrstyle_mod15.C"


void drawMultijetJER(string mode = "EC1") {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  //TFile *fdt = new TFile("rootfiles/output-DATA-2b-legacy16-GH.root","READ");
  //TFile *fdt = new TFile("rootfiles/output-DATA-2b-17nov17-DE.root","READ");
  //TFile *fdt = new TFile("rootfiles/output-DATA-2b-Fall18-RunD.root","READ");
  //TFile *fdt = new TFile("rootfiles/output-DATA-2b-Fall18V8-A.root","READ");
  TFile *fdt = new TFile("rootfiles/output-DATA-2b-Fall18V8-B.root","READ");
  //TFile *fdt = new TFile("rootfiles/output-DATA-2b-Fall18V8-C.root","READ");
  //TFile *fdt = new TFile("rootfiles/output-DATA-2b-Fall18V8-D.root","READ");
  assert(fdt && !fdt->IsZombie());

  //TFile *fmc = new TFile("rootfiles/output-MC-2b-legacy16-GH.root","READ");
  //TFile *fmc = new TFile("rootfiles/output-P8CP5-2b-17nov17-DE.root","READ");
  //TFile *fmc = new TFile("rootfiles/output-MC-2b-Fall18-RunD.root","READ");
  //TFile *fmc = new TFile("rootfiles/output-MCNU-2b-Fall18V8-A.root","READ");
  TFile *fmc = new TFile("rootfiles/output-MCNU-2b-Fall18V8-B.root","READ");
  //TFile *fmc = new TFile("rootfiles/output-MCNU-2b-Fall18V8-C.root","READ");
  //TFile *fmc = new TFile("rootfiles/output-MCNU-2b-Fall18V8-D.root","READ");
  assert(fmc && !fmc->IsZombie());

  //TFile *fhw = new TFile("rootfiles/output-HW-2b-legacy16-GH.root","READ");
  //TFile *fhw = new TFile("rootfiles/output-HW-2b-17nov17-DE.root","READ");
  TFile *fhw = 0;
  //assert(fhw && !fhw->IsZombie());

  TFile *fm1 = 0;
  //TFile *fm1 = new TFile("rootfiles/output-P8M1-2b-17nov17-DE.root","READ");
  //TFile *fm1 = 0;
  //assert(fm1 && !fm1->IsZombie());

  string sp = "";
  if (mode=="BB") sp = "Standard/Eta_0.0-1.3";
  if (mode=="EC0") sp = "Standard/Eta_1.5-2.0";
  if (mode=="EC1") sp = "Standard/Eta_2.0-2.5";
  if (mode=="EC2") sp = "Standard/Eta_2.5-3.0";
  if (mode=="HF") sp = "Standard/Eta_3.2-4.7";
  const char *cp = sp.c_str();
  const char *cm = mode.c_str();

  //const char *cp = "Standard/Eta_0.0-1.3"; //BB
  //const char *cp = "Standard/Eta_2.0-2.5"; // EC1
  //const char *cp = "Standard/Eta_2.5-3.0"; // EC2

  // jetphys/HistosFill.C:
  //   double metstuff = met1 * cos(DPhi(metphi1, phiprobe));
  //   h->pmpfx->Fill(pttag, 1 + metstuff / pttag, _w);
  //   h->pmpfy->Fill(ptprobe, 1 + metstuff / ptprobe, _w);
  //   h->pmpfz->Fill(ptave, 1 + metstuff / ptave, _w);

  TProfile *pdttag = (TProfile*)fdt->Get(Form("%s/pmpfx",cp)); assert(pdttag);
  TProfile *pdtpro = (TProfile*)fdt->Get(Form("%s/pmpfy",cp)); assert(pdtpro);
  TProfile *pdtave = (TProfile*)fdt->Get(Form("%s/pmpfz",cp)); assert(pdtave);

  TProfile *pmctag = (TProfile*)fmc->Get(Form("%s/pmpfx",cp)); assert(pmctag);
  TProfile *pmcpro = (TProfile*)fmc->Get(Form("%s/pmpfy",cp)); assert(pmcpro);
  TProfile *pmcave = (TProfile*)fmc->Get(Form("%s/pmpfz",cp)); assert(pmcave);

  TProfile *pdt0 = (TProfile*)fdt->Get(Form("%s/pmpf",cp)); assert(pdt0);
  TProfile *pdt1 = (TProfile*)fdt->Get(Form("%s/pmpf1",cp)); assert(pdt1);
  TProfile *pdt2 = (TProfile*)fdt->Get(Form("%s/pmpf2",cp)); assert(pdt2);

  TProfile *pmc0 = (TProfile*)fmc->Get(Form("%s/pmpf",cp)); assert(pmc0);
  TProfile *pmc1 = (TProfile*)fmc->Get(Form("%s/pmpf1",cp)); assert(pmc1);
  TProfile *pmc2 = (TProfile*)fmc->Get(Form("%s/pmpf2",cp)); assert(pmc2);

  TProfile *phwtag(0), *phwpro(0), *phwave(0), *phw0(0), *phw1(0), *phw2(0);
  if (fhw) {
    phwtag = (TProfile*)fhw->Get(Form("%s/pmpfx",cp)); //assert(phwtag);
    phwpro = (TProfile*)fhw->Get(Form("%s/pmpfy",cp)); //assert(phwpro);
    phwave = (TProfile*)fhw->Get(Form("%s/pmpfz",cp)); //assert(phwave);

    phw0 = (TProfile*)fhw->Get(Form("%s/pmpf",cp)); //assert(phw0);
    phw1 = (TProfile*)fhw->Get(Form("%s/pmpf1",cp)); //assert(phw1);
    phw2 = (TProfile*)fhw->Get(Form("%s/pmpf2",cp)); //assert(phw2);
  }
  TProfile *pm1tag(0), *pm1pro(0), *pm1ave(0), *pm10(0), *pm11(0), *pm12(0);
  if (fm1) {
    pm1tag = (TProfile*)fm1->Get(Form("%s/pmpfx",cp)); //assert(pm1tag);
    pm1pro = (TProfile*)fm1->Get(Form("%s/pmpfy",cp)); //assert(pm1pro);
    pm1ave = (TProfile*)fm1->Get(Form("%s/pmpfz",cp)); //assert(pm1ave);

    pm10 = (TProfile*)fm1->Get(Form("%s/pmpf",cp)); //assert(pm10);
    pm11 = (TProfile*)fm1->Get(Form("%s/pmpf1",cp)); //assert(pm11);
  pm12 = (TProfile*)fm1->Get(Form("%s/pmpf2",cp)); //assert(pm12);
  }

  TH1D *h = new TH1D("h",";p_{T} (GeV);MPF",3970,30,4000);
  //if (mode=="BB")  { h->SetMinimum(0.93); h->SetMaximum(1.08); }
  if (mode=="BB")  { h->SetMinimum(0.80); h->SetMaximum(1.20); }
  //if (mode=="EC0") { h->SetMinimum(0.80); h->SetMaximum(1.25); }
  if (mode=="EC0") { h->SetMinimum(0.80); h->SetMaximum(1.20); }
  if (mode=="EC1") { h->SetMinimum(0.80); h->SetMaximum(1.20); }
  if (mode=="EC2") { h->SetMinimum(0.00); h->SetMaximum(2.00); }
  if (mode=="HF")  { h->SetMinimum(0.00); h->SetMaximum(2.00); }
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  
  //lumi_13TeV = "Run2016GH, 16.8 fb^{-1}";
  //lumi_13TeV = "Run2017DE, 13.5 fb^{-1}";
  //lumi_13TeV = "Run2018D, 29.2 fb^{-1}";
  //lumi_13TeV = "Run2018A, 14.0 fb^{-1}";
  lumi_13TeV = "Run2018B, 7.1 fb^{-1}";
  //lumi_13TeV = "Run2018C, 6.9 fb^{-1}";
  //lumi_13TeV = "Run2018D, 29.2 fb^{-1}";
  TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);
  gPad->SetLogx();

  if (phwave) {
    tdrDraw(phwave,"HISTL",kNone,kGreen+2,kSolid,-1,kNone); // MPF vs pTave
    tdrDraw(phwtag,"HISTL",kNone,kRed,kSolid,-1,kNone); // MPF vs pTtag
    tdrDraw(phwpro,"HISTL",kNone,kBlue,kSolid,-1,kNone); // MPF vs pTprobe
  }
  if (pm1ave) {
    tdrDraw(pm1ave,"HISTL",kNone,kGreen+2,kDotted,-1,kNone); // MPF vs pTave
    tdrDraw(pm1tag,"HISTL",kNone,kRed,kDotted,-1,kNone); // MPF vs pTtag
    tdrDraw(pm1pro,"HISTL",kNone,kBlue,kDotted,-1,kNone); // MPF vs pTprobe
  }    

  tdrDraw(pmcave,"PzE0X0",kOpenDiamond,kGreen+2); // MPF vs pTave
  tdrDraw(pmctag,"PzE0X0",kOpenTriangleDown,kRed); // MPF vs pTtag
  tdrDraw(pmcpro,"PzE0X0",kOpenTriangleUp,kBlue); // MPF vs pTprobe

  tdrDraw(pdtave,"PzE0X0",kFullDiamond,kGreen+2); // MPF vs pTave
  tdrDraw(pdttag,"PzE0X0",kFullTriangleDown,kRed); // MPF vs pTtag
  tdrDraw(pdtpro,"PzE0X0",kFullTriangleUp,kBlue); // MPF vs pTprobe

  TLatex *tex = new TLatex();
  tex->SetTextSize(0.045);
  tex->SetNDC();
  tex->DrawLatex(0.20,0.75,cm);

  TLegend *legmc = tdrLeg(0.48,0.70,0.83,0.90);
  legmc->SetHeader("MC");
  legmc->AddEntry(pmctag,"","PL");
  legmc->AddEntry(pmcave,"","PL");
  legmc->AddEntry(pmcpro,"","PL");

  TLegend *legdt = tdrLeg(0.55,0.70,0.90,0.90);
  legdt->SetHeader("Data");
  legdt->AddEntry(pdttag,"vs tag","P");
  legdt->AddEntry(pdtave,"vs average","P");
  legdt->AddEntry(pdtpro,"vs probe","P");

  if (phwave) {
    phwtag->GetXaxis()->SetRangeUser(30,4000);
    phwave->GetXaxis()->SetRangeUser(30,4000);
    phwpro->GetXaxis()->SetRangeUser(30,4000);
  }
  if (pm1ave) {
    pm1tag->GetXaxis()->SetRangeUser(30,4000);
    pm1ave->GetXaxis()->SetRangeUser(30,4000);
    pm1pro->GetXaxis()->SetRangeUser(30,4000);
  }
  pmctag->GetXaxis()->SetRangeUser(30,4000);
  pmcave->GetXaxis()->SetRangeUser(30,4000);
  pmcpro->GetXaxis()->SetRangeUser(30,4000);
  pdttag->GetXaxis()->SetRangeUser(30,4000);
  pdtave->GetXaxis()->SetRangeUser(30,4000);
  pdtpro->GetXaxis()->SetRangeUser(30,4000);

  //c1->SaveAs(Form("pdf/drawMultijetJER_binvars_16GH_%s.pdf",cm));
  //c1->SaveAs(Form("pdf/drawMultijetJER_binvars_17DE_%s.pdf",cm));
  //c1->SaveAs(Form("pdf/drawMultijetJER_binvars_18D_%s.pdf",cm));
  //c1->SaveAs(Form("pdf/drawMultijetJER_binvars_18AV8_%s.pdf",cm));
  c1->SaveAs(Form("pdf/drawMultijetJER_binvars_18BV8_%s.pdf",cm));
  //c1->SaveAs(Form("pdf/drawMultijetJER_binvars_18CV8_%s.pdf",cm));
  //c1->SaveAs(Form("pdf/drawMultijetJER_binvars_18DV8_%s.pdf",cm));

  TH1D *h2 = new TH1D("h2",";p_{T} (GeV);MPF",3900,100,4000);
  h2->SetMinimum(0.90);
  h2->SetMaximum(1.05);
  h2->GetXaxis()->SetMoreLogLabels();
  h2->GetXaxis()->SetNoExponent();
  
  TCanvas *c2 = tdrCanvas("c2",h2,4,11,kSquare);
  gPad->SetLogx();

  tdrDraw(pmcpro,"Pz",kOpenTriangleUp,kBlack);
  tdrDraw(pmc0,"Pz",kOpenTriangleUp,kRed);
  tdrDraw(pmc1,"Pz",kOpenTriangleUp,kGreen+2);
  tdrDraw(pmc2,"Pz",kOpenTriangleUp,kBlue);

  tdrDraw(pdtpro,"Pz",kFullTriangleUp,kBlack);
  tdrDraw(pdt0,"Pz",kFullTriangleUp,kRed);
  tdrDraw(pdt1,"Pz",kFullTriangleUp,kGreen+2);
  tdrDraw(pdt2,"Pz",kFullTriangleUp,kBlue);

  tex->DrawLatex(0.20,0.75,cm);

  TLegend *legmc2 = tdrLeg(0.48,0.65,0.83,0.90);
  legmc2->SetHeader("MC");
  legmc2->AddEntry(pmcpro,"","P");
  legmc2->AddEntry(pmc0,"","P");
  legmc2->AddEntry(pmc1,"","P");
  legmc2->AddEntry(pmc2,"","P");

  TLegend *legdt2 = tdrLeg(0.55,0.65,0.90,0.90);
  legdt2->SetHeader("Data");
  legdt2->AddEntry(pdtpro,"Type-1 (dijet)","P");
  legdt2->AddEntry(pdt0,"Type-0","P");
  legdt2->AddEntry(pdt1,"Type-1","P");
  legdt2->AddEntry(pdt2,"Type-2","P");

  pmcpro->GetXaxis()->SetRangeUser(100,4000);
  pmc0->GetXaxis()->SetRangeUser(100,4000);
  pmc1->GetXaxis()->SetRangeUser(100,4000);
  pmc2->GetXaxis()->SetRangeUser(100,4000);
  pdtpro->GetXaxis()->SetRangeUser(100,4000);
  pdt0->GetXaxis()->SetRangeUser(100,4000);
  pdt1->GetXaxis()->SetRangeUser(100,4000);
  pdt2->GetXaxis()->SetRangeUser(100,4000);

  //c2->SaveAs(Form("pdf/drawMultijetJER_mettypes_16GH_%s.pdf",cm));
  //c2->SaveAs(Form("pdf/drawMultijetJER_mettypes_17DE_%s.pdf",cm));
  //c2->SaveAs(Form("pdf/drawMultijetJER_mettypes_18D_%s.pdf",cm));
  //c2->SaveAs(Form("pdf/drawMultijetJER_mettypes_18AV8_%s.pdf",cm));
  c2->SaveAs(Form("pdf/drawMultijetJER_mettypes_18BV8_%s.pdf",cm));
  //c2->SaveAs(Form("pdf/drawMultijetJER_mettypes_18CV8_%s.pdf",cm));
  //c2->SaveAs(Form("pdf/drawMultijetJER_mettypes_18DV8_%s.pdf",cm));
}
