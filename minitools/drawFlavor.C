// Purpose: Draw results using Flavor.C/h
// Run with 'root -l -b -q minitools/mk_drawFlavor.C'
#include "TFile.h"
#include "TLine.h"

#include "../Flavor.h"
#include "../tdrstyle_mod15.C"
#include "tools.C"

// Data/MC false (i.e. MC) is good for checking HDM vs MPF vs MJB
bool doDataMC = true;//false;

void drawFlavor() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  Flavor *f = new Flavor;
  f->doRobin = false;

  TFile *fmj = new TFile("rootfiles/multijet_Rebin2_20201207_UL2018ABCD_jecV5_jerV2.root","READ");
  //TFile *fmj = new TFile("rootfiles/multijet_Rebin2_20201209_UL2017BCDEF_jecV6_jerV3.root","READ");
  assert(fmj && !fmj->IsZombie());
  TH1D *hmj = (TH1D*)fmj->Get("MC/Pt30/MPF_ptave_MG"); assert(hmj);
  TH1D *hjb = (TH1D*)fmj->Get("MC/Pt30/MJB_ptave_MG"); assert(hjb);
  TH1D *hmu = (TH1D*)fmj->Get("MC_JER-down/Pt30/MPF_ptave_MG"); assert(hmu);
  //
  if (doDataMC) {
    TH1D *hmjd = (TH1D*)fmj->Get("Data/Pt30/MPF_ptave_L1L2L3Res"); assert(hmjd);
    TH1D *hjbd = (TH1D*)fmj->Get("Data/Pt30/MJB_ptave_L1L2L3Res"); assert(hjbd);
    TH1D *hmud = (TH1D*)hmjd->Clone("hmud");
    hmjd->Divide(hmj);
    hjbd->Divide(hjb);
    hmud->Divide(hmu);
    hmj = hmjd;
    hjb = hjbd;
    hmu = hmud;
  }

  TFile *fgf = new TFile("rootfiles/jecdataRun2Test.root","READ");
  //TFile *fgf = new TFile("rootfiles/jecdata2018ABCD.root","READ");
  // TFile *fgf = new TFile("rootfiles/jecdata2017BCDEF.root","READ");
  assert(fgf && !fgf->IsZombie());
  TH1D *hdm = (TH1D*)fgf->Get("mc/eta00-13/hdm_mpfchs1_multijet"); assert(hdm);
  if (doDataMC) {
    hdm = (TH1D*)fgf->Get("ratio/eta00-13/hdm_mpfchs1_multijet"); assert(hdm);
  }
  TH1D *hdu = (TH1D*)hdm->Clone("hdu");
  hdu->Multiply(hmu);
  hdu->Divide(hmj);

  // Retrieve results before/after flavor corrections. Only data+MC, no ratio
  TGraphErrors *gmjbdb(0), *gmjbda(0), *gmjbm(0);
  TGraphErrors *gmpfdb(0), *gmpfda(0), *gmpfm(0);
  fgf->cd("data/eta00-13/orig");
  TDirectory *db = gDirectory;
  fgf->cd("data/eta00-13");
  TDirectory *da = gDirectory;
  fgf->cd("mc/eta00-13/orig");
  TDirectory *dm = gDirectory;
  gmjbdb = (TGraphErrors*)db->Get("ptchs_multijet_a30"); assert(gmjbdb);
  gmjbda = (TGraphErrors*)da->Get("ptchs_multijet_a30"); assert(gmjbda);
  gmjbm  = (TGraphErrors*)dm->Get("ptchs_multijet_a30"); assert(gmjbm);
  gmpfdb = (TGraphErrors*)db->Get("mpfchs1_multijet_a30"); assert(gmpfdb);
  gmpfda = (TGraphErrors*)da->Get("mpfchs1_multijet_a30"); assert(gmpfda);
  gmpfm  = (TGraphErrors*)dm->Get("mpfchs1_multijet_a30"); assert(gmpfm);
  if (doDataMC) {
    gmjbda = tools::ratioGraphs(gmjbda, gmjbm);
    gmjbdb = tools::ratioGraphs(gmjbdb, gmjbm);
    gmpfda = tools::ratioGraphs(gmpfda, gmpfm);
    gmpfdb = tools::ratioGraphs(gmpfdb, gmpfm);
  }
  else { // useMC 
    gmjbdb = gmjbm;
    gmjbda = gmjbm;
    gmpfdb = gmpfm;
    gmpfda = gmpfm;
  }
  TGraphErrors *gaob = tools::ratioGraphs(gmpfda,gmpfdb);
  TGraphErrors *ghdmda = new TGraphErrors(hdm);
  TGraphErrors *ghdmdb = tools::ratioGraphs(ghdmda,gaob);

  curdir->cd();

  //lumi_13TeV = "Flavor.C / Autumn18_V3";
  lumi_13TeV = "Flavor.C / Run2Test";
  TH1D *h = tdrHist("h","Flavor response (MC)",0.96,1.06);
  if (doDataMC) h->SetYTitle("Flavor response (Data/MC)");
  if (f->doRobin) {
    h = tdrHist("hr","Correction",0.96,1.34,"p_{T}^{reco} (GeV)",7,6100);
  }
  TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);
  gPad->SetLogx();
  
  TLine *l = new TLine();
  l->SetLineStyle(kDotted);
  l->DrawLine(15,1,3500,1);

  TH1D *hl = tdrHist("hl","Leading jet");
  TH1D *hr = tdrHist("hr","Recoil");
  TH1D *hmjb = tdrHist("hmjb","Multijet balance");

  TH1D *hzl = tdrHist("hzl","Z+jet (pT,lead)");
  TH1D *hzr = tdrHist("hzr","Z+jet (pT,recoil)");
  TH1D *hzmjb = tdrHist("hzmjb","Multijet balance wrt Z+jet");

  TH1D *hq = tdrHist("hq","Quark");
  TH1D *hs = tdrHist("hs","Strange");
  TH1D *hc = tdrHist("hc","Charm");
  TH1D *hb = tdrHist("hb","Bottom");
  TH1D *hg = tdrHist("hg","Gluon");

  // Robin's plot style
  if (f->doRobin) {
    hq = tdrHist("hqr","Quark",0,1,"p_{T} (GeV)",9,5100);
    hs = tdrHist("hsr","Strange",0,1,"p_{T} (GeV)",9,5100);
    hc = tdrHist("hcr","Charm",0,1,"p_{T} (GeV)",9,5100);
    hb = tdrHist("hbr","Bottom",0,1,"p_{T} (GeV)",9,5100);
    hg = tdrHist("hgr","Gluon",0,1,"p_{T} (GeV)",9,5100);
  }

  for (int i = 1; i != hl->GetNbinsX()+1; ++i) {

    double pt = hl->GetBinCenter(i);
    if (f->doRobin) {
      double offset = 0.5*0.5*20;
      pt = hl->GetBinCenter(i) - offset;
    }
    double pta = pt;
    double ptl = pt;
    double ptr = 0.45*pt;

    double eta = 0.;
    double absetamin = 0;
    double absetamax = 1.3;
    double w = 0; // Pythia8
    if (doDataMC) {
      w = -1; // data
    }
    //double w = 1; // Herwig7
    if (f->doRobin) { 
      w = 0;
      absetamin = 0;
      absetamax = 0.783;
    }
    //if (pt>1650) continue;

    //double rl = f->getResp(pt, eta, "MultijetLeading13", 0);
    double rl = f->getResp(ptl, absetamin, absetamax, "MultijetLeading13",
			   w, "13");
    hl->SetBinContent(i, rl);
    // Error used in hl,hr ratio and f1mjb fit
    hl->SetBinError(i, 0.0002);
    //
    double rzl = f->getResp(ptl, absetamin, absetamax, "ZJet13", w, "13");
    hzl->SetBinContent(i, rzl);
    hzl->SetBinError(i, 0.0002);

    //cout << "pt="<<pt<<" rl="<<r<<endl;
    //double rr = f->getResp(pt, eta, "MultijetRecoil25", 0); // NOPE!
    // => Retrieve fractions at pT,recoil, but response at Crecoil*pT,recoil
    double ff[5]; f->getFracs(pta, eta, "MultijetRecoil25", ff);
    double rr = f->getResp(ptr, absetamin, absetamax, ff, w, "25");
    hr->SetBinContent(i, rr);
    hr->SetBinError(i, 0);
    //
    double rzr = f->getResp(ptr, absetamin, absetamax, "ZJet13", w, "13");
    hzr->SetBinContent(i, rzr);
    hzr->SetBinError(i, 0.0002);

    double rq = f->getResp(pt, absetamin, absetamax, "ud", w);
    hq->SetBinContent(i, rq);
    double rs = f->getResp(pt, absetamin, absetamax, "s", w);
    hs->SetBinContent(i, rs);
    double rc = f->getResp(pt, absetamin, absetamax, "c", w);
    hc->SetBinContent(i, rc);
    double rb = f->getResp(pt, absetamin, absetamax, "b", w);
    hb->SetBinContent(i, rb);
    double rg = f->getResp(pt, absetamin, absetamax, "g", w);
    hg->SetBinContent(i, rg);
  } // for i

  // Robin Aggleton's colors (update_14_6_19.pdf)
  if (f->doRobin) {
    tdrDraw(hq,"PL",kNone,kRed,kSolid,-1,kNone);
    tdrDraw(hg,"PL",kNone,kAzure,kSolid,-1,kNone);
    tdrDraw(hs,"PL",kNone,kBlue,kSolid,-1,kNone);
    tdrDraw(hc,"PL",kNone,kGreen+2,kSolid,-1,kNone);
    tdrDraw(hb,"PL",kNone,kOrange+1,kSolid,-1,kNone);
  }
  else if (false) {
    // Btag colors
    tdrDraw(hq,"PL",kNone,kBlue,kDotted,-1,kNone);
    tdrDraw(hs,"PL",kNone,kCyan+2,kDotted,-1,kNone);
    tdrDraw(hc,"PL",kNone,kGreen+2,kDotted,-1,kNone);
    tdrDraw(hb,"PL",kNone,kRed,kDotted,-1,kNone);
    tdrDraw(hg,"PL",kNone,kOrange+2,kDotted,-1,kNone);
  }

  // Multijet responses
  tdrDraw(hl,"PL",kNone,kGreen+2,kSolid,-1,kNone);
  tdrDraw(hr,"PL",kNone,kRed,kSolid,-1,kNone);
  hmjb->Divide(hl,hr);
  tdrDraw(hmjb,"PL",kNone,kBlack,kSolid,-1,kNone);

  tdrDraw(hzl,"PL",kNone,kOrange+2,kSolid,-1,kNone);
  tdrDraw(hzr,"PL",kNone,kOrange+1,kSolid,-1,kNone);
  hzmjb->Divide(hl,hr);
  hzmjb->Divide(hzl);
  hzmjb->Multiply(hzr);
  tdrDraw(hzmjb,"PL",kNone,kOrange+3,kSolid,-1,kNone);

  //tdrDraw(hjb,"PL",kOpenDiamond,kBlack);
  //tdrDraw(hmj,"PL",kOpenCircle,kBlack);
  tdrDraw(hdm,"PL",kFullCircle,kGray);//kBlack);
  //tdrDraw(hdu,"PL",kFullCircle,kBlue); hdu->SetMarkerSize(0.7);
  
  // New graphs
  tdrDraw(gmjbdb,"Pz",kOpenDiamond,kBlack);//kGray+2);
  tdrDraw(gmpfdb,"Pz",kOpenCircle,kBlack);//kGray+2);
  tdrDraw(ghdmdb,"Pz",kFullCircle,kBlack);//kGray+2);

  // Newly fixed HDM
  TGraphErrors *glor = new TGraphErrors(hmjb);
  TGraphErrors *ghdmda2 = tools::ratioGraphs(ghdmdb,glor);
  tdrDraw(ghdmda2,"PL",kFullSquare,kRed);
  ghdmda2->SetMarkerSize(0.7);

  double dx = (doDataMC ? +0.05 : 0.);
  TLegend *leg = tdrLeg(0.4+dx,0.90-6*0.05,0.6+dx,0.90);
  leg->AddEntry(hdm,"Multijet HDM","PLE");
  //leg->AddEntry(hdu,"Multijet HDM (JER down)","PLE");
  //leg->AddEntry(hmj,"Multijet MPF","PLE");
  leg->AddEntry(gmpfdb,"Multijet MPF","PLE");
  //leg->AddEntry(hjb,"Multijet MJB","PLE");
  leg->AddEntry(gmjbdb,"Multijet MJB","PLE");
  leg->AddEntry(hmjb,"Leading jet / recoil","L");
  leg->AddEntry(hl,"Leading jet","L");
  leg->AddEntry(hr,"Recoil","L");

  TLegend *leg2 = tdrLeg(0.4+dx,0.20,0.6+dx,0.20+3*0.05);
  leg2->AddEntry(hzmjb,"Leading jet / recoil wrt Z+jet","L");
  leg2->AddEntry(hzl,"Z+jet (pT,lead)","L");
  leg2->AddEntry(hzr,"Z+jet (pT,recoil)","L");

  TF1 *f1mjb = new TF1("f1mjb","[0]+[1]*pow(x,[2])",114,2500);  
  f1mjb->SetParameters(1,1,-1);
  //hmjb->Fit(f1mjb,"QRNW");
  hmjb->Fit(f1mjb,"QRN");
  f1mjb->SetLineWidth(2);
  f1mjb->SetLineStyle(kDotted);
  f1mjb->SetLineColor(kOrange+2);
  f1mjb->DrawClone("SAME");
  f1mjb->SetLineWidth(1);
  f1mjb->SetRange(15,3500);
  f1mjb->DrawClone("SAME");

  cout << "  // MultijetRecoilScale from minitools/drawFlavor.C" << endl;
  cout << Form("  TF1 *%s = new TF1(\"%s\",\"%s\",114,2500);",
	       f1mjb->GetName(),f1mjb->GetName(),
	       f1mjb->GetExpFormula().Data()) << endl;
  cout << Form("  %s->SetParameters(",f1mjb->GetName());
  for (int i = 0; i != f1mjb->GetNpar(); ++i) {
    cout << Form("%1.4g%s",f1mjb->GetParameter(i),
		 i==f1mjb->GetNpar()-1 ? ");" : ", ");
  }
  cout << Form(" // chi2/NDF=%1.1f/%d",f1mjb->GetChisquare(),
	       f1mjb->GetNDF()) << endl;


  TF1 *f1zmjb = new TF1("f1zmjb","[0]+[1]*pow(x,[2])",114,2500);  
  f1zmjb->SetParameters(1,1,-1);
  hzmjb->Fit(f1zmjb,"QRN");
  f1zmjb->SetLineWidth(2);
  f1zmjb->SetLineStyle(kDotted);
  f1zmjb->SetLineColor(kOrange+3);
  f1zmjb->DrawClone("SAME");
  f1zmjb->SetLineWidth(1);
  f1zmjb->SetRange(15,3500);
  f1zmjb->DrawClone("SAME");

  cout << "  // MultijetRecoilScale from minitools/drawFlavor.C" << endl;
  cout << Form("  TF1 *%s = new TF1(\"%s\",\"%s\",114,2500);",
	       f1zmjb->GetName(),f1zmjb->GetName(),
	       f1zmjb->GetExpFormula().Data()) << endl;
  cout << Form("  %s->SetParameters(",f1zmjb->GetName());
  for (int i = 0; i != f1zmjb->GetNpar(); ++i) {
    cout << Form("%1.4g%s",f1zmjb->GetParameter(i),
		 i==f1zmjb->GetNpar()-1 ? ");" : ", ");
  }
  cout << Form(" // chi2/NDF=%1.1f/%d",f1zmjb->GetChisquare(),
	       f1zmjb->GetNDF()) << endl;

  gPad->RedrawAxis();
  c1->SaveAs("pdf/drawFlavor.pdf");
}
