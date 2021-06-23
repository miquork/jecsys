// Purpose: Use minitools/hadW.C outputs to estimate semileptonic BR in HF jets
//          Provided are electron and muon fractions and jets with e/mu
// Assumption is that muons are well correlated with SL BR, while electrons
// have some gamma>e+e- contamination, but otherwise flavors are symmetric for B
// Charm may have some flavor asymmetry due to higher muon mass relative to D's
//
// run with 'root -l -b -q minitools/drawSemilepBR.C+g'
#include "TFile.h"
#include "TProfile.h"
#include "TF1.h"
#include "TH2D.h"

#include "../tdrstyle_mod15.C"

//string mode = "e";
string _mode = "mu";
//string mode = "zmu";
//string _mode = "ze";

void drawSemilepBR(string mode = _mode) {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  //TFile *fd = new TFile("rootfiles/hadWUL1718V5_EMUF.root","READ");
  //TFile *fd = new TFile("rootfiles/hadWUL1718V5_Glu_v6.root","READ");
  TFile *fd = new TFile("rootfiles/hadWUL17V5New_JEC.root","READ");
  assert(fd && !fd->IsZombie());
  //TFile *fm = new TFile("rootfiles/hadWMC1718V5_EMUF.root","READ");
  //TFile *fm = new TFile("rootfiles/hadWMC1718V5_Glu_v6.root","READ");
  TFile *fm = new TFile("rootfiles/hadWMC17V5New_JEC.root","READ");
  assert(fm && !fm->IsZombie());

  //TFile *fz = new TFile("rootfiles/jme_bplusZ_merged_v34.root","READ");
  //TFile *fz = new TFile("rootfiles/jme_bplusZ_merged_v35.root","READ");
  //TFile *fz = new TFile("rootfiles/jme_bplusZ_merged_v35.root","READ");
  TFile *fz = new TFile("rootfiles/jme_bplusZ_merged_v35_2017_emu_wTTJets.root","READ");
  assert(fz && !fz->IsZombie());

  curdir->cd();

  vector<string> vf;
  if (mode=="e") {
    vf.push_back("belf");
    vf.push_back("qelf");
  }
  if (mode=="mu") {
    vf.push_back("bmuf");
    vf.push_back("qmuf");
  }
  if (mode=="ze") {
    vf.push_back("zbelf");
    vf.push_back("zqelf");
  }
  if (mode=="zmu") {
    vf.push_back("zbmuf");
    vf.push_back("zqmuf");
  }
  vector<string> va;
  //va.push_back("a");
  va.push_back("b");
  va.push_back("n");

  map<string,int> color;
  color["qelf"] = kBlue;
  color["qmuf"] = kGreen+2;
  color["belf"] = kOrange+2;
  color["bmuf"] = kRed;
  color["zqmuf"] = kGreen+3;
  color["zbmuf"] = kRed+1;
  color["zqelf"] = kBlue+1;
  color["zbelf"] = kOrange+3;

  map<string,int> marker;
  marker["qelf"] = kOpenSquare;
  marker["qmuf"] = kFullSquare;
  marker["belf"] = kOpenCircle;
  marker["bmuf"] = kFullCircle;
  marker["zqmuf"] = kFullStar;
  marker["zbmuf"] = kFullDiamond;
  marker["zqelf"] = kOpenStar;
  marker["zbelf"] = kOpenDiamond;
  //
  map<string,int> omarker;
  omarker["qelf"] = kOpenSquare;
  omarker["qmuf"] = kOpenSquare;
  omarker["belf"] = kOpenCircle;
  omarker["bmuf"] = kOpenCircle;
  omarker["zqmuf"] = kOpenStar;
  omarker["zbmuf"] = kOpenDiamond;
  omarker["zqelf"] = kOpenStar;
  omarker["zbelf"] = kOpenDiamond;
  //
  map<string,int> fmarker;
  fmarker["qelf"] = kFullSquare;
  fmarker["qmuf"] = kFullSquare;
  fmarker["belf"] = kFullCircle;
  fmarker["bmuf"] = kFullCircle;
  fmarker["zqmuf"] = kFullStar;
  fmarker["zbmuf"] = kFullDiamond;
  fmarker["zqelf"] = kFullStar;
  fmarker["zbelf"] = kFullDiamond;

  map<string,int> style;
  style["a"] = kDotted;//kSolid;
  style["b"] = kDashed;
  style["n"] = kSolid;//kDotted;

  map<string,map<string, string> > label;
  /*
  label["bmuf"]["b"] = "#mu E frac. in #mu-tagged b jets";
  label["qmuf"]["b"] = "#mu E frac. in #mu-tagged W>qq' jets ";
  label["bmuf"]["n"] = "#mu-tagged frac. of b tagged jets";
  label["qmuf"]["n"] = "#mu-tagged frac. of W>qq' tagged jets";
  */
  label["bmuf"]["b"] = "#mu E frac. (#mu-tag b jet)";
  label["qmuf"]["b"] = "#mu E frac. (#mu-tag W>qq') ";
  label["bmuf"]["n"] = "#mu-tag frac. (b jet)";
  label["qmuf"]["n"] = "#mu-tag frac. (W>qq')";

  label["zbmuf"]["b"] = "#mu E frac. (#mu-tag b jet)";
  label["zqmuf"]["b"] = "#mu E frac. (#mu-tag Z+jet) ";
  label["zbmuf"]["n"] = "#mu-tag frac. (b jet)";
  label["zqmuf"]["n"] = "#mu-tag frac. (Z+jet)";

  label["belf"]["b"] = "e E frac. (e-tag b jet)";
  label["qelf"]["b"] = "e E frac. (e-tag W>qq') ";
  label["belf"]["n"] = "e-tag frac. (b jet)";
  label["qelf"]["n"] = "e-tag frac. (W>qq')";

  //TH1D *h = tdrHist("h","Fraction",0.002,2,"p_{T} (GeV)",30,230);
  //TH1D *h = tdrHist("h","Fraction",0.04,0.4,"p_{T} (GeV)",30,230); // w/o a 
  //TH1D *h = tdrHist("h","Fraction",0.08,0.4,"p_{T} (GeV)",30,230); // w/o q, el 

  TH1D *hu = tdrHist("hu","Fraction",0.04+1e-5,0.4-1e-5,"p_{T} (GeV)",30,230);
  //TH1D *hd = tdrHist("hd","Data-MC (0.01)",-2,+0.5,"p_{T} (GeV)",30,230); //mu
  TH1D *hd = tdrHist("hd","Data-MC (0.01)",-3.0,+1.0,"p_{T,jet} (GeV)",30,230); //el

  hu->GetYaxis()->SetMoreLogLabels();
  hu->GetYaxis()->SetNoExponent();

  //lumi_13TeV = "TT lepton+jet, UL17+18, 101.4 fb^{-1}";
  lumi_13TeV = "[TT lepton+jet] 2017, 41.5 fb^{-1}";
  if (mode=="zmu"||mode=="ze")
    lumi_13TeV = "[Z+jet] 2017, 41.5 fb^{-1}";
  //TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);
  TCanvas *c1 = tdrDiCanvas("c1",hu,hd,4,11);


  c1->cd(1);
  gPad->SetLogx();
  gPad->SetLogy();

  TLegend *legm = tdrLeg(0.40,0.90-5*0.05,0.70,0.90);
  //legm->SetTextSize(0.040);
  legm->SetHeader("MC");
  TLegend *legd = tdrLeg(0.48,0.90-5*0.05,0.78,0.90);
  //legd->SetTextSize(0.040);
  legd->SetHeader("Data");

  c1->cd(2);
  gPad->SetLogx();

  for (int ia = 0; ia != va.size(); ++ia) {
    for (int jf = 0; jf != vf.size(); ++jf) {
      
      string &sa = va[ia];
      string &sf = vf[jf];
      const char *ca = sa.c_str();
      const char *cf = sf.c_str();

      TH1D *hd(0), *hm(0);
      if (sf=="zqmuf" || sf=="zbmuf" || sf=="zqelf" || sf=="zbelf") {
	TH2D *h2d = (TH2D*)fz->Get(Form("data/eta_00_25/"
					"h_JetPt_%sEF%s_alpha100",
					(sf=="zqmuf" || sf=="zbmuf") ?
					"mu" : "chEm",
					(sf=="zbmuf" || sf=="zbelf") ?
					"_btagDeepBtight" : ""));
	assert(h2d);
	TH2D *h2m = (TH2D*)fz->Get(Form("mc/eta_00_25/"
					"h_JetPt_%sEF%s_alpha100",
					(sf=="zqmuf" || sf=="zbmuf") ?
					"mu" : "chEm",
					(sf=="zbmuf" || sf=="zbelf") ?
					"_btagDeepBtight" : ""));
	assert(h2m);

	if (sa=="n") {
	  hd = h2d->ProjectionX(Form("hd%s%s",cf,ca),2,
				h2d->GetNbinsY());
	  TH1D *h1d = h2d->ProjectionX(Form("h1d%s%s",cf,ca),1,
					h2d->GetNbinsY());
	  hd->Divide(hd,h1d,1,1,"B");

	  hm = h2m->ProjectionX(Form("hm%s%s",cf,ca),2,
				h2m->GetNbinsY());
	  TH1D *h1m  = h2m->ProjectionX(Form("h1m%s%s",cf,ca),1,
					h2m->GetNbinsY());
	  hm->Divide(hm,h1m,1,1,"B");
	}
	else if (sa=="b") {
	  TProfile *pd = h2d->ProfileX(Form("pd%s%s",cf,ca),2,
				       h2d->GetNbinsY());
	  TProfile *pm = h2m->ProfileX(Form("pm%s%s",cf,ca),2,
				       h2m->GetNbinsY());
	  hd = pd->ProjectionX(Form("hd%s%s",cf,ca));
	  hm = pm->ProjectionX(Form("hm%s%s",cf,ca));
	}
	else
	  assert(false);
      }
      else {
	TProfile *pd = (TProfile*)fd->Get(Form("p%s%s",cf,ca)); assert(pd);
	TProfile *pm = (TProfile*)fm->Get(Form("p%s%s",cf,ca)); assert(pm);
	hd = pd->ProjectionX(Form("hd%s%s",cf,ca));
	hm = pm->ProjectionX(Form("hm%s%s",cf,ca));
      }
      assert(hd);
      assert(hm);

      // Scale quark fraction from 27% charm to effective 100%
      //if ((sf=="qmuf" || sf=="qelf") && sa=="n") {
      //if ((sf=="qmuf") && sa=="n") {
      if (false) {
	hd->Scale(100./27.);
	hm->Scale(100./27.);
      }

      // Calculate data-MC difference in percent-units
      //TH1D *hdm = pd->ProjectionX(Form("h%s%s",cf,ca));
      TH1D *hdm = (TH1D*)hd->Clone(Form("h%s%s",cf,ca));
      hdm->Add(hd,hm,100,-100);
      
      c1->cd(1);
      tdrDraw(hm,"HIST][",kNone,color[sf],style[sa],-1,kNone);
      tdrDraw(hd,"Pz",sa=="b" ? omarker[sf] : fmarker[sf],color[sf],kSolid,
	      -1, kNone);
      
      legm->AddEntry(hm," ","L");
      legd->AddEntry(hd,label[sf][sa].c_str(),"PLE");
      //legd->AddEntry(hd,label["bmuf"]["b"].c_str(),"PLE");
      //legd->AddEntry(hd,"test","PLE");
      
      c1->cd(2);
      tdrDraw(hdm,"Pz",sa=="b" ? omarker[sf] : fmarker[sf],color[sf],kSolid);

      TF1 *f1 = new TF1(Form("f1%s%s",cf,ca),"[0]+[1]*log(0.01*x)",30,230);
      f1->FixParameter(1,0);
      hdm->Fit(f1,"QRN");
      f1->SetLineColor(color[sf]);
      f1->SetLineStyle(style[sa]);
      f1->Draw("SAME");

    } // for jf
  } // for ia

  c1->SaveAs(Form("pdf/drawSemilepBR_2017_%s.pdf",mode.c_str()));
} // drawSemilepBR
