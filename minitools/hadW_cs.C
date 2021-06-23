// Purpose: test separating W>qq' into ud and cs enriched regions,
//          to estimate c and s quark JES relative to ud
#include "TFile.h"
#include "TH2D.h"
#include "TChain.h"
#include "TCut.h"
#include "TF1.h"
#include "TProfile2D.h"
#include "TBox.h"

#include "../tdrstyle_mod15.C"

// Fit function for ABCD fractions
// x=bin index in TH2D
// p[0]=light c-mistag rate modifier (worse)
// p[1]=charm c-efficiency modifier (worse)
// p[2]=light+charm g-mistag rate modified (worse)
// p[3]=gluon (others) rate modifier (worse), e.g. due to FSR combinatorics
// Input MC counts are in h2m
// Output data counts are in h2d
TH2D *h2nm(0), *h2nd(0), *h2nf(0);
Double_t fnABCD(Double_t *x, Double_t *p);

TH2D *h2mm(0), *h2md(0), *h2mf(0);
Double_t frABCD(Double_t *x, Double_t *p);

bool fixGluonJES = true;

void hadW_cs() {
 
  setTDRStyle();
  TDirectory *curdir = gDirectory;

  //TFile *fd = new TFile("rootfiles/hadWUL18V5_MPDGcorrNoW.root","READ");
  TFile *fd = new TFile("rootfiles/hadWUL1718V5_MPDGcorrNoW.root","READ");
  assert(fd && !fd->IsZombie());
  //TFile *fm = new TFile("rootfiles/hadWMC18V5_MPDGcorrNoW.root","READ");
  TFile *fm = new TFile("rootfiles/hadWMC1718V5_MPDGcorrNoW.root","READ");
  assert(fm && !fm->IsZombie());

  curdir->cd();

  h2nd = (TH2D*)fd->Get("h2n"); assert(h2nd); h2nd->SetName("h2nd");
  h2nm = (TH2D*)fm->Get("h2n"); assert(h2nm); h2nm->SetName("h2nm");

  TProfile2D *p2md = (TProfile2D*)fd->Get("p2m"); assert(p2md);
  p2md->SetName("p2md");
  TProfile2D *p2mm = (TProfile2D*)fm->Get("p2m"); assert(p2mm);
  p2mm->SetName("p2mm");
  h2md = (TH2D*)p2md->ProjectionXY("h2md");
  h2mm = (TH2D*)p2mm->ProjectionXY("h2mm");

  // Map TH2D to TH1D for fitting, include parameter constraints if preferred
  const int nx = h2nd->GetNbinsX();
  const int np = 4;
  TF1 *fN = new TF1("fN",fnABCD,0.5-np,nx+0.5,np);
  //fN->SetParameters(0.18,0.15,0.15,0.05); // initial guess
  fN->SetParameters(1.88848e-01,1.65628e-01,1.43455e-01,2.45679e-02); // 1st fit
  //const int np = fN->GetNpar();
  TH1D *h1d = new TH1D("h1d","",nx+np,0.5-np,nx+0.5);
  // Set 5% uncertainty on parameters around starting value
  for (int i=0; i != np; ++i) {
    h1d->SetBinContent(h1d->FindBin(-i),fN->GetParameter(i));
    h1d->SetBinError(h1d->FindBin(-i),0.05);
  }
  for (int i = 1; i != h2nd->GetNbinsX()+1; ++i) {
    int j = 1;
    int k = nx*(j-1) + (i-1) + 1;
    assert(k<nx+1+np);
    int kk = h1d->FindBin(k);
    h1d->SetBinContent(kk, h2nd->GetBinContent(i,j));
    h1d->SetBinError(kk, h2nd->GetBinError(i,j));
  } // for i
  h2nf = (TH2D*)h2nd->Clone("h2nf");
  //fN->Eval(1); // set bins to default without fitting
  h1d->Fit(fN,"RN");

  TH1D *h1 = tdrHist("h1","Fraction of events",1e-2,3,//1.25,
		     "Tag region",-0.5,4.5);
  h1->GetXaxis()->SetBinLabel(1,"All");
  h1->GetXaxis()->SetBinLabel(2,"A");
  h1->GetXaxis()->SetBinLabel(3,"B");
  h1->GetXaxis()->SetBinLabel(4,"C");
  h1->GetXaxis()->SetBinLabel(5,"D");
  
  //lumi_13TeV = "2018, 59.9 fb^{-1}";
  //lumi_13TeV = "2017, 41.5 fb^{-1}";
  lumi_13TeV = "UL17+18, 101.4 fb^{-1}";
  TCanvas *c1 = tdrCanvas("c1",h1,4,11,kSquare);
  gPad->SetLogy();

  TLegend *leg1 = tdrLeg(0.70,0.90-5*0.05,1.00,0.90);

  int color[] = {kGray+2,kBlue,kGreen+2,kOrange+1,kCyan+2};
  int ncolor = sizeof(color)/sizeof(color[0]);
  assert(ncolor==h2nm->GetNbinsY());

  const char* name[] = {"MC","ud+us","cs+cd","gx","others"};
  int nname = sizeof(name)/sizeof(name[0]);
  assert(nname==h2nm->GetNbinsY());

  // Data
  double nd = h2nd->GetBinContent(1,1);
  double nm = h2nm->GetBinContent(1,1);
  TH1D *hnd = h2nd->ProjectionX("hnd",1,1);
  hnd->Scale(1./nd);
  //hnd->Scale(1./nm);
  tdrDraw(hnd,"HP",kFullCircle,kBlack,kSolid,-1,kNone);
  leg1->AddEntry(hnd,"Data","PL");

  // MC and fit
  for (int j = 1; j != h2nm->GetNbinsY()+1; ++j) {

    TH1D *hn = h2nm->ProjectionX(Form("hnm%d",j),j,j);
    hn->Scale(1./nm);
    //hn->SetBinContent(1,hn->GetBinContent(1)/4.);
    tdrDraw(hn,"HP",j==1 ? kOpenCircle : kNone,color[j-1],kSolid,-1,kNone);

    leg1->AddEntry(hn,name[j-1],"PL");

    TH1D *hnf = h2nf->ProjectionX(Form("hnf%d",j),j,j);
    hnf->Scale(1./nd);
    hnf->SetLineWidth(2);
    hnf->SetMarkerSize(0.5);
    tdrDraw(hnf,"HP",j==1 ? kOpenCircle : kNone,color[j-1],kDotted,-1,kNone);
  } // for j

  gPad->RedrawAxis();

  //gPad->SetLogy(kFALSE);
  //h2f->Draw("BOX");
  //h1d->Draw();

  //p2f = p2d->ProjectionXY("p2f"); // TH2D to be filled
  h2mf = (TH2D*)h2md->Clone("h2mf");
  TH1D *hmd = h2md->ProjectionX("p1d",1,1); // TH1D to be fit for data
  //TH1D *hmm = h2mm->ProjectionX("p1m",1,1); // TH1D to be fit for MC (test)
  TF1 *fr = new TF1("fr",frABCD,-0.5,4.5,4); // Fit function
  //fr->SetParameters(1,1,1,1); h2nf = h2nm; // MC reproduction
  double kd = h2md->GetBinContent(1,1)/h2mm->GetBinContent(1,1);
  double kg = 0.5*(0.994+1.000);
  fr->SetParameters(kd,kd,kd,kd*kg,kd);
  fr->FixParameter(3,kd);
  if (fixGluonJES) fr->FixParameter(2,kd*kg);
  //fr->Eval(0); // Fill in p2f
  //h2f = h2m; // Use MC purities instead of data
  hmd->Fit(fr,"RN");
  //hmm->Fit(fr,"RN"); // test fit to MC


  TH1D *h2 = tdrHist("h2","Mass-based response",0.98,1.05,
		     "Tag region",-0.5,4.5);
  h2->GetXaxis()->SetBinLabel(1,"All");
  h2->GetXaxis()->SetBinLabel(2,"A");
  h2->GetXaxis()->SetBinLabel(3,"B");
  h2->GetXaxis()->SetBinLabel(4,"C");
  h2->GetXaxis()->SetBinLabel(5,"D");
  
  TCanvas *c2 = tdrCanvas("c2",h2,4,11,kSquare);

  TLegend *leg2 = tdrLeg(0.70,0.90-5*0.05,1.00,0.90);

  // Data
  double md = h2md->GetBinContent(1,1);
  double mm = h2mm->GetBinContent(1,1);
  TH1D *hmd2 = h2md->ProjectionX("hd2",1,1);
  hmd2->Scale(1./mm);
  //hmd2->Scale(1./md);
  //hmd2->Scale(1./hd2->GetBinContent(1));
  tdrDraw(hmd2,"HP",kFullCircle,kBlack,kSolid,-1,kNone);
  leg2->AddEntry(hmd2,"Data","PL");


  // MC and fit
  for (int j = 1; j != h2mm->GetNbinsY()+1; ++j) {

    TH1D *hmm = h2mm->ProjectionX(Form("hmm%d2",j),j,j);
    hmm->Scale(1./mm);
    tdrDraw(hmm,"HP",j==1 ? kOpenCircle : kNone,color[j-1],kSolid,-1,kNone);

    leg2->AddEntry(hmm,name[j-1],"PL");

    TH1D *hmf = h2mf->ProjectionX(Form("hmf%d2",j),j,j);
    //hmf->Scale(1./md);
    hmf->Scale(1./mm);
    hmf->SetLineWidth(2);
    hmf->SetMarkerSize(0.5);
    tdrDraw(hmf,"HP",j==1 ? kOpenCircle : kNone,color[j-1],kDotted,-1,kNone);
  } // for j

  gPad->RedrawAxis();



  TH1D *h3 = tdrHist("h3","Flavor response vs ud (%)",-0.4,0.4,//-0.35,0.35,
		     "Jet flavor",-0.5,4.5);
  h3->GetXaxis()->SetBinLabel(1,"All");
  h3->GetXaxis()->SetBinLabel(2,"ud");
  h3->GetXaxis()->SetBinLabel(3,"cs");
  h3->GetXaxis()->SetBinLabel(4,"gx");
  h3->GetXaxis()->SetBinLabel(5,"others");

  TH1D *hr = (TH1D*)h3->Clone("hr"); hr->Reset();
  /*
  double dm = p2f->GetBinContent(1,2)/md-p2m->GetBinContent(1,2)/mm;
  hr->SetBinContent(1,p2f->GetBinContent(1,1)/md-p2m->GetBinContent(1,1)/mm-dm);
  hr->SetBinContent(2,p2f->GetBinContent(1,2)/md-p2m->GetBinContent(1,2)/mm-dm);
  hr->SetBinContent(3,p2f->GetBinContent(1,3)/md-p2m->GetBinContent(1,3)/mm-dm);
  hr->SetBinContent(4,p2f->GetBinContent(1,4)/md-p2m->GetBinContent(1,4)/mm-dm);
  hr->SetBinError(1,p2f->GetBinError(1,1)/md);
  hr->SetBinError(2,p2f->GetBinError(1,2)/md+0.0001);
  hr->SetBinError(3,p2f->GetBinError(1,3)/md);
  hr->SetBinError(4,p2f->GetBinError(1,4)/md);
  hr->Scale(100.);
  */
  //double dm = fr->GetParameter(0);
  double dm0 = (fr->GetParameter(0)*h2nf->GetBinContent(1,2) +
		fr->GetParameter(1)*h2nf->GetBinContent(1,3) +
		fr->GetParameter(2)*h2nf->GetBinContent(1,4) +
		fr->GetParameter(3)*h2nf->GetBinContent(1,5)) /
    (h2nf->GetBinContent(1,2) + h2nf->GetBinContent(1,3) +
     h2nf->GetBinContent(1,4) + h2nf->GetBinContent(1,5));
  double dm = fr->GetParameter(0); // ud as reference
  hr->SetBinContent(1,100.*(dm0-dm));
  hr->SetBinContent(2,100.*(fr->GetParameter(0)-dm));
  hr->SetBinContent(3,100.*(fr->GetParameter(1)-dm));
  hr->SetBinContent(4,100.*(fr->GetParameter(2)-dm));
  hr->SetBinContent(5,100.*(fr->GetParameter(3)-dm));
  hr->SetBinError(1,100.*fr->GetParError(0));
  hr->SetBinError(2,100.*fr->GetParError(0));
  hr->SetBinError(3,100.*sqrt(pow(fr->GetParError(0),2)+
			      pow(fr->GetParError(1),2)));
  hr->SetBinError(4,100.*sqrt(pow(fr->GetParError(0),2)+
			      pow(fr->GetParError(2),2)));
  hr->SetBinError(5,100.*sqrt(pow(fr->GetParError(0),2)+
			      pow(fr->GetParError(3),2)));

  TCanvas *c3 = tdrCanvas("c3",h3,4,11,kSquare);

  tdrDraw(hr,"P",kFullCircle);

  TBox box1(0.5,-0.40,1.5,0.40);
  box1.SetFillColorAlpha(kGray,0.30);
  box1.Draw("SAME");
  TBox box2(2.5,-0.40,4.5,0.40);
  box2.SetFillColorAlpha(kGray,0.30);
  box2.Draw("SAME");
  TLine *l = new TLine();
  l->SetLineStyle(kDotted);
  l->DrawLine(-0.5,0.,4.5,0.);

  gPad->RedrawAxis();

  // Observations on fractions: B,D high, A,C low
  // Hypotheses for regions
  // A: higher c-tag mistag rate in data, L down
  // B: higher c-tag mistag, lower c-tag eff, L+H both up
  // C: lower c-tag eff, H down
  // D: higher others (FSR), higher g-tag mistag rate, L+H up, O plus/minus

  //TH1D *h4 = tdrHist("h2","Leading jet flavor",-5.5,6.5,
  //		     "Subleading jet flavor",-5.5,6.5);
  //TH2D *h4 = new TH2D("h2",";Leading jet flavor;Subleading jet flavor",
  //		      12,-5.5,6.5,12,-5.5,6.5);
  //TCanvas *c4 = tdrCanvas("c4",(TH1D*)h4,4,0,kSquare);

  TH2D *h2LHGO = (TH2D*)fm->Get("h2LHGO"); assert(h2LHGO);
  TH2D *h2L    = (TH2D*)fm->Get("h2L");    assert(h2L);
  TH2D *h2H    = (TH2D*)fm->Get("h2H");    assert(h2H);
  TH2D *h2G    = (TH2D*)fm->Get("h2G");    assert(h2G);
  TH2D *h2O    = (TH2D*)fm->Get("h2O");    assert(h2O);

  h2LHGO->UseCurrentStyle();
  TCanvas *c4 = tdrCanvas("c4",(TH1D*)h2LHGO,4,0,kSquare);
  h2LHGO->SetXTitle("Leading jet flavor");
  h2LHGO->SetYTitle("Subleading jet flavor");

  h2LHGO->GetXaxis()->SetBinLabel(1,"#bar{b}");
  h2LHGO->GetXaxis()->SetBinLabel(2,"#bar{c}");
  h2LHGO->GetXaxis()->SetBinLabel(3,"#bar{s}");
  h2LHGO->GetXaxis()->SetBinLabel(4,"#bar{u}");
  h2LHGO->GetXaxis()->SetBinLabel(5,"#bar{d}");
  h2LHGO->GetXaxis()->SetBinLabel(7,"d");
  h2LHGO->GetXaxis()->SetBinLabel(8,"u");
  h2LHGO->GetXaxis()->SetBinLabel(9,"s");
  h2LHGO->GetXaxis()->SetBinLabel(10,"c");
  h2LHGO->GetXaxis()->SetBinLabel(11,"b");
  h2LHGO->GetXaxis()->SetBinLabel(12,"g");
  //
  h2LHGO->GetYaxis()->SetBinLabel(1,"#bar{b}");
  h2LHGO->GetYaxis()->SetBinLabel(2,"#bar{c}");
  h2LHGO->GetYaxis()->SetBinLabel(3,"#bar{s}");
  h2LHGO->GetYaxis()->SetBinLabel(4,"#bar{u}");
  h2LHGO->GetYaxis()->SetBinLabel(5,"#bar{d}");
  h2LHGO->GetYaxis()->SetBinLabel(7,"d");
  h2LHGO->GetYaxis()->SetBinLabel(8,"u");
  h2LHGO->GetYaxis()->SetBinLabel(9,"s");
  h2LHGO->GetYaxis()->SetBinLabel(10,"c");
  h2LHGO->GetYaxis()->SetBinLabel(11,"b");
  h2LHGO->GetYaxis()->SetBinLabel(12,"g");

  h2LHGO->SetFillStyle(1001); h2LHGO->SetFillColor(kBlack);
  h2LHGO->Draw("BOX SAME");
  h2L->SetFillStyle(1001); h2L->SetFillColor(kBlue);
  h2L->Draw("BOX SAME");
  h2H->SetFillStyle(1001); h2H->SetFillColor(kGreen+2);
  h2H->Draw("BOX SAME");
  h2G->SetFillStyle(1001); h2G->SetFillColor(kOrange+1);
  h2G->Draw("BOX SAME");
  h2O->SetFillStyle(1001); h2O->SetFillColor(kRed);
  h2O->Draw("BOX SAME");

  
  TH1D *h5 = tdrHist("h5","Momentum-based response",0.65,1.35,
		     "Flavor category",-0.5,2.5);
  TCanvas *c5 = tdrCanvas("c5",h5,4,11,kSquare);

  double ptref = 66.;
  //double ptrefm = 66.5;
  //double ptrefd = 65.;

  h5->GetXaxis()->SetBinLabel(1,"Up-type");
  h5->GetXaxis()->SetBinLabel(2,"Down-type");
  h5->GetXaxis()->SetBinLabel(3,"Both");

  TProfile *plnum = (TProfile*)fm->Get("plnu"); assert(plnum);
  TH1D *hlnum = plnum->ProjectionX("hlnum");
  TProfile *plnud = (TProfile*)fd->Get("plnu"); assert(plnud);
  TH1D *hlnud = plnud->ProjectionX("hlnud");

  // Scale lepton pT to unity (though slightly different phase space)
  double ptrefm = hlnum->GetBinContent(1);
  double ptrefd = hlnud->GetBinContent(1);

  hlnum->Scale(1./ptrefm);
  hlnud->Scale(1./ptrefd);

  TProfile *pudm = (TProfile*)fm->Get("pud"); assert(pudm);
  TH1D *hudm = pudm->ProjectionX("hudm");
  hudm->Scale(1./ptrefm);
  TProfile *pudd = (TProfile*)fd->Get("pud"); assert(pudd);
  TH1D *hudd = pudd->ProjectionX("hudd");
  hudd->Scale(1./ptrefm);

  TProfile *pcsm = (TProfile*)fm->Get("pcs"); assert(pcsm);
  TH1D *hcsm = pcsm->ProjectionX("hcsm");
  hcsm->Scale(1./ptrefm);
  TProfile *pcsd = (TProfile*)fd->Get("pcs"); assert(pcsd);
  TH1D *hcsd = pcsd->ProjectionX("hcsd");
  hcsd->Scale(1./ptrefd);

  tdrDraw(hudm,"HPz",kOpenSquare,kBlue,kSolid,-1,kNone);
  tdrDraw(hcsm,"HPz",kOpenSquare,kGreen+2,kSolid,-1,kNone);
  tdrDraw(hlnum,"HPz",kOpenSquare,kBlack,kSolid,-1,kNone);
  tdrDraw(hlnum,"HPz",kOpenSquare,kBlack,kSolid,-1,kNone);
  tdrDraw(hudd,"HPz",kFullCircle,kBlue,kSolid,-1,kNone);
  tdrDraw(hcsd,"HPz",kFullCircle,kGreen+2,kSolid,-1,kNone);
  tdrDraw(hlnud,"HPz",kFullCircle,kBlack,kSolid,-1,kNone);

  gPad->RedrawAxis();

  c1->SaveAs("pdf/hadW_cs_fractions.pdf");
  c2->SaveAs("pdf/hadW_cs_responses.pdf");
  c3->SaveAs("pdf/hadW_cs_flavorjes.pdf");
  c4->SaveAs("pdf/hadW_cs_2Dflavor.pdf");
  c5->SaveAs("pdf/hadW_cs_pairbalance.pdf");
} // hadW_cs


// Early testing code
void test_hadW_cs() {

  TChain *c = new TChain("tree");
  c->AddFile("rootfiles/HadW/UL18V5/WMass_Muo18_PowHeg.root");

  // No easy way to select first permutation?
  TCut cb = "abs(eta1)<1.3 && abs(eta2)<1.3 && fitProb>0.2"; // barrel
  //TCut cb = "abs(eta1)<1.3 && abs(eta2)<1.3 && fitProb>0.01"; // barrel
  TCut cm = "recoWMass>60 && recoWMass<110"; // mass

  // tight ctag
  // https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL18
  /*
  TCut ct1 = "ctag1/(ctag1+udstag1+gtag1)>0.282 &&"
    "ctag1/(ctag1+btag1+bleptag1)>0.267"; // ctag1
  TCut ct2 = "ctag2/(ctag2+udstag2+gtag2)>0.282 &&"
    "ctag2/(ctag2+btag2+bleptag2)>0.267"; //ctag2
  TCut ct = (ct1||ct2); // ctag
  TCut cu1 = "ctag1/(ctag1+udstag1+gtag1)<0.282 ||"
    "ctag1/(ctag1+btag1+bleptag1)<0.267"; // !ctag1
  TCut cu2 = "ctag2/(ctag2+udstag2+gtag2)<0.282 ||"
    "ctag2/(ctag2+btag2+bleptag2)<0.267"; // !ctag2
  TCut cu = (cu1&&cu2); // utag
  */


  // medium ctag
  // https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL18
  TCut ct1 = "ctag1/(ctag1+udstag1+gtag1)>0.099 &&"
    "ctag1/(ctag1+btag1+bleptag1)>0.325"; // ctag1
  TCut ct2 = "ctag2/(ctag2+udstag2+gtag2)>0.099 &&"
    "ctag2/(ctag2+btag2+bleptag2)>0.325"; //ctag2
  TCut ct = (ct1||ct2); // ctag
  TCut cu1 = "ctag1/(ctag1+udstag1+gtag1)<0.099 ||"
    "ctag1/(ctag1+btag1+bleptag1)<0.325"; // !ctag1
  TCut cu2 = "ctag2/(ctag2+udstag2+gtag2)<0.099 ||"
    "ctag2/(ctag2+btag2+bleptag2)<0.325"; // !ctag2
  TCut cu = (cu1&&cu2); // utag

  /*
  // loose ctag
  // https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL18
  TCut ct1 = "ctag1/(ctag1+udstag1+gtag1)>0.038 &&"
    "ctag1/(ctag1+btag1+bleptag1)>0.246"; // ctag1
  TCut ct2 = "ctag2/(ctag2+udstag2+gtag2)>0.038 &&"
    "ctag2/(ctag2+btag2+bleptag2)>0.246"; //ctag2
  TCut ct = (ct1||ct2); // ctag
  TCut cu1 = "ctag1/(ctag1+udstag1+gtag1)<0.038 ||"
    "ctag1/(ctag1+btag1+bleptag1)<0.246"; // !ctag1
  TCut cu2 = "ctag2/(ctag2+udstag2+gtag2)<0.038 ||"
    "ctag2/(ctag2+btag2+bleptag2)<0.246"; // !ctag2
  TCut cu = (cu1&&cu2); // utag
  */

  // Gluon antitag
  TCut ng = "gtag1<0.8 && gtag2<0.8";
  cb = (cb&&ng);

  // truth tags
  TCut cs = "(abs(flav1)==4&&abs(flav2)==3||abs(flav1)==3&&abs(flav2)==4)";
  TCut cd = "(abs(flav1)==4&&abs(flav2)==1||abs(flav1)==1&&abs(flav2)==4)";
  TCut su = "(abs(flav1)==3&&abs(flav2)==2||abs(flav1)==2&&abs(flav2)==3)";
  TCut ud = "(abs(flav1)==2&&abs(flav2)==1||abs(flav1)==1&&abs(flav2)==2)";

  /*
  c->Draw("abs(flav2):abs(flav1)>>h2cs(22,-0.5,21.5,22,-0.5,21.5)",cb&&cm&&ct,"box");
  c->Draw("abs(flav2):abs(flav1)>>h2ud(22,-0.5,21.5,22,-0.5,21.5)",cb&&cm&&cu,"box");
	
  TH2D *h2cs = (TH2D*)gROOT->FindObject("h2cs"); assert(h2ds);
  h2cs->SetFillStyle(1001);
  h2cs->SetFillColor(kGreen+2);
  TH2D *h2ud = (TH2D*)gROOT->FindObject("h2ud"); assert(h2ud);
  h2ud->SetFillStyle(1001);
  h2ud->SetFillColor(kBlue);

  h2ud->Draw("BOX");
  h2cs->Draw("BOX SAME");
  */

  c->Draw("recoWMass>>h1cs(50,60,110)",cb&&cm&&ct);
  c->Draw("recoWMass>>h1cs_cs(50,60,110)",cb&&cm&&ct&&cs);
  c->Draw("recoWMass>>h1cs_cd(50,60,110)",cb&&cm&&ct&&cd);
  c->Draw("recoWMass>>h1cs_su(50,60,110)",cb&&cm&&ct&&su);
  c->Draw("recoWMass>>h1cs_ud(50,60,110)",cb&&cm&&ct&&ud);
  c->Draw("recoWMass>>h1cs_ot(50,60,110)",cb&&cm&&ct&&!(cs||ud||cd||su));
  c->Draw("recoWMass>>h1ud(50,60,110)",cb&&cm&&cu);
  c->Draw("recoWMass>>h1ud_cs(50,60,110)",cb&&cm&&cu&&cs);
  c->Draw("recoWMass>>h1ud_cd(50,60,110)",cb&&cm&&cu&&cd);
  c->Draw("recoWMass>>h1ud_su(50,60,110)",cb&&cm&&cu&&su);
  c->Draw("recoWMass>>h1ud_ud(50,60,110)",cb&&cm&&cu&&ud);
  c->Draw("recoWMass>>h1ud_ot(50,60,110)",cb&&cm&&cu&&!(cs||ud||cd||su));

  TH2D *h1ud = (TH2D*)gROOT->FindObject("h1ud"); assert(h1ud);
  h1ud->SetFillStyle(1001);
  h1ud->SetFillColor(kBlue-9);
  TH2D *h1cs = (TH2D*)gROOT->FindObject("h1cs"); assert(h1cs);
  //h1cs->SetFillStyle(1001);
  //h1cs->SetFillColor(kGreen+2);
  h1cs->SetLineColor(kGreen+2);

  TH2D *h1ud_cs = (TH2D*)gROOT->FindObject("h1ud_cs"); assert(h1ud_cs);
  TH2D *h1ud_cd = (TH2D*)gROOT->FindObject("h1ud_cd"); assert(h1ud_cd);
  TH2D *h1ud_su = (TH2D*)gROOT->FindObject("h1ud_su"); assert(h1ud_su);
  TH2D *h1ud_ud = (TH2D*)gROOT->FindObject("h1ud_ud"); assert(h1ud_ud);
  TH2D *h1ud_ot = (TH2D*)gROOT->FindObject("h1ud_ot"); assert(h1ud_ot);
  
  TH2D *h1cs_cs = (TH2D*)gROOT->FindObject("h1cs_cs"); assert(h1cs_cs);
  TH2D *h1cs_cd = (TH2D*)gROOT->FindObject("h1cs_cd"); assert(h1cs_cd);
  TH2D *h1cs_su = (TH2D*)gROOT->FindObject("h1cs_su"); assert(h1cs_su);
  TH2D *h1cs_ud = (TH2D*)gROOT->FindObject("h1cs_ud"); assert(h1cs_ud);
  TH2D *h1cs_ot = (TH2D*)gROOT->FindObject("h1cs_ot"); assert(h1cs_ot);
  
  h1cs_cs->SetLineColor(kGreen+3);
  h1cs_cd->SetLineColor(kGreen+4);
  h1cs_su->SetLineColor(kCyan+3);
  h1cs_ud->SetLineColor(kCyan+2);
  h1cs_ot->SetLineColor(kGray+1);
  
  h1ud_cs->SetLineColor(kMagenta+2);
  h1ud_cd->SetLineColor(kMagenta+3);
  h1ud_su->SetLineColor(kBlue+3);
  h1ud_ud->SetLineColor(kBlue+2);
  h1ud_ot->SetLineColor(kOrange+2);
  
  h1ud_cs->Scale(1./h1ud->Integral());
  h1ud_cd->Scale(1./h1ud->Integral());
  h1ud_su->Scale(1./h1ud->Integral());
  h1ud_ud->Scale(1./h1ud->Integral());
  h1ud_ot->Scale(1./h1ud->Integral());
  
  h1cs_cs->Scale(1./h1cs->Integral());
  h1cs_cd->Scale(1./h1cs->Integral());
  h1cs_su->Scale(1./h1cs->Integral());
  h1cs_ud->Scale(1./h1cs->Integral());
  h1cs_ot->Scale(1./h1cs->Integral());
  
  h1ud->Scale(1./h1ud->Integral());
  h1cs->Scale(1./h1cs->Integral());
  
  h1ud->Draw("H");
  h1cs->Draw("H SAMES");
  
  h1ud_cs->Draw("H SAMES");
  h1ud_cd->Draw("H SAMES");
  h1ud_su->Draw("H SAMES");
  h1ud_ud->Draw("H SAMES");
  h1ud_ot->Draw("H SAMES");
  

  h1cs_cs->Draw("H SAMES");
  h1cs_cd->Draw("H SAMES");
  h1cs_su->Draw("H SAMES");
  h1cs_ud->Draw("H SAMES");
  h1cs_ot->Draw("H SAMES");
} // test_hadW_cs


// Fit to mass-based response in ABCD regions
// Allow others response to differ in ABCD vs D
// by estimating different gluon fraction
Double_t frABCD(Double_t *x, Double_t *p) {

  assert(h2mm); // masses in MC
  assert(h2nf); // tag fractions (fitted) in data
  assert(h2mf); // masses fit to data

  int k = int(*x)+1;
  assert(k>=1);
  assert(k<=5);

  // Create full response matrix for easier plotting later
  // although this is a bit of extra work
  for (int i = 1; i != h2nf->GetNbinsX()+1; ++i) {

    double sumw(0), sumrw(0);
    for (int j = 2; j != h2nf->GetNbinsY()+1; ++j) {

      double w = h2nf->GetBinContent(i,j);
      double r0 = h2mm->GetBinContent(i,j);
      double k = p[j-2];
      double r = k * r0;

      sumrw += r * w;
      sumw += w;
      h2mf->SetBinContent(i, j, r);
    } // for j>1

    h2mf->SetBinContent(i, 1, sumrw / sumw);
  } // for i

  return (h2mf->GetBinContent(k,1));
}

// Fit to event numbers in ABCD
Double_t fnABCD(Double_t *x, Double_t *p) {

  int k = int(*x);

  if (k<=0) {
    assert(-k<4);
    return p[-k];
  }
  // int k = nx*(j-1) + (i-1) + 1;
  int nx = h2nd->GetNbinsX();
  int i0 = ((k-1) % nx) + 1;
  int j0 = ((k-1) / nx) + 1;

  assert(h2nm);
  assert(h2nd);
  assert(h2nf);

  // First redistribute MC events to ABCD for each of LHGO
  for (int i = 2; i != h2nm->GetNbinsX()+1; ++i) {
    for (int j = 2; j != h2nm->GetNbinsY()+1; ++j) {

      bool isA = (i==2); int iA = 2;
      bool isB = (i==3); int iB = 3;
      bool isC = (i==4); int iC = 4;
      bool isD = (i==5); int iD = 5;

      bool isL = (j==2); int jL = 2;
      bool isH = (j==3); int jH = 3;
      bool isG = (j==4); int jG = 4;
      bool isO = (j==5); int jO = 5;

      double n(1);
      
      // Light mistag rate (eff in A vs BC, D vs ABCD)
      if (isL) {
	double nL = h2nm->GetBinContent(1, jL);
	double pA = h2nm->GetBinContent(iA,jL) / nL;
	double pB = h2nm->GetBinContent(iB,jL) / nL;
	double pC = h2nm->GetBinContent(iC,jL) / nL;
	double pD = h2nm->GetBinContent(iD,jL) / nL;
	double pX = h2nm->GetBinContent(i, jL) / nL;
	double kLA = p[0]; // A relative decrease (charm mistag)
	double kLD = p[2]; // D relative increase (gluon mistag)
	double mLA = kLA*pA/(pB+pC); // BC relative increase
	double mLD = kLD*pD/(pA+pB+pC); // ABC relative decrease
	if (isA)      n = (1-kLA)*(1-mLD)*pX*nL;
	if (isB||isC) n = (1+kLA)*(1-mLD)*pX*nL;
	if (isD)      n = (1+kLD)*pX*nL;
      }

      // Charm efficiency (eff in C vs AB)
      if (isH) {
	double nH = h2nm->GetBinContent(1, jH);
	double pA = h2nm->GetBinContent(iA,jH) / nH;
	double pB = h2nm->GetBinContent(iB,jH) / nH;
	double pC = h2nm->GetBinContent(iC,jH) / nH;
	double pD = h2nm->GetBinContent(iD,jH) / nH;
	double pX = h2nm->GetBinContent(i, jH) / nH;
	double kHC = p[1]; // C relative decrease (charm efficiency)
	double kHD = p[2]; // D relative increase (gluon mistag)
	double mHC = kHC*pC/(pA+pB); // AB relative increase
	double mHD = kHD*pD/(pA+pB+pC); // ABC relative increase
	if (isC)      n = (1-kHC)*(1-mHD)*pX*nH;
	if (isA||isB) n = (1+mHC)*(1-mHD)*pX*nH;
	if (isD)      n = (1+kHD)*pX*nH;
      }
      
      // gluon jet rate
      if (isG) {
	double nG = h2nm->GetBinContent(1, jG);
	double pX = h2nm->GetBinContent(i, jG) / nG;
	double kG = p[3]; // gluon jet rate
	n = (1+kG)*pX*nG;
      }
      
      // other jet rate
      if (isO) {
	double nO = h2nm->GetBinContent(1, jO);
	double pX = h2nm->GetBinContent(i, jO) / nO;
	double kO = 0;//p[3];
	n = (1+kO)*pX*nO;
      }

      h2nf->SetBinContent(i,j,n);
    } // for j
  } // for i


  double nd = h2nd->GetBinContent(1,1);
  double nf = h2nf->Integral(2,5,2,5);
  double nm = h2nm->GetBinContent(1,1);
  if (!((nf/nm-1)<1e-3)) {
    cout << "nf="<<nf<<" nm="<<nm<<" nd="<<nd<<endl<<flush;
    assert((nf/nm-1)<1e-3);
  }
  //double kd = nd/nm;
  double kd = nd/nf;

  // Normalize exclusive L,H,G,O x A,B,C,D to data integral
  for (int i = 2; i != h2nf->GetNbinsX()+1; ++i) {
    for (int j = 2; j != h2nf->GetNbinsY()+1; ++j) {
      h2nf->SetBinContent(i,j,kd*h2nf->GetBinContent(i,j));
    } // for j
  } // for i

  // Sum up LHGO for exclusive A,B,C,D
  for (int i = 2; i != h2nf->GetNbinsX()+1; ++i) {
    h2nf->SetBinContent(i,1,h2nf->Integral(i,i,2,5));
  } // for i

  // Sum up ABCD for exclusive L,H,G,O
  for (int j = 2; j != h2nf->GetNbinsY()+1; ++j) {
    h2nf->SetBinContent(1,j,h2nf->Integral(2,5,j,j));
  } // for j

  // Set total to data
  //h2f->SetBinContent(1,1,nd);
  h2nf->SetBinContent(1,1,h2nf->Integral(2,5,2,5));
  
  return (h2nf->GetBinContent(i0,j0));
} // fnABCD
