// Purpose: draw data-MC difference in PF composition
//          and try to use it to estimate JEC shift
//          and/or source (tracker, ECAL, HCAL) of likely JEC shift

// Other tools: jecsys2018/drawHerwigVsPythia.C
//              jecsys2018/drawProbePerTag.C
//
// References to detector performance:
// * Tracking *
// https://twiki.cern.ch/twiki/bin/view/CMSPublic/PixelOfflinePlotsOctober2018#Bad_Components (one sector in Layer 3 and two sectors in Layer 4 off; overlaps with dead module in Layer 1? (module 2, ladder 3))
// https://twiki.cern.ch/twiki/bin/view/CMSPublic/PixelOfflinePlotsAugust2017#Dynamic_Inefficiencies (Layer 1 efficiency drops from 99.0-99.5% at PU<30 (RunBCDE=) to 97.5% at PU=50 (RunF))
// * HCAL (IsoTrack), from Marina Chadeeva through Markus Seidel *
// https://indico.cern.ch/event/804390/contributions/3346405/attachments/1820286/2977230/IsoTrackN84.pdf
// https://indico.cern.ch/event/804390/contributions/3346405/attachments/1820286/2977086/chadeeva_29mar2019.pdf
#include "TFile.h"
#include "TProfile.h"
#include "TF1.h"

#include <vector>
#include <map>

#include "../tdrstyle_mod15.C"

using namespace std;

const  bool doscale = false;//true;
const  bool dodelta = false;//true;
//bool doscale = true;

void drawPFcompMC(string era="E", string tp="");
void drawPFcompIOV(string sf, string seta, string tp="tp");

void drawPFcomp() {
  /*
  drawPFcompMC("B","");
  drawPFcompMC("C","");
  drawPFcompMC("D","");
  drawPFcompMC("E","");
  drawPFcompMC("F","");
  drawPFcompMC("BCDEF","");
  drawPFcompMC("B","tp");
  drawPFcompMC("C","tp");
  drawPFcompMC("D","tp");
  drawPFcompMC("E","tp");
  drawPFcompMC("F","tp");
  drawPFcompMC("BCDEF","tp");
  */
  drawPFcompIOV("chf","barrel","tp");
  drawPFcompIOV("nef","barrel","tp");
  drawPFcompIOV("nhf","barrel","tp");
  drawPFcompIOV("puf","barrel","tp");
  drawPFcompIOV("chf","barrel","");
  drawPFcompIOV("nef","barrel","");
  drawPFcompIOV("nhf","barrel","");
  drawPFcompIOV("puf","barrel","");

  /*
  drawPFcompIOV("chf","endcap");
  drawPFcompIOV("nef","endcap");
  drawPFcompIOV("nhf","endcap");
  */
}

void drawPFcompMC(string era, string tp) {

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  const char *cera = era.c_str();
  const char *ctp = tp.c_str();

  //TFile *fd = new TFile("rootfiles/output-DATA-2b-UL17V2_BCDEF.root","READ");
  //TFile *fd = new TFile("rootfiles/output-DATA-2b-UL17V4_BCDEF.root","READ");
  //TFile *fd = new TFile("rootfiles/output-DATA-2b-UL17V4_E.root","READ");
  //TFile *fd = new TFile("rootfiles/output-DATA-2b-UL17V4_BCDEF.root","READ");
  TFile *fd = new TFile(Form("rootfiles/output-DATA-2b-UL17V4_%s.root",cera),
			"READ");
  assert(fd && !fd->IsZombie());
  //TFile *fm = new TFile("rootfiles/output-MCNU-2b-UL17V2_BCDEF.root","READ");
  //TFile *fm = new TFile("rootfiles/output-MCNU-2b-UL17V4_BCDEF.root","READ");
  //TFile *fm = new TFile("rootfiles/output-MC80NU-2b-UL17V4_E.root","READ");
  //TFile *fm = new TFile("rootfiles/output-MC80NU-2b-UL17V4_BCDEF.root","READ");
  TFile *fm = new TFile(Form("rootfiles/output-MC80NU-2b-UL17V4_%s.root",cera),
			"READ");
  assert(fm && !fm->IsZombie());

  vector<string> vf;
  vf.push_back("chf");
  vf.push_back("nhf");
  vf.push_back("nef");
  vf.push_back("cef");
  vf.push_back("muf");
  vf.push_back("puf");

  map<string, int> color;
  color["chf"] = kRed;
  color["puf"] = kRed+2;
  color["nhf"] = kGreen+2;
  color["nef"] = kBlue;
  color["muf"] = kMagenta+2;
  color["cef"] = kCyan+2;
  map<string, int> markerd;
  markerd["chf"] = kFullCircle;
  markerd["puf"] = kOpenCircle;
  markerd["nhf"] = kFullDiamond;
  markerd["nef"] = kFullSquare;
  markerd["muf"] = kDot;
  markerd["cef"] = kDot;
  map<string, int> markerm;
  markerm["chf"] = kOpenCircle;
  markerm["puf"] = kOpenDiamond;
  markerm["nhf"] = kOpenDiamond;
  markerm["nef"] = kOpenSquare;
  markerm["muf"] = kDot;
  markerm["cef"] = kDot;
  map<string, string> legendm;
  legendm["chf"] = "Charged hadrons";
  legendm["nhf"] = "Neutral hadrons";
  legendm["nef"] = "Photons";
  legendm["muf"] = "Muons";
  legendm["cef"] = "Electrons";
  legendm["puf"] = "Charged PU";

  /*
  map<string, double> scale;
  scale["chf"] = 1-1.0/62.-0.5/62.;
  scale["nhf"] = 1;
  scale["nef"] = 1-0.2/25.-0.4/25.-0.3/62.;
  scale["muf"] = 1;
  scale["cef"] = 1;
  */

  // Scaling component fractions
  map<string, TF1*> fscale;
  //fscale["chf"] = new TF1("chf","1-1.0/62.-0.5/62.",15,3500);
  //fscale["chf"] = new TF1("chf","1-1.3/62.",15,3500);
  //fscale["chf"] = new TF1("chf","1-1./62.",15,3500);
  fscale["chf"] = new TF1("chf","1-0.015",15,3500);
  //fscale["chf"] = new TF1("chf","1",15,3500);
  //fscale["nhf"] = new TF1("nhf","1+([0]+[1]*log(x/208.))/10.",15,3500);
  //fscale["nhf"]->SetParameters(0.017408025,0.25303863);
  fscale["nhf"] = new TF1("nhf","1",15,3500);
  //fscale["nef"] = new TF1("chf","1-0.2/25.-0.4/25.-0.15/25.",15,3500);
  //fscale["nef"] = new TF1("nef","1-0.65/25.",15,3500);
  fscale["nef"] = new TF1("nef","1",15,3500);
  //fscale["muf"] = new TF1("muf","1.2",15,3500);
  fscale["muf"] = new TF1("muf","1",15,3500);
  fscale["cef"] = new TF1("cef","1",15,3500);

  // Applying delta shift to component fractions
  map<string, TF1*> fdelta;
  fdelta["chf"] = new TF1("dchf","0",15,3500);
  //fdelta["nhf"] = new TF1("dnhf","0",15,3500);
  //fdelta["nef"] = new TF1("dnef","0",15,3500);
  fdelta["nef"] = new TF1("dnef","-0.0050",15,3500);
  fdelta["muf"] = new TF1("dmuf","+0.0005",15,3500);
  fdelta["cef"] = new TF1("dcef","+0.0005",15,3500);

  // Fits from minitools/varPlots.C
  TF1 *fx = new TF1("fc","0.01*1.0862*([p0]+[p1]*pow(x/[p2],[p3])/(1+pow(x/[p2],[p3]))*(1-pow(x/[p2],-[p3])))",15,4500);
  fx->SetParameters(1.184, 1.428, 1402, 1.225); // toyPF
  fdelta["nhf"] = fx;

 
  //map<string, TProfile*> vd;
  //map<string, TProfile*> vm;
  map<string, TH1D*> vd;
  map<string, TH1D*> vm;
  map<string, TH1D*> vm2;
  map<string, TH1D*> vr;
  map<string, TH1D*> vs;

  TH1D *hmsum(0);
  const char *cd = "Standard/Eta_0.0-1.3/";
  for (unsigned int i = 0; i != vf.size(); ++i) {

    string sf = vf[i];
    const char *cf = sf.c_str();
    TProfile *pd = (TProfile*)fd->Get(Form("%sp%s%s",cd,cf,ctp)); assert(pd);
    TProfile *pm = (TProfile*)fm->Get(Form("%sp%s%s",cd,cf,ctp)); assert(pm);

    TH1D *hd = pd->ProjectionX(Form("hd_%s",cf));
    TH1D *hm = pm->ProjectionX(Form("hm_%s",cf));
    TH1D *hm2 = pm->ProjectionX(Form("hm2_%s",cf));

    // Difference of unscaled fractions
    TH1D *hr = pd->ProjectionX(Form("hr_%s",cf));
    hr->Add(hd,hm,+100,-100);

    // Scale MC components
    //pm->Scale(scale[sf]);
    if (doscale) hm2->Multiply(fscale[sf]);
    if (dodelta) hm2->Add(fdelta[sf]);
      /*
   {
      for (int i = 1; i != pm->GetNbinsX()+1; ++i) {
	double pt = pm->GetBinCenter(i);
	double d = fdelta[sf]->Eval(pt);
	pm->SetBinContent(i, pm->GetBinContent(i)+d)
      }
    }
      */

    //if (!hmsum) hmsum = pm->ProjectionX("hmsum");
    if (!hmsum) hmsum = (TH1D*)hm2->Clone("hmsum");
    else if (sf!="puf") hmsum->Add(hm2);//pm);

    vd[sf] = hd;
    vm[sf] = hm;
    vm2[sf] = hm2;
    vr[sf] = hr;
  }

  // Renormalize MC back to unity, calculate new fractions
  for (unsigned int i = 0; i != vf.size(); ++i) {

    string sf = vf[i];
    const char *cf = sf.c_str();
    TH1D *hm2 = vm2[sf];//->ProjectionX(Form("hm_%s",sf.c_str()));
    hm2->Divide(hmsum);

    TH1D *hs = (TH1D*)hm2->Clone(Form("hs_%s",cf));
    hs->Add(vd[sf],hm2,+100,-100);

    vm2[sf] = hm2;
    vs[sf] = hs;
  }

  // Sum everything up to make sure it's at unity
  TH1D *hr = (TH1D*)vr["chf"]->Clone("hr");
  hr->Reset();
  for (unsigned int i = 0; i != vf.size(); ++i) {
    string sf = vf[i];
    //hr->Add(vr[sf]);
    if (sf!="puf") hr->Add(vs[sf]);
  }

  curdir->cd();

  TH1D *hu = new TH1D("hu",";p_{T,tag} (GeV);Probe E fraction",3485,15,3500);
  if (tp=="") hu->SetXTitle("p_{T,probe}");
  hu->GetXaxis()->SetMoreLogLabels();
  hu->GetXaxis()->SetNoExponent();
  hu->SetMaximum(0.75);
  hu->SetMinimum(0.00);

  TH1D *hd = new TH1D("hd",";p_{T,tag} (GeV);Data-MC (10^{-2})",3485,15,3500);
  if (tp=="") hd->SetXTitle("p_{T,probe}");
  hd->GetXaxis()->SetMoreLogLabels();
  hd->GetXaxis()->SetNoExponent();
  hd->SetMaximum(+1.2);
  hd->SetMinimum(-1.2);

  //lumi_13TeV = "UL17 BCDEF";
  string lum = Form("UL17 %s",cera);
  lumi_13TeV = lum.c_str();//"UL17 E";
  TCanvas *c1 = tdrDiCanvas("c1",hu,hd,4,11);

  c1->cd(1);
  gPad->SetLogx();

  TLegend *leg1 = tdrLeg(0.20,0.38,0.40,0.68);
  leg1->SetHeader("MC");
  TLegend *leg2 = tdrLeg(0.27,0.38,0.47,0.68);
  leg2->SetHeader("Data");

  for (unsigned int i = 0; i != vf.size(); ++i) {
    string sf = vf[i];
    //tdrDraw(vm[sf],"Pz",markerm[sf],color[sf]);
    tdrDraw(vm2[sf],"Pz",markerm[sf],color[sf]);
    tdrDraw(vd[sf],"Pz",markerd[sf],color[sf]);

    leg1->AddEntry(vm2[sf]," ","PLE");
    leg2->AddEntry(vd[sf],legendm[sf].c_str(),"PLE");
  }

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);

  if (!(doscale || dodelta))
    tex->DrawLatex(0.50,0.84,"Uncorrected composition");
  else
    tex->DrawLatex(0.50,0.84,"Corrected composition");

  c1->cd(2);
  gPad->SetLogx();

  tdrDraw(hr,"HIST",kBlack);
  for (unsigned int i = 0; i != vf.size(); ++i) {
    string sf = vf[i];
    //tdrDraw(vr[sf],"Pz",markerm[sf],color[sf]); vr[sf]->SetMarkerSize(0.7);
    tdrDraw(vs[sf],"Pz",markerd[sf],color[sf]); vs[sf]->SetMarkerSize(0.7);
  }

  tex->SetTextSize(0.090);
  if (!(doscale || dodelta)) {

    tex->DrawLatex(0.35,0.80,"Uncorrected composition");
    c1->SaveAs("pdf/drawPFComp_uncorr.pdf");
  }
  else {
    tex->DrawLatex(0.35,0.80,"Corrected composition");
    c1->SaveAs("pdf/drawPFComp_corr.pdf");
  }

  c1->SaveAs(Form("pdf/drawPFCompMC%s_%s.pdf",ctp,cera));
} // drawPFcompMC

// Evolution in fractions over time
void drawPFcompIOV(string sf, string seta, string tp) {
  const char *cf = sf.c_str();
  const char *ceta = seta.c_str();
  const char *ctp = tp.c_str();
    
  setTDRStyle();
  TDirectory *curdir = gDirectory;
  
  string aIOV[] = {"BCDEF","B","C","D","E","F"};
  int nIOV = sizeof(aIOV)/sizeof(aIOV[0]);

  map<string, int> markerd;
  markerd["BCDEF"] = kFullSquare;
  /*
  markerd["B"] = kOpenCircle;
  markerd["C"] = kOpenCircle;
  markerd["D"] = kOpenSquare;
  markerd["E"] = kOpenSquare;
  markerd["F"] = kOpenDiamond;
  */
  markerd["B"] = kFullCircle;
  markerd["C"] = kFullDiamond;
  markerd["D"] = kFullSquare;
  markerd["E"] = kFullCross;
  markerd["F"] = kFullCrossX;

  map<string, int> colord;
  colord["BCDEF"] = kBlack;
  /*
  colord["B"] = kRed;
  colord["C"] = kRed+1;
  colord["D"] = kGreen+2;
  colord["E"] = kGreen+3;
  colord["F"] = kBlue;
  */

  colord["B"] = kBlue;
  colord["C"] = kRed;
  colord["D"] = kGreen+2;
  colord["E"] = kMagenta+1;
  colord["F"] = kOrange+2;

  const double ptmin = 15;
  const double ptmax = 1700;
  TH1D *hu = new TH1D("huI",Form(";p_{T,tag} (GeV);Probe %s E fraction in data",
				 cf),
		      int(ptmax-ptmin),ptmin,ptmax);
  hu->GetXaxis()->SetMoreLogLabels();
  hu->GetXaxis()->SetNoExponent();
  hu->SetMaximum(1.00);//0.75);
  hu->SetMinimum(0.00);

  TH1D *hd = new TH1D("hdI",";p_{T,tag} (GeV);Data-MC (10^{-2})",
		      int(ptmax-ptmin),ptmin,ptmax);
  hd->GetXaxis()->SetMoreLogLabels();
  hd->GetXaxis()->SetNoExponent();
  if (tp=="") hd->SetXTitle("p_{T,probe} (GeV)");
  hd->SetMaximum(seta=="barrel" ? +1.5 : +3.0);
  hd->SetMinimum(seta=="barrel" ? -1.5 : -3.0);

  if (sf=="nef") hu->SetMaximum(0.55);
  if (sf=="puf") hu->SetMaximum(0.40);//0.35);
  if (sf=="nhf") hu->SetMaximum(0.30);
  if (sf=="puf") hd->SetMinimum(-3.0);
  if (sf=="puf") hd->SetMaximum(+2.0);

  TH1D *h = new TH1D("hI",Form(";p_{T,tag} (GeV);IOV-BCDEF for %s (10^{-2})",
			       cf),
		     int(ptmax-ptmin),ptmin,ptmax);
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  if (tp=="") h->SetXTitle("p_{T,probe} (GeV)");
  h->SetMaximum(seta=="barrel" ? +1.2 : +2.4);
  h->SetMinimum(seta=="barrel" ? -1.2 : -2.4);

  lumi_13TeV = "UL17 BCDEF";
  TCanvas *c2 = tdrDiCanvas(Form("c2_%s_iov",cf),hu,hd,4,11);
  TCanvas *c1 = tdrCanvas(Form("c1_%s_iov",cf),h,4,11,kSquare);

  c2->cd(1);
  gPad->SetLogx();
  TLegend *leg2 = tdrLeg(0.70,0.62,0.90,0.87);
  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  if (seta=="barrel") tex->DrawLatex(0.39,0.84,"Barrel |#eta|<1.3");
  if (seta=="endcap") tex->DrawLatex(0.39,0.84,"Endcap 2<|#eta|<2.5");
  c2->cd(2);
  gPad->SetLogx();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(ptmin,0,ptmax,0);

  c1->cd();
  gPad->SetLogx();
  TLegend *leg1 = tdrLeg(0.70,0.70,0.90,0.90);
  if (seta=="barrel") tex->DrawLatex(0.38,0.87,"Barrel |#eta|<1.3");
  if (seta=="endcap") tex->DrawLatex(0.38,0.87,"Endcap 2<|#eta|<2.5");

  TH1D *hr0(0), *hd0(0);
  for (int iov = 0; iov != nIOV; ++iov) {

    string sera = aIOV[iov];
    const char *cera = sera.c_str();

    // For puf reference
    TF1 *fpu = new TF1(Form("fpu_%s",cera),"[0]/x",ptmin,ptmax);
    fpu->SetLineStyle(kDotted);
    fpu->SetParameter(0,tp=="tp" ? 5.5/0.9 : 5.5);

    TFile *fd = new TFile(Form("rootfiles/output-DATA-2b-UL17V4_%s.root",
			       cera), "READ");
    assert(fd && !fd->IsZombie());
    //TFile *fm = new TFile(Form("rootfiles/output-MCNU-2b-UL17V4_%s.root",
    TFile *fm = new TFile(Form("rootfiles/output-MC80NU-2b-UL17V4_%s.root",
			       cera), "READ");
    assert(fm && !fm->IsZombie());
    
    curdir->cd();

    //string sf = "chf";
    //const char *cf = sf.c_str();
    string sd = "";
    if (seta=="barrel") sd = "Standard/Eta_0.0-1.3/";
    if (seta=="endcap") sd = "Standard/Eta_2.0-2.5/";
    const char *cd = sd.c_str();
    //TProfile *pd = (TProfile*)fd->Get(Form("%sp%stp",cd,cf)); assert(pd);
    //TProfile *pm = (TProfile*)fm->Get(Form("%sp%stp",cd,cf)); assert(pm);
    //TProfile *pd = (TProfile*)fd->Get(Form("%sp%s",cd,cf)); assert(pd);
    //TProfile *pm = (TProfile*)fm->Get(Form("%sp%s",cd,cf)); assert(pm);
    TProfile *pd = (TProfile*)fd->Get(Form("%sp%s%s",cd,cf,ctp)); assert(pd);
    TProfile *pm = (TProfile*)fm->Get(Form("%sp%s%s",cd,cf,ctp)); assert(pm);

    c2->cd(1);
    fpu->SetLineColor(colord[sera]);
    fpu->SetParameter(0,pd->GetBinContent(pd->GetXaxis()->FindBin(20.))*20.);
    if (sf=="puf") fpu->Draw("SAME");

    tdrDraw(pd,"Pz",markerd[sera],colord[sera]);
    //pd->SetMarkerSize(0.5);
    leg2->AddEntry(pd,cera,"PLE");

    TH1D *hd = pd->ProjectionX(Form("hd_%s_%s",cf,cera));
    TH1D *hm = pm->ProjectionX(Form("hm_%s_%s",cf,cera));

    // Difference of unscaled fractions
    TH1D *hr = pd->ProjectionX(Form("hr_%s_%s",cf,cera));
    hr->Add(hd,hm,+100,-100);

    c2->cd(2);
    tdrDraw(hr,"Pz",markerd[sera],colord[sera]);
    //hr->SetMarkerSize(0.5);

    c1->cd();

    if (hr0==0) hr0 = hr;
    if (hd0==0) hd0 = hd;
    TH1D *hr2 = (TH1D*)hr->Clone(Form("hr2_%s_%s",cf,cera));
    hr2->Add(hr0,-1);

    // Fix statistical uncertainty: MC zero, data correlated
    for (int i = 1; i != hr2->GetNbinsX()+1; ++i) {
      double s0 = hd0->GetBinError(i);
      double s1 = hd->GetBinError(i);
      double s2 = sqrt(fabs(s1*s1-s0*s0));
      hr2->SetBinError(i, 100.*s2);
    }

    tdrDraw(hr2,"Pz",markerd[sera],colord[sera]);
    //hr2->SetMarkerSize(0.5);
    if (markerd[sera]==kFullDiamond) hr2->SetMarkerSize(1.3);
    if (markerd[sera]==kFullCross || markerd[sera]==kFullCrossX)
      hr2->SetMarkerSize(1.2);
    if (sera!="BCDEF") leg1->AddEntry(hr2,cera,"PLE");
    
    TF1 *f1 = new TF1(Form("f1_%s_%s",cf,cera),
		      //"[0]+[1]*(pow(x/208.,[2])-1)",
		      "[0]+[1]*log(x/208)+[2]/x",
		      ptmin,ptmax);
    f1->SetParameters(0,0,0);
    f1->FixParameter(2,0);
    hr2->Fit(f1,"QRN");
    f1->SetLineColor(colord[sera]);
    f1->SetLineStyle(kDotted);
    f1->Draw("SAME");

  } // for iov

  c2->SaveAs(Form("pdf/drawPFCompIOV2%s_%s_%s.pdf",ctp,ceta,cf));
  c1->SaveAs(Form("pdf/drawPFCompIOV1%s_%s_%s.pdf",ctp,ceta,cf));
} // drawPFcompIOV
