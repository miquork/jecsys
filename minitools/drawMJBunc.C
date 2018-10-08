#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TF1.h"

#include "../tdrstyle_mod15.C"

#include <vector>
#include <string>

using namespace std;

void drawMJBunc() {

  TDirectory *curdir = gDirectory;

  setTDRStyle();

  TH1D *h = new TH1D("h",";p_{T} (GeV);Uncertainty (%)",100,200,3500);
  h->SetMinimum(-1.0);
  h->SetMaximum(+1.5);
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();

  TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);

  TLegend *leg = tdrLeg(0.4,0.60,0.7,0.86);

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);
  tex->DrawLatex(0.3,0.17,"Full: MPF, Open: MJB");

  map<string, map<string, int> > marker;
  marker["MPF"]["JEC"] = kFullSquare;
  marker["MPF"]["JER"] = kFullCircle;
  marker["MPF"]["PU"] = kFullDiamond;
  marker["MPF"]["FSR"] = kOpenTriangleDown;//kFullTriangleDown;
  marker["MJB"]["JEC"] = kOpenSquare;
  marker["MJB"]["JER"] = kOpenCircle;
  marker["MJB"]["PU"] = kOpenDiamond;
  marker["MJB"]["FSR"] = kOpenTriangleUp;//kOpenTriangleDown;

  map<string, map<string, int> > color;
  color["MPF"]["JEC"] = kRed;
  color["MPF"]["JER"] = kBlue;
  color["MPF"]["PU"] = kGreen+2;
  color["MPF"]["FSR"] = kOrange+2;
  color["MJB"]["JEC"] = kRed;
  color["MJB"]["JER"] = kBlue;
  color["MJB"]["PU"] = kGreen+2;
  color["MJB"]["FSR"] = kOrange+2;

  // === CODE FROM MULTIJET.C === (start)

  //double ptbins[] = {200, 250, 300, 330, 370, 410, 450, 510, 530, 550, 575,
  //                 600, 650, 700, 800, 900, 1000, 1100, 1200, 1400, 1600,
  //                 1800, 2000, 2200};// 2016Early (hass 3000 also?)  
  double ptbins[] = {200, 250, 300, 330, 370, 420, 465, 505, 550,
                     620, 700, 800, 1000, 1200, 1500,
                     2000,3500}; 
  const int npt = sizeof(ptbins)/sizeof(ptbins[0])-1;
  
  // Load Anne-Laure's MJB systematics                                          
  TFile *fs = new TFile("../rootfiles/compareMC_signedValues_recoilPtHLTBin.root",
                        "READ");
  assert(fs && !fs->IsZombie());

  vector<TH1D*> vsrc;
  vector<string> src;
  src.push_back("SystTot");
  src.push_back("JEC");
  src.push_back("JER");
  src.push_back("PU");

  const int nmethods = 2;
  //const int nmethods = 1; // MJB only
  for (unsigned int isrc = 0; isrc != src.size(); ++isrc) {

    for (int im = 0; im != nmethods; ++im) {

      //const char *cm = (im==0 ? "MJB" : "MPF");
      const char *cm = (im==0 ? "MPF" : "MJB");

      assert(fs->cd(cm));
      TDirectory *ds1 = fs->GetDirectory(cm); assert(ds1);
      assert(ds1->cd("PtBin"));
      TDirectory *ds2 = ds1->GetDirectory("PtBin"); assert(ds2);
      const char *cs = src[isrc].c_str();
      TGraphErrors *gjec = (TGraphErrors*)ds2->Get(Form("g%s_%s",cm,cs));
      assert(gjec);
      string s = Form("%s_multijet_src%d", im==0 ? "ptchs" : "mpfchs1", isrc);
      TH1D *hs = new TH1D(Form("bm%d_%s",1<<im, s.c_str()),
                          Form("%s;p_{T}^{recoil};%s unc.",s.c_str(),cs),
                          npt, ptbins);

      for (int i = 0; i != gjec->GetN(); ++i) {

        double pt = gjec->GetX()[i];
        int ipt = hs->FindBin(pt);
        hs->SetBinContent(ipt, gjec->GetY()[i]);
        hs->SetBinError(ipt, gjec->GetY()[i]);
	
	// Set to percentage
	gjec->SetPoint(i, gjec->GetX()[i], 100.*gjec->GetY()[i]);
	gjec->SetPointError(i, gjec->GetEX()[i], 100.*gjec->GetEY()[i]);
      } // for i     
      
      vsrc.push_back(hs);

      c1->cd();
      hs->Scale(100);

      if (string(cs)=="SystTot") {
	if (string(cm)=="MPF") {
	  tdrDraw(hs,"HIST][",kNone);
	  hs->DrawClone("HISTSAME][");
	  hs->Scale(-1);
	  hs->DrawClone("HISTSAME][");
	  leg->AddEntry(hs,"Total (ex. FSR)","F");
	}
      }
      else {
	tdrDraw(gjec,"P",marker[cm][cs],color[cm][cs]);
	if (string(cm)=="MPF") leg->AddEntry(gjec,cs,"P");

	TF1 *f1 = new TF1(Form("f1_%s%s",cm,cs),
			  "[0]+[1]*log(x/200.)+[2]*log(x/200.)*log(x/200.)",
			  200,1500);
	f1->SetParameters(0.1,0.01,0.001);
	gjec->Fit(f1,"QRNW");
	f1->SetLineColor(color[cm][cs]);
	f1->Draw("SAME");

	//cout << endl;
	cout << Form("  %s_%sFunc (\"%s_%sFunc\","
		     "\"0.01*(%+6.3f + %6.4f*log(x/200.)"
		     " %+8.5f*log(x/200.)*log(x/200.))\",10,7000),",
		     cm,cs,cm,cs,
		     f1->GetParameter(0), f1->GetParameter(1),
		     f1->GetParameter(2)) << endl;
	//cout << endl;
      }

    } // for isrc                                                               
  } // for im         

  // === CODE FROM MULTIJET.C === (end)

  // === CODE FROM GLOBAL FIT === (start)

  // New uncertainty sources for multijets
  //if (false && string (samples[0])=="multijet") {
  TH1D *hfsr(0);
  if (true) {

    //const int isample = 0;
    const int nmethods = 2;
    for (int imethod = 0; imethod != nmethods; ++imethod) {

      const char *cm = (imethod==0 ? "MPF" : "MJB"); //methods[imethod];
      const char *cs = "FSR";//"multijet";

      //const int npt = 25;
      double vpt[]= {200, 250, 300, 330, 370, 410, 450, 510, 530, 550, 575,
		     600, 650, 700, 800, 900, 1000, 1100, 1200, 1400, 1600,
		     1800, 2000, 2200, 3500};
      const int npt = sizeof(vpt)/sizeof(vpt[0])-1;
      TH1D *h = new TH1D(Form("hmj_%s",cm),"",npt,vpt);

      // Pt10/Pt30 fits with drawMultijet.C                                     
      // Using 10/30 rather than 20/30, because this shows more effect          
      // at low pT where MPF and MJB differ for Pt30                            
      //MJB:                                                                    
      double vmjb[4] = {1.8054e+06, -3.2847, -1.0056e+09, -4.3954};
      //MPF:                                                                    
      double vmpf[4] = {1.221e+05, -2.7955, -2.5519e+08, -4.1888};
      TF1 *f1 = new TF1(Form("f1_%s",cm),
                        "1+[0]*pow(x,[1])+[2]*pow(x,[3])",200,1600);
      //f1->SetParameters(string(cm)=="ptchs" ? vmjb : vmpf);
      f1->SetParameters(string(cm)=="MJB" ? vmjb : vmpf);
      for (int i = 1; i != h->GetNbinsX()+1; ++i) {
        h->SetBinContent(i, f1->Eval(h->GetBinCenter(i))-1);
        h->SetBinError(i, 0);
      }

      int ibm = 0;//isample + nsamples*imethod;
      h->SetName(Form("bm%d_%s",(1<<ibm),h->GetName()));

      //hs.push_back(h);

      c1->cd();
      h->Scale(100);
      tdrDraw(h,"P",marker[cm][cs],color[cm][cs]);
      //if (string(cm)=="MPF") leg->AddEntry(h,cs,"P");

      if (hfsr)    {
	tdrDraw(hfsr,"P",kFullTriangleDown,color[cm][cs]);
	leg->AddEntry(hfsr,cs,"P");
      }
      if (hfsr==0) hfsr = (TH1D*)h->Clone("hfsr");
      else         {

	hfsr->Add(h,-1);

	TF1 *f1 = new TF1(Form("f1_%s%s",cm,cs),
			  "[0]+[1]*pow(x/200.,[2])",
			  200,3500);
	f1->SetParameters(0.,1,-0.5);
	hfsr->Fit(f1,"QRNW");
	f1->SetLineColor(color[cm][cs]);
	f1->DrawClone("SAME");

	cout << endl;
	cout << Form("  MPF_FSRFunc (\"MPF_FSRFunc\","
		     "\"0.01*(%5.3f + %5.3f*pow(x/200.,%+7.4f))\",10,7000),",
		     f1->GetParameter(0), f1->GetParameter(1),
		     f1->GetParameter(2)) << endl;
	cout << Form("  MJB_FSRFunc (\"MJB_FSRFunc\","
		     "\"0.01*(%5.3f + %5.3f*pow(x/200.,%+7.4f))\",10,7000)",
		     2.*f1->GetParameter(0), 2.*f1->GetParameter(1),
		     f1->GetParameter(2)) << endl;
	cout << endl;

	f1->SetParameters(f1->GetParameter(0)*2,
			  f1->GetParameter(1)*2,
			  f1->GetParameter(2));
	f1->DrawClone("SAME");
      }

    } // for imethod
  } // for ieig              

  // (NB: may have had MJB and MPF uncertainties swapped...)
  // === CODE FROM GLOBAL FIT === (start)
 

  TFile *fm = new TFile("../rootfiles/jecdataGH_Multijet.root");
  assert(fm && !fm->IsZombie());

  //TH1D *hmjb = (TH1D*)fm->Get("MJB_RatioPostFit"); assert(hmjb);
  //TH1D *hmpf = (TH1D*)fm->Get("MPF_RatioPostFit"); assert(hmpf);
  TH1D *hmjb = (TH1D*)fm->Get("MJB_RatioRaw"); assert(hmjb);
  TH1D *hmpf = (TH1D*)fm->Get("MPF_RatioRaw"); assert(hmpf);

  TH1D *h1 = (TH1D*)hmjb->Clone("h1");
  h1->Divide(hmjb);
  hmjb->Add(h1,-1); hmjb->Scale(100.); 
  hmpf->Add(h1,-1); hmpf->Scale(100.);

  c1->cd();
  hmjb->GetXaxis()->SetRangeUser(200,1000);//500);
  hmpf->GetXaxis()->SetRangeUser(200,1000);//500);
  //tdrDraw(hmpf,"P",kFullCircle);
  //tdrDraw(hmjb,"P",kOpenCircle);
  hmjb->Add(hmpf,-1);
  tdrDraw(hmjb,"HIST",kNone,kBlack,kDashed,-1,kNone);

  gPad->SetLogx();
  gPad->RedrawAxis();
  leg->Draw();

  c1->SaveAs("../pdf/drawMJBunc.pdf");

  // Load Andrey's new systematics
  TFile *f = new TFile("../../jecsys/rootfiles/"
		       "multijet_20161202_Run2016FlateG.root","READ");
  // Data/Pt10,Pt15,Pt20,Pt30
  // MC, MC_JER-up, MC_JER-down

}
