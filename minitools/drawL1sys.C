// Purpose: draw L1 uncertainties stored in rootfiles/jecdataBCDEF.root
#include "TFile.h"
#include "TLine.h"
#include "TF1.h"

#include "../tdrstyle_mod15.C"

#include <string>
#include <vector>
#include <map>

using namespace std;

void drawL1sys() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TFile *f = new TFile("rootfiles/jecdataBCDEF.root","READ");
  assert(f && !f->IsZombie());

  string as[] = {"hsys1_rcmpf","hsys2a_l1sf","hsys2b_l1",
		 "hsys3a_l1sfr","hsys3b_l1r"};
  const int ns = sizeof(as)/sizeof(as[0]);
  map<string,Style_t> marker;
  marker["hsys1_rcmpf"] = kOpenDiamond;
  marker["hsys2a_l1sf"] = kFullSquare;
  marker["hsys2b_l1"] = kFullCircle;
  marker["hsys3a_l1sfr"] = kOpenSquare;
  marker["hsys3b_l1r"] = kOpenCircle;
  map<string,Color_t> color;
  color["hsys1_rcmpf"] = kBlack;
  color["hsys2a_l1sf"] = kRed;
  color["hsys2b_l1"] = kBlue;
  color["hsys3a_l1sfr"] = kRed-9;
  color["hsys3b_l1r"] = kBlue-9;
  map<string,const char*> label;
  label["hsys1_rcmpf"] = "L1RC_MED (1)";
  label["hsys2a_l1sf"] = "SF_MED (2a)";
  label["hsys2b_l1"] = "L1SemiSimple (2b)";
  label["hsys3a_l1sfr"] = "SF_MED (3a)";
  label["hsys3b_l1r"] = "L1SemiSimple (3b)";


  double ptmin = 15;
  double ptmax = 4500;
  TH1D *h = tdrHist("h","Offset uncertainty (%)",
                     -1.,+2,"p_{T} (GeV)",ptmin,ptmax);

  lumi_13TeV = "2017UL";
  TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);
  gPad->SetLogx();

  TLegend *leg = tdrLeg(0.50,0.65,0.80,0.90);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(ptmin,0,ptmax,0);

  
  for (int i = 0; i != ns; ++i) {
    string s = as[i];
    const char *c = s.c_str();
    TH1D *hs = (TH1D*)f->Get(Form("ratio/eta00-13/%s",c)); assert(hs);
    for (int i = 1; i != hs->GetNbinsX()+1; ++i) {
      hs->SetBinContent(i, 100.*(hs->GetBinContent(i)-1));
    }
    tdrDraw(hs,"P",marker[s],color[s],kSolid,-1,kNone);
    leg->AddEntry(hs,label[s],"LP");

    TF1 *f1 = new TF1(Form("f1_%s",hs->GetName()),
		      "100.*([0]/x+[1]*log(x)/x+[2]/(x*log(x)))",ptmin,ptmax);
    f1->SetParameters(0.1,0.01,0.1);
    f1->SetLineColor(color[s]);
    f1->SetLineWidth(2);
    hs->Fit(f1,"QRN");
    f1->Draw("SAME");

    cout << Form("TF1 *f1_%s = new TF1(\"f1_%s\","
		 "\"[0]/x+[1]*log(x)/x+[2]/(x*log(x))\",%d,%d);\n",
		 c,c,int(ptmin+0.5),int(ptmax+0.5));
  } //  for i

  gPad->RedrawAxis();

  c1->SaveAs("pdf/drawL1sys.pdf");
}
