// compare results from two jecdata files
#include "TFile.h"
#include "TGraphErrors.h"

#include "tdrstyle_mod15.C"

#include <string>

void compareJECdata(string dir = "mc", 
		    string ch = "zmmjet",
		    string var = "mpfchs1",
		    string alpha = "a30",
		    string eta = "eta00-13",
		    //string file1 = "rootfiles/jecdataG_80XSum16V5.root",
		    ///string label1 = "Gv5",
		    //string file1 = "rootfiles/jecdataG.root",
		    //string file1 = "rootfiles/jecdataG_GammaJetV6.root",
		    //string file1 = "rootfiles/jecdataG_80XSum16V6.root",
		    //string file1 = "rootfiles/jecdataG_RSv6.root",
		    string file1 = "rootfiles/jecdataH-JetMET100-fix.root",
		    //string label1 = "Gv6",
		    //string label1 = "GRSv6",
		    string label1 = "JME100",
		    //string file2 = "rootfiles/jecdataG_TBv6.root",
		    //string file2 = "rootfiles/jecdataG_80XSum16V4.root",
		    //string file2 = "rootfiles/jecdataG_GammaJetV4.root",
		    string file2 = "rootfiles/jecdataH-Summer16_23Sep2016V4.root",
		    //string label2 = "Gv4"
		    //string label2 = "GTBv6"
		    string label2 = "Sum16V4"
		    ) {
  //string file2 = "rootfiles/jecdataG.root",
  //		    string label2 = "Gv5_08") {

  TDirectory *curdir = gDirectory;
  setTDRStyle();

  TFile *fin1 = new TFile(file1.c_str(),"READ");
  assert(fin1 && !fin1->IsZombie());

  TFile *fin2 = new TFile(file2.c_str(),"READ");
  assert(fin2 && !fin2->IsZombie());
  
  string s = Form("%s/%s/%s_%s_%s", dir.c_str(), eta.c_str(),
		  var.c_str(), ch.c_str(), alpha.c_str());
  cout << s << endl;
  TGraphErrors *g1 = (TGraphErrors*)fin1->Get(s.c_str());
  assert(g1);
  TGraphErrors *g2 = (TGraphErrors*)fin2->Get(s.c_str());
  assert(g2);

  //assert(g1->GetN()==g2->GetN());
  // Check if numbers of points match
  if (!(g1->GetN()==g2->GetN())) {

    // Print out numbers of points and their ranges
    cout << label1 << ": " << g1->GetN()
	 << "  (" << g1->GetX()[0]
	 << " - " << g1->GetX()[g1->GetN()-1] << ")" << endl;
    cout << label2 << ": " << g2->GetN()
	 << "  (" << g2->GetX()[0]
	 << " - " << g2->GetX()[g2->GetN()-1] << ")" << endl;
    cout << flush;

    //g1->RemovePoint(g1->GetN()-1);
    // Remove points of g1 not in g2
    for (int i = g1->GetN()-1; i != -1; --i) {
      double x = g1->GetX()[i];
      double dxmin(fabs(g1->GetX()[0]-g1->GetX()[g1->GetN()-1]));
      double dx(dxmin);
      for (int j = 0; j != g2->GetN(); ++j) {
	double dx = fabs(x-g2->GetX()[j]);
	if (dx < dxmin) dxmin = dx;
      } // for j
      
      if (dxmin > 0.5*fabs(g1->GetX()[max(0,i+1)]-x)  &&
	  dxmin > 0.5*fabs(g1->GetX()[min(g1->GetN()-1,i-1)]-x))
	g1->RemovePoint(i);
    } // for i
    // Remove points of g2 not in g1
    for (int i = g2->GetN()-1; i != -1; --i) {
      double x = g2->GetX()[i];
      double dxmin(fabs(g2->GetX()[0]-g2->GetX()[g2->GetN()-1]));
      double dx(dxmin);
      for (int j = 0; j != g1->GetN(); ++j) {
	double dx = fabs(x-g1->GetX()[j]);
	if (dx < dxmin) dxmin = dx;
      } // for j
      
      if (dxmin > 0.5*fabs(g2->GetX()[max(0,i+1)]-x)  &&
	  dxmin > 0.5*fabs(g2->GetX()[min(g2->GetN()-1,i-1)]-x))
	g2->RemovePoint(i);
    } // for i

    // Print out cleaned points
    cout << label1 << ": " << g1->GetN()
	 << "  (" << g1->GetX()[0]
	 << " - " << g1->GetX()[g1->GetN()-1] << ")" << endl;
    cout << label2 << ": " << g2->GetN()
	 << "  (" << g2->GetX()[0]
	 << " - " << g2->GetX()[g2->GetN()-1] << ")" << endl;
    cout << flush;

    assert(g1->GetN()==g2->GetN());
  }
  TGraphErrors *gr = (TGraphErrors*)g1->Clone("gr");

  for (int i = 0; i != gr->GetN(); ++i) {
    gr->SetPoint(i, 0.5*(g1->GetX()[i]+g2->GetX()[i]),
		 g1->GetY()[i]/g2->GetY()[i]);
    gr->SetPointError(i, 0.5*fabs(g1->GetX()[i]-g2->GetX()[i]),
		      0.5*(g1->GetEY()[i] + g2->GetEY()[i]));
  }
  
  const double ptmin = (ch=="multijet" ? 200 : 30);
  const double ptmax = (ch=="multijet" ? 3000 : 1300);
  
  TH1D *hup = new TH1D("hup", Form(";p_{T} (GeV);Response (%s)",dir.c_str()),
		       1270,ptmin,ptmax);
  hup->GetXaxis()->SetMoreLogLabels();
  hup->GetXaxis()->SetNoExponent();
  hup->SetMinimum(ch=="multijet" ? 0.975 : 0.84);
  hup->SetMaximum(ch=="multijet" ? 1.095 : 1.09);

  TH1D *hdw = new TH1D("hdw", Form(";p_{T} (GeV);%s / %s",
				   label1.c_str(), label2.c_str()),
		       1270,ptmin,ptmax);
  hdw->GetXaxis()->SetMoreLogLabels();
  hdw->GetXaxis()->SetNoExponent();
  hdw->SetMinimum(ch=="multijet" ? 0.99-1e-4 : 0.95);
  hdw->SetMaximum(ch=="multijet" ? 1.01+1e-4 : 1.05);
  
  //lumi_13TeV = "Run2016FG re-reco, 8.0 fb^{-1}";
  lumi_13TeV = "Run2016H re-mAOD, 8.8 fb^{-1}";
  TCanvas *c1 = tdrDiCanvas("c1",hup,hdw,4,11);

  c1->cd(1);
  gPad->SetLogx();
  tdrDraw(g2,"Pz",kOpenCircle);
  tdrDraw(g1,"Pz",kFullCircle);

  TLegend *leg = tdrLeg(0.50,0.72,0.70,0.90);
  leg->SetHeader(Form("%s %s %s",var.c_str(),ch.c_str(), alpha.c_str()));
  leg->AddEntry(g1,label1.c_str(),"PL");
  leg->AddEntry(g2,label2.c_str(),"PL");
  
  c1->cd(2);
  gPad->SetLogx();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(ptmin,1,ptmax,1);

  tdrDraw(gr,"Pz",kFullCircle);

  c1->SaveAs(Form("pdf/compareJECdata_%svs%s_%s_%s_%s_%s.pdf",
		  label1.c_str(), label2.c_str(),
		  dir.c_str(), ch.c_str(), var.c_str(), alpha.c_str()));
} // compareJECdata
