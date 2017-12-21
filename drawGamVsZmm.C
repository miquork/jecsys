#include "TFile.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TMatrixD.h"

#include "tdrstyle_mod15.C"

void drawGamVsZmm() {

  setTDRStyle();

  //TFile *f = new TFile("rootfiles/jecdataBCDEFGH_lowpt.root","READ");
  TFile *f = new TFile("rootfiles/jecdataBCD.root","READ");
  //TFile *f = new TFile("rootfiles/jecdataGH.root","READ");
  //TFile *f = new TFile("rootfiles/jecdataBCDEFGH.root","READ");
  assert(f && !f->IsZombie());

  string dir = "data/eta00-13";
  const char *cd = dir.c_str();
  TGraphErrors *gp = (TGraphErrors*)f->Get(Form("%s/mpfchs1_gamjet_a30",cd));
  assert(gp);
  TGraphErrors *gz = (TGraphErrors*)f->Get(Form("%s/mpfchs1_zmmjet_a30",cd));
  assert(gz);
  TGraphErrors *ge = (TGraphErrors*)f->Get(Form("%s/mpfchs1_zeejet_a30",cd));
  assert(ge);

  TH1D *h = (TH1D*)f->Get(Form("%s/h",cd)); assert(h);
  h->SetMaximum(1.01-1e-4);
  h->SetMinimum(0.95+1e-4);
  h->GetXaxis()->SetRangeUser(30.,1500.);

  TH1D *h2 = (TH1D*)h->Clone("h2");
  h2->SetYTitle("#gamma+jet / Z#mu#mu+jet");
  h2->SetMaximum(1.00+1e-5);
  h2->SetMinimum(0.98-1e-5);

  //TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);
  lumi_13TeV = "Run2016BCDEFGH Legacy, 36.5 fb^{-1}";
  TCanvas *c1 = tdrDiCanvas("c1",h,h2,4,11);//,kSquare);

  c1->cd(1);
  gPad->SetLogx();

  //ge->Draw("SAMEPz");
  tdrDraw(ge,"Pz",kFullDiamond,kGreen+2);
  //gp->Draw("SAMEP");
  tdrDraw(gp,"Pz",kFullSquare,kBlue); 
  //gz->Draw("SAMEP");
  tdrDraw(gz,"Pz",kFullCircle,kRed); 

  TLegend *leg = tdrLeg(0.65,0.12,0.90,0.30);
  leg->AddEntry(gz,"Z#mu#mu+jet","P");
  leg->AddEntry(ge,"Zee+jet","P");
  leg->AddEntry(gp,"#gamma+jet","P");

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);
  tex->DrawLatex(0.20,0.70,"Type-I MPF+CHS");
  tex->DrawLatex(0.20,0.64,"#alpha < 0.30");


  c1->cd(2);
  gPad->SetLogx();

  for (int i = gz->GetN()-1; i != -1; --i) {
    //if (gz->GetX()[i]<40. || gz->GetX()[i]>500.)
    if (gz->GetX()[i]<gp->GetX()[0]*0.9 || gz->GetX()[i]>500.)
      gz->RemovePoint(i);
  }
  for (int i = ge->GetN()-1; i != -1; --i) {
    //if (ge->GetX()[i]<40. || ge->GetX()[i]>500.)
    if (ge->GetX()[i]<gp->GetX()[0]*0.9 || ge->GetX()[i]>500.)
      ge->RemovePoint(i);
  }
  //gz->RemovePoint(0);

  assert(gp->GetN()>=gz->GetN());
  TGraphErrors *gr = (TGraphErrors*)gz->Clone("gr");
  TGraphErrors *gr2 = (TGraphErrors*)ge->Clone("gr2");
  for (int i = 0; i != gz->GetN(); ++i) {
    if (fabs(gz->GetX()[i]/gp->GetX()[i]-1)<0.1 && gp->GetX()[i]>40.) {
      gr->SetPoint(i, 0.5*(gz->GetX()[i]+gp->GetX()[i]),
		   gp->GetY()[i]/gz->GetY()[i]);
      gr->SetPointError(i, 0.5*(gz->GetX()[i]-gp->GetX()[i]),
			gp->GetY()[i]/gz->GetY()[i]*
			sqrt(pow(gz->GetEY()[i]/gz->GetY()[i],2) +
			     pow(gp->GetEY()[i]/gp->GetY()[i],2)));
      //
      gr2->SetPoint(i, 0.5*(gz->GetX()[i]+ge->GetX()[i]),
		   ge->GetY()[i]/gz->GetY()[i]);
      gr2->SetPointError(i, 0.5*(gz->GetX()[i]-ge->GetX()[i]),
			ge->GetY()[i]/gz->GetY()[i]*
			sqrt(pow(gz->GetEY()[i]/gz->GetY()[i],2) +
			     pow(ge->GetEY()[i]/ge->GetY()[i],2)));
    } // if matching x
  } // for i

  tdrDraw(gr2,"Pz",kFullDiamond,kGreen+2);
  tdrDraw(gr,"Pz",kFullCircle,kBlue);
  gr->SetLineWidth(2);

  TF1 *f1 = new TF1("f1","[0]+[1]*log(x)",60,500);
  f1->SetLineColor(kBlue);
  f1->SetLineWidth(2);
  gr->Fit(f1,"QRN");
  f1->Draw("SAME");

  TMatrixD emat(2,2);
  gMinuit->mnemat(&emat[0][0],2);

  TF1 *f1eup = new TF1("f1eup","[0]+[1]*log(x)+"
		       "[2]*sqrt([3]*1+2*[4]*log(x)+[5]*log(x)*log(x))",60,500);
  TF1 *f1edw = new TF1("f1edw","[0]+[1]*log(x)+"
		       "[2]*sqrt([3]*1+2*[4]*log(x)+[5]*log(x)*log(x))",60,500);
  f1eup->SetLineColor(kBlue);
  f1eup->SetLineStyle(kDashed);
  f1eup->SetParameters(f1->GetParameter(0),f1->GetParameter(1), +1,
		     emat[0][0],emat[0][1],emat[1][1]);
  f1eup->Draw("SAME");

  f1edw->SetLineColor(kBlue);
  f1edw->SetLineStyle(kDashed);
  f1edw->SetParameters(f1->GetParameter(0),f1->GetParameter(1), -1,
		       emat[0][0],emat[0][1],emat[1][1]);
  f1edw->Draw("SAME");

  c1->SaveAs("pdf/drawGamVsZmm.pdf");
}
