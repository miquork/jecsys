#include "../tdrstyle_mod15.C"
#include "../tools.C"

#include "TFile.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TLine.h"

bool ispr = true;

void drawZmass() {

  setTDRStyle();

  TFile *f = new TFile("../rootfiles/jecdataBCDEFGH.root","READ");
  //TFile *f = new TFile("../rootfiles/jecdataGH.root","READ");
  //TFile *f = new TFile("../rootfiles/jecdataBCD.root","READ");
  assert(f && !f->IsZombie());

  TGraphErrors *gzmmd = (TGraphErrors*)f->Get("data/eta00-13/mass_zmmjet_a30");
  assert(gzmmd);
  TGraphErrors *gzmmm = (TGraphErrors*)f->Get("mc/eta00-13/mass_zmmjet_a30");
  assert(gzmmm);
  TGraphErrors *gzmmr = (TGraphErrors*)f->Get("ratio/eta00-13/mass_zmmjet_a30");
  assert(gzmmr);

  TGraphErrors *gzeed = (TGraphErrors*)f->Get("data/eta00-13/mass_zeejet_a30");
  assert(gzeed);
  TGraphErrors *gzeem = (TGraphErrors*)f->Get("mc/eta00-13/mass_zeejet_a30");
  assert(gzeem);
  TGraphErrors *gzeer = (TGraphErrors*)f->Get("ratio/eta00-13/mass_zeejet_a30");
  assert(gzeer);

  double ptmax = 1200.*2.;
  TH1D *hup = new TH1D("hup",";p_{T,Z} (GeV);m_{Z} (GeV)",670,30,ptmax);
  hup->SetMinimum(90.5);//88.5);
  hup->SetMaximum(92.5);//94.5);
  hup->GetXaxis()->SetMoreLogLabels();
  hup->GetXaxis()->SetNoExponent();

  TH1D *hdw = new TH1D("hdw",";p_{T,Z} (GeV);Data / MC - 1 (%)",670,30,ptmax);
  hdw->SetMinimum((0.990+1e-5-1)*100);
  //hdw->SetMaximum(1.010+1e-5);
  hdw->SetMaximum((1.025+1e-5-1)*100);
  hdw->GetXaxis()->SetMoreLogLabels();
  hdw->GetXaxis()->SetNoExponent();

  //lumi_13TeV = "Run2016BCDEFGH 36.5 fb^{-1}";
  lumi_13TeV = "Run2016 36.5 fb^{-1}";// Private Work";
  extraText = "Private Work";
  //TCanvas *c1 = tdrDiCanvas("c1",hup,4,11),kSquare);
  //TCanvas *c1 = tdrDiCanvas("c1",hup,hdw,4,11);
  TCanvas *c1 = tdrDiCanvas("c1",hdw,hup,4,11);

  //c1->cd(1);
  c1->cd(2);
  gPad->SetLogx();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  double mzpdg = 91.2;
  l->DrawLine(30,mzpdg,ptmax,mzpdg);

  tdrDraw(gzeed,"Pz",kFullSquare,kGreen+2);
  tdrDraw(gzeem,"Pz",kOpenSquare,kGreen+2);

  tdrDraw(gzmmd,"Pz",kFullCircle,kRed);
  tdrDraw(gzmmm,"Pz",kOpenCircle,kRed);


  //c1->cd(2);
  c1->cd(1);
  gPad->SetLogx();

  // Plot originals on the back to see the default fits
  //tdrDraw(gzeer,"Pz",kFullSquare,kGreen+2);
  //tdrDraw(gzmmr,"Pz",kFullCircle,kRed);

  // Plot the fit used in reprocess.C
  TF1 *f1jec = new TF1("f1jec","(([0]+[1]*log(0.01*x)"
			"+[2]*pow(log(0.01*x),2))"
			"-1)*100", 30, ptmax);
  //TF1 *f1ezee = new TF1("f1ezee","sqrt([0]+pow(log(2*x),2)*[1]"
  //			"+pow(log(2*x),4)*[2]+2*log(2*x)*[3]"
  //			"+2*pow(log(2*x),2)*[4]+2*pow(log(2*x),3)*[5])"
  //			,30,ptmax);
  // BCDEFGH fit with minitools/drawZmass.C
  //f1mzee->SetParameters(1.01732, -0.00909, 0.00116, 0.99855);
  f1jec->SetParameters(0.99885, 0.00176, 0.00135);//EGM1
  //f1jec->SetParameters(1.00017, 0.00166, 0.00114);//EGM2
  //f1jec->SetParameters(1.00279, 0.00166, 0.00112);//EGM3
  //f1ezee->SetParameters(7.54e-05, 1.41e-05, 1.63e-07,
  //			-3.26e-05, 3.47e-06, -1.51e-06);
  f1jec->SetLineColor(kBlack);
  f1jec->SetLineWidth(2);
  f1jec->Draw("SAME");
  

  l->DrawLine(30,0,ptmax,0);

  // Recreate ratios, because GetListOfFuctions segfaults
  // Also have small shift in pT,Z values within bins
  TGraphErrors *gzeer2 = tools::ratioGraphs(gzeed,gzeem); 
  TGraphErrors *gzmmr2 = tools::ratioGraphs(gzmmd,gzmmm);

  // Increase uncertainty for Zmm mZ kink at pT,Z~70 GeV
  int ik = 3;
  //double ek = 0.00000; // no extra
  //double ek = 0.0004; // 0.04% extra
  //double ek = 0.0002; // 0.02% extra
  // same as ee * mumu(100)/ee(100):
  double ek = gzeer2->GetEY()[3] * gzmmr2->GetEY()[4]/gzeer2->GetEY()[4];
  gzmmr2->SetPointError(ik, gzmmr2->GetEX()[ik],
			sqrt(pow(gzmmr2->GetEY()[ik],2)+pow(ek,2)));

  // Change result to (data/MC-1)*100
  for (int i = 0; i != gzeer2->GetN(); ++i) {
    gzeer2->SetPoint(i, gzeer2->GetX()[i], (gzeer2->GetY()[i]-1)*100);
    gzeer2->SetPointError(i, gzeer2->GetEX()[i], (gzeer2->GetEY()[i])*100);
  }
  for (int i = 0; i != gzmmr2->GetN(); ++i) {
    gzmmr2->SetPoint(i, gzmmr2->GetX()[i], (gzmmr2->GetY()[i]-1)*100);
    gzmmr2->SetPointError(i, gzmmr2->GetEX()[i], (gzmmr2->GetEY()[i])*100);
  }

  TF1 *f1mzee = new TF1("f1zee","([0]+[1]*log(0.01*x)+[2]*pow(log(0.01*x),2)-1)*100",
			30,ptmax);
  f1mzee->SetParameters(1,0.001,0.0001);
  f1mzee->SetLineColor(kGreen+2);
  gzeer2->Fit(f1mzee,"QRN");
  f1mzee->Draw("SAME");

  TMatrixD emat(f1mzee->GetNpar(),f1mzee->GetNpar());
  gMinuit->mnemat(&emat[0][0],f1mzee->GetNpar());

  // d/dp0=1, d/dp1=log(x), d/dp2=log(x)^2
  // eps =   (d/dp0)^2 m_00 +   (d/dp1)^2 m_11   +   (d/dp2)^2 m_22
  //     + 2*(d/dp0/p1)m_01 + 2*(d/dp0/p2)m_02 + 2*(d/dp1/p2)m_12
  TF1 *f1e = new TF1("f1e","([0]+[1]*log(0.01*x)+[2]*pow(log(0.01*x),2)"
		     "+ [3]*sqrt([4] + pow(log(0.01*x),2)*[5]"
		     "+ pow(log(0.01*x),4)*[6]"
		     "+ 2*log(0.01*x)*[7]+2*pow(log(0.01*x),2)*[8]"
		     "+ 2*pow(log(0.01*x),3)*[9])"
		     " - 1)*100",30,ptmax);
  f1e->SetParameters(f1mzee->GetParameter(0),f1mzee->GetParameter(1),
		     f1mzee->GetParameter(2), +1,
		     emat[0][0], emat[1][1], emat[2][2],
		     emat[0][1], emat[0][2], emat[1][2]);

  f1e->SetLineColor(kGreen-8);
  f1e->DrawClone("SAME");
  f1e->SetParameter(3, -1);
  f1e->DrawClone("SAME");

  // Print out corrections for reprocess.C
  cout <<
    "  // Smoothen mass corrections\n"
    "  TF1 *f1mzee = new TF1(\"f1mzee\",\"[0]+[1]*log(0.01*x)+[2]*pow(log(0.01*x),2)\",\n"
    "			     fzeeptmin, fzeeptmax);\n"
    "  TF1 *f1ezee = new TF1(\"f1ezee\",\"sqrt([0]+pow(log(0.01*x),2)*[1]\"\n"
    "                        \"+pow(log(0.01*x),4)*[2]\"\n"
    "                        \"+2*log(0.01*x)*[3]+2*pow(log(0.01*x),2)*[4]\"\n"
    "                        \"+2*pow(log(0.01*x),3)*[5])\"\n"
    "			     ,fzeeptmin, fzeeptmax);\n"
    "  if (correctZeeMass || correctGamMass) {\n"
    "    if (useFixedFit) {\n"
    "      // BCDEFGH fit with minitools/drawZmass.C\n";
  cout << Form(""
	       "      f1mzee->SetParameters(%1.5f, %1.5f, %1.5f);\n"
	       "      f1ezee->SetParameters(%+9.3g, %+9.3g, %+9.3g,\n"
	       "			    %+9.3g, %+9.3g, %+9.3g);\n"
	       "    }\n"
	       "    else\n"
	       "      hmzee->Fit(f1mzee);\n"
	       "\n"
	       "  }\n",
	       f1mzee->GetParameter(0),f1mzee->GetParameter(1),
	       f1mzee->GetParameter(2),
	       emat[0][0], emat[1][1], emat[2][2],
	       emat[0][1], emat[0][2], emat[1][2]);


  //TF1 *f1mzee0 = new TF1("f1zee0","1+[0]*log(x)+[1]*log(x)*log(x)",30,ptmax);
  TF1 *f1mzee0 = new TF1("f1zee0","(1+[0]*log(x)+[1]*log(x)*log(x)-1)*100",
			 30,ptmax);
  f1mzee0->SetParameters(0.001,0.0001);
  //TF1 *f1mzee0 = new TF1("f1zee0","[0]+[1]*pow(x,[2])",30,ptmax);
  //f1mzee0->SetParameters(1,0.01,0.5);
  f1mzee0->SetLineColor(kGreen+3);
  gzeer2->Fit(f1mzee0,"QRN");
  //f1mzee0->Draw("SAME");

  TF1 *f1mzmm = new TF1("f1zmm","([0]+[1]*log(0.01*x)+[2]*pow(log(0.01*x),2)-1)*100",
			30,ptmax);
  //f1mzmm->SetParameters(1,-0.001,0.00015);
  //f1mzmm->FixParameter(1,1);
  f1mzmm->SetParameters(1,0,0);
  f1mzmm->FixParameter(2,0);
  f1mzmm->FixParameter(1,0);
  f1mzmm->SetLineColor(kRed);
  gzmmr2->Fit(f1mzmm,"QRN");
  f1mzmm->Draw("SAME");

  TMatrixD emat2(f1mzmm->GetNpar(),f1mzmm->GetNpar());
  gMinuit->mnemat(&emat2[0][0],f1mzmm->GetNpar());
  
  f1e->SetLineStyle(kSolid);
  f1e->SetLineColor(kRed-9);
  f1e->SetParameters(f1mzmm->GetParameter(0),0,0,+1,
		     pow(f1mzmm->GetParError(0),2),0,0, 0,0,0);
  f1e->DrawClone("SAME");
  f1e->SetParameter(3,-1);
  f1e->DrawClone("SAME");

  // Print out corrections for reprocess.C
  cout <<
    "\n"
    "  TF1 *f1mzmm = new TF1(\"f1mzmm\",\"[0]+[1]*log(0.01*x)+[2]*pow(log(0.01*x),2)\",\n"
    "			     fzmmptmin, fzmmptmax);\n"
    "  TF1 *f1ezmm = new TF1(\"f1ezmm\",\"sqrt([0]+pow(log(0.01*x),2)*[1]\"\n"
    "                        \"+pow(log(0.01*x),4)*[2]\"\n"
    "                        \"+2*log(0.01*x)*[3]+2*pow(log(0.01*x),2)*[4]\"\n"
    "                        \"+2*pow(log(0.01*x),3)*[5])\"\n"
    "			     ,fzeeptmin, fzeeptmax);\n"
    "  if (correctZmmMass) {\n"
    "    if (useFixedFit) {\n"
    "      // BCDEFGH fit with minitools/drawZmass.C\n";
  cout << Form(""
	       "      f1mzmm->SetParameters(%1.5f, %1.5f, %1.5f);\n"
	       "      f1ezmm->SetParameters(%+9.3g, %+9.3g, %+9.3g,\n"
	       "			    %+9.3g, %+9.3g, %+9.3g);\n"
	       "    }\n"
	       "    else\n"
	       "      hmzmm->Fit(f1mzmm);\n"
	       "\n"
	       "  }\n",
	       f1mzmm->GetParameter(0),f1mzmm->GetParameter(1),
	       f1mzmm->GetParameter(2),
	       emat2[0][0], emat2[1][1], emat2[2][2],
	       emat2[0][1], emat2[0][2], emat2[1][2]);

  // Draw points last so they stay on top of fits
  tdrDraw(gzeer2,"Pz",kFullSquare,kGreen+2); gzeer2->SetMarkerSize(0.5);
  tdrDraw(gzmmr2,"Pz",kFullCircle,kRed);     gzmmr2->SetMarkerSize(0.5);

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);

  // Don't print all the internal information to the PR plots
  if (!ispr) {

    tex->SetTextColor(kGreen+2);
    tex->DrawLatex(0.20,0.40,Form("#chi^{2} / NDF = %1.1f / %d (p3)",
				  f1mzee->GetChisquare(), f1mzee->GetNDF()));
    tex->SetTextColor(kGreen+2);
    tex->DrawLatex(0.20,0.70,Form("p_{0} = %1.5f #pm %1.5f",
				  f1mzee->GetParameter(0),
				  f1mzee->GetParError(0)));
    tex->DrawLatex(0.20,0.65,Form("p_{1} = %1.5f #pm %1.5f",
				  f1mzee->GetParameter(1),
				  f1mzee->GetParError(1)));
    tex->DrawLatex(0.20,0.60,Form("p_{2} = %1.5f #pm %1.5f",
				  f1mzee->GetParameter(2),
				  f1mzee->GetParError(2)));
    
    tex->SetTextSize(0.035);
    tex->DrawLatex(0.20,0.54,Form("%10.3g %10.3g %10.3g",
				  emat[0][0], emat[0][1], emat[0][2]));
    tex->DrawLatex(0.20,0.50,Form("%10.3g %10.3g %10.3g",
				  emat[1][0], emat[1][1], emat[1][2]));
    tex->DrawLatex(0.20,0.46,Form("%10.3g %10.3g %10.3g",
				  emat[2][0], emat[2][1], emat[2][2]));
    tex->SetTextSize(0.045);
    
    //tex->DrawLatex(0.20,1.09-dy,Form("p_{0} = %1.5f #pm %1.5f",
    //tex->DrawLatex(0.20,0.73-dy,Form("#chi^{2}/NDF = %1.1f/%d (p1)",

    tex->SetTextColor(kRed);
    tex->DrawLatex(0.20,0.05,Form("#chi^{2} / NDF = %1.1f / %d (p1)",
				  f1mzmm->GetChisquare(), f1mzmm->GetNDF()));
    tex->DrawLatex(0.20,0.11,Form("p_{0} = %1.5f #pm %1.5f",
				  f1mzmm->GetParameter(0),
				  f1mzmm->GetParError(0)));
  } // !ispr
  if (ispr) {

    TLegend *leg = tdrLeg(0.20,0.50,0.40,0.74);
    leg->AddEntry(gzeed,"Z(#rightarrow e^{+}e^{-}) data","PL");
    leg->AddEntry(gzmmd,"Z(#rightarrow #mu^{+}#mu^{-}) data","PL");
    leg->AddEntry(gzeem,"Z(#rightarrow e^{+}e^{-}) MC","PL");
    leg->AddEntry(gzmmm,"Z(#rightarrow #mu^{+}#mu^{-}) MC","PL");


    tex->SetTextColor(kGreen+2);
    tex->DrawLatex(0.2,0.40,Form("#chi^{2} / NDF = %1.1f / %d (e^{+}e^{-})",
				  f1mzee->GetChisquare(), f1mzee->GetNDF()));

    tex->SetTextColor(kRed);
    tex->DrawLatex(0.2,0.05,Form("#chi^{2} / NDF = %1.1f / %d (#mu^{+}#mu^{-})",
				 f1mzmm->GetChisquare(), f1mzmm->GetNDF()));
    tex->DrawLatex(0.2,0.11,Form("#DeltaM(#mu^{+}#mu^{-}) = %1.3f #pm %1.3f%%",
				 (f1mzmm->GetParameter(0)-1)*100,
				 f1mzmm->GetParError(0)*100));
  }


  c1->SaveAs("../pdf/drawZmass.pdf");
}
