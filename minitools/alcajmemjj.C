// Purpose: fit AlCaRaw hmjji13 (or hmjj13) to find W>qq' and Z>qqbar
#include "../tdrstyle_mod22.C"

void alcajmemjj() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  //TFile *f = new TFile("rootfiles/alcajme_out_50M_v6_mjj.root","READ");
  //TFile *f = new TFile("rootfiles/alcajme_out_50M_v8_doak4.root","READ");
  TFile *f = new TFile("rootfiles/alcajme_out_50M_v10_4p86fb.root","READ");
  assert(f && !f->IsZombie());

  curdir->cd();

  TH1D *hd = (TH1D*)f->Get("hmjji13"); assert(hd);

  //TF1 *f1 = new TF1("f1","[0]*pow(x,[1])*exp([2]/x)",40,300);
  //f1->SetParameters(1e16,-5, 0.);// 1e7,80.4,10.);
  //TF1 *f1 = new TF1("f1","[0]*pow(x,[1])",40,300);
  //f1->SetParameters(1e16,-5);

  TF1 *f1 = new TF1("f1","[0]*pow(x,[1]+[2]*log(x/100.)"
		    "+[3]*pow(log(x/100.),2)"
		    "+[4]*pow(log(x/100.),3)"
		    "+[5]*pow(log(x/100.),4)"
		    "+[6]*pow(log(x/100.),5)"
		    "+[7]*pow(log(x/100.),6)"
		    "+[8]*pow(log(x/100.),7))",
		    40,300);
  f1->SetParameters(1.05e16,-5.46,-0.057,0.099,0.108,-0.082,-0.025,0.01,0.01);
  f1->SetNpx(2600);

  // Blind 70-90 GeV range
  TH1D *hb = (TH1D*)hd->Clone("hb");
  for (int i = 1; i != hb->GetNbinsX()+1; ++i) {
    double pt = hb->GetBinCenter(i);
    //if (pt>70. && pt<105.) {
    if (pt>70. && pt<110.) {
      hb->SetBinContent(i, 0.);
      hb->SetBinError(i, 0.);
    }
  }

  TH1D *h1 = tdrHist("h1","N_{pairs}",3e2,1e9,
		     "Dijet mass m_{jj} (GeV)",30,300);
  //lumi_136TeV = "AlCaRaw DCSOnly, X fb^{-1}";
  lumi_136TeV = "AlCaRaw 1.79, Golden 4.86 fb^{-1}";
  TCanvas *c1 = tdrCanvas("c1",h1,8,11,kSquare);

  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextSize(0.035);
  t->DrawLatex(0.50,0.72,"All jet pairs");
  t->DrawLatex(0.50,0.67,"with |#Delta#eta|<1.3");

  hd->Draw("SAME");
  hb->SetLineColor(kGreen+2);
  hb->Draw("SAME");
  //hd->Fit(f1,"RN");
  hb->Fit(f1,"RN");
  f1->Draw("SAME");
  //hd->GetXaxis()->SetRangeUser(30.,400.);
  //gPad->SetLogx();
  gPad->SetLogy();

  cout << Form("Chi2/NDF = %1.1f / %d",
	       f1->GetChisquare(), f1->GetNDF()) << endl;

  TH1D *hdf = (TH1D*)hd->Clone("hdf");
  TH1D *hdfa = (TH1D*)hd->Clone("hdfa");
  for (int i = 1; i != hd->GetNbinsX()+1; ++i) {
    double pt = hd->GetBinCenter(i);
    double n = hd->GetBinContent(i);
    double en = hd->GetBinError(i);
    if (pt>32.) {
      hdf->SetBinContent(i, n>1 ? (n-f1->Eval(pt))/n : 0.);
      hdf->SetBinError(i, n>1 ? en/n : 0.);
      hdfa->SetBinContent(i, n>1 ? (n-f1->Eval(pt)) : 0.);
      hdfa->SetBinError(i, n>1 ? en : 0.);
    }
    else {
      hdf->SetBinContent(i, 0.);
      hdf->SetBinError(i, 0.);
      hdfa->SetBinContent(i, 0.);
      hdfa->SetBinError(i, 0.);
    }
  }


  TH1D *h2 = tdrHist("h2","(Data - fit) / fit",-0.07,0.07,
		    "Dijet mass m_{jj} (GeV)",30,300);
  //lumi_136TeV = "AlCaRaw DCSOnly, X fb^{-1}";
  lumi_136TeV = "AlCaRaw 1.79, Golden 4.86 fb^{-1}";
  TCanvas *c2 = tdrCanvas("c2",h2,8,11,kSquare);

  TLine *l = new TLine();
  l->DrawLine(30,0,300,0);
  l->SetLineColor(kRed);
  l->DrawLine(40,0,70,0);
  //l->DrawLine(105,0,300,0);
  l->DrawLine(110,0,300,0);

  hdf->Draw("SAME");
  //hd->SetYTitle("(Data-fit)/fit");
  //hd->SetXTitle("Dijet mass m_{jj} (GeV)");
  //hd->GetYaxis()->SetRangeUser(-0.05,0.05);
  //hd->GetYaxis()->SetMoreLogLabels();
  //hd->GetYaxis()->SetNoExponent();
  gPad->SetLogx();
  //gPad->SetLogy(kNone);
  
  //TLatex *t = new TLatex();
  //t->SetNDC();
  t->SetTextSize(0.035);
  t->DrawLatex(0.20,0.72,"All jet pairs");
  t->DrawLatex(0.20,0.67,"with |#Delta#eta|<1.3");
  t->SetTextSize(0.035);
  t->DrawLatex(0.20,0.27,"Fit range 40-300 GeV");
  t->DrawLatex(0.20,0.23,"Fit blinded at 70-110 GeV");
  t->DrawLatex(0.20,0.19,"Fit shape N#times(p_{T})^{logpol7}");
  t->DrawLatex(0.20,0.15,Form("Fit #chi^{2} / NDF = %1.1f / %d",
			      f1->GetChisquare(), f1->GetNDF()));

  c1->SaveAs("pdf/alcajme/alcajmemjj_mjj_v10.pdf");
  c2->SaveAs("pdf/alcajme/alcajmemjj_fit_v10.pdf");

  // Mass resolution estimated as sqrt(m/2.)
  TF1 *fg1 = new TF1("fg1","[0]*TMath::Gaus(x,[1]*80.4,[2]*6.34,0)",60,120);
  TF1 *fg2 = new TF1("fg2","[0]*TMath::Gaus(x,[1]*91.2,[2]*6.75,0)",60,120);
  TF1 *fwz = new TF1("fwz","[0]*TMath::Gaus(x,[1]*80.4,[2]*6.34,0)"
		     "+[3]*[0]*TMath::Gaus(x,[4]*[1]*91.2,[5]*[2]*6.75,0)",
		     60,120);
  fwz->SetParameters(0.025,1,1,1,1.,1);
  //fwz->SetParameters(0.026,1.024,0.82,0.034,1.019,0.99);
  //fwz->SetParameters(0.026,1.024,0.82,0.034,1.019,0.99);
  hdf->Fit(fwz,"RN");
  fwz->Draw("SAME");
  
  fg1->SetParameters(fwz->GetParameter(0),fwz->GetParameter(1),
		     fwz->GetParameter(2));
  fg1->SetLineColor(kGreen+2);
  fg1->Draw("SAME");

  fg2->SetParameters(fwz->GetParameter(0)*fwz->GetParameter(3),
		     fwz->GetParameter(1)*fwz->GetParameter(4),
		     fwz->GetParameter(2)*fwz->GetParameter(5));
  fg2->SetLineColor(kMagenta+2);
  fg2->Draw("SAME");
  


  t->SetTextColor(kRed);
  t->SetTextSize(0.030);
  t->DrawLatex(0.2,0.55,Form("Fit #chi^{2}/NDF = %1.1f/%d",
			     fwz->GetChisquare(), fwz->GetNDF()));
  t->SetTextColor(kGreen+2);
  t->DrawLatex(0.3,0.50,Form("N_{1} = %1.3g#pm%1.3g",
			     fwz->GetParameter(0),fwz->GetParError(0)));
  t->DrawLatex(0.3,0.47,Form("#mu_{1} = m_{W}#times%1.03f#pm%1.03f",
			      fwz->GetParameter(1),fwz->GetParError(1)));
  t->DrawLatex(0.3,0.43,Form("#sigma_{1} = #sqrt{m_{W}/2}#times%1.03f#pm%1.03f",
			      fwz->GetParameter(2),fwz->GetParError(2)));
  t->SetTextColor(kMagenta+2);
  t->DrawLatex(0.3,0.39,Form("N_{2} = N_{1}#times%1.3g#pm%1.3g",
			     fwz->GetParameter(3),fwz->GetParError(3)));
  t->DrawLatex(0.3,0.36,Form("#mu_{2} = #m_{Z}#timesm_{1}"
			     "#times%1.03f#pm%1.03f",
			     fwz->GetParameter(4),fwz->GetParError(4)));
  t->DrawLatex(0.3,0.32,Form("#sigma_{2} = #sqrt{m_{Z}/2}"
			     "#times#sigma_{1}#times%1.03f#pm%1.03f",
			     fwz->GetParameter(5),fwz->GetParError(5)));

  c2->SaveAs("pdf/alcajme/alcajme_fit2_v10.pdf");


  TH1D *h3 = tdrHist("h3","Data - fit",-2e3,2e4,
		     "Dijet mass m_{jj} (GeV)",30,300);
  //lumi_136TeV = "AlCaRaw DCSOnly, X fb^{-1}";
  lumi_136TeV = "AlCaRaw 1.79, Golden 4.86 fb^{-1}";
  TCanvas *c3 = tdrCanvas("c3",h3,8,11,kSquare);
  h3->GetYaxis()->SetNoExponent(kFALSE);
  
  hdfa->Draw("SAME");
  gPad->SetLogx();

  // Cross sections
  double xw = 20508.9; // pb, BR(W>lnu) included, l=mu
  double xz = 1928.0; // pb, Z->ll, l=mu
  double totxw = xw * 9.;
  double hadxw = totxw * 2./3.;
  double totxz = xz / 0.034;
  double hadxz = totxz * 0.692;
  double hadx = hadxw+hadxz;
  double lumi = 1.79e3 / 9277; // max estimate, could be half this
  double nw = hadxw * lumi;
  double rz = hadxz/hadxw;
  double nz = nw * rz;
  cout << "Total W xsec = " << totxw << " pb" << endl;
  cout << "Total Z xsec = " << totxz << " pb" << endl;
  cout << "Total hadronic W+Z xsec = " << hadx << " pb" << endl;
  cout << "Expected Z/W ratio = " << rz << endl;
  cout << "Expected Nw = " << nw << endl;
  cout << "Expected Nz = " << nz << endl;

  // Mass resolution estimated as sqrt(m/2.)
  TF1 *fga1 = new TF1("fga1",Form("[0]*%1.3g*TMath::Gaus(x,[1]*80.4,[2]*6.34,1)",nw),60,120);
  TF1 *fga2 = new TF1("fga2",Form("[0]*%1.3g*TMath::Gaus(x,[1]*91.2,[2]*6.75,1)",nz),60,120);
  TF1 *fwza = new TF1("fwza",Form("[0]*%1.3g*TMath::Gaus(x,[1]*80.4,[2]*6.34,1)"
				  "+[3]*[0]*%1.3g*TMath::Gaus(x,[4]*[1]*91.2,[5]*[2]*6.75,1)",nw,nz),
		      60,120);
  //fwza->SetParameters(5e5*0.025,1,1,1,1.,1);
  fwza->SetParameters(1,1,1,1,1.,1);

  hdfa->Fit(fwza,"RN");
  fwza->Draw("SAME");
  
  fga1->SetParameters(fwza->GetParameter(0),fwza->GetParameter(1),
		      fwza->GetParameter(2));
  fga1->SetLineColor(kGreen+2);
  fga1->Draw("SAME");

  fga2->SetParameters(fwza->GetParameter(0)*fwza->GetParameter(3),
		      fwza->GetParameter(1)*fwza->GetParameter(4),
		      fwza->GetParameter(2)*fwza->GetParameter(5));
  fga2->SetLineColor(kMagenta+2);
  fga2->Draw("SAME");
  
  t->SetTextColor(kBlack);
  t->SetTextSize(0.035);
  t->DrawLatex(0.20,0.72,"All jet pairs");
  t->DrawLatex(0.20,0.67,"with |#Delta#eta|<1.3");

  t->SetTextColor(kRed);
  t->SetTextSize(0.030);
  t->DrawLatex(0.6,0.55,Form("#chi^{2} / NDF = %1.1f / %d",
			     fwza->GetChisquare(), fwza->GetNDF()));
  t->SetTextColor(kGreen+2);
  t->DrawLatex(0.6,0.50,Form("N_{1} = N_{W,13 TeV}#times%1.3g#pm%1.3g",
			     fwza->GetParameter(0),fwza->GetParError(0)));
  t->DrawLatex(0.6,0.47,Form("#mu_{1} = m_{W}#times%1.03f#pm%1.03f",
			      fwza->GetParameter(1),fwza->GetParError(1)));
  t->DrawLatex(0.6,0.43,Form("#sigma_{1} = #sqrt{m_{W}/2}#times%1.03f#pm%1.03f",
			      fwza->GetParameter(2),fwza->GetParError(2)));
  t->SetTextColor(kMagenta+2);
  t->DrawLatex(0.6,0.39,Form("N_{2} = R_{Z/W}#timesn_{1}#times%1.3g#pm%1.3g",
			     fwza->GetParameter(3),fwza->GetParError(3)));
  t->DrawLatex(0.6,0.36,Form("#mu_{2} = m_{Z}#timesm_{1}"
			     "#times%1.03f#pm%1.03f",
			     fwza->GetParameter(4),fwza->GetParError(4)));
  t->DrawLatex(0.6,0.32,Form("#sigma_{2} = #sqrt{m_{Z}/2}"
			     "#timess_{1}#times%1.03f#pm%1.03f",
			     fwza->GetParameter(5),fwza->GetParError(5)));

  c3->SaveAs("pdf/alcajme/alcajme_fita_v10.pdf");
}
