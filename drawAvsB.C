#include "tdrstyle_mod15.C"

TF1 *fhb(0);
Double_t jesFit(Double_t *x, Double_t *p) {

  double pt = *x;

  // Initialize SinglePionHCAL and PileUpPt shapes
  if (!fhb) fhb = new TF1("fhb","max(0.,[0]+[1]*pow(x,[2]))",10,3500);
  fhb->SetParameters(1.03091e+00, -5.11540e-02, -1.54227e-01); // SPRH
 
  // p[0]: overall scale shift, p[1]: HCAL shift in % (full band +3%)
  double ptref = 208; // pT that minimizes correlation in p[0] and p[1]
  
  return (p[0] + p[1]/3.*100*(fhb->Eval(pt)-fhb->Eval(ptref)));
}

// Determine sensitivity to tracker dynamic inefficiency
// by studying ratio of jet responses in Runs G and F (and BCD / F, E / F)
void drawAvsB() {

  setTDRStyle();

  string epocha = "H";//"H";//"F";//"BCD";//"F";//"E";//"BCD";//"F";
  string epochb = "G";//"G";//"BCD";//"G";//"E";//"E";//"F";//"G";

  string type = "data";

  vector<string> methods;
  methods.push_back("mpfchs1");
  methods.push_back("ptchs");
  bool nogjmpf = false;
  bool nogjptb = false;//true;
  bool mjvsjes = false;
  
  vector<string> samples;
  samples.push_back("zeejet");
  samples.push_back("zmmjet");
  samples.push_back("gamjet");
  samples.push_back("multijet");

  cout << "draw"<<epocha<<"vs"<<epochb<<endl;
  const char *ct = type.c_str();
  const char *pa = epocha.c_str();
  const char *pb = epochb.c_str();

  TFile *fg = new TFile(Form("rootfiles/jecdata%s.root",pb),"READ");
  assert(fg && !fg->IsZombie());

  TFile *ff = new TFile(Form("rootfiles/jecdata%s.root",pa),"READ");
  assert(ff && !ff->IsZombie());

  TH1D *h = new TH1D("h",
		     Form(";p_{T,ref} (GeV);%s ratio (%s / %s)",
			  (type=="ratio" ? "Data/MC" :
			   type=="data" ? "Data/data" : "MC/MC"),
			  epocha.c_str(),epochb.c_str()),
		     //870,30,900);
		     //2170,30,2200);
		     3470,30,3500);
  h->SetMinimum(0.90);//epocha=="BCD" ? 0.95 : 0.92);
  h->SetMaximum(1.15);//epocha=="BCD" ? 1.11 : 1.08);
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();

  if (epocha=="F" && epochb=="G")
    lumi_13TeV = "Run2016F+G, 3.1+7.1 fb^{-1}";
  if (epocha=="BCD" && epochb=="G")
    //lumi_13TeV = "Run2016BCD+G, 13+7.1 fb^{-1}";
    lumi_13TeV = "Run2016BCD+FearlyGH, 12.9+16.8 fb^{-1}";
  if (epocha=="BCD" && epochb=="F")
    lumi_13TeV = "Run2016BCD+F, 13+3.1 fb^{-1}";
  if (epocha=="BCD" && epochb=="E")
    lumi_13TeV = "Run2016BCD+E, 13+4.0 fb^{-1}";
  if (epocha=="E" && epochb=="F")
    lumi_13TeV = "Run2016E+F, 4.0+3.1 fb^{-1}";
  if (epocha=="F" && epochb=="E")
    lumi_13TeV = "Run2016E+F, 4.0+3.1 fb^{-1}";

  if (epocha=="BCDEF" && epochb=="GH")
    lumi_13TeV = "Run2016BCDEF+GH, 19.7+16.8 fb^{-1}";
  if (epocha=="EF" && epochb=="BCD")
    lumi_13TeV = "Run2016BCD+EF, 12.9+6.8 fb^{-1}";
  if (epocha=="H" && epochb=="G")
    lumi_13TeV = "Run2016G+H, 8.0+8.8 fb^{-1}";

  TCanvas *c1 = tdrCanvas("c1",h,4,11,true);
  c1->SetLogx();

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);

  TMultiGraph *mg = new TMultiGraph();
  string s = "draw"+epocha+"vs"+epochb;

  TGraphErrors *gmjb(0), *gmpf(0);

  for (unsigned int im = 0; im != methods.size(); ++im) {
    const char *cm = methods[im].c_str();

    tex->DrawLatex(0.20,0.75-0.06*im,cm);
    s += "_" + methods[im];

  for (unsigned int is = 0; is != samples.size(); ++is) {

    //if (samples[is]=="gamjet" && methods[im]=="mpfchs1" && nogjmpf) continue;
    
    const char *cs = samples[is].c_str();
    TGraphErrors *gg = (TGraphErrors*)fg->Get(Form("%s/eta00-13/%s_%s_a30",ct,cm,cs));
    cout << cm << " " << cs << endl << flush;
    assert(gg);
    
    TGraphErrors *gf = (TGraphErrors*)ff->Get(Form("%s/eta00-13/%s_%s_a30",ct,cm,cs));
    assert(gf);
    
    //gg->Draw("AP");
    //gf->Draw("SAMEP");

    TGraphErrors *g = (TGraphErrors*)gg->Clone(Form("ge_%s_%s",cm,cs));
    if (!(gf->GetN()==gg->GetN())) {
      cout << "sample " << samples[is] << " method " << methods[im]
	   << " gf->N: " << gf->GetN() << " gg->N: " << gg->GetN() << endl;
      assert(gf->GetN()==gg->GetN());
    }
    for (int i = 0; i != g->GetN(); ++i) {
      double yg = gg->GetY()[i];
      double yf = gf->GetY()[i];
      g->SetPoint(i, gg->GetX()[i], yf / yg);
      double ex = gg->GetEX()[i];
      double eg = gg->GetEY()[i];
      double ef = gf->GetEY()[i];
      g->SetPointError(i, ex, yf/yg*sqrt(pow(eg/yg,2)+pow(ef/yf,2)));
    }
    //g->Draw(is==0 ? "AP" : "SAMEP");
    g->SetLineWidth(1+is);
    g->Draw("SAMEPZ");

    if (samples[is]=="gamjet" && methods[im]=="mpfchs1" && nogjmpf) {
      tex->SetTextColor(kBlue);
      tex->DrawLatex(0.20,0.63,"#gamma+jet MPF excl. from fit");
      tex->SetTextColor(kBlack);
    }
    else if (samples[is]=="gamjet" && methods[im]=="ptchs" && nogjptb) {
      tex->SetTextColor(kBlue);
      tex->DrawLatex(0.20,0.63,"#gamma+jet p_{T}^{bal} excl. from fit");
      tex->SetTextColor(kBlack);
    }
    else if (samples[is]=="multijet") {
      g->SetMarkerColor(kGray+1);
      g->SetLineColor(kGray+1);
      if (methods[im]=="ptchs") gmjb = g;
      if (methods[im]=="mpfchs1") gmpf = g;
    }
    else
      mg->Add(g);
  } // for is
  } // for im
  
  if (nogjmpf) s += "_nogjmpf";
  if (nogjptb) s += "_nogptb";
  if (mjvsjes) {
    s += "_mjvsjes";
    tex->SetTextColor(kBlack);
    tex->DrawLatex(0.20,0.58,"Multijet vs JES fit");
  }

  TF1 *fjes = new TF1("fjes",jesFit,30,2200,2);
  fjes->SetParameters(0.99,0.05);
  mg->Fit(fjes,"RN");
  fjes->SetLineColor(kBlack);
  fjes->SetLineStyle(kDashed);
  fjes->SetLineWidth(2);
  fjes->SetRange(10.,3500.);
  fjes->Draw("SAME");
  
  //TF1 *ft = new TF1("ft","0.85-[0] + "
  //		    "(0.15+[0])*([3]-[1]*pow(x,[2])) + [4]/x",30,1300);
  //ft->SetParameters(0.05,1,0.8-1,1,0.5);
  // Standard fit 58.7/60
  TF1 *ft = new TF1("ft","1-[0]-[1]*pow(x,[2]) + ([3]+[4]*log(x))/x",30,2200);
  ft->SetParameters(0,0.05,-0.5,1,0.1); // 58.8/60
  //ft->FixParameter(0,0); // 59.4/61
  //ft->FixParameter(4,0); // 59.5/62
  ft->FixParameter(3,0);

  //ft->FixParameter(3,1);
  mg->Fit(ft,"RN");
  ft->SetLineColor(kBlue);
  ft->SetLineWidth(2);
  ft->SetRange(10.,3500.);
  ft->Draw("SAME");

  // Map multijet with response ratio
  TGraphErrors *gmpf2 = (TGraphErrors*)gmpf->Clone("gmpf2");
  gmpf2->SetMarkerColor(kBlack);//kGray+1);
  gmpf2->SetLineColor(kBlack);//kGray+1);
  for (int i = 0; i != gmpf->GetN(); ++i) {
    if (mjvsjes) {
      gmpf2->SetPoint(i, 0.4*gmpf->GetX()[i],
		      fjes->Eval(gmpf->GetX()[i])/gmpf->GetY()[i]);
      gmpf2->SetPointError(i, 0.4*gmpf->GetEX()[i],
			   gmpf->GetEY()[i]);
    }
    else {
      gmpf2->SetPoint(i, 0.4*gmpf->GetX()[i],
		      ft->Eval(gmpf->GetX()[i])/gmpf->GetY()[i]);
      gmpf2->SetPointError(i, 0.4*gmpf->GetEX()[i],
			   gmpf->GetEY()[i]);
    }
  }
  gmpf2->Draw("SAMEPz");

  tex->SetTextColor(kBlue);
  tex->DrawLatex(0.50,0.85,Form("#chi^{2} / NDF = %1.1f / %d",
				ft->GetChisquare(),
				ft->GetNDF()));
  tex->SetTextColor(kBlack);
  tex->SetTextSize(0.040);
  tex->DrawLatex(0.50,0.80,Form("(#chi^{2} / NDF = %1.1f / %d)",
				fjes->GetChisquare(),
				fjes->GetNDF()));


  tex->SetTextColor(kBlue-9);
  tex->SetTextSize(0.030);
  tex->DrawLatex(0.20,0.25,ft->GetExpFormula());
  tex->DrawLatex(0.20,0.20,
		 Form("p_{0}=%1.3f#pm%1.3f"
		      ", p_{1}=%1.3f#pm%1.3f"
		      ", p_{2}=%1.3f#pm%1.3f",
		      ft->GetParameter(0),ft->GetParError(0),
		      ft->GetParameter(1),ft->GetParError(1),
		      ft->GetParameter(2),ft->GetParError(2)));
  tex->DrawLatex(0.20,0.17,
		 Form("p_{3}=%1.3f#pm%1.3f"
		      ", p_{4}=%1.3f#pm%1.3f",
		      ft->GetParameter(3),ft->GetParError(3),
		      ft->GetParameter(4),ft->GetParError(4)));

  c1->SaveAs(Form("pdf/%s.pdf",s.c_str()));

  for (int i = 0; i != ft->GetNpar(); ++i) {
    cout << Form("%s%1.4g",i==0 ? "{" : ", ",ft->GetParameter(i));
  }
  cout << "}" << endl;
    

}
