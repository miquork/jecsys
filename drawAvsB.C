#include "TFile.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"

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

TGraphErrors* addGraph(TGraphErrors *g1, TGraphErrors *g2) {
  
  assert(g1->GetN()==g2->GetN());
  TGraphErrors *g = (TGraphErrors*)g1->Clone();//new TGraphErrors(g1->GetN());
  for (int i = 0; i != g1->GetN(); ++i) {

    double x1 = g1->GetX()[i];
    double x2 = g2->GetX()[i];
    assert(fabs(x1/x2-1)<0.1);

    double ex1 = g1->GetEX()[i];
    double ex2 = g2->GetEX()[i];

    double y1 = g1->GetY()[i];
    double y2 = g2->GetY()[i];

    double ey1 = g1->GetEY()[i];
    double ey2 = g2->GetEY()[i];

    double x = (ey2*ey2*x1 + ey1*ey1*x2) / (ey2*ey2+ey1*ey1);
    double y = (ey2*ey2*y1 + ey1*ey1*y2) / (ey2*ey2+ey1*ey1);
    double ex = (ey2*ey2*ex1 + ey1*ey1*ex2) / (ey2*ey2+ey1*ey1);
    double ey = (ey2*ey2*ey1 + ey1*ey1*ey2) / (ey2*ey2+ey1*ey1);

    g->SetPoint(i, x, y);
    g->SetPointError(i, ex, ey);
  }

  return g;
}

// Determine sensitivity to tracker dynamic inefficiency
// by studying ratio of jet responses in Runs G and F (and BCD / F, E / F)
void drawAvsB() {

  setTDRStyle();

  string epocha = "BCD";//"BCD";//"H";//"F";//"BCD";//"F";//"E";//"BCD";//"F";
  string epochb = "GH";//"G";//"BCD";//"G";//"E";//"E";//"F";//"G";

  // Add the rest as well
  string epocha2 = "";//"EF";
  string epochb2 = "";//"G";

  string type = "data";

  vector<string> methods;
  methods.push_back("mpfchs1");
  methods.push_back("ptchs");
  bool nozjptb = false;
  bool nogjmpf = false;
  bool nogjptb = true;
  bool mjvsjes = false;
  
  vector<string> samples;
  samples.push_back("zeejet");
  samples.push_back("zmmjet");
  samples.push_back("gamjet");
  //samples.push_back("multijet");

  cout << "draw"<<epocha<<"vs"<<epochb<<endl;
  const char *ct = type.c_str();
  const char *pa = epocha.c_str();
  const char *pb = epochb.c_str();

  const char *pa2 = epocha2.c_str();
  const char *pb2 = epochb2.c_str();

  TFile *fg = new TFile(Form("rootfiles/jecdata%s.root",pb),"READ");
  assert(fg && !fg->IsZombie());

  TFile *ff = new TFile(Form("rootfiles/jecdata%s.root",pa),"READ");
  assert(ff && !ff->IsZombie());

  TFile *fg2(0), *ff2(0);
  if (epochb2!="") fg2 = new TFile(Form("rootfiles/jecdata%s.root",pb2),"READ");
  if (epocha2!="") ff2 = new TFile(Form("rootfiles/jecdata%s.root",pa2),"READ");

  TH1D *h = new TH1D("h",
		     Form(";p_{T,ref} (GeV);%s ratio (%s / %s)",
			  (type=="ratio" ? "Data/MC" :
			   type=="data" ? "Data/data" : "MC/MC"),
			  (epocha + (epocha2!="" ? "+"+epocha2 : "")).c_str(),
			  (epochb + (epochb2!="" ? "+"+epochb2 : "")).c_str()),
		     3470,30,3500);
  h->SetMinimum(0.90);
  h->SetMaximum(1.15);
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();

  if (epocha=="F" && epochb=="G")
    lumi_13TeV = "Run2016F+G, 3.1+7.1 fb^{-1}";
  if (epocha=="BCD" && epochb=="G")
    lumi_13TeV = "Run2016BCD+H, 12.9+8.8 fb^{-1}";
  if (epocha=="BCD" && epochb=="G")
    lumi_13TeV = "Run2016BCD+FearlyGH, 12.9+16.8 fb^{-1}";
  if (epocha=="BCD" && epochb=="F")
    lumi_13TeV = "Run2016BCD+F, 13+3.1 fb^{-1}";
  if (epocha=="BCD" && epochb=="E")
    lumi_13TeV = "Run2016BCD+E, 13+4.0 fb^{-1}";
  if (epocha=="E" && epochb=="F")
    lumi_13TeV = "Run2016E+F, 4.0+3.1 fb^{-1}";
  if (epocha=="F" && epochb=="E")
    lumi_13TeV = "Run2016E+F, 4.0+3.1 fb^{-1}";

  if ((epocha=="BCDEF" && epochb=="GH") ||
      (epocha=="BCD" && epocha2=="EF" && epochb=="H" && epochb2=="G")) 
    lumi_13TeV = "Run2016BCDEF+GH, 19.7+16.8 fb^{-1}";
  if (epocha=="EF" && epochb=="BCD")
    lumi_13TeV = "Run2016BCD+EF, 12.9+6.8 fb^{-1}";
  if (epocha=="H" && epochb=="G")
    lumi_13TeV = "Run2016G+H, 8.0+8.8 fb^{-1}";

  if ((epocha=="BCD" && epocha2=="EF" && epochb=="G" && epochb2=="H")) 
    lumi_13TeV = "Run2016BCDFearly+FlateGH, 19.7+16.8 fb^{-1}";

  if ((epocha=="BCD" && epocha2=="" && ((epochb=="GH" && epochb2=="") ||
					(epochb=="G" && epochb2=="H"))))
    lumi_13TeV = "Run2016BCD+FlateGH, 12.9+16.8 fb^{-1}";
  if ((epocha=="EF" && epocha2=="" && ((epochb=="GH" && epochb2=="") ||
				       (epochb=="G" && epochb2=="H"))))
    lumi_13TeV = "Run2016EF+FlateGH, 6.8+16.8 fb^{-1}";

  if ((epocha=="EF" && epocha2=="" && epochb=="G" && epochb2=="H")) 
    lumi_13TeV = "Run2016EFearly+FlateGH, 6.8+16.8 fb^{-1}";


  TCanvas *c1 = tdrCanvas("c1",h,4,11,true);
  c1->SetLogx();

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);

  TMultiGraph *mg = new TMultiGraph();
  string s = "draw"+epocha+(epocha2!="" ? "p" + epocha2 : "")
    +"vs"+epochb+(epochb2!="" ? "p" + epochb2 : "");

  TGraphErrors *gmjb(0), *gmpf(0);

  for (unsigned int im = 0; im != methods.size(); ++im) {
    const char *cm = methods[im].c_str();

    tex->DrawLatex(0.20,0.75-0.06*im,cm);
    s += "_" + methods[im];

  for (unsigned int is = 0; is != samples.size(); ++is) {

    const char *cs = samples[is].c_str();
    TGraphErrors *gg = (TGraphErrors*)fg->Get(Form("%s/eta00-13/%s_%s_a30",ct,cm,cs));
    cout << cm << " " << cs << endl << flush;
    assert(gg);
    if (fg2) {
      TGraphErrors *gg2 = (TGraphErrors*)fg2->Get(Form("%s/eta00-13/%s_%s_a30",ct,cm,cs));
      assert(gg2);
      gg = addGraph(gg,gg2);
    }
    
    TGraphErrors *gf = (TGraphErrors*)ff->Get(Form("%s/eta00-13/%s_%s_a30",ct,cm,cs));
    assert(gf);
    if (ff2) {
      TGraphErrors *gf2 = (TGraphErrors*)ff2->Get(Form("%s/eta00-13/%s_%s_a30",ct,cm,cs));
      assert(gf2);
      gf = addGraph(gf,gf2);
    }
    
    if (!(gf->GetN()==gg->GetN())) {

      // Remove highest pT point is that is the offender (BCD vs GH)
      if (gg->GetN()>gf->GetN() &&
	  fabs(gg->GetX()[gg->GetN()-1]/gf->GetX()[gf->GetN()-1]-1)>0.1 &&
	  fabs(gg->GetX()[gg->GetN()-2]/gf->GetX()[gf->GetN()-1]-1)<0.1) {
	cout << "Remove point B(N-1)" << endl;
	gg->RemovePoint(gg->GetN()-1);
      }
      else {
	cout << "sample " << samples[is] << " method " << methods[im]
	     << " gf->N: " << gf->GetN() << " gg->N: " << gg->GetN() << endl;
	cout << " x_gf(N-1)=" << gf->GetX()[gf->GetN()-1]
	     << " x_gg(N-1)=" << gg->GetX()[gg->GetN()-1]
	     << " x_gg(N-2)=" << gg->GetX()[gg->GetN()-2] << endl;
      }

      assert(gf->GetN()==gg->GetN());
    }

    TGraphErrors *g = (TGraphErrors*)gg->Clone(Form("ge_%s_%s",cm,cs));
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
    else if ((samples[is]=="zmmjet" || samples[is]=="zeejet") &&
	     methods[im]=="ptchs" && nozjptb) {
      tex->SetTextColor(kRed);
      tex->DrawLatex(0.20,0.63,"Z+jet p_{T}^{bal} excl. from fit");
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
  if (nozjptb) s += "_nozptb";
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
  
  //TF1 *ft = new TF1("ft","1-[0]-[1]*pow(x,[2]) + ([3]+[4]*log(x))/x",30,2200);
  //ft->SetParameters(0,0.05,-0.5,1,0.1);
  //ft->FixParameter(3,0);

  // Logarithmic sigmoid
  //TF1 *ft = new TF1("ft","[0]+(1-[0])/(1. + exp(-(log(x)-log(abs([1])))"
  //	       "/(log(abs([2])+abs([1]))-log(abs([1])))))", 30,2200);
  //ft->SetParameters(0.98, 150, 50);
  TF1 *ft = new TF1("ft","[0]+(1-[0])/(1. + exp(-(log(x)-[1])/[2]))",30,2200);
  //ft->SetParameters(0.98,log(145),log(190)-log(145));
  //ft->SetParameters(0.982,4.967,0.271);
  //ft->SetParameters(0.976,5.040,0.370); // ENDCAP
  //ft->SetParameters(0.985,5.0,0.3);
  ft->SetParameters(0.985,5.025,0.3);
  //ft->FixParameter(1,5.03); // semi-weighted average of BCD and EF
  //ft->FixParameter(2,0.395); // combined fit to BCD+EF / G+H 

  // ( 12.9*5.055+6.8*5.000)/(12.9+6.8)
  ft->FixParameter(1,5.036); // semi-weighted average of BCD/GH and EF/GH
  // ( 12.9*0.344 + 6.8*0.455)/(12.9+6.8)
  ft->FixParameter(2,0.391); // combined fit to BCD+EF / GH 

  // Log-sigmoid + powerlaw
  //TF1 *ft = new TF1("ft","[0]+(1-[0])/(1. + exp(-(log(x)-[1])/[2]))"
  //	       "*(1-[3]*pow(x,[4]))",30,2200);
  //ft->SetParameters(0.982,4.967,0.271,0.1,-0.2);
  // Double powerlaw
  //TF1 *ft = new TF1("ft","[4]-[0]*pow(x,[1])-[2]*pow(x,[3])",30,2200);
  //ft->SetParameters(0.05,-0.15,0.01,-0.3,1);
  

  mg->Fit(ft,"RN");
  ft->SetLineColor(kBlue);
  ft->SetLineWidth(2);
  ft->SetRange(10.,3500.);
  ft->Draw("SAME");

  // Map multijet with response ratio
  if (gmpf) { // we have multijet available
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
  } // multijet

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
  if (ft->GetNpar()>3)
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
