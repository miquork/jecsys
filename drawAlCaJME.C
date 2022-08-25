// Purpose: produce plots of AlCaJME results
#include "TFile.h"
#include "TProfile.h"
#include "TLine.h"
#include "TH2D.h"
#include "TF1.h"
#include "TString.h"

#include "tdrstyle_mod22.C"
#include <map>
#include <string>

TF1 *F(0), *f1(0);
Double_t _F(Double_t *x, Double_t *p) {
  double pt = x[0];
  double eta = p[0];
  double sqrts = 13600.;
  double alpha = -5;
  //double beta = 10.-2.*fabs(eta); // lhcewwg_jets_2018_06_13.pdf
  double beta = 10.-1.*fabs(eta); // 
  return (1e11*pow(pt,alpha)*pow(1.-2.*pt*cosh(eta)/sqrts,beta));
} // _F

TH1D *fixB(TProfile* p) {
  TH1D *h = p->ProjectionX(Form("h_%s",p->GetName()));
  if (!f1) f1 = new TF1("f1","[0]",-1.3,+1.3);
  if (!F) F = new TF1("F",_F,5,0.5*13600,1);
  h->Fit(f1,"QRN");
  double p0 = f1->GetParameter(0); // MPF
  if (p0>0.35) cout << p->GetName() << ": " << p0 << endl;
  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    if (p0>0.35) {
      // Stage 1: constant spectrum and JER
      //h->SetBinContent(i, h->GetBinContent(i)/p0);

      // Stage 2: add eta-dependent spectrum
      // voutilainen_thesis.pdf p.259 Eq.(C.8):
      // <pT,ptcl> = pT - alpha*sigma^2 for spectrum N0*exp(-alpha*pT)
      // Assume sigma constant (holds for tag at least), but alpha variable
      // Locally, alpha = dF/dpT * 1/F
      // #1 MPF = 1+(pT,probe-pT,tag)/pT,tag = pT,ptcl/(pT,ptcl+a*s^2)
      // => MPF*pT+MPF*as2 = pT => as2 = (1-MPF)/MPF * pT
      // #1b MPF = 1+2*(pT,probe-pT,tag)/(pTprobe+pT,tag)
      // =>  MPF = 1-2*(a*s^2)/(2*pT,ptcl+a*s^2)
      // => 2as2 = (1-MPF)*2*pT +(1-MPF)*as2
      // => as2 = 2*(1-MPF)/(1+MPF) * pT
      //
      // #2 MPF = 1+(pT,probe-pT,tag)/pT,probe = 2-pT,ptcl/(pT,ptcl+a*s^2)
      // => (2-MPF)*pT+(2-MPF)*as2 = pT => as2 = (1-(2-MPF))/(2-MPF) * pT
      // => as2 = (MPF-1)/(2-MPF) * pT
      //
      // #3 MPF = 1+(pT,probe-pT,tag)/pTave = 1+pTdiff/pTave = 1+dpT/(pT+as^2/2)
      // => (MPF-1)*pT+(MPF-1)*as2/2 = dpT
      // => as2 = (d-(MPF-1))/(MPF-1)*2 * pT
      // (just that MPF=1 for barrel so unsolvable => use as2 from tag, probe)
      double pt = 40.;
      double eta = h->GetBinCenter(i);
      F->SetParameter(0,0.);
      double df0 = F->Derivative(pt) / F->Eval(pt);
      F->SetParameter(0,eta);
      double df = F->Derivative(pt) / F->Eval(pt);
      //double p0new = 1 + (p0-1) * df / df0;
      //h->SetBinContent(i, h->GetBinContent(i)/p0new);
      if (p0<0.97) { // pttag
	double as2ref = (1-p0)/p0*pt;
	double as2new = as2ref * df / df0;
	double p0new = pt/(pt+as2new);
	h->SetBinContent(i, h->GetBinContent(i)/p0new);

	// Post-processing for central value
	// MPF = 1 + (pTp-pTt)/pTt
	// pTp=c*pTt => MPF = 1 + (c-1)/1 = c
	// => MPF=c and no post-processing needed
      }
      else if (p0>1.03) { // ptprobe
	double as2ref = (p0-1)/(2-p0)*pt;
	double as2new = as2ref * df / df0;
	double p0new = 2-pt/(pt+as2new);
	double mpf = h->GetBinContent(i)/p0new;
	h->SetBinContent(i, mpf);

	// Post-processing for central value
	// MPF = 1 + (pTp-pTt)/pTp
	// pTp=c*pTt => MPF = 1 + (c-1)/c
	// => (MPF-1)*c = (c-1) => (MPF-2)*c = -1 => c = 1/(2-MPF)
	double c = 1./(2.-mpf);
	h->SetBinContent(i, c);
      }
      else { // ptave
	double as2ref = 0.5*( (1.07-1.)/(2.-1.07)*pt + (1.-0.93)/0.93*pt);
	double as2new =  as2ref * df / df0;
	double d = (h->GetBinContent(i)-1)*(pt+as2new/2.)/pt;
	double mpf = 1+d;
	h->SetBinContent(i, mpf);
	
	// Post-processing for central value
	// MPF = 1 + (pTp-pTt)/((pTp+pTt)/2)
	// pTp=c*pTt => MPF = 1 + 2*(c-1)/(c+1)
	// => MPF = (3*c-1)/(c+1) => MPF*c+MPF=3*c-1
	// => (MPF-3)*c = -(1+MPF) => c = (1+MPF)/(3-MPF)
	double c = (1.+mpf)/(3.-mpf);
	h->SetBinContent(i, c);
      }
    } // if p0>0.35
    else {
      h->SetBinContent(i, h->GetBinContent(i)-p0);
    }
  } // for i

  return h;
} // fixB

// difference band
TH1D *diffB(TH1D *h1, TH1D *h2) {
  assert(h1->GetNbinsX()==h2->GetNbinsX());
  TH1D *hd = (TH1D*)h1->Clone(Form("hd_%s_%s",h1->GetName(),h2->GetName()));
  for (int i = 1; i != hd->GetNbinsX()+1; ++i) {
    hd->SetBinError(i, 0.5*fabs(h2->GetBinContent(i) - h1->GetBinContent(i)));
    hd->SetBinContent(i, 0.5*(h2->GetBinContent(i) + h1->GetBinContent(i)));
  }
  return hd;
} // diffB

void drawJER(const char *cs, const char *cs2);
void drawComp(const char *cs);

void drawAlCaJME(string set) {

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  const char *cs = set.c_str();
  TString s = cs;
  cout << "Running over set " << cs << endl << flush;
  cout << "(set="<<set<<", s="<<s<<")"<<endl<<flush;

  TFile *f(0);
  if (s.Contains("jmenanodt")) {
    //f = new TFile("../dijet/rootfiles/jmenano_data_out_v2.root","READ");
    //f = new TFile("../dijet/rootfiles/jmenano_data_out_v2_fwd40.root","READ");
    f = new TFile("../dijet/rootfiles/jmenano_data_out_v5_4p86fb.root","READ");
    if (f) { f->cd("HLT_DiPFJetAve40/Pt_40_50"); f = (TFile*)gDirectory; }
  }
  else if (s.Contains("jmenanomc"))
    f = new TFile("../dijet/rootfiles/jmenano_mc_out_v2.root","READ");
  else {
    //f = new TFile("rootfiles/alcajme_out_50M_v5_denom.root","READ");
    //f = new TFile("rootfiles/alcajme_out_50M_v8_doak4.root","READ");
    f = new TFile("rootfiles/alcajme_out_50M_v10_4p86fb.root","READ");
    if (f) { f->cd("Pt_40_50"); f = (TFile*)gDirectory; }
  }
  assert(f && !f->IsZombie());
  curdir->cd();

  TProfile *p(0);
  map<string, TProfile*> mp;
  mp["m0a"] = p = (TProfile*)f->Get("pm0ab"); assert(p);
  mp["m0t"] = p = (TProfile*)f->Get("pm0tc"); assert(p);
  mp["m0p"] = p = (TProfile*)f->Get("pm0pf"); //assert(p);
  if (!p) { mp["m0p"] = p = (TProfile*)f->Get("pm0pb"); assert(p); }

  mp["m2a"] = p = (TProfile*)f->Get("pm2ab"); assert(p);
  mp["m2t"] = p = (TProfile*)f->Get("pm2tc"); assert(p);
  mp["m2p"] = p = (TProfile*)f->Get("pm2pf"); assert(p);

  mp["mna"] = p = (TProfile*)f->Get("pmnab"); assert(p);
  mp["mnt"] = p = (TProfile*)f->Get("pmntc"); assert(p);
  mp["mnp"] = p = (TProfile*)f->Get("pmnpf"); assert(p);

  mp["mua"] = p = (TProfile*)f->Get("pmuab"); assert(p);
  mp["mut"] = p = (TProfile*)f->Get("pmutc"); assert(p);
  mp["mup"] = p = (TProfile*)f->Get("pmupf"); assert(p);

  mp["moa"] = p = (TProfile*)f->Get("pmoab"); assert(p);
  mp["mot"] = p = (TProfile*)f->Get("pmotc"); assert(p);
  mp["mop"] = p = (TProfile*)f->Get("pmopf"); assert(p);


  if (true) {
    //lumi_136TeV = "AlCaRaw DCSonly, X fb^{-1}";
    lumi_136TeV = "AlCaRaw 1.79, Golden 4.86 fb^{-1}";
    if (s.Contains("jmenanodt")) {
      //lumi_136TeV = "JMENANO DCSonly, X fb^{-1}";
      lumi_136TeV = "JMENANO 0.450, Golden 4.86 fb^{-1}";
    }
    if (s.Contains("jmenanomc")) {
      lumi_136TeV = "JMENANO FlatQCD";
    }
    extraText = "Preliminary";
    TH1D *h = tdrHist("h1","MPF",-0.30,1.60,"#eta",-5.191,5.191);
    TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);

    TLine *l = new TLine();
    l->DrawLine(-5.191,1,5.191,1);
    l->DrawLine(-5.191,0,5.191,0);

    tdrDraw(mp["m2a"],"H",kNone,kGreen+2,kDotted,-1,kNone);
    tdrDraw(mp["m2t"],"H",kNone,kBlue,kDotted,-1,kNone);
    tdrDraw(mp["m2p"],"H",kNone,kRed,kDotted,-1,kNone);
    
    tdrDraw(mp["mna"],"H",kNone,kGreen+2,kDotted,-1,kNone);
    tdrDraw(mp["mnt"],"H",kNone,kBlue,kDotted,-1,kNone);
    tdrDraw(mp["mnp"],"H",kNone,kRed,kDotted,-1,kNone);
    
    tdrDraw(mp["mua"],"H",kNone,kGreen+2,kDashed,-1,kNone);
    tdrDraw(mp["mut"],"H",kNone,kBlue,kDashed,-1,kNone);
    tdrDraw(mp["mup"],"H",kNone,kRed,kDashed,-1,kNone);

    tdrDraw(mp["moa"],"H",kNone,kGreen+2,kSolid,-1,kNone);
    tdrDraw(mp["mot"],"H",kNone,kBlue,kSolid,-1,kNone);
    tdrDraw(mp["mop"],"H",kNone,kRed,kSolid,-1,kNone);
    
    tdrDraw(mp["m0a"],"H",kNone,kGreen+2,kSolid,-1,kNone);
    tdrDraw(mp["m0t"],"H",kNone,kBlue,kSolid,-1,kNone);
    tdrDraw(mp["m0p"],"H",kNone,kRed,kSolid,-1,kNone);
    
    TLegend *t1a = tdrLeg(0.40,0.90-2*0.045,0.65,0.90);
    t1a->AddEntry(mp["m2p"],"DB","L");
    t1a->AddEntry(mp["m0p"],"MPF","L");
    TLegend *t1b = tdrLeg(0.57,0.90-3*0.045,0.82,0.90);
    t1b->AddEntry(mp["m0p"],"PtProbe","L");
    t1b->AddEntry(mp["m0a"],"PtAve","L");
    t1b->AddEntry(mp["m0t"],"PtTag","L");

    TLegend *t1c = tdrLeg(0.40,0.32,0.65,0.32+3*0.045);
    t1c->AddEntry(mp["mnt"],"n-jet","L");
    t1c->AddEntry(mp["mut"],"uncl.","L");
    t1c->AddEntry(mp["mot"],"both","L");
    TLegend *t1d = tdrLeg(0.57,0.32,0.82,0.32+3*0.045);
    t1d->AddEntry(mp["mut"],"PtTag","L");
    t1d->AddEntry(mp["mua"],"PtAve","L");
    t1d->AddEntry(mp["mup"],"PtProbe","L");

    gPad->RedrawAxis();
    c1->SaveAs(Form("pdf/alcajme/drawAlCaJME_loose_%s.pdf",cs));
  } // c1


  if (true) {
    //lumi_136TeV = "AlCaRaw DCSonly, X fb^{-1}";
    lumi_136TeV = "AlCaRaw 1.79, Golden 4.86 fb^{-1}";
    if (s.Contains("jmenanodt")) {
      //lumi_136TeV = "JMENANO DCSonly, X fb^{-1}";
      lumi_136TeV = "JMENANO 0.450, Golden 4.86 fb^{-1}";
    }
    if (s.Contains("jmenanomc")) {
      lumi_136TeV = "JMENANO FlatQCD";//, X fb^{-1}";
    }
    extraText = "Preliminary";
    TH1D *h = tdrHist("h2","MPF",-0.30,1.60,"#eta",-5.191,5.191);
    TCanvas *c1 = tdrCanvas("c2",h,8,11,kSquare);

    TLine *l = new TLine();
    l->DrawLine(-5.191,1,5.191,1);
    l->DrawLine(-5.191,0,5.191,0);
    
    tdrDraw(fixB(mp["m2a"]),"H",kNone,kGreen+2,kDotted,-1,kNone);
    tdrDraw(fixB(mp["m2t"]),"H",kNone,kBlue,kDotted,-1,kNone);
    tdrDraw(fixB(mp["m2p"]),"H",kNone,kRed,kDotted,-1,kNone);
    
    tdrDraw(fixB(mp["mna"]),"H",kNone,kGreen+2,kDotted,-1,kNone);
    tdrDraw(fixB(mp["mnt"]),"H",kNone,kBlue,kDotted,-1,kNone);
    tdrDraw(fixB(mp["mnp"]),"H",kNone,kRed,kDotted,-1,kNone);
    
    tdrDraw(fixB(mp["mua"]),"H",kNone,kGreen+2,kDotted,-1,kNone);
    tdrDraw(fixB(mp["mut"]),"H",kNone,kBlue,kDotted,-1,kNone);
    tdrDraw(fixB(mp["mup"]),"H",kNone,kRed,kDotted,-1,kNone);

    tdrDraw(fixB(mp["moa"]),"H",kNone,kGreen+2,kSolid,-1,kNone);
    tdrDraw(fixB(mp["mot"]),"H",kNone,kBlue,kSolid,-1,kNone);
    tdrDraw(fixB(mp["mop"]),"H",kNone,kRed,kSolid,-1,kNone);
    
    tdrDraw(fixB(mp["m0a"]),"H",kNone,kGreen+2,kSolid,-1,kNone);
    tdrDraw(fixB(mp["m0t"]),"H",kNone,kBlue,kSolid,-1,kNone);
    tdrDraw(fixB(mp["m0p"]),"H",kNone,kRed,kSolid,-1,kNone);
    
    TLegend *t1a = tdrLeg(0.40,0.90-2*0.045,0.65,0.90);
    t1a->AddEntry(mp["m2p"],"DB","L");
    t1a->AddEntry(mp["m0p"],"MPF","L");
    TLegend *t1b = tdrLeg(0.57,0.90-3*0.045,0.82,0.90);
    t1b->AddEntry(mp["m0p"],"PtProbe","L");
    t1b->AddEntry(mp["m0a"],"PtAve","L");
    t1b->AddEntry(mp["m0t"],"PtTag","L");

    TLegend *t1c = tdrLeg(0.40,0.32,0.65,0.32+3*0.045);
    t1c->AddEntry(mp["mnt"],"n-jet","L");
    t1c->AddEntry(mp["mut"],"uncl.","L");
    t1c->AddEntry(mp["mot"],"both","L");
    TLegend *t1d = tdrLeg(0.57,0.32,0.82,0.32+3*0.045);
    t1d->AddEntry(mp["mut"],"PtTag","L");
    t1d->AddEntry(mp["mua"],"PtAve","L");
    t1d->AddEntry(mp["mup"],"PtProbe","L");

    gPad->RedrawAxis();
    c1->SaveAs(Form("pdf/alcajme/drawAlCaJME_tight_%s.pdf",cs));
  } // c2

  if (true) {
    //lumi_136TeV = "AlCaRaw DCSonly, X fb^{-1}";
    lumi_136TeV = "AlCaRaw 1.79, Golden 4.86 fb^{-1}";
    if (s.Contains("jmenanodt")) {
      //lumi_136TeV = "JMENANO DCSonly, X fb^{-1}";
      lumi_136TeV = "JMENANO 0.450, Golden 4.86 fb^{-1}";
    }
    if (s.Contains("jmenanomc")) {
      lumi_136TeV = "JMENANO FlatQCD";//, X fb^{-1}";
    }
    extraText = "Preliminary";
    TH1D *h = tdrHist("h3","MPF",0.65,1.20,"#eta",-5.191,5.191);
    TCanvas *c1 = tdrCanvas("c3",h,8,11,kSquare);

    TLine *l = new TLine();
    l->DrawLine(-5.191,1,5.191,1);
    l->DrawLine(-5.191,0,5.191,0);
    
    tdrDraw(fixB(mp["m2a"]),"H",kNone,kGreen+2,kDotted,-1,kNone);
    tdrDraw(fixB(mp["m2t"]),"H",kNone,kBlue,kDotted,-1,kNone);
    tdrDraw(fixB(mp["m2p"]),"H",kNone,kRed,kDotted,-1,kNone);
    
    tdrDraw(fixB(mp["m0a"]),"H",kNone,kGreen+2,kSolid,-1,kNone);
    tdrDraw(fixB(mp["m0t"]),"H",kNone,kBlue,kSolid,-1,kNone);
    tdrDraw(fixB(mp["m0p"]),"H",kNone,kRed,kSolid,-1,kNone);
    
    TLegend *t1a = tdrLeg(0.40,0.90-2*0.045,0.65,0.90);
    t1a->AddEntry(mp["m2p"],"DB","L");
    t1a->AddEntry(mp["m0p"],"MPF","L");
    TLegend *t1b = tdrLeg(0.57,0.90-3*0.045,0.82,0.90);
    t1b->AddEntry(mp["m0p"],"PtProbe","L");
    t1b->AddEntry(mp["m0a"],"PtAve","L");
    t1b->AddEntry(mp["m0t"],"PtTag","L");

    gPad->RedrawAxis();
    c1->SaveAs(Form("pdf/alcajme/drawAlCaJME_zoom_%s.pdf",cs));
  } // c3

    drawJER("a",cs);
    drawJER("t",cs);
    drawJER("p",cs);

    drawComp(cs);
}


void drawJER(const char *cs, const char *cs2) {

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  TString s(cs2);
  
  TFile *f(0);
  if (s.Contains("jmenanodt")) {
    //f = new TFile("../dijet/rootfiles/jmenano_data_out_v2.root","READ");
    //f = new TFile("../dijet/rootfiles/jmenano_data_out_v2_fwd40.root","READ");
    f = new TFile("../dijet/rootfiles/jmenano_data_out_v5_4p86db.root","READ");
    if (f) { f->cd("HLT_DiPFJetAve40/Pt_40_50"); f = (TFile*)gDirectory; }
  }
  else if (s.Contains("jmenanomc"))
    f = new TFile("../dijet/rootfiles/jmenano_mc_out_v2.root","READ");
  //else
    //f = new TFile("rootfiles/alcajme_out_50M_v5_denom.root","READ");
  else {
    //f = new TFile("rootfiles/alcajme_out_50M_v8_doak4.root","READ");
    f = new TFile("rootfiles/alcajme_out_50M_v10_4p86fb.root","READ");
    if (f) { f->cd("Pt_40_50"); f = (TFile*)gDirectory; }
  }
  assert(f && !f->IsZombie());
  curdir->cd();

  /*
  TH2D *h2m = (TH2D*)f->Get("h2m0aeta"); assert(h2m);
  TH2D *h2x = (TH2D*)f->Get("h2m0xaeta"); assert(h2x);
  TH2D *h2a = (TH2D*)f->Get("h2m2aeta"); assert(h2a);
  TH2D *h2ax = (TH2D*)f->Get("h2m2xaeta"); assert(h2ax);
  */
  /*
  TH2D *h2m = (TH2D*)f->Get("h2m0peta"); assert(h2m);
  TH2D *h2x = (TH2D*)f->Get("h2m0xpeta"); assert(h2x);
  TH2D *h2a = (TH2D*)f->Get("h2m2peta"); assert(h2a);
  TH2D *h2ax = (TH2D*)f->Get("h2m2xpeta"); assert(h2ax);
  */
  /*
  TH2D *h2m = (TH2D*)f->Get("h2m0teta"); assert(h2m);
  TH2D *h2x = (TH2D*)f->Get("h2m0xteta"); assert(h2x);
  TH2D *h2a = (TH2D*)f->Get("h2m2teta"); assert(h2a);
  TH2D *h2ax = (TH2D*)f->Get("h2m2xteta"); assert(h2ax);
  */
  TH2D *h2m = (TH2D*)f->Get(Form("h2m0%sb",cs)); assert(h2m);
  TH2D *h2x = (TH2D*)f->Get(Form("h2m0%sbx",cs)); assert(h2x);
  TH2D *h2a = (TH2D*)f->Get(Form("h2m2%sb",cs)); assert(h2a);
  TH2D *h2ax = (TH2D*)f->Get(Form("h2m2%sbx",cs)); 
  if (!h2ax) h2ax = (TH2D*)f->Get(Form("h2m2x%sb",cs)); 
  assert(h2ax);

  TH1D *hm = h2m->ProjectionX(Form("hm_%s",cs));
  TH1D *hx = h2x->ProjectionX(Form("hx_%s",cs));
  TH1D *ha = h2a->ProjectionX(Form("ha_%s",cs));
  TH1D *hax = h2ax->ProjectionX(Form("hax_%s",cs));
  TH1D *hjer = h2m->ProjectionX(Form("hjer_%s",cs));
  TH1D *hjera = h2m->ProjectionX(Form("hjera_%s",cs));
  for (int i = 1; i != h2m->GetNbinsX()+1; ++i) {
    TH1D *h1m = h2m->ProjectionY("h1m",i,i);
    TH1D *h1x = h2x->ProjectionY("h1x",i,i);
    TH1D *h1a = h2a->ProjectionY("h1a",i,i);
    TH1D *h1ax = h2ax->ProjectionY("h1ax",i,i);
    double sm = h1m->GetRMS();
    double sx = h1x->GetRMS();
    double sa = h1a->GetRMS();
    double sax = h1ax->GetRMS();
    double em = h1m->GetRMSError();
    double ex = h1x->GetRMSError();
    double ea = h1a->GetRMSError();
    double eax = h1ax->GetRMSError();
    hm->SetBinContent(i, sm);
    hm->SetBinError(i, em);
    hx->SetBinContent(i, sx);
    hx->SetBinError(i, ex);
    ha->SetBinContent(i, sa);
    ha->SetBinError(i, ea);
    hax->SetBinContent(i, sax);
    hax->SetBinError(i, eax);
    double jer = sqrt(max(sm*sm - sx*sx,0.)) / sqrt(2.);
    double err = 0.5*(sqrt(max(pow(sm+em,2) - pow(sx-ex,2),0.))
		      - sqrt(max(pow(sm-em,2) - pow(sx+ex,2),0.))) / sqrt(2);
    hjer->SetBinContent(i, jer);
    hjer->SetBinError(i, err);
    double jera = sqrt(max(sa*sa - sax*sax,0.)) / sqrt(2);
    double erra = 0.5*(sqrt(max(pow(sa+ea,2) - pow(sax-eax,2),0.))
		       - sqrt(max(pow(sa-ea,2) - pow(sax+eax,2),0.))) / sqrt(2);
    hjera->SetBinContent(i, jera);
    hjera->SetBinError(i, erra);
    delete h1m;
    delete h1x;
    delete h1a;
    delete h1ax;
  } // for i

  TH1D *h = tdrHist(Form("h3_%s",cs),"RMS",0,1.0,"#eta",-5.191,5.191);
  //lumi_136TeV = "AlCaRaw, 1 fb^{-1}";
  lumi_136TeV = "AlCaRaw 1.79, Golden 4.86 fb^{-1}";
  TCanvas *c1 = tdrCanvas(Form("c1jer_%s",cs),h,8,11,kSquare);
  tdrDraw(hx,"H",kNone,kRed,kSolid,-1,kNone);
  tdrDraw(hm,"H",kNone,kBlue,kSolid,-1,kNone);
  tdrDraw(ha,"H",kNone,kBlack,kSolid,-1,kNone);
  tdrDraw(hax,"H",kNone,kGray+2,kSolid,-1,kNone);
  tdrDraw(hjer,"H",kNone,kGreen+2,kSolid,-1,kNone);
  tdrDraw(hjera,"H",kNone,kMagenta+2,kSolid,-1,kNone);

  c1->SaveAs(Form("pdf/alcajme/drawAlCaJME_drawJER_pt%s_%s.pdf",cs,cs2));
}

// Compare AlCaRaw to JMENANO data (and MC?)
void drawComp(const char *cs) {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  //TFile *fd = new TFile("../dijet/rootfiles/jmenano_data_out_v2.root","READ");
  TFile *fd = new TFile("../dijet/rootfiles/jmenano_data_out_v5_4p86fb.root",
			"READ");
  assert(fd && !fd->IsZombie());
  fd->cd("HLT_DiPFJetAve40/Pt_40_50"); fd = (TFile*)gDirectory;
  //fd->cd("HLT_DiPFJetAve40/Pt_50_60"); fd = (TFile*)gDirectory;

  //TFile *fa = new TFile("rootfiles/alcajme_out_50M_v5_denom.root","READ");
  //TFile *fa = new TFile("rootfiles/alcajme_out_50M_v8_doak4.root","READ");
  TFile *fa = new TFile("rootfiles/alcajme_out_50M_v10_4p86fb.root","READ");
  assert(fa && !fa->IsZombie());
  fa->cd("Pt_40_50"); fa = (TFile*)gDirectory;
  //fa->cd("Pt_50_60"); fa = (TFile*)gDirectory;
  //fa->cd("Pt_60_70"); fa = (TFile*)gDirectory;
  //fa->cd("Pt_70_85"); fa = (TFile*)gDirectory;
  curdir->cd();

  TProfile *p(0);
  map<string, TProfile*> mp;
  mp["m0a"] = p = (TProfile*)fa->Get("pm0ab"); assert(p);
  mp["m2a"] = p = (TProfile*)fa->Get("pm2ab"); assert(p);
  mp["m0d"] = p = (TProfile*)fd->Get("pm0ab"); assert(p);
  mp["m2d"] = p = (TProfile*)fd->Get("pm2ab"); assert(p);

  lumi_136TeV = "AlCaRaw 1.79, JMENANO 0.450, Golden 4.86 fb^{-1}";
  extraText = "Preliminary";
  TH1D *h = tdrHist("hc","MPF PtAve",0.50,1.2,"#eta",-5.191,5.191);
  TCanvas *c1 = tdrCanvas("cc",h,8,11,kSquare);

  TLine *l = new TLine();
  l->DrawLine(-5.191,1,5.191,1);
  l->DrawLine(-5.191,0,5.191,0);
  
  TH1D *hm2d = fixB(mp["m2d"]);
  TH1D *hm0d = fixB(mp["m0d"]);
  TH1D *hm2a = fixB(mp["m2a"]);
  TH1D *hm0a = fixB(mp["m0a"]);

  TH1D *hdiffd = diffB(hm2d,hm0d);
  TH1D *hdiffa = diffB(hm2a,hm0a);

  tdrDraw(hdiffd,"E2",kNone,kGreen+2,kNone,-1,1001,kGray+2);
  hdiffd->SetFillColorAlpha(kGray+2,0.7);
  tdrDraw(hdiffa,"E2",kNone,kGreen+2,kNone,-1,1001,kGreen+2);
  hdiffa->SetFillColorAlpha(kGreen+2,0.7);
  
  tdrDraw(hm2d,"H",kNone,kBlack,kDotted,-1,kNone);
  tdrDraw(hm0d,"H",kNone,kBlack,kSolid,-1,kNone);
  hm0d->SetLineWidth(2);

  tdrDraw(hm2a,"H",kNone,kGreen+2,kDotted,-1,kNone);
  tdrDraw(hm0a,"H",kNone,kGreen+2,kSolid,-1,kNone);
  hm0a->SetLineWidth(2);    

  TLegend *t1 = tdrLeg(0.40,0.90-2*0.045,0.65,0.90);
  t1->AddEntry(hm0a,"AlCaRaw","L");
  t1->AddEntry(hm0d,"JMENANO","L");
  TLegend *t2 = tdrLeg(0.70,0.90-2*0.045,0.95,0.90);
  t2->AddEntry(hdiffa,"MPF-DB","F");
  t2->AddEntry(hdiffd,"MPF-DB","F");

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  tex->DrawLatex(0.40,0.40,"AK4PUPPI");
  tex->DrawLatex(0.40,0.35,"40<p_{T,ave}<50 GeV");
  //tex->DrawLatex(0.40,0.35,"50<p_{T,ave}<60 GeV");
  tex->SetTextColor(kGreen+2);
  tex->DrawLatex(0.40,0.30,"40<p_{T,ave}<50 GeV");
  //tex->DrawLatex(0.40,0.30,"50<p_{T,ave}<60 GeV");
  //tex->DrawLatex(0.40,0.30,"70<p_{T,ave}<85 GeV");
  tex->DrawLatex(0.40,0.23,"AlCa_PFJet40_v21");
  tex->SetTextColor(kGray+2);
  tex->DrawLatex(0.40,0.17,"HLT_DiPFJetAve40");

  gPad->RedrawAxis();
  c1->SaveAs(Form("pdf/alcajme/drawAlCaJME_drawComp_%s.pdf",cs));
} // drawComp
