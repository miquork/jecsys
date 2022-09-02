// Purpose: draw L2Res from JMENANO at various pTs, merging triggers
//#include "drawAlCaJME.C" // get fixB
#include "TProfile.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

#include "tdrstyle_mod22.C"

#include <fstream>

class trg {
public:
  string name;
  double ptmin;
  double ptmax;
  double etamin;
  double etamax;
  trg(string _name, double _ptmin, double _ptmax,
      double _etamin, double _etamax) {
    name = _name;
    ptmin = _ptmin;
    ptmax = _ptmax;
    etamin = _etamin;
    etamax = _etamax;
  }
};

// fixB and _F originally implented in drawAlCaJME.C
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
  //if (p0>0.35) cout << p->GetName() << ": " << p0 << endl;
  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    if (p0>0.35) {

      // Analytical approximation for JER bias on pT spectrum from
      // voutilainen_thesis.pdf p.259 Eq.(C.8):
      // <pT,ptcl> = pT - alpha*sigma^2 for spectrum N0*exp(-alpha*pT),
      // where alpha = dF/dpT * 1/F
      //
      // #1 MPF = 1+(pT,probe-pT,tag)/pT,tag = pT,ptcl/(pT,ptcl+a*s^2)
      // => MPF*pT+MPF*as2 = pT => as2 = (1-MPF)/MPF * pT
      //
      // #2 MPF = 1+(pT,probe-pT,tag)/pT,probe = 2-pT,ptcl/(pT,ptcl+a*s^2)
      // => (2-MPF)*pT+(2-MPF)*as2 = pT => as2 = (MPF-1)/(2-MPF) * pT
      //
      // #3 MPF = 1+(pT,probe-pT,tag)/pTave = 1+pTdiff/pTave = 1+dpT/(pT+as^2/2)
      // => (MPF-1)*pT+(MPF-1)*as2/2 = dpT => as2 = (d-(MPF-1))/(MPF-1)*2 * pT
      // (just that MPF=1 for barrel so unsolvable => use as2 from tag, probe)

      double pt = 40.;
      double eta = h->GetBinCenter(i);
      F->SetParameter(0,0.);
      double df0 = F->Derivative(pt) / F->Eval(pt);
      F->SetParameter(0,eta);
      double df = F->Derivative(pt) / F->Eval(pt);

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

// Symmetrize TH2D vs |eta|
bool doSymmetrize = true;
TH2D *symmetrize(TH2D *h2) {

  TH2D *h2s = (TH2D*)h2->Clone(Form("%s_symm",h2->GetName()));
  for (int i1 = 1; i1 != h2s->GetNbinsX()+1; ++i1) {
    double eta = h2s->GetXaxis()->GetBinCenter(i1);
    int i2 = h2s->GetXaxis()->FindBin(-eta);
    for (int j = 1; j != h2s->GetNbinsY()+1; ++j) {
      
      double y1 = h2s->GetBinContent(i1, j);
      double ey1 = h2s->GetBinError(i1, j);
      double y2 = h2s->GetBinContent(i2, j);
      double ey2 = h2s->GetBinError(i2, j);
      double y(0), ey(0);
      if (eta<0) { // calculate average of plus and minus
	y = (ey1*ey2!=0 ? 0.5*(y1+y2) : (ey1>0 ? y1 : (ey2>0 ? y2 : 0.)));
	ey = (ey1*ey2!=0 ? 0.5*(ey1+ey2) : (ey1>0 ? ey1 : (ey2>0 ? ey2 : 0.)));
      }
      else { // replace plus with previously calculated average
	y = y2;
	ey = ey2;
      }
      int i = i1;
      h2s->SetBinContent(i, j, y);
      h2s->SetBinError(i, j, ey);
    } // for k
  } // for i

  return h2s;
} // void symmetrize


void drawJMENANOs(int iFile, TFile *fout); // forward declaration
void drawPF(string trg, int ptmin, int ptmax);

void drawJMENANO() {

  cout << "drawJMENANO.C" << endl << flush;
  /*
  TFile *fout = new TFile("rootfiles/drawJMENANO.root","RECREATE");
  drawJMENANOs(0, fout);
  drawJMENANOs(1, fout);
  drawJMENANOs(2, fout);
  fout->Close();
  */

  drawPF("HLT_DiPFJetAve40",40,50); // bad mc
  drawPF("HLT_DiPFJetAve40",50,60);
  drawPF("HLT_DiPFJetAve40",60,70);
  //drawPF("HLT_DiPFJetAve60",60,70); // bad mc
  drawPF("HLT_DiPFJetAve40",70,85); // mc ok, data same as 60?
  //drawPF("HLT_DiPFJetAve60",70,85); // less mc, data same as 40?
  //drawPF("HLT_DiPFJetAve40",85,100); // too little data for 80-85
  drawPF("HLT_DiPFJetAve60",85,100); // goldilocks
  //drawPF("HLT_DiPFJetAve80",85,100); // too little mc for 80-85
  //drawPF("HLT_DiPFJetAve60",100,125); // mc ok, data same as 80?
  drawPF("HLT_DiPFJetAve80",100,125); // mc ok, data same as 60?
  drawPF("HLT_DiPFJetAve80",125,155);
  //drawPF("HLT_DiPFJetAve80",155,180); // data bad
  drawPF("HLT_DiPFJetAve140",155,180); // data,mc ok
  drawPF("HLT_DiPFJetAve140",180,210);
  //drawPF("HLT_DiPFJetAve140",210,250); // less data
  drawPF("HLT_DiPFJetAve200",210,250); // data fair, mc ok
  //drawPF("HLT_DiPFJetAve140",250,300); // bad data
  drawPF("HLT_DiPFJetAve200",250,300);
  drawPF("HLT_DiPFJetAve260",300,350);
  drawPF("HLT_DiPFJetAve260",350,400); // data chf lower?
  //drawPF("HLT_DiPFJetAve320",350,400); // data chf higher?
  drawPF("HLT_DiPFJetAve320",400,500);
  //drawPF("HLT_DiPFJetAve400",400,500);
  drawPF("HLT_DiPFJetAve400",500,600);
  drawPF("HLT_DiPFJetAve500",500,600);
  drawPF("HLT_DiPFJetAve500",600,800);
  drawPF("HLT_DiPFJetAve500",800,1000);
  drawPF("HLT_DiPFJetAve500",1200,1500);
  drawPF("HLT_DiPFJetAve500",1500,1800);
  drawPF("HLT_DiPFJetAve500",1800,2100);
}

void drawJMENANOs(int iFile, TFile *fout) {

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  bool isMC = (iFile==0);
  bool isDT = (iFile==1);
  bool isRatio = (iFile==2);
  assert(isMC || isDT || isRatio);

  const char *cs = (isMC ? "MC" : (isDT ? "DT" : "DM"));

  TFile *f(0);
  if (isMC)
    f = new TFile("../dijet/rootfiles/jmenano_mc_out_v6b.root","READ");
  //fd = new TFile("../dijet/rootfiles/jmenano_mc_out_v6.root","READ");
  //fd = new TFile("../dijet/rootfiles/jmenano_mc_out_v5.root","READ");
  else
    f = new TFile("../dijet/rootfiles/jmenano_data_out_v6.root","READ");
  //fd = new TFile("../dijet/rootfiles/jmenano_data_out_v5_4p86fb.root","READ");
  assert(f && !f->IsZombie());

  TFile *fa = new TFile("rootfiles/alcajme_out_50M_v10_4p86fb.root","READ");
  assert(fa && !fa->IsZombie());
  TProfile *pam = (TProfile*)fa->Get("Pt_40_50/pm0ab"); assert(pam);
  TProfile *pad = (TProfile*)fa->Get("Pt_40_50/pm2ab"); assert(pad);
  TProfile *pjd = (TProfile*)f->Get("HLT_DiPFJetAve40/Pt_40_50/pm2ab"); assert(pjd);
  TH1D *ham = fixB(pam);
  TH1D *had = fixB(pad);
  TH1D *hjd = fixB(pjd);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);

  TH1D *h = tdrHist("h","p_{T,ave}",15,3400,"#eta_{jet}",-5.2,5.2);
  lumi_136TeV = (isMC ? "FlatQCD" : "RunC, 2.90 of 4.86 fb^{-1}");
  extraText = "WIP";
  TCanvas *c1 = tdrCanvas("c1",h,8,0,kSquare);

  gPad->SetRightMargin(0.20);
  gPad->SetLeftMargin(0.20);
  h->GetYaxis()->SetTitleOffset(1.6);
  gPad->SetLogy();
  h->GetYaxis()->SetMoreLogLabels();
  h->GetYaxis()->SetNoExponent();

   // Copy L2Res histograms for multiple pT bins
   const double vpt[] = {40, 50, 60, 70, 85, 100, 125, 155, 180, 210, 250, 300,
			 350, 400, 500, 600, 800, 1000, 1200, 1500,
			 1800, 2100, 2400, 2700, 3000};
   const int npt = sizeof(vpt)/sizeof(vpt[0])-1;

  // Listing of available triggers
   vector<trg> vtrg;
   //vtrg.push_back(trg("HLT_DiPFJetAve40",40,70,0,5.2));//40,60
   //vtrg.push_back(trg("HLT_DiPFJetAve40",60,70,0,2.853));//40,60
   //vtrg.push_back(trg("HLT_DiPFJetAve60",70,100,0,2.853)); //60,85
   vtrg.push_back(trg("HLT_DiPFJetAve40",40,60,0,5.2)); // v6f
   vtrg.push_back(trg("HLT_DiPFJetAve60",60,100,0,2.853)); // v6f
   vtrg.push_back(trg("HLT_DiPFJetAve80",100,155,0,2.853)); //85,155
   vtrg.push_back(trg("HLT_DiPFJetAve140",155,250,0,2.853)); //180,210
   vtrg.push_back(trg("HLT_DiPFJetAve200",250,300,0,2.853)); // 250,300
   vtrg.push_back(trg("HLT_DiPFJetAve260",300,350,0,2.853));
   vtrg.push_back(trg("HLT_DiPFJetAve320",350,400,0,2.853));
   vtrg.push_back(trg("HLT_DiPFJetAve400",400,500,0,2.853));
   vtrg.push_back(trg("HLT_DiPFJetAve500",500,3000,0,2.853));
   /*
   vtrg.push_back("HLT_PFJet40");
   vtrg.push_back("HLT_PFJet60");
   vtrg.push_back("HLT_PFJet80");
   //vtrg.push_back("HLT_PFJet110");
   vtrg.push_back("HLT_PFJet140");
   vtrg.push_back("HLT_PFJet200");
   vtrg.push_back("HLT_PFJet260");
   vtrg.push_back("HLT_PFJet320");
   vtrg.push_back("HLT_PFJet450");
   vtrg.push_back("HLT_PFJet500");
   vtrg.push_back("HLT_PFJet550");
   */
   vtrg.push_back(trg("HLT_DiPFJetAve60_HFJEC",60,100,2.853,5.2)); // v6f
   vtrg.push_back(trg("HLT_DiPFJetAve80_HFJEC",100,125,2.853,5.2));
   vtrg.push_back(trg("HLT_DiPFJetAve100_HFJEC",125,180,2.853,5.2));
   vtrg.push_back(trg("HLT_DiPFJetAve160_HFJEC",180,250,2.853,5.2));
   vtrg.push_back(trg("HLT_DiPFJetAve220_HFJEC",250,350,2.853,5.2));
   vtrg.push_back(trg("HLT_DiPFJetAve300_HFJEC",350,3000,2.853,5.2));
   /*
   vtrg.push_back("HLT_PFJetFwd15");
   vtrg.push_back("HLT_PFJetFwd25");
   vtrg.push_back("HLT_PFJetFwd40");
   vtrg.push_back("HLT_PFJetFwd60");
   vtrg.push_back("HLT_PFJetFwd80");
   vtrg.push_back("HLT_PFJetFwd140");
   vtrg.push_back("HLT_PFJetFwd200");
   vtrg.push_back("HLT_PFJetFwd260");
   vtrg.push_back("HLT_PFJetFwd320");
   vtrg.push_back("HLT_PFJetFwd400");
   vtrg.push_back("HLT_PFJetFwd450");
   vtrg.push_back("HLT_PFJetFwd500");
   */
   int ntrg = vtrg.size();
  
   // Regular L2Relative and L2Res eta binning
   double vx[] =
     {-5.191,
      -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, -3.489,
      -3.314, -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172, -2.043,
      -1.93, -1.83, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218,
      -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, -0.435, 
      -0.348, -0.261, -0.174, -0.087, 0, 0.087, 0.174, 0.261, 0.348, 0.435,
      0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305,
      1.392, 1.479, 1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5,
      2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191,
      4.363, 4.538, 4.716, 4.889, 5.191};
   const int nx = sizeof(vx)/sizeof(vx[0])-1;

   TH2D *h2 = new TH2D(Form("h2_%s",cs),
		       ";#eta_{jet};p_{T,ave};JES",nx,vx,npt,vpt);
   //TH2D *h2s = new TH2D("h2s",";#eta_{jet};p_{T,ave};MPF",nx,vx,npt,vpt);
   //TH2D *h2d = new TH2D("h2d",";#eta_{jet};p_{T,ave};DB",nx,vx,npt,vpt);
   TH2D *h2m = (isRatio ? (TH2D*)fout->Get("h2s_MC") : 0);
   assert(h2m || !isRatio);

   for (int itrg = 0; itrg != ntrg; ++itrg) {

     trg &t = vtrg[itrg];
     f->cd(t.name.c_str());
     TDirectory *d = gDirectory;

     for (int ipt = 0; ipt != npt; ++ipt) {
       
       double ptmin = vpt[ipt];
       double ptmax = vpt[ipt+1];

       if (ptmin >= t.ptmin && ptmax <= t.ptmax) {
   
	 d->cd(Form("Pt_%d_%d",int(ptmin+0.5),int(ptmax+0.5)));
	 TDirectory *dd = gDirectory;
	 TProfile *pm = (TProfile*)dd->Get("pm0ab"); assert(pm);
	 //TProfile *pm = (TProfile*)dd->Get("pm2ab"); assert(pm);
	 //TProfile *pd = (TProfile*)dd->Get("pm2ab"); assert(pd);
	 
	 TH1D *hm = fixB(pm); hm->SetName(Form("hm_%d",ipt));
	 //TH1D *hd = fixB(pd); hd->SetName(Form("hd_%d",ipt));

	 for (int ieta = 0; ieta != nx; ++ieta) {

	   double etamin = vx[ieta];
	   double etamax = vx[ieta+1];

	   if ((etamin>=0 && (etamin >= t.etamin && etamax <= t.etamax)) ||
	       (etamin<0 && (etamin >= -t.etamax && etamax <= -t.etamin))) {

	     /*
	     h2->SetBinContent(ieta+1,ipt+1, pm->GetBinContent(ieta+1));
	     h2->SetBinError(ieta+1,ipt+1, pm->GetBinError(ieta+1));
	     h2->SetBinContent(ieta+1,ipt+1, pd->GetBinContent(ieta+1));
	     h2->SetBinError(ieta+1,ipt+1, pd->GetBinError(ieta+1));
	     */
	     double ref = (h2m ? h2m->GetBinContent(ieta+1,ipt+1) : 1);
	     if (ref!=0) {
	       h2->SetBinContent(ieta+1,ipt+1, hm->GetBinContent(ieta+1)/ref);
	       h2->SetBinError(ieta+1,ipt+1, hm->GetBinError(ieta+1)/ref);
	     }
	     else {
	       h2->SetBinContent(ieta+1,ipt+1, 0.);
	       h2->SetBinError(ieta+1,ipt+1, 0.);
	     }
	   } // eta
	 } // ieta
       } // pt
     } // ipt
   } // itrg 

   //if (doSymmetrize) 
   h2 = symmetrize(h2);

   h2->Draw("SAME COLZ");
   h2->GetZaxis()->SetRangeUser(0.5+1e-4,1.5-1e-4);
   h2->GetZaxis()->SetTitleOffset(1.2);
   gPad->RedrawAxis();
   gPad->Update();

   c1->SaveAs(Form("pdf/alcajme/drawJMENANO_2D_%s.pdf",cs));
   h2->GetZaxis()->SetRangeUser(0.8+1e-4,1.2-1e-4);
   c1->SaveAs(Form("pdf/alcajme/drawJMENANO_2D_zoom_%s.pdf",cs));
   h2->GetZaxis()->SetRangeUser(0.9+1e-4,1.1-1e-4);
   c1->SaveAs(Form("pdf/alcajme/drawJMENANO_2D_zoom2_%s.pdf",cs));

   // Run smoothing vs pT, keep track of chi2/NDF
   TH2D *h2s = (TH2D*)h2->Clone(Form("h2s%s",cs));
   h2s->Reset();
   TH1D *hchi20 = h2s->ProjectionX(Form("hchi20_%s",cs));
   TH1D *hchi21 = h2s->ProjectionX(Form("hchi21_%s",cs));
   TH1D *hchi22 = h2s->ProjectionX(Form("hchi22_%s",cs));
   TH1D *hjes40 = h2s->ProjectionX(Form("hjes40_%s",cs));
   TH1D *hjes15b = h2s->ProjectionX(Form("hjes15b_%s",cs));
   TH1D *hjesmaxb = h2s->ProjectionX(Form("hjesmaxb_%s",cs));
   TH1D *hjes40orig = h2->ProjectionX(Form("hjes40orig_%s",cs),1,1);
   double i60 = h2->GetYaxis()->FindBin(60.);
   TH1D *hjes60orig = h2->ProjectionX(Form("hjes100orig_%s",cs),i60,i60);
   double i100 = h2->GetYaxis()->FindBin(100.);
   TH1D *hjes100orig = h2->ProjectionX(Form("hjes100orig_%s",cs),i100,i100);
   TH1D *hdpt = h2s->ProjectionX(Form("hdpt_%s",cs));
   TH1D *hdpt0 = h2s->ProjectionX(Form("hdpt0_%s",cs));
   for (int i = 1; i != h2->GetNbinsX()+1; ++i) {

     double eta = h2->GetXaxis()->GetBinCenter(i);
     double emax1 = 6800.*0.5;
     double emax2 = 4500.;
     double ptmax = (fabs(eta)>2.5 ? emax2 : emax1) / cosh(eta);
     TH1D *h = h2->ProjectionY("htmp",i,i);
     
     // Quick cleaning
     for (int i = 1; i != h->GetNbinsX()+1; ++i) {
       if (h->GetBinError(i)<0.001) {
	 h->SetBinContent(i, 0.);
	 h->SetBinError(i, 0.);
       }
     }

     TF1 *f0 = new TF1(Form("f0_%d",i),"[0]",40,ptmax);
     h->Fit(f0,"QRN");
     TF1 *f1 = new TF1(Form("f1_%d",i),"[0]+[1]*log(x)",40,ptmax);
     f1->SetParameters(f0->GetParameter(0), 0.);
     h->Fit(f1,"QRN");
     TF1 *f2 = new TF1(Form("f2_%d",i),"[0]+[1]*log(x)+[2]/x",40,ptmax);
     f2->SetParameters(f0->GetParameter(0), f1->GetParameter(1), 0.);
     h->Fit(f2,"QRN");

     hchi20->SetBinContent(i, f0->GetNDF()>0 ? f0->GetChisquare()/f0->GetNDF() : 0);
     hchi20->SetBinError(i, f0->GetNDF()>0 ? f0->GetChisquare()/f0->GetNDF()*1./sqrt(f0->GetNDF()) : 0);

     hchi21->SetBinContent(i, f1->GetNDF()>0 ? f1->GetChisquare()/f1->GetNDF() : 0);
     hchi21->SetBinError(i, f1->GetNDF()>0 ? f1->GetChisquare()/f1->GetNDF()*1./sqrt(f1->GetNDF()) : 0);

     hchi22->SetBinContent(i, f2->GetNDF()>0 ? f2->GetChisquare()/f2->GetNDF() : 0);
     hchi22->SetBinError(i, f2->GetNDF()>0 ? f2->GetChisquare()/f2->GetNDF()*1./sqrt(f2->GetNDF()) : 0);


     double nsig = 2.5;//1.5;
     //TF1 *fjes = (fabs(f1->GetParameter(1))>nsig*f1->GetParError(1) ? f1 : f0);
     TF1 *fjes = f0; // default choice
     if (fabs(f1->GetParameter(1))>nsig*f1->GetParError(1)) fjes = f1;
     //if (fabs(f2->GetParameter(2))>nsig*f2->GetParError(2)) fjes = f2;
     if (isMC && (fabs(eta)>2.5 || fabs(eta)<1.0)) fjes = f0;
     //if (isRatio && (fabs(eta)>4.191 || fabs(eta)<1.0)) fjes = f0;
     //if (!isMC && (fabs(eta)<1.0)) fjes = f0;
  
     hjes40->SetBinContent(i, fjes->GetNDF()>=0 ? fjes->Eval(40.) : 1);

     hjes15b->SetBinContent(i, fjes->GetNDF()>=0 ?
			    0.5*(fjes->Eval(40.)+fjes->Eval(15.)) : 1);
     hjes15b->SetBinError(i, fjes->GetNDF()>=0 ?
			  0.5*fabs(fjes->Eval(40.)-fjes->Eval(15.)) : 1);

     hjesmaxb->SetBinContent(i, fjes->GetNDF()>=0 ?
			     0.5*(fjes->Eval(40.)+fjes->Eval(ptmax)) : 1);
     hjesmaxb->SetBinError(i, fjes->GetNDF()>=0 ?
			   0.5*fabs(fjes->Eval(40.)-fjes->Eval(ptmax)) : 1);

     hdpt->SetBinContent(i, f1->GetNDF()>=0 ? 
			 100.*f1->GetParameter(1) : 0);
     hdpt->SetBinError(i, f1->GetNDF()>=0 ? 
		       100.*f1->GetParError(1) : 0);

     hdpt0->SetBinContent(i, fjes==f1 ? 100.*f1->GetParameter(1) : 0);
     hdpt0->SetBinError(i, fjes==f1 ? 100.*f1->GetParError(1) : 0);

     for (int j = 1; j != h2->GetNbinsY()+1; ++j) {
       if (//h2->GetBinContent(i,j)!=0 &&
	   h2->GetYaxis()->GetBinLowEdge(j)<ptmax) {
	 double pt = h2->GetYaxis()->GetBinCenter(j);
	 h2s->SetBinContent(i, j, fjes->Eval(pt));
       }
       else {
	 h2s->SetBinContent(i, j, 0.);
       }
     }

     delete h;
   } // for i


   TH1D *hs = tdrHist("hs","p_{T,ave}",15,3400,"#eta_{jet}",-5.2,5.2);
   TCanvas *c1s = tdrCanvas("c1s",hs,8,0,kSquare);

   gPad->SetRightMargin(0.20);
   gPad->SetLeftMargin(0.20);
   hs->GetYaxis()->SetTitleOffset(1.6);
   gPad->SetLogy();
   hs->GetYaxis()->SetMoreLogLabels();
   hs->GetYaxis()->SetNoExponent();

   h2s->Draw("SAME COLZ");
   h2s->GetZaxis()->SetRangeUser(0.5+1e-4,1.5-1e-4);
   gPad->RedrawAxis();
   gPad->Update();

   c1s->SaveAs(Form("pdf/alcajme/drawJMENANO_2DS_%s.pdf",cs));
   h2s->GetZaxis()->SetRangeUser(0.8+1e-4,1.2-1e-4);
   c1s->SaveAs(Form("pdf/alcajme/drawJMENANO_2DS_zoom_%s.pdf",cs));
   h2s->GetZaxis()->SetRangeUser(0.9+1e-4,1.1-1e-4);
   c1s->SaveAs(Form("pdf/alcajme/drawJMENANO_2DS_zoom2_%s.pdf",cs));


   TH1D *h2c = tdrHist("h2c","#chi^{2} / NDF",0,10,"#eta_{jet}",-5.2,5.2);
   TCanvas *c2 = tdrCanvas("c2",h2c,8,0,kSquare);
   tdrDraw(hchi20,"HISTE",kNone,kOrange+2,kSolid,-1,1001,kOrange+1);
   tdrDraw(hchi21,"HISTE",kNone,kYellow+2,kSolid,-1,1001,kYellow+1);
   //tdrDraw(hchi22,"HISTE",kNone,kGreen+2,kSolid,-1,1001,kGreen-9);
   hchi21->SetFillColorAlpha(kYellow+1,0.7);
   hchi22->SetFillColorAlpha(kGreen-9,0.7);
   gPad->RedrawAxis();
   l->DrawLine(-5.2,1,+5.2,1);
   c2->SaveAs(Form("pdf/alcajme/drawJMENANO_chi2_%s.pdf",cs));

   //TH1D *h3 = tdrHist("h3","JES at p_{T}=40 GeV (p_{T,min} to E_{max})",
   TH1D *h3 = tdrHist("h3","JES from p_{T,min} to E_{max}",
		      0.38+1e-4,1.28-1e-4,"#eta_{jet}",-5.2,5.2);
   TCanvas *c3 = tdrCanvas("c3",h3,8,0,kSquare);
   h3->GetXaxis()->SetTitleOffset(0.8);
   h3->GetYaxis()->SetTitleOffset(1.18);
   tdrDraw(hjes40,"HIST",kNone,kYellow+2,kSolid,-1,1001,kYellow+1);
   tdrDraw(hjes15b,"E2",kNone,kCyan+2,kSolid,-1,1001,kCyan+1);
   hjes15b->SetFillColorAlpha(kCyan+1,0.7);
   tdrDraw(hjesmaxb,"E2",kNone,kOrange+2,kSolid,-1,1001,kOrange+1);
   hjesmaxb->SetFillColorAlpha(kOrange+1,0.7);
   if (!isMC) tdrDraw(hjes40orig,"P",kFullCircle,kBlack,kSolid,-1,kNone,0);
   //if (!isMC) tdrDraw(hjes60orig,"P",kFullCircle,kBlack,kSolid,-1,kNone,0);
   if (isMC)  tdrDraw(hjes100orig,"P",kFullCircle,kBlack,kSolid,-1,kNone,0);
   hjes40orig->SetMarkerSize(0.5);
   hjes100orig->SetMarkerSize(0.5);
   gPad->RedrawAxis();

   TLegend *leg = tdrLeg(0.35,0.90-4*0.035,0.55,0.90);
   leg->SetTextSize(0.035);
   leg->AddEntry(hjes40orig,isMC ? "JES @ 100 GeV" : "JES @ 40 GeV","PLE");
   //leg->AddEntry(hjes40orig,"JES @ 40 GeV","PLE");
   //leg->AddEntry(hjes40orig,isMC ? "MPF @ 100 GeV" : "MPF @ 60 GeV","PLE");
   leg->AddEntry(hjes40,"Fit @ 40 GeV","F");
   leg->AddEntry(hjes15b,"Fit @ 15 GeV","F");
   leg->AddEntry(hjesmaxb,"Fit @ E_{max}","F");

   l->DrawLine(-5.2,1,+5.2,1);

   c3->SaveAs(Form("pdf/alcajme/drawJMENANO_jes40_%s.pdf",cs));

   tdrDraw(hjd,"Pz",kOpenDiamond,kGreen+2,kSolid,-1,kNone,0);
   tdrDraw(had,"Pz",kOpenSquare,kBlue,kSolid,-1,kNone,0);
   tdrDraw(ham,"Pz",kOpenCircle,kBlack,kSolid,-1,kNone,0);

   hjd->SetMarkerSize(0.7);
   had->SetMarkerSize(0.7);
   ham->SetMarkerSize(0.7);

   TLegend *leg3 = tdrLeg(0.35,0.20,0.65,0.20+3*0.045);
   leg3->AddEntry(hjd,"DB @ 40 GeV","PLE");
   leg3->AddEntry(had,"HLT-DB @ 40 GeV","PLE");
   leg3->AddEntry(ham,"HLT-MHPF @ 40 GeV","PLE");

   c3->SaveAs(Form("pdf/alcajme/drawJMENANO_jes40_v2_%s.pdf",cs));

   TH1D *h4 = tdrHist("h4","dJES/dlog(p_{T}) (%)",-30,30,"#eta_{jet}",-5.2,5.2);
   TCanvas *c4 = tdrCanvas("c4",h4,8,0,kSquare);
   tdrDraw(hdpt,"HISTE",kNone,kYellow+2,kSolid,-1,1001,kYellow+1);
   tdrDraw(hdpt0,"HISTE",kNone,kOrange+2,kSolid,-1,1001,kOrange+1);
   gPad->RedrawAxis();

   l->SetLineStyle(kSolid);
   l->DrawLine(-5.2,0,+5.2,0);

   c4->SaveAs(Form("pdf/alcajme/drawJMENANO_dpt_%s.pdf",cs));

   // Save combined and smoothed results
   fout->cd();
   h2->Write(Form("h2_%s",cs),TObject::kOverwrite);
   h2s->Write(Form("h2s_%s",cs),TObject::kOverwrite);
   curdir->cd();

   // Map <PtAve> to <PtRaw> = JES(eta, PtAve)*<PtAve>
   // and JES to correction=1./JES for L2Residuals
   for (int i = 1; i != h2->GetNbinsX()+1; ++i) {

   }

   // Write out text file
   bool writeTextFile = false;
   if (writeTextFile && iFile==2) {
     ofstream fout("textFiles/drawJMENANO_Run2022C-PromptReco_DATA_L2Res.txt",
		   ios::out);
     fout << "{1 JetEta 1 JetPt [0]+[1]*log(x) Correction L2Relative}" << endl;
   }
} // drawJMENANO


void drawPF(string trg, int ptmin, int ptmax) {

  TDirectory *curdir = gDirectory;
  setTDRStyle();

  TFile *fm = new TFile("../dijet/rootfiles/jmenano_mc_out_v7.root","READ");
  TFile *fd = new TFile("../dijet/rootfiles/jmenano_data_out_v7.root","READ");
  assert(fm && !fm->IsZombie());
  assert(fd && !fd->IsZombie());

  fd->cd(trg.c_str());
  gDirectory->cd(Form("Pt_%d_%d",ptmin,ptmax));
  TDirectory *dd = gDirectory;
  fm->cd(trg.c_str());
  gDirectory->cd(Form("Pt_%d_%d",ptmin,ptmax));
  TDirectory *dm = gDirectory;
  
  string vpf[] = {"muf","cef","chf","nhf","nef","hfhf","hfef"};
  const int npf = sizeof(vpf)/sizeof(vpf[0]);
  map<string,int> color;
  color["muf"] = kMagenta+1;
  color["cef"] = kCyan+1;
  color["chf"] = kRed;
  color["nhf"] = kGreen+2;
  color["nef"] = kBlue;
  color["hfhf"] = kMagenta+1;
  color["hfef"] = kCyan+1;
  map<string,const char*> name;
  name["muf"] = "Muons";
  name["cef"] = "Electrons";
  name["chf"] = "Charged hadrons";
  name["nhf"] = "Neutral hadrons";
  name["nef"] = "Photons";
  name["hfhf"] = "HF hadrons";
  name["hfef"] = "HF EM";
  
  TH1D *h = tdrHist("h","PF energy fraction",0,1.2,"#eta_{jet}",-5.2,5.2);
  TH1D *h2 = tdrHist("h","Data-MC (%)",-10,10,"#eta_{jet}",-5.2,5.2);
  lumi_136TeV = "RunC, 2.90 of 4.86 fb^{-1}";
  TCanvas *c1 = tdrDiCanvas("c1",h,h2,8,11);

  TLegend *leg = tdrLeg(0.43,0.92-npf*0.030,0.68,0.92);
  leg->SetTextSize(0.035);
  TLatex *tex = new TLatex(); tex->SetNDC(); tex->SetTextSize(0.035);
  tex->DrawLatex(0.37,0.55,"Markers data, histos MC");
  tex->DrawLatex(0.65,0.90,Form("%4d<p_{T}<%d GeV",ptmin,ptmax));
  
  for (int i = 0; i != npf; ++i) {
    const char *cpf = vpf[i].c_str();
    TProfile *pd = (TProfile*)dd->Get(Form("p%s",cpf)); assert(pd);
    TProfile *pm = (TProfile*)dm->Get(Form("p%s",cpf)); assert(pm);

    c1->cd(1);
    tdrDraw(pm,"HISTE",kNone,color[cpf],kSolid,-1,kNone,kNone);
    tdrDraw(pd,"Pz",kFullCircle,color[cpf],kSolid,-1,kNone,kNone);
    pd->SetMarkerSize(0.5);
    gPad->RedrawAxis();

    leg->AddEntry(pd,name[cpf],"PLE");
    
    c1->cd(2);
    TH1D *hd = pd->ProjectionX();
    hd->Add(pm,-1);
    hd->Scale(100.);
    tdrDraw(hd,"Pz",kFullCircle,color[cpf],kSolid,-1,kNone,kNone);
    hd->SetMarkerSize(0.5);
    gPad->RedrawAxis();
  } // for i

  c1->SaveAs(Form("pdf/alcajme/PF/%s_%d_%d.pdf",trg.c_str(),ptmin,ptmax));
} // drawPF
