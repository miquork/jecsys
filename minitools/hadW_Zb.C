// Purpose: Use hadW_cs.C inspired code for Z+b
//          to estimate b, c and g JES relative to uds
#include "TFile.h"
#include "TH2D.h"
#include "TChain.h"
#include "TCut.h"
#include "TF1.h"
#include "TProfile2D.h"
#include "TBox.h"
#include "TGraphErrors.h"
#include "TLatex.h"

#include "../tdrstyle_mod15.C"

// Fit function for ABCD fractions
// x=bin index in TH2D
// p[0]=light c-mistag rate modifier (worse)
// p[1]=charm c-efficiency modifier (worse)
// p[2]=light+charm g-mistag rate modified (worse)
// p[3]=gluon (others) rate modifier (worse), e.g. due to FSR combinatorics
// Input MC counts are in h2m
// Output data counts are in h2d
TH2D *h2nm(0), *h2nd(0), *h2nf(0);
Double_t fnABCD(Double_t *x, Double_t *p);

TH2D *h2mm(0), *h2md(0), *h2mf(0);
Double_t frABCD(Double_t *x, Double_t *p);

bool debug = false;

double minpt = 30;//15;//30
double fitminpt = 35;

bool useIncRef = true;//false;

bool testMC = false;//true;
bool isEta25 = true;
bool showNDcomp = false; // Zoom data/MC for N out to show decomposition
//bool fixParsNeff = true;
bool fixParNgeff = true;//false;//true;
bool fixParNbeff = true;
bool fixParNceff = true;
bool fixParNg = true;//false;//true;
bool invertNgNeff = false; // invert SF for Ng and Neff
bool fixParNb = true;//false;//true;
bool fixParNc = true;//false;//true;
bool fixParNo = true;//false;//true;
//bool fixParsNall = false;
bool fixParsMCN = false;//false;//true;
bool constrainParsN = true;//false;

struct zbFit {

  double nd, nde;
  double nm, nme;


  TH2D *h2nd;
  TH2D *h2nm;

  //double nq, nqe;
  //double ng, nge;
  //double nc, nce;
  //double nb, nbe;
  //double no, noe;

  double nchi2;
  double nndf;

  TH2D *h2rd;
  TH2D *h2rm;
  
  double rq, rqe;
  double rg, rge;
  double rc, rce;
  double rb, rbe;
  double ro, roe;

  double rchi2;
  double rndf;  
};

zbFit hadW_Zbs(double pt);

void hadW_Zb() {

  TFile *f = new TFile("rootfiles/jecdata2018ABCD.root","READ");
  assert(f && !f->IsZombie());

  // MC/data ratios of event counts for each region
  TGraphErrors *gmi = new TGraphErrors(0);
  TGraphErrors *gmq = new TGraphErrors(0);
  TGraphErrors *gmg = new TGraphErrors(0);
  TGraphErrors *gmc = new TGraphErrors(0);
  TGraphErrors *gmb = new TGraphErrors(0);
  // ...and same for fit/data ratio
  TGraphErrors *gfi = new TGraphErrors(0);
  TGraphErrors *gfq = new TGraphErrors(0);
  TGraphErrors *gfg = new TGraphErrors(0);
  TGraphErrors *gfc = new TGraphErrors(0);
  TGraphErrors *gfb = new TGraphErrors(0);

  // Flavor fractions in data
  TGraphErrors *gndi = new TGraphErrors(0);
  TGraphErrors *gndq = new TGraphErrors(0);
  TGraphErrors *gndg = new TGraphErrors(0);
  TGraphErrors *gndc = new TGraphErrors(0);
  TGraphErrors *gndb = new TGraphErrors(0);
  // ...and their data/MC ratio
  TGraphErrors *gni = new TGraphErrors(0);
  TGraphErrors *gnq = new TGraphErrors(0);
  TGraphErrors *gng = new TGraphErrors(0);
  TGraphErrors *gnc = new TGraphErrors(0);
  TGraphErrors *gnb = new TGraphErrors(0);

  // Flavor responses in data
  TGraphErrors *grdi = new TGraphErrors(0);
  TGraphErrors *grdq = new TGraphErrors(0);
  TGraphErrors *grdg = new TGraphErrors(0);
  TGraphErrors *grdc = new TGraphErrors(0);
  TGraphErrors *grdb = new TGraphErrors(0);
  // ...in MC
  TGraphErrors *grmi = new TGraphErrors(0);
  TGraphErrors *grmq = new TGraphErrors(0);
  TGraphErrors *grmg = new TGraphErrors(0);
  TGraphErrors *grmc = new TGraphErrors(0);
  TGraphErrors *grmb = new TGraphErrors(0);
  // ...and their data/MC ratio
  TGraphErrors *gri = new TGraphErrors(0);
  TGraphErrors *grq = new TGraphErrors(0);
  TGraphErrors *grg = new TGraphErrors(0);
  TGraphErrors *grc = new TGraphErrors(0);
  TGraphErrors *grb = new TGraphErrors(0);

  TH1D *hz = (TH1D*)f->Get("data/eta00-25/hdm_mpfchs1_zi"); assert(hz);
  for (int i = 1; i != hz->GetNbinsX(); ++i) {
    if (hz->GetBinContent(i)!=0) {
      zbFit zb = hadW_Zbs(hz->GetBinCenter(i));

      // Counts per region in data
      double ndi = h2nd->GetBinContent(1,1);
      double ndq = h2nd->GetBinContent(2,1) / ndi;
      double ndg = h2nd->GetBinContent(3,1) / ndi;
      double ndc = h2nd->GetBinContent(4,1) / ndi;
      double ndb = h2nd->GetBinContent(5,1) / ndi;
      double ndie = h2nd->GetBinError(1,1);
      double ndqe = h2nd->GetBinError(2,1) / ndi;
      double ndge = h2nd->GetBinError(3,1) / ndi;
      double ndce = h2nd->GetBinError(4,1) / ndi;
      double ndbe = h2nd->GetBinError(5,1) / ndi;
      ndie /= ndi;
      ndi  /= ndi;
      // ...in MC,
      double nmi = zb.h2nm->GetBinContent(1,1);
      double nmq = zb.h2nm->GetBinContent(2,1) / nmi;
      double nmg = zb.h2nm->GetBinContent(3,1) / nmi;
      double nmc = zb.h2nm->GetBinContent(4,1) / nmi;
      double nmb = zb.h2nm->GetBinContent(5,1) / nmi;
      double nmie = zb.h2nm->GetBinError(1,1);
      double nmqe = zb.h2nm->GetBinError(2,1) / nmi;
      double nmge = zb.h2nm->GetBinError(3,1) / nmi;
      double nmce = zb.h2nm->GetBinError(4,1) / nmi;
      double nmbe = zb.h2nm->GetBinError(5,1) / nmi;
      nmie /= nmi;
      nmi  /= nmi;
      // ...and in fit
      double nfi = zb.h2nd->GetBinContent(1,1);
      double nfq = zb.h2nd->GetBinContent(2,1) / nfi;
      double nfg = zb.h2nd->GetBinContent(3,1) / nfi;
      double nfc = zb.h2nd->GetBinContent(4,1) / nfi;
      double nfb = zb.h2nd->GetBinContent(5,1) / nfi;
      double nfie = zb.h2nd->GetBinError(1,1);
      double nfqe = zb.h2nd->GetBinError(2,1) / nfi;
      double nfge = zb.h2nd->GetBinError(3,1) / nfi;
      double nfce = zb.h2nd->GetBinError(4,1) / nfi;
      double nfbe = zb.h2nd->GetBinError(5,1) / nfi;
      nfie /= nfi;
      nfi  /= nfi;

      // Inferred counts per flavor in data (from fit)
      double nid  = zb.h2nd->GetBinContent(1,1);
      double nide = zb.h2nd->GetBinError(1,1);
      double nqd  = zb.h2nd->GetBinContent(1,2) / nid;
      double nqde = zb.h2nd->GetBinError(1,2) / nid;
      double ngd  = zb.h2nd->GetBinContent(1,3) / nid;
      double ngde = zb.h2nd->GetBinError(1,3) / nid;
      double ncd  = zb.h2nd->GetBinContent(1,4) / nid;
      double ncde = zb.h2nd->GetBinError(1,4) / nid;
      double nbd  = zb.h2nd->GetBinContent(1,5) / nid;
      double nbde = zb.h2nd->GetBinError(1,5) / nid;
      nide /= nid;
      nid  /= nid;
      // ...and in original MC
      double nim  = zb.h2nm->GetBinContent(1,1);
      double nime = zb.h2nm->GetBinError(1,1);
      double nqm  = zb.h2nm->GetBinContent(1,2) / nim;
      double nqme = zb.h2nm->GetBinError(1,2) / nim;
      double ngm  = zb.h2nm->GetBinContent(1,3) / nim;
      double ngme = zb.h2nm->GetBinError(1,3) / nim;
      double ncm  = zb.h2nm->GetBinContent(1,4) / nim;
      double ncme = zb.h2nm->GetBinError(1,4) / nim;
      double nbm  = zb.h2nm->GetBinContent(1,5) / nim;
      double nbme = zb.h2nm->GetBinError(1,5) / nim;
      nime /= nim;
      nim  /= nim;
      // ...plus data/MC ratio for counts
      double ni  = nid / nim;
      double nie = sqrt(pow(nide/nid,2)+pow(nime/nim,2));
      double nq  = nqd / nqm;
      double nqe = sqrt(pow(nqde/nqd,2)+pow(nqme/nqm,2));
      double ng  = ngd / ngm;
      double nge = sqrt(pow(ngde/ngd,2)+pow(ngme/ngm,2));
      double nc  = ncd / ncm;
      double nce = sqrt(pow(ncde/ncd,2)+pow(ncme/ncm,2));
      double nb  = nbd / nbm;
      double nbe = sqrt(pow(nbde/nbd,2)+pow(nbme/nbm,2));

      // Inferred responses per flavor in data (from fit)
      double rid  = zb.h2rd->GetBinContent(1,1);
      double ride = zb.h2rd->GetBinError(1,1);
      double rqd  = zb.h2rd->GetBinContent(1,2);
      double rqde = zb.h2rd->GetBinError(1,2);
      double rgd  = zb.h2rd->GetBinContent(1,3);
      double rgde = zb.h2rd->GetBinError(1,3);
      double rcd  = zb.h2rd->GetBinContent(1,4);
      double rcde = zb.h2rd->GetBinError(1,4);
      double rbd  = zb.h2rd->GetBinContent(1,5);
      double rbde = zb.h2rd->GetBinError(1,5);
      // ...in original MC
      double rim  = zb.h2rm->GetBinContent(1,1);
      double rime = zb.h2rm->GetBinError(1,1);
      double rqm  = zb.h2rm->GetBinContent(1,2);
      double rqme = zb.h2rm->GetBinError(1,2);
      double rgm  = zb.h2rm->GetBinContent(1,3);
      double rgme = zb.h2rm->GetBinError(1,3);
      double rcm  = zb.h2rm->GetBinContent(1,4);
      double rcme = zb.h2rm->GetBinError(1,4);
      double rbm  = zb.h2rm->GetBinContent(1,5);
      double rbme = zb.h2rm->GetBinError(1,5);
      // ...and in MC double tags for reference and uncertainties
      //double rim  = zb.h2rm->GetBinContent(1,1);
      //double rime = zb.h2rm->GetBinError(1,1);
      double rqqm  = zb.h2rm->GetBinContent(2,2);
      double rqqme = zb.h2rm->GetBinError(2,2);
      double rggm  = zb.h2rm->GetBinContent(3,3);
      double rggme = zb.h2rm->GetBinError(3,3);
      double rgqm  = zb.h2rm->GetBinContent(2,3);
      double rgqme = zb.h2rm->GetBinError(2,3);
      double rccm  = zb.h2rm->GetBinContent(4,4);
      double rccme = zb.h2rm->GetBinError(4,4);
      double rbbm  = zb.h2rm->GetBinContent(5,5);
      double rbbme = zb.h2rm->GetBinError(5,5);
      // ...plus data/MC ratio
      double ri  = rid/rim;
      double rie = ri*sqrt(pow(ride/rid,2)+pow(rime/rim,2));

      // Flavor response data/MC ratios directly from fit
      // Uncertainties need some more thought, on how to incorporate
      // statistical fluctuations from reference MC side
      double rq  = zb.rq;
      double rqe = sqrt(pow(zb.rqe,2) + pow(rqqme/rqqm,2));
      double rg  = zb.rg;
      double rge = sqrt(pow(zb.rge,2) + pow(rggme/rggm,2) +
			+ pow(rgqme/rgqm,2));
      double rc  = zb.rc;
      double rce = sqrt(pow(zb.rce,2) + pow(rccme/rccm,2));
      double rb  = zb.rb;
      double rbe = sqrt(pow(zb.rbe,2) + pow(rbbme/rbbm,2));
      
      double pt = hz->GetBinCenter(i);
      //if (rbde<0.10 && fabs(rb-rq)>-0.06 && fabs(rb-rq)<0.04) {
      {
	int n = grdb->GetN();

	gmi->SetPoint(n, pt, nmi/ndi);
	gmi->SetPointError(n, 0, nmi/ndi*sqrt(pow(nmie/nmi,2)+pow(ndie/ndi,2)));
	gmq->SetPoint(n, pt, nmq/ndq);
	gmq->SetPointError(n, 0, nmq/ndq*sqrt(pow(nmqe/nmq,2)+pow(ndqe/ndq,2)));
	gmg->SetPoint(n, pt, nmg/ndg);
	gmg->SetPointError(n, 0, nmg/ndg*sqrt(pow(nmge/nmg,2)+pow(ndge/ndg,2)));
	gmc->SetPoint(n, pt, nmc/ndc);
	gmc->SetPointError(n, 0, nmc/ndc*sqrt(pow(nmce/nmc,2)+pow(ndce/ndc,2)));
	gmb->SetPoint(n, pt, nmb/ndb);
	gmb->SetPointError(n, 0, nmb/ndb*sqrt(pow(nmbe/nmb,2)+pow(ndbe/ndb,2)));

	gfi->SetPoint(n, pt, nfi/ndi);
	gfi->SetPointError(n, 0, nfi/ndi*sqrt(pow(nfie/nfi,2)+pow(ndie/ndi,2)));
	gfq->SetPoint(n, pt, nfq/ndq);
	gfq->SetPointError(n, 0, nfq/ndq*sqrt(pow(nfqe/nfq,2)+pow(ndqe/ndq,2)));
	gfg->SetPoint(n, pt, nfg/ndg);
	gfg->SetPointError(n, 0, nfg/ndg*sqrt(pow(nfge/nfg,2)+pow(ndge/ndg,2)));
	gfc->SetPoint(n, pt, nfc/ndc);
	gfc->SetPointError(n, 0, nfc/ndc*sqrt(pow(nfce/nfc,2)+pow(ndce/ndc,2)));
	gfb->SetPoint(n, pt, nfb/ndb);
	gfb->SetPointError(n, 0, nfb/ndb*sqrt(pow(nfbe/nfb,2)+pow(ndbe/ndb,2)));
	
	//double refd = rqd;
	double refd = (useIncRef ? rid : rqd);
	grdi->SetPoint(n, pt, 100.*(rid-refd));
	grdi->SetPointError(n, 0, 100.*ride);
	grdq->SetPoint(n, pt, 100.*(rqd-refd));
	grdq->SetPointError(n, 0, 100.*rqde);
	grdg->SetPoint(n, pt, 100.*(rgd-refd));
	//grdg->SetPointError(n, 0, 100.*rgde);
	grdg->SetPointError(n, 0, 100.*sqrt(pow(rggme/rggm,2)
					    +pow(rgqme/rgqm,2)));
	grdc->SetPoint(n, pt, 100.*(rcd-refd));
	grdc->SetPointError(n, 0, 100.*rcde);
	grdb->SetPoint(n, pt, 100.*(rbd-refd));
	grdb->SetPointError(n, 0, 100.*rbde);

	//double refm = rqm;
	double refm = (useIncRef ? rim : rqm);
	double k = 0.98; // shift pT for error bar visibility
	grmi->SetPoint(n, k*pt, 100.*(rim-refm));
	grmi->SetPointError(n, 0, 100.*rime);
	grmq->SetPoint(n, k*pt, 100.*(rqm-refm));
	grmq->SetPointError(n, 0, 100.*rqme);
	grmg->SetPoint(n, k*pt, 100.*(rgm-refm));
	//grmg->SetPointError(n, 0, 100.*rgde);
	grmg->SetPointError(n, 0, 100.*sqrt(pow(rggme/rggm,2)
					    +pow(rgqme/rgqm,2)));
	grmc->SetPoint(n, k*pt, 100.*(rcm-refm));
	grmc->SetPointError(n, 0, 100.*rcme);
	grmb->SetPoint(n, k*pt, 100.*(rbm-refm));
	grmb->SetPointError(n, 0, 100.*rbme);
	
	//double ref = rq;
	double ref = (useIncRef ? ri : rq);
	gri->SetPoint(n, pt, 100.*(ri-ref));
	gri->SetPointError(n, 0, 100.*rie);
	grq->SetPoint(n, pt, 100.*(rq-ref));
	grq->SetPointError(n, 0, 100.*rqe);
	grg->SetPoint(n, pt, 100.*(rg-ref));
	grg->SetPointError(n, 0, 100.*rge);
	grc->SetPoint(n, pt, 100.*(rc-ref));
	grc->SetPointError(n, 0, 100.*rce);
	grb->SetPoint(n, pt, 100.*(rb-ref));
	grb->SetPointError(n, 0, 100.*rbe);


	gndi->SetPoint(n, pt, 100.*nid);
	gndi->SetPointError(n, 0, 100.*nide);
	gndq->SetPoint(n, pt, 100.*nqd);
	gndq->SetPointError(n, 0, 100.*nqde);
	gndg->SetPoint(n, pt, 100.*ngd);
	gndg->SetPointError(n, 0, 100.*ngde);
	gndc->SetPoint(n, pt, 100.*ncd);
	gndc->SetPointError(n, 0, 100.*ncde);
	gndb->SetPoint(n, pt, 100.*nbd);
	gndb->SetPointError(n, 0, 100.*nbde);

	gni->SetPoint(n, pt, 100.*(nid/nim-1));
	gni->SetPointError(n, 0, 100.*nide/nid);
	gnq->SetPoint(n, pt, 100.*(nqd/nqm-1));
	gnq->SetPointError(n, 0, 100.*nqde/nqd);
	gng->SetPoint(n, pt, 100.*(ngd/ngm-1));
	gng->SetPointError(n, 0, 100.*ngde/ngd);
	gnc->SetPoint(n, pt, 100.*(ncd/ncm-1));
	gnc->SetPointError(n, 0, 100.*ncde/ncd);
	gnb->SetPoint(n, pt, 100.*(nbd/nbm-1));
	gnb->SetPointError(n, 0, 100.*nbde/nbd);
      }
    }
  } // for i

  const int nf = 5;
  int markers[nf] = {kFullSquare,kFullCircle,kFullCircle,
		     kFullCircle,kFullCircle};
  int markersm[nf] = {kOpenSquare,kOpenCircle,kOpenCircle,
		      kOpenCircle,kOpenCircle};
  int colors[nf] = {kBlack,kBlue,kOrange+2,kGreen+2,kRed};
  string labels[nf] = {"Z+jet","q-jet","g-jet","c-jet","b-jet"};
  TGraphErrors *gms[nf] = {gmi, gmq, gmg, gmc, gmb};
  TGraphErrors *gfs[nf] = {gfi, gfq, gfg, gfc, gfb};
  TGraphErrors *grds[nf] = {grdi, grdq, grdg, grdc, grdb};
  TGraphErrors *grms[nf] = {grmi, grmq, grmg, grmc, grmb};
  TGraphErrors *grs[nf] = {gri, grq, grg, grc, grb};
  TGraphErrors *gnds[nf] = {gndi, gndq, gndg, gndc, gndb};
  TGraphErrors *gns[nf] = {gni, gnq, gng, gnc, gnb};
  const int ic = 3;
  const int ib = 4;

  TH1D *hup = new TH1D("hup",";p_{T,Z} (GeV);Z+f - Z+q (%)",
		       2000-minpt,minpt,2000);
  if (useIncRef) hup->SetYTitle("Z+f - Z+jet (%)");
  hup->SetMinimum(-13);//-16);//-13);
  hup->SetMaximum(+6.5);
  hup->GetXaxis()->SetMoreLogLabels();
  hup->GetXaxis()->SetNoExponent();

  TH1D *hdw = new TH1D("hdw",";p_{T,Z} (GeV);Data-MC (%)",
		       2000-minpt,minpt,2000);
  hdw->SetMinimum(-6+1e-4);//-12+1e-4);//-3
  hdw->SetMaximum(+6-1e-4);//+4
  hdw->GetXaxis()->SetMoreLogLabels();
  hdw->GetXaxis()->SetNoExponent();

  TCanvas *c0 = tdrDiCanvas("c0b",hup,hdw,4,11);
  hup->GetYaxis()->SetTitleOffset(1.00);

  c0->cd(1);
  gPad->SetLogx();

  TLegend *leg(0);
  leg = tdrLeg(0.70,0.88-nf*0.05,0.90,0.88);
  //leg->SetHeader("Data");

  TLatex *tex = new TLatex();
  tex->SetTextSize(0.044); tex->SetNDC();
  tex->DrawLatex(0.61,0.37,"Anti-k_{T} R=0.4 CHS");
  tex->DrawLatex(0.61,0.31,"HDM 15 GeV");
  if (isEta25) tex->DrawLatex(0.61,0.25,"|#eta|<2.5, #alpha<1.0");
  else         tex->DrawLatex(0.61,0.25,"|#eta|<1.3, #alpha<1.0");

  c0->cd(2);
  gPad->SetLogx();

  TLegend *legb = tdrLeg(0.585,0.40,0.785,0.90);
  legb->SetTextSize(0.035*2);
  legb->SetHeader("Flavor response residual");

  for (int i = 0; i != nf; ++i) {

    c0->cd(1);

    TGraphErrors *grm = (TGraphErrors*)grms[i];
    TGraphErrors *grd = (TGraphErrors*)grds[i];
    tdrDraw(grm,"Pz",markersm[i],colors[i]);
    tdrDraw(grd,"Pz",markers[i],colors[i]);
    leg->AddEntry(grd,labels[i].c_str(),"PL");

    // Fit Z+g and Z+jet responses in MC
    if (labels[i]=="g-jet" || labels[i]=="Z+jet") {
      TF1 *f1 = new TF1(Form("f0_%d",i),
			"[2]+[0]*pow(0.01*x,[1])+[3]/(0.01*x)",
			fitminpt,250);
      f1->SetParameters(-1,-0.5,-1,+0.5);
      f1->FixParameter(3,0);
      grd->Fit(f1,"QRN");
      f1->SetLineColor(colors[i]);
      f1->Draw("SAME");
      
      TLatex *tf1 = new TLatex();
      tf1->SetNDC(); tf1->SetTextSize(0.03);
      tf1->SetTextColor(colors[i]);
      tf1->DrawLatex(0.61,labels[i]=="Z+jet" ? 0.53-0.03 : 0.37+0.06,
		     Form("(%+1.2f%+1.2f(p_{T}/100)^{%+1.3f})%%",
			  f1->GetParameter(2),f1->GetParameter(0),
			  f1->GetParameter(1)));
    }

    c0->cd(2);

    TGraphErrors *grr = (TGraphErrors*)grs[i];

    TF1 *f0 = new TF1(Form("fr0_%d",i),"[2]+[0]*pow(0.01*x,[1])",fitminpt,300);
    //f0->SetParameters(0.01,-0.5,0);
    f0->SetParameters(0.5,-0.5,0);
    f0->FixParameter(2,0);
    f0->SetParLimits(1,-1,0);
    if (i==ic) f0->FixParameter(1,0); // only fix C to constant
    if (i==ib) f0->FixParameter(1,0); // b as well
    //f0->FixParameter(1,0); // b as well
    //if (i==ii) f0->FixParameter(1,0); // and why not inclusive
    grr->Fit(f0,"QRN");
    f0->SetLineColor(colors[i]);
    //if (i!=iq) f0->Draw("SAME");
    if ((useIncRef && i!=0) || (!useIncRef && i!=1))
      f0->Draw("SAME");

    // De-stress inclusive and Z+c
    //if (i==ii || i==ic) {
    //for (int j = 0; j != gr->GetN() && i==ic; ++j)
    //	gr->SetPoint(j, gr->GetX()[j]*0.98, gr->GetY()[j]);
    //gr->SetMarkerSize(0.5);
    //f0b->SetLineStyle(kDashed);
    //}

    //if (i!=iq) tdrDraw(gr,"Pz",markers[i],colors[i]);
    if ((useIncRef && i!=0) || (!useIncRef && i!=1))
      tdrDraw(grr,"Pz",markers[i],colors[i]);

    //if (i!=iq)
    if ((useIncRef && i!=0) || (!useIncRef && i!=1))
      legb->AddEntry(grr,Form("%+1.2f#pm%1.2f%%,"
			      " #chi^{2}=%4.1f/%d",
			      f0->GetParameter(0),f0->GetParError(0),
			      f0->GetChisquare(),f0->GetNDF()),"PL");
    //if (i==ig) {
    //tex->SetTextColor(colors[i]);
    //tex->SetTextSize(0.03*2);
    //tex->DrawLatex(0.65,0.40-0.05,Form("(1-p_{0}(p_{T}/100)^{%1.3f})%%",
    //					 f0->GetParameter(1)));
    //}
  } // for i

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(minpt,0,300,0);

  gPad->Update();
  c0->SaveAs("pdf/hadW_Zb_ZfMinusZq.pdf");


  TH1D *hup3 = new TH1D("hup3",";p_{T,Z} (GeV);Flavor fraction (%)",
			2000-minpt,minpt,2000);
  hup3->SetMinimum(0);
  hup3->SetMaximum(110);
  hup3->GetXaxis()->SetMoreLogLabels();
  hup3->GetXaxis()->SetNoExponent();

  TH1D *hdw3 = new TH1D("hdw3",";p_{T,Z} (GeV);Data/MC (%)",
			2000-minpt,minpt,2000);
  hdw3->SetMinimum(-50);
  hdw3->SetMaximum(+50);
  hdw3->GetXaxis()->SetMoreLogLabels();
  hdw3->GetXaxis()->SetNoExponent();

  TCanvas *c3 = tdrDiCanvas("c3b",hup3,hdw3,4,11);
  hup3->GetYaxis()->SetTitleOffset(1.00);

  c3->cd(1);
  gPad->SetLogx();

  TLegend *leg3 = tdrLeg(0.70,0.88-nf*0.05,0.90,0.88);

  tex->SetTextSize(0.045);
  tex->DrawLatex(0.61,0.37,"Anti-k_{T} R=0.4 CHS");
  tex->DrawLatex(0.61,0.31,"HDM 15 GeV");
  if (isEta25) tex->DrawLatex(0.61,0.25,"|#eta|<2.5, #alpha<1.0");
  else         tex->DrawLatex(0.61,0.25,"|#eta|<1.3, #alpha<1.0");

  c3->cd(2);
  gPad->SetLogx();

  //TLegend *leg3b = tdrLeg(0.585,0.40,0.785,0.90);
  //leg3b->SetTextSize(0.035*2);
  //leg3b->SetHeader("Flavor response residual");

  l->SetLineStyle(kDashed);
  l->DrawLine(minpt,0,300,0);

  for (int i = 0; i != nf; ++i) {

    c3->cd(1);

    TGraphErrors *gnd = (TGraphErrors*)gnds[i];
    tdrDraw(gnd,"Pz",markers[i],colors[i]);
    leg3->AddEntry(gnd,labels[i].c_str(),"PL");

    c3->cd(2);

    TGraphErrors *gn = (TGraphErrors*)gns[i];
    tdrDraw(gn,"Pz",markers[i],colors[i]);
  } // for i


  gPad->Update();
  c3->SaveAs("pdf/hadW_Zb_FlavorFractions.pdf");

  TH1D *h4 = new TH1D("h4",";p_{T,Z} (GeV);MC / Data (Event fraction)",
		      2000-minpt,minpt,2000);
  h4->SetMinimum(0.8);//0.7);
  h4->SetMaximum(1.3);//1.5);
  h4->GetXaxis()->SetMoreLogLabels();
  h4->GetXaxis()->SetNoExponent();

  TCanvas *c4 = tdrCanvas("c4",h4,4,11,kSquare);
  //h4->GetYaxis()->SetTitleOffset(1.00);
  gPad->SetLogx();

  TLegend *leg4 = tdrLeg(0.70,0.88-nf*0.05,0.90,0.88);

  tex->SetTextSize(0.045);
  tex->DrawLatex(0.61,0.37,"Anti-k_{T} R=0.4 CHS");
  tex->DrawLatex(0.61,0.31,"HDM 15 GeV");
  if (isEta25) tex->DrawLatex(0.61,0.25,"|#eta|<2.5, #alpha<1.0");
  else         tex->DrawLatex(0.61,0.25,"|#eta|<1.3, #alpha<1.0");

  for (int i = 0; i != nf; ++i) {

    TGraphErrors *gm = (TGraphErrors*)gms[i];
    tdrDraw(gm,"Pz",kOpenCircle, i!=0 ? colors[i]-9 : colors[i]);
    if (i!=0) gm->SetMarkerSize(0.7);
    if (i!=0) gm->SetLineStyle(kDotted);
    TGraphErrors *gf = (TGraphErrors*)gfs[i];
    tdrDraw(gf,"Pz",kFullCircle,colors[i]);
    gf->SetMarkerSize(0.7);

    leg4->AddEntry(gm,labels[i].c_str(),"PL");
  } // for

  c4->SaveAs("pdf/hadW_Zb_counts.pdf");

} // for hadW_Zb

zbFit hadW_Zbs(double pt) {
 
  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TFile *f = new TFile("rootfiles/jecdata2018ABCD.root","READ");
  assert(f && !f->IsZombie());

  if (isEta25) f->cd("data/eta00-25");
  else         f->cd("data/eta00-13");
  TDirectory *fd = gDirectory;
  if (isEta25) f->cd("mc/eta00-25");
  else         f->cd("mc/eta00-13");
  TDirectory *fm = gDirectory;

  curdir->cd();

  // Copy B, C, Q, G into similar TH2D as with hadW_cs, for each pT bin in turn
  TH1D *hind = (TH1D*)fd->Get("counts_zjet_a100"); assert(hind);
  TH1D *hqnd = (TH1D*)fd->Get("counts_zq_a100");   assert(hqnd);
  TH1D *hgnd = (TH1D*)fd->Get("counts_zg_a100");   assert(hgnd);
  TH1D *hcnd = (TH1D*)fd->Get("counts_zc_a100");   assert(hcnd);
  TH1D *hbnd = (TH1D*)fd->Get("counts_zb_a100");   assert(hbnd);
  TH1D *hond = (TH1D*)hind->Clone("counts_zo_a100_dt0");
  hond->Add(hqnd,-1);   hond->Add(hgnd,-1);
  hond->Add(hcnd,-1);   hond->Add(hbnd,-1);

  TH1D *hinm = (TH1D*)fm->Get("counts_zjet_a100"); assert(hinm);
  TH1D *hqnm = (TH1D*)fm->Get("counts_zq_a100");   assert(hqnm);
  TH1D *hgnm = (TH1D*)fm->Get("counts_zg_a100");   assert(hgnm);
  TH1D *hcnm = (TH1D*)fm->Get("counts_zc_a100");   assert(hcnm);
  TH1D *hbnm = (TH1D*)fm->Get("counts_zb_a100");   assert(hbnm);
  TH1D *honm = (TH1D*)hinm->Clone("counts_zo_a100_mc");
  honm->Add(hqnm,-1);   honm->Add(hgnm,-1);
  honm->Add(hcnm,-1);   honm->Add(hbnm,-1);

  const int nf = 5;
  const char *vf[nf+1] = {"i","q","g","c","b","o"};
  TH1D* vhnd[nf+1];
  TH1D* vhnm[nf+1];
  TH1D* ahnm[nf+1][nf+1];
  TH1D* vhmd[nf+1];
  TH1D* vhmm[nf+1];
  TH1D* ahmm[nf+1][nf+1];
  for (int i = 0; i != nf; ++i) {
    if (debug) cout << "Reading in z"<<vf[i]<<endl<<flush;
    vhnd[i] = (TH1D*)fd->Get(Form("counts_z%s_a100",vf[i]));
    assert(vhnd[i]);
    vhnm[i] = (TH1D*)fm->Get(Form("counts_z%s_a100",vf[i]));
    assert(vhnm[i]);
    //
    vhmd[i] = (TH1D*)fd->Get(Form("hdm_mpfchs1_z%s",vf[i]));
    assert(vhmd[i]);
    vhmm[i] = (TH1D*)fm->Get(Form("hdm_mpfchs1_z%s",vf[i]));
    assert(vhmm[i]);

    for (int j = 0; j != nf; ++j) {
      if (debug) cout << "Reading in z"<<vf[i]<<vf[j]<<endl<<flush;
      ahnm[i][j] = (TH1D*)fm->Get(Form("counts_z%s%s_a100",vf[i],vf[j]));
      //ahnm[j][i] = (TH1D*)fm->Get(Form("counts_z%s%s_a100",vf[i],vf[j]));
      assert(ahnm[i][j]);
      ahmm[i][j] = (TH1D*)fm->Get(Form("hdm_mpfchs1_z%s%s",vf[i],vf[j]));
      //ahmm[j][i] = (TH1D*)fm->Get(Form("hdm_mpfchs1_z%s%s",vf[i],vf[j]));
      assert(ahmm[i][j]);
    } // for j
  } // for i
  vhnd[nf] = (TH1D*)vhnd[0]->Clone("counts_zo_a100_dt");
  vhnm[nf] = (TH1D*)vhnm[0]->Clone("counts_zo_a100_mc");
  vhmd[nf] = (TH1D*)vhmd[0]->Clone("hdm_mpfchs1_zo_dt");
  vhmm[nf] = (TH1D*)vhmm[0]->Clone("hdm_mpfchs1_zo_mc");
  for (int i = 1; i != nf; ++i) {
    vhnd[nf]->Add(vhnd[i],-1);
    vhnm[nf]->Add(vhnm[i],-1);
  } // for i
  for (int i = 0; i != nf; ++i) {
    if (debug) cout << "Processing o for " << vf[i] << endl << flush;
    ahnm[nf][i] = (TH1D*)ahnm[0][i]->Clone(Form("counts_zo%s_a100_mc",vf[i]));
    ahnm[i][nf] = (TH1D*)ahnm[i][0]->Clone(Form("counts_z%so_a100_mc",vf[i]));
    if (i==0) ahnm[nf][nf] = (TH1D*)ahnm[0][0]->Clone("counts_zoo_a100_mc");
    //
    ahmm[nf][i] = (TH1D*)ahmm[0][i]->Clone(Form("hdm_zo%s_a100_mc",vf[i]));
    ahmm[i][nf] = (TH1D*)ahmm[i][0]->Clone(Form("hdm_z%so_a100_mc",vf[i]));
    if (i==0) ahmm[nf][nf] = (TH1D*)ahmm[0][0]->Clone("hmd_zoo_a100_mc");
    for (int j = 1; j != nf; ++j) {
      ahnm[nf][i]->Add(ahnm[j][i],-1);
      ahnm[i][nf]->Add(ahnm[i][j],-1);
      ahnm[nf][nf]->Add(ahnm[j][i],-1);
      if (j!=i) ahnm[nf][nf]->Add(ahnm[i][j],-1);
    } // for j
  } // for i


  TH1D *himd = (TH1D*)fd->Get("hdm_mpfchs1_zjet"); assert(himd);
  TH1D *hqmd = (TH1D*)fd->Get("hdm_mpfchs1_zq");   assert(hqmd);
  TH1D *hgmd = (TH1D*)fd->Get("hdm_mpfchs1_zg");   assert(hgmd);
  TH1D *hcmd = (TH1D*)fd->Get("hdm_mpfchs1_zc");   assert(hcmd);
  TH1D *hbmd = (TH1D*)fd->Get("hdm_mpfchs1_zb");   assert(hbmd);

  TH1D *himm = (TH1D*)fm->Get("hdm_mpfchs1_zjet"); assert(himm);
  TH1D *hqmm = (TH1D*)fm->Get("hdm_mpfchs1_zq");   assert(hqmm);
  TH1D *hgmm = (TH1D*)fm->Get("hdm_mpfchs1_zg");   assert(hgmm);
  TH1D *hcmm = (TH1D*)fm->Get("hdm_mpfchs1_zc");   assert(hcmm);
  TH1D *hbmm = (TH1D*)fm->Get("hdm_mpfchs1_zb");   assert(hbmm);

  h2nd = new TH2D(Form("h2nd_%1.0f",pt),
		  ";Tag region;N_{events}",6,-0.5,5.5,6,-0.5,5.5);
  h2nm = new TH2D(Form("h2nm_%1.0f",pt),
		  ";Tag region;N_{events}",6,-0.5,5.5,6,-0.5,5.5);
  
  h2md = new TH2D(Form("h2md_%1.0f",pt),
		  ";Tag region;HDM response",6,-0.5,5.5,6,-0.5,5.5);
  h2mm = new TH2D(Form("h2mm_%1.0f",pt),
		  ";Tag region;HDM response",6,-0.5,5.5,6,-0.5,5.5);

  //int k = hind->FindBin(35.);
  //int k = hind->FindBin(65.);
  //int k = hind->FindBin(80.);
  //int k = hind->FindBin(95.);
  int k = hind->FindBin(pt);
  //int k = hind->FindBin(110.);
  int kmax = hind->FindBin(300.-1);
  //double pt = hind->GetBinCenter(k);
  double ptmin = hind->GetBinLowEdge(k);
  double ptmax = hind->GetBinLowEdge(k+1);
  bool isIntegral = false;
  TF1 *f1 = new TF1("f1","[0]",fitminpt,300);
  for (int i = 0; i != nf+1; ++i) {
    for (int j = 0; j != nf+1; ++j) {
      double errnd(0), errnm(0);
      assert(vhnd[i]);
      if (j==0) {
	h2nd->SetBinContent(i+1, j+1, vhnd[i]->GetBinContent(k));
	h2nd->SetBinError(i+1, j+1, vhnd[i]->GetBinError(k));
	if (isIntegral) {
	  h2nd->SetBinContent(i+1,j+1,vhnd[i]->IntegralAndError(1,kmax,errnd));
	  h2nd->SetBinError(i+1,j+1,errnd);
	}
      }
      assert(vhmd[i]);
      if (j==0) {
	h2md->SetBinContent(i+1, j+1, vhmd[i]->GetBinContent(k));
	h2md->SetBinError(i+1, j+1, vhmd[i]->GetBinError(k));
	if (isIntegral) {
	  vhmd[i]->Fit(f1,"QRN");
	  h2md->SetBinContent(i+1, j+1, f1->GetParameter(0));
	  h2md->SetBinError(i+1, j+1, f1->GetParError(0));
	}
      }
      //
      assert(vhnm[i]);
      if (j==0) {
	h2nm->SetBinContent(i+1, j+1, vhnm[i]->GetBinContent(k));
	h2nm->SetBinError(i+1, j+1, vhnm[i]->GetBinError(k));
	if (isIntegral) {
	  h2nm->SetBinContent(i+1,j+1,vhnm[i]->IntegralAndError(1,kmax,errnm));
	  h2nm->SetBinError(i+1, j+1, errnm);
	}
      }
      assert(ahnm[i][j]);
      if (j!=0) {
	h2nm->SetBinContent(i+1, j+1, ahnm[i][j]->GetBinContent(k));
	h2nm->SetBinError(i+1, j+1, ahnm[i][j]->GetBinError(k));
	if (isIntegral) {
	  h2nm->SetBinContent(i+1,j+1,ahnm[i][j]->IntegralAndError(1,kmax,errnm));
	  h2nm->SetBinError(i+1, j+1, errnm);
	}
      }
      assert(ahmm[i][j]);
      //if (j==0 && false) {
      if (j==0) {
	h2mm->SetBinContent(i+1, j+1, vhmm[i]->GetBinContent(k));
	h2mm->SetBinError(i+1, j+1, vhmm[i]->GetBinError(k));
	if (isIntegral) {
	  vhmm[i]->Fit(f1,"QRN");
	  h2mm->SetBinContent(i+1, j+1, f1->GetParameter(0));
	  h2mm->SetBinError(i+1, j+1, f1->GetParError(0));
	}
      }
      //if (j!=0 || j==0) {
      if (j!=0) {
	h2mm->SetBinContent(i+1, j+1, ahmm[i][j]->GetBinContent(k));
	h2mm->SetBinError(i+1, j+1, ahmm[i][j]->GetBinError(k));
	if (isIntegral) {
	  ahmm[i][j]->Fit(f1,"QRN");
	  h2mm->SetBinContent(i+1, j+1, f1->GetParameter(0));
	  h2mm->SetBinError(i+1, j+1, f1->GetParError(0));
	}
      }
    } // for j
  } // for i

  // Check internal consistency by replacing data with MC
  if (testMC) {
    h2nd = h2nm;
    h2md = h2mm;
  }

  // Map TH2D to TH1D for fitting, include parameter constraints if preferred
  const int nx = h2nd->GetNbinsX();
  const int np = 7;//4;
  TF1 *fN = new TF1(Form("fN_%1.0f",ptmin),fnABCD,0.5-np,nx+0.5,np);
  // LD>0.5 for gluons 0.34 MC, 0.24 data
  // (1-0.34)/(1-0.24)-1 = 0.13
  // 0.34/0.24-1 = 0.42
  // something around 0.25 close to data, if gluon xsec modified 0.25
  //fN->SetParameters(0.25,0.12,0.15, 0.25,0.2,0.2,0.3); // initial guess
  //fN->SetParameters(0.13,0.12,0.15, 0.15,0.15,0.15,0.3); // initial guess
  fN->SetParameters(0.1745,0.12,0.15, 0.0,0.,0.0,0.3); // initial guess
  if (testMC || fixParsMCN) {
    fN->SetParameters(0,0,0, 0,0,0,0);
  }
  //if (fixParsNeff) {
  //fN->FixParameter(0,0.13);//0.25); // gluon eff
  //fN->FixParameter(1,0.12);
  //fN->FixParameter(2,0.15);
  //}
  {//if (fixParNgeff) {
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/QuarkGluonLikelihood#Systematics
    // Pythia8 gluon (|eta| < 2.0, inclusive pT) = 2.5626*x^3 - 3.2240*x^2 + 1.8687*x + 0.6770,
    // where x is the value of the qgl output (i.e. central value of the QGL bin).
    TF1 *fqgl = new TF1(Form("fqll_%1.0f",ptmin),
		       "2.5626*x^3 - 3.2240*x^2 + 1.8687*x + 0.6770",0,1);
    // Use rootfiles/output-MC-1-UL17V4_BCDEF.root/Standard/Eta_0.0-1.3/mc
    // and TH2D hqgl2_g_g for estimating efficiency scale factor
    TFile *f = new TFile("rootfiles/output-MC-1-UL17V4_BCDEF.root","READ");
    assert(f && !f->IsZombie());
    f->cd("Standard/Eta_0.0-1.3/mc");
    TH2D *h2qgl = (TH2D*)gDirectory->Get("hqgl2_g_g"); assert(h2qgl);
    int i1 = h2qgl->GetXaxis()->FindBin(ptmin);
    int i2 = h2qgl->GetXaxis()->FindBin(ptmax-1);
    TH1D *hqgl = h2qgl->ProjectionY(Form("hqgl_%1.0f",ptmin),i1,i2);
    TH1D *hqglw = (TH1D*)hqgl->Clone(Form("hqglw_%1.0f",ptmin));
    hqglw->Multiply(fqgl);
    int i3 = hqgl->FindBin(0.5-1e-5);
    double effm = hqgl->Integral(1,i3) / hqgl->Integral();
    double effd = hqglw->Integral(1,i3)/hqglw->Integral();
    double effsf = 1 - effd / effm;
    fN->SetParameter(0, effsf); // 0.1745 for 85-105 bin
    if (testMC) fN->SetParameter(0,0.);
    else if (fixParNgeff) fN->FixParameter(0, effsf);
  }
  {//if (fixParNbeff) {
    // https://twiki.cern.ch/twiki/bin/view/CMSPublic/BTV13TeV2017DeepJet
    fN->SetParameter(1,0.12);
    if (testMC) fN->SetParameter(1,0.);
    else if (fixParNbeff) fN->FixParameter(1,0.12);
  }
  {//if (fixParNceff) {
    // https://twiki.cern.ch/twiki/bin/view/CMSPublic/BTV13TeV2017DeepJet
    fN->SetParameter(2,0.15);
    if (testMC) fN->SetParameter(2,0.);
    else if (fixParNceff) fN->FixParameter(2,0.15);
  }
  {//if (fixParNg) {
    // Any reference on this? QGL assumes Q/G ratio from MC...
    fN->SetParameter(3, 0.00);
    if (testMC) fN->SetParameter(3,0.);
    else if (fixParNg) fN->FixParameter(3, 0.00);
    // But, QGL did not separate b/c from uds, which are also low
    //fN->FixParameter(3, 0.30);
    //fN->FixParameter(3, 0.45);
    if (invertNgNeff && !testMC) {
      double keff = fN->GetParameter(0);
      double kn = fN->GetParameter(3);
      fN->ReleaseParameter(0);
      fN->ReleaseParameter(3);
      fN->SetParameter(0, kn);
      fN->SetParameter(3,keff);
      if (fixParNgeff) fN->FixParameter(3, keff);
      if (fixParNg)    fN->FixParameter(0, kn);
    }
  }
  {//if (fixParNb) {
    // https://cms.cern.ch/iCMS/analysisadmin/cadilines?line=SMP-20-015
    // Figure 3c
    // Z+b is is 20-30% high for MG5_aMC[NLO,FxFx]
    // Z+b is is about right for MG5_aMC[LO,MLM] => use this 
    // MLM used for Z+jet JEC productions
    //fN->FixParameter(4, 0.25);
    fN->SetParameter(4, 0.);//0.20);
    if (testMC) fN->SetParameter(4,0.);
    else if (fixParNb) fN->FixParameter(4, 0.);//0.20);
  }
  {//if (fixParNc) {
    // https://cms.cern.ch/iCMS/analysisadmin/cadilines?line=SMP-19-011
    // Figure 6b
    // Z+c is 20% high in MC for MG5_aMC+PY8 (<=2j NLO + PS)
    // Z+c is 5%  low  in MC for MG5_aMC+PY8 (<=4j LO + PS) => use this (0%)
    // MLM used for Z+jet JEC productions
    fN->SetParameter(5, 0);//0.20);
    if (testMC) fN->SetParameter(5,0.);
    else if (fixParNc) fN->FixParameter(5, 0);//0.20);
  }
  {//if (fixParNo) {
    fN->SetParameter(6, 0.30);
    if (testMC) fN->SetParameter(6,0.);
    else if (fixParNo) fN->FixParameter(6, 0.30);
  }
  /*
  if (fixParsNall) {
    fN->FixParameter(0,0.13);//0.25); // gluon eff
    fN->FixParameter(1,0.12);
    fN->FixParameter(2,0.15);
    fN->FixParameter(3,0.15);//0.25); // gluon xsec
    fN->FixParameter(4,0.15);//0.2);
    fN->FixParameter(5,0.15);//2);
    fN->FixParameter(6,0.3);
  }
  */
  if (fixParsMCN) {
    fN->FixParameter(0,0.); // gluon eff
    fN->FixParameter(1,0.);
    fN->FixParameter(2,0.);
    fN->FixParameter(3,0.); // gluon xsec
    fN->FixParameter(4,0.);
    fN->FixParameter(5,0.);
    fN->FixParameter(6,0.);
  }
  

  TH1D *h1nd = new TH1D(Form("h1nd_%1.0f",pt),"",nx+np,0.5-np,nx+0.5);
  // Set 5% uncertainty on parameters around starting value
  for (int i=0; i != np; ++i) {
    h1nd->SetBinContent(h1nd->FindBin(-i),fN->GetParameter(i));
    h1nd->SetBinError(h1nd->FindBin(-i),constrainParsN ? 0.05 : 10);
  }
  for (int i = 1; i != h2nd->GetNbinsX()+1; ++i) {
    int j = 1;
    int k = nx*(j-1) + (i-1) + 1;
    assert(k<nx+1+np);
    int kk = h1nd->FindBin(k);
    h1nd->SetBinContent(kk, h2nd->GetBinContent(i,j));
    h1nd->SetBinError(kk, h2nd->GetBinError(i,j));
  } // for i
  h2nf = (TH2D*)h2nd->Clone("h2nf");
  fN->Eval(1); // set bins to default without fitting
  h1nd->Fit(fN,debug ? "RN" : "QRN");

  TH1D *h1 = tdrHist(Form("h1_%1.0f",ptmin),
		     "Fraction of events",2e-5,10-1e-3,//1.25,
		     "Tag region",-0.5,5.5);
  h1->GetXaxis()->SetBinLabel(1,"All");
  h1->GetXaxis()->SetBinLabel(2,"QGL>0.5");//"Quark");
  h1->GetXaxis()->SetBinLabel(3,"QGL<0.5");//Gluon");
  h1->GetXaxis()->SetBinLabel(4,"DeepCT");//"Charm");
  h1->GetXaxis()->SetBinLabel(5,"DeepBT");//"Bottom");
  h1->GetXaxis()->SetBinLabel(6,"None");
  TH1D *h1d = (TH1D*)h1->Clone(Form("h1d_%1.0f",ptmin));
  if (showNDcomp) h1d->GetYaxis()->SetRangeUser(0.00,1.50); // composition also
  else            h1d->GetYaxis()->SetRangeUser(0.80,1.30); // only inclusive
  h1d->SetYTitle("MC / Data");

  //lumi_13TeV = "2018, 59.9 fb^{-1}";
  //lumi_13TeV = "2017, 41.5 fb^{-1}";
  lumi_13TeV = "Z+jet, UL17+18, 101.4 fb^{-1}";
  //TCanvas *c1 = tdrCanvas(Form("c1_%1.0f",ptmin),h1,4,11,kSquare);
  TCanvas *c1 = tdrDiCanvas(Form("c1_%1.0f",ptmin),h1,h1d,4,11);

  c1->cd(1);
  gPad->SetLogy();

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  tex->DrawLatex(0.40,0.86,Form("%1.0f<p_{T}<%1.0f GeV",ptmin,ptmax));
  if (isEta25) tex->DrawLatex(0.40,0.80,"|#eta|<2.5, #alpha<1.0");
  else         tex->DrawLatex(0.40,0.80,"|#eta|<1.3, #alpha<1.0");

  TLegend *leg1 = tdrLeg(0.70,0.90-6*0.05,1.00,0.90);

  c1->cd(2);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(-0.5,1,5.5,1);
 
  //int color[] = {kGray+2,kRed,kGreen+2,kBlue,kOrange+1,kCyan+2};
  int color[] = {kGray+2,kBlue,kOrange+2,kGreen+2,kRed,kCyan+2};
  int ncolor = sizeof(color)/sizeof(color[0]);
  assert(ncolor==h2nm->GetNbinsY());

  const char* name[] = {"MC","q","g","c","b","others"};
  int nname = sizeof(name)/sizeof(name[0]);
  assert(nname==h2nm->GetNbinsY());

  c1->cd(1);

  // Data
  double nd  = h2nd->GetBinContent(1,1);
  double nde = h2nd->GetBinError(1,1);
  double nm  = h2nm->GetBinContent(1,1);
  double nme = h2nm->GetBinError(1,1);
  TH1D *hnd = h2nd->ProjectionX("hnd",1,1);
  hnd->Scale(1./nd);
  //hnd->Scale(1./nm); // >50% off?
  tdrDraw(hnd,"HP",kFullCircle,kBlack,kSolid,-1,kNone);
  leg1->AddEntry(hnd,testMC ? "MC'" : "Data","PL");

  c1->cd(2);

  TH1D *hndr = (TH1D*)hnd->Clone(Form("hndr_%1.0f",ptmin));
  hndr->Divide(hnd);
  tdrDraw(hndr,"HP",kFullCircle,kBlack,kSolid,-1,kNone);
  hndr->SetMarkerSize(0.5);

  // MC and fit
  assert(nm!=0);
  for (int j = 1; j != h2nm->GetNbinsY()+1; ++j) {

    c1->cd(1);

    TH1D *hnm = h2nm->ProjectionX(Form("hnm%d",j),j,j);
    hnm->Scale(nm!=0 ? 1./nm : 1);
    //hnm->SetLineWidth(2);
    //hnm->SetBinContent(1,hnm->GetBinContent(1)/4.);
    //tdrDraw(hnm,"HP",j==1 ? kOpenCircle : kNone,color[j-1],kSolid,-1,kNone);
    tdrDraw(hnm,"HP",j==1 ? kOpenCircle : kNone,color[j-1],kDotted,-1,kNone);

    //leg1->AddEntry(hnm,name[j-1],"PL");
    if (j==1) leg1->AddEntry(hnm,name[j-1],"PL");

    if (!h2nf) continue;
    TH1D *hnf = h2nf->ProjectionX(Form("hnf%d",j),j,j);
    hnf->Scale(nd!=0 ? 1./nd : 1);
    //hnf->SetLineWidth(2);
    hnf->SetMarkerSize(0.5);
    //tdrDraw(hnf,"HP",j==1 ? kOpenCircle : kNone,color[j-1],kDotted,-1,kNone);
    tdrDraw(hnf,"HP",j==1 ? kOpenCircle : kNone,color[j-1],kSolid,-1,kNone);

    if (j!=1) leg1->AddEntry(hnf,name[j-1],"PL");

    c1->cd(2);

    TH1D *hnr = (TH1D*)hnm->Clone(Form("hnr_%d_%1.0f",j,ptmin));
    hnr->Divide(hnd);
    //tdrDraw(hnr,"HP",j==1 ? kOpenCircle : kNone,color[j-1],kSolid,-1,kNone);
    if (showNDcomp || j==1)
      tdrDraw(hnr,"HP",j==1 ? kOpenCircle : kNone,color[j-1],kDotted,-1,kNone);

    TH1D *hnfr = (TH1D*)hnf->Clone(Form("hnfr_%d_%1.0f",j,ptmin));
    hnfr->Divide(hnd);
    //tdrDraw(hnfr,"HP",j==1 ? kOpenCircle : kNone,color[j-1],kDotted,-1,kNone);
    if (showNDcomp || j==1)
      tdrDraw(hnfr,"HP",j==1 ? kOpenCircle : kNone,color[j-1],kSolid,-1,kNone);
  } // for j

  gPad->RedrawAxis();

  //gPad->SetLogy(kFALSE);
  //h2f->Draw("BOX");
  //h1d->Draw();
  //p2f = p2d->ProjectionXY("p2f"); // TH2D to be filled
  //h2mf = h2md; //
  h2mf = (TH2D*)h2md->Clone("h2mf");
  TH1D *hmd = h2md->ProjectionX("p1d",1,1); // TH1D to be fit for data
  //TH1D *hmm = h2mm->ProjectionX("p1m",1,1); // MC for adding uncertainty
  // Add MC uncertainty to data points so fit is not overconstrained
  //for (int i = 1; i != hmd->GetNbinsX()+1; ++i) {
  //if (hmd->GetBinContent(i)!=0) {
  //  hmd->SetBinError(i,sqrt(pow(hmd->GetBinError(i)/hmd->GetBinContent(i),2)+
  //			      pow(hmm->GetBinError(i)/hmm->GetBinContent(i),2))*
  //		       hmd->GetBinContent(i));
  //}
  //}
  // Fit full range [-0.5,5.5] with inclusive and other, or only [0.5,4.5]?
  //TF1 *fr = new TF1(Form("fr_%1.0f",ptmin),frABCD,-0.5,4.5,4); // Fit function
  //const int nf = 5; // q,g,c,b, o [not i]
  const int ns = nf*(nf+1); // also i in reco regions
  const int npr = nf + ns;
  //TF1 *fr = new TF1(Form("fr_%1.0f",ptmin),frABCD,-0.5,5.5,5); // Fit function
  TF1 *fr = new TF1(Form("fr_%1.0f",ptmin),frABCD,-ns-0.5,nf+0.5,npr);
  // seems much too sensitive to including "Other"
  //fr->SetParameters(1,1,1,1); h2nf = h2nm; // MC reproduction
  //fr->SetParameters(1.05,0.98,1,1,1,1);
  fr->SetParameters(1.005,0.980,1,1,1); // q,g,c,b,o
  for (int i = nf; i != np; ++i) fr->SetParameter(i, 0);
  TH1D *h1md = new TH1D(Form("h1md_%1.0f",ptmin),"",
			npr+1,-ns-0.5,nf+0.5);
  for (int i = 1; i != h1md->GetNbinsX()+1; ++i) {
    if (h1md->GetBinCenter(i)>-0.5) {
      int j = hmd->FindBin(h1md->GetBinCenter(i));
      h1md->SetBinContent(i, hmd->GetBinContent(j));
      h1md->SetBinError(i, hmd->GetBinError(j));
    }
    else {
      h1md->SetBinContent(i, 0.);
      h1md->SetBinError(i, 1.);
    }
  }
  //fr->FixParameter(4,1); // others
  //fr->FixParameter(5,1); // others
  //fr->SetParLimits(0,0.85,1.25);
  //fr->SetParLimits(1,0.85,1.25);
  //fr->SetParLimits(2,0.85,1.25);
  //fr->SetParLimits(3,0.85,1.25);
  //fr->SetParLimits(5,0.85,1.25);
  //double kd = h2md->GetBinContent(1,1)/h2mm->GetBinContent(1,1);
  //double kg = 0.5*(0.994+1.000);
  //fr->SetParameters(kd,kd,kd,kd*kg,kd);
  //fr->FixParameter(3,kd);
  fr->Eval(0); // Fill in p2f
  //h2f = h2m; // Use MC purities instead of data
  hmd->Fit(fr,debug ? "RN" : "QRN");
  //hmm->Fit(fr,"RN"); // test fit to MC

  c1->cd(2);

  tex->SetTextSize(2*0.030);
  tex->DrawLatex(0.17,0.60,Form("%1.2f#pm%1.2f",nm/nd,
				nm/nd*sqrt(pow(nde/nd,2)+pow(nme/nm,2)))); 

  TH1D *h2 = tdrHist(Form("h2_%1.0f",ptmin),
		     "HDM response",0.85,1.15,//0.90,1.05,
		     "Tag region",-0.5,5.5);
  h2->GetXaxis()->SetBinLabel(1,"All");
  h2->GetXaxis()->SetBinLabel(2,"QGL>0.5");//"Quark");
  h2->GetXaxis()->SetBinLabel(3,"QGL<0.5");//"Gluon");
  h2->GetXaxis()->SetBinLabel(4,"DeepCT");//"Charm");
  h2->GetXaxis()->SetBinLabel(5,"DeepBT");//"Bottom");
  h2->GetXaxis()->SetBinLabel(6,"Others");
  
  TCanvas *c2 = tdrCanvas(Form("c2_%1.0f",ptmin),h2,4,11,kSquare);

  tex->SetTextSize(0.045);
  tex->DrawLatex(0.65,0.87,Form("%1.0f<p_{T}<%1.0f GeV",ptmin,ptmax));
  if (isEta25) tex->DrawLatex(0.65,0.82,"|#eta|<2.5, #alpha<1.0");
  else         tex->DrawLatex(0.65,0.82,"|#eta|<1.3, #alpha<1.0");

  TLegend *leg2 = tdrLeg(0.40,0.90-6*0.045,0.70,0.90);

  // Data
  double md = h2md->GetBinContent(1,1);
  double mm = h2mm->GetBinContent(1,1);
  TH1D *hmd2 = h2md->ProjectionX("hd2",1,1);
  hmd2->Scale(1./md);
  //hmd2->Scale(1./mm);
  //hmd2->Scale(1./hd2->GetBinContent(1));
  tdrDraw(hmd2,"HP",kFullCircle,kBlack,kSolid,-1,kNone);
  leg2->AddEntry(hmd2,testMC ? "MC'" : "Data","PL");


  // MC and fit
  for (int j = 1; j != h2mm->GetNbinsY()+1; ++j) {

    TH1D *hmm = h2mm->ProjectionX(Form("hmm%d2",j),j,j);
    hmm->Scale(1./mm);
    tdrDraw(hmm,"HP",j==1 ? kOpenCircle : kNone,color[j-1],kDotted,-1,kNone);

    leg2->AddEntry(hmm,name[j-1],"PL");

    if (!h2mf) continue;
    TH1D *hmf = h2mf->ProjectionX(Form("hmf%d2",j),j,j);
    hmf->Scale(1./md);
    //hmf->Scale(1./mm);
    //hmf->SetLineWidth(2);
    hmf->SetMarkerSize(0.5);
    tdrDraw(hmf,"HP",j==1 ? kOpenCircle : kNone,color[j-1],kSolid,-1,kNone);
  } // for j

  gPad->RedrawAxis();

  zbFit zb;

  zb.nd  = nd;
  zb.nde = nde;
  zb.nm  = nm;
  zb.nme = nme;

  zb.h2nd = h2nf;
  zb.h2nm = h2nm;

  zb.nchi2 = fN->GetChisquare();
  zb.nndf = fN->GetNDF();

  zb.h2rd = h2mf;
  zb.h2rm = h2mm;

  zb.rq  = fr->GetParameter(0);
  zb.rqe = fr->GetParError(0);
  zb.rg  = fr->GetParameter(1);
  zb.rge = fr->GetParError(1);
  zb.rc  = fr->GetParameter(2);
  zb.rce = fr->GetParError(2);
  zb.rb  = fr->GetParameter(3);
  zb.rbe = fr->GetParError(3);
  zb.ro  = fr->GetParameter(5);
  zb.roe = fr->GetParError(5);

  zb.rchi2 = fr->GetChisquare();
  zb.rndf = fr->GetNDF();

  //c1->SaveAs("pdf/hadW_Zb_fractions.pdf");
  //c2->SaveAs("pdf/hadW_Zb_responses.pdf");
  if ((ptmin>=30 && ptmax<=35) ||
      (ptmin>=85 && ptmax<=105) ||
      (ptmin>=60 && ptmax<=70)) {
    c1->SaveAs(Form("pdf/hadW_Zb_fractions_pt_%1.0f_%1.0f.pdf",ptmin,ptmax));
    c2->SaveAs(Form("pdf/hadW_Zb_responses_pt_%1.0f_%1.0f.pdf",ptmin,ptmax));
  }

  return zb;
} // hadW_Zbs


// Fit to mass-based response in ABCD regions
// Allow others response to differ in ABCD vs D
// by estimating different gluon fraction
Double_t frABCD(Double_t *x, Double_t *p) {

  assert(h2mm); // masses in MC
  assert(h2nf); // tag fractions (fitted) in data
  assert(h2mf); // masses fit to data

  const int nf = 5;
  int k = int(*x)+1;
  //assert(k>=1);
  assert(k<=nf+1);
  assert(h2nf->GetNbinsX()==nf+1);
  assert(h2nf->GetNbinsY()==nf+1);

  // Create full response matrix for easier plotting later
  // although this is a bit of extra work
  for (int i = 1; i != h2nf->GetNbinsX()+1; ++i) { // reco-tag

    double sumw(0), sumrw(0), sumre2w(0), sumwe2(0);
    for (int j = 2; j != h2nf->GetNbinsY()+1; ++j) { // gen-flavor

      double w = h2nf->GetBinContent(i,j);
      double ew = h2nf->GetBinError(i,j);
      double r0  = h2mm->GetBinContent(i,j);
      double re0 = h2mm->GetBinError(i,j);
      int is = nf + (i-1)*nf + (j-2); 
      //double r1 = r0 + p[is]*re0;
      //double k = (j>5 ? 1 : p[j-2]);
      double rj = p[j-2];
      //double r  = rj * r1;
      double r  = rj * r0 + p[is]*re0;
      double re = rj * re0;
      // Linearize parameters for more stable solution?
      //double r  = (rj-1) + r1;
      //double re = rj * re0;
      // Scale response of all flavors in None category
      //if (i==6) r  *= p[5];
      //if (i==6) re *= p[5];

      sumw    += w;
      sumrw   += r * w;
      sumre2w += re * re * w;
      sumwe2 += ew * ew;
      h2mf->SetBinContent(i, j, r);
      h2mf->SetBinError(i, j, re);
    } // for j>1

    // non-tagged true flavor response
    h2mf->SetBinContent(i, 1, sumrw / sumw);
    h2mf->SetBinError(i, 1, sqrt(sumre2w / sumw));
    //h2mf->SetBinError(i, 1, sqrt(sumwe2)/sumw);
    //h2mf->SetBinError(i, 1, sqrt(sumre2w / sumw + sumwe2 / sumw));
  } // for i

  // Return nuisance parameter
  if (k<0) {
    return p[(-k)+nf];
  }

  return (h2mf->GetBinContent(k,1));
}

// Fit to event numbers in ABCD
Double_t fnABCD(Double_t *x, Double_t *p) {

  int k = int(*x);
  
  const int nfl = 6; // should be 6, but...
  if (k<=0) {
    //assert(-k<4);
    assert(-k<7);
    return p[-k];
  }
  // int k = nx*(j-1) + (i-1) + 1;
  int nx = h2nd->GetNbinsX();
  int i0 = ((k-1) % nx) + 1;
  int j0 = ((k-1) / nx) + 1;

  assert(h2nm);
  assert(h2nd);
  assert(h2nf);

  // First redistribute MC events to reco QGCBO for each of gen QGCBO
  for (int i = 2; i != h2nm->GetNbinsX()+1; ++i) {
    for (int j = 2; j != h2nm->GetNbinsY()+1; ++j) {

      // reco QGCBO
      bool isQ = (i==2); int iQ = 2;
      bool isG = (i==3); int iG = 3;
      bool isC = (i==4); int iC = 4;
      bool isB = (i==5); int iB = 5;
      bool isO = (i==6); int iO = 6;

      // gen QGCBO
      bool jsQ = (j==2); int jQ = 2;
      bool jsG = (j==3); int jG = 3;
      bool jsC = (j==4); int jC = 4;
      bool jsB = (j==5); int jB = 5;
      bool jsO = (j==6); int jO = 6;

      //double n(1);
      double n  = max(0., h2nm->GetBinContent(i, j));
      double ne = max(0., h2nm->GetBinError  (i, j));

      //if (jsQ) {
      //double nQ = h2nm->GetBinContent(1, jQ);
	//double pQ = h2nm->GetBinContent(iQ,jQ) / nQ;
	//double pG = h2nm->GetBinContent(iG,jQ) / nQ;
	//double pB = h2nm->GetBinContent(iB,jQ) / nQ;
	//double pC = h2nm->GetBinContent(iC,jQ) / nQ;
	//double pO = h2nm->GetBinContent(iO,jQ) / nQ;
	//double pX = h2nm->GetBinContent(i ,jQ) / nQ;
	//n = pX * nQ;
      //}

      // Gluon tagging efficiency scale factor, estimated e.g. from
      // https://cds.cern.ch/record/2254861/files/1512294_73-78.pdf?#page=5
      // Fig.7a, LD>0.5 MC 0.24 data 0.34, so LD<0.5 MC 0.76 data 0.66
      // kGG ~ 1-0.66/0.766 = 0.13
      // Also change gluon cross section
      if (jsG) {
	double nG = h2nm->GetBinContent(1, jG);
	double pQ = h2nm->GetBinContent(iQ,jG) / nG;
	double pG = h2nm->GetBinContent(iG,jG) / nG;
	double pB = h2nm->GetBinContent(iB,jG) / nG;
	double pC = h2nm->GetBinContent(iC,jG) / nG;
	double pO = h2nm->GetBinContent(iO,jG) / nG;
	double pX = h2nm->GetBinContent(i ,jG) / nG;
	double kGG = p[0]; // reco G relative decrease for gen G
	double kG  = p[3]; // gluon fraction
	nG = (1-kG)*nG;
	// Lost reco gluon tags reappear as reco light quark tags
	// Reduced true gluon fraction scaled later
	if (isQ)      n = pX*nG + kGG*pG*nG;
	if (isG)      n = (1-kGG)*pX*nG;
	if (isC||isB) n = pX*nG;
	if (isO)      n = max(0.,pX*nG);

	double epG = h2nm->GetBinError(iG,jG) / nG;
	double epX = h2nm->GetBinError(i ,jG) / nG;
	// Lost reco gluon tags reappear as reco light quark tags
	// Reduced true gluon fraction scaled later
	if (isQ)      ne = sqrt(pow(epX*nG,2) + pow(kGG*epG*nG,2));
	if (isG)      ne = (1-kGG)*epX*nG;
	if (isC||isB) ne = epX*nG;
	if (isO)      ne = max(0.,epX*nG);
      }

      // B-tagging efficiency scale factors around 0.88
      // https://twiki.cern.ch/twiki/bin/view/CMSPublic/BTV13TeV2017DeepJet
      // Also change B cross section
      if (jsB) {
	double nB = h2nm->GetBinContent(1, jB);
	double pQ = h2nm->GetBinContent(iQ,jB) / nB;
	double pG = h2nm->GetBinContent(iG,jB) / nB;
	double pB = h2nm->GetBinContent(iB,jB) / nB;
	double pC = h2nm->GetBinContent(iC,jB) / nB;
	double pO = h2nm->GetBinContent(iO,jB) / nB;
	double pX = h2nm->GetBinContent(i ,jB) / nB;
	double kBB = p[1]; // reco B relative decrease for gen B
	double kB  = p[4]; // decrease of b jet cross section
	nB *= (1-kB);
	// Lost b-tags reappear as quarks, gluons and charm
	if (isQ||isG||isC) n = pX*nB + kBB*pB*nB*pX/(pQ+pG+pC);
	if (isB)           n = (1-kBB)*pX*nB;
	if (isO)           n = max(0.,pX*nB);


	double epB = h2nm->GetBinError(iB,jB) / nB;
	double epX = h2nm->GetBinError(i ,jB) / nB;
	// Lost b-tags reappear as quarks, gluons and charm
	if (isQ||isG||isC) ne = sqrt(pow(epX*nB,2) +
				     pow(kBB*epB*nB*pX/(pQ+pG+pC),2) +
				     pow(kBB*pB*nB*epX/(pQ+pG+pC),2));
	if (isB)           ne = (1-kBB)*epX*nB;
	if (isO)           ne = max(0.,epX*nB);
      }

      // C-tagging efficiency scale factors around 0.85
      // https://twiki.cern.ch/twiki/bin/view/CMSPublic/BTV13TeV2017DeepJet
      // Also change C cross section
      if (jsC) {
	double nC = h2nm->GetBinContent(1, jC);
	double pQ = h2nm->GetBinContent(iQ,jC) / nC;
	double pG = h2nm->GetBinContent(iG,jC) / nC;
	double pB = h2nm->GetBinContent(iB,jC) / nC;
	double pC = h2nm->GetBinContent(iC,jC) / nC;
	double pO = h2nm->GetBinContent(iO,jC) / nC;
	double pX = h2nm->GetBinContent(i ,jC) / nC;
	double kCC = p[2]; // reco C relative decrease for gen C
	double kC = p[5];
	nC = (1-kC)*nC;
	// Lost c-tags reappear as quarks and gluons
	if (isQ||isG) n = pX*nC + kCC*pC*nC*pX/(pQ+pG);
	if (isC)      n = (1-kCC)*pX*nC;
	if (isB)      n = pX*nC;
	if (isO)      n = max(0.,pX*nC);

	double epC = h2nm->GetBinError(iC,jC) / nC;
	double epX = h2nm->GetBinError(i ,jC) / nC;
	// Lost c-tags reappear as quarks and gluons
	if (isQ||isG) ne = sqrt(pow(epX*nC,2) +
				pow(kCC*epC*nC*pX/(pQ+pG),2) +
				pow(kCC*pC*nC*epX/(pQ+pG),2));
	if (isC)      ne = (1-kCC)*epX*nC;
	if (isB)      ne = epX*nC;
	if (isO)      ne = max(0.,epX*nC);
      }

      if (jsO) {
	//double nO = h2nm->GetBinContent(1, jO);
	//double pQ = h2nm->GetBinContent(iQ,jO) / nO;
	//double pG = h2nm->GetBinContent(iG,jO) / nO;
	//double pB = h2nm->GetBinContent(iB,jO) / nO;
	//double pC = h2nm->GetBinContent(iC,jO) / nO;
	//double pO = h2nm->GetBinContent(iO,jO) / nO;
	//double pX = h2nm->GetBinContent(i ,jO) / nO;
	//double kO = p[6];
	//n = (nO!=0 ? pX * nO : 0.);
	//n = max(0.,(1+kO)*n);
	n = max(0.,n);
	ne = max(0., ne);
      }

      // Increased none rate of all flavors
      // changes cross section slightly up
      double kXO = p[6];
      if (isO) n  = (1+kXO)*n;
      if (isO) ne = (1+kXO)*ne;

      h2nf->SetBinContent(i,j,n);
      h2nf->SetBinError(i,j,ne);
    } // for j
  } // for i


  double nd = h2nd->GetBinContent(1,1);
  //double nfe(0);
  //double nf = h2nf->IntegralAndError(2,nfl,2,nfl,nfe);
  double nf = h2nf->Integral(2,nfl,2,nfl);
  //double nm = h2nm->GetBinContent(1,1);
  // not keeping MC cross section since g,b,c xsec are varied
  //if (!((nf/nm-1)<1e-3)) {
  //cout << "nf="<<nf<<" nm="<<nm<<" nd="<<nd<<endl<<flush;
  //assert((nf/nm-1)<1e-3);
  //}
  //double kd = nd/nm;
  double kd = nd/nf;

  // Normalize exclusive recoQGBCO x genQGBCO to data integral
  for (int i = 2; i != h2nf->GetNbinsX()+1; ++i) {
    for (int j = 2; j != h2nf->GetNbinsY()+1; ++j) {
      h2nf->SetBinContent(i,j,kd*h2nf->GetBinContent(i,j));
      h2nf->SetBinError(i,j,kd*h2nf->GetBinError(i,j));
    } // for j
  } // for i

  // Sum up reco QGBCO for exclusive reco QGBCO
  for (int i = 2; i != h2nf->GetNbinsX()+1; ++i) {
    //h2nf->SetBinContent(i,1,h2nf->Integral(i,i,2,nfl));
    double nfle(0);
    h2nf->SetBinContent(i,1,h2nf->IntegralAndError(i,i,2,nfl,nfle));
    h2nf->SetBinError(i,1,nfle);
  } // for i

  // Sum up reco QGBCO for exclusive gen QGBCO
  for (int j = 2; j != h2nf->GetNbinsY()+1; ++j) {
    //h2nf->SetBinContent(1,j,h2nf->Integral(2,nfl,j,j));
    double nfle(0);
    h2nf->SetBinContent(1,j,h2nf->IntegralAndError(2,nfl,j,j,nfle));
    h2nf->SetBinError(1,j,nfle);
  } // for j

  // Set total to data
  //h2f->SetBinContent(1,1,nd);
  //assert(h2nf->GetNbinsX()==nfl+1);
  //assert(h2nf->GetNbinsY()==nfl+1);
  //assert(h2nf->GetNbinsX()==nfl);
  //assert(h2nf->GetNbinsY()==nfl);
  //h2nf->SetBinContent(1,1,h2nf->Integral(2,nfl,2,nfl));
  double nfle(0);
  h2nf->SetBinContent(1,1,h2nf->IntegralAndError(2,nfl,2,nfl,nfle));
  h2nf->SetBinError(1,1,nfle);  

  return (h2nf->GetBinContent(i0,j0));
} // fnABCD
