// Purpose: Rewrite drawZflavor.C and hadW_Zb.C for more efficient processing
//          Store event fractions in TH2D for easy access and manipulation
//          x-axis: pT x tag(i,uds,g,c,b)
//          y-axis: true flavor (i,uds,g,c,b) 
//
// Fundamental math:
// Nt = sum_f Etf * Nf,      R't = sum_f (Etf * Nf / Nt) * R'tf
// N = sum_t Nt = sum_f Nf,  R   = sum_t (Nt / N) * R't
// Writing wt = Nt / N, wf = Nf / N, we get...
// wt = sum_f Etf * wf,      
// 1 = sum_t wt = sum_f wf,  1 = sum_t Rt * wt
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TLine.h"
#include "TMatrixD.h"
//#include "TVectorD.h"

#include "../tdrstyle_mod15.C"

#include <vector>
#include <map>
#include <string>

TH2D *h2nm(0), *h2nd(0), *h2nf(0);
TH2D *h2rm(0), *h2rd(0), *h2rf(0);

// Notes: Smooth MC responses, maybe also event counts, to create smooth matrix.
//        This could stabilize inversion. Check condidion number?

bool drawFit = false;
double quarkScale = 1;//0.998;//1.001;//1.005;
double gluonScale = 1;//0.995;//0.970;
// Could check stability of quark scale with Qtag from W>qq'
double quarkScaleQtag = 1;//1.001;//1.005;
double quarkScaleGtag = 1;//0.997;
// Could gluon scale change in data as well?
double gluonScaleQtag = 1;//1.002;
double gluonScaleGtag = 1;//0.995;

// Model of gluon efficiency based on dijet MC and official shape SF
double gluonEffSF(double ptmin, double ptmax, double &ineffsf);
double quarkEffSF(double ptmin, double ptmax, double &ineffsf);
double heavyEffSF(double pt, double &mistag);

// 2D model of tagged flavor response over inclusive
double *_ptbins; int _npt;
Double_t fR(const Double_t *x, const Double_t *p) {

  //cout << (*x) << ", " << flush;

  // Decode flavor and pt
  // double x = i*log10(ptbins[npt]/pt0) + log10(pt/pt0);
  double pt0 = _ptbins[0];
  int itag = int((*x) / log10(_ptbins[_npt]/pt0));
  if (!(itag>=0 && itag<5)) return 0;
  assert(itag>=0);
  assert(itag<5);
  double pt = pt0*pow(10., (*x)-itag*log10(_ptbins[_npt]/pt0));
  if (!(pt>_ptbins[0] && pt<_ptbins[_npt])) return 0; // excludes nan also
  assert(pt>_ptbins[0] && pt<_ptbins[_npt]);

  double r = p[2*itag] + p[2*itag+1]*(pow(0.01*pt,-0.3)-1);

  //cout << r << ", " << flush;

  return max(0.85,min(1.15,r));
}

// 2D model of tagged flavor fraction over tag
const int _npw = 4;
Double_t fW(const Double_t *x, const Double_t *p) {

  double pt0 = _ptbins[0];
  int itag = int((*x) / log10(_ptbins[_npt]/pt0));
  if (!(itag>=0 && itag<5)) return 0;
  assert(itag>=0);
  assert(itag<5);
  double pt = pt0*pow(10., (*x)-itag*log10(_ptbins[_npt]/pt0));
  if (!(pt>_ptbins[0] && pt<_ptbins[_npt])) return 0; // excludes nan also
  assert(pt>_ptbins[0] && pt<_ptbins[_npt]);

  double w(0);
  // Fixed power law and offset term
  if (_npw==3) {
    w = p[3*itag] + p[3*itag+1]*(pow(0.01*pt,-0.3)-1)
      + p[3*itag+2]*(1./(0.01*pt)-1);
  }
  // Free power law and offset term => not really improving
  if (_npw==4 && false) {
    w = p[4*itag] + p[4*itag+1]*(pow(0.01*pt,p[4*itag+2])-1)
      + p[4*itag+3]*(1./(0.01*pt)-1);
  }
  // Fixed power and log offset => does help especially q and g
  if (_npw==4) {
    w = p[4*itag] + p[4*itag+1]*(pow(0.01*pt,-0.3)-1)
      + p[4*itag+2]*(1./(0.01*pt)-1)
      + p[4*itag+3]*log(0.01*pt)/(0.01*pt);
  }

  //double w = p[3*itag] + p[3*itag+1]*log(0.01*pt)
  //+ p[3*itag+2]*(1./(0.01*pt)-1);

  return w;
}



void Zflavor() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  
  // pT range can be adjusted here to avoid empty bins
  double ptbins[] =
  // Region with sensible gluon behavior (Rq>Rg for q-tag)
    {40, 45, 50, 60, 70, 85, 105, 130, 175, 230, 300};
  // Avoid non-physics region
    //{30, 35, 40, 45, 50, 60, 70, 85, 105, 130, 175, 230, 300};
  // Full range
    //{15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 85, 105, 130, 175, 230, 300};
  //400, 500, 700};
  const int npt = sizeof(ptbins)/sizeof(ptbins[0])-1;
  _ptbins = &ptbins[0]; // warning: ptbins destroyed at the end, should clone it
  _npt = npt;

  // Ordering of tag and flavor bins can be changed here
  // Naming scheme corresponds to rootfiles/jecdata[X].root
  const int nt = 5;
  string tagbins[] = {"i","q","g","c","b"};
  // Flavor bins separately from tag to enable differentiating ud, s,or other
  const int nf = 5;
  string flvbins[] = {"i","q","g","c","b"};

  // Map tags and pT bins on a single x-axis for easy handling and plotting
  vector<double> vx;
  for (int i = 0; i != nt; ++i) {
    for (int j = 0; j != npt; ++j) {
      double pt0 = ptbins[0];
      double x = i*log10(ptbins[npt]/pt0) + log10(ptbins[j]/pt0);
      vx.push_back(x);
    } // for j
  } // for i
  vx.push_back(nt*log10(ptbins[npt]/ptbins[0]));

  // Axis labels to be used on plots
  map<string,const char*> taglabel;
  taglabel["i"] = "Inclusive";
  taglabel["q"] = "q tag";
  taglabel["g"] = "g tag";
  taglabel["c"] = "c tag";
  taglabel["b"] = "b tag";
  // Legen labels to be used on plots
  map<string,const char*> flvlabel;
  flvlabel["i"] = "Any jet";
  flvlabel["q"] = "uds jet";
  flvlabel["g"] = "g jet";
  flvlabel["c"] = "c jet";
  flvlabel["b"] = "b jet";
  
  // Reference 2D histogram
  //TH2D *h2ref = new TH2D("h2ref",";Tag region and p_{T} bin;True flavor",
  TH2D *h2ref = new TH2D("h2ref",";Global bin(tag, p_{T});True flavor",
			 vx.size()-1,&vx[0], nf,-0.5,nf-0.5);
  for (int i = 0; i != nt; ++i)
    h2ref->GetXaxis()->SetBinLabel(int((i+0.5)*npt), taglabel[tagbins[i]]);
  for (int i = 0; i != nf; ++i)
    h2ref->GetYaxis()->SetBinLabel(i+1, flvlabel[flvbins[i]]);
  h2ref->GetXaxis()->SetTitleOffset(1.5);
  h2ref->GetZaxis()->SetTitleOffset(1.3);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);

  TFile *f = new TFile("rootfiles/jecdata2018ABCD.root","READ");
  assert(f && !f->IsZombie());
  f->cd("mc/eta00-25");
  TDirectory *dm = gDirectory;
  f->cd("data/eta00-25");
  TDirectory *dd = gDirectory;

  TH2D *h2n = (TH2D*)h2ref->Clone("h2n"); // event counts
  TH2D *h2r = (TH2D*)h2ref->Clone("h2r"); // responses
  //
  TH2D *h2nd = (TH2D*)h2ref->Clone("h2nd"); // event counts
  TH2D *h2rd = (TH2D*)h2ref->Clone("h2rd"); // responses
  
  // Load results from jecdata.root
  for (int i = 0; i != nt; ++i) {
    for (int j = 0; j != nf; ++j) {

      // Get input histograms
      const char *ct = tagbins[i].c_str();
      const char *cf = flvbins[j].c_str();
      
      // MC simulation
      TH1D *hn = (TH1D*)dm->Get(Form("counts_z%s%s_a100",ct,cf));
      assert(hn);
      TH1D *hr = (TH1D*)dm->Get(Form("hdm_mpfchs1_z%s%s",ct,cf));
      assert(hr);

      // Data
      TH1D *hnd = (TH1D*)dd->Get(Form("counts_z%s%s_a100",ct,""));
      assert(hnd);
      // Careful: Data z[t][f] are clones of MC, z[t] not
      TH1D *hrd = (TH1D*)dd->Get(Form("hdm_mpfchs1_z%s%s",ct,""));
      assert(hrd);
      
      // Copy results to 2D histograms with new x-axis
      for (int ipt = 1; ipt != hn->GetNbinsX()+1; ++ipt) {

	// Check valid pT range
	double pt = hn->GetBinCenter(ipt);
	if (pt>ptbins[0] && pt<ptbins[npt]) {

	  // Map to new axis
	  double pt0 = ptbins[0];
	  double dpt = hn->GetBinWidth(ipt);
	  double x = i*log10(ptbins[npt]/pt0) + log10(pt/pt0);
	  double ix = h2ref->GetXaxis()->FindBin(x);
	  double iy = h2ref->GetYaxis()->FindBin(j);
	  assert(h2ref->GetBinContent(ix,iy)==0);

	  h2n->SetBinContent(ix,iy, hn->GetBinContent(ipt) / dpt);
	  h2n->SetBinError(ix,iy, hn->GetBinError(ipt) / dpt);

	  h2r->SetBinContent(ix,iy, hr->GetBinContent(ipt));
	  h2r->SetBinError(ix,iy, hr->GetBinError(ipt));

	  h2nd->SetBinContent(ix,iy, hnd->GetBinContent(ipt) / dpt);
	  h2nd->SetBinError(ix,iy, hnd->GetBinError(ipt) / dpt);

	  h2rd->SetBinContent(ix,iy, hrd->GetBinContent(ipt));
	  h2rd->SetBinError(ix,iy, hrd->GetBinError(ipt));
	}
      } // for ipt
    } // for j
  } // for i

  // PATCH 'no' tags by recreating inclusive from exclusive flavor tags
  for (int ipt = 0; ipt != npt; ++ipt) {

    int i0 = 0*npt + ipt + 1;
    double sumwf(0), sumwfr(0);
    for (int iq = 1; iq != nf; ++iq) {
      int j = iq+1;
      double sumwt(0), sumwtr(0);
      for (int it = 1; it != nt; ++it) {
	int i = it*npt + ipt + 1;
	sumwt += h2n->GetBinContent(i,j);
	sumwf += h2n->GetBinContent(i,j);
	sumwtr += h2n->GetBinContent(i,j) * h2r->GetBinContent(i,j);
	sumwfr += h2n->GetBinContent(i,j) * h2r->GetBinContent(i,j);
      } // for it
      // Exclusive flavor tags->inclusive flavor
      h2n->SetBinContent(i0,j,sumwt);
      h2r->SetBinContent(i0,j,sumwtr/sumwt);
    } // for iq

    // Check magnitude of change
    if (!(fabs(sumwf/h2n->GetBinContent(i0,1)-1)<1e-3))
      cout << "ipt="<<ipt<<" sumwf="<<sumwf
	   << " h2n="<<h2n->GetBinContent(i0,1)<<endl;

    // Exclusive flavor tags->fully inclusive Z+jet
    h2n->SetBinContent(i0,1,sumwf);
    h2r->SetBinContent(i0,1,sumwfr/sumwf);
    //assert(fabs(sumwf-1)<1e-3);
    //assert(fabs(sumwf-1)<1e-1);
  } // for ipt

  // For data, can only patch by summing up tagged regions
  // Sanity check MC on the same go
  for (int ipt = 0; ipt != npt; ++ipt) {

    int i0 = 0*npt + ipt + 1;
    int j = 1;
    double sumwt(0), sumwtr(0);
    double sumwtm(0), sumwtrm(0);
    for (int it = 1; it != nt; ++it) {
      int i = it*npt + ipt + 1;
      sumwt += h2nd->GetBinContent(i,j);
      sumwtr += h2nd->GetBinContent(i,j) * h2rd->GetBinContent(i,j);
      //
      sumwtm += h2n->GetBinContent(i,j);
      sumwtrm += h2n->GetBinContent(i,j) * h2r->GetBinContent(i,j);
    } // for it

    //if (!(h2nd->GetBinContent(i0,j)==sumwt))
    if (fabs(h2nd->GetBinContent(i0,j)/sumwt-1)>1e-3)
      cout << "ipt="<<ipt<< " h2nd=="<<h2nd->GetBinContent(i0,j)
	   <<" sumwt="<<sumwt<<endl<<flush;
    //assert(h2rd->GetBinContent(i0,j)==sumwtr/sumwt);
    //if (!(h2rd->GetBinContent(i0,j)==sumwtr/sumwt))
    if (fabs(h2rd->GetBinContent(i0,j)-sumwtr/sumwt)>1e-3)
      cout << "ipt="<<ipt<< " h2rd=="<<h2rd->GetBinContent(i0,j)
	   <<" <r>="<<sumwtr/sumwt<<endl<<flush;
    h2nd->SetBinContent(i0,j,sumwt);
    h2rd->SetBinContent(i0,j,sumwtr/sumwt);

    //assert(h2n->GetBinContent(i0,j)==sumwtm);
    //if (!(h2n->GetBinContent(i0,j)==sumwtm))
    if (fabs(h2n->GetBinContent(i0,j)/sumwtm-1)>1e-3)
      cout << "ipt="<<ipt<< " h2n="<<h2n->GetBinContent(i0,j)
	   <<" sumwtm="<<sumwtm<<endl<<flush;
    //assert(h2r->GetBinContent(i0,j)==sumwtrm/sumwtm);
    //if (!(h2r->GetBinContent(i0,j)==sumwtrm/sumwtm))
    if (fabs(h2r->GetBinContent(i0,j)-sumwtrm/sumwtm)>1e-3)
      cout << "ipt="<<ipt<< " h2r="<<h2r->GetBinContent(i0,j)
	   <<" <r>="<<sumwtrm/sumwtm<<endl<<flush;
  }

  // Plot original event counts
  h2n->Draw("COLZ");
  //h2n->GetZaxis()->SetRangeUser(10,3e6);
  h2n->GetZaxis()->SetRangeUser(0.1,3e6);
  //h2n->Scale(1./3e6);
  //h2n->GetZaxis()->SetRangeUser(1./3e5,1);
  gPad->SetLogz();
  gPad->SetRightMargin(0.20);
  gPad->SetBottomMargin(0.20);

  // Normalize event counts to inclusive sample or to tagged sample
  TH2D *h2wi = (TH2D*)h2ref->Clone("h2wi"); //h2wi->Reset();
  TH2D *h2wf = (TH2D*)h2ref->Clone("h2wf"); //h2wf->Reset();
  TH2D *h2wF = (TH2D*)h2ref->Clone("h2wF"); //h2wf->Reset();
  //
  TH2D *h2wid = (TH2D*)h2ref->Clone("h2wid"); //h2wid->Reset();
  TH2D *h2wfd = (TH2D*)h2ref->Clone("h2wfd"); //h2wfd->Reset();
  TH2D *h2wFd = (TH2D*)h2ref->Clone("h2wFd"); //h2wfd->Reset();
  for (int i = 1; i != h2ref->GetNbinsX()+1; ++i) {
    for (int j = 1; j != h2ref->GetNbinsY()+1; ++j) {
      double k = (i-1) % npt + 1;
      h2wi->SetBinContent(i,j,h2n->GetBinContent(i,j)/h2n->GetBinContent(k,1));
      h2wi->SetBinError(i,j,h2n->GetBinError(i,j)/h2n->GetBinContent(k,1));
      h2wf->SetBinContent(i,j,h2n->GetBinContent(i,j)/h2n->GetBinContent(i,1));
      h2wf->SetBinError(i,j,h2n->GetBinError(i,j)/h2n->GetBinContent(i,1));
      //
      h2wid->SetBinContent(i,j,h2nd->GetBinContent(i,j) /
			   h2nd->GetBinContent(k,1));
			   //h2n->GetBinContent(k,1));
      h2wid->SetBinError(i,j,h2nd->GetBinError(i,j) / 
			 h2nd->GetBinContent(k,1));
      h2wfd->SetBinContent(i,j,h2nd->GetBinContent(i,j) / 
			   h2nd->GetBinContent(i,1));
			   //h2n->GetBinContent(i,1));
      h2wfd->SetBinError(i,j,h2nd->GetBinError(i,j) / 
			 h2nd->GetBinContent(i,1));
      // scale c,b fractions to better represent their contribution to
      // HDM response differences (due to 5-10% of neutrinos)
      int t = (i-1)/npt+1;
      double F = 1;
      if (j==2 && t>=4) F = (t==5 ? 10 : 5);
      if (j==3 && t>=4) F = (t==5 ? 10 : 5);
      if (j==4 && t<4)  F = 5;
      if (j==4 && t==5) F = 2;
      if (j==5 && t<4)  F = 10;
      if (j==5 && t==4) F = 2;
      h2wF->SetBinContent(i,j,F*h2n->GetBinContent(i,j) / 
			  h2n->GetBinContent(i,1));
      h2wF->SetBinError(i,j,F*h2n->GetBinError(i,j) / 
			h2n->GetBinContent(i,1));
    } // for j
  } // for i
  
  // Normalize responses to inclusive sample or to tagged sample
  TH2D *h2ri = (TH2D*)h2ref->Clone("h2ri"); //h2ri->Reset();
  TH2D *h2rf = (TH2D*)h2ref->Clone("h2rf"); //h2rf->Reset();
  //
  TH2D *h2rid = (TH2D*)h2ref->Clone("h2rid"); //h2rid->Reset();
  TH2D *h2rfd = (TH2D*)h2ref->Clone("h2rfd"); //h2rfd->Reset();
  for (int i = 1; i != h2ref->GetNbinsX()+1; ++i) {
    for (int j = 1; j != h2ref->GetNbinsY()+1; ++j) {
      double k = (i-1) % npt + 1;
      if (h2r->GetBinContent(i,j)>0) { // filter out nan in one q bin
	h2ri->SetBinContent(i,j,h2r->GetBinContent(i,j) /
			    h2r->GetBinContent(k,1));
	h2ri->SetBinError(i,j,h2r->GetBinError(i,j)/h2r->GetBinContent(k,1));
	h2rf->SetBinContent(i,j,h2r->GetBinContent(i,j) /
			    h2r->GetBinContent(i,1));
	h2rf->SetBinError(i,j,h2r->GetBinError(i,j)/h2r->GetBinContent(i,1));

	h2rid->SetBinContent(i,j,h2rd->GetBinContent(i,j) /
			    h2rd->GetBinContent(k,1));
	h2rid->SetBinError(i,j,h2rd->GetBinError(i,j)/h2rd->GetBinContent(k,1));
	h2rfd->SetBinContent(i,j,h2rd->GetBinContent(i,j) /
			    h2rd->GetBinContent(i,1));
	h2rfd->SetBinError(i,j,h2rd->GetBinError(i,j)/h2rd->GetBinContent(i,1));
      }
    } // for j
  } // for i


  // Inclusive weights
  h2wi->Draw("COLZ");
  h2wi->GetZaxis()->SetRangeUser(1e-4,1.);
  h2wi->SetZTitle("Per inclusive weight");

  // Per flavor weights
  /*
  h2wf->Draw("COLZ");
  h2wf->GetZaxis()->SetRangeUser(1e-4,1.);
  h2wf->SetZTitle("Per flavor weight");
  */

  // Flavor responses per tag
  /*
  h2r->Draw("COLZ");
  h2r->GetZaxis()->SetRangeUser(0.80,1.20);
  gPad->SetLogz(kFALSE);
  */

  TH1D *href = h2wi->ProjectionX("href",1,1); href->Reset();
  for (int i = 1; i != href->GetNbinsX()+1; ++i)
    href->GetXaxis()->SetBinLabel(i, " ");
  for (int i = 0; i != nt; ++i)
    href->GetXaxis()->SetBinLabel(int((i+0.5)*npt), taglabel[tagbins[i]]);
  double xmin  = href->GetXaxis()->GetBinLowEdge(1)+1e-4;
  double xmin2 = href->GetXaxis()->GetBinLowEdge(npt+1)+1e-4;
  double dx = xmin2-xmin;
  double xmax = href->GetXaxis()->GetBinLowEdge(href->GetNbinsX()+1)-1e-4;

  lumi_13TeV = "[Z+jet] UL17+18, 101.4 fb^{-1}";
  TH1D *h0 = (TH1D*)href->Clone("h0");
  h0->SetYTitle("Data / MC");
  
  // Counts: data/MC and flavor/inclusive double ratio
  TCanvas *c0 = tdrCanvas("c0",h0,4,11,kSquare);
  h0->GetXaxis()->SetTitleOffset(1.2);//1.5);
  h0->SetMinimum(0.60);
  h0->SetMaximum(1.50); 
  gPad->SetBottomMargin(0.15);//0.20);
  
  TH1D *hnm = h2n->ProjectionX("hnm",1,1);
  TH1D *hnd = h2nd->ProjectionX("hnd",1,1);
  TH1D *hwm = h2wi->ProjectionX("hwim",1,1);
  TH1D *hwd = h2wid->ProjectionX("hwid",1,1);
  //TH1D *hrm = h2ri->ProjectionX("hrim",1,1);
  //TH1D *hrd = h2rid->ProjectionX("hrid",1,1);
  hnd->Divide(hnm);
  hwd->Divide(hwm);
  /*
  hrd->Add(hrm,-1);
  hrd->Divide(hrm);
  hrd->Scale(10);
  hrm->Divide(hrm);
  hrd->Add(hrm);
  tdrDraw(hrd,"Pz",kOpenSquare,kRed);
  */
  l->SetLineStyle(kDotted);
  l->DrawLine(xmin,1.1,xmax,1.1);
  l->DrawLine(xmin,0.9,xmax,0.9);
  l->SetLineStyle(kDashed);
  l->DrawLine(xmin,1,xmax,1);
  tdrDraw(hnd,"Pz",kOpenCircle);
  tdrDraw(hwd,"Pz",kFullCircle);

  TLegend *leg0 = tdrLeg(0.50,0.90-0.05*2,0.80,0.90);
  leg0->AddEntry(hnd,"Counts","PLE");
  leg0->AddEntry(hwd,"Counts / inclusive","PLE");

  // Draw tag region boundaries
  l->SetLineStyle(kSolid);
  l->DrawLine(xmin+1*dx,0.6,xmin+1*dx,1.3);
  l->DrawLine(xmin+2*dx,0.6,xmin+2*dx,1.3);
  l->DrawLine(xmin+3*dx,0.6,xmin+3*dx,1.3);
  l->DrawLine(xmin+4*dx,0.6,xmin+4*dx,1.3);

  // Fit b rate to estimate tagging efficiency vs pT
  TF1 *fb = new TF1("fb",fW,xmin,xmax,_npw*nt);
  hwd->Fit(fb,"QRN");
  fb->Draw("SAME");

  gPad->Update();
  c0->SaveAs("pdf/Zflavor_counts_v1.pdf");


  //lumi_13TeV = "UL17+18, 101.4 fb^{-1}";
  TH1D *h1 = (TH1D*)href->Clone("h1");

  // Flavor fractions vs tag (smoothing)
  TCanvas *c1 = tdrCanvas("c1b",h1,4,11,kSquare);
  h1->SetYTitle("Counts / tag");
  h1->GetXaxis()->SetTitleOffset(1.2);//1.5);
  h1->SetMinimum(0.0);
  h1->SetMaximum(1.30);//1.25); 
  gPad->SetBottomMargin(0.15);//0.20);

  TLegend *leg1d = tdrLeg(0.45,0.90-3*0.05,0.75,0.90);
  TLegend *leg1f = tdrLeg(0.76,0.90-(nt+1-3)*0.05,1.06,0.90);
  TLegend *leg1d2 = tdrLeg(0.45,0.90-3*0.05+0.01,0.75,0.90+0.01);
  TLegend *leg1f2 = tdrLeg(0.76,0.90-(nt+1-3)*0.05+0.01,1.05,0.90+0.01);

  // Create smoothed variant of flavor fractions vs Z+tag and Z+jet
  TH2D *h2wfs = (TH2D*)h2ref->Clone("h2wfs");
  TH2D *h2wis = (TH2D*)h2ref->Clone("h2wis");
  int color[nf] = {kBlack,kBlue,kOrange+2,kGreen+2,kRed};
  for (int iq = 0; iq != nf; ++iq) {
    int j = iq+1;
    //TH1D *hw = h2wi->ProjectionX(Form("hw%d",j),j,j);
    TH1D *hw = h2wf->ProjectionX(Form("hw%d",j),j,j);
    //TH1D *hw = h2wF->ProjectionX(Form("hw%d",j),j,j); // SCALED
    //tdrDraw(hw,"P",kNone,color[iq],kSolid,-1,kNone);
    tdrDraw(hw,"P",kNone,color[iq],kDotted,-1,kNone);
    
    // draw data/MC ratio of counts on top
    if (true && iq==0) { // now in separate Zflavor_counts.pdf
      TH1D *hwim = h2wi->ProjectionX(Form("hwim%d",j),j,j);
      TH1D *hwid = h2wid->ProjectionX(Form("hwid%d",j),j,j);
      //hwid->Divide(hwim);
      //if (iq==0) tdrDraw(hwid,"P",kOpenCircle,color[iq],kSolid,-1,kNone);
      //hwid->SetMarkerSize(0.5);
      if (iq==0) tdrDraw(hwim,"P",kOpenCircle,color[iq],kDotted,-1,kNone);
      if (iq==0) tdrDraw(hwid,"P",kFullCircle,color[iq],kSolid,-1,kNone);
      hwim->SetMarkerSize(0.4);
      hwid->SetMarkerSize(0.5);

      leg1d->AddEntry(hwid,"Data / incl.","PL");
      leg1d->AddEntry(hwim,"MC / incl.","PL");
      leg1d2->AddEntry(hwid," ","");
      leg1d2->AddEntry(hwim," ","");
    }

    if (iq!=0 && iq<2) {
      leg1d->AddEntry(hw,Form("%s / tag",flvlabel[flvbins[iq]]),"L");
      //leg1d2->AddEntry(hw0," ","L");
    }
    if (iq!=0 && iq>=2) {
      leg1f->AddEntry(hw,Form("%s",flvlabel[flvbins[iq]]),"L");
      //leg1f2->AddEntry(hw0," ","L");
    }
    
    //double xmin = hw->GetBinLowEdge(1)+1e-3;
    //double xmin2 = hw->GetBinLowEdge(npt+1)+1e-3;
    //double xmax = hw->GetBinLowEdge(hw->GetNbinsX()+1)-1e-3;
    TF1 *f1 = new TF1(Form("f1w_%d",j),fW,xmin,xmax,_npw*nt);
    f1->SetLineColor(color[iq]);
    // First free fit
    //f1->SetParameters(0.5,0,0, 0.5,0,0, 0.5,0,0, 0.5,0,0, 0.5,0,0);
    for (int ipar = 0; ipar != _npw*nt; ++ipar) {
      f1->SetParameter(ipar, ipar%_npw==0 ? 0.5 : 0);
      if (_npw==4 && ipar%_npw==1 && false) f1->SetParameter(ipar, -0.3);
    }
    // Reduce degrees of freedom when statistics low and no clear evidence
    if (_npw==4) {
      if (iq==4) f1->FixParameter(_npw*1+3,0.);
      if (iq==4) f1->FixParameter(_npw*2+3,0.);
      if (iq==4) f1->FixParameter(_npw*3+3,0.);
      if (iq==4) f1->FixParameter(_npw*4+3,0.); // b-in-b
      //
      if (iq==3) f1->FixParameter(_npw*1+3,0.);
      if (iq==3) f1->FixParameter(_npw*2+3,0.);
      if (iq==3) f1->FixParameter(_npw*3+3,0.); // c-in-c
      if (iq==3) f1->FixParameter(_npw*4+3,0.);
      //
      if (iq==2) f1->FixParameter(_npw*3+3,0.);
      if (iq==2) f1->FixParameter(_npw*4+3,0.);
      //
      if (iq==1) f1->FixParameter(_npw*3+3,0.);
      if (iq==1) f1->FixParameter(_npw*4+3,0.);
    }
    hw->Fit(f1,"QRN");
    f1->SetLineStyle(kDotted);
    if (drawFit) f1->DrawClone("SAME");
    cout << "j="<<j<<", chi2/NDF="<<f1->GetChisquare()<<"/"<<f1->GetNDF()<<endl;

    // Copy fitted values for smoothed 2D
    // Edges will need to be recalculated with smoother counts for consistency
    for (int i = 1; i != h2wfs->GetNbinsX()+1; ++i) {
      double txpt = h2wfs->GetXaxis()->GetBinCenter(i);
      //double sf = h2wf->GetBinContent(i,j) / h2wF->GetBinContent(i,j); // SCALED
      double sf = 1;
      h2wfs->SetBinContent(i,j,sf*f1->Eval(txpt));
      double fi = h2wi->GetBinContent(i,1);
      h2wis->SetBinContent(i,j,fi*sf*f1->Eval(txpt));
    } // for j

    // Refit only flavor-tagged regions, which matter for the matrix
    f1->SetRange(xmin2,xmax);
    hw->Fit(f1,"QRN");
    f1->SetLineStyle(kSolid);
    if (drawFit) f1->Draw("SAME");
    cout << "  *, chi2/NDF="<<f1->GetChisquare()<<"/"<<f1->GetNDF()+4<<endl;
  } // for i
  
  // 1) Ensure that exclusive flavor tag fractions add up to 100%
  for (int ipt = 0; ipt != npt; ++ipt) {
    for (int it = 0; it != nt; ++it) {

      int i = it*npt + ipt + 1;
      
      // Calculate fraction sum
      double sumw(0);
      for (int iq = 1; iq != nf; ++iq) {
	int j = iq+1;
	sumw += h2wfs->GetBinContent(i,j);
      } // for iq
      // Normalize fraction sum to 100%
      for (int iq = 1; iq != nf; ++iq) {
	int j = iq+1;
	h2wfs->SetBinContent(i,j, h2wfs->GetBinContent(i,j)/sumw);
	h2wis->SetBinContent(i,j, h2wis->GetBinContent(i,j)/sumw);
      } // for iq
    } // for it
  } // for ipt

  // 2) Ensure that exclusive flavor tags add up to inclusive flavor and Z+jet
  for (int ipt = 0; ipt != npt; ++ipt) {

    int i0 = 0*npt + ipt + 1;
    double sumwf(0);
    for (int iq = 1; iq != nf; ++iq) {
      int j = iq+1;
      double sumwt(0);
      for (int it = 1; it != nt; ++it) {
	int i = it*npt + ipt + 1;
	sumwt += h2wis->GetBinContent(i,j);
	sumwf += h2wis->GetBinContent(i,j);
      } // for it
      // Exclusive flavor tags->inclusive flavor
      h2wis->SetBinContent(i0,j,sumwt);
      h2wfs->SetBinContent(i0,j,sumwt);
    } // for iq
    // Exclusive flavor tags->fully inclusive Z+jet
    h2wis->SetBinContent(i0,1,sumwf);
    h2wfs->SetBinContent(i0,1,sumwf);
    //assert(fabs(sumwf-1)<1e-3);
    if (!(fabs(sumwf-1)<1e-3))
      cout << " sumwf-1="<<sumwf-1<<endl;
  } // for ipt

  // Draw smoother plots on top
  for (int iq = 0; iq != nf; ++iq) {
    int j = iq+1;
    TH1D *hw = h2wfs->ProjectionX(Form("hwfs%d",j),j,j);
    TH1D *hwi = h2wis->ProjectionX(Form("hwis%d",j),j,j);
    //tdrDraw(hw,"HIST",kNone,color[iq],kSolid,-1,kNone);   
    if (iq==0) tdrDraw(hwi,"HIST",kNone,color[iq],kDotted,-1,kNone);     
    tdrDraw(hw,"HIST",kNone,color[iq],kDotted,-1,kNone);    
  }

  // Draw tag region boundaries
  l->SetLineStyle(kSolid);
  l->DrawLine(xmin+1*dx,0.,xmin+1*dx,1.);
  l->DrawLine(xmin+2*dx,0.,xmin+2*dx,1.);
  l->DrawLine(xmin+3*dx,0.,xmin+3*dx,1.);
  l->DrawLine(xmin+4*dx,0.,xmin+4*dx,1.);

  gPad->Update();
  c1->SaveAs("pdf/Zflavor_mcfrac_v1.pdf");


  //TH1D *h2 = h2ri->ProjectionX("h2",1,1); h2->Reset();
  //for (int i = 1; i != h2->GetNbinsX()+1; ++i)
  //h2->GetXaxis()->SetBinLabel(i, " ");
  //for (int i = 0; i != nt; ++i)
  //h2->GetXaxis()->SetBinLabel(int((i+0.5)*npt), taglabel[tagbins[i]]);
  TH1D *h2 = (TH1D*)href->Clone("h2");

  // Response vs inclusive (smoothing)
  TCanvas *c2 = tdrCanvas("c2",h2,4,11,kSquare);
  h2->SetYTitle("Ratio to inclusive HDM");
  h2->GetXaxis()->SetTitleOffset(1.2);//1.5);
  //h2->SetMinimum(0.1); // for h2n
  //h2->SetMaximum(3e6); // for h2n
  //h2->SetMinimum(1e-4); // for hwi
  //h2->SetMaximum(10-1e-4); // for hwi
  //h2->SetMinimum(2e-3); // for hwf
  //h2->SetMaximum(5.); // for hwf
  //gPad->SetLogy();
  //h2->SetMinimum(0.9); // for hr
  //h2->SetMaximum(1.3); // for hr
  h2->SetMinimum(0.85); // for hri
  h2->SetMaximum(1.10); // for hri
  gPad->SetBottomMargin(0.15);//0.20);
  
  // Create smoothed variant of flavor response vs Z+jet
  TH2D *h2ris = (TH2D*)h2ref->Clone("h2ris"); //h2ris->Reset();
  //int color[nf] = {kBlack,kBlue,kOrange+2,kGreen+2,kRed};
  for (int iq = 0; iq != nf; ++iq) {
    int j = iq+1;
    //TH1D *hw = h2n->ProjectionX(Form("hw%d",j),j,j);
    //TH1D *hw = h2wi->ProjectionX(Form("hw%d",j),j,j);
    //TH1D *hw = h2wf->ProjectionX(Form("hw%d",j),j,j);
    //tdrDraw(hw,"HP",kNone,color[iq],kSolid,-1,kNone);
    //TH1D *hr = h2r->ProjectionX(Form("hr%d",j),j,j);
    TH1D *hr = h2ri->ProjectionX(Form("hr%d",j),j,j);
    //TH1D *hr = h2rf->ProjectionX(Form("hr%d",j),j,j);
    tdrDraw(hr,"P",kNone,color[iq],kSolid,-1,kNone);
    if (iq==0) hr->SetLineWidth(2);

    TH1D *hrd = h2rid->ProjectionX(Form("hrd%d",j),j,j);
    if (iq==0) tdrDraw(hrd,"P",kOpenCircle,color[iq],kSolid,-1,kNone);
    hrd->SetMarkerSize(0.5);

    //cout << hr->GetBinLowEdge(1) << endl << flush;
    //cout << hr->GetBinCenter(1) << endl << flush;
    //cout << hr->GetBinCenter(hr->GetNbinsX()+1) << endl << flush;
    //cout << hr->GetBinLowEdge(hr->GetNbinsX()+1) << endl << flush;

    TF1 *f1 = new TF1(Form("f1r_%d",iq),fR,xmin,xmax,
		      //hr->GetBinLowEdge(1)+1e-3,
		      //hr->GetBinLowEdge(hr->GetNbinsX()+1)-1e-3,
		      2*nt);
    f1->SetLineColor(color[iq]);
    // First free fit
    f1->SetParameters(1,0,1,0,1,0,1,0,1,0);
    hr->Fit(f1,"QRN");
    // Next, fix slope of each tag region to inclusive (high statistics)
    if (iq!=0) {
      f1->FixParameter(3,f1->GetParameter(1));
      f1->FixParameter(5,f1->GetParameter(1));
      f1->FixParameter(7,f1->GetParameter(1));
      f1->FixParameter(9,f1->GetParameter(1));
    }
    hr->Fit(f1,"QRN");
    // Finally, fix value to inclusive if not significant (>2 sigma)
    if (iq!=0) {
      double k = f1->GetChisquare()/f1->GetNDF();
      for (int it = 0; it != nt; ++it) {
	if (fabs(f1->GetParameter(2*it)-f1->GetParameter(0)) < 
	    2*k*f1->GetParError(2*it) && !(iq==1&&it==1))
	  f1->FixParameter(2*it, f1->GetParameter(0));
      }
    }
    hr->Fit(f1,"QRN");
    if (drawFit) 
      f1->Draw("SAME");

    // For reference, draw inclusive flavor responses to each tag region
    TF1 *fi = new TF1(Form("fi_%d",iq),fR,xmin,xmax,
		      //hr->GetBinLowEdge(1)+1e-3,
		      //hr->GetBinLowEdge(hr->GetNbinsX()+1)-1e-3,
		      2*nt);
    for (int k = 0; k != 2*nt; ++k) {
      fi->SetParameter(k, f1->GetParameter(k%2));
    }
    fi->SetLineColor(color[iq]);
    fi->SetLineStyle(kDotted);
    if (drawFit || true) fi->Draw("SAME");
    // => 1) gluon-tagged quarks jets have lower HDM response than average
    // => 2) b-tagged b-jets have higher HDM response than average
    //       (and q/c-tagged are lower, while gluons unbiased)

    // Copy fitted values for smoothed 2D
    // Edges will need to be recalculated with smoother counts for consistency
    for (int i = 1; i != h2ris->GetNbinsX()+1; ++i) {
      double txpt = h2ris->GetXaxis()->GetBinCenter(i);
      h2ris->SetBinContent(i,j,f1->Eval(txpt));
    }
  } // for i

  // Recalculate edges of smoothed response based on smoothed weights
  for (int ipt = 0; ipt != npt; ++ipt) {

    // First, "left" edge for flavor response
    int i0 = 0*npt + ipt + 1;
    double sumwf(0), sumwfr(0);
    for (int iq = 1; iq != nf; ++iq) {
      int j = iq+1;
      double sumwt(0), sumwtr(0);
      for (int it = 1; it != nt; ++it) {
	int i = it*npt + ipt + 1;
	sumwt += h2wis->GetBinContent(i,j);
	sumwf += h2wis->GetBinContent(i,j);
	sumwtr += h2wis->GetBinContent(i,j) * h2ris->GetBinContent(i,j);
	sumwfr += h2wis->GetBinContent(i,j) * h2ris->GetBinContent(i,j);
      } // for it
      // Exclusive flavor tags->inclusive flavor
      h2ris->SetBinContent(i0,j,sumwtr/sumwt);
      h2ris->SetBinContent(i0,j,sumwtr/sumwt);
    } // for iq
    // Exclusive flavor tags->fully inclusive Z+jet
    h2ris->SetBinContent(i0,1,sumwfr/sumwf);
    h2ris->SetBinContent(i0,1,sumwfr/sumwf);
    //assert(fabs(sumwf-1)<1e-3);
    if(!(fabs(sumwf-1)<1e-3))
      cout << "  sumwf-1=="<<sumwf-1<<endl;
    
    // Then, "bottom" edge for tagged jet response
    for (int it = 1; it != nt; ++it) {
      
      int i = it*npt + ipt + 1;
      double sumwq(0), sumwqr(0);
      for (int iq = 1; iq != nf; ++iq) {
	int j = iq+1;
	sumwq += h2wis->GetBinContent(i,j);
	sumwq += h2wis->GetBinContent(i,j);
	sumwqr += h2wis->GetBinContent(i,j) * h2ris->GetBinContent(i,j);
	sumwqr += h2wis->GetBinContent(i,j) * h2ris->GetBinContent(i,j);
      } // for it
      // Exclusive flavor tags->inclusive flavor
      h2ris->SetBinContent(i,1,sumwqr/sumwq);
      h2ris->SetBinContent(i,1,sumwqr/sumwq);
    } // for iq
  } // for ipt

  // Draw smoother plots on top
  for (int iq = 0; iq != nf; ++iq) {
    int j = iq+1;
    TH1D *hr = h2ris->ProjectionX(Form("hris%d",j),j,j);
    tdrDraw(hr,"HIST",kNone,color[iq],kSolid,-1,kNone);    
  }

  // Draw tag region boundaries
  l->SetLineStyle(kSolid);
  l->DrawLine(xmin+1*dx,0.85,xmin+1*dx,1.04);
  l->DrawLine(xmin+2*dx,0.85,xmin+2*dx,1.04);
  l->DrawLine(xmin+3*dx,0.85,xmin+3*dx,1.04);
  l->DrawLine(xmin+4*dx,0.85,xmin+4*dx,1.04);
  
  gPad->Update();
  c2->SaveAs("pdf/Zflavor_mcresp_v1.pdf");


  // Apply scaling to event fractions, then update responses
  TH2D *h2wiss = (TH2D*)h2wis->Clone("h2wiss");
  TH2D *h2wfss = (TH2D*)h2wfs->Clone("h2wfss");
  for (int ipt = 0; ipt != npt; ++ipt) {

    double pt = 0.5*(ptbins[ipt]+ptbins[ipt+1]);

    // Mapping to global bins in 2D
    int ii = 0*npt + ipt + 1;
    int iq = 1*npt + ipt + 1;
    int ig = 2*npt + ipt + 1;
    int ic = 3*npt + ipt + 1;
    int ib = 4*npt + ipt + 1;
    int ji = 0+1;
    int jq = 1+1;
    int jg = 2+1;
    int jc = 2+1;
    int jb = 4+1;

    // Correct for gluon efficiency
    bool fixGluon = true;
    if (fixGluon) {
      double ngi = h2wiss->GetBinContent(ii,jg);
      double ngq = h2wiss->GetBinContent(iq,jg);
      double ngg = h2wiss->GetBinContent(ig,jg);
      double effg = ngg/ngi;
      // Update this based on better calculation of gluon efficiency
      double ineffgsf;
      //double effgnew = effg*0.55; // 0.55 is a good guesstimate (or not)
      double effgsf = gluonEffSF(ptbins[ipt],ptbins[ipt+1],ineffgsf);
      double effgnew = effg*effgsf;
      double nggnew = ngi * effgnew; // = ngg*0.8
      double ngqnew = ngq + (ngg-nggnew);
      
      // Correct tagged event counts
      double niq = h2wiss->GetBinContent(iq,ji);
      double nig = h2wiss->GetBinContent(ig,ji);
      double niqnew = niq + (ngqnew-ngq);
      double nignew = nig + (nggnew-ngg);
      
      // Propagate changed counts to matrix
      h2wiss->SetBinContent(iq,jg, ngqnew);
      h2wiss->SetBinContent(ig,jg, nggnew);
      h2wiss->SetBinContent(iq,ji, niqnew);
      h2wiss->SetBinContent(ig,ji, nignew);
    } // fixGluon

    // Correct for quark efficiency
    bool fixQuark = true;
    if (fixQuark) {
      double nqi = h2wiss->GetBinContent(ii,jq);
      double nqq = h2wiss->GetBinContent(iq,jq);
      double nqg = h2wiss->GetBinContent(ig,jq);
      double effq = nqq/nqi;
      // Update this based on better calculation of gluon efficiency
      double ineffqsf;
      double effqsf = quarkEffSF(ptbins[ipt],ptbins[ipt+1],ineffqsf);
      double effqnew = effq*effqsf;
      double nqqnew = nqi * effqnew;
      double nqgnew = nqg + (nqq-nqqnew);
      
      // Correct tagged event counts
      double niq = h2wiss->GetBinContent(iq,ji);
      double nig = h2wiss->GetBinContent(ig,ji);
      double niqnew = niq + (nqqnew-nqq);
      double nignew = nig + (nqgnew-nqg);
      
      // Propagate changed counts to matrix
      h2wiss->SetBinContent(iq,jq, nqqnew);
      h2wiss->SetBinContent(ig,jq, nqgnew);
      h2wiss->SetBinContent(iq,ji, niqnew);
      h2wiss->SetBinContent(ig,ji, nignew);
    } // fixQuark
    
    // correct HF tagging efficiency
    bool fixHF = true;
    //bool fixB = true;
    //bool fixC = true;
    if (fixHF) {
      for (int jh = jc; jh != nf+1; ++jh)  {
      
	//int ih = jh*npt + ipt + 1;
	int ih = (jh==jb ? ib : ic);

	double nhi = h2wiss->GetBinContent(ii,jh);
	double nhq = h2wiss->GetBinContent(iq,jh);
	double nhg = h2wiss->GetBinContent(ig,jh);
	double nhc = h2wiss->GetBinContent(ic,jh);
	double nhb = h2wiss->GetBinContent(ib,jh);
	double nhh =  (jh==jb ? nhb : nhc);
	double effh = nhh/nhi;//(jh==jb ? nhb/nhi : nhc/nhi);
	double rq = nhq/(nhq+nhg);
	// Update this based on HF tagging efficiency. Guess 10% at low pT
	//double effhnew = effh * (pt < 100 ? 0.9 : 1.1);
	double x = h2wiss->GetXaxis()->GetBinCenter(ih);
	//double effsf = 1 + (fb->Eval(x)-1) / (jh==jb ? 0.93 : 0.60);
	double mistag;
	double effsf = heavyEffSF(pt,mistag);
	double effhnew = effh * effsf;
	double nhhnew = nhi * effhnew;
	double nhqnew = nhq + (nhh-nhhnew) * rq;
	double nhgnew = nhg + (nhh-nhhnew) * (1-rq);
	
	// Correct tagged event counts
	double niq = h2wiss->GetBinContent(iq,ji);
	double nig = h2wiss->GetBinContent(ig,ji);
	double nih = h2wiss->GetBinContent(ih,ji);
	double niqnew = niq + (nhqnew-nhq);
	double nignew = nig + (nhgnew-nhg);
	double nihnew = nih + (nhhnew-nhh);
	
	// Propagate changed counts to matrix
	h2wiss->SetBinContent(iq,jh, nhqnew);
	h2wiss->SetBinContent(ig,jh, nhgnew);
	h2wiss->SetBinContent(ih,jh, nhhnew);
	//
	h2wiss->SetBinContent(iq,ji, niqnew);
	h2wiss->SetBinContent(ig,ji, nignew);
	h2wiss->SetBinContent(ih,ji, nihnew);
      } // iq
    } // fixHF

    // Propagate changed fractions to matrix
    for (int it = 0; it != nt; ++it) {
      for (int iq = 0; iq != nf; ++iq) {
	int i = it*npt + ipt + 1;
	int j = iq+1;
	h2wfss->SetBinContent(i,j,h2wiss->GetBinContent(i,j) / 
			      h2wiss->GetBinContent(i,1));
      } // for iq
    } // for it
  } // for ipt

  // Draw updated counts
  c0->cd();

  TH1D *hwmv2 = h2wiss->ProjectionX("hwimv2",1,1);
  TH1D *hwdv2 = h2wid->ProjectionX("hwidv2",1,1);
  hwdv2->Divide(hwmv2);  

  tdrDraw(hwdv2,"HIST",kNone,color[0],kDotted,-1,kNone);    

  // Fit b rate to estimate tagging efficiency vs pT
  TF1 *fbv2 = new TF1("fbv2",fW,xmin,xmax,_npw*nt);
  fbv2->SetLineColor(kBlue);
  hwdv2->Fit(fbv2,"QRN");
  fbv2->Draw("SAME");

  c0->SaveAs("pdf/Zflavor_counts_v2.pdf");

  // Draw updated fractions
  c1->cd();

  // Draw smoother plots on top
  for (int iq = 0; iq != nf; ++iq) {
    int j = iq+1;
    TH1D *hwi = h2wiss->ProjectionX(Form("hwiss%d",j),j,j);
    TH1D *hwf = h2wfss->ProjectionX(Form("hwfss%d",j),j,j);
    if (iq==0) tdrDraw(hwi,"HIST",kOpenSquare,color[iq],kSolid,-1,kNone);    
    //if (iq==0) tdrDraw(hwi,"HIST",kOpenSquare,color[iq],kDotted,-1,kNone);    
    //tdrDraw(hwf,"HIST",kNone,color[iq],kDotted,-1,kNone);    
    tdrDraw(hwf,"HIST",kNone,color[iq],kSolid,-1,kNone);    

    if (iq!=0 && iq<2) {
      //leg1d->AddEntry(hw,Form("%s / tag",flvlabel[flvbins[iq]]),"L");
      leg1d2->AddEntry(hwf," ","L");
    }
    if (iq!=0 && iq>=2) {
      //leg1f->AddEntry(hw,Form("%s",flvlabel[flvbins[iq]]),"L");
      leg1f2->AddEntry(hwf," ","L");
    }
  }
  c1->SaveAs("pdf/Zflavor_mcfrac_v2.pdf");


  // Solve inclusive flavor fractions based on tagged counts
  // Ntag = [Epsilon] * Nflavor => Nflavor = [Epsilon]^{-1} * Ntag
  // ft = Ntag/Npt, ff = Nflavor/Npt
  // fd is ft for data

  // First task is to expand data into square matrix that can be inverted
  TMatrixD MEm((nt-1)*npt,(nt-1)*npt), MEd((nt-1)*npt,(nt-1)*npt);
  TMatrixD ftm((nt-1)*npt,1), ffm((nt-1)*npt,1), ff((nt-1)*npt,1);
  TMatrixD ftd((nt-1)*npt,1), ffd((nt-1)*npt,1);
  TH2D *h2mem = h2wis; // use: h2wi, h2wis (default)
  TH2D *h2med = h2wiss; // use: h2wi, h2wis, h2wiss (default)

  for (int ipt = 0; ipt != npt; ++ipt) {
    for (int it = 1; it != nt; ++it) { // skip inclusive '0'
      for (int iq = 0; iq != nf; ++iq) { // inclusive '0' separate treatment
	
	// Indeces in matrix, excludes inclusive flavor and tag
	int kt = npt*(it-1) + ipt;
	int kf = npt*(iq-1) + ipt;
	// Indeces in TH2D, offset by 1 and no duplicated pT dimension in y
	int i = npt*it + ipt + 1;
	int j = iq + 1;
	int k = ipt + 1;
	// value in x-axis
	double x = h2ref->GetXaxis()->GetBinCenter(i);

	// Tag fractions ft, smoothed for both data and MC
	if (j==1) {
	  ftm[kt][0] =
	    h2mem->GetBinContent(i,1) /
	    h2mem->GetBinContent(k,1);
	  ftd[kt][0] =
	    h2med->GetBinContent(i,1) * fbv2->Eval(x) /
	    h2med->GetBinContent(k,1);
	}
	// Use smoothed MC for ff and ME matrix
	if (j!=1) {
	  ff[kf][0] =
	    h2mem->GetBinContent(k,j) /
	    h2mem->GetBinContent(k,1);
	  MEm[kt][kf] =
	    h2mem->GetBinContent(i,  j) /
	    h2mem->GetBinContent(k,j);
	  MEd[kt][kf] =
	    h2med->GetBinContent(i,  j) /
	    h2med->GetBinContent(k,j);
	}
      } // for iq
    } // for it
  } // for ipt

  // Sanity check of corrected index order, ftm == MEm * ff
  if (true) {
    TMatrixD fttm((nt-1)*_npt,1);
    fttm.Mult(MEm,ff);
    assert(fttm.GetNcols()==1);
    for (int i = 0; i != fttm.GetNrows(); ++i) {
      if (fabs(fttm[i][0]/ftm[i][0]-1)>1e-4)
	cout << "i="<<i<<", ftm="<<ftm[i][0]<<", ratio="
	     << fttm[i][0]/ftm[i][0] << endl;
    }
  } // sanity check

  // Invert efficiency matrix to solve true fractions
  // Check that inversion succeeded (det!=0)
  Double_t detMEm, detMEd;
  TMatrixD MEmi = MEm;
  TMatrixD MEdi = MEd;
  MEmi.Invert(&detMEm);
  MEdi.Invert(&detMEd);
  std::cout << "  Determinant(MEm)     = " << detMEm << std::endl;
  std::cout << "  Determinant(MEd)     = " << detMEd << std::endl;
  //assert(det1m!=0);

  // Solve flavor fractions in (smoothed) data
  //TMatrixD ffm((nt-1)*npt,1);
  ffm.Mult(MEmi,ftm);
  //TMatrixD ffd((nt-1)*npt,1);
  ffd.Mult(MEdi,ftd);
  TMatrixD fttd((nt-1)*npt,1); // remultiplication check
  fttd.Mult(MEd,ffd);

  // Sanity check matrix inversion for MC where we know the answer,
  // and for data, where we can remultiply solved fraction back to tagged
  if (true) {
    assert(ffm.GetNcols()==1);
    for (int i = 0; i != ffm.GetNrows(); ++i) {
      if (fabs(ffm[i][0]/ff[i][0]-1)>1e-4)
	cout << "i="<<i<<", ffm="<<ffm[i][0]<<", ratio="
	     << ffm[i][0]/ff[i][0] << endl;
      if (fabs(fttd[i][0]/ftd[i][0]-1)>1e-4)
	cout << "i="<<i<<", ftd="<<ftd[i][0]<<", ratio="
	     << fttd[i][0]/ftd[i][0] << endl;
    }
  } // sanity checks

  // Create matrices mapping responses to tag reponse
  TMatrixD MEWRm((nt-1)*npt,(nt-1)*npt);
  TMatrixD MEWRd((nt-1)*npt,(nt-1)*npt);
  TMatrixD rtm((nt-1)*npt,1), rfm((nt-1)*npt,1), rf((nt-1)*npt,1);
  TMatrixD rtd((nt-1)*npt,1), rfd((nt-1)*npt,1);
  
  // Create TH2D based on new flavor fractions solved from data
  TH2D *h2wissd = (TH2D*)h2ref->Clone("h2wissd"); h2wissd->Reset();
  TH2D *h2wfssd = (TH2D*)h2ref->Clone("h2wfssd"); h2wfssd->Reset();

  for (int ipt = 0; ipt != npt; ++ipt) {
    for (int it = 1; it != nt; ++it) { // skip inclusive '0'
      for (int iq = 1; iq != nf; ++iq) { // inclusive '0' separate treatment
	
	// Indeces in matrix, excludes inclusive flavor and tag
	int kt = npt*(it-1) + ipt;
	int kf = npt*(iq-1) + ipt;
	// Indeces in TH2D, offset by 1 and no duplicated pT dimension in y
	int i = npt*it + ipt + 1;
	int j = iq + 1;
	int k = ipt + 1;
	// value in x-axis
	double x = h2ref->GetXaxis()->GetBinCenter(i);

	double Ntf = MEd[kt][kf]*ffd[kf][0];
	h2wissd->SetBinContent(i,j,Ntf); 
	h2wissd->SetBinContent(k,j,h2wissd->GetBinContent(k,j)+Ntf);
	h2wissd->SetBinContent(i,1,h2wissd->GetBinContent(i,1)+Ntf);
	h2wissd->SetBinContent(k,1,h2wissd->GetBinContent(k,1)+Ntf);

	// Tag fractions ft, smoothed for both data and MC
	if (j==1) {
	  rtm[kt][0] = 
	    h2ris->GetBinContent(i,1) /
	    h2ris->GetBinContent(k,1);
	    //h2ris->GetBinContent(i,1) /
	    //h2ris->GetBinContent(k,1);
	  rtd[kt][0] =
	    h2rid->GetBinContent(i,1) /
	    h2rid->GetBinContent(k,1);
	    //h2ris->GetBinContent(i,1) * f1rd->Eval(x) /
	    //h2ris->GetBinContent(k,1);
	}
	// Use smoothed MC for ff and ME matrix
	if (j!=1) {
	  rf[kf][0] =
	    h2ris->GetBinContent(k,j) /
	    h2ris->GetBinContent(k,1);
	  MEWRm[kt][kf] = MEm[kt][kf] * ffm[kf][0] / ftm[kt][0] *
	    h2ris->GetBinContent(i,j) /
	    h2ris->GetBinContent(k,j);
	  MEWRd[kt][kf] = MEd[kt][kf] * ffd[kf][0] / ftd[kt][0] *
	    h2ris->GetBinContent(i,  j) /
	    h2ris->GetBinContent(k,j);
	}

      } // for iq
    } // for it

    // Propagate changed fractions to matrix
    for (int it = 0; it != nt; ++it) {
      for (int iq = 0; iq != nf; ++iq) {
	int i = it*npt + ipt + 1;
	int j = iq+1;
	h2wfssd->SetBinContent(i,j,h2wissd->GetBinContent(i,j) / 
			       h2wissd->GetBinContent(i,1));
      } // for iq
    } // for it
  } // for ipt


  // Update data/MC with data-customised smoothed and scaled matrix
  c0->cd();

  TH1D *hwmv3 = h2wissd->ProjectionX("hwimv3",1,1);
  TH1D *hwdv3 = h2wid->ProjectionX("hwidv3",1,1);
  hwdv3->Divide(hwmv3);  

  tdrDraw(hwdv3,"HIST",kViolet+2,color[0],kDotted,-1,kNone);    
  hwdv3->SetLineWidth(2);

  // Fit b rate to estimate tagging efficiency vs pT
  TF1 *fbv3 = new TF1("fbv3",fW,xmin,xmax,_npw*nt);
  fbv3->SetLineColor(kViolet+2);
  hwdv3->Fit(fbv3,"QRN");
  fbv3->Draw("SAME");

  c0->SaveAs("pdf/Zflavor_counts_v3.pdf");


  // Final v3 update for mcfrac with MC normalized to data in two days
  TH1D *h3 = (TH1D*)h1->Clone("h3");
  TCanvas *c3 = tdrCanvas("c3",h3,4,11,kSquare);
  h3->GetXaxis()->SetTitleOffset(1.20);//1.5);
  gPad->SetBottomMargin(0.15);//0.20);

  //TLegend *leg3d = tdrLeg(0.45,0.90-3*0.05,0.75,0.90);
  //leg3d->SetHeader("Data");
  //TLegend *leg3m = tdrLeg(0.37,0.90-3*0.05,0.67,0.90);
  //leg3m->SetHeader("MC");
  TLegend *leg3d = tdrLeg(0.45,0.90-3*0.05,0.75,0.90);
  TLegend *leg3f = tdrLeg(0.76,0.90-(nt+1-3)*0.05,1.06,0.90);
  TLegend *leg3d2 = tdrLeg(0.45,0.90-3*0.05+0.01,0.75,0.90+0.01);
  TLegend *leg3f2 = tdrLeg(0.76,0.90-(nt+1-3)*0.05+0.01,1.05,0.90+0.01);
  
  if (true) { // now in separate Zflavor_counts.pdf
    TH1D *hwim0 = h2wis->ProjectionX(Form("hwim0s%d",1),1,1);
    TH1D *hwim1 = h2wiss->ProjectionX(Form("hwim1ss%d",1),1,1);
    TH1D *hwim2 = h2wissd->ProjectionX(Form("hwim2ssd%d",1),1,1);
    TH1D *hwim = h2wi->ProjectionX(Form("hwimm%d",1),1,1);
    TH1D *hwid = h2wid->ProjectionX(Form("hwidss%d",1),1,1);
    tdrDraw(hwim0,"HIST",kNone,color[0],kDotted,-1,kNone);
    //tdrDraw(hwim1,"HIST",kNone,color[0],kDashed,-1,kNone);
    tdrDraw(hwim2,"HIST",kNone,color[0],kSolid,-1,kNone);
    tdrDraw(hwim,"P",kOpenCircle,color[0],kDotted,-1,kNone);
    tdrDraw(hwid,"P",kFullCircle,color[0],kSolid,-1,kNone);
    hwid->SetMarkerSize(0.5);
    hwim->SetMarkerSize(0.5);

    //leg3d->AddEntry(hwid,"Tag / incl.","PL");
    //leg3m->AddEntry(hwim," ","PL");
    leg3d->AddEntry(hwid,"Data / incl.","PL");
    leg3d->AddEntry(hwim,"MC / incl.","PL");
    leg3d2->AddEntry(hwid," ","");
    leg3d2->AddEntry(hwim," ","");
  }

  for (int iq = 0; iq != nf; ++iq) {
    int j = iq+1;
    TH1D *hw0 = h2wfs->ProjectionX(Form("hw0s%d",j),j,j);
    TH1D *hw1 = h2wfss->ProjectionX(Form("hw1ss%d",j),j,j);
    TH1D *hw2 = h2wfssd->ProjectionX(Form("hw2ssd%d",j),j,j);
    tdrDraw(hw0,"HIST",kNone,color[iq],kDotted,-1,kNone);
    //tdrDraw(hw1,"HIST",kNone,color[iq],kDashed,-1,kNone);
    tdrDraw(hw2,"HIST",kNone,color[iq],kSolid,-1,kNone);

    if (iq!=0 && iq<2) {
      leg3d->AddEntry(hw2,Form("%s / tag",flvlabel[flvbins[iq]]),"L");
      leg3d2->AddEntry(hw0," ","L");
      //leg3m->AddEntry(hw0," ","L");
    }
    if (iq!=0 && iq>=2) {
      leg3f->AddEntry(hw2,Form("%s",flvlabel[flvbins[iq]]),"L");
      leg3f2->AddEntry(hw0," ","L");
      //leg3f->AddEntry(hw0," ","L");
    }
  }

  // Draw tag region boundaries
  l->SetLineStyle(kSolid);
  l->DrawLine(xmin+1*dx,0.,xmin+1*dx,1.);
  l->DrawLine(xmin+2*dx,0.,xmin+2*dx,1.);
  l->DrawLine(xmin+3*dx,0.,xmin+3*dx,1.);
  l->DrawLine(xmin+4*dx,0.,xmin+4*dx,1.);

  gPad->Update();
  c3->SaveAs("pdf/Zflavor_mcfrac_v3.pdf");


  // Combine smoothed response and fraction matrices into single matrix
  // To be used for unfolding. Update with more smoothed versions
  TH2D *h2rnis = (TH2D*)h2ref->Clone("h2rnis");
  for (int i = 1; i != h2ref->GetNbinsX()+1; ++i) {
    for (int j = 1; j != h2ref->GetNbinsY()+1; ++j) {
      h2rnis->SetBinContent(i, j, h2ris->GetBinContent(i,j) *
			    h2wfs->GetBinContent(i,j));
     } // for j
  } // for i

  // c3 originally here

  // Recalculate edges of smoothed response based on smoothed and scaled weights
  TH2D *h2riss = (TH2D*)h2ris->Clone("h2riss");
  // Combine with (quark and) gluon scale shift
  TH2D *h2rix = (TH2D*)h2ris->Clone("h2rix");
  // ..and with flavor fractions from data
  TH2D *h2rixd = (TH2D*)h2ris->Clone("h2rixd");
  for (int ipt = 0; ipt != npt; ++ipt) {

    // Before starting, scale quark and gluon responses
    for (int it = 0; it != nt; ++it) {
      int i = it*npt + ipt + 1;
      int jq = 1+1;
      int jg = 2+1;
      //h2rix->SetBinContent(i,jq, h2ris->GetBinContent(i,jq)*1.010);
      //h2rix->SetBinContent(i,jg, h2ris->GetBinContent(i,jg)*0.970);
      // Inclusive shift of quark and gluon response

      h2rix->SetBinContent(i,jq, h2ris->GetBinContent(i,jq)*quarkScale);
      h2rix->SetBinContent(i,jg, h2ris->GetBinContent(i,jg)*gluonScale);
      h2rixd->SetBinContent(i,jq, h2ris->GetBinContent(i,jq)*quarkScale);
      h2rixd->SetBinContent(i,jg, h2ris->GetBinContent(i,jg)*gluonScale);
      // Extra shifts for Q/G tagged quarks and gluons
      if (it==1) { // quark tag
	h2rix->SetBinContent(i,jq, h2rix->GetBinContent(i,jq)*quarkScaleQtag);
	h2rix->SetBinContent(i,jg, h2rix->GetBinContent(i,jg)*gluonScaleQtag);
	h2rixd->SetBinContent(i,jq, h2rix->GetBinContent(i,jq)*quarkScaleQtag);
	h2rixd->SetBinContent(i,jg, h2rix->GetBinContent(i,jg)*gluonScaleQtag);
      }
      if (it==2) { // gluon tag
	h2rix->SetBinContent(i,jq, h2rix->GetBinContent(i,jq)*quarkScaleGtag);
	h2rix->SetBinContent(i,jg, h2rix->GetBinContent(i,jg)*gluonScaleGtag);
	h2rixd->SetBinContent(i,jq, h2rix->GetBinContent(i,jq)*quarkScaleGtag);
	h2rixd->SetBinContent(i,jg, h2rix->GetBinContent(i,jg)*gluonScaleGtag);
      }

    }

    // First, "left" edge for flavor response
    int i0 = 0*npt + ipt + 1;
    double sumwf(0), sumwfr(0), sumwfx(0), sumwfd(0), sumwfxd(0);
    for (int iq = 1; iq != nf; ++iq) {
      int j = iq+1;
      double sumwt(0), sumwtr(0), sumwtx(0), sumwtd(0),sumwtxd(0);
      for (int it = 1; it != nt; ++it) {
	int i = it*npt + ipt + 1;
	sumwt += h2wiss->GetBinContent(i,j);
	sumwtd += h2wissd->GetBinContent(i,j);
	sumwf += h2wiss->GetBinContent(i,j);
	sumwfd += h2wissd->GetBinContent(i,j);
	sumwtr += h2wiss->GetBinContent(i,j) * h2riss->GetBinContent(i,j);
	sumwfr += h2wiss->GetBinContent(i,j) * h2riss->GetBinContent(i,j);
	sumwtx += h2wiss->GetBinContent(i,j) * h2rix->GetBinContent(i,j);
	sumwfx += h2wiss->GetBinContent(i,j) * h2rix->GetBinContent(i,j);
	sumwtxd += h2wissd->GetBinContent(i,j) * h2rixd->GetBinContent(i,j);
	sumwfxd += h2wissd->GetBinContent(i,j) * h2rixd->GetBinContent(i,j);
      } // for it
      // Exclusive flavor tags->inclusive flavor
      h2riss->SetBinContent(i0,j,sumwtr/sumwt);
      h2riss->SetBinContent(i0,j,sumwtr/sumwt);
      h2rix->SetBinContent(i0,j,sumwtx/sumwt);
      h2rix->SetBinContent(i0,j,sumwtx/sumwt);
      h2rixd->SetBinContent(i0,j,sumwtxd/sumwtd);
      h2rixd->SetBinContent(i0,j,sumwtxd/sumwtd);
    } // for iq

    // Bottom left corner: Exclusive flavor tags->fully inclusive Z+jet
    h2riss->SetBinContent(i0,1,sumwfr/sumwf);
    h2riss->SetBinContent(i0,1,sumwfr/sumwf);
    h2rix->SetBinContent(i0,1,sumwfx/sumwf);
    h2rix->SetBinContent(i0,1,sumwfx/sumwf);
    h2rixd->SetBinContent(i0,1,sumwfxd/sumwfd);
    h2rixd->SetBinContent(i0,1,sumwfxd/sumwfd);
    //assert(fabs(sumwf-1)<1e-3);
    if(!(fabs(sumwf-1)<1e-3))
      cout << "  sumwf-1=="<<sumwf-1<<endl;
    
    // Then, "bottom" edge for tagged jet response
    for (int it = 1; it != nt; ++it) {
      
      int i = it*npt + ipt + 1;
      double sumwq(0), sumwqr(0), sumwqx(0), sumwqd(0), sumwqxd(0);
      for (int iq = 1; iq != nf; ++iq) {
	int j = iq+1;
	sumwq += h2wiss->GetBinContent(i,j);
	sumwqd += h2wissd->GetBinContent(i,j);
	sumwqr += h2wiss->GetBinContent(i,j) * h2riss->GetBinContent(i,j);
	sumwqx += h2wiss->GetBinContent(i,j) * h2rix->GetBinContent(i,j);
	sumwqxd += h2wissd->GetBinContent(i,j) * h2rixd->GetBinContent(i,j);
      } // for it
      // Exclusive flavor tags->inclusive flavor
      h2riss->SetBinContent(i,1,sumwqr/sumwq);
      h2rix->SetBinContent(i,1,sumwqx/sumwq);
      h2rixd->SetBinContent(i,1,sumwqxd/sumwqd);
    } // for iq
  } // for ipt

  c2->cd();

  // Draw smoother plots on top
  for (int iq = 0; iq != nf; ++iq) {
    int j = iq+1;
    TH1D *hr = h2riss->ProjectionX(Form("hriss%d",j),j,j);
    TH1D *hx = h2rix->ProjectionX(Form("hrix%d",j),j,j);
    TH1D *hxd = h2rixd->ProjectionX(Form("hrixd%d",j),j,j);
    //tdrDraw(hr,"HIST",kNone,color[iq],kDashed,-1,kNone);
    tdrDraw(hx,"HIST",kNone,color[iq],kDotted,-1,kNone);
    tdrDraw(hxd,"HIST",kNone,color[iq],kDashed,-1,kNone);
  }

  c2->SaveAs("pdf/Zflavor_mcresp_v2.pdf");



  TH1D *h4 = (TH1D*)href->Clone("h4");
  h4->SetYTitle("Ratio to inclusive HDM");
  
  // Ratio to inclusive HDM for data, MC, fits
  TCanvas *c4 = tdrCanvas("c4",h4,4,11,kSquare);
  h4->GetXaxis()->SetTitleOffset(1.2);//1.5);
  h4->SetMinimum(0.93);
  h4->SetMaximum(1.03); 
  gPad->SetBottomMargin(0.15);//0.20);

  TH1D *hrim = h2ri->ProjectionX(Form("hrim%d",10),1,1);
  TH1D *hrid = h2rid->ProjectionX(Form("hrid%d",10),1,1);
  TH1D *hrims = h2ris->ProjectionX(Form("hrims%d",10),1,1);
  TH1D *hrimss = h2riss->ProjectionX(Form("hrimss%d",10),1,1);
  TH1D *hrimx = h2rix->ProjectionX(Form("hrimx%d",10),1,1);
  TH1D *hrimxd = h2rixd->ProjectionX(Form("hrimxd%d",10),1,1);

  l->SetLineStyle(kDotted);
  l->DrawLine(xmin,1.1,xmax,1.1);
  l->DrawLine(xmin,0.9,xmax,0.9);
  l->SetLineStyle(kDashed);
  l->DrawLine(xmin,1,xmax,1);
  
  tdrDraw(hrimxd,"HIST",kNone,color[0],kSolid,-1,kNone);
  //tdrDraw(hrimx,"HIST",kNone,color[0],kDashed,-1,kNone);
  //tdrDraw(hrimx,"HIST",kNone,color[0],kSolid,-1,kNone);
  tdrDraw(hrimss,"HIST",kNone,color[0],kDashed,-1,kNone);
  tdrDraw(hrims,"HIST",kNone,color[0],kDotted,-1,kNone);
  tdrDraw(hrim,"P",kOpenCircle,color[0],kSolid,-1,kNone);
  tdrDraw(hrid,"P",kFullCircle,color[0],kSolid,-1,kNone);
  hrims->SetLineWidth(2);
  hrim->SetMarkerSize(0.6);
  hrid->SetMarkerSize(0.5);
  hrid->SetLineWidth(2);
  
  TLegend *leg4 = tdrLeg(0.50,0.90-0.040*4,0.80,0.90);
  leg4->SetTextSize(0.040);
  leg4->AddEntry(hrid,"Data","PLE");
  leg4->AddEntry(hrim,"MC","PLE");
  leg4->AddEntry(hrims,"smoothed","L");
  //leg4->AddEntry(hrimss,"& scaled fractions","L");
  leg4->AddEntry(hrimss,"+ scaled efficiency","L");
  //leg4->AddEntry(hrimss,"& scaled efficiency","L");
  //leg4->AddEntry(hrimx,"& response","L");
  //leg4->AddEntry(hrimx,"+ scaled eff.&R","L");
  //leg4->AddEntry(hrimxd,"+ data fractions","L");
  if (quarkScale==1 && gluonScale==1)
      leg4->AddEntry(hrimxd,"+ data fractions","L");
  else
    leg4->AddEntry(hrimxd,"+ data fractions & R","L");


  // Draw tag region boundaries
  l->SetLineStyle(kSolid);
  l->DrawLine(xmin+1*dx,0.93,xmin+1*dx,1.012);
  l->DrawLine(xmin+2*dx,0.93,xmin+2*dx,1.012);
  l->DrawLine(xmin+3*dx,0.93,xmin+3*dx,1.00);
  l->DrawLine(xmin+4*dx,0.93,xmin+4*dx,1.00);
 

  gPad->Update();
  c4->SaveAs("pdf/Zflavor_resp.pdf");


  // Plot efficiency matrix
  //TH1D *h5 = (TH1D*)href->Clone("h5");
  //TCanvas *c5 = tdrCanvas("c5",h5,4,11,kSquare);
  

} // Zflavor


// gluonEffSF copied from drawZflavor.C
// hadW_Zb.C also has in-line version of this
//TF1 *fqgl(0);
//TFile *fi(0);
TF1 *f1mg(0), *f1dg(0);
double gluonEffSF(double ptmin, double ptmax, double &ineffsf) {

  // hadW_QGL.C QGL>0.5 efficiency fits
  // For Z+jet, multiply MC results by 1.05
  TF1 *f1mg = new TF1("f1mg","[0]+[1]*pow(x,[2])",34,200);
  TF1 *f1dg = new TF1("f1dg","[0]+[1]*pow(x,[2])",34,200);
  // Fit with Z+jet included
  //f1mg->SetParameters(0.2468,0.4722,-0.5116);
  //f1dg->SetParameters(0.4225,0.7519,-0.4858);
  // Fit without Z+jet (only dijet and TT+jet)
  f1mg->SetParameters(0.2534,0.3997,-0.5167);
  f1dg->SetParameters(0.4329,0.65,-0.4916);

  double pt = 0.5*(ptmin+ptmax);
  //double effm = 1.05*f1mg->Eval(pt);
  double effm = 1.02*f1mg->Eval(pt);
  double effd = f1dg->Eval(pt);
  
  // Efficiencies above are for QGL>0.5, so 1-eff(g)
  double effsf = (1-effd) / (1-effm);
  ineffsf = effd / effm;

  return effsf;

  /*
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/QuarkGluonLikelihood#Systematics
//  Pythia8 gluon (|eta| < 2.0, inclusive pT)
//  = 2.5626*x^3 - 3.2240*x^2 + 1.8687*x + 0.6770,
//  where x is the value of the qgl output (i.e. central value of the QGL bin).

  if (!fqgl) fqgl = new TF1(Form("fqll_%1.0f",ptmin),
			    "2.5626*x^3 - 3.2240*x^2 + 1.8687*x + 0.6770",0,1);

  // Use rootfiles/output-MC-1-UL17V4_BCDEF.root/Standard/Eta_0.0-1.3/mc
  // and TH2D hqgl2_g_g for estimating efficiency scale factor
  if (!fi) fi = new TFile("rootfiles/output-MC-1-UL17V4_BCDEF.root","READ");
  assert(fi && !fi->IsZombie());
  fi->cd("Standard/Eta_0.0-1.3/mc");
  TH2D *h2qgl = (TH2D*)gDirectory->Get("hqgl2_g_g"); assert(h2qgl);
  int i1 = h2qgl->GetXaxis()->FindBin(ptmin);
  int i2 = h2qgl->GetXaxis()->FindBin(ptmax-1);
  TH1D *hqgl = h2qgl->ProjectionY(Form("hqgl_%1.0f",ptmin),i1,i2);
  TH1D *hqglw = (TH1D*)hqgl->Clone(Form("hqglw_%1.0f",ptmin));
  hqglw->Multiply(fqgl);
  int i3 = hqgl->FindBin(0.5-1e-5);
  double effm = hqgl->Integral(1,i3) / hqgl->Integral();
  double effd = hqglw->Integral(1,i3)/hqglw->Integral();

  double effsf = (effd / effm);
  ineffsf = (1-effd) / (1-effm);

  return effsf;
  */
}
// Gluoneffsf

// Custom quark SF, otherwise same as gluon
//TF1 *fqgl_q(0);
TF1 *f1mq(0), *f1dq(0);
double quarkEffSF(double ptmin, double ptmax, double &ineffsf) {

  // hadW_QGL.C QGL>0.5 efficiency fits
  // For Z+jet, multiply MC results by 1.05
  if (!f1mq) f1mq = new TF1("f1mq","[0]+[1]*pow(x,[2])",34,200);
  if (!f1dq) f1dq = new TF1("f1dq","[0]+[1]*pow(x,[2])",34,200);
  // Fit with Z+jet included
  //f1mq->SetParameters(0.7315,-15.79,-1.346);
  //f1dq->SetParameters(0.7926,-17.77,-1.357);
  // Fit without Z+jet (W>qq' and dijet only)
  f1mq->SetParameters(0.7325,-12.78,-1.287);
  f1dq->SetParameters(0.7936,-14.48,-1.299);

  double pt = 0.5*(ptmin+ptmax);
  //double effm = 1.05*f1mq->Eval(pt);
  double effm = 1.02*f1mq->Eval(pt);
  double effd = f1dq->Eval(pt);
  
  double effsf = (effd / effm);
  ineffsf = (1-effd) / (1-effm);

  return effsf;

  /*
  // Extracted using rootfiles/hadW[UL,MC]1718V5_EMUF.root/h2qgl
  if (!fqgl_q) fqgl_q = new TF1(Form("fql_%1.0f",ptmin),"(x<0.99)*(2.204*x^3 - 3.527*x^2 + 2.021*x + 0.5724)+(x>=0.99)*1.1395",0,1);

  // Use rootfiles/output-MC-1-UL17V4_BCDEF.root/Standard/Eta_0.0-1.3/mc
  // and TH2D hqgl2_g_g for estimating efficiency scale factor
  if (!fi) fi = new TFile("rootfiles/output-MC-1-UL17V4_BCDEF.root","READ");
  assert(fi && !fi->IsZombie());
  fi->cd("Standard/Eta_0.0-1.3/mc");
  TH2D *h2qgl = (TH2D*)gDirectory->Get("hqgl2_q_g"); assert(h2qgl);
  int i1 = h2qgl->GetXaxis()->FindBin(ptmin);
  int i2 = h2qgl->GetXaxis()->FindBin(ptmax-1);
  TH1D *hqgl = h2qgl->ProjectionY(Form("hqglq_%1.0f",ptmin),i1,i2);
  TH1D *hqglw = (TH1D*)hqgl->Clone(Form("hqglqw_%1.0f",ptmin));
  hqglw->Multiply(fqgl_q);
  //int i3 = hqgl->FindBin(0.5-1e-5);
  //double effm = hqgl->Integral(1,i3) / hqgl->Integral();
  //double effd = hqglw->Integral(1,i3)/hqglw->Integral();
  int i3 = hqgl->FindBin(0.5+1e-5);
  double effm = hqgl->Integral(i3,hqgl->GetNbinsX()) / hqgl->Integral();
  double effd = hqglw->Integral(i3,hqgl->GetNbinsX())/hqglw->Integral();

  double effsf = (effd / effm);
  ineffsf = (1-effd) / (1-effm);

  return effsf;
  */
} // quarkEffSF

TF1 *fbmu18(0), *fbcm18(0), *fbmu17(0), *fbcm17(0);
TF1 *finc18(0), *finc17(0);
double heavyEffSF(double pt, double &mistag) {
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/BTV13TeV2017DeepJet
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation106XUL17
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation106XUL18

  // 2, mujets, central, 1, -2.5, 2.5, 20, 1000, 0, 1, "0.885007+0.00159965*log(x+19)*log(x+18)*(3-0.320398*log(x+18))" 
  //2, comb, central, 1, -2.5, 2.5, 20, 1000, 0, 1, "0.716648+0.00833545*log(x+19)*log(x+18)*(3-0.370069*log(x+18))" 

  // Tight working point, UL18
  if (!fbmu18) fbmu18 = new TF1("fbmu18","0.885007+0.00159965*log(x+19)*log(x+18)*(3-0.320398*log(x+18))",20,1000);
  if (!fbcm18) fbcm18 = new TF1("fbcm18","0.716648+0.00833545*log(x+19)*log(x+18)*(3-0.370069*log(x+18))",20,1000);

  // Medium and loose WP, UL18
  //TF1 *fbcm2 = new TF1("fbcm2","0.763354+0.0081767*log(x+19)*log(x+18)*(3-0.399925*log(x+18))",20,1000);
  //TF1 *fbcm3 = new TF1("fbcm2","0.882297+0.00426142*log(x+19)*log(x+18)*(3-0.411195*log(x+18))",20,1000);


  // 2, mujets, central, 1, -2.5, 2.5, 20, 1000, 0, 1, "1.28462+(-(0.0149197*(log(x+19)*(log(x+18)*(3-(0.392618*log(x+18)))))))" 
  //2, comb, central, 1, -2.5, 2.5, 20, 1000, 0, 1, "0.88695+(0.000743492*(log(x+19)*(log(x+18)*(3-(-(0.0812416*log(x+18)))))))" 

  // Tight working point, UL17
  if (!fbmu17) fbmu17 = new TF1("fbmu17","1.28462+(-(0.0149197*(log(x+19)*(log(x+18)*(3-(0.392618*log(x+18)))))))",20,1000); // weird shape?
  if (!fbcm17) fbcm17 = new TF1("fbcm17","0.88695+(0.000743492*(log(x+19)*(log(x+18)*(3-(-(0.0812416*log(x+18)))))))",20,1000); 

  
  // Mistag rates UL18 and UL17
  //2, incl, central, 2, 0, 2.4, 20, 1000, 0, 1, "0.864506+2.79354/sqrt(x)" 
  //2, incl, central, 2, 0, 2.4, 20, 1000, 0, 1, "0.850069+1.99726/sqrt(x)" 
  
  if (!finc18) finc18 = new TF1("finc18","0.864506+2.79354/sqrt(x)",20,1000);
  if (!finc17) finc17 = new TF1("finc17","0.850069+1.99726/sqrt(x)",20,1000);

  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopSystematics#Luminosity
  double lumi18 = 59.74;
  double lumi17 = 41.53;

  mistag =  (lumi17*finc17->Eval(pt)+lumi18*finc18->Eval(pt)) / (lumi17+lumi18);
  
  return (lumi17*fbcm17->Eval(pt) + lumi18*fbcm18->Eval(pt)) / (lumi17+lumi18);
  //return (lumi17*fbmu17->Eval(pt) + lumi18*fbmu18->Eval(pt)) / (lumi17+lumi18);
} // heavyEffSF
