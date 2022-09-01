// Some helpful tools
#ifndef __TOOLSC__
#define __TOOLSC__
#include "tools.h"

#include "TChain.h"
#include "TMath.h"

#include <fstream>
#include <iostream>
#include <string>

using namespace std;
using namespace tools;

// Add a list of files to a TChain
int tools::addFiles(TChain *c, string filelistname) {
  
  ifstream filelist(filelistname.c_str(), ios::in);
  string filename;
  int igood(0), ibad(0);
  
  while (filelist >> filename) {
    
    if (c->AddFile(filename.c_str())) ++igood;
    else ++ibad;
  }
  
  cout << "Loaded " << igood << " files." << endl; 
  if (ibad!=0)
    cerr << "Warning: failed to load " << ibad << " files  out of "
	 << igood << endl;
  
  return igood;
}

double tools::delta_phi(double a, double b) {
  double dphi = a - b;
  if (dphi > TMath::Pi()) dphi -= TMath::TwoPi();
  if (dphi < -TMath::Pi()) dphi += TMath::TwoPi();
  return fabs(dphi);
}

double tools::delta_eta(double a, double b) {
  return (a - b);
}

double tools::oplus(double a, double b) {
  return sqrt(a*a + b*b);
}

void tools::swap(double &a, double &b) {
  double c = a;
  a = b;
  b = c;
}

// Vector manipulations
vector<double> tools::make_vector(double *a, int na) {
  
  vector<double> *vx = new vector<double>(na);
  for (int i = 0; i != na; ++i) (*vx)[i] = a[i];

  return (*vx);
} // make_vector
// Interpolation between bin centers
double tools::interpolate(double x, vector<double> const& vx,
			  vector<double> const& vy) {

  assert(vx.size()==vy.size());
  assert(vx.size()>=2);
  int i = TMath::BinarySearch(vx.size(), &vx[0], x);
  // Fixed value at the edges
  if (i==-1) return vy[0];
  if (i==int(vx.size())-1) return vy[vx.size()-1];

  double x1 = vx[i];
  double x2 = vx[i+1];
  double y1 = vy[i];
  double y2 = vy[i+1];
  assert(x1!=x2);

  return (y1 + (y2-y1)/(x2-x1)*(x-x1));
} // interpolate

// Ratio of two graphs
void tools::GetPoint(TGraphErrors *g, int n, double &x, double &y,
		     double &ex, double &ey) {
  g->GetPoint(n, x, y);
  ex = g->GetErrorX(n);
  ey = g->GetErrorY(n);
}
void tools::SetPoint(TGraphErrors *g, int n, double x, double y,
	      double ex, double ey) {
  g->SetPoint(n, x, y);
  g->SetPointError(n, ex, ey);
}
TGraphErrors *tools::makeGraph(TH1 *hx, TH1 *hy, double scale ) {

  assert(hx->GetNbinsX()==hy->GetNbinsX());

  TGraphErrors *g = new TGraphErrors(0);
  for (int i = 1; i != hx->GetNbinsX()+1; ++i) {

    assert(hx->GetBinLowEdge(i)==hy->GetBinLowEdge(i));
    int n = g->GetN();
    SetPoint(g, n, hx->GetBinContent(i), scale*hy->GetBinContent(i),
	     hx->GetBinError(i), scale*hy->GetBinError(i));
  }

  return g;
}
TGraphErrors *tools::diffGraphs(TGraphErrors *g1, TGraphErrors *g2,
				double c1, double c2) {

  assert(g1); assert(g2);
  TGraphErrors *g = new TGraphErrors(0);
  g->SetName(Form("diff_%s_%s",g1->GetName(),g2->GetName()));
  g->SetMarkerStyle(g1->GetMarkerStyle());
  g->SetMarkerColor(g1->GetMarkerColor());
  g->SetLineColor(g1->GetLineColor());

  // take diff only of points that are closer to each other than any others
  for (int i = 0, j = 0; (i != g1->GetN() && j != g2->GetN());) {
    
    double x1m, x1p, x2m, x2p, y, ex, ey;
    GetPoint(g1, max(i-1,0), x1m, y, ex, ey);
    GetPoint(g1, min(i+1,g1->GetN()-1), x1p, y, ex, ey);
    GetPoint(g2, max(j-1,0), x2m, y, ex, ey);
    GetPoint(g2, min(j+1,g2->GetN()-1), x2p, y, ex, ey);
    double x1, y1, ex1, ey1;
    GetPoint(g1, i, x1, y1, ex1, ey1);
    double x2, y2, ex2, ey2;
    GetPoint(g2, j, x2, y2, ex2, ey2);

    if (fabs(x1-x2)<=fabs(x1-x2m) && fabs(x1-x2)<=fabs(x1-x2p) &&
	fabs(x1-x2)<=fabs(x2-x1m) && fabs(x1-x2)<=fabs(x2-x1p)) {

      int n = g->GetN();
      SetPoint(g, n, 0.5*(x1+x2), c1*y1-c2*y2, 0.5*fabs(x1-x2),
	       oplus(c1*ey1, c2*ey2));
      ++i, ++j;
    }
    else
      (x1<x2 ? ++i : ++j);
  } // for i, j

  return g;
} // ratioGraphs
TGraphErrors *tools::ratioGraphs(TGraphErrors *g1, TGraphErrors *g2,
				 double erry) {

  assert(g1); assert(g2);
  TGraphErrors *g = new TGraphErrors(0);
  g->SetName(Form("ratio_%s_%s",g1->GetName(),g2->GetName()));
  g->SetMarkerStyle(g1->GetMarkerStyle());
  g->SetMarkerColor(g1->GetMarkerColor());
  g->SetLineColor(g1->GetLineColor());

  // take ratio only of points that are closer to each other than any others
  for (int i = 0, j = 0; (i != g1->GetN() && j != g2->GetN());) {
    
    double x1m, x1p, x2m, x2p, y, ex, ey;
    GetPoint(g1, max(i-1,0), x1m, y, ex, ey);
    GetPoint(g1, min(i+1,g1->GetN()-1), x1p, y, ex, ey);
    GetPoint(g2, max(j-1,0), x2m, y, ex, ey);
    GetPoint(g2, min(j+1,g2->GetN()-1), x2p, y, ex, ey);
    double x1, y1, ex1, ey1;
    GetPoint(g1, i, x1, y1, ex1, ey1);
    double x2, y2, ex2, ey2;
    GetPoint(g2, j, x2, y2, ex2, ey2);

    if (fabs(x1-x2)<=fabs(x1-x2m) && fabs(x1-x2)<=fabs(x1-x2p) &&
	fabs(x1-x2)<=fabs(x2-x1m) && fabs(x1-x2)<=fabs(x2-x1p)) {

      int n = g->GetN();
      SetPoint(g, n, 0.5*(x1+x2), (y2 ? y1/y2 : 0.), 0.5*fabs(x1-x2),
	       oplus(ey1/y1, ey2/y2) * fabs(y2 ? y1/y2 : 0.) * erry); 
      ++i, ++j;
    }
    else
      (x1<x2 ? ++i : ++j);
  } // for i, j

  return g;
} // ratioGraphs
TGraphErrors *tools::ratioGraphs(TGraphErrors *g1, TF1 *f2) {

  assert(g1); assert(f2);
  TGraphErrors *g = new TGraphErrors(0);
  g->SetMarkerStyle(g1->GetMarkerStyle());

  for (int i = 0; i != g1->GetN(); ++i) {

    double x, y, ex, ey;
    GetPoint(g1, i, x, y, ex, ey);
    double v = f2->Eval(x);
    SetPoint(g, i, x, v ? y/v : 0, ex, v ? ey/v : 0.);
  }

  return g;
} // ratioGraphs (TF1*)
int tools::findPoint(TGraph *g, double x) {

  int k = 0;
  double dxmin = -1;
  for (int i = 0; i != g->GetN(); ++i) {
    
    double dx = fabs(g->GetX()[i]-x);
    if (dx<=dxmin || dxmin<0) {
      dxmin = dx;
      k = i;
    }
  } // for i
  
  return k;
} // findPoint(TGraph*)

// combine points statistically
TGraphErrors *tools::mergeGraphs(TGraphErrors *g1, TGraphErrors *g2) {
  // find j closests to i, and i closest to j
  // when both are the same, pair points
  map<int,int> itoj, jtoi;
  map<int,double> idxmin, jdxmin;
  for (int i = 0; i != g1->GetN(); ++i) {
    for (int j = 0; j != g2->GetN(); ++j) {

      double dx = fabs(g1->GetX()[i]-g2->GetX()[j]);
      if (dx<idxmin[i] || j==0) { itoj[i] = j; idxmin[i] = dx; }
      if (dx<jdxmin[j] || i==0) { jtoi[j] = i; jdxmin[j] = dx; }
    }
  }

  TGraphErrors *g = new TGraphErrors(0);
  for (int i = 0; i != g1->GetN(); ++i) {
    int j = itoj[i];
    //cout <<Form("i %i j %i itoj[i] %i",i,j, itoj[i]) << endl;
    //cout << Form("g1->GetN() %i g2->GetN() %i", g1->GetN(), g2->GetN()) << endl;
    if (jtoi[j]==i && g2->GetN()>0) {
      int n = g->GetN();
      // 2021-03-11: logic may have been wrong here
      //double w1 = g2->GetEY()[j] / (g2->GetEY()[j] + g1->GetEY()[i]);
      //double w2 = g1->GetEY()[i] / (g2->GetEY()[j] + g1->GetEY()[i]);
      // Updated logic: 1/err^2 ~ N, since err ~ 1/sqrt(N)
      double w1(0.5), w2(0.5);
      // 2021-06-14: w1==1 and w2==0 change to w1=1 and w2=0 etc.
      if      (g2->GetN()==0) { w1=1; w2=0; }
      else if (g1->GetN()==0) { w1=0; w2=1; }
      else {
	double n1 = 1./pow(g1->GetEY()[i],2);
	double n2 = 1./pow(g2->GetEY()[j],2);
	w1 = n1 / (n1+n2);
	w2 = n2 / (n1+n2);
      }
      double y = (g1->GetY()[i]*w1 + g2->GetY()[j]*w2);
      double dy = sqrt(pow(g1->GetEY()[i]*w1,2)+pow(g2->GetEY()[j]*w2,2));
      double x = (g1->GetX()[i]*w1 + g2->GetX()[j]*w2);
      double dx = sqrt(pow(g1->GetEX()[i]*w1,2)+pow(g2->GetEX()[j]*w2,2));
      g->SetPoint(n, x, y);
      g->SetPointError(n, dx, dy);
    }
  }

  return g;
} // mergeGraphs

// Divide histograms, invoking rebin if needed
TH1D *tools::Divide(const TH1D *h1, const TH1D *h2, double c1, double c2,
		    const char *opt) {

  TH1D *h1r(0), *h2r(0);
  if (h1->GetNbinsX()>h2->GetNbinsX()) {
    h1r = Rebin(h1, h2); h1 = h1r;
  }
  if (h2->GetNbinsX()>h1->GetNbinsX()) {
    h2r = Rebin(h2, h1); h2 = h2r;
  }
  
  TH1D *h3 = (TH1D*)h1->Clone(Form("ratio_%s_%s",h1->GetName(),h2->GetName()));
  h3->Divide(h1, h2, c1, c2, opt);

  // delete temporary copies
  if (h1r) h1r->Delete();
  if (h2r) h2r->Delete();
  
  return h3;
} // Divide
// Rebin first histogram to match the second
TH1D *tools::Rebin(const TH1D *h, const TH1D* href) {

  //assert(href->GetNbinsX()<=h->GetNbinsX());
  if (!(href->GetNbinsX()<=h->GetNbinsX())) {
    cout << "Histo has less bins than ref: "
	 << h->GetNbinsX() << " vs " << href->GetNbinsX()
	 << " for " << h->GetName() << endl;
  }
  
  // First, we need to rebin inclusive jets to match b-tagged jets
  TH1D *hre = (TH1D*)href->Clone(Form("%s_rebin",h->GetName()));
  hre->Reset();

  for (int i = 1; i != h->GetNbinsX()+1; ++i) {

    double x = h->GetBinLowEdge(i);
    int j = hre->FindBin(x);
    // Check that h is fully contained within href bin
    if (h->GetBinContent(i)!=0) {

      if (!(h->GetBinLowEdge(i)>=hre->GetBinLowEdge(j) - 1e-5 &&
	    h->GetBinLowEdge(i+1)<=hre->GetBinLowEdge(j+1) + 1e-5)) {
	cerr << Form("Warning, bin edges overlapping: h=[%1.0f,%1.0f],"
		     " hre=[%1.0f,%1.0f] (%s)",
		     h->GetBinLowEdge(i), h->GetBinLowEdge(i+1),
		     hre->GetBinLowEdge(j), hre->GetBinLowEdge(j+1),
		     h->GetName()) << endl;
      }

      double y = ( hre->GetBinContent(j)*hre->GetBinWidth(j)
		   + h->GetBinContent(i)*h->GetBinWidth(i) )
	/ hre->GetBinWidth(j);
      //double ey = ( hre->GetBinError(j)*hre->GetBinWidth(j)
      //	    + h->GetBinError(i)*h->GetBinWidth(i) )
      // / hre->GetBinWidth(j);
      double ey = sqrt( pow(hre->GetBinError(j)*hre->GetBinWidth(j),2)
			+ pow(h->GetBinError(i)*h->GetBinWidth(i),2) )
	/ hre->GetBinWidth(j);
      hre->SetBinContent(j, y);
      hre->SetBinError(j, ey);
    }
  } // for i

  return hre;
} // Rebin


// Add two histograms by averaging
void tools::Hadd(TH1 *h1, TH1 *h2, double ptmax, bool syserr) {

  //assert(h1->GetNbinsX()==h2->GetNbinsX());

  for (int i = 1; i != h1->GetNbinsX()+1; ++i) {

    double y1 = h1->GetBinContent(i);
    double ey1 = h1->GetBinError(i);
    int j = h2->FindBin(h1->GetBinCenter(i));
    if (j>0 && j<h2->GetNbinsX()+1 &&
	(ptmax==0 || h1->GetBinLowEdge(i+1)<=ptmax)) {

      double y2 = h2->GetBinContent(j);
      double ey2 = h2->GetBinError(j);
      
      //if (y1!=0 && y2!=0 && ey1!=0 && ey2!=0) {
      
      double n1 = (ey1 ? 1./pow(ey1,2) : 0.);
      double n2 = (ey2 ? 1./pow(ey2,2) : 0.);
      if (y1!=0 && ey1==0) n1 = 1;
      if (y2!=0 && ey2==0) n2 = 1;
      double y = (n1+n2 ? (y1*n1 + y2*n2)/(n1+n2) : 0.);
      double ey = (n1+n2 ? y/sqrt(n1+n2) : 0.);
      if (syserr) ey = (n1+n2 ? (ey1*n1 + ey2*n2)/(n1+n2) : 0.);
      h1->SetBinContent(i, y);
      h1->SetBinError(i, ey);
      //}
    }
  } // for i

} // Hadd


// Generic error calculation using differentials
void tools::drawErrBand(TF1 *f1, TFitResultPtr &fp, double xmin, double xmax) {

  TMatrixDSym emat = fp->GetCovarianceMatrix();

  //TH1D *he = new TH1D(Form("he_%d",(int)fp),"",1,xmin,xmax);
  TGraphErrors *ge = new TGraphErrors(int(xmax-xmin)+1);
  TGraph *ge_up = new TGraph(int(xmax-xmin)+1);
  TGraph *ge_dw = new TGraph(int(xmax-xmin)+1);
  vector<double> dfdp(f1->GetNpar());
  //vector<double> dp(f1->GetNpar());
  for (int i = 0; i != ge->GetN(); ++i) {

    double x = xmin + (xmax-xmin)*i/(ge->GetN()-1);
    double y = f1->Eval(x);
    for (int j = 0; j != f1->GetNpar(); ++j) {

      double p = f1->GetParameter(j);
      double ep = f1->GetParError(j);
      double dp = 0.1*ep;
      f1->SetParameter(j, p+0.1*ep);
      double yup = f1->Eval(x);
      f1->SetParameter(j, p-0.1*ep);
      double ydw = f1->Eval(x);
      f1->SetParameter(j, p);
      dfdp[j] = (dp!=0 ? 0.5*(yup-ydw)/dp : 0);
      //dp[j] = ep;
    } // for j
    
    double e2(0);
    for (int j = 0; j != f1->GetNpar(); ++j) {
      for (int k = 0; k != f1->GetNpar(); ++k) {
	e2 += dfdp[j]*dfdp[k]*emat[j][k];
      } // for k
    } //  for j
    
    double err = sqrt(e2);
    ge->SetPoint(i, x, y);
    ge->SetPointError(i, 0, err);
    ge_up->SetPoint(i, x, y+err);
    ge_dw->SetPoint(i, x, y-err);
  } // for i

  ge->SetFillColorAlpha(kBlue-9,0.5);
  ge->SetFillStyle(1001);
  ge->Draw("SAMEE3");
  ge_up->SetLineStyle(kDotted);
  ge_up->Draw("SAMEL");
  ge_dw->SetLineStyle(kDotted);
  ge_dw->Draw("SAMEL");

} // drawErr
#endif
