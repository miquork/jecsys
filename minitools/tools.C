// Some helpful tools
#include "tools.h"

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

// Divide histograms, invoking rebin if needed
TH1D *tools::Divide(const TH1D *_h1, const TH1D *_h2, double c1, double c2, const char *opt) {

  TH1D *h1 = dynamic_cast<TH1D*>(_h1->Clone());
  TH1D *h2 = dynamic_cast<TH1D*>(_h2->Clone());
  TH1D *h1r(0), *h2r(0);
  if (_h1->GetNbinsX()>_h2->GetNbinsX()) {
    h1r = Rebin(h1,h2);
    h1->Delete();
    h1 = h1r;
  }
  if (_h2->GetNbinsX()>_h1->GetNbinsX()) {
    h2r = Rebin(h2,h1);
    h2->Delete();
    h2 = h2r;
  }

  TH1D *h3 = dynamic_cast<TH1D*>(h1->Clone(Form("ratio_%s_%s",h1->GetName(),h2->GetName())));
  h3->Divide(h1, h2, c1, c2, opt);

  return h3;
} // Divide

TH1D *tools::Rebin(TH1D *h, TH1D* href) {
  TH1D *hr(0);
  vector<double> hrefx;
  for (int i = 1; i <= href->GetNbinsX(); ++i)
    hrefx.push_back(href->GetBinLowEdge(i));
  hrefx.push_back(href->GetBinLowEdge(href->GetNbinsX()+1));
  hr = dynamic_cast<TH1D*>(h->Rebin(href->GetNbinsX(),h->GetName(),&hrefx[0]));
  return hr;
}


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
