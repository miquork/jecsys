#ifndef __tools_h__
#define __tools_h__

#include "TChain.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TF1.h"

#include <string>
#include <vector>

namespace tools {

  int addFiles(TChain *c, std::string filelistname);
  double delta_phi(double a, double b);
  double delta_eta(double a, double b);
  double oplus(double a, double b);
  void swap(double &a, double &b);

  // vector manipulation
  std::vector<double> make_vector(double *a, int na);
  double interpolate(double x, std::vector<double> const& vx,
		     std::vector<double> const& vy); 

  // Graph manipulation
  void GetPoint(TGraphErrors *g, int n, double &x, double &y,
		double &ex, double &ey);
  void SetPoint(TGraphErrors *g, int n, double x, double y,
		double ex, double ey);
  TGraphErrors *makeGraph(TH1 *hx, TH1 *hy, double scale = 1.);
  TGraphErrors *ratioGraphs(TGraphErrors *g1, TGraphErrors *g2, double erry=1.);
  TGraphErrors *ratioGraphs(TGraphErrors *g1, TF1 *f2);
  int findPoint(TGraph *g, double x);

  // Histogram manipulation
  TH1D *Divide(const TH1D *h1, const TH1D *h2, double c1=1, double c2=1,
	       const char *opt="");
  TH1D *Rebin(const TH1D *h, const TH1D* href);

  void Hadd(TH1 *h1, TH1 *h2, double ptmax=0, bool syserr = false);
} // namespace tools

#endif
// __tools_h__
