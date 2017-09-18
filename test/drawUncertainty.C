#include "TF1.h" 
#include "TLine.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "tdrstyle_mod14.C"

JetCorrectionUncertainty *_unc(0);
Double_t uncert(Double_t *x, Double_t *par) {

  assert(_unc);

  double pt = *x;
  double etamin = par[0]; assert(etamin==-0.8);
  double etamax = par[1]; assert(etamax==+0.8);

  // Average uncertainty over these etabins (flat weight)
  const double etabins[] = {-0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8};
  const int neta = sizeof(etabins)/sizeof(etabins[0])-1;

  /*
  // Average uncertainty over these pT in these bins
  double ptbins[] =
    {20, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
     97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 362, 430,
     507, 592, 686, 790, 905, 1032, 1172, 1327, 1497, 1684, 1999};
  const int npt = sizeof(ptbins)/sizeof(ptbins[0])-1;
  */

  // Calculate average uncertainty over etamin<eta<etamax
  double sumerr(0);
  for (int ieta = 0; ieta != neta; ++ieta) {
    
    double eta = 0.5*(etabins[ieta]+etabins[ieta+1]);
    _unc->setJetPt(pt);
    _unc->setJetEta(eta);
    double err = _unc->getUncertainty(true);
    double w = 1./neta;
    
    sumerr += w * err;
  } // for ieta
  
  return sumerr;
}

void drawUncertainty(string subset = "TotalNoFlavorNoTime",
		     string algo = "AK5PFchs.txt") {

  setTDRStyle();

  // Path and uncertainty source file name 
  const char *a = algo.c_str();
  string ssf = Form("CondFormats/JetMETObjects/data/Winter14_V8_DATA_UncertaintySources_%s.txt",a);
  const char *sf = ssf.c_str();
    
  // Source name
  const char *src = subset.c_str();
  cout << sf << ":" << src << endl << flush;
  JetCorrectorParameters *p = new JetCorrectorParameters(sf, src);
  JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
  _unc = unc;

  TF1 *f1 = new TF1("f1",uncert,10,2000,2);
  f1->SetParameters(-0.8, 0.8);
  f1->SetNpx(1000);
  f1->Draw();

  cout << " Min(0.32%): " << f1->GetX(0.0032,100,250)
       << " Max(0.32%): " << f1->GetX(0.0032,250,400) << endl;

}
