// Purpose: Recreate Sifu's offset vs pTgen plots from L1RC files
// (e.g. https://indico.cern.ch/event/759372/contributions/3149363/attachments/1721390/2783446/20180924_JERC_Status_2016_2017.pdf, slides 3 and 5)
// Compare 2016 to either RC offset or to 2017 offset
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TF1.h"
#include "TMultiGraph.h"
#include "TLatex.h"
#include "TMinuit.h"
#include "TMatrixD.h"
#include "TMath.h"
#include "TLine.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "tools.h"
#include "tdrstyle_mod14.C"

#include <string>

using namespace std;

FactorizedJetCorrector *getJEC(string version) {

  string s;
  const char *cd = "CondFormats/JetMETObjects/data";
  const char *gt = "Summer16_07Aug2017";

  FactorizedJetCorrector *jecl1;
  s = Form("%s/%s%s.txt",cd,gt,version.c_str());
  cout << s << endl << flush;
  JetCorrectorParameters *l1 = new JetCorrectorParameters(s);
  vector<JetCorrectorParameters> v;
  v.push_back(*l1);
  jecl1 = new FactorizedJetCorrector(v);

  return jecl1;
}

// Calculate corrected jet pT for given pTreco
TF1 *_fCorrPt(0);
FactorizedJetCorrector *_fjec(0);
Double_t fCorrPt(Double_t *x, Double_t *p) {
  double ptreco = x[0];
  double eta = p[0];
  double rho = p[1];
  _fjec->setJetPt(ptreco);
  _fjec->setJetEta(eta);
  _fjec->setJetA(TMath::Pi()*0.4*0.4);
  _fjec->setRho(rho);

  return (_fjec->getCorrection() * ptreco);
}

// Solve JEC (which is a function of pTreco) for given pTgen
double getCorrGen(FactorizedJetCorrector *jec, double ptgen, double eta, double rho) {

  // First invert JEC to get ptreco that corresponds to given ptcorr=ptgen
  _fjec = jec;
  if (_fCorrPt==0) _fCorrPt = new TF1("_fCorrPt",fCorrPt,0,6500.,2);
  _fCorrPt->SetParameters(eta, rho);
  double ptreco = _fCorrPt->GetX(ptgen, max(5.,0.5*ptgen), min(6500./cosh(eta),2*ptgen),
				 1e-4);
  // Then recalculate ptgen=ptcorr in case ptreco hit boundaries above
  double ptcorr = _fCorrPt->Eval(ptreco);
  double corr = ptcorr / ptreco;

  return corr;
}

double getRho(double mu) {
  // Derived from Bahareh's 2016 legacy neutrino gun MC at
  // /Volumes/LaCie_10TB/data/offset/Legacy_tuples/SingleNeutrino_MC.root
  //
  // T->Draw("rho:mu>>hrho(100,0,100)","","");//,100000)
  // TProfile *prho = hrho->ProfileX()
  // TF1 *f1 = new TF1("f1","[0]+[1]*x+[2]*x*x",0,100)
  // prho->Fit(f1,"QRN"); prho->Draw("SAME"); f1->Draw("SAME");
  // for (int i=0; i!=f1->GetNpar(); ++i) cout << ", "<<f1->GetParameters()[i]; cout << endl;
  const double p[3] = {-0.334599, 0.527231, 0.000433922};
  
  return (p[0]+p[1]*mu+p[2]*mu*mu);
}


void drawOffset() {

  setTDRStyle();

  FactorizedJetCorrector *jecrc = getJEC("_V15_MC_L1RC_AK4PFchs");
  FactorizedJetCorrector *jecl1 = getJEC("_V15_MC_L1FastJet_AK4PFchs");

  const double vx[] = //{1, 5, 6, 8, 
    {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832, 4037};
  //, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000};
  const int nx = sizeof(vx)/sizeof(vx[0])-1;
  
  const double veta[] = {0., 0.25, 0.5, 0.75, 1.0, 1.25};
  const int neta = sizeof(veta)/sizeof(veta[0]);

  const double vmu[] = {5, 15, 25, 35, 45, 55};
  const int nmu = sizeof(vmu)/sizeof(vmu[0]); 
  const int color[nmu] = {kRed+1, kOrange+1, kYellow+2, kGreen-9+2, kGreen+2, kCyan+1};

  TH1D *h = new TH1D("h",";p_{T}^{GEN};#LToffset#GT (GeV)",nx,vx);
  h->SetMinimum(-10.);
  h->SetMaximum(+25.);
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();

  extraText = "Simulation";
  lumi_13TeV = "Summer16_07Aug2017_V15_MC";
  TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);
  c1->SetLogx();

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  tex->DrawLatex(0.50,0.87,"Anti-k_{T} R=0.4 PF+CHS");
  tex->DrawLatex(0.50,0.82,"|#eta|<1.3");

  TLegend *leg = tdrLeg(0.40,0.80-0.04*nmu,0.60,0.80);

  for (int imu = 0; imu != nmu; ++imu) {

    double mu = vmu[imu];//40;
    double rho = getRho(mu)+2.00;
      
    TH1D *hrc = new TH1D(Form("hrc_%d",imu),";p_{T}^{GEN} (GeV);#LToffset#GT (GeV)",nx,vx);
    TH1D *hl1 = new TH1D(Form("hl1_%d",imu),";p_{T}^{GEN} (GeV);#LToffset#GT (GeV)",nx,vx);
    
    for (int i = 1 ; i != h->GetNbinsX()+1; ++i) {
      
      double ptgen = h->GetBinCenter(i);
      //double eta = 1.3;
      
      double sumcorrl1(0), sumcorrrc(0), sumw(0);
      for (int ieta = 0; ieta != neta; ++ieta) {
	
	double eta = veta[ieta];
	
	double w = 1;
	double l1p = getCorrGen(jecl1, ptgen, eta, rho);
	double rcp = getCorrGen(jecrc, ptgen, eta, rho);
	sumcorrl1 += l1p*w;
	sumcorrrc += rcp*w;
	sumw += w;
	
	double l1m = getCorrGen(jecl1, ptgen, -eta, rho);
	double rcm = getCorrGen(jecrc, ptgen, -eta, rho);
	sumcorrl1 += l1m*w;
	sumcorrrc += rcm*w;
	sumw += w;
      } // for ieta
      double corrl1 = sumcorrl1/sumw;
      double corrrc = sumcorrrc/sumw;
      
      // corr = (ptreco-off)/ptreco = ptgen/(ptgen+off) => off = ptgen/corr-ptgen
      double offl1 = ptgen*(1./corrl1-1.);
      hl1->SetBinContent(i, offl1);
      double offrc = ptgen*(1./corrrc-1.);
      hrc->SetBinContent(i, offrc);
    } // for i

    tdrDraw(hl1,"HISTP",kFullCircle,color[imu]);//kGreen+2);
    tdrDraw(hrc,"HISTP",kOpenDiamond,color[imu]);//kGreen+2);
    leg->AddEntry(hl1,Form("N_{PU} = %1.0f",mu),"PL");
  } // for imu
  //h->Draw();

  c1->RedrawAxis();

  c1->SaveAs("pdf/drawOffset_AK4PFchs.pdf");
}
