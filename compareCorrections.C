
#include <cassert>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include "TROOT.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"
#include "Math/IFunction.h"
#include "Math/BrentRootFinder.h"
#include "Math/RootFinderAlgorithms.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

const double x_pt[] =
    {8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
     97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 362, 430,
     507, 592, 686, 790, 905, 1032, 1172, 1327, 1497, 1684, 1890, //1999};
     2000, 2238, 2500, 2787, 3103, 3450};
const int ndiv_pt = sizeof(x_pt)/sizeof(x_pt[0])-1;
const double x_eta[] =  {-5.4,-5.0,-4.4,-4,-3.5,-3,-2.8,-2.6,-2.4,-2.2,-2.0,-1.8,-1.6,-1.4,-1.2,-1.0, -0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3,3.5,4,4.4,5.0,5.4};
const int ndiv_eta = sizeof(x_eta)/sizeof(x_eta[0])-1;

int npv = 14;

bool plotraw = true;//false

// For JEC central value
/*
gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+");
gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+");
gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+");
*/


double rho(const double npvmean) {

  // For data 2012-10-24 (9.2/fb 53X Fall12 V1; range 0.5-40.5)
  //const float _uepf = 1.068; // from 2010 studies
  const float _ootpf = 1.606; // UE+OOT GeV/A 2012
  const float _itpf1 = 0.745; // IT p1 GeV/A 2012
  const float _itpf2 = +0.0016; // IT p2 GeV/A 2012
  double rho = _ootpf + (npvmean-1) * (_itpf1 + (npvmean-1) * _itpf2);

  return rho;
}

class ResponseFunc : public ROOT::Math::IBaseFunctionOneDim
{
public:
  ResponseFunc(double pTprime, FactorizedJetCorrector *jec, int npv, double eta, double rho, double jeta) :
    ROOT::Math::IBaseFunctionOneDim(),_pTprime(pTprime),_jec(jec),_npv(npv),_eta(eta),_rho(rho),_jeta(jeta) {}

  double DoEval(double pTraw) const {
    _jec->setJetPt(pTraw);
    _jec->setJetE(pTraw*cosh(_eta));
    _jec->setJetEta(_eta);
    _jec->setNPV(_npv);
    _jec->setRho(_rho); _jec->setJetA(_jeta);
    double cor = _jec->getCorrection();
    //std::cout << "in Brent: pTraw:" << pTraw << "  cor:" << cor << "  dist:" << pTraw * cor - _pTprime << '\n';
    return pTraw * cor - _pTprime;
  }

  ROOT::Math::IBaseFunctionOneDim* Clone() const {
    return new ResponseFunc(_pTprime,_jec,_npv,_eta,_rho,_jeta);
  }
private:
  double _pTprime;
  FactorizedJetCorrector *_jec;
  int _npv;
  double _eta;
  double _rho;
  double _jeta;
};

// Solve pTraw from pTprime = pTraw / R(pTraw) using Brent's method
// We want to provide JEC uncertainties vs pTprime, not pTraw, but JEC
// is only available as a function of pTraw
double findPtraw(double ptref, double eta, double jetarea, FactorizedJetCorrector* jec) {


  ResponseFunc f(ptref,jec,npv,eta,rho(npv),jetarea);

  ROOT::Math::Roots::Brent brf;
  //brf.SetLogScan(true);
  // Set parameters of the method
  brf.SetFunction(f,std::max(2.0,0.25*ptref),std::min(4*ptref,3500.0));
  bool found_root = brf.Solve(50,1e-4,1e-5);
  double ptraw = brf.Root();
  double rjet = ptraw /ptref;
  jec->setJetE(ptraw*cosh(eta));
  jec->setJetEta(eta);
  jec->setNPV(npv);
  jec->setRho(rho(npv));
  jec->setJetA(jetarea);
  jec->setJetPt(ptraw);
  double corr = jec->getCorrection();
  if(std::abs((ptraw * corr - ptref)/ptref) > 0.001) {
    std::cout << "NPV:" << npv << "   Brent: status:" << brf.Status() << " " << brf.Iterations() << '\n';
    std::cout << "R: pTprime:" << ptref << "  eta:" << eta << " pTraw:" << ptraw << " corr:" << corr
              << " pTcor:" << corr * ptraw << " rjet*pTprime:" << rjet * ptref << '\n';
  }
  assert(std::abs((ptraw * corr - ptref)/ptref) < 0.001);
  assert(found_root);
  assert(brf.Status() == 0);
  /*
  if((eta > 2.6) && (eta < 2.8)) {
    _jec->setJetE(pTprime*cosh(eta));
    _jec->setJetEta(eta); _jec->setNPV(_npv);
    _jec->setJetPt(pTprime);
    std::cout << "NPV:" << _npv << " cor:" << _jec->getCorrection() << "   Brent: status:" << brf.Status() << " " << brf.Iterations() << '\n';
    std::cout << "R: pTprime:" << pTprime << "  eta:" << eta << " pTraw:" << pTraw << " corr:" << corr
              << " pTcor:" << corr * pTraw << " rjet*pTprime:" << _rjet * pTprime << '\n';
  }
  */
  return ptraw;
} // _Rjet

TProfile *GetProfile(TH2D* hd2, double eta1, double eta2) {
  int bin1 = hd2->GetYaxis()->FindBin(eta1);
  int bin2 = hd2->GetYaxis()->FindBin(eta2);

  //TProfile *hp = new TProfile("hp","",hd2->GetXaxis()->GetNbins(),hd2->GetXaxis()->GetXbins()->GetArray());
  TProfile* hp = hd2->ProfileX();
  hp->Reset();
  for(int i = 0 ; i < hp->GetNbinsX()+1 ; ++i) {
    double pt = hp->GetBinCenter(i);
    for(int j = bin1 ; j <=bin2 ; ++j) { 
      double eta = hd2->GetYaxis()->GetBinCenter(j);
      if(pt*cosh(eta) <= 3500) {
	hp->Fill(pt,hd2->GetBinContent(i,j));
      }
    }
  }
  hp->SetStats(0);
  return hp;
}

TProfile *GetProfile2(TH2D* hd2, double pt1, double pt2) {
  int bin1 = hd2->GetXaxis()->FindBin(pt1);
  int bin2 = hd2->GetXaxis()->FindBin(pt2);

  //TProfile *hp = new TProfile("hp","",hd2->GetXaxis()->GetNbins(),hd2->GetXaxis()->GetXbins()->GetArray());
  TProfile* hp = hd2->ProfileY();
  hp->Reset();
  for(int i = 0 ; i < hp->GetNbinsX()+1 ; ++i) {
    double eta = hp->GetBinCenter(i);
    for(int j = bin1 ; j <=bin2 ; ++j) { 
	hp->Fill(eta,hd2->GetBinContent(j,i));
    }
  }
  hp->SetStats(0);
  return hp;
}


TH2D* getHist(const std::string& jec,const std::string& algo, double offsetPt = 0) {
  double eta1, eta2,pt, val1, val2;
  int nbins;
  double r = 0;
  assert((algo[2]=='5') || (algo[2]=='7'));
  if(algo[2] == '5') r = 0.5;
  if(algo[2] == '7') r = 0.7;
  std::vector<JetCorrectorParameters> params;
  //std::cout << jec+"_L1FastJet_"+algo+".txt" << '\n';
  params.push_back(jec+"_L1FastJet_"+algo+".txt");
  params.push_back(jec+"_L2Relative_"+algo+".txt");
  params.push_back(jec+"_L3Absolute_"+algo+".txt");
  params.push_back(jec+"_L2L3Residual_"+algo+".txt");
  FactorizedJetCorrector* cor = new FactorizedJetCorrector(params);
  std::string title = jec+" "+algo;
  int ptstartbin = plotraw ? 2 : 0;
  TH2D* hist = new TH2D(title.c_str(),title.c_str(),ndiv_pt-ptstartbin, x_pt+ptstartbin,ndiv_eta,x_eta);

  if(plotraw) {
      hist->SetXTitle("p_{T}^{raw} [GeV]");
  } else {
      hist->SetXTitle("p_{T}^{ref} [GeV]");
  }

  hist->SetYTitle("#eta");
  hist->SetStats(0);
  hist->GetXaxis()->SetMoreLogLabels();
  hist->GetXaxis()->SetNoExponent();
  hist->SetMinimum(0.5);
  //hist->SetMaximum(1.8);
  std::cout << "initialzing jet corrections " << jec << "for " << algo << '\n';

  for(unsigned int i = 0 ; i < ndiv_eta ; ++i) {
   double eta = (x_eta[i]+x_eta[i+1])/2;
   if(eta >= 5.2) eta = 5.1;
   if(eta <= -5.2) eta = -5.1;
   double cosheta = cosh(eta);
    for(unsigned int j = ptstartbin ; j < ndiv_pt ; ++j) {
      double pt = (x_pt[j]+x_pt[j+1])/2;
      if( pt*cosheta <= 3500) {
	double ptraw = plotraw ? pt : findPtraw(pt, eta, TMath::Pi()*r*r, cor );
	cor->setJetE(ptraw*cosheta);
	cor->setJetPt(ptraw);
	cor->setJetEta(eta);
	cor->setJetPhi(0);
	cor->setNPV(npv);
	cor->setRho(rho(npv));
	cor->setJetA(TMath::Pi()*r*r);
	double c = cor->getCorrection();
	c+= offsetPt/ptraw;
	//std::cout << pt << ", " << eta << ", " << c << '\n';
	hist->Fill(pt,eta,c);
      } else {
	//hist->Fill(pt,eta,0);
      }
    }
  }

  delete cor;

  return hist;

}


void compareCorrections() {
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadLeftMargin(0.10);
  gStyle->SetPadRightMargin(0.10);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetTitleOffset(1.10,"y");
  gStyle->SetOptTitle(1);

  std::vector<std::string> list;

  list.push_back("AK5PFchs");
  //list.push_back("AK7PFchs");
  //list.push_back("AK5PF");
  //list.push_back("AK7PF");
  //list.push_back("AK5Calo");
  //list.push_back("AK7Calo");
  //list.push_back("AK5JPT");

  std::string newJEC = "/afs/desy.de/user/s/stadie/xxl/jetcorrections/Winter14_V7_DATA";
  //std::string oldJEC = "/afs/desy.de/user/s/stadie/xxl/jetcorrections/Winter14_VDENIS_DATA";
  //std::string oldJEC = "/afs/desy.de/user/s/stadie/xxl/jetcorrections/Winter14_V2PT_DATA";
  std::string oldJEC = "/afs/desy.de/user/s/stadie/xxl/jetcorrections/FT53_V21A_AN6";
  //std::string oldJEC = "/afs/desy.de/user/s/stadie/xxl/jetcorrections/PTFIXV2_FT_53_V21_AN5_private";
  //std::string oldJEC = "/afs/desy.de/user/s/stadie/xxl/jetcorrections/Winter14_V1_MC";
  //std::string newJEC = "/afs/desy.de/user/s/stadie/xxl/jetcorrections/START53_V27";
  for(int i = 0 ; i < list.size() ; ++i) {
    TCanvas* c = new TCanvas();
    TH2D* h2 = getHist(newJEC,list[i]);
    h2->DrawCopy("COLZ");
    c->SetLogx();
    //c->SetLogz();
    std::string name = list[i];
    std::string pname = "new"+list[i]+".eps";
    c->Print(pname.c_str());


    std::string oldname = "old"+list[i];
    TH2D* h2old = getHist(oldJEC,list[i]);
    if(h2old) {
      TCanvas* c = new TCanvas();
      h2old->DrawCopy("COLZ");
      c->SetLogx();
      //c->SetLogz();
      pname = oldname+".eps";
      c->Print(pname.c_str());
      c = new TCanvas();
      pname = name+"diff";
      TH2D* hdiff = (TH2D*)h2->Clone(pname.c_str());
      pname = name+" (new-old)/old";
      hdiff->SetTitle(pname.c_str());
      hdiff->Add(h2old,-1);
      hdiff->Divide(h2old);
      //hdiff->SetMaximum(1.10);
      //hdiff->SetMinimum(0.90);
      hdiff->DrawCopy("COLZ");
      c->SetLogx();
      pname = name+"ratio"+".eps";
      c->Print(pname.c_str());
      TProfile* hpx = GetProfile(hdiff,-1.3,1.3);
      c = new TCanvas();
      c->SetLogx();
      hpx->DrawCopy();
      pname = name+"ratiobarrel"+".eps";
      c->Print(pname.c_str()); 
      TProfile* hpy = GetProfile2(hdiff,100,100);
      c = new TCanvas();
      //c->SetLogx();
      hpy->SetMaximum(0.2);
      hpy->SetMinimum(-0.2);
      hpy->DrawCopy(); 
      pname = name+"ratiobarrel2"+".eps";
      c->Print(pname.c_str());
      delete hpx;
      delete hpy;
      delete h2old;
      delete hdiff;

    }
    delete h2;
  }
  /*
  TH2D* hold = readHist("txt/Summer13_V1_DATA_Uncertainty_AK5Calo.txt","AK5Calo old");
  TH2D* hnew = readHist("txt/Summer13_V4_DATA_Uncertainty_AK5Calo.txt","differences (old-new)");

  TCanvas* c = new TCanvas();
  hold->DrawClone("COLZ");
  c->SetLogx();
  c->Print("oldAK5Calo.eps");

  TH2D* hdiff = (TH2D*)hold->Clone("hdiffcalo");
  hdiff->SetTitle("AK5Calo differences (old-new)");
  hdiff->Add(hnew, -1.0);
  //hold->Scale(-1.0);
  c = new TCanvas();
  hdiff->SetMaximum(0.2);
  hdiff->SetMinimum(-0.0);
  hdiff->DrawClone("COLZ");
  c->SetLogx();
  //c->SetLogz();
  c->Print("diffchsv3calo.eps");
  c = new TCanvas();
  hdiff->Divide(hold);
  hdiff->SetTitle("AK5Calo differences (old-new)/old");
  hdiff->SetMaximum(0.2);
  hdiff->SetMinimum(-0.1);
  hdiff->DrawClone("COLZ");
  c->SetLogx();
  c->Print("diffchsv3calorel.eps");
  */
}
