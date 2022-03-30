// Purpose: Create new input L2L3Res text file with weighted inputs
//          Uses simple effective parameterization for the output file
//          same way as done in minitools/createL2L3ResTextFile.C
//          Prototype is only setup to merge either years or all Run2 IOVs
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>

using namespace std;

const bool debug = true;//false;

void mergeL2L3ResTextFiles(string run = "Run2", double dz = 0.002) {

  // List IOVs to be merged and their luminosities
  vector<pair<string,double> > iovs;
  if (run=="2016APV" || run=="Run2") {
    iovs.push_back(make_pair<string,double>("2016BCD",12.9));
    iovs.push_back(make_pair<string,double>("2016EF",6.8));
  }
  if (run=="2016GH" || run=="Run2") {
    iovs.push_back(make_pair<string,double>("2016GH",16.8));
  }

  if (run=="2017" || run=="Run2") {
    iovs.push_back(make_pair<string,double>("2017B",4.8));
    iovs.push_back(make_pair<string,double>("2017C",9.6));
    iovs.push_back(make_pair<string,double>("2017D",4.2));
    iovs.push_back(make_pair<string,double>("2017E",9.3));
    iovs.push_back(make_pair<string,double>("2017F",13.4));
  }

  if (run=="2018" || run=="Run2") {
    iovs.push_back(make_pair<string,double>("2018A",14.0));
    iovs.push_back(make_pair<string,double>("2018B",7.1));
    iovs.push_back(make_pair<string,double>("2018C",6.9));
    iovs.push_back(make_pair<string,double>("2018D",31.9));
  }

  double totlumi(0);
  for (int iov = 0; iov != iovs.size(); ++iov) totlumi += iovs[iov].second;
  
  // Helpers to map to file name
  map<string,const char*> mrun;
  mrun["2016BCD"] = "BCD";
  mrun["2016EF"] = "EF";
  mrun["2016GH"] = "FGH";

  mrun["2017B"] = "B";
  mrun["2017C"] = "C";
  mrun["2017D"] = "D";
  mrun["2017E"] = "E";
  mrun["2017F"] = "F";

  mrun["2018A"] = "A";
  mrun["2018B"] = "B";
  mrun["2018C"] = "C";
  mrun["2018D"] = "D";

  map<string,const char*> mver;
  mver["2016BCD"] = "V7";
  mver["2016EF"] = "V7";
  mver["2016GH"] = "V7";

  mver["2017B"] = "V5";
  mver["2017C"] = "V5";
  mver["2017D"] = "V5";
  mver["2017E"] = "V5";
  mver["2017F"] = "V5";

  mver["2018A"] = "V5";
  mver["2018B"] = "V5";
  mver["2018C"] = "V5";
  mver["2018D"] = "V5";


  map<string,const char*> mera;
  mera["2016BCD"] = "Summer19UL16APV";
  mera["2016EF"] = "Summer19UL16APV";
  mera["2016GH"] = "Summer19UL16";

  mera["2017B"] = "Summer19UL17";
  mera["2017C"] = "Summer19UL17";
  mera["2017D"] = "Summer19UL17";
  mera["2017E"] = "Summer19UL17";
  mera["2017F"] = "Summer19UL17";

  mera["2018A"] = "Summer19UL18";
  mera["2018B"] = "Summer19UL18";
  mera["2018C"] = "Summer19UL18";
  mera["2018D"] = "Summer19UL18";

  // pt bins from inclusive jets
  const double vpt[] =
  //1, 5, 6, 8,
    {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133,
     153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592,
     638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410,
     1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787,
     2941, 3103, 3273, 3450, 3637, 3832, 4037, 4252, 4477, 4713, 4961, 5220,
     5492, 5777, 6076, 6389, 6717};
  //, 7000};
  const int npt = sizeof(vpt)/sizeof(vpt[0])-1;

  // eta bins with 'cat X.txt  | awk '{print $1", "}''
  double veta[] =
    {-5.191, -3.839, -3.489, -3.139, -2.964, -2.853, -2.650, -2.500,
     -2.322, -2.172, -1.930, -1.653, -1.479, -1.305, -1.044,
     -0.783, -0.522, -0.261, 0.000, 0.261, 0.522, 0.783, 1.044,
     1.305, 1.479, 1.653, 1.930, 2.172, 2.322, 2.500,
     2.650, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};
  const int neta = sizeof(veta)/sizeof(veta[0])-1;
  // pt range from text file, pt2max with awk '{print $7", "}'
  double minpt3 = 10;
  double maxpt3 = 6500;
  double minpt2 = 30;
  // min and max range of L2Res fits. Problem is that this varies
  double vmaxpt2[] =
    {140, 190, 260, 320, 370, 410, 490, 570, 680, 800,
     940, 950, 1020, 1030, 1070, 1120, 1140, 1130,
     1130, 1140, 1120, 1070, 1030, 1020, 950, 940,
     800, 680, 570, 490, 410, 370, 320, 260, 190, 140};
  const int nptmax = sizeof(vmaxpt2)/sizeof(vmaxpt2[0]);
  assert(nptmax==neta);

  TH2D *h2 = new TH2D("h2",";#eta;p_{T} (GeV)",neta,veta,npt,vpt);
  TH2D *h3 = new TH2D("h3",";#eta;p_{T} (GeV)",neta,veta,npt,vpt);
  TH2D *h23 = new TH2D("h23",";#eta;p_{T} (GeV)",neta,veta,npt,vpt);

  TH2D *h2f = new TH2D("h2f",";#eta;p_{T} (GeV)",neta,veta,npt,vpt);
  TH2D *h3f = new TH2D("h3f",";#eta;p_{T} (GeV)",neta,veta,npt,vpt);
  TH2D *h23f = new TH2D("h23f",";#eta;p_{T} (GeV)",neta,veta,npt,vpt);

  TH2D *h2d = new TH2D("h2d",";#eta;p_{T} (GeV)",neta,veta,npt,vpt);
  TH2D *h3d = new TH2D("h3d",";#eta;p_{T} (GeV)",neta,veta,npt,vpt);
  TH2D *h23d = new TH2D("h23d",";#eta;p_{T} (GeV)",neta,veta,npt,vpt);

  TF1 *f2 = new TF1("f2","[2]+[3]*TMath::Log(max([0],min([1],x)))",10,7000);
  f2->SetParameters(30,140,1.32,-0.048);
  TF1 *f3 = new TF1("f3","1./([0]+[1]/x+[2]*log(x)/x+"
		    "[3]*(pow(x/[4],[5])-1)/(pow(x/[4],[5])+1)+"
		    "[6]*pow(x,-0.3051))",10,7000);
  f3->SetParameters(2.3875, 0.5544, -1.55835, -1.32828, 1.6, 0.6532, -0.68140);
  TF1 *f23 = new TF1("f23","[2]+[3]*TMath::Log(max([0],min([1],x)))*"
		     "1./([4]+[5]/x+[6]*log(x)/x+"
		     "[7]*(pow(x/[8],[9])-1)/(pow(x/[8],[9])+1)+"
		     "[10]*pow(x,-0.3051))",10,7000);

  // Fill averaged JECs
  for (int iov = 0; iov != iovs.size(); ++iov) {

    string siov = iovs[iov].first;
    double lumi = iovs[iov].second;
    const char *crun = mrun[siov];
    const char *cver = mver[siov];
    const char *cera = mera[siov];
    string s2 = Form("%s/%s_Run%s_%s_DATA/"
		     "%s_Run%s_%s_DATA_L2Residual_AK4PFchs.txt",
		     "../JECDatabase/textFiles",cera,crun,cver,
		     cera,crun,cver);
    string s23 = Form("%s/%s_Run%s_%s_DATA/"
		      "%s_Run%s_%s_DATA_L2L3Residual_AK4PFchs.txt",
		      "../JECDatabase/textFiles",cera,crun,cver,
		      cera,crun,cver);

    ifstream fin2(s2.c_str());
    string header2;
    getline(fin2, header2);
    if (debug) cout << "Old L2LResidual header:" << endl;
    if (debug) cout << header2 << endl;
    
    ifstream fin23(s23.c_str());
    string header23;
    getline(fin23, header23);
    if (debug) cout << "Old L2L3Residual header:" << endl;
    if (debug) cout << header23 << endl;

    cout << s2 << endl << flush;
    cout << s23 << endl << flush;

    vector<JetCorrectorParameters> v2;
    JetCorrectorParameters *p2 = new JetCorrectorParameters(s2);
    v2.push_back(*p2);
    FactorizedJetCorrector *jec2 = new FactorizedJetCorrector(v2);

    vector<JetCorrectorParameters> v23;
    JetCorrectorParameters *p23 = new JetCorrectorParameters(s23);
    v23.push_back(*p23);
    FactorizedJetCorrector *jec23 = new FactorizedJetCorrector(v23);
    
    for (int ieta = 1; ieta != h2->GetNbinsX()+1; ++ieta) {

      double eta = h2->GetXaxis()->GetBinCenter(ieta);
      
      for (int ipt = 1; ipt != h2->GetNbinsY()+1; ++ipt) {
	double pt  = h2->GetYaxis()->GetBinCenter(ipt);
	
	jec2->setJetPt(pt);
	jec2->setJetEta(eta);
	double c2 = jec2->getCorrection();

	jec23->setJetPt(pt);
	jec23->setJetEta(eta);
	double c23 = jec23->getCorrection();
	double c3 = c23 / c2;
	
	double w = lumi/totlumi;
	h2->SetBinContent(ieta, ipt, h2->GetBinContent(ieta,ipt) + w * c2);
	h3->SetBinContent(ieta, ipt, h3->GetBinContent(ieta,ipt) + w * c3);
	h23->SetBinContent(ieta, ipt, h23->GetBinContent(ieta,ipt) + w * c23);

	double k = 1;
	h2->SetBinError(ieta, ipt, 0.001*k);
	h3->SetBinError(ieta, ipt, 0.001*k);
	h23->SetBinError(ieta, ipt, 0.001*k);
      } // for ipt
    } // for ieta
  } // for iov
    

  string sout2 = Form("textFiles/mergeL2ResTextFiles_%s.txt",run.c_str());
  ofstream fout2(sout2.c_str());

  //string header2 = "{ 1 JetEta 1 JetPt [2]+[3]*TMath::Log(max([0],min([1],x))) Correction L2Relative}";
  string header2 = "{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[8]*0.021*(-1.+1./(1.+exp(-(TMath::Log(x)-5.030)/0.395))))) Correction L2Relative}";
  if (debug) cout << "New L2Residual header:" << endl;
  if (debug) cout << header2 << endl;
  fout2 << header2 << endl;

  string sout23 = Form("textFiles/mergeL2L3ResTextFiles_%s.txt",run.c_str());
  ofstream fout23(sout23.c_str());
  
  string header23 = "{ 1 JetEta 1 JetPt ([2]+[3]*log(max([0],min([1],x))))*1./([4]+[5]/x+[6]*log(x)/x+[7]*(pow(x/[8],[9])-1)/(pow(x/[8],[9])+1)+[10]*pow(x,-0.3051)) Correction L2Relative}";
  if (debug) cout << "New L2L3Residual header:" << endl;
  if (debug) cout << header23 << endl;
  fout23 << header23 << endl;

  // Fit results
  for (int ieta = 1; ieta != h2->GetNbinsX()+1; ++ieta) {

    double eta = h2->GetXaxis()->GetBinCenter(ieta);
    double etamin = veta[ieta-1];
    double etamax = veta[ieta];
    assert(fabs(eta-0.5*(etamin+etamax))<0.5*(etamax-etamin));

    // Fit results with simple parameterization
    double maxpt2 = vmaxpt2[ieta-1];
    f2->FixParameter(0,minpt2);
    f2->FixParameter(1,maxpt2);
    f2->SetRange(minpt2,maxpt2);
    f2->SetParameter(2,1.);
    f2->SetParameter(3,0.);
    TH1D *h2t = h2->ProjectionY("h2t",ieta,ieta);
    h2t->Fit(f2,"QRNM");
    delete h2t;
    f2->SetRange(10,7000);

    // Print L2Residuals
    double q[4]; assert(f2->GetNpar()==4);
    for (int i = 0; i != f2->GetNpar(); ++i) {
      q[i] = f2->GetParameter(i);
    }
    
    int nparnew2 = f2->GetNpar()+2;
    int xmin3 = int(minpt3+0.5);
    int xmax3 = int(maxpt3+0.5);
    int ptmin2 = int(minpt2+0.5);
    int ptmax2 = int(maxpt2+0.5);
    assert(fabs(q[0]-ptmin2)<1);
    assert(fabs(q[1]-ptmax2)<1);

    //fout2 << Form("  %9.6f %9.6f   %d   %d %d   %d   %d   %8.6f %8.6f ",
    fout2 << Form("  %9.3f %9.3f  %d  %d %d  %d %d   %1.0f %1.0f  "
		   "%8.6f %8.6f  %1.0f %1.0f %1.0f",
		  etamin, etamax, nparnew2, xmin3, xmax3, ptmin2, ptmax2,
		  //q[2], q[3]) 
		  1., 1., q[2], q[3], 1., 0., 0.)
	   << endl;


    f3->SetRange(minpt3,maxpt3);
    TH1D *h3t = h3->ProjectionY("h3t",ieta,ieta);
    h3t->Fit(f3,"QRNM");
    delete h3t;
    f3->SetRange(10,7000);

    // Print L3Residuals
    f23->SetParameters(f2->GetParameter(0),f2->GetParameter(1),
		       f2->GetParameter(2),f2->GetParameter(3),
		       f3->GetParameter(0),f3->GetParameter(1),
		       f3->GetParameter(2),f3->GetParameter(3),
		       f3->GetParameter(4),f3->GetParameter(5),
		       f3->GetParameter(6));
    double p[11]; assert(f23->GetNpar()==11);
    for (int i = 0; i != f23->GetNpar(); ++i) {
      p[i] = f23->GetParameter(i);
    }
    
    int nparnew23 = f23->GetNpar()+2;
    //int xmin3 = int(minpt3+0.5);
    //int xmax3 = int(maxpt3+0.5);
    //int ptmin2 = int(minpt2+0.5);
    //int ptmax2 = int(maxpt2+0.5);
    assert(fabs(p[0]-ptmin2)<1);
    assert(fabs(p[1]-ptmax2)<1);

    fout23 << Form("  %9.6f %9.6f   %d   %d %d   %d   %d   %8.6f %8.6f "
		   "%5.4f %5.4f %5.5f %5.5f %5.1f %5.4f %5.5f",
		   etamin, etamax, nparnew23, xmin3, xmax3, ptmin2, ptmax2,
		   p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10])
	   << endl;

    for (int ipt = 1; ipt != h2->GetNbinsY()+1; ++ipt) {
      double pt  = h2->GetYaxis()->GetBinCenter(ipt);
      
      double c2 = f2->Eval(pt);
      double c3 = f3->Eval(pt);
      double c23 = c2 * c3;
      
      h2f->SetBinContent(ieta, ipt, c2);
      h3f->SetBinContent(ieta, ipt, c3);
      h23f->SetBinContent(ieta, ipt, c23);
    } // for ipt
  } // for ieta
  
  
  h2d->Add(h2,h2f,1,-1);
  h3d->Add(h3,h3f,1,-1);
  h23d->Add(h23,h23f,1,-1);


  gStyle->SetOptStat(0);

  //double dz = 0.005;//0.002;

  TCanvas *c23 = new TCanvas("c23","c23",600,600);
  gPad->SetLogy();
  gPad->SetRightMargin(0.15);
  //h23->Draw("COLZ");
  //h23f->Draw("COLZ");
  h23d->Draw("COLZ");
  h23->GetZaxis()->SetRangeUser(0.93,1.23);
  h23f->GetZaxis()->SetRangeUser(0.93,1.23);
  h23d->GetZaxis()->SetRangeUser(-dz,+dz);//-0.001,0.001);

  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  gPad->SetLogy();
  gPad->SetRightMargin(0.15);
  //h2->Draw("COLZ");
  //h2f->Draw("COLZ");
  h2d->Draw("COLZ");
  h2->GetZaxis()->SetRangeUser(0.93,1.23);
  h2f->GetZaxis()->SetRangeUser(0.93,1.23);
  h2d->GetZaxis()->SetRangeUser(-dz,+dz);//-0.001,0.001);

  TCanvas *c3 = new TCanvas("c3","c3",600,600);
  gPad->SetLogy();
  gPad->SetRightMargin(0.15);
  //h3->Draw("COLZ");
  //h3f->Draw("COLZ");
  h3d->Draw("COLZ");
  //h3->GetZaxis()->SetRangeUser(0.93,1.23);
  h3->GetZaxis()->SetRangeUser(0.95,1.05);
  h3f->GetZaxis()->SetRangeUser(0.95,1.05);
  h3d->GetZaxis()->SetRangeUser(-dz,+dz);//-0.001,0.001);

} // mergeL2L3ResTextFiles
