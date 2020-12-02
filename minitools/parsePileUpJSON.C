// Parse pileup_latest.txt and return true pileup for (run,LS) pair
// Input from /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/pileup_latest.txt
// Re-formatted with minitools/JSONtoASCII.py
// Run with 'root -l -b -q minitools/parsePileupJSON.C+g'
#ifndef __parsePileUpJSON_C__
#define __parsePileUpJSON_C__

#include "TH1D.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TChain.h"
#include "TFile.h"

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <utility> // std::pair

using namespace std;

bool debug = false;//true;

struct lumiInfo {
  double lum;
  double muavg;
  double murms;
  lumiInfo(double lum_=0, double muavg_=0, double murms_=0)
    : lum(lum_), muavg(muavg_), murms(murms_) {};
};

map<int, map<int, lumiInfo> > _mus;
map<int, map<int, lumiInfo> > _mus2;

// Read pileup_latest.txt into memory
// NB: first reformat pileup_latext.txt into C++ -friedly format:
// minitools/JSONtoASCII.py rootfiles/pileup_latest_2018.txt > pileup_ASCII_2018.txt
void parsePileUpJSON(string filename="rootfiles/pileup_ASCII_2018.txt", 
		     string filename2="", 
		     string json="",
		     double minbXsec = 69200) {

  TDirectory *curdir = gDirectory;

  cout << Form("parsePileUpJSON(minbXsec=%1.3g):",minbXsec) << endl;
  cout << "Opening " << filename << "...";
  ifstream f(filename.c_str(),ios::in);
  if (f.is_open()) cout << "ok" << endl;
  else { cout << "failure" << endl; return; }

  if (filename2!="") cout << "Opening " << filename2 << "...";
  ifstream f2(filename2.c_str(),ios::in);
  if (f2.is_open()) cout << "ok" << endl;
  else if (filename2!="") { cout << "failure" << endl; return; }

  if (json!="") cout << "Opening " << json << "...";
  ifstream fj(json.c_str(),ios::in);
  if (fj.is_open()) cout << "ok" << endl;
  else if (json!="") { cout << "failure" << endl; return; }

  TH1D *hpu = new TH1D("hpu",";true PU;",600,0,60);
  TH1D *hpux = new TH1D("hpux",";true PU (cleaned);",600,0,60);
  TH1D *hpus = new TH1D("hpus",";true PU (smeared);",600,0,60);
  TH1D *hpuo = new TH1D("hpuo",";observed PU;",60,0,60);

  TH1D *hruns = new TH1D("hruns",";run;",20000,190000,210000);
  TH1D *hrunsx = new TH1D("hrunsx",";run;",20000,190000,210000);

  TH2D *h2pu = new TH2D("h2pu",";true PU (px); true PU (hf);",
			600,0,60, 600,0,60);

  set<int> hiruns;
  hiruns.insert(199703); // part
  hiruns.insert(200190); // part
  hiruns.insert(203830);
  hiruns.insert(203832);
  hiruns.insert(203833);
  hiruns.insert(203834);
  hiruns.insert(203835); // part
  hiruns.insert(208509); // part

  // Load list of good LS from re-reco JSON
  set<pair<int, int> > goodls;
  if (fj.is_open()) {
    int run, ls1, ls2;
    while (fj >> run >> ls1 >> ls2) {

      for (int i = ls1; i != ls2+1; ++i) {
	//goodls.insert(make_pair<int,int>(run, ls1+i));
	//goodls.insert(make_pair<int,int>(run, i));
	goodls.insert(pair<int,int>(run, i));
      }
    }
  }

  // Read information from file
  double maxmu = 60;
  int nrun(0), nls(0);

  if (f.is_open()) {
    int run, ls, cnt(0), cntMax(-1);
    double lum, xsavg, xsrms;
    while ( f >> run >> ls >> lum >> xsrms >> xsavg && ++cnt!=cntMax) {
      
      if (debug && cnt<10)
	cout << run<<" "<<ls<<" "<<lum<<" "<<xsavg<<" "<<xsrms << endl;
      
      //if (json=="" || goodls.find(make_pair<int,int>(run,ls))!=goodls.end()) {
      if (json=="" || goodls.find(pair<int,int>(run,ls))!=goodls.end()) {
	
	double mu = xsavg * minbXsec;
	double murms = xsrms * minbXsec;
	//if (murms>mu) mu = min(mu,maxmu);
	//murms = max(0.,min(mu,murms));
	
	if (_mus.find(run)==_mus.end()) ++nrun;
	if (_mus[run].find(ls)==_mus[run].end()) ++nls;

	_mus[run][ls] = lumiInfo(lum, mu, murms);
	
	lumiInfo const& m = _mus[run][ls];
	if (debug && cnt<10)
	cout << "=>"<<m.muavg<<" +/- "<<m.murms<<endl;
	
	if (debug) {
	  hruns->Fill(run, lum);
	  if (mu>maxmu) hrunsx->Fill(run, lum);
	  
	  hpu->Fill(mu, lum);
	  if (hiruns.find(run)!=hiruns.end()) hpux->Fill(run, lum);
	  double murnd = gRandom->Gaus(mu, murms); 
	  hpus->Fill(murnd, lum);
	  int muobs = gRandom->Poisson(murnd);
	  hpuo->Fill(muobs, lum);
	}
      }
    }
  } // read from file

  // Second opinion on the luminosity using another luminosity monitor
  if (f2.is_open()) {

    int run, ls, cnt(0), cntMax(-1);
    double lum, xsavg, xsrms;
    while ( f2 >> run >> ls >> lum >> xsrms >> xsavg && ++cnt!=cntMax) {

      if (debug && cnt<10)
	cout << run<<" "<<ls<<" "<<lum<<" "<<xsavg<<" "<<xsrms << endl;
      
      //if (json=="" || goodls.find(make_pair<int,int>(run,ls))!=goodls.end()) {
      if (json=="" || goodls.find(pair<int,int>(run,ls))!=goodls.end()) {

	double mu = xsavg * minbXsec;
	double murms = xsrms * minbXsec;
	
	_mus2[run][ls] = lumiInfo(lum, mu, murms);
	
	lumiInfo const& m = _mus2[run][ls];
	if (debug && cnt<10)
	  cout << "=>"<<m.muavg<<" +/- "<<m.murms<<endl;
	
	double mupx = _mus[run][ls].muavg;
	h2pu->Fill(mupx, mu, lum);      
      }
    }
  }

  // Store lumi information in ntuple for easy access
  if (f2.is_open()) {

    TFile *lumfile = new TFile("lumtree.root","RECREATE");
    TTree *lumtree = new TTree("t","t");
    Int_t t_run, t_ls;
    Float_t t_lpx, t_lhf, t_mpx, t_mhf, t_spx, t_shf;
    lumtree->Branch("run",&t_run,"run/I");
    lumtree->Branch("ls",&t_ls,"ls/I");
    lumtree->Branch("lpx",&t_lpx,"lpx/F");
    lumtree->Branch("lhf",&t_lhf,"lhf/F");
    lumtree->Branch("mpx",&t_mpx,"mpx/F");
    lumtree->Branch("mhf",&t_mhf,"mhf/F");
    lumtree->Branch("spx",&t_spx,"spx/F");
    lumtree->Branch("shf",&t_shf,"shf/F");

    typedef map<int, map<int, lumiInfo> >::const_iterator IT;
    typedef map<int, lumiInfo>::const_iterator JT;
    for (IT it = _mus.begin(); it != _mus.end(); ++it) {
      for (JT jt = it->second.begin(); jt != it->second.end(); ++jt) {
	
	int run = it->first;
	int ls = jt->first;

	t_run = run;
	t_ls = ls;
	t_lpx = _mus[run][ls].lum;
	t_mpx = _mus[run][ls].muavg;
	t_spx = _mus[run][ls].murms;
	
	if (_mus2.find(run)==_mus2.end()) {
	  cout << "Run "<<run<<" of "<<_mus[run].size()<<" LS missing from "
	       << filename2 << endl;
	}
	if (_mus2[run].find(ls)==_mus2[run].end()) {
	  t_lhf = t_mhf = t_shf = -1;
	}
	else {
	  t_lhf = _mus2[run][ls].lum;
	  t_mhf = _mus2[run][ls].muavg;
	  t_shf = _mus2[run][ls].murms;
	}
	
	lumtree->Fill();
      } // for jt
    } // for it
    
    // Check that all the LS are also in the first file
    for (IT it = _mus2.begin(); it != _mus2.end(); ++it) {
      for (JT jt = it->second.begin(); jt != it->second.end(); ++jt) {
	
	int run = it->first;
	int ls = jt->first;
	
	assert(_mus.find(run)!=_mus.end());
	assert(_mus[run].find(ls)!=_mus[run].end());
      } // for jt
    } // for it

    if (debug) lumtree->Print();
    lumtree->AutoSave();
    cout << "Saved pileup data to " << lumfile->GetName() << endl;
    delete lumfile;
    curdir->cd();

    TCanvas *c3 = new TCanvas("c3","c3",600,600);
    h2pu->Draw("COLZ");
    c3->SaveAs("pdf/parsePileUpJSON_HFvsPX.pdf");
  } // second file open

  if (debug) {
    hpu->SetLineColor(kBlue);
    hpu->Draw();
    hpux->SetLineColor(kBlue+2);
    hpux->SetLineStyle(kDotted);
    hpux->Draw("SAME");
    hpus->SetLineColor(kBlack);
    hpus->Draw("SAME");
    hpuo->Scale(0.1); // for bin size
    hpuo->SetLineColor(kRed);
    hpuo->Draw("SAME");

    TCanvas *c2 = new TCanvas("c2","c2",600,600); c2->cd();
    hruns->Draw();
    hrunsx->SetLineColor(kRed);
    hrunsx->Draw("SAME");
  }

  cout << "Read in "<<nls<<" lumisections in "<<nrun<<" runs"<<endl;
} // parsePileUpJSON

double getTruePU(int run, int ls) {

  if (_mus.find(run)==_mus.end()) return 0;
  if (_mus[run].find(ls)==_mus[run].end()) return 0;
  return _mus[run][ls].muavg;
}

#endif __parsePileUpJSON_C__
