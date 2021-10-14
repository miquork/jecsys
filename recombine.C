// Purpose: recombine IOVs or epochs for JEC global fit
// Nominally setup to combine 2016BCDEF, 2016GH, 2017BCDEF and 2018ABCD
// into Run 2 full data set
// Take jecdataX.root as input and create new jecdataY.root
//
// NB: need to be careful about combining Zee and Zmm, applying corrections,
//     dealing with missing (sniped) data points
// orig/ directory is safest for consistent data sets, but
// misses all corrections, combinations, posterior ratios etc.
//
// reprocess.C is currently not setup to reused jecdataX.root as input,
// although it could make sense to implement this option
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TKey.h"

#include <string>
#include <vector>
#include <iostream>

using namespace std;

const bool _debug = true;
const bool skipEmpty = true; // leave bin empty if any input empty
const bool skipEmptyEnd = false; // leave bin empty if any input end bin missing
const bool patchHist = true; // copy overlapping contents if bin # mismatch
const bool skipDifferent = true; // leav bin empty if any bins mismatch

void recurseJECDataFile(std::vector<pair<double, TDirectory*> > &indirs,
			TDirectory *outdir, int lvl = 0);

void recombine(string type="Run2") {

  vector<pair<double,string> > finns;
  //finns.push_back(make_pair<double,string>(19.7,"2016BCDEF"));
  //finns.push_back(make_pair<double,string>(16.8,"2016GH"));
  //finns.push_back(make_pair<double,string>(41.48,"2017BCDEF"));
  //finns.push_back(make_pair<double,string>(59.82,"2018ABCD"));
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiLUM
  finns.push_back(make_pair<double,string>(36.33*19.7/36.5,"2016BCDEF"));
  finns.push_back(make_pair<double,string>(36.33*16.8/36.5,"2016GH"));
  finns.push_back(make_pair<double,string>(41.48,"2017BCDEF"));
  finns.push_back(make_pair<double,string>(59.83,"2018ABCD"));

  vector<pair<double,TDirectory*> > fins(finns.size());
  for (unsigned int i = 0; i != finns.size(); ++i) {
    fins[i].first = finns[i].first;
    string s = finns[i].second;
    TFile *f = new TFile(Form("rootfiles/jecdata%s.root",s.c_str()),"READ");
    assert(f && !f->IsZombie());
    TDirectory *d = f->GetDirectory("");
    assert(d);
    fins[i].second = d;
  } // for i in finns

  TFile *fout = new TFile("rootfiles/jecdataRun2TestData.root","RECREATE");

  cout << "Calling recombine("<<type<<");" << endl;
  cout << "Input files ";
  for (unsigned int i = 0; i != fins.size(); ++i) {
    cout << (i==0 ? "":", ") << fins[i].second->GetName();
  }
  cout << endl;
  cout << "Output file " << fout->GetName() << endl;
  cout << "Starting recursive loop. This may take a minute" << endl << flush;

  // Loop over all the directories recursively                                 
  recurseJECDataFile(fins, fout);

  cout << endl;
  cout << "Recursive loop done." << endl;
  cout << "Output written in " << fout->GetName() << endl;
  fout->Close();
  cout << "Output file closed" << endl;
  fout->Delete();
  cout << "Output file pointer deleted" << endl << flush;

  for (unsigned int i = 0; i != fins.size(); ++i) {
    fins[i].second->Close();
    fins[i].second->Delete();
  }
} // recombine

void recurseJECDataFile(std::vector<pair<double, TDirectory*> > &indirs,
			TDirectory *outdir, int lvl) {
  
  // Automatically go through the list of keys (directories)                   
  assert(indirs.size()!=0);
  TDirectory *indir = indirs[0].second;
  TList *keys = indir->GetListOfKeys();
  TListIter itkey(keys);
  TObject *obj;
  TKey *key;

  while ( (key = dynamic_cast<TKey*>(itkey.Next())) ) {
    const char *c = Form("%%%ds%%s/%%s",lvl+1);
    if (_debug) cout << Form(c," ",indir->GetName(),key->GetName())
		     << endl << flush;
    obj = key->ReadObj(); assert(obj);

    // Found a subdirectory: copy it to output and go deeper                   
    if (obj->InheritsFrom("TDirectory")) {
      int loclvl = lvl;
      outdir->mkdir(obj->GetName());
      bool enteroutdir = outdir->cd(obj->GetName());
      assert(enteroutdir);
      TDirectory *outdir2 = outdir->GetDirectory(obj->GetName());
      assert(outdir2);
      outdir2->cd();

      vector<pair<double, TDirectory*> > indirs2(indirs.size());
      for (unsigned int i = 0; i != indirs.size(); ++i) {
	TDirectory *indir = indirs[i].second;
	bool enterindir = indir->cd(obj->GetName());
	assert(enterindir);
	TDirectory *indir2 = indir->GetDirectory(obj->GetName());
	indir2->cd();
	indirs2[i].first = indirs[i].first;
	indirs2[i].second = indir2;
      }

      if (loclvl>=0) loclvl++;

      //if (loclvl==1) 
      cout << endl << "Entering: " << indirs2[0].second->GetName() << endl;

      recurseJECDataFile(indirs2, outdir2, loclvl);
    } // inherits from TDirectory  
    else if (obj->InheritsFrom("TH1")) { // Combine histograms

      outdir->cd();
      //string name = obj->GetName();
      //string trgname = indir->GetName();

      TObject *obj2 = obj->Clone(key->GetName()); // Copy input histo to output
      TH1 *hout(0);
      if (obj->InheritsFrom("TH1D")) hout = dynamic_cast<TH1D*>(obj2);
      if (obj->InheritsFrom("TH1F")) hout = dynamic_cast<TH1F*>(obj2);
      //if (obj->InheritsFrom("TH2D")) hout = dynamic_cast<TH2D*>(obj2);
      if (!hout) { // TH2D, TH3D etc.
	obj2->Delete();
	continue;
      }

      // Start adding histograms from files
      hout->Reset();
      double sumw(0);
      for (unsigned int i = 0; i != indirs.size(); ++i) {

	double w1 = sumw;
	double w2 = indirs[i].first;
	sumw = (w1 + w2);

	TH1D *hin = (TH1D*)indirs[i].second->Get(key->GetName());
	//assert(hin);
	if (!hin) {
	  cout << outdir->GetName()<<"/"<<key->GetName()
	       << ": missing from input #"<<(i+1)<<", "
	       << (skipEmpty ? "skipping" : "beware") << endl;
	  if (skipEmpty) {
	    hout->Reset();
	    break;
	  }
	  else
	    continue;
	}
	
	//assert(hin->GetNbinsX()==hout->GetNbinsX()+1);
	if (hin->GetNbinsX()!=hout->GetNbinsX()) {//+1) {
	  cout << outdir->GetName()<<"/"<<key->GetName()
	       << ": bin # mismatch ("
	       << hin->GetNbinsX() << " vs " << hout->GetNbinsX() <<"), " 
	       << (skipEmptyEnd ? "skip end" : 
		   (patchHist ? "patch histo" : "beware")) << endl;
	  //hout->Reset();
	  //break;
	  if (patchHist) {
	    TH1D *htmp = hin;
	    hin = (TH1D*)hout->Clone(Form("%s_patch",htmp->GetName()));
	    hin->Reset();
	    for (int i = 1; i != hin->GetNbinsX()+1; ++i) {
	      int j = htmp->FindBin(hin->GetBinCenter(i));
	      assert(htmp->GetBinLowEdge(j)==hin->GetBinLowEdge(i));
	      assert(htmp->GetBinLowEdge(j+1)==hin->GetBinLowEdge(i+1));
	      hin->SetBinContent(i, htmp->GetBinContent(j));
	      hin->SetBinError(i, htmp->GetBinError(j));
	    } // for i
	  } // patchHist
	} // bin # mistmatch
	

	for (int j = 1; j != min(hout->GetNbinsX(),hin->GetNbinsX())+1; ++j) {

	  if (obj->InheritsFrom("TH1D") || obj->InheritsFrom("TH1F")) {
	    //assert(hin->GetBinLowEdge(j)==hout->GetBinLowEdge(j));
	    if (hin->GetBinLowEdge(j)!=hout->GetBinLowEdge(j)) {
	      cout << outdir->GetName()<<"/"<<key->GetName()
		   << Form(": bin mismatch ([%1.0f,%1.0f] vs [%1.0f,%1.0f]), ",
			   hin->GetBinLowEdge(j), hin->GetBinLowEdge(j+1),
			   hout->GetBinLowEdge(j), hout->GetBinLowEdge(j+1))
		   << (skipDifferent ? "skipping" : "beware") << endl;
	      if (skipDifferent) {
		hout->SetBinContent(j, 0);
		hout->SetBinError(j, 0);
	      }
	      continue;
	    }
	  }
	  
	  // Patch 'sniped' bins
	  if (hin->GetBinContent(j)==0 && hin->GetBinError(j)==0) {
	    hout->SetBinContent(j, skipEmpty ? 0 : hout->GetBinContent(j));
	    hout->SetBinError(j, skipEmpty ? 0 : hout->GetBinError(j));
	  }
	  else if (hout->GetBinContent(j)==0 && hout->GetBinError(j)==0) {
	    hout->SetBinContent(j, skipEmpty && i ? 0 : hin->GetBinContent(j));
	    hout->SetBinError(j, skipEmpty && i ? 0 : hin->GetBinError(j));
	  }
	  else {
	    hout->SetBinContent(j, (w1 * hout->GetBinContent(j) +
				    w2 * hin->GetBinContent(j)) / sumw);
	    hout->SetBinError(j, sqrt(pow(w1 * hout->GetBinError(j),2) +
				      pow(w2 * hin->GetBinError(j),2)) / sumw);
	  }
	} // for j in hout
	for (int j = min(hout->GetNbinsX(),hin->GetNbinsX())+1;
	     j != hout->GetNbinsX()+1; ++j) {
	  if (skipEmptyEnd) {
	    hout->SetBinContent(j, 0);
	    hout->SetBinError(j, 0);
	  }
	}
      } // for i in indirs
      
      // Save the stuff into an identical directory                            
      outdir->cd();
      obj2->Write();
      obj2->Delete();
      indir->cd();
    } // inherits from TH1
    else if (obj->InheritsFrom("TGraphErrors")) { // Combine graphs

      outdir->cd();
      //string name = obj->GetName();
      //string trgname = indir->GetName();

      TObject *obj2 = obj->Clone(key->GetName()); // Copy input graph to output
      TGraphErrors *gout(0), *gin(0);
      if (obj->InheritsFrom("TGraphErrors"))
	gout = dynamic_cast<TGraphErrors*>(obj2);
      if (!gout) {
	continue;
      }

      // Start adding histograms from files
      double sumw(indirs[0].first);
      for (unsigned int i = 1; i != indirs.size(); ++i) {

	double w1 = sumw;
	double w2 = indirs[i].first;
	sumw = (w1 + w2);

	gin = (TGraphErrors*)indirs[i].second->Get(key->GetName());
	if (!gin) {
	  cout << outdir->GetName()<<"/"<<key->GetName()
	       << ": missing from input #"<<(i+1)<<", "
	       << (skipEmpty ? "skipping" : "beware") << endl;
	  if (skipEmpty) {
	    //hout->Reset();
	    break;
	  }
	  else
	    continue;
	} // !gin

	// For every point in gout, find nearest point in gin
	// Then check that gout is also the nearest point to gin
	// If yes, remove point from bin to see that all points get matched
	gin = (TGraphErrors*)gin->Clone("tmpgin");
	for (int j = 0; j != gout->GetN(); ++j) {
	  double drmin(-1); int kk(-1);
	  for (int k = 0; k != gin->GetN(); ++k) {
	    double dr = fabs(gout->GetX()[j]-gin->GetX()[k]);
	    if (drmin<0 || dr<drmin) { drmin = dr; kk = k; }
	  } // for k in gin
	  double drmin2(-1); int jj(-1);
	  for (int k = 0; k != gout->GetN() && kk != -1; ++k) {
	    double dr = fabs(gout->GetX()[k]-gin->GetX()[kk]);
	    if (drmin2<0 || dr<drmin2) { drmin2 = dr; jj = k; }
	  } // for k in gout
	  // Nearest match both ways
	  if (jj == j) {
	   
	    if (jj!=-1 && kk != -1) {
	      // Merge points
	      double y = (w1*gout->GetY()[jj] + w2*gin->GetY()[kk]) / sumw;
	      double x = (w1*gout->GetX()[jj] + w2*gin->GetX()[kk]) / sumw;
	      double ey = sqrt(pow(w1*gout->GetEY()[jj],2) +
			       pow(w2*gin->GetEY()[kk],2)) / sumw;
	      double ex = sqrt(pow(w1*gout->GetEX()[jj],2) +
			       pow(w2*gin->GetEX()[kk],2)) / sumw;

	      gout->SetPoint(jj, x, y);
	      gout->SetPointError(jj, ex, ey);
	      
	      // To keep track of unmatched points, remove matches
	      gin->RemovePoint(kk);
	    }
	    else { // no points in gin
	      // Keep existing point in gout
	    }
	  }
	  else {
	    cout << outdir->GetName()<<"/"<<key->GetName()
		 << " point "<<j<<" pt "<<gout->GetX()[j]
		 << " near pt "<< (kk!=-1 ? gin->GetX()[kk] : 0)
		 << " missing match in input #"<<(i+1) << endl;
	  }
	} // for j in gout
	if (gin->GetN()!=0) {
	  cout << outdir->GetName()<<"/"<<key->GetName()
	       << " "<<gin->GetN()<<" points from pt "<<gin->GetX()[0]
	       << " missing match in input #"<<(i+1) << endl;
	}

	// Reset pointer for next rounds
	gin->Delete(); gin = 0;
      } // for i in indirs

      // Save the stuff into an identical directory                            
      outdir->cd();
      obj2->Write();
      obj2->Delete();
      indir->cd();
    }
    else if (false) { // Just skip others

      // Save the stuff into an identical directory                            
      TObject *obj2 = obj->Clone(key->GetName()); // Copy the input to output      outdir->cd();
      obj2->Write();
      obj2->Delete();
      indir->cd();
    }
    obj->Delete();
  } // while key      

} // recurseJECDataFile
