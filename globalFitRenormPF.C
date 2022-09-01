// Purpose: Renormalize Data and MC CHF, NEF and NHF to add up to 1.000
//          Update Data-MC difference accordingly
//          Use '_a100' as input and store result as '_ren'
//          Explicit loop for easy readability
//          Bonus: pre-combine PF inputs for better visualization
//          Bonus: also pre-combine response inputs
#include "TFile.h"
#include "TGraphErrors.h"

#include "tools.C"

#include <map>
#include <string>

using namespace std;

void globalFitPrecombinePF(string s);
void globalFitPrecombineHDM(string s, bool addMJ, bool scaleRes);

void globalFitRenormPF(string s = "Run2Test") {

  TDirectory *curdir = gDirectory;

  const char *c = s.c_str();
  TFile *f = new TFile(Form("rootfiles/jecdata%s.root",c),"UPDATE");
  assert(f && !f->IsZombie());

  curdir->cd();

  // List all datasets providing PF composition inputs
  vector<string> sets;
  if (s=="2017H") {
    sets.push_back("zjet");
    sets.push_back("pfjet");
  }
  else {
    sets.push_back("gamjet");
    sets.push_back("zeejet");
    sets.push_back("zmmjet");
    sets.push_back("zlljet");
    sets.push_back("zjet");
    sets.push_back("pfjet");
  }

  // List inputs to be normalized
  vector<string> vars;
  vars.push_back("chf");
  vars.push_back("nef");
  vars.push_back("nhf");

  // Directories to loop for reading in data
  vector<string> dirs;
  dirs.push_back("data");
  dirs.push_back("mc");
  dirs.push_back("ratio");

  // Read in data and calculate sums
  TGraphErrors *g(0), *gs(0);
  map<string, map<string, map<string, TGraphErrors*> > > m;
  for (unsigned int id = 0; id != dirs.size(); ++id) {
    for (unsigned int iv = 0; iv != vars.size(); ++iv) {
      for (unsigned int is = 0; is != sets.size(); ++is) {
    
	cout << "." << flush;

	// Index to string for maps
	const char *cd = dirs[id].c_str();
	const char *cv = vars[iv].c_str();
	const char *cs = sets[is].c_str();

	// Read in data
	g = (TGraphErrors*)f->Get(Form("%s/eta00-13/%s_%s_a100",cd,cv,cs));
	if (!g) // for pfjet_a30
	  g = (TGraphErrors*)f->Get(Form("%s/eta00-13/%s_%s_a30",cd,cv,cs));
	assert(g);
	m[cs][cv][cd] = g;

	// Calculate sums
	gs = m[cs]["sum"][cd];
	if (!gs) gs = (TGraphErrors*)g->Clone(Form("%s_sum",g->GetName()));
	else gs = tools::diffGraphs(g, gs, 1, -1);
	m[cs]["sum"][cd] = gs;
      } // for is in sets
    } // for iv in vars
  } // for id in dirs

  // Renormalize data and MC, recalculate ratio
  TGraphErrors *gr(0), *gd(0), *gm(0);
  for (unsigned int id = 0; id != dirs.size(); ++id) {
    for (unsigned int is = 0; is != sets.size(); ++is) {
      for (unsigned int iv = 0; iv != vars.size(); ++iv) {

	// Index to string for maps
	string sd = dirs[id];
	string ss = sets[is];
	string sv = vars[iv];

	// Renormalize data
	if (sd=="data" || sd=="mc") {
	  g = m[ss][sv][sd];      assert(g);
	  gs = m[ss]["sum"][sd];  assert(gs);
	  gr = tools::ratioGraphs(g, gs);
	}
	// Recalculate Data-MC difference
	else if (sd=="ratio") {
	  gd = m[ss][sv+"_ren"]["data"]; assert(gd);
	  gm = m[ss][sv+"_ren"]["mc"];   assert(gm);
	  gr = tools::diffGraphs(gd, gm);
	}
	assert(gr);
	m[ss][sv+"_ren"][sd] = gr;
	  
      } // for is in sets
    } // for iv in vars
  } // for id in dirs

  // Write out results
  for (unsigned int id = 0; id != dirs.size(); ++id) {
    for (unsigned int iv = 0; iv != vars.size(); ++iv) {
      for (unsigned int is = 0; is != sets.size(); ++is) {
    
	// Index to string for maps
	string sv = vars[iv];
	const char *cd = dirs[id].c_str();
	const char *cv = vars[iv].c_str();
	const char *cs = sets[is].c_str();

	g = m[cs][sv+"_ren"][cd]; assert(g);
	f->cd(Form("%s/eta00-13",cd));
	g->SetName(Form("%s_%s_ren",cv,cs));
	g->Write(Form("%s_%s_ren",cv,cs),TObject::kOverwrite);
      } // for is in sets
    } // for iv in vars
  } // for id in dirs

  f->Write();
  f->Close();
  curdir->cd();

  globalFitPrecombinePF(s);
  globalFitPrecombineHDM(s,false,false);
  if (s!="2017H") {
    globalFitPrecombineHDM(s,true,false); // addMJ
    globalFitPrecombineHDM(s,true,true); // addMJ+scaleRes
  }
  else {
    globalFitPrecombineHDM(s,false,true); // scaleRes
  }
} // globalFitRenormPF.C


// Purpose: pre-combine PF compositions for easier plotting / fitting
void globalFitPrecombinePF(string s) {

  TDirectory *curdir = gDirectory;

  const char *c = s.c_str();
  TFile *f = new TFile(Form("rootfiles/jecdata%s.root",c),"UPDATE");
  assert(f && !f->IsZombie());

  curdir->cd();

  // List all datasets providing PF composition inputs
  vector<string> sets;
  if (s=="2017H") {
    sets.push_back("zjet");
    sets.push_back("pfjet");
  }
  else {
    sets.push_back("gamjet");
    sets.push_back("zlljet");
    sets.push_back("zjet");
    sets.push_back("pfjet");
  }

  // List inputs to be pre-combined
  vector<string> vars;
  vars.push_back("chf");
  vars.push_back("nef");
  vars.push_back("nhf");

  // Directories to loop for reading in data
  vector<string> dirs;
  dirs.push_back("data");
  dirs.push_back("mc");
  dirs.push_back("ratio");

  // zlljet: 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 85, 105, 130, 175, 230, 300, 400, 500, 700, 1000, 1500
  // gamjet: 15, 20, 25, 30, 35, 40, 50, 60, 70, 85, 105, 130, 175, 230, 300, 400, 500, 600, 700, 850, 1000, 1200, 1450, 1750,
  // multijet: 10, 15, 21, 28, 37, 49, 64, 84, 114, 153, 196, 245, 300, 362, 430, 507, 592, 686, 790, 905, 1032, 1172, 1327, 1497, 1684, 1890, 2116, 2366, 2640, 2941, 3273, 5220,
  // pfjet?: 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890,
  const double x[] =
    {15, 20, 25, 30, 35, 40, 45, 50, 60, // Z/W only
     70, 85, 130, // +gamjet
     175, 230, 300, 400, 500, // still Z/gamma domain
     592, 686, 790, 905, 1032, 1172, 1327, 1497, // multijet
     1684, 1890, 2116, 2366, 2640, 2941, 3273, 5220}; // multijet
  const int nx = sizeof(x)/sizeof(x[0])-1;
  TH1D *h = new TH1D(Form("h_%s",c),"",nx,x);
  TF1 *f1 = new TF1(Form("f1_%s",c),"[0]",15,7000);

  // Read in data and calculate sums
  TGraphErrors *g(0);
  map<string, map<string, map<string, TGraphErrors*> > > m;
  TH1D *hs(0), *hs2(0), *hsw(0); 
  map<string, map<string, map<string, TH1D*> > > ms;
  for (unsigned int id = 0; id != dirs.size(); ++id) {
    for (unsigned int iv = 0; iv != vars.size(); ++iv) {
      for (unsigned int is = 0; is != sets.size(); ++is) {
    
	cout << "." << flush;

	// Index to string for maps
	const char *cd = dirs[id].c_str();
	const char *cv = vars[iv].c_str();
	const char *cs = sets[is].c_str(); string ss = cs;

	// Read in data
	g = (TGraphErrors*)f->Get(Form("%s/eta00-13/%s_%s_ren",cd,cv,cs));
	assert(g);
	m[cs][cv][cd] = g;

	// Calculate sums
	hs  = ms["sum"][cv][cd];
	hs2 = ms["sum2"][cv][cd];
	hsw = ms["sumw"][cv][cd];
	if (!hs)  hs  = (TH1D*)h->Clone(Form("%s_sum",g->GetName()));
	if (!hs2) hs2 = (TH1D*)h->Clone(Form("%s_sum2",g->GetName()));
	if (!hsw) hsw = (TH1D*)h->Clone(Form("%s_sumw",g->GetName()));

	hs->SetLineColor(g->GetLineColor());
	hs2->SetLineColor(g->GetLineColor());
	hsw->SetLineColor(g->GetLineColor());

	hs->SetMarkerColor(g->GetMarkerColor());
	hs2->SetMarkerColor(g->GetMarkerColor());
	hsw->SetMarkerColor(g->GetMarkerColor());

	for (int i = 1; i != hs->GetNbinsX()+1; ++i) {
	  double ptmin = hs->GetBinLowEdge(i);
	  double ptmax = hs->GetBinLowEdge(i+1);
	  f1->SetRange(ptmin,ptmax);
	  f1->SetParameter(0,0); // in case no points to fit
	  f1->SetParError(0,0);  // in case no points to fit
	  g->Fit(f1,"QRN");
	  if (f1->GetParError(0)>0 && // check something was fit
	      //!(ptmin<40) &&
	      !(ptmin<30 && ss=="pfjet") && // exclude low pT PF jets
	      !(ptmin<70 && ss=="gamjet")) { // exclude low pT photon+jet
	    double w0 = hsw->GetBinContent(i);
	    double yw0 = hs->GetBinContent(i) * w0;
	    double y2w0 = pow(hs->GetBinContent(i),2) * w0;
	    double eyw0 = hs->GetBinError(i) * w0;
	    double w = 1./pow(f1->GetParError(0),2);
	    double yw = f1->GetParameter(0) * w;
	    double y2w = pow(f1->GetParameter(0),2) * w;
	    double eyw = f1->GetParError(0) * w;
	    hs->SetBinContent(i, (yw0 + yw) / (w0 + w));
	    hs->SetBinError(i, (eyw0 + eyw) / (w0 + w));
	    hs2->SetBinContent(i, (y2w0 + y2w) / (w0 + w));
	    hsw->SetBinContent(i, w0 + w);
	  }
	} // for i
	ms["sum"][cv][cd] = hs;
	ms["sum2"][cv][cd] = hs2;
	ms["sumw"][cv][cd] = hsw;
      } // for is in sets
    } // for iv in vars
  } // for id in dirs

  // Write out results
  for (unsigned int id = 0; id != dirs.size(); ++id) {
    for (unsigned int iv = 0; iv != vars.size(); ++iv) {
    
      // Index to string for maps
      string sv = vars[iv];
      const char *cd = dirs[id].c_str();
      const char *cv = vars[iv].c_str();

      hs = ms["sum"][cv][cd]; assert(hs);
      hs2 = ms["sum2"][cv][cd]; assert(hs2);
      curdir->cd();
      TH1D *hs3 = (TH1D*)hs2->Clone(Form("%st",hs2->GetName()));
      hs3->SetLineColor(hs2->GetLineColor());
      hs3->SetMarkerColor(hs2->GetMarkerColor());

      // Calculate hs2 with RMS as error
      for (int i = 1; i != hs->GetNbinsX()+1; ++i) {
	double y = hs->GetBinContent(i);
	double y2 = hs2->GetBinContent(i);
	double ey = hs->GetBinError(i);
	hs2->SetBinContent(i, y);
	hs2->SetBinError(i, sqrt(fabs(y2 - y*y)));
	hs3->SetBinContent(i, y);
	hs3->SetBinError(i, sqrt(fabs(y2 - y*y) + ey*ey));
      } // for i

      TGraphErrors *g = new TGraphErrors(hs3);
      g->SetLineColor(hs3->GetLineColor());
      g->SetMarkerColor(hs3->GetMarkerColor());
      for (int i = g->GetN()-1; i!=-1; --i) {
	if (g->GetEY()[i]==0) g->RemovePoint(i);
      } // for i

      f->cd(Form("%s/eta00-13",cd));
      g->SetName(Form("%s_%s_ren",cv,"cmb"));
      g->Write(Form("%s_%s_ren",cv,"cmb"),TObject::kOverwrite);
      hs->SetName(Form("%s_%s_ren_stat",cv,"cmb"));
      hs->Write(Form("%s_%s_ren_stat",cv,"cmb"),TObject::kOverwrite);
      hs2->SetName(Form("%s_%s_ren_syst",cv,"cmb"));
      hs2->Write(Form("%s_%s_ren_syst",cv,"cmb"),TObject::kOverwrite);
      hs3->SetName(Form("%s_%s_ren_tot",cv,"cmb"));
      hs3->Write(Form("%s_%s_ren_tot",cv,"cmb"),TObject::kOverwrite);
      curdir->cd();
    } // for iv in vars
  } // for id in dirs

  f->Write();
  f->Close();
  curdir->cd();

} // globalFitPrecombinePF


// Purpose: pre-combine HD compositions for easier plotting / fitting
void globalFitPrecombineHDM(string s, bool addMJ, bool scaleRes) {

  TDirectory *curdir = gDirectory;

  const char *c = s.c_str();
  TFile *f = new TFile(Form("rootfiles/jecdata%s.root",c),"UPDATE");
  assert(f && !f->IsZombie());

  curdir->cd();

  // List all datasets providing PF composition inputs
  vector<string> sets;
  //if (s=="2017H") return;
  if (s=="2017H") {
    sets.push_back("zjet");
  }
  else {
    sets.push_back("gamjet");
    sets.push_back("zlljet");
    sets.push_back("zjet");
    sets.push_back("hadw");
  }
  if (s!="Run2Test") sets.push_back("incjet");
  // keep multijet last, as referenced to all those preceding it
  if (addMJ) sets.push_back("multijet");

  // Directories to loop for reading in data
  vector<string> dirs;
  //dirs.push_back("data");
  //dirs.push_back("mc");
  dirs.push_back("ratio");

  // zlljet: 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 85, 105, 130, 175, 230, 300, 400, 500, 700, 1000, 1500
  // gamjet: 15, 20, 25, 30, 35, 40, 50, 60, 70, 85, 105, 130, 175, 230, 300, 400, 500, 600, 700, 850, 1000, 1200, 1450, 1750,
  // multijet: 10, 15, 21, 28, 37, 49, 64, 84, 114, 153, 196, 245, 300, 362, 430, 507, 592, 686, 790, 905, 1032, 1172, 1327, 1497, 1684, 1890, 2116, 2366, 2640, 2941, 3273, 5220,
  // pfjet?: 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890,
  const double x[] =
    {15, 20, 25, 30, 35, 40, 45, 50, 60, // Z/W only
     70, 85, 130, // +gamjet
     175, 230, 300, 400, 500, // still Z/gamma domain
     592, 686, 790, 905, 1032, 1172, 1327, 1497, // multijet
     1684, 1890, 2116, 2366, 2640, 2941, 3273, 5220}; // multijet
  const int nx = sizeof(x)/sizeof(x[0])-1;
  TH1D *h1 = new TH1D(Form("h1_hdm_%s",c),"",nx,x);
  TF1 *f1 = new TF1(Form("f1hdm_%s",c),"[0]",15,7000);

  // Read in data and calculate sums
  TH1D *h(0), *hc(0), *hs(0), *hs2(0), *hsw(0);
  map<string, map<string, TH1D*> > m;
  for (unsigned int id = 0; id != dirs.size(); ++id) {
    for (unsigned int is = 0; is != sets.size(); ++is) {
    
      cout << "." << flush;
      
      // Index to string for maps
      const char *cd = dirs[id].c_str();
      const char *cs = sets[is].c_str();
      string ss = cs;

      // Read in data
      h = (TH1D*)f->Get(Form("%s/eta00-13/hdm_mpfchs1_%s",cd,cs));
      if (!h) h = (TH1D*)f->Get(Form("%s/eta00-13/hdm_%s",cd,cs));
      if (!h) cout << "Missing "<<cd<<"/hdm[X]_"<<cs<<endl<<flush;
      assert(h);
      m[cs][cd] = h;

      // Calculate sums
      hs  = m["sum"][cd];
      hs2 = m["sum2"][cd];
      hsw = m["sumw"][cd];
      if (!hs)  hs  = (TH1D*)h1->Clone(Form("%s_sum",h->GetName()));
      if (!hs2) hs2 = (TH1D*)h1->Clone(Form("%s_sum2",h->GetName()));
      if (!hsw) hsw = (TH1D*)h1->Clone(Form("%s_sumw",h->GetName()));
      
      hs->SetLineColor(h->GetLineColor());
      hs2->SetLineColor(h->GetLineColor());
      hsw->SetLineColor(h->GetLineColor());

      hs->SetMarkerColor(h->GetMarkerColor());
      hs2->SetMarkerColor(h->GetMarkerColor());
      hsw->SetMarkerColor(h->GetMarkerColor());

      for (int i = 1; i != hs->GetNbinsX()+1; ++i) {
	double ptmin = hs->GetBinLowEdge(i);
	double ptmax = hs->GetBinLowEdge(i+1);
	f1->SetRange(ptmin,ptmax);
	f1->SetParameter(0,0); // in case no points to fit
	f1->SetParError(0,0);  // in case no points to fit
	h->Fit(f1,"QRN");
	if (f1->GetParError(0)>0) { // check that something was fit
	  double w0 = hsw->GetBinContent(i);
	  double yw0 = hs->GetBinContent(i) * w0;
	  double y2w0 = pow(hs->GetBinContent(i),2) * w0;
	  double eyw0 = hs->GetBinError(i) * w0;
	  //
	  double w = 1./pow(f1->GetParError(0),2);
	  double yw = f1->GetParameter(0) * w;
	  double y2w = pow(f1->GetParameter(0),2) * w;
	  double eyw = f1->GetParError(0) * w;
	  if (ss=="multijet") {
	    TGraphErrors *gc =
	      (TGraphErrors*)f->Get(Form("%s/eta00-13/crecoil_%s_a30",cd,cs));
	    assert(gc);
	    double ptlead = hs->GetBinCenter(i);
	    double ptrec = gc->Eval(ptlead) * ptlead;
	    int j = hs->FindBin(ptrec);
	    double yref = hs->GetBinContent(j);
	    double eyref = hs->GetBinError(j);
	    double ey = sqrt(pow(f1->GetParError(0)*yref,2) + pow(eyref,2));
	    w = 1./pow(ey,2);
	    yw = f1->GetParameter(0) * yref * w;
	    y2w = pow(f1->GetParameter(0) * yref,2) * w;
	    eyw = ey * w;
	  }
	  hs->SetBinContent(i, (yw0 + yw) / (w0 + w));
	  hs->SetBinError(i, (eyw0 + eyw) / (w0 + w));
	  hs2->SetBinContent(i, (y2w0 + y2w) / (w0 + w));
	  hsw->SetBinContent(i, w0 + w);
	}
      } // for i
      m["sum"][cd] = hs;
      m["sum2"][cd] = hs2;
      m["sumw"][cd] = hsw;
    } // for is in sets
  } // for id in dirs

  // Write out results
  for (unsigned int id = 0; id != dirs.size(); ++id) {
    
    // Index to string for maps
    const char *cd = dirs[id].c_str();

    hs = m["sum"][cd]; assert(hs);
    hs2 = m["sum2"][cd]; assert(hs2);
    curdir->cd();

    // Optional L2L3Res scaling
    if (scaleRes) {
      hc = (TH1D*)f->Get("ratio/eta00-13/herr_l2l3res"); assert(hc);
      //hs->Multiply(hc);
      //hs2->Multiply(hc);
      //hs2->Multiply(hc);
      for (int i = 1; i != hs->GetNbinsX()+1; ++i) {
	double pt = hs->GetBinCenter(i);
	int j = hc->FindBin(pt);
	hs->SetBinContent(i,hs->GetBinContent(i)*hc->GetBinContent(j));
	hs2->SetBinContent(i,hs2->GetBinContent(i)*pow(hc->GetBinContent(j),2));
      } // for i
    } // scaleRes

    TH1D *hs3 = (TH1D*)hs2->Clone(Form("%s_tot",hs2->GetName()));
    hs3->SetLineColor(hs2->GetLineColor());
    hs3->SetMarkerColor(hs2->GetMarkerColor());

    // Calculate hs2 with RMS as error
    for (int i = 1; i != hs->GetNbinsX()+1; ++i) {
      double y = hs->GetBinContent(i);
      double y2 = hs2->GetBinContent(i);
      double ey = hs->GetBinError(i);
      hs2->SetBinContent(i, y);
      hs2->SetBinError(i, sqrt(fabs(y2 - y*y)));
      hs3->SetBinContent(i, y);
      hs3->SetBinError(i, sqrt(fabs(y2 - y*y) + ey*ey));
    } // for i

    f->cd(Form("%s/eta00-13",cd));
    //const char *cs = (addMJ ? "cmb_mj" : "cmb");
    string ss = "cmb";
    if (addMJ)  ss += "_mj";
    if (scaleRes) ss += "_res";
    const char *cs = ss.c_str();
    hs->SetName(Form("%s_%s","hdm",cs));
    hs->Write(Form("%s_%s","hdm",cs),TObject::kOverwrite);
    hs2->SetName(Form("%s_%s_syst","hdm",cs));
    hs2->Write(Form("%s_%s_syst","hdm",cs),TObject::kOverwrite);
    hs3->SetName(Form("%s_%s_tot","hdm",cs));
    hs3->Write(Form("%s_%s_tot","hdm",cs),TObject::kOverwrite);
    curdir->cd();
  } // for id in dirs

  f->Write();
  f->Close();
  curdir->cd();

} // globalFitPrecombineHDM
