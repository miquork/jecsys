// Purpose: Create systematic variations for global fit.
//          These are stored in the sys/ folder in jecdataX.root
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TProfile.h"
#include "TGraphErrors.h"

bool doMultijetJER = true;
// Which binning used for multijets ("leading" or "recoil");
//string multijetModeS = "leading"; // defined in reprocess.C => check consistency
string multijetModeS = "ptave"; // defined in reprocess.C => check consistency

bool doGamMass = true; // photon scale from Zee
bool doGamGains = true; // photon scale with gain1, gain6 path

bool doHadWfitProb = true; // hadronic W fitProb
bool doHadWfsr = true; // hadronic W FSR uncertainty

bool doIncjetHDM = true; // TH1D "hdm_incjet" for globalFitRenormPF.C

void globalFitSyst(string run = "Run2Test") {

  TDirectory *curdir = gDirectory;

  TFile *f = new TFile(Form("rootfiles/jecdata%s.root",run.c_str()),
		       "UPDATE");
  assert(f && !f->IsZombie());

  if (run=="2017H") { 
    doMultijetJER = doGamMass = doGamGains = doHadWfitProb = doHadWfsr = false;
  }
  
  assert(f->cd("ratio"));
  TDirectory *d1 = f->GetDirectory("ratio");
  assert(d1->cd("eta00-13"));
  TDirectory *d2 = d1->GetDirectory("eta00-13");
  assert(d2->cd("fsr"));
  TDirectory *dfsr = d2->GetDirectory("fsr");
  //
  if (!d2->FindObject("sys")) d2->mkdir("sys");
  assert(d2->cd("sys"));
  TDirectory *dsys = d2->GetDirectory("sys"); dsys->cd();

  curdir->cd();

  if (doMultijetJER) {
    // Systematics fit done in minitools/systMultijet.C
    // Very similar sensitivity in all variants of
    // (recoil, leading) x (MPF, MPFtype2, MJB) as expected
    // NB: assumes flat JER SF uncertainty. To account for difference in pT,
    // add the systematic from comparing recoil vs lead (necessary for high pT?)
    // So far: leading without systematics (except FSR) is chi2/NDF=125.4/52
    
    // Fits done in minitools/systMultijet.C on UL17BCDEF
    TF1 *frecoil = new TF1("frecoil","[0]+[1]*pow(0.01*x,[2])",114,2000);
    frecoil->SetParameters(0.138, -0.993, -0.561);
    
    TF1 *fleading = new TF1("fleading","[0]+[1]*pow(0.01*x,[2])",114,2000);
    fleading->SetParameters(0.192, -0.384, -3.260);
    
    // histogram to get correct binning
    //TH1D *h = (TH1D*)dfsr->Get("hkfsr_ptchs_multijet"); assert(h);
    TH1D *h = (TH1D*)dfsr->Get("hkfsr3_ptchs_multijet"); assert(h);

    TH1D *h1 = (TH1D*)h->Clone("mpfchs1_multijet_src0");
    TH1D *h2 = (TH1D*)h->Clone("ptchs_multijet_src0");
    h1->Reset();
    h2->Reset();
    
    for (int i = 1; i != h->GetNbinsX()+1; ++i) {
      double pt = h->GetBinCenter(i);
      if (pt>114 && pt<2000) { // NB, hard-coded for now (use counts later?)
	if (multijetModeS=="leading") {
	  h1->SetBinContent(i, 0.01*fleading->Eval(pt));
	  h2->SetBinContent(i, 0.01*fleading->Eval(pt));
	}
	else if (multijetModeS=="recoil") {
	  h1->SetBinContent(i, 0.01*frecoil->Eval(pt));
	  h2->SetBinContent(i, 0.01*frecoil->Eval(pt));
	}
	// placeholder for actual updated JER
	else if (multijetModeS=="ptave") {
	  h1->SetBinContent(i, 0.01*0.5*(fleading->Eval(pt)+frecoil->Eval(pt)));
	  h2->SetBinContent(i, 0.01*0.5*(fleading->Eval(pt)+frecoil->Eval(pt)));
	}
      }
    } // for i
    
    dsys->cd();
    //h1->Write(h1->GetName(), TObject::kOverwrite);
    //h2->Write(h2->GetName(), TObject::kOverwrite);
    h1->Write("jer_multijet", TObject::kOverwrite);
    curdir->cd();
  } // do MultijetJER


  if (doGamMass) {

    // histogram to get correct binning
    //TH1D *hk = (TH1D*)dfsr->Get("hkfsr_mpfchs1_gamjet"); assert(hk);      
    TH1D *hk = (TH1D*)dfsr->Get("hkfsr3_mpfchs1_gamjet"); assert(hk);      
    
    // Copied from reprocess.C
    // Extra functions for gamma+jet, modified from Zee+jet
    TF1 *f1mgam = new TF1("f1mgam","[0]+[1]*log(0.01*x)+[2]*pow(log(0.01*x),2)",
			  15,3000);
    TF1 *f1egam = new TF1("f1egam","sqrt([0]+pow(log(0.01*x),2)*[1]"
			  "+pow(log(0.01*x),4)*[2]"
			  "+2*log(0.01*x)*[3]+2*pow(log(0.01*x),2)*[4]"
			  "+2*pow(log(0.01*x),3)*[5])",
			  15,3000);
    
    // Copied from minitools/drawZmass.C
    //"      f1ezee->SetParameters(%+9.3g, %+9.3g, %+9.3g,\n"
    //"			    %+9.3g, %+9.3g, %+9.3g);\n"
    //       emat[0][0], emat[1][1], emat[2][2],
    //       emat[0][1], emat[0][2], emat[1][2]);

    // Copied from reprocess.C
    // 17Nov17_V10 BCDEF (EOY2017)
    f1egam->SetParameters(+3.73e-08, +1.17e-07, +2.02e-07,
    			  +1.5e-08, -5.07e-08, -6.39e-08);

    // Create error matrix
    const int n = 3;
    TMatrixD emat(n,n);
    emat[0][0] = f1egam->GetParameter(0);
    emat[1][1] = f1egam->GetParameter(1);
    emat[2][2] = f1egam->GetParameter(2);
    emat[0][1] = f1egam->GetParameter(3);
    emat[1][0] = f1egam->GetParameter(3);
    emat[0][2] = f1egam->GetParameter(4);
    emat[2][0] = f1egam->GetParameter(4);
    emat[1][2] = f1egam->GetParameter(5);
    emat[2][1] = f1egam->GetParameter(5);

    // Factorize error matrix into eigenvectors
    // Remember: A = Q*Lambda*Q^-1, where
    // A is emat, Q is eigmat, and Lambda is a diagonal matrix with
    // eigenvalues from eigvec on the diagonal. For eigenmatrix
    // Q^-1 = Q^T, i.e. inverse matrix is the original transposed
    TVectorD eigvec(n);
    TMatrixD eigmat = emat.EigenVectors(eigvec);

    // Eigenvectors are the columns and sum of eigenvectors squared
    // equals original uncertainty. Calculate histograms from the
    // eigenvectors and store them
    TF1 *fkeig = (TF1*)f1mgam->Clone("fkeig");
    for (int ieig = 0; ieig != n; ++ieig) {

      // Eigenvector functions
      for (int i = 0; i != n; ++i) {
	//fkeig->SetParameter(i, fk->GetParameter(i)
	//		    + eigmat[i][ieig] * sqrt(eigvec[ieig]));
	fkeig->SetParameter(i, eigmat[i][ieig] * sqrt(eigvec[ieig]));
      }
    
      // Eigenvector histograms evaluated at bin mean pT
      //TH1D *hke = (TH1D*)hk->Clone(Form("%s_eig%d",hk->GetName(),ieig));
      TH1D *hke = (TH1D*)hk->Clone(Form("zee_gamjet_eig%d",ieig));
      hke->Reset();
      
      for (int i = 1; i != hk->GetNbinsX()+1; ++i) {

	double pt = hk->GetBinCenter(i);
	// uncertainty sources are signed
	hke->SetBinContent(i, fkeig->Eval(2*pt));
	//hke->SetBinError(i, fabs(fkeig->Eval(pt)));
      }
      dsys->cd();
      hke->Write(hke->GetName(), TObject::kOverwrite);
      curdir->cd();
    } // for ieig
  } // doGamMass

  // Estimate impact on photon+jet from gain1 path
  if (doGamGains) {
    TFile *f = new TFile("../gamjet/files/GamHistosFill_data_2018ABCD_v19.root");
    assert(f && !f->IsZombie());

    TProfile *p1 = (TProfile*)f->Get("control/pgain1vspt"); assert(p1);
    TH1D *h1 = p1->ProjectionX("gamjet_gain1");
    h1->Scale(0.01);

    TProfile *p6 = (TProfile*)f->Get("control/pgain6vspt"); assert(p6);
    TH1D *h6 = p6->ProjectionX("gamjet_gain6");
    h6->Scale(0.01);

    dsys->cd();
    h1->Write(h1->GetName(), TObject::kOverwrite);
    h6->Write(h6->GetName(), TObject::kOverwrite);
    curdir->cd();

    // Create an extra "pure HDM variant" for globalFitRun2.root
    TH1D *hg = (TH1D*)d2->Get("hdm_mpfchs1_gamjet"); assert(hg);
    hg = (TH1D*)hg->Clone("hdm_gamjet");
    //hg->Add(h1,-0.5);
    // Do by hand to avoid empty bins
    for (int i = 1; i != hg->GetNbinsX()+1; ++i) {
      if (hg->GetBinContent(i)!=0) {
	int j = h1->FindBin(hg->GetBinCenter(i));
	// 0.5 to match zjet, 0 to match Zll+jet
	hg->SetBinContent(i, hg->GetBinContent(i)
			  - 0.0*h1->GetBinContent(j)
			  - 0.0*h6->GetBinContent(j));
			  //- 0.5*h1->GetBinContent(j)
			  //- 0.5*h6->GetBinContent(j));
      }
    } // for i in hg

    d2->cd();
    hg->Write(hg->GetName(), TObject::kOverwrite);
    curdir->cd();
  }

  if (doHadWfitProb) {
    // fitProb scanning done in minitools/hadW.C::DrawFP()

    // histogram to get correct binning
    TH1D *hk = (TH1D*)d2->Get("counts_hadw_a30"); assert(hk);

    // Fit from minitools/hadW.C::DrawFP()
    //TF1 *fhadw = new TF1("fse","[0]+[1]/(log(x/0.218) * x)",30,200);
    //fhadw->SetParameters(0, 125.54);

    TF1 *fhadw_a = new TF1("fhadw_ptave","[0]+[1]/(log(x/0.218) * x)",30,200);
    //fhadw_a->SetParameters(0.20652, 86.857); // UL17
    fhadw_a->SetParameters(0.95564, -139); // UL18
    TF1 *fhadw_b = new TF1("fhadw_ptboth","[0]+[1]/(log(x/0.218) * x)",30,200);
    //fhadw_b->SetParameters(0, 125.54); // UL17
    fhadw_b->SetParameters(0.99556, -128.37); // UL18

    TF1 *fhadw2_a = new TF1("fhadw2_ptave","[0]+fabs([1])/(log(x/0.218))",
			    30,200);
    fhadw2_a->SetParameters(0.57464, -2.104e-10);
    TF1 *fhadw2_b = new TF1("fhadw2_ptboth","[0]+fabs([1])/(log(x/0.218))",
			    30,200);
    fhadw2_b->SetParameters(0.65814, -4.1327e-07);
    
    TH1D *hkea = (TH1D*)hk->Clone("hadw_ptave_fitprob");
    hkea->Reset();
    TH1D *hkeb = (TH1D*)hk->Clone("hadw_ptboth_fitprob");
    hkeb->Reset();

    TH1D *hkea2 = (TH1D*)hk->Clone("hadw_ptave_fitprob2");
    hkea2->Reset();
    TH1D *hkeb2 = (TH1D*)hk->Clone("hadw_ptboth_fitprob2");
    hkeb2->Reset();
      
    // uncertainty sources are signed, and in %'s
    for (int i = 1; i != hkea->GetNbinsX()+1; ++i) {
      double pt = hkea->GetBinCenter(i);
      hkea->SetBinContent(i, 0.01*fhadw_a->Eval(pt));
      hkea2->SetBinContent(i, 0.01*(fhadw2_a->Eval(pt)-fhadw_a->Eval(pt)));
    }
    for (int i = 1; i != hkeb->GetNbinsX()+1; ++i) {
      double pt = hkeb->GetBinCenter(i);
      hkeb->SetBinContent(i, 0.01*fhadw_b->Eval(pt));
      hkeb2->SetBinContent(i, 0.01*(fhadw2_b->Eval(pt)-fhadw_b->Eval(pt)));
    }

    dsys->cd();
    //hke->Write("hadw_fitprob", TObject::kOverwrite);
    hkea->Write("hadw_ptave_fitprob", TObject::kOverwrite);
    hkeb->Write("hadw_ptboth_fitprob", TObject::kOverwrite);
    hkea2->Write("hadw_ptave_fitprob2", TObject::kOverwrite);
    hkeb2->Write("hadw_ptboth_fitprob2", TObject::kOverwrite);
    curdir->cd();
  } // doHadWfitProb

  if (doHadWfsr) {
    // FSR uncertainty using PSWeights (only in UL16V7)
    // Variants: FSRup (mu=0.5) and FSRdw (mu=2.0)
    TFile *ffsr = new TFile("rootfiles/hadWMC16GHV7_JEC.root","READ");
    assert(ffsr && !ffsr->IsZombie());

    TProfile *pmc = (TProfile*)ffsr->Get("pmm13bptave"); assert(pmc);
    TProfile *pup = (TProfile*)ffsr->Get("pmm13bptave_FSRup"); assert(pup);
    TProfile *pdw = (TProfile*)ffsr->Get("pmm13bptave_FSRdw"); assert(pdw);

    // Expected about 1-22 sigma shift on FSR => 1 sigma shift
    // Uncertainty for ptave
    TH1D *hfsr = pup->ProjectionX("hadw_fsr");
    hfsr->Add(pdw,-1);
    //hfsr->Scale(1.0);//2.0);
    //hfsr->Scale(1.5); // good match to both zjet and zlljet
    hfsr->Scale(1.0); // better match to zjet and gamjet
    hfsr->Add(pmc,+1);
    hfsr->Divide(pmc);

    // Store for use in globalFitL3Res.C (adapt later in softrad3.C?)
    dfsr->cd();
    hfsr->Write("hadw_fsr", TObject::kOverwrite);
    curdir->cd();

    // Correct this for an "HDM" variant of W>qq' as well
    TGraphErrors *gw = (TGraphErrors*)d2->Get("mpfchs1_hadw_a30"); assert(gw);
    gw = (TGraphErrors*)gw->Clone("hdm_mpfchs1_hadw");
    TH1D *hw = (TH1D*)hfsr->Clone("hdm_hadw"); hw->Reset();
    for (int i = 0; i != gw->GetN(); ++i) {
      double pt = gw->GetX()[i];
      double fsr = hfsr->GetBinContent(hfsr->FindBin(pt));
      double newhadw = gw->GetY()[i]/fsr;
      gw->SetPoint(i, pt, newhadw);
      int j = hw->FindBin(pt);
      hw->SetBinContent(j, newhadw);
      hw->SetBinError(j, gw->GetEY()[i]/fsr);
      // Don't add hfsr uncertainty (yet): likely correlated for fsrup & fsrdw
    }

    // Store for use in globalFitRun2.C
    d2->cd();
    gw->Write("ghdm_hadw", TObject::kOverwrite);
    hw->Write("hdm_hadw", TObject::kOverwrite);
    curdir->cd();

    // Uncertainty for ptave
    TH1D *hfsra = pup->ProjectionX("hadw_ptave_fsr");
    hfsra->Add(pdw,-1);
    hfsra->Divide(pmc);

    // ptboth is a clone of ptave
    TH1D *hfsrb = pup->ProjectionX("hadw_ptboth_fsr");
    hfsrb->Add(pdw,-1);
    hfsrb->Divide(pmc);

    dsys->cd();
    hfsra->Write("hadw_ptave_fsr", TObject::kOverwrite);
    hfsrb->Write("hadw_ptboth_fsr", TObject::kOverwrite);
    curdir->cd();
  } // doHadWfsr

  if (doIncjetHDM) {
    // Create TH1D out for TGraphErrors for incjet and factor in Run2Test
    TFile *f2 = new TFile("rootfiles/jecdataRun2Test.root","READ");
    TH1D *href(0);
    //if (f2 && !f2->IsZombie()) {
    assert(f2 && !f2->IsZombie());
    href = (TH1D*)f2->Get("ratio/eta00-13/hdm_cmb_mj"); assert(href);
    //}

    TGraphErrors *g = (TGraphErrors*)d2->Get("mpfchs1_incjet_a100");
    assert(g);
    g = (TGraphErrors*)g->Clone("ghdm_incjet");
    vector<double> x;
    x.push_back(g->GetX()[0]-g->GetEX()[0]);
    for (int i = 0; i != g->GetN(); ++i) {
      x.push_back(g->GetX()[i]+g->GetEX()[i]);
    }
    TH1D *h = new TH1D("hdm_incjet","",x.size()-1,&x[0]);
    for (int i = 0; i != g->GetN(); ++i) {
      int j = h->FindBin(g->GetX()[i]);
      //h->SetBinContent(j, g->GetY()[i]);
      //h->SetBinError(j, g->GetEY()[i]);
      //if (href) {
      int k = href->FindBin(g->GetX()[i]);
      double cref = href->GetBinContent(k);
      h->SetBinContent(j, g->GetY()[i] * cref);
      h->SetBinError(j, g->GetEY()[i] * cref);
      g->SetPoint(i, g->GetX()[i], g->GetY()[i] * cref);
      //}
    }

    d2->cd();
    g->Write("ghdm_incjet",TObject::kOverwrite);
    h->Write("hdm_incjet",TObject::kOverwrite);
    curdir->cd();

    if (f2) f2->Close();
  } // doIncJetHDM
  
  f->Close();
} // globalFitSyst

