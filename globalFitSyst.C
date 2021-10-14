// Purpose: Create systematic variations for global fit.
//          These are stored in the sys/ folder in jecdataX.root
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMatrixD.h"
#include "TVectorD.h"

const bool doMultijetJER = true;
// Which binning used for multijets ("leading" or "recoil");
//string multijetModeS = "leading"; // defined in reprocess.C => check consistency
string multijetModeS = "ptave"; // defined in reprocess.C => check consistency

const bool doGamMass = true; // photon scale from Zee

const bool doHadWfitProb = true; // hadronic W fitProb

void globalFitSyst(string run = "BCDEF") {

  TFile *f = new TFile(Form("rootfiles/jecdata%s.root",run.c_str()),
		       "UPDATE");
  assert(f && !f->IsZombie());
  
  assert(f->cd("ratio"));
  TDirectory *d1 = f->GetDirectory("ratio");
  assert(d1->cd("eta00-13"));
  TDirectory *d2 = d1->GetDirectory("eta00-13");
  assert(d2->cd("fsr"));
  TDirectory *dfsr = d2->GetDirectory("fsr");
  //
  if (!d2->FindObject("sys")) d2->mkdir("sys");
  assert(d2->cd("sys"));
  TDirectory *d3 = d2->GetDirectory("sys"); d3->cd();
  
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
    TH1D *h = (TH1D*)dfsr->Get("hkfsr_ptchs_multijet"); assert(h);

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
    
    d3->cd();
    //h1->Write(h1->GetName(), TObject::kOverwrite);
    //h2->Write(h2->GetName(), TObject::kOverwrite);
    h1->Write("jer_multijet", TObject::kOverwrite);
  } // do MultijetJER


  if (doGamMass) {

    // histogram to get correct binning
    TH1D *hk = (TH1D*)dfsr->Get("hkfsr_mpfchs1_gamjet"); assert(hk);      
    
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
      d3->cd();
      hke->Write(hke->GetName(), TObject::kOverwrite);
    } // for ieig
  } // doGamMass


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

    d3->cd();
    //hke->Write("hadw_fitprob", TObject::kOverwrite);
    hkea->Write("hadw_ptave_fitprob", TObject::kOverwrite);
    hkeb->Write("hadw_ptboth_fitprob", TObject::kOverwrite);
    hkea2->Write("hadw_ptave_fitprob2", TObject::kOverwrite);
    hkeb2->Write("hadw_ptboth_fitprob2", TObject::kOverwrite);
  } // doHadWfitProb
  
  f->Close();
} // globalFitSyst

