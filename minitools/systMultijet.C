// Purpose: Evaluate multijet systematic uncertainties from
//          FSR) Already done in softrad.C and incorporated in globalFitL3Res.C
//               but cross-check with Pt30 vs Pt15
//          JER) JERup vs JERdown
//               + cross-check with recoil vs leading binning (opposite bias?)
//          JEC) Some estimate of L2Res MPF vs pTbal vs eta for recoil jets?
//          PU)  Some estimate of L1Res vs eta for recoil jets?
//          Flavor) More gluon jets in recoil. Use FlavorPureX sources?
//          Unclustered pT) metType2 vs regular variant, for data/MC
//
// run with 'root -l -b -q systMultijet.C+g'
#include "TFile.h"
#include "TF1.h"
#include "TLine.h"

#include "../tdrstyle_mod15.C"

using namespace std;

const bool doFSR = true;
const bool doJER = true;
//const bool doTeV = true;
const double ptmin = 30;//114;

void systMultijets(string mode="", string spt="Pt30", string spt2="") {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TFile *f = new TFile("rootfiles/multijet_Rebin2_20200320_UL2017BCDEF_jecV1_jerV3.root","READ");
  assert(f && !f->IsZombie());

  TFile *freco = new TFile("rootfiles/jecdataBCDEF_MJrecoil.root","READ");
  assert(freco && !freco->IsZombie());
  freco->cd("ratio/eta00-13/fsr");
  TDirectory *dreco = gDirectory;

  TFile *flead = new TFile("rootfiles/jecdataBCDEF_MJleading_v2.root","READ");
  assert(flead && !flead->IsZombie());
  flead->cd("ratio/eta00-13/fsr");
  TDirectory *dlead = gDirectory;

  curdir->cd();

  if (doFSR) {

    TH1D *hmpfl = (TH1D*)dlead->Get("hkfsr_mpfchs1_multijet"); assert(hmpfl);
    TH1D *hmjbl = (TH1D*)dlead->Get("hkfsr_ptchs_multijet");   assert(hmjbl);
    TH1D *hmpfr = (TH1D*)dreco->Get("hkfsr_mpfchs1_multijet"); assert(hmpfr);
    TH1D *hmjbr = (TH1D*)dreco->Get("hkfsr_ptchs_multijet");   assert(hmjbr);

    TH1D *h = new TH1D("h_fsr",";p_{T,lead} or p_{T,recoil} (GeV);dR / d#alpha",
		       100,ptmin,2000);
    h->GetXaxis()->SetMoreLogLabels();
    h->GetXaxis()->SetNoExponent();
    h->SetMaximum(+0.06);
    h->SetMinimum(-0.03);

    TLine *l = new TLine();
    l->SetLineStyle(kDashed);

    TCanvas *c1 = tdrCanvas("c1_fsr",h,4,11,kSquare);
    c1->SetLogx();

    l->DrawLine(ptmin,0,2000,0);
    tdrDraw(hmjbr,"",kOpenSquare,kRed);
    tdrDraw(hmjbl,"",kFullSquare,kRed);  hmjbl->SetLineWidth(2);
    tdrDraw(hmpfr,"",kOpenCircle,kBlue); 
    tdrDraw(hmpfl,"",kFullCircle,kBlue); hmpfl->SetLineWidth(2);

    /*
    // Show original quad-log fit
    TF1 *fqlmjb = new TF1("fqlmjb","[0]+[1]*log(0.01*x)+[2]*pow(log(0.01*x),2)",
			  ptmin,2000);
    fqlmjb->SetParameters(0.05,-0.03,+0.01);
    fqlmjb->SetLineColor(kRed);
    hmjbl->Fit(fqlmjb,"QRN");
    //fqlmjb->Draw("SAME");

    TF1 *fqlmpf = new TF1("fqlmpf","[0]+[1]*log(0.01*x)+[2]*pow(log(0.01*x),2)",
			  ptmin,2000);
    fqlmpf->SetParameters(0.02,-0.01,+0.001);
    fqlmpf->SetLineColor(kBlue);
    hmpfl->Fit(fqlmpf,"QRN");
    //fqlmpf->Draw("SAME");

    // Show improved shapes:
    // MJB lead: proportional to alpha_s times 15 GeV / pT
    // MPF lead: proportional to MJB times unclustered pT response difference
    TF1 *ffsrmjb = new TF1("ffsrmjb",
			   "[0]*(1./log(x*x/([1]*[1]))+[1]*[1]/([1]*[1] - x*x))"
			   "/x",ptmin,2000);
    ffsrmjb->SetParameters(1, 0.218);
    ffsrmjb->SetLineColor(kRed);
    hmjbl->Fit(ffsrmjb,"QRN");
    ffsrmjb->Draw("SAME");
    TF1 *ffsrmpf = new TF1("ffsrmpf",
			   "[0]*(1./log(x*x/([1]*[1]))+[1]*[1]/([1]*[1] - x*x))"
			   "/x",ptmin,2000);
    //ffsrmpf->SetParameters(0.05,-0.03,+0.01);
    ffsrmpf->SetParameters(1, 0.218);
    ffsrmpf->FixParameter(1, 0.218);
    //ffsrmpf->SetParameters(ffsrmjb->GetParameter(0),ffsrmjb->GetParameter(1));
    //ffsrmpf->FixParameter(1,ffsrmjb->GetParameter(1));
    ffsrmpf->SetLineColor(kBlue);
    hmpfl->Fit(ffsrmpf,"QRN");
    ffsrmpf->Draw("SAME");
    */

    // Simplified version
    TF1 *ffsrmjb2 = new TF1("ffsrmjb2",
			    "[0] / (log(x/0.218) * x)", ptmin,2000);
    ffsrmjb2->SetParameter(0, 45);
    ffsrmjb2->SetLineColor(kRed);
    hmjbl->Fit(ffsrmjb2,"RN");
    ffsrmjb2->Draw("SAME");
    TF1 *ffsrmpf2 = new TF1("ffsrmpf2",
			    "[0] / (log(x/0.218) * x)", ptmin,2000);
    ffsrmpf2->SetParameter(0, 37*0.2);
    ffsrmpf2->SetLineColor(kBlue);
    hmpfl->Fit(ffsrmpf2,"RN");
    ffsrmpf2->Draw("SAME");

    c1->SaveAs("pdf/systMultijet_FSR.pdf");
  }
    
  if (doJER) {

    TH1D *hup(0), *hmd(0), *hdw(0), *hupl(0), *hmdl(0), *hdwl(0);
    TH1D *hdt(0), *hdtl(0);
    string sm("");
    const char *cpt = spt.c_str();
    if (spt2=="") spt2 = spt;
    const char *cpt2 = spt2.c_str();

    if (mode=="MPF") {
      sm = "MPF";
      hup = (TH1D*)f->Get(Form("MC_JER-up/%s/MPF_recoil_MG",cpt2));
      hmd = (TH1D*)f->Get(Form("MC/%s/MPF_recoil_MG",cpt));
      hdw = (TH1D*)f->Get(Form("MC_JER-down/%s/MPF_recoil_MG",cpt2));
      //
      hupl = (TH1D*)f->Get(Form("MC_JER-up/%s/MPF_leading_MG",cpt2));
      hmdl = (TH1D*)f->Get(Form("MC/%s/MPF_leading_MG",cpt));
      hdwl = (TH1D*)f->Get(Form("MC_JER-down/%s/MPF_leading_MG",cpt2));
      //
      hdt = (TH1D*)f->Get(Form("Data/%s/MPF_recoil_L1L2Res",cpt));
      hdtl = (TH1D*)f->Get(Form("Data/%s/MPF_leading_L1L2Res",cpt));
    }
    if (mode=="MPFtype2") {
      hup = (TH1D*)f->Get(Form("MC_JER-up/%s/MPF_recoil_MG_metType2",cpt2));
      hmd = (TH1D*)f->Get(Form("MC/%s/MPF_recoil_MG_metType2",cpt));
      hdw = (TH1D*)f->Get(Form("MC_JER-down/%s/MPF_recoil_MG_metType2",cpt2));
      //
      hupl = (TH1D*)f->Get(Form("MC_JER-up/%s/MPF_leading_MG_metType2",cpt2));
      hmdl = (TH1D*)f->Get(Form("MC/%s/MPF_leading_MG_metType2",cpt));
      hdwl = (TH1D*)f->Get(Form("MC_JER-down/%s/MPF_leading_MG_metType2",cpt2));
      //
      hdt = (TH1D*)f->Get(Form("Data/%s/MPF_recoil_L1L2Res_metType2",cpt));
      hdtl = (TH1D*)f->Get(Form("Data/%s/MPF_leading_L1L2Res_metType2",cpt));
    }
    if (mode=="MJB") {
      hup = (TH1D*)f->Get(Form("MC_JER-up/%s/MJB_recoil_MG",cpt2));
      hmd = (TH1D*)f->Get(Form("MC/%s/MJB_recoil_MG",cpt));
      hdw = (TH1D*)f->Get(Form("MC_JER-down/%s/MJB_recoil_MG",cpt2));
      //
      hupl = (TH1D*)f->Get(Form("MC_JER-up/%s/MJB_leading_MG",cpt2));
      hmdl = (TH1D*)f->Get(Form("MC/%s/MJB_leading_MG",cpt));
      hdwl = (TH1D*)f->Get(Form("MC_JER-down/%s/MJB_leading_MG",cpt2));
      //
      hdt = (TH1D*)f->Get(Form("Data/%s/MJB_recoil_L1L2Res",cpt));
      hdtl = (TH1D*)f->Get(Form("Data/%s/MJB_leading_L1L2Res",cpt));
    }

    assert(hup);
    assert(hmd);
    assert(hdw);
    assert(hupl);
    assert(hmdl);
    assert(hdwl);
    assert(hdt);
    assert(hdtl);

    hup->GetXaxis()->SetRangeUser(ptmin,2000);
    hmd->GetXaxis()->SetRangeUser(ptmin,2000);
    hdw->GetXaxis()->SetRangeUser(ptmin,2000);
    //
    hupl->GetXaxis()->SetRangeUser(ptmin,2000);
    hmdl->GetXaxis()->SetRangeUser(ptmin,2000);
    hdwl->GetXaxis()->SetRangeUser(ptmin,2000);
    //
    hdt->GetXaxis()->SetRangeUser(ptmin,2000);
    hdtl->GetXaxis()->SetRangeUser(ptmin,2000);

    TH1D *hup2 = (TH1D*)hup->Clone("hup2_jer");
    TH1D *hmd2 = (TH1D*)hmd->Clone("hmd2_jer");
    TH1D *hdw2 = (TH1D*)hdw->Clone("hdw2_jer");
    //
    TH1D *hupl2 = (TH1D*)hupl->Clone("hupl2_jer");
    TH1D *hmdl2 = (TH1D*)hmdl->Clone("hmdl2_jer");
    TH1D *hdwl2 = (TH1D*)hdwl->Clone("hdwl2_jer");
    //
    TH1D *hdt2 = (TH1D*)hdt->Clone("hdt2_jer");
    TH1D *hdtl2 = (TH1D*)hdtl->Clone("hdtl2_jer");

    hup2->Divide(hmd); 
    hmd2->Divide(hmd);
    hdw2->Divide(hmd);
    //
    hupl2->Divide(hmdl); 
    hmdl2->Divide(hmdl);
    hdwl2->Divide(hmdl);
    //
    hdt2->Divide(hmd);
    hdtl2->Divide(hmdl);

    TH1D *hup3 = (TH1D*)hup2->Clone("hup3_jer"); hup3->Scale(100);
    TH1D *hmd3 = (TH1D*)hmd2->Clone("hmd3_jer"); hmd3->Scale(100);
    TH1D *hdw3 = (TH1D*)hdw2->Clone("hdw3_jer"); hdw3->Scale(100);
    hup3->Add(hmd2,-100);
    hdw3->Add(hmd2,-100);
    hmd3->Add(hmd2,-100);
    //
    TH1D *hupl3 = (TH1D*)hupl2->Clone("hupl3_jer"); hupl3->Scale(100);
    TH1D *hmdl3 = (TH1D*)hmdl2->Clone("hmdl3_jer"); hmdl3->Scale(100);
    TH1D *hdwl3 = (TH1D*)hdwl2->Clone("hdwl3_jer"); hdwl3->Scale(100);
    hupl3->Add(hmdl2,-100);
    hdwl3->Add(hmdl2,-100);
    hmdl3->Add(hmdl2,-100);
    //
    TH1D *hdt3 = (TH1D*)hdt2->Clone("hdt3_jer"); hdt3->Scale(100);
    TH1D *hdtl3 = (TH1D*)hdtl2->Clone("hdtl3_jer"); hdtl3->Scale(100);
    hdt3->Add(hmd2,-100);
    hdtl3->Add(hmdl2,-100);

    TH1D *hu1 = new TH1D(Form("hu1_jer_%s",mode.c_str()),
			Form(";p_{T,recoil} or p_{T,lead} (GeV);"
			     "%s ( #LTp_{T,lead}#GT / #LTp_{T,recoil}#GT )",
			     mode.c_str()),
			100,ptmin,2000);
    hu1->GetXaxis()->SetMoreLogLabels();
    hu1->GetXaxis()->SetNoExponent();
    //hu1->SetMaximum(1.01-1e-5); // MPF_recoil only
    //hu1->SetMaximum(1.06-1e-5); // MPF recoil and lead
    hu1->SetMaximum(1.12-1e-5); // MPF and MJB
    //hu1->SetMinimum(0.93+1e-5); // MPF only 
    hu1->SetMinimum(0.89+1e-5); // MPF and MJB

    TH1D *hd1 = new TH1D(Form("hd1_jer_%s",mode.c_str()),
			";p_{T,recoil} or p_{T,lead} (GeV);"
			"Ratio-1 (%)",100,ptmin,2000);
    hd1->GetXaxis()->SetMoreLogLabels();
    hd1->GetXaxis()->SetNoExponent();
    hd1->SetMaximum(+1.1);//+0.8);
    hd1->SetMinimum(-1.1);//-0.8);

    lumi_13TeV = "UL2017, 41.5 fb^{-1}";
    TCanvas *c1 = tdrDiCanvas("c1_jer",hu1,hd1,4,11);

    c1->cd(1);
    gPad->SetLeftMargin(0.17);
    hu1->GetYaxis()->SetTitleOffset(1.2);
    gPad->SetLogx();

    tdrDraw(hmdl,"",kNone,kBlack,kSolid);//kDotted);
    tdrDraw(hupl,"",kNone,kRed,kSolid);//kDotted);
    tdrDraw(hdwl,"",kNone,kBlue,kSolid);//kDotted);
    tdrDraw(hdtl,"",kNone,kGreen+3,kSolid); hdtl->SetLineWidth(2);
    hmdl->SetLineWidth(2);
    hupl->SetLineWidth(2);
    hdwl->SetLineWidth(2);

    tdrDraw(hmd,"",kNone,kBlack);
    tdrDraw(hup,"",kNone,kRed);
    tdrDraw(hdw,"",kNone,kBlue);
    tdrDraw(hdt,"",kNone,kGreen+3); //hdt->SetLineWidth(2);

    //TLegend *leg = tdrLeg(0.20,0.60,0.50,0.75);
    //TLegend *leg = tdrLeg(0.20,0.35,0.50,0.55);
    TLegend *leg = tdrLeg(0.18,0.31,0.43,0.51); // MJB and MPF
    leg->AddEntry(hdt,"JER data","L");
    leg->AddEntry(hdw,"JER up","L");
    leg->AddEntry(hmd,"JER nom.","L");
    leg->AddEntry(hup,"JER down","L");

    c1->cd(2);
    gPad->SetLeftMargin(0.17);
    hd1->GetYaxis()->SetTitleOffset(0.62);//1.2);
    gPad->SetLogx();

    tdrDraw(hmd3,"",kNone,kBlack);
    tdrDraw(hup3,"",kNone,kRed);
    tdrDraw(hdw3,"",kNone,kBlue);
    tdrDraw(hdt3,"",kNone,kGreen+3); //hdt3->SetLineWidth(2);
    //hup3->Scale(-1);
    //
    tdrDraw(hmdl3,"",kNone,kBlack,kSolid);//kDotted);
    tdrDraw(hupl3,"",kNone,kRed,kSolid);//kDotted);
    tdrDraw(hdwl3,"",kNone,kBlue,kSolid);//xkDotted);
    tdrDraw(hdtl3,"",kNone,kGreen+3,kSolid); hdtl3->SetLineWidth(2);
    //hupl3->Scale(-1);
    hmdl3->SetLineWidth(2);
    hupl3->SetLineWidth(2);
    hdwl3->SetLineWidth(2);

    TF1 *f1dw = new TF1("f1dw","[0]+[1]*pow(0.01*x,[2])",ptmin,2000);
    f1dw->SetParameters(0.1,0.6,-0.5);
    hdw3->Fit(f1dw,"QRN");

    f1dw->SetLineColor(kBlue);
    f1dw->Draw("SAME");

    TF1 *f1up = new TF1("f1up","[0]+[1]*pow(0.01*x,[2])",ptmin,2000);
    f1up->SetParameters(-f1dw->GetParameter(0), -f1dw->GetParameter(1),
			f1dw->GetParameter(2));
    f1up->SetLineColor(kRed);
    f1up->Draw("SAME");

    TF1 *f1dwl = new TF1("f1dwl","[0]+[1]*pow(0.01*x,[2])",ptmin,2000);
    f1dwl->SetParameters(0.1,0.6,-0.5);
    hdwl3->Fit(f1dwl,"QRN");

    f1dwl->SetLineColor(kBlue);
    f1dwl->SetLineStyle(kSolid);//kDotted);
    f1dwl->SetLineWidth(2);
    f1dwl->Draw("SAME");

    TF1 *f1upl = new TF1("f1upl","[0]+[1]*pow(0.01*x,[2])",ptmin,2000);
    f1upl->SetParameters(-f1dwl->GetParameter(0), -f1dwl->GetParameter(1),
			f1dwl->GetParameter(2));
    f1upl->SetLineColor(kRed);
    f1upl->SetLineStyle(kSolid);//kDotted);
    f1upl->SetLineWidth(2);
    f1upl->Draw("SAME");

    if (mode=="MPF") {
      cout << "  // Fits done in minitools/systMultijet.C on UL17BCDEF" << endl;
      cout << Form("  TF1 *frecoil = new TF1(\"frecoil\","
		   "\"[0]+[1]*pow(0.01*x,[2])\",114,2000);\n"
		   "  frecoil->SetParameters(%1.3f, %1.3f, %1.3f);\n",
		   f1up->GetParameter(0), f1up->GetParameter(1),
		   f1up->GetParameter(2)) << endl;
      cout << Form("  TF1 *fleading = new TF1(\"fleading\","
		   "\"[0]+[1]*pow(0.01*x,[2])\",114,2000);\n"
		   "  fleading->SetParameters(%1.3f, %1.3f, %1.3f);\n",
		   f1upl->GetParameter(0), f1upl->GetParameter(1),
		   f1upl->GetParameter(2)) << endl;
    }

    c1->Update();
    c1->SaveAs(Form("pdf/systMultijet_JER_%s.pdf",mode.c_str()));


    // High pT separately
    TH1D *hu2 = new TH1D(Form("hu2_jer_%s",mode.c_str()),
			Form(";p_{T,recoil} or p_{T,lead} (GeV);"
			     "%s ( lead / recoil - 1 ) (%%)",
			     mode.c_str()),
			100,ptmin,2000);
    hu2->GetXaxis()->SetMoreLogLabels();
    hu2->GetXaxis()->SetNoExponent();
    hu2->SetMaximum(+1.7-1e-5);
    hu2->SetMinimum(-0.5+1e-5);

    TH1D *hd2 = new TH1D(Form("hd2_jer_%s",mode.c_str()),
			";p_{T,recoil} or p_{T,lead} (GeV);"
			"Ratio-1 (%)",100,ptmin,2000);
    hd2->GetXaxis()->SetMoreLogLabels();
    hd2->GetXaxis()->SetNoExponent();
    hd2->SetMaximum(+1.1);
    hd2->SetMinimum(-0.5);

    // Create hybrid of recoil and lead that cancels out JER uncertainty
    TH1D *hdth3 = (TH1D*)hdt3->Clone("hdth3"); hdth3->Reset();
    for (int i = 1; i != hdth3->GetNbinsX()+1; ++i) {

      double x = hdt3->GetBinCenter(i);
      double yr = hdt3->GetBinContent(i);
      double er = hdt3->GetBinError(i);
      double jr = f1up->Eval(x);
      //
      double yl = hdtl3->GetBinContent(i);
      double el = hdtl3->GetBinError(i);
      double jl = f1upl->Eval(x);
      //
      // jh=w*jl+(1-w)*jr=0 => w*(jl-jr)+jr=0 => w=jr/(jr-jl)
      double w = (jr!=jl ? jr/(jr-jl) : 0.5);
      double yh = w*yl + (1-w)*yr;
      double eh = w*el + (1-w)*er;
      hdth3->SetBinContent(i, yh);
      hdth3->SetBinError(i, eh);
    } // for i

    // Relative weight of C term in NSC JER
    // (try to estimate shape of pT-dependent JER uncertainty)
    TF1 *fcw = new TF1("fcw","(sqrt([0]*[3]/(x*x)+[1]*[4]/x+[2]*[5])"
		       "/sqrt([3]/(x*x)+[4]/x+[5])-1)/([6])",
		       ptmin,2000);
    fcw->SetParameters(pow(0.9,2),pow(0.9,2),pow(1.1,2),
		       pow(8,2),pow(0.8,2),pow(0.05,2),
		       0.1);



    TH1D *hdt3r = (TH1D*)hdt3->Clone("hdt3r");
    hdt3r->Reset();
    TH1D *hup3r = (TH1D*)hup3->Clone("hup3r");
    hup3r->Reset();
    TH1D *hup3rb = (TH1D*)hup3->Clone("hup3rb");
    hup3rb->Reset();

    for (int i = 1; i != hdt3r->GetNbinsX()+1; ++i) {
      hdt3r->SetBinContent(i, ((1+0.01*hdtl3->GetBinContent(i)) /
			       (1+0.01*hdt3->GetBinContent(i)) - 1)*100);
      hdt3r->SetBinError(i, sqrt(pow(hdtl3->GetBinError(i),2)+
				 pow(hdt3->GetBinError(i),2)));
      //
      hup3r->SetBinContent(i, ((1+0.01*hupl3->GetBinContent(i)) /
			       (1+0.01*hup3->GetBinContent(i)) - 1)*100);
      hup3r->SetBinError(i, sqrt(pow(hupl3->GetBinError(i),2)+
				 pow(hup3->GetBinError(i),2)));
      //
      double x = hup3r->GetBinCenter(i);
      hup3rb->SetBinContent(i, hup3r->GetBinContent(i)*fcw->Eval(x));
      hup3rb->SetBinError(i, hup3r->GetBinError(i)*fcw->Eval(x));
    }
    


    lumi_13TeV = "UL2017, 41.5 fb^{-1}";
    TCanvas *c2 = tdrDiCanvas("c2_jer",hu2,hd2,4,11);

    c2->cd(1);
    gPad->SetLogx();

    TLine *l = new TLine();
    l->SetLineStyle(kDashed);
    l->DrawLine(ptmin,0,2000,0);

    tdrDraw(hdt3,"",kNone,kGreen+3,kSolid);
    tdrDraw(hdtl3,"",kNone,kGreen+3,kSolid); hdtl3->SetLineWidth(2);
    tdrDraw(hdth3,"",kNone,kBlack,kSolid); //hdth3->SetLineWidth(2);

    c2->cd(2);
    gPad->SetLogx();
    
    l->DrawLine(ptmin,0,2000,0);
    
    tdrDraw(hup3r,"",kNone,kRed,kSolid);
    tdrDraw(hup3rb,"",kNone,kRed+2,kSolid);
    tdrDraw(hdt3r,"",kNone,kGreen+3,kSolid); hdt3r->SetLineWidth(2);

    fcw->Draw("SAME");

    c2->SaveAs(Form("pdf/systMultijet_JER_%s_TeV.pdf",mode.c_str()));
  } // doJER


} // systMultijet

void systMultijet() {

  /*
  systMultijets("MPF");
  systMultijets("MPFtype2");
  systMultijets("MJB");
  */
  
  // Check how methods compare with lower pT threshold
  // (but JET not yet implemented for 15, 20)
  systMultijets("MPF","Pt15","Pt30");
  systMultijets("MPFtype2","Pt15","Pt30");
  systMultijets("MJB","Pt15","Pt30");


} // sysMultijets
