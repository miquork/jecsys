// Purpose: draw mass of jets from W to estimate jet mass scale bias
// Uses output from minitools/mk_hadW.C
// run with 'root -l minitools/drawWjetmass.C+g'
#include "TFile.h"
#include "TProfile.h"
#include "TF1.h"
#include "TGraphErrors.h"

#include "../tdrstyle_mod15.C"

void drawBjetMass(string mode); // Check also b-jet mass etc

void drawWjetMass(string mode="18V5") {

  //drawBjetMass(mode);
  //continue;

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  const char *cm = mode.c_str();

  TFile *fm = new TFile(Form("rootfiles/hadWMC%s_MPDGcorrNoW.root",cm),"READ");
  //TFile *fm = new TFile("rootfiles/hadWMC18V4_MPDcorrNoW.root","READ");
  //TFile *fm = new TFile("rootfiles/hadWMC18V4_MPcorrNoW.root","READ");
  //TFile *fm = new TFile("rootfiles/hadWMC18V4_McorrNoW.root","READ");
  assert(fm && !fm->IsZombie());

  TFile *fd = new TFile(Form("rootfiles/hadWUL%s_MPDGcorrNoW.root",cm),"READ");
  //TFile *fd = new TFile("rootfiles/hadWUL18V4_MPDcorrNoW.root","READ");
  //TFile *fd = new TFile("rootfiles/hadWUL18V4_MPcorrNoW.root","READ");
  //TFile *fd = new TFile("rootfiles/hadWUL18V4_McorrNoW.root","READ");
  assert(fd && !fd->IsZombie());

  curdir->cd();

  TProfile *pmreco = (TProfile*)fm->Get("pmreco"); assert(pmreco);
  TProfile *pmgen = (TProfile*)fm->Get("pmgen");   assert(pmgen);
  TProfile *pmm = (TProfile*)fm->Get("pmjet");     assert(pmm);
  TProfile *pmd = (TProfile*)fd->Get("pmjet");     assert(pmd);

  TProfile *pavereco = (TProfile*)fm->Get("pavereco"); assert(pavereco);
  TProfile *pavegen = (TProfile*)fm->Get("pavegen");   assert(pavegen);

  double ptave = pavereco->GetMean();

  TH1D *hmreco = pmreco->ProjectionX("hmreco"); // <mjet>/<ptave>
  hmreco->Divide(pmgen);  // divide by <mgen>/<ptave> => <mjet>/<mgen>
  TH1D *hmmc = pmm->ProjectionX("hmmc"); // <mjet>/<ptave>
  hmmc->Divide(pmgen); // divide by <mgen>/<ptave> => <mjet>/<mgen>
  TH1D *hmdata = pmd->ProjectionX("hmdata");
  hmdata->Divide(pmgen);

  TH1D *havegen = pavegen->ProjectionX("havegen");
  //haver->Divide(pavegen);
  TGraphErrors *gavegen = new TGraphErrors(0);
  for (int i = 1; i != havegen->GetNbinsX()+1; ++i) {
    double xmid = havegen->GetBinCenter(i);
    if (xmid<30 || xmid>230) continue;
    double x = pavereco->GetBinContent(pavereco->FindBin(xmid));
    double ex = pavereco->GetBinError(pavereco->FindBin(xmid));
    int n = gavegen->GetN();
    gavegen->SetPoint(n, x, havegen->GetBinContent(i));
    gavegen->SetPointError(n, ex, havegen->GetBinError(i));
  }

  TH1D *hmreco2 = (TH1D*)hmreco->Clone("hmreco2");
  hmreco2->Multiply(havegen);
  TH1D *hmmc2 = (TH1D*)hmmc->Clone("hmmc2");
  hmmc2->Multiply(havegen); // divide by <ptgen>/<ptave>
                            //  => (<mjet>/<ptave>) / (<mgen>/<ptgen>)
  TH1D *hmdata2 = (TH1D*)hmdata->Clone("hmdata2");
  hmdata2->Multiply(havegen);

  TH1D *hmgen = pmgen->ProjectionX("hmgen");
  hmgen->Divide(pavegen);
  TH1D *hmrecorc = pmreco->ProjectionX("hmrecorc");
  //hmreco->Divide(pavereco); // already divided
  TH1D *hmrecomc = pmm->ProjectionX("hmrecomc");  
  TH1D *hmrecodt = pmd->ProjectionX("hmrecodt");  

  TH1D *hm = (TH1D*)hmdata2->Clone("hm");
  hm->Add(hmmc2,-1); // are stat. unc. correct with add+divide?
  hm->Divide(hmmc2);
  hm->Scale(100.);

  TH1D *hm2 = (TH1D*)hmrecodt->Clone("hm2");
  hm2->Add(hmrecomc,-1); // are stat. unc. correct with add+divide?
  hm2->Divide(hmrecomc);
  hm2->Scale(100.);

  // Correction for reco mass/pt ratio
  TF1 *fmj = new TF1("fmj","1+[1]*log(x/[0])",30,230);
  fmj->SetParameters(ptave,0.05);
  fmj->SetLineColor(kRed);
  hmmc2->Fit(fmj,"QRN");
  //hmreco2->Fit(fmj,"QRN");
  cout << Form("  fmj->SetParameters(%1.3g, %1.4g); // ptave = %1.3g GeV",
	       fmj->GetParameter(0), fmj->GetParameter(1), ptave) << endl;

  // Correction for pTgen/pTave ratio
  TF1 *fpt = new TF1("fpt","[0]+[1]*exp(-[2]*x)",30,230);
  fpt->SetParameters(0.958, 6.664, 0.1424); // chi2/NDF=276.3/8
  fpt->SetLineColor(kBlack);
  gavegen->Fit(fpt,"QRN");
  cout << Form("  fpt->SetParameters(%1.4g, %1.4g, %1.4g);"
  	       " // chi2/NDF=%1.1f/%d",
  	       fpt->GetParameter(0), fpt->GetParameter(1),
  	       fpt->GetParameter(2), 
  	       fpt->GetChisquare(), fpt->GetNDF()) << endl;
  
  // Correction for mjdata/mjmc
  TF1 *fmd = new TF1("fmd","[0]",30,230);
  fmd->SetParameter(0, 3.5);
  fmd->SetLineColor(kBlue);
  hm->Fit(fmd,"QRN");
  cout << Form("  fmd->SetParameter(0, %1.4g);"
  	       " // chi2/NDF=%1.1f/%d",
  	       fmd->GetParameter(0),
  	       fmd->GetChisquare(), fmd->GetNDF()) << endl;


  //TH1D *h1 = tdrHist("h1","#LTm_{jet,X}#GT / #LTm_{jet,Y/Z}#GT",
  //TH1D *h1 = tdrHist("h1","(#LTm_{jet}#GT/#LTp_{T,ave}#GT) /"
  //		     " (#LTm_{gen}#GT/#LTp_{T,gen}#GT)",
  TH1D *h1 = tdrHist("h1","#LTm_{jet}#GT / #LTp_{T,ave}#GT"
  		     " over gen",
		     0.90+1e-4,1.20-1e-4,"p_{T,ave} (GeV)",30,230);
  TH1D *h1d = tdrHist("h1","Data/MC-1 (%)",-2,+2,//+2,+5,
		      "p_{T,ave} (GeV)",30,230);

  //TCanvas *c1 = tdrCanvas("c1",h1,4,11,kSquare);
  if (mode=="17V5") lumi_13TeV = "2017, 43.9 fb^{-1}";
  if (mode=="18V5") lumi_13TeV = "2018, 59.9 fb^{-1}";
  TCanvas *c1 = tdrDiCanvas("c1",h1,h1d,4,11);
  c1->cd(1);
  gPad->SetLogx();

  tdrDraw(hmreco,"Pz",kFullSquare,kGray);
  tdrDraw(hmreco2,"Pz",kFullSquare,kGray+1);
  tdrDraw(hmmc2,"Pz",kFullSquare,kRed);
  tdrDraw(hmdata2,"Pz",kFullCircle,kBlue);
  //tdrDraw(havereco,"Pz",kOpenCircle,kBlack);
  //tdrDraw(havegen,"Pz",kOpenSquare,kBlack);
  tdrDraw(havegen,"Pz",kNone,kBlack);
  tdrDraw(gavegen,"Pz",kOpenSquare,kBlack);
  //gavegen->Draw("AP");
  fmj->Draw("SAME");
  fpt->Draw("SAME");

  
  TLine *l = new TLine();
  l->SetLineStyle(kDotted);
  l->DrawLine(30,1,230,1);
  l->SetLineColor(kGray+2);
  l->DrawLine(ptave,0.93,ptave,1.07);

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  tex->SetTextColor(kGray+1);
  tex->DrawLatex(0.41,0.06,"#LTp_{T,ave}#GT");

  TLegend *leg1 = tdrLeg(0.40,0.69,0.60,0.69+0.040*5);
  leg1->SetTextSize(0.035);
  leg1->AddEntry(hmdata2,"DATA","PLE");
  leg1->AddEntry(hmmc2,"MC","PLE");
  //leg1->AddEntry(hmreco2,"X=reco,Y=gen,Z=gen","PLE");
  //leg1->AddEntry(hmreco,"X=reco,Y=gen,Z=reco","PLE");
  //leg1->AddEntry(havegen,"gen/reco (p_{T,ave})","PLE");
  leg1->AddEntry(hmreco2,"Reco-gen match","PLE");
  //leg1->AddEntry(hmreco,"#LTp_{T,gen}#GT #rightarrow #LTp_{T,ave}#GT "
  //		 " (#LTm_{jet}#GT / #LTm_{gen}#GT)","PLE");
  leg1->AddEntry(hmreco,"#LTm_{jet}#GT / #LTm_{gen}#GT","PLE");
  //leg1->AddEntry(havegen,"#LTp_{T,gen}#GT / #LTp_{T,ave}#GT","PLE");
  leg1->AddEntry(gavegen,"#LTp_{T,gen}#GT / #LTp_{T,ave}#GT","PLE");


  c1->cd(2);
  gPad->SetLogx();

  tdrDraw(hm,"Pz",kFullCircle,kBlue);
  fmd->Draw("SAME");


  // Non-perturbatively
  //	m^2 = mu_NP*R*pT   (pT^0.5)
  // Perturbatively
  // (Eq. (6.1) of https://arxiv.org/abs/1901.10342
  //  and Eq. (4) of https://arxiv.org/abs/1704.05066)
  //	<m^2> = 1/2*R^2*pT^2 * alpha_s * CF / pi.   (pT^1)
  // Here multiplier is 0.045^2 for alphaS(MZ)=0.118, although I guess
  // Q^2 should be (mW/2)^2. But would we then always have pT~mW/2 as well?
  // UE effects should be pT^0, but some high exponent of R

  //TH1D *h2 = tdrHist("h2","#LTm_{jet,gen}#GT / #LTp_{T,ave,gen}#GT",
  TH1D *h2 = tdrHist("h2","#LTm_{jet}#GT / #LTp_{T,ave}#GT",
		     0.05,0.25,"p_{T,ave} (GeV)",30,230);
  TH1D *h2d = tdrHist("h2d","Data/MC-1 (%)",-2,+2,//+2,+5,
		      "p_{T,ave} (GeV)",30,230);
  
  //TCanvas *c2 = tdrCanvas("c2",h2,4,11,kSquare);
  TCanvas *c2 = tdrDiCanvas("c2",h2,h2d,11);
  c2->cd(1);
  gPad->SetLogx();

  tdrDraw(hmgen,"Pz",kFullSquare,kGreen+2);
  tdrDraw(hmrecorc,"Pz",kFullSquare,kGray+1);
  tdrDraw(hmrecomc,"Pz",kFullSquare,kRed);
  tdrDraw(hmrecodt,"Pz",kFullCircle,kBlue);

  TLegend *leg2 = tdrLeg(0.60,0.69,0.80,0.69+0.040*5);
  leg2->SetTextSize(0.035);
  leg2->AddEntry(hmrecodt,"DATA","PLE");
  leg2->AddEntry(hmrecomc,"MC","PLE");
  leg2->AddEntry(hmrecorc,"Reco-gen match","PLE");
  leg2->AddEntry(hmgen,"Gen","PLE");

  c2->cd(2);
  gPad->SetLogx();

  tdrDraw(hm,"Pz",kFullCircle,kBlue);

  c1->SaveAs(Form("pdf/drawWjetMass_RecoOverGen_MPDcorrNoW_%s.pdf",cm));
  c2->SaveAs(Form("pdf/drawWjetMass_MassOverPt_MPDcorrNoW_%s.pdf",cm));
  //c1->SaveAs("pdf/drawWjetMass_RecoOverGen_MPcorrNoW.pdf");
  //c2->SaveAs("pdf/drawWjetMass_MassOverPt_MPcorrNoW.pdf");
  //c1->SaveAs("pdf/drawWjetMass_RecoOverGen_McorrNoW.pdf");
  //c2->SaveAs("pdf/drawWjetMass_MassOverPt_McorrNoW.pdf");
} // drawWjetMass


void drawBjetMass(string mode) {

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  const char *cm = mode.c_str();

  TFile *fm = new TFile(Form("rootfiles/hadWMC%s_MPDcorrNoW.root",cm),"READ");
  assert(fm && !fm->IsZombie());

  TFile *fd = new TFile(Form("rootfiles/hadWUL%s_MPDcorrNoW.root",cm),"READ");
  assert(fd && !fd->IsZombie());

  curdir->cd();

  TProfile *pmbreco = (TProfile*)fm->Get("pmbreco"); assert(pmbreco);
  TProfile *pmbgen = (TProfile*)fm->Get("pmbgen");   assert(pmbgen);
  TProfile *pmbm = (TProfile*)fm->Get("pmbjet");     assert(pmbm);
  TProfile *pmbd = (TProfile*)fd->Get("pmbjet");     assert(pmbd);

  TProfile *pbreco = (TProfile*)fm->Get("pbreco"); assert(pbreco);
  TProfile *pbgen = (TProfile*)fm->Get("pbgen");   assert(pbgen);

  double ptb = pbreco->GetMean();

  TH1D *hmbreco = pmbreco->ProjectionX("hmbreco"); // <mbjet>/<ptb>
  hmbreco->Divide(pmbgen);  // divide by <mbgen>/<ptb> => <mbjet>/<mbgen>
  TH1D *hmbmc = pmbm->ProjectionX("hmbmc"); // <mbjet>/<ptb>
  hmbmc->Divide(pmbgen); // divide by <mbgen>/<ptb> => <mbjet>/<mbgen>
  TH1D *hmbdata = pmbd->ProjectionX("hmbdata");
  hmbdata->Divide(pmbgen);

  TH1D *hbgen = pbgen->ProjectionX("hbgen");
  TGraphErrors *gbgen = new TGraphErrors(0);
  for (int i = 1; i != hbgen->GetNbinsX()+1; ++i) {
    double xmid = hbgen->GetBinCenter(i);
    if (xmid<30 || xmid>230) continue;
    double x = pbreco->GetBinContent(pbreco->FindBin(xmid));
    double ex = pbreco->GetBinError(pbreco->FindBin(xmid));
    int n = gbgen->GetN();
    gbgen->SetPoint(n, x, hbgen->GetBinContent(i));
    gbgen->SetPointError(n, ex, hbgen->GetBinError(i));
  }

  TH1D *hmbreco2 = (TH1D*)hmbreco->Clone("hmbreco2");
  hmbreco2->Multiply(hbgen);
  TH1D *hmbmc2 = (TH1D*)hmbmc->Clone("hmbmc2");
  hmbmc2->Multiply(hbgen); // divide by <ptbgen>/<ptb>
                          //  => (<mbjet>/<ptb>) / (<mbgen>/<ptbgen>)
  TH1D *hmbdata2 = (TH1D*)hmbdata->Clone("hmbdata2");
  hmbdata2->Multiply(hbgen);

  TH1D *hmbgen = pmbgen->ProjectionX("hmbgen");
  hmbgen->Divide(pbgen);
  TH1D *hmbrecorc = pmbreco->ProjectionX("hmbrecorc");
  //hmbreco->Divide(pbreco); // already divided
  TH1D *hmbrecomc = pmbm->ProjectionX("hmbrecomc");  
  TH1D *hmbrecodt = pmbd->ProjectionX("hmbrecodt");  

  TH1D *hmb = (TH1D*)hmbdata2->Clone("hmb");
  hmb->Add(hmbmc2,-1); // are stat. unc. correct with add+divide?
  hmb->Divide(hmbmc2);
  hmb->Scale(100.);

  TH1D *hmb2 = (TH1D*)hmbrecodt->Clone("hmb2");
  hmb2->Add(hmbrecomc,-1); // are stat. unc. correct with add+divide?
  hmb2->Divide(hmbrecomc);
  hmb2->Scale(100.);

  // Correction for reco mass/pt ratio
  TF1 *fmbj = new TF1("fmbj","1+[1]*log(x/[0])",30,230);
  fmbj->SetParameters(ptb,0.05);
  fmbj->SetLineColor(kRed);
  hmbmc2->Fit(fmbj,"QRN");
  //hmbreco2->Fit(fmbj,"QRN");
  cout << Form("  fmbj->SetParameters(%1.3g, %1.4g); // ptb = %1.3g GeV",
	       fmbj->GetParameter(0), fmbj->GetParameter(1), ptb) << endl;

  // Correction for pTgen/pTreco ratio
  TF1 *fbpt = new TF1("fbpt","[0]+[1]*exp(-[2]*x)",30,230);
  fbpt->SetParameters(0.958, 6.664, 0.1424); // chi2/NDF=276.3/8
  fbpt->SetLineColor(kBlack);
  gbgen->Fit(fbpt,"QRN");
  cout << Form("  fbpt->SetParameters(%1.4g, %1.4g, %1.4g);"
  	       " // chi2/NDF=%1.1f/%d",
  	       fbpt->GetParameter(0), fbpt->GetParameter(1),
  	       fbpt->GetParameter(2), 
  	       fbpt->GetChisquare(), fbpt->GetNDF()) << endl;
  
  // Correction for mjdata/mjmc
  TF1 *fmbd = new TF1("fmbd","[0]",30,230);
  fmbd->SetParameter(0, 3.5);
  fmbd->SetLineColor(kBlue);
  hmb->Fit(fmbd,"QRN");
  cout << Form("  fmbd->SetParameter(0, %1.4g);"
  	       " // chi2/NDF=%1.1f/%d",
  	       fmbd->GetParameter(0),
  	       fmbd->GetChisquare(), fmbd->GetNDF()) << endl;


  TH1D *hb1 = tdrHist("hb1","#LTm_{b-jet}#GT / #LTp_{T,b}#GT"
		      " over gen",
		     0.90+1e-4,1.20-1e-4,"p_{T,b} (GeV)",30,230);
  TH1D *hb1d = tdrHist("hb1","Data/MC-1 (%)",
		       +2,5,"p_{T,b} (GeV)",30,230);

  TCanvas *cb1 = tdrDiCanvas("cb1",hb1,hb1d,4,11);
  cb1->cd(1);
  gPad->SetLogx();

  tdrDraw(hmbreco,"Pz",kFullSquare,kGray);
  tdrDraw(hmbreco2,"Pz",kFullSquare,kGray+1);
  tdrDraw(hmbmc2,"Pz",kFullSquare,kRed);
  tdrDraw(hmbdata2,"Pz",kFullCircle,kBlue);
  //tdrDraw(hbreco,"Pz",kOpenCircle,kBlack);
  //tdrDraw(hbgen,"Pz",kOpenSquare,kBlack);
  tdrDraw(hbgen,"Pz",kNone,kBlack);
  tdrDraw(gbgen,"Pz",kOpenSquare,kBlack);
  //gbgen->Draw("AP");
  fmbj->Draw("SAME");
  fbpt->Draw("SAME");

  
  TLine *l = new TLine();
  l->SetLineStyle(kDotted);
  l->DrawLine(30,1,230,1);
  l->SetLineColor(kGray+2);
  l->DrawLine(ptb,0.93,ptb,1.07);

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  tex->SetTextColor(kGray+1);
  tex->DrawLatex(0.41,0.06,"#LTp_{T,b}#GT");

  TLegend *leg1 = tdrLeg(0.40,0.69,0.60,0.69+0.040*5);
  leg1->SetTextSize(0.035);
  leg1->AddEntry(hmbdata2,"DATA","PLE");
  leg1->AddEntry(hmbmc2,"MC","PLE");
  //leg1->AddEntry(hmbreco2,"X=reco,Y=gen,Z=gen","PLE");
  //leg1->AddEntry(hmbreco,"X=reco,Y=gen,Z=reco","PLE");
  //leg1->AddEntry(hbgen,"gen/reco (p_{T,b})","PLE");
  leg1->AddEntry(hmbreco2,"Reco-gen match","PLE");
  //leg1->AddEntry(hmbreco,"#LTp_{T,gen}#GT #rightarrow #LTp_{T,b}#GT "
  //		 " (#LTm_{bjet}#GT / #LTm_{gen}#GT)","PLE");
  leg1->AddEntry(hmbreco,"#LTm_{bjet}#GT / #LTm_{bgen}#GT","PLE");
  //leg1->AddEntry(hbgen,"#LTp_{T,bgen}#GT / #LTp_{T,b}#GT","PLE");
  leg1->AddEntry(gbgen,"#LTp_{T,bgen}#GT / #LTp_{T,b}#GT","PLE");


  cb1->cd(2);
  gPad->SetLogx();

  tdrDraw(hmb,"Pz",kFullCircle,kBlue);
  fmbd->Draw("SAME");

  TH1D *hb2 = tdrHist("hb2","#LTm_{bjet}#GT / #LTp_{T,b}#GT",
		      0.05,0.25,"p_{T,b} (GeV)",30,230);
  TH1D *hb2d = tdrHist("h2d","Data/MC-1 (%)",
		     +2,+5,"p_{T,b} (GeV)",30,230);
  
  TCanvas *cb2 = tdrDiCanvas("cb2",hb2,hb2d,11);
  cb2->cd(1);
  gPad->SetLogx();

  tdrDraw(hmbgen,"Pz",kFullSquare,kGreen+2);
  tdrDraw(hmbrecorc,"Pz",kFullSquare,kGray+1);
  tdrDraw(hmbrecomc,"Pz",kFullSquare,kRed);
  tdrDraw(hmbrecodt,"Pz",kFullCircle,kBlue);

  TLegend *leg2 = tdrLeg(0.60,0.69,0.80,0.69+0.040*5);
  leg2->SetTextSize(0.035);
  leg2->AddEntry(hmbrecodt,"DATA","PLE");
  leg2->AddEntry(hmbrecomc,"MC","PLE");
  leg2->AddEntry(hmbrecorc,"Reco-gen match","PLE");
  leg2->AddEntry(hmbgen,"Gen","PLE");

  cb2->cd(2);
  gPad->SetLogx();

  tdrDraw(hmb,"Pz",kFullCircle,kBlue);

  cb1->SaveAs(Form("pdf/drawWjetMass_BRecoOverBGen_MPDGcorrNoW_%s.pdf",cm));
  cb2->SaveAs(Form("pdf/drawWjetMass_BMassOverBPt_MPDGcorrNoW_%s.pdf",cm));
} // drawBjetmass
