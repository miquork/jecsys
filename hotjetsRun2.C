#include "TFile.h"
#include "TH2D.h"
#include "TMath.h"
#include "TBox.h"

#include "tdrstyle_mod15.C"

bool excludeHEP17 = false;

// Set negative values to zero in TH2D for better "BOX" drawing
void rezero(TH2D* h2, double thr=0, double max=10) {
  
  assert(h2);

  for (int i = 1; i != h2->GetNbinsX()+1; ++i) {
    for (int j = 1; j != h2->GetNbinsY()+1; ++j) {
      if (h2->GetBinContent(i,j)<thr) h2->SetBinContent(i,j,0);
      if (h2->GetBinContent(i,j)>max) h2->SetBinContent(i,j,max);
    }
  }
} // void rezero


void hotjetsRun2() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TFile *f16 = new TFile("rootfiles/hotjets-UL16.root","READ");
  assert(f16 && !f16->IsZombie());

  TFile *f17 = new TFile("rootfiles/hotjets-UL17_V2.root","READ");
  assert(f17 && !f17->IsZombie());

  TFile *f18 = new TFile("rootfiles/hotjets-UL18.root","READ");
  assert(f18 && !f18->IsZombie());


  curdir->cd();

  TH2D *h216 = (TH2D*)f16->Get("h2hot_ul16"); assert(h216);
  TH2D *h217 = (TH2D*)f17->Get("h2hot_ul17"); assert(h217);
  TH2D *h217b = (TH2D*)f17->Get("h2hot_ul17_plus_hep17_plus_hbpw89"); assert(h217b);
  TH2D *h218 = (TH2D*)f18->Get("h2hot_ul18"); assert(h218);
  TH2D *h218b = (TH2D*)f18->Get("h2hot_ul18_plus_hem1516_and_h2bp2m1"); assert(h218);
  //TH2D *h2all = (TH2D*)f2->Get("h2hotfilter"); assert(h2all);
  TH2D *h2em = (TH2D*)fem->Get("all/h2hole"); assert(h2em);
  TH2D *h2mc = (TH2D*)fmc->Get("h2hotfilter"); assert(h2mc);

  TH1D *h = new TH1D("h",";#eta_{jet};#phi_{jet}",100,-4.7,4.7);
  h->SetMaximum(+TMath::Pi());
  h->SetMinimum(-TMath::Pi());

  // 36.5+43.9+59.9
  lumi_13TeV = "Run2, 140.3 fb^{-1}";
  TCanvas *c1 = tdrCanvas("c1",h,4,0,kRectangular);

  TLine *l = new TLine();

  l->SetLineStyle(kSolid);
  double etahf = 2.964;
  l->DrawLine(-etahf,-TMath::Pi(),-etahf,+TMath::Pi());
  l->DrawLine(+etahf,-TMath::Pi(),+etahf,+TMath::Pi());
  l->SetLineStyle(kDashed);
  double etatr = 2.5;
  l->DrawLine(-etatr,-TMath::Pi(),-etatr,+TMath::Pi());
  l->DrawLine(+etatr,-TMath::Pi(),+etatr,+TMath::Pi());
  l->SetLineStyle(kDotted);
  double etaec = 1.305;
  l->DrawLine(-etaec,-TMath::Pi(),-etaec,+TMath::Pi());
  l->DrawLine(+etaec,-TMath::Pi(),+etaec,+TMath::Pi());

  rezero(h2mc);
  h2mc->GetZaxis()->SetRangeUser(-10,10);
  h2mc->SetLineColor(kMagenta);
  h2mc->SetFillStyle(1001);
  h2mc->SetFillColor(kNone);
  h2mc->SetFillColorAlpha(kMagenta-9, 0.35); // 35% transparent
  h2mc->Draw("SAMEBOX");

  rezero(h2all);
  h2all->GetZaxis()->SetRangeUser(-10,10);
  h2all->SetLineColor(kRed);
  h2all->SetLineStyle(kNone);
  h2all->SetFillStyle(1001);
  //h2all->DrawClone("SAMEBOX");
  h2all->SetFillColor(kRed);
  h2all->SetFillColorAlpha(kRed, 0.35); // 35% transparent
  //h2all->DrawClone("SAMEBOX");

  rezero(h2bcd);
  h2bcd->GetZaxis()->SetRangeUser(-10,10);
  h2bcd->SetFillStyle(1001);
  h2bcd->SetLineColor(kYellow+1);
  h2bcd->SetFillColorAlpha(kYellow, 0.35); // 35% transparent
  h2bcd->Draw("SAMEBOX");

  rezero(h2ef);
  h2ef->GetZaxis()->SetRangeUser(-10,10);
  h2ef->SetFillStyle(1001);
  h2ef->SetLineColor(kOrange+1);
  h2ef->SetFillColorAlpha(kOrange, 0.35); // 35% transparent
  h2ef->Draw("SAMEBOX");

  rezero(h2gh);
  h2gh->GetZaxis()->SetRangeUser(-10,10);
  h2gh->SetFillStyle(1001);
  h2gh->SetLineColor(kRed+1);
  h2gh->SetFillColorAlpha(kRed, 0.35); // 35% transparent
  h2gh->Draw("SAMEBOX");

  //h2all->SetLineColor(kRed);
  //h2all->SetLineStyle(kNone);
  h2all->SetFillStyle(1001);
  //h2all->SetFillColor(kNone);
  h2all->SetFillColor(kRed);
  h2all->SetFillColorAlpha(kRed, 0.35); // 35% transparent
  //h2all->DrawClone("SAMEBOX");

  rezero(h2em);
  h2em->GetZaxis()->SetRangeUser(-10,10);
  h2em->SetLineColor(kBlue);
  h2em->SetFillStyle(1001);
  h2em->SetFillColor(kNone);
  //h2em->DrawClone("SAMEBOX");
  //h2em->SetFillColorAlpha(kBlue, 0.35); // 35% transparent
  //h2em->SetFillColorAlpha(kAzure-9, 0.35); // 35% transparent
  //h2em->Draw("SAMEBOX");

  // combination of regions
  TH2D *h2sum = (TH2D*)h2bcd->Clone("h2hot_ul16");
  h2sum->Add(h2bcd);
  h2sum->Add(h2ef);
  h2sum->Add(h2gh);
  //h2sum->Add(h2all);
  rezero(h2sum,20,10); // overlap min. 2
  //rezero(h2sum); // no overlap needed

  /*
  // Remove also HEP17
  TBox HEP17(1.31,-0.5236,2.96,-0.8727); // centered at 28*dphi+/-2
//TBox HBPw89(0,2.793,1.4835,3.1416); // centered at 8.5*4*dphi+/-2, wide barrel
  TBox HBPw89(0,2.705,1.4835,3.1416); // centered at 8.5*4*dphi+/-2, v2
  TH2D *h2hep17 = (TH2D*)h2sum->Clone("h2hot_ul17_plus_hep17");
  TH2D *h2hbpw89 = (TH2D*)h2sum->Clone("h2hot_ul17_plus_hbpw89");
  TH2D *h2both = (TH2D*)h2sum->Clone("h2hot_ul17_plus_hep17_plus_hbpw89");
  for (int i = 1; i != h2sum->GetNbinsX()+1; ++i) {
    for (int j = 1; j != h2sum->GetNbinsY()+1; ++j) {
      double eta = h2sum->GetXaxis()->GetBinCenter(i);
      double phi = h2sum->GetYaxis()->GetBinCenter(j);
      if (eta>HEP17.GetX1() && eta<HEP17.GetX2() &&
	  phi>HEP17.GetY1() && phi<HEP17.GetY2())
	h2hep17->SetBinContent(i, j, 10);
      if (eta>HBPw89.GetX1() && eta<HBPw89.GetX2() &&
	  phi>HBPw89.GetY1() && phi<HBPw89.GetY2())
	h2hbpw89->SetBinContent(i, j, 10);
      if ((eta>HEP17.GetX1() && eta<HEP17.GetX2() &&
	   phi>HEP17.GetY1() && phi<HEP17.GetY2()) ||
	  (eta>HBPw89.GetX1() && eta<HBPw89.GetX2() &&
	   phi>HBPw89.GetY1() && phi<HBPw89.GetY2()))
	h2both->SetBinContent(i, j, 10);
    } // for i
  } // for j

  h2both->SetLineColor(kMagenta-9);
  h2both->SetLineStyle(kNone);
  h2both->SetFillStyle(1001);
  h2both->SetFillColor(kNone);
  h2both->DrawClone("SAMEBOX");

  h2hbpw89->SetLineColor(kBlue-9);
  h2hbpw89->SetLineStyle(kNone);
  h2hbpw89->SetFillStyle(1001);
  h2hbpw89->SetFillColor(kNone);
  h2hbpw89->DrawClone("SAMEBOX");

  h2hep17->SetLineColor(kRed-9);
  h2hep17->SetLineStyle(kNone);
  h2hep17->SetFillStyle(1001);
  h2hep17->SetFillColor(kNone);
  h2hep17->DrawClone("SAMEBOX");
  */

  //h2em->Draw("SAMEBOX");
  h2em->DrawClone("SAMEBOX");
  //h2em->SetFillColorAlpha(kBlue, 0.35); // 35% transparent
  h2em->SetFillColorAlpha(kAzure-9, 0.35); // 35% transparent
  h2em->Draw("SAMEBOX");

  h2sum->SetLineColor(kBlack);
  h2sum->SetLineStyle(kNone);
  h2sum->SetFillStyle(1001);
  h2sum->SetFillColor(kNone);
  h2sum->DrawClone("SAMEBOX");

  //TLegend *leg = tdrLeg(0.15,0.45,0.35,0.90);
  TLegend *leg = tdrLeg(0.15,0.90-6*0.05,0.35,0.90);
  leg->AddEntry(h2em,"Cold (GH)","F");
  leg->AddEntry(h2bcd,"BCD","F");
  leg->AddEntry(h2ef,"EF","F");
  leg->AddEntry(h2gh,"GH","F");
  //leg->AddEntry(h2all,"BCDEFGH","F");
  //leg->AddEntry(h2em,"EM mask","F");
  leg->AddEntry(h2sum,"UL (min. 2)","F");
  //leg->AddEntry(h2sum,"UL (min. 1)","F");
  leg->AddEntry(h2mc,"MC hot","F");

  // Count fraction of towers in the veto map
  /*
  double sum_bb(0), sum_ec1(0), sum_ec2(0), sum_hf(0);
  double em_bb(0), em_ec1(0), em_ec2(0), em_hf(0);
  double hot_bb(0), hot_ec1(0), hot_ec2(0), hot_hf(0);
  for (int ieta = 1; ieta != h2em->GetNbinsX()+1; ++ieta) {

    double eta = h2em->GetXaxis()->GetBinCenter(ieta);
    for (int iphi = 1; iphi != h2em->GetNbinsY()+1; ++iphi) {

      if (fabs(eta)<etaec) {
	++sum_bb;
	if (h2em->GetBinContent(ieta,iphi)>0) ++em_bb;
	if (h2all->GetBinContent(ieta,iphi)>0) ++hot_bb;
      }
      if (fabs(eta)>etaec && fabs(eta)<etatr) {
	++sum_ec1;
	if (h2em->GetBinContent(ieta,iphi)>0) ++em_ec1;
	if (h2all->GetBinContent(ieta,iphi)>0) ++hot_ec1;
      }
      if (fabs(eta)>etatr && fabs(eta)<etahf) {
	++sum_ec2;
	if (h2em->GetBinContent(ieta,iphi)>0) ++em_ec2;
	if (h2all->GetBinContent(ieta,iphi)>0) ++hot_ec2;
      }
      if (fabs(eta)>etahf) {
	++sum_hf;
	if (h2em->GetBinContent(ieta,iphi)>0) ++em_hf;
	if (h2all->GetBinContent(ieta,iphi)>0) ++hot_hf;
      }
      
    } // for iphi
  } // for ieta

  // Draw fraction of towers in veto map
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.030);

  // For EM masked towers (cold jets)
  tex->SetTextColor(kBlue);
  tex->DrawLatex(0.439,0.656,Form("%1.1f%% (%1.1f / %1.1f / %1.1f)",
				100.*(em_bb+em_ec1+em_ec2) /
				  (sum_bb+sum_ec1+sum_ec2),
				  100.*em_bb/sum_bb, 100.*em_ec1/sum_ec1,
				  100.*em_ec2/sum_ec2));

  // For EM IC issues and other excesses (hot jets)
  tex->SetTextColor(kRed);
  tex->DrawLatex(0.668,0.817,Form("%1.1f%% (%1.1f / %1.1f / %1.1f / %1.1f)",
				  100.*(hot_bb+hot_ec1+hot_ec2+hot_hf) /
				  (sum_bb+sum_ec1+sum_ec2+sum_hf),
				  100.*hot_bb/sum_bb, 100.*hot_ec1/sum_ec1,
				  100.*hot_ec2/sum_ec2, 100.*hot_hf/sum_hf));
  */
  /*
  // Add HEP17
  // https://indico.cern.ch/event/619421/contributions/2500270/attachments/1424160/2184465/HE-plan1-xpog-feedback.pdf (slide 9 for coordinates)
  // https://cds.cern.ch/record/357153/files/CMS_HCAL_TDR.pdf (18 wedges, 4dphi)
  double dphi = TMath::TwoPi()/72.;
  //TBox HEP17(1.31,-0.7,2.96,-0.9); // guess
  //TBox HEP17(1.31,-0.32,2.96,-1.07);  // slide 9 above
  //TBox HEP17(1.31,-0.3491,2.96,-1.047); // adjusted to -pi+N*dphi edges
  //TBox HEP17(1.31,-0.5236,2.96,-0.8727); // centered at 28*dphi+/-2
  HEP17.SetFillStyle(kNone);
  HEP17.SetLineColor(kRed+2);
  HEP17.Draw("SAME");

  HBPw89.SetFillStyle(kNone);
  HBPw89.SetLineColor(kBlue+1);
  HBPw89.Draw("SAME");

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  tex->SetTextColor(kRed+2);
  tex->DrawLatex(0.66,0.46,"HEP17");
  tex->SetTextColor(kBlue+1);
  tex->DrawLatex(0.51,0.83,"HBPw8/9: BPIX?");
  */
  gPad->Paint();

  c1->SaveAs("pdf/hotjets2016UL.pdf");

  TFile *fout = new TFile("rootfiles/hotjets-UL16.root","RECREATE");
  h2sum->Write();
  //h2hep17->Write();
  //h2hbpw89->Write();
  //h2both->Write();
  fout->Close();
}
