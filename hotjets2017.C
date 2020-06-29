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


void hotjets2017() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TFile *fem = new TFile("../jecsys2017/rootfiles/coldjets-17runBCDEF.root","READ");
  assert(fem && !fem->IsZombie());

  TFile *f0 = new TFile("../jecsys2017/rootfiles/hotjets-17runBCDEF.root","READ");
  assert(f0 && !f0->IsZombie());


  TFile *f1 = new TFile("rootfiles/hotjets-17runB.root","READ");
  assert(f1 && !f1->IsZombie());

  TFile *f2 = new TFile("rootfiles/hotjets-17runC.root","READ");
  assert(f2 && !f2->IsZombie());

  TFile *f3 = new TFile("rootfiles/hotjets-17runD.root","READ");
  assert(f3 && !f3->IsZombie());

  TFile *f4 = new TFile("rootfiles/hotjets-17runE.root","READ");
  assert(f4 && !f4->IsZombie());

  TFile *f5 = new TFile("rootfiles/hotjets-17runF.root","READ");
  assert(f5 && !f5->IsZombie());


  curdir->cd();

  TH2D *h2b = (TH2D*)f1->Get("h2hotfilter"); assert(h2b);
  TH2D *h2c = (TH2D*)f2->Get("h2hotfilter"); assert(h2c);
  TH2D *h2d = (TH2D*)f3->Get("h2hotfilter"); assert(h2d);
  TH2D *h2e = (TH2D*)f4->Get("h2hotfilter"); assert(h2e);
  TH2D *h2f = (TH2D*)f5->Get("h2hotfilter"); assert(h2f);
  TH2D *h2all = (TH2D*)f0->Get("h2hotfilter"); assert(h2all);
  TH2D *h2em = (TH2D*)fem->Get("h2hole"); assert(h2em);

  TH1D *h = new TH1D("h",";#eta_{jet};#phi_{jet}",100,-4.7,4.7);
  h->SetMaximum(+TMath::Pi());
  h->SetMinimum(-TMath::Pi());

  lumi_13TeV = "2017 UL, 43.9 fb^{-1}";
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


  rezero(h2all);
  h2all->GetZaxis()->SetRangeUser(-10,10);
  h2all->SetLineColor(kRed);
  h2all->SetLineStyle(kNone);
  h2all->SetFillStyle(1001);
  //h2all->DrawClone("SAMEBOX");
  h2all->SetFillColor(kRed);
  h2all->SetFillColorAlpha(kRed, 0.35); // 35% transparent
  //h2all->DrawClone("SAMEBOX");

  rezero(h2b);
  h2b->GetZaxis()->SetRangeUser(-10,10);
  h2b->SetFillStyle(1001);
  h2b->SetLineColor(kBlue);
  h2b->SetFillColorAlpha(kBlue, 0.35); // 35% transparent
  h2b->Draw("SAMEBOX");

  rezero(h2c);
  h2c->GetZaxis()->SetRangeUser(-10,10);
  h2c->SetFillStyle(1001);
  h2c->SetLineColor(kCyan+1);
  h2c->SetFillColorAlpha(kCyan, 0.35); // 35% transparent
  h2c->Draw("SAMEBOX");

  rezero(h2d);
  h2d->GetZaxis()->SetRangeUser(-10,10);
  h2d->SetFillStyle(1001);
  h2d->SetLineColor(kGreen+1);
  h2d->SetFillColorAlpha(kGreen, 0.35); // 35% transparent
  h2d->Draw("SAMEBOX");

  rezero(h2e);
  h2e->GetZaxis()->SetRangeUser(-10,10);
  h2e->SetFillStyle(1001);
  h2e->SetLineColor(kYellow+1);
  h2e->SetFillColorAlpha(kYellow, 0.35); // 35% transparent
  h2e->Draw("SAMEBOX");

  rezero(h2f);
  h2f->GetZaxis()->SetRangeUser(-10,10);
  h2f->SetFillStyle(1001);
  h2f->SetLineColor(kOrange);
  h2f->SetFillColorAlpha(kOrange, 0.35); // 35% transparent
  h2f->Draw("SAMEBOX");

  h2all->SetLineColor(kGray+1);//kRed);
  h2all->SetLineStyle(kNone);
  h2all->SetFillStyle(1001);
  h2all->SetFillColor(kNone);
  h2all->DrawClone("SAMEBOX");

  rezero(h2em);
  h2em->GetZaxis()->SetRangeUser(-10,10);
  h2em->SetLineColor(kBlue);
  h2em->SetFillStyle(1001);
  h2em->SetFillColor(kNone);
  //h2em->SetFillColorAlpha(kBlue, 0.35); // 35% transparent
  //h2em->Draw("SAMEBOX");

  // combination of regions
  TH2D *h2sum = (TH2D*)h2b->Clone("h2hot_ul17");
  h2sum->Add(h2c);
  h2sum->Add(h2d);
  h2sum->Add(h2e);
  h2sum->Add(h2f);
  rezero(h2sum,20,10); // overlap min. 2
  //rezero(h2sum); // no overlap needed

  // Remove also HEP17
  TBox HEP17(1.31,-0.5236,2.96,-0.8727); // centered at 28*dphi+/-2
  TH2D *h2hep17 = (TH2D*)h2sum->Clone("h2hot_ul17_plus_hep17");
  for (int i = 1; i != h2sum->GetNbinsX()+1; ++i) {
    for (int j = 1; j != h2sum->GetNbinsY()+1; ++j) {
      double eta = h2sum->GetXaxis()->GetBinCenter(i);
      double phi = h2sum->GetYaxis()->GetBinCenter(j);
      if (eta>HEP17.GetX1() && eta<HEP17.GetX2() &&
	  phi>HEP17.GetY1() && phi<HEP17.GetY2())
	h2hep17->SetBinContent(i, j, 10);
    } // for i
  } // for j

  h2hep17->SetLineColor(kRed-9);
  h2hep17->SetLineStyle(kNone);
  h2hep17->SetFillStyle(1001);
  h2hep17->SetFillColor(kNone);
  h2hep17->DrawClone("SAMEBOX");

  h2sum->SetLineColor(kBlack);
  h2sum->SetLineStyle(kNone);
  h2sum->SetFillStyle(1001);
  h2sum->SetFillColor(kNone);
  h2sum->DrawClone("SAMEBOX");

  //TLegend *leg = tdrLeg(0.43,0.64,0.63,0.84);
  TLegend *leg = tdrLeg(0.43,0.60,0.63,0.90);
  //leg->AddEntry(h2all,"BCDEF","F");
  leg->AddEntry(h2all,"non-UL","F");
  leg->AddEntry(h2b,"B","F");
  leg->AddEntry(h2c,"C","F");
  leg->AddEntry(h2d,"D","F");
  leg->AddEntry(h2e,"E","F");
  leg->AddEntry(h2f,"F","F");
  //leg->AddEntry(h2em,"EM mask","F");
  leg->AddEntry(h2sum,"UL (min. 2)","F");
  //leg->AddEntry(h2sum,"UL (min. 1)","F");

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

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  tex->SetTextColor(kRed+2);
  tex->DrawLatex(0.66,0.46,"HEP17");
  
  gPad->Paint();

  c1->SaveAs("pdf/hotjets2017UL.pdf");

  TFile *fout = new TFile("rootfiles/hotjets-UL17.root","RECREATE");
  h2sum->Write();
  h2hep17->Write();
  fout->Close();
}
