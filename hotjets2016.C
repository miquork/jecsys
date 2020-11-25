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


void hotjets2016() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  //TFile *fem = new TFile("rootfiles/coldjets-UL16runGH.root","READ");
  TFile *fem = new TFile("rootfiles/coldjets-UL16runBCDEFGH.root","READ");
  assert(fem && !fem->IsZombie());

  TFile *f0 = new TFile("../jecsys2016/rootfiles/hotjets-runBCDEFGH.root","READ");
  assert(f0 && !f0->IsZombie());

  TFile *fmc = new TFile("rootfiles/hotjetsmc-UL16runGH.root","READ");
  assert(fmc && !fmc->IsZombie());


  TFile *f1 = new TFile("rootfiles/hotjets-UL16runBCD.root","READ");
  assert(f1 && !f1->IsZombie());

  TFile *f2 = new TFile("rootfiles/hotjets-UL16runEF.root","READ");
  assert(f2 && !f2->IsZombie());

  TFile *f3 = new TFile("rootfiles/hotjets-UL16runGH.root","READ");
  assert(f3 && !f3->IsZombie());

  TFile *f4 = new TFile("rootfiles/hotjets-UL16runBCDEFGH.root","READ");
  //TFile *f4 = new TFile("rootfiles/hotjets-UL16runGH.root","READ");
  assert(f4 && !f4->IsZombie());

  TFile *fecal = new TFile("rootfiles/UL2018_badchannels.root","READ");
  assert(fecal && !fecal->IsZombie());


  curdir->cd();

  TH2D *h2bcd = (TH2D*)f1->Get("h2hotfilter"); assert(h2bcd);
  TH2D *h2ef = (TH2D*)f2->Get("h2hotfilter"); assert(h2ef);
  TH2D *h2gh = (TH2D*)f3->Get("h2hotfilter"); assert(h2gh);
  TH2D *h2all = (TH2D*)f4->Get("h2hotfilter"); assert(h2all);
  //TH2D *h2old = (TH2D*)f0->Get("h2hotfilter"); assert(h2old);
  TH2D *h2old = (TH2D*)f0->Get("h2jet"); assert(h2old);
  //TH2D *h2em = (TH2D*)fem->Get("h2hole"); assert(h2em);
  TH2D *h2em = (TH2D*)fem->Get("all/h2hole"); assert(h2em);
  TH2D *h2mc = (TH2D*)fmc->Get("h2hotfilter"); assert(h2mc);

  TH2D *h2ecal = (TH2D*)fecal->Get("ecalmap"); assert(h2ecal);

  TH1D *h = new TH1D("h",";#eta_{jet};#phi_{jet}",100,-4.7,4.7);
  h->SetMaximum(+TMath::Pi());
  h->SetMinimum(-TMath::Pi());

  lumi_13TeV = "2016 UL, 36.5 fb^{-1}";
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
  h2mc->SetLineColor(kNone);//kMagenta);
  h2mc->SetLineStyle(kNone);
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
  h2gh->SetLineColor(kOrange+2);//kRed+1);
  //h2gh->SetFillColorAlpha(kRed, 0.35); // 35% transparent
  h2gh->SetFillColorAlpha(kOrange+1, 0.35); // 35% transparent
  h2gh->Draw("SAMEBOX");

  //h2all->SetLineColor(kRed);
  //h2all->SetLineStyle(kNone);
  h2all->SetFillStyle(1001);
  //h2all->SetFillColor(kNone);
  h2all->SetFillColor(kRed);
  h2all->SetFillColorAlpha(kRed, 0.35); // 35% transparent
  h2all->DrawClone("SAMEBOX");

  rezero(h2em);
  h2em->GetZaxis()->SetRangeUser(-10,10);
  h2em->SetLineColor(kBlue);
  h2em->SetFillStyle(1001);
  h2em->SetFillColor(kNone);
  //h2em->DrawClone("SAMEBOX");
  //h2em->SetFillColorAlpha(kBlue, 0.35); // 35% transparent
  //h2em->SetFillColorAlpha(kAzure-9, 0.35); // 35% transparent
  //h2em->Draw("SAMEBOX");

  rezero(h2ecal);
  h2ecal->GetZaxis()->SetRangeUser(-10,10);
  h2ecal->SetLineColor(kNone);//kBlue);
  h2ecal->SetLineStyle(kNone);
  h2ecal->SetFillStyle(1001);
  h2ecal->SetFillColor(kNone);

  // combination of regions
  TH2D *h2sum = (TH2D*)h2bcd->Clone("h2hot_ul16");
  h2sum->Add(h2bcd);
  h2sum->Add(h2ef);
  h2sum->Add(h2gh);
  h2sum->Add(h2all);
  h2sum->Add(h2all); // h2all or any two others
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

  //h2ecal->DrawClone("SAMEBOX");
  h2ecal->SetFillColorAlpha(kBlue,0.50);//kAzure-9, 0.35); // 35% transparent
  h2ecal->Draw("SAMEBOX");

  h2sum->SetLineColor(kBlack);
  h2sum->SetLineStyle(kNone);
  h2sum->SetFillStyle(1001);
  h2sum->SetFillColor(kNone);
  h2sum->DrawClone("SAMEBOX");

  for (int i = 1; i != h2old->GetNbinsX()+1; ++i) {
    for (int j = 1; j != h2old->GetNbinsY()+1; ++j) {
      if (h2old->GetBinContent(i,j)<0) h2old->SetBinContent(i,j,0);
    } // for j
  } // for i
  h2old->SetLineColor(kGray+1);
  h2old->SetLineStyle(kDotted);
  h2old->SetFillStyle(1001);
  h2old->SetFillColor(kNone);
  //h2old->DrawClone("SAMEBOX");

  //TLegend *leg = tdrLeg(0.15,0.45,0.35,0.90);
  //TLegend *leg = tdrLeg(0.15,0.90-7*0.05,0.35,0.90);
  TLegend *leg = tdrLeg(0.15,0.90-8*0.05,0.35,0.90);
  //leg->AddEntry(h2em,"Cold (GH)","F");
  leg->AddEntry(h2bcd,"BCD","F");
  leg->AddEntry(h2ef,"EF","F");
  leg->AddEntry(h2gh,"GH","F");
  leg->AddEntry(h2all,"B-H","F");
  //leg->AddEntry(h2em,"EM mask","F");
  leg->AddEntry(h2sum,"Hot (min. 2)","F");
  //leg->AddEntry(h2sum,"UL (min. 1)","F");
  leg->AddEntry(h2mc,"MC hot","F");
  leg->AddEntry(h2em,"Cold","F");
  leg->AddEntry(h2ecal,"ECAL map","F");

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
  // Layer 2, Modules 1-4, Ladders -7 to -6 out of 2x14
  // Layer 3, Modules 1-4, Ladders -11 to -9 out of 2x22 
  // double dy0p = -TMath::Pi()/3.;
  //double dy0p = -TMath::Pi()/6.;
  //double dy0p = -TMath::Pi()/4.;
  double dy0p = -9*TMath::TwoPi()/72.;
  double dy2 = TMath::Pi()/14.;
  double dy3 = TMath::Pi()/22.;
  TBox HBP2(0,  (-7)*dy2+dy0p, 1.31, (-6+1)*dy2+dy0p);
  TBox HBP3(0, (-11)*dy3+dy0p, 1.31, (-9+1)*dy3+dy0p);
  TBox HBP2x(0, (+6-1)*dy2-dy0p, 1.31, (+7)*dy2-dy0p);
  TBox HBP3x(0, (+9-1)*dy3-dy0p, 1.31, (+11)*dy3-dy0p);
  // Layer 2, Modules -3 to -2, Ladder 5 out of 2x14
  // Layer 4, Module -2, Ladder 10 out of 2x32
  //double dy0m = -TMath::Pi()/6.;
  double dy0m = -9*TMath::TwoPi()/72.;
  double dx = 1.31 / 4.;
  double dy4 = TMath::Pi()/32.;
  TBox HBM2((-3)*dx,  (+5-1)*dy2+dy0m, (-2+1)*dx,  (+5)*dy2+dy0m);
  TBox HBM4((-2)*dx, (+10-1)*dy4+dy0m, (-2+1)*dx, (+10)*dy4+dy0m);

  HBP2.SetFillStyle(kNone);
  HBP2.SetLineColor(kRed+2);
  HBP2.Draw("SAME");
  HBP3.SetFillStyle(kNone);
  HBP3.SetLineColor(kOrange+2);
  HBP3.SetLineStyle(kDashed);
  HBP3.Draw("SAME");

  HBP2x.SetFillStyle(kNone);
  HBP2x.SetLineColor(kRed+2);
  HBP2x.SetLineStyle(kDashDotted);
  //HBP2x.Draw("SAME");
  HBP3x.SetFillStyle(kNone);
  HBP3x.SetLineColor(kOrange+2);
  HBP3x.SetLineStyle(kDotted);
  //HBP3x.Draw("SAME");

  HBM2.SetFillStyle(kNone);
  HBM2.SetLineColor(kRed+2);
  HBM2.Draw("SAME");
  HBM4.SetFillStyle(kNone);
  HBM4.SetLineColor(kOrange+2);
  HBM4.SetLineStyle(kDashed);
  HBM4.Draw("SAME");
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

  
  // HO kuumat:
  // https://twiki.cern.ch/twiki/pub//CMSPublic/HcalDPGResultsCMSDPS2017017/const_term_2d.png:
  // - HBM iphi 55, -6 to -5, plus -10 (hot)
  // - HBP iphi 23, +5 to 10 and 14 (hot)
  // - HBM iphi 57, -9 to -7 (cold)
  // https://cds.cern.ch/record/357153/files/CMS_HCAL_TDR.pdf (18 wedges, 4dphi)
  // "wedges are numbered 1 through 18 for each half-barrel, starting with 1 at the x axis and proceeding counterclockwise towards positive y"
  // "The CMS experiment uses a right-handed coordinate system, with the origin at the nominal collision point, the x-axis pointing to the centre of the LHC ring, the y-axis pointing up (perpendicular to the LHC plane), and the z-axis along the anticlockwise beam direction."
  //
  // From Salavat Abdoulline, 18 Nov 2020
  // (1) iphi=1 : phi=0.0437 (HO, HB, HE(S) with N_phi=72)
  //     : phi=0.0873 (HE(D) and HF with N_phi=36, except extreme ieta=40,41).
  //
  // NB: iphi=36 : phi=3.1 ->  iphi=37 : phi=-3.1 ->  iphi=72 : phi=-0.0447

  //TBox HBM2x1(-1.31,0.3927,0,0.4800); // wedge 2, tower 1 (iphi=5)
  TBox HBM2x1(-1.31,0.4363,0,0.5236); // wedge 2, tower 2 (iphi=6)
  HBM2x1.SetFillStyle(kNone);
  HBM2x1.SetLineColor(kRed+2);
  HBM2x1.Draw("SAME");

  //TBox HBP12x2(0,-2.4000,1.31,-2.225); // wedge 12, towers 1-2 (iphi=45-46)
  TBox HBP12x2(0,-2.4435,1.31,-2.2689); // wedge 12, towers 1-2 (iphi=45-46)
  HBP12x2.SetFillStyle(kNone);
  HBP12x2.SetLineColor(kRed+2);
  HBP12x2.Draw("SAME");

  //TBox HBP6x2(0,1.9199,1.31,2.0944); // wedge 6 towers 3-4 (iphi=23-24)
  //HBP6x2.SetFillStyle(kNone);
  //HBP6x2.SetLineColor(kRed-9);
  //HBP6x2.Draw("SAME");
  //TBox HBP6x1(0,1.9199,1.31,2.0071); // wedge 6 tower 3 (iphi=23)
  TBox HBP6x1(0,2.0071,1.31,2.0944); // wedge 6 tower 4 (iphi=24)
  HBP6x1.SetFillStyle(kNone);
  HBP6x1.SetLineColor(kRed-9);
  HBP6x1.Draw("SAME");

  double dphi = TMath::TwoPi()/72.;
  double dphi2 = dphi/2.;
  //double dphiHF = TMath::TwoPi()/36.;
  //double dphiHF2 = dphiHF/2.;
  double deta = TMath::TwoPi()/72.;
  double deta2 = deta/2.;
  double detaHF = 2*deta;
  double detaHF2 = detaHF/2.;
  //double phi0 = -TMath::Pi()-TMath::Pi()/2.5;
  double phi0 = 0;//dphi*23;

  double phi55 = 55*dphi - dphi2 + phi0;
  if (phi55>+TMath::Pi()) phi55 -= TMath::TwoPi();
  if (phi55<-TMath::Pi()) phi55 += TMath::TwoPi();
  double etam10 = -10*deta + deta2;
  double etam5 = -5*deta + deta2;

  double phi23 = 23*dphi - dphi2 + phi0;
  if (phi23>TMath::Pi()) phi23 -= TMath::TwoPi();
  double etap5 = 5*deta - deta2;
  double etap14 =14*deta - deta2;

  double phi57 = 57*dphi - dphi2 + phi0;
  if (phi57>TMath::Pi()) phi57 -= TMath::TwoPi();
  double etam9 = -9*deta + deta2;
  double etam7 = -7*deta + deta2;

  double phi39 = 39*dphi - dphi2;
  if (phi39>TMath::Pi()) phi39 -= TMath::TwoPi();
  double etap30 = min(30,25)*deta + max(30-25,0)*detaHF - deta2;
  double etap34 = min(34,25)*deta + max(34-25,0)*detaHF - deta2;

  TBox HBM55(etam10-deta2, phi55-dphi2, etam5+deta2, phi55+dphi2);
  HBM55.SetFillStyle(kNone);
  HBM55.SetLineColor(kRed+2);
  //HBM55.Draw("SAME");
  
  TBox HBP23(etap5-deta2, phi23-dphi2, etap14+deta2, phi23+dphi2);
  HBP23.SetFillStyle(kNone);
  HBP23.SetLineColor(kRed+2);
  //HBP23.Draw("SAME");
  
  TBox HBM57(etam9-deta2, phi57-dphi2, etam7+deta2, phi57+dphi2);
  HBM57.SetFillStyle(kNone);
  HBM57.SetLineColor(kBlue+2);
  //HBM57.Draw("SAME");
  
  //TBox QIE11(etap30-detaHF2, phi39-dphi2, etap34+detaHF2, phi39+dphi2); 
  TBox QIE11(etap30-detaHF2, phi39-dphi2, etap34+detaHF2, phi39+3*dphi2); 
  QIE11.SetFillStyle(kNone);
  QIE11.SetLineColor(kMagenta+1);
  QIE11.Draw("SAME");


  TH2D *h2hbm2 = (TH2D*)h2sum->Clone("h2hot_hbm2"); h2hbm2->Reset();
  TH2D *h2hbp12 = (TH2D*)h2sum->Clone("h2hot_hbp12"); h2hbp12->Reset();
  TH2D *h2qie11 = (TH2D*)h2sum->Clone("h2cold_qie11"); h2qie11->Reset();
  TH2D *h2mchot = (TH2D*)h2sum->Clone("h2hot_mc"); h2mchot->Reset();
  TH2D *h2ref = (TH2D*)h2sum->Clone("h2hot_ul16_plus_hbm2_hbp12_qie11");
  for (int i = 1; i != h2sum->GetNbinsX()+1; ++i) {
    for (int j = 1; j != h2sum->GetNbinsY()+1; ++j) {
      double eta = h2sum->GetXaxis()->GetBinCenter(i);
      double phi = h2sum->GetYaxis()->GetBinCenter(j);
      if (eta>HBM2x1.GetX1() && eta<HBM2x1.GetX2() &&
	  phi>HBM2x1.GetY1() && phi<HBM2x1.GetY2()) {
	h2hbm2->SetBinContent(i, j, 10);
	h2ref->SetBinContent(i, j, 10);
      }
      if (eta>HBP12x2.GetX1() && eta<HBP12x2.GetX2() &&
	  phi>HBP12x2.GetY1() && phi<HBP12x2.GetY2()) {
	h2hbp12->SetBinContent(i, j, 10);
	h2ref->SetBinContent(i, j, 10);
      }
      if (eta>QIE11.GetX1() && eta<QIE11.GetX2() &&
	  phi>QIE11.GetY1() && phi<QIE11.GetY2()) {
	h2qie11->SetBinContent(i, j, 10);
	h2ref->SetBinContent(i, j, 10);
      }
      if (eta>HBP6x1.GetX1() && eta<HBP6x1.GetX2() &&
	  phi>HBP6x1.GetY1() && phi<HBP6x1.GetY2()) {
	h2mchot->SetBinContent(i, j, 10);
	//h2ref->SetBinContent(i, j, 10);
      }
    } // for i
  } // for j



  gPad->Paint();

  c1->SaveAs("pdf/hotjets2016UL.pdf");

  TFile *fout = new TFile("rootfiles/hotjets-UL16.root","RECREATE");
  h2sum->Write();
  h2ref->Write();
  h2hbm2->Write();
  h2hbp12->Write();
  h2qie11->Write();
  h2mchot->Write();
  fout->Close();
}
