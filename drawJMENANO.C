// Purpose: draw L2Res from JMENANO at various pTs, merging triggers
#include "drawAlCaJME.C" // get fixB

#include "tdrstyle_mod22.C"

class trg {
public:
  string name;
  double ptmin;
  double ptmax;
  double etamin;
  double etamax;
  trg(string _name, double _ptmin, double _ptmax,
      double _etamin, double _etamax) {
    name = _name;
    ptmin = _ptmin;
    ptmax = _ptmax;
    etamin = _etamin;
    etamax = _etamax;
  }
};

void drawJMENANO(bool isMC = true) {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TFile *fd(0);
  if (isMC)
    fd = new TFile("../dijet/rootfiles/jmenano_mc_out_v5.root","READ");
  else
    fd = new TFile("../dijet/rootfiles/jmenano_data_out_v6.root","READ");
  //fd = new TFile("../dijet/rootfiles/jmenano_data_out_v5_4p86fb.root","READ");
  assert(fd && !fs->IsZombie());

  const char *cs = (isMC ? "MC" : "DT");

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);

  //TH1D *h = tdrHist("h","MPF",0.5,1.2,"#eta_{jet}",-5.2,5.2);
  TH1D *h = tdrHist("h","p_{T,ave}",15,3400,"#eta_{jet}",-5.2,5.2);
  //lumi_136TeV = "JMENANO 0.450 out of 4.86 fb^{-1}";
  //lumi_136TeV = "RunC, 0.450 of 4.86 fb^{-1}";
  lumi_136TeV = (isMC ? "Flat2018QCD" : "RunC, 2.90 of 4.86 fb^{-1}");
  extraText = "WIP";
  TCanvas *c1 = tdrCanvas("c1",h,8,0,kSquare);

  gPad->SetRightMargin(0.20);
  gPad->SetLeftMargin(0.20);
  h->GetYaxis()->SetTitleOffset(1.6);
  gPad->SetLogy();
  h->GetYaxis()->SetMoreLogLabels();
  h->GetYaxis()->SetNoExponent();

   // Copy L2Res histograms for multiple pT bins
   const double vpt[] = {40, 50, 60, 70, 85, 100, 125, 155, 180, 210, 250, 300,
			 350, 400, 500, 600, 800, 1000, 1200, 1500,
			 1800, 2100, 2400, 2700, 3000};
   const int npt = sizeof(vpt)/sizeof(vpt[0])-1;

  // Listing of available triggers
   vector<trg> vtrg;
   vtrg.push_back(trg("HLT_DiPFJetAve40",40,70,0,5.2));//40,60
   vtrg.push_back(trg("HLT_DiPFJetAve60",70,100,0,5.2)); //60,85
   vtrg.push_back(trg("HLT_DiPFJetAve80",100,155,0,5.2)); //85,155
   vtrg.push_back(trg("HLT_DiPFJetAve140",155,250,0,5.2)); //180,210
   vtrg.push_back(trg("HLT_DiPFJetAve200",250,300,0,5.2)); // 250,300
   vtrg.push_back(trg("HLT_DiPFJetAve260",300,350,0,5.2));
   vtrg.push_back(trg("HLT_DiPFJetAve320",350,400,0,5.2));
   vtrg.push_back(trg("HLT_DiPFJetAve400",400,500,0,5.2));
   vtrg.push_back(trg("HLT_DiPFJetAve500",500,3000,0,5.2));
   /*
   vtrg.push_back("HLT_PFJet40");
   vtrg.push_back("HLT_PFJet60");
   vtrg.push_back("HLT_PFJet80");
   //vtrg.push_back("HLT_PFJet110");
   vtrg.push_back("HLT_PFJet140");
   vtrg.push_back("HLT_PFJet200");
   vtrg.push_back("HLT_PFJet260");
   vtrg.push_back("HLT_PFJet320");
   vtrg.push_back("HLT_PFJet450");
   vtrg.push_back("HLT_PFJet500");
   vtrg.push_back("HLT_PFJet550");

   vtrg.push_back("HLT_DiPFJetAve60_HFJEC");
   vtrg.push_back("HLT_DiPFJetAve80_HFJEC");
   vtrg.push_back("HLT_DiPFJetAve100_HFJEC");
   vtrg.push_back("HLT_DiPFJetAve160_HFJEC");
   vtrg.push_back("HLT_DiPFJetAve220_HFJEC");
   vtrg.push_back("HLT_DiPFJetAve300_HFJEC");

   vtrg.push_back("HLT_PFJetFwd15");
   vtrg.push_back("HLT_PFJetFwd25");
   vtrg.push_back("HLT_PFJetFwd40");
   vtrg.push_back("HLT_PFJetFwd60");
   vtrg.push_back("HLT_PFJetFwd80");
   vtrg.push_back("HLT_PFJetFwd140");
   vtrg.push_back("HLT_PFJetFwd200");
   vtrg.push_back("HLT_PFJetFwd260");
   vtrg.push_back("HLT_PFJetFwd320");
   vtrg.push_back("HLT_PFJetFwd400");
   vtrg.push_back("HLT_PFJetFwd450");
   vtrg.push_back("HLT_PFJetFwd500");
   */
   int ntrg = vtrg.size();
  
   // Regular L2Relative and L2Res eta binning
   double vx[] =
     {-5.191,
      -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, -3.489,
      -3.314, -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172, -2.043,
      -1.93, -1.83, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218,
      -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, -0.435, 
      -0.348, -0.261, -0.174, -0.087, 0, 0.087, 0.174, 0.261, 0.348, 0.435,
      0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305,
      1.392, 1.479, 1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5,
      2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191,
      4.363, 4.538, 4.716, 4.889, 5.191};
   const int nx = sizeof(vx)/sizeof(vx[0])-1;

   TH2D *h2m = new TH2D("h2m",";#eta_{jet};p_{T,ave};MPF",nx,vx,npt,vpt);
   //TH2D *h2ms = new TH2D("h2ms",";#eta_{jet};p_{T,ave};MPF",nx,vx,npt,vpt);
   TH2D *h2d = new TH2D("h2d",";#eta_{jet};p_{T,ave};DB",nx,vx,npt,vpt);

   for (int itrg = 0; itrg != ntrg; ++itrg) {

     trg &t = vtrg[itrg];
     fd->cd(t.name.c_str());
     TDirectory *d = gDirectory;

     for (int ipt = 0; ipt != npt; ++ipt) {
       
       double ptmin = vpt[ipt];
       double ptmax = vpt[ipt+1];

       if (ptmin >= t.ptmin && ptmax <= t.ptmax) {
   
	 d->cd(Form("Pt_%d_%d",int(ptmin+0.5),int(ptmax+0.5)));
	 TDirectory *dd = gDirectory;
	 TProfile *pm = (TProfile*)dd->Get("pm0ab"); assert(pm);
	 TProfile *pd = (TProfile*)dd->Get("pm2ab"); assert(pd);
	 
	 for (int ieta = 0; ieta != nx; ++ieta) {

	   double etamin = vx[ieta];
	   double etamax = vx[ieta+1];

	   if ((etamin>=0 && (etamin >= t.etamin && etamax <= t.etamax)) ||
	       (etamin<0 && (etamin >= -t.etamax && etamax <= -t.etamin))) {

	     h2m->SetBinContent(ieta+1,ipt+1, pm->GetBinContent(ieta+1));
	     h2m->SetBinError(ieta+1,ipt+1, pm->GetBinError(ieta+1));
	   } // eta
	 } // ieta
       } // pt
     } // ipt
   } // itrg 

   h2m->Draw("SAME COLZ");
   h2m->GetZaxis()->SetRangeUser(0.5+1e-4,1.5-1e-4);
   h2m->GetZaxis()->SetTitleOffset(1.2);
   gPad->RedrawAxis();
   gPad->Update();

   c1->SaveAs(Form("pdf/alcajme/drawJMENANO_2D_%s.pdf",cs));
   h2m->GetZaxis()->SetRangeUser(0.8+1e-4,1.2-1e-4);
   c1->SaveAs(Form("pdf/alcajme/drawJMENANO_2D_zoom_%s.pdf",cs));
   h2m->GetZaxis()->SetRangeUser(0.9+1e-4,1.1-1e-4);
   c1->SaveAs(Form("pdf/alcajme/drawJMENANO_2D_zoom2_%s.pdf",cs));

   // Run smoothing vs pT, keep track of chi2/NDF
   TH2D *h2ms = (TH2D*)h2m->Clone("h2ms");
   h2ms->Reset();
   TH1D *hchi20 = h2ms->ProjectionX("hchi20");
   TH1D *hchi21 = h2ms->ProjectionX("hchi21");
   TH1D *hjes40 = h2ms->ProjectionX("hjes40");
   TH1D *hjes15b = h2ms->ProjectionX("hjes15b");
   TH1D *hjesmaxb = h2ms->ProjectionX("hjesmaxb");
   TH1D *hjes40orig = h2m->ProjectionX("hjes40orig",1,1);
   double i100 = h2m->GetYaxis()->FindBin(100.);
   TH1D *hjes100orig = h2m->ProjectionX("hjes100orig",i100,i100);
   TH1D *hdpt = h2ms->ProjectionX("hdpt");
   TH1D *hdpt0 = h2ms->ProjectionX("hdpt0");
   for (int i = 1; i != h2m->GetNbinsX()+1; ++i) {

     double eta = h2m->GetXaxis()->GetBinCenter(i);
     double emax1 = 6800.*0.5;
     double emax2 = 4500.;
     double ptmax = (fabs(eta)>2.5 ? emax2 : emax1) / cosh(eta);
     TH1D *h = h2m->ProjectionY("htmp",i,i);
     
     // Quick cleaning
     for (int i = 1; i != h->GetNbinsX()+1; ++i) {
       if (h->GetBinError(i)<0.001) {
	 h->SetBinContent(i, 0.);
	 h->SetBinError(i, 0.);
       }
     }

     TF1 *f0 = new TF1(Form("f0_%d",i),"[0]",40,ptmax);
     h->Fit(f0,"QRN");
     TF1 *f1 = new TF1(Form("f1_%d",i),"[0]+[1]*log(x)",40,ptmax);
     f1->SetParameters(f0->GetParameter(0), 0.);
     h->Fit(f1,"QRN");

     hchi20->SetBinContent(i, f0->GetNDF()>0 ? f0->GetChisquare()/f0->GetNDF() : 0);
     hchi20->SetBinError(i, f0->GetNDF()>0 ? f0->GetChisquare()/f0->GetNDF()*1./sqrt(f0->GetNDF()) : 0);

     hchi21->SetBinContent(i, f1->GetNDF()>0 ? f1->GetChisquare()/f1->GetNDF() : 0);
     hchi21->SetBinError(i, f1->GetNDF()>0 ? f1->GetChisquare()/f1->GetNDF()*1./sqrt(f1->GetNDF()) : 0);

     double nsig = 2.5;//1.5;
     TF1 *fjes = (fabs(f1->GetParameter(1))>nsig*f1->GetParError(1) ? f1 : f0);
     if (isMC && (fabs(eta)>2.5 || fabs(eta)<1.0)) fjes = f0;
     if (!isMC && (fabs(eta)<1.0)) fjes = f0;
  
     hjes40->SetBinContent(i, fjes->GetNDF()>=0 ? fjes->Eval(40.) : 1);

     hjes15b->SetBinContent(i, fjes->GetNDF()>=0 ?
			    0.5*(fjes->Eval(40.)+fjes->Eval(15.)) : 1);
     hjes15b->SetBinError(i, fjes->GetNDF()>=0 ?
			  0.5*fabs(fjes->Eval(40.)-fjes->Eval(15.)) : 1);

     hjesmaxb->SetBinContent(i, fjes->GetNDF()>=0 ?
			     0.5*(fjes->Eval(40.)+fjes->Eval(ptmax)) : 1);
     hjesmaxb->SetBinError(i, fjes->GetNDF()>=0 ?
			   0.5*fabs(fjes->Eval(40.)-fjes->Eval(ptmax)) : 1);

     hdpt->SetBinContent(i, f1->GetNDF()>=0 ? 
			 100.*f1->GetParameter(1) : 0);
     hdpt->SetBinError(i, f1->GetNDF()>=0 ? 
		       100.*f1->GetParError(1) : 0);

     hdpt0->SetBinContent(i, fjes==f1 ? 100.*f1->GetParameter(1) : 0);
     hdpt0->SetBinError(i, fjes==f1 ? 100.*f1->GetParError(1) : 0);

     for (int j = 1; j != h2m->GetNbinsY()+1; ++j) {
       if (//h2m->GetBinContent(i,j)!=0 &&
	   h2m->GetYaxis()->GetBinLowEdge(j)<ptmax) {
	 double pt = h2m->GetYaxis()->GetBinCenter(j);
	 h2ms->SetBinContent(i, j, f1->Eval(pt));
       }
       else {
	 h2ms->SetBinContent(i, j, 0.);
       }
     }

     delete h;
   } // for i


   TH1D *hs = tdrHist("hs","p_{T,ave}",15,3400,"#eta_{jet}",-5.2,5.2);
   TCanvas *c1s = tdrCanvas("c1s",hs,8,0,kSquare);

   gPad->SetRightMargin(0.20);
   gPad->SetLeftMargin(0.20);
   hs->GetYaxis()->SetTitleOffset(1.6);
   gPad->SetLogy();
   hs->GetYaxis()->SetMoreLogLabels();
   hs->GetYaxis()->SetNoExponent();

   h2ms->Draw("SAME COLZ");
   h2ms->GetZaxis()->SetRangeUser(0.5+1e-4,1.5-1e-4);
   gPad->RedrawAxis();
   gPad->Update();

   c1s->SaveAs(Form("pdf/alcajme/drawJMENANO_2DS_%s.pdf",cs));
   h2ms->GetZaxis()->SetRangeUser(0.8+1e-4,1.2-1e-4);
   c1s->SaveAs(Form("pdf/alcajme/drawJMENANO_2DS_zoom_%s.pdf",cs));
   h2ms->GetZaxis()->SetRangeUser(0.9+1e-4,1.1-1e-4);
   c1s->SaveAs(Form("pdf/alcajme/drawJMENANO_2DS_zoom2_%s.pdf",cs));


   TH1D *h2 = tdrHist("h2","#chi^{2} / NDF",0,10,"#eta_{jet}",-5.2,5.2);
   TCanvas *c2 = tdrCanvas("c2",h2,8,0,kSquare);
   tdrDraw(hchi20,"HISTE",kNone,kOrange+2,kSolid,-1,1001,kOrange+1);
   tdrDraw(hchi21,"HISTE",kNone,kYellow+2,kSolid,-1,1001,kYellow+1);
   hchi21->SetFillColorAlpha(kYellow+1,0.7);
   gPad->RedrawAxis();
   l->DrawLine(-5.2,1,+5.2,1);
   c2->SaveAs(Form("pdf/alcajme/drawJMENANO_chi2_%s.pdf",cs));

   TH1D *h3 = tdrHist("h3","JES at p_{T}=40 GeV (p_{T,min} to E_{max})",
		      0.38+1e-4,1.28-1e-4,"#eta_{jet}",-5.2,5.2);
   TCanvas *c3 = tdrCanvas("c3",h3,8,0,kSquare);
   h3->GetXaxis()->SetTitleOffset(0.8);
   h3->GetYaxis()->SetTitleOffset(1.18);
   tdrDraw(hjes40,"HIST",kNone,kYellow+2,kSolid,-1,1001,kYellow+1);
   tdrDraw(hjes15b,"E2",kNone,kCyan+2,kSolid,-1,1001,kCyan+1);
   hjes15b->SetFillColorAlpha(kCyan+1,0.7);
   tdrDraw(hjesmaxb,"E2",kNone,kOrange+2,kSolid,-1,1001,kOrange+1);
   hjesmaxb->SetFillColorAlpha(kOrange+1,0.7);
   if (!isMC) tdrDraw(hjes40orig,"P",kFullCircle,kBlack,kSolid,-1,kNone,0);
   if (isMC)  tdrDraw(hjes100orig,"P",kFullCircle,kBlack,kSolid,-1,kNone,0);
   hjes40orig->SetMarkerSize(0.5);
   hjes100orig->SetMarkerSize(0.5);
   gPad->RedrawAxis();

   TLegend *leg = tdrLeg(0.35,0.90-4*0.035,0.55,0.90);
   leg->SetTextSize(0.035);
   leg->AddEntry(hjes40orig,isMC ? "JES @ 100 GeV" : "JES @ 40 GeV","PLE");
   leg->AddEntry(hjes40,"Fit @ 40 GeV","F");
   leg->AddEntry(hjes15b,"Fit @ 15 GeV","F");
   leg->AddEntry(hjesmaxb,"Fit @ E_{max}","F");

   l->DrawLine(-5.2,1,+5.2,1);

   c3->SaveAs(Form("pdf/alcajme/drawJMENANO_jes40_%s.pdf",cs));

   TH1D *h4 = tdrHist("h4","dJES/dlog(p_{T}) (%)",-30,30,"#eta_{jet}",-5.2,5.2);
   TCanvas *c4 = tdrCanvas("c4",h4,8,0,kSquare);
   tdrDraw(hdpt,"HISTE",kNone,kYellow+2,kSolid,-1,1001,kYellow+1);
   tdrDraw(hdpt0,"HISTE",kNone,kOrange+2,kSolid,-1,1001,kOrange+1);
   gPad->RedrawAxis();

   l->SetLineStyle(kSolid);
   l->DrawLine(-5.2,0,+5.2,0);

   c4->SaveAs(Form("pdf/alcajme/drawJMENANO_dpt_%s.pdf",cs));
} // drawJMENANO
