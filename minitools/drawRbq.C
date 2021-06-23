// Purpose: explore the R_bq variables introduced by ATLAS
// Run with 'root -l -b -q minitools/drawRbq.C+g'
#include "TFile.h"
#include "TProfile.h"
#include "TH2D.h"
#include "TF1.h"

#include "../tdrstyle_mod15.C"

void invertY(TH1D *h) {
  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    double y = h->GetBinContent(i);
    double ey = h->GetBinError(i);
    h->SetBinContent(i, (y!=0 ? 1./y : 0.));
    h->SetBinError(i, (y!=0 ? 1./y * ey/y : 0.));
  } // for i
} // invertY

void solveY(TH1D *h, double c) {
  // Y = 1+((b1+b2)-c*(p1+p2))/(b1+b2+c*(p1+p2)) vs (b1+b2+c*(p1+p2))
  // => (b1+b2)/(p1+p2) = c*Y / (2-Y)
  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    double y = h->GetBinContent(i);
    double ey = h->GetBinError(i);
    double F = ((2-y)!=0 ? c*y/(2-y) : 0.);
    h->SetBinContent(i, F); 
    //h->SetBinError(i, (2-y!=0 ? c*y/(2-y) * ey/y : 0.));
    h->SetBinError(i, (F+1./c*F*F) * ey/y);
  } // for i
} // solveY

TH1D *diffY(TH1D *h1, TH1D *h2, const char* name) {
  
  TH1D *h = (TH1D*)h1->Clone(name);
  h->Add(h1,h2,1,-1);
  h->Divide(h2);
  h->Scale(100.);
  
  return h;
}

void drawRbq(string mode="18V5") {

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  const char *cm = mode.c_str();

  //TFile *fm = new TFile("rootfiles/hadWMC18V4.root","READ");
  TFile *fm = new TFile(Form("rootfiles/hadWMC%s_MPDGcorrNoW.root",cm),"READ");
  assert(fm && !fm->IsZombie());
  //TFile *fd = new TFile("rootfiles/hadWUL18V4.root","READ");
  TFile *fd = new TFile(Form("rootfiles/hadWUL%s_MPDGcorrNoW.root",cm),"READ");
  assert(fd && !fd->IsZombie());

  // MC truth bJES
  TProfile *pg1 = (TProfile*)fm->Get("pb1"); assert(pg1);
  TProfile *pb1 = (TProfile*)fm->Get("pbq1"); assert(pb1);
  TProfile *pq1 = (TProfile*)fm->Get("pbb1"); assert(pq1);
  TProfile *pa1 = (TProfile*)fm->Get("pba1"); assert(pa1);

  // Differential R_bq estimators
  TH1D *hrbqm = (TH1D*)fm->Get("hrbq"); assert(hrbqm);
  TH2D *h2bqm = (TH2D*)fm->Get("h2bq"); assert(h2bqm);
  TProfile *prqm = (TProfile*)fm->Get("prbqq"); assert(prqm);
  TProfile *prbm = (TProfile*)fm->Get("prbqb"); assert(prbm);
  TProfile *pram = (TProfile*)fm->Get("prbqa"); assert(pram);
  TProfile *prcm = (TProfile*)fm->Get("prbqc"); assert(prcm);
  //
  TH1D *hrbqd = (TH1D*)fd->Get("hrbq"); assert(hrbqd);
  TH2D *h2bqd = (TH2D*)fd->Get("h2bq"); assert(h2bqd);
  TProfile *prqd = (TProfile*)fd->Get("prbqq"); assert(prqd);
  TProfile *prbd = (TProfile*)fd->Get("prbqb"); assert(prbd);
  TProfile *prad = (TProfile*)fd->Get("prbqa"); assert(prad);
  TProfile *prcd = (TProfile*)fd->Get("prbqc"); assert(prcd);

  // NB: need to add mapping to <b1+b2> on x-axis
  // (b1+b2)/(p1+p2) vs (p1+p2) => y-axis ok
  TH1D *hrqm = prqm->ProjectionX("hrqm");
  TH1D *hrqd = prqd->ProjectionX("hrqd");
  // (p1+p2)/(b1+b2) vs (b1+b2) => swap y-axis
  TH1D *hrbm = prbm->ProjectionX("hrbm");
  TH1D *hrbd = prbd->ProjectionX("hrbd");
  invertY(hrbm);
  invertY(hrbd);
  // 1+((b1+b2)-(p1+p2))/(b1+b2+p1+p2) vs (b1+b2+p1+p2) => solve for Y
  // => 2*(b1+b2) / (b1+b2+p1+p2) = Y
  // => -2*(p1+p2) / (b1+b2+p1+p2) = Y-2
  // => (b1+b2)/(p1+p2) = Y / (2-Y)
  TH1D *hram = pram->ProjectionX("hram");
  TH1D *hrad = prad->ProjectionX("hrad");
  solveY(hram,1.0);
  solveY(hrad,1.0);
  // 1+((b1+b2)-c*(p1+p2))/(b1+b2+c*(p1+p2)) vs (b1+b2+c*(p1+p2)) => solve for Y
  // => 2*(b1+b2) / (b1+b2+c*(p1+p2)) = Y
  // => -2*c*(p1+p2) / (b1+b2+c*(p1+p2)) = Y-2
  // => (b1+b2)/(p1+p2) = c*Y / (2-Y)
  TH1D *hrcm = prcm->ProjectionX("hrcm");
  TH1D *hrcd = prcd->ProjectionX("hrcd");
  solveY(hrcm,1.2);
  solveY(hrcd,1.2);

  TH1D *hrq = diffY(hrqd,hrqm,"hrq");
  TH1D *hrb = diffY(hrbd,hrbm,"hrb");
  TH1D *hra = diffY(hrad,hram,"hra");
  TH1D *hrc = diffY(hrcd,hrcm,"hrc");
  
  //////////////////////////////
  // MC truth b jet response  //
  //////////////////////////////

  // NB1: Need to fix <1/pT,gen> to 1/<pT,gen> for pb1,pq1,pa1
  // NB2: Also need to map p_{T,ref} to <pT,gen>
  TH1D *h1 = tdrHist("h1","#LTp_{T,b,rec}/p_{T,b,gen}#GT",
		     0.90,1.15,"p_{T,ref} (GeV)",30,230);

  TCanvas *c1 = tdrCanvas("c1",h1,4,11,kSquare);
  gPad->SetLogx();

  tdrDraw(pg1,"Pz",kFullCircle,kBlack);
  tdrDraw(pb1,"Pz",kOpenTriangleUp,kRed);
  tdrDraw(pq1,"Pz",kOpenTriangleDown,kBlue);
  tdrDraw(pa1,"Pz",kOpenDiamond,kGreen+2);

  TLegend *leg1 = tdrLeg(0.65,0.66,0.85,0.90);
  leg1->AddEntry(pg1,"vs b,gen","PLE");
  leg1->AddEntry(pb1,"vs b,ave","PLE");
  leg1->AddEntry(pq1,"vs q,ave","PLE");
  leg1->AddEntry(pa1,"vs b+q,ave","PLE");

  TLatex *tex = new TLatex();
  tex->SetTextSize(0.045); tex->SetNDC();
  tex->DrawLatex(0.25,0.73,"Anti-k_{T} R=0.4");
  tex->DrawLatex(0.25,0.67,"|#eta|<1.3, p_{T}>30 GeV");



  /////////////////////////////////
  // Differential R_bq, data/MC  //
  ////////////////////////////////

  TH1D *h2up = tdrHist("h2up","#LTp_{T,b1}+p_{T,b2}#GT /"
		       " #LTp_{T,q1}+p_{T,q2}#GT",
		       0.00,3.50,"p_{T,ref} (GeV)",30,230);
  TH1D *h2dw = tdrHist("h2dw","Data/MC-1 (%)",
		       /*-7,+7,*/-1.5,1.5,"p_{T,ref} (GeV)",30,230);

  //TCanvas *c2 = tdrCanvas("c2",h2,4,11,kSquare);
  TCanvas *c2 = tdrDiCanvas("c2",h2up,h2dw,4,11);
  c2->cd(1);
  gPad->SetLogx();

  tdrDraw(hrqm,"Pz",kOpenTriangleUp,kBlue);
  tdrDraw(hrbm,"Pz",kOpenTriangleDown,kRed);
  tdrDraw(hram,"Pz",kOpenDiamond,kGreen+2);
  tdrDraw(hrcm,"Pz",kOpenDiamond,kCyan+2);
  //
  tdrDraw(hrqd,"Pz",kFullTriangleUp,kBlue);
  tdrDraw(hrbd,"Pz",kFullTriangleDown,kRed);
  tdrDraw(hrad,"Pz",kFullDiamond,kGreen+2);
  tdrDraw(hrcd,"Pz",kFullDiamond,kCyan+2);

  TLegend *leg2 = tdrLeg(0.65,0.66,0.85,0.90);
  //TLegend *leg2 = tdrLeg(0.65,0.60,0.85,0.90);
  leg2->AddEntry(hrqm,"vs q,ave","PLE");
  leg2->AddEntry(hrbm,"vs b,ave","PLE");
  leg2->AddEntry(hram,"vs b+q,ave","PLE");
  leg2->AddEntry(hrcm,"vs b+q,wgt","PLE");
  //leg2->AddEntry(,"R_{bq}","PLE");


  tex->DrawLatex(0.25,0.73,"Anti-k_{T} R=0.4");
  tex->DrawLatex(0.25,0.67,"|#eta|<1.3, p_{T}>30 GeV");

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  double rbq = hrbqm->GetMean();
  l->DrawLine(30,rbq,230,rbq);
  l->SetLineStyle(kDotted);
  double ptq = h2bqm->GetMean(1);
  double ptb = h2bqm->GetMean(2);
  double pta = 0.5*(ptq+ptb);
  double rptbq = ptb/ptq;
  l->DrawLine(30,rptbq,230,rptbq);
  l->SetLineColor(kBlue);
  l->DrawLine(ptq,0,ptq,rptbq);
  l->SetLineColor(kRed);
  l->DrawLine(ptb,0,ptb,rptbq);
  l->SetLineColor(kGreen+2);
  l->DrawLine(pta,0,pta,rptbq);

  TF1 *fq = new TF1("fq","[0]+log(x)*([1]+log(x)*([2]+[3]*log(x)))",30,230);
  fq->SetParameters(1.2,0,0,0);
  fq->SetLineColor(kBlue);
  hrqm->Fit(fq,"QRN");
  fq->Draw("SAME");
  TF1 *fb = new TF1("fb","[0]+log(x)*([1]+log(x)*([2]+[3]*log(x)))",30,230);
  fb->SetParameters(1.2,0,0,0);
  fb->SetLineColor(kRed);
  hrbm->Fit(fb,"QRN");
  fb->Draw("SAME");
  TF1 *fa = new TF1("fa","[0]+log(x)*([1]+log(x)*([2]+[3]*log(x)))",30,230);
  fa->SetParameters(1.2,0,0,0);
  fa->SetLineColor(kGreen+2);
  hram->Fit(fa,"QRN");
  fa->Draw("SAME");
  TF1 *fc = new TF1("fc","[0]+log(x)*([1]+log(x)*([2]+[3]*log(x)))",30,230);
  fc->SetParameters(1.2,0,0,0);
  fc->SetLineColor(kCyan+2);
  hrcm->Fit(fc,"QRN");
  fc->Draw("SAME");

  c2->cd(2);
  gPad->SetLogx();

  double rbqm = hrbqm->GetMean();
  double rbqd = hrbqd->GetMean();
  double drbq = 100.*(rbqd/rbqm-1);
  double erbqm = hrbqm->GetMeanError();
  double erbqd = hrbqd->GetMeanError();
  double edrbq = 100.*(1+0.01*drbq)*sqrt(pow(erbqm/rbqm,2)+pow(erbqd/rbqd,2));

  double ptqm = h2bqm->GetMean(1);
  double ptbm = h2bqm->GetMean(2);
  double rptbqm = ptbm/ptqm;
  double ptqd = h2bqd->GetMean(1);
  double ptbd = h2bqd->GetMean(2);
  double rptbqd = ptbd/ptqd;
  double drptbq = 100.*(rptbqd/rptbqm-1);
  //
  double eptqm = h2bqm->GetMeanError(1);
  double eptbm = h2bqm->GetMeanError(2);
  double erptbqm = rptbqm*sqrt(pow(eptqm/ptqm,2)+pow(eptbm/ptbm,2));
  double eptqd = h2bqd->GetMeanError(1);
  double eptbd = h2bqd->GetMeanError(2);
  double erptbqd = rptbqm*sqrt(pow(eptqd/ptqd,2)+pow(eptbd/ptbd,2));
  double edrptbq = 100.*(1+0.01*drptbq)*sqrt(pow(erptbqm/rptbqm,2)+pow(erptbqd/rptbqd,2));

  l->SetLineColor(kBlack);
  l->SetLineStyle(kDashed);
  l->DrawLine(30,drbq,230,drbq);
  l->SetLineStyle(kDotted);
  l->DrawLine(30,drptbq,230,drptbq);

  tdrDraw(hrq,"Pz",kFullTriangleUp,kBlue);
  tdrDraw(hrb,"Pz",kFullTriangleDown,kRed);
  tdrDraw(hra,"Pz",kFullDiamond,kGreen+2);
  tdrDraw(hrc,"Pz",kFullDiamond,kCyan+2);

  TF1 *fda = new TF1("fda","[0]",30,230);
  fda->SetLineColor(kGreen+2);
  hra->Fit(fda,"QRN");
  fda->Draw("SAME");
  //TF1 *fdc = new TF1("fdc","[0]+log(x)*([1]+log(x)*([2]+[3]*log(x)))",30,230);
  //fdc->SetParameters(0.5,0,0,0);
  TF1 *fdc = new TF1("fdc","[0]",30,230);
  fdc->SetLineColor(kCyan+2);
  hrc->Fit(fdc,"QRN");
  fdc->Draw("SAME");

  tex->SetTextSize(0.045);
  tex->SetTextColor(kGreen+2);
  tex->DrawLatex(0.30,0.42,Form("%1.2f#pm%1.2f%%",fda->GetParameter(0),fda->GetParError(0)));
  tex->SetTextColor(kCyan+2);
  tex->DrawLatex(0.30,0.35,Form("%1.2f#pm%1.2f%%",fdc->GetParameter(0),fdc->GetParError(0)));
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.50,0.42,Form("%1.2f#pm%1.2f%%",drptbq,edrptbq));
  tex->DrawLatex(0.50,0.35,Form("%1.2f#pm%1.2f%%",drbq,edrbq));

  c1->SaveAs("pdf/drawRbq_bJES.pdf");
  c2->SaveAs("pdf/drawRbq.pdf");
} // drawRbq
