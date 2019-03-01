// Purpose: Compare L3res parameterizations from 2018 to 2017 and 2016
//          (For first 2018 iteration, guesses based on 2017 and 2016)
//          Code uses custom TF1, take parameters from JECdatabase .txt
// Run with 'root -l minitools/compareL3Res.C' (i.e. not within /minitools)
#include "TF1.h"

#include "../tdrstyle_mod15.C"

void compareL3Res() {

  setTDRStyle();

  // 2016: JEC 07Aug17 V18 (DATA) and V15 (MC). P8M1 (+HS1)
  // => Summer16_07Aug2017BCD_V18_DATA
  string s16 = "([0]+[1]*(100./3.)*(TMath::Max(0.,1.03091-0.051154*TMath::Power(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[2]*((1-(1./x)*(3.906-1.652*TMath::Log(x)+0.2257*TMath::Power(TMath::Log(x),2)))-(1-(1./208.)*(3.906-1.652*TMath::Log(208.)+0.2257*TMath::Power(TMath::Log(208.),2)))))";
  TF1 *f16bcd = new TF1("f16bcd",s16.c_str(),10,4000);
  TF1 *f16ef = new TF1("f16ef",s16.c_str(),10,4000);
  TF1 *f16gh = new TF1("f16gh",s16.c_str(),10,4000);

  // File locations
  // https://github.com/cms-jet/JECDatabase/blob/master/textFiles/Summer16_07Aug2017BCD_V18_DATA/Summer16_07Aug2017BCD_V18_DATA_L2L3Residual_AK4PFchs.txt
  // https://github.com/cms-jet/JECDatabase/blob/master/textFiles/Summer16_07Aug2017EF_V18_DATA/Summer16_07Aug2017EF_V18_DATA_L2L3Residual_AK4PFchs.txt
  // https://github.com/cms-jet/JECDatabase/blob/master/textFiles/Summer16_07Aug2017GH_V18_DATA/Summer16_07Aug2017GH_V18_DATA_L2L3Residual_AK4PFchs.txt
  f16bcd->SetParameters(0.9851, 0.1604, -2.268); // BCD
  f16ef->SetParameters(0.9796, 0.1953, -2.513); // EF
  f16gh->SetParameters(0.9909, 0.0840, -1.380); // GH

  // 2017: JEC 17Nov17 V31 (DATA) and V24 (MC). P8CP5 (+P1M1, +HS1)
  // => Fall17_17Nov2017B_V32_DATA
  string s17 = "([0]+[1]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[2]*0.021*(-1.+1./(1.+exp(-(TMath::Log(x)-5.030)/0.395))))";
  TF1 *f17b = new TF1("f17b",s17.c_str(),10,4000);
  TF1 *f17c = new TF1("f17c",s17.c_str(),10,4000);
  TF1 *f17de = new TF1("f17de",s17.c_str(),10,4000);
  TF1 *f17f = new TF1("f17de",s17.c_str(),10,4000);

  // File locations:
  // https://github.com/cms-jet/JECDatabase/blob/master/textFiles/Fall17_17Nov2017B_V32_DATA/Fall17_17Nov2017B_V32_DATA_L2L3Residual_AK4PFchs.txt (17B)
  // https://github.com/cms-jet/JECDatabase/blob/master/textFiles/Fall17_17Nov2017C_V32_DATA/Fall17_17Nov2017C_V32_DATA_L2L3Residual_AK4PFchs.txt
  // https://github.com/cms-jet/JECDatabase/blob/master/textFiles/Fall17_17Nov2017DE_V32_DATA/Fall17_17Nov2017DE_V32_DATA_L2L3Residual_AK4PFchs.txt
  // https://github.com/cms-jet/JECDatabase/blob/master/textFiles/Fall17_17Nov2017F_V32_DATA/Fall17_17Nov2017F_V32_DATA_L2L3Residual_AK4PFchs.txt
  f17b->SetParameters (0.980, 0.068, 0.0); //B
  f17c->SetParameters (0.984, 0.037, 0.0); //C
  f17de->SetParameters(0.978, 0.026, 0.0); //DE
  f17f->SetParameters (0.968, 0.031, 0.0); //F
  
  TF1 *f18a = new TF1("f18a",s17.c_str(),10,4000);
  TF1 *f18b = new TF1("f18b",s17.c_str(),10,4000);
  TF1 *f18c = new TF1("f18c",s17.c_str(),10,4000);
  TF1 *f18d = new TF1("f18d",s17.c_str(),10,4000);

  // Hannu's composition plots to guess 2018 vs 2017 L3Res:
  // http://hsiikone.web.cern.ch/hsiikone/
  // 2018AB similar, both between 2017DE and 2017F for CH+NH, but NE 1% lower
  // 2018C is similar to 2017DE, but NH bit lower at pT<200 GeV (-1% at 30 GeV)
  // 2018D is again similar to AB, but bit lower CH and higher NE
  // Guess: 18AB +1% for NH wrt 17DE
  //        18C at similar scale to 17DE
  //        18D -1% compared to 17 DE, or avefage of 17DE and 17F
  // In 2016-2017 had about -1% per 15/fb
  f18a->SetParameters(0.978+0.013, 0.026, 0.0); //17DE->18A (15/fb)
  f18b->SetParameters(0.978+0.008, 0.026, 0.0); //17DE->18B (7/fb)
  f18c->SetParameters(0.978+0.003, 0.026, 0.0); //17DE->18C (3->7/fb?)
  f18d->SetParameters(0.978-0.012, 0.026, 0.0); //17DE->18D (30/fb)
  
  lumi_13TeV = "Run II (2016-2017)";
  TH1D *h = new TH1D("h",";p_{T} (GeV);Jet response (L3Res)",399,10,4000);
  h->SetMinimum(0.96);//0.94);//0.95);
  h->SetMaximum(1.02);//1.01);
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  h->GetXaxis()->SetRangeUser(30,4000.);

  TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);
  gPad->SetLogx();

  // 2016
  f16bcd->SetLineWidth(2); f16bcd->SetLineStyle(kDotted);
  f16bcd->SetLineColor(kGray);//kBlue);
  f16bcd->DrawClone("SAME");
  f16ef->SetLineWidth(2); f16ef->SetLineStyle(kDotted);
  f16ef->SetLineColor(kGray);//kRed);
  f16ef->Draw("SAME");
  f16gh->SetLineWidth(5); f16gh->SetLineStyle(kDotted);
  f16gh->SetLineColor(kBlack);//kGreen+2);
  f16gh->Draw("SAME");

  // 2017
  f17b->SetLineWidth(2); f17b->SetLineStyle(kDashed);
  f17b->SetLineColor(kGray);//kBlack);
  f17b->Draw("SAME");
  f17c->SetLineWidth(5); f17c->SetLineStyle(kDashed);
  f17c->SetLineColor(kBlue);
  f17c->Draw("SAME");
  f17de->SetLineWidth(5); f17de->SetLineStyle(kDashed);
  f17de->SetLineColor(kGreen+2);
  f17de->Draw("SAME");
  f17f->SetLineWidth(5); f17f->SetLineStyle(kDashed);
  f17f->SetLineColor(kRed);
  f17f->Draw("SAME");

  // 2018
  f18a->SetLineWidth(5); f18a->SetLineStyle(kSolid);
  f18a->SetLineColor(kBlack);
  f18a->Draw("SAME");
  f18b->SetLineWidth(5); f18b->SetLineStyle(kSolid);
  f18b->SetLineColor(kBlue);
  f18b->Draw("SAME");
  f18c->SetLineWidth(5); f18c->SetLineStyle(kSolid);
  f18c->SetLineColor(kGreen+2);
  f18c->Draw("SAME");
  f18d->SetLineWidth(5); f18d->SetLineStyle(kSolid);
  f18d->SetLineColor(kRed);
  f18d->Draw("SAME");

  TLatex *tex = new TLatex();
  tex->SetTextSize(0.045);
  tex->SetTextAlign(31); // adjust right
  tex->SetTextColor(kGray);
  //tex->DrawLatex(3000,1.014,"16BCD+EF");
  tex->SetTextColor(kBlack);
  //tex->DrawLatex(3000,1.005,"16GH");
  tex->DrawLatex(3000,0.998,"18A");
  tex->SetTextColor(kGray);
  //tex->DrawLatex(3000,0.9975,"17B");
  tex->SetTextColor(kBlue);
  //tex->DrawLatex(3000,0.9935,"17C");
  tex->DrawLatex(3000,0.9935,"18B");
  tex->SetTextColor(kGreen+2);
  //tex->DrawLatex(3000,0.985,"17DE");
  tex->DrawLatex(3000,0.988,"18C");
  tex->SetTextColor(kRed);
  //tex->DrawLatex(3000,0.977,"17F");
  tex->DrawLatex(3000,0.974,"17D");


  c1->SaveAs("pdf/compareL3Res.pdf");
}
