// Purpose: Produce plots of <rho> vs <Npv> for Zero bias
//          These are used as input in RhoVsNPV.C to determine UE
// run with 'root -l -b -q minitools/drawZBrhoVsNpv.C'

void drawZBrhoVsNpv() {

  TDirectory *curdir = gDirectory;

  TChain *tma = new TChain("T","T"); // APV / preVFP / BCDEF
  tma->Add("rootfiles/L1OffsetUL16_VFP/MCUL16/Offset_MC_UL2016_FlatPU0to75_106X_preVFP_v8_total.root");

  /*
  TChain *tmb = new TChain("T","T"); // non-APV / postVFP / GH
  tmb->Add("rootfiles/L1OffsetUL16_VFP/MCUL16/Offset_DATA_UL2016_postVFP_total.root");
  */
  TChain *tda = new TChain("T","T"); // APV / preVFP / BCDEF
  tda->Add("rootfiles/L1OffsetUL16_VFP/DataUL16/Offset_Data_UL2016B_preVFP_reduced.root");
  tda->Add("rootfiles/L1OffsetUL16_VFP/DataUL16/Offset_Data_UL2016C_preVFP_reduced.root");
  tda->Add("rootfiles/L1OffsetUL16_VFP/DataUL16/Offset_Data_UL2016D_preVFP_reduced.root");
  tda->Add("rootfiles/L1OffsetUL16_VFP/DataUL16/Offset_Data_UL2016E_preVFP_reduced.root");
  tda->Add("rootfiles/L1OffsetUL16_VFP/DataUL16/Offset_Data_UL2016F_preVFP_reduced.root");

  TChain *tdb = new TChain("T","T"); // non-APV / postVFP / GH
  tdb->Add("rootfiles/L1OffsetUL16_VFP/DataUL16/Offset_Data_UL2016F_postVFP_reduced.root");
  tdb->Add("rootfiles/L1OffsetUL16_VFP/DataUL16/Offset_Data_UL2016G_postVFP_reduced.root");
  tdb->Add("rootfiles/L1OffsetUL16_VFP/DataUL16/Offset_Data_UL2016H_postVFP_reduced.root");


  TFile *fout = new TFile("rootfiles/L1RC_RhoVsNpv.root","RECREATE");

  tma->Draw("rho:mu>>p_rho_nPU_MC16BCDEF(100,0,100)","","prof");
  tma->Draw("nPV:mu>>p_nPV_nPU_MC16BCDEF(100,0,100)","","prof");
  //tmb->Draw("rho:mu>>p_rho_nPU_MC16GH(100,0,100)","","prof");
  //tmb->Draw("nPV:mu>>p_nPV_nPU_MC16GH(100,0,100)","","prof");

  tda->Draw("rho:mu>>p_rho_nPU_UL16BCDEF(100,0,100)","","prof");
  tda->Draw("nPV:mu>>p_nPV_nPU_UL16BCDEF(100,0,100)","","prof");
  tdb->Draw("rho:mu>>p_rho_nPU_UL16GH(100,0,100)","","prof");
  tdb->Draw("nPV:mu>>p_nPV_nPU_UL16GH(100,0,100)","","prof");

  fout->Write();
} // drawZBRhoVsNpv
