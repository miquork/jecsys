#! /usr/bin/env python 
# compare dijet and Z+jets responses in various eta bins
# A.Karavdina 14.07.2017
from ROOT import *
import sys
import numpy
import os 


pathGL = "/afs/desy.de/user/k/karavdia/xxl/af-cms/GlobalFit_2016_dijet_Zjets_data/"
filelist = {'Dijet_BCD':'/Dijets_20170712/RunBCD/output/JEC_L2_Dijet_AK4PFchs_pythia8.root','Zeejet_BCD':'/Zjets_20170710/combination_Zee_BCD_remini_madgraph_NJ_2017-07-10.root','Zmumujet_BCD':'/Zjets_20170710/combination_Zmm_BCD_remini_madgraph_NJ_2017-07-10.root'}
channellist = {'Zeejet_BCD':'Zjets','Zmumujet_BCD':'Zjets','Dijet_BCD':'Dijet'}
alphalist = {'0.1':'a10','0.2':'a20','0.3':'a30'}
etalist = {'00_03':'00-026','03_05':'026-052','05_08':'052-078','08_10':'078-104','10_13':'130-147','13_15':'130-147','15_17':'147-165',
           '17_19':'165-193','19_22':'193-217','22_23':'217-232','23_25':'232-250','25_27':'250-265','27_29':'265-285','29_30':'285-296',
           '30_31':'296-313','31_35':'313-348','35_38':'348-383','38_52':'383-519'}
etalistnames = {'00_03':'0.00<|#eta|<0.26','03_05':'0.26<|#eta|<0.52','05_08':'0.52<|#eta|<0.78','08_10':'0.78<|#eta|<1.04','10_13':'1.30<|#eta|<1.47',
                '13_15':'1.30<|#eta|<1.47','15_17':'1.47<|#eta|<1.65','17_19':'1.65<|#eta|<1.93','19_22':'1.93<|#eta|<2.17','22_23':'2.17<|#eta|<2.32','23_25':'2.32<|#eta|<2.50',
                '25_27':'2.50<|#eta|<2.65','27_29':'2.65<|#eta|<2.85','29_30':'2.85<|#eta|<2.96',
                '30_31':'2.96<|#eta|<3.13','31_35':'3.13<|#eta|<3.48','35_38':'3.48<|#eta|<3.83','38_52':'3.83<|#eta|<5.19'}

hist_response_mpf_data = {}   
hist_response_mpf_dataDraw = {}
hist_response_mpf_mc = {}                
hist_response_mpf_ratio = {}                
hist_response_mpf_dataSum = {} 
cCompare = {}
legend = {}
for eta_key in etalist:
    for alpha_key in alphalist:
        cCompare[eta_key+alpha_key] = TCanvas("cCompare_Responses","Responses",800,600)                                                           
        legend[eta_key+alpha_key] = TLegend(.74,.75,.99,.95)
        hist_response_mpf_dataSum[eta_key+alpha_key] = TMultiGraph()
        gStyle.SetOptStat(0)                                                                                                                                      
        gStyle.SetOptTitle(0)
        i = 1
        for key_file in filelist:
            file = TFile(pathGL+filelist[key_file])
            if(channellist[key_file]=='Dijet'):
                #hist_response_mpf_data[key_file+alpha_key] = file.Get('data/eta'+etalist[eta_key]+'/mpfchs_dijet_'+alphalist[alpha_key]) 
                #hist_response_mpf_data[key_file+alpha_key] = file.Get('mc/eta'+etalist[eta_key]+'/mpfchs_dijet_'+alphalist[alpha_key])
                hist_response_mpf_data[key_file+alpha_key] = file.Get('data/eta'+etalist[eta_key]+'/ptchs_dijet_'+alphalist[alpha_key])  
                #hist_response_mpf_data[key_file+alpha_key] = file.Get('mc/eta'+etalist[eta_key]+'/ptchs_dijet_'+alphalist[alpha_key]) 

            if(channellist[key_file]=='Zjets'):
                #hist_response_mpf_data[key_file+alpha_key] = file.Get('Data_MPF_CHS_'+alphalist[alpha_key]+'_eta_'+eta_key+'_L1L2L3') 
                #hist_response_mpf_data[key_file+alpha_key] = file.Get('MC_MPF_CHS_'+alphalist[alpha_key]+'_eta_'+eta_key+'_L1L2L3') 
                hist_response_mpf_data[key_file+alpha_key] = file.Get('Data_PtBal_CHS_'+alphalist[alpha_key]+'_eta_'+eta_key+'_L1L2L3') 
                #hist_response_mpf_data[key_file+alpha_key] = file.Get('MC_PtBal_CHS_'+alphalist[alpha_key]+'_eta_'+eta_key+'_L1L2L3') 


            hist_response_mpf_data[key_file+alpha_key].SetMarkerStyle(20)
            hist_response_mpf_data[key_file+alpha_key].SetMarkerColor(i)
            hist_response_mpf_data[key_file+alpha_key].SetLineColor(i)
            hist_response_mpf_dataSum[eta_key+alpha_key].Add(hist_response_mpf_data[key_file+alpha_key])

            legend[eta_key+alpha_key].AddEntry(hist_response_mpf_data[key_file+alpha_key],key_file,"lp")     
            i = i+1
            #legend[eta_key+alpha_key].SetHeader('MPF, #alpha<'+alpha_key+', '+etalistnames[eta_key])
            legend[eta_key+alpha_key].SetHeader('PtBal, #alpha<'+alpha_key+', '+etalistnames[eta_key])
            hist_response_mpf_dataSum[eta_key+alpha_key].Draw('AP')
            hist_response_mpf_dataSum[eta_key+alpha_key].GetYaxis().SetTitle('response data')
           # hist_response_mpf_dataSum[eta_key+alpha_key].GetYaxis().SetTitle('response MC')
            hist_response_mpf_dataSum[eta_key+alpha_key].GetXaxis().SetTitle('p_{T}')
            hist_response_mpf_dataSum[eta_key+alpha_key].GetXaxis().SetRangeUser(10,10000)
            hist_response_mpf_dataSum[eta_key+alpha_key].GetYaxis().SetRangeUser(0.7,1.2)
            legend[eta_key+alpha_key].Draw()
                
    #   cCompare.SetLogx()
#        cCompare[eta_key+alpha_key].SaveAs('Dijet_Zjet_response_data_'+alphalist[alpha_key]+'_eta_'+eta_key+'.pdf')
#        cCompare[eta_key+alpha_key].SaveAs('Dijet_Zjet_response_MC_'+alphalist[alpha_key]+'_eta_'+eta_key+'.pdf')
#        cCompare[eta_key+alpha_key].SaveAs('Dijet_Zjet_response_MC_PtBal'+alphalist[alpha_key]+'_eta_'+eta_key+'.pdf')
        cCompare[eta_key+alpha_key].SaveAs('Dijet_Zjet_response_data_PtBal'+alphalist[alpha_key]+'_eta_'+eta_key+'.pdf')

