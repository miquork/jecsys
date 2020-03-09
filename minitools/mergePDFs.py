#! /usr/bin/python
import os
import collections
from string import ascii_uppercase

#IOV_list= ['A','B','C','D','ABC','ABCD']
IOV_list= ['A','B','D','ABCD']
#IOV_list= ['D','ABCD']
channel_list = ["DJ", "gam_zll", "MJDJ_gam_zll"]

##SampleDict=collections.OrderedDict()
##SampleDict["V16"]         ="Old"
##SampleDict["V17L1L2"]     ="NoClosure"
##SampleDict["V17Res"]="WithClosure"
##SampleDict["V17ResPtDep"]="WithClosurePtDepRel"
##SampleDict["V17ResPtDepFunc3"]="WithClosurePtDepRelFunc3"
##
##print("Current Working Directory ", os.getcwd())
##
##for iov in IOV_list:
##    for chann in channel_list:
##        inp = ""
##        shuffle = ""
##        out = ""
##        for name,path in SampleDict.iteritems():
##            idx = SampleDict.keys().index(name)
##            inp+= "{}=L1L2L3Res_2Par_{}Uncertainties/CollectL2Output_{}_PtBalMPF/jecslides_FineEta_Autumn18_V16_{}.pdf ".format(ascii_uppercase[idx], path, chann, iov)
##            shuffle+= "{} ".format( ascii_uppercase[idx])
##            out+= "_{}_{}{}".format( ascii_uppercase[idx],name, path)
##        print("pdftk {} shuffle {} output jecslides_Autumn18_V16_{}_{}{}.pdf".format(inp, shuffle, iov, chann, out))
##        os.system("pdftk {} shuffle {} output jecslides_Autumn18_V16_{}_{}_{}.pdf".format(inp, shuffle, iov, chann, out))
#
#
##SampleDict=collections.OrderedDict()
##SampleDict["JERConst"]         =""
##SampleDict["JERFunc1"]="Func1"
##SampleDict["JERFunc3"]="Func3"
##
##print("Current Working Directory ", os.getcwd())
##
##### 
##for iov in IOV_list:
##    for chann in channel_list:
##        inp = ""
##        shuffle = ""
##        out = ""
##        for name,path in SampleDict.iteritems():
##            idx = SampleDict.keys().index(name)
##            inp+= "{}=L1L2L3Res_2Par_{}NoClosureUncertainties/CollectL2Output_{}_PtBalMPF/jecslides_FineEta_Autumn18_V16_{}.pdf ".format(ascii_uppercase[idx], path, chann, iov)
##            shuffle+= "{} ".format( ascii_uppercase[idx])
##            out+= "_{}_{}{}".format( ascii_uppercase[idx],name, path)
##        print("pdftk {} shuffle {} output jecslides_Autumn18_V16_{}_{}{}.pdf".format(inp, shuffle, iov, chann, out))
##        os.system("pdftk {} shuffle {} output jecslides_Autumn18_V16_{}_{}_{}.pdf".format(inp, shuffle, iov, chann, out))
##
##
##for iov in IOV_list:
##    for name,path in SampleDict.iteritems():
##        inp = ""
##        shuffle = ""
##        out = ""
##        for idx,chann in enumerate(channel_list):
##            #idx = SampleDict.keys().index(name)
##            inp+= "{}=L1L2L3Res_2Par_{}NoClosureUncertainties/CollectL2Output_{}_PtBalMPF/jecslides_FineEta_Autumn18_V16_{}.pdf ".format(ascii_uppercase[idx], path, chann, iov)
##            shuffle+= "{} ".format( ascii_uppercase[idx])
##            out+= "_{}_{}".format( ascii_uppercase[idx],chann)
##        print("pdftk {} shuffle {} output jecslides_Autumn18_V16_L1L2L3Res{}_{}{}.pdf".format(inp, shuffle, iov, name, out))
##        os.system("pdftk {} shuffle {} output jecslides_Autumn18_V16_L1L2L3Res_{}_{}_{}.pdf".format(inp, shuffle, iov, name, out))
##        
##
#SampleDict=collections.OrderedDict()
#SampleDict["JERConst"]         =""
#SampleDict["JERFunc1"]="_Func1"
#SampleDict["JERFunc3"]="_Func3"
#
#for chann in channel_list:
#    for name,path in SampleDict.iteritems():
#        inp = ""
#        shuffle = ""
#        out = ""
#        for idx,iov in enumerate(IOV_list):
#            #            idx = SampleDict.keys().index(name)
#            inp+= "{}=L1L2_2Par{}/CollectL2Output_{}_PtBalMPF/jecslides_FineEta_Autumn18_V16_{}.pdf ".format(ascii_uppercase[idx], path, chann, iov)
#            shuffle+= "{} ".format( ascii_uppercase[idx])
#            out+= "_{}_{}".format( ascii_uppercase[idx],iov)
#        print("pdftk {} shuffle {} output jecslides_Autumn18_V16_L1L2_{}_{}{}.pdf".format(inp, shuffle, name, chann, out))
#        os.system("pdftk {} shuffle {} output jecslides_Autumn18_V16_L1L2_{}_{}_{}.pdf".format(inp, shuffle, name, chann, out))
#
#
#
##for iov in IOV_list:
##    for chann in channel_list:
##        inp = ""
##        shuffle = ""
##        out = ""
##        for name,path in SampleDict.iteritems():
##            idx = SampleDict.keys().index(name)
##            inp+= "{}=L1L2_2Par{}/CollectL2Output_{}_PtBalMPF/jecslides_FineEta_Autumn18_V16_{}.pdf ".format(ascii_uppercase[idx], path, chann, iov)
##            shuffle+= "{} ".format( ascii_uppercase[idx])
##            out+= "_{}_{}".format( ascii_uppercase[idx],name)
##        print("pdftk {} shuffle {} output jecslides_Autumn18_V16_L1L2_{}_{}{}.pdf".format(inp, shuffle, iov, chann, out))
##        os.system("pdftk {} shuffle {} output jecslides_Autumn18_V16_L1L2_{}_{}_{}.pdf".format(inp, shuffle, iov, chann, out))
##
##
##for iov in IOV_list:
##    for name,path in SampleDict.iteritems():
##        inp = ""
##        shuffle = ""
##        out = ""
##        for idx,chann in enumerate(channel_list):
##            #idx = SampleDict.keys().index(name)
##            inp+= "{}=L1L2_2Par{}/CollectL2Output_{}_PtBalMPF/jecslides_FineEta_Autumn18_V16_{}.pdf ".format(ascii_uppercase[idx], path, chann, iov)
##            shuffle+= "{} ".format( ascii_uppercase[idx])
##            out+= "_{}_{}".format( ascii_uppercase[idx],chann)
##        print("pdftk {} shuffle {} output jecslides_Autumn18_V16_L1L2_{}_{}{}.pdf".format(inp, shuffle, iov, name, out))
##        os.system("pdftk {} shuffle {} output jecslides_Autumn18_V16_L1L2_{}_{}_{}.pdf".format(inp, shuffle, iov, name, out))
##
#


IOV_list= ['A','B','C','D','ABC','ABCD']
channel_list = ["DJ", "gam_zll", "MJDJ_gam_zll"]

for chann in channel_list:
    inp = ""
    shuffle = ""
    out = ""
    for idx, iov in enumerate(IOV_list):
#        for name,path in SampleDict.iteritems():
        inp+= "{}=Autumn18_V17_L1L2_2Par/CollectL2Output_{}_PtBalMPF/jecslides_FineEta_Autumn18_V16_{}.pdf ".format(ascii_uppercase[idx], chann, iov)
        shuffle+= "{} ".format( ascii_uppercase[idx])
        out+= "_{}_{}".format( ascii_uppercase[idx],iov)
    print("pdftk {} shuffle {} output jecslides_Autumn18_V17_L1L2_{}_{}{}.pdf".format(inp, shuffle, iov, chann, out))
    os.system("pdftk {} shuffle {} output jecslides_Autumn18_V17_L1L2_{}_{}_{}.pdf".format(inp, shuffle, iov, chann, out))






for iov in IOV_list:
    inp = ""
    shuffle = ""
    out = ""
    for idx,chann in enumerate(channel_list):
        inp+= "{}=Autumn18_V17_L1L2L3Res_2Par/CollectL2Output_{}_PtBalMPF/jecslides_FineEta_Autumn18_V16_{}.pdf ".format(ascii_uppercase[idx], chann, iov)
        shuffle+= "{} ".format( ascii_uppercase[idx])
        out+= "_{}_{}".format( ascii_uppercase[idx],chann)
    print("pdftk {} shuffle {} output jecslides_Autumn18_V17_fixedIOV_{}_{}.pdf".format(inp, shuffle, iov, out))
    os.system("pdftk {} shuffle {} output jecslides_Autumn18_V17_fixedIOV_{}_{}.pdf".format(inp, shuffle, iov, out))


for chann in channel_list:
    inp = ""
    shuffle = ""
    out = ""
    for idx,iov in enumerate(IOV_list):
        inp+= "{}=Autumn18_V17_L1L2L3Res_2Par/CollectL2Output_{}_PtBalMPF/jecslides_FineEta_Autumn18_V16_{}.pdf ".format(ascii_uppercase[idx], chann, iov)
        shuffle+= "{} ".format( ascii_uppercase[idx])
        out+= "_{}_{}".format( ascii_uppercase[idx],chann)
    print("pdftk {} shuffle {} output jecslides_Autumn18_V17_fixedChann_{}_{}.pdf".format(inp, shuffle, chann, out))
    os.system("pdftk {} shuffle {} output jecslides_Autumn18_V17_fixedChann_{}_{}.pdf".format(inp, shuffle, chann, out))
    

        

    
