import ROOT
import os
from array import array
import math
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import numpy as np

x_pt = array('d',
             [8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
              97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 362, 430, 507,
              592, 686, 790, 905, 1032, 1172, 1327, 1497, 1684, 1890, 2116,
              2366, 2640, 2941, 3273, 3637, 4037, 4477, 4961, 5492, 6076, 7000])
np_pt = np.array(x_pt[3:])
print np_pt

x_eta = array('d',
    [-5.4,-5.0,-4.4,-4,-3.5,-3,-2.8,-2.6,-2.4,-2.2,-2.0,
     -1.8,-1.6,-1.4,-1.2,-1.0, -0.8,-0.6,-0.4,-0.2,0.,
     0.2,0.4,0.6,0.8,1.0,1.2,1.4,
     1.6,1.8,2.0,2.2,2.4,2.6,2.8,3,3.5,4,4.4,5.0,5.4])

np_eta = np.array(x_eta[1:-1])
print(np_eta)

ROOT.gROOT.ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc++");
ROOT.gROOT.ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectorParameters.cc++");
ROOT.gROOT.ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc++");
ROOT.gROOT.ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc++");
#//
ROOT.gROOT.ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc++");
ROOT.gROOT.ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc++");

verbose = True

class jetmetUncertaintiesReader():
    def __init__(self, era="2018", globalTag="Autumn18_V18_DATA", jesUncertainties = [ "All" ], jetType = "AK4PFchs"):
        self.era=era
        self.globalTag=globalTag
        self.jetType=jetType
        self.jesUncertainties = jesUncertainties
        if not "AK4PFchs" in jetType : 
            raise ValueError("ERROR: Invalid jet type = '%s'!" % jetType)

        # read jet energy scale (JES) uncertainties
        self.jesInputFilePath = "CondFormats/JetMETObjects/data"
        self.jesUncertaintyInputFileName = globalTag + "_UncertaintySources_" + jetType + ".txt"
        print("use for uncertainties: {}/{}".format(self.jesInputFilePath ,self.jesUncertaintyInputFileName))

        # read all uncertainty source names from the loaded file
        if jesUncertainties[0] == "All":
            with open(self.jesInputFilePath+'/'+self.jesUncertaintyInputFileName) as f:
                lines = f.read().split("\n")
                sources = filter(lambda x: x.startswith("[") and x.endswith("]"), lines)
                sources = map(lambda x: x[1:-1], sources)
                self.jesUncertainties = sources
        if verbose:
            print("List of contained sources:")
            self.jesUncertainties.sort()
            print(self.jesUncertainties)
        
        self.jesUncertainty = {} 

        for jesUncertainty in self.jesUncertainties:
            jesUncertainty_label = jesUncertainty
            if jesUncertainty == 'Total' and len(self.jesUncertainties) == 1:
                jesUncertainty_label = ''
            pars = ROOT.JetCorrectorParameters(os.path.join(self.jesInputFilePath, self.jesUncertaintyInputFileName),jesUncertainty_label)
            self.jesUncertainty[jesUncertainty_label] = ROOT.JetCorrectionUncertainty(pars)    
        
        self.KeepSingleSources = ['FlavorQCD','RelativeBal']#,'Total']
        
        self.DictOfCommonSources = {}
        self.DictOfCommonSources["Absolute"] = ['AbsoluteMPFBias', 'AbsoluteScale', 'Fragmentation', 'PileUpDataMC', 'PileUpPtRef', 'RelativeFSR', 'SinglePionECAL', 'SinglePionHCAL']
        self.DictOfCommonSources["BBEC1"] = ['PileUpPtBB','PileUpPtEC1','RelativePtBB']
        self.DictOfCommonSources["EC2"] = ['PileUpPtEC2']
        self.DictOfCommonSources["HF"] = ['PileUpPtHF','RelativeJERHF','RelativePtHF']
        #[AbsoluteCommonSources, BBEC1CommonSources, EC2CommonSources, HFCommonSources ]
        
        self.DictOfPerYearSources = {}
        PerYearSourceSuffix = self.era
        #PerYearSourceSuffix = "YEAR"
        self.DictOfPerYearSources["Absolute_"+PerYearSourceSuffix] =['AbsoluteStat', 'RelativeStatFSR','TimePtEta']
        self.DictOfPerYearSources["BBEC1_"+PerYearSourceSuffix] = ['RelativeJEREC1','RelativePtEC1','RelativeStatEC']
        self.DictOfPerYearSources["EC2_"+PerYearSourceSuffix] = ['RelativeJEREC2','RelativePtEC2']
        self.DictOfPerYearSources["HF_"+PerYearSourceSuffix] = ['RelativeStatHF']
        self.DictOfPerYearSources["RelativeSample_"+PerYearSourceSuffix] = ['RelativeSample']

        self.AllSourcesConsideredForMerging = list(self.KeepSingleSources)
        for v in self.DictOfCommonSources.values() :
            self.AllSourcesConsideredForMerging+=v
        for v in self.DictOfPerYearSources.values() :
            self.AllSourcesConsideredForMerging+=v
        #print(self.AllSourcesConsideredForMerging)
        #print(self.DictOfPerYearSources)
        #print(self.DictOfCommonSources)
        
        
    def printTestValuesForSource(self):
        for jesUnc in self.jesUncertainties:
            # (cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetCorUncertainties )
            self.jesUncertainty[jesUnc].setJetPt(pt)
            self.jesUncertainty[jesUnc].setJetEta(eta)
            delta = self.jesUncertainty[jesUnc].getUncertainty(True)
            print(jesUnc, delta)
    def writeOutSources(self):
        f = open(self.jesUncertaintyInputFileName,"w") 
        f.write('#Uncertainty sources for {}_DATA_{}\n'.format(self.globalTag,self.jetType))
        for jesUnc in self.jesUncertainties:
            f.write("[{}]\n".format(jesUnc))
            f.write("{1 JetEta 1 JetPt \"\" Correction JECSource}\n")
            for ieta in range(len(x_eta)-1):
                etamin = x_eta[ieta]
                etamax = x_eta[ieta+1]
                eta= 0.5*(etamin+etamax)
                #print(etamin, etamax)
                etaline = "{:1.1f} {:1.1f} {:d} ".format(etamin,etamax,(len(x_pt)-1)*3)
                #                etaline=""
                for ipt in range(len(x_pt)-1):
                    pt = 0.5*(x_pt[ipt]+x_pt[ipt+1])
                    self.jesUncertainty[jesUnc].setJetPt(pt)
                    self.jesUncertainty[jesUnc].setJetEta(eta)
                    err = float(self.jesUncertainty[jesUnc].getUncertainty(True))
                    etaline+= "{:1.1f} {:1.4f} {:1.4f} ".format(pt, err, err)
                    #                    etaline+= ("%1.1f  " % {pt})
                    #                    etaline+= ("%1.1f %1.4f %1.4f " % {pt, err, err})
                    #                    print(jesUnc, delta)
                f.write(etaline+"\n") 

        f.close()
        
        
        
    def combinedSourceSection(self,sourcesToCombine=[],combinedName=""):
        if combinedName=="" or len(sourcesToCombine)==0:
            raise ValueError("ERROR: No way to combine this... '%s' with [possibly] empty source list" %combinedName)
        if verbose: print("combined name: {} will merge ".format(combinedName),sourcesToCombine)
        collectLine='[{}]\n'.format(combinedName)
        collectLine+="{1 JetEta 1 JetPt \"\" Correction JECSource}\n"
        for ieta in range(len(x_eta)-1):
            etamin = x_eta[ieta]
            etamax = x_eta[ieta+1]
            eta= 0.5*(etamin+etamax)
            etaline = "{:1.1f} {:1.1f} {:d} ".format(etamin,etamax,(len(x_pt)-1)*3)
            for ipt in range(len(x_pt)-1):
                pt = 0.5*(x_pt[ipt]+x_pt[ipt+1])
                err2=0
                for source in sourcesToCombine:
                    if source in self.jesUncertainties:
                        self.jesUncertainty[source].setJetPt(pt)
                        self.jesUncertainty[source].setJetEta(eta)
                        err2 += math.pow(float(self.jesUncertainty[source].getUncertainty(True)),2)
                    else:
                        raise ValueError('"{}" not in self.jesUncertainties of era {}. Probably this era needs special treatment (because source not present)'.format(source,self.era))
                err = math.sqrt(err2)
                etaline+= "{:1.1f} {:1.4f} {:1.4f} ".format(pt, err, err)
            collectLine+=(etaline+"\n") 
        return collectLine

    def writeCorrelationGrouping(self):
        print("Now write out...: Regrouped_"+self.jesUncertaintyInputFileName)
        f = open("Regrouped_"+self.jesUncertaintyInputFileName,"w")
        if verbose:
            print("Now copying ", self.KeepSingleSources)
        for s in self.KeepSingleSources: f.write(self.combinedSourceSection([s],s))
        for k,v in self.DictOfCommonSources.items():
            if verbose: print(k,v)
            f.write(self.combinedSourceSection(v,k))
        for k,v in self.DictOfPerYearSources.items():
            f.write(self.combinedSourceSection(v,k))
        #Add Total (as a copy) and calculate it:
        f.write(self.combinedSourceSection(["Total"],"Total"))
        #f.write(self.combinedSourceSection(self.AllSourcesConsideredForMerging,"TotalRecalculated"))

            

pt=50.
eta=2.

def jesHelper(pt,eta,jes):
    #print(pt, eta)
    jes.setJetPt(pt)
    jes.setJetEta(eta)
    return jes.getUncertainty(True)




jmeCorrections2018 = jetmetUncertaintiesReader(era="2018", globalTag="Autumn18_V19_MC")
#jmeCorrections2018.printTestValuesForSource()
#jmeCorrections2018.writeOutSources()
#print jmeCorrections2018.combinedSourceSection(["AbsoluteMPFBias","AbsoluteScale","Fragmentation","PileUpDataMC"],"Absolute")
jmeCorrections2018.writeCorrelationGrouping()
#
#
jmeCorrections2017 = jetmetUncertaintiesReader(era="2017", globalTag="Fall17_17Nov2017_V32_MC")
jmeCorrections2017.writeCorrelationGrouping()
#                    
jmeCorrections2016 = jetmetUncertaintiesReader(era="2016", globalTag="Summer16_07Aug2017_V11_MC")
jmeCorrections2016.writeCorrelationGrouping()

jmeCorrections2018RG = jetmetUncertaintiesReader(era="2018RG", globalTag="Regrouped_Autumn18_V19_MC")
jmeCorrections2017RG = jetmetUncertaintiesReader(era="2017RG", globalTag="Regrouped_Fall17_17Nov2017_V32_MC")
jmeCorrections2016RG = jetmetUncertaintiesReader(era="2016RG", globalTag="Regrouped_Summer16_07Aug2017_V11_MC")


plot_etas = [0.0,2.7,3.5]
plot_pts = [30.,100.,1000.]

print np_pt


def writeSummaryFig(jeslist=[], sourceList=[], sourceListShortName="", perJESPlots = True, perSourcePlots=True, Prefix=""):
    KeepArrays = {}
    for jes in jeslist:
        for s in sourceList:
            KeepArrays[(jes,s)] = [
#                ["$p_{T}$", np.array([jesHelper(30.,x,jes.jesUncertainty[s]) for x in np_eta])],
#                ["", np.array([jesHelper(100.,x,jes.jesUncertainty[s]) for x in np_eta])],
#                ["", np.array([jesHelper(1000.,x,jes.jesUncertainty[s]) for x in np_eta])],
#                ["eta", np.array([jesHelper(x,0.,jes.jesUncertainty[s]) for x in np_pt])],
#                ["eta", np.array([jesHelper(x,2.7,jes.jesUncertainty[s]) for x in np_pt])],
#                ["eta", np.array([jesHelper(x,3.5,jes.jesUncertainty[s]) for x in np_pt])],
                ["$p_{T}=30\\,GeV$", np.array([jesHelper(30.,x,jes.jesUncertainty[s]) for x in np_eta])],
                ["$p_{T}=100\\,GeV$", np.array([jesHelper(100.,x,jes.jesUncertainty[s]) for x in np_eta])],
                ["$p_{T}=1000\\,GeV$", np.array([jesHelper(1000.,x,jes.jesUncertainty[s]) for x in np_eta])],
                ["$\\eta^{jet}=0$", np.array([jesHelper(x,0.,jes.jesUncertainty[s]) for x in np_pt])],
                ["$\\eta^{jet}=2.7$", np.array([jesHelper(x,2.7,jes.jesUncertainty[s]) for x in np_pt])],
                ["$\\eta^{jet}=3.5$", np.array([jesHelper(x,3.5,jes.jesUncertainty[s]) for x in np_pt])],
            ]
            


        #per jes: plot sourcelist
        if perJESPlots:
            fig, axs = plt.subplots(2, 3, sharey=True, sharex=False)#, clear=True)
            for s in sourceList:
                for i in range(0,3):
                    axs[0,i].plot(np_eta,KeepArrays[(jes,s)][i][1],label=s)
                    axs[0,i].set_title(KeepArrays[(jes,s)][i][0],loc='left')
                    axs[0,i].set_xlabel("$\\eta^{jet}$")#,loc="right")
                    axs[1,i].semilogx(np_pt,KeepArrays[(jes,s)][i+3][1],label=s)
                    axs[1,i].set_title(KeepArrays[(jes,s)][i+3][0],loc='left')
                    axs[1,i].set_xlabel("$p_{T}$")
            niceName = ""
            #FigTitle = sourceListShortName
            if sourceListShortName=="":
                niceName="_".join(sourceList)
            else:
                niceName=sourceListShortName
            #FigTitle+="_".join(sourceList)
            handles, labels = axs[0,0].get_legend_handles_labels()
            fig.legend(handles, labels, loc="lower center", ncol=3)
            fig.suptitle("JES era: "+jes.era.replace("_","\_"))
            plt.tight_layout()
            #plt.subplots_adjust(bottom=0.15)
            plt.subplots_adjust(bottom=0.2,top=0.9)
            plt.savefig("pdf/UncSources/"+Prefix+"PerJes_"+jes.era+"_"+niceName+'.pdf')
            plt.savefig("pdf/UncSources/"+Prefix+"PerJes_"+jes.era+"_"+niceName+'.png')
            plt.close(fig)

    #per source: plot jeslist
    if perSourcePlots:
        for s in sourceList:
            fig, axs = plt.subplots(2, 3, sharey=True, sharex=False)#, clear=True)
            jesEras=""
            for jes in jeslist:
                jesEras+=jes.era+"_"
                for i in range(0,3):
                    axs[0,i].plot(np_eta,KeepArrays[(jes,s)][i][1],label=jes.era)
                    axs[0,i].set_title(KeepArrays[(jes,s)][i][0],loc='left')
                    axs[0,i].set_xlabel("$\\eta^{jet}$")#,loc="right")
                    axs[1,i].semilogx(np_pt,KeepArrays[(jes,s)][i+3][1],label=jes.era)
                    axs[1,i].set_title(KeepArrays[(jes,s)][i+3][0],loc='left')
                    axs[1,i].set_xlabel("$p_{T}$")
            handles, labels = axs[0,0].get_legend_handles_labels()
            fig.legend(handles, labels, loc="lower center", ncol=4)
            fig.suptitle("Source: "+s.replace("_","\_"))
            plt.tight_layout()
            plt.subplots_adjust(bottom=0.2,top=0.9)
            plt.savefig("pdf/UncSources/"+Prefix+"PerSource_"+s+"_"+jesEras+'.pdf')
            plt.savefig("pdf/UncSources/"+Prefix+"PerSource_"+s+"_"+jesEras+'.png')
            plt.close(fig)
        

#writeSummaryFig([jmeCorrections2016, jmeCorrections2017, jmeCorrections2018, ],["Total","SubTotalRelative"])



for k,v in jmeCorrections2018.DictOfCommonSources.items():
    writeSummaryFig([jmeCorrections2016, jmeCorrections2017, jmeCorrections2018, ],v,sourceListShortName="CommonSources"+k,perSourcePlots=False)


for k,v in jmeCorrections2018.DictOfPerYearSources.items():
    writeSummaryFig([jmeCorrections2016, jmeCorrections2017, jmeCorrections2018, ],v,sourceListShortName="PerYearSources"+k,perSourcePlots=False)

writeSummaryFig([jmeCorrections2016, jmeCorrections2017, jmeCorrections2018, ],jmeCorrections2018.KeepSingleSources,sourceListShortName="SingleSources",perSourcePlots=False)

writeSummaryFig([jmeCorrections2016, jmeCorrections2017, jmeCorrections2018, ],jmeCorrections2018.AllSourcesConsideredForMerging,perJESPlots=False)


#writeSummaryFig([jmeCorrections2016RG, jmeCorrections2017RG, jmeCorrections2018RG, ],jmeCorrections2018RG.jesUncertainties,perJESPlots=False,Prefix="Regrouped_")
#writeSummaryFig([jmeCorrections2016, jmeCorrections2017, jmeCorrections2018, ],['AbsoluteMPFBias', 'AbsoluteScale', 'Fragmentation', 'PileUpDataMC', 'PileUpPtRef', 'RelativeFSR', 'SinglePionECAL', 'SinglePionHCAL'],sourceListShortName="AbsCommonSources",perSourcePlots=False)



