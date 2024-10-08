import matplotlib
matplotlib.use('Agg')

import ROOT as r
from optparse import OptionParser
from time import sleep
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import mplhep as hep
from root_numpy import hist2array
import ctypes
from pathlib import Path
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import os 
import misc.plotRPFWithUnc as plotRPFWithUnc
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
matplotlib.use('Agg')
r.gROOT.SetBatch(True)
r.gStyle.SetOptFit(111)
matplotlib.rc('axes.formatter', use_mathtext=True)

def plotDeltaNLLComp(baseDir,odir,outFile,xRange=[0.,30.0],yRange=[0.,15.],signalFactor=10.,extraText="HZy coupling"):
    rss     = []
    dNLLs   = []
    hbbBR   = 0.58

    for fName in ["higgsCombinefreezeAll.MultiDimFit.mH125.root","higgsCombinefloatAll.MultiDimFit.mH125.root","higgsCombinefreezeAllPlusRPF.MultiDimFit.mH125.root"]:
        dNLL    = []
        rs      = []
        f       = r.TFile.Open(baseDir+"/"+fName)
        ttree   = f.Get("limit")
        for i,event in enumerate(ttree):
            if i==0:
                continue
            rs.append(event.r*signalFactor/hbbBR)
            #rs.append(event.r)
            dNLL.append(event.deltaNLL*2)
        dNLLs.append(dNLL)
        rss.append(rs)
        f.Close()

    plt.style.use([hep.style.CMS])
    f, ax = plt.subplots()
    plt.plot(rss[1],dNLLs[1],'bo-',label="Total uncert.")
    plt.plot(rss[0],dNLLs[0],'ro-',label="Frozen constrained nuis.")
    plt.plot(rss[2],dNLLs[2],'go-',label="Frozen constrained nuis.+ RPF")

    hep.cms.text("WiP",loc=0)
    lumiText = "138 $fb^{-1} (13 TeV)$"    
    hep.cms.lumitext(lumiText)
    ax.set_xlim(xRange)
    ax.set_ylim(yRange)

    plt.xlabel(r"$\sigma (pp\rightarrow H\gamma)$ [fb]",horizontalalignment='right', x=1.0)
    #plt.xlabel("Signal strength [r]",horizontalalignment='right', x=1.0)
    plt.ylabel("2$\Delta$ NLL",horizontalalignment='right', y=1.0)
    ax.yaxis.set_tick_params(which='minor', left=False)    
    ax.yaxis.set_tick_params(which='minor', right=False)    


    ax.axhline(y=1.0, color="grey",linestyle="--")
    ax.axhline(y=4.0, color="grey",linestyle="--")
    ax.axhline(y=9.0, color="grey",linestyle="--")

    plt.text(1.0,1.2,"1$\sigma$",color="grey")
    plt.text(1.0,4.2,"2$\sigma$",color="grey")
    plt.text(1.0,9.2,"3$\sigma$",color="grey")
    plt.text(1.0,14,extraText,color="black")

    plt.legend(loc="best")

    plt.savefig("test.png", bbox_inches='tight')

    print("Saving "+odir+"{0}.pdf".format(outFile))
    plt.savefig(odir+"/{0}.pdf".format(outFile), bbox_inches='tight')
    plt.savefig(odir+"/{0}.png".format(outFile), bbox_inches='tight')
    plt.cla()
    plt.clf()


def plotDeltaNLL(baseDir,odir,outFile,xRange=[0.,30.0],yRange=[0.,15.],signalFactor=10.,extraText="HZy coupling"):
    f       = r.TFile.Open(baseDir+"/higgsCombineTest.MultiDimFit.mH125.root")
    ttree   = f.Get("limit")
    rs      = []
    dNLL    = []
    hbbBR   = 0.58
    for i,event in enumerate(ttree):
        if i==0:
            continue
        rs.append(event.r*signalFactor/hbbBR)
        #rs.append(event.r)
        dNLL.append(event.deltaNLL*2)

    plt.style.use([hep.style.CMS])
    f, ax = plt.subplots()
    plt.plot(rs,dNLL,'ro-')

    hep.cms.text("WiP",loc=0)
    lumiText = "138 $fb^{-1} (13 TeV)$"    
    hep.cms.lumitext(lumiText)
    ax.set_xlim(xRange)
    ax.set_ylim(yRange)

    plt.xlabel(r"$\sigma (pp\rightarrow H\gamma)$ [fb]",horizontalalignment='right', x=1.0)
    #plt.xlabel("Signal strength [r]",horizontalalignment='right', x=1.0)
    plt.ylabel("2$\Delta$ NLL",horizontalalignment='right', y=1.0)
    ax.yaxis.set_tick_params(which='minor', left=False)    
    ax.yaxis.set_tick_params(which='minor', right=False)    


    ax.axhline(y=1.0, color="grey",linestyle="--")
    ax.axhline(y=4.0, color="grey",linestyle="--")
    ax.axhline(y=9.0, color="grey",linestyle="--")

    plt.text(1.0,1.2,"1$\sigma$",color="grey")
    plt.text(1.0,4.2,"2$\sigma$",color="grey")
    plt.text(1.0,9.2,"3$\sigma$",color="grey")
    plt.text(1.0,14,extraText,color="black")

    plt.savefig("test.png", bbox_inches='tight')

    print("Saving "+odir+"{0}.pdf".format(outFile))
    plt.savefig(odir+"/{0}.pdf".format(outFile), bbox_inches='tight')
    plt.savefig(odir+"/{0}.png".format(outFile), bbox_inches='tight')
    plt.cla()
    plt.clf()


def plotRPFMC(tplFile,odir,passTag="M",failTag="F",xRange=[60,200],yTitle="$R_{M/F}$",rebinX=1,CRFlag=False):
    f               = r.TFile.Open(tplFile)
    foutName        = "R{0}{1}_MC".format(passTag,failTag)
    histos          = []
    histosErrs      = []
    projectionsLo   = [1,2]
    projectionsHi   = [1,-1]
    maxRpf          = 0.
    if(CRFlag):
        hPass        = f.Get("QCD_H_m_pT_{0}__nominal".format(passTag)).ProjectionX("hPass_temp")
        hFail        = f.Get("QCD_H_m_pT_{0}__nominal".format(failTag)).ProjectionX("hFail_temp")

        hPass.RebinX(rebinX)
        hFail.RebinX(rebinX)


        hRpf    = hPass.Clone("hRPF")
        hRpf.Divide(hFail)
        hRpf.Scale(1000)
        maxRpf  = max(maxRpf,hRpf.GetMaximum())

        hRpfErr = []
        for i in range(1,hRpf.GetNbinsX()+1):
            hRpfErr.append(hRpf.GetBinError(i))
        hRpf, edges = hist2array(hRpf,return_edges=True)
        histos.append(hRpf)
        histosErrs.append(hRpfErr)

    else:
        for i in range(2):
            projectionLo = projectionsLo[i]
            projectionHi = projectionsHi[i]
            hPass        = f.Get("GJets_H_m_pT_{0}__nominal".format(passTag)).ProjectionX("hPass_temp",projectionLo,projectionHi)
            hFail        = f.Get("GJets_H_m_pT_{0}__nominal".format(failTag)).ProjectionX("hFail_temp",projectionLo,projectionHi)

            hPass.RebinX(rebinX)
            hFail.RebinX(rebinX)


            hRpf    = hPass.Clone("hRPF")
            hRpf.Divide(hFail)
            hRpf.Scale(1000)
            maxRpf  = max(maxRpf,hRpf.GetMaximum())

            hRpfErr = []
            for i in range(1,hRpf.GetNbinsX()+1):
                hRpfErr.append(hRpf.GetBinError(i))
            hRpf, edges = hist2array(hRpf,return_edges=True)
            histos.append(hRpf)
            histosErrs.append(hRpfErr)


    f.Close()

    plt.style.use([hep.style.CMS])
    f, ax = plt.subplots()
    if(CRFlag):
        hep.histplot(histos,yerr=histosErrs, bins=edges[0],label=["$p_T > 450$ GeV"])
    else:
        hep.histplot(histos,yerr=histosErrs, bins=edges[0],label=["$300 < p_T < 400$ GeV", "$p_T > 400$ GeV"])
    plt.legend()

    hep.cms.text("Work in progress",loc=0)
    
    ax.set_xlim(xRange)
    ax.set_ylim([0.,maxRpf*1.3])

    plt.xlabel("$M_{PNet}$ (GeV)",horizontalalignment='right', x=1.0)
    plt.ylabel(yTitle+" x $10^{3}$",horizontalalignment='right', y=1.0)
    #ax.yaxis.set_tick_params(which='minor', left=False)    
    #ax.yaxis.set_tick_params(which='minor', right=False)    

    print("Saving "+odir+"{0}.pdf".format(foutName))
    plt.savefig(odir+"/{0}.pdf".format(foutName), bbox_inches='tight')
    plt.savefig(odir+"/{0}.png".format(foutName), bbox_inches='tight')
    plt.cla()
    plt.clf()

def plotRPF(postfitShapesFile,odir,qcdTag,passTag="M",failTag="F",xRange=[60,200],yTitle="$R_{M/F}$"):
    hPass2D = get2DPostfitPlot(postfitShapesFile,qcdTag,passTag)
    hFail2D = get2DPostfitPlot(postfitShapesFile,"qcd",failTag)
    nBins   = hPass2D.GetNbinsY()
    foutName = "R{0}{1}".format(passTag,failTag)
    maxRpf  = 0.
    histos  = []
    if(nBins==2):
        legendFlag  = True
        labels      = ["300<$p_T$<400 GeV","$p_T$>400 GeV"]
    else:
        labels      = ["_dummy"]
        legendFlag  = False

    for i in range(1,nBins+1):
        hPass = hPass2D.ProjectionX("hPass_temp",i,i)
        hFail = hFail2D.ProjectionX("hFail_temp",i,i)
        hRpf  = hPass.Clone("hRPF")
        hRpf.Divide(hFail)
        hRpf.Scale(1000)
        maxRpf = max(maxRpf,hRpf.GetMaximum())
        hRpf, edges = hist2array(hRpf,return_edges=True)
        histos.append(hRpf)

    plt.style.use([hep.style.CMS])
    f, ax = plt.subplots()
    hep.histplot(histos, bins=edges[0],label=labels)
    if legendFlag:
        plt.legend()

    hep.cms.text("Work in progress",loc=0)
    
    ax.set_xlim(xRange)
    ax.set_ylim([0.,maxRpf*1.3])

    plt.xlabel("$M_{PNet}$ (GeV)",horizontalalignment='right', x=1.0)
    plt.ylabel(yTitle+" x $10^{3}$",horizontalalignment='right', y=1.0)
    #ax.yaxis.set_tick_params(which='minor', left=False)    
    #ax.yaxis.set_tick_params(which='minor', right=False)    

    print("Saving "+odir+"{0}.pdf".format(foutName))
    plt.savefig(odir+"/{0}.pdf".format(foutName), bbox_inches='tight')
    plt.savefig(odir+"/{0}.png".format(foutName), bbox_inches='tight')
    plt.cla()
    plt.clf()


def plotTriggerEff(hPass,hTotal,year,luminosity,outFile,xlabel="Photon $p_{T}$ (GeV)",ylabel="Trigger efficiency / 10 GeV"):
    TEff   = r.TEfficiency(hPass,hTotal)
    effs   = []
    errsUp = []
    errsDn = []
    binsX  = []
    for i in range(1,hPass.GetNbinsX()+1):
        binX  = hPass.GetBinCenter(i)
        eff   = TEff.GetEfficiency(i)
        errUp = TEff.GetEfficiencyErrorUp(i)
        errDn = TEff.GetEfficiencyErrorLow(i)

        binsX.append(binX)
        effs.append(eff)
        errsDn.append(errDn)
        errsUp.append(errUp)


    plt.style.use([hep.style.CMS])
    f, ax = plt.subplots()
    plt.errorbar(binsX, effs, yerr=[errsDn,errsUp], fmt='o',color="black",label="{0} Data trigger efficiency".format(year))

    plt.xlabel(xlabel, horizontalalignment='right', x=1.0)
    plt.ylabel(ylabel,horizontalalignment='right', y=1.0)

    lumiText = luminosity + " $fb^{-1}\ (13 TeV)$"
    hep.cms.lumitext(text=lumiText, ax=ax, fontname=None, fontsize=None)
    hep.cms.text("WiP",loc=0)
    ax.set_ylim([0.5,1.02])
    plt.legend(loc="best")#loc = 'best'
    plt.tight_layout()

    print("Saving {0}".format(outFile))
    plt.savefig(outFile)
    plt.savefig(outFile.replace("pdf","png"))
    plt.clf()

def mergeLoMassBins(hist,edges):
    #Merges first 3 pairs of mass bins into 10 GeV bins to combat ParticleNet mass regression sculpting
    #Converts hist to Events / GeV

    binsToRemove    = [1,3,5]

    newEdges        = np.delete(edges,binsToRemove)
    newHist         = np.delete(hist,binsToRemove)

    for binIdx in binsToRemove:
        newBinIdx   = int((binIdx-1)/2)
        newHist[newBinIdx]+=hist[binIdx]

    for i in range(len(newHist)):
        binWidth    = newEdges[i+1]-newEdges[i]
        newHist[i]  = newHist[i]/binWidth

    return newHist, [newEdges]


def mergeMassBinsErrors(histErrs,edges):
    #Merges errors in first 3 pairs of mass bins into 10 GeV bins to combat ParticleNet mass regression sculpting
    #Converts hist to Events / GeV
    #edges should already be merged (output of mergeLoMassBins(hist,edges))
    binsToRemove    = [1,3,5]
    
    newErrs         = np.delete(histErrs,binsToRemove)
    
    for binIdx in binsToRemove:
        newBinIdx   = int((binIdx-1)/2)
        newErrs[newBinIdx]=np.sqrt(newErrs[newBinIdx]**2+histErrs[binIdx]**2)
    
    for i in range(len(newErrs)):
        binWidth    = edges[i+1]-edges[i]
        newErrs[i]  = newErrs[i]/binWidth

    return newErrs
    
def blindHiggsMass(hist):
    for i in range(1,hist.GetNbinsX()+1):
        binCenter = hist.GetBinCenter(i)
        if binCenter>110 and binCenter<145:
            hist.SetBinContent(i,0)
    return hist


def errListFromTH1(hist):
    errs = []
    for i in range(1,hist.GetNbinsX()+1):
        errs.append(hist.GetBinError(i))
    return errs

def plotVarStack(data,var,outFile,xTitle="",yTitle="",yRange=[],xRange=[],log=True,rebinX=1,luminosity="36.3",mergeMassBins=False,blind=True,kFactor=1.25,projection="",ratioRange=[0.,2.]):
    histos  = []
    labels  = []
    edges   = []
    colors  = []
    histosData = []#we're assuming only one data_obs dataset
    edgesData  = []#it's still kept in array (with one element) to be similar to other processes
    labelsData = []
    histosSig  = []#we're assuming only one signal dataset
    edgesSig   = []#it's still kept in array (with one element) to be similar to other processes
    labelsSig  = []
    colorsSig  = []

    totalMCYields = []
    totalMCErrors = []

    data = sorted(data.items(), key=lambda x: x[1]["order"])
    for sample, sample_cfg in data:
        if("singlephoton" in sample.lower() or "data" in sample.lower() or "jetht" in sample.lower()):
            dataFlag = True
        else:
            dataFlag = False
        #These are workaround because we now have all signals in the same file
        #The histogram names also have awkwardly placed signal coupling name
        if(sample=="HgammaHZy"):
            sigFlag  = True
            hName    = "{0}_{1}".format(sample.replace("HZy",""),var.replace("__nominal","__weight_ghza2__nominal"))
        elif(sample=="HgammaHyy"):
            sigFlag  = True
            hName    = "{0}_{1}".format(sample.replace("Hyy",""),var.replace("__nominal","__weight_gha2__nominal"))
        else:
            sigFlag  = False
            hName    = "{0}_{1}".format(sample,var)

        tempFile = r.TFile.Open(sample_cfg["file"])
        #print("Opening ", hName)
        h = tempFile.Get(hName)

        #Normalization adjustments
        if(sample=="SMHiggs"):
            sigFlag  = True
            h.Scale(10.)
        if(sample=="HgammaHZy"): 
            h.Scale(10./(0.58*1.33*10**4))#10fb/(BR(H->bb)*xsec(pp->Hgamma for g2hza=1))
        if(sample=="HgammaHyy"): 
            h.Scale(10./(0.58*0.735*10**4))#10fb/(BR(H->bb)*xsec(pp->Hgamma for g2ha=1))

        if(projection=="x"):
            h = h.ProjectionX(h.GetName()+"_x")
        elif(projection=="y"):
            h = h.ProjectionX(h.GetName()+"_y")

        if("GJets" in sample or "QCD" in sample):
            h.Scale(kFactor)
        h.RebinX(rebinX)
        if(dataFlag and blind):
            h = blindHiggsMass(h)

        hist, edges = hist2array(h,return_edges=True)
        histErrs    = errListFromTH1(h)
        if(mergeMassBins):
            hist, edges = mergeLoMassBins(hist,edges[0])
            histErrs    = mergeMassBinsErrors(histErrs,edges[0])
        if dataFlag:
            histosData.append(hist)
            dataMax = max(hist)
            edgesData.append(edges[0])
            labelsData.append(sample_cfg["label"])
            continue
        elif sigFlag:
            histosSig.append(hist)
            edgesSig.append(edges[0])
            labelsSig.append(sample_cfg["label"])
            colorsSig.append(sample_cfg["color"])
            continue            
        else:
            histos.append(hist)
            edges.append(edges[0])
            labels.append(sample_cfg["label"])
            colors.append(sample_cfg["color"])
            if(sample_cfg["label"]=="ttbar+Gamma"):
                labels[-1]=r"t$\bar{t}+Gamma$"#json restrictions workaround

            if len(totalMCYields)==0:
                totalMCYields = hist
                totalMCErrors = histErrs
            else:
                totalMCYields = [x + y for x, y in zip(totalMCYields, hist)]
                totalMCErrors = [np.sqrt(x**2 + y**2) for x, y in zip(totalMCErrors, histErrs)]

    plt.style.use([hep.style.CMS])
    f, axs = plt.subplots(2,1, sharex=True, sharey=False,gridspec_kw={'height_ratios': [4, 1],'hspace': 0.07})
    axs = axs.flatten()
    plt.sca(axs[0])

    centresData = (edgesData[0][:-1] + edgesData[0][1:])/2.
    if(mergeMassBins):
        errorsData = []
        for i, histCount in enumerate(histosData[0]):
            width       = edges[0][i+1]-edges[0][i]
            errorData   = np.sqrt(histCount*width)/width
            errorsData.append(errorData)
    else:
        errorsData  = np.sqrt(histosData[0])
    xerrorsData = []
    #if(mergeMassBins):
    for i in range(len(edges[0])-1):
        xerror = (edges[0][i+1]-edges[0][i])/2.
        xerrorsData.append(xerror)
    #--------------------------#
    hep.histplot(histos,edges[0],stack=True,ax=axs[0],label=labels,linewidth=1,histtype="fill",facecolor=colors,edgecolor='black')
    if(histosSig):
        hep.histplot(histosSig,edges[0],stack=False,ax=axs[0],label=labelsSig,linewidth=3,histtype="step",edgecolor=colorsSig)
    if mergeMassBins:
        plt.errorbar(centresData,histosData[0], yerr=errorsData, xerr=xerrorsData, fmt='o',color="k",label = labelsData[0])    
    else:
        plt.errorbar(centresData,histosData[0], yerr=errorsData, fmt='o',color="k",label = labelsData[0])    
    if(log):
        axs[0].set_yscale("log")
    axs[0].legend()
    axs[0].set_ylabel(yTitle)
    axs[1].set_xlabel(xTitle)
    axs[1].set_ylabel("Data/MC")
    #plt.xlabel(xTitle, horizontalalignment='right', x=1.0)
    plt.ylabel(yTitle,horizontalalignment='right', y=1.0)

    if(yRange):
        if(yRange[1]==None):
            yRange[1] = dataMax*2.0
        axs[0].set_ylim(yRange)
    else:
        axs[0].set_ylim([0,dataMax*2.0])
    if(xRange):
        axs[0].set_xlim(xRange)
    lumiText = luminosity + " $fb^{-1}\ (13 TeV)$"
    hep.cms.lumitext(text=lumiText, ax=axs[0], fontname=None, fontsize=None)
    hep.cms.text("WiP",loc=0)
    plt.legend(loc="best",ncol=2,handletextpad=0.3)#loc = 'best'
    
    #Ratio
    ratioVals = []
    ratioErrs = []
    for i in range(len(totalMCYields)):
        data      = histosData[0][i]
        dataSq    = data**2
        dataErrSq = errorsData[i]**2
        mc        = totalMCYields[i]
        mcSq      = mc**2
        mcErrSq   = totalMCErrors[i]**2

        ratioVal     = data/(mc+0.000001)#Protect division by zero
        relDataErrSq = dataErrSq/(dataSq+0.000001)
        relmcErrSq   = mcErrSq/(mcSq+0.000001)
        ratioErr     = ratioVal*np.sqrt(relDataErrSq+relmcErrSq)
        ratioVals.append(ratioVal)
        ratioErrs.append(ratioErr)

    plt.sca(axs[1])#switch to lower pad
    axs[1].set_ylim([ratioRange[0],ratioRange[1]])

    plt.errorbar(centresData,ratioVals,yerr=ratioErrs, xerr=xerrorsData, fmt='.',color="lightcoral") 
     
    axs[1].axhline(y=1.0, color="grey",linestyle="--")

    plt.tight_layout()
    print("Saving {0}".format(outFile))
    plt.savefig(outFile)
    plt.savefig(outFile.replace(".png",".pdf"))
    plt.cla()
    plt.clf()

def plotVarStackMC(data,var,outFile,xTitle="",yTitle="",yRange=[],xRange=[],log=True,rebinX=1,luminosity="36.3",mergeMassBins=False,projection=""):
    histos  = []
    labels  = []
    edges   = []
    colors  = []
    histosSig  = []#we're assuming only one signal dataset
    edgesSig   = []#it's still kept in array (with one element) to be similar to other processes
    labelsSig  = []
    colorsSig  = []

    data = sorted(data.items(), key=lambda x: x[1]["order"])#VERY HANDY, reads samples in order
    for sample, sample_cfg in data:
        if("singlephoton" in sample.lower() or "data" in sample.lower()):
            continue
        if("hgamma" in sample.lower()):
            sigFlag  = True
        else:
            sigFlag  = False

        tempFile = r.TFile.Open(sample_cfg["file"])
        #print("Opening ", sample_cfg["file"], "{0}_{1}".format(sample,var))
        h = tempFile.Get("{0}_{1}".format(sample,var))

        if(projection=="x"):
            h = h.ProjectionX(h.GetName()+"_x")
        elif(projection=="y"):
            h = h.ProjectionX(h.GetName()+"_y")

        h.RebinX(rebinX)
        hist, edges = hist2array(h,return_edges=True)
        if(mergeMassBins):
            hist, edges = mergeLoMassBins(hist,edges[0])
        if sigFlag:
            histosSig.append(hist)
            edgesSig.append(edges[0])
            labelsSig.append(sample_cfg["label"])
            colorsSig.append(sample_cfg["color"])
            continue            
        else:
            histos.append(hist)
            edges.append(edges[0])
            labels.append(sample_cfg["label"])
            colors.append(sample_cfg["color"])
            if(sample_cfg["label"]=="ttbar+Gamma"):
                labels[-1]=r"t$\bar{t}+Gamma$"#json restrictions workaround

    plt.style.use([hep.style.CMS])
    f, ax = plt.subplots()
    plt.sca(ax)

    if(mergeMassBins):
        for i in range(len(edges[0])-1):
            xerror = (edges[0][i+1]-edges[0][i])/2.
            xerrorsData.append(xerror)
    #--------------------------#
    hep.histplot(histos,edges[0],stack=True,ax=ax,label=labels,linewidth=1,histtype="fill",facecolor=colors,edgecolor='black')
    hep.histplot(histosSig,edges[0],stack=False,ax=ax,label=labelsSig,linewidth=3,histtype="step",edgecolor=colorsSig)
    if(log):
        ax.set_yscale("log")
    ax.legend()
    ax.set_ylabel(yTitle)
    ax.set_xlabel(xTitle)
    plt.ylabel(yTitle,horizontalalignment='right', y=1.0)
    
    if(yRange):
        ax.set_ylim(yRange)
    if(xRange):
        ax.set_xlim(xRange)
    lumiText = luminosity + " $fb^{-1}\ (13 TeV)$"
    hep.cms.lumitext(text=lumiText, ax=ax, fontname=None, fontsize=None)
    hep.cms.text("WiP",loc=0)
    plt.legend(loc="best",ncol=2,handletextpad=0.3)#loc = 'best'
    
    plt.tight_layout()
    print("Saving {0}".format(outFile))
    plt.savefig(outFile)
    plt.savefig(outFile.replace(".png",".pdf"))
    plt.cla()
    plt.clf()

def getPoissonErrors(hist,binWidthDivision=False):
    hist.SetBinErrorOption(1)

    #This is needed due to some nasty float precision inaccuracy causing some data content to be 0.9999998
    #The kPoisson error does not get calculated correctly in that case for some reason
    tempHist   = hist.Clone("tempHist_forErrs")
    tempHist.Reset()
    tempHist.SetBinErrorOption(1)

    errors_low = []
    errors_hi  = []
    for i in range(1,hist.GetNbinsX()+1):
        tempHist.SetBinContent(i,int(round(hist.GetBinContent(i))))
        #print(int(hist.GetBinContent(i)),tempHist.GetBinErrorLow(i),tempHist.GetBinErrorUp(i))
        if(binWidthDivision):
            errors_low.append(tempHist.GetBinErrorLow(i)/tempHist.GetBinWidth(i))
            errors_hi.append(tempHist.GetBinErrorUp(i)/tempHist.GetBinWidth(i))
        else:
            errors_low.append(tempHist.GetBinErrorLow(i))
            errors_hi.append(tempHist.GetBinErrorUp(i))

    return [errors_low,errors_hi]

def rebinHisto(hModel,hToRebin,name,scale=1.0):
    hRes = hModel.Clone(name)
    hRes.Reset()
    xaxis = hToRebin.GetXaxis()
    yaxis = hToRebin.GetYaxis()
    xaxis_re = hRes.GetXaxis()
    yaxis_re = hRes.GetYaxis()
    for i in range(1,hToRebin.GetNbinsX()+1):
        for j in range(1,hToRebin.GetNbinsY()+1):
            x = xaxis.GetBinCenter(i)
            y = yaxis.GetBinCenter(j)
            i_re = xaxis_re.FindBin(x)
            j_re = yaxis_re.FindBin(y)
            value = hToRebin.GetBinContent(i,j)
            if(value<0.):
                value = 0.
            err = hToRebin.GetBinError(i,j)
            err_re = np.sqrt(hRes.GetBinError(i_re,j_re)*hRes.GetBinError(i_re,j_re)+err*err)
            hRes.Fill(x,y,value)
            hRes.SetBinError(i_re,j_re,err_re)
    hRes.Scale(scale)
    hRes.SetDirectory(0)
    return hRes

def getUncBand(totalHistos):
    yLo = []
    yUp = []
    for i in range(1,totalHistos.GetNbinsX()+1):
        errLo  = totalHistos.GetBinErrorLow(i)
        errUp  = totalHistos.GetBinErrorUp(i)
        mcPred = totalHistos.GetBinContent(i)
        yLo.append(mcPred-errLo)
        yUp.append(mcPred+errUp)
    return np.array(yLo), np.array(yUp)

def calculatePull(hData,dataErrors,hTotBkg,uncBand):
    pulls = []
    for i,dataYield in enumerate(hData):
        mcYield     = hTotBkg[i]
        diff        = dataYield-mcYield
        dataErr     = np.sqrt(dataYield)
        if(dataYield>=mcYield):
            dataErr = dataErrors[0][i]#ErrorLo
            mcErr   = uncBand[1][i]-mcYield#ErrorUp
        elif(dataYield==0):
            pull = 0.
            pulls.append(pull)
            continue
        else:
            dataErr = dataErrors[1][i]#ErrorUp
            mcErr   = uncBand[0][i]-mcYield#ErrorLo

        sigma =  np.sqrt(dataErr*dataErr+mcErr*mcErr)
        pull        = diff/sigma
        pulls.append(pull)

    return np.array(pulls)

def divideByBinWidth(h,divideFlag=True):
    hist, edges = hist2array(h,return_edges=True)

    if not divideFlag:
        return hist,edges

    else:
        newHist     = []
        for i in range(len(hist)):
            binWidth    = edges[0][i+1]-edges[0][i]
            newHist.append(hist[i]/binWidth)

        return newHist,edges

def plotShapesWithRatio(hData,hMC,uncBand,labelsMC,colorsMC,xlabel,outputFile,xRange=[],yRange=[],projectionText="",binWidthDivision=True):
    dataErrors      = getPoissonErrors(hData,binWidthDivision=binWidthDivision)
    hData, edges    = divideByBinWidth(hData,divideFlag=binWidthDivision)

    centresData     = (edges[0][:-1] + edges[0][1:])/2.#Bin centres
    xerrorsData     = []

    for i in range(len(edges[0])-1):
        xerror = (edges[0][i+1]-edges[0][i])/2.
        xerrorsData.append(xerror)

    histosMC        = []
    for h in hMC:
        histosMC.append(divideByBinWidth(h,divideFlag=binWidthDivision)[0])


    plt.style.use([hep.style.CMS])
    #matplotlib.rcParams.update({'font.size': 30})
    f, axs = plt.subplots(2,1, sharex=True, sharey=False,gridspec_kw={'height_ratios': [4, 1],'hspace': 0.08})
    axs = axs.flatten()
    plt.sca(axs[0])

    if any("H+" in label for label in labelsMC):
        #last entry of histosMC = bkg-only
        #next-to-lasta entry of histosMC = bkg+sig
        hep.histplot(histosMC[:-3],edges[0],stack=True,label = labelsMC[:-1], histtype="fill",facecolor=colorsMC[:-1])
        hep.histplot(histosMC[-2],edges[0],stack=True,label = labelsMC[-1], histtype="step")
    else:
        hep.histplot(histosMC[:-2],edges[0],stack=True,label = labelsMC, histtype="fill",facecolor=colorsMC)
    plt.errorbar(centresData,hData, yerr=dataErrors, xerr=xerrorsData, fmt='o',color="k",label = "Data",zorder=11)

    def calcRatioAndErrors(hData,hMC,dataErrs,MCUncBand):
        ratioVals=[]
        ratioErrs=[]
        for i in range(len(hData)):
            data      = hData[i]
            mc        = hMC[i]
            dataSq    = data**2
            mcSq      = mc**2
            if data<mc:
                dataErr = dataErrs[1][i] #dataErrs[0]/[1] = errors low/high
                MCErr = MCUncBand[0][i]-mc
            else:
                dataErr = dataErrs[0][i]
                MCErr = mc-MCUncBand[1][i]
            
            mcErrSq   = MCErr**2
            dataErrSq = dataErr**2

            ratioVal     = data/(mc+0.000001)#Protect division by zero
            relDataErrSq = dataErrSq/(dataSq+0.000001)
            relmcErrSq   = mcErrSq/(mcSq+0.000001)
            ratioErr     = ratioVal*np.sqrt(relDataErrSq+relmcErrSq)
            ratioVals.append(ratioVal)
            ratioErrs.append(ratioErr)

        return ratioVals, ratioErrs

    for i in range(len(uncBand[0])):
        binWidth            = edges[0][i+1]-edges[0][i]
        if(binWidthDivision):
            uncBand[0][i]       = uncBand[0][i]/binWidth
            uncBand[1][i]       = uncBand[1][i]/binWidth
    
    ratioVals, ratioErrs = calcRatioAndErrors(hData,histosMC[-1],dataErrors,uncBand)

    uncBandLow = np.append(uncBand[0],[0],axis=0)
    uncBandHi  = np.append(uncBand[1],[0],axis=0)#Hack to get the last bin uncertainty to plot, since we're using step="post"

    plt.fill_between(edges[0],uncBandLow,uncBandHi,facecolor="none", hatch="xxx", edgecolor="grey", linewidth=0.0,step="post",zorder=10)

    axs[0].legend()
    if(binWidthDivision):
        plt.ylabel("Events/GeV",horizontalalignment='right', y=1.0)
    else:
        plt.ylabel("Events/bin",horizontalalignment='right', y=1.0)
    axs[1].set_ylabel("Ratio")

    if(xRange):
        axs[0].set_xlim(xRange)
    if(yRange):
        axs[0].set_ylim(yRange)
    else:
        yMaximum = max(hData)*2.0
        axs[0].set_ylim([0.,yMaximum])

    if("16APV" in outputFile):
        lumi = 19.5
    elif("16" in outputFile):
        lumi = 16.8
    elif("17" in outputFile):
        lumi = 41.5
    elif("18" in outputFile):
        lumi = 59.8
    else:
        lumi = 138

    lumiText = str(lumi)+ " $fb^{-1} (13 TeV)$"    
    hep.cms.lumitext(lumiText)
    hep.cms.text("WiP",loc=0)
    plt.legend(loc=1,ncol=2)

    if(projectionText):
        plt.text(0.60, 0.60, projectionText, horizontalalignment='center',verticalalignment='center',transform=axs[0].transAxes)


    plt.sca(axs[1])#switch to lower pad
    axs[1].axhline(y=1.0, xmin=0, xmax=1, color="grey",linestyle="--",alpha=0.5)
    if "CR" in outputFile.split("/")[-1]:
        axs[1].set_ylim([0.95,1.05])
    else:
        axs[1].set_ylim([0.0,2.0])
    plt.xlabel(xlabel,horizontalalignment='right', x=1.0)

    plt.errorbar(centresData,ratioVals,yerr=ratioErrs, xerr=xerrorsData, fmt='.',color="lightcoral") 

    print("Saving ", outputFile)
    #plt.tight_layout()
    plt.savefig(outputFile,bbox_inches="tight")
    plt.savefig(outputFile.replace("png","pdf"))

    plt.clf()
    plt.cla()

def plotShapes(hData,hMC,uncBand,labelsMC,colorsMC,xlabel,outputFile,xRange=[],yRange=[],projectionText="",binWidthDivision=True):
    dataErrors      = getPoissonErrors(hData,binWidthDivision=binWidthDivision)
    hData, edges    = divideByBinWidth(hData,divideFlag=binWidthDivision)

    centresData     = (edges[0][:-1] + edges[0][1:])/2.#Bin centres
    xerrorsData     = []

    for i in range(len(edges[0])-1):
        xerror = (edges[0][i+1]-edges[0][i])/2.
        xerrorsData.append(xerror)

    histosMC        = []
    for h in hMC:
        histosMC.append(divideByBinWidth(h,divideFlag=binWidthDivision)[0])


    plt.style.use([hep.style.CMS])
    #matplotlib.rcParams.update({'font.size': 30})
    f, axs = plt.subplots(2,1, sharex=True, sharey=False,gridspec_kw={'height_ratios': [4, 1],'hspace': 0.08})
    axs = axs.flatten()
    plt.sca(axs[0])

    hep.histplot(histosMC[:-1],edges[0],stack=True,ax=axs[0],label = labelsMC, histtype="fill",facecolor=colorsMC,edgecolor='black')
    plt.errorbar(centresData,hData, yerr=dataErrors, xerr=xerrorsData, fmt='o',color="k",label = "Data",zorder=11)

    for i in range(len(uncBand[0])):
        binWidth            = edges[0][i+1]-edges[0][i]
        if(binWidthDivision):
            uncBand[0][i]       = uncBand[0][i]/binWidth
            uncBand[1][i]       = uncBand[1][i]/binWidth
        else:
            uncBand[0][i]       = uncBand[0][i]
            uncBand[1][i]       = uncBand[1][i]

    uncBandLow = np.append(uncBand[0],[0],axis=0)
    uncBandHi  = np.append(uncBand[1],[0],axis=0)#Hack to get the last bin uncertainty to plot, since we're using step="post"

    plt.fill_between(edges[0],uncBandLow,uncBandHi,facecolor="none", hatch="xxx", edgecolor="grey", linewidth=0.0,step="post",zorder=10)

    axs[0].legend()
    if(binWidthDivision):
        plt.ylabel("Events/GeV",horizontalalignment='right', y=1.0)
    else:
        plt.ylabel("Events/bin",horizontalalignment='right', y=1.0)
    axs[1].set_ylabel("Pulls")

    if(xRange):
        axs[0].set_xlim(xRange)
    if(yRange):
        axs[0].set_ylim(yRange)
    else:
        yMaximum = max(hData)*2.0
        axs[0].set_ylim([0.,yMaximum])

    if("16APV" in outputFile):
        lumi = 19.5
    elif("16" in outputFile):
        lumi = 16.8
    elif("17" in outputFile):
        lumi = 41.5
    elif("18" in outputFile):
        lumi = 59.8
    else:
        lumi = 138

    lumiText = str(lumi)+ " $fb^{-1} (13 TeV)$"    
    hep.cms.lumitext(lumiText)
    hep.cms.text("WiP",loc=0)
    plt.legend(loc=1,ncol=2)

    if(projectionText):
        plt.text(0.60, 0.60, projectionText, horizontalalignment='center',verticalalignment='center',transform=axs[0].transAxes)


    plt.sca(axs[1])#switch to lower pad
    #axs[1].axhline(y=0.0, xmin=0, xmax=1, color="r")
    axs[1].axhline(y=1.0, xmin=0, xmax=1, color="grey",linestyle="--",alpha=0.5)
    axs[1].axhline(y=-1.0, xmin=0, xmax=1, color="grey",linestyle="--",alpha=0.5)
    axs[1].axhline(y=2.0, xmin=0, xmax=1, color="grey",linestyle="-.",alpha=0.5)
    axs[1].axhline(y=-2.0, xmin=0, xmax=1, color="grey",linestyle="-.",alpha=0.5)
    axs[1].set_ylim([-2.5,2.5])
    plt.xlabel(xlabel,horizontalalignment='right', x=1.0)

    pulls = calculatePull(hData,dataErrors,histosMC[-1],uncBand)
    if(max(pulls)>2.5 or min(pulls)<-2.5):
        axs[1].set_ylim([-4.,4.])
    else:
        axs[1].set_ylim([-2.5,2.5])
    hep.histplot(pulls,edges[0],ax=axs[1],linewidth=1,histtype="fill",facecolor="grey",edgecolor='black')

    print("Saving ", outputFile)
    #plt.tight_layout()
    plt.savefig(outputFile,bbox_inches="tight")
    plt.savefig(outputFile.replace("png","pdf"))

    plt.clf()
    plt.cla()

def plotShapesNoPull(hData,hMC,uncBand,labelsMC,colorsMC,xlabel,outputFile,xRange=[],yRange=[],projectionText="",binWidthDivision=True,twoSignal=False):
    dataErrors      = getPoissonErrors(hData,binWidthDivision=binWidthDivision)
    hData, edges    = divideByBinWidth(hData,divideFlag=binWidthDivision)

    centresData     = (edges[0][:-1] + edges[0][1:])/2.#Bin centres
    xerrorsData     = []

    for i in range(len(edges[0])-1):
        xerror = (edges[0][i+1]-edges[0][i])/2.
        xerrorsData.append(xerror)

    histosMC        = []
    for h in hMC:
        histosMC.append(divideByBinWidth(h,divideFlag=binWidthDivision)[0])

    if "CR_F" in outputFile:
        regionText="Untagged, fail"
    elif "CR_M" in outputFile:
        regionText="Untagged, medium"
    elif "CR_T" in outputFile:
        regionText="Untagged, tight"
    elif "F." in outputFile:
        regionText="$\gamma$-tagged, fail"
    elif "M." in outputFile:
        regionText="$\gamma$-tagged,\nmedium"
    elif "T." in outputFile:
        regionText="$\gamma$-tagged, tight"
    else:
        regionText=""

    plt.style.use([hep.style.CMS])
    matplotlib.rcParams.update({'font.size': 30})
    f, ax = plt.subplots()

    if any("gamma H" in label for label in labelsMC):
        if twoSignal:
            hep.histplot(histosMC[:-2],edges[0],stack=True,label = labelsMC[:-2], histtype="fill",facecolor=colorsMC[:-2])
            hep.histplot(histosMC[-1],edges[0],stack=False,label = "_nolegend_", histtype="step",color=colorsMC[-1])            
            hep.histplot(histosMC[-2],edges[0],stack=False,label = "_nolegend_", histtype="step",color=colorsMC[-2],linestyle="--")
            custom_handle1 = Patch(facecolor="none", edgecolor=colorsMC[-1], label=labelsMC[-1])  # Solid rectangle
            custom_handle2 = Patch(facecolor="none", edgecolor=colorsMC[-2], label=labelsMC[-2], linestyle='--')  # Dashed rectangle

        else:
            hep.histplot(histosMC[:-2],edges[0],stack=True,label = labelsMC[:-1], histtype="fill",facecolor=colorsMC[:-1])
            hep.histplot(histosMC[-1],edges[0],stack=True,label = "_nolegend_", histtype="step")
            custom_handle1 = Patch(facecolor="none", edgecolor=colorsMC[-1], label=labelsMC[-1])  # Solid rectangle
            
    else:
        hep.histplot(histosMC[:-1],edges[0],stack=True,label = labelsMC, histtype="fill",facecolor=colorsMC)
    plt.errorbar(centresData,hData, yerr=dataErrors, fmt='o',color="k",label = "Observed",zorder=11,markersize=10)

    # for i in range(len(uncBand[0])):
    #     binWidth            = edges[0][i+1]-edges[0][i]
    #     if(binWidthDivision):
    #         uncBand[0][i]       = uncBand[0][i]/binWidth
    #         uncBand[1][i]       = uncBand[1][i]/binWidth
    #     else:
    #         uncBand[0][i]       = uncBand[0][i]
    #         uncBand[1][i]       = uncBand[1][i]

    # uncBandLow = np.append(uncBand[0],[0],axis=0)
    # uncBandHi  = np.append(uncBand[1],[0],axis=0)#Hack to get the last bin uncertainty to plot, since we're using step="post"

    #plt.fill_between(edges[0],uncBandLow,uncBandHi,facecolor="none", hatch="xxx", edgecolor="grey", linewidth=0.0,step="post",zorder=10)

    ax.legend()
    if(binWidthDivision):
        #plt.ylabel("Events/GeV",horizontalalignment='right', y=1.0)
        plt.ylabel("Events/GeV",horizontalalignment='center', y=0.5)
    else:
        #plt.ylabel("Events / bin",horizontalalignment='right', y=1.0)
        plt.ylabel("Events / 10 GeV",horizontalalignment='center', y=0.5)

    if(xRange):
        ax.set_xlim(xRange)
    if(yRange):
        ax.set_ylim(yRange)
    else:
        yMaximum = max(hData)*1.4
        ax.set_ylim([0.,yMaximum])

    if("16APV" in outputFile):
        lumi = 19.5
    elif("16" in outputFile):
        lumi = 16.8
    elif("17" in outputFile):
        lumi = 41.5
    elif("18" in outputFile):
        lumi = 59.8
    else:
        lumi = 138

    lumiText = str(lumi)+ " $fb^{-1} (13 TeV)$"    
    hep.cms.lumitext(lumiText)
    #hep.cms.text("Preliminary",loc=0)
    hep.cms.text("",loc=0)
    if any("gamma H" in label for label in labelsMC):
        handles, labels = plt.gca().get_legend_handles_labels()
        if twoSignal:
            handles.extend([custom_handle1, custom_handle2])
            labels.extend([labelsMC[-1], labelsMC[-2]])
        else:
            handles.extend([custom_handle1])
            labels.extend([labelsMC[-1]])
        plt.legend(handles, labels,loc=1,ncol=1)
    else:
        plt.legend(loc=1,ncol=1)

    if(projectionText):
        plt.text(0.60, 0.60, projectionText, horizontalalignment='center',verticalalignment='center',transform=ax.transAxes)

    #plt.xlabel(xlabel,horizontalalignment='right', x=1.0)
    plt.xlabel(xlabel,horizontalalignment='center', x=0.5)

    plt.text(0.05, 0.95, regionText, transform=plt.gca().transAxes, verticalalignment='top')


    print("Saving ", outputFile)
    plt.tight_layout()
    plt.savefig(outputFile,bbox_inches="tight")
    plt.savefig(outputFile.replace("png","pdf"))

    plt.clf()
    plt.cla()


def printMCYields(data,region,year,debug=False):
    data = sorted(data.items(), key=lambda x: x[1]["order"])
    print("----------")
    print("Yields in {0} - {1}".format(year,region))
    for sample, sample_cfg in data:
        f       = r.TFile.Open(sample_cfg["file"])
        hist    = f.Get("{0}_H_m_pT_{1}__nominal".format(sample,region))
        if(debug):
            print(sample_cfg["file"], "{0}_H_m_{1}__nominal".format(sample,region))
        print("{0:10}\t{1:.0f}".format(sample,hist.Integral()))

def get2DPostfitPlot(file,process,region,prefit=False):
    f       = r.TFile.Open(file)
    if prefit:
        fitStatus = "prefit"
    else:
        fitStatus = "postfit"
    hLow    = f.Get("{0}_LOW_{2}/{1}".format(region,process,fitStatus))
    hSig    = f.Get("{0}_SIG_{2}/{1}".format(region,process,fitStatus))
    hHigh   = f.Get("{0}_HIGH_{2}/{1}".format(region,process,fitStatus))
    h2      = merge_low_sig_high(hLow,hSig,hHigh,hName="h2_{0}_{1}".format(process,region))
    h2.SetDirectory(0)
    return h2
    
def get_binning_x(hLow,hSig,hHigh):
    bins = []
    for i in range(1,hLow.GetNbinsX()+1):
        bins.append(hLow.GetXaxis().GetBinLowEdge(i))
    for i in range(1,hSig.GetNbinsX()+1):
        bins.append(hSig.GetXaxis().GetBinLowEdge(i))
    for i in range(1,hHigh.GetNbinsX()+2):#low edge of overflow is high edge of last bin
        bins.append(hHigh.GetXaxis().GetBinLowEdge(i))
    bins = np.array(bins,dtype='float64')
    return bins

def get_binning_y(hLow,hSig,hHigh):
    #histos should have same binning in Y
    bins = []
    for i in range(1,hLow.GetNbinsY()+2):
        bins.append(hLow.GetYaxis().GetBinLowEdge(i))
    bins = np.array(bins,dtype='float64')
    return bins

def merge_low_sig_high(hLow,hSig,hHigh,hName="temp"):
    n_x_low     = hLow.GetNbinsX()
    n_x_sig     = hSig.GetNbinsX()
    n_x_high    = hHigh.GetNbinsX()
    n_x         = n_x_low + n_x_sig + n_x_high
    n_y         = hLow.GetNbinsY()#assumes Y bins are the same
    bins_x      = get_binning_x(hLow,hSig,hHigh)
    bins_y      = get_binning_y(hLow,hSig,hHigh)
    h_res       = r.TH2F(hName,"",n_x,bins_x,n_y,bins_y)
    for i in range(1,n_x_low+1):
        for j in range(1,n_y+1):
            h_res.SetBinContent(i+0,j,hLow.GetBinContent(i,j))
            h_res.SetBinError(i+0,j,hLow.GetBinError(i,j))

    for i in range(1,n_x_sig+1):
        for j in range(1,n_y+1):
            h_res.SetBinContent(i+n_x_low,j,hSig.GetBinContent(i,j))
            h_res.SetBinError(i+n_x_low,j,hSig.GetBinError(i,j))

    for i in range(1,n_x_high+1):
        for j in range(1,n_y+1):
            h_res.SetBinContent(i+n_x_sig+n_x_low,j,hHigh.GetBinContent(i,j))
            h_res.SetBinError(i+n_x_sig+n_x_low,j,hHigh.GetBinError(i,j))
    return h_res

def plotPostfit(postfitShapesFile,region,odir,prefitTag=False,blind=True,binWidthDivision=True,debug=False,signal="HgammaHZy",plotSlices=False,paper=True):

    if(region=="CR_T" or region=="CR_F" or region=="CR_M"):
        CRFlag = True
    else:
        CRFlag = False

    if not paper:
        if(region=="pass" or region=="T" or region=="M"):
            labels              = ["Data","Non-resonant","W+Gamma","Z+Gamma","H+Gamma (HZy coupling)"]
            tags                = ["data_obs","qcd","WGamma","ZGamma",signal]
            colors              = ["black","gold","mistyrose","blue","red"]
        elif(CRFlag and region=="CR_F"):
            labels              = ["Data","Non-resonant","W+jets","Z+jets"]
            tags                = ["data_obs","qcd","WJets","ZJets"]
            colors              = ["black","gold","mistyrose","blue"]
        elif(CRFlag and region!="CR_F"):
            labels              = ["Data","Non-resonant","W+jets","Z+jets","SM H"]
            tags                = ["data_obs","qcd","WJets","ZJets","SMHiggs"]
            colors              = ["black","gold","mistyrose","blue","green"]
        else:
            labels              = ["Data","Non-resonant","W+Gamma","Z+Gamma"]
            tags                = ["data_obs","qcd","WGamma","ZGamma"]
            colors              = ["black","gold","mistyrose","blue"]
    else:
        if(region=="pass" or region=="T" or region=="M"):
            labels              = ["Observed","Non-resonant","W+$\gamma$","Z+$\gamma$",r"$\gamma H \; (\sigma \times B = 10 fb)$"]
            tags                = ["data_obs","qcd","WGamma","ZGamma",signal]
            colors              = ["black","#99ccff","lightgreen","#669966","#669966"]
        elif(CRFlag and region=="CR_F"):
            labels              = ["Observed","Non-resonant","W+jets","Z+jets"]
            tags                = ["data_obs","qcd","WJets","ZJets"]
            colors              = ["black","#99ccff","lightgreen","#669966"]
        elif(CRFlag and region!="CR_F"):
            labels              = ["Observed","Non-resonant","W+jets","Z+jets","SM H"]
            tags                = ["data_obs","qcd","WJets","ZJets","SMHiggs"]
            colors              = ["black","#99ccff","lightgreen","#669966","#ff9b9b"]
        else:
            labels              = ["Observed","Non-resonant","W+$\gamma$","Z+$\gamma$"]
            tags                = ["data_obs","qcd","WGamma","ZGamma"]
            colors              = ["black","#99ccff","lightgreen","#669966"]

    if(prefitTag):
        outFile = "prefit"
    else:
        outFile = "postfit"


    twoDShapes  = []

    dirRegion   = region
    if CRFlag:
        polyOrder = odir.split("/")[-2].split("CR")[0][-1]
    else:
        polyOrder = odir.split("/")[-2].split("SR")[0][-1]
    #polyOrder   = odir.split("/")[-2]

    #Merge sliced histograms
    for tag in tags:
        if(tag=="qcd" and (region=="pass" or region=="T" or region=="M")):
            tag     = tag+"_"+polyOrder
        elif(tag=="qcd" and (region=="CR_T" or region=="CR_M")):
            tag     = tag+"_CR_"+polyOrder
        if debug:
            print(postfitShapesFile,tag,dirRegion,prefitTag)
        twoDShape   = get2DPostfitPlot(postfitShapesFile,tag,dirRegion,prefit=prefitTag)
        twoDShapes.append(twoDShape)
        
        totalBkg  = get2DPostfitPlot(postfitShapesFile,"TotalBkg",dirRegion,prefit=prefitTag)
        totalProcs  = get2DPostfitPlot(postfitShapesFile,"TotalProcs".format(region),dirRegion,prefit=prefitTag)

    if CRFlag:
        xRange = [60,150]
    else:
        xRange = [60,200]


    if(plotSlices):
        #Plot mass
        for i in range(1,totalProcs.GetNbinsY()+1):
            projections         = []
            for j,twoDShape in enumerate(twoDShapes):
                proj        = twoDShape.ProjectionX(tags[j]+"_projx_{0}".format(i),i,i)
                if(tags[j]=="data_obs" and blind):
                    proj = blindHiggsMass(proj)
                projections.append(proj)
            totalProcs_proj     = totalProcs.ProjectionX("totalprocs_projx_{0}".format(i),i,i)
            totalBkg_proj     = totalBkg.ProjectionX("totalbkg_projx_{0}".format(i),i,i)
            projections.append(totalProcs_proj)
            projections.append(totalBkg_proj)
            uncBand_proj        = getUncBand(totalBkg_proj)
            projLowEdge         = int(totalProcs.GetYaxis().GetBinLowEdge(i))
            projUpEdge          = int(totalProcs.GetYaxis().GetBinUpEdge(i))
            projectionText      = "{0}".format(projLowEdge)+"<$p_{T}$<"+"{0} GeV".format(projUpEdge)
            if(projUpEdge==1000):
                projectionText  = "$p_{T}$>"+"{0} GeV".format(projLowEdge)

            if binWidthDivision:
                plotShapesWithRatio(projections[0],projections[1:],uncBand_proj,labels[1:],colors[1:],"$M_{PNet}_rescaled$ (GeV)","{0}/{1}_{2}_{3}.png".format(odir,outFile,region,i),xRange=xRange,projectionText=projectionText,binWidthDivision=binWidthDivision)
            else:
                plotShapesWithRatio(projections[0],projections[1:],uncBand_proj,labels[1:],colors[1:],"$M_{PNet}$ (GeV)","{0}/{1}_{2}_{3}.png".format(odir,outFile,region,i),xRange=xRange,projectionText=projectionText,binWidthDivision=binWidthDivision)

    projections         = []
    for j,twoDShape in enumerate(twoDShapes):
        proj            = twoDShape.ProjectionX(tags[j]+"_projx")
        if(tags[j]=="data_obs" and blind):
            proj = blindHiggsMass(proj)
        #if(tags[j]=="SMHiggs"):
        #    proj.Scale(10.)
        projections.append(proj)
    totalProcs_proj     = totalProcs.ProjectionX("totalprocs_projx")
    projections.append(totalProcs_proj)
    uncBand_proj        = getUncBand(totalProcs_proj)

    if binWidthDivision:
        plotName = "{0}/{1}_{2}_rescaled.png".format(odir,outFile,region)
    else:
        plotName = "{0}/{1}_{2}.png".format(odir,outFile,region)
    if binWidthDivision:
        if not paper:
            plotShapes(projections[0],projections[1:],uncBand_proj,labels[1:],colors[1:],"$M_{PNet}$ (GeV)",plotName,xRange=xRange,binWidthDivision=binWidthDivision)
        else:
            plotShapesNoPull(projections[0],projections[1:],uncBand_proj,labels[1:],colors[1:],"$M_{PNet}$ (GeV)",plotName,xRange=xRange,binWidthDivision=binWidthDivision)
    else:
        if not paper:
            plotShapes(projections[0],projections[1:],uncBand_proj,labels[1:],colors[1:],"$M_{PNet}$ (GeV)",plotName,xRange=xRange,binWidthDivision=binWidthDivision)
        else:
            #plotShapesNoPull(projections[0],projections[1:],uncBand_proj,labels[1:],colors[1:],"$M_{PNet}$ (GeV)",plotName,xRange=xRange,binWidthDivision=binWidthDivision)
            plotShapesWithRatioAndBand(projections[0],projections[1:],totalProcs_proj,uncBand_proj,labels[1:],colors[1:],"$M_{PNet}$ (GeV)",plotName,xRange=xRange,binWidthDivision=False,twoSignal=True)

    f = r.TFile.Open("{0}/{1}_{2}.root".format(odir,outFile,region),"RECREATE")
    f.cd()
    totalProcs.Write()
    for h in twoDShapes:
        h.Write()
    f.Close()

def plotPostfitTwoSignalsForPaper(postfitShapesFile,region,odir):
    signal = "HgammaHZy"
    signal_alt = "HgammaHyy"
    prefitTag = False

    if(region=="T" or region=="M"):
        labels              = ["Observed","Non-resonant","W+$\gamma$","Z+$\gamma$",r"$\gamma H \; (c_{z\gamma },\;\sigma \times B = 10 fb)$",r"$\gamma H \; (c_{\gamma \gamma },\;\sigma \times B = 10 fb)$"]
        tags                = ["data_obs","qcd","WGamma","ZGamma",signal,signal_alt]
        colors              = ["black","#99ccff","lightgreen","#669966","#FF0000","#FF0000"]
    else:
        print("ERROR: Region needs to be either M or T in plotPostfitTwoSignalsForPaper()")
        exit()

    outFile = "postfit"


    twoDShapes  = []

    dirRegion   = region
    polyOrder = odir.split("/")[-2].split("SR")[0][-1]
    

    #Merge sliced histograms
    for tag in tags:
        if(tag=="qcd"):
            tag     = tag+"_"+polyOrder
        if(tag==signal_alt):
            twoDShape   = get2DPostfitPlot(postfitShapesFile.replace("HZy","Hyy"),tag,dirRegion,prefit=prefitTag)
            twoDShapes.append(twoDShape)
        else:
            twoDShape   = get2DPostfitPlot(postfitShapesFile,tag,dirRegion,prefit=prefitTag)
            twoDShapes.append(twoDShape)        
    totalBkg  = get2DPostfitPlot(postfitShapesFile,"TotalBkg",dirRegion,prefit=prefitTag)
    totalBkg_X = totalBkg.ProjectionX("bkg_x")
    totalProcs  = get2DPostfitPlot(postfitShapesFile,"TotalProcs".format(region),dirRegion,prefit=prefitTag)

    xRange = [60,200]

    projections         = []
    for j,twoDShape in enumerate(twoDShapes):
        if("Hgamma" in tags[j]):
            proj = twoDShape.ProjectionX(tags[j]+"_projx")
            proj.Add(totalBkg_X)
        else:
            proj            = twoDShape.ProjectionX(tags[j]+"_projx")
        projections.append(proj)
    # totalProcs_proj     = totalProcs.ProjectionX("totalprocs_projx")
    # projections.append(totalProcs_proj)
    uncBand_proj        = [] #We don't plot unc band for paper

    plotName = "{0}/{1}_{2}.png".format(odir,outFile,region)
    plotShapesNoPull(projections[0],projections[1:],uncBand_proj,labels[1:],colors[1:],"$M_{PNet}$ (GeV)",plotName,xRange=xRange,binWidthDivision=False,twoSignal=True)

    f = r.TFile.Open("{0}/{1}_{2}.root".format(odir,outFile,region),"RECREATE")
    f.cd()
    totalProcs.Write()
    for h in twoDShapes:
        h.Write()
    f.Close()

def plotPostfitTwoSignalsForPaperWithRatio(postfitShapesFile,region,odir):
    signal = "HgammaHZy"
    signal_alt = "HgammaHyy"
    prefitTag = False

    if(region=="T" or region=="M"):
        labels              = ["Observed","Non-resonant","W+$\gamma$","Z+$\gamma$",r"$\gamma H \; (c_{z\gamma },\;\sigma \times B = 10 fb)$",r"$\gamma H \; (c_{\gamma \gamma },\;\sigma \times B = 10 fb)$"]
        tags                = ["data_obs","qcd","WGamma","ZGamma",signal,signal_alt]
        colors              = ["black","#99ccff","lightgreen","#669966","#FF0000","#FF0000"]
    else:
        print("ERROR: Region needs to be either M or T in plotPostfitTwoSignalsForPaper()")
        exit()

    outFile = "postfit"


    twoDShapes  = []

    dirRegion   = region
    polyOrder = odir.split("/")[-2].split("SR")[0][-1]
    

    #Merge sliced histograms
    for tag in tags:
        if(tag=="qcd"):
            tag     = tag+"_"+polyOrder
        if(tag==signal_alt):
            twoDShape   = get2DPostfitPlot(postfitShapesFile.replace("HZy","Hyy"),tag,dirRegion,prefit=prefitTag)
            twoDShapes.append(twoDShape)
        else:
            twoDShape   = get2DPostfitPlot(postfitShapesFile,tag,dirRegion,prefit=prefitTag)
            twoDShapes.append(twoDShape)        
    totalBkg  = get2DPostfitPlot(postfitShapesFile,"TotalBkg",dirRegion,prefit=prefitTag)
    totalBkg_X = totalBkg.ProjectionX("bkg_x")
    totalProcs  = get2DPostfitPlot(postfitShapesFile,"TotalProcs".format(region),dirRegion,prefit=prefitTag)

    xRange = [60,200]

    projections         = []
    for j,twoDShape in enumerate(twoDShapes):
        if("Hgamma" in tags[j]):
            proj = twoDShape.ProjectionX(tags[j]+"_projx")
            proj.Add(totalBkg_X)
        else:
            proj            = twoDShape.ProjectionX(tags[j]+"_projx")
        projections.append(proj)
    # totalProcs_proj     = totalProcs.ProjectionX("totalprocs_projx")
    # projections.append(totalProcs_proj)
    uncBand_proj        = getUncBand(totalBkg_X)

    plotName = "{0}/{1}_{2}.png".format(odir,outFile,region)
    plotShapesWithRatioAndBand(projections[0],projections[1:],totalBkg_X,uncBand_proj,labels[1:],colors[1:],"$M_{PNet}$ (GeV)",plotName,xRange=xRange,binWidthDivision=False,twoSignal=True)

    f = r.TFile.Open("{0}/{1}_{2}.root".format(odir,outFile,region),"RECREATE")
    f.cd()
    totalProcs.Write()
    for h in twoDShapes:
        h.Write()
    f.Close()

def plotShapesWithRatioAndBand(hData,hMC,bkg_proj,uncBand,labelsMC,colorsMC,xlabel,outputFile,xRange=[],yRange=[],projectionText="",binWidthDivision=True,twoSignal=False):
    dataErrors      = getPoissonErrors(hData,binWidthDivision=binWidthDivision)
    hData, edges    = divideByBinWidth(hData,divideFlag=binWidthDivision)

    centresData     = (edges[0][:-1] + edges[0][1:])/2.#Bin centres
    xerrorsData     = []

    for i in range(len(edges[0])-1):
        xerror = (edges[0][i+1]-edges[0][i])/2.
        xerrorsData.append(xerror)

    histosMC        = []
    for h in hMC:
        histosMC.append(divideByBinWidth(h,divideFlag=binWidthDivision)[0])

    bkg_proj = divideByBinWidth(bkg_proj,divideFlag=binWidthDivision)[0]

    if "CR_F" in outputFile:
        regionText="Untagged, fail"
    elif "CR_M" in outputFile:
        regionText="Untagged, medium"
    elif "CR_T" in outputFile:
        regionText="Untagged, tight"
    elif "F." in outputFile:
        regionText="$\gamma$-tagged, fail"
    elif "M." in outputFile:
        regionText="$\gamma$-tagged,\nmedium"
    elif "T." in outputFile:
        regionText="$\gamma$-tagged, tight"
    else:
        regionText=""

    plt.style.use([hep.style.CMS])
    matplotlib.rcParams.update({'font.size': 30})
    f, axs = plt.subplots(2,1, sharex=True, sharey=False,gridspec_kw={'height_ratios': [4, 1],'hspace': 0.00}, figsize=(10, 11),constrained_layout=True)
    axs = axs.flatten()
    plt.sca(axs[0])

    if any("gamma H" in label for label in labelsMC):
        if twoSignal:
            hep.histplot(histosMC[:-2],edges[0],stack=True,label = labelsMC[:-2], histtype="fill",facecolor=colorsMC[:-2])
            hep.histplot(histosMC[-1],edges[0],stack=False,label = "_nolegend_", histtype="step",color=colorsMC[-1])            
            hep.histplot(histosMC[-2],edges[0],stack=False,label = "_nolegend_", histtype="step",color=colorsMC[-2],linestyle="--")
            custom_handle1 = Patch(facecolor="none", edgecolor=colorsMC[-1], label=labelsMC[-1])  # Solid rectangle
            custom_handle2 = Patch(facecolor="none", edgecolor=colorsMC[-2], label=labelsMC[-2], linestyle='--')  # Dashed rectangle

        else:
            hep.histplot(histosMC[:-2],edges[0],stack=True,label = labelsMC[:-1], histtype="fill",facecolor=colorsMC[:-1])
            hep.histplot(histosMC[-1],edges[0],stack=True,label = "_nolegend_", histtype="step")
            custom_handle1 = Patch(facecolor="none", edgecolor=colorsMC[-1], label=labelsMC[-1])  # Solid rectangle
            
    else:
        hep.histplot(histosMC[:-1],edges[0],stack=True,label = labelsMC, histtype="fill",facecolor=colorsMC)
    plt.errorbar(centresData,hData, yerr=dataErrors, fmt='o',color="k",label = "Observed",zorder=11,markersize=10)

    # for i in range(len(uncBand[0])):
    #     binWidth            = edges[0][i+1]-edges[0][i]
    #     if(binWidthDivision):
    #         uncBand[0][i]       = uncBand[0][i]/binWidth
    #         uncBand[1][i]       = uncBand[1][i]/binWidth
    #     else:
    #         uncBand[0][i]       = uncBand[0][i]
    #         uncBand[1][i]       = uncBand[1][i]

    # uncBandLow = np.append(uncBand[0],[0],axis=0)
    # uncBandHi  = np.append(uncBand[1],[0],axis=0)#Hack to get the last bin uncertainty to plot, since we're using step="post"

    #plt.fill_between(edges[0],uncBandLow,uncBandHi,facecolor="none", hatch="xxx", edgecolor="grey", linewidth=0.0,step="post",zorder=10)

    def calcRatio(hData,hMC,dataErrs):
        ratioVals=[]
        ratioErrs=dataErrs
        for i in range(len(hData)):
            data      = hData[i]
            mc        = hMC[i]

            ratioVal     = data/(mc+0.000001)#Protect division by zero
            ratioErrs[0][i] = ratioErrs[0][i]/(mc+0.000001)
            ratioErrs[1][i] = ratioErrs[1][i]/(mc+0.000001)
            ratioVals.append(ratioVal)

        return ratioVals, ratioErrs

    def calcSystBand(hMC,uncBand):
        systBand = [[],[]]
        for i in range(len(hMC)):
            systBandUp = uncBand[0][i]/hMC[i]
            systBandDn = uncBand[1][i]/hMC[i]
            systBand[0].append(systBandUp)
            systBand[1].append(systBandDn)
        return systBand

    ratioVals, ratioErrs = calcRatio(hData,bkg_proj,dataErrors)
    systBand = calcSystBand(bkg_proj,uncBand)

    axs[0].legend()
    if(binWidthDivision):
        #plt.ylabel("Events/GeV",horizontalalignment='right', y=1.0)
        plt.ylabel("Events/GeV",horizontalalignment='center', y=0.5)
    else:
        #plt.ylabel("Events / bin",horizontalalignment='right', y=1.0)
        plt.ylabel("Events / 10 GeV",horizontalalignment='center', y=0.5)
    axs[1].set_ylabel("Data/bkg.")

    if(xRange):
        axs[0].set_xlim(xRange)
    if(yRange):
        axs[0].set_ylim(yRange)
    else:
        if "CR_M" in outputFile:
            yMaximum = max(hData)*1.7
        else:
            yMaximum = max(hData)*1.6
        axs[0].set_ylim([0.,yMaximum])

    if("16APV" in outputFile):
        lumi = 19.5
    elif("16" in outputFile):
        lumi = 16.8
    elif("17" in outputFile):
        lumi = 41.5
    elif("18" in outputFile):
        lumi = 59.8
    else:
        lumi = 138

    lumiText = str(lumi)+ " $fb^{-1} (13 TeV)$"    
    hep.cms.lumitext(lumiText)
    #hep.cms.text("Preliminary",loc=0)
    hep.cms.text("",loc=0)
    if any("gamma H" in label for label in labelsMC):
        handles, labels = plt.gca().get_legend_handles_labels()
        if twoSignal:
            handles.extend([custom_handle1, custom_handle2])
            labels.extend([labelsMC[-1], labelsMC[-2]])
        else:
            handles.extend([custom_handle1])
            labels.extend([labelsMC[-1]])
        plt.legend(handles, labels,loc=1,ncol=1)
    else:
        plt.legend(loc=1,ncol=1)

    if(projectionText):
        plt.text(0.60, 0.60, projectionText, horizontalalignment='center',verticalalignment='center',transform=axs[0].transAxes)
    

    plt.text(0.05, 0.95, regionText, transform=plt.gca().transAxes, verticalalignment='top')

    plt.sca(axs[1])#switch to lower pad
    axs[1].axhline(y=1.0, xmin=0, xmax=1, color="grey",linestyle="--",alpha=0.5)
    if "CR" in outputFile.split("/")[-1]:
        axs[1].set_ylim([0.95,1.05])
    else:
        axs[1].set_ylim([0.0,2.3])
    plt.xlabel(xlabel,horizontalalignment='center', x=0.5)
    plt.errorbar(centresData,ratioVals,yerr=ratioErrs, fmt='o',color="k",markersize=10) 

    #Hacky way to plot shaded region, assumes bin width = 10
    left_edges = np.array(centresData) - 5
    right_edges = np.array(centresData) + 5
    for i in range(len(centresData)):
        if i==0:
            label = "Bkg. uncertainty"
        else:
            label="_nolegend_"
        plt.fill_between([left_edges[i], right_edges[i]], systBand[0][i], systBand[1][i], color='lightgrey', alpha=0.8,label=label)
    plt.legend(bbox_to_anchor=(-0.01, 1.1),loc='upper left')
    axs[1].tick_params(axis='x', pad=10)  #Fix overlap between numbers on x and y axis of ratio plot
    #plt.fill_between(centresData, systBand[0], systBand[1], color='lightgrey', alpha=0.5)


    print("Saving ", outputFile)
    #plt.tight_layout()
    plt.savefig(outputFile)#,bbox_inches="tight")
    plt.savefig(outputFile.replace("png","pdf"))#,bbox_inches="tight")

    plt.clf()
    plt.cla()

def plotVJets(data,var,outFile,xTitle="",yTitle="",yRange=[],xRange=[],log=True,rebinX=1,luminosity="36.3",proj=""):
    histos = []
    labels  = []
    edges   = []
    colors  = []
    histosData = []#we're assuming only one data_obs dataset
    edgesData  = []#it's still kept in array (with one element) to be similar to other processes
    labelsData = []
    data = sorted(data.items(), key=lambda x: x[1]["order"])#VERY HANDY, reads samples in order
    for sample, sample_cfg in data:
        if not(("ZJets" in sample) or ("WJets" in sample)):
            continue

        tempFile = r.TFile.Open(sample_cfg["file"])
        if(proj=="X"):
            h = tempFile.Get("{0}_{1}".format(sample,var)).ProjectionX()
        elif(proj=="Y"):
            h = tempFile.Get("{0}_{1}".format(sample,var)).ProjectionY()
        else:
            h = tempFile.Get("{0}_{1}".format(sample,var))
        h.RebinX(rebinX)
        hist, edges = hist2array(h,return_edges=True)
        histos.append(hist)
        edges.append(edges[0])
        labels.append(sample_cfg["label"])
        colors.append(sample_cfg["color"])

    plt.style.use([hep.style.CMS])
    f, ax = plt.subplots()

    hep.histplot(histos,edges[0],stack=False,ax=ax,label=labels,linewidth=3,histtype="step",color=colors)
    if(log):
        ax.set_yscale("log")
    ax.legend()
    ax.set_ylabel(yTitle)
    ax.set_xlabel(xTitle)
    plt.xlabel(xTitle, horizontalalignment='right', x=1.0)
    plt.ylabel(yTitle,horizontalalignment='right', y=1.0)
    if(yRange):
        ax.set_ylim(yRange)
    if(xRange):
        ax.set_xlim(xRange)
    lumiText = luminosity + " $fb^{-1}\ (13 TeV)$"
    hep.cms.lumitext(text=lumiText, ax=ax, fontname=None, fontsize=None)
    hep.cms.text("Simulation WiP",loc=0)
    plt.legend(loc="best",ncol=2)#loc = 'best'
    plt.tight_layout()

    print("Saving {0}".format(outFile))

    plt.savefig(outFile)
    plt.savefig(outFile.replace(".png",".pdf"))

    plt.clf()


def plotSigTemplates(iFile,hName,outFile,xTitle="",yTitle="",yRange=[],xRange=[],log=True,rebinX=1,luminosity="36.3",projection="",slices=[0,-1],normalize=False):
    histos  = []
    labels  = []
    couplings = ["gha4","ghza4","ghza4_gha4","gha2","ghza2","ghza2_gha2"]
    tempFile = r.TFile.Open(iFile)
    linestyles = ["solid","solid","dotted","dotted","dashed","dashed"]

    for coupling in couplings:
        h = tempFile.Get("{0}__weight_{1}__nominal".format(hName,coupling))
        if normalize:
            h.Scale(1./h.Integral())

        if(projection=="x"):
            h = h.ProjectionX(h.GetName()+"_x",slices[0],slices[1])
        elif(projection=="y"):
            h = h.ProjectionY(h.GetName()+"_y",slices[0],slices[1])

        h.RebinX(rebinX)
        hist, edges = hist2array(h,return_edges=True)
        histos.append(hist)
        labels.append(coupling)

    plt.style.use([hep.style.CMS])
    f, ax = plt.subplots()
    plt.sca(ax)

    #--------------------------#
    hep.histplot(histos,edges[0],stack=False,ax=ax,label=labels,linestyle=linestyles,linewidth=3,histtype="step")
    if(log):
        ax.set_yscale("log")
    ax.legend()
    ax.set_ylabel(yTitle)
    ax.set_xlabel(xTitle)
    plt.ylabel(yTitle,horizontalalignment='right', y=1.0)
    
    if(yRange):
        ax.set_ylim(yRange)
    else:
        ax.set_ylim([None,max(histos[1])*1.5])
    if(xRange):
        ax.set_xlim(xRange)
    if(luminosity):
        lumiText = luminosity + " $fb^{-1}\ (13 TeV)$"
        hep.cms.lumitext(text=lumiText, ax=ax, fontname=None, fontsize=None)
    hep.cms.text("WiP",loc=0)
    plt.legend(loc="best",ncol=2,handletextpad=0.3)#loc = 'best'
    
    plt.tight_layout()
    print("Saving {0}".format(outFile))
    plt.savefig(outFile)
    plt.savefig(outFile.replace(".png",".pdf"))
    plt.cla()
    plt.clf()
    tempFile.Close()

def plotPostfitCombined(postfitShapesFile,odir,prefitTag=False,blind=True,binWidthDivision=True,debug=False,signal="HgammaHZy",plotSlices=False):

    labels              = ["Data","Non-resonant","W+Gamma","Z+Gamma","H+Gamma (HZy coupling)"]
    tags                = ["data_obs","qcd","WGamma","ZGamma",signal]
    colors              = ["black","gold","mistyrose","blue","red"]

    twoDShapes  = {}
    totalProcs  = {}

    region      = "M"
    for region in ["M","T"]:
        dirRegion   = region
        polyOrder = odir.split("/")[-2].split("SR")[0][-1]
        twoDShapes[region]=[]
        #Merge sliced histograms
        for tag in tags:
            if(tag=="qcd"):
                tag     = tag+"_"+polyOrder
            twoDShape   = get2DPostfitPlot(postfitShapesFile,tag,dirRegion,prefit=prefitTag)
            twoDShapes[region].append(twoDShape)
            totalProcs[region]  = get2DPostfitPlot(postfitShapesFile,"TotalProcs".format(region),dirRegion,prefit=prefitTag)

    #print(twoDShapes)
    twoDShapesMerged = []
    for i in range(len(tags)):
        shapeM = twoDShapes["M"][i]
        shapeT = twoDShapes["T"][i]
        shapeM.Add(shapeT)
        twoDShapesMerged.append(shapeM)

    twoDShapes = twoDShapesMerged

    totalProcsTemp = totalProcs["M"]
    totalProcsTemp.Add(totalProcs["T"])
    totalProcs     = totalProcsTemp

    xRange = [60,200]

    projections         = []
    for j,twoDShape in enumerate(twoDShapes):
        proj            = twoDShape.ProjectionX(tags[j]+"_projx")
        if(tags[j]=="data_obs" and blind):
            proj = blindHiggsMass(proj)
        #if(tags[j]=="SMHiggs"):
        #    proj.Scale(10.)
        projections.append(proj)
    totalProcs_proj     = totalProcs.ProjectionX("totalprocs_projx")
    projections.append(totalProcs_proj)
    uncBand_proj        = getUncBand(totalProcs_proj)

    plotName = "{0}/Combined_SR.png".format(odir)
    plotShapes(projections[0],projections[1:],uncBand_proj,labels[1:],colors[1:],"$M_{PNet}$ (GeV)",plotName,xRange=xRange,binWidthDivision=binWidthDivision)

    # f = r.TFile.Open("{0}/{1}_{2}.root".format(odir,outFile,region),"RECREATE")
    # f.cd()
    # totalProcs.Write()
    # for h in twoDShapes:
    #     h.Write()
    # f.Close()


if __name__ == '__main__':

    wp = "tight_medium"
    # #for year in ["2016","2016APV","2017","2018","RunII"]:
    # for year in ["RunII"]:
    #     odir = "results/plots/{0}/{1}/".format(wp,year)
    #     Path(odir).mkdir(parents=True, exist_ok=True)
    #     # if(year=="RunII"):
    #     #     plotRPFMC("results/templates/tight_medium/RunII/scaled/GJets.root",odir,passTag="M",failTag="F",xRange=[60,200],yTitle="$R_{M/F}$",rebinX=4)
    #     #     plotRPFMC("results/templates/tight_medium/RunII/scaled/GJets.root",odir,passTag="T",failTag="F",xRange=[60,200],yTitle="$R_{T/F}$",rebinX=4)
    #     #     plotRPFMC("results/templates/tight_medium/RunII/scaled/GJets.root",odir,passTag="T",failTag="M",xRange=[60,200],yTitle="$R_{T/M}$",rebinX=4)

        
    #     if(year=="2016APV"):
    #         luminosity="19.5"
    #     elif(year=="2016"):
    #         luminosity="16.8"
    #     elif(year=="2017"):
    #         luminosity="41.5"
    #     elif(year=="2018"):
    #         luminosity="59.8"
    #     elif(year=="RunII"):
    #         luminosity="138"

    #     with open("plotConfigs/{0}_{1}.json".format(year,wp)) as json_file:
    #         data = json.load(json_file)
    #         plotVarStack(data,"H_m_T__nominal","{0}/m_lin_T_data.png".format(odir),xTitle="$M_{PNet}$ (GeV)",yTitle="Events / 10 GeV",yRange=[0,None],log=False,xRange=[60,200],rebinX=2,luminosity=luminosity,mergeMassBins=False)
    #         plotVarStack(data,"H_m_M__nominal","{0}/m_lin_M_data.png".format(odir),xTitle="$M_{PNet}$ (GeV)",yTitle="Events / 10 GeV",yRange=[0,None],log=False,xRange=[60,200],rebinX=2,luminosity=luminosity,mergeMassBins=False)
    #         plotVarStack(data,"H_m_F__nominal","{0}/m_lin_F_data.png".format(odir),xTitle="$M_{PNet}$ (GeV)",yTitle="Events / 10 GeV",yRange=[0,None],log=False,xRange=[60,200],rebinX=2,luminosity=luminosity,mergeMassBins=False,blind=False)
    #         plotVarStack(data,"H_m_F__nominal","{0}/m_lin_F_data.png".format(odir),xTitle="$M_{PNet}$ (GeV)",yTitle="Events / 10 GeV",yRange=[0,None],log=False,xRange=[60,200],rebinX=2,luminosity=luminosity,mergeMassBins=False,blind=False)
            
    #         #This plotting function is not yet adjusted for new signal template naming
    #         plotVarStackMC(data,"Gamma_pT_T_nom","{0}/gamma_pT_T_data.png".format(odir),xTitle="Photon $p_{T}$ (GeV)",yTitle="Events / 50 GeV",yRange=[1.,None],log=True,xRange=[300,1000],rebinX=1,luminosity=luminosity,mergeMassBins=False)
    #         plotVarStackMC(data,"Gamma_pT_M_nom","{0}/gamma_pT_M_data.png".format(odir),xTitle="Photon $p_{T}$ (GeV)",yTitle="Events / 50 GeV",yRange=[1.,None],log=True,xRange=[300,1000],rebinX=1,luminosity=luminosity,mergeMassBins=False)
    #         plotVarStackMC(data,"Gamma_pT_F_nom","{0}/gamma_pT_F_data.png".format(odir),xTitle="Photon $p_{T}$ (GeV)",yTitle="Events / 50 GeV",yRange=[1.,None],log=True,xRange=[300,1000],rebinX=1,luminosity=luminosity,mergeMassBins=False)

    # #for year in ["2016APV","2016","2017","2018","RunII"]:
    for year in ["RunII"]:

        odir = "results/plots/tight_medium_CR/{0}/".format(year)
        Path(odir).mkdir(parents=True, exist_ok=True)

        #plotRPFMC("results/templates_CR/tight_medium/RunII/scaled/QCD.root",odir,passTag="CR_M",failTag="CR_F",xRange=[60,200],yTitle="$R_{M/F}^{0\gamma}$",rebinX=4,CRFlag=True)
        #plotRPFMC("results/templates_CR/tight_medium/RunII/scaled/QCD.root",odir,passTag="CR_T",failTag="CR_F",xRange=[60,200],yTitle="$R_{T/F}^{0\gamma}$",rebinX=4,CRFlag=True)

        
        if(year=="2016APV"):
            luminosity="19.5"
        elif(year=="2016"):
            luminosity="16.8"
        elif(year=="2017"):
            luminosity="41.5"
        elif(year=="2018"):
            luminosity="59.8"
        elif(year=="RunII"):
            luminosity="138"

        # with open("plotConfigs/{0}_tight_CR.json".format(year)) as json_file:
        #     data = json.load(json_file)
        #     plotVarStack(data,"H_m_pT_CR_T__nominal","{0}/m_lin_CR_T_data.png".format(odir),xTitle="$M_{PNet}$ (GeV)",yTitle="Events / 10 GeV",yRange=[0,None],log=False,xRange=[40,200],rebinX=2,luminosity=luminosity,mergeMassBins=False,projection="x",kFactor=1.15,blind=False,ratioRange=[0.8,1.2])
        #     plotVarStack(data,"H_m_pT_CR_M__nominal","{0}/m_lin_CR_M_data.png".format(odir),xTitle="$M_{PNet}$ (GeV)",yTitle="Events / 10 GeV",yRange=[0,None],log=False,xRange=[40,200],rebinX=2,luminosity=luminosity,mergeMassBins=False,projection="x",kFactor=1.15,blind=False,ratioRange=[0.8,1.2])
        #     plotVarStack(data,"H_m_pT_CR_F__nominal","{0}/m_lin_CR_F_data.png".format(odir),xTitle="$M_{PNet}$ (GeV)",yTitle="Events / 10 GeV",yRange=[0,None],log=False,xRange=[40,200],rebinX=2,luminosity=luminosity,mergeMassBins=False,projection="x",kFactor=0.95,blind=False,ratioRange=[0.8,1.2])
        #     plotVarStack(data,"H_m_pT_CR_F__nominal","{0}/m_CR_F_data.png".format(odir),xTitle="$M_{PNet}$ (GeV)",yTitle="Events / 10 GeV",yRange=[10**2,10**10],log=True,xRange=[40,200],rebinX=2,luminosity=luminosity,mergeMassBins=False,projection="x",kFactor=0.95,blind=False,ratioRange=[0.8,1.2])


        # with open("plotConfigs/Vjets{0}.json".format(year)) as json_file:
        #     data = json.load(json_file)
        #     plotVJets(data,"VpT_CR_F_nom","{0}/vpT_CR_F_lin.png".format(odir),xTitle="$Gen V p_{T}$ (GeV)",yTitle="Events / 10 GeV",log=False,xRange=[0,1500],yRange=[1.,5*10**4],rebinX=1,luminosity=luminosity)
        #     plotVJets(data,"VpT_CR_T_nom","{0}/vpT_CR_T_lin.png".format(odir),xTitle="$Gen V p_{T}$ (GeV)",yTitle="Events / 10 GeV",log=False,xRange=[0,1500],yRange=[1.,600],rebinX=1,luminosity=luminosity)
        #     plotVJets(data,"VpT_CR_F_nom","{0}/vpT_CR_F.png".format(odir),xTitle="$Gen V p_{T}$ (GeV)",yTitle="Events / 10 GeV",log=True,xRange=[0,1500],yRange=[1.,10**6],rebinX=1,luminosity=luminosity)
        #     plotVJets(data,"VpT_CR_T_nom","{0}/vpT_CR_T.png".format(odir),xTitle="$Gen V p_{T}$ (GeV)",yTitle="Events / 10 GeV",log=True,xRange=[0,1500],yRange=[1.,10**6],rebinX=1,luminosity=luminosity)
        #     plotVJets(data,"H_m_pT_CR_F__nominal","{0}/mVJets_CR_F_lin.png".format(odir),xTitle="$M_{PNet}$ (GeV)",yTitle="Events / 5 GeV",log=False,xRange=[40,150],yRange=[0,10**5],rebinX=1,luminosity=luminosity,proj="X")
        #     plotVJets(data,"H_m_pT_CR_T__nominal","{0}/mVJets_CR_T_lin.png".format(odir),xTitle="$M_{PNet}$ (GeV)",yTitle="Events / 5 GeV",log=False,xRange=[40,150],yRange=[0,2500],rebinX=1,luminosity=luminosity,proj="X")


    # Postfit
    cmsswArea       = "../CMSSW_10_6_14/src/"
    polyOrders      = ["1","1"]
    workingAreas    = ["SR_CR_standalone_HZy"]


    for workingArea in workingAreas:
        if("HZy" in workingArea):
            signal      = "HgammaHZy"
            coupling    = "HZy coupling"
        else:
            signal      = "HgammaHyy"
            coupling    = "Hyy coupling"
            
        baseDir         = cmsswArea + workingArea + "/{0}SR_{1}CR_area/".format(polyOrders[0],polyOrders[1])
        odir = "results/plots/{0}/{1}SR_{2}CR_area/".format(workingArea,polyOrders[0],polyOrders[1])
        Path(odir).mkdir(parents=True, exist_ok=True)
        fitFile         = baseDir+"postfitshapes_b.root"
        # if polyOrders == ["1","1"]:
        #     iFile  = baseDir+"/fitDiagnosticsTest.root"
        #     rpfTag = "qcd_rpfT"
        #     plotRPFWithUnc.getAndPlotRPF(iFile,rpfTag,odir,CRflag=False)
        #     rpfTag = "qcd_rpfM"
        #     plotRPFWithUnc.getAndPlotRPF(iFile,rpfTag,odir,CRflag=False)
        #     rpfTag = "qcd_CR_rpfT"
        #     plotRPFWithUnc.getAndPlotRPF(iFile,rpfTag,odir,CRflag=True)
        #     rpfTag = "qcd_CR_rpfM"
        #     plotRPFWithUnc.getAndPlotRPF(iFile,rpfTag,odir,CRflag=True)

        # plotPostfitCombined(fitFile,odir,binWidthDivision=False,signal=signal)
        # plotRPF(fitFile,odir,"qcd_{0}".format(polyOrders[0]),yTitle="$R_{M/F}$")
        # plotRPF(fitFile,odir,"qcd_{0}".format(polyOrders[0]),passTag="T",yTitle="$R_{T/F}$")
        # plotRPF(fitFile,odir,"qcd_CR_{0}".format(polyOrders[1]),passTag="CR_T",failTag="CR_F",xRange=[60,150],yTitle="$R_{T/F}^{0\gamma}$")
        # plotRPF(fitFile,odir,"qcd_CR_{0}".format(polyOrders[1]),passTag="CR_M",failTag="CR_F",xRange=[60,150],yTitle="$R_{M/F}^{0\gamma}$")
        # plotDeltaNLL(baseDir,odir,"DeltaNLL",extraText=coupling,xRange=[-40.,40.0])
        # plotDeltaNLL(baseDir,odir,"DeltaNLL",extraText="",xRange=[-2.,2.0])
        #plotDeltaNLLComp(baseDir,"./","DeltaNLL",extraText="",xRange=[0.,40.0])
        try:
            #For paper
            #plotPostfitTwoSignalsForPaper(fitFile,"T",odir)
            #plotPostfitTwoSignalsForPaper(fitFile,"M",odir)
            plotPostfitTwoSignalsForPaperWithRatio(fitFile,"T",odir)
            plotPostfitTwoSignalsForPaperWithRatio(fitFile,"M",odir)
            #plotPostfit(fitFile,"T",odir,binWidthDivision=False,signal=signal,blind=False,plotSlices=False)
            #plotPostfit(fitFile,"M",odir,binWidthDivision=False,signal=signal,blind=False,plotSlices=False)
            plotPostfit(fitFile,"F",odir,blind=False,binWidthDivision=False,signal=signal,plotSlices=False)
            plotPostfit(fitFile,"CR_T",odir,binWidthDivision=False,signal=signal,blind=False,plotSlices=False)
            plotPostfit(fitFile,"CR_M",odir,binWidthDivision=False,signal=signal,blind=False,plotSlices=False)
            plotPostfit(fitFile,"CR_F",odir,binWidthDivision=False,signal=signal,blind=False,plotSlices=False)
            #For AN
            # plotPostfit(fitFile,"T",odir,binWidthDivision=False,signal=signal,plotSlices=True,paper=False)
            # plotPostfit(fitFile,"M",odir,binWidthDivision=False,signal=signal,plotSlices=True,paper=False)
            # plotPostfit(fitFile,"F",odir,blind=False,binWidthDivision=False,signal=signal,plotSlices=True,paper=False)
            # plotPostfit(fitFile,"CR_T",odir,binWidthDivision=False,signal=signal,blind=False,paper=False)
            # plotPostfit(fitFile,"CR_M",odir,binWidthDivision=False,signal=signal,blind=False,paper=False)
            # plotPostfit(fitFile,"CR_F",odir,binWidthDivision=False,signal=signal,blind=False,paper=False)
        except:
            print("Couldn't plot for {0} {1}".format(workingArea,polyOrders))

    # #Plot unblinded RPF
    # baseDir         = "/users/mrogul/Work/Hgamma/test_datacards/CMSSW_11_3_4/src/datacards/"
    # odir = "results/plots/combination/"
    # Path(odir).mkdir(parents=True, exist_ok=True)
    # iFile  = baseDir+"/fitDiagnosticsTest.root"
    # rpfTag = "qcd_rpfT"
    # plotRPFWithUnc.getAndPlotRPF(iFile,rpfTag,odir,CRflag=False)
    # rpfTag = "qcd_rpfM"
    # plotRPFWithUnc.getAndPlotRPF(iFile,rpfTag,odir,CRflag=False)
    # rpfTag = "qcd_CR_rpfT"
    # plotRPFWithUnc.getAndPlotRPF(iFile,rpfTag,odir,CRflag=True)
    # rpfTag = "qcd_CR_rpfM"
    # plotRPFWithUnc.getAndPlotRPF(iFile,rpfTag,odir,CRflag=True)

# # plotSigTemplates("test.root","Hgamma_H_m_pT_T","MPNet_proj1_norm.png",xTitle="MPNet (GeV)",yTitle="Events / bin [a.u.]",yRange=[],log=False,rebinX=2,luminosity="",projection="x",slices=[1,1],normalize=True,xRange=[50,200])
# # plotSigTemplates("test.root","Hgamma_H_m_pT_T","MPNet_proj2_norm.png",xTitle="MPNet (GeV)",yTitle="Events / bin [a.u.]",yRange=[],log=False,rebinX=2,luminosity="",projection="x",slices=[2,-1],normalize=True,xRange=[50,200])
# # plotSigTemplates("test.root","Hgamma_H_m_pT_T","MPNet_proj1.png",xTitle="MPNet (GeV)",yTitle="Events / bin",yRange=[],log=False,rebinX=2,luminosity="",projection="x",slices=[1,1],xRange=[50,200])
# # plotSigTemplates("test.root","Hgamma_H_m_pT_T","MPNet_proj2.png",xTitle="MPNet (GeV)",yTitle="Events / bin",yRange=[],log=False,rebinX=2,luminosity="",projection="x",slices=[2,-1],xRange=[50,200])
