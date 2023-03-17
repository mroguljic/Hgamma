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

matplotlib.use('Agg')
r.gROOT.SetBatch(True)
r.gStyle.SetOptFit(111)

def plotRPF(postfitShapesFile,odir,polyOrder,passTag="M",failTag="F"):
    hPass = get2DPostfitPlot(postfitShapesFile,"qcd_{0}".format(polyOrder),passTag)
    hFail = get2DPostfitPlot(postfitShapesFile,"qcd",failTag)

    hPass = hPass.ProjectionX("hPass_temp")
    hFail = hFail.ProjectionX("hFail_temp")

    foutName = "R{0}{1}".format(passTag,failTag)

    hRpf  = hPass.Clone("hRPF")
    hRpf.Divide(hFail)
    hRpf.Scale(1000)
    maxRpf = hRpf.GetMaximum()

    hRpf, edges = hist2array(hRpf,return_edges=True)

    plt.style.use([hep.style.CMS])
    f, ax = plt.subplots()
    hep.histplot(hRpf, bins=edges[0])

    hep.cms.text("Work in progress",loc=0)
    
    ax.set_xlim([60,200])
    ax.set_ylim([0.,maxRpf*1.3])

    plt.xlabel("$M_{PNet}$ [GeV]",horizontalalignment='right', x=1.0)
    plt.ylabel("$R_{{{0}/{1}}}$".format(passTag,failTag)+" x $10^{3}$",horizontalalignment='right', y=1.0)
    #ax.yaxis.set_tick_params(which='minor', left=False)    
    #ax.yaxis.set_tick_params(which='minor', right=False)    

    print("Saving "+odir+"{0}.pdf".format(foutName))
    plt.savefig(odir+"/{0}.pdf".format(foutName), bbox_inches='tight')
    plt.savefig(odir+"/{0}.png".format(foutName), bbox_inches='tight')
    plt.cla()
    plt.clf()


def plotTriggerEff(hPass,hTotal,year,luminosity,outFile,xlabel="Photon $p_{T}$ [GeV]",ylabel="Trigger efficiency / 10 GeV"):
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

def blindHiggsMass(hist):
    for i in range(1,hist.GetNbinsX()+1):
        binCenter = hist.GetBinCenter(i)
        if binCenter>100 and binCenter<150:
            hist.SetBinContent(i,0)
    return hist

def plotVarStack(data,var,outFile,xTitle="",yTitle="",yRange=[],xRange=[],log=True,rebinX=1,luminosity="36.3",mergeMassBins=False,blind=True,kFactor=1.25):
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

    data = sorted(data.items(), key=lambda x: x[1]["order"])#VERY HANDY, reads samples in order
    for sample, sample_cfg in data:
        if("singlephoton" in sample.lower() or "data" in sample.lower()):
            dataFlag = True
        else:
            dataFlag = False
        if("hgamma" in sample.lower()):
            sigFlag  = True
        else:
            sigFlag  = False

        tempFile = r.TFile.Open(sample_cfg["file"])
        #print("Opening ", sample_cfg["file"])
        h = tempFile.Get("{0}_{1}".format(sample,var))
        if("GJets" in sample or "QCD" in sample):
            h.Scale(kFactor)
        h.RebinX(rebinX)
        if(dataFlag and blind):
            h = blindHiggsMass(h)
        hist, edges = hist2array(h,return_edges=True)
        if(mergeMassBins):
            hist, edges = mergeLoMassBins(hist,edges[0])
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

    plt.style.use([hep.style.CMS])
    f, axs = plt.subplots(2,1, sharex=True, sharey=False,gridspec_kw={'height_ratios': [4, 1],'hspace': 0.05})
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
    if(mergeMassBins):
        for i in range(len(edges[0])-1):
            xerror = (edges[0][i+1]-edges[0][i])/2.
            xerrorsData.append(xerror)
    #--------------------------#
    hep.histplot(histos,edges[0],stack=True,ax=axs[0],label=labels,linewidth=1,histtype="fill",facecolor=colors,edgecolor='black')
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
            yRange[1] = dataMax*1.5
        axs[0].set_ylim(yRange)
    else:
        axs[0].set_ylim([0,dataMax*1.5])
    if(xRange):
        axs[0].set_xlim(xRange)
    lumiText = luminosity + " $fb^{-1}\ (13 TeV)$"
    hep.cms.lumitext(text=lumiText, ax=axs[0], fontname=None, fontsize=None)
    hep.cms.text("WiP",loc=0)
    plt.legend(loc="best",ncol=2,handletextpad=0.3)#loc = 'best'
    


    totalMC     = [sum(x) for x in zip(*histos)]
    dataMCRatio = [a/b for a,b in zip(histosData[0],totalMC)]
    plt.sca(axs[1])#switch to lower pad
    axs[1].set_ylim([0,2])
    hep.histplot(dataMCRatio,edges[0],ax=axs[1],linewidth=1,histtype="step",edgecolor='red')
    axs[1].axhline(y=1.0, color="grey",linestyle="--")

    plt.tight_layout()
    print("Saving {0}".format(outFile))
    plt.savefig(outFile)
    plt.savefig(outFile.replace(".png",".pdf"))
    plt.cla()
    plt.clf()


def plotVarStackMC(data,var,outFile,xTitle="",yTitle="",yRange=[],xRange=[],log=True,rebinX=1,luminosity="36.3",mergeMassBins=False):
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
        #print("Opening ", sample_cfg["file"])
        h = tempFile.Get("{0}_{1}".format(sample,var))
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
    f, axs = plt.subplots(2,1, sharex=True, sharey=False,gridspec_kw={'height_ratios': [4, 1],'hspace': 0.05})
    axs = axs.flatten()
    plt.sca(axs[0])

    hep.histplot(histosMC[:-1],edges[0],stack=True,ax=axs[0],label = labelsMC, histtype="fill",facecolor=colorsMC,edgecolor='black')
    plt.errorbar(centresData,hData, yerr=dataErrors, xerr=xerrorsData, fmt='o',color="k",label = "Data")

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

    plt.fill_between(edges[0],uncBandLow,uncBandHi,facecolor="none", hatch="xxx", edgecolor="grey", linewidth=0.0,step="post")

    axs[0].legend()
    if(binWidthDivision):
        plt.ylabel("Events/GeV",horizontalalignment='right', y=1.0)
    else:
        plt.ylabel("Events",horizontalalignment='right', y=1.0)
    axs[1].set_ylabel("Pulls")

    if(xRange):
        axs[0].set_xlim(xRange)
    if(yRange):
        axs[0].set_ylim(yRange)
    else:
        yMaximum = max(hData)*1.5+10.
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
        plt.text(0.60, 0.15, projectionText, horizontalalignment='center',verticalalignment='center',transform=axs[0].transAxes)


    plt.sca(axs[1])#switch to lower pad
    #axs[1].axhline(y=0.0, xmin=0, xmax=1, color="r")
    axs[1].axhline(y=1.0, xmin=0, xmax=1, color="grey",linestyle="--",alpha=0.5)
    axs[1].axhline(y=-1.0, xmin=0, xmax=1, color="grey",linestyle="--",alpha=0.5)
    axs[1].axhline(y=2.0, xmin=0, xmax=1, color="grey",linestyle="-.",alpha=0.5)
    axs[1].axhline(y=-2.0, xmin=0, xmax=1, color="grey",linestyle="-.",alpha=0.5)
    axs[1].set_ylim([-2.5,2.5])
    plt.xlabel(xlabel,horizontalalignment='right', x=1.0)

    pulls = calculatePull(hData,dataErrors,histosMC[-1],uncBand)
    hep.histplot(pulls,edges[0],ax=axs[1],linewidth=1,histtype="fill",facecolor="grey",edgecolor='black')

    print("Saving ", outputFile)
    #plt.tight_layout()
    plt.savefig(outputFile,bbox_inches="tight")
    plt.savefig(outputFile.replace("png","pdf"))

    plt.clf()
    plt.cla()

def printMCYields(data,region,year):
    data = sorted(data.items(), key=lambda x: x[1]["order"])
    print("----------")
    print("Yields in {0} - {1}".format(year,region))
    for sample, sample_cfg in data:
        f       = r.TFile.Open(sample_cfg["file"])
        hist    = f.Get("{0}_H_m_{1}__nominal".format(sample,region))
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

def plotPostfit(postfitShapesFile,region,odir,prefitTag=False,blind=True,binWidthDivision=True):

    if(region=="pass" or region=="T" or region=="M"):
        labels              = ["Data","Non-resonant","Z+Gamma","H+Gamma"]
        tags                = ["data_obs","qcd","ZGamma","Hgamma"]
        colors              = ["black","deepskyblue","blue","red"]
    else:
        labels              = ["Data","Non-resonant","Z+Gamma"]
        tags                = ["data_obs","qcd","ZGamma"]
        colors              = ["black","deepskyblue","blue"]

    if(prefitTag):
        outFile = "prefit"
    else:
        outFile = "postfit"


    twoDShapes  = []

    dirRegion   = region
    polyOrder   = odir.split("/")[-2]

    #Merge sliced histograms
    for tag in tags:
        if(tag=="qcd" and (region=="pass" or region=="T" or region=="M")):
            tag     = tag+"_"+polyOrder
        twoDShape   = get2DPostfitPlot(postfitShapesFile,tag,dirRegion,prefit=prefitTag)
        twoDShapes.append(twoDShape)
        totalProcs  = get2DPostfitPlot(postfitShapesFile,"TotalProcs".format(region),dirRegion,prefit=prefitTag)
    
    projections         = []
    for j,twoDShape in enumerate(twoDShapes):
        proj            = twoDShape.ProjectionX(tags[j]+"_projx")
        if(tags[j]=="data_obs" and blind):
            proj = blindHiggsMass(proj)
        projections.append(proj)
    totalProcs_proj     = totalProcs.ProjectionX("totalprocs_projx")
    projections.append(totalProcs_proj)
    uncBand_proj        = getUncBand(totalProcs_proj)

    if binWidthDivision:
        plotName = "{0}/{1}_{2}_rescaled.png".format(odir,outFile,region)
    else:
        plotName = "{0}/{1}_{2}.png".format(odir,outFile,region)
    if binWidthDivision:
        plotShapes(projections[0],projections[1:],uncBand_proj,labels[1:],colors[1:],"$M_{PNet}$ [GeV]",plotName,xRange=[60,200],yRange=[0,20],binWidthDivision=binWidthDivision)
    else:
        plotShapes(projections[0],projections[1:],uncBand_proj,labels[1:],colors[1:],"$M_{PNet}$ [GeV]",plotName,xRange=[60,200],binWidthDivision=binWidthDivision)

if __name__ == '__main__':

    wp = "tight_medium"
    # for year in ["2017"]:
    #     odir = "results_AK15/plots/{0}/{1}/".format(wp,year)
    #     Path(odir).mkdir(parents=True, exist_ok=True)        
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
    #         # plotVarStack(data,"H_m_T__nominal","{0}/m_lin_T_data.png".format(odir),xTitle="$M_{PNet}$ [GeV]",yTitle="Events / GeV",yRange=[0,None],log=False,xRange=[60,200],rebinX=1,luminosity=luminosity,mergeMassBins=True)
    #         # plotVarStack(data,"H_m_M__nominal","{0}/m_lin_M_data.png".format(odir),xTitle="$M_{PNet}$ [GeV]",yTitle="Events / GeV",yRange=[0,None],log=False,xRange=[60,200],rebinX=1,luminosity=luminosity,mergeMassBins=True)
    #         # plotVarStack(data,"H_m_F__nominal","{0}/m_lin_F_data.png".format(odir),xTitle="$M_{PNet}$ [GeV]",yTitle="Events / GeV",yRange=[0,None],log=False,xRange=[60,200],rebinX=1,luminosity=luminosity,mergeMassBins=True,blind=False)
    #         # plotVarStack(data,"H_m_F__nominal","{0}/m_lin_F_data.png".format(odir),xTitle="$M_{PNet}$ [GeV]",yTitle="Events / GeV",yRange=[0,None],log=False,xRange=[60,200],rebinX=1,luminosity=luminosity,mergeMassBins=True,blind=False)
    #         # plotVarStackMC(data,"Gamma_pT_T_nom","{0}/gamma_pT_T_data.png".format(odir),xTitle="Photon $p_{T}$ [GeV]",yTitle="Events / 50 GeV",yRange=[1.,None],log=True,xRange=[200,1000],rebinX=1,luminosity=luminosity,mergeMassBins=False)
    #         # plotVarStackMC(data,"Gamma_pT_M_nom","{0}/gamma_pT_M_data.png".format(odir),xTitle="Photon $p_{T}$ [GeV]",yTitle="Events / 50 GeV",yRange=[1.,None],log=True,xRange=[200,1000],rebinX=1,luminosity=luminosity,mergeMassBins=False)
    #         # plotVarStackMC(data,"Gamma_pT_F_nom","{0}/gamma_pT_F_data.png".format(odir),xTitle="Photon $p_{T}$ [GeV]",yTitle="Events / 50 GeV",yRange=[1.,None],log=True,xRange=[200,1000],rebinX=1,luminosity=luminosity,mergeMassBins=False)
    #         f = r.TFile.Open(data["data_obs"]["file"])
    #         print(data["data_obs"]["file"])
    #         hTotal = f.Get("data_obs_GammapTnoTriggers_nom")
    #         hPass  = f.Get("data_obs_GammapTtriggersAll_nom")
    #         hPass.RebinX(5)
    #         hTotal.RebinX(5)
    #         eff = r.TEfficiency(hPass,hTotal)
    #         eff.SetName("trig_eff")
    #         #g   = r.TFile.Open("data/trig_eff_{0}.root".format(year),"RECREATE")
    #         g   = r.TFile.Open("trig_eff_{0}.root".format(year),"RECREATE")
    #         g.cd()
    #         eff.Write()
    #         g.Close()

    #         plotTriggerEff(hPass,hTotal,year,luminosity,"{0}/Trig_eff_{1}.pdf".format(odir,year),ylabel="Trigger efficiency / 50 GeV")

    #Postfit
    cmsswArea       = "../CMSSW_10_6_14/src/"
    polyOrder       = "2"
    workingAreas    = ["AK15_2017_tight_medium","AK8_2017_tight_medium"]
    for workingArea in workingAreas:
        baseDir         = cmsswArea + workingArea + "/" + polyOrder + "_area/"
        fitFile         = baseDir+"postfitshapes_b.root"
        Path("results_AK15/plots/{0}/{1}/".format(workingArea,polyOrder)).mkdir(parents=True, exist_ok=True)
        #plotRPF(fitFile,"results_AK15/plots/{0}/{1}/".format(workingArea,polyOrder),polyOrder)
        #plotRPF(fitFile,"results_AK15/plots/{0}/{1}/".format(workingArea,polyOrder),polyOrder,passTag="T")

        try:
            plotPostfit(fitFile,"T","results_AK15/plots/{0}/{1}/".format(workingArea,polyOrder),binWidthDivision=False)
            plotPostfit(fitFile,"M","results_AK15/plots/{0}/{1}/".format(workingArea,polyOrder),binWidthDivision=False)
            plotPostfit(fitFile,"F","results_AK15/plots/{0}/{1}/".format(workingArea,polyOrder),blind=False,binWidthDivision=False)
            plotPostfit(fitFile,"T","results_AK15/plots/{0}/{1}/".format(workingArea,polyOrder),binWidthDivision=True)
            plotPostfit(fitFile,"M","results_AK15/plots/{0}/{1}/".format(workingArea,polyOrder),binWidthDivision=True)
            plotPostfit(fitFile,"F","results_AK15/plots/{0}/{1}/".format(workingArea,polyOrder),blind=False,binWidthDivision=True)
        except:
           print("Couldn't plot for {0} {1}".format(workingArea,polyOrder))
           