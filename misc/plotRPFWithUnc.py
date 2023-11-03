import ROOT as r
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mplhep as hep


plt.style.use([hep.style.CMS])

def getA(M,pT,CRflag=False):
    #Assumes RPF = 0.002*(a+b*M+c*pT)
    #Used to calculate unc.**2 = A*covMat*A^T
    if CRflag:
        A = np.matrix([[0.002,0.002*M]])
    else:
        A = np.matrix([[0.002,0.002*M,0.002*pT]])
    return A

def getUnc(M,pT,covMat,CRflag=False):
    if CRflag:
        A    = getA(M,pT,CRflag)
        AT   = A.transpose()
    else:
        A    = getA(M,pT,CRflag)
        AT   = A.transpose()
    return np.sqrt(np.linalg.multi_dot([A,covMat,AT]))

def getVal(M,pT,rpfParams,CRflag=False):
    if CRflag:
        val = 0.002*(rpfParams[0]+rpfParams[1]*M)
    else:
        val = 0.002*(rpfParams[0]+rpfParams[1]*M+rpfParams[2]*pT)
    return val

def plot(rpfTag,rpfParams,rpfCovMat,foutName="test",odir=".",CRflag=False):
    if CRflag:
        masses = np.linspace(65,195,num=14)
        massMax = 200.
    else:
        masses = np.linspace(65,145,num=9)
        massMax = 150.
    massMin = 60.
    pTMin   = 300.
    pTMax   = 2000.
    f, ax   = plt.subplots()
    hep.cms.text("Work in progress",loc=0)
    pTLabels={350.:"300<$p_T$<400 GeV",1200.:"$p_T$>400 GeV"}
    rpfMax  = 0.
    xRange  = [60,massMax]

    if rpfTag == "qcd_rpfT":
        ylabel  = "$R_{T/F}$ x $10^{3}$"
    elif rpfTag =="qcd_rpfM":
        ylabel  = "$R_{M/F}$ x $10^{3}$"
    elif rpfTag == "qcd_CR_rpfT":
        ylabel  = "$R_{T/F}^{0\gamma}$ x $10^{3}$"
    elif rpfTag == "qcd_CR_rpfM":
        ylabel  = "$R_{M/F}^{0\gamma}$ x $10^{3}$"
    else:
        raise ValueError('Unknown rpfTag: {0}'.format(rpfTag))

    if(CRflag):
        plt_x       = []
        plt_y       = []
        plt_y_err   = []
        for mass in masses:
            massRescaled = (mass-massMin)/(massMax-massMin)    
            rpfVal = getVal(massRescaled,0,rpfParams,CRflag)
            rpfUnc = getUnc(massRescaled,0,rpfCovMat,CRflag).item()
            plt_x.append(mass)
            plt_y.append(rpfVal*1000)
            plt_y_err.append(rpfUnc*1000)
            if rpfVal*1000>rpfMax:
                rpfMax = rpfVal*1000
        plt.errorbar(plt_x, plt_y,yerr=plt_y_err,xerr=5,ls='none')
    else:
        for pT in [350., 1200.]:
            pTRescaled  = (pT-pTMin)/(pTMax-pTMin)#For RPF evaluation, we scale bin coordinates to [0,1] range
            plt_x       = []
            plt_y       = []
            plt_y_err   = []
            for mass in masses:
                massRescaled = (mass-massMin)/(massMax-massMin)    
                rpfVal = getVal(massRescaled,pTRescaled,rpfParams,CRflag)
                rpfUnc = getUnc(massRescaled,pTRescaled,rpfCovMat,CRflag).item()
                plt_x.append(mass)
                plt_y.append(rpfVal*1000)
                plt_y_err.append(rpfUnc*1000)
                if rpfVal*1000>rpfMax:
                    rpfMax = rpfVal*1000

    #print(plt_x,plt_y,plt_y_err)
            plt.errorbar(plt_x, plt_y,yerr=plt_y_err,xerr=5,ls='none',label=pTLabels[pT])
        plt.legend()

    ax.set_xlim(xRange)
    ax.set_ylim([0.,rpfMax*1.3])

    plt.xlabel("$M_{PNet}$ [GeV]",horizontalalignment='right', x=1.0)
    plt.ylabel(ylabel,horizontalalignment='right', y=1.0)
    #ax.yaxis.set_tick_params(which='minor', left=False)    
    #ax.yaxis.set_tick_params(which='minor', right=False)    

    print("Saving "+odir+"/{0}.pdf".format(foutName))
    plt.savefig(odir+"/{0}.pdf".format(foutName), bbox_inches='tight')
    plt.savefig(odir+"/{0}.png".format(foutName), bbox_inches='tight')
    plt.cla()
    plt.clf()

def getAndPlotRPF(iFile,rpfTag,odir="/.",CRflag=False):
    #Assumes RPF = 0.002*(a+b*M+c*pT)
    #Assumes CR RPF = 0.002*(a+b*M)
    f = r.TFile.Open(iFile)

    fit_b       = f.Get("fit_b")
    covMat      = fit_b.covarianceMatrix()
    finalPars   = fit_b.floatParsFinal()
    nParams     = covMat.GetNcols()

    rpfIndices  = []
    rpfParams   = []
    for cm_index in range(nParams):
        parName = finalPars.at(cm_index).GetName()
        if (rpfTag in parName):
            rpfIndices.append(cm_index)
            rpfParams.append(finalPars.at(cm_index).getValV())
            #print(parName, cm_index, covMat[cm_index][cm_index])
        else:
            continue

    nRPFParams = len(rpfIndices)
    rpfCovMat  = np.zeros((nRPFParams,nRPFParams))

    for i,rpfIndexI in enumerate(rpfIndices):
        for j,rpfIndexJ in enumerate(rpfIndices):
            rpfCovMat[i][j] = covMat[rpfIndexI][rpfIndexJ]

    #print(rpfCovMat)
    #print("RPF parameters: ", rpfParams)
    plot(rpfTag,rpfParams,rpfCovMat,foutName=rpfTag,odir=odir,CRflag=CRflag)
    f.Close()


if __name__ == '__main__':

    iFile  = "/users/mrogul/Work/Hgamma/CMSSW_10_6_14/src/SR_CR_HZy_ANv5/1SR_1CR_area/fitDiagnosticsTest.root"
    odir   = "."
    rpfTag = "qcd_rpfT"
    getAndPlotRPF(iFile,rpfTag,odir,CRflag=False)
    rpfTag = "qcd_rpfM"
    getAndPlotRPF(iFile,rpfTag,odir,CRflag=False)
    rpfTag = "qcd_CR_rpfT"
    getAndPlotRPF(iFile,rpfTag,odir,CRflag=True)
    rpfTag = "qcd_CR_rpfM"
    getAndPlotRPF(iFile,rpfTag,odir,CRflag=True)