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

def compareNonResShapes():
    tplDir      = "/users/mrogul/Work/Hgamma/Hgamma/results/templates/tight_medium/RunII/scaled/"
    fileGJets   = "GJets.root"
    fileGJets   = r.TFile.Open(tplDir+fileGJets)
    hGJets      = fileGJets.Get("GJets_H_m_pT_F__nominal")
    hGJets      = hGJets.ProjectionX()
    hGJets.RebinX(2)
    hGJets.Scale(1./hGJets.Integral())
    hNumpyGJets, edges = hist2array(hGJets,return_edges=True)

    fileQCD   = "QCD.root"
    fileQCD   = r.TFile.Open(tplDir+fileQCD)
    hQCD      = fileQCD.Get("QCD_H_m_pT_F__nominal")
    hQCD      = hQCD.ProjectionX()
    hQCD.RebinX(2)
    hQCD.Scale(1./hQCD.Integral())
    hNumpyQCD, edges = hist2array(hQCD,return_edges=True)

    hRatio = hGJets.Clone("hRatio")
    hRatio.Divide(hQCD)

    nBinsX      = hGJets.GetNbinsX()
    hGJetsErrs  = []
    hQCDErrs    = []
    hRatioErrs     = []
    for i in range(1,nBinsX+1):
        err1   = hGJets.GetBinError(i)
        err2   = hQCD.GetBinError(i)

        hGJetsErrs.append(err1)
        hQCDErrs.append(err2)
        hRatioErrs.append(hRatio.GetBinError(i))

    hRatio, edges = hist2array(hRatio,return_edges=True)

    fileQCD.Close()
    fileGJets.Close()


    plt.style.use([hep.style.CMS])
    f, axs = plt.subplots(2,1, sharex=True, sharey=False,gridspec_kw={'height_ratios': [4, 1],'hspace': 0.05})
    axs = axs.flatten()
    plt.sca(axs[0])
    hep.histplot([hNumpyGJets,hNumpyQCD],yerr=[hGJetsErrs,hQCDErrs], bins=edges[0],label=["$\gamma$+jets", "QCD"])


    hep.cms.text("Work in progress",loc=0)
    axs[0].set_xlim([60.,200.])
    #ax.set_ylim([0.,maxRpf*1.3])

    plt.legend()

    lumiText = "138 $fb^{-1}\ (13 TeV)$"
    hep.cms.lumitext(text=lumiText)

    axs[1].set_xlabel("$M_{PNet}$ [GeV]",horizontalalignment='right', x=1.0)
    axs[0].set_ylabel("Event fraction / 10 GeV",horizontalalignment='right', y=1.0)


    plt.sca(axs[1])#switch to lower pad
    axs[1].set_ylim([0.5,1.5])
    hep.histplot(hRatio,edges[0],yerr=hRatioErrs,ax=axs[1],linewidth=1)#,histtype="step",edgecolor='red')
    axs[1].axhline(y=1.0, color="grey",linestyle="--")


    foutName = "nonResShapeComp"
    print("Saving {0}.pdf".format(foutName))
    plt.savefig("{0}.pdf".format(foutName), bbox_inches='tight')
    plt.savefig("{0}.png".format(foutName), bbox_inches='tight')
    plt.cla()
    plt.clf()



compareNonResShapes()
