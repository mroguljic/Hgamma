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


dir1 = "/users/mrogul/Work/Hgamma/Hgamma/results/plots/SR_CR_HZy_powerLaw/1SR_1CR_area"
dir2 = "/users/mrogul/Work/Hgamma/Hgamma/results/plots/SR_CR_HZy/1SR_1CR_area"

for region in ["T","M"]:
    f1   = r.TFile.Open("{0}/postfit_{1}.root".format(dir1,region))
    f2   = r.TFile.Open("{0}/postfit_{1}.root".format(dir2,region))

    h1   = f1.Get("h2_qcd_1_{0}".format(region))
    h1.Add(f1.Get("h2_WGamma_{0}".format(region)))
    h1.Add(f1.Get("h2_ZGamma_{0}".format(region)))
    h2   = f2.Get("h2_qcd_1_{0}".format(region))
    h2.Add(f2.Get("h2_WGamma_{0}".format(region)))
    h2.Add(f2.Get("h2_ZGamma_{0}".format(region)))

    h1   = h1.ProjectionX("powerLaw")
    h2   = h2.ProjectionX("poly")

    hRatio = h1.Clone("hRatio")
    hRatio.Divide(h2)

    # nBinsX = h1.GetNbinsX()
    # nBinsY = h1.GetNbinsY()

    # for i in range(1,nBinsX+1):
    #   for j in range(1,nBinsY+1):
    #       yield1 = h1.GetBinContent(i,j)
    #       yield2 = h2.GetBinContent(i,j)
    #       err1   = h1.GetBinError(i,j)
    #       err2   = h2.GetBinError(i,j)

    #       print("Bin {0}-{1}: {2:.2f}+/-{3:.2f} | {4:.2f}+/-{5:.2f}".format(i,j,yield1,err1,yield2,err2))
     
    nBinsX = h1.GetNbinsX()
    h1Errs = []
    h2Errs = []
    hRatioErrs = []
    for i in range(1,nBinsX+1):
        yield1 = h1.GetBinContent(i)
        yield2 = h2.GetBinContent(i)
        err1   = h1.GetBinError(i)
        err2   = h2.GetBinError(i)

        h1Errs.append(err1)
        h2Errs.append(err2)
        hRatioErrs.append(hRatio.GetBinError(i))

        print("Bin {0}: {2:.2f}+/-{3:.2f} | {4:.2f}+/-{5:.2f}".format(i,i,yield1,err1,yield2,err2))


    h1, edges = hist2array(h1,return_edges=True)
    h2        = hist2array(h2,return_edges=False)
    hRatio    = hist2array(hRatio,return_edges=False)

    h1[0]=h1[0]/2.#First bin is twice as large as others
    h2[0]=h2[0]/2.
    h1Errs[0]=h1Errs[0]/2.
    h2Errs[0]=h2Errs[0]/2.

    h1 = h1/5.#Display Events / GeV
    h2 = h2/5.
    h1Errs = np.array(h1Errs)
    h2Errs = np.array(h2Errs)
    h1Errs = h1Errs/5.
    h2Errs = h2Errs/5.

    plt.style.use([hep.style.CMS])
    f, axs = plt.subplots(2,1, sharex=True, sharey=False,gridspec_kw={'height_ratios': [4, 1],'hspace': 0.05})
    axs = axs.flatten()
    plt.sca(axs[0])
    hep.histplot([h1,h2],yerr=[h1Errs,h2Errs], bins=edges[0],label=["Power law", "Polynomial"])


    hep.cms.text("Work in progress",loc=0)
    axs[0].set_xlim([60.,200.])
    #ax.set_ylim([0.,maxRpf*1.3])

    plt.legend()

    lumiText = "138 $fb^{-1}\ (13 TeV)$"
    hep.cms.lumitext(text=lumiText, ax=axs[0], fontname=None, fontsize=None)

    axs[1].set_xlabel("$M_{PNet}$ [GeV]",horizontalalignment='right', x=1.0)
    axs[0].set_ylabel("Events / GeV",horizontalalignment='right', y=1.0)
    #ax.yaxis.set_tick_params(which='minor', left=False)    
    #ax.yaxis.set_tick_params(which='minor', right=False)    

    plt.sca(axs[1])#switch to lower pad
    axs[1].set_ylim([0.5,1.5])
    hep.histplot(hRatio,edges[0],yerr=hRatioErrs,ax=axs[1],linewidth=1)#,histtype="step",edgecolor='red')
    axs[1].axhline(y=1.0, color="grey",linestyle="--")




    foutName = "RpfSyst_{0}".format(region)
    print("Saving {0}.pdf".format(foutName))
    plt.savefig("{0}.pdf".format(foutName), bbox_inches='tight')
    plt.savefig("{0}.png".format(foutName), bbox_inches='tight')
    plt.cla()
    plt.clf()

    f1.Close()
    f2.Close()