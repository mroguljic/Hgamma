import ROOT as r
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mplhep as hep
from root_numpy import hist2array


def plotRPF(iFile,oFile):
    iFile       = r.TFile.Open(iFile)
    hMedium     = iFile.Get("GJets_H_m_M__nominal")
    hTight      = iFile.Get("GJets_H_m_T__nominal")
    hPass       = hMedium.Clone("hPass")
    hPass.Add(hTight)
    hFail       = iFile.Get("GJets_H_m_F__nominal")
    hRpf        = hPass.Clone("hRpf")
    hRpf.Divide(hFail)
    hRpf.Scale(100)

    maxFail     = hFail.GetMaximum()
    minRpf      = hRpf.GetMinimum()
    maxRpf      = hRpf.GetMaximum()


    hPass, edges = hist2array(hPass,return_edges=True)
    hFail, edges = hist2array(hFail,return_edges=True)
    hRpf, edges  = hist2array(hRpf,return_edges=True)


    plt.style.use([hep.style.CMS])
    f, axs = plt.subplots(2,1, sharex=True, sharey=False,gridspec_kw={'height_ratios': [1, 1],'hspace': 0.05})
    axs = axs.flatten()
    plt.sca(axs[0])
    hep.histplot([hPass,hFail],edges[0],stack=False,ax=axs[0],label=["Pass","Fail"],linewidth=1,histtype="step",color=["red","blue"])
    hep.cms.text("WiP",loc=0)
    plt.sca(axs[1])#switch to lower pad
    hep.histplot(hRpf,edges[0],ax=axs[1],linewidth=2,histtype="step",edgecolor='darkgreen')


    axs[0].set_ylim([1,maxFail*10])
    axs[1].set_ylim([0,maxRpf*1.3])
    axs[1].set_xlim([50,200])

    axs[0].set_yscale("log")
    axs[0].legend()
    axs[0].set_ylabel("Events / bin")
    axs[1].set_xlabel("$M_{PNet}$")
    axs[1].set_ylabel("$R_{P/F}$ x 100")

    plt.tight_layout()
    print("Saving {0}".format(oFile))
    plt.savefig(oFile)
    plt.savefig(oFile.replace(".png",".pdf"))
    plt.cla()
    plt.clf()


origFile = "../results/templates/tight_medium_orig/RunII/scaled/GJets.root"
testFile = "../results/templates/tight_medium/RunII/scaled/GJets.root"
plotRPF(origFile,"inclusiveFail.png")
plotRPF(testFile,"failAbove0p2.png")