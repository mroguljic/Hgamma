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
r.gStyle.SetOptStat(0000)

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
    plt.errorbar(binsX, effs, yerr=[errsDn,errsUp], fmt='.',color="black",label="{0} Data trigger efficiency".format(year))

    plt.xlabel(xlabel, horizontalalignment='right', x=1.0)
    plt.ylabel(ylabel,horizontalalignment='right', y=1.0)

    if(luminosity):
        lumiText = luminosity + " $fb^{-1}\ (13 TeV)$"
    else:
        lumiText = year + " (13 TeV)"
    hep.cms.lumitext(text=lumiText, ax=ax, fontname=None, fontsize=None)
    hep.cms.text("WiP",loc=0)
    ax.set_ylim([0.5,1.02])
    plt.legend(loc="best")#loc = 'best'
    plt.tight_layout()

    print("Saving {0}".format(outFile))
    plt.savefig(outFile)
    plt.savefig(outFile.replace("pdf","png"))
    plt.clf()

def trigEffFromFile(iFile,year,lumi,writeToRoot=False):
        f = r.TFile.Open(iFile)
        hTotal = f.Get("data_pT_denominator")
        hPass  = f.Get("data_pT_numerator")
        eff = r.TEfficiency(hPass,hTotal)
        eff.SetName("trig_eff")
        if(writeToRoot):
            g   = r.TFile.Open("trig_eff_{0}.root".format(year),"RECREATE")
            g.cd()
            eff.Write()
            g.Close()

        plotTriggerEff(hPass,hTotal,year,lumi,"plots/Trig_eff_{0}.pdf".format(year))

def rwt2DNum(h2,eff):
    hRes = h2.Clone(h2.GetName()+"_eff_rwt")
    for i in range(1,h2.GetNbinsX()+1):
        for j in range(1,h2.GetNbinsY()+1):
            eff_pt      = eff.GetEfficiency(j)
            yield_rwt   =  hRes.GetBinContent(i,j)/eff_pt
            hRes.SetBinContent(i,j,yield_rwt)
    return hRes

def testMassDependency(iFile):
    f           = r.TFile.Open(iFile)
    m_pt_den    = f.Get("data_m_pT_denominator")
    m_pt_num    = f.Get("data_m_pT_numerator")

    m_pt_den.RebinX(5)
    m_pt_num.RebinX(5)
    m_pt_den.RebinY(2)
    m_pt_num.RebinY(2)

    m_pt_num.GetXaxis().SetTitle("M_{PNet} [GeV] / 5 GeV")
    m_pt_num.GetYaxis().SetTitle("p_{T} [GeV] / 20 GeV")
    h2_for_plotting = m_pt_num.Clone("h2_plotting")
    h2_for_plotting.GetXaxis().SetRangeUser(60,200)
    h2_for_plotting.GetZaxis().SetRangeUser(0,1)
    h2_for_plotting.Reset()
    eff2D           = r.TEfficiency(m_pt_num,m_pt_den)
    #Not needed
    # pt_num      = m_pt_num.ProjectionY("pt_num")
    # pt_den      = m_pt_den.ProjectionY("pt_den")
    # eff_pt      = r.TEfficiency(pt_num,pt_den)

    # m_num       = m_pt_num.ProjectionX("m_num")
    # m_den       = m_pt_den.ProjectionX("m_den")
    # m_pt_num_rwt= rwt2DNum(m_pt_num,eff_pt)
    # m_num_rwt   = m_pt_num_rwt.ProjectionX("m_num_rwt")

    #eff_m       = r.TEfficiency(m_num,m_den)
    #eff_m_rwt   = r.TEfficiency(m_num_rwt,m_den)


    #Plotting
    c = r.TCanvas("c","",1500,1500)
    c.cd()
    h2_for_plotting.Draw("colz")
    eff2D.Draw("colz same")
    c.SetLeftMargin(0.15);
    c.SetBottomMargin(0.10);
    r.gPad.RedrawAxis()
    c.SaveAs("2Deff.pdf")
    c.SaveAs("2Deff.png")

if __name__ == '__main__':

    for year in ["2016","2016APV","2017","2018"]:
        
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

    #     iFile = "mergedOutputs/SingleMuon{0}Total.root".format(year)
    #     trigEffFromFile(iFile,year,luminosity,writeToRoot=True)

    # for year in ["2017B","2017C","2017D","2017E","2017F"]:
    #     iFile = "mergedOutputs/SingleMuon{0}.root".format(year)
    #     trigEffFromFile(iFile,year,"",writeToRoot=False)

    testMassDependency("mergedOutputs/SingleMuonRunII.root")