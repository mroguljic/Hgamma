import ROOT as r
import os
import csv
import numpy as np


datasets            = ["ZGamma","WGamma","GJets","TTGJets","QCD","HgammaHyy","HgammaHZy"]
tplDir              = "/users/mrogul/Work/Hgamma/Hgamma/results/templates/tight_medium/RunII/scaled/"
values              = []
labelFlag           = False
additionalRegions   = ["T","M","F"]
totalBkg            = []
labels              = ["Sample"]
for sample in datasets:
    if ("Hgamma" in sample):
        datasetType = "sig"
    else:
        datasetType = "bkg"
        

    f           = r.TFile.Open(tplDir+sample+".root")
    h           = f.Get("{0}_cutflow_nom".format(sample))
    nBin        = h.GetNbinsX()
    totalBkgRow = []
    latexRow    = sample+" &"
    for i in range(1,nBin-2):
        if(h.GetXaxis().GetBinLabel(i)=="Trigger"):
            continue
        if(labelFlag==False):
            labels.append(h.GetXaxis().GetBinLabel(i))
        nCut     = h.GetBinContent(i)
        #nCut = nCut/(h.GetBinContent(1))
        if(nCut>10000):
            latexRow = latexRow +" {0:1.2e} &".format(nCut)
        else:   
            latexRow = latexRow +" {0:.1f} &".format(nCut)
        totalBkgRow.append(nCut)

    if(labelFlag==False):
        labelFlag=True
        print(" & ".join(labels), " \\\\")

    #Additional regions
    for region in additionalRegions:
        nCut = f.Get("{0}_H_m_pT_{1}__nominal".format(sample,region)).Integral()
        #nCut = nCut/(h.GetBinContent(1))
        if(nCut>10000):
            latexRow = latexRow +" {0:1.2e} &".format(nCut)
        else:   
            latexRow = latexRow +" {0:.1f} &".format(nCut)
        totalBkgRow.append(nCut)
    #-----#
    

    latexRow = latexRow[:-2]+" \\\\"
    print(latexRow)
    if(datasetType=="bkg"):
        if len(totalBkg)==0:
            totalBkg=totalBkgRow
        else:
            totalBkg = np.add(totalBkg,totalBkgRow)
    # a = latexRow.split(" & ")
    # print("{0} {1} {2}".format(a[0],a[-4],a[-5]))

latexTotalBkg = "Total_bkg "
for val in totalBkg:
    if(val>10000):
        latexTotalBkg += "& {0:1.2e} ".format(val)
    else:
        latexTotalBkg += "& {0:.2f} ".format(val)
print(latexTotalBkg+"\\\\")
