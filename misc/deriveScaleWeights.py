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

f 				= r.TFile.Open("zpt_scales_nom.root")
hNom  		 	= f.Get("ZGamma_pT__nominal")
nomYield        = hNom.Integral()
hNom, edges  	= hist2array(hNom,return_edges=True)


#nBins 			= hNom.GetNbinsX() 
# for i in range(1,nBins+1):
# 	nomVal  = hNom.GetBinContent(i)
# 	minDiff = 0.
# 	maxDiff = 0.
# 	for j in range(8):
# 		hTemp   = f.Get("ZGamma_pT_{0}".format(j))
# 		tempVal = hTemp.GetBinContent(i)
# 		relDiff = (tempVal - nomVal)/nomVal
# 		if relDiff>maxDiff:
# 			maxDiff = relDiff
# 		if relDiff<minDiff:
# 			minDiff = relDiff
# 	print(minDiff,maxDiff,(maxDiff-minDiff)/2.)

histos = [hNom]
labels = ["(1.0,1.0)","(0.5,0.5)","(1.0,0.5)","(2.0,0.5)","(0.5,1.0)","(2.0,1.0)","(0.5,2.0)","(1.0,2.0)","(2.0,2.0)"]

for j in range(8):
	hTemp   	= f.Get("ZGamma_pT_{0}".format(j))

	if(j==0 or j==7):
		print("{0:.2f}".format(hTemp.Integral()/nomYield))

	hTemp		= hist2array(hTemp,return_edges=False)
	histos.append(hTemp)


plt.style.use([hep.style.CMS])
f, ax = plt.subplots()
hep.histplot(histos, bins=edges[0],label=labels)
plt.legend(title="$(\mu_{F},\mu_{R})$ factors")

hep.cms.text("Work in progress",loc=0)
hep.cms.lumitext(text="2018 (13 TeV)", ax=ax, fontname=None, fontsize=None)

ax.set_xlim([300.,1000.])
#ax.set_ylim([0.,maxRpf*1.3])

plt.xlabel("H candidate $p_{T}$ [GeV]",horizontalalignment='right', x=1.0)
plt.ylabel("Events / 50 GeV [a.u.]",horizontalalignment='right', y=1.0)
#ax.yaxis.set_tick_params(which='minor', left=False)    
#ax.yaxis.set_tick_params(which='minor', right=False)    

foutName = "scaleVariations"
print("Saving {0}.pdf".format(foutName))
plt.savefig("{0}.pdf".format(foutName), bbox_inches='tight')
plt.savefig("{0}.png".format(foutName), bbox_inches='tight')

ax.set_yscale("log")
plt.savefig("{0}_log.pdf".format(foutName), bbox_inches='tight')
plt.savefig("{0}_log.png".format(foutName), bbox_inches='tight')

plt.cla()
plt.clf()