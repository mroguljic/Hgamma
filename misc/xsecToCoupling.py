import ROOT as r
from numpy import sqrt

def xsecToHyy(xsec):
	#Eq. 32 in https://arxiv.org/pdf/2109.13363.pdf
	#Returns either limit on g_2^yy or g_4^yy (they give the same result)
	#Other three couplings set to zero
	xsecRef = 1.33*10**4 #in fb
	Hyy     = sqrt(xsec/(xsecRef*0.553))
	return Hyy

def xsecToHZy(xsec):
	#Eq. 32 in https://arxiv.org/pdf/2109.13363.pdf
	#Returns either limit on g_2^Zy or g_4^Zy (they give the same result)
	#Other three couplings set to zero
	xsecRef = 1.33*10**4 #in fb
	HZy     = sqrt(xsec/xsecRef)
	return HZy

HZylimitFileName 	= "limitsHZy.root"
HyylimitFileName 	= "limitsHyy.root"
hbbBR           	= 0.58

limitFile 			= r.TFile.Open(HZylimitFileName)
limitTree  			= limitFile.Get("limit")
print("HZy")
for entry in limitTree:
	xsecLimitInFb = entry.limit*10./hbbBR #signal is normalized to xsec = 10fb, BR is included
	HZy           = xsecToHZy(xsecLimitInFb)	
	print("{0:.1f} {1:.4f}".format(xsecLimitInFb,HZy))
limitFile.Close()

limitFile 			= r.TFile.Open(HyylimitFileName)
limitTree  			= limitFile.Get("limit")
print("Hyy")
for entry in limitTree:
	xsecLimitInFb = entry.limit*10./hbbBR #signal is normalized to xsec = 10fb, BR is included
	Hyy           = xsecToHyy(xsecLimitInFb)	
	print("{0:.1f} {1:.4f}".format(xsecLimitInFb,Hyy))
limitFile.Close()