import ROOT as r

#ZGamma
f = r.TFile.Open("/users/mrogul/Work/Hgamma/Hgamma/results/templates/tight_medium/RunII/scaled/ZGamma.root")
yieldT = f.Get("ZGamma_H_m_pT_T__nominal").Integral()
yieldM = f.Get("ZGamma_H_m_pT_M__nominal").Integral()
yieldF = f.Get("ZGamma_H_m_pT_F__nominal").Integral()
yieldTot = yieldF+yieldM+yieldT

effT = yieldT/yieldTot
effM = yieldM/yieldTot

print("rate_F rateParam F_* ZGamma (1-{0:.3f}*@0-{1:.3f}*@1)/(1-{0:.3f}-{1:.3f}) SF_T,SF_M".format(effT,effM))
f.Close()

#ZGamma
f = r.TFile.Open("/users/mrogul/Work/Hgamma/Hgamma/results/templates_CR/tight_medium/RunII/scaled/ZJets.root")
yieldT = f.Get("ZJets_H_m_pT_CR_T__nominal").Integral()
yieldM = f.Get("ZJets_H_m_pT_CR_M__nominal").Integral()
yieldF = f.Get("ZJets_H_m_pT_CR_F__nominal").Integral()
yieldTot = yieldF+yieldM+yieldT

effT = yieldT/yieldTot
effM = yieldM/yieldTot

print("rate_CR_F rateParam CR_F_* ZJets (1-{0:.3f}*@0-{1:.3f}*@1)/(1-{0:.3f}-{1:.3f}) SF_T,SF_M".format(effT,effM))