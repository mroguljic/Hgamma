#To be used with trees from event selection
import ROOT as r
import time, os
import sys
from optparse import OptionParser
from collections import OrderedDict

from TIMBER.Tools.Common import *
from TIMBER.Analyzer import *

def getNweighted(analyzer,isData):
    if not isData:
        nWeighted = analyzer.DataFrame.Sum("genWeight").GetValue()
    else:
        nWeighted = analyzer.DataFrame.Count().GetValue()
    return nWeighted

def getNCut(analyzer,cut,cutName):
    beforeNode = analyzer.GetActiveNode()
    analyzer.Cut(cutName,cut)
    nCut = analyzer.DataFrame.Sum("genWeight").GetValue()
    analyzer.SetActiveNode(beforeNode)
    return nCut

def getTaggingEfficiencies(analyzer,wpM,wpT):
    beforeNode = analyzer.GetActiveNode()
    nTot = analyzer.DataFrame.Sum("genWeight").GetValue()
    analyzer.Cut("Eff_M_cut","pnetHiggs>{0} && pnetHiggs<{1}".format(wpM,wpT))
    nM   = analyzer.DataFrame.Sum("genWeight").GetValue()
    analyzer.SetActiveNode(beforeNode)
    analyzer.Cut("Eff_T_cut","pnetHiggs>{0}".format(wpT))
    nT   = analyzer.DataFrame.Sum("genWeight").GetValue()
    effM = nM/nTot
    effT = nT/nTot
    analyzer.SetActiveNode(beforeNode)
    return effM, effT


parser = OptionParser()

parser.add_option('-i', '--input', metavar='IFILE', type='string', action='store',
                default   =   '',
                dest      =   'input',
                help      =   'A root file to analyze')
parser.add_option('-o', '--output', metavar='OFILE', type='string', action='store',
                default   =   'output.root',
                dest      =   'output',
                help      =   'Output file name.')
parser.add_option('-y', '--year', metavar='year', type='string', action='store',
                default   =   '2016',
                dest      =   'year',
                help      =   'Dataset year')
parser.add_option('-p', '--process', metavar='PROCESS', type='string', action='store',
                default   =   'X1600_Y100',
                dest      =   'process',
                help      =   'Process in the given file')
parser.add_option('-v','--var', metavar='variation', type='string', action='store',
                default   =   "nom",
                dest      =   'variation',
                help      =   'nom; jer, jes, jmr, jms + Up/Down')
parser.add_option('-m', metavar='mode', type='string', action='store',
                default   =   "RECREATE",
                dest      =   'mode',
                help      =   'RECREATE or UPDATE outputfile')
parser.add_option('-w', '--wp', metavar='working point', action="append", type=float,
                default   =   [],
                dest      =   'wps',
                help      =   'Working point')



(options, args) = parser.parse_args()
variation = options.variation
iFile = options.input

if not variation in iFile:
    if("je" in variation or "jm" in variation):
        iFile = iFile.replace(".root","_{0}.root".format(variation))
        print("{0} not in {1}, swapping input to {2}".format(variation,options.input,iFile))
    else:
        if not("nom" in iFile):
            iFile = iFile.replace(".root","_nom.root")


if variation=="nom":
    #Nominal tree processing will run all variations which do not require separate selection
    nomTreeFlag = True
else:
    nomTreeFlag = False


a               = analyzer(iFile)
print(options.wps)
taggerWpUp   = options.wps[0]
taggerWpLo   = options.wps[1]
year            = options.year
histos          = []
histGroups      = []
taggerBranch    = "pnetHiggs"
massVar         = "HiggsPnetMass"

if("data" in options.process.lower() or "singlephoton" in options.process.lower() or "egamma" in options.process.lower()):
    print("Running on data")
    isData=True
else:
    print("Running on MC")
    isData=False

#PNet SF part
CompileCpp('TIMBER/Framework/src/btagSFHandler.cc')
if(variation=="pnetUp"):
    pnetVar=2
elif(variation=="pnetDown"):
    pnetVar=1
else:
    pnetVar = 0

if("Hgamma" in options.process or "ZGamma" in options.process):
    eff_M, eff_T = getTaggingEfficiencies(a,taggerWpLo,taggerWpUp)
    print("{0} ParticleNet (M,T) efficiencies: ({1:.2f},{2:.2f})".format(options.process,eff_M,eff_T))
else:
    eff_M = 0.1
    eff_T = 0.1#placeholders, not needed if not applying pnet SF

CompileCpp('btagSFHandler btagHandler = btagSFHandler({%f,%f},{%f,%f},%s,%i);' %(taggerWpLo,taggerWpUp,eff_M,eff_T,'"{0}"'.format(year),pnetVar))#wps, efficiencies, year, var
a.Define("TaggerCat","btagHandler.createTaggingCategories(pnetHiggs)")

if("Hgamma" in options.process or "ZGamma" in options.process):
    a.Define("ScaledPnet","btagHandler.updateTaggingCategories(TaggerCat,Higgs_pt)")
else:
    a.Define("ScaledPnet","TaggerCat")

if isData:
    a.Define("genWeight","1")
#Workaround to include genWeight into template maker
genWCorr    = Correction('genW',"TIMBER/Framework/Zbb_modules/BranchCorrection.cc",corrtype='corr',mainFunc='evalCorrection')
a.AddCorrection(genWCorr, evalArgs={'val':'genWeight'})

if not isData:
    photIDFile   = "data/photIDSFs.root"
    photIDName   = "UL{0}_sf".format(year)

    a.Define("Gamma_pt_CutID","TMath::Min(Double_t(Gamma_pt),999.)")#Cut based SF are measured up to 1000 GeV
    photIDCorr   = Correction('photID',"TIMBER/Framework/src/HistLoader.cc",constructor=[photIDFile,photIDName],corrtype='weight',mainFunc='eval')
    a.AddCorrection(photIDCorr, evalArgs={'xval':'Gamma_eta','yval':'Gamma_pt_CutID','zval':0})

    if not "Hgamma" in options.process:
        #Private signal does not have proper pdf nor pu correction
        pdfCorr     = genWCorr.Clone("pdfUnc",newMainFunc="evalUncert",newType="uncert")
        puCorr      = genWCorr.Clone("puUnc",newMainFunc="evalWeight",newType="weight")
        a.AddCorrection(pdfCorr, evalArgs={'valUp':'Pdfweight__up','valDown':'Pdfweight__down'})
        a.AddCorrection(puCorr, evalArgs={'val':'Pileup__nom','valUp':'Pileup__up','valDown':'Pileup__down'})

    if "ZGamma" in options.process:
        NLOfile    = "data/ewk_corr.root"
        ewkName    = "zgamma_ewk"
        a.Define("genGammaPt_rescaled","TMath::Max(200.,TMath::Min(Double_t(genGammaPt),1000.))")#Weights applied in 200-2000 GeV gen V pt range
        NLOewkCorr = Correction('ewk_nlo',"TIMBER/Framework/src/HistLoader.cc",constructor=[NLOfile,ewkName],corrtype='corr',isClone=True,cloneFuncInfo=photIDCorr._funcInfo,isNewConstr=True)
        ISRcorr    = genWCorr.Clone("ISRunc",newMainFunc="evalUncert",newType="uncert")
        FSRcorr    = genWCorr.Clone("FSRunc",newMainFunc="evalUncert",newType="uncert")
        a.AddCorrection(NLOewkCorr, evalArgs={'xval':'genGammaPt_rescaled','yval':0,'zval':0})
        a.AddCorrection(ISRcorr, evalArgs={'valUp':'ISR__up','valDown':'ISR__down'})
        a.AddCorrection(FSRcorr, evalArgs={'valUp':'FSR__up','valDown':'FSR__down'})

    if(year=="2018"):
        hemCorr = genWCorr.Clone("hemCorrection")
        a.AddCorrection(hemCorr, evalArgs={'val':'HEM_drop__nom'})
    if(year!="2018"):
        prefireCorr = genWCorr.Clone("prefireUnc",newMainFunc="evalWeight",newType="weight")
        a.AddCorrection(prefireCorr, evalArgs={'val':'Prefire__nom','valUp':'Prefire__up','valDown':'Prefire__down'})

    #Trigger is negligible
    #trigFile   = "data/trig_eff_{0}.root".format(year)
    #a.Define("pt_for_trig","TMath::Min(Double_t(FatJet_pt0),999.)")#Trigger eff, measured up to 1000 GeV (well withing 100% eff. regime)
    #triggerCorr = Correction('trig',"TIMBER/Framework/Zbb_modules/TrigEff.cc",constructor=["{0}".format(trigFile),"trig_eff"],corrtype='weight')
    #a.AddCorrection(triggerCorr, evalArgs={'xval':'pt_for_trig','yval':0,'zval':0})

a.MakeWeightCols()

regionDefs      = [("T","ScaledPnet==2"),("F","ScaledPnet==0"),("M","ScaledPnet==1")]
regionYields    = {}

for region,cut in regionDefs:
    a.Define(region,cut)

checkpoint = a.GetActiveNode()

for region,cut in regionDefs:
    a.SetActiveNode(checkpoint)
    a.Cut("{0}_cut".format(region),cut)

    if nomTreeFlag:
        HmassHist   = r.TH1F('{0}_H_m_{1}'.format(options.process,region),';{0} [GeV];'.format(massVar),32,40,200)
        HmassPtHist = r.TH2F('{0}_H_m_pT_{1}'.format(options.process,region),';{0} [GeV];pT [GeV]'.format(massVar),32,40,200,17,300,2000)
        HptHist     = a.DataFrame.Histo1D(('{0}_H_pT_{1}'.format(options.process,region),';pT [GeV];',31,300,2000),"Higgs_pt","weight__nominal")
        GammaptHist = a.DataFrame.Histo1D(('{0}_Gamma_pT_{1}'.format(options.process,region),';pT [GeV];',31,300,2000),"Gamma_pt","weight__nominal")
        templates   = a.MakeTemplateHistos(HmassHist,[massVar])
        templates2D = a.MakeTemplateHistos(HmassPtHist,[massVar,"Higgs_pt"])
        histGroups.append(templates)
        histGroups.append(templates2D)
        histos.append(HptHist)
        histos.append(GammaptHist)
    else:
        #For jms/jmr/jes/jer trees, we don't need to calculate uncertainties on nominal trees
        template   = a.DataFrame.Histo1D(('{0}_H_m_{1}'.format(options.process,region),';{0} [GeV];'.format(massVar),32,40,200),massVar,"weight__nominal")
        template2D = a.DataFrame.Histo2D(('{0}_H_m_pT_{1}'.format(options.process,region),';{0} [GeV];pT [GeV]'.format(massVar),32,40,200,17,300,2000),massVar,"Higgs_pt","weight__nominal")
        histos.append(template)
        histos.append(template2D)

    regionYields[region] = getNweighted(a,isData)

#include histos from evt sel in the template file for nominal template
if(options.variation=="nom"):
    in_f = ROOT.TFile(iFile)
    for key in in_f.GetListOfKeys():
        h = key.ReadObj()
        hName = h.GetName()
        if(hName=="Events"):
            continue
        h.SetDirectory(0)
        if("cutflow" in hName.lower()):
            h.SetBinContent(h.GetNbinsX()-1,regionYields["T"])
            h.SetBinContent(h.GetNbinsX(),regionYields["F"])
        histos.append(h)

out_f = ROOT.TFile(options.output,options.mode)
out_f.cd()
for h in histos:
    h.SetName(h.GetName()+"_"+options.variation)
    h.Write()

for hGroup in histGroups:
    hGroup.Do("Write")

out_f.Close()

#a.PrintNodeTree('test.ps')
