import ROOT as r
import json
import sys
import re
import os

def normalizeProcess(process,year,inFile,outFile):
    h_dict = {}
    f = r.TFile.Open(inFile)
    print(process,inFile)
    json_file = open("xsecs.json")
    config = json.load(json_file)
    xsec    = config[year][process]["xsec"]
    luminosity  = config[year]["lumi"]
    sumGenW     = f.Get("{0}_cutflow_nom".format(process)).GetBinContent(1)
    nLumi       = xsec*luminosity
    scaling     = nLumi/sumGenW
    print("Scale: {0:.6f}".format(scaling))
    for key in f.GetListOfKeys():
        h = key.ReadObj()
        hName = h.GetName()
        h.Scale(scaling)
        h.SetDirectory(0)
        h_dict[hName] = h
    f.Close()

    f = r.TFile(outFile,"recreate")
    f.cd()
    for key in h_dict:
        histo = h_dict[key]
        histo.Write()
    f.Close()

def mergeSamples(inFiles,outFile,regexMatch,regexReplace):
    h_dict = {}
    print("Merging to {0}".format(outFile))
    for inFile in inFiles:
        print(inFile)
        f        = r.TFile.Open(inFile) 
        for key in f.GetListOfKeys():
            h = key.ReadObj()
            hName = h.GetName()
            h.SetDirectory(0)
            hKey = re.sub(regexMatch,regexReplace,hName,count=1)
            if not hKey in h_dict:
                h.SetName(hKey)
                h_dict[hKey] = h
            else:
                h_dict[hKey].Add(h)
        f.Close()
    f = r.TFile(outFile,"recreate")
    f.cd()
    for key in h_dict:
        histo = h_dict[key]
        histo.Write()
    f.Close()
    print("\n")

def lumiNormalization(wp="tight",tagger="ParticleNet"):

    #processes = ["ZGamma","WGamma","GJets200","GJets400","GJets600","GJetsHT200","GJetsHT400","GJetsHT600","TTGJets","Hgamma"
    #,"QCD500","QCD700","QCD1000","QCD1500","QCD2000"]
    processes = ["ZGamma","WGamma","GJets200","GJets400","GJets600","TTGJets","Hgamma","Hgamma_HZy","QCD500","QCD700","QCD1000","QCD1500","QCD2000"]
    for year in ['2016','2016APV','2017','2018']:
        print(year)
        nonScaledDir = "results/templates/{2}/{0}/{1}/nonScaled/".format(wp,year,tagger)
        lumiScaledDir = "results/templates/{2}/{0}/{1}/scaled/".format(wp,year,tagger)

        for proc in processes:
            nonScaledFile = "{0}/{1}.root".format(nonScaledDir,proc)
            if(os.path.isfile(nonScaledFile)):
                try:                 
                    normalizeProcess(proc,year,"{0}/{1}.root".format(nonScaledDir,proc),"{0}/{1}.root".format(lumiScaledDir,proc))
                except:
                    print("Couldn't normalize {0}".format(proc))
            else:
                print("{0} does not exist, skipping!".format(nonScaledFile))
        
        GJetsSamples = ["GJets200.root","GJets400.root","GJets600.root"]
        GJetsSamples = [lumiScaledDir+f for f in GJetsSamples if (os.path.isfile(os.path.join(lumiScaledDir, f)))]
        mergeSamples(GJetsSamples,"{0}/GJets{1}.root".format(lumiScaledDir,year[2:]),"GJets\d+_","GJets_")

        # GJetsHTSamples = ["GJetsHT200.root","GJetsHT400.root","GJetsHT600.root"]
        # GJetsHTSamples = [lumiScaledDir+f for f in GJetsHTSamples if (os.path.isfile(os.path.join(lumiScaledDir, f)))]
        # mergeSamples(GJetsHTSamples,"{0}/GJetsHT{1}.root".format(lumiScaledDir,year[2:]),"GJetsHT\d+_","GJetsHT_")

        QCDSamples = ["QCD500.root","QCD700.root","QCD1000.root","QCD1500.root","QCD2000.root"]
        QCDSamples = [lumiScaledDir+f for f in QCDSamples if (os.path.isfile(os.path.join(lumiScaledDir, f)))]
        mergeSamples(QCDSamples,"{0}/QCD{1}.root".format(lumiScaledDir,year[2:]),"QCD\d+_","QCD_")
        
        # ttSamples = ["TTHadronic.root","TTSemileptonic.root"]
        # ttSamples = [lumiScaledDir+f for f in ttSamples if (os.path.isfile(os.path.join(lumiScaledDir, f)))]
        # mergeSamples(ttSamples,"{0}/TTbar{1}.root".format(lumiScaledDir,year[2:]),"TTSemileptonic|TTHadronic","TTbar")


        SinglePhotonSamples = [nonScaledDir+f for f in os.listdir(nonScaledDir) if (os.path.isfile(os.path.join(nonScaledDir, f)) and "SinglePhoton" in f)]
        mergeSamples(SinglePhotonSamples,"{0}/SinglePhoton{1}.root".format(lumiScaledDir,year[2:]),"SinglePhoton201[0-9][a-zA-Z0-9]+_","data_obs_")

def lumiNormalizationCR(wp="tight",tagger="ParticleNet"):
    processes = ["ZJets400","ZJets600","ZJets800","WJets400","WJets600","WJets800"]
    for year in ['2016','2016APV','2017','2018']:
        print(year)
        nonScaledDir = "results/templates_CR/{2}/{0}/{1}/nonScaled/".format(wp,year,tagger)
        lumiScaledDir = "results/templates_CR/{2}/{0}/{1}/scaled/".format(wp,year,tagger)

        for proc in processes:
            nonScaledFile = "{0}/{1}.root".format(nonScaledDir,proc)
            if(os.path.isfile(nonScaledFile)):
                try:                 
                    normalizeProcess(proc,year,"{0}/{1}.root".format(nonScaledDir,proc),"{0}/{1}.root".format(lumiScaledDir,proc))
                except:
                    print("Couldn't normalize {0}".format(proc))
            else:
                print("{0} does not exist, skipping!".format(nonScaledFile))
        
        WJetsSamples = ["WJets400.root","WJets600.root","WJets800.root"]
        WJetsSamples = [lumiScaledDir+f for f in WJetsSamples if (os.path.isfile(os.path.join(lumiScaledDir, f)))]
        mergeSamples(WJetsSamples,"{0}/WJets{1}.root".format(lumiScaledDir,year[2:]),"[A-Z]Jets\d+_","WJets_")

        ZJetsSamples = ["ZJets400.root","ZJets600.root","ZJets800.root"]
        ZJetsSamples = [lumiScaledDir+f for f in ZJetsSamples if (os.path.isfile(os.path.join(lumiScaledDir, f)))]
        mergeSamples(ZJetsSamples,"{0}/ZJets{1}.root".format(lumiScaledDir,year[2:]),"[A-Z]Jets\d+_","ZJets_")

        JetHTSamples = [nonScaledDir+f for f in os.listdir(nonScaledDir) if (os.path.isfile(os.path.join(nonScaledDir, f)) and "JetHT" in f)]
        mergeSamples(JetHTSamples,"{0}/JetHT{1}.root".format(lumiScaledDir,year[2:]),"JetHT201[0-9][a-zA-Z0-9]+_","data_obs_")

def mergeRunII(wp,tagger):
    lumiScaledDir = "results/templates/{0}/{1}".format(tagger,wp)
    runIIDir      = "results/templates/{0}/{1}/RunII/scaled/".format(tagger,wp)
    if not os.path.exists(runIIDir):
        os.makedirs(runIIDir)
    os.system("hadd -f {0}/SinglePhoton.root {1}/201*/scaled/SinglePhoton*.root".format(runIIDir,lumiScaledDir))
    os.system("hadd -f {0}/ZGamma.root {1}/201*/scaled/ZGamma*.root".format(runIIDir,lumiScaledDir))
    os.system("hadd -f {0}/WGamma.root {1}/201*/scaled/WGamma*.root".format(runIIDir,lumiScaledDir))
    os.system("hadd -f {0}/GJets.root {1}/201*/scaled/GJets1*.root".format(runIIDir,lumiScaledDir))
    os.system("hadd -f {0}/TTGJets.root {1}/201*/scaled/TTGJets.root".format(runIIDir,lumiScaledDir))
    os.system("hadd -f {0}/QCD.root {1}/201*/scaled/QCD1?.root {1}/201*/scaled/QCD1?APV.root".format(runIIDir,lumiScaledDir))
    os.system("hadd -f {0}/Hgamma.root {1}/201*/scaled/Hgamma.root".format(runIIDir,lumiScaledDir))
    os.system("hadd -f {0}/Hgamma_HZy.root {1}/201*/scaled/Hgamma_HZy.root".format(runIIDir,lumiScaledDir))

def mergeRunIICR(wp,tagger):
    lumiScaledDir = "results/templates_CR/{0}/{1}".format(tagger,wp)
    runIIDir      = "results/templates_CR/{0}/{1}/RunII/scaled/".format(tagger,wp)
    if not os.path.exists(runIIDir):
        os.makedirs(runIIDir)
    os.system("hadd -f {0}/JetHT.root {1}/201*/scaled/JetHT*.root".format(runIIDir,lumiScaledDir))
    os.system("hadd -f {0}/ZJets.root {1}/201*/scaled/ZJets1*.root".format(runIIDir,lumiScaledDir))
    os.system("hadd -f {0}/WJets.root {1}/201*/scaled/WJets1*.root".format(runIIDir,lumiScaledDir))

if __name__ == '__main__':
    #lumiNormalization(wp="tight_medium",tagger="/")
    #mergeRunII("tight_medium","/")
    lumiNormalizationCR(wp="tight_medium",tagger="/")
    mergeRunIICR(wp="tight_medium",tagger="/")
