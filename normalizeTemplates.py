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

    processes = ["ZGamma","WGamma","GJets200","GJets400","GJets600","TTGJets","Hgamma"]
    #for year in ['2016','2016APV','2017','2018']:
    for year in ['2017']:
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

        SinglePhotonSamples = [nonScaledDir+f for f in os.listdir(nonScaledDir) if (os.path.isfile(os.path.join(nonScaledDir, f)) and "SinglePhoton" in f)]
        mergeSamples(SinglePhotonSamples,"{0}/SinglePhoton{1}.root".format(lumiScaledDir,year[2:]),"SinglePhoton201[0-9][a-zA-Z]+_","data_obs_")


if __name__ == '__main__':
    lumiNormalization(wp="tight",tagger="/")