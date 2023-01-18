import ROOT

import time, os
from optparse import OptionParser
from collections import OrderedDict

from TIMBER.Tools.Common import *
from TIMBER.Analyzer import *
import sys


def getNweighted(analyzer,isData):
    if not isData:
        nWeighted = analyzer.DataFrame.Sum("genWeight").GetValue()
    else:
        nWeighted = analyzer.DataFrame.Count().GetValue()
    return nWeighted

def ptCorrectorString(variation,isData):
    if isData:
        jetCorrector  = '{FatJet_JES_nom}'
    else:
        if(variation=="jesUp"):
            jetCorrector  = '{FatJet_JES_up,FatJet_JER_nom}'
        elif(variation=="jesDown"):
            jetCorrector  = '{FatJet_JES_down,FatJet_JER_nom}'
        elif(variation=="jerUp"):
            jetCorrector  = '{FatJet_JES_nom,FatJet_JER_up}'
        elif(variation=="jerDown"):
            jetCorrector  = '{FatJet_JES_nom,FatJet_JER_down}'
        else:
            jetCorrector  = '{FatJet_JES_nom,FatJet_JER_nom}'
    return jetCorrector

def massCorrectorString(variation,isData):
    if isData:
        jetCorrector  = '{FatJet_JES_nom}'
    else:
        if(variation=="jesUp"):
            jetCorrector  = '{FatJet_JES_up,FatJet_JER_nom,FatJet_JMS_nom,FatJet_JMR_nom}'
        elif(variation=="jesDown"):
            jetCorrector  = '{FatJet_JES_down,FatJet_JER_nom,FatJet_JMS_nom,FatJet_JMR_nom}'
        elif(variation=="jerUp"):
            jetCorrector  = '{FatJet_JES_nom,FatJet_JER_up,FatJet_JMS_nom,FatJet_JMR_nom}'
        elif(variation=="jerDown"):
            jetCorrector  = '{FatJet_JES_nom,FatJet_JER_down,FatJet_JMS_nom,FatJet_JMR_nom}'
        elif(variation=="jmsUp"):
            jetCorrector  = '{FatJet_JES_nom,FatJet_JER_nom,FatJet_JMS_up,FatJet_JMR_nom}'
        elif(variation=="jmsDown"):
            jetCorrector  = '{FatJet_JES_nom,FatJet_JER_nom,FatJet_JMS_down,FatJet_JMR_nom}'
        elif(variation=="jmrUp"):
            jetCorrector  = '{FatJet_JES_nom,FatJet_JER_nom,FatJet_JMS_nom,FatJet_JMR_up}'
        elif(variation=="jmrDown"):
            jetCorrector  = '{FatJet_JES_nom,FatJet_JER_nom,FatJet_JMS_nom,FatJet_JMR_down}'
        else:
            jetCorrector  = '{FatJet_JES_nom,FatJet_JER_nom,FatJet_JMS_nom,FatJet_JMR_nom}'

    return jetCorrector

def eventSelection(options):

    #----Initialize RDataFrame-----#
    start_time = time.time()
    year = options.year
    a = analyzer(options.input)

    #For testing only, faster execution
    # nToRun = 10000
    # small_rdf = a.GetActiveNode().DataFrame.Range(nToRun) 
    # small_node = Node('small',small_rdf)
    # a.SetActiveNode(small_node)

    runNumber = a.DataFrame.Range(1).AsNumpy(["run"])#check run number to see if data
    if(runNumber["run"][0]>10000):
        isData=True
        a.Define("genWeight","1")
        print("Running on data")
    else:
        isData=False
        print("Running on MC")
    nProc = a.genEventSumw

    CompileCpp("TIMBER/Framework/Hgamma_modules/Hgamma_Functions.cc") 
    #CompileCpp("TIMBER/Framework/Zbb_modules/helperFunctions.cc") 

    ptCorrector     = ptCorrectorString(options.variation,isData) 
    massCorrector   = massCorrectorString(options.variation,isData) 

    histos          = []
    nSkimmed = getNweighted(a,isData)

    #-------MET filters------------#
    MetFilters = ["Flag_BadPFMuonFilter","Flag_EcalDeadCellTriggerPrimitiveFilter","Flag_HBHENoiseIsoFilter","Flag_HBHENoiseFilter",
    "Flag_globalSuperTightHalo2016Filter","Flag_goodVertices","Flag_BadPFMuonDzFilter","Flag_eeBadScFilter","Flag_ecalBadCalibFilter"]
    if(year=="2017" or year=="2018"):
        MetFilters.append("Flag_ecalBadCalibFilter")
    MetFiltersString = a.GetFlagString(MetFilters)
    if MetFiltersString:#RDF crashes if METstring is empty
        a.Cut("MET_Filters",MetFiltersString)
    #------------------------------#


    #----------Triggers------------#
    beforeTrigCheckpoint    = a.GetActiveNode()
    if(year=="2016" or year=="2016APV"):
        triggerList         = ["HLT_Photon200"]
    elif(year=="2017"):
        triggerList         =["HLT_Photon200"]
    elif(year=="2018"):
       triggerList          =["HLT_Photon200"]
    triggersStringAll       = a.GetTriggerString(triggerList)    
    if(isData): #Only applying trigger to data, will apply trigger turn-on to MC
        a.Cut("Triggers",triggersStringAll)
    nTrig                   = getNweighted(a,isData)
    #------------------------------#


    #----------Selection-----------#
    a.Cut("JetAndPhoton","nFatJet>0 && nPhoton>0")
    a.Define("Hidx","leadingNonGammaAK8Idx(nFatJet,FatJet_eta,FatJet_phi,Photon_eta[0],Photon_phi[0])")
    a.Cut("HidxCut","Hidx>-1")
    nJetGamma = getNweighted(a,isData)
    a.Cut("EtaCut","abs(FatJet_eta[Hidx])<2.4 && abs(Photon_eta[0])<2.4")
    nEta = getNweighted(a,isData)

    evtColumns = VarGroup("Event columns")
    evtColumns.Add('FatJet_pt_corr','hardware::MultiHadamardProduct(FatJet_pt,%s)'%ptCorrector)
    evtColumns.Add('FatJet_msoftdrop_corr','hardware::MultiHadamardProduct(FatJet_msoftdrop,%s)'%massCorrector)
    evtColumns.Add('FatJet_mpnet_corr','hardware::MultiHadamardProduct(FatJet_particleNet_mass,%s)'%massCorrector)
    evtColumns.Add("Higgs_pt","FatJet_pt_corr[Hidx]")
    evtColumns.Add("Higgs_eta","FatJet_eta[Hidx]")
    evtColumns.Add("Higgs_phi","FatJet_phi[Hidx]")
    evtColumns.Add("HiggsSDMass",'FatJet_msoftdrop_corr[Hidx]')
    evtColumns.Add("HiggsPnetMass",'FatJet_mpnet_corr[Hidx]')
    evtColumns.Add("Gamma_pt","Photon_pt[0]")
    evtColumns.Add("Gamma_eta","Photon_eta[0]")
    evtColumns.Add("Gamma_phi","Photon_phi[0]")

    a.Apply([evtColumns])

    a.Cut("pT","Higgs_pt>300 && Gamma_pt>300")
    npT = getNweighted(a,isData)

    a.Cut("JetPnetMassCut","HiggsPnetMass>50")
    nJetMass = getNweighted(a,isData)

    a.Define("pnetHiggs","FatJet_particleNetMD_Xbb[Hidx]/(FatJet_particleNetMD_Xbb[Hidx]+FatJet_particleNetMD_QCD[Hidx])")
    #--------Store Output----------#

    snapshotColumns = ["pnetHiggs","Higgs_pt","Gamma_pt","HiggsSDMass","HiggsPnetMass","PV_npvsGood","nFatJet","nPhoton"]
    outputFile      = options.output.replace(".root","_{0}.root".format(options.variation))

    if not isData:
        snapshotColumns.extend(['Pileup__nom','Pileup__up','Pileup__down','Pdfweight__up','Pdfweight__down','Pileup_nTrueInt','genWeight'])
        a.Define("ISR__up","PSWeight[2]")
        a.Define("ISR__down","PSWeight[0]")
        a.Define("FSR__up","PSWeight[3]")
        a.Define("FSR__down","PSWeight[1]")
        snapshotColumns.append("ISR__up")
        snapshotColumns.append("ISR__down")
        snapshotColumns.append("FSR__up")
        snapshotColumns.append("FSR__down")

        if year=="2018":
            snapshotColumns.append("HEM_drop__nom")

        if "2016" in year or "2017" in year:
            snapshotColumns.extend(['Prefire__nom','Prefire__up','Prefire__down'])


    a.Snapshot(snapshotColumns,outputFile,'Events',saveRunChain=False)

    cutFlowVars         = [nProc,nSkimmed,nTrig,nJetGamma,nEta,npT,nJetMass]
    cutFlowLabels       = ["Processed","Skimmed","Trigger","JetPlusGamma","Eta","pT","JetMass","pass","fail"]#tagging bins will be filled out in template making
    nCutFlowlabels      = len(cutFlowLabels)
    hCutFlow            = ROOT.TH1F('{0}_cutflow'.format(options.process),"Number of events after each cut",nCutFlowlabels,0.5,nCutFlowlabels+0.5)
    for i,label in enumerate(cutFlowLabels):
        hCutFlow.GetXaxis().SetBinLabel(i+1, label)

    for i,var in enumerate(cutFlowVars):
        hCutFlow.AddBinContent(i+1,var)

    histos.append(hCutFlow)


    out_f = ROOT.TFile(outputFile,"UPDATE")
    out_f.cd()
    for h in histos:
        h.Write()
    out_f.Close()
    #------------------------------#


    #a.PrintNodeTree('node_tree.dot',verbose=True)
    print("Total time: "+str((time.time()-start_time)/60.) + ' min')



parser = OptionParser()

parser.add_option('-i', '--input', metavar='IFILE', type='string', action='store',
                default   =   '',
                dest      =   'input',
                help      =   'A root file or text file with multiple root file locations to analyze')
parser.add_option('-o', '--output', metavar='OFILE', type='string', action='store',
                default   =   'output.root',
                dest      =   'output',
                help      =   'Output file name.')
parser.add_option('-p', '--process', metavar='PROCESS', type='string', action='store',
                default   =   'ZJets400',
                dest      =   'process',
                help      =   'Process in the given MC file')
parser.add_option('-y', '--year', metavar='year', type='string', action='store',
                default   =   '2016',
                dest      =   'year',
                help      =   'Dataset year')
parser.add_option('--var', metavar='variation', type='string', action='store',
                default   =   "nom",
                dest      =   'variation',
                help      =   'jmrUp/Down, jmsUp/Down, jesUp/Down, jerUp/Down')

(options, args) = parser.parse_args()
eventSelection(options)
