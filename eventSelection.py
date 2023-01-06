import ROOT

import time, os
from optparse import OptionParser
from collections import OrderedDict

from TIMBER.Tools.Common import *
from TIMBER.Analyzer import *
import sys




def eventSelection(options,isData):

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
    histos          = []

    #----------Triggers------------#
    beforeTrigCheckpoint    = a.GetActiveNode()
    if(year=="2016" or year=="2016APV"):
        triggerList         = ["HLT_Photon200"]
    elif(year=="2017"):
        triggerList         =["HLT_Photon200"]
    elif(year=="2018"):
       triggerList          =["HLT_Photon200"]
    a.Cut("Triggers","HLT_Photon200==1")
    nTrig                   = getNweighted(a,isData)
    #------------------------------#



    #----------Selection-----------#
    a.Cut("nFatJet","nFatJet>0 && nPhoton>0")
    a.Cut("ID","FatJet_jetId[0]>1")#bit 1 is loose, bit 2 is tight, bit3 is tightlepVeto, we select tight
    nJetID = getNweighted(a,isData)


    a.Cut("Eta","abs(FatJet_eta[0])<2.4 && abs(Photon_eta[0])<2.4")
    nEta = getNweighted(a,isData)

    evtColumns = VarGroup("Event columns")
    evtColumns.Add('H_pt','FatJet_pt[0]')
    evtColumns.Add('H_mass','FatJet_particleNet_mass[0]')
    evtColumns.Add('gamma_pt','Photon_pt[0]')

    a.Apply([evtColumns])

    # a.Cut("H_pt","H_pt>300")
    # a.Cut("gamma_pt","H_pt>200")
    # npT = getNweighted(a,isData)

    # a.Cut("H_mass","JetPnetMass>50")
    # nJetMass = getNweighted(a,isData)

    a.Define("H_pnet","FatJet_particleNetMD_Xbb[0]/(FatJet_particleNetMD_Xbb[0]+FatJet_particleNetMD_QCD[0])")
    #------------------------------#


    h_H_mass    = a.DataFrame.Histo1D(('{0}_H_mass'.format(options.process),';H Pnet mass',60,0.,300.),"H_mass")
    h_H_pt      = a.DataFrame.Histo1D(('{0}_H_pt'.format(options.process),';H pT',100,0.,1000.),"H_pt")
    h_gamma_pt  = a.DataFrame.Histo1D(('{0}_gamma_pt'.format(options.process),';Photon pT',100,0.,1000.),"gamma_pt")


    histos.append(h_gamma_pt)
    histos.append(h_H_mass)
    histos.append(h_H_pt)
    #--------Store Output----------#

    snapshotColumns = ["pnet0","H_pt","gamma_pt","H_mass","H_pnet","nPhoton","nFatJet"]
    outputFile      = options.output.replace(".root","_{0}.root".format(options.variation))



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
eventSelection(options,False)
