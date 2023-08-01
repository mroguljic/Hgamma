import ROOT

import time, os
from optparse import OptionParser
from collections import OrderedDict

from TIMBER.Tools.Common import *
from TIMBER.Analyzer import *
import sys

def eventSelection(options):

    #----Initialize RDataFrame-----#
    start_time = time.time()
    year = options.year
    a = analyzer(options.input)

    # #For testing only, faster execution
    # nToRun = 1000
    # small_rdf = a.GetActiveNode().DataFrame.Range(nToRun) 
    # small_node = Node('small',small_rdf)
    # a.SetActiveNode(small_node)

    #nProc = a.DataFrame.Count().GetValue()
    CompileCpp("TIMBER/Framework/Zbb_modules/Zbb_Functions.cc") 
    CompileCpp("TIMBER/Framework/Zbb_modules/helperFunctions.cc") 

    histos          = []

    if(options.year=="2016APV"):
       deepJetM    = 0.2598
    elif(options.year=="2016"):
       deepJetM    = 0.2489
    elif(options.year=="2017"):
       deepJetM    = 0.3040
    elif(options.year=="2018"):
       deepJetM    = 0.2783
    else:
        print("Year not supported")
        sys.exit()


    #----------Triggers------------#
    #beforeTrigCheckpoint    = a.GetActiveNode()
    if(year=="2016" or year=="2016APV"):
        refTrigger          = ["HLT_IsoMu24","HLT_IsoTkMu24","HLT_Mu50"]
        triggerList         = ["HLT_AK8DiPFJet280_200_TrimMass30","HLT_AK8PFHT650_TrimR0p1PT0p03Mass50",
            "HLT_AK8PFHT700_TrimR0p1PT0p03Mass50","HLT_PFHT800","HLT_PFHT900","HLT_AK8PFJet360_TrimMass30","HLT_AK8PFJet450"]
    elif(year=="2017"):
        refTrigger          = ["HLT_IsoMu27","HLT_Mu50"]
        triggerList         =["HLT_PFHT1050","HLT_AK8PFJet400_TrimMass30","HLT_AK8PFJet420_TrimMass30",
        "HLT_AK8PFHT800_TrimMass50","HLT_PFJet500","HLT_AK8PFJet500"]
    elif(year=="2018"):
        refTrigger          = ["HLT_IsoMu24","HLT_Mu50"]
        triggerList          = ["HLT_PFHT1050","HLT_AK8PFJet400_TrimMass30","HLT_AK8PFJet420_TrimMass30","HLT_AK8PFHT800_TrimMass50","HLT_PFJet500","HLT_AK8PFJet500"]
    triggersStringTarget       = a.GetTriggerString(triggerList)    
    triggersStringRef       = a.GetTriggerString(refTrigger)    
    #------------------------------#


    #----------Selection-----------#
    a.Cut("nFatJet","nFatJet>1")

    evtColumns = VarGroup("Event columns")
    evtColumns.Add("FatJet_pt0","FatJet_pt[0]")
    evtColumns.Add("FatJet_pt1","FatJet_pt[1]")
    evtColumns.Add("FatJet_eta0","FatJet_eta[0]")
    evtColumns.Add("FatJet_phi0","FatJet_phi[0]")
    evtColumns.Add("JetPnetMass",'FatJet_particleNet_mass[0]')


    evtColumns.Add("topVetoFlag","topVeto(FatJet_eta0,FatJet_phi0,nJet,Jet_eta,{0},Jet_phi,Jet_pt,Jet_btagDeepFlavB,{1})".format(2.4,deepJetM))

    a.Apply([evtColumns])
    
    a.Cut("pT","FatJet_pt0>450 && FatJet_pt1>200")
    a.Cut("Eta","abs(FatJet_eta[0])<2.4 && abs(FatJet_eta[1])<2.4")
    a.Cut("JetPnetMassCut","JetPnetMass>40")
    a.Cut("topVeto","topVetoFlag==0")
    a.Cut("photonVeto","!(nPhoton>0 && Photon_cutBased[0]==3 && abs(Photon_eta[0])<2.4 && Photon_pt[0]>300)")


    #-------MET filters------------#
    MetFilters = ["Flag_BadPFMuonFilter","Flag_EcalDeadCellTriggerPrimitiveFilter","Flag_HBHENoiseIsoFilter","Flag_HBHENoiseFilter",
    "Flag_globalSuperTightHalo2016Filter","Flag_goodVertices","Flag_BadPFMuonDzFilter","Flag_eeBadScFilter","Flag_ecalBadCalibFilter"]
    if(year=="2017" or year=="2018"):
        MetFilters.append("Flag_ecalBadCalibFilter")
    MetFiltersString = a.GetFlagString(MetFilters)
    if MetFiltersString:#RDF crashes if METstring is empty
        a.Cut("MET_Filters",MetFiltersString)
    #------------------------------#

    a.Cut("Reference Triggers",triggersStringRef) #Only applying reference trigger for denominator


    h2_denom   = a.GetActiveNode().DataFrame.Histo2D(('data_m_pT_denominator',';M_{PNet} [GeV] / 1 GeV;p_{T} [GeV] / 10 GeV;',160,40,200,55,450,1000),'JetPnetMass','FatJet_pt0')
    h_pT_denom = a.GetActiveNode().DataFrame.Histo1D(('data_pT_denominator',';Leading jet pT [GeV]; Events/10 GeV;',55,450,1000),"FatJet_pt0")
    
    a.Cut("JetHT triggers",triggersStringTarget)
        
    h2_numer   = a.GetActiveNode().DataFrame.Histo2D(('data_m_pT_numerator',';M_{PNet} [GeV] / 1 GeV;p_{T} [GeV] / 10 GeV;',160,40,200,55,450,1000),'JetPnetMass','FatJet_pt0')
    h_pT_numer = a.GetActiveNode().DataFrame.Histo1D(('data_pT_numerator',';Leading jet pT [GeV]; Events/10 GeV;',55,450,1000),"FatJet_pt0")

    histos.append(h2_denom)
    histos.append(h_pT_denom)
    histos.append(h2_numer)
    histos.append(h_pT_numer)
    #------------------------------#
    outputFile      = options.output
    out_f = ROOT.TFile(outputFile,"RECREATE")
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
parser.add_option('-y', '--year', metavar='year', type='string', action='store',
                default   =   '2016',
                dest      =   'year',
                help      =   'Dataset year')

(options, args) = parser.parse_args()
eventSelection(options)
