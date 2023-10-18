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
    a = analyzer(options.input)

    # #For testing only, faster execution
    # nToRun = 10000
    # small_rdf = a.GetActiveNode().DataFrame.Range(nToRun) 
    # small_node = Node('small',small_rdf)
    # a.SetActiveNode(small_node)


    histos          = []
    import correctionlib
    correctionlib.register_pyroot_binding()
    CompileCpp("TIMBER/Framework/Hgamma_modules/Hgamma_Functions.cc") 
    CompileCpp("TIMBER/Framework/Zbb_modules/Zbb_Functions.cc") 
    CompileCpp("TIMBER/Framework/Zbb_modules/helperFunctions.cc") 

    #----------Selection-----------#
    a.Cut("JetAndPhoton","nFatJet>0 && nPhoton>0")
    a.Define("Hidx","leadingNonGammaAK8Idx(nFatJet,FatJet_eta,FatJet_phi,Photon_eta[0],Photon_phi[0])")
    a.Cut("HidxCut","Hidx>-1")
    a.Cut("EtaCut","abs(FatJet_eta[Hidx])<2.4 && abs(Photon_eta[0])<2.4")

    a.Cut("GammaID","Photon_cutBased[0]==3")#cut-based ID Fall17V2 Tight
    #a.Cut("GammaID","Photon_mvaID_WP80[0]==1")#MVA ID Fall17V2

    evtColumns = VarGroup("Event columns")
    evtColumns.Add("Higgs_pt","FatJet_pt[Hidx]")
    evtColumns.Add("Higgs_mass","FatJet_particleNet_mass[Hidx]")
    a.Apply([evtColumns])

    a.Define("Gamma_pt_up","Photon_pt[0]*(1+Photon_dEsigmaUp[0])")
    a.Define("Gamma_pt_dn","Photon_pt[0]*(1+Photon_dEsigmaDown[0])")
    a.Define("Gamma_pt_nom","Photon_pt[0]")

    h_nom  = a.GetActiveNode().DataFrame.Histo1D(('pt_precut_nom',';Leading photon pT [GeV]; Events/10 GeV;',100,0,1000),"Gamma_pt_nom")
    h_up   = a.GetActiveNode().DataFrame.Histo1D(('pt_precut_up',';Leading photon pT [GeV]; Events/10 GeV;',100,0,1000),"Gamma_pt_up")
    h_down = a.GetActiveNode().DataFrame.Histo1D(('pt_precut_down',';Leading photon pT [GeV]; Events/10 GeV;',100,0,1000),"Gamma_pt_dn")
    histos.extend([h_nom,h_up,h_down])

    a.Cut("pT","Higgs_pt>300 && Higgs_mass>50")

    h_nom  = a.GetActiveNode().DataFrame.Histo1D(('pt_nom',';Leading photon pT [GeV]; Events/10 GeV;',100,0,1000),"Gamma_pt_nom")
    h_up   = a.GetActiveNode().DataFrame.Histo1D(('pt_up',';Leading photon pT [GeV]; Events/10 GeV;',100,0,1000),"Gamma_pt_up")
    h_down = a.GetActiveNode().DataFrame.Histo1D(('pt_down',';Leading photon pT [GeV]; Events/10 GeV;',100,0,1000),"Gamma_pt_dn")

    histos.extend([h_nom,h_up,h_down])

    out_f = ROOT.TFile("testPhotonEr.root","RECREATE")
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

(options, args) = parser.parse_args()
eventSelection(options)
