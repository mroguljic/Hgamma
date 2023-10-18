import ROOT

import time, os
from optparse import OptionParser
from collections import OrderedDict

from TIMBER.Tools.Common import *
from TIMBER.Analyzer import *
import sys

def add4VecToSnapshot(options):

    #----Initialize RDataFrame-----#
    start_time = time.time()
    a = analyzer(options.input)
    import correctionlib
    correctionlib.register_pyroot_binding()
    CompileCpp("TIMBER/Framework/Hgamma_modules/Hgamma_Functions.cc") 
    
    evtColumns = VarGroup("Event columns")
    evtColumns.Add("bVec","LHEfourVector(nGenPart,GenPart_pdgId,GenPart_pt,GenPart_eta,GenPart_phi,GenPart_mass,5)")
    evtColumns.Add("bbarVec","LHEfourVector(nGenPart,GenPart_pdgId,GenPart_pt,GenPart_eta,GenPart_phi,GenPart_mass,-5)")
    evtColumns.Add("gammaVec","LHEfourVector(nGenPart,GenPart_pdgId,GenPart_pt,GenPart_eta,GenPart_phi,GenPart_mass,22)")
    evtColumns.Add("hVec","LHEfourVector(nGenPart,GenPart_pdgId,GenPart_pt,GenPart_eta,GenPart_phi,GenPart_mass,25)")


    a.Apply([evtColumns])
    a.Snapshot("all",options.output,'Events',saveRunChain=True)
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

(options, args) = parser.parse_args()
add4VecToSnapshot(options)
