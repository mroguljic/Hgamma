#!/usr/bin/env python

import numpy as np
import gzip
from lhereader import LHEReader
import ROOT
import os

ROOT.gROOT.SetBatch(ROOT.kTRUE)

def MakePlots(proc, vals):

  lhef = os.path.join(os.getcwd(), vals['fin'])
  fout = vals['fin'].replace('lhe', 'root')

  fout = ROOT.TFile(fout, 'recreate')
  fout.cd()

  h_pt_h     = ROOT.TH1D('h_pt_h', ';p_{T} (h) [GeV];;',50, 0.,500)
  h_pt_gamma = ROOT.TH1D('h_pt_gamma', ';p_{T} (h) [GeV];;',50, 0.,500)

  reader = LHEReader(lhef)
  
  for iev, event in enumerate(reader):
  
      if iev > 10000:
        break
  
      hs     = list(filter(lambda x: abs(x.pdgid)== 25, event.particles))
      gammas = list(filter(lambda x: abs(x.pdgid)== 22, event.particles))
  
      if len(hs) > 0 and len(gammas) > 0:
        h_pt_h.Fill(hs[0].p4().pt)
        h_pt_gamma.Fill(gammas[0].p4().pt)
      else:
        print("No gamma and Higgs in event")
  
  fout.Write()
  fout.Close()
  

def main():

  procs = {
      'HGamma_SM': {
        'fin': 'Hgamma_SM_LO.lhe'
        }
      }

  for proc, vals in procs.items():
    MakePlots(proc, vals)

if __name__ == "__main__":
  main()
