import ROOT as r
import csv

hNom        = r.TH1F("nominal","",23,80,1000)
hCorrUp     = r.TH1F("up","",23,80,1000)
hCorrDn     = r.TH1F("down","",23,80,1000)
hUncertUp   = r.TH1F("uncUp","",23,80,1000)
hUncertDn   = r.TH1F("uncDn","",23,80,1000)

with open('corrections.csv', 'r') as csvfile:
    csvreader = csv.reader(csvfile, delimiter=',', quotechar='|')
    for i,row in enumerate(csvreader):
        value  = float(row[1])/100. #percentage
        uncert = value**2
        hNom.SetBinContent(i+1,1.+value)        
        hCorrUp.SetBinContent(i+1,1.+value+uncert)
        hCorrDn.SetBinContent(i+1,1.+value-uncert)
        hUncertUp.SetBinContent(i+1,1.+uncert)
        hUncertDn.SetBinContent(i+1,1.-uncert)

f = r.TFile.Open("ewk_zgamma.root","RECREATE")
f.cd()
hNom.Write()
hCorrUp.Write()
hCorrDn.Write()
hUncertUp.Write()
hUncertDn.Write()
f.Close()