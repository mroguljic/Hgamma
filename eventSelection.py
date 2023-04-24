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
            #jetCorrector  = '{FatJet_JES_nom,FatJet_JER_nom,FatJet_JMS_up,FatJet_JMR_nom}'
            #Applying a custom jms uncertainty in the event loop
            jetCorrector  = '{FatJet_JES_nom,FatJet_JER_nom,FatJet_JMS_nom,FatJet_JMR_nom}'
        elif(variation=="jmsDown"):
            #Applying a custom jms uncertainty in the event loop
            #jetCorrector  = '{FatJet_JES_nom,FatJet_JER_nom,FatJet_JMS_down,FatJet_JMR_nom}'
            jetCorrector  = '{FatJet_JES_nom,FatJet_JER_nom,FatJet_JMS_nom,FatJet_JMR_nom}'
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
    CompileCpp("TIMBER/Framework/Zbb_modules/Zbb_Functions.cc") 
    CompileCpp("TIMBER/Framework/Zbb_modules/helperFunctions.cc") 

    if(options.year=="2016APV"):
       deepCsvM    = 0.6001
    elif(options.year=="2016"):
       deepCsvM    = 0.5847
    elif(options.year=="2017"):
       deepCsvM    = 0.4506
    elif(options.year=="2018"):
       deepCsvM    = 0.4168 
    else:
        print("Year not supported")
        sys.exit()

    ptCorrector     = ptCorrectorString(options.variation,isData) 
    massCorrector   = massCorrectorString(options.variation,isData) 

    histos          = []
    nSkimmed = getNweighted(a,isData)

    #-------MET filters------------#
    MetFilters = ["Flag_BadPFMuonFilter","Flag_EcalDeadCellTriggerPrimitiveFilter","Flag_HBHENoiseIsoFilter","Flag_HBHENoiseFilter",
    "Flag_globalSuperTightHalo2016Filter","Flag_goodVertices","Flag_BadPFMuonDzFilter","Flag_eeBadScFilter"]
    if(year=="2017" or year=="2018"):
        MetFilters.append("Flag_ecalBadCalibFilter")
    MetFiltersString = a.GetFlagString(MetFilters)
    if MetFiltersString:#RDF crashes if METstring is empty
        a.Cut("MET_Filters",MetFiltersString)
    #------------------------------#
    
    #b-tag reshaping
    if(options.variation == "sfDown"):
        sfVar = 1
    elif(options.variation=="sfUp"):
        sfVar = 2
    else:
        sfVar = 0 

    #----------Triggers------------#
    beforeTrigCheckpoint    = a.GetActiveNode()
    if(year=="2016" or year=="2016APV"):
        triggerList         = ["HLT_Photon175"]
    elif(year=="2017"):
        triggerList         =["HLT_Photon200"]
    elif(year=="2018"):
       triggerList          =["HLT_Photon200"]
    triggersStringAll       = a.GetTriggerString(triggerList)    
    if(isData): #Only applying trigger to data, will apply trigger turn-on to MC
        a.Cut("Triggers",triggersStringAll)
    nTrig                   = getNweighted(a,isData)
    #------------------------------#

    if("ZGamma" in options.process):
        a.Define("genGammaPt","genGammaPt(nGenPart,GenPart_pdgId,GenPart_pt,GenPart_statusFlags)")


    #----------Selection-----------#
    a.Cut("JetAndPhoton","nFatJet>0 && nPhoton>0")
    a.Define("Hidx","leadingNonGammaAK8Idx(nFatJet,FatJet_eta,FatJet_phi,Photon_eta[0],Photon_phi[0])")
    a.Cut("HidxCut","Hidx>-1")
    nJetGamma = getNweighted(a,isData)
    a.Cut("EtaCut","abs(FatJet_eta[Hidx])<2.4 && abs(Photon_eta[0])<2.4")
    nEta = getNweighted(a,isData)

    a.Cut("GammaID","Photon_cutBased[0]==3")#cut-based ID Fall17V2 Tight
    #a.Cut("GammaID","Photon_mvaID_WP80[0]==1")#MVA ID Fall17V2
    nID = getNweighted(a,isData)

    evtColumns = VarGroup("Event columns")
    evtColumns.Add('FatJet_pt_corr','hardware::MultiHadamardProduct(FatJet_pt,%s)'%ptCorrector)
    evtColumns.Add('FatJet_msoftdrop_corr','hardware::MultiHadamardProduct(FatJet_msoftdrop,%s)'%massCorrector)
    evtColumns.Add('FatJet_mpnet_corr','hardware::MultiHadamardProduct(FatJet_particleNet_mass,%s)'%massCorrector)
    evtColumns.Add("Higgs_pt","FatJet_pt_corr[Hidx]")
    evtColumns.Add("Higgs_eta","FatJet_eta[Hidx]")
    evtColumns.Add("Higgs_phi","FatJet_phi[Hidx]")
    evtColumns.Add("HiggsSDMass",'FatJet_msoftdrop_corr[Hidx]')
    
    #Adding custom jms uncertainty of 2%
    if(options.variation=="jmsUp"):
        evtColumns.Add("HiggsPnetMass",'FatJet_mpnet_corr[Hidx]*1.02')
    elif(options.variation=="jmsDown"):
        evtColumns.Add("HiggsPnetMass",'FatJet_mpnet_corr[Hidx]*0.98')
    else:
        evtColumns.Add("HiggsPnetMass",'FatJet_mpnet_corr[Hidx]')
    evtColumns.Add("Gamma_eta","Photon_eta[0]")
    evtColumns.Add("Gamma_phi","Photon_phi[0]")
    evtColumns.Add("nEle","nElectrons(nElectron,Electron_cutBased,0,Electron_pt,20,Electron_eta)")
    #0:fail,1:veto,2:loose,3:medium,4:tight
    #condition is, cutBased>cut
    evtColumns.Add("nMu","nMuons(nMuon,Muon_looseId,Muon_pfIsoId,0,Muon_pt,20,Muon_eta)")
    #1=PFIsoVeryLoose, 2=PFIsoLoose, 3=PFIsoMedium, 4=PFIsoTight, 5=PFIsoVeryTight, 6=PFIsoVeryVeryTight

    if(isData):
        evtColumns.Add("btagDisc",'Jet_btagDeepB') 
    else:
        CompileCpp('TIMBER/Framework/Zbb_modules/AK4Btag_SF.cc')
        print('AK4Btag_SF ak4SF = AK4Btag_SF("{0}", "DeepCSV", "reshaping");'.format(options.year))
        CompileCpp('AK4Btag_SF ak4SF = AK4Btag_SF("{0}", "DeepCSV", "reshaping");'.format(options.year))
        evtColumns.Add("btagDisc",'ak4SF.evalCollection(nJet,Jet_pt, Jet_eta, Jet_hadronFlavour,Jet_btagDeepB,{0})'.format(sfVar)) 
    evtColumns.Add("topVetoFlag","topVeto(Higgs_eta,Higgs_phi,nJet,Jet_eta,{0},Jet_phi,Jet_pt,btagDisc,{1})".format(2.4,deepCsvM))


    a.Apply([evtColumns])

    #Apply photon energy scale/resolution unc. scale factors
    if("photonEs" in options.variation):
        import correctionlib
        correctionlib.register_pyroot_binding()
        if(year=="2016APV"):
            yearAlt = "2016preVFP"
        elif(year=="2016"):
            yearAlt = "2016postVFP"
        else:
            yearAlt = year
        ROOT.gInterpreter.Declare('auto photon_scale_set = correction::CorrectionSet::from_file("/users/mrogul/Work/Hgamma/Hgamma/data/EGM_ScaleUnc_{0}.json");'.format(year))
        ROOT.gInterpreter.Declare('auto photon_scale_unc = photon_scale_set->at("UL-EGM_ScaleUnc");')
        if(options.variation=="photonEsUp"):
            a.Define("Gamma_gain","Photon_seedGain[0]")
            a.Define('scale_corr','photon_scale_unc->evaluate({"%s","scaleup",Gamma_eta,Gamma_gain})'%yearAlt)
            a.Define("Gamma_pt","Photon_pt[0]*scale_corr")
        elif(options.variation=="photonEsDown"):
            a.Define("Gamma_gain","Photon_seedGain[0]")
            a.Define('scale_corr','photon_scale_unc->evaluate({"%s","scaledown",Gamma_eta,Gamma_gain})'%yearAlt)
            a.Define("Gamma_pt","Photon_pt[0]*scale_corr")
        else:
            print("WARNING: Unknown photon es variation")
    elif(options.variation=="photonErUp"):
        a.Define("Gamma_pt","Photon_pt[0]*(1+Photon_dEsigmaUp[0])")
    elif(options.variation=="photonErDown"):
        a.Define("Gamma_pt","Photon_pt[0]*(1+Photon_dEsigmaDown[0])")
    else:
        a.Define("Gamma_pt","Photon_pt[0]")

    a.Cut("pT","Higgs_pt>300 && Gamma_pt>300")
    npT = getNweighted(a,isData)

    a.Cut("JetPnetMassCut","HiggsPnetMass>50")
    nJetMass = getNweighted(a,isData)

    a.Cut("LeptonVeto","nMu==0 && nEle==0")
    nLeptonVeto = getNweighted(a,isData)


    a.Cut("topVeto","topVetoFlag==0")
    nTopVeto = getNweighted(a,isData)

    a.Define("pnetHiggs","FatJet_particleNetMD_Xbb[Hidx]/(FatJet_particleNetMD_Xbb[Hidx]+FatJet_particleNetMD_QCD[Hidx])")

    checkpoint  = a.GetActiveNode()
    #-----Trigger study part-------#
    if(options.variation=="nom"):
        baselineTrigger="HLT_PFJet260"
        a.SetActiveNode(beforeTrigCheckpoint)
        a.Cut("Baseline",baselineTrigger)
        if(MetFiltersString):
            a.Cut("MET For Trigger",MetFiltersString)
        #need to change names to create nodes with different names than already existing
        a.Cut("JetAndPhotonForTrig","nFatJet>0 && nPhoton>0")
        a.Define("Hidx","leadingNonGammaAK8Idx(nFatJet,FatJet_eta,FatJet_phi,Photon_eta[0],Photon_phi[0])")
        a.Cut("HidxCutForTrig","Hidx>-1")
        nJetGamma = getNweighted(a,isData)
        a.Cut("EtaCutForTrig","abs(FatJet_eta[Hidx])<2.4 && abs(Photon_eta[0])<2.4")
        a.Cut("GammaIDForTrig","Photon_cutBased[0]==3")
        evtColumns.name = "Event Columns For Trigger"
        a.Apply([evtColumns])
        a.Define("Gamma_pt_trig","Photon_pt[0]")
        a.Cut("pT_ForTrigger","Higgs_pt>300 && Gamma_pt_trig>300")
        a.Cut("JetPnetMassCut_ForTrigger","HiggsPnetMass>50")

        triggersStringAll   = a.GetTriggerString(triggerList)  
        h_pTnoTriggers      = a.GetActiveNode().DataFrame.Histo1D(('{0}_GammapTnoTriggers'.format(options.process),';Gamma pT [GeV]; Events/10 GeV;',70,300,1000),"Gamma_pt_trig","genWeight")
        a.Cut("Triggers for trig measurement",triggersStringAll)
        h_pTtriggersAll      = a.GetActiveNode().DataFrame.Histo1D(('{0}_GammapTtriggersAll'.format(options.process),';Gamma pT [GeV]; Events/10 GeV;',70,300,1000),"Gamma_pt_trig","genWeight")

        histos.append(h_pTnoTriggers)
        histos.append(h_pTtriggersAll)
        a.SetActiveNode(checkpoint)
    #------------------------------#



    #--------Store Output----------#

    snapshotColumns = ["pnetHiggs","Higgs_pt","Gamma_pt","Gamma_eta","Gamma_phi","HiggsSDMass","HiggsPnetMass","PV_npvsGood","nFatJet","nPhoton"]
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

    if("ZGamma" in options.process):
        snapshotColumns.append("genGammaPt")


    a.Snapshot(snapshotColumns,outputFile,'Events',saveRunChain=False)

    cutFlowVars         = [nProc,nSkimmed,nTrig,nJetGamma,nEta,nID,npT,nJetMass,nLeptonVeto,nTopVeto]
    cutFlowLabels       = ["Processed","Skimmed","Trigger","JetPlusGamma","Eta","Gamma ID","pT","JetMass","Lepton Veto","Top Veto","pass","fail"]#tagging bins will be filled out in template making
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
