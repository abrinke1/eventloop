#! /usr/bin/env python

## Plots for ttbar --> mu + ele + AK8 jet control region
## https://github.com/abrinke1/cmssw/blob/AWB_displaced_dev_v1_macros/L1Trigger/L1TNtuples/macros/SkimNTuples.py

import os
import sys
import subprocess
import numpy as np
import ROOT as R

R.gROOT.SetBatch(True)  ## Don't display histograms or canvases when drawn

# ## Load local python modules
# sys.path.insert(0, '%s/python' % os.getcwd())
# import EventWeights
# EVT_WGT_BPH = EventWeights.GetWeights()
# import EventWeightsMuEG
# EVT_WGT_MUEG = EventWeightsMuEG.GetWeightsMuEG()
# import EventSelection
# EVT_SEL = EventSelection.GetSelection()

## User configuration
PRT_EVT = 100  ## Print every Nth event
MAX_EVT = 10000     ## Number of events to process
LUMI    = 59830  ## 2018 integrated lumi (pb-1), certified "Good"

CUT_GEN_BB = True ## Require AK8 jet to match GEN bb pair from the ttbar decay
PRT_MVA_IN = False ## Print inputs to MVA, and use looser selection to boost training statistics

LABEL = 'TTTo2L2Nu'
CAT   = 'DoubleLep'

if len(sys.argv) > 1:
    print('\nLABEL changed from %s to %s' % (LABEL, str(sys.argv[1])))
    LABEL = str(sys.argv[1])
if len(sys.argv) > 2:
    print('\nCAT changed from %s to %s' % (CAT, str(sys.argv[2])))
    CAT = str(sys.argv[2])
if len(sys.argv) > 3:
    print('\nCUT_GEN_BB changed from %d to %d' % (CUT_GEN_BB, int(sys.argv[3])))
    CUT_GEN_BB = int(sys.argv[3])


def main():

    print('\nInside ControlTTBarMuEleFat3\n')

    in_file_names = []
    in_dir = '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/'
    if LABEL == 'MuonEG':        in_dir += 'data/MuonEG/Nano25Oct2019/2018/SkimsTTBar/'
    if LABEL == 'TTTo2L2Nu':     in_dir += 'MC/PNet_v1_2023_10_06/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/r1/'
    if 'TTJets_Incl' in LABEL:   in_dir += 'MC/TTJets/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/Nano25Oct2019/SkimsTTBar/'
    if LABEL == 'ZJets_MG':      in_dir += 'MC/ZJets/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/Nano25Oct2019/SkimsTTBar/'
    if LABEL == 'ZJets_AMC':     in_dir += 'MC/ZJets/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/Nano25Oct2019/SkimsTTBar/'
    if LABEL == 'ST_tW_top':     in_dir += 'MC/SingleTop/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/Nano25Oct2019/SkimsTTBar/'
    if LABEL == 'ST_tW_antitop': in_dir += 'MC/SingleTop/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/Nano25Oct2019/SkimsTTBar/'
    if LABEL == 'WW':            in_dir += 'MC/WW/WW_TuneCP5_PSweights_13TeV-pythia8/Nano25Oct2019/SkimsTTBar/'
    if LABEL == 'WZ':            in_dir += 'MC/WZ/WZ_TuneCP5_PSweights_13TeV-pythia8/Nano25Oct2019/SkimsTTBar/'
    if LABEL == 'ZZ':            in_dir += 'MC/ZZ/ZZ_TuneCP5_13TeV-pythia8/Nano25Oct2019/SkimsTTBar/'

    isTTBar = (('TTTo' in LABEL) or ('TTJets' in LABEL)) 

    if not LABEL.split('_')[0] in in_dir:
        print('\n\n***  TRIED TO APPLY LABEL %s TO DIRECTORY %s  ***' % (LABEL, in_dir))
        print('***  CORRECT? QUITTING!  ***')
        sys.exit()
    if not isTTBar and CUT_GEN_BB:
        print('\n\n***  TRIED TO REQUIRE CUT_GEN_BB %s TO DIRECTORY %s  ***' % (CUT_GEN_BB, in_dir))
        print('***  CORRECT? QUITTING!  ***')
        sys.exit()

    for f_name in subprocess.check_output(['ls', in_dir], encoding='UTF-8').splitlines():
        if not '.root' in str(f_name): continue
        # if not 'Skim_Mu_Ele_Fat.root' in f_name: continue
        # if not 'noAK4_bMed.root' in f_name: continue
        if not '2-4_small' in f_name: continue
        in_file_names.append(in_dir+f_name)
        print('Appending file: %s' % in_file_names[-1])

    out_dir = 'plots/'
    out_file_str = out_dir+'ControlTTBarMuEleFat_%s_%s' % (CAT, LABEL)
    # out_file_str = out_dir+'ControlTTBarMuEleFat_MuonEG_2018'
    if CUT_GEN_BB: out_file_str += '_GenBB'
    if CUT_GEN_BB and PRT_MVA_IN: out_file_str += '_Pt140'
    else:          out_file_str += '_noGenBB'
    if PRT_MVA_IN: out_file_str += '_mvaIn'
    if MAX_EVT > 0: out_file_str += '_%dk' % (MAX_EVT / 1000)
    out_file_str += '.root'
    out_file = R.TFile(out_file_str,'recreate')

    isData = False
    if 'MuonEG' in in_dir: isData = True

    chains = {}
    chains['Events'] = 0

    for i in range(len(in_file_names)):
        print('Adding file %s' % in_file_names[i])
        for key in chains.keys():
            if i == 0: chains[key] = R.TChain(key)
            chains[key].Add( in_file_names[i] )
            print('Added TChain %s' % key)


    ## Histogram binning for each variable
    bins = {}
    bins['pt']   = [200,  0, 1000]
    bins['eta']  = [200, -5,    5]
    bins['phi']  = [128, -3.2, 3.2]
    bins['sum']  = [100,  0, 1000]
    bins['dR']   = [100,  0,   10]
    bins['dEta'] = [100,  0,   10]
    bins['dPhi'] = [ 32,  0,  3.2]
    bins['mass'] = [200,  0, 2000]
    bins['MT']   = [200,  0, 2000]
    bins['HT']   = [200,  0, 2000]


    ## Sets of objects to be plotted
    Objs = ['lep1','lep2','fat','MET','ISR','ISRs','ll','llv','llJ','llvJ','jet1','jet2','jet3']
    Vars = ['pt','eta','phi','dR','dEta','dPhi','mass','MT','HT']

    ## Book histograms (most will not be used, and will be deleted at the end)
    hst = {}
    ## All combinations of one and two objects for each variable
    for var in Vars:
        for obj1 in Objs:
            hst['%s_%s' % (var, obj1)] = R.TH1D( 'h_%s_%s' % (var, obj1), '%s %s' % (obj1, var),
                                                bins[var][0], bins[var][1], bins[var][2] )
            for obj2 in Objs:
                hst['%s_%s_%s' % (var, obj1, obj2)] = R.TH1D( 'h_%s_%s_%s' % (var, obj1, obj2),
                                                              '%s(%s, %s)' % (var, obj1, obj2),
                                                              bins[var][0], bins[var][1], bins[var][2] )


    hst['nLep'] = R.TH1D('h_nLep', 'h_nLep',  6, -0.5,  5.5)
    hst['nMu']  = R.TH1D('h_nMu',  'h_nMu',   6, -0.5,  5.5)
    hst['nEle'] = R.TH1D('h_nEle', 'h_nEle',  6, -0.5,  5.5)
    hst['nFat'] = R.TH1D('h_nFat', 'h_nFat',  6, -0.5,  5.5)
    hst['nJet'] = R.TH1D('h_nJet', 'h_nJet', 17, -0.5, 16.5)
    hst['nMuInFat']  = R.TH1D('h_nMuInFat',  'h_nMuInFat',   6, -0.5,  5.5)
    hst['nEleInFat'] = R.TH1D('h_nEleInFat', 'h_nEleInFat',  6, -0.5,  5.5)

    del hst['mass_ll']
    hst['mass_ll'] = R.TH1D('h_mass_ll', 'h_mass_ll', 200, 0, 400)

    hst['mu_mvaTTH']      = R.TH1D('h_mu_mvaTTH',      'h_mu_mvaTTH',      40, -1, 1)
    hst['ele_mvaTTH']     = R.TH1D('h_ele_mvaTTH',     'h_ele_mvaTTH',     40, -1, 1)
    hst['min_lep_mvaTTH'] = R.TH1D('h_min_lep_mvaTTH', 'h_min_lep_mvaTTH', 40, -1, 1)

    hst['mu_miniIso']      = R.TH1D('h_mu_miniIso',      'h_mu_miniIso',      11, 0, 0.22)
    hst['ele_miniIso']     = R.TH1D('h_ele_miniIso',     'h_ele_miniIso',     16, 0, 0.8)
    hst['max_lep_miniIso'] = R.TH1D('h_max_lep_miniIso', 'h_max_lep_miniIso', 16, 0, 0.8)

    hst['mu_ID']  = R.TH1D('h_mu_ID',  'h_mu_ID',   4, -0.5,  3.5)
    hst['ele_ID'] = R.TH1D('h_ele_ID', 'h_ele_ID', 32, -0.5, 31.5)

    hst['mu_SIP']      = R.TH1D('h_mu_SIP',      'h_mu_SIP',      16, 0, 8.0)
    hst['ele_SIP']     = R.TH1D('h_ele_SIP',     'h_ele_SIP',     16, 0, 8.0)
    hst['max_lep_SIP'] = R.TH1D('h_max_lep_SIP', 'h_max_lep_SIP', 16, 0, 8.0)

    hst['mass_min_lep_fat'] = R.TH1D('h_mass_min_lep_fat', 'min mass(lep, fat)', 100, 0, 1000)
    hst['mass_max_lep_fat'] = R.TH1D('h_mass_max_lep_fat', 'max mass(lep, fat)', 100, 0, 1000)
    hst['dR_min_lep_fat']   = R.TH1D('h_dR_min_lep_fat',   'min dR(lep, fat)',   100, 0, 10)
    hst['dR_max_lep_fat']   = R.TH1D('h_dR_max_lep_fat',   'max dR(lep, fat)',   100, 0, 10)
    hst['dPhi_min_lep_fat'] = R.TH1D('h_dPhi_min_lep_fat', 'min dPhi(lep, fat)',  32, 0, 3.2)
    hst['dPhi_max_lep_fat'] = R.TH1D('h_dPhi_max_lep_fat', 'max dPhi(lep, fat)',  32, 0, 3.2)
    hst['dEta_min_lep_fat'] = R.TH1D('h_dEta_min_lep_fat', 'min dEta(lep, fat)', 100, 0, 10)
    hst['dEta_max_lep_fat'] = R.TH1D('h_dEta_max_lep_fat', 'max dEta(lep, fat)', 100, 0, 10)

    hst['MT_min_lep_MET'] = R.TH1D('h_MT_min_lep_MET', 'min MT(lep, MET)', 100, 0, 500)
    hst['MT_max_lep_MET'] = R.TH1D('h_MT_max_lep_MET', 'max MT(lep, MET)', 100, 0, 500)

    hst['msoft_fat'] = R.TH1D('h_msoft_fat', 'fat msoft', 100, 0, 500)

    hst['deepB_max_jet'] = R.TH1D('h_deepB_max_jet', 'Maximum AK4 jet deepB tag score', 55, -0.1, 1.0)
    hst['flavB_max_jet'] = R.TH1D('h_flavB_max_jet', 'Maximum AK4 jet deepFlavB tag score', 55, -0.1, 1.0)
    hst['nSV_max_jet']   = R.TH1D('h_nSV_max_jet',   'Maximum AK4 jet # of secondary vertices', 7, -0.5, 6.5)
    hst['doubB_fat']     = R.TH1D('h_doubB_fat',     'AK8 jet btagDDBvLV2 score', 55, -0.1, 1.0)

    hst['nSV_ISR']       = R.TH1D('h_nSV_ISR',  'ISR jet # of secondary vertices',   7, -0.5,  6.5)
    hst['nSV_ISRs']      = R.TH1D('h_nSV_ISRs', 'ISR jets # of secondary vertices', 13, -0.5, 12.5)

    hst['wgt']      = R.TH1D('h_wgt',      'Event weights (unweighted)', 2000, -10, 10.)
    hst['wgt_wgt']  = R.TH1D('h_wgt_wgt',  'Event weights (weighted)',   2000, -10, 10.)
    hst['nPV']      = R.TH1D('h_nPV',      '# of PVs',      101, -0.5, 100.5)
    hst['nPV_good'] = R.TH1D('h_nPV_good', '# of good PVs', 101, -0.5, 100.5)


    ## 2018 DeepCSV and DeepFlavB cuts: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
    DeepB = {'L': 0.1241, 'M': 0.4184, 'T': 0.7527}
    FlavB = {'L': 0.0494, 'M': 0.2770, 'T': 0.7264}


    ## Loop through events, select, and plot
    nEntries = chains['Events'].GetEntries()
    print('\nEntering loop over %d events\n' % (nEntries))

    # print('\nEvent, isSignal, msoft_fat, mass_fat, deepB_max_jet, flavB_max_jet, pt_llJ, pt_llvJ, MT_llvJ, mass_llJ, mass_llvJ, mass_max_lep_fat, mass_min_lep_fat, dR_max_lep_fat, dR_min_lep_fat, dEta_max_lep_fat, dEta_min_lep_fat, dR_ll_fat, dEta_ll_fat\n')

    print('\nEvent, isSignal, mu_ele, msoft_fat, mass_fat, deepB_max_jet, flavB_max_jet, mass_llvJ, MT_llvJ, pt_llvJ, eta_llvJ, mass_llJ, pt_jet1, pt_ISRs, dR_max_lep_fat, dPhi_max_lep_fat, mass_max_lep_fat, dR_min_lep_fat, dR_ll_fat, dPhi_ll_fat, dR_lep1_fat, dPhi_lep1_fat, dPhi_llJ_MET, dPhi_llvJ_ISRs, pt_llvJ_ISRs, dPhi_llJ_ISRs, MT_llJ_ISR\n')


    ch = chains['Events']  ## Shortcut expression
    ## Count passing events
    nPassPre   = 0
    nPassBtag  = 0
    nPassMass  = 0
    nPassLep   = 0
    nPassMET   = 0
    nPassKin   = 0
    nPassTight = 0
    nPassData  = 0
    nPassTrg   = 0
    nPassGen   = 0
    for iEvt in range(nEntries):

        ch.GetEntry(iEvt)

        if iEvt > MAX_EVT and MAX_EVT > 0: break
        if (iEvt % PRT_EVT) == 0: print('Looking at event #%d / %d' % (iEvt, nEntries))

        if ch.nFatJet   < 1: continue

        if ch.nMuon + ch.nElectron            < 2: continue
        if CAT == 'DoubleMu' and ch.nMuon     < 2: continue
        if CAT == 'DoubleEG' and ch.nElectron < 2: continue
        if CAT == 'MuonEG'   and min(ch.nMuon, ch.nElectron) < 1: continue
        
        ## Require qualifying AK8 jet off the bat
        hasGoodFat = False
        for iFat in range(ch.nFatJet):
            if      ch.FatJet_pt [iFat]  < 170: continue
            if  abs(ch.FatJet_eta[iFat]) > 2.4: continue
            if ch.FatJet_msoftdrop[iFat] < 20: continue
            hasGoodFat = True
            break
        if not hasGoodFat: continue

        # ## Remove data events not in Golden JSON
        # if isData:
        #     pass_JSON = EVT_SEL.PassJSON(ch.run, ch.luminosityBlock)
        #     if not pass_JSON:
        #         continue
        ## Remove events with no good collision vertices
        if ch.PV_npvsGood < 1:
            continue

        muIdxs  = []  ## Indices of selected muon
        eleIdxs = []  ## Indices of selected electron
        fatIdxs = []  ## Indices of selected AK8 jets
        bIdxs   = []  ## Indices of b-quarks from the ttbar decay

        muVecs   = []  ## TLorentzVectors of selected muons
        muVecTs  = []  ## TLorentzVectors of selected muons
        eleVecs  = []  ## TLorentzVectors of selected electrons
        eleVecTs = []  ## TLorentzVectors of selected electrons
        lepVecs  = []  ## TLorentzVectors of two selected leptons
        lepVecTs = []  ## TLorentzVectors of two selected leptons
        fatVecs  = []  ## TLorentzVectors of selected AK8 jets
        fatVecTs = []  ## TLorentzVectors of selected AK8 jets
        bVecs    = []  ## TLorentzVectors of b-quarks from ttbar decay

        ## In ttbar events, look for leptons from W decay
        nGenLep = 0
        lastPdgID = -999
        genLepMoms = []
        for iGen in range(-1 if not isTTBar else ch.nGenPart):
            pdgID = ch.GenPart_pdgId[iGen]
            ## Particle must be a charged lepton
            if not abs(pdgID) in [11,13,15]: continue
            ## Particle must have W boson mother
            iMom = ch.GenPart_genPartIdxMother[iGen]
            if iMom < 0: continue
            if abs(ch.GenPart_pdgId[iMom]) != 24: continue
            if iMom in genLepMoms:
                ## Special case of W radiating a gamma --> l+l- pair
                if iMom == genLepMoms[-1] and pdgID == -1*lastPdgID:
                    nGenLep -= 1  ## Delete previous lepton, and skip this one
                continue
            else: genLepMoms.append(iMom)
            nGenLep += 1
            lastPdgID = pdgID
        ## End loop: for iGen in range(-1 if not isTTBar else ch.nGenPart)
        if nGenLep > 2:
            print('\n*** Super-weird event!!! LS = %d, event = %d has %d leptons from W decay! ***\n' % (ch.luminosityBlock, ch.event, nGenLep))
        if 'TTTo2L2Nu' in LABEL and nGenLep < 2:
            print('\n*** Super-weird event!!! LS = %d, event = %d has %d leptons from W decay! ***\n' % (ch.luminosityBlock, ch.event, nGenLep))
            sys.exit()
        ## For inclusive ttbar sample, separate into di-lepton and non di-lepton (hadronic) samples
        if 'TTJets_Incl_DiLept' in LABEL or 'TTTo2L2Nu' in LABEL and nGenLep < 2:
            continue
        if 'TTJets_Incl_Had' in LABEL and nGenLep >= 2:
            continue


        ## Look for final-state b-quarks from ttbar decay
        for iGen in range(-1 if isData else ch.nGenPart):
            ## Particle must be a b-quark
            if abs(ch.GenPart_pdgId[iGen]) != 5: continue

            ## Particle's mom must be the final-state top quark or initial-state b-quark
            iMom = ch.GenPart_genPartIdxMother[iGen]
            if iMom < 0: continue
            if not ( (abs(ch.GenPart_pdgId[iMom]) == 6 and ch.GenPart_status[iMom] == 62) or \
                     (abs(ch.GenPart_pdgId[iMom]) == 5 and ch.GenPart_status[iMom] == 23) ): continue

            ## If particle's mom is a b, it must come directly from top decay
            if abs(ch.GenPart_pdgId[iMom]) == 5:
                iGrand = ch.GenPart_genPartIdxMother[iMom]
                if iGrand < 0: continue
                if (abs(ch.GenPart_pdgId[iGrand]) != 6 or abs(ch.GenPart_status[iGrand]) != 62): continue

            ## If nearly identical to previous GEN particle, skip
            if len(bIdxs) > 0 and \
               ch.GenPart_pdgId[iGen]  == ch.GenPart_pdgId[bIdxs[-1]]  and \
               ch.GenPart_status[iGen] == ch.GenPart_status[bIdxs[-1]] and \
               abs((ch.GenPart_pt[iGen] / bVecs[-1].Pt()) - 1.0) < 0.1 and \
               abs(ch.GenPart_eta[iGen] - bVecs[-1].Eta())       < 0.1: continue

            ## Set the 4-vector for this GEN b
            bVec = R.TLorentzVector()
            bVec.SetPtEtaPhiM( ch.GenPart_pt  [iGen],
                               ch.GenPart_eta [iGen],
                               ch.GenPart_phi [iGen],
                               ch.GenPart_mass[iGen] )

            ## If particle's mom is a b, and the mom exists in the current set of bVecs, replace it
            if abs(ch.GenPart_pdgId[iMom]) == 5:
                for iIdx in range(len(bIdxs)):
                    if ch.GenPart_genPartIdxMother[iMom] == bIdxs[iIdx]:
                        print('\nIn LS %d, event %d, replacing idx %d with %d' % (ch.luminosityBlock, ch.event, bIdxs[iIdx], iGen))
                        bIdxs[iIdx] = iGen
                        bVecs[iIdx] = bVec
            else:
                bIdxs.append(iGen)
                bVecs.append(bVec)

        ## End loop: for iGen in range(nGen)


        ## Check that event *has* at least 2 GEN b-quarks
        if CUT_GEN_BB and len(bVecs) != 2:
            print('\n\n*** Interesting LS %d, event %d with %d GEN b-quarks' % (ch.luminosityBlock, ch.event, len(bVecs)))
            if len(bVecs) > 2:
                for bVec in bVecs:
                    print('  * pT = %.2f, eta = %.2f, phi = %.2f' % (bVec.Pt(), bVec.Eta(), bVec.Phi()))
                    for iGen in range(ch.nGenPart):
                        if ( abs(ch.GenPart_pdgId[iGen]) == 5 and
                             abs(ch.GenPart_pt[iGen]  - bVec.Pt()) < 0.1 and
                             abs(ch.GenPart_eta[iGen] - bVec.Eta()) < 0.1 ):
                            print('    - ID = %d, Idx = %d, MomIdx = %d, Status = %d' % (ch.GenPart_pdgId[iGen], iGen, ch.GenPart_genPartIdxMother[iGen], ch.GenPart_status[iGen]))
                print('')
            if len(bVecs) < 2:
                print('\n\nZERO OR ONE???! WHAT KIND OF TTBAR EVENT IS THIS?!!!')
                for iGen in range(ch.nGenPart):
                    if abs(ch.GenPart_pdgId[iGen]) == 5:
                        print('  * ID = %d, Idx = %d, MomIdx = %d, Status = %d' % (ch.GenPart_pdgId[iGen], iGen, ch.GenPart_genPartIdxMother[iGen], ch.GenPart_status[iGen]))
                print('')
        ## End conditional: if CUT_GEN_BB and len(bVecs) != 2

        ## Skip "signal" events if GEN b-jets are too far away from each other, or have too-low pT
        if CUT_GEN_BB and bVecs[0].DeltaR(bVecs[1]) > 1.6: continue
        if CUT_GEN_BB and PRT_MVA_IN and (bVecs[0]+bVecs[1]).Pt() < 140: continue

        ## Find muon(s) passing cuts
        for iMu in range(ch.nMuon):
            if              ch.Muon_pt [iMu]  < 10.: continue
            if          abs(ch.Muon_eta[iMu]) > 2.4: continue
            if   ch.Muon_mediumPromptId[iMu]  != 1 : continue
            ## if     ch.Muon_miniIsoId[iMu]  <  2 : continue  ## Doesn't work!!! - AWB 2021.06.10
            ## Equivalent: https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatAlgos/plugins/PATMuonProducer.cc#L699
            if ch.Muon_miniPFRelIso_all[iMu] > 0.20: continue
            if           ch.Muon_mvaTTH[iMu] < -0.4: continue

            ## Save the selected muon
            muVec  = R.TLorentzVector()
            muVecT = R.TLorentzVector()
            muVec .SetPtEtaPhiM( ch.Muon_pt[iMu], ch.Muon_eta[iMu], ch.Muon_phi[iMu], 0.106 )
            muVecT.SetPtEtaPhiM( ch.Muon_pt[iMu],                0, ch.Muon_phi[iMu], 0.106 )
            muVecs .append(muVec)
            muVecTs.append(muVecT)
            muIdxs .append(iMu)
        ## End loop: for iMu in range(ch.nMuon)


        ## Find electron(s) passing cuts
        for iEle in range(ch.nElectron):
            if                 ch.Electron_pt [iEle]  < 10.: continue
            if             abs(ch.Electron_eta[iEle]) > 2.5: continue
            ## https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2#Recommended_MVA_Recipe_V2_for_re
            if ch.Electron_mvaFall17V2Iso_WP90[iEle]  != 1 : continue
            if              ch.Electron_mvaTTH[iEle] < -0.4: continue

            ## Save the selected electron
            eleVec  = R.TLorentzVector()
            eleVecT = R.TLorentzVector()
            eleVec .SetPtEtaPhiM( ch.Electron_pt[iEle], ch.Electron_eta[iEle], ch.Electron_phi[iEle], 0.0005 )
            eleVecT.SetPtEtaPhiM( ch.Electron_pt[iEle],                     0, ch.Electron_phi[iEle], 0.0005 )
            eleVecs .append(eleVec)
            eleVecTs.append(eleVecT)
            eleIdxs .append(iEle)
        ## End loop: for iEle in range(ch.nElectron)


        ## Only keep events with exactly two selected leptons of opposite charge
        if len(muIdxs) + len(eleIdxs) != 2: continue

        ## Define a new lepton category event-by-event
        LepCat = CAT

        if CAT == 'DoubleLep':
            if min(len(muIdxs), len(eleIdxs)) > 0 and     ch.Muon_charge[ muIdxs[0]] + ch.Electron_charge[eleIdxs[0]] == 0: LepCat = 'MuonEG'
            elif   len(muIdxs)                > 1 and     ch.Muon_charge[ muIdxs[0]] +     ch.Muon_charge[ muIdxs[1]] == 0: LepCat = 'DoubleMu'
            elif                len(eleIdxs)  > 1 and ch.Electron_charge[eleIdxs[0]] + ch.Electron_charge[eleIdxs[1]] == 0: LepCat = 'DoubleEG'
            else: continue
        if LepCat == 'DoubleMu' and (    len(muIdxs)  < 2 or               ch.Muon_charge[ muIdxs[0]] +     ch.Muon_charge[ muIdxs[1]] != 0): continue
        if LepCat == 'DoubleEG' and (    len(eleIdxs) < 2 or           ch.Electron_charge[eleIdxs[0]] + ch.Electron_charge[eleIdxs[1]] != 0): continue
        if LepCat == 'MuonEG'   and (min(len(muIdxs), len(eleIdxs)) < 1 or ch.Muon_charge[ muIdxs[0]] + ch.Electron_charge[eleIdxs[0]] != 0): continue

        if LepCat == 'DoubleMu':
            xMu1 = muIdxs[0]
            xMu2 = muIdxs[1]
        if LepCat == 'DoubleEG':
            xEle1 = eleIdxs[0]
            xEle2 = eleIdxs[1]
        if LepCat == 'MuonEG':
            xMu  = muIdxs[0]
            xEle = eleIdxs[0]

        ## Set lepton vectors by high-pT, low-pT
        if LepCat == 'DoubleMu':
            lepVecsTmp  = [muVecs [0], muVecs [1]]
            lepVecTsTmp = [muVecTs[0], muVecTs[1]]
        if LepCat == 'DoubleEG':
            lepVecsTmp  = [eleVecs [0], eleVecs [1]]
            lepVecTsTmp = [eleVecTs[0], eleVecTs[1]]
        if LepCat == 'MuonEG':
            lepVecsTmp  = [muVecs [0], eleVecs [0]]
            lepVecTsTmp = [muVecTs[0], eleVecTs[0]]

        if lepVecsTmp[0].Pt() < lepVecsTmp[1].Pt():
            lepVecs  = [lepVecsTmp [1], lepVecsTmp [0]]
            lepVecTs = [lepVecTsTmp[1], lepVecTsTmp[0]]
        else:
            lepVecs  = lepVecsTmp
            lepVecTs = lepVecTsTmp
        if lepVecs[0].Pt() < lepVecs[1].Pt():
            print('\n*** Super-weird event!!! LS = %d, event = %d has pT mis-ordered (%f, %f). ***\n' % (ch.luminosityBlock, ch.event, lepVecs[0].Pt(), lepVecs[1].Pt()))


        ## Find AK8 jet(s) passing cuts
        for iFat in range(ch.nFatJet):
            if      ch.FatJet_pt [iFat]  < 170: continue
            if  abs(ch.FatJet_eta[iFat]) > 2.4: continue
            if ch.FatJet_msoftdrop[iFat] < 20:  continue
            ## Save 4-vectors of passing AK8 jets
            fatVec  = R.TLorentzVector()
            fatVecT = R.TLorentzVector()
            fatVec .SetPtEtaPhiM( ch.FatJet_pt[iFat], ch.FatJet_eta[iFat], ch.FatJet_phi[iFat], ch.FatJet_mass[iFat] )
            fatVecT.SetPtEtaPhiM( ch.FatJet_pt[iFat],                   0, ch.FatJet_phi[iFat], ch.FatJet_mass[iFat] )
            ## Skip jets which overlap selected leptons
            if fatVec.DeltaR(lepVecs[0]) < 0.8 or fatVec.DeltaR(lepVecs[1]) < 0.8: continue
            ## Save selected AK8 jet
            fatIdxs .append(iFat)
            fatVecs .append(fatVec)
            fatVecTs.append(fatVecT)
        ## End loop: for iFat in range(ch.nFatJet)

        ## Only keep event if there is at least one selected fat jet
        if len(fatVecs) == 0: continue


        ## Store secondary vertex 4-vectors
        svVecs = []
        for iSV in range(ch.nSV):
            svVec = R.TLorentzVector()
            svVec.SetPtEtaPhiM( ch.SV_pt[iSV], ch.SV_eta[iSV], ch.SV_phi[iSV], ch.SV_mass[iSV])
            svVecs.append(svVec)


        ##############################################
        ###  Big loop over all fat jet candidates  ###
        ##############################################
        hasCand = False  ## Track whether we've already found an AK8 jet passing all selection cuts
        for jFat in range(len(fatVecs)):
            if hasCand: break

            ## Select a single AK8 jet candidate
            fatIdx  = fatIdxs[jFat]
            fatVec  = fatVecs[jFat]
            fatVecT = fatVecTs[jFat]

            ## Count selected leptons overlapping selected AK8 jet
            muInFat,eleInFat = (0,0)
            for muVec in muVecs:
                if fatVec.DeltaR(muVec) < 0.8: muInFat += 1
            for eleVec in eleVecs:
                if fatVec.DeltaR(eleVec) < 0.8: eleInFat += 1

            ## Shortcut 4-vectors for multi-object quantities
            ll_vec    = lepVecs [0] + lepVecs [1]
            ll_vecT   = lepVecTs[0] + lepVecTs[1]
            llJ_vec   = ll_vec  + fatVec
            llJ_vecT  = ll_vecT + fatVecT

            ## Store MET vector (eta depends on selected AK8 jet)
            metVec  = R.TLorentzVector()
            metVecT = R.TLorentzVector()
            metVec .SetPtEtaPhiM( ch.MET_pt, llJ_vec.Eta(), ch.MET_phi, 0)
            metVecT.SetPtEtaPhiM( ch.MET_pt,             0, ch.MET_phi, 0)

            ## Shortcut 4-vector for entire ttbar system
            llv_vec   = ll_vec   + metVec
            llv_vecT  = ll_vecT  + metVecT
            llvJ_vec  = llJ_vec  + metVec
            llvJ_vecT = llJ_vecT + metVecT


            jetIdxs  = []  ## Indices of selected AK4 jets
            jetVecs  = []  ## TLorentzVectors of selected AK4 jets
            jetVecTs = []  ## Transverse component only (eta = 0)
            isrVec   = R.TLorentzVector()  ## First selected AK4 jet
            isrVecT  = R.TLorentzVector()
            isrsVec  = R.TLorentzVector()  ## All selected AK4 jets
            isrsVecT = R.TLorentzVector()
            isrs_HT  = 0   ## Sum of pT of selected AK4 jets

            nSV_max_jet = 0  ## Maximum number of secondary vertices in one AK4 jet
            nSV_ISR     = 0  ## Number of secondary vertices in highest-pT AK4 jet
            nSV_ISRs    = 0  ## Total number of secondary vertices in all AK4 jets

            ## Store AK4 jets not matched to AK8 jet
            max_deepB = -9.99
            max_flavB = -9.99
            for iJet in range(ch.nJet):
                if  abs(ch.Jet_eta[iJet]) > 4.7: continue
                if      ch.Jet_pt [iJet]  <  25: continue
                ## Tight jet ID: https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2018
                if    ch.Jet_jetId[iJet]  <=  1: continue
                ## Loose pileup ID: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID
                if      ch.Jet_pt [iJet]  <  50:
                    if ch.Jet_puId[iJet]  <=  3: continue
                ## Save 4-vectors of loose b-tagged AK4 jets
                jetVec  = R.TLorentzVector()
                jetVecT = R.TLorentzVector()
                jetVec .SetPtEtaPhiM( ch.Jet_pt[iJet], ch.Jet_eta[iJet], ch.Jet_phi[iJet], ch.Jet_mass[iJet] )
                jetVecT.SetPtEtaPhiM( ch.Jet_pt[iJet],                0, ch.Jet_phi[iJet], ch.Jet_mass[iJet] )
                ## Don't allow overlap with muon, electron, or AK8 jet
                if jetVec.DeltaR(lepVecs[0]) < 0.4 or jetVec.DeltaR(lepVecs[1]) < 0.4 or jetVec.DeltaR(fatVec) < 0.8: continue
                jetIdxs .append(iJet)
                jetVecs .append(jetVec)
                jetVecTs.append(jetVecT)

                ## Store un-matched jets as "ISR" (for "signal" events)
                isrsVec  = isrsVec +jetVec
                isrsVecT = isrsVecT+jetVecT
                isrs_HT += jetVec.Pt()
                if len(jetVecs) == 1:
                    isrVec  = jetVec
                    isrVecT = jetVecT

                ## Store b-tag values of un-matched jets
                if abs(jetVec.Eta()) < 2.4:
                    max_deepB = max(max_deepB, ch.Jet_btagDeepB[iJet])
                    max_flavB = max(max_flavB, ch.Jet_btagDeepFlavB[iJet])

                ## Count matched secondary vertices
                nSV_jet = 0
                for svVec in svVecs:
                    if jetVec.DeltaR(svVec) < 0.4:
                        nSV_jet += 1
                
                nSV_max_jet = max(nSV_max_jet, nSV_jet)
                nSV_ISRs   += nSV_jet
                if len(jetVecs) == 1:
                    nSV_ISR = nSV_jet

            ## End loop: for iJet in range(ch.nJet)


            ##########################
            ## Apply selection cuts ##
            ##########################
            nPassPre += 1
            ## No medium b-tagged AK4 jet
            if max_deepB > DeepB['M'] or max_flavB > FlavB['M']: continue
            nPassBtag += 1
            ## Require soft-drop mass > 20 GeV (duplicates selection above)
            if ch.FatJet_msoftdrop[fatIdx] < 20: continue
            nPassMass += 1
            ## Lepton pT cuts to match double lepton triggers, plus invariant mass cut
            if LepCat != 'MuonEG':
                if ll_vec.M() < 4 or abs(ll_vec.M() - 91.2) < 10: continue
            if not PRT_MVA_IN:
                if ll_vec.M() < 12: continue
                if lepVecs[0].Pt() < 25 or lepVecs[1].Pt() < 15: continue

            # print('\n* Muon pT = %.1f, eta = %.1f, phi = %.1f, miniPFRelIso = %.3f, pfRelIso03 (04) = %.3f (%.3f), dxy = %.6f, dz = %.6f, IP = %.6f, SIP = %.2f, mvaTTH = %.3f' % (muVecs[0].Pt(), muVecs[0].Eta(), muVecs[0].Phi(), ch.Muon_miniPFRelIso_all[xMu], ch.Muon_pfRelIso03_all[xMu], ch.Muon_pfRelIso04_all[xMu], ch.Muon_dxy[xMu], ch.Muon_dz[xMu], ch.Muon_ip3d[xMu], ch.Muon_sip3d[xMu], ch.Muon_mvaTTH[xMu]))
            # print('* Electron pT = %.1f, eta = %.1f, phi = %.1f, miniPFRelIso = %.2f, pfRelIso03 = %.3f, dxy = %.6f, dz = %.6f, IP = %.6f, SIP = %.2f, mvaTTH = %.3f' % (eleVecs[0].Pt(), eleVecs[0].Eta(), eleVecs[0].Phi(), ch.Electron_miniPFRelIso_all[xEle], ch.Electron_pfRelIso03_all[xEle], ch.Electron_dxy[xEle], ch.Electron_dz[xEle], ch.Electron_ip3d[xEle], ch.Electron_sip3d[xEle], ch.Electron_mvaTTH[xEle]))

            nPassLep += 1
            # ## MET cut to suppress background
            # if metVec.Pt() < 30: continue
            nPassMET += 1

            # ## MT(llvJ) < 550 GeV, mass(llvJ) < 650 GeV
            # if llvJ_vecT.M() > 550: continue
            # if llvJ_vec.M() > 650: continue
            nPassKin += 1

            # ## Tighter selection cuts
            # if isrVec .Pt() < 60: continue
            # if isrsVec.Pt() < 60: continue
            # if llvJ_vecT.M() > 450: continue
            # if llvJ_vec .M() > 500: continue
            nPassTight += 1

            # ## Selection cuts to remove background in data
            # if lepVecs[0].Pt() < 20: continue
            # if lepVecs[1].Pt() < 20: continue
            # if metVec.Pt() < 30: continue
            # if llv_vec.M() < 70: continue
            # if min( (lepVecTs[0]+metVecT).M(), (lepVecTs[1]+metVecT).M() ) < 10: continue
            if isData and ch.run > 319077:  ## HEM veto
                if fatVec.Eta() < -1.17 and fatVec.Phi() > -1.97 and fatVec.Phi() < -0.47:
                    continue
            nPassData += 1

            ## Trigger selection in data and MC
            if LepCat == 'DoubleMu' and not PRT_MVA_IN:
                if not ( ch.L1_DoubleMu_15_7 and ch.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 ): continue

            if LepCat == 'DoubleEG' and not PRT_MVA_IN:
                L1_SingleEG = (ch.L1_SingleEG36er2p5 or ch.L1_SingleEG38er2p5 or ch.L1_SingleEG40er2p5 or
                               ch.L1_SingleIsoEG30er2p1 or ch.L1_SingleIsoEG32er2p5) ##  or ch.L1_SingleLooseIsoEG28er2p5)  Not in MC?!? -- AWB 2021.06.11
                L1_DoubleEG = (ch.L1_DoubleEG_25_12_er2p5 or ch.L1_DoubleEG_25_14_er2p5 or
                               ch.L1_DoubleEG_LooseIso22_12_er2p5)
                if not ( (  L1_SingleEG                 and ch.HLT_Ele32_WPTight_Gsf ) or
                         ( (L1_SingleEG or L1_DoubleEG) and ch.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL ) ): continue

            if LepCat == 'MuonEG' and not PRT_MVA_IN:
                ## L1_Mu7_EG20er2p5 was unprescaled for most of the year, but didn't exist in early runs, so leave it out for now
                if not ( ( (ch.L1_Mu7_EG23er2p5 or ch.L1_Mu7_LooseIsoEG20er2p5) and ch.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ ) or
                         ( (ch.L1_Mu20_EG10er2p5 or ch.L1_SingleMu22)           and ch.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL )   or
                         (  ch.L1_SingleMu25                                    and ch.HLT_Mu27_Ele37_CaloIdL_MW )                         or
                         (  ch.L1_SingleMu25                                    and ch.HLT_Mu37_Ele27_CaloIdL_MW ) ):
                    continue
            nPassTrg += 1

            ## AK8 jet has passsed all event selection cuts
            hasCand = True
            ## Signal events must have AK8 jet matching two b-quarks; background events must not
            if CUT_GEN_BB and (fatVec.DeltaR(bVecs[0]) > 0.8 or fatVec.DeltaR(bVecs[1]) > 0.8): continue
            if CUT_GEN_BB and PRT_MVA_IN and (bVecs[0]+bVecs[1]).Pt() < 140: continue
            if isTTBar:
                if (not CUT_GEN_BB) and (fatVec.DeltaR(bVecs[0]) < 0.8 and fatVec.DeltaR(bVecs[1]) < 0.8): continue
            nPassGen += 1

            ## *** GEN pT > 140 GeV --> + mass_llvJ < 650, MT_llvJ < 550 --> + mass/MT < 500/450 + ISR > 60 *** ##
            ## dR_max_lep_fat  , SigSq2 = 3.075 (2.987) --> 2.176 (2.123)  [# 2] --> 1.489 (1.460)  [# 1]  {A+}
            ## mass_llvJ       , SigSq2 = 2.972 (2.965) --> 2.110 (2.096)  [# 1] --> 1.370 (1.362)  [# 3]  {A+}
            ## MT_llvJ         , SigSq2 = 2.813 (2.804) --> 1.886 (1.881)  [# 5] --> 1.238 (1.237)  [#14]  {A+}
            ## dR_ll_fat       , SigSq2 = 2.635 (2.572) --> 1.983 (1.938)  [# 4] --> 1.423 (1.394)  [# 2]  {A}
            ## dPhi_max_lep_fat, SigSq2 = 2.366 (2.256) --> 1.818 (1.742)  [# 6] --> 1.335 (1.264)  [# 5]  {A}
            ## pt_ISRs         , SigSq2 = 2.129 (1.997) --> 1.940 (1.746)  [# 3] --> 1.364 (1.211)  [# 4]  {A}
            ## dR_lep1_fat     , SigSq2 = 2.114 (2.114) --> 1.686 (1.686)  [# 9] --> 1.329 (1.329)  [# 6]  {A}
            ## mass_max_lep_fat, SigSq2 = 2.192 (2.186) --> 1.533 (1.529)  [#13] --> 1.243 (1.241)  [#13]  {C}
            ## mass_llJ        , SigSq2 = 2.094 (2.088) --> 1.478 (1.474)  [#15] --> 1.182 (1.178)  [#23]  {E}
            ## dPhi_ll_fat     , SigSq2 = 2.063 (2.063) --> 1.679 (1.679)  [#10] --> 1.293 (1.293)  [# 8]  {B}
            ## pt_jet1         , SigSq2 = 1.850 (1.730) --> 1.713 (1.579)  [# 8] --> 1.293 (1.150)  [# 9]  {A+}
            ## dPhi_lep1_fat   , SigSq2 = 1.807 (1.807) --> 1.529 (1.529)  [#12] --> 1.251 (1.251)  [#11]  {C}
            ## pt_llvJ         , SigSq2 = 1.776 (1.662) --> 1.706 (1.490)  [# 7] --> 1.286 (1.106)  [#10]  {B}
            ## dR_min_lep_fat  , SigSq2 = 1.685 (1.685) --> 1.472 (1.472)  [#16] --> 1.213 (1.213)  [#16]  {D}
            ## dPhi_llvJ_ISRs  , SigSq2 = 1.676 (1.676) --> 1.564 (1.564)  [#11] --> 1.305 (1.305)  [# 7]  {B}
            ## dPhi_llJ_MET    , SigSq2 = 1.661 (1.661) --> 1.402 (1.402)  [#20] --> 1.146 (1.146)         {E}
            ## dPhi_llJ_ISRs   , SigSq2 = 1.659 (1.659) --> 1.447 (1.447)  [#17] --> 1.180 (1.179)  [#20]  {E}
            ## HT_ISRs         , SigSq2 = 1.545 (1.488) --> 1.515 (1.411)  [#14] --> 1.143 (1.116)         {E}
            ## MT_llJ_ISR      , SigSq2 = 1.078 (1.078) --> 1.403 (1.288)  [#19] --> 1.238 (1.126)  [#15]  {C}
            ## HT_llJ_ISR      , SigSq2 = 1.025 (1.025) --> 1.318 (1.225)        --> 1.194 (1.078)  [#19]  {D}
            ## MT_llvJ_ISR     , SigSq2 = 1.012 (1.012) --> 1.236 (1.159)        --> 1.199 (1.050)  [#18]  {D}
            ## pt_llJ          , SigSq2 = 1.180 (1.177) --> 1.343 (1.279)        --> 1.180 (1.102)  [#21]  {E}
            ## dPhi_llvJ_ISR   , SigSq2 = 1.416 (1.416) --> 1.349 (1.349)        --> 1.176 (1.175)  [#24]  {E}
            ## pt_llvJ_ISRs    , SigSq2 = 1.238 (1.238) --> 1.219 (1.218)        --> 1.173 (1.172)  [#25]  {E}
            ## HT_llvJ_ISR     , SigSq2 = 1.013 (1.011) --> 1.215 (1.189)        --> 1.180 (1.070)  [#22]  {E}
            ## mass_llJ_ISR    , SigSq2 = 1.032 (1.032) --> 1.189 (1.188)        --> 1.168 (1.059)  [#26]  {E}


            ## *** Store highly discriminating output variables for MVA algorithm ***
            ## Event, isSignal, msoft_fat, mass_fat, deepB_max_jet, flavB_max_jet, pt_llJ, pt_llvJ, MT_llvJ, mass_llJ, mass_llvJ,
            ## mass_max_lep_fat, mass_min_lep_fat, dR_max_lep_fat, dR_min_lep_fat, dEta_max_lep_fat, dEta_min_lep_fat, dR_ll_fat, dEta_ll_fat'

            ## Event, isSignal, mu_ele, msoft_fat, mass_fat, deepB_max_jet, flavB_max_jet, mass_llvJ, MT_llvJ, pt_llvJ, eta_llvJ, mass_llJ, pt_jet1, pt_ISRs, dR_max_lep_fat, dR_min_lep_fat,
            ## dPhi_max_lep_fat, mass_max_lep_fat, dR_ll_fat, dPhi_ll_fat, dR_lep1_fat, dPhi_lep1_fat, dPhi_llJ_MET, dPhi_llvJ_ISRs, pt_llvJ_ISRs, dPhi_llJ_ISRs, MT_llJ_ISR'

            if PRT_MVA_IN:
                print('%16d, %d, %d, %6.1f, %6.1f, %6.3f, %6.3f, %6.1f, %6.1f, %6.1f, %5.3f, %6.1f, %6.1f, %6.1f, %5.3f, %5.3f,' \
                    ' %5.3f, %6.1f, %5.3f, %5.3f, %5.3f, %5.3f, %5.3f, %5.3f, %6.1f, %5.3f, %6.1f' \
                    % ( ch.luminosityBlock*1000000000 + ch.event, CUT_GEN_BB, (LepCat == 'MuonEG'), ## Event, isSignal, mu_ele
                        ch.FatJet_msoftdrop[fatIdx], fatVec.M(),                                    ## msoft_fat, mass_fat
                        max(-0.099, max_deepB), max(-0.099, max_flavB),                             ## deepB_max_jet, flavB_max_jet
                        llvJ_vec.M(), llvJ_vecT.M(), llvJ_vec.Pt(), abs(llvJ_vec.Eta()),            ## mass_llvJ, MT_llvJ, pt_llvJ, eta_llvJ
                        llJ_vec.M(), isrVec.Pt(), isrsVec.Pt(),                                     ## mass_llJ, pt_jet1, pt_ISRs
                        max(lepVecs[0].DeltaR(fatVec), lepVecs[1].DeltaR(fatVec)),                  ## dR_max_lep_fat
                        min(lepVecs[0].DeltaR(fatVec), lepVecs[1].DeltaR(fatVec)),                  ## dR_min_lep_fat
                        max(abs(lepVecs[0].DeltaPhi(fatVec)), abs(lepVecs[1].DeltaPhi(fatVec))),    ## dPhi_max_lep_fat
                        max((lepVecs[0]+fatVec).M(), (lepVecs[1]+fatVec).M()),                      ## mass_max_lep_fat
                        ll_vec.DeltaR(fatVec), abs(ll_vec.DeltaPhi(fatVec)),                        ## dR_ll_fat, dPhi_ll_fat
                        lepVecs[0].DeltaR(fatVec), abs(lepVecs[0].DeltaPhi(fatVec)),                ## dR_lep1_fat, dPhi_lep1_fat
                        abs(llJ_vec.DeltaPhi(metVec)),                                              ## dPhi_llJ_MET
                        abs(llvJ_vec.DeltaPhi(isrsVec)), (llvJ_vec+isrsVec).Pt(),                   ## dPhi_llvJ_ISRs, pt_llvJ_ISRs
                        abs(llJ_vec.DeltaPhi(isrsVec)), (llJ_vecT+isrVecT).M() ) )                  ## dPhi_llJ_ISRs, MT_llJ_ISR


            #######################
            ## Get event weights ##
            #######################
            
            WGT         = 1.0  ## Overall event weight
            WGT_NO_LUMI = 1.0  ## Event weight before luminosity scaling

            if not isData:
                # if LepCat == 'MuonEG' and not PRT_MVA_IN:
                #     WGT *= EVT_WGT_MUEG.GetTrigSF   (  muVecs[0].Pt(), eleVecs[0].Pt() )
                #     WGT *= EVT_WGT_MUEG.GetMuonIDSF (  muVecs[0].Pt(),  muVecs[0].Eta() )
                #     WGT *= EVT_WGT_MUEG.GetMuonIsoSF(  muVecs[0].Pt(),  muVecs[0].Eta() )
                #     WGT *= EVT_WGT_MUEG.GetEleIDSF  ( eleVecs[0].Pt(), eleVecs[0].Eta() )

                # # WGT *= EVT_WGT_MUEG.GetBtagSF   ( bFlavs[0], bVecs[0].Eta(), bVecs[0].Pt() )
                # # WGT *= EVT_WGT_MUEG.GetBtagSF   ( bFlavs[1], bVecs[1].Eta(), bVecs[1].Pt() )

                # ## Top pT weighting
                # if isTTBar:
                #     WGT *= EVT_WGT_MUEG.GetTopPtWgt( top_pt[0], top_pt[1] )
                ## HEM veto weight
                if fatVec.Eta() < -1.17 and fatVec.Phi() > -1.97 and fatVec.Phi() < -0.47:
                    WGT *= (21090. / LUMI)
                ## GEN negative weight
                if ch.genWeight < 0:
                    WGT *= -1
            ## End conditional: if not isData

            WGT_NO_LUMI = WGT  ## Track event weight before cross-section and luminosity scaling
            if not isData:
                # WGT *= EVT_WGT_MUEG.GetXsecPerEvt( in_dir )
                WGT *= LUMI
            if MAX_EVT > 0 and MAX_EVT < nEntries:
                WGT *= (1.0*nEntries / MAX_EVT)


            #####################
            ## Fill histograms ##
            #####################

            ## Objs = ['lep1','lep2','fat','MET','ISR','ll','llv','llJ','llvJ','jet1','jet2','jet3']
            ## Vars = ['pt','eta','dR','dEta','dPhi','mass','MT','HT']

            ## Whole event variables
            hst['nPV']     .Fill( min( ch.PV_npvs, 100 ), WGT )
            hst['nPV_good'].Fill( min( ch.PV_npvsGood, 100), WGT )
            hst['wgt']     .Fill( max( 0.001, min( 1.999, WGT_NO_LUMI ) ) )
            hst['wgt_wgt'] .Fill( max( 0.001, min( 1.999, WGT_NO_LUMI ) ), WGT_NO_LUMI )

            hst['nLep'].Fill( len(muVecs) + len(eleVecs) - muInFat - eleInFat, WGT)
            hst['nMu'] .Fill( len(muVecs)  -  muInFat, WGT)
            hst['nEle'].Fill( len(eleVecs) - eleInFat, WGT)
            hst['nFat'].Fill( len(fatVecs), WGT)
            hst['nJet'].Fill( len(jetVecs), WGT)
            hst['nMuInFat'] .Fill(  muInFat, WGT)
            hst['nEleInFat'].Fill( eleInFat, WGT)

            ## Lepton ID
            if LepCat == 'MuonEG':
                hst['mu_mvaTTH']     .Fill( ch.Muon_mvaTTH[xMu], WGT)
                hst['ele_mvaTTH']    .Fill( ch.Electron_mvaTTH[xEle], WGT)
                hst['min_lep_mvaTTH'].Fill( min(ch.Muon_mvaTTH[xMu], ch.Electron_mvaTTH[xEle]), WGT)
            
                hst['mu_miniIso']     .Fill( min(ch.Muon_miniPFRelIso_all[xMu], 0.21), WGT)
                hst['ele_miniIso']    .Fill( min(ch.Electron_miniPFRelIso_all[xEle], 0.79), WGT)
                hst['max_lep_miniIso'].Fill( min( max(ch.Muon_miniPFRelIso_all[xMu], ch.Electron_miniPFRelIso_all[xEle]), 0.79), WGT)

                hst['mu_ID'] .Fill(   ch.Muon_mediumPromptId[xMu] + \
                                    2*ch.Muon_tightId       [xMu], WGT)
                hst['ele_ID'].Fill(    ch.Electron_mvaFall17V2Iso_WPL   [xEle] + \
                                     2*ch.Electron_mvaFall17V2noIso_WP90[xEle] + \
                                     4*ch.Electron_mvaFall17V2Iso_WP90  [xEle] + \
                                     8*ch.Electron_mvaFall17V2noIso_WP80[xEle] + \
                                    16*ch.Electron_mvaFall17V2Iso_WP80  [xEle], WGT)

                hst['mu_SIP']     .Fill( min( abs(ch.Muon_sip3d[xMu]), 7.99), WGT)
                hst['ele_SIP']    .Fill( min( abs(ch.Electron_sip3d[xEle]), 7.99), WGT)
                hst['max_lep_SIP'].Fill( min( max(abs(ch.Muon_sip3d[xMu]), abs(ch.Electron_sip3d[xEle])), 7.99), WGT)

            ## pT
            hst['pt_lep1'].Fill( lepVecs[0].Pt(), WGT)
            hst['pt_lep2'].Fill( lepVecs[1].Pt(), WGT)
            hst['pt_fat'] .Fill( fatVec.Pt(), WGT)
            hst['pt_MET'] .Fill( metVec.Pt(), WGT)
            hst['pt_ll']  .Fill( ll_vec.Pt(), WGT)
            hst['pt_llv'] .Fill( llv_vec.Pt(), WGT)
            hst['pt_llJ'] .Fill( llJ_vec.Pt(), WGT)
            hst['pt_llvJ'].Fill( llvJ_vec.Pt(), WGT)
            for i in range(min(len(jetVecs),3)):
                hst['pt_jet%d' % (i+1)].Fill( jetVecs[i].Pt(), WGT)

            ## eta
            hst['eta_lep1'].Fill( lepVecs[0].Eta(), WGT)
            hst['eta_lep2'].Fill( lepVecs[1].Eta(), WGT)
            hst['eta_fat'] .Fill( fatVec.Eta(), WGT)
            hst['eta_MET'] .Fill( metVec.Eta(), WGT)
            hst['eta_ll']  .Fill( ll_vec.Eta(), WGT)
            hst['eta_llv'] .Fill( llv_vec.Eta(), WGT)
            hst['eta_llJ'] .Fill( llJ_vec.Eta(), WGT)
            hst['eta_llvJ'].Fill( llvJ_vec.Eta(), WGT)
            for i in range(min(len(jetVecs),3)):
                hst['eta_jet%d' % (i+1)].Fill( jetVecs[i].Eta(), WGT)

            hst['phi_lep1'].Fill( lepVecs[0].Phi(), WGT)
            hst['phi_lep2'].Fill( lepVecs[1].Phi(), WGT)
            hst['phi_fat'] .Fill( fatVec.Phi(), WGT)
            hst['phi_MET'] .Fill( metVec.Phi(), WGT)
            hst['phi_ll']  .Fill( ll_vec.Phi(), WGT)

            ## mass
            hst['mass_fat'] .Fill( fatVec.M(), WGT)
            hst['mass_ll']  .Fill( min(ll_vec.M(), 399.9), WGT)
            hst['mass_llv'] .Fill( llv_vec.M(), WGT)
            hst['mass_llJ'] .Fill( llJ_vec.M(), WGT)
            hst['mass_llvJ'].Fill( llvJ_vec.M(), WGT)

            ## MT
            hst['MT_llv']        .Fill( llv_vecT.M(), WGT)
            hst['MT_llvJ']       .Fill( llvJ_vecT.M(), WGT)
            hst['MT_lep1_MET']   .Fill( (lepVecTs[0]+metVecT).M(), WGT)
            hst['MT_lep2_MET']   .Fill( (lepVecTs[1]+metVecT).M(), WGT)
            hst['MT_min_lep_MET'].Fill( min( (lepVecTs[0]+metVecT).M(), (lepVecTs[1]+metVecT).M() ), WGT)
            hst['MT_max_lep_MET'].Fill( max( (lepVecTs[0]+metVecT).M(), (lepVecTs[1]+metVecT).M() ), WGT)

            ## HT
            hst['HT_ll']  .Fill( lepVecs[0].Pt()+lepVecs[1].Pt(), WGT)
            hst['HT_llv'] .Fill( lepVecs[0].Pt()+lepVecs[1].Pt()+metVec.Pt(), WGT)
            hst['HT_llJ'] .Fill( lepVecs[0].Pt()+lepVecs[1].Pt()+fatVec.Pt(), WGT)
            hst['HT_llvJ'].Fill( lepVecs[0].Pt()+lepVecs[1].Pt()+fatVec.Pt()+metVec.Pt(), WGT)

            ## Dilepton
            hst['dR_ll']  .Fill( abs( lepVecs[0].DeltaR  (lepVecs[1]) ), WGT)
            hst['dPhi_ll'].Fill( abs( lepVecs[0].DeltaPhi(lepVecs[1]) ), WGT)
            hst['dEta_ll'].Fill( abs( lepVecs[0].Eta() -  lepVecs[1].Eta() ), WGT)

            ## lep1-jet
            hst['dR_lep1_fat']  .Fill( abs( lepVecs[0].DeltaR  (fatVec) ), WGT)
            hst['dPhi_lep1_fat'].Fill( abs( lepVecs[0].DeltaPhi(fatVec) ), WGT)
            hst['dEta_lep1_fat'].Fill( abs( lepVecs[0].Eta() -  fatVec.Eta() ), WGT)

            ## lep2-jet
            hst['dR_lep2_fat']  .Fill( abs( lepVecs[1].DeltaR  (fatVec) ), WGT)
            hst['dPhi_lep2_fat'].Fill( abs( lepVecs[1].DeltaPhi(fatVec) ), WGT)
            hst['dEta_lep2_fat'].Fill( abs( lepVecs[1].Eta() -  fatVec.Eta() ), WGT)

            ## lepton-jet
            hst['mass_min_lep_fat'].Fill( min(     (lepVecs[0]         +fatVec).M(),   (lepVecs[1]         +fatVec).M() ), WGT)
            hst['mass_max_lep_fat'].Fill( max(     (lepVecs[0]         +fatVec).M(),   (lepVecs[1]         +fatVec).M() ), WGT)
            hst['dR_min_lep_fat']  .Fill( min( abs( lepVecs[0].DeltaR  (fatVec) ), abs( lepVecs[1].DeltaR  (fatVec) ) ), WGT)
            hst['dR_max_lep_fat']  .Fill( max( abs( lepVecs[0].DeltaR  (fatVec) ), abs( lepVecs[1].DeltaR  (fatVec) ) ), WGT)
            hst['dPhi_min_lep_fat'].Fill( min( abs( lepVecs[0].DeltaPhi(fatVec) ), abs( lepVecs[1].DeltaPhi(fatVec) ) ), WGT)
            hst['dPhi_max_lep_fat'].Fill( max( abs( lepVecs[0].DeltaPhi(fatVec) ), abs( lepVecs[1].DeltaPhi(fatVec) ) ), WGT)
            hst['dEta_min_lep_fat'].Fill( min( abs( lepVecs[0].Eta() -  fatVec.Eta() ), abs( lepVecs[1].Eta() - fatVec.Eta() ) ), WGT)
            hst['dEta_max_lep_fat'].Fill( max( abs( lepVecs[0].Eta() -  fatVec.Eta() ), abs( lepVecs[1].Eta() - fatVec.Eta() ) ), WGT)

            ## dilepton-jet
            hst['dR_ll_fat']  .Fill( abs( ll_vec.DeltaR  (fatVec) ), WGT)
            hst['dEta_ll_fat'].Fill( abs( ll_vec.Eta() -  fatVec.Eta() ), WGT)
            hst['dPhi_ll_fat'].Fill( abs( ll_vec.DeltaPhi(fatVec) ), WGT)

            ## dPhi with MET
            hst['dPhi_ll_MET'] .Fill( abs( ll_vec       .DeltaPhi(metVec) ), WGT)
            hst['dPhi_llJ_MET'].Fill( abs( llJ_vec.DeltaPhi(metVec) ), WGT)

            ## ISR plots
            hst['pt_ISRs'].Fill( isrsVec.Pt(), WGT)
            hst['HT_ISRs'].Fill( isrs_HT, WGT)

            if len(jetVecs) > 0:
                hst['dR_llJ_ISR']   .Fill( llJ_vec .DeltaR  (isrVec), WGT)
                hst['dPhi_llvJ_ISR'].Fill( llvJ_vec.DeltaPhi(isrVec), WGT)
                hst['dPhi_llJ_ISR'] .Fill( llJ_vec .DeltaPhi(isrVec), WGT)
                hst['pt_llvJ_ISR']  .Fill( (llvJ_vec+isrVec).Pt(), WGT)
                hst['pt_llJ_ISR']   .Fill( (llJ_vec +isrVec).Pt(), WGT)
                hst['mass_llvJ_ISR'].Fill( (llvJ_vec+isrVec).M(), WGT)
                hst['mass_llJ_ISR'] .Fill( (llJ_vec +isrVec).M(), WGT)
                hst['MT_llvJ_ISR']  .Fill( (llvJ_vecT+isrVecT).M(), WGT)
                hst['MT_llJ_ISR']   .Fill( (llJ_vecT +isrVecT).M(), WGT)
                hst['HT_llvJ_ISR']  .Fill(  lepVecs[0].Pt()+lepVecs[1].Pt()+fatVec.Pt()+metVec.Pt()+isrVec.Pt(), WGT)
                hst['HT_llJ_ISR']   .Fill(  lepVecs[0].Pt()+lepVecs[1].Pt()+fatVec.Pt()            +isrVec.Pt(), WGT)

                hst['dR_llJ_ISRs']   .Fill( llJ_vec .DeltaR  (isrsVec), WGT)
                hst['dPhi_llvJ_ISRs'].Fill( llvJ_vec.DeltaPhi(isrsVec), WGT)
                hst['dPhi_llJ_ISRs'] .Fill( llJ_vec .DeltaPhi(isrsVec), WGT)
                hst['pt_llvJ_ISRs']  .Fill( (llvJ_vec+isrsVec).Pt(), WGT)
                hst['pt_llJ_ISRs']   .Fill( (llJ_vec +isrsVec).Pt(), WGT)
                hst['mass_llvJ_ISRs'].Fill( (llvJ_vec+isrsVec).M(), WGT)
                hst['mass_llJ_ISRs'] .Fill( (llJ_vec +isrsVec).M(), WGT)
                hst['MT_llvJ_ISRs']  .Fill( (llvJ_vecT+isrsVecT).M(), WGT)
                hst['MT_llJ_ISRs']   .Fill( (llJ_vecT +isrsVecT).M(), WGT)
                hst['HT_llvJ_ISRs']  .Fill(  lepVecs[0].Pt()+lepVecs[1].Pt()+fatVec.Pt()+metVec.Pt()+isrs_HT, WGT)
                hst['HT_llJ_ISRs']   .Fill(  lepVecs[0].Pt()+lepVecs[1].Pt()+fatVec.Pt()            +isrs_HT, WGT)


            hst['deepB_max_jet'].Fill( max(-0.099, min(0.999, max_deepB) ), WGT)
            hst['flavB_max_jet'].Fill( max(-0.099, min(0.999, max_flavB) ), WGT)
            hst['doubB_fat']    .Fill( max(-0.099, min(0.999, ch.FatJet_btagDDBvLV2[fatIdx]) ), WGT)
            hst['msoft_fat']    .Fill( ch.FatJet_msoftdrop[fatIdx], WGT)

            hst['nSV_max_jet'].Fill(nSV_max_jet, WGT)
            hst['nSV_ISR']    .Fill(nSV_ISR, WGT)
            hst['nSV_ISRs']   .Fill(nSV_ISRs, WGT)

        ## End loop: for jFat in range(len(fatVecs))
    ## End loop: for iEvt in range(nEntries):

    print('\nFinished loop over events')
    
    print('\nOut of %d pre-selected events, %d pass b-tag veto, %d mass cut, %d lep, %d MET, %d (%d, %d) pass (tight, data) kinematic cuts, %d trigger, %d GEN' % (nPassPre, nPassBtag, nPassMass, nPassLep, nPassMET, nPassKin, nPassTight, nPassData, nPassTrg, nPassGen))

    out_file.cd()

    keys_to_delete = []
    for key in hst.keys():
        if hst[key].Integral() == 0:
            keys_to_delete.append(key)
            continue
        hst[key].SetLineColor(R.kBlue if CUT_GEN_BB else R.kBlack)
        hst[key].SetLineWidth(2)

    for key in keys_to_delete:
        del hst[key]
        
    print('\nSaved histograms to output file %s\n' % out_file_str)

    out_file.Write()
    out_file.Close()


if __name__ == '__main__':
    main()
