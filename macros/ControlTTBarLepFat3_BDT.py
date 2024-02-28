#! /usr/bin/env python

## Plots for ttbar --> lepton + AK8 jet control region

import os
import sys
import subprocess
import numpy as np
import ROOT as R

R.gROOT.SetBatch(True)  ## Don't display histograms or canvases when drawn

## User configuration
PRT_EVT = 1000     ## Print every Nth event
MAX_EVT = 100000   ## Number of events to process
VERBOSE = False    ## Verbose print-outs
LUMI    = 59830    ## 2018 integrated lumi (pb-1), certified "Good"

PRT_BDT_IN   = False   ## Print input variables for BDT training in CSV format  
SKIM_H4B     = True    ## Use input files from skims/Hto4b_0p8/
CUT_GEN_BBQQ = 0       ## 4 for BBQQ, 3 for BBQ, 2 for BB, 1 for B, 0 for none, -1 for not-BB
CUT_MASS     = 140     ## Minimum AK8 jet PF mass
CUT_MSOFT    = 20      ## Minimum AK8 jet soft-drop mass
CUT_DR_LJ    = -1      ## Maximum dR(AK8, lepton)
CUT_PT_LVJ   = 0       ## Minimum pT(AK8, lep, MET)
CUT_MT_LVJ   = 700     ## Maximum MT(AK8, lep, MET)
CUT_HT_ISR   = 300     ## Minimum HT(AK4 jets)
CUT_PT_ISR   = 250     ## Minimum pT(AK4 jets)
CUT_LFRAC    = 0.6     ## Minimum lepton/AK4 jet pT ratio for events with dR(lep,jet) < 0.4
LABEL = 'TTToSemiLep'
#LABEL = 'TTTo2L2Nu'
#LABEL = 'WToLNu_HT-600'
CAT   = 'SingleLep'

if len(sys.argv) > 1:
    print('\nLABEL changed from %s to %s' % (LABEL, str(sys.argv[1])))
    LABEL = str(sys.argv[1])
if len(sys.argv) > 2:
    print('\nCAT changed from %s to %s' % (CAT, str(sys.argv[2])))
    CAT = str(sys.argv[2])
if len(sys.argv) > 3:
    print('\nCUT_GEN_BBQQ changed from %d to %d' % (CUT_GEN_BBQQ, int(sys.argv[3])))
    CUT_GEN_BBQQ = int(sys.argv[3])

## Small function to both fill histogram and print out value for BDT
#BDT_VARS = ['pt_lep','pt_fat','pt_MET','pt_ISRs','pt_lJ','pt_lv','pt_lvJ','pt_jet1','pt_jet2','pt_jet3','dR_lep_fat','dR_lJ_ISR','dEtaS_lep_fat','dPhi_lep_fat','dPhi_lep_MET','dPhi_lvJ_ISR','mass_fat','mass_lJ','HT_ISRs','HT_lvJ','msoft_fat','deepB_max_jet','flavB_max_jet','dR_fat_jet_min']
#BDT_VARS = ['msoft_fat','mass_fat','mass_lJ','dR_lep_fat','dEtaS_lep_fat','dPhi_lep_fat','deepB_max_jet','flavB_max_jet','dR_lJ_deepB_max','dR_lJ_flavB_max']
BDT_VARS = ['pt_lep','pt_fat','pt_MET','pt_ISRs','pt_lJ','pt_lJ_ISR','pt_lJ_ISRs','pt_lv','pt_lvJ','pt_lvJ_ISR','pt_lvJ_ISRs','pt_jet1','pt_jet2','pt_jet3','eta_lep','eta_fat','eta_lJ','eta_jet1','eta_jet2','eta_jet3','etaA_lep','etaA_fat','etaA_lJ','etaA_jet1','etaA_jet2','etaA_jet3','etaS_lep','etaS_fat','etaS_jet1','etaS_jet2','etaS_jet3','lepjet_frac','lepjet_dR','dR_lep_fat','dR_lJ_ISR','dR_lJ_ISRs','dEta_lep_fat','dEtaS_lep_fat','dPhi_lep_fat','dPhi_lep_MET','dPhi_lJ_MET','dPhi_lJ_ISR','dPhi_lJ_ISRs','dPhi_lv_fat','dPhi_lvJ_ISR','dPhi_lvJ_ISRs','mass_fat','mass_lJ','mass_lJ_ISR','mass_lJ_ISRs','mass_lvJ','mass_lvJ_ISR','mass_lvJ_ISRs','MT_lJ_ISR','MT_lJ_ISRs','MT_lv','MT_lvJ','MT_lvJ_ISR','MT_lvJ_ISRs','HT_ISRs','HT_lJ','HT_lJ_ISR','HT_lJ_ISRs','HT_lv','HT_lvJ','HT_lvJ_ISR','HT_lvJ_ISRs','nFat','nJet','msoft_fat','deepB_max_jet','flavB_max_jet','dR_fat_jet_min','dR_lJ_jet_min','dR_fat_deepB_max','dR_lJ_deepB_max','dR_fat_flavB_max','dR_lJ_flavB_max','jet_pt_near_fat','jet_pt_near_lJ','deepB_near_fat','deepB_near_lJ','flavB_near_fat','flavB_near_lJ','bb_DDBvLV2_fat','bb_PNet_Xbb_fat','bb_deep_ZHbb_fat','bb_batg_Hbb_fat','Hto4b_Haa4b_fat','Hto4b_Haa34b_fat','wgt']

def FillX(csv, cat, hists, name, val, wgt):
    hists['%s_0'  % (name     )].Fill(val, wgt)  ## Fill inclusive histogram
    hists['%s_%d' % (name, cat)].Fill(val, wgt)  ## Fill histogram for specific category
    if PRT_BDT_IN and name in BDT_VARS:
        if val < 0: csv[name] = (', '+('{:.99f}'.format(val))[:max(7,len(name))])
        else:       csv[name] = (',  '+('{:.99f}'.format(val))[:max(6,len(name)-1)])
    return csv
## End function: FillX(csv, hists, name, val, wgt)

def main():

    print('\nInside ControlTTBarLepFat3\n')

    in_file_names = []
    in_dir = '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/'
    if LABEL == 'SingleMuon':  in_dir += 'data/PNet_v1_2023_10_06/Run2018D-UL2018_MiniAODv2_GT36-v2/SingleMuon/'
    if LABEL == 'EGamma':      in_dir += 'data/PNet_v1_2023_10_06/Run2018D-UL2018_MiniAODv2_GT36-v3/EGamma/'
    if LABEL == 'TTToSemiLep': in_dir += 'MC/PNet_v1_2023_10_06/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/'
    if LABEL == 'TTTo2L2Nu':   in_dir += 'MC/PNet_v1_2023_10_06/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/'
    if LABEL == 'ST_s_lep':    in_dir += 'MC/PNet_v1_2023_10_06/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/'
    if LABEL == 'ST_tW_top':   in_dir += 'MC/PNet_v1_2023_10_06/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/'
    if LABEL == 'ST_tW_anti':  in_dir += 'MC/PNet_v1_2023_10_06/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/'
    if LABEL == 'ST_t_top':    in_dir += 'MC/PNet_v1_2023_10_06/ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/'
    if LABEL == 'ST_t_anti':   in_dir += 'MC/PNet_v1_2023_10_06/ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/'
    if LABEL == 'WToLNu_NLO':  in_dir += 'MC/PNet_v1_2023_10_06/WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/'
    if LABEL == 'WToLNu_HT-600': in_dir += 'MC/PNet_v1_2023_10_06/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/'
    if LABEL == 'WW':          in_dir += 'MC/PNet_v1_2023_10_06/WW_TuneCP5_13TeV-pythia8/'
    if LABEL == 'WZ':          in_dir += 'MC/PNet_v1_2023_10_06/WZ_TuneCP5_13TeV-pythia8/'
    ## TTJets_Incl, TTJets_Incl_Had, TTJets_Incl_SemiLep, TTJets_Incl_2L2Nu
    ## ZJets_MG, ZJets_AMC, ST_tW_top, ST_tW_antitop, WW, WZ, ZZ
    if SKIM_H4B: in_dir += 'skims/Hto4b_0p8/'
    else:        in_dir += 'r1/'

    print('\nGetting input files from %s' % in_dir)

    isData  = ('/data/' in in_dir)
    isTTBar = (LABEL.startswith('TTTo') or LABEL.startswith('TTJets')) 

    if not LABEL.split('_')[0] in in_dir and not LABEL.split('_')[0] in in_dir.replace('Jets',''):
        print('\n\n***  TRIED TO APPLY LABEL %s TO DIRECTORY %s  ***' % (LABEL, in_dir))
        print('***  CORRECT? QUITTING!  ***')
        sys.exit()
    if not isTTBar and CUT_GEN_BBQQ >= 2:
        print('\n\n***  TRIED TO REQUIRE CUT_GEN_BBQQ %s TO DIRECTORY %s  ***' % (CUT_GEN_BBQQ, in_dir))
        print('***  CORRECT? QUITTING!  ***')
        sys.exit()

    for f_name in subprocess.check_output(['ls', in_dir], encoding='UTF-8').splitlines():
        if not '.root' in str(f_name): continue
        in_file_names.append(in_dir+f_name)
        print('Appending file: %s' % in_file_names[-1])
        if MAX_EVT > 0 and len(in_file_names)*50001*(1 + 19*(SKIM_H4B==0)) > MAX_EVT: break

    out_dir = 'plots/'
    out_file_str = out_dir+'ControlTTBarLepFat_%s_%s' % (CAT, LABEL)
    if   SKIM_H4B:           out_file_str += '_Hto4b_0p8' 
    if   CUT_GEN_BBQQ == -1: out_file_str += '_noGenBB'
    elif CUT_GEN_BBQQ ==  1: out_file_str += '_GenB'
    elif CUT_GEN_BBQQ ==  2: out_file_str += '_GenBB'
    elif CUT_GEN_BBQQ ==  3: out_file_str += '_GenBBQ'
    elif CUT_GEN_BBQQ ==  4: out_file_str += '_GenBBQQ'
    else:                     out_file_str += '_GenAll'
    if CUT_MASS   > 0: out_file_str += ('_mass_%d'  % CUT_MASS)
    if CUT_MSOFT  > 0: out_file_str += ('_msoft_%d' % CUT_MSOFT)
    if CUT_DR_LJ  > 0: out_file_str += ('_dRLJ_%d'  % (CUT_DR_LJ*100))
    if CUT_PT_LVJ > 0: out_file_str += ('_ptLVJ_%d' % CUT_PT_LVJ)
    if CUT_MT_LVJ > 0: out_file_str += ('_MtLVJ_%d' % CUT_MT_LVJ)
    if CUT_HT_ISR > 0: out_file_str += ('_HtISR_%d' % CUT_HT_ISR)
    if CUT_PT_ISR > 0: out_file_str += ('_ptISR_%d' % CUT_PT_ISR)
    if CUT_LFRAC  > 0: out_file_str += ('_lfrac_%d' % (CUT_LFRAC*100))
    if MAX_EVT    > 9999999: out_file_str += '_%dM' % (MAX_EVT / 1000000)
    elif MAX_EVT  >       0: out_file_str += '_%dk' % (MAX_EVT / 1000)
    # else:           out_file_str += '_1188k'
    out_file_str += '.root'
    out_file = R.TFile(out_file_str,'recreate')
    if PRT_BDT_IN:
        bdt_file = open(out_file_str.replace('plots','csv').replace('.root','.csv'), 'w')

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
    bins['pt']    = [200,  0, 1000]
    bins['eta']   = [200, -5,    5]
    bins['etaA']  = [100,  0,    5]  ## Absolute value of eta
    bins['etaS']  = [200, -5,    5]  ## "Signed" eta (relative to lep + AK8 jet vector)
    bins['phi']   = [128, -3.2, 3.2]
    bins['sum']   = [100,  0, 1000]
    bins['dR']    = [100,  0,   10]
    bins['dEta']  = [100,  0,   10]  ## Absolute value of dEta
    bins['dEtaS'] = [200, -10,  10]  ## "Signed" dEta (positive if abs(eta) of obj1 > obj2)
    bins['dPhi']  = [ 32,  0,  3.2]
    bins['mass']  = [200,  0, 2000]
    bins['MT']    = [200,  0, 2000]
    bins['HT']    = [200,  0, 2000]


    ## Sets of objects to be plotted
    Objs = ['lep','fat','MET','ISR','ISRs','lJ','lv','lvJ','jet1','jet2','jet3']
    Vars = ['pt','eta','etaA','etaS','phi','dR','dEta','dEtaS','dPhi','mass','MT','HT']

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
    hst['nQC_vs_B']  = R.TH2D('h_nQC_vs_B',  'h_nQC_vs_B',  4, -0.5, 3.5, 6, -0.5, 5.5)
    hst['nCh_vs_Bh'] = R.TH2D('h_nCh_vs_Bh', 'h_nCh_vs_Bh', 6, -0.5, 5.5, 9, -0.5, 8.5)
    hst['lepjet_frac_vs_dR'] = R.TH2D('h_lepjet_frac_vs_dR', 'h_lepjet_frac_vs_dR', 50, 0, 0.5, 55, 0, 1.1)

    hst['mu_mvaTTH']  = R.TH1D('h_mu_mvaTTH',  'h_mu_mvaTTH',  40, -1, 1)
    hst['ele_mvaTTH'] = R.TH1D('h_ele_mvaTTH', 'h_ele_mvaTTH', 40, -1, 1)
    hst['lep_mvaTTH'] = R.TH1D('h_lep_mvaTTH', 'h_lep_mvaTTH', 40, -1, 1)

    hst['mu_miniIso']  = R.TH1D('h_mu_miniIso',  'h_mu_miniIso',  11, 0, 0.22)
    hst['ele_miniIso'] = R.TH1D('h_ele_miniIso', 'h_ele_miniIso', 16, 0, 0.8)
    hst['lep_miniIso'] = R.TH1D('h_lep_miniIso', 'h_lep_miniIso', 16, 0, 0.8)

    hst['mu_ID']  = R.TH1D('h_mu_ID',  'h_mu_ID',   4, -0.5,  3.5)
    hst['ele_ID'] = R.TH1D('h_ele_ID', 'h_ele_ID', 32, -0.5, 31.5)

    hst['mu_SIP']  = R.TH1D('h_mu_SIP',  'h_mu_SIP',  16, 0, 8.0)
    hst['ele_SIP'] = R.TH1D('h_ele_SIP', 'h_ele_SIP', 16, 0, 8.0)
    hst['lep_SIP'] = R.TH1D('h_lep_SIP', 'h_lep_SIP', 16, 0, 8.0)

    hst['lepjet_frac'] = R.TH1D('h_lepjet_frac', 'h_lepjet_frac', 110, 0, 1.1)
    hst['lepjet_dR']   = R.TH1D('h_lepjet_dR',   'h_lepjet_dR',   50, 0, 0.5)
    
    hst['msoft_fat'] = R.TH1D('h_msoft_fat', 'fat msoft', 100, 0, 500)

    hst['deepB_max_jet'] = R.TH1D('h_deepB_max_jet', 'Maximum AK4 jet deepB tag score', 110, -0.1, 1.0)
    hst['flavB_max_jet'] = R.TH1D('h_flavB_max_jet', 'Maximum AK4 jet deepFlavB tag score', 110, -0.1, 1.0)
    hst['nSV_max_jet']   = R.TH1D('h_nSV_max_jet',   'Maximum AK4 jet # of secondary vertices', 7, -0.5, 6.5)

    hst['dR_fat_jet_min']   = R.TH1D('h_dR_fat_jet_min', 'Minimum dR(AK4, AK8)', 100, 0, 10)
    hst['dR_lJ_jet_min']    = R.TH1D('h_dR_lJ_jet_min', 'Minimum dR(AK4, lep+AK8)', 100, 0, 10)
    hst['dR_fat_deepB_max'] = R.TH1D('h_dR_fat_deepB_max', 'dR(high deepB, AK8)', 100, 0, 10)
    hst['dR_lJ_deepB_max']  = R.TH1D('h_dR_lJ_deepB_max', 'dR(high deepB, lep+AK8)', 100, 0, 10)
    hst['dR_fat_flavB_max'] = R.TH1D('h_dR_fat_flavB_max', 'dR(high flavB, AK8)', 100, 0, 10)
    hst['dR_lJ_flavB_max']  = R.TH1D('h_dR_lJ_flavB_max', 'dR(high flavB, lep+AK8)', 100, 0, 10)

    hst['jet_pt_near_fat'] = R.TH1D('h_jet_pt_near_fat', 'pT of AK4 nearest AK8', 60, 0, 300)
    hst['jet_pt_near_lJ']  = R.TH1D('h_jet_pt_near_lJ', 'pT of AK4 nearest lep+AK8', 60, 0, 300)
    hst['deepB_near_fat'] = R.TH1D('h_deepB_near_fat', 'deepB nearest AK8', 110, -0.1, 1.0)
    hst['deepB_near_lJ']  = R.TH1D('h_deepB_near_lJ', 'deepB nearest lep+AK8', 110, -0.1, 1.0)
    hst['flavB_near_fat'] = R.TH1D('h_flavB_near_fat', 'flavB nearest AK8', 110, -0.1, 1.0)
    hst['flavB_near_lJ']  = R.TH1D('h_flavB_near_lJ', 'flavB nearest lep+AK8', 110, -0.1, 1.0)
    
    hst['nSV_ISR']       = R.TH1D('h_nSV_ISR',  'ISR jet # of secondary vertices',   7, -0.5,  6.5)
    hst['nSV_ISRs']      = R.TH1D('h_nSV_ISRs', 'ISR jets # of secondary vertices', 13, -0.5, 12.5)

    hst['bb_DDBvLV2_fat']   = R.TH1D('h_bb_DDBvLV2_fat',   'AK8 jet btagDDBvLV2 score', 110, -0.1, 1.0)
    hst['bb_PNet_Xbb_fat']  = R.TH1D('h_bb_PNet_Xbb_fat',  'AK8 jet ParticleNetMD Xbb/QCD score', 110, -0.1, 1.0)
    hst['bb_deep_ZHbb_fat'] = R.TH1D('h_bb_deep_ZJbb_fat', 'AK8 jet DeepTagMD ZHbb score', 110, -0.1, 1.0)
    hst['bb_batg_Hbb_fat']  = R.TH1D('h_bb_btag_Hbb_fat',  'AK8 jet btagHbb score', 110, -0.1, 1.0)

    hst['qq_PNet_Xqq_fat']  = R.TH1D('h_qq_PNet_Xqq_fat',  'AK8 jet ParticleNetMD Xqq/QCD score', 110, -0.1, 1.0)
    hst['4q_PNet_H4q_fat']  = R.TH1D('h_4q_PNet_H4q_fat',  'AK8 jet ParticleNetMD H4q/QCD score', 110, -0.1, 1.0)
    
    hst['TvsQCD_deepMD_fat'] = R.TH1D('h_TvsQCD_deepMD_fat', 'AK8 jet DeepTagMD top/QCD score', 110, -0.1, 1.0)
    hst['TvsQCD_deep_fat']   = R.TH1D('h_TvsQCD_deep_fat',   'AK8 jet DeepTag top/QCD score', 110, -0.1, 1.0)
    hst['TvsQCD_PNet_fat']   = R.TH1D('h_TvsQCD_PNet_fat',   'AK8 jet ParticleNet top/QCD score', 110, -0.1, 1.0)
    hst['WvsQCD_deepMD_fat'] = R.TH1D('h_WvsQCD_deepMD_fat', 'AK8 jet DeepTagMD W/QCD score', 110, -0.1, 1.0)
    hst['WvsQCD_deep_fat']   = R.TH1D('h_WvsQCD_deep_fat',   'AK8 jet DeepTag W/QCD score', 110, -0.1, 1.0)
    hst['WvsQCD_PNet_fat']   = R.TH1D('h_WvsQCD_PNet_fat',   'AK8 jet ParticleNet W/QCD score', 110, -0.1, 1.0)
    
    hst['Hto4b_Haa4b_fat']    = R.TH1D('h_Hto4b_Haa4b_fat',    'AK8 jet Hto4b Haa4b score', 101, -0.01, 1.0)
    hst['Hto4b_Haa3b_fat']    = R.TH1D('h_Hto4b_Haa3b_fat',    'AK8 jet Hto4b Haa3b score', 101, -0.01, 1.0)
    hst['Hto4b_Haa34b_fat']   = R.TH1D('h_Hto4b_Haa34b_fat',   'AK8 jet Hto4b Haa3b+Haa4b score', 101, -0.01, 1.0)
    hst['Hto4b_binary_fat']   = R.TH1D('h_Hto4b_binary_fat',   'AK8 jet Hto4b binary score', 101, -0.01, 1.0)
    hst['Hto4b_binaryLF_fat'] = R.TH1D('h_Hto4b_binaryLF_fat', 'AK8 jet Hto4b binary LF score', 101, -0.01, 1.0)

    hst['wgt']      = R.TH1D('h_wgt',      'Event weights (unweighted)', 2000, -10, 10.)
    hst['wgt_wgt']  = R.TH1D('h_wgt_wgt',  'Event weights (weighted)',   2000, -10, 10.)
    hst['nPV']      = R.TH1D('h_nPV',      '# of PVs',      101, -0.5, 100.5)
    hst['nPV_good'] = R.TH1D('h_nPV_good', '# of good PVs', 101, -0.5, 100.5)

    ## Duplicate histograms for different categories
    new_hst = {}
    for key in hst.keys():
        for BBQQ_cat in range(10):  ## See X2 variable below 
            new_hst['%s_%d' % (key, BBQQ_cat)] = hst[key].Clone('%s_%d' % (hst[key].GetName(), BBQQ_cat))
    hst.update(new_hst)
    del new_hst

    ## 2018 DeepCSV and DeepFlavB cuts: https://btv-wiki.docs.cern.ch/ScaleFactors/UL2018/
    DeepB = {'L': 0.1208, 'M': 0.4168, 'T': 0.7665}
    FlavB = {'L': 0.0490, 'M': 0.2783, 'T': 0.7100}

    ## Print header for BDT input CSV file
    bcsv_head = 'BDT_IN, BBQQ'
    if PRT_BDT_IN:
        for BDT_var in BDT_VARS:
            bcsv_head += (', '+'{0: >7}'.format(BDT_var))
        bdt_file.write('\n%s' % bcsv_head)
        
    ## Loop through events, select, and plot
    nEntries = chains['Events'].GetEntries()
    print('\nEntering loop over %d events\n' % (nEntries))

    ch = chains['Events']  ## Shortcut expression
    ## Count passing events
    nPassPre   = 0
    nPassBtag  = 0
    nPassLep   = 0
    nPassMET   = 0
    nPassKin   = 0
    nPassTight = 0
    nPassData  = 0
    nPassTrg   = 0
    nPassGen   = 0
    nPassLepFrac = 0
    for iEvt in range(nEntries):

        ch.GetEntry(iEvt)

        if iEvt > MAX_EVT and MAX_EVT > 0: break
        if (iEvt % PRT_EVT) == 0: print('Looking at event #%d / %d' % (iEvt, nEntries))

        if ch.nFatJet < 1: continue

        # if VERBOSE: print('  * Pass ch.nFatJet < 1 cut')

        if ch.nMuon + ch.nElectron            < 1: continue
        if CAT == 'SingleMu' and ch.nMuon     < 1: continue
        if CAT == 'SingleEG' and ch.nElectron < 1: continue

        # if VERBOSE: print('  * Pass ch.nMuon + ch.nElectron < 1 cut')
        
        ## Require qualifying AK8 jet off the bat
        hasGoodFat = False
        for iFat in range(ch.nFatJet):
            if ch.FatJet_mass[iFat]      < CUT_MASS:  continue
            if ch.FatJet_msoftdrop[iFat] < CUT_MSOFT: continue
            if      ch.FatJet_pt [iFat]  < 170: continue
            if  abs(ch.FatJet_eta[iFat]) > 2.4: continue
            if ch.FatJet_jetId[iFat]     < 6:   continue
            hasGoodFat = True
            break
        if not hasGoodFat: continue

        if VERBOSE: print('  * Pass hasGoodFat cut')

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
        cIdxs   = []  ## Indices of c-quarks from the ttbar decay
        qIdxs   = []  ## Indices of uds-quarks from the ttbar decay

        muVecs   = []  ## TLorentzVectors of selected muons
        muVecTs  = []  ## TLorentzVectors of selected muons
        eleVecs  = []  ## TLorentzVectors of selected electrons
        eleVecTs = []  ## TLorentzVectors of selected electrons
        lepVecs  = []  ## TLorentzVectors of two selected leptons
        lepVecTs = []  ## TLorentzVectors of two selected leptons
        fatVecs  = []  ## TLorentzVectors of selected AK8 jets
        fatVecTs = []  ## TLorentzVectors of selected AK8 jets
        bVecs    = []  ## TLorentzVectors of b-quarks from ttbar decay
        cVecs    = []  ## TLorentzVectors of c-quarks from ttbar decay
        qVecs    = []  ## TLorentzVectors of uds-quarks from ttbar decay

        ## In ttbar events, look for leptons and quarks from W decay
        nGenLep = 0
        nGenC   = 0
        nGenQ   = 0
        lastPdgID = -999
        genLepMoms = []
        for iGen in range(-1 if not isTTBar else ch.nGenPart):
            pdgID = ch.GenPart_pdgId[iGen]
            ## Particle must be a charged lepton or quark
            if not abs(pdgID) in [1,2,3,4,5,11,13,15]: continue
            ## Particle must have W boson mother
            iMom = ch.GenPart_genPartIdxMother[iGen]
            if iMom < 0: continue
            if abs(ch.GenPart_pdgId[iMom]) != 24: continue
            ## Get 4-vector for particle
            iVec = R.TLorentzVector()
            iVec.SetPtEtaPhiM( ch.GenPart_pt  [iGen],
                               ch.GenPart_eta [iGen],
                               ch.GenPart_phi [iGen],
                               ch.GenPart_mass[iGen] )
            ## Count leptons from W decay
            if abs(pdgID) in [11,13,15]:
                if  iMom in genLepMoms:
                    ## Special case of W radiating a gamma --> l+l- pair
                    if iMom == genLepMoms[-1] and pdgID == -1*lastPdgID:
                        nGenLep -= 1  ## Delete previous lepton, and skip this one
                    continue
                else : genLepMoms.append(iMom)
                nGenLep += 1
                lastPdgID = pdgID
            ## Count quarks from W decay
            elif abs(pdgID) == 4:
                nGenC += 1
                cIdxs.append(iGen)
                cVecs.append(iVec)
            elif abs(pdgID) <= 3:
                nGenQ += 1
                qIdxs.append(iGen)
                qVecs.append(iVec)
        ## End loop: for iGen in range(-1 if not isTTBar else ch.nGenPart)
        if nGenLep > 2:
            print('\n*** Super-weird event!!! LS = %d, event = %d has %d leptons from W decay! ***\n' % (ch.luminosityBlock, ch.event, nGenLep))
        if 'TTToSemiLep' in LABEL and nGenLep != 1:
            print('\n*** Super-weird event!!! LS = %d, event = %d has %d leptons from W decay! ***\n' % (ch.luminosityBlock, ch.event, nGenLep))
            sys.exit()
        if 'TTTo2L2Nu' in LABEL and nGenLep != 2:
            print('\n*** Super-weird event!!! LS = %d, event = %d has %d leptons from W decay! ***\n' % (ch.luminosityBlock, ch.event, nGenLep))
            sys.exit()
        # print('Event has %d leptons, %d charm, and %d light-flavor quarks from W decay' % (nGenLep, nGenC, nGenQ))

        ## Check that event *has* GEN quarks from W decay
        if 'TTToSemiLep' in LABEL and len(cVecs) + len(qVecs) < (CUT_GEN_BBQQ - 2):
            print('\n*** Super-weird event!!! LS = %d, event = %d has %d quarks from W decay! ***\n' % (ch.luminosityBlock, ch.event, nGenC+nGenQ))

        ## For inclusive ttbar sample, separate into di-lepton, semi-leptonic, and hadronic samples
        if 'TTJets_Incl_2L2Nu' in LABEL and nGenLep < 2:
            continue
        if 'TTJets_Incl_SemiLep' in LABEL and nGenLep != 1:
            continue
        if 'TTJets_Incl_Had' in LABEL and nGenLep > 0:
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
               abs(ch.GenPart_eta[iGen] - bVecs[-1].Eta())       < 0.1:
                print('\nIn LS %d, event %d, apparent duplicate b-hadrons with pT %.2f/%.2f, eta %.3f/%.3f, phi %.3f/%.3f' % \
                      (ch.luminosityBlock, ch.event, bVecs[-1].Pt(), ch.GenPart_pt[iGen], bVecs[-1].Eta(), ch.GenPart_eta[iGen], \
                       bVecs[-1].Phi(), ch.GenPart_phi[iGen]))
                continue

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


        ## Check that event *has* exactly 2 GEN b-quarks
        if CUT_GEN_BBQQ >= 2 and len(bVecs) != 2:
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
        ## End conditional: if CUT_GEN_BBQQ >= 2 and len(bVecs) != 2

        ## Skip "signal" events if GEN quarks are too far away from each other
        if CUT_GEN_BBQQ >= 2 and bVecs[0].DeltaR(bVecs[1]) > 1.6: continue
        if CUT_GEN_BBQQ >= 3:
            nQinBB = 0
            for xVec in cVecs+qVecs:
                if max(bVecs[0].DeltaR(xVec), bVecs[1].DeltaR(xVec)) < 1.6:
                    nQinBB += 1
            if nQinBB < (CUT_GEN_BBQQ - 2): continue

        if VERBOSE: print('  * Pass GEN_BBQQ cuts')

        ## Find muon(s) passing cuts (following HIG-21-002, AN-2020/032 Table 7 on p. 18)
        for iMu in range(ch.nMuon):
            if           ch.Muon_mvaTTH[iMu]  < 0.5 : continue
            if              ch.Muon_pt [iMu]  < 10. : continue
            if          abs(ch.Muon_eta[iMu]) > 2.4 : continue
            if         ch.Muon_mediumId[iMu]  != 1  : continue
            if ch.Muon_miniPFRelIso_all[iMu]  > 0.40: continue
            if        abs(ch.Muon_sip3d[iMu]) > 8.0 : continue
            if          abs(ch.Muon_dxy[iMu]) > 0.05: continue
            if           abs(ch.Muon_dz[iMu]) > 0.10: continue

            ## Save the selected muon
            muVec  = R.TLorentzVector()
            muVecT = R.TLorentzVector()
            muVec .SetPtEtaPhiM( ch.Muon_pt[iMu], ch.Muon_eta[iMu], ch.Muon_phi[iMu], 0.106 )
            muVecT.SetPtEtaPhiM( ch.Muon_pt[iMu],                0, ch.Muon_phi[iMu], 0.106 )
            muVecs .append(muVec)
            muVecTs.append(muVecT)
            muIdxs .append(iMu)
        ## End loop: for iMu in range(ch.nMuon)

        ## Find electron(s) passing cuts (following HIG-21-002, AN-2020/032 Table 6 on p. 17)
        for iEle in range(ch.nElectron):
            if                ch.Electron_mvaTTH[iEle]  < 0.3 : continue
            if                   ch.Electron_pt [iEle]  < 10. : continue
            if               abs(ch.Electron_eta[iEle]) > 2.5 : continue
            if ch.Electron_mvaFall17V2noIso_WP90[iEle]  != 1  : continue
            if      ch.Electron_miniPFRelIso_all[iEle]  > 0.40: continue
            if             abs(ch.Electron_sip3d[iEle]) > 8.0 : continue
            if               abs(ch.Electron_dxy[iEle]) > 0.05: continue
            if                abs(ch.Electron_dz[iEle]) > 0.10: continue

            ## Save the selected electron
            eleVec  = R.TLorentzVector()
            eleVecT = R.TLorentzVector()
            eleVec .SetPtEtaPhiM( ch.Electron_pt[iEle], ch.Electron_eta[iEle], ch.Electron_phi[iEle], 0.0005 )
            eleVecT.SetPtEtaPhiM( ch.Electron_pt[iEle],                     0, ch.Electron_phi[iEle], 0.0005 )
            eleVecs .append(eleVec)
            eleVecTs.append(eleVecT)
            eleIdxs .append(iEle)
        ## End loop: for iEle in range(ch.nElectron)


        ## Only keep events with exactly one selected lepton
        if len(muIdxs) + len(eleIdxs) != 1: continue
        if VERBOSE: print('  * Pass len(muIdxs) + len(eleIdxs) != 1 cut')

        ## Define lepton flavor category event-by-event
        if   len(muIdxs)  > 0:
            LepCat = 'SingleMu'
            xMu = muIdxs[0]
            lepVecs = muVecs
            lepVecTs = muVecTs
        elif len(eleIdxs) > 0:
            LepCat = 'SingleEG'
            xEle = eleIdxs[0]
            lepVecs  = eleVecs
            lepVecTs = eleVecTs

        ## Store MET vector (transverse component only)
        metVecT = R.TLorentzVector()
        metVecT.SetPtEtaPhiM(ch.MET_pt, 0, ch.MET_phi, 0)
        
        ## Find AK8 jet(s) passing cuts
        for iFat in range(ch.nFatJet):
            if      ch.FatJet_mass[iFat]  < CUT_MASS: continue
            if ch.FatJet_msoftdrop[iFat]  < CUT_MSOFT: continue
            if        ch.FatJet_pt[iFat]  < 170: continue
            if   abs(ch.FatJet_eta[iFat]) > 2.4: continue
            if     ch.FatJet_jetId[iFat]  <   6: continue
            ## Save 4-vectors of passing AK8 jets
            if VERBOSE: print('    - iFat = %d' % iFat)
            fatVec  = R.TLorentzVector()
            fatVecT = R.TLorentzVector()
            fatVec .SetPtEtaPhiM( ch.FatJet_pt[iFat], ch.FatJet_eta[iFat], ch.FatJet_phi[iFat], ch.FatJet_mass[iFat] )
            fatVecT.SetPtEtaPhiM( ch.FatJet_pt[iFat],                   0, ch.FatJet_phi[iFat], ch.FatJet_mass[iFat] )
            ## Skip jets which overlap selected leptons
            if fatVec.DeltaR(lepVecs[0]) < 0.8: continue
            if VERBOSE: print('    - Passed fatVec.DeltaR(lepVecs[0]) < 0.8 cut')
            ## Skip jets which are too far from selected lepton
            if fatVec.DeltaR(lepVecs[0]) > CUT_DR_LJ and CUT_DR_LJ > 0: continue
            if VERBOSE: print('    - Passed fatVec.DeltaR(lepVecs[0]) > CUT_DR_LJ cut')
            ## Save selected AK8 jet
            fatIdxs .append(iFat)
            fatVecs .append(fatVec)
            fatVecTs.append(fatVecT)
        ## End loop: for iFat in range(ch.nFatJet)

        ## Only keep event if there is at least one selected fat jet
        if len(fatVecs) == 0: continue

        if VERBOSE: print('  * Pass len(fatVecs) == 0 cut')

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
            lJ_vec  = lepVecs [0] + fatVec
            lJ_vecT = lepVecTs[0] + fatVecT

            ## Store MET vector with eta equal to lep+AK8 vector
            metVec = R.TLorentzVector()
            metVec.SetPtEtaPhiM( ch.MET_pt, lJ_vec.Eta(), ch.MET_phi, 0)

            ## Shortcut 4-vector for entire ttbar system
            lv_vec   = lepVecs [0]  + metVec
            lv_vecT  = lepVecTs[0]  + metVecT
            lvJ_vec  = lJ_vec  + metVec
            lvJ_vecT = lJ_vecT + metVecT

            jetIdxs  = []  ## Indices of selected AK4 jets
            jetVecs  = []  ## TLorentzVectors of selected AK4 jets
            jetVecTs = []  ## Transverse component only (eta = 0)
            jlepVecs = []  ## TLorentzVectors of selected AK4 jets overlapping selected lepton
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
            deepB_vec = R.TLorentzVector()
            flavB_vec = R.TLorentzVector()
            near_fat_idx = -99
            near_lJ_idx  = -99
            near_fat_vec = R.TLorentzVector()
            near_lJ_vec  = R.TLorentzVector()
            for iJet in range(ch.nJet):
                if      ch.Jet_pt [iJet]  <  25: continue
                if  abs(ch.Jet_eta[iJet]) > 4.7: continue
                ## Tight jet ID: https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVUL
                ## PhysicsTools/NanoAOD/python/jets_cff.py
                ## Jet ID flags bit1 is loose, bit2 is tight, bit3 is tightLepVeto
                if    ch.Jet_jetId[iJet]  <   6: continue
                ## Loose pileup ID: https://twiki.cern.ch/twiki/bin/view/CMS/PileupJetIDUL
                ## puId==4 means 100: pass loose ID, fail medium, fail tight
                ## PhysicsTools/NanoAOD/python/jets_cff.py
                ## userInt('puId106XUL18Id') : Pileup ID flags with 106X (2018) training
                if      ch.Jet_pt [iJet]  <  50:
                    if ch.Jet_puId[iJet]  <   4: continue
                ## Save 4-vectors of AK4 jets
                jetVec  = R.TLorentzVector()
                jetVecT = R.TLorentzVector()
                jetVec .SetPtEtaPhiM( ch.Jet_pt[iJet], ch.Jet_eta[iJet], ch.Jet_phi[iJet], ch.Jet_mass[iJet] )
                jetVecT.SetPtEtaPhiM( ch.Jet_pt[iJet],                0, ch.Jet_phi[iJet], ch.Jet_mass[iJet] )
                ## Don't allow overlap with muon, electron, or AK8 jet
                if jetVec.DeltaR(fatVec) < 0.8: continue
                if jetVec.DeltaR(lepVecs[0]) < 0.4:
                    jlepVecs.append(jetVec)
                    continue
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
                    if ch.Jet_btagDeepB[iJet] > max_deepB:
                        max_deepB = ch.Jet_btagDeepB[iJet]
                        deepB_vec = jetVec
                    if ch.Jet_btagDeepFlavB[iJet] > max_flavB:
                        max_flavB = ch.Jet_btagDeepFlavB[iJet]
                        flavB_vec = jetVec

                ## Store jets closest to fat jet and lep+fat system
                if near_fat_idx < 0 or jetVec.DeltaR(fatVec) < near_fat_vec.DeltaR(fatVec):
                    near_fat_idx = iJet
                    near_fat_vec = jetVec
                if near_lJ_idx < 0 or jetVec.DeltaR(lJ_vec) < near_lJ_vec.DeltaR(lJ_vec):
                    near_lJ_idx = iJet
                    near_lJ_vec = jetVec

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

            ## No double-medium b-tagged AK4 jet
            if max_deepB > DeepB['M'] and max_flavB > FlavB['M']: continue
            nPassBtag += 1

            ## Lepton pT cuts to match single lepton triggers
            if lepVecs[0].Pt() < 25: continue
            nPassLep += 1

            # ## MET cut to suppress background
            # if metVec.Pt() < 30: continue
            nPassMET += 1

            # ## Full system kinematic selection cuts
            if lvJ_vec.Pt() < CUT_PT_LVJ and CUT_PT_LVJ > 0: continue
            if lvJ_vecT.M() > CUT_MT_LVJ and CUT_MT_LVJ > 0: continue
            if isrs_HT      < CUT_HT_ISR and CUT_HT_ISR > 0: continue
            if isrsVec.Pt() < CUT_PT_ISR and CUT_PT_ISR > 0: continue
            nPassKin += 1

            # ## Tighter selection cuts
            nPassTight += 1

            if isData and ch.run > 319077:  ## HEM veto
                if fatVec.Eta() < -1.17 and fatVec.Phi() > -1.97 and fatVec.Phi() < -0.47:
                    continue
            nPassData += 1

            ## Trigger selection in data and MC
            if LepCat == 'SingleMu':
                if not ( ch.L1_SingleMu22 and (ch.HLT_IsoMu24 or ch.HLT_Mu50)): continue
            if LepCat == 'SingleEG':
                if not ( (ch.L1_SingleIsoEG28er2p5 or ch.L1_SingleEG36er2p5) and ch.HLT_Ele32_WPTight_Gsf): continue
            nPassTrg += 1

            ## AK8 jet has passsed all event selection cuts
            hasCand = True
            ## Signal events must have AK8 jet matching two b-quarks; background events must not
            nBinJ = 0
            nQinJ = 0
            nCinJ = 0
            for xVec in bVecs:
                if fatVec.DeltaR(xVec) < 0.8: nBinJ += 1
            for xVec in qVecs:
                if fatVec.DeltaR(xVec) < 0.8: nQinJ += 1
            for xVec in cVecs:
                if fatVec.DeltaR(xVec) < 0.8: nCinJ += 1
            if       nBinJ < min(CUT_GEN_BBQQ,2): continue
            if nQinJ+nCinJ <    (CUT_GEN_BBQQ-2): continue
            if CUT_GEN_BBQQ < 0 and nBinJ >= 2:   continue
            nPassGen += 1

            ## Remove events where lepton/jet pT fraction is small for dR(lep,jet) < 0.4
            lepjet_frac = 999.
            lepjet_dR   = 999.
            for ii in range(len(jlepVecs)):
                if (lepVecs[0].Pt() / jlepVecs[ii].Pt()) < lepjet_frac:
                    lepjet_frac = lepVecs[0].Pt() / jlepVecs[ii].Pt()
                    lepjet_dR   = lepVecs[0].DeltaR(jlepVecs[ii])

            if CUT_LFRAC > 0 and lepjet_dR < 0.4 and lepjet_frac < 0.6: continue
            nPassLepFrac += 1
            
            # if nBinJ >= 2 and nQinJ+nCinJ >= 1:
            #     print('lepjet_frac = %.3f' % lepjet_frac)
            
            # for jj in range(len(jlepVecs)):
            #     if lepVecs[0].Pt() / jlepVecs[jj].Pt() < 0.7 and nBinJ >= 2 and nQinJ+nCinJ >= 1:
            #         print('\n\n*** Weird case where lepVecs[0] pT / jlepVecs[%d] pT = %.2f  ***' % (jj, lepVecs[0].Pt() / jlepVecs[jj].Pt()))
            #         print('run = %d, LS = %d, event = %d' % (ch.run, ch.luminosityBlock, ch.event))
            #         print('nBinJ = %d, nCinJ = %d, nQinJ = %d' % (nBinJ, nCinJ, nQinJ))
            #         for ii in range(len(muVecs)):
            #             print('muVecs[%d] pT = %d, phi = %.2f, eta = %.2f' % (ii, muVecs[ii].Pt(), muVecs[ii].Phi(), muVecs[ii].Eta()))
            #         for ii in range(len(eleVecs)):
            #             print('eleVecs[%d] pT = %d, phi = %.2f, eta = %.2f' % (ii, eleVecs[ii].Pt(), eleVecs[ii].Phi(), eleVecs[ii].Eta()))
            #         print('lepVecs[0] pT = %d, phi = %.2f, eta = %.2f' % (lepVecs[0].Pt(), lepVecs[0].Phi(), lepVecs[0].Eta()))
            #         for ii in range(len(jlepVecs)):
            #             print('jlepVec[%d] pT = %d, phi = %.2f, eta = %.2f, dR = %.2f, pT frac = %.2f' % (ii, jlepVecs[ii].Pt(), jlepVecs[ii].Phi(), jlepVecs[ii].Eta(), jlepVecs[ii].DeltaR(lepVecs[0]), lepVecs[0].Pt() / jlepVecs[ii].Pt()))

            
            # if (lvJ_vec+isrsVec).Pt() > 200 and lepjet_frac > 0.45:
            #     print('\n\n*** Weird case where pt_lvJ_ISRs = %d  ***' % (lvJ_vec+isrsVec).Pt())
            #     print('run = %d, LS = %d, event = %d' % (ch.run, ch.luminosityBlock, ch.event))
            #     print('nBinJ = %d, nCinJ = %d, nQinJ = %d' % (nBinJ, nCinJ, nQinJ))
            #     print('lvJ_ISRs   pT = %d, phi = %.2f, eta = %.2f' % ((lvJ_vec+isrsVec).Pt(), (lvJ_vec+isrsVec).Phi(), (lvJ_vec+isrsVec).Eta()))
            #     print('lepVecs[0] pT = %d, phi = %.2f, eta = %.2f' % (lepVecs[0].Pt(), lepVecs[0].Phi(), lepVecs[0].Eta()))
            #     print('metVec     pT = %d, phi = %.2f, eta = %.2f' % (metVec.Pt(), metVec.Phi(), metVec.Eta()))
            #     print('fatVec     pT = %d, phi = %.2f, eta = %.2f' % (fatVec.Pt(), fatVec.Phi(), fatVec.Eta()))
            #     print('lvJ_vec    pT = %d, phi = %.2f, eta = %.2f' % (lvJ_vec.Pt(), lvJ_vec.Phi(), lvJ_vec.Eta()))
            #     print('isrsVec    pT = %d, phi = %.2f, eta = %.2f' % (isrsVec.Pt(), isrsVec.Phi(), isrsVec.Eta()))
            #     for ii in range(len(jetVecs)):
            #         print('jetVec[%d] pT = %d, phi = %.2f, eta = %.2f' % (ii, jetVecs[ii].Pt(), jetVecs[ii].Phi(), jetVecs[ii].Eta()))
            #     for ii in range(len(jlepVecs)):
            #         print('jlepVec[%d] pT = %d, phi = %.2f, eta = %.2f, dR = %.2f, pT frac = %.2f' % (ii, jlepVecs[ii].Pt(), jlepVecs[ii].Phi(), jlepVecs[ii].Eta(), jlepVecs[ii].DeltaR(lepVecs[0]), lepVecs[0].Pt() / jlepVecs[ii].Pt()))

                    
            #######################
            ## Get event weights ##
            #######################
            
            WGT         = 1.0  ## Overall event weight
            WGT_NO_LUMI = 1.0  ## Event weight before luminosity scaling

            if not isData:
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

            # WGT_NO_LUMI = WGT  ## Track event weight before cross-section and luminosity scaling
            # if not isData:
            #     # WGT *= EVT_WGT_MUEG.GetXsecPerEvt( in_dir )
            #     WGT *= LUMI
            # if MAX_EVT > 0 and MAX_EVT < nEntries:
            #     WGT *= (1.0*nEntries / MAX_EVT)


            #####################
            ## Fill histograms ##
            #####################

            ## Construct string for BDT input CSV file
            X1 = {}
            ## Numerical representation of GEN level AK8 jet contents
            X2 = nQinJ + nCinJ + 3*nBinJ + 1
            
            ## Objs = ['lep','fat','MET','ISR','ISRs','lJ','lv','lvJ','jet1','jet2','jet3']
            ## Vars = ['pt','eta','phi','dR','dEta','dEtaS','dPhi','mass','MT','HT']

            ## Whole event variables
            X1=FillX(X1,X2,hst, 'nPV'     , min( ch.PV_npvs, 100 ), WGT )
            X1=FillX(X1,X2,hst, 'nPV_good', min( ch.PV_npvsGood, 100), WGT )
            X1=FillX(X1,X2,hst, 'wgt'     , max( 0.001, min( 1.999, WGT_NO_LUMI ) ), 1.0 )
            X1=FillX(X1,X2,hst, 'wgt_wgt' , max( 0.001, min( 1.999, WGT_NO_LUMI ) ), WGT_NO_LUMI )

            X1=FillX(X1,X2,hst, 'nLep', len(muVecs) + len(eleVecs) - muInFat - eleInFat, WGT)
            X1=FillX(X1,X2,hst, 'nMu' , len(muVecs)  -  muInFat, WGT)
            X1=FillX(X1,X2,hst, 'nEle', len(eleVecs) - eleInFat, WGT)
            X1=FillX(X1,X2,hst, 'nFat', len(fatVecs), WGT)
            X1=FillX(X1,X2,hst, 'nJet', len(jetVecs), WGT)
            X1=FillX(X1,X2,hst, 'nMuInFat' ,  muInFat, WGT)
            X1=FillX(X1,X2,hst, 'nEleInFat', eleInFat, WGT)
            hst['nQC_vs_B_0']       .Fill( min(nBinJ, 3),
                                           min(nQinJ + 3*nCinJ, 5), WGT)
            hst['nQC_vs_B_%d' % X2] .Fill( min(nBinJ, 3),
                                           min(nQinJ + 3*nCinJ, 5), WGT)
            hst['nCh_vs_Bh_0']      .Fill( min(ch.FatJet_nBHadrons[fatIdx], 5),
                                           min(ch.FatJet_nCHadrons[fatIdx], 8), WGT)
            hst['nCh_vs_Bh_%d' % X2].Fill( min(ch.FatJet_nBHadrons[fatIdx], 5),
                                           min(ch.FatJet_nCHadrons[fatIdx], 8), WGT)
            hst['lepjet_frac_vs_dR_0']      .Fill( min(lepjet_dR, 0.499), min(lepjet_frac, 1.099), WGT)
            hst['lepjet_frac_vs_dR_%d' % X2].Fill( min(lepjet_dR, 0.499), min(lepjet_frac, 1.099), WGT)

            ## Lepton ID
            if LepCat == 'SingleMu':
                X1=FillX(X1,X2,hst, 'mu_mvaTTH'  , ch.Muon_mvaTTH[xMu], WGT)
                X1=FillX(X1,X2,hst, 'lep_mvaTTH' , ch.Muon_mvaTTH[xMu], WGT)
                X1=FillX(X1,X2,hst, 'mu_miniIso' , min(ch.Muon_miniPFRelIso_all[xMu], 0.21), WGT)
                X1=FillX(X1,X2,hst, 'lep_miniIso', min(ch.Muon_miniPFRelIso_all[xMu], 0.21), WGT)
                X1=FillX(X1,X2,hst, 'mu_SIP'     , min( abs(ch.Muon_sip3d[xMu]), 7.99), WGT)
                X1=FillX(X1,X2,hst, 'lep_SIP'    , min( abs(ch.Muon_sip3d[xMu]), 7.99), WGT)
                X1=FillX(X1,X2,hst, 'mu_ID'      , ch.Muon_mediumPromptId[xMu] + \
                                         2*ch.Muon_tightId     [xMu], WGT)
            if LepCat == 'SingleEG':
                X1=FillX(X1,X2,hst, 'ele_mvaTTH' , ch.Electron_mvaTTH[xEle], WGT)
                X1=FillX(X1,X2,hst, 'lep_mvaTTH' , ch.Electron_mvaTTH[xEle], WGT)
                X1=FillX(X1,X2,hst, 'ele_miniIso', min(ch.Electron_miniPFRelIso_all[xEle], 0.79), WGT)
                X1=FillX(X1,X2,hst, 'lep_miniIso', min(ch.Electron_miniPFRelIso_all[xEle], 0.79), WGT)
                X1=FillX(X1,X2,hst, 'ele_SIP'    , min( abs(ch.Electron_sip3d[xEle]), 7.99), WGT)
                X1=FillX(X1,X2,hst, 'lep_SIP'    , min( abs(ch.Electron_sip3d[xEle]), 7.99), WGT)
                X1=FillX(X1,X2,hst, 'ele_ID'     ,   ch.Electron_mvaFall17V2Iso_WPL   [xEle] + \
                                         2*ch.Electron_mvaFall17V2noIso_WP90[xEle] + \
                                         4*ch.Electron_mvaFall17V2Iso_WP90  [xEle] + \
                                         8*ch.Electron_mvaFall17V2noIso_WP80[xEle] + \
                                        16*ch.Electron_mvaFall17V2Iso_WP80  [xEle], WGT)

            ## pT
            X1=FillX(X1,X2,hst, 'pt_lep', lepVecs[0].Pt(), WGT)
            X1=FillX(X1,X2,hst, 'pt_fat', fatVec.Pt(), WGT)
            X1=FillX(X1,X2,hst, 'pt_MET', metVec.Pt(), WGT)
            X1=FillX(X1,X2,hst, 'pt_lJ' , lJ_vec.Pt(), WGT)
            X1=FillX(X1,X2,hst, 'pt_lv' , lv_vec.Pt(), WGT)
            X1=FillX(X1,X2,hst, 'pt_lvJ', lvJ_vec.Pt(), WGT)
            for i in range(min(len(jetVecs),3)):
                X1=FillX(X1,X2,hst, 'pt_jet%d' % (i+1), jetVecs[i].Pt(), WGT)

            ## eta and phi
            X1=FillX(X1,X2,hst, 'eta_lep', lepVecs[0].Eta(), WGT)
            X1=FillX(X1,X2,hst, 'eta_fat', fatVec.Eta(), WGT)
            X1=FillX(X1,X2,hst, 'eta_lJ' , lJ_vec.Eta(), WGT)
            for i in range(min(len(jetVecs),3)):
                X1=FillX(X1,X2,hst, 'eta_jet%d' % (i+1), jetVecs[i].Eta(), WGT)

            X1=FillX(X1,X2,hst, 'etaA_lep', abs(lepVecs[0].Eta()), WGT)
            X1=FillX(X1,X2,hst, 'etaA_fat', abs(fatVec.Eta()), WGT)
            X1=FillX(X1,X2,hst, 'etaA_lJ' , abs(lJ_vec.Eta()), WGT)
            for i in range(min(len(jetVecs),3)):
                X1=FillX(X1,X2,hst, 'etaA_jet%d' % (i+1), abs(jetVecs[i].Eta()), WGT)

            X1=FillX(X1,X2,hst, 'etaS_lep', lepVecs[0].Eta()*(1 - 2*(lJ_vec.Eta() < 0)), WGT)
            X1=FillX(X1,X2,hst, 'etaS_fat', fatVec.Eta()*(1 - 2*(lJ_vec.Eta() < 0)), WGT)
            for i in range(min(len(jetVecs),3)):
                X1=FillX(X1,X2,hst, 'etaS_jet%d' % (i+1), jetVecs[i].Eta()*(1 - 2*(lJ_vec.Eta() < 0)), WGT)

            X1=FillX(X1,X2,hst, 'phi_lep', lepVecs[0].Phi(), WGT)
            X1=FillX(X1,X2,hst, 'phi_fat' , fatVec.Phi(), WGT)
            X1=FillX(X1,X2,hst, 'phi_MET' , metVec.Phi(), WGT)

            ## lepton vs. enclosing jet
            X1=FillX(X1,X2,hst, 'lepjet_frac', min(lepjet_frac, 1.099), WGT)
            X1=FillX(X1,X2,hst, 'lepjet_dR',   min(lepjet_dR,   0.499), WGT)

            ## mass
            X1=FillX(X1,X2,hst, 'mass_fat', fatVec.M(), WGT)
            X1=FillX(X1,X2,hst, 'mass_lJ' , lJ_vec.M(), WGT)
            X1=FillX(X1,X2,hst, 'mass_lvJ', lvJ_vec.M(), WGT)

            ## MT
            X1=FillX(X1,X2,hst, 'MT_lv' , lv_vecT.M(), WGT)
            X1=FillX(X1,X2,hst, 'MT_lvJ', lvJ_vecT.M(), WGT)

            ## HT
            X1=FillX(X1,X2,hst, 'HT_lv' , lepVecs[0].Pt()+metVec.Pt(), WGT)
            X1=FillX(X1,X2,hst, 'HT_lJ' , lepVecs[0].Pt()+fatVec.Pt(), WGT)
            X1=FillX(X1,X2,hst, 'HT_lvJ', lepVecs[0].Pt()+fatVec.Pt()+metVec.Pt(), WGT)

            ## lep-jet
            X1=FillX(X1,X2,hst, 'dR_lep_fat'  ,  abs( lepVecs[0].DeltaR  (fatVec) ), WGT)
            X1=FillX(X1,X2,hst, 'dPhi_lep_fat',  abs( lepVecs[0].DeltaPhi(fatVec) ), WGT)
            X1=FillX(X1,X2,hst, 'dEta_lep_fat',  abs( lepVecs[0].Eta() -  fatVec.Eta() ), WGT)
            X1=FillX(X1,X2,hst, 'dEtaS_lep_fat',     (lepVecs[0].Eta() -  fatVec.Eta())*(1 - 2*(fatVec.Eta() < 0)), WGT)

            ## dPhi with MET
            X1=FillX(X1,X2,hst, 'dPhi_lep_MET', abs( lepVecs[0].DeltaPhi(metVec) ), WGT)
            X1=FillX(X1,X2,hst, 'dPhi_lJ_MET' , abs( lJ_vec    .DeltaPhi(metVec) ), WGT)
            X1=FillX(X1,X2,hst, 'dPhi_lv_fat' , abs( lv_vec    .DeltaPhi(fatVec) ), WGT)

            ## ISR plots
            X1=FillX(X1,X2,hst, 'pt_ISRs', isrsVec.Pt(), WGT)
            X1=FillX(X1,X2,hst, 'HT_ISRs', isrs_HT, WGT)

            if len(jetVecs) > 0:
                X1=FillX(X1,X2,hst, 'dR_fat_jet_min'  , fatVec.DeltaR(near_fat_vec), WGT)
                X1=FillX(X1,X2,hst, 'dR_lJ_jet_min'   , lJ_vec.DeltaR(near_lJ_vec), WGT)
                X1=FillX(X1,X2,hst, 'dR_fat_deepB_max', fatVec.DeltaR(deepB_vec), WGT)
                X1=FillX(X1,X2,hst, 'dR_lJ_deepB_max' , lJ_vec.DeltaR(deepB_vec), WGT)
                X1=FillX(X1,X2,hst, 'dR_fat_flavB_max', fatVec.DeltaR(flavB_vec), WGT)
                X1=FillX(X1,X2,hst, 'dR_lJ_flavB_max' , lJ_vec.DeltaR(flavB_vec), WGT)

                X1=FillX(X1,X2,hst, 'jet_pt_near_fat', min(299, ch.Jet_pt[near_fat_idx]), WGT)
                X1=FillX(X1,X2,hst, 'jet_pt_near_lJ' , min(299, ch.Jet_pt[near_lJ_idx] ), WGT)
                X1=FillX(X1,X2,hst, 'deepB_near_fat', max(-0.099, min(0.999, ch.Jet_btagDeepB[near_fat_idx]) ), WGT)
                X1=FillX(X1,X2,hst, 'deepB_near_lJ' , max(-0.099, min(0.999, ch.Jet_btagDeepB[near_lJ_idx] ) ), WGT)
                X1=FillX(X1,X2,hst, 'flavB_near_fat', max(-0.099, min(0.999, ch.Jet_btagDeepFlavB[near_fat_idx]) ), WGT)
                X1=FillX(X1,X2,hst, 'flavB_near_lJ' , max(-0.099, min(0.999, ch.Jet_btagDeepFlavB[near_lJ_idx] ) ), WGT)

                X1=FillX(X1,X2,hst, 'dR_lJ_ISR'   ,     lJ_vec .DeltaR  (isrVec), WGT)
                X1=FillX(X1,X2,hst, 'dPhi_lvJ_ISR', abs(lvJ_vec.DeltaPhi(isrVec)), WGT)
                X1=FillX(X1,X2,hst, 'dPhi_lJ_ISR' , abs(lJ_vec .DeltaPhi(isrVec)), WGT)
                X1=FillX(X1,X2,hst, 'pt_lvJ_ISR'  ,    (lvJ_vec+isrVec).Pt(), WGT)
                X1=FillX(X1,X2,hst, 'pt_lJ_ISR'   ,    (lJ_vec +isrVec).Pt(), WGT)
                X1=FillX(X1,X2,hst, 'mass_lvJ_ISR',    (lvJ_vec+isrVec).M(), WGT)
                X1=FillX(X1,X2,hst, 'mass_lJ_ISR' ,    (lJ_vec +isrVec).M(), WGT)
                X1=FillX(X1,X2,hst, 'MT_lvJ_ISR'  ,    (lvJ_vecT+isrVecT).M(), WGT)
                X1=FillX(X1,X2,hst, 'MT_lJ_ISR'   ,    (lJ_vecT +isrVecT).M(), WGT)
                X1=FillX(X1,X2,hst, 'HT_lvJ_ISR'  , lepVecs[0].Pt()+fatVec.Pt()+metVec.Pt()+isrVec.Pt(), WGT)
                X1=FillX(X1,X2,hst, 'HT_lJ_ISR'   , lepVecs[0].Pt()+fatVec.Pt()            +isrVec.Pt(), WGT)

                X1=FillX(X1,X2,hst, 'dR_lJ_ISRs'   ,     lJ_vec .DeltaR  (isrsVec), WGT)
                X1=FillX(X1,X2,hst, 'dPhi_lvJ_ISRs', abs(lvJ_vec.DeltaPhi(isrsVec)), WGT)
                X1=FillX(X1,X2,hst, 'dPhi_lJ_ISRs' , abs(lJ_vec .DeltaPhi(isrsVec)), WGT)
                X1=FillX(X1,X2,hst, 'pt_lvJ_ISRs'  ,    (lvJ_vec+isrsVec).Pt(), WGT)
                X1=FillX(X1,X2,hst, 'pt_lJ_ISRs'   ,    (lJ_vec +isrsVec).Pt(), WGT)
                X1=FillX(X1,X2,hst, 'mass_lvJ_ISRs',    (lvJ_vec+isrsVec).M(), WGT)
                X1=FillX(X1,X2,hst, 'mass_lJ_ISRs' ,    (lJ_vec +isrsVec).M(), WGT)
                X1=FillX(X1,X2,hst, 'MT_lvJ_ISRs'  ,    (lvJ_vecT+isrsVecT).M(), WGT)
                X1=FillX(X1,X2,hst, 'MT_lJ_ISRs'   ,    (lJ_vecT +isrsVecT).M(), WGT)
                X1=FillX(X1,X2,hst, 'HT_lvJ_ISRs'  , lepVecs[0].Pt()+fatVec.Pt()+metVec.Pt()+isrs_HT, WGT)
                X1=FillX(X1,X2,hst, 'HT_lJ_ISRs'   , lepVecs[0].Pt()+fatVec.Pt()            +isrs_HT, WGT)
            ## End if len(jetVecs) > 0:

            X1=FillX(X1,X2,hst, 'deepB_max_jet', max(-0.099, min(0.999, max_deepB) ), WGT)
            X1=FillX(X1,X2,hst, 'flavB_max_jet', max(-0.099, min(0.999, max_flavB) ), WGT)

            X1=FillX(X1,X2,hst, 'msoft_fat'    , ch.FatJet_msoftdrop[fatIdx], WGT)

            X1=FillX(X1,X2,hst, 'bb_DDBvLV2_fat'  , max(-0.099, min(0.999, ch.FatJet_btagDDBvLV2[fatIdx]) ), WGT)
            X1=FillX(X1,X2,hst, 'bb_PNet_Xbb_fat' , max(-0.099, min(0.999, ch.FatJet_particleNetMD_Xbb[fatIdx]
                                                          / max( ch.FatJet_particleNetMD_Xbb[fatIdx] +
                                                                 ch.FatJet_particleNetMD_QCD[fatIdx] , 0.001) ) ), WGT)
            X1=FillX(X1,X2,hst, 'bb_deep_ZHbb_fat', max(-0.099, min(0.999, ch.FatJet_deepTagMD_ZHbbvsQCD[fatIdx]) ), WGT)
            X1=FillX(X1,X2,hst, 'bb_batg_Hbb_fat' , max(-0.099, min(0.999, ch.FatJet_btagHbb[fatIdx]) ), WGT)

            X1=FillX(X1,X2,hst, 'qq_PNet_Xqq_fat' , max(-0.099, min(0.999, ch.FatJet_particleNetMD_Xqq[fatIdx]
                                                          / max( ch.FatJet_particleNetMD_Xqq[fatIdx] +
                                                                 ch.FatJet_particleNetMD_QCD[fatIdx] , 0.001) ) ), WGT)
            X1=FillX(X1,X2,hst, '4q_PNet_H4q_fat'  , max(-0.099, min(0.999, ch.FatJet_particleNet_H4qvsQCD[fatIdx]) ), WGT)
            
            X1=FillX(X1,X2,hst, 'TvsQCD_deepMD_fat' , max(-0.099, min(0.999, ch.FatJet_deepTagMD_TvsQCD[fatIdx]) ), WGT)
            X1=FillX(X1,X2,hst, 'TvsQCD_deep_fat'   , max(-0.099, min(0.999, ch.FatJet_deepTag_TvsQCD[fatIdx]) ), WGT)
            X1=FillX(X1,X2,hst, 'TvsQCD_PNet_fat'   , max(-0.099, min(0.999, ch.FatJet_particleNet_TvsQCD[fatIdx]) ), WGT)
            X1=FillX(X1,X2,hst, 'WvsQCD_deepMD_fat' , max(-0.099, min(0.999, ch.FatJet_deepTagMD_WvsQCD[fatIdx]) ), WGT)
            X1=FillX(X1,X2,hst, 'WvsQCD_deep_fat'   , max(-0.099, min(0.999, ch.FatJet_deepTag_WvsQCD[fatIdx]) ), WGT)
            X1=FillX(X1,X2,hst, 'WvsQCD_PNet_fat'   , max(-0.099, min(0.999, ch.FatJet_particleNet_WvsQCD[fatIdx]) ), WGT)
            
            X1=FillX(X1,X2,hst, 'Hto4b_Haa4b_fat'   , max(-0.009, min(0.9999, ch.FatJet_particleNetMD_Hto4b_Haa4b[fatIdx]) ), WGT)
            X1=FillX(X1,X2,hst, 'Hto4b_Haa3b_fat'   , max(-0.009, min(0.9999, ch.FatJet_particleNetMD_Hto4b_Haa3b[fatIdx]) ), WGT)
            X1=FillX(X1,X2,hst, 'Hto4b_Haa34b_fat'  , max(-0.009, min(0.9999, ch.FatJet_particleNetMD_Hto4b_Haa3b[fatIdx] +
                                                                    ch.FatJet_particleNetMD_Hto4b_Haa4b[fatIdx]) ), WGT)
            X1=FillX(X1,X2,hst, 'Hto4b_binary_fat'  , max(-0.009, min(0.9999, ch.FatJet_particleNetMD_Hto4b_binary_Haa4b[fatIdx]) ), WGT)
            X1=FillX(X1,X2,hst, 'Hto4b_binaryLF_fat', max(-0.009, min(0.9999, ch.FatJet_particleNetMD_Hto4b_binaryLF_Haa4b[fatIdx]) ), WGT)

            X1=FillX(X1,X2,hst, 'nSV_max_jet', nSV_max_jet, WGT)
            X1=FillX(X1,X2,hst, 'nSV_ISR'    , nSV_ISR, WGT)
            X1=FillX(X1,X2,hst, 'nSV_ISRs'   , nSV_ISRs, WGT)

            ## Write out line with BDT input variables to CSV file
            bcsv_line = 'BDT_IN,    %d' % X2
            if PRT_BDT_IN:
                for BDT_var in BDT_VARS:
                    if BDT_var in X1.keys():
                        bcsv_line += X1[BDT_var]
                    else:
                        bcsv_line += (', '+('{:.99f}'.format(-99))[:max(7,len(BDT_var))])
                bdt_file.write('\n%s' % bcsv_line)

        ## End loop: for jFat in range(len(fatVecs))
    ## End loop: for iEvt in range(nEntries):

    print('\nFinished loop over events')
    if PRT_BDT_IN:
        bdt_file.close()
    
    print('\nOut of %d pre-selected events, %d pass b-tag veto, %d lep, %d MET, %d (%d, %d) pass (tight, data) kinematic cuts, %d trigger, %d GEN, %d lep pT frac' % (nPassPre, nPassBtag, nPassLep, nPassMET, nPassKin, nPassTight, nPassData, nPassTrg, nPassGen, nPassLepFrac))

    out_file.cd()

    keys_to_delete = []
    key_list_str = ""
    for key in hst.keys():
        if hst[key].Integral() == 0:
            keys_to_delete.append(key)
            continue
        ## Strip off number at the end to get the histogram name
        key_base = ("'%s', " % key.replace('_'+key.split('_')[-1], ""))
        if not key_base in key_list_str:
            key_list_str += key_base
        hst[key].SetLineWidth(2)
        hst[key].SetLineColor(R.kBlack)
        if int(key.split('_')[-1]) == 1:
            hst[key].SetLineColor(R.kGray+2)
        if abs(int(key.split('_')[-1]) - 3) <= 1:
            hst[key].SetLineColor(R.kViolet + int(key.split('_')[-1]))
        if int(key.split('_')[-1]) == 5:
            hst[key].SetLineColor(R.kBlue)
        if int(key.split('_')[-1]) == 6:
            hst[key].SetLineColor(R.kGreen+1)
        if int(key.split('_')[-1]) == 7:
            hst[key].SetLineColor(R.kOrange+8)
        if int(key.split('_')[-1]) == 8:
            hst[key].SetLineColor(R.kMagenta)
        if int(key.split('_')[-1]) == 9:
            hst[key].SetLineColor(R.kRed)

    for key in keys_to_delete:
        del hst[key]
        
    print('\nSaved histograms to output file %s\n' % out_file_str)

    # if PRT_BDT_IN:
    print('\n\n *** FULL LIST OF POSSIBLE VARIABLES FOR BDT INPUT ***')
    print(key_list_str)

    out_file.Write()
    out_file.Close()


if __name__ == '__main__':
    main()
