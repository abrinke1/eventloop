#! /usr/bin/env python

## Compute trigger acceptance

import os
import sys
import subprocess
import numpy as np
import ROOT as R

R.gROOT.SetBatch(True)  ## Don't display histograms or canvases when drawn

## User configuration
PRT_EVT = 100000  ## Print every Nth event
MAX_EVT = -1      ## Number of events to process


def main():

    print('\nInside HtoAA_trigger\n')

    in_file_names = []
    in_dir = '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/'

    # samp = 'sig'
    # in_dir += 'MC/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-All_TuneCP5_13TeV_madgraph_pythia8/'
    # in_file_names.append(in_dir+'hadd/GluGluH_01J_HToAATo4B_Pt150_M-All_NanoAODv9.root')

    samp = 'data'
    in_dir += 'data/PNet_v1_2023_10_06/Run2018-UL2018_MiniAODv2_GT36-v123/JetHT/'
    # in_file_names.append(in_dir+'skims/Hto4b_0p8/PNet_v1.root')
    in_file_names.append(in_dir+'skims/Hto4b_0p8_new/PNet_v1.root')
    
    print('Appending file: %s' % in_file_names[-1])

    out_dir = 'plots/'
    out_file_str = out_dir+'trigger_%s' % samp
    if MAX_EVT > 0: out_file_str += '_%dk' % (MAX_EVT / 1000)
    out_file_str += '_new.root'
    out_file = R.TFile(out_file_str,'recreate')

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
    bins['pt']    = [ 83, 170, 1000]
    bins['logPt'] = [ 26, 7.4, 10.0]
    bins['eta']   = [200,  -5,    5]
    bins['phi']   = [128, -3.2, 3.2]
    bins['mass']  = [ 60,   0,  300]
    bins['msoft'] = [ 60,   0,  300]

    ## Sets of objects to be plotted
    Objs = ['fat','MET','pMET']
    Vars = ['pt','logPt','eta','phi','mass','msoft']

    ## Book histograms (most will not be used, and will be deleted at the end)
    hst = {}
    ## All combinations of objects and variables
    for trg in ['num','den']:
        for var in Vars:
            for obj in Objs:
                hst['%s_%s_%s' % (var, obj, trg)] = \
                R.TH1D( 'h_%s_%s_%s' % (var, obj, trg),
                        '%s %s %s' % (obj, var, trg),
                        bins[var][0], bins[var][1], bins[var][2] )

        hst['nPV_%s' % trg] = R.TH1D('h_nPV_%s' % trg,      \
                                     '# of PVs (%s)' % trg, \
                                     101, -0.5, 100.5)
        hst['nPV_good_%s' % trg] = R.TH1D('h_nPV_good_%s' % trg,      \
                                          '# of good PVs (%s)' % trg, \
                                          101, -0.5, 100.5)

    ## Set plot graphics
    for key in hst.keys():
        hst[key].SetLineWidth(2)
        if '_num' in key: hst[key].SetLineColor(R.kBlue)
        if '_den' in key: hst[key].SetLineColor(R.kBlack)

    ## Loop through events, select, and plot
    nEntries = chains['Events'].GetEntries()
    print('\nEntering loop over %d events\n' % (nEntries))
    nPass    = 0
    nPassTrg = 0

    ch = chains['Events']  ## Shortcut expression

    for iEvt in range(nEntries):

        ch.GetEntry(iEvt)

        if iEvt > MAX_EVT and MAX_EVT > 0: break
        if (iEvt % PRT_EVT) == 0: print('Looking at event #%d / %d' % (iEvt, nEntries))

        fatVecM = R.TLorentzVector()  ## TLorentzVector of selected AK8 jets with full mass
        fatVecS = R.TLorentzVector()  ## TLorentzVector of selected AK8 jets with soft-drop mass

        ## Find AK8 jet(s) passing cuts and matching HtoAA decay
        hasGoodFat = False
        for iFat in range(ch.nFatJet):
            ## For signal MC, cut on GEN b hadron count
            if samp == 'sig':
                if ch.FatJet_nBHadrons[iFat] < 4: continue
            ## For MC and data, cut on ParticleNet score
            if max(max(ch.FatJet_particleNetMD_Hto4b_Haa4b[iFat] + \
                       ch.FatJet_particleNetMD_Hto4b_Haa3b[iFat],  \
                       ch.FatJet_particleNetMD_Hto4b_binary_Haa4b[iFat]), \
                   ch.FatJet_particleNetMD_Hto4b_binaryLF_Haa4b[iFat]) < 0.8: continue
            ## Basic AK8 jet cuts
            if      ch.FatJet_pt [iFat]  < 170: continue
            if  abs(ch.FatJet_eta[iFat]) > 2.4: continue
            if ch.FatJet_msoftdrop[iFat] <  50: continue
            if ch.FatJet_jetId[iFat]     <   6: continue
            hasGoodFat = True
            ## Save 4-vector
            fatVecM.SetPtEtaPhiM(ch.FatJet_pt[iFat], ch.FatJet_eta[iFat], \
                                 ch.FatJet_phi[iFat], ch.FatJet_mass[iFat] )
            fatVecS.SetPtEtaPhiM(ch.FatJet_pt[iFat], ch.FatJet_eta[iFat], \
                                 ch.FatJet_phi[iFat], ch.FatJet_msoftdrop[iFat] )
            break
        if not hasGoodFat: continue

        ## Remove events with no good collision vertices
        if ch.PV_npvsGood < 1:
            continue
        ## Remove events failing MET noise filters
        if not (ch.Flag_goodVertices and ch.Flag_globalSuperTightHalo2016Filter and \
                ch.Flag_HBHENoiseFilter and ch.Flag_HBHENoiseIsoFilter and \
                ch.Flag_EcalDeadCellTriggerPrimitiveFilter and ch.Flag_BadPFMuonFilter and \
                ch.Flag_BadPFMuonDzFilter and ch.Flag_eeBadScFilter and ch.Flag_ecalBadCalibFilter):
            continue

        nPass += 1
        
        ## Get MET 4-vectors
        METvec  = R.TLorentzVector()
        pMETvec = R.TLorentzVector()
        METvec .SetPtEtaPhiM(ch.MET_pt, 0, ch.MET_phi, 0)
        pMETvec.SetPtEtaPhiM(ch.PuppiMET_pt, 0, ch.PuppiMET_phi, 0)

        ## Evaluate passing triggers
        pass_ggH_HLTa = (ch.HLT_PFJet500 or ch.HLT_AK8PFJet500 or \
                         ch.HLT_AK8PFJet400_TrimMass30 or \
                         ch.HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4)
        pass_ggH_HLTb = (ch.HLT_PFHT1050 or ch.HLT_AK8PFHT800_TrimMass50)
        pass_ggH_HLT  = (pass_ggH_HLTa or pass_ggH_HLTb)
        pass_ggH      = (pass_ggH_HLTa and ch.L1_SingleJet180) or \
                        (pass_ggH_HLTb and (ch.L1_SingleJet180 or ch.L1_HTT360er))

        pass_Vjj_HLTa = (ch.HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5)
        pass_Vjj_HLTb = (ch.HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71)
        pass_Vjj_HLT  = (pass_Vjj_HLTa or pass_Vjj_HLTb)
        pass_Vjj      = (pass_Vjj_HLTa and \
                         (ch.L1_HTT320er or ch.L1_HTT360er or ch.L1_HTT400er or ch.L1_ETT2000 or \
                          ch.L1_HTT320er_QuadJet_70_55_40_40_er2p4 or \
                          ch.L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3)) or \
                        (pass_Vjj_HLTb and (ch.L1_DoubleJet150er2p5 or \
                                            ch.HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71))

        pass_VBF_HLT = (ch.HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1 or \
                        ch.HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2)
        pass_VBF     = pass_VBF_HLT and (ch.L1_SingleJet180 or \
                                         ch.L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5)

        if ((pass_ggH != pass_ggH_HLT) or \
            (pass_Vjj != pass_Vjj_HLT) or \
            (pass_VBF != pass_VBF_HLT)):
            print('\nRun %d, LS %d, event %d:' % (ch.run, ch.luminosityBlock, ch.event))
            print('ggH, Vjj, VBF HLT = %d, %d, %d' % (pass_ggH_HLT, pass_Vjj_HLT, pass_VBF_HLT))
            print('ggH, Vjj, VBF L1T = %d, %d, %d' % (pass_ggH, pass_Vjj, pass_VBF))

        pass_trig = (pass_ggH or pass_Vjj or pass_VBF)
        if pass_trig: nPassTrg += 1
        
        
        #####################
        ## Fill histograms ##
        #####################
        
        ## Objs = ['fat','MET','pMET']
        ## Vars = ['pt','logPt','eta','phi','mass','msoft']

        WGT = 1.0

        for trg in ['num','den']:
            if trg == 'num' and not pass_trig: continue
            
            ## Whole event variables
            hst['nPV_%s' % trg]     .Fill( min( ch.PV_npvs, 100 ), WGT )
            hst['nPV_good_%s' % trg].Fill( min( ch.PV_npvsGood, 100), WGT )

            ## mass
            hst['mass_fat_%s' % trg] .Fill(fatVecM.M(), WGT)
            hst['msoft_fat_%s' % trg].Fill(fatVecS.M(), WGT)

            ## pT
            hst['pt_fat_%s' % trg].Fill( fatVecM.Pt(), WGT)
            hst['logPt_fat_%s' % trg].Fill( np.log2(fatVecM.Pt()), WGT)

            ## eta
            hst['eta_fat_%s' % trg].Fill( fatVecM.Eta(), WGT)

            ## phi
            hst['phi_fat_%s' % trg].Fill( fatVecM.Phi(), WGT)

    ## End loop: for iEvt in range(nEntries):

    print('\nFinished loop over %d events (%d passed)' % (nEntries if MAX_EVT < 0 else min(nEntries, MAX_EVT), nPass))

    out_file.cd()

    keys_to_delete = []
    for key in hst.keys():
        if hst[key].Integral() == 0:
            keys_to_delete.append(key)
            continue
    for key in keys_to_delete:
        del hst[key]
    for key in hst.keys():
        hst[key].Write()
        
    print('\nSaved histograms to output file %s\n' % out_file_str)

    out_file.Write()
    out_file.Close()


if __name__ == '__main__':
    main()
