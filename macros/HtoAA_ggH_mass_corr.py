#! /usr/bin/env python

## Compute corrections to Higgs candidate AK8 mass based on balancing jet(s)

import os
import sys
import subprocess
import numpy as np
import ROOT as R

R.gROOT.SetBatch(True)  ## Don't display histograms or canvases when drawn

## User configuration
PRT_EVT = 10000   ## Print every Nth event
MAX_EVT = 1000000 ## Number of events to process


def main():

    print('\nInside HtoAA_ggH_mass_corr\n')

    in_file_names = []
    in_dir = '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/'
    in_dir += 'MC/SUSY_GluGluH_01J_HToAATo4B_Pt150_M-All_TuneCP5_13TeV_madgraph_pythia8/'
    in_file_names.append(in_dir+'hadd/GluGluH_01J_HToAATo4B_Pt150_M-All_NanoAODv9.root')
    print('Appending file: %s' % in_file_names[-1])

    out_dir = 'plots/'
    out_file_str = out_dir+'ggH_mass_corr'
    if MAX_EVT > 0: out_file_str += '_%dk' % (MAX_EVT / 1000)
    out_file_str += '.root'
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
    bins['pt']   = [200,  0, 1000]
    bins['eta']  = [200, -5,    5]
    bins['phi']  = [128, -3.2, 3.2]
    bins['dR']   = [100,  0,   10]
    bins['dEta'] = [100,  0,   10]
    bins['dPhi'] = [ 32,  0,  3.2]
    bins['mass'] = [ 30,  50, 200]
    bins['msoft'] = [30,  50, 200]

    ## Sets of objects to be plotted
    Objs = ['fat','jets','jet1','jet2','jet3','MET','pMET']
    Vars = ['pt','eta','phi','dR','dEta','dPhi','mass','msoft']

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

    hst['nJets']    = R.TH1D('h_nJets',    'Number of AK4 jets outside AK8', 11, -0.5, 10.5)
    hst['nJetsCen'] = R.TH1D('h_nJetsCen', 'Number of AK4 jets outside AK8, |eta| < 2.4', 11, -0.5, 10.5)
    hst['nJetsFwd'] = R.TH1D('h_nJetsFwd', 'Number of AK4 jets outside AK8, |eta| > 2.4', 11, -0.5, 10.5)
    hst['nPV']      = R.TH1D('h_nPV',      '# of PVs',      101, -0.5, 100.5)
    hst['nPV_good'] = R.TH1D('h_nPV_good', '# of good PVs', 101, -0.5, 100.5)

    ## Book corrected mass histograms
    for mass in ['mass', 'msoft']:
        for cor in ['cor1', 'cors', 'corM', 'corP']:
            hst[cor] = R.TH1D('h_%s' % cor, 'Correction %s' % cor, 100, 0.9, 1.9)
            for lim in ['', '_1p1', '_1p2', '_1p3']:
                hst['%s_fat_%s%s' % (mass, cor, lim)] = R.TH1D( 'h_%s_fat_%s%s' % (mass, cor, lim),
                                                                '%s %s %s' % (mass, cor, lim), 30, 50, 200)

    ## Loop through events, select, and plot
    nEntries = chains['Events'].GetEntries()
    print('\nEntering loop over %d events\n' % (nEntries))
    nPass = 0

    ch = chains['Events']  ## Shortcut expression

    for iEvt in range(nEntries):

        ch.GetEntry(iEvt)

        if iEvt > MAX_EVT and MAX_EVT > 0: break
        if (iEvt % PRT_EVT) == 0: print('Looking at event #%d / %d' % (iEvt, nEntries))

        ## Remove events with no good collision vertices
        if ch.PV_npvsGood < 1:
            continue
        
        fatVecM = R.TLorentzVector()  ## TLorentzVector of selected AK8 jets with full mass
        fatVecS = R.TLorentzVector()  ## TLorentzVector of selected AK8 jets with soft-drop mass

        ## Find AK8 jet(s) passing cuts and matching HtoAA decay
        hasGoodFat = False
        for iFat in range(ch.nFatJet):
            if ch.FatJet_nBHadrons[iFat] <   4: continue
            if      ch.FatJet_pt [iFat]  < 170: continue
            if  abs(ch.FatJet_eta[iFat]) > 2.4: continue
            if ch.FatJet_msoftdrop[iFat] <  50: continue
            if ch.FatJet_jetId[iFat]     <   6: continue
            hasGoodFat = True
            ## Save 4-vector
            fatVecM.SetPtEtaPhiM(ch.FatJet_pt[iFat], ch.FatJet_eta[iFat], ch.FatJet_phi[iFat], ch.FatJet_mass[iFat] )
            fatVecS.SetPtEtaPhiM(ch.FatJet_pt[iFat], ch.FatJet_eta[iFat], ch.FatJet_phi[iFat], ch.FatJet_msoftdrop[iFat] )
            break
        if not hasGoodFat: continue
        nPass += 1
        
        jetVecs = []  ## TLorentzVectors of selected AK4 jets
        nJetsCen = 0
        nJetsFwd = 0
        ## Find AK4 jet(s) passing cuts and not overlapping AK8 jet
        for iJet in range(ch.nJet):
            if  abs(ch.Jet_eta[iJet]) > 4.7: continue
            if      ch.Jet_pt [iJet]  <  25: continue
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
            jetVec = R.TLorentzVector()
            jetVec.SetPtEtaPhiM( ch.Jet_pt[iJet], ch.Jet_eta[iJet], ch.Jet_phi[iJet], ch.Jet_mass[iJet] )
            ## Don't allow overlap with AK8 jet
            if jetVec.DeltaR(fatVecM) < 0.8: continue
            if abs(jetVec.Eta()) < 2.4:
                nJetsCen += 1
            else:
                nJetsFwd += 1
            jetVecs.append(jetVec)
        ## End loop: for iJet in range(ch.nJet)

        ## Sum up all AK4 4-vectors
        jetsVec = R.TLorentzVector()
        if len(jetVecs) > 0:
            jetsVec = jetVecs[0]
        for ii in range(1, len(jetVecs)):
            jetsVec = jetsVec + jetVecs[ii]

        ## Get MET 4-vectors
        METvec  = R.TLorentzVector()
        pMETvec = R.TLorentzVector()
        METvec .SetPtEtaPhiM(ch.MET_pt, 0, ch.MET_phi, 0)
        pMETvec.SetPtEtaPhiM(ch.PuppiMET_pt, 0, ch.PuppiMET_phi, 0)

        
        #####################
        ## Fill histograms ##
        #####################
        
        ## Objs = ['fat','jets','jet1','jet2','jet3','MET','pMET']
        ## Vars = ['pt','eta','phi','dR','dEta','dPhi','mass','msoft']

        WGT = 1.0
        
        ## Whole event variables
        hst['nPV']     .Fill( min( ch.PV_npvs, 100 ), WGT )
        hst['nPV_good'].Fill( min( ch.PV_npvsGood, 100), WGT )

        hst['nJets']   .Fill( len(jetVecs), WGT)
        hst['nJetsCen'].Fill( nJetsCen,     WGT)
        hst['nJetsFwd'].Fill( nJetsFwd,     WGT)

        ## mass with scalings
        hst['mass_fat'] .Fill(fatVecM.M(), WGT)
        hst['msoft_fat'].Fill(fatVecS.M(), WGT)

        jet1_scale = 1.0
        jets_scale = 1.0
        MET_scale  = 1.0
        pMET_scale = 1.0
        if len(jetVecs) > 0:
            if abs(fatVecM.DeltaPhi(jetVecs[0])) > 2.75:
                jet1_scale = max(1.0, jetVecs[0].Pt() / fatVecM.Pt())
            if abs(fatVecM.DeltaPhi(jetsVec)) > 2.75:
                jets_scale = max(1.0, jetsVec.Pt() / fatVecM.Pt())

        dPhi_MET  = abs(fatVecM.DeltaPhi(METvec))
        dPhi_pMET = abs(fatVecM.DeltaPhi(pMETvec))
        if abs(dPhi_MET < 0.8):
            MET_scale = 1 + ( (METvec.Pt() / fatVecM.Pt()) * np.cos(dPhi_MET) )
        if abs(dPhi_pMET < 0.8):
            pMET_scale = 1 + ( (pMETvec.Pt() / fatVecM.Pt()) * np.cos(dPhi_pMET) )

        hst['cor1'].Fill(min(jet1_scale, 1.899), WGT)
        hst['cors'].Fill(min(jets_scale, 1.899), WGT)
        hst['corM'].Fill(min(MET_scale, 1.899),  WGT)
        hst['corP'].Fill(min(pMET_scale, 1.899), WGT)

        hst['mass_fat_cor1'] .Fill(fatVecM.M()*jet1_scale, WGT)
        hst['msoft_fat_cor1'].Fill(fatVecS.M()*jet1_scale, WGT)
        hst['mass_fat_cors'] .Fill(fatVecM.M()*jets_scale, WGT)
        hst['msoft_fat_cors'].Fill(fatVecS.M()*jets_scale, WGT)
        hst['mass_fat_corM'] .Fill(fatVecM.M()*MET_scale,  WGT)
        hst['msoft_fat_corM'].Fill(fatVecS.M()*MET_scale,  WGT)
        hst['mass_fat_corP'] .Fill(fatVecM.M()*pMET_scale, WGT)
        hst['msoft_fat_corP'].Fill(fatVecS.M()*pMET_scale, WGT)

        hst['mass_fat_cor1_1p1'] .Fill(fatVecM.M()*min(1.1, jet1_scale), WGT)
        hst['msoft_fat_cor1_1p1'].Fill(fatVecS.M()*min(1.1, jet1_scale), WGT)
        hst['mass_fat_cors_1p1'] .Fill(fatVecM.M()*min(1.1, jets_scale), WGT)
        hst['msoft_fat_cors_1p1'].Fill(fatVecS.M()*min(1.1, jets_scale), WGT)
        hst['mass_fat_corM_1p1'] .Fill(fatVecM.M()*min(1.1, MET_scale ), WGT)
        hst['msoft_fat_corM_1p1'].Fill(fatVecS.M()*min(1.1, MET_scale ), WGT)
        hst['mass_fat_corP_1p1'] .Fill(fatVecM.M()*min(1.1, pMET_scale), WGT)
        hst['msoft_fat_corP_1p1'].Fill(fatVecS.M()*min(1.1, pMET_scale), WGT)

        hst['mass_fat_cor1_1p2'] .Fill(fatVecM.M()*min(1.2, jet1_scale), WGT)
        hst['msoft_fat_cor1_1p2'].Fill(fatVecS.M()*min(1.2, jet1_scale), WGT)
        hst['mass_fat_cors_1p2'] .Fill(fatVecM.M()*min(1.2, jets_scale), WGT)
        hst['msoft_fat_cors_1p2'].Fill(fatVecS.M()*min(1.2, jets_scale), WGT)
        hst['mass_fat_corM_1p2'] .Fill(fatVecM.M()*min(1.2, MET_scale ), WGT)
        hst['msoft_fat_corM_1p2'].Fill(fatVecS.M()*min(1.2, MET_scale ), WGT)
        hst['mass_fat_corP_1p2'] .Fill(fatVecM.M()*min(1.2, pMET_scale), WGT)
        hst['msoft_fat_corP_1p2'].Fill(fatVecS.M()*min(1.2, pMET_scale), WGT)

        hst['mass_fat_cor1_1p3'] .Fill(fatVecM.M()*min(1.3, jet1_scale), WGT)
        hst['msoft_fat_cor1_1p3'].Fill(fatVecS.M()*min(1.3, jet1_scale), WGT)
        hst['mass_fat_cors_1p3'] .Fill(fatVecM.M()*min(1.3, jets_scale), WGT)
        hst['msoft_fat_cors_1p3'].Fill(fatVecS.M()*min(1.3, jets_scale), WGT)
        hst['mass_fat_corM_1p3'] .Fill(fatVecM.M()*min(1.3, MET_scale ), WGT)
        hst['msoft_fat_corM_1p3'].Fill(fatVecS.M()*min(1.3, MET_scale ), WGT)
        hst['mass_fat_corP_1p3'] .Fill(fatVecM.M()*min(1.3, pMET_scale), WGT)
        hst['msoft_fat_corP_1p3'].Fill(fatVecS.M()*min(1.3, pMET_scale), WGT)

        ## pT
        hst['pt_fat'] .Fill( fatVecM.Pt(), WGT)
        if len(jetVecs) > 0:
            hst['pt_jets'].Fill( jetsVec.Pt(), WGT)
            for i in range(min(len(jetVecs),3)):
                hst['pt_jet%d' % (i+1)].Fill( jetVecs[i].Pt(), WGT)

        ## eta
        hst['eta_fat'] .Fill( fatVecM.Eta(), WGT)
        if len(jetVecs) > 0:
            hst['eta_jets'].Fill( jetsVec.Eta(), WGT)
            for i in range(min(len(jetVecs),3)):
                hst['eta_jet%d' % (i+1)].Fill( jetVecs[i].Eta(), WGT)

        ## dEta
        if len(jetVecs) > 0:
            hst['dEta_fat_jets'].Fill( abs(fatVecM.Eta() - jetsVec.Eta()), WGT)
            for i in range(min(len(jetVecs),3)):
                hst['dEta_fat_jet%d' % (i+1)].Fill( abs(fatVecM.Eta() - jetVecs[i].Eta()), WGT)

        ## dPhi
        if len(jetVecs) > 0:
            hst['dPhi_fat_jets'].Fill( abs(fatVecM.DeltaPhi(jetsVec)), WGT)
            for i in range(min(len(jetVecs),3)):
                hst['dPhi_fat_jet%d' % (i+1)].Fill( abs(fatVecM.DeltaPhi(jetVecs[i])), WGT)

        ## dR
        if len(jetVecs) > 0:
            hst['dR_fat_jets'].Fill( abs(fatVecM.DeltaR(jetsVec)), WGT)
            for i in range(min(len(jetVecs),3)):
                hst['dPhi_fat_jet%d' % (i+1)].Fill( abs(fatVecM.DeltaR(jetVecs[i])), WGT)

    ## End loop: for iEvt in range(nEntries):

    print('\nFinished loop over %d events (%d passed)' % (nEntries if MAX_EVT < 0 else min(nEntries, MAX_EVT), nPass))

    out_file.cd()

    keys_to_delete = []
    for key in hst.keys():
        if hst[key].Integral() == 0:
            keys_to_delete.append(key)
            continue
        hst[key].SetLineWidth(2)
        hst[key].SetLineColor(R.kBlack)
        if key.endswith('cor1'):
            hst[key].SetLineColor(R.kViolet)
        if key.endswith('cors'):
            hst[key].SetLineColor(R.kBlue)
        if key.endswith('corM'):
            hst[key].SetLineColor(R.kTeal-1)
        if key.endswith('corP'):
            hst[key].SetLineColor(R.kMagenta+2)
        if key.endswith('1p1'):
            hst[key].SetLineColor(R.kMagenta)
        if key.endswith('1p2'):
            hst[key].SetLineColor(R.kRed)
        if key.endswith('1p3'):
            hst[key].SetLineColor(R.kGreen)

    for key in keys_to_delete:
        del hst[key]
        
    print('\nSaved histograms to output file %s\n' % out_file_str)

    out_file.Write()
    out_file.Close()


if __name__ == '__main__':
    main()
