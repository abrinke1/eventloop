#! /usr/bin/env python

## Script to plot "resolved" and "boosted" topologies vs. Higgs pT and "a" mass
import os
import sys
import math
import subprocess
import ROOT as R
from array import array

R.gROOT.SetBatch(True)  ## Don't display histograms or canvases when drawn
R.gStyle.SetOptStat(0)  ## Don't display stat boxes

## User configuration
MAX_EVT  = -1
PRT_EVT  = 10000
VERBOSE  = False
IN_DIR   = '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/'
OUT_FILE = 'plots/HtoAA_resolved_boosted_sel2.root'
OUT_DIR  = 'plots/pdf/HtoAA_resolved_boosted_sel2/'
MASSES   = [12, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]
#MASSES   = [55]


def main():

    print('\nInside HtoAA_resolved_boosted')

    print('\nGetting input files %s' % IN_DIR)

    in_file_names = []
    for MASS in MASSES:
        subdir = 'SUSY_GluGluH_01J_HToAATo4B_M-%d_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18NanoAODv9/' % MASS
        for f_name in subprocess.check_output(['ls', IN_DIR+subdir], encoding='UTF-8').splitlines():
            if not '.root' in str(f_name): continue
            in_file_names.append(IN_DIR+subdir+f_name)
            print('Appending file: %s' % in_file_names[-1])

    chains = {}
    chains['Events'] = 0
    for i in range(len(in_file_names)):
        # print('Adding file %s' % in_file_names[i])
        print('Adding file #%d' % (i+1))
        for key in chains.keys():
            if i == 0: chains[key] = R.TChain(key)
            chains[key].Add( in_file_names[i] )
            print('Added TChain %s' % key)

    out_file = R.TFile(OUT_FILE,'recreate')
    if not os.path.exists(OUT_DIR):
        os.makedirs(OUT_DIR)

    hst = {}
    for sel in ['all','sel']:
        hst[sel+'_resolved'] = R.TH2D(sel+'_resolved', sel+'_resolved', 40, 2, 10, 66, -0.5, 65.5)
        hst[sel+'_overlap1'] = R.TH2D(sel+'_overlap1', sel+'_overlap1', 40, 2, 10, 66, -0.5, 65.5)
        hst[sel+'_overlap2'] = R.TH2D(sel+'_overlap2', sel+'_overlap2', 40, 2, 10, 66, -0.5, 65.5)
        hst[sel+'_merged1']  = R.TH2D(sel+'_merged1',  sel+'_merged1',  40, 2, 10, 66, -0.5, 65.5)
        hst[sel+'_merged2']  = R.TH2D(sel+'_merged2',  sel+'_merged2',  40, 2, 10, 66, -0.5, 65.5)
        hst[sel+'_boosted3'] = R.TH2D(sel+'_boosted3', sel+'_boosted3', 40, 2, 10, 66, -0.5, 65.5)
        hst[sel+'_boosted4'] = R.TH2D(sel+'_boosted4', sel+'_boosted4', 40, 2, 10, 66, -0.5, 65.5)
    hst['other'] = R.TH2D('other', 'other', 40, 2, 10, 66, -0.5, 65.5)

    

    ch = chains['Events']
    nEntries = ch.GetEntries()
    for iEvt in range(nEntries):

        ch.GetEntry(iEvt)

        if iEvt > MAX_EVT and MAX_EVT > 0: break
        if (iEvt % PRT_EVT) == 0: print('Looking at event #%d / %d' % (iEvt, nEntries))

        ## First find 'a' bosons
        iMom1 = -99
        iMom2 = -99
        aVecs = []
        for iGen in range(ch.nGenPart):
            pdgID = ch.GenPart_pdgId[iGen]
            if pdgID == 36:
                if   iMom1 < 0: iMom1 = iGen
                elif iMom2 < 0: iMom2 = iGen
                else: break
                ## Get 4-vector for particle 
                iVec = R.TLorentzVector()
                iVec.SetPtEtaPhiM( ch.GenPart_pt  [iGen],
                                   ch.GenPart_eta [iGen],
                                   ch.GenPart_phi [iGen],
                                   ch.GenPart_mass[iGen] )
                aVecs.append(iVec)
            if iMom2 > 0: break

        ## Now set 4-vectors for b-quarks
        bVecs1 = []  ## TLorentzVectors of b-quarks from 1st Higgs to aa decay
        bVecs2 = []  ## TLorentzVectors of b-quarks from 2nd Higgs to aa decay
        for iGen in range(ch.nGenPart):
            pdgID = ch.GenPart_pdgId[iGen]
            if abs(pdgID) != 5: continue
            iMom = ch.GenPart_genPartIdxMother[iGen]
            if iMom < 0: continue
            if ch.GenPart_pdgId[iMom] != 36: continue
            ## Get 4-vector for particle 
            iVec = R.TLorentzVector()
            iVec.SetPtEtaPhiM( ch.GenPart_pt  [iGen],
                               ch.GenPart_eta [iGen],
                               ch.GenPart_phi [iGen],
                               ch.GenPart_mass[iGen] )
            if iMom == iMom1:
                bVecs1.append(iVec)
            if iMom == iMom2:
                bVecs2.append(iVec)
            if len(bVecs1) + len(bVecs2) == 4: break

        hVec   = bVecs1[0]+bVecs1[1]+bVecs2[0]+bVecs2[1]
        ptH    = hVec.Pt()
        massA  = (aVecs[0].M() + aVecs[1].M())/2.0
        dR_bb1 = bVecs1[0].DeltaR(bVecs1[1])
        dR_bb2 = bVecs2[0].DeltaR(bVecs2[1])
        dR_min1 = min(bVecs1[0].DeltaR(bVecs2[0]), bVecs1[0].DeltaR(bVecs2[1]))
        dR_min2 = min(bVecs1[1].DeltaR(bVecs2[0]), bVecs1[1].DeltaR(bVecs2[1]))
        dR_max1 = max(bVecs1[0].DeltaR(bVecs2[0]), bVecs1[0].DeltaR(bVecs2[1]))
        dR_max2 = max(bVecs1[1].DeltaR(bVecs2[0]), bVecs1[1].DeltaR(bVecs2[1]))

        labels = []

        if (dR_bb1 < 0.4) + (dR_bb2 < 0.4) + (dR_min1 < 0.4) + (dR_min2 < 0.4) == 0:
            labels.append('all_resolved')
            if bVecs1[0].Pt() > 30 and bVecs1[1].Pt() >	30 and bVecs2[0].Pt() >	30 and bVecs2[1].Pt() > 30:
                labels.append('sel_resolved')

        if (dR_min1 < 0.4) + (dR_min2 < 0.4) == 1:
            labels.append('all_overlap1')
            if bVecs1[0].Pt() > 30 and bVecs1[1].Pt() >	30 and bVecs2[0].Pt() >	30 and bVecs2[1].Pt() > 30:
                labels.append('sel_overlap1')

        if (dR_min1 < 0.4) + (dR_min2 < 0.4) == 2:
            labels.append('all_overlap2')
            if bVecs1[0].Pt() > 30 and bVecs1[1].Pt() >	30 and bVecs2[0].Pt() >	30 and bVecs2[1].Pt() > 30:
                labels.append('sel_overlap2')

        if (dR_bb1 < 0.4) + (dR_bb2 < 0.4) == 1 and (dR_min1 < 0.4) + (dR_min2 < 0.4) == 0:
            labels.append('all_merged1')
            if bVecs1[0].DeltaR(bVecs1[1]) < 0.4 and (bVecs1[0] + bVecs1[1]).Pt() > 40 and \
               bVecs2[0].Pt() > 30 and bVecs2[1].Pt() >	30:
                labels.append('sel_merged1')
            elif bVecs2[0].DeltaR(bVecs2[1]) < 0.4 and (bVecs2[0] + bVecs2[1]).Pt() > 40 and \
               bVecs1[0].Pt() > 30 and bVecs1[1].Pt() >	30:
                labels.append('sel_merged1')

        if (dR_bb1 < 0.4) + (dR_bb2 < 0.4) == 2 and (dR_min1 < 0.4) + (dR_min2 < 0.4) == 0:
            labels.append('all_merged2')
            if (bVecs1[0] + bVecs1[1]).Pt() > 40 and (bVecs2[0] + bVecs2[1]).Pt() > 40:
                labels.append('sel_merged2')

        if (hVec.DeltaR(bVecs1[0]) < 0.8) + (hVec.DeltaR(bVecs1[1]) < 0.8) + \
           (hVec.DeltaR(bVecs2[0]) < 0.8) + (hVec.DeltaR(bVecs2[1]) < 0.8) == 3:
            labels.append('all_boosted3')
            if hVec.Pt() > 200:
                labels.append('sel_boosted3')
        if (hVec.DeltaR(bVecs1[0]) < 0.8) + (hVec.DeltaR(bVecs1[1]) < 0.8) + \
           (hVec.DeltaR(bVecs2[0]) < 0.8) + (hVec.DeltaR(bVecs2[1]) < 0.8) == 4:
            labels.append('all_boosted4')
            if hVec.Pt() > 200:
                labels.append('sel_boosted4')
        if len(labels) == 0:
            labels.append('other')

        if VERBOSE: print(labels)

        if VERBOSE: print('log2(ptH) = %.2f, massA = %.2f' % (min(max(math.log2(ptH), 0.01), 9.99), massA))
        for label in labels:
            hst[label].Fill( min(max(math.log2(ptH), 2.01), 9.99), massA )

    ## End loop: for iEvt in range(nEntries)

    print('\nFinished loop over events')

    out_file.cd()

    can = R.TCanvas()

    for key in hst.keys():
        can.cd()
        hst[key].Draw('colz')
        can.SaveAs(OUT_DIR+hst[key].GetName()+'.pdf')
        hst[key].Write()
        ## Make selection efficiency plots
        if key.startswith('sel_'):
            all_name = key.replace('sel_','all_')
            eff_name = key.replace('sel_','eff_')
            hst[key].Divide(hst[all_name])
            hst[key].SetName(eff_name)
            hst[key].SetTitle(eff_name)
            hst[key].GetZaxis().SetRangeUser(0,1)
            hst[key].Draw('colz')
            can.SaveAs(OUT_DIR+eff_name+'.pdf')
            hst[key].Write()

    out_file.Write()
    
## End function: def main()


if __name__ == '__main__':
    main()
