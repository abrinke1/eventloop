#! /usr/bin/env python
## Develop cut-based selection for AK8 jets with 3 or 4 b-hadrons
## to calibrate PNet X4b tagger efficiency

import os
import sys
import subprocess
import numpy as np
import ROOT as R

R.gROOT.SetBatch(True)  ## Don't display histograms or canvases when drawn

MAX_EVT = 1000000  ## Maximum number of events to process per MC sample
PRT_EVT = 10000    ## Print every Nth event while processing
DEBUG   = False
JET_PT  = 15     ## Minimum AK4 jet pT
#JET_TAG = 'M'    ## AK4 b-tag threshold (Tight, Medium, Loose)


## Location of postprocessed input files
IN_DIR = '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/PNet_v2_2024_11_22/'
## MC samples to test
SAMPS = ['QCD_bEnr_HT700to1000','QCD_BGen_HT700to1000']
DSETS = {'QCD_bEnr_HT700to1000':'QCD_bEnriched_HT700to1000_TuneCP5_13TeV-madgraph-pythia8',
         'QCD_BGen_HT700to1000':'QCD_HT700to1000_BGenFilter_TuneCP5_13TeV-madgraph-pythia8'}

## 2018 DeepFlavB cuts: https://btv-wiki.docs.cern.ch/ScaleFactors/UL2018/
FlavB = {'L': 0.0490, 'M': 0.2783, 'T': 0.7100}

btagM = [0,0,0,0,0,0,0]
btagT = [0,0,0,0,0,0,0]
bHadr = [0,0,0,0,0,0,0]
bPart = [0,0,0,0,0,0,0]
minPts = []

for samp in SAMPS:

    print('\n\n*** Beginning to look at sample %s ***\n' % samp)

    in_dir_name = IN_DIR+DSETS[samp]+'/r1/'
    in_file_names = [in_dir_name+fn for fn in os.listdir(in_dir_name) if fn.endswith('.root')]
    chains = {}
    chains['Events'] = 0

    ## Combine input files into a single "chain"
    for i in range(len(in_file_names)):
        print('Adding file %s' % in_file_names[i])
        for key in chains.keys():
            if i == 0: chains[key] = R.TChain(key)
            chains[key].Add( in_file_names[i] )

    ## Loop through events, select, make plots
    ch = chains['Events']  ## Shortcut expression
    nEntries = ch.GetEntries()
    print('\nEntering loop over %d events\n' % (nEntries))

    for iEvt in range(nEntries):

        if iEvt > MAX_EVT and MAX_EVT > 0: break
        if (iEvt % PRT_EVT) == 0: print('Looking at event #%d / %d' % (iEvt, nEntries))

        ch.GetEntry(iEvt)

        if DEBUG: print('tagHaa4b_v1 = %.4f, cat_idx = %d, passFilters = %d' % (ch.Haa4b_FatH_tagHaa4b_v1,
                                                                                ch.Haa4b_cat_idx,
                                                                                ch.Haa4b_passFilters))

        if ch.Haa4b_iFatH < 0: continue  ## Require good AK8(H-->4b) candidate, including Xbb > 0.75
        iH = ch.Haa4b_iFatH

        ## Loop over AK4 jets to find jets overlapping candidate AK8 jet
        ijms = []
        ijts = []
        jpts = []
        for ij in range(ch.nJet):
            if ch.Jet_btagDeepFlavB[ij] < FlavB['M']: continue  ## b-tag cut
            if ch.Jet_Haa4b_ovlp_iFat[ij] != iH: continue  ## dR < 0.8 to AK8
            if ch.Jet_Haa4b_presel[ij] < 1: continue  ## AK4 selection not including pT or eta cuts
            if abs(ch.Jet_eta[ij]) > 2.4:   continue  ## b-tagging only valid for |eta| < 2.4
            if ch.Jet_pt[ij] < JET_PT:      continue
            ijms.append(ij)
            jpts.append(ch.Jet_pt[ij])
            if ch.Jet_btagDeepFlavB[ij] > FlavB['T']:
                ijts.append(ij)
        ## Select events with at least 4 medium or 3 tight b-tagged AK4 jets overlapping AK8
        if len(ijms) < 4 and len(ijts) < 3: continue
        if DEBUG: print('Passed initial selection!')

        print('* Event %d, %d(%d) M(T) b-tags, %d b-hadrons, %d b-partons, min pT %.1f' % (iEvt, len(ijms), len(ijts),
                                                                                           ch.FatJet_nBHadrons[iH],
                                                                                           ch.FatJet_nBPartons[iH], min(jpts)))
        btagM[len(ijms)] += 1
        btagT[len(ijts)] += 1
        bHadr[ch.FatJet_nBHadrons[iH]] += 1
        bPart[ch.FatJet_nBPartons[iH]] += 1
        minPts.append(min(jpts))
    ## End loop: for iEvt in range(nEntries)

    print('\n*** Finished looking at sample %s ***\n\n' % samp)

## End loop: for samp in SAMPS

print('\n*** Finished looking at all samples! ***\n\n')

print('\n      Counts:     0,     1,     2,     3,     4,     5,     6')
print(' Med. b-tags: '+', '.join([('%5d' % ii) for (ii) in btagM]))
print('Tight b-tags: '+', '.join([('%5d' % ii) for (ii) in btagT]))
print('GEN bHadrons: '+', '.join([('%5d' % ii) for (ii) in bHadr]))
print('GEN bPartons: '+', '.join([('%5d' % ii) for (ii) in bPart]))
print('\nMinimum jet pT:')
for ipt in range(7):
    minpt = 15+5*ipt
    maxpt = 20+5*ipt
    print('%d <= pT < %d: %5d' % (minpt, maxpt, sum(pt >= minpt and pt < maxpt for pt in minPts)))
print('     pT >= 50: %5d' % sum(pt >= 50 for pt in minPts))

print('\n*** All done!!! ***\n\n')
