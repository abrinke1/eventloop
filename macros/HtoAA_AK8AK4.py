#! /usr/bin/env python

## Study HtoAAto4b in events with AK8+AK4 yet

import os
import sys
import subprocess
import numpy as np
import ROOT as R

R.gROOT.SetBatch(True)  ## Don't display histograms or canvases when drawn

## User configuration
PRT_EVT = 10000   ## Print every Nth event
MAX_EVT = -1      ## Number of events to process
LUMI    = 59830   ## 2018 integrated lumi (pb-1), certified "Good"
TRIG    = False   ## Require ggH triggers for gg0l category
PSCALE  = 10      ## Prescale events
DEBUG   = False

CAT   = 'gg0l'     ## gg0l, tt
SAMP  = 'ggHtoAA'  ## JetHT, QCD_BGen, QCD_bEnr, ggHtoAA, ttbar_had, ttHtoAA
MASSA = ['mALo','mAMed','mAHi']
GENB  = ['all0j','440j','330j','all1j','44','44','43','33','XX']  ## dR-matched b-quarks and nBHadrons in AK8 jet
PTCUT = 170
MSCUT = 20
AK4PT = 15
BBTAG = 'PNetXbb'  ## PNetXbb, deepTagXbb, avgXbb
BBCUT = 0.75
DRMIN = 0.8  ## Minimum dR(AK8, AK4)
DRMAX = 1.5  ## Maximum dR(AK8, AK4)
H34BWP = 'WP60'  ## WP for Hto(3)4b
H3BCUT = (0.92 if H34BWP == 'WP80' else (0.960 if H34BWP == 'WP60' else None))
H4BCUT = (0.92 if H34BWP == 'WP80' else (0.975 if H34BWP == 'WP60' else None))


## Open relevant files for a given sample
def load_files(samp):
    chain = R.TChain('Events')
    in_dirs  = []
    in_files = []
    top_dir = '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/PNet_v1_2023_10_06/'
    suff = '/skims/Hto4b_0p8/'
    if samp.startswith('ggHtoAA'):
        in_dirs.append(top_dir+'SUSY_GluGluH_01J_HToAATo4B_Pt150_M-All_TuneCP5_13TeV_madgraph_pythia8')
    if samp.startswith('ttHtoAA'):
        in_dirs.append(top_dir+'SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-All_TuneCP5_13TeV_madgraph_pythia8')
    if samp == 'QCD_BGen':
        for HT in ['300to500','500to700','700to1000','1000to1500','1500to2000','2000toInf']:
            in_dirs.append(top_dir+'QCD_HT%s_BGenFilter_TuneCP5_13TeV-madgraph-pythia8' % HT)
    if samp == 'QCD_bEnr':
        for HT in ['300to500','500to700','700to1000','1000to1500','1500to2000','2000toInf']:
            in_dirs.append(top_dir+'QCD_bEnriched_HT%s_TuneCP5_13TeV-madgraph-pythia8' % HT)
    if samp == 'ttbar_had':
        in_dirs.append(top_dir+'TTToHadronic_TuneCP5_13TeV-powheg-pythia8')
    if samp == 'JetHT':
        in_dirs.append(top_dir.replace('/MC/','/data/')+'Run2018-UL2018_MiniAODv2_GT36-v123/JetHT')

    for in_dir in in_dirs:
        print(in_dir)
        for f_name in subprocess.check_output(['ls', in_dir+suff], encoding='UTF-8').splitlines():
            if not str(f_name).endswith('.root'): continue
            in_files.append(in_dir+suff+f_name)

    for i in range(len(in_files)):
        print('Adding file %s' % in_files[i])
        chain.Add( in_files[i] )

    return chain
## End function: load_files(samp)


## Compute cross section weights for MC, currently just using numbers from Google spreadsheet
## https://docs.google.com/spreadsheets/d/1xDLsr3ikLJxuMPNiSRs79YjTzbN64RetXL3A-tL6-hY
def get_weight(samp, genPtH, mA, LHE_HT):
    WGT = 1.0
    if samp == 'JetHT':
        return WGT
    else:
        WGT = LUMI
    
    if samp.startswith('ggHtoAA'):
        if genPtH < 0.01:
            print('\n*** WEIRD ERROR!!! genPtH = %.6f! Quitting. ***')
            sys.exit()
        WGT *= (2.77077 / 500000)
        WGT *= (0.2 if mA == 'mAMed' else 0.333) ## 5 samples in medium mass(a), 3 in low/high
        ## getHiggsPtRewgtForGGToHToAATo4B from https://github.com/siddhesh86/htoaa/blob/ana_SS/htoaa_CommonTools.py
        a0 = 1.45849
        a1 = -0.00400668
        a2 = 0.00000402577
        a3 = -0.00000000138804
        WGT *= min(max(a0 + a1*genPtH + a2*pow(genPtH,2) + a3*pow(genPtH,3), 0.09), 1.02)
    if samp.startswith('ttHtoAA'):
        WGT *= (0.1445 / 200000)
    if samp == 'ttbar_had':
        WGT *= (380.133 / 331506194)
    if samp == 'QCD_BGen':
        if   LHE_HT <  500: WGT *= (27360 /   14144826)
        elif LHE_HT <  700: WGT *= ( 2991 /    8004808)
        elif LHE_HT < 1000: WGT *= (  731.8 /  4642245)
        elif LHE_HT < 1500: WGT *= (  139.3 /  1537452)
        elif LHE_HT < 2000: WGT *= (   14.74 / 1263157)
        else:               WGT *= (    3.09 / 1300672)
    if samp == 'QCD_bEnr':
        if   LHE_HT <  500: WGT *= (16600 /     11197722)
        elif LHE_HT <  700: WGT *= ( 1503 /      9246898)
        elif LHE_HT < 1000: WGT *= (  297.4 /    1844165)
        elif LHE_HT < 1500: WGT *= (   48.08 /   1330829)
        elif LHE_HT < 2000: WGT *= (    3.951 /  1431254)
        else:               WGT *= (    0.6957 / 1357334)

    return WGT
## End function: get_weight(samp, genPtH, mA, LHE_HT)


def main():

    print('\nInside HtoAA_AK8AK4\n')
    isSig = ('HtoAA' in SAMP)

    ## Dictionary for all histograms
    hst = {}
    clr = {}
    clr['all0j'] = R.kGray+1  if isSig else R.kBlack 
    clr['440j']  = R.kMagenta if isSig else R.kBlack 
    clr['330j']  = R.kGreen   if isSig else R.kBlack 
    clr['all1j'] = R.kBlue    if isSig else R.kBlack 
    clr['44'] = R.kViolet
    clr['43'] = R.kCyan-7
    clr['33'] = R.kGreen+1
    clr['XX'] = R.kRed

    ## Loop over mass(a) and GEN-matched b's to book histograms
    for mA in MASSA:
        hst[mA] = {}
        for gB in GENB:
            hst[mA][gB] = {}

            hst[mA][gB]['AK84_logPt'] = R.TH1D('b', 'AK8+AK4 log_{2}(p_{T})', 52, 7.4, 10.)
            hst[mA][gB]['AK8_logPt']  = R.TH1D('c', 'AK8 log_{2}(p_{T})',     52, 7.4, 10.)
            hst[mA][gB]['AK4_logPt']  = R.TH1D('d', 'AK4 log_{2}(p_{T})',     51, 3.9,  9.)
            hst[mA][gB]['AK84_mass']  = R.TH1D('f', 'AK8+AK4 mass', 60, 0, 300)
            hst[mA][gB]['AK8_mass']   = R.TH1D('g', 'AK8 mass',     60, 0, 300)
            hst[mA][gB]['AK8_msoft']  = R.TH1D('h', 'AK8 soft-drop mass', 40, 0, 200)
            hst[mA][gB]['AK8_massA']  = R.TH1D('i', 'AK8 mass(a)', 140, 0,  70)
            hst[mA][gB]['AK84_dR']    = R.TH1D('k', '#DeltaR(AK8, AK4)', int((DRMAX-DRMIN)*50), DRMIN, DRMAX)
            hst[mA][gB]['AK8_dRb']    = R.TH1D('l', '#DeltaR(AK8, outer GEN b)', int(DRMAX*100),  0.0, DRMAX)
            hst[mA][gB]['AK4_dRb']    = R.TH1D('m', '#DeltaR(AK4, outer GEN b)',  80, 0.0, 0.8)
            hst[mA][gB]['AK8_Xbb']    = R.TH1D('n', 'AK8 %s score' % BBTAG, 100, 0, 1.0)
            hst[mA][gB]['AK4_flavB']  = R.TH1D('o', 'AK4 deepFlavB score',  100, 0, 1.0)
            hst[mA][gB]['AK8_Haa4b']  = R.TH1D('p', 'AK8 Haa4b vs. QCD',    200, 0, 1.0)
            hst[mA][gB]['AK8_Haa34b'] = R.TH1D('q', 'AK8 Haa3b+4b vs. QCD',  40, 0.92, 1.0)

            DRB = int((DRMAX-DRMIN)*20)
            hst[mA][gB]['AK84_logPt_vs_dR84'] = R.TH2D('bb', 'AK8+AK4 log_{2}(p_{T})', 52, 7.4, 10., DRB, DRMIN, DRMAX)
            hst[mA][gB]['AK8_logPt_vs_dR84']  = R.TH2D('cc', 'AK8 log_{2}(p_{T})',     52, 7.4, 10., DRB, DRMIN, DRMAX)
            hst[mA][gB]['AK4_logPt_vs_dR84']  = R.TH2D('dd', 'AK4 log_{2}(p_{T})',     51, 3.9,  9., DRB, DRMIN, DRMAX)
            hst[mA][gB]['AK84_mass_vs_dR84']  = R.TH2D('ff', 'AK8+AK4 mass', 60, 0, 300, DRB, DRMIN, DRMAX)
            hst[mA][gB]['AK8_mass_vs_dR84']   = R.TH2D('gg', 'AK8 mass',     60, 0, 300, DRB, DRMIN, DRMAX)
            hst[mA][gB]['AK8_msoft_vs_dR84']  = R.TH2D('hh', 'AK8 soft-drop mass', 40, 0, 200, DRB, DRMIN, DRMAX)
            hst[mA][gB]['AK8_massA_vs_dR84']  = R.TH2D('ii', 'AK8 mass(a)', 140, 0,  70, DRB, DRMIN, DRMAX)
            hst[mA][gB]['AK8_dRb_vs_dR84']    = R.TH2D('jj', '#DeltaR(AK8, outer GEN b)', int(DRMAX*100),  0.0, DRMAX, DRB, DRMIN, DRMAX)
            hst[mA][gB]['AK4_dRb_vs_dR84']    = R.TH2D('kk', '#DeltaR(AK4, outer GEN b)',  80, 0.0, 0.8, DRB, DRMIN, DRMAX)
            hst[mA][gB]['AK8_Haa4b_vs_dR84']  = R.TH2D('ll', 'AK8 Haa4b vs. QCD',    200, 0, 1.0, DRB, DRMIN, DRMAX)
            hst[mA][gB]['AK8_Haa34b_vs_dR84'] = R.TH2D('mm', 'AK8 Haa3b+4b vs. QCD',  40, 0.92, 1.0, DRB, DRMIN, DRMAX)

            ## Special plots for merged AK8+AK4 candidates with different dR cuts
            for iR in range(10, int(DRMAX*10)):  ## Maximum dR cut x 0.1
            # for iR in range(60,86):  ## Minimum dR cut x 0.01
                hst[mA][gB]['AKr%02d_logPt' % iR] = R.TH1D('a%02d' % iR, 'AK8+AK4(dR < %.1f) log_{2}(p_{T})' % (iR*0.1), 52, 7.4, 10.)
                hst[mA][gB]['AKr%02d_mass'  % iR] = R.TH1D('b%02d' % iR, 'AK8+AK4(dR < %.1f) mass' % (iR*0.1), 60, 0, 300)
                hst[mA][gB]['AKr%02d_logPt_vs_dR84' % iR] = R.TH2D('aa%02d' % iR, 'AK8+AK4(dR < %.1f) log_{2}(p_{T})' % (iR*0.1), 50, 7.4, 10., DRB, DRMIN, DRMAX)
                hst[mA][gB]['AKr%02d_mass_vs_dR84'  % iR] = R.TH2D('bb%02d' % iR, 'AK8+AK4(dR < %.1f) mass' % (iR*0.1), 60, 0, 300, DRB, DRMIN, DRMAX)

            for key in hst[mA][gB].keys():
                hst[mA][gB][key].SetName('h_%s_%s_%s' % (mA, gB, key))
                hst[mA][gB][key].GetXaxis().SetTitle(hst[mA][gB][key].GetTitle())
                if key.endswith('_vs_dR84'):
                    hst[mA][gB][key].GetYaxis().SetTitle('#DeltaR(AK8, AK4)')
                else:
                    hst[mA][gB][key].SetLineWidth(2)
                    hst[mA][gB][key].SetLineColor(clr[gB])
                hst[mA][gB][key].SetDirectory(0)

        ## End loop: for gB in GENB
    ## End loop: for mA in MASSA

    ## --------------------------------------------------------------------
    ## 2018 DeepCSV and DeepFlavB cuts: https://btv-wiki.docs.cern.ch/ScaleFactors/UL2018/
    DeepB = {'L': 0.1208, 'M': 0.4168, 'T': 0.7665}
    FlavB = {'L': 0.0490, 'M': 0.2783, 'T': 0.7100}

    ## Process events, filling histograms
    ch = load_files(SAMP)
    nEntries = ch.GetEntries()
    print('\nFor %s, entering loop over %d events\n' % (SAMP, nEntries))
    nPass = 0

    for iEvt in range(nEntries):
        if iEvt > MAX_EVT and MAX_EVT > 0: break
        if PSCALE > 1 and (iEvt % PSCALE) != 0: continue
        if (iEvt % PRT_EVT) == 0: print('Looking at event #%d / %d' % (iEvt, nEntries))
        ch.GetEntry(iEvt)

        ## Find AK8 jet(s) passing cuts
        xFat = -99
        xH4bVsQCD = -99
        xH34bVsQCD = -99
        fatVec = R.TLorentzVector()
        bbX = {}
        for iFat in range(ch.nFatJet):
            ## Basic AK8 jet cuts
            if ch.FatJet_pt[iFat]        < PTCUT: continue
            if ch.FatJet_msoftdrop[iFat] < MSCUT: continue
            if  abs(ch.FatJet_eta[iFat]) > 2.4: continue
            if ch.FatJet_jetId[iFat]     <   6: continue
            ## Cut on double b-tag score
            bbPN_num = ch.FatJet_particleNetMD_Xbb[iFat]
            bbPN_den = bbPN_num + ch.FatJet_particleNetMD_QCD[iFat]
            ibbX = {}
            ibbX['PNetXbb'] = (0 if bbPN_den <= 0 else (bbPN_num / bbPN_den))
            ibbX['deepTagXbb'] = max(ch.FatJet_deepTagMD_ZHbbvsQCD[iFat], 0.0)
            ibbX['avgXbb'] = (ibbX['PNetXbb'] + ibbX['deepTagXbb']) / 2.0
            if ibbX[BBTAG] < BBCUT: continue
            ## Cut on ParticleNet score (v1)
            PNet_QCD = \
                ch.FatJet_particleNetMD_Hto4b_QCD0b[iFat] + \
                ch.FatJet_particleNetMD_Hto4b_QCD1b[iFat] + \
                ch.FatJet_particleNetMD_Hto4b_QCD2b[iFat] + \
                ch.FatJet_particleNetMD_Hto4b_QCD3b[iFat] + \
                ch.FatJet_particleNetMD_Hto4b_QCD4b[iFat]
            PNet_Hto4b  = ch.FatJet_particleNetMD_Hto4b_Haa4b[iFat] 
            PNet_Hto34b = PNet_Hto4b + ch.FatJet_particleNetMD_Hto4b_Haa3b[iFat] 
            if PNet_Hto4b + PNet_QCD <= 0: continue
            H34bVsQCD = PNet_Hto34b / (PNet_Hto34b + PNet_QCD)
            H4bVsQCD  = PNet_Hto4b / (PNet_Hto4b + PNet_QCD)
            if H34bVsQCD < H3BCUT and H4bVsQCD < H4BCUT: continue
            ## Pick the highest-score AK8 jet
            if H34bVsQCD > xH34bVsQCD or (H34bVsQCD < H3BCUT and H4bVsQCD > xH4bVsQCD):
                xFat = iFat
                xH34bVsQCD = H34bVsQCD
                xH4bVsQCD  = H4bVsQCD
                for key in ibbX.keys():
                    bbX[key] = ibbX[key]
        ## End loop: for iFat in range(ch.nFatJet)
        if xFat < 0: continue
        ## Save 4-vector
        fatVec.SetPtEtaPhiM(ch.FatJet_pt[xFat], ch.FatJet_eta[xFat], \
                            ch.FatJet_phi[xFat], ch.FatJet_mass[xFat] )
        fat_msoft = ch.FatJet_msoftdrop[xFat]

        ## Find AK4 b-tagged jet with dR(AK8, AK4) closest to 0.8 (preferably outside)
        xJet = -99  ## b-tagged AK4 with lowest dR > 0.8
        yJet = -99  ## b-tagged AK4 with highest dR < 0.8
        xdR84 = 99.
        ydR84 = -99.
        xJetVec = R.TLorentzVector()
        yJetVec = R.TLorentzVector()
        nBTags = 0
        for iJet in range(ch.nJet):
            if ch.Jet_btagDeepFlavB[iJet]  < FlavB['M']: continue
            if            ch.Jet_pt[iJet]  < AK4PT: continue
            if       abs(ch.Jet_eta[iJet]) > 2.4: continue
            ## Tight jet ID: https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVUL
            ## PhysicsTools/NanoAOD/python/jets_cff.py
	    ## Jet ID flags bit1 is loose, bit2 is tight, bit3 is tightLepVeto
            if ch.Jet_jetId[iJet] < 6: continue
	    ## Loose pileup ID: https://twiki.cern.ch/twiki/bin/view/CMS/PileupJetIDUL
            ## puId==4 means 100: pass loose ID, fail medium, fail tight
	    ## PhysicsTools/NanoAOD/python/jets_cff.py
            ## userInt('puId106XUL18Id') : Pileup ID flags with 106X (2018) training
            if     ch.Jet_pt  [iJet] <  50:
                if ch.Jet_puId[iJet] <   4: continue
	    ## Set 4-vectors of this AK4 jet
            iJetVec = R.TLorentzVector()
            iJetVec.SetPtEtaPhiM( ch.Jet_pt[iJet], ch.Jet_eta[iJet], ch.Jet_phi[iJet], ch.Jet_mass[iJet] )
            if iJetVec.DeltaR(fatVec) > 0.8 and iJetVec.Pt() > 30:
                nBTags += 1
            if iJetVec.DeltaR(fatVec) < DRMIN: continue
            if iJetVec.DeltaR(fatVec) > DRMAX: continue
            if iJetVec.DeltaR(fatVec) > 0.8 and iJetVec.DeltaR(fatVec) < xdR84:
                xJet    = iJet
                xJetVec = iJetVec
                xdR84   = iJetVec.DeltaR(fatVec)
                ## Count low-pT AK4 b-tags only if they can form AK8+AK4
                if iJetVec.Pt() < 30:
                    nBTags += 1
            if iJetVec.DeltaR(fatVec) < 0.8 and iJetVec.DeltaR(fatVec) > ydR84:
                yJet    = iJet
                yJetVec = iJetVec
                ydR84   = iJetVec.DeltaR(fatVec)
        ## End loop: for iJet in range(ch.nJet)

        ## Select additional AK4 jet with dR(AK8, AK4) > DRMIN (> 0.8 preferred)
        if xJet >= 0:
            jetVec = xJetVec
        elif yJet >= 0 and ydR84 > DRMIN:
            xJet   = yJet
            jetVec = yJetVec
            xdR84  = ydR84
        ## Distinguish between "1j" and "0j" events
        if xJet >= 0:
            if xH34bVsQCD < H3BCUT: continue
            AK84Vec = fatVec+jetVec
        else:  ## If no extra AK4 jet, require Hto4b vs. QCD > H4BCUT
            if xH4bVsQCD < H4BCUT: continue

        ## Require no extra b-tagged AK4 jets for gg0l category
        if CAT == 'gg0l':
            if nBTags > 1*(xJet >= 0 and xdR84 > 0.8): continue


        ## Remove events with no good collision vertices
        if ch.PV_npvsGood < 1:
            continue
        ## Remove events failing MET noise filters
        if not (ch.Flag_goodVertices and ch.Flag_globalSuperTightHalo2016Filter and \
                ch.Flag_HBHENoiseFilter and ch.Flag_HBHENoiseIsoFilter and \
                ch.Flag_EcalDeadCellTriggerPrimitiveFilter and ch.Flag_BadPFMuonFilter and \
                ch.Flag_BadPFMuonDzFilter and ch.Flag_eeBadScFilter and ch.Flag_ecalBadCalibFilter):
            continue
        ## Require relevant trigger for ggH
        if TRIG and cat == 'gg0l' and not \
           ( ( (ch.HLT_PFJet500 or ch.HLT_AK8PFJet500 or ch.HLT_AK8PFJet400_TrimMass30 or \
                ch.HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4) and ch.L1_SingleJet180) or \
             ( (ch.HLT_PFHT1050 or ch.HLT_AK8PFHT800_TrimMass50) and \
               (ch.L1_SingleJet180 or ch.L1_HTT360er) ) or \
             ( (ch.HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71) and \
               (ch.L1_DoubleJet112er2p3_dEta_Max1p6 or ch.L1_DoubleJet150er2p5) ) ): continue

        ## -----------------------------------------------------------------
        ## Done with basic selection cuts, fill WP yield and mass histograms
        nPass += 1
        
        ## Get GEN H, a, and b 4-vectors
        
        if isSig:
            xB    = -99  ## GEN b furthest from AK8 jet
            xdR8b = -99
            xdR4b = -99
            bVec  = R.TLorentzVector()
            bVecs1 = []  ## TLorentzVectors of b-quarks from Higgs to aa decay (a1)
            bVecs2 = []  ## TLorentzVectors of b-quarks from Higgs to aa decay (a2)
            aIdxs  = []  ## Indices of a bosons
            nBinAK8 = 0
            for iGen in range(ch.nGenPart):
                pdgID = ch.GenPart_pdgId[iGen]
                if abs(pdgID) != 5: continue
                iMom = ch.GenPart_genPartIdxMother[iGen]
                if iMom < 0: continue
                if ch.GenPart_pdgId[iMom] != 36: continue
                ## Get 4-vector for particle
                iVec = R.TLorentzVector()
                iVec.SetPtEtaPhiM( ch.GenPart_pt[iGen], ch.GenPart_eta[iGen], ch.GenPart_phi[iGen], ch.GenPart_mass[iGen] )
                if not iMom in aIdxs:
                    aIdxs.append(iMom)
                if iMom == aIdxs[0]:
                    bVecs1.append(iVec)
                else:
                    bVecs2.append(iVec)
                ## Count GEN b's with dR(AK8, b) < 0.8
                if fatVec.DeltaR(iVec) < 0.8:
                    nBinAK8 += 1
                ## Keep GEN b furthest from AK8 jet
                if fatVec.DeltaR(iVec) > xdR8b:
                    xB = iGen
                    bVec = iVec
                    xdR8b = fatVec.DeltaR(iVec)
                    if xJet >= 0:
                        xdR4b = jetVec.DeltaR(iVec)
                if len(bVecs1) + len(bVecs2) == 4: break
            ## End loop: for iGen in range(ch.nGenPart)

            if len(aIdxs) != 2 or len(bVecs1) != 2 or len(bVecs2) != 2:
                print('\n\n*** SUPER-WEIRD ERROR!!! %d a bosons --> %d + %d b-quarks! Quitting. ***' % (len(aIdxs), len(bVecs1), len(bVecs2)))
                sys.exit()
            HVec  = bVecs1[0]+bVecs1[1]+bVecs2[0]+bVecs2[1]
            massGenA = ch.GenPart_mass[aIdxs[0]]
            if abs(HVec.M() - 115) > 15.0 or massGenA < 8 or massGenA > 62.5:
                print('\n\n*** SUPER-WEIRD ERROR!!! mass(H) = %.3f, mass(a) = %.3f! ***' % ((HVec.M(), massGenA)))
        ## End conditional: if isSig


        ## Get mass(a) and set category
        massA = (ch.FatJet_particleNet_massA_Hto4b_v0[xFat] + \
                 ch.FatJet_particleNet_massA_Hto4b_v1[xFat] + \
                 ch.FatJet_particleNet_massA_Hto4b_v3[xFat]) / 3.0
        if isSig:
            mA = ('mALo' if massGenA < 22. else ('mAHi' if massGenA > 47 else 'mAMed'))
        else:
            mA = ('mALo' if massA < 22. else ('mAHi' if massA > 47 else 'mAMed'))

        ## Set GEN b category
        if isSig:
            gBx = str(nBinAK8)+str(min(ch.FatJet_nBHadrons[xFat],4))
        else:
            gBx = str(min(max(ch.FatJet_nCHadrons[xFat], \
                              ch.FatJet_nBHadrons[xFat]),4)) + \
                              str(min(ch.FatJet_nBHadrons[xFat],4))
        if gBx == '34':
            gBx = '43'  ## Combine two rare situations
        if not gBx in GENB:
            gBx = 'XX'

        if xJet >= 0:
            gBs = ['all1j', gBx]
        else:
            gBs = ['all0j']
            if gBx in ['33','44']:
                gBs.append(gBx+'0j')
        
        ## Compute weights
        WGT = get_weight(SAMP, (HVec.Pt() if isSig else 0), mA, 0 if SAMP == 'JetHT' else ch.LHE_HT)
        if PSCALE > 1: WGT *= PSCALE

        for gB in gBs:
            hst[mA][gB]['AK8_logPt'] .Fill( min(np.log2( fatVec.Pt()), 9.99), WGT)
            hst[mA][gB]['AK8_mass']  .Fill( min( fatVec.M(), 299), WGT)
            hst[mA][gB]['AK8_msoft'] .Fill( min( fat_msoft,  199), WGT)
            hst[mA][gB]['AK8_massA'] .Fill( min(     massA, 69.9), WGT)
            hst[mA][gB]['AK8_dRb']   .Fill( min(xdR8b, DRMAX-0.001) if isSig else 0, WGT)
            hst[mA][gB]['AK8_Xbb']   .Fill( bbX[BBTAG], WGT)
            hst[mA][gB]['AK8_Haa4b'] .Fill( xH4bVsQCD,  WGT)
            hst[mA][gB]['AK8_Haa34b'].Fill( xH34bVsQCD, WGT)
            
            if xJet < 0: continue  ## Fill only AK8 histograms if there is no AK4 jet

            ## Combined AK8+AK4 histograms
            hst[mA][gB]['AK84_logPt'].Fill( min(np.log2(AK84Vec.Pt()), 9.99), WGT)
            hst[mA][gB]['AK4_logPt'] .Fill( min(np.log2( jetVec.Pt()), 8.99), WGT)
            hst[mA][gB]['AK84_mass'] .Fill( min(AK84Vec.M(), 299), WGT)
            hst[mA][gB]['AK84_dR']   .Fill( min(xdR84, DRMAX-0.001), WGT)
            hst[mA][gB]['AK4_dRb']   .Fill( min(xdR4b, 0.799) if isSig else 0, WGT)
            hst[mA][gB]['AK4_flavB'] .Fill( ch.Jet_btagDeepFlavB[xJet], WGT)
            
            hst[mA][gB]['AK84_logPt_vs_dR84'].Fill( min(np.log2(AK84Vec.Pt()), 9.99), min(xdR84, DRMAX-0.001), WGT)
            hst[mA][gB]['AK8_logPt_vs_dR84'] .Fill( min(np.log2( fatVec.Pt()), 9.99), min(xdR84, DRMAX-0.001), WGT)
            hst[mA][gB]['AK4_logPt_vs_dR84'] .Fill( min(np.log2( jetVec.Pt()), 8.99), min(xdR84, DRMAX-0.001), WGT)
            hst[mA][gB]['AK84_mass_vs_dR84'] .Fill( min(AK84Vec.M(), 299), min(xdR84, DRMAX-0.001), WGT)
            hst[mA][gB]['AK8_mass_vs_dR84']  .Fill( min( fatVec.M(), 299), min(xdR84, DRMAX-0.001), WGT)
            hst[mA][gB]['AK8_msoft_vs_dR84'] .Fill( min( fat_msoft,  199), min(xdR84, DRMAX-0.001), WGT)
            hst[mA][gB]['AK8_massA_vs_dR84'] .Fill( min(     massA, 69.9), min(xdR84, DRMAX-0.001), WGT)
            hst[mA][gB]['AK8_dRb_vs_dR84']   .Fill( min(xdR8b, DRMAX-0.001) if isSig else 0, min(xdR84, DRMAX-0.001), WGT)
            hst[mA][gB]['AK4_dRb_vs_dR84']   .Fill( min(xdR4b, 0.799) if isSig else 0, min(xdR84, DRMAX-0.001), WGT)
            hst[mA][gB]['AK8_Haa4b_vs_dR84'] .Fill( xH4bVsQCD,  min(xdR84, DRMAX-0.001))
            hst[mA][gB]['AK8_Haa34b_vs_dR84'].Fill( xH34bVsQCD, min(xdR84, DRMAX-0.001), WGT)

            ## Special plots for merged AK8+AK4 candidates with different dR cuts
            for iR in range(10, int(DRMAX*10)):  ## Maximum dR cut x 0.1
            # for iR in range(60,86):  ## Minimum dR cut x 0.01
                if xJet < 0: break
                if xdR84 > iR*0.1: continue  ## For maximum dR, just drop events
                AKrVec = AK84Vec
                # AKrVec = (AK84Vec if xdR84 > 0.01*iR else fatVec)  ## Previous version for minimum dR
                hst[mA][gB]['AKr%d_logPt' % iR].Fill( min(np.log2(AKrVec.Pt()), 9.99), WGT)
                hst[mA][gB]['AKr%d_mass' % iR] .Fill( min(AKrVec.M(), 299), WGT)
                hst[mA][gB]['AKr%d_logPt_vs_dR84' % iR].Fill( min(np.log2(AKrVec.Pt()), 9.99), min(xdR84, DRMAX-0.001), WGT)
                hst[mA][gB]['AKr%d_mass_vs_dR84' % iR] .Fill( min(AKrVec.M(), 299), min(xdR84, DRMAX-0.001), WGT)

            
            
            

    ## End loop: for iEvt in range(nEntries)
    print('Finished with event loop in %s %s: %d passed' % (CAT, SAMP, nPass))
                
    ## End event processing and histogram filling
    ## --------------------------------------------------------------------


    ## --------------------------------------------------------------------
    ## Write out histograms
    out_dir = 'plots/'
    out_file_str = out_dir+'HtoAA_AK8AK4_%s_%s' % (CAT, SAMP)
    if TRIG: out_file_str += '_trig'
    out_file_str += '_dR_%02d_%02d' % (int(DRMIN*10), int(DRMAX*10))
    out_file_str += '_%s' % H34BWP
    out_file_str += '_AK4_pt%d' % AK4PT
    if MAX_EVT > 0: out_file_str += '_%dk' % (MAX_EVT / 1000)
    if PSCALE  > 1: out_file_str += '_ps%d' % PSCALE
    out_file_str += '.root'
    print('\nWriting to %s' % out_file_str)

    out_file = R.TFile(out_file_str,'recreate')
    out_file.cd()

    nHist = 0
    for mA in hst.keys():
        for gB in hst[mA].keys():
            for key in hst[mA][gB].keys():
                if hst[mA][gB][key].Integral() > 0:
                    hst[mA][gB][key].Write()
                    nHist += 1
    print('All done!  Wrote out %d histograms.' % nHist)

    ## End event processing and histogram writing
    ## --------------------------------------------------------------------
    

if __name__ == '__main__':
    main()
