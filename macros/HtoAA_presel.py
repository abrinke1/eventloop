#! /usr/bin/env python

## Optimize HtoAAto4b pre-selection acceptance
## Run on "skimmed" datasets
## Use ggH signal vs. JetHT data (or QCD MC?) and
##     ttH signal vs. ttbar MC
## (For data, would need to compute TF from sideband to signal)
## Include trigger selection?
## Separate into 3 mass(a) ranges
## Plot 1D mass of signal and background, sum S/sqrt(B)
## Test mass algos: mass, msoft, PNet, scaled mass, scaled msoft
## Test double b-tag cuts: PNet Xbb, deepTag Xbb, avg.
## Test Higgs pT cuts
## Test msoft / mass cuts
## Store histograms in giant dictionary, compute S/sqrt(B)

import os
import sys
import subprocess
import numpy as np
import ROOT as R

R.gROOT.SetBatch(True)  ## Don't display histograms or canvases when drawn

## User configuration
PRT_EVT = 10000   ## Print every Nth event
MAX_EVT = -1     ## Number of events to process
DO_EVT  = False   ## Re-generate histograms from event loop
LUMI    = 59830   ## 2018 integrated lumi (pb-1), certified "Good"
TRIG    = False    ## Require ggH triggers for gg0l category
PSCALE  = -1      ## Prescale events
DEBUG   = False

#CATS  = ['gg0l','Zvv','Wlv','tt']
CATS  = ['Wlv']
SAMPS = {}
#SAMPS['gg0l'] = ['JetHT','QCD_BGen','QCD_bEnr','ggHtoAA','ggHtoAAGen']
#SAMPS['gg0l'] = ['ggHtoAA','ggHtoAAGen3','ggHtoAAGen4']
#SAMPS['gg0l'] = ['QCD_BGen','QCD_bEnr']
#SAMPS['gg0l'] = ['QCD_bEnr']
#SAMPS['gg0l'] = ['JetHT']
#SAMPS['tt']   = ['ttbar_had','ttHtoAA','ttHtoAAGen3','ttHtoAAGen4']
#SAMPS['Zvv']  = ['ZHtoAA','ZHtoAAGen3','ZHtoAAGen4']
#SAMPS['Zvv']  = ['QCD_BGen','QCD_bEnr']
#SAMPS['Zvv']  = ['ttbar_2l']
#SAMPS['Zvv']  = ['MET']
SAMPS['Wlv']  = ['WHtoAA','WHtoAAGen3','WHtoAAGen4']
#SAMPS['Wlv']  = ['WtoLNu']
#SAMPS['Wlv']  = ['EGamma']
MASSA = ['mALo','mAMed','mAHi']
#MASSH = ['logPt','mass','msoft','mPNet','massPt','msoftPt']
MASSH = ['logPt','mass','mass1j','mass1j1p6','mass1j1p5','mass1j1p4','mass1j1p3','mass1j1p2','mass1j1p1']
#MASSH = ['logPt','mass','mass1j','mass1j20','mass1j25','mass1j30']
PTCUT = [170,200,225,250]
#MRCUT = [0.0,0.4,0.5,0.6,0.7,0.8,0.9]
#MSCUT = [20,30,40,50,60,70]
MSCUT = [20,40,50]
#BBTAG = ['PNetXbb','deepTagXbb','avgXbb']
#BBCUT = [0.0,0.6,0.7,0.75,0.8,0.85,0.9]
BBTAG = ['PNetXbb']
BBCUT = [0.75]
PTMAX = -1
PTMIN = -1
#PTCUT = [400]
AK4PT = 15
AK4DR = 1.7

# # # out_file_str += '_msoft_scan'
# # # if DEBUG:
# # SAMPS['gg0l'] = ['JetHT','ggHtoAA']
# # SAMPS['gg0l'] = ['QCD_BGen','QCD_bEnr','ggHtoAA']
# SAMPS['Zvv'] = ['MET','ZHtoAA']
# SAMPS['Zvv'] = ['QCD_BGen','QCD_bEnr','ZHtoAA']
# SAMPS['Zvv'] = ['ttbar_had','ttbar_1l','ttbar_2l','ZHtoAA']
SAMPS['Wlv'] = ['WtoLNu','ttbar_1l','ttbar_2l','WHtoAA']
# SAMPS['Wlv'] = ['SingleMuon','EGamma','WHtoAA']
MASSA = ['mAMedHi']
# # MASSH = ['mass','mass1j']
PTCUT = [170]
MSCUT = [20]
# # BBTAG = ['PNetXbb']
# # BBCUT = [0.75]

## Open relevant files for a given sample
def load_files(samp):
    ## Don't need to separately process signal requiring 4 GEN b's in AK8
    if 'HtoAA' in samp and 'Gen' in samp:
        return None
    chain = R.TChain('Events')
    in_dirs  = []
    in_files = []
    top_dir = '/eos/cms/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/PNet_v1_2023_10_06/'
    suff = '/skims/Hto4b_0p8/'
    genMA = 'All' if MAX_EVT < 0 else '45'
    if samp.startswith('ggHtoAA'):
        in_dirs.append(top_dir+'SUSY_GluGluH_01J_HToAATo4B_Pt150_M-%s_TuneCP5_13TeV_madgraph_pythia8' % genMA)
    if samp.startswith('WHtoAA'):
        in_dirs.append(top_dir+'SUSY_WH_WToAll_HToAATo4B_Pt150_M-%s_TuneCP5_13TeV_madgraph_pythia8' % genMA)
    if samp.startswith('ZHtoAA'):
        in_dirs.append(top_dir+'SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-%s_TuneCP5_13TeV_madgraph_pythia8' % genMA)
    if samp.startswith('ttHtoAA'):
        in_dirs.append(top_dir+'SUSY_TTH_TTToAll_HToAATo4B_Pt150_M-%s_TuneCP5_13TeV_madgraph_pythia8' % genMA)
    if samp == 'QCD_BGen':
        for HT in ['300to500','500to700','700to1000','1000to1500','1500to2000','2000toInf']:
            in_dirs.append(top_dir+'QCD_HT%s_BGenFilter_TuneCP5_13TeV-madgraph-pythia8' % HT)
    if samp == 'QCD_bEnr':
        for HT in ['300to500','500to700','700to1000','1000to1500','1500to2000','2000toInf']:
            in_dirs.append(top_dir+'QCD_bEnriched_HT%s_TuneCP5_13TeV-madgraph-pythia8' % HT)
    if samp == 'ttbar_had':
        in_dirs.append(top_dir+'TTToHadronic_TuneCP5_13TeV-powheg-pythia8')
    if samp == 'ttbar_1l':
        in_dirs.append(top_dir+'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8')
    if samp == 'ttbar_2l':
        in_dirs.append(top_dir+'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8')
    if samp == 'WtoLNu':
        for HT in ['70To100','100To200','200To400','400To600','600To800','800To1200','1200To2500','2500ToInf']:
            in_dirs.append(top_dir+'WJetsToLNu_HT-%s_TuneCP5_13TeV-madgraphMLM-pythia8' % HT)
    if samp in ['JetHT','MET','SingleMuon','EGamma']:
        in_dirs.append(top_dir.replace('/MC/','/data/')+'Run2018-UL2018_MiniAODv2_GT36-v123/'+samp)

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
    elif samp.startswith('WHtoAA'):
        WGT *= (0.179256 / 100000)
    elif samp.startswith('ZHtoAA'):
        WGT *= (0.11352 / 100000)
    elif samp.startswith('ttHtoAA'):
        WGT *= (0.1445 / 200000)
    elif samp == 'ttbar_had':
        WGT *= (380.133 / 331506194)
    elif samp == 'ttbar_1l':
        WGT *= (364.328 / 472557630)
    elif samp == 'ttbar_2l':
        WGT *= (87.339 / 143848848)
    elif samp == 'QCD_BGen':
        if   LHE_HT <  500: WGT *= (27360 /   14144826)
        elif LHE_HT <  700: WGT *= ( 2991 /    8004808)
        elif LHE_HT < 1000: WGT *= (  731.8 /  4642245)
        elif LHE_HT < 1500: WGT *= (  139.3 /  1537452)
        elif LHE_HT < 2000: WGT *= (   14.74 / 1263157)
        else:               WGT *= (    3.09 / 1300672)
    elif samp == 'QCD_bEnr':
        if   LHE_HT <  500: WGT *= (16600 /     11197722)
        elif LHE_HT <  700: WGT *= ( 1503 /      9246898)
        elif LHE_HT < 1000: WGT *= (  297.4 /    1844165)
        elif LHE_HT < 1500: WGT *= (   48.08 /   1330829)
        elif LHE_HT < 2000: WGT *= (    3.951 /  1431254)
        else:               WGT *= (    0.6957 / 1357334)
    elif samp == 'WtoLNu':
        if LHE_HT < 100:
            WGT *= (1440 /      66220256)
        elif LHE_HT < 200:
            WGT *= (1431 /      51408967)
        elif LHE_HT < 400:
            WGT *= ( 382.1 /    58225632)
        elif LHE_HT < 600:
            WGT *= (  51.54 /    7444030)
        elif LHE_HT < 800:
            WGT *= (  12.49 /    7718765)
        elif LHE_HT < 1200:
            WGT *= (   5.619 /   7306187)
        elif LHE_HT < 2500:
            WGT *= (   1.321 /   6481518)
        else:
            WGT *= (   0.02992 / 2097648)
    else:
        print('\n\n*** FATAL ERROR! %s has no get_weight! Quitting.' % samp)
        sys.exit()

    return WGT
## End function: get_weight(samp, genPtH, mA, LHE_HT)


## Find significance given signal and background histograms
def find_significance(hSigX, hBkgX, wBkgX, mH, mA, cat):
    isData = ('JetHT' in hBkgX.GetName() or 'MET' in hBkgX.GetName() or \
              'SingleMu' in hBkgX.GetName() or 'EGam' in hBkgX.GetName())
    sumSSB  = 0  ## S^2/B
    sumSSBe = 0  ## S^2/(B+err)

    hSig = hSigX.Clone('hSig')
    hBkg = hBkgX.Clone('hBkg')
    wBkg = wBkgX.Clone('wBkg')

    ## Scale signal by 0.1 for 10% branching
    hSig.Scale(0.1)
    ## Scale background by ratio of WP60 / WP80
    iw = 4*(mH.startswith('mass1j'))  ## Use last 4 bins for mass1j (Haa34b WP)
    if isData:  ## Data excludes WP40 events
        WPscl = wBkg.Integral(iw+3,iw+4) / wBkg.Integral(iw+2,iw+3)
    else:
        WPscl = wBkg.Integral(iw+3,iw+4) / wBkg.Integral(iw+2,iw+4)
    if DEBUG:
        print('WP scaling bkg by %.3f' % WPscl)
    hBkg.Scale(WPscl)

    ## Rebin to 10 GeV for mAMed and mAHi, 20 GeV for mALo or Zvv
    if mH != 'logPt':
        rebin = 2
        if mA == 'mALo' or 'Wlv' in cat  or 'Zvv' in cat:
            rebin = 4
        hSig.Rebin(rebin)
        hBkg.Rebin(rebin)

    ## Scale background to S+B to represent floating background shape
    scl = (hSig.Integral() + hBkg.Integral()) / hBkg.Integral()
    if DEBUG: print('Sig scaling bkg by %.3f' % scl)

    maxSig = hSig.GetMaximum()
    pSig   = 0  ## Carry over previous signal from empty-background bins
    for ii in range(1, hSig.GetNbinsX()+1):
        if mH.startswith('mass')  and hSig.GetBinLowEdge(ii) < 49.9: continue
        if mH.startswith('mPNet') and hSig.GetBinLowEdge(ii) < 49.9: continue
        if mH.startswith('msoft') and hSig.GetBinLowEdge(ii) < 19.9: continue
        nSig = hSig.GetBinContent(ii)
        nBkg = hBkg.GetBinContent(ii)
        eBkg = hBkg.GetBinError(ii)
        if max(nBkg, eBkg) <= 0:
            pSig += nSig
            continue
        pSig = 0
        dSig    = nSig+nBkg - nBkg*scl
        dSigUp  = nSig+nBkg+eBkg - (nBkg+eBkg)*scl
        dSigDn  = nSig+nBkg-eBkg - (nBkg-eBkg)*scl
        dSigMin = min(abs(dSigUp), abs(dSigDn))
        sumSSB  += pow(dSig, 2) / (max(nBkg, eBkg)*scl)
        sumSSBe += pow(dSigMin, 2) / ((nBkg+eBkg)*scl)
        if DEBUG: print('Bin %d : %.1f sig, %.1f +/- %.1f bkg, dSig(Min) = %.1f (%.1f)' % \
                        (ii, nSig, nBkg, eBkg, dSig, dSigMin))
    ## End loop: for ii in range(1, hSig.GetNbinsX()+1)

    sigInt = hSig.Integral()
    bkgInt = hBkg.Integral()
    del hSig,hBkg,wBkg
    
    return [sigInt, bkgInt, np.sqrt(sumSSB), np.sqrt(sumSSBe)]
## End function: find_significance(hSigX, hBkgX, wBkgX, mH, mA, cat)

## Get selected muons
def get_muons(ch):
    ## Following https://indico.cern.ch/event/1430644/#2-update-on-higgs-aa-4b-booste
    muVecs = []
    muTrgs = []  ## List of vector indices for "triggerable" muons
    for iMu in range(ch.nMuon):
        if               ch.Muon_pt[iMu]  < 10. : continue
        if          abs(ch.Muon_eta[iMu]) > 2.4 : continue
        if ch.Muon_miniPFRelIso_all[iMu]  > 0.10: continue  ## Manual implementation of miniIsoId >= 3
        if   ch.Muon_mediumPromptId[iMu]  < 1   :
            if           ch.Muon_pt[iMu]  < 53. : continue
            if     ch.Muon_highPtId[iMu]  < 1   : continue
        if          abs(ch.Muon_dxy[iMu]) > 0.02: continue
        if           abs(ch.Muon_dz[iMu]) > 0.10: continue
        muVec = R.TLorentzVector()
        muVec.SetPtEtaPhiM(ch.Muon_pt[iMu], ch.Muon_eta[iMu], ch.Muon_phi[iMu], 0.106)
        muVecs.append(muVec)
        if ch.Muon_pt[iMu] < 26: continue
        muTrgs.append(len(muVecs)-1)
    return muVecs, muTrgs
## End function: get_muons(ch)

## Get selected electrons
def get_electrons(ch):
    ## Find electron(s) passing cuts
    ## Following https://indico.cern.ch/event/1430644/#2-update-on-higgs-aa-4b-booste
    eleVecs = []
    eleTrgs = []  ## List of vector indices for "triggerable" electrons
    for iEle in range(ch.nElectron):
        if                   ch.Electron_pt [iEle]  < 10. : continue
        if               abs(ch.Electron_eta[iEle]) > 2.5 : continue
        if abs(abs(ch.Electron_eta[iEle]) - 1.505) < 0.065: continue  ## Veto 1.44 - 1.57                                
        if    ch.Electron_mvaFall17V2Iso_WPL[iEle]  < 1   : continue
        if   ch.Electron_mvaFall17V2Iso_WP90[iEle]  < 1   :
            if                ch.Electron_pt[iEle]  < 35. : continue
            if     ch.Electron_cutBased_HEEP[iEle]  < 1   : continue
        ## Save the selected electron
        eleVec = R.TLorentzVector()
        eleVec.SetPtEtaPhiM(ch.Electron_pt[iEle], ch.Electron_eta[iEle], ch.Electron_phi[iEle], 0.0005)
        eleVecs.append(eleVec)
        if                  ch.Electron_pt[iEle] < 35. : continue
        if ch.Electron_mvaFall17V2Iso_WP90[iEle] < 1   : continue
        if ch.Electron_mvaFall17V2Iso_WP80[iEle] < 1   :
            if   ch.Electron_cutBased_HEEP[iEle] < 1   : continue
        if            abs(ch.Electron_dxy[iEle]) > 0.02: continue
        if            abs(ch.Electron_dz[iEle])  > 0.10: continue
        eleTrgs.append(len(eleVecs)-1)
    ## End loop: for iEle in range(ch.nElectron)
    return eleVecs, eleTrgs
## End function: get_electrons(ch)


def main():

    print('\nInside HtoAA_presel\n')

    chn = {}
    
    ## Dictionary for all histograms
    hst = {}

    ## Loop over categories and samples to create array of histograms
    for cat in CATS:
        hst[cat] = {}
        chn[cat] = {}
        for samp in SAMPS[cat]:
            isSig = ('HtoAA' in samp)
            hst[cat][samp] = {}
            if DO_EVT:
                ## Load input files and TChain 'Events' for this sample
                chn[cat][samp] = load_files(samp)
            ## Loop over mA, mH, pt, ms, bbt, bb to book histograms
            for mA in MASSA:
                hst[cat][samp][mA] = {}
                for mH in MASSH:
                    hst[cat][samp][mA][mH] = {}
                    for ptc in PTCUT:
                        pt = str(ptc)
                        hst[cat][samp][mA][mH][pt] = {}
                        for msc in MSCUT:
                            ## MRCUT on ratio of msoft / mass, so not needed if mH = msoft
                            ## MSCUT on msoft, so not needed if mH = msoft
                            if 'msoft' in mH and msc > MSCUT[0]: continue
                            # mr = str(mrc).replace('.','p')
                            ms = str(msc)
                            hst[cat][samp][mA][mH][pt][ms] = {}
                            for bbt in BBTAG:
                                hst[cat][samp][mA][mH][pt][ms][bbt] = {}
                                for bbc in BBCUT:
                                    bb = str(bbc).replace('.','p')
                                    hst[cat][samp][mA][mH][pt][ms][bbt][bb] = [0,0]
                                    ## If just reading histograms, array is now set; don't need to book TH1D
                                    if not DO_EVT: continue
                                    hname = 'h_%s_%s_%s_%s_pt%s_ms%s_%s_%s' % \
                                        (cat, samp, mA, mH, pt, ms, bbt, bb)
                                    htitle = '%s(AK8) : %s, %s, %s, pT > %d, ms > %d, %s > %.2f' % \
                                        (mH, cat, samp, mA, ptc, msc, bbt, bbc)
                                    hnameWP  = hname.replace(mH, 'WP')
                                    htitleWP = htitle.replace(mH+'(AK8)', 'Hto4b WPs')
                                    bins = ([30,7.4,10.4] if mH == 'logPt' else [60,0,300])
                                    ## Store mass distribution in [0], yields for different WPs in [1]
                                    hst[cat][samp][mA][mH][pt][ms][bbt][bb][0] = \
                                        R.TH1D(hname, htitle, bins[0], bins[1], bins[2])
                                    hst[cat][samp][mA][mH][pt][ms][bbt][bb][0].SetLineWidth(2)
                                    hst[cat][samp][mA][mH][pt][ms][bbt][bb][0].SetLineColor(R.kBlue if isSig else R.kBlack)
                                    hst[cat][samp][mA][mH][pt][ms][bbt][bb][0].SetDirectory(0)
                                    ## Store WP yield histogram just once
                                    if mH == MASSH[0]:
                                        hst[cat][samp][mA][mH][pt][ms][bbt][bb][1] = \
                                            R.TH1D(hnameWP, htitleWP, 8, -0.5, 7.5)
                                        hst[cat][samp][mA][mH][pt][ms][bbt][bb][1].SetLineWidth(2)
                                        hst[cat][samp][mA][mH][pt][ms][bbt][bb][1].SetLineColor(R.kBlue if isSig else R.kBlack)
                                        hst[cat][samp][mA][mH][pt][ms][bbt][bb][1].SetDirectory(0)
                                ## End loop: for bbc in BBCUT
                            ## End loop: for bbt in BBTAG
                        ## End loop: for msc in MSCUT
                    ## End loop: for ptc in PTCUT
                ## End loop: for mH in MASSH
            ## End loop: for mA in MASSA
        ## End loop: for samp in SAMPS[CAT]
    ## End loop: for cat in CATS

    ## --------------------------------------------------------------------
    ## Process events, filling histograms based on  mA, mH, pt, ms, bbt, bb

    ## 2018 DeepCSV and DeepFlavB cuts: https://btv-wiki.docs.cern.ch/ScaleFactors/UL2018/
    DeepB = {'L': 0.1208, 'M': 0.4168, 'T': 0.7665}
    FlavB = {'L': 0.0490, 'M': 0.2783, 'T': 0.7100}
    
    for cat in CATS:
        if not DO_EVT: break
        for samp in SAMPS[cat]:
            isSig = ('HtoAA' in samp)
            isData = ('JetHT' in samp or 'MET' in samp or \
                      'SingleMu' in samp or 'EGam' in samp)
            
            ## Don't need to separately process signal requiring 4 GEN b's in AK8
            if isSig and 'Gen' in samp: continue
            ## Loop through events, select, and plot
            ch = chn[cat][samp]
            nEntries = ch.GetEntries()
            print('\nFor %s %s, entering loop over %d events\n' % (cat, samp, nEntries))
            nPass = 0

            for iEvt in range(nEntries):
                if iEvt > MAX_EVT and MAX_EVT > 0: break
                if PSCALE > 1 and (iEvt % PSCALE) != 0: continue
                if (iEvt % PRT_EVT) == 0: print('Looking at event #%d / %d' % (iEvt, nEntries))
                ch.GetEntry(iEvt)

                ## If Zvv category, require MET
                if ('Zvv' in cat) and ch.MET_pt < 200: continue
                ## If Wlv category, require lepton
                if ('Wlv' in cat) and ch.nMuon + ch.nElectron < 1: continue

                ## Find AK8 jet(s) passing cuts
                xFat = -99
                xH34bVsQCD = -99
                xH4bVsQCD = -99
                fatVec = R.TLorentzVector()
                for iFat in range(ch.nFatJet):
                    ## Basic AK8 jet cuts
                    if PTMIN > 0 and ch.FatJet_pt[iFat] < PTMIN: continue
                    if PTMAX > 0 and ch.FatJet_pt[iFat] > PTMAX: continue
                    if ch.FatJet_msoftdrop[iFat] <  20: continue
                    ## if ch.FatJet_mass[iFat]      <  50: continue
                    if  abs(ch.FatJet_eta[iFat]) > 2.4: continue
                    if ch.FatJet_jetId[iFat]     <   6: continue
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
                    ## Cut at PNet (v1) WPL initially to store WP yields
                    if H34bVsQCD < 0.8 and H4bVsQCD < 0.80: continue
                    ## Pick the highest-score AK8 jet
                    if H34bVsQCD > xH34bVsQCD or (H34bVsQCD < 0.80 and H4bVsQCD > xH4bVsQCD):
                        xFat = iFat
                        xH34bVsQCD = H34bVsQCD
                        xH4bVsQCD  = H4bVsQCD
                        # break  ## Temporarily suspend sorting for speed - AWB 2024.08.08
                ## End loop: for iFat in range(ch.nFatJet)
                if xFat < 0: continue
                ## Save 4-vector
                fatVec.SetPtEtaPhiM(ch.FatJet_pt[xFat], ch.FatJet_eta[xFat], \
                                    ch.FatJet_phi[xFat], ch.FatJet_mass[xFat] )

                ## Save MET 4-vector
                metVec = R.TLorentzVector()
                metVec.SetPtEtaPhiM(ch.MET_pt, 0, ch.MET_phi, 0)
                
                ## For Zvv category, require dPhi(MET, AK8) > pi/2
                if ('Zvv' in cat) and abs(metVec.DeltaPhi(fatVec)) < (np.pi/2.0): continue

                ## Get selected leptons to "clean" jets and veto events
                muVecs,  muTrgs  = get_muons(ch)
                eleVecs, eleTrgs = get_electrons(ch)
                lepVecs = muVecs + eleVecs
                ## No trigger leptons allowed in non-lepton categories
                if (('gg0l' in cat) or ('Zvv' in cat)) and len(muTrgs)+len(eleTrgs) > 0: continue
                ## Require trigger lepton for lepton categories
                if ('Wlv' in cat) and len(muTrgs)+len(eleTrgs) < 1: continue
                
                ## Require exactly one lepton for Wlv category, with dPhi(AK8,lv) > 3*pi/4
                if ('Wlv' in cat):
                    if len(lepVecs) > 1: continue
                    if abs(fatVec.DeltaPhi(metVec+lepVecs[0])) < (3*np.pi/4.0): continue
                    if 'SingleMu' in samp and len(muVecs) < 1: continue
                    if ('EGam' in samp) and len(eleVecs) < 1: continue

                ## "Clean" AK8 jet w.r.t. all selected leptons
                lepVeto = False
                for lepVec in muVecs+eleVecs:
                    if lepVec.DeltaR(fatVec) < 0.8:
                        lepVeto = True
                        if not 'ZHtoAA' in samp and not 'ttbar_' in samp:
                            print('*** dR(lep,AK8) = %.3f, lep pT = %.1f. Vetoing! ***' % \
                                  (lepVec.DeltaR(fatVec), lepVec.Pt()))
                        break
                if lepVeto: continue

                ## Find AK4 b-tagged jet with dR(AK8, AK4) closest to 0.8
                xJet = -99  ## b-tagged AK4 with lowest dR > 0.8
                xdR84 = 99.
                xJetVec = R.TLorentzVector()
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
                    ## Clean w.r.t. AK8 jet and leptons
                    if iJetVec.DeltaR(fatVec) < 0.8: continue
                    for lepVec in lepVecs:
                        if iJetVec.DeltaR(lepVec) < 0.4: continue
                    if iJetVec.Pt() > 30:
                        nBTags += 1
                    if iJetVec.DeltaR(fatVec) > AK4DR: continue
                    if iJetVec.DeltaR(fatVec) < xdR84:
                        xJet    = iJet
                        xJetVec = iJetVec
                        xdR84   = iJetVec.DeltaR(fatVec)
                        ## Count low-pT AK4 b-tags only if they can form AK8+AK4
                        if iJetVec.Pt() < 30:
                            nBTags += 1
                ## End loop: for iJet in range(ch.nJet)

                ## Distinguish between "1j" and "0j" events
                if xJet >= 0:
                    if xH34bVsQCD < 0.80: continue
                else:
                    if xH4bVsQCD < 0.80: continue

                ## Require no extra b-tagged AK4 jets for all non-ttbar categories
                if (not ('tt' in cat)) and nBTags > 1*(xJet >= 0): continue

                
                ## Remove events with no good collision vertices
                if ch.PV_npvsGood < 1:
                    continue
                ## Remove events failing MET noise filters
                if not (ch.Flag_goodVertices and ch.Flag_globalSuperTightHalo2016Filter and \
                        ch.Flag_HBHENoiseFilter and ch.Flag_HBHENoiseIsoFilter and \
                        ch.Flag_EcalDeadCellTriggerPrimitiveFilter and ch.Flag_BadPFMuonFilter and \
                        ch.Flag_BadPFMuonDzFilter and ch.Flag_eeBadScFilter and ch.Flag_ecalBadCalibFilter):
                    continue
                ## Require relevant triggers
                if TRIG and ('gg0l' in cat) and not \
                   ( ( (ch.HLT_PFJet500 or ch.HLT_AK8PFJet500 or ch.HLT_AK8PFJet400_TrimMass30 or \
                        ch.HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4) and ch.L1_SingleJet180) or \
                     ( (ch.HLT_PFHT1050 or ch.HLT_AK8PFHT800_TrimMass50) and \
                       (ch.L1_SingleJet180 or ch.L1_HTT360er) ) or \
                     ( (ch.HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71) and \
                       (ch.L1_DoubleJet112er2p3_dEta_Max1p6 or ch.L1_DoubleJet150er2p5) ) ): continue
                if TRIG and ('Zvv' in cat) and not \
                   ( ((ch.HLT_PFMET120_PFMHT120_IDTight_PFHT60 or ch.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60) and \
                      (ch.L1_ETMHF90_HTT60er or ch.L1_ETMHF100_HTT60er or ch.L1_ETMHF110_HTT60er)) or \
                     ((ch.HLT_PFMETTypeOne140_PFMHT140_IDTight or ch.HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned or \
                       ch.HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1) and \
                      (ch.L1_ETMHF100 or ch.L1_ETMHF110 or ch.L1_ETMHF120 or ch.L1_ETMHF130)) ): continue

                ## -----------------------------------------------------------------
                ## Done with basic selection cuts, fill WP yield and mass histograms
                nPass += 1

                ## Get GEN pT(H) for ggH reweighting
                genPtH = 0
                if samp.startswith('ggHtoAA'):
                    for iGen in range(ch.nGenPart):
                        if ch.GenPart_pdgId[iGen] == 25 and ch.GenPart_status[iGen] == 62:
                            genPtH = ch.GenPart_pt[iGen]
                            break
                # ## Temporary hack to reduce processing time - AWB 2024.08.08
                # genPtH = ch.FatJet_pt[xFat]
                
                ## Get mass(a)
                massA = (ch.FatJet_particleNet_massA_Hto4b_v0[xFat] + \
                         ch.FatJet_particleNet_massA_Hto4b_v1[xFat] + \
                         ch.FatJet_particleNet_massA_Hto4b_v3[xFat]) / 3.0
                mA = ('mALo' if massA < 22. else ('mAHi' if massA > 47 else 'mAMed'))
                if not mA in MASSA: continue

                ## Get pT(H) and mass(H)
                ptH = ch.FatJet_pt[xFat]
                MH = {}
                MH['logPt'] = min(max(np.log2(ptH), 7.41), 10.39)
                MH['mass']  = ch.FatJet_mass[xFat] if xJet < 0 else 1
                for xMH in MASSH:
                    if xMH.startswith('mass1j'):
                        MH[xMH] = (fatVec+xJetVec).M() if xJet >= 0 else 1
                MH['msoft'] = ch.FatJet_msoftdrop[xFat]
                MH['mPNet'] = (ch.FatJet_particleNet_massH_Hto4b_v0[xFat]*1.01 + \
                               ch.FatJet_particleNet_massH_Hto4b_v1[xFat] + \
                               ch.FatJet_particleNet_massH_Hto4b_v2[xFat] + \
                               ch.FatJet_particleNet_massH_Hto4b_v3[xFat]) / 4.0
                ## See last slide of https://indico.cern.ch/event/1430644/contributions/6018705
                scalePt = 1 + (8.644 / 15) * (1 - pow(2, -1*pow(8.644 - np.log2(ptH), 2)))
                MH['massPt']  = MH['mass'] * (1 if ptH >= 400 else scalePt)
                MH['msoftPt'] = MH['msoft']* (1 if ptH >= 400 else scalePt)

                ## Get double b-tag scores
                bbX = {}
                bbPN_num = ch.FatJet_particleNetMD_Xbb[xFat]
                bbPN_den = bbPN_num + ch.FatJet_particleNetMD_QCD[xFat]
                bbX['PNetXbb'] = (0 if bbPN_den <= 0 else (bbPN_num / bbPN_den))
                bbX['deepTagXbb'] = max(ch.FatJet_deepTagMD_ZHbbvsQCD[xFat], 0.0)
                bbX['avgXbb'] = (bbX['PNetXbb'] + bbX['deepTagXbb']) / 2.0

                ## Compute weights
                WGT = 1.0 if isData else get_weight(samp, genPtH, mA, ch.LHE_HT)
                if PSCALE > 1: WGT *= PSCALE

                ## Loop to fill histograms requiring particular cuts
                for mH in MASSH:
                    if mH == 'mass' and xJet >= 0: continue
                    if mH.startswith('mass1j'):
                        if xJet < 0: continue
                        if '1j1p' in mH:
                            if xdR84 > float((mH.replace('mass1j','')).replace('p','.')):
                                continue
                        if mH.replace('mass1j','') in ['20','25','30']:
                            if xJetVec.Pt() < float(mH.replace('mass1j','')):
                                continue
                    for ptc in PTCUT:
                        pt = str(ptc)
                        if ptH < ptc: continue
                        for msc in MSCUT:
                            # mr = str(mrc).replace('.','p')
                            ms = str(msc)
                            ## MRCUT on ratio of msoft / mass, so not needed if mH = msoft
                            ## MSCUT on msoft, so not needed if mH = msoft
                            if 'msoft' in mH and msc > MSCUT[0]: continue
                            if (MH[mH] if mH != 'logPt' else MH['mass']) <= 0:
                                print('\n\n*** Super weird error in LS = %d, event = %d!!!\n\n' % (ch.luminosityBlock, ch.event))
                                print('MH[msoft] = %.6f, MH[%s] = %.6f, MH[mass] = %.6f' % (MH['msoft'], MH[mH], MH['mass']))
                                continue
                            # if (MH['msoft'] / (MH[mH] if mH != 'logPt' else MH['mass'])) < mrc: continue
                            if MH['msoft'] < msc: continue
                            for bbt in BBTAG:
                                for bbc in BBCUT:
                                    bb = str(bbc).replace('.','p')
                                    if bbX[bbt] < bbc: continue

                                    ## Fill PNet (v1) WP80/60/40
                                    if mH == MASSH[0]:
                                        WP4 = (xH4bVsQCD > 0.92) + (xH4bVsQCD > 0.975) + (xH4bVsQCD > 0.992)
                                        WP34 = (xH34bVsQCD > 0.92) + (xH34bVsQCD > 0.96) + (xH34bVsQCD > 0.98)
                                        WPx = (xJet < 0)*WP4 + (xJet >= 0)*(4 + WP34)
                                        hst[cat][samp][mA][mH][pt][ms][bbt][bb][1].Fill(WPx, WGT)
                                        #if ((not isSig) and mH == MASSH[0] and pt == '170' and ms == '20' and bbt == BBTAG[0] and bb == '0p0'):
                                            #print('Filled %s with %d x %.6f' % (hst[cat][samp][mA][mH][pt][ms][bbt][bb][1].GetName(), WPx, WGT))
                                        ## Plot signal with 4 GEN b's in AK8 separately
                                        # if isSig and ( (xJet <  0 and ch.FatJet_nBHadrons[xFat] >= 4) or \
                                        #                (xJet >= 0 and ch.FatJet_nBHadrons[xFat] == 3) ):
                                        if isSig and ch.FatJet_nBHadrons[xFat] >= 4:
                                            hst[cat][samp+'Gen4'][mA][mH][pt][ms][bbt][bb][1].Fill(WPx, WGT)
                                        if isSig and ch.FatJet_nBHadrons[xFat] == 3:
                                            hst[cat][samp+'Gen3'][mA][mH][pt][ms][bbt][bb][1].Fill(WPx, WGT)

                                    ## For mass plots, cut at PNet (v1) WP60 for signal, WP80 for background
                                    ## Veto WP40 for data
                                    if xJet <  0 and xH4bVsQCD < (0.975 if isSig else 0.92): continue
                                    if xJet >= 0 and xH34bVsQCD < (0.96 if isSig else 0.92): continue
                                    if isData and xJet <  0 and xH4bVsQCD > 0.992: continue
                                    if isData and xJet >= 0 and xH34bVsQCD > 0.98: continue
                                    hst[cat][samp][mA][mH][pt][ms][bbt][bb][0].Fill(MH[mH], WGT)
                                    #if ((not isSig) and mH == 'mass' and pt == '170' and ms == '20' and bbt == BBTAG[0] and bb == '0p0'):
                                        #print('Filled %s with %d x %.6f' % (hst[cat][samp][mA][mH][pt][ms][bbt][bb][0].GetName(), MH[mH], WGT))
                                    ## Plot signal with 4 GEN b's in AK8(+AK4) separately
                                    # if isSig and ( (xJet <  0 and ch.FatJet_nBHadrons[xFat] >= 4) or \
                                    #                (xJet >= 0 and ch.FatJet_nBHadrons[xFat] == 3) ):
                                    if isSig and ch.FatJet_nBHadrons[xFat] >= 4:
                                        hst[cat][samp+'Gen4'][mA][mH][pt][ms][bbt][bb][0].Fill(MH[mH], WGT)
                                    if isSig and ch.FatJet_nBHadrons[xFat] == 3:
                                        hst[cat][samp+'Gen3'][mA][mH][pt][ms][bbt][bb][0].Fill(MH[mH], WGT)
                                ## End loop: for bbc in BBCUT
                            ## End loop: for bbt in BBTAG
                        ## End loop: for msc in MSCUT
                    ## End loop: for ptc in PTCUT
                ## End loop: for mH in MASSH
                
            ## End loop: for iEvt in range(nEntries)
            print('Finished with event loop in %s %s: %d passed' % (cat, samp, nPass))
        ## End loop: for samp in SAMPS[CAT]
    ## End loop: for cat in CATS
                
    ## End event processing and histogram filling
    ## --------------------------------------------------------------------


    ## ------------------------------------------------------------------------
    ## Loop over cat, samp, mA, mH, pt, ms, bbt, bb to write or load histograms
    for cat in CATS:
        out_file = None
        for samp in SAMPS[cat]:
            isSig = ('HtoAA' in samp)
            ## Create or load output file containing histograms for this sample
            out_dir = 'plots/'
            out_file_str = out_dir+'HtoAA_presel_%s_%s' % (cat, samp)
            if TRIG: out_file_str += '_trig'
            if 'mass1j' in MASSH or not DO_EVT:
                out_file_str += '_01j'
                out_file_str += '_ptj%d' % AK4PT
                out_file_str += '_dRj%s' % str(AK4DR).replace('.','p')
            if 'mass1j1p4' in MASSH: out_file_str += '_dRs'
            if 'mass1j20' in MASSH: out_file_str += '_pTs'
            if PTMIN > 0: out_file_str += '_ptMin%d' % PTMIN
            if PTMAX > 0: out_file_str += '_ptMax%d' % PTMAX
            if MAX_EVT > 0: out_file_str += '_%dk' % (MAX_EVT / 1000)
            if PSCALE  > 1: out_file_str += '_ps%d' % PSCALE
            if PSCALE < 1 and not DO_EVT and samp.startswith('ggHtoAA'):
                out_file_str += '_ps10'
                print('\n\n*** Using PS = 10 for ggHtoAA signal ONLY! ***')
            out_file_str += '.root'
            nHist = 0
            if DO_EVT:
                out_file = R.TFile(out_file_str,'recreate')
                out_file.cd()
            else:
                out_file = R.TFile(out_file_str,'read')
                print('Opened '+out_file_str)
                
            for mA in MASSA:
                for mH in MASSH:
                    for ptc in PTCUT:
                        pt = str(ptc)
                        for msc in MSCUT:
                            ## MRCUT on ratio of msoft / mass, so not needed if mH = msoft
                            ## MSCUT on msoft, so not needed if mH = msoft
                            if 'msoft' in mH and msc > MSCUT[0]: continue
                            # mr = str(mrc).replace('.','p')
                            ms = str(msc)
                            for bbt in BBTAG:
                                for bbc in BBCUT:
                                    bb = str(bbc).replace('.','p')
                                    if DO_EVT:
                                        hst[cat][samp][mA][mH][pt][ms][bbt][bb][0].Write()
                                        if mH == MASSH[0]:
                                            hst[cat][samp][mA][mH][pt][ms][bbt][bb][1].Write()
                                    else:
                                        mAx = ('mAMed' if mA == 'mAMedHi' else mA)
                                        hname = 'h_%s_%s_%s_%s_pt%s_ms%s_%s_%s' % \
                                            (cat, samp, mAx, mH, pt, ms, bbt, bb)
                                        hst[cat][samp][mA][mH][pt][ms][bbt][bb][0] = out_file.Get(hname)
                                        hst[cat][samp][mA][mH][pt][ms][bbt][bb][0].SetDirectory(0)
                                        if mA == 'mAMedHi':
                                            hnameHi = hname.replace('mAMed','mAHi')
                                            htitleX = ((out_file.Get(hname)).GetTitle()).replace('mAMed','mAMedHi')
                                            hst[cat][samp][mA][mH][pt][ms][bbt][bb][0].Add(out_file.Get(hnameHi))
                                            hst[cat][samp][mA][mH][pt][ms][bbt][bb][0].SetTitle(htitleX)
                                        if mH == MASSH[0]:
                                            hnameWP  = hname.replace(mH, 'WP')
                                            hst[cat][samp][mA][mH][pt][ms][bbt][bb][1] = out_file.Get(hnameWP)
                                            hst[cat][samp][mA][mH][pt][ms][bbt][bb][1].SetDirectory(0)
                                            if mA == 'mAMedHi':
                                                hnameWPHi = hnameWP.replace('mAMed','mAHi')
                                                htitleWPX = ((out_file.Get(hnameWP)).GetTitle()).replace('mAMed','mAMedHi')
                                                hst[cat][samp][mA][mH][pt][ms][bbt][bb][1].Add(out_file.Get(hnameWPHi))
                                                hst[cat][samp][mA][mH][pt][ms][bbt][bb][1].SetTitle(htitleWPX)
                                    ## End conditional: if DO_EVT / else
                                    nHist += (1 + (mH == MASSH[0]))
                                ## End loop: for bbc in BBCUT
                            ## End loop: for bbt in BBTAG
                        ## End loop: for msc in MSCUT
                    ## End loop: for ptc in PTCUT
                ## End loop: for mH in MASSH
            ## End loop: for mA in MASSA
            print('\n%s %d histograms for %s %s\n' % (('Wrote' if DO_EVT else 'Read'), nHist, cat, samp))
            out_file.Close()
        ## End loop: for samp in SAMPS[CAT]
    ## End loop: for cat in CATS

    ## End event processing and histogram writing
    ## --------------------------------------------------------------------

    
    ## Loop over cat, mA, mH, pt, ms, bbt, bb to estimate sensitivity
    print('Mass, Category, Signal, Mass(a), pT cut, msoft cut, Xbb cut, nSig, nBkg, SBB, SBBe')
    for cat in CATS:
        if DO_EVT: break
        for mA in MASSA:
            for mH in MASSH:
                for ptc in PTCUT:
                    pt = str(ptc)
                    for msc in MSCUT:
                        ## MRCUT on ratio of msoft / mass, so not needed if mH = msoft
                        ## MSCUT on msoft, so not needed if mH = msoft
                        if 'msoft' in mH and msc > MSCUT[0]: continue
                        # mr = str(mrc).replace('.','p')
                        ms = str(msc)
                        for bbt in BBTAG:
                            for bbc in BBCUT:
                                bb = str(bbc).replace('.','p')
                                hSig = None
                                for sampS in SAMPS[cat]:
                                    if not 'HtoAA' in sampS: continue
                                    hSig = hst[cat][sampS][mA][mH][pt][ms][bbt][bb][0]
                                    if DEBUG: print('Initialize sig %s' % hSig.GetName())
                                    hBkg = None
                                    wBkg = None
                                    for sampB in SAMPS[cat]:
                                        if 'HtoAA' in sampB: continue
                                        if not hBkg:
                                            hBkg = hst[cat][sampB][mA][mH][pt][ms][bbt][bb][0]
                                            wBkg = hst[cat][sampB][mA][MASSH[0]][pt][ms][bbt][bb][1]
                                            if DEBUG: print('Initialize bkg %s' % hBkg.GetName())
                                        else:
                                            hBkg.Add(hst[cat][sampB][mA][mH][pt][ms][bbt][bb][0])
                                            wBkg.Add(hst[cat][sampB][mA][MASSH[0]][pt][ms][bbt][bb][1])
                                            if DEBUG: print('Add bkg %s' % hst[cat][sampB][mA][mH][pt][ms][bbt][bb][0].GetName())
                                    ## End loop: for sampB in SAMPS
                                    ssb = find_significance(hSig, hBkg, wBkg, mH, mA, cat)
                                    # print('%s, S = %.2f, B = %.2f, SSB = %.4f, SSBe = %.4f' % \
                                        #       (hSig.GetTitle(), hSig.Integral(), hBkg.Integral(), ssb[0], ssb[1]))
                                    print('%s, %.2f, %.2f, %.4f, %.4f' % \
                                          (hSig.GetTitle().replace('(AK8) :',',').replace('pT > ','').replace('ms > ',''), \
                                           ssb[0], ssb[1], ssb[2], ssb[3]))
                                ## End loop: for sampS in SAMPS
                            ## End loop: for bbc in BBCUT
                        ## End loop: for bbt in BBTAG
                    ## End loop: for msc in MSCUT
                ## End loop: for ptc in PTCUT
            ## End loop: for mH in MASSH
        ## End loop: for mA in MASSA
    ## End loop: for cat in CATS
                                
    ## End sensitivity estimation
    ## --------------------------------------------------------------------
    

if __name__ == '__main__':
    main()
