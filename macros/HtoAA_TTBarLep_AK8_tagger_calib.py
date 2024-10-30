#! /usr/bin/env python
## Script to generate input ROOT files for Combine fit of AK8 tagger efficiencies
## Before running, set up same CMSSW environment as Combine release:
## cd /afs/cern.ch/work/a/abrinke1/public/HiggsToAA/2DAlphabet/CMSSW_11_3_4/src/; cmssw-cc7
## cmsenv; cd /afs/cern.ch/work/a/abrinke1/public/HiggsToAA/coffea/eventloop/
## python3 macros/HtoAA_TTBarLep_AK8_tagger_calib.py

import os
import sys
import math
import ROOT as R
import ctypes
from array import array

R.gROOT.SetBatch(True)  ## Don't display histograms or canvases when drawn
R.gStyle.SetOptStat(0)  ## Don't display stat boxes

## User configuration
VERBOSE  = True
IN_DIR   = '/afs/cern.ch/work/c/csutanta/public/rootfiles/unskimmed_singlemuon_'
CATS     = ['bdtVeto',
            'bdtLo',
            'bdtMed',
            'bdtHi']
SELS     = ['0b_BBQ']
TAGGERS  = {'particleNetMD_XbbOverQCD':[0.1,0.5,0.75]}
#TAGGERS  = {'FatJet_PNetMD_Hto4b_Htoaa34bOverQCD':[0.8]}
TAGNM    = {'particleNetMD_XbbOverQCD':'Xbb',
            'FatJet_PNetMD_Hto4b_Htoaa34bOverQCD':'Htoaa34b'}
SVAR  = 2.0  ## Systematic factor of variation in tagging efficiency

YEAR     = '2018'
DATA     = 'SingleMuon'
ERAS     = ['Run'+YEAR+er for er in ['A','B','C','D']]

MC_2bq   = ['TTToSemiLeptonic_powheg_bbqq',
            'TTToSemiLeptonic_powheg_bbq']
MC_2b    = ['TTToSemiLeptonic_powheg_bb',
            'TTTo2L2Nu_powheg_bb']
MC_bqq   = ['TTToSemiLeptonic_powheg_bqq']
MC_1b    = ['TTToSemiLeptonic_powheg_1b',
            'TTTo2L2Nu_powheg_1b',
            'SingleTop']
MC_0b    = ['TTToSemiLeptonic_powheg_0b',
            'TTTo2L2Nu_powheg_0b',
            'WJetsToLNu_HT_LO',
            'DYJets_M-50_Incl_NLO',
            'DYJets_M-10to50_Incl_NLO']
MCNM = {'TTToSemiLeptonic_powheg_bbqq':'TT1L_bbqq',
        'TTToSemiLeptonic_powheg_bbq' :'TT1L_bbq',
        'TTToSemiLeptonic_powheg_bb'  :'TT1L_bb',
        'TTTo2L2Nu_powheg_bb'         :'TT2L_bb',
        'TTToSemiLeptonic_powheg_bqq' :'TT1L_bqq',
        'TTToSemiLeptonic_powheg_1b'  :'TT1L_1b',
        'TTTo2L2Nu_powheg_1b'         :'TT2L_2b',
        'SingleTop'                   :'SingleT',
        'TTToSemiLeptonic_powheg_0b'  :'TT1L_0b',
        'TTTo2L2Nu_powheg_0b'         :'TT2L_0b',
        'WJetsToLNu_HT_LO'            :'WToLNu',
        'DYJets_M-50_Incl_NLO'        :'DY_M50',
        'DYJets_M-10to50_Incl_NLO'    :'DY_M10'}

OUT_DIR  = 'plots/HtoAA_TTBarLep_AK8_tagger_calib/'


## Loop over bins in original histogram to find boundary bin
def find_boundary(h_in, cut):
    nBin = h_in.GetNbinsX()
    iBin = -99
    for iX in range(1, nBin+1):
        if h_in.GetBinCenter(iX) < cut:
            iBin = iX
        else:
            break
    return iBin
## End function: find_boundary(h_in, cut)


## Fill rebinned pass-fail histogram with integrals from original histogram
def fill_pass_fail(h_in, h_out, cuts):
    ## Add low- and high-score boundaries
    cutsX = [-999.]+cuts+[999.]
    nBin = h_in.GetNbinsX()
    for iCut in range(1, len(cutsX)):
        cutLo = cutsX[iCut-1]
        cutHi = cutsX[iCut]
        iBinLo = max(find_boundary(h_in, cutLo), 0) + 1
        iBinHi = find_boundary(h_in, cutHi)
        inA  = h_out.GetBinContent(iCut)
        inAe = h_out.GetBinError(iCut)
        inBe = ctypes.c_double()
        inB  = h_in.IntegralAndError(iBinLo, iBinHi, inBe)
        inBe = float(str(inBe).replace('c_double(','').replace(')',''))
        h_out.SetBinContent(iCut, inA+inB)
        h_out.SetBinError(iCut, math.sqrt(pow(inAe,2) + pow(float(inBe),2)))
        if VERBOSE: print('  * %s(%.1f, %.1f)(%d, %d) = %.2f +/- %.2f' % (h_in.GetName(), cutLo, cutHi, iBinLo, iBinHi, inA+inB, h_out.GetBinError(iCut)))
    ## End loop: for iCut in range(1, len(cutsX))
    return
## End function: fill_pass_fail(h_in, h_out, cut)


## Make systematic variations
def make_syst_hists(h_out, syst):
    if VERBOSE: print('  * Integral %s = %.1f' % (h_out.GetName(), h_out.Integral()))
    h_outs_up_dn = []
    for iBin in range(1, h_out.GetNbinsX()+1):
        bStr = '_%dbin' % iBin
        h_out_up = h_out.Clone(h_out.GetName()+'_s'+syst+bStr+'Up')
        h_out_dn = h_out.Clone(h_out.GetName()+'_s'+syst+bStr+'Down')
        iN   = h_out.GetBinContent(iBin)
        iNe  = h_out.GetBinError(iBin)
        tot  = h_out.Integral()
        oN   = tot - iN
        ## No shape variation if bin contains all or no events
        if iN == 0 or oN == 0:
            h_outs_up_dn.append(h_out_up)
            h_outs_up_dn.append(h_out_dn)
            continue
        iNup = min(iN*(SVAR-1.0), oN*(SVAR-1.0)/SVAR)  ## How much to add for up variation
        oFup = (oN - iNup) / oN                        ## Factor of reduction for other bins
        iNdn = min(iN*(SVAR-1.0)/SVAR, oN*(SVAR-1.0))  ## How much to subtract for down variation
        oFdn = (oN + iNdn) / oN                        ## Factor of increase for other bins
        if VERBOSE: print('    - Bin %d %.2f --> +%.2f/-%.2f' % (iBin, iN, iNup, iNdn))
        h_out_up.SetBinContent(iBin, iN+iNup)
        h_out_dn.SetBinContent(iBin, iN-iNdn)
        h_out_up.SetBinError(iBin, iNe*(iN+iNup)/iN)
        h_out_dn.SetBinError(iBin, iNe*(iN-iNdn)/iN)
        for xBin in range(1, h_out.GetNbinsX()+1):
            if xBin == iBin: continue
            xN  = h_out.GetBinContent(xBin)
            xNe = h_out.GetBinError(xBin)
            h_out_up.SetBinContent(xBin, xN*oFup)
            h_out_dn.SetBinContent(xBin, xN*oFdn)
            h_out_up.SetBinError(xBin, xNe*oFup)
            h_out_dn.SetBinError(xBin, xNe*oFdn)
        ## End loop: for xBin in range(1, h_out.GetNbinsX()+1)
        if abs((h_out_up.Integral() / h_out.Integral()) - 1.0) > 0.001 or \
           abs((h_out_dn.Integral() / h_out.Integral()) - 1.0) > 0.001:
            print('\n\nBAD ERROR!!! Integral = %.3f, up = %.3f, down = %.3f.\n\n' % (h_out.Integral(), h_out_up.Integral(), h_out_dn.Integral()))
            sys.exit()
        h_outs_up_dn.append(h_out_up)
        h_outs_up_dn.append(h_out_dn)
    ## End loop: for iBin in range(1, h_out.GetNbinsX()+1)
    return h_outs_up_dn
## End function: make_syst_hists(h_out, syst)

    
def main():

    print('\nInside HtoAA_TTBarLep_AK8_tagger_calib\n')

    if not os.path.exists(OUT_DIR):
        os.mkdir(OUT_DIR)

    ## Initialize dictionary of output histograms
    h_outs = {}
    for sel in SELS:
        h_outs[sel] = {}
        
    for cat in CATS:
        if VERBOSE: print('\nNow looking at category %s' % cat)
        in_file_str = IN_DIR+cat+'.root'
        in_file = R.TFile(in_file_str, 'open')
        print('\n*******\nReading from %s' % in_file_str)

        for sel in SELS:
            if VERBOSE: print('\nNow looking at selection %s' % sel)
            for tag in TAGGERS.keys():
                nCuts = len(TAGGERS[tag])
                if VERBOSE: print('\nNow looking at tagger %s, will rebin to %d bins' % (tag, nCuts+1))

                ####################################
                ## Get data histograms from each era
                for era in ERAS:
                    if VERBOSE: print('\nNow looking at era %s' % era)
                    h_in_name = 'evt/'+DATA+'_'+era+'/'+tag+'_sel_'+sel+'_noweight'
                    h_in = in_file.Get(h_in_name)
                    if VERBOSE: print('\nGot histogram %s' % h_in_name)
                    if VERBOSE: print('  * Integral = %.1f' % h_in.Integral())

                    ## Fill new rebinned histogram with events in tagger ranges
                    h_out_name = DATA+'_'+YEAR+'_'+cat+'_'+TAGNM[tag]
                    if not h_out_name in h_outs[sel].keys():
                        h_outs[sel][h_out_name] = R.TH1D(h_out_name, h_out_name, nCuts+1, 0, nCuts+1)
                        h_outs[sel][h_out_name].SetDirectory(0) ## Save locally
                    fill_pass_fail(h_in, h_outs[sel][h_out_name], TAGGERS[tag])
                ## End loop: for era in ERAS

                ####################
                ## Get MC histograms
                for mc in MC_2bq+MC_2b+MC_bqq+MC_1b+MC_0b:
                    if VERBOSE: print('\nNow looking at MC sample %s' % mc)
                    h_in_name = 'evt/'+mc+'/'+tag+'_sel_'+sel+'_central'
                    h_in = in_file.Get(h_in_name)
                    if VERBOSE: print('\nGot histogram %s' % h_in_name)
                    if VERBOSE: print('  * Integral = %.1f' % h_in.Integral())

                    ## Fill new rebinned histogram with events failing and passing cuts
                    mcm = MCNM[mc]
                    h_out_name = mcm+'_'+cat+'_'+TAGNM[tag]
                    if not h_out_name in h_outs[sel].keys():
                        h_outs[sel][h_out_name] = R.TH1D(h_out_name, h_out_name, nCuts+1, 0, nCuts+1)
                        h_outs[sel][h_out_name].SetDirectory(0) ## Save locally
                    else:
                        print('\n\nHow is %s already in h_outs[%s]??? Quitting.' % (h_out_name, sel))
                        sys.exit()
                    fill_pass_fail(h_in, h_outs[sel][h_out_name], TAGGERS[tag])
                    
                    if mc in MC_2bq: systs = ['2bq','2B']
                    if mc in MC_2b:  systs = ['2b', '2B']
                    if mc in MC_bqq: systs = ['bqq','1B']
                    if mc in MC_1b:  systs = ['1b', '1B']
                    if mc in MC_0b:  systs = ['0b']
                    for syst in systs:
                        ## Generate additional histograms with sum of MC
                        h_MC_name = 'Sum'+syst+'_'+cat+'_'+TAGNM[tag]
                        if not h_MC_name in h_outs[sel].keys():
                            h_outs[sel][h_MC_name] = R.TH1D(h_MC_name, h_MC_name, nCuts+1, 0, nCuts+1)
                            h_outs[sel][h_MC_name].SetDirectory(0) ## Save locally
                        h_outs[sel][h_MC_name].Add(h_outs[sel][h_out_name])

                        ## Perform systematic variations
                        for h_syst in make_syst_hists(h_outs[sel][h_out_name], syst):
                            h_outs[sel][h_syst.GetName()] = h_syst
                            h_outs[sel][h_syst.GetName()].SetDirectory(0) ## Save locally
                            ## Generate additional systematic histograms with sum of MC
                            h_MC_name_syst = h_MC_name+(h_syst.GetName().replace(h_out_name,''))
                            if not h_MC_name_syst in h_outs[sel].keys():
                                h_outs[sel][h_MC_name_syst] = R.TH1D(h_MC_name_syst, h_MC_name_syst, nCuts+1, 0, nCuts+1)
                                h_outs[sel][h_MC_name_syst].SetDirectory(0) ## Save locally
                            h_outs[sel][h_MC_name_syst].Add(h_syst)
                        ## End loop: for h_syst in make_syst_hists(h_outs[sel][h_out_name], syst)
                    ## End loop: for syst in systs
                ## End loop: for mc in MC_2bq+MC_2b+MC_bqq+MC_1b+MC_0b
                
            ## End loop: for tag in TAGGERS.keys()
        ## End loop: for sel in SELS
        in_file.Close()
    ## End loop: for cat in CATS

    
    ## Create a separate output ROOT file for each selection option
    for sel in SELS:
        out_file_str = OUT_DIR+'AK8_tagger_calib_%s_%s_slc7.root' % (sel, str(SVAR).replace('.','p'))
        out_file = R.TFile(out_file_str, 'recreate')
        print('\n*******\nWriting to %s' % out_file_str)
        for h_out_name in h_outs[sel].keys():
            if VERBOSE: print('Writing out %s' % h_out_name)
            if VERBOSE: print('  * Integral = %.3f' % h_outs[sel][h_out_name].Integral())
            if VERBOSE: print('  * Bins = '+', '.join('%.3f' % h_outs[sel][h_out_name].GetBinContent(iX)
                                                      for iX in range(1, h_outs[sel][h_out_name].GetNbinsX()+1)))
            if not (h_out_name.endswith('Up') or h_out_name.endswith('Down')):
                h_outs[sel][h_out_name].SetLineWidth(2)
            if DATA in h_out_name:
                h_outs[sel][h_out_name].SetLineWidth(3)
                h_outs[sel][h_out_name].SetLineColor(R.kBlack)
            elif '_bbq' in h_out_name or '2bq' in h_out_name or '2B' in h_out_name:
                h_outs[sel][h_out_name].SetLineColor(R.kRed)
            elif '_bb' in h_out_name or '2b' in h_out_name:
                h_outs[sel][h_out_name].SetLineColor(R.kOrange+8)
            elif '_bqq' in h_out_name:
                h_outs[sel][h_out_name].SetLineColor(R.kGreen+1)
            elif '1b' in h_out_name or '1B' in h_out_name:
                h_outs[sel][h_out_name].SetLineColor(R.kBlue)
            else:
                h_outs[sel][h_out_name].SetLineColor(R.kViolet)
            h_outs[sel][h_out_name].Write()

            # ## Add "extended" histograms failing selection cut
            # if sel == 'sel_JetID': continue
            # if 'bdtHi' in h_out_name: continue
            # ## Pick looser selection criteria and clone histogram
            # selX = ('sel_JetID' if sel == 'sel_Mass140' else 'sel_Mass140')
            # h_name_X = h_out_name.replace('bdtLo','ext')
            # h_out_X  = h_outs[selX][h_out_name].Clone(h_name_X)
            # ## Add in bdtHi histogram
            # h_out_X.Add(h_outs[sel][h_out_name.replace('bdtLo','bdtHi')])
            # ## Subtract tighter selection from looser selection
            # h_out_X.Add(h_outs[sel][h_out_name], -1)
            # h_out_X.Add(h_outs[sel][h_out_name.replace('bdtLo','bdtHi')], -1)
            # h_out_X.SetTitle(h_name_X)
            # h_out_X.Write()
        
        ## End loop: for h_out_name in h_outs[sel].keys()
    ## End loop: for sel in SELS:

    print('\n\nAll done!')
    
## End function: def main()


if __name__ == '__main__':
    main()
