#! /usr/bin/env python
## Script to generate input ROOT files for Combine fit of AK8 tagger efficiencies
## Before running, set up same CMSSW environment as Combine release:
## cd /afs/cern.ch/work/a/abrinke1/public/HiggsToAA/2DAlphabet/CMSSW_11_3_4/src/; cmssw-cc7
## cmsenv; cd /afs/cern.ch/work/a/abrinke1/public/HiggsToAA/coffea/eventloop/
## python3 macros/HtoAA_QCD_AK8_tagger_calib.py
## In Slack, see https://baylorhep.slack.com/archives/C013B0LRAEA/p1742318447152739?thread_ts=1741042320.395499&cid=C013B0LRAEA

import os
import sys
import math
import ROOT as R
import ctypes
from array import array

R.gROOT.SetBatch(True)  ## Don't display histograms or canvases when drawn
R.gStyle.SetOptStat(0)  ## Don't display stat boxes

## User configuration
VERBOSE  = False
IN_DIR   = '/eos/cms/store/user/ssawant/htoaa/analysis/20250317_CR_QCD4b_FullSyst/2018/CombineFit_inputFiles/'
CATS     = ['3M2T','3M3T','4M3T','4M4T']
#TAGGERS  = {'X4b_v2ab_Haa4b_score':[0.40,0.66,0.93]}  ## For WP60
TAGGERS  = {'X4b_v2ab_Haa4b_score':[0.40,0.84,0.96]}  ## For WP40
TAGNM    = {'X4b_v2ab_Haa4b_score':'X4b_v2'}
WP = 'WPX'
if TAGGERS['X4b_v2ab_Haa4b_score'][2] == 0.96: WP = 'WP40'
if TAGGERS['X4b_v2ab_Haa4b_score'][2] == 0.93: WP = 'WP60'

SVAR  = 2.0  ## Systematic factor of variation in tagging efficiency

YEAR     = '2018'
DATA     = 'Data'

MC_4b   = ['QCD_4bAndMoreCat']
MC_3b   = ['QCD_3bGen']
MC_012b = ['QCD_012bCat']
MC_TWZ  = ['TT0l','TT1l','Wqq','Zqq']
MC_SUMS = ['4B','3B','012B','TWZ']

MCNM = {'QCD_4bAndMoreCat':'QCD4b',
        'QCD_3bGen'       :'QCD3b',
        'QCD_012bCat'     :'QCD012b',
        'TT0l'            :'TT0l',
        'TT1l'            :'TT1l',
        'Wqq'             :'Wqq',
        'Zqq'             :'Zqq'}

OUT_DIR  = 'plots/HtoAA_QCD_AK8_tagger_calib/'


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

    ## Vary bins N,N-1,...2. When varying bin iN, adjust bins 1...iN-1 to keep total normalization.
    for iBin in range(h_out.GetNbinsX(), 1, -1):
        bStr = '_%dbin' % iBin
        h_out_up = h_out.Clone(h_out.GetName()+'_s'+syst+bStr+'Up')
        h_out_dn = h_out.Clone(h_out.GetName()+'_s'+syst+bStr+'Down')
        iN  = h_out.GetBinContent(iBin)
        iNe = h_out.GetBinError(iBin)
        oN  = h_out.Integral(1, iBin-1)  ## Integrate bins *below* iBin
        ## No shape variation if bin contains all or no events
        if iN == 0 or oN == 0:
            h_outs_up_dn.append(h_out_up)
            h_outs_up_dn.append(h_out_dn)
            continue
        iNup = min(iN*(SVAR-1.0), oN*(SVAR-1.0)/SVAR)  ## How much to add for up variation
        iNdn = min(iN*(SVAR-1.0)/SVAR, oN*(SVAR-1.0))  ## How much to subtract for down variation
        iFsy = min((iNup+iN)/iN, iN/(iN-iNdn))         ## Pick smaller of up/down for symmetric varation
        oFup = (oN - (iN*iFsy - iN)) / oN              ## Factor of reduction for other bins
        oFdn = (oN - (iN/iFsy - iN)) / oN              ## Factor of increase for other bins

        if VERBOSE: print('    - Bin %d %.2f --> +%.2f/-%.2f' % (iBin, iN, iN*(iFsy-1), iN*(1 - (1/iFsy))))
        if VERBOSE: print('    - Bin %d  iN = %.2f, oN = %.2f, iFsy = %.3f, oFup = %.3f, oFdn = %.3f.' % (iBin, iN, oN, iFsy, oFup, oFdn))
        if (iFsy <= 0 or oFup <= 0 or oFdn <= 0):
            print('\n\nBAD ERROR! Bin %d iFsy = %.3f, oFup = %.3f, oFdn = %.3f.\n\n' % (iBin, iFsy, oFup, oFdn))
            sys.exit()
        h_out_up.SetBinContent(iBin, iN*iFsy)
        h_out_dn.SetBinContent(iBin, iN/iFsy)
        h_out_up.SetBinError(iBin, iNe*iFsy)
        h_out_dn.SetBinError(iBin, iNe/iFsy)
        for xBin in range(1, iBin):
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

    print('\nInside HtoAA_QCD_AK8_tagger_calib\n')

    if not os.path.exists(OUT_DIR):
        os.mkdir(OUT_DIR)

    ## Initialize dictionary of output histograms
    h_outs = {}
        
    for cat in CATS:
        if VERBOSE: print('\nNow looking at category %s' % cat)

        for tag in TAGGERS.keys():
            nCuts = len(TAGGERS[tag])
            if VERBOSE: print('\nNow looking at tagger %s, will rebin to %d bins' % (tag, nCuts+1))

            ######################
            ## Get data histograms
            if VERBOSE: print('\nNow looking at data')
            in_file_str = IN_DIR+'CR4b_'+cat+'/CR4b_'+cat+'_'+DATA+'_'+YEAR+'.root'
            in_file = R.TFile(in_file_str, 'open')
            print('\n*******\nReading from %s' % in_file_str)
            h_in_name = 'CR4b_'+cat+'_'+DATA+'_'+YEAR+'_'+tag+'_Nom'
            h_in = in_file.Get(h_in_name)
            if VERBOSE: print('\nGot histogram %s' % h_in_name)
            if VERBOSE: print('  * Integral = %.1f' % h_in.Integral())

            ## Fill new rebinned histogram with events in tagger ranges
            h_out_name = DATA+'_'+YEAR+'_'+cat+'_'+TAGNM[tag]
            if not h_out_name in h_outs.keys():
                h_outs[h_out_name] = R.TH1D(h_out_name, h_out_name, nCuts+1, 0, nCuts+1)
                h_outs[h_out_name].SetDirectory(0) ## Save locally
            fill_pass_fail(h_in, h_outs[h_out_name], TAGGERS[tag])
            in_file.Close()

            ####################
            ## Get MC histograms
            for mc in MC_4b+MC_3b+MC_012b+MC_TWZ:
                if VERBOSE: print('\nNow looking at MC sample %s' % mc)
                in_file_str = IN_DIR+'CR4b_'+cat+'/CR4b_'+cat+'_'+mc+'_'+YEAR+'.root'
                in_file = R.TFile(in_file_str, 'open')
                h_in_name = 'CR4b_'+cat+'_'+mc+'_'+YEAR+'_'+tag+'_Nom'
                h_in = in_file.Get(h_in_name)
                if VERBOSE: print('\nGot histogram %s' % h_in_name)
                if VERBOSE: print('  * Integral = %.1f' % h_in.Integral())

                ## Fill new rebinned histogram with events failing and passing cuts
                mcm = MCNM[mc]
                h_out_name = mcm+'_'+cat+'_'+TAGNM[tag]
                if not h_out_name in h_outs.keys():
                    h_outs[h_out_name] = R.TH1D(h_out_name, h_out_name, nCuts+1, 0, nCuts+1)
                    h_outs[h_out_name].SetDirectory(0) ## Save locally
                else:
                    print('\n\nHow is %s already in h_outs??? Quitting.' % h_out_name)
                    sys.exit()
                fill_pass_fail(h_in, h_outs[h_out_name], TAGGERS[tag])
                    
                if mc in MC_4b:   systs = ['4B']
                if mc in MC_3b:   systs = ['3B']
                if mc in MC_012b: systs = ['012B']
                if mc in MC_TWZ:  systs = ['TWZ']
                for syst in systs:
                    ## Generate additional histograms with sum of MC
                    h_MC_name = 'Sum'+syst+'_'+cat+'_'+TAGNM[tag]
                    if not h_MC_name in h_outs.keys():
                        h_outs[h_MC_name] = R.TH1D(h_MC_name, h_MC_name, nCuts+1, 0, nCuts+1)
                        h_outs[h_MC_name].SetDirectory(0) ## Save locally
                    h_outs[h_MC_name].Add(h_outs[h_out_name])

                    ## Perform systematic variations
                    for h_syst in make_syst_hists(h_outs[h_out_name], syst):
                        h_outs[h_syst.GetName()] = h_syst
                        h_outs[h_syst.GetName()].SetDirectory(0) ## Save locally
                        ## Generate additional systematic histograms with sum of MC
                        h_MC_name_syst = h_MC_name+(h_syst.GetName().replace(h_out_name,''))
                        if not h_MC_name_syst in h_outs.keys():
                            h_outs[h_MC_name_syst] = R.TH1D(h_MC_name_syst, h_MC_name_syst, nCuts+1, 0, nCuts+1)
                            h_outs[h_MC_name_syst].SetDirectory(0) ## Save locally
                        h_outs[h_MC_name_syst].Add(h_syst)
                    ## End loop: for h_syst in make_syst_hists(h_outs[h_out_name], syst)
                ## End loop: for syst in systs
                in_file.Close()
            ## End loop: for mc in MC_4b+MC_3b+MC_012b+MC_TWZ
                
        ## End loop: for tag in TAGGERS.keys()
    ## End loop: for cat in CAT

    
    ## Create summed histograms across 3M3T, 4M3T, and 4M4T categories,
    ##   then use as shape templates for 4M3T and 4M4T categories
    h_summs = {}
    ## Loop over taggers
    for tag in TAGGERS.keys():
        ## Loop over summed MC samples
        for mcSum in MC_SUMS:
            base = 'Sum'+mcSum+'_3M3T_'+TAGNM[tag]
            ## Loop over histograms
            for h_out_name in h_outs.keys():
                if not h_out_name.startswith(base): continue
                suff = h_out_name.replace(base, '')
                summ = base.replace('3M3T','34T')+suff
                ## Sanity check: summ should not already exist
                if summ in h_summs.keys():
                    print('\n\nWEIRD ERROR!!! %s already exists, but now revisiting. Will quit instead.\n\n' % summ)
                    sys.exit()
                h_summs[summ] = h_outs[h_out_name].Clone(summ)
                n_4M3T = h_out_name.replace('3M3T','4M3T')
                n_4M4T = h_out_name.replace('3M3T','4M4T')
                h_summs[summ].Add(h_outs[n_4M3T])
                h_summs[summ].Add(h_outs[n_4M4T])
                yields = {'3M3T': h_outs[base].Integral(),
                          '4M3T': h_outs[base.replace('3M3T','4M3T')].Integral(),
                          '4M4T': h_outs[base.replace('3M3T','4M4T')].Integral(),
                          '34T':  h_summs[summ].Integral()}
                ## Sanity check: 34T should add up to sum of other categories
                diff = yields['34T'] - yields['4M4T'] - yields['4M3T'] - yields['3M3T']
                if abs(diff) > 0.01:
                    print('\n\nWEIRD ERROR!!! For %s, 34T = %.4f (%.f higher than sum). Quitting.' % (h_out_name, yields['34T'], diff))
                    for yld in yields.keys():
                        print('%s = %.4f' % (yld, yields[yld]))
                    sys.exit()
                ## Set 4M3T and 4M4T shapes to equal 34T shape
                scale_4M3T = yields['4M3T'] / yields['34T']
                scale_4M4T = yields['4M4T'] / yields['34T']
                for iBin in range(1, h_summs[summ].GetNbinsX()+1):
                    h_outs[n_4M3T].SetBinContent(iBin, h_summs[summ].GetBinContent(iBin)*scale_4M3T)
                    h_outs[n_4M3T].SetBinError  (iBin, h_summs[summ].GetBinError  (iBin)*scale_4M3T)
                    h_outs[n_4M4T].SetBinContent(iBin, h_summs[summ].GetBinContent(iBin)*scale_4M4T)
                    h_outs[n_4M4T].SetBinError  (iBin, h_summs[summ].GetBinError  (iBin)*scale_4M4T)
                ## End loop: for iBin in range(1, h_summs[summ].GetNbinsX()+1)
                h_summs[summ].SetDirectory(0) ## Save locally
                h_outs[n_4M3T].SetDirectory(0) ## Save locally
                h_outs[n_4M4T].SetDirectory(0) ## Save locally
    
            ## End loop: for h_out_name in h_outs.keys()
        ## End loop: for mcSum in MC_SUMS
    ## End loop: for tag in TAGGERS.keys()
    h_outs.update(h_summs)  ## Merge two dictionaries of histograms


    ## Create an output ROOT file
    tag_str = '%s'.join(TAGNM[tag] for tag in TAGGERS.keys())
    out_file_str = OUT_DIR+'AK8_tagger_calib_%s_%s_%s_slc7.root' % (tag_str, WP, str(SVAR).replace('.','p'))
    out_file = R.TFile(out_file_str, 'recreate')
    print('\n*******\nWriting to %s' % out_file_str)
    for h_out_name in h_outs.keys():
        if VERBOSE: print('Writing out %s' % h_out_name)
        if VERBOSE: print('  * Integral = %.3f' % h_outs[h_out_name].Integral())
        if VERBOSE: print('  * Bins = '+', '.join('%.3f' % h_outs[h_out_name].GetBinContent(iX)
                                                  for iX in range(1, h_outs[h_out_name].GetNbinsX()+1)))
        if not (h_out_name.endswith('Up') or h_out_name.endswith('Down')):
            h_outs[h_out_name].SetLineWidth(2)
        if DATA in h_out_name:
            h_outs[h_out_name].SetLineWidth(3)
            h_outs[h_out_name].SetLineColor(R.kBlack)
        elif 'QCD4b' in h_out_name or 'Sum4B' in h_out_name:
            h_outs[h_out_name].SetLineColor(R.kRed)
        elif 'QCD3b' in h_out_name or 'Sum3B' in h_out_name:
            h_outs[h_out_name].SetLineColor(R.kGreen+1)
        elif 'QCD012b' in h_out_name or 'Sum012B' in h_out_name:
            h_outs[h_out_name].SetLineColor(R.kBlue)
        else:
            h_outs[h_out_name].SetLineColor(R.kViolet)
        h_outs[h_out_name].Write()

    ## End loop: for h_out_name in h_outs.keys()

    print('\n\nAll done!')
    
## End function: def main()


if __name__ == '__main__':
    main()
