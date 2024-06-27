#! /usr/bin/env python

## Script to estimate S/B ratio for various categories
import os
import sys
import math
import ROOT as R
from array import array

R.gROOT.SetBatch(True)  ## Don't display histograms or canvases when drawn
R.gStyle.SetOptStat(0)  ## Don't display stat boxes

## User configuration
VERBOSE  = False
YEAR     = '2018'
MASSH    = 'mass'
MASSWIN  = [100, 140]
MASSESA  = ['15', '30', '55']

# # IN_DIR  = '/eos/cms/store/user/ssawant/htoaa/analysis/20240530_ggH0l_for2DAlphabet/2018/2DAlphabet_inputFiles/'
# IN_DIR  = '/eos/cms/store/user/ssawant/htoaa/analysis/20240627_gg0l_1/2018/2DAlphabet_inputFiles/gg0lLo/'
# SIGS    = ['ggHtoaato4b']
# BKGS    = []
# CATS    = ['gg0lLo']
# WP_CUTS = ['WP40', 'WP60']

IN_DIR  = '/afs/cern.ch/user/m/moanwar/public/forAndrew/2DAlphabetfiles/'
SIGS    = ['VBFHtoaato4b']
BKGS    = []
CATS    = ['VBFjj']
WP_CUTS = ['WP40', 'WP60']

# IN_DIR   = '/afs/cern.ch/user/h/hboucham/public/2D_062324/' ## Single lepton
# SIGS     = ['WHtoaato4b', 'ttHtoaato4b']
# BKGS     = []  # ['Wlv', 'TT1l']
# CATS     = ['WlvLo', 'WlvHi', 'ttblv', 'ttbblv']
# WP_CUTS  = ['WP60']

# IN_DIR   = '/afs/cern.ch/user/h/hboucham/public/2D_Alphabet_root_052424/'
# SIGS     = ['ZHtoaato4b']
# BKGS     = []  # ['Wlv', 'TT1l']
# CATS     = ['Zll']
# WP_CUTS  = ['WP60']

# # IN_DIR  = '/eos/cms/store/user/ssawant/htoaa/analysis/20240530_VHHadronicMode_for2DAlphabet_1/2018/2DAlphabet_inputFiles/'
# IN_DIR  = '/eos/cms/store/user/ssawant/htoaa/analysis/20240626_Vjj/2018/2DAlphabet_inputFiles/Vjj/'
# SIGS    = ['WHtoaato4b', 'ZHtoaato4b']
# BKGS    = []
# CATS    = ['Vjj']
# WP_CUTS = ['WP40', 'WP60']

# # IN_DIR  = '/eos/cms/store/user/ssawant/htoaa/analysis/20240529_ZH_4b2nu_for2DAlphabet/2018/2DAlphabet_inputFiles/'
# IN_DIR  = '/eos/cms/store/user/ssawant/htoaa/analysis/20240626_Zvv_METDataset_2/2018/2DAlphabet_inputFiles/ZvvLo/'
# SIGS    = ['ZHtoaato4b']
# BKGS    = []
# CATS    = ['ZvvLo']
# WP_CUTS = ['WP40', 'WP60']


# OUT_FILE = 'plots/HtoAA_StoB_est.root'
# OUT_DIR  = 'plots/pdf/HtoAA_StoB_est/'


def main():

    print('\nInside HtoAA_StoB_est\n')

    samps = ['Data']
    for bkg in BKGS:
        samps.append(bkg)
    for sig in SIGS:
        for mA in MASSESA:
            samps.append(sig+'_mA_'+mA)

    yields = {}

    for wp in WP_CUTS:
        yields[wp] = {}
        for cat in CATS:
            yields[wp][cat] = {}
            for samp in samps:
                if cat.startswith('W') and samp.startswith('ttHtoaa'): continue
                if cat.startswith('ttb') and samp.startswith('WHtoaa'): continue
                yields[wp][cat][samp] = {}

                if 'hboucham' in IN_DIR:
                    in_file_str = IN_DIR+'%s/%s_%s_%s.root' % (wp, cat, samp, YEAR)
                else:
                    in_file_str = IN_DIR+'%s_%s_%s.root' % (cat, samp, YEAR)
                in_file = R.TFile(in_file_str, 'open')
                if VERBOSE: print('\nOpened %s' % in_file_str)

                for pf in ['Pass', 'Fail']:
                    if pf == 'Fail' and 'Htoaato4b' in samp: continue
                    
                    h_name = '%s_%s_%s_%s_%s_%s_Nom' % (cat, samp, YEAR, MASSH, wp, pf)
                    hist = in_file.Get(h_name)
                    if VERBOSE: print('\nGot histogram %s' % h_name)
                    if VERBOSE: print('  * Integral = %.1f' % hist.Integral())

                    yields[wp][cat][samp][pf] = {}
                    yields[wp][cat][samp][pf]['mH_in'] = 0

                    nY = hist.GetNbinsY()
                    for iX in range(1, hist.GetNbinsX()+1):
                        if hist.GetXaxis().GetBinCenter(iX) > MASSWIN[0] and \
                           hist.GetXaxis().GetBinCenter(iX) < MASSWIN[1]:
                            yields[wp][cat][samp][pf]['mH_in'] += hist.Integral(iX, iX, 1, nY)
                    yields[wp][cat][samp][pf]['mH_out'] = hist.Integral() - yields[wp][cat][samp][pf]['mH_in']

                    if VERBOSE: print('%s %s %s %s total = %.1f, mH_in = %.1f, mH_out = %.1f' % \
                                      ( wp, cat, samp, pf, hist.Integral(), \
                                        yields[wp][cat][samp][pf]['mH_in'], \
                                        yields[wp][cat][samp][pf]['mH_out'] ) )
                    del hist
                
                in_file.Close()
                
            ## End loop: for samp in samps

            nSig = 0
            for sig in SIGS:
                if cat.startswith('W') and sig.startswith('ttHtoaa'): continue
                if cat.startswith('ttb') and sig.startswith('WHtoaa'): continue
                for mA in MASSESA:
                    nSig += yields[wp][cat]['%s_mA_%s' % (sig, mA)]['Pass']['mH_in']
            nSig /= (10.0*len(MASSESA))  ## Scale down to 10% branching fraction
            if cat == 'gg0l':
                print('\nWARNING!!! In gg0l from 20240530, signal over-estimated by roughly a factor of 4!')
                print('Scaling signal by 0.25 to compensate. - AWB 2024.06.26')
                nSig *= 0.25

            nBkg  = yields[wp][cat]['Data']['Fail']['mH_in']
            nBkg *= yields[wp][cat]['Data']['Pass']['mH_out']
            nBkg /= yields[wp][cat]['Data']['Fail']['mH_out']

            eBkg  = (1.0 / yields[wp][cat]['Data']['Fail']['mH_in'])
            eBkg += (1.0 / yields[wp][cat]['Data']['Pass']['mH_out'])
            eBkg += (1.0 / yields[wp][cat]['Data']['Fail']['mH_out'])
            eBkg = nBkg * math.sqrt(eBkg)

            print('%s %s S/10 = %.3f, B = %.2f +/- %.2f, S/B = %.3f +/- %.3f, S/sqrt(B+err) = %.3f, limit = %.4f' % \
                  ( wp, cat, nSig, nBkg, eBkg, nSig/nBkg, (nSig/nBkg)*(eBkg/nBkg), \
                    nSig/math.sqrt(nBkg+eBkg), (0.2/nSig) * (1 + math.sqrt(nBkg+eBkg + 1)) ) )

        ## End loop: for cat in CATS
    ## End loop: for wp in WP_CUTS

    print('\nAll done!')
    
## End function: def main()


if __name__ == '__main__':
    main()
