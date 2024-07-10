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
MASSESH  = ['mass','msoft','pnet']
MASSESA  = ['15', '30', '55']
WP_CUTS  = ['WP40', 'WP60', 'WP80']
DATE     = '2024_07_09'

# BKGS = ['Wlv', 'TT1l']

# CAT_OUT = 'LepHi'
# CAT_INS = ['WlvHi', 'ttbblv', 'Zll']

CAT_OUT = 'XXHi'
CAT_INS = ['WlvHi', 'ttbblv', 'ttbll', 'ttbbll', 'Zll', 'ZvvHi']
CATS_IN = {}
for cat in CAT_INS:
    CATS_IN[cat] = {}
CATS_IN['WlvHi']['dir']   = '/afs/cern.ch/user/h/hboucham/public/2D_1L_070924/'
CATS_IN['WlvHi']['sigs']  = ['WHtoaato4b']
CATS_IN['ttbblv']['dir']  = '/afs/cern.ch/user/h/hboucham/public/2D_1L_070924/'
CATS_IN['ttbblv']['sigs'] = ['ttHtoaato4b']
CATS_IN['ttbll']['dir']   = '/afs/cern.ch/user/h/hboucham/public/2D_2Ltt_070924/'
CATS_IN['ttbll']['sigs']  = ['ttHtoaato4b']
CATS_IN['ttbbll']['dir']  = '/afs/cern.ch/user/h/hboucham/public/2D_2Ltt_070924/'
CATS_IN['ttbbll']['sigs'] = ['ttHtoaato4b']
CATS_IN['Zll']['dir']     = '/afs/cern.ch/user/h/hboucham/public/2D_2LZ_070924/'
CATS_IN['Zll']['sigs']    = ['ZHtoaato4b']
CATS_IN['ZvvHi']['dir']   = '/eos/cms/store/user/ssawant/htoaa/analysis/20240626_Zvv_METDataset_2/2018/2DAlphabet_inputFiles/ZvvHi/'
CATS_IN['ZvvHi']['sigs']  = ['ZHtoaato4b']

# CAT_OUT = 'XXLo'
# CAT_INS = ['WlvLo', 'ttblv', 'ZvvLo']
# CATS_IN = {}
# for cat in CAT_INS:
#     CATS_IN[cat] = {}
# CATS_IN['WlvLo']['dir']   = '/afs/cern.ch/user/h/hboucham/public/2D_1L_070924/'
# CATS_IN['WlvLo']['sigs']  = ['WHtoaato4b']
# CATS_IN['ttblv']['dir']  = '/afs/cern.ch/user/h/hboucham/public/2D_1L_070924/'
# CATS_IN['ttblv']['sigs'] = ['ttHtoaato4b']
# CATS_IN['ZvvLo']['dir']   = '/eos/cms/store/user/ssawant/htoaa/analysis/20240626_Zvv_METDataset_2/2018/2DAlphabet_inputFiles/ZvvLo/'
# CATS_IN['ZvvLo']['sigs']  = ['ZHtoaato4b']


OUT_DIR  = 'plots/HtoAA_2DAlphabet_merge_inputs_%s/%s/' % (CAT_OUT, DATE)

def main():

    print('\nInside HtoAA_2DAlphabet_merge_inputs\n')

    print('\n\nWARNING!!! There is a hacky conditional that will fail if')
    print('signal samples are shared between categories!!! - AWB 2024.06.24\n')
    print('\n\nWARNING!!! Should run from inside 2DAlphabet/CMSSW_11_3_4/src/2DAlphabet/')
    print('to avoid "list is accessing an object already deleted" error! - AWB 2024.06.24\n')
    print('See https://root-forum.cern.ch/t/error-in-tlist-clear-a-list-is-accessing-an-object-already-deleted-list-name-tlist-when-opening-a-file-created-by-root-6-30-using-root-6-14-09/57588/1')
    

    if not os.path.exists(OUT_DIR):
        os.mkdir(OUT_DIR)

    h_outs = {}
    for cat in CAT_INS:
        samps = ['Data']
        for sig in CATS_IN[cat]['sigs']:
            for mA in MASSESA:
                samps.append(sig+'_mA_'+mA)
        for samp in samps:
            for wp in WP_CUTS:
                if 'hboucham' in CATS_IN[cat]['dir']:
                    in_file_str = CATS_IN[cat]['dir']+'%s/%s_%s_%s.root' % (wp, cat, samp, YEAR)
                else:
                    in_file_str = CATS_IN[cat]['dir']+'%s_%s_%s.root' % (cat, samp, YEAR)
                in_file = R.TFile(in_file_str, 'open')
                print('\n*******\nReading from %s' % in_file_str)

                for mH in MASSESH:
                    for pf in ['Pass', 'Fail']:
                        h_in_name = '%s_%s_%s_%s_%s_%s_Nom' % (cat, samp, YEAR, mH, wp, pf)
                        h_in = in_file.Get(h_in_name)
                        if VERBOSE: print('\nGot histogram %s' % h_in_name)
                        if VERBOSE: print('  * Integral = %.1f' % h_in.Integral())

                        h_out_name = '%s_%s_%s_%s_%s_%s_Nom' % (CAT_OUT, samp, YEAR, mH, wp, pf)
                        if not h_out_name in h_outs.keys():
                            h_outs[h_out_name] = h_in.Clone(h_out_name)
                            if VERBOSE: print('Created %s' % h_out_name)
                            if VERBOSE: print('  * Integral = %.1f' % h_outs[h_out_name].Integral())
                            h_outs[h_out_name].SetDirectory(0) ## Save locally
                        else:
                            if VERBOSE: print('Adding %d to %s (with %d)' % (h_in.Integral(), h_out_name, \
                                                                             h_outs[h_out_name].Integral()))
                            nXo = h_outs[h_out_name].GetNbinsX()
                            nYo = h_outs[h_out_name].GetNbinsY()
                            nXi = h_in.GetNbinsX()
                            nYi = h_in.GetNbinsY()
                            xLo = h_outs[h_out_name].GetXaxis().GetBinLowEdge(1)
                            xHo = h_outs[h_out_name].GetXaxis().GetBinLowEdge(nXo+1)
                            yLo = h_outs[h_out_name].GetYaxis().GetBinLowEdge(1)
                            yHo = h_outs[h_out_name].GetYaxis().GetBinLowEdge(nYo+1)
                            xLi = h_in.GetXaxis().GetBinLowEdge(1)
                            xHi = h_in.GetXaxis().GetBinLowEdge(nXo+1)
                            yLi = h_in.GetYaxis().GetBinLowEdge(1)
                            yHi = h_in.GetYaxis().GetBinLowEdge(nYo+1)
                            if (nXo != nXi or nYo != nYi or xLo != xLi or xHo != xHi or yLo != yLi or yHo != yHi):
                                print('\nMAJOR ERROR!!! %s is %d x %d, %s is %d x %d' % (h_out_name, nXo, nYo,
                                                                                         nYi, nYi, h_in.GetName()))
                                print('Spanning [%.1f-%.1f] x [%.1f-%.1f] vs. [%.1f-%.1f]' % (xLo, xHo, yLo, yHo,
                                                                                              xLi, xHi,	yLi, yHi))
                                sys.exit()
                            else:    
                                h_outs[h_out_name].Add(h_in)
                                if VERBOSE: print('  * Integral = %.1f' % h_outs[h_out_name].Integral())

                    ## End loop: for pf in ['Pass', 'Fail']
                ## End loop: for mH in MASSESH
                in_file.Close()

                ## Hacky conditional that will fail if signal samples are shared between categories!!! - AWB 2024.06.24
                if cat == CAT_INS[-1] or samp.split('_mA_')[0] in CATS_IN[cat]['sigs']:
                    out_file_str = OUT_DIR+'%s_%s_%s.root' % (CAT_OUT, samp, YEAR)
                    if wp == WP_CUTS[0]:
                        out_file = R.TFile(out_file_str, 'recreate')
                    else:
                        out_file = R.TFile(out_file_str, 'update')
                    print('\n*******\nWriting to %s' % out_file_str)
                    for mH in MASSESH:
                        for pf in ['Pass', 'Fail']:
                            h_out_name = '%s_%s_%s_%s_%s_%s_Nom' % (CAT_OUT, samp, YEAR, mH, wp, pf)
                            h_outs[h_out_name].Write()
                            if VERBOSE: print('Wrote out %s' % h_out_name)
                            if VERBOSE: print('  * Integral = %.1f' % h_outs[h_out_name].Integral())
                        ## End loop: for pf in ['Pass', 'Fail']
                    ## End loop: for mH in MASSESH
                    out_file.Write()
                    out_file.Close()
            ## End loop: for wp in WP_CUTS
        ## End loop: for samp in samps
    ## End loop: for cat in CAT_INS

    print('\n\nAll done!')
    
## End function: def main()


if __name__ == '__main__':
    main()
