#! /usr/bin/env python
## Script to compute efficiency at different tagger working points
## python3 macros/HtoAA_TTBarLep_AK8_tagger_eff.py

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
DSETS    = ['SingleMuon','EGamma','SingleLep']
IN_DIR   = '/afs/cern.ch/work/c/csutanta/public/rootfiles/v1/lepjet/unskimmed_'
CATS     = ['bdtHi']
SELS     = ['0b_BBQQ']
TAGGERS  = {'Xbb':'particleNetMD_XbbOverQCD'}
            #'Hto34b':'FatJet_PNetMD_Hto4b_Htoaa34bOverQCD'}

MC_2B = ['TT1L_bbqq','TT1L_bbq','TT1L_bb','TT2L_bb','Sum2B']
samps = {'TT1L_bbqq': 'TTToSemiLeptonic_powheg_bbqq',
         'TT1L_bbq':  'TTToSemiLeptonic_powheg_bbq',
         'TT1L_bb':   'TTToSemiLeptonic_powheg_bb',
         'TT2L_bb':   'TTTo2L2Nu_powheg_bb',
         'Sum2B':     'Sum2B'}
colors = {'TT1L_bbqq': R.kGreen+1,
          'TT1L_bbq':  R.kBlue,
          'TT1L_bb':   R.kRed,
          'TT2L_bb':   R.kMagenta,
          'Sum2B':     R.kBlack}

OUT_DIR  = 'plots/HtoAA_TTBarLep_AK8_tagger_eff/'


def main():

    print('\nInside HtoAA_TTBarLep_AK8_tagger_eff\n')

    if not os.path.exists(OUT_DIR):
        os.mkdir(OUT_DIR)

    ## Histograms and yields by category, dataset, selection, tagger
    hst = {}
    evt = {}
        
    for cat in CATS:
        hst[cat] = {}
        evt[cat] = {}
        for dset in DSETS:
            hst[cat][dset] = {}
            evt[cat][dset] = {}
            if VERBOSE: print('\nNow looking at category %s, dataset %s' % (cat, dset))
            in_file_str = IN_DIR+dset+'_'+cat+'.root'
            if not dset == 'SingleLep':
                in_file = R.TFile(in_file_str, 'open')
                print('\n*******\nReading from %s' % in_file_str)
            ## End conditional: if not dset == 'SingleLep'

            for sel in SELS:
                if VERBOSE: print('\nNow looking at selection %s' % sel)
                hst[cat][dset][sel] = {}
                evt[cat][dset][sel] = {}
                for tag in TAGGERS.keys():
                    if VERBOSE: print('\nNow looking at tagger %s' % tag)
                    hst[cat][dset][sel][tag] = {}
                    evt[cat][dset][sel][tag] = {}
                    for mc in MC_2B:
                        if VERBOSE: print('\nNow looking at MC %s (%s)' % (mc, samps[mc]))
                        if not dset == 'SingleLep' and not mc == 'Sum2B':
                            h_in_name = 'evt/'+samps[mc]+'/'+TAGGERS[tag]+'_sel_'+sel+'_central'
                            h_in = in_file.Get(h_in_name)
                            if VERBOSE: print('%s %s %s %s %s = %.2f' % (cat, dset, sel, tag, mc, h_in.Integral()))
                            hst[cat][dset][sel][tag][mc] = h_in.Clone()
                            hst[cat][dset][sel][tag][mc].SetDirectory(0)
                            hst[cat][dset][sel][tag][mc].SetName('%s_%s_%s_%s_%s' % (cat, dset, sel, tag, mc))
                            evt[cat][dset][sel][tag][mc] = [h_in.Integral(), h_in.Integral(151,200)]
                            if not 'Sum2B' in hst[cat][dset][sel][tag].keys():
                                hst[cat][dset][sel][tag]['Sum2B'] = h_in.Clone()
                                hst[cat][dset][sel][tag]['Sum2B'].SetDirectory(0)
                                hst[cat][dset][sel][tag]['Sum2B'].SetName('%s_%s_%s_%s_%s' % (cat, dset, sel, tag, 'Sum2B'))
                                evt[cat][dset][sel][tag]['Sum2B'] = [h_in.Integral(), h_in.Integral(151,200)]
                            else:
                                hst[cat][dset][sel][tag]['Sum2B'].Add(h_in)
                                evt[cat][dset][sel][tag]['Sum2B'][0] += h_in.Integral()
                                evt[cat][dset][sel][tag]['Sum2B'][1] += h_in.Integral(151,200)
                        ## End conditional: if not dset == 'SingleLep' and not mc == 'Sum2B'
                        if dset == 'SingleLep':
                            hst[cat][dset][sel][tag][mc] = hst[cat]['SingleMuon'][sel][tag][mc].Clone()
                            hst[cat][dset][sel][tag][mc].SetDirectory(0)
                            hst[cat][dset][sel][tag][mc].SetName('%s_%s_%s_%s_%s' % (cat, dset, sel, tag, mc))
                            hst[cat][dset][sel][tag][mc].Add(hst[cat]['EGamma'][sel][tag][mc])
                            evt[cat][dset][sel][tag][mc] = [hst[cat][dset][sel][tag][mc].Integral(),
                                                            hst[cat][dset][sel][tag][mc].Integral(151,200)]
                        ## End conditional: if dset == 'SingleLep'
                            
                    ## End loop: for mc in MC_2B

                    if dset == 'SingleLep' or VERBOSE:
                        print('\nYields and efficiencies for different cut values:')
                        print('Category = %s, dataset = %s, selection = %s, tagger = %s, cut = %.3f' % (cat, dset, sel, tag, h_in.GetBinLowEdge(151)))
                        numEvt = evt[cat][dset][sel][tag]['Sum2B'][1]
                        denEvt = evt[cat][dset][sel][tag]['Sum2B'][0]
                        print('Nominal efficiency = %5.1f / %5.1f = %.2f%%' % (numEvt, denEvt, 100.0*numEvt/denEvt))
                        print('------------------------------------------------')
                        den = hst[cat][dset][sel][tag]['Sum2B'].Integral()
                        for ii in range(100,201):
                            num = hst[cat][dset][sel][tag]['Sum2B'].Integral(ii,200)
                            cut = hst[cat][dset][sel][tag]['Sum2B'].GetBinLowEdge(ii)
                            print('Eff. at %.3f = %5.1f / %5.1f = %.2f%%' % (cut, num, den, 100*num/den))
                        print('------------------------------------------------\n')
                    ## End conditional: if dset == 'SingleLep'

                ## End loop: for tag in TAGGERS.keys()
            ## End loop: for sel in SELS
            if not dset == 'SingleLep':
                in_file.Close()
        ## End loop: for dset in DSETS
    ## End loop: for cat in CATS

    ## Create a separate output ROOT file for each selection option
    for cat in CATS:
        tag_str = '%s'.join(TAGGERS[tag] for tag in TAGGERS.keys())
        out_file_str = OUT_DIR+'AK8_tagger_eff_SingleLep_%s_%s.root' % (tag_str, cat)
        out_file = R.TFile(out_file_str, 'recreate')
        out_file.cd()
        print('\n*******\nWriting to %s' % out_file_str)
        for sel in SELS:
            for tag in TAGGERS.keys():
                can = R.TCanvas('can_%s_%s_%s_%s' % (cat, sel, tag, mc))
                can.cd()
                leg = R.TLegend(0.12,0.58,0.36,0.88)
                for mc in reversed(MC_2B):
                    hst[cat]['SingleLep'][sel][tag][mc].SetLineWidth(3 if mc == 'Sum2B' else 2)
                    hst[cat]['SingleLep'][sel][tag][mc].SetLineColor(colors[mc])
                    hst[cat]['SingleLep'][sel][tag][mc].Rebin(10)
                    hst[cat]['SingleLep'][sel][tag][mc].SetTitle('%s %s %s %s' % (cat, sel, tag, mc))
                    hst[cat]['SingleLep'][sel][tag][mc].GetXaxis().SetTitle('PNet %s score' % tag)
                    hst[cat]['SingleLep'][sel][tag][mc].Write()
                    hst[cat]['SingleLep'][sel][tag][mc].Draw('hist' if mc == 'Sum2B' else 'histsame')
                    leg.AddEntry(hst[cat]['SingleLep'][sel][tag][mc], mc)
                ## End loop: for mc in MC_2B.reverse()
                leg.Draw()
                can.Write()
                can.SaveAs(out_file_str.replace('.root','.pdf'))
                can.SaveAs(out_file_str.replace('.root','.png'))
                del can
            ## End loop: for tag in TAGGERS.keys()
        ## End loop: for sel in SELS
        out_file.Write()
        out_file.Close()

    print('\n\nAll done!')
    
## End function: def main()


if __name__ == '__main__':
    main()
