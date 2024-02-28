#! /usr/bin/env python

## Script to compute mass scale for different ParticleNet "a" boson mass regressions
import os
import sys
import math
import ROOT as R
from array import array

R.gROOT.SetBatch(True)  ## Don't display histograms or canvases when drawn

## User configuration
VERBOSE  = False
REBIN    = True
MASSES   = [12, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]
#ALGOS    = ['v0','v1','v2','v3']
#ALGOS    = ['avg_v01','avg_v02','avg_v03','avg_v12','avg_v13','avg_v23']
#ALGOS    = ['avg_v012','avg_v013','avg_v023','avg_v123','avg_v0123']
ALGOS     = ['v0','avg_v01','avg_v03','avg_v013','avg_v0123']
#MASSES   = [60]
#ALGOS    = ['v1']
WP_CUT   = 'WP80'  ## WP40, WP60, or WP80
IN_DIR   = '/eos/cms/store/user/ssawant/htoaa/analysis/'
#IN_FILE  = IN_DIR+'20231110_Signal_PNetHtoaa4bOverQCDWP40_woTrg/2018/analyze_htoaa_stage1.root' ## Wider binning
#IN_FILE  = IN_DIR+'20231110_Signal_PNetHtoaa4bOverQCDWP40_woTrg_1/2018/analyze_htoaa_stage1.root' ## Finer binning
IN_FILE  = IN_DIR+'20231117_Signal_woTrg/2018/analyze_htoaa_stage1.root' ## Include combined algos
OUT_FILE = 'plots/HtoAA_massA_scale_%s.root' % WP_CUT
OUT_DIR  = 'plots/pdf/HtoAA_massA_scale_%s/' % WP_CUT
SEL_CUTS = [0.1, 0.3, 0.5]  ## Height of bin w.r.t. central bin to be included in peak-finding


def find_peak(hstX, mass, cut):
    hmax = hstX.GetMaximum()
    hint = hstX.Integral()
    hbin = hstX.GetBinWidth(1)
    bins = []
    cent = []
    evts = []
    evt_sq  = 0
    mass_sq = 0
    diff_sq = 0
    inv_evt_sq = 0


    ## Find first and last bins passing threshold
    for i in range(1,hstX.GetNbinsX()+1):
        if hstX.GetBinContent(i) > cut*hmax:
            bins.append(i)

    ## For all bins in peak, compute peak and width
    for i in range(min(bins), max(bins)+1):
        center = hstX.GetBinCenter(i)
        events = hstX.GetBinContent(i)
        cent.append(center)
        evts.append(events)
        evt_sq  += pow(events, 2)
        mass_sq += center*pow(events, 2)
        inv_evt_sq += pow(events/hint,2)
    peak_val   = mass_sq/evt_sq
    peak_width = hbin/inv_evt_sq
    peak_bins  = max(bins) - min(bins) + 1

    if VERBOSE:
        print(bins)
        print(cent)
        print(evts)
        print('For M-%d cut %.2f, peak found at %.4f' % (mass, cut, peak_val))

    return [peak_val, peak_width, peak_bins]

## end function: find_peak(hstX, mass, cut)


def main():

    print('\nInside HtoAA_massA_scale\n')

    print('\nGetting input file %s' % IN_FILE)
    infile = R.TFile(IN_FILE,'open')
    outfile = R.TFile(OUT_FILE,'recreate')
    if not os.path.exists(OUT_DIR):
        os.makedirs(OUT_DIR)

    hist_beg = 'evt/SUSY_GluGluH_01J_HToAATo4B_M-'
    hist_mid = '_HPtAbv150/hLeadingFatJetParticleNet_massA_Hto4b_'
    # hist_end = '_SR_central'  ## Initial file (20231110)
    hist_end = '_SR%s_central' % WP_CUT
    hist = {}

    ## Store peak and with values
    masses_val = {}
    masses_err = {}
    peaks_val  = {}
    peaks_hi   = {}
    peaks_lo   = {}
    widths_mass = {}
    widths_val  = {}
    widths_hi   = {}
    widths_lo   = {}
    for algo in ALGOS:
        masses_val[algo] = array('d')
        masses_err[algo] = array('d')
        peaks_val [algo] = array('d')
        peaks_hi  [algo] = array('d')
        peaks_lo  [algo] = array('d')
        widths_val[algo] = array('d')
        widths_hi [algo] = array('d')
        widths_lo [algo] = array('d')

    ## Loop over mass points
    for mass in MASSES:
        if VERBOSE:
            print('\nStarting to look at M-%d' % mass)

        ## Loop over algos
        for algo in ALGOS:
            if VERBOSE:
                print('\nGetting histogram '+hist_beg+str(mass)+hist_mid+algo+hist_end)
            hst = infile.Get(hist_beg+str(mass)+hist_mid+algo+hist_end)
            if VERBOSE:
                print('\nM-%d %s has integral %.1f, peak %.2f' % (mass, algo, hst.Integral(), hst.GetMaximum()))

            ## Rebin for larger "a" boson masses, with wider distributions
            rebin = int(math.floor(mass/3)) if (REBIN and mass != 12) else 1
            nbins = hst.GetNbinsX()
            drop  = (nbins % rebin)
            if VERBOSE or True:
                print('Rebinning M-%d by a factor of %d' % (mass, rebin))

            ## Drop last bins to get an even division for the rebinning
            print('Cloning M-%d histogram and dropping %d bins' % (mass, drop))
            hst_new = R.TH1D(hst.GetName()+'_copy', hst.GetTitle()+'_copy', nbins-drop,
                             hst.GetBinLowEdge(1), hst.GetBinLowEdge(nbins-drop+1))
            for ii in range(1,nbins-drop+1):
                hst_new.SetBinContent(ii, hst.GetBinContent(ii))
                hst_new.SetBinError  (ii, hst.GetBinError(ii))
            ## Sanity check that the histograms are identical
            for ii in range(1, hst_new.GetNbinsX()):
                if abs(hst_new.GetBinContent(ii) - hst.GetBinContent(ii)) > 0.0001 or \
                   abs(hst_new.GetBinError(ii)   - hst.GetBinError(ii))   > 0.0001 or \
                   abs(hst_new.GetBinLowEdge(ii) - hst.GetBinLowEdge(ii)) > 0.001:
                    print('Mismatch in bin %d!!! Quitting.' % ii)
                    print(hst_new.GetBinContent(ii), hst.GetBinContent(ii))
                    print(hst_new.GetBinError(ii), hst.GetBinError(ii))
                    print(hst_new.GetBinLowEdge(ii), hst.GetBinLowEdge(ii))
                    sys.exit()

            ## Rebin and find the peak!
            hst_new.Rebin(rebin)
            peaks = []
            widths = []
            for cut in SEL_CUTS:
                peak = find_peak(hst_new, mass, cut)
                print('For M-%d %s cut %.2f, peak found at %.4f (width = %.5f, %d bins)' % (mass, algo, cut, peak[0], peak[1], peak[2]))
                #print('For M-%d %s cut %.2f, peak shift %.2f%% (width = %.2f%%, %d bins)' % (mass, algo, cut, 100*((peak[0]/mass) - 1), 100*peak[1]/mass, peak[2]))
                peaks.append(100*((peak[0]/mass) - 1))
                widths.append(100*peak[1]/mass)


            ## Save values for TGraphAsymmErrors
            masses_val[algo].append(mass)
            masses_err[algo].append(0)
            peaks_val [algo].append(sum(peaks)/len(peaks))
            peaks_hi  [algo].append(max(peaks) - (sum(peaks)/len(peaks)))
            peaks_lo  [algo].append((sum(peaks)/len(peaks)) - min(peaks))
            widths_val [algo].append(sum(widths)/len(widths))
            widths_hi  [algo].append(max(widths) - (sum(widths)/len(widths)))
            widths_lo  [algo].append((sum(widths)/len(widths)) - min(widths))

            del hst
            del hst_new

        ## End loop: for algo in ALGOS
    ## End loop: for mass in MASSES


    ## Generate TGraphAsymmErrors
    outfile.cd()
    grp = {}
    color = 2
    for algo in ALGOS:
        grp[algo+'_peak'] = R.TGraphAsymmErrors(len(MASSES), masses_val[algo], peaks_val[algo],
                                                masses_err[algo], masses_err[algo],
                                                peaks_lo[algo], peaks_hi[algo])
        grp[algo+'_width'] = R.TGraphAsymmErrors(len(MASSES), masses_val[algo], widths_val[algo],
                                                 masses_err[algo], masses_err[algo],
                                                 widths_lo[algo], widths_hi[algo])

        # color = R.kViolet  if algo.endswith('0') else \
        #         (R.kBlue   if algo.endswith('1') else \
        #          (R.kGreen if algo.endswith('2') else \
        #           (R.kRed  if algo.endswith('3') else R.kBlack)))

        grp[algo+'_peak'] .SetName('g_'+algo+'_peak')
        grp[algo+'_peak'] .SetTitle('Algo %s mass peak scale shift (%%)' % algo)
        grp[algo+'_width'].SetName('g_'+algo+'_width')
        grp[algo+'_width'].SetTitle('Algo %s mass peak width (%%)' % algo)
        grp[algo+'_peak'] .SetMarkerColor(color)
        grp[algo+'_peak'] .SetMarkerStyle(21)
        grp[algo+'_width'].SetMarkerColor(color)
        grp[algo+'_width'].SetMarkerStyle(21)
        grp[algo+'_peak'].Write()
        grp[algo+'_width'].Write()
        color += 1

    ## Plot overlays into pdf files
    for plot in ['peak','width']:
        can = R.TCanvas()
        can.cd()
        leg = R.TLegend()

        for algo in ALGOS:
            leg.AddEntry(grp[algo+'_'+plot], algo)
            if algo == ALGOS[0]:
                grp[algo+'_'+plot].Draw('ALP')
                if plot == 'peak':
                    grp[algo+'_'+plot].GetYaxis().SetRangeUser(-10,5)
                if plot == 'width':
                    grp[algo+'_'+plot].GetYaxis().SetRangeUser(0,50)
            else:
                grp[algo+'_'+plot].Draw('LPsame')
        leg.Draw('same')
        can.SaveAs(OUT_DIR+plot+'.pdf')
        del leg
        del can

    del grp
    outfile.Write()
    outfile.Close()
    infile.cd()
    del hist
    infile.Close()

## End function: def main()


if __name__ == '__main__':
    main()
