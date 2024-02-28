#! /usr/bin/env python

## Script to compute signal vs. background for different ParticleNet "a" boson mass regressions
import os
import sys
import math
import ROOT as R
from array import array

R.gROOT.SetBatch(True)  ## Don't display histograms or canvases when drawn
R.gStyle.SetOptStat(0)  ## Don't display stat boxes

## User configuration
VERBOSE  = False
REBIN    = True
MASSES   = [12, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]
#ALGOS    = ['v0','v1','v2','v3']  ## vX
#ALGOS    = ['avg_v01','avg_v02','avg_v03','avg_v12','avg_v13','avg_v23']  ## vXY
#ALGOS    = ['avg_v012','avg_v013','avg_v023','avg_v123','avg_v0123']  ## vXYZ
#ALGOS     = ['v0','avg_v03','avg_v13','avg_v0123']  ## best8A
#ALGOS     = ['v3','avg_v01','avg_v013','avg_v023']  ## best8B
#ALGOS     = ['v0','avg_v01','avg_v03','avg_v013']  ## best4
ALGOS     = ['avg_v01','avg_v03','avg_v013']  ## best3
#MASSES   = [55]
#ALGOS    = ['v3']
WP_CUT   = 'WP40'  ## WP40, WP60, or WP80
BKGS     = ['QCD_0bCat','QCD_1bCat','QCD_2bCat','QCD_3bCat','QCD_4bCat','QCD_5bAndMoreCat','TTToHadronic_powheg','SingleTop','ZJetsToQQ_HT','WJetsToQQ_HT']
WINDOWS  = [1.0,2.0,4.0]  ## Size of averaging window in GeV

IN_DIR   = '/eos/cms/store/user/ssawant/htoaa/analysis/'
IN_FILE  = IN_DIR+'20231122_MC_woTrg/2018/analyze_htoaa_stage1.root' ## Include backgrounds
OUT_FILE = 'plots/HtoAA_massA_sig_bkg_%s.root' % WP_CUT
OUT_DIR  = 'plots/pdf/HtoAA_massA_sig_bkg_%s_best3_FWHM/' % WP_CUT


def find_significance(hstS, hstB, wind):
    sumSSB = 0
    for i in range(1, hstS.GetNbinsX()+1):
        nSig = hstS.GetBinContent(i)
        # if nSig < 0.1*hstS.GetMaximum(): continue
        if nSig < 0.5*hstS.GetMaximum(): continue
        bin_mass = hstS.GetBinCenter(i)
        ## Shift mass window if it includes masses < 11.5
        if bin_mass - (wind/2.0) < 11.5:
            if VERBOSE:
                print('  * Shifting mass center from %.2f to %.2f' % (bin_mass, min(11.5, bin_mass)+(wind/2.0)))
            bin_mass = min(11.5, bin_mass) + (wind/2.0)
        ## Find average background value within mass window
        total = 0
        nBins = 0
        for j in range(1,hstB.GetNbinsX()+1):
            #print('  - Bin %d, center %.3f vs. mass %.3f' % (j,hstB.GetBinCenter(j),bin_mass))
            if abs(bin_mass - hstB.GetBinCenter(j)) <= (wind/2.0):
                #print('    ** KEEP!')
                total += hstB.GetBinContent(j)
                nBins += 1
        ## Use wider window if there are no background events
        fact = 1
        #if total == 0: sys.exit()
        while total == 0:
            print('\nBin %d, mass %.2f, window %.2f (%d) has %.2f signal but no background! Widening.' % (i, bin_mass, wind*fact, nBins, nSig))
            fact *= 2
            nBins = 0
            for jj in range(1,hstB.GetNbinsX()+1):
                if abs(bin_mass - hstB.GetBinCenter(jj)) <= (fact*wind/2.0):
                    total += hstB.GetBinContent(jj)
                    nBins += 1
        ## Compute average background per signal bin width
        nBkg = (total/nBins)*(hstS.GetBinWidth(1)/hstB.GetBinWidth(1))
        #try:
        sumSSB += (pow(nSig,2)/nBkg)
        #except:
        if VERBOSE:
            print('Bin %d, mass %.2f, window %.2f (%d), nSig = %.2f, nBkg = %.2f' % (i, bin_mass, wind, nBins, nSig, nBkg))

    return math.sqrt(sumSSB)
## end function: find_significance(hstS, hstB, wind)


def main():

    print('\nInside HtoAA_massA_sig_bkg\n')

    print('\nGetting input file %s' % IN_FILE)
    infile = R.TFile(IN_FILE,'open')
    outfile = R.TFile(OUT_FILE,'recreate')
    if not os.path.exists(OUT_DIR):
        os.makedirs(OUT_DIR)
    infile.cd()

    sig_beg  = 'evt/SUSY_GluGluH_01J_HToAATo4B_M-'
    sig_end  = '_HPtAbv150/'
    hist_beg = 'hLeadingFatJetParticleNet_massA_Hto4b_'
    hist_end = '_SR%s_central' % WP_CUT
    hist_bkg = {}
    hist_sig = {}
    colors   = {}

    ## Store summed S/sqrt(B) values
    masses_val = {}
    masses_err = {}
    srbs_val   = {}
    srbs_hi    = {}
    srbs_lo    = {}

    color = 2
    for algo in ALGOS:
        masses_val[algo] = array('d')
        masses_err[algo] = array('d')
        srbs_val  [algo] = array('d')
        srbs_hi   [algo] = array('d')
        srbs_lo   [algo] = array('d')
        hist_bkg  [algo] = 0
        hist_sig  [algo] = {}
        colors    [algo] = (R.kViolet if color == 5 else color)
        color += 1
        for mass in MASSES:
            hist_sig[algo][str(mass)] = 0

    ## Loop over algos
    for algo in ALGOS:
        if VERBOSE:
            print('\nStarting to look at algo %s' % algo)
            print('\nGetting background histogram '+'evt/'+BKGS[0]+'/'+hist_beg+algo+hist_end)
        hst_bkg = infile.Get('evt/'+BKGS[0]+'/'+hist_beg+algo+hist_end)
        if VERBOSE:
            print('\n%s %s has integral %.1f' % (BKGS[0], algo, hst_bkg.Integral()))
        for bkg in BKGS:
            if bkg == BKGS[0]: continue
            hst_bkg.Add(infile.Get('evt/'+bkg+'/'+hist_beg+algo+hist_end))
            if VERBOSE:
                print('  * With %s, integral %.1f' % (bkg, hst_bkg.Integral()))

        ## Loop over mass points
        for mass in MASSES:
            if VERBOSE:
                print('\nStarting to look at M-%d' % mass)
                print('\nGetting signal histogram '+sig_beg+str(mass)+sig_end+hist_beg+algo+hist_end)
            hst_sig = infile.Get(sig_beg+str(mass)+sig_end+hist_beg+algo+hist_end)
            if VERBOSE:
                print('\nM-%d %s has integral %.1f, peak %.2f' % (mass, algo, hst_sig.Integral(), hst_sig.GetMaximum()))
            ## Rebin for larger "a" boson masses, with wider distributions
            rebin = 1
            if REBIN:
                binW  = hst_sig.GetBinWidth(1)
                if mass == 12:
                    rebin = int(math.floor(0.05/binW))
                else:
                    rebin = int(math.floor((0.05/binW)*(mass/3.01)))
                rebin = max(rebin, 1)
            nbins = hst_sig.GetNbinsX()
            drop  = (nbins % rebin)
            if VERBOSE or True:
                print('Rebinning M-%d by a factor of %d' % (mass, rebin))

            ## Drop last bins to get an even division for the rebinning
            print('Cloning M-%d histogram and dropping %d bins' % (mass, drop))
            hst_sig_new = R.TH1D(hst_sig.GetName()+'_sig_copy', hst_sig.GetTitle()+'_sig_copy', nbins-drop,
                                 hst_sig.GetBinLowEdge(1), hst_sig.GetBinLowEdge(nbins-drop+1))
            for ii in range(1,nbins-drop+1):
                hst_sig_new.SetBinContent(ii, hst_sig.GetBinContent(ii))
                hst_sig_new.SetBinError  (ii, hst_sig.GetBinError(ii))
            ## Sanity check that the histograms are identical
            for ii in range(1, hst_sig_new.GetNbinsX()):
                if abs(hst_sig_new.GetBinContent(ii) - hst_sig.GetBinContent(ii)) > 0.0001 or \
                   abs(hst_sig_new.GetBinError(ii)   - hst_sig.GetBinError(ii))   > 0.0001 or \
                   abs(hst_sig_new.GetBinLowEdge(ii) - hst_sig.GetBinLowEdge(ii)) > 0.001:
                    print('Mismatch in signal bin %d!!! Quitting.' % ii)
                    print(hst_sig_new.GetBinContent(ii), hst_sig.GetBinContent(ii))
                    print(hst_sig_new.GetBinError(ii), hst_sig.GetBinError(ii))
                    print(hst_sig_new.GetBinLowEdge(ii), hst_sig.GetBinLowEdge(ii))
                    sys.exit()

            ## Rebin and find sqrt of sum of S^2/B
            hst_sig_new.Rebin(rebin)
            srbs = []
            widths = []
            for window in WINDOWS:
                srb = find_significance(hst_sig_new, hst_bkg, window)
                print('For M-%d %s window %.2f, S/sqrt(B) = %.3f' % (mass, algo, window, srb))
                srbs.append(srb)

            ## Save values for TGraphAsymmErrors
            masses_val[algo].append(mass)
            masses_err[algo].append(0)
            srbs_val [algo].append(sum(srbs)/len(srbs))
            srbs_hi  [algo].append(max(srbs) - (sum(srbs)/len(srbs)))
            srbs_lo  [algo].append((sum(srbs)/len(srbs)) - min(srbs))

            ## Save signal (rebinned) mass distribution
            outfile.cd()
            hist_sig[algo][str(mass)] = hst_sig_new.Clone('ggH_mA-%d_PNet_%s' % (mass, algo))
            hist_sig[algo][str(mass)].SetLineWidth(2)
            hist_sig[algo][str(mass)].SetLineColor(colors[algo])
            infile.cd()

            del hst_sig
            del hst_sig_new
        ## End loop: for mass in MASSES

        ## Save background mass distribution
        outfile.cd()
        hist_bkg[algo] = hst_bkg.Clone('bkg_PNet_%s' % algo)
        rebin_bkg = max(int(math.floor(0.5/hist_bkg[algo].GetBinWidth(1))), 1)
        hist_bkg[algo].Rebin(rebin_bkg)
        hist_bkg[algo].SetLineWidth(2)
        hist_bkg[algo].SetLineColor(colors[algo])
        infile.cd()
        del hst_bkg
    ## End loop: for algo in ALGOS


    ## Generate TGraphAsymmErrors
    outfile.cd()
    grp = {}
    for algo in ALGOS:
        grp[algo+'_srb'] = R.TGraphAsymmErrors(len(MASSES), masses_val[algo], srbs_val[algo],
                                                masses_err[algo], masses_err[algo],
                                                srbs_lo[algo], srbs_hi[algo])

        # color = R.kViolet  if algo.endswith('0') else \
        #         (R.kBlue   if algo.endswith('1') else \
        #          (R.kGreen if algo.endswith('2') else \
        #           (R.kRed  if algo.endswith('3') else R.kBlack)))

        grp[algo+'_srb'] .SetName('g_'+algo+'_srb')
        grp[algo+'_srb'] .SetTitle('Algo %s S/sqrt(B)' % algo)
        grp[algo+'_srb'] .SetMarkerColor(colors[algo])  ## Use violet instead of yellow
        grp[algo+'_srb'] .SetMarkerStyle(21)
        grp[algo+'_srb'].Write()

    ## Plot overlays into pdf files
    yScl = 2.0 if WP_CUT == 'WP40' else (1.5 if WP_CUT == 'WP60' else 1.0) 
    for plot in ['srb']:
        can = R.TCanvas()
        can_bkg = R.TCanvas()
        can_sig = {}
        leg = R.TLegend(0.58,0.68,0.88,0.88)
        leg_bkg = R.TLegend(0.12,0.68,0.42,0.88)
        leg_sig = {}
        for mass in MASSES:
            can_sig[str(mass)] = R.TCanvas()
            if mass > 37 or mass == 20:
                leg_sig[str(mass)] = R.TLegend(0.12,0.68,0.42,0.88)
            else:
                leg_sig[str(mass)] = R.TLegend(0.58,0.68,0.88,0.88)

        for algo in ALGOS:
            ## Significance trends
            can.cd()
            leg.AddEntry(grp[algo+'_'+plot], algo)
            if algo == ALGOS[0]:
                grp[algo+'_'+plot].Draw('ALP')
            else:
                grp[algo+'_'+plot].Draw('LPsame')
            ## Background mass distributions
            can_bkg.cd()
            leg_bkg.AddEntry(hist_bkg[algo], algo)
            if algo == ALGOS[0]:
                hist_bkg[algo].Draw('hist')
            else:
                hist_bkg[algo].Draw('histsame')
            ## Signal mass distributions
            for mass in MASSES:
                can_sig[str(mass)].cd()
                leg_sig[str(mass)].AddEntry(hist_sig[algo][str(mass)], algo)
                if algo == ALGOS[0]:
                    hist_sig[algo][str(mass)].Draw('hist')
                else:
                    hist_sig[algo][str(mass)].Draw('histsame')
            ## End loop: for mass in MASSES
        ## End loop: for algo in ALGOS

        ## Overlay significance graphs
        can.cd()
        leg.Draw('same')
        grp[ALGOS[0]+'_'+plot].GetYaxis().SetRangeUser(0,120*yScl)
        can.SaveAs(OUT_DIR+plot+'_lin.pdf')
        grp[ALGOS[0]+'_'+plot].GetYaxis().SetRangeUser(0,40*yScl)
        can.SaveAs(OUT_DIR+plot+'_zoom.pdf')
        grp[ALGOS[0]+'_'+plot].GetYaxis().SetRangeUser(0,20*yScl)
        can.SaveAs(OUT_DIR+plot+'_zoom2.pdf')
        grp[ALGOS[0]+'_'+plot].GetYaxis().SetRangeUser(1*yScl,200*yScl)
        can.SetLogy()
        can.SaveAs(OUT_DIR+plot+'_log.pdf')

        ## Overlay background mass distributions
        can_bkg.cd()
        leg_bkg.Draw('same')
        hist_bkg[ALGOS[0]].GetXaxis().SetRangeUser(10,30)
        can_bkg.SaveAs(OUT_DIR+'bkg_zoom.pdf')
        hist_bkg[ALGOS[0]].GetXaxis().SetRangeUser(0,70)
        for algo in ALGOS:
            hist_bkg[algo].Rebin(5)
        hist_bkg[ALGOS[0]].GetYaxis().SetRangeUser(0, hist_bkg[ALGOS[0]].GetMaximum()*1.3)
        can_bkg.SaveAs(OUT_DIR+'bkg_lin.pdf')

        ## Overlay signal mass distributions
        for mass in MASSES:
            can_sig[str(mass)].cd()
            leg_sig[str(mass)].Draw('same')
            hist_sig[ALGOS[0]][str(mass)].GetXaxis().SetRangeUser(mass*(2.0/3), mass*(4.0/3))
            ## Custom mass window for low masses, include possible "spike" at 12 GeV
            if mass < 21:
                hist_sig[ALGOS[0]][str(mass)].GetXaxis().SetRangeUser(11,16+5*(mass>14)+5*(mass>19))
            hist_sig[ALGOS[0]][str(mass)].GetYaxis().SetRangeUser(0, hist_sig[ALGOS[0]][str(mass)].GetMaximum()*1.5)
            can_sig[str(mass)].SaveAs(OUT_DIR+'ggH_mA-%d.pdf' % mass)
            if mass < 21:
                hist_sig[ALGOS[0]][str(mass)].GetYaxis().SetRangeUser(0, hist_sig[ALGOS[0]][str(mass)].GetMaximum()*0.2)
                can_sig[str(mass)].SaveAs(OUT_DIR+'ggH_mA-%d_zoom.pdf' % mass)

        del leg_sig
        del can_sig
        del leg_bkg
        del can_bkg
        del leg
        del can
    ## End loop: for plot in ['srb']

    del grp
    outfile.Write()
    outfile.Close()
    infile.cd()
    del hist_bkg
    infile.Close()

## End function: def main()


if __name__ == '__main__':
    main()
