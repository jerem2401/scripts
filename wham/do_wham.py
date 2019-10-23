#!/usr/bin/env python
#This version have not been tested with bCycl=True => careful or use do_whampy2.py

import numpy as np
import os, sys, re, math


def do_wham(z, histos, pos, kappa, RT = 2.49434, tol = 1e-6, bCycl = False,
                zmin = float('nan'), zmax = float('nan')):

    maxchangeZ = 1
    nHist      = len(histos)
    nz_input   = len(z)
    dz         = z[1] - z[0]


    if (bCycl):
        if (math.isnan(zmin) or math.isnan(zmax)):
            print('Error, with bcycl, you must provide zmin and zmax')
    else:
        if (not math.isnan(zmin) or not math.isnan(zmax)):
            print('Error, provide zmin and zmax only with bCycl')
        zmin = z[0]
        zmax = z[-1]
            
    izmin = int(round((zmin-z[0])/dz))
    izmax = int(round((zmax-z[0])/dz))
    if (izmin < 0):
        izmin = 0
    if (izmax > nz_input-1):
        izmax = nz_input-1

    nzused = izmax-izmin+1

    zused = z[izmin:izmin+nzused]
    
    print('\nDoing WHAM:')
    print('\t# of histograms = %d' % nHist)
    print('\t# of bins given = %d' % nz_input)
    print('\t# of bins used  = %d' % nzused)
    print('\tz range given   = %g to %g' % (z    [0], z    [-1]+dz))
    print('\tz range used    = %g to %g' % (zused[0], zused[-1]+dz))
    print('\tTolerance       = %g' % tol)

    histos_used = np.zeros([nHist, nzused])
    
    if (bCycl == False):
        histos_used = histos
    else:
        for i in range(nHist):
            nBinsOutside_lower = izmin
            nBinsOutside_upper = nz_input - izmax
            # Map histos outside boundaryies [izmin:izmax] back into the [izmin:izmax]
            histos_used[i]                                    = histos[i][izmin:izmin+nzused]
            if nBinsOutside_upper > 0:
                histos_used[i][nzused-nBinsOutside_lower:] += histos[i][0:nBinsOutside_lower]
            if nBinsOutside_upper:
                histos_used[i][0:nBinsOutside_upper]         += histos[i][nz_input-nBinsOutside_upper:]
    
    expZ         = np.ones(nHist)   # = e^Z[i]
    nPointsHist  = np.zeros(nHist)
    nPoints_times_exp_minus_U_over_RT = []
    denom        = np.zeros(nzused)
    numer        = np.zeros(nzused)
    exp_minus_U_over_RT = []              # exp( -V[umb]/RT)

    # Write all histos into one xvg file
    if (bCycl):
        fp = open('histo_1D_cycl.xvg', 'w')
        for i in range(len(zused)):
            #print >> fp, '%10g ' % (zused[i]),
            print("{:%10g}".format(zused[i]), file=fp)
            for w in range(nHist):
                #print >> fp, ' %10g' % (histos_used[w][i]),
                print("{:%10g}".format(histos_used[w][i]), file=fp)
            print >> fp #what is this ? How to translate to python 3.X ?
        fp.close()

    
    # Store e^(-U/RT) and Npoints * e^(-U/RT) for all histos, so
    # to avoid that we need to compute a lot of exp() in the loop
    for i in range(nHist):
        nPointsHist[i] = np.sum(histos_used[i])

        if (bCycl == False):
            sqr_distance = ((zused+dz/2) - pos[i])**2
        else:
            # Get shortest distance, possibly over the periodic boundaries
            L             = zused[-1] - zused[0] + dz
            deltaz        = (zused+dz/2) - pos[i]
            sqr_distance1 = (deltaz   )**2
            sqr_distance2 = (deltaz + L)**2
            sqr_distance3 = (deltaz - L)**2
            sqr_distance  = np.min([sqr_distance1, sqr_distance2, sqr_distance3], axis=0)            
        
        expU            = np.exp(-0.5 * kappa[i] * sqr_distance / RT)
        exp_minus_U_over_RT.append( expU )
        nPoints_times_exp_minus_U_over_RT.append( nPointsHist[i] * expU)
    
    # The main WHAM loop
    N = 1
    while (maxchangeZ > tol):

        bCheckChange = ((N % 100) == 0)
        if (bCheckChange):
            Zold = np.log(expZ)

        # Compute profile
        # E.g.: Roux, Comp Phys Comm, 91, 275-282 (1995), eq. 8
        numer.fill(0)
        denom.fill(0)
        for i in range(nHist):
            numer += histos_used[i]
            denom += nPoints_times_exp_minus_U_over_RT[i] * expZ[i]

        rho = numer / denom

        # Update Z
        for i in range(nHist):
            expZ[i] = 1.0 / (np.sum(exp_minus_U_over_RT[i] * rho))
            
        # Get change in Z
        if (bCheckChange):
            Znew       = np.log(expZ)
            maxchangeZ = (np.fabs(np.log(expZ) - Zold)).max()
            print('WHAM iteration %5d, maxchange of Z\'s: %g' % (N, maxchangeZ))
            
        N += 1
    
    return zused, rho, np.log(expZ)
