#!/usr/bin/env python

import numpy as np
import os, sys, re, math
from do_wham import do_wham
import matplotlib.pyplot as plt
import matplotlib
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f', help='give histo files for wham', action='store', nargs='*', dest='f')
args = parser.parse_args()


RT = 2.49434

#histo_files = os.popen('ls histo* | sort -t _ -k 2 -n').read().split()
histo_files = args.f
N = len(histo_files)

hist  = []
pos   = np.zeros(N)
kappa = np.zeros(N)


n = 0
for h in histo_files:

    print('Reading %s' % h)
    fp = open(h, 'r')

    l      = fp.readline()
    l      = fp.readline()
    pos1   = float(l.split()[3])
    l      = fp.readline()    
    kappa1 = float(l.split()[3])
    fp.close()

    z, thishist = np.loadtxt(h, unpack=True)
    hist.append( thishist )
    
    print ('\tpos = %10g  kappa = %10g   %d bins' % (pos1, kappa1, len(z)))

    pos  [n] = pos1
    kappa[n] = kappa1
    n       += 1

    

zused, rho, Z = do_wham(z, hist, pos, kappa, RT = 2.49434, tol = 1e-6, bCycl = False)

#The main WHAM loop
#    N = 1
#    while (maxchangeZ > tol):
#
#        bCheckChange = ((N % 100) == 0)
#        if (bCheckChange):
#            Zold = np.log(expZ)
#
#        # Compute profile
#        # E.g.: Roux, Comp Phys Comm, 91, 275-282 (1995), eq. 8
#        numer.fill(0)
#        denom.fill(0)
#        for i in range(nHist):
#            numer += histos_used[i]
#            denom += nPoints_times_exp_minus_U_over_RT[i] * expZ[i]
#
#        rho = numer / denom
#
#        # Update Z
#        for i in range(nHist):
#            expZ[i] = 1.0 / (np.sum(exp_minus_U_over_RT[i] * rho))
#            
#        # Get change in Z
#        if (bCheckChange):
#            Znew       = np.log(expZ)
#            maxchangeZ = (np.fabs(np.log(expZ) - Zold)).max()
#            print('WHAM iteration %5d, maxchange of Z\'s: %g' % (N, maxchangeZ))
#            
#        N += 1
#    
#    return zused, rho, np.log(expZ)


pmf = np.zeros_like(rho)
for i in range(len(rho)):
    if (rho[i] > 0):
        pmf[i] = -RT * np.log(rho[i])
    else:
        pmf[i] = 0

fp = open('test_wham_pmf.out', 'w')
for i in range(len(zused)):
    print(str(zused[i])+'    '+str(pmf[i]), file=fp)
    #print("{:.3f}    {:.3f}".format(zused[i],pmf[i]), file=fp) #To do: adding format to remove digits
fp.close()

fp = open('test_wham_z.out', 'w')
for i in range(N):
    print(str(pos[i])+'    '+str(Z[i]), file=fp)
fp.close()

#plot PMF

matplotlib.rcParams.update({'font.size': 25})

x, y = np.loadtxt('test_wham_pmf.out', unpack=True, delimiter='    ')

fig = plt.figure(figsize=(12, 8))
ax1 = fig.add_subplot(1, 1, 1)

ax1.plot(x, y, '-')
ax1.set_xlabel(r'$\xi_{new}$')
ax1.set_ylabel('Free energy (kJ/mol)')

id=0
while os.path.exists("pmf.%s.jpeg" % id):
    id+=1
plt.savefig("pmf."+str(id)+".jpeg")
plt.close()
