#!/usr/bin/env python
'''
Author: Jeremy Lapierre <jeremy.lapierre@uni-saarland.de>
Description:

'''

import argparse
import numpy as np
from scipy import stats
import glob
import uncertainties as u


def bywin(tdir):
    toth = []

    out = open(f"{tdir}/analyzehb.txt", 'a')
    out.write(f"[ Prot-nuc H-bonds for nuc in bubble ]\n")
    out.write(f"#nuc_resid hbond\n")

    xvgs = glob.glob(f"{tdir}/DNA*/*.xvg")
    xvgs = sorted(xvgs, key=lambda x: x.split('/')[-2].split('_')[-1])

    for xvg in xvgs:
        hbond = []
        with open(xvg, 'r') as reader:
            for line in reader:
                if ("#" in line or "@" in line):
                    continue
                else:
                    hbond.append(int(line.split()[1]))
        hbond = np.array(hbond)
        hmean = round(np.mean(hbond), 2)
        hsem = round(stats.sem(hbond), 2)

        hunc = u.ufloat(hmean, hsem)
        toth.append(hunc)
        tmp = xvg.split('/')
        nucname = tmp[-2]
        out.write(f"{nucname:<35} {str(toth[-1]):>11}\n")

    usumt = u.ufloat(0, 0)
    for i in toth:
        usumt += i
    out.write(f"\n[ Prot-nuc H-bonds for the whole bubble ]\n")
    out.write(f"{usumt}\n")

    out.close()
#    d = dict(zip(name, toth))
#    tnuc = list(set(tnuc))
#    tnuc = [int(x.strip('t')) for x in tnuc]
#    tnuc.sort()
#
#    out.write(f"\n[ H-bonds/nuc in bubble ]\n")
#    out.write(f"#tnuc_nbr hbonds\n")
#
#    usumt = u.ufloat(0, 0)
#    for i in tnuc:
#        usum = u.ufloat(0, 0)
#        for key, val in d.items():
#            pat = 't'+str(i)
#            if pat in key:
#                usum += val
#        out.write(f"{pat} {usum}\n")
#        usumt += usum
#
#    out.write(f"\n[ H-bonds 4 the whole bubble ]\n")
#    out.write(f"{usumt}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='give target directories containing hbond.xvg files',
                        action='store', dest='f', required=True)
    args = parser.parse_args()

    for i in glob.glob(f"{args.f}/hbondP"):
        bywin(i)
        print(f"hbonds between DNA-prot done in {i} directory")


if __name__ == "__main__":
    main()
