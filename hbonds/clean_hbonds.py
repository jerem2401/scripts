#!/usr/bin/env python
'''
Author: Jeremy Lapierre <jeremy.lapierre@uni-saarland.de>
Description:

'''

import argparse
import glob
import uncertainties as u
import numpy as np


Vdef bywin(tdir):
    toth = []

    out = open(f"{tdir}/analyzehb.txt", 'a')
    out.write(f"[ Prot-nuc H-bonds for nuc in bubble ]\n")
    out.write(f"#nuc_resid hbond\n")

    xvgs = glob.glob(f"{tdir}/DNA*/hbond.xvg")
    xvgs = sorted(xvgs, key=lambda x: x.split('/')[-2].split('_')[-1])

    for xvg in xvgs:
        with open(xvg, 'r') as reader:
            for line in reader:
                if "#@ s0" in line:
                    hmean = round(float(line.split('"')[1].split()[-1]), 2)
                elif "#@ s1" in line:
                    hsem = round(float(line.split('"')[1].split()[-1]), 2)
                elif "#std" in line:
                    if hsem > 10:
                        std = round(float(line.split()[1]))
                        hsem = std/np.sqrt(2)
                else:
                    break

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

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='give target directories containing hbond.xvg files',
                        action='store', dest='f', required=True)
    args = parser.parse_args()

    with open(f"{args.f}/hbond/analyzehb.txt", 'r'):
        for line in reader:
            if "E_49.360_t49_nt45_O6_N4_H41_H42" in line:

        print(f"hbonds between DNA-prot done in {i} directory")


if __name__ == "__main__":
    main()
