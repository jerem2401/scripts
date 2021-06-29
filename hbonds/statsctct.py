#!/usr/bin/env python
'''
Author: Jeremy Lapierre <jeremy.lapierre@uni-saarland.de>
Description:

'''

import argparse
import glob
import uncertainties as u
import numpy as np


def bywin(tdir):
    out = open(f"{tdir}/analyzehb.txt", 'a')
    out.write(f"[ Contact ]\n")

    xvgs = glob.glob(f"{tdir}/numcount.xvg")
    xvgs = sorted(xvgs, key=lambda x: x.split('/')[-2].split('_')[-1])

    for xvg in xvgs:
        with open(xvg, 'r') as reader:
            for line in reader:
                if "#@ s0" in line:
                    hmean = round(float(line.split('"')[1].split()[-1]), 2)
                elif "#@ s1" in line:
                    hsem = round(float(line.split('"')[1].split()[-1]), 2)
                elif "#std" in line:
                    if hsem > 10000:
                        std = round(float(line.split()[1]))
                        hsem = std/np.sqrt(2)
                else:
                    break

        hunc = u.ufloat(hmean, hsem)
        out.write(f"{hunc}\n")

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

    for i in glob.glob(f"{args.f}"):
        bywin(i)
        print(f"hbonds between DNA-prot done in {i} directory")


if __name__ == "__main__":
    main()
