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
    name = []
    tnuc = []

    out = open(f"{tdir}/analyzehb.txt", 'a')
    out.write(f"[ induvidual H bond for nuc in bubble ]\n")
    out.write(f"#nuc_and_atm_name hbond\n")

    for xvg in glob.glob(f"{tdir}/*.xvg"):
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
        tmp = xvg.rstrip('.xvg')
        name.append(tmp.split("/")[-1])

        tnuc.append(name[-1].split("_")[2])

        out.write(f"{name[-1]:<35} {str(toth[-1]):>11}\n")

    d = dict(zip(name, toth))
    tnuc = list(set(tnuc))
    tnuc = [int(x.strip('t')) for x in tnuc]
    tnuc.sort()

    out.write(f"\n[ H-bonds/nuc in bubble ]\n")
    out.write(f"#tnuc_nbr hbonds\n")

    usumt = u.ufloat(0, 0)
    for i in tnuc:
        usum = u.ufloat(0, 0)
        for key, val in d.items():
            pat = 't'+str(i)
            if pat in key:
                usum += val
        out.write(f"{pat} {usum}\n")
        usumt += usum

    out.write(f"\n[ H-bonds 4 the whole bubble ]\n")
    out.write(f"{usumt}\n")

    out.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='give target directories containing hbond.xvg files',
                        action='store', dest='f', required=True)
    args = parser.parse_args()

    for i in glob.glob(f"{args.f}/hbond"):
        bywin(i)
        print(f"{i} direcotry is done")


if __name__ == "__main__":
    main()
