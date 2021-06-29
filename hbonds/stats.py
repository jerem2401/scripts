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
    toth = []
    name = []
    tnuc = []

    out = open(f"{tdir}/analyzehb.txt", 'a')
    out.write(f"[ induvidual H bond for nuc in bubble ]\n")
    out.write(f"#nuc_and_atm_name hbond\n")

    for xvg in glob.glob(f"{tdir}/E_*.xvg"):
        with open(xvg, 'r') as reader:
            hmean = 0
            hsem = 0
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
