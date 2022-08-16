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


def pull(tdir):
    from itertools import dropwhile
    import re

    print('doing pull MF')
    xvgs = glob.glob(f"{tdir}/DNA*/hbond.xvg")
    xvgs = sorted(xvgs, key=lambda x: x.split('/')[-2].split('_')[-1])
    out = open(f"{tdir}/hbondsP_pull.txt", 'a')
#    #reduce bubble size
    r = re.compile('.*(DNAt_5[2-9]|DNAt_6[0-3]|DNAnt_3[1-9]|DNAnt_4[0-2])')
    files = list(filter(r.match, xvgs))

    tot = []
    time = []
    for file in files:
        with open(file, 'r') as reader:
            hb = []
            for line in dropwhile(lambda x: '#' in x or '@' in x, reader):
                hb.append(int(line.split()[1]))
                if file == files[0]:
                    time.append(float(line.split()[0]))
            tot.append(hb)

    tot = np.array(tot)
    tot_sum = np.sum(tot, axis=0)

    for i in range(len(time)):
        out.write(f"{time[i]:<10}{tot_sum[i]:>5}\n")
    out.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='give target directories containing hbond.xvg files',
                        action='store', dest='f', required=True)
    parser.add_argument('-pull', help='If pulling sim instead of US',
                        action='store_true', dest='pull')
    args = parser.parse_args()

    for i in glob.glob(f"{args.f}/hbondP"):
        if args.pull:
            pull(i)
            print("doing pull")
            print(f"hbonds between DNA-prot done in {i} directory")
        else:
            bywin(i)
            print(f"hbonds between DNA-prot done in {i} directory")


if __name__ == "__main__":
    main()
