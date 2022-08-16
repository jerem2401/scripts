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


def pull(tdir):
    from itertools import dropwhile
    import re

    out = open(f"{tdir}/hbonds_pull.txt", 'a')
    files = glob.glob(f"{tdir}/_*")
#    #reduce bubble size
    r = re.compile('.*_t(5[2-9]|6[0-3])')
    files = list(filter(r.match, files))

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

    for i in glob.glob(f"{args.f}/hbond"):
        if args.pull:
            pull(i)
            print("doing pull")
            print(f"{i} direcotry is done")
        else:
            bywin(i)
            print(f"{i} direcotry is done")


if __name__ == "__main__":
    main()
