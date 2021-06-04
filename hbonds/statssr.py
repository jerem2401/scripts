#!/usr/bin/env python
'''
Author: Jeremy Lapierre <jeremy.lapierre@uni-saarland.de>
Description:

'''

import argparse
import glob
from uncertainties import ufloat


def bywin(tdir):
    out = open(f"{tdir}/analyzehb.txt", 'a')
    out.write(f"[ Prot-DNA nb SR interactions ]\n")
    out.write(f"#coul lj\n")

    xvgs = glob.glob(f"{tdir}/outnrg.txt")
    xvgs = sorted(xvgs, key=lambda x: x.split('/')[-2].split('_')[-1])

    for xvg in xvgs:
        with open(xvg, 'r') as reader:
            for line in reader:
                if "Coul-SR:Protein-DNA" in line and 'LJ-SR:Protein-DNA' not in line:
                    mean_coul = float(line.split()[1])
                    err_coul = float(line.split()[2])
                elif "LJ-SR:Protein-DN" in line and 'Coul-SR:Protein-DNA' not in line:
                    mean_lj = float(line.split()[1])
                    err_lj = float(line.split()[2])

        coul = ufloat(mean_coul, err_coul)
        lj = ufloat(mean_lj, err_lj)
        out.write(f"{coul}     {lj}\n")
        out.write(f"\n[ Prot-DNA Coul+LJ ]\n")
        sr_sum = coul+lj
        out.write(f"{sr_sum}\n")

    out.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='give target directories containing hbond.xvg files',
                        action='store', dest='f', required=True)
    args = parser.parse_args()

    for i in glob.glob(f"{args.f}"):
        bywin(i)
        print(f"SR between DNA-prot done in {i} directory")


if __name__ == "__main__":
    main()
