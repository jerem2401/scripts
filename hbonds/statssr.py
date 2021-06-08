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

    xvgs = glob.glob(f"{tdir}/outnrg*.txt")
    xvgs = sorted(xvgs, key=lambda x: x.split('/')[-2].split('_')[-1])

    for xvg in xvgs:
        with open(xvg, 'r') as reader:
            for line in reader:
                if "Coul-SR:Protein-DNA" in line and 'LJ-SR:Protein-DNA' not in line:
                    tag = 'Protein-DNA'
                    mean_coul = float(line.split()[1])
                    err_coul = float(line.split()[2])
                elif "LJ-SR:Protein-DNA" in line and 'Coul-SR:Protein-DNA' not in line:
                    mean_lj = float(line.split()[1])
                    err_lj = float(line.split()[2])
                elif "Coul-SR:DNA-DNA" in line and 'LJ-SR:DNA-DNA' not in line:
                    tag = 'DNA-DNA'
                    mean_coul = float(line.split()[1])
                    err_coul = float(line.split()[2])
                elif "LJ-SR:DNA-DNA" in line and 'Coul-SR:DNA-DNA' not in line:
                    mean_lj = float(line.split()[1])
                    err_lj = float(line.split()[2])
                if "Coul-SR:DNA-Water" in line and 'LJ-SR:DNA-Water' not in line:
                    tag = 'DNA-Water'
                    mean_coul = float(line.split()[1])
                    err_coul = float(line.split()[2])
                elif "LJ-SR:DNA-Water" in line and 'Coul-SR:DNA-Water' not in line:
                    mean_lj = float(line.split()[1])
                    err_lj = float(line.split()[2])

        out.write(f"[ {tag} nb SR interactions ]\n")
        out.write(f"#coul lj\n")
        coul = ufloat(mean_coul, err_coul)
        lj = ufloat(mean_lj, err_lj)
        out.write(f"{coul}     {lj}\n")
        out.write(f"\n[ {tag} Coul+LJ ]\n")
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
