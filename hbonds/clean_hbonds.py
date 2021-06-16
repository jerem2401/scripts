#!/usr/bin/env python
'''
Author: Jeremy Lapierre <jeremy.lapierre@uni-saarland.de>
Description:

'''

import argparse
import glob
from uncertainties import ufloat


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='give target directories containing hbond.xvg files',
                        action='store', dest='f', required=True)
    args = parser.parse_args()

    files = glob.glob(f"{args.f}")
    for i in files:
        vals = []
        with open(i, 'r') as reader:
            for line in reader:
                if '[' not in line and '#' not in line and line != "\n":
                    vals.append(line.split()[1])

        tot = ufloat(0, 0)
        for j in vals:
            mean = float(j.split('+/-')[0])
            err = float(j.split('+/-')[1])
            tot += ufloat(mean, err)

        with open(i, 'a') as writer:
            writer.write(f"[ Total ]\n{tot}")

        print(f"{i} is cleaned")


if __name__ == "__main__":
    main()
