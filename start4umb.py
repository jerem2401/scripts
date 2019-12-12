#!/usr/bin/env python

import argparse
import numpy as np
import pandas as pd
import plumed_pandas

def main():
    '''
    parsing the colvar.txt, find the closest nCV value to target_val
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='give colvar_files', action='store', dest='f', type=str)
    parser.add_argument('-v', help='tnCV', action='store', dest='v', type=float)
    args = parser.parse_args()

    #read COLVAR.txt
    df=plumed_pandas.read_as_pandas(args.f)

    #extract time of frame closest to the target_value and cnt of gksi at this frame 
    if args.v == 1.00 or args.v == -1.00:
        tvalue=df.iloc[(df['nCV']-float(args.v)).abs().idxmin()]['time']
        gksi=df.iloc[(df['nCV']-float(args.v)).abs().idxmin()]['restraint2.RMSDMID_cntr']
    else:
        tvalue=df.iloc[(df['restraint.nCV_cntr']-float(args.v)).abs().idxmin()]['time']
        gksi=df.iloc[(df['restraint.nCV_cntr']-float(args.v)).abs().idxmin()]['restraint2.RMSDMID_cntr']

    return (tvalue,gksi)

# If called from the command line... 
if __name__ == "__main__":
    print(main())
