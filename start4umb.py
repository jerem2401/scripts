#!/usr/bin/env python

import sys
import numpy as np
import pandas as pd
import plumed_pandas

def main():
    '''
    parsing the colvar.txt, find the closest nCV value to target_val
    '''

    #read COLVAR.txt
    df=plumed_pandas.read_as_pandas(sys.argv[1])

    #extract time of frame closest to the target_value and cnt of gksi at this frame
    tvalue=df.iloc[(df['nCV']-float(sys.argv[2])).abs().idxmin()]['time']
    gksi=df.iloc[(df['nCV']-float(sys.argv[2])).abs().idxmin()]['restraint2.RMSDMID_cntr']
    return (tvalue,gksi)

# If called from the command line... 
if __name__ == "__main__":
    print(main())
