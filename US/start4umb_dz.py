#!/usr/bin/env python

import plumed_pandas

def main(colvf, ksi):
    '''
    parsing the colvar.txt, find the closest nCV value to target_val
    '''

    #read COLVAR.txt
    df = plumed_pandas.read_as_pandas(colvf)

    mindiff = (df['d.z']-float(ksi)).abs().idxmin()
    found = df.iloc[mindiff]['d.z']
    tvalue = df.iloc[mindiff]['time']

    return found, tvalue

# If called from the command line...
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='give colvar_files', action='store', dest='f', type=str)
    parser.add_argument('-v', help='tnCV', action='store', dest='v', type=float)
    args = parser.parse_args()

    F, TV = main(args.f, args.v)
    print(F, TV)
