#!/usr/bin/env python

import argparse
import plumed_pandas

def main():
    '''
    parsing the colvar.txt, find the closest nCV value to target_val
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='give colvar_files', action='store', dest='f', type=str)
    parser.add_argument('-v', help='tnCV', action='store', dest='v', type=float)
    parser.add_argument('-wm', help="""give min THEN max value of cv,
                        in case of adding windows this will change a little bit the behavior
                        but this should not be problematic""",
                        nargs='*', action='store', dest='w', type=float)
    args = parser.parse_args()

    #read COLVAR.txt
    df = plumed_pandas.read_as_pandas(args.f)

    #extract time of frame closest to the target_value and cnt of gksi at this frame
    if args.v == args.w[0] or args.v == args.w[1]:
        tvalue = df.iloc[(df['nCV']-float(args.v)).abs().idxmin()]['time']
        gksi = df.iloc[(df['nCV']-float(args.v)).abs().idxmin()]['restraint2.RMSDMID_cntr']
        gkap = df.iloc[(df['nCV']-float(args.v)).abs().idxmin()]['restraint2.RMSDMID_kappa']
    else:
        #get idx of min diff btw args.v & restraint.ncv_cnt for all the pulling simulation in df
        mindiff = (df['restraint.nCV_cntr']-float(args.v)).abs().idxmin()
        tvalue = df.iloc[mindiff]['time']
        gksi = df.iloc[mindiff]['restraint2.RMSDMID_cntr']
        gkap = df.iloc[mindiff]['restraint2.RMSDMID_kappa']

    return (tvalue, gksi, gkap)

# If called from the command line...
if __name__ == "__main__":
    print(main())
