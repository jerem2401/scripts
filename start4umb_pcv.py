#!/usr/bin/env python

import plumed_pandas

def main(colvf, ksi, minmax):
    '''
    parsing the colvar.txt, find the closest nCV value to target_val
    '''

    #read COLVAR.txt
    df = plumed_pandas.read_as_pandas(colvf)

    #extract time of frame closest to the target_value and cnt of gksi at this frame
    if ksi == minmax[0] or ksi == minmax[1]:
        mindiff = (df['p1.sss']-float(ksi)).abs().idxmin()
        found = df.iloc[mindiff]['p1.sss']
        tvalue = df.iloc[mindiff]['time']
    else:
        #get idx of min diff btw args.v & path.p1.sss_cntr for all the pulling simulation in df
        #mindiff = (df['path.p1.sss_cntr']-float(ksi)).abs().idxmin()
        mindiff = (df['p1.sss']-float(ksi)).abs().idxmin()
        found = df.iloc[mindiff]['p1.sss']
        tvalue = df.iloc[mindiff]['time']

    return found, tvalue

# If called from the command line...
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='give colvar_files', action='store', dest='f', type=str)
    parser.add_argument('-v', help='tnCV', action='store', dest='v', type=float)
    parser.add_argument('-wm', help="""give min THEN max value of cv,
                        in case of adding windows this will change a little bit the behavior
                        but this should not be problematic""",
                        nargs='*', action='store', dest='w', type=float)
    args = parser.parse_args()

    F, TV = main(args.f, args.v, args.w)
    print(F, TV)
