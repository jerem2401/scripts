#!/usr/bin/env python

import pandas as pd
import plumed_pandas
import argparse
import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='give 2dpmf_clean.txt file', action='store', dest='f', type=str)
    parser.add_argument('-o', help='give path for projected data output', action='store', dest='o')
    args = parser.parse_args()

    #open 2dwham output
    df = pd.read_table(args.f, delim_whitespace=True, dtype={'#X': np.float64, 'Y': np.float64, 'Free': np.float64, 'Pro': np.float64})

    #computing Prob from free nrg because not enought decimals in Pro column
    df['Pro2'] = np.exp((df['Free'])/-2.494)

    #Projection on nCV: adding all proba in same nCV values
    aggregation_functions = {'Pro2': 'sum'}
    df_new = df.groupby(df['#X']).aggregate(aggregation_functions)
    df_new=df_new.reset_index()

    #recomputing free from Prob
    df_new['Free2'] = -2.494*np.log(df_new['Pro2'])

    #save projected Free_nrg and Prob
    if args.o == None:
        np.savetxt(r'./1dpmf.txt', df_new.values, header="nCV   Prob   Free", fmt='%.3f')
    else:
        np.savetxt(args.o, df_new.values, header="nCV   Prob   Free", fmt='%.3f')

# If called from the command line...
if __name__ == "__main__":
    main()
