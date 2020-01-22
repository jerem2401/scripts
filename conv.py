#!/usr/bin/env python

import pandas as pd
import plumed_pandas
import argparse
import matplotlib.pyplot as plt
import matplotlib
import os

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='give colvar_file, for plotting f(var)=t', action='store', nargs='*', dest='f')
    parser.add_argument('-v', help='give colvar_file, for building weithed histo', action='store', nargs='*', dest='var')
    args = parser.parse_args()

    args.var.append('time')

    for i in args.f:
        df=plumed_pandas.read_as_pandas(i)
        df=df.drop(list(set(df.columns)-set(args.var)), axis='columns')
        fig=plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        for j in args.var[0:-1]:
            ax.plot(df[args.var[-1]][::1000], df[j][::1000], label=j)
            ax.legend()
        plt.savefig(os.path.basename(i).split('.txt')[0]+'_conv.jpeg', quality=30)
        print(os.path.basename(i).split('.txt')[0]+'_conv.jpeg done')

# If called from the command line... 
if __name__ == "__main__":
    main()
