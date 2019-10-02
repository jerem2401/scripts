#!/usr/bin/env python

import plumed_pandas
import numpy as np
import pandas as pd 
import argparse

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='give colvar_file', nargs='1', action='store', dest='f')        

    df=plumed_pandas.read_as_pandas(args.f)
    df = df[['nCV','guide_restraint.bias']]

    bstep=(df['nCV'].max()-df['nCV'].min())/200
    print(df['nCV'].max(), df['nCV'].min())
    print(bstep)
    print(round(bstep,6))
    bval=[]
    for i in range(1,201):
        bval.append(round(df['nCV'].min()+bstep*i,6))
    print(bval)
    hist=pd.DataFrame(bval)
    hist.columns = ['bval']
    hist['counts'] = 0
    for c, (i, j) in enumerate(zip(df['nCV'], df['guide_restraint.bias'])):
        if c in [i for i in range(1000,500000,1000)]:
            print(c)
        for index, row in hist.iterrows():
            if index == 0:
                N=df['nCV'].min()
            if N <= i < row['bval']:
                hist.loc[index,'counts'] = hist.loc[index,'counts'] + np.exp(j/2.5)
            N=row['bval']
    #print(hist.loc[0,'counts'])
    display(hist)
