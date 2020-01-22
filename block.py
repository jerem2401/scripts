#!/usr/bin/env python

import pandas as pd
import plumed_pandas
import argparse
import subprocess
import os
import numpy as np
import re

parser = argparse.ArgumentParser()
parser.add_argument('-f', help='give colvar_file of each umb windows (give posta colvar files if posta as been performed)', action='store', nargs='*', dest='f')
parser.add_argument('-c', help='give nbr of chunk you want', action='store', dest='c', type=int)
parser.add_argument('-k', help='force constant used for each windows', action='store', dest='k', required=True)
args = parser.parse_args()

#opening colvar files, chunking them, do rew_histo for each chunk
for i in args.f:
    print('reading '+i)
    df=plumed_pandas.read_as_pandas(i)
    dflen=len(df.index)
    print('dflen='+str(dflen))
    chlen=round(dflen/args.c)
    print('chlen='+str(chlen))
    for n,k in enumerate(range(0,dflen,chlen)):
        chunk=df[['nCV','guide_restraint.bias']][k:k+chlen]

        print('head of chunk '+str(n)+' from '+str(k)+' to '+str(k+chlen)+' :', chunk.head(), chunk.tail(), sep='\n')

        #Put next lines as a function in rew_histo, to make it versatile
        hist2, bin_edges = np.histogram(chunk['nCV'], bins=[i for i in np.around(np.arange(-1.25,1.25,0.01), decimals=3)], weights=[np.exp(i/2.49434) for i in chunk['guide_restraint.bias']],density=True)
        hist2 = hist2 / hist2.sum()
        d = {'z': bin_edges[0:-1],'hist': hist2}
        hdf=pd.DataFrame(d)
        
        # /!\ saved in working directory
        pos=re.search('(?<=colvar).*(?=_)', i).group()
        np.savetxt(r'./histo_'+str(pos)+'_c'+str(n)+'.dat', hdf.values, header="col1=z col2=hist\n#1 #2 "+str(pos)+"\n#1 #2 "+str(args.k), fmt='%.6f')
        print('./histo_'+str(pos)+'_c'+str(n)+'.dat is done')

#do wham, put next lines as function in wham.pmf, to make it versatile
print('here comes the wham loop')
for n,k in enumerate(range(0,dflen,chlen)):
    opt='-f histo*_c'+str(n)+'.dat'
    subprocess.call(['wham.py', opt])
    os.rename('test_wham_z.out', 'test_wham_z_'+str(n)+'.out')
    os.rename('test_wham_pmf.out', 'test_wham_pmf_'+str(n)+'.out')

    print('pmf of chunk '+str(n)+' is done')
