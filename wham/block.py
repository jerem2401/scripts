#!/usr/bin/env python

import pandas as pd
import plumed_pandas
import argparse
import os
import numpy as np
import re
import glob
from wham_pmf import whamloop
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib


def rew_chunk(f,c,k):
    #opening colvar files, chunking them, do rew_histo for each chunk
    for i in f:
        pos=re.search('(?<=colvar).*(?=_)', i).group()
        check=glob.glob('./histo_'+str(pos)+'_c*.dat')
        if check != []:
            print('./histo_'+str(pos)+'_c*.dat  exists !')
        else:
            print('reading '+i)
            df=plumed_pandas.read_as_pandas(i)
            dflen=len(df.index)
            print('dflen='+str(dflen))
            chlen=round(dflen/c)
            print('chlen='+str(chlen))

            for n,j in enumerate(range(0,dflen,chlen)):
                chunk=df[['nCV','guide_restraint.bias']][j:j+chlen]

                print('head of chunk '+str(n)+' from '+str(j)+' to '+str(j+chlen)+' :', chunk.head(), chunk.tail(), sep='\n')

                #Put next lines as a function in rew_histo, to make it versatile
                hist2, bin_edges = np.histogram(chunk['nCV'], bins=[i for i in np.around(np.arange(-1.25,1.25,0.01), decimals=3)], weights=[np.exp(i/2.49434) for i in chunk['guide_restraint.bias']],density=True)
                hist2 = hist2 / hist2.sum()
                d = {'z': bin_edges[0:-1],'hist': hist2}
                hdf=pd.DataFrame(d)

                # /!\ saved in working directory
                np.savetxt(r'./histo_'+str(pos)+'_c'+str(n)+'.dat', hdf.values, header="col1=z col2=hist\n#1 #2 "+str(pos)+"\n#1 #2 "+str(k), fmt='%.6f')
                print('./histo_'+str(pos)+'_c'+str(n)+'.dat is done')
    return(dflen,chlen)



def wham_chunk(c):
    #do wham with the whamloop function of wham_pmf.py
    print('here comes the wham loop')

    for i in range(c):
        histo_files=glob.glob('hist*_c'+str(i)+'.dat')
        whamloop(histo_files)

        os.rename('test_wham_z.out', 'test_wham_z_'+str(i)+'.out')
        os.rename('test_wham_pmf.out', 'test_wham_pmf_'+str(i)+'.out')

        print('pmf of chunk '+str(i)+' is done')



def block_avg(f2):
    #computing block avergaes. std and sem plotting the final block_averaged pmf with error bars
    print('Here comes the block averaging')

    files=glob.glob(f2)
    files.sort()
    lf=len(files)
    ref=files[round(lf/2)]
    ksir, pror, freer = np.loadtxt(ref, unpack=True, delimiter=' ')

    min_freer=np.amin(freer)

    allf=[]
    for i in files:
        ksi, pro, free = np.loadtxt(i, unpack=True, delimiter=' ')
        #delta=min_freer-free[np.where(freer == min_freer)]
        #shift every pmf to the ref (min of middle block)
        #free=free+delta
        #shift ervery pmf so the ref free nrg point is shifted to 0
        #free=free-min_freer
        free=free-free[np.where(freer == min_freer)]
        allf.append(free)

    arr=np.vstack(allf)

    avg=np.average(arr, axis=0)
    #std=np.std(arr, axis=0, ddof=1)
    sem=stats.sem(arr, axis=0)

    fp = open('block_pmf.out', 'w')
    for i in range(len(ksir)):
        print(str(ksir[i])+'    '+str(avg[i])+'    '+str(sem[i]), file=fp)
    fp.close()
    np.savetxt(r'block_pmf.out', np.c_[ksir,avg,sem],fmt='%.3f',delimiter='    ')

    return(ksir,avg,sem)

def plot_bpmf(ksir,avg,sem,c):
    #plotting block_pmf, todo: make plotpmf function in wham_pmf more felixible so it can be used here
    print('Here comes the plotting of the block_pmf.jpeg')

    matplotlib.rcParams.update({'font.size': 25})
    
    fig = plt.figure(figsize=(12, 8))
    ax1 = fig.add_subplot(1, 1, 1)
    
    ax1.errorbar(ksir, avg, yerr=2*sem)
    ax1.set_xlabel(r'$\xi_{new}$')
    ax1.set_ylabel('<G> (with '+str(c)+' blocks), kJ/mol')
    
    plt.savefig("block_pmf.jpeg")
    plt.close()

def main(raw_args=None):
    '''defining what to do when script called from the command line or main() called from another script'''

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='give colvar_file of each umb windows (give posta colvar files if posta as been performed)', action='store', nargs='*', dest='f')
    parser.add_argument('-c', help='give nbr of chunk you want', action='store', dest='c', type=int)
    parser.add_argument('-k', help='force constant used for each windows', action='store', dest='k', required=True)
    args = parser.parse_args(raw_args)
    
    dflen, chlen = rew_chunk(args.f,args.c,args.k)
    wham_chunk(args.c)
    ksir, avg, sem = block_avg('test_wham_pmf_*')
    plot_bpmf(ksir,avg,sem,args.c)

if __name__ == "__main__":
    '''if called from the command line do both functions, using argparse'''
    main()
