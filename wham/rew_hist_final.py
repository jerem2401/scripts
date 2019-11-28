#!/usr/bin/env python

import plumed_pandas
import numpy as np
import pandas as pd 
import argparse
import re
import sys

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='give colvar_file, for building weithed histo', action='store', dest='f', required=True)        
    parser.add_argument('--plot', help='plot histograms', action='store_true', dest='plot')
    parser.add_argument('-k', help='force constant used for each windows', action='store', dest='k', required=True) 
    args = parser.parse_args()


    if not args.plot:
        import os
        pos=re.search('(-|[0-9]|\.)+(?=_)', args.f).group()    
        
        if os.path.exists('./histo_'+str(pos)+'_.dat'):
            print('./histo_'+str(pos)+'_.dat exists !')
            sys.exit()
        else:
            df=plumed_pandas.read_as_pandas(args.f)
            df = df[['nCV','guide_restraint.bias']]
            
            hist2, bin_edges = np.histogram(df['nCV'], bins=[i for i in np.around(np.arange(-1.25,1.25,0.01), decimals=3)], weights=[np.exp(i/2.49434) for i in df['guide_restraint.bias']],density=True)
            
            d = {'z': bin_edges[0:-1],'hist': hist2}
            
            hdf=pd.DataFrame(d)
            
            # /!\ saved in working directory
            print('./histo_'+str(pos)+'_.dat is done')
            np.savetxt(r'./histo_'+str(pos)+'_.dat', hdf.values, header="col1=z col2=hist\n#1 #2 "+str(pos)+"\n#1 #2 "+str(args.k), fmt='%.6f')

    if args.plot:

        import matplotlib.pyplot as plt
        import os       
        import matplotlib
        matplotlib.rcParams.update({'font.size': 23})

 
        histo_files = os.popen('ls histo* | sort -t _ -k 2 -n').read().split()
        
        fig=plt.figure(figsize=(15,8))
        ax = fig.add_subplot(1, 1, 1)
        
        for h in histo_files:
        
            z, thishist = np.loadtxt(h, unpack=True)
            #ax.bar(z, thishist, width = z[0]-z[1])
            ax.plot(z, thishist)

        plt.xlim(min(z), max(z))
        plt.savefig('hist.jpeg')

# If called from the command line... 
if __name__ == "__main__":
    main()
