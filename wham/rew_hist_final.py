#!/usr/bin/env python

import plumed_pandas
import numpy as np
import pandas as pd 
import argparse
import re

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='give colvar_file, for building weithed histo', action='store', dest='f')        
    args = parser.parse_args()

    df=plumed_pandas.read_as_pandas(args.f)
    df = df[['nCV','guide_restraint.bias']]
    
    hist2, bin_edges = np.histogram(df['nCV'], bins=[i for i in np.around(np.arange(-1.0,1.0,0.01), decimals=3)], weights=[np.exp(i/2.49434) for i in df['guide_restraint.bias']],density=True)
    
    d = {'z': bin_edges[0:-1],'hist': hist2}
    
    hdf=pd.DataFrame(d)

    pos=re.search('(-|[0-9]|\.)+(?=.txt)', args.f).group()    

    np.savetxt(r'./histo_'+str(pos)+'_.dat', hdf.values, header="col1=z col2=hist\n#1 #2 "+str(pos)+"\n#1 #2 5500", fmt='%.6f')

# If called from the command line... 
if __name__ == "__main__":
    main()
