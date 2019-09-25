#!/usr/users/jlapier/.conda/envs/env1/bin/python

import sys
import numpy as np
import pandas as pd

def main():
    '''
    parsing the colvar.txt, extracting RMSDs and angles data, find the closest nCV value to target_val
    '''
    
    time = []
    nCV = []
    
    #read file and extract RMSDs and angle data
    with open(sys.argv[1]) as f:
        fl = f.readline()
        fls=fl.split(' ')
        fls = [x.replace("\n", "") for x in fls]
        fls=[x.replace('.', '_') for x in fls]
        fls=fls[2:]
        #print(fls)
        timeid=fls.index('time')
        nCVid=fls.index('nCV')
        #print(timeid, nCVid)
        for line in f:
            t=line.split()
            if t[0] == '#!':
                next(f)
            else:
                time.append(round(float(t[timeid]), 2))
                nCV.append(float(t[nCVid]))
    
    #extract time step of frame closest to the target_value
    time=pd.DataFrame(time)
    nCV=pd.DataFrame(nCV)
    time.columns = ['time']
    nCV.columns = ['nCV']
    df=time.join(nCV)
    #neareast = df.iloc[(df['nCV']-float(sys.argv[2])).abs().idxmin()]
    value=df.iloc[(df['nCV']-float(sys.argv[2])).abs().idxmin()]['time']
    return value

# If called from the command line... 
if __name__ == "__main__":
    print(main())
