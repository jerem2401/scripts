#!/usr/bin/env python

import sys
import numpy as np
import pandas as pd

def rdata(colvar):
    '''
    parsing the colvar.txt, extracting RMSDs and angles data
    '''
    
    time = []
    nCV = []
    RMSDeq = []
    RMSDax = []
    phi=[]
    psi=[]
    lwall_bias=[]
    nCVcnt=[]
    
    #read file and extract RMSDs and angle data
    with open(colvar) as f:
        fl = f.readline()
        fls=fl.split(' ')
        fls = [x.replace("\n", "") for x in fls]
        fls=[x.replace('.', '_') for x in fls]
        fls=fls[2:] #get all the names of the values printed by plumed from the first line
        #print(fls)
        timeid=fls.index('time')
        #nCVid=fls.index('nCV')
        RMSDaxid=fls.index('RMSDAX')
        RMSDeqid=fls.index('RMSDEQ')
        phiid=fls.index('phi')
        psiid=fls.index('psi')
        #lwall_biasid=fls.index('lwall.bias') #get indexes of all these names, will be used to allocate column indexes to object
        #print(timeid, nCVid, RMSDaxid, RMSDeqid, phiid, psiid)
        for line in f:
            t=line.split()
            if t[0] == '#!':
                next(f)#jump the header to go to data
            else:
                time.append(round(float(t[timeid]), 2))
                try:
                    nCVid=fls.index('nCV')
                    nCV.append(float(t[nCVid]))
                except ValueError:
                    pass
                try:
                    nCVcntid=fls.index('restraint_nCV_cntr')
                    nCVcnt.append(float(t[nCVcntid]))
                except ValueError:
                    pass
                try:
                    lwall_biasid=fls.index('lwall_bias')
                    lwall_bias.append(float(t[lwall_biasid]))
                except ValueError:
                    pass
                RMSDax.append(float(t[RMSDaxid]))
                RMSDeq.append(float(t[RMSDeqid]))
                phi.append(float(t[phiid]))
                psi.append(float(t[psiid]))#each line are split and each element are assigned to corresponding object
    f.close()
    return(RMSDeq, RMSDax, phi, psi, lwall_bias, nCV, nCVcnt)
