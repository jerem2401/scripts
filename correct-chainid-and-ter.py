#!/usr/bin/env python
#script given by Jochen Hub

import sys, string

try:
        f=open(sys.argv[1],"r")
except:
        print('Usage:',sys.argv[0],'pdbfile')

labels='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz123456789'
labels=labels+labels+labels+labels+labels+labels+labels+labels+labels

last=-1
chainnr=0

for line in f:
    if (not (line[0:6]=='ATOM  ' or line[0:6]=='HETATM')):
        print(line, end = '')
    else:
        res     = int(line[22:26])
        # chainID = labels[chainnr]
        
        if (res<last):
            print('TER')
            chainnr += 1
        last = res

        line = line[0:21] + labels[chainnr] + line[22:]
        print(line,end='')
