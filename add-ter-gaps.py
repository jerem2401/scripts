#!/usr/bin/env python
#script given by Jochen Hub

import sys, string

try:
        f=open(sys.argv[1],"r")
except:
        print('Usage:',sys.argv[0],'pdbfile')

last=-1

for line in f:
    if (not (line[0:6]=='ATOM  ' or line[0:6]=='HETATM')):
        print(line,end='')
    else:
        res     = int(line[22:26])
        
        if (res>last+1):
            print('TER')
        last = res

        print(line,end='')
