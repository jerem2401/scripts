#!/usr/bin/env python

import argparse
import sys
import os
import numpy as np
import pandas as pd
import plumed_pandas

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='give colvar_file, for building weithed histo', action='store', dest='f', required=True)
    parser.add_argument('--plot', help='plot histograms together: all1 (all2 for separating plots in the middle of the CV in case of cyclic CV), or separately: sep', action='store', dest='plot')
    parser.add_argument('-k', help='force constant used for each windows', action='store', dest='k', required=True)
    parser.add_argument('-s', help='binning step, deflt=0.01 for ext -1.25 and 1.25', action='store', dest='s', default=0.01, type=float)
    parser.add_argument('-min', help='min boundary for hist, deflt= %(default)s', action='store', dest='min', default=-1.25, type=float)
    parser.add_argument('-max', help='max boundary for hist, deflt= %(default)s', action='store', dest='max', default=1.25, type=float)
    parser.add_argument('-rew', help='if present, allows to reweight histo by the value of the guide_restraint.bias', action='store_true', dest='rew')
    parser.add_argument('-col', help='from which column of the colvar file do you want to do an histogram, deflt= %(default)s', default='nCV', action='store', dest='col', type=str)
    parser.add_argument('--o', help='name of hist file', action='store', dest='o', type=str)
    parser.add_argument('-pos', help='center of haromic potential', action='store', dest='pos')
    parser.add_argument('-nopld', help='tell taht file are not in plumed format',
                        action='store_true', dest='nopld')
    args = parser.parse_args()

    print('careful: put the -rew option added in the last maj for reweighting !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    print(args.f, args.o)
    if args.plot == None:
        pos=args.pos
        #pos=re.search('(-|[0-9]|\.)+(?=_)', args.f).group()

        if os.path.exists('./histo_'+str(pos)+'_.txt'):
            print('./histo_'+str(pos)+'_.txt exists !')
            sys.exit()
        elif args.o != None and os.path.exists(args.o):
            print(args.o+' exists !')
            sys.exit()
        else:
            if not args.nopld:
                df=plumed_pandas.read_as_pandas(args.f)
            else:
                df=pd.read_csv(args.f, sep=' ', header=None)
                df.rename(columns={ df.columns[1]: "d.z" }, inplace = True)

            if not args.rew:
                df = df[args.col]
                hist2, bin_edges = np.histogram(df, bins=[i for i in np.around(np.arange(args.min,args.max,args.s), decimals=3)],density=True)
                hist2 = hist2 / hist2.sum()
            else:
                df = df[['dz', 'zrest.bias']]
                hist2, bin_edges = np.histogram(df['dz'], bins=[i for i in np.around(np.arange(args.min,args.max,args.s), decimals=3)], weights=[np.exp(i/2.49434) for i in df['zrest.bias']],density=True)
                hist2 = hist2 / hist2.sum()

            d = {'z': bin_edges[0:-1],'hist': hist2}

            hdf = pd.DataFrame(d)

            # /!\ saved in working directory
            if args.o == None:
                np.savetxt(r'./histo_'+str(pos)+'_.txt', hdf.values, header="col1=z col2=hist\n#1 #2 "+str(pos)+"\n#1 #2 "+str(args.k), fmt='%.6f')
                print('./histo_'+str(pos)+'_.txt is done')
            else:
                np.savetxt(args.o, hdf.values, fmt='%.6f')
                print(args.o+' is done')

    if args.plot == 'all2':

        import matplotlib.pyplot as plt
        import matplotlib

        matplotlib.rcParams.update({'font.size': 23})


        histo_files = os.popen('ls histo*.txt | sort -t _ -k 2 -n').read().split()
        h1 = histo_files[:len(histo_files)//2]
        h2 = histo_files[len(histo_files)//2:]

        for i in [h1,h2]:

            fig=plt.figure(figsize=(15,8))
            ax = fig.add_subplot(1, 1, 1)

            for j in i:

                z, thishist = np.loadtxt(j, unpack=True)
                #ax.bar(z, thishist, width = z[0]-z[1])
                ax.plot(z, thishist)

            plt.xlim(min(z), max(z))
            plt.savefig('hist_'+str([h1, h2].index(i))+'.jpeg')
            plt.close()

    if args.plot == 'all1':
        import statistics
        import matplotlib.pyplot as plt
        import matplotlib
        matplotlib.rcParams.update({'font.size': 15})
        histo_files = os.popen('ls histo*.txt | sort -t _ -k 2 -n').read().split()
        fig = plt.figure(figsize=(15,8))
        fig, axs = plt.subplots(figsize=(15,8))
        n = 1
        allz = []
        for i in  histo_files:
            label = str(n)
            xs = i.split('_')[1]
            x = float(xs)
            k = os.popen(f"grep -oP '(?<=KAPPA=)[0-9]*$' ../E_{xs}/plumed_{xs}.dat").read()
            z, thishist = np.loadtxt(i, unpack=True)
            y = max(thishist)
            axs.plot(z, thishist, label=label+":"+str(x)+':'+k)
            axs.annotate(s=label, xy=(x, y), fontsize=8)
            n+=1
            allz = np.concatenate((allz, z), axis=None)
        axs.legend(loc='upper left', ncol=7, fontsize=7, handlelength=0.4)
        plt.xlim(np.amin(allz), np.amax(allz))
        plt.savefig('hist_all1.jpeg')
        plt.close()


    if args.plot == 'sep':

        import matplotlib.pyplot as plt
        import matplotlib

        histo_files = os.popen('ls histo* | sort -t _ -k 2 -n').read().split()

        for h in histo_files:

            fig, ax = plt.subplots(figsize=(15,8))

            z, thishist = np.loadtxt(h, unpack=True)
            Eg = h.split('_')[2]
            Eg = Eg.strip('.txt')

            ax.plot(z, thishist)
            ax.plot(np.full((2,), float(Eg)), np.asarray([min(thishist), max(thishist)]))

            plt.title(h.strip('.txt'))
            plt.savefig(h.strip('.txt')+'.jpeg')

# If called from the command line...
if __name__ == "__main__":
    main()
