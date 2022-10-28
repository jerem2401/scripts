#!/usr/bin/env python
'''
Author: Jeremy Lapierre <jeremy.lapierre@uni-saarland.de>
Description:

'''
import sys
import numpy as np
from itertools import dropwhile, islice, takewhile
import matplotlib.pyplot as plt
import matplotlib as mpl
import os


def log_parser(file):
    """useful function when average probabilities of exchange are not printed, for example in a
    continued run that was actually finished"""

    time = []
    ex = []

    with open(file, 'r') as reader:
        start = dropwhile(lambda x: 'Command line:' not in x, reader)
        end = takewhile(lambda x: 'GROMACS version:' not in x, start)
        for line in islice(end, 1, 2):
            tmp = line.split(' ')
        tmp = [x for x in tmp if 'E' in x or 'TMP' in x]
        ksis = [x.strip('E_').strip('TMP_').strip('/').strip('./E_').strip('./TMP_') for x in tmp]

        start = dropwhile(lambda x: 'Replica exchange at step 5700' not in x, reader)
        print('in parsing function, lets parse')
        for line in start:
            if 'Replica exchange at step ' in line:
                time.append(line.split()[-1])
            if 'Repl pr ' in line:
                tmp = line.strip('Repl pr')
                ex.append(tmp.split())
    print('I got all proba from file')
#    if len(time) != len(ex):
#        sys.exit('exiting, time not the same length as ex')

    time = [float(x) for x in time]

    len_list = [len(x) for x in ex]
    len_list = list(set(len_list))
    if len(len_list) > 2:
        sys.exit('some proba arrays are not complete, exiting')
    ex1 = []
    ex2 = []
    for i in ex:
        if len(i) == len_list[0]:
            ex1.append(i)
        else:
            ex2.append(i)
    #n = 1
    #for i in ex:
    #    if n % 2 == 1:
    #        ex1.append(i)
    #    else:
    #        ex2.append(i)
    #    n += 1

    print('proba are split')

    for l in ex1:
        for k in range(0, len(l)):
            l[k] = float(l[k])

    for l in ex2:
        for k in range(0, len(l)):
            l[k] = float(l[k])

    if len(ksis) != len(ex1[0])+len(ex2[0])+1:
        sys.exit('exiting, len(ksis) should be equal to len(ex1)+len(ex2)+1')

    #ex1 = [ex1[i] for i in range(0,len(ex1)-1,10)]
    #ex2 = [ex2[i] for i in range(0,len(ex2)-1,10)]

    ex1_ar = np.asarray(ex1)
    ex2_ar = np.asarray(ex2)

    print('computing mean')

    mex1 = np.mean(ex1_ar, axis=0)
    mex2 = np.mean(ex2_ar, axis=0)

    print('I m done computing means')

    fex = [item for sublist in zip(mex1, mex2) for item in sublist]
    if len(mex1) != len(mex2):
        fex.append(mex1[-1])
    if len(ksis) != len(fex)+1:
        sys.exit('exiting, len(ksis) should be equal to len(ex1)+len(ex2)+1')

    print('quiting parsing function.')
    return time, fex, ksis


def main(n, files, getavrp):
    n = int(n)
    getavrp = int(getavrp)
    times = []
    fexs = []
    ksiss = []
    for file in files:
        if getavrp == 1:
            print(f"starting parsing {file}")
            a = log_parser(file)
            times.append(a[0])
            fexs.append(a[1])
            ksiss.append(a[2])
            print(f"finished parsing {file}")
        else:
            with open(file, 'r') as reader:
                start = dropwhile(lambda x: 'Command line:' not in x, reader)
                end = takewhile(lambda x: 'GROMACS version:' not in x, start)
                for line in islice(end, 1, 2):
                    tmp = line.split(' ')
                tmp = [x for x in tmp if 'E' in x or 'TMP' in x]
                ksis = [x.strip('E_').strip('TMP_').strip('/').strip('./E_').strip('./TMP_') for x in tmp]
                ksiss.append(ksis)

                start = dropwhile(lambda x: 'Repl  average probabilities' not in x, reader)
                end  = takewhile(lambda x: 'Repl  number of exchanges:' not in x, start)
                for line in islice(end, 1, None):
                    tmp.append(line)
                tmp = tmp[1].split()[1::]
                fexs.append([float(i) for i in tmp])

    ksiss = [[int(float(x) * 10**2) / 10.0**2 for x in item] for item in ksiss]
    print(fexs, ksiss, sep='\n\n')
    mpl.rcParams.update(plt.rcParamsDefault)
    mpl.rcParams['figure.dpi'] = 300
    mpl.rcParams['axes.linewidth'] = 0.1*n
    mpl.rcParams['grid.linewidth'] = 0.05*n
    mpl.rcParams['lines.linewidth'] = 0.33*n
    mpl.rcParams['lines.markersize'] = 1*n
    mpl.rcParams['lines.markeredgewidth'] = 0.2*n
    mpl.rcParams['xtick.major.size'] = mpl.rcParams['ytick.major.size'] = 0.6*n
    mpl.rcParams['xtick.major.width'] = mpl.rcParams['ytick.major.width'] = 0.1*n
    mpl.rcParams['axes.labelpad'] = 4*n
    mpl.rcParams['xtick.major.pad'] = mpl.rcParams['ytick.major.pad'] = 1.2*n
    mpl.rcParams['mathtext.default'] = 'regular'
    plt.rcParams.update({'font.sans-serif':'Helvetica'})
    mpl.rcParams['font.size'] = 0.1*n
    cm = 1/2.54
    label_size = 4*n
    tick_size = 4*n
    fontax = label_size+3

    plt.close()

    fig, axs = plt.subplots(figsize=(18*cm*n,4*cm*n))
    axs.set_xlabel('dz (nm)', fontsize=fontax)
    axs.set_ylabel('Probability of\nexchange acceptance', fontsize=fontax)

    n = 0
    for x, y in zip(ksiss, fexs):
        name = files[n].strip('.log')
        axs.plot(x, [0]+y, label=name, linestyle='-', marker='o')
        #axs.plot(x, [0.18]*(len(y)+1))
        #axs.set_xticks(x)
        n+=1
    x_tic = [item for items in ksiss for item in items]
    x_tic = list(set(x_tic))
    print(x_tic, min(x_tic), max(x_tic))
    axs.plot([min(x_tic),max(x_tic)],[0.25,0.25], color='black')
    axs.plot([min(x_tic),max(x_tic)],[0.47,0.47], color='black')
    axs.set_xticks(x_tic)
    axs.set_xlim([min(x_tic)-0.05, max(x_tic)+0.05])
    plt.xticks(fontsize=tick_size*0.55, rotation=90)
    plt.yticks(fontsize=tick_size)
    plt.savefig('./all.png', dpi=150, bbox_inches="tight")
    plt.close()
    print(f"finished ploting all.png")

    n = 0
    for x,y in zip(ksiss, fexs):
        fig = plt.figure(figsize=(40,6))
        x = [int(float(f) * 10**2) / 10.0**2 for f in x]
        plt.plot(x, [0]+y)
        plt.plot(x, [0.18]*(len(y)+1))
        plt.xticks(x)
        name = files[n].strip('.log')
        plt.savefig(name+'.png')
        plt.close()
        n+=1
        print(f"finished ploting {name}.png")

# If called from the command line...
if __name__ == "__main__":
    """ex: getavgP_andplot.py 100 mid*.log 1 -> getavrp of mid*.log and plot
           getavgP_andplot.py 100 mid*.log 0 ->  just plot"""
    n = sys.argv[1]
    files = sys.argv[2:-1]
    getavrp = sys.argv[-1]
    main(n, files, getavrp)
