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
        for line in start:
            if 'Replica exchange at step ' in line:
                time.append(line.split()[-1])
            if 'Repl pr ' in line:
                tmp = line.strip('Repl pr')
                ex.append(tmp.split())

    if len(time) != len(ex):
        sys.exit('exiting, time not the same length as ex')

    time = [float(x) for x in time]

    ex1 = []
    ex2 = []
    n = 1
    for i in ex:
        if n % 2 == 1:
            ex1.append(i)
        else:
            ex2.append(i)
        n += 1

    for l in ex1:
        for k in range(0, len(l)):
            l[k] = float(l[k])

    for l in ex2:
        for k in range(0, len(l)):
            l[k] = float(l[k])

    if len(ksis) != len(ex1[0])+len(ex2[0])+1:
        sys.exit('exiting, len(ksis) should be equal to len(ex1)+len(ex2)+1')

    ex1_ar = np.asarray(ex1)
    ex2_ar = np.asarray(ex2)

    mex1 = np.mean(ex1_ar, axis=0)
    mex2 = np.mean(ex2_ar, axis=0)

    fex = [item for sublist in zip(mex1, mex2) for item in sublist]
    if len(mex1) != len(mex2):
        fex.append(mex1[-1])
    if len(ksis) != len(fex)+1:
        sys.exit('exiting, len(ksis) should be equal to len(ex1)+len(ex2)+1')

    return time, fex, ksis


def main(n, files, getavrp):
    n = int(n)
    getavrp = int(getavrp)
    times = []
    fexs = []
    ksiss = []
    for file in files:
        if getavrp == 1:
            times.append(log_parser(file)[0])
            fexs.append(log_parser(file)[1])
            ksiss.append(log_parser(file)[2])
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

    print(fexs, ksiss, sep='\n\n')
    mpl.rcParams.update(plt.rcParamsDefault)
    mpl.rcParams['font.size'] = 0.1*n
    tick_size = 1.7*n
    plt.xticks(fontsize=tick_size)
    plt.yticks(fontsize=tick_size)

    plt.close()

    fig, axs = plt.subplots(figsize=(40,6))

    n = 0
    for x, y in zip(ksiss, fexs):
        name = files[n].strip('.log')
        x = [int(float(f) * 10**2) / 10.0**2 for f in x]
        axs.plot(x, [0]+y, label=name)
        axs.plot(x, [0.18]*(len(y)+1))
        axs.set_xticks(x)
        n+=1
    axs.legend()
    plt.savefig('./all.png')
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
