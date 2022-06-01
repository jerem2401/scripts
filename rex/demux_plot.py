#!/usr/bin/env python
'''
Author: Jeremy Lapierre <jeremy.lapierre@uni-saarland.de>
Description:

'''
import pandas as pd
from itertools import dropwhile, islice, takewhile
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse


def demux_parser(path):
    #get windows ksis
    with open(path+'/md.log', 'r') as reader:
        start = dropwhile(lambda x: 'Command line:' not in x, reader)
        end = takewhile(lambda x: 'GROMACS version:' not in x, start)
        for line in islice(end, 1, 2):
            tmp = line.split(' ')
        tmp = [x for x in tmp if 'E' in x or 'TMP' in x]
        ksis = [x.strip('E_').strip('TMP_').strip('/').strip('./E_').strip('./TMP_') for x in tmp]

#   make df from replica_temp.xvg
    df = pd.read_csv(path+'/replica_temp.xvg', delim_whitespace=True)
    df = pd.concat([df.columns.to_frame().T, df])
    df.columns = ['Time']+list(range(0, len(ksis)))
    return df, ksis


def plot(n, curve_per_plot, col, skip, df, ksis, path):
    mpl.rcParams.update(plt.rcParamsDefault)
#   style.use('dark_background')
    mpl.rcParams['font.size'] = 0.1*n
    tick_size = 1.7*n
    plt.xticks(fontsize=tick_size)
    plt.yticks(fontsize=tick_size)

    plt.close()
#   for 56 curves we typically need figsize=(35,70),
#   the next line evaluate figsize on this basis
    yfigsiz = (len(ksis)/56)*70

    rest = (len(ksis)/curve_per_plot) % col
    if rest != 0:
        x = (len(ksis)/curve_per_plot)//col + 1
        fig, axs = plt.subplots(int(x), col, figsize=(35, yfigsiz))
    else:
        x = 0
        fig, axs = plt.subplots((len(ksis)/curve_per_plot)/col, col, figsize=(35, yfigsiz))

    axs = axs.flatten()

    dic = dict(zip(list(range(0, len(ksis))), ksis))

    for i in range(0, len(ksis), curve_per_plot):
        print(f'making {curve_per_plot} plots on canva nbr {int(i/curve_per_plot)}')
        for j in range(i, i+curve_per_plot):
            try:
                x = df['Time'].tolist()[::skip]
                x = [float(i) for i in x]
                y = df[j].tolist()[::skip]
                y = [int(i) for i in y]
                axs[int(i/curve_per_plot)].plot(x, y, label=str(dic[j]))
            except KeyError:
                print('leaving loop, more canvas than data to plot')
                break

            axs[int(i/curve_per_plot)].legend()
            axs[int(i/curve_per_plot)].set_ylim(-1, len(ksis)+1)
            axs[int(i/curve_per_plot)].set_yticks(list(range(0, len(ksis), 5)))
            axs[int(i/curve_per_plot)].set_yticklabels(ksis[::5])
#       title=re.search('(?<=plumed_).*(?=.dat)',i).group()
#       titlel.append(title)
#       axs[n].set_title(title)
    name1 = path.split('demux/')[-1].replace('/','_')+'demuxed_plot.png'
    name2 = path.split('demux/')[-1].replace('/','_')+'ksis.txt'
    plt.savefig(path+name1)
    with open(path+name2, 'w') as writer:
        writer.write('\n'.join(ksis))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='path to directory to process', action='store', dest='f'
                        , nargs='+', required=True)
    parser.add_argument('-n', help='control tick and label sizes', action='store', dest='n'
                        , default=200)
    parser.add_argument('-cpp', help='curves per plot', action='store', default=3, dest='cpp')
    parser.add_argument('-col', help='number of column', action='store', default=2, dest='col')
    parser.add_argument('-skip', help='read every *skip* data', action='store', default=2000
                        , dest='skip')
    args = parser.parse_args()

    print(f'variables are: {args.f}, {args.n}, {args.cpp}, {args.col}, {args.skip}')
    for i in args.f:
        df, ksis = demux_parser(i)
        print(f'demux_parser returned as df:\n{df}\nand as ksis:\n{ksis}')
        plot(args.n, args.cpp, args.col, args.skip, df, ksis, i)
        print(f'plots for {i} are done.')


# If called from the command line...
if __name__ == "__main__":
    main()
