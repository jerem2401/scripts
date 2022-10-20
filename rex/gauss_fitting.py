#!/usr/bin/env python
'''
Author: Jeremy Lapierre <jeremy.lapierre@uni-saarland.de>
Description:

'''
import numpy as np


def func(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))


def make_histo(data, process_col, binning, histoname):
    import pandas as pd
    df = pd.read_csv(data, delimiter=r"\s+", header=None, dtype='float')
    df = df.iloc[:, 1]
    hist2, bin_edges = np.histogram(df, bins=[i for i in np.around(np.arange(min(df), max(df),
                                    binning), decimals=3)], density=True)
    hist2 = hist2 / hist2.sum()
    d = {'z': bin_edges[0:-1], 'hist': hist2}
    hdf = pd.DataFrame(d)
    np.savetxt(histoname, hdf.values, fmt='%.6f')
    print(histoname+' is done')
    return hdf


def plot_lbd(files):
    import matplotlib.pyplot as plt
    import matplotlib
    import matplotlib.font_manager
    import matplotlib as mpl
    from scipy.optimize import curve_fit
    from matplotlib.pyplot import cm
    import re

    dpi = 300
    n = 4

    mpl.rcParams.update(plt.rcParamsDefault)
    mpl.rcParams['figure.dpi'] = dpi  # Frame width
    mpl.rcParams['axes.linewidth'] = 0.1*n  # edge line width
    mpl.rcParams['axes.grid.which'] = 'both'
    mpl.rcParams['grid.linewidth'] = 0.05*n
    mpl.rcParams['grid.alpha'] = 1
    mpl.rcParams['lines.linewidth'] = 0.15*n
    mpl.rcParams['lines.markersize'] = 0.5*n
    mpl.rcParams['lines.markeredgewidth'] = 0*n
    # major tick size in points
    mpl.rcParams['xtick.major.size'] = mpl.rcParams['ytick.major.size'] = 0.6*n
    # major tick width in points
    mpl.rcParams['xtick.major.width'] = mpl.rcParams['ytick.major.width'] = 0.1*n
    mpl.rcParams['axes.labelpad'] = 1.4*n
    mpl.rcParams['xtick.major.pad'] = mpl.rcParams['ytick.major.pad'] = 0.8*n
    plt.rcParams.update({'font.sans-serif': 'Helvetica'})
    mpl.rcParams['font.size'] = 1.7*n
    mpl.rcParams['legend.edgecolor'] = 'white'
    mpl.rcParams['legend.framealpha'] = 1
    mpl.rcParams['legend.handlelength'] = n*0.2
    mpl.rcParams['legend.fontsize'] = n*1.1

    centi = 1/2.54  # centimeters in inches

    fig, axs = plt.subplots(1, 1, figsize=(1.75*centi*n, 1.75*centi*n))
    color = iter(cm.rainbow(np.linspace(0, 1, 16)))

    for n, i in enumerate(files):
        hdf = make_histo(i, ['time', 'energy'], 200, i.replace('.xvg', '_histo.xvg'))

        x = np.asarray(hdf.iloc[:, 0])
        y = np.asarray(hdf.iloc[:, 1])

        popt, pcov = curve_fit(func, x, y, [max(y), np.mean(x), np.std(x)])
        ym = func(x, popt[0], popt[1], popt[2])

        m = re.search('(?<=ener_).*(?=.xvg)', i)
        label = 1-float(m.group(0))
        c = next(color)
        axs.plot(x, y, '.', c=c)
        axs.plot(x, ym, '-', c=c, label=str(label))

    plt.legend()
    m = re.search('.*(?=/ener_*)', i)
    plt.savefig(m.group(0)+'/potential_lbd.jpeg')


def main():
    import sys
    import os
    print(sys.argv[1:])
    for i in sys.argv[1:]:
        files = os.popen(f"ls {i}/energy/*xvg | sort -t _ -k 3 -n").read().split()
        plot_lbd(files)


if __name__ == "__main__":
    main()
