#!/usr/bin/env python

#for i in ../../E_*/colv*; do var2=$(basename "$i"); awk '(NR>5) {print $1,$7,$6}' $i > ./$var2; done
#for i in ./colv*; do cnt=$(echo "$i" | grep -oP '(?<=colvar).*(?=_)'); cnt2=$(echo "$i" | grep -oP '(?<=_).*(?=.txt)'); echo "$i   $cnt   $cnt2   5000   10000" >> metd.txt; done
#column -t metd.txt > temp.txt
#sed '/inf/d' 2dpmf.txt > 2dpmf_clean.txt

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='give 2dpmf_clean.txt file', action='store', dest='f', type=str)
    parser.add_argument('-o', help='give path for plot output', action='store', dest='o')
    args = parser.parse_args()

    #data
    x, y, free, prob = np.loadtxt(args.f, unpack=True)
    
    #fig
    fig = plt.figure(figsize=(20,20))
    
    cmap = plt.cm.get_cmap('rainbow')
    norm = matplotlib.colors.Normalize(vmin=min(free), vmax=max(free))
    colors = [cmap(norm(value)) for value in free]

    ax1 = fig.add_subplot(1, 1, 1)
    plot = ax1.scatter(x, y, color=colors, marker="s")

    cax, _ = matplotlib.colorbar.make_axes(ax1, shrink=0.5, anchor=(0.9,1.2))
    cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)
 
    plt.tight_layout()
    plt.savefig(str(args.o)+'/2dpmf.jpeg', quality=60)

# If called from the command line... 
if __name__ == "__main__":
    main()

