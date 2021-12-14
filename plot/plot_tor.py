#!/usr/bin/env python
'''plotting output of torcv tries'''

import argparse
import matplotlib.pyplot as plt
import plumed_pandas

def main():
    '''plotting output of torcv tries'''

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='give colvar_file', action='store', dest='f', type=str)
    parser.add_argument('-start', help='give start position of column', action='store',
                        dest='start', type=int)
    parser.add_argument('-end', help='give end position +1 of column', action='store',
                        dest='end', type=int)
    parser.add_argument('-col', help='list of column names to plot', nargs='*', dest='col')
    parser.add_argument('-o', help='output image name / path', action='store', dest='out',
                        type=str)
    args = parser.parse_args()

    #def skip(index):
    #    if index % 100 == 0:
    #        return True
    #    return False

    dataf = plumed_pandas.read_as_pandas(args.f, skiprows=lambda x: x % 8000)

    if args.col is not None:
        ax = dataf.plot(x='time', y=args.col)
        ax.legend([args.f])
    else:
        ax = dataf.plot(x='time', y=range(args.start, args.end))
    plt.savefig(f'{args.out}.png', format='png')
    plt.close()
    print(f'{args.out}.png has been created')

if __name__ == '__main__':
    main()
