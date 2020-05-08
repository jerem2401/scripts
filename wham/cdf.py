#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('-hist', help='taking counts and x from histo files', nargs='*', action='store', dest='hist')
parser.add_argument('-plot', help='plot the cdf', action='store_true', dest='plot')
parser.add_argument('-val', help='give y (could be many separated by " "), ruturn x', nargs='*', action='store', dest='val', type=float)
args = parser.parse_args()

if args.hist != None:
    #files=os.popen('ls '+args.hist+'*').read().split()
    v_given = [args.val[0]]*len(args.hist)
    v_closest = []
    torm = []

    for i in args.hist:
        f = open(i, "r")
        lines = f.readlines()
        x = []
        y = []
        for j in lines:
            x.append(float(j.split(' ')[0]))
            y.append(float(j.split(' ')[1]))
        f.close()

        cdf = np.cumsum(y)
        dic = {k:v for k, v in zip(cdf, x)}
        cdfa = np.asarray(cdf)
        idx = (np.abs(cdfa - args.val)).argmin()

        print('value given: '+str(args.val)+'\nclosest value found: '+str(cdf[idx])+'\nmust remove everything > '+str(dic[cdf[idx]])+' in file: '+str(i))
        v_closest.append(float(cdf[idx]))
        torm.append(float(dic[cdf[idx]]))

        if args.plot:
            fig, ax = plt.subplots(figsize=(10, 7))
            ax.plot(x, cdf)

            plt.title(i.strip('.txt'))
            plt.savefig('./'+i.split('/')[-1].strip('.txt')+'.jpeg')
print(v_given[0], v_closest[0], args.hist[0], torm[0], sep='    ')

f = open('./hist_curing.txt', 'w')
#f.write('#v_given v_closest file 2rm')
for i in range(len(args.hist)):
    f.write("%f %f %s %f\n" % (v_given[i], v_closest[i], args.hist[i], torm[i]))
f.close()

#np.savetxt(r'./hist_curing.txt', np.transpose([v_given,v_closest,args.hist,torm]), header='#v_given v_closest file 2rm', fmt="%f %f %s %f")
#else:
#not supported yet

#num_bins = 20
#counts, bin_edges = np.histogram (data, bins=num_bins, normed=True)
#cdf = np.cumsum (counts)
#plt.plot (bin_edges[1:], cdf/cdf[-1])
