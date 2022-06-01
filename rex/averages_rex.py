#!/usr/bin/env python
'''
Author: Jeremy Lapierre <jeremy.lapierre@uni-saarland.de>
Description:

'''
import os
import pandas as pd
import glob
from itertools import dropwhile, islice, takewhile, zip_longest
import numpy as np

e = os.popen('export LC_ALL=C; command ls E* -d | sort -t _ -k 2 -n').read().split()
tormp = os.popen('cat latepos.txt 2>/dev/null').read()
torm = os.popen('cat late.txt 2>/dev/null').read()

if len(torm) > 0:
    e = e[0:len(e)-int(torm)]
elif len(tormp) > 0:
    e = e[0:len(e)-int(tormp)]

e = [float(x[2::]) for x in e]

files = glob.glob('./log*.txt')
data = []

for i in files:
    tmp = []
    with open(i, 'r') as reader:
        start = dropwhile(lambda x: 'Repl  average probabilities' not in x, reader)
        end = takewhile(lambda x: 'Repl  number of exchanges:' not in x, start)
        for line in islice(end, 1, None):
            tmp.append(line)
    tmp = tmp[1].split()[1::]
    data.append([float(i) for i in tmp])

if len(torm) > 0:
    data = [i[0:len(i)-int(torm)] for i in data]
print(data)

arr = np.vstack(data)
avg = np.around(np.average(arr, axis=0), decimals=5)
std = np.around(np.std(arr, axis=0, ddof=1), decimals=5)

avg2 = np.average(avg, axis=0)
std2 = np.std(avg, axis=0, ddof=1)
print(f"avg:\n{avg}\nstd:\n{std}\navg2:\n{avg2}\nstde:\n{std2}\n")


p = avg
d = []
prev = e[0]
for i in e[1::]:
    d.append(round(i-prev, 5))
    prev = i

scale_key = [(e[i], e[i+1]) for i in range(0, len(e)-1)]
scale_func = [(round((1-p[i])/d[i], 5)) for i in range(0, len(p))]
scale = pd.DataFrame(list(zip(scale_key, scale_func, p, [1-x for x in p], d)),
                     columns=['scale_interval', 'scale_funcparam', 'proba_acc', 'proba_rejec', 'distance'])


############

ksis = [e[0]]      #IMPORTANT
tp_reject = 0.64   #
step = 0.00001     #
ruler = []         #
pruler = []        #

for j in range(0, len(scale)):
    inter = list(np.arange(scale.iloc[j, 0][0], scale.iloc[j, 0][1], step))[0:-1]
    ruler += [round(x, 5) for x in inter]
    pruler += [(scale.iloc[j, 3]/len(inter))]*len(inter)
ruler += [e[-1]]

benchs = glob.glob('./bench*/')
if benchs == []:
    maxne = e[-1]
    prev_pruler = pruler
else:
    benchs = sorted(benchs)
    maxne = os.popen('export LC_ALL=C; command ls bench1/E* -d | sort -t _ -k 2 -n | tail -n1').read()
    maxne = float(maxne[9::])
    with open(benchs[-1]+'rex_opti.out', 'r') as reader:
        start = dropwhile(lambda x: 'pruler' not in x, reader)
        end   = takewhile(lambda x: 'results:' not in x, start)
        for line in islice(end, 1, None):
            prev_pruler = line[1:-2]
    prev_pruler = [float(x) for x in prev_pruler.split(', ')]
    

print(f'previous pruler: {prev_pruler[0]}, {prev_pruler[-1]}, {len(prev_pruler)}')
print(f'curent pruler: {pruler[0]}, {pruler[-1]}, {len(pruler)}')

concat = np.array(list(zip_longest(*[prev_pruler, pruler])),dtype=float).transpose()
mean_pruler = np.nanmean(concat, axis=0)
print(f'mean pruler: {mean_pruler[0]}, {mean_pruler[-1]}, {len(mean_pruler)}')
fmean_pruler = mean_pruler[0:len(prev_pruler)]
print(f'final pruler: {fmean_pruler[0]}, {fmean_pruler[-1]}, {len(fmean_pruler)}')

nbrruler2add = len(fmean_pruler)-len(pruler)
if nbrruler2add != 0:
    ruler = list(np.arange(e[0], maxne, step))
    ruler = [round(x, 5) for x in ruler]

j = 0
while (j<= len(scale)-1) and (tp_reject+0.1 >= scale.iloc[j, 3] >= tp_reject-0.10):
    ksis += [scale.iloc[j, 0][1]]
    j += 1
tmp = []
[tmp.append(x) for x in ksis if x not in tmp]
# remove last ksi because this was the one responsible for getting out of the while loop
if len(tmp) != 1:
    tmp = tmp[0:-1]
    ksis = tmp
    index = ruler.index(ksis[-1])
else:
    index = 1

p = 0
for i in range(index-1,len(fmean_pruler)):
    p += fmean_pruler[i]
    if p >= tp_reject:
        ksis += [ruler[i+1]]
        p = 0
    else:
        continue

p = 0
while ksis[-1] < ruler[-1]:
    p += fmean_pruler[-1]
    if p >= tp_reject:
        ksis += [round(ksis[-1]+(step*(p/fmean_pruler[-1])), 5)]
        p = 0

print('\n')
print(p, scale_key, scale_func, scale, sep='\n\n')
print(ruler[0], ruler[-1], len(ruler), pruler[0], pruler[-1], len(pruler))
print('\n')
print(f'pruler:\n{pruler}')
print('results:')
print(*(f'tmp_{x:.5f}' for x in ksis), sep=' ')
