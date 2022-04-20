#!/usr/bin/env python
'''
Author: Jeremy Lapierre <jeremy.lapierre@uni-saarland.de>
Description:

'''
import matplotlib.pyplot as plt
import matplotlib

c1 = 0
c2 = 0
avg = []
std = []

with open('avg.txt', 'r') as file:
    for line in file:
        avg.append(float(line))
        if float(line) <= 0.18:
            c1 += 1
        elif float(line) >= 0.5:
            c2 += 1
print(f"there is {c1} values below 0.18\nthere is {c2} values above 0.5")

with open('std.txt', 'r') as file:
    for line in file:
        std.append(float(line))

matplotlib.rcParams.update({'font.size': 25})

fig = plt.figure(figsize=(12, 8))
ax1 = fig.add_subplot(1, 1, 1)
ax1.errorbar(range(0,len(avg)), avg, yerr=std, marker='o', ls='none')
ax1.plot(range(0, len(avg)), [0.5]*(len(avg)))
ax1.plot(range(0, len(avg)), [0.18]*(len(avg)))

plt.savefig("bench_processed.jpeg")
