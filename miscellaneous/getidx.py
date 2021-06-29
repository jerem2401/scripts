#!/usr/bin/env python
'''
Author: Jeremy Lapierre <jeremy.lapierre@uni-saarland.de>
Description:

'''
from pymol import cmd
from pymol import stored

cmd.load("0_chains.pdb")
cmd.select("rpb1_2", "(chain A or chain B) and 0_chains")
stored.idw = []
cmd.iterate("rpb1_2", "stored.idw.append(ID)")
with open("id.txt", "w") as writer:
    writer.write('[ prb1_2 ]\n')
    for i in range(0, len(stored.idw), 15):
        out = stored.idw[i:i+15]
        outs = ' '.join(f'{str(x):>5}' for x in out)
        writer.write(f'{outs:15}\n')
