#!/usr/bin/env python
'''
Author: Jeremy Lapierre <jeremy.lapierre@uni-saarland.de>
Description:

'''
from pymol import cmd
from pymol import stored
from chempy import cpv 
import centroid

def centroid(selection, center=0, quiet=1):

    selection = 'sele'
    model = cmd.get_model(selection)
    nAtom = len(model.atom)

    centroid = cpv.get_null()

    for a in model.atom:
        centroid = cpv.add(centroid, a.coord)
    centroid = cpv.scale(centroid, 1. / nAtom)

    if not int(quiet):
        print(' centroid: [%8.3f,%8.3f,%8.3f]' % tuple(centroid))

    if int(center):
        cmd.alter_state(1, selection, "(x,y,z)=sub((x,y,z), centroid)",
                        space={'centroid': centroid, 'sub': cpv.sub})

    return centroid


def comfcen(n):
    xyz=centroid("sele")
    cmd.pseudoatom(object=n,selection='sele',pos=xyz)

cmd.extend("comfcen", comfcen)
#cmd.iterate("rpb1_2", "stored.idw.append(ID)")
#with open("id.txt", "w") as writer:
#    writer.write('[ prb1_2 ]\n')
#    for i in range(0, len(stored.idw), 15):
#        out = stored.idw[i:i+15]
#        outs = ' '.join(f'{str(x):>5}' for x in out)
#        writer.write(f'{outs:15}\n')
