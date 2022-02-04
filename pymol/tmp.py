#!/usr/bin/env python
from pymol import cmd
from chempy import cpv


def centroid(selection='sele', center=0, quiet=1):

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


def tmp(n):
    xyz = centroid('sele')
    print(xyz)
    cmd.pseudoatom(object=n, pos=xyz)


cmd.extend("tmp", tmp)
