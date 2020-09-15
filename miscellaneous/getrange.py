#!/usr/bin/env python
from pymol import cmd


def getrange(sel):
    minid = min(cmd.identify(sel, 0))
    maxid = max(cmd.identify(sel, 0))
    print(str(minid)+'-'+str(maxid))

cmd.extend("getrange", getrange)
