#!/usr/bin/env python
'''
Author: Jeremy Lapierre <jeremy.lapierre@uni-saarland.de>
Description:

'''

import pymol
from pymol import cmd
from pymol import stored

import __main__
__main__.pymol_argv = ['pymol', '-qc']  # Quiet and no GUI


pymol.finish_launching()

cmd.load("1_ref_structures.pse")


#stored.RPB1=[]
#stored.RPB2=[]
#stored.RPB3=[]
#stored.RPB4=[]
#stored.RPB5=[]
#stored.RPB6=[]
#stored.RPB7=[]
#stored.RPB8=[]
#stored.RPB9=[]
#stored.RPB10=[]
#stored.RPB11=[]
#stored.RPB12=[]
#stored.TFIIB=[]
#stored.TFIIAa=[]
#stored.TFIIAb=[]
#stored.TATAbbp=[]
#stored.TFIIEa=[]
#stored.TFIIEb=[]
#stored.TFIIFa=[]
#stored.TFIIFb=[]
#stored.TFIIS=[]
#stored.DNAnt=[]
#stored.DNAt=[]

print("dic = {")
for i in (cmd.get_object_list('(all)')):
    if "_cc" in i:
        print(f'"{i}": {cmd.identify(i)},')
print("}")

pymol.cmd.quit()
