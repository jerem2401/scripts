#!/usr/bin/env python
'''
Author(s): Jeremy Lapierre <jeremy.lapierre@uni-saarland.de>
\nDescription: This script modifies gromacs topology files to simulate with heavy hydrogens. This
allows to decrease angular and out-of-plane motion involving hydrogens and therefore allowing to
increase the time step. To conserve the total mass, the heavy atoms connected to hydrogens will
also have a proportional reduced mass. Example for methyl group with factor of 3 (default, see
argument descriptions at the end of the help) :\n
           H (1.008)                     H (3.024)
           |                             |
(1.008) H--C--H (1.008)  ==>  (3.024) H--C--H (3.024)
        (12.01)                       (5.9620)

(FLOAT): Atom mass\n

If you are simulating proteins, don't forget to constrain all bonds.\n

Some script launch examples:\n
1. heavyh_final.py : a .top file is looked in current dir and all .itp files (except *posre*.itp and
forcefield.itp) included in the .top will be processed as well (if in current directory). All
hydrogen masses are multiply by default factor (i.e. 3).
2. heavyh_final.py -p topol.top : same as 1. but .itp also looked in directory of the .top file.
3. heavyh_final.py -p topol1.itp topol2.itp: only processes the .itp files given, does not look
for .top.
4. heavyh_final.py -p topol.top -factor 4 : same as 2. but all hydrogen masses are multiply by a
factor of 4.

TODO/Notes:\n1. It is not possible to have several Hydrogens with different masses in the same topolgy
file yet. If those topologies are in different .itp files included in the .top, this will be all
fine though.\n2. If included .itp files in .top are not in the same directory as the .top,
the script will not process it.

'''

import argparse
from glob import glob
import re
from os import path
from shutil import move


def file_checker(top, itp, suff):
    """This function checks for .itp (not containing 'posre' or '/') included in .top and try to
    look for them in the current directory or in the directory name of the .top file given by -p.\n
    Check if .itp/.top files have an [ atoms ] block to modifiy, else remove it from input files.\n
    Finally, builds output names for the topologies"""


    if top != []:
        print(f'Working on following .top file: {top}\n')
        regxitp = re.compile(r'(?<=#include ")(?!.*(\b\/\b|posre)).*(?=")')
        with open(top[0], 'r') as intop:
            for line in intop:
                if regxitp.search(line):
                    if regxitp.search(line).group() not in itp:
                        itp.append(regxitp.search(line).group())

        #Check if those .itp are in current directory or in the dirname of .top given by -p
        for i in itp:
            if not path.isfile(i):
                dirname_of_top = path.dirname(top[0])
                print(f"{i} included in .top not found in './', "
                      f"let's look in dirname of the -p argument (i.e. {dirname_of_top})")
                guess_itp_path = f'{dirname_of_top}/{i}'
                if path.isfile(guess_itp_path):
                    print(f'{i} found in {dirname_of_top}')
                    itp[itp.index(i)] = guess_itp_path
                else:
                    print(f"{i} included in .top also not present in {dirname_of_top}:"
                          f" exiting script")
                    raise SystemExit()

    topols = top + itp

    regx_atmblock = re.compile(r'^\[ atoms \]\n')
    topol2rm = []

    for i in topols:
        keep = 0
        with open(i, 'r') as topfile:
            for line in topfile:
                if regx_atmblock.search(line):
                    keep = 1
        if keep == 0:
            print(f'\nNo [ atoms ] block in {i}, this file will be ignored, make sure to have all'
                  f' itp and top in the same directory in case full path of itp is not in'
                  f' the top file though.')
            topol2rm.append(i)

    topols = [topfile for topfile in topols if topfile not in topol2rm]

    print(f'\nThe following .itp and/or .top will be used for generating heavy H topologies:'
          f'\n{topols}\n')

    topolout = [re.sub(r'(\.itp|\.top)', '', path.basename(name))
                +suff
                +'.'
                +name.split('.')[-1] for name in topols]

    topolin_topolout = {topols[i]: topolout[i] for i in range(len(topols))}

    return topolin_topolout

def top_parser(topol):
    """This function is parsing the topol file: 1. Pass lines until the [ atoms ] block is reached,
    2. Check format of [ atoms ] block,
    3.Read the atom block (but skip commments) until it is finished to create hydrogen and
    heavy atom sets (based on hydrogen_mass), 4. Read the bond block and make a list of col 1 and 2.
    Also get the hydrogen mass of the topology."""

    hydrogen_mass = set()
    hydrogenid = set()
    heavyatmid = set()
    bondlist = []

    #Sanity check variable
    massin = 0
    error = []

    with open(topol, 'r') as top:
        for atom_line in iter(top.readline, '[ atoms ]\n'):
            pass
        for index, atom_line in enumerate(iter(top.readline, '\n')):
            atom_line_list = atom_line.split()
            if index == 0:
                if 'mass' not in atom_line_list:
                    error.append('no_mass_error') #'No mass in [ atoms ] block'
                    #raise SystemExit(error)
                if 'mass' not in atom_line_list[8] and 'charge' not in atom_line_list[7]:
                    error.append('atm_format_error') #'unusual [ atoms ] block format input, improve the code.'
                    #raise SystemExit(error)
                string_untilmass = re.search(rf".*(?=\bcharge\b)", atom_line).group()
                nbrchar_b4_mass = len(string_untilmass) + 6 #6 is len(charge)
                string_withmass = re.search(rf".*(?=\bmass\b)", atom_line).group()
                nbrchar4mass = len(string_withmass) + 4 - nbrchar_b4_mass #4 is len(mass)
            elif atom_line.startswith(';'):
                pass
            elif 0.99 <= float(atom_line_list[7]) <= 1.1:
                hydrogen_mass.add(atom_line_list[7])
                hydrogenid.add(atom_line_list[0])
                massin += float(atom_line_list[7])
            else:
                heavyatmid.add(atom_line_list[0])
                massin += float((atom_line_list[7]))
        iterbond = iter(top.readline, '\n')
        next(iterbond)
        next(iterbond)
        for atom_line in iterbond:
            bondlist.append(atom_line.split()[0:2])

    if len(hydrogen_mass) > 1:
        #error = ('Not possible to have different hydrogen masses in the same topology, see --help')
        error.append('multi_h')
        #raise SystemExit(error)

    try:
        hydrogen_mass = list(hydrogen_mass)[0]
    except IndexError:
        error.append('no_hydrogen')

    return (nbrchar_b4_mass, nbrchar4mass, float(hydrogen_mass),
            hydrogenid, heavyatmid, bondlist, massin, error)


def h_heavyatm_matcher(bondlist, hydrogenid):
    """This function creates a dictionary:\n
       { heavy atom id being connected to hydrogens : nbr of associated H }."""

    print('ok')
    h_connected_heavy = {}
    #Sanity check variable
    hbond_nbr = 0

    for col1col2 in bondlist:
        if col1col2[1] in hydrogenid:
            h_connected_heavy[col1col2[0]] = h_connected_heavy.get(col1col2[0], 0) +1
            hbond_nbr += 1
        elif col1col2[0] in hydrogenid:
            h_connected_heavy[col1col2[1]] = h_connected_heavy.get(col1col2[1], 0) +1
            hbond_nbr += 1

    if hbond_nbr != len(hydrogenid):
        error = ('Error: the number of parsed hydrogens in hydrogenid variable is not the same as '
                 'the number of hydrogen bonds in hbond_nbr')
        raise SystemExit(error)

    print('ok')
    return h_connected_heavy


def writer(hydrogen_mass, massfactor, nbrchar_b4_mass, nbrchar4mass, topol, out, h_connected_heavy):
    """This function writes the output topology (and/or .itps)"""

    new_hydrogen_mass = hydrogen_mass*massfactor
    mass2remove = new_hydrogen_mass-hydrogen_mass
    new_hydrogen_mass = f'{new_hydrogen_mass:.4f}'
    new_hydrogen_mass = f'{new_hydrogen_mass:>11}'
    #Sanity check variable
    massout = 0
    regx4mass = re.compile(r"(?<=^.{%s}).{%s}" % (nbrchar_b4_mass, nbrchar4mass))

    with open(topol, 'r') as intop, open(out, 'w') as outtop:
        for atom_line in iter(intop.readline, '[ atoms ]\n'):
            outtop.write(atom_line)
        outtop.write(f'[ atoms ]\n')
        for atom_line in iter(intop.readline, '\n'):
            atom_line_list = atom_line.split()
            if atom_line.startswith(';'):
                outtop.write(atom_line)
            elif 0.99 <= float(atom_line_list[7]) <= 1.1:
                atom_line = regx4mass.sub(new_hydrogen_mass, atom_line)
                outtop.write(atom_line)
                massout += float(new_hydrogen_mass)
            elif atom_line_list[0] in h_connected_heavy.keys():
                new_heavyatm_mass = float(atom_line_list[7])                            \
                                    - (mass2remove*h_connected_heavy[atom_line_list[0]])
                new_heavyatm_mass = f'{new_heavyatm_mass:.4f}'
                new_heavyatm_mass = f'{new_heavyatm_mass:>11}'
                atom_line = regx4mass.sub(new_heavyatm_mass, atom_line)
                outtop.write(atom_line)
                massout += float(new_heavyatm_mass)
            else:
                massout += float(atom_line_list[7])
                outtop.write(atom_line)
        outtop.write('\n')
        for atom_line in iter(intop.readline, ''):
            outtop.write(atom_line)

    return massout


def main():
    """This is just the main function of the script executing parser, h_heavyatm_matcher
    and writer functions, see __doc__ for details"""

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-p',
                        default=glob('./*.top'),
                        help='Give input topology (.top and/or itp), default: ./*.top, '
                             'if several .top in current dir or given, error will be raised',
                        dest='topol',
                        nargs='*')
    parser.add_argument('-factor',
                        default=3,
                        help='Give the factor with which each hydrogen mass should be multiplied by'
                             ' (default: %(default)s)',
                        dest='massfactor',
                        type=int)
    parser.add_argument('-o',
                        default='_heavyH',
                        help=f'Suffix for output topology file name (default: %(default)s)',
                        action='store',
                        dest='suff',
                        type=str)
    args = parser.parse_args()

    top = [top for top in args.topol if ".top" in top]
    itp = [itp for itp in args.topol if ".itp" in itp]

    #Manage only one .top/inputs, this simplify the script. If several .top just loop over the
    if len(top) > 1:
        error = ('Several .top files in the current directory or given, please use -p to precise'
                 ' which .top should be used (you can however give as many .itp as you want or even'
                 ' just one .itp as input if you want). Tip: if you want to process several .top,'
                 ' just use a bash for loop over your .top using this scripti.')
        raise SystemExit(error)

    #see file_checker function
    topolin_topolout = file_checker(top, itp, args.suff)

    #see top_parser, h_heavyatm_matcher and writer function
    #variable needed to not modifiy the name of .itp without H in the .top
    no_hydrogen = []
    for items in topolin_topolout.items():
        print(f"processing {items[0]}...\n")

        try:
            
            nbrchar_b4_mass, nbrchar4mass, hydrogen_mass, hydrogenid, heavyatmid, bondlist, massin = top_parser(items[0])
            h_connected_heavy = h_heavyatm_matcher(bondlist, hydrogenid)
            massout = writer(hydrogen_mass, args.massfactor, nbrchar_b4_mass, nbrchar4mass,
                             items[0], items[1], h_connected_heavy)
            print(f'Difference between overall mass in input '
                  f'and overall mass in output is: {massin-massout} (should be 0)'
                  f'\n(massin is: {massin}, massout is: {massout})')
        #catch IndexEerror raised by top_parser because hydrogen_mass is empty list
        except IndexError:
            print(f'***CAREFUL*** No hydrogens were found in: {items[0]}. If you think it makes '
                  f'sens for this molecular entity to not have hydrogens, first make sure that (in '
                  f'case of .itp file) this file is copied in the working directory and then you'
                  f' can ignore what follows.\n'
                  f'Else possible reasons could be:\n1. The hydrogen masses in the input are not '
                  f'between 0.99 and 1.1, this means that the script needs to be optimized.'
                  f'\n2. Because the top formatting is unusual'
                  f', you should review the code in this case too.\n')
            no_hydrogen.append(items[0])
    #This replaces the name of the modified .itp in the .top file if .itp contains H
    if top != []:
        #case where .top does not contain H
        if top[0] not in topolin_topolout:
            new_name = top[0].replace('.top', args.suff+'.top')
            with open(top[0], 'r') as intop, open(new_name, 'w') as outtop:
                for line in intop:
                    for i in topolin_topolout:
                        if i not in no_hydrogen and i in line:
                            line = line.replace(i, topolin_topolout[i])
                    outtop.write(line)
        #case where .top does contain H
        else:
            with open(topolin_topolout[top[0]], 'r') as outtop, open('temp.top', 'w') as temp:
                for line in outtop:
                    for i in topolin_topolout:
                        if i not in no_hydrogen and i in line:
                            line = line.replace(i, topolin_topolout[i])
                    temp.write(line)
            #print('moving ./temp.top to', "./"+topolin_topolout[top[0]])
            move('./temp.top', "./"+topolin_topolout[top[0]])


# If called from the command line...
if __name__ == "__main__":
    main()
