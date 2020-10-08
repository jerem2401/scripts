#!/usr/bin/env python
'''
Author(s): Jeremy Lapierre <jeremy.lapierre@uni-saarland.de>
\nDescription: This script modify gromacs topology files to simulate with heavy hydrogens. This
allows to decrease angular and out-of-plane motion involving hydrogens and therefore allowing to
increase the time step. To conserve the total mass, the heavy atoms conntected to hydrogens will
also have a proportional reduced mass. Example for methyl group with factor of 3 (default, see
argument descriptions at the end of the help) :\n
           H (1.008)                     H (3.024)
           |                             |
(1.008) H--C--H (1.008)  ==>  (3.024) H--C--H (3.024)
        (12.01)                       (5.9620)

(FLOAT): Atom mass\n

'''

import argparse
from glob import glob
import re
from os import path

def top_parser(topol):
    """This function is parsing the topol file: 1. Pass until the [ atoms ] block is reached,
    2.Read the atom block (but skip commments) until it is finished to create hydrogen and
     heavy atom sets (based on H mass), 3. Read the bond block and make a list of col 1 and 2"""

    hydrogenid = set()
    heavyatmid = set()
    bondlist = []

    #Sanity check variable
    massin = 0

    with open(topol, 'r') as top:
        for atom_line in iter(top.readline, '[ atoms ]\n'):
            pass
        for atom_line in iter(top.readline, '\n'):
            if atom_line.startswith(';'):
                continue
            atom_line_list = atom_line.split()
            if 0.99 <= float(atom_line_list[7]) <= 1.1:
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

    return hydrogenid, heavyatmid, bondlist, massin


def h_heavyatm_matcher(bondlist, hydrogenid):
    """This function creates a dictionary:\n
       { heavy atom id being connected to hydrogens : nbr of associated H }."""

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

    return h_connected_heavy


def writer(massfactor, topol, out, h_connected_heavy):
    """This function writes the output topology (and/or .itps)"""

    new_hydrogen_mass = 1.008*massfactor
    mass2remove = new_hydrogen_mass-1.008
    new_hydrogen_mass = f'{new_hydrogen_mass:.4f}'
    new_hydrogen_mass = f'{new_hydrogen_mass:>11}'
    #Sanity check variable
    massout = 0
    regx4mass = re.compile(r"(?<=^.{56}).{11}")

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
                #heavy_h = atom_line.replace(" 1.008 ", f" {new_hydrogen_mass} ")
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
                        help='Give input topology (default: ./*.top, '
                             'if several .top in current dir, error will be raised',
                        dest='topol',
                        nargs='*')
    parser.add_argument('-factor',
                        default=3,
                        help='Give the factor with which each hydrogen mass should be multiply to'
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

    if len(args.topol) > 1:
        error = 'Several .top files in the current directory, please use -p to select one'
        raise SystemExit(error)

    print(f'Working on following .top file: {args.topol[0]}\n')

    #making list of .itp in .top (but which does not contain 'posre' or '/' in name)
    topol = [args.topol[0]]
    regxitp = re.compile(r'(?<=#include ")(?!.*(\b\/\b|posre)).*(?=")')

    with open(args.topol[0], 'r') as intop:
        for line in intop:
            if regxitp.search(line):
                topol.append(regxitp.search(line).group())

    #Check if those .itp are in current directory or in the dirname of .top given by -p
    for i in topol:
        if not path.isfile(i):
            dirname_of_top = path.dirname(args.topol[0])
            print(f"{i} included in .top not found in './', "
                  f"let's look in dirname of the -p argument (i.e. {dirname_of_top})")
            guess_itp_path = f'{dirname_of_top}/{i}'
            if path.isfile(guess_itp_path):
                print(f'{i} found in {dirname_of_top}')
                topol[topol.index(i)] = guess_itp_path
            else:
                print(f"{i} included in .top also not present in {dirname_of_top}: exiting script")
                raise SystemExit()

    #checking if .itp files have an [ atoms ] block to modifiy, else remove it from input files
    regx_atmblock = re.compile(r'^\[ atoms \]\n')
    topol2rm = []

    for i in topol:
        keep = 0
        with open(i, 'r') as top:
            for line in top:
                if regx_atmblock.search(line):
                    keep = 1
        if keep == 0:
            print(f'\nNo [ atoms ] block in {i}, this file will be ignored, make sure to have all'
                  f' itp and top in the same directory in case full path of itp is not in'
                  f' the top file though.')
            topol2rm.append(i)

    topol = [top for top in topol if top not in topol2rm]

    print(f'\nThe following .itp and/or .top will be used for generating heavy H topologies:'
          f'\n{topol}\n')

    topolout = [re.sub(r'(\.itp|\.top)', '', path.basename(name))
                +args.suff
                +'.'
                +name.split('.')[-1] for name in topol]

    for topin, topout in zip(topol, topolout):
        print(f"processing {topin}...\n")

        hydrogenid, heavyatmid, bondlist, massin = top_parser(topin)
        if len(hydrogenid) != 0:
            h_connected_heavy = h_heavyatm_matcher(bondlist, hydrogenid)
            massout = writer(args.massfactor, topin, topout, h_connected_heavy)
            print(f'Difference between overall mass in input '
                  f'and overall mass in output is: {massin-massout} (should be 0)'
                  f'\n(massin is: {massin}, massout is: {massout})')
        else:
            print(f'***CAREFUL*** No hydrogens were found in: {topin}. If you think it makes sens '
                  f'for this molecular entity to not have hydrogens, first make sure that (in case '
                  f'of .itp file) this file is copied in the working directory and then you can'
                  f' ignore what follows.\n'
                  f'Else possible reasons could be:\n1. The hydrogen masses in the input are not '
                  f'between 0.99 and 1.1, this means that the script needs to be optimized.'
                  f'\n2. Because the top formatting is unusual (probably mass not in '
                  f'8th column), you should review the code in this case too.\n')


# If called from the command line...
if __name__ == "__main__":
    main()
