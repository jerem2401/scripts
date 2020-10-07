#!/usr/bin/env python
'''
Author(s): Jeremy Lapierre <jeremy.lapierre@uni-saarland.de>
\nDescription:\n
           H (1.008)                     H (3.024)
           |                             |
(1.008) H--C--H (1.008)  ==>  (3.024) H--C--H (3.024)
        (12.01)                       (5.9620)

(FLOAT): Atom mass
'''

import argparse

def main():
    """This is just the main function of the script, see __doc__ for details"""

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-p',
                        help='Give input topology',
                        dest='topol',
                        required=True,
                        type=str)
    parser.add_argument('-factor',
                        default=3,
                        help='Give the factor with which each hydrogen mass should be multiply to'
                             ' (default: %(default)s)',
                        dest='massfactor',
                        type=int)
    parser.add_argument('-o',
                        default='topol_heavyH.top',
                        help=f'Output topology file name (default: %(default)s)',
                        action='store',
                        dest='out',
                        type=str)
    args = parser.parse_args()

    ####Variable initializations
    new_hydrogen_mass = 1.008*args.massfactor
    mass2remove = new_hydrogen_mass-1.008

    hydrogenid = set()
    heavyatmid = set()
    bondlist = []

    #Sanity check variables
    massin = 0
    massout = 0
    hbond_nbr = 0

    ####Parsing
    with open(args.topol, 'r') as top:
        # Read line by line and do nothing until the line '[ atoms ]\n' is reached
        for atom_line in iter(top.readline, '[ atoms ]\n'):
            pass
        # THEN, read line by line until the [ atoms ] block is finished...
        for atom_line in iter(top.readline, '\n'):
            # ...but skip lines starting with ";"
            if atom_line.startswith(';'):
                continue
            # Split the line in a list
            atom_line_list = atom_line.split()
            # If mass of atom is 1.008, it's an H, so add the id to hydrogenid
            if atom_line_list[7] == "1.008":
                hydrogenid.add(atom_line_list[0])
                massin += float(atom_line_list[7])
            # Else it's an heavy atom, so add the id to heavyatmid
            else:
                heavyatmid.add(atom_line_list[0])
                massin += float((atom_line_list[7]))
        # THEN, read line by line until the [ bonds ] block is finished but skip the first 2 lines
        iterbond = iter(top.readline, '\n')
        next(iterbond)
        next(iterbond)
        for atom_line in iterbond:
            # Put [col1,col2] of each line of bond sections into bondlist
            bondlist.append(atom_line.split()[0:2])


    #Create a dic:  { heavy atom id being connected to hydrogens : nbr of associated H }
    h_connected_heavy = {}
    for col1col2 in bondlist:
        #if col 2 of bonds section is in hydrogenid (assuming that H can only be in 2nd col)...
        if col1col2[1] in hydrogenid:
            #...add to dic the heavy atom id (from 1st col of bonds) and associate +1 hydrogen
            h_connected_heavy[col1col2[0]] = h_connected_heavy.get(col1col2[0], 0) +1
            hbond_nbr += 1
        elif col1col2[0] in hydrogenid:
            h_connected_heavy[col1col2[1]] = h_connected_heavy.get(col1col2[1], 0) +1
            hbond_nbr += 1
    #print(f"h_connected_heavy is: {h_connected_heavy}")

    if hbond_nbr != len(hydrogenid):
        print('Error: the number of parsed hydrogens in hydrogenid variable is not the same as '
              'the number of hydrogen bonds in hbond_nbr')

    with open(args.topol, 'r') as intop, open(args.out, 'w') as outtop:
        for atom_line in iter(intop.readline, '[ atoms ]\n'):
            outtop.write(atom_line)
        outtop.write(f'[ atoms ]\n')
        for atom_line in iter(intop.readline, '\n'):
            atom_line_list = atom_line.split()
            print(atom_line_list)
            if atom_line_list[7] == "1.008":
                massout += float(atom_line_list[7])
                heavy_h = atom_line.replace(" 1.008 ", f" {new_hydrogen_mass} ")
                outtop.write(heavy_h)
                continue
            if atom_line_list[0] in h_connected_heavy.keys():
                massout += float(atom_line_list[7])
                new_heavyatm_mass = float(atom_line_list[7])                                   \
                                    - (mass2remove*h_connected_heavy[atom_line_list[0]])
                new_heavyatm_mass = f'{new_heavyatm_mass:.4f}'
                outtop.write(f'{atom_line_list[0]:>6}{atom_line_list[1]:>11}{atom_line_list[2]:>7}'
                             f'{atom_line_list[3]:>7}{atom_line_list[4]:>7}{atom_line_list[5]:>7}'
                             f'{atom_line_list[6]:>11}{new_heavyatm_mass:>11}{atom_line_list[8]:>4}'
                             f'{atom_line_list[9]:>5} {atom_line_list[10]}\n')
                continue
            if atom_line.startswith(';'):
                outtop.write(atom_line)
            else:
                massout += float(atom_line_list[7])
                outtop.write(atom_line)
        outtop.write('\n')
        for atom_line in iter(intop.readline, ''):
            outtop.write(atom_line)
#    print(f'ninth is: {atom_line_list[9]}\n10th is: {atom_line_list[10]}')
#    print(f"len is: {len(hydrogenid)}, c1 is: {c1}")

    print(f'Difference between overall mass in input '
          f'and overall mass in output is: {massin-massout} (should be 0)'
          f'\n(massin is: {massin}, massout is: {massout})')
# If called from the command line...
if __name__ == "__main__":
    main()
