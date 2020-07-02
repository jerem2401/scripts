#! /bin/bash

for i in rep*/*/; do
    echo $i
    echo 1 | gmx trjconv -f $i/traj_comp.xtc -s $i/md.tpr -o $i/nopbc1.xtc -pbc atom -ur compact -dt 4000
    wait $!
    echo 1 | gmx trjconv -f $i/nopbc1.xtc -s $i/md.tpr -o $i/nopbc2.xtc -pbc whole
    wait $!
    echo 1 | gmx trjconv -f $i/nopbc2.xtc -s $i/md.tpr -dump 0 -o $i/0.gro
done
