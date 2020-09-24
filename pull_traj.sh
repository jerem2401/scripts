#! /bin/bash

#for i in rep*/*/; do
#    echo $i
#    echo 1 | gmx trjconv -f $i/traj_comp.xtc -s $i/md.tpr -o $i/nopbc1.xtc -pbc atom -ur compact -dt 600
#    wait $!
#    echo 1 | gmx trjconv -f $i/nopbc1.xtc -s $i/md.tpr -o $i/nopbc2.xtc -pbc whole
#    wait $!
#    echo 1 | gmx trjconv -f $i/nopbc2.xtc -s $i/md.tpr -dump 0 -o $i/0.gro
#done
echo $1
echo 24 | gmx trjconv -f $1/traj_comp.xtc -s $1/md.tpr -o $1/nopbc1.xtc -ur compact -pbc atom -n /data/users/jeremy/simulation/syncsim/pol/meta/index.ndx -skip 2
wait
echo 24 | gmx trjconv -f $1/nopbc1.xtc -s $1/md.tpr -o $1/nopbc2.xtc -ur compact -pbc whole -n /data/users/jeremy/simulation/syncsim/pol/meta/index.ndx
wait
echo 24 | gmx trjconv -f $1/nopbc2.xtc -s $1/md.tpr -o $1/0.pdb -dump 0 -n /data/users/jeremy/simulation/syncsim/pol/meta/index.ndx
