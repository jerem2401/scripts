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

if [ -z "$2" ]; then
	skip=4
else
	skip=$2
fi

echo $2

if [ -z "$3" ]; then
	group=26
else
	group=$3
fi

echo $3

read -p "which index.ndx, 0 (old) or 1(new) ?" ival

if (( "$ival" == 0 )); then
	index="/data/users/jeremy/simulation/syncsim/pol/meta/index.ndx"
elif (( "$ival" == 1 )); then
	index="/data/users/jeremy/simulation/syncsim/pol/ref/cc_ZN_noTFIIS_params/prep/nvt/index.ndx"
else
	echo "error: no index file choosen"
	exit 1
fi

echo $group | gmx trjconv -f $1/traj_comp.xtc -s $1/md.tpr -o $1/nopbc1.xtc -ur compact -pbc atom -n $index -skip $skip
wait
echo $group | gmx trjconv -f $1/nopbc1.xtc -s $1/md.tpr -o $1/nopbc2.xtc -ur compact -pbc whole -n $index
wait
echo $group | gmx trjconv -f $1/nopbc2.xtc -s $1/md.tpr -o $1/0.pdb -dump 0 -n $index

correct-chainid-and-ter.py $1/0.pdb > $1/0_chains.pdb
