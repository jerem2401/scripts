#! /bin/bash

PWD=$(pwd)

if [[ "$PWD" =~ .*t4.* ]]; then
	i=$1
	echo $i
	echo 1 | gmx trjconv -f $i/traj_comp.xtc -s $i/md.tpr -o $i/nopbc1.xtc -pbc atom -ur compact
	wait $!
	echo 1 | gmx trjconv -f $i/nopbc1.xtc -s $i/md.tpr -o $i/nopbc2.xtc -pbc whole
	wait $!
	echo 1 | gmx trjconv -f $i/nopbc2.xtc -s $i/md.tpr -dump 0 -o $i/0.gro
fi


if [[ "$PWD" =~ .*pol.* ]]; then

	skip=4
	group=26
	index="$HOME/simulation/syncsim/pol/ref/cc_ZN_noTFIIS_params/prep/nvt/index.ndx"
	proto=1

	while [ $# -gt 0 ]; do
		case "$1" in
			-d)
			  shift
			  dir=$1
			  ;;
			-skip)
			  shift
			  skip=$1
			  ;;
			-group)
			  shift
			  group=$1
			  ;;
			-protocol)
			  shift
			  proto=$1
			  ;;
		esac
		shift
	done
	echo "dir is: $dir, skip is: $skip, group is: $group, protocol is: $proto"

	if [[ "$dir" == '' ]]; then
		echo "no dir given, exiting"
		exit
	fi

	if (("$proto" == 1)); then
		echo $group | gmx trjconv -f $dir/traj_comp.xtc -s $dir/md.tpr -o $dir/nopbc1.xtc -ur compact -pbc atom -n $index -skip $skip
		wait
		echo $group | gmx trjconv -f $dir/nopbc1.xtc -s $dir/md.tpr -o $dir/nopbc2.xtc -ur compact -pbc whole -n $index
		wait
	else
		echo $group | gmx trjconv -f $dir/traj_comp.xtc -s $dir/md.tpr -o $dir/nopbc1.xtc -ur compact -pbc mol -n $index -skip $skip
		wait
		echo $group | gmx trjconv -f $dir/nopbc1.xtc -s $dir/md.tpr -o $dir/nopbc2.xtc -pbc nojump -n $index
		wait
	fi

	module load conda
	source activate env1

	echo $group | gmx trjconv -f $dir/nopbc2.xtc -s $dir/md.tpr -o $dir/0.pdb -dump 0 -n $index
	correct-chainid-and-ter.py $dir/0.pdb > $dir/0_chains.pdb
fi
