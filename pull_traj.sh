#! /bin/bash

PWD=$(pwd)
if [ "$HOSTNAME" == smaug ]; then
	base=$(echo ${PWD%/simulation*})
else
	module load anaconda3/2020.07
	source activate env1
	base=$HOME
fi

if [[ "$PWD" =~ .*t4.* ]]; then
	i=$1
	echo $i
	echo 1 | gmx trjconv -nice 0 -f $i/traj_comp.xtc -s $i/md.tpr -o $i/nopbc1.xtc -pbc atom -ur compact
	wait $!
	echo 1 | gmx trjconv -nice 0 -f $i/nopbc1.xtc -s $i/md.tpr -o $i/nopbc2.xtc -pbc whole
	wait $!
	echo 1 | gmx trjconv -nice 0 -f $i/nopbc2.xtc -s $i/md.tpr -dump 0 -o $i/0.gro
fi


if [[ "$PWD" =~ .*pol.* ]]; then

	skip=4
	group=26
	index="$base/simulation/syncsim/pol/ref/cc_ZN_noTFIIS_params/prep/nvt/index.ndx"
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
			-traj) shift
			  traj=$1;;
			-tpr) shift
			  tpr=$1;;
		esac
		shift
	done

	if [ -z "$traj" ]; then
		traj="$dir/traj_comp.xtc"
	fi
	if [ -z "$tpr" ]; then
		tpr="$tpr"
	fi

	echo "dir is: $dir, skip is: $skip, group is: $group, protocol is: $proto, traj is: $traj, tpr is: $tpr"

	if [[ "$dir" == '' ]]; then
		echo "no dir given, exiting"
		exit
	fi

	if (("$proto" == 1)); then
		echo $group | gmx trjconv -nice 0 -f $traj -s $tpr -o $dir/nopbc1.xtc -ur compact -pbc atom -n $index -skip $skip
		wait
		echo $group | gmx trjconv -nice 0 -f $dir/nopbc1.xtc -s $tpr -o $dir/nopbc2.xtc -ur compact -pbc whole -n $index
		wait
	elif [ "$proto" = path ]; then
		echo $group | gmx trjconv -nice 0 -f $traj -s $tpr -o $dir/nopbc1.xtc -ur compact -pbc mol -n $index -skip $skip
		wait
		echo $group | gmx trjconv -nice 0 -f $dir/nopbc1.xtc -s $tpr -o $dir/nopbc2.xtc -pbc nojump -n $index
		wait
	elif [ "$proto" = ocfree ]; then
		if [ $dir = rep3 ]; then
			echo $group | gmx trjconv -nice 0 -f $dir/traj_comp.xtc -s $tpr -o $dir/nopbc1.xtc -ur compact -pbc atom -n $index -skip $skip -trans 0 0 3
		elif [ $dir = rep5 ]; then
			echo $group | gmx trjconv -nice 0 -f $dir/traj_comp.xtc -s $tpr -o $dir/nopbc1.xtc -ur compact -pbc atom -n $index -skip $skip -trans -1 0 1
		elif [ $dir = rep8 ]; then
			echo $group | gmx trjconv -nice 0 -f $dir/traj_comp.xtc -s $tpr -o $dir/nopbc1.xtc -ur compact -pbc atom -n $index -skip $skip -trans 1 -4 -2
		fi
		wait
		echo $group | gmx trjconv -nice 0 -f $dir/nopbc1.xtc -s $tpr -o $dir/nopbc2.xtc -pbc whole -n $index
		wait
	fi

	echo $group | gmx trjconv -nice 0 -f $dir/nopbc2.xtc -s $tpr -o $dir/0.pdb -dump 0 -n $index
	correct-chainid-and-ter.py $dir/0.pdb > $dir/0_chains.pdb
fi
