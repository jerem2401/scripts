#! /bin/bash

PWD=$(pwd)
if [ "$HOSTNAME" == smaug ]; then
	base=$(echo ${PWD%/simulation*})
else
	module load conda
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
		esac
		shift
	done

	xtc=$(command ls $dir/*.xtc)
	nxtc=$(echo $xtc | wc -w)

	if (($nxtc > 1)); then
		echo "more than 1 xtc in $dir, please chose 1 file among: $xtc"
		read traj
	else
		traj=$xtc
	fi

	echo "dir is: $dir, skip is: $skip, group is: $group, protocol is: $proto, traj is: $traj"

	if [[ "$dir" == '' ]]; then
		echo "no dir given, exiting"
		exit
	fi

	if (("$proto" == 1)); then
		echo $group | gmx trjconv -nice 0 -f $traj -s $dir/md.tpr -o $dir/nopbc1.xtc -ur compact -pbc atom -n $index -skip $skip
		wait
		echo $group | gmx trjconv -nice 0 -f $dir/nopbc1.xtc -s $dir/md.tpr -o $dir/nopbc2.xtc -ur compact -pbc whole -n $index
		wait
	elif [ "$proto" = path ]; then
		echo $group | gmx trjconv -nice 0 -f $traj -s $dir/md.tpr -o $dir/nopbc1.xtc -ur compact -pbc mol -n $index -skip $skip
		wait
		echo $group | gmx trjconv -nice 0 -f $dir/nopbc1.xtc -s $dir/md.tpr -o $dir/nopbc2.xtc -pbc nojump -n $index
		wait
	elif [ "$proto" = ocfree ]; then
		if [ $dir = rep3 ]; then
			echo $group | gmx trjconv -nice 0 -f $dir/traj_comp.xtc -s $dir/md.tpr -o $dir/nopbc1.xtc -ur compact -pbc atom -n $index -skip $skip -trans 0 0 3
		elif [ $dir = rep5 ]; then
			echo $group | gmx trjconv -nice 0 -f $dir/traj_comp.xtc -s $dir/md.tpr -o $dir/nopbc1.xtc -ur compact -pbc atom -n $index -skip $skip -trans -1 0 1
		elif [ $dir = rep8 ]; then
			echo $group | gmx trjconv -nice 0 -f $dir/traj_comp.xtc -s $dir/md.tpr -o $dir/nopbc1.xtc -ur compact -pbc atom -n $index -skip $skip -trans 1 -4 -2
		fi
		wait
		echo $group | gmx trjconv -nice 0 -f $dir/nopbc1.xtc -s $dir/md.tpr -o $dir/nopbc2.xtc -pbc whole -n $index
		wait
	fi

	echo $group | gmx trjconv -nice 0 -f $dir/nopbc2.xtc -s $dir/md.tpr -o $dir/0.pdb -dump 0 -n $index
	correct-chainid-and-ter.py $dir/0.pdb > $dir/0_chains.pdb
fi
