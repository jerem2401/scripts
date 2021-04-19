#! /bin/bash

PWD=$(pwd)
if [ "$HOSTNAME" == smaug ]; then
	base=$(echo ${PWD%/simulation*})
else
	module load anaconda3/2020.07
	source activate env1
	base=$HOME
fi

while [ $# -gt 0 ]; do
	case "$1" in
		-d)
		  shift
		  dir=$1
		  ;;
		-mol)
		  shift
		  mol=$1;;
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
	tpr="$dir/md.tpr"
fi

if [[ "$dir" == '' ]]; then
	echo "no dir given, exiting"
	exit
fi

echo "dir is: $dir, skip is: $skip, group is: $group, protocol is: $proto, traj is: $traj, tpr is: $tpr"
sleep 1

if [[ $mol == 't4' ]]; then
	echo 1 | gmx trjconv -nice 0 -f $traj -s $tpr -o $dir/nopbc1.xtc -pbc atom -ur compact
	wait $!
	echo 1 | gmx trjconv -nice 0 -f $traj -s $tpr -o $dir/nopbc2.xtc -pbc whole
	wait $!
	echo 1 | gmx trjconv -nice 0 -f $dir/nopbc2.xtc -s $tpr -dump 0 -o $dir/0.gro
fi

#TODO put -dump and correct-chainid-and-ter.py lines in function as it is used several times
if [[ $mol == 'pol' ]]; then
	index="$base/simulation/syncsim/pol/ref/cc_ZN_noTFIIS_params/prep/nvt/index.ndx"
	[ -z "$group" ] && group=26
	[ -z "$skip" ] && skip=4
	[ -z "$proto" ] && proto=1

	if (("$proto" == 1)); then
		echo $group | gmx trjconv -nice 0 -f $traj -s $tpr -o $dir/nopbc1.xtc -ur compact -pbc atom -n $index -skip $skip
		wait
		echo $group | gmx trjconv -nice 0 -f $dir/nopbc1.xtc -s $tpr -o $dir/nopbc2.xtc -ur compact -pbc whole -n $index
		wait
		echo $group | gmx trjconv -nice 0 -f $dir/nopbc2.xtc -s $tpr -o $dir/0.pdb -dump 0 -n $index
		correct-chainid-and-ter.py $dir/0.pdb > $dir/0_chains.pdb
	elif [ "$proto" = path ]; then
		echo $group | gmx trjconv -nice 0 -f $traj -s $tpr -o $dir/nopbc1.xtc -ur compact -pbc mol -n $index -skip $skip
		wait
		echo $group | gmx trjconv -nice 0 -f $dir/nopbc1.xtc -s $tpr -o $dir/nopbc2.xtc -pbc nojump -n $index
		wait
		echo $group | gmx trjconv -nice 0 -f $dir/nopbc2.xtc -s $tpr -o $dir/0.pdb -dump 0 -n $index
		correct-chainid-and-ter.py $dir/0.pdb > $dir/0_chains.pdb
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
		echo $group | gmx trjconv -nice 0 -f $dir/nopbc2.xtc -s $tpr -o $dir/0.pdb -dump 0 -n $index
		correct-chainid-and-ter.py $dir/0.pdb > $dir/0_chains.pdb
	elif [ "$proto" = umb ]; then
		index="$base/simulation/syncsim/pol/heavy_h/ref/protwithin3ofdna_dna_protref.ndx"
		echo $group | gmx trjconv -nice 0 -f $traj -s $tpr -o $dir/umb1.xtc -ur compact -pbc mol -n $index -skip $skip ; wait
		echo $group | gmx convert-tpr -s $tpr -n $index -o $dir/subset.tpr
		echo 0 | gmx trjconv -nice 0 -f $dir/umb1.xtc -s $dir/subset.tpr -o $dir/umb2.xtc -pbc nojump -ur compact; wait
		rm -f $dir/umb1.xtc
		echo 0 | gmx trjconv -nice 0 -f $dir/umb2.xtc -s $dir/subset.tpr -o $dir/0umb.pdb -ur compact -dump 0; wait
		correct-chainid-and-ter.py $dir/0umb.pdb > $dir/0umb_chains.pdb
	fi

fi

if [[ $mol == 'opro' ]]; then
	echo "1 $group" | gmx trjconv -nice 0 -s $tpr -f $traj -o $dir/nopbc.xtc -pbc mol -center
	echo $group | gmx trjconv -nice 0 -f $dir/nopbc.xtc -s $tpr -o $dir/0.pdb -dump 0
fi
