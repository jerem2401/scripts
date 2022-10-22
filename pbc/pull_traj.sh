#! /bin/bash

PWD=$(pwd)
if [ "$HOSTNAME" == smaug ]; then
	base=$(echo ${PWD%/simulation*})
else
	module load anaconda3/2020.07
	source activate env1
	base=$HOME
fi

nstep=''

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
		-group2)
		  shift
		  group2=$1
		  ;;
		-protocol)
		  shift
		  proto=$1
		  ;;
		-traj) shift
		  traj=$1;;
		-tpr) shift
		  tpr=$1;;
		-n) shift
		  index=$1;;
		-nstep) shift
		  nstep=' -nsteps 1';;
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
	[ -z "$index" ] && index="$base/simulation/syncsim/pol/ref/cc_ZN_noTFIIS_params/prep/nvt/index.ndx"
	[ -z "$group" ] && group=26
	[ -z "$skip" ] && skip=4
	[ -z "$proto" ] && proto=1

	if (("$proto" == 1)); then
		echo $group | gmx trjconv -nice 0 -f $traj -s $tpr -o $dir/nopbc1.xtc -ur compact -pbc atom -n $index -skip $skip
		echo $group | gmx convert-tpr $nstep -s $tpr -n $index -o $dir/convert.tpr
		echo $group | gmx trjconv -nice 0 -f $dir/nopbc1.xtc -s $dir/convert.tpr -o $dir/nopbc2.xtc -ur compact -pbc whole -n $index
		echo $group | gmx trjconv -nice 0 -f $dir/nopbc2.xtc -s $dir/convert.tpr -o $dir/0.pdb -dump 0 -n $index
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
		echo $group | gmx convert-tpr $nstep -s $tpr -n $index -o $dir/subset.tpr
		echo 0 | gmx trjconv -nice 0 -f $dir/umb1.xtc -s $dir/subset.tpr -o $dir/umb2.xtc -pbc nojump -ur compact; wait
		rm -f $dir/umb1.xtc
		echo 0 | gmx trjconv -nice 0 -f $dir/umb2.xtc -s $dir/subset.tpr -o $dir/0umb.pdb -ur compact -dump 0; wait
		correct-chainid-and-ter.py $dir/0umb.pdb > $dir/0umb_chains.pdb
	elif [ "$proto" = test  ]; then
		index="$base/simulation/syncsim/pol/heavy_h/ref/master_index.ndx"
		echo $group | gmx trjconv -f $traj -s $tpr -ur compact -pbc atom -skip $skip -o gold1.xtc -n $index -trans 5.5 7 -4
		echo $group | gmx convert-tpr $nstep -s $tpr -n $index -o $dir/subset.tpr
		echo 0 | gmx trjconv -nice 0 -f $dir/gold1.xtc -s $dir/subset.tpr -o $dir/0.pdb -ur compact -dump 0
		correct-chainid-and-ter.py $dir/0.pdb > $dir/0_chains.pdb
		#echo "13 31" | gmx trjconv -f $traj -s $tpr -ur compact -pbc atom -skip $skip -o $dir/gold1.xtc -n ./E_87.870/test/master_index.ndx -center
		#echo 31 | gmx convert-tpr $nstep -s $tpr -n ./E_87.870/test/master_index.ndx -o $dir/subset.tpr
		#echo 0 | gmx trjconv -f $dir/gold1.xtc -s $dir/subset.tpr -ur compact -pbc whole -o $dir/gold2.xtc
	elif [ "$proto" = test2 ]; then
		index="$base/simulation/syncsim/pol/heavy_h/ref/index_ions.ndx"
		echo "29 27" | gmx trjconv -nice 0 -f $traj -s $tpr -o $dir/test2.xtc -ur compact -pbc atom -n $index -skip $skip -center
		echo 27 | gmx convert-tpr $nstep -s $tpr -n $index -o $dir/test2.tpr
		echo 0 | gmx trjconv -nice 0 -f $dir/test2.xtc -s $dir/test2.tpr -o $dir/test2_2.xtc -ur compact -pbc whole
		echo 0 | gmx trjconv -nice 0 -f $dir/test2_2.xtc -s $dir/test2.tpr -o $dir/test2.pdb -dump 0
	elif [ "$proto" = ions ]; then
		index="$base/simulation/syncsim/pol/heavy_h/ref/index_ions.ndx"
		echo "29 0" | gmx trjconv -nice 0 -f $traj -s $tpr -o $dir/solvant.xtc -ur compact -pbc atom -n $index -skip $skip -center
		echo 0 | gmx trjconv -nice 0 -f $dir/solvant.xtc -s $tpr -o $dir/solvant2.xtc -ur compact -pbc whole -n $index
		echo 0 | gmx trjconv -nice 0 -f $dir/solvant2.xtc -s $tpr -o $dir/solvant.pdb -dump 0 -n $index
	fi

fi

if [[ $mol == 'opro' ]]; then
	[ -z "$group" ] && group='Protein'
	[ -z "$proto" ] && proto=1

	if (("$proto" == 1)); then
		echo "$group $group2" | gmx trjconv -nice 0 -s $tpr -f $traj -o $dir/nopbc.xtc -pbc mol -center -skip $skip -n $index
		echo $group2 | gmx convert-tpr -s $tpr -n $index -o $dir/convert.tpr
		echo 0 | gmx trjconv -nice 0 -f $dir/nopbc.xtc -s $dir/convert.tpr -o $dir/0.gro -dump 0
	elif [ "$proto" = rot ]; then
		index="${base}/simulation/syncsim/opro/finalsetup/umb_wrex/po3/ref/index_prot-fom.ndx"
		[ -z "$group2" ] && group2='FOM'
		echo "$group" "$group2" | gmx trjconv -nice 0 -s $tpr -f $traj -o $dir/fom_nopbc.xtc -pbc mol -center -skip $skip -n $index
		echo $group2 | gmx convert-tpr -s $tpr -n $index -o $dir/convert.tpr
		echo 0 | gmx trjconv -nice 0 -f $dir/fom_nopbc.xtc -s $dir/convert.tpr -o $dir/fom.gro -dump 0
	fi
fi
