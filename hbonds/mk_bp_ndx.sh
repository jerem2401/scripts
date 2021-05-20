#!/bin/bash
###############################################################################
#Author: Jeremy Lapierre <jeremy.lapierre@uni-saarland.de>
#Option description:
#
###############################################################################
#set -o errexit   # abort on nonzero exitstatus
#set -o nounset   # abort on unbound variable
#set -o pipefail  # dont hide errors within pipes

PWD=$(pwd)
if [ "$HOSTNAME" == smaug ]; then
	base=$(echo ${PWD%/simulation*})
else
	module load anaconda3/2020.07 && source activate env1
	base=$HOME
fi

traj=''
dir=''
mknd=0
traj_hb=1
ana=''
index="$base/simulation/syncsim/pol/heavy_h/ref/hbond.ndx"

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
		-traj) shift
		  traj=$1;;
		-tpr) shift
		  tpr=$1;;
		-mknd)
		  mknd=1
		  traj_hb=0;;
		-ndx) shift
		  index=$1;;
		-ana) shift
		  ana=$1
		  traj_hb=0;;
	esac
	shift
done


if (($mknd == 1)); then
	while read line; do
		f=$(echo $line | awk '{print $1}')
		fn=$(echo $line | awk '{print $2}')
		fa3=''
		sa3=''
		if [ $fn == "DC" ]; then
			fa1="atomname N3"
			fa2="atomname O2"
			fa3="(atomname N4 or atomname H41 or atomname H42)"
			sa1="(atomname N1 or atomname H1)"
			sa2="(atomname N2 or atomname H21 or atomname H22)"
			sa3="atomname O6"
		elif [ $fn == "DG" ]; then
			fa1="(atomname N1 or atomname H1)"
			fa2="(atomname N2 or atomname H21 or atomname H22)"
			fa3="atomname O6"
			sa1="atomname N3"
			sa2="atomname O2"
			sa3="(atomname N4 or atomname H41 or atomname H42)"
		elif [ $fn == "DA" ]; then
			fa1="(atomname N6 or atomname H61 or atomname H62)"
			fa2="atomname N1"
			sa1="atomname O4"
			sa2="(atomname N3 or atomname H3)"
		elif [ $fn == "DT" ]; then
			fa1="atomname O4"
			fa2="(atomname N3 or atomname H3)"
			sa1="(atomname N6 or atomname H61 or atomname H62)"
			sa2="atomname N1"
		else
			echo "error, exiting"
			exit
		fi
		s=$(echo $line | awk '{print $3}')
		#sn=$(echo $line | awk '{print $4}')

		sel=$(echo "(group \"DNAt\" and resid $f and $fa1) or (group \"DNAnt\" and resid $s and $sa1);")
		echo $sel >> selection.dat
		sel=$(echo "(group \"DNAt\" and resid $f and $fa2) or (group \"DNAnt\" and resid $s and $sa2);")
		echo $sel >> selection.dat
		if [[ ! -z $fa3 ]]; then
			sel=$(echo "(group \"DNAt\" and resid $f and $fa3) or (group \"DNAnt\" and resid $s and $sa3);")
			echo $sel >> selection.dat
		fi
	done < $base/simulation/syncsim/pol/heavy_h/ref/pairs2.txt

	gmx select -sf selection.dat -s md.tpr -n $base/simulation/syncsim/pol/heavy_h/ref/index_dna.ndx -on hbindex.ndx

elif (($traj_hb == 1)); then
	if [ -z $dir ]; then
		echo "you should give target dir, exiting"
		exit
	fi

	if [ -z "$traj" ]; then
		traj="$dir/traj_comp.xtc"
	fi

	[[ -d ${dir}/hbond ]] && echo "${dir}/hbond already exists, exiting" && exit || mkdir ${dir}/hbond

	index_1="$base/simulation/syncsim/pol/heavy_h/ref/index.ndx"
	echo DNA | gmx trjconv -nice 0 -f $traj -s $tpr -o $dir/nopbc1.xtc -ur compact -pbc mol \
	-n $index_1 -skip $skip
	echo DNA | gmx convert-tpr -s $tpr -n $index_1 -o $dir/subset.tpr
	tpr=$dir/subset.tpr
	wait
	echo 0 | gmx trjconv -nice 0 -f $dir/nopbc1.xtc -s $tpr -o $dir/nopbc2.xtc -pbc nojump
	wait
	echo 0 | gmx trjconv -nice 0 -f $dir/nopbc2.xtc -s $tpr -o $dir/0.pdb -dump 0
	correct-chainid-and-ter.py $dir/0.pdb > $dir/0_chains.pdb

	echo "trajconv of $dir is done" >> hbout.txt

	for i in $(seq 77 1 118); do
		line=$(($i + 1))
		resn=$(grep '\[' $index | sed "${line}q;d" | grep -oP '(?<=resid_)[^_\]\ \)]*')
		temp=$(echo $resn | awk '{print $1}')
		ntemp=$(echo $resn | awk '{print $2}')
		atn=$(grep '\[' $index | sed "${line}q;d" | grep -oP '(?<=atomname_)[^_\]\ \)]*')
		atnb=$(echo $atn | sed 's/ /_/g')
		ksi=$(echo $dir | grep -oP 'E_.*(?=\/)')
		name="${ksi}_t${temp}_nt${ntemp}_${atnb}"
		#echo "$i $i" | gmx hbond -nthreads 1 -f "$dir/nopbc2.xtc" -s "$tpr" -n "$index"  \
		#-num "${dir}/${name}.xvg" -r 0.358 &
		echo "$i $i" | gmx hbond -nthreads 1 -f "$dir/nopbc2.xtc" -s "$tpr" -n "$index"  \
		-num "${dir}/${name}.xvg"
		wait
		gmx analyze -f "${dir}/${name}.xvg" -ee "${dir}/err_${name}.xvg" &> "${dir}/out.txt"
		wait
		av=$(grep "\@ s0" "${dir}/err_${name}.xvg")
		std=$(cat "${dir}/out.txt" | grep SS1 | awk '{print $3}')
		err=$(grep "\@ s1" "${dir}/err_${name}.xvg")
		sed  -i "1i \#std $std" "${dir}/${name}.xvg"
		sed  -i "1i \#$err" "${dir}/${name}.xvg"
		sed  -i "1i \#$av " "${dir}/${name}.xvg"
	done

	echo "hbond of traj of $dir is done" >> hbout.txt

elif [[ ! -z $ana ]]; then
	set -o noglob
	stats.py -f $ana
	set +o noglob
fi

#for i in E*/colv*; do dir=$(dirname $i); echo ${dir:2}; mk_bp_ndx.sh -d $dir -skip 5 -tpr "${dir}/conf_${dir:2}.tpr" -ndx /data/users/jeremy/simulation/syncsim/pol/heavy_h/ref/hbond_dnaonly.ndx; done
#for i in E*/hbond/analyzehb.txt; do a=$(dirname $i); E=$(echo ${a:0:-6}); hb=$(tail -n1 $i ); real=$(awk '{ total += $2 } END { print total/NR }' $E/colvar_${E:2}.txt); printf "%-8s %20s %20s \n" $E $real $hb >> hbout.txt; done
