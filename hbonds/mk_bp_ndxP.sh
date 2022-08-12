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
mknd=1
ana=''
index="$base/simulation/syncsim/pol/heavy_h/ref/hbondP.ndx"

while [ $# -gt 0 ]; do
	case "$1" in
		-d)
		  shift
		  dir="$1/hbondP"
		  dirumb=$1
		  ;;
		-ndx) shift
		  index=$1;;
		-ana) shift
		  ana=$1
		  mknd=0;;
	esac
	shift
done


if (($mknd == 1)); then

	#[[ -d $dir ]] && echo "$dir already exists, exiting" && exit || mkdir $dir
	#if ls $dir/DNA* &> /dev/null; then echo "DNAs dir exist"; exit; fi

	#while read line; do
	#	f=$(echo $line | awk '{print $1}')
	#	fn=$(echo $line | awk '{print $2}')
	#	fa3=''
	#	sa3=''
	#	if [ $fn == "DC" ]; then
	#		fa1="atomname N3"
	#		fa2="atomname O2"
	#		fa3="(atomname N4 or atomname H41 or atomname H42)"
	#		sa1="(atomname N1 or atomname H1)"
	#		sa2="(atomname N2 or atomname H21 or atomname H22)"
	#		sa3="atomname O6"
	#	elif [ $fn == "DG" ]; then
	#		fa1="(atomname N1 or atomname H1)"
	#		fa2="(atomname N2 or atomname H21 or atomname H22)"
	#		fa3="atomname O6"
	#		sa1="atomname N3"
	#		sa2="atomname O2"
	#		sa3="(atomname N4 or atomname H41 or atomname H42)"
	#	elif [ $fn == "DA" ]; then
	#		fa1="(atomname N6 or atomname H61 or atomname H62)"
	#		fa2="atomname N1"
	#		sa1="atomname O4"
	#		sa2="(atomname N3 or atomname H3)"
	#	elif [ $fn == "DT" ]; then
	#		fa1="atomname O4"
	#		fa2="(atomname N3 or atomname H3)"
	#		sa1="(atomname N6 or atomname H61 or atomname H62)"
	#		sa2="atomname N1"
	#	else
	#		echo "error, exiting"
	#		exit
	#	fi
	#	s=$(echo $line | awk '{print $3}')
	#	#sn=$(echo $line | awk '{print $4}')


	#	if [[ ! -z $fa3 ]]; then
	#		mkdir "${dir}/DNAt_${f}"
	#		sel=$(echo "group \"Protein\" and within \
	#		0.4 of (group \"DNAt\" and resid $f and ($fa1 or $fa2 or $fa3));")
	#		echo $sel > "${dir}/DNAt_${f}/selection.dat"
	#		gmx select -sf "${dir}/DNAt_${f}/selection.dat" -f ${dirumb}/umb2.xtc \
	#		-s ${dirumb}/subset.tpr -n $index -on ${dir}/DNAt_${f}/hbindex.ndx

	#		mkdir "${dir}/DNAnt_${s}"
	#		sel=$(echo "group \"Protein\" and within \
	#		0.4 of (group \"DNAnt\" and resid $s and ($sa1 or $sa2 or $sa3));")
	#		echo $sel > "${dir}/DNAnt_${s}/selection.dat"
	#		gmx select -sf "${dir}/DNAnt_${s}/selection.dat" -f ${dirumb}/umb2.xtc \
	#		-s ${dirumb}/subset.tpr -n $index -on ${dir}/DNAnt_${s}/hbindex.ndx
	#	else
	#		mkdir "${dir}/DNAt_${f}"
	#		sel=$(echo "group \"Protein\" and within \
	#		0.4 of (group \"DNAt\" and resid $f and ($fa1 or $fa2));")
	#		echo $sel > "${dir}/DNAt_${f}/selection.dat"
	#		gmx select -sf "${dir}/DNAt_${f}/selection.dat" -f ${dirumb}/umb2.xtc \
	#		-s ${dirumb}/subset.tpr -n $index -on ${dir}/DNAt_${f}/hbindex.ndx

	#		mkdir "${dir}/DNAnt_${s}"
	#		sel=$(echo "group \"Protein\" and within \
	#		0.4 of (group \"DNAnt\" and resid $s and ($sa1 or $sa2));")
	#		echo $sel > "${dir}/DNAnt_${s}/selection.dat"
	#		gmx select -sf "${dir}/DNAnt_${s}/selection.dat" -f ${dirumb}/umb2.xtc \
	#		-s ${dirumb}/subset.tpr -n $index -on ${dir}/DNAnt_${s}/hbindex.ndx
	#	fi


	#done < $base/simulation/syncsim/pol/heavy_h/ref/pairs2P.txt

	for i in ${dir}/DNA*; do
		atominrange=$(egrep -c -v '(\[|^$)' $i/hbindex.ndx)
		if (($atominrange == 0)); then
			echo "no atom in range for $i"
			echo -e "#@ s0 legend \"av 0.0\"\n#@ s1 legend \"ee  0.0\"" > "${i}/hbond.xvg"
		else
			nbndx=$(grep -c '\[' $i/hbindex.ndx)
			nbndxf=$(($nbndx-1))
			mergegrp=$(for i in $(seq 0 $nbndxf); do echo -n "group $i or "; done)
			pgrp=$(echo ${mergegrp:0: -4})
			echo "$pgrp;" > "${i}/merge_sel.dat"
			gmx select -sf "${i}/merge_sel.dat" -f ${dirumb}/umb2.xtc \
			-s ${dirumb}/subset.tpr -n $i/hbindex.ndx -on ${i}/hbindexm.ndx
			sed '1 s/^.*$/\[ my_custom_group_along_traj \]/' hbindexm.ndx

			grep -oP '(?<=of \().*(?=\);)' $i/selection.dat > ${i}/selectionnuc.dat
			gmx select -sf ${i}/selectionnuc.dat -f ${dirumb}/umb2.xtc \
			-s ${dirumb}/subset.tpr -n $index -on ${i}/hbindexdna.ndx

			cat ${i}/hbindexm.ndx ${i}/hbindexdna.ndx > ${i}/final.ndx
			echo "0 1" | gmx hbond -nthreads 1 -f ${dirumb}/umb2.xtc \
			-s ${dirumb}/subset.tpr -n ${i}/final.ndx -num "${i}/hbond.xvg"
			wait
			gmx analyze -f ${i}/hbond.xvg -ee ${i}/err_hbond.xvg &> ${i}/out.txt
			wait
			av=$(grep "\@ s0" "${i}/err_hbond.xvg")
			std=$(cat ${i}/out.txt | grep SS1 | awk '{print $3}')
			err=$(grep "\@ s1" "${i}/err_hbond.xvg")
			sed  -i "1i \#std $std" "${i}/hbond.xvg"
			sed  -i "1i \#$err" "${i}/hbond.xvg"
			sed  -i "1i \#$av " "${i}/hbond.xvg"
		fi
	done
elif [[ ! -z $ana ]]; then
	set -o noglob
	statsP.py -f $ana
	set +o noglob
fi
#for i in E*/hbond/analyzehb.txt; do a=$(dirname $i); E=$(echo ${a:0:-6}); hb=$(tail -n1 $i ); real=$(awk '{ total += $2 } END { print total/NR }' $E/colvar_${E:2}.txt); printf "%-8s %20s %20s \n" $E $real $hb >> hbout.txt; done
