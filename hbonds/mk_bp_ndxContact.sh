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
basic=0
acidic=0
hydrophob=0
hydrophil=0
ndx="-n ${base}/simulation/syncsim/pol/heavy_h/ref/moreindex/umb_dnabubble_protein_2.ndx"

while [ $# -gt 0 ]; do
	case "$1" in
		-d)
		  shift
		  dirumb=$1
		  ;;
		-basic)
		  basic=1;;
		-acidic)
		  acidic=1;;
		-hydrophob)
		  hydrophob=1;;
		-hydrophil)
		  hydrophil=1;;
		-ana) shift
		  ana=$1
		  mknd=0;;
	esac
	shift
done

if (($basic == 1)); then
	dir="${dirumb}/contact4_basic"
	grp1=0
	grp2=1
elif (($acidic == 1)); then
	dir="${dirumb}/contact4_acidic"
	grp1=0
	grp2=1
elif (($hydrophob == 1)); then
	dir="${dirumb}/contact4_hydrophob"
	grp1=0
	grp2=1
elif (($hydrophil == 1)); then
	dir="${dirumb}/contact4_hydrophil"
	grp1=0
	grp2=1
else
	dir="${dirumb}/contact4"
	grp1=1
	grp2=3
fi

traj='umb2.xtc'
tpr='subset.tpr'


if (($mknd == 1)); then

	[[ -d $dir ]] && echo "$dir already exists, exiting" && exit || mkdir $dir

	for i in ${dir}; do
		if (($basic==1)); then
			gmx select -nice 0 -select \
			"resname ARG or resname LYS or resname HIS; group 3" \
			-f ${dirumb}/${traj} -s ${dirumb}/${tpr} -on ${i}/basic.ndx $ndx
			ndx="-n ${i}/basic.ndx"
		elif (($acidic==1)); then
			gmx select -nice 0 -select \
			"resname ASP or resname GLU; group 3" \
			-f ${dirumb}/${traj} -s ${dirumb}/${tpr} -on ${i}/acidic.ndx $ndx
			ndx="-n ${i}/acidic.ndx"
		elif (($hydrophob == 1)); then
			gmx select -nice 0 -select \
			"resname GLY or resname PRO or resname ALA or resname ILE or resname LEU or resname MET or resname PHE or resname TRP or resname TYR or resname VAL; group 3" \
			-f ${dirumb}/${traj} -s ${dirumb}/${tpr} -on ${i}/hydrophob.ndx $ndx
			ndx="-n ${i}/hydrophob.ndx"
		elif (($hydrophil == 1)); then
			gmx select -nice 0 -select \
			"resname SER or resname THR or resname ASN or resname GLN or resname CYS; group 3" \
			-f ${dirumb}/${traj} -s ${dirumb}/${tpr} -on ${i}/hydrophil.ndx $ndx
			ndx="-n ${i}/hydrophil.ndx"
		fi

		echo "$grp1 $grp2" | gmx mindist -f ${dirumb}/${traj} -s ${dirumb}/${tpr} -d 0.3 \
		-on "${i}/numcount.xvg" -od "${i}/mindist.xvg" $ndx
		gmx analyze -nice 0 -f ${i}/numcount.xvg -ee ${i}/err_numcount.xvg &> ${i}/out.txt
		wait
		av=$(grep "\@ s0" "${i}/err_numcount.xvg")
		std=$(cat ${i}/out.txt | grep SS1 | awk '{print $3}')
		err=$(grep "\@ s1" "${i}/err_numcount.xvg")
		sed  -i "1i \#std $std" "${i}/numcount.xvg"
		sed  -i "1i \#$err" "${i}/numcount.xvg"
		sed  -i "1i \#$av " "${i}/numcount.xvg"
	done

elif [[ ! -z $ana ]]; then
	set -o noglob
	statsctct.py -f $ana
	set +o noglob
fi
#for i in E*/hbond/analyzehb.txt; do a=$(dirname $i); E=$(echo ${a:0:-6}); hb=$(tail -n1 $i ); real=$(awk '{ total += $2 } END { print total/NR }' $E/colvar_${E:2}.txt); printf "%-8s %20s %20s \n" $E $real $hb >> hbout.txt; done
