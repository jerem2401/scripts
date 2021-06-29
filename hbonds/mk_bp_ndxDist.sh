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

while [ $# -gt 0 ]; do
	case "$1" in
		-d)
		  shift
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

dir="${dirumb}/dist"
traj='umb2.xtc'
tpr='subset.tpr'


if (($mknd == 1)); then

	[[ -d $dir ]] && echo "$dir already exists, exiting" && exit || mkdir $dir
	dir=$(realpath $dir)

	for i in ${dir}; do
		sed "s=DIR=${dir}=g" $base/simulation/syncsim/pol/heavy_h/ref/dist.dat \
		> ${dir}/dist.dat
		wait
		plumed driver --plumed ${dir}/dist.dat --mf_xtc ${dirumb}/${traj}
		gmx analyze -f ${i}/colvar_dist.xvg -ee ${i}/err_dist.xvg &> ${i}/out.txt
		av=$(grep "\@ s0" "${i}/err_dist.xvg")
		std=$(cat ${i}/out.txt | grep SS1 | awk '{print $3}')
		err=$(grep "\@ s1" "${i}/err_dist.xvg")
		sed  -i "1i \#std $std" "${i}/colvar_dist.xvg"
		sed  -i "1i \#$err" "${i}/colvar_dist.xvg"
		sed  -i "1i \#$av " "${i}/colvar_dist.xvg"
	done

elif [[ ! -z $ana ]]; then
	set -o noglob
	statsdist.py -f $ana
	set +o noglob
fi
#for i in E*/hbond/analyzehb.txt; do a=$(dirname $i); E=$(echo ${a:0:-6}); hb=$(tail -n1 $i ); real=$(awk '{ total += $2 } END { print total/NR }' $E/colvar_${E:2}.txt); printf "%-8s %20s %20s \n" $E $real $hb >> hbout.txt; done
