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

dir="${dirumb}/contact"
traj='umb2.xtc'
tpr='subset.tpr'


if (($mknd == 1)); then

	[[ -d $dir ]] && echo "$dir already exists, exiting" && exit || mkdir $dir

	for i in ${dir}; do
		echo "1 12" | gmx mindist -nice 0 -f ${dirumb}/${traj} -s ${dirumb}/${tpr} -d 0.5 \
		-on "${i}/numcount.xvg" -od "${i}/mindist.xvg"
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
