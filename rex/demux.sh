#!/bin/bash
###############################################################################
#Author: Jeremy Lapierre <jeremy.lapierre@uni-saarland.de>
#Option description:
#
###############################################################################
set -o errexit   # abort on nonzero exitstatus
set -o nounset   # abort on unbound variable
set -o pipefail  # dont hide errors within pipes

while [ $# -gt 0 ]; do
    case "$1" in
        -f)
            shift
            rpath="$1";;
    esac
    shift
done

mkdir ${rpath}/analysis/demux/{neg,mid,pos}
logd=$(ls -d ${rpath}/mid/E* | head -n1)
cp ${logd}/md.log ${rpath}/analysis/demux/mid
(cd ${rpath}/analysis/demux/mid; echo 0.004 | ~/gitrepo/demux/demux.pl ./md.log) &
for i in neg pos; do
    dirs=$(ls -d ${rpath}/${i}/*/ | grep -oP "(?<=${i}/).*(?=/)")
    for j in ${dirs}; do
        mkdir ${rpath}/analysis/demux/${i}/${j}
	logd=$(ls -d ${rpath}/${i}/${j}/E* | head -n1)
	cp ${logd}/md.log ${rpath}/analysis/demux/${i}/${j}
	(cd ${rpath}/analysis/demux/${i}/${j}; echo 0.004 | ~/gitrepo/demux/demux.pl ./md.log) &
    done
done

wait
module load anaconda3/2020.07 && source activate env1
demux_plot.py -f ${rpath}/analysis/demux/{neg/*/,mid/,pos/*/}
#rm -f plot.txt; for i in $(seq 2 $((n+2))); do echo ${mid_arr[$i]} >> plot.txt; echo $((i-2)) >> plot.txt; gnuplot -e "set term dumb size 70,${n_2}; set yr [0:${n}]; set ytics 1; p 'replica_temp.xvg' u 1:${i} pt '*'" >> plot.txt; done
