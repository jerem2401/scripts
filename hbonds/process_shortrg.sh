#!/bin/bash
###############################################################################
#Author: Jeremy Lapierre <jeremy.lapierre@uni-saarland.de>
#Option description:
#
###############################################################################
set -o errexit   # abort on nonzero exitstatus
set -o nounset   # abort on unbound variable
set -o pipefail  # dont hide errors within pipes

#for i in srint_dna-dna.xvg; do
#    sed '1,25d' $i > ${i%.xvg}_noheader.xvg
#    awk '{ print $1, $2 + $3 }' ${i%.xvg}_noheader.xvg > ${i%.xvg}_noheader_sum.xvg
#    awk '{ if (NR == 1 ) { Value = $2 } { print $1, $2 - Value } }' ${i%.xvg}_noheader_sum.xvg  > ${i%.xvg}_zeroed.xvg
#done

#for i in srint_dna-water.xvg; do
#    #sed '1,30d' $i > ${i%.xvg}_noheader.xvg
#    sed '1,25d' $i > ${i%.xvg}_noheader.xvg
#    awk '{ print $1, $2 + $3 }' ${i%.xvg}_noheader.xvg > ${i%.xvg}_noheader_sum.xvg
#    awk '{ if (NR == 1 ) { Value = $2 } { print $1, $2 - Value } }' ${i%.xvg}_noheader_sum.xvg  > ${i%.xvg}_zeroed.xvg
#done

#for i in srint_prot-dna.xvg; do
#    #sed '1,27d' $i > ${i%.xvg}_noheader.xvg
#    #sed '1,25d' $i > ${i%.xvg}_noheader.xvg
#    #awk '{ print $1, $2 + $3 }' ${i%.xvg}_noheader.xvg > ${i%.xvg}_noheader_sum.xvg
#    awk '{ if (NR == 1 ) { Value = $2 } { print $1, $2 - Value } }' ${i%.xvg}_noheader_sum.xvg  > ${i%.xvg}_zeroed.xvg
#done

#for i in srint_dna-dna.xvg; do
#    sed '1,30d' $i > ${i%.xvg}_noheader.xvg
#    awk 'NR % 5 == 0' ${i%.xvg}_noheader.xvg > ${i%.xvg}_halved.xvg
#    awk '{ print $1, $2 + $3 }' ${i%.xvg}_halved.xvg > ${i%.xvg}_sum.xvg
#    awk '{ if (NR == 1 ) { Value = $2 } { print $1, $2 - Value } }' ${i%.xvg}_sum.xvg > ${i%.xvg}_zeroed.xvg
#done
#
#for i in srint_dna-water.xvg; do
#    sed '1,30d' $i > ${i%.xvg}_noheader.xvg
#    awk 'NR % 5 == 0' ${i%.xvg}_noheader.xvg > ${i%.xvg}_halved.xvg
#    awk '{ print $1, $2 + $3 }' ${i%.xvg}_halved.xvg > ${i%.xvg}_sum.xvg
#    awk '{ if (NR == 1 ) { Value = $2 } { print $1, $2 - Value } }' ${i%.xvg}_sum.xvg > ${i%.xvg}_zeroed.xvg
#done
#
#for i in srint_prot-dna.xvg; do
#    sed '1,30d' $i > ${i%.xvg}_noheader.xvg
#    awk 'NR % 5 == 0' ${i%.xvg}_noheader.xvg > ${i%.xvg}_halved.xvg
#    awk '{ print $1, $2 + $3 }' ${i%.xvg}_halved.xvg > ${i%.xvg}_sum.xvg
#    awk '{ if (NR == 1 ) { Value = $2 } { print $1, $2 - Value } }' ${i%.xvg}_sum.xvg > ${i%.xvg}_zeroed.xvg
#done
#done
#for j in $(seq 0 1 12); do
#        for i in srint_dna-water.xvg srint_prot-dna.xvg srint_dna-dna.xvg; do
#	    nbr=$((25+$j))
#	    sed "1,${nbr}d" $i > ${i%.xvg}_noheader.xvg
#	    awk 'NR % 10 == 0' ${i%.xvg}_noheader.xvg > ${i%.xvg}_halved.xvg
#	    awk '{ print $1, $2 + $3 }' ${i%.xvg}_halved.xvg > ${i%.xvg}_sum.xvg
#	    awk '{ if (NR == 1 ) { Value = $2 } { print $1, $2 - Value } }' ${i%.xvg}_sum.xvg > ${i%.xvg}_${j}_zeroed.xvg
#	done
#done
for i in rep*/short_range_2/{srint_prot-dna_noheader.xvg,srint_dna-dna_noheader.xvg,srint_dna-water_noheader.xvg};do tmp=$(basename $i); awk -v sig=14000 -f ~/gitrepo/scripts/hbonds/gauss1.awk < ${i} > $(dirname ${i})/${tmp%.xvg}_smoothed.xvg; echo $i; done
for i in rep*/short_range_2/{srint_prot-dna_noheader_smoothed.xvg,srint_dna-dna_noheader_smoothed.xvg,srint_dna-water_noheader_smoothed.xvg};do tmp=$(basename $i); awk '{ if (NR == 1 ) { Value = $2 } { print $1, $2 - Value } }' $i > $(dirname ${i})/${tmp%.xvg}_zeroed.xvg; echo $i; done
