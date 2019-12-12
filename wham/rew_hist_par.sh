#!/bin/bash

dirl=$(command ls -d ../E_*/colv* | sort -n -t _ -k 2)
#t=$(echo $dirl)
#IFS=' ' read -r -a ADDR <<< "$t" #dirl is read into an array as tokens separated by IFS

#echo ${ADDR[*]}
#echo ${ADDR[@]:0:10}
#echo "${#ADDR[@]}"
#for i in ${ADDR[@]:0:10}; do
#	echo $i
#done

for i in $dirl; do
	rew_hist_final.py $i;
done


