#!/bin/bash

touch fname.dat
>fname.dat

for dir in "$1"E_*
do
	echo "$dir"
	val=$(echo "$dir" | grep -oP '(?<=E_).*')
	echo "$val"
	kappa=$(echo "$dir" | grep -oP '(?<=_k).*(?=_a)')
	time=$(cat "${dir}/colvar_${val}.txt" | awk 'FNR > 5 {print $1}')
	distance=$(cat "${dir}/colvar_${val}.txt" | awk -v v="$val" 'FNR > 5 {print ($6 - v)}')
	paste <(echo "$time") <(echo "$distance") > "./E_${val}.pdo"
        echo "# UMBRELLA      3.0
# Component selection: 0 0 1
# nSkip 1
# Ref. Group TestAtom
# Nr. of pull groups 1
# Group 1 GR1 Umb. Pos. $val Umb. Cons. $kappa
#####" | cat - "./E_${val}.pdo" > temp && mv temp "./E_${val}.pdo"
	gzip "./E_${val}.pdo"
	echo "E_${val}.pdo.gz" >> fname.dat
done


