#!/bin/bash

#touch fname.dat
#>fname.dat

rm ./fname.dat

while [ $# -gt 0  ]; do
    case "$1" in
        -colv) shift
	  colv=$1;;
    esac
    shift
done

if [[ -z $colv ]]; then
    echo "no colvar files given, exiting"
    exit
fi

for f in $colv; do
    echo "processing $f"
    val=$(echo "$f" | grep -oP '(?<=colvar_).*(?=.txt)')
    echo "ksi cnt is $val" #TODO add var for pld dir
    [[ -f ../../E_${val}/plumed.dat ]] && pld="../../E_${val}/plumed.dat" || \
    ( echo "no ../../E_${val}/plumed.dat, exiting" && exit ) #TODO add var for pld
    echo "check if $pld is the pld you want"
    kappa=$(grep -oP '(?<=   KAPPA=).*' "$pld")
    time=$(awk '{print $1}' $f)
    distance=$(awk -v v="$val" '{print ($2 -v)}' $f)
    paste <(echo "$time") <(echo "$distance") > "./E_${val}.pdo"
    echo "# UMBRELLA      3.0
# Component selection: 0 0 1
# nSkip 1
# Ref. Group 'TestAtom'
# Nr. of pull groups 1
# Group 1 'GR1'  Umb. Pos. $val Umb. Cons. $kappa
#####" | cat - "./E_${val}.pdo" > temp && mv temp "./E_${val}.pdo"
    gzip "./E_${val}.pdo"
    echo "E_${val}.pdo.gz" >> fname.dat
done
