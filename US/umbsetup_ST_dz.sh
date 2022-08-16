#!/bin/bash

dir=$(d=$(pwd); basename $d)

nbOfWin=$(echo "$dir" | grep -oP '(?<=N).*(?=_t.*_k)')

ttot=$(echo "$dir" | grep -oP '(?<=_t).*(?=_k)')

kappa=$(echo "$dir" | grep -oP '(?<=_k).*(?=_rep)')

scale=5

tempering=0
wrex=0

while [ $# -gt 0 ]; do
    case "$1" in
	--file|-f)
            shift
	    tmd_plum_out="$1"
            ;;
        -mdp)
            shift
            mdp="$1"
            ;;
	-top)
            shift
            top="$1"
            ;;
	-n)
	    shift
	    index="$1"
	    ;;
	-ts)
	    shift
	    ts="$1" ;;
        --ttot|-t)
            shift
            ttot="$1"
	    tstep=$(echo "scale=0; ${ttot}/${ts}" | bc -l)
            ;;
	--kappa|-k)
	    shift
	    kappa="$1"
	    ;;
	-scale) shift
	    scale=$1;;
	-lbd) shift
	    LBD=$1;;
	-weights) shift
	    WEIGHTS=$1;;
	-Tmin) shift
	    TMIN=$1;;
	-Tmax) shift
	    TMAX=$1;;
	#because we could also use this script for annealing or equil
	-tempering)
	    tempering=1;;
	-wrex)
	    wrex=1;;
    esac
    shift
done


for i in E*; do
	if (($tempering==1)); then
		sed "s/KSI/${i:2}/g;s/NSTEPS/${tstep}/g;s/_ST_/${ts}/g;s/LBD/$LBD/g;s/WEIGHTS/$WEIGHTS/g;s/TMIN/$TMIN/g;s/TMAX/$TMAX/g" $mdp > $i/md.mdp
	elif (($wrex==1)); then
		:
	else
		sed "s/KSI/${i:2}/g;s/NSTEPS/${tstep}/g;s/_ST_/${ts}/g" $mdp > $i/md.mdp
	fi
done

ksis=$(command ls -d E*)

ksis=($ksis)
l=${#ksis[@]}
for k in $(seq 0 10 $l); do
    for ksi in "${ksis[@]:$k:10}"; do
        echo "ksi is $ksi"
	if (($wrex==1)); then
	    (gmx grompp -nice 0 -f $mdp -c "./equil/${ksi}/confout.gro" -p "$top" -o "./${ksi}/md.tpr" -maxwarn 1 -n "$index" &> ./${ksi}/grompp.out) &
	else
            mdp="$ksi/md.mdp"
            (gmx grompp -nice 0 -f $mdp -c "./${ksi}/equil.gro" -p "$top" -o "./${ksi}/md.tpr" -maxwarn 1 -n "$index" -r "./${ksi}/equil.gro" &> ./${ksi}/grompp.out) &
            touch "./${ksi}/plumed_${ksi:2}.dat"
	fi
    done
    wait
done
