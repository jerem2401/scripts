#!/bin/bash

dir=$(d=$(pwd); basename $d)

nbOfWin=$(echo "$dir" | grep -oP '(?<=N).*(?=_t.*_k)')

ttot=$(echo "$dir" | grep -oP '(?<=_t).*(?=_k)')
tstep=$(echo "scale=0; $ttot/0.000004" | bc -l) #assuming that the time step is 0.004ps, add timestep as arguments ?

kappa=$(echo "$dir" | grep -oP '(?<=_k).*(?=_rep)')

scale=5

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
        --ttot|-t)
            shift
            ttot="$1"
	    tstep=$(echo "scale=0; $ttot/0.000002" | bc -l)
            ;;
	--kappa|-k)
	    shift
	    kappa="$1"
	    ;;
	-scale) shift
	    scale=$1;;
    esac
    shift
done

for i in E*; do
	sed "s/KSI/${i:2}/g;s/NSTEPS/${tstep}/g" $mdp > $i/md.mdp
done

ksis=$(command ls -d E*)

ksis=($ksis)
l=${#ksis[@]}
for k in $(seq 0 10 $l); do
    for ksi in "${ksis[@]:$k:10}"; do
        echo "ksi is $ksi"
	mdp="$ksi/md.mdp"
        (gmx grompp -nice 0 -f $mdp -c "./${ksi}/equil.gro" -p "$top" -o "./${ksi}/md.tpr" -maxwarn 1 -n "$index" -r "./${ksi}/equil.gro" &> ./${ksi}/grompp.out) &
	touch "./${ksi}/plumed_${ksi:2}.dat"
    done
    wait
done
