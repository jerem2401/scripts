#!/bin/bash

colf() {
	for i in 
}


while [ $# -gt 0 ]; do
    case "$1" in
	--file|-f)
            shift
	    tmd_plum_out="$1"
            ;;
        -h)
	    echo "easy peasy args detection with following dir format: N*nbOfWin*_t*ttot*_k*kappa*_gk*gkappa*_a*alpha*_deq*deq*"
            echo "-f tmd.txt -tpr tmd.tpr -traj tmd.traj -mdp temp.mdp -top mol.top -ref1 refpdb1.pdb -ref2 refpdb2.pdb -pltmp plumed_tmp.dat"
	    echo "other options: -w=nbOfWin, -wmax=winmax, -wmin=winmin, -refmid=ref pdb for gCV"
	    exit
            ;;
        -tpr)


#for i in ../../E_*/colv*; do var2=$(basename "$i"); awk '(NR>5) {print $1,$7,$6}' $i > ./$var2; done
#for i in ./colv*; do cnt=$(echo "$i" | grep -oP '(?<=colvar).*(?=_)'); cnt2=$(echo "$i" | grep -oP '(?<=_).*(?=.txt)'); echo "$i   $cnt   $cnt2   5000   10000" >> metd.txt; done


if [ $col -eq 1 ]; then
   for i in
