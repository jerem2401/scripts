#! /bin/bash

wat=0

while [ $# -gt 0 ]; do
    case "$1" in
	-wat)
            shift
	    wat=1
            ;;
    esac
    shift
done

if [ "$wat" -eq 1 ]; then
    c=2
    for i in $(command ls -d ../../E_*); do
        sub=${i#../../E_}
        if (( $(echo "$sub > 0.9" | bc -l) )) && (( $c%2 == 0 )); then
            echo $sub
            echo 0 | gmx trjconv -f "${i}/traj_comp.xtc" -s $i/*.tpr -o "nopbc1_${sub}_wat.xtc" -pbc atom -ur compact -dt 200
            wait $!
            echo 1 0 | gmx trjconv -f "nopbc1_${sub}_wat.xtc" -s $i/*.tpr -o "nopbc2_${sub}_wat.xtc" -center -pbc mol -ur compact
            #echo 0 | gmx trjconv -f "nopbc2_${sub}.xtc" -s $i/*.tpr -o "nopbc3_${sub}.xtc" -pbc whole
            lsub=$sub
            li=$i
        fi
    #let c=c+1
    done
    #gmx trjcat -f nopbc2_* -o cat_wat.xtc -cat
    echo 0 | gmx trjconv -f "nopbc2_${lsub}_wat.xtc" -s $li/*.tpr -o 0_wat.pdb -dump 0
    wait $!
    #rm -f ./nopbc*

else
    c=2
    for i in $(command ls -d ../../E_*); do
        sub=${i#../../E_}
	if (( $(echo "$sub > -0.25" | bc -l) )); then
        #if (( $(echo "$sub > 0.8" | bc -l) )) && (( $(echo "$sub < 1.000" | bc -l) )) && (( $c%2 == 0 )); then
            echo $sub
            echo 1 | gmx trjconv -f "${i}/traj_comp.xtc" -s $i/*.tpr -o "nopbc1_${sub}.xtc" -pbc atom -ur compact -dt 2000
            wait $!
            #echo 1 | gmx trjconv -f "nopbc1_${sub}.xtc" -s $i/*.tpr -o "nopbc2_${sub}.xtc" -pbc nojump
            echo 1 | gmx trjconv -f "nopbc1_${sub}.xtc" -s $i/*.tpr -o "nopbc2_${sub}.xtc" -pbc mol -ur compact
            lsub=$sub
            li=$i
        fi
    let c=c+1
    done
    gmx trjcat -f $(command ls nopbc2_* | sort -t '_' -k2 -n) -o cat_prot.xtc -cat
    echo 1 | gmx trjconv -f "nopbc2_${lsub}.xtc" -s $li/*.tpr -o 0_prot.pdb -dump 0
    wait $!
    rm -f ./nopbc*
fi
