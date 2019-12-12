#!/bin/bash

dir=$(pwd)

nbOfWin=$(echo "$dir" | grep -oP '(?<=_N).*(?=_t.*_k)')
winmax=1
winmin=-1

ttot=$(echo "$dir" | grep -oP '(?<=_t).*(?=_k)')
tstep=$(echo "scale=0; $ttot/0.000002" | bc -l) #assuming that the time step is 0.002ps, add timestep as arguments ?

alpha=$(echo "$dir" | grep -oP '(?<=_a).*(?=_deq)')

deq=$(echo "$dir" | grep -oP '(?<=_deq)([0-9]|\.)*')

kappa=$(echo "$dir" | grep -oP '(?<=_k).*(?=_gk)')
gkappa=$(echo "$dir" | grep -oP '(?<=_gk).*(?=_a)')

#tmd_tpr='md.tpr'
#tmd_traj='traj_comp.xtc'
temp_mdp='../temp/temp.mdp'
top='../temp/c7eq.top'
refpdb1='../temp/c7eqf.pdb'
refpdb2='../temp/c7axf.pdb'
plumed_tmp='../temp/plumed_tmp.dat'

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
            shift
            tmd_tpr="$1"
            ;;
        -traj)
            shift
            tmd_traj="$1"
            ;;
        -mdp)
            shift
            temp_mdp="$1"
            ;;
	-top)
            shift
            top="$1"
            ;;
	-ref1)
            shift
            refpdb1="$1"
	    ;;
	-ref2)
            shift
            refpdb2="$1"
            ;;
	-refmid)
	    shift
	    refmid="$1"
	    ;;
	-pltmp)
            shift
            plumed_tmp="$1"
            ;;
        --win|-w)
            shift
            nbOfWin="$1"
	    ;;
	--winmax|-wmax)
	    shift
	    winmax="$1"
	    ;;
        --winmin|-wmin)
            shift
            winmin="$1"
            ;;
        --ttot|-t)
            shift
            ttot="$1"
	    tstep=$(echo "scale=0; $ttot/0.000002" | bc -l)
            ;;
        --alpha|-a)
            shift
            alpha="$1"
            ;;
	--deq|-d)
	    shift
	    deq="$1"
	    ;;
	--kappa|-k)
	    shift
	    kappa="$1"
	    ;;
	-ext)
	    shift
	    tmd_tpr='../tmd.tpr'
	    tmd_traj='../tmd.trr'
	    temp_mdp='../temp.mdp'
	    top='../c7eq.top'
	    refpdb1='../genc7eq_2.pdb'
	    refpdb2='../genc7ax_3.pdb'
	    plumed_tmp='../plumed_temp.dat'
	    ;;
    esac
    shift
done

max_min=$(echo "$winmax- $winmin" | bc)
dist="${max_min#-}"
winStep=$(echo "$dist/($nbOfWin- 1)" | bc -l)

sed "s/NSTEPS/${tstep}/" "$temp_mdp" > ./md.mdp

echo "arg1:${tmd_plum_out} nbOfWin:${nbOfWin} wmax:${winmax} wmin:${winmin} dist:${dist} winStep:${winStep} ttot${ttot} tstep:${tstep} alpha:${alpha} deq:${deq} kappa:${kappa}"

echo -n 'on which group do you want to extract frames ?'
read group

for ksi in $(seq -f "%.2f" "$winmin" "$winStep" "$winmax"); do
   if mkdir "E_${ksi}"; then
	values=$(start4umb.py -f $tmd_plum_out -v $ksi)
	value=$(echo "$values" | grep -oP '(?<=\().*(?=,)')
	gksi=$(echo "$values" | grep -oP '(?<=, ).*(?=\))')
        gmx trjconv -s "$tmd_tpr" -f "$tmd_traj" -dump "$value" -o "./E_${ksi}/conf_${value}.gro" <<EOF
	$group
EOF
	gmx grompp -f md.mdp -c "./E_${ksi}/conf_${value}.gro" -p "$top" -o "./E_${ksi}/conf_${ksi}.tpr" -maxwarn 1	
	cp "$refpdb1" "./E_${ksi}"
	cp "$refpdb2" "./E_${ksi}"
	cp "$refmid" "./E_${ksi}"
	sed "s/_POSI_/${ksi}/g;s/_gPOSI_/${gksi}/g;s/_DEQ_/${deq}/g;s/_ALPHA_/${alpha}/g;s/_KAPPA_/${kappa}/g;s/_gKAPPA_/${gkappa}/g;" "$plumed_tmp" > "./E_${ksi}/plumed_${ksi}.dat"
	echo -e "\nused tmd directory: ${tmd_plum_out}" >> mdout.mdp
	echo 'done'
   else
	echo "window already exists"
   fi
done


