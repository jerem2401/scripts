#!/bin/bash

module load anaconda3/2020.07 && source activate env1
dir=$(d=$(pwd); basename $d)

nbOfWin=$(echo "$dir" | grep -oP '(?<=N).*(?=_t.*_k)')
#winmax=1
#winmin=-1

ttot=$(echo "$dir" | grep -oP '(?<=_t).*(?=_k)')
tstep=$(echo "scale=0; $ttot/0.000004" | bc -l) #assuming that the time step is 0.004ps, add timestep as arguments ?

kappa=$(echo "$dir" | grep -oP '(?<=_k).*(?=_rep)')

#tmd_tpr='md.tpr'
#tmd_traj='traj_comp.xtc'
temp_mdp='../temp/temp.mdp'
top=$(echo ../temp/*.top)
plumed_tmp='../temp/plumed_tmp.dat'
index='../temp/index.ndx'

lbd=''
scale=2

while [ $# -gt 0 ]; do
    case "$1" in
	--file|-f)
            shift
	    tmd_plum_out="$1"
            ;;
        -h)
	    echo "easy peasy args detection with following dir format: N*nbOfWin*_t*ttot*_k*kappa*_wk*wkappa*_wcnt*wallcnt*_wof*walloffset*"
            echo "-f tmd.txt -tpr tmd.tpr -traj tmd.traj -mdp temp.mdp -top mol.top -pltmp plumed_tmp.dat"
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
	-pltmp)
            shift
            plumed_tmp="$1"
            ;;
	-n)
	    shift
	    index="$1"
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
	    tstep=$(echo "scale=0; $ttot/0.000004" | bc -l)
            ;;
	--kappa|-k)
	    shift
	    kappa="$1"
	    ;;
	--wkappa|-wk)
	    shift
	    wkappa="$1"
	    ;;
	--wcnt)
	    shift
	    wcnt="$1"
	    ;;
	-path)
	    shift
	    path=$1
	    ;;
	-lam)
	    shift
	    lam=$1
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
	-hrex)
	    shift
	    lbd=$1
	    ;;
	-listwin)
	    shift
	    listwin=$1;;
	-scale) shift
	    scale=$1;;
    esac
    shift
done

if [[ ! -z $lbd ]]; then
	for i in $lbd; do
		if [[ ! -f ../temp/topol_${i}.top ]]; then
			echo "../temp/topol_${i}.top not present, exiting"
			exit
		fi
	done
fi


max_min=$(echo "$winmax- $winmin" | bc)
dist="${max_min#-}"
winStep=$(echo "scale=$scale; $dist/($nbOfWin- 1)" | bc -l)

sed "s/NSTEPS/${tstep}/" "$temp_mdp" > ./md.mdp

echo "arg1:${tmd_plum_out} nbOfWin:${nbOfWin} wmax:${winmax} wmin:${winmin} dist:${dist} winStep:${winStep} ttot:${ttot} tstep:${tstep} kappa:${kappa} wkappa:${wkappa} wcnt=${wcnt} wof=${wof} path: ${path} lambda: ${lam}"


if [[ ! -z $listwin ]]; then
	ksis=$(echo $listwin | sed 's/tmp_//g')
	echo -e "\nksis are:\n$ksis"
else
	ksis=$(seq -f "%.3f" "$winmin" "$winStep" "$winmax")
fi

ksis=($ksis)
l=${#ksis[@]}
for k in $(seq 0 10 $l); do
    for ksi in "${ksis[@]:$k:10}"; do
        echo "ksisi $ksi"
        tmp=$(find ./bench* -type d -name "E_${ksi}")
        echo "tmp $tmp"
        if [[ ! -z ${tmp} ]]; then
    	copy=$(echo $tmp | awk '{print $1;}')
            cp -r ${copy} ./
    	echo "${ksi} in ${tmp}, copying."
        elif mkdir "E_${ksi}"; then
            values=$(start4umb_dz.py -f $tmd_plum_out -v $ksi)
            read -ra ADDR <<< $(echo "${values//[\(\)\,]}")
            closestksi=${ADDR[0]}
            value=${ADDR[1]}
            echo "given: $ksi, found: $closestksi"
            b=$(echo "$value - 2" | bc -l)
            echo 0 | gmx trjconv -nice 0 -s "$tmd_tpr" -f "$tmd_traj" -b "$b" -dump "$value" -o "./E_${ksi}/conf_${value}.gro" -n "$index"
            if [[ ! -z $lbd ]]; then
                for i in $lbd; do
                    mkdir "E_${ksi}/${i}"
                    gmx grompp -nice 0 -f ./md.mdp -c "./E_${ksi}/conf_${value}.gro" \
                    -p ../temp/topol_${i}.top -o "./E_${ksi}/${i}/md.tpr" \
                    -maxwarn 1 -n "$index"
                done
    	    index_in_plumed='../index.ndx'
            else
                (gmx grompp -nice 0 -f ./md.mdp -c "./E_${ksi}/conf_${value}.gro" -p "$top" -o "./E_${ksi}/conf_${ksi}.tpr" -maxwarn 1 -n "$index") &
    	    index_in_plumed='index.ndx'
    	fi
    
            sed "s=_POSI_=${ksi}=g;s=_KAPPA_=${kappa}=g;s=_wKAPPA_=${wkappa}=g;s=index.ndx=${index_in_plumed}=g" "$plumed_tmp" > "./E_${ksi}/plumed_${ksi}.dat"
            cp "$index" "./E_${ksi}/"
            echo 'done'
        else
            echo "window already exists"
        fi
    done
    wait
done
