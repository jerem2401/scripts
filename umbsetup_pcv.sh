#!/bin/bash

dir=$(d=$(pwd); basename $d)

nbOfWin=$(echo "$dir" | grep -oP '(?<=N).*(?=_t.*_k)')
#winmax=1
#winmin=-1

ttot=$(echo "$dir" | grep -oP '(?<=_t).*(?=_k)')
tstep=$(echo "scale=0; $ttot/0.000004" | bc -l) #assuming that the time step is 0.004ps, add timestep as arguments ?

wcnt=$(echo "$dir" | grep -oP '(?<=_wcnt).*(?=_wof)')
wof=$(echo "$dir" | grep -oP '(?<=_wof).*')

kappa=$(echo "$dir" | grep -oP '(?<=_k).*(?=_wk)')
wkappa=$(echo "$dir" | grep -oP '(?<=_wk).*(?=_wcnt)')

#tmd_tpr='md.tpr'
#tmd_traj='traj_comp.xtc'
temp_mdp='../temp/temp.mdp'
top=$(echo ../temp/*.top)
plumed_tmp='../temp/plumed_tmp.dat'
index='../temp/index.ndx'

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
    esac
    shift
done

if [ "$path" = '' ]; then
	echo "you should give the path.pdb path, exiting"
	exit
else
	path=$(readlink -f $path)
fi

max_min=$(echo "$winmax- $winmin" | bc)
dist="${max_min#-}"
winStep=$(echo "$dist/($nbOfWin- 1)" | bc -l)

sed "s/NSTEPS/${tstep}/" "$temp_mdp" > ./md.mdp

echo "arg1:${tmd_plum_out} nbOfWin:${nbOfWin} wmax:${winmax} wmin:${winmin} dist:${dist} winStep:${winStep} ttot:${ttot} tstep:${tstep} kappa:${kappa} wkappa:${wkappa} wcnt=${wcnt} wof=${wof} path: ${path} lambda: ${lam}"

echo -n 'on which group do you want to extract frames ?'
read group

for ksi in $(seq -f "%.3f" "$winmin" "$winStep" "$winmax"); do
   if mkdir "E_${ksi}"; then
	values=$(start4umb_pcv.py -f $tmd_plum_out -v $ksi -wm $winmin $winmax)
	read -ra ADDR <<< $(echo "${values//[\(\)\,]}")
	closestksi=${ADDR[0]}
	value=${ADDR[1]}
	echo "given: $ksi, found: $closestksi"
	#b=$(echo "$value - 1" | bc -l)
        #gmx trjconv -s "$tmd_tpr" -f "$tmd_traj" -b "$b" -dump "$value" -o "./E_${ksi}/conf_${value}.gro" -n "$index" <<EOF
	#$group
#EOF
	#gmx grompp -f md.mdp -c "./E_${ksi}/conf_${value}.gro" -p "$top" -o "./E_${ksi}/conf_${ksi}.tpr" -maxwarn 1 -n "$index"
	#sed "s=_POSI_=${ksi}=g;s=_KAPPA_=${kappa}=g;s=_wKAPPA_=${wkappa}=g;s=_OFFSET_=${wof}=g;s=_wPOSI_=${wcnt}=g;s=_PATH_=${path}=g;s=_LAMBDA_=${lam}=g;" "$plumed_tmp" > "./E_${ksi}/plumed_${ksi}.dat"
	#echo -e "\nused tmd directory: ${tmd_plum_out}" >> mdout.mdp
	#echo 'done'
   else
	echo "window already exists"
   fi
done


