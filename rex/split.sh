#!/bin/bash
###############################################################################
#Author: Jeremy Lapierre <jeremy.lapierre@uni-saarland.de>
#Option description:
#
###############################################################################
#set -o errexit   # abort on nonzero exitstatus
set -o nounset   # abort on unbound variable
set -o pipefail  # dont hide errors within pipes

ndir=$(command ls E* -d | sort -t _ -k 2 | wc -l)
ed=$(( $ndir / 14 ))
ed_1=$(($ed-1))
rem=$(( $ndir-(14*$ed) ))
echo $ed
echo $rem

if [[ $(basename $PWD) = pos ]]; then
	for i in $(seq 0 $ed_1); do
		mkdir $i
		for j in $(command ls E* -d | sort -t _ -k 2 | head -n 14 ); do
			mv $j $i
		done
	done
	tmpmid=$(command ls ../mid/E* -d | sort -t _ -k 2 -n | tail -n1)
	tmp=$(basename $tmpmid)
	mkdir ./0/TMP_${tmp:2}
	cp $tmpmid/{md.tpr,index.ndx,plumed*,conf_*.gro} ./0/TMP_${tmp:2}
	tmp=$(command ls 1/E* -d | sort -t _ -k 2 -n | head -n1)
	cp -r $tmp ./0/TMP_${tmp:4}

	mkdir $ed
	mv E* $ed
	for i in $(seq 1 $ed); do
		echo $i
		k=$(($i-1))
		k2=$(($i+1))
		tmp=$(command ls $k/E* -d | sort -t _ -k 2 -n | tail -n1)
		cp -r $tmp ./$i/TMP_${tmp:4}
		if (($i!=$ed)); then
			tmp=$(command ls $k2/E* -d | sort -t _ -k 2 -n | head -n1)
			cp -r $tmp ./$i/TMP_${tmp:4}
		else
			echo "entering dir 4"
			tmp=$((8-$rem-1))
			if (($tmp>0)); then
				tmp2=$(command ls $k/E* -d | sort -t _ -k 2 -n | head -n 13 | tail -n $tmp)
				echo $tmp2
				for l in $tmp2; do
					echo "cp -r $l ./$i/TMP_${tmp2:4}"
					cp -r $l ./$i/TMP_${l:4}
				done
			elif (($tmp<0)); then
				echo "$ed contains more than 8 and less than 14 ksis, have a look"
			fi
		fi
	done
elif [[ $(basename $PWD) = neg ]]; then
	for i in $(seq 0 $ed_1); do
		mkdir $i
		for j in $(command ls E* -d | sort -t _ -k 2 | head -n 14 ); do
			mv $j $i
		done
	done
	tmpmid=$(command ls ../mid/E* -d | sort -t _ -k 2 -n | head -n1)
	tmp=$(basename $tmpmid)
	mkdir ./0/TMP_${tmp:2}
	cp $tmpmid/{md.tpr,index.ndx,plumed*,conf_*.gro} ./0/TMP_${tmp:2}
	tmp=$(command ls 1/E* -d | sort -t _ -k 2 -n | tail -n1)
	cp -r $tmp ./0/TMP_${tmp:4}

	mkdir $ed
	mv E* $ed
	for i in $(seq 1 $ed); do
		echo $i
		k=$(($i-1))
		k2=$(($i+1))
		tmp=$(command ls $k/E* -d | sort -t _ -k 2 -n | head -n1)
		cp -r $tmp ./$i/TMP_${tmp:4}
		if (($i!=$ed)); then
			tmp=$(command ls $k2/E* -d | sort -t _ -k 2 -n | tail -n1)
			cp -r $tmp ./$i/TMP_${tmp:4}
		else
			tmp=$((8-$rem-1))
			if (($tmp>0)); then
				tmp2=$(command ls $k/E* -d | sort -t _ -k 2 -n | tail -n 13 | head -n $tmp)
				for l in $tmp2; do
					cp -r $l ./$i/TMP_${l:4}
				done
			elif (($tmp<0)); then
				echo "$ed contains more than 8 and less than 14 ksis, have a look"
			fi
		fi
	done
fi
