#!/bin/bash

while read line; do 
	E=$(echo $line | grep -oP '(?<=histo_).*(?=_)'); 
	echo E_$E;
	v=$(echo $line | awk '{print $4}'); 
	echo $v;
	f=$(echo $line | grep -oP '(?<=histo_).*(?= )'); 
	echo "colvar${f}.txt";
	Elist+=( "E_${E}" )
	vlist+=( "${v}" )
	flist+=( "colvar${f}" )
done < ../histo_RMSDMID/hist_curing.txt

echo "var attribution done"

slice=$1
len=$(echo "${#flist[@]}")

#lenT=$(echo "$((${len}-1))")

PWD=$(pwd)

for i in $(seq 0 "$slice" "$len")
do
	echo "i= ${i}"
	endi=$(echo "$((${i}+${slice}))")
	echo "endi= ${endi}"
	ib=$i
	for file in "${flist[@]:${i}:${slice}}"
	do
		if [ ! -f $file ]
		then
			v="${vlist[${ib}]}"
			E="${Elist[${ib}]}"
			echo "will execute awk on ${file} with v= ${v} and dir = ${E}" &
			awk -v v=$v 'NR<=5{print}; (NR>5) && ($6 < v) {print}' "../../${E}/${file}" > "./${file}" &
			PID="$!"
			ib=$(echo "$((${ib}+1))")
			#PID_LIST+="$PID "
		else
			echo "${file} already exists"
		fi
	done
	wait $PID
done
