#!/bin/bash
echo $(date)

dir=$1
slice=$2
cmd=$3

for i in $dir
do
	dirl+=( "$i" )
done

len=$(echo "${#dirl[@]}")
echo "total files: $len"
echo "slice chosen : $slice"

for i in $(seq 0 "$slice" "$len")
do
	for file in "${dirl[@]:${i}:${slice}}"
	do
		name=$(echo ${file} | grep -oP '(?<=colvar).*(?=.txt)')
		echo "will execute rew on ${file} with command: ${cmd} and output as histo_${name}.txt"
		$(rew_hist_final.py -f $file -k 20000 -min 0 -max 0.18 -s 0.001 -col RMSDMID --o histo_$name.txt) &
		PID="$!"
		#PID_LIST+="$PID "
	done
	wait $PID
done

echo $(date)
