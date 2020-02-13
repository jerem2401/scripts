#!/bin/bash
echo $(date)
echo '!!!!!!!!!!!!!!!!!!! CAREFUL: k=8000 and rew=True, think about argument parsing to make this script versatile !!!!!!!!!!!!!!!!'
dir=$1
for i in $dir
do
	dirl+=( "$i" )
done

slice=$2
len=$(echo "${#dirl[@]}")
echo "total files: $len"
echo "slice chosen : $slice"

for i in $(seq 0 "$slice" "$len")
do
	for file in "${dirl[@]:${i}:${slice}}"
	do
		echo "will execute rew on ${file} with k =8000 and rew=True" &
		$(rew_hist_final.py -f $file -k 8000 -rew) &
		PID="$!"
		#PID_LIST+="$PID "
	done
	wait $PID
done

echo $(date)
