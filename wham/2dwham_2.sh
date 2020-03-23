#!/usr/bin/bash

while [ $# -gt 0 ]; do
    case "$1" in
        -f)
            shift
            files="$1"
            ;;
        -c)
            shift
            c="$1"
	    ;;
        -h)
            echo "easy peasy args detection with following dir format: N*nbOfWin*_t*ttot*_k*kappa*_gk*gkappa*_a*alpha*_deq*deq*"
            ;;
    esac
    shift
done

file1=$(echo $files | grep -o '^\S*')
collen=$(<$file1 wc -l)
chlen=$(($collen / $c))

echo $collen $chlen

seq -f "c_%.0f" 1 $c | xargs mkdir

for i in ${files}; do
    startl=1
    endl=$chlen
    for k in `seq 1 $c`; do
        file=$(basename "$i" .txt);
	echo "file=${file}";
        echo "k=${k}";
	echo "startl = $startl";
	echo "endl   = $endl";
	$(awk -v a=$startl -v b=$endl '(NR>=a) && (NR<=b)' $i > "c_${k}/${file}_${k}.txt") &
	PID="$!"
	startl=$((startl+chlen));
	endl=$((endl+chlen));
    done
    wait $PID
done
#for i in ./colv*; do cnt=$(echo "$i" | grep -oP '(?<=colvar).*(?=_)'); cnt2=$(echo "$i" | grep -oP '(?<=_).*(?=.txt)'); echo "$i   $cnt   $cnt2   5000   10000" >> metd.txt; done
#for i in ./colv*; do cnt=$(echo "$i" | grep -oP '(?<=colvar).*(?=_)'); cnt2=$(echo "$i" | grep -oP '(?<=_).*(?=.txt)'); read cntkappa cnt2kappa <<< $(grep -oP '(?<=KAPPA=).*(?= )' ../../../E_$cnt/plumed_* ); echo "$i   $cnt   $cnt2   $cntkappa   $cnt2kappa" >> metd.txt; done
#column -t metd.txt > temp.txt
#sed '/inf/d' 2dpmf.txt > 2dpmf_clean.txt

