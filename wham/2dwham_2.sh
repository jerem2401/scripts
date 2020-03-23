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
tmp=$(<$file1 wc -l)
collen=$(($tmp - 6))
chlen=$(($collen / $c))

echo $collen $chlen

for i in ${files}; do
    for k in `seq 0 $c`; do
        var2=$(basename "$i" .txt);
        echo "var2=${var2}";
        echo "k=${k}"i
	echo
        #awk '(NR>5) {print $1,$7,$6}' $i > ./$var2;
    done
done
#for i in ./colv*; do cnt=$(echo "$i" | grep -oP '(?<=colvar).*(?=_)'); cnt2=$(echo "$i" | grep -oP '(?<=_).*(?=.txt)'); echo "$i   $cnt   $cnt2   5000   10000" >> metd.txt; done
#for i in ./colv*; do cnt=$(echo "$i" | grep -oP '(?<=colvar).*(?=_)'); cnt2=$(echo "$i" | grep -oP '(?<=_).*(?=.txt)'); read cntkappa cnt2kappa <<< $(grep -oP '(?<=KAPPA=).*(?= )' ../../../E_$cnt/plumed_* ); echo "$i   $cnt   $cnt2   $cntkappa   $cnt2kappa" >> metd.txt; done
#column -t metd.txt > temp.txt
#sed '/inf/d' 2dpmf.txt > 2dpmf_clean.txt

