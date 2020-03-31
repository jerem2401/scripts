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

#args needed: -f, -c
mk_chunk() {
   seq -f "c_%.0f" 1 $c | xargs mkdir

   for i in ${files}; do
       collen=$(wc -l < $i)
       chlen=$(($collen / $c))
       echo "colvar len = ${collen} and chunk lenght = ${chlen}"
       startl=1
       endl=$chlen
       for k in `seq 1 $c`; do
           file=$(basename "$i" .txt);
   	   echo "file=${file}";
           echo "k=${k}";
   	   echo "startl = $startl";
	   if [ $k == $c ]; then
              endl=$collen
           else
   	      endl=$((endl+chlen))
           fi
   	   echo "endl   = $endl";
   	   $(awk -v a=$startl -v b=$endl '(NR>=a) && (NR<=b)' $i > "c_${k}/${file}_${k}.txt") &
   	   startl=$((startl+chlen));
   	   PID="$!"
       done
       wait $PID
   done
}

#args needed: -c
mk_metd() {
   for k in `seq 1 $c`; do
      echo "making metd.txt for c_${k}"
      for i in c_$k/colv*; do
         cnt=$(echo "$i" | grep -oP '(?<=colvar)[^_]*')
         cnt2="$(cut -d'_' -f3 <<<$i)"
	 read cntkappa cnt2kappa <<< $(grep -oPh '(?<=KAPPA=).*(?= )' ../../../../E_$cnt/plumed_* | tr "\n" " ")
         echo "$(basename $i)   $cnt   $cnt2   $cntkappa   $cnt2kappa" >> "c_${k}/metd.txt"
      done
   done
}

#args needed: -c, E1min,E1max,bin,E2min,E2max
do_wham2d() {
   for i in `seq 0 4 $c`; do
      for k in `seq $(($i+1)) $(($i+4))`; do
         if cd c_$k; then
            echo "launching 2dwham in c_$k with: E1min=$E1min; E1max=$E1max; bin=$bin; E2min=$E2min; E2max=$E2max"
            $(~/bin/wham/wham-2d/wham-2d Px=0 $E1min $E1max $bin Py=0 $E2min $E2max $bin 0.000001 300 0 metd.txt 2dpmf.txt 0) &
            PID="$!"
            cd ..
         fi
      done
      wait $PID
   done
}

clean_wham2d() {
    for k in `seq 1 $c`; do
        sed '/inf/d' "c_${k}/2dpmf.txt" > "c_${k}/2dpmf_clean.txt"
        $(projection.py -f "c_${k}/2dpmf_clean.txt" -o "c_${k}/1dpmf_${k}.txt")
    done
}

#mk_chunk
#mk_metd

#column -t metd.txt > temp.txt
#sed '/inf/d' 2dpmf.txt > 2dpmf_clean.txt

