#!/usr/bin/bash

prep2dwham() {
   echo "careful, prep2dwham is aimed to work in a specific directory tree fromat"
   for i in ../colv*; do
      var2=$(basename "$i")
      awk '(NR>5) {print $1,$7,$6}' $i > ./$var2
   done
   echo "prep2dwham done"
}

mk_chunk() {
   echo "starting mk_chunk"
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
   echo "mk_chunk done"
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
   echo "mk_metd done"
}

#args needed: -c, E1min,E1max,bin,E2min,E2max
do_wham2d() {
   echo "starting do_wham2d"
   if [ $(hostname) = dema9 ]; then
      wham='/localdisk/bin/wham/wham-2d/wham-2d'
   else
      wham='~/bin/wham/wham-2d/wham-2d'
   fi

   for i in `seq 0 4 $c`; do
      for k in `seq $(($i+1)) $(($i+4))`; do
         if cd c_$k; then
            echo "launching 2dwham in c_$k with: E1min=$E1min; E1max=$E1max; bin=$bin; E2min=$E2min; E2max=$E2max"
            ($wham Px=0 $E1min $E1max $bin Py=0 $E2min $E2max $bin 0.000001 300 0 metd.txt 2dpmf.txt 0 &> 2dwham.out) &
            PID="$!"
            cd ..
         fi
      done
      wait $PID
   done
   sleep 10
   echo "do_wham2d done"
}

clean_wham2d() {
   echo "starting clean_wham2d"
   for k in `seq 1 $c`; do
       sed '/inf/d' "c_${k}/2dpmf.txt" > "c_${k}/2dpmf_clean.txt"
       $(projection.py -f "c_${k}/2dpmf_clean.txt" -o "c_${k}/1dpmf_${k}.txt")
   done
   python -c "from block import block_avg; print(block_avg('c_*/1dpmf*'))"
   echo "clean_wham2d done"
}

#mk_chunk
#mk_metd

#column -t metd.txt > temp.txt
#sed '/inf/d' 2dpmf.txt > 2dpmf_clean.txt

