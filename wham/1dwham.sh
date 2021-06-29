prep2dwham() {

    #check if we need to clean the data according to dir name (if clean > 0 => cure else not)
    clean=$(basename $(pwd) | grep -o '[0-9]*')

    if (("$clean" == 0)); then

        echo "careful, prep2dwham is aimed to work in a specific directory tree fromat"

        unset coll
        declare -a coll

        unset PID
        declare -a PID

        for i in ../../E*/colv*; do
            coll+=( "$i" )
        done
        colllen=$(echo "${#coll[@]}")

        time5ns=$(awk '$1~/^5000.*/{print NR;exit}' "${coll[0]}")
        echo "time5ns = $time5ns"
	echo "column chosen: $column"

        for index in $(seq 0 4 "$colllen"); do
            for i in "${coll[@]:${index}:4}"; do
                var2=$(basename "$i")
                if [ -e "$var2" ]; then
                    echo "colvar file $var2 exists, skipping this file"
                    #just to play around with continue, could be replace by a else
                    continue
                fi
        	echo "preparing colvar: $var2"
                (awk -v var=$time5ns -v col=$column '(NR>=var) && ($0 !~ /^#.*/) {print $1,$col}' $i > ./$var2) &
                PID+=( "$!" )
            done
            for pid in ${PID[*]}; do
                wait $pid
            done
        done
        echo "prep2dwham done"

    else
        echo "option not supported for 1dwham"
        exit
    fi
}


mk_chunk() {
   echo "starting mk_chunk"
   mkdir "c${c}"
   cd "c${c}"
   seq -f "c_%.0f" 1 $c | xargs mkdir

   for i in ${files}; do
       collen=$(wc -l < $i)
       chlen=$(($collen / $c))
       echo "colvar len = ${collen} and chunk lenght = ${chlen}"
       startl=1
       endl=$chlen

       unset PID
       declare -a PID

       for k in `seq 1 $c`; do
           file=$(basename "$i" .txt)
   	   echo "file=${file}"
           echo "k=${k}"
   	   echo "startl = $startl"

	   if [ $k == $c ]; then
              endl=$collen
           elif [ $k == 1 ]; then
              endl=$chlen
           else
   	      endl=$((endl+chlen))
           fi

   	   echo "endl   = $endl"
   	   $(awk -v a=$startl -v b=$endl '(NR>=a) && (NR<=b)' $i > "c_${k}/${file}_${k}.txt") &
   	   startl=$((startl+chlen))
   	   PID+=( "$!" )
       done
       for pid in ${PID[*]}; do
           wait $pid
       done
   done
   echo "mk_chunk done"
}

#args needed: -c
mk_metd() {
   for k in `seq 1 $c`; do
      [ -e "c${c}/c_${k}/metd.txt" ] \
      && echo "c${c}/c_${k}/metd.txt exists, removing file first" \
      && rm -f "c${c}/c_${k}/metd.txt"

      echo "making metd.txt for c_${k}"
      for i in c$c/c_$k/colv*; do
         cnt=$(echo "$i" | grep -oP '(?<=colvar_).*(?=_)')
	 cntkappa=$(grep -oPh '(?<=KAPPA=)[0-9]*' ../../E_$cnt/plumed_* | tail -1)
         echo "$(basename $i)   $cnt   $cntkappa" >> "c${c}/c_${k}/metd.txt"
      done
   done
   echo "mk_metd done"
}

#args needed: -c, E1min,E1max,bin,E2min,E2max
do_wham2d() {
   echo "starting do_wham2d"
   if [ $(hostname) = dema9 ]; then
      wham='/localdisk/bin/wham/wham-2d/wham-2d'
   elif [[ $(hostname) == gwdu* ]]; then
      wham='/usr/users/jlapier/bin/wham/wham/wham'
   elif [[ $(hostname) == fullmetal ]]; then
      wham='/home/jeremy/opt/wham/wham/wham'
   fi

   for i in `seq 0 4 $c`; do
      unset PID
      declare -a PID
      for k in `seq $(($i+1)) $(($i+4))`; do
         if cd "c${c}/c_${k}"; then
            echo "launching 2dwham in c_$k with: E1min=$E1min; E1max=$E1max; bin=$bin"
            ($wham $E1min $E1max $bin 0.000001 300 0 metd.txt 1dpmf.txt &> 1dwham.out) &
            PID+=( "$!" )
            cd ../..
         fi
      done
      for pid in ${PID[*]}; do
          wait $pid
      done
   done
   echo "do_wham2d done"
}

clean_wham2d() {
   echo "starting clean_wham2d"
   for k in `seq 1 $c`; do
       awk '/^[^#]/ { print $1 " " $4  " " $2 }' "c${c}/c_${k}/1dpmf.txt" > "c${c}/c_${k}/tmp.txt"
       wait $!
       sed '/inf/d' "c${c}/c_${k}/tmp.txt" > "c${c}/c_${k}/1dpmf_clean_${k}.txt"
       wait $!
   done
   python -c "from block import block_avg; print(block_avg('c*/c_*/1dpmf_clean*'))"
   rm -f c*/c*/tmp.txt
   echo "clean_wham2d done"
}
