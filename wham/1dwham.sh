prep2dwham() {

    #check if we need to clean the data according to dir name (if clean > 0 => cure else not)
    clean=$(basename $(pwd) | grep -o '[0-9]*')

    if (("$clean" == 0)); then

        echo "careful, prep2dwham is aimed to work in a specific directory tree fromat"

        unset coll
        declare -a coll

        unset PID
        declare -a PID

	unset tosed
	declare -a tosed

	#for i in "${rpath}/E*/colv*"; do
	#echo "first element of rpath is: ${rpath[0]} ###################"
	#for i in "$rpath"; do
        #    coll+=( "$i" )
        #done
	read -r -a coll <<< "$rpath"
        colllen=$(echo "${#coll[@]}")

	echo "timeinit is $timeinit, timefinal is $timefinal"
        time5ns=$(awk -v var=$timeinit '$1~"^"var".*" {print NR;exit}' "${coll[0]}")

        echo "timeinit line = $time5ns"
	echo "column chosen: $column"

        for index in $(seq 0 8 "$colllen"); do
            for i in "${coll[@]:${index}:8}"; do
                var2=$(basename "$i")
                if [ -e "$var2" ]; then
                    echo "colvar file $var2 exists, skipping this file"
                    #just to play around with continue, could be replace by a else
                    continue
                fi
        	echo "preparing colvar: $var2"
                (awk -v var=$time5ns -v col=$column '(NR>=var) {print $1,$col}' $i > ./$var2 \
        	 && TMPFILE=$(mktemp ./foo-XXXXX) \
        	 && awk '!seen[$1]++' $var2 > $TMPFILE \
        	 && mv $TMPFILE $var2 \
		 && echo "timefinal line = $timefinal" \
		 && timefinal=$(awk -v var=$timefinal '$1~"^"var".*" {print NR;exit}' $var2) \
		 && TMPFILE2=$(mktemp ./foo-XXXXX) \
		 && awk -v var=$timefinal '(NR<=var)' $var2 > $TMPFILE2 \
		 && mv $TMPFILE2 $var2 \
		 && TMPFILE3=$(mktemp ./foo-XXXXX) \
        	 && grep -vi '[a-z]' $var2 | grep -vi '\#' > $TMPFILE3 \
		 && mv $TMPFILE3 $var2 \
		 && echo "removing [a-z]|\# from $var2") &
                PID+=( "$!" )
            done
            for pid in ${PID[*]}; do
                wait $pid
            done
        done

	if (($st==1)); then
	    for index in $(seq 0 8 "$colllen"); do
	        for i in "${coll[@]:${index}:8}"; do
		    var2=$(basename "$i")
		    echo "preparing colvar for umbST: $var2"
		    tmp=$(dirname $i)
		    cnt=$(echo "$i" | grep -oP '(?<=colvar_).*(?=.txt)')
		    (awk '($0 !~ /^#.*/) && ($0 !~ /^@.*/)' $tmp/dhdl.xvg > ./dhdl_${cnt}.txt && awk 'NR % 2 == 1' ./dhdl_${cnt}.txt > ./dhdl_${cnt}_every2nd.txt && awk '($4 == '0.0000000') {print NR}' ./dhdl_${cnt}_every2nd.txt > linenumber_${cnt}.txt && awk 'NR==FNR{linesToPrint[$0];next} FNR in linesToPrint' linenumber_${cnt}.txt ./$var2 > out_${cnt}.txt && mv out_${cnt}.txt ./$var2) &
		    PID+=( "$!" )
		done
		for pid in ${PID[*]}; do
		    wait $pid
		done
	    done
	    mkdir cleaning_st
	    mv dhdl* linenumber_* cleaning_st
	fi

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

   for i in ../colv*; do
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
      #for i in c$c/c_$k/colv*; do
      #   cnt=$(echo "$i" | grep -oP '(?<=colvar_).*(?=_)')
      #   #cnt=$(echo "$i" | grep -oP '(?<=colvar_).*(?=\.0_[0-9]\.txt)')
      #   #read cntkappa cnt2kappa <<< $(grep -oPh '(?<=KAPPA=)[0-9]*' ../../E_$cnt/plumed_* | tail -1)
      #   ########################
      #   #read cntkappa cnt2kappa <<< $(grep -oPh '(?<=KAPPA=)[0-9]*' ${rpath}/E_$cnt/plumed_* | head -1)
      #   read cntkappa cnt2kappa <<< $(grep -oPh '(?<=KAPPA=)[0-9]*' ${rpath}/plumed_files/E_$cnt/plumed_* | head -1)
      #   echo "$(basename $i)   $cnt   $cntkappa" >> "c${c}/c_${k}/metd.txt"
      #done
      read -r -a coll <<< "$rpath"
      for i in c$c/c_$k; do
          for i in "${coll[@]}"; do
              cnt=$(echo "$i" | grep -oP '(?<=colvar_).*(?=\.[0-9]*.txt)')
              tmp=$(dirname $i)
              tmp2=$(basename $i)
	      plumed="$tmp/plumed.dat"
	      read cntkappa cnt2kappa <<< $(grep -oPh '(?<=KAPPA=)[0-9]*' $plumed | head -1)
              echo "${tmp2%.txt}_${k}.txt   $cnt   $cntkappa" >> "c${c}/c_${k}/metd.txt"
          done
      done
   done
   echo "mk_metd done"
}

#args needed: -c, E1min,E1max,bin,E2min,E2max
do_wham2d() {
   echo "starting do_wham2d"
   if [[ $(hostname) = gwdu* ]]; then
      wham='/usr/users/jlapier/bin/wham/wham/wham'
   elif [[ $(hostname) == fullmetal ]]; then
      wham='/home/jeremy/opt/wham/wham/wham'
   elif [[ $(hostname) == smaug ]]; then
      wham='/home/users/jeremy/bin/wham/wham/wham'
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
