#!/usr/bin/bash

#prep2dwham() {
#   echo "careful, prep2dwham is aimed to work in a specific directory tree fromat"
#   for f in ./colv*; do
#      [ -e "$f" ] && echo "colvar files exist, exiting prep2dwham" && return || echo "colvar files do not exist, continuing prep2dwham"
#   break
#   done
#
#   for i in ../../../colv*; do
#      var2=$(basename "$i")
#      awk '(NR>5) {print $1,$7,$6}' $i > ./$var2
#   done
#   echo "prep2dwham done"
#}

#prep2dwham() {
#    echo "careful, prep2dwham is aimed to work in a specific directory tree fromat"
#    for i in ../../../E_*/colv*; do
#        var2=$(basename "$i")
#        if [ -e "$var2" ]; then
#            echo "colvar file $var2 exists, skipping this file"
#            continue
#        fi
#	echo "preparing colvar: $var2"
#        time5ns=$(awk '$1~/^5000.*/{print NR;exit}' $i)
#        awk -v var=$time5ns '(NR>=var) {print $1,$7,$6}' $i > ./$var2
#    done
#    echo "prep2dwham done"
#}

prep2dwham() {

    #check if we need to clean the data according to dir name (if clean > 0 => cure else not)
    clean=$(basename $(pwd) | grep -o '[0-9]*')

    if (("$clean" == 0)); then

        echo "careful, prep2dwham is aimed to work in a specific directory tree fromat"

        unset coll
        declare -a coll

        unset PID
        declare -a PID

        for i in ../../colv*; do
            coll+=( "$i" )
        done
        colllen=$(echo "${#coll[@]}")

        time5ns=$(awk '$1~/^5000.*/{print NR;exit}' "${coll[0]}")
        echo "time5ns = $time5ns"

        for index in $(seq 0 4 "$colllen"); do
            for i in "${coll[@]:${index}:4}"; do
                var2=$(basename "$i")
                if [ -e "$var2" ]; then
                    echo "colvar file $var2 exists, skipping this file"
                    #just to play around with continue, could be replace by a else
                    continue
                fi
        	echo "preparing colvar: $var2"
                (awk -v var=$time5ns '(NR>=var) {print $1,$7,$6}' $i > ./$var2) &
                PID+=( "$!" )
            done
            for pid in ${PID[*]}; do
                wait $pid
            done
        done
        echo "prep2dwham done"

     else
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

        time5ns=$(awk '$1~/^5000.*/{print NR;exit}' "../../${Elist[0]}/${flist[0]}")
        echo "time5ns = $time5ns"

        echo "var attribution done"

        #slice=$1
        len=$(echo "${#flist[@]}")

        #lenT=$(echo "$((${len}-1))")

        for i in $(seq 0 4 "$len"); do
        	echo "i= ${i}"
        	endi=$(echo "$((${i}+4))")
        	echo "endi= ${endi}"
        	ib=$i

                unset PID
                declare -a PID

        	for file in "${flist[@]:${i}:4}"; do
                    if [ ! -f $file ]; then
                    	v="${vlist[${ib}]}"
                    	E="${Elist[${ib}]}"
                    	echo "will execute awk on ${file} with v= ${v} and dir = ${E}"
                    	awk -v var=$time5ns -v v=$v '(NR>var) && ($6 < v) {print $1,$7,$6}' "../../${E}/${file}" > "./${file}" &
                    	PID+=( "$!" )
                    	ib=$(echo "$((${ib}+1))")
                    else
                    	echo "${file} already exists"
                    fi
        	done
                for pid in ${PID[*]}; do
        	    wait $pid
                done
        done
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
      [ -e "c${c}/c_${k}/metd.txt" ] && echo "c${c}/c_${k}/metd.txt exists, removing file first"
      echo "making metd.txt for c_${k}"
      for i in c$c/c_$k/colv*; do
         cnt=$(echo "$i" | grep -oP '(?<=colvar)[^_]*')
         cnt2="$(cut -d'_' -f3 <<<$i)"
	 read cntkappa cnt2kappa <<< $(grep -oPh '(?<=KAPPA=).*(?= )' ../../E_$cnt/plumed_* | tr "\n" " ")
         echo "$(basename $i)   $cnt   $cnt2   $cntkappa   $cnt2kappa" >> "c${c}/c_${k}/metd.txt"
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
      wham='/usr/users/jlapier/bin/wham/wham-2d/wham-2d'
   fi

   for i in `seq 0 4 $c`; do
      unset PID
      declare -a PID
      for k in `seq $(($i+1)) $(($i+4))`; do
         if cd "c${c}/c_${k}"; then
            echo "launching 2dwham in c_$k with: E1min=$E1min; E1max=$E1max; bin=$bin; E2min=$E2min; E2max=$E2max"
            ($wham Px=0 $E1min $E1max $bin Py=0 $E2min $E2max $bin 0.000001 300 0 metd.txt 2dpmf.txt 0 &> 2dwham.out) &
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
       sed '/inf/d' "c${c}/c_${k}/2dpmf.txt" > "c${c}/c_${k}/2dpmf_clean.txt"
       wait $!
       $(projection.py -f "c${c}/c_${k}/2dpmf_clean.txt" -o "c${c}/c_${k}/1dpmf_${k}.txt")
   done
   python -c "from block import block_avg; print(block_avg('c*/c_*/1dpmf*'))"
   echo "clean_wham2d done"
}

#column -t metd.txt > temp.txt
