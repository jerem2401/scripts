#!/usr/bin/bash

while [ $# -gt 0 ]; do
    case "$1" in
        -f)
            shift
            rpath="$1"
	    rpath=$(eval echo ${rpath}/colv*)
            ;;
        -c)
            shift
            c="$1"
	        ;;
        -h)
            echo "TODO: fill up the help message"
             ;;
        -E1min)
            shift
            E1min="$1"
            ;;
        -E1max)
            shift
            E1max="$1"
            ;;
        -bin)
            shift;
            bin="$1";;
	-column) shift
	    column="$1";;
    esac
    shift
done

module load anaconda3/2020.07 && source activate env1

source 1dwham.sh

prep2dwham 2>&1 | tee -a prep1dwham.log
mk_chunk 2>&1 | tee -a mk_chunk.log
mk_metd 2>&1 | tee -a mk_metd.log
do_wham2d | tee -a do_wham1d.log
clean_wham2d | tee -a clean_wham2d.log
