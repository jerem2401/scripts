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
            bin="$1"
            ;;
        -E2min)
            shift;
            E2min="$1"
            ;;
        -E2max)
            shift;
            E2max="$1"
            ;;
    esac
    shift
done

module load conda
source activate env1

source 2dwham_2.sh

prep2dwham 2>&1 | tee -a prep2dwham.log
mk_chunk 2>&1 | tee -a mk_chunk.log
mk_metd 2>&1 | tee -a mk_metd.log
do_wham2d | tee -a do_wham2d.log
clean_wham2d | tee -a clean_wham2d.log
