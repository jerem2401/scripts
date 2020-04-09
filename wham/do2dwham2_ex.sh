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


source 2dwham_2.sh

prep2dwham
mk_chunk
mk_metd
do_wham2d
clean_wham2d
