#!/bin/bash
#set -x

REP_vals=(1000 10000 100000 1000000 10000000)
run_vals=( {1..1000} )
bash_file=$1
cpp_file=$2
n_cores=$3

for REP in ${REP_vals[@]}
do
    for run in ${run_vals[@]}
    do  
        ./$bash_file $cpp_file $n_cores $REP $run 
    done
done

exit $rc

