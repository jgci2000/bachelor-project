#!/bin/bash
#set -x

exp_f_vals=( {0..13} )
run_vals=( {1..1000} )
bash_file=$1
cpp_file=$2

for REP in ${exp_f_vals[@]}
do
    for run in ${run_vals[@]}
    do  
        ./$bash_file $cpp_file $REP $run 
    done
done

exit $rc

