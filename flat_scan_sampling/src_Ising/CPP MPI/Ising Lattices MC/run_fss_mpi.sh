#!/bin/bash
#set -x

CC=mpic++
input_file=$1
shift # pull off first arg
n_cores=$1
shift
args="$*"
filename=${input_file%%.cpp}


$CC -o $filename $CFLAGS $input_file -O3 -m64 -march=native -msse4.2 -mavx2
rc=$?

if [[ $rc == 0 ]]; then
   mpirun -np $n_cores ./$filename $args
   exit $?
fi

exit $rc
