#!/bin/bash
#set -x

CC=g++
input_file=$1
shift # pull off first arg
args="$*"
filename=${input_file%%.cpp}

$CC -o $filename $input_file -O3 -m64 -march=native -msse4.2 -mavx2 -fopenmp
rc=$?

if [[ $rc == 0 ]]; then
   ./$filename $args
   exit $?
fi

exit $rc
