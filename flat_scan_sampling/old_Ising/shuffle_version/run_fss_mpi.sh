#!/bin/bash

#
# Script to compile and run FSS method
# João Inácio, Mar. 25th, 2021
#


# Compiler
CC=mpic++

# File names
input_file=$1
filename=${input_file%%.cpp}

# Number of processes
shift
n_cores=$1

# Runtime arguments
shift
args="$*"

# Compilation
$CC -o $filename -Ofast -std=c++17 -march=native -m64 $input_file Ising.cpp Fss_Functions.cpp
rc=$?

# Run the program
if [[ $rc == 0 ]]; then
   mpiexec -np $n_cores ./$filename $args
   exit $?
fi

exit $rc
