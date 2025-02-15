#!/bin/bash

#
# Script to compile and run FSS method
# João Inácio, Mar. 25th, 2021
#


# Compiler
CC=g++

# File names
input_file=$1
filename=${input_file%%.cpp}

# Runtime arguments
shift
args="$*"

# Compilation
$CC -o $filename -Ofast -std=c++17 -m64 -march=native $input_file Ising.cpp Fss_Functions.cpp
rc=$?

# Run the program
if [[ $rc == 0 ]]; then
   ./$filename $args
   exit $?
fi

exit $rc
