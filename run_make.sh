#! /bin/sh
source activate sharpy

# This flag forces gfortran to print exactly when the
# print statement is, but has some computational overhead.
# Use only for debugging the library
# export GFORTRAN_UNBUFFERED_ALL=1

make install
