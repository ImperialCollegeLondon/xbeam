#! /bin/sh

# This flag forces gfortran to print exactly when the
# print statement is, but has some computational overhead.
# Use only for debugging the library
# echo "GFORTRAN_UNBUFFERED_ALL is 1"
# export GFORTRAN_UNBUFFERED_ALL=1

make install
