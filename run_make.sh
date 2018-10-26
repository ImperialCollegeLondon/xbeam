#! /bin/sh

# This flag forces gfortran to print exactly when the
# print statement is, but has some computational overhead.
# Use only for debugging the library
# echo "GFORTRAN_UNBUFFERED_ALL is 1"
# export GFORTRAN_UNBUFFERED_ALL=1

# lapack flag
export LAPACK_LIB_DIR=$(conda info --json | python -c "import sys, json; print(json.load(sys.stdin)['active_prefix'])")/lib

make install
