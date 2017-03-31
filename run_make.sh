#! /bin/sh

# This flag forces gfortran to print exactly when the
# print statement is, but has some computational overhead.
# Use only for debugging the library
# export GFORTRAN_UNBUFFERED_ALL=1

cd src/xbeam_base/
mkdir -p obj
cd obj
rm *
cd ..
make
cd ..
mkdir -p obj
cd obj
rm *
cd ..
make
cd ..
make
cp lib/* ../sharpy/lib/
