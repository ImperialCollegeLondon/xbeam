#! /bin/sh
# export GFORTRAN_UNBUFFERED_ALL=1
cd src/xbeam_base/
cd obj
rm *
cd ..
make
cd ..
cd obj
rm *
cd ..
make
cd ..
make
cp lib/* ../SHARPy/lib/
