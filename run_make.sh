#! /bin/sh

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

