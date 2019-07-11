## GFORTRAN SETTINGS
## ===========================================================================
export FCOMP = gfortran
## GFORTRAN RELEASE
export FCFLAGS = -fPIC -O3 -funroll-loops -ftree-parallelize-loops=4 -march=native -fopenmp # -fno-underscoring -g -fno-omit-frame-pointer
export LDFLAGS = -L$(LAPACK_LIB_DIR) -llapack -shared -fopenmp
## GFORTRAN DEBUG
# export FCFLAGS = -g -fPIC -fcheck=all -fbacktrace -pedantic -fno-omit-frame-pointer  -ffpe-trap=invalid,zero,overflow,underflow,precision,denormal -std=f2008 -W -Wtabs -O -fbacktrace -fbounds-check -fstack-arrays -fno-underscoring
# export LDFLAGS = -L$(LAPACK_LIB_DIR) -llapack -shared -fopenmp
## ===========================================================================

## INTEL FORTRAN SETTINGS
## ===========================================================================
export FCOMP = ifort
### IFORT RELEASE
### no MKL
#export FCFLAGS = -I${MKLROOT}/include -fPIC -O3 -funroll-loops -march=native -fopenmp -heap-arrays -xHost -wrap-margin- -parallel -mkl
### with MKL
# export FCFLAGS = -fPIC -O3 -funroll-loops -march=native -heap-arrays -xHost -wrap-margin- -mkl=parallel -parallel -qopenmp
# export LDFLAGS = -lmatmul -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl -mkl -shared -parallel -qopenmp
#`export LDFLAGS = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl -shared

## IFORT DEBUG
export FCFLAGS = -fpic -O0 -g3 -stand f08 -traceback  -fstack-protector  -assume protect_parens  -implicitnone -check bounds -ftrapuv -debug all -fpe-all=0 -no-ftz -wrap-margin-
export LDFLAGS = -L$(LAPACK_LIB_DIR) -llapack -shared
## ===========================================================================

FOLDER = src/
LIBRARY_NAME =  libxbeam
UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
EXT=so
endif
ifeq ($(UNAME), Darwin)
EXT=dylib
endif

LIBRARY = $(LIBRARY_NAME:=.$(EXT))
export LIBRARY

all: $(LIBRARY)

install: $(LIBRARY)
	mv ./lib/$(LIBRARY) ../sharpy/lib/$(LIBRARY)

$(LIBRARY):
	$(MAKE) -C src
	mkdir -p ./lib
	mv $(FOLDER)$(LIBRARY) ./lib/$(LIBRARY)

.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD
	rm -f ./lib/$(LIBRARY)
	$(MAKE) -C src clean

veryclean: clean
	rm -f *~ $(LIBRARY)
	$(MAKE) -C src veryclean
