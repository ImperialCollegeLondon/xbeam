## GFORTRAN SETTINGS
## ===========================================================================
### GFORTRAN RELEASE - FAST RUN, SLOW COMPILATION
export FCOMP = gfortran
export FCFLAGS = -fPIC -O3 -march=x86-64 -mtune=generic -funroll-loops -ftree-parallelize-loops=4 -fopenmp  -m64
export LDFLAGS = -shared  -Wl,--no-as-needed -lpthread -lm -ldl -fopenmp

### GFORTRAN DEBUG
#export FCOMP = gfortran
#export FCFLAGS = -g3 -fPIC -fcheck=all -fbacktrace -pedantic -fno-omit-frame-pointer  -ffpe-trap=invalid,zero,overflow,underflow,precision,denormal -std=f2008 -W -Wtabs -O0 -fbacktrace -fbounds-check -fstack-arrays -fno-underscoring 
#export LDFLAGS = -lblas -llapack -shared -lm -ldl
## ===========================================================================

## INTEL FORTRAN SETTINGS
## ===========================================================================
### IFORT RELEASE
###     with MKL
#export FCOMP = ifort
#export FCFLAGS = -fPIC -O3 -funroll-loops -heap-arrays -xHost -wrap-margin- -mkl=parallel -parallel -qopenmp
#export LDFLAGS = -lmatmul -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl -mkl -shared -parallel -qopenmp

### IFORT DEBUG
#export FCOMP = ifort
#export FCFLAGS = -fPIC -O0 -g3 -stand f08 -traceback -assume protect_parens  -implicitnone -check bounds -ftrapuv -debug all -fpe-all=0 -no-ftz -wrap-margin- -fp-model source
#export LDFLAGS = -llapack -shared
## ===========================================================================

## FLANG SETTINGS
## ===========================================================================
### FLANG DEBUG
#export FCOMP = flang
#export FCFLAGS = -O0 -g3 -fPIC
#export LDFLAGS = -llapack -shared

### FLANG RELEASE
#export FCOMP = flang
#export FCFLAGS = -O3 -fPIC -funroll-loops
#export LDFLAGS = -llapack -shared

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
