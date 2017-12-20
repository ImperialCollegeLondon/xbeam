export FCOMP = gfortran
export FCFLAGS = -fPIC -O3 -funroll-loops -ftree-parallelize-loops=4 -march=native -fomit-frame-pointer -ffast-math
# export FCFLAGS = -g -Og -fPIC -fcheck=all -fbacktrace -pedantic # -Wall

# export FCOMP = ifort
# export FCFLAGS = -fpic -O0 -g -stand f08    -assume realloc_lhs  -check bounds -traceback  -fstack-protector  -assume protect_parens  -implicitnone -check uninit -ftrapuv -debug all -fpe-all=0 -no-ftz
# export FCFLAGS = -fPIC -O3 -funroll-loops -march=native -fomit-frame-pointer -xHost -fast

export LDFLAGS = -llapack -shared -fopenmp
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
