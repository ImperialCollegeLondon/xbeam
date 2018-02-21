export FCOMP = gfortran
export FCFLAGS = -fPIC -O3 -funroll-loops -ftree-parallelize-loops=4 -march=native -fomit-frame-pointer -ffast-math
# export FCFLAGS = -fPIC -O3 -funroll-loops -march=native -fcheck=all
# export FCFLAGS = -g -fPIC -fcheck=all -fbacktrace -pedantic -fno-omit-frame-pointer  -ffpe-trap=invalid,zero,overflow,underflow,precision,denormal -std=f2008 -W -Wtabs -O -fbacktrace -fbounds-check -fstack-arrays # -fsanitize=address -static-libasan # -Wall

export GFORTRAN_UNBUFFERED_ALL = 1
# export FCOMP = ifort
# export FCFLAGS = -fpic -O0 -g -stand f08 -traceback  -fstack-protector  -assume protect_parens  -implicitnone -check bounds -ftrapuv -debug all -fpe-all=0 -no-ftz -wrap-margin-
# export FCFLAGS = -fPIC -O3 -funroll-loops -march=native -fomit-frame-pointer -wrap-margin- # -Xhost -fast

export LDFLAGS = -llapack -shared # -lasan -lubsan -static-libasan -fsanitize=address
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
