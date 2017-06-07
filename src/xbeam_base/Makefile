# FFLAGS = -g3 -fPIC -O1 -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow
FFLAGS = -fPIC -O3 -funroll-loops -ftree-parallelize-loops=4 -march=native -fomit-frame-pointer -I/usr/include 
LINKER_FLAGS = -shared -llapack -fopenmp
UNAME := $(shell uname)
#ifeq ($(UNAME), Linux)
FC = gfortran
#endif
#ifeq ($(UNAME), Darwin)
#FC = ifort
#endif
LINKER = ar cr
object_dir = obj/
lib_dir = lib/

ifeq ($(UNAME), Linux)
target = libxbeam_base.so
endif
ifeq ($(UNAME), Darwin)
target = libxbeam_base.dylib
endif
files = lib_fem\
lib_rot\
lib_rotvect\
lib_bgeom\
lib_sparse\
lib_lu\
lib_cbeam3\
lib_solv\
lib_out\
lib_xbeam\
lib_perturb

objects = $(files:=.o)
sources = $(files:=.f90)
modules = $(files:=.mod)

default: $(target)

$(target): $(objects)
	$(FC) -o $(target) $(addprefix $(object_dir), $^) $(LINKER_FLAGS)
	mkdir -p $(lib_dir)
	mv $(target) $(lib_dir)$(target)
	mv *.mod $(object_dir)

%.o: %.f90
	$(FC) $(FFLAGS) -c $<
	mkdir -p $(object_dir)
	mv $@ $(object_dir)$@

clean:
	rm -rf $(object_dir)

clobber: clean
	rm -rf $(lib_dir)

.PHONY: default clean clobber
