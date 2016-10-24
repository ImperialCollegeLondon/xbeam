FFLAGS = -g -fPIC -O3
FC = gfortran
LINKER = ar cr
object_dir = obj/
lib_dir = lib/

target = libxbeam.a
base_dir = xbeam_base/
base_lib = libxbeam_base.a

files = xbeam_shared\
cbeam3_asbly\
xbeam_shared\
cbeam3_solv\
interface_lapack_v3\
xbeam_asbly\
xbeam_fdiff\
xbeam_perturb\
xbeam_solv\
xbeam_undef

objects = $(files:=.o)
sources = $(files:=.f90)
modules = $(files:=.mod)

libs = xbeam_base
libs_dir = $(base_dir)lib/
include_dir = $(base_dir)obj/

default: $(target)

$(target): $(base_lib) $(objects)
	$(LINKER) $(target) $(addprefix $(object_dir), $(objects))
	mkdir -p $(lib_dir)
	mv $(target) $(lib_dir)$(target)
	mv *.mod $(object_dir)

%.o: %.f90
	$(FC) $(FFLAGS) -I$(include_dir) -c $<
	mkdir -p $(object_dir)
	mv $@ $(object_dir)$@

$(base_lib):
	$(MAKE) -C $(base_dir)

clean:
	rm -rf $(object_dir)

clobber: clean
	rm -rf $(lib_dir)

.PHONY: default clean clobber

