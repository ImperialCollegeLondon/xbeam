INCLUDE = -Ixbeam_base

FILES = xbeam_shared\
debug_utils\
cbeam3_asbly\
cbeam3_solv\
xbeam_asbly\
xbeam_fdiff\
xbeam_perturb\
xbeam_solv\
xbeam_undef\
cbeam3_interface\
xbeam_interface

OBJS = $(FILES:=.o)
SRCS = $(FILES:=.f90)
MODS = $(FILES:=.mod)

all: $(LIBRARY)

$(LIBRARY): $(OBJS)
	$(FCOMP) $(FCFLAGS) $(INCLUDE) -o $@ $^ $(LDFLAGS)

%.o: %.f90 XBEAM_BASE
	$(FCOMP) $(FCFLAGS) $(INCLUDE) -c $<

.PHONY: XBEAM_BASE

XBEAM_BASE:
	$(MAKE) -C xbeam_base

.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD
	$(MAKE) -C xbeam_base clean


veryclean: clean
	rm -f *~ $(LIBRARY)

