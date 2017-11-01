TARGET = base
FILES = lib_fem\
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

OBJS = $(FILES:=.o)
SRCS = $(FILES:=.f90)
MODS = $(FILES:=.mod)

all: $(TARGET)

$(TARGET): $(OBJS)

%.o: %.f90
	$(FCOMP) $(FCFLAGS) -c $<

.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~

