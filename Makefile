UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
EXT = so
endif
ifeq ($(UNAME), Darwin)
EXT = dylib
endif
default:
	$(MAKE) -C src/
	mkdir -p lib
	mv src/lib/*.$(EXT) lib/
	mv src/xbeam_base/lib/*.$(EXT) lib/
	mkdir -p include
	mv src/obj/*.mod include/
	mv src/xbeam_base/obj/*.mod include/

clean:
	cd src/xbeam_base; make clean
	cd src; make clean
