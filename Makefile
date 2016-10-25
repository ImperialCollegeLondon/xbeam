default:
	$(MAKE) -C src/
	mkdir -p lib
	mv src/lib/*.so lib/
	mv src/xbeam_base/lib/*.so lib/
	mkdir -p include
	mv src/obj/*.mod include/
	mv src/xbeam_base/obj/*.mod include/

clean:
	cd src/xbeam_base; make clean
	cd src; make clean
