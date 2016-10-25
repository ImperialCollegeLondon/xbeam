default:
	$(MAKE) -C src/

clean:
	cd src/xbeam_base; make clean
	cd src; make clean
