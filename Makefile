
all : 
	cd boost; ./bootstrap.sh
	boost/bjam -j4 install

install:
	cp bin/* /usr/bin/

distclean:
	boost/bjam -j4 clean
	rm -Rf build-dir