#A simple wrapper Makefile around the boost build system, also downloads and unpacks boost if required.
BOOST_FILE=boost_1_46_1.tar.bz2
BOOST_DL="http://downloads.sourceforge.net/project/boost/boost/1.46.1/$(BOOST_FILE)"

all : build_deps
	./boost/bjam -j4 install

#Make sure we've downloaded boost and built the bjam executable inside
build_deps:
	if [ ! -d ./boost ]; then wget $(BOOST_DL) && tar xf $(BOOST_FILE) && rm $(BOOST_FILE) && mv boost_* boost; fi
	if [ ! -x ./boost/bjam ]; then cd boost; ./bootstrap.sh; fi

install: build_deps
	mkdir -p $(DESTDIR)/usr/bin/
	cp bin/* $(DESTDIR)/usr/bin/

distclean: build_deps
	rm -Rf build-dir

.PHONY : build_deps all install distclean