#A simple wrapper Makefile around the boost build system, also downloads and unpacks boost if required.
BOOST_DIR=boost_1_46_1
BOOST_FILE=$(BOOST_DIR).tar.bz2
BOOST_DL="http://downloads.sourceforge.net/project/boost/boost/1.46.1/$(BOOST_FILE)"
BJAM="./src/boost/bjam"

all : build_deps
	$(BJAM) -j4 install

debug : build_deps
	$(BJAM) -j4 debug

test : build_deps
	$(BJAM) -j4 test

docs :
	doxygen

#Make sure we've downloaded boost and built the bjam executable inside
build_deps:
	if [ ! -d ./src/boost ]; then wget $(BOOST_DL) && tar xf $(BOOST_FILE) && rm $(BOOST_FILE) && mv $(BOOST_DIR) src/boost; fi
	if [ ! -x $(BJAM) ]; then cd src/boost; ./bootstrap.sh; fi

install: build_deps all
	if [ -d ./bin ]; then mkdir -p $(DESTDIR)/usr/bin/; cp bin/* $(DESTDIR)/usr/bin/; fi
	if [ -d ./lib ]; then mkdir -p $(DESTDIR)/usr/lib/; cp lib/* $(DESTDIR)/usr/lib/; fi
	if [ -d ./include ]; then mkdir -p $(DESTDIR)/usr/include/; cp include/* $(DESTDIR)/usr/include/; fi

distclean: build_deps
	rm -Rf build-dir

.PHONY : build_deps all install distclean test docs