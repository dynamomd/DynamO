#A simple wrapper Makefile around the boost build system, also downloads and unpacks boost if required.
BOOST_DIR=boost_1_47_0
BOOST_FILE=$(BOOST_DIR).tar.bz2
BOOST_DL="http://sourceforge.net/projects/boost/files/boost/1.47.0/boost_1_47_0.tar.bz2"
BJAM="./src/boost/bjam"

all : build_deps
	$(BJAM) -j4 install

debug : build_deps
	$(BJAM) -j4 debug

test : build_deps
	$(BJAM) -j4 test

docs :
	doxygen

clean: build_deps
	$(BJAM) clean

#Make sure we've downloaded boost and built the bjam executable inside
build_deps:
	if [ ! -d ./src/boost ]; then wget $(BOOST_DL) && tar xf $(BOOST_FILE) && rm $(BOOST_FILE) && mv $(BOOST_DIR) src/boost; fi
	if [ ! -x $(BJAM) ]; then cd src/boost; ./bootstrap.sh; fi

install: build_deps all
	if [ -d ./bin ]; then mkdir -p $(DESTDIR)/usr/bin/; cp -R bin/* $(DESTDIR)/usr/bin/; fi
	if [ -d ./lib ]; then mkdir -p $(DESTDIR)/usr/lib/; cp -R lib/* $(DESTDIR)/usr/lib/; fi
	if [ -d ./include ]; then mkdir -p $(DESTDIR)/usr/include/; cp -R include/* $(DESTDIR)/usr/include/; fi

distclean: build_deps
	rm -Rf build-dir
	rm -Rf src/boost

.PHONY : build_deps all install distclean test docs