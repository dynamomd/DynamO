#A simple wrapper Makefile around the boost build system, also downloads and unpacks boost if required.
BOOST_VER=1.49.0
BOOST_DIR="boost_"$(subst .,_,$(BOOST_VER))
BOOST_FILE=$(BOOST_DIR)".tar.bz2"
BOOST_DL="http://sourceforge.net/projects/boost/files/boost/"$(BOOST_VER)"/"$(BOOST_FILE)
BJAM="./src/boost/bjam"

all : build_deps
	$(BJAM) install

debug : build_deps
	$(BJAM) debug

test : build_deps
	$(BJAM) test

docs :
	doxygen

clean: build_deps
	rm -Rf build-dir

#Make sure we've downloaded boost and built the bjam executable inside
build_deps:
	if [ ! -d ./src/boost ]; then wget $(BOOST_DL) && tar xf $(BOOST_FILE) && rm $(BOOST_FILE) && mv $(BOOST_DIR) src/boost; fi
	if [ ! -x $(BJAM) ]; then cd src/boost; ./bootstrap.sh; fi

install: build_deps all
	if [ -d ./bin ]; then mkdir -p $(DESTDIR)/usr/bin/; cp -R bin/* $(DESTDIR)/usr/bin/; fi
	if [ -d ./lib ]; then mkdir -p $(DESTDIR)/usr/lib/; cp -R lib/* $(DESTDIR)/usr/lib/; fi
	if [ -d ./include ]; then mkdir -p $(DESTDIR)/usr/include/; cp -R include/* $(DESTDIR)/usr/include/; fi

distclean:
	rm -Rf build-dir
	rm -Rf src/boost


.PHONY : build_deps all install distclean test docs