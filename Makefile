#A simple wrapper Makefile around the boost build system, also downloads and unpacks boost if required.
BOOST_VER=1.50.0
BOOST_DIR=boost_$(subst .,_,$(BOOST_VER))
BOOST_FILE=$(BOOST_DIR).tar.bz2
BOOST_DL=http://sourceforge.net/projects/boost/files/boost/$(BOOST_VER)/$(BOOST_FILE)
BJAM=./src/boost/bjam

all: boost_bjam
	echo "### Building release version of DynamO"
	$(BJAM) install toolset=gcc

debug: boost_bjam
	echo "### Building debug version of DynamO"
	$(BJAM) debug toolset=gcc

test: boost_bjam
	echo "### Testing DynamO software"
	$(BJAM) test toolset=gcc

docs:
	echo "### Building DynamO documentation"
	doxygen

clean:
	rm -Rf build-dir

boost_dl:
	echo "### Testing if the boost sources have been downloaded:"
	if [ ! -f ./$(BOOST_FILE) ]; then echo "### Boost sources are missing, downloading boost sources:"; curl -L -o $(BOOST_FILE) $(BOOST_DL) || wget $(BOOST_DL); fi
	if [ ! -f ./$(BOOST_FILE) ]; then echo "### Failed to download the boost sources"; exit 1; fi

boost_untar: boost_dl
	echo "### Testing if the boost sources have been extracted:"
	if [ -f ./$(BOOST_FILE) ] && [ ! -d ./src/boost ]; then echo "### Extracting boost sources"; tar xf $(BOOST_FILE) && mv $(BOOST_DIR) src/boost; fi
	if [ -f ./$(BOOST_FILE) ] && [ ! -d ./src/boost ]; then echo "### Failed to untar the boost sources"; exit 1; fi

#Make sure we've downloaded boost and built the bjam executable inside
boost_bjam: boost_untar
	echo "### Testing if the build system has been setup"
	if [ ! -x $(BJAM) ]; then echo "### Setting up the build system"; cd src/boost; ./bootstrap.sh; else echo "### Build system setup"; fi

install: all
	if [ -d ./bin ]; then mkdir -p $(DESTDIR)/usr/bin/; cp -R bin/* $(DESTDIR)/usr/bin/; fi
	if [ -d ./lib ]; then mkdir -p $(DESTDIR)/usr/lib/; cp -R lib/* $(DESTDIR)/usr/lib/; fi
	if [ -d ./include ]; then mkdir -p $(DESTDIR)/usr/include/; cp -R include/* $(DESTDIR)/usr/include/; fi

distclean:
	rm -Rf build-dir src/boost lib/ include/ bin/dynarun bin/dynamod bin/dynahist_rw


.PHONY: boost_untar boost_dl boost_bjam all install distclean test docs
.SILENT: boost_untar boost_dl boost_bjam install all debug test docs clean distclean
