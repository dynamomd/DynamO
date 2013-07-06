#A simple wrapper Makefile around the boost build system
all: 
	echo "### Building release version of DynamO"
	mkdir -p build-dir
	bjam install toolset=gcc

debug: 
	echo "### Building debug version of DynamO"
	bjam debug toolset=gcc

test:
	echo "### Testing DynamO software"
	bjam test toolset=gcc

docs:
	echo "### Building DynamO documentation"
	doxygen

clean:
	rm -Rf build-dir

install: all
	if [ -d ./bin ]; then mkdir -p $(DESTDIR)/usr/bin/; cp -R bin/* $(DESTDIR)/usr/bin/; fi
	if [ -d ./lib ]; then mkdir -p $(DESTDIR)/usr/lib/; cp -R lib/* $(DESTDIR)/usr/lib/; fi
	if [ -d ./include ]; then mkdir -p $(DESTDIR)/usr/include/; cp -R include/* $(DESTDIR)/usr/include/; fi

distclean:
	rm -Rf build-dir lib/ include/ bin/


.PHONY: all install distclean test docs
.SILENT: install all debug test docs clean distclean
