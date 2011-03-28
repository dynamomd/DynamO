
all : 
	cd boost; ./bootstrap.sh
	boost/bjam -j4 install

install:
	cp bin/* /usr/bin/