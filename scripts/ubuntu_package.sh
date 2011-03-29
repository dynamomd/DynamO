#!/bin/bash
BUILD_DIR=/tmp/dynamo_build/
BOOST_FILE=boost_1_46_1.tar.bz2
DYNAMO_VER=1.0
GIT_BRANCH=dynamo-1-0
MAINTAINER="dynamo@marcusbannerman.co.uk"
MAINTAINER_NAME="Marcus Bannerman"
PACKAGE_REVISION="0"
LICENCE=gpl3
DEPENDS=", libbz2-1.0, zlib1g"
SRC_DEPENDS=", libbz2-dev, zlib1g-dev"
URL="http://www.marcusbannerman.co.uk/dynamo"

############################################################################
################################ FUNCTIONS #################################
############################################################################
function files {
cd $PACKAGENAME-$DYNAMO_VER
DEBFULLNAME=$MAINTAINER_NAME dh_make \
    -e $MAINTAINER \
    -f ../$PACKAGENAME"_"$DYNAMO_VER.orig.tar.gz \
    --multi \
    -c $LICENCE
rm debian/*.ex debian/*.EX

mkdir -p debian/source
echo "3.0 (native)" > debian/source/format
> debian/info #Blank any info files from being created

echo $PACKAGENAME" ("$DYNAMO_VER-$PACKAGE_REVISION") "$(lsb_release -c -s)"; urgency=low

  * Direct build of the DynamO simulator, see the git repository or
    website http://www.marcusbannerman.co.uk for change information.

 -- "$MAINTAINER_NAME" <"$MAINTAINER">  "$(date -R) > debian/changelog

echo "Source: "$PACKAGENAME"
Section: science
Priority: extra
Maintainer: "$MAINTAINER_NAME" <"$MAINTAINER">
Build-Depends: debhelper (>= 7)"$SRC_DEPENDS"
Standards-Version: 3.8.3
Homepage: "$URL"

Package: "$PACKAGENAME"
Architecture: any
Depends: "'${shlibs:Depends}'", "'${misc:Depends}'$DEPENDS"
Conflicts: "$CONFLICTS"
Description: A general event-driven particle simulator.
 A general event-driven simulator capable of simulating millions of
 particles in real time." > debian/control

echo "This work was packaged for Debian by:

    "$MAINTAINER_NAME" <"$MAINTAINER"> on "$(date -R)"

It was downloaded from:

    <url://"$URL">

Upstream Author(s):

    Marcus Bannerman <marcus@marcusbannerman.co.uk>

Copyright:

    <Copyright (C) 2011 Marcus Bannerman>

License:

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This package is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

On Debian systems, the complete text of the GNU General
Public License version 3 can be found in '/usr/share/common-licenses/GPL-3'.

The Debian packaging is:

    Copyright (C) 2011 "$MAINTAINER_NAME" <"$MAINTAINER">

and is licensed under the GPL version 3, see above.
" > debian/copyright
cd ../
}

############################################################################
################################ Build System ##############################
############################################################################

###Make/Clean the build directory
mkdir -p $BUILD_DIR
cd $BUILD_DIR
rm -Rf *

###Download the source, clean off the git cruft and add in the boost libraries
git clone git://marcusbannerman.co.uk/dynamo.git
cd dynamo
git checkout $GIT_BRANCH
cd ..
rm -Rf dynamo/.git
mv dynamo dynamo-$DYNAMO_VER
cp -R dynamo-$DYNAMO_VER dynaview-$DYNAMO_VER
#Patch the dynamo package to not build the visualizer
echo "WIBBLE" > dynamo-$DYNAMO_VER/src/coil/tests/openclpp.cpp

#Make the source tarballs
tar czf dynamo_$DYNAMO_VER.orig.tar.gz dynamo-$DYNAMO_VER
tar czf dynaview_$DYNAMO_VER.orig.tar.gz dynaview-$DYNAMO_VER

#Source file is created! Begin building a package
############################################################################
####WITHOUT VISUALIZATION SUPPORT
############################################################################
PACKAGENAME="dynamo"
CONFLICTS="dynaview"
files

############################################################################
####WITH VISUALIZATION SUPPORT
############################################################################
PACKAGENAME="dynaview"
CONFLICTS="dynamo"

DEPENDS=$DEPENDS", libgtkmm-2.4-1c2a, freeglut3, libftgl2 "
SRC_DEPENDS=$SRC_DEPENDS", libgtkmm-2.4-dev, libgl1-mesa-dev, freeglut3-dev, libftgl-dev"

files

echo " ." >> $PACKAGENAME-$DYNAMO_VER/debian/control
echo " This version of dynamo is built with visualization support enabled." >> $PACKAGENAME-$DYNAMO_VER/debian/control
