#!/bin/bash
VERSION=1.4-rc1
PACKAGE_NAME=dynamo
GIT_BRANCH=master
REPOSITORY=
export DEBFULLNAME="Marcus Bannerman"
export DEBEMAIL="support@dynamomd.org"
BUILDDEPENDS=", curl, libbz2-dev, make, g++"
DEPENDS=", libbz2"
URL="http://www.dynamomd.org"
ARCHIVE_NAME=$PACKAGE_NAME-$VERSION.tar.gz
BUILD_DIR=/tmp/dynamo-build/"$PACKAGE_NAME"-$VERSION
DISTRIBUTION=precise


###Download the sources and prepare the upstream package
echo "Making the build dir ("$BUILD_DIR")"
mkdir -p $BUILD_DIR || exit 1
cd $BUILD_DIR
git clone https://github.com/toastedcrumpets/DynamO.git $BUILD_DIR || exit 1
cd $BUILD_DIR
rm -Rf .git
###Start building the package
dh_make -e $DEBEMAIL --createorig -p $PACKAGE_NAME --multi -c gpl3 || exit 1

echo "Source: "$PACKAGE_NAME"
Section: science
Priority: optional
Maintainer: Marcus Bannerman <support@dynamomd.org>
Build-Depends: debhelper (>= 8.0.0) "$BUILDDEPENDS"
Standards-Version: 3.9.2
Homepage: http://dynamomd.org
Vcs-Git: git://github.com/toastedcrumpets/DynamO.git

Package: dynamo
Architecture: any
Depends: "$DEPENDS"
Description: An event-driven molecular dynamics simulator 
 DynamO is a modern event-driven molecular dynamics simulator with a
 wide range of tools and models built in. It also includes a world
 class visualisation package for particulate systems." > debian/control

rm -f debian/README.source debian/README.Debian debian/*.ex debian/dynamo?doc* || exit 1

#Prevent the test from running when building the package
echo "override_dh_auto_test: " >> debian/rules

#Test building the package!
debuild -us -uc -S