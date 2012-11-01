#!/bin/bash
VERSION=1.4-rc1
PACKAGE_NAME=dynamo
export DEBFULLNAME="Marcus Bannerman"
export DEBEMAIL="support@dynamomd.org"
OldWD=$(pwd)
BUILD_DIR=$OldWD/dynamo-build
SRC_DIR=$BUILD_DIR/"$PACKAGE_NAME"-$VERSION
###This requries quite a bit of setting up, install the dependencies
#sudo apt-get install pbuilder debootstrap devscripts
###Then set up the distributions to build for
#pbuilder-dist precise amd64 create
#pbuilder-dist precise i386 create

###Download the sources and prepare the upstream package
echo "Making the build dir ("$BUILD_DIR")"
mkdir -p $BUILD_DIR || exit 1
cd $BUILD_DIR; git clone https://github.com/toastedcrumpets/DynamO.git $SRC_DIR || exit 1
cd $SRC_DIR; rm -Rf .git
cd $SRC_DIR; make boost_dl

cp -R $OldWD/debian $SRC_DIR/

#Create the orig file
cd $BUILD_DIR; tar -czf $PACKAGE_NAME"_"$VERSION.orig.tar.gz $SRC_DIR/*
#Build the signed source package
cd $SRC_DIR; debuild -S

echo "Do you wish to test the build on amd64 and i386?"
select yn in "Yes" "No"; do
    case $yn in
        Yes ) cd $BUILD_DIR; pbuilder-dist precise i386 build $PACKAGE_NAME"_"$VERSION-1.dsc; pbuilder-dist precise amd64 build $PACKAGE_NAME"_"$VERSION-1.dsc; break;;
        No ) exit;;
    esac
done

echo "Do you wish to upload the build?"
select yn in "Yes" "No"; do
    case $yn in
        Yes ) cd $BUILD_DIR; dput ppa:dynamomd/ppa $PACKAGE_NAME"_"$VERSION-1_source.changes; break;;
        No ) exit;;
    esac
done

#############THE COMMANDS USED TO GENERATE THE PACKAGE
###Start building the package
#BUILDDEPENDS="debhelper (>= 8.0.0), curl, libbz2-dev, make, g++"
#DEPENDS="libbz2"
#dh_make -e $DEBEMAIL --createorig -p $PACKAGE_NAME --multi -c gpl3 || exit 1
#echo "Source: "$PACKAGE_NAME"
#Section: science
#Priority: optional
#Maintainer: Marcus Bannerman <support@dynamomd.org>
#Build-Depends: "$BUILDDEPENDS"
#Standards-Version: 3.9.2
#Homepage: http://dynamomd.org
#Vcs-Git: git://github.com/toastedcrumpets/DynamO.git
#
#Package: dynamo
#Architecture: any
#Depends: "$DEPENDS"
#Description: An event-driven molecular dynamics simulator 
# DynamO is a modern event-driven molecular dynamics simulator with a
# wide range of tools and models built in. It also includes a world
# class visualisation package for particulate systems." > debian/control
#
#cat $OldWD/changelog > debian/changelog
#
#rm -f debian/README.source debian/README.Debian debian/*.ex debian/dynamo?doc* || exit 1
#
##Prevent the test from running when building the package
#echo "override_dh_auto_test: " >> debian/rules
