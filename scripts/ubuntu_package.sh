#!/bin/bash
VERSION=1.4-rc4
PACKAGE_REV=2
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
echo "Cloning a clean copy of the DynamO sources"
cd $BUILD_DIR; git clone $OldWD/../ $SRC_DIR || exit 1
cd $SRC_DIR; rm -Rf .git || exit 1

#Determine if the boost sources are already downloaded in the parent
echo "Looking for a copy of the boost libraries"
regex="BOOST_VER=([0-9.]*)"
[[ $(cat $SRC_DIR/Makefile) =~ $regex ]]
BOOST_VER=${BASH_REMATCH[1]}
BOOST_FILE=boost_${BOOST_VER//./_}.tar.bz2
if [ -f $OldWD/../$BOOST_FILE ]; then
    echo "Found a local copy of the boost sources ("$OldWD/../$BOOST_FILE")"
    cp $OldWD/../$BOOST_FILE $SRC_DIR/ || exit 1
else
    echo "Asking the makefile to download a copy of the boost sources"
    cd $SRC_DIR; make boost_dl || exit 1
fi

#Copy the package configuration across
echo "Copying the package configuration across"
cp -R $OldWD/debian $SRC_DIR/ || exit 1

#Create the orig file
echo "Creating a tar of the sources"
cd $SRC_DIR; tar -czf ../$PACKAGE_NAME"_"$VERSION.orig.tar.gz . || exit 1
#Build the signed source package
echo "Creating the source package"
cd $SRC_DIR; debuild -S || exit 1

#echo "Do you wish to test the build on amd64 and i386?"
#select yn in "Yes" "No"; do
#    case $yn in
#        Yes ) cd $BUILD_DIR; pbuilder-dist precise i386 build $PACKAGE_NAME"_"$VERSION-$PACKAGE_REV.dsc || exit 1; cd $BUILD_DIR pbuilder-dist precise amd64 build $PACKAGE_NAME"_"$VERSION-$PACKAGE_REV.dsc || exit 1; break;;
#        No ) exit;;
#    esac
#done
#
#echo "Do you wish to upload the build?"
#select yn in "Yes" "No"; do
#    case $yn in
#        Yes ) cd $BUILD_DIR; dput ppa:dynamomd/ppa $PACKAGE_NAME"_"$VERSION-$PACKAGE_REV"_"source.changes; break;;
#        No ) exit;;
#    esac
#done
