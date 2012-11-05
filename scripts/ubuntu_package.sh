#!/bin/bash
VERSION=1.4.0
PACKAGE_REV=1
PACKAGE_NAME=dynamo
export DEBFULLNAME="Marcus Bannerman"
export DEBEMAIL="support@dynamomd.org"
OldWD=$(pwd)
BUILD_DIR=$OldWD/dynamo-build
SRC_DIR=$BUILD_DIR/"$PACKAGE_NAME"-$VERSION

###Download the sources and prepare the upstream package
echo "Making the build dir ("$BUILD_DIR")"
mkdir -p $BUILD_DIR || exit 1
echo "Cloning a clean copy of the DynamO sources"
cd $BUILD_DIR; git clone $OldWD/../ $SRC_DIR || exit 1
cd $SRC_DIR; rm -Rf .git || exit 1

#Create the orig file
echo "Creating a tar of the sources"
cd $SRC_DIR; tar -czf ../$PACKAGE_NAME"_"$VERSION.orig.tar.gz . || exit 1
#Build the signed source package
echo "Creating the source package"
cd $SRC_DIR; debuild -S || exit 1
