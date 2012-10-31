#!/bin/bash
BUILD_DIR=/tmp/dynamo_build/
DYNAMO_VER=1.4-rc1
GIT_BRANCH=master
REPOSITORY=
MAINTAINER="support@dynamomd.org"
MAINTAINER_NAME="Marcus Bannerman"
PACKAGE_REVISION="0"
LICENCE=gpl3
DEPENDS=", libbz2-1.0, zlib1g"
SRC_DEPENDS=", libbz2-dev, zlib1g-dev"
URL="http://www.dynamomd.org"

echo "Making the build dir ("$BUILD_DIR")"
mkdir -p $BUILD_DIR || exit 1
cd $BUILD_DIR
git clone https://github.com/toastedcrumpets/DynamO.git || exit 1
cd DynamO
make boost_dl || exit 1
rm -Rf .git
cd ../
tar -czf DynamO-$DYNAMO_VER.tar.gz DynamO || exit 1