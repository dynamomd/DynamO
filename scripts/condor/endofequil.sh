#!/bin/bash
if [ ! -e $2/Equilibration/Run$1/config.out.xml ]; then
    exit 1
fi

cp $2/Equilibration/Run$1/config.out.xml \
    $2/Production/Run0/config.out.xml