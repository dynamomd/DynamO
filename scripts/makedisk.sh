#!/bin/bash

bin/dynamod -m0 --i1 2 -C 25 -x 4 -d 1 --text --rectangular-box 
bzcat config.out.xml.bz2 | xmlstarlet ed \
    -d "//Pt/P[@x>=0.5]/.."  \
    -d "//Pt/P[@x < '-0.5']/.." \
    -u "//Pt/V/@x" -v 0 \
    -u "//Pt/P/@x" -v 0 \
    | bzip2 > config.plane.xml.bz2

MaxR=$(bzcat config.plane.xml.bz2 | xmlstarlet sel -t -v '//Units/@BoxLength' | gawk '{print $1*0.5}')

R2Outer=$(echo $MaxR | gawk '{print ($1)**2}')

R2Inner=$(echo $MaxR | gawk '{print (0.5*$1)**2}')

bzcat config.plane.xml.bz2 | xmlstarlet ed \
    -d "//Pt/P[@x*@x+@y*@y+@z*@z >= $R2Outer]/.." \
    -d "//Pt/P[@x*@x+@y*@y+@z*@z < $R2Inner]/.." \
    | bzip2 > config.annuli.xml.bz2


