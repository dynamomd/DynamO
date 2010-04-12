#!/bin/bash

bzcat config.out.xml.bz2 \
    | xmlstarlet sel -t -m '//PotentialDeformation/Entry' -v '@Energy' -o " " \
    -v '@Shift' -n | gawk '{print $1+0.0,$2+0.0}' | head -n-1  > tmp1.dat

bzcat output.xml.bz2 | xmlstarlet sel -t -v '//EnergyHist/WeightHistogram' \
    | head -n-1 | tail -n+2 > tmp2.dat


data=`gawk 'pass == 1 {Emc[$1+0.0] = $2} pass == 2 {if ($1+0.0 in Emc) Emc[$1+0.0] += $2; else Emc[$1+0.0] = $2} END {for (x in Emc) {print "<Entry Energy=\""x"\" Shift=\""Emc[x]"\"/>"}}' pass=1 tmp1.dat pass=2 tmp2.dat`
bzcat config.out.xml.bz2 | xmlstarlet ed -P -d '//Entry' -u '//PotentialDeformation' -v "$data" | xmlstarlet unesc | bzip2 > config.opt.xml.bz2