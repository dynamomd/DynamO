#!/bin/bash

dynamod="../bin/dynamod"
dynarun="../bin/dynarun"

NUMRUN=4
NCOLL=5000000
function testrun {
    > speedvals
    > memvals
    $dynamod -m 0 -d $dens -C $C --f1 $elas #--i2 2000000

    bzcat config.out.xml.bz2 | xmlstarlet ed -u '//Scheduler/Sorter/@Type' -v $sorter | bzip2 > config.tmp.xml.bz2
    mv config.tmp.xml.bz2 config.out.xml.bz2
    
    for i in $(seq 0 $NUMRUN); do
	echo -n "Running test $i for $C cells, $elas elasticity, and $dens density...."
	val=$($dynarun config.out.xml.bz2 -c $NCOLL | grep "Avg Coll" | gawk '{print $4}')
	echo $val >> speedvals
	bzcat output.xml.bz2  | xmlstarlet sel -t -v '//MemoryUsage/@ResidentSet' >> memvals
    done
    echo $C $elas $(cat speedvals | gawk 'BEGIN {sum=0; sqrsum=0} { sum += $1; sqrsum += $1*$1} END {print "Colls Avg "sum/NR" Dev "sqrt((sqrsum - sum * sum /NR) / NR)}') \
	$(cat memvals | gawk 'BEGIN {sum=0; sqrsum=0} { sum += $1; sqrsum += $1*$1} END {print "Mem Avg "sum/NR" Dev "sqrt((sqrsum - sum * sum /NR) / NR)}')
    echo $sorter $(cat speedvals | gawk 'BEGIN {sum=0; sqrsum=0} { sum += $1; sqrsum += $1*$1} END {print sum/NR, sqrt((sqrsum - sum * sum /NR) / NR)}') \
	$(cat memvals | gawk 'BEGIN {sum=0; sqrsum=0} { sum += $1; sqrsum += $1*$1} END {print sum/NR, sqrt((sqrsum - sum * sum /NR) / NR)}') \
	>> sorter.dat
    
    > speedvals
    > memvals
}

C=25
dens=0.5
elas=1.0
for sorter in BoundedPQMinMax2 BoundedPQMinMax3 BoundedPQMinMax4 BoundedPQMinMax5 BoundedPQMinMax6 BoundedPQMinMax7 BoundedPQMinMax8 BoundedPQ; do
    testrun
done
