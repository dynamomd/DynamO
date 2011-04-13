#!/bin/bash
dynamod="../bin/dynamod"
dynarun="../bin/dynarun"

NUMRUN=4
NCOLL=5000000

function testrun {
    > speedvals
    > memvals
    $dynamod -m 0 -d $dens -C $C --f1 $elas #--i2 2000000
    
    for i in $(seq 0 $NUMRUN); do
	echo -n "Running test $i for $C cells, $elas elasticity, and $dens density...."
	val=$($dynarun config.out.xml.bz2 -c $NCOLL | grep "Avg Coll" | gawk '{print $4}')
	echo $val >> speedvals
	bzcat output.xml.bz2  | xmlstarlet sel -t -v '//MemoryUsage/@ResidentSet' >> memvals
    done
    echo $C $elas $(cat speedvals | gawk 'BEGIN {sum=0; sqrsum=0} { sum += $1; sqrsum += $1*$1} END {print "Colls Avg "sum/NR" Dev "sqrt((sqrsum - sum * sum /NR) / NR)}') \
	$(cat memvals | gawk 'BEGIN {sum=0; sqrsum=0} { sum += $1; sqrsum += $1*$1} END {print "Mem Avg "sum/NR" Dev "sqrt((sqrsum - sum * sum /NR) / NR)}')
    echo $C $elas $(cat speedvals | gawk 'BEGIN {sum=0; sqrsum=0} { sum += $1; sqrsum += $1*$1} END {print sum/NR, sqrt((sqrsum - sum * sum /NR) / NR)}') \
	$(cat memvals | gawk 'BEGIN {sum=0; sqrsum=0} { sum += $1; sqrsum += $1*$1} END {print sum/NR, sqrt((sqrsum - sum * sum /NR) / NR)}') \
	>> dens$dens.morton.dat
    
    > speedvals
    > memvals
    
    $dynamod -m 0 -d $dens --b2 -C $C --f1 $elas #--i2 2000000
    
    for i in $(seq 0 $NUMRUN); do
	echo -n "Running test $i for $C cells, $elas elasticity, and $dens density...."
	val=$($dynarun config.out.xml.bz2 -c $NCOLL | grep "Avg Coll" | gawk '{print $4}')
	echo $val >> speedvals
	bzcat output.xml.bz2  | xmlstarlet sel -t -v '//MemoryUsage/@ResidentSet' >> memvals
    done
    echo $C $elas $(cat speedvals | gawk 'BEGIN {sum=0; sqrsum=0} { sum += $1; sqrsum += $1*$1} END {print "Colls Avg "sum/NR" Dev "sqrt((sqrsum - sum * sum /NR) / NR)}') \
	$(cat memvals | gawk 'BEGIN {sum=0; sqrsum=0} { sum += $1; sqrsum += $1*$1} END {print "Mem Avg "sum/NR" Dev "sqrt((sqrsum - sum * sum /NR) / NR)}')
    echo $C $elas $(cat speedvals | gawk 'BEGIN {sum=0; sqrsum=0} { sum += $1; sqrsum += $1*$1} END {print sum/NR, sqrt((sqrsum - sum * sum /NR) / NR)}') \
	$(cat memvals | gawk 'BEGIN {sum=0; sqrsum=0} { sum += $1; sqrsum += $1*$1} END {print sum/NR, sqrt((sqrsum - sum * sum /NR) / NR)}') \
	>> dens$dens.dat
}

for elas in 1.0; do # 0.4 0.5 0.6 0.7 0.8 0.9 
    for dens in 0.8; do
	#> dens$dens.dat
	for C in 75 80 85 90 95 100 105 110 115 120 125 130 135 ; do 
	    testrun
	done
    done
done

for elas in 1.0; do # 0.4 0.5 0.6 0.7 0.8 0.9 
    for dens in 0.5 1.1; do
	#> dens$dens.dat
	for C in  5 6 7 8 9 10 12 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 105 110 115 120 125 130 135 ; do 
	    testrun
	done
    done
done