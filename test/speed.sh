#!/bin/bash
dynamod="/home/mjki2mb2/dynamo/bin/dynamod"
dynarun="/home/mjki2mb2/dynamo/bin/dynarun"

NUMRUN=4
NCOLL=10000000

for dens in 0.5; do
    > dens$dens.dat
    for C in  75 100 125 150 160 175; do #5 6 7 8 9 10 15 20 25 30 40 50
	> speedvals
	> memvals
	$dynamod -m 0 -d $dens -C $C

	for i in $(seq 0 $NUMRUN); do
	    echo -n "Running test $i for $C cells and $dens density...."
	    val=$($dynarun config.out.xml.bz2 -c $NCOLL | grep "Avg Coll" | gawk '{print $4}')
	    echo $val >> speedvals
	    bzcat output.xml.bz2  | xmlstarlet sel -t -v '//MemoryUsage/@ResidentSet' >> memvals
	done
	echo $C $(cat speedvals | gawk 'BEGIN {sum=0; sqrsum=0} { sum += $1; sqrsum += $1*$1} END {print "Colls Avg "sum/NR" Dev "sqrt((sqrsum - sum * sum /NR) / NR)}') \
	    $(cat memvals | gawk 'BEGIN {sum=0; sqrsum=0} { sum += $1; sqrsum += $1*$1} END {print "Mem Avg "sum/NR" Dev "sqrt((sqrsum - sum * sum /NR) / NR)}')
	echo $C $(cat speedvals | gawk 'BEGIN {sum=0; sqrsum=0} { sum += $1; sqrsum += $1*$1} END {print sum/NR, sqrt((sqrsum - sum * sum /NR) / NR)}') \
	    $(cat memvals | gawk 'BEGIN {sum=0; sqrsum=0} { sum += $1; sqrsum += $1*$1} END {print sum/NR, sqrt((sqrsum - sum * sum /NR) / NR)}') \
	    >> dens$dens.morton.dat

	> speedvals
	> memvals
	$dynamod -m 0 -d $dens --b2 -C $C

	mv config.tmp.xml.bz2 config.out.xml.bz2
	for i in $(seq 0 $NUMRUN); do
	    echo -n "Running test $i for $C cells and $dens density...."
	    val=$($dynarun config.out.xml.bz2 -c $NCOLL | grep "Avg Coll" | gawk '{print $4}')
	    echo $val >> speedvals
	    bzcat output.xml.bz2  | xmlstarlet sel -t -v '//MemoryUsage/@ResidentSet' >> memvals
	done
	echo $C $(cat speedvals | gawk 'BEGIN {sum=0; sqrsum=0} { sum += $1; sqrsum += $1*$1} END {print "Colls Avg "sum/NR" Dev "sqrt((sqrsum - sum * sum /NR) / NR)}') \
	    $(cat memvals | gawk 'BEGIN {sum=0; sqrsum=0} { sum += $1; sqrsum += $1*$1} END {print "Mem Avg "sum/NR" Dev "sqrt((sqrsum - sum * sum /NR) / NR)}')
	echo $C $(cat speedvals | gawk 'BEGIN {sum=0; sqrsum=0} { sum += $1; sqrsum += $1*$1} END {print sum/NR, sqrt((sqrsum - sum * sum /NR) / NR)}') \
	    $(cat memvals | gawk 'BEGIN {sum=0; sqrsum=0} { sum += $1; sqrsum += $1*$1} END {print sum/NR, sqrt((sqrsum - sum * sum /NR) / NR)}') \
	    >> dens$dens.non_morton.dat

    done
done
